program fracture
	!A fast fortran program to record every bond breaking event in the simulation of fracture.
	!Input is parameters: L,f_s,p,strain_max,strain_step
	!Output is event logs: the format is # broken bond , strain , stress drop
	
	implicit none
	
	double precision f_s,p,strain_max,strain_step
	integer L,N
	double precision,allocatable::r(:),k(:,:),lt(:,:),rm(:)
	double precision strain
	double precision stress(0:1,0:1)
	double precision eps,sp,sn,sd,lts(1000),xb(1000),yb(1000),xc,yc
	integer flag,bc,sc,cc(1000),adjc,nc
	double precision sdp
	character*250 filename,filename2
	integer i,j,c,m,seed,ii,jj
	double precision b,x,y,z(2)
	double precision fortn
	
	
	read(*,*) L,f_s,p,strain_max,strain_step
	read(*,*) filename
	read(*,*) filename2
	read(*,*) seed
	read(*,*) eps
	
	!Assign the initial values of position,spring const,et. al.
	Do i = 1, seed
		call random_number(b)
	enddo
	call init_random_seed(seed)
	flag = 0
	nc = 0
	sd = strain_step
	N = L*L
	allocate(k(0:N-1,0:2))
	allocate(lt(0:N-1,0:2))
	allocate(r(0:2*N-1))
	allocate(rm(0:2*N-1))
	strain = 0
	DO 10 i = 0 , L - 1
		DO 20 j = 0 , L - 1
			m = i + j * L
			x = dble(i - 0.5*j)
			y = dble(sqrt(3.)/2 * j)
			r(m) = x
			r(m+N) = y
			call random_number(b)
			if(b < p) k(m,0) = 1.
			if( j < (L-1)) then
				call random_number(b)
				if(b < p) k(m,1) = 1.
				call random_number(b)
				if(b < p) k(m,2) = 1.
			endif
		20 continue
	10 continue
	Do 30 i = N / 2 + L / 2 - 1 , N / 2 + L / 2 + 2
		k(i,1) = 0.
		k(i,2) = 0.
	30 continue
	rm = r			
	!Do an optimization to minimize energy

	35 continue 
	call relax(N,L,rm,k,strain,stress,1.d-5)
	call bonds(rm,k,L,strain,lt,N)
! 	print *,'relax,bondortn:',fortn,',strain:',strain,',maxlt:',maxval(lt)
	!Compute bond elongation, break bonds, 
	!eps = 1.d-4
	if( flag .ne. 1) then
		z=maxloc(lt)
		i=z(1)-1
		j=z(2)-1
		b = 1.d0 + f_s 
		x = b + eps
		if(lt(i,j) > x .and. k(i,j) > 0.1 ) then
			!choose smaller step and redo relaxation
			flag = 2
			adjc = adjc + 1
			strain = sp + (sn - sp) / 2.d0
			sn = strain 
			print *,'adjust,strain:',strain,sp
			sd = (sn - sp) / 4.5d0
			if(sd < 1.e-10) then
				stop
			endif
			rm = r
			go to 35
		endif
	endif
	b = 1.d0 + f_s 
	z=maxloc(lt)
	i=z(1)-1
	j=z(2)-1
	if(lt(i,j) > b .and. k(i,j) > 0.1) then
		!break the bond and redo relaxation
		if(flag .ne. 1) then
			sdp = stress(1,1)
		endif
		k(i,j) = 0.
		r=rm
		flag = 1
		bc=bc+1
		nc = nc + 1
		if(bc .eq. 1) print *,'new event'
		print *,'break at strain',strain,'max_length',lt(i,j)
		lts(bc) = lt(i,j)
		x =dble(mod(i,L))
		y = dble(int(i/L))
		xb(bc) = - 1.0 * dble(y) * 0.5 + dble(x) 
		yb(bc) = sqrt(3.)/2. * dble(y)
		if(j .eq. 0) then
			xb(bc) = xb(bc) + 0.5
		elseif(j .eq. 1) then
			xb(bc) = xb(bc) - 0.25
			yb(bc) = yb(bc) + sqrt(3./4.)
		elseif(j .eq. 2) then
			xb(bc) = xb(bc) + 0.25
			yb(bc) = yb(bc) + sqrt(3./4.)
		endif
			
		sc = 0
		do ii = 0,N-1
			do jj = 0,2
				if(abs(lt(ii,jj)-1.d0) > 0.5 * f_s .and. k(ii,jj) > 0.1) then
					sc = sc + 1
				endif
			enddo
		enddo
		cc(bc) = sc
		go to 35
	endif
	!no bond break, increase strain
	
	if(flag .eq. 0) then
		sp = strain
		sd = sd * 2.2d0
		strain = strain + sd
	else
		sd = strain_step
		sp = strain
		strain = strain + sd
		if(flag .eq. 1) then
! 			x = 0
! 			y = 0
! 			xc = 0
! 			yc = 0
! 			do i = 1, bc
! 				do j = 1, bc
! ! 			y = y + (xb(i)-xb(j) -dble(L)*dnint(xb(i)/dble(L)-xb(j)/dble(L)))**2.d0
! ! 			y = y + (yb(i)-yb(j))**2.d0
! 					if( j .eq. (i + 1)) then
! 			y = (xb(i)-xb(j) -dble(L)*dnint(xb(i)/dble(L)-xb(j)/dble(L)))**2.d0
! 			y = y + (yb(i)-yb(j))**2.d0
! 			x = x + dsqrt(y)
! 					endif
! 				enddo
! 			enddo
! 			if(bc .ge. 2) then
! ! 				y = dsqrt(y / dble(bc - 1)/ dble(bc))
! 				x = x / dble(bc - 1)
! 			else
! 				x = 0.d0
! 				y = 0.d0
! 			endif
! 			xc = xb(1)
! 			yc = yb(1)
			open(file=filename,unit=32,access='append')
			write(32,*)bc,strain,sdp
			close(unit=32) 
			open(file=filename2,unit=33,access='append')
			do i = 1, bc
				write(33,'(2F9.3,3I8)')xb(i),yb(i),cc(i),bc,nc
			enddo
			close(unit=33) 
		endif
	endif
	sn = strain
	flag = 0
	r=rm
	bc = 0
	sdp = 0.d0
	if(strain < strain_max) go to 35
	print *,'end'
end program

subroutine init_random_seed(a)
	implicit none
	integer, allocatable :: seed(:)
	integer a,t,n
	integer :: i,dt(8)
	
	call random_seed(size = n)
	allocate(seed(n))
	
	call date_and_time(values = dt)
	t = (dt(5)) * 60 * 60 + dt(7)  + dt(8) 
	t = t + dt(6) * 60 
	t = t + a
	do i = 1, n
		seed(i) = mod(t * i,49973)
	enddo
	call random_seed(put=seed)
end subroutine
	
	
	
	
