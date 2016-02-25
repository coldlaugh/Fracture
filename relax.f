C======================================================================
C     energy.F
C     Computational function that takes a scalar and doubles it.
C     This is a MEX-file for MATLAB.
C  [v lengths E dE] = opti_lbfgs(v,k,strain,kappa,L,N,f_c,eps)
C======================================================================


      subroutine relax(N,L,v,k,strain,sts,eps)
      implicit none

      integer::N,L
	  double precision,dimension(0:2*N-1)::v
	  double precision, dimension(0:N-1,0:2) ::k

      real(8),dimension(0:2*N-1)::dE
      real(8),dimension(:,:),allocatable::w
      real(8)::strain,kappa,E,sts(0:1,0:1)
      real(8)::eps,xtol,f_c,eps_small
      integer::iprint(2),iflag,icall,i,j,imax,jmax,fmax
      integr::ir

      real(8)::dt,vv,ff,vf,f,acoef,mass
      real(8),dimension(:),allocatable::md_velocity
      real(8)::u
      integer::nstep

	    kappa = 0.d0
      allocate(md_velocity(0:2*N-1))
      md_velocity(:) = 0.d0
      dt = 1.d-4
      f = 0.5
      acoef = f
      nstep = 0
      mass = 0.1
      xtol = 1.d-15;
      eps_small = 1.d-10;
      NN = 2*N
      dE(:) =0.d0

  20  CONTINUE
      moveBoundary(v,L,strain,N)
      do ir = 1,N
          call random_number(u)
          i = int(N * u)
          v(i) = v(i) + dt * md_velocity(ir) - 0.5 * dE(i) * dt * dt / mass
          md_velocity(i) = md_velocity(i) - 0.5 * dE(i) * dt / mass
          v(i+N)=v(i+N)+dt*md_velocity(i+N)-0.5*dE(i+N)*dt*dt/mass
          md_velocity(i+N) = md_velocity(i+N) - 0.5 * dE(i+N) * dt / mass
          call singleDEnergy(v,k,L,N,nsx,nsy,dE(i),dE(i+N))
          md_velocity(i)=md_velocity(i)-0.5*dE(i)*dt/mass
          md_velocity(i+N) = md_velocity(i+N) - 0.5*dE(i+N)*dt/mass
      enddo
      if(maxval(dabs(dE(:))) .le. eps) go to 50
      vf = - DOT_PRODUCT(dE , md_velocity)
      vv = DOT_PRODUCT(md_velocity , md_velocity)
      ff = DOT_PRODUCT(dE , dE)
      if(vf < 0) then
           md_velocity(:) = 0.d0
           dt = dt * 0.5
           acoef = f
           nstep = 0
      else
           nstep = nstep + 1
           if(nstep > 5) then
                dt = min(dt*1.1,5.d-3)
                acoef = acoef * 0.99
           endif
      endif
        do i = 0,2*N-1
        md_velocity(i)=md_velocity(i)*(1.d0-acoef)+
     +                  (-acoef)*dsqrt(vv/ff)*dE(i)
        enddo

      icall = icall +1
      if(icall .gt. 800000) then
        print *, 'not finished'
		stop
        go to 50
      endif
      go to 20
  50  CONTINUE

      call DEnergy(v,k,kappa,L,strain,dE,N,sts)



      return
      end
