C======================================================================
C     energy.F
C     Computational function that takes a scalar and doubles it.
C     This is a MEX-file for MATLAB.
C  [v lengths E dE] = opti_lbfgs(v,k,strain,kappa,L,N,f_c,eps)
C======================================================================


C     Gateway routine
      subroutine relax(N,L,v,k,strain,sts,eps)

C     Declarations
      implicit none



C     Arguments for computational routine:
      integer::N,L
	  double precision,dimension(0:2*N-1)::v
	  double precision, dimension(0:N-1,0:2) ::k
	  
      real(8),dimension(0:2*N-1)::dE
      real(8),dimension(:,:),allocatable::w
      real(8)::strain,kappa,E,sts(0:1,0:1)
      real(8)::eps,xtol,f_c,eps_small

      real(8),dimension(:),allocatable::work,diag
      integer::iprint(2),iflag,icall,i,j,imax,jmax,fmax

      DOUBLE PRECISION GTOL,STPMIN,STPMAX
      integer :: NN,M,MP,LP,indicator
      logical::diagco
      EXTERNAL LB2
      COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX

      real(8)::dt,vv,ff,vf,f,acoef,mass
      real(8),dimension(:),allocatable::md_velocity
      integer::nstep

	  kappa = 0.d0
        
      iprint(1)=0
      iprint(2)=0
      icall = 0
      iflag = 0
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
      M = 5
	  
      call Energy(v,k,kappa,L,strain,E,N)
      call DEnergy(v,k,kappa,L,strain,dE,N,sts)
	  
	  	 

  51  CONTINUE
  20  CONTINUE
      v = v + dt * md_velocity - 0.5 * dE * dt * dt / mass
      md_velocity = md_velocity - 0.5 * dE * dt / mass
      call DEnergy(v,k,kappa,L,strain,dE,N,sts)
      md_velocity(0:2*N-1)=md_velocity(0:2*N-1)-0.5*dE(0:2*N-1)*dt/mass
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

    
C      call LBFGS(NN,M,v,E,dE,diagco,diag,iprint,+
C     +              eps_small,xtol,work,iflag)
      icall = icall +1   
C      if(iflag .lt. 0) go to 50   
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
    
     
  
