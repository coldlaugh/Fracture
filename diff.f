      subroutine Energy(v,k,kappa,L,strain,E,N)
C      spring energy, spring const k,  and bending energy, const kappa
      implicit none
      integer ::N
      double precision, dimension(0:2*N-1) :: v
      integer ::L
      double precision, dimension(0:N-1,0:2) ::k
      double precision ::strain,kappa
      double precision :: E
      integer :: i,j,i1,j1,m,m0,m1,m2,ll,m3,m4,m5,im1,jm1
      double precision :: l0,l1,l2,dx,dy,c,c0,c1,eps,dpp
	  double precision dxm,dym,l3,l4,l5,theta
      double precision, dimension(0:N-1) :: ux,uy
      double precision :: Ebend

      c0=sqrt(3.)/2.
      c=c0*(1.+strain)*dble(L-1)
      c1=c0*(1.+strain)*dble(L)
      ll=L*(L-1)
      E=0.
      Ebend=0.
      eps=0.00001
      do m=0,N-1
      ux(m)=v(m)
      uy(m)=v(m+N)
      enddo
      do i=0,L-1
      do j=0,1
         m=i+j*L*(L-1)
         E=E+1.d3*(v(m+N)-dble(j)*c)*(v(m+N)-dble(j)*c)
      enddo
      enddo

      do i=0,L-1
      do j=0,L-1
        m=i+j*L
        i1=i+1
        if (i1==L) i1=0
        j1=j+1
        if (j1==L) j1=0
        im1=i-1
        if (im1.lt.0) im1=L-1
        jm1=j-1
        if (jm1.lt.0) jm1=L-1
        m0=i1+j*L
        m1=i+j1*L
        m2=i1+j1*L
        m3=im1+j*L
        m4=i+jm1*L
        m5=im1+jm1*L
        l0=1.
        if (k(m,0).gt.eps) then
           dx=ux(m0)-ux(m)
           dy=uy(m0)-uy(m)
C           dx=dx +0.5*dble(L)*dnint(dy/c1)
           dx=dx - L * dnint(dx/L)
C           dy= dy - c1*dnint(dy/c1)
           l0=sqrt(dx*dx+dy*dy)
           E=E+k(m,0)*((l0)-1.)**2
           if (k(m3,0).gt.2.0) then
              dxm=ux(m)-ux(m3)
              dym=uy(m)-uy(m3)
C              dxm=dxm +0.5*dble(L)*dnint(dym/c1)
              dxm=dxm - L * dnint(dxm/L)
C              dym=dym - c1 * dnint(dym/c1)
              l3=sqrt(dxm*dxm+dym*dym)
              dpp=min((dxm*dx+dym*dy)/(l0*l3),1.d0)
              theta=dacos(dpp)
              E=E+kappa*theta*theta
           endif
        endif
        l1=1.
        if (k(m,1).gt.eps) then
           dx=ux(m1)-ux(m)
           dy=uy(m1)-uy(m)
C           dx=dx +0.5*dble(L)*dnint(dy/c1)
           dx=dx - L * dnint(dx/L)
C          added rounding by leyou
C           dy= dy - c1*dnint(dy/c1)
           l1=sqrt(dx*dx+dy*dy)
           E=E+k(m,1)*((l1)-1.)**2
           if (k(m4,1).gt.2.0) then
              dxm=ux(m)-ux(m4)
              dym=uy(m)-uy(m4)
C              dxm=dxm +0.5*dble(L)*dnint(dym/c1)
              dxm=dxm - L * dnint(dxm/L)
C              dym= dym - c1*dnint(dym/c1)
              l4=sqrt(dxm*dxm+dym*dym)
              dpp=min((dxm*dx+dym*dy)/(l1*l4),1.d0)
              theta=dacos(dpp)

              E=E+kappa*theta*theta
           endif
        endif
        l2=1.
        if (k(m,2).gt.eps) then
           dx=ux(m2)-ux(m)
           dy=uy(m2)-uy(m)
C           dx=dx +0.5*dble(L)*dnint(dy/c1)
           dx=dx - L * dnint(dx/L)
C          added rounding by leyou
C           dy= dy - c1*dnint(dy/c1)
           l2=sqrt(dx*dx+dy*dy)
           E=E+k(m,2)*((l2)-1.)**2
           if (k(m5,2).gt.2.0) then
              dxm=ux(m)-ux(m5)
              dym=uy(m)-uy(m5)
C              dxm=dxm +0.5*dble(L)*dnint(dym/c1)
              dxm=dxm - L * dnint(dxm/L)
C              dym= dym - c1*dnint(dym/c1)
              l5=sqrt(dxm*dxm+dym*dym)
              dpp=min((dxm*dx+dym*dy)/(l2*l5),1.d0)
              theta=dacos(dpp)

              E=E+kappa*theta*theta
           endif
        endif
      enddo
      enddo
      do i=0,L-1
      	do j=0,1
         	m=i+j*L*(L-1)
         	E=E+1.d3*(v(m+N)-(dble(j)*c))**2
      	enddo
      enddo

      E=E/2.
      end subroutine Energy



      subroutine DEnergy(v,k,kappa,L,strain,dE,N,sts)
        implicit none
        integer, intent(in)::N
        real(8), intent(inout), dimension(0:2*N-1) :: v
        integer, intent(in) ::L
        real(8), intent(in), dimension(0:N-1,0:2) ::k
        real(8), intent(in) ::strain, kappa
        real(8), intent(inout), dimension(0:2*N-1) :: dE
        integer :: i,j,i1,j1,im1,jm1,m,m0,m1,m2,m3,m4,m5,ll
        real(8) :: dxp,dyp,dEx,dEy,dEE,c0,c,c1,eps,den,dxm,dym
        real(8)::  dpp,theta,dott,lp,lm,dEbx1,dEby1,dEbx2,dEby2
        real(8)::  sts(0:1,0:1)

        real(8), dimension(0:N-1) :: ux,uy
        c0=sqrt(3.)/2.
        c=c0*(1.+strain)*(L-1)
        c1=c0*(1.+strain)*dble(L)
        ll=L*(L-1)
        eps=0.0000000001d0
        sts(:,:)=0.d0
        dE(:)=0.d0
        do m=0,N-1
           ux(m)=v(m)
           uy(m)=v(m+N)
        enddo

        do j=0,L-1
        do i=0,L-1
C           !i=mod(m,L)
C           !j=(m/L)
           m=i+j*L
           i1=i+1
           if (i1==L) i1=0
           j1=j+1
           if (j1==L) j1=0
           im1=i-1
           if (im1.lt.0) im1=L-1
           jm1=j-1
           if (jm1.lt.0) jm1=L-1


           m0=i1+j*L
           m1=i+j1*L
           m2=i1+j1*L
           m3=im1+j*L
           m4=i+jm1*L
           m5=im1+jm1*L


           dEx=0.
           dEy=0.
           lp=1.
           if (k(m,0).gt.eps) then
              dxp=ux(m0)-ux(m)
              dyp=uy(m0)-uy(m)
C              dxp=dxp +0.5*dble(L)*dnint(dyp/c1)
              dxp=dxp - L * dnint(dxp/L)
C              dyp= dyp - c1*dnint(dyp/c1)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m,0)*(lp-1.)/lp
              sts(0,0)=sts(0,0) + dEE*dxp*dxp/2.0
              sts(0,1)=sts(0,1) + dEE*dxp*dyp/2.0
              sts(1,0)=sts(1,0) + dEE*dxp*dyp/2.0
              sts(1,1)=sts(1,1) + dEE*dyp*dyp/2.0
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
              lm=1.
              if (k(m3,0).gt.2.0) then
                 dxm=ux(m)-ux(m3)
                 dym=uy(m)-uy(m3)
C                 dxm=dxm +0.5*dble(L)*dnint(dym/c1)
                 dxm=dxm - L * dnint(dxm/L)
C                 dym= dym-c1*dnint(dym/c1)
                 lm=sqrt(dxm*dxm+dym*dym)
                 dott=(dxm*dxp+dym*dyp)
                 den=sqrt(1.-(dott/lp/lm)**2)
                 if (den.ge.eps) then
                    dpp=min(dott/lp/lm,1.d0)
                    theta=dacos(dpp)
                    dEbx1= -theta*kappa*( (dxp)/lp/lm+
     +                      (-dxm/lm**3/lp)*dott )/den
                    dEbx2= -theta*kappa*( (-dxm)/lp/lm+
     +                      (dxp/lp**3/lm)*dott )/den
                    dEx=dEx +dEbx1+dEbx2
                    dEby1= -theta*kappa*( (dyp)/lp/lm+
     +                      (-dym/lm**3/lp)*dott )/den
                    dEby2= -theta*kappa*( (-dym)/lp/lm+
     +                      (dyp/lp**3/lm)*dott )/den
                    dEy=dEy+dEby1+dEby2
                    dE(m3) = dE(m3) - dEbx1
                    dE(m3+N) = dE(m3+N) - dEby1
                    dE(m0) = dE(m0) -dEbx2
                    dE(m0+N) = dE(m0+N) -dEby2

                    sts(0,0)=sts(0,0)+(-dEbx2)*dxp+(dEbx1)*dxm
                    sts(0,1)=sts(0,1)+(-dEbx2)*dyp+(dEbx1)*dym
                    sts(1,0)=sts(1,0)+(-dEby2)*dxp+(dEby1)*dxm
                    sts(1,1)=sts(1,1)+(-dEby2)*dyp+(dEby1)*dym

                 endif
              endif



           endif

           lp=1.
           if (k(m,1).gt.eps) then
              dxp=ux(m1)-ux(m)
              dyp=uy(m1)-uy(m)
C              dxp=dxp +0.5*dble(L)*dnint(dyp/c1)
              dxp=dxp - L * dnint(dxp/L)
C              dyp= dyp - c1*dnint(dyp/c1)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m,1)*(lp-1.)/lp
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
              sts(0,0)=sts(0,0) + dEE*dxp*dxp/2.0
              sts(0,1)=sts(0,1) + dEE*dxp*dyp/2.0
              sts(1,0)=sts(1,0) + dEE*dxp*dyp/2.0
              sts(1,1)=sts(1,1) + dEE*dyp*dyp/2.0
              lm=1.
              if (k(m4,1).gt.2.0) then
                 dxm=ux(m)-ux(m4)
                 dym=uy(m)-uy(m4)
C                 dxm=dxm +0.5*dble(L)*dnint(dym/c1)
                 dxm=dxm - L * dnint(dxm/L)

C                 dym= dym-c1*dnint(dym/c1)
                 lm=sqrt(dxm*dxm+dym*dym)
                 dott=(dxm*dxp+dym*dyp)
                 den=sqrt(1.-dott**2/lp**2/lm**2)
                 if (den.ge.eps) then
                    dpp=min(dott/lp/lm,1.d0)
                    theta=dacos(dpp)
                    dEbx1= -theta*kappa*( (dxp)/lp/lm+
     +                      (-dxm/lm**3/lp)*dott )/den
                    dEbx2= -theta*kappa*( (-dxm)/lp/lm+
     +                      (dxp/lp**3/lm)*dott )/den
                    dEx=dEx +dEbx1+dEbx2
                    dEby1= -theta*kappa*( (dyp)/lp/lm+
     +                      (-dym/lm**3/lp)*dott )/den
                    dEby2= -theta*kappa*( (-dym)/lp/lm+
     +                      (dyp/lp**3/lm)*dott )/den
                    dEy=dEy+dEby1+dEby2

                    dE(m4) = dE(m4) - dEbx1
                    dE(m4+N) = dE(m4+N) - dEby1
                    dE(m1) = dE(m1) -dEbx2
                    dE(m1+N) = dE(m1+N) -dEby2

                    sts(0,0)=sts(0,0)+(-dEbx2)*dxp+(dEbx1)*dxm
                    sts(0,1)=sts(0,1)+(-dEbx2)*dyp+(dEbx1)*dym
                    sts(1,0)=sts(1,0)+(-dEby2)*dxp+(dEby1)*dxm
                    sts(1,1)=sts(1,1)+(-dEby2)*dyp+(dEby1)*dym

                 endif
              endif

           endif

           lp=1.
           if (k(m,2).gt.eps) then
              dxp=ux(m2)-ux(m)
              dyp=uy(m2)-uy(m)
C              dxp=dxp +0.5*dble(L)*dnint(dyp/c1)
              dxp=dxp - L * dnint(dxp/L)
C              dyp= dyp - c1*dnint(dyp/c1)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m,2)*(lp-1.)/lp
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
              sts(0,0)=sts(0,0) + dEE*dxp*dxp/2.0
              sts(0,1)=sts(0,1) + dEE*dxp*dyp/2.0
              sts(1,0)=sts(1,0) + dEE*dxp*dyp/2.0
              sts(1,1)=sts(1,1) + dEE*dyp*dyp/2.0
              lm=1.
              if (k(m5,2).gt.2.0) then
                 dxm=ux(m)-ux(m5)
                 dym=uy(m)-uy(m5)
C                 dxm=dxm +0.5*dble(L)*dnint(dym/c1)
                 dxm=dxm - L * dnint(dxm/L)
C                 dym= dym-c1*dnint(dym/c1)
                 lm=sqrt(dxm*dxm+dym*dym)
                 dott=(dxm*dxp+dym*dyp)
                 den=sqrt(1.-dott**2/lp**2/lm**2)
                 if (den.ge.eps) then
                    dpp=min(dott/lp/lm,1.d0)
                    theta=dacos(dpp)
                    dEbx1= -theta*kappa*( (dxp)/lp/lm+
     +                      (-dxm/lm**3/lp)*dott )/den
                    dEbx2= -theta*kappa*( (-dxm)/lp/lm+
     +                      (dxp/lp**3/lm)*dott )/den
                    dEx=dEx +dEbx1+dEbx2
                    dEby1= -theta*kappa*( (dyp)/lp/lm+
     +                      (-dym/lm**3/lp)*dott )/den
                    dEby2= -theta*kappa*( (-dym)/lp/lm+
     +                      (dyp/lp**3/lm)*dott )/den
                    dEy=dEy+dEby1+dEby2

                    dE(m5) = dE(m5) - dEbx1
                    dE(m5+N) = dE(m5+N) - dEby1
                    dE(m2) = dE(m2) -dEbx2
                    dE(m2+N) = dE(m2+N) -dEby2

                    sts(0,0)=sts(0,0)+(-dEbx2)*dxp+(dEbx1)*dxm
                    sts(0,1)=sts(0,1)+(-dEbx2)*dyp+(dEbx1)*dym
                    sts(1,0)=sts(1,0)+(-dEby2)*dxp+(dEby1)*dxm
                    sts(1,1)=sts(1,1)+(-dEby2)*dyp+(dEby1)*dym


                 endif
              endif

           endif

           lp=1.
           if (k(m3,0).gt.eps) then
              dxp=ux(m3)-ux(m)
              dyp=uy(m3)-uy(m)
C              dxp=dxp +0.5*dble(L)*dnint(dyp/c1)
              dxp=dxp - L * dnint(dxp/L)
C              dyp= dyp - c1*dnint(dyp/c1)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m3,0)*(lp-1.)/lp
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
              sts(0,0)=sts(0,0) + dEE*dxp*dxp/2.0
              sts(0,1)=sts(0,1) + dEE*dxp*dyp/2.0
              sts(1,0)=sts(1,0) + dEE*dxp*dyp/2.0
              sts(1,1)=sts(1,1) + dEE*dyp*dyp/2.0
           endif

           lp=1.
           if (k(m4,1).gt.eps) then
              dxp=ux(m4)-ux(m)
              dyp=uy(m4)-uy(m)
C              dxp=dxp +0.5*dble(L)*dnint(dyp/c1)
              dxp=dxp - L * dnint(dxp/L)
C              dyp= dyp - c1*dnint(dyp/c1)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m4,1)*(lp-1.)/(lp)
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
              sts(0,0)=sts(0,0) + dEE*dxp*dxp/2.0
              sts(0,1)=sts(0,1) + dEE*dxp*dyp/2.0
              sts(1,0)=sts(1,0) + dEE*dxp*dyp/2.0
              sts(1,1)=sts(1,1) + dEE*dyp*dyp/2.0

           endif

           lp=1.
           if (k(m5,2).gt.eps) then
              dxp=ux(m5)-ux(m)
              dyp=uy(m5)-uy(m)
C              dxp=dxp +0.5*dble(L)*dnint(dyp/c1)
              dxp=dxp - L * dnint(dxp/L)
C              dyp= dyp - c1*dnint(dyp/c1)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m5,2)*(lp-1.)/lp
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
              sts(0,0)=sts(0,0) + dEE*dxp*dxp/2.0
              sts(0,1)=sts(0,1) + dEE*dxp*dyp/2.0
              sts(1,0)=sts(1,0) + dEE*dxp*dyp/2.0
              sts(1,1)=sts(1,1) + dEE*dyp*dyp/2.0
           endif
           dE(m)=dE(m)+dEx
           dE(m+N)=dE(m+N)+dEy
        enddo
        enddo

        do i=0,L-1
        do j=0,1
           m=i+j*L*(L-1)
		   dE(m+N) = 0.d0
           v(m+N)=dble(j)*c
        enddo
        enddo

       sts(0,0)=sts(0,0)/dble(L)/c1
       sts(0,1)=sts(0,1)/dble(L)/c1
       sts(1,0)=sts(1,0)/dble(L)/c1
       sts(1,1)=sts(1,1)/dble(L)/c1


      end subroutine DEnergy

      subroutine bonds(v,k,L,strain,lengths,N)
        implicit none
        integer, intent(in)::N
        real(8), intent(in), dimension(0:2*N-1) :: v
        integer, intent(in) ::L
        real(8), intent(in), dimension(0:N-1,0:2) ::k
        real(8), intent(in) ::strain
        real(8), intent(inout), dimension(0:N-1,0:2) ::lengths
        integer :: i,j,i1,j1,m,m0,m1,m2,ll
        real(8) :: l0,l1,l2,dx,dy,c,c0,c1
        real(8), dimension(0:N-1) :: ux,uy
        c0=sqrt(3.)/2.
        c=c0*(1.+strain)*(L-1)
        c1=c0*(1.+strain)*dble(L)
        ll=L*(L-1)

        do m=0,N-1
           ux(m)=v(m)
           uy(m)=v(m+N)
        enddo

        do m = 0, N-1
           i=mod(m,L)
           j=(m/L)
           i1=mod(i+1,L)
           j1=j+1
           if (j1.ge.L) j1=0
           m0=i1+j*L
           m1=i+j1*L
           m2=i1+j1*L
           dx=ux(m0)-ux(m)
           dy=uy(m0)-uy(m)
           dx=dx - L * dnint(dx/L)
           l0=dsqrt(dx*dx+dy*dy)
           lengths(m,0)=l0*k(m,0)
           dx=ux(m1)-ux(m)
           dy=uy(m1)-uy(m)
           dx=dx - L * dnint(dx/L)
           l1=dsqrt(dx*dx+dy*dy)
           lengths(m,1)=l1*k(m,1)
           dx=ux(m2)-ux(m)
           dy=uy(m2)-uy(m)
           dx=dx - L * dnint(dx/L)
           l2=dsqrt(dx*dx+dy*dy)
           lengths(m,2)=l2*k(m,2)
        enddo

      end subroutine bonds

      subroutine singleDEnergy(v,k,L,N,nsx,nsy,dEx,dEy)
        implicit none
        integer, intent(in)::N
        real(8), intent(inout), dimension(0:2*N-1) :: v
        integer, intent(in) ::L
        real(8), intent(in), dimension(0:N-1,0:2) ::k
        real(8), intent(inout) :: dEx,dEy
        integer :: nsx,nsy
        integer :: i,j,i1,j1,im1,jm1,m,m0,m1,m2,m3,m4,m5,ll
        real(8) :: dxp,dyp,dEE,c0,c,c1,eps,den,dxm,dym
        real(8)::  dpp,theta,dott,lp,lm,dEbx1,dEby1,dEbx2,dEby2

        real(8), dimension(0:N-1) :: ux,uy

        c0=sqrt(3.)/2.
        eps = 0.01
        ll=L*(L-1)

        i = mod(nsx,L)
        j = int(nsx/L)

        if ( j .eq. 0) then
            dEx = 0.d0
            dEy = 0.d0
            return
        endif


        m=i+j*L
        i1=i+1
        if (i1==L) i1=0
        j1=j+1
        if (j1==L) j1=0
        im1=i-1
        if (im1.lt.0) im1=L-1
        jm1=j-1
        if (jm1.lt.0) jm1=L-1


        m0=i1+j*L
        m1=i+j1*L
        m2=i1+j1*L
        m3=im1+j*L
        m4=i+jm1*L
        m5=im1+jm1*L


        dEx=0.
        dEy=0.
        lp=1.
           if (k(m,0).gt.eps) then
              dxp=v(m0)-v(m)
              dyp=v(m0+N)-v(m+N)
              dxp=dxp - L * dnint(dxp/L)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m,0)*(lp-1.)/lp
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
           endif

           lp=1.
           if (k(m,1).gt.eps) then
              dxp=v(m1)-v(m)
              dyp=v(m1+N)-v(m+N)
              dxp=dxp - L * dnint(dxp/L)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m,1)*(lp-1.)/lp
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
           endif

           lp=1.
           if (k(m,2).gt.eps) then
              dxp=v(m2)-v(m)
              dyp=v(m2+N)-v(m+N)
              dxp=dxp - L * dnint(dxp/L)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m,2)*(lp-1.)/lp
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
           endif

           lp=1.
           if (k(m3,0).gt.eps) then
              dxp=v(m3)-v(m)
              dyp=v(m3)-v(m+N)
              dxp=dxp - L * dnint(dxp/L)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m3,0)*(lp-1.)/lp
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
           endif

           lp=1.
           if (k(m4,1).gt.eps) then
              dxp=v(m4)-v(m)
              dyp=v(m4+N)-v(m+N)
              dxp=dxp - L * dnint(dxp/L)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m4,1)*(lp-1.)/(lp)
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
           endif

           lp=1.
           if (k(m5,2).gt.eps) then
              dxp=v(m5)-v(m)
              dyp=v(m5+N)-v(m+N)
              dxp=dxp - L * dnint(dxp/L)
              lp=sqrt(dxp*dxp+dyp*dyp)
              dEE=k(m5,2)*(lp-1.)/lp
              dEx=dEx-dxp*dEE
              dEy=dEy-dyp*dEE
           endif


        if (j .eq. (L-1)) then
            dEy = 0.d0
            return
        endif

      end subroutine singleDEnergy

      subroutine moveBoundary(v,L,strain,N)
          real(8), intent(inout), dimension(0:2*N-1) :: v
          integer, intent(in) ::L
          double precision ::strain
          integer ::N
          integer ::i,j,m
          double precision::c
          c=0.5d0*sqrt(3.)*(1.+strain)*(L-1)
          do i=0,L-1
          do j=0,1
             m=i+j*L*(L-1)
             v(m+N)=dble(j)*c
          enddo
          enddo
      end subroutine moveBoundary
