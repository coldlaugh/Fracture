      subroutine f_ortn(v,k,N,L,E)
C      spring energy, spring const k,  and bending energy, const kappa
      implicit none
      integer ::N,L
      double precision, dimension(0:2*N-1) :: v
      double precision, dimension(0:N-1,0:2) ::k
	  double precision :: E,NE
      integer :: i,j,i1,j1,m,m0,m1,m2,ll,m3,m4,m5,im1,jm1
      double precision :: l0,l1,l2,dx,dy,eps
      double precision, dimension(0:N-1) :: ux,uy

      eps=0.00001
	  E = 0
	  
      do m=0,N-1
      ux(m)=v(m)
      uy(m)=v(m+N)
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
           E=E+dabs(l0-1.)*(2.d0*dy*dy/l0/l0 -1.d0)    
		   NE = NE + dabs(l0-1.)     
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
           E=E+dabs(l1-1.)*(2.d0*dy*dy/l1/l1-1.d0)
		   NE = NE + dabs(l1-1.)
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
           E=E+dabs(l2-1.)*(2.d0*dy*dy/l2/l2-1.d0)
		   NE = NE + dabs(l2-1.)
        endif
      enddo
      enddo
	  E = E/NE
      end subroutine f_ortn
