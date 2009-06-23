module spline_util

  use vpi_defines

  implicit none
  real (kind=b8) pi=acos(-1.0d0)

contains

subroutine spline(x,y,n,yp1,ypn,y2)
  integer, parameter :: nmax=10000
  integer :: n
  real (kind=b8) :: yp1, ypn
  real (kind=b8), dimension(n) :: x,y,y2
  real (kind=b8), dimension(nmax) u

  if (yp1.gt..99d30) then
    y2a(1)=0.d0
    u(1)=0.d0
  else
    y2a(1)=-0.5d0
    u(1)=(3.d0/(xa(2)-xa(1)))*((ya(2)-ya(1))/(xa(2)-xa(1))-yp1)
  endif

  do i=2,n-1
    sig=(xa(i)-xa(i-1))/(xa(i+1)-xa(i-1))
    p=sig*y2a(i-1)+2.d0
    y2a(i)=(sig-1.d0)/p
    u(i)=(6.d0*((ya(i+1)-ya(i))/(xa(i+1)-xa(i))-(ya(i)-ya(i-1)) &
         (xa(i)-xa(i-1)))/(xa(i+1)-xa(i-1))-sig*u(i-1))/p
  end do

  if (ypn.gt..99d30) then
    qn=0.d0
    un=0.d0
  else
    qn=0.5d0
    un=(3.d0/(xa(n)-xa(n-1)))*(ypn-(ya(n)-ya(n-1))/(xa(n)-xa(n-1)))
  endif

  y2a(n)=(un-qn*u(n-1))/(qn*y2a(n-1)+1.)
  do k=n-1,1,-1
    y2a(k)=y2a(k)*y2a(k+1)+u(k)
  end do
end subroutine spline

subroutine splint(xa,ya,y2a,n,rjac,y,maxy,ny,nlam)
  integer :: n,ny,nlam
  real (kind=b8), dimension(n) :: xa,ya 
  real (kind=b8), dimension(n,nlam) :: xa,ya,y2a
  real (kind=b8), dimension(maxy,nlam) :: y
  real (kind=b8), dimension(ny) :: rjac

  do i=1,ny
    x=rjac(i)

    h=(xa(2)-xa(1))
    hinv=1.d0/h
    klo=(x-xa(1))*hinv+1
    khi=klo+1

    a=(xa(khi)-x)*hinv
    b=(x-xa(klo))*hinv

    do j=1,nlam
      y(i,j)=a*ya(klo,j)+b*ya(khi,j)+
    & ((a**3-a)*y2a(klo,j)+(b**3-b)*y2a(khi,j))*(h**2)/6.d0
    end do
  end do
end subroutine splint

end module spline_util

