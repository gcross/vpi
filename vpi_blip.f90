module vpi_blip
  use vpi_defines
  implicit none

contains 

subroutine  blip(x,x0,b0,y,dy,ddy)
  real(kind=b8), intent(in) :: x,x0,b0
  real(kind=b8), intent(out) :: y,dy,ddy

  real(kind=b8) :: b2i,xx,ax,ax2
  integer :: sgnax

  b2i = 2.0_b8/b0
  b2i2 = b2i**2

  xx = 2.0*(x-x0)/b0
  ax = abs(xx)
  ax2 = ax**2
  sgnax = int(xx/ax)

  y = 0.0_b8
  dy = 0.0_b8
  ddy = 0.0_b8
  if (ax < 1) then
    y = 1.0_b8 - 1.5_b8*ax2 + 0.75_b8*ax2*ax
    dy = b2i*sgnax*(-3.0_b8*ax + 2.25_b8*ax2)
    ddy = b2i2*(-3.0_b8 + 4.5_b8*ax)
  else
    if (ax < 2) then
      tmp1 = (2.0_b8_-ax)
      tmp2 = tmp1**2
      tmp3 = tmp2*tmp3
      y = 0.25*tmp3
      dy = b2i*sgnax*-0.75_b8*tmp2
      ddy = b2i2*1.5_b8*tmp1
    end if
  end if
end  subroutine blip


end module vpi_blip
