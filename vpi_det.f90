module vpi_potential
  use vpi_defines
  use vpi_aziz
  implicit none

contains 

function vpi_det( x, slice, ip, nslice, np, ndim ) result ( det )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( np , np ) :: mat

  real(kind=b8) :: Usp

  do i = 1, np
    do j = 1, np-1, 2
      mat(i,j) = cos(M_PI_2*j*x(slice,ip,1))
      mat(i,j+1) = sin(M_PI_2*j*x(slice,ip,1))
    end do
  end do



end function vpi_Usp_NULL

end module vpi_potential
