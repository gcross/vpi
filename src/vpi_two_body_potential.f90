!@+leo-ver=4-thin
!@+node:gcross.20090629153134.1740:@thin vpi_two_body_potential.f90
!@@language fortran90

module vpi_two_body_potential

use kinds

implicit none

contains

subroutine init_tb_potential ()
end subroutine init_tb_potential

function Uij_func( x, xij2, slice, ip, nslice, np, ndim, reject_flag ) result ( Uij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, ip, nslice, np, ndim
  logical :: reject_flag
  real(kind=b8) :: Uij
end function Uij_func

function gUij_func( x, xij2, slice, nslice, np, ndim ) result ( gUij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, nslice, np, ndim
  real(kind = b8), dimension ( np , ndim ) :: gUij
end function gUij_func

end module vpi_two_body_potential
!@-node:gcross.20090629153134.1740:@thin vpi_two_body_potential.f90
!@-leo
