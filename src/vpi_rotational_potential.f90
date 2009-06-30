!@+leo-ver=4-thin
!@+node:gcross.20090629153134.1743:@thin vpi_rotational_potential.f90
!@@language fortran90

module vpi_rotational_potential

use kinds
use vpi_defines

implicit none

contains

subroutine init_rot_potential ()
end subroutine init_rot_potential

function Uij_rot_func( x, x_rot, xij2, slice, ip, nslice, np, ndim ) result ( Uij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, ip, nslice, np, ndim
  real(kind=b8) :: Uij
end function Uij_rot_func

function gUij_rot_func( x, x_rot, xij2, slice, nslice, np, ndim ) result ( gUij )
  implicit none
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, nslice, np, ndim
  real(kind = b8), dimension ( np , ndim ) :: gUij
end function gUij_rot_func

end module

!@-node:gcross.20090629153134.1743:@thin vpi_rotational_potential.f90
!@-leo
