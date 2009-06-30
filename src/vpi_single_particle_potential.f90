!@+leo-ver=4-thin
!@+node:gcross.20090629153134.1738:@thin vpi_single_particle_potential.f90
!@@language fortran90

module vpi_single_particle_potential

use kinds

implicit none

contains

subroutine init_sp_potential ()
end subroutine init_sp_potential

function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
  integer :: nslice, np, ndim
  integer :: slice, ip
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8) :: Usp
end function Usp_func

function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
  real(kind=b8), dimension ( nslice, np , ndim ) ::  x
  integer :: slice, nslice, np, ndim
  real(kind=b8), dimension ( np, ndim ) :: gUsp
end function gUsp_func

end module vpi_single_particle_potential
!@-node:gcross.20090629153134.1738:@thin vpi_single_particle_potential.f90
!@-leo
