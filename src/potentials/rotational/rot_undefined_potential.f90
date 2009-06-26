!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1772:@thin rot_undefined_potential.f90
!@@language fortran90
module rot_undefined_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1773:<< Imported modules >>
  use kinds
  use vpi_defines
  !@-node:gcross.20090624144408.1773:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1774:<< Variables >>
  real (kind=b8), private :: coefficient = 1e-4_b8

  !@-node:gcross.20090624144408.1774:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1775:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1776:init_rot_potential
  subroutine init_rot_potential ()
    stop "The rotational potential is undefined and should not be initialized."
  end subroutine init_rot_potential
  !@-node:gcross.20090624144408.1776:init_rot_potential
  !@+node:gcross.20090624144408.1777:RUij
  function Uij_rot_func( x, x_rot, xij2, slice, ip, nslice, np, ndim ) result ( Uij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, ip, nslice, np, ndim

    real(kind=b8) :: Uij

    stop "The rotational potential is undefined and should not be evaluated."

  end function Uij_rot_func
  !@-node:gcross.20090624144408.1777:RUij
  !@+node:gcross.20090624144408.1778:gRUij
  function gUij_rot_func( x, x_rot, xij2, slice, nslice, np, ndim ) result ( gUij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, nslice, np, ndim

    real(kind = b8), dimension ( np , ndim ) :: gUij

    stop "The rotational potential is undefined and so its derivative should not be evaluated."

  end function gUij_rot_func
  !@-node:gcross.20090624144408.1778:gRUij
  !@-others
  !@-node:gcross.20090624144408.1775:<< Subroutines >>
  !@nl

end module rot_undefined_potential
!@-node:gcross.20090624144408.1772:@thin rot_undefined_potential.f90
!@-leo
