!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1411:@thin sp_null_potential.f90
!@@language fortran90
module vpi_single_particle_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1412:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1412:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1413:<< Variables >>
  !@-node:gcross.20090624144408.1413:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1414:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1415:init_sp_potential
  subroutine init_sp_potential ()
    write(*,*) "Using no single particular potential."
  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624144408.1415:init_sp_potential
  !@+node:gcross.20090624144408.1416:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp

    Usp = 0

  end function Usp_func
  !@-node:gcross.20090624144408.1416:Usp
  !@+node:gcross.20090624144408.1417:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    gUsp = 0

  end function gUsp_func

  !@-node:gcross.20090624144408.1417:gUsp
  !@-others
  !@-node:gcross.20090624144408.1414:<< Subroutines >>
  !@nl

end module vpi_single_particle_potential
!@-node:gcross.20090624144408.1411:@thin sp_null_potential.f90
!@-leo
