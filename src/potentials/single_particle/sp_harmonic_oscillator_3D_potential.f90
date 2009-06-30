!@+leo-ver=4-thin
!@+node:gcross.20090624094338.1410:@thin sp_harmonic_oscillator_3D_potential.f90
!@@language fortran90
module vpi_single_particle_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624094338.1411:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624094338.1411:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624094338.1412:<< Variables >>
  real (kind=b8), private :: x_coefficient = 1.0_b8
  real (kind=b8), private :: y_coefficient = 1.0_b8
  real (kind=b8), private :: z_coefficient = 1.0_b8
  !@-node:gcross.20090624094338.1412:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624094338.1413:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624094338.1414:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ x_coefficient, y_coefficient, z_coefficient

    read(unit=10,nml=single_particle_potential_parameters)

    write(*,*) "Using single particle harmonic oscillator 3D potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@-node:gcross.20090624094338.1414:init_sp_potential
  !@+node:gcross.20090624094338.1415:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp

    Usp = ( x_coefficient*x(slice,ip,1)**2 &
          + y_coefficient*x(slice,ip,2)**2 &
          + z_coefficient*x(slice,ip,3)**2 &
          ) / 2.0_b8

  end function Usp_func
  !@-node:gcross.20090624094338.1415:Usp
  !@+node:gcross.20090624094338.1416:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    gUsp(:,1) = x_coefficient*x(slice,:,1)
    gUsp(:,2) = y_coefficient*x(slice,:,2)
    gUsp(:,3) = z_coefficient*x(slice,:,3)

  end function gUsp_func
  !@nonl
  !@-node:gcross.20090624094338.1416:gUsp
  !@-others
  !@-node:gcross.20090624094338.1413:<< Subroutines >>
  !@nl

end module vpi_single_particle_potential
!@-node:gcross.20090624094338.1410:@thin sp_harmonic_oscillator_3D_potential.f90
!@-leo
