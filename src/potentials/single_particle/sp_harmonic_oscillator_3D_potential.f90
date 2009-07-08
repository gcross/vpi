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
  real (kind=b8), dimension(3), private :: coefficients = (/1.0_b8,1.0_b8,1.0_b8/)

  !@-node:gcross.20090624094338.1412:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624094338.1413:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624094338.1414:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ coefficients

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

    Usp = dot_product(coefficients,x(slice,ip,:)**2) / 2.0_b8

  end function Usp_func
  !@-node:gcross.20090624094338.1415:Usp
  !@+node:gcross.20090624094338.1416:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim
    integer :: i

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    do i=1,3
      gUsp(:,i) = coefficients(i)*x(slice,:,i)
    end do

  end function gUsp_func
  !@-node:gcross.20090624094338.1416:gUsp
  !@-others
  !@-node:gcross.20090624094338.1413:<< Subroutines >>
  !@nl

end module vpi_single_particle_potential
!@-node:gcross.20090624094338.1410:@thin sp_harmonic_oscillator_3D_potential.f90
!@-leo
