!@+leo-ver=4-thin
!@+node:gcross.20090624094338.1396:@thin sp_harmonic_oscillator_1D_potential.f90
!@@language fortran90
module sp_harmonic_oscillator_1D_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624094338.1397:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624094338.1397:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624094338.1398:<< Variables >>
  real (kind=b8), private :: coefficient = 1.0_b8
  !@-node:gcross.20090624094338.1398:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624094338.1399:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624094338.1400:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ coefficient

    read(unit=10,nml=single_particle_potential_parameters)

    write(*,*) "Using single particle harmonic oscillator 1D potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624094338.1400:init_sp_potential
  !@+node:gcross.20090624094338.1401:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp

    Usp = coefficient*x(slice,ip,3)**2 

  end function Usp_func
  !@-node:gcross.20090624094338.1401:Usp
  !@+node:gcross.20090624094338.1402:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    gUsp(:,1) = 0.0_b8
    gUsp(:,2) = 0.0_b8
    gUsp(:,3) = coefficient*x(slice,:,3)

  end function gUsp_func
  !@-node:gcross.20090624094338.1402:gUsp
  !@-others
  !@-node:gcross.20090624094338.1399:<< Subroutines >>
  !@nl

end module sp_harmonic_oscillator_1D_potential
!@-node:gcross.20090624094338.1396:@thin sp_harmonic_oscillator_1D_potential.f90
!@-leo