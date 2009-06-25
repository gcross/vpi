!@+leo-ver=4-thin
!@+node:gcross.20090624094338.1396:@thin sp_HO_1D_potential.f90
!@@language fortran90
module sp_harmonic_os_1D_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624094338.1397:<< Imported modules >>
  use vpi_defines
  !@-node:gcross.20090624094338.1397:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624094338.1398:<< Variables >>
  real (kind=b8), private :: cylinder_radius
  !@-node:gcross.20090624094338.1398:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624094338.1399:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624094338.1400:init_sp
  subroutine init_sp ()
    namelist /single_particle_configuration/ cylinder_radius
    read(unit=10,nml=single_particle_configuration)
    write(*,*) "Using single particle box potential with"
    write(*,nml=single_particle_configuration)
  end subroutine init_sp
  !@-node:gcross.20090624094338.1400:init_sp
  !@+node:gcross.20090624094338.1401:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp,r

    Usp = 0
    r = sqrt( x(slice,ip,1)**2 + x(slice,ip,2)**2 ) 
    if(r > r_cylinder) then
      Usp = realbignumber*r**2
    end if

  end function Usp_func
  !@-node:gcross.20090624094338.1401:Usp
  !@+node:gcross.20090624094338.1402:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    integer :: slice, nslice, np, ndim
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    real(kind=b8), dimension ( np, ndim ) :: gUsp
    real(kind=b8), dimension ( np ) :: r

    gUsp = 0

  end function gUsp_func
  !@-node:gcross.20090624094338.1402:gUsp
  !@-others
  !@-node:gcross.20090624094338.1399:<< Subroutines >>
  !@nl

end module sp_harmonic_oscillator_1D_potential
!@-node:gcross.20090624094338.1396:@thin sp_HO_1D_potential.f90
!@-leo
