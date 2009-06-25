!@+leo-ver=4-thin
!@+node:gcross.20090624094338.1424:@thin sp_harmonic_oscillator_overlapping_potential.f90
!@@language fortran90
module sp_harmonic_oscillator_overlapping_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624094338.1425:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624094338.1425:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624094338.1426:<< Variables >>
  real(kind=b8), private :: domega = 1.0001_b8
  real(kind=b8), private :: lam_ho = 1.0_b8
  !@-node:gcross.20090624094338.1426:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624094338.1427:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624094338.1428:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ domega, lam_ho

    read(unit=10,nml=single_particle_potential_parameters)

    write(*,*) "Using single particle harmonic oscillator overlapping potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624094338.1428:init_sp_potential
  !@+node:gcross.20090624094338.1429:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp

    if(slice .gt. CSLICE) then
      Usp = domega*( x(slice,ip,1)**2 + x(slice,ip,2)**2 + lam_ho*x(slice,ip,3)**2)
    else
      Usp = ( x(slice,ip,1)**2 + x(slice,ip,2)**2 + lam_ho*x(slice,ip,3)**2)
    end if

  end function Usp_func
  !@-node:gcross.20090624094338.1429:Usp
  !@+node:gcross.20090624094338.1430:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    if(slice .gt. CSLICE) then
      gUsp(:,1) = 2.0_b8*domega*x(slice,:,1)
      gUsp(:,2) = 2.0_b8*domega*x(slice,:,2)
      gUsp(:,3) = 2.0_b8*domega*lam_ho*x(slice,:,3)
    else
      gUsp(:,1) = 2.0_b8*x(slice,:,1)
      gUsp(:,2) = 2.0_b8*x(slice,:,2)
      gUsp(:,3) = 2.0_b8*lam_ho*x(slice,:,3)
    end if

  end function gUsp_func
  !@-node:gcross.20090624094338.1430:gUsp
  !@-others
  !@-node:gcross.20090624094338.1427:<< Subroutines >>
  !@nl

end module sp_harmonic_oscillator_overlapping_potential
!@-node:gcross.20090624094338.1424:@thin sp_harmonic_oscillator_overlapping_potential.f90
!@-leo
