!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1426:@thin sp_wlink_potential.f90
!@@language fortran90
module sp_wlink_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1427:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1427:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1428:<< Variables >>
  real (kind=b8), private :: a1wlink, a2wlink, vbwlink, w0wlink
  !@-node:gcross.20090624144408.1428:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1429:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1430:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ a1wlink, a2wlink, vbwlink, w0wlink

    read(unit=10,nml=single_particle_potential_parameters)

    stop "The wlink potential has not been fully implemented yet."

    write(*,*) "Using single particle 'wlink' potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624144408.1430:init_sp_potential
  !@+node:gcross.20090624144408.1431:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8) :: Usp

    real(kind=b8) :: term1, term2, r2

    stop "The wlink potential has not been fully implemented yet."

    r2 = sum(x(slice,ip,:)**2)
    term1 = exp(-(x(slice,ip,3)/a1wlink)**12)
    term2 = exp(-((x(slice,ip,1)**2+x(slice,ip,2)**2)/a2wlink**2)**6)
    Usp = vbwlink*term1*(1.0_b8-term2) + w0wlink*r2
  end function Usp_func
  !@-node:gcross.20090624144408.1431:Usp
  !@+node:gcross.20090624144408.1432:gUsp
  function gUsp_func( x, ip, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: ip, slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp
    real(kind=b8) :: term1, term2, r2

    r2 = sum(x(slice,ip,:)**2)
    term1 = exp(-(x(slice,ip,3)/a1wlink)**12)
    term2 = exp(-((x(slice,ip,1)**2+x(slice,ip,2)**2)/a2wlink**2)**6)

    gUsp = 0

    stop "The wlink potential has not been fully implemented yet."

  !@+at
  !   gUsp(:,1) = 2.0_b8*x(slice,:,1)
  !   gUsp(:,2) = 2.0_b8*x(slice,:,2)
  !   gUsp(:,3) = 2.0_b8*lam_ho*x(slice,:,3)
  !@-at
  !@@c
  end function gUsp_func
  !@-node:gcross.20090624144408.1432:gUsp
  !@-others
  !@-node:gcross.20090624144408.1429:<< Subroutines >>
  !@nl

end module sp_wlink_potential
!@-node:gcross.20090624144408.1426:@thin sp_wlink_potential.f90
!@-leo
