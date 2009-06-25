!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1485:@thin sp_N_well_potential.f90
!@@language fortran90
module sp_N_well_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1486:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.1486:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1487:<< Variables >>
  real(kind=b8), private :: lam_ho = 1.0_b8
  real(kind=b8), private :: p_nw_vb = 10.0_b8
  real(kind=b8), private :: p_nw_l = M_PI
  !@-node:gcross.20090624144408.1487:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1488:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1489:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ lam_ho, p_nw_vb, p_nw_l

    read(unit=10,nml=single_particle_potential_parameters)

    write(*,*) "Using single particle N well potential with"
    write(*,nml=single_particle_potential_parameters)

  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624144408.1489:init_sp_potential
  !@+node:gcross.20090624144408.1490:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp

    Usp = ( x(slice,ip,1)**2 + x(slice,ip,2)**2 + lam_ho*x(slice,ip,3)**2) + p_nw_vb*cos(x(slice,ip,3)*p_nw_l)

  end function Usp_func
  !@-node:gcross.20090624144408.1490:Usp
  !@+node:gcross.20090624144408.1491:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    gUsp(:,1) = 2.0_b8*x(slice,:,1)
    gUsp(:,2) = 2.0_b8*x(slice,:,2)
    gUsp(:,3) = 2.0_b8*lam_ho*x(slice,:,3) - p_nw_vb*p_nw_l*sin(x(slice,:,3)*p_nw_l)

  end function gUsp_func
  !@-node:gcross.20090624144408.1491:gUsp
  !@-others
  !@-node:gcross.20090624144408.1488:<< Subroutines >>
  !@nl

end module sp_N_well_potential
!@-node:gcross.20090624144408.1485:@thin sp_N_well_potential.f90
!@-leo
