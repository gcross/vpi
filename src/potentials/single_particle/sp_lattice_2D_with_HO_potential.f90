!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1587:@thin sp_lattice_2D_with_HO_potential.f90
!@@language fortran90
module vpi_single_particle_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1588:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.1588:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1589:<< Variables >>
  real(kind=b8), private :: p_lattice_vb = 80
  real(kind=b8), private :: lam_ho = 1
  real(kind=b8), private :: p_lattice_ax = M_PI
  real(kind=b8), private :: p_lattice_az = M_PI
  !@-node:gcross.20090624144408.1589:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1590:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1591:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ p_lattice_vb, p_lattice_ax, p_lattice_az, lam_ho

    read(unit=10,nml=single_particle_potential_parameters)

    write(*,*) "Using 2D lattice potential plus a harmonic oscillator potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624144408.1591:init_sp_potential
  !@+node:gcross.20090624144408.1592:Usp
  function vpi_Usp_lattice( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp

    Usp = 0.0_b8

    Usp = ( x(slice,ip,1)**2 + lam_ho*x(slice,ip,2)**2 + x(slice,ip,3)**2) &
         + p_lattice_vb*cos(x(slice,ip,1)*p_lattice_ax)**2 &
         + p_lattice_vb*cos(x(slice,ip,3)*p_lattice_az)**2

  end function vpi_Usp_lattice
  !@-node:gcross.20090624144408.1592:Usp
  !@+node:gcross.20090624144408.1593:gUsp
  function vpi_gUsp_lattice( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    gUsp(:,1) = 2.0_b8*x(slice,:,1) - 2.0_b8*p_lattice_vb*p_lattice_ax*cos(x(slice,:,1)*p_lattice_ax)*sin(x(slice,:,1)*p_lattice_ax)
    gUsp(:,2) = 2.0_b8*lam_ho*x(slice,:,2)
    gUsp(:,3) = 2.0_b8*x(slice,:,3) - 2.0_b8*p_lattice_vb*p_lattice_ax*cos(x(slice,:,3)*p_lattice_az)*sin(x(slice,:,3)*p_lattice_az)

  end function vpi_gUsp_lattice
  !@-node:gcross.20090624144408.1593:gUsp
  !@-others
  !@-node:gcross.20090624144408.1590:<< Subroutines >>
  !@nl

end module vpi_single_particle_potential
!@-node:gcross.20090624144408.1587:@thin sp_lattice_2D_with_HO_potential.f90
!@-leo
