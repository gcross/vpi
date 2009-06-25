!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1601:@thin sp_lattice_3D_with_HO_potential.f90
!@@language fortran90
module sp_lattice_3D_with_HO_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1602:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.1602:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1603:<< Variables >>
  real(kind=b8), private :: p_lattice_vb = 80
  real(kind=b8), private :: lam_ho = 1
  real(kind=b8), private :: p_lattice_ax = M_PI
  real(kind=b8), private :: p_lattice_ay = M_PI
  real(kind=b8), private :: p_lattice_az = M_PI
  real(kind=b8), private :: p_lattice_phase_x = M_PI/2.0_b8
  real(kind=b8), private :: p_lattice_phase_y = M_PI/2.0_b8
  real(kind=b8), private :: p_lattice_phase_z = M_PI/2.0_b8
  !@-node:gcross.20090624144408.1603:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1604:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1605:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ &
      p_lattice_vb, lam_ho, &
      p_lattice_ax, p_lattice_ay, p_lattice_az, &
      p_lattice_phase_x, p_lattice_phase_y, p_lattice_phase_z

    read(unit=10,nml=single_particle_potential_parameters)

    write(*,*) "Using 2D lattice potential plus a harmonic oscillator potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624144408.1605:init_sp_potential
  !@+node:gcross.20090624144408.1606:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp

    Usp = 0.0_b8

    Usp = 0.5_b8*( x(slice,ip,1)**2 + x(slice,ip,2)**2 + lam_ho*x(slice,ip,3)**2) +  &
       p_lattice_vb*cos(x(slice,ip,1)*p_lattice_ax+p_lattice_phase_x)**2 + &
       p_lattice_vb*cos(x(slice,ip,2)*p_lattice_ay+p_lattice_phase_y)**2 + & 
       p_lattice_vb*cos(x(slice,ip,3)*p_lattice_az+p_lattice_phase_z)**2

  end function Usp_func
  !@nonl
  !@-node:gcross.20090624144408.1606:Usp
  !@+node:gcross.20090624144408.1607:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    gUsp(:,1) = x(slice,:,1) - &
       2.0_b8*p_lattice_vb*p_lattice_ax*cos(x(slice,:,1)*p_lattice_ax + p_lattice_phase_x)* &
       sin(x(slice,:,1)*p_lattice_ax+p_lattice_phase_x)
    gUsp(:,2) = lam_ho*x(slice,:,2) - &
       2.0_b8*p_lattice_vb*p_lattice_ay*cos(x(slice,:,2)*p_lattice_ay+p_lattice_phase_y)* &
       sin(x(slice,:,2)*p_lattice_ay+p_lattice_phase_y)
    gUsp(:,3) = x(slice,:,3) - &
       2.0_b8*p_lattice_vb*p_lattice_az*cos(x(slice,:,3)*p_lattice_az+p_lattice_phase_z)* &
       sin(x(slice,:,3)*p_lattice_az+p_lattice_phase_z)

  end function gUsp_func
  !@-node:gcross.20090624144408.1607:gUsp
  !@-others
  !@-node:gcross.20090624144408.1604:<< Subroutines >>
  !@nl

end module sp_lattice_3D_with_HO_potential
!@-node:gcross.20090624144408.1601:@thin sp_lattice_3D_with_HO_potential.f90
!@-leo
