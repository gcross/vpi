!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1981:@thin sp_lattice_trial.f90
!@@language fortran90
module sp_lattice_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1982:<< Imported modules >>
  use kinds
  use constants
  use sp_trial_numeric_differentiator
  !@-node:gcross.20090624144408.1982:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1983:<< Variables >>
  real(kind=b8), private :: x_harmonic_coefficient = 1_b8
  real(kind=b8), private :: y_harmonic_coefficient = 1_b8
  real(kind=b8), private :: z_harmonic_coefficient = 1_b8
  real(kind=b8), private :: lattice_ax = M_PI
  real(kind=b8), private :: lattice_ay = M_PI
  real(kind=b8), private :: lattice_az = M_PI
  real(kind=b8), private :: lattice_weight = 1
  real(kind=b8), private :: lattice_phase_x = M_PI/2.0_b8
  real(kind=b8), private :: lattice_phase_y = M_PI/2.0_b8
  real(kind=b8), private :: lattice_phase_z = M_PI/2.0_b8
  !@-node:gcross.20090624144408.1983:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1984:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1985:init_sp_tfunc
  subroutine init_sp_tfunc ()
    namelist /single_particle_trial_function_parameters/ &
      x_harmonic_coefficient, y_harmonic_coefficient, z_harmonic_coefficient, &
      lattice_ax, lattice_ay, lattice_az, &
      lattice_phase_x, lattice_phase_y, lattice_phase_z, &
      lattice_weight

    read(unit=10,nml=single_particle_trial_function_parameters)

    write(*,*) "Using anharmonic single particle trial function with"
    write(*,nml=single_particle_trial_function_parameters)
  end subroutine init_sp_tfunc
  !@-node:gcross.20090624144408.1985:init_sp_tfunc
  !@+node:gcross.20090624144408.1986:tfunc
  function tfunc( x, sl, nslice, np, ndim ) result( y )
    integer :: sl, nslice, np, ndim
    real(kind=b8), dimension( nslice, np , ndim ) :: x
    real(kind=b8) :: y
    real, dimension( np ) :: psi_t

    psi_t(:)  = -(  x_harmonic_coefficient*x(sl,:,1)**2 &
                  + y_harmonic_coefficient*x(sl,:,2)**2 &
                  + z_harmonic_coefficient*x(sl,:,3)**2 &
                 )/2.0 + &
                lattice_weight*sin(x(sl,:,1)*lattice_ax+lattice_phase_x)**2 + &
                lattice_weight*sin(x(sl,:,2)*lattice_ay+lattice_phase_y)**2 + &
                lattice_weight*sin(x(sl,:,3)*lattice_az+lattice_phase_z)**2
    y = sum(psi_t)

  end function tfunc
  !@-node:gcross.20090624144408.1986:tfunc
  !@+node:gcross.20090624144408.1987:grad & lapacian of tfunc
  function grad_lap_sp_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
    implicit none

    real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, slice
    integer :: y

    y = numeric_grad_lap_spf( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn, tfunc )

  end function grad_lap_sp_tfun
  !@-node:gcross.20090624144408.1987:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1984:<< Subroutines >>
  !@nl

end module sp_lattice_trial
!@-node:gcross.20090624144408.1981:@thin sp_lattice_trial.f90
!@-leo
