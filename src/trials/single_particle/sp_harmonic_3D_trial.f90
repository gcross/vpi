!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1852:@thin sp_harmonic_3D_trial.f90
!@@language fortran90
module vpi_single_particle_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1853:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1853:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1854:<< Variables >>
  real (kind=b8), private :: x_coefficient = 1.0_b8
  real (kind=b8), private :: y_coefficient = 1.0_b8
  real (kind=b8), private :: z_coefficient = 1.0_b8
  !@-node:gcross.20090624144408.1854:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1855:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1856:init_sp_tfunc
  subroutine init_sp_tfunc ()
    namelist /single_particle_trial_function_parameters/ x_coefficient, y_coefficient, z_coefficient

    read(unit=10,nml=single_particle_trial_function_parameters)

    write(*,*) "Using harmonic trial function with"
    write(*,nml=single_particle_trial_function_parameters)
  end subroutine init_sp_tfunc
  !@-node:gcross.20090624144408.1856:init_sp_tfunc
  !@+node:gcross.20090624144408.1857:tfunc
  function tfunc( x, sl, nslice, np, ndim ) result( y )
    integer :: sl, nslice, np, ndim
    real(kind=b8), dimension( nslice, np , ndim ) :: x
    real(kind=b8) :: y
    real, dimension( np ) :: psi_t

    psi_t(:)  = -( x_coefficient*x(sl,:,1)**2 &
                 + y_coefficient*x(sl,:,2)**2 &
                 + z_coefficient*x(sl,:,3)**2 &
                 )/2.0_b8
    y = sum(psi_t)

  end function tfunc
  !@-node:gcross.20090624144408.1857:tfunc
  !@+node:gcross.20090624144408.1858:grad & lapacian of tfunc
  function grad_lap_sp_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
    real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, slice
    integer :: y

    grad_lntfn(:,1) = -x_coefficient*x(slice,:,1)
    grad_lntfn(:,2) = -y_coefficient*x(slice,:,2)
    grad_lntfn(:,3) = -z_coefficient*x(slice,:,3)

    lap_lntfn = -(x_coefficient + y_coefficient + z_coefficient)*dble(np) 
    y = 1

  end function grad_lap_sp_tfun
  !@-node:gcross.20090624144408.1858:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1855:<< Subroutines >>
  !@nl

end module vpi_single_particle_trial
!@-node:gcross.20090624144408.1852:@thin sp_harmonic_3D_trial.f90
!@-leo
