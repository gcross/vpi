!@+leo-ver=4-thin
!@+node:gcross.20090901084550.2618:@thin harmonic_oscillator.f95
!@@language fortran90

module harmonic_oscillator
  implicit none

contains

!@+others
!@+node:gcross.20090901084550.2619:compute_trial_weight
pure function compute_trial_weight( &
    x, &
    trial_coefficients, &
    n_particles, n_dimensions &
  ) result( weight )
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_dimensions,n_particles), intent(in) :: x
  double precision, dimension(n_dimensions), intent(in) :: trial_coefficients
  double precision  :: weight

  weight = -dot_product(sum(x*x,2),trial_coefficients)/2d0
end function
!@-node:gcross.20090901084550.2619:compute_trial_weight
!@+node:gcross.20090901084550.2621:accumulate_trial_derivatives
pure subroutine accumulate_trial_derivatives( &
    x, &
    trial_coefficients, &
    n_particles, n_dimensions, &
    gradient_of_log_trial_fn, laplacian_of_log_trial_fn &
  )
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_dimensions,n_particles), intent(in) :: x
  double precision, dimension(n_dimensions), intent(in) :: trial_coefficients
  double precision, dimension(n_dimensions,n_particles), intent(inout) :: gradient_of_log_trial_fn
  double precision, intent(inout) :: laplacian_of_log_trial_fn

  integer :: particle

  do particle = 1, n_particles
    gradient_of_log_trial_fn(:,particle) = gradient_of_log_trial_fn(:,particle) - trial_coefficients(:)*x(:,particle)
  end do

  laplacian_of_log_trial_fn = laplacian_of_log_trial_fn - sum(trial_coefficients)*n_particles
end subroutine
!@-node:gcross.20090901084550.2621:accumulate_trial_derivatives
!@+node:gcross.20090901084550.2626:accumulate_potential
subroutine accumulate_potential( &
    x, &
    potential_coefficients, &
    n_slices, n_particles, n_dimensions, &
    U, &
    gradU2 &
  )
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension(n_dimensions,n_particles,n_slices), intent(in) :: x
  double precision, dimension(n_dimensions), intent(in) :: potential_coefficients
  double precision, intent(inout) :: &
    U(n_particles,n_slices), &
    gradU2(n_slices)

  integer :: slice, particle

  do slice = 1, n_slices
    do particle = 1, n_particles
      U(particle,slice) = U(particle,slice) + dot_product(potential_coefficients,x(:,particle,slice)**2)/2d0
      gradU2(slice) = gradU2(slice) + dot_product(potential_coefficients,x(:,particle,slice))**2
    end do
  end do

end subroutine
!@-node:gcross.20090901084550.2626:accumulate_potential
!@-others

end module
!@-node:gcross.20090901084550.2618:@thin harmonic_oscillator.f95
!@-leo
