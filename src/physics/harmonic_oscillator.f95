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
  double precision, dimension( n_particles, n_dimensions ), intent(in) :: x
  double precision, dimension( n_dimensions ), intent(in) :: trial_coefficients
  double precision  :: weight

  double precision, dimension( n_dimensions ) :: temp
  integer :: i

  forall (i=1:n_dimensions) &
    temp(i) = sum(x(:,i)**2)

  weight = -dot_product(temp,trial_coefficients)/2d0
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
  double precision, dimension( n_particles, n_dimensions ), intent(in) :: x
  double precision, dimension( n_dimensions ), intent(in) :: trial_coefficients
  double precision, dimension( n_particles, n_dimensions ), intent(inout) :: gradient_of_log_trial_fn
  double precision, intent(inout) :: laplacian_of_log_trial_fn

  integer :: i

  forall (i=1:n_dimensions) &
    gradient_of_log_trial_fn(:,i) = gradient_of_log_trial_fn(:,i) - trial_coefficients(i)*x(:,i)

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
  double precision, dimension( n_slices, n_particles, n_dimensions ), intent(in) :: x
  double precision, dimension( n_dimensions ), intent(in) :: potential_coefficients
  double precision, intent(inout) :: &
    U(n_slices,n_particles), &
    gradU2(n_slices)

  integer :: i,j

  forall (i=1:n_slices,j=1:n_particles) &
    U(i,j) = U(i,j) + dot_product(potential_coefficients,x(i,j,:)**2)/2d0

  do i = 1, n_slices
    do j = 1, n_particles
      gradU2(i) = gradU2(i) + dot_product(potential_coefficients,x(i,j,:))**2
    end do
  end do

end subroutine
!@-node:gcross.20090901084550.2626:accumulate_potential
!@-others

end module
!@-node:gcross.20090901084550.2618:@thin harmonic_oscillator.f95
!@-leo
