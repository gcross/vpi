!@+leo-ver=4-thin
!@+node:gcross.20090818081913.1354:@thin observables.f95
!@@language fortran90
module observables

implicit none

contains

!@+others
!@+node:gcross.20090818081913.1355:compute_local_energy_estimate
pure function compute_local_energy_estimate( &
  U, &
  gradient_of_log_trial_fn, laplacian_of_log_trial_fn, &
  lambda, &
  n_particles, n_dimensions &
  ) result ( energy )
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_particles), intent(in) :: U
  double precision, dimension(n_particles,n_dimensions), intent(in) :: gradient_of_log_trial_fn
  double precision, intent(in) :: laplacian_of_log_trial_fn, lambda
  double precision :: energy

  energy = sum(U) - lambda * ( sum(gradient_of_log_trial_fn(:,:)**2) + laplacian_of_log_trial_fn )
end function
!@-node:gcross.20090818081913.1355:compute_local_energy_estimate
!@-others

end module
!@nonl
!@-node:gcross.20090818081913.1354:@thin observables.f95
!@-leo
