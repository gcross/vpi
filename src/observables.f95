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
!@+node:gcross.20090825141639.1529:accumulate_1d_position_averages
pure subroutine accumulate_1d_position_averages( &
  x, &
  n_particles, n_dimensions, &
  position_averages &
  )
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_particles, n_dimensions), intent(in) :: x
  double precision, dimension(n_dimensions), intent(inout) :: position_averages

  integer :: axis

  forall (axis = 1:n_dimensions) &
    position_averages(axis) = position_averages(axis) + sum(x(:,axis))/n_particles

end subroutine
!@-node:gcross.20090825141639.1529:accumulate_1d_position_averages
!@+node:gcross.20090825141639.1531:compute_radius_average
pure function compute_radius_average( &
  x, &
  n_particles, n_dimensions &
  ) result (radius_average)
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_particles, n_dimensions), intent(in) :: x

  integer :: i
  double precision :: radius_average

  radius_average = 0
  do i = 1, n_particles
    radius_average = radius_average + sqrt(sum(x(i,:)**2))
  end do
  radius_average = radius_average/n_particles

end function
!@-node:gcross.20090825141639.1531:compute_radius_average
!@+node:gcross.20090825141639.1533:compute_plane_radius_average
pure function compute_plane_radius_average( &
  x, &
  plane_axis_1, plane_axis_2, &
  n_particles, n_dimensions &
  ) result (radius_average)
  integer, intent(in) :: plane_axis_1, plane_axis_2
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_particles, n_dimensions), intent(in) :: x

  integer :: i
  double precision :: radius_average

  radius_average = 0
  do i = 1, n_particles
    radius_average = radius_average + sqrt(x(i,plane_axis_1)**2+x(i,plane_axis_2)**2)
  end do
  radius_average = radius_average/n_particles

end function
!@-node:gcross.20090825141639.1533:compute_plane_radius_average
!@+node:gcross.20090825141639.1535:compute_recip_plane_r_sq_average
pure function compute_recip_plane_r_sq_average( &
  x, &
  plane_axis_1, plane_axis_2, &
  n_particles, n_dimensions &
  ) result (radius_average)
  integer, intent(in) :: plane_axis_1, plane_axis_2
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_particles, n_dimensions), intent(in) :: x

  integer :: i
  double precision :: radius_average

  radius_average = 0
  do i = 1, n_particles
    radius_average = radius_average + 1.0/(x(i,plane_axis_1)**2+x(i,plane_axis_2)**2)
  end do
  radius_average = radius_average/n_particles

end function
!@-node:gcross.20090825141639.1535:compute_recip_plane_r_sq_average
!@-others

end module
!@nonl
!@-node:gcross.20090818081913.1354:@thin observables.f95
!@-leo
