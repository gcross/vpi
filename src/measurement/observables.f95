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
  double precision, dimension(n_dimensions,n_particles), intent(in) :: gradient_of_log_trial_fn
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
  double precision, dimension(n_dimensions,n_particles), intent(in) :: x
  double precision, dimension(n_dimensions), intent(inout) :: position_averages

  double precision :: current_average(n_dimensions)
  integer :: particle

  current_average = 0
  do particle = 1, n_particles
    current_average(:) = current_average(:) + x(:,particle)
  end do
  position_averages = current_average/n_particles

end subroutine
!@-node:gcross.20090825141639.1529:accumulate_1d_position_averages
!@+node:gcross.20090825141639.1531:compute_radius_average
pure function compute_radius_average( &
  x, &
  n_particles, n_dimensions &
  ) result (radius_average)
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_dimensions,n_particles), intent(in) :: x

  integer :: particle
  double precision :: radius_average

  radius_average = 0
  do particle = 1, n_particles
    radius_average = radius_average + sqrt(sum(x(:,particle)**2))
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
  double precision, dimension(n_dimensions,n_particles), intent(in) :: x

  integer :: particle
  double precision :: radius_average

  radius_average = 0
  do particle = 1, n_particles
    radius_average = radius_average + sqrt(x(plane_axis_1,particle)**2+x(plane_axis_2,particle)**2)
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

  integer :: particle
  double precision :: radius_average

  radius_average = 0
  do particle = 1, n_particles
    radius_average = radius_average + 1.0/(x(plane_axis_1,particle)**2+x(plane_axis_2,particle)**2)
  end do
  radius_average = radius_average/n_particles

end function
!@-node:gcross.20090825141639.1535:compute_recip_plane_r_sq_average
!@+node:gcross.20090826091347.1420:compute_average_angular_separation
pure function compute_average_angular_separation(angles,n_particles) result (average)
  integer, intent(in) :: n_particles
  double precision, dimension(n_particles), intent(in) :: angles

  integer :: particle_1, particle_2
  double precision :: average

  average = 0
  do particle_1 = 1, n_particles
    do particle_2 = particle_1+1, n_particles
      average = average + abs(angles(particle_1)-angles(particle_2))
    end do
  end do
  average = average / dble(n_particles * (n_particles-1) / 2)

end function
!@-node:gcross.20090826091347.1420:compute_average_angular_separation
!@+node:gcross.20090826091347.1422:compute_avg_neighbor_angular_sep
pure function compute_avg_neighbor_angular_sep(angles,n_particles) result (average)
  integer, intent(in) :: n_particles
  double precision, dimension(n_particles), intent(in) :: angles

  integer :: i
  double precision :: average

  average = 0
  do i = 1, n_particles
    average = average + min(minval(abs(angles(i)-angles(:i-1))),minval(abs(angles(i)-angles(i+1:))))
  end do
  average = average / dble(n_particles)

end function
!@-node:gcross.20090826091347.1422:compute_avg_neighbor_angular_sep
!@+node:gcross.20090830224709.2068:compute_average_particle_separation
pure function compute_particle_separation_average( &
  xij2, &
  n_particles &
  ) result (separation_average)
  integer, intent(in) :: n_particles
  double precision, dimension(n_particles,n_particles), intent(in) :: xij2

  integer :: particle_1, particle_2
  double precision :: separation_average

  separation_average = 0
  do particle_1 = 1, n_particles
    do particle_2 = particle_1+1, n_particles
      separation_average = separation_average + sqrt(xij2(particle_2,particle_1))
    end do
  end do
  separation_average = separation_average/(n_particles*(n_particles-1))

end function
!@-node:gcross.20090830224709.2068:compute_average_particle_separation
!@-others

end module
!@nonl
!@-node:gcross.20090818081913.1354:@thin observables.f95
!@-leo
