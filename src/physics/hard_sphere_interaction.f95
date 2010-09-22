!@+leo-ver=4-thin
!@+node:gcross.20090901084550.2616:@thin hard_sphere_interaction.f95
!@@language fortran90

module hard_sphere_interaction

  use gfn

  implicit none

contains

!@+others
!@+node:gcross.20090901084550.2617:compute_trial_weight
subroutine compute_trial_weight( &
    xij2, &
    hard_sphere_radius, &
    n_particles, &
    weight, &
    reject_flag &
  )
  integer, intent(in) ::  n_particles
  double precision, dimension(n_particles,n_particles), intent(in) :: xij2
  double precision, intent(in) :: hard_sphere_radius
  double precision, intent(out) :: weight
  logical, intent(out) :: reject_flag
  double precision  :: r2, hard_sphere_radius_squared

  integer :: particle_1, particle_2

  hard_sphere_radius_squared = hard_sphere_radius ** 2

  weight = 0d0
  reject_flag = .false.
  do particle_1 = 1, n_particles
    do particle_2 = particle_1 + 1, n_particles
      r2 = xij2(particle_2,particle_1)
      if ( r2 > hard_sphere_radius_squared ) then
        weight = weight + log(1.0d0 - hard_sphere_radius/sqrt(r2))
      else
        reject_flag = .true.
        return
      end if
    end do
  end do

end subroutine
!@-node:gcross.20090901084550.2617:compute_trial_weight
!@+node:gcross.20090901084550.2622:accumulate_trial_derivatives
pure subroutine accumulate_trial_derivatives( &
    x, xij2, &
    hard_sphere_radius, &
    n_particles, n_dimensions, &
    grad_lntfn, lap_lntfn &
  )
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_dimensions,n_particles), intent(in) :: x
  double precision, dimension(n_particles,n_particles), intent(in) :: xij2
  double precision, intent(in) :: hard_sphere_radius
  double precision, dimension(n_dimensions,n_particles), intent(inout) :: grad_lntfn 
  double precision, intent(inout) :: lap_lntfn 

  double precision :: fi,fi2,ri,ri2,ri3,ri4
  double precision, dimension(n_dimensions) :: gtmp
  double precision :: hard_sphere_radius_squared
  integer :: particle_1, particle_2

  hard_sphere_radius_squared = hard_sphere_radius ** 2

  do particle_1 = 1, n_particles
    do particle_2 = particle_1 + 1, n_particles
      ri2 = 1.0d0/xij2(particle_2,particle_1)
      ri = sqrt(ri2)
      ri3 = ri*ri2
      ri4 = ri2*ri2
      fi = 1.0d0/(1.0d0 - hard_sphere_radius*ri)
      fi2 = fi*fi
      gtmp(:) = ri3*fi*(x(:,particle_1) - x(:,particle_2))
      grad_lntfn(:,particle_1) = grad_lntfn(:,particle_1) + hard_sphere_radius*gtmp(:)
      grad_lntfn(:,particle_2) = grad_lntfn(:,particle_2) - hard_sphere_radius*gtmp(:)
      lap_lntfn = lap_lntfn - 2.0d0*hard_sphere_radius_squared*fi2*ri4
    end do
  end do
end subroutine
!@-node:gcross.20090901084550.2622:accumulate_trial_derivatives
!@+node:gcross.20090828201103.2127:has_collision
pure function has_collision(xij2,hard_sphere_radius_squared,n_slices,n_particles)
  integer, intent(in) :: n_slices, n_particles
  double precision, intent(in) :: hard_sphere_radius_squared
  double precision, dimension(n_particles,n_particles,n_slices), intent(in) :: xij2
  logical :: has_collision

  integer :: slice, particle

  has_collision = .false.
  do slice = 1, n_slices
    do particle = 1, n_particles-1
      if (any(xij2(particle+1:,particle,slice) <= hard_sphere_radius_squared)) then
        has_collision = .true.
        return
      end if
    end do
  end do
end function
!@-node:gcross.20090828201103.2127:has_collision
!@+node:gcross.20090828095451.1455:compute_greens_function
pure function compute_greens_function( &
    xij2, &
    dt, hard_sphere_radius, &
    slice_start, slice_end, &
    particle, &
    n_slices, n_particles &
  ) result ( ln_gfn )
  integer, intent(in) :: slice_start, slice_end, particle
  integer, intent(in) :: n_slices, n_particles
  double precision, intent(in) :: dt, hard_sphere_radius
  double precision, dimension (n_particles,n_particles,n_slices), intent(in) :: xij2

  double precision :: gfn, ln_gfn
  integer :: other_particle

  gfn = 1.0d0

  do other_particle = 1, particle-1
    gfn = gfn * compute_contribution_from_particle(other_particle)
  end do

  do other_particle = particle+1, n_particles
    gfn = gfn * compute_contribution_from_particle(other_particle)
  end do

  ln_gfn = log(gfn)

contains

  pure function compute_contribution_from_particle(other_particle) result (gfn)
    integer, intent(in) :: other_particle
    double precision :: gfn

    double precision, dimension(n_slices) :: distances

    distances(slice_start:slice_end) = sqrt(xij2(other_particle,particle,slice_start:slice_end)) - hard_sphere_radius

    gfn = compute_green_fn_from_distances( &
            distances, &
            dt, &
            slice_start, slice_end, &
            n_slices &
          )
  end function

end function
!@nonl
!@-node:gcross.20090828095451.1455:compute_greens_function
!@-others

end module
!@-node:gcross.20090901084550.2616:@thin hard_sphere_interaction.f95
!@-leo
