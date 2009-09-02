!@+leo-ver=4-thin
!@+node:gcross.20090901084550.2616:@thin hard_sphere_interaction.f95
!@@language fortran90

module hard_sphere_interaction
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
  double precision, dimension( n_particles, n_particles ), intent(in) :: xij2
  double precision, intent(in) :: hard_sphere_radius
  double precision, intent(out) :: weight
  logical, intent(out) :: reject_flag
  double precision  :: r2, hard_sphere_radius_squared

  integer :: i, j

  hard_sphere_radius_squared = hard_sphere_radius ** 2

  weight = 0d0
  reject_flag = .false.
  do i = 1, n_particles
    do j = i + 1, n_particles
      r2 = xij2(i,j)
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
  double precision, dimension( n_particles, n_dimensions ), intent(in) :: x
  double precision, dimension( n_particles, n_particles ), intent(in) :: xij2
  double precision, intent(in) :: hard_sphere_radius
  double precision, dimension( n_particles, n_dimensions ), intent(inout) :: grad_lntfn 
  double precision, intent(inout) :: lap_lntfn 

  double precision :: fi,fi2,ri,ri2,ri3,ri4
  double precision, dimension( n_dimensions ) :: gtmp
  double precision :: hard_sphere_radius_squared
  integer :: i, j

  hard_sphere_radius_squared = hard_sphere_radius ** 2

  grad_lntfn = 0.0d0
  lap_lntfn = 0.0d0
  do i = 1, n_particles
    do j = i + 1, n_particles
      ri2 = 1.0d0/xij2(i,j)
      ri = sqrt(ri2)
      ri3 = ri*ri2
      ri4 = ri2*ri2
      fi = 1.0d0/(1.0d0 - hard_sphere_radius*ri)
      fi2 = fi*fi
      gtmp(:) = ri3*fi*(x(i,:) - x(j,:))
      grad_lntfn(i,:) = grad_lntfn(i,:) + hard_sphere_radius*gtmp(:)
      grad_lntfn(j,:) = grad_lntfn(j,:) - hard_sphere_radius*gtmp(:)
      lap_lntfn = lap_lntfn - 2.0d0*hard_sphere_radius_squared*fi2*ri4
    end do
  end do
end subroutine
!@-node:gcross.20090901084550.2622:accumulate_trial_derivatives
!@+node:gcross.20090828201103.2127:has_collision
pure function has_collision(xij2,hard_sphere_radius_squared,n_slices,n_particles)
  integer, intent(in) :: n_slices, n_particles
  double precision, intent(in) :: hard_sphere_radius_squared
  double precision, dimension(n_slices,n_particles,n_particles), intent(in) :: xij2
  logical :: has_collision

  integer :: i, j

  has_collision = .false.
  do i = 1, n_particles
    do j = i+1, n_particles
      if (any(xij2(:,i,j) <= hard_sphere_radius_squared)) then
        has_collision = .true.
        return
      end if
    end do
  end do
end function
!@-node:gcross.20090828201103.2127:has_collision
!@-others

end module
!@-node:gcross.20090901084550.2616:@thin hard_sphere_interaction.f95
!@-leo
