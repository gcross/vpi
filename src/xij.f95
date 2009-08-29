!@+leo-ver=4-thin
!@+node:gcross.20090805093617.1832:@thin xij.f95
!@@language fortran90
module xij

implicit none

contains 

!@+others
!@+node:gcross.20090805153643.1842:wrap_around
elemental function wrap_around(x,period_length)
  double precision, intent(in) :: x, period_length
  double precision :: wrap_around
  wrap_around = x - period_length*(floor(x/period_length - 0.5D0) + 1.0D0)
end function
!@-node:gcross.20090805153643.1842:wrap_around
!@+node:gcross.20090805093617.1833:update_xij
pure subroutine update_xij( xij2, x, N_SLICES, N_PARTICLES, N_DIMENSIONS )
  integer, intent(in) :: N_SLICES, N_PARTICLES, N_DIMENSIONS
  double precision, dimension ( N_SLICES, N_PARTICLES , N_PARTICLES ), intent(inout) :: xij2
  double precision, dimension ( N_SLICES, N_PARTICLES , N_DIMENSIONS ), intent(in) :: x

  integer :: s, i, j

  forall (s=1:n_slices, i=1:N_PARTICLES, j=1:N_PARTICLES) &
    xij2(s,i,j) = sum( (x(s,i,:) - x(s,j,:))**2 )

end subroutine update_xij
!@-node:gcross.20090805093617.1833:update_xij
!@+node:gcross.20090805093617.1834:update_xij_pbc
pure subroutine update_xij_pbc( xij2, x, period_length, N_SLICES, N_PARTICLES, N_DIMENSIONS )
  integer, intent(in) :: N_SLICES, N_PARTICLES, N_DIMENSIONS
  double precision, dimension ( N_SLICES, N_PARTICLES , N_PARTICLES ), intent(inout) :: xij2
  double precision, dimension ( N_SLICES, N_PARTICLES , N_DIMENSIONS ), intent(in) :: x
  double precision, intent(in) :: period_length

  integer :: s, i, j

  forall (s=1:n_slices, i=1:N_PARTICLES, j=1:N_PARTICLES) &
    xij2(s,i,j) = sum( wrap_around( (x(s,i,:) - x(s,j,:)), period_length )**2 )

end subroutine update_xij_pbc
!@-node:gcross.20090805093617.1834:update_xij_pbc
!@+node:gcross.20090828201103.2127:hard_sphere_violation
pure function hard_sphere_violation(xij2,hard_sphere_radius_squared,n_slices,n_particles) result (is_violated)
  integer, intent(in) :: n_slices, n_particles
  double precision, intent(in) :: hard_sphere_radius_squared
  double precision, dimension(n_slices,n_particles,n_particles), intent(in) :: xij2
  logical :: is_violated

  integer :: i, j

  is_violated = .false.
  do i = 1, n_particles
    do j = i+1, n_particles
      if (any(xij2(:,i,j) <= hard_sphere_radius_squared)) then
        is_violated = .true.
        return
      end if
    end do
  end do
end function
!@-node:gcross.20090828201103.2127:hard_sphere_violation
!@-others

end module xij
!@-node:gcross.20090805093617.1832:@thin xij.f95
!@-leo
