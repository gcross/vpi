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
pure subroutine update_xij( xij2, x, n_slices, n_particles, n_dimensions )
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension ( n_slices, n_particles , n_particles ), intent(inout) :: xij2
  double precision, dimension ( n_slices, n_particles , n_dimensions ), intent(in) :: x

  integer :: s, i, j

  forall (s=1:n_slices, i=1:n_particles, j=1:n_particles) &
    xij2(s,i,j) = sum( (x(s,i,:) - x(s,j,:))**2 )

end subroutine update_xij
!@nonl
!@-node:gcross.20090805093617.1833:update_xij
!@+node:gcross.20090805093617.1834:update_xij_pbc
pure subroutine update_xij_pbc( xij2, x, period_length, n_slices, n_particles, n_dimensions )
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension ( n_slices, n_particles , n_particles ), intent(inout) :: xij2
  double precision, dimension ( n_slices, n_particles , n_dimensions ), intent(in) :: x
  double precision, intent(in) :: period_length

  integer :: s, i, j

  forall (s=1:n_slices, i=1:n_particles, j=1:n_particles) &
    xij2(s,i,j) = sum( wrap_around( (x(s,i,:) - x(s,j,:)), period_length )**2 )

end subroutine update_xij_pbc
!@nonl
!@-node:gcross.20090805093617.1834:update_xij_pbc
!@-others

end module xij
!@-node:gcross.20090805093617.1832:@thin xij.f95
!@-leo
