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
subroutine update_xij( xij2, x, sl_start, sl_end, N_SLICES, N_PARTICLES, N_DIMENSIONS )
  integer, intent(in) :: N_SLICES, N_PARTICLES, N_DIMENSIONS
  integer, intent(in) :: sl_start, sl_end
  double precision, dimension ( N_SLICES, N_PARTICLES , N_PARTICLES ), intent(inout) :: xij2
  double precision, dimension ( N_SLICES, N_PARTICLES , N_DIMENSIONS ), intent(in) :: x

  integer :: s, i, j

  forall (s=sl_start:sl_end, i=1:N_PARTICLES, j=1:N_PARTICLES) &
    xij2(s,i,j) = sum( (x(s,i,:) - x(s,j,:))**2 )

end subroutine update_xij
!@-node:gcross.20090805093617.1833:update_xij
!@+node:gcross.20090805093617.1834:update_xij_pbc
subroutine update_xij_pbc( xij2, x, period_length, sl_start, sl_end, N_SLICES, N_PARTICLES, N_DIMENSIONS )
  integer, intent(in) :: N_SLICES, N_PARTICLES, N_DIMENSIONS
  integer, intent(in) :: sl_start, sl_end
  double precision, dimension ( N_SLICES, N_PARTICLES , N_PARTICLES ), intent(inout) :: xij2
  double precision, dimension ( N_SLICES, N_PARTICLES , N_DIMENSIONS ), intent(in) :: x
  double precision, intent(in) :: period_length

  integer :: s, i, j

  forall (s=sl_start:sl_end, i=1:N_PARTICLES, j=1:N_PARTICLES) &
    xij2(s,i,j) = sum( wrap_around( (x(s,i,:) - x(s,j,:)), period_length )**2 )

end subroutine update_xij_pbc
!@-node:gcross.20090805093617.1834:update_xij_pbc
!@-others

end module xij
!@-node:gcross.20090805093617.1832:@thin xij.f95
!@-leo
