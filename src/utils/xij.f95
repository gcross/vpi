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
!@+node:gcross.20090805093617.1833:compute_xij
pure subroutine compute_xij( x, n_slices, n_particles, n_dimensions, xij2 )
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension (n_dimensions,n_particles,n_slices), intent(in) :: x
  double precision, dimension (n_particles,n_particles,n_slices), intent(out) :: xij2

  integer :: slice, particle_1, particle_2
  double precision :: p1x(n_dimensions), r_sq
  do slice = 1, n_slices
    do particle_1 = 1, n_particles
      p1x = x(:,particle_1,slice)
      do particle_2 = particle_1+1, n_particles
        r_sq = sum( (p1x - x(:,particle_2,slice))**2 )
        xij2(particle_2,particle_1,slice) = r_sq
        xij2(particle_1,particle_2,slice) = r_sq
      end do
      xij2(particle_1,particle_1,slice) = 0
    end do
  end do

end subroutine compute_xij
!@-node:gcross.20090805093617.1833:compute_xij
!@+node:gcross.20090805093617.1834:compute_xij_pbc
pure subroutine compute_xij_pbc( x, period_length, n_slices, n_particles, n_dimensions, xij2 )
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension (n_dimensions,n_particles,n_slices), intent(in) :: x
  double precision, dimension (n_particles,n_particles,n_slices), intent(out) :: xij2
  double precision, intent(in) :: period_length

  integer :: slice, particle_1, particle_2
  double precision :: p1x(n_dimensions), r_sq
  do slice = 1, n_slices
    do particle_1 = 1, n_particles
      p1x = x(:,particle_1,slice)
      do particle_2 = particle_1+1, n_particles
        r_sq = sum(wrap_around(p1x - x(:,particle_2,slice),period_length)**2)
        xij2(particle_2,particle_1,slice) = r_sq
        xij2(particle_1,particle_2,slice) = r_sq
      end do
      xij2(particle_1,particle_1,slice) = 0
    end do
  end do

contains

  ! The stupid gfortran inliner will not inline this function unless it is nested...
  elemental function wrap_around(x,period_length)
    double precision, intent(in) :: x, period_length
    double precision :: wrap_around
    wrap_around = x - period_length*(floor(x/period_length - 0.5D0) + 1.0D0)
  end function

end subroutine compute_xij_pbc
!@-node:gcross.20090805093617.1834:compute_xij_pbc
!@-others

end module xij
!@-node:gcross.20090805093617.1832:@thin xij.f95
!@-leo
