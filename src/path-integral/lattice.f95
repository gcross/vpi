!@+leo-ver=4-thin
!@+node:gcross.20090819093822.1396:@thin lattice.f95
!@@language fortran90

module lattice

  use rand_utils

  implicit none

contains 

!@+others
!@+node:gcross.20090819093822.1397:make_lattice
subroutine make_lattice( size, n_slices, n_particles, n_dimensions, q )
  double precision, dimension(n_dimensions,n_particles,n_slices), intent(out) :: q
  double precision, intent(in) :: size
  integer, intent(in) :: n_slices, n_particles, n_dimensions

  integer :: particle, slice
  double precision, dimension( n_dimensions ) :: nu

  do particle = 1, n_particles
    call random_number( nu )
    q(:,particle,1) = (nu(:)-0.5)*size
  end do

  do slice = 2, n_slices
    q(:,:,slice) = q(:,:,1)
  end do

end subroutine make_lattice
!@-node:gcross.20090819093822.1397:make_lattice
!@-others

end module lattice
!@-node:gcross.20090819093822.1396:@thin lattice.f95
!@-leo
