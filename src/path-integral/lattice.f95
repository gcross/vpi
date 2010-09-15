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
  double precision, dimension( n_slices, n_particles, n_dimensions ), intent(out) :: q
  double precision, intent(in) :: size
  integer, intent(in) :: n_slices, n_particles, n_dimensions

  integer :: i,k
  double precision, dimension( n_dimensions ) :: nu

  do i = 1, n_particles
    call random_number( nu )

    q(1,i,:) = (nu(:)-0.5)*size

    do k = 1, n_dimensions
      q(:,i,k) = q(1,i,k)
    end do

  end do 

end subroutine make_lattice
!@nonl
!@-node:gcross.20090819093822.1397:make_lattice
!@-others

end module lattice
!@-node:gcross.20090819093822.1396:@thin lattice.f95
!@-leo
