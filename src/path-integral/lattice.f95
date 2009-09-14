!@+leo-ver=4-thin
!@+node:gcross.20090819093822.1396:@thin lattice.f95
!@@language fortran90

module lattice

  use rand_utils

  implicit none

contains 

!@+others
!@+node:gcross.20090819093822.1397:make_lattice
subroutine make_lattice( size, nslice, np, ndim, q )
  double precision, dimension( nslice, np, ndim ), intent(out) :: q
  double precision, intent(in) :: size
  integer, intent(in) :: nslice, np, ndim

  integer :: i,k
  double precision, dimension( ndim ) :: nu

  do i = 1, np
    call random_number( nu )

    q(1,i,:) = (nu(:)-0.5)*size

    do k = 1, ndim
      q(:,i,k) = q(1,i,k)
    end do

  end do 

end subroutine make_lattice
!@-node:gcross.20090819093822.1397:make_lattice
!@-others

end module lattice
!@-node:gcross.20090819093822.1396:@thin lattice.f95
!@-leo
