!@+leo-ver=4-thin
!@+node:gcross.20090805153643.1843:@thin thermalize.f90
!@@language fortran90
module thermalize

use xij
use sample

implicit none

contains 

!@+others
!@+node:gcross.20090623152316.32:accept_path
function accept_path( lngfn0, lngfn1 ) result( accept )
  double precision, intent(in) :: lngfn0, lngfn1
  logical :: accept

  real :: Pa, Ptest 

  Pa = exp(lngfn1 - lngfn0) 
  call random_number( Ptest )
  accept = ( Pa > Ptest )

end function accept_path
!@-node:gcross.20090623152316.32:accept_path
!@-others

end module thermalize
!@-node:gcross.20090805153643.1843:@thin thermalize.f90
!@-leo
