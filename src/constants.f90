!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1462:@thin constants.f90
!@@language fortran90

module constants
  use kinds

  real(kind=b8), parameter :: M_PI = 3.141592653589793239_b8
  real(kind=b8), parameter :: M_4PI = 4.0_b8*M_PI
  real(kind=b8), parameter :: M_PI_2 = M_PI/2.0_b8
  real(kind=b8), parameter :: M_3PI_2 = 3.0_b8*M_PI/2.0_b8
  real(kind=b8), parameter :: M_2PI = 2.0_b8*M_PI
  real(kind=b8), parameter :: M_SQRT2PI = 2.506628274631_b8

  real(kind=b8), parameter :: realbignumber = 1e30

end module constants
!@-node:gcross.20090624144408.1462:@thin constants.f90
!@-leo
