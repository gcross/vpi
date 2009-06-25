!@+leo-ver=4-thin
!@+node:gcross.20090623152316.34:@thin vpi_potential.f90
!@@language fortran90
module vpi_potential
  !@  << Imported modules >>
  !@+node:gcross.20090623152316.35:<< Imported modules >>
  use vpi_defines
  use vpi_xij
  use vpi_aziz
  !@-node:gcross.20090623152316.35:<< Imported modules >>
  !@nl
  implicit none

contains 

!@+others
!@+node:gcross.20090623184136.3812:Rotation potentials
!@+node:gcross.20090623152316.90:undefined
function vpi_RUij_undefined( x, x_rot, xij2, slice, ip, nslice, np, ndim ) result ( Uij )
  real(kind=b8), dimension ( nslice, np , ndim ) :: x
  real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  integer :: slice, ip, nslice, np, ndim
  logical :: acc_flag

  real(kind=b8) :: Uij

  acc_flag = .true.

  STOP "Code has attempted to evaluate an undefined potential!"

end function vpi_RUij_undefined
!@-node:gcross.20090623152316.90:undefined
!@-node:gcross.20090623184136.3812:Rotation potentials
!@-others

end module vpi_potential
!@-node:gcross.20090623152316.34:@thin vpi_potential.f90
!@-leo
