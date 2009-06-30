!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1630:@thin tb_null_potential.f90
!@@language fortran90
module vpi_two_body_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1631:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1631:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1632:<< Variables >>
  !@-node:gcross.20090624144408.1632:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1633:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1634:init_tb_potential
  subroutine init_tb_potential ()
    write(*,*) "Using null two-body potential."
  end subroutine init_tb_potential
  !@-node:gcross.20090624144408.1634:init_tb_potential
  !@+node:gcross.20090624144408.1635:Uij
  function Uij_func( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, ip, nslice, np, ndim
    logical :: acc_flag

    real(kind=b8) :: Uij

    acc_flag = .true.

    Uij = 0

  end function Uij_func
  !@-node:gcross.20090624144408.1635:Uij
  !@+node:gcross.20090624144408.1636:gUij
  function gUij_func( x, xij2, slice, nslice, np, ndim ) result ( gUij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, nslice, np, ndim

    real(kind = b8), dimension ( np , ndim ) :: gUij

    gUij(:,:) = 0

  end function gUij_func
  !@-node:gcross.20090624144408.1636:gUij
  !@-others
  !@-node:gcross.20090624144408.1633:<< Subroutines >>
  !@nl

end module vpi_two_body_potential
!@-node:gcross.20090624144408.1630:@thin tb_null_potential.f90
!@-leo
