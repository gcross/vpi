!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1714:@thin tb_hard_sphere_potential.f90
!@@language fortran90
module tb_hard_sphere_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1715:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.1715:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1716:<< Variables >>
  real (kind=b8), private :: radius
  real (kind=b8), private :: radius_squared
  !@-node:gcross.20090624144408.1716:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1717:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1718:init_tb_potential
  subroutine init_tb_potential ()
    namelist /two_body_potential_parameters/ radius

    read(unit=10,nml=two_body_potential_parameters)

    radius_squared = radius**2

    write(*,*) "Using two-body hard-sphere potential with"
    write(*,nml=two_body_potential_parameters)
  end subroutine init_tb_potential
  !@-node:gcross.20090624144408.1718:init_tb_potential
  !@+node:gcross.20090624144408.1719:Uij
  function Uij_func( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, ip, nslice, np, ndim
    logical, intent(inout):: acc_flag

    integer k
    real(kind=b8) :: Uij

    Uij = 0.0_b8

    do k=1,ip-1
      if(xij2(slice, ip, k) .le. radius_squared) then
        acc_flag = .false.
        Uij = realbignumber
        goto 111
      end if
    end do

    do k=ip+1,np
      if(xij2(slice, ip, k) .le. radius_squared) then
        acc_flag = .false.
        Uij = realbignumber
        goto 111
      end if
    end do

  111 return

  end function Uij_func

  !@-node:gcross.20090624144408.1719:Uij
  !@+node:gcross.20090624144408.1720:gUij
  function gUij_func( x, xij2, slice, nslice, np, ndim ) result ( gUij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, nslice, np, ndim

    real(kind = b8), dimension ( np , ndim ) :: gUij

    gUij = 0

  end function gUij_func
  !@-node:gcross.20090624144408.1720:gUij
  !@-others
  !@-node:gcross.20090624144408.1717:<< Subroutines >>
  !@nl

end module tb_hard_sphere_potential
!@-node:gcross.20090624144408.1714:@thin tb_hard_sphere_potential.f90
!@-leo
