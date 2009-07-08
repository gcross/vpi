!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1644:@thin tb_charge_potential.f90
!@@language fortran90
module vpi_two_body_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1645:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1645:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1646:<< Variables >>
  real (kind=b8), private :: coefficient
  !@-node:gcross.20090624144408.1646:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1647:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1648:init_tb_potential
  subroutine init_tb_potential ()
    namelist /two_body_potential_parameters/ coefficient

    read(unit=10,nml=two_body_potential_parameters)

    write(*,*) "Using two-body charge potential with"
    write(*,nml=two_body_potential_parameters)
  end subroutine init_tb_potential
  !@-node:gcross.20090624144408.1648:init_tb_potential
  !@+node:gcross.20090624144408.1649:Uij
  function Uij_func( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, ip, nslice, np, ndim
    logical :: acc_flag

    real(kind=b8) :: Uij

    acc_flag = .true.

    Uij = coefficient * ( sum(xij2(slice,ip,1:ip-1)**(-0.5)) + sum(xij2(slice,ip,ip+1:np)**(-0.5)) )

  end function Uij_func
  !@-node:gcross.20090624144408.1649:Uij
  !@+node:gcross.20090624144408.1650:gUij
  function gUij_func( x, xij2, slice, nslice, np, ndim ) result ( gUij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, nslice, np, ndim

    real(kind = b8), dimension ( np , ndim ) :: gUij

    real(kind = b8), dimension ( ndim ) :: tg_ee
    real(kind = b8), dimension ( np , ndim ) :: g_ee
    real(kind=b8) :: rij3i
    integer :: i,j

    g_ee = 0

    do i = 1, np
      do j = i+1, np
        rij3i =  xij2(slice,i,j)**(-3.0/2)
        tg_ee(:) = ( x(slice,i,:) - x(slice,j,:) ) * rij3i
        g_ee(i,:) = g_ee(i,:) + tg_ee(:)
        g_ee(j,:) = g_ee(j,:) - tg_ee(:)
      end do
    end do

    gUij(:,1) = -coefficient*g_ee(:,1)
    gUij(:,2) = -coefficient*g_ee(:,2)
    gUij(:,3) = -coefficient*g_ee(:,3)

  end function gUij_func
  !@-node:gcross.20090624144408.1650:gUij
  !@-others
  !@-node:gcross.20090624144408.1647:<< Subroutines >>
  !@nl

end module vpi_two_body_potential
!@-node:gcross.20090624144408.1644:@thin tb_charge_potential.f90
!@-leo
