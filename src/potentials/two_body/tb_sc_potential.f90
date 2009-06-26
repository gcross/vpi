!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1700:@thin tb_sc_potential.f90
!@@language fortran90
module tb_sc_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1701:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.1701:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1702:<< Variables >>
  real (kind=b8), private :: coefficient = 2352.5_b8
  real (kind=b8), private :: length_scale = 0.05_b8
  real (kind=b8), private :: length_scaling_factor
  real (kind=b8), private :: derivative_coefficient
  !@-node:gcross.20090624144408.1702:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1703:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1704:init_tb_potential
  subroutine init_tb_potential ()
    namelist /two_body_potential_parameters/ coefficient, length_scale

    read(unit=10,nml=two_body_potential_parameters)

    length_scaling_factor = 1./(2*length_scale**2)
    derivative_coefficient = coefficient/(length_scale*M_SQRT2PI)

    write(*,*) "Using two-body S.C. potential with"
    write(*,nml=two_body_potential_parameters)
  end subroutine init_tb_potential
  !@-node:gcross.20090624144408.1704:init_tb_potential
  !@+node:gcross.20090624144408.1705:Uij
  function vpi_Uij_Sc( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, ip, nslice, np, ndim
    logical :: acc_flag

    integer k
    real(kind=b8) :: Uij

    acc_flag = .true.

    Uij = 0

    do k = 1,ip-1
      Uij = Uij + exp(-length_scaling_factor*xij2(slice,ip,k))
    end do
    do k = ip+1,np
      Uij = Uij + exp(-length_scaling_factor*xij2(slice,ip,k))
    end do

    Uij = Uij * coefficient

  end function vpi_Uij_Sc
  !@-node:gcross.20090624144408.1705:Uij
  !@+node:gcross.20090624144408.1706:gUij
  function vpi_gUij_Sc( x, xij2, slice, nslice, np, ndim ) result ( gUij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, nslice, np, ndim

    real(kind = b8), dimension ( np , ndim ) :: gUij

    real(kind = b8), dimension ( ndim ) :: t_gUij
    real(kind = b8) :: Ui
    real(kind = b8) :: gc1,gc2
    integer :: i,j,k

    gUij = 0

    do i = 1, np

      Ui = 0
      do k = 1,i-1
        Ui = Ui + exp(-length_scaling_factor*xij2(slice,i,k))
      end do
      do k = i+1,np
        Ui = Ui + exp(-length_scaling_factor*xij2(slice,i,k))
      end do

      Ui = derivative_coefficient * Ui

      do j = i+1, np
        t_gUij(:) = ( x(slice,i,:) - x(slice,j,:) ) * Ui 
        gUij(i,:) = gUij(i,:) - 2.0*length_scaling_factor*t_gUij(:)
        gUij(j,:) = gUij(j,:) + 2.0*length_scaling_factor*t_gUij(:)
      end do

    end do

  end function vpi_gUij_Sc
  !@-node:gcross.20090624144408.1706:gUij
  !@-others
  !@-node:gcross.20090624144408.1703:<< Subroutines >>
  !@nl

end module tb_sc_potential
!@-node:gcross.20090624144408.1700:@thin tb_sc_potential.f90
!@-leo
