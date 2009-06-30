!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1729:@thin tb_dimer_Z_polarized_potential.f90
!@@language fortran90
module vpi_two_body_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1730:<< Imported modules >>
  use kinds
  use constants
  use vpi_defines
  !@-node:gcross.20090624144408.1730:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1731:<< Variables >>
  real (kind=b8), private :: lj_coefficient
  !@-node:gcross.20090624144408.1731:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1732:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1733:init_tb_potential
  subroutine init_tb_potential ()
    namelist /two_body_potential_parameters/ lj_coefficient

    read(unit=10,nml=two_body_potential_parameters)

    stop "Incomplete implementation."

    if(.not. use_pbc) then
      stop "Must be using periodic boundary conditions to employ the dimer Z-polarized potential."
    end if

    write(*,*) "Using two-body dimer Z-polarized potential with"
    write(*,nml=two_body_potential_parameters)
  end subroutine init_tb_potential
  !@-node:gcross.20090624144408.1733:init_tb_potential
  !@+node:gcross.20090624144408.1734:Uij
  function Uij_func( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, ip, nslice, np, ndim
    logical :: acc_flag

    real(kind=b8) :: Uij
    real(kind=b8) :: t3,tmp

    integer :: i,k

    acc_flag = .true.
    t3 = 0.0_b8
    do i=1,ip-1
      tmp =  lj_coefficient*xij2(slice,ip,i)**(-6) &
           + xij2(slice,ip,i)**(-1.5)*(1.0_b8-3.0_b8*xij2(slice,ip,i)**(-1)*(x(slice,ip,3)-x(slice,i,3))**2)
      t3 = t3 + exp(-(xij2(slice,ip,i)/(period_length*.99))**50)*tmp
    end do

    do i=ip+1,np
      tmp =  lj_coefficient*xij2(slice,ip,i)**(-6) &
           + xij2(slice,ip,i)**(-1.5)*(1.0_b8-3.0_b8*xij2(slice,ip,i)**(-1)*(x(slice,ip,3)-x(slice,i,3))**2)
      t3 = t3 + exp(-(xij2(slice,ip,i)/(period_length*.99))**50)*tmp
    end do

    Uij = t3 
  end function Uij_func
  !@-node:gcross.20090624144408.1734:Uij
  !@+node:gcross.20090624144408.1735:gUij
  function vpi_gUij_z_polarized_dimer( x, xij2, sl, nslice, np, ndim ) result ( gUij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: sl, nslice, np, ndim

    real(kind = b8), dimension ( np , ndim ) :: gUij

    real(kind=b8), dimension( nslice, np, ndim ) :: tx
    real(kind=b8), dimension( nslice, np, np ) :: txij2
    real(kind=b8), dimension( ndim ) :: fhi,flo
    integer :: i, k
    logical :: acc_flag

    gUij = 0.0
    do i = 1, np
      tx(sl,:,:) = x(sl,:,:)
      do k = 1, ndim
        txij2(sl,:,:) = xij2(sl,:,:)
        tx(sl,i,k) = x(sl,i,k) + ntol_eps
        if(use_pbc) then
          call vpi_update_xij_pbc( txij2, tx, sl, sl, i, nslice, np, ndim )
        else
          call vpi_update_xij( txij2, tx, sl, sl, i, nslice, np, ndim )
        end if
        fhi(k) = vpi_Uij_z_polarized_dimer( tx, txij2, sl, i, nslice, np, ndim, acc_flag )

        tx(sl,i,k) = x(sl,i,k) - ntol_eps
        if(use_pbc) then
          call vpi_update_xij_pbc( txij2, tx, sl, sl, i, nslice, np, ndim )
        else
          call vpi_update_xij( txij2, tx, sl, sl, i, nslice, np, ndim )
        end if
        flo(k) =  vpi_Uij_z_polarized_dimer( tx, txij2, sl, i, nslice, np, ndim, acc_flag )
        tx(sl,i,k) = x(sl,i,k)
      end do
      gUij(i,:) = fhi(:) - flo(:)
    end do
    gUij(:,:) = gUij(:,:)/(2.0_b8*ntol_eps)

  end function vpi_gUij_z_polarized_dimer
  !@-node:gcross.20090624144408.1735:gUij
  !@-others
  !@-node:gcross.20090624144408.1732:<< Subroutines >>
  !@nl

end module vpi_two_body_potential
!@-node:gcross.20090624144408.1729:@thin tb_dimer_Z_polarized_potential.f90
!@-leo
