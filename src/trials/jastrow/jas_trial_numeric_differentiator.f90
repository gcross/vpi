!@+leo-ver=4-thin
!@+node:gcross.20090623152316.102:@thin jas_trial_numeric_differentiator.f90
!@@language fortran90
module jas_trial_numeric_differentiator

use kinds
use vpi_defines
use vpi_xij

contains

function numeric_grad_lap_jas( x, xij2, sl, np, ndim, nslice, grad_lntfn, lap_lntfn, jfunc ) result (y)
  integer, intent(in) :: sl, np, ndim, nslice
  real(kind=b8), dimension( nslice, np, ndim ) :: x
  real(kind=b8), dimension( nslice, np, np ) :: xij2
  real(kind=b8), dimension( np, ndim ), intent(out) :: grad_lntfn
  real(kind=b8), intent(out) :: lap_lntfn

  interface
    function jfunc( x, xij2, sl, nslice, np, ndim) result(y)
      use kinds
      implicit none
      integer :: sl, nslice, np, ndim
      real(kind=b8), dimension( nslice, np, ndim ), intent(in) :: x
      real(kind=b8), dimension( nslice, np, np ), intent(in) :: xij2
      real(kind=b8)  :: y
    end function jfunc
  end interface

  integer :: y

  real(kind=b8), dimension( nslice, np, ndim ) :: tx
  real(kind=b8), dimension( nslice, np, np ) :: txij2
  real(kind=b8), dimension( ndim ) :: fhi,flo
  real(kind=b8) :: f0
  integer :: i, k

  real(kind=b8), dimension( np, ndim ) :: tgrad_lntfn
  real(kind=b8) :: tlap_lntfn, tmp
  integer :: pass

  y = 1

  grad_lntfn = 0.0_b8
  lap_lntfn = 0.0_b8
  f0 = jfunc( x, xij2, sl, nslice, np, ndim )
  do i = 1, np
    do k = 1, ndim
      tx(sl,:,:) = x(sl,:,:)
      tx(sl,i,k) = x(sl,i,k) + ntol_eps
      txij2(sl,:,:) = xij2(sl,:,:)
      if(use_pbc) then
        call vpi_update_xij_pbc( txij2, tx, sl, sl, i, nslice, np, ndim )
      else
        call vpi_update_xij( txij2, tx, sl, sl, i, nslice, np, ndim )
      end if
      fhi(k) = jfunc( tx, txij2, sl, nslice, np, ndim )

      txij2(sl,:,:) = xij2(sl,:,:)
      tx(sl,i,k) = x(sl,i,k) - ntol_eps
      if(use_pbc) then
        call vpi_update_xij_pbc( txij2, tx, sl, sl, i, nslice, np, ndim )
      else
        call vpi_update_xij( txij2, tx, sl, sl, i, nslice, np, ndim )
      end if
      flo(k) = jfunc( tx, txij2, sl, nslice, np, ndim )
    end do
    grad_lntfn(i,:) = fhi(:) - flo(:)
    tmp = sum( fhi(:) + flo(:) ) -2.0*ndim*f0
    lap_lntfn = lap_lntfn + tmp
!    write(1001,"(i10, 12g20.12)") i, flo,f0,fhi
!    write(1002,"(i10, 12g20.12)") i, grad_lntfn(i,:)
!    write(1003,"(i10, 12g20.12)") i, tmp
  end do
  grad_lntfn(:,:) = grad_lntfn(:,:)/(2.0_b8*ntol_eps)
  lap_lntfn = lap_lntfn / ntol_eps**2
!  write(1004,"(12g20.12)") grad_lntfn
!  write(1005,"(12g20.12)") lap_lntfn

!  pass = grad_lap_lj_tfun( x, xij2, sl, N_PARTICLE, N_DIM, N_SLICE, tgrad_lntfn, tlap_lntfn, jfunc_params, jfunc )

!  write(1006,"(12g20.12)") grad_lntfn(:,:) - tgrad_lntfn(:,:)
!  write(1007,"(12g20.12)") lap_lntfn - tlap_lntfn


end function numeric_grad_lap_jas

end module jas_trial_numeric_differentiator
!@-node:gcross.20090623152316.102:@thin jas_trial_numeric_differentiator.f90
!@-leo
