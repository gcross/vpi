!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1973:@thin sp_trial_numeric_differentiator.f90
!@@language fortran90

module sp_trial_numeric_differentiator

  use kinds
  use vpi_defines

  implicit none

contains

function numeric_grad_lap_spf( x, sl, np, ndim, nslice, grad_lntfn, lap_lntfn, spf_func ) result( y )
  real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
  real(kind=b8),intent(out) :: lap_lntfn 
  integer, intent(in) :: np, ndim, nslice, sl

  interface
    function spf_func( x, sl, nslice, np, ndim ) result( y )
      use kinds
      implicit none
      integer :: sl, nslice, np, ndim
      real(kind=b8), dimension( nslice, np , ndim ) :: x
      real(kind=b8) :: y
    end function spf_func
  end interface
  integer :: y

  real(kind=b8), dimension( nslice, np, ndim ) :: tx
  real(kind=b8), dimension( ndim ) :: fhi,flo,f
  integer :: i, k

  y = 0
  tx(sl,:,:) = x(sl,:,:)

  grad_lntfn = 0
  lap_lntfn = 0
  do i = 1, np
    do k = 1, ndim
      f(k) = spf_func( x, sl, nslice, np, ndim )
      tx(sl,i,k) = x(sl,i,k) + ntol_eps
      fhi(k) = spf_func( tx, sl, nslice, np, ndim )
      tx(sl,i,k) = x(sl,i,k) - ntol_eps
      flo(k) = spf_func( tx, sl, nslice, np, ndim )
      tx(sl,i,k) = x(sl,i,k)
    end do
!@@raw
#ifdef DEBUG
!@@end_raw
    write(1000,*) fhi
    write(1000,*) f
    write(1000,*) flo
!@@raw
#endif
!@@end_raw
    grad_lntfn(i,:) = (fhi(:)-flo(:))/(2.0_b8*ntol_eps)
    lap_lntfn = lap_lntfn + sum((fhi(:)-2.0_b8*f(:)+flo(:)))/ntol_eps**2
  end do

!  write(1000,*) grad_lntfn
!  write(1000,*) lap_lntfn

  y = 1

end function numeric_grad_lap_spf

end module sp_trial_numeric_differentiator
!@-node:gcross.20090624144408.1973:@thin sp_trial_numeric_differentiator.f90
!@-leo
