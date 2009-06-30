!@+leo-ver=4-thin
!@+node:gcross.20090629153134.1747:@thin vpi_jastrow_trial.f90
!@@language fortran90

module vpi_jastrow_trial

use kinds

implicit none

contains

subroutine init_jas_tfunc ()
end subroutine init_jas_tfunc

function jas_tfun( x, xij2, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( nslice , np, ndim ), intent(in) :: x
  real(kind=b8), dimension( nslice , np, np ), intent(in) :: xij2
  real(kind=b8)  :: y
end function jas_tfun


function grad_lap_jas_tfun( x, xij2, sl, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
  integer :: sl, np, ndim, nslice
  real(kind=b8), dimension( nslice, np, ndim ) :: x
  real(kind=b8), dimension( nslice, np, np ) :: xij2
  real(kind=b8), dimension( np, ndim ), intent(out) :: grad_lntfn 
  real(kind=b8), intent(out) :: lap_lntfn 
  integer :: y
end function grad_lap_jas_tfun

end module vpi_jastrow_trial
!@-node:gcross.20090629153134.1747:@thin vpi_jastrow_trial.f90
!@-leo
