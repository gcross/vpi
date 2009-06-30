!@+leo-ver=4-thin
!@+node:gcross.20090629153134.1745:@thin vpi_single_particle_trial.f90
!@@language fortran90

module vpi_single_particle_trial

use kinds

implicit none

contains

subroutine init_sp_tfunc ()
end subroutine init_sp_tfunc

function tfunc( x, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( nslice, np , ndim ) :: x
  real(kind=b8) :: y
end function tfunc

function grad_lap_sp_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
  real(kind=b8), dimension( : , : , : ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn
  real(kind=b8),intent(out) :: lap_lntfn
  integer, intent(in) :: np, ndim, nslice, slice
  integer :: y
end function grad_lap_sp_tfun

end module vpi_single_particle_trial
!@-node:gcross.20090629153134.1745:@thin vpi_single_particle_trial.f90
!@-leo
