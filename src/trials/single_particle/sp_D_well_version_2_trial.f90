!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1953:@thin sp_D_well_version_2_trial.f90
!@@language fortran90
module sp_D_well_version_2_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1954:<< Imported modules >>
  use kinds
  use sp_D_well_trial_common
  !@-node:gcross.20090624144408.1954:<< Imported modules >>
  !@nl

  implicit none

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1955:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1956:tfunc
  function tfunc( x, sl, param, nslice, np, ndim ) result( y )
    use kinds
    implicit none
    integer :: sl, nslice, np, ndim
    real(kind=b8), dimension( nslice, np , ndim ) :: x
    real(kind=b8), dimension(:) :: param
    real(kind=b8) :: y

    y =  sum(-(x_harmonic_coefficient*x(sl,:,1)**2 + y_harmonic_coefficient*x(sl,:,2)**2)/2.0_b8 &
             - p_dw_f0*((x(sl,:,3)/p_dw_f1)**2-1.0_b8)**2)
  end function tfunc
  !@-node:gcross.20090624144408.1956:tfunc
  !@+node:gcross.20090624144408.1958:grad & lapacian of tfunc
  function grad_lap_sp_tfun( x, sl, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
    real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, sl
    integer :: y
    y = 1
    grad_lntfn(:,1) = -x_harmonic_coefficient*x(sl,:,1)
    grad_lntfn(:,2) = -y_harmonic_coefficient*x(sl,:,2)
    grad_lntfn(:,3) = -4.0*p_dw_f0*x(sl,:,3)*((x(sl,:,3)/p_dw_f1)**2-1)/p_dw_f1**2
    lap_lntfn = sum(-8.0*p_dw_f0*x(sl,:,3)**2/p_dw_f1**4-4.0*p_dw_f0*((x(sl,:,3)/p_dw_f1)**2-1)/p_dw_f1**2)
    lap_lntfn = lap_lntfn -np*(x_harmonic_coefficient+y_harmonic_coefficient)

  end function grad_lap_sp_tfun
  !@-node:gcross.20090624144408.1958:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1955:<< Subroutines >>
  !@nl

end module sp_D_well_version_2_trial
!@-node:gcross.20090624144408.1953:@thin sp_D_well_version_2_trial.f90
!@-leo
