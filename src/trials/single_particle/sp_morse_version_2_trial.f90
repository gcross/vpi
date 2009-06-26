!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1917:@thin sp_morse_version_2_trial.f90
!@@language fortran90
module sp_morse_version_2_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1918:<< Imported modules >>
  use kinds
  use sp_morse_trial_common
  !@-node:gcross.20090624144408.1918:<< Imported modules >>
  !@nl

  implicit none

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1919:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1921:tfunc
  function tfunc( x, np, ndim ) result( y )
    integer :: np, ndim
    real(kind=b8), dimension( np , ndim ) :: x
    real(kind=b8) :: y
    real, dimension( np ) :: psi_t

    real(kind=b8), dimension( np ) :: r

    integer :: i, j


    r(:) = sqrt(x(:,1)**2 + x(:,2)**2 + x(:,3)**2)
    psi_t(:)  = -p_MO_vpa*(r-p_MO_vpb)**2
    y = sum(psi_t)
  end function tfunc
  !@-node:gcross.20090624144408.1921:tfunc
  !@+node:gcross.20090624144408.1922:grad & lapacian of tfunc
  !
  !> q := -a*(x-x0)**2;                                        
  !                                                                      2
  !                                                      q := -a (r - x0)
  !
  !>  simplify((grad(q,[r,theta,phi],coords=spherical)[1]));   
  !                                                        -2 a (r - x0)
  !
  !>  simplify((laplacian(q,[r,theta,phi],coords=spherical))); 
  !                                                         a (3 r - 2 x0)
  !                                                      -2 --------------
  !                                                               r
  !
  function grad_lap_sp_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
    real(kind=b8), dimension( : , : , : ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, slice
    integer :: y 

    real, dimension( np ) :: t_grad
    real(kind=b8), dimension( np ) :: r

    r(:) = sqrt(x(slice,:,1)**2 + x(slice,:,2)**2 + x(slice,:,3)**2)
    t_grad(:) = -2.0_b8*p_MO_vpa*(1 - p_MO_vpb/r(:))
    grad_lntfn(:,1) = x(slice,:,1)*t_grad(:)
    grad_lntfn(:,2) = x(slice,:,2)*t_grad(:)
    grad_lntfn(:,3) = x(slice,:,3)*t_grad(:)

    lap_lntfn = -2.0_b8*p_MO_vpa*sum(3.0_b8 - 2.0_b8*p_MO_vpb/r(:))

    y = 1

  end function grad_lap_sp_tfun
  !@-node:gcross.20090624144408.1922:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1919:<< Subroutines >>
  !@nl

end module sp_morse_version_2_trial
!@-node:gcross.20090624144408.1917:@thin sp_morse_version_2_trial.f90
!@-leo
