!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1882:@thin sp_morse_version_1_trial.f90
!@@language fortran90
module sp_morse_version_1_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1883:<< Imported modules >>
  use kinds
  use sp_morse_trial_common
  !@-node:gcross.20090624144408.1883:<< Imported modules >>
  !@nl

  implicit none

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1885:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1887:tfunc
  function tfunc( x, sl, nslice, np, ndim ) result( y )
    integer :: np, ndim, nslice,sl
    real(kind=b8), dimension( nslice, np , ndim ) :: x
    real(kind=b8) :: y
    real, dimension( np ) :: psi_t

    real(kind=b8), dimension( np ) :: r

    integer :: i, j

  !  print *,"IN morse_tfun"

    r(:) = sqrt(x(sl,:,1)**2 + x(sl,:,2)**2 + x(sl,:,3)**2)
  !  print *,r(:)
    psi_t(:)  = -( p_MO_vpa*r(:) +  p_MO_vpb/r(:)**3 )
    y = sum(psi_t)
  !  print *,y

  end function tfunc
  !@nonl
  !@-node:gcross.20090624144408.1887:tfunc
  !@+node:gcross.20090624144408.1888:grad & lapacian of tfunc
  !> q;
  !                                           b
  !                                  -a rr - ---
  !                                            3
  !                                          rr
  !
  !> simplify((grad(q,[rr,theta,phi],coords=spherical)[1]));
  !                                       4
  !                                   a rr  - 3 b
  !                                 - -----------
  !                                         4
  !                                       rr
  !
  !> simplify((laplacian(q,[rr,theta,phi],coords=spherical)));
  !                                       4
  !                                   a rr  + 3 b
  !                                -2 -----------
  !                                         5
  !                                       rr
  function grad_lap_sp_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
    real(kind=b8), dimension( : , : , : ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, slice
    integer :: y 

    real, dimension( np ) :: t_grad

    real(kind=b8), dimension( np ) :: r
    real(kind=b8), dimension( np ) :: r3
    real(kind=b8), dimension( np ) :: r5

    integer :: i, j

    r(:) = sqrt(x(slice,:,1)**2 + x(slice,:,2)**2 + x(slice,:,3)**2)
    r3(:) = r(:)**3
    r5(:) = r3(:)*r(:)**2
    t_grad(:) = -p_MO_vpa/r(:) + 3*p_MO_vpb/r5(:)
    grad_lntfn(:,1) = x(slice,:,1)*t_grad(:)
    grad_lntfn(:,2) = x(slice,:,2)*t_grad(:)
    grad_lntfn(:,3) = x(slice,:,3)*t_grad(:)

    lap_lntfn = -2.0_b8*sum( p_MO_vpa/r(:) + 3.0_b8*p_MO_vpb/r5(:) )

    y = 1

  end function grad_lap_sp_tfun

  !@-node:gcross.20090624144408.1888:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1885:<< Subroutines >>
  !@nl

end module sp_morse_version_1_trial
!@-node:gcross.20090624144408.1882:@thin sp_morse_version_1_trial.f90
!@-leo
