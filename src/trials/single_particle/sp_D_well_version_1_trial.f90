!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1941:@thin sp_D_well_version_1_trial.f90
!@@language fortran90
module sp_D_well_version_1_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1942:<< Imported modules >>
  use kinds
  use sp_D_well_trial_common
  !@-node:gcross.20090624144408.1942:<< Imported modules >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1943:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1944:tfunc
  function tfunc( x, sl, nslice, np, ndim ) result( y )
    integer :: sl, nslice, np, ndim
    real(kind=b8), dimension( nslice, np , ndim ) :: x
    real(kind=b8) :: y
    real(kind=b8), dimension( np ) :: psi_t

    y = 0

    psi_t(:)  = log( exp( -( x_harmonic_coefficient*x(sl,:,1)**2 + y_harmonic_coefficient*x(sl,:,2)**2 )/2.0 ) * &
                   ( p_dw_f0*exp(-p_dw_f1*(x(sl,:,3)-p_dw_f2)**2) + &
                   p_dw_f0*exp(-p_dw_f1*(x(sl,:,3)+p_dw_f2)**2) - &
                   p_dw_f3*exp(-p_dw_f4*(x(sl,:,3)**2))) )

  !  psi_t(:)  = -( x(:,1)**2 + x(:,2)**2 + x(:,3)**2 )/2.0
    y = sum(psi_t)

  end function tfunc
  !@-node:gcross.20090624144408.1944:tfunc
  !@+node:gcross.20090624144408.1946:grad
  function grad_sp_tfun( x, ndim, grad_lntfn ) result( y )
    real(kind=b8), dimension( ndim ), intent(in) :: x
    real(kind=b8), dimension( ndim ), intent(out) :: grad_lntfn 
    integer, intent(in) :: ndim 
    integer :: y

    real(kind=b8), dimension( ndim ) :: tx
    real(kind=b8), dimension( ndim ) :: fhi,flo
    integer :: k

    grad_lntfn = 0
    tx(:) = x(:)
    do k = 1, ndim
      tx(k) = x(k) + ntol_eps
      fhi(k)  = log( exp( -( x_harmonic_coefficient*tx(1)**2 + y_harmonic_coefficient*tx(2)**2 )/2.0 ) &
                              * ( p_dw_f0*exp(-p_dw_f1*(tx(3)-p_dw_f2)**2) &
               + p_dw_f0*exp(-p_dw_f1*(tx(3)+p_dw_f2)**2) &
               + p_dw_f3*exp(-p_dw_f4*tx(3)**2)))

      tx(k) = x(k) - ntol_eps
      flo(k)  = log( exp( -( x_harmonic_coefficient*tx(1)**2 + y_harmonic_coefficient*tx(2)**2 )/2.0 ) &
                              * ( p_dw_f0*exp(-p_dw_f1*(tx(3)-p_dw_f2)**2) &
                + p_dw_f0*exp(-p_dw_f1*(tx(3)+p_dw_f2)**2) &
                + p_dw_f3*exp(-p_dw_f4*tx(3)**2)))

      tx(k) = x(k)
    end do

    grad_lntfn(:) = (fhi(:)-flo(:))/(2.0_b8*ntol_eps)

    y = 1

  end function grad_sp_tfun
  !@-node:gcross.20090624144408.1946:grad
  !@+node:gcross.20090624144408.1945:grad & lapacian of tfunc
  function grad_lap_sp_tfun( x, sl, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
    real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, sl
    integer :: y

    real(kind=b8), dimension( ndim ) :: tgrad_lntfn 
    real(kind=b8), dimension( ndim ) :: tx
    real(kind=b8), dimension( ndim ) :: fhi,flo
    integer :: i, j, k

    lap_lntfn = 0
    grad_lntfn = 0
    do i = 1, np
      tx(:) = x(sl,i,:)
      do k = 1, ndim
        tx(k) = x(sl,i,k) + ntol_eps
        y = grad_sp_tfun( tx, ndim, tgrad_lntfn )
        fhi(k) = tgrad_lntfn(k)
        tx(k) = x(sl,i,k) - ntol_eps
        y  = grad_sp_tfun( tx, ndim, tgrad_lntfn )
        flo(k) = tgrad_lntfn(k)
        tx(k) = x(sl,i,k)
      end do
      y = grad_sp_tfun( tx, ndim, tgrad_lntfn )
      grad_lntfn(i,:) = tgrad_lntfn(:)
      lap_lntfn = lap_lntfn + sum((fhi(:)-flo(:))/(2.0_b8*ntol_eps))
    end do

  end function grad_lap_sp_tfun
  !@-node:gcross.20090624144408.1945:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1943:<< Subroutines >>
  !@nl

end module sp_D_well_version_1_trial
!@-node:gcross.20090624144408.1941:@thin sp_D_well_version_1_trial.f90
!@-leo
