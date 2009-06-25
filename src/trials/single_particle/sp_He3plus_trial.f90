!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1809:@thin sp_He3plus_trial.f90
!@@language fortran90
module sp_He3plus_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1810:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1810:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1811:<< Variables >>
  !@-node:gcross.20090624144408.1811:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1812:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1813:init_sp_tfunc
  subroutine init_sp_tfunc ()
    write(*,*) "Using He3+ single particle trial function."
  end subroutine init_sp_tfunc
  !@-node:gcross.20090624144408.1813:init_sp_tfunc
  !@+node:gcross.20090624144408.1814:tfunc
  function He3_fn_plus_tfun( x, islice, nslice, np, ndim ) result( y )
    integer :: islice, nslice, np, ndim
    real(kind=b8), dimension( : , : , : ) :: x

    real(kind=b8) :: y
    real(kind=b8)  :: r1,r1_sq,r2,r2_sq
    real(kind=b8)  :: rc
    integer :: i, j

    r1_sq =  dot_product(x(islice,1,:),x(islice,1,:))
    r1 = sqrt(r1_sq)
    r2_sq =  dot_product(x(islice,2,:),x(islice,2,:))
    r2 = sqrt(r2_sq)

    rc = r1-r2

    if( rc .gt. 0 ) then
      y = log(r1-r2)
    else 
      y = -realbignumber
    end if

  end function He3_fn_plus_tfun
  !@-node:gcross.20090624144408.1814:tfunc
  !@+node:gcross.20090624144408.1815:grad & lapacian of tfunc
  function grad_lap_sp_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
    real(kind=b8), dimension( : , : , : ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, slice
    integer :: y

    integer :: i, j

    grad_lntfn(:,1) = -2.0_b8*p_hox*x(slice,:,1)
    grad_lntfn(:,2) = -2.0_b8*p_hoy*x(slice,:,2)
    grad_lntfn(:,3) = -2.0_b8*p_hoz*x(slice,:,3)

    lap_lntfn = -2.0_b8*(p_hox+p_hoy+p_hoz)*np 
    y = 1

  end function grad_lap_sp_tfun


  !@-node:gcross.20090624144408.1815:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1812:<< Subroutines >>
  !@nl

end module sp_He3plus_trial
!@-node:gcross.20090624144408.1809:@thin sp_He3plus_trial.f90
!@-leo
