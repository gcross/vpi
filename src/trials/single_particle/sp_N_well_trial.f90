!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1966:@thin sp_N_well_trial.f90
!@@language fortran90
module sp_N_well_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1967:<< Imported modules >>
  use kinds
  use sp_trial_numeric_differentiator
  !@-node:gcross.20090624144408.1967:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1968:<< Variables >>
  real (kind=b8), private :: coefficient_2nd_order
  real (kind=b8), private :: coefficient_4th_order
  !@-node:gcross.20090624144408.1968:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1969:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1970:init_sp_tfunc
  subroutine init_sp_tfunc ()
    namelist /single_particle_trial_function_parameters/ coefficient_2nd_order, coefficient_4th_order

    read(unit=10,nml=single_particle_trial_function_parameters)

    write(*,*) "Using anharmonic single particle trial function with"
    write(*,nml=single_particle_trial_function_parameters)
  end subroutine init_sp_tfunc
  !@-node:gcross.20090624144408.1970:init_sp_tfunc
  !@+node:gcross.20090624144408.1971:tfunc
  function tfunc( x, np, ndim ) result( y )
    integer :: np, ndim
    real(kind=b8), dimension( np , ndim ) :: x
    real(kind=b8) :: y
    real(kind=b8), dimension( np ) :: psi_t
    real(kind=b8) :: lam_nw = 1./10

    integer :: i, j

    psi_t(:)  = -( x(:,1)**2 + x(:,2)**2 + x(:,3)**2 )/2.0
    y = sum(psi_t)

  end function tfunc
  !@-node:gcross.20090624144408.1971:tfunc
  !@+node:gcross.20090624144408.1972:grad & lapacian of tfunc
  function grad_lap_sp_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
    implicit none

    real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, slice
    integer :: y

    y = numeric_grad_lap_spf( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn, tfunc )

  end function grad_lap_sp_tfun
  !@-node:gcross.20090624144408.1972:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1969:<< Subroutines >>
  !@nl

end module sp_N_well_trial
!@-node:gcross.20090624144408.1966:@thin sp_N_well_trial.f90
!@-leo
