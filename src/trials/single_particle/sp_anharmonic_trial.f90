!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1866:@thin sp_anharmonic_trial.f90
!@@language fortran90
module vpi_single_particle_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1867:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1867:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1868:<< Variables >>
  real (kind=b8), private :: coefficient_2nd_order
  real (kind=b8), private :: coefficient_4th_order
  !@-node:gcross.20090624144408.1868:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1869:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1870:init_sp_tfunc
  subroutine init_sp_tfunc ()
    namelist /single_particle_trial_function_parameters/ coefficient_2nd_order, coefficient_4th_order

    read(unit=10,nml=single_particle_trial_function_parameters)

    write(*,*) "Using anharmonic single particle trial function with"
    write(*,nml=single_particle_trial_function_parameters)
  end subroutine init_sp_tfunc
  !@-node:gcross.20090624144408.1870:init_sp_tfunc
  !@+node:gcross.20090623152316.110:tfunc
  function tfunc( x, sl, nslice, np, ndim ) result( y )
    integer :: sl, nslice, np, ndim
    real(kind=b8), dimension( nslice, np , ndim ) :: x
    real(kind=b8) :: y
    real, dimension( np ) :: psi_t

    psi_t(:)  = -coefficient_2nd_order*x(sl,:,3)**2 - coefficient_4th_order*x(sl,:,3)**4
    y = sum(psi_t)

  end function tfunc
  !@-node:gcross.20090623152316.110:tfunc
  !@+node:gcross.20090624144408.1872:grad & lapacian of tfunc
  function grad_lap_sp_tfun( x, sl, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
    real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, sl
    integer :: y

    grad_lntfn(:,1) = 0
    grad_lntfn(:,2) = 0
    grad_lntfn(:,3) = - 2*coefficient_2nd_order*x(sl,:,3) - 4*coefficient_4th_order*x(sl,:,3)**3

    lap_lntfn = - 2*coefficient_2nd_order*dble(np) - 12*coefficient_4th_order*sum(x(sl,:,3)**2)
    y = 1

  end function grad_lap_sp_tfun
  !@-node:gcross.20090624144408.1872:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1869:<< Subroutines >>
  !@nl

end module vpi_single_particle_trial
!@-node:gcross.20090624144408.1866:@thin sp_anharmonic_trial.f90
!@-leo
