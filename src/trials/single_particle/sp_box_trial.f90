!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1838:@thin sp_box_trial.f90
!@@language fortran90
module sp_box_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1839:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1839:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1840:<< Variables >>
  real (kind=b8), private :: box_size
  !@-node:gcross.20090624144408.1840:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1841:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1842:init_sp_tfunc
  subroutine init_sp_tfunc ()
    namelist /single_particle_trial_function_parameters/ box_size

    read(unit=10,nml=single_particle_trial_function_parameters)

    write(*,*) "Using box trial function with"
    write(*,nml=single_particle_trial_function_parameters)
  end subroutine init_sp_tfunc
  !@-node:gcross.20090624144408.1842:init_sp_tfunc
  !@+node:gcross.20090624144408.1843:tfunc
  function tfunc( x, islice, nslice, np, ndim ) result( y )
    integer :: islice, nslice, np, ndim
    real(kind=b8), dimension( nslice, np , ndim ) :: x
    real(kind=b8) :: y
    real, dimension( np ) :: psi_t

    integer :: i, j

    psi_t(:)  = sum(log(cos( M_PI*x(islice,:,:)/(2.0_b8*box_size) )))
    y = sum(psi_t)

  end function tfunc
  !@-node:gcross.20090624144408.1843:tfunc
  !@+node:gcross.20090624144408.1844:grad & lapacian of tfunc
  function grad_lap_sp_tfun( x, sl, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
    real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, sl

    integer :: i, j

    grad_lntfn(:,1) = -tan(x(sl,:,1))
    grad_lntfn(:,2) = -tan(x(sl,:,2))
    grad_lntfn(:,3) = -tan(x(sl,:,3))

    lap_lntfn = -3.0*np - sum(grad_lntfn(:,:)**2)
    y = 1

  end function grad_lap_sp_tfun
  !@-node:gcross.20090624144408.1844:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1841:<< Subroutines >>
  !@nl

end module sp_box_trial
!@-node:gcross.20090624144408.1838:@thin sp_box_trial.f90
!@-leo
