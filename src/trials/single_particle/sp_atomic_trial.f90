!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1791:@thin sp_atomic_trial.f90
!@@language fortran90
module vpi_single_particle_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1792:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1792:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1793:<< Variables >>
  real (kind=b8), private :: coefficient_0, coefficient_1
  !@-node:gcross.20090624144408.1793:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1794:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1795:init_sp_tfunc
  subroutine init_sp_tfunc ()
    namelist /single_particle_trial_function_parameters/ coefficient_0, coefficient_1

    read(unit=10,nml=single_particle_trial_function_parameters)

    write(*,*) "Using atomic single particle trial function with"
    write(*,nml=single_particle_trial_function_parameters)
  end subroutine init_sp_tfunc
  !@-node:gcross.20090624144408.1795:init_sp_tfunc
  !@+node:gcross.20090624144408.1796:tfunc
  function tfunc( x, sl, nslice, np, ndim ) result( y )
    integer :: sl, nslice, np, ndim
    real(kind=b8), dimension( nslice, np , ndim ) :: x

    real(kind=b8) :: y

    real(kind=b8)  :: r2,r1
    integer :: i

    y = 0.0
    do i = 1, np
      r2 = dot_product(x(sl,i,:),x(sl,i,:))
      r1 = sqrt(r2)
      y = y - ( coefficient_0*r1 + coefficient_1*r2 )/( 1 + coefficient_1*r1 )
    end do

  end function tfunc
  !@-node:gcross.20090624144408.1796:tfunc
  !@+node:gcross.20090624144408.1797:grad & lapacian of tfunc
  function grad_lap_sp_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
    real(kind=b8), dimension( : , : , : ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn
    real(kind=b8),intent(out) :: lap_lntfn
    integer, intent(in) :: np, ndim, nslice, slice
    integer :: y

    real, dimension( np ) :: t_grad

    real(kind=b8), dimension( np ) :: r
    real(kind=b8), dimension( np ) :: r2
    real(kind=b8), dimension( np ) :: r3

    r2(:) = x(slice,:,1)**2 + x(slice,:,2)**2 + x(slice,:,3)**2
    r(:) = sqrt(r2(:))
    r3(:) = r(:)*r2(:)
    t_grad(:) = -(coefficient_0+coefficient_1*(2*r(:)+coefficient_1*r2(:)))/(r(:)*(1+coefficient_1*r2(:))**2)
    grad_lntfn(:,1) = x(slice,:,1)*t_grad(:)
    grad_lntfn(:,2) = x(slice,:,2)*t_grad(:)
    grad_lntfn(:,3) = x(slice,:,3)*t_grad(:)

    lap_lntfn = sum( &
                  (-2.0_b8*(coefficient_0 + coefficient_1*(3.0_b8*r(:) + 3.0_b8*coefficient_1*r2(:) + coefficient_1**2*r3(:) ))) &
                  / (r(:)*(1.0_b8+coefficient_1*r2(:))**3) &
              )

    y = 1

  end function grad_lap_sp_tfun
  !@-node:gcross.20090624144408.1797:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1794:<< Subroutines >>
  !@nl

end module vpi_single_particle_trial
!@-node:gcross.20090624144408.1791:@thin sp_atomic_trial.f90
!@-leo
