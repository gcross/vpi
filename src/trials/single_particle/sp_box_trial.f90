!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1838:@thin sp_box_trial.f90
!@@language fortran90
!@@tabwidth -2
module vpi_single_particle_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1839:<< Imported modules >>
  use kinds
  use constants
  !@nonl
  !@-node:gcross.20090624144408.1839:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1840:<< Variables >>
  real (kind=b8), private, dimension (3) :: box_dimensions = (/1_b8,1_b8,1_b8/)
  real (kind=b8), private :: x_length, x_wall_location
  real (kind=b8), private :: y_length, y_wall_location
  real (kind=b8), private :: z_length, z_wall_location
  !@-node:gcross.20090624144408.1840:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1841:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1842:init_sp_tfunc
  subroutine init_sp_tfunc ()
    namelist /single_particle_trial_function_parameters/ box_dimensions

    read(unit=10,nml=single_particle_trial_function_parameters)

    x_length = box_dimensions(1)
    y_length = box_dimensions(2)
    z_length = box_dimensions(3)

    x_wall_location = x_length/2  
    y_wall_location = y_length/2  
    z_wall_location = z_length/2

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
    integer :: i

    y = 0

    do i = 1, np
      if( (abs(x(islice,i,1)) > x_wall_location) .or. &
          (abs(x(islice,i,2)) > y_wall_location) .or. &
          (abs(x(islice,i,3)) > z_wall_location) ) then
          y = -realbignumber
          return
      endif

      y = y + log(sum(cos( M_PI * x(islice,i,:) / box_dimensions(:) )))

    end do

  end function tfunc
  !@-node:gcross.20090624144408.1843:tfunc
  !@+node:gcross.20090624144408.1844:grad & lapacian of tfunc
  function grad_lap_sp_tfun( x, sl, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
    real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
    real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
    real(kind=b8),intent(out) :: lap_lntfn 
    integer, intent(in) :: np, ndim, nslice, sl
    integer :: y

    integer :: i, j

    grad_lntfn(:,:) = -tan(x(sl,:,:))

    lap_lntfn = -3.0*np - sum(grad_lntfn(:,:)**2)
    y = 1

  end function grad_lap_sp_tfun
  !@-node:gcross.20090624144408.1844:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1841:<< Subroutines >>
  !@nl

end module vpi_single_particle_trial
!@-node:gcross.20090624144408.1838:@thin sp_box_trial.f90
!@-leo
