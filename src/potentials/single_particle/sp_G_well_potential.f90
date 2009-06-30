!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1455:@thin sp_G_well_potential.f90
!@@language fortran90
module vpi_single_particle_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1456:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.1456:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1457:<< Variables >>
  real (kind=b8), private :: x_HO_coefficient = 1.0_b8
  real (kind=b8), private :: y_HO_coefficient = 1.0_b8
  real (kind=b8), private :: z_HO_coefficient = 1.0_b8
  real (kind=b8), private :: exponential_coefficient = 20.0_b8
  real (kind=b8), private :: exponential_characteristic_length = 0.15_b8
  real (kind=b8), private :: exponential_center = 0.0_b8

  real (kind=b8), private :: exp_coeff_with_normalization
  real (kind=b8), private :: exp_characteristic_length_squared
  real (kind=b8), private :: exp_characteristic_length_squared_times_2

  !@-node:gcross.20090624144408.1457:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1458:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1459:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ &
          x_HO_coefficient, y_HO_coefficient, z_HO_coefficient, &
          exponential_coefficient, exponential_characteristic_length, exponential_center

    read(unit=10,nml=single_particle_potential_parameters)

    exp_characteristic_length_squared = exponential_characteristic_length * exponential_characteristic_length
    exp_characteristic_length_squared_times_2 = 2.0_b8 * exponential_characteristic_length
    exp_coeff_with_normalization = exponential_coefficient / (M_SQRT2PI * exponential_characteristic_length)

    write(*,*) "Using single particle G well potential with"
    write(*,nml=single_particle_potential_parameters)

  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624144408.1459:init_sp_potential
  !@+node:gcross.20090624144408.1460:Usp
  function vpi_Usp_Gwell( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp

    Usp =   x_HO_coefficient*x(slice,ip,1)**2 &
          + y_HO_coefficient*x(slice,ip,2)**2 &
          + z_HO_coefficient*x(slice,ip,3)**2 &
          + exp_coeff_with_normalization &
              * exp(-(x(slice,ip,3)-exponential_center)**2/exp_characteristic_length_squared_times_2)

  end function vpi_Usp_Gwell
  !@-node:gcross.20090624144408.1460:Usp
  !@+node:gcross.20090624144408.1461:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    gUsp(:,1) = 2.0_b8*x_HO_coefficient*x(slice,:,1)
    gUsp(:,2) = 2.0_b8*y_HO_coefficient*x(slice,:,2)
    gUsp(:,3) = 2.0_b8*z_HO_coefficient*x(slice,:,3) - &
      (x(slice,:,3)-exponential_center) / exp_characteristic_length_squared &
        * exp_coeff_with_normalization &
           * exp(-(x(slice,:,3)-exponential_center)**2/exp_characteristic_length_squared_times_2)
  end function gUsp_func

  !@-node:gcross.20090624144408.1461:gUsp
  !@-others
  !@-node:gcross.20090624144408.1458:<< Subroutines >>
  !@nl

end module vpi_single_particle_potential
!@-node:gcross.20090624144408.1455:@thin sp_G_well_potential.f90
!@-leo
