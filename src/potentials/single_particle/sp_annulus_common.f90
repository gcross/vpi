!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1500:@thin sp_annulus_common.f90
!@@language fortran90
module sp_annulus_common

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1501:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.1501:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1502:<< Variables >>
  !@+at
  ! Input quantities.
  !@-at
  !@@c
  real (kind=b8), dimension(3) :: harmonic_oscillator_coefficients
  real (kind=b8) :: hump_coefficient, hump_characteristic_radius
  real (kind=b8) :: extra_hump_x_threshold_hack
  real (kind=b8) :: species_2_attraction_coefficient

  !@+at
  ! Derived quantities.
  !@-at
  !@@c

  real (kind=b8) :: normalized_hump_coefficient
  real (kind=b8) :: hump_radius_squared, hump_radius_squared_times_2
  !@-node:gcross.20090624144408.1502:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1503:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1504:init_sp_potential
  subroutine common_init_sp_potential ()
    namelist /single_particle_potential_parameters/ &
      harmonic_oscillator_coefficients, &
      hump_coefficient, hump_characteristic_radius, &
      extra_hump_x_threshold_hack

    read(unit=10,nml=single_particle_potential_parameters)

    hump_radius_squared = hump_characteristic_radius * hump_characteristic_radius
    hump_radius_squared_times_2 = 2.0_b8 * hump_radius_squared

    normalized_hump_coefficient = hump_coefficient/(M_SQRT2PI*hump_characteristic_radius)

    write(*,*) "Using single particle annulus potential with"
    write(*,nml=single_particle_potential_parameters)

  end subroutine common_init_sp_potential
  !@-node:gcross.20090624144408.1504:init_sp_potential
  !@-others
  !@-node:gcross.20090624144408.1503:<< Subroutines >>
  !@nl

end module sp_annulus_common
!@-node:gcross.20090624144408.1500:@thin sp_annulus_common.f90
!@-leo
