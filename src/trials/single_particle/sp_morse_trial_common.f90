!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1904:@thin sp_morse_trial_common.f90
!@@language fortran90
module sp_morse_trial_common

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1905:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1905:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1906:<< Variables >>
  ! variational parameters for Morse potential trial function
  real(kind=b8) :: p_MO_vpa = 15.29_b8
  real(kind=b8) :: p_MO_vpb = 6.82_b8
  !@-node:gcross.20090624144408.1906:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1907:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1908:init_sp_tfunc
  subroutine init_sp_tfunc ()
    namelist /single_particle_trial_function_parameters/ p_MO_vpa, p_MO_vpb

    read(unit=10,nml=single_particle_trial_function_parameters)

    write(*,*) "Using Morse trial function with"
    write(*,nml=single_particle_trial_function_parameters)
  end subroutine init_sp_tfunc
  !@-node:gcross.20090624144408.1908:init_sp_tfunc
  !@-others
  !@-node:gcross.20090624144408.1907:<< Subroutines >>
  !@nl

end module sp_morse_trial_common
!@-node:gcross.20090624144408.1904:@thin sp_morse_trial_common.f90
!@-leo
