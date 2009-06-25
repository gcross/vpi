!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1931:@thin sp_D_well_trial_common.f90
!@@language fortran90
module sp_D_well_trial_common

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1932:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1932:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1933:<< Variables >>
  real(kind=b8) :: x_harmonic_coefficient, y_harmonic_coefficient
  real(kind=b8) :: p_dw_f0 =  0.877013259465906
  real(kind=b8) :: p_dw_f1 =  3.7487604667038
  real(kind=b8) :: p_dw_f2 =  1.4485743834981
  real(kind=b8) :: p_dw_f3 =  0.0
  real(kind=b8) :: p_dw_f4 =  0.0
  !@-node:gcross.20090624144408.1933:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1934:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1935:init_sp_tfunc
  subroutine init_sp_tfunc ()
    namelist /single_particle_trial_function_parameters/ &
      x_harmonic_coefficient, y_harmonic_coefficient, &
      p_dw_f0, p_dw_f1, p_dw_f2, p_dw_f3, p_dw_f4

    read(unit=10,nml=single_particle_trial_function_parameters)

    write(*,*) "Using D well trial function with"
    write(*,nml=single_particle_trial_function_parameters)
  end subroutine init_sp_tfunc
  !@-node:gcross.20090624144408.1935:init_sp_tfunc
  !@-others
  !@-node:gcross.20090624144408.1934:<< Subroutines >>
  !@nl

end module sp_D_well_trial_common
!@-node:gcross.20090624144408.1931:@thin sp_D_well_trial_common.f90
!@-leo
