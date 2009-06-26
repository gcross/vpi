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
  real :: ab_e = 0
  real :: ab_a = 0.1
  real :: ab_asq 
  real :: ab_2asq
  real :: ab_e_norm
  real :: ab_x0 = -100

  real :: ap_e = 20
  real :: ap_a = 1.5
  real :: ap_asq 
  real :: ap_2asq
  real :: ap_e_norm
  real :: ap_lam = 1.0
  real :: ae3 = 4  ! extra harmonic confinment for second species in Annulus3
  !@-node:gcross.20090624144408.1502:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1503:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1504:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ ab_e, ab_a, ab_x0, ap_e, ap_a, ap_lam

    read(unit=10,nml=single_particle_potential_parameters)

    ab_asq = ab_a*ab_a
    ab_2asq = 2.0*ab_asq
    ab_e_norm = ab_e/(M_SQRT2PI*ab_a)

    ap_asq = ap_a*ap_a
    ap_2asq = 2.0*ap_asq
    ap_e_norm = ap_e/(M_SQRT2PI*ap_a)

    write(*,*) "Using single particle annulus potential with"
    write(*,nml=single_particle_potential_parameters)

  end subroutine init_sp_potential
  !@-node:gcross.20090624144408.1504:init_sp_potential
  !@-others
  !@-node:gcross.20090624144408.1503:<< Subroutines >>
  !@nl

end module sp_annulus_common
!@-node:gcross.20090624144408.1500:@thin sp_annulus_common.f90
!@-leo
