!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1512:@thin sp_annulus_version_1_potential.f90
!@@language fortran90
module vpi_single_particle_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1513:<< Imported modules >>
  use kinds
  use constants
  use sp_annulus_common
  !@-node:gcross.20090624144408.1513:<< Imported modules >>
  !@nl

  implicit none

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1515:<< Subroutines >>
  !@+others
  !@+node:gcross.20090707121444.1782:init_sp_potential
  subroutine init_sp_potential ()
    call common_init_sp_potential
  end subroutine init_sp_potential
  !@-node:gcross.20090707121444.1782:init_sp_potential
  !@+node:gcross.20090624144408.1516:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    integer :: slice, ip, nslice, np, ndim
    real(kind=b8) :: Usp

    Usp = dot_product(harmonic_oscillator_coefficients,x(slice,ip,:)**2) &
          + normalized_hump_coefficient*exp(-( x(slice,ip,1)**2 + x(slice,ip,3)**2 )/hump_radius_squared_times_2)

  !@+at
  !   if ( x(slice,ip,1) .gt. extra_hump_x_threshold_hack ) then
  !     Usp = Usp + normalized_hump_coefficient*exp(-( x(slice,ip,3)**2 
  ! )/hump_radius_squared_times_2)
  !   end if
  !@-at
  !@@c

  end function Usp_func
  !@-node:gcross.20090624144408.1516:Usp
  !@+node:gcross.20090624144408.1517:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    integer :: i

    gUsp(:,1) = 2.0_b8*harmonic_oscillator_coefficients(1)*x(slice,:,1) &
                 - normalized_hump_coefficient/hump_radius_squared*x(slice,:,1) &
                      * exp(-(x(slice,:,1)**2+x(slice,:,3)**2)/hump_radius_squared_times_2)
    gUsp(:,2) = 2.0_b8*harmonic_oscillator_coefficients(2)*x(slice,:,2)
    gUsp(:,3) = 2.0_b8*harmonic_oscillator_coefficients(3)*x(slice,:,3) &
                 - normalized_hump_coefficient/hump_radius_squared*x(slice,:,3) &
                      * exp(-(x(slice,:,1)**2+x(slice,:,3)**2)/hump_radius_squared_times_2)

  !@+at
  !   do i = 1, np
  !     if ( x(slice,i,1) > extra_hump_x_threshold_hack ) then
  !       gUsp(i,3) = gUsp(i,3) - x(slice,i,3) &
  !                                     * 
  ! normalized_hump_coefficient/hump_radius_squared &
  !                                     * 
  ! exp(-(x(slice,i,3)**2)/hump_radius_squared_times_2)
  !     end if
  !   end do
  !@-at
  !@@c
  end function gUsp_func
  !@-node:gcross.20090624144408.1517:gUsp
  !@-others
  !@-node:gcross.20090624144408.1515:<< Subroutines >>
  !@nl

end module vpi_single_particle_potential
!@-node:gcross.20090624144408.1512:@thin sp_annulus_version_1_potential.f90
!@-leo
