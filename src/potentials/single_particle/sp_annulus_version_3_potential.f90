!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1533:@thin sp_annulus_version_3_potential.f90
!@@language fortran90
module vpi_single_particle_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1534:<< Imported modules >>
  use kinds
  use constants
  use vpi_defines
  use sp_annulus_common
  !@-node:gcross.20090624144408.1534:<< Imported modules >>
  !@nl

  implicit none

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1535:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1536:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    integer :: slice, ip, nslice, np, np2, ndim
    real(kind=b8) :: Usp

    Usp = dot_product(harmonic_oscillator_coefficients,x(slice,ip,:)**2) &
          + normalized_hump_coefficient*exp(-( x(slice,ip,1)**2 + x(slice,ip,3)**2 )/hump_radius_squared_times_2)

    if ( x(slice,ip,1) .gt. extra_hump_x_threshold_hack ) then 
      Usp = Usp + normalized_hump_coefficient*exp(-( x(slice,ip,3)**2 )/hump_radius_squared_times_2)
    end if

    if ( ip .gt. np-N_PARTICLE2 ) then 
      Usp = Usp + species_2_attraction_coefficient*x(slice,ip,3)**2 
    end if

  end function Usp_func
  !@-node:gcross.20090624144408.1536:Usp
  !@+node:gcross.20090624144408.1537:gUsp
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

    do i = 1, np
      if ( x(slice,i,1) > extra_hump_x_threshold_hack ) then 
        gUsp(i,3) = gUsp(i,3) - x(slice,i,3) &
                                      * normalized_hump_coefficient/hump_radius_squared &
                                      * exp(-(x(slice,i,3)**2)/hump_radius_squared_times_2)
      end if
    end do

    gUsp(np-N_PARTICLE2+1:np,1) = gUsp(np-N_PARTICLE2+1:np,3) + &
      2.0_b8*species_2_attraction_coefficient*x(slice,np-N_PARTICLE2+1:np,3)

  end function gUsp_func
  !@-node:gcross.20090624144408.1537:gUsp
  !@-others
  !@-node:gcross.20090624144408.1535:<< Subroutines >>
  !@nl

end module vpi_single_particle_potential
!@-node:gcross.20090624144408.1533:@thin sp_annulus_version_3_potential.f90
!@-leo
