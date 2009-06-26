!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1533:@thin sp_annulus_version_3_potential.f90
!@@language fortran90
module sp_annulus_version_3_potential

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

    Usp = ( x(slice,ip,2)**2 + ap_lam*( x(slice,ip,1)**2  + x(slice,ip,3)**2) ) + &
          ap_e_norm*exp(-( x(slice,ip,1)**2 + x(slice,ip,3)**2 )/ap_2asq)

    if ( x(slice,ip,1) .gt. ab_x0 ) then 
      Usp = Usp + ab_e_norm*exp(-( x(slice,ip,3)**2 )/ab_2asq)
    end if

    if ( ip .gt. np-N_PARTICLE2 ) then 
      Usp = Usp + ae3*x(slice,ip,3)**2 
    end if

  end function Usp_func
  !@-node:gcross.20090624144408.1536:Usp
  !@+node:gcross.20090624144408.1537:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    integer :: i

    gUsp(:,1) = 2.0_b8*ap_lam*x(slice,:,1) - &
      x(slice,:,1)*ap_e_norm*exp(-(x(slice,:,1)**2 + &
      x(slice,:,3)**2)/ap_2asq)/ap_asq 
    gUsp(np-N_PARTICLE2+1:np,1) = gUsp(np-N_PARTICLE2+1:np,3) + &
      2.0_b8*ae3*x(slice,np-N_PARTICLE2+1:np,3)
    gUsp(:,2) = 2.0_b8*x(slice,:,2)
    gUsp(:,3) = 2.0_b8*ap_lam*x(slice,:,3) - &
      x(slice,:,3)*ap_e_norm*exp(-(x(slice,:,2)**2+x(slice,:,3)**2)/ap_2asq)/ap_asq 

    do i = 1, np
      if ( x(slice,i,1) .gt. ab_x0 ) then 
        gUsp(i,3) = gUsp(i,3) + x(slice,i,3)*ab_e_norm*exp(-(x(slice,i,3)**2)/ab_2asq)/ab_asq 
      end if
    end do

  end function gUsp_func
  !@-node:gcross.20090624144408.1537:gUsp
  !@-others
  !@-node:gcross.20090624144408.1535:<< Subroutines >>
  !@nl

end module sp_annulus_version_3_potential
!@-node:gcross.20090624144408.1533:@thin sp_annulus_version_3_potential.f90
!@-leo
