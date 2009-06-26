!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1471:@thin sp_D_well_potential.f90
!@@language fortran90
module sp_D_well_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1472:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.1472:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1473:<< Variables >>
  real(kind=b8), private :: length_scale = 1.5_b8
  real(kind=b8), private :: length_scale_squared
  real(kind=b8), private :: ep = 40.0_b8
  !@nonl
  !@-node:gcross.20090624144408.1473:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1474:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1475:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ ep, length_scale

    read(unit=10,nml=single_particle_potential_parameters)

    length_scale_squared = length_scale * length_scale

    write(*,*) "Using single particle D well potential with"
    write(*,nml=single_particle_potential_parameters)

  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624144408.1475:init_sp_potential
  !@+node:gcross.20090624144408.1476:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp

    Usp = ( x(slice,ip,1)**2 + x(slice,ip,2)**2 )/2.0_b8 + ep*((x(slice,ip,3)/length_scale)**2 - 1.0_b8)**2

  end function Usp_func
  !@-node:gcross.20090624144408.1476:Usp
  !@+node:gcross.20090624144408.1477:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    gUsp(:,1) = x(slice,:,1)
    gUsp(:,2) = x(slice,:,2)
    gUsp(:,3) = 4.0_b8*ep*((x(slice,:,3)/length_scale)**2 - 1.0_b8)*(x(slice,:,3)/length_scale_squared) 

  end function gUsp_func
  !@-node:gcross.20090624144408.1477:gUsp
  !@-others
  !@-node:gcross.20090624144408.1474:<< Subroutines >>
  !@nl

end module sp_D_well_potential
!@-node:gcross.20090624144408.1471:@thin sp_D_well_potential.f90
!@-leo
