!@+leo-ver=4-thin
!@+node:gcross.20090624094338.1381:@thin sp_cylindrical_box_potential.f90
!@@language fortran90
module sp_cylindrical_box_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624094338.1382:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624094338.1382:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624094338.1383:<< Variables >>
  real (kind=b8), private :: cylinder_radius = 1.0_b8
  !@-node:gcross.20090624094338.1383:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624094338.1384:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624094338.1385:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ cylinder_radius

    read(unit=10,nml=single_particle_potential_parameters)

    write(*,*) "Using single particle cylindrical box potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624094338.1385:init_sp_potential
  !@+node:gcross.20090624094338.1386:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp,r

    Usp = 0
    r = sqrt( x(slice,ip,1)**2 + x(slice,ip,2)**2 ) 
    if(r > cylinder_radius) then
      Usp = realbignumber*r**2
    end if

  end function Usp_func
  !@-node:gcross.20090624094338.1386:Usp
  !@+node:gcross.20090624094338.1387:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    integer :: slice, nslice, np, ndim
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    real(kind=b8), dimension ( np, ndim ) :: gUsp
    real(kind=b8), dimension ( np ) :: r

    gUsp = 0

  end function gUsp_func
  !@-node:gcross.20090624094338.1387:gUsp
  !@-others
  !@-node:gcross.20090624094338.1384:<< Subroutines >>
  !@nl

end module sp_cylindrical_box_potential
!@-node:gcross.20090624094338.1381:@thin sp_cylindrical_box_potential.f90
!@-leo
