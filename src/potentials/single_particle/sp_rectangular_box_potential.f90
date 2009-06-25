!@+leo-ver=4-thin
!@+node:gcross.20090624094338.1367:@thin sp_rectangular_box_potential.f90
!@@language fortran90
module sp_rectangular_box_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624094338.1373:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624094338.1373:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624094338.1372:<< Variables >>
  real (kind=b8), private :: x_length = 1_b8, x_wall_location
  real (kind=b8), private :: y_length = 1_b8, y_wall_location
  real (kind=b8), private :: z_length = 1_b8, z_wall_location
  !@-node:gcross.20090624094338.1372:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624094338.1371:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624094338.1370:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ x_length, y_length, z_length

    read(unit=10,nml=single_particle_potential_parameters)

    x_wall_location = x_length/2  
    y_wall_location = y_length/2  
    z_wall_location = z_length/2

    write(*,*) "Using single particle rectangular box potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@-node:gcross.20090624094338.1370:init_sp_potential
  !@+node:gcross.20090624094338.1368:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp,r

    Usp = 0
  !  r = sqrt(dot_product(x(slice,ip,:),x(slice,ip,:))) 
    if (  (abs(x(slice,ip,1)) > x_wall_location) &
     .or. (abs(x(slice,ip,2)) > y_wall_location) &
     .or. (abs(x(slice,ip,3)) > z_wall_location)  ) & 
    then
      Usp = realbignumber
    end if

  end function Usp_func
  !@-node:gcross.20090624094338.1368:Usp
  !@+node:gcross.20090624094338.1369:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    gUsp = 0

  end function gUsp_func

  !@-node:gcross.20090624094338.1369:gUsp
  !@-others
  !@-node:gcross.20090624094338.1371:<< Subroutines >>
  !@nl

end module sp_rectangular_box_potential
!@-node:gcross.20090624094338.1367:@thin sp_rectangular_box_potential.f90
!@-leo
