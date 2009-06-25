!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1573:@thin sp_atomic_overlapping_potential.f90
!@@language fortran90
module sp_atomic_overlapping_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1574:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1574:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1575:<< Variables >>
  real (kind=b8), private :: coefficient_on_overlap = 1.0_b8
  real (kind=b8), private :: coefficient_otherwise = 1.0_b8

  !@-node:gcross.20090624144408.1575:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1576:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1577:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ coefficient_on_overlap, coefficient_otherwise

    read(unit=10,nml=single_particle_potential_parameters)

    write(*,*) "Using single particle atomic (1/r) potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624144408.1577:init_sp_potential
  !@+node:gcross.20090624144408.1578:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8) :: Usp

    real(kind=b8) :: ri

    ri = sqrt( dot_product(x(slice,ip,:),x(slice,ip,:)) )

    if(slice .gt. CSLICE) then
      Usp = coefficient_on_overlap * ri
    else
      Usp = coefficient_otherwise * ri
    end if

  end function Usp_func
  !@-node:gcross.20090624144408.1578:Usp
  !@+node:gcross.20090624144408.1579:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    real(kind=b8) :: r3i
    integer :: i

    if(slice .gt. CSLICE) then
      do i = 1, np
        r3i = dot_product(x(slice,i,:),x(slice,i,:))**(-3.0/2)
        gUsp(i,1) = -coefficient_on_overlap*x(slice,i,1)*r3i
        gUsp(i,2) = -coefficient_on_overlap*x(slice,i,2)*r3i
        gUsp(i,3) = -coefficient_on_overlap*x(slice,i,3)*r3i
      end do
    else
      do i = 1, np
        r3i = dot_product(x(slice,i,:),x(slice,i,:))**(-3.0/2)
        gUsp(i,1) = -coefficient_otherwise*x(slice,i,1)*r3i
        gUsp(i,2) = -coefficient_otherwise*x(slice,i,2)*r3i
        gUsp(i,3) = -coefficient_otherwise*x(slice,i,3)*r3i
      end do
    end if

  end function gUsp_func
  !@-node:gcross.20090624144408.1579:gUsp
  !@-others
  !@-node:gcross.20090624144408.1576:<< Subroutines >>
  !@nl

end module sp_atomic_overlapping_potential
!@-node:gcross.20090624144408.1573:@thin sp_atomic_overlapping_potential.f90
!@-leo
