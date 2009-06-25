!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1559:@thin sp_atomic_potential.f90
!@@language fortran90
module sp_atomic_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1560:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1560:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1561:<< Variables >>
  real (kind=b8), private :: coefficient = 1.0_b8
  !@-node:gcross.20090624144408.1561:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1562:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1563:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ coefficient

    read(unit=10,nml=single_particle_potential_parameters)

    write(*,*) "Using single particle atomic (1/r) potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624144408.1563:init_sp_potential
  !@+node:gcross.20090624144408.1564:Usp
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8) :: Usp

    real(kind=b8) :: ri

    ri = dot_product(x(slice,ip,:),x(slice,ip,:))**(-0.5_b8)

    Usp = coefficient * ri

  end function Usp_func
  !@nonl
  !@-node:gcross.20090624144408.1564:Usp
  !@+node:gcross.20090624144408.1565:gUsp
  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    real(kind=b8) :: r3i
    integer :: i

    do i = 1, np
      r3i = dot_product(x(slice,i,:),x(slice,i,:))**(-3.0_b8/2)
      gUsp(i,1) = -coefficient*x(slice,i,1)*r3i
      gUsp(i,2) = -coefficient*x(slice,i,2)*r3i
      gUsp(i,3) = -coefficient*x(slice,i,3)*r3i
    end do

  end function gUsp_func
  !@-node:gcross.20090624144408.1565:gUsp
  !@-others
  !@-node:gcross.20090624144408.1562:<< Subroutines >>
  !@nl

end module sp_atomic_potential
!@-node:gcross.20090624144408.1559:@thin sp_atomic_potential.f90
!@-leo
