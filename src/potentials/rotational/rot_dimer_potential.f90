!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1758:@thin rot_dimer_potential.f90
!@@language fortran90
module rot_dimer_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1759:<< Imported modules >>
  use kinds
  use vpi_defines
  !@-node:gcross.20090624144408.1759:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1760:<< Variables >>
  real (kind=b8), private :: coefficient = 1e-4_b8

  !@-node:gcross.20090624144408.1760:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1761:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1762:init_rot_potential
  subroutine init_rot_potential ()
    namelist /rotational_potential_parameters/ rotational

    read(unit=10,nml=rotational_potential_parameters)

    write(*,*) "Using dimer rotational potential with"
    write(*,nml=rotational_potential_parameters)
  end subroutine init_rot_potential
  !@-node:gcross.20090624144408.1762:init_rot_potential
  !@+node:gcross.20090624144408.1763:RUij
  function Uij_rot_func( x, x_rot, xij2, slice, ip, nslice, np, ndim ) result ( Uij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, ip, nslice, np, ndim

    real(kind=b8) :: Uij

    real(kind=b8) :: h1,h2
    real(kind=b8) :: t3,t6,costh
    real(kind=b8), dimension (3) :: dr
    integer :: i

    t3 = 0d0
    do i=1,ip-1
      costh = dot_product(x_rot(slice,ip,:),x_rot(slice,i,:))
      dr = x(slice,ip,:)-x(slice,i,:)
      h1 = dot_product(dr,x_rot(slice,i ,:))
      h2 = dot_product(dr,x_rot(slice,ip,:))
      t3 = t3 + (costh-3d0*h1*h2/xij2(slice,ip,i))*xij2(slice,ip,i)**(-1.5d0)
    end do
    do i=ip+1,np
      costh = dot_product(x_rot(slice,ip,:),x_rot(slice,i,:))
      dr = x(slice,ip,:)-x(slice,i,:)
      h1 = dot_product(dr,x_rot(slice,i ,:))
      h2 = dot_product(dr,x_rot(slice,ip,:))
      t3 = t3 + (costh-3d0*h1*h2/xij2(slice,ip,i))*xij2(slice,ip,i)**(-1.5d0)
    end do

    Uij = coefficient*t3

  end function Uij_rot_func

  !@-node:gcross.20090624144408.1763:RUij
  !@+node:gcross.20090624144408.1764:gRUij
  function gUij_rot_func( x, x_rot, xij2, slice, nslice, np, ndim ) result ( gUij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , N_DIM_ROT ) :: x_rot
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, nslice, np, ndim

    real(kind = b8), dimension ( np , ndim ) :: gUij

    gUij = 0
  end function gUij_rot_func
  !@-node:gcross.20090624144408.1764:gRUij
  !@-others
  !@-node:gcross.20090624144408.1761:<< Subroutines >>
  !@nl

end module rot_dimer_potential
!@-node:gcross.20090624144408.1758:@thin rot_dimer_potential.f90
!@-leo
