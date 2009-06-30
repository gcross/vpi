!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1672:@thin tb_lj_potential.f90
!@@language fortran90
module vpi_two_body_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1673:<< Imported modules >>
  use kinds
  !@-node:gcross.20090624144408.1673:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1674:<< Variables >>
  real (kind=b8), private :: coefficient = 2352.5_b8
  real (kind=b8), private :: length_scale = 0.05_b8
  real (kind=b8), private :: length_scale_squared
  !@-node:gcross.20090624144408.1674:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1675:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1676:init_tb_potential
  subroutine init_tb_potential ()
    namelist /two_body_potential_parameters/ coefficient, length_scale

    read(unit=10,nml=two_body_potential_parameters)

    length_scale_squared = length_scale * length_scale

    write(*,*) "Using two-body L.J. potential with"
    write(*,nml=two_body_potential_parameters)
  end subroutine init_tb_potential
  !@-node:gcross.20090624144408.1676:init_tb_potential
  !@+node:gcross.20090624144408.1677:Uij
  function Uij_func( x, xij2, slice, ip, nslice, np, ndim, acc_flag ) result ( Uij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, ip, nslice, np, ndim
    logical :: acc_flag

    real(kind=b8) :: Uij

    real(kind=b8) :: t3,t6

    acc_flag = .true.

    t6 = length_scale_squared**6 * (sum( xij2(slice,ip,1:ip-1)**(-6) ) + sum( xij2(slice,ip,ip+1:np)**(-6) ))
    t3 = length_scale_squared**3 * (sum( xij2(slice,ip,1:ip-1)**(-3) ) + sum( xij2(slice,ip,ip+1:np)**(-3) ))

    Uij = coefficient * ( t6 - t3 )/2.0_b8

  end function Uij_func
  !@-node:gcross.20090624144408.1677:Uij
  !@+node:gcross.20090624144408.1678:gUij
  !> f := e*((r/s)^(-12) - (r/s)^(-6))/2;
  !> v := [r, theta, phi]:
  !> simplify(grad(f,v,coords=spherical)[1]/r);
  !                                              6      6    6
  !                                         3 e s  (-2 s  + r )
  !                                         -------------------
  !                                                  14
  !                                                 r


  function gUij_func( x, xij2, slice, nslice, np, ndim ) result ( gUij )
    real(kind=b8), dimension ( nslice, np , ndim ) :: x
    real(kind=b8), dimension ( nslice, np , np ) :: xij2
    integer :: slice, nslice, np, ndim

    real(kind = b8), dimension ( np , ndim ) :: gUij

    real(kind = b8), dimension ( ndim ) :: t_gUij
    real(kind = b8) :: rij, tij
    real(kind = b8) :: gc1,gc2
    integer :: i,j

    gc1 =  -6.0_b8*length_scale_squared**6
    gc2  =  3.0_b8*length_scale_squared**3

    gUij = 0

    do i = 1, np
      do j = i+1, np
        rij =  xij2(slice,i,j)
        tij = gc1*rij**(-7) + gc2*rij**(-4)
        t_gUij(:) = ( x(slice,i,:) - x(slice,j,:) ) * tij
        gUij(i,:) = gUij(i,:) + coefficient*t_gUij(:)
        gUij(j,:) = gUij(j,:) - coefficient*t_gUij(:)
      end do
    end do

  end function gUij_func
  !@-node:gcross.20090624144408.1678:gUij
  !@-others
  !@-node:gcross.20090624144408.1675:<< Subroutines >>
  !@nl

end module vpi_two_body_potential
!@-node:gcross.20090624144408.1672:@thin tb_lj_potential.f90
!@-leo
