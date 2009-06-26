!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1441:@thin sp_morse_oscillator_potential.f90
!@@language fortran90
module sp_morse_oscillator_potential

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1442:<< Imported modules >>
  use kinds
  use vpi_defines
  !@-node:gcross.20090624144408.1442:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.1443:<< Variables >>
  !Morse oscillator parameters where
  !H = -0.5T + p_MO_De*( exp(-2*p_MO_a*(r-p_MO_r0)) -2*exp(-p_MO_a*(r-p_MO_r0)) )
  real(kind=b8), private :: p_MO_a = 10.0_b8
  !real(kind=b8), parameter :: p_MO_ap = 1.010971279201_b8 ! for evaluation of overlap integral
  real(kind=b8), private :: p_MO_ap = 10.0_b8 ! for evaluation of overlap integral
  real(kind=b8), private :: p_MO_r0 = 1.0_b8
  real(kind=b8), private :: p_MO_De = 50.0_b8
  real(kind=b8), private :: p_MO_Dep = 50.0_b8! for evaluation of overlap integral
  !@-node:gcross.20090624144408.1443:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1444:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.1445:init_sp_potential
  subroutine init_sp_potential ()
    namelist /single_particle_potential_parameters/ p_MO_a, p_MO_ap, p_MO_r0, p_MO_De, p_MO_Dep

    read(unit=10,nml=single_particle_potential_parameters)

    write(*,*) "Using single particle Morse oscillator potential with"
    write(*,nml=single_particle_potential_parameters)
  end subroutine init_sp_potential
  !@nonl
  !@-node:gcross.20090624144408.1445:init_sp_potential
  !@+node:gcross.20090624144408.1446:Usp
  !dimensionless 3D morse oscilator H = -0.5T +De[exp(-2a(r-r0))-2*exp(-a(r-r0))]
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer :: nslice, np, ndim
    integer :: slice, ip
    real(kind=b8), dimension ( nslice, np , ndim ) :: x

    real(kind=b8) :: Usp
    real(kind=b8) :: r,t1,t2

    r = sqrt(x(slice,ip,1)**2 + x(slice,ip,2)**2 + x(slice,ip,3)**2)

    if(slice .le. CSLICE) then
      t1 = exp(-2.0_b8*p_MO_a*(r-p_MO_r0))
      t2 = exp(-p_MO_a*(r-p_MO_r0))
      Usp = p_MO_De*(t1-2.0_b8*t2)
    else
      t1 = exp(-2.0_b8*p_MO_ap*(r-p_MO_r0))
      t2 = exp(-p_MO_ap*(r-p_MO_r0))
      Usp = p_MO_Dep*(t1-2.0_b8*t2)
    end if

  end function Usp_func
  !@-node:gcross.20090624144408.1446:Usp
  !@+node:gcross.20090624144408.1447:gUsp
  !                                   2    2    2 1/2
  !-2 p_MO_De p_MO_a x exp(p_MO_a (-(x  + y  + z )    + p_MO_r0))
  !
  !                    2    2    2 1/2                    /   2    2    2 1/2
  !    (exp(p_MO_a (-(x  + y  + z )    + p_MO_r0)) - 1)  /  (x  + y  + z )

  function gUsp_func( x, slice, nslice, np, ndim ) result ( gUsp )
    real(kind=b8), dimension ( nslice, np , ndim ) ::  x
    integer :: slice, nslice, np, ndim

    real(kind=b8), dimension ( np, ndim ) :: gUsp

    real(kind=b8), dimension ( np ) :: r,t1,t2

    r(:) = sqrt(x(slice,:,1)**2 + x(slice,:,2)**2 + x(slice,:,3)**2)
    if(slice .le. CSLICE) then
      t1(:) = exp(-2.0_b8*p_MO_a*(r(:)-p_MO_r0))
      t2(:) = exp(-p_MO_a*(r(:)-p_MO_r0))
      gUsp(:,1) = 2.0_b8*p_MO_De*p_MO_a*(t2(:) - t1(:))*x(slice,:,1)/r(:)
      gUsp(:,2) = 2.0_b8*p_MO_De*p_MO_a*(t2(:) - t1(:))*x(slice,:,2)/r(:)
      gUsp(:,3) = 2.0_b8*p_MO_De*p_MO_a*(t2(:) - t1(:))*x(slice,:,3)/r(:)
    else
      t1(:) = exp(-2.0_b8*p_MO_ap*(r(:)-p_MO_r0))
      t2(:) = exp(-p_MO_ap*(r(:)-p_MO_r0))
      gUsp(:,1) = 2.0_b8*p_MO_Dep*p_MO_ap*(t2(:) - t1(:))*x(slice,:,1)/r(:)
      gUsp(:,2) = 2.0_b8*p_MO_Dep*p_MO_ap*(t2(:) - t1(:))*x(slice,:,2)/r(:)
      gUsp(:,3) = 2.0_b8*p_MO_Dep*p_MO_ap*(t2(:) - t1(:))*x(slice,:,3)/r(:)
    end if

  end function gUsp_func
  !@-node:gcross.20090624144408.1447:gUsp
  !@-others
  !@-node:gcross.20090624144408.1444:<< Subroutines >>
  !@nl

end module sp_morse_oscillator_potential
!@-node:gcross.20090624144408.1441:@thin sp_morse_oscillator_potential.f90
!@-leo
