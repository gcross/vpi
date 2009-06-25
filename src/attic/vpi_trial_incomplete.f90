!@+leo-ver=4-thin
!@+node:gcross.20090623152316.91:@thin vpi_trial_incomplete.f90
!@@language fortran90
module vpi_trial_func

  !@  << Import modules >>
  !@+node:gcross.20090623152316.128:<< Import modules >>
    use vpi_defines
    use vpi_xij
    use vpi_potential
  !@-node:gcross.20090623152316.128:<< Import modules >>
  !@nl

  implicit none

  contains

!@+others
!@+node:gcross.20090623152316.103:Trial functions
!@+node:gcross.20090623152316.104:Single particle functions
!@+node:gcross.20090623152316.127:fcc (?)
function fcc_tfun( x, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( nslice, np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np, ndim ) :: psi_t

  psi_t(:,:)  = (x(sl,:,:)-p_lat_r0(:,:))**2
  y = -p_lat_w*sum(psi_t)/2.0_b8

end function fcc_tfun
!@-node:gcross.20090623152316.127:fcc (?)
!@+node:gcross.20090623152316.120:(functions w/ incomplete interfaces)
!@+node:gcross.20090623152316.121:annulus
function annulus_tfun( x, np, ndim ) result( y )
  integer :: np, ndim
  real(kind=b8), dimension( np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np ) :: psi_t

  psi_t(:)  = -( p_hox*x(:,1)**2 + p_hoy*x(:,2)**2 + p_hoz*x(:,3)**2 )/2.0 -&
              p_aw/( x(:,1)**2 + x(:,3)**2 )
  y = sum(psi_t)

end function annulus_tfun
!@-node:gcross.20090623152316.121:annulus
!@+node:gcross.20090623152316.122:3D dwell
function tfunc_3D_dwell( x ) result( y )
  real(kind=b8), dimension( : , : ) :: x
  real(kind=b8) :: y

  real :: f0, f1, f2
  integer :: i, j
  integer :: np

  f0 = 0.758993
  f1 = 2.31308
  f2 = 0.807643   
  y = 1.0
  np = size(x,1)

  do i = 1, np
    y =  y * (f0*exp(-f1*(x(i,3)-f2)**2) + f0*exp(-f1*(x(i,3)+f2)**2))
    y =  y * exp( -( x(i,1)**2 + x(i,2)**2 ) )
  end do

end function tfunc_3D_dwell
!@nonl
!@-node:gcross.20090623152316.122:3D dwell
!@+node:gcross.20090623152316.125:poor
function tfunc_poor( x ) result( y )
  real(kind=b8), dimension( : , : ) :: x
  real(kind=b8) :: y

  integer :: i, j

  y = 1.0

  do i = 1, size(x,1)
    do j = 1, size(x,2)
      y = y * ( exp( -abs( x(i,j) - 2 ) ) + exp( -abs( x(i,j) + 1 ) ) )
    end do
  end do
end function tfunc_poor
!@-node:gcross.20090623152316.125:poor
!@+node:gcross.20090623152316.126:lattice, old version
function lattice_tfun0( x, np, ndim ) result( y )
  integer :: np, ndim
  real(kind=b8), dimension( np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np ) :: psi_t

  integer :: i, j

  psi_t(:)  = -( p_hox*x(:,1)**2 + p_hoy*x(:,2)**2 + p_hoz*x(:,3)**2 )/2.0 + p_lwx*sin(x(:,1))**2 + p_lwz*sin(x(:,3))**2
  y = sum(psi_t)

end function lattice_tfun0
!@nonl
!@-node:gcross.20090623152316.126:lattice, old version
!@-node:gcross.20090623152316.120:(functions w/ incomplete interfaces)
!@-node:gcross.20090623152316.104:Single particle functions
!@+node:gcross.20090623152316.92:Jastrow functions
!@@language fortran90
!@+node:gcross.20090623152316.94:test
function test_tfun( x, xij2, p, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8) , dimension( nslice , np, ndim ), intent(in) :: x
  real(kind=b8) , dimension( nslice , np, np ), intent(in) :: xij2
  real(kind=b8) , dimension(:) :: p
  real(kind=b8)  :: y
  real(kind=b8)  :: r5, r1

  integer :: i, j

  y = 0.0_b8
  do i = 1, np
    do j = i + 1, np
      r1 = xij2(sl,i,j)
      y = y - xij2(sl,i,j)
    end do
  end do

end function test_tfun
!@-node:gcross.20090623152316.94:test
!@+node:gcross.20090623152316.95:dimer
function dimer_tfun( x, xij2, p, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8) , dimension( nslice , np, ndim ), intent(in) :: x
  real(kind=b8) , dimension( nslice , np, np ), intent(in) :: xij2
  real(kind=b8) , dimension(:) :: p
  real(kind=b8)  :: y

  integer :: i
  logical :: acc_flag

  y = 0.0
  do i = 1, np
    y = y - p(1)*vpi_Uij_z_polarized_dimer( x, xij2, sl, i, nslice, np, ndim, acc_flag )
  end do

end function dimer_tfun
!@-node:gcross.20090623152316.95:dimer
!@+node:gcross.20090623152316.97:aziz
function aziz_tfun( x, xij2, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8) , dimension( nslice , np, ndim ), intent(in) :: x
  real(kind=b8) , dimension( nslice , np, np ), intent(in) :: xij2
  real(kind=b8)  :: y
  real(kind=b8)  :: r5, r1

  real(kind=b8)  :: alpha, beta

  integer :: i, j

  alpha = 19.0_b8
  beta = 0.12_b8

  y = 0.0
  do i = 1, np
    do j = i + 1, np
      r1 = sqrt(xij2(sl,i,j))
      r5 = r1*xij2(sl,i,j)**2
      y = y - alpha/(1.0_b8+beta*r5)
    end do
  end do
  y = y*0.5_b8

end function aziz_tfun
!@-node:gcross.20090623152316.97:aziz
!@-node:gcross.20090623152316.92:Jastrow functions
!@-node:gcross.20090623152316.103:Trial functions
!@-others

end module vpi_trial_func
!@-node:gcross.20090623152316.91:@thin vpi_trial_incomplete.f90
!@-leo
