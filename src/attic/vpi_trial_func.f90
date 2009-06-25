!@+leo-ver=4-thin
!@+node:gcross.20090623152316.91:@thin vpi_trial_func.f90
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
!@+node:gcross.20090623152316.115:eval_E_local
subroutine eval_E_local(np,ndim,grad_lnspf, lap_lnspf, grad_lnjas, lap_lnjas, U, E_l)
  integer, intent(in) :: np, ndim
  real(kind=b8), dimension( np , ndim ), intent(in) :: grad_lnjas 
  real(kind=b8), dimension( np , ndim ), intent(in) :: grad_lnspf 
  real(kind=b8), intent(in) :: lap_lnspf
  real(kind=b8), intent(in) :: lap_lnjas 
  real(kind=b8), intent(in) :: U
  real(kind=b8), intent(out) :: E_l

  real(kind=b8) :: t0,t1,T

  integer i,j

  t0 = 0.0_b8
!  do i = 1, np
!    do j = 1, ndim
!      t1  = grad_lnspf(i,j) +  grad_lnjas(i,j) 
!      t0 = t0 + t1*t1
!    end do
!  end do

  t0 = sum((grad_lnspf(:,:) +  grad_lnjas(:,:))**2)

  T = t0 + lap_lnspf + lap_lnjas
!  do i = 1, np
!    do j = 1, ndim
!      write(111,"(1a,2g18.6)") "#",grad_lnspf(i,j), grad_lnjas(i,j) 
!    end do
!  end do

!  write (111,"(7a18)") "#grad^2","lap_lnspf","lap_lnjas","T","U","E_l","lambda"
!  write (111,"(7g18.6)") t0,lap_lnspf,lap_lnjas,T,U,E_l,lambda
  E_l = U - lambda*T


end subroutine  eval_E_local
!@-node:gcross.20090623152316.115:eval_E_local
!@+node:gcross.20090623152316.100:numeric differentiation
!@+node:gcross.20090623152316.102:Jastrow function
function numeric_grad_lap_jas( x, xij2, sl, np, ndim, nslice, grad_lntfn, lap_lntfn, jfunc_params, jfunc ) result (y)
  integer, intent(in) :: sl, np, ndim, nslice
  real(kind=b8), dimension( nslice, np, ndim ) :: x
  real(kind=b8), dimension( nslice, np, np ) :: xij2
  real(kind=b8), dimension( np, ndim ), intent(out) :: grad_lntfn
  real(kind=b8), intent(out) :: lap_lntfn
  real(kind=b8) , dimension(:) :: jfunc_params

  interface
    function jfunc( x, xij2, params, sl, nslice, np, ndim) result(y)
      use kinds
      implicit none
      integer :: sl, nslice, np, ndim
      real(kind=b8), dimension(:) :: params
      real(kind=b8), dimension( nslice, np, ndim ), intent(in) :: x
      real(kind=b8), dimension( nslice, np, np ), intent(in) :: xij2
      real(kind=b8)  :: y
    end function jfunc
  end interface

  integer :: y

  real(kind=b8), dimension( nslice, np, ndim ) :: tx
  real(kind=b8), dimension( nslice, np, np ) :: txij2
  real(kind=b8), dimension( ndim ) :: fhi,flo
  real(kind=b8) :: f0
  integer :: i, k

  real(kind=b8), dimension( np, ndim ) :: tgrad_lntfn
  real(kind=b8) :: tlap_lntfn, tmp
  integer :: pass

  y = 1

  grad_lntfn = 0.0_b8
  lap_lntfn = 0.0_b8
  f0 = jfunc( x, xij2, jfunc_params, sl, nslice, np, ndim )
!  f0 = test_tfun( txij2, jfunc_params, sl, nslice, np )
  do i = 1, np
    do k = 1, ndim
      tx(sl,:,:) = x(sl,:,:)
      tx(sl,i,k) = x(sl,i,k) + ntol_eps
      txij2(sl,:,:) = xij2(sl,:,:)
      if(use_pbc) then
        call vpi_update_xij_pbc( txij2, tx, sl, sl, i, nslice, np, ndim )
      else
        call vpi_update_xij( txij2, tx, sl, sl, i, nslice, np, ndim )
      end if
      fhi(k) = jfunc( tx, txij2, jfunc_params, sl, nslice, np, ndim )

      txij2(sl,:,:) = xij2(sl,:,:)
      tx(sl,i,k) = x(sl,i,k) - ntol_eps
      if(use_pbc) then
        call vpi_update_xij_pbc( txij2, tx, sl, sl, i, nslice, np, ndim )
      else
        call vpi_update_xij( txij2, tx, sl, sl, i, nslice, np, ndim )
      end if
      flo(k) = jfunc( tx, txij2, jfunc_params, sl, nslice, np, ndim )
    end do
    grad_lntfn(i,:) = fhi(:) - flo(:)
    tmp = sum( fhi(:) + flo(:) ) -2.0*ndim*f0
    lap_lntfn = lap_lntfn + tmp
!    write(1001,"(i10, 12g20.12)") i, flo,f0,fhi
!    write(1002,"(i10, 12g20.12)") i, grad_lntfn(i,:)
!    write(1003,"(i10, 12g20.12)") i, tmp
  end do
  grad_lntfn(:,:) = grad_lntfn(:,:)/(2.0_b8*ntol_eps)
  lap_lntfn = lap_lntfn / ntol_eps**2
!  write(1004,"(12g20.12)") grad_lntfn
!  write(1005,"(12g20.12)") lap_lntfn

!  pass = grad_lap_lj_tfun( x, xij2, sl, N_PARTICLE, N_DIM, N_SLICE, tgrad_lntfn, tlap_lntfn, jfunc_params, jfunc )

!  write(1006,"(12g20.12)") grad_lntfn(:,:) - tgrad_lntfn(:,:)
!  write(1007,"(12g20.12)") lap_lntfn - tlap_lntfn


end function numeric_grad_lap_jas
!@-node:gcross.20090623152316.102:Jastrow function
!@+node:gcross.20090623152316.106:???
!function numeric_grad_lap(np,ndim,func,x,xij2,eps) result(grad_f)
!  integer :: np, ndim
!  real(kind=b8), dimension( np , ndim ) :: x
!end function
!@-node:gcross.20090623152316.106:???
!@-node:gcross.20090623152316.100:numeric differentiation
!@+node:gcross.20090623152316.103:Trial functions
!@+node:gcross.20090623152316.104:Single particle functions
!@+node:gcross.20090623152316.124:lattice
!@-node:gcross.20090623152316.124:lattice
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
!@+node:gcross.20090623152316.129:independent (w/ optimized differentiator)
function independent_tfun( x, xij2, params, sl, nslice, np, ndim ) result( y )
  use kinds
  implicit none
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( nslice , np, ndim ), intent(in) :: x
  real(kind=b8), dimension( nslice , np, np ), intent(in) :: xij2
  real(kind=b8), dimension(:) :: params
  real(kind=b8)  :: r2
  real(kind=b8)  :: y

  y = 0

end function independent_tfun

function grad_lap_independent_tfun( x, xij2, sl, np, ndim, nslice, grad_lntfn, lap_lntfn, jfunc_params, jfunc ) result (y)
  integer :: sl, np, ndim, nslice
  real(kind=b8), dimension( nslice, np, ndim ) :: x
  real(kind=b8), dimension( nslice, np, np ) :: xij2
  real(kind=b8), dimension( np, ndim ), intent(out) :: grad_lntfn 
  real(kind=b8), intent(out) :: lap_lntfn 
  real(kind=b8), dimension(:) :: jfunc_params
  interface
    function jfunc( x, xij2, params, sl, nslice, np, ndim) result(y)
      use kinds
      implicit none
      integer :: sl, nslice, np, ndim
      real(kind=b8), dimension(:) :: params
      real(kind=b8), dimension( nslice, np, ndim ), intent(in) :: x
      real(kind=b8), dimension( nslice, np, np ), intent(in) :: xij2
      real(kind=b8)  :: y
    end function jfunc
  end interface
  integer :: y


  real(kind=b8)  :: fi,fi2,ri,ri2,ri3,ri4
  real(kind=b8), dimension( ndim ) :: gtmp
  integer :: i, j

  grad_lntfn = 0.0_b8
  lap_lntfn = 0.0_b8
  y = 1
end function grad_lap_independent_tfun
!@-node:gcross.20090623152316.129:independent (w/ optimized differentiator)
!@+node:gcross.20090623152316.93:hard sphere (w/ optimized differentiator)
function hs_tfun( x, xij2, params, sl, nslice, np, ndim ) result( y )
  use kinds
  implicit none
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( nslice , np, ndim ), intent(in) :: x
  real(kind=b8), dimension( nslice , np, np ), intent(in) :: xij2
  real(kind=b8), dimension(:) :: params
  real(kind=b8)  :: r2
  real(kind=b8)  :: y

  integer :: i, j

  a_hs = params(1)
  a_hs2 = params(2)
  y = 0.0_b8
  do i = 1, np
    do j = i + 1, np
      r2 = xij2(sl,i,j)
      if ( r2 .gt. a_hs2 ) then
        y = y + log(1.0_b8 - a_hs/sqrt(r2))
      else
        y = -realbignumber
      end if
    end do
  end do

end function hs_tfun

function grad_lap_hs_tfun( x, xij2, sl, np, ndim, nslice, grad_lntfn, lap_lntfn, jfunc_params, jfunc ) result (y)
  integer :: sl, np, ndim, nslice
  real(kind=b8), dimension( nslice, np, ndim ) :: x
  real(kind=b8), dimension( nslice, np, np ) :: xij2
  real(kind=b8), dimension( np, ndim ), intent(out) :: grad_lntfn 
  real(kind=b8), intent(out) :: lap_lntfn 
  real(kind=b8), dimension(:) :: jfunc_params
  interface
    function jfunc( x, xij2, params, sl, nslice, np, ndim) result(y)
      use kinds
      implicit none
      integer :: sl, nslice, np, ndim
      real(kind=b8), dimension(:) :: params
      real(kind=b8), dimension( nslice, np, ndim ), intent(in) :: x
      real(kind=b8), dimension( nslice, np, np ), intent(in) :: xij2
      real(kind=b8)  :: y
    end function jfunc
  end interface
  integer :: y


  real(kind=b8)  :: fi,fi2,ri,ri2,ri3,ri4
  real(kind=b8), dimension( ndim ) :: gtmp
  integer :: i, j

  grad_lntfn = 0.0_b8
  lap_lntfn = 0.0_b8
  do i = 1, np
    do j = i + 1, np
      if ( xij2(sl,i,j) .gt. a_hs2 ) then
        ri2 = 1.0_b8/xij2(sl,i,j)
        ri = sqrt(ri2)
        ri3 = ri*ri2
        ri4 = ri2*ri2
        fi = 1.0_b8/(1.0_b8 - a_hs*ri)
        fi2 = fi*fi
        gtmp(:) = ri3*fi*(x(sl,i,:) - x(sl,j,:))
!        write (222,"(3g18.6)") gtmp(1),gtmp(2),gtmp(3)
        grad_lntfn(i,:) = grad_lntfn(i,:) + gtmp(:)
        grad_lntfn(j,:) = grad_lntfn(j,:) - gtmp(:)
        lap_lntfn = lap_lntfn - fi2*ri4 
      else
        y = -1
        return
      end if
    end do
  end do
  lap_lntfn = 2.0_b8*a_hs2*lap_lntfn
  grad_lntfn(:,:) = a_hs*grad_lntfn(:,:)
  y = 1
end function grad_lap_hs_tfun
!@-node:gcross.20090623152316.93:hard sphere (w/ optimized differentiator)
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
!@+node:gcross.20090623152316.96:lj (w/ optimized differentiator)
function lj_tfun( x, xij2, p, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8) , dimension( nslice , np, ndim ), intent(in) :: x
  real(kind=b8) , dimension( nslice , np, np ), intent(in) :: xij2
  real(kind=b8) , dimension(:) :: p
  real(kind=b8)  :: y
  real(kind=b8)  :: r5, r1

  integer :: i, j

  y = 0.0
  do i = 1, np
    do j = i + 1, np
      r1 = xij2(sl,i,j)**(-0.5_b8)
      r5 = r1*xij2(sl,i,j)**(-2)
      y = y - p_ljc5*r5 - p_ljc1*r1
    end do
  end do

end function lj_tfun

function grad_lap_lj_tfun( x, xij2, sl, np, ndim, nslice, grad_lntfn, lap_lntfn, jfunc_params, jfunc ) result (y)
  integer, intent(in) :: sl, np, ndim, nslice
  real(kind=b8), dimension( nslice, np, ndim ) :: x
  real(kind=b8), dimension( nslice, np, np ) :: xij2
  real(kind=b8), dimension( np, ndim ), intent(out) :: grad_lntfn
  real(kind=b8), intent(out) :: lap_lntfn
  real(kind=b8) , dimension(:) :: jfunc_params

  interface
    function jfunc( xij2, params, sl, nslice, np) result(y)
      use kinds
      implicit none
      integer :: sl, nslice, np
      real(kind=b8), dimension(:) :: params
      real(kind=b8), dimension( nslice, np, np ) :: xij2
      real(kind=b8)  :: y
    end function jfunc
  end interface

  real(kind=b8)  :: r,ri,ri2,ri3,ri4,ri6,ri7,ri8,ri12
  real(kind=b8), dimension( ndim ) :: gtmp
  real(kind=b8) :: p5sq,p1sq,tmp
  integer :: y
  integer :: i, j

  p5sq = p_ljc5*p_ljc5
  p1sq = p_ljc1*p_ljc1

  grad_lntfn = 0
  lap_lntfn = 0
  do i = 1, np
    do j = i + 1, np
      r = sqrt(xij2(sl,i,j))
      ri2 = 1.0_b8/xij2(sl,i,j)
      ri = 1.0/r
      ri3 = ri*ri2
      ri4 = ri2*ri2
      ri6 = ri3*ri3
      ri8 = ri4*ri4
      ri7 = ri*ri6
      ri12 = ri6*ri6

      gtmp(:) =  (5.0*p_ljc5*ri7 + p_ljc1*ri3)*(x(sl,i,:) - x(sl,j,:))
      grad_lntfn(i,:) = grad_lntfn(i,:) + gtmp(:)
      grad_lntfn(j,:) = grad_lntfn(j,:) - gtmp(:)
      tmp = 25.0*p5sq*ri12 + 10.0*p_ljc5*(p_ljc1 - 2*r)*ri8 + p1sq*ri4 
      lap_lntfn = lap_lntfn + tmp
    end do
  end do
  y = 1
end function grad_lap_lj_tfun
!@-node:gcross.20090623152316.96:lj (w/ optimized differentiator)
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
!@+node:gcross.20090623152316.98:charge
function charge_tfun( x, xij2, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8) , dimension( nslice , np, ndim ), intent(in) :: x
  real(kind=b8) , dimension( nslice , np, np ), intent(in) :: xij2
  real(kind=b8)  :: y

  real(kind=b8)  :: r2,r1

  integer :: i, j

  y = 0.0
  do i = 1, np
    do j = i + 1, np
      r2 = xij2(sl,i,j)
      r1 = sqrt(r2)
      y = y - ( p_cc0*r1 + p_cc1*r2 )/( 1 + p_cc2*r1 )
    end do
  end do

end function charge_tfun
!@-node:gcross.20090623152316.98:charge
!@-node:gcross.20090623152316.92:Jastrow functions
!@-node:gcross.20090623152316.103:Trial functions
!@-others

end module vpi_trial_func
!@-node:gcross.20090623152316.91:@thin vpi_trial_func.f90
!@-leo
