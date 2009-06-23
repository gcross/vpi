module vpi_trial_func 

  use vpi_defines
  use vpi_xij
  use vpi_potential
  implicit none

  contains


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

function numeric_grad_lap_spf( x, sl, np, ndim, nslice, grad_lntfn, lap_lntfn, sp_param, spf_func ) result( y )
  real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
  real(kind=b8), dimension( : ) :: sp_param
  real(kind=b8),intent(out) :: lap_lntfn 
  integer, intent(in) :: np, ndim, nslice, sl

  interface
    function spf_func( x, sl, param, nslice, np, ndim ) result( y )
      use kinds
      implicit none
      integer :: sl, nslice, np, ndim
      real(kind=b8), dimension( nslice, np , ndim ) :: x
      real(kind=b8), dimension( : ) :: param
      real(kind=b8) :: y
    end function spf_func
  end interface
  integer :: y

  real(kind=b8), dimension( nslice, np, ndim ) :: tx
  real(kind=b8), dimension( ndim ) :: fhi,flo,f
  integer :: i, k

  y = 0
  tx(sl,:,:) = x(sl,:,:)

  grad_lntfn = 0
  lap_lntfn = 0
  do i = 1, np
    do k = 1, ndim
      f(k) = spf_func( x, sl, sp_param, nslice, np, ndim )
      tx(sl,i,k) = x(sl,i,k) + ntol_eps
      fhi(k) = spf_func( tx, sl, sp_param, nslice, np, ndim )
      tx(sl,i,k) = x(sl,i,k) - ntol_eps
      flo(k) = spf_func( tx, sl, sp_param, nslice, np, ndim )
      tx(sl,i,k) = x(sl,i,k)
    end do
#ifdef DEBUG
    write(1000,*) fhi
    write(1000,*) f
    write(1000,*) flo
#endif
    grad_lntfn(i,:) = (fhi(:)-flo(:))/(2.0_b8*ntol_eps)
    lap_lntfn = lap_lntfn + sum((fhi(:)-2.0_b8*f(:)+flo(:)))/ntol_eps**2
  end do

!  write(1000,*) grad_lntfn
!  write(1000,*) lap_lntfn

  y = 1

end function numeric_grad_lap_spf

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

!function numeric_grad_lap(np,ndim,func,x,xij2,eps) result(grad_f)
!  integer :: np, ndim
!  real(kind=b8), dimension( np , ndim ) :: x
!end function

function test_tfun( xij2, p, sl, nslice, np ) result( y )
  integer :: sl, nslice, np
  real(kind=b8) , dimension( nslice , np, np ) :: xij2
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




function aziz_tfun( xij2, sl, nslice, np ) result( y )
  integer :: sl, nslice, np
  real(kind=b8) , dimension( nslice , np, np ) :: xij2
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

function charge_tfun( xij2, sl, nslice, np ) result( y )
  integer :: sl, nslice, np
  real(kind=b8) , dimension( nslice , np, np ) :: xij2
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

function atom_tfun( x, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( nslice, np , ndim ) :: x
  
  real(kind=b8) :: y

  real(kind=b8)  :: r2,r1
  integer :: i

  y = 0.0
  do i = 1, np
    r2 = dot_product(x(sl,i,:),x(sl,i,:))
    r1 = sqrt(r2)
    y = y - ( p_ac0*r1 + p_ac1*r2 )/( 1 + p_ac1*r1 )
  end do

end function atom_tfun

function grad_lap_atom_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
  real(kind=b8), dimension( : , : , : ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn
  real(kind=b8),intent(out) :: lap_lntfn
  integer, intent(in) :: np, ndim, nslice, slice
  integer :: y

  real, dimension( np ) :: t_grad

  real(kind=b8), dimension( np ) :: r
  real(kind=b8), dimension( np ) :: r2
  real(kind=b8), dimension( np ) :: r3

  r2(:) = x(slice,:,1)**2 + x(slice,:,2)**2 + x(slice,:,3)**2
  r(:) = sqrt(r2(:))
  r3(:) = r(:)*r2(:)
  t_grad(:) = -(p_ac0+p_ac1*(2*r(:)+p_ac1*r2(:)))/(r(:)*(1+p_ac1*r2(:))**2)
  grad_lntfn(:,1) = x(slice,:,1)*t_grad(:)
  grad_lntfn(:,2) = x(slice,:,2)*t_grad(:)
  grad_lntfn(:,3) = x(slice,:,3)*t_grad(:)

  lap_lntfn = sum((-2.0_b8*(p_ac0 + p_ac1*(3.0_b8*r(:) + 3.0_b8*p_ac1*r2(:) + p_ac1**2*r3(:) )))/(r(:)*(1.0_b8+p_ac1*r2(:))**3))

  y = 1

end function grad_lap_atom_tfun


function He3_fn_plus_tfun( x, islice, nslice, np, ndim ) result( y )
  integer :: islice, nslice, np, ndim
  real(kind=b8), dimension( : , : , : ) :: x
  
  real(kind=b8) :: y
  real(kind=b8)  :: r1,r1_sq,r2,r2_sq
  real(kind=b8)  :: rc
  integer :: i, j

  r1_sq =  dot_product(x(islice,1,:),x(islice,1,:))
  r1 = sqrt(r1_sq)
  r2_sq =  dot_product(x(islice,2,:),x(islice,2,:))
  r2 = sqrt(r2_sq)

  rc = r1-r2

  if( rc .gt. 0 ) then
    y = log(r1-r2)
  else 
    y = -realbignumber
  end if

end function He3_fn_plus_tfun


function grad_lap_He3_plus_tfun( x, slice, grad_lntfn, lap_lntfn, np, ndim ) result( y )
  real(kind=b8), dimension( : , : , : ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
  real(kind=b8),intent(out) :: lap_lntfn 
  integer, intent(in) :: np, ndim, slice
  integer :: y

  integer :: i, j

  grad_lntfn(:,1) = -2.0_b8*p_hox*x(slice,:,1)
  grad_lntfn(:,2) = -2.0_b8*p_hoy*x(slice,:,2)
  grad_lntfn(:,3) = -2.0_b8*p_hoz*x(slice,:,3)

  lap_lntfn = -2.0_b8*(p_hox+p_hoy+p_hoz)*np 
  y = 1

end function grad_lap_He3_plus_tfun

function He3_fn_neg_tfun( x, islice, nslice, np, ndim ) result( y )
  integer :: islice, nslice, np, ndim
  real(kind=b8), dimension( :,:,: ) :: x
  
  real(kind=b8) :: y
  real(kind=b8)  :: r1,r1_sq,r2,r2_sq
  real(kind=b8)  :: rc
  integer :: i, j

  r1_sq =  dot_product(x(islice,1,:),x(islice,1,:))
  r1 = sqrt(r1_sq)
  r2_sq =  dot_product(x(islice,2,:),x(islice,2,:))
  r2 = sqrt(r2_sq)

  rc = r1-r2

  if( rc .gt. 0 ) then
    y = log(r1-r2)
  else 
    y = -realbignumber
  end if

end function He3_fn_neg_tfun

function ho_tfun( x, sl, param, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( : ) :: param
  real(kind=b8), dimension( nslice, np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np ) :: psi_t


  psi_t(:)  = -( p_hox*x(sl,:,1)**2 + p_hoy*x(sl,:,2)**2 + p_hoz*x(sl,:,3)**2 )/2.0_b8
  y = sum(psi_t)

end function ho_tfun

function hoz_tfun( x, sl, param, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( : ) :: param
  real(kind=b8), dimension( nslice, np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np ) :: psi_t


  psi_t(:)  = -p_hoz*x(sl,:,3)**2 -p_hoy*x(sl,:,3)**4
  y = sum(psi_t)

end function hoz_tfun

function box_tfun( x, np, ndim ) result( y )
  integer :: np, ndim
  real(kind=b8), dimension( np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np ) :: psi_t

  integer :: i, j

  psi_t(:)  = sum(log(cos( M_PI*x(:,:)/(2.0_b8*abox) )))
  y = sum(psi_t)

end function box_tfun

function grad_lap_box_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
  real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
  real(kind=b8),intent(out) :: lap_lntfn 
  integer, intent(in) :: np, ndim, nslice, slice
  integer :: y

  integer :: i, j

  grad_lntfn(:,1) = -tan(x(slice,:,1))
  grad_lntfn(:,2) = -tan(x(slice,:,2))
  grad_lntfn(:,3) = -tan(x(slice,:,3))

  lap_lntfn = -3.0*np - sum(grad_lntfn(:,:)**2)
  y = 1

end function grad_lap_box_tfun

function grad_lap_ho_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
  real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
  real(kind=b8),intent(out) :: lap_lntfn 
  integer, intent(in) :: np, ndim, nslice, slice
  integer :: y

  grad_lntfn(:,1) = -p_hox*x(slice,:,1)
  grad_lntfn(:,2) = -p_hoy*x(slice,:,2)
  grad_lntfn(:,3) = -p_hoz*x(slice,:,3)

  lap_lntfn = -(p_hox+p_hoy+p_hoz)*dble(np) 
  y = 1

end function grad_lap_ho_tfun

function morse_tfun2( x, np, ndim ) result( y )
  integer :: np, ndim
  real(kind=b8), dimension( np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np ) :: psi_t

  real(kind=b8), dimension( np ) :: r

  integer :: i, j


  r(:) = sqrt(x(:,1)**2 + x(:,2)**2 + x(:,3)**2)
  psi_t(:)  = -p_MO_vpa*(r-p_MO_vpb)**2
  y = sum(psi_t)
end function morse_tfun2

! 
!> q := -a*(x-x0)**2;                                        
!                                                                      2
!                                                      q := -a (r - x0)
!
!>  simplify((grad(q,[r,theta,phi],coords=spherical)[1]));   
!                                                        -2 a (r - x0)
!
!>  simplify((laplacian(q,[r,theta,phi],coords=spherical))); 
!                                                         a (3 r - 2 x0)
!                                                      -2 --------------
!                                                               r
!
function grad_lap_morse_tfun2( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
  real(kind=b8), dimension( : , : , : ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
  real(kind=b8),intent(out) :: lap_lntfn 
  integer, intent(in) :: np, ndim, nslice, slice
  integer :: y 

  real, dimension( np ) :: t_grad
  real(kind=b8), dimension( np ) :: r

  r(:) = sqrt(x(slice,:,1)**2 + x(slice,:,2)**2 + x(slice,:,3)**2)
  t_grad(:) = -2.0_b8*p_MO_vpa*(1 - p_MO_vpb/r(:))
  grad_lntfn(:,1) = x(slice,:,1)*t_grad(:)
  grad_lntfn(:,2) = x(slice,:,2)*t_grad(:)
  grad_lntfn(:,3) = x(slice,:,3)*t_grad(:)

  lap_lntfn = -2.0_b8*p_MO_vpa*sum(3.0_b8 - 2.0_b8*p_MO_vpb/r(:))

  y = 1

end function grad_lap_morse_tfun2

function morse_tfun( x, sl,nslice,np, ndim ) result( y )
  integer :: np, ndim, nslice,sl
  real(kind=b8), dimension( nslice, np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np ) :: psi_t

  real(kind=b8), dimension( np ) :: r

  integer :: i, j

!  print *,"IN morse_tfun"

  r(:) = sqrt(x(sl,:,1)**2 + x(sl,:,2)**2 + x(sl,:,3)**2)
!  print *,r(:)
  psi_t(:)  = -( p_MO_vpa*r(:) +  p_MO_vpb/r(:)**3 )
  y = sum(psi_t)
!  print *,y

end function morse_tfun


!> q;
!                                           b
!                                  -a rr - ---
!                                            3
!                                          rr
!
!> simplify((grad(q,[rr,theta,phi],coords=spherical)[1]));
!                                       4
!                                   a rr  - 3 b
!                                 - -----------
!                                         4
!                                       rr
!
!> simplify((laplacian(q,[rr,theta,phi],coords=spherical)));
!                                       4
!                                   a rr  + 3 b
!                                -2 -----------
!                                         5
!                                       rr
function grad_lap_morse_tfun( x, slice, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
  real(kind=b8), dimension( : , : , : ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
  real(kind=b8),intent(out) :: lap_lntfn 
  integer, intent(in) :: np, ndim, nslice, slice
  integer :: y 

  real, dimension( np ) :: t_grad

  real(kind=b8), dimension( np ) :: r
  real(kind=b8), dimension( np ) :: r3
  real(kind=b8), dimension( np ) :: r5

  integer :: i, j

  r(:) = sqrt(x(slice,:,1)**2 + x(slice,:,2)**2 + x(slice,:,3)**2)
  r3(:) = r(:)**3
  r5(:) = r3(:)*r(:)**2
  t_grad(:) = -p_MO_vpa/r(:) + 3*p_MO_vpb/r5(:)
  grad_lntfn(:,1) = x(slice,:,1)*t_grad(:)
  grad_lntfn(:,2) = x(slice,:,2)*t_grad(:)
  grad_lntfn(:,3) = x(slice,:,3)*t_grad(:)

  lap_lntfn = -2.0_b8*sum( p_MO_vpa/r(:) + 3.0_b8*p_MO_vpb/r5(:) )

  y = 1

end function grad_lap_morse_tfun

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

function annulus_tfun( x, np, ndim ) result( y )
  integer :: np, ndim
  real(kind=b8), dimension( np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np ) :: psi_t

  psi_t(:)  = -( p_hox*x(:,1)**2 + p_hoy*x(:,2)**2 + p_hoz*x(:,3)**2 )/2.0 -&
              p_aw/( x(:,1)**2 + x(:,3)**2 )
  y = sum(psi_t)

end function annulus_tfun

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

function dw_tfun( x, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( nslice, np , ndim ) :: x
  real(kind=b8) :: y
  real(kind=b8), dimension( np ) :: psi_t

  y = 0

  psi_t(:)  = log( exp( -( p_hox*x(sl,:,1)**2 + p_hoy*x(sl,:,2)**2 )/2.0 ) * &
                 ( p_dw_f0*exp(-p_dw_f1*(x(sl,:,3)-p_dw_f2)**2) + &
                 p_dw_f0*exp(-p_dw_f1*(x(sl,:,3)+p_dw_f2)**2) - &
                 p_dw_f3*exp(-p_dw_f4*(x(sl,:,3)**2))) )

!  psi_t(:)  = -( x(:,1)**2 + x(:,2)**2 + x(:,3)**2 )/2.0
  y = sum(psi_t)

end function dw_tfun

function grad_dw_tfun( x, ndim, grad_lntfn ) result( y )
  real(kind=b8), dimension( ndim ), intent(in) :: x
  real(kind=b8), dimension( ndim ), intent(out) :: grad_lntfn 
  integer, intent(in) :: ndim 
  integer :: y

  real(kind=b8), dimension( ndim ) :: tx
  real(kind=b8), dimension( ndim ) :: fhi,flo
  integer :: k

  grad_lntfn = 0
  tx(:) = x(:)
  do k = 1, ndim
    tx(k) = x(k) + ntol_eps
    fhi(k)  = log( exp( -( p_hox*tx(1)**2 + p_hoy*tx(2)**2 )/2.0 ) * ( p_dw_f0*exp(-p_dw_f1*(tx(3)-p_dw_f2)**2) + p_dw_f0*exp(-p_dw_f1*(tx(3)+p_dw_f2)**2) + p_dw_f3*exp(-p_dw_f4*tx(3)**2)))
  
    tx(k) = x(k) - ntol_eps
    flo(k)  = log( exp( -( p_hox*tx(1)**2 + p_hoy*tx(2)**2 )/2.0 ) * ( p_dw_f0*exp(-p_dw_f1*(tx(3)-p_dw_f2)**2) + p_dw_f0*exp(-p_dw_f1*(tx(3)+p_dw_f2)**2) + p_dw_f3*exp(-p_dw_f4*tx(3)**2)))

    tx(k) = x(k)
  end do

  grad_lntfn(:) = (fhi(:)-flo(:))/(2.0_b8*ntol_eps)

  y = 1

end function grad_dw_tfun

function grad_lap_dw_tfun( x, sl, np, ndim, nslice, grad_lntfn, lap_lntfn ) result( y )
  real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
  real(kind=b8),intent(out) :: lap_lntfn 
  integer, intent(in) :: np, ndim, nslice, sl
  integer :: y

  real(kind=b8), dimension( ndim ) :: tgrad_lntfn 
  real(kind=b8), dimension( ndim ) :: tx
  real(kind=b8), dimension( ndim ) :: fhi,flo
  integer :: i, j, k

  lap_lntfn = 0
  grad_lntfn = 0
  do i = 1, np
    tx(:) = x(sl,i,:)
    do k = 1, ndim
      tx(k) = x(sl,i,k) + ntol_eps
      y = grad_dw_tfun( tx, ndim, tgrad_lntfn )
      fhi(k) = tgrad_lntfn(k)
      tx(k) = x(sl,i,k) - ntol_eps
      y  = grad_dw_tfun( tx, ndim, tgrad_lntfn )
      flo(k) = tgrad_lntfn(k)
      tx(k) = x(sl,i,k)
    end do
    y = grad_dw_tfun( tx, ndim, tgrad_lntfn )
    grad_lntfn(i,:) = tgrad_lntfn(:)
    lap_lntfn = lap_lntfn + sum((fhi(:)-flo(:))/(2.0_b8*ntol_eps))
  end do

end function grad_lap_dw_tfun

function nw_tfun( x, np, ndim ) result( y )
  integer :: np, ndim
  real(kind=b8), dimension( np , ndim ) :: x
  real(kind=b8) :: y
  real(kind=b8), dimension( np ) :: psi_t
  real(kind=b8) :: lam_nw = 1./10

  integer :: i, j

  psi_t(:)  = -( x(:,1)**2 + x(:,2)**2 + x(:,3)**2 )/2.0
  y = sum(psi_t)

end function nw_tfun

!> f := (f0*exp(-f1*(z-f2)**2) + f0*exp(-f1*(z+f2)**2))*exp(-f3*(x**2+y**2));
!                             2                       2             2    2
!    f := (f0 exp(-f1 (z - f2) ) + f0 exp(-f1 (z + f2) )) exp(-f3 (x  + y ))
!
!> simplify(laplacian(f,v));
!                               2  2          2  2          2  2
!2 f0 (-2 f3 %2 - 2 f3 %1 + 2 f3  x  %2 + 2 f3  x  %1 + 2 f3  y  %2
!
!           2  2                     2  2          2               2   2
!     + 2 f3  y  %1 - %2 f1 + 2 %2 f1  z  - 4 %2 f1  z f2 + 2 %2 f1  f2
!
!                      2  2          2               2   2
!     - %1 f1 + 2 %1 f1  z  + 4 %1 f1  z f2 + 2 %1 f1  f2 )
!
!               2       2       2                    2
!%1 := exp(-f3 x  - f3 y  - f1 z  - 2 f1 z f2 - f1 f2 )
!
!               2       2       2                    2
!%2 := exp(-f3 x  - f3 y  - f1 z  + 2 f1 z f2 - f1 f2 )

function dw_tfun_lap( x, np, ndim ) result( y )
  integer :: np, ndim
  real, dimension( np , ndim ) :: x
  real :: y

  y = 0


end function dw_tfun_lap

function grad_lap_dw_tfun2( x, sl, np, ndim, nslice, grad_lntfn, lap_lntfn, sp_param, spf_func ) result( y )
  real(kind=b8), dimension( nslice , np , ndim ), intent(in) :: x
  real(kind=b8), dimension( np , ndim ), intent(out) :: grad_lntfn 
  real(kind=b8), dimension( : ) :: sp_param
  real(kind=b8),intent(out) :: lap_lntfn 
  integer, intent(in) :: np, ndim, nslice, sl

  interface
    function spf_func( x, sl, param, nslice, np, ndim ) result( y )
      use kinds
      implicit none
      integer :: sl, nslice, np, ndim
      real(kind=b8), dimension( nslice, np , ndim ) :: x
      real(kind=b8), dimension( : ) :: param
      real(kind=b8) :: y
    end function spf_func
  end interface
  integer :: y

  y = 1
  grad_lntfn(:,1) = -p_hox*x(sl,:,1)
  grad_lntfn(:,2) = -p_hoy*x(sl,:,2)
  grad_lntfn(:,3) = -4.0*p_dw_f0*x(sl,:,3)*((x(sl,:,3)/p_dw_f1)**2-1)/p_dw_f1**2
  lap_lntfn = sum(-8.0*p_dw_f0*x(sl,:,3)**2/p_dw_f1**4-4.0*p_dw_f0*((x(sl,:,3)/p_dw_f1)**2-1)/p_dw_f1**2)
  lap_lntfn = lap_lntfn -np*(p_hox+p_hoy)

end function grad_lap_dw_tfun2

function dw_tfun2( x, sl, param, nslice, np, ndim ) result( y )
  use kinds
  implicit none
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( nslice, np , ndim ) :: x
  real(kind=b8), dimension(:) :: param
  real(kind=b8) :: y

  y =  sum(-(p_hox*x(sl,:,1)**2 + p_hoy*x(sl,:,2)**2)/2.0_b8 - p_dw_f0*((x(sl,:,3)/p_dw_f1)**2-1.0_b8)**2)
end function dw_tfun2

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

function lattice_tfun0( x, np, ndim ) result( y )
  integer :: np, ndim
  real(kind=b8), dimension( np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np ) :: psi_t

  integer :: i, j

  psi_t(:)  = -( p_hox*x(:,1)**2 + p_hoy*x(:,2)**2 + p_hoz*x(:,3)**2 )/2.0 + p_lwx*sin(x(:,1))**2 + p_lwz*sin(x(:,3))**2
  y = sum(psi_t)

end function lattice_tfun0

function lattice_tfun( x, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( nslice, np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np ) :: psi_t

  psi_t(:)  = -( p_hox*x(sl,:,1)**2 + p_hoy*x(sl,:,2)**2 + p_hoz*x(sl,:,3)**2 )/2.0 + &
              p_lat_w*sin(x(sl,:,1)*p_lattice_ax+p_lattice_phase_x)**2 + &
              p_lat_w*sin(x(sl,:,2)*p_lattice_ay+p_lattice_phase_y)**2 + &
              p_lat_w*sin(x(sl,:,3)*p_lattice_az+p_lattice_phase_z)**2
  y = sum(psi_t)

end function lattice_tfun

function fcc_tfun( x, sl, nslice, np, ndim ) result( y )
  integer :: sl, nslice, np, ndim
  real(kind=b8), dimension( nslice, np , ndim ) :: x
  real(kind=b8) :: y
  real, dimension( np, ndim ) :: psi_t

  psi_t(:,:)  = (x(sl,:,:)-p_lat_r0(:,:))**2
  y = -p_lat_w*sum(psi_t)/2.0_b8

end function fcc_tfun

end module vpi_trial_func
