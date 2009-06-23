module vpi_gfn

  use vpi_defines

  implicit none

contains 


function vpi_gfn2_sp( sl_start, sl_end, ip, U, nslice, np, ndim, dt ) result ( ln_gfn )
  real(kind=b8), dimension( nslice, np ), intent(in) :: U
  integer :: sl_start, sl_end, ip
  integer :: nslice, np, ndim
  real(kind=b8) :: dt
  real(kind=b8) :: ln_gfn

  integer :: ii,ioffset

  ln_gfn = -dt*sum( U(sl_start:sl_end,ip) )

!  print *,ln_gfn
end function vpi_gfn2_sp


function vpi_gfn4_sp( sl_start, sl_end, ip, U, gradU2, U_weight, gU2_weight, nslice, np, ndim, lambda, dt ) result ( ln_gfn )
  real(kind=b8), dimension( nslice, np ), intent(in) :: U
  real(kind=b8), dimension( nslice ), intent(in) :: gradU2
  real(kind=b8), dimension( nslice ), intent(in) :: U_weight
  real(kind=b8), dimension( nslice ), intent(in) :: gU2_weight
  integer :: sl_start, sl_end, ip
  integer :: nslice, np, ndim
  real(kind=b8) :: lambda
  real(kind=b8) :: dt
  real(kind=b8) :: ln_gfn

  integer :: ii,ioffset

  ln_gfn = -2.0_b8*dt*sum( U(sl_start:sl_end,ip)*U_weight(sl_start:sl_end) )/3.0_b8 &
           -2.0_b8*lambda*(dt**3)*sum( gradU2(sl_start:sl_end)*gU2_weight(sl_start:sl_end) )/(9.0_b8)  

!  print *,ln_gfn
end function vpi_gfn4_sp

! image approximation 
function vpi_hw_gfn( sl_start, sl_end, ip, q, nslice, np, ndim, dt ) result ( hw_gfn )
  integer :: sl_start, sl_end, ip
  integer :: nslice, np, ndim
  real(kind=b8) :: dt
  real(kind=b8), dimension ( nslice, np, ndim ) :: q
  real(kind=b8) :: hw_gfn,tgfn

  real(kind=b8) :: x,x0,xr,d
  integer :: i,j

  hw_gfn = 1.0_b8
  tgfn = 0.0_b8

  do i = 1, ndim 
    do j = sl_start+1, sl_end
      x0 = q(j-1,ip,i)
      x = q(j,ip,i)
      d = (abox - x0)
      xr = x0 + d*2.0_b8
      tgfn = exp(-( (xr-x)**2 )/(2*dt))
      d = (x0+abox)
      xr = x0 - d*2.0_b8
      tgfn = tgfn + exp(-( (xr-x)**2 )/(2*dt))
      hw_gfn = hw_gfn*(1.0_b8-tgfn*exp((x0-x)**2/(2*dt)))
    end do
  end do

end function vpi_hw_gfn

! image approximation 
function vpi_hs_gfn( sl_start, sl_end, ip, xij2, nslice, np, ndim, dt ) result ( hs_gfn )
  integer :: sl_start, sl_end, ip
  integer :: nslice, np, ndim
  real(kind=b8) :: dt
  real(kind=b8), dimension ( nslice, np , np ) :: xij2
  real(kind=b8) :: hs_gfn

  real(kind=b8) :: r,rp,d,dp
  integer :: i,j

  hs_gfn = 1.0_b8

  do i = 1, np
    if(i .ne. ip) then
      do j = sl_start+1, cslice,2
        r = sqrt(xij2(j-1,ip,i))
        rp = sqrt(xij2(j,ip,i))
        d = (r-a_hs)
        dp = (rp-a_hs)
        hs_gfn = hs_gfn*(1.0_b8 - exp(-d*dp/dt))
      end do
      do j = cslice+2,sl_end,2
        r = sqrt(xij2(j-1,ip,i))
        rp = sqrt(xij2(j,ip,i))
        d = (r-a_hs)
        dp = (rp-a_hs)
        hs_gfn = hs_gfn*(1.0_b8 - exp(-d*dp/dt))
      end do
    end if
  end do

end function vpi_hs_gfn

! A new quantum propagator for hard sphere and cavity systems
!J. Cao ad B.J. Berne
!JCP 97 (4) 2382  (1992)
function vpi_hs_gfn2( sl_start, sl_end, ip, q, nslice, np, ndim, dt ) result ( hs_gfn )
  integer :: sl_start, sl_end, ip
  integer :: nslice, np, ndim
  real(kind=b8) :: dt
  real(kind=b8), dimension ( nslice, np , ndim ) :: q
  real(kind=b8) :: hs_gfn

  real(kind=b8) :: tmp

  real(kind=b8), dimension ( ndim ) :: dx,dxp,dx12
  real(kind=b8) :: u,r,rp,drp,D,r12
  integer :: i,j

  hs_gfn = 1.0_b8

  do j = sl_start+1, cslice, 2
    do i = 1, np
      if(i .ne. ip) then
        dx(:) = q(j-1,ip,:) - q(j-1,i,:)
        r = sqrt(sum(dx(:)**2))
        dxp(:) = q(j,ip,:) - q(j,i,:)
        rp = sqrt(sum(dxp(:)**2))
        dx12(:) = dx(:) - dxp(:)
        r12 = sum(dx12(:)**2)
        u = r+rp
        tmp = a_hs*(u-a_hs)/(r*rp)
        D = u**2 - 2.0_b8*(a_hs*u - a_hs2 + tmp*dot_product(dx,dxp))
        hs_gfn = hs_gfn * (1.0_b8 - tmp*exp(-(D-r12)/(2.0_b8*lambda*dt)))
      end if
    end do
  end do

  do j = cslice+2, sl_end, 2
    do i = 1, np
      if(i .ne. ip) then
        dx(:) = q(j-1,ip,:) - q(j-1,i,:)
        r = sqrt(sum(dx(:)**2))
        dxp(:) = q(j,ip,:) - q(j,i,:)
        rp = sqrt(sum(dxp(:)**2))
        dx12(:) = dx(:) - dxp(:)
        r12 = sum(dx12(:)**2)
        u = r+rp
        tmp = a_hs*(u-a_hs)/(r*rp)
        D = u**2 - 2.0_b8*(a_hs*u - a_hs2 + tmp*dot_product(dx,dxp))
        hs_gfn = hs_gfn * (1.0_b8 - tmp*exp(-(D-r12)/(2.0_b8*lambda*dt)))
      end if
    end do
  end do

end function vpi_hs_gfn2

end module vpi_gfn
