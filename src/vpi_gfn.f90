!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1798:@thin vpi_gfn.f90
!@@language fortran90
module vpi_gfn

  use vpi_defines

  implicit none

contains 

!@+others
!@+node:gcross.20090624144408.1799:2nd order
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
!@nonl
!@-node:gcross.20090624144408.1799:2nd order
!@+node:gcross.20090624144408.1800:4th order
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
!@-node:gcross.20090624144408.1800:4th order
!@+node:dubois9.20090625101848.1680:hard wall
! image approximation
function vpi_hw_gfn( sl_start, sl_end, ip, q, nslice, np, ndim, dt ) result ( hw_gfn )
  integer :: sl_start, sl_end, ip
  integer :: nslice, np, ndim
  real(kind=b8) :: dt
  real(kind=b8), dimension ( nslice, np, ndim ) :: q
  real(kind=b8) :: hw_gfn,tgfn

  real(kind=b8) :: x,x0,xr,d
  integer :: i,j

  hw_gfn = 0.0_b8
  tgfn = 0.0_b8

  do i = 1, ndim
    do j = sl_start+1, sl_end
      x0 = q(j-1,ip,i)
      x = q(j,ip,i)
      d = (hard_wall_locations(i) - x)
      xr = hard_wall_locations(i) + d
      tgfn = exp(( -(xr-x0)**2  + (x-x0)**2)/(2.0_b8*lambda*dt))
      hw_gfn = hw_gfn + log(1.0_b8 - tgfn)

      d = (hard_wall_locations(i)+ x)
      xr = -hard_wall_locations(i) - d
      tgfn = exp(( -(xr-x0)**2  + (x-x0)**2)/(2.0_b8*lambda*dt))
      hw_gfn = hw_gfn+ log(1.0_b8 - tgfn)
    end do
  end do

  hw_gfn = exp(hw_gfn)

end function vpi_hw_gfn
!@-node:dubois9.20090625101848.1680:hard wall
!@+node:gcross.20090701101352.1736:(disabled)
!@+at
! ! image approximation
! function vpi_hs_gfn( sl_start, sl_end, ip, xij2, nslice, np, ndim, dt ) 
! result ( hs_gfn )
!   integer :: sl_start, sl_end, ip
!   integer :: nslice, np, ndim
!   real(kind=b8) :: dt
!   real(kind=b8), dimension ( nslice, np , np ) :: xij2
!   real(kind=b8) :: hs_gfn
! 
!   real(kind=b8) :: r,rp,d,dp
!   integer :: i,j
! 
!   hs_gfn = 1.0_b8
! 
!   do i = 1, np
!     if(i .ne. ip) then
!       do j = sl_start+1, cslice,2
!         r = sqrt(xij2(j-1,ip,i))
!         rp = sqrt(xij2(j,ip,i))
!         d = (r-hard_sphere_radius)
!         dp = (rp-hard_sphere_radius)
!         hs_gfn = hs_gfn*(1.0_b8 - exp(-d*dp/dt))
!       end do
!       do j = cslice+2,sl_end,2
!         r = sqrt(xij2(j-1,ip,i))
!         rp = sqrt(xij2(j,ip,i))
!         d = (r-hard_sphere_radius)
!         dp = (rp-hard_sphere_radius)
!         hs_gfn = hs_gfn*(1.0_b8 - exp(-d*dp/dt))
!       end do
!     end if
!   end do
! 
! end function vpi_hs_gfn
! 
! ! A new quantum propagator for hard sphere and cavity systems
! !J. Cao ad B.J. Berne
! !JCP 97 (4) 2382  (1992)
! function vpi_hs_gfn2( sl_start, sl_end, ip, q, nslice, np, ndim, dt ) result 
! ( hs_gfn )
!   integer :: sl_start, sl_end, ip
!   integer :: nslice, np, ndim
!   real(kind=b8) :: dt
!   real(kind=b8), dimension ( nslice, np , ndim ) :: q
!   real(kind=b8) :: hs_gfn
! 
!   real(kind=b8) :: tmp
! 
!   real(kind=b8), dimension ( ndim ) :: dx,dxp,dx12
!   real(kind=b8) :: u,r,rp,drp,D,r12
!   integer :: i,j
! 
!   hs_gfn = 1.0_b8
! 
!   do j = sl_start+1, cslice, 2
!     do i = 1, np
!       if(i .ne. ip) then
!         dx(:) = q(j-1,ip,:) - q(j-1,i,:)
!         r = sqrt(sum(dx(:)**2))
!         dxp(:) = q(j,ip,:) - q(j,i,:)
!         rp = sqrt(sum(dxp(:)**2))
!         dx12(:) = dx(:) - dxp(:)
!         r12 = sum(dx12(:)**2)
!         u = r+rp
!         tmp = hard_sphere_radius*(u-hard_sphere_radius)/(r*rp)
!         D = u**2 - 2.0_b8*(hard_sphere_radius*u - hard_sphere_radius_squared 
! + tmp*dot_product(dx,dxp))
!         hs_gfn = hs_gfn * (1.0_b8 - tmp*exp(-(D-r12)/(2.0_b8*lambda*dt)))
!       end if
!     end do
!   end do
! 
!   do j = cslice+2, sl_end, 2
!     do i = 1, np
!       if(i .ne. ip) then
!         dx(:) = q(j-1,ip,:) - q(j-1,i,:)
!         r = sqrt(sum(dx(:)**2))
!         dxp(:) = q(j,ip,:) - q(j,i,:)
!         rp = sqrt(sum(dxp(:)**2))
!         dx12(:) = dx(:) - dxp(:)
!         r12 = sum(dx12(:)**2)
!         u = r+rp
!         tmp = hard_sphere_radius*(u-hard_sphere_radius)/(r*rp)
!         D = u**2 - 2.0_b8*(hard_sphere_radius*u - hard_sphere_radius_squared 
! + tmp*dot_product(dx,dxp))
!         hs_gfn = hs_gfn * (1.0_b8 - tmp*exp(-(D-r12)/(2.0_b8*lambda*dt)))
!       end if
!     end do
!   end do
! 
! end function vpi_hs_gfn2
! 
!@-at
!@@c
!@-node:gcross.20090701101352.1736:(disabled)
!@-others

end module vpi_gfn
!@-node:gcross.20090624144408.1798:@thin vpi_gfn.f90
!@-leo
