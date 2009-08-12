!@+leo-ver=4-thin
!@+node:gcross.20090812093015.1722:@thin gfn.f95
!@@language fortran90
!@@tabwidth -2

module gfn

  implicit none

contains 

!@+others
!@+node:gcross.20090812093015.1753:initialize weights
subroutine initialize_4th_order_weights(n_slices,U_weight,gU2_weight)
  integer, intent(in) :: n_slices
  double precision, dimension(n_slices), intent(out) :: U_weight, gU2_weight
  integer :: cslice, ii

  cslice = n_slices / 2

  if (mod(CSLICE,2) .ne. 1) then
    print *,"ERROR: CSLICE = N_SLICE/2 must be odd"
    stop
  end if

!@+at
! These "slice weights" are used to give the proper weight to even
! and odd terms when calculating the Green's function.
! In the fourth order short time approximation for the Green's funciton,
! odd slices contribute twice the potential energy of even slices
! the gradient of the potential energy only contributes for odd slides.
!@-at
!@@c
  U_weight = 0.0
  gU2_weight = 0.0
!@+at
! The fourth order propagator works in sets of 3 like this
!   0.5 1 0.5 , 0.5 1 0.5
! We concatentate the adjacent half weighted steps (unless we are at the end 
! of a path)
! in order for paths to be broken properly at the center we need CSLICE to be 
! odd.
!@-at
!@@c
  do ii = 1, CSLICE
    U_weight(ii) = mod(ii+1,2) + 1 
    gU2_weight(ii) = mod(ii+1,2)
  end do
  ! the center slice is doubled to handle broken paths 
  ! (probably a less kludgy way to do this but... )
  do ii = CSLICE+2,n_slices
    U_weight(ii) = mod(ii,2) + 1 
    gU2_weight(ii) = mod(ii,2)
  end do
  U_weight(1) = 0.5d0
  U_weight(N_SLICES) = 0.5d0
  U_weight(CSLICE) = 0.5d0
  U_weight(CSLICE+1) = 0.5d0
end subroutine
!@-node:gcross.20090812093015.1753:initialize weights
!@+node:gcross.20090624144408.1799:2nd order
pure function gfn2_sp( sl_start, sl_end, ip, U, nslice, np, dt ) result ( ln_gfn )
  double precision, dimension( nslice, np ), intent(in) :: U
  integer, intent(in) :: sl_start, sl_end, ip
  integer, intent(in) :: nslice, np
  double precision, intent(in) :: dt
  double precision :: ln_gfn

  ln_gfn = -dt*sum( U(sl_start:sl_end,ip) )

end function gfn2_sp
!@-node:gcross.20090624144408.1799:2nd order
!@+node:gcross.20090812093015.1845:4th order
function gfn4_sp( sl_start, sl_end, ip, U, gradU2, U_weight, gU2_weight, nslice, np, lambda, dt ) result ( ln_gfn )
  double precision, dimension( nslice, np ), intent(in) :: U
  double precision, dimension( nslice ), intent(in) :: gradU2
  double precision, dimension( nslice ), intent(in) :: U_weight
  double precision, dimension( nslice ), intent(in) :: gU2_weight
  integer, intent(in) :: sl_start, sl_end, ip
  integer, intent(in) :: nslice, np
  double precision, intent(in) :: lambda, dt
  double precision :: ln_gfn

  integer :: slice_length
  double precision, external :: ddot

  slice_length = sl_end-sl_start+1

  ln_gfn = -2.0d0*dt*ddot(slice_length,U(sl_start,ip),1,U_weight(sl_start),1)/3.0d0 &
           -2.0d0*lambda*(dt**3)*ddot(slice_length,gradU2(sl_start),1,gU2_weight(sl_start),1)/9.0d0

end function gfn4_sp
!@-node:gcross.20090812093015.1845:4th order
!@+node:dubois9.20090625101848.1680:hard wall
!@+at
! ! image approximation
! function vpi_hw_gfn( sl_start, sl_end, ip, q, nslice, np, ndim, dt ) result 
! ( hw_gfn )
!   integer :: sl_start, sl_end, ip
!   integer :: nslice, np, ndim
!   real(kind=b8) :: dt
!   real(kind=b8), dimension ( nslice, np, ndim ) :: q
!   real(kind=b8) :: hw_gfn,tgfn
! 
!   real(kind=b8) :: x,x0,xr,d
!   integer :: i,j
! 
!   hw_gfn = 0.0_b8
!   tgfn = 0.0_b8
! 
!   do i = 1, ndim
!     do j = sl_start+1, sl_end
!       x0 = q(j-1,ip,i)
!       x = q(j,ip,i)
!       d = (hard_wall_locations(i) - x)
!       xr = hard_wall_locations(i) + d
!       tgfn = exp(( -(xr-x0)**2  + (x-x0)**2)/(2.0_b8*lambda*dt))
!       hw_gfn = hw_gfn + log(1.0_b8 - tgfn)
! 
!       d = (hard_wall_locations(i)+ x)
!       xr = -hard_wall_locations(i) - d
!       tgfn = exp(( -(xr-x0)**2  + (x-x0)**2)/(2.0_b8*lambda*dt))
!       hw_gfn = hw_gfn+ log(1.0_b8 - tgfn)
!     end do
!   end do
! 
!   hw_gfn = exp(hw_gfn)
! 
! end function vpi_hw_gfn
!@-at
!@@c

!@-node:dubois9.20090625101848.1680:hard wall
!@+node:gcross.20090624144408.1800:(broken 4th order)
!@+at
! function gfn4_sp( sl_start, sl_end, ip, U, gradU2, nslice, np, lambda, dt ) 
! result ( ln_gfn )
!   double precision, dimension( nslice, np ), intent(in) :: U
!   double precision, dimension( nslice ), intent(in) :: gradU2
!   integer, intent(in) :: sl_start, sl_end, ip
!   integer, intent(in) :: nslice, np
!   double precision, intent(in) :: lambda
!   double precision, intent(in) :: dt
! 
!   integer :: k, cslice, start_boundary, end_boundary
!   double precision :: coefficient
!   double precision :: ln_gfn
!   double precision :: ln_gfn_U, ln_gfn_gU2
! 
!   cslice = nslice/2
! 
!   ln_gfn_U = 0d0
!   ln_gfn_gU2 = 0d0
! 
!   if (sl_start <= 1) then
!     ln_gfn_U   = ln_gfn_U   + U(1,ip)  *0.5d0
!   end if
! 
!   start_boundary = max(2,sl_start)
!   end_boundary   = min(cslice-1,sl_end)
!   if (start_boundary <= end_boundary .and. end_boundary >= 2 .and. 
! start_boundary <= cslice-1) then
!     coefficient = dble(2d0-mod(start_boundary,2))
!     do k = start_boundary, end_boundary
!       ln_gfn_U   = ln_gfn_U   + U(k,ip)  *coefficient
!       coefficient = 3.0d0 - coefficient
!     end do
!     ln_gfn_gU2 = ln_gfn_gU2 + sum(gradU2(start_boundary:end_boundary:2))
!   end if
! 
!   if (sl_start < cslice .and. cslice < sl_end) then
!     ln_gfn_U   = ln_gfn_U   + U(cslice,ip)  *0.5d0
!   end if
! 
!   if (sl_start < cslice+1 .and. cslice+1 < sl_end) then
!     ln_gfn_U   = ln_gfn_U   + U(cslice+1,ip)*0.5d0
!   end if
! 
! 
!   start_boundary = max(cslice+2,sl_start)
!   end_boundary   = min(nslice-1,sl_end)
!   if (start_boundary <= end_boundary .and. end_boundary >= cslice+2 .and. 
! start_boundary <= nslice-1) then
!     coefficient = dble(mod(start_boundary,2)+1d0)
!     do k = start_boundary, end_boundary
!       ln_gfn_U   = ln_gfn_U   + U(k,ip)  *coefficient
!       coefficient = 3.0d0 - coefficient
!     end do
!     ln_gfn_gU2 = ln_gfn_gU2 + sum(gradU2(start_boundary:end_boundary:2))
!   end if
! 
!   if (sl_end == nslice ) then
!     ln_gfn_U   = ln_gfn_U   + U(nslice,ip)  *0.5d0
!   end if
! 
!   print *, ln_gfn_U, ln_gfn_gU2
! 
!   ln_gfn = -2.0d0*dt*ln_gfn_U/3.0d0 &
!            -2.0d0*lambda*(dt**3)*ln_gfn_gU2/(9.0d0)
! 
! end function gfn4_sp
!@-at
!@@c
!@-node:gcross.20090624144408.1800:(broken 4th order)
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

end module gfn
!@-node:gcross.20090812093015.1722:@thin gfn.f95
!@-leo
