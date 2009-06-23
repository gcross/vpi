module vpi_rotation
  
  use vpi_rand_utils
  use vpi_defines
  implicit none

  contains

function gfn_rot(q_rot,ip,sl_start,sl_end,dtau) result(gfn)
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM_ROT ), intent(in):: q_rot
  integer, intent(in) :: ip, sl_start, sl_end
  real(kind=b8), intent(in):: dtau
  real(kind=b8) :: gfn

  real(kind=b8) :: cx, gfn_i
  integer :: i


  gfn = 1.0_b8
!  write (2000,*) 
!  write (2000,*) sl_start," :: ",sl_end
  do i=sl_start,sl_end-1
     if((i .eq. CSLICE) .and. (ip .eq. od_pnum)) then
       cx = 1.0
     else
       cx = dot_product( q_rot(i,ip,:),q_rot(i+1,ip,:) )
       gfn_i = lookup_rot_gfn (cx,gfn_v)
       gfn = gfn*gfn_i
     end if
!     write (2000,*,ADVANCE="no") i,cx,gfn_i,gfn
  end do
!  write (3000,*) gfn

end function gfn_rot


function init_rot_gfn (nbins, lmax, dtau) result(gfn)
  integer, intent (IN) :: nbins, lmax
  real(kind=b8), dimension (nbins) :: gfn
  real(kind=b8), intent (IN) :: dtau

  real(kind=b8), dimension (0:lmax) :: P
  real(kind=b8), dimension (0:lmax) :: wght

  integer :: L,i
  real(kind=b8) :: x

  do L = 0, lmax
    wght(L) = (2.0_b8*L+1.0_b8)*exp(-dtau*p_B*L*(L+1))
  end do
  wght(:) = wght(:) / M_4PI

  P(0) = 1.0_b8
  do i = 1,nbins
    x = 2.0_b8*(real(i-1)/real(nbins-1)) - 1.0_b8
    P(1) = x
    do L = 1, lmax-1
      P(L+1) = ((2.0_b8*L+1.0_b8)*x*P(L) - L*P(L-1)) / (L+1)
    end do

    gfn(i) = abs(sum(P(:)*wght(:)))
  end do

end function init_rot_gfn


function lookup_rot_gfn (x,gfn_v) result(gfn)
  real(kind=b8), intent (IN) :: x
  real(kind=b8), intent (IN), dimension (:) :: gfn_v
  real(kind=b8) :: gfn

  integer :: nbins,ix
  real(kind=b8) :: rx,di

  nbins = size(gfn_v)

  rx = (x+1.0)*(nbins-1)/2.0 + 1
  ix = int(floor(rx))
  di = rx - ix
  
  if ( ix .ge. 1 ) then
    if (  ix .lt. nbins  ) then
      gfn = gfn_v(ix) + (gfn_v(ix+1)-gfn_v(ix))*di
    else
      gfn = gfn_v(nbins)
    end if
  else
    gfn = gfn_v(1)
  end if

end function lookup_rot_gfn

subroutine write_rot_gfn (gfn_v)
  real(kind=b8), intent (IN), dimension (:) :: gfn_v

  integer :: i,nsteps
  real(kind=b8) :: x
  
  nsteps = size(gfn_v)*10

  do i=1,size(gfn_v)
    x = 2.0_b8*real(i)/size(gfn_v)-1.0_b8
    write(5000,*) x, gfn_v(i)
  end do

  do i=0,nsteps
    x = 2.0_b8*real(i)/nsteps-1.0_b8
    write(4000,*) x, lookup_rot_gfn(x,gfn_v)
  end do

end subroutine write_rot_gfn

end module vpi_rotation

