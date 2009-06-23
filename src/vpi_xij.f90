module vpi_xij

use vpi_defines
implicit none

contains 

subroutine vpi_update_xij( xij2, x, sl_start, sl_end, ip, nslice, np, ndim )
  integer, intent(in) :: nslice, np, ndim
  integer, intent(in) :: sl_start, sl_end, ip
  real(kind=b8), dimension ( nslice, np , np ), intent(inout) :: xij2
  real(kind=b8), dimension ( nslice, np , ndim ), intent(in) :: x

  integer j

!  print *,x(slice,ip,:)

  do j = 1, ip-1
    xij2(sl_start:sl_end,ip,j) =  ( x(sl_start:sl_end,ip,1) - x(sl_start:sl_end,j,1) )**2 +&
                                  ( x(sl_start:sl_end,ip,2) - x(sl_start:sl_end,j,2) )**2 +&
                                  ( x(sl_start:sl_end,ip,3) - x(sl_start:sl_end,j,3) )**2
    xij2(sl_start:sl_end,j,ip) =  xij2(sl_start:sl_end,ip,j)
  end do

  do j = ip+1, np
    xij2(sl_start:sl_end,ip,j) =  ( x(sl_start:sl_end,ip,1) - x(sl_start:sl_end,j,1) )**2 +&
                                  ( x(sl_start:sl_end,ip,2) - x(sl_start:sl_end,j,2) )**2 +&
                                  ( x(sl_start:sl_end,ip,3) - x(sl_start:sl_end,j,3) )**2
    xij2(sl_start:sl_end,j,ip) =  xij2(sl_start:sl_end,ip,j)
  end do
end subroutine vpi_update_xij

subroutine vpi_update_xij_pbc( xij2, x, sl_start, sl_end, ip, nslice, np, ndim )
  integer, intent(in) :: nslice, np, ndim
  integer, intent(in) :: sl_start, sl_end, ip
  real(kind=b8), dimension ( nslice, np , np ), intent(inout) :: xij2
  real(kind=b8), dimension ( nslice, np , ndim ), intent(in) :: x

  real(kind=b8), dimension ( ndim ) :: dxij,pbc_dxij
  integer j,k

!  print *,x(slice,ip,:)

  do j = 1, ip-1
    do k = sl_start, sl_end
      dxij(:) =  x(k,ip,:) - x(k,j,:)
      pbc_dxij(:) = dxij(:) - p_pbc_L*floor(dxij(:)/p_pbc_L - 0.5_b8) - p_pbc_L
      xij2(k,ip,j) = sum(pbc_dxij(:)**2)
      xij2(k,j,ip) =  xij2(k,ip,j)
!      print *,xij2(k,j,ip)
    end do
  end do

  do j = ip+1, np
    do k = sl_start, sl_end
      dxij(:) =  x(k,ip,:) - x(k,j,:)
      pbc_dxij(:) = dxij(:) - p_pbc_L*floor(dxij(:)/p_pbc_L - 0.5_b8) - p_pbc_L
      xij2(k,ip,j) = sum(pbc_dxij(:)**2)
      xij2(k,j,ip) =  xij2(k,ip,j)
!      print *,xij2(k,j,ip)
    end do
  end do

end subroutine vpi_update_xij_pbc

end module vpi_xij
