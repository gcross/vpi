module vpi_bisect
  
  use vpi_rand_utils
  implicit none

  contains

recursive subroutine bisect( q , i0 , i1 , dtau ) 
  real(kind=b8), dimension( :, :, : ), intent(inout):: q
  integer, intent(in):: i0 , i1 
  real, intent(in):: dtau

  integer i,di
  integer :: j,k
  integer :: npart, ndim
  real(kind=b8), dimension( 2 )::  nu
  real  Pa

  Pa = 0

  npart = size(q,2)
  ndim = size(q,3)

  di = i1 - i0

  if ( di < 2 ) then
    return
  endif

  i = ( i0 + i1 ) / 2

  !print *,"bisect> [",i0,i,i1,"]"
  !print *,"bisect> shape(q) = ",shape(q)

  do j = 1, npart
    do k = 1, ndim
      call ru_gasdev( nu )
      q(i,j,k) = ( q(i0,j,k) + q(i1,j,k) ) / 2.0 + nu(1) * sqrt(di*dtau*2.0_b8)
      !print *, j, k, q(i0,j,k),  q(i1,j,k)
    end do
  end do

  call bisect( q, i0, i, dtau )
  call bisect( q, i, i1, dtau )

end subroutine bisect

recursive subroutine bisect_sp( q , i0 , i1 , dtau ) 
  real(kind=b8), dimension( :, : ), intent(inout):: q
  integer, intent(in):: i0 , i1 
  real, intent(in):: dtau

  integer :: i,di
  integer :: j
  integer :: s0,s1
  integer :: ndim, nslice
  real(kind=b8), dimension( 2 )::  nu
  real :: Pa

  Pa = 0

  nslice = size(q,1)
  ndim = size(q,2)
  
  !print *, nslice, ndim

  di = i1 - i0

  if ( di < 2 ) then
    return
  endif

  i = ( i0 + i1 ) / 2

  !print *,"bisect_sp> [",i0,i,i1,"]"
  !print *,"bisect_sp> shape(q) = ",shape(q)

  do j = 1, ndim
    call ru_gasdev( nu )
    q(i,j) = ( q(i0,j) + q(i1,j) ) / 2.0 + nu(1) * sqrt(di*dtau*2.0_b8)
!    print *, j, q(i0,j),  q(i1,j)
!    print *,"bisect_sp> ",j,": [",i0,i,i1,"]"
!    print *,"bisect_sp> ",j,": [",q(i0,j),q(i,j),q(i1,j),"]"
  end do

  call bisect_sp( q, i0, i, dtau )
  call bisect_sp( q, i, i1, dtau )

end subroutine bisect_sp

subroutine bisect_endpoint( q , idx_end , idx_start , dtau ) 
  real(kind=b8), dimension( :, : ), intent(inout):: q
  integer, intent(in):: idx_end , idx_start
  real, intent(in):: dtau

  integer :: i,di
  integer :: j
  integer :: s1
  integer :: ndim, nslice
  real(kind=b8), dimension( 2 )::  nu
  real :: Pa

  Pa = 0

  nslice = size(q,1)
  ndim = size(q,2)
  
  !print *, nslice, ndim

  di = abs(idx_end - idx_start)

  if ( di .lt. 1 ) then
    return
  endif

!  print *,"bisect_endpoint> [",idx_end,idx_start,"]"
  !print *,"bisect_sp> shape(q) = ",shape(q)

  do j = 1, ndim
    call ru_gasdev( nu )
    q(idx_end,j) = q(idx_start,j)  + nu(1) * sqrt(di*dtau*2.0_b8)
!    print *, j, q(i0,j),  q(i1,j)
!    print *,"bisect_sp> ",j,": [",i0,i,i1,"]"
!    print *,"bisect_sp> ",j,": [",q(i0,j),q(i,j),q(i1,j),"]"
  end do

  if( idx_end .gt. idx_start ) then 
    call bisect_sp( q, idx_start, idx_end, dtau )
  else
    call bisect_sp( q, idx_end, idx_start, dtau )
  end if
end subroutine bisect_endpoint

recursive subroutine bisect_sp_1d( q , id, i0 , i1 , dtau ) 
  real(kind=b8), dimension( :, : ), intent(inout):: q
  integer, intent(in):: id, i0 , i1 
  real, intent(in):: dtau

  integer :: i, j, di
  real(kind=b8), dimension( 2 )::  nu

  di = i1 - i0

  if ( di < 2 ) then
    return
  endif

  i = ( i0 + i1 ) / 2


  call ru_gasdev( nu )
  q(i,id) = ( q(i0,id) + q(i1,id) ) / 2.0 + nu(1) * sqrt(di*dtau*2.0_b8)

  call bisect_sp_1d( q, id, i0, i, dtau )
  call bisect_sp_1d( q, id, i, i1, dtau )

end subroutine bisect_sp_1d

subroutine bisect_endpoint_1d( q , id, idx_end , idx_start , dtau ) 
  real(kind=b8), dimension( :, : ), intent(inout):: q
  integer, intent(in):: id, idx_end , idx_start
  real, intent(in):: dtau

  integer :: i, j, di
  integer :: s1
  real(kind=b8), dimension( 2 )::  nu

  di = abs(idx_end - idx_start)

  if ( di .lt. 1 ) then
    return
  endif

  call ru_gasdev( nu )
  q(idx_end,id) = q(idx_start,id)  + nu(1) * sqrt(di*dtau*2.0_b8)

  if( idx_end .gt. idx_start ) then 
    call bisect_sp_1d( q, id, idx_start, idx_end, dtau )
  else
    call bisect_sp_1d( q, id, idx_end, idx_start, dtau )
  end if
end subroutine bisect_endpoint_1d

end module vpi_bisect

