module vpi_obdm
  
  use vpi_defines
  implicit none

contains 

subroutine vpi_eval_obdm_cut( obdm, x1, x2, nbins, ndim, xsize, dndx)
  integer :: id, nbins, ndim
  real(kind=b8), dimension ( nbins, nbins ) :: obdm
  real(kind=b8) :: x1, x2
  real :: xsize, dndx

  integer :: bx1, bx2

!  print *, xsize, dndx
!  print *, x1, x2
  bx1 = floor((x1 + xsize)*dndx)+1
  bx2 = floor((x2 + xsize)*dndx)+1
!  print *, bx1,bx2
  if ( (bx1 .le. nbins) .and. (bx2 .le. nbins) ) then
    if ( (bx1 .ge. 1) .and. (bx2 .ge. 1) ) then
       obdm(bx1,bx2) = obdm(bx1,bx2) + 1
    end if
  end if

end subroutine vpi_eval_obdm_cut

subroutine vpi_eval_obdm_ring( obdm, x1, x2, nbins, ndim)
  integer :: id, nbins, ndim
  real(kind=b8), dimension ( nbins, nbins ) :: obdm
  real(kind=b8), dimension( 3 ) :: x1, x2
  real :: xsize, dndx
  real :: t1,t2

  integer :: bx1, bx2

!  print *, xsize, dndx
!  print *, x1, x2
  if( x1(1) .gt. 0 ) then
    t1 = (atan(x1(3)/x1(1))/M_2PI) + 0.25
  else
    t1 = (atan(x1(3)/x1(1))/M_2PI) + 0.75
  end if

  if( x2(1) .gt. 0 ) then
    t2 = (atan(x2(3)/x2(1))/M_2PI) + 0.25
  else
    t2 = (atan(x2(3)/x2(1))/M_2PI) + 0.75
  end if

  bx1 = floor(t1*nbins)+1
  bx2 = floor(t2*nbins)+1
!  print *, bx1,bx2
  if ( (bx1 .le. nbins) .and. (bx2 .le. nbins) ) then
    if ( (bx1 .ge. 1) .and. (bx2 .ge. 1) ) then
       obdm(bx1,bx2) = obdm(bx1,bx2) + 1
    end if
  end if

end subroutine vpi_eval_obdm_ring

subroutine vpi_eval_obdm_full( obdmf, x1, x2, nbins, ndim, xsize, dndx)
  integer :: nbins, ndim
  real(kind=b8), dimension ( nbins**ndim, nbins**ndim ) :: obdmf
  real(kind=b8), dimension( ndim ) :: x1, x2
  real,dimension(ndim) :: xsize, dndx

  integer, dimension(ndim)  :: bx1, bx2
  integer :: ii,i1,i2

  bx1(:) = floor((x1(:) + xsize(:))*dndx)+1
  bx2(:) = floor((x2(:) + xsize(:))*dndx)+1
  i1 = 0
  i2 = 0
  do ii=1,ndim
    i1 = i1 + bx1(ii)*nbins**(ii-1)
    i2 = i2 + bx2(ii)*nbins**(ii-1)
  enddo
  print *, bx1(:)
  print *, bx2(:)
  print *, i1, i2
  if ( (i1 .ge. 1) .and. (i2 .ge. 1) ) then
    if ( (i1 .le. nbins**ndim) .and. (i2 .le. nbins**ndim) ) then
       obdmf(i1,i2) = obdmf(i1,i2) + 1
    end if
  end if

end subroutine vpi_eval_obdm_full

subroutine vpi_eval_nrdm(nrdm, x1, x2)
  integer :: id, nbins, ndim
  real(kind=b8), dimension ( 2**N_OD_PARTICLE, 2**N_OD_PARTICLE ) :: nrdm
  real(kind=b8), dimension ( N_PARTICLE ) :: x1, x2
  real :: xsize, dndx

  integer :: i,j,bx1, bx2

  bx1 = 1
  bx2 = 1
  do i = 0, N_OD_PARTICLE-1
    if(x1(i+1) > 0) then
      bx1 = bx1 + 2**i
    end if
    if(x2(i+1) > 0) then
      bx2 = bx2 + 2**i
    end if
  end do

  nrdm(bx1,bx2) = nrdm(bx1,bx2) + 1

end subroutine vpi_eval_nrdm


end module vpi_obdm
