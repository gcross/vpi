module vpi_bbridge
  
  use vpi_rand_utils
  use vpi_defines
  implicit none

  contains

! given the free particle hamiltonian
! $$
! H  = -\lambda \nabla^2
! $$
! the exact propagator is
! $$
! (4 \pi \lambda \tau)^(-3N /2) \exp( (r - r')^2/(4 \lambda \tau) )
! $$

!/*
! * Form a Brownian Bridge from q(i0) to q(i1) with time step dtau
! */
subroutine bbridge(q0, q1, pnum, low_dim, high_dim, btype, i0, i1, lambda, dtau)
  real(kind=b8), dimension( :, :, : ), intent(in):: q0
  real(kind=b8), dimension( :, :, : ), intent(inout):: q1
  integer, intent(in):: pnum, low_dim, high_dim, btype
  integer, intent(inout) :: i0
  integer, intent(inout) :: i1
  real(kind=b8), intent(in):: lambda
  real(kind=b8), intent(in):: dtau

  real(kind=b8) :: dt

  integer :: i
  integer :: j,cnt
  integer :: s0,s1
  integer :: ndim
  real(kind=b8) :: t1,t2,mean,sigma,sigmasq
  real(kind=b8) :: a,b

  real(kind=b8), dimension( 2*ceiling((i1-i0+2.)*(high_dim-low_dim+1.)/2.) )::  nu

  dt = lambda*dtau*2.0_b8

  call ru_gasdev( nu )

  !print *, 'i0= ',i0,'i1 = ',i1,' btype= ',btype 
  !write (111,*) 'i0= ',i0,'i1 = ',i1,' btype= ',btype 
  q1(i0, pnum, low_dim:high_dim) = q0(i0, pnum, low_dim:high_dim)
  q1(i1, pnum, low_dim:high_dim) = q0(i1, pnum, low_dim:high_dim)

  cnt = 1

  if(btype .eq. LEFT_BRIDGE) then
    sigma= sqrt(dt*(i1-1))
    do j=low_dim, high_dim
      mean = q1(i1, pnum, j)
      q1(i0, pnum, j) = nu(cnt)*sigma + mean
      cnt = cnt + 1
    end do
  else if(btype .eq. RIGHT_BRIDGE) then
    sigma= sqrt(dt*(N_SLICE-i0))
    do j=low_dim, high_dim
      mean = q1(i0, pnum, j)
      q1(i1, pnum, j) = nu(cnt)*sigma + mean
      cnt = cnt + 1
    end do
  end if

  if( (i0 .le. CSLICE) .and. (i1 .gt. CSLICE) ) then
    i1 = i1 + 1
    if(i1 .gt. N_SLICE) then
      i1 = N_SLICE
    end if

    if( ((pnum .eq. OD_PNUM) .and. eval_off_diagonal) .or. ((pnum .le. N_OD_PARTICLE) .and. (eval_nrdm)) ) then
      sigma= sqrt(dt*(CSLICE-i0))
      do j = OD_DIM_LOW,OD_DIM_HIGH
        mean = q1(i0, pnum, j)
        q1(CSLICE, pnum, j) = nu(cnt)*sigma + mean
        cnt = cnt + 1
      end do
      sigma= sqrt(dt*(i1-(CSLICE+1)))
      do j = OD_DIM_LOW,OD_DIM_HIGH
        mean = q1(i1, pnum, j)
        q1(CSLICE+1, pnum, j) = nu(cnt)*sigma + mean
        cnt = cnt + 1
      end do
    else
      t1 = dt*(CSLICE-i0)
      t2 = dt*(i1-(CSLICE+1))
      sigmasq=1./( 1./t1 + 1./t2)
      sigma= sqrt(sigmasq)
      do j = low_dim, high_dim
        a = q1(i0, pnum, j)
        b = q1(i1, pnum, j)
        mean=(t2*a  + t1*b)/(t1 + t2)
        q1(CSLICE, pnum, j) = nu(cnt)*sigma + mean
        q1(CSLICE+1, pnum, j) = q1(CSLICE, pnum, j)
        cnt = cnt + 1
      end do
    end if

    do i=i0+1, CSLICE-1
      t1 = dt
      t2 = dt*(CSLICE-i)
      sigmasq=1./( 1./t1 + 1./t2)
      sigma= sqrt(sigmasq)
      do j=low_dim, high_dim
        a = q1(i-1, pnum, j)
        b = q1(CSLICE, pnum, j)
        mean=(t2*a  + t1*b)/(t1 + t2)
        q1(i, pnum, j) = nu(cnt)*sigma + mean
        cnt = cnt + 1
      end do
    end do

    do i=CSLICE+2, i1-1
      t1 = dt
      t2 = dt*(i1-i)
      sigmasq=1./( 1./t1 + 1./t2)
      sigma= sqrt(sigmasq)
      do j=low_dim, high_dim
        a = q1(i-1, pnum, j)
        b = q1(i1, pnum, j)
        mean=(t2*a  + t1*b)/(t1 + t2)
        q1(i, pnum, j) = nu(cnt)*sigma + mean
        cnt = cnt + 1
      end do
    end do
  else
    do i=i0+1, i1-1
      t1 = dt
      t2 = dt*(i1-i)
      sigmasq=1./( 1./t1 + 1./t2)
      sigma= sqrt(sigmasq)
      do j=low_dim, high_dim
        a = q1(i-1, pnum, j)
        b = q1(i1, pnum, j)
        mean=(t2*a  + t1*b)/(t1 + t2)
        q1(i, pnum, j) = nu(cnt)*sigma + mean
        cnt = cnt + 1
      end do
    end do
  end if
end subroutine bbridge

end module vpi_bbridge

