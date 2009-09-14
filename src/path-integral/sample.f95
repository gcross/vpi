!@+leo-ver=4-thin
!@+node:gcross.20090623152316.73:@thin sample.f95
!@@language fortran90
module sample

 !@ << Imported modules >>
 !@+node:gcross.20090623152316.74:<< Imported modules >>
 use constants
 use rand_utils
 !@-node:gcross.20090623152316.74:<< Imported modules >>
 !@nl

 implicit none

 !@ << Move type constants >>
 !@+node:gcross.20090623152316.75:<< Move type constants >>
 integer, parameter :: MT_BBRIDGE = 1
 integer, parameter :: MT_RIGID = 2
 integer, parameter :: MT_SWAP = 3
 integer, parameter :: N_MOVE_TYPES = 3

 integer, parameter :: LEFT_BRIDGE = 1
 integer, parameter :: RIGHT_BRIDGE = 2
 integer, parameter :: S_BRIDGE = 3
 !@-node:gcross.20090623152316.75:<< Move type constants >>
 !@nl

 contains

!@+others
!@+node:gcross.20090623152316.76:Sampling subroutines
!@+node:gcross.20090623152316.77:scheme 1
subroutine sample_scheme1( &
    q0, q1, &
    move_start, move_end, &
    move_type_probabilities, move_type_differentials, &
    dM, lambda, &
    low_swap_dim, high_swap_dim, &
    part_num, mtype, &
    N_SLICE, N_PARTICLE, N_DIM, &
    od_pnum, PROB_OD_PNUM, &
    n_od_particle, &
    od_dim_low, od_dim_high &
  )
  integer, intent(in) :: N_SLICE, N_PARTICLE, N_DIM
  integer, parameter :: N_MOVE_TYPES = 3
  double precision, dimension ( N_MOVE_TYPES ), intent(in) :: move_type_probabilities, move_type_differentials
  double precision, dimension ( N_SLICE , N_PARTICLE, N_DIM ), intent(in) :: q0
  double precision, dimension ( N_SLICE , N_PARTICLE, N_DIM ), intent(inout) :: q1
  integer, intent(in) :: low_swap_dim, high_swap_dim
  integer, intent(out) :: move_start,move_end,part_num,mtype
  double precision, intent(in), optional :: PROB_OD_PNUM
  integer, intent(in), optional :: od_pnum, n_od_particle, od_dim_low, od_dim_high
  double precision, intent(in) :: dM, lambda

  double precision :: rnd
  integer :: i0,i1,cslice
  integer :: swap_dim
  cslice = n_slice/2

  call random_number( rnd )
  if (present(od_pnum) .and. ( rnd .lt. PROB_od_pnum ) ) then
    part_num = od_pnum
  else
    part_num = ceiling(rnd*N_PARTICLE)
  end if

  call random_number( rnd )
  if (rnd < move_type_probabilities(MT_BBRIDGE)) then
    mtype = MT_BBRIDGE
    call random_number( rnd )
    i0 = floor(rnd*((N_SLICE-1)-(2-dM))) - dM + 2
    i1 = i0 + dM
! print *, "i0, i1 ",i0,i1
    if( ((i0 .ge. 1) .and. (i1 .le. N_SLICE)) ) then
      call bbridge(q0, q1, part_num, 1, N_DIM, S_BRIDGE, i0, i1, lambda, &
            move_type_differentials(MT_BBRIDGE), N_SLICE, N_PARTICLE, N_DIM)
      move_start = i0+1
      move_end = i1-1
    else
      if(i0 .lt. 1) then
        i0 = 1
        call bbridge(q0, q1, part_num, 1, N_DIM, LEFT_BRIDGE, i0, i1, lambda, &
            move_type_differentials(MT_BBRIDGE), N_SLICE, N_PARTICLE, N_DIM)
        move_start = 1
        move_end = i1-1
      else
        if(i1 .gt. N_SLICE) then
          i1 = N_SLICE
          call bbridge(q0, q1, part_num, 1, N_DIM, RIGHT_BRIDGE, i0, i1, lambda, &
            move_type_differentials(MT_BBRIDGE), N_SLICE, N_PARTICLE, N_DIM, &
            od_pnum, n_od_particle, &
            od_dim_low, od_dim_high &
            )
          move_start = i0+1
          move_end = N_SLICE
        end if
      end if
    end if
  else 
! move all time steps of a single particle at once
    if (rnd < move_type_probabilities(MT_BBRIDGE) + move_type_probabilities(MT_RIGID)) then
      mtype = MT_RIGID
      call rigid_move(q0, q1, part_num, 1, N_DIM, 1, N_SLICE, &
              move_type_differentials(MT_RIGID), N_SLICE, N_PARTICLE, N_DIM)
      if( (present(od_pnum) .and. (part_num .eq. od_pnum)) .or. (present(n_od_particle) .and. (part_num .le. N_OD_PARTICLE)) ) then
        call rigid_move(q0, q1, part_num, OD_DIM_LOW, OD_DIM_HIGH, 1, CSLICE, &
                move_type_differentials(MT_RIGID), N_SLICE, N_PARTICLE, N_DIM)
        call rigid_move(q0, q1, part_num, OD_DIM_LOW, OD_DIM_HIGH, CSLICE+1, N_SLICE, &
                move_type_differentials(MT_RIGID), N_SLICE, N_PARTICLE, N_DIM)
      end if
      move_start = 1
      move_end = N_SLICE
    else 
! attempt a rotation around a symmetry axis.  right now only does moves like q1 = -q0
    if (rnd < move_type_probabilities(MT_BBRIDGE) + move_type_probabilities(MT_RIGID) + move_type_probabilities(MT_SWAP) ) then
        mtype = MT_SWAP
        call random_number( rnd )
        swap_dim = int(floor( rnd * (high_swap_dim-low_swap_dim+1) ))+low_swap_dim
        if( (present(od_pnum) .and. (part_num == OD_PNUM)) .or. (present(N_OD_PARTICLE) .and. (part_num <= N_OD_PARTICLE)) ) then
          call random_number( rnd )
          if(rnd .ge. 0.5) then
            call swap_move(q0, q1, part_num, swap_dim, swap_dim, 1, CSLICE, &
                    move_type_differentials(MT_SWAP), N_SLICE, N_PARTICLE, N_DIM)
          else
            call swap_move(q0, q1, part_num, swap_dim, swap_dim, CSLICE+1, N_SLICE, &
                    move_type_differentials(MT_SWAP), N_SLICE, N_PARTICLE, N_DIM)
          end if
        else
          call cascading_swap_move(q0, q1, part_num, swap_dim, N_SLICE, N_PARTICLE, N_DIM)
        end if
        move_start = 1
        move_end = N_SLICE
      end if
    end if
  end if

!  print *, "move_start, move_end ",move_start, move_end,part_num
end subroutine sample_scheme1
!@-node:gcross.20090623152316.77:scheme 1
!@-node:gcross.20090623152316.76:Sampling subroutines
!@+node:gcross.20090623152316.84:Utility subroutines
!@+node:gcross.20090623152316.85:ipush / ipop
subroutine ipush(list,val,ip)
  integer, dimension(:),intent(inout) :: list
  integer, intent(in) :: val
  integer,intent(inout) :: ip

  if(ip < size(list)) then
    ip = ip + 1
    list(ip) = val 
  end if
end subroutine ipush!}}}

subroutine ipop(list,val,ip)
  integer, dimension(:),intent(inout) :: list
  integer, intent(inout) :: val,ip

  if(ip > 0) then
    val = list(ip)
    ip = ip -1 
  else
    val = -1
    ip = 0
  end if
end subroutine ipop!}}}
!@-node:gcross.20090623152316.85:ipush / ipop
!@+node:gcross.20090623152316.88:get_collisions
!@+at
! subroutine get_collisions(q,ip,clist,n_coll,N_SLICE,N_PARTICLE,N_DIM)
!   integer, intent(in) :: N_SLICE, N_PARTICLE, N_DIM
!   double precision, dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: q
!   integer, dimension ( N_PARTICLE ) :: clist
!   integer :: ip,n_coll
!   integer :: i,j,cflag
!   double precision :: xij
! 
!   n_coll = 0
!   do j = 1, n_particle
!     if(j .ne. ip) then
!       cflag = 0
!       do i = 1, n_slice
!         xij = sum((q(i,j,:) - q(i,ip,:))**2)
!         if(xij < hard_sphere_radius_squared) then
!           cflag = 1
!           exit
!         end if
!       end do
!       if(cflag == 1) call ipush(clist,j,n_coll)
!     end if
!   end do
! 
! end subroutine get_collisions!}}}
!@-at
!@@c

!@-node:gcross.20090623152316.88:get_collisions
!@-node:gcross.20090623152316.84:Utility subroutines
!@+node:gcross.20090623152316.80:Move subroutines
!@+node:gcross.20090623152316.81:rigid
subroutine rigid_move(qin, qout, part_num, low_dim, high_dim, i0, i1, dnu, N_SLICE, N_PARTICLE, N_DIM)
  integer, intent(in) :: N_SLICE, N_PARTICLE, N_DIM
  double precision, dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qin
  double precision, dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qout
  integer :: part_num, low_dim, high_dim, i0, i1
  double precision :: dnu

  real, dimension( N_DIM )  :: nu
  integer :: j

  call random_number( nu )
  nu = ( nu - 0.5 ) * dnu

  do j = low_dim, high_dim
    qout(i0:i1,part_num,j) = qin(i0:i1,part_num,j) + nu(j)
  end do
end subroutine rigid_move!}}}
!@-node:gcross.20090623152316.81:rigid
!@+node:gcross.20090623152316.82:swap
subroutine swap_move(qin, qout, part_num, low_dim, high_dim, i0, i1, dnu, N_SLICE, N_PARTICLE, N_DIM)
  integer, intent(in) :: N_SLICE, N_PARTICLE, N_DIM
  double precision, dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qin
  double precision, dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qout
  integer :: part_num, low_dim, high_dim, i0, i1
  double precision :: dnu

  double precision, dimension( N_DIM )  :: nu
  integer :: j

  call random_number( nu )
  nu = ( nu - 0.5 ) * dnu

  do j = low_dim, high_dim
    qout(i0:i1,part_num,j) = -qin(i0:i1,part_num,j) + nu(j)
  end do
end subroutine swap_move!}}}
!@-node:gcross.20090623152316.82:swap
!@+node:gcross.20090623152316.83:cascading swap
subroutine cascading_swap_move(qin, qout, part_num, swap_dim, N_SLICE, N_PARTICLE, N_DIM, pbc_period_length)
  integer, intent(in) :: N_SLICE, N_PARTICLE, N_DIM
  double precision, dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qin
  double precision, dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qout
  integer :: part_num, swap_dim
  double precision, intent(in), optional :: pbc_period_length
  integer, dimension(N_PARTICLE) :: plist,alist,clist
  integer :: pcnt,acnt,t_acnt,ccnt
  integer :: j,k,tmp
  integer :: ip_swap
  double precision :: rnd

  call random_number( rnd )
  pcnt = 0
  acnt = 0
  call ipush(plist,part_num,pcnt)
  do j = 1, n_particle
    if(j .ne. part_num) then
      call ipush(alist,j,acnt)
    end if
  end do
  ip_swap = part_num
  do while(pcnt > 0)
    if( present(pbc_period_length) ) then
      qout(:,ip_swap,swap_dim) = -qin(:,ip_swap,swap_dim)+ rnd*pbc_period_length
    else
      qout(:,ip_swap,swap_dim) = -qin(:,ip_swap,swap_dim)
    end if
!@+at
!     ccnt = 0
!     call get_collisions(qout,ip_swap,clist,ccnt)
! #ifdef DEBUG_CASCADE_SWAP
!     print *, "n_coll = ",ccnt
! #endif
!     do j=1,ccnt
!       t_acnt = acnt
!       do k = t_acnt,1,-1
!         if(clist(j) == alist(k)) then
!           alist(k) = alist(acnt)
!           call ipop(alist,tmp,acnt)
!           call ipush(plist,clist(j),pcnt)
! #ifdef DEBUG_CASCADE_SWAP
!           print *, j,k,tmp,acnt
!           print *, j,k,clist(j),pcnt
! #endif
!         end if
!       end do
!     end do
!@-at
!@@c
    call ipop(plist,ip_swap,pcnt)
  end do
end subroutine cascading_swap_move!}}}
!@-node:gcross.20090623152316.83:cascading swap
!@+node:gcross.20090805153643.1851:brownian bridge
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
subroutine bbridge(q0, q1, pnum, low_dim, high_dim, btype, i0, i1, lambda, dtau, &
    N_SLICE, N_PARTICLE, N_DIM, &
    od_pnum, n_od_particle, &
    od_dim_low, od_dim_high &
  )
  integer, intent(in) :: N_SLICE, N_PARTICLE, N_DIM
  double precision, dimension( :, :, : ), intent(in):: q0
  double precision, dimension( :, :, : ), intent(inout):: q1
  integer, intent(in):: pnum, low_dim, high_dim, btype
  integer, intent(inout) :: i0
  integer, intent(inout) :: i1
  double precision, intent(in):: lambda
  double precision, intent(in):: dtau
  integer, intent(in), optional :: od_pnum, n_od_particle, od_dim_low, od_dim_high

  double precision :: dt

  integer :: i
  integer :: j,cnt
  integer :: s0,s1
  integer :: ndim
  double precision :: t1,t2,mean,sigma,sigmasq
  double precision :: a,b

  double precision, dimension( 2*ceiling((i1-i0+2.)*(high_dim-low_dim+1.)/2.) )::  nu

  integer :: cslice
  cslice = n_slice/2

  dt = lambda*dtau*2.0D0

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

  if( (i0 <= CSLICE) .and. (i1 >= CSLICE) ) then
    i1 = i1 + 1
    if(i1 .gt. N_SLICE) then
      i1 = N_SLICE
    end if

    if( (present(od_pnum) .and. (pnum == OD_PNUM)) .or. (present(N_OD_PARTICLE) .and. (pnum <= N_OD_PARTICLE)) ) then
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
!@-node:gcross.20090805153643.1851:brownian bridge
!@-node:gcross.20090623152316.80:Move subroutines
!@-others

end module sample
!@-node:gcross.20090623152316.73:@thin sample.f95
!@-leo
