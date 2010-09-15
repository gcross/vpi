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
    swap_dimension_low, swap_dimension_high, &
    particle_number, move_type, &
    n_slices, n_particles, n_dimensions, &
    off_diagonal_particle_number, prb_off_diagonal_particle_number, &
    n_off_diagonal_particles, &
    off_diagonal_dimension_low, off_diagonal_dimension_high &
  )
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  integer, parameter :: N_MOVE_TYPES = 3
  double precision, dimension ( N_MOVE_TYPES ), intent(in) :: move_type_probabilities, move_type_differentials
  double precision, dimension ( n_slices , n_particles, n_dimensions ), intent(in) :: q0
  double precision, dimension ( n_slices , n_particles, n_dimensions ), intent(inout) :: q1
  integer, intent(in) :: swap_dimension_low, swap_dimension_high
  integer, intent(out) :: move_start,move_end,particle_number,move_type
  double precision, intent(in), optional :: prb_off_diagonal_particle_number
  integer, intent(in), optional :: &
    off_diagonal_particle_number, &
    n_off_diagonal_particles, &
    off_diagonal_dimension_low, off_diagonal_dimension_high
  double precision, intent(in) :: dM, lambda

  double precision :: rnd
  integer :: i0,i1,center_slice
  integer :: swap_dimension
  center_slice = n_slices/2

  call random_number( rnd )
  if (present(off_diagonal_particle_number) .and. ( rnd .lt. prb_off_diagonal_particle_number ) ) then
    particle_number = off_diagonal_particle_number
  else
    particle_number = ceiling(rnd*n_particles)
  end if

  call random_number( rnd )
  if (rnd < move_type_probabilities(MT_BBRIDGE)) then
    move_type = MT_BBRIDGE
    call random_number( rnd )
    i0 = floor(rnd*((n_slices-1)-(2-dM))) - dM + 2
    i1 = i0 + dM
! print *, "i0, i1 ",i0,i1
    if( ((i0 .ge. 1) .and. (i1 .le. n_slices)) ) then
      call bbridge(q0, q1, particle_number, 1, n_dimensions, S_BRIDGE, i0, i1, lambda, &
            move_type_differentials(MT_BBRIDGE), n_slices, n_particles, n_dimensions)
      move_start = i0+1
      move_end = i1-1
    else
      if(i0 .lt. 1) then
        i0 = 1
        call bbridge(q0, q1, particle_number, 1, n_dimensions, LEFT_BRIDGE, i0, i1, lambda, &
            move_type_differentials(MT_BBRIDGE), n_slices, n_particles, n_dimensions)
        move_start = 1
        move_end = i1-1
      else
        if(i1 .gt. n_slices) then
          i1 = n_slices
          call bbridge(q0, q1, particle_number, 1, n_dimensions, RIGHT_BRIDGE, i0, i1, lambda, &
            move_type_differentials(MT_BBRIDGE), n_slices, n_particles, n_dimensions, &
            off_diagonal_particle_number, n_off_diagonal_particles, &
            off_diagonal_dimension_low, off_diagonal_dimension_high &
            )
          move_start = i0+1
          move_end = n_slices
        end if
      end if
    end if
  else 
! move all time steps of a single particle at once
    if (rnd < move_type_probabilities(MT_BBRIDGE) + move_type_probabilities(MT_RIGID)) then
      move_type = MT_RIGID
      call rigid_move(q0, q1, particle_number, 1, n_dimensions, 1, n_slices, &
              move_type_differentials(MT_RIGID), n_slices, n_particles, n_dimensions)
      if( (present(off_diagonal_particle_number) .and. &
              (particle_number .eq. off_diagonal_particle_number)) .or. &
          (present(n_off_diagonal_particles) .and. &
              (particle_number .le. n_off_diagonal_particles)) &
      ) then
        call rigid_move( &
                q0, q1, &
                particle_number, &
                off_diagonal_dimension_low, off_diagonal_dimension_high, &
                1, center_slice, &
                move_type_differentials(MT_RIGID), &
                n_slices, n_particles, n_dimensions &
              )
        call rigid_move( &
                q0, q1, &
                particle_number, &
                off_diagonal_dimension_low, off_diagonal_dimension_high, &
                center_slice+1, n_slices, &
                move_type_differentials(MT_RIGID), &
                n_slices, n_particles, n_dimensions &
              )
      end if
      move_start = 1
      move_end = n_slices
    else 
! attempt a rotation around a symmetry axis.  right now only does moves like q1 = -q0
    if (rnd < move_type_probabilities(MT_BBRIDGE) + move_type_probabilities(MT_RIGID) + move_type_probabilities(MT_SWAP) ) then
        move_type = MT_SWAP
        call random_number( rnd )
        swap_dimension = int(floor( rnd * (swap_dimension_high-swap_dimension_low+1) ))+swap_dimension_low
        if( (present(off_diagonal_particle_number) .and. &
                (particle_number == off_diagonal_particle_number)) .or. &
            (present(n_off_diagonal_particles) .and. &
                (particle_number <= n_off_diagonal_particles)) &
        ) then
          call random_number( rnd )
          if(rnd .ge. 0.5) then
            call swap_move(q0, q1, particle_number, swap_dimension, swap_dimension, 1, center_slice, &
                    move_type_differentials(MT_SWAP), n_slices, n_particles, n_dimensions)
          else
            call swap_move(q0, q1, particle_number, swap_dimension, swap_dimension, center_slice+1, n_slices, &
                    move_type_differentials(MT_SWAP), n_slices, n_particles, n_dimensions)
          end if
        else
          call cascading_swap_move(q0, q1, particle_number, swap_dimension, n_slices, n_particles, n_dimensions)
        end if
        move_start = 1
        move_end = n_slices
      end if
    end if
  end if

!  print *, "move_start, move_end ",move_start, move_end,particle_number
end subroutine sample_scheme1
!@nonl
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
! subroutine 
! get_collisions(q,ip,clist,n_coll,n_slices,n_particles,n_dimensions)
!   integer, intent(in) :: n_slices, n_particles, n_dimensions
!   double precision, dimension ( n_slices , n_particles, n_dimensions ) :: q
!   integer, dimension ( n_particles ) :: clist
!   integer :: ip,n_coll
!   integer :: i,j,cflag
!   double precision :: xij
! 
!   n_coll = 0
!   do j = 1, n_particles
!     if(j .ne. ip) then
!       cflag = 0
!       do i = 1, n_slices
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
!@nonl
!@-node:gcross.20090623152316.88:get_collisions
!@-node:gcross.20090623152316.84:Utility subroutines
!@+node:gcross.20090623152316.80:Move subroutines
!@+node:gcross.20090623152316.81:rigid
subroutine rigid_move(qin, qout, particle_number, low_dim, high_dim, i0, i1, dnu, n_slices, n_particles, n_dimensions)
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension ( n_slices , n_particles, n_dimensions ) :: qin
  double precision, dimension ( n_slices , n_particles, n_dimensions ) :: qout
  integer :: particle_number, low_dim, high_dim, i0, i1
  double precision :: dnu

  real, dimension( n_dimensions )  :: nu
  integer :: j

  call random_number( nu )
  nu = ( nu - 0.5 ) * dnu

  do j = low_dim, high_dim
    qout(i0:i1,particle_number,j) = qin(i0:i1,particle_number,j) + nu(j)
  end do
end subroutine rigid_move!}}}
!@nonl
!@-node:gcross.20090623152316.81:rigid
!@+node:gcross.20090623152316.82:swap
subroutine swap_move(qin, qout, particle_number, low_dim, high_dim, i0, i1, dnu, n_slices, n_particles, n_dimensions)
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension ( n_slices , n_particles, n_dimensions ) :: qin
  double precision, dimension ( n_slices , n_particles, n_dimensions ) :: qout
  integer :: particle_number, low_dim, high_dim, i0, i1
  double precision :: dnu

  double precision, dimension( n_dimensions )  :: nu
  integer :: j

  call random_number( nu )
  nu = ( nu - 0.5 ) * dnu

  do j = low_dim, high_dim
    qout(i0:i1,particle_number,j) = -qin(i0:i1,particle_number,j) + nu(j)
  end do
end subroutine swap_move!}}}
!@nonl
!@-node:gcross.20090623152316.82:swap
!@+node:gcross.20090623152316.83:cascading swap
subroutine cascading_swap_move(qin, qout, particle_number, swap_dimension, n_slices, n_particles, n_dimensions, pbc_period_length)
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension ( n_slices , n_particles, n_dimensions ) :: qin
  double precision, dimension ( n_slices , n_particles, n_dimensions ) :: qout
  integer :: particle_number, swap_dimension
  double precision, intent(in), optional :: pbc_period_length
  integer, dimension(n_particles) :: plist,alist
  integer :: pcnt,acnt
  integer :: j
  integer :: ip_swap
  double precision :: rnd

  call random_number( rnd )
  pcnt = 0
  acnt = 0
  call ipush(plist,particle_number,pcnt)
  do j = 1, n_particles
    if(j .ne. particle_number) then
      call ipush(alist,j,acnt)
    end if
  end do
  ip_swap = particle_number
  do while(pcnt > 0)
    if( present(pbc_period_length) ) then
      qout(:,ip_swap,swap_dimension) = -qin(:,ip_swap,swap_dimension)+ rnd*pbc_period_length
    else
      qout(:,ip_swap,swap_dimension) = -qin(:,ip_swap,swap_dimension)
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
!@nonl
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
subroutine bbridge(q0, q1, particle_number, low_dim, high_dim, btype, i0, i1, lambda, dtau, &
    n_slices, n_particles, n_dimensions, &
    off_diagonal_particle_number, n_off_diagonal_particles, &
    off_diagonal_dimension_low, off_diagonal_dimension_high &
  )
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension( n_slices, n_particles, n_dimensions ), intent(in) :: q0
  double precision, dimension( n_slices, n_particles, n_dimensions ), intent(inout) :: q1
  integer, intent(in):: particle_number, low_dim, high_dim, btype
  integer, intent(inout) :: i0
  integer, intent(inout) :: i1
  double precision, intent(in):: lambda
  double precision, intent(in):: dtau
  integer, intent(in), optional :: &
    off_diagonal_particle_number, &
    n_off_diagonal_particles, &
    off_diagonal_dimension_low, off_diagonal_dimension_high

  double precision :: dt

  integer :: i
  integer :: j,cnt
  double precision :: t1,t2,mean,sigma,sigmasq
  double precision :: a,b

  double precision, dimension( 2*ceiling((i1-i0+2.)*(high_dim-low_dim+1.)/2.) )::  nu

  integer :: center_slice
  center_slice = n_slices/2

  dt = lambda*dtau*2.0D0

  call ru_gasdev( nu )

  !print *, 'i0= ',i0,'i1 = ',i1,' btype= ',btype 
  !write (111,*) 'i0= ',i0,'i1 = ',i1,' btype= ',btype 
  q1(i0, particle_number, low_dim:high_dim) = q0(i0, particle_number, low_dim:high_dim)
  q1(i1, particle_number, low_dim:high_dim) = q0(i1, particle_number, low_dim:high_dim)

  cnt = 1

  if(btype .eq. LEFT_BRIDGE) then
    sigma= sqrt(dt*(i1-1))
    do j=low_dim, high_dim
      mean = q1(i1, particle_number, j)
      q1(i0, particle_number, j) = nu(cnt)*sigma + mean
      cnt = cnt + 1
    end do
  else if(btype .eq. RIGHT_BRIDGE) then
    sigma= sqrt(dt*(n_slices-i0))
    do j=low_dim, high_dim
      mean = q1(i0, particle_number, j)
      q1(i1, particle_number, j) = nu(cnt)*sigma + mean
      cnt = cnt + 1
    end do
  end if

  if( (i0 <= center_slice) .and. (i1 >= center_slice) ) then
    i1 = i1 + 1
    if(i1 .gt. n_slices) then
      i1 = n_slices
    end if

    if( (present(off_diagonal_particle_number) .and. (particle_number == off_diagonal_particle_number)) .or. &
        (present(n_off_diagonal_particles) .and. (particle_number <= n_off_diagonal_particles)) &
    ) then
      sigma= sqrt(dt*(center_slice-i0))
      do j = off_diagonal_dimension_low,off_diagonal_dimension_high
        mean = q1(i0, particle_number, j)
        q1(center_slice, particle_number, j) = nu(cnt)*sigma + mean
        cnt = cnt + 1
      end do
      sigma= sqrt(dt*(i1-(center_slice+1)))
      do j = off_diagonal_dimension_low,off_diagonal_dimension_high
        mean = q1(i1, particle_number, j)
        q1(center_slice+1, particle_number, j) = nu(cnt)*sigma + mean
        cnt = cnt + 1
      end do
    else
      t1 = dt*(center_slice-i0)
      t2 = dt*(i1-(center_slice+1))
      sigmasq=1./( 1./t1 + 1./t2)
      sigma= sqrt(sigmasq)
      do j = low_dim, high_dim
        a = q1(i0, particle_number, j)
        b = q1(i1, particle_number, j)
        mean=(t2*a  + t1*b)/(t1 + t2)
        q1(center_slice, particle_number, j) = nu(cnt)*sigma + mean
        q1(center_slice+1, particle_number, j) = q1(center_slice, particle_number, j)
        cnt = cnt + 1
      end do
    end if

    do i=i0+1, center_slice-1
      t1 = dt
      t2 = dt*(center_slice-i)
      sigmasq=1./( 1./t1 + 1./t2)
      sigma= sqrt(sigmasq)
      do j=low_dim, high_dim
        a = q1(i-1, particle_number, j)
        b = q1(center_slice, particle_number, j)
        mean=(t2*a  + t1*b)/(t1 + t2)
        q1(i, particle_number, j) = nu(cnt)*sigma + mean
        cnt = cnt + 1
      end do
    end do

    do i=center_slice+2, i1-1
      t1 = dt
      t2 = dt*(i1-i)
      sigmasq=1./( 1./t1 + 1./t2)
      sigma= sqrt(sigmasq)
      do j=low_dim, high_dim
        a = q1(i-1, particle_number, j)
        b = q1(i1, particle_number, j)
        mean=(t2*a  + t1*b)/(t1 + t2)
        q1(i, particle_number, j) = nu(cnt)*sigma + mean
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
        a = q1(i-1, particle_number, j)
        b = q1(i1, particle_number, j)
        mean=(t2*a  + t1*b)/(t1 + t2)
        q1(i, particle_number, j) = nu(cnt)*sigma + mean
        cnt = cnt + 1
      end do
    end do
  end if
end subroutine bbridge
!@nonl
!@-node:gcross.20090805153643.1851:brownian bridge
!@-node:gcross.20090623152316.80:Move subroutines
!@-others

end module sample
!@-node:gcross.20090623152316.73:@thin sample.f95
!@-leo
