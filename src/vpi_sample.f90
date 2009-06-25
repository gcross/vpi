!@+leo-ver=4-thin
!@+node:gcross.20090623152316.73:@thin vpi_sample.f90
!@@language fortran90
module vpi_sample

 !@ << Imported modules >>
 !@+node:gcross.20090623152316.74:<< Imported modules >>
 use vpi_rand_utils
 use vpi_bbridge
 use vpi_bisect
 use vpi_defines
 !@-node:gcross.20090623152316.74:<< Imported modules >>
 !@nl

 implicit none

 !@ << Move type constants >>
 !@+node:gcross.20090623152316.75:<< Move type constants >>
 integer, parameter :: MT_BBRIDGE = 0
 integer, parameter :: MT_RIGID = 1
 integer, parameter :: MT_SWAP = 2
 !@nonl
 !@-node:gcross.20090623152316.75:<< Move type constants >>
 !@nl

 contains

!@+others
!@+node:gcross.20090623152316.76:Sampling subroutines
!@+node:gcross.20090623152316.77:scheme 1
subroutine sample_scheme1(q0, q1, move_start, move_end, part_num, mtype)
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM ), intent(in) :: q0
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM ), intent(inout) :: q1
  integer, intent(out) :: move_start,move_end,part_num,mtype

  real(kind=b8) :: rnd
  integer :: i0,i1
  integer :: swap_dim

  call random_number( rnd )
  if (eval_off_diagonal .and. ( rnd .lt. PROB_od_pnum ) ) then
    part_num = od_pnum
  else
    part_num = ceiling(rnd*N_PARTICLE)
  end if

  call random_number( rnd )
  if (rnd < PROB_BBRIDGE_MOVE) then
    mtype = MT_BBRIDGE
    call random_number( rnd )
    i0 = floor(rnd*((N_SLICE-1)-(2-dM))) - dM + 2
    i1 = i0 + dM
! print *, "i0, i1 ",i0,i1
    if( ((i0 .ge. 1) .and. (i1 .le. N_SLICE)) ) then
      call bbridge(q0, q1, part_num, 1, N_DIM, S_BRIDGE, i0, i1, lambda, dtau)
      move_start = i0+1
      move_end = i1-1
    else
      if(i0 .lt. 1) then
        i0 = 1
        call bbridge(q0, q1, part_num, 1, N_DIM, LEFT_BRIDGE, i0, i1, lambda, dtau)
        move_start = 1
        move_end = i1-1
      else
        if(i1 .gt. N_SLICE) then
          i1 = N_SLICE
          call bbridge(q0, q1, part_num, 1, N_DIM, RIGHT_BRIDGE, i0, i1, lambda, dtau)
          move_start = i0+1
          move_end = N_SLICE
        end if
      end if
    end if
  else 
! move all time steps of a single particle at once
    if (rnd < PROB_BBRIDGE_MOVE + PROB_RIGID_MOVE) then
      mtype = MT_RIGID
      call rigid_move(q0, q1, part_num, 1, N_DIM, 1, N_SLICE, dnu)
      if( ((part_num .eq. od_pnum) .and. (eval_off_diagonal)) .or. ((part_num .le. N_OD_PARTICLE) .and. (eval_nrdm)) ) then
        call rigid_move(q0, q1, part_num, OD_DIM_LOW, OD_DIM_HIGH, 1, CSLICE, dnu)
        call rigid_move(q0, q1, part_num, OD_DIM_LOW, OD_DIM_HIGH, CSLICE+1, N_SLICE, dnu)
      end if
      move_start = 1
      move_end = N_SLICE
    else 
! attempt a rotation around a symmetry axis.  right now only does moves like q1 = -q0
      if (rnd < PROB_BBRIDGE_MOVE + PROB_RIGID_MOVE + PROB_SWAP_MOVE) then
        mtype = MT_SWAP
        call random_number( rnd )
        if( swap_in_12 ) then
          swap_dim = int(floor( rnd * 2 ))+1
        else
          if( swap_in_123 ) then
            swap_dim = int(floor( rnd * 3 ))+1
          else
            swap_dim = 3
          end if
        end if
        if( ((part_num .eq. od_pnum) .and. (eval_off_diagonal)) .or. ((part_num .le. N_OD_PARTICLE) .and. (eval_nrdm)) ) then
          call random_number( rnd )
          if(rnd .ge. 0.5) then
            call swap_move(q0, q1, part_num, swap_dim, swap_dim, 1, CSLICE, dswap)
          else
            call swap_move(q0, q1, part_num, swap_dim, swap_dim, CSLICE+1, N_SLICE, dswap)
          end if
        else
          call cascading_swap_move(q0, q1, part_num, swap_dim)
        end if
        move_start = 1
        move_end = N_SLICE
      end if
    end if
  end if

!  print *, "move_start, move_end ",move_start, move_end,part_num
end subroutine sample_scheme1
!@-node:gcross.20090623152316.77:scheme 1
!@+node:gcross.20090623152316.78:rotation
subroutine sample_rotation(q_rot0, q_rot1, move_start, move_end, part_num)
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM_ROT ), intent(in) :: q_rot0
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM_ROT ), intent(inout) :: q_rot1
  integer, intent(in) :: move_start,move_end,part_num

  real(kind=b8), dimension(N_DIM_ROT) :: n_vect
  real(kind=b8) :: nu
  real(kind=b8) :: ctheta,phi,stheta,sint0,cost0,sinp0,cosp0
  real(kind=b8), dimension(3) :: e1
  real(kind=b8), dimension(3,3) :: Ryz,Ry,Rz

  integer :: i

  do i=move_start,move_end
    cost0=q_rot0(i,part_num,3)
    if(abs(1d0-abs(cost0))<1d-10) then
      sint0=0d0; sinp0=0d0; cosp0=1d0
    else
      sint0=sqrt(abs(1d0-cost0**2))
      sinp0=q_rot0(i,part_num,2)/sint0
      cosp0=q_rot0(i,part_num,1)/sint0
    endif
    Rz=reshape( (/ cosp0,sinp0,0.0_b8 , -sinp0,cosp0,0.0_b8 , 0.0_b8,0.0_b8,1.0_b8/), (/3,3/) )
    Ry=reshape( (/ cost0,0.0_b8,-sint0, 0.0_b8,1.0_b8,0.0_b8 , sint0,0.0_b8,cost0/) , (/3,3/) )
!    Ryz = matmul(Ry,Rz)
    Ryz = matmul(Rz,Ry)
e1=matmul(q_rot0(i,part_num,:),Ryz)
if( abs(e1(1))>1d-6 .or. abs(e1(2))>1d-6 .or. abs(1d0-e1(3))>1d-6 )then
write(*,*)'sample_rotation()::ERROR ',e1
stop
endif
    call random_number( nu )
    ctheta = 1.0_b8-nu*drot
    stheta = sqrt(abs(1d0-ctheta**2))
    call random_number( nu )
    phi = M_2PI*nu 
    e1(1) = cos(phi)*stheta
    e1(2) = sin(phi)*stheta
    e1(3) = ctheta
!    q_rot1(i,part_num,:) = matmul(e1,Ryz)
    q_rot1(i,part_num,:) = matmul(Ryz,e1)
  end do

  if ( (eval_off_diagonal .eqv. .false.) .or. (part_num .ne. od_pnum ) )then
    q_rot1(CSLICE+1,part_num,:) = q_rot1(CSLICE,part_num,:)
  end if

end subroutine sample_rotation!}}}
!@-node:gcross.20090623152316.78:rotation
!@+node:gcross.20090623152316.79:rotation xx (?)
subroutine xxsample_rotation(q_rot0, q_rot1, move_start, move_end, part_num)
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM_ROT ), intent(in) :: q_rot0
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM_ROT ), intent(inout) :: q_rot1
  integer, intent(in) :: move_start,move_end,part_num

  real(kind=b8), dimension(N_DIM_ROT) :: n_vect
  real(kind=b8), dimension(N_DIM_ROT) :: nu
  real(kind=b8) :: norm

  integer :: i

  do i=move_start,move_end
    call ru_gasdev( nu )
    norm = -1
    do while( norm .le. 0.0 )
      n_vect = q_rot0(i,part_num,:) + drot*nu(:)
      norm =  sqrt(dot_product(n_vect,n_vect))
    end do
    q_rot1(i,part_num,:) = n_vect(:)/norm
  end do

  if ( (eval_off_diagonal .eqv. .false.) .or. (part_num .ne. od_pnum ) )then
    q_rot1(CSLICE+1,part_num,:) = q_rot1(CSLICE,part_num,:)
  end if

end subroutine xxsample_rotation!}}}
!@-node:gcross.20090623152316.79:rotation xx (?)
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
subroutine get_collisions(q,ip,clist,n_coll)
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: q
  integer, dimension ( N_PARTICLE ) :: clist
  integer :: ip,n_coll
  integer :: i,j,cflag
  real(kind=b8) :: xij

  n_coll = 0
  do j = 1, n_particle
    if(j .ne. ip) then
      cflag = 0
      do i = 1, n_slice
        xij = sum((q(i,j,:) - q(i,ip,:))**2) 
        if(xij < hard_sphere_radius_squared) then
          cflag = 1
          exit
        end if
      end do
      if(cflag == 1) call ipush(clist,j,n_coll)
    end if
  end do

end subroutine get_collisions!}}}
!@-node:gcross.20090623152316.88:get_collisions
!@-node:gcross.20090623152316.84:Utility subroutines
!@+node:gcross.20090623152316.80:Move subroutines
!@+node:gcross.20090623152316.81:rigid
subroutine rigid_move(qin, qout, part_num, low_dim, high_dim, i0, i1, dnu)
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qin
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qout
  integer :: part_num, low_dim, high_dim, i0, i1
  real(kind=b8) :: dnu

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
subroutine swap_move(qin, qout, part_num, low_dim, high_dim, i0, i1, dnu)
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qin
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qout
  integer :: part_num, low_dim, high_dim, i0, i1
  real(kind=b8) :: dnu

  real(kind=b8), dimension( N_DIM )  :: nu
  integer :: j

  call random_number( nu )
  nu = ( nu - 0.5 ) * dnu

  do j = low_dim, high_dim
    qout(i0:i1,part_num,j) = -qin(i0:i1,part_num,j) + nu(j)
  end do
end subroutine swap_move!}}}
!@-node:gcross.20090623152316.82:swap
!@+node:gcross.20090623152316.83:cascading swap
subroutine cascading_swap_move(qin, qout, part_num, swap_dim)
  implicit none
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qin
  real(kind=b8), dimension ( N_SLICE , N_PARTICLE, N_DIM ) :: qout
  integer :: part_num, swap_dim
  integer, dimension(N_PARTICLE) :: plist,alist,clist
  integer :: pcnt,acnt,t_acnt,ccnt
  integer :: j,k,tmp
  integer :: ip_swap
  real(kind=b8) :: rnd

  call random_number( rnd )
!@@raw
#ifdef DEBUG_CASCADE_SWAP
!@@end_raw
  print *, "** starting cascading swap"
!@@raw
#endif
!@@end_raw
  pcnt = 0
  acnt = 0
  call ipush(plist,part_num,pcnt)
  do j = 1, n_particle
    if(j .ne. part_num) then
      call ipush(alist,j,acnt)
!@@raw
#ifdef DEBUG_CASCADE_SWAP
!@@end_raw
      print *, "alist = ",j," acnt = ",acnt
!@@raw
#endif
!@@end_raw
    end if
  end do
  ip_swap = part_num
  do while(pcnt > 0)
!@@raw
#ifdef DEBUG_CASCADE_SWAP
!@@end_raw
    print *, "p_swap = ",ip_swap," pcnt = ",pcnt
!@@raw
#endif
!@@end_raw
    if(use_pbc) then
      qout(:,ip_swap,swap_dim) = -qin(:,ip_swap,swap_dim)+ rnd*p_pbc_l
    else
      qout(:,ip_swap,swap_dim) = -qin(:,ip_swap,swap_dim)
    end if
    ccnt = 0
    call get_collisions(qout,ip_swap,clist,ccnt)
!@@raw
#ifdef DEBUG_CASCADE_SWAP
!@@end_raw
    print *, "n_coll = ",ccnt
!@@raw
#endif
!@@end_raw
    do j=1,ccnt
      t_acnt = acnt
      do k = t_acnt,1,-1
        if(clist(j) == alist(k)) then
          alist(k) = alist(acnt)
          call ipop(alist,tmp,acnt)
          call ipush(plist,clist(j),pcnt)
!@@raw
#ifdef DEBUG_CASCADE_SWAP
!@@end_raw
          print *, j,k,tmp,acnt
          print *, j,k,clist(j),pcnt
!@@raw
#endif
!@@end_raw
        end if
      end do
    end do
    call ipop(plist,ip_swap,pcnt)
  end do
end subroutine cascading_swap_move!}}}
!@-node:gcross.20090623152316.83:cascading swap
!@-node:gcross.20090623152316.80:Move subroutines
!@-others

end module vpi_sample
!@-node:gcross.20090623152316.73:@thin vpi_sample.f90
!@-leo
