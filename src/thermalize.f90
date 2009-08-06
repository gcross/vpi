!@+leo-ver=4-thin
!@+node:gcross.20090805153643.1843:@thin thermalize.f90
!@@language fortran90
module thermalize

use vpi_xij
use vpi_sample

implicit none

contains 

!@+others
!@+node:gcross.20090623152316.32:accept_path
function accept_path( lngfn0, lngfn1 ) result( accept )
  double precision, intent(in) :: lngfn0, lngfn1
  logical :: accept

  real :: Pa, Ptest 

  Pa = exp(lngfn1 - lngfn0) 
  call random_number( Ptest )
  accept = ( Pa > Ptest )

end function accept_path
!@-node:gcross.20090623152316.32:accept_path
!@+node:gcross.20090721121051.1746:compute_log_acceptance_weight
function compute_log_acceptance_weight (&
! INPUT: particle position / rotation information
    x, xij2, x_rot, &
! INPUT: path slice to consider
    move_start, move_end, &
! INPUT: 
    nslice, np, ndim, &
! OUTPUT: potential and square gradients
    U, gradU2, &
! OUTPUT: whether the move should be rejected outright
    reject_flag &
! OUTPUT: computed weight
    ) result ( weight )

!@+at
! Array dimensions / slicing
!@-at
!@@c
  integer, intent(in) :: move_start, move_end, nslice, np, ndim

!@+at
! Function input
!@-at
!@@c

  real(kind=b8), dimension ( nslice, np , ndim ), intent(in) :: x
  real(kind=b8), dimension ( nslice, np , np ), intent(in) :: xij2
  real(kind=b8), dimension ( nslice, np , N_DIM_ROT ), intent(in) :: x_rot

!@+at
! Function output
!@-at
!@@c
  real(kind=b8), dimension( nslice, np ), target, intent(out) :: U
  real(kind=b8), dimension( nslice ), intent(out) :: gradU2
  real(kind=b8) :: weight
  logical, intent(out) :: reject_flag

!@+at
! Temporary variables.
!@-at
!@@c
  real (kind=b8), dimension( nslice, np ), target :: U_effective
  real (kind=b8), dimension( :, : ), pointer :: U_total
  real (kind=b8) :: lngfn, rotgfn, hsgfn, lntfn
  integer :: sl_start, sl_end
!@+at
! Code begins:
!@-at
!@@c
  !@  << Compute contribution from potentials >>
  !@+node:gcross.20090626112946.1696:<< Compute contribution from potentials >>
  call compute_physical_potential (&
      x, xij2, x_rot, &
      move_start, move_end, &
      nslice, np, ndim, &
      U, gradU2, &
      reject_flag &
      )

  if (reject_flag) then
    return
  endif

  if (fixed_angular_momentum > 0) then
    call compute_effective_rotational_potential (&
      x, &
      move_start, move_end, &
      nslice, np, ndim, &
      U &
      )
  end if

  U_total => U

  !@-node:gcross.20090626112946.1696:<< Compute contribution from potentials >>
  !@nl

  !@  << Compute contribution from propagators >>
  !@+node:gcross.20090721121051.1752:<< Compute contribution from propagators >>
  if(move_start .le. 1) then
    sl_start = 1
  else
    sl_start = move_start-1
  end if
  if(move_end .ge. N_SLICE) then
    sl_end = N_SLICE
  else
    sl_end = move_end+1
  end if

  if( n_slice > 2 ) then 
    if ( use_gfn4 ) then
      lngfn = vpi_gfn4_sp( sl_start, sl_end, part_num, U_total, gradU2, U_weight, gU2_weight, &
                           N_SLICE, N_PARTICLE, N_DIM, lambda, dtau ) 
    else 
      lngfn = vpi_gfn2_sp( sl_start, sl_end, part_num, U_total, N_SLICE, N_PARTICLE, N_DIM, dtau ) 
    end if

  !@+at
  !   if ( use_HS_gfn ) then
  ! #ifdef USE_IMAGE_HSGFN
  !     hsgfn = vpi_hs_gfn( sl_start, sl_end, part_num, xij2, N_SLICE, 
  ! N_PARTICLE, N_DIM, dtau*lambda )
  ! #else
  !     hsgfn = vpi_hs_gfn2( sl_start, sl_end, part_num, x, N_SLICE, 
  ! N_PARTICLE, N_DIM, dtau*lambda )
  ! #endif
  !     if( hsgfn <= 0 .and. hsgfn1 > 0 ) then
  !       print *, "****ERROR***** hsgfn < 0 ", hsgfn0, hsgfn1
  !     else
  !       lngfn = lngfn + log(hsgfn)
  !     else
  !   end if
  !@-at
  !@@c

    if ( use_HW_gfn ) then
      hsgfn = vpi_hw_gfn( sl_start, sl_end, part_num, x, N_SLICE, N_PARTICLE, N_DIM, dtau )
      lngfn = lngfn + log(hsgfn)
    end if

    if ( eval_rotation ) then
      rotgfn =  gfn_rot(x_rot,part_num,sl_start,sl_end,dtau)
  !          write (1000+my_rank,*) rotgfn
      lngfn = lngfn + log(rotgfn)
    end if
  end if

  !@-node:gcross.20090721121051.1752:<< Compute contribution from propagators >>
  !@nl

  !@  << Compute contribution from trial functions >>
  !@+node:gcross.20090721121051.1754:<< Compute contribution from trial functions >>
  lntfn = tfunc(x, 1, N_SLICE, N_PARTICLE, N_DIM) + &
          tfunc(x, N_SLICE, N_SLICE, N_PARTICLE, N_DIM) + &
          jas_tfun(x, xij2, 1, N_SLICE, N_PARTICLE, N_DIM) + &
          jas_tfun(x, xij2, N_SLICE, N_SLICE, N_PARTICLE, N_DIM)
  !@-node:gcross.20090721121051.1754:<< Compute contribution from trial functions >>
  !@nl

  weight = lngfn + lntfn

  return

end function
!@-node:gcross.20090721121051.1746:compute_log_acceptance_weight
!@+node:gcross.20090626112946.1694:thermalize_path
subroutine thermalize_path( &
  q, xij2, U, gradU2, &
  N_SLICES, N_PARTICLES, N_DIMENSIONS, &
  n_trials, &
  slice_move_attempted_counts, move_type_attempted_counts, &
  slice_move_accepted_counts, move_type_accepted_counts, &  
  q_rot, N_DIMENSIONS_ROTATION, &
  pbc_period_length &
  )

  ! Input variables
  integer, intent(in) :: N_SLICES, N_PARTICLES, N_DIMENSIONS, n_trials
  integer, intent(in), optional :: N_DIMENSIONS_ROTATION
  double precision, intent(in), optional :: pbc_period_length

  ! Updated variables
  double precision, dimension( N_SLICES, N_PARTICLES, N_DIMENSIONS ), intent(inout) :: q
  double precision, dimension( :, :, : ), intent(inout), optional :: q_rot
  double precision, dimension( N_SLICES, N_PARTICLES, N_PARTICLES ), intent(inout) :: xij2
  double precision, dimension( N_SLICES, N_PARTICLES ), intent(inout) :: U
  double precision, dimension( N_SLICES ), intent(inout) :: gradU2
  integer, dimension( N_SLICES ), intent(inout) :: slice_move_attempted_counts, slice_move_accepted_counts
  integer, dimension( N_MOVE_TYPES ), intent(inout) :: move_type_attempted_counts, move_type_accepted_counts

  ! Local variables
  integer :: i, mtype, move_start, move_end, part_num
  double precision :: old_weight, new_weight
  logical :: reject_flag
  double precision, dimension( N_SLICES, N_PARTICLES, N_DIMENSIONS ) :: q_trial
  double precision, dimension( :, :, : ), allocatable :: q_rot_trial
  double precision, dimension( N_SLICES, N_PARTICLES, N_PARTICLES ), intent(inout) :: xij2_trial
  double precision, dimension( N_SLICES, N_PARTICLES ) :: U_trial
  double precision, dimension( N_SLICES ) :: gradU2_trial

  !@  << Initialize trial functions >>
  !@+node:gcross.20090805153643.1848:<< Initialize trial functions >>
  q_trial = q
  if( present(q_rot) ) then
    allocate( q_rot_trial( N_SLICES, N_PARTICLES, N_DIMENSIONS ) )
    q_rot_trial = q_rot
  end if
  xij2_trial = xij2
  U_trial = U
  gradU2_trial = gradU2
  !@nonl
  !@-node:gcross.20090805153643.1848:<< Initialize trial functions >>
  !@nl

  do i = 1, n_trials
    !@    << Randomly choose a move to make and apply it to the system >>
    !@+node:gcross.20090626112946.1688:<< Randomly choose a move to make and apply it to the system >>
    call sample_scheme1(q, q_trial, move_start, move_end, part_num, mtype)
    if( present(q_rot) ) then
      call sample_rotation(q_rot, q_rot_trial, move_start, move_end, part_num)
    end if

    !@<< Impose periodic boundary conditions >>
    !@+node:gcross.20090626112946.1691:<< Impose periodic boundary conditions >>
    if( present(pbc_period_length) ) then
      q_trial(move_start:move_end,:,:) = wrap_around(q_trial(move_start:move_end,:,:))
    end if
    !@-node:gcross.20090626112946.1691:<< Impose periodic boundary conditions >>
    !@nl

    !@<< Update move type statistics >>
    !@+node:gcross.20090626112946.1690:<< Update move type statistics >>
    slice_move_attempted_counts(move_start:move_end) =  slice_move_attempted_counts(move_start:move_end) + 1
    move_type_attempted_counts(mtype) = move_type_attempted_counts(mtype) + 1
    !@-node:gcross.20090626112946.1690:<< Update move type statistics >>
    !@nl

    !@<< Update the displacement matrix >>
    !@+node:gcross.20090626112946.1689:<< Update the displacement matrix >>
    if( present(pbc_period_length) ) then
      call vpi_update_xij_pbc( xij2_trial, q_trial, move_start, move_end, N_SLICES, N_PARTICLES, N_DIMENSIONS  )
    else
      call vpi_update_xij( xij2_trial, q_trial, move_start, move_end, N_SLICES, N_PARTICLES, N_DIMENSIONS  )
    end if
    !@-node:gcross.20090626112946.1689:<< Update the displacement matrix >>
    !@nl
    !@-node:gcross.20090626112946.1688:<< Randomly choose a move to make and apply it to the system >>
    !@nl

    !@    << Determine whether the move should be accepted >>
    !@+node:gcross.20090626112946.1693:<< Determine whether the move should be accepted >>
    !@<< Compute logarithmic probability of acceptance >>
    !@+node:gcross.20090721121051.1755:<< Compute logarithmic probability of acceptance >>
    old_weight = compute_log_acceptance_weight (&
        q, xij2, q_rot_trial, &
        move_start, move_end, &
        n_slices, n_particles, n_dimensions, &
        U, gradU2, &
        reject_flag &
        )

    if(reject_flag) then
      stop "Error:  A path that had been okay before is somehow invalid now."
    end if

    new_weight = compute_log_acceptance_weight (&
        q_trial, xij2_trial, q_rot_trial, &
        move_start, move_end, &
        n_slices, n_particles, n_dimensions, &
        U_trial, gradU2_trial, &
        reject_flag &
        )

    if(reject_flag) then
      call revert_move
      cycle
    end if
    !@-node:gcross.20090721121051.1755:<< Compute logarithmic probability of acceptance >>
    !@nl

    !@<< Accept or reject the move >>
    !@+node:gcross.20090626112946.1706:<< Accept or reject the move >>
    if ( accept_path( old_weight, new_weight, 1.0D0 ) ) then
      !@  << Update the degrees of freedom >>
      !@+node:gcross.20090626112946.1701:<< Update the degrees of freedom >>
      q(move_start:move_end,:,:) = q_trial(move_start:move_end,:,:)
      if( present(q_rot) ) then
        q_rot(move_start:move_end,:,:) = q_rot_trial(move_start:move_end,:,:)
      end if
      !@-node:gcross.20090626112946.1701:<< Update the degrees of freedom >>
      !@nl
      !@  << Update derived quantities >>
      !@+node:gcross.20090805093617.1831:<< Update derived quantities >>
      xij2(move_start:move_end,:,:) = xij2_trial(move_start:move_end,:,:) 
      U(move_start:move_end,:) = U_trial(move_start:move_end,:)
      gradU2(move_start:move_end) = gradU2_trial(move_start:move_end)
      !@-node:gcross.20090805093617.1831:<< Update derived quantities >>
      !@nl
      !@  << Update statistics >>
      !@+node:gcross.20090626112946.1703:<< Update statistics >>
      slice_move_accepted_counts(move_start:move_end) =  slice_move_accepted_counts(move_start:move_end) + 1
      move_type_accepted_counts(mtype) = move_type_accepted_counts(mtype) + 1
      !@-node:gcross.20090626112946.1703:<< Update statistics >>
      !@nl
    else
      call revert_move
    end if
    !@-node:gcross.20090626112946.1706:<< Accept or reject the move >>
    !@nl

    !@-node:gcross.20090626112946.1693:<< Determine whether the move should be accepted >>
    !@nl
  end do

  if( allocated( q_rot_trial ) ) then
    deallocate( q_rot_trial )
  end if

  contains

  !@  @+others
  !@+node:gcross.20090805093617.1835:subroutine revert_move
  subroutine revert_move
    q_trial(move_start:move_end,:,:) = q(move_start:move_end,:,:)
    if( eval_rotation ) then
      q_rot_trial(move_start:move_end,:,:) = q_rot(move_start:move_end,:,:)
    end if
    xij2_trial(move_start:move_end,:,:) = xij2(move_start:move_end,:,:)
    U_trial(move_start:move_end,:) = U(move_start:move_end,:)
    gradU2_trial(move_start:move_end) =  gradU2(move_start:move_end)
  end subroutine revert_move
  !@-node:gcross.20090805093617.1835:subroutine revert_move
  !@-others

end subroutine thermalize_path
!@-node:gcross.20090626112946.1694:thermalize_path
!@-others

end module thermalize
!@-node:gcross.20090805153643.1843:@thin thermalize.f90
!@-leo
