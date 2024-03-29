!@+leo-ver=4-thin
!@+node:gcross.20090805153643.1843:@thin thermalize.f95
!@@language fortran90
module thermalize

use xij
use sample
use gfn
use angular_momentum

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
!@+node:gcross.20090626112946.1694:thermalize_path
subroutine thermalize_path( &
  q, xij2, U, gradU2, &
  n_slices, n_particles, n_dimensions, &
  n_trials, &
  move_type_probabilities, move_type_differentials, &
  dM, lambda, &
  swap_dimension_low, swap_dimension_high, &
  slice_move_attempted_counts, move_type_attempted_counts, &
  slice_move_accepted_counts, move_type_accepted_counts, &
  compute_potential, trial_function, greens_function, &
  pbc_period_length, &
  off_diagonal_particle_number, prb_off_diagonal_particle_number, &
  n_off_diagonal_particles, &
  off_diagonal_dimension_low, off_diagonal_dimension_high &
  )

!@+at
! Input variables
!@-at
!@@c
  integer, intent(in) :: n_slices, n_particles, n_dimensions, n_trials
  double precision, dimension(N_MOVE_TYPES), intent(in) :: &
    move_type_probabilities, move_type_differentials
  double precision, intent(in) :: dM, lambda
  integer, intent(in) :: swap_dimension_low, swap_dimension_high
  double precision, intent(in), optional :: pbc_period_length
  double precision, intent(in), optional :: prb_off_diagonal_particle_number
  integer, intent(in), optional :: &
    off_diagonal_particle_number, &
    n_off_diagonal_particles, &
    off_diagonal_dimension_low, off_diagonal_dimension_high

!@+at
! Potential and trial functions
!@-at
!@@c
!f2py external, intent(callback) :: compute_potential, trial_function, greens_function
interface
  !@  << Potential callback interface >>
  !@+middle:gcross.20090817102318.2271:Interface
  !@+node:gcross.20090809223137.1724:<< Potential callback interface >>
  subroutine compute_potential( x, xij2, n_slices, n_particles, n_dimensions, U, gradU2, reject_flag )
    double precision, dimension ( n_dimensions, n_particles, n_slices ), intent(in) :: x
    double precision, dimension ( n_particles, n_particles, n_slices ), intent(in) :: xij2
    integer, intent(in) :: n_slices, n_particles, n_dimensions
    double precision, dimension( n_particles, n_slices ), intent(out) :: U
    double precision, dimension( n_slices ), intent(out) :: gradU2
    logical, intent(out) :: reject_flag
  end subroutine compute_potential
  !@nonl
  !@-node:gcross.20090809223137.1724:<< Potential callback interface >>
  !@-middle:gcross.20090817102318.2271:Interface
  !@nl
  !@  << Trial callback interface >>
  !@+middle:gcross.20090817102318.2271:Interface
  !@+node:gcross.20090812093015.1848:<< Trial callback interface >>
  function trial_function( x, xij2, n_particles, n_dimensions, reject_flag ) result ( log_probability )
    integer, intent(in) :: n_particles, n_dimensions
    double precision, dimension( n_dimensions, n_particles ), intent(in) :: x
    double precision, dimension( n_particles, n_particles ), intent(in) :: xij2
    logical, intent(out) :: reject_flag
    double precision  :: log_probability
  end function trial_function
  !@-node:gcross.20090812093015.1848:<< Trial callback interface >>
  !@-middle:gcross.20090817102318.2271:Interface
  !@nl
  !@  << Green's function callback interface >>
  !@+middle:gcross.20090817102318.2271:Interface
  !@+node:gcross.20090828095451.1676:<< Green's function callback interface >>
  function greens_function( &
      x, xij2, U, gradU2, &
      lambda, dt, &
      slice_start, slice_end, &
      particle_number, &
      n_slices, n_particles, n_dimensions &
    ) result ( log_probability )
    integer, intent(in) :: n_slices, n_particles, n_dimensions
    double precision, dimension( n_dimensions, n_particles, n_slices ), intent(in) :: x
    double precision, dimension( n_particles, n_particles, n_slices ), intent(in) :: xij2
    double precision, dimension( n_particles, n_slices ), intent(in) :: U
    double precision, dimension( n_slices ), intent(in) :: gradU2
    double precision, intent(in) :: lambda, dt
    integer, intent(in) :: slice_start, slice_end
    integer, intent(in) :: particle_number
    double precision  :: log_probability
  end function
  !@nonl
  !@-node:gcross.20090828095451.1676:<< Green's function callback interface >>
  !@-middle:gcross.20090817102318.2271:Interface
  !@nl
end interface

!@+at
! Updated trial functions
!@-at
!@@c
  double precision, dimension(n_dimensions,n_particles,n_slices), intent(inout) :: q
  double precision, dimension(n_particles,n_particles,n_slices), intent(inout) :: xij2
  double precision, dimension(n_particles,n_slices), intent(inout) :: U
  double precision, dimension(n_slices), intent(inout) :: gradU2
  integer, dimension(n_slices), intent(inout) :: slice_move_attempted_counts, slice_move_accepted_counts
  integer, dimension(N_MOVE_TYPES), intent(inout) :: move_type_attempted_counts, move_type_accepted_counts

!@+at
! Local variables
!@-at
!@@c
  integer :: i, move_type, move_start, move_end, particle_number
  double precision :: old_weight, new_weight, dtau
  logical :: reject_flag, system_is_periodic, some_particles_are_off_diagonal
  double precision, dimension(n_dimensions,n_particles,n_slices) :: q_trial
  double precision, dimension(n_particles,n_particles,n_slices) :: xij2_trial
  double precision, dimension(n_particles,n_slices) :: U_trial
  double precision, dimension(n_slices) :: gradU2_trial

  system_is_periodic = present(pbc_period_length) .and. pbc_period_length > 0
  some_particles_are_off_diagonal = present(off_diagonal_particle_number) .and. off_diagonal_particle_number > 0

  !@  << Initialize trial functions >>
  !@+node:gcross.20090805153643.1848:<< Initialize trial functions >>
  q_trial = q
  xij2_trial = xij2
  U_trial = U
  gradU2_trial = gradU2
  !@nonl
  !@-node:gcross.20090805153643.1848:<< Initialize trial functions >>
  !@nl

  ! for convenience
  dtau = move_type_differentials(MT_BBRIDGE)

  do i = 1, n_trials
    !@    << Randomly choose a move to make and apply it to the system >>
    !@+node:gcross.20090626112946.1688:<< Randomly choose a move to make and apply it to the system >>
    if( some_particles_are_off_diagonal ) then
      call sample_scheme1(q, q_trial, move_start, move_end, &
          move_type_probabilities, move_type_differentials, &
          dM, lambda, &
          swap_dimension_low, swap_dimension_high, &
          particle_number, move_type, &
          n_slices, n_particles, n_dimensions, &
          off_diagonal_particle_number, prb_off_diagonal_particle_number, &
          n_off_diagonal_particles, &
          off_diagonal_dimension_low, off_diagonal_dimension_high &
        )
    else
      call sample_scheme1(q, q_trial, move_start, move_end, &
          move_type_probabilities, move_type_differentials, &
          dM, lambda, &
          swap_dimension_low, swap_dimension_high, &
          particle_number, move_type, &
          n_slices, n_particles, n_dimensions &
        )
    end if

    if( move_start < 1 .or. move_start > move_end .or. move_end > n_slices ) then
      call revert_move
      cycle
    end if

    !@<< Impose periodic boundary conditions >>
    !@+node:gcross.20090626112946.1691:<< Impose periodic boundary conditions >>
    if( system_is_periodic ) then
      q_trial(:,:,move_start:move_end) = wrap_around(x=q_trial(:,:,move_start:move_end),period_length=pbc_period_length)
    end if
    !@-node:gcross.20090626112946.1691:<< Impose periodic boundary conditions >>
    !@nl

    !@<< Update move type statistics >>
    !@+node:gcross.20090626112946.1690:<< Update move type statistics >>
    slice_move_attempted_counts(move_start:move_end) =  slice_move_attempted_counts(move_start:move_end) + 1
    move_type_attempted_counts(move_type) = move_type_attempted_counts(move_type) + 1
    !@-node:gcross.20090626112946.1690:<< Update move type statistics >>
    !@nl

    !@<< Update the displacement matrix >>
    !@+node:gcross.20090626112946.1689:<< Update the displacement matrix >>
    if( system_is_periodic ) then
      call update_xij_pbc( &
        q_trial(:,:,move_start:move_end), particle_number, &
        pbc_period_length, (move_end-move_start+1), n_particles, n_dimensions, &
        xij2_trial(:,:,move_start:move_end) &
        )
    else
      call update_xij( &
        q_trial(:,:,move_start:move_end), particle_number, &
        (move_end-move_start+1), n_particles, n_dimensions, &
        xij2_trial(:,:,move_start:move_end) &
        )
    end if
    !@nonl
    !@-node:gcross.20090626112946.1689:<< Update the displacement matrix >>
    !@nl
    !@nonl
    !@-node:gcross.20090626112946.1688:<< Randomly choose a move to make and apply it to the system >>
    !@nl

    !@    << Determine whether the move should be accepted >>
    !@+node:gcross.20090626112946.1693:<< Determine whether the move should be accepted >>
    !@<< Compute logarithmic probability of acceptance >>
    !@+node:gcross.20090721121051.1755:<< Compute logarithmic probability of acceptance >>
    reject_flag = .false.

    old_weight = compute_log_acceptance_weight (&
        q, xij2, &
        move_start, move_end, &
        U, gradU2, &
        reject_flag &
        )

    if(reject_flag) then
      stop "Error:  A path that had been okay before is somehow invalid now."
    end if

    reject_flag = .false.

    new_weight = compute_log_acceptance_weight (&
        q_trial, xij2_trial, &
        move_start, move_end, &
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
    if ( accept_path( old_weight, new_weight ) ) then
      !@  << Update the degrees of freedom >>
      !@+node:gcross.20090626112946.1701:<< Update the degrees of freedom >>
      q(:,:,move_start:move_end) = q_trial(:,:,move_start:move_end)
      !@-node:gcross.20090626112946.1701:<< Update the degrees of freedom >>
      !@nl
      !@  << Update derived quantities >>
      !@+node:gcross.20090805093617.1831:<< Update derived quantities >>
      xij2(:,:,move_start:move_end) = xij2_trial(:,:,move_start:move_end) 
      U(:,move_start:move_end) = U_trial(:,move_start:move_end)
      gradU2(move_start:move_end) = gradU2_trial(move_start:move_end)
      !@-node:gcross.20090805093617.1831:<< Update derived quantities >>
      !@nl
      !@  << Update statistics >>
      !@+node:gcross.20090626112946.1703:<< Update statistics >>
      slice_move_accepted_counts(move_start:move_end) =  slice_move_accepted_counts(move_start:move_end) + 1
      move_type_accepted_counts(move_type) = move_type_accepted_counts(move_type) + 1
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

  contains

  !@  @+others
  !@+node:gcross.20090817102318.2271:Interface
  !@-node:gcross.20090817102318.2271:Interface
  !@+node:gcross.20090817102318.2263:Subroutines
  !@+node:gcross.20090805093617.1835:subroutine revert_move
  subroutine revert_move
    q_trial(:,:,move_start:move_end) = q(:,:,move_start:move_end)
    xij2_trial(:,:,move_start:move_end) = xij2(:,:,move_start:move_end)
    U_trial(:,move_start:move_end) = U(:,move_start:move_end)
    gradU2_trial(move_start:move_end) =  gradU2(move_start:move_end)
  end subroutine revert_move
  !@-node:gcross.20090805093617.1835:subroutine revert_move
  !@+node:gcross.20090817102318.2264:subroutine compute_log_acceptance_weight
  function compute_log_acceptance_weight (&
      x, xij2, &
      move_start, move_end, &
      U, gradU2, &
      reject_flag &
      ) result ( weight )

  !@+at
  ! Array dimensions / slicing
  !@-at
  !@@c
    integer, intent(in) :: move_start, move_end

  !@+at
  ! Function input
  !@-at
  !@@c

    double precision, dimension (n_dimensions,n_particles,n_slices), intent(in) :: x
    double precision, dimension (n_particles,n_particles,n_slices), intent(in) :: xij2

  !@+at
  ! Function output
  !@-at
  !@@c
    double precision, dimension(n_particles,n_slices), intent(out) :: U
    double precision, dimension(n_slices), intent(out) :: gradU2
    logical, intent(out) :: reject_flag
    double precision :: weight

  !@+at
  ! Local variables.
  !@-at
  !@@c
    double precision :: lngfn, lntfn
    integer :: slice_start, slice_end

  !@+at
  ! Code begins:
  !@-at
  !@@c

    !@  << Compute contribution from potentials >>
    !@+node:gcross.20090817102318.2266:<< Compute contribution from potentials >>
    call compute_potential (&
      x(:,:,move_start:move_end), xij2(:,:,move_start:move_end), &
      (move_end-move_start+1), n_particles, n_dimensions, &
      U(:,move_start:move_end), gradU2(move_start:move_end), &
      reject_flag &
      )

    if (reject_flag) then
      return
    end if

    if(move_start .le. 1) then
      slice_start = 1
    else
      slice_start = move_start-1
    end if
    if(move_end .ge. n_slices) then
      slice_end = n_slices
    else
      slice_end = move_end+1
    end if

    if( n_slices > 2 ) then
      lngfn = greens_function( &
                x, xij2, U, gradU2, &
                lambda, dtau, &
                slice_start, slice_end, &
                particle_number, &
                n_slices, n_particles, n_dimensions &
              )
    else
      lngfn = 0
    end if
    !@nonl
    !@-node:gcross.20090817102318.2266:<< Compute contribution from potentials >>
    !@nl

    !@  << Compute contribution from trial functions >>
    !@+node:gcross.20090817102318.2270:<< Compute contribution from trial functions >>
    lntfn = 0d0

    lntfn = lntfn + trial_function( &
                      x(:,:,1), xij2(:,:,1), &
                      n_particles, n_dimensions, &
                      reject_flag &
                    )
    if(reject_flag) then
      return
    end if

    lntfn = lntfn + trial_function( &
                      x(:,:,n_slices), xij2(:,:,n_slices), &
                      n_particles, n_dimensions, &
                      reject_flag &
                    )
    if(reject_flag) then
      return
    end if
    !@nonl
    !@-node:gcross.20090817102318.2270:<< Compute contribution from trial functions >>
    !@nl

    weight = lntfn + lngfn

  end function compute_log_acceptance_weight
  !@nonl
  !@-node:gcross.20090817102318.2264:subroutine compute_log_acceptance_weight
  !@-node:gcross.20090817102318.2263:Subroutines
  !@-others

end subroutine thermalize_path
!@-node:gcross.20090626112946.1694:thermalize_path
!@-others

end module thermalize
!@-node:gcross.20090805153643.1843:@thin thermalize.f95
!@-leo
