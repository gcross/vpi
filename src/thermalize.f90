!@+leo-ver=4-thin
!@+node:gcross.20090805153643.1843:@thin thermalize.f90
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
!@+node:gcross.20090721121051.1764:compute_physical_potential
subroutine compute_physical_potential (&
    x, xij2, &
    Usp_func, gUsp_func, &
    Uij_func, gUij_func, &
    move_start, move_end, &
    n_slices, n_particles, n_dimensions, &
    U, gradU2, &
    reject_flag &
    )
!@+at
! Array dimensions / slicing
!@-at
!@@c
  integer, intent(in) :: move_start, move_end, n_slices, n_particles, n_dimensions
!@+at
! Function input
!@-at
!@@c
  double precision, dimension ( n_slices, n_particles, n_dimensions ), intent(in) :: x
  double precision, dimension ( n_slices, n_particles, n_particles ), intent(in) :: xij2
!@+at
! Potential functions
!@-at
!@@c
!f2py external, intent(callback) :: Usp_func, gUsp_func, Uij_func, gUij_func
interface
  !@  << Potential callback interface >>
  !@+node:gcross.20090809223137.1724:<< Potential callback interface >>
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer, intent(in) :: nslice, np, ndim
    integer, intent(in) :: slice, ip
    double precision, dimension ( nslice, np , ndim ), intent(in) :: x
    double precision :: Usp
  end function Usp_func

  subroutine gUsp_func( x, slice, nslice, np, ndim, gUsp )
    double precision, dimension ( nslice, np , ndim ), intent(in) ::  x
    integer, intent(in) :: slice, nslice, np, ndim
    double precision, dimension ( np, ndim ), intent(out) :: gUsp
  end subroutine gUsp_func

  function Uij_func( x, xij2, slice, ip, nslice, np, ndim, reject_flag ) result ( Uij )
    double precision, dimension ( nslice, np , ndim ), intent(in) :: x
    double precision, dimension ( nslice, np , np ), intent(in) :: xij2
    integer, intent(in) :: slice, ip, nslice, np, ndim
    logical, intent(out) :: reject_flag
    double precision :: Uij
  end function Uij_func

  subroutine gUij_func( x, xij2, slice, nslice, np, ndim, gUij )
    double precision, dimension ( nslice, np , ndim ), intent(in) :: x
    double precision, dimension ( nslice, np , np ), intent(in) :: xij2
    integer, intent(in) :: slice, nslice, np, ndim
    double precision, dimension ( np , ndim ), intent(out) :: gUij
  end subroutine gUij_func
  !@-node:gcross.20090809223137.1724:<< Potential callback interface >>
  !@nl
end interface
!@+at
! Function output
!@-at
!@@c
  double precision, dimension( n_slices, n_particles ), intent(out) :: U
  double precision, dimension( n_slices ), intent(out) :: gradU2
  logical, intent(out) :: reject_flag
!@+at
! Local variables.
!@-at
!@@c
  double precision, dimension( n_particles, n_dimensions ) :: grad_Usp, grad_Uij, grad_Uij_rot, grad_U
  double precision :: U_sp, U_ij
  integer :: ii, jj
!@+at
! Code begins:
!@-at
!@@c

  reject_flag = .false.

  do ii = move_start, move_end
    U(ii,:) = 0.0d0
    do jj = 1, n_particles
      U_sp = Usp_func( x, ii, jj, n_slices, n_particles, n_dimensions )
      U_ij = Uij_func( x, xij2, ii, jj, n_slices, n_particles, n_dimensions, reject_flag )
      if(reject_flag) then
          return
      end if
      U(ii,jj) = U(ii,jj) + U_sp + U_ij
    end do

    call gUsp_func( x, ii, n_slices, n_particles, n_dimensions, grad_Usp )
    call gUij_func( x, xij2, ii, n_slices, n_particles, n_dimensions, grad_Uij )
    grad_U = grad_Usp + grad_Uij
    gradU2(ii) = sum( grad_U(:,:)**2 )
  end do

end subroutine
!@-node:gcross.20090721121051.1764:compute_physical_potential
!@+node:gcross.20090721121051.1746:compute_log_acceptance_weight
subroutine compute_log_acceptance_weight (&
    x, xij2, &
    move_start, move_end, &
    n_slices, n_particles, n_dimensions, &
    particle_number, &
    Usp_func, gUsp_func, &
    Uij_func, gUij_func, &
    U, gradU2, &
    U_weights, gU2_weights, &
    fixed_rotation_axis, frame_angular_velocity, fixed_angular_momentum, &
    lambda, dtau, &
    use_4th_order_green_function, &
    sp_trial_function, jastrow_trial_function, &
    reject_flag, &
    weight &
    )

!@+at
! Array dimensions / slicing
!@-at
!@@c
  integer, intent(in) :: move_start, move_end, n_slices, n_particles, n_dimensions

!@+at
! Function input
!@-at
!@@c

  double precision, dimension ( n_slices, n_particles , n_dimensions ), intent(in) :: x
  double precision, dimension ( n_slices, n_particles , n_particles ), intent(in) :: xij2
  logical, intent(in) :: use_4th_order_green_function
  double precision, dimension ( n_slices ), intent(in) :: U_weights, gU2_weights
  integer, intent(in) :: fixed_rotation_axis, fixed_angular_momentum
  double precision, intent(in) :: frame_angular_velocity
  integer, intent(in) :: particle_number
  double precision, intent(in) :: lambda, dtau

!@+at
! Potential and trial functions
!@-at
!@@c
!f2py external, intent(callback) :: Usp_func, gUsp_func, Uij_func, gUij_func
!f2py external, intent(callback) :: sp_trial_function, jastrow_trial_function
interface
  !@  << Potential callback interface >>
  !@+node:gcross.20090809223137.1724:<< Potential callback interface >>
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer, intent(in) :: nslice, np, ndim
    integer, intent(in) :: slice, ip
    double precision, dimension ( nslice, np , ndim ), intent(in) :: x
    double precision :: Usp
  end function Usp_func

  subroutine gUsp_func( x, slice, nslice, np, ndim, gUsp )
    double precision, dimension ( nslice, np , ndim ), intent(in) ::  x
    integer, intent(in) :: slice, nslice, np, ndim
    double precision, dimension ( np, ndim ), intent(out) :: gUsp
  end subroutine gUsp_func

  function Uij_func( x, xij2, slice, ip, nslice, np, ndim, reject_flag ) result ( Uij )
    double precision, dimension ( nslice, np , ndim ), intent(in) :: x
    double precision, dimension ( nslice, np , np ), intent(in) :: xij2
    integer, intent(in) :: slice, ip, nslice, np, ndim
    logical, intent(out) :: reject_flag
    double precision :: Uij
  end function Uij_func

  subroutine gUij_func( x, xij2, slice, nslice, np, ndim, gUij )
    double precision, dimension ( nslice, np , ndim ), intent(in) :: x
    double precision, dimension ( nslice, np , np ), intent(in) :: xij2
    integer, intent(in) :: slice, nslice, np, ndim
    double precision, dimension ( np , ndim ), intent(out) :: gUij
  end subroutine gUij_func
  !@-node:gcross.20090809223137.1724:<< Potential callback interface >>
  !@nl
  !@  << Trial callback interface >>
  !@+node:gcross.20090812093015.1848:<< Trial callback interface >>
  function sp_trial_function( x, sl, nslice, np, ndim ) result( y )
    integer, intent(in) :: sl, nslice, np, ndim
    double precision, dimension( nslice, np , ndim ), intent(in) :: x
    double precision :: y
  end function sp_trial_function

  function jastrow_trial_function( x, xij2, sl, nslice, np, ndim ) result( y )
    integer :: sl, nslice, np, ndim
    double precision, dimension( nslice , np, ndim ), intent(in) :: x
    double precision, dimension( nslice , np, np ), intent(in) :: xij2
    double precision  :: y
  end function jastrow_trial_function
  !@-node:gcross.20090812093015.1848:<< Trial callback interface >>
  !@nl
end interface

!@+at
! Function output
!@-at
!@@c
  double precision, dimension( n_slices, n_particles ), intent(out) :: U
  double precision, dimension( n_slices ), intent(out) :: gradU2
  double precision, intent(out) :: weight
  logical, intent(out) :: reject_flag

!@+at
! Local variables.
!@-at
!@@c
  double precision :: lngfn, rotgfn, hsgfn, lntfn
  integer :: sl_start, sl_end
!@+at
! Code begins:
!@-at
!@@c
  !@  << Compute contribution from potentials >>
  !@+node:gcross.20090626112946.1696:<< Compute contribution from potentials >>
  call compute_physical_potential (&
    x, xij2, &
    Usp_func, gUsp_func, &
    Uij_func, gUij_func, &
    move_start, move_end, &
    n_slices, n_particles, n_dimensions, &
    U, gradU2, &
    reject_flag &
    )

  if (reject_flag) then
    return
  endif

  if (fixed_angular_momentum > 0) then
    call compute_effective_rotational_potential (&
      x, &
      fixed_rotation_axis, frame_angular_velocity, fixed_angular_momentum, &
      move_start, move_end, &
      n_slices, n_particles, n_dimensions, &
      U &
      )
  end if

  !@-node:gcross.20090626112946.1696:<< Compute contribution from potentials >>
  !@nl

  !@  << Compute contribution from propagators >>
  !@+node:gcross.20090721121051.1752:<< Compute contribution from propagators >>
  if(move_start .le. 1) then
    sl_start = 1
  else
    sl_start = move_start-1
  end if
  if(move_end .ge. n_slices) then
    sl_end = n_slices
  else
    sl_end = move_end+1
  end if

  if( n_slices > 2 ) then 
    if ( use_4th_order_green_function ) then
      lngfn = gfn4_sp( sl_start, sl_end, particle_number, U, gradU2, U_weights, gU2_weights, &
                           n_slices, n_particles, lambda, dtau ) 
    else 
      lngfn = gfn2_sp( sl_start, sl_end, particle_number, U, n_slices, n_particles, dtau ) 
    end if
  end if

  !@-node:gcross.20090721121051.1752:<< Compute contribution from propagators >>
  !@nl

  !@  << Compute contribution from trial functions >>
  !@+node:gcross.20090721121051.1754:<< Compute contribution from trial functions >>
  lntfn = sp_trial_function(x, 1, n_slices, n_particles, n_dimensions) + &
          sp_trial_function(x, n_slices, n_slices, n_particles, n_dimensions) + &
          jastrow_trial_function(x, xij2, 1, n_slices, n_particles, n_dimensions) + &
          jastrow_trial_function(x, xij2, n_slices, n_slices, n_particles, n_dimensions)
  !@-node:gcross.20090721121051.1754:<< Compute contribution from trial functions >>
  !@nl

  weight = lngfn + lntfn

end subroutine

!@-node:gcross.20090721121051.1746:compute_log_acceptance_weight
!@+node:gcross.20090626112946.1694:thermalize_path
subroutine thermalize_path( &
  q, xij2, U, gradU2, &
  n_slices, n_particles, n_dimensions, &
  n_trials, &
  move_type_probabilities, move_type_differentials, &
  dM, lambda, &
  low_swap_dim, high_swap_dim, &
  slice_move_attempted_counts, move_type_attempted_counts, &
  slice_move_accepted_counts, move_type_accepted_counts, &
  Usp_func, gUsp_func, &
  Uij_func, gUij_func, &
  U_weights, gU2_weights, &
  fixed_rotation_axis, frame_angular_velocity, fixed_angular_momentum, &
  use_4th_order_green_function, &
  sp_trial_function, jastrow_trial_function, &
  pbc_period_length, &
  od_pnum, PROB_OD_PNUM, &
  n_od_particle, &
  od_dim_low, od_dim_high &
  )

!@+at
! Input variables
!@-at
!@@c
  integer, intent(in) :: n_slices, n_particles, n_dimensions, n_trials
  double precision, dimension( N_MOVE_TYPES ), intent(in) :: move_type_probabilities, move_type_differentials
  double precision, intent(in) :: dM, lambda
  integer, intent(in) :: low_swap_dim, high_swap_dim
  double precision, intent(in), optional :: pbc_period_length
  logical, intent(in) :: use_4th_order_green_function
  double precision, dimension ( n_slices ), intent(in) :: U_weights, gU2_weights
  integer, intent(in) :: fixed_rotation_axis, fixed_angular_momentum
  double precision, intent(in) :: frame_angular_velocity
  double precision, intent(in), optional :: PROB_OD_PNUM
  integer, intent(in), optional :: od_pnum, n_od_particle, od_dim_low, od_dim_high

!@+at
! Potential and trial functions
!@-at
!@@c
!f2py external, intent(callback) :: Usp_func, gUsp_func, Uij_func, gUij_func
!f2py external, intent(callback) :: sp_trial_function, jastrow_trial_function
interface
  !@  << Potential callback interface >>
  !@+node:gcross.20090809223137.1724:<< Potential callback interface >>
  function Usp_func( x, slice, ip, nslice, np, ndim ) result ( Usp )
    integer, intent(in) :: nslice, np, ndim
    integer, intent(in) :: slice, ip
    double precision, dimension ( nslice, np , ndim ), intent(in) :: x
    double precision :: Usp
  end function Usp_func

  subroutine gUsp_func( x, slice, nslice, np, ndim, gUsp )
    double precision, dimension ( nslice, np , ndim ), intent(in) ::  x
    integer, intent(in) :: slice, nslice, np, ndim
    double precision, dimension ( np, ndim ), intent(out) :: gUsp
  end subroutine gUsp_func

  function Uij_func( x, xij2, slice, ip, nslice, np, ndim, reject_flag ) result ( Uij )
    double precision, dimension ( nslice, np , ndim ), intent(in) :: x
    double precision, dimension ( nslice, np , np ), intent(in) :: xij2
    integer, intent(in) :: slice, ip, nslice, np, ndim
    logical, intent(out) :: reject_flag
    double precision :: Uij
  end function Uij_func

  subroutine gUij_func( x, xij2, slice, nslice, np, ndim, gUij )
    double precision, dimension ( nslice, np , ndim ), intent(in) :: x
    double precision, dimension ( nslice, np , np ), intent(in) :: xij2
    integer, intent(in) :: slice, nslice, np, ndim
    double precision, dimension ( np , ndim ), intent(out) :: gUij
  end subroutine gUij_func
  !@-node:gcross.20090809223137.1724:<< Potential callback interface >>
  !@nl
  !@  << Trial callback interface >>
  !@+node:gcross.20090812093015.1848:<< Trial callback interface >>
  function sp_trial_function( x, sl, nslice, np, ndim ) result( y )
    integer, intent(in) :: sl, nslice, np, ndim
    double precision, dimension( nslice, np , ndim ), intent(in) :: x
    double precision :: y
  end function sp_trial_function

  function jastrow_trial_function( x, xij2, sl, nslice, np, ndim ) result( y )
    integer :: sl, nslice, np, ndim
    double precision, dimension( nslice , np, ndim ), intent(in) :: x
    double precision, dimension( nslice , np, np ), intent(in) :: xij2
    double precision  :: y
  end function jastrow_trial_function
  !@-node:gcross.20090812093015.1848:<< Trial callback interface >>
  !@nl
end interface

!@+at
! Updated trial functions
!@-at
!@@c
  double precision, dimension( n_slices, n_particles, n_dimensions ), intent(inout) :: q
  double precision, dimension( n_slices, n_particles, n_particles ), intent(inout) :: xij2
  double precision, dimension( n_slices, n_particles ), intent(inout) :: U
  double precision, dimension( n_slices ), intent(inout) :: gradU2
  integer, dimension( n_slices ), intent(inout) :: slice_move_attempted_counts, slice_move_accepted_counts
  integer, dimension( N_MOVE_TYPES ), intent(inout) :: move_type_attempted_counts, move_type_accepted_counts

!@+at
! Local variables
!@-at
!@@c
  integer :: i, move_type, move_start, move_end, particle_number
  double precision :: old_weight, new_weight
  logical :: reject_flag
  double precision, dimension( N_SLICES, N_PARTICLES, N_DIMENSIONS ) :: q_trial
  double precision, dimension( N_SLICES, N_PARTICLES, N_PARTICLES ) :: xij2_trial
  double precision, dimension( N_SLICES, N_PARTICLES ) :: U_trial
  double precision, dimension( N_SLICES ) :: gradU2_trial

  !@  << Initialize trial functions >>
  !@+node:gcross.20090805153643.1848:<< Initialize trial functions >>
  q_trial = q
  xij2_trial = xij2
  U_trial = U
  gradU2_trial = gradU2
  !@nonl
  !@-node:gcross.20090805153643.1848:<< Initialize trial functions >>
  !@nl

  do i = 1, n_trials
    !@    << Randomly choose a move to make and apply it to the system >>
    !@+node:gcross.20090626112946.1688:<< Randomly choose a move to make and apply it to the system >>
    if( present(od_pnum) ) then
      call sample_scheme1(q, q_trial, move_start, move_end, &
          move_type_probabilities, move_type_differentials, &
          dM, lambda, &
          low_swap_dim, high_swap_dim, &
          particle_number, move_type, &
          n_slices, n_particles, n_dimensions, &
          od_pnum, PROB_OD_PNUM, &
          n_od_particle, &
          od_dim_low, od_dim_high &
        )
    else
      call sample_scheme1(q, q_trial, move_start, move_end, &
          move_type_probabilities, move_type_differentials, &
          dM, lambda, &
          low_swap_dim, high_swap_dim, &
          particle_number, move_type, &
          n_slices, n_particles, n_dimensions &
        )
    end if
    !@<< Impose periodic boundary conditions >>
    !@+node:gcross.20090626112946.1691:<< Impose periodic boundary conditions >>
    if( present(pbc_period_length) ) then
      q_trial(move_start:move_end,:,:) = wrap_around(x=q_trial(move_start:move_end,:,:),period_length=pbc_period_length)
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
    if( present(pbc_period_length) ) then
      call update_xij_pbc( xij2_trial, q_trial, pbc_period_length, move_start, move_end, N_SLICES, N_PARTICLES, N_DIMENSIONS  )
    else
      call update_xij( xij2_trial, q_trial, move_start, move_end, N_SLICES, N_PARTICLES, N_DIMENSIONS  )
    end if
    !@-node:gcross.20090626112946.1689:<< Update the displacement matrix >>
    !@nl
    !@-node:gcross.20090626112946.1688:<< Randomly choose a move to make and apply it to the system >>
    !@nl

    !@    << Determine whether the move should be accepted >>
    !@+node:gcross.20090626112946.1693:<< Determine whether the move should be accepted >>
    !@<< Compute logarithmic probability of acceptance >>
    !@+node:gcross.20090721121051.1755:<< Compute logarithmic probability of acceptance >>
    call compute_log_acceptance_weight (&
        q, xij2, &
        move_start, move_end, &
        n_slices, n_particles, n_dimensions, &
        particle_number, &
        Usp_func, gUsp_func, &
        Uij_func, gUij_func, &
        U, gradU2, &
        U_weights, gU2_weights, &
        fixed_rotation_axis, frame_angular_velocity, fixed_angular_momentum, &
        lambda, move_type_differentials(MT_BBRIDGE), &
        use_4th_order_green_function, &
        sp_trial_function, jastrow_trial_function, &
        reject_flag, &
        old_weight &
        )

    if(reject_flag) then
      stop "Error:  A path that had been okay before is somehow invalid now."
    end if

    call compute_log_acceptance_weight (&
        q_trial, xij2_trial, &
        move_start, move_end, &
        n_slices, n_particles, n_dimensions, &
        particle_number, &
        Usp_func, gUsp_func, &
        Uij_func, gUij_func, &
        U_trial, gradU2_trial, &
        U_weights, gU2_weights, &
        fixed_rotation_axis, frame_angular_velocity, fixed_angular_momentum, &
        lambda, move_type_differentials(MT_BBRIDGE), &
        use_4th_order_green_function, &
        sp_trial_function, jastrow_trial_function, &
        reject_flag, &
        new_weight &
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
      q(move_start:move_end,:,:) = q_trial(move_start:move_end,:,:)
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
  !@+node:gcross.20090805093617.1835:subroutine revert_move
  subroutine revert_move
    q_trial(move_start:move_end,:,:) = q(move_start:move_end,:,:)
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
