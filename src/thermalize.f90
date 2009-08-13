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
! Potential functions
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
  double precision, dimension( n_slices, n_particles ), target, intent(out) :: U
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
!@-others

end module thermalize
!@-node:gcross.20090805153643.1843:@thin thermalize.f90
!@-leo
