!@+leo-ver=4-thin
!@+node:gcross.20090805153643.1843:@thin thermalize.f90
!@@language fortran90
module thermalize

use xij
use sample

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
!@-others

end module thermalize
!@-node:gcross.20090805153643.1843:@thin thermalize.f90
!@-leo
