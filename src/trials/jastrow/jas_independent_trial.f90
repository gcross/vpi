!@+leo-ver=4-thin
!@+node:gcross.20090624144408.2038:@thin jas_independent_trial.f90
!@@language fortran90
module jas_independent_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.2039:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.2039:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.2040:<< Variables >>
  !@-node:gcross.20090624144408.2040:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.2041:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.2042:init_jas_tfunc
  subroutine init_jas_tfunc ()
    write(*,*) "Using independent Jastrow trial function."
  end subroutine init_jas_tfunc
  !@-node:gcross.20090624144408.2042:init_jas_tfunc
  !@+node:gcross.20090624144408.2043:jas_tfunc
  function jas_tfun( x, xij2, sl, nslice, np, ndim ) result( y )
    use kinds
    implicit none
    integer :: sl, nslice, np, ndim
    real(kind=b8), dimension( nslice , np, ndim ), intent(in) :: x
    real(kind=b8), dimension( nslice , np, np ), intent(in) :: xij2
    real(kind=b8)  :: r2
    real(kind=b8)  :: y

    y = 0

  end function jas_tfun
  !@-node:gcross.20090624144408.2043:jas_tfunc
  !@+node:gcross.20090624144408.2044:grad & lapacian of tfunc
  function grad_lap_jas_tfun( x, xij2, sl, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
    integer :: sl, np, ndim, nslice
    real(kind=b8), dimension( nslice, np, ndim ) :: x
    real(kind=b8), dimension( nslice, np, np ) :: xij2
    real(kind=b8), dimension( np, ndim ), intent(out) :: grad_lntfn 
    real(kind=b8), intent(out) :: lap_lntfn 
    integer :: y


    real(kind=b8)  :: fi,fi2,ri,ri2,ri3,ri4
    real(kind=b8), dimension( ndim ) :: gtmp
    integer :: i, j

    grad_lntfn = 0.0_b8
    lap_lntfn = 0.0_b8
    y = 1
  end function grad_lap_jas_tfun

  !@-node:gcross.20090624144408.2044:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.2041:<< Subroutines >>
  !@nl

end module jas_independent_trial
!@-node:gcross.20090624144408.2038:@thin jas_independent_trial.f90
!@-leo
