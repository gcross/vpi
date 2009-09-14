!@+leo-ver=4-thin
!@+node:gcross.20090903090230.2039:@thin numeric_differentiation.f95
!@@language fortran90

module numeric_differentiation

  use xij

  implicit none

contains 

!@+others
!@+node:gcross.20090903090230.2040:numerically_differentiate_jastrow
subroutine numerically_differentiate_jastrow( &
    x, xij2, &
    n_particles, n_dimensions, &
    jfunc, &
    grad_lntfn, lap_lntfn &
  )
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension( n_particles, n_dimensions ), intent(in) :: x
  double precision, dimension( n_particles, n_particles ), intent(in) :: xij2
  double precision, dimension( n_particles, n_dimensions ), intent(out) :: grad_lntfn
  double precision, intent(out) :: lap_lntfn

!f2py external, intent(callback) :: jfunc
  interface
    pure function jfunc( x, xij2, n_particles, n_dimensions ) result(y)
      integer, intent(in) :: n_particles, n_dimensions
      double precision, dimension( n_particles, n_dimensions ), intent(in) :: x
      double precision, dimension( n_particles, n_particles ), intent(in) :: xij2
      double precision  :: y
    end function
  end interface

  double precision, parameter ::  epsilon = 1d-2

  double precision, dimension( n_particles, n_dimensions ) :: tx
  double precision, dimension( n_particles, n_particles ) :: txij2
  double precision, dimension( n_dimensions ) :: fhi,flo
  double precision :: f0
  integer :: i, k
  double precision :: tmp

  grad_lntfn = 0.0d0
  lap_lntfn = 0.0d0
  f0 = jfunc( x, xij2, n_particles, n_dimensions )
  do i = 1, n_particles
    do k = 1, n_dimensions
      tx(:,:) = x(:,:)
      tx(i,k) = x(i,k) + epsilon
      txij2(:,:) = xij2(:,:)
      call update_xij( txij2, tx, 1, n_particles, n_dimensions )
      fhi(k) = jfunc( tx, txij2, n_particles, n_dimensions )

      txij2(:,:) = xij2(:,:)
      tx(i,k) = x(i,k) - epsilon
      call update_xij( txij2, tx, 1, n_particles, n_dimensions )
      flo(k) = jfunc( tx, txij2, n_particles, n_dimensions )
    end do
    grad_lntfn(i,:) = fhi(:) - flo(:)
    tmp = sum( fhi(:) + flo(:) ) -2.0d0*n_dimensions*f0
    lap_lntfn = lap_lntfn + tmp
  end do
  grad_lntfn(:,:) = grad_lntfn(:,:)/(2.0d0*epsilon)
  lap_lntfn = lap_lntfn / epsilon**2

end subroutine
!@-node:gcross.20090903090230.2040:numerically_differentiate_jastrow
!@-others

end module
!@-node:gcross.20090903090230.2039:@thin numeric_differentiation.f95
!@-leo
