!@+leo-ver=4-thin
!@+node:gcross.20090624144408.2010:@thin jas_lj_trial.f90
!@@language fortran90
module jas_lj_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.2011:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.2011:<< Imported modules >>
  !@nl

  implicit none

  !@  << Variables >>
  !@+node:gcross.20090624144408.2012:<< Variables >>
  real(kind=b8) :: p_ljc5 = 1.14083e-07_b8
  real(kind=b8) :: p_ljc1 = 0.0123049_b8
  !@-node:gcross.20090624144408.2012:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.2013:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.2014:init_jas_tfunc
  subroutine init_jas_tfunc ()
    namelist /jastrow_trial_parameters/ p_ljc1, p_ljc5

    read(unit=10,nml=jastrow_trial_parameters)

    write(*,*) "Using L.J. Jastrow trial function with"
    write(*,nml=jastrow_trial_parameters)
  end subroutine init_jas_tfunc
  !@-node:gcross.20090624144408.2014:init_jas_tfunc
  !@+node:gcross.20090624144408.2015:tfunc
  function jas_tfun( x, xij2, p, sl, nslice, np, ndim ) result( y )
    integer :: sl, nslice, np, ndim
    real(kind=b8) , dimension( nslice , np, ndim ), intent(in) :: x
    real(kind=b8) , dimension( nslice , np, np ), intent(in) :: xij2
    real(kind=b8) , dimension(:) :: p
    real(kind=b8)  :: y
    real(kind=b8)  :: r5, r1

    integer :: i, j

    y = 0.0
    do i = 1, np
      do j = i + 1, np
        r1 = xij2(sl,i,j)**(-0.5_b8)
        r5 = r1*xij2(sl,i,j)**(-2)
        y = y - p_ljc5*r5 - p_ljc1*r1
      end do
    end do

  end function jas_tfun
  !@-node:gcross.20090624144408.2015:tfunc
  !@+node:gcross.20090624144408.2016:grad & lapacian of tfunc
  function grad_lap_jas_tfun( x, xij2, sl, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
    integer, intent(in) :: sl, np, ndim, nslice
    real(kind=b8), dimension( nslice, np, ndim ) :: x
    real(kind=b8), dimension( nslice, np, np ) :: xij2
    real(kind=b8), dimension( np, ndim ), intent(out) :: grad_lntfn
    real(kind=b8), intent(out) :: lap_lntfn

    real(kind=b8)  :: r,ri,ri2,ri3,ri4,ri6,ri7,ri8,ri12
    real(kind=b8), dimension( ndim ) :: gtmp
    real(kind=b8) :: p5sq,p1sq,tmp
    integer :: y
    integer :: i, j

    p5sq = p_ljc5*p_ljc5
    p1sq = p_ljc1*p_ljc1

    grad_lntfn = 0
    lap_lntfn = 0
    do i = 1, np
      do j = i + 1, np
        r = sqrt(xij2(sl,i,j))
        ri2 = 1.0_b8/xij2(sl,i,j)
        ri = 1.0/r
        ri3 = ri*ri2
        ri4 = ri2*ri2
        ri6 = ri3*ri3
        ri8 = ri4*ri4
        ri7 = ri*ri6
        ri12 = ri6*ri6

        gtmp(:) =  (5.0*p_ljc5*ri7 + p_ljc1*ri3)*(x(sl,i,:) - x(sl,j,:))
        grad_lntfn(i,:) = grad_lntfn(i,:) + gtmp(:)
        grad_lntfn(j,:) = grad_lntfn(j,:) - gtmp(:)
        tmp = 25.0*p5sq*ri12 + 10.0*p_ljc5*(p_ljc1 - 2*r)*ri8 + p1sq*ri4 
        lap_lntfn = lap_lntfn + tmp
      end do
    end do
    y = 1
  end function grad_lap_jas_tfun
  !@-node:gcross.20090624144408.2016:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.2013:<< Subroutines >>
  !@nl

end module jas_lj_trial
!@-node:gcross.20090624144408.2010:@thin jas_lj_trial.f90
!@-leo
