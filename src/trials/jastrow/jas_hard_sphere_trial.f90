!@+leo-ver=4-thin
!@+node:gcross.20090624144408.1996:@thin jas_hard_sphere_trial.f90
!@@language fortran90
module jas_hard_sphere_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.1997:<< Imported modules >>
  use kinds
  use constants
  !@-node:gcross.20090624144408.1997:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.1998:<< Variables >>
  real (kind=b8), private :: radius
  real (kind=b8), private :: radius_squared
  !@-node:gcross.20090624144408.1998:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.1999:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.2000:init_jas_tfunc
  subroutine init_jas_tfunc ()
    namelist /jastrow_trial_parameters/ radius

    read(unit=10,nml=jastrow_trial_parameters)

    radius_squared = radius**2

    write(*,*) "Using hard-sphere Jastrow trial function with"
    write(*,nml=jastrow_trial_parameters)
  end subroutine init_jas_tfunc
  !@-node:gcross.20090624144408.2000:init_jas_tfunc
  !@+node:gcross.20090624144408.2001:jas_tfun
  function jas_tfun( x, xij2, sl, nslice, np, ndim ) result( y )
    use kinds
    implicit none
    integer :: sl, nslice, np, ndim
    real(kind=b8), dimension( nslice , np, ndim ), intent(in) :: x
    real(kind=b8), dimension( nslice , np, np ), intent(in) :: xij2
    real(kind=b8)  :: r2
    real(kind=b8)  :: y

    integer :: i, j

    y = 0.0_b8
    do i = 1, np
      do j = i + 1, np
        r2 = xij2(sl,i,j)
        if ( r2 .gt. radius_squared ) then
          y = y + log(1.0_b8 - radius/sqrt(r2))
        else
          y = -realbignumber
        end if
      end do
    end do

  end function jas_tfun
  !@nonl
  !@-node:gcross.20090624144408.2001:jas_tfun
  !@+node:gcross.20090624144408.2002:grad & lapacian of tfunc
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
    do i = 1, np
      do j = i + 1, np
        if ( xij2(sl,i,j) .gt. radius_squared ) then
          ri2 = 1.0_b8/xij2(sl,i,j)
          ri = sqrt(ri2)
          ri3 = ri*ri2
          ri4 = ri2*ri2
          fi = 1.0_b8/(1.0_b8 - radius*ri)
          fi2 = fi*fi
          gtmp(:) = ri3*fi*(x(sl,i,:) - x(sl,j,:))
  !        write (222,"(3g18.6)") gtmp(1),gtmp(2),gtmp(3)
          grad_lntfn(i,:) = grad_lntfn(i,:) + gtmp(:)
          grad_lntfn(j,:) = grad_lntfn(j,:) - gtmp(:)
          lap_lntfn = lap_lntfn - fi2*ri4 
        else
          y = -1
          return
        end if
      end do
    end do
    lap_lntfn = 2.0_b8*radius_squared*lap_lntfn
    grad_lntfn(:,:) = radius*grad_lntfn(:,:)
    y = 1
  end function grad_lap_jas_tfun

  !@-node:gcross.20090624144408.2002:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.1999:<< Subroutines >>
  !@nl

end module jas_hard_sphere_trial
!@-node:gcross.20090624144408.1996:@thin jas_hard_sphere_trial.f90
!@-leo
