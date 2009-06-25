!@+leo-ver=4-thin
!@+node:gcross.20090624144408.2024:@thin jas_charge_trial.f90
!@@language fortran90
module jas_charge_trial

  !@  << Imported modules >>
  !@+node:gcross.20090624144408.2025:<< Imported modules >>
  use kinds
  use constants
  use jas_trial_numeric_differentiator
  !@-node:gcross.20090624144408.2025:<< Imported modules >>
  !@nl

  !@  << Variables >>
  !@+node:gcross.20090624144408.2026:<< Variables >>
  real(kind=b8) :: coefficient_0
  real(kind=b8) :: coefficient_1
  real(kind=b8) :: coefficient_2
  !@-node:gcross.20090624144408.2026:<< Variables >>
  !@nl

contains

  !@  << Subroutines >>
  !@+node:gcross.20090624144408.2027:<< Subroutines >>
  !@+others
  !@+node:gcross.20090624144408.2028:init_jas_tfunc
  subroutine init_jas_tfunc ()
    namelist /jastrow_trial_parameters/ coefficient_0, coefficient_1, coefficient_2

    read(unit=10,nml=jastrow_trial_parameters)

    radius_squared = radius**2

    write(*,*) "Using charge Jastrow trial function with"
    write(*,nml=jastrow_trial_parameters)
  end subroutine init_jas_tfunc
  !@-node:gcross.20090624144408.2028:init_jas_tfunc
  !@+node:gcross.20090624144408.2029:jas_tfunc
  function jas_tfun( x, xij2, sl, nslice, np, ndim ) result( y )
    integer :: sl, nslice, np, ndim
    real(kind=b8) , dimension( nslice , np, ndim ), intent(in) :: x
    real(kind=b8) , dimension( nslice , np, np ), intent(in) :: xij2
    real(kind=b8)  :: y

    real(kind=b8)  :: r2,r1

    integer :: i, j

    y = 0.0
    do i = 1, np
      do j = i + 1, np
        r2 = xij2(sl,i,j)
        r1 = sqrt(r2)
        y = y - ( coefficient_0*r1 + coefficient_1*r2 )/( 1 + coefficient_2*r1 )
      end do
    end do

  end function jas_tfun

  !@-node:gcross.20090624144408.2029:jas_tfunc
  !@+node:gcross.20090624144408.2030:grad & lapacian of tfunc
  function grad_lap_jas_tfun( x, xij2, sl, np, ndim, nslice, grad_lntfn, lap_lntfn ) result (y)
    integer, intent(in) :: sl, np, ndim, nslice
    real(kind=b8), dimension( nslice, np, ndim ) :: x
    real(kind=b8), dimension( nslice, np, np ) :: xij2
    real(kind=b8), dimension( np, ndim ), intent(out) :: grad_lntfn
    real(kind=b8), intent(out) :: lap_lntfn

    y = numeric_grad_lap_jas( x, xij2, sl, np, ndim, nslice, grad_lntfn, lap_lntfn, jas_tfun )

  end function grad_lap_jas_tfun
  !@-node:gcross.20090624144408.2030:grad & lapacian of tfunc
  !@-others
  !@-node:gcross.20090624144408.2027:<< Subroutines >>
  !@nl

end module jas_charge_trial
!@-node:gcross.20090624144408.2024:@thin jas_charge_trial.f90
!@-leo
