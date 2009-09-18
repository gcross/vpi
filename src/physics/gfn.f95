!@+leo-ver=4-thin
!@+node:gcross.20090812093015.1722:@thin gfn.f95
!@@language fortran90
!@@tabwidth -2

module gfn

  implicit none

contains 

!@+others
!@+node:gcross.20090812093015.1753:initialize weights
subroutine initialize_4th_order_weights(n_slices,U_weight,gU2_weight)
  integer, intent(in) :: n_slices
  double precision, dimension(n_slices), intent(out) :: U_weight, gU2_weight
  integer :: cslice, ii

  cslice = n_slices / 2

  if (mod(CSLICE,2) .ne. 1) then
    print *,"ERROR: CSLICE = N_SLICE/2 must be odd"
    stop
  end if

!@+at
! These "slice weights" are used to give the proper weight to even
! and odd terms when calculating the Green's function.
! In the fourth order short time approximation for the Green's funciton,
! odd slices contribute twice the potential energy of even slices
! the gradient of the potential energy only contributes for odd slides.
!@-at
!@@c
  U_weight = 0.0
  gU2_weight = 0.0
!@+at
! The fourth order propagator works in sets of 3 like this
!   0.5 1 0.5 , 0.5 1 0.5
! We concatentate the adjacent half weighted steps (unless we are at the end 
! of a path)
! in order for paths to be broken properly at the center we need CSLICE to be 
! odd.
!@-at
!@@c
  do ii = 1, CSLICE
    U_weight(ii) = mod(ii+1,2) + 1 
    gU2_weight(ii) = mod(ii+1,2)
  end do
  ! the center slice is doubled to handle broken paths 
  ! (probably a less kludgy way to do this but... )
  do ii = CSLICE+2,n_slices
    U_weight(ii) = mod(ii,2) + 1 
    gU2_weight(ii) = mod(ii,2)
  end do
  U_weight(1) = 0.5d0
  U_weight(N_SLICES) = 0.5d0
  U_weight(CSLICE) = 0.5d0
  U_weight(CSLICE+1) = 0.5d0
end subroutine
!@-node:gcross.20090812093015.1753:initialize weights
!@+node:gcross.20090624144408.1799:2nd order
pure function gfn2_sp( sl_start, sl_end, ip, U, nslice, np, dt ) result ( ln_gfn )
  double precision, dimension( nslice, np ), intent(in) :: U
  integer, intent(in) :: sl_start, sl_end, ip
  integer, intent(in) :: nslice, np
  double precision, intent(in) :: dt
  double precision :: ln_gfn

  ln_gfn = -dt*sum( U(sl_start:sl_end,ip) )

end function gfn2_sp
!@-node:gcross.20090624144408.1799:2nd order
!@+node:gcross.20090812093015.1845:4th order
pure function gfn4_sp( sl_start, sl_end, ip, U, gradU2, U_weight, gU2_weight, nslice, np, lambda, dt ) result ( ln_gfn )
  double precision, dimension( nslice, np ), intent(in) :: U
  double precision, dimension( nslice ), intent(in) :: gradU2
  double precision, dimension( nslice ), intent(in) :: U_weight
  double precision, dimension( nslice ), intent(in) :: gU2_weight
  integer, intent(in) :: sl_start, sl_end, ip
  integer, intent(in) :: nslice, np
  double precision, intent(in) :: lambda, dt
  double precision :: ln_gfn

  integer :: slice_length

  interface
    pure function ddot(n,x,incx,y,incy)
      integer, intent(in) :: n, incx, incy
      double precision, intent(in), dimension(n*incx) :: x
      double precision, intent(in), dimension(n*incy) :: y
      double precision :: ddot
    end function ddot
  end interface

  slice_length = sl_end-sl_start+1

  ln_gfn = -2.0d0*dt*ddot(slice_length,U(sl_start,ip),1,U_weight(sl_start),1)/3.0d0 &
           -2.0d0*lambda*(dt**3)*ddot(slice_length,gradU2(sl_start),1,gU2_weight(sl_start),1)/9.0d0

end function gfn4_sp
!@-node:gcross.20090812093015.1845:4th order
!@+node:gcross.20090828095451.1453:hard wall
pure function gfn_hard_wall_contribution( &
    q, &
    hard_wall_locations, &
    lambda, dt, &
    slice_start, slice_end,  &
    particle_number, &
    n_slices, n_particles, n_dimensions &
  ) result ( ln_gfn )
  integer, intent(in) :: slice_start, slice_end, particle_number
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, intent(in) :: dt, lambda
  double precision, dimension ( n_slices, n_particles, n_dimensions ), intent(in) :: q
  double precision, dimension ( n_dimensions ), intent(in) :: hard_wall_locations

  double precision :: ln_gfn, tgfn

  double precision :: x,x0,xr,d
  integer :: i,j

  ln_gfn = 0.0d0
  tgfn = 0.0d0

  do i = 1, n_dimensions
    do j = slice_start+1, slice_end
      x0 = q(j-1,particle_number,i)
      x = q(j,particle_number,i)
      d = (hard_wall_locations(i) - x)
      xr = hard_wall_locations(i) + d
      tgfn = exp(( -(xr-x0)**2  + (x-x0)**2)/(2.0d0*lambda*dt))
      ln_gfn = ln_gfn + log(1.0d0 - tgfn)

      d = (hard_wall_locations(i)+ x)
      xr = -hard_wall_locations(i) - d
      tgfn = exp(( -(xr-x0)**2  + (x-x0)**2)/(2.0d0*lambda*dt))
      ln_gfn = ln_gfn + log(1.0d0 - tgfn)
    end do
  end do

end function
!@-node:gcross.20090828095451.1453:hard wall
!@+node:gcross.20090916153857.1828:compute_green_fn_from_distances
! image approximation
pure function compute_green_fn_from_distances( &
    distances, &
    denominator, &
    slice_start, slice_end, &
    n_slices &
  ) result ( gfn )
  integer, intent(in) :: slice_start, slice_end, n_slices
  double precision, intent(in) :: denominator
  double precision, dimension ( n_slices ), intent(in) :: distances
  double precision :: gfn

  integer :: center_slice_number

  center_slice_number = n_slices / 2

  if(slice_start <= center_slice_number .and. slice_end >= center_slice_number) then  
    gfn = 1d0 &
      * product(compute_slice_gfn( &
          distances(slice_start:center_slice_number-1), &
          distances(slice_start+1:center_slice_number) &
        )) &
      * product(compute_slice_gfn( &
          distances(center_slice_number+1:slice_end-1), &
          distances(center_slice_number+2:slice_end) &
        ))
  else
    gfn = &
        product(compute_slice_gfn( &
          distances(slice_start:slice_end-1), &
          distances(slice_start+1:slice_end) &
        ))
  end if

contains

  elemental function compute_slice_gfn(d1,d2) result (slice_gfn)
    double precision, intent(in) :: d1, d2
    double precision :: slice_gfn
    slice_gfn = 1d0 - exp(-d1*d2/denominator)
  end function

end function
!@-node:gcross.20090916153857.1828:compute_green_fn_from_distances
!@-others

end module gfn
!@-node:gcross.20090812093015.1722:@thin gfn.f95
!@-leo
