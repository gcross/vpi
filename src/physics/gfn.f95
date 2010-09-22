!@+leo-ver=4-thin
!@+node:gcross.20090812093015.1722:@thin gfn.f95
!@@language fortran90
!@@tabwidth -2

module gfn

  implicit none

contains 

!@+others
!@+node:gcross.20100106123346.1699:2nd order
!@+node:gcross.20100106123346.1698:initialize 2nd order weights
subroutine initialize_2nd_order_weights(n_slices,U_weight)
  integer, intent(in) :: n_slices
  double precision, dimension(n_slices), intent(out) :: U_weight
  integer :: center_slice, ii

  center_slice = n_slices / 2

  if (mod(center_slice,2) .ne. 1) then
    print *,"ERROR: center_slice = n_slices/2 must be odd"
    stop
  end if

  U_weight = 1.0

  U_weight(1) = 0.5d0
  U_weight(n_slices) = 0.5d0
  U_weight(center_slice) = 0.5d0
  U_weight(center_slice+1) = 0.5d0
end subroutine
!@nonl
!@-node:gcross.20100106123346.1698:initialize 2nd order weights
!@+node:gcross.20090624144408.1799:2nd order
pure function gfn2_sp( sl_start, sl_end, U, U_weight, n_slices, n_particles, dt ) result ( ln_gfn )
  double precision, dimension(n_particles,n_slices), intent(in) :: U
  double precision, dimension(n_slices), intent(in) :: U_weight
  integer, intent(in) :: sl_start, sl_end
  integer, intent(in) :: n_slices, n_particles
  double precision, intent(in) :: dt
  double precision :: ln_gfn
  integer :: slice

  ln_gfn = 0
  do slice = sl_start, sl_end
    ln_gfn = ln_gfn + sum(U(:,slice)) * U_weight(slice)
  end do
  ln_gfn = -dt*ln_gfn

end function gfn2_sp
!@nonl
!@-node:gcross.20090624144408.1799:2nd order
!@-node:gcross.20100106123346.1699:2nd order
!@+node:gcross.20100106123346.1700:4th order
!@+node:gcross.20090812093015.1753:initialize 4th order weights
subroutine initialize_4th_order_weights(n_slices,U_weight,gU2_weight)
  integer, intent(in) :: n_slices
  double precision, dimension(n_slices), intent(out) :: U_weight, gU2_weight
  integer :: center_slice, ii

  center_slice = n_slices / 2

  if (mod(center_slice,2) .ne. 1) then
    print *,"ERROR: center_slice = n_slices/2 must be odd"
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
! in order for paths to be broken properly at the center we need center_slice 
! to be odd.
!@-at
!@@c
  do ii = 1, center_slice
    U_weight(ii) = mod(ii+1,2) + 1 
    gU2_weight(ii) = mod(ii+1,2)
  end do
  ! the center slice is doubled to handle broken paths 
  ! (probably a less kludgy way to do this but... )
  do ii = center_slice+2,n_slices
    U_weight(ii) = mod(ii,2) + 1 
    gU2_weight(ii) = mod(ii,2)
  end do
  U_weight(1) = 0.5d0
  U_weight(n_slices) = 0.5d0
  U_weight(center_slice) = 0.5d0
  U_weight(center_slice+1) = 0.5d0
end subroutine
!@nonl
!@-node:gcross.20090812093015.1753:initialize 4th order weights
!@+node:gcross.20090812093015.1845:4th order
pure function gfn4_sp( &
  sl_start, sl_end, &
  U, gradU2, &
  U_weight, gU2_weight, &
  n_slices, n_particles, &
  lambda, dt &
) result ( ln_gfn )
  double precision, dimension(n_particles,n_slices), intent(in) :: U
  double precision, dimension(n_slices), intent(in) :: gradU2
  double precision, dimension(n_slices), intent(in) :: U_weight
  double precision, dimension(n_slices), intent(in) :: gU2_weight
  integer, intent(in) :: sl_start, sl_end
  integer, intent(in) :: n_slices, n_particles
  double precision, intent(in) :: lambda, dt
  double precision :: ln_gfn, factor
  integer :: slice

  ln_gfn = 0

  do slice = sl_start, sl_end
    ln_gfn = ln_gfn + sum(U(:,slice)) * U_weight(slice)
  end do
  ln_gfn = -2.0d0*dt/3.0d0 * ln_gfn - 2.0d0*lambda*(dt**3)*dot_product(gradU2(sl_start:sl_end),gU2_weight(sl_start:sl_end))/9.0d0

end function gfn4_sp
!@nonl
!@-node:gcross.20090812093015.1845:4th order
!@-node:gcross.20100106123346.1700:4th order
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
!@nonl
!@-node:gcross.20090916153857.1828:compute_green_fn_from_distances
!@-others

end module gfn
!@-node:gcross.20090812093015.1722:@thin gfn.f95
!@-leo
