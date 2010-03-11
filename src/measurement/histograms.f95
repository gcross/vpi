!@+leo-ver=4-thin
!@+node:gcross.20090819083142.1363:@thin histograms.f95
!@@language fortran90
module histograms

use constants

implicit none

contains

!@+others
!@+node:gcross.20090819083142.1365:within_bins
pure function within_bins(index,nbins)
  integer, intent(in) :: index, nbins
  logical :: within_bins

  within_bins = (index >= 1) .and. (index <= nbins)
end function within_bins
!@-node:gcross.20090819083142.1365:within_bins
!@+node:gcross.20090819083142.1371:place_in_bin
pure subroutine place_in_bin(value,offset,dndx,n_bins,histogram)
  double precision, intent(in) :: value, offset, dndx
  integer, intent(in) :: n_bins
  integer, dimension(n_bins), intent(inout) :: histogram

  integer bin

  bin = floor((value+offset)*dndx)+1
  if ( within_bins(bin,n_bins) ) then
    histogram(bin) = histogram(bin) + 1
  end if

end subroutine place_in_bin
!@-node:gcross.20090819083142.1371:place_in_bin
!@+node:gcross.20090826112349.1420:accumulate
pure subroutine accumulate(x,left_x,right_x,n_values,n_bins,histogram)
  double precision, dimension(n_values), intent(in) :: x
  double precision, intent(in) :: left_x, right_x
  integer, intent(in) :: n_values, n_bins
  integer, dimension(n_bins), intent(inout) :: histogram

  integer :: i
  double precision :: offset, dndx

  offset = -left_x
  dndx = n_bins/(right_x-left_x)
  do i = 1, n_values
    call place_in_bin(x(i),offset,dndx,n_bins,histogram)
  end do

end subroutine
!@-node:gcross.20090826112349.1420:accumulate
!@+node:gcross.20090819083142.1364:accumulate_1d_densities
pure subroutine accumulate_1d_densities(x,left_x,right_x,n_particles,n_dimensions,n_bins,histogram)
  double precision, dimension(n_particles,n_dimensions), intent(in) :: x
  double precision, dimension(n_dimensions), intent(in) :: left_x, right_x
  integer, intent(in) :: n_particles, n_dimensions, n_bins
  integer, dimension(n_dimensions,n_bins), intent(inout) :: histogram

  integer :: i,j
  double precision :: offset, dndx

  do j = 1, n_dimensions
    offset = -left_x(j)
    dndx = n_bins/(right_x(j)-left_x(j))
    do i = 1, n_particles
      call place_in_bin(x(i,j),offset,dndx,n_bins,histogram(j,:))
    end do
  end do

end subroutine accumulate_1d_densities
!@-node:gcross.20090819083142.1364:accumulate_1d_densities
!@+node:gcross.20100226131523.1656:accumulate_2d_density
subroutine accumulate_2d_density(positions,left_x,left_y,right_x,right_y,n_particles,n_bins,histogram)
  double precision, dimension(n_particles,2), intent(in) :: positions
  double precision, intent(in) :: left_x, left_y, right_x, right_y
  integer, intent(in) :: n_particles, n_bins
  integer, dimension(n_bins,n_bins), intent(inout) :: histogram

  integer :: i, bin_x, bin_y
  double precision :: offset, dndx, dndy
  dndx = n_bins/(right_x-left_x)
  dndy = n_bins/(right_y-left_y)

  do i = 1, n_particles
    bin_x = floor((positions(i,1)-left_x)*dndx)+1
    bin_y = floor((positions(i,2)-left_y)*dndy)+1
    if ( within_bins(bin_x,n_bins) .and. within_bins(bin_y,n_bins) ) then
      histogram(bin_x,bin_y) = histogram(bin_x,bin_y) + 1
    end if
  end do

end subroutine accumulate_2d_density
!@-node:gcross.20100226131523.1656:accumulate_2d_density
!@+node:gcross.20100226131523.1658:accumulate_2d_density_matrix
pure subroutine accumulate_2d_density_matrix(positions_1,positions_2,left_x,left_y,right_x,right_y,n_particles,n_bins,histogram)
  double precision, dimension(n_particles,2), intent(in) :: positions_1, positions_2
  double precision, intent(in) :: left_x, left_y, right_x, right_y
  integer, intent(in) :: n_particles, n_bins
  integer, dimension(n_bins,n_bins,n_bins,n_bins), intent(inout) :: histogram

  integer :: i, bin_1_x, bin_1_y, bin_2_x, bin_2_y
  double precision :: offset, dndx, dndy
  dndx = n_bins/(right_x-left_x)
  dndy = n_bins/(right_y-left_y)

  do i = 1, n_particles
    bin_1_x = floor((positions_1(i,1)-left_x)*dndx)+1
    bin_1_y = floor((positions_1(i,2)-left_y)*dndy)+1
    bin_2_x = floor((positions_2(i,1)-left_x)*dndx)+1
    bin_2_y = floor((positions_2(i,2)-left_y)*dndy)+1
    if ( within_bins(bin_1_x,n_bins) .and. within_bins(bin_1_y,n_bins) .and. &
         within_bins(bin_2_x,n_bins) .and. within_bins(bin_2_y,n_bins) &
       ) then
      histogram(bin_1_x,bin_1_y,bin_2_x,bin_2_y) = histogram(bin_1_x,bin_1_y,bin_2_x,bin_2_y) + 1
    end if
  end do

end subroutine accumulate_2d_density_matrix
!@-node:gcross.20100226131523.1658:accumulate_2d_density_matrix
!@+node:gcross.20090819083142.1370:accumulate_radial_densities
pure subroutine accumulate_radial_densities(x,maximum_radius,n_particles,n_dimensions,n_bins,histogram)
  double precision, dimension(n_particles,n_dimensions), intent(in) :: x
  double precision, intent(in) :: maximum_radius
  integer, intent(in) :: n_particles, n_dimensions, n_bins
  integer, dimension(n_bins), intent(inout) :: histogram

  integer :: i
  double precision :: dndx

  dndx = n_bins/maximum_radius
  do i = 1, n_particles
    call place_in_bin(sqrt(sum(x(i,:)**2)),0d0,dndx,n_bins,histogram)
  end do

end subroutine
!@-node:gcross.20090819083142.1370:accumulate_radial_densities
!@+node:gcross.20090825141639.1523:accumulate_recip_r_sq_densities
pure subroutine accumulate_recip_r_sq_densities(x,maximum_value,n_particles,n_dimensions,n_bins,histogram)
  double precision, dimension(n_particles,n_dimensions), intent(in) :: x
  double precision, intent(in) :: maximum_value
  integer, intent(in) :: n_particles, n_dimensions, n_bins
  integer, dimension(n_bins), intent(inout) :: histogram

  integer :: i
  double precision :: dndx

  dndx = n_bins/maximum_value
  do i = 1, n_particles
    call place_in_bin(1.0/sum(x(i,:)**2),0d0,dndx,n_bins,histogram)
  end do

end subroutine
!@-node:gcross.20090825141639.1523:accumulate_recip_r_sq_densities
!@+node:gcross.20090825141639.1525:accumulate_plane_radial_densities
pure subroutine accumulate_plane_radial_densities( &
    x,maximum_radius, &
    plane_axis_1,plane_axis_2, &
    n_particles,n_dimensions,n_bins, &
    histogram &
    )
  double precision, dimension(n_particles,n_dimensions), intent(in) :: x
  double precision, intent(in) :: maximum_radius
  integer, intent(in) :: plane_axis_1, plane_axis_2, n_particles, n_dimensions, n_bins
  integer, dimension(n_bins), intent(inout) :: histogram

  integer :: i
  double precision :: dndx

  dndx = n_bins/maximum_radius
  do i = 1, n_particles
    call place_in_bin(sqrt(x(i,plane_axis_1)**2+x(i,plane_axis_2)**2),0d0,dndx,n_bins,histogram)
  end do

end subroutine
!@-node:gcross.20090825141639.1525:accumulate_plane_radial_densities
!@+node:gcross.20090826091347.1416:accumulate_angular_separation_densities
pure subroutine accumulate_angular_separation_densities( &
    angles, &
    n_particles,n_bins, &
    histogram &
    )
  double precision, dimension(n_particles), intent(in) :: angles
  integer, intent(in) :: n_particles, n_bins
  integer, dimension(n_bins), intent(inout) :: histogram

  integer :: i, j
  double precision :: dndtheta

  dndtheta = n_bins/m_2pi
  do i = 1, n_particles
    do j = i+1, n_particles
      call place_in_bin(abs(angles(i)-angles(j)),0d0,dndtheta,n_bins,histogram)
    end do
  end do

end subroutine
!@-node:gcross.20090826091347.1416:accumulate_angular_separation_densities
!@+node:gcross.20090826091347.1424:accumulate_neighbor_angular_separation_densities
pure subroutine accumulate_neighbor_angular_separation_densities( &
    angles, &
    n_particles,n_bins, &
    histogram &
    )
  double precision, dimension(n_particles), intent(in) :: angles
  integer, intent(in) :: n_particles, n_bins
  integer, dimension(n_bins), intent(inout) :: histogram

  integer :: i
  double precision :: dndtheta

  dndtheta = n_bins/m_2pi
  do i = 1, n_particles
    call place_in_bin( &
      min(minval(abs(angles(i)-angles(:i-1))),minval(abs(angles(i)-angles(i+1:)))), &
      0d0,dndtheta,n_bins,histogram &
    )
  end do

end subroutine
!@-node:gcross.20090826091347.1424:accumulate_neighbor_angular_separation_densities
!@+node:gcross.20090825141639.1527:accumulate_recip_plane_r_sq_densities
pure subroutine accumulate_recip_plane_r_sq_densities( &
    x,maximum_value, &
    plane_axis_1,plane_axis_2, &
    n_particles,n_dimensions,n_bins, &
    histogram &
    )
  double precision, dimension(n_particles,n_dimensions), intent(in) :: x
  double precision, intent(in) :: maximum_value
  integer, intent(in) :: plane_axis_1, plane_axis_2, n_particles, n_dimensions, n_bins
  integer, dimension(n_bins), intent(inout) :: histogram

  integer :: i
  double precision :: dndx

  dndx = n_bins/maximum_value
  do i = 1, n_particles
    call place_in_bin(1.0/(x(i,plane_axis_1)**2+x(i,plane_axis_2)**2),0d0,dndx,n_bins,histogram)
  end do

end subroutine
!@-node:gcross.20090825141639.1527:accumulate_recip_plane_r_sq_densities
!@+node:gcross.20090830224709.2035:accumulate_particle_separation_densities
pure subroutine accumulate_particle_separation_densities( &
    xij2, &
    maximum_value, &
    n_particles,n_bins, &
    histogram &
    )
  double precision, dimension(n_particles,n_particles), intent(in) :: xij2
  double precision, intent(in) :: maximum_value
  integer, intent(in) :: n_particles, n_bins
  integer, dimension(n_bins), intent(inout) :: histogram

  integer :: i, j
  double precision :: dndx

  dndx = n_bins/maximum_value
  do i = 1, n_particles
    do j = i+1, n_particles
      call place_in_bin(sqrt(xij2(j,i)),0d0,dndx,n_bins,histogram)
    end do
  end do

end subroutine
!@-node:gcross.20090830224709.2035:accumulate_particle_separation_densities
!@-others

end module
!@nonl
!@-node:gcross.20090819083142.1363:@thin histograms.f95
!@-leo
