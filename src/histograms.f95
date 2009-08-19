!@+leo-ver=4-thin
!@+node:gcross.20090819083142.1363:@thin histograms.f95
!@@language fortran90
module histograms

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
!@+node:gcross.20090819083142.1370:accumulate_radial_densities
pure subroutine accumulate_radial_densities(x,maximum_radius,n_particles,n_dimensions,n_bins,histogram)
  double precision, dimension(n_particles,n_dimensions), intent(in) :: x
  double precision, intent(in) :: maximum_radius
  integer, intent(in) :: n_particles, n_dimensions, n_bins
  integer, dimension(n_bins), intent(inout) :: histogram

  integer :: i
  double precision :: radius, dndx

  dndx = n_bins/maximum_radius
  do i = 1, n_particles
    call place_in_bin(sqrt(sum(x(i,:)**2)),0d0,dndx,n_bins,histogram)
  end do

end subroutine accumulate_radial_densities
!@-node:gcross.20090819083142.1370:accumulate_radial_densities
!@-others

end module
!@nonl
!@-node:gcross.20090819083142.1363:@thin histograms.f95
!@-leo
