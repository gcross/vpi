!@+leo-ver=4-thin
!@+node:gcross.20090807144330.2245:@thin angular_momentum.f90
!@@language fortran90
module angular_momentum

 implicit none

 !@ << Constants >>
 !@+node:gcross.20090807144330.2252:<< Constants >>
 integer, parameter :: x_axis_label = 1, y_axis_label = 2, z_axis_label = 3
 !@-node:gcross.20090807144330.2252:<< Constants >>
 !@nl

 contains

!@+others
!@+node:gcross.20090721121051.1757:get_rotation_plane_axes
subroutine get_rotation_plane_axes(rotation_axis,plane_axis_1,plane_axis_2)
  integer, intent(in) :: rotation_axis
  integer, intent(out) :: plane_axis_1, plane_axis_2

  select case (rotation_axis)
    case (x_axis_label)
      plane_axis_1 = y_axis_label
      plane_axis_2 = z_axis_label
    case (y_axis_label)
      plane_axis_1 = x_axis_label
      plane_axis_2 = z_axis_label
    case (z_axis_label)
      plane_axis_1 = x_axis_label
      plane_axis_2 = y_axis_label
    case default
      print *, "Rotation axis ", rotation_axis, " is invalid;  must be X (1), Y (2) or Z (3)."
      stop
  end select

end subroutine get_rotation_plane_axes
!@-node:gcross.20090721121051.1757:get_rotation_plane_axes
!@+node:gcross.20090803153449.1835:perform_special_matmul
subroutine perform_special_matmul(vector,amplitudes,size)
  double complex, dimension(size), intent(inout) :: vector
  double complex, dimension(size), intent(in) :: amplitudes
  integer, intent(in) :: size

  integer :: i
  double complex :: partial_sum

  partial_sum = 0D0
  do i = 1,size
    partial_sum = partial_sum + vector(i)
    vector(i) = partial_sum * amplitudes(i)
  end do

end subroutine perform_special_matmul
!@-node:gcross.20090803153449.1835:perform_special_matmul
!@-others

end module angular_momentum
!@-node:gcross.20090807144330.2245:@thin angular_momentum.f90
!@-leo
