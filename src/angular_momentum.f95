!@+leo-ver=4-thin
!@+node:gcross.20090807144330.2245:@thin angular_momentum.f95
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
pure subroutine get_rotation_plane_axes(rotation_axis,plane_axis_1,plane_axis_2)
  integer, intent(in) :: rotation_axis
  integer, intent(out) :: plane_axis_1, plane_axis_2

  integer, dimension(3), parameter :: &
    plane_axis_1_list = (/ 2, 1, 1/), &
    plane_axis_2_list = (/ 3, 3, 2/)

  plane_axis_1 = plane_axis_1_list(rotation_axis)
  plane_axis_2 = plane_axis_2_list(rotation_axis)

end subroutine get_rotation_plane_axes
!@-node:gcross.20090721121051.1757:get_rotation_plane_axes
!@+node:gcross.20090803153449.1835:perform_special_matmul
pure subroutine perform_special_matmul(vector,amplitudes,size)
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
!@+node:gcross.20090803153449.1836:sum_over_symmetrizations
pure function sum_over_symmetrizations(amplitudes,N_particles,N_excited)
  ! Input variables
  integer, intent(in) :: N_particles, N_excited
  double complex, dimension(N_particles), intent(in) :: amplitudes

  ! Function Result
  double complex :: sum_over_symmetrizations

  ! Local variables
  double complex, dimension(N_particles-N_excited+1) :: vector
  integer :: i


  vector = amplitudes(1:N_particles-N_excited+1)
  do i = 2, N_excited
    call perform_special_matmul(vector,amplitudes(i:i+N_particles-N_excited),N_particles-N_excited+1)
  end do
  sum_over_symmetrizations = sum(vector)

end function sum_over_symmetrizations
!@-node:gcross.20090803153449.1836:sum_over_symmetrizations
!@+node:gcross.20090803153449.1837:compute_angular_derivatives
pure subroutine compute_angular_derivatives(&
    x, &
    fixed_rotation_axis, fixed_angular_momentum, &
    N_particles,N_dimensions, &
    derivatives &
  )

  ! Input variables
  integer, intent(in) :: N_particles, N_dimensions
  double precision, dimension(N_particles,N_dimensions), intent(in) :: x
  integer, intent(in) :: fixed_rotation_axis, fixed_angular_momentum

  ! Output variables
  double precision, dimension(N_particles), intent(out) :: derivatives

  ! Local variables
  integer :: plane_axis_1, plane_axis_2
  double complex, dimension(N_particles) :: amplitudes, partial_sums
  double complex :: amplitude
  integer :: i

  if(fixed_angular_momentum == 0) then
    derivatives(:) = 0
    return
  end if

  if(N_particles == fixed_angular_momentum) then
    derivatives(:) = 1
    return
  end if

  call get_rotation_plane_axes(fixed_rotation_axis,plane_axis_1,plane_axis_2)

  amplitudes(:) = &
    (x(:,plane_axis_1)*(1.0d0,0) + x(:,plane_axis_2)*(0,1.0d0))/sqrt(x(:,plane_axis_1)**2 + x(:,plane_axis_2)**2)

!@+at
! To get the partial sums, i.e. the the sums which are required to have a 
! selected coordinate "i" appear, we compute the full sum but set the 
! amplitude of the coordinate "i" to zero in order to get a sum over all sets 
! of (m-1) coordinates which don't include "i"'s amplitude, then we multiply 
! by "i"'s amplitude.
!@-at
!@@c
  do i = 1, N_particles
    amplitude = amplitudes(i)
    amplitudes(i) = (0d0,0d0)
    partial_sums(i) = sum_over_symmetrizations(amplitudes,N_particles,fixed_angular_momentum)*amplitude
    amplitudes(i) = amplitude
  end do
!@+at
! The full sum over all sets of m coordinates is equal to the sum of all of 
! the partial sums divided by m to account for the fact that each term appears 
! m times, since the term for any given subset of m coordinates appears in m 
! of the partial sums (one for each coordinate).
!@-at
!@@c
  amplitude = sum(partial_sums)/fixed_angular_momentum
  derivatives(:) = real(amplitude*conjg(partial_sums(:)))/(abs(amplitude)**2)

end subroutine compute_angular_derivatives
!@-node:gcross.20090803153449.1837:compute_angular_derivatives
!@+node:gcross.20090721121051.1756:compute_effective_rotational_potential
pure subroutine compute_effective_rotational_potential (&
    x, &
    fixed_rotation_axis, frame_angular_velocity, fixed_angular_momentum, &
    move_start, move_end, &
    n_slices, n_particles, n_dimensions, &
    U_rot &
  )

  ! Input variables
  integer, intent(in) :: move_start, move_end, n_slices, n_particles, n_dimensions
  double precision, dimension ( n_slices, n_particles , n_dimensions ), intent(in) :: x
  double precision, intent(in) :: frame_angular_velocity
  integer, intent(in) :: fixed_rotation_axis, fixed_angular_momentum

  ! Output variables
  double precision, dimension( n_slices, n_particles ), intent(inout) :: U_rot

  ! Local variables
  integer :: i

  forall (i = move_start:move_end) &
    U_rot(i,:) = U_rot(i,:) + compute_for_slice(i)  

  contains

    pure function compute_for_slice (i) result (potential)
      integer, intent(in) :: i
      double precision, dimension( n_particles ) :: derivatives, potential
      call compute_angular_derivatives( &
              x(i,:,:), &
              fixed_rotation_axis, fixed_angular_momentum, &
              n_particles, n_dimensions, &
              derivatives &
            )
      potential = derivatives(:)*(derivatives(:)-frame_angular_velocity)
    end function compute_for_slice

end subroutine compute_effective_rotational_potential
!@-node:gcross.20090721121051.1756:compute_effective_rotational_potential
!@-others

end module angular_momentum
!@-node:gcross.20090807144330.2245:@thin angular_momentum.f95
!@-leo
