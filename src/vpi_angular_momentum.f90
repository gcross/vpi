!@+leo-ver=4-thin
!@+node:gcross.20090724102123.1771:@thin vpi_angular_momentum.f90
!@@language fortran90
module vpi_angular_momentum

 !@ << Imported modules >>
 !@+node:gcross.20090724102123.1772:<< Imported modules >>
 use vpi_defines
 !@-node:gcross.20090724102123.1772:<< Imported modules >>
 !@nl

 implicit none

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
  complex(kind=b8), dimension(size), intent(inout) :: vector
  complex(kind=b8), dimension(size), intent(in) :: amplitudes
  integer, intent(in) :: size

  integer :: i
  complex(kind=b8) :: partial_sum

  partial_sum = 0_b8  
  do i = 1,size
    partial_sum = partial_sum + vector(i)
    vector(i) = partial_sum * amplitudes(i)
  end do

end subroutine perform_special_matmul
!@-node:gcross.20090803153449.1835:perform_special_matmul
!@+node:gcross.20090803153449.1836:sum_over_symmetrizations
function sum_over_symmetrizations(amplitudes,N_particles,N_excited)
  ! Input variables
  integer, intent(in) :: N_particles, N_excited
  complex(kind=b8), dimension(N_particles), intent(in) :: amplitudes

  ! Function Result
  complex(kind=b8) :: sum_over_symmetrizations

  ! Local variables
  complex(kind=b8), dimension(N_particles-N_excited+1) :: vector
  integer :: i


  vector = amplitudes(1:N_particles-N_excited+1)
  do i = 2, N_excited
    call perform_special_matmul(vector,amplitudes(i:i+N_particles-N_excited),N_particles-N_excited+1)
  end do
  sum_over_symmetrizations = sum(vector)

end function sum_over_symmetrizations
!@-node:gcross.20090803153449.1836:sum_over_symmetrizations
!@+node:gcross.20090803153449.1837:compute_angular_derivatives
subroutine compute_angular_derivatives(&
! INPUT: particle position information
    x, &
! INPUT: array dimensions
    N_particles,N_dimensions, &
! OUTPUT: partial derivatives of the phase function with respect to each coordinate
    derivatives &
  )

  ! Input variables
  integer, intent(in) :: N_particles, N_dimensions
  real(kind=b8), dimension(N_particles,N_dimensions), target, intent(in) :: x

  ! Output variables
  real(kind=b8), dimension(N_particles), intent(out) :: derivatives

  ! Local variables
  integer :: plane_axis_1, plane_axis_2
  real(kind=b8), dimension(:), pointer :: x_slice_1, x_slice_2
  complex(kind=b8), dimension(N_particles) :: amplitudes, partial_sums
  complex(kind=b8) :: amplitude
  integer :: i


  call get_rotation_plane_axes(fixed_rotation_axis,plane_axis_1,plane_axis_2)
  x_slice_1 => x(:,plane_axis_1)
  x_slice_2 => x(:,plane_axis_2)

  amplitudes(:) = &
    (x_slice_1(:)*(1.0_b8,0) + x_slice_2(:)*(0,1.0_b8))/sqrt(x_slice_1(:)**2 + x_slice_2(:)**2)

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
    amplitudes(i) = 0
    partial_sums(i) = sum_over_symmetrizations(amplitudes,N_particles,fixed_angular_momentum-1)*amplitude
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
subroutine compute_effective_rotational_potential (&
! INPUT: particle position information
    x, &
! INPUT: path slice to consider
    move_start, move_end, &
! INPUT: array dimensions
    nslice, np, ndim, &
! OUTPUT: effective rotational potential
    U_rot &
  )

  ! Input variables
  integer, intent(in) :: move_start, move_end, nslice, np, ndim
  real(kind=b8), dimension ( nslice, np , ndim ), intent(in) :: x

  ! Output variables
  real(kind=b8), dimension( nslice, np ), intent(out) :: U_rot

  ! Local variables
  integer :: i
  real(kind=b8), dimension( np ) :: derivatives

  do i = move_start, move_end
    call compute_angular_derivatives(x(i,:,:),np,ndim,derivatives)
    U_rot(i,:) = derivatives(:)*(derivatives(:)-frame_angular_velocity)
  end do
end subroutine compute_effective_rotational_potential
!@-node:gcross.20090721121051.1756:compute_effective_rotational_potential
!@-others

end module vpi_angular_momentum
!@-node:gcross.20090724102123.1771:@thin vpi_angular_momentum.f90
!@-leo
