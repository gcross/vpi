!@+leo-ver=4-thin
!@+node:gcross.20090807144330.2245:@thin angular_momentum.f95
!@@language fortran90
module angular_momentum

 use numeric_differentiation

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
pure function sum_over_symmetrizations(amplitudes,N_particles,N_rotating_particles)
  ! Input variables
  integer, intent(in) :: N_particles, N_rotating_particles
  double complex, dimension(N_particles), intent(in) :: amplitudes

  ! Function Result
  double complex :: sum_over_symmetrizations

  ! Local variables
  double complex, dimension(N_particles-N_rotating_particles+1) :: vector
  integer :: i

  if(n_rotating_particles == 0) then
    sum_over_symmetrizations = 1
    return
  end if

  vector = amplitudes(1:N_particles-N_rotating_particles+1)
  do i = 2, N_rotating_particles
    call perform_special_matmul(vector,amplitudes(i:i+N_particles-N_rotating_particles),N_particles-N_rotating_particles+1)
  end do
  sum_over_symmetrizations = sum(vector)

end function
!@-node:gcross.20090803153449.1836:sum_over_symmetrizations
!@+node:gcross.20090827135657.1422:accumulate_gradient_fancy
pure subroutine accumulate_gradient_fancy (&
    x, &
    n_rotating_particles, &
    rotation_plane_axis_1, rotation_plane_axis_2, &
    n_slices, n_particles, n_dimensions, &
    gradient_phase &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension(n_slices, n_particles, n_dimensions), intent(in) :: x
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2, n_rotating_particles

  ! Output variables
  double precision, dimension(n_slices, n_particles, n_dimensions), intent(inout) :: gradient_phase

  ! Local variables
  double complex, dimension(n_particles) :: amplitudes
  double complex :: full_sum, full_sum_conj, partial_sum
  double precision :: full_sum_amplitude_squared
  integer :: s, i

  !@  << Special cases >>
  !@+node:gcross.20090915142144.1647:<< Special cases >>
  if(n_rotating_particles == 0) then
    return
  end if

  if(mod(n_particles,2) == 0 .and. n_rotating_particles == N_particles / 2) then
    do i = 1, n_particles
      do s = 1, n_slices
        call accumulate_gradient(0.5d0,gradient_phase(s,i,:))
      end do
    end do
    return
  end if

  if(n_rotating_particles == n_particles) then
    do i = 1, n_particles
      do s = 1, n_slices
        call accumulate_gradient(1.0d0,gradient_phase(s,i,:))
      end do
    end do
    return
  end if
  !@-node:gcross.20090915142144.1647:<< Special cases >>
  !@nl

  do s = 1, n_slices

    amplitudes(:) = &
      (x(s,:,rotation_plane_axis_1)*(1.0d0,0) + x(s,:,rotation_plane_axis_2)*(0,1.0d0)) &
    / sqrt(x(s,:,rotation_plane_axis_1)**2 + x(s,:,rotation_plane_axis_2)**2)
    full_sum = sum_over_symmetrizations(amplitudes,n_particles,n_rotating_particles)
    full_sum_conj = conjg(full_sum)
    full_sum_amplitude_squared = real(full_sum*full_sum_conj)

    if(n_rotating_particles == 1) then
      do i = 1, n_particles
        call compute_and_accumulate_gradient(amplitudes(i),gradient_phase(s,i,:))
      end do
      cycle
    else
      do i = 1, n_particles
        call compute_partial_sum(amplitudes,i,n_particles,n_rotating_particles,partial_sum)
        call compute_and_accumulate_gradient(partial_sum,gradient_phase(s,i,:))
      end do
    end if
  end do

contains

  !@  << Helper routines >>
  !@+node:gcross.20090915142144.1648:<< Helper routines >>
  pure subroutine compute_and_accumulate_gradient(amplitude,gradient_phase)
    double complex, intent(in) :: amplitude
    double precision, dimension(n_dimensions), intent(inout) :: gradient_phase
    double precision :: phase_derivative

    phase_derivative = real(amplitude*full_sum_conj)/full_sum_amplitude_squared

    call accumulate_gradient(phase_derivative,gradient_phase)
  end subroutine

  pure subroutine accumulate_gradient(phase_derivative,gradient_phase)
    double precision, intent(in) :: phase_derivative
    double precision, dimension(n_dimensions), intent(inout) :: gradient_phase
    call accum_angle_drv_into_gradient( &
          phase_derivative, x(s,i,:), &
          rotation_plane_axis_1, rotation_plane_axis_2, &
          n_dimensions, &
          gradient_phase &
    )
  end subroutine
  !@-node:gcross.20090915142144.1648:<< Helper routines >>
  !@nl

end subroutine
!@-node:gcross.20090827135657.1422:accumulate_gradient_fancy
!@+node:gcross.20090915142144.1652:compute_gradient_fancy_amplitude
pure subroutine compute_gradient_fancy_amplitude (&
    x, &
    n_rotating_particles, &
    rotation_plane_axis_1, rotation_plane_axis_2, &
    n_particles, n_dimensions, &
    gradient_amplitude &
  )

  ! Input variables
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_particles, n_dimensions), intent(in) :: x
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2, n_rotating_particles

  ! Output variables
  double precision, dimension(n_particles, n_dimensions), intent(out) :: gradient_amplitude

  ! Local variables
  double complex, dimension(n_particles) :: amplitudes
  double complex :: full_sum, full_sum_conj, partial_sum
  double precision :: full_sum_amplitude_squared
  integer :: i

  !@  << Special cases >>
  !@+node:gcross.20090915142144.1653:<< Special cases >>
  if(mod(n_rotating_particles,n_particles) == 0) then
    gradient_amplitude = 0
    return
  end if
  !@-node:gcross.20090915142144.1653:<< Special cases >>
  !@nl

  amplitudes(:) = &
    (x(:,rotation_plane_axis_1)*(1.0d0,0) + x(:,rotation_plane_axis_2)*(0,1.0d0)) &
  / sqrt(x(:,rotation_plane_axis_1)**2 + x(:,rotation_plane_axis_2)**2)
  full_sum = sum_over_symmetrizations(amplitudes,n_particles,n_rotating_particles)
  full_sum_conj = conjg(full_sum)
  if(n_rotating_particles == 1) then
    forall (i=1:n_particles) &
      gradient_amplitude(i,:) = compute_gradient(amplitudes(i),i)
  else
    do i = 1, n_particles
      call compute_partial_sum(amplitudes,i,n_particles,n_rotating_particles,partial_sum)
      gradient_amplitude(i,:) = compute_gradient(partial_sum,i)
    end do
  end if

contains

  !@  << Helper routines >>
  !@+node:gcross.20090915142144.1654:<< Helper routines >>
  pure function compute_gradient(amplitude,i) result (gradient)
    double complex, intent(in) :: amplitude
    integer, intent(in) :: i
    double precision, dimension(n_dimensions) :: gradient

    double precision :: derivative

    derivative = -2d0*imag(amplitude*full_sum_conj)

    gradient = 0
    call accumulate_gradient(derivative,i,gradient)
  end function

  pure subroutine accumulate_gradient(derivative,i,gradient_phase)
    double precision, intent(in) :: derivative
    integer, intent(in) :: i
    double precision, dimension(n_dimensions), intent(inout) :: gradient_phase
    call accum_angle_drv_into_gradient( &
          derivative, x(i,:), &
          rotation_plane_axis_1, rotation_plane_axis_2, &
          n_dimensions, &
          gradient_phase &
    )
  end subroutine
  !@-node:gcross.20090915142144.1654:<< Helper routines >>
  !@nl

end subroutine
!@-node:gcross.20090915142144.1652:compute_gradient_fancy_amplitude
!@+node:gcross.20090915142144.1645:compute_partial_sum
pure subroutine compute_partial_sum( &
    amplitudes, &
    excluded_index, &
    n_particles, n_rotating_particles, &
    partial_sum &
  )

  double complex, dimension(n_particles), intent(inout) :: amplitudes
  integer, intent(in) :: excluded_index, n_particles, n_rotating_particles
  double complex, intent(out) :: partial_sum

  double complex :: saved_amplitude

!@+at
! To get the partial sums, i.e. the the sums which are required to have a
! selected coordinate "i" appear, we compute the full sum but set the
! amplitude of the coordinate "i" to zero in order to get a sum over all sets
! of (m-1) coordinates which don't include "i"'s amplitude, then we multiply
! by "i"'s amplitude.
!@-at
!@@c

  saved_amplitude = amplitudes(excluded_index)
  amplitudes(excluded_index) = (0d0,0d0)
  partial_sum = sum_over_symmetrizations(amplitudes,n_particles,n_rotating_particles-1)*saved_amplitude
  amplitudes(excluded_index) = saved_amplitude

end subroutine
!@-node:gcross.20090915142144.1645:compute_partial_sum
!@+node:gcross.20090915142144.1644:accum_angle_drv_into_gradient
pure subroutine accum_angle_drv_into_gradient(&
    derivative_of_fn_by_angle, x, &
    rotation_plane_axis_1, rotation_plane_axis_2, &
    n_dimensions, &
    gradient &
  )
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2
  integer, intent(in) :: n_dimensions
  double precision, intent(in) :: derivative_of_fn_by_angle
  double precision, dimension(n_dimensions), intent(in) :: x
  double precision, dimension(n_dimensions), intent(inout) :: gradient
  double precision :: rot_x, rot_y, rot_r_squared

  rot_x = x(rotation_plane_axis_1)
  rot_y = x(rotation_plane_axis_2)
  rot_r_squared = rot_x**2 + rot_y**2

  gradient(rotation_plane_axis_1) = gradient(rotation_plane_axis_1) &
    + derivative_of_fn_by_angle * rot_y/rot_r_squared
  gradient(rotation_plane_axis_2) = gradient(rotation_plane_axis_2) &
    - derivative_of_fn_by_angle * rot_x/rot_r_squared
end subroutine
!@nonl
!@-node:gcross.20090915142144.1644:accum_angle_drv_into_gradient
!@+node:gcross.20090903090230.2072:accumulate_effective_potential
subroutine accumulate_effective_potential (&
    gradient_phase, &
    frame_angular_velocity, lambda, &
    n_slices, n_particles, n_dimensions, &
    U &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension ( n_slices, n_particles, n_dimensions ), intent(in) :: gradient_phase
  double precision, intent(in) :: frame_angular_velocity, lambda

  ! Output variables
  double precision, dimension( n_slices, n_particles ), intent(inout) :: U

  ! Local variables
  integer :: s, i

  forall (s=1:n_slices, i=1:n_particles) &
    U(s,i) = U(s,i) + sum(gradient_phase(s,i,:)*(lambda*gradient_phase(s,i,:)-frame_angular_velocity))

end subroutine
!@-node:gcross.20090903090230.2072:accumulate_effective_potential
!@+node:gcross.20090908085435.1633:accumulate_gradient_feynman
subroutine accumulate_gradient_feynman (&
    x, &
    rotation_rate, &
    rotation_plane_axis_1, rotation_plane_axis_2, &
    n_slices, n_particles, n_dimensions, &
    gradient_phase &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension( n_slices, n_particles, n_dimensions ), intent(in) :: x
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2
  double precision, intent(in) :: rotation_rate

  ! Output variables
  double precision, dimension ( n_slices, n_particles, n_dimensions ), intent(inout) :: gradient_phase

  ! Local variables
  integer :: s, i

  forall (s=1:n_slices, i=1:n_particles)
    gradient_phase(s,i,rotation_plane_axis_1) = gradient_phase(s,i,rotation_plane_axis_1) &
        + rotation_rate * x(s,i,rotation_plane_axis_2) / (x(s,i,rotation_plane_axis_1)**2 + x(s,i,rotation_plane_axis_2)**2)
    gradient_phase(s,i,rotation_plane_axis_2) = gradient_phase(s,i,rotation_plane_axis_2) &
        - rotation_rate * x(s,i,rotation_plane_axis_1) / (x(s,i,rotation_plane_axis_1)**2 + x(s,i,rotation_plane_axis_2)**2)
  end forall

end subroutine
!@-node:gcross.20090908085435.1633:accumulate_gradient_feynman
!@-others

end module angular_momentum
!@-node:gcross.20090807144330.2245:@thin angular_momentum.f95
!@-leo
