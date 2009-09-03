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
!@+node:gcross.20090827135657.1422:compute_angular_derivatives
pure subroutine compute_angular_derivatives(&
    x, &
    rotation_plane_axis_1, rotation_plane_axis_2, N_rotating_particles, &
    N_particles,N_dimensions, &
    derivatives &
  )

  ! Input variables
  integer, intent(in) :: N_particles, N_dimensions
  double precision, dimension(N_particles,N_dimensions), intent(in) :: x
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2, N_rotating_particles

  ! Output variables
  double precision, dimension(N_particles), intent(out) :: derivatives

  ! Local variables
  double complex, dimension(N_particles) :: amplitudes
  double complex :: saved_amplitude, full_sum, full_sum_conj, partial_sum
  double precision :: full_sum_amplitude_squared
  integer :: i

  if(N_rotating_particles == 0) then
    derivatives = 0d0
    return
  end if

  if(mod(N_particles,2) == 0 .and. N_rotating_particles == N_particles / 2) then
    derivatives = 0.5d0
    return
  end if

  if(N_rotating_particles == N_particles) then
    derivatives = 1d0
    return
  end if

  amplitudes(:) = &
    (x(:,rotation_plane_axis_1)*(1.0d0,0) + x(:,rotation_plane_axis_2)*(0,1.0d0)) &
  / sqrt(x(:,rotation_plane_axis_1)**2 + x(:,rotation_plane_axis_2)**2)

  if(N_rotating_particles == 1) then
    full_sum = sum(amplitudes)
    full_sum_conj = conjg(full_sum)
    full_sum_amplitude_squared = real(full_sum*conjg(full_sum))
    derivatives(:) = real(amplitudes(:)*full_sum_conj)/full_sum_amplitude_squared
    return
  end if

!@+at
! To get the partial sums, i.e. the the sums which are required to have a
! selected coordinate "i" appear, we compute the full sum but set the
! amplitude of the coordinate "i" to zero in order to get a sum over all sets
! of (m-1) coordinates which don't include "i"'s amplitude, then we multiply
! by "i"'s amplitude.
!@-at
!@@c
  full_sum = sum_over_symmetrizations(amplitudes,N_particles,N_rotating_particles)
  full_sum_conj = conjg(full_sum)
  full_sum_amplitude_squared = real(full_sum*full_sum_conj)
  do i = 1, N_particles
    saved_amplitude = amplitudes(i)
    amplitudes(i) = (0d0,0d0)
    partial_sum = sum_over_symmetrizations(amplitudes,N_particles,N_rotating_particles-1)*saved_amplitude
    derivatives(i) = real(partial_sum*full_sum_conj)/full_sum_amplitude_squared
    amplitudes(i) = saved_amplitude
  end do

end subroutine
!@-node:gcross.20090827135657.1422:compute_angular_derivatives
!@+node:gcross.20090827135657.1421:compute_rotation_potential_term
pure function compute_rotation_potential_term (&
    x, lambda, &
    first_angular_derivative, &
    rotation_plane_axis_1, rotation_plane_axis_2, frame_angular_velocity, &
    n_dimensions &
  ) result ( U )

  ! Input variables
  integer, intent(in) :: n_dimensions
  double precision, dimension ( n_dimensions ), intent(in) :: x
  double precision, intent(in) :: frame_angular_velocity, lambda
  double precision, intent(in) :: first_angular_derivative
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2

  ! Output variables
  double precision :: U

  ! Local variables
  double precision :: rho_squared

  rho_squared = x(rotation_plane_axis_1)**2 + x(rotation_plane_axis_2)**2
  U = first_angular_derivative*(first_angular_derivative*lambda/rho_squared - frame_angular_velocity)

end function
!@-node:gcross.20090827135657.1421:compute_rotation_potential_term
!@+node:gcross.20090827135657.1424:accumulate_rotation_potential
pure subroutine accumulate_rotation_potential (&
    x, lambda, &
    first_angular_derivatives, &
    rotation_plane_axis_1, rotation_plane_axis_2, frame_angular_velocity, &
    n_particles, n_dimensions, &
    U &
  )

  ! Input variables
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension ( n_particles, n_dimensions ), intent(in) :: x
  double precision, intent(in) :: frame_angular_velocity, lambda
  double precision, dimension( n_particles ), intent(in) :: first_angular_derivatives
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2

  ! Output variables
  double precision, dimension( n_particles ), intent(inout) :: U

  ! Local variables
  integer :: i

  forall (i=1:n_particles) &
    U(i) = U(i) + &
           compute_rotation_potential_term (&
              x(i,:), lambda, &
              first_angular_derivatives(i), &
              rotation_plane_axis_1, rotation_plane_axis_2, frame_angular_velocity, &
              n_dimensions &
            )

end subroutine
!@-node:gcross.20090827135657.1424:accumulate_rotation_potential
!@+node:gcross.20090827135657.1431:accumulate_effective_potential
pure subroutine accumulate_effective_potential (&
    x, lambda, &
    rotation_plane_axis_1, rotation_plane_axis_2, &
    frame_angular_velocity, n_rotating_particles, &
    n_slices, n_particles, n_dimensions, &
    U &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions, n_rotating_particles
  double precision, dimension ( n_slices, n_particles, n_dimensions ), intent(in) :: x
  double precision, intent(in) :: frame_angular_velocity, lambda
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2

  ! Output variables
  double precision, dimension( n_slices, n_particles ), intent(inout) :: U

  ! Local variables
  integer :: i
  double precision, dimension( n_particles ) :: first_angular_derivatives

  do i = 1, n_slices
    call compute_angular_derivatives( &
      x(i,:,:), &
      rotation_plane_axis_1, rotation_plane_axis_2, &
      n_rotating_particles, &
      n_particles, n_dimensions, &
      first_angular_derivatives &
    )
    call accumulate_rotation_potential ( &
      x(i,:,:), lambda, &
      first_angular_derivatives, &
      rotation_plane_axis_1, rotation_plane_axis_2, frame_angular_velocity, &
      n_particles, n_dimensions, &
      U(i,:) &
    )
  end do

end subroutine
!@-node:gcross.20090827135657.1431:accumulate_effective_potential
!@+node:gcross.20090827135657.1433:accumulate_effective_potential2
pure subroutine accumulate_effective_potential2 (&
    x, lambda, &
    rotation_plane_axis_1, rotation_plane_axis_2, &
    frame_angular_velocity, n_rotating_particles, &
    n_slices, n_particles, n_dimensions, &
    U &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions, n_rotating_particles
  double precision, dimension ( n_slices, n_particles, n_dimensions ), intent(in) :: x
  double precision, intent(in) :: frame_angular_velocity, lambda
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2

  ! Output variables
  double precision, dimension( n_slices, n_particles ), intent(inout) :: U

  ! Local variables
  integer :: i
  double precision, dimension( n_particles ) :: first_angular_derivatives

  first_angular_derivatives = dble(n_rotating_particles)/dble(n_particles)
  do i = 1, n_slices
    call accumulate_rotation_potential ( &
      x(i,:,:), lambda, &
      first_angular_derivatives, &
      rotation_plane_axis_1, rotation_plane_axis_2, frame_angular_velocity, &
      n_particles, n_dimensions, &
      U(i,:) &
    )
  end do

end subroutine
!@-node:gcross.20090827135657.1433:accumulate_effective_potential2
!@+node:gcross.20090903090230.2068:compute_feynman_phase
pure function compute_feynman_phase( &
    x, n_rotating_particles, &
    rotation_axis_1, rotation_axis_2, &
    n_particles, n_dimensions &
  ) result (phase)
  integer, intent(in) :: rotation_axis_1, rotation_axis_2, n_rotating_particles
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_particles, n_dimensions), intent(in) :: x

  double precision :: phase

  phase = dble(n_rotating_particles)/n_particles & 
        * sum(atan2(x(:,rotation_axis_2),x(:,rotation_axis_1)))

end function
!@-node:gcross.20090903090230.2068:compute_feynman_phase
!@+node:gcross.20090903090230.2070:compute_feynman_backflow
pure function compute_feynman_backflow( &
    x, xij2, gradient_of_log_trial_fn, &
    rotation_axis_1, rotation_axis_2, &
    n_particles, n_dimensions &
  ) result (correction)
  integer, intent(in) :: rotation_axis_1, rotation_axis_2
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_particles, n_dimensions), intent(in) :: x, gradient_of_log_trial_fn
  double precision, dimension(n_particles, n_particles), intent(in) :: xij2

  double precision :: correction

  correction = sum(compute_term( &
    gradient_of_log_trial_fn(:,rotation_axis_1), & ! gx
    gradient_of_log_trial_fn(:,rotation_axis_2), & ! gy
    x(:,rotation_axis_1),x(:,rotation_axis_2) &
  ))

contains

  elemental function compute_term(gx,gy,x,y) result (term)
    double precision, intent(in) :: gx, gy, x, y
    double precision :: term

    double precision :: r_squared

    r_squared = x**2+y**2

    term = gx*y/r_squared - gy*x/r_squared
  end function

end function
!@-node:gcross.20090903090230.2070:compute_feynman_backflow
!@+node:gcross.20090903090230.2073:compute_feynman_phase_with_correction
pure function compute_feynman_phase_with_correction( &
    x, xij2, gradient_of_log_trial_fn, &
    n_rotating_particles, &
    rotation_axis_1, rotation_axis_2, &
    n_particles, n_dimensions &
  ) result (phase)
  integer, intent(in) :: rotation_axis_1, rotation_axis_2, n_rotating_particles
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension(n_particles, n_dimensions), intent(in) :: x, gradient_of_log_trial_fn
  double precision, dimension(n_particles, n_particles), intent(in) :: xij2

  double precision :: phase

  phase = compute_feynman_phase( &
            x, n_rotating_particles, &
            rotation_axis_1, rotation_axis_2, &
            n_particles, n_dimensions &
          ) &
        + compute_feynman_backflow( &
            x, xij2, gradient_of_log_trial_fn, &
            rotation_axis_1, rotation_axis_2, &
            n_particles, n_dimensions &
          )

end function
!@-node:gcross.20090903090230.2073:compute_feynman_phase_with_correction
!@+node:gcross.20090903090230.2072:accum_potential_via_numeric_diffs
subroutine accum_potential_via_numeric_diffs (&
    x, xij2, &
    frame_angular_velocity, lambda, &
    n_slices, n_particles, n_dimensions, &
    compute_phase, &
    U &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension ( n_slices, n_particles, n_dimensions ), intent(in) :: x
  double precision, dimension ( n_slices, n_particles, n_particles ), intent(in) :: xij2
  double precision, intent(in) :: frame_angular_velocity, lambda

!f2py external, intent(callback) :: compute_phase
  interface
    pure function compute_phase(x,xij2,n_particles,n_dimensions) result (phase)
      integer, intent(in) :: n_particles, n_dimensions
      double precision, dimension ( n_particles, n_dimensions ), intent(in) :: x
      double precision, dimension ( n_particles, n_particles ), intent(in) :: xij2
      double precision :: phase
    end function
  end interface

  ! Output variables
  double precision, dimension( n_slices, n_particles ), intent(inout) :: U

  ! Local variables
  integer :: i, j
  double precision, dimension( n_particles, n_dimensions ) :: phase_gradient
  double precision :: phase_laplacian

  do i = 1, n_slices
    call numerically_differentiate_jastrow( &
          x,xij2, &
          n_particles,n_dimensions, &
          compute_phase, &
          phase_gradient,phase_laplacian &
        )
    forall(j=1:n_particles) &
      U(i,j) = U(i,j) + & 
        sum(phase_gradient(j,:)*(lambda*phase_gradient(j,:)-frame_angular_velocity))
  end do

end subroutine
!@-node:gcross.20090903090230.2072:accum_potential_via_numeric_diffs
!@-others

end module angular_momentum
!@-node:gcross.20090807144330.2245:@thin angular_momentum.f95
!@-leo
