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
function sum_over_symmetrizations(amplitudes,N_particles,N_rotating_particles)
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

end function sum_over_symmetrizations
!@-node:gcross.20090803153449.1836:sum_over_symmetrizations
!@+node:gcross.20090803153449.1837:compute_angular_derivatives
subroutine compute_angular_derivatives(&
    x, &
    rotation_plane_axis_1, rotation_plane_axis_2, N_rotating_particles, &
    N_particles,N_dimensions, &
    first_derivatives, second_derivatives &
  )

  ! Input variables
  integer, intent(in) :: N_particles, N_dimensions
  double precision, dimension(N_particles,N_dimensions), intent(in) :: x
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2, N_rotating_particles

  ! Output variables
  double precision, dimension(N_particles), intent(out) :: first_derivatives
  double precision, dimension(N_particles,N_particles), intent(out) :: second_derivatives

  ! Local variables
  double complex, dimension(N_particles) :: amplitudes, one_index_restricted_sums
  double complex, dimension(N_particles,N_particles) :: two_index_restricted_sums
  double complex :: full_sum, full_sum_conj, ai, aj
  double precision :: full_sum_amplitude_squared
  integer :: i, j

  if(N_rotating_particles == 0) then
    first_derivatives = 0
    second_derivatives = 0
    return
  end if

  if(N_particles == N_rotating_particles) then
    first_derivatives = 1
    second_derivatives = 0
    return
  end if

  amplitudes(:) = &
    (x(:,rotation_plane_axis_1)*(1.0d0,0) + x(:,rotation_plane_axis_2)*(0,1.0d0)) &
  / sqrt(x(:,rotation_plane_axis_1)**2 + x(:,rotation_plane_axis_2)**2)

  if(N_rotating_particles == 1) then
    full_sum = sum(amplitudes)
    full_sum_conj = conjg(full_sum)
    full_sum_amplitude_squared = real(full_sum*conjg(full_sum))
    first_derivatives(:) = real(amplitudes(:)*full_sum_conj)/full_sum_amplitude_squared
    do i = 1, n_particles
      do j = 1, n_particles
        if(i /= j) then
          second_derivatives(i,j) = &
                imag(amplitudes(i)*conjg(amplitudes(j)))/full_sum_amplitude_squared &
            + 2*real(amplitudes(i)*full_sum)*imag(amplitudes(j)*full_sum)/full_sum_amplitude_squared**2
        else
          second_derivatives(i,j) = &
                imag(amplitudes(i)*full_sum)/full_sum_amplitude_squared &
            + 2*real(amplitudes(i)*full_sum)*imag(amplitudes(j)*full_sum)/full_sum_amplitude_squared**2
        end if
      end do
    end do
    return
  end if

! Compute sums which require the appearance of two specific indicies
  do i = 1, N_particles
    do j = i+1, N_particles
      ai = amplitudes(i)
      aj = amplitudes(j)
      amplitudes(i) = (0d0,0d0)
      amplitudes(j) = (0d0,0d0)
      two_index_restricted_sums(i,j) = sum_over_symmetrizations(amplitudes,N_particles,N_rotating_particles-2)*ai*aj
      two_index_restricted_sums(j,i) = two_index_restricted_sums(i,j)
      amplitudes(i) = ai
      amplitudes(j) = aj
    end do
  end do

! Compute sums which require the appearance of one specific index
  do i = 1, N_particles
    ai = amplitudes(i)
    amplitudes(i) = (0d0,0d0)
    one_index_restricted_sums(i) = sum_over_symmetrizations(amplitudes,N_particles,N_rotating_particles-1)*ai
    two_index_restricted_sums(i,i) = one_index_restricted_sums(i)
    amplitudes(i) = ai
  end do

! Compute the full sum
  full_sum = sum_over_symmetrizations(amplitudes,N_particles,N_rotating_particles)

! Compute first and second derivatives
  full_sum_amplitude_squared = real(full_sum*conjg(full_sum))
  full_sum_conj = conjg(full_sum)
  first_derivatives(:) = real(one_index_restricted_sums(:)*full_sum_conj)/full_sum_amplitude_squared

  do i = 1, N_particles
    do j = 1, N_particles
      second_derivatives(i,j) = &
        ( &
            imag(one_index_restricted_sums(i)*conjg(one_index_restricted_sums(j))) &
          - imag(two_index_restricted_sums(i,j)*full_sum_conj) &
        ) / full_sum_amplitude_squared &
      +2* ( &
            real(one_index_restricted_sums(i)*full_sum_conj) &
          * imag(one_index_restricted_sums(j)*full_sum_conj) &
        ) / full_sum_amplitude_squared**2
    end do
  end do

end subroutine compute_angular_derivatives
!@-node:gcross.20090803153449.1837:compute_angular_derivatives
!@+node:gcross.20090721121051.1756:compute_effective_rotational_potential
subroutine compute_effective_rotational_potential (&
    x, lambda, &
    rotation_plane_axis_1, rotation_plane_axis_2, frame_angular_velocity, N_rotating_particles, &
    n_slices, n_particles, n_dimensions, &
    U, gradU &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension ( n_slices, n_particles, n_dimensions ), intent(in) :: x
  double precision, intent(in) :: frame_angular_velocity, lambda
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2, N_rotating_particles

  ! Output variables
  double precision, dimension( n_slices, n_particles ), intent(inout) :: U
  double precision, dimension( n_slices, n_particles, n_dimensions ), intent(inout) :: gradU
  double precision, dimension( n_particles ) :: first_derivatives
  double precision, dimension( n_particles, n_particles ) :: second_derivatives

  ! Local variables
  integer :: slice, i, j
  double precision :: rho_i_squared, rho_j, rho_j_squared, C, S, F

  do slice = 1, n_slices
    call compute_angular_derivatives( &
            x(slice,:,:), &
            rotation_plane_axis_1, rotation_plane_axis_2, N_rotating_particles, &
            n_particles, n_dimensions, &
            first_derivatives, second_derivatives &
          )
    do j = 1, n_particles
      rho_j_squared = x(slice,j,rotation_plane_axis_1)**2 + x(slice,j,rotation_plane_axis_2)**2
      U(slice,j) = U(slice,j) + first_derivatives(j)*(first_derivatives(j)*lambda/rho_j_squared - frame_angular_velocity)
      rho_j = sqrt(rho_j_squared)
      C = x(slice,j,rotation_plane_axis_1)/rho_j
      S = x(slice,j,rotation_plane_axis_2)/rho_j
      do i = 1, n_particles
        rho_i_squared = x(slice,i,rotation_plane_axis_1)**2 + x(slice,i,rotation_plane_axis_2)**2
        F = second_derivatives(i,j) * &
          ( &
            first_derivatives(i) * 2d0*lambda/(rho_j*rho_i_squared) &
          - frame_angular_velocity / rho_j &
          )
        gradU(slice,j,rotation_plane_axis_1) = &
          gradU(slice,j,rotation_plane_axis_1) - S*F
        gradU(slice,j,rotation_plane_axis_2) = &
          gradU(slice,j,rotation_plane_axis_2) + C*F
      end do
      F = first_derivatives(j)**2 * 2d0*lambda/(rho_j*rho_j_squared)
      gradU(slice,j,rotation_plane_axis_1) = &
        gradU(slice,j,rotation_plane_axis_1) - C*F
      gradU(slice,j,rotation_plane_axis_2) = &
        gradU(slice,j,rotation_plane_axis_2) - S*F
    end do
  end do

end subroutine compute_effective_rotational_potential
!@-node:gcross.20090721121051.1756:compute_effective_rotational_potential
!@-others

end module angular_momentum
!@-node:gcross.20090807144330.2245:@thin angular_momentum.f95
!@-leo
