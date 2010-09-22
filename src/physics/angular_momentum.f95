!@+leo-ver=4-thin
!@+node:gcross.20090807144330.2245:@thin angular_momentum.f95
!@@language fortran90
module angular_momentum

 use numeric_differentiation
 use gfn

 implicit none

 !@ << Constants >>
 !@+node:gcross.20090807144330.2252:<< Constants >>
 integer, parameter :: x_axis_label = 1, y_axis_label = 2, z_axis_label = 3
 !@-node:gcross.20090807144330.2252:<< Constants >>
 !@nl

 !@ << Interfaces >>
 !@+node:gcross.20090916153857.1668:<< Interfaces >>
 interface
   !@  @+others
   !@+node:gcross.20090916153857.1669:dnrm2
   pure function dnrm2(n,x,incx) result(norm)
     integer, intent(in) :: n, incx
     double precision, dimension(n), intent(in) :: x
     double precision :: norm
   end function
   !@-node:gcross.20090916153857.1669:dnrm2
   !@+node:gcross.20090916153857.1670:dgels
   pure subroutine dgelss(m,n,nrhs,a,lda,b,ldb,s,rcond,rank,work,lwork,info)
     integer, intent(in) :: m,n,nrhs,lda,ldb,lwork
     integer, intent(out) :: info, rank
     double precision, dimension(lda,n), intent(inout) :: a
     double precision, dimension(ldb,nrhs), intent(inout) :: b
     double precision, intent(inout) :: work
     double precision, intent(in) :: rcond
     double precision, dimension(min(m,n)), intent(out) :: s
   end subroutine
   !@-node:gcross.20090916153857.1670:dgels
   !@-others
 end interface
 !@-node:gcross.20090916153857.1668:<< Interfaces >>
 !@nl

 contains

!@+others
!@+node:gcross.20100117204224.1699:compute_sum_and_its_derivatives
subroutine compute_sum_and_its_derivatives( &
  number_of_amplitudes, amplitudes, &
  number_of_amplitudes_to_include, &
  sum_over_symmetrizations, gradient_of_sum &
  )
  integer, intent(in) :: number_of_amplitudes, number_of_amplitudes_to_include
  double complex, intent(in) :: amplitudes(number_of_amplitudes)
  double complex, intent(out) :: sum_over_symmetrizations, gradient_of_sum(number_of_amplitudes)

  integer :: m

  gradient_of_sum = (0d0,0d0)
  sum_over_symmetrizations = (1d0,0d0)

  do m = 1, number_of_amplitudes_to_include
    gradient_of_sum(:) = sum_over_symmetrizations - amplitudes(:)*gradient_of_sum(:)
    sum_over_symmetrizations = sum(amplitudes(:)*gradient_of_sum(:))/m
  end do
end subroutine
!@-node:gcross.20100117204224.1699:compute_sum_and_its_derivatives
!@+node:gcross.20100113111641.1972:compute_gradient_ho_phase_power
subroutine compute_gradient_ho_phase_power (&
    x, &
    n_total_rotating_particles, &
    n_slices, n_particles, n_dimensions, &
    gradient_phase &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension(n_dimensions,n_particles,n_slices), intent(in) :: x
  integer, intent(in) :: n_total_rotating_particles

  ! Output variables
  double precision, intent(out) :: gradient_phase(n_dimensions,n_particles,n_slices)

  ! Local variables
  integer :: slice
  double complex :: &
    full_sum, full_sum_abs_sq, partial_sum, &
    amplitudes(n_particles), &
    x_sum, y_sum, r_sq_sum

  gradient_phase = 0

  if (n_total_rotating_particles == 0) then
    return
  end if

  do slice = 1, n_slices
    x_sum = sum(x(1,:,slice))
    y_sum = sum(x(2,:,slice))
    r_sq_sum = x_sum**2 + y_sum**2
    gradient_phase(1,:,slice) = n_total_rotating_particles * (-y_sum / r_sq_sum)
    gradient_phase(2,:,slice) = n_total_rotating_particles * ( x_sum / r_sq_sum)
  end do

end subroutine
!@nonl
!@-node:gcross.20100113111641.1972:compute_gradient_ho_phase_power
!@+node:gcross.20091210143551.1691:compute_gradient_ho_phase_choice
subroutine compute_gradient_ho_phase_choice (&
    x, &
    n_total_rotating_particles, &
    n_slices, n_particles, n_dimensions, &
    gradient_phase &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension(n_slices,n_particles,n_dimensions), intent(in) :: x
  integer, intent(in) :: n_total_rotating_particles

  ! Output variables
  double precision, intent(out) :: gradient_phase(n_slices,n_particles,n_dimensions)

  ! Local variables
  integer :: particle, slice, n_excess_rotating_particles, feynman_factor
  double complex :: &
    full_sum, full_sum_abs_sq, &
    partial_sums(n_particles), amplitudes(n_particles)
  double precision :: &
    x_sum, y_sum, r_sq_sum, &
    factor


  if (n_total_rotating_particles >= n_particles) then
    feynman_factor = n_total_rotating_particles / n_particles
    do slice = 1, n_slices
      do particle = 1, n_particles
        factor = feynman_factor / (x(1,particle,slice)**2+x(2,particle,slice)**2)
        gradient_phase(1,particle,slice) = -x(2,particle,slice) * factor
        gradient_phase(2,particle,slice) =  x(1,particle,slice) * factor
      end do
    end do
    gradient_phase(3:,:,:) = 0
  else
    gradient_phase = 0
  end if

  n_excess_rotating_particles = mod(n_total_rotating_particles,n_particles)

  if (n_excess_rotating_particles == 0) then
    return
  end if

  if (n_excess_rotating_particles == 1) then
    do slice = 1, n_slices
      x_sum = sum(x(1,:,slice))
      y_sum = sum(x(2,:,slice))
      r_sq_sum = x_sum**2 + y_sum**2
      gradient_phase(1,:,slice) = gradient_phase(1,:,slice) - y_sum / r_sq_sum
      gradient_phase(2,:,slice) = gradient_phase(2,:,slice) + x_sum / r_sq_sum
    end do
    return
  end if

  do slice = 1, n_slices
    amplitudes = x(2,:,slice)*(1,0)-x(1,:,slice)*(0,1)
    call compute_sum_and_its_derivatives( &
      n_particles,amplitudes, &
      n_excess_rotating_particles, &
      full_sum, partial_sums &
    )
    full_sum_abs_sq = real(full_sum * conjg(full_sum))
    gradient_phase(1,:,slice) = gradient_phase(1,:,slice) - &
      (imag(full_sum)*imag(partial_sums(:)) + real(full_sum)*real(partial_sums(:)))/full_sum_abs_sq
    gradient_phase(2,:,slice) = gradient_phase(2,:,slice) - &
      (imag(full_sum)*real(partial_sums(:)) - real(full_sum)*imag(partial_sums(:)))/full_sum_abs_sq
  end do

end subroutine
!@nonl
!@-node:gcross.20091210143551.1691:compute_gradient_ho_phase_choice
!@+node:gcross.20090903090230.2072:accumulate_rotating_frame_potential
pure subroutine accumulate_rotating_frame_potential (&
    x, gradient_phase, &
    frame_angular_velocity, lambda, &
    rotation_plane_axis_1, rotation_plane_axis_2, &
    n_slices, n_particles, n_dimensions, &
    U &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2
  double precision, dimension (n_dimensions,n_particles,n_slices), intent(in) :: x, gradient_phase
  double precision, intent(in) :: frame_angular_velocity, lambda

  ! Output variables
  double precision, dimension(n_particles,n_slices), intent(inout) :: U

  ! Local variables
  integer :: slice, particle

  do slice = 1, n_slices
    do particle = 1, n_particles
      U(particle,slice) = U(particle,slice) + lambda*sum(gradient_phase(:,particle,slice)**2) &
                         - frame_angular_velocity*( &
                             x(rotation_plane_axis_1,particle,slice)*gradient_phase(rotation_plane_axis_2,particle,slice) &
                           - x(rotation_plane_axis_2,particle,slice)*gradient_phase(rotation_plane_axis_1,particle,slice) &
                           )
    end do
  end do

end subroutine
!@nonl
!@-node:gcross.20090903090230.2072:accumulate_rotating_frame_potential
!@+node:gcross.20090919132620.2307:accumulate_effective_potential
pure subroutine accumulate_effective_potential (&
    gradient_phase, &
    lambda, &
    n_slices, n_particles, n_dimensions, &
    U &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension (n_dimensions,n_particles,n_slices), intent(in) :: gradient_phase
  double precision, intent(in) :: lambda

  ! Output variables
  double precision, dimension(n_particles,n_slices), intent(inout) :: U

  ! Local variables
  integer :: slice, particle

  do slice = 1, n_slices
    do particle = 1, n_particles
      U(particle,slice) = U(particle,slice) + lambda*sum(gradient_phase(:,particle,slice)**2)
    end do
  end do

end subroutine
!@nonl
!@-node:gcross.20090919132620.2307:accumulate_effective_potential
!@+node:gcross.20090919132620.2310:accumulate_magnetic_field_phase
pure subroutine accumulate_magnetic_field_phase (&
    x, &
    magnetic_field_strength, &
    rotation_plane_axis_1, rotation_plane_axis_2, &
    n_slices, n_particles, n_dimensions, &
    gradient_phase &
  )

  ! Input variables
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension(n_dimensions,n_particles,n_slices), intent(in) :: x
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2
  double precision, intent(in) :: magnetic_field_strength

  ! Output variables
  double precision, dimension (n_dimensions,n_particles,n_slices), intent(inout) :: gradient_phase

  ! Local variables
  integer :: slice, particle

  do slice = 1, n_slices
    do particle = 1, n_particles
      gradient_phase(rotation_plane_axis_1,particle,slice) = gradient_phase(rotation_plane_axis_1,particle,slice) &
        - magnetic_field_strength * x(rotation_plane_axis_2,particle,slice)
      gradient_phase(rotation_plane_axis_2,particle,slice) = gradient_phase(rotation_plane_axis_2,particle,slice) &
        + magnetic_field_strength * x(rotation_plane_axis_1,particle,slice)
    end do
  end do

end subroutine
!@nonl
!@-node:gcross.20090919132620.2310:accumulate_magnetic_field_phase
!@-others

end module angular_momentum
!@-node:gcross.20090807144330.2245:@thin angular_momentum.f95
!@-leo
