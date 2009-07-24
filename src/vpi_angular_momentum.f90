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
!@+node:gcross.20090721121051.1756:compute_effective_rotational_potential
subroutine compute_effective_rotational_potential (&
! INPUT: particle position information
    x, &
! INPUT: path slice to consider
    move_start, move_end, &
! INPUT: 
    nslice, np, ndim, &
! OUTPUT: effective rotational potential
    U_rot &
  )
!@+at
! Array dimensions / slicing
!@-at
!@@c
  integer, intent(in) :: move_start, move_end, nslice, np, ndim
!@+at
! Function input
!@-at
!@@c
  real(kind=b8), dimension ( nslice, np , ndim ), intent(in) :: x
!@+at
! Function output
!@-at
!@@c
  real(kind=b8), dimension( nslice, np ), intent(out) :: U_rot
!@+at
! Temporary variables
!@-at
!@@c
  integer :: plane_axis_1, plane_axis_2
!@+at
! Code begins
!@-at
!@@c
  call get_rotation_plane_axes(fixed_rotation_axis,plane_axis_1,plane_axis_2)
  U_rot(move_start:move_end,:) = abs(float(fixed_angular_momentum)) / &
                                    ( x(move_start:move_end,:,plane_axis_1)**2 & 
                                     +x(move_start:move_end,:,plane_axis_2)**2 )
end subroutine compute_effective_rotational_potential
!@-node:gcross.20090721121051.1756:compute_effective_rotational_potential
!@+node:gcross.20090723093414.1759:compute_angular_interference
!@+at
! This function computers the effect of interference due to the effect that 
! some particles have angular momentum and others do not, but which particles 
! have this momentum is non-deterministic due to boson symmetry.  To see why 
! this matters, consider two bosons in a periodic box, with one excited and 
! the other in its ground state.  The (unnormalized) amplitude of the wave 
! function is then $e^{x_1} + e^{x_2}$.
! 
! To generalize this, if there are m particles with one unit of angular 
! momentum in our system, then we need to sum over all terms 
! $e^{x_{i_1}+x_{i_2}+\dots+x_{i_m}}$ for all $i_1\ne i_2\ne \dots \ne i_m$.
! 
! This function employs a trick to carry out this.
!@-at
!@@c
function compute_angular_interference_amp (&
! INPUT: particle position information
    x, &
! INPUT: path slice to consider
    move_start, move_end, &
! INPUT: array dimensions
    nslice, np, ndim &
! OUTPUT: complex amplitude resulting from the interference
  ) result ( weight )
!@+at
! Array dimensions / slicing
!@-at
!@@c
  integer, intent(in) :: move_start, move_end, nslice, np, ndim
!@+at
! Function input
!@-at
!@@c
  real(kind=b8), dimension ( nslice, np , ndim ), intent(in), target :: x
!@+at
! Function output
!@-at
!@@c
  real(kind=b8) :: weight
!@+at
! Temporary variables
!@-at
!@@c
  integer :: plane_axis_1, plane_axis_2
  complex(kind=b8), dimension( move_end-move_start+1, np ) :: single_particle_amplitudes
  complex(kind=b8), dimension( move_end-move_start+1, np-fixed_angular_momentum+1 ) :: interference_vector
  real(kind=b8), dimension(:,:), pointer :: x_slice_1, x_slice_2
  complex(kind=b8) :: amplitude
  integer :: i, j
!@+at
! Code begins
!@-at
!@@c 
  if( fixed_angular_momentum == 0 .or. &
      move_end == move_start .or. &
      size(interference_vector,1) == 0 &
      ) then
    weight = 1.0_b8
    return
  end if

  call get_rotation_plane_axes(fixed_rotation_axis,plane_axis_1,plane_axis_2)

  x_slice_1 => x(move_start:move_end,:,plane_axis_1)
  x_slice_2 => x(move_start:move_end,:,plane_axis_2)

  single_particle_amplitudes(:,:) = &
    (x_slice_1(:,:)*(1.0_b8,0) + x_slice_2(:,:)*(0,1.0_b8))/sqrt(x_slice_1(:,:)**2 + x_slice_2(:,:)**2)

  interference_vector = single_particle_amplitudes(:,:size(interference_vector,2))

  do i = 1, fixed_angular_momentum-1
    do j = size(interference_vector,2), 1, -1
      interference_vector(:,j) = sum(interference_vector(:,j:),dim=2)*single_particle_amplitudes(:,i+j)
    end do
  end do

  weight = 1_b8
  do i = 1, size(interference_vector,1)  
    amplitude = sum(interference_vector(i,:))/float(np)
    weight = weight + (real(amplitude)**2 + imag(amplitude)**2)
  end do

end function
!@-node:gcross.20090723093414.1759:compute_angular_interference
!@-others

end module vpi_angular_momentum
!@-node:gcross.20090724102123.1771:@thin vpi_angular_momentum.f90
!@-leo
