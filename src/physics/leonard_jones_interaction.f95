!@+leo-ver=4-thin
!@+node:gcross.20090901084550.2627:@thin leonard_jones_interaction.f95
!@@language fortran90

module leonard_jones_interaction
  implicit none

contains

!@+others
!@+node:gcross.20090901084550.2628:accumulate_potential
pure subroutine accumulate_potential( &
    xij2, &
    coefficient, cross_over_point_squared, &
    n_slices, n_particles, n_dimensions, &
    U &
  )
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension ( n_slices, n_particles, n_particles ), intent(in) :: xij2
  double precision, dimension( n_slices, n_particles ), intent(inout) :: U
  double precision, intent(in) :: coefficient, cross_over_point_squared

  integer :: i, j

  forall (i=1:n_slices,j=1:n_particles) &
    U(i,j) = U(i,j) + compute_U(xij2(i,:,j))

contains

  pure function compute_U(displacements) result ( U )
    double precision, dimension( n_particles ), intent(in) :: displacements
    double precision :: U, t3, t6
    t6 = cross_over_point_squared**6 * sum(displacements)**(-6)
    t3 = cross_over_point_squared**3 * sum(displacements)**(-3)
    U = coefficient * ( t6 - t3 )/2.0d0
  end function

end subroutine
!@-node:gcross.20090901084550.2628:accumulate_potential
!@-others

end module
!@-node:gcross.20090901084550.2627:@thin leonard_jones_interaction.f95
!@-leo
