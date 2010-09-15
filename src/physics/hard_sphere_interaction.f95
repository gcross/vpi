!@+leo-ver=4-thin
!@+node:gcross.20090901084550.2616:@thin hard_sphere_interaction.f95
!@@language fortran90

module hard_sphere_interaction

  use gfn

  implicit none

contains

!@+others
!@+node:gcross.20090901084550.2617:compute_trial_weight
subroutine compute_trial_weight( &
    xij2, &
    hard_sphere_radius, &
    n_particles, &
    weight, &
    reject_flag &
  )
  integer, intent(in) ::  n_particles
  double precision, dimension( n_particles, n_particles ), intent(in) :: xij2
  double precision, intent(in) :: hard_sphere_radius
  double precision, intent(out) :: weight
  logical, intent(out) :: reject_flag
  double precision  :: r2, hard_sphere_radius_squared

  integer :: i, j

  hard_sphere_radius_squared = hard_sphere_radius ** 2

  weight = 0d0
  reject_flag = .false.
  do i = 1, n_particles
    do j = i + 1, n_particles
      r2 = xij2(i,j)
      if ( r2 > hard_sphere_radius_squared ) then
        weight = weight + log(1.0d0 - hard_sphere_radius/sqrt(r2))
      else
        reject_flag = .true.
        return
      end if
    end do
  end do

end subroutine
!@-node:gcross.20090901084550.2617:compute_trial_weight
!@+node:gcross.20090901084550.2622:accumulate_trial_derivatives
pure subroutine accumulate_trial_derivatives( &
    x, xij2, &
    hard_sphere_radius, &
    n_particles, n_dimensions, &
    grad_lntfn, lap_lntfn &
  )
  integer, intent(in) :: n_particles, n_dimensions
  double precision, dimension( n_particles, n_dimensions ), intent(in) :: x
  double precision, dimension( n_particles, n_particles ), intent(in) :: xij2
  double precision, intent(in) :: hard_sphere_radius
  double precision, dimension( n_particles, n_dimensions ), intent(inout) :: grad_lntfn 
  double precision, intent(inout) :: lap_lntfn 

  double precision :: fi,fi2,ri,ri2,ri3,ri4
  double precision, dimension( n_dimensions ) :: gtmp
  double precision :: hard_sphere_radius_squared
  integer :: i, j

  hard_sphere_radius_squared = hard_sphere_radius ** 2

  do i = 1, n_particles
    do j = i + 1, n_particles
      ri2 = 1.0d0/xij2(i,j)
      ri = sqrt(ri2)
      ri3 = ri*ri2
      ri4 = ri2*ri2
      fi = 1.0d0/(1.0d0 - hard_sphere_radius*ri)
      fi2 = fi*fi
      gtmp(:) = ri3*fi*(x(i,:) - x(j,:))
      grad_lntfn(i,:) = grad_lntfn(i,:) + hard_sphere_radius*gtmp(:)
      grad_lntfn(j,:) = grad_lntfn(j,:) - hard_sphere_radius*gtmp(:)
      lap_lntfn = lap_lntfn - 2.0d0*hard_sphere_radius_squared*fi2*ri4
    end do
  end do
end subroutine
!@-node:gcross.20090901084550.2622:accumulate_trial_derivatives
!@+node:gcross.20090908085435.1631:compute_gradient_backflow
pure subroutine compute_gradient_backflow( &
    x, xij2, &
    hard_sphere_radius, rotation_rate, &
    rotation_plane_axis_1, rotation_plane_axis_2, &
    n_slices, n_particles, n_dimensions, &
    gradient_backflow &
  )
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, dimension( n_slices, n_particles, n_dimensions ), intent(in) :: x
  double precision, dimension( n_slices, n_particles, n_particles ), intent(in) :: xij2
  double precision, intent(in) :: hard_sphere_radius, rotation_rate
  integer, intent(in) :: rotation_plane_axis_1, rotation_plane_axis_2
  double precision, dimension( n_slices, n_particles, n_dimensions ), intent(out) :: gradient_backflow

  double precision :: recipricol_rho_ip_squared, recipricol_rho_ip_4th, r_ip_j_squared, r_ip_j, recipricol_rho_j_squared
  double precision :: difference_recipical_rho_squared, common_factor, C_ip_j, r_ip_j_cubed, D_ip_j
  double precision, dimension( n_dimensions ) :: term
  integer :: s, ip, j, k1, k2

  ! Rename parameters to match notes
  k1 = rotation_plane_axis_1
  k2 = rotation_plane_axis_2

  gradient_backflow = 0d0

  do ip = 1, n_particles
    do s = 1, n_slices
      recipricol_rho_ip_squared = 1d0/(x(s,ip,k1)**2 + x(s,ip,k2)**2)
      recipricol_rho_ip_4th = recipricol_rho_ip_squared**2
      do j = 1, n_particles
        if (j == ip) then
          cycle
        end if
        r_ip_j_squared = xij2(s,ip,j)
        r_ip_j = sqrt(r_ip_j_squared)
        recipricol_rho_j_squared = 1d0/(x(s,j,k1)**2 + x(s,j,k2)**2)
        difference_recipical_rho_squared = recipricol_rho_ip_squared - recipricol_rho_j_squared
        common_factor = (1.0/r_ip_j_squared)*difference_recipical_rho_squared &
                            *(3d0+hard_sphere_radius/(r_ip_j*(1-hard_sphere_radius/r_ip_j)))
        term(:) = common_factor * (x(s,ip,:)-x(s,j,:))

        term(k1) = term(k1) + 2*x(s,ip,k1)*recipricol_rho_ip_4th
        term(k2) = term(k2) + 2*x(s,ip,k2)*recipricol_rho_ip_4th

        C_ip_j = x(s,ip,k1)*x(s,j,k2)-x(s,ip,k2)*x(s,j,k1)
        term(:) = (-C_ip_j)*term(:)

        term(k1) = term(k1) + x(s,j,k2)*difference_recipical_rho_squared
        term(k2) = term(k2) - x(s,j,k1)*difference_recipical_rho_squared

        r_ip_j_cubed = r_ip_j * r_ip_j_squared
        D_ip_j = r_ip_j_cubed*(1d0-hard_sphere_radius/r_ip_j)    
        term(:) = hard_sphere_radius*rotation_rate/D_ip_j*term(:)

        gradient_backflow(s,ip,:) = gradient_backflow(s,ip,:) + term(:)
      end do
    end do
  end do
end subroutine
!@nonl
!@-node:gcross.20090908085435.1631:compute_gradient_backflow
!@+node:gcross.20090828201103.2127:has_collision
pure function has_collision(xij2,hard_sphere_radius_squared,n_slices,n_particles)
  integer, intent(in) :: n_slices, n_particles
  double precision, intent(in) :: hard_sphere_radius_squared
  double precision, dimension(n_slices,n_particles,n_particles), intent(in) :: xij2
  logical :: has_collision

  integer :: i, j

  has_collision = .false.
  do i = 1, n_particles
    do j = i+1, n_particles
      if (any(xij2(:,i,j) <= hard_sphere_radius_squared)) then
        has_collision = .true.
        return
      end if
    end do
  end do
end function
!@nonl
!@-node:gcross.20090828201103.2127:has_collision
!@+node:gcross.20090902085220.2790:Green's functions
!@+node:gcross.20090828095451.1455:strategy 1 (image approximation)
! image approximation
pure function compute_greens_function( &
    xij2, &
    dt, hard_sphere_radius, &
    slice_start, slice_end, &
    particle_number, &
    n_slices, n_particles &
  ) result ( ln_gfn )
  integer, intent(in) :: slice_start, slice_end, particle_number
  integer, intent(in) :: n_slices, n_particles
  double precision, intent(in) :: dt, hard_sphere_radius
  double precision, dimension ( n_slices, n_particles, n_particles ), intent(in) :: xij2

  double precision :: gfn, ln_gfn
  integer :: i

  gfn = 1.0d0

  do i = 1, particle_number-1
    gfn = gfn * compute_contribution_from_particle(i)
  end do

  do i = particle_number+1, n_particles
    gfn = gfn * compute_contribution_from_particle(i)
  end do

  ln_gfn = log(gfn)

contains

  pure function compute_contribution_from_particle(i) result (gfn)
    integer, intent(in) :: i
    double precision :: gfn

    double precision, dimension(n_slices) :: distances

    distances(slice_start:slice_end) = sqrt(xij2(slice_start:slice_end,particle_number,i)) - hard_sphere_radius

    gfn = compute_green_fn_from_distances( &
            distances, &
            dt, &
            slice_start, slice_end, &
            n_slices &
          )
  end function

end function
!@nonl
!@-node:gcross.20090828095451.1455:strategy 1 (image approximation)
!@+node:gcross.20090828095451.1456:strategy 2
! A new quantum propagator for hard sphere and cavity systems
!J. Cao ad B.J. Berne
!JCP 97 (4) 2382  (1992)
pure function compute_greens_function2( &
    q, &
    lambda, dt, hard_sphere_radius, &
    slice_start, slice_end, &
    particle_number, &
    n_slices, n_particles, n_dimensions &
  ) result ( ln_gfn )
  integer, intent(in) :: slice_start, slice_end, particle_number
  integer, intent(in) :: n_slices, n_particles, n_dimensions
  double precision, intent(in) :: dt, lambda, hard_sphere_radius
  double precision, dimension ( n_slices, n_particles, n_dimensions ), intent(in) :: q

  double precision :: hs_gfn, ln_gfn

  double precision :: tmp

  double precision, dimension ( n_dimensions ) :: dx,dxp,dx12
  double precision :: u,r,rp,D,r12
  double precision :: hard_sphere_radius_squared
  integer :: i,j
  integer :: center_slice_number

  center_slice_number = n_slices / 2
  hard_sphere_radius_squared = hard_sphere_radius**2

  hs_gfn = 1.0d0

  do j = slice_start+1, center_slice_number, 2
    do i = 1, n_particles
      if(i .ne. particle_number) then
        dx(:) = q(j-1,particle_number,:) - q(j-1,i,:)
        r = sqrt(sum(dx(:)**2))
        dxp(:) = q(j,particle_number,:) - q(j,i,:)
        rp = sqrt(sum(dxp(:)**2))
        dx12(:) = dx(:) - dxp(:)
        r12 = sum(dx12(:)**2)
        u = r+rp
        tmp = hard_sphere_radius*(u-hard_sphere_radius)/(r*rp)
        D = u**2 - 2.0d0*(hard_sphere_radius*u - hard_sphere_radius_squared + tmp*dot_product(dx,dxp))
        hs_gfn = hs_gfn * (1.0d0 - tmp*exp(-(D-r12)/(2.0d0*lambda*dt)))
      end if
    end do
  end do

  do j = center_slice_number+2, slice_end, 2
    do i = 1, n_particles
      if(i .ne. particle_number) then
        dx(:) = q(j-1,particle_number,:) - q(j-1,i,:)
        r = sqrt(sum(dx(:)**2))
        dxp(:) = q(j,particle_number,:) - q(j,i,:)
        rp = sqrt(sum(dxp(:)**2))
        dx12(:) = dx(:) - dxp(:)
        r12 = sum(dx12(:)**2)
        u = r+rp
        tmp = hard_sphere_radius*(u-hard_sphere_radius)/(r*rp)
        D = u**2 - 2.0d0*(hard_sphere_radius*u - hard_sphere_radius_squared + tmp*dot_product(dx,dxp))
        hs_gfn = hs_gfn * (1.0d0 - tmp*exp(-(D-r12)/(2.0d0*lambda*dt)))
      end if
    end do
  end do

  ln_gfn = log(hs_gfn)

end function
!@nonl
!@-node:gcross.20090828095451.1456:strategy 2
!@-node:gcross.20090902085220.2790:Green's functions
!@-others

end module
!@-node:gcross.20090901084550.2616:@thin hard_sphere_interaction.f95
!@-leo
