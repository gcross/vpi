module vpi_phase_maps
  use vpi_defines
  use vpi_coordinate_utils
  implicit none

contains 

subroutine vortex( phase, sl_start, sl_end, x, xij2 )
  real(kind=b8), dimension( N_SLICE ), intent(inout) :: phase
  integer :: sl_start, sl_end
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_PARTICLE ) :: xij2

  real(kind=b8) :: qx, qy
  real(kind=b8) :: tph

  integer :: i,j

  do i = sl_start, sl_end
    tph = 0
    do j = 1, N_PARTICLE
      qx = x(i,j,1)
      qy = x(i,j,2)
      tph = tph + eval_polar_angle(qx,qy)
    end do
    phase(i) = cos(tph)
  end do

end subroutine vortex


subroutine nvortex( phase, sl_start, sl_end, x, xij2 )
  real(kind=b8), dimension( N_SLICE ), intent(inout) :: phase
  integer :: sl_start, sl_end
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_PARTICLE ) :: xij2

  real(kind=b8) :: cx, cy
  real(kind=b8) :: qx, qy
  real(kind=b8) :: tph

  integer :: i,j,k

  do i = sl_start, sl_end
    tph = 0
    do j = 1, N_PARTICLE
      cx = x(i,j,1)
      cy = x(i,j,2)
      do k=1, N_VORTEX
        qx = cx + p_vcoords(k,1)
        qy = cy + p_vcoords(k,2)
        tph = tph + eval_polar_angle(qx,qy)
      end do
    end do
    phase(i) = cos(tph)
  end do

end subroutine nvortex

! simple superposition
subroutine vortex_sup( phase, sl_start, sl_end, x, xij2 )
  real(kind=b8), dimension( N_SLICE ), intent(inout) :: phase
  integer :: sl_start, sl_end
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_PARTICLE ) :: xij2

  real(kind=b8) :: qx, qy, qr
  real(kind=b8) :: angle
  complex(kind=b8) :: tph

  integer :: i,j

  do i = sl_start, sl_end
    tph = (1.0_b8,0)
    do j = 1, N_PARTICLE
      qx = x(i,j,1)
      qy = x(i,j,3)
      angle = eval_polar_angle(qx,qy)
      tph = tph*cmplx(1+cos(angle),sin(angle))
    end do
    phase(i) = atan(aimag(tph)/real(tph))
!    write (5000,*) phase(i)
  end do

end subroutine vortex_sup


end module vpi_phase_maps
