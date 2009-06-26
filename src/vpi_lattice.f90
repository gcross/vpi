!@+leo-ver=4-thin
!@+node:gcross.20090624144408.2046:@thin vpi_lattice.f90
!@@language fortran90

module vpi_lattice

  use vpi_defines
  use vpi_rand_utils

  implicit none

  integer, parameter :: LATTICE_FILL = 0
  integer, parameter :: GAUSSIAN_FILL = 1
  integer, parameter :: DWELL_GAUSSIAN_FILL = 2


contains 

subroutine vpi_make_lattice( q, a, fill_method) 
  real(kind=b8), dimension( :, :, : ), intent(inout) :: q
  real :: a

  integer :: nslice, np, ndim, fill_method
  integer :: i,j,k
  real :: x,y,z,r2
  real(kind=b8), dimension( 3 ) :: nu

  nslice = size(q,1)
  np = size(q,2)
  ndim = size(q,3)

  select case (fill_method)
   case (LATTICE_FILL)

    i = 1
    do while ( i .le. np )
10  call random_number( nu )
      q(1,i,1) = (nu(1)-0.5)*a
      q(1,i,2) = (nu(2)-0.5)*a
      q(1,i,3) = (nu(3)-0.5)*a

      do j = 1, i-1
        x = get_pbc(q(1,j,1)-q(1,i,1))
        y = get_pbc(q(1,j,2)-q(1,i,2))
        z = get_pbc(q(1,j,3)-q(1,i,3))
        r2 = x*x + y*y + z*z
        if ( r2 .le. hard_sphere_radius_squared ) then
          goto 10
        end if
      end do

      do k = 1, ndim
        q(:,i,k) = q(1,i,k)
      end do

      i = i + 1
    end do 

!@+at
!    case (GAUSSIAN_FILL)
! 
!     i = 1
!     do while ( i .le. np )
! 20  call ru_gasdev( nu )
!       q(1,i,1) = nu(1)/p_hox
!       q(1,i,2) = nu(2)/p_hoy
!       q(1,i,3) = nu(3)/p_hoz
! 
!       do j = 1, i-1
!         x = get_pbc(q(1,j,1)-q(1,i,1))
!         y = get_pbc(q(1,j,2)-q(1,i,2))
!         z = get_pbc(q(1,j,3)-q(1,i,3))
!         r2 = x*x + y*y + z*z
!         if ( r2 .le. hard_sphere_radius_squared ) then
!           goto 20
!         end if
!       end do
! 
!       do k = 1, ndim
!         q(:,i,k) = q(1,i,k)
!       end do
! 
!       i = i + 1
!     end do
! 
!    case (DWELL_GAUSSIAN_FILL)
! 
!     i = 1
!     do while ( i .le. np )
! 30  call ru_gasdev( nu )
!       q(1,i,:) = nu(:)/2
!       if(i .le. np/2) then
!         q(1,i,3) = q(1,i,3)+a_dw
!       else
!         q(1,i,3) = q(1,i,3)-a_dw
!       end if
! 
!       do j = 1, i-1
!         x = get_pbc(q(1,j,1)-q(1,i,1))
!         y = get_pbc(q(1,j,2)-q(1,i,2))
!         z = get_pbc(q(1,j,3)-q(1,i,3))
!         r2 = x*x + y*y + z*z
!         if ( r2 .le. 2.0*hard_sphere_radius_squared ) then
!           goto 30
!         end if
!       end do
! 
!       do k = 1, ndim
!         q(:,i,k) = q(1,i,k)
!       end do
! 
!       i = i + 1
!     end do
!@-at
!@@c
    case default
      stop "Unsupported lattice fill strategy requested."
  end select

end subroutine vpi_make_lattice

function get_pbc(x) result (y)
  real(kind=b8) :: x,y

  if(use_pbc) then
    y = x - p_pbc_l*(floor(x/p_pbc_l-0.5_b8)+1.0_b8)
  else
    y = x
  endif
end function get_pbc


end module vpi_lattice
!@-node:gcross.20090624144408.2046:@thin vpi_lattice.f90
!@-leo
