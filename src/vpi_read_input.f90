!@+leo-ver=4-thin
!@+node:gcross.20090624144408.2047:@thin vpi_read_input.f90
!@@language fortran90

module vpi_read_input
  use vpi_defines
  implicit none

contains

subroutine read_input_file(pid)
  integer :: pid

  logical :: iexist
  character (len=30) :: param_file
  character (len=50) :: pstr

  write(param_file,"(a8)") "param.in"

  inquire(file=param_file, exist=iexist)
  if (iexist) then
     ! Seed file exists, read seed from seedfile
     open(10, file=param_file, status="old")
     read(10, nml=configuration)
     close(10)
  endif

end subroutine read_input_file

subroutine read_lattice_file(pid)
  integer :: pid

  logical :: iexist
  character (len=30) :: lattice_file
  character (len=50) :: pstr

  integer :: i

  write(lattice_file,"(a8)") "lattice.in"

  allocate( p_lat_r0(N_PARTICLE, N_DIM+1) )

  inquire(file=lattice_file, exist=iexist)
  if (iexist) then
     ! Seed file exists, read seed from seedfile
     open(10, file=lattice_file, status="old")
     do i=1, N_PARTICLE
       read(10, *) pstr, p_lat_r0(i,1),  p_lat_r0(i,2),  p_lat_r0(i,3), p_lat_r0(i,4)
       write(*,*) pstr, p_lat_r0(i,1),  p_lat_r0(i,2),  p_lat_r0(i,3), p_lat_r0(i,4)
     enddo
     close(10)
  endif

end subroutine read_lattice_file

subroutine read_expot_file(pid)
  integer :: pid

  logical :: iexist
  character (len=30) :: expot_file
  character (len=50) :: pstr

  integer :: i

  stop "Reading expot file not presently supported."

!@+at
!   write(expot_file,"(a8)") "expot.in"
! 
!   inquire(file=expot_file, exist=iexist)
!   if (iexist) then
!      ! Seed file exists, read seed from seedfile
!      open(12, file=expot_file, status="old")
!      read(12, *) pstr, p_sa_a(1), p_sa_a(2), p_sa_a(3)
!      write(*,*) pstr, p_sa_a(1), p_sa_a(2), p_sa_a(3)
!      read(12, *) pstr, p_sa_b(1), p_sa_b(2), p_sa_b(3)
!      write(*,*) pstr, p_sa_b(1), p_sa_b(2), p_sa_b(3)
!      read(12, *) pstr, p_sa_c(1), p_sa_c(2), p_sa_c(3)
!      write(*,*) pstr, p_sa_c(1), p_sa_c(2), p_sa_c(3)
!      read(12, *) pstr, p_sa_d(1), p_sa_d(2), p_sa_d(3)
!      write(*,*) pstr, p_sa_d(1), p_sa_d(2), p_sa_d(3)
!      close(12)
!   endif
!@-at
!@@c
end subroutine read_expot_file

end module vpi_read_input
!@-node:gcross.20090624144408.2047:@thin vpi_read_input.f90
!@-leo
