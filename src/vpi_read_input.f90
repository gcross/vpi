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
     read(10, *) pstr, N_PARTICLE
     write(*,*) pstr, N_PARTICLE
     read(10, *) pstr, N_SLICE
     write(*,*) pstr, N_SLICE
     read(10, *) pstr, dM
     write(*,*) pstr, dM
     read(10, *) pstr, dtau
     write(*,*) pstr, dtau
     read(10, *) pstr, dnu
     write(*,*) pstr, dnu
     read(10, *) pstr, dswap
     write(*,*) pstr, dswap
     read(10, *) pstr, ntol_eps
     write(*,*) pstr, ntol_eps
     read(10, *) pstr, PROB_BBRIDGE_MOVE
     write(*,*) pstr, PROB_BBRIDGE_MOVE
     read(10, *) pstr, PROB_RIGID_MOVE
     write(*,*) pstr, PROB_RIGID_MOVE
     read(10, *) pstr, PROB_SWAP_MOVE
     write(*,*) pstr, PROB_SWAP_MOVE
     read(10, *) pstr, PROB_OD_PNUM
     write(*,*) pstr, PROB_OD_PNUM
     read(10, *) pstr, OD_PNUM
     write(*,*) pstr, OD_PNUM
     read(10, *) pstr, nblocks
     write(*,*) pstr, nblocks
     read(10, *) pstr, nsteps
     write(*,*) pstr, nsteps

     read(10, *) pstr, N_BINS
     write(*,*) pstr, N_BINS
     read(10, *) pstr, N_BINS_ROT
     write(*,*) pstr, N_BINS_ROT
     read(10, *) pstr, N_BINS_full
     write(*,*) pstr, N_BINS_full
     read(10, *) pstr, N_BINS_GOFR
     write(*,*) pstr, N_BINS_GOFR
     
     read(10, *) pstr, size_x
     write(*,*) pstr, size_x
     read(10, *) pstr, size_r
     write(*,*) pstr, size_r
     read(10, *) pstr, size_x_gofr
     write(*,*) pstr, size_x_gofr
     read(10, *) pstr, size_rot
     write(*,*) pstr, size_rot

     read(10, *) pstr, eval_off_diagonal
     write(*,*) pstr, eval_off_diagonal
     read(10, *) pstr, write_paths
     write(*,*) pstr, write_paths
     read(10, *) pstr, eval_column_density
     write(*,*) pstr, eval_column_density
     read(10, *) pstr, eval_full_density
     write(*,*) pstr, eval_full_density
     read(10, *) pstr, eval_qsq_sum
     write(*,*) pstr, eval_qsq_sum
     read(10, *) pstr, use_eval_cfn
     write(*,*) pstr, use_eval_cfn
     read(10, *) pstr, use_HS_gfn
     write(*,*) pstr, use_HS_gfn
     read(10, *) pstr, use_gfn4
     write(*,*) pstr, use_gfn4
     read(10, *) pstr, use_pbc
     write(*,*) pstr, use_pbc
     read(10, *) pstr, p_pbc_L
     write(*,*) pstr, p_pbc_L
     read(10, *) pstr, use_lattice_file
     write(*,*) pstr, use_lattice_file
     read(10, *) pstr, use_expot_file
     write(*,*) pstr, use_expot_file

     read(10, *) pstr, a_hs
     write(*,*) pstr, a_hs

     read(10, *) pstr, p_ljc5
     write(*,*) pstr, p_ljc5
     read(10, *) pstr, p_ljc1
     write(*,*) pstr, p_ljc1

     read(10, *) pstr, a_lj
     write(*,*) pstr, a_lj
     a_lj2 = a_lj*a_lj
     read(10, *) pstr, e_lj
     write(*,*) pstr, e_lj

     read(10, *) pstr, lam_ho
     write(*,*) pstr, lam_ho
     read(10, *) pstr, a_dw
     write(*,*) pstr, a_dw
     read(10, *) pstr, ep_dw
     write(*,*) pstr, ep_dw
     read(10, *) pstr, atom_qn
     write(*,*) pstr, atom_qn
     read(10, *) pstr, p_ac0
     write(*,*) pstr, p_ac0
     read(10, *) pstr, p_ac1
     write(*,*) pstr, p_ac1
     read(10, *) pstr, p_lattice_vb
     write(*,*) pstr, p_lattice_vb
     read(10, *) pstr, p_hox
     write(*,*) pstr, p_hox
     read(10, *) pstr, p_hoy
     write(*,*) pstr, p_hoy
     read(10, *) pstr, p_hoz
     write(*,*) pstr, p_hoz
     read(10, *) pstr, p_lat_w
     write(*,*) pstr, p_lat_w
     read(10, *) pstr, p_dw_f0
     write(*,*) pstr, p_dw_f0
     read(10, *) pstr, p_dw_f1
     write(*,*) pstr, p_dw_f1
     read(10, *) pstr, p_dw_f2
     write(*,*) pstr, p_dw_f2
     read(10, *) pstr, e_dimer
     write(*,*) pstr, e_dimer
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

  write(expot_file,"(a8)") "expot.in"

  inquire(file=expot_file, exist=iexist)
  if (iexist) then
     ! Seed file exists, read seed from seedfile
     open(12, file=expot_file, status="old")
     read(12, *) pstr, p_sa_a(1), p_sa_a(2), p_sa_a(3)
     write(*,*) pstr, p_sa_a(1), p_sa_a(2), p_sa_a(3)
     read(12, *) pstr, p_sa_b(1), p_sa_b(2), p_sa_b(3)
     write(*,*) pstr, p_sa_b(1), p_sa_b(2), p_sa_b(3)
     read(12, *) pstr, p_sa_c(1), p_sa_c(2), p_sa_c(3)
     write(*,*) pstr, p_sa_c(1), p_sa_c(2), p_sa_c(3)
     read(12, *) pstr, p_sa_d(1), p_sa_d(2), p_sa_d(3)
     write(*,*) pstr, p_sa_d(1), p_sa_d(2), p_sa_d(3)
     close(12)
  endif

end subroutine read_expot_file

end module vpi_read_input
