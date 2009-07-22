!@+leo-ver=4-thin
!@+node:gcross.20090626112946.1710:@thin vpi_defines.f90
!@@language fortran90
module vpi_defines

  use kinds
  use constants
  implicit none

  integer,save :: N_SLICE = 1000
  integer,save :: N_SLICE_SAMPLES = 15
  integer,save :: CSLICE
  integer,save :: N_PARTICLE = 1
  integer, parameter :: N_OD_PARTICLE = 2
  integer, parameter :: N_PARTICLE2 = 0
!  integer, parameter :: OD_PNUM = N_PARTICLE - N_PARTICLE2 + 1
  integer :: OD_PNUM = 1
  real(kind=b8) :: PROB_OD_PNUM = 0.0
  integer, parameter :: N_DIM = 3
  integer, parameter :: N_DIM_ROT = 3
  integer, parameter :: OD_DIM_LOW = 1
  integer, parameter :: OD_DIM_HIGH = 3
  integer, parameter :: ODDX = 1
  integer, parameter :: ODDY = 2
  integer :: N_BINS = 53
  integer :: N_BINS_FULL = 11
  integer :: N_BINS_ROT = 21
  integer :: N_BINS_ROT_XYZ = 10
  integer :: N_BINS_GOFR = 100
  integer :: dM = 128
  real(kind=b8) :: dnu = 0.05_b8
  real(kind=b8) :: dswap = 0.50_b8
  real(kind=b8) :: dtau = 5.0e-4_b8

  real(kind=b8) :: ntol_eps = 1e-2

  real(kind=b8), parameter :: drot = 2.0e-2_b8

  real(kind=b8) :: PROB_BBRIDGE_MOVE = 1.00_b8
  real(kind=b8) :: PROB_RIGID_MOVE = 0.00_b8
  real(kind=b8) :: PROB_SWAP_MOVE = 0.00_b8

  real(kind=b8), parameter :: PROB_ROT_MOVE = 0.00_b8

  integer :: nblocks = 10000
  integer :: nsteps = 10000
  integer, parameter :: nend_steps = 0
  integer, parameter :: nb_discard = 10
  real, dimension(N_DIM) :: size_x = (/6.0,6.0,6.0/)

  real, dimension(N_DIM) :: dxdn
  real, dimension(N_DIM) :: dndx
  real :: dthetadn
  real :: dndtheta
  real, dimension(N_DIM) :: dxdn_rot_xyz
  real, dimension(N_DIM) :: dndx_rot_xyz
  real :: size_rot = 2.0
  real :: dxdn_rot
  real :: dndx_rot
  real :: size_r = 6.0
  real :: drdn
  real :: dndr
  real :: size_x_GOFR = 2.0
  real :: dndx_GOFR
  real :: dxdn_GOFR

  real :: hard_sphere_radius
  real :: hard_sphere_radius_squared

  integer, parameter :: LEFT_BRIDGE = 1
  integer, parameter :: RIGHT_BRIDGE = 2
  integer, parameter :: S_BRIDGE = 3

  logical :: eval_off_diagonal = .false.
  logical :: eval_obdm_full = .false.
  logical :: eval_obdm_rz = .false.
  logical, parameter :: xz_oda = .false.
  logical, parameter :: eval_nrdm = .false.
  logical :: eval_rdm = .false.
  logical :: write_paths = .false.

  logical :: eval_column_density = .false.
  logical :: eval_full_density = .false.
  logical :: eval_correlations = .false.

  logical :: use_gfn4 = .true.
  logical, parameter :: force_space_flip = .false.
  logical, parameter :: swap_in_123 = .false.
  logical, parameter :: swap_in_12 = .true.
  logical, parameter :: swap_in_13 = .false.

  logical, parameter :: use_eval_phase = .false.
  logical, parameter :: use_eval_corr_phase = .false.
  logical :: use_eval_cfn = .false.
  logical :: eval_qsq_sum = .false.
  logical, parameter :: impose_fixed_phase = .false.
  logical, parameter :: impose_fixed_node = .false.

  logical, parameter :: use_bbridge = .true.
  logical :: use_HS_gfn = .false.

  logical :: use_HW_gfn = .false.
  real(kind=b8), dimension(3) :: hard_wall_locations = (/ 0.5_b8, 0.5_b8, 0.5_b8 /)

  logical, parameter :: eval_rdm2_ring = .false.
  logical, parameter :: eval_rdm22_ring = .false.
  logical :: eval_obdm_angle_in_XZ_plane = .false.
  logical, parameter :: use_eval_N_R = .false.

  logical :: try_split_move = .false.
  logical,parameter :: wrap_path = .false.

  logical, parameter :: eval_rotation = .false.

  character, parameter, dimension(100) :: run_type = ""

  integer, parameter :: N_VORTEX = 4
  real, dimension( 4,2 ), parameter :: p_vcoords = reshape( (/ -0.75, 0.75, -0.75, 0.75, -.75, -.75, 0.75, 0.75/), (/ 4, 2 /) )
  real, parameter :: p_M = 0.5

  real(kind=b8), parameter :: LAMBDA = 0.5_b8

  logical :: use_lattice_file
  logical :: use_expot_file
  real(kind=b8), allocatable, dimension(:,:) :: p_lat_r0 

  real(kind=b8),parameter :: p_B = 10.0_b8

  integer, parameter :: L_MAX = 300
  real(kind=b8), dimension(10000) :: gfn_v

  logical :: use_pbc
  real(kind=b8) :: p_pbc_L

  integer, parameter :: x_axis_label = 1, y_axis_label = 2, z_axis_label = 3
  integer :: fixed_rotation_axis
  logical :: eval_2particle_angle_correlation = .false.

  namelist /configuration/ &
    N_PARTICLE, &
    N_SLICE, &
    dM, &
    dtau, &
    dnu, &
    dswap, &
    ntol_eps, &
    PROB_BBRIDGE_MOVE, &
    PROB_RIGID_MOVE, &
    PROB_SWAP_MOVE, &
    PROB_OD_PNUM, &
    OD_PNUM, &
    nblocks, &
    nsteps, &
    N_BINS, &
    N_BINS_ROT, &
    N_BINS_full, &
    N_BINS_GOFR, &
    size_x, &
    size_r, &
    size_x_gofr, &
    size_rot, &
    eval_off_diagonal, &
    write_paths, &
    eval_column_density, &
    eval_full_density, &
    eval_qsq_sum, &
    eval_rdm, &
    eval_obdm_angle_in_XZ_plane, &
    eval_obdm_full, &
    eval_obdm_rz, &
    use_eval_cfn, &
    use_HS_gfn, &
    use_HW_gfn, &
    hard_wall_locations, &
    use_gfn4, &
    use_pbc, &
    p_pbc_L, &
    use_lattice_file, &
    use_expot_file, &
    fixed_rotation_axis, &
    eval_2particle_angle_correlation

contains

subroutine init_global_parameters ()
  read(10, nml=configuration)
  write(*, nml=configuration)
end subroutine

end module vpi_defines
!@-node:gcross.20090626112946.1710:@thin vpi_defines.f90
!@-leo
