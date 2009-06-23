module vpi_defines

  use kinds
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

  real(kind=b8), parameter :: realbignumber = 1e30

  integer, parameter :: LEFT_BRIDGE = 1
  integer, parameter :: RIGHT_BRIDGE = 2
  integer, parameter :: S_BRIDGE = 3

  logical :: eval_off_diagonal = .false.
  logical :: eval_obdm_full = .false.
  logical :: eval_obdm_rz = .true.
  logical, parameter :: xz_oda = .false.
  logical, parameter :: eval_nrdm = .false.
  logical, parameter :: eval_rdm = .false.
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

  logical, parameter :: eval_rdm2_ring = .false.
  logical, parameter :: eval_rdm22_ring = .false.
  logical, parameter :: eval_obdm_ring = .false.
  logical, parameter :: use_eval_N_R = .true.

  logical :: try_split_move = .false.
  logical,parameter :: wrap_path = .false.

  logical, parameter :: eval_rotation = .false.

  character, parameter, dimension(100) :: run_type = ""

  real(kind=b8), parameter :: M_PI = 3.141592653589793239_b8
  real(kind=b8), parameter :: M_4PI = 4.0_b8*M_PI
  real(kind=b8), parameter :: M_PI_2 = M_PI/2.0_b8
  real(kind=b8), parameter :: M_3PI_2 = 3.0_b8*M_PI/2.0_b8
  real(kind=b8), parameter :: M_2PI = 2.0_b8*M_PI
  real(kind=b8), parameter :: M_SQRT2PI = 2.506628274631_b8

  real(kind=b8) :: e_lj = 2352.5_b8
  real(kind=b8) :: a_lj = 0.05_b8
  real(kind=b8) :: a_lj2
  real(kind=b8) :: lam_ho = 1.0_b8

  real(kind=b8) :: p_hox = 0.66_b8
  real(kind=b8) :: p_hoy = 0.66_b8
  real(kind=b8) :: p_hoz = 0.66_b8
  real(kind=b8), parameter :: p_aw = 1.0_b8
  real(kind=b8) :: p_ljc5 = 1.14083e-07_b8
  real(kind=b8) :: p_ljc1 = 0.0123049_b8

  real(kind=b8), parameter :: p_cc0 = -3.32997978580414
  real(kind=b8), parameter :: p_cc1 = 0.513167460966539
  real(kind=b8), parameter :: p_cc2 = 0.067030102564702
  real(kind=b8) :: p_ac0
  real(kind=b8) :: p_ac1

  real, parameter :: e0_nua = 8
  real, parameter :: m_nua = 1.0

  real, parameter :: ep_nw = 40.0
  real, parameter :: a_nw = 0.5/3.14159265358979
  real, parameter :: lam_nw = 1.0

  real, parameter :: gw_e = 20.0_b8
  real, parameter :: gw_a = 0.15_b8
  real, parameter :: gw_z0 = 0.0_b8
  real, parameter :: gw_asq = gw_a*gw_a
  real, parameter :: gw_2asq = 2.0_b8*gw_asq
  real, parameter :: gw_e_norm = gw_e/(M_SQRT2PI*gw_a)

  real, parameter :: ap_e = 20
  real, parameter :: ap_a = 1.5
  real, parameter :: ap_asq = ap_a*ap_a
  real, parameter :: ap_2asq = 2.0*ap_asq
  real, parameter :: ap_e_norm = ap_e/(M_SQRT2PI*ap_a)
  real, parameter :: ap_lam = 1.0
  real, parameter :: ae3 = 4  ! extra harmonic confinment for second species in Annulus3

  real, parameter :: ab_e = 0
  real, parameter :: ab_a = 0.1
  real, parameter :: ab_asq = ab_a*ab_a
  real, parameter :: ab_2asq = 2.0*ab_asq
  real, parameter :: ab_e_norm = ab_e/(M_SQRT2PI*ab_a)
  real, parameter :: ab_x0 = -100

  real(kind=b8) :: a_dw = 1.5_b8
  real(kind=b8) :: a_dw2 
  real(kind=b8) :: ep_dw = 40.0_b8

  real, parameter :: e_hs = 100
  real :: a_hs = 0.01_b8
  real,save  :: a_hs2

  real(kind=b8) :: atom_qn
  real(kind=b8) :: atom_qn_p
  real(kind=b8), parameter :: atom_qe = 2.0

  integer, parameter :: N_VORTEX = 4
  real, dimension( 4,2 ), parameter :: p_vcoords = reshape( (/ -0.75, 0.75, -0.75, 0.75, -.75, -.75, 0.75, 0.75/), (/ 4, 2 /) )
  real, parameter :: p_M = 0.5

  real(kind=b8), parameter :: domega = 1.0001_b8

!Morse oscillator parameters where
!H = -0.5T + p_MO_De*( exp(-2*p_MO_a*(r-p_MO_r0)) -2*exp(-p_MO_a*(r-p_MO_r0)) )
  real(kind=b8), parameter :: p_MO_a = 10.0_b8
  !real(kind=b8), parameter :: p_MO_ap = 1.010971279201_b8 ! for evaluation of overlap integral
  real(kind=b8), parameter :: p_MO_ap = 10.0_b8 ! for evaluation of overlap integral
  real(kind=b8), parameter :: p_MO_r0 = 1.0_b8
  real(kind=b8), parameter :: p_MO_De = 50.0_b8
  real(kind=b8), parameter :: p_MO_Dep = 50.0_b8! for evaluation of overlap integral
! variational parameters for Morse potential trial function
  real(kind=b8), parameter :: p_MO_vpa = 15.29_b8
  real(kind=b8), parameter :: p_MO_vpb = 6.82_b8

  real(kind=b8), parameter :: tscw = 5.0e-2_b8
  real(kind=b8), parameter :: tsca = 1e9_b8
  real(kind=b8), parameter :: p_sc_w = 1./(2*tscw**2)
  real(kind=b8), parameter :: p_sc_a = tsca/(tscw*M_SQRT2PI)

  real(kind=b8), parameter :: LAMBDA = 0.5_b8

  real(kind=b8) :: p_dw_f0 =  0.877013259465906
  real(kind=b8) :: p_dw_f1 =  3.7487604667038
  real(kind=b8) :: p_dw_f2 =  1.4485743834981
  real(kind=b8) :: p_dw_f3 =  0.0
  real(kind=b8) :: p_dw_f4 =  0.0

  real(kind=b8),parameter :: p_nw_vb = 10.0_b8
  real(kind=b8),parameter :: p_nw_l = M_PI

  logical :: use_lattice_file
  logical :: use_expot_file
  real(kind=b8) :: p_lattice_vb = 80
  real(kind=b8),parameter :: p_lattice_ax = M_PI
  real(kind=b8),parameter :: p_lattice_ay = M_PI
  real(kind=b8),parameter :: p_lattice_az = M_PI
  real(kind=b8),parameter :: p_lattice_phase_x = M_PI/2.0_b8
  real(kind=b8),parameter :: p_lattice_phase_y = M_PI/2.0_b8
  real(kind=b8),parameter :: p_lattice_phase_z = M_PI/2.0_b8
  real(kind=b8),parameter :: p_lwx = 1.0_b8
  real(kind=b8),parameter :: p_lwz = 1.0_b8
  real(kind=b8) :: p_lat_w = 3.0_b8
  real(kind=b8), allocatable, dimension(:,:) :: p_lat_r0 

  real(kind=b8),parameter :: p_B = 10.0_b8
  real(kind=b8) :: e_dimer = 1e-4_b8

  real(kind=b8),parameter :: a1wlink =  1.0_b8
  real(kind=b8),parameter :: a2wlink =  4.0_b8
  real(kind=b8),parameter :: vbwlink =  100.0_b8
  real(kind=b8),parameter :: w0wlink =  0.01_b8

  integer, parameter :: L_MAX = 300
  real(kind=b8), dimension(10000) :: gfn_v

  real(kind=b8),parameter :: abox =  0.5_b8


  real(kind=b8), dimension(3) :: p_sa_a
  real(kind=b8), dimension(3) :: p_sa_b
  real(kind=b8), dimension(3) :: p_sa_c
  real(kind=b8), dimension(3) :: p_sa_d

  real(kind=b8) :: r_cylinder = 1.0_b8

  logical :: use_pbc
  real(kind=b8) :: p_pbc_L

  integer, parameter :: PNUM_ljc1 = 1
  integer, parameter :: PNUM_ljc5 = 2
  integer, parameter :: num_jas_params = 2

end module vpi_defines
