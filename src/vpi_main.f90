!@+leo-ver=4-thin
!@+node:gcross.20090623152316.3:@thin vpi_main.f90
!@@language fortran90
!@@tabwidth -2
program Test_VPI

  !@  << Imported modules >>
  !@+node:gcross.20090623152316.5:<< Imported modules >>
  use vpi_rand_utils
  use vpi_coordinate_utils
  use vpi_gfn
  use vpi_phase_maps, update_phase => vortex
  use vpi_bisect
  use vpi_bbridge
  use vpi_sample
  use vpi_rotation
  use vpi_defines

  use vpi_lattice
  use vpi_xij
  use timers

  use vpi_single_particle_potential
  use vpi_two_body_potential
  use vpi_rotational_potential
  use vpi_single_particle_trial
  use vpi_jastrow_trial

  use kinds
  use constants
  !@-node:gcross.20090623152316.5:<< Imported modules >>
  !@nl

  implicit none

!@@raw
#ifdef USE_MPI
  include 'mpif.h'
#endif
!@@end_raw

  !@  << Global Variables >>
  !@+node:gcross.20090623152316.6:<< Global variables >>
  INTEGER :: short, medium, long, vlong
  PARAMETER (short = SELECTED_INT_KIND(2), &
             medium= SELECTED_INT_KIND(4), &
             long  = SELECTED_INT_KIND(10),&
             vlong = SELECTED_INT_KIND(100))

  real (kind=b8), allocatable, dimension( :, :, : ) :: q0, q1, qobs
  real (kind=b8), allocatable, dimension( :, :, : ) :: q_opt
  real (kind=b8), allocatable, dimension( :, :, : ) :: q_rot0, q_rot1
  real (kind=b8), allocatable, dimension( :, : ) :: qsq_sum
  real (kind=b8), allocatable, dimension( : ) :: mcount
  real (kind=b8), allocatable, dimension( : ) :: accr_v
  real (kind=b8), allocatable, dimension( :, : ) :: sp_qtest
  real (kind=b8), allocatable, dimension( : ) :: tmp_sp_chain
  real (kind=b8), allocatable, dimension( :, : ) :: ss_qtest

  real (kind=b8), allocatable, dimension( :, :, : ) :: xij2_0
  real (kind=b8), allocatable, dimension( :, :, : ) :: xij2_1

  real (kind=b8), allocatable, dimension( :, : ) :: U_0
  real (kind=b8), allocatable, dimension( :, : ) :: U_1
  real (kind=b8), allocatable, dimension( : ) :: U_cum
  real (kind=b8), allocatable, dimension( : ) :: U_avg
  real (kind=b8), allocatable, dimension( : ) :: gradU2_0
  real (kind=b8), allocatable, dimension( : ) :: gradU2_1

  real (kind=b8), allocatable, dimension( :, : ) :: grad_lntfn
  real (kind=b8), allocatable, dimension( :, : ) :: grad_lnjas
  real (kind=b8) :: lap_lntfn
  real (kind=b8) :: lap_lnjas

  real (kind=b8), allocatable, dimension( : ) :: phase_0
  real (kind=b8), allocatable, dimension( : ) :: phase_sum_0
  real (kind=b8), allocatable, dimension( : ) :: phase_1
  real (kind=b8) :: dphase

  real (kind=b8), allocatable, dimension( : , : ) :: corr_x
  real (kind=b8), allocatable, dimension( : , : ) :: corr_xz
  real (kind=b8), allocatable, dimension( : ) :: corr_r
  real (kind=b8), allocatable, dimension( : ) :: corr_phase
  real (kind=b8), allocatable, dimension( : ) :: U_weight, gU2_weight

  real (kind=b8), allocatable, dimension( : ) :: N_R
  real (kind=b8), allocatable, dimension( : ) :: N_zs
  real (kind=b8), allocatable, dimension( : ) :: N_za
  real (kind=b8), allocatable, dimension( :, : ) :: corr_N_R
  real (kind=b8), allocatable, dimension( :, : ) :: nrdm_z

  real (kind=b8), allocatable, dimension( : , : ) :: rho
  real (kind=b8), allocatable, dimension( : , : ) :: rho_red
  real (kind=b8), allocatable, dimension( : ) :: rho_r
  real (kind=b8), allocatable, dimension( : ) :: gofr
  real (kind=b8), allocatable, dimension( : ) :: gofz
  real (kind=b8), allocatable, dimension( : ) :: gofrho
  integer, allocatable, dimension( : ) :: goftheta
  integer, allocatable, dimension( : ) :: goftheta_red
  real (kind=b8), allocatable, dimension( : ) :: gofr_red
  real (kind=b8), allocatable, dimension( : ) :: gofr_rot
  real (kind=b8), allocatable, dimension( : ) :: gof_rot
  real (kind=b8), allocatable, dimension( : , : , : , : ) :: gof_rot_xyz
  real (kind=b8), allocatable, dimension( : , : , : ) :: gof_rot_xyz_avg
  real (kind=b8), allocatable, dimension( : , : ) :: trial_rho_L
  real (kind=b8), allocatable, dimension( : , : ) :: trial_rho_R
  real (kind=b8), allocatable, dimension( : , : , : ) :: full_density
  real (kind=b8), allocatable, dimension( : , : ) :: gofzrho
  real (kind=b8), allocatable, dimension( : , : ) :: odxy
  real (kind=b8), allocatable, dimension( : , : ) :: odxy_avg
  real (kind=b8), allocatable, dimension( : , : ) :: odxy_sig
  real (kind=b8), allocatable, dimension( : , : ) :: odxz
  real (kind=b8), allocatable, dimension( : , : ) :: odxz_avg
  real (kind=b8), allocatable, dimension( : , : ) :: odxz_sig
  real (kind=b8), allocatable, dimension( : , : ) :: odyz
  real (kind=b8), allocatable, dimension( : , : ) :: odyz_avg
  real (kind=b8), allocatable, dimension( : , : ) :: odyz_sig
  real (kind=b8), allocatable, dimension( : , : ) :: od2
  real (kind=b8), allocatable, dimension( : , : ) :: od2_avg
  real (kind=b8), allocatable, dimension( : , : ) :: od2_sig
  real (kind=b8), allocatable, dimension( : , : ) :: dxlxr
  real (kind=b8), allocatable, dimension( : , : ) :: rdm2_z
  real (kind=b8), allocatable, dimension( : , : ) :: rdm2_theta
  real (kind=b8), allocatable, dimension( : , : ) :: rdm22_theta
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_x
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_y
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_z
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_full
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_rz
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_tmp
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_ftmp
  !  real (kind=b8), allocatable, dimension( : , : , : , : , : , : ) :: obdm_xyz
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_rot_z
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_theta

  real (kind=b8) :: t0,t1,t2,t3

  real (kind=b8) :: accr = 0.0
  real (kind=b8) :: ccnt = 0.0
  real (kind=b8) :: sp_accr = 0.0
  real (kind=b8) :: swap_accr = 0.0
  real (kind=b8) :: bi_accr = 0.0
  real (kind=b8) :: end_accr = 0.0
  real (kind=b8) :: split_accr = 0.0
  real (kind=b8) :: E_left
  real (kind=b8) :: E_right
  real (kind=b8) :: E_center
  real (kind=b8) :: e_diff
  real (kind=b8) :: e_red
  real (kind=b8) :: E_l
  real (kind=b8) :: dE = 0
  real (kind=b8) :: dn_rms = 0
  real (kind=b8) :: dn_rms_red = 0

  integer :: i,j,ik,k,l,m,ii,jj,kk,ll,pass
  integer :: coll_iter,corpcnt,corper,opt_nconfig_total,opt_nconfig_max

  integer(long) :: n_moves, n_moves_pb, sp_moves, swap_moves, bi_moves, end_moves, split_moves
  integer(long) :: left_moves, right_moves
  integer :: mtype
  integer :: sp_try, swap_try, bi_try, end_try, split_try 
  integer :: part_num, move_start, move_end
  integer :: tmp_move_start, tmp_move_end
  integer :: sl_start, sl_end
  integer :: islice
  integer :: dmfix = 1
  integer :: maxfix_iter = 10

  real(kind=b8), dimension(2) :: tx0,tx1
  real, dimension(2) :: tsize,tdxdn,tdndx

  logical :: reject_flag = .true.
  real(kind=b8) :: weight_0, weight_1

  integer :: ierr
  integer :: my_rank = 0
  logical :: ex
  double precision time0, time1
  character (len=50) :: my_fname
  INTEGER no_fn_evals,no_fn_evals_nl2sol
  !@-node:gcross.20090623152316.6:<< Global variables >>
  !@nl

  !@  << Initialization >>
  !@+node:gcross.20090623152316.8:<< Initialization >>
  !@@raw
#ifdef USE_MPI
  !@@end_raw
  !@<< MPI Initialization >>
  !@+node:gcross.20090623152316.9:<< MPI Initialization >>
    integer :: num_procs
    integer :: source
    integer :: dest
    integer :: tag=0
    character (len = 100) :: message
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer :: instr_size, outstr_size

    call MPI_Init(ierr)

    call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  !   write(6,*)'hi, this is my rank',my_rank
    call MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)

  !@-node:gcross.20090623152316.9:<< MPI Initialization >>
  !@nl
  !@@raw
#endif
  !@@end_raw

  print *, "*** INITIALIZE ***"

  !@<< Read in configuration >>
  !@+node:gcross.20090623152316.10:<< Read in configuration >>
  !@+at
  ! Read in the configuration parameters.
  !@-at
  !@@c
    call read_configuration_file(my_rank)
  !@+at
  ! Read in optional specifications of the lattice.
  !@-at
  !@@c
  if(use_lattice_file) then
    call read_lattice_file(my_rank)
  endif
  if(use_expot_file) then
    call read_expot_file(my_rank)
  endif
  !@-node:gcross.20090623152316.10:<< Read in configuration >>
  !@nl

  !@<< Initialize variables >>
  !@+node:gcross.20090623152316.11:<< Initialize variables >>
    dMfix = 16

    corpcnt = 0
    corper = N_PARTICLE*n_slice/dM
    opt_nconfig_max = n_slice/corper
    opt_nconfig_total = 0

    n_slice_samples = n_slice

    CSLICE = N_SLICE/2
    hard_sphere_radius_squared = hard_sphere_radius*hard_sphere_radius
    dxdn = size_x / N_BINS
    dndx = 1.0/dxdn
    dxdn_rot_xyz = size_x / N_BINS_ROT_XYZ
    dndx_rot_xyz = 1.0/dxdn_rot_xyz
    dxdn_rot = size_rot / N_BINS_ROT
    dndx_rot = 1.0 / dxdn_rot
    drdn = size_r / N_BINS
    dndr = 1.0/drdn
    dndx_GOFR = (N_BINS_GOFR-1) / size_x_GOFR
    dxdn_GOFR = size_x_GOFR / (N_BINS_GOFR-1)
    dthetadn = M_PI / N_BINS
    dndtheta = 1.0 / dthetadn

    if(nsteps < corper) then
      print *, "WARNING nsteps < corper : increase nsteps to collect stats in each block"
    end if
  !@-node:gcross.20090623152316.11:<< Initialize variables >>
  !@nl

  !@<< Allocate memory for arrays >>
  !@+node:gcross.20090623152316.12:<< Allocate memory for arrays >>
    allocate( q0(N_SLICE, N_PARTICLE, N_DIM) )
    allocate( q1(N_SLICE, N_PARTICLE, N_DIM) )
    allocate( qobs(N_SLICE, N_PARTICLE, N_DIM) )
    allocate( q_opt(opt_nconfig_max, N_PARTICLE, N_DIM) )
    allocate( q_rot0(N_SLICE, N_PARTICLE, N_DIM_ROT) )
    allocate( q_rot1(N_SLICE, N_PARTICLE, N_DIM_ROT) )
    allocate( qsq_sum( N_SLICE_SAMPLES, N_DIM ) )
    allocate( mcount( N_SLICE ) )
    allocate( accr_v( N_SLICE ) )
    allocate( sp_qtest( N_SLICE, N_DIM ) )
    allocate( tmp_sp_chain( N_SLICE ) )
    allocate( ss_qtest( N_PARTICLE, N_DIM ) )
    allocate( xij2_0( N_SLICE, N_PARTICLE, N_PARTICLE ) )
    allocate( xij2_1( N_SLICE, N_PARTICLE, N_PARTICLE ) )
    allocate( U_0( N_SLICE, N_PARTICLE ) )
    allocate( U_1( N_SLICE, N_PARTICLE ) )
    allocate( U_cum( N_SLICE_SAMPLES ) )
    allocate( U_avg( N_SLICE_SAMPLES ) )
    allocate( gradU2_0( N_SLICE ) )
    allocate( gradU2_1( N_SLICE ) )
    allocate( grad_lntfn( N_PARTICLE, N_DIM ) )
    allocate( grad_lnjas( N_PARTICLE, N_DIM ) )
    allocate( phase_0( N_SLICE ) )
    allocate( phase_sum_0( N_SLICE ) )
    allocate( phase_1( N_SLICE ) )
    if(use_eval_cfn) then
      allocate( corr_x( N_SLICE, N_DIM ) )
      allocate( corr_xz( N_SLICE, 2 ) )
      allocate( corr_r( N_SLICE ) )
      allocate( corr_phase( N_SLICE ) )
      allocate( corr_N_R( N_SLICE, N_PARTICLE + 1 ) )
    end if
    allocate( N_R( N_PARTICLE + 1 ) )
    allocate( N_zs( N_PARTICLE2 + 1 ) )
    allocate( N_za( N_PARTICLE2 + 1 ) )
    allocate( nrdm_z( N_OD_PARTICLE**2, N_OD_PARTICLE**2 ) )
    allocate( U_weight( N_SLICE ) )
    allocate( gU2_weight( N_SLICE ) )

    allocate( rho( N_BINS, N_DIM ) )
    allocate( rho_red( N_BINS, N_DIM ) )
    allocate( rho_r( N_BINS_GOFR ) )
    allocate( gofr( N_BINS_GOFR ) )
    allocate( gofz( N_BINS_GOFR ) )
    allocate( gofrho( N_BINS_GOFR ) )
    allocate( gofr_red( N_BINS_GOFR ) )
    if(eval_rotation) then
      allocate( gofr_rot( N_BINS_GOFR ) )
      allocate( gof_rot( N_BINS_ROT ) )
      allocate( gof_rot_xyz( N_BINS_ROT_XYZ, N_BINS_ROT_XYZ, N_BINS_ROT_XYZ, N_BINS_ROT ) )
      allocate( gof_rot_xyz_avg(N_BINS_ROT_XYZ, N_BINS_ROT_XYZ, N_BINS_ROT_XYZ ) )
    end if
    if(eval_2particle_angle_correlation) then
      allocate( goftheta( N_BINS ) )
      if(my_rank .eq. 0) then
          allocate( goftheta_red( N_BINS ) )
      end if
    end if
    allocate( trial_rho_L( N_BINS, N_DIM ) )
    allocate( trial_rho_R( N_BINS, N_DIM ) )
    if(eval_full_density) then
      allocate( full_density( N_BINS, N_BINS, N_BINS ) )
    end if
    if(eval_column_density) then
      allocate( odxy( N_BINS, N_BINS ) )
      allocate( odxz( N_BINS, N_BINS ) )
      allocate( odyz( N_BINS, N_BINS ) )
      allocate( odxy_avg( N_BINS, N_BINS ) )
      allocate( odxz_avg( N_BINS, N_BINS ) )
      allocate( odyz_avg( N_BINS, N_BINS ) )
      allocate( odxy_sig( N_BINS, N_BINS ) )
      allocate( odxz_sig( N_BINS, N_BINS ) )
      allocate( odyz_sig( N_BINS, N_BINS ) )
    end if
    allocate( gofzrho( N_BINS, N_BINS ) )
    allocate( od2( N_BINS, N_BINS ) )
    allocate( od2_avg( N_BINS, N_BINS ) )
    allocate( od2_sig( N_BINS, N_BINS ) )
    allocate( rdm2_z(N_BINS, N_BINS ) )
    allocate( rdm2_theta( N_BINS, N_BINS ) )
    allocate( rdm22_theta( N_BINS, N_BINS ) )
    if(eval_off_diagonal) then
      allocate( dxlxr( N_BINS, N_DIM ) )
      allocate( obdm_x(N_BINS, N_BINS ) )
      allocate( obdm_y(N_BINS, N_BINS ) )
      allocate( obdm_z(N_BINS, N_BINS ) )
      if(eval_obdm_full) then
        allocate( obdm_full(N_BINS_FULL**n_dim, N_BINS_FULL**n_dim ) )
        allocate( obdm_ftmp(N_BINS_FULL**n_dim, N_BINS_FULL**n_dim ) )
      else
        if(eval_obdm_rz) then
          allocate( obdm_rz(N_BINS_FULL**2, N_BINS_FULL**2 ) )
          allocate( obdm_ftmp(N_BINS_FULL**2, N_BINS_FULL**2 ) )
        end if
      end if
      allocate( obdm_tmp(N_BINS, N_BINS ) )
  !  allocate( obdm_xyz(N_BINS, N_BINS, N_BINS, N_BINS, N_BINS, N_BINS ) )
      allocate( obdm_rot_z(N_BINS, N_BINS ) )
      if(eval_obdm_angle_in_XZ_plane) then
        allocate( obdm_theta(N_BINS, N_BINS ) )
      end if
    end if
  !@nonl
  !@-node:gcross.20090623152316.12:<< Allocate memory for arrays >>
  !@nl

  !@<< Set statistical counters to zero >>
  !@+node:gcross.20090626131222.1708:<< Set statistical counters to zero >>
  mcount = 0
  accr_v = 0
  U_cum = 0
  if(use_eval_cfn) then
    corr_x = 0.0
    corr_xz = 0.0
    corr_phase = 0.0
    corr_N_R = 0.0
  endif
  rho = 0.0
  rho_r = 0.0
  gofr = 0.0
  gofz = 0.0
  gofrho = 0.0
  if(eval_rotation) then
    gofr_rot = 0.0
    gof_rot = 0.0
  end if
  if(eval_2particle_angle_correlation) then
    goftheta = 0
  end if
  trial_rho_L = 0.0
  trial_rho_R = 0.0
  N_R = 0.0
  N_zs = 0.0
  N_za = 0.0
  U_0 = 0.0
  U_1 = 0.0

  accr = 0
  sp_accr = 0
  swap_accr = 0
  bi_accr = 0
  end_accr = 0
  split_accr = 0

  n_moves = 0
  n_moves_pb = 0
  sp_moves = 0
  swap_moves = 0
  bi_moves = 0
  end_moves = 0
  left_moves = 0
  right_moves = 0 
  split_moves = 0

  grad_lnjas = 0
  grad_lntfn = 0
  lap_lnjas = 0
  lap_lntfn = 0
  !@nonl
  !@-node:gcross.20090626131222.1708:<< Set statistical counters to zero >>
  !@nl

  !@<< Initialize "slice weights" >>
  !@+node:gcross.20090626131222.1709:<< Initialize "slice weights" >>
  ! these "slice weights" are used to give the proper weight to even
  ! and odd terms when calculating the Green's function.  
  ! in the fourth order short time approximation for the Green's funciton,
  ! odd slices contribute twice the potential energy of even slices
  ! the gradient of the potential energy only contributes for odd slides.
  U_weight = 0.0
  gU2_weight = 0.0
  ! the fourth order propagator works in sets of 3 like this
  ! 0.5 1 0.5 , 0.5 1 0.5
  ! we concatentate the adjacet half weighted steps (unless we are at the end of a path)
  ! in order for paths to be broken properly at the center we need CSLICE to be odd
  if (mod(CSLICE,2) .ne. 1) then
    print *,"ERROR: CSLICE = N_SLICE/2 must be odd"
    stop
  end if
  if (dM > cslice-1 .and. eval_off_diagonal) then
    print *,"ERROR: dM > N_slice/2"
    stop
  end if

  if (use_gfn4) then
    do ii = 1, CSLICE
      U_weight(ii) = mod(ii+1,2) + 1 
      gU2_weight(ii) = mod(ii+1,2)
    end do
    ! the center slice is doubled to handle broken paths 
    ! (probably a less kludgy way to do this but... )
    do ii = CSLICE+2,n_slice
      U_weight(ii) = mod(ii,2) + 1 
      gU2_weight(ii) = mod(ii,2)
    end do
  else
    U_weight = 1.0
    gU2_weight = 0.0
  end if
  U_weight(1) = 0.5_b8
  U_weight(N_SLICE) = 0.5_b8
  U_weight(CSLICE) = 0.5_b8
  U_weight(CSLICE+1) = 0.5_b8
  !@nonl
  !@-node:gcross.20090626131222.1709:<< Initialize "slice weights" >>
  !@nl

    !@  << Seed random number generator >>
    !@+node:gcross.20090626131222.1710:<< Seed random number generator >>
    ! read in seed file / create a new seed for random number generator
    call ru_init_seed_file(my_rank)
    if (eval_rotation) then
      gfn_v = init_rot_gfn (size(gfn_v), L_MAX, dtau)
      call write_rot_gfn(gfn_v)
    end if
    !@-node:gcross.20090626131222.1710:<< Seed random number generator >>
    !@nl

    !@  << Set the physical system to its initial state >>
    !@+node:gcross.20090626131222.1715:<< Set the physical system to its initial state >>
    !@<< Set up the positional degrees of freedom >>
    !@+node:gcross.20090626131222.1711:<< Set up the positional degrees of freedom >>
    q0 = 0.0
    q1 = 0.0

    !print *, q0(1,:,:)
    !print *, "pfunc = ", pfunc(q0(1,:,:),N_PARTICLE, N_DIM)
    !print *, "test: ", test_fvar(q0,pfunc)

    write(my_fname,"(a14,i4.4)")"last_path.dat.",my_rank
    inquire(file=my_fname, exist=ex)
    if ( ex ) then
      !@  << Read initial path from file with same rank as myself >>
      !@+node:gcross.20090626131222.1712:<< Read initial path from file with same rank as myself >>
      open(10, file=my_fname, action="read")
      do i = 1,N_SLICE
        read(UNIT=10 ,END=111 ,FMT=* ) q0(i,:,:)
        !print *, q0(i,:,:)
      end do
      111   close(10)
      if(i .lt. N_SLICE) then
        do j = 1, N_SLICE
          k = j*float(i-1)/float(N_SLICE)
          if ( k .le. i-1 ) then
            if( k .ge. 1 ) then
              q1(j,:,:) = q0(k,:,:) 
            else
              q1(j,:,:) = q0(1,:,:)
            end if
          else
            q1(j,:,:) = q0(i-1,:,:)
          end if
        end do
        if((i-1)/2+1 < i-1) then
          q1(CSLICE,:,:) = q0((i-1)/2,:,:)
          q1(CSLICE+1,:,:) = q0((i-1)/2+1,:,:)
        end if
        do ii =1, n_particle
          if( .not.((ii .eq. od_pnum) .and. eval_off_diagonal) ) then  
            q1(cslice+1,ii,:) = q1(cslice,ii,:)
          endif
        end do
        q0 = q1
      endif

      !@-node:gcross.20090626131222.1712:<< Read initial path from file with same rank as myself >>
      !@nl
    else
      write(my_fname,"(a14,i4.4)")"last_path.dat.",0
      inquire(file=my_fname, exist=ex)
      if ( ex ) then
        !@    << Read initial path from file tagged with rank 0 >>
        !@+node:gcross.20090626131222.1714:<< Read initial path from file tagged with rank 0 >>
        open(10, file=my_fname, action="read")
        do i = 1,N_SLICE
          read(UNIT=10 ,END=113 ,FMT=* ) q0(i,:,:)
          !print *, q0(i,:,:)
        end do
        113     close(10)
        if(i .lt. N_SLICE) then
          do j = 1, N_SLICE
            k = j*float(i-1)/float(N_SLICE)
            if ( k .le. i-1 ) then
              if( k .ge. 1 ) then
                q1(j,:,:) = q0(k,:,:) 
              else
                q1(j,:,:) = q0(1,:,:)
              end if
            else
              q1(j,:,:) = q0(i-1,:,:)
            end if
          end do
          do ii =0, n_particle
            if( .not.((ii .eq. od_pnum) .and. eval_off_diagonal) ) then  
              q1(cslice+1,ii,:) = q1(cslice,ii,:)
            endif
          end do
          q0 = q1
        endif 
        !@nonl
        !@-node:gcross.20090626131222.1714:<< Read initial path from file tagged with rank 0 >>
        !@nl
      else
        !@    << Generate initial path >>
        !@+node:gcross.20090626131222.1713:<< Generate initial path >>
        call vpi_make_lattice( q0, real(p_pbc_l), LATTICE_FILL)
        do ii =1, n_particle
          if( .not.((ii .eq. od_pnum) .and. eval_off_diagonal) ) then  
            q0(cslice+1,ii,:) = q0(cslice,ii,:)
          endif
        end do
        !@-node:gcross.20090626131222.1713:<< Generate initial path >>
        !@nl
      end if
    end if
    !@-node:gcross.20090626131222.1711:<< Set up the positional degrees of freedom >>
    !@nl

    !@<< Set up the rotational degrees of freedom >>
    !@+node:gcross.20090626131222.1716:<< Set up the rotational degrees of freedom >>
    if (eval_rotation) then
      write(my_fname,"(a18,i4.4)")"last_rot_path.dat.",my_rank
      inquire(file=my_fname, exist=ex)
      if ( ex ) then
        open(10, file=my_fname, action="read")
        do i = 1,N_SLICE
          read(UNIT=10, END=112, FMT=*) q_rot0(i,:,:)
          !print *, q0(i,:,:)
        end do
    112     close(10)

      else
        q_rot0(:,:,1) = 0.0_b8
        q_rot0(:,:,2) = 0.0_b8
        q_rot0(:,:,3) = 1.0_b8
      end if
      q_rot1 = q_rot0
    end if
    !@-node:gcross.20090626131222.1716:<< Set up the rotational degrees of freedom >>
    !@nl

    if(force_space_flip) then
      q0(CSLICE+2:N_SLICE,:,3) = -q0(CSLICE+2:N_SLICE,:,3)
    end if

    q1 = q0
    if(use_pbc) then
      qobs(:,:,:) = q0(:,:,:) - p_pbc_L*floor(q0(:,:,:)/p_pbc_L - 0.5_b8) - p_pbc_L
      q0(:,:,:) = qobs(:,:,:)
    else
      qobs(:,:,:) = q0(:,:,:)
    end if

    !@<< Update displacement matrix >>
    !@+node:gcross.20090626131222.1717:<< Update displacement matrix >>
    do ii = 1, N_PARTICLE
      if(use_pbc) then
        call vpi_update_xij_pbc( xij2_0, q0, 1, N_SLICE, ii, N_SLICE, N_PARTICLE, N_DIM  )
      else
        call vpi_update_xij( xij2_0, q0, 1, N_SLICE, ii, N_SLICE, N_PARTICLE, N_DIM  )
      end if
    end do
    !@nonl
    !@-node:gcross.20090626131222.1717:<< Update displacement matrix >>
    !@nl

    !@<< Update the phase >>
    !@+node:gcross.20090626131222.1718:<< Update the phase >>
    if( use_eval_phase ) then
      call update_phase( phase_0, 1, N_SLICE, qobs, xij2_0 )
    end if
    phase_1 = phase_0
    !@nonl
    !@-node:gcross.20090626131222.1718:<< Update the phase >>
    !@nl
    !@-node:gcross.20090626131222.1715:<< Set the physical system to its initial state >>
    !@nl

    !@  << Update the phase >>
    !@+middle:gcross.20090626131222.1715:<< Set the physical system to its initial state >>
    !@+node:gcross.20090626131222.1718:<< Update the phase >>
    if( use_eval_phase ) then
      call update_phase( phase_0, 1, N_SLICE, qobs, xij2_0 )
    end if
    phase_1 = phase_0
    !@nonl
    !@-node:gcross.20090626131222.1718:<< Update the phase >>
    !@-middle:gcross.20090626131222.1715:<< Set the physical system to its initial state >>
    !@nl

    !@  << Compute initial values of some physical quantities >>
    !@+node:gcross.20090626131222.1719:<< Compute initial values of some physical quantities >>
    call compute_potential (&
        qobs, xij2_0, q_rot0, &
        1, N_SLICE, &
        N_SLICE, N_PARTICLE, N_DIM, &
        U_0, gradU2_0, &
        reject_flag &
        )

    if (reject_flag) then
      print *,"Unfortunately, the generated initial state is invalid because there are"
      print *,"particles which overlap in a region forbidden by the 2-body potential."
      stop
    end if

    U_1(:,:) = U_0(:,:)
    gradU2_1(:) = gradU2_0(:)
    xij2_1(:,:,:) = xij2_0(:,:,:)
    !@-node:gcross.20090626131222.1719:<< Compute initial values of some physical quantities >>
    !@nl

    qsq_sum = 0.0
  !@-node:gcross.20090623152316.8:<< Initialization >>
  !@nl

    !@    << Main iteration >>
    !@+node:gcross.20090623152316.33:<< Main iteration >>
    do i = 1,nblocks
      !@  << Initialize accumulating variables to zero >>
      !@+node:gcross.20090626131222.1720:<< Initialize accumulating variables to zero >>
      n_moves_pb = 0
      dE = 0.0_b8
      dn_rms = 0.0_b8
      e_diff = 0.0_b8
      E_left = 0.0_b8
      E_center = 0.0_b8
      E_right = 0.0_b8
      !@-node:gcross.20090626131222.1720:<< Initialize accumulating variables to zero >>
      !@nl

      !@  << Start the timer >>
      !@+node:gcross.20090626091217.1689:<< Start the timer >>
      !@@raw
#ifdef USE_MPI
      !@@end_raw
      time0 = MPI_Wtime()
      !@@raw
#else
      !@@end_raw
      time0 = day_timer()
      !@@raw
#endif
      !@@end_raw
      !@-node:gcross.20090626091217.1689:<< Start the timer >>
      !@nl

    !      write(my_fname,"(a10,i4.4)")"xl-xr.dat.",my_rank
    !      open(13, file=my_fname,  POSITION="APPEND")
      do j = 1,nsteps
        !@    << Perform Monte Carlo iteration >>
        !@+node:gcross.20090626112946.1694:<< Perform Monte Carlo iteration >>
        !@<< If debugging, then dump out all particle positions >>
        !@+node:gcross.20090626112946.1708:<< If debugging, then dump out all particle positions >>
        !@@raw
#ifdef DEBUG_L0
        !@@end_raw
        do k = 1, N_SLICE
          do l = 1, N_PARTICLE
            do m = 1, N_DIM
              if(q1(k,l,m) .ne. q0(k,l,m)) then
                print *, "q0 != q1",k,l,m
              end if
            end do
          end do
        end do
        !@@raw
#endif
        !@@end_raw
        !@-node:gcross.20090626112946.1708:<< If debugging, then dump out all particle positions >>
        !@nl

        !@<< Randomly choose a move to make and apply it to the system >>
        !@+node:gcross.20090626112946.1688:<< Randomly choose a move to make and apply it to the system >>
        call sample_scheme1(q0, q1, move_start, move_end, part_num, mtype)
        if( eval_rotation ) then
          call sample_rotation(q_rot0, q_rot1, move_start, move_end, part_num)
        end if

        !@<< Impose periodic boundary conditions >>
        !@+node:gcross.20090626112946.1691:<< Impose periodic boundary conditions >>
        if(use_pbc) then
          q1(move_start:move_end,:,:) = q1(move_start:move_end,:,:) &
                                       - p_pbc_L*floor(q1(move_start:move_end,:,:)/p_pbc_L - 0.5_b8) &
                                       - p_pbc_L
        end if
        !@-node:gcross.20090626112946.1691:<< Impose periodic boundary conditions >>
        !@nl

        !        if( .not.((part_num .eq. od_pnum) .and. eval_off_diagonal) ) then  
        !          q1(cslice+1,part_num,:) = q1(cslice,part_num,:)
        !        endif

        !@<< Update move type statistics >>
        !@+node:gcross.20090626112946.1690:<< Update move type statistics >>
        mcount(move_start:move_end) =  mcount(move_start:move_end) + 1
        bi_try = 0 
        sp_try = 0 
        swap_try = 0 
        select case (mtype)
          case (MT_BBRIDGE)
            bi_try = 1 
          case (MT_RIGID)
            sp_try = 1 
          case (MT_SWAP)
            swap_try = 1 
        end select
        split_moves = split_moves + split_try
        bi_moves = bi_moves + bi_try
        end_moves = end_moves + end_try
        sp_moves = sp_moves + sp_try
        swap_moves = swap_moves + swap_try

        !        print *,"move_start, move_end", move_start, move_end
        !        do ii=1,N_SLICE
        !          print *, "q0: ",ii,  q0(ii,:,:)
        !          print *, "q1: ",ii,  q1(ii,:,:)
        !        end do
        !@nonl
        !@-node:gcross.20090626112946.1690:<< Update move type statistics >>
        !@nl

        !@<< Update the displacement matrix >>
        !@+node:gcross.20090626112946.1689:<< Update the displacement matrix >>
        do ii = 1, N_PARTICLE
          if(use_pbc) then
            call vpi_update_xij_pbc( xij2_1, q1, move_start, move_end, ii, N_SLICE, N_PARTICLE, N_DIM  )
          else
            call vpi_update_xij( xij2_1, q1, move_start, move_end, ii, N_SLICE, N_PARTICLE, N_DIM  )
          end if
        end do
        !@-node:gcross.20090626112946.1689:<< Update the displacement matrix >>
        !@nl
        !@-node:gcross.20090626112946.1688:<< Randomly choose a move to make and apply it to the system >>
        !@nl

        !@<< (Inactive experimental stuff that I don't understand) >>
        !@+node:gcross.20090626112946.1707:<< (Inactive experimental stuff that I don't understand) >>
        !@@raw
#ifdef EXPERIMENTAL
        !@@end_raw
        if(mtype == MT_SWAP .or. mtype == MT_RIGID ) then
          acc_flag = .true.
          do ii = dMfix, N_SLICE,dMfix*2
            acc_flag = .true.
            do jj = ii-dMfix,ii+dMfix
              U_ij = Uij_func( q1, xij2_1, jj, part_num, N_SLICE, N_PARTICLE, N_DIM, acc_flag )
            end do
            coll_iter = 0
        crepair:    do while ((coll_iter < maxfix_iter) .and. (acc_flag .eqv. .false.))
              tmp_move_start = ii-dMfix-1
              tmp_move_end = ii+dMfix+1
              if(tmp_move_start .le. 1) then
                tmp_move_start = 1
                call bbridge(q1, q1, part_num, 1, N_DIM, LEFT_BRIDGE, tmp_move_start,tmp_move_end, lambda, dtau)
              else 
                if(tmp_move_end .ge. N_SLICE) then
                  tmp_move_end = N_SLICE
                  call bbridge(q1, q1, part_num, 1, N_DIM, RIGHT_BRIDGE, tmp_move_start,tmp_move_end, lambda, dtau)
                else
                  call bbridge(q1, q1, part_num, 1, N_DIM, S_BRIDGE, tmp_move_start,tmp_move_end, lambda, dtau)
                end if
              end if
              call vpi_update_xij( xij2_1, q1, tmp_move_start, tmp_move_end, part_num, N_SLICE, N_PARTICLE, N_DIM  )
              acc_flag = .true.
              do jj = tmp_move_start,tmp_move_end
                U_ij = Uij_func( q1, xij2_1, jj, part_num, N_SLICE, N_PARTICLE, N_DIM, acc_flag )
              end do
              coll_iter = coll_iter + 1
        !@@raw
#ifdef DEBUG_EXP
        !@@end_raw
              if(acc_flag .eqv. .true.) then 
                print *, my_rank,'fixed step: ', j, 'slice: ', ii, 'iter: ', coll_iter
              end if
        !@@raw
#endif DEBUG_EXP
        !@@end_raw
            end do crepair
            if(acc_flag .eqv. .false.) then 
        !@@raw
#ifdef DEBUG_EXP
        !@@end_raw
              print *, my_rank,'failed step: ', j, 'slice: ', ii, 'iter: ', coll_iter
        !@@raw
#endif DEBUG_EXP
        !@@end_raw
              exit
            end if
          end do 
        end if
        !@@raw
#endif
        !@@end_raw
        !@nonl
        !@-node:gcross.20090626112946.1707:<< (Inactive experimental stuff that I don't understand) >>
        !@nl

        if ( use_eval_phase ) then
          call update_phase( phase_1, move_start, move_end, q1, xij2_1 )
        end if 

        !@<< Determine whether the move should be accepted >>
        !@+node:gcross.20090626112946.1693:<< Determine whether the move should be accepted >>
        !@<< Compute probability of acceptance >>
        !@+node:gcross.20090721121051.1755:<< Compute probability of acceptance >>
        weight_0 = compute_log_acceptance_weight (&
            q0, xij2_0, q_rot0, &
            move_start, move_end, &
            N_SLICE, N_PARTICLE, N_DIM, &
            U_0, gradU2_0, &
            reject_flag &
            )

        if(reject_flag) then
          stop "ERROR, the path that was accepted before has just been rejected!"
        end if

        weight_1 = compute_log_acceptance_weight (&
            q1, xij2_1, q_rot1, &
            move_start, move_end, &
            N_SLICE, N_PARTICLE, N_DIM, &
            U_1, gradU2_1, &
            reject_flag &
            )

        if(reject_flag) then
          goto 10
        end if
        !@-node:gcross.20090721121051.1755:<< Compute probability of acceptance >>
        !@nl

        !@<< Optionally impose the fixed phase condition >>
        !@+node:gcross.20090626112946.1698:<< Optionally impose the fixed phase condition >>
        if (impose_fixed_phase) then
          dphase = product(phase_1(move_start:move_end))/product(phase_0(move_start:move_end))
        else
          dphase = 1.0_b8
        end if
        !@-node:gcross.20090626112946.1698:<< Optionally impose the fixed phase condition >>
        !@nl

        !@<< Accept or reject the move >>
        !@+node:gcross.20090626112946.1706:<< Accept or reject the move >>
        10 if ( (.not. reject_flag) .and. vpi_accept_path( weight_0, weight_1, dphase)) then
          !@  << Accept the move >>
          !@+node:gcross.20090626112946.1699:<< Accept the move >>
          !          print *, "!!! ACCEPTED !!!"
          !@<< Update the degrees of freedom and derived quantities >>
          !@+node:gcross.20090626112946.1701:<< Update the degrees of freedom and derived quantities >>
          !@+at
          ! Update degrees of freedom
          !@-at
          !@@c

          q0(move_start:move_end,:,:) = q1(move_start:move_end,:,:)
          if(use_pbc) then
            qobs(move_start:move_end,:,:) = q0(move_start:move_end,:,:) &
                                            - p_pbc_L*floor(q0(move_start:move_end,:,:)/p_pbc_L - 0.5_b8) &
                                            - p_pbc_L
            q0(move_start:move_end,:,:) = qobs(move_start:move_end,:,:)
          else
            qobs(move_start:move_end,:,:) = q0(move_start:move_end,:,:)
          end if
          if( eval_rotation ) then
            q_rot0(move_start:move_end,:,:) = q_rot1(move_start:move_end,:,:)
          end if

          !@+at
          ! Update derived quantities
          !@-at
          !@@c

          xij2_0(move_start:move_end,:,:) = xij2_1(move_start:move_end,:,:) 
          U_0(move_start:move_end,:) = U_1(move_start:move_end,:)
          gradU2_0(move_start:move_end) =  gradU2_1(move_start:move_end)
          if( use_eval_phase ) then
            phase_0(move_start:move_end) = phase_1(move_start:move_end)
          end if
          !@-node:gcross.20090626112946.1701:<< Update the degrees of freedom and derived quantities >>
          !@nl

          !@<< Update statistics >>
          !@+node:gcross.20090626112946.1703:<< Update statistics >>
          accr = accr + 1;
          accr_v(move_start:move_end) = accr_v(move_start:move_end) + 1
          bi_accr = bi_accr + bi_try
          sp_accr = sp_accr + sp_try
          swap_accr = swap_accr + swap_try
          end_accr = end_accr + end_try
          split_accr = split_accr + split_try
          !@-node:gcross.20090626112946.1703:<< Update statistics >>
          !@nl
          !@nonl
          !@-node:gcross.20090626112946.1699:<< Accept the move >>
          !@nl
        else
          !@  << Reject the move >>
          !@+node:gcross.20090626112946.1704:<< Reject the move >>
          !          print *, "!!! REJECTED !!!"
          !@+at
          ! Reset all degrees of freedom and derived quantities.
          !@-at
          !@@c
          q1(move_start:move_end,:,:) = q0(move_start:move_end,:,:)
          if( eval_rotation ) then
            q_rot1(move_start:move_end,:,:) = q_rot0(move_start:move_end,:,:)
          end if
          xij2_1(move_start:move_end,:,:) = xij2_0(move_start:move_end,:,:)
          U_1(move_start:move_end,:) = U_0(move_start:move_end,:)
          gradU2_1(move_start:move_end) =  gradU2_0(move_start:move_end)
          !@-node:gcross.20090626112946.1704:<< Reject the move >>
          !@nl
        end if
        !@-node:gcross.20090626112946.1706:<< Accept or reject the move >>
        !@nl

        !@-node:gcross.20090626112946.1693:<< Determine whether the move should be accepted >>
        !@nl

        !@<< Compute quantities of interest >>
        !@+node:gcross.20090626112946.1692:<< Compute quantities of interest >>
        corpcnt = corpcnt + 1
        if(corpcnt >= corper) then
          call vpi_eval_gfn_sep_in_full_space( gofr, xij2_0(CSLICE,:,:), N_BINS_GOFR, dndx_GOFR )
          call vpi_eval_gfn_sep_along_Z_axis( gofz, qobs, cslice, N_BINS_GOFR, dndx_GOFR )
          call vpi_eval_gfn_sep_in_XY_plane( gofrho, qobs, cslice, N_BINS_GOFR, dndx_GOFR )
          call vpi_eval_gfn_sep_in_Z_and_XY( gofzrho, qobs, cslice, N_BINS, dndx )
          if(eval_2particle_angle_correlation) then
            call vpi_eval_gfn_sep_angle( goftheta, qobs, cslice, N_BINS, dndtheta )
          end if
          corpcnt = 0
          n_moves = n_moves + 1
          n_moves_pb = n_moves_pb + 1
          pass = grad_lap_sp_tfun( qobs, 1, N_PARTICLE, N_DIM, N_SLICE, grad_lntfn, lap_lntfn )
          pass = grad_lap_jas_tfun( qobs, xij2_0, 1, N_PARTICLE, N_DIM, N_SLICE, grad_lnjas, lap_lnjas )
          call eval_E_local(N_PARTICLE, N_DIM, grad_lntfn, lap_lntfn, grad_lnjas, lap_lnjas, sum(U_0(1,:)), E_l)
          E_left = E_left + E_l

          pass = grad_lap_sp_tfun( qobs, N_SLICE, N_PARTICLE, N_DIM, N_SLICE, grad_lntfn, lap_lntfn )
          pass = grad_lap_jas_tfun( qobs, xij2_0, N_SLICE, N_PARTICLE, N_DIM, N_SLICE, grad_lnjas, lap_lnjas )
          call eval_E_local(N_PARTICLE, N_DIM, grad_lntfn, lap_lntfn, grad_lnjas, lap_lnjas, sum(U_0(N_SLICE,:)), E_l)
          E_right = E_right + E_l

          pass = grad_lap_sp_tfun( qobs, CSLICE, N_PARTICLE, N_DIM, N_SLICE, grad_lntfn, lap_lntfn )
          pass = grad_lap_jas_tfun( qobs, xij2_0, CSLICE, N_PARTICLE, N_DIM, N_SLICE, grad_lnjas, lap_lnjas )
          call eval_E_local(N_PARTICLE, N_DIM, grad_lntfn, lap_lntfn, grad_lnjas, lap_lnjas, sum(U_0(CSLICE,:)), E_l)
          E_center = E_center + E_l
          call MPI_ALLreduce(E_l,t1,1,mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,ierr)
          e_diff = e_diff + abs(E_l-t1/num_procs)
          if(n_moves_pb < opt_nconfig_max ) then
            q_opt(n_moves_pb,:,:) = qobs(cslice,:,:)
          end if

        !         print "(56g30.15)", my_rank, e_diff,E_l,e_red/num_procs,abs(E_l-e_red/num_procs)

          t1 =  sum(U_0(CSLICE+1,:)) - sum(U_0(CSLICE,:))
          if(eval_off_diagonal .eqv. .false.) then
            if(t1 .gt. 1e-12) then
              print *,n_moves,part_num,"error: de=",t1,q0(CSLICE,part_num,:)-q0(CSLICE+1,part_num,:)
              print *,n_moves,part_num,"U_0(CSLICE) : ",sum(U_0(CSLICE,:))
              print *,n_moves,part_num,"U_0(CSLICE+1) : ",sum(U_0(CSLICE+1,:))
            end if
          end if
          dE =  dE + t1
          if(eval_qsq_sum) then
            do k = 1, N_SLICE_SAMPLES
        !            ik =  int((k-1)*float(N_SLICE)/(N_SLICE_SAMPLES-1)) + 1
              ik = k
              U_cum(k) = U_cum(k) + sum(U_0(ik,:))/(N_PARTICLE)
              U_avg(k) = U_cum(k)/n_moves
              phase_sum_0(k) = phase_sum_0(k) + phase_0(ik)
              do m = 1, N_DIM
                qsq_sum(k,m) = qsq_sum(k,m) + sum(qobs(ik,1:N_PARTICLE,m)**2)
              end do
            end do
          end if

          call vpi_eval_density( rho, qobs(CSLICE,:,:), N_BINS, size_x/2, dndx )
          call vpi_eval_radial_density( rho_r, qobs(CSLICE,:,:), n_bins_gofr, dndr ) 
        !        if ( eval_off_diagonal ) then
        !          call vpi_eval_density_single_particle( dxlxr, (qobs(CSLICE,od_pnum,:)-qobs(CSLICE+1,od_pnum,:)), N_BINS, size_x/2, dndx )
        !        end if
          call vpi_eval_density( trial_rho_L, qobs(1,:,:), N_BINS, size_x/2, dndx )
          call vpi_eval_density( trial_rho_R, qobs(N_SLICE,:,:), N_BINS, size_x/2, dndx )
          if(eval_full_density) then
            call vpi_eval_full_density(full_density, qobs(CSLICE,:,:), N_BINS, size_x/2, dndx )
          end if
          if(eval_column_density) then
            odxy_avg = odxy/( n_moves*(N_PARTICLE-N_PARTICLE2)*dxdn(1)*dxdn(2) )
            odxz_avg = odxz/( n_moves*(N_PARTICLE-N_PARTICLE2)*dxdn(1)*dxdn(3) )
            odyz_avg = odyz/( n_moves*(N_PARTICLE-N_PARTICLE2)*dxdn(2)*dxdn(3) )
            if(N_PARTICLE2 > 0 ) then
              od2_avg = od2/( n_moves*N_PARTICLE2*dxdn(1)*dxdn(3) )
            endif
            call vpi_eval_density_in_plane( odxy, odxy_avg, odxy_sig, qobs(CSLICE,:,:), 1, 2, N_BINS, size_x/2, dndx )
            call vpi_eval_density_in_plane( odxz, odxz_avg, odxz_sig, qobs(CSLICE,:,:), 1, 3, N_BINS, size_x/2, dndx )
            call vpi_eval_density_in_plane( odyz, odyz_avg, odyz_sig, qobs(CSLICE,:,:), 2, 3, N_BINS, size_x/2, dndx )
            if( N_PARTICLE2 .gt. 0 ) then 
              call vpi_eval_density_in_plane2( od2, od2_avg, od2_sig, qobs(CSLICE,:,:), N_BINS, size_x/2, dndx )
            end if
            if (i .lt. nb_discard) then
              odxy_sig = 0 
              odxz_sig = 0 
              odyz_sig = 0 
            end if
          end if


          if( eval_rotation ) then
            call vpi_eval_gof_rot( gof_rot, gof_rot_xyz, gof_rot_xyz_avg, gofr_rot, qobs, q_rot0, xij2_0, &
                                   CSLICE, N_BINS_GOFR, dndx_GOFR,  N_BINS_ROT_XYZ, size_x/2, dndx_rot_xyz, N_BINS_ROT, dndx_rot)
          end if
          if(use_eval_cfn) then
            call vpi_eval_corr_x( corr_x, 1, N_SLICE, q0 )
            call vpi_eval_corr_xz( corr_xz, 1, N_SLICE, q0 )
            call vpi_eval_corr_r( corr_r, 1, N_SLICE, q0 )
          end if
          if (use_eval_corr_phase) then
            call vpi_eval_corr_phase( corr_phase, phase_0, 1, N_SLICE )
          end if
          if( use_eval_N_R ) then
            call vpi_eval_N_R( N_R, dn_rms, qobs(:,:,:), CSLICE, 0 )
          end if
        !        call eval_N_z2( N_zs,N_za, q0(CSLICE,:,:) )
        !       do ii = 1, N_SLICE
        !         call eval_N_R( corr_N_R(ii,:), q0(ii,:,:) )
        !       end do
          if(eval_rdm) then
            do ii = 1, N_PARTICLE
              do jj = ii+1, N_PARTICLE
                call vpi_eval_obdm_single_axis( rdm2_z, qobs(CSLICE,ii,3), qobs(CSLICE,jj,3), N_BINS, N_DIM, size_x(3)/2, dndx(3))
                if(eval_rdm2_ring) then
                  call vpi_eval_obdm_angle_in_XZ_plane( rdm2_theta, qobs(CSLICE,ii,:), qobs(CSLICE,jj,:), N_BINS, N_DIM)
                end if
              end do
            end do
          end if
          if(eval_rdm22_ring) then
            do ii = N_PARTICLE-N_PARTICLE2+1, N_PARTICLE
              do jj = ii+1, N_PARTICLE
                call vpi_eval_obdm_angle_in_XZ_plane( rdm22_theta, qobs(CSLICE,ii,:), qobs(CSLICE,jj,:), N_BINS, N_DIM)
              end do
            end do
          end if
          if( eval_off_diagonal ) then
            call vpi_eval_obdm_single_axis( obdm_x, qobs(CSLICE,od_pnum,1), qobs(CSLICE+1,od_pnum,1), &
                                            N_BINS, N_DIM, size_x(1)/2, dndx(1))
            call vpi_eval_obdm_single_axis( obdm_y, qobs(CSLICE,od_pnum,2), qobs(CSLICE+1,od_pnum,2), &
                                            N_BINS, N_DIM, size_x(2)/2, dndx(2))
            call vpi_eval_obdm_single_axis( obdm_z, qobs(CSLICE,od_pnum,3), qobs(CSLICE+1,od_pnum,3), &
                                            N_BINS, N_DIM, size_x(3)/2, dndx(3))
            if(eval_obdm_full) then
              call vpi_eval_obdm_full( obdm_full, qobs(CSLICE,od_pnum,:), qobs(CSLICE+1,od_pnum,:), &
                                       N_BINS_full, N_DIM, size_x(:)/2, dndx(:))
            else
              if(eval_obdm_rz) then
                tx0(1) = qobs(CSLICE,od_pnum,3)
                tx0(2) = -sqrt( (qobs(CSLICE+1,od_pnum,1) &
                         -qobs(CSLICE,od_pnum,1))**2 &
                         + (qobs(CSLICE+1,od_pnum,2) &
                         -qobs(CSLICE,od_pnum,2))**2)/2
                tx1(1) = qobs(CSLICE+1,od_pnum,3)
                tx1(2) =  -tx0(2)
                tsize(1) = size_x(3) 
                tsize(2) = sqrt(size_x(1)**2+size_x(2)**2)
                tdxdn = tsize / N_BINS_full
                tdndx = 1.0/tdxdn
                call vpi_eval_obdm_full( obdm_rz, tx0, tx1, N_BINS_full, 2, tsize(:)/2, tdndx(:))
              end if
            end if
            if(eval_obdm_angle_in_XZ_plane) then
              call vpi_eval_obdm_angle_in_XZ_plane( obdm_theta, qobs(CSLICE,od_pnum,:), qobs(CSLICE+1,od_pnum,:), N_BINS, N_DIM)
            endif
            if( eval_rotation ) then
              call vpi_eval_obdm_single_axis( obdm_rot_z, q_rot0(CSLICE,od_pnum,2), q_rot0(CSLICE+1,od_pnum,2), &
                                      N_BINS, N_DIM, 1.0, real(N_BINS)/2.0)
            end if
          end if
          if( eval_nrdm ) then
            call vpi_eval_nrdm( nrdm_z, qobs(CSLICE,:,3), qobs(CSLICE+1,:,3) )
          end if
        end if
        !@nonl
        !@-node:gcross.20090626112946.1692:<< Compute quantities of interest >>
        !@nl

        !@-node:gcross.20090626112946.1694:<< Perform Monte Carlo iteration >>
        !@nl
      end do

      !@  << Compute and output current simulation results >>
      !@+node:gcross.20090626091217.1690:<< Compute and output current simulation results >>
      !      if( eval_off_diagonal ) then
      !        write(13, "(3g18.9)", ADVANCE="NO") qobs(CSLICE,1,:)
      !        write(13, "(3g18.9)", ADVANCE="NO") qobs(CSLICE+1,1,:)
      !        write(13, "(3g)", ADVANCE="NO") q0(CSLICE,2,:)
      !        write(13, "(3g)") q0(CSLICE+1,2,:)
      !      end if
      !      close(13)


      call MPI_REDUCE(rho,rho_red,n_bins*n_dim,mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(my_rank .eq. 0) then
        write(my_fname,"(a11)")"rho_cum.dat"
        open(12, file=my_fname, status="replace")
        do ii = 1,N_BINS
      !          write(12, "(6g30.21)") ii, rho_red(ii,:)/( num_procs*n_moves*dxdn(:)*N_PARTICLE )
          write(12, "(6g30.21)") dxdn(:)*(ii-0.5) - size_x(:)/2.0, rho_red(ii,:)/( num_procs*n_moves*dxdn(:)*N_PARTICLE )
        end do
        close(12)
      end if

      if(eval_off_diagonal) then

        call MPI_REDUCE(obdm_x,obdm_tmp,n_bins*n_bins,mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(my_rank .eq. 0) then
          write(my_fname,"(a10)")"obdm_x.dat"
          open(12, file=my_fname, status="replace")
          do ii = 1, N_BINS
            do jj = 1, N_BINS
              write(12, *) dxdn(1)*(ii-0.5) - size_x(1)/2.0, dxdn(1)*(jj-0.5) - size_x(1)/2.0, &
                           obdm_tmp(ii,jj)/( num_procs*2*n_moves*dxdn(1) )
            end do
            write(12, *)
          end do
          close(12)
        end if

        call MPI_REDUCE(obdm_y,obdm_tmp,n_bins*n_bins,mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(my_rank .eq. 0) then
          write(my_fname,"(a10)")"obdm_y.dat"
          open(12, file=my_fname, status="replace")
          do ii = 1, N_BINS
            do jj = 1, N_BINS
              write(12, *) dxdn(2)*(ii-0.5) - size_x(2)/2.0, dxdn(2)*(jj-0.5) - size_x(2)/2.0, &
                           obdm_tmp(ii,jj)/( num_procs*2*n_moves*dxdn(2) )
            end do
            write(12, *)
          end do
          close(12)
        end if

        call MPI_REDUCE(obdm_z,obdm_tmp,n_bins*n_bins,mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        if(my_rank .eq. 0) then
          write(my_fname,"(a10)")"obdm_z.dat"
          open(12, file=my_fname, status="replace")
          do ii = 1, N_BINS
            do jj = 1, N_BINS
              write(12, *) dxdn(3)*(ii-0.5) - size_x(3)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0, &
                           obdm_tmp(ii,jj)/( num_procs*2*n_moves*dxdn(3) )
            end do
            write(12, *)
          end do
          close(12)
        end if

        if(eval_obdm_full) then
          call MPI_REDUCE(obdm_full,obdm_ftmp,(n_bins_full**n_dim)*(n_bins_full**n_dim), &
                          mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if(my_rank .eq. 0) then
            write(my_fname,"(a13)")"obdm_full.dat"
            open(12, file=my_fname, status="replace")
            do ii = 1, N_BINS_full**n_dim
              do jj = 1, N_BINS_full**n_dim
                write(12, *) dxdn(3)*(ii-0.5) - size_x(3)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0, &
                             obdm_ftmp(ii,jj)/( num_procs*2*n_moves*dxdn(3) )
              end do
              write(12, *)
            end do
            close(12)
          end if
        else
          if(eval_obdm_rz) then
            call MPI_REDUCE(obdm_rz,obdm_ftmp,(n_bins_full**2)*(n_bins_full**2),mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(my_rank .eq. 0) then
              write(my_fname,"(a11)")"obdm_rz.dat"
              open(12, file=my_fname, status="replace")
              do ii = 1, N_BINS_full**2
                do jj = 1, N_BINS_full**2
                  write(12, *) tdxdn(1)*(ii-0.5) - tsize(1)/2.0, tdxdn(1)*(jj-0.5) - tsize(1)/2.0, &
                               obdm_ftmp(ii,jj)/( num_procs*2*n_moves*tdxdn(1) )
                end do
                write(12, *)
              end do
              close(12)
            end if
          end if
        end if

        if(eval_obdm_angle_in_XZ_plane) then
          call MPI_REDUCE(obdm_theta,obdm_tmp,n_bins*n_bins,mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if(my_rank .eq. 0) then
            write(my_fname,"(a26)")"obdm_angle_in_XZ_plane.dat"
            open(12, file=my_fname, status="replace")
            do ii = 1, N_BINS
              do jj = 1, N_BINS
                write(12, *) dxdn(3)*(ii-0.5) - size_x(3)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0, &
                             obdm_tmp(ii,jj)/( num_procs*2*n_moves*dxdn(3) )
              end do
              write(12, *)
            end do
            close(12)
          end if
        end if

      end if ! ( eval_off_diagonal )

      write(my_fname,"(a6,i4.4)")"U.dat.",my_rank
      open(12, file=my_fname,  POSITION="APPEND")
      t1 = sum(U_0(CSLICE,:))
      t2 = sum(U_0(CSLICE+1,:))
      write(12, "(2g20.10)")  t1/N_PARTICLE,t2/N_PARTICLE
      close(12)

      call MPI_ALLreduce(e_diff,e_red,1,mpi_double_precision,MPI_SUM,MPI_COMM_WORLD,ierr)
      write(my_fname,"(a11,i4.4)")"Energy.dat.",my_rank
      open(12, file=my_fname,  POSITION="APPEND")
      write(12, "(56g30.20)")  i,E_left/n_moves_pb,E_center/n_moves_pb,E_right/n_moves_pb,dE/n_moves_pb,e_red/(n_moves_pb*num_procs)
      close(12)

      call MPI_REDUCE(rho_r,gofr_red,n_bins_gofr,mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(my_rank .eq. 0) then
        open(12, file="rho_r.dat", status="replace")
        do ii = 1, size(rho_r)
          t1 = drdn*(ii)
          t2 = drdn*(ii-1)
          t3 = 4*M_PI*(t1**3-t2**3)/3
          write(12, "(4g12.3)") drdn*(ii-0.5), gofr_red(ii)/( t3*n_moves*N_PARTICLE*num_procs )
        end do
        close(12)
      end if

      call MPI_REDUCE(gofr,gofr_red,n_bins_gofr,mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(my_rank .eq. 0) then
        open(12, file="gofr.dat", status="replace")
        do ii = 1,size(gofr)
          t0 = dxdn_GOFR*(ii-0.5)
          t1 = dxdn_GOFR*(ii-1)
          t2 = dxdn_GOFR*(ii)
          t3 = (4.*M_PI/3)*(t2**3 - t1**3)*N_PARTICLE*(N_PARTICLE+1)/2
          write(12, "(2g12.3)") t0, gofr_red(ii)/( n_moves*t3*num_procs )
        end do
        close(12)
      end if

      call MPI_REDUCE(gofz,gofr_red,n_bins_gofr,mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(my_rank .eq. 0) then
        open(12, file="gofz.dat", status="replace")
        do ii = 1,size(gofz)
          t0 = dxdn_GOFR*(ii-0.5)
          t1 = dxdn_GOFR*N_PARTICLE*(N_PARTICLE+1)/2
          write(12, "(2g12.3)") t0, gofr_red(ii)/( n_moves*t1*num_procs )
        end do
        close(12)
      end if

      call MPI_REDUCE(gofrho,gofr_red,n_bins_gofr,mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(my_rank .eq. 0) then
        open(12, file="gofrho.dat", status="replace")
        do ii = 1,size(gofr_red)
          t0 = dxdn_GOFR*(ii-0.5)
          t1 = dxdn_GOFR*(ii-1)
          t2 = dxdn_GOFR*(ii)
          t3 = (2.*M_PI)*(t2**2 - t1**2)*N_PARTICLE*(N_PARTICLE+1)/2
          write(12, "(2g12.3)") t0, gofr_red(ii)/( n_moves*t3*num_procs )
        end do
        close(12)
      end if

      call MPI_REDUCE(gofzrho,od2,(n_bins)*(n_bins),mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(my_rank .eq. 0) then
        open(12, file="gofzrho.dat", status="replace")
        do ii = 1,n_bins
          do jj = 1,n_bins
            t0 = dxdn(1)*(jj-0.5)/2
            t1 = dxdn(1)*(jj-1)/2
            t2 = dxdn(1)*(jj)/2
            t3 = (2.*M_PI)*(t2**2 - t1**2)*N_PARTICLE*(N_PARTICLE+1)/2
            write(12, "(3g12.3)")  dxdn(3)*(ii-0.5)/2, t0, od2(ii,jj)/( n_moves*t3*num_procs*dxdn(3) )
          end do
          write(12, *)
        end do
        close(12)
      end if

      if (eval_2particle_angle_correlation) then
          call MPI_REDUCE(goftheta,goftheta_red,n_bins,mpi_integer,MPI_SUM,0,MPI_COMM_WORLD,ierr)
          if(my_rank .eq. 0) then
            open(12, file="goftheta.dat", status="replace")
            do ii = 1,n_bins
              write(12, "(2g12.3)") dthetadn*ii, float(goftheta_red(ii))/( n_moves*num_procs*n_particle )
            end do
            close(12)
          end if
      end if

      if ( eval_rotation ) then
        write(my_fname,"(a12,i4.4)")"gof_rot.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 1,N_BINS_ROT
          write(12, "(2g12.3)")  2.0*(ii-1)/(N_BINS_ROT-1)-1.0, gof_rot(ii)/( n_moves )
        end do
        close(12)
      end if

      if ( eval_rotation ) then
        write(my_fname,"(a13,i4.4)")"gofr_rot.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 1,size(gofr)
          t0 = dxdn_GOFR*(ii-0.5)
          t1 = dxdn_GOFR*(ii-1)
          t2 = dxdn_GOFR*(ii)
          t3 = (4.*M_PI/3)*(t2**3 - t1**3)*N_PARTICLE*(N_PARTICLE+1)/2
          write(12, "(2g12.3)") t0, gofr_rot(ii)/( n_moves*t3 )
        end do
        close(12)
        write(my_fname,"(a16,i4.4)")"gof_rot_xyz.dat.",my_rank
        open(12, file=my_fname, status="replace")
        write(my_fname,"(a20,i4.4)")"gof_rot_xyz_avg.dat.",my_rank
        open(13, file=my_fname, status="replace")
        do ii = 1,N_BINS_ROT_XYZ
          do jj = 1,N_BINS_ROT_XYZ
            do kk = 1,N_BINS_ROT_XYZ
              write(13, "(5g12.3)", ADVANCE="NO") dxdn_rot_xyz(1)*(ii-0.5) - size_x(1)/2.0, &
                                                  dxdn_rot_xyz(2)*(jj-0.5) - size_x(2)/2.0, &
                                                  dxdn_rot_xyz(3)*(kk-0.5) - size_x(3)/2.0, &
                                                  gof_rot_xyz_avg(ii,jj,kk)/( n_moves )
              write(13, *)
              do ll = 1,N_BINS_ROT
                write(12, "(5g12.3)", ADVANCE="NO") dxdn_rot_xyz(1)*(ii-0.5) - size_x(1)/2.0, &
                                                    dxdn_rot_xyz(2)*(jj-0.5) - size_x(2)/2.0, &
                                                    dxdn_rot_xyz(3)*(kk-0.5) - size_x(3)/2.0, &
                                                    2.0*(ll-1)/(N_BINS_ROT-1) - 1.0, &
                                                    gof_rot_xyz(ii,jj,kk,ll)/( n_moves )
                write(12, *)
              end do
              write(12, *)
            end do
            write(13, *)
          end do
        end do
        close(12)
        close(13)
      end if

      if(use_eval_cfn) then
        write(my_fname,"(a11,i4.4)")"corr_x.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 2, N_SLICE-1
          if ( ii .le. CSLICE ) then
            write(12, "(5g24.16)") dtau*(ii-1), corr_x(ii,:)/n_moves, dtau*(ii-1)
          else
            write(12, "(5g24.16)") (ii-1)*dtau, corr_x(ii,:)/n_moves, (N_SLICE-ii)*dtau
          end if
        end do
        close(12)

        write(my_fname,"(a12,i4.4)")"corr_xz.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 2, N_SLICE-1
          if ( ii .le. CSLICE ) then
            write(12, "(4g24.16)") dtau*(ii-1), corr_xz(ii,:)/n_moves, dtau*(ii-1)
          else
            write(12, "(4g24.16)") (ii-1)*dtau, corr_xz(ii,:)/n_moves, (N_SLICE-ii)*dtau
          end if
        end do
        close(12)

        if( use_eval_corr_phase ) then
          write(my_fname,"(a15,i4.4)")"corr_phase.dat.",my_rank
          open(12, file=my_fname, status="replace")
          do ii = 2, N_SLICE-1
            if ( ii .le. CSLICE ) then
              write(12, "(5g24.16)") ii*dtau, corr_phase(ii)/n_moves, ii*dtau
            else
              write(12, "(5g24.16)") ii*dtau, corr_phase(ii)/n_moves, (N_SLICE+1-ii)*dtau
            end if
          end do
          close(12)
        end if 

        write(my_fname,"(a11,i4.4)")"corr_r.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 2, N_SLICE-1
          if ( ii .le. CSLICE ) then
            write(12, "(5g24.16)") dtau*(ii-1), corr_r(ii)/n_moves, dtau*(ii-1)
          else
            write(12, "(5g24.16)") (ii-1)*dtau, corr_r(ii)/n_moves, (N_SLICE-ii)*dtau
          end if
        end do
        close(12)
      end if

      !      write(my_fname,"(a10,i4.4)")"dxlxr.dat.",my_rank
      !      open(12, file=my_fname, status="replace")
      !      do ii = 1,size(dxlxr,1)
      !        write(12, "(6g24.16)", ADVANCE="NO") dxdn*(ii-1) - size_x/2.0, dxlxr(ii,:)/( n_moves*dxdn*N_PARTICLE )
      !        write(12, *)
      !      end do
      !      close(12)

      if(eval_rdm) then
        write(my_fname,"(a11,i4.4)")"rdm2_z.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 1, N_BINS
          do jj = 1, N_BINS
            write(12, *) dxdn(3)*(ii-0.5) - size_x(3)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0,  rdm2_z(ii,jj)/( n_moves*dxdn(3)**2 )
          end do
          write(12, *)
        end do
        close(12)
      end if

      if( eval_nrdm ) then
        write(my_fname,"(a11,i4.4)")"nrdm_z.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 1, 2**N_OD_PARTICLE
          do jj = 1, 2**N_OD_PARTICLE
            write(12, *) ii,jj,nrdm_z(ii,jj)/( n_moves )
          end do
          write(12, *)
        end do
        close(12)
      end if

      if( eval_off_diagonal ) then
        if( eval_rotation ) then
          write(my_fname,"(a15,i4.4)")"obdm_rot_z.dat.",my_rank
          open(12, file=my_fname, status="replace")
          do ii = 1, N_BINS
            do jj = 1, N_BINS
              write(12, *) 2.0_b8*(ii-0.5)/N_BINS - 1.0_b8, 2.0_b8*(jj-0.5)/N_BINS - 1.0_b8,  obdm_rot_z(ii,jj)/( n_moves)
            end do
            write(12, *)
          end do
          close(12)
        end if
      end if

      if( eval_off_diagonal .and. eval_obdm_angle_in_XZ_plane) then
        write(my_fname,"(a27,i4.4)")"obdm_angle_in_XZ_plane.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 1, N_BINS
          do jj = 1, N_BINS
            write(12, *) (ii-0.5)/(N_BINS), (jj-0.5)/(N_BINS),  obdm_theta(ii,jj)/( 2*n_moves )
          end do
          write(12, *)
        end do
        close(12)
      end if

      if(eval_rdm2_ring) then
        write(my_fname,"(a14,i4.4)")"rdm2_ring.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 1, N_BINS
          do jj = 1, N_BINS
            write(12, *) (ii-0.5)/(N_BINS), (jj-0.5)/(N_BINS),  rdm2_theta(ii,jj)/( 2*n_moves )
          end do
          write(12, *)
        end do
        close(12)

        write(my_fname,"(a15,i4.4)")"rdm22_ring.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 1, N_BINS
          do jj = 1, N_BINS
            write(12, *) (ii-0.5)/(N_BINS), (jj-0.5)/(N_BINS),  rdm22_theta(ii,jj)/( 2*n_moves )
          end do
          write(12, *)
        end do
        close(12)
      end if

      write(my_fname,"(a11,i4.4)")"mcount.dat.",my_rank
      open(12, file=my_fname, status="replace")
      do ii = 1, N_SLICE
        write(12, *) ii,mcount(ii)/n_moves, mcount(ii)
      end do
      close(12)

      write(my_fname,"(a11,i4.4)")"accr_v.dat.",my_rank
      open(12, file=my_fname, status="replace")
      do ii = 1, N_SLICE
        write(12, *) ii,accr_v(ii)/mcount(ii)
      end do
      close(12)

      if( use_eval_N_R ) then
        call MPI_REDUCE(dn_rms,dn_rms_red,1,mpi_double_precision,MPI_SUM,0,MPI_COMM_WORLD,ierr)

        if(my_rank .eq. 0) then
          open(12, file="dn_rms.dat",  POSITION="APPEND")
          write(12, "(5g30.15)")  i,sqrt(dn_rms_red / (n_moves_pb*num_procs))
          close(12)
        end if

        write(my_fname,"(a8,i4.4)")"N_R.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 1, N_PARTICLE+1
          write(12, *) (ii-1), N_R(ii)/( n_moves )
        end do
        close(12)

        write(my_fname,"(a9,i4.4)")"N_zs.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 1, N_PARTICLE2+1
          write(12, *) (ii-1), N_zs(ii)/( n_moves )
        end do
        close(12)

        write(my_fname,"(a9,i4.4)")"N_za.dat.",my_rank
        open(12, file=my_fname, status="replace")
        do ii = 1, N_PARTICLE2+1
          write(12, *) (ii-1), N_za(ii)/( n_moves )
        end do
        close(12)
      end if 

      !      write(my_fname,"(a13,i4.4)")("corr_N_R.dat.",my_rank)
      !      open(12, file=my_fname, status="replace")
      !      do ii = 1, N_PARTICLE+1
      !        write(12, "(i12)", ADVANCE="NO") (ii-1)
      !        do jj = 1, N_SLICE
      !          write(12, "(g20.12)", ADVANCE="NO") corr_N_R(jj,ii)/( n_moves )
      !        end do
      !        write(12, *)
      !      end do
      !      close(12)


      if(eval_full_density) then
        write(my_fname,"(a9,i4.4)")"odxyz.dat.",my_rank
        open(23, file=my_fname, status="replace")
        do ii = 1, N_BINS
          do jj = 1, N_BINS
            do kk = 1, N_BINS
              write(23, *) dxdn(1)*(ii-0.5) - size_x(1)/2.0, &
                           dxdn(2)*(jj-0.5) - size_x(2)/2.0, &
                           dxdn(3)*(jj-0.5) - size_x(3)/2.0, &
                           full_density(ii,jj,kk)/( n_moves*N_PARTICLE*dxdn(1)*dxdn(2)*dxdn(3) )
            end do
            write(23, *)
          end do
          write(23, *)
        end do
        close(23)
      end if

      if(eval_column_density) then
        write(my_fname,"(a9,i4.4)")"odxy.dat.",my_rank
        open(23, file=my_fname, status="replace")
        do ii = 1, N_BINS
          do jj = 1, N_BINS
            write(23, *) dxdn(1)*(ii-0.5) - size_x(1)/2.0, dxdn(2)*(jj-0.5) - size_x(2)/2.0, odxy_avg(ii,jj)
          end do
          write(23, *)
        end do
        close(23)

        write(my_fname,"(a9,i4.4)")"odxz.dat.",my_rank
        open(23, file=my_fname, status="replace")
        do ii = 1, N_BINS
          do jj = 1, N_BINS
            write(23, *) dxdn(1)*(ii-0.5) - size_x(1)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0, odxz_avg(ii,jj)
          end do
          write(23, *)
        end do
        close(23)

        write(my_fname,"(a9,i4.4)")"odyz.dat.",my_rank
        open(23, file=my_fname, status="replace")
        do ii = 1, N_BINS
          do jj = 1, N_BINS
            write(23, *) dxdn(2)*(ii-0.5) - size_x(2)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0, odyz_avg(ii,jj)
          end do
          write(23, *)
        end do
        close(23)

        write(my_fname,"(a8,i4.4)")"od2.dat.",my_rank
        open(23, file=my_fname, status="replace")
        do ii = 1, N_BINS
          do jj = 1, N_BINS
            write(23, *) dxdn(1)*(ii-0.5) - size_x(1)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0, od2_avg(ii,jj)
          end do
          write(23, *)
        end do
        close(23)
      end if

      !      write(my_fname,"(a10,i4.4)")"od_sig.dat.",my_rank
      !      open(24, file=my_fname, status="replace")
      !      do ii = 1, N_BINS
      !        do jj = 1, N_BINS
      !          write(24, *) dxdn(1)*(ii-0.5) - size_x(1)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0, od_sig(ii,jj)/( n_moves*N_PARTICLE )
      !        end do
      !        write(24, *)
      !      end do
      !      close(24)

      if(eval_qsq_sum) then
        write(my_fname,"(a8,i4.4)")"qsq.dat.",my_rank
        open(10, file=my_fname, status="replace")
        write(10, *) "#", accr/(n_moves)
        do ii = 1,N_SLICE_SAMPLES
      !          ik =  int((ii-1)*float(N_SLICE)/(N_SLICE_SAMPLES-1)) + 1
          ik = ii
          write(10, "(4g20.12)") (ik-1)*dtau, qsq_sum(ii,:)/(N_PARTICLE*n_moves)
        end do
        close(10)

        write(my_fname,"(a10,i4.4)")"U_avg.dat.",my_rank
        open(12, file=my_fname,  status="replace")
        do ii = 1, size(U_avg)
      !          ik =  int((ii-1)*float(N_SLICE)/(N_SLICE_SAMPLES-1)) + 1
          ik = ii
          write(12, "(2g30.15)")  (ik-1)*dtau, U_avg(ii)
        end do
        close(12)

      end if

      if(use_eval_phase) then
        write(my_fname,"(a10,i4.4)")"phase.dat.",my_rank
        open(30, file=my_fname, status="replace")
        write(30, *) "#", accr/(n_moves)
        do ii = 1,N_SLICE
          write(30, "(i12,2g20.12)") ii, phase_0(ii), phase_sum_0(ii)/mcount(ii)
        end do
        close(30)
      end if

      write(my_fname,"(a14,i4.4)")"last_path.dat.",my_rank
      open(10, file=my_fname, status="replace")
      do ii = 1,N_SLICE
        write(10, *) q0(ii,:,:)
      end do
      close(10)

      if( eval_rotation ) then
        write(my_fname,"(a18,i4.4)")"last_rot_path.dat.",my_rank
        open(10, file=my_fname, status="replace")
        do ii = 1,N_SLICE
          write(10, *) q_rot0(ii,:,:)
        end do
        close(10)
      end if
      if(write_paths) then
        write(my_fname,"(a10,i4.4)")"paths.dat.",my_rank
        open(10, file=my_fname, status="replace")
        do ii = 1,CSLICE
          write(10, "(i12)", ADVANCE="NO") ii
          jj = 1
          do jj = 1,N_PARTICLE
            write(10, "(3g20.12)", ADVANCE="NO") q0(ii,jj,1),q0(ii,jj,2),q0(ii,jj,3)
          end do
          write(10, *) 
        end do
        write(10, *) 
        do ii = CSLICE+1, N_SLICE
          write(10, "(i12)", ADVANCE="NO") ii
          do jj = 1,N_PARTICLE
            write(10, "(3g20.12)", ADVANCE="NO") q0(ii,jj,1),q0(ii,jj,2),q0(ii,jj,3)
          end do
          write(10, *) 
        end do
        close(10)

        write(my_fname,"(a14,i4.4)")"rot_paths.dat.",my_rank
        open(10, file=my_fname, status="replace")
        do ii = 1,CSLICE
          write(10, "(i12,3g20.12)") ii, q_rot0(ii,od_pnum,1),q_rot0(ii,od_pnum,2),q_rot0(ii,od_pnum,3)
        end do
          write(10, *) 
        do ii = CSLICE+1, N_SLICE
          write(10, "(i12,3g20.12)") ii-1, q_rot0(ii,od_pnum,1),q_rot0(ii,od_pnum,2),q_rot0(ii,od_pnum,3)
        end do
        close(10)
      end if

      write(my_fname,"(a11,i4.4)")"trho_L.dat.",my_rank
      open(10, file=my_fname, status="replace")
      write(my_fname,"(a11,i4.4)")"trho_R.dat.",my_rank
      open(11, file=my_fname, status="replace")
      do ii = 1,size(trial_rho_L,1)
        write(10, "(6g12.3)") dxdn*(ii-0.5)-size_x/2.0, trial_rho_L(ii,:)/( n_moves*N_PARTICLE*dxdn )
        write(11, "(6g12.3)") dxdn*(ii-0.5)-size_x/2.0, trial_rho_R(ii,:)/( n_moves*N_PARTICLE*dxdn )
      end do
      close(10)
      close(11)

      !@-node:gcross.20090626091217.1690:<< Compute and output current simulation results >>
      !@nl

      !@  << Stop the timer and output iteration statistics for this block >>
      !@+node:gcross.20090626112946.1687:<< Stop the timer and output iteration statistics for this block >>
      !@@raw
#ifdef USE_MPI
      !@@end_raw
        time1 = MPI_Wtime()
      !@@raw
#else
      !@@end_raw
        time1 = day_timer()
      !@@raw
#endif
      !@@end_raw

      write(my_fname,"(a10,i4.4)")"stats.dat.",my_rank
      open(27, file=my_fname,  POSITION="APPEND")
      write(27,"(a10, 7a12)") "#block ","dtime ","bi_accr ", "sp_accr ", "swap accr", "coll count"
      if(bi_moves .eq. 0) bi_moves = 1
      if(sp_moves .eq. 0) sp_moves = 1
      if(swap_moves .eq. 0) swap_moves = 1
      write(27,"(i10, 5g12.3)")  i, time1-time0, &
        bi_accr/bi_moves, sp_accr/sp_moves, swap_accr/swap_moves,ccnt/n_moves
      !     write(27,"(a4, 2i12)") "#", left_moves, right_moves
      close(27)

      inquire(file="terminate_all", exist=ex)
      if ( ex ) then
        goto 13
      end if
      !@-node:gcross.20090626112946.1687:<< Stop the timer and output iteration statistics for this block >>
      !@nl

    end do
    !@-node:gcross.20090623152316.33:<< Main iteration >>
    !@nl

!    call nl2sno_driver()

    ! Write random number seed
13  call ru_write_seed_file(my_rank)
#ifdef USE_MPI
    call MPI_Finalize(ierr)
#endif

  contains

!@+others
!@+node:gcross.20090623152316.13:Subroutines
!@+others
!@+node:gcross.20090624144408.2048:Evaluators of physical quantities
!@+node:gcross.20090626170642.1720:Number density
!@+node:gcross.20090623152316.15:vpi_eval_density
!@+at
! For each coordinate axis in position space, vpi_eval_density evaluates the 
! number density as a function of the coordinate on that axis, integrating 
! over the other axes.
! 
! That is to say, rho(:,n) is the number density as a function of the n-th 
! coordinate, integrating over the other coordinates.
!@-at
!@@c
subroutine vpi_eval_density( rho, x, nbins, size_x, dndx )
  real(kind=b8), dimension( nbins , N_DIM ), intent(out):: rho
  real(kind=b8), dimension( N_PARTICLE , N_DIM ) :: x
  integer nbins
  real, dimension(N_DIM) :: size_x, dndx

  integer bin
  integer :: i,j

  do i = 1, N_PARTICLE
    do j = 1, N_DIM
      bin = floor((x(i,j)+size_x(j))*dndx(j))+1
      if ( within_bins(bin,nbins) ) then
        rho(bin,j) = rho(bin,j) + 1
      end if
    end do
  end do
end subroutine vpi_eval_density
!@-node:gcross.20090623152316.15:vpi_eval_density
!@+node:gcross.20090623152316.18:vpi_eval_density_single_particle
!@+at
! Same as vpi_eval_density, but only computes the density for a single given 
! particle.
!@-at
!@@c
subroutine vpi_eval_density_single_particle( rho, x, nbins, size_x, dndx )
  real(kind=b8), dimension( nbins , N_DIM ), intent(out):: rho
  real(kind=b8), dimension( N_DIM ) :: x
  integer nbins
  real, dimension(N_DIM) :: size_x, dndx
  integer bin

  integer :: j

  do j = 1, N_DIM
    bin = anint((x(j)+size_x(j))*dndx(j))+1
    if ( within_bins(bin,nbins) ) then
      rho(bin,j) = rho(bin,j) + 1
    end if 
  end do
end subroutine vpi_eval_density_single_particle
!@nonl
!@-node:gcross.20090623152316.18:vpi_eval_density_single_particle
!@+node:gcross.20090626170642.1721:vpi_eval_full_density
!@+at
! Evaluates the number density as a function of the first three position-space 
! coordinates, integrating over the others.
!@-at
!@@c
subroutine vpi_eval_full_density( od, x, nbins, size_x, dndx )
  real(kind=b8), dimension( : , : , : ), intent(out):: od
  real(kind=b8), dimension( : , : ) :: x
  integer :: nbins
  real, dimension(N_DIM) :: size_x, dndx

  integer :: xb, yb, zb

  integer :: i

  do i = 1, N_PARTICLE
    xb = anint((x(i,1)+size_x(1))*dndx(1))+1
    if ( within_bins(xb,nbins) ) then
      yb = anint((x(i,2)+size_x(2))*dndx(2))+1
      if ( within_bins(yb,nbins) ) then
        zb = anint((x(i,3)+size_x(3))*dndx(3))+1
        if ( within_bins(zb,nbins) ) then
          od(xb,yb,zb) = od(xb,yb,zb) + 1
        end if
      end if
    end if 
  end do
end subroutine vpi_eval_full_density
!@-node:gcross.20090626170642.1721:vpi_eval_full_density
!@+node:gcross.20090623152316.17:vpi_eval_radial_density
!@+at
! Evaluates number density of particles as a function of the radius.
!@-at
!@@c
subroutine vpi_eval_radial_density( rho, x, nbins, dndx )
  real(kind=b8), dimension( : ), intent(out):: rho
  real(kind=b8), dimension( : , : ) :: x
  integer nbins
  real :: dndx

  real(kind=b8) :: r
  integer bin

  integer :: i,j

  do i = 1, size(x,1)
    r = sqrt(dot_product(x(i,:),x(i,:)))
    bin = floor(r*dndx)+1
    if ( within_bins(bin,nbins) ) then
      rho(bin) = rho(bin) + 1
    end if 
  end do
end subroutine vpi_eval_radial_density
!@-node:gcross.20090623152316.17:vpi_eval_radial_density
!@+node:gcross.20090623152316.14:vpi_eval_density_in_chebyshev_basis
!@+at
! Calculate number density independently along each dimension, like 
! vpi_eval_density, but the result (in each dimension) is expressed in a 
! Chebychev polynomial basis.
!@-at
!@@c
subroutine vpi_eval_density_in_chebyshev_basis( rho, x, nbins, size_x, dndx )
  real(kind=b8), dimension( nbins , N_DIM ), intent(out):: rho
  real(kind=b8), dimension( N_PARTICLE , N_DIM ) :: x
  integer nbins
  real, dimension(N_DIM) :: size_x, dndx

  real, dimension(nbins) :: yy

  integer bin
  integer :: i,j,n


  do i = 1, N_PARTICLE
   do j = 1, N_DIM
     yy(1) = 1
     yy(2) = x(i,j)
     do n = 3, nbins
       yy(n) = 2.0*x(i,j)*yy(n-1)/size_x(j) - yy(n-2)
     end do
     rho(:,j) = rho(:,j) + yy(:)
   end do
 end do
end subroutine vpi_eval_density_in_chebyshev_basis
!@nonl
!@-node:gcross.20090623152316.14:vpi_eval_density_in_chebyshev_basis
!@+node:gcross.20090623152316.19:vpi_eval_density_in_plane
!@+at
! Calculate number density as a function of position in a plane spanning the 
! "id1" and "id2" axes, integrating over any remaining axes.
!@-at
!@@c
subroutine vpi_eval_density_in_plane( od, od_avg, od_sig, x, id1, id2, nbins, size_x, dndx )
  real(kind=b8), dimension( nbins , nbins ), intent(out):: od
  real(kind=b8), dimension( nbins , nbins ), intent(in):: od_avg
  real(kind=b8), dimension( nbins , nbins ), intent(out):: od_sig
  real(kind=b8), dimension( N_PARTICLE , N_DIM ) :: x
  real(kind=b8), dimension( nbins , nbins ) :: tod
  integer nbins, id1, id2
  real, dimension(N_DIM) :: size_x, dndx
  integer xbin, zbin

  integer :: i

  tod = 0

  do i = 1, N_PARTICLE-N_PARTICLE2
    xbin = anint((x(i,id1)+size_x(id1))*dndx(id1))+1
    if ( within_bins(xbin,nbins) ) then
      zbin = anint((x(i,id2)+size_x(id2))*dndx(id2))+1
      if ( within_bins(zbin,nbins) ) then
        tod(xbin,zbin) = tod(xbin,zbin) + 1
      end if
    end if 
  end do
  od = od + tod
  od_sig(:,:) = od_sig(:,:) + (tod(:,:) - od_avg(:,:))**2
end subroutine vpi_eval_density_in_plane
!@nonl
!@-node:gcross.20090623152316.19:vpi_eval_density_in_plane
!@+node:gcross.20090623152316.20:vpi_eval_density_in_plane2
!@+at
! Calculate number density as a function of position in a plane spanning the X 
! and Z axes, integrating over any remaining axes, for the second particle 
! species?
!@-at
!@@c
subroutine vpi_eval_density_in_plane2( od, od_avg, od_sig, x, nbins, size_x, dndx )
  real(kind=b8), dimension( nbins , nbins ), intent(out):: od
  real(kind=b8), dimension( nbins , nbins ), intent(in):: od_avg
  real(kind=b8), dimension( nbins , nbins ), intent(out):: od_sig
  real(kind=b8), dimension( N_PARTICLE , N_DIM ) :: x
  real(kind=b8), dimension( nbins , nbins ) :: tod
  integer nbins
  real, dimension(N_DIM) :: size_x, dndx
  integer xbin, zbin

  integer :: i

  tod = 0

  do i = N_PARTICLE - N_PARTICLE2 + 1, N_PARTICLE
    xbin = anint((x(i,1)+size_x(1))*dndx(1))+1
    if ( within_bins(xbin,nbins) ) then
      zbin = anint((x(i,3)+size_x(3))*dndx(3))+1
      if ( within_bins(zbin,nbins) ) then
        tod(xbin,zbin) = tod(xbin,zbin) + 1
      end if
    end if 
  end do
  od = od + tod
  od_sig(:,:) = od_sig(:,:) + (tod(:,:) - od_avg(:,:))**2
end subroutine vpi_eval_density_in_plane2
!@nonl
!@-node:gcross.20090623152316.20:vpi_eval_density_in_plane2
!@-node:gcross.20090626170642.1720:Number density
!@+node:gcross.20090626170642.1727:Correlators
!@+node:gcross.20090623152316.23:vpi_eval_corr_x
subroutine vpi_eval_corr_x( corr_x, move_start, move_end, x )
  real(kind=b8), dimension( N_SLICE, N_DIM ), intent(inout):: corr_x 
  integer :: move_start, move_end
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x

  integer :: i,j

  do i = move_start+1, move_end-1
    do j = 1, N_DIM
      corr_x(i,j) = corr_x(i,j) + sum( x(1,:,j)*x(i,:,j) )
      corr_x(N_SLICE-i+1,j) = corr_x(N_SLICE-i+1,j) + sum( x(N_SLICE,:,j)*x(N_SLICE-i+1,:,j) )
    end do
  end do

end subroutine vpi_eval_corr_x
!@-node:gcross.20090623152316.23:vpi_eval_corr_x
!@+node:gcross.20090623152316.24:vpi_eval_corr_phase
subroutine vpi_eval_corr_phase( corr_phase, phase, move_start, move_end )
  real(kind=b8), dimension( N_SLICE ), intent(inout):: corr_phase
  real(kind=b8), dimension( N_SLICE ), intent(in) :: phase
  integer :: move_start, move_end

  integer :: i

  do i = move_start+1, move_end-1
    corr_phase(i) = corr_phase(i) + phase(1)*phase(i)
    corr_phase(N_SLICE-i+1) = corr_phase(N_SLICE-i+1) + phase(N_SLICE)*phase(N_SLICE-i+1)
  end do

end subroutine vpi_eval_corr_phase
!@-node:gcross.20090623152316.24:vpi_eval_corr_phase
!@+node:gcross.20090623152316.25:vpi_eval_corr_r
subroutine vpi_eval_corr_r( corr_r, move_start, move_end, x )
  real(kind=b8), dimension( N_SLICE ), intent(inout):: corr_r
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x
  integer :: move_start, move_end

  integer :: i
  real(kind=b8) :: r1,rM,tmp

  r1 = sum(x(1,:,1:3)**2)
  rM = sum(x(N_SLICE,:,1:3)**2)
  do i = move_start+1, move_end-1
    tmp = sum(x(i,:,1:3)**2)
    corr_r(i) = corr_r(i) + r1*tmp
    tmp = sum(x(N_SLICE-i+1,:,1:3)**2)
    corr_r(N_SLICE-i+1) = corr_r(N_SLICE-i+1) + rM*tmp
  end do

end subroutine vpi_eval_corr_r
!@-node:gcross.20090623152316.25:vpi_eval_corr_r
!@+node:gcross.20090623152316.26:vpi_eval_corr_xz
subroutine vpi_eval_corr_xz( corr_xz, move_start, move_end, x )
  real(kind=b8), dimension( N_SLICE, 2 ), intent(inout):: corr_xz
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x
  integer :: move_start, move_end

  integer :: i,j
  real(kind=b8) :: r1,rM,tmp

  do i = move_start+1, move_end-1
    do j = 1, 2
      corr_xz(i,j) = corr_xz(i,j) + sum( x(1,:,j)*x(1,:,3)*x(i,:,j)*x(i,:,3) )
      corr_xz(N_SLICE-i+1,j) = corr_xz(N_SLICE-i+1,j) + sum( x(N_SLICE,:,j)*x(N_SLICE,:,3)*x(N_SLICE-i+1,:,j)*x(N_SLICE-i+1,:,3) )
    end do
  end do

end subroutine vpi_eval_corr_xz
!@-node:gcross.20090623152316.26:vpi_eval_corr_xz
!@-node:gcross.20090626170642.1727:Correlators
!@+node:gcross.20090706131953.1744:Green's functions
!@+at
! These functions compute distributions that relate to the correlations 
! between particles.
!@-at
!@@c
!@+node:gcross.20090623152316.27:vpi_eval_gfn_sep_in_full_space
!@+at
! Computes the number density as a function of the distance between 
! particles;  i.e. given any particle this distribution gives information 
! about how many particles you'd expect to see as a function of distance from 
! said particle.
! 
! The "in_full_space" suffix indicates that this function considers the 
! distance in the full N-dimensional physical space, in contrast to many 
! others which consider the distance restricted to various axes and planes.
!@-at
!@@c
subroutine vpi_eval_gfn_sep_in_full_space( gofr, xij2, nbins, dndx )
  integer :: nbins
  real :: dndx
  real(kind=b8), dimension( nbins ), intent(inout):: gofr 
  real(kind=b8), dimension( N_PARTICLE , N_PARTICLE ) :: xij2

  integer :: i,j
  integer :: rbin
  real :: r

  do i = 1, N_PARTICLE
    do j = i + 1, N_PARTICLE
      r = sqrt(xij2(i,j))
      rbin = floor( r * dndx ) + 1
      if( within_bins(rbin,nbins) ) then
        gofr(rbin) = gofr(rbin) + 1
      end if
    end do
  end do

end subroutine vpi_eval_gfn_sep_in_full_space
!@nonl
!@-node:gcross.20090623152316.27:vpi_eval_gfn_sep_in_full_space
!@+node:gcross.20090623152316.30:vpi_eval_gfn_sep_in_XY_plane
subroutine vpi_eval_gfn_sep_in_XY_plane( gofr, x, islice, nbins, dndx )
  integer :: nbins, islice
  real :: dndx
  real(kind=b8), dimension( nbins ), intent(inout):: gofr 
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x

  integer :: i,j
  integer :: rbin
  real :: r

  real(kind=b8), dimension ( 2 ) :: dxij

  do i = 1, N_PARTICLE
    do j = i + 1, N_PARTICLE
      dxij(:) =  x(islice,i,1:2) - x(islice,j,1:2)
      if (use_pbc) then
        dxij(:) = dxij(:) - p_pbc_L*floor(dxij(:)/p_pbc_L - 0.5_b8) - p_pbc_L
      end if
      r = sqrt(sum(dxij(1:2)**2))
      rbin = floor( r * dndx ) + 1
      if( within_bins(rbin,nbins) ) then
        gofr(rbin) = gofr(rbin) + 1
      end if
    end do
  end do

end subroutine vpi_eval_gfn_sep_in_XY_plane
!@-node:gcross.20090623152316.30:vpi_eval_gfn_sep_in_XY_plane
!@+node:gcross.20090623152316.29:vpi_eval_gfn_sep_along_Z_axis
!@+at
! Computes the number density as a function of the separation between 
! particles along the Z axis;  i.e., like vpi_eval_gfn_sep_in_full_space, but 
! expressed as a separate function of the distance along each axis rather than 
! as the distance in the full N-dimensional space.
!@-at
!@@c
subroutine vpi_eval_gfn_sep_along_Z_axis( gofr, x, islice, nbins, dndx )
  integer :: nbins, islice
  real :: dndx
  real(kind=b8), dimension( nbins ), intent(inout):: gofr 
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x

  integer :: i,j
  integer :: rbin
  real :: r

  do i = 1, N_PARTICLE
    do j = i + 1, N_PARTICLE
      r = abs(x(islice,i,3) - x(islice,j,3))
      rbin = floor( r * dndx ) + 1
      if( within_bins(rbin,nbins) ) then
        gofr(rbin) = gofr(rbin) + 1
      end if
    end do
  end do

end subroutine vpi_eval_gfn_sep_along_Z_axis
!@nonl
!@-node:gcross.20090623152316.29:vpi_eval_gfn_sep_along_Z_axis
!@+node:gcross.20090623152316.28:vpi_eval_gfn_sep_in_Z_and_XY
subroutine vpi_eval_gfn_sep_in_Z_and_XY( gofr, x, islice, nbins, dndx )
  integer :: nbins, islice
  real :: dndx(3)
  real(kind=b8), dimension( nbins,nbins ), intent(inout):: gofr 
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x

  integer :: i,j
  integer :: zbin,rbin
  real :: z,r

  real(kind=b8), dimension ( 3 ) :: dxij

  do i = 1, N_PARTICLE
    do j = i + 1, N_PARTICLE
      dxij(:) =  x(islice,i,:) - x(islice,j,:)
      if(use_pbc) then
        dxij(:) = dxij(:) - p_pbc_L*floor(dxij(:)/p_pbc_L - 0.5_b8) - p_pbc_L
      end if
      ! separation in the XY plane
      r = sqrt(sum(dxij(1:2)**2))
      rbin = floor( r * dndx(1)*2 ) + 1
      ! separation along the Z axis
      z = abs(dxij(3))*2
      zbin = floor( z * dndx(3) ) + 1
      if( within_bins(rbin,nbins) .and. within_bins(zbin,nbins) ) then
        gofr(zbin,rbin) = gofr(zbin,rbin) + 1
      end if
    end do
  end do

end subroutine vpi_eval_gfn_sep_in_Z_and_XY
!@nonl
!@-node:gcross.20090623152316.28:vpi_eval_gfn_sep_in_Z_and_XY
!@+node:gcross.20090721121051.1759:vpi_eval_gfn_sep_angle
subroutine vpi_eval_gfn_sep_angle( goftheta, x, islice, nbins, dndtheta)
  integer, intent(in) :: nbins, islice
  integer, dimension( nbins ), intent(inout) :: goftheta
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ), intent(in) :: x
  real, intent(in) :: dndtheta

  integer :: plane_axis_1, plane_axis_2
  integer :: zbin,rbin
  real :: z,r

  real(kind=b8), dimension( n_particle ) :: theta

  real(kind=b8) :: dtheta
  integer :: dthetabin

  integer :: i, j

  call get_rotation_plane_axes(fixed_rotation_axis,plane_axis_1,plane_axis_2)

  theta(:) = atan2(x(islice,:,plane_axis_2),x(islice,:,plane_axis_1))

  do i = 1, N_PARTICLE
    do j = i + 1, N_PARTICLE
      dtheta = abs(theta(i)-theta(j))
      if(dtheta > m_pi) then
        dtheta = 2*m_pi - dtheta
      end if
      dthetabin = floor( dtheta * dndtheta ) + 1
      goftheta(dthetabin) = goftheta(dthetabin) + 1
    end do
  end do

end subroutine vpi_eval_gfn_sep_angle
!@-node:gcross.20090721121051.1759:vpi_eval_gfn_sep_angle
!@-node:gcross.20090706131953.1744:Green's functions
!@+node:gcross.20090706131953.1750:Off-diagonal elements
!@+node:gcross.20090706131953.1751:vpi_eval_obdm_single_axis
!@+at
! Evaluates the single-particle reduced density matrix (not to be confused 
! with the number density) along a single axis in position space.
!@-at
!@@c
subroutine vpi_eval_obdm_single_axis( obdm, x1, x2, nbins, ndim, xsize, dndx)
  integer :: id, nbins, ndim
  real(kind=b8), dimension ( nbins, nbins ) :: obdm
  real(kind=b8) :: x1, x2
  real :: xsize, dndx

  integer :: bx1, bx2

!  print *, xsize, dndx
!  print *, x1, x2
  bx1 = floor((x1 + xsize)*dndx)+1
  bx2 = floor((x2 + xsize)*dndx)+1
!  print *, bx1,bx2
  if ( within_bins(bx1,nbins) .and. within_bins(bx2,nbins) ) then
    obdm(bx1,bx2) = obdm(bx1,bx2) + 1
  end if

end subroutine vpi_eval_obdm_single_axis
!@-node:gcross.20090706131953.1751:vpi_eval_obdm_single_axis
!@+node:gcross.20090706131953.1752:vpi_eval_obdm_angle_in_XZ_plane
!@+at
! Evaluates the single-particle reduced density matrix (not to be confused 
! with the number density) for the angle in the XZ plane.
!@-at
!@@c
subroutine vpi_eval_obdm_angle_in_XZ_plane( obdm, x1, x2, nbins, ndim)
  integer :: id, nbins, ndim
  real(kind=b8), dimension ( nbins, nbins ) :: obdm
  real(kind=b8), dimension( 3 ) :: x1, x2
  real :: xsize, dndx
  real :: t1,t2

  integer :: bx1, bx2

!  print *, xsize, dndx
!  print *, x1, x2
  if( x1(1) .gt. 0 ) then
    t1 = (atan(x1(3)/x1(1))/M_2PI) + 0.25
  else
    t1 = (atan(x1(3)/x1(1))/M_2PI) + 0.75
  end if

  if( x2(1) .gt. 0 ) then
    t2 = (atan(x2(3)/x2(1))/M_2PI) + 0.25
  else
    t2 = (atan(x2(3)/x2(1))/M_2PI) + 0.75
  end if

  bx1 = floor(t1*nbins)+1
  bx2 = floor(t2*nbins)+1
!  print *, bx1,bx2
  if ( within_bins(bx1,nbins) .and. within_bins(bx2,nbins) ) then
    obdm(bx1,bx2) = obdm(bx1,bx2) + 1
  end if

end subroutine vpi_eval_obdm_angle_in_XZ_plane
!@-node:gcross.20090706131953.1752:vpi_eval_obdm_angle_in_XZ_plane
!@+node:gcross.20090706131953.1753:vpi_eval_obdm_full
!@+at
! Compute the full N-dimensional single-particle reduced density matrix, 
! flattening the multi-dimensional coordinates into a single index so that the 
! result can be expressed as a matrix instead of a 2N-dimensional tensor.
!@-at
!@@c
subroutine vpi_eval_obdm_full( obdmf, x1, x2, nbins, ndim, xsize, dndx)
  integer :: nbins, ndim
  real(kind=b8), dimension ( nbins**ndim, nbins**ndim ) :: obdmf
  real(kind=b8), dimension( ndim ) :: x1, x2
  real,dimension(ndim) :: xsize, dndx

  integer, dimension(ndim)  :: bx1, bx2
  integer :: ii,i1,i2

  bx1(:) = floor((x1(:) + xsize(:))*dndx)+1
  bx2(:) = floor((x2(:) + xsize(:))*dndx)+1
  i1 = 0
  i2 = 0
  do ii=1,ndim
    i1 = i1 + bx1(ii)*nbins**(ii-1)
    i2 = i2 + bx2(ii)*nbins**(ii-1)
  enddo
  print *, bx1(:)
  print *, bx2(:)
  print *, i1, i2
  if ( within_bins(i1,nbins**ndim) .and. within_bins(i2,nbins**ndim) ) then
    obdmf(i1,i2) = obdmf(i1,i2) + 1
  end if

end subroutine vpi_eval_obdm_full
!@-node:gcross.20090706131953.1753:vpi_eval_obdm_full
!@-node:gcross.20090706131953.1750:Off-diagonal elements
!@+node:gcross.20090626170642.1722:No idea what these do yet
!@+node:gcross.20090623152316.31:vpi_eval_gof_rot
subroutine vpi_eval_gof_rot( gof_rot, gof_rot_xyz, gof_rot_xyz_avg, gofr_rot, q, q_rot, xij2, &
                             slice, nbins_rr, dndx_rr, nbins_xyz, xsize, dndx, nbins_rot, dndx_rot )
  integer :: nbins_xyz, nbins_rot, nbins_rr
  integer :: slice
  real(kind=b8), dimension( : ), intent(inout):: gof_rot 
  real(kind=b8), dimension( : ), intent(inout):: gofr_rot 
  real(kind=b8), dimension( : , : , : , : ), intent(inout):: gof_rot_xyz
  real(kind=b8), dimension( : , : , : ), intent(inout):: gof_rot_xyz_avg
  real(kind=b8), dimension( : , : , : ) :: q
  real(kind=b8), dimension( : , : , : ) :: q_rot
  real(kind=b8), dimension( N_SLICE, N_PARTICLE , N_PARTICLE ) :: xij2
  real, dimension(3) :: xsize, dndx
  real :: dndx_rot,dndx_rr

  integer :: i,j
  integer :: br, brr, bx, by, bz
  real(kind=b8) :: r,rr

  do i = 1, N_PARTICLE
    do j = i + 1, N_PARTICLE
      rr = sqrt(xij2(slice,i,j))
      brr = floor( rr * dndx_rot ) + 1
      if((brr .ge. 1) .and. (brr .le. nbins_rot) ) then
        gofr_rot(brr) = gofr_rot(brr) + r
      end if
      bx = floor((q(slice,i,1) + xsize(1))*dndx(1))+1
      by = floor((q(slice,i,2) + xsize(2))*dndx(2))+1
      bz = floor((q(slice,i,3) + xsize(3))*dndx(3))+1
      r = dot_product(q_rot(slice,i,:),q_rot(slice,j,:))
      br = floor((r + 1.0)*nbins_rot/2.0)+1
      if( (br .ge. 1) .and. (br .le. nbins_rot) ) then
        gof_rot(br) = gof_rot(br) + 1
        if( within_bins(bx,nbins_xyz) .and. &
            within_bins(by,nbins_xyz) .and. &
            within_bins(bz,nbins_xyz) &
          ) then
              gof_rot_xyz(bx,by,bz,br) = gof_rot_xyz(bx,by,bz,br) + 1
              gof_rot_xyz_avg(bx,by,bz) = gof_rot_xyz_avg(bx,by,bz) + r
        end if
      end if
    end do
  end do

end subroutine vpi_eval_gof_rot
!@-node:gcross.20090623152316.31:vpi_eval_gof_rot
!@+node:gcross.20090706131953.1754:vpi_eval_nrdm
subroutine vpi_eval_nrdm(nrdm, x1, x2)
  integer :: id, nbins, ndim
  real(kind=b8), dimension ( 2**N_OD_PARTICLE, 2**N_OD_PARTICLE ) :: nrdm
  real(kind=b8), dimension ( N_PARTICLE ) :: x1, x2
  real :: xsize, dndx

  integer :: i,j,bx1, bx2

  bx1 = 1
  bx2 = 1
  do i = 0, N_OD_PARTICLE-1
    if(x1(i+1) > 0) then
      bx1 = bx1 + 2**i
    end if
    if(x2(i+1) > 0) then
      bx2 = bx2 + 2**i
    end if
  end do

  nrdm(bx1,bx2) = nrdm(bx1,bx2) + 1

end subroutine vpi_eval_nrdm
!@-node:gcross.20090706131953.1754:vpi_eval_nrdm
!@+node:gcross.20090623152316.16:vpi_eval_dphase
subroutine vpi_eval_dphase( rho, x, nbins, size_x, dndx )
  real(kind=b8), dimension( : , : ), intent(out):: rho
  real(kind=b8), dimension( : , : ) :: x
  integer nbins
  real :: size_x, dndx
  integer bin

  integer :: i,j

  do i = 1, size(x,1)
    do j = 1, size(x,2)
      bin = floor((x(i,j)+size_x)*dndx)+1
      if ( within_bins(bin,nbins) ) then
        rho(bin,j) = rho(bin,j) + 1
      end if 
    end do
  end do
end subroutine vpi_eval_dphase
!@-node:gcross.20090623152316.16:vpi_eval_dphase
!@+node:gcross.20090623152316.21:vpi_eval_N_R
subroutine vpi_eval_N_R( N_R, dn_rms, x, islice, dslice )
  real(kind=b8), dimension( : ), intent(inout):: N_R 
  real(kind=b8), intent(inout):: dn_rms
  integer, intent(in):: islice, dslice
  real(kind=b8), dimension( : , : , : ), intent(in) :: x

  integer :: i,j, cnt

  stop "Evaluation of N_R has been disabled."

!@+at
!   do j = islice-dslice, islice+dslice
!     cnt = 1
!     do i = 1, N_PARTICLE
!       if ( x(j,i,3) .gt. gw_z0 ) then
!         cnt = cnt + 1
!       end if
!     end do
!     N_R(cnt) = N_R(cnt) + 1.0_b8/(2.0_b8*dslice+1.0_b8)
!   end do
!   dn_rms = dn_rms + (cnt-N_Particle/2)**2
!@-at
!@@c

end subroutine vpi_eval_N_R
!@-node:gcross.20090623152316.21:vpi_eval_N_R
!@+node:gcross.20090623152316.22:vpi_eval_B_z2
subroutine vpi_eval_N_z2( N_zs,N_za, x )
  real(kind=b8), dimension( : ), intent(inout):: N_zs, N_za
  real(kind=b8), dimension( : , : ) :: x

  integer :: i, tr_cnt, tl_cnt, br_cnt, bl_cnt

  tr_cnt = 0
  tl_cnt = 0
  br_cnt = 0
  bl_cnt = 0

  do i = N_PARTICLE-N_PARTICLE2+1, N_PARTICLE
    if ( x(i,3) .gt. 0 ) then
      if ( x(i,1) .gt. 0 ) then
        tr_cnt = tr_cnt + 1
      else 
        tl_cnt = tl_cnt + 1
      end if
    else 
      if ( x(i,1) .gt. 0 ) then
        br_cnt = br_cnt + 1
      else 
        bl_cnt = bl_cnt + 1
      end if
    end if 
  end do

  N_zs(tr_cnt+tl_cnt+1) = N_zs(tr_cnt+tl_cnt+1) + 1
  N_za((tr_cnt-tl_cnt)+N_PARTICLE2+1) = N_za((tr_cnt-tl_cnt)+N_PARTICLE2+1) + 1
end subroutine vpi_eval_N_z2
!@-node:gcross.20090623152316.22:vpi_eval_B_z2
!@-node:gcross.20090626170642.1722:No idea what these do yet
!@-node:gcross.20090624144408.2048:Evaluators of physical quantities
!@+node:gcross.20090624144408.2049:Input
!@+node:gcross.20090624144408.2051:read_configuration_file
subroutine read_configuration_file(pid)
  integer :: pid

  logical :: iexist
  character (len=30) :: param_file
  character (len=50) :: pstr

  write(param_file,"(a8)") "param.in"

  inquire(file=param_file, exist=iexist)
  if (iexist) then
     ! Seed file exists, read seed from seedfile
     open(10, file=param_file, status="old")

     call init_global_parameters

     call init_sp_potential
     call init_tb_potential     
     call init_sp_tfunc
     call init_jas_tfunc

     close(10)
  endif

end subroutine read_configuration_file
!@-node:gcross.20090624144408.2051:read_configuration_file
!@+node:gcross.20090624144408.2052:read_lattice_file
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
!@-node:gcross.20090624144408.2052:read_lattice_file
!@+node:gcross.20090624144408.2053:read_expot_file
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
!@-node:gcross.20090624144408.2053:read_expot_file
!@-node:gcross.20090624144408.2049:Input
!@+node:gcross.20090624144408.2050:Misc
!@+node:gcross.20090706131953.1749:within_bins
function within_bins(index,nbins)
  integer :: index, nbins
  logical :: within_bins

  within_bins = (index >= 1) .and. (index <= nbins)
end function within_bins
!@-node:gcross.20090706131953.1749:within_bins
!@+node:gcross.20090623152316.32:vpi_accept_path
function vpi_accept_path( lngfn0, lngfn1, dphase ) result( accept )
  real(kind = b8) :: lngfn0, lngfn1
  real(kind = b8) :: dphase
  logical :: accept

  real :: Pa, Ptest 

  Pa = dphase*exp(lngfn1 - lngfn0) 
  call random_number( Ptest )
!  print "(2a20)", "Pa", "Ptest"
!  print *, Pa, Ptest
  if ( Pa .gt. Ptest ) then
    accept = .true.
  else
    accept = .false.
  end if

end function vpi_accept_path
!@-node:gcross.20090623152316.32:vpi_accept_path
!@+node:gcross.20090623152316.115:eval_E_local
subroutine eval_E_local(np,ndim,grad_lnspf, lap_lnspf, grad_lnjas, lap_lnjas, U, E_l)
  integer, intent(in) :: np, ndim
  real(kind=b8), dimension( np , ndim ), intent(in) :: grad_lnjas 
  real(kind=b8), dimension( np , ndim ), intent(in) :: grad_lnspf 
  real(kind=b8), intent(in) :: lap_lnspf
  real(kind=b8), intent(in) :: lap_lnjas 
  real(kind=b8), intent(in) :: U
  real(kind=b8), intent(out) :: E_l

  real(kind=b8) :: t0,t1,T

  integer i,j

  t0 = 0.0_b8
!  do i = 1, np
!    do j = 1, ndim
!      t1  = grad_lnspf(i,j) +  grad_lnjas(i,j) 
!      t0 = t0 + t1*t1
!    end do
!  end do

  t0 = sum((grad_lnspf(:,:) +  grad_lnjas(:,:))**2)

  T = t0 + lap_lnspf + lap_lnjas
!  do i = 1, np
!    do j = 1, ndim
!      write(111,"(1a,2g18.6)") "#",grad_lnspf(i,j), grad_lnjas(i,j) 
!    end do
!  end do

!  write (111,"(7a18)") "#grad^2","lap_lnspf","lap_lnjas","T","U","E_l","lambda"
!  write (111,"(7g18.6)") t0,lap_lnspf,lap_lnjas,T,U,E_l,lambda
  E_l = U - lambda*T


end subroutine  eval_E_local
!@-node:gcross.20090623152316.115:eval_E_local
!@+node:gcross.20090721121051.1756:compute_gradU2
subroutine compute_potential (&
! INPUT: particle position / rotation information
    x, xij2, x_rot, &
! INPUT: path slice to consider
    move_start, move_end, &
! INPUT: 
    nslice, np, ndim, &
! OUTPUT: potential and square gradients
    U, gradU2, &
! OUTPUT: whether the move should be rejected outright
    reject_flag &
    )
  implicit none

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
  real(kind=b8), dimension ( nslice, np , np ), intent(in) :: xij2
  real(kind=b8), dimension ( nslice, np , N_DIM_ROT ), intent(in) :: x_rot

!@+at
! Function output
!@-at
!@@c
  real(kind=b8), dimension( nslice, np ), intent(out) :: U
  real(kind=b8), dimension( nslice ), intent(out) :: gradU2
  logical, intent(out) :: reject_flag

!@+at
! Temporary variables.
!@-at
!@@c
  real (kind=b8), dimension( np, ndim ) :: grad_Usp, grad_Uij, grad_U
  real (kind=b8) :: U_sp, U_ij, U_ij_rot

  reject_flag = .false.
  do ii = move_start, move_end
    U(ii,:) = 0.0_b8
    do jj = 1,n_particle
      U_sp = Usp_func( x, ii, jj, N_SLICE, N_PARTICLE, N_DIM )
      U_ij = Uij_func( x, xij2, ii, jj, N_SLICE, N_PARTICLE, N_DIM, reject_flag )
      if(reject_flag) then
          return
      end if
      U_ij_rot=0
      if( eval_rotation ) then
        U_ij_rot = Uij_rot_func( x, x_rot, xij2, ii, jj, N_SLICE, N_PARTICLE, N_DIM )
      end if
      U(ii,jj) = U(ii,jj) + ( U_sp + U_ij + U_ij_rot )
    end do

    if ( gU2_weight(ii) .ne. 0 ) then
      grad_Usp = gUsp_func( x, ii, N_SLICE, N_PARTICLE, N_DIM )
      grad_Uij = gUij_func( x, xij2, ii, N_SLICE, N_PARTICLE, N_DIM )
      grad_U = grad_Usp + grad_Uij
      gradU2(ii) = sum( grad_U(:,1)**2 + grad_U(:,2)**2 + grad_U(:,3)**2 )
    else
      gradU2(ii) = 0
    end if
  end do

  !        print *, "U_1:", U_1
  !        print *, "U_0:", U_0
  !        print *, "gradU2_1:", gradU2_1
  !        print *, "gradU2_0:", gradU2_0

end subroutine
!@-node:gcross.20090721121051.1756:compute_gradU2
!@+node:gcross.20090721121051.1746:compute_log_acceptance_weight
function compute_log_acceptance_weight (&
! INPUT: particle position / rotation information
    x, xij2, x_rot, &
! INPUT: path slice to consider
    move_start, move_end, &
! INPUT: 
    nslice, np, ndim, &
! OUTPUT: potential and square gradients
    U, gradU2, &
! OUTPUT: whether the move should be rejected outright
    reject_flag &
! OUTPUT: computed weight
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

  real(kind=b8), dimension ( nslice, np , ndim ), intent(in) :: x
  real(kind=b8), dimension ( nslice, np , np ), intent(in) :: xij2
  real(kind=b8), dimension ( nslice, np , N_DIM_ROT ), intent(in) :: x_rot

!@+at
! Function output
!@-at
!@@c
  real(kind=b8), dimension( nslice, np ), intent(out) :: U
  real(kind=b8), dimension( nslice ), intent(out) :: gradU2
  real(kind=b8) :: weight
  logical, intent(out) :: reject_flag

!@+at
! Temporary variables.
!@-at
!@@c
  real (kind=b8) :: lngfn, rotgfn, hsgfn, lntfn
  integer :: sl_start, sl_end

  !@  << Compute contribution from potentials >>
  !@+node:gcross.20090626112946.1696:<< Compute contribution from potentials >>
  call compute_potential (&
      x, xij2, x_rot, &
      move_start, move_end, &
      nslice, np, ndim, &
      U, gradU2, &
      reject_flag &
      )

  if (reject_flag) then
    return
  endif
  !@-node:gcross.20090626112946.1696:<< Compute contribution from potentials >>
  !@nl

  !@  << Compute contribution from propagators >>
  !@+node:gcross.20090721121051.1752:<< Compute contribution from propagators >>
  if(move_start .le. 1) then
    sl_start = 1
  else
    sl_start = move_start-1
  end if
  if(move_end .ge. N_SLICE) then
    sl_end = N_SLICE
  else
    sl_end = move_end+1
  end if

  if( n_slice > 2 ) then 
    if ( use_gfn4 ) then
      lngfn = vpi_gfn4_sp( sl_start, sl_end, part_num, U, gradU2, U_weight, gU2_weight, &
                           N_SLICE, N_PARTICLE, N_DIM, lambda, dtau ) 
    else 
      lngfn = vpi_gfn2_sp( sl_start, sl_end, part_num, U, N_SLICE, N_PARTICLE, N_DIM, dtau ) 
    end if

  !@+at
  !   if ( use_HS_gfn ) then
  ! #ifdef USE_IMAGE_HSGFN
  !     hsgfn = vpi_hs_gfn( sl_start, sl_end, part_num, xij2, N_SLICE, 
  ! N_PARTICLE, N_DIM, dtau*lambda )
  ! #else
  !     hsgfn = vpi_hs_gfn2( sl_start, sl_end, part_num, x, N_SLICE, 
  ! N_PARTICLE, N_DIM, dtau*lambda )
  ! #endif
  !     if( hsgfn <= 0 .and. hsgfn1 > 0 ) then
  !       print *, "****ERROR***** hsgfn < 0 ", hsgfn0, hsgfn1
  !     else
  !       lngfn = lngfn + log(hsgfn)
  !     else
  !   end if
  !@-at
  !@@c

    if ( use_HW_gfn ) then
      hsgfn = vpi_hw_gfn( sl_start, sl_end, part_num, x, N_SLICE, N_PARTICLE, N_DIM, dtau )
      lngfn = lngfn + log(hsgfn)
    end if

    if ( eval_rotation ) then
      rotgfn =  gfn_rot(x_rot,part_num,sl_start,sl_end,dtau)
  !          write (1000+my_rank,*) rotgfn
      lngfn = lngfn + log(rotgfn)
    end if
  end if
  !@-node:gcross.20090721121051.1752:<< Compute contribution from propagators >>
  !@nl

  !@  << Compute contribution from trial functions >>
  !@+node:gcross.20090721121051.1754:<< Compute contribution from trial functions >>
  lntfn = tfunc(x, 1, N_SLICE, N_PARTICLE, N_DIM) + &
          tfunc(x, N_SLICE, N_SLICE, N_PARTICLE, N_DIM) + &
          jas_tfun(x, xij2, 1, N_SLICE, N_PARTICLE, N_DIM) + &
          jas_tfun(x, xij2, N_SLICE, N_SLICE, N_PARTICLE, N_DIM)
  !@-node:gcross.20090721121051.1754:<< Compute contribution from trial functions >>
  !@nl

  weight = lngfn + lntfn

  return

end function
!@-node:gcross.20090721121051.1746:compute_log_acceptance_weight
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
!@-node:gcross.20090624144408.2050:Misc
!@-others
!@-node:gcross.20090623152316.13:Subroutines
!@-others

end program Test_VPI

!@-node:gcross.20090623152316.3:@thin vpi_main.f90
!@-leo
