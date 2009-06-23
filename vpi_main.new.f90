program Test_VPI

  
  use vpi_rand_utils
  use vpi_coordinate_utils
  use vpi_gfn
  use vpi_potential, Usp_func => vpi_Usp_dwell, &
                     Uij_func => vpi_Uij_hs, &
                     gUsp_func => vpi_gUsp_dwell, &
                     gUij_func => vpi_gUij_hs, &
                     Uij_rot_func => vpi_Uij_dimer
  use vpi_trial_func, tfunc => dw_tfun, &
                      jas_tfun => hs_tfun, &
                      grad_lap_sp_tfun => numeric_grad_lap_spf, &
                      grad_lap_jas_tfun => grad_lap_hs_tfun
  use vpi_phase_maps, update_phase => vortex
  use vpi_bisect
  use vpi_bbridge
  use vpi_sample
  use vpi_rotation
  use vpi_defines

  use vpi_lattice
  use vpi_xij
  use vpi_obdm
  use vpi_read_input
  use timers

  implicit none

#ifdef USE_MPI
  include 'mpif.h'
#endif
  INTEGER :: short, medium, long, vlong
  PARAMETER (short = SELECTED_INT_KIND(2), &
             medium= SELECTED_INT_KIND(4), &
             long  = SELECTED_INT_KIND(10),&
             vlong = SELECTED_INT_KIND(100))

  real (kind=b8), allocatable, dimension( :, :, : ) :: q0, q1
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
  real (kind=b8), allocatable, dimension( :, : ) :: grad_Usp, grad_Uij, grad_U
  real (kind=b8), allocatable, dimension( : ) :: U_weight, gU2_weight

  real (kind=b8), allocatable, dimension( : ) :: N_R
  real (kind=b8), allocatable, dimension( : ) :: N_zs
  real (kind=b8), allocatable, dimension( : ) :: N_za
  real (kind=b8), allocatable, dimension( :, : ) :: corr_N_R
  real (kind=b8), allocatable, dimension( :, : ) :: nrdm_z


  real (kind=b8), allocatable, dimension( : , : ) :: rho
  real (kind=b8), allocatable, dimension( : ) :: rho_r
  real (kind=b8), allocatable, dimension( : ) :: gofr
  real (kind=b8), allocatable, dimension( : ) :: gofr_rot
  real (kind=b8), allocatable, dimension( : ) :: gof_rot
  real (kind=b8), allocatable, dimension( : , : , : , : ) :: gof_rot_xyz
  real (kind=b8), allocatable, dimension( : , : , : ) :: gof_rot_xyz_avg
  real (kind=b8), allocatable, dimension( : , : ) :: trial_rho_L
  real (kind=b8), allocatable, dimension( : , : ) :: trial_rho_R
  real (kind=b8), allocatable, dimension( : , : , : ) :: full_density
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
!  real (kind=b8), allocatable, dimension( : , : , : , : , : , : ) :: obdm_xyz
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_rot_z
  real (kind=b8), allocatable, dimension( : , : ) :: obdm_theta

  real (kind=b8) :: U_sp, U_ij,Uij_rot

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
  real (kind=b8) :: E_l
  real (kind=b8) :: dE = 0

  integer :: i,j,k,l,m,ii,jj,kk,ll,pass
  integer :: coll_iter

  integer(long) :: n_moves, sp_moves, swap_moves, bi_moves, end_moves, split_moves
  integer(long) :: left_moves, right_moves
  integer :: mtype
  integer :: sp_try, swap_try, bi_try, end_try, split_try 
  integer :: part_num, move_start, move_end
  integer :: tmp_move_start, tmp_move_end
  integer :: sl_start, sl_end
  integer :: islice

  real(kind=b8) :: lngfn0, lngfn1 
  real(kind=b8) :: rotgfn0, rotgfn1 
  real(kind=b8) :: hsgfn0, hsgfn1 
  real(kind=b8) :: lntfn0, lntfn1 

  logical :: acc_flag = .true.
  integer :: ierr
  integer :: my_rank = 0
  logical :: ex
  double precision time0, time1
  character (len=50) :: my_fname


#ifdef USE_MPI
  integer :: p
  integer :: source
  integer :: dest
  integer :: tag=0
  character (len = 100) :: message
  integer, dimension(MPI_STATUS_SIZE) :: status
  integer :: instr_size, outstr_size

  call MPI_Init(ierr)

  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
!   write(6,*)'hi, this is my rank',my_rank
  call MPI_Comm_size(MPI_COMM_WORLD, p, ierr)
#endif

  call read_input_file(my_rank)

  CSLICE = N_SLICE/2
  a_hs2 = a_hs*a_hs
  a_dw2 = a_dw**2
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


  allocate( q0(N_SLICE, N_PARTICLE, N_DIM) )
  allocate( q1(N_SLICE, N_PARTICLE, N_DIM) )
  allocate( q_rot0(N_SLICE, N_PARTICLE, N_DIM_ROT) )
  allocate( q_rot1(N_SLICE, N_PARTICLE, N_DIM_ROT) )
  allocate( qsq_sum( N_SLICE, N_DIM ) )
  allocate( mcount( N_SLICE ) )
  allocate( accr_v( N_SLICE ) )
  allocate( sp_qtest( N_SLICE, N_DIM ) )
  allocate( tmp_sp_chain( N_SLICE ) )
  allocate( ss_qtest( N_PARTICLE, N_DIM ) )
  allocate( xij2_0( N_SLICE, N_PARTICLE, N_PARTICLE ) )
  allocate( xij2_1( N_SLICE, N_PARTICLE, N_PARTICLE ) )
  allocate( U_0( N_SLICE, N_PARTICLE ) )
  allocate( U_1( N_SLICE, N_PARTICLE ) )
  allocate( U_cum( N_SLICE ) )
  allocate( U_avg( N_SLICE ) )
  allocate( gradU2_0( N_SLICE ) )
  allocate( gradU2_1( N_SLICE ) )
  allocate( grad_lntfn( N_PARTICLE, N_DIM ) )
  allocate( grad_lnjas( N_PARTICLE, N_DIM ) )
  allocate( phase_0( N_SLICE ) )
  allocate( phase_sum_0( N_SLICE ) )
  allocate( phase_1( N_SLICE ) )
  allocate( corr_x( N_SLICE, N_DIM ) )
  allocate( corr_xz( N_SLICE, 2 ) )
  allocate( corr_r( N_SLICE ) )
  allocate( corr_phase( N_SLICE ) )
  allocate( N_R( N_PARTICLE + 1 ) )
  allocate( N_zs( N_PARTICLE2 + 1 ) )
  allocate( N_za( N_PARTICLE2 + 1 ) )
  allocate( corr_N_R( N_SLICE, N_PARTICLE + 1 ) )
  allocate( nrdm_z( N_OD_PARTICLE**2, N_OD_PARTICLE**2 ) )
  allocate( grad_Usp( N_PARTICLE , N_DIM ) )
  allocate( grad_Uij( N_PARTICLE , N_DIM ) )
  allocate( grad_U( N_PARTICLE , N_DIM ) )
  allocate( U_weight( N_SLICE ) )
  allocate( gU2_weight( N_SLICE ) )

  allocate( rho( N_BINS, N_DIM ) )
  allocate( rho_r( N_BINS ) )
  allocate( gofr( N_BINS_GOFR ) )
  allocate( gofr_rot( N_BINS_GOFR ) )
  allocate( gof_rot( N_BINS_ROT ) )
  allocate( gof_rot_xyz( N_BINS_ROT_XYZ, N_BINS_ROT_XYZ, N_BINS_ROT_XYZ, N_BINS_ROT ) )
  allocate( gof_rot_xyz_avg(N_BINS_ROT_XYZ, N_BINS_ROT_XYZ, N_BINS_ROT_XYZ ) )
  allocate( trial_rho_L( N_BINS, N_DIM ) )
  allocate( trial_rho_R( N_BINS, N_DIM ) )
  allocate( full_density( N_BINS, N_BINS, N_BINS ) )
  allocate( odxy( N_BINS, N_BINS ) )
  allocate( odxz( N_BINS, N_BINS ) )
  allocate( odyz( N_BINS, N_BINS ) )
  allocate( odxy_avg( N_BINS, N_BINS ) )
  allocate( odxz_avg( N_BINS, N_BINS ) )
  allocate( odyz_avg( N_BINS, N_BINS ) )
  allocate( odxy_sig( N_BINS, N_BINS ) )
  allocate( odxz_sig( N_BINS, N_BINS ) )
  allocate( odyz_sig( N_BINS, N_BINS ) )
  allocate( od2( N_BINS, N_BINS ) )
  allocate( od2_avg( N_BINS, N_BINS ) )
  allocate( od2_sig( N_BINS, N_BINS ) )
  allocate( dxlxr( N_BINS, N_DIM ) )
  allocate( rdm2_z(N_BINS, N_BINS ) )
  allocate( rdm2_theta( N_BINS, N_BINS ) )
  allocate( rdm22_theta( N_BINS, N_BINS ) )
  allocate( obdm_x(N_BINS, N_BINS ) )
  allocate( obdm_y(N_BINS, N_BINS ) )
  allocate( obdm_z(N_BINS, N_BINS ) )
!  allocate( obdm_xyz(N_BINS, N_BINS, N_BINS, N_BINS, N_BINS, N_BINS ) )
  allocate( obdm_rot_z(N_BINS, N_BINS ) )
  allocate( obdm_theta(N_BINS, N_BINS ) )

  mcount = 0
  accr_v = 0
  U_cum = 0
  corr_x = 0.0
  corr_xz = 0.0
  corr_phase = 0.0
  rho = 0.0
  rho_r = 0.0
  gofr = 0.0
  gofr_rot = 0.0
  gof_rot = 0.0
  trial_rho_L = 0.0
  trial_rho_R = 0.0
  N_R = 0.0
  N_zs = 0.0
  N_za = 0.0
  corr_N_R = 0.0
  U_0 = 0.0
  U_1 = 0.0



  if ( my_rank .eq. 0 ) then
    open(12, file="params.dat", status="replace")
    write(12,*)  "N_SLICE",N_SLICE
    write(12,*)  "CSLICE",CSLICE
    write(12,*)  "N_PARTICLE",N_PARTICLE
    write(12,*)  "N_DIM",N_DIM
    write(12,*)  "N_BINS", N_BINS
    write(12,*)  "dM", dM
    write(12,*)  "nblocks", nblocks
    write(12,*)  "nsteps", nsteps
    write(12,*)  "dnu", dnu
    write(12,*)  "dswap", dswap
    write(12,*)  "drot", drot
    write(12,*)  "dtau", dtau
    write(12,*)  "size_x", size_x
    write(12,*)  "eval_off_diagonal", eval_off_diagonal
    write(12,*)  "od_pnum", od_pnum
    write(12,*)  "use_gfn4",use_gfn4
    write(12,*)  "lambda",lambda
    write(12,*)  "a_dw", a_dw
    write(12,*)  "ep_dw", ep_dw
    write(12,*)  "e_lj", e_lj
    write(12,*)  "a_lj", a_lj
    write(12,*)  "lam_ho", lam_ho
    write(12,*)  "p_hox", p_hox
    write(12,*)  "p_hoy", p_hoy
    write(12,*)  "p_hoz", p_hoz
    write(12,*)  "p_ljc5", p_ljc5
    write(12,*)  "p_ljc1", p_ljc1
    write(12,*)  "e0_nua", e0_nua
    write(12,*)  "m_nua", m_nua
    write(12,*)  "lam_nw", m_nua
    write(12,*)  "gw_e", gw_e
    write(12,*)  "gw_asq", gw_asq
    write(12,*)  "gw_z0", gw_z0
    write(12,*)  "run_type:", run_type
    write(12,*)  "ap_e::", ap_e 
    write(12,*)  "ap_a::", ap_a
    write(12,*)  "ap_e_norm::", ap_e_norm
    write(12,*)  "ap_lam::", ap_lam
    write(12,*)  "ae3::", ae3
    write(12,*)  "ab_e::", ab_e
    write(12,*)  "ab_a::", ab_a
    write(12,*)  "ab_e_norm::", ab_e_norm
    write(12,*)  "ab_x0::", ab_x0
    write(12,*)  "domega::", domega
    write(12,*)  "p_MO_a::", p_MO_a
    write(12,*)  "p_MO_ap::", p_MO_ap
    write(12,*)  "p_MO_r0::", p_MO_r0
    write(12,*)  "p_MO_De::", p_MO_De
    write(12,*)  "p_sc_a::", p_sc_a
    write(12,*)  "p_sc_w::", p_sc_w
    write(12,*)  "a_hs::", a_hs

    write(12,*)  "p_lattice_vb::", p_lattice_vb
    write(12,*)  "p_lattice_ax::", p_lattice_ax
    write(12,*)  "p_lattice_az::", p_lattice_az
    write(12,*)  "p_lwx::", p_lwx
    write(12,*)  "p_lwz::", p_lwz

    write(12,*)  "eval_rotation", eval_rotation
    write(12,*)  "p_B", p_B
    write(12,*)  "e_dimer", e_dimer

    write(12,*)  "PROB_BBRIDGE_MOVE", PROB_BBRIDGE_MOVE
    write(12,*)  "PROB_RIGID_MOVE", PROB_RIGID_MOVE
    write(12,*)  "PROB_SWAP_MOVE", PROB_SWAP_MOVE
    write(12,*)  "PROB_ROT_MOVE", PROB_ROT_MOVE
    write(12,*)  "force_space_flip",  force_space_flip
    do i = 1, N_VORTEX
      write(12,*)  p_vcoords(i,1), p_vcoords(i,2)
    end do

    close(12)
  end if


  accr = 0
  sp_accr = 0
  swap_accr = 0
  bi_accr = 0
  end_accr = 0
  split_accr = 0

  n_moves = 0
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

  ! these "slice weights" are used to give the proper weight to even
  ! and odd terms when calculating the Green's function.  
  ! in the fourth order short time approximation for the Green's funciton,
  ! odd slices contribute twice the potential energy of even slices
  ! the gradient of the potential energy only contributes for odd slides.
  U_weight = 0.0
  gU2_weight = 0.0
  if (use_gfn4) then
    do ii = 1, CSLICE
      U_weight(ii) = mod(ii+1,2) + 1 
      gU2_weight(ii) = mod(ii+1,2)
    end do
    do ii = CSLICE+2, N_SLICE
      U_weight(ii) = mod(ii,2) + 1
      gU2_weight(ii) = mod(ii,2)
    end do
  else
    U_weight = 1.0
    gU2_weight = 0.0
  end if

  if ( eval_off_diagonal ) then
    U_weight(CSLICE+1) = U_weight(CSLICE)
    gU2_weight(CSLICE+1) = gU2_weight(CSLICE)
  else
    U_weight(CSLICE+1) = 0
    gU2_weight(CSLICE+1) = 0
  end if


  ! read in seed file / create a new seed for random number generator
    call ru_init_seed_file(my_rank)
    if (eval_rotation) then
      gfn_v = init_rot_gfn (size(gfn_v), L_MAX, dtau)
      call write_rot_gfn(gfn_v)
    end if

    q0 = 0.0
    q1 = 0.0
    !print *, q0(1,:,:)
    !print *, "pfunc = ", pfunc(q0(1,:,:),N_PARTICLE, N_DIM)
    !print *, "test: ", test_fvar(q0,pfunc)

    print *, "*** INITIALIZE ***"
     
    write(my_fname,"(a14,i4.4)")"last_path.dat.",my_rank
    inquire(file=my_fname, exist=ex)
    if ( ex ) then
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
        q0 = q1
      endif 
    else
      call vpi_make_lattice( q0, size_x(1), GAUSSIAN_FILL) !NOT SURE OF THIS MODIF
    end if

 
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

    if(force_space_flip) then
      q0(CSLICE+2:N_SLICE,:,3) = -q0(CSLICE+2:N_SLICE,:,3)
    end if

    do ii = 1, N_PARTICLE
      call vpi_update_xij( xij2_0, q0, 1, N_SLICE, ii, N_SLICE, N_PARTICLE, N_DIM  )
    end do

    if( use_eval_phase ) then
      call update_phase( phase_0, 1, N_SLICE, q0, xij2_0 )
    end if
    phase_1 = phase_0

    acc_flag = .true.
    do ii = 1, N_SLICE
      do jj = 1, N_PARTICLE
        U_sp = Usp_func( q0, ii, jj, N_SLICE, N_PARTICLE, N_DIM )
        U_ij = Uij_func( q0, xij2_0, ii, jj, N_SLICE, N_PARTICLE, N_DIM, acc_flag )
        if( eval_rotation ) then
          Uij_rot =  Uij_rot_func( q0, q_rot0, xij2_0, ii, jj, N_SLICE, N_PARTICLE, N_DIM )
          U_ij = U_ij + Uij_rot
        end if
        if(acc_flag .eqv. .false.) then
          print *,"ERROR: particle overlap"
          stop
        end if
        U_0(ii,jj) = ( U_sp + U_ij )
      end do
      grad_Usp = gUsp_func( q0, ii, N_SLICE, N_PARTICLE, N_DIM )
      grad_Uij = gUij_func( q0, xij2_0, ii, N_SLICE, N_PARTICLE, N_DIM )
      grad_U = grad_Usp + grad_Uij
      gradU2_0(ii) = sum( grad_U(:,1)**2 + grad_U(:,2)**2 + grad_U(:,3)**2 )
    end do

#ifdef USE_SWFN 
    lntfn0 = 2.0_b8*tfunc(q0, CSLICE, N_SLICE, N_PARTICLE, N_DIM) + 2.0_b8*jas_tfun(xij2_0, CSLICE, N_SLICE, N_PARTICLE)
#else
    lntfn0 = tfunc(q0, 1, N_SLICE, N_PARTICLE, N_DIM) + tfunc(q0, N_SLICE, N_SLICE, N_PARTICLE, N_DIM) + jas_tfun(xij2_0, 1, N_SLICE, N_PARTICLE) + jas_tfun(xij2_0, N_SLICE, N_SLICE, N_PARTICLE)
#endif

    q1 = q0
    U_1(:,:) = U_0(:,:)
    gradU2_1(:) = gradU2_0(:)
    xij2_1(:,:,:) = xij2_0(:,:,:)

    qsq_sum = 0.0

    do i = 1,nblocks
      dE = 0.0_b8
      E_left = 0.0_b8
      E_center = 0.0_b8
      E_right = 0.0_b8
#ifdef USE_MPI
      time0 = MPI_Wtime()
#else
      time0 = day_timer()
#endif
      write(my_fname,"(a10,i4.4)")"xl-xr.dat.",my_rank
      open(13, file=my_fname,  POSITION="APPEND")
      do j = 1,nsteps
#ifdef DEBUG_L0
        do k = 1, N_SLICE
          do l = 1, N_PARTICLE
            do m = 1, N_DIM
              if(q1(k,l,m) .ne. q0(k,l,m)) then
                print *, "q0 != q1",k,l,m
              end if
            end do
          end do
        end do
#endif

        call sample_scheme1(q0, q1, move_start, move_end, part_num, mtype)
        if( eval_rotation ) then
          call sample_rotation(q_rot0, q_rot1, move_start, move_end, part_num)
        end if

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

        call vpi_update_xij( xij2_1, q1, move_start, move_end, part_num, N_SLICE, N_PARTICLE, N_DIM  )
#ifdef EXPERIMENTAL
        acc_flag = .true.
        do ii = move_start, move_end
          acc_flag = .true.
          U_ij = Uij_func( q1, xij2_1, ii, part_num, N_SLICE, N_PARTICLE, N_DIM, acc_flag )
          coll_iter = 0
          do while ((coll_iter < 20) .and. (acc_flag .eqv. .false.))
            tmp_move_start = ii-1
            tmp_move_end = ii+1
            if(tmp_move_start < 1) then
              tmp_move_start = 1
              call bbridge(q1, q1, part_num, 1, N_DIM, LEFT_BRIDGE, tmp_move_start,tmp_move_end, lambda, dtau)
            else 
              if(tmp_move_end > N_SLICE) then
                tmp_move_end = N_SLICE
                call bbridge(q1, q1, part_num, 1, N_DIM, RIGHT_BRIDGE, tmp_move_start,tmp_move_end, lambda, dtau)
              else
                call bbridge(q1, q1, part_num, 1, N_DIM, S_BRIDGE, tmp_move_start,tmp_move_end, lambda, dtau)
              end if
            end if
            call vpi_update_xij( xij2_1, q1, ii, ii, part_num, N_SLICE, N_PARTICLE, N_DIM  )
            acc_flag = .true.
            U_ij = Uij_func( q1, xij2_1, ii, part_num, N_SLICE, N_PARTICLE, N_DIM, acc_flag )
            coll_iter = coll_iter + 1
          end do
          if(acc_flag .eqv. .false.) then 
            print *, coll_iter
          end if
        end do 
#endif

        if ( use_eval_phase ) then
          call update_phase( phase_1, move_start, move_end, q1, xij2_1 )
        end if 


        acc_flag = .true.
        do ii = move_start, move_end
          U_sp = Usp_func( q1, ii, part_num, N_SLICE, N_PARTICLE, N_DIM )
          U_ij = Uij_func( q1, xij2_1, ii, part_num, N_SLICE, N_PARTICLE, N_DIM, acc_flag )
          if(acc_flag .eqv. .false.) then
!            print *,"overlap", mtype
            exit 
          end if
          if( eval_rotation ) then
            Uij_rot =  Uij_rot_func( q1, q_rot1, xij2_1, ii, part_num, N_SLICE, N_PARTICLE, N_DIM )
            U_ij = U_ij + Uij_rot
          end if
          U_1(ii,part_num) = ( U_sp + U_ij )

          U_sp = Usp_func( q0, ii, part_num, N_SLICE, N_PARTICLE, N_DIM )
          U_ij = Uij_func( q0, xij2_0, ii, part_num, N_SLICE, N_PARTICLE, N_DIM, acc_flag )
          if(acc_flag .eqv. .false.) then
            print *,"ERROR, overlap in supposedly clean path"
            stop
          end if
          if( eval_rotation ) then
            Uij_rot =  Uij_rot_func( q0, q_rot0, xij2_0, ii, part_num, N_SLICE, N_PARTICLE, N_DIM )
            U_ij = U_ij + Uij_rot
          end if
          U_0(ii,part_num) = ( U_sp + U_ij )

          if ( gU2_weight(ii) .ne. 0 ) then
            grad_Usp = gUsp_func( q1, ii, N_SLICE, N_PARTICLE, N_DIM )
            grad_Uij = gUij_func( q1, xij2_1, ii, N_SLICE, N_PARTICLE, N_DIM )
            grad_U = grad_Usp + grad_Uij
            gradU2_1(ii) = sum( grad_U(:,1)**2 + grad_U(:,2)**2 + grad_U(:,3)**2 )

            grad_Usp = gUsp_func( q0, ii, N_SLICE, N_PARTICLE, N_DIM )
            grad_Uij = gUij_func( q0, xij2_0, ii, N_SLICE, N_PARTICLE, N_DIM )
            grad_U = grad_Usp + grad_Uij
            gradU2_0(ii) = sum( grad_U(:,1)**2 + grad_U(:,2)**2 + grad_U(:,3)**2 )
          else
            gradU2_1(ii) = 0
            gradU2_0(ii) = 0
          end if
        end do

        if(acc_flag .eqv. .false.) then
          ccnt = ccnt + 1
          goto 10
        end if

!        print *, "U_1:", U_1
!        print *, "U_0:", U_0
!        print *, "gradU2_1:", gradU2_1
!        print *, "gradU2_0:", gradU2_0

        
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
       

        if ( use_gfn4 ) then
          lngfn1 = vpi_gfn4_sp( sl_start, sl_end, part_num, U_1, gradU2_1, U_weight, gU2_weight, N_SLICE, N_PARTICLE, N_DIM, lambda, dtau ) 
          lngfn0 = vpi_gfn4_sp( sl_start, sl_end, part_num, U_0, gradU2_0, U_weight, gU2_weight, N_SLICE, N_PARTICLE, N_DIM, lambda, dtau ) 
        else 
          lngfn1 = vpi_gfn2_sp( sl_start, sl_end, part_num, U_1, N_SLICE, N_PARTICLE, N_DIM, dtau ) 
          lngfn0 = vpi_gfn2_sp( sl_start, sl_end, part_num, U_0, N_SLICE, N_PARTICLE, N_DIM, dtau ) 
        end if

        if ( use_HS_gfn ) then
#ifdef USE_IMAGE_HSGFN
          hsgfn1 = vpi_hs_gfn( sl_start, sl_end, part_num, xij2_1, N_SLICE, N_PARTICLE, N_DIM, dtau*lambda )
          hsgfn0 = vpi_hs_gfn( sl_start, sl_end, part_num, xij2_0, N_SLICE, N_PARTICLE, N_DIM, dtau*lambda )
#else
          hsgfn1 = vpi_hs_gfn2( sl_start, sl_end, part_num, q1, N_SLICE, N_PARTICLE, N_DIM, dtau*lambda )
          hsgfn0 = vpi_hs_gfn2( sl_start, sl_end, part_num, q0, N_SLICE, N_PARTICLE, N_DIM, dtau*lambda )
#endif
          if( hsgfn0 > 0 .and. hsgfn1 > 0 ) then
            lngfn1 = lngfn1 + log(hsgfn1)
            lngfn0 = lngfn0 + log(hsgfn0)
          else
            print *, "****ERROR***** hsgfn < 0 ", hsgfn0, hsgfn1
            lngfn1 = 0
          endif
        end if

        if ( use_HW_gfn ) then
          hsgfn1 = vpi_hw_gfn( sl_start, sl_end, part_num, q1, N_SLICE, N_PARTICLE, N_DIM, dtau*lambda )
          hsgfn0 = vpi_hw_gfn( sl_start, sl_end, part_num, q0, N_SLICE, N_PARTICLE, N_DIM, dtau*lambda )
          lngfn1 = lngfn1 + log(hsgfn1)
          lngfn0 = lngfn0 + log(hsgfn0)
        end if

        if ( eval_rotation ) then
          rotgfn1 =  gfn_rot(q_rot1,part_num,sl_start,sl_end,dtau)
          rotgfn0 =  gfn_rot(q_rot0,part_num,sl_start,sl_end,dtau)
!          write (1000+my_rank,*) rotgfn0, rotgfn1
          lngfn1 = lngfn1 + log(rotgfn1)
          lngfn0 = lngfn0 + log(rotgfn0)
        end if


        lntfn1 = tfunc(q1, 1, N_SLICE, N_PARTICLE, N_DIM) + tfunc(q1, N_SLICE, N_SLICE, N_PARTICLE, N_DIM) + jas_tfun(xij2_1, 1, N_SLICE, N_PARTICLE) + jas_tfun(xij2_1, N_SLICE, N_SLICE, N_PARTICLE)
        lntfn0 = tfunc(q0, 1, N_SLICE, N_PARTICLE, N_DIM) + tfunc(q0, N_SLICE, N_SLICE, N_PARTICLE, N_DIM) + jas_tfun(xij2_0, 1, N_SLICE, N_PARTICLE) + jas_tfun(xij2_0, N_SLICE, N_SLICE, N_PARTICLE)


!        print "(2a20)", "lngfn0", "lngfn1"
!        print *, lngfn0, lngfn1
!        print "(2a20)", "lntfn0", "lntfn1"
!        print *, lntfn0, lntfn1
        if (impose_fixed_phase) then
          dphase = product(phase_1(move_start:move_end))/product(phase_0(move_start:move_end))
        else
          dphase = 1.0_b8
        end if
10      if ( vpi_accept_path( lngfn0 + lntfn0, lngfn1 + lntfn1, dphase) .and. acc_flag) then
!          print *, "!!! ACCEPTED !!!"
          q0(move_start:move_end,part_num,:) = q1(move_start:move_end,part_num,:)
          if( eval_rotation ) then
            q_rot0(move_start:move_end,part_num,:) = q_rot1(move_start:move_end,part_num,:)
          end if
          xij2_0(move_start:move_end,:,:) = xij2_1(move_start:move_end,:,:) 
          U_0(move_start:move_end,part_num) = U_1(move_start:move_end,part_num)
          gradU2_0(move_start:move_end) =  gradU2_1(move_start:move_end)
          if( use_eval_phase ) then
            phase_0(move_start:move_end) = phase_1(move_start:move_end)
          end if
          accr = accr + 1;
          accr_v(move_start:move_end) = accr_v(move_start:move_end) + 1
          lntfn0 = lntfn1
          bi_accr = bi_accr + bi_try
          sp_accr = sp_accr + sp_try
          swap_accr = swap_accr + swap_try
          end_accr = end_accr + end_try
          split_accr = split_accr + split_try
        else
!          print *, "!!! REJECTED !!!"
          q1(move_start:move_end,part_num,:) = q0(move_start:move_end,part_num,:)
          if( eval_rotation ) then
            q_rot1(move_start:move_end,part_num,:) = q_rot0(move_start:move_end,part_num,:)
          end if
          xij2_1(move_start:move_end,:,:) = xij2_0(move_start:move_end,:,:)
          U_1(move_start:move_end,part_num) = U_0(move_start:move_end,part_num)
          gradU2_1(move_start:move_end) =  gradU2_0(move_start:move_end)
        end if

        n_moves = n_moves + 1

        pass = grad_lap_sp_tfun( q0, CSLICE, N_PARTICLE, N_DIM, N_SLICE, grad_lntfn, lap_lntfn )
        pass = grad_lap_jas_tfun( q0, xij2_0, CSLICE, N_PARTICLE, N_DIM, N_SLICE, grad_lnjas, lap_lnjas )
        call eval_E_local(N_PARTICLE, N_DIM, grad_lntfn, lap_lntfn, grad_lnjas, lap_lnjas, sum(U_0(CSLICE,:)), E_l)
        E_center = E_center + E_l

        pass = grad_lap_sp_tfun( q0, 1, N_PARTICLE, N_DIM, N_SLICE, grad_lntfn, lap_lntfn )
        pass = grad_lap_jas_tfun( q0, xij2_0, 1, N_PARTICLE, N_DIM, N_SLICE, grad_lnjas, lap_lnjas )
        call eval_E_local(N_PARTICLE, N_DIM, grad_lntfn, lap_lntfn, grad_lnjas, lap_lnjas, sum(U_0(1,:)), E_l)
        E_left = E_left + E_l

        pass = grad_lap_sp_tfun( q0, N_SLICE, N_PARTICLE, N_DIM, N_SLICE, grad_lntfn, lap_lntfn )
        pass = grad_lap_jas_tfun( q0, xij2_0, N_SLICE, N_PARTICLE, N_DIM, N_SLICE, grad_lnjas, lap_lnjas )
        call eval_E_local(N_PARTICLE, N_DIM, grad_lntfn, lap_lntfn, grad_lnjas, lap_lnjas, sum(U_0(N_SLICE,:)), E_l)
        E_right = E_right + E_l

        dE = dE + U_0(CSLICE+1,part_num) - U_0(CSLICE,part_num)
        if(eval_qsq_sum) then
          do k = 1, N_SLICE
            U_cum(k) = U_cum(k) + sum(U_0(k,:))/(N_PARTICLE)
            U_avg(k) = U_cum(k)/n_moves
            phase_sum_0(k) = phase_sum_0(k) + phase_0(k)
            do m = 1, N_DIM
              qsq_sum(k,m) = qsq_sum(k,m) + sum(q0(k,1:N_PARTICLE,m)**2)
            end do
          end do
        end if

        call vpi_eval_density( rho, q0(CSLICE,:,:), N_BINS, size_x/2, dndx )
        call vpi_eval_radial_density( rho_r, q0(CSLICE,:,:), N_BINS, dndr ) 
        if ( eval_off_diagonal ) then
          call vpi_eval_density_sp( dxlxr, (q0(CSLICE,od_pnum,:)-q0(CSLICE+1,od_pnum,:)), N_BINS, size_x/2, dndx )
        end if
        call vpi_eval_density( trial_rho_L, q0(1,:,:), N_BINS, size_x/2, dndx )
        call vpi_eval_density( trial_rho_R, q0(N_SLICE,:,:), N_BINS, size_x/2, dndx )
        if(compute_full_density) then
          call vpi_eval_full_density(full_density, q0(CSLICE,:,:), N_BINS, size_x/2, dndx )
        end if
        if(eval_column_density) then
          odxy_avg = odxy/( n_moves*(N_PARTICLE-N_PARTICLE2)*dxdn(1)*dxdn(2) )
          odxz_avg = odxz/( n_moves*(N_PARTICLE-N_PARTICLE2)*dxdn(1)*dxdn(3) )
          odyz_avg = odyz/( n_moves*(N_PARTICLE-N_PARTICLE2)*dxdn(2)*dxdn(3) )
          if(N_PARTICLE2 > 0 ) then
            od2_avg = od2/( n_moves*N_PARTICLE2*dxdn(1)*dxdn(3) )
          endif
          call vpi_eval_optical_density( odxy, odxy_avg, odxy_sig, q0(CSLICE,:,:), 1, 2, N_BINS, size_x/2, dndx )
          call vpi_eval_optical_density( odxz, odxz_avg, odxz_sig, q0(CSLICE,:,:), 1, 3, N_BINS, size_x/2, dndx )
          call vpi_eval_optical_density( odyz, odyz_avg, odyz_sig, q0(CSLICE,:,:), 2, 3, N_BINS, size_x/2, dndx )
          if( N_PARTICLE2 .gt. 0 ) then 
            call vpi_eval_optical_density2( od2, od2_avg, od2_sig, q0(CSLICE,:,:), N_BINS, size_x/2, dndx )
          end if
        end if
        if (i .lt. nb_discard) then
          odxy_sig = 0 
          odxz_sig = 0 
          odyz_sig = 0 
        end if
        call vpi_eval_gofr( gofr, xij2_0(CSLICE,:,:), N_BINS_GOFR, dndx_GOFR )

        if( eval_rotation ) then
          call vpi_eval_gof_rot( gof_rot, gof_rot_xyz, gof_rot_xyz_avg, gofr_rot, q0, q_rot0, xij2_0, CSLICE, N_BINS_GOFR, dndx_GOFR,  N_BINS_ROT_XYZ, size_x/2, dndx_rot_xyz, N_BINS_ROT, dndx_rot)
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
          call vpi_eval_N_R( N_R, q0(:,:,:), CSLICE, 0 )
        end if
!        call eval_N_z2( N_zs,N_za, q0(CSLICE,:,:) )
!       do ii = 1, N_SLICE
!         call eval_N_R( corr_N_R(ii,:), q0(ii,:,:) )
!       end do
        if(eval_rdm) then
          do ii = 1, N_PARTICLE
            do jj = ii+1, N_PARTICLE
              call vpi_eval_obdm_cut( rdm2_z, q0(CSLICE,ii,3), q0(CSLICE,jj,3), N_BINS, N_DIM, size_x(3)/2, dndx(3))
              if(eval_rdm2_ring) then
                call vpi_eval_obdm_ring( rdm2_theta, q0(CSLICE,ii,:), q0(CSLICE,jj,:), N_BINS, N_DIM)
              end if
            end do
          end do
        end if
        if(eval_rdm22_ring) then
          do ii = N_PARTICLE-N_PARTICLE2+1, N_PARTICLE
            do jj = ii+1, N_PARTICLE
              call vpi_eval_obdm_ring( rdm22_theta, q0(CSLICE,ii,:), q0(CSLICE,jj,:), N_BINS, N_DIM)
            end do
          end do
        end if
        if( eval_off_diagonal ) then
          call vpi_eval_obdm_cut( obdm_x, q0(CSLICE,od_pnum,1), q0(CSLICE+1,od_pnum,1), N_BINS, N_DIM, size_x(1)/2, dndx(1))
          call vpi_eval_obdm_cut( obdm_y, q0(CSLICE,od_pnum,2), q0(CSLICE+1,od_pnum,2), N_BINS, N_DIM, size_x(2)/2, dndx(2))
          call vpi_eval_obdm_cut( obdm_z, q0(CSLICE,od_pnum,3), q0(CSLICE+1,od_pnum,3), N_BINS, N_DIM, size_x(3)/2, dndx(3))
          if(eval_obdm_ring) then
            call vpi_eval_obdm_ring( obdm_theta, q0(CSLICE,od_pnum,:), q0(CSLICE+1,od_pnum,:), N_BINS, N_DIM)
          endif
          if( eval_rotation ) then
            call vpi_eval_obdm_cut( obdm_rot_z, q_rot0(CSLICE,od_pnum,2), q_rot0(CSLICE+1,od_pnum,2), N_BINS, N_DIM, 1.0, real(N_BINS)/2.0)
          end if
        end if
        if( eval_nrdm ) then
          call vpi_eval_nrdm( nrdm_z, q0(CSLICE,:,3), q0(CSLICE+1,:,3) )
        end if
      end do

      if( eval_off_diagonal ) then
        write(13, "(3g12.3)", ADVANCE="NO") q0(CSLICE,1,:)
        write(13, "(3g12.3)", ADVANCE="NO") q0(CSLICE+1,1,:)
!        write(13, "(3g)", ADVANCE="NO") q0(CSLICE,2,:)
!        write(13, "(3g)") q0(CSLICE+1,2,:)
      end if
      close(13)

      write(my_fname,"(a8,i4.4)")"rho.dat.",my_rank
      open(12, file=my_fname, status="replace")
      do ii = 1,N_BINS
        write(12, "(6g12.3)") dxdn(:)*(ii-0.5) - size_x(:)/2.0, rho(ii,:)/( n_moves*dxdn(:)*N_PARTICLE )
      end do
      close(12)

      write(my_fname,"(a6,i4.4)")"U.dat.",my_rank
      open(12, file=my_fname,  POSITION="APPEND")
      t1 = sum(U_0(CSLICE,:))
      write(12, "(g20.10)")  t1/N_PARTICLE
      close(12)

      write(my_fname,"(a10,i4.4)")"U_avg.dat.",my_rank
      open(12, file=my_fname,  status="replace")
      do ii = 1, size(U_avg)
        write(12, "(2g30.15)")  ii*dtau, U_avg(ii)
      end do
      close(12)

      write(my_fname,"(a11,i4.4)")"Energy.dat.",my_rank
      open(12, file=my_fname,  POSITION="APPEND")
      write(12, "(5g30.15)")  i,E_left/nsteps,E_center/nsteps,E_right/nsteps,dE/nsteps
      close(12)

      write(my_fname,"(a10,i4.4)")"rho_r.dat.",my_rank
      open(12, file=my_fname, status="replace")
      do ii = 1, size(rho_r)
        t1 = drdn*(ii)
        t2 = drdn*(ii-1)
        t3 = 4*M_PI*(t1**3-t2**3)/3
        write(12, "(4g12.3)") drdn*(ii-0.5), rho_r(ii)/( t3*n_moves*N_PARTICLE )
      end do
      close(12)

      write(my_fname,"(a9,i4.4)")"gofr.dat.",my_rank
      open(12, file=my_fname, status="replace")
      do ii = 1,size(gofr)
        t0 = dxdn_GOFR*(ii-0.5)
        t1 = dxdn_GOFR*(ii-1)
        t2 = dxdn_GOFR*(ii)
        t3 = (4.*M_PI/3)*(t2**3 - t1**3)*N_PARTICLE*(N_PARTICLE+1)/2
        write(12, "(2g12.3)") t0, gofr(ii)/( n_moves*t3 )
      end do
      close(12)

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
              write(13, "(5g12.3)", ADVANCE="NO") dxdn_rot_xyz(1)*(ii-0.5) - size_x(1)/2.0, dxdn_rot_xyz(2)*(jj-0.5) - size_x(2)/2.0, dxdn_rot_xyz(3)*(kk-0.5) - size_x(3)/2.0, gof_rot_xyz_avg(ii,jj,kk)/( n_moves )
              write(13, *)
              do ll = 1,N_BINS_ROT
                write(12, "(5g12.3)", ADVANCE="NO") dxdn_rot_xyz(1)*(ii-0.5) - size_x(1)/2.0,  dxdn_rot_xyz(2)*(jj-0.5) - size_x(2)/2.0,  dxdn_rot_xyz(3)*(kk-0.5) - size_x(3)/2.0,  2.0*(ll-1)/(N_BINS_ROT-1) - 1.0, gof_rot_xyz(ii,jj,kk,ll)/( n_moves )
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

      write(my_fname,"(a10,i4.4)")"dxlxr.dat.",my_rank
      open(12, file=my_fname, status="replace")
      do ii = 1,size(dxlxr,1)
        write(12, "(6g24.16)", ADVANCE="NO") dxdn*(ii-1) - size_x/2.0, dxlxr(ii,:)/( n_moves*dxdn*N_PARTICLE )
        write(12, *)
      end do
      close(12)

      write(my_fname,"(a11,i4.4)")"rdm2_z.dat.",my_rank
      open(12, file=my_fname, status="replace")
      do ii = 1, N_BINS
        do jj = 1, N_BINS
          write(12, *) dxdn(3)*(ii-0.5) - size_x(3)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0,  rdm2_z(ii,jj)/( n_moves*dxdn(3)**2 )
        end do
        write(12, *)
      end do
      close(12)

      write(my_fname,"(a11,i4.4)")"obdm_z.dat.",my_rank
      open(12, file=my_fname, status="replace")
      do ii = 1, N_BINS
        do jj = 1, N_BINS
          write(12, *) dxdn(3)*(ii-0.5) - size_x(3)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0,  obdm_z(ii,jj)/( 2*n_moves*dxdn(3) )
        end do
        write(12, *)
      end do
      close(12)

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

      if(eval_obdm_ring) then
        write(my_fname,"(a14,i4.4)")"obdm_ring.dat.",my_rank
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

      write(my_fname,"(a11,i4.4)")"obdm_x.dat.",my_rank
      open(12, file=my_fname, status="replace")
      do ii = 1, N_BINS
        do jj = 1, N_BINS
          write(12, *) dxdn(1)*(ii-0.5) - size_x(1)/2.0, dxdn(1)*(jj-0.5) - size_x(1)/2.0,  obdm_x(ii,jj)/( 2*n_moves*dxdn(1) )
        end do
        write(12, *)
      end do
      close(12)

      write(my_fname,"(a11,i4.4)")"obdm_y.dat.",my_rank
      open(12, file=my_fname, status="replace")
      do ii = 1, N_BINS
        do jj = 1, N_BINS
          write(12, *) dxdn(2)*(ii-0.5) - size_x(2)/2.0, dxdn(2)*(jj-0.5) - size_x(2)/2.0,  obdm_y(ii,jj)/( 2*n_moves*dxdn(2) )
        end do
        write(12, *)
      end do
      close(12)

      if( use_eval_N_R ) then
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


      if(compute_full_density) then
        write(my_fname,"(a9,i4.4)")"odxy.dat.",my_rank
        open(23, file=my_fname, status="replace")
        do ii = 1, N_BINS
          do jj = 1, N_BINS
            do kk = 1, N_BINS
              write(23, *) dxdn(1)*(ii-0.5) - size_x(1)/2.0, dxdn(2)*(jj-0.5) - size_x(2)/2.0, dxdn(3)*(jj-0.5) - size_x(3)/2.0, full_density(ii,jj,kk)/( n_moves*N_PARTICLE*dxdn(1)*dxdn(2)*dxdn(3) )
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

      write(my_fname,"(a11,i4.4)")"bisect.dat.",my_rank
      open(10, file=my_fname, status="replace")
      write(10, *) "#", accr/(n_moves)
      do ii = 1,N_SLICE
        write(10, "(4g20.12)") (ii-1)*dtau, qsq_sum(ii,:)/(N_PARTICLE*n_moves)
      end do
      close(10)

      write(my_fname,"(a10,i4.4)")"phase.dat.",my_rank
      open(30, file=my_fname, status="replace")
      write(30, *) "#", accr/(n_moves)
      do ii = 1,N_SLICE
        write(30, "(i12,2g20.12)") ii, phase_0(ii), phase_sum_0(ii)/mcount(ii)
      end do
      close(30)

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

      write(my_fname,"(a10,i4.4)")"paths.dat.",my_rank
      open(10, file=my_fname, status="replace")
      do ii = 1,CSLICE
        write(10, "(i12)", ADVANCE="NO") ii
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

#ifdef USE_MPI
      time1 = MPI_Wtime()
#else
      time1 = day_timer()
#endif

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
        stop
      end if
  real(kind=b8), dimension( nbins , N_DIM ), intent(out):: rho

    ! Write random number seed
    call ru_write_seed_file(my_rank)
#ifdef USE_MPI
    call MPI_Finalize(ierr)
#endif

  contains

subroutine vpi_eval_density( rho, x, nbins, size_x, dndx )
  real(kind=b8), dimension( nbins , N_DIM ), intent(out):: rho
  real(kind=b8), dimension( N_PARTICLE , N_DIM ) :: x
  integer nbins
  real, dimension(N_DIM) :: size_x, dndx

end subroutine vpi_eval_density
  integer :: i,j

  real(kind=b8), dimension( : , : ), intent(out):: rho
    do j = 1, N_DIM
      bin = floor((x(i,j)+size_x(j))*dndx(j))+1
      if ( (bin .le. nbins) .and. (bin .ge. 1) ) then
        rho(bin,j) = rho(bin,j) + 1
      end if 
    end do
  end do
end subroutine vpi_eval_density

subroutine vpi_eval_dphase( rho, x, nbins, size_x, dndx )
  real(kind=b8), dimension( : , : ), intent(out):: rho
  real(kind=b8), dimension( : , : ) :: x
  integer nbins
  real :: size_x, dndx
  integer bin
end subroutine vpi_eval_dphase!}}}
  integer :: i,j

  real(kind=b8), dimension( : ), intent(out):: rho
    do j = 1, size(x,2)
      bin = floor((x(i,j)+size_x)*dndx)+1
      if ( (bin .le. nbins) .and. (bin .ge. 1) ) then
        rho(bin,j) = rho(bin,j) + 1
      end if 
    end do
  end do
end subroutine vpi_eval_dphase

subroutine vpi_eval_radial_density( rho, x, nbins, dndx )
  real(kind=b8), dimension( : ), intent(out):: rho
  real(kind=b8), dimension( : , : ) :: x
  integer nbins
  real :: dndx

  real(kind=b8) :: r
end subroutine vpi_eval_radial_density!}}}

  integer :: i,j
  real(kind=b8), dimension( nbins , N_DIM ), intent(out):: rho
  do i = 1, size(x,1)
    r = sqrt(dot_product(x(i,:),x(i,:)))
    bin = floor(r*dndx)+1
    if ( (bin .le. nbins) .and. (bin .ge. 1) ) then
      rho(bin) = rho(bin) + 1
    end if 
  end do
end subroutine vpi_eval_radial_density

subroutine vpi_eval_density_sp( rho, x, nbins, size_x, dndx )
  real(kind=b8), dimension( nbins , N_DIM ), intent(out):: rho
  real(kind=b8), dimension( N_DIM ) :: x
  integer nbins
end subroutine vpi_eval_density_sp!}}}
  integer bin

  real(kind=b8), dimension( : , : ), intent(out):: od

  do j = 1, N_DIM
    bin = anint((x(j)+size_x(j))*dndx(j))+1
    if ( (bin .le. nbins) .and. (bin .ge. 1) ) then
      rho(bin,j) = rho(bin,j) + 1
    end if 
  end do
end subroutine vpi_eval_density_sp

subroutine vpi_eval_optical_density( od, od_avg, od_sig, x, id1, id2, nbins, size_x, dndx )
  real(kind=b8), dimension( : , : ), intent(out):: od
  real(kind=b8), dimension( : , : ), intent(in):: od_avg
  real(kind=b8), dimension( : , : ), intent(out):: od_sig
  real(kind=b8), dimension( : , : ) :: x
  real(kind=b8), dimension( nbins , nbins ) :: tod
  integer nbins, id1, id2
  real, dimension(N_DIM) :: size_x, dndx
  integer xbin, zbin

  integer :: i

  tod = 0

end subroutine vpi_eval_optical_density!}}}
    xbin = anint((x(i,id1)+size_x(id1))*dndx(id1))+1
    if ( (xbin .le. nbins) .and. (xbin .ge. 1) ) then
  real(kind=b8), dimension( : , : , : ), intent(out):: od
      if ( (zbin .le. nbins) .and. (zbin .ge. 1) ) then
        tod(xbin,zbin) = tod(xbin,zbin) + 1
      end if
    end if 
  end do
  od = od + tod
  od_sig(:,:) = od_sig(:,:) + (tod(:,:) - od_avg(:,:))**2
end subroutine vpi_eval_optical_density

subroutine vpi_eval_full_density( od, x, nbins, size_x, dndx )
  real(kind=b8), dimension( : , : , : ), intent(out):: od
  real(kind=b8), dimension( : , : ) :: x
  integer :: nbins
  real, dimension(N_DIM) :: size_x, dndx

  integer :: xb, yb, zb

  integer :: i

  do i = 1, N_PARTICLE
end subroutine vpi_eval_full_density!}}}
    if ( (xb .le. nbins) .and. (xb .ge. 1) ) then
      yb = anint((x(i,2)+size_x(2))*dndx(2))+1
  real(kind=b8), dimension( : , : ), intent(out):: od
        zb = anint((x(i,3)+size_x(3))*dndx(3))+1
        if ( (zb .le. nbins) .and. (zb .ge. 1) ) then
          od(xb,yb,zb) = od(xb,yb,zb) + 1
        end if
      end if
    end if 
  end do
end subroutine vpi_eval_full_density

subroutine vpi_eval_optical_density2( od, od_avg, od_sig, x, nbins, size_x, dndx )
  real(kind=b8), dimension( : , : ), intent(out):: od
  real(kind=b8), dimension( : , : ), intent(in):: od_avg
  real(kind=b8), dimension( : , : ), intent(out):: od_sig
  real(kind=b8), dimension( : , : ) :: x
  real(kind=b8), dimension( nbins , nbins ) :: tod
  integer nbins
  real, dimension(N_DIM) :: size_x, dndx
  integer xbin, zbin

  integer :: i

  tod = 0

end subroutine vpi_eval_optical_density2!}}}
    xbin = anint((x(i,1)+size_x(1))*dndx(1))+1
    if ( (xbin .le. nbins) .and. (xbin .ge. 1) ) then
  real(kind=b8), dimension( : ), intent(inout):: N_R 
      if ( (zbin .le. nbins) .and. (zbin .ge. 1) ) then
        tod(xbin,zbin) = tod(xbin,zbin) + 1
      end if
    end if 
  end do
  od = od + tod
  od_sig(:,:) = od_sig(:,:) + (tod(:,:) - od_avg(:,:))**2
end subroutine vpi_eval_optical_density2

subroutine vpi_eval_N_R( N_R, x, islice, dslice )
  real(kind=b8), dimension( : ), intent(inout):: N_R 
  integer, intent(in):: islice, dslice
  real(kind=b8), dimension( : , : , : ), intent(in) :: x

  integer :: i,j, cnt
end subroutine vpi_eval_N_R!}}}

  do j = islice-dslice, islice+dslice
  real(kind=b8), dimension( : ), intent(inout):: N_zs, N_za
    do i = 1, N_PARTICLE
      if ( x(j,i,3) .gt. gw_z0 ) then
        cnt = cnt + 1
      end if 
    end do
    N_R(cnt) = N_R(cnt) + 1.0_b8/(2.0_b8*dslice+1.0_b8)
  end do
end subroutine vpi_eval_N_R

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
end subroutine vpi_eval_N_z2!}}}
      if ( x(i,1) .gt. 0 ) then
        br_cnt = br_cnt + 1
  real(kind=b8), dimension( N_SLICE, N_DIM ), intent(inout):: corr_x 
        bl_cnt = bl_cnt + 1
      end if
    end if 
  end do

  N_zs(tr_cnt+tl_cnt+1) = N_zs(tr_cnt+tl_cnt+1) + 1
  N_za((tr_cnt-tl_cnt)+N_PARTICLE2+1) = N_za((tr_cnt-tl_cnt)+N_PARTICLE2+1) + 1
end subroutine vpi_eval_N_z2

subroutine vpi_eval_corr_x( corr_x, move_start, move_end, x )
  real(kind=b8), dimension( N_SLICE, N_DIM ), intent(inout):: corr_x 
  integer :: move_start, move_end
end subroutine vpi_eval_corr_x!}}}

  integer :: i,j
  real(kind=b8), dimension( N_SLICE ), intent(inout):: corr_phase
  do i = move_start+1, move_end-1
    do j = 1, N_DIM
      corr_x(i,j) = corr_x(i,j) + sum( x(1,:,j)*x(i,:,j) )
      corr_x(N_SLICE-i+1,j) = corr_x(N_SLICE-i+1,j) + sum( x(N_SLICE,:,j)*x(N_SLICE-i+1,:,j) )
    end do
  end do

end subroutine vpi_eval_corr_x

subroutine vpi_eval_corr_phase( corr_phase, phase, move_start, move_end )
end subroutine vpi_eval_corr_phase
  real(kind=b8), dimension( N_SLICE ), intent(in) :: phase
  integer :: move_start, move_end
  real(kind=b8), dimension( N_SLICE ), intent(inout):: corr_r
  integer :: i

  do i = move_start+1, move_end-1
    corr_phase(i) = corr_phase(i) + phase(1)*phase(i)
    corr_phase(N_SLICE-i+1) = corr_phase(N_SLICE-i+1) + phase(N_SLICE)*phase(N_SLICE-i+1)
  end do

end subroutine vpi_eval_corr_phase!}}}

subroutine vpi_eval_corr_r( corr_r, move_start, move_end, x )
  real(kind=b8), dimension( N_SLICE ), intent(inout):: corr_r
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x
  integer :: move_start, move_end

  integer :: i
end subroutine vpi_eval_corr_r!}}}

  r1 = sum( (x(1,:,1)**2 + x(1,:,2)**2) ) - N_PARTICLE
  real(kind=b8), dimension( N_SLICE, 2 ), intent(inout):: corr_xz
  do i = move_start+1, move_end-1
    tmp = sum( ( x(i,:,1)**2 + x(i,:,2)**2 ) ) - N_PARTICLE
    corr_r(i) = corr_r(i) + r1*tmp
    tmp = sum( ( x(N_SLICE-i+1,:,1)**2 + x(N_SLICE-i+1,:,2)**2 ) ) - N_PARTICLE
    corr_r(N_SLICE-i+1) = corr_r(N_SLICE-i+1) + rM*tmp
  end do

end subroutine vpi_eval_corr_r

subroutine vpi_eval_corr_xz( corr_xz, move_start, move_end, x )
  real(kind=b8), dimension( N_SLICE, 2 ), intent(inout):: corr_xz
  real(kind=b8), dimension( N_SLICE, N_PARTICLE, N_DIM ) :: x
  integer :: move_start, move_end
end subroutine vpi_eval_corr_xz!}}}
  integer :: i,j
  real(kind=b8) :: r1,rM,tmp
  integer :: nbins
  do i = move_start+1, move_end-1
    do j = 1, 2
      corr_xz(i,j) = corr_xz(i,j) + sum( x(1,:,j)*x(1,:,3)*x(i,:,j)*x(i,:,3) )
      corr_xz(N_SLICE-i+1,j) = corr_xz(N_SLICE-i+1,j) + sum( x(N_SLICE,:,j)*x(N_SLICE,:,3)*x(N_SLICE-i+1,:,j)*x(N_SLICE-i+1,:,3) )
    end do
  end do

end subroutine vpi_eval_corr_xz

subroutine vpi_eval_gofr( gofr, xij2, nbins, dndx )
  integer :: nbins
  real :: dndx
  real(kind=b8), dimension( nbins ), intent(inout):: gofr 
  real(kind=b8), dimension( N_PARTICLE , N_PARTICLE ) :: xij2

  integer :: i,j
  integer :: rbin
  real :: r
end subroutine vpi_eval_gofr!}}}
  do i = 1, N_PARTICLE
    do j = i + 1, N_PARTICLE
  integer :: nbins_xyz, nbins_rot, nbins_rr
      rbin = floor( r * dndx ) + 1
      if( (rbin .ge. 1) .and. (rbin .le. nbins) ) then
        gofr(rbin) = gofr(rbin) + 1
      end if
    end do
  end do

end subroutine vpi_eval_gofr

subroutine vpi_eval_gof_rot( gof_rot, gof_rot_xyz, gof_rot_xyz_avg, gofr_rot, q, q_rot, xij2, slice, nbins_rr, dndx_rr, nbins_xyz, xsize, dndx, nbins_rot, dndx_rot )
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
        if( (bx .ge. 1) .and. (bx .le. nbins_xyz) ) then
end subroutine vpi_eval_gof_rot!}}}
            if( (bz .ge. 1) .and. (bz .le. nbins_xyz) ) then
              gof_rot_xyz(bx,by,bz,br) = gof_rot_xyz(bx,by,bz,br) + 1
  real(kind = b8) :: lngfn0, lngfn1
            end if
          end if
        end if
      end if
    end do
  end do

end subroutine vpi_eval_gof_rot

function vpi_accept_path( lngfn0, lngfn1, dphase ) result( accept )
  real(kind = b8) :: lngfn0, lngfn1
  real(kind = b8) :: dphase
  logical :: accept

  real :: Pa, Ptest 
end function vpi_accept_path!}}}
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

!function eval_local_energy(q, tfun, tlap, pfun, np, ndim ) result ( E )
end program Test_VPI!}}}
!  real, dimension( : , : ) :: q
!  real :: tfun, pfun
!  real :: E
!
!  real :: hbsq2m = 1
!
!  E = (-hbsq2m*tlap(q, np,ndim) + pfun(q, np,ndim))/tfun(q)
!  
!end function eval_local_energy

end program Test_VPI

