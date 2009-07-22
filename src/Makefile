MODULE_PATH = potentials/single_particle:potentials/two_body:potentials/rotational:trials/single_particle:trials/jastrow:.

################################################################################
# AIX

#FC = xlf90
#FC = xlf90_r
#SPRNG_INCDIR = /home/1/huang/sprng2.0/include
#SPRNG_LIBDIR = /home/1/huang/sprng2.0/lib
#F90FLAGS = -qsuffix=cpp=f90 -Q -O3 -qarch=auto -qtune=auto -qhot -qfloat=hssngl
#F90FLAGS = -qsuffix=cpp=f90 -Q -O3 -qarch=auto -qtune=auto -qhot -qfloat=hssngl -WF,-DQMC_SMP -I${SPRNG_INCDIR} -L${SPRNG_LIBDIR} -qsmp=noauto
################################################################################
# Blue Horizon

#FC = mpxlf90
#FC = mpxlf90_r
#SPRNG_INCDIR = /paci/ucb/ux450650/sprng2.0/include
#SPRNG_LIBDIR = /paci/ucb/ux450650/sprng2.0/lib
#F90FLAGS = -qsuffix=cpp=f90 -Q -O3 -qarch=pwr3 -qtune=pwr3 -qhot -qfloat=hssngl -WF,-DQMC_MPI -I${SPRNG_INCDIR} -L${SPRNG_LIBDIR} -L/usr/local/apps/mass
#F90FLAGS = -qsuffix=cpp=f90 -Q -O3 -qarch=pwr3 -qtune=pwr3 -qhot -qfloat=hssngl -WF,-DQMC_MIX -I${SPRNG_INCDIR} -L${SPRNG_LIBDIR} -qsmp=noauto -L/usr/local/apps/mass -bmaxdata:0x80000000 -bmaxstack:0x10000000
#LIBS = -lmass
################################################################################
# Univac
# Note that for PGF90, -Msse and -Mprefetch causes problems on aluminum...

#FC = mpif77
#MYFFLAGS = -Mpreprocess -Minline -Mcache_align -Minfo=inline,loop,opt -Minform=inform -Munroll=c:6
### BEGIN NOTE!!! fastsse flag BREAKS THE CODE on UNIVAC!!!!
#MYFFLAGS = -Mprof=mpi -qp -Mpreprocess -fastsse -Minline -Mcache_align -Minfo=inline,loop,opt -Minform=inform -Munroll=c:6
#### END NOTE!!!
#MYFFLAGS = -g -pg -Mbounds  -Minform=inform
#MYFFLAGS = -g -Minform=inform
################################################################################
# DataStar
#
#FC = mpxlf90_r
#MYFFLAGS = -q64 -qsuffix=cpp=f90 -O3 -qstrict -qarch=pwr4 -qtune=pwr4 

#gfortran
#FC = mpigfortran
#F90FLAGS =  -ffree-line-length-none -x f95-cpp-input
#MYFFLAGS =  -O6 -funroll-loops -frerun-loop-opt -fomit-frame-pointer -finline-functions  -D USE_IMAGE_HSGFN
#MYFFLAGS =  -O6 -funroll-loops -frerun-loop-opt -fomit-frame-pointer -finline-functions  -D USE_MPI  -D USE_IMAGE_HSGFN
#MYFFLAGS =  -O6 -g -pg -funroll-loops -frerun-loop-opt -finline-functions  -D USE_MPI  -D USE_IMAGE_HSGFN
#MYFFLAGS =  -O6 -pg -funroll-loops -frerun-loop-opt -finline-functions  -D USE_MPI
#MYFFLAGS =   -g -fbounds-check -ffpe-trap=invalid -D USE_MPI  -D USE_IMAGE_HSGFN

#g95
#FC = g95
#F90FLAGS =  -ffree-line-length-huge -cpp
#MYFFLAGS =  -O6 -funroll-loops -frerun-loop-opt -fomit-frame-pointer -finline-functions 
#MYFFLAGS =  -g -DUSE_IMAGE_HSGFN
#MYLIBS = 

#ifort
FC = mpif90
F90FLAGS = -cpp
#MYFFLAGS =  -O3 -funroll-loops -fomit-frame-pointer -finline-functions -D USE_MPI  -D USE_IMAGE_HSGFN -D EXPERIMENTAL
#MYFFLAGS =  -p -O3 -funroll-loops -fomit-frame-pointer -finline-functions -D USE_MPI  -D USE_IMAGE_HSGFN -prof-genx 


MYFFLAGS =  -O3 -funroll-loops -fomit-frame-pointer -finline-functions -fimplicit-none -D USE_MPI -D USE_IMAGE_HSGFN -Jmods -Imods -x f95-cpp-input
#MYFFLAGS =  -g -pg -fimplicit-none -D USE_MPI -D USE_IMAGE_HSGFN -Jmods -Imods -x f95-cpp-input


#MYFFLAGS =  -g -D USE_MPI -D USE_IMAGE_HSGFN
#zeus
#MYLIBS =  -Wl,-rpath=/usr/local/tools/mkl/lib/ -L/usr/local/tools/mkl/lib/ -lmkl -lmkl_lapack64 -lguide -lpthread
MYLIBS =  -Wl -llapack -lpthread
#MYLIBS = -p -L/usr/local/intel/mkl/lib/em64t -lmkl -lmkl_lapack64 -lguide -lpthread

#pathscale
#FC = mpipathf90
#F90FLAGS = -cpp 
#MYFFLAGS =   -O3 -m64 -align128 -apo -D USE_MPI  -D USE_IMAGE_HSGFN -D EXPERIMENTAL
#MYFFLAGS =   -O3 -m64 -align128 -D USE_MPI  -D USE_IMAGE_HSGFN
#MYFFLAGS =   -g -C
#MYLIBS =   -I/usr/local/include/acml/pathscale -L/usr/local/lib/acml/pathscale -lacml

#thunder
#MYLIBS = -L/usr/local/intel/mkl/lib/64 -lmkl -lmkl_lapack64 -lguide -lpthread

PROGRAM_OBJS = kinds.o constants.o vpi_defines.o vpi_rand_utils.o vpi_xij.o vpi_aziz.o vpi_gfn.o vpi_bisect.o vpi_bbridge.o vpi_sample.o vpi_lattice.o vpi_coordinate_utils.o vpi_phase_maps.o vpi_rotation.o timers.o vpi_main.o

SINGLE_PARTICLE_POTENTIAL_OBJS = \
	potentials/single_particle/sp_null_potential.o \
	potentials/single_particle/sp_rectangular_box_potential.o \
	potentials/single_particle/sp_cylindrical_box_potential.o \
	potentials/single_particle/sp_harmonic_oscillator_1D_potential.o \
	potentials/single_particle/sp_harmonic_oscillator_3D_potential.o \
	potentials/single_particle/sp_harmonic_oscillator_overlapping_potential.o \
	potentials/single_particle/sp_wlink_potential.o \
	potentials/single_particle/sp_morse_oscillator_potential.o \
	potentials/single_particle/sp_G_well_potential.o \
	potentials/single_particle/sp_D_well_potential.o \
	potentials/single_particle/sp_N_well_potential.o \
	potentials/single_particle/sp_annulus_common.o \
	potentials/single_particle/sp_annulus_version_1_potential.o \
	potentials/single_particle/sp_annulus_version_2_potential.o \
	potentials/single_particle/sp_annulus_version_3_potential.o \
	potentials/single_particle/sp_atomic_potential.o \
	potentials/single_particle/sp_atomic_overlapping_potential.o \
	potentials/single_particle/sp_lattice_2D_with_HO_potential.o \
	potentials/single_particle/sp_lattice_3D_with_HO_potential.o \


TWO_BODY_POTENTIAL_OBJS = \
	potentials/two_body/tb_null_potential.o \
	potentials/two_body/tb_charge_potential.o \
	potentials/two_body/tb_lj_potential.o \
	potentials/two_body/tb_sc_potential.o \
	potentials/two_body/tb_hard_sphere_potential.o \
#	potentials/two_body/tb_sa_potential.o \
	potentials/two_body/tb_dimer_Z_polarized_potential.o \
	potentials/two_body/tb_dimer_Z_polarized_with_hard_sphere_potential.o \

ROTATIONAL_POTENTIAL_OBJS = \
	potentials/rotational/rot_dimer_potential.o \
	potentials/rotational/rot_undefined_potential.o \

POTENTIAL_OBJS = \
	${SINGLE_PARTICLE_POTENTIAL_OBJS} \
	${TWO_BODY_POTENTIAL_OBJS} \
	${ROTATIONAL_POTENTIAL_OBJS} \

SINGLE_BODY_TRIAL_OBJS = \
	trials/single_particle/sp_trial_numeric_differentiator.o \
	trials/single_particle/sp_atomic_trial.o \
	trials/single_particle/sp_He3plus_trial.o \
	trials/single_particle/sp_box_trial.o \
	trials/single_particle/sp_harmonic_3D_trial.o \
	trials/single_particle/sp_anharmonic_trial.o \
	trials/single_particle/sp_morse_trial_common.o \
	trials/single_particle/sp_morse_version_1_trial.o \
	trials/single_particle/sp_morse_version_2_trial.o \
	trials/single_particle/sp_D_well_trial_common.o \
	trials/single_particle/sp_D_well_version_1_trial.o \
	trials/single_particle/sp_D_well_version_2_trial.o \
	trials/single_particle/sp_lattice_trial.o \
#	trials/single_particle/sp_N_well_trial.o \

JASTROW_TRIAL_OBJS = \
	trials/jastrow/jas_trial_numeric_differentiator.o \
	trials/jastrow/jas_hard_sphere_trial.o \
	trials/jastrow/jas_lj_trial.o \
	trials/jastrow/jas_charge_trial.o \
	trials/jastrow/jas_independent_trial.o \

TRIAL_OBJS = ${SINGLE_BODY_TRIAL_OBJS} ${JASTROW_TRIAL_OBJS}

OBJS = ${POTENTIAL_OBJS} ${TRIAL_OBJS} ${PROGRAM_OBJS}

INTERFACE_OBJS = \
	vpi_single_particle_potential.o \
	vpi_two_body_potential.o \
	vpi_rotational_potential.o \
	vpi_single_particle_trial.o \
	vpi_jastrow_trial.o

%.o: %.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c $< -o $@

all: kinds.o constants.o vpi_defines.o vpi_xij.o \
	potentials/single_particle/sp_annulus_common.o \
	trials/single_particle/sp_morse_trial_common.o \
	trials/single_particle/sp_D_well_trial_common.o \
	trials/single_particle/sp_trial_numeric_differentiator.o \
	trials/jastrow/jas_trial_numeric_differentiator.o \
	../lib/libvpi.a ${OBJS} \
	../lib/vpi_main.o \

vpi_main.f90 : ${INTERFACE_OBJS}

.Default: all

vpi_rand_utils.f90: timers.o

# Final executables
vpi :: ${OBJS}
	${FC} ${OBJS}  -o ../bin/vpi ${MYLIBS}

../lib/libvpi.a : ${PROGRAM_OBJS} ${INTERFACE_OBJS}
	ar r ../lib/libvpi.a ${PROGRAM_OBJS}

../lib/vpi_main.o: vpi_main.o
	cp vpi_main.o ../lib

hnode :: vpi_hnode_defines.o vpi_rand_utils.o vpi_xij.o vpi_hnode_potential.o vpi_gfn.o vpi_hnode_trial_func.o vpi_bisect.o vpi_hnode_lattice.o hnode_main.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} vpi_hnode_defines.o vpi_rand_utils.o vpi_xij.o vpi_hnode_potential.o vpi_gfn.o vpi_hnode_trial_func.o vpi_bisect.o vpi_hnode_lattice.o vpi_obdm.o hnode_main.f90 -o ../bin/hnode-test

clean : 
	rm -f *.mod
	rm -f */*/*.mod
	rm -f *.o
	rm -f */*/*.o
	rm -f *.lst
	rm -f *.i
	rm -f *.inc
	rm -f *.vo
	rm -f *~
	rm -f *.stb
