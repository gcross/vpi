all: vpi

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
MYFFLAGS =  -O3 -funroll-loops -fomit-frame-pointer -finline-functions -D USE_MPI -D USE_IMAGE_HSGFN
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


# Modules
kinds.o : kinds.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c kinds.f90
vpi_defines.o : vpi_defines.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_defines.f90
vpi_timers.o : timers.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c timers.f90
vpi_hnode_defines.o : vpi_hnode_defines.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_hnode_defines.f90
vpi_rand_utils.o : vpi_rand_utils.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c timers.f90 vpi_rand_utils.f90
vpi_aziz.o : vpi_aziz.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_aziz.f90
vpi_potential.o : vpi_potential.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_potential.f90
vpi_hnode_potential.o : vpi_hnode_potential.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_hnode_potential.f90
vpi_trial_func.o : vpi_trial_func.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_trial_func.f90
vpi_hnode_trial_func.o : vpi_hnode_trial_func.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_hnode_trial_func.f90
vpi_gfn.o : vpi_gfn.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_gfn.f90
vpi_bisect.o : vpi_bisect.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_bisect.f90
vpi_bbridge.o : vpi_bbridge.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_bbridge.f90
vpi_sample.o : vpi_sample.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_sample.f90
vpi_lattice.o : vpi_lattice.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_lattice.f90
vpi_hnode_lattice.o : vpi_hnode_lattice.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_hnode_lattice.f90
vpi_xij.o : vpi_xij.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_xij.f90
vpi_obdm.o : vpi_obdm.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_obdm.f90
vpi_coordinate_utils.o : vpi_coordinate_utils.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_coordinate_utils.f90
vpi_phase_maps.o : vpi_phase_maps.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_phase_maps.f90
vpi_rotation.o : vpi_rotation.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_rotation.f90
vpi_read_input.o : vpi_read_input.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c vpi_read_input.f90
nl2sol.o : nl2sol.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} -c nl2sol.f90

#main
vpi_main.o : vpi_main.f90
	 ${FC} ${F90FLAGS} ${MYFFLAGS} -c  vpi_main.f90

# Final executables
vpi :: kinds.o vpi_defines.o vpi_rand_utils.o vpi_xij.o vpi_aziz.o vpi_potential.o vpi_gfn.o vpi_trial_func.o vpi_bisect.o vpi_bbridge.o vpi_sample.o vpi_lattice.o vpi_obdm.o vpi_coordinate_utils.o vpi_phase_maps.o vpi_rotation.o vpi_read_input.o timers.o vpi_main.o
	${FC} kinds.o vpi_defines.o vpi_rand_utils.o vpi_xij.o vpi_aziz.o vpi_potential.o vpi_gfn.o vpi_trial_func.o vpi_bisect.o vpi_bbridge.o vpi_sample.o vpi_lattice.o vpi_obdm.o vpi_coordinate_utils.o vpi_phase_maps.o vpi_rotation.o vpi_read_input.o timers.o vpi_main.o  -o ../bin/vpi-hoz-independent ${MYLIBS}

hnode :: vpi_hnode_defines.o vpi_rand_utils.o vpi_xij.o vpi_hnode_potential.o vpi_gfn.o vpi_hnode_trial_func.o vpi_bisect.o vpi_hnode_lattice.o vpi_obdm.o hnode_main.f90
	${FC} ${F90FLAGS} ${MYFFLAGS} vpi_hnode_defines.o vpi_rand_utils.o vpi_xij.o vpi_hnode_potential.o vpi_gfn.o vpi_hnode_trial_func.o vpi_bisect.o vpi_hnode_lattice.o vpi_obdm.o hnode_main.f90 -o ../bin/hnode-test

clean : 
	rm -f *.mod
	rm -f *.o
	rm -f *.lst
	rm -f *.i
	rm -f *.inc
	rm -f *.vo
	rm -f *~
	rm -f *.stb
