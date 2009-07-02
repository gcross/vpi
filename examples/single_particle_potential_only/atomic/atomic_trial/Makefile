VPI_ROOT = ../../../..

OBJS =  \
	${VPI_ROOT}/src/potentials/single_particle/sp_atomic_potential.o \
	${VPI_ROOT}/src/potentials/two_body/tb_null_potential.o \
	${VPI_ROOT}/src/potentials/rotational/rot_undefined_potential.o \
	${VPI_ROOT}/src/trials/single_particle/sp_atomic_trial.o \
	${VPI_ROOT}/src/trials/jastrow/jas_independent_trial.o \
	${VPI_ROOT}/lib/vpi_main.o \
	${VPI_ROOT}/lib/libvpi.a \

vpi: ${OBJS}
	mpif90 ${OBJS} -o vpi
