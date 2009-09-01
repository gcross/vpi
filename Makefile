#@+leo-ver=4-thin
#@+node:gcross.20090806151612.1842:@thin Makefile
#@@language Makefile
#@@tabwidth 4
all: lib/vpi.so

include paths.mk

CC = gcc
FC = gfortran
FLAGS = -fPIC -O3
CFLAGS = -I ${NUMPYDIR}/include -I ${PYTHONINCDIR} -I src
FFLAGS = -cpp -fbounds-check -M src

SOURCES = \
  src/angular_momentum.f95 \
  src/timers.f90 \
  src/rand_utils.f95 \
  src/gfn.f95 \
  src/sample.f90 \
  src/xij.f95 \
  src/observables.f95 \
  src/thermalize.f90 \
  src/histograms.f95 \
  src/lattice.f90 \
  src/harmonic_oscillator_3d.f95 \
  src/hard_sphere_interaction.f95 \
  src/leonard_jones_interaction.f95

OBJS = \
  obj/constants.o \
  obj/angular_momentum.o \
  obj/timers.o \
  obj/rand_utils.o \
  obj/gfn.o \
  obj/sample.o \
  obj/xij.o \
  obj/observables.o \
  obj/thermalize.o \
  obj/histograms.o \
  obj/lattice.o \
  obj/harmonic_oscillator_3d.o \
  obj/hard_sphere_interaction.o \
  obj/leonard_jones_interaction.o \
  obj/vpimodule.o \
  obj/vpi-f2pywrappers2.o \
  obj/fortranobject.o

src/vpimodule.c src/vpi-f2pywrappers2.f90: ${SOURCES} Makefile
	f2py ${SOURCES} -m vpi
	mv vpimodule.c vpi-f2pywrappers2.f90 src

obj/vpimodule.o: src/vpimodule.c
	${CC} ${FLAGS} ${CFLAGS} -c $< -o $@

obj/%.o: src/%.c Makefile
	${CC} ${FLAGS} ${CFLAGS} -c $< -o $@
obj/%.o: src/%.f95 Makefile
	${FC} ${FLAGS} ${FFLAGS} -c $< -o $@
obj/%.o: src/%.f90 Makefile
	${FC} ${FLAGS} ${FFLAGS} -c $< -o $@

lib/vpi.so: ${OBJS}
	gfortran ${PYTHONLINK} -shared -o lib/vpi.so -lgfortran -lblas ${OBJS}

clean:
	rm -f obj/* lib/* src/*.mod src/vpimodule.c src/vpi-f2pywrappers2.f90
#@-node:gcross.20090806151612.1842:@thin Makefile
#@-leo
