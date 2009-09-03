#@+leo-ver=4-thin
#@+node:gcross.20090806151612.1842:@thin Makefile
#@@language Makefile
#@@tabwidth 4
all: lib/vpif.so

include paths.mk

CC = gcc
FC = gfortran
FLAGS = -fPIC -O3
CFLAGS = -I ${NUMPYDIR}/include -I ${PYTHONINCDIR} -I src
FFLAGS = -cpp -fbounds-check -M src

SOURCES = \
  src/numeric_differentiation.f95 \
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
  src/harmonic_oscillator.f95 \
  src/hard_sphere_interaction.f95 \
  src/leonard_jones_interaction.f95

OBJS = \
  obj/constants.o \
  obj/xij.o \
  obj/numeric_differentiation.o \
  obj/angular_momentum.o \
  obj/timers.o \
  obj/rand_utils.o \
  obj/gfn.o \
  obj/sample.o \
  obj/observables.o \
  obj/thermalize.o \
  obj/histograms.o \
  obj/lattice.o \
  obj/harmonic_oscillator.o \
  obj/hard_sphere_interaction.o \
  obj/leonard_jones_interaction.o \
  obj/vpifmodule.o \
  obj/vpif-f2pywrappers2.o \
  obj/fortranobject.o

src/vpifmodule.c src/vpif-f2pywrappers2.f90: ${SOURCES} Makefile
	f2py ${SOURCES} -m vpif
	mv vpifmodule.c vpif-f2pywrappers2.f90 src

obj/vpifmodule.o: src/vpifmodule.c
	${CC} ${FLAGS} ${CFLAGS} -c $< -o $@

obj/%.o: src/%.c Makefile
	${CC} ${FLAGS} ${CFLAGS} -c $< -o $@
obj/%.o: src/%.f95 Makefile
	${FC} ${FLAGS} ${FFLAGS} -c $< -o $@
obj/%.o: src/%.f90 Makefile
	${FC} ${FLAGS} ${FFLAGS} -c $< -o $@

lib/vpif.so: ${OBJS}
	gfortran ${PYTHONLINK} -shared -o lib/vpif.so -lgfortran -lblas ${OBJS}

clean:
	rm -f obj/* lib/* src/*.mod src/vpifmodule.c src/vpif-f2pywrappers2.f90
#@nonl
#@-node:gcross.20090806151612.1842:@thin Makefile
#@-leo
