#@+leo-ver=4-thin
#@+node:gcross.20090806151612.1842:@thin Makefile
#@@language Makefile
#@@tabwidth 4
all: lib/vpif.so lib/vpi.py

include paths.mk
include options.mk

SOURCES = \
  src/utils/erfn.f95 \
  src/utils/timers.f95 \
  src/utils/rand_utils.f95 \
  src/utils/xij.f95 \
  src/utils/numeric_differentiation.f95 \
  src/measurement/observables.f95 \
  src/measurement/histograms.f95 \
  src/physics/angular_momentum.f95 \
  src/physics/gfn.f95 \
  src/physics/hard_sphere_interaction.f95 \
  src/physics/harmonic_oscillator.f95 \
  src/physics/lennard_jones_interaction.f95 \
  src/path-integral/lattice.f95 \
  src/path-integral/sample.f95 \
  src/path-integral/thermalize.f95

OBJS = \
  obj/utils/constants.o \
  obj/utils/erf.o \
  obj/utils/erfn.o \
  obj/utils/timers.o \
  obj/utils/rand_utils.o \
  obj/utils/xij.o \
  obj/utils/numeric_differentiation.o \
  obj/measurement/observables.o \
  obj/measurement/histograms.o \
  obj/physics/gfn.o \
  obj/physics/angular_momentum.o \
  obj/physics/hard_sphere_interaction.o \
  obj/physics/harmonic_oscillator.o \
  obj/physics/lennard_jones_interaction.o \
  obj/path-integral/lattice.o \
  obj/path-integral/sample.o \
  obj/path-integral/thermalize.o \
  obj/wrappers/vpifmodule.o \
  obj/wrappers/vpif-f2pywrappers2.o \
  obj/wrappers/fortranobject.o

src/wrappers/vpifmodule.c src/wrapers/vpif-f2pywrappers2.f90: ${SOURCES} Makefile
	${F2PY} ${SOURCES} -m vpif --build-dir src/wrappers

obj/%.o: src/%.c Makefile
	${CC} ${FLAGS} ${CFLAGS} -c $< -o $@
obj/%.o: src/%.f95 Makefile
	${FC} ${FLAGS} ${FFLAGS} -c $< -o $@
obj/%.o: src/%.f90 Makefile
	${FC} ${FLAGS} ${FFLAGS} -c $< -o $@
obj/%.o: src/%.f Makefile
	${FC} ${FLAGS} ${FFLAGS} -c $< -o $@

lib/vpif.so: ${OBJS}
	gfortran ${PYTHONLINK} ${FLAGS} -shared -o lib/vpif.so -lgfortran -lblas -llapack ${OBJS}

lib/vpi.py: src/vpi.py
	cp $< $@

clean:
	rm -f obj/*/* lib/* src/*/*.mod src/wrappers/vpifmodule.c src/wrappers/vpif-f2pywrappers2.f90
#@-node:gcross.20090806151612.1842:@thin Makefile
#@-leo
