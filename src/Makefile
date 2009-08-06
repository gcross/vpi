#@+leo-ver=4-thin
#@+node:gcross.20090806151612.1842:@thin Makefile
#@@language Makefile
all: vpi.so

NUMPYDIR = ${HOME}/local/lib/python2.6/site-packages/numpy/core

CC = gcc
FC = gfortran
FLAGS = -fPIC
CFLAGS = -I ${NUMPYDIR}/include -I ${HOME}/local/include/python2.6
FFLAGS = -cpp

SOURCES = \
	timers.f90 \
	rand_utils.f90 \
	sample.f90 \
	xij.f95

OBJS = \
  constants.o \
	timers.o \
	rand_utils.o \
	sample.o \
	xij.o \
	vpimodule.o \
	vpi-f2pywrappers2.o \
	fortranobject.o

vpimodule.c vpi-f2pywrappers2.f90: ${SOURCES} Makefile
	f2py ${SOURCES} -m vpi

%.o: %.c Makefile
	${CC} ${FLAGS} ${CFLAGS} -c $< -o $@
%.o: %.f95 Makefile
	${FC} ${FLAGS} ${FFLAGS} -c $< -o $@
%.o: %.f90 Makefile
	${FC} ${FLAGS} ${FFLAGS} -c $< -o $@

vpi.so: ${OBJS}
	gfortran -shared -o vpi.so -lgfortran ${OBJS}
#@-node:gcross.20090806151612.1842:@thin Makefile
#@-leo
