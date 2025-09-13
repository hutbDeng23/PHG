FILES =	Makefile functions.h poisson.c simplest.c .cvsignore *-plot.sh jobs.sh \
	maxwell*.c heat.c non-smooth.c eigen.c elastic.c navier-stokes.[ch] \
	navier-stokes-unsteady.c robin.c petsc-asm.options phg-asm.options \
	asp_gen_options specht.c ipdg.c interface.c

default: lib simplest poisson maxwell maxwell-complex heat eigen elastic \
	maxwell-eigen  non-smooth navier-stokes robin

all: lib poisson simplest maxwell maxwell-complex heat non-smooth eigen \
	elastic maxwell-eigen navier-stokes robin

examples.zip: ${FILES}
	@zip -9 -u -y $@ $^

clean:
	-/bin/rm -f *.o *.tmp plot.eps PI[0-9]* \
		*.G.* *.A.* *.Aalpha.* *.Abeta.* \
		*.Gx.* *.Gy.* *.Gz.* *.x.* *.y.* *.z.* *.b.* *.x0.*

distclean: clean
	-/bin/rm -f poisson simplest maxwell heat mesh*.* eigen elastic \
		maxwell-complex maxwell-eigen maxwell-eigen1 non-smooth \
		navier-stokes navier-stokes-unsteady robin specht \
		robin ipdg interface

lib:
	@(cd ..; $(MAKE) lib)

include ../Makefile.inc

poisson.o:	../src/libphg${LIB_SUFFIX} functions.h poisson.c
simplest.o:	../src/libphg${LIB_SUFFIX} functions.h simplest.c
robin.o:	../src/libphg${LIB_SUFFIX} functions.h robin.c
maxwell.o:	../src/libphg${LIB_SUFFIX} maxwell.c
heat.o:		../src/libphg${LIB_SUFFIX} heat.c
non-smooth.o:	../src/libphg${LIB_SUFFIX} non-smooth.c
eigen.o:	../src/libphg${LIB_SUFFIX} eigen.c
maxwell-eigen.o: ../src/libphg${LIB_SUFFIX} maxwell-eigen.c
elastic.o:	../src/libphg${LIB_SUFFIX} elastic.c
maxwell-complex.o: ../src/libphg${LIB_SUFFIX} maxwell-complex.c
navier-stokes.o: ../src/libphg${LIB_SUFFIX} navier-stokes.c navier-stokes.h
navier-stokes-unsteady.o: ../src/libphg${LIB_SUFFIX} navier-stokes-unsteady.c navier-stokes.h
robin: robin.o
ipdg: ipdg.o
interface: interface.o

.PHONY: default all clean distclean lib
