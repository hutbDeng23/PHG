# $Id: Makefile,v 1.92 2022/09/13 01:53:45 zlb Exp $

default:
all:

quadrule-dist:
	sh quadrule-dist.sh

lib:
	@(cd ../src; $(MAKE))

CLEAN_FILES = abaqus cpp-filter density dof_info mirror quadrule tahoe2medit \
		grid2vtk *.a *.so *.o *.bak core* *.sage.py grid2vtk
clean:
	-/bin/rm -rf ${CLEAN_FILES}

distclean: clean
	-/bin/rm -f phgdoc

cpp-filter: cpp-filter.c
	${CC} ${CFLAGS} -o $@ $^

abaqus: abaqus.o ../src/elem-info.c ../include/phg/elem-info.h
	${LINKER} ${LDFLAGS} -I../include -o $@ abaqus.o ../src/elem-info.c ${LIBS}
abaqus.o: abaqus.c ../include/phg/elem-info.h

###############################################################################
# Note: leave the line below unchanged for quadrule-dist.sh
# Begin Makefile for quadrule
#-----------------------------------------------------------

# NONLINEAR_SOLVER={NEWTON,MINPACK,PETSC,TENSOLVE,CCMNF,GSLSIMAN}
CDEFS0=-DNONLINEAR_SOLVER=MINPACK
CDEFS = -DUSE_MPI=1 ${CDEFS0}

MPFR_CPPFLAGS = -DUSE_MPFR=1
MPFR_LIB = -lmpfr -lgmp

MINPACK_LIB = minpack/libminpack.a
GSL_CPPFLAGS =
GSL_LIB = #-lgsl

GAUL_CPPFLAGS =
GAUL_LIB = #-lgaul -lgaul_util

# PETSc
# Headers (-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include)
PETSC_CPPFLAGS = $(shell sh petsc-configs.sh CPPFLAGS)
# Libs (-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -llapack -lblas -lpthread)
PETSC_LIBS = $(shell sh petsc-configs.sh LDFLAGS RPATHS LIBS)

CC = mpicc
CFLAGS = -Wall -O2 -g
CPPFLAGS = ${PETSC_CPPFLAGS}
F77 = mpif77
# Note: with GCC 10 '-std=legacy' is needed, otherwise 'Rank mismatch in argument' errors are produced when compiling the original toms-768/Src/Dp/src.f file.
FFLAGS = -g -O2 #-std=legacy
LDFLAGS =
LINKER = mpif77
LIBS = ${PETSC_LIBS} #-lquadmath

quadrule: quadrule.o
	@(cd minpack; ${MAKE} F77="$(F77)" FFLAGS="$(FFLAGS)" libminpack.a)
	@(cd ccmnf; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS)" LINKER="$(LINKER)" ccmNF.o)
	${LINKER} -o $@ $^ ccmnf/ccmNF.o ${LDFLAGS} ${MINPACK_LIB} ${GAUL_LIB} ${GSL_LIB} ${MPFR_LIB} ${LIBS}

quadrule.o: quadrule.c quadrule.h Makefile
	@if ! test -r ccmnf; then ln -s ../misc/ccmnf .; fi
	${CC} ${CFLAGS} ${CPPFLAGS} ${GAUL_CPPFLAGS} ${GSL_CPPFLAGS} ${MPFR_CPPFLAGS} ${CDEFS} -c quadrule.c

# End Makefile for quadrule
# Note: leave the line above unchanged for quadrule-dist.sh
###############################################################################

include ../Makefile.inc

dof_info.o: ../src/libphg${LIB_SUFFIX} dof_info.c
mirror.o: ../src/libphg${LIB_SUFFIX} mirror.c
grid2vtk.o: ../src/libphg${LIB_SUFFIX} grid2vtk.c

tahoe-expat: tahoe-expat.c
	$(CC) $(CFLAGS) -o $@ $^ -lexpat

FILES = .cvsignore gridview.tcl Makefile makedeps \
	?etgen2medit readvtk.c abaqus.c phgdoc.in \
	functions opendx.net opendx.sh ASCII.tcl \
	valgrind.sh valgrind.supp chown.sh \
	.cvsignore indent.sh dof_info.c quad_table \
	quadrule-plot.tex quadrule-plot.sty quadrule-plot.eps \
	quadrule.[ch] quadrule-ml quadrule-ml.res quadrule*.sh quadrule*.txt \
	quadrule*.sage quadrule-ml-bases.c unisolvence_test.c README.quadrule \
	minpack/*.f* minpack/README minpack/Makefile \
	minpack/toms-768/Src/Dp/*.f density.c \
	newton read_rule.func quad2ml newton.bc common.bc [123]Dfuncs.bc \
	[123]Dpermu.sh [123]Dpermu.bc rules-[123]d.data \
	cpp-filter.c dlamch.c f2c.h dolfin2medit \
	phg_debug tahoe2medit.c tahoe-expat.c \
	mirror.c grid2vtk.c lagrange_gen quad_symmetry.c \
	quadrule-dump-csv quadrule-from-csv quadrule-from-csv.py \
	quadrule-retrieve medit2vtk

utils.zip: $(FILES)
	@zip -9 -u -y $@ $^
