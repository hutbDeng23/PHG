# $Id: Makefile,v 1.147 2022/09/23 06:31:13 zlb Exp $
#
# Note: 'make clean' only deletes refine.o and coarsen.o whose src files exist.
#
# To use the binary distribution of a .o file:
# 	1. delete corresponding .c file
# 	2. comment out the line 'include Makefile.dep', or delete corresponding
#	   dependency in the file Makefile.dep

.SUFFIXES: .cxx .inc .y .l

default: lib
all: lib

include ../Makefile.inc
include Makefile.dep

PHG_BUILD	= 20220923		# updated by "make dist"
CPPFLAGS	:= ${CPPFLAGS} ${BUILD_CPPFLAGS} ${USER_CPPFLAGS} \
		   -DPHG_BUILD="\"$(strip ${PHG_BUILD})\""
CFLAGS		:= ${CFLAGS} ${BUILD_CFLAGS} ${USER_CFLAGS}
CXXFLAGS	:= ${CXXFLAGS} ${BUILD_CXXFLAGS} ${USER_CXXFLAGS}

# Note: the next line allows refine.c and coarsen.c to not exist
.SECONDARY: refine.c coarsen.c

OBJS =	albert.o expr.o utils.o refine.o grid.o coarsen.o dof.o dof-utils.o \
	mapping.o geom.o quad.o quad-gauss.o option.o perf.o phg-mpi.o \
	distribute.o sturm.o mpi-utils.o mpi-2level.o medit.o gambit.o \
	elem-info.o io.o vtk.o opendx.o ensight.o \
	lagrange.o lagrange1.o mass-lumping.o specht.o \
	nedelec.o dg.o hfeb.o hierarchical-basis.o \
	solver.o solver-diagonal.o solver-phg.o solver-petsc.o solver-hypre.o \
	solver-trilinos.o solver-laspack.o solver-superlu.o \
	solver-mumps.o solver-spooles.o solver-ams.o solver-ssparse.o \
	solver-pardiso.o solver-x9amg.o solver-gather.o \
	solver-pastix.o solver-asm.o solver-smoother.o solver-hips.o \
	solver-minres.o solver-asp.o solver-2grid.o solver-hsparse.o \
	map.o matvec.o mark.o \
	eigen.o eigen-arpack.o eigen-jdbsym.o eigen-blopex.o eigen-trilinos.o \
	eigen-primme.o eigen-slepc.o eigen-gwg.o \
	surface-cut.o partition-metis.o partition-rtk.o hamilton.o \
	partition-zoltan.o hp.o partition-sfc.o partition-utils.o \
	intergrid.o periodic.o petsc-utils.o checkpoint.o oct-search.o \
	quad-interface.o quad-interface-triangle.o nonlinear-solver.o \
	quad-interface-cuboid.o halo.o quad-cache.o xfem.o xfem-utils.o

p4estdir = ../p4est
PHG_TO_P4EST_OBJS = phg-to-p4est.o xfem-p4est.o xfem-utils-p4est.o

phg-to-p4est:
	(cd ${p4estdir}; ${MAKE} ${PHG_TO_P4EST_OBJS})

lib: libphg${LIB_SUFFIX}
libs: libphg${LIB_SUFFIX}
libphg.a: ${OBJS} phg-to-p4est
	${AR} -ruv $@ ${OBJS} \
		`for f in ${PHG_TO_P4EST_OBJS}; do echo ${p4estdir}/$$f; done`
	${RANLIB} $@
libphg.dll: ${OBJS} phg-to-p4est
	@make libphg.so
libphg.so: ${OBJS} phg-to-p4est
	@echo 1>&2 "libphg${LIB_SUFFIX}: refreshing $?"
	@if ! ${BUILD_SHARED} ${OBJS} \
		`for f in ${PHG_TO_P4EST_OBJS}; do echo ${p4estdir}/$$f; done` \
		${BUILD_SHARED_LIBS} -o libphg${LIB_SUFFIX}; then \
    echo 1>&2 "*************************************************************";\
    echo 1>&2 "* Error: building libphg.so failed (see errors above).      *";\
    echo 1>&2 "*                                                           *";\
    echo 1>&2 "* If the error is caused by an incompatible static library, *";\
    echo 1>&2 "* you may consider reconfigure PHG using either the option  *";\
    echo 1>&2 "* '--disable-link' to skip linking libphg.so, or the option *";\
    echo 1>&2 "* '--disable-shared' to create a static library libphg.a    *";\
    echo 1>&2 "*************************************************************";\
	    false; \
	fi
	@#chcon -t textrel_shlib_t $@

clean:
	@sh -c 'for f in *.o; do					\
		    if test "$$f" != "refine.o" -a "$$f" != "coarsen.o";\
		    then						\
			/bin/rm -f $$f;					\
		    elif test						\
			-r `echo $$f | sed -e "s/\.o/.c/g"` -o		\
			-r `echo $$f | sed -e "s/\.o/.cxx/g"` -o	\
			-r `echo $$f | sed -e "s/\.o/.cpp/g"` -o	\
			-r `echo $$f | sed -e "s/\.o/.f/g"` -o		\
			-r `echo $$f | sed -e "s/\.o/.F/g"` -o		\
			-r `echo $$f | sed -e "s/\.o/.f90/g"`;		\
		    then						\
			/bin/rm -f $$f;					\
		    fi;							\
	done'

distclean: clean
	-/bin/rm -f libphg.a libphg.so

mass-lumping.inc: ../utils/mass-lumping-gen
	@echo 1>&2 "Regenerating mass-limping.inc ..."
	bash ../utils/mass-lumping-gen >$@
expression.o: expression.c expression.inc
expression0.o: expression0.c expression0.inc

tetrahedralize: elem-info.c ../include/elem-info.h
	${CC} ${CPPFLAGS} -DTEST=1 -o $@ elem-info.c
	rm -f $@.o

.y.c:
	if yacc $*.y; then \
	    mv -f y.tab.c $*.c; \
	else \
	    touch $*.c; \
	fi

.l.inc:
	if lex -P _phg_yy $*.l; then \
	    mv -f lex._phg_yy.c $*.inc; \
	else \
	    touch $*.inc; \
	fi

# Hacky: perform partial linking to speed up later linkings
#.cxx.o:
#	${CXX} ${CXXFLAGS} ${CPPFLAGS} -c $*.cxx
#	${CXX} -r -o tmp.o $*.o \
#		-L/usr/local/trilinos/lib -lteuchos -ltriutils -lepetra
#	mv -f tmp.o $*.o

.PHONY: default all clean distclean phg-to-p4est lib libs

FILES = Makefile .cvsignore Makefile.dep \
	expr.c expression*.[ylc] expression.inc \
	refine.c grid.c phg-mpi.c distribute.c utils.c \
	mpi-utils.c mpi-2level.c geom.c quad.c dof.c dof-utils.c mapping.c \
	coarsen.c quad-gauss.c sturm.c \
	medit.c gambit.c albert.c io.c vtk.c opendx.c ensight.c \
	lagrange.c nedelec.c dg.c option.c mark.c \
	solver*.c solver-*.cxx eigen*.c eigen-*.cxx \
	hierarchical-basis.c hfeb.c perf.c hp.c \
	elem-info.c map.c matvec.c surface-cut.c \
	partition-metis.c hamilton.c partition-rtk.c \
	partition-zoltan.c partition-utils.c partition-sfc.c \
	intergrid.c periodic.c lagrange1.c mass-lumping.c petsc-utils.c \
	checkpoint.c specht.c oct-search.c \
	quad-interface.c quad-interface-triangle.c nonlinear-solver.c \
	quad-interface-cuboid.cxx halo.c quad-cache.c xfem.c

src.zip: $(FILES)
	@zip -9 -u -y $@ $^
