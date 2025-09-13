#!/bin/bash
#
# Note: to run the program
#	cd examples
#	open.mpirun -np 1 ./poisson -solver petsc -partitioner metis

set -x

#--------------- Get list of Makefile macros
if false; then
    cat <<END >Makefile.tmp
default:
	true
END
    gmake -f Makefile.tmp -p default >macros.log 2>&1
    rm -f Makefile.tmp
    exit 0
fi

#--------------- configure
cp obj/OpenMPI-x86_64/*.o src/.
rm -f src/{refine,coarsen}.c
test -r src/libPEPCF90.so || ln -s /usr/lib64/libc.so src/libPEPCF90.so
env LDFLAGS="$LDFLAGS -L`pwd`/src" ./configure \
  --with-hypre-dir=/home_soft/soft/x86_64/lib/OpenSourceLib/hypre2.0 \
  --with-trilinos-dir=/home_soft/soft/x86_64/lib/OpenSourceLib/trilinos-9.0.3/LINUX_MPI \
  --with-mumps-dir=/home_soft/soft/x86_64/lib/OpenSourceLib/MUMPS/MUMPS_4.7.3 \
  --with-mumps-optlib="-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64"

#--------------- make
set -e
gmake lib
cd examples
gmake poisson
