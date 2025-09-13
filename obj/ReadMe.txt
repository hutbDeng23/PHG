The object files 'refine.o' and 'coarsen.o' in this directory support parallel
mesh refinement and coarsening.

Usage:	Copy the files 'refine.o' and 'coarsen.o' from a subdirectory which
	matches your OS type, MPI and data types to the 'src/' directory
	(overwriting existing ones), and (optionally) remove the files
	'refine.c' and 'coarsen.c' from the 'src/' directory.

	For example, if you are using MPICH2 on an x86_64 machine, then
	you can use the following commands to build PHG with full parallel
	refinement/coarsening support (suppose 'mpicc', 'mpicxx', etc., are
	in the PATH. 'double_int' in the pathname corresponds to data types,
	which means use 'double' and 'int' for 'FLOAT' and 'INT' respectively):
		tar xjpvf phg-x.x.x.tar.bz2
		cd phg-x.x.x
		cp obj/mpich-x86_64/double_int/*.o src/.
		rm -f src/{refine,coarsen}.c src/libphg.*
		./configure
		make

Note:	Object files for MPICH can also be used with other MPI implementations
	which are based on or derived from MPICH, like MVAPICH and Intel MPI.

-------------------------------------------------------------------------------

List of subdirectories:

	Subdirectory name		Target system
	-----------------		-------------
	mpich-i386/			MPICH, Linux, i386
	mpich-x86_64/			MPICH, Linux, x86_64
	openmpi-i386/			OpenMPI, Linux, i386
	openmpi-x86_64/			OpenMPI, Linux, x86_64

	intelmpi-mic/			Intel MPI, Intel Xeon Phi (MIC)

	macports-mpich-x86_64/		MPICH, MacOS, x86_64
	macports-openmpi-x86_64/	OpenMPI, MacOS, x86_64

	sgimpt-x86_64/			MPT, SGI Altix UV, x86_64
