/* Parallel Hierarchical Grid -- an adaptive finite element library.
 *
 * Copyright (C) 2005-2010 State Key Laboratory of Scientific and
 * Engineering Computing, Chinese Academy of Sciences. */

/* This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA */

/* $Id: petsc-utils.h,v 1.2 2017/02/28 08:02:33 zlb Exp $ */

#ifndef PHG_PETSC_UTILS_H
#define PHG_PETSC_UTILS_H

#if USE_PETSC

#include <petscksp.h>

#ifdef __cplusplus
extern "C" {
#endif

/* PETSc may define some MPI functions as macros, we undefine them here to
 * avoid some compiler warnings (e.g., in matvec.c) */
#ifdef MPI_Alltoall
# undef MPI_Alltoall
#endif
#ifdef MPI_Alltoallv
# undef MPI_Alltoallv
#endif
#ifdef MPI_Allreduce
# undef MPI_Allreduce
#endif
#ifdef MPI_Recv
# undef MPI_Recv
#endif
#ifdef MPI_Send
# undef MPI_Send
#endif

/* for compatibility with PETSc <= 2.2.0 */
#if PETSC_VERSION_MAJOR < 2 || \
    (PETSC_VERSION_MAJOR == 2 && (PETSC_VERSION_MINOR < 2 || \
        (PETSC_VERSION_MINOR == 2 && PETSC_VERSION_SUBMINOR == 0)))
typedef int PetscInt;
typedef int PetscMPIInt;
typedef int PetscErrorCode;
#endif

#if PETSC_VERSION_MAJOR > 3 || \
    (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR >= 2)
# define phgPetscKSPDestroy(a)	KSPDestroy(a)
# define phgPetscMatDestroy(a)	MatDestroy(a)
# define phgPetscVecDestroy(a)	VecDestroy(a)
typedef PetscBool PetscTruth;
#else	/* #if PETSC_VERSION_MAJOR ... */
# define phgPetscKSPDestroy(a)	KSPDestroy(*(a))
# define phgPetscMatDestroy(a)	MatDestroy(*(a))
# define phgPetscVecDestroy(a)	VecDestroy(*(a))
#endif	/* #if PETSC_VERSION_MAJOR ... */

/* for PETSc 3.3 changes */
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=3) || PETSC_VERSION_MAJOR>3
# define MatCreateMPIAIJ        MatCreateAIJ
# define MatCreateMPIBAIJ       MatCreateBAIJ
# define MatCreateMPISBAIJ      MatCreateSBAIJ
# define MatCreateMPIDense      MatCreateDense
#endif	/* PETSC_VERSION_MAJOR==3 && ... */

Mat phgPetscCreateMatAIJ(MAT *A);
Mat phgPetscCreateMatShell(MAT *A);
Mat phgPetscCreateMat(MAT *A);

#ifdef __cplusplus
}
#endif

#endif	/* USE_PETSC */

#endif	/* !defined(PHG_PETSC_UTILS_H) */
