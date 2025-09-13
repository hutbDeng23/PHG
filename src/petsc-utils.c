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

/* $Id: petsc-utils.c,v 1.10 2015/06/19 09:00:26 zlb Exp $ */

/* Note: "petscksp.h" must be included in "phg.h" before "mpi-utils.h" */
#define NEED_PETSCKSP_H

#include "phg.h"

#if USE_PETSC

/* for PETSc 3.3 changes */
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=3) || PETSC_VERSION_MAJOR>3
# define MatCreateMPIAIJ	MatCreateAIJ
# define MatCreateMPIBAIJ	MatCreateBAIJ
# define MatCreateMPISBAIJ	MatCreateSBAIJ
# define MatCreateMPIDense	MatCreateDense
#endif

Mat
phgPetscCreateMatAIJ(MAT *A)
/* Setup PETSc AIJ matrix */
/* By copying data from MAT A to PETSc Mat a */
{
    Mat a;
    const MAT_ROW *row;
    PetscInt prealloc = PETSC_DECIDE, *d_nnz, *o_nnz;
    PetscInt gi, i, j, maxcols, *cols;
    PetscScalar *data;

    assert(A != NULL);
    assert(A->type != PHG_DESTROYED);

    if (!A->assembled)
	phgMatAssemble(A);

    d_nnz = phgAlloc(A->rmap->nlocal * sizeof(*d_nnz));
    o_nnz = phgAlloc(A->rmap->nlocal * sizeof(*o_nnz));
    maxcols = 0;
    for (i = 0; i < A->rmap->nlocal; i++) {
	if (sizeof(INT) == sizeof(PetscInt)) {
	    phgMatGetNnzInRow(A, i, (INT *)d_nnz + i, (INT *)o_nnz + i);
	}
	else {
	    INT dn, on;
	    phgMatGetNnzInRow(A, i, &dn, &on);
	    d_nnz[i] = dn;
	    o_nnz[i] = on;
	}
	if (maxcols < d_nnz[i] + o_nnz[i])
	    maxcols = d_nnz[i] + o_nnz[i];
    }
    if (A->rmap->nlocal == 0 || A->cmap->nlocal == 0)
	phgInfo(-1, "%s: WARNING: local size is 0 (PETSc may hang).\n",
			__func__);
    MatCreateMPIAIJ(A->rmap->comm, A->rmap->nlocal, A->cmap->nlocal,
		    A->rmap->nglobal, A->cmap->nglobal,
		    prealloc, d_nnz, prealloc, o_nnz, &a);
    /*MatSetFromOptions(a);*/
    phgFree(d_nnz);
    phgFree(o_nnz);

    gi = A->rmap->partition[A->rmap->rank];
    cols = phgAlloc(maxcols * sizeof(*cols));
    data = phgAlloc(maxcols * sizeof(*data));
    for (i = 0; i < A->rmap->nlocal; i++, gi++) {
	row = phgMatGetRow(A, i);
	assert(row->ncols <= maxcols);
	if (row->ncols > 0) {
	    for (j = 0; j < row->ncols; j++) {
		cols[j] = row->cols[j];
		data[j] = (PetscScalar)row->data[j];
	    }
	    MatSetValues(a, 1, &gi, row->ncols, cols, data, INSERT_VALUES);
	}
    }
    phgFree(cols);
    phgFree(data);

    MatAssemblyBegin(a, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(a, MAT_FINAL_ASSEMBLY);

    return a;
}

static int
to_petsc(VEC *x, Vec v)
/* copy PHG VEC 'x' to Petsc Vec 'v' */
{
    PetscScalar *vec;
    PetscInt size;
    PetscErrorCode code;
    INT i;

    if (v == PETSC_NULL)
	phgError(1, "PETSc CopySolution: invalid PETSc vector.\n");

    code = VecGetLocalSize(v, &size); CHKERRQ(code);
    code = VecGetArray(v, &vec); CHKERRQ(code);
    for (i = 0; i < size; i++)
	vec[i] = x->data[i];
    code = VecRestoreArray(v, &vec); CHKERRQ(code);

    return 0;
}

static int
from_petsc(VEC *x, Vec v)
/* copy Petsc Vec 'v' to PHG VEC 'x' */
{
    PetscScalar *vec;
    PetscInt size;
    PetscErrorCode code;
    INT i;

    if (v == PETSC_NULL)
	phgError(1, "PETSc CopySolution: invalid PETSc vector.\n");

    code = VecGetLocalSize(v, &size); CHKERRQ(code);
    code = VecGetArray(v, &vec); CHKERRQ(code);
    for (i = 0; i < size; i++)
	x->data[i] = vec[i];
    code = VecRestoreArray(v, &vec); CHKERRQ(code);

    return 0;
}

static VEC *_wx = NULL;
static VEC *_wy = NULL;

static int
matop_mult(Mat mat, Vec x, Vec y)
{
    MAT *phg_mat;

    /* Note: ctx points to PHG MAT struct */
    MatShellGetContext(mat, (void *)&phg_mat);
    from_petsc(_wx, x);
    phgMatVec(MAT_OP_N, 1.0, phg_mat, _wx, 0.0, &_wy);
    to_petsc(_wy, y);

    return 0;
}

static int
matop_mult_transpose(Mat mat, Vec x, Vec y)
{
    MAT *phg_mat;

    /* Note: ctx points to PHG MAT struct */
    MatShellGetContext(mat, (void *)&phg_mat);
    from_petsc(_wx, x);
    phgMatVec(MAT_OP_T, 1.0, phg_mat, _wx, 0.0, &_wy);
    to_petsc(_wy, y);
    
    return 0;
}

static int
matop_get_diagonal(Mat mat, Vec x)
{
    MAT *phg_mat;
    PetscScalar *vec;
    INT i;

    /* Note: ctx points to PHG MAT struct */
    MatShellGetContext(mat, (void *)&phg_mat);
    if (phg_mat->diag == NULL)
	phgMatSetupDiagonal(phg_mat);
    VecGetArray(x, &vec);
    for (i = 0; i < phg_mat->rmap->nlocal; i++)
	vec[i] = phg_mat->diag[i];
    VecRestoreArray(x, &vec);

    return 0;
}

Mat
phgPetscCreateMatShell(MAT *A)
/* Setup PHG's MatrixFreeMatrix for PETSc's Shell matrix */
{
    Mat a;

    assert(A != NULL);
    assert(A->type != PHG_DESTROYED);

    MatCreateShell(A->rmap->comm, A->rmap->nlocal, A->cmap->nlocal,
		   A->rmap->nglobal, A->cmap->nglobal, (void *)A, &a);
    MatShellSetOperation(a, MATOP_MULT, (void(*)(void))&matop_mult);
    MatShellSetOperation(a, MATOP_MULT_TRANSPOSE,
				(void(*)(void))&matop_mult_transpose);
    MatShellSetOperation(a, MATOP_GET_DIAGONAL,
				(void(*)(void))&matop_get_diagonal);
    MatSetFromOptions(a);

    return a;
}

Mat
phgPetscCreateMat(MAT *A)
{
#if 0
    /* This will force using shell matrix for testing. To test:
     *	eigen -eigen_solver slepc -oem_options '-st_pc_type jacobi' */
    if (TRUE)
#else	/* 0 */
    if (A->type == PHG_MATRIX_FREE && A->blocks == NULL)
#endif	/* 0 */
	return phgPetscCreateMatShell(A);
    else
	return phgPetscCreateMatAIJ(A);
}

#endif	/* USE_PETSC */
