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

/* $Id: eigen-primme.c,v 1.27 2011/04/14 06:34:17 zlb Exp $ */

#include "phg.h"

#if USE_PRIMME

#include <primme.h>

static MAT *ma, *mb;
static VEC *va, *vb;
static FLOAT *diag;

/* Note: the current version of PRIMME does not allow the case mb != NULL.
 * A simple workaround consists of first decompose B = L L^T, then compute
 *	L^(-1) A L^(-T) y = lambda y
 * with y = L^T x */

static void
mat_vec(MAT *a, FLOAT *x, FLOAT *y, int nvec)
/* computes y := Ax */
{
    va->data = x;
    vb->data = y;
    va->nvec = vb->nvec = nvec;
    phgMatVec(MAT_OP_N, 1.0, a, va, 0.0, &vb);

    va->data = NULL;
    vb->data = NULL;
}

static void
mat_vec_a(void *x, void *y, int *blockSize, struct primme_params *primme)
/* computes y := Ax */
{
    mat_vec(ma, x, y, *blockSize);
}

static void
mat_vec_b(void *x, void *y, int *blockSize, struct primme_params *primme)
/* computes y := Bx */
{
    mat_vec(mb, x, y, *blockSize);
}

/* diagonal preconditioner */
static void
precon_diag(void *x0, void *y0, int *blockSize, struct primme_params *primme)
{ 
    FLOAT *x = x0, *y = y0;
    INT i;

    assert(*blockSize == 1);
    for (i = 0; i < ma->rmap->nlocal; i++)
	y[i] = diag[i]*x[i];
}
#if USE_MPI
static void
globalSumDouble(void *sbuf, void *rbuf, int *count,
				struct primme_params *primme)
{
    MPI_Allreduce(sbuf, rbuf, *count, PHG_MPI_FLOAT, PHG_SUM,
					*(MPI_Comm *)primme->commInfo);
}
#endif	/* USE_MPI */

EigenSolverPrototype(phgEigenSolverPRIMME)
/* (A, B, n, which, tau, evals, evecs, eps, itmax, nit) */
{
    static BOOLEAN initialized = FALSE;
    INT i, j, n0;
    MAT_ROW *row;
    VEC *vec;
    int ret;
    FLOAT *rnorms;
    primme_params primme;

    if (!initialized) {
	/* register options */
	initialized = TRUE;
	/*phgOptionsRegisterTitle("\nPRIMME options:", "\n", "primme");*/
	return 0;
    }

    assert(FT_PHG == FT_DOUBLE);
    assert(A->type == PHG_UNPACKED);

    if (B != NULL)
	phgError(1, "%s can't compute generalized eigenvalues\n", __func__);

    if (n <= 0 || A->cmap->nglobal <= 0)
	return 0;

    n0 = A->cmap->partition[A->cmap->rank];
    vec = *evecs;
    MagicCheck(VEC, vec);

    rnorms = phgAlloc(n * sizeof(*rnorms));

    primme_initialize(&primme);
#if 0
    primme_set_method(DYNAMIC, &primme);	/* divergent! */
#else
    primme_set_method(/*JDQR*/JDQMR_ETol, &primme);
#endif

    primme.maxOuterIterations = itmax * n;
    primme.correctionParams.convTest = primme_full_LTolerance;
    /*primme.correctionParams.convTest = primme_decreasing_LTolerance;*/
    primme.eps = eps;
    primme.n = A->cmap->nglobal;
#if USE_MPI
    /* parallel PRIMME */
    primme.numProcs = A->cmap->nprocs;
    primme.procID = A->cmap->rank;
    primme.commInfo = &A->cmap->comm;
    primme.nLocal = A->cmap->nlocal;
    primme.globalSumDouble = globalSumDouble;
#endif	/* USE_MPI */

    primme.numEvals = n;
    primme.matrixMatvec = mat_vec_a;

    primme.massMatrixMatvec = mat_vec_b;

    primme.applyPreconditioner = precon_diag;
    primme.correctionParams.precondition = 0;	/* 0 ==> don't do precon. */

    primme.initSize = n;
    primme.printLevel = phgVerbosity;

    /* primme.basisSize = ;*/
    /* primme.blockSize = ;*/

    ma = A;
    mb = B;
    va = phgMapCreateVec(A->cmap, 1);
    vb = phgMapCreateVec(A->cmap, 1);
    phgFree(va->data);
    va->data = NULL;
    phgFree(vb->data);
    vb->data = NULL;

    diag = phgAlloc(A->rmap->nlocal * sizeof(FLOAT));
    for (i = 0, row = A->rows; i < A->rmap->nlocal; i++, row++) {
	for (j = 0; j < row->ncols; j++)
	    if (row->cols[j] == i)
		break;
	if (j >= row->ncols || row->data[j] == 0.0)
	    phgError(1, "diagonal entry of row %d is zero!\n", i + n0);
	diag[i] = 1.0 / row->data[j];
    } 

    if (vec == NULL) {
	*evecs = vec = phgMapCreateVec(A->cmap, n);
	primme.initSize = 0;
    }
    else if (vec->nvec != n) {
	vec->nvec = n;
	phgFree(vec->data);
	vec->data = phgAlloc(n * A->cmap->nlocal * sizeof(*vec->data));
    }

    assert(A->rmap->nglobal > 1);

    ret = dprimme(evals, vec->data, rnorms, &primme);
    if (ret != 0)
	phgError(1, "PRIMME ret value = %d\n", ret);

    *nit = primme.stats.numOuterIterations;

    /* FIXME: free primme */
    primme_Free(&primme);
    phgFree(rnorms);
    phgFree(diag);
    phgVecDestroy(&va);
    phgVecDestroy(&vb);

    return primme.initSize;
}

#else	/* USE_PRIMME */
EigenSolverVoid(phgEigenSolverPRIMME)
#endif	/* USE_PRIMME */
