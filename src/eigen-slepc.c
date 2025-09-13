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

/* $Id: eigen-slepc.c,v 1.29 2014/09/19 05:21:34 zlb Exp $ */

/* Note: "petscksp.h" must be included in "phg.h" before "mpi-utils.h" */
#define NEED_PETSCKSP_H

#include "phg.h"
#include "phg/petsc-utils.h"

#if USE_SLEPC

#include <slepceps.h>
#include <stdlib.h>

static PetscScalar	*array;
static EPSWhich		whch;

static int
ecomp(const void *p0, const void *p1)
{
    PetscScalar *v0 = array + *(const int *)p0;
    PetscScalar *v1 = array + *(const int *)p1;
    int ret;

    ret = (*v0 < *v1 ? -1 : (*v0 > *v1 ? 1 : 0));

    return whch == EPS_SMALLEST_MAGNITUDE ? ret : -ret;
}

EigenSolverPrototype(phgEigenSolverSLEPc)
/* (A, B, n, which, tau, evals, evecs, eps, itmax, nit) */
{
    static BOOLEAN	initialized = FALSE;
    VEC			*vec;
    FLOAT		*data;
    int			i, j, nev, ncv, its, *ind;
    Mat			mA, mB;		/* matrices */
    Vec			evec;
    EPS			solver;		/* eigenproblem solver context */
    PetscScalar		eval;
    PetscErrorCode	ierr;

    if (!initialized) {
	/* register options */
	initialized = TRUE;
	/*phgOptionsRegisterTitle("\nSLEPc options:", "\n", "slepc");*/
	return 0;
    }

    if (n <= 0 || A->cmap->nglobal <= 0)
	return 0;

    vec = *evecs;
    MagicCheck(VEC, vec);

    if (n > A->cmap->nglobal)
	n = A->cmap->nglobal;

    if (vec == NULL) {
	*evecs = vec = phgMapCreateVec(A->cmap, n);
    }
    else {
	if (vec->nvec < n) {
	    phgFree(vec->data);
	    vec->data = phgAlloc(n * A->cmap->nlocal * sizeof(*vec->data));
	}
	vec->nvec = n;
	/* TODO: set initial vector */
    }

    assert(A->rmap->nglobal > 1);

    mA = phgPetscCreateMat(A);
    mB = phgPetscCreateMat(B);

    Unused(ierr);
    ierr = EPSCreate(A->cmap->comm, &solver);
    ierr = EPSSetOperators(solver, mA, mB);
    /* FIXME: why EPS_GHEP does not work! */
    /*ierr = EPSSetProblemType(solver, EPS_GHEP);*/
    ierr = EPSSetFromOptions(solver);
    ierr = EPSSetTolerances(solver, eps, itmax * n);
    ncv = 2 * n;
#if SLEPC_VERSION_MAJOR < 3
    ierr = EPSSetDimensions(solver, n, ncv);
#else	/* SLEPC_VERSION_MAJOR < 3 */
    ierr = EPSSetDimensions(solver, n, ncv, PETSC_DECIDE);
#endif	/* SLEPC_VERSION_MAJOR < 3 */
    whch = (which <= 0 ? EPS_SMALLEST_MAGNITUDE : EPS_LARGEST_MAGNITUDE);
    ierr = EPSSetWhichEigenpairs(solver, whch);
    ierr = EPSSolve(solver);
    ierr = EPSGetConverged(solver, &nev);
    ind = phgAlloc(nev * sizeof(*ind));
    if (nev > 0) {
	array = phgAlloc(nev * sizeof(*array));
	for (i = 0; i < nev; i++) {
	    ind[i] = i;
	    EPSGetEigenpair(solver, i, array + i, PETSC_NULL,
					PETSC_NULL, PETSC_NULL);
	}
#if 0
	EPSSortEigenvalues(nev, array, array, whch, n, ind);
#else
	qsort(ind, nev, sizeof(*ind), ecomp);
#endif
	phgFree(array);
    }
    if (nev > n)
	nev = n;
    vec->nvec = nev;
    if (nev == 0) {
	phgFree(vec->data);
	vec->data = NULL;
    }
    ierr = EPSGetIterationNumber(solver, &its);
    *nit = its;

    ierr = MatGetVecs(mA, PETSC_NULL, &evec);
    data = vec->data;
#if DEBUG
    VecGetLocalSize(evec, &j);
    assert(j == A->cmap->nlocal);
#endif	/* DEBUG */
    for (i = 0; i < nev; i++) {
	/*ierr = EPSComputeRelativeError(solver, i, &error);*/
	EPSGetEigenpair(solver, ind[i], &eval, PETSC_NULL, evec, PETSC_NULL);
	evals[i] = eval;
	VecGetArray(evec, &array);
	for (j = 0; j < A->cmap->nlocal; j++)
	    *(data++) = array[j];
	VecRestoreArray(evec, &array);
    }
    phgFree(ind);

#if (SLEPC_VERSION_MAJOR==3 && SLEPC_VERSION_MINOR<3) || SLEPC_VERSION_MAJOR<3
    ierr = EPSDestroy(solver);	/* for version < 3.3 */
#else
    ierr = EPSDestroy(&solver);	/* for version >= 3.3 */
#endif
    ierr = phgPetscVecDestroy(&evec);
    ierr = phgPetscMatDestroy(&mA);
    if (mB != NULL)
	ierr = phgPetscMatDestroy(&mB);

    return nev;
}

#else	/* USE_SLEPC */
EigenSolverVoid(phgEigenSolverSLEPc)
#endif	/* USE_SLEPC */
