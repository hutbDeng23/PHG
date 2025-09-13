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

/* $Id: eigen.c,v 1.51 2012/12/14 04:59:57 zlb Exp $
 *
 * The eigensolver interface. */

#include "phg.h"
#include <stdlib.h>
#include <string.h>

static const char *eigen_solver_names[] = {
#if USE_ARPACK
    "arpack",
#endif
#if USE_JDBSYM
    "jdbsym",
#endif
#if USE_BLOPEX
    "blopex",
#endif
#if USE_SLEPC
    "slepc",
#endif
#if USE_TRILINOS_ANASAZI
    "trilinos",
#endif
#if USE_PRIMME
    "primme",
#endif
#if USE_GWLOBPCG
    "gwlobpcg",
#endif
    NULL
};

static EIGEN_SOLVER eigen_solver_list[] = {
#if USE_ARPACK
    phgEigenSolverARPACK,
#endif
#if USE_JDBSYM
    phgEigenSolverJDBSYM,
#endif
#if USE_BLOPEX
    phgEigenSolverBLOPEX,
#endif
#if USE_SLEPC
    phgEigenSolverSLEPc,
#endif
#if USE_TRILINOS_ANASAZI
    phgEigenSolverTrilinos,
#endif
#if USE_PRIMME
    phgEigenSolverPRIMME,
#endif
#if USE_GWLOBPCG
    phgEigenSolverGWLOBPCG,
#endif
    NULL	/* this will make the list always non-empty */
};

/* variables for command-line options */
static int eigen_solver = 0;
static INT maxit = 200;
static FLOAT eps = 1e-6;
static FLOAT tau_opt = 0.0;	/* Note: it overwrites the argument 'tau' */

static FLOAT *tmp;
static int
eigen_comp(const void *p0, const void *p1)
{
    FLOAT a = tmp[*(const int *)p0] - tmp[*(const int *)p1];
    return (a == 0. ? 0 : (a < 0. ? -1 : 1));
}

int
phgEigenSolve_(MAT *A, MAT *B, int n, int which, FLOAT tau,
	       FLOAT *evals, VEC **evecs, int *nit, void *ssd)
{
    int i, *ordering, nit0;
    VEC *vec;
    INT m;
    static BOOLEAN initialized = FALSE;

    FunctionEntry;

    if (!initialized) {
	/* register options */
	initialized = TRUE;

	phgOptionsRegisterTitle("\nEigen solver options:", "\n", "eigen");
	phgOptionsRegisterKeyword("-eigen_solver", "Eigen solver",
		eigen_solver_names, &eigen_solver);
	phgOptionsRegisterFloat("-eigen_tau", "The shift value (overwriting "
				"the argument 'tau' of phgEigenSolve)",
				&tau_opt);
	phgOptionsRegisterFloat("-eigen_tol", "Convergence tolerance", &eps);
	phgOptionsRegisterInt("-eigen_maxit",
				"Maxi. nbr of outer iterations per eigenvalue",
				&maxit);
	/* register options for OEM eigensolvers */
	for (i = 0; eigen_solver_list[i] != NULL; i++)
	    eigen_solver_list[i](NULL, NULL, 0, 0, 0., NULL, NULL, 0., 0,
								NULL, NULL);

	Return 0;
    }

    if (n <= 0)
	Return 0;

    assert(evecs != NULL);
    assert(A->rmap->nlocal == A->cmap->nlocal);
    assert(A->type != PHG_DESTROYED);
    assert(B == NULL || B->type != PHG_DESTROYED);
    assert(B == NULL || (A->rmap->nlocal == B->rmap->nlocal &&
			 A->cmap->nlocal == B->cmap->nlocal));
    assert(*evecs == NULL || (*evecs)->map->nlocal == A->cmap->nlocal);

    if (phgOptionsIfUsed("-eigen_tau"))
	tau = tau_opt;

#if 0
    if (which != EIGEN_CLOSEST)
	tau = 0.0;
#endif

    if (nit != NULL)
	*nit = 0;

    if (A->cmap->nglobal == 1) {
	/* special case: 1x1 matrices */
	n = 1;
	if ((vec = *evecs) == NULL) {
	    *evecs = vec = phgMapCreateVec(A->cmap, 1);
	}
	else {
	    vec->nvec = 1;
	    phgFree(vec->data);
	    vec->data = phgAlloc(vec->map->nlocal * sizeof(*vec->data));
	}
	if (A->rmap->nlocal == 1) {
	    *vec->data = 1.0;
	    tau = (A->rows[0].ncols == 1 ? A->rows[0].data[0] : 0.0);
	    if (B != NULL) {
		assert(B->rmap->nlocal == 1);
		assert(B->rows[0].ncols == 1 && B->rows[0].data[0] != 0.0);
		tau /= B->rows[0].data[0];
	    }
	}
	else {
	    assert(A->rmap->nlocal == 0);
	    tau = 0.0;
	}
#if USE_MPI
	MPI_Allreduce(&tau, evals, 1, PHG_MPI_FLOAT, PHG_SUM, A->rmap->comm);
#else	/* USE_MPI */
	*evals = tau;
#endif	/* USE_MPI */
    }
    else {
	if (eigen_solver_names[0] == NULL)
	    phgError(1, "%s: no available eigensolver!.\n", __func__);
	phgInfo(1, "using eigensolver %s\n", eigen_solver_names[eigen_solver]);

	/* check whether we have nlocal == 0 on some process */
	m = 0;
	for (i = 0; i < A->rmap->nprocs; i++)
	    if (A->rmap->partition[i] == A->rmap->partition[i + 1])
		m++;
	if (m > 0 && A->rmap->rank == 0)
	    phgWarning("The local matrix size is 0 on some processes!\n");

	/* call external eigensolver */
	n = eigen_solver_list[eigen_solver](A, B, n, which, tau, evals, evecs,
				eps, maxit, nit == NULL ? &nit0 : nit, ssd);
    }

    if (((*evecs)->nvec = n) > 1) {
	/* sort the eigenvalues into increasing order */
	ordering = phgAlloc(n * sizeof(*ordering));
	for (i = 0; i < n; i++)
	    ordering[i] = i;
	tmp = evals;
	qsort(ordering, n, sizeof(*ordering), eigen_comp);
	for (i = 0; i < n; i++)
	    if (ordering[i] != i)
		break;
	if (i >= n) {
	    phgFree(ordering);
	    Return n;
	}
	/* reorder data */
	vec = *evecs;
	tmp = phgAlloc(n * (m = vec->map->nlocal) * sizeof(*tmp));
	for (i = 0; i < n; i++)
	    memcpy(tmp + i * m, vec->data + ordering[i] * m, m * sizeof(*tmp));
	phgFree(vec->data);
	vec->data = tmp;
	tmp = phgAlloc(n * sizeof(*tmp));
	for (i = 0; i < n; i++)
	    tmp[i] = evals[ordering[i]];
	memcpy(evals, tmp, n * sizeof(*tmp));
	phgFree(tmp);
	phgFree(ordering);
    }

    Return n;
    return n;
}
