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

/* $Id: solver.c,v 1.397 2022/09/23 09:06:49 zlb Exp $ */

#define NEED_GetVariableArgs

#include "phg.h"
#include "phg/blas-lapack.h"

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdarg.h>
#include <math.h>
#include <limits.h>	/* INT_MAX */

#if !defined(PHG_TO_P4EST)	/* avoid compile/link following funcs twice */

/* the workaround when the local size of the matrix is 0 on some processes.
 * Typical test case:
 *      mpirun -np 4 poisson -fem_type=P1 -mesh_file=tetra.dat -pre_refines=2
 */

/* list of solvers */
#define NSolvers	(sizeof(phgSolverList_)/sizeof(phgSolverList_[0]) - 1)

/* Solver names */
const char *phgSolverNames_[] = {
    "none",		/* dummy entry used by other interfaces (BLOPEX) */
    "gmres",
#if USE_MINRES
    "minres",
#endif
    "pcg",
    "bicgstab",
    "preonly",
    "asp",
    "2grid",
    "ams",
    "gather",
    "asm",
    "smoother",
#if USE_HYPRE
    "hypre",
#endif
#if USE_X9AMG
    "x9amg",
#endif
#if USE_PETSC
    "petsc",
#endif
#if USE_TRILINOS
    "trilinos",
#endif
#if USE_MUMPS
    "mumps",
#endif
#if USE_PASTIX
    "pastix",
#endif
#if USE_HIPS
    "hips",
#endif
#if USE_PARDISO
    "pardiso",
#endif
#if USE_SUPERLU
    "superlu",
#endif
#if USE_SPOOLES
    "spooles",
#endif
#if USE_SSPARSE
    "ssparse",
#endif
#if USE_LASPACK
    "laspack",
#endif
#if USE_HSPARSE
    "hsparse",
#endif
    "diagonal",
    NULL
};

/* Solvers list */
OEM_SOLVER *phgSolverList_[] = {
    NULL,		/* dummy entry reserved for other interfaces */
    SOLVER_GMRES,
#if USE_MINRES
    SOLVER_MINRES,
#endif
    SOLVER_PCG,
    SOLVER_BiCGstab,
    SOLVER_PreOnly,
    SOLVER_ASP,
    SOLVER_2Grid,
    SOLVER_AMS,
    SOLVER_GATHER,
    SOLVER_ASM,
    SOLVER_SMOOTHER,
#if USE_HYPRE
    SOLVER_HYPRE,
#endif
#if USE_X9AMG
    SOLVER_X9AMG,
#endif
#if USE_PETSC
    SOLVER_PETSC,
#endif
#if USE_TRILINOS
    SOLVER_TRILINOS,
#endif
#if USE_MUMPS
    SOLVER_MUMPS,
#endif
#if USE_PASTIX
    SOLVER_PASTIX,
#endif
#if USE_HIPS
    SOLVER_HIPS,
#endif
#if USE_PARDISO
    SOLVER_PARDISO,
#endif
#if USE_SUPERLU
    SOLVER_SUPERLU,
#endif
#if USE_SPOOLES
    SOLVER_SPOOLES,
#endif
#if USE_SSPARSE
    SOLVER_SSPARSE,
#endif
#if USE_LASPACK
    SOLVER_LASPACK,
#endif
#if USE_HSPARSE
    SOLVER_HSPARSE,
#endif
    SOLVER_DIAGONAL
};
static int default_solver = 0;

/* Parameters */
static INT maxit = 10000;
static BOOLEAN warn_maxit = TRUE;
static FLOAT rtol = 1e-10;
static FLOAT atol_ = 0.;	/* avoid conflicts with stdlib.h */
static FLOAT btol = 1e-10;
static BOOLEAN monitor = FALSE;	/* solver->monitor */
static BOOLEAN cond = FALSE;	/* whether estimate condition number */
static int scaling = 0;	/* solver->scaling */
static const char *scaling_args[] = {"off", "left", "right", "sym", NULL};
static int symmetry = M_UNSYM;	/* solver->symmetry */
static const char *symmetry_args[] =
	{"unsym", "ssym", "sym", "spd", /*"auto",*/ NULL};

static int matrix_free = ALLOW;	/* enum {ALLOW, DISALLOW, ENFORCE} */
static const char *matrix_free_args[] = {"allow", "disallow", "enforce", NULL};

int
phgSolverRegisterOptions(void)
{
    int i;

    phgOptionsRegisterTitle("\nLinear solver options:", "\n", "solver");
    phgOptionsRegisterKeyword("-solver", "Default solver",
			      phgSolverNames, &default_solver);
    phgOptionsRegisterInt("-solver_maxit",
			  "Default maxi. number of iterations", &maxit);
    phgOptionsRegisterNoArg("-solver_cond",
			    "Estimate condition number (solver->cond)", &cond);
    phgOptionsRegisterNoArg("-solver_monitor",
			    "Print convergence history", &monitor);
    phgOptionsRegisterNoArg("-solver_warn_maxit",
			    "Warn if maxi. number of iterations attained",
				&warn_maxit);
    phgOptionsRegisterFloat("-solver_rtol", "Relative tolerance "
			    "(||residual|| <= ||residual0|| * rtol)", &rtol);
    phgOptionsRegisterFloat("-solver_btol", "B tolerance "
			    "(||residual|| <= ||RHS|| * btol)", &btol);
    phgOptionsRegisterFloat("-solver_atol", "Absolute tolerance "
			    "(||residual|| <= atol)", &atol_);
    phgOptionsRegisterKeyword("-solver_symmetry", "Symmetry of the matrix",
			symmetry_args, &symmetry);
    phgOptionsRegisterKeyword("-solver_scaling", "Perform diagonal scaling",
			scaling_args, &scaling);
    phgOptionsRegisterKeyword("-solver_matrix_free",
			"Use of matrix-free form in solver when supported",
			matrix_free_args, &matrix_free);

    /* Register eigensolver options */
    phgEigenSolve(NULL, NULL, 0, 0, 0.0, NULL, NULL, NULL);

    for (i = 0; i < NSolvers; i++)
	if (phgSolverList[i] != NULL &&
	    phgSolverList[i]->RegisterOptions != NULL)
	    phgSolverList[i]->RegisterOptions();

    return 0;
}

static int
send_matrix_to_oem_solver(SOLVER *solver)
/* sends data to OEM solver */
{
    MAT *mat = solver->mat;
    INT i, j;
    const MAT_ROW *row;
    FLOAT total_cols;

    FunctionEntry;

    if (solver->sent_to_oem)
	Return 0;

    if (!solver->oem_created) {
	solver->oem_created = TRUE;
	solver->oem_solver->Create(solver);
    }

    /* check whether the solver has AddMatrixEntries function */
    if (solver->oem_solver->AddMatrixEntries == NULL)
	Return 0;

    if (solver->oem_solver->allow_matrix_free)
	if ((solver->matrix_free == ALLOW && mat->type == PHG_MATRIX_FREE)
	    || solver->matrix_free == ENFORCE)
	    Return 0;

    assert(mat->type != PHG_DESTROYED);

    phgInfo(1, "sending matrix to OEM solver...\n");

    /* send matrix to OEM solver and free PHG's matrix */
    j = mat->rmap->partition[mat->rmap->rank];
    total_cols = 0.0;
    for (i = 0; i < mat->rmap->nlocal; i++, j++) {
	row = phgMatGetRow(mat, i);
	if (row->ncols <= 0)
	    phgError(1, "%s: matrix row %d is empty!\n", __func__, i);
	total_cols += row->ncols;
	solver->oem_solver->AddMatrixEntries(solver, 1, &j, row->ncols,
					     row->cols, row->data);

	if (mat->refcount == 0 && mat->rows != NULL) {
	    /* free PHG matrix row */
	    phgFree(mat->rows[i].cols);
	    phgFree(mat->rows[i].data);
	    mat->rows[i].cols = NULL;
	    mat->rows[i].data = NULL;
	    mat->rows[i].ncols = mat->rows[i].alloc = 0;
	}
    }

    phgInfo(1, "average cols = %lg\n", solver->mat->rmap->nlocal == 0 ?
		0 : ((double)total_cols) / solver->mat->rmap->nlocal);

    if (mat->refcount == 0)
	    phgMatFreeMatrix(mat);

    solver->sent_to_oem = TRUE;

    Return 0;
    return 0;		/* to avoid possible gcc warnings */
}

BOOLEAN
phgSolverDenseLU(int n, FLOAT *a0, int pvt[])
/* LU factorization (partial pivoting). TODO: call LAPACK for large n */
{
#if USE_LAPACK && (FT_PHG == FT_DOUBLE || FT_PHG == FT_FLOAT) \
	&& SIZEOF_BLAS_INT == SIZEOF_INT
    /* use LAPACK XGETRF. Note the matrix here is in the C order, thus
     * TRANS='T' will be used when calling XGETRS in phgSolverDenseSV */
    int info;

#if FT_PHG == FT_DOUBLE
    F77_FUNC(dgetrf,DGETRF)(&n, &n, a0, &n, pvt, &info);
#elif FT_PHG == FT_FLOAT
    F77_FUNC(sgetrf,SGETRF)(&n, &n, a0, &n, pvt, &info);
#endif
    return info == 0;

#else	/* USE_LAPACK */
    int i, j, k;
    FLOAT d, r;
#define a(i,j)	(a0[(i) * n + (j)])

    for (i = 0; i < n; i++) {
	d = Fabs(a(i,i));
	k = i;
	for (j = i + 1; j < n; j++) {
	    if ((r = Fabs(a(j,i))) > d) {
		d = r;
		k = j;
	    }
	}
	if (d == 0.0) {
#if 0 * DEBUG
	    printf("Dense solver failed: i = %d\n", i);
	    for (j = 0; j < n; j++) {
		for (k = 0; k < n; k++)
		    printf(" %9.2le", (k >= i || k >= j) ? (double)a(j,k) : 0.);
	    }
#endif
	    return FALSE;
	}
	pvt[i] = k;
	if (k != i) {
	    /* exchange row i and row k */
	    for (j = i; j < n; j++) {
		d = a(i,j);
		a(i,j) = a(k,j);
		a(k,j) = d;
	    }
	}
	if ((d = a(i,i)) != 1.0) {
	    a(i,i) = d = 1.0 / d;
	    for (j = i + 1; j < n; j++)
		a(i,j) *= d;
	}
	for (j = i + 1; j < n; j++) {
	    if ((d = a(j,i)) == 0.0)
		continue;
	    for (k = i + 1; k < n; k++)
		a(j,k) -= d * a(i,k);
	}
    }
#undef a

    return TRUE;
#endif	/* USE_LAPACK */
}

void
phgSolverDenseSV(int n, FLOAT *a0, int pvt[], int m, FLOAT *b0)
/* solves L*U*X = B. TODO: call LAPACK for large n */
{
#if USE_LAPACK && (FT_PHG == FT_DOUBLE || FT_PHG == FT_FLOAT) \
	&& SIZEOF_BLAS_INT == SIZEOF_INT
    int i, j, info;
    FLOAT *b = b0;

    if (m > 1) {
	/* Need to set up B in Fortran order */
	b = phgAlloc(n * m * sizeof(*b));
	for (i = 0; i < n; i++)
	    for (j = 0; j < m; j++)
		b[j * n + i] = b0[i * m + j];
    }
#if FT_PHG == FT_DOUBLE
    F77_FUNC(dgetrs,DGETRS)("Trans", &n, &m, a0, &n, pvt, b, &n, &info, 1);
#elif FT_PHG == FT_FLOAT
    F77_FUNC(sgetrs,SGETRS)("Trans", &n, &m, a0, &n, pvt, b, &n, &info, 1);
#endif
    if (m > 1) {
	/* Restore B from Fortran order */
	for (i = 0; i < n; i++)
	    for (j = 0; j < m; j++)
		b0[i * m + j] = b[j * n + i];
	phgFree(b);
    }

#else	/* USE_LAPACK */
    int i, j, k;
    FLOAT d;
#define a(i,j)	(a0[(i) * n + (j)])
#define b(i,j)	(b0[(i) * m + (j)])

    for (i = 0; i < n; i++) {
	k = pvt[i];
	if (k != i) {
	    /* exchange row i with row k */
	    for (j = 0; j < m; j++) {
		d = b(i,j);
		b(i,j) = b(k,j);
		b(k,j) = d;
	    }
	}
	if ((d = a(i,i)) != 1.0) {
	    for (j = 0; j < m; j++)
		b(i,j) *= d;
	}
	for (j = i + 1; j < n; j++) {
	    if ((d = a(j,i)) == 0.0)
		continue;
	    for (k = 0; k < m; k++)
		b(j,k) -= d * b(i,k);
	}
    }

    for (i = n - 2; i >= 0; i--)
	for (j = i + 1; j < n; j++) {
	    d = a(i,j);
	    for (k = 0; k < m; k++)
		b(i,k) -= d * b(j,k);
	}
#undef a
#undef b
#endif	/* USE_LAPACK */

    return;
}

BOOLEAN
phgSolverDenseSolver(int n, int m, FLOAT *a0, FLOAT *b0)
/* solves AX = B, returns TRUE if successful and FALSE if A is singular.
 * TODO: call LAPACK for large n */
{
#if 1
    int i, j, k;
    FLOAT d, r;
#define a(i,j)	(a0[(i) * n + (j)])
#define b(i,j)	(b0[(i) * m + (j)])

    for (i = 0; i < n; i++) {
	d = Fabs(a(i,i));
	k = i;
	for (j = i + 1; j < n; j++) {
	    if ((r = Fabs(a(j,i))) > d) {
		d = r;
		k = j;
	    }
	}
	if (d == 0.0) {
#if 0 * DEBUG
	    printf("Dense solver failed: i = %d\n", i);
	    for (j = 0; j < n; j++) {
		for (k = 0; k < n; k++)
		    printf(" %9.2le", (k >= i || k >= j) ? (double)a(j,k) : 0.);
		for (k = 0; k < m; k++)
		    printf(" = %9.2le\n", (double)b(j,k));
	    }
#endif
	    return FALSE;
	}
	if (k != i) {
	    /* exchange row i with row k */
	    for (j = i; j < n; j++) {
		d = a(i,j);
		a(i,j) = a(k,j);
		a(k,j) = d;
	    }
	    for (j = 0; j < m; j++) {
		d = b(i,j);
		b(i,j) = b(k,j);
		b(k,j) = d;
	    }
	}
	if ((d = a(i,i)) != 1.0) {
	    d = 1.0 / d;
	    for (j = i + 1; j < n; j++)
		a(i,j) *= d;
	    for (j = 0; j < m; j++)
		b(i,j) *= d;
	}
	for (j = i + 1; j < n; j++) {
	    if ((d = a(j,i)) == 0.0)
		continue;
	    for (k = i + 1; k < n; k++)
		a(j,k) -= d * a(i,k);
	    for (k = 0; k < m; k++)
		b(j,k) -= d * b(i,k);
	}
    }

    for (i = n - 2; i >= 0; i--)
	for (j = i + 1; j < n; j++) {
	    d = a(i,j);
	    for (k = 0; k < m; k++)
		b(i,k) -= d * b(j,k);
	}
#undef a
#undef b
#else
    int *pvt = phgAlloc(n * sizeof(*pvt));

    if (phgSolverDenseLU(n, a0, pvt) == FALSE)
	return FALSE;
    phgSolverDenseSV(n, a0, pvt, m, b0);
    phgFree(pvt);
#endif

    return TRUE;
}

static void
update_rhs(MAT *smat, VEC *srhs, MAT *dmat, VEC *drhs, FLOAT coeff,
	   BOOLEAN diag)
/* This function computes:
 *	drhs -= coeff * B(dmat) * D(smat)^(-1) * srhs
 * where:
 *	A(mat) = mat with all bdry columns set to 0,
 *	B(mat) = the matrix composed of the bdry columns of mat,
 *	D(mat) = diag matrix with non-bdry part = I and bdry part = diag(mat)
 *
 * Note: if mat->handle_bdry_eqns == TRUE, then after MAT assembly, only the
 * diagonal entries in the bdry rows can be nonzero, and A is stored in
 * mat->rows, B+D is stored in mat->bdry_rows (note both mat->rows and
 * mat->bdry_rows contain D). Since mat->rows might have been destroyed when
 * sending data to OEM solvers, we get D from mat->bdry_rows.
 *
 * TODO: replace bdry_rows with a MAT object, and use mat-vec operations. */
{
    MAT *mat;
    VEC *rhs;
    MAT_ROW *row;
    INT i, j, k;
#if USE_MPI
    INT O2Gsize, n, *rlist, *slist;
    FLOAT *rbuff, *sbuff;
    BYTE *flags;
    int rank, rsize, *sdsps, *scnts, *rdsps, *rcnts;
#endif	/* USE_MPI */

    if (smat == NULL || !smat->handle_bdry_eqns)
	return;

    assert(srhs->nvec == drhs->nvec);

    if (smat->blocks != NULL) {
	FLOAT *sdata, *ddata, c0, c;
	INT p = smat->blocks->p, q = smat->blocks->q;
	assert(smat->type == PHG_MATRIX_FREE || smat->type == PHG_DESTROYED);
	assert(dmat->blocks != NULL && p == q);
	if (p == 1) {
	    assert(dmat->blocks->p == 1 && dmat->blocks->q == 1);
	    if ((c = dmat->blocks->coeff[0]) == 0.)
		return;
	    c0 = smat->blocks->coeff[0];
	    assert(c0 != 0.);
	    update_rhs(smat->blocks->pmat[0], srhs,
		       dmat->blocks->pmat[0], drhs, c / c0, TRUE);
	    return;
	}
	sdata = srhs->data;
	for (j = 0; j < q; j++) {
	    VEC *rhs0;
	    MAT *mat0;
	    /* update column j */
	    mat0 = smat->blocks->pmat[j * q + j];	/* diagonal block */
	    if (mat0 == NULL) {
		/* zero diagonal, block which shouldn't have bdry entries */
		sdata += srhs->nvec * smat->blocks->ncols[j];
		continue;
	    }
	    assert(mat0 != NULL && mat0->type != PHG_MATRIX_FREE);
	    rhs0 = phgMapCreateVec(mat0->cmap, srhs->nvec);
	    phgFree(rhs0->data);
	    rhs0->data = sdata;
	    ddata = drhs->data;
	    c0 = smat->blocks->coeff[j * q + j];
	    assert(c0 != 0.);
	    c0 = 1. / c0;
	    for (i = 0; i < p; i++) {
		c = dmat->blocks->coeff[i * q + j];
		mat = dmat->blocks->pmat[i * q + j];
		if (c != 0. && mat != NULL) {
		    assert(mat->type != PHG_MATRIX_FREE);
		    rhs = phgMapCreateVec(mat->rmap, drhs->nvec);
		    phgFree(rhs->data);
		    rhs->data = ddata;
		    assert((mat->handle_bdry_eqns == TRUE) == (i == j));
		    mat->handle_bdry_eqns = TRUE;
		    update_rhs(mat0, rhs0, mat, rhs, c * c0, i == j);
		    mat->handle_bdry_eqns = (i == j);
		    rhs->data = NULL;
		    phgVecDestroy(&rhs);
		}
		ddata += drhs->nvec * smat->blocks->nrows[i];
	    }
	    rhs0->data = NULL;
	    phgVecDestroy(&rhs0);
	    sdata += srhs->nvec * smat->blocks->ncols[j];
	}

	return;
    }

    assert(smat->assembled && dmat->assembled);
    assert(dmat->handle_bdry_eqns);
    /* The diagonal block must be square */
    assert(smat->rmap->nlocal == smat->cmap->nlocal);
    /* The number of columns in the src and dest blocks must match */
    assert(smat->cmap->nlocal == dmat->cmap->nlocal);

    if (!srhs->assembled)
	phgVecAssemble(srhs);

    if (!drhs->assembled)
	phgVecAssemble(drhs);

    mat = dmat;
    rhs = drhs;

#if USE_MPI
    /* prepare communications for getting off-proc bdry entries */
    O2Gsize = mat->localsize - mat->cmap->nlocal;
    rsize = mat->bdry_localsize - mat->cmap->bdry_nlocal;
    rlist = phgAlloc(rsize * sizeof(*rlist));
    rbuff = phgAlloc(rsize * sizeof(*rbuff));
    flags = phgCalloc(O2Gsize, sizeof(*flags));
    scnts = phgAlloc(4 * mat->cmap->nprocs * sizeof(*scnts));
    sdsps = scnts + mat->cmap->nprocs;
    rcnts = sdsps + mat->cmap->nprocs;
    rdsps = rcnts + mat->cmap->nprocs;

    n = mat->cmap->partition[mat->cmap->rank];

    /* expand off-proc bdry entries for faster accessing */
    for (i = mat->cmap->bdry_nlocal; i < mat->bdry_localsize; i++) {
	assert(mat->bdry[i] >= mat->cmap->nlocal);
	flags[mat->bdry[i] - mat->cmap->nlocal] = 1;
    }

    rank = 0;
    bzero(rcnts, mat->cmap->nprocs * sizeof(*rcnts));
    for (i = 0, j = 0; i < O2Gsize; i++) {
	k = mat->ordering[i];
	if (flags[k] == 0)
	    continue;
	k = rlist[j++] = mat->O2Gmap[k];
	assert(j == 1 || rlist[j - 1] > rlist[j - 2]);
	while (rank < mat->cmap->nprocs &&
	       mat->cmap->partition[rank + 1] <= k)
	    rank++;
	assert(rank < mat->cmap->nprocs && rank != mat->cmap->rank);
	rcnts[rank]++;
    }
    assert(j == rsize);
    MPI_Alltoall(rcnts, 1, MPI_INT, scnts, 1, MPI_INT, mat->cmap->comm);
    i = j = 0;
    for (rank = 0; rank < mat->cmap->nprocs; rank++) {
	rdsps[rank] = i;
	sdsps[rank] = j;
	i += rcnts[rank];
	j += scnts[rank];
    }
    /* Note: i/j = total number of entries to recv/send */
    assert(i == rsize);
    slist = phgAlloc(j * sizeof(*slist));
    sbuff = phgAlloc(j * sizeof(*sbuff));
    MPI_Alltoallv(rlist, rcnts, rdsps, PHG_MPI_INT,
		  slist, scnts, sdsps, PHG_MPI_INT, mat->cmap->comm);
    /* collect entries to send */
    for (i = 0; i < j; i++) {
	k = slist[i] - n;
	assert(k >= 0 && k < smat->cmap->nlocal);
#if DEBUG
	/* Note: k must be a boundary entry, otherwise it means incompatibility
	 * of boundary masks in different submeshes */
	if (mat->cmap->bdry_nlocal == 0 ||
		bsearch(&k, mat->cmap->bdry, mat->cmap->bdry_nlocal, sizeof(k),
			phgCompINT) == NULL) {
	    /* identify the source process */
	    j = 0;
	    while (j + 1 < mat->cmap->nprocs && sdsps[j + 1] < i)
		j++;
	    phgError(1, "%s:%d, row %d from p%d is not a bdry row\n",
			__FILE__, __LINE__, k + n, j);
	}
	if (smat->bdry_rows[k].ncols != 1 ||
	    smat->bdry_rows[k].cols[0] != k)
	    phgError(1, "%s:%d, unexpected error: row %d, ncols = %d\n",
			__FILE__, __LINE__, k, smat->bdry_rows[k].ncols);
#endif	/* DEBUG */
	sbuff[i] = srhs->data[k] / smat->bdry_rows[k].data[0];
    }
    phgFree(rlist);
    phgFree(slist);
    MPI_Alltoallv(sbuff, scnts, sdsps, PHG_MPI_FLOAT,
		  rbuff, rcnts, rdsps, PHG_MPI_FLOAT, mat->cmap->comm);
    phgFree(sbuff);
    phgFree(scnts);
    /* save received RHS's in srhs->offp_data */
    phgFree(srhs->offp_data);
    srhs->offp_data = phgAlloc(O2Gsize * sizeof(*srhs->offp_data));
    for (i = 0, j = 0; i < O2Gsize; i++) {
	k = mat->ordering[i];
	if (flags[k] == 0)
	    continue;
	srhs->offp_data[k] = rbuff[j++];
    }
    phgFree(rbuff);
#endif	/* USE_MPI */

    for (i = 0, row = mat->bdry_rows; i < mat->rmap->nlocal; i++, row++) {
	for (j = 0; j < row->ncols; j++) {
	    if ((k = row->cols[j]) == i && diag)
		continue;		/* skip column on diagonal */
#if USE_MPI
	    assert(k < mat->cmap->nlocal || mat->cmap->nprocs > 1);
	    if (k >= mat->cmap->nlocal) {
		k -= mat->cmap->nlocal;
		assert(flags[k] != 0);
		rhs->data[i] -= coeff * srhs->offp_data[k] * row->data[j];
		continue;
	    }
#endif	/* USE_MPI */
	    assert(smat->bdry_rows[k].ncols == 1 &&
		   smat->bdry_rows[k].cols[0] == k);
	    rhs->data[i] -= coeff * srhs->data[k] * row->data[j] /
					smat->bdry_rows[k].data[0];
	}
    }
#if USE_MPI
    phgFree(flags);
    if (O2Gsize > 0)
	bzero(srhs->offp_data, O2Gsize * sizeof(*srhs->offp_data));
#endif	/* USE_MPI */
}

void
phgSolverUpdateRHS(SOLVER *solver)
{
    MAT *A;

    if (solver->rhs_updated)
	return;
    if (!solver->sent_to_oem && !solver->mat->assembled)
	phgMatAssemble(solver->mat);
    A = solver->rhs->mat;
    update_rhs(A, solver->rhs, A, solver->rhs, 1.0, TRUE);

    assert(solver->scaling >= 0 && solver->scaling <= 3);
    if (solver->D != NULL && solver->scaling != 2) {
	INT i;
	/* multiply RHS by sqrt(|diag(solver->mat)|) */
	for (i = 0; i < solver->rhs->map->nlocal; i++)
	    solver->rhs->data[i] *= solver->D[i];
    }

    solver->rhs_updated = TRUE;
}

int
phgSolverResetRHS(SOLVER *solver)
{
    assert(solver->rhs != NULL || solver->mat != NULL);
    if (solver->rhs == NULL) {
	solver->rhs = phgMapCreateVec(solver->mat->rmap, 1);
    }
    else {
	bzero(solver->rhs->data, solver->rhs->nvec *
			solver->rhs->map->nlocal * sizeof(*solver->rhs->data));
#if USE_MPI
	if (solver->rhs->localsize > solver->rhs->map->nlocal)
	    bzero(solver->rhs->offp_data, sizeof(*solver->rhs->offp_data) *
			(solver->rhs->localsize - solver->rhs->map->nlocal));
#endif	/* USE_MPI */
    }
    phgVecDisassemble(solver->rhs);
    solver->rhs_updated = FALSE;
    return 0;
}

/*-------------------------- User interface --------------------------------*/

int
phgSolverInitialize(int *argc, char ***argv)
{
    int i = 0;

    for (i = 0; i < NSolvers; i++)
	if (phgSolverList[i] != NULL && phgSolverList[i]->Initialize != NULL)
	    phgSolverList[i]->Initialize(argc, argv);

    return 0;
}

int
phgSolverFinalize(void)
{
    int i = 0;

    for (i = 0; i < NSolvers; i++)
	if (phgSolverList[i] != NULL && phgSolverList[i]->Finalize != NULL)
	    phgSolverList[i]->Finalize();

    return 0;
}

MAT *
phgSolverGetMat(SOLVER *solver)
{
    solver->mat->refcount++;
    return solver->mat;
}

void
phgSolverSetDefaultSuboptions(void)
/* set default options for all sub-solvers */
{
    phgOptionsSetOptions("+solver_monitor +solver_cond -solver_scaling=off "
			 /*"-solver_symmetry=unsym "*/);
}

SOLVER *
phgSolverMat2Solver(OEM_SOLVER *oem_solver, MAT *mat)
{
    SOLVER *solver = phgCalloc(1, sizeof(*solver));
    int j;

    solver->monitor = monitor;
    solver->scaling = scaling;
    solver->symmetry = scaling == 1 || scaling == 2 ? M_UNSYM : symmetry;
    if (cond) {
	/* prepare for condition number estimate */
	solver->cond = -1.0;
	if (solver->symmetry < M_SYM) {
	    /* create solver->St with the same solver options */
	    MAT dummy;
	    memset(&dummy, 0, sizeof(dummy));
	    dummy.rmap = mat->cmap;
	    dummy.cmap = mat->rmap;
	    phgOptionsPush();
	    phgSolverSetDefaultSuboptions();
	    /* solver->St is either unsymmetric or unneeded */
	    phgOptionsSetKeyword("-solver_symmetry", "unsym");
	    solver->St = phgSolverMat2Solver(oem_solver, &dummy);
	    solver->St->mat = NULL;
	    solver->St->rhs->mat = NULL;
	    phgOptionsPop();
	}
    }

    solver->oem_solver = oem_solver;
    if (solver->oem_solver == NULL && default_solver >= 0) {
	static BOOLEAN warned = FALSE;
	solver->oem_solver = phgSolverList[default_solver];
	if (!warned && mat->rmap->rank == 0 && solver->oem_solver == NULL) {
	    warned = TRUE;
	    phgWarning("solver \"%s\" not available, will use the first "
			"available solver.\n", phgSolverNames[default_solver]);
	}
    }

    /* fall back to first available solver if the requested solver is not
     * available or is SOLVER_DEFAULT */
    if (solver->oem_solver == NULL) {
	for (j = 0; j < NSolvers; j++)
	    if (phgSolverList[j] != NULL) {
		solver->oem_solver = phgSolverList[j];
		break;
	    }
    }

    if (solver->oem_solver == NULL) {
	phgInfo(-1, "phgSolverCreate: no solver available\n");
	phgError(1, "phgSolverCreate: abort.\n");
    }

    if (solver->oem_solver->parallel == FALSE && mat->rmap->nprocs > 1)
	phgError(1, "solver \"%s\" can't be used for parallel systems.\n",
					solver->oem_solver->name);

    phgInfo(1, "using linear solver \"%s\".\n", solver->oem_solver->name);

    mat->refcount++;
    solver->mat = mat;
    solver->rhs = phgMapCreateVec(mat->rmap, 1);
    solver->rhs->mat = mat;	/* reverse pointer to the matrix */
    solver->sent_to_oem = FALSE;
    solver->oem_created = FALSE;
    solver->oem_data = NULL;

    solver->rtol = rtol;
    solver->atol = atol_;
    solver->btol = btol;
    solver->maxit = maxit;
    solver->warn_maxit = warn_maxit;
    solver->matrix_free = matrix_free;
    solver->pc = NULL;

    if (solver->oem_solver->Init != NULL)
	solver->oem_solver->Init(solver);

    return solver;
}

SOLVER *
phgSolverCreate(OEM_SOLVER *oem_solver, DOF *u, ...)
{
    MAP *map;
    MAT *mat;
    int ndof;
    DOF **dofs;
    SOLVER *solver;

    GetVariableArgs(ndof, dofs, u, DOF *);

    map = phgMapCreateN(ndof, dofs);
    phgFree(dofs);
    mat = phgMapCreateMat(map, map);
    phgMapDestroy(&map);

    mat->handle_bdry_eqns = TRUE;

    solver = phgSolverMat2Solver(oem_solver, mat);
    phgMatDestroy(&mat);
    phgVecDisassemble(solver->rhs);
    solver->rhs_updated = FALSE;

    return solver;
}

int
phgSolverDestroy(SOLVER **solver_ptr)
{
    SOLVER *solver = *solver_ptr;

    if (solver == NULL)
	return 0;

    if (solver->D != NULL)
	phgFree(solver->D);

    if (solver->St != NULL)
	phgSolverDestroy(&solver->St);

    solver->oem_solver->Destroy(solver);

    phgMatDestroy(&solver->mat_bak);
    phgMatDestroy(&solver->mat);

    phgVecDestroy(&solver->rhs);

    if (solver->pc != NULL) {
	phgSolverDestroy(&solver->pc->solver_owned);
	phgFree(solver->pc);
    }

    phgFree(solver);

    *solver_ptr = NULL;

    return 0;
}

int
phgSolverAssemble(SOLVER *solver)
{
    MAT *A = NULL;
    int ret = 0;

    FunctionEntry;

    assert(solver->mat != NULL);

    if (solver->assembled)
	Return ret;
    solver->assembled = TRUE;

    if (!solver->mat->assembled)
	phgMatAssemble(solver->mat);

    if (solver->scaling) {
	/* perform diagonal scaling */
	INT i, i0;
	MAT *D, *tmp;
	FLOAT d, dmi = 0.;	/* = |diag(|solver->mat|)| - I|_oo */
	if (solver->monitor)
	    phgPrintf("*** %s: applying %s diagonal scaling.\n", __func__,
			solver->scaling == 1 ? "left" :
			solver->scaling == 2 ? "right" : "symmetric");
	phgMatSetupDiagonal(A = solver->mat);
	D = phgMapCreateMat(A->rmap, A->cmap);
	i0 = D->rmap->partition[D->rmap->rank];
	solver->D = phgAlloc(D->rmap->nlocal * sizeof(*solver->D));
	switch (solver->scaling) {
	    case 3:	/* "sym", A := Sqrt(D)*A*Sqrt(D) */
		for (i = 0; i < D->rmap->nlocal; i++) {
		    solver->D[i] = 1. / Sqrt(Fabs(A->diag[i]));
		    if ((d = Fabs(1. - solver->D[i])) > dmi)
			dmi = d;
		    phgMatAddGlobalEntry(D, i + i0, i + i0, solver->D[i]);
		}
		phgMatAssemble(D);
		tmp = phgMatMat(MAT_OP_N, MAT_OP_N, 1., D, A, 0.0, NULL);
		A = phgMatMat(MAT_OP_N, MAT_OP_N, 1., tmp, D, 0.0, NULL);
		phgMatDestroy(&tmp);
		break;
	    case 2:	/* "right", A := A*D */
		for (i = 0; i < D->rmap->nlocal; i++) {
		    solver->D[i] = 1. / Fabs(A->diag[i]);
		    if ((d = Fabs(1. - solver->D[i])) > dmi)
			dmi = d;
		    phgMatAddGlobalEntry(D, i + i0, i + i0, solver->D[i]);
		}
		tmp = A;
		A = phgMatMat(MAT_OP_N, MAT_OP_N, 1., tmp, D, 0.0, NULL);
		break;
	    case 1:	/* "left", A := D*A */
		for (i = 0; i < D->rmap->nlocal; i++) {
		    solver->D[i] = 1. / Fabs(A->diag[i]);
		    if ((d = Fabs(1. - solver->D[i])) > dmi)
			dmi = d;
		    phgMatAddGlobalEntry(D, i + i0, i + i0, solver->D[i]);
		}
		tmp = A;
		A = phgMatMat(MAT_OP_N, MAT_OP_N, 1., D, tmp, 0.0, NULL);
		break;
	    default:
		phgError(1, "%s:%d: unexpected.\n", __FILE__, __LINE__);
	}
#if 0
#warning
/* save mat before scaling */
phgMatDumpMATLAB(solver->mat, "A0", "A0.m");
/* save D */
phgMatDumpMATLAB(D, "D", "D.m");
#endif
	phgMatDestroy(&D);

	/* only apply scaling when D != I */
	if (dmi > FLOAT_EPSILON * 100.0) {
	    /* Note: we will need solver->mat->*bdry* and solver->mat->O2G* in
	     * update_rhs(), and we cannot simply copy solver->mat->*bdry* to A
	     * since they are not compatible with A->O2Gmap and A->ordering,
	     * so we save solver->mat in solver->mat_bak (freeing its data) */
	    phgMatFreeMatrix(solver->mat);
	    solver->mat_bak = solver->mat;
	    solver->mat = A;
	    if (solver->scaling == 1 || solver->scaling == 2)
		solver->symmetry = M_UNSYM;
	    else
		phgSolverDestroy(&solver->St);
	}
	else {
	    solver->scaling = 0;
	    phgFree(solver->D);
	    solver->D = NULL;
	}
    }

    if (solver->symmetry < solver->oem_solver->req_symmetry &&
	solver->mat->rmap->rank == 0)
	phgWarning("specified symmetry (%s) doesn't match requirement (%s) of "
		   "the %s solver, results may be wrong!!!\n",
		   symmetry_args[solver->symmetry],
		   symmetry_args[solver->oem_solver->req_symmetry],
		   solver->oem_solver->name);

    if (solver->cond < 0.0) {
	/* protect solver->mat (=> solver->mat->ref_count++) */
	A = phgSolverGetMat(solver);
	/* setup solver->St */
	if (solver->St != NULL)
	    solver->St->mat = solver->St->rhs->mat = phgMatTranspose(A);
    }

    if (!solver->sent_to_oem)
	send_matrix_to_oem_solver(solver);

    ret = solver->oem_solver->Assemble == NULL ?
			0 : solver->oem_solver->Assemble(solver);

    if (solver->cond < 0.0) {
	/* estimate condition number */
	double t = phgGetTime(NULL);
	solver->cond = 0.0;	/* prevent dead-loop */
#if 1
	/* share solvers in phgMatConditionNumber */
	solver->cond = phgMatConditionNumber_(A, solver->symmetry,
							solver, solver->St);
#else
#warning
	/* create new solvers in phgMatConditionNumber */
	solver->cond = phgMatConditionNumber_(A, solver->symmetry, NULL, NULL);
#endif
	phgSolverDestroy(&solver->St);
	phgMatDestroy(&A);
	if (solver->monitor)
	    phgPrintf("*** %s: condition number: %0.2le, time: %0.2lg\n",
			__func__, (double)solver->cond, phgGetTime(NULL) - t);
    }

    Return ret;
}

int
phgSolverVecSolve(SOLVER *solver, BOOLEAN destroy, VEC *x)
/* returns number of iterations, < 0 ==> error */
{
    int ret;
    INT i, j;
    FLOAT *rhs;
    size_t mem_peak = 0;
    double time0 = 0.;

    assert(solver->rhs != NULL && solver->mat != NULL);

    if (solver->mat->cmap->nlocal != x->map->nlocal) {
	phgError(1, "(%s:%d) inconsistent RHS!\n", __FILE__, __LINE__);
    }

    if (solver->mat->handle_bdry_eqns == TRUE &&
	   solver->rhs->mat != solver->mat)
	phgError(1, "(%s:%d) solver->rhs->mat != solver->mat.\n",
			__FILE__, __LINE__);

    phgSolverAssemble(solver);

    if (!solver->rhs->assembled)
	phgVecAssemble(solver->rhs);

    if (!solver->rhs_updated)
	phgSolverUpdateRHS(solver);

    /* send RHS to OEM solver and free PHG's RHS */
    if (solver->oem_solver->AddRHSEntries != NULL) {
	j = solver->mat->rmap->partition[solver->mat->rmap->rank];
	rhs = solver->rhs->data;
	for (i = 0; i < solver->mat->rmap->nlocal; i++, j++)
	    solver->oem_solver->AddRHSEntries(solver, 1, &j, rhs++);
    }

    assert(solver->scaling >= 0 && solver->scaling <= 3);
    if (solver->scaling > 1) {	/* "right" or "sym" */
	assert(solver->D != NULL);
	/* divide initial solution by solver->D */
	for (i = 0; i < solver->rhs->map->nlocal; i++)
	    x->data[i] /= solver->D[i];
    }

    if (solver->monitor) {
	mem_peak = phgMemoryPeakReset(0);
	time0 = phgGetTime(NULL);
    }
    ret = solver->oem_solver->Solve(solver, x, destroy);
    if (solver->monitor) {
	mem_peak = phgMemoryPeakRestore(mem_peak);
	phgPrintf("%s: solution time: %0.4lg, memory: %0.4lgMB.\n",
		solver->oem_solver->name,
		phgGetTime(NULL) - time0, mem_peak / (1024.0 * 1024.0));
    }

    if (x->map->rank == 0) {
	if (ret < 0)
	    phgError(1, "error occured in OEM solver \"%s\".\n",
		     solver->oem_solver->name);
	if (solver->warn_maxit && ret >= solver->maxit && solver->rtol > 0.0 &&
	    solver->oem_solver->iterative)
	    phgWarning("maxit attained in OEM solver, "
		       "the result may be inaccurate.\n");
    }

    assert(solver->scaling >= 0 && solver->scaling <= 3);
    if (solver->scaling > 1) {	/* "right" or "sym" */
	assert(solver->D != NULL);
	/* multiply x by solver->D */
	for (i = 0; i < solver->rhs->map->nlocal; i++)
	    x->data[i] *= solver->D[i];
    }

    return ret;
}

int
phgSolverSolve(SOLVER *solver, BOOLEAN destroy, DOF *u, ...)
/* returns number of iterations, < 0 ==> error */
{
    int i, ndof;
    DOF **dofs;
    VEC *x;

    GetVariableArgs(ndof, dofs, u, DOF *);

    /* make sure the DOFs match those used in creating the solver */
    assert(solver->rhs->map->dofs != NULL);
    assert(solver->rhs->map->ndof == ndof);
    for (i = 0; i < ndof; i++) {
	if (dofs[i]->type != solver->rhs->map->dofs[i]->type)
	    phgError(1, "phgSolveSolve: DOF %d (%s), types mismatch: %s<=>%s\n",
		     i, dofs[i]->name, dofs[i]->type,
		     solver->rhs->map->dofs[i]->type);
	if (dofs[i]->dim != solver->rhs->map->dofs[i]->dim)
	    phgError(1, "phgSolveSolve: DOF %d (%s), dims mismatch: %d<=>%d\n",
		     i, dofs[i]->name, dofs[i]->dim,
		     solver->rhs->map->dofs[i]->dim);
    }

    x = phgMapCreateVec(solver->rhs->map, 1);
    /* copy initial solution to x */
    phgMapDofToLocalData(x->map, ndof, dofs, x->data);
    /* call solver */
    i = phgSolverVecSolve(solver, destroy, x);
    /* copy solution from x */
    phgMapLocalDataToDof(x->map, ndof, dofs, x->data);
    phgFree(dofs);
    phgVecDestroy(&x);

    return i;
}

/* The functions below are maintained only for backward compatibility */

void
phgSolverSetTol(SOLVER *solver, FLOAT user_tol)
{
#if USE_MPI
    FLOAT tol = user_tol;

    MPI_Bcast(&tol, 1, PHG_MPI_FLOAT, 0, solver->mat->rmap->comm);
    user_tol = tol;
#endif
    solver->rtol = user_tol;
}

void
phgSolverSetMaxIt(SOLVER *solver, int user_maxit)
{
#if USE_MPI
    int maxit = user_maxit;

    MPI_Bcast(&maxit, 1, MPI_INT, 0, solver->mat->rmap->comm);
    user_maxit = maxit;
#endif
    solver->maxit = user_maxit;
}

#endif	/* !defined(PHG_TO_P4EST) */

#ifndef DEBUG_RHS
# define DEBUG_RHS 0	/* 0 or 1 */
#endif

static FLOAT (*coords)[Dim + 1 + DEBUG_RHS] = NULL;
static int
comp_coord(const void *p0, const void *p1)
{
    FLOAT *c0 = coords[*(const INT *)p0], *c1 = coords[*(const INT *)p1];
    int i;

    for (i = Dim; i >= 0; i--)
	if (Fabs(c0[i] - c1[i]) > FLOAT_EPSILON * 100.0)
	    return c0[i] < c1[i] ? -1 : 1;

    return 0;
}

void
phgSolverDumpMATLAB_(SOLVER *solver, const char *mat_name, const char *rhs_name,
		     const char *perm_name, BOOLEAN reorder)
/* writes .m files for the matrix and vector.
 *
 * "mat_name" and "rhs_name" are respectively the var name for the mat and rhs.
 *
 * "perm_name", if not NULL, is the var name for the permutation vector.
 * It is only used when "reorder"==TRUE.
 *
 * If "reorder"==TRUE then the unknowns will be reordered in a way which is
 * independent of grid partitioning and element indices, for parallel debugging.
 */
{
    FLOAT (*c)[Dim + 1 + DEBUG_RHS];
    MAP *map, *m;
    GRID *g;
    ELEMENT *e;
    MAT *A = NULL;
    VEC *b = NULL, *p = NULL;
    int ii, dof_id, *counts = NULL, *displs = NULL;
    INT i, *perm;
    char *flags;

    if (mat_name == NULL && rhs_name == NULL &&
	(reorder == FALSE || perm_name == NULL))
	return;		/* nothing to do */

    map = solver->rhs->map;
    assert(map->dofs != NULL);
    g = ((DOF *)map->dofs[0])->g;

    if (rhs_name != NULL)
	phgVecAssemble(solver->rhs);
    if (mat_name != NULL)
	phgMatAssemble(solver->mat);

    if (!reorder) {
	if (mat_name != NULL) {
	    phgMatDumpMATLAB(solver->mat, mat_name, mat_name);
	    phgPrintf("%s: matrix saved in \"%s.m\".\n", __func__, mat_name);
	}
	if (rhs_name != NULL) {
	    phgVecDumpMATLAB(solver->rhs, rhs_name, rhs_name);
	    phgPrintf("%s: RHS saved in \"%s.m\".\n", __func__, rhs_name);
	}
	return;
    }

    assert((size_t)map->nglobal * (size_t)(Dim + 1 + DEBUG_RHS) <= INT_MAX);

    c = phgCalloc(map->nlocal + 2, sizeof(*c));	/* 2 entries for work space */
    coords = c;		/* for comp_cord() */
    flags = phgCalloc(map->nlocal, sizeof(*flags));	/* TRUE => processed */
    ForAllElements(g, e) {
	for (dof_id = 0; dof_id < map->ndof; dof_id++) {
	    DOF *u = (DOF *)map->dofs[dof_id];
	    int n = DofNBas(u, e), bindex;
#ifndef PHG_TO_P4EST
	    static BOOLEAN warned = FALSE;
	    int gindex;
	    GTYPE gtype;
#endif	/* !defined(PHG_TO_P4EST) */

	    assert(u->type != NULL);
	    for (ii = 0; ii < n; ii++) {
		i = phgSolverMapE2G(solver, dof_id, e, ii)
					- map->partition[g->rank];
		if (i < 0 || i >= map->nlocal || flags[i])
		    continue;
		flags[i] = 1;	/* process each DOF only once */
#if defined(PHG_TO_P4EST)	/*--------------------------------------------*/
		/* DG only for p4est element, use the lower corner */
		c[i][1] = e->corners[0][0];
		c[i][2] = e->corners[0][1];
		c[i][3] = e->corners[0][2];
		bindex = ii;
#else	/* defined(PHG_TO_P4EST) ---------------------------------------------*/
		phgDofGetElementBasisInfo_(u, e, ii / u->dim,
				&gtype, &gindex, &bindex,
				!u->type->is_nodal ? NULL : c[i] + 1, NULL);
		if (u->type->np_vert+u->type->np_edge+u->type->np_face == 0) {
		    /* force gtype to VOLUME */
		    gtype = VOLUME;
		    gindex = 0;
		    bindex = ii;
		}
		else if (!u->type->is_nodal && !warned) {
		    /* TODO: need to orient faces/edges according to physical
		     * coordinates, instead of global numbering, of vertices
		     * for edge and face DOFs.
		     *
		     * For example, for the following 2 cases (with cube5.dat):
		     *	    mpirun -np 1 ./ipdg -refine0 1 -dof_type HB3
		     *	    mpirun -np 2 ./ipdg -refine0 1 -dof_type HB3
		     * the element with global index 3 gets different global
		     * vertex numbers: {1,6,9,8}, for the 1st case and
		     * {1,6,8,9} for the 2nd case. And the basis function 15
		     * (1 on edge 5) for the 2 cases has opposite sign. */
		    warned = TRUE;
		    phgPrintf("*** %s: WARNING: reorder=TRUE may not work "
			      "with \"%s\".\n", __func__, u->type->name);
		}
		if (gtype == VOLUME || !u->type->is_nodal) {
		    FLOAT lam[] = {0.,0.,0.,0.};
		    if (gtype == VOLUME) {
			assert((g->types_elem[e->index] & OWNER));
			lam[0] = lam[1] = lam[2] = lam[3] = 0.25;
		    }
		    else if (gtype == FACE) {
			assert((g->types_face[e->faces[gindex]] & OWNER));
			lam[0] = lam[1] = lam[2] = lam[3] = 1. / 3.;
			lam[gindex] = 0.;
			/* TODO: map bindex */
		    }
		    else if (gtype == EDGE) {
			int v0, v1;
			assert((g->types_edge[e->edges[gindex]] & OWNER));
			GetEdgeVertices(e, gindex, v0, v1);
			lam[v0] = 0.5;
			lam[v1] = 0.5;
#if 0
			/* map bindex */
			INT i0 = i, i1 = map->nlocal;
			c[i0][0] = c[i1][0] = 0.;
			memcpy(c[i0]+1, g->verts[e->verts[v0]], sizeof(COORD));
			memcpy(c[i1]+1, g->verts[e->verts[v1]], sizeof(COORD));
#if 0
#warning
INT j = i + map->partition[g->rank];
/* Test case: ipdg -refine0 -dof_type HB3 */
if ((g->nprocs == 1 && (j == 76)) || (g->nprocs == 2 && (j == 96))) {
    phgInfo(-1, "--- %d: e = %d, ii = %d, gt = %d, gi = %d\n", j, GlobalElement(g, e->index), ii, gtype, gindex);
    phgInfo(-1, "*** %d: v0=(%g,%g,%g), v1=(%g,%g,%g), bi=%d, cmp=%d, rhs=%g\n",
	j, c[i0][1], c[i0][2], c[i0][3], c[i1][1], c[i1][2], c[i1][3], bindex,
	comp_coord(&i0, &i1), solver->rhs->data[i]);
}
#endif
			if (comp_coord(&i0, &i1) > 0)
			    bindex = bindex % u->dim + (u->type->np_edge - 1
						- bindex / u->dim) * u->dim;
#endif
		    }
		    else {
			assert(gtype == VERTEX);
			assert((g->types_vert[e->verts[gindex]] & OWNER));
			lam[gindex] = 1.0;
		    }
		    phgGeomLambda2XYZ(g, e, lam, &c[i][1], &c[i][2], &c[i][3]);
		}
#endif	/* defined(PHG_TO_P4EST) ---------------------------------------------*/
		c[i][0] = 1000000.0 * dof_id + bindex;
#if DEBUG_RHS
		c[i][Dim + 1] = solver->rhs->data[i];
#endif
	    }
	}
    }
    phgFree(flags);

    if (g->rank == 0) {
	coords = phgAlloc(map->nglobal * sizeof(*coords));
	counts = phgAlloc(2 * g->nprocs * sizeof(*counts));
	displs = counts + g->nprocs;
	for (ii = 0; ii < g->nprocs; ii++) {
	    displs[ii] = map->partition[ii] * (Dim + 1 + DEBUG_RHS);
	    counts[ii] = (map->partition[ii + 1]-map->partition[ii])
			 * (Dim + 1 + DEBUG_RHS);
	}
    }
#if USE_MPI
    MPI_Gatherv(c, map->nlocal * (Dim + 1 + DEBUG_RHS), PHG_MPI_FLOAT,
		coords, counts, displs, PHG_MPI_FLOAT, 0, g->comm);
    phgFree(c);
#else	/* USE_MPI */
    phgFree(coords);
    coords = c;
#endif	/* USE_MPI */
    perm = phgAlloc(map->nglobal * sizeof(*perm));
    if (g->rank == 0) {
	INT *tmp = phgAlloc(map->nglobal * sizeof(*tmp));
	phgFree(counts);
	for (i = 0; i < map->nglobal; i++)
	    tmp[i] = i;
	qsort(tmp, map->nglobal, sizeof(*tmp), comp_coord);
#if 0
#warning
for (i = 74; i < 75/*map->nglobal*/; i++)
phgInfo(-1, "%"dFMT": coords[%"dFMT"] = %d %g %g %g (%g)\n", i, tmp[i], (int)(coords[tmp[i]][0]+0.5), (double)coords[tmp[i]][1], (double)coords[tmp[i]][2], (double)coords[tmp[i]][3], DEBUG_RHS ? (double)coords[tmp[i]][4] : 0.0);
#endif
	/* perm := tmp^-1 */
	for (i = 0; i < map->nglobal; i++)
	    perm[tmp[i]] = i;
	phgFree(tmp);
	phgFree(coords);
	coords = NULL;
    }
#if USE_MPI
    MPI_Bcast(perm, map->nglobal, PHG_MPI_INT, 0, g->comm);
#endif	/* USE_MPI */

    m = phgMapCreateSimpleMap(g->comm, map->nlocal, map->nglobal);
    if (rhs_name != NULL) {
	b = phgMapCreateVec(m, 1);
	phgVecDisassemble(b);
    }
    if (mat_name != NULL) {
	A = phgMapCreateMat(m, m);
    }
    if (perm_name != NULL) {
	p = phgMapCreateVec(m, 1);
	phgVecDisassemble(p);
    }
    phgMapDestroy(&m);

    for (i = 0; i < map->nlocal; i++) {
	INT I = perm[map->partition[g->rank] + i];
	if (mat_name != NULL) {
	    const MAT_ROW *row = phgMatGetRow(solver->mat, i);
	    for (ii = 0; ii < row->ncols; ii++) {
		phgMatAddGlobalEntry(A, I, perm[row->cols[ii]], row->data[ii]);
	    }
	}
	if (rhs_name != NULL)
	    phgVecAddGlobalEntry(b, 0, I, solver->rhs->data[i]);
	if (perm_name != NULL)
	    phgVecAddEntry(p, 0, i, (FLOAT)(I + 1));
    }
    phgFree(perm);

    if (mat_name != NULL) {
	phgMatDumpMATLAB(A, mat_name, mat_name);
	phgPrintf("%s: permuted matrix saved in \"%s.m\".\n",
					__func__, mat_name);
	phgMatDestroy(&A);
    }

    if (rhs_name != NULL) {
	phgVecDumpMATLAB(b, rhs_name, rhs_name);
	phgPrintf("%s: permuted RHS vector saved in \"%s.m\".\n",
					__func__, rhs_name);
	phgVecDestroy(&b);
    }

    if (perm_name != NULL) {
	phgVecDumpMATLAB(p, perm_name, perm_name);
	phgPrintf("%s: permutation [1..%d]->[1..%d] saved in \"%s.m\".\n",
			__func__, p->map->nglobal, p->map->nglobal, perm_name);
	phgVecDestroy(&p);
    }
}
