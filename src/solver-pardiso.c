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

/* $Id: solver-pardiso.c,v 1.26 2022/09/21 02:35:28 zlb Exp $
 *
 * Interface to the PARDISO solver */

#include "phg.h"

#if USE_PARDISO		/* whole file */

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

#include <stdlib.h>
#include <string.h>

#define PARDISOINIT	F77_FUNC(pardisoinit,PARDISOINIT)
#define PARDISO		F77_FUNC(pardiso,PARDISO)

/* Note: according to the manual pt and iparm are int32 arrays on 32-bit
 * architectures and int64 arrays on 64-bit architectures, so they should
 * be casted to size_t which works for both 32 and 64 bit architectures.
 * But apperently this is not true on em64t so we cast them to 'int'. */
typedef int PARM_t;

#ifdef __cplusplus
extern "C" {
#endif
  void PARDISOINIT(PARM_t pt[], int *mtype, PARM_t iparm[]);
  void PARDISO(PARM_t pt[], int *maxfct, int *mnum, int *mtype, int *phase,
	       int *n, double *a, int *ia, int *ja, int *perm, int *nrhs,
	       PARM_t iparm[], int *msglvl, double *b, double *x, int *error);
#ifdef __cplusplus
}
#endif

#if 0
/* debugging PARDISO calls */
# define PARDISOINIT_(pt, mtype, iparm) \
	phgInfo(-1, "calling PARDISOINIT with mtype = %d\n", *(mtype)); \
	PARDISOINIT(pt, mtype, iparm);
# define PARDISO_(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, \
		  iparm, msglvl, b, x, error) \
	phgInfo(-1, "calling PARDISO with phase = %d, mtype = %d\n", \
			*(phase), *(mtype)); \
	PARDISO(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, \
                iparm, msglvl, b, x, error);
#else
# define PARDISOINIT_	PARDISOINIT
# define PARDISO_	PARDISO
#endif

static INT refine = 0;

typedef struct {
    PARM_t	pt[64], iparm[64];
    int		*Ap, *Ai;
    double	*Ad;
    int		mtype;		/* matrix type */
    int		refine;
    BOOLEAN	factored;
} OEM_DATA;

#define	_t	((OEM_DATA *)solver->oem_data)

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nPARDISO options:", "\n", "pardiso");
    phgOptionsRegisterInt("-pardiso_refine", "Max number of iterative "
				"refinements", &refine);

    return 0;
}

static int
Init(SOLVER *solver)
{
    /* Symmetry:
	 1 real and structurally symmetric
	 2 real and symmetric positive definite
	-2 real and symmetric indefinite
	 3 complex and structurally symmetric
	 4 complex and Hermitian positive definite
	-4 complex and Hermitian indefinite
	 6 complex and symmetric
	11 real and nonsymmetric
	13 complex and nonsymmetric */
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));

    _t->refine = refine;    

    return 0;
}

static int
Create(SOLVER *solver)
{
#define IS_SYMMETRIC(mtype)	(mtype == 2 || mtype == -2)
    switch (solver->symmetry) {
	case M_UNSYM:	_t->mtype = 11; break;
	case M_SSYM:	_t->mtype = 1; break;
	case M_SYM:	_t->mtype = -2; break;
	case M_SPD:	_t->mtype = 2; break;
	default:	phgError(1, "invalid symmetry %d\n", solver->symmetry);
    }

    PARDISOINIT_(_t->pt, &_t->mtype, _t->iparm);

    _t->factored = FALSE;

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    int maxfct = 1, mnum = 1, msglvl = solver->monitor, phase;
    int nrhs = 1, error = 0, idummy, nlocal = solver->mat->rmap->nlocal;
    double ddummy;

    if (solver->oem_data == NULL)
	return 0;

    if (_t->Ap != NULL) {
	phase = -1;
	PARDISO_(_t->pt, &maxfct, &mnum, &_t->mtype, &phase,
		 &nlocal, _t->Ad, _t->Ap, _t->Ai,
		 &idummy, &nrhs, _t->iparm, &msglvl, &ddummy, &ddummy, &error);
	phgFree(_t->Ap);
	phgFree(_t->Ai);
	phgFree(_t->Ad);
    }

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static INT *sort_key = NULL;

static int
comp_cols(const void *p0, const void *p1)
{
    return sort_key[*(const INT *)p0] - sort_key[*(const INT *)p1];
}

static int
Assemble(SOLVER *solver)
{
    MAP *map = solver->mat->rmap;
    INT i, j, k, m, *index;
    int n = map->nglobal, nnz = solver->mat->nnz_d, omp_threads;
    const MAT_ROW *row;
    double time0;
    char *env;

    int maxfct = 1, mnum = 1, msglvl = solver->monitor, phase;
    int nrhs = 1, error = 0, idummy;
    double ddummy;

    if (_t->factored)
	return 0;

    _t->factored = TRUE;

    if (solver->sent_to_oem == TRUE)
	phgError(1, "%s:%d, data already sent to OEM solver!\n",
						__FILE__, __LINE__);
    solver->sent_to_oem = TRUE;

    assert(sizeof(INT) == sizeof(int));

    assert(solver->mat->type != PHG_DESTROYED);
    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }

    phgInfo(2, "send data to PARDISO solver.\n");

    /* allocate index array */
    index = phgAlloc(n * sizeof(*index));
    _t->Ap = phgCalloc(n + 1, sizeof(*_t->Ap));
    _t->Ai = phgAlloc(nnz * sizeof(*_t->Ai));
    _t->Ad = phgAlloc(nnz * sizeof(*_t->Ad));

    _t->Ap[0] = 1;
    j = 0;
    for (i = 0; i < n; i++) {
	row = phgMatGetRow(solver->mat, i);
	m = row->ncols;
	if (m == 1) {
	    if (!IS_SYMMETRIC(_t->mtype) || row->cols[0] >= i) {
		_t->Ai[j] = row->cols[0] + 1;
		_t->Ad[j++] = row->data[0];
	    }
	}
	else if (m > 0) {
	    for (k = 0; k < m; k++)
		index[k] = k;
	    sort_key = row->cols;
	    qsort(index, m, sizeof(*index), comp_cols);
	    k = 0;
	    if (IS_SYMMETRIC(_t->mtype)) {
		/* Note: may use binary search */
		while (row->cols[index[k]] < i)
		    k++;
	    }
	    for (; k < m; k++) {
		_t->Ai[j] = row->cols[index[k]] + 1;
		_t->Ad[j++] = row->data[index[k]];
	    }
	}
	_t->Ap[i + 1] = j + 1;
	if (solver->mat->refcount == 0 && solver->mat->rows != NULL) {
	    /* free matrix row */
	    phgFree(solver->mat->rows[i].cols);
	    phgFree(solver->mat->rows[i].data);
	    solver->mat->rows[i].cols = NULL;
	    solver->mat->rows[i].data = NULL;
	    solver->mat->rows[i].ncols = solver->mat->rows[i].alloc = 0;
	}
    }
    if (IS_SYMMETRIC(_t->mtype)) {
	nnz = j;
	_t->Ai = phgRealloc_(_t->Ai, nnz * sizeof(*_t->Ai),
				     nnz * sizeof(*_t->Ai));
	_t->Ad = phgRealloc_(_t->Ad, nnz * sizeof(*_t->Ad),
				     nnz * sizeof(*_t->Ad));
    }
    assert(j == nnz);
    phgFree(index);

    if (solver->mat->refcount == 0)
	phgMatFreeMatrix(solver->mat);

    /* determine number of threads to use */
    if ((env = getenv("OMP_NUM_THREADS")) != NULL)
	sscanf(env, "%d", &omp_threads);
    _t->iparm[2] = omp_threads <= 0 ? 1 : omp_threads;

    phgInfo(2, "analyse and factor matrix.\n");

    /* symbolic factorization */
    phase = 11;
    time0 = phgGetTime(NULL);
    PARDISO_(_t->pt, &maxfct, &mnum, &_t->mtype, &phase,
	     &n, _t->Ad, _t->Ap, _t->Ai, &idummy, &nrhs,
	     _t->iparm, &msglvl, &ddummy, &ddummy, &error);
    if (error != 0)
	phgError(1, "error during symbolic factorization: %d", error);
    phgInfo(1, "Symbolic factorization time: %lg\n",
			(double)(phgGetTime(NULL) - time0));
#if 0
    phgInfo(1, "Number of nonzeros in factors  = %d\n", (int)_t->iparm[17]);
    phgInfo(1, "Number of factorization MFLOPS = %d\n", (int)_t->iparm[18]);
#endif

    /* numeric factorization */
    phase = 22;
    time0 = phgGetTime(NULL);
    PARDISO_(_t->pt, &maxfct, &mnum, &_t->mtype, &phase,
	     &n, _t->Ad, _t->Ap, _t->Ai, &idummy, &nrhs,
	     _t->iparm, &msglvl, &ddummy, &ddummy, &error);
    if (error != 0)
	phgError(1, "error during numeric factorization: %d", error);
    phgInfo(1, "Numeric factorization time: %lg\n",
			(double)(phgGetTime(NULL) - time0));

    return 0;
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    double time0;

    int maxfct = 1, mnum = 1, msglvl = solver->monitor, phase;
    int nrhs = 1, error = 0, idummy, nlocal = solver->mat->rmap->nlocal;
    double ddummy;

    Assemble(solver);

    phase = 33;
    _t->iparm[7] = _t->refine;		/* # iterative refinements */
    time0 = phgGetTime(NULL);
#if FT_PHG == FT_DOUBLE
    PARDISO_(_t->pt, &maxfct, &mnum, &_t->mtype, &phase,
	    &nlocal, _t->Ad, _t->Ap, _t->Ai, &idummy, &nrhs,
	    _t->iparm, &msglvl, solver->rhs->data, x->data, &error);
#else	/* FT_PHG == FT_DOUBLE */
    {
	int i, n = solver->mat->rmap->nlocal;
	double *rhs, *sol;
	rhs = phgAlloc(2 * n * sizeof(*rhs));
	sol = rhs + n;
	for (i = 0; i < n; i++) {
	    rhs[i] = (double)solver->rhs->data[i];
	    sol[i] = (double)x->data[i];
	}
	PARDISO_(_t->pt, &maxfct, &mnum, &_t->mtype, &phase,
	    &nlocal, _t->Ad, _t->Ap, _t->Ai, &idummy, &nrhs,
	    _t->iparm, &msglvl, rhs, sol, &error);
	for (i = 0; i < n; i++)
	    x->data[i] = (FLOAT)sol[i];
	phgFree(rhs);
    }
#endif	/* FT_PHG == FT_DOUBLE */
    if (error != 0)
	phgError(1, "error during solution: %d", error);
    phgInfo(1, "Time for solution: %lg\n", (double)(phgGetTime(NULL) - time0));

    if (destroy) {
	phase = -1;
	PARDISO_(_t->pt, &maxfct, &mnum, &_t->mtype, &phase,
		 &nlocal, _t->Ad, _t->Ap, _t->Ai, &idummy,
		 &nrhs, _t->iparm, &msglvl, &ddummy, &ddummy, &error);
	phgFree(_t->Ap);
	_t->Ap = NULL;
	phgFree(_t->Ai);
	_t->Ai = NULL;
	phgFree(_t->Ad);
	_t->Ad = NULL;
    }

    solver->residual = 0.0;	/* fake residual norm */
    return solver->nits = 1;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverPARDISO_ = {
    "PARDISO", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy,
    AddMatrixEntries, AddRHSEntries,
    Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, FALSE, FALSE, FALSE, FALSE
};

#endif	/* USE_PARDISO */
