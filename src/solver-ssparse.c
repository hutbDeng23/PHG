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

/* $Id: solver-ssparse.c,v 1.14 2022/09/21 02:35:28 zlb Exp $ */

/* Interface to the SuiteSparse solvers:
 *	http://www.cise.ufl.edu/research/sparse/SuiteSparse
 *
 * TODO: more solvers (currently only UMFPACK is implemented)
 * */

#include "phg.h"

#if USE_SSPARSE		/* whole file */

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

#include <umfpack.h>

#include <stdlib.h>
#include <string.h>

typedef struct {
    void	*symbolic, *numeric;
    int		solver_id;		/* for SuiteSParse solver selection */
    int		*Ap, *Ai;
    double	*Ad;
    BOOLEAN	factored;
} OEM_DATA;

#define	_t	((OEM_DATA *)solver->oem_data)

static int solver_id = 0;

static int
RegisterOptions(void)
{
#if 0
    phgOptionsRegisterTitle("\nThe SuiteSparse solvers options:", "\n", "ssparse");
    phgOptionsRegisterKeyword("ssparse_solver", "SuiteSparse solver",
						&solver_names, solver_id);
#endif

    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));
    _t->solver_id = solver_id;

    return 0;
}

static int
Create(SOLVER *solver)
{
    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }
    _t->factored = FALSE;

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    if (_t->symbolic != NULL)
	umfpack_di_free_symbolic(&_t->symbolic);

    if (_t->numeric != NULL)
	umfpack_di_free_numeric(&_t->numeric);

    phgFree(_t->Ap);
    phgFree(_t->Ai);
    phgFree(_t->Ad);

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    MAP *map = solver->mat->rmap;
    INT i, j, m;
    int n = map->nglobal, nnz = solver->mat->nnz_d;
    int status;
    MAT_ROW *row;
    double time0;

    if (_t->factored)
	return 0;

    _t->factored = TRUE;

    if (solver->sent_to_oem == TRUE)
	phgError(1, "%s:%d, data already sent to OEM solver!\n",
						__FILE__, __LINE__);
    solver->sent_to_oem = TRUE;

    assert(FT_PHG == FT_DOUBLE);
    assert(sizeof(INT) == sizeof(int));
    assert(solver->mat->type == PHG_PACKED ||
	   solver->mat->type == PHG_UNPACKED);
    if (solver->mat->type == PHG_PACKED)
	phgMatUnpack(solver->mat);

    phgInfo(2, "send data to SuiteSparse solver.\n");

    /* Note: UMFPACK uses CSC storage with sorted row indices in each row */

    /* allocate index array */
    _t->Ai = phgAlloc(nnz * sizeof(*_t->Ai));
    _t->Ad = phgAlloc(nnz * sizeof(*_t->Ad));
    _t->Ap = phgCalloc(n + 1, sizeof(*_t->Ap));

    /* count # of nonzeros in the columns */
    for (i = 0, row = solver->mat->rows; i < n; i++, row++) {
	m = row->ncols;
	for (j = 0; j < m; j++) {
	    _t->Ap[row->cols[j] + 1]++;
	}
    }

    for (i = 0; i < n; i++)
	_t->Ap[i + 1] += _t->Ap[i];

    for (i = 0, row = solver->mat->rows; i < n; i++, row++) {
	m = row->ncols;
	for (j = 0; j < m; j++) {
	    _t->Ai[_t->Ap[row->cols[j]]] = i;
	    _t->Ad[_t->Ap[row->cols[j]]] = row->data[j];
	    _t->Ap[row->cols[j]]++;
	}
	if (solver->mat->refcount == 0) {
	    /* free matrix row */
	    phgFree(row->cols);
	    phgFree(row->data);
	    row->cols = NULL;
	    row->data = NULL;
	    row->ncols = row->alloc = 0;
	}
    }

    assert(_t->Ap[n] == nnz);
    for (i = n; i > 0; i--)
	_t->Ap[i] = _t->Ap[i - 1];
    _t->Ap[0] = 0;

    if (solver->mat->refcount == 0)
	phgMatFreeMatrix(solver->mat);

    phgInfo(2, "analyse and factor matrix.\n");
    time0 = phgGetTime(NULL);
    status = umfpack_di_symbolic(n, n, _t->Ap, _t->Ai, _t->Ad, &_t->symbolic,
					NULL, NULL) ;
    if (status != UMFPACK_OK)
	phgError(1, "error calling umfpack_di_symbolic(), status = %d.\n",
					status);
    phgInfo(1, "Time for symbolic factorize: %lg\n",
			(double)(phgGetTime(NULL) - time0));
    time0 = phgGetTime(NULL);
    status = umfpack_di_numeric(_t->Ap, _t->Ai, _t->Ad, _t->symbolic,
					&_t->numeric, NULL, NULL);
    if (status != UMFPACK_OK)
	phgError(1, "error calling umfpack_di_numeric(), status = %d.\n",
					status);
    umfpack_di_free_symbolic(&_t->symbolic) ;
    _t->symbolic = NULL;
    phgInfo(1, "Time for factorize: %lg\n", (double)(phgGetTime(NULL) - time0));

    return 0;
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    int status;
    double time0;

    Assemble(solver);

    time0 = phgGetTime(NULL);
    status = umfpack_di_solve(UMFPACK_A, _t->Ap, _t->Ai, _t->Ad, x->data,
				solver->rhs->data, _t->numeric, NULL, NULL);
    if (status != UMFPACK_OK)
	phgError(1, "error calling umfpack_di_solve(), status = %d.\n",
					status);
    phgInfo(1, "Time for solve: %lg\n", (double)(phgGetTime(NULL) - time0));

    if (destroy) {
	if (_t->symbolic != NULL) {
	    umfpack_di_free_symbolic(&_t->symbolic);
	    _t->symbolic = NULL;
	}
	if (_t->numeric != NULL) {
	    umfpack_di_free_numeric(&_t->numeric);
	    _t->numeric = NULL;
	}
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

OEM_SOLVER phgSolverSSPARSE_ = {
    "SuiteSparse", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy,
    AddMatrixEntries, AddRHSEntries,
    Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, FALSE, FALSE, FALSE, FALSE
};

#endif	/* USE_SSPARSE */
