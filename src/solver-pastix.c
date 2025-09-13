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

/* $Id: solver-pastix.c,v 1.30 2022/09/21 02:35:28 zlb Exp $
 *
 * The PaStiX interface (https://gforge.inria.fr/projects/pastix/)
 */

#include "phg.h"

/* Note: USE_PASTIX implies USE_MPI (no need to test USE_MPI in this file) */

#if USE_PASTIX	/* whole file */

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

#include <pastix.h>

/* undefine conflicting macros defined in "pastix.h" */
#undef INT	/* pastix_int_t */
#undef FLOAT	/* pastix_float_t */

#include <stdlib.h>
#include <string.h>

typedef struct {
    pastix_data_t	*pastix_data;
    pastix_int_t	ncol, *colptr, *loc2glb, *row, *perm, *invp;
    pastix_float_t	*vals, *rhs;
    pastix_int_t	iparm[IPARM_SIZE];
    pastix_float_t	dparm[DPARM_SIZE];
    pastix_int_t	*global_map;
    int			*cnts0, *dsps0;	/* for gathering data with nlocal */
    int			*cnts1, *dsps1; /* for gathering data with _t->ncol */

    BOOLEAN		factored;

    int			ordering;
    int			refine;
} OEM_DATA;

#define	_t	((OEM_DATA *)solver->oem_data)
#define call_dpastix(solver)						\
    dpastix(&_t->pastix_data, map->comm, _t->ncol, _t->colptr, _t->row,	\
	    _t->vals, _t->loc2glb, _t->perm, _t->invp, _t->rhs,	1,	\
	    _t->iparm, _t->dparm);

/*enum {MAT_SYM = 0, MAT_SSYM = 1, MAT_UNSYM = 2};*/
#define IsSymmetric     (solver->symmetry >= M_SYM)

enum {ORDER_AUTO = 0, ORDER_SCOTCH = 1, ORDER_METIS = 2, ORDER_NONE = 3};
static const char *ordering_keywords[] = {
    "auto", "scotch", "metis", "none", NULL};
static int ordering = 0;	/* ordering scheme */

static INT refine = 0;		/* max nbr of iterative refinements */

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nPaStiX options:", "\n", "pastix");
    phgOptionsRegisterKeyword("-pastix_ordering", "Ordering scheme",
				ordering_keywords, &ordering);
    phgOptionsRegisterInt("-pastix_refine", "Max number of iterative "
				"refinements", &refine);
    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));

    _t->ordering	= ordering;
    _t->refine		= refine;

    return 0;
}

static int
Create(SOLVER *solver)
{
    _t->factored = FALSE;
    phgInfo(2, "create PaStiX solver.\n");

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    MAP *map = solver->mat->rmap;

    if (solver->oem_data == NULL)
	return 0;

    phgInfo(2, "destroy PaStiX solver.\n");

    if (_t->pastix_data != NULL) {
	_t->iparm[IPARM_START_TASK]=API_TASK_CLEAN;
	_t->iparm[IPARM_END_TASK]=API_TASK_CLEAN;
	call_dpastix(solver);
    }

    phgFree(_t->colptr);
    phgFree(_t->loc2glb);
    phgFree(_t->row);
    phgFree(_t->perm);
    phgFree(_t->invp);
    phgFree(_t->vals);
    phgFree(_t->rhs);

    phgFree(_t->global_map);
    phgFree(_t->cnts0);
    phgFree(_t->dsps0);
    phgFree(_t->cnts1);
    phgFree(_t->dsps1);

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    MAP *map = solver->mat->rmap;
    int i, j, k, nz;
    int Nlocal = map->nlocal, n0 = map->partition[map->rank];
    const MAT_ROW *row;
    double time0;
    FLOAT a;

    assert(IsSymmetric);	/* need to transpose matrix if unsymmetric */

    if (_t->factored)
	return 0;

    if (solver->sent_to_oem == TRUE)
	phgError(1, "%s:%d, data already sent to OEM solver!\n",
						__FILE__, __LINE__);
    solver->sent_to_oem = TRUE;

    assert(solver->mat->type != PHG_DESTROYED);
    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }

    phgInfo(2, "send data to PaStiX solver.\n");

    /* count # of nonzeros */
    nz = 0;
    if (!IsSymmetric) {
	for (i = 0; i < Nlocal; i++)
	    nz += phgMatGetNnzInRow(solver->mat, i, NULL, NULL);
    }
    else {
	for (i = 0; i < Nlocal; i++) {
	    row = phgMatGetRow(solver->mat, i);
	    for (j = 0; j < row->ncols; j++) {
		if (row->cols[j] <= i + n0)
		    nz++;
	    }
	}
    }

    /* allocate buffers */
#define ALLOC(ptr, count) (ptr) = phgAlloc((count) * sizeof(*(ptr)))
    _t->ncol = Nlocal;
    ALLOC(_t->colptr, Nlocal + 1);
    ALLOC(_t->loc2glb, Nlocal);
    ALLOC(_t->row, nz);
    ALLOC(_t->perm, Nlocal);
    ALLOC(_t->vals, nz);
    ALLOC(_t->rhs, Nlocal);

    /* setup PaStiX matrix */
    _t->colptr[0] = 1;
    nz = 0;
    for (i = 0; i < Nlocal; i++) {
	_t->perm[i] = 1 + i + n0;
	_t->loc2glb[i] = 1 + i + n0;
	row = phgMatGetRow(solver->mat, i);
	for (j = 0; j < row->ncols; j++) {
	    k = row->cols[j];
	    if (IsSymmetric && k > n0 + i)
		continue;
	    a = row->data[j];
	    if (a > DBL_MAX)
		a = DBL_MAX;
	    _t->vals[nz] = a;
	    _t->row[nz++] = 1 + k;
	}
	_t->colptr[i + 1] = 1 + nz;
	/* free matrix row in solver */
	if (solver->mat->refcount == 0 && solver->mat->rows != NULL) {
	    phgFree(solver->mat->rows[i].cols);
	    phgFree(solver->mat->rows[i].data);
	    solver->mat->rows[i].cols = NULL;
	    solver->mat->rows[i].data = NULL;
	    solver->mat->rows[i].ncols = solver->mat->rows[i].alloc = 0;
	}
    }

    if (solver->mat->refcount == 0)
	phgMatFreeMatrix(solver->mat);

    _t->factored = TRUE;

    phgInfo(2, "initialize dpastix().\n");

    _t->iparm[IPARM_MODIFY_PARAMETER] = API_NO;
    call_dpastix(solver);

    phgInfo(2, "analyse and factor the matrix.\n");

    /*_t->iparm[IPARM_THREAD_NBR] = nbthread;*/
    if (IsSymmetric) {
	_t->iparm[IPARM_SYM]		= API_SYM_YES;
	_t->iparm[IPARM_FACTORIZATION]	= API_FACT_LDLT;
    }
    else {
	_t->iparm[IPARM_SYM]		= API_SYM_NO;
	_t->iparm[IPARM_FACTORIZATION]	= API_FACT_LU;
    }
    _t->iparm[IPARM_FREE_CSCUSER]	= API_NO;
    _t->iparm[IPARM_FREE_CSCPASTIX]	= API_YES;
    _t->iparm[IPARM_RHS_MAKING]         = API_RHS_B;
    _t->iparm[IPARM_MATRIX_VERIFICATION]= API_NO;
    if (!solver->monitor)
	_t->iparm[IPARM_VERBOSE]	= API_VERBOSE_NOT;	/* no mesg */
    else
	_t->iparm[IPARM_VERBOSE]	= API_VERBOSE_YES;

    _t->iparm[IPARM_GRAPHDIST]		= API_YES;
    if (_t->ordering == ORDER_SCOTCH) {
	_t->iparm[IPARM_ORDERING]	= API_ORDER_SCOTCH;
    }
    else if (_t->ordering == ORDER_METIS) {
	_t->iparm[IPARM_GRAPHDIST]	= API_NO; /* non distr. graph only */
	_t->iparm[IPARM_LEVEL_OF_FILL]	= -1;
	_t->iparm[IPARM_ORDERING]	= API_ORDER_METIS;
    }
    else if (_t->ordering == ORDER_NONE) {
	_t->iparm[IPARM_LEVEL_OF_FILL]	= -1;
	_t->iparm[IPARM_ORDERING]	= API_ORDER_PERSONAL;
    }
    _t->iparm[IPARM_START_TASK]		= API_TASK_ORDERING;
    _t->iparm[IPARM_END_TASK]		= API_TASK_ORDERING;
    time0 = phgGetTime(NULL);
    call_dpastix(solver);
    if (_t->iparm[IPARM_ERROR_NUMBER] != 0 && map->rank == 0)
	phgWarning("PaStiX returns a non-zero code in ORDERING step: %d\n",
			_t->iparm[IPARM_ERROR_NUMBER]);
    phgInfo(1, "Time for ORDERING: %lg\n", phgGetTime(NULL) - time0);

    _t->iparm[IPARM_START_TASK]		= API_TASK_SYMBFACT;
    _t->iparm[IPARM_END_TASK]		= API_TASK_SYMBFACT;
    time0 = phgGetTime(NULL);
    call_dpastix(solver);
    if (_t->iparm[IPARM_ERROR_NUMBER] != 0 && map->rank == 0)
	phgWarning("PaStiX returns a non-zero code in SYMBFACT step: %d\n",
			_t->iparm[IPARM_ERROR_NUMBER]);
    phgInfo(1, "Time for SYMBFACT: %lg\n", phgGetTime(NULL) - time0);

    _t->iparm[IPARM_START_TASK]		= API_TASK_ANALYSE;
    _t->iparm[IPARM_END_TASK]		= API_TASK_ANALYSE;
    time0 = phgGetTime(NULL);
    call_dpastix(solver);
    if (_t->iparm[IPARM_ERROR_NUMBER] != 0 && map->rank == 0)
	phgWarning("PaStiX returns a non-zero code in ANALYSE step: %d\n",
			_t->iparm[IPARM_ERROR_NUMBER]);
    phgInfo(1, "Time for ANALYSE: %lg\n", phgGetTime(NULL) - time0);

    if (map->nprocs > 1) {
	/* redispatch the matrix according to new distribution */
	pastix_int_t *colptr0, *loc2glb0, *row0;
	pastix_float_t *vals0, *rhs0;
	extern int cscd_redispatch(
	    pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja,
	    pastix_float_t *a, pastix_float_t *rhs, pastix_int_t *l2g,
	    pastix_int_t dn, pastix_int_t **dia, pastix_int_t **dja,
	    pastix_float_t **da, pastix_float_t **drhs, pastix_int_t *dl2g,
	    MPI_Comm comm
	);
	colptr0 = _t->colptr;
	row0 = _t->row;
	vals0 = _t->vals;
	rhs0 = _t->rhs;
	loc2glb0 = _t->loc2glb;
	_t->ncol = pastix_getLocalNodeNbr(&_t->pastix_data);
	phgInfo(1, "redispatch matrix: old nlocal = %d, new nlocal = %d\n",
			Nlocal, _t->ncol);
	_t->loc2glb = phgAlloc(_t->ncol * sizeof(*_t->loc2glb));
	pastix_getLocalNodeLst(&_t->pastix_data, _t->loc2glb);
	time0 = phgGetTime(NULL);
	cscd_redispatch(Nlocal, colptr0, row0, vals0, rhs0, loc2glb0,
			_t->ncol, &_t->colptr, &_t->row, &_t->vals, &_t->rhs,
			_t->loc2glb, map->comm);
	phgInfo(1, "Time for redispatching: %lg\n", phgGetTime(NULL) - time0);
	phgFree(colptr0);
	phgFree(row0);
	phgFree(vals0);
	phgFree(rhs0);
	phgFree(loc2glb0);

	ALLOC(_t->global_map, map->nglobal);
	ALLOC(_t->cnts0, map->nprocs);
	ALLOC(_t->dsps0, map->nprocs);
	ALLOC(_t->cnts1, map->nprocs);
	ALLOC(_t->dsps1, map->nprocs);

	MPI_Allgather(&_t->ncol, 1, MPI_PASTIX_INT,
		      _t->cnts1, 1, MPI_PASTIX_INT, map->comm);
	k = 0;
	for (i = 0; i < map->nprocs; i++) {
	    _t->cnts0[i] = map->partition[i + 1] - map->partition[i];
	    _t->dsps0[i] = map->partition[i];
	    _t->dsps1[i] = k;
	    k += _t->cnts1[i];
	}
	assert(k == map->nglobal);

	MPI_Allgatherv(_t->loc2glb, _t->ncol, MPI_PASTIX_INT,
		      	_t->global_map, _t->cnts1, _t->dsps1, MPI_PASTIX_INT,
			map->comm);
    }

#define FREE(ptr) phgFree(ptr); ptr = NULL

    /* perm[] is no more needed */
    FREE(_t->perm);

    /* FIXME: bug in PaStiX 5.1.2: NUMFACT step includes the SOLVE step.
     * 
     * To test: commenting out the next block and change 'SOLVE' to 'NUMFACT'
     *		in Solve(), "./poisson -solver pastix" still runs correctly! */
    time0 = phgGetTime(NULL);
    _t->iparm[IPARM_RHS_MAKING]	= API_RHS_B;
    if (_t->ncol > 0)
	bzero(_t->rhs, _t->ncol * sizeof(*_t->rhs));
    _t->iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
    _t->iparm[IPARM_END_TASK] = API_TASK_NUMFACT;
    call_dpastix(solver);
    /* check return code */
    if (_t->iparm[IPARM_ERROR_NUMBER] != 0 && map->rank == 0)
	phgWarning("PaStiX returns a non-zero code in NUMFACT step: %d\n",
			_t->iparm[IPARM_ERROR_NUMBER]);
    phgInfo(1, "Time for NUMFACT: %lg\n", phgGetTime(NULL) - time0);

    return 0;
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    MAP *map = solver->mat->rmap;
    INT i, k, Nlocal = map->nlocal;
    double time0;
    FLOAT a;
    pastix_float_t *buf0, *buf1;

    Assemble(solver);

    if (map->nprocs == 1) {
	buf0 = _t->rhs;
	buf1 = NULL;
    }
    else {
	buf0 = phgAlloc(Nlocal * sizeof(*buf0));
	buf1 = phgAlloc(map->nglobal * sizeof(*buf1));
    }

    for (i = 0; i < Nlocal; i++) {
	a = solver->rhs->data[i];
	if (a > DBL_MAX)
	    a = DBL_MAX;
	buf0[i] = (pastix_float_t)a;
    }

    if (map->nprocs > 1) {
	/* redispatch RHS (FIXME: better to use MPI_Alltoallv) */
	MPI_Allgatherv(buf0, Nlocal, MPI_PASTIX_FLOAT,
		       buf1, _t->cnts0, _t->dsps0, MPI_PASTIX_FLOAT, map->comm);
	for (i = 0; i < _t->ncol; i++)
	    _t->rhs[i] = buf1[_t->loc2glb[i] - 1];
    }

    time0 = phgGetTime(NULL);
    _t->iparm[IPARM_RHS_MAKING]	= API_RHS_B;
    _t->iparm[IPARM_START_TASK]	= API_TASK_SOLVE;
    _t->iparm[IPARM_END_TASK]	= API_TASK_SOLVE;
    call_dpastix(solver);
    if (_t->iparm[IPARM_ERROR_NUMBER] != 0 && map->rank == 0)
	phgWarning("PaStiX returns a non-zero code in SOLVE step: %d\n",
			_t->iparm[IPARM_ERROR_NUMBER]);
    phgInfo(1, "Time for SOLVE: %lg\n", phgGetTime(NULL) - time0);

    if (map->nprocs > 1) {
	/* redispatch solution (FIXME: better to use MPI_Alltoallv) */
	INT n0 = map->partition[map->rank], n1 = map->partition[map->rank + 1];
	MPI_Allgatherv(_t->rhs, _t->ncol, MPI_PASTIX_FLOAT,
		       buf1, _t->cnts1, _t->dsps1, MPI_PASTIX_FLOAT, map->comm);
	for (i = 0; i < map->nglobal; i++) {
	    k = _t->global_map[i] - 1;
	    if (k < n0 || k >= n1)
		continue;
	    x->data[k - n0] = buf1[i];
	}

	phgFree(buf0);
	phgFree(buf1);
    }
    else {
	for (i = 0; i < Nlocal; i++)
	    x->data[i] = buf0[i];
    }

    if (destroy) {
	if (_t->pastix_data != NULL) {
	    _t->iparm[IPARM_START_TASK]=API_TASK_CLEAN;
	    _t->iparm[IPARM_END_TASK]=API_TASK_CLEAN;
	    call_dpastix(solver);
	    _t->pastix_data = NULL;
	}
	FREE(_t->colptr);
	FREE(_t->loc2glb);
	FREE(_t->row);
	FREE(_t->perm);
	FREE(_t->invp);
	FREE(_t->vals);
	FREE(_t->rhs);

	FREE(_t->global_map);
	FREE(_t->cnts0);
	FREE(_t->dsps0);
	FREE(_t->cnts1);
	FREE(_t->dsps1);
    }

    solver->residual = 0.0;	/* fake residual norm */

    return solver->nits = 1;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverPaStiX_ = {
    "PaStiX", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy,
    AddMatrixEntries, AddRHSEntries,
    Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, FALSE, TRUE, FALSE, FALSE
};

#endif	/* USE_PASTIX */
