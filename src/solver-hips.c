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

/* $Id: solver-hips.c,v 1.23 2022/09/21 02:35:28 zlb Exp $ */

/* Interface to HIPS (http://hips.gforge.inria.fr/) */

#include "phg.h"
#include <limits.h>	/* INT_MAX */

#if USE_HIPS		/* whole file */

#define Initialize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

#include <hips.h>

#include <stdlib.h>
#include <string.h>

static BOOLEAN initialized = FALSE;

/* parameters with cmdline options */

static const char *verbose_keywords[] = {"0", "1", "2", "3", "4", NULL};
static int verbose = 0;

static const char *strategy_keywords[] = {
    "direct", "ilut", "iterative", "hybrid", "block", NULL};
static const INTS strategy_values[] = {
    HIPS_DIRECT, HIPS_ILUT, HIPS_ITERATIVE, HIPS_HYBRID, HIPS_BLOCK};
static int strategy = 2;

/* Options for HIPS parameters (< 0 means unchanged) */
static INT dof = -1;		/* HIPS_DOF */
static FLOAT prec = -1.0;	/* HIPS_PREC */
static FLOAT droptol0 = -1.0;	/* HIPS_DROPTOL0 */
static FLOAT droptol1 = -1.0;	/* HIPS_DROPTOL1 */
static FLOAT dropschur = -1.0;	/* HIPS_DROPSCHUR */
static FLOAT droptole = -1.0;	/* HIPS_DROPTOLE */
static FLOAT amalg = -1.0;	/* HIPS_AMALG */

static BOOLEAN overlap = FALSE;		/* whether use overlap */

typedef struct {
    INTS	id, verbose, strategy, dof;
    INTS	*nodelist;
    COEF	*sol, *rhs;
    REAL	prec, droptol0, droptol1, dropschur, droptole, amalg;
    BOOLEAN	overlap, factored;
} OEM_DATA;

#define	_t		((OEM_DATA *)solver->oem_data)
#define CheckErr	HIPS_ExitOnError(ierr)

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nThe HIPS solver options:", "\n", "hips");
    phgOptionsRegisterKeyword("-hips_verbose", "HIPS verbosity (0-4)",
			verbose_keywords, &verbose);
    phgOptionsRegisterKeyword("-hips_strategy", "HIPS strategy",
			strategy_keywords, &strategy);
    phgOptionsRegisterNoArg("-hips_overlap", "Overlap", &overlap);
    phgOptionsRegisterInt("-hips_dof", "HIPS_DOF", &dof);
    phgOptionsRegisterFloat("-hips_prec", "HIPS_PREC", &prec);
    phgOptionsRegisterFloat("-hips_droptol0", "HIPS_DROPTOL0", &droptol0);
    phgOptionsRegisterFloat("-hips_droptol1", "HIPS_DROPTOL1", &droptol1);
    phgOptionsRegisterFloat("-hips_dropschur", "HIPS_DROPSCHUR", &dropschur);
    phgOptionsRegisterFloat("-hips_droptole", "HIPS_DROPTOLE", &droptole);
    phgOptionsRegisterFloat("-hips_amalg", "HIPS_AMALG", &amalg);

    return 0;
}

static int
Finalize(void)
{
    INTS ierr;

    if (initialized) {
	ierr = HIPS_Finalize(); CheckErr;
	initialized = FALSE;
    }

    return 0;
}

static int
Init(SOLVER *solver)
{
    INTS ierr;

    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));
    _t->id = 0;
    _t->factored = FALSE;
    _t->verbose = verbose;
    _t->strategy = strategy_values[strategy];
#define T(v)	_t->v = v
    T(overlap);
    T(dof);
    T(prec);
    T(droptol0);
    T(droptol1);
    T(dropschur);
    T(droptole);
    T(amalg);
#undef T

    if (!initialized) {
	ierr = HIPS_Initialize(1); CheckErr;
	initialized = TRUE;
    }

    return 0;
}

static int
Create(SOLVER *solver)
{
    INTS ierr;

    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }

    ierr = HIPS_SetDefaultOptions(_t->id, _t->strategy); CheckErr;
    ierr = HIPS_SetOptionINT(_t->id, HIPS_VERBOSE, _t->verbose); CheckErr;

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    INTS ierr;

    if (solver->oem_data == NULL)
	return 0;

    /* FIXME: how to clean up a given HIPS solver instance? */
    /*ierr = HIPS_FreePrecond(_t->id); CheckErr;
    ierr = HIPS_MatrixReset(_t->id); CheckErr;*/
    ierr = HIPS_Clean(_t->id); CheckErr;

    phgFree(_t->nodelist);
    phgFree(_t->sol);
    phgFree(_t->rhs);
    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    MAT *mat = solver->mat;
    MAP *map = mat->rmap;
    INT i, j, m, n0 = map->partition[map->rank];
    INT n = map->nlocal, N = map->nglobal;
    const MAT_ROW *row;

    INTS ierr;
    INTL nnz, nnz0 = solver->mat->nnz_d + solver->mat->nnz_o;
    INTL *ia;
    INTS *ja;
    COEF *la;
    INTS *mapptr = NULL, *mapp = NULL;

    if (_t->factored)
	return 0;

    _t->factored = TRUE;

    if (solver->sent_to_oem == TRUE)
	phgError(1, "%s:%d, data already sent to OEM solver!\n",
						__FILE__, __LINE__);
    solver->sent_to_oem = TRUE;

    assert(solver->mat->type==PHG_PACKED || solver->mat->type==PHG_UNPACKED);

    phgInfo(2, "send data to HIPS solver.\n");

    /* allocate index array */
    if (solver->symmetry >= M_SYM) {
	nnz = nnz0 - n;
	if (nnz < 0)
	    nnz = 0;
	nnz = n + (nnz + 1) / 2;
    }
    else {
	nnz = nnz0;
    }
    if (_t->overlap)
	_t->nodelist = phgAlloc(mat->localsize * sizeof(*_t->nodelist));
    else
	_t->nodelist = phgAlloc(n * sizeof(*_t->nodelist));
    _t->sol = phgAlloc(n * sizeof(*_t->sol));
    _t->rhs = phgAlloc(n * sizeof(*_t->rhs));
    ia = phgAlloc((n + 1) * sizeof(*ia));
    ja = phgAlloc(nnz * sizeof(*ja));
    la = phgAlloc(nnz * sizeof(*la));

    for (i = 0, nnz = 0; i < n; i++) {
	_t->nodelist[i] = n0 + i;
	row = phgMatGetRow(solver->mat, i);
	ia[i] = nnz;
	for (m = j = 0; j < row->ncols; j++) {
	    if (solver->symmetry >= M_SYM && row->cols[j] < n0 + i)
		continue;
	    ja[nnz + m] = row->cols[j];
	    la[nnz + m] = row->data[j];
	    m++;
	}
	nnz += m;
	if (solver->mat->refcount == 0 && solver->mat->rows != NULL) {
	    MAT_ROW *r = solver->mat->rows + i;
	    /* free matrix row */
	    phgFree(r->cols);
	    phgFree(r->data);
	    r->cols = NULL;
	    r->data = NULL;
	    r->ncols = r->alloc = 0;
	}
    }
    ia[n] = nnz;
    assert(nnz == nnz0 || solver->symmetry >= M_SYM);

    if (solver->mat->refcount == 0)
	phgMatFreeMatrix(solver->mat);

    /* Set parameters */
    ierr = HIPS_SetCommunicator(_t->id, map->comm); CheckErr;

#define SET(value, parm, func)	\
	if (_t->value >= 0) {ierr = func(_t->id, parm, _t->value); CheckErr;}
    SET(dof,		HIPS_DOF,	HIPS_SetOptionINT);
    SET(prec,		HIPS_PREC,	HIPS_SetOptionREAL);
    SET(droptol0,	HIPS_DROPTOL0,	HIPS_SetOptionREAL);
    SET(droptol1,	HIPS_DROPTOL1,	HIPS_SetOptionREAL);
    SET(dropschur,	HIPS_DROPSCHUR,	HIPS_SetOptionREAL);
    SET(droptole,	HIPS_DROPTOLE,	HIPS_SetOptionREAL);
    SET(amalg,		HIPS_AMALG,	HIPS_SetOptionREAL);
#undef SET

    /* Enter graph */
    ierr = HIPS_SetOptionINT(_t->id, HIPS_SYMMETRIC,
					solver->symmetry < M_SYM ? 0 : 1);
    CheckErr;
    ierr = HIPS_SetOptionINT(_t->id, HIPS_FORTRAN_NUMBERING, 0); CheckErr;

    /*ierr = HIPS_SetOptionINT(_t->id, HIPS_DOMSIZE, n); CheckErr;*/
    ierr = HIPS_SetOptionINT(_t->id, HIPS_DOMNBR, map->nprocs); CheckErr;

#if 0
    /* This code requires calling Initialize/Finalize HIPS for each solver */
    ierr = HIPS_GraphDistrCSR(_t->id, N, n, _t->nodelist, ia, ja); CheckErr; 
#else
    /* Note: this code is faster than using HIPS_GraphDistrCSR(), and can
     *	     use HIPS_MatrixReset() instead of Initialize/Finalize */
    ierr = HIPS_GraphBegin(_t->id, N, nnz); CheckErr;
    for (i = 0; i < n; i++) {
	INTL j;
	for (j = ia[i]; j < ia[i+1]; j++) {
	    ierr = HIPS_GraphEdge(_t->id, n0 + i, ja[j]); CheckErr;
	}
    }
    ierr = HIPS_GraphEnd(_t->id); CheckErr;
#endif

    ierr = HIPS_SetOptionINT(_t->id, HIPS_GRAPH_SYM,
				solver->symmetry < M_SYM ? 0 : 1);
    CheckErr;

    /* Setup partition */
    if (_t->overlap) {
#if 0
	assert(map->nprocs == 1);	/* requires the whole graph on proc 0 */
	ierr = HIPS_GraphPartition(map->nprocs, _t->overlap, 0,
				   n, ia, ja, 0, &mapptr, &mapp); CheckErr;
#else
	MPI_Datatype type;
	INTS m;
	int size, *cnts = NULL, *dsps = NULL;
	if (map->rank == 0) {
	    cnts = phgAlloc(2 * map->nprocs * sizeof(*cnts));
	    dsps = cnts + map->nprocs;
	}
	assert(mat->localsize <= INT_MAX);
	size = (int)(mat->localsize - map->nlocal);
	MPI_Gather(&size, 1, MPI_INT, cnts, 1, MPI_INT, 0, map->comm);
	if (map->rank == 0) {
	    mapptr = phgAlloc((map->nprocs + 1) * sizeof(*mapptr));
	    for (m = 0, i = 0; i < map->nprocs; i++) {
		assert(m <= INT_MAX);
		dsps[i] = m;
		mapptr[i] = m + map->partition[i];
		m += cnts[i];
	    }
	    mapptr[i] = m + map->partition[i];
	    mapp = phgAlloc((N + m) * sizeof(*mapp));
	}
	for (i = 0; i < size; i++)
	    _t->nodelist[n + i] = solver->mat->O2Gmap[i];
	MPI_Type_contiguous(sizeof(INTS), MPI_BYTE, &type);
	MPI_Type_commit(&type);
	MPI_Gatherv(_t->nodelist + n, size, type, mapp, cnts, dsps, type,
			0, map->comm);
	MPI_Type_free(&type);
	if (map->rank == 0) {
	    for (i = map->nprocs - 1; i >= 0; i--) {
		m = map->partition[i + 1] - map->partition[i];
		/* move off-proc indices to place */
		if (cnts[i] > 0) {
		    memcpy(mapp + mapptr[i] + m, mapp + dsps[i],
					cnts[i] * sizeof(*mapp));
		}
		/* fill in-proc indices */
		for (j = 0; j < m; j++)
		    mapp[mapptr[i] + j] = map->partition[i] + j;
	    }
	    phgFree(cnts);
	}
#endif
    }
    else if (map->rank == 0) {	/* _t->overlap */
	mapptr = phgAlloc((map->nprocs + 1) * sizeof(*mapptr));
	mapp = phgAlloc(N * sizeof(*mapp));
	for (i = 0; i < map->nprocs + 1; i++)
		mapptr[i] = map->partition[i];
	for (i = 0; i < N; i++)
		mapp[i] = i;
    }

    if  (map->rank == 0) {
	ierr = HIPS_SetPartition(_t->id, map->nprocs, mapptr, mapp); CheckErr;
	phgFree(mapptr);
	phgFree(mapp);
    }

    /* Construct the local CSR matrix  */
#warning FIXME: need extra CSR data?
    ierr = HIPS_MatrixDistrCSR(_t->id, n, _t->nodelist, ia, ja, la,
		HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_FOOL,
		solver->symmetry < M_SYM ? 0 : 1); CheckErr;

    phgFree(ia);
    phgFree(ja);
    phgFree(la);

    return 0;
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    INT i, n = solver->mat->rmap->nlocal;
    double time0;

    INTS ierr;
    INTL ival;
    REAL rval;

    Assemble(solver);

    /* Enter RHS */
#warning FIXME: fetch off-proc RHS data
    for (i = 0; i < n; i++)
	_t->rhs[i] = solver->rhs->data[i];

    ierr = HIPS_SetOptionINT(_t->id, HIPS_ITMAX, solver->maxit); CheckErr;
    ierr = HIPS_SetOptionREAL(_t->id, HIPS_PREC, solver->rtol); CheckErr;
    ierr = HIPS_RHSReset(_t->id); CheckErr;

    time0 = phgGetTime(NULL);
    ierr = HIPS_SetRHS(_t->id, n, _t->nodelist, _t->rhs, HIPS_ASSEMBLY_OVW,
		HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_FOOL); CheckErr;
    ierr = HIPS_GetSolution(_t->id, n, _t->nodelist, _t->sol,
		HIPS_ASSEMBLY_FOOL); CheckErr;
    phgInfo(1, "Time for solve: %lg\n", (double)(phgGetTime(NULL) - time0));

    /* Copy solution */
    for (i = 0; i < n; i++)
	x->data[i] = _t->sol[i];

    ierr = HIPS_GetInfoINT(_t->id, HIPS_INFO_ITER, &ival); CheckErr;
    ierr = HIPS_GetInfoREAL(_t->id, HIPS_INFO_RES_NORM, &rval); CheckErr;
    solver->nits = (int)ival;
    solver->residual = (FLOAT)rval;

    if (destroy) {
	ierr = HIPS_MatrixReset(_t->id); CheckErr;
	phgFree(_t->nodelist);
	_t->nodelist = NULL;
	phgFree(_t->sol);
	_t->sol = NULL;
	phgFree(_t->rhs);
	_t->rhs = NULL;
    }

    return solver->nits;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverHIPS_ = {
    "HIPS", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy,
    AddMatrixEntries, AddRHSEntries,
    Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, FALSE, FALSE
};

#endif	/* USE_HIPS */
