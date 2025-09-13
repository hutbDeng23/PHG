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

/* $Id: solver-superlu.c,v 1.63 2022/09/21 02:35:28 zlb Exp $ */

/* The SuperLU interface for PHG */

/* This controls fill-reducing code using ParMETIS (about 10% improvement when
   enabled, currently only implemented for uniprocess case) */
#define FILL_REDUCING 0

#include "phg.h"

#if USE_SUPERLU

#if !USE_PARMETIS
# undef FILL_REDUCING
# define FILL_REDUCING 0
#endif

/* Hack: 'FLOAT' in util_dist.h conflicts with PHG's type */
#define FLOAT	SUPERLU_FLOAT
#include <superlu_ddefs.h>
#undef FLOAT

#if defined(SUPERLU_DIST_MAJOR_VERSION) && SUPERLU_DIST_MAJOR_VERSION >= 6
# if SUPERLU_DIST_MAJOR_VERSION > 6 || SUPERLU_DIST_MINOR_VERSION > 2
# define ScalePermstruct_t	dScalePermstruct_t
# define LUstruct_t		dLUstruct_t
# define SOLVEstruct_t		dSOLVEstruct_t
# define ScalePermstructFree	dScalePermstructFree
# define Destroy_LU		dDestroy_LU
# define LUstructFree		dLUstructFree
# define ScalePermstructInit	dScalePermstructInit
# define LUstructInit		dLUstructInit
# endif
#endif

#if SUPERLU_VERSION_MAJOR >= 5
    /* 'superlu_options_t' is renamed to 'superlu_dist_options_t' in v5.0 */
    typedef superlu_dist_options_t superlu_options_t;
#endif

#if FILL_REDUCING
# include <stdlib.h>
# include <math.h>
# include <metis.h>
# include <parmetis.h>
#ifndef PARMETIS_MAJOR_VERSION
# define PARMETIS_MAJOR_VERSION 3
#endif
# if PARMETIS_MAJOR_VERSION <= 3
typedef idxtype idx_t;
# endif	/* PARMETIS_MAJOR_VERSION <= 3 */
#endif

typedef struct {
    double	*B;		/* B=RHS and solution */
    /* SuperLU matrix data */
    SuperMatrix	*A;		/* A=matrix */
    int		nnz, firstrow;	/* # of local nonzeros and first row index */
    int		*rowptr;	/* indices to array data and cols
				   (beginning of rows) [nprocs + 1] */
    int		*cols;		/* array of column indices [nnz] */
    double	*data;		/* nonzero values [nnz] */
#if FILL_REDUCING
    idx_t	*iorder, *order;
#endif

    /* SuperLU data */
    BOOLEAN		superlu_has_data;
    superlu_options_t	options;
    SuperLUStat_t	stat;
    ScalePermstruct_t	ScalePermstruct;
    LUstruct_t		LUstruct;
    SOLVEstruct_t	SOLVEstruct;
    gridinfo_t		grid;

    int		refine_index;
    BOOLEAN	print_stat;
    BOOLEAN	par_symb_fact;
} OEM_DATA;

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)

#define Initialize	NULL
#define Finalize	NULL
#define SetPC		NULL

/* SuperLU options (more options will be added when needed) */

/* Note: must match the definition of IterRefine_t in superlu_defs.h (2.x) or
 * superlu_enum_consts.h (3.x) */
enum {slu_norefine, slu_single=1, slu_double, slu_extra};
static int refine_index	  = 1;
static const char *refine_names[] = {"norefine", "double"/*, "extra"*/, NULL};
static IterRefine_t refine_types[] = {slu_norefine, slu_double/*, slu_extra*/};

static BOOLEAN print_stat = FALSE;
/* FIXME: the program "./poisson -solver superlu" crashes if parallel symbolic
 * factorization is enabled (and it currently requires nprocs to be 2^k) */
static BOOLEAN par_symb_fact = FALSE;

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nSuperLU_DIST options:", "\n", "superlu");
    phgOptionsRegisterKeyword("superlu_iter_refine", "Iterative refinement",
		refine_names, &refine_index);
    phgOptionsRegisterNoArg("superlu_print_stat", "Print statistics",
		&print_stat);
    phgOptionsRegisterNoArg("superlu_par_symb_fact",
		"Perform parallel symbolic factorization", &par_symb_fact);

    return 0;
}

#if FILL_REDUCING

static void
fill_reducing_ordering(SOLVER *solver)
/* calls ParMETIS to redistribute the grid into nparts processes.
   The result of the partitioning is stored in the 'mark' field, i.e.,
   e->mark is the rank of the process the element is assigned to. */
{
    MAP *map = solver->mat->rmap;
    idx_t *vtxdist;
    idx_t *xadj, *adjncy;
    idx_t *sizes;
    idx_t numflag = 0, n;
    idx_t options[3] = {1, /*verbosity=*/ 0, /*seed=*/15};
    int i, j, k, l, index0;
    const MAT_ROW *row;
    INT nlocal = map->nlocal;

    /* Note: ParMETIS requires map->nprocs to be power of 2!!! */
#ifdef log2
# undef log2
#endif
#define log2(x) (log(x)/log(2.0))
    assert(solver->mat != NULL);
    assert(map->nprocs == (1 << round(log2((double)map->nprocs))));

    /* the fill-reducing code currently only works with map->nprocs == 1 */
    if (map->nprocs != 1)
	return;
 
    vtxdist = phgAlloc((map->nprocs + 1) * sizeof(*vtxdist));
    for (i = 0; i <= map->nprocs; i++)
	vtxdist[i] = map->partition[i];
    index0 = vtxdist[map->rank];

    /* compute size of the adjacency array (non zero off-diag entries) */
    n = 0;	/* used to count # of edges in the subgraph */
    for (i = 0; i < nlocal; i++) {
	row = phgMatGetRow(solver->mat, i);
	for (j = 0; j < row->ncols; j++) {
	    if (row->cols[j] != i + index0)
		n++;
	}
    }

    xadj    = phgAlloc((nlocal + 1) * sizeof(*xadj));
    adjncy  = phgAlloc(n * sizeof(*adjncy));
    _t->order   = phgAlloc(nlocal * sizeof(*_t->order));
    _t->iorder  = phgAlloc(nlocal * sizeof(*_t->iorder));
    sizes   = phgAlloc(2 * map->nprocs * sizeof(*sizes));

    /* set up the distributed CSR arrays */
    k = 0;
    for (i = 0; i < nlocal; i++) {
	xadj[i] = k;
	row = phgMatGetRow(solver->mat, i);
	for (j = 0; j < row->ncols; j++) {
	    if ((l = row->cols[j]) != i + index0) {
		adjncy[k++] = l;
	    }
	}
    }
    xadj[i] = k;

#if DEBUG
    for (i = 0; i < nlocal; i++) {
	for (j = xadj[i]; j < xadj[i + 1]; j++) {
	    n = adjncy[j] - index0;
	    for (k = xadj[n]; k < xadj[n + 1]; k++)
		if (adjncy[k] == n + index0)
		    break;
	    /*if (k >= xadj[n + 1])
		phgWarning("A(%d,%d) is present, but A(%d,%d) is not.\n",
			i + index0, n + index0, n + index0, i + index0);*/
	}
    }
#endif	/* DEBUG */

    if (map->nprocs == 1) {
	options[0] = 0;
	n = nlocal;
	METIS_NodeND(&n, xadj, adjncy, &numflag, options, _t->order,_t->iorder);
    }
    else {
	MPI_Comm comm = map->comm;
	ParMETIS_V3_NodeND(vtxdist, xadj, adjncy, &numflag, options,
				_t->iorder, sizes, &comm);
	for (i = 0; i < nlocal; i++) {
	    _t->order[_t->iorder[i]] = i;
	}
    }

#if DEBUG
    /* print the ordering computed */
    for (i = 0; i < nlocal; i++) {
	/*phgInfo(-1, "    %d ==> %d\n", index0 + i, order[i]);*/
	if (map->nprocs == 1) {
	    assert((k = _t->order[i] - index0) >= 0 && k < nlocal);
	    assert((k = _t->iorder[i] - index0) >= 0 && k < nlocal);
	    assert(_t->iorder[_t->order[i]] == i
			&& _t->order[_t->iorder[i]] == i);
	}
    }
#endif	/* DEBUG */

    phgFree(vtxdist);
    phgFree(xadj);
    phgFree(adjncy);
    phgFree(sizes);
}
#endif	/* FILL_REDUCING */

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgAlloc(sizeof(OEM_DATA));

    _t->refine_index	= refine_index;
    _t->print_stat	= print_stat;
    _t->par_symb_fact	= par_symb_fact;

    return 0;
}

static int
Create(SOLVER *solver)
{
    MAP *map = solver->mat->rmap;
    int i, j, Nlocal = map->nlocal;

    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }
    assert(solver->mat->assembled);	/* only works with assembled matrix */

    _t->superlu_has_data = FALSE;
    _t->A = NULL;
    _t->B = phgAlloc(Nlocal * sizeof(*_t->B));
    _t->firstrow = map->partition[map->rank];
    if ((_t->rowptr = intMalloc_dist(Nlocal + 1)) == NULL)
	phgError(1, "%s:%d, failed malloc.\n", __FILE__, __LINE__);

#if FILL_REDUCING
    /* compute fill-reducing ordering by ParMETIS */
    _t->order = NULL;
    _t->iorder = NULL;
    fill_reducing_ordering(solver);
#endif	/* FILL_REDUCING */

    j = 0;
    for (i = 0; i < Nlocal; i++) {
	_t->B[i] = 0.0;
	_t->rowptr[i] = j;
#if FILL_REDUCING
	if (_t->order != NULL)
	    j += phgMatGetNnzInRow(solver->mat, _t->order[i], NULL, NULL);
	else
#endif
	j += phgMatGetNnzInRow(solver->mat, i, NULL, NULL);
    }
    _t->rowptr[Nlocal] = _t->nnz = j;

    /* allocate arrays */
    if ((_t->cols = intMalloc_dist(_t->nnz)) == NULL ||
	(_t->data = doubleMalloc_dist(_t->nnz)) == NULL)
	phgError(1, "%s:%d, failed malloc.\n", __FILE__, __LINE__);

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    if (_t->superlu_has_data) {
	PStatFree(&_t->stat);
	ScalePermstructFree(&_t->ScalePermstruct);
	Destroy_LU(solver->mat->rmap->nglobal, &_t->grid, &_t->LUstruct);
	LUstructFree(&_t->LUstruct);
	if (_t->options.SolveInitialized)
	    dSolveFinalize(&_t->options, &_t->SOLVEstruct);
	superlu_gridexit(&_t->grid);
	_t->superlu_has_data = FALSE;
    }

    if (_t->A != NULL) {
	Destroy_CompRowLoc_Matrix_dist(_t->A);
	phgFree(_t->A);
    }
    else {
	if (_t->rowptr != NULL)
	    SUPERLU_FREE(_t->rowptr);
	if (_t->cols != NULL)
	    SUPERLU_FREE(_t->cols);
	if (_t->data != NULL)
	    SUPERLU_FREE(_t->data);
    }
    phgFree(_t->B);
#if FILL_REDUCING
    phgFree(_t->order);
    phgFree(_t->iorder);
#endif	/* FILL_REDUCING */
    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
AddMatrixEntries(SOLVER *solver, INT nrows, INT *rows, INT ncols, INT *cols,
			FLOAT *data)
{
    INT i;

    /* note: assuming this function is only called once for each row */
#if FILL_REDUCING
    assert(solver->mat->rmap->nprocs == 1);
    for (i = 0; i < nrows; i++) {
	int j, k = _t->rowptr[_t->iorder[rows[i]] - _t->firstrow];
	for (j = 0; j < ncols; j++) {
	    _t->cols[k + j] = _t->iorder[cols[j]];
	    _t->data[k + j] = data[j];
	}
	data += ncols;
    }
#else	/* FILL_REDUCING */
    for (i = 0; i < nrows; i++) {
#if FT_PHG == FT_DOUBLE && IT_PHG == IT_INT
	memcpy(_t->cols + _t->rowptr[rows[i] - _t->firstrow], cols,
		ncols * sizeof(*cols));
	memcpy(_t->data + _t->rowptr[rows[i] - _t->firstrow], data,
		ncols * sizeof(*data));
	data += ncols;
#else	/* FT_PHG == FT_DOUBLE */
	int j, k = _t->rowptr[rows[i] - _t->firstrow];
	for (j = 0; j < ncols; j++) {
	    _t->cols[k + j] = cols[j];
	    _t->data[k + j] = data[j];
	}
	data += ncols;
#endif	/* FT_PHG == FT_DOUBLE */
    }
#endif	/* FILL_REDUCING */

    return 0;
}

static int
AddRHSEntries(SOLVER *solver, INT n, INT *ni, FLOAT *values)
{
    INT i;

#if FILL_REDUCING
    assert(solver->mat->rmap->nprocs == 1);
    for (i = 0; i < n; i++)
	_t->B[_t->iorder[ni[i]] - _t->firstrow] = *(values++);
#else	/* FILL_REDUCING */
    for (i = 0; i < n; i++)
	_t->B[ni[i] - _t->firstrow] = *(values++);
#endif	/* FILL_REDUCING */

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    if (_t->A != NULL)	/* already assembled */
	return 0;

    _t->A = phgAlloc(sizeof(*_t->A));

    dCreate_CompRowLoc_Matrix_dist(_t->A, solver->mat->rmap->nglobal,
						solver->mat->rmap->nglobal,
	_t->nnz, solver->mat->rmap->nlocal, _t->firstrow, _t->data,
						_t->cols, _t->rowptr,
	SLU_NR_loc, SLU_D, SLU_GE);

#if 0
    /* print matrix rows */
    int i, j, ncols, row;
    for (i = 0; i < solver->mat->rmap->nlocal; i++) {
	row = i + _t->firstrow;
	j = _t->rowptr[i];
	ncols = _t->rowptr[i + 1] - j;
	phgInfo(-1, "row %d, ncols = %d, rhs = %lg\n", row, ncols,
			(double)_t->B[i]);
	for (; j < _t->rowptr[i + 1]; j++)
	    phgInfo(-1, "\tA[%d][%d] = %lg\n", row, _t->cols[j],
			(double)_t->data[j]);
    }
#endif

    return 0;
}
 
static void
copy_solution(SOLVER *solver, VEC *x)
/* copy solution from SuperLU's vector 'B' to local DOFs */
{
    double *vec;
    size_t size;
    INT i;

    size = solver->mat->rmap->nlocal;
#if FILL_REDUCING
    assert(solver->mat->rmap->nprocs == 1);
    vec = phgAlloc(size * sizeof(*vec));
    for (i = 0; i < size; i++)
	vec[i] = _t->B[_t->iorder[i]];
#else	/* FILL_REDUCING */
    vec = _t->B;
#endif	/* FILL_REDUCING */

    for (i = 0; i < size; i++)
	x->data[i] = vec[i];

#if FILL_REDUCING
    phgFree(vec);
#endif

    return;
}

static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    double error;
    int info, dims[2] = {0, 0};

    Assemble(solver);

    if (!_t->superlu_has_data) {
	_t->superlu_has_data = TRUE;
	MPI_Dims_create(solver->mat->rmap->nprocs, 2, dims);
	superlu_gridinit(solver->mat->rmap->comm, dims[0], dims[1],
			 &_t->grid);
	ScalePermstructInit(solver->mat->rmap->nglobal,
			    solver->mat->rmap->nglobal,
			    &_t->ScalePermstruct);
#if SUPERLU_VERSION_MAJOR <= 3
	/* SuperLU-dist <= 3.x */
	LUstructInit(solver->mat->rmap->nglobal, solver->mat->rmap->nglobal,
			&_t->LUstruct);
#else
	/* SuperLU-dist >= 4.x */
	LUstructInit(solver->mat->rmap->nglobal, &_t->LUstruct);
#endif
    	PStatInit(&_t->stat);
	/* Set the default input options
	    options.Fact = DOFACT;
	    options.Equil = YES;
	    options.ColPerm = MMD_AT_PLUS_A;
	    options.RowPerm = LargeDiag;
	    options.ReplaceTinyPivot = YES;
	    options.Trans = NOTRANS;
	    options.IterRefine = slu_double;
	    options.SolveInitialized = NO;
	    options.RefineInitialized = NO;
	    options.PrintStat = YES;
	*/
	set_default_options_dist(&_t->options);
#if 0
	_t->options.RowPerm = NOROWPERM;
	_t->options.ColPerm = NATURAL;
	_t->options.Equil = NO; 
	_t->options.ReplaceTinyPivot = NO;
#endif
	_t->options.IterRefine = refine_types[_t->refine_index];
	_t->options.PrintStat = (yes_no_t)_t->print_stat;
	_t->options.Fact = DOFACT;
#if USE_SUPERLU_PARSYMBFACT	/* note: available in SuperLU_Dist >= 2.1 */
	if (_t->par_symb_fact == TRUE) {
	    _t->options.ParSymbFact = YES;
	    _t->options.ColPerm = PARMETIS;
	}
	else {
	    _t->options.ParSymbFact = NO;
	}
#endif	/* USE_SUPERLU_PARSYMBFACT */
    }
    else {
	_t->options.Fact = FACTORED;
    }

    pdgssvx(&_t->options, _t->A, &_t->ScalePermstruct, _t->B,
	    solver->mat->rmap->nlocal,
	    1, &_t->grid, &_t->LUstruct, &_t->SOLVEstruct, &error, &_t->stat,
	    &info);
    if (_t->print_stat)
	PStatPrint(&_t->options, &_t->stat, &_t->grid);

    if (destroy) {
	PStatFree(&_t->stat);
	ScalePermstructFree(&_t->ScalePermstruct);
	Destroy_LU(solver->mat->rmap->nglobal, &_t->grid, &_t->LUstruct);
	LUstructFree(&_t->LUstruct);
	if (_t->options.SolveInitialized)
	    dSolveFinalize(&_t->options, &_t->SOLVEstruct);
	superlu_gridexit(&_t->grid);

	_t->superlu_has_data = FALSE;

	Destroy_CompRowLoc_Matrix_dist(_t->A);
	phgFree(_t->A);
	_t->A = NULL;
#if 0
	/* Note: the arrays are freed by Destroy_CompRowLoc_Matrix_dist() */
	if (_t->rowptr != NULL)
	    SUPERLU_FREE(_t->rowptr);
	if (_t->cols != NULL)
	    SUPERLU_FREE(_t->cols);
	if (_t->data != NULL)
	    SUPERLU_FREE(_t->data);
#endif
	_t->rowptr = NULL;
	_t->cols = NULL;
	_t->data = NULL;
    }

    copy_solution(solver, x);

    solver->residual = /*error*/0.0;
    return solver->nits = 1;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverSuperLU_ = {
    "SuperLU", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, FALSE, TRUE, FALSE, FALSE
};

#endif		/* USE_SUPERLU */
