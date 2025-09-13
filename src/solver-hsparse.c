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

/* $Id: solver-hsparse.c,v 1.13 2022/09/21 02:35:28 zlb Exp $ */

/* Interface for swSparse solver */

#include "phg.h"

#if USE_HSPARSE

#include "swSparse.h"

#if USE_METIS
# include <metis.h>
# ifndef METIS_VER_MAJOR
#  define METIS_VER_MAJOR 4
# endif
# if METIS_VER_MAJOR <= 4
typedef idxtype idx_t;
# endif	/* if METIS_VER_MAJOR <= 4 */
static idx_t metis_options[METIS_NOPTIONS];
#endif	/* USE_METIS */

typedef struct {
    swSparseMatrix	l_trsv, u_trsv;
    INT			levels;		/* k for ILU(k), <0 => full LU */
    int			ordering;	/* ordering type */
    BOOLEAN		assembled;
} OEM_DATA;

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)

#define Initialize		NULL
#define Finalize		NULL
#define SetPC			NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL

static INT levels = 0;

enum {ORDERING_NONE = 0, ORDERING_METIS = 1};
static const char *ordering_keywords[] = { "none", "metis", NULL};
#if USE_METIS
static int ordering = ORDERING_METIS;
#else	/* USE_METIS */
static int ordering = ORDERING_NONE;
#endif	/* USE_METIS */

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nhSparse solver options:", "\n", "hsparse");

    /* number of levels of fill for ILU(k) */
    phgOptionsRegisterInt("-hsparse_levels", "Fill levels (negative->full LU)",
				&levels);
    phgOptionsRegisterKeyword("-hsparse_ordering", "Ordering scheme",
				ordering_keywords, &ordering);

    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));

    _t->levels = levels;
    _t->ordering = ordering;

    return 0;
}

static int
Create(SOLVER *solver)
{
    /* do nothing */
    Unused(solver);
    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    /* TODO: free allocated swSparse objects */
    swSparse_destroy(_t->l_trsv);
    swSparse_destroy(_t->u_trsv);

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

#ifndef SORT_COLS
# define SORT_COLS	0
#endif

static const MAT_ROW *row;

#if SORT_COLS
static int *index = NULL, allocated = 0;

static
int comp(const void *c0, const void *c1)
{
    INT d = row->cols[*((const INT *)c0)] - row->cols[*((const INT *)c1)];
    return d == 0 ? 0 : (d < 0 ? -1 : 1);
}
#endif	/* SORT_COLS */

static int
Assemble(SOLVER *solver)
{
    swSparseInt i, j, k;
    swSparseInt n = solver->mat->rmap->nglobal, nnz = solver->mat->nnz_d;
    swSparseInt *ciperm, *riperm, *row_off, *col_idx;
    double *values;
    swSparseMatrix a, p, l, u;
    BOOLEAN flag, is_symmetric;

    assert(solver->mat->rmap->nprocs == 1 && solver->mat->cmap->nprocs == 1);
    if (_t->assembled)
	return 0;

    _t->assembled = TRUE;

    /* create swSparse matrices */
    row_off = phgAlloc((n + 1) * sizeof(*row_off));
    col_idx = phgAlloc(nnz * sizeof(*col_idx));
    values = phgAlloc(nnz * sizeof(*values));

    k = 0;
    for (i = 0; i < n; i++) {
	row_off[i] = k;
	row = phgMatGetRow(solver->mat, i);
	if (row->ncols == 0)
	    continue;
#if SORT_COLS
	/* sort row->cols[] into increasing order */
	if (allocated < row->ncols) {
	    phgFree(index);
	    index = phgAlloc((allocated = row->ncols + 16) * sizeof(*index));
	}
	for (j = 0; j < row->ncols; j++)
	    index[j] = j;
	qsort(index, row->ncols, sizeof(*index), comp);
	for (j = 0; j < row->ncols; j++, k++) {
	    /* make sure columns are ordered */
	    assert(j == 0 || row->cols[index[j]] > row->cols[index[j - 1]]);
	    col_idx[k] = row->cols[index[j]];
	    values[k] = row->data[index[j]];
	}
#else	/* SORT_COLS */
	for (j = 0; j < row->ncols; j++, k++) {
	    col_idx[k] = row->cols[j];
	    values[k] = row->data[j];
	}
#endif	/* SORT_COLS */
    }
    row_off[n] = k;
    assert(k == nnz);

#if SORT_COLS
    phgFree(index);
    index = NULL;
    allocated = 0;
#endif	/* SORT_COLS */

    /* Create CSR matrix */
    swSparse_d_create_csr(&p, n, n, nnz, row_off, col_idx, values, 0);

    phgFree(values);

    /* compute ordering */
    ciperm = phgAlloc(n * sizeof(*ciperm));
    riperm = phgAlloc(n * sizeof(*riperm));
    flag = FALSE;
    switch (_t->ordering) {
	case ORDERING_NONE:
	    break;
	case ORDERING_METIS:
#if USE_METIS
	    /* Note: METIS_NodeND() requires the graph data to be strictly
	     *	     symmetric and without diagonal entries */
	    assert(sizeof(idx_t) == sizeof(swSparseInt));
#if 1
	    /* check symmetry of the graph */
	    for (i = k = 0; i < n; i++) {
		for (; k < row_off[i + 1]; k++) {
		    swSparseInt l;
		    j = col_idx[k];
		    if (j == i)
			continue;
		    /* check whether col-i exists in row-j */
		    for (l = row_off[j]; l < row_off[j + 1]; l++)
			if (col_idx[l] == i)
			    break;
		    if (l == row_off[j + 1])
			break;	/* not found -> unsymmetric graph */
		}
		if (k < row_off[i + 1])
		    break;
	    }
	    is_symmetric = (k == row_off[n]);
#else	/* 1 */
	    is_symmetric = FALSE;
#endif	/* 1 */
	    /*printf("\nis_symmetric = %d\n", is_symmetric);*/
	    if (is_symmetric) {
	 	/* the graph is symmetric, only remove diagonal entries */
	 	for (i = j = k = 0; i < n; i++) {
		row_off[i] = k;
		for (; j < row_off[i + 1]; j++)
		    if (col_idx[j] != i)
			col_idx[k++] = col_idx[j];
	 	}
	 	row_off[i] = k;
	    }
	    else {	/* is_symmetric */
	 	/* remove diagonal entries and force symmetry */
		swSparseInt c, r;
		INT nnz = 0;
		MAT_ROW *rows = phgCalloc(n, sizeof(*rows)), *row;
		/* pass 1: copy rows */
		for (i = k = 0; i < n; i++) {
		    row = rows + i;
		    row->alloc = row_off[i + 1] - row_off[i];
		    row->cols = phgAlloc(row->alloc * sizeof (*row->cols));
		    row->ncols = 0;
		    for (; k < row_off[i + 1]; k++) {
			if (col_idx[k] == i)
			    continue;
			row->cols[row->ncols++] = col_idx[k];
			nnz++;
		    }
		}
		/* pass 2: force symmetry */
		for (i = k = 0; i < n; i++) {
		    for (; k < row_off[i + 1]; k++) {
			if (col_idx[k] == i)
			    continue;
			c = i;
			row = rows + (r = col_idx[k]);
			for (j = 0; j < row->ncols; j++)
			    if (row->cols[j] == c)
				break;
			if (j < row->ncols)
			    continue;
			if (row->ncols >= row->alloc) {
			    j = row->alloc;
			    row->alloc += 8;
			    row->cols = phgRealloc_(row->cols,
						row->alloc * sizeof(*row->cols),
						j * sizeof(*row->cols));
			}
			row->cols[row->ncols++] = c;
			nnz++;
		    }
		}
		/* pass 3: copy results back to row_off/col_idx */
		phgFree(col_idx);
		col_idx = phgAlloc(nnz * sizeof(col_idx));
		for (i = k = 0; i < n; i++) {
		    row = rows + i;
		    row_off[i] = k;
		    for (j = 0; j < row->ncols; j++)
			col_idx[k++] = row->cols[j];
		    phgFree(row->cols);
		}
		row_off[i] = k;
		phgFree(rows);
	    }		/* is_symmetric */
	    METIS_SetDefaultOptions(metis_options);
	    metis_options[METIS_OPTION_NUMBERING] = 0;	/* C numbering */
	    /* metis_options[METIS_OPTION_DBGLVL] = 31; */
	    i = METIS_NodeND(&n, row_off, col_idx, /*vwgt=*/NULL, metis_options,
				/*perm=*/ciperm, /*iperm=*/riperm);
	    if (i != METIS_OK) {
		phgWarning("Error in METIS_NodeND (ret=%d).\n", i);
		break;
	    }
	    for (i = 0; i < n; i++)
		ciperm[i] = riperm[i];
	    flag = TRUE;
#else	/* USE_METIS */
	    phgWarning("METIS is unavailable, can't do METIS ordering.\n");
#endif	/* USE_METIS */
	    break;
    }
    if (flag != TRUE) {
	for (i = 0; i < n; i++)
	    ciperm[i] = riperm[i] = i;
    }

    phgFree(row_off);
    phgFree(col_idx);

    /* Convert to CSC */
    swSparse_convert_csc(p, &a);	/* transpose */
    swSparse_destroy(p);

    /* Set ordering */
    swSparse_reorder_rows(a, SWSPARSE_IN_PLACE, NULL, riperm);
    swSparse_reorder_columns(a, SWSPARSE_IN_PLACE, NULL, ciperm);

    phgFree(ciperm);
    phgFree(riperm);

    /* factorization */
    if (_t->levels >= 0)
	swSparse_spitrf(a, &l, &u, _t->levels, 0);
    else
	swSparse_sptrf(a, &l, &u, 0);
    swSparse_destroy(a);

    swSparse_sptrsv_inspect(l, &_t->l_trsv, SWSPARSE_LOWER_TRI |
					    SWSPARSE_UNIT_DIAG);
    swSparse_sptrsv_inspect(u, &_t->u_trsv, SWSPARSE_UPPER_TRI);

    swSparse_destroy(l);
    swSparse_destroy(u);

    return 0;
}

static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    INT i, n = solver->mat->rmap->nglobal;
    double *a, *b;

    a = phgAlloc(n * sizeof(*a));
    b = phgAlloc(n * sizeof(*b));

    /* set up RHS */
    for (i = 0; i < solver->mat->rmap->nglobal; i++)
	a[i] = solver->rhs->data[i];

    /* solve */
    swSparse_d_sptrsv_execute(_t->l_trsv, a, b);
    swSparse_d_sptrsv_execute(_t->u_trsv, b, a);

    solver->residual = 0.0;
    solver->nits = 1;

    if (destroy)
	Destroy(solver);

    /* copy solution to DOFs */
    for (i = 0; i < solver->mat->rmap->nglobal; i++)
	x->data[i] = a[i];

    phgFree(a);
    phgFree(b);

    return solver->nits;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverHSparse_ = {
    "hsparse", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, FALSE, FALSE, FALSE, FALSE
};

#endif		/* USE_HSPARSE */
