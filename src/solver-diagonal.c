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

/* $Id: solver-diagonal.c,v 1.7 2022/09/21 02:35:28 zlb Exp $ */

/* The diagonal solver. */

#include "phg.h"

#define Initialize		NULL
#define Finalize		NULL
#define Init			NULL
#define RegisterOptions		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

static int
Create(SOLVER *solver)
{
    return 0;
}

static int
Destroy(SOLVER *solver)
{
    phgMatDestroy(&solver->mat);
    solver->sent_to_oem = TRUE;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    return 0;
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    MAT *mat = solver->mat;
    MAP *map = mat->rmap;
    FLOAT *rhs, *p;
    INT i, n0;

    if (map->nlocal <= 0)
	return 0;

    p = x->data;
    n0 = map->partition[map->rank];
    rhs = solver->rhs->data;
#if 0
    if (mat->type == PHG_PACKED)
	phgMatUnpack(mat);
    assert(mat->type == PHG_UNPACKED);
    for (i = 0; i < map->nlocal; i++, rhs++, p++) {
	MAT_ROW *row = mat->rows + i;
	if (row->ncols == 0)
	    phgError(1, "empty row, global index = %d\n", i + n0);
	if (row->ncols != 1 || *row->cols != i)
	    phgError(1, "the diagonal solver only solves diagonal systems.\n");
	if (*row->data == 0.0)
	    phgError(1, "%s:%d, zero diagonal entry found in the matrix.\n",
			__FILE__, __LINE__);
	*p = *rhs / *row->data;
    }
#else
    if (mat->diag == NULL)
	phgMatSetupDiagonal(mat);
    for (i = 0; i < map->nlocal; i++) {
	if (mat->diag[i] == 0.0)
	    phgError(1, "%s:%d, zero diagonal entry, global index = %d\n",
			__FILE__, __LINE__, i + n0);
	p[i] = rhs[i] / mat->diag[i];
    }
#endif

    if (destroy) {
	phgMatDestroy(&solver->mat);
	solver->sent_to_oem = TRUE;
    }

    solver->residual = 0.0;	/* fake residual norm */
    return solver->nits = 1;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverDiagonal_ = {
    "Diagonal", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, FALSE, TRUE, TRUE, FALSE
};
