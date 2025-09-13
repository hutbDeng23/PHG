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

/* $Id: solver-laspack.c,v 1.34 2022/09/21 02:35:28 zlb Exp $ */

/* The LASPack interface for PHG */

#include "phg.h"

#if USE_LASPACK

#include <string.h>

/* undefine some conflicting macros */
#undef	Dim
#undef	IsZero
#undef	IsOne

/* LASPack includes */
#include <laspack/vector.h>
#include <laspack/operats.h>

#include <laspack/qmatrix.h>

#include <laspack/itersolv.h>
#include <laspack/rtc.h>
#include <laspack/errhandl.h>
#undef	BOOLEAN				/* defined in lastypes.h */

typedef struct {
    Vector		*B;		/* RHS */
    Vector		*U;		/* solution */
    QMatrix		*A;		/* LASPack matrix */
    int			ksp_type;
    int			pc_type;
    INT			gmres_restart;
} OEM_DATA;

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)

#define Initialize	NULL
#define Finalize	NULL
#define SetPC			NULL

/* Options */

static const char *ksp_names[] = {"cg", "cgn", "cgs", "bicg",
	"bicgstab", "qmr", "ssor", "jacobi", "gmres", "chebyshev", NULL};
static const IterProcType ksp_table[] = {CGIter, CGNIter, CGSIter, BiCGIter,
	BiCGSTABIter, QMRIter, SSORIter, JacobiIter, GMRESIter, ChebyshevIter};
static int ksp_type = 0;

static const char *pc_names[] = {"none", "ilu", "jacobi", "ssor", NULL};
static const PrecondProcType pc_table[] = 
	{NULL, ILUPrecond, JacobiPrecond, SSORPrecond};
static int pc_type = 2;

static INT gmres_restart = 30;

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nLASPack options:", "\n", "laspack");

    /* Options for KSP method */
    phgOptionsRegisterKeyword("laspack_ksp_type", "KSP type",
				ksp_names, &ksp_type);
    phgOptionsRegisterInt("laspack_restart", "GMRES restart size",
			  &gmres_restart);
    /* Options for preconditioner */
    phgOptionsRegisterKeyword("laspack_pc_type", "Preconditioner type",
				pc_names, &pc_type);

    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgAlloc(sizeof(OEM_DATA));

    _t->ksp_type	= ksp_type;
    _t->pc_type		= pc_type;
    _t->gmres_restart	= gmres_restart;

    return 0;
}

static int
Create(SOLVER *solver)
{
    int N = solver->mat->rmap->nglobal;

    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }
    assert(solver->mat->rmap->nlocal == N);
    assert(solver->mat->assembled);	/* only works with assembled matrix */

    /* create matrix */
    _t->A = phgAlloc(sizeof(*_t->A));
    Q_Constr(_t->A, "LASPack matrix", N, FALSE, Rowws, Normal, TRUE);

    /* create vectors */
    _t->B = phgAlloc(sizeof(*_t->B));
    V_Constr(_t->B, "RHS", N, Normal, TRUE);
    _t->U = phgAlloc(sizeof(*_t->U));
    V_Constr(_t->U, "solution", N, Normal, TRUE);

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    if (_t->A != NULL) {
	Q_Destr(_t->A);
	phgFree(_t->A);
    }
    if (_t->B != NULL) {
	V_Destr(_t->B);
	phgFree(_t->B);
    }
    if (_t->U != NULL) {
	V_Destr(_t->U);
	phgFree(_t->U);
    }

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    Unused(solver);
    /* do nothing */
    return 0;
}

static int
AddMatrixEntries(SOLVER *solver, INT nrows, INT *rows, INT ncols, INT *cols,
			FLOAT *values)
{
    INT i, j;

    for (i = 0; i < nrows; i++, rows++) {
	/* note: this function is only called once for each row */
	Q_SetLen(_t->A, 1 + *rows, ncols);
	for (j = 0; j < ncols; j++, cols++, values++)
	    Q__SetEntry(_t->A, 1 + *rows, j, 1 + *cols, *values);
    }

    return 0;
}

static int
AddRHSEntries(SOLVER *solver, INT n, INT *ni, FLOAT *values)
{
    INT i;

    for (i = 0; i < n; i++, ni++, values++)
	V__SetCmp(_t->B, 1 + *ni, *values);

    return 0;
}

static void
laspack_callback(int it, double rNorm, double bNorm, IterIdType id)
{
    phgPrintf("    iter %d, residual %le\n", it,
			(double)(rNorm / (bNorm == 0.0 ? 1.0 : bNorm)));
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    INT i;

    /* set up initial solution */
    for (i = 0; i < solver->mat->rmap->nlocal; i++)
	V__SetCmp(_t->U, 1 + i, x->data[i]);

    assert(solver->rtol > 0. || solver->btol > 0.);
    SetRTCAccuracy(solver->rtol > solver->btol ? solver->rtol : solver->btol);
    SetGMRESRestart(_t->gmres_restart);
    if (solver->monitor)
	SetRTCAuxProc(laspack_callback);
    ksp_table[_t->ksp_type](_t->A, _t->U, _t->B, solver->maxit,
			pc_table[_t->pc_type], 1.);
    solver->residual = GetLastAccuracy();
    solver->nits = GetLastNoIter();

    if (destroy) {
	Q_Destr(_t->A);
	phgFree(_t->A);
	_t->A = NULL;
	V_Destr(_t->B);
	phgFree(_t->B);
	_t->B = NULL;
    }

    /* copy solution to DOFs */
    assert(sizeof(*_t->U->Cmp) == sizeof(FLOAT));
    memcpy(x->data, _t->U->Cmp + 1, x->map->nlocal * sizeof(*x->data));

    if (destroy) {
	V_Destr(_t->U);
	phgFree(_t->U);
	_t->U = NULL;
    }

    return solver->nits;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverLASPack_ = {
    "LASPack", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, FALSE, FALSE, FALSE
};

#endif		/* USE_LASPACK */
