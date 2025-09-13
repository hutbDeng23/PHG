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

/* $Id: solver-x9amg.c,v 1.15 2022/09/21 02:35:28 zlb Exp $ */

/* The XTU-9Suo AMG interface */

#include "phg.h"

#if USE_X9AMG	/* whole file */

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

#include <x9amg.h>

typedef XJ_AMG_HO_CONTEXT OEM_DATA;

#define	_t	((OEM_DATA *)solver->oem_data)

static const char *solver_list[] = {"amg", "pcg-amg", "gmres-amg","pj-amg","pcg-pj","gmres-pj", NULL};
static int solver_id = 2;
static const char *relax_types[] = {"bgs", "jacobi", NULL};
static int relax_type = 0;
static INT relax_sweeps = 1;
static FLOAT linear_tol = 1.0e-6;
static INT linear_max_iter = 1;
static FLOAT THETA = 0.2;

/* X9 AMG interface functions:
 *
 * void *XJ_AMG_HO_CONTEXTCreate();
 * int XJ_AMG_HO_CONTEXTDestroy(XJ_AMG_HO_CONTEXT *XJ_AMG_HO_vdata);
 * int XJ_AMG_HO_CONTEXTtol_maxit_print(XJ_AMG_HO_CONTEXT *XJ_AMG_HO_data,
 * 		int max_iter,double rtol,int print_lever );
 * int xj_ho_amg_setup(MPI_Comm comm, XJ_AMG_HO_CONTEXT *XJ_AMG_HO_vdata,
 * 		int nfree_in, int num_freedom_linear_in,
 *		int *ia,int *ja, double *a,const int MYINFO);
 * int xj_ho_amg_bgs_solve(MPI_Comm  comm, XJ_AMG_HO_CONTEXT *XJ_AMG_HO_vdata,
 * 		double *rhs, double **xx_data_out, double **residual_out,
 * 		int *bgs_iters_out);
 * xj_ho_amg_bgs_solve(MPI_Comm  comm, XJ_AMG_HO_CONTEXT *XJ_AMG_HO_vdata,
 *		double *rhs, double *xx_data, double *residual_data,
 *		int *bgs_iters_out);
 * int xj_ho_amg_solve(MPI_Comm comm, XJ_AMG_HO_CONTEXT *XJ_AMG_HO_vdata,
 * 		double *rhs, double *xx_data, double *residual_data,
 * 		int *bgs_iters_out);
*/

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nx9AMG options:", "\n", "x9amg");
    phgOptionsRegisterKeyword("-x9amg_solver", "x9AMG solver",
				solver_list, &solver_id);
    phgOptionsRegisterKeyword("-x9amg_relax_type", "Relaxation type",
				relax_types, &relax_type);
    phgOptionsRegisterInt("-x9amg_relax_sweeps", "Relaxation sweeps",
				&relax_sweeps);
    phgOptionsRegisterInt("-x9amg_coarse_maxit",
			  "Max # of iterations for the coarse solver",
			  &linear_max_iter);
    phgOptionsRegisterFloat("-x9amg_coarse_tol",
			  "Tolerance for the for the coarse solver",
			  &linear_tol);
    phgOptionsRegisterFloat("-x9amg_strength", "?????????????", &THETA);

    return 0;
}

static int
Init(SOLVER *solver)
{
phgPrintf("=============== solver_id = %d\n", solver_id);
    solver->oem_data = XJ_AMG_HO_CONTEXTCreate(solver_id + 1);

    _t->bgs_smooth = relax_type;
    _t->linear_tol = linear_tol;
    _t->linear_max_iter = linear_max_iter;
    _t->THETA = THETA;

    return 0;
}

static int
Create(SOLVER *solver)
{
    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    XJ_AMG_HO_CONTEXTDestroy(_t);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    MAT *mat = solver->mat;
    MAP *map = mat->rmap;
    int i, j, *index, *cols, *pi;
    int Nlocal = map->nlocal, n0 = map->partition[map->rank];
    MAT_ROW *row;
    double *a, *pd;

    assert(map->dofs != NULL && map->dofs[0] != NULL);
    assert(solver->mat->type != PHG_DESTROYED);
    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }

    if (solver->sent_to_oem == TRUE)
	return 0;
    solver->sent_to_oem = TRUE;

    if (mat->type == PHG_PACKED)
	phgMatUnpack(mat);

    index = phgAlloc((Nlocal + 1) * sizeof(*index));
    cols = phgAlloc((mat->nnz_o + mat->nnz_d) * sizeof(*cols));
    a = phgAlloc((mat->nnz_o + mat->nnz_d) * sizeof(*a));

    pi = cols;
    pd = a;
    index[0] = 0;
    for (i = 0; i < Nlocal; i++) {
        if (mat->blocks == NULL) {
	    row = mat->rows + i;
	    for (j = 0; j < row->ncols; j++) {
		INT r;
		*(pd++) = row->data[j];
		r = row->cols[j];
#if USE_MPI
		*(pi++) = (r < Nlocal ? r + n0 : mat->O2Gmap[r - Nlocal]);
#else	/* USE_MPI */
		*(pi++) = r;
#endif	/* USE_MPI */
	    }
	    if (mat->refcount == 0) {
		phgFree(row->cols);
		phgFree(row->data);
		row->cols = NULL;
		row->data = NULL;
		row->ncols = row->alloc = 0;
	    }
	}
        else {
	    phgError(1, "unimplemented.\n");
            /*phgMatGetBlockMatrixRow(mat, i, row);*/
	}
	index[i + 1] = (int)(pi - cols);
    }

    if (mat->refcount == 0)
	phgMatFreeMatrix(mat);

    assert(sizeof(FLOAT) == sizeof(double));

    _t->comm = mat->rmap->comm;

    xj_ho_amg_setup(_t->comm, _t, Nlocal, map->dofs[0]->g->nvert_owned,
		    index, cols, a, 0);

    phgFree(index);
    phgFree(cols);
    phgFree(a);

    return 0;
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    Assemble(solver);

    assert(sizeof(INT) == sizeof(int));
   // _t->bgs_tol = solver->rtol;
    _t->rtol = solver->rtol;
    _t->max_iter = solver->maxit;
    _t->print_or_no = solver->monitor;
   
     XJ_AMG_HO_CONTEXTtol_maxit_print(_t, solver-> maxit,
						solver->rtol, solver->monitor);
//    xj_ho_amg_bgs_solve(_t, &x->data, &solver->nits);
double *dummy = phgAlloc(x->map->nlocal * sizeof(*dummy));
    xj_ho_amg_solve(solver->mat->rmap->comm, _t, solver->rhs->data, x->data, dummy, &solver->nits);
phgFree(dummy);

    if (destroy) {
	XJ_AMG_HO_CONTEXTDestroy(_t);
	solver->oem_data = NULL;
    }

    solver->residual = 0.0;	/* fake residual norm */

    return solver->nits;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverX9AMG_ = {
    "x9AMG", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy,
    AddMatrixEntries, AddRHSEntries,
    Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, FALSE, TRUE, FALSE, FALSE
};

#endif	/* USE_X9AMG */
