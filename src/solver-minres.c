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

/* $Id: solver-minres.c,v 1.26 2022/09/21 06:37:44 zlb Exp $ */

/* Interface to the MINRES solver:
 *	http://www.stanford.edu/group/SOL/software/minres.html
 *
 * Note: must use package ftp://ftp.cc.ac.cn/pub/RPMS/minres*.rpm */

#include "phg.h"

#if USE_MINRES	/* whole file */

#include <string.h>	/* memcpy */

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL

typedef struct {
    int		pc_type;
    char        *pc_opts;
    BOOLEAN	verbose;

    /* state and scratch variables */
    VEC		*vy, *vr;
    double	rb_norm, r0_norm, tol;
} OEM_DATA;

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)

enum {PC_TYPE_NONE = 0, PC_TYPE_DIAG = 1, PC_TYPE_SOLVER = 2};
static const char *pc_names[] = {"none", "diag", "solver", NULL};
int pc_type = 1;

static char *pc_opts = NULL;
static BOOLEAN verbose = FALSE;		/* MINRES verbose */

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nThe MINRES solver options:", "\n", "minres");
    phgOptionsRegisterKeyword("-minres_pc_type", "Preconditioner type",
				pc_names, &pc_type);
    phgOptionsRegisterString("-minres_pc_opts", "Options passed to the "
				"preconditioner", &pc_opts);
    phgOptionsRegisterNoArg("-minres_verbose", "MINRES verbose", &verbose);

    return 0;
}

static int
Init(SOLVER *solver)
{
    phgInfo(2, "Create %s solver\n", solver->oem_solver->name);

    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));
    _t->pc_type = pc_type;
    _t->pc_opts = (pc_opts == NULL ? NULL : strdup(pc_opts));
    solver->pc_option_flag = phgOptionsIfUsed("-minres_pc_type");
    _t->verbose = verbose;

    return 0;
}

static int
Create(SOLVER *solver)
{
    /* Set up PC if not already set */

    if (solver->pc == NULL && _t->pc_type == PC_TYPE_SOLVER) {
	/* Create PC solver */
	SOLVER *pc;
	phgOptionsPush();
	phgSolverSetDefaultSuboptions();
	phgOptionsSetOptions("-solver_rtol 0 -solver_atol 0 -solver_btol 0"
			     " -solver_maxit 1 +solver_warn_maxit");
	phgOptionsSetOptions(_t->pc_opts);
	pc = phgSolverMat2Solver(SOLVER_DEFAULT, solver->mat);
	phgOptionsPop();
	if (pc->oem_solver == SOLVER_MINRES) {
	    /* avoid regression */
	    ((OEM_DATA *)pc->oem_data)->pc_type = PC_TYPE_NONE;
	}
	phgSolverSetPC_(solver, pc, NULL);
	solver->pc->solver_owned = pc;	/* destroy pc afterwards */
    }

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    phgMatDestroy(&solver->mat);
    solver->sent_to_oem = TRUE;

    if (_t->pc_opts != NULL)
	phgFree(_t->pc_opts);
    _t->pc_opts = NULL;

    phgFree(solver->oem_data);

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    return 0;
}

/*--------------------------- callback functions ------------------------*/
static SOLVER *solver = NULL;

static void
Aprod(int *n, double x[], double y[])
/* computes y := Ax */
{
    assert(*n == solver->mat->rmap->nlocal);
    assert(FT_PHG == FT_DOUBLE);
    _t->vy->data = (FLOAT *)x;
    _t->vr->data = (FLOAT *)y;
    phgMatVec(MAT_OP_N, 1.0, solver->mat, _t->vy, 0.0, &_t->vr);
    _t->vy->data = _t->vr->data = NULL;
}

static void
Msolve(int *n, double x[], double y[])
/* solves My = x */
{
    int i, pc_type = _t->pc_type;

    assert(*n == solver->mat->rmap->nlocal);

    switch (pc_type) {
	case PC_TYPE_NONE:	/* no preconditioner */
	    if (*n > 0)
		memcpy(y, x, *n * sizeof(y[0]));
	    break;
	case PC_TYPE_DIAG:	/* Diagonal preconditioner */
	    if (solver->mat->diag == NULL)
		phgMatSetupDiagonal(solver->mat);
	    for (i = 0; i < solver->mat->rmap->nlocal; i++)
		y[i] = x[i] / solver->mat->diag[i];
	    break;
	case PC_TYPE_SOLVER:	/* External solver: calling solver->pc */
	    assert(solver->pc != NULL && solver->pc->pc_proc != NULL);
	    assert(FT_PHG == FT_DOUBLE);
	    _t->vy->data = (FLOAT *)x;
	    _t->vr->data = (FLOAT *)y;
	    solver->pc->pc_proc(solver->pc->ctx, _t->vy, &_t->vr);
	    _t->vy->data = _t->vr->data = NULL;
	    break;
    }
}

static double
dot_product(double x[], double y[])
{
    int i;
    double d = 0.0;

    for (i = 0; i < solver->mat->rmap->nlocal; i++)
	d += x[i] * y[i];

#if USE_MPI
    if (solver->mat->rmap->nprocs > 1) {
	double d0;
	MPI_Allreduce(&d, &d0, 1, MPI_DOUBLE, MPI_SUM, solver->mat->rmap->comm);
	d = d0;
    }
#endif	/* USE_MPI */

    return d;
}

static int
check_res(int *n, double *x, double *rnorm, int *itn)
/* Note: returns !0 will stop MINRES iterations */
{
    assert(*n == solver->mat->rmap->nlocal);

    if (*itn == 1) {
	double d;

	if (solver->monitor) {
	    phgPrintf("======================= MINRES ====================\n");
	    phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
			(double)solver->atol, (double)solver->rtol,
			(double)solver->btol);
	    phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
	    phgPrintf("-----   ------------   -------------  -------------\n");
	}
	if (*rnorm == 0.0) {
	    /* RHS = 0: set x = 0 and exit */
	    if (*n > 0)
		memset(x, 0, *n * sizeof(*x));
	    if (solver->monitor)
		phgPrintf("% 5d   %12le   %12le   %12le\n", *itn, 0., 0., 0.);
	    return 1;	/* convergence achieved */
	}
	_t->r0_norm = *rnorm;
	_t->tol = solver->atol;
	d = _t->r0_norm * solver->rtol;
	if (_t->tol <= d)
	    _t->tol = d;
	d = _t->rb_norm * solver->btol;
	if (_t->tol <= d)
	    _t->tol = d;
	/* inverse r0_norm and rb_norm to avoid a few divisions */
	if (solver->monitor) {
	    _t->r0_norm = 1.0 / _t->r0_norm;
	    _t->rb_norm = 1.0 / _t->rb_norm;
	}
    }

    if (solver->monitor)
	phgPrintf("% 5d   %12le   %12le   %12le\n", *itn,
			(double)*rnorm,
			(double)(*rnorm * _t->r0_norm),
			(double)(*rnorm * _t->rb_norm));

    return *rnorm <= _t->tol ? 1 : 0;
}

static int
SetPC(SOLVER *solver)
{
    _t->pc_type = PC_TYPE_SOLVER;

    return 0;
}

static int
Solve(SOLVER *solver0, VEC *x, BOOLEAN destroy)
{
    extern void FC_FUNC(minreswrapper,MINRESWRAPPER)(
		int *n,
		void (*Aprod)(int *, double *, double *),
		void (*Msolve)(int *, double *, double *),
		double *b, double *shift, int *checkA, int *precon,
		double *x, int *itnlim, int *nout, double *rtol,
		int *istop, int *itn, double *Anorm, double *Acond,
		double *rnorm, double *Arnorm, double *ynorm,
		double (*dot_product)(double *, double *),
		int (*check_res)(int *, double *, double *, int *) );
    /* input variables */
    double shift = 0.0, rtol = 0.0;
    int n = solver0->mat->rmap->nlocal, checkA = 0, precon = 1,
	itnlim = solver0->maxit, nout = 0;
    /* output variables */
    double Anorm, Acond, rnorm, Arnorm, ynorm;
    int istop, itn;
    FLOAT *y, *r;
    SOLVER *solver_save;

    assert(FT_PHG == FT_DOUBLE);
    solver_save = solver;
    solver = solver0;

    if (solver->mat->rmap->nglobal <= 0 ||
	(_t->rb_norm = phgVecNorm2(solver->rhs, 0, NULL)) == 0.0) {
	/* empty system */
	solver->nits = 0;
	solver->residual = 0;
	solver = solver_save;
	return 0;
    }

    /* Note: this version of MINRES always sets the initial solution to 0.
     * as a work around we solve the residual equation, then update x using
     * the solution */

    /* compute initial residual */
    _t->vy = phgMapCreateVec(x->map, 1);
    if (phgVecNormInfty(x, 0, NULL) == 0.0) {
	/* the initial solution is zero */
	_t->vr = phgVecCopy(solver->rhs, NULL);
    }
    else {
	_t->vr = phgMatVec(MAT_OP_N, 1.0, solver->mat, x, 0.0, NULL);
	phgVecAXPBY(1.0, solver->rhs, -1.0, &_t->vr);
	if (phgVecNormInfty(_t->vr, 0, NULL) == 0.0) {
	    /* initial residual is zero: direct return */
	    phgVecDestroy(&_t->vy);
	    phgVecDestroy(&_t->vr);
	    solver->nits = 0;
	    solver->residual = 0;
	    solver = solver_save;
	    return 0;
	}
    }

    /* save the data pointers of the vectors for use in MINRES */
    r = _t->vr->data;	/* the residual */
    y = _t->vy->data;
    _t->vy->data = _t->vr->data = NULL;

    if (solver->mat->rmap->rank == 0 && _t->verbose)
	nout = 6;

    assert(FT_PHG == FT_DOUBLE);
    FC_FUNC(minreswrapper,MINRESWRAPPER)(
		   &n, Aprod, Msolve, (double *)r, &shift, &checkA, &precon,
		   (double *)y, &itnlim, &nout, &rtol, &istop, &itn, &Anorm,
		   &Acond, &rnorm, &Arnorm, &ynorm, dot_product, check_res);
    solver->residual = rnorm;
    solver->nits = itn;

    /* check the return value (see minresModule.f90 in the MINRES source) */
    switch (istop) {
	case -1:	/* got an eigen vector: Ax = lambda Mx */
	    break;
	case 0:		/* b = 0, should not happen here */
	    phgError(1, "%s:%d: unexpected.\n", __FILE__, __LINE__);
	    break;
	case 1:		/* normal convergence */
	    break;
	case 2:		/* convergence to machine precision */
	    break;
	case 3:		/* convergence to an eigenvector of A: A*x=shift*x */
	    break;
	case 4:		/* A very ill conditioned, convergence unassured */
	    break;
	case 5:		/* maxit achieved */
	    break;
	case 6:		/* A not symmetric: no iteration performed */
	    phgPrintf("(%s:%d) Warning: unsymmetric matrix, solution "
		      "unconverged.\n", __FILE__, __LINE__);
	    solver->residual = 1e20;
	    break;
	case 7:		/* M unsymmetric: no iteration performed */
	    phgPrintf("(%s:%d) Warning: unsymmetric preconditioner, solution "
		      "unconverged.\n" __FILE__, __LINE__);
	    solver->residual = 1e20;
	    break;
	case 8:		/* M not positive definite */
	    phgPrintf("(%s:%d) Warning: non positive preconditioner, solution "
		      "unconverged.\n", __FILE__, __LINE__);
	    solver->residual = 1e20;
	    break;
	default:	/* b = 0, should not happen here */
	    phgError(1, "%s:%d: unexpected.\n", __FILE__, __LINE__);
    }

    /* restore data pointers of the vectors */
    _t->vy->data = y;
    _t->vr->data = r;

    /* update x */
    phgVecAXPBY(1.0, _t->vy, 1.0, &x);

    phgVecDestroy(&_t->vy);
    phgVecDestroy(&_t->vr);

    if (destroy) {
	phgMatDestroy(&solver->mat);
	solver->sent_to_oem = TRUE;
    }

    solver = solver_save;
    return solver0->nits;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverMINRES_ = {
    "MINRES", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_SYM, TRUE, TRUE, TRUE, TRUE
};

#endif	/* USE_MINRES */
