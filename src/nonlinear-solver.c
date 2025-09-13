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

/* $Id: nonlinear-solver.c,v 1.7 2020/07/02 12:51:47 zlb Exp $
 *
 * Nonlinear Equations Solver. */

#ifndef NS_TEST	/*------------------------------------------------------------*/

/* force phg.h to include petscksp.h before mpi-utils.h */
#define NEED_PETSCKSP_H
#include "phg.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <setjmp.h>
#include <signal.h>

#if HAVE_FEENABLEEXCEPT || USE_PETSC
static jmp_buf jmp_env;
#if USE_OMP
# pragma omp threadprivate(jmp_env)
#endif  /* USE_OMP */
#endif	/* HAVE_FEENABLEEXCEPT || USE_PETSC */

#if HAVE_FEENABLEEXCEPT
static void
fpe_handler(int signum, siginfo_t *siginfo, void *dummy)
{
    phgInfo(1, "%s:%d: FPE trapped.\n", __FILE__, __LINE__);
    longjmp(jmp_env, 1);
}
#endif	/* HAVE_FEENABLEEXCEPT */

#if USE_PETSC

#include <petscsnes.h>

typedef struct {
    void		*ctx;	/* original user context */
    FLOAT		*x, *y, *jac;	/* work vector */
    PetscInt		*idx;
    PetscScalar		*petsc_jac;
    NS_USER_FUNC	f, fjac;
    FLOAT		tolf, tolx;
    int			n;
} CTX;

static PetscErrorCode
petsc_func(SNES snes, Vec px, Vec pf, void *pc)
{
    PetscErrorCode ierr;
    PetscScalar *ptr;
    CTX *pctx = pc;
    int i;

    ierr = VecGetArray(px, &ptr); CHKERRQ(ierr);
    for (i = 0; i < pctx->n; i++)
	pctx->x[i] = ptr[i];
    ierr = VecRestoreArray(px, &ptr); CHKERRQ(ierr);

    /* call the user function */
    if (pctx->f(pctx->x, pctx->y, pctx->ctx) == FALSE)
	longjmp(jmp_env, 1);

    ierr = VecGetArray(pf, &ptr); CHKERRQ(ierr);
    for (i = 0; i < pctx->n; i++)
	ptr[i] = pctx->y[i];
    ierr = VecRestoreArray(pf, &ptr); CHKERRQ(ierr);

    return 0;
}

static PetscErrorCode
petsc_jac(SNES snes, Vec x, Mat jac, Mat B, void *pc)
/* ref: http://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tutorials/ex1.c.html */
{
    PetscErrorCode ierr;
    PetscScalar *ptr;
    CTX *pctx = pc;
    int i;

    ierr = VecGetArray(x, &ptr); CHKERRQ(ierr);
    for (i = 0; i < pctx->n; i++)
	pctx->y[i] = ptr[i];
    ierr = VecRestoreArray(x, &ptr); CHKERRQ(ierr);

    if (pctx->fjac(pctx->y, pctx->jac, pctx->ctx) == FALSE)
	longjmp(jmp_env, 1);

    if (sizeof(FLOAT) != sizeof(PetscScalar))
	for (i = 0; i < pctx->n * pctx->n; i++)
	    pctx->petsc_jac[i] = pctx->jac[i];
    MatSetValues(B, pctx->n, pctx->idx, pctx->n, pctx->idx, pctx->petsc_jac,
			INSERT_VALUES);	
    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
    if (jac != B) {
	MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY);
    }

    return 0;
}

static PetscErrorCode
conv_test(SNES snes, PetscInt its, PetscReal xnorm, PetscReal gnorm,
	  PetscReal fnorm, SNESConvergedReason *reason, void *pc)
{
    CTX *pctx = pc;

    if (phgVerbosity > 0)
	phgInfo(-1, "Nonlinear Solver:     its = %d, residual = %le\n",
			(int)its, (double)fnorm);
    /* TODO: check pctx->tolx */
    if (fnorm <= pctx->tolf)
	*reason = SNES_CONVERGED_FNORM_ABS;
    else
	*reason = SNES_CONVERGED_ITERATING;

    return 0;
}

static int
petsc_snes(int n, NS_USER_FUNC f, NS_USER_FUNC fjac, FLOAT x[],
	   int maxits, FLOAT tolf, FLOAT tolx,
	   FLOAT *residual, void *user_ctx, FLOAT delta)
{
    /* We temporarily use PETSc with finite difference Jacobian to solve the
     * 2x2 nonlinear system of equations.
     *
     * TODO: implement Newton iterations using analytic Jacobian */
    SNES		snes;
    Vec			sol, res;	/* solution and residual */
    Mat			jac;		/* optional Jacobian */
    PetscScalar		*ptr;
    PetscErrorCode	ierr, ret;
    PetscInt		nits = 0;
    CTX			ctx;
    int i;

    /* save parameters in pctx */
    ctx.ctx = user_ctx;
    ctx.x = x;
    ctx.y = phgAlloc(n * sizeof(x[0]));
    ctx.f = f;
    ctx.fjac = fjac;
    ctx.tolf = tolf;
    ctx.tolx = tolx;
    ctx.n = n;

    if (fjac != NULL) {
	ctx.jac = phgAlloc(n * n * sizeof(FLOAT));
	ctx.idx = phgAlloc(n * sizeof(PetscInt));
	for (i = 0; i < n; i++)
	    ctx.idx[i] = i;
	if (sizeof(FLOAT) == sizeof(PetscScalar))
	    ctx.petsc_jac = (void *)ctx.jac;
	else
	    ctx.petsc_jac = phgAlloc(n * n *sizeof(PetscScalar));
    }

    /* create nonlinear solver */
    ierr = SNESCreate(PETSC_COMM_SELF, &snes); CHKERRQ(ierr);

    /* create solution and residual vectors */
    ierr = VecCreate(PETSC_COMM_SELF, &sol); CHKERRQ(ierr);
    ierr = VecSetSizes(sol, PETSC_DECIDE, n); CHKERRQ(ierr);
    ierr = VecSetFromOptions(sol); CHKERRQ(ierr);
    ierr = VecDuplicate(sol, &res); CHKERRQ(ierr);

    ierr = SNESSetFunction(snes, res, petsc_func, &ctx); CHKERRQ(ierr);
    ierr = SNESSetTolerances(snes, tolf, 0., tolx, maxits, 999999);
    CHKERRQ(ierr);
    ierr = SNESSetType(snes, "newtonls");

    /* Jacobian */
    MatCreate(PETSC_COMM_SELF, &jac);
    MatSetSizes(jac, PETSC_DECIDE, PETSC_DECIDE, n, n);
    MatSetFromOptions(jac);
    MatSetUp(jac);
    if (fjac != NULL)
	SNESSetJacobian(snes, jac, jac, petsc_jac, &ctx);
    else
	SNESSetJacobian(snes, jac, jac, SNESComputeJacobianDefault, NULL);

#if 0
    /* set linear solver */
    KSP			ksp;
    PC			pc;
    ierr = SNESGetKSP(snes, &ksp); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 0., 0., 0., 1); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPPREONLY); CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCLU); CHKERRQ(ierr);
#endif

    /* Set SNES/KSP/KSP/PC runtime options */
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

    /* Set initial solution */
    ierr = VecGetArray(sol, &ptr); CHKERRQ(ierr);
    for (i = 0; i < n; i++)
	ptr[i] = x[i];
    ierr = VecRestoreArray(sol, &ptr); CHKERRQ(ierr);

    SNESSetConvergenceTest(snes, conv_test, &ctx, NULL);

    /* Solve nonlinear system */
    if (setjmp(jmp_env) != 0)
	ret = -9999;
    else
	ret = SNESSolve(snes, PETSC_NULL, sol); 

    ierr = MatDestroy(&jac); CHKERRQ(ierr);
    phgFree(ctx.y);
    if (fjac != NULL) {
	phgFree(ctx.jac);
	if (ctx.petsc_jac != (void *)ctx.jac)
	    phgFree(ctx.petsc_jac);
	phgFree(ctx.idx);
    }

    /* Get solution, # iterations */
    ierr = SNESGetIterationNumber(snes, &nits); CHKERRQ(ierr);
    ierr = VecGetArray(sol, &ptr); CHKERRQ(ierr);
    for (i = 0; i < n; i++)
	x[i] = ptr[i];
    ierr = VecRestoreArray(sol, &ptr); CHKERRQ(ierr);

    /* Compute residual */
    if (residual != NULL) {
	FLOAT r = 0.0;
	ierr = VecGetArray(res, &ptr); CHKERRQ(ierr);
	for (i = 0; i < n; i++)
	    r += ptr[i] * ptr[i];
	ierr = VecRestoreArray(res, &ptr); CHKERRQ(ierr);
	*residual = Sqrt(r);
    }

    /* clean up */
    ierr = VecDestroy(&sol); CHKERRQ(ierr);
    ierr = VecDestroy(&res); CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);

    return ret == 0 ? nits : -1;
}

#else	/* USE_PETSC */

static int
petsc_snes(int n, NS_USER_FUNC f, NS_USER_FUNC fjac, FLOAT x[],
	   int maxits, FLOAT tolf, FLOAT tolx,
	   FLOAT *residual, void *ctx, FLOAT delta)
{
    phgError(1, "%s (%s:%d): PETSc unavailable, abort.\n",
		__func__, __FILE__, __LINE__);
    return -1;
}

#endif	/* USE_PETSC */

static struct {
    NS_USER_FUNC f;
    int n;
    FLOAT *f0, *f1, delta, delta2;
} parms;
#if USE_OMP
# pragma omp threadprivate(parms)
#endif  /* USE_OMP */

static BOOLEAN
fd_jac(FLOAT x[], FLOAT A[], void *ctx)
{
    FLOAT a;
    int i, j;

    for (i = 0; i < parms.n; i++) {
	a = x[i];
	x[i] = a - parms.delta;
	if (parms.f(x, parms.f0, ctx) != TRUE)
	    return FALSE;
	x[i] = a  + parms.delta;
	if (parms.f(x, parms.f1, ctx) != TRUE)
	    return FALSE;
	x[i] = a;
	for (j = 0; j < parms.n; j++)
	    A[j * parms.n + i] = (parms.f1[j] - parms.f0[j]) * parms.delta2;
    }

    return TRUE;
}

static int
newton(int n, NS_USER_FUNC f, NS_USER_FUNC fjac, FLOAT x[],
       int maxits, FLOAT tolf, FLOAT tolx,
       FLOAT *residual, void *ctx, FLOAT delta)
{
    int its, i;
    FLOAT *A, *y, res = 0.0, xres;


    y = phgAlloc((n + 1) * n * sizeof(*y));
    A = y + n;

    if (fjac == NULL) {
	parms.f = f;
	parms.n = n;
	parms.f0 = phgAlloc(2 * n * sizeof(FLOAT));
	parms.f1 = parms.f0 + n;
	fjac = (NS_USER_FUNC)fd_jac;	/* force type to override 'const' */
	parms.delta  = 0.5 * delta;
	parms.delta2  = 1.0 / delta;
    }
    else {
	parms.f0 = NULL;
    }

    tolf *= tolf;
    tolx *= tolx;
    for (its = 0; its <= maxits; its++) {
	/* evaluate function values */
	if (f(x, y, ctx) != TRUE) {
	    /* user break */
	    its = -1;
	    break;
	}

	/* compute residual */
	res = 0.0;
	for (i = 0; i < n; i++)
	    res += y[i] * y[i];

	if (phgVerbosity > 0)
	    phgInfo(-1, "Nonlinear Solver:     its = %d, residual = %le\n",
			(int)its, (double)Sqrt(res));

	/* check converegence */
	if (res <= tolf)
	    break;		/* tolf met */

	/* check number of iteration */
	if (its == maxits) {
	    /* maxits attained */
	    its = -1;
	    break;
	}

	/* evaluate Jacobian */
	if (fjac(x, A, ctx) != TRUE) {
	    /* user break */
	    its = -1;
	    break;
	}

	/* solve linear system of equations */
	if (phgSolverDenseSolver(n, 1, A, y) != TRUE) {
	    its = -1;
	    break;
	}

	/* update approximate solution */
	xres = 0.0;
	for (i = 0; i < n; i++) {
	    x[i] -= y[i];
	    xres += y[i] * y[i];
	}

	if (xres <= tolx)
	    break;		/* tolx met */
    }

    phgFree(y);
    phgFree(parms.f0);

    if (residual != NULL)
	*residual = Sqrt(res);

    return its;
}

int
phgNonlinearSolver(int n, NS_USER_FUNC f, NS_USER_FUNC fjac, FLOAT x[],
		   int maxits, FLOAT tolf, FLOAT tolx,
		   FLOAT *residual, void *ctx)
/* Solves a system of nonlinear equations and returns the number of iterations
 * (-1 => failure, including FPE and exceeding 'maxits')
 *
 * Note:
 *   1. 'n' is the number of equations and unknowns (dimension of f and x).
 *   2. 'x' contains the initial (input) and final (output) approx. solution.
 *   3. if 'fjac' == NULL then finite difference Jacobian is used.
 *   4. if 'residual' != NULL then it returns the l2 norm of the residual.
 */
{
    static BOOLEAN initialized = FALSE;
    static int ns_id = 0;
    static const char *ns_names[] = {"newton", "petsc", NULL};
    static FLOAT delta = -1.0;

    int nits = 0;

    if (initialized != TRUE) {
	/* register options */
	initialized = TRUE;
	phgOptionsRegisterTitle("\nNonlinear solver options:", "\n",
			"nonlinear_solver");
	phgOptionsRegisterKeyword("-ns_solver", "Nonlinear solver",
			ns_names, &ns_id);
	delta = Sqrt(FLOAT_EPSILON);
	phgOptionsRegisterFloat("-ns_delta",
			"Delta for computing finite difference Jacobian",
			&delta);
	return 0;
    }

#if HAVE_FEENABLEEXCEPT
    extern void (*__phg_user_fpe_handler)(int, siginfo_t *, void *);
    void (*fpe_handler_bak)(int, siginfo_t *, void *) = __phg_user_fpe_handler;
    __phg_user_fpe_handler = fpe_handler;
#endif	/* HAVE_FEENABLEEXCEPT */

    if (phgVerbosity > 0) {
	phgInfo(-1, "Nonlinear Solver: --- entering nonlinear solver.\n");
	phgInfo(-1, "Nonlinear Solver: *** %s Jacobian.\n",
			fjac == NULL ? "finite difference" : "analytic");
    }

    switch (ns_id) {
	case 0:	/* Newton */
	    if (phgVerbosity > 0)
		phgInfo(-1, "Nonlinear Solver: *** built-in Newton method.\n");
	    nits = newton(n, f, fjac, x, maxits, tolf, tolx,
			  residual, ctx, delta);
	    break;
	case 1:	/* PETSc SNES */
	    if (phgVerbosity > 0)
		phgInfo(-1, "Nonlinear Solver: *** PETSc SNES.\n");
	    nits = petsc_snes(n, f, fjac, x, maxits, tolf, tolx,
			      residual, ctx, delta);
	    break;
    }

#if HAVE_FEENABLEEXCEPT
    __phg_user_fpe_handler = fpe_handler_bak;
#endif	/* HAVE_FEENABLEEXCEPT */

    if (phgVerbosity > 0) {
	if (nits < 0)
	    phgInfo(-1, "Nonlinear Solver: *** unconverged or interrupted.\n");
	else
	    phgInfo(-1, "Nonlinear Solver: *** converged after %d it%s.\n",
			nits, nits > 1 ? "s" : "");
    }

    return nits;
}

#else	/*------------------------- defined(NS_TEST) -------------------------*/

#include "phg.h"

static BOOLEAN
func(const FLOAT x[], FLOAT y[], void *ctx)
{
    y[0] = x[0] * x[0] + 2.0 * x[0] * x[1] + x[1] * x[1] - 4.0;
    y[1] = x[0] * x[0] - 2.0 * x[0] * x[1] + x[1] * x[1] - 9.0;
    return TRUE;
}

static BOOLEAN
fjac(const FLOAT x[], FLOAT y[], void *ctx)
{
    y[0] = 2. * (x[0] + x[1]); y[1] = 2. * (x[0] + x[1]);
    y[2] = 2. * (x[0] - x[1]); y[3] = 2. * (-x[0] + x[1]);
    return TRUE;
}

int
main(int argc, char *argv[])
{
    BOOLEAN fd_jac = FALSE;
    FLOAT x[2] = {2.0, -0.25};	/* exact sol: x=(2+3)/2=2.5, y=(2-3)/2=-0.5 */

    phgOptionsRegisterNoArg("-fd_jac", "Use finite difference Jacobian",
				&fd_jac);
    phgInit(&argc, &argv);

    phgNonlinearSolver(2, func, fd_jac == FALSE ? fjac : NULL, x,
						20, 1e-10, 1e-10, NULL, NULL);
    printf("solution: %lg %lg\n", x[0], x[1]);

    phgFinalize();
    return 0;
}

#endif	/*------------------------- defined(NS_TEST) -------------------------*/
