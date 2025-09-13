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

/* $Id: eigen-arpack.c,v 1.81 2022/09/16 22:41:49 zlb Exp $ */

#include "phg.h"

#if USE_ARPACK

#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>

static VEC *va, *vb;

/* cmdline options */
enum {SHIFT_INVERT, REGULAR_INVERT};
static const int mode_list[] = {SHIFT_INVERT, REGULAR_INVERT};
static const char *mode_name[] = {"shift-invert", "regular-invert", NULL};
static int mode_index = 0;	/* default = shift-invert */
static BOOLEAN keep_unconverged = FALSE; /* whether keep unconverged evals */
static int mode;		/* will be set to mode_list[mode_index] */

static INT solver_maxit = 2000;
static FLOAT solver_rtol = 1e-15;
static int solver_id = 0;
static char *solver_opts = NULL;
static INT cmdline_ncv = 0;

static int matvec_count, solve_count;

typedef unsigned int LOGICAL;	/* FORTRAN LOGICAL */
#define F77_True	1
#define F77_False	0

/* TODO: use the following code to figure out the Fortran logical type and
   the values for .true. and .false.:
   ------------ tmp.f:
      logical a(2)
      data a/.true., .false./
      call tmpc(a(1))
      call tmpc(a(2))
      stop
      end
  ------------- tmpc.c:
  #include <stdio.h>
  void tmpc_(unsigned char *ptr)
  {
    static unsigned char *ptr0 = NULL;
    int i, n;
    if (ptr0 == NULL) {
	ptr0 = ptr;
	return;
    }
    n = (int)(ptr - ptr0);
    printf("LOGICAL size = %d\n .TRUE.  =", n);
    for (i = 0; i < n; i++)
	printf(" 0x%02x", ptr0[i]);
    printf("\n .FALSE. =");
    for (i = 0; i < n; i++)
	printf(" 0x%02x", ptr[i]);
    printf("\n");
  }
 */

extern void F77_FUNC(pdsaupd,PDSAUPD) (
    int   *comm,	/* the communicator */
    int   *ido,		/* what to do next */
    const char *bmat,	/* type of the matrix B ('I' or 'G') */
    int   *n,		/* dimension of the eigenproblem */
    const char *which,	/* 'LA' (algebraicly largest),
			 * 'SA' (algebraicly smallest),
			 * 'LM' (largest in magnitude),
			 * 'SM' (smallest in magnitude),
			 * 'BE' (half from each end of the spectrum) */
    int    *nev,	/* number of eigenvalues (Ritz values) to compute */
    double *tol,	/* relative accuracy (<=0 => use default) */
    double resid[],	/* residual vector, array of length N */
    int    *ncv,	/* number of columns of the matrix V (<= N) */
    double v[],		/* the Lanczos basis vectors, N by NCV array (out) */
    int    *ldv,	/* leading dimension of V */
    int    iparam[11],	/* parameters */
    int    ipntr[11],	/* starting locations in the WORKD and WORKL */
    double workd[],	/* work array of length 3*N */
    double workl[],	/* work array of length LWORKL */
    int    *lworkl,	/* sizeof WORKL, at least NCV**2 + 8*NCV */
    int    *info,	/* Input:
			 *	0:  a randomly initial residual vector is used
			 *	!0: RESID contains the initial residual vector
			 * Output:
			 *	0: normal exit
			 *	1: maximum number of iterations taken.
			 *	... ... */
    /* lengths of character strings */
    int bmat_len, int which_len);

extern void F77_FUNC(pdseupd,PDSEUPD) (
    int    *comm,
    LOGICAL *rvec,	/* FALSE: Ritz values only, TRUE: Ritz vectors */
    const char *howmny,	/* 'A': all, 'S': some (specified by SELECT[]) */
    LOGICAL select[],	/* logical array of dimension NCV */
    double d[],		/* eigenvalues, array of dimension NEV (output) */
    double z[],		/* Ritz vectors, N by NEV array (may use V) */
    int    *ldz,	/* leading dimension of Z */
    double *sigma,	/* represents the shift */
    /* the following arguments are the same as for pdsaupd */
    const char *bmat,
    int    *n,
    const char *which,
    int    *nev,
    double *tol,
    double resid[],
    int    *ncv,
    double v[],
    int    *ldv,
    int    iparam[11],
    int    ipntr[11],
    double workd[],
    double workl[],
    int    *lworkl,
    int    *info,
    /* sizes of strings */
    int howmny_len, int bmat_len, int which_len);

static int
mat_vec(MAT *A, SOLVER *solver, FLOAT *x, FLOAT *y)
/* computes y := inv[B] * A * x */
{
    va->data = x;
    vb->data = y;
    va->nvec = vb->nvec = 1;

    if (A != NULL) {
	matvec_count++;
	phgMatVec(MAT_OP_N, 1.0, A, va, 0.0, &vb);
    }
    else {	/* A == NULL ==> A == I */
	memcpy(y, x, va->map->nlocal * sizeof(*x));
    }

    /* If solver != NULL solve with RHS = A*x and set x <-- A*x */
    if (solver != NULL) {
	FLOAT res;
	int nits;
	INT verb = phgVerbosity;

	if (mode == REGULAR_INVERT) {
	    /* copy vb (A*x) to va as required by the regular invert mode */
	    memcpy(va->data, vb->data, vb->map->nlocal * sizeof(*vb->data));
	    solver->rhs->data = va->data;
	}
	else {
	    memcpy(solver->rhs->data, vb->data,
			vb->map->nlocal * sizeof(*vb->data));
	}

	if (phgVerbosity != 0)
	    phgVerbosity--;
#if 0
	/* call PHG's built-in PCG solver (slow!) */
	nits = solver_pcg(B, va, &vb, NULL, solver_maxit, 1e-17, 1e-30, 1e-30, &res, FALSE, 0);
#else
	/* call an external solver */
	solver->rhs->assembled = TRUE;
	solve_count++;
	nits = phgSolverVecSolve(solver, FALSE, vb);
	if (mode == REGULAR_INVERT)
	    solver->rhs->data = NULL;
	res = solver->residual;
#endif
	phgInfo(2, "linear solver: nits = %d, res = %lg\n", nits, (double)res);
	phgVerbosity = verb;
    }

    va->data = NULL;
    vb->data = NULL;

    return 0;
}

static int
comp_double(const void *p0, const void *p1)
{
    if (*(const double *)p0 == *(const double *)p1)
	return 0;
    else
	return (*(const double *)p0 < *(const double *)p1) ? -1 : 1;
}

EigenSolverPrototype(phgEigenSolverARPACK)
/* (A, B, n, which, tau, evals, evecs, eps, itmax, nit) */
{
    static BOOLEAN initialized = FALSE;
    VEC *vec, *vec0 = NULL;

    int comm;
    const char *bmat, *whch, *howmny;
    int i, it, ido, nloc, nev, ncv, info;
    int ldv, iparam[11], ipntr[11], lworkl;
    double sigma, tol, res_max, *resid, *v, *workd, *workl, *bounds;
    LOGICAL rvec, *select;
    BOOLEAN flag;
    SOLVER *solver;

#define IPARAM(i)	iparam[(i) - 1]
#define IPNTR(i)	ipntr[(i) - 1]
#define WORKD(i)	workd[(i) - 1]
#define WORKL(i)	workl[(i) - 1]

    if (!initialized) {
	/* register options */
	initialized = TRUE;
	phgOptionsRegisterTitle("\nARPACK options:", "\n", "arpack");
	phgOptionsRegisterKeyword("-arpack_mode", "ARPACK mode", mode_name,
					&mode_index);
	phgOptionsRegisterInt("-arpack_ncv", "ARPACK ncv (0 means use default)",
					&cmdline_ncv);

	/* set default solver to the first available parallel direct solver */
	for (solver_id = 0;
	     phgSolverNames[solver_id] != NULL
	     && (phgSolverList[solver_id]->iterative == TRUE
#if USE_MPI
		 || phgSolverList[solver_id]->parallel == FALSE
#endif	/* USE_MPI */
		 || phgSolverList[solver_id] == SOLVER_DIAGONAL);
	     solver_id++);
	if (phgSolverNames[solver_id] == NULL)
	    solver_id = 0;		/* fall back to the first solver */
	phgOptionsRegisterKeyword("arpack_solver", "ARPACK linear solver",
					phgSolverNames, &solver_id);
	phgOptionsRegisterInt("arpack_solver_maxit", "ARPACK linear solver "
					"maxit", &solver_maxit);
	phgOptionsRegisterFloat("arpack_solver_rtol", "ARPACK linear solver "
					"rtol", &solver_rtol);
	phgOptionsRegisterString("-arpack_solver_opts",
					"Options for the linear solver",
					&solver_opts);
	phgOptionsRegisterNoArg("-arpack_keep_unconverged",
					"Keep unconverged eigen pairs",
					&keep_unconverged);
	return 0;
    }

    assert(FT_PHG == FT_DOUBLE);

    if (n <= 0 || A->cmap->nglobal <= 0)
	return 0;

    vec = *evecs;
    MagicCheck(VEC, vec);

    mode = mode_list[mode_index];

    va = phgMapCreateVec(A->cmap, 1);
    vb = phgMapCreateVec(A->cmap, 1);
    phgFree(va->data);
    va->data = NULL;
    phgFree(vb->data);
    vb->data = NULL;

    if (vec == NULL) {
	vec0 = NULL;
	*evecs = vec = phgMapCreateVec(A->cmap, n);
    }
    else {
	/* save initial input vectors */
	vec0 = phgMapCreateVec(vec->map, 0);
	vec0->data = vec->data;
	vec0->nvec = vec->nvec;
	vec->nvec = n;
	vec->data = phgAlloc(n * vec->map->nlocal * sizeof(*vec->data));
    }

    assert(A->cmap->nglobal > 1);

    comm = MPI_Comm_c2f(A->cmap->comm);

    bmat = (B == NULL ? "I" : "G");
    ido = 0;
    nloc = A->cmap->nlocal;
    nev = n;
    ncv = 2 * nev;
    info = 0;
    ldv = nloc;
    tol = eps;

    howmny = "All";
    rvec = F77_True;
    sigma = tau;
    flag = TRUE;
    solver = NULL;

    if (nev >= A->cmap->nglobal - 1)
	nev = A->cmap->nglobal - 2;

    if (nev <= 0) {
	*nit = 0;
	return 0;
    }

    /*if (mode == REGULAR_INVERT && ncv < nev + 40)*/
	ncv = nev + 40;
    if (cmdline_ncv > 0)
	ncv = cmdline_ncv;
    if (ncv <= nev)
	ncv = nev + 1;
    if (ncv > A->cmap->nglobal) {
	ncv = A->cmap->nglobal;
    }
    lworkl = ncv * (ncv + 8);

    resid = phgAlloc(nloc * sizeof(*resid));
    v = phgAlloc(nloc * ncv * sizeof(*v));
    workd = phgAlloc(3 * nloc * sizeof(*workd));
    workl = phgAlloc((lworkl + ncv) * sizeof(*workl));
    bounds = workl + lworkl;

    select = phgAlloc(ncv * sizeof(*select));

    /* set up initial residual vector */
#if 1
    if (vec0 != NULL) {
	double *lambda = phgAlloc(vec0->nvec * sizeof(*lambda));
	/* use a random linear combination of input vectors */
	phgInfo(1, "PARPACK: constructing initial residual vector\n");
	info = 1;   /* tell ARPACK to use the initial residual vector */
	if (A->cmap->rank == 0) {
	    for (it = 0; it < vec0->nvec; it++)
		lambda[it] = drand48();
	}
#if USE_MPI
	MPI_Bcast(lambda, vec0->nvec, MPI_DOUBLE, 0, A->cmap->comm);
#endif	/* USE_MPI */ 
	bzero(resid, nloc * sizeof(*resid));
	for (it = 0; it < vec0->nvec; it++) {
	    res_max = lambda[it];
	    for (i = 0; i < nloc; i++)
		resid[i] += res_max * vec0->data[it * nloc + i];
	}
	phgFree(lambda);
    }
#endif
    phgVecDestroy(&vec0);

    IPARAM(1) = 1;		    /* use exact ishfts */
    IPARAM(3) = itmax * nev;	    /* maxitr */
    phgOptionsPush();
    phgSolverSetDefaultSuboptions();
    phgOptionsSetKeyword("-solver", phgSolverNames[solver_id]);
    phgOptionsSetInt("-solver_maxit", solver_maxit);
    phgOptionsSetFloat("-solver_rtol", solver_rtol);
    phgOptionsSetOptions(solver_opts);
    if (mode == SHIFT_INVERT) {
	/* shift-invert mode: IPARAM(7) = 3 */
	whch = (which <= 0 ? "LM" : "SM");
	IPARAM(7) = 3;
	if (tau != 0.0) {
	    MAT *T = NULL;
	    phgMatCopy(A, &T);
	    phgMatAXPBY(-tau, B, 1.0, &T);
	    solver = phgMat2Solver(SOLVER_DEFAULT, T);
	    phgMatDestroy(&T);
	}
	else {
	    solver = phgMat2Solver(SOLVER_DEFAULT, A);
	}
    }
    else {
	whch = (which <= 0 ? "SM" : "LM");
	/* regular-invert mode: IPARAM(7) = 1 or 2 */
	/* See dsdrv2.f if bmat='I', dsdrv4.f if bmat='G' */
	if (B == NULL) {
	    IPARAM(7) = 1;  /* regular mode (standard problem) */
	}
	else {
	    IPARAM(7) = 2;  /* regular invert mode (generalized problem) */
	    solver = phgMat2Solver(SOLVER_DEFAULT, B);
	    phgFree(solver->rhs->data);
	    solver->rhs->data = NULL;
	}
    }
    phgOptionsPop();

    /* Note: ssd, if not NULL, must point to an user pc_proc function */
    if (solver != NULL && ssd != NULL)
	phgSolverSetPC(solver, solver, ssd);

    /* repeatedly call PARPACK's pdsaupd subroutine */
    if (phgVerbosity > 0)
	phgPrintf("PARPACK: starting Arnoldi update iterations.\n");
    for (it = 1; flag; ) {
	F77_FUNC(pdsaupd,PDSAUPD)(
			&comm, &ido, bmat, &nloc, whch, &nev,
			&tol, resid, &ncv, v, &ldv, iparam,
			ipntr, workd, workl, &lworkl, &info, 1, 2);

	/* Note: it seems that in regular-invert mode, when computing
		 y <-- inv[M]*A*x, the vector x should be overwritten with
		 A*x, and this is not clearly documented in the
		 documentation (see, e.g., dsdrv3.f) */

	/* IPNTR(4) points to iter number (NOTE: NEED PATCHED pdsaup2.f) */
	i = (int)(WORKL(IPNTR(4)));
	if (phgVerbosity > 0 && it < i) {
	    it = i;
	    /* get the ncv error bounds (scaled by the Ritz values) */
	    for (i = 0; i < ncv; i++) {
		res_max = fabs(WORKL(IPNTR(6) + i));
		if (res_max == 0.0)
		    bounds[i] = WORKL(IPNTR(7) + i);
		else
		    bounds[i] = WORKL(IPNTR(7) + i) / res_max;
	    }
	    /* sort the error bounds */
	    qsort(bounds, ncv, sizeof(*bounds), comp_double);
	    phgPrintf("    iteration %d, estimate = %0.4le\n", it - 1,
			    (double)bounds[nev]);
	}
	switch (ido) {
	    case -1:
		/* Perform y <-- OP(x), where OP(x) is
		    inv[A-sigma*M]*M*x for shift-invert mode, or
		    inv[M]*A*x for regular-invert mode,
		    with x = WORKD(IPNTR(1)), y = WORKD(IPNTR(2))
		 */
		if (mode == SHIFT_INVERT)
		    mat_vec(B, solver, &WORKD(IPNTR(1)), &WORKD(IPNTR(2)));
		else
		    mat_vec(A, solver, &WORKD(IPNTR(1)), &WORKD(IPNTR(2)));
		break;
	    case 1:
		/* Perform y <-- OP(x), where OP(x) is
		    inv[A-sigma*M]*M*x for shift-invert mode, or
		    inv[M]*A*x for regular-invert mode,
		    with x = WORKD(IPNTR(1)), y = WORKD(IPNTR(2))
		    Note: M*x is already stored at WORKD(IPNTR(3))
		    and need not be recomputed for shift-invert mode
		 */
		if (mode == SHIFT_INVERT)
		    mat_vec(NULL,solver,&WORKD(IPNTR(3)),&WORKD(IPNTR(2)));
		else
		    mat_vec(A,solver,&WORKD(IPNTR(1)),&WORKD(IPNTR(2)));
		break;
	    case 2:
		/* Perform y <-- M*x with
		    x = WORKD(IPNTR(1)), y = WORKD(IPNTR(2)) */
		mat_vec(B, NULL, &WORKD(IPNTR(1)), &WORKD(IPNTR(2)));
		break;
	    case 3:
		phgError(1, "%s:%d, unexpected.\n", __FILE__, __LINE__);
		break;
	    default:
		/* Terminate the loop */
		assert(ido == 99);
		flag = FALSE;
	}
    }

    /* destroy linear solver */
    if (solver != NULL)
	phgSolverDestroy(&solver);

    if (info < 0)
	phgError(1, "error in PDSAUPD, info = %d\n", info);

    if (keep_unconverged && IPARAM(5) < n) {
	phgPrintf("WARNING: only %d eigenvalue%s converged.\n",
			IPARAM(5), IPARAM(5) < 2 ? " has" : "s have");
	/* Hacky: try to fool PARPACK on the # of converged eigenvalues */
	tol = 1e10;
	IPARAM(5) = n;
    }

    /* Note (from ./PARPACK/SRC/MPI/pdsaupd.f):
     *     -------------------------------------------------------------
     *     IPNTR(1): pointer to the current operand vector X in WORKD.
     *     IPNTR(2): pointer to the current result vector Y in WORKD.
     *     IPNTR(3): pointer to the vector B * X in WORKD when used in
     *               the shift-and-invert mode.
     *     IPNTR(4): pointer to the next available location in WORKL
     *               that is untouched by the program.
     *     IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
     *     IPNTR(6): pointer to the NCV RITZ values array in WORKL.
     *     IPNTR(7): pointer to the Ritz estimates in array WORKL associated
     *               with the Ritz values located in RITZ in WORKL.
     *     IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
     *
     *     Note: IPNTR(8:10) is only referenced by pdseupd . See Remark 2.
     *     IPNTR(8): pointer to the NCV RITZ values of the original system.
     *     IPNTR(9): pointer to the NCV corresponding error bounds.
     *     IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
     *                of the tridiagonal matrix T. Only referenced by
     *                pdseupd  if RVEC = .TRUE. See Remarks. */

    F77_FUNC(pdseupd,PDSEUPD)(
		&comm, &rvec, howmny, select, evals, v, &ldv, &sigma,
		bmat, &nloc, whch, &nev, &tol, resid, &ncv, v, &ldv,
		iparam, ipntr, workd, workl, &lworkl, &info, 1, 1, 2);
    if (info != 0)
	phgError(1, "error in PDSEUPD, info = %d\n", info);

    if (phgVerbosity > 0)
	phgPrintf("mat-vec: %d, linsolve: %d\n", matvec_count, solve_count);

    nev = IPARAM(5);		    /* number of converged eigenvalues */
    if (n > nev)
	n = nev;
    memcpy(vec->data, v, nloc * n * sizeof(*v));
    *nit = IPARAM(3);	    /* number of Arnoldi update iterations taken */

    /* Note: IPNTR(9) points to the NCV error bounds (in workl) */
    if (phgVerbosity > 0) {
	for (i = 0, res_max = 0.; i < n; i++) {
	    tol = fabs(WORKL(IPNTR(8) + i));
	    if (tol == 0.0)
		tol = WORKL(IPNTR(9) + i);
	    else
		tol = WORKL(IPNTR(9) + i) / tol;
	    if (phgVerbosity > 1)
		phgPrintf("    residual[%d] = %lg\n", i, (double)tol);
	    if (res_max < tol)
		res_max = tol;
	}
	phgPrintf("PARPACK: %d iterations, unscaled residual = %0.4le\n",
		    it, (double)res_max);
    }

    phgFree(select);
    phgFree(resid);
    phgFree(v);
    phgFree(workd);
    phgFree(workl);

    phgVecDestroy(&vec0);
    phgVecDestroy(&va);
    phgVecDestroy(&vb);

    return n;
}

#else	/* USE_ARPACK */
EigenSolverVoid(phgEigenSolverARPACK)
#endif	/* USE_ARPACK */
