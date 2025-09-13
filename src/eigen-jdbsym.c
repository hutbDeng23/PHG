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

/* $Id: eigen-jdbsym.c,v 1.50 2011/06/08 06:31:34 zlb Exp $ */

#include "phg.h"

#if USE_JDBSYM

#include <stdlib.h>
#include <jdbsym.h>

static FLOAT *diag = NULL;

static MAT *ma, *mb;
static VEC *va, *vb;

static void
mat_vec(MAT *a, FLOAT *x, FLOAT *y)
/* computes y := Ax */
{
    va->data = x;
    vb->data = y;
    phgMatVec(MAT_OP_N, 1.0, a, va, 0.0, &vb);
    va->data = NULL;
    vb->data = NULL;
}

static void
mat_vec_a(FLOAT *x, FLOAT *y)
/* computes y := Ax */
{
    mat_vec(ma, x, y);
}

static void
mat_vec_b(FLOAT *x, FLOAT *y)
/* computes y := Bx */
{
    mat_vec(mb, x, y);
}

/* diagonal preconditioner: inverse diag of A - tau*B */
static void
precon_diag(FLOAT *x, FLOAT *y)
{ 
    INT i, nlocal;
    nlocal = ma->cmap->nlocal;

    assert(diag != NULL);

    for (i = 0; i < nlocal; i++)
	y[i] = diag[i] * x[i];
}

#if USE_JDBSYM_PARALLEL && USE_MPI
static void
gsum(FLOAT *sbuf, FLOAT *rbuf, int count, void *user_data)
{ 
    (void)MPI_Allreduce(sbuf, rbuf, count, PHG_MPI_FLOAT, PHG_SUM,
						*(MPI_Comm *)user_data);
}
#endif	/* USE_JDBSYM_PARALLEL && USE_MPI */

EigenSolverPrototype(phgEigenSolverJDBSYM)
/* (A, B, n, which, tau, evals, evecs, eps, itmax, nit) */
{
    static BOOLEAN initialized = FALSE;
    INT i, n0, nlocal_min;
    VEC *vec;
    void *prec = NULL;

    /* JDBSYM parameters */
    static FLOAT eps_tr = 1e-4;
    static FLOAT toldecay = 2.0;
    static INT blksize0 = 1;
    static INT blkwise = 0;
    static INT mvmax = 200;
    static INT jmax0 = -1;
    static INT jmin0 = -1;

    int clvl = phgVerbosity > 0 ? phgVerbosity + 1 : 0;	/* verbose level */
    int kmax, jmax, jmin, blksize, k_conv = 0, v0dim = 0; 
    FLOAT *v0 = NULL;

    static const char *linsolver_names[] = {"cgs",
#if USE_JDBSYM_PARALLEL
					    "cg",
#endif
					    "minres", "symmlq", "qmr", NULL};
    static int linsolver_list[] = {LINSOLV_CGS,
#if USE_JDBSYM_PARALLEL
				   LINSOLV_CG,
#endif
				   LINSOLV_MINRES, LINSOLV_SYMMLQ, LINSOLV_QMR};
    static int linsolver = 0;

    static const char *precon_names[] = {"none", "diag", NULL};
    static void (*precon_list[])(FLOAT *x, FLOAT *y) = {NULL, precon_diag};
    static int precon = 0;

    static const char *optype_names[] = {"unsym1", "unsym2", "sym", NULL};
    static int optype_list[] = {OP_UNSYM1, OP_UNSYM2, OP_SYM};
    static int optype = 0;

    static const char *strategy_names[] = {"fixed_tau", "variable_tau", NULL};
    static int strategy = 0;

    assert(FT_PHG == FT_DOUBLE);

    if (!initialized) {
	/* register options */
	initialized = TRUE;
	phgOptionsRegisterTitle("\nJDBSYM options:", "\n", "jdbsym");
	phgOptionsRegisterKeyword("jdbsym_precon", "Preconditioner",
					precon_names, &precon);
	phgOptionsRegisterKeyword("jdbsym_linsolver", "Linear solver",
					linsolver_names, &linsolver);
	phgOptionsRegisterKeyword("jdbsym_optype", "Type of linear operator",
					optype_names, &optype);
	phgOptionsRegisterKeyword("jdbsym_strategy", "How to set TAU after "
				  "convergence of an eigenpair",
				  strategy_names, &strategy);
	phgOptionsRegisterFloat("jdbsym_eps_tr", "Tracking parameter", &eps_tr);
	phgOptionsRegisterFloat("jdbsym_toldecay", "Stopping criterion for the"
				" correction equation", &toldecay);
	phgOptionsRegisterInt("jdbsym_blksize","Blocking parameter",&blksize0);
	phgOptionsRegisterInt("jdbsym_blkwise", "Blocking method (0 or 1)",
				&blkwise);
	phgOptionsRegisterInt("jdbsym_mvmax", "Maximum number of iterations "
			      "for solving the correction equation", &mvmax);
	phgOptionsRegisterInt("jdbsym_jmax", "Maximum dimension of search "
			      "subspace (default 2*n)", &jmax0);
	phgOptionsRegisterInt("jdbsym_jmin", "Minimum dimension of search "
			      "subspace (default n)", &jmin0);
	return 0;
    }

    if (n <= 0 || A->cmap->nglobal <= 0)
	return 0;

    /* Note: JDBSYM requires jmax <= nlocal and kmax <= nlocal */
#if USE_MPI
    (void) /* to avoid MPICH warning */
    MPI_Allreduce(&A->cmap->nlocal, &nlocal_min, 1, PHG_MPI_INT, MPI_MIN,
			A->cmap->comm);
#else	/* USE_MPI */
    nlocal_min = A->cmap->nlocal;
#endif	/* USE_MPI */

    if (nlocal_min <= 0)
	return 0;

    if (n > nlocal_min)
	n = nlocal_min;

    kmax = n;
    jmax = (jmax0 <= 0 ? 2 * kmax : jmax0);

    if (jmax > nlocal_min)
	jmax = nlocal_min;

    jmin = (jmin0 <= 0 ? jmax - blksize0 : jmin0);

    if (jmin <= 0)
	jmin = 1;

    blksize = blksize0;
    if (blksize <= 0)
	blksize = 1;
    else if (blksize > jmin)
	blksize = jmin;

    if (A->cmap->nprocs > 1 && linsolver_list[linsolver] != LINSOLV_CGS
#if USE_JDBSYM_PARALLEL
	&& linsolver_list[linsolver] != LINSOLV_CG
#endif
	)
	phgError(1, "%s: linear solver '%s' is not yet parallel.\n",
			linsolver_names[linsolver]);

    vec = *evecs;
    MagicCheck(VEC, vec);

    n0 = A->cmap->partition[A->cmap->rank];
    Unused(n0);

    ma = A;
    mb = B;
    va = phgMapCreateVec(A->cmap, 1);
    vb = phgMapCreateVec(A->cmap, 1);
    phgFree(va->data);
    va->data = NULL;
    phgFree(vb->data);
    vb->data = NULL;

    prec = precon_list[precon];
#if USE_JDBSYM_PARALLEL
    if (prec != NULL) {
	phgWarning("jdbsym-parallel does not support preconditioner!\n");
	prec = NULL;
    }
#endif	/* USE_JDBSYM_PARALLEL */
    if (prec == precon_diag) {
	diag = phgAlloc(A->cmap->nlocal * sizeof(FLOAT));
	phgMatSetupDiagonal(A);
	if (B == NULL) {
	    for (i = 0; i < A->cmap->nlocal; i++)
		diag[i] = 1.0 / (A->diag[i] - tau);
	}
	else {
	    phgMatSetupDiagonal(B);
	    for (i = 0; i < A->cmap->nlocal; i++)
		diag[i] = 1.0 / (A->diag[i] - tau * B->diag[i]);
	}
    } 

    if (vec == NULL) {
	*evecs = vec = phgMapCreateVec(A->cmap, n);
    }
    else {
	v0dim = vec->nvec;
	v0 = vec->data;
	if (v0dim >= jmax)
	    v0dim = jmax - 1;
	vec->nvec = n;
	vec->data = phgAlloc(n * A->cmap->nlocal * sizeof(*vec->data));
    }

    assert(A->rmap->nglobal > 1);

    if (A->cmap->nprocs > 1) {
#if USE_JDBSYM_PARALLEL && USE_MPI
	jdbsym_set_parms(&A->cmap->comm, gsum, A->cmap->nglobal,
			    A->cmap->nprocs, A->cmap->rank);
#else	/* USE_JDBSYM_PARALLEL */
	phgError(1, "This version of JDBSYM requires nprocs == 1.\n");
#endif	/* USE_JDBSYM_PARALLEL */
    }

    *nit = 0;

    jdbsym(A->cmap->nlocal,
       tau, eps, kmax, jmax, jmin, itmax * kmax, blksize, blkwise,
       v0dim, v0, linsolver_list[linsolver], optype_list[optype],
       mvmax, eps_tr, toldecay, strategy, clvl, &k_conv, vec->data, evals,
       nit, mat_vec_a, B == NULL ? NULL : mat_vec_b, prec);

    if (*nit <= 0)
	k_conv = 0;

    phgFree(v0);
    phgFree(diag);
    diag = NULL;
    phgVecDestroy(&va);
    phgVecDestroy(&vb);

    return k_conv;
}

#else	/* USE_JDBSYM */
EigenSolverVoid(phgEigenSolverJDBSYM)
#endif	/* USE_JDBSYM */
