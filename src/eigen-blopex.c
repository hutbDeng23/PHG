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

/* $Id: eigen-blopex.c,v 1.67 2022/09/16 22:41:49 zlb Exp $ */

#include "phg.h"

#if USE_BLOPEX

#include <blopex.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

/* define some macros for older versions of BLOPEX for convenience */
#ifndef BLOPEX_VERSION
# define BLOPEX_VERSION		0.0
# define BLOPEX_VERSION_MAJOR	0
# define BLOPEX_VERSION_MINOR	0
# define BLOPEX_RELEASE		0
# define BlopexInt		int
# define FLOATP			FLOAT *
#else	/* !defined(BLOPEX_VERSION) */
# define FLOATP			void *
#endif	/* !defined(BLOPEX_VERSION) */

/* Prototypes of needed LAPACK functions */
extern int F77_FUNC(dsygv,DSYGV)(int *itype, char *jobz, char *uplo, int *
                    n, FLOAT *a, int *lda, FLOAT *b, int *ldb,
                    FLOAT *w, FLOAT *work, int *lwork, int *info);
extern int F77_FUNC(dpotrf,DPOTRF)(char *uplo, int *n, FLOAT *a, int * lda,
		    int *info);

/*------------------------ Interface functions ------------------------*/

#if 0
# define CHECK(x, y) \
  check(__func__, nlocal,  x->nact, x->list, x->data, y->nact, y->list, y->data)
static void
check(const char *func, int nlocal,  int nx, int *xlist, FLOAT *xdata,
				int ny, int *ylist, FLOAT *ydata)
{
  FLOAT sx = 0., sy = 0., *p;
  int i, j;

  for (i = 0; i < nx; i++) {
    p = xdata + xlist[i] * nlocal;
    for (j = 0; j < nlocal; j++) {
      sx += *(p++);
    }
  }
  for (i = 0; i < ny; i++) {
    p = ydata + ylist[i] * nlocal;
    for (j = 0; j < nlocal; j++) {
      sy += *(p++);
    }
  }
  printf("%-25s nx=%d, %12.5le, ny=%d, %12.5le\n",
		func, nx, (double)sx, ny, (double)sy);
}
#else
# define CHECK(x, y)
#endif

static INT nlocal;
static int rank, nprocs;
#if USE_MPI
static MPI_Comm comm;
#endif	/* USE_MPI */

typedef struct {
    FLOAT	*data;
    int		*list;
    int		nvec, nact;	/* number of total/active vectors */
} MULTI_VEC;

/*---------------------------------------------------------------------------*/
/* Command-line options variables */

/* The global PC */
static int pc_id = 1;
static char *pc_opts = NULL;
static FLOAT pc_tau = 0.0;

/* The local PC (pc_id1 == 0 disables the local PC) */
static int pc_id1 = 0;
static char *pc_opts1 = NULL;

/* Other options */
static INT verbosity = 0;
/*---------------------------------------------------------------------------*/

static void *
#if BLOPEX_VERSION_MAJOR == 0
CreateMultiVector(BlopexInt nvec)
#else	/* BLOPEX_VERSION_MAJOR == 0 */
CreateMultiVector(void *InterfaceInterpreter, BlopexInt nvec, void *sample)
#endif	/* BLOPEX_VERSION_MAJOR == 0 */
{
  MULTI_VEC *x0 = phgAlloc(sizeof(*x0));
  BlopexInt i;

  x0->nact = x0->nvec = nvec;
  x0->data = phgAlloc(nvec * nlocal * sizeof(*x0->data));
  x0->list = phgAlloc(nvec * sizeof(*x0->list));
  for (i = 0; i < nvec; i++)
    x0->list[i] = i;
  return x0;
}

static void *
CreateCopyMultiVector(void *src_, BlopexInt copyValues)
{
   MULTI_VEC *src = src_;
   MULTI_VEC *dest;
   BlopexInt i, nvec;
  
   dest = phgAlloc(sizeof(*dest));
   dest->nact = dest->nvec = nvec = src->nvec;
   dest->data = phgAlloc(nvec * nlocal * sizeof(*dest->data));
   dest->list = phgAlloc(nvec * sizeof(*dest->list));
   for (i = 0; i < nvec; i++)
	dest->list[i] = i;
   if (copyValues)
      memcpy(dest->data, src->data, nvec * nlocal * sizeof(*dest));

   return dest;
}

static void
DestroyMultiVector(void *vec_)
{
  MULTI_VEC *vec = vec_;

  phgFree(vec->data);
  phgFree(vec->list);
  phgFree(vec);
}

static BlopexInt
MultiVectorWidth(void *v)
{
  return ((MULTI_VEC *)v)->nvec;
}

static void
MultiSetMask(void *vec_, BlopexInt *mask)
{
  MULTI_VEC *vec = vec_;
  BlopexInt i, m = 0;

  for (i = 0; i < vec->nvec; i++)
    if (mask == NULL || mask[i] != 0)
      vec->list[m++] = i;
  vec->nact = m;
}

static void
CopyMultiVector(void *x_, void *y_)
{
  MULTI_VEC *x = x_, *y = y_;
  BlopexInt i;

  assert(x->nact == y->nact);
  for (i = 0; i < x->nact; i++)
    memcpy(y->data + y->list[i] * nlocal, x->data + x->list[i] * nlocal,
		nlocal * sizeof(*x->data));
  CHECK(x, y);
}

static void
ClearMultiVector(void *x_)
{
  MULTI_VEC *x = x_;
  FLOAT *p;
  BlopexInt i;
  INT j;

  for (i = 0; i < x->nact; i++) {
    p = x->data + x->list[i] * nlocal;
    for (j = 0; j < nlocal; j++)
      *(p++) = 0.0;
  }
}

static void
SetMultiVectorRandomValues(void *v_, BlopexInt seed)
{
  MULTI_VEC *v = v_;
  FLOAT *start, *end;
  BlopexInt i;

  srand48(seed + rank);
  for (i = 0; i < v->nact; i++) {
    start = v->data + v->list[i] * nlocal;
    end = start + nlocal;
    while (start < end) {
       *(start++) = 2.0 * drand48() - 1.0;
    }
  }
}

static void
MultiInnerProd(void *x_, void *y_, BlopexInt gh, BlopexInt h, BlopexInt w,
	       FLOATP v_)
{
  FLOAT *v = v_;
  MULTI_VEC *x = x_, *y = y_;
  BlopexInt i, j;
  INT k;
  FLOAT d, *p, *q, *r, *v0;

  assert(x->nact == h && y->nact == w);
  r = v0 = phgAlloc(h * w * sizeof(*v0));
  for (j = 0; j < w; j++) {
    p = y->data + y->list[j] * nlocal;
    for (i = 0; i < h; i++) {
      q = x->data + x->list[i] * nlocal;
      d = 0.0;
      for (k = 0; k < nlocal; k++)
	d += p[k] * *(q++);
      *(r++) = d;
    }
  }

#if USE_MPI
  if (nprocs > 1) {
    r = phgAlloc(h * w * sizeof(*r));
    MPI_Allreduce(v0, r, w * h, PHG_MPI_FLOAT, PHG_SUM, comm);
    phgFree(v0);
    v0 = r;
  }
#endif	/* USE_MPI */
  for (j = 0, r = v0; j < w; j++, r += h, v += gh)
    memcpy(v, r, h * sizeof(*r));
  phgFree(v0);
  CHECK(x, y);
}

static void
MultiInnerProdDiag(void *x_, void *y_, BlopexInt *mask, BlopexInt n0,
		   FLOATP diag_)
{
  FLOAT *diag = diag_;
  MULTI_VEC *x = x_, *y = y_;
  BlopexInt *list = phgAlloc(n0 * sizeof(*list));
  BlopexInt i;
  INT k;
  FLOAT d, *p, *q, *r, *v0;

  assert(x->nact == y->nact && x->nact <= n0);

  for (i = 0, k = 0; i < x->nvec; i++)
    if (mask == NULL || mask[i] != 0)
      list[k++] = i;
  assert (k == x->nact);

  r = v0 = phgAlloc(x->nact * sizeof(*v0));
  for (i = 0; i < x->nact; i++) {
    p = y->data + y->list[i] * nlocal;
    q = x->data + x->list[i] * nlocal;
    d = 0.0;
    for (k = 0; k < nlocal; k++)
      d += *(p++) * *(q++);
    *(r++) = d;
  }

#if USE_MPI
  if (nprocs > 1) {
    r = phgAlloc(x->nact * sizeof(*r));
    MPI_Allreduce(v0, r, x->nact, PHG_MPI_FLOAT, PHG_SUM, comm);
    phgFree(v0);
    v0 = r;
  }
#endif	/* USE_MPI */

  for (i = 0, r = v0; i < x->nact; i++)
    diag[list[i]] = *(r++);

  phgFree(v0);
  phgFree(list);
  CHECK(x, y);
}

static void
MultiVectorByDiagonal(void *x_, BlopexInt *mask, BlopexInt n0, FLOATP alpha_,
		      void *y_)
/* y(<y_mask>) = alpha(<mask>) .* x(<x_mask>) */
{
  FLOAT *alpha = alpha_;
  MULTI_VEC *x = x_, *y = y_;
  BlopexInt *list = phgAlloc(n0 * sizeof(*list));
  BlopexInt i, j;
  INT k;
  FLOAT d, *p, *q;

  assert(x->nact == y->nact && x->nact <= n0);

  for (i = 0, j = 0; i < x->nvec; i++)
    if (mask == NULL || mask[i] != 0)
      list[j++] = i;
  assert (j == x->nact);

  for (i = 0; i < x->nact; i++) {
    p = y->data + y->list[i] * nlocal;
    q = x->data + x->list[i] * nlocal;
    d = alpha[list[i]];
    for (k = 0; k < nlocal; k++)
      *(p++) = *(q++) * d;
  }

  phgFree(list);
  CHECK(x, y);
}

static void 
MultiVectorByMatrix(void *x_, BlopexInt gh, BlopexInt h, BlopexInt w,
		    FLOATP v_, void *y_)
{
  FLOAT *v = v_;
  MULTI_VEC *x = x_, *y = y_;
  BlopexInt i, j;
  INT k;
  FLOAT d, *p, *q;

  assert (h == x->nact && w == y->nact);
 
  for (j = 0; j < w; j++) {
    p = y->data + y->list[j] * nlocal;
    q = x->data + x->list[0] * nlocal;
    d = *(v++);
    for (k = 0; k < nlocal; k++)
      p[k] = d * *(q++);
    for (i = 1; i < h; i++) {
      q = x->data + x->list[i] * nlocal;
      d = *(v++);
      for (k = 0; k < nlocal; k++)
        p[k] += d * *(q++);
    }
    v += gh - h;
  }

  CHECK(x, y);
}

static void
MultiVectorAxpy(FLOAT alpha, void *x_, void *y_)
{
  MULTI_VEC *x = x_, *y = y_;
  BlopexInt i;
  INT j;
  FLOAT *p, *q;

  assert (x->nact == y->nact);
  for (i = 0; i < x->nact; i++) {
    p = y->data + y->list[i] * nlocal;
    q = x->data + x->list[i] * nlocal;
    for (j = 0; j < nlocal; j++)
      *(p++) += alpha * *(q++);
  }
  CHECK(x, y);
}

static VEC *va = NULL, *vb = NULL;
static SOLVER *solver = NULL;		/* The global preconditioner */
static SOLVER *solver1 = NULL;		/* The local preconditioner */
static MAT *matA, *matB;
static double *eigen_values;

static void
setup_solver1(FLOAT pc_tau1)
{
    if (pc_id1 <= 0)
	return;

    if (solver1 != NULL && !solver1->oem_solver->iterative)
	phgSolverDestroy(&solver1);

    if (solver1 == NULL) {
	MAT *T = NULL;
	phgOptionsPush();
	phgOptionsSetKeyword("-solver", phgSolverNames_[pc_id1]);
	phgOptionsSetInt("-verbosity", (INT)0);
	phgOptionsSetInt("-solver_maxit", 10);
	phgOptionsSetFloat("-solver_atol", (FLOAT)0.);
	phgOptionsSetFloat("-solver_btol", (FLOAT)0.);
	phgOptionsSetFloat("-solver_rtol", (FLOAT)0.1);
	phgOptionsSetNoArg("-solver_warn_maxit", FALSE);
	phgOptionsSetOptions(pc_opts1);
	phgMatAXPBY_(1.0, matA, 0.0, &T, TRUE /* force T to be matrix-free */);
	phgMatAXPBY(-1.0 /* must be nonzero here */, matB, 1.0, &T);
	solver1 = phgMat2Solver(SOLVER_DEFAULT, T);
	phgMatDestroy(&T);
	phgFree(solver1->rhs->data);
	solver1->rhs->data = NULL;
	phgOptionsPop();
	return;
    }

    assert(solver1->oem_solver->iterative);

    /* The next lines are tricky and depend on the creation of solver1 */
    assert(solver1->mat->type == PHG_MATRIX_FREE);
    assert(solver1->mat->mv_x == matB);
    solver1->mat->mv_a = -pc_tau1;
    phgFree(solver1->mat->diag);
    solver1->mat->diag = NULL;
}

static void
MatMultiVec(void *data, void *x_, void *y_)
{
  MULTI_VEC *x = x_, *y = y_;
  MAT *mat = data;
  BlopexInt i;

  assert (x->nact == y->nact);
  for (i = 0; i < x->nact; i++) {
    va->data = x->data + x->list[i] * nlocal;
    vb->data =  y->data + y->list[i] * nlocal;
    if (mat != matA && mat != matB) { /* apply the preconditioner */
	assert(data == solver);

	if (solver == NULL && pc_id1 <= 0) {
	    memcpy(vb->data, va->data, nlocal * sizeof(*vb->data));
	    continue;
	}

#if 0
	bzero(vb->data, nlocal * sizeof(*vb->data));
#else
	memcpy(vb->data, va->data, nlocal * sizeof(*vb->data));
#endif

	/* Apply the global preconditioner (A - pc_tau * B)^(-1) */
	if (solver != NULL) {
	    solver->rhs->data = va->data;
	    solver->rhs->assembled = TRUE;
	    phgSolverVecSolve(solver, FALSE, vb);
	    solver->rhs->data = NULL;
	}

	/* Apply the local preconditioner (A - eigen_values[i] * B)^(-1) */
	if (pc_id1 > 0) {
#if 0
static FLOAT eigen_values[] = {
    2.961001761849e+01,
    5.922810928096e+01,
    5.922810928096e+01,
    5.922810928096e+01,
    8.886094398911e+01,
    8.886094398912e+01,
    8.886094398912e+01,
    1.086317278746e+02,
    1.086317278746e+02,
    1.086317278746e+02
};
#endif
	    setup_solver1(eigen_values[i]);
	    solver1->rhs->data = va->data;
	    solver1->rhs->assembled = TRUE;
	    phgSolverVecSolve(solver1, FALSE, vb);
	    solver1->rhs->data = NULL;
	}

	continue;
    }
    phgMatVec(MAT_OP_N, 1.0, mat, va, 0.0, &vb);
    va->data = NULL;
    vb->data = NULL;
  }

  CHECK(x, y);
}

static void
initInterfaceInterpreter(mv_InterfaceInterpreter *ii)
{
    /* Vector part */
    ii->CreateVector = NULL;
    ii->DestroyVector = NULL;
    ii->InnerProd = NULL; 
    ii->CopyVector = NULL;
    ii->ClearVector = NULL;
    ii->SetRandomValues = NULL;
    ii->ScaleVector = NULL;
    ii->Axpy = NULL;

    /* Multivector part */
    ii->CopyCreateMultiVector = CreateCopyMultiVector;
    ii->DestroyMultiVector = DestroyMultiVector;
    ii->Width = MultiVectorWidth;
    ii->Height = NULL;
    ii->SetMask = MultiSetMask;
    ii->CopyMultiVector = CopyMultiVector;
    ii->ClearMultiVector = ClearMultiVector;
    ii->SetRandomVectors = SetMultiVectorRandomValues;
    ii->MultiInnerProd = MultiInnerProd;
    ii->MultiInnerProdDiag = MultiInnerProdDiag;
    ii->MultiVecMat = MultiVectorByMatrix;
    ii->MultiVecMatDiag = MultiVectorByDiagonal;
    ii->MultiAxpy = MultiVectorAxpy;
    ii->MultiXapy = NULL;
    ii->Eval = NULL;
}

EigenSolverPrototype(phgEigenSolverBLOPEX)
/* (A, B, n, which, tau, evals, evecs, eps, itmax, nit) */
{
    static BOOLEAN initialized = FALSE;
    MULTI_VEC *x0;
    mv_MultiVectorPtr x;
    FLOAT *resid = NULL;
    lobpcg_Tolerance lobpcg_tol;
    mv_InterfaceInterpreter ii;
    lobpcg_BLASLAPACKFunctions blap_fn;
    VEC *vec;
    MAT *m = (A != NULL ? A : B);
    int i, ret;
 
    if (!initialized) {
	/* register options */
	initialized = TRUE;

	phgOptionsRegisterTitle("\nBLOPEX (LOBPCG) options:", "\n", "blopex");

	phgOptionsRegisterKeyword("-blopex_pc_solver",
				"Solver for the global preconditioner",
				phgSolverNames_, &pc_id);
	phgOptionsRegisterString("-blopex_pc_opts",
				"Options for the global preconditioner",
				&pc_opts);
	phgOptionsRegisterFloat("-blopex_pc_tau", "Tau for the preconditioner",
				&pc_tau);

	phgOptionsRegisterKeyword("-blopex_pc_solver1",
				"Solver for the local preconditioner",
				phgSolverNames_, &pc_id1);
	phgOptionsRegisterString("-blopex_pc_opts1",
				"Options for the local preconditioner",
				&pc_opts1);

	phgOptionsRegisterInt("-blopex_verbosity", "BLOPEX verbosity level",
				&verbosity);
	return 0;
    }

    assert(FT_PHG == FT_DOUBLE);

    *nit = 0;
    vec = *evecs;
    MagicCheck(VEC, vec);

    if (n <= 0 || A->cmap->nglobal <= 0)
	return 0;

    if (which > 0) {
	MAT *m = A;
	A = B;
	B = m;
    }

    if (n >= m->cmap->nglobal)
	n = m->cmap->nglobal;

    nlocal = m->cmap->nlocal;
    rank = m->cmap->rank;
    nprocs = m->cmap->nprocs;
#if USE_MPI
    comm = m->cmap->comm;
#endif

    matA = A;
    matB = B;
    if (tau != 0.0) {
#if 0
	matA = NULL;
	phgMatCopy(A, &matA);
	phgMatAXPBY(-tau, B, 1.0, &matA);
#else
	if (m->cmap->rank == 0)
            phgWarning("%s: the value of 'tau' (%lg) is ignored.\n",
						__func__, (double)tau);
#endif
    }

    initInterfaceInterpreter(&ii);
#if BLOPEX_VERSION_MAJOR == 0
    x0 = CreateMultiVector(n);
#else	/* BLOPEX_VERSION_MAJOR == 0 */
    x0 = CreateMultiVector(&ii, n, NULL);
#endif	/* BLOPEX_VERSION_MAJOR == 0 */
    /* FIXME: this causes a memory leak of 24 bytes, to be investigated */
    x = mv_MultiVectorWrap(&ii, x0, 0/*ownsData, FIXME: what's its meaning?*/);

    if (vec == NULL) {
	*evecs = vec = phgMapCreateVec(m->cmap, n);
	vec->nvec = 0;
    }
    if (vec->nvec != n) {
	phgFree(vec->data);
	SetMultiVectorRandomValues(x0, 1);
	vec->nvec = n;
    }
    else {
	phgFree(x0->data);
	x0->data = vec->data;
    }
    vec->data = NULL;

    lobpcg_tol.absolute = eps;
    lobpcg_tol.relative = 1e-50;

    assert(sizeof(BLAS_INT) == sizeof(int));
    blap_fn.dpotrf = F77_FUNC(dpotrf,DPOTRF);
    blap_fn.dsygv = F77_FUNC(dsygv,DSYGV);

    va = phgMapCreateVec(m->cmap, 1);
    vb = phgMapCreateVec(m->cmap, 1);
    phgFree(va->data);
    va->data = NULL;
    phgFree(vb->data);
    vb->data = NULL;    

    resid = phgAlloc(n * sizeof(*resid));

    /* Setup the global PC */
    if (pc_id > 0) {
	phgOptionsPush();
	phgSolverSetDefaultSuboptions();
	phgOptionsSetKeyword("-solver", phgSolverNames_[pc_id]);
	phgOptionsSetInt("-solver_maxit", 1);
	phgOptionsSetFloat("-solver_atol", (FLOAT)0.);
	phgOptionsSetFloat("-solver_btol", (FLOAT)0.);
	phgOptionsSetFloat("-solver_rtol", (FLOAT)0.);
	phgOptionsSetNoArg("-solver_warn_maxit", FALSE);
	phgOptionsSetOptions(pc_opts);
	if (pc_tau == 0.0) {
	    solver = phgMat2Solver(SOLVER_DEFAULT, matA);
	}
	else {
	    MAT *T = NULL;
	    phgMatCopy(matA, &T);
	    phgMatAXPBY(-pc_tau, B, 1.0, &T);
	    solver = phgMat2Solver(SOLVER_DEFAULT, T);
	    phgMatDestroy(&T);
	}
	/* ssd points to an user pc_proc function (subject to change!) */
	if (ssd != NULL)
	    phgSolverSetPC(solver, solver, ssd);
	phgFree(solver->rhs->data);
	solver->rhs->data = NULL;
	phgOptionsPop();
    }

    assert(FT_PHG == FT_DOUBLE);
    eigen_values = evals;	/* for referencing evals[] in MatMultiVec() */
    ret =
#if BLOPEX_VERSION_MAJOR == 0
	lobpcg_solve
#else	/* BLOPEX_VERSION_MAJOR == 0 */
	lobpcg_solve_double
#endif	/* BLOPEX_VERSION_MAJOR == 0 */
	       (x,		/* eigenvectors */
		matA,		/* A data */
		MatMultiVec,	/* A*x */
		B,		/* B data */
		MatMultiVec,	/* B*x */
		solver,		/* T data */
		MatMultiVec,	/* T*x (precon matrix: (A - tau * B)^(-1)?) */
		NULL,		/* blockVectorY (constraints) */
		blap_fn,	/* BLAS/LAPACK functions */
		lobpcg_tol,	/* tolerances */
		itmax * n,	/* MAXIT */
		verbosity,	/* verbose level */
		nit,		/* iteration nunmber */
		evals,		/* lambda's */
		NULL,		/* lambda history: neig * (it + 1) */
		0,		/* lamda history LDA */
		resid,		/* residual norms */
		NULL,		/* residual norms history: neig * (it + 1) */
		0);		/* residual norms history LDA */
    if (solver != NULL)
	phgSolverDestroy(&solver);
    if (solver1 != NULL)
	phgSolverDestroy(&solver1);

    phgFree(resid);

    /* Note: possible return values of lobpcg_solve() are defined in the file
     * 'blopex_abstract/krylov/lobpcg.c':
     *		 1: problem size is too small,
     *		 2: block size < 1,            
     *		 3: linearly dependent constraints,
     *		-1: requested accuracy not achieved. */
#if 0
    if (ret < 0) {
	for (i = 0; i < n; i++)
	    phgPrintf("%d: resid = %le, evals = %lg\n",
				i, (double)resid[i], (double)evals[i]);
    }
#endif
    if (ret != 0 || *nit <= 0 || *nit > itmax * n + 1)
	n = 0;

    vec->data = x0->data;
    x0->data = NULL;
    DestroyMultiVector(x0);

    phgVecDestroy(&va);
    phgVecDestroy(&vb);

    if (which > 0) {
	for (i = 0; i < n; i++)
	    evals[i] = 1.0 / evals[i];
    }

    return n;
}

#else	/* USE_BLOPEX */
EigenSolverVoid(phgEigenSolverBLOPEX)
#endif	/* USE_BLOPEX */
