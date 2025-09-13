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

/* $Id: solver-phg.c,v 1.120 2022/09/21 02:35:28 zlb Exp $ */

/* PHG's built-in solvers. */

#include "phg.h"
#include <string.h>
#include <strings.h>
#include <math.h>

#define TRUNK_SIZE 8

#define Initialize		NULL
#define Finalize		NULL
#define Assemble		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL

static int solver_pcg(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc, 
	int max_its, FLOAT rtol, FLOAT btol, FLOAT atol,
	FLOAT *res, BOOLEAN *warned, INT verb);
static int solver_BiCGstab(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc,
	int max_its, FLOAT rtol, FLOAT btol, FLOAT atol, FLOAT *res,
	INT verb);

static void void_pc(void *ctx, VEC *b, VEC **x);
static void jacobi_pc(void *ctx, VEC *b, VEC **x);
static void bjacobi_pc(void *ctx, VEC *b, VEC **x);

static const char *pc_names[] = {
    "none", "jacobi", "bjacobi", "solver", NULL
};
static const PC_PROC pc_list[] = {
    void_pc, jacobi_pc, bjacobi_pc, NULL
};
static int pcg_pc_type = 1;	/* default to Jacobi PC */
static char *pcg_pc_opts = NULL;	/* PC options for PCG */

/* Prototypes of various GMRES solvers */
#define GMRES_FUNC(foo) int \
    foo(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc, int restart, \
	int augk, int max_its, FLOAT rtol, FLOAT btol, FLOAT atol, \
	FLOAT *res, INT verb)
static GMRES_FUNC(solver_GMRES_l);
static GMRES_FUNC(solver_GMRES_r);
static GMRES_FUNC(solver_FGMRES_r);
static GMRES_FUNC(solver_LGMRES_l);
static GMRES_FUNC(solver_LGMRES_r);
static GMRES_FUNC(solver_FLGMRES_r);

typedef GMRES_FUNC((*gmres_func_t));
static const char *gmres_names[] = {
    "gmres_l", "gmres_r", "fgmres_r",
    "lgmres_l", "lgmres_r", "flgmres_r", NULL};
static gmres_func_t gmres_funcs[] = {
    solver_GMRES_l, solver_GMRES_r, solver_FGMRES_r,
    solver_LGMRES_l, solver_LGMRES_r, solver_FLGMRES_r};
static int gmres_type = 5;	/* default to FLGMRES_r */
static int gmres_pc_type = 1;	/* default to Jacobi PC */
static char *gmres_pc_opts = NULL;	/* PC options for GMRES */
static INT gmres_restart = 10;
static INT lgmres_augk = 2;

static int bjacobi_solver_id = -1;
static char *bjacobi_solver_opts = NULL;

static int preonly_pc_type = 1;	/* default to Jacobi PC */
static char *preonly_pc_opts = NULL;	/* PC options for PreOnly solver */

typedef struct {
    SOLVER *bjacobi_solver;	/* local solver for block Jacobi PC */
    char *bjacobi_solver_opts;
    BOOLEAN warned;		/* used to warn non SPD matrix only once */
    int pc_type;
    char *pc_opts;
    INT gmres_restart;
    int bjacobi_solver_id;
    int gmres_type;
    INT lgmres_augk;
} OEM_DATA;

#define _t      ((OEM_DATA *)solver->oem_data)

static FLOAT
get_tol(FLOAT r0_tol, FLOAT b_tol, FLOAT a_tol)
{
    FLOAT tol;

    tol = r0_tol;
    if (tol < b_tol)
	tol = b_tol;
    if (tol < a_tol)
	tol = a_tol;

    return tol;
}

#if 0
static int
solver_Jacobi(MAT *A, VEC *b, VEC **px, int max_its, FLOAT rtol,
		  FLOAT btol, FLOAT atol, FLOAT *res, INT verb)
{
    VEC *x, *r;
    FLOAT residual = 1.0, b_norm, tol;
    INT i, nits;

    if (max_its <= 0)
	return 0;

    if (*px == NULL || b == NULL)
	phgError(1, "%s:%d invalid vector handler\n", __FILE__, __LINE__);

    assert(b->nvec == 1);

    x = *px;
    r = phgMapCreateVec(b->map, 1);
    if (IsZero(b_norm = phgVecNorm2(b, 0, NULL))) {
	bzero(x->data, sizeof(*x->data) * x->map->nlocal);
	return 1;
    }

    if (A->diag == NULL)
	phgMatSetupDiagonal(A);

    nits = 0;
    tol = -1.0;

    while (TRUE) {
	phgVecCopy(b, &r);
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
	for (i = 0; i < A->cmap->nlocal; i++)
	    x->data[i] += r->data[i] / A->diag[i];
	if (++nits >= max_its)
	    break;
	residual = phgVecNorm2(r, 0, NULL);
	if (tol < 0.)
	    tol = get_tol(residual * rtol, b_norm * btol, atol);
	if (residual <= tol)
	    break;
    }				/* end of while-loop */
    phgVecDestroy(&r);
    if (res != NULL)
	*res = residual;

    return nits;
}
#endif

static void
default_pc_proc(void *ctx, VEC *b, VEC **x)
{
    SOLVER *pc_solver = ctx;
    FLOAT *old_rhs;

    phgSolverAssemble(pc_solver);	/* assemble matrix and RHS */

    if (pc_solver->mat->rmap->nlocal != b->map->nlocal) {
	phgPrintf("Invalid PC!\n");
	return;
    }

    /* updating RHS */
    old_rhs = pc_solver->rhs->data;
    pc_solver->rhs->data = b->data;
    pc_solver->rhs->assembled = TRUE;
    pc_solver->rhs_updated = TRUE;
    bzero((*x)->data, (*x)->map->nlocal * sizeof(*(*x)->data));

    /* solving M x = b */
    phgSolverVecSolve(pc_solver, FALSE, *x);
    pc_solver->rhs->data = old_rhs;

    return;
}

static void
void_pc(void *ctx, VEC *b, VEC **x)
{
    SOLVER *pc_solver = ctx;

    phgVecCopy(b, x);
    pc_solver->nits = 0;

    return;
}

static void
jacobi_pc(void *ctx, VEC *b, VEC **x)
{
    SOLVER *pc_solver = ctx;
    MAT *A;
    VEC *v;
    FLOAT *rhs, *p;
    FLOAT value;
    INT i;

    phgSolverAssemble(pc_solver);

#if 0
    /* test the Jacobi solver */
    bzero((*u)->data, B->map->nlocal * sizeof(*(*u)->data));
    /*phgVecCopy(b, x); */
    solver_Jacobi(pc_solver->mat, b, x, 1, 0.0, 1.0e-30, 1.0e-30, NULL, 0);
    return;
#endif

    A = pc_solver->mat;
    if (A->diag == NULL)
	phgMatSetupDiagonal(A);

    v = *x;
    p = v->data;
    rhs = b->data;
    for (i = 0; i < A->rmap->nlocal; i++, rhs++, p++) {
	if ((value = A->diag[i]) == 0.0)
	    phgError(1,
		     "%s:%d, zero diagonal (loc. %" dFMT ", glob. %" dFMT
		     ")\n", __FILE__, __LINE__, i,
		     i + A->rmap->partition[A->rmap->rank]);
	*p = *rhs / value;
    }

    pc_solver->nits = 1;

    return;
}

static void
bjacobi_pc(void *ctx, VEC *b, VEC **x)
/* Submesh-based block jacobi PC */
{
    SOLVER *pc_solver = ctx;
    MAT *A;
    VEC *y;
    SOLVER *solver;

    if (pc_solver->mat->rmap->nlocal == 0)
	return;

    phgSolverAssemble(pc_solver);

    solver = pc_solver;

    A = pc_solver->mat;
    assert(A->type != PHG_DESTROYED);
    if (A->type == PHG_MATRIX_FREE && A->blocks == NULL) {
	/* unimplemented yet */
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
		 __FILE__, __LINE__);
    }

    if (_t->bjacobi_solver == NULL) {
	MAT *B;
	/* set up local solver */
	if (A->rmap->nprocs > 1) {
	    /* create matrix for diagonal block (note: can use phgMatDup) */
	    MAP *map;
	    size_t size;
	    map = phgMapCreateSimpleMap(MPI_COMM_SELF, A->rmap->nlocal,
					A->rmap->nlocal);
	    B = phgMapCreateMat(map, map);

	    /* copy diagonal part of A to B */
	    B->assembled = TRUE;
	    B->handle_bdry_eqns = FALSE;
	    B->nnz_o = 0;
	    B->nnz_d = A->nnz_d;

	    if (A->type == PHG_PACKED) {
		B->type = PHG_PACKED;

		size = B->nnz_d * sizeof(*B->packed_cols);
		B->packed_cols = phgAlloc(size);
		memcpy(B->packed_cols, A->packed_cols, size);

		size = B->nnz_d * sizeof(*B->packed_data);
		B->packed_data = phgAlloc(size);
		memcpy(B->packed_data, A->packed_data, size);

		size = (B->rmap->nlocal + 1) * sizeof(*B->packed_ind);
		B->packed_ind = phgAlloc(size);
		memcpy(B->packed_ind, A->packed_ind, size);
	    }
	    else {
		MAT_ROW *row, *row0;
		INT i, j, k, ncols0, ncols;
		assert(A->type == PHG_UNPACKED);
		if (B->rows == NULL)
		    B->rows = phgCalloc(map->nlocal, sizeof(*row));
		row0 = A->rows;
		row = B->rows;
		for (i = 0; i < map->nlocal; i++, row++, row0++) {
		    ncols0 = row0->ncols;
		    for (ncols = 0, j = 0; j < ncols0; j++) {
			if (row0->cols[j] >= map->nlocal)
			    continue;
			ncols++;
		    }
		    row->ncols = row->alloc = ncols;
		    phgFree(row->cols);
		    row->cols = phgAlloc(ncols * sizeof(*row->cols));
		    phgFree(row->data);
		    row->data = phgAlloc(ncols * sizeof(*row->data));
		    for (k = 0, j = 0; j < ncols0; j++) {
			if (row0->cols[j] >= map->nlocal)
			    continue;
			row->cols[k] = row0->cols[j];
			row->data[k] = row0->data[j];
			if (++k >= ncols)
			    break;
		    }
		}
	    }

	    phgMapDestroy(&map);
	}
	else {
	    B = A;
	}

	phgOptionsPush();
	phgSolverSetDefaultSuboptions();
	/* set some default options for block Jacobi solver */
	phgOptionsSetKeyword("-solver",
			     phgSolverNames[_t->bjacobi_solver_id]);
	phgOptionsSetInt("-solver_maxit", (INT)1);
	phgOptionsSetOptions(_t->bjacobi_solver_opts);
	_t->bjacobi_solver = phgMat2Solver(SOLVER_DEFAULT, B);
	phgOptionsPop();
	_t->bjacobi_solver->warn_maxit = FALSE;

	if (_t->bjacobi_solver->oem_solver == solver->oem_solver) {
	    /* avoid regression */
	    assert(solver->oem_solver == SOLVER_PCG ||
		   solver->oem_solver == SOLVER_BiCGstab ||
		   solver->oem_solver == SOLVER_GMRES);
	    ((OEM_DATA *) _t->bjacobi_solver->oem_data)->pc_type = 0;
	}

	if (B != A)
	    phgMatDestroy(&B);
	phgFree(_t->bjacobi_solver->rhs->data);
	_t->bjacobi_solver->rhs->data = NULL;
    }

    y = phgMapCreateVec(_t->bjacobi_solver->rhs->map, 0);
    y->nvec = (*x)->nvec;
    y->data = (*x)->data;
    _t->bjacobi_solver->rhs->data = b->data;
    y->assembled = _t->bjacobi_solver->rhs->assembled = TRUE;
    _t->bjacobi_solver->rhs_updated = TRUE;
    phgSolverVecSolve(_t->bjacobi_solver, FALSE, y);
    y->data = _t->bjacobi_solver->rhs->data = NULL;
    phgVecDestroy(&y);

    return;
}

/*---------------------------------------------------------------------------*/

void
phgSolverSetPC_(SOLVER *solver, void *ctx, PC_PROC pc_proc)
/* internal function */
{
    SOLVER_PC *pc;

    if (solver->pc == NULL)
	solver->pc = phgCalloc(1, sizeof(SOLVER_PC));
    else if (solver->pc->solver_owned != NULL)
	phgSolverDestroy(&solver->pc->solver_owned);
    pc = solver->pc;
    pc->ctx = ctx;
    pc->pc_proc = (pc_proc == NULL ? default_pc_proc : pc_proc);

    if (pc->pc_proc != NULL && solver->oem_solver->SetPC != NULL)
	solver->oem_solver->SetPC(solver);

    return;
}

void
phgSolverSetPC(SOLVER *solver, void *ctx, PC_PROC pc_proc)
/* user callable function */
{
    assert(solver != NULL);

    if (solver->pc_option_flag) {
	/* commandline options override phgSolverSetPC */
	return;
    }

    phgSolverSetPC_(solver, ctx, pc_proc);
}

/* -------------------------------------------------------------*/

static int
BJacobiRegisterOptions(void)
{
    phgOptionsRegisterTitle("\nThe block Jacobi preconditioner options:",
			    "\n", "bjacobi");

    /* set default solver to the first direct solver if available */
    for (bjacobi_solver_id = 0;
	 phgSolverNames[bjacobi_solver_id] != NULL
	 && (phgSolverList[bjacobi_solver_id]->iterative == TRUE ||
	     phgSolverList[bjacobi_solver_id] == SOLVER_DIAGONAL);
	 bjacobi_solver_id++);
    if (phgSolverNames[bjacobi_solver_id] == NULL)
	bjacobi_solver_id = 0;	/* fall back to the first solver */
    phgOptionsRegisterKeyword("-bjacobi_solver",
			      "Local solver for block Jacobi preconditioner",
			      phgSolverNames, &bjacobi_solver_id);
    phgOptionsRegisterString("-bjacobi_solver_opts",
			     "Options passed to block " "Jacobi solver",
			     &bjacobi_solver_opts);
    return 0;
}

static int
PCGRegisterOptions(void)
{
    BJacobiRegisterOptions();	/* FIXME: put it elsewhere? */

    phgOptionsRegisterTitle("\nThe PCG and BiCGstab solvers options:", "\n",
			    "pcg");
    phgOptionsRegisterKeyword("-pcg_pc_type", "Preconditioner type",
			      pc_names, &pcg_pc_type);
    phgOptionsRegisterString("-pcg_pc_opts", "Options passed to the "
			     "preconditioner", &pcg_pc_opts);
    return 0;
}

static int
GMRESRegisterOptions(void)
{
    phgOptionsRegisterTitle("\nThe GMRES solver options:", "\n",
			    "gmres");
    phgOptionsRegisterKeyword("-gmres_type",
	"Type of GMRES solver ('_l': left precond., '_r': right precond.)",
			      gmres_names, &gmres_type);
    phgOptionsRegisterKeyword("-gmres_pc_type", "Preconditioner type",
			      pc_names, &gmres_pc_type);
    phgOptionsRegisterString("-gmres_pc_opts",
			     "Options passed to the " "preconditioner",
			     &gmres_pc_opts);
    phgOptionsRegisterInt("-gmres_restart", "Restart steps", &gmres_restart);
    phgOptionsRegisterInt("-lgmres_augk",
			   "The parameter augk for LGMRES/FLGMRES",
			   &lgmres_augk);
    return 0;
}

static int
PreOnlyRegisterOptions(void)
{
    phgOptionsRegisterTitle("\nThe PreOnly solver options:", "\n", "preonly");
    phgOptionsRegisterKeyword("-preonly_pc_type", "Preconditioner type",
			      pc_names, &preonly_pc_type);
    phgOptionsRegisterString("-preonly_pc_opts", "Options passed to the "
			     "preconditioner", &preonly_pc_opts);
    return 0;
}

static int
solver_pcg(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc, int max_its,
	       FLOAT rtol, FLOAT btol, FLOAT atol, FLOAT *res,
	       BOOLEAN *warned, INT verb)
{
    int nits;
    FLOAT residual = 0., r0_norm, rb_norm, tol;
    FLOAT rz_1 = 0.0, rz_2 = 0.0, alpha;
    VEC *r = NULL, *r_1 = NULL, *z = NULL, *Ap = NULL;
    VEC *p = NULL, *p_1 = NULL, *x;

    assert(b->nvec == 1);
    x = *x_ptr;
    nits = 0;
    if (IsZero(rb_norm = phgVecNorm2(b, 0, NULL))) {
	bzero(x->data, sizeof(*x->data) * x->map->nlocal);
	return 0;
    }

    r = phgMapCreateVec(x->map, 1);
    r_1 = phgMapCreateVec(x->map, 1);
    z = phgMapCreateVec(x->map, 1);
    Ap = phgMapCreateVec(x->map, 1);
    p = phgMapCreateVec(x->map, 1);
    p_1 = phgMapCreateVec(x->map, 1);

    /* r = b - A * _t_U  */
    if (!IsZero(phgVecNorm1(x, 0, NULL) / A->cmap->nglobal)) {
	phgMatVec(MAT_OP_N, 1.0, A, x, 0.0, &r_1);
	phgVecCopy(r_1, &r);
	phgVecAXPBY(1.0, b, -1.0, &r);
    }
    else {
	phgVecCopy(b, &r);
    }

    r0_norm = residual = phgVecNorm2(r, 0, NULL);

    tol = get_tol(r0_norm * rtol, rb_norm * btol, atol);

    if (verb > 0) {
	phgPrintf("======================== PCG =====================\n");
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
		  (double)atol, (double)rtol, (double)btol);
	phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
	phgPrintf("-----   ------------   -------------  ------------\n");
	phgPrintf("%5d   %12le   %12le   %12le\n", 0,
		  (double)residual,
		  (double)(residual / r0_norm), (double)(residual / rb_norm));
    }

    while (residual > tol && nits < max_its) {
	assert(pc != NULL);
	pc->pc_proc(pc->ctx, r, &z);
	nits++;
	if (nits == 1) {
	    phgVecCopy(z, &p);
	    rz_1 = phgVecDot(r, 0, z, 0, NULL);
	}
	else {
	    rz_2 = rz_1;
	    rz_1 = phgVecDot(r, 0, z, 0, NULL);
	    phgVecCopy(p, &p_1);
	    if (rz_2 == 0.0) {
		/*phgError(1, "%s:%d, the matrix may not be SPD.\n",
		   __FILE__, __LINE__); */
		if (warned == NULL || !*warned) {
		    phgPrintf("%s:%d, Warning! the matrix may not be SPD.\n",
			      __FILE__, __LINE__);
		    if (warned != NULL)
			*warned = TRUE;
		}
		break;
	    }
	    phgVecAXPBY(1.0, z, rz_1 / rz_2, &p);
	}
	phgMatVec(MAT_OP_N, 1.0, A, p, 0.0, &Ap);
	if ((alpha = phgVecDot(p, 0, Ap, 0, NULL)) == 0.0) {
	    /*phgError(1, "%s:%d, the matrix may not be SPD.\n",
	       __FILE__, __LINE__); */
	    if (warned == NULL || !*warned) {
		phgPrintf("%s:%d, Warning! the matrix may not be SPD.\n",
			  __FILE__, __LINE__);
		if (warned != NULL)
		    *warned = TRUE;
	    }
	    break;
	}
	alpha = rz_1 / alpha;
	phgVecAXPBY(alpha, p, 1.0, &x);
	phgVecCopy(r, &r_1);
	phgVecAXPBY(-1.0 * alpha, Ap, 1.0, &r);
	residual = phgVecNorm2(r, 0, NULL);
	if (verb > 0) {
	    phgPrintf("%5d   %12e   %12e   %12e", nits, (double)residual,
		      (double)(residual / r0_norm),
		      (double)(residual / rb_norm));
	    if (pc->solver_owned != NULL &&
#if 0
		pc->solver_owned != solver &&
#endif
		pc->solver_owned->nits > 0)
		phgPrintf(" (%d it%s in PC)\n", pc->solver_owned->nits,
			  pc->solver_owned->nits > 1 ? "s" : ".");
	    else
		phgPrintf("\n");
	}
    }

    phgVecDestroy(&r);
    phgVecDestroy(&r_1);
    phgVecDestroy(&z);
    phgVecDestroy(&Ap);
    phgVecDestroy(&p_1);
    phgVecDestroy(&p);

    if (res != NULL)
	*res = residual;

    return nits;
}

static int
solver_BiCGstab(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc, int max_its,
		    FLOAT rtol, FLOAT btol, FLOAT atol, FLOAT *res, INT verb)
{
    int nits;
    FLOAT residual = 0., r0_norm, rb_norm, tol;
    VEC *r_1 = NULL;
    VEC *x;

    VEC *rg;
    VEC *rh;
    VEC *pg;
    VEC *ph;
    VEC *sg;
    VEC *sh;
    VEC *tg;
    VEC *vg;

    FLOAT rho0 = 0, rho1 = 0;
    FLOAT alpha = 0, beta = 0, omega = 0;

    assert(b->nvec == 1);
    x = *x_ptr;
    nits = 0;
    if (IsZero(rb_norm = phgVecNorm2(b, 0, NULL))) {
	bzero(x->data, sizeof(*x->data) * x->map->nlocal);
	return 0;
    }

    r_1 = phgMapCreateVec(x->map, 1);
    rg = phgMapCreateVec(x->map, 1);
    rh = phgMapCreateVec(x->map, 1);
    pg = phgMapCreateVec(x->map, 1);
    ph = phgMapCreateVec(x->map, 1);
    sg = phgMapCreateVec(x->map, 1);
    sh = phgMapCreateVec(x->map, 1);
    tg = phgMapCreateVec(x->map, 1);
    vg = phgMapCreateVec(x->map, 1);

    /* r = b - A * _t_U  */
    if (!IsZero(phgVecNorm1(x, 0, NULL) / A->cmap->nglobal)) {
	phgMatVec(MAT_OP_N, 1.0, A, x, 0.0, &r_1);
	phgVecCopy(r_1, &rg);
	phgVecAXPBY(1.0, b, -1.0, &rg);
    }
    else {
	phgVecCopy(b, &rg);
    }

    r0_norm = residual = phgVecNorm2(rg, 0, NULL);

    tol = get_tol(r0_norm * rtol, rb_norm * btol, atol);

    if (verb > 0) {
	phgPrintf
	    ("======================== BiCGstab =====================\n");
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n", (double)atol,
		  (double)rtol, (double)btol);
	phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
	phgPrintf("-----   ------------   -------------  ------------\n");
	phgPrintf("%5d   %12le   %12le   %12le\n", 0,
		  (double)residual,
		  (double)(residual / r0_norm), (double)(residual / rb_norm));
    }

    /* set rh = r0, init sh, ph */
    phgVecCopy(rg, &rh);
    bzero(sh->data, sh->map->nlocal * sizeof(*sh->data));
    bzero(ph->data, ph->map->nlocal * sizeof(*ph->data));

    while (residual > tol && nits < max_its) {
	nits++;

	rho1 = phgVecDot(rg, 0, rh, 0, NULL);

	/* check if fails */
	if (rho1 == 0.0) {
	    phgPrintf("phg: bicgstab: method failed.!\n");
	    break;
	}

	if (nits == 1) {
	    /* p = r */
	    phgVecCopy(rg, &pg);
	}
	else {			/* update beta, p */
	    beta = (rho1 * alpha) / (rho0 * omega);

	    phgVecCopy(pg, &sg);
	    phgVecAXPBY(-omega, vg, 1., &sg);

	    phgVecCopy(rg, &pg);
	    phgVecAXPBY(beta, sg, 1., &pg);
	}

	/* save rho1 to rho0 */
	rho0 = rho1;

	/* precondition, solve ph */
	assert(pc != NULL);
	pc->pc_proc(pc->ctx, pg, &ph);

	/* v = A ph */
	phgMatVec(MAT_OP_N, 1.0, A, ph, 0.0, &vg);

	/* update alpha, s */
	alpha = rho1 / phgVecDot(rh, 0, vg, 0, NULL);
	phgVecCopy(rg, &sg);
	phgVecAXPBY(-alpha, vg, 1., &sg);

	/* check if s is small enough */
	omega = phgVecNorm2(sg, 0, NULL);
	if (omega <= 1e-60) {
	    phgPrintf("phg: bicgstab: ||s|| is too small: %e, terminated.\n",
		      omega);

	    phgVecAXPBY(alpha, ph, 1., &x);

	    phgMatVec(MAT_OP_N, 1.0, A, x, 0.0, &r_1);
	    phgVecCopy(b, &rg);
	    phgVecAXPBY(-1, r_1, 1., &rg);

	    residual = phgVecNorm2(rg, 0, NULL);
	    break;
	}

	/* precondition, solve sh, M sh = sg */
	assert(pc != NULL);
	pc->pc_proc(pc->ctx, sg, &sh);

	/* t = A sh */
	phgMatVec(MAT_OP_N, 1.0, A, sh, 0.0, &tg);

	/* update omega, x, r */
	omega = phgVecDot(tg, 0, sg, 0, NULL) / phgVecDot(tg, 0, tg, 0, NULL);

	phgVecAXPBY(omega, sh, alpha, &ph);
	phgVecAXPBY(1, ph, 1, &x);
	phgVecCopy(sg, &rg);
	phgVecAXPBY(-omega, tg, 1, &rg);

	/* check error */
	residual = phgVecNorm2(rg, 0, NULL);

	if (verb > 0) {
	    phgPrintf("%5d   %12e   %12e   %12e\n", nits,
		      residual, residual / r0_norm, residual / rb_norm);
	}

	if (residual <= tol)
	    break;
    }

#if 0
    phgMatVec(MAT_OP_N, 1.0, A, x, 0.0, &r_1);
    phgVecCopy(r_1, &rg);
    phgVecAXPBY(1.0, b, -1.0, &rg);
    phgPrintf("*** final residual %12e\n", phgVecNorm2(rg, 0, NULL));
#endif

    phgVecDestroy(&r_1);
    phgVecDestroy(&rg);
    phgVecDestroy(&rh);
    phgVecDestroy(&pg);
    phgVecDestroy(&ph);
    phgVecDestroy(&sg);
    phgVecDestroy(&sh);
    phgVecDestroy(&tg);
    phgVecDestroy(&vg);

    if (res != NULL)
	*res = residual;

    return nits;
}

/* GMRES, left preconditioned */
static int
solver_GMRES_l(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc, int restart,
		int augk, int max_its, FLOAT rtol, FLOAT btol, FLOAT atol,
		FLOAT *res, INT verb)
    /*
     *  for details to the algorithm GMRES(m)-Aya see
     *
     *  R. Barrett, M. Berry, T. F. Chan, J. Demmel, J. Donato, J. Dongarra,
     *  V. Eijkhout, R. Pozo, Ch. Romine, H. van der Vorst:
     *  Templates for the Solution of Linear Systems: Building Blocks
     *  for Iterative Solvers;
     *  SIAM, Philadelphia, 1994
     *
     */
{
    int nits, i, j, k;
    BOOLEAN break_flag = FALSE, false_conv = FALSE;
    INT nlocal;
    FLOAT residual = 0., r0_norm, r_norm, b_norm, tol, v_norm;
    FLOAT tmp, h1, h2;
    FLOAT **H = NULL, *h_data = NULL, **HH = NULL, *hh_data = NULL;
    FLOAT *y = NULL, *s = NULL, *c1 = NULL, *c2 = NULL;
    VEC *r = NULL, *z = NULL, *v1 = NULL;
    VEC **v_ptr = NULL, *x;

    Unused(augk);
    assert(b->nvec == 1);

    if (restart <= 0)
	restart = max_its + 1;

    x = *x_ptr;
    nlocal = x->map->nlocal;

    b_norm = phgVecNorm2(b, 0, NULL);
    if (IsZero(b_norm)) {
	bzero(x->data, sizeof(*x->data) * nlocal);
	return 0;
    }
    /*  allocation of matrix H and Vec v */
    v_ptr = (VEC **)phgAlloc((restart + 2) * sizeof(*v_ptr));
    H = (FLOAT **)phgAlloc((restart + 2) * sizeof(*H));
    HH = (FLOAT **)phgAlloc((restart + 2) * sizeof(*HH));
    h_data =
	(FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));
    hh_data =
	(FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));
    for (i = 0; i < restart + 2; i++) {
	v_ptr[i] = phgMapCreateVec(x->map, 1);
	H[i] = h_data + i * (restart + 2);
	HH[i] = hh_data + i * (restart + 2);
    }

    y = (FLOAT *)phgAlloc((restart + 2) * 4 * sizeof(*y));
    s = y + restart + 2;
    c1 = y + (restart + 2) * 2;
    c2 = y + (restart + 2) * 3;
    z = phgMapCreateVec(x->map, 1);

    /* loop for 'MaxIter' GMRES cycles */
    /* v_ptr[0]= r = b - A * x */
    r = v_ptr[0];
    phgVecCopy(b, &r);
    if (phgVecNorm1(x, 0, NULL) != 0.)
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);

    r0_norm = r_norm = residual = phgVecNorm2(r, 0, NULL);

    tol = get_tol(r0_norm * rtol, b_norm * btol, atol);

    if (verb > 0) {
	phgPrintf("======================= GMRES ====================\n");
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
		  (double)atol, (double)rtol, (double)btol);
	phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
	phgPrintf("-----   ------------   -------------  ------------\n");
	phgPrintf("%5d   %12le   %12le   %12le\n", 0,
		  (double)r_norm,
		  (double)(r_norm / r0_norm), (double)(r_norm / b_norm));
    }
    nits = 0;

    /* M*z = b - A*x */
    assert(pc != NULL);
    pc->pc_proc(pc->ctx, r, &z);

    s[0] = 1.;
    while (nits < max_its) {
	v_norm = s[1] = r_norm = phgVecNorm2(z, 0, NULL);
	phgMatVec(MAT_OP_N, 1.0 / r_norm, NULL, z, 0.0, &v_ptr[1]);
	/* GMRES iteration */
	i = 0;
	bzero(h_data, (restart + 2) * (restart + 2) * sizeof(*h_data));
	bzero(hh_data, (restart + 2) * (restart + 2) * sizeof(*h_data));
	/* v_ptr[1...i]= (M^-1 * A)^(i-1) * (M^-1)(b-A*x)  */
	while (i < restart && (r_norm > tol || false_conv) && nits < max_its) {
	    i++, nits++;
	    v1 = v_ptr[i + 1];
	    /* M*v[j+1] = A*v[j] */
	    phgMatVec(MAT_OP_N, 1.0, A, v_ptr[i], 0.0, &r);
	    assert(pc != NULL);
	    pc->pc_proc(pc->ctx, r, &v1);

	    /* modified Gram-Schmidt orthogonalization */
	    for (k = 1; k <= i; k++) {
		/* H[i][k]=<v[i+1],v[k]> */
		HH[k][i] = H[i][k] = phgVecDot(v1, 0, v_ptr[k], 0, NULL);
		/* v[i+1]=v[i+1]- H[i][k]*v[k] */
		phgMatVec(MAT_OP_N, -H[i][k], NULL, v_ptr[k], 1.0, &v1);
	    }			/* end for k-iter */

	    HH[i + 1][i] = H[i][i + 1] = phgVecNorm2(v1, 0, NULL);
	    if (IsZero(H[i][i + 1])) {
		phgInfo(2, "%le happy break i=%d ", (double)H[i][i + 1], i);
		i--;
		break_flag = TRUE;
		break;
	    }
	    phgMatVec(MAT_OP_N, 0.0, NULL, NULL, 1.0 / H[i][i + 1], &v1);

	    /* Q-R algoritm */
	    for (j = 1; j < i; j++) {
		h1 = c1[j] * H[i][j] + c2[j] * H[i][j + 1];
		h2 = -c2[j] * H[i][j] + c1[j] * H[i][j + 1];
		H[i][j] = h1;
		H[i][j + 1] = h2;
	    }

	    tmp = Sqrt(H[i][i] * H[i][i] + H[i][i + 1] * H[i][i + 1]);
	    if (IsZero(tmp)) {
		tmp = 1.e-16;
	    }
	    c1[i] = H[i][i] / tmp;
	    c2[i] = H[i][i + 1] / tmp;
	    s[i + 1] = -c2[i] * s[i];
	    s[i] = c1[i] * s[i];
	    H[i][i] = c1[i] * H[i][i] + c2[i] * H[i][i + 1];
	    r_norm = Fabs(s[i + 1]);
	    /* check whether real converged */
	    if (false_conv == TRUE && r_norm * s[0] <= tol) {
		false_conv = FALSE;
	    }

	    if (verb > 0) {
		phgPrintf("%5d   %12e   %12e   %12e", nits,
			  (double)r_norm, (double)(r_norm / r0_norm),
			  (double)(r_norm / b_norm));
		if (pc->solver_owned != NULL &&
#if 0
		    pc->solver_owned != solver &&
#endif
		    pc->solver_owned->nits > 0)
		    phgPrintf(" (%d it%s in PC)\n", pc->solver_owned->nits,
			      pc->solver_owned->nits > 1 ? "s" : ".");
		else
		    phgPrintf("\n");
	    }
	}			/* end-for-while(i)-sub-space */

	/* solve H * y = s */
	for (j = i; j > 0; j--) {
	    y[j] = s[j] / H[j][j];
	    for (k = j - 1; k > 0; k--) {
		s[k] -= H[j][k] * y[j];
	    }
	}
	/* Update x */
	/*  r =  V * y   */
#if 0
	phgMatVec(MAT_OP_N, y[i], NULL, v_ptr[i], 0.0, &z);
	for (j = i - 1; j >= 1; j--) {
	    phgMatVec(MAT_OP_N, y[j], NULL, v_ptr[j], 1.0, &z);
	}
#else
	if (i > 0) {
	    phgMatVec(MAT_OP_N, y[1], NULL, v_ptr[1], 0.0, &z);
	    for (j = 2; j <= i; j++) {
		phgMatVec(MAT_OP_N, y[j], NULL, v_ptr[j], 1.0, &z);
	    }
	}
#endif
	/*  x = x +  r  */
	phgMatVec(MAT_OP_N, 1.0, NULL, z, 1.0, &x);
	false_conv = FALSE;

	if (r_norm <= tol) {
	    /* r = b - A * x */
	    phgVecCopy(b, &r);
	    phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
	    residual = phgVecNorm2(r, 0, NULL);
	    if (residual <= tol) {
		if (verb > 0)
		    phgPrintf
			("GMRES iteration %d ||b-Ax||_2 = %e  rel.res %e\n",
			 nits, (double)residual,
			 (double)(residual / r0_norm));
		break;
	    }
	    else if (break_flag) {
		if (verb > 0) {
		    phgPrintf("You reached the happy break down.\n");
		    phgPrintf("GMRES iteration %d ||b-Ax||_2 = %e  "
			      "rel.res %e\n", nits,
			      (double)residual, (double)(residual / r0_norm));
		}
		break;
	    }
	    else if (i < restart) {
		if (verb > 0) {
		    phgPrintf("false convergence, restart iterations\n");
		}
		false_conv = TRUE;
		/* compute the factor */
		s[0] = residual / r_norm;
	    }

	    assert(pc != NULL);
	    pc->pc_proc(pc->ctx, r, &z);
	    r_norm = phgVecNorm2(z, 0, NULL);
	}
	else {
/* computer residual vector and continue loop 
 * z_{i+1}=M^-1 * r_{i+1} =M^-1 * (b  - A*x_{i+1})
 *         = M^-1 *(b  - A(x0 + V_i y))
 *         = z0 - M^-1 *A*V_i* y 
 *         = z0 - V_{i+1} * {HH_i}^T * y
*/
#if 0
	    for (j = 1; j <= i + 1; j++) {
		s[j] = 0.;
		for (k = 1; k <= i; k++)
		    s[j] -= HH[j][k] * y[k];
	    }

	    if (i > 0) {
		phgMatVec(MAT_OP_N, s[1] + v_norm, NULL, v_ptr[1], 0.0, &z);
		for (j = 2; j <= i + 1; j++)
		    phgMatVec(MAT_OP_N, s[j], NULL, v_ptr[j], 1.0, &z);
	    }
#else
	    Unused(v_norm);
	    phgVecCopy(b, &r);
	    phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
	    assert(pc != NULL);
	    pc->pc_proc(pc->ctx, r, &z);
#endif
	}
    }				/* main-loop */

    if (res != NULL)
	*res = residual;

    phgVecDestroy(&z);

    for (i = 0; i < restart + 2; i++)
	phgVecDestroy(&v_ptr[i]);
    phgFree(v_ptr);
    phgFree(H);
    phgFree(HH);
    phgFree(h_data);
    phgFree(hh_data);
    phgFree(y);

    return nits;
}

/* GMRES, right preconditioned */
static int
solver_GMRES_r(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc, int restart,
		int augk, int max_its, FLOAT rtol, FLOAT btol, FLOAT atol,
		FLOAT *res, INT verb)
{
    int nits, i, j, k;
    INT nlocal;
    FLOAT residual = 0., r0_norm, r_norm, b_norm, tol;
    FLOAT tmp, h1, h2;
    FLOAT **H = NULL, *h_data = NULL, **HH = NULL, *hh_data = NULL;
    FLOAT *y = NULL, *s = NULL, *c1 = NULL, *c2 = NULL;
    VEC *r = NULL, *z = NULL, *v1 = NULL;
    VEC **v_ptr = NULL, *x;

    Unused(augk);
    assert(b->nvec == 1);

    if (restart <= 0)
	restart = max_its + 1;

    x = *x_ptr;
    nlocal = x->map->nlocal;

    b_norm = phgVecNorm2(b, 0, NULL);
    if (IsZero(b_norm)) {
	bzero(x->data, sizeof(*x->data) * nlocal);
	return 0;
    }
    /*  allocation of matrix H and Vec v */
    v_ptr = (VEC **)phgAlloc((restart + 2) * sizeof(*v_ptr));
    H = (FLOAT **)phgAlloc((restart + 2) * sizeof(*H));
    HH = (FLOAT **)phgAlloc((restart + 2) * sizeof(*HH));
    h_data = (FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));
    hh_data = (FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));
    for (i = 0; i < restart + 2; i++) {
	v_ptr[i] = phgMapCreateVec(x->map, 1);
	H[i] = h_data + i * (restart + 2);
	HH[i] = hh_data + i * (restart + 2);
    }

    y = (FLOAT *)phgAlloc((restart + 2) * 4 * sizeof(*y));
    s = y + restart + 2;
    c1 = y + (restart + 2) * 2;
    c2 = y + (restart + 2) * 3;
    z = phgMapCreateVec(x->map, 1);

    /* loop for 'MaxIter' GMRES cycles */
    /* v_ptr[0]= r = b - A * x */
    r = v_ptr[0];
    phgVecCopy(b, &z);
    if (phgVecNorm1(x, 0, NULL) != 0.)
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &z);

    r0_norm = r_norm = residual = phgVecNorm2(z, 0, NULL);

    tol = get_tol(r0_norm * rtol, b_norm * btol, atol);

    if (verb > 0) {
	phgPrintf("======================= GMRES ====================\n");
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
		  (double)atol, (double)rtol, (double)btol);
	phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
	phgPrintf("-----   ------------   -------------  ------------\n");
	phgPrintf("%5d   %12le   %12le   %12le\n", 0,
		  (double)r_norm,
		  (double)(r_norm / r0_norm), (double)(r_norm / b_norm));
    }
    nits = 0;

    s[0] = 1.;
    while (nits < max_its) {
	s[1] = r_norm = phgVecNorm2(z, 0, NULL);
	phgMatVec(MAT_OP_N, 1.0 / r_norm, NULL, z, 0.0, &v_ptr[1]);

	/* GMRES iteration */
	i = 0;
	bzero(h_data, (restart + 2) * (restart + 2) * sizeof(*h_data));
	bzero(hh_data, (restart + 2) * (restart + 2) * sizeof(*h_data));

	/* Gram-Schmidt */
	while (i < restart && r_norm > tol && nits < max_its) {
	    i++, nits++;
	    v1 = v_ptr[i + 1];

	    /* preconditioning */
	    assert(pc != NULL);
	    pc->pc_proc(pc->ctx, v_ptr[i], &r);
	    phgMatVec(MAT_OP_N, 1.0, A, r, 0.0, &v1);

	    /* modified Gram-Schmidt orthogonalization */
	    for (k = 1; k <= i; k++) {
		/* H[i][k]=<v[i+1],v[k]> */
		HH[k][i] = H[i][k] = phgVecDot(v1, 0, v_ptr[k], 0, NULL);

		/* v[i+1]=v[i+1]- H[i][k]*v[k] */
		phgMatVec(MAT_OP_N, -H[i][k], NULL, v_ptr[k], 1.0, &v1);
	    }			/* end for k-iter */

	    HH[i + 1][i] = H[i][i + 1] = phgVecNorm2(v1, 0, NULL);
	    if (IsZero(H[i][i + 1])) {
		phgInfo(2, "%le happy break i=%d ", (double)H[i][i + 1], i);
		i--;
		break;
	    }
	    phgMatVec(MAT_OP_N, 0.0, NULL, NULL, 1.0 / H[i][i + 1], &v1);

	    /* Q-R algoritm */
	    for (j = 1; j < i; j++) {
		h1 = c1[j] * H[i][j] + c2[j] * H[i][j + 1];
		h2 = -c2[j] * H[i][j] + c1[j] * H[i][j + 1];
		H[i][j] = h1;
		H[i][j + 1] = h2;
	    }

	    tmp = Sqrt(H[i][i] * H[i][i] + H[i][i + 1] * H[i][i + 1]);
	    if (IsZero(tmp)) {
		tmp = 1.e-16;
	    }
	    c1[i] = H[i][i] / tmp;
	    c2[i] = H[i][i + 1] / tmp;
	    s[i + 1] = -c2[i] * s[i];
	    s[i] = c1[i] * s[i];
	    H[i][i] = c1[i] * H[i][i] + c2[i] * H[i][i + 1];
	    r_norm = Fabs(s[i + 1]);

	    if (verb > 0) {
		phgPrintf("%5d   %12e   %12e   %12e", nits,
			  (double)r_norm, (double)(r_norm / r0_norm),
			  (double)(r_norm / b_norm));
		if (pc->solver_owned != NULL &&
#if 0
		    pc->solver_owned != solver &&
#endif
		    pc->solver_owned->nits > 0)
		    phgPrintf(" (%d it%s in PC)\n", pc->solver_owned->nits,
			      pc->solver_owned->nits > 1 ? "s" : ".");
		else
		    phgPrintf("\n");
	    }

	    /* check whether converged */
	    if (r_norm <= tol) break; 
	}			/* end-for-while(i)-sub-space */

	/* solve H * y = s */
	for (j = i; j > 0; j--) {
	    y[j] = s[j] / H[j][j];
	    for (k = j - 1; k > 0; k--) {
		s[k] -= H[j][k] * y[j];
	    }
	}

	/* Update x */
	/*  z =  V * y   */
	if (i > 0) {
	    phgMatVec(MAT_OP_N, y[1], NULL, v_ptr[1], 0.0, &z);
	    for (j = 2; j <= i; j++) {
		phgMatVec(MAT_OP_N, y[j], NULL, v_ptr[j], 1.0, &z);
	    }
	}

        /* Mz = r */
        assert(pc != NULL);
        pc->pc_proc(pc->ctx, z, &r);

	/*  x = x +  r  */
	phgMatVec(MAT_OP_N, 1.0, NULL, r, 1.0, &x);

	if (r_norm <= tol) {

            /* verification, output should be 1. */
#if 0
            /* r = b - A * x */
            phgVecCopy(b, &r);
            phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
            residual = phgVecNorm2(r, 0, NULL);

            phgPrintf("r_norm / residual %e\n",  r_norm / residual);
#endif

            residual = r_norm;

            if (verb > 0)
                phgPrintf("GMRES iteration %d ||b-Ax||_2 = %e  rel.res %e\n", nits, (double)residual,
                        (double)(residual / r0_norm));

            break;
        }
        else {
            phgVecCopy(b, &z);
            phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &z);
	}
    }				/* main-loop */

    if (res != NULL) *res = residual;

    phgVecDestroy(&z);

    for (i = 0; i < restart + 2; i++)
	phgVecDestroy(&v_ptr[i]);

    phgFree(v_ptr);
    phgFree(H);
    phgFree(HH);
    phgFree(h_data);
    phgFree(hh_data);
    phgFree(y);

    return nits;
}

/* Flexible GMRES, right preconditioned */
static int
solver_FGMRES_r(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc, int restart,
		int augk, int max_its, FLOAT rtol, FLOAT btol, FLOAT atol,
		FLOAT *res, INT verb)
{
    int nits, i, j, k;
    INT nlocal;
    FLOAT residual = 0., r0_norm, r_norm, b_norm, tol;
    FLOAT tmp, h1, h2;
    FLOAT **H = NULL, *h_data = NULL, **HH = NULL, *hh_data = NULL;
    FLOAT *y = NULL, *s = NULL, *c1 = NULL, *c2 = NULL;
    VEC *r = NULL, *z = NULL, *v1 = NULL;
    VEC **v_ptr = NULL, *x;
    VEC **flv;

    Unused(augk);
    assert(b->nvec == 1);

    if (restart <= 0)
	restart = max_its + 1;

    x = *x_ptr;
    nlocal = x->map->nlocal;

    b_norm = phgVecNorm2(b, 0, NULL);
    if (IsZero(b_norm)) {
	bzero(x->data, sizeof(*x->data) * nlocal);
	return 0;
    }
    /*  allocation of matrix H and Vec v */
    v_ptr = (VEC **)phgAlloc((restart + 2) * sizeof(*v_ptr));
    flv = (VEC **)phgAlloc((restart + 2) * sizeof(*v_ptr));

    H = (FLOAT **)phgAlloc((restart + 2) * sizeof(*H));
    HH = (FLOAT **)phgAlloc((restart + 2) * sizeof(*HH));
    h_data = (FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));
    hh_data = (FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));
    for (i = 0; i < restart + 2; i++) {
	v_ptr[i] = phgMapCreateVec(x->map, 1);
	H[i] = h_data + i * (restart + 2);
	HH[i] = hh_data + i * (restart + 2);

        if (i > 0) flv[i] = phgMapCreateVec(x->map, 1);
    }

    y = (FLOAT *)phgAlloc((restart + 2) * 4 * sizeof(*y));
    s = y + restart + 2;
    c1 = y + (restart + 2) * 2;
    c2 = y + (restart + 2) * 3;
    z = phgMapCreateVec(x->map, 1);

    /* loop for 'MaxIter' GMRES cycles */
    /* v_ptr[0]= r = b - A * x */
    r = v_ptr[0];
    phgVecCopy(b, &z);
    if (phgVecNorm1(x, 0, NULL) != 0.)
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &z);

    r0_norm = r_norm = residual = phgVecNorm2(z, 0, NULL);

    tol = get_tol(r0_norm * rtol, b_norm * btol, atol);

    if (verb > 0) {
	phgPrintf("======================= GMRES ====================\n");
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
		  (double)atol, (double)rtol, (double)btol);
	phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
	phgPrintf("-----   ------------   -------------  ------------\n");
	phgPrintf("%5d   %12le   %12le   %12le\n", 0,
		  (double)r_norm,
		  (double)(r_norm / r0_norm), (double)(r_norm / b_norm));
    }
    nits = 0;

    s[0] = 1.;
    while (nits < max_its) {
	s[1] = r_norm = phgVecNorm2(z, 0, NULL);
	phgMatVec(MAT_OP_N, 1.0 / r_norm, NULL, z, 0.0, &v_ptr[1]);

	/* GMRES iteration */
	i = 0;
	bzero(h_data, (restart + 2) * (restart + 2) * sizeof(*h_data));
	bzero(hh_data, (restart + 2) * (restart + 2) * sizeof(*h_data));

	/* v_ptr[1...i]= (M^-1 * A)^(i-1) * (M^-1)(b-A*x)  */
	while (i < restart && r_norm > tol && nits < max_its) {
	    i++, nits++;
	    v1 = v_ptr[i + 1];

	    /* preconditioning */
	    assert(pc != NULL);
	    pc->pc_proc(pc->ctx, v_ptr[i], &r);
	    phgMatVec(MAT_OP_N, 1.0, A, r, 0.0, &v1);

            /* back up z */
            phgVecCopy(r, &flv[i]);

	    /* modified Gram-Schmidt orthogonalization */
	    for (k = 1; k <= i; k++) {
		/* H[i][k]=<v[i+1],v[k]> */
		HH[k][i] = H[i][k] = phgVecDot(v1, 0, v_ptr[k], 0, NULL);
		/* v[i+1]=v[i+1]- H[i][k]*v[k] */
		phgMatVec(MAT_OP_N, -H[i][k], NULL, v_ptr[k], 1.0, &v1);
	    }			/* end for k-iter */

	    HH[i + 1][i] = H[i][i + 1] = phgVecNorm2(v1, 0, NULL);
	    if (IsZero(H[i][i + 1])) {
		phgInfo(2, "%le happy break i=%d ", (double)H[i][i + 1], i);
		i--;
		break;
	    }
	    phgMatVec(MAT_OP_N, 0.0, NULL, NULL, 1.0 / H[i][i + 1], &v1);

	    /* Q-R algoritm */
	    for (j = 1; j < i; j++) {
		h1 = c1[j] * H[i][j] + c2[j] * H[i][j + 1];
		h2 = -c2[j] * H[i][j] + c1[j] * H[i][j + 1];
		H[i][j] = h1;
		H[i][j + 1] = h2;
	    }

	    tmp = Sqrt(H[i][i] * H[i][i] + H[i][i + 1] * H[i][i + 1]);
	    if (IsZero(tmp)) {
		tmp = 1.e-16;
	    }
	    c1[i] = H[i][i] / tmp;
	    c2[i] = H[i][i + 1] / tmp;
	    s[i + 1] = -c2[i] * s[i];
	    s[i] = c1[i] * s[i];
	    H[i][i] = c1[i] * H[i][i] + c2[i] * H[i][i + 1];
	    r_norm = Fabs(s[i + 1]);

	    if (verb > 0) {
		phgPrintf("%5d   %12e   %12e   %12e", nits,
			  (double)r_norm, (double)(r_norm / r0_norm),
			  (double)(r_norm / b_norm));
		if (pc->solver_owned != NULL &&
#if 0
		    pc->solver_owned != solver &&
#endif
		    pc->solver_owned->nits > 0)
		    phgPrintf(" (%d it%s in PC)\n", pc->solver_owned->nits,
			      pc->solver_owned->nits > 1 ? "s" : ".");
		else
		    phgPrintf("\n");
	    }

	    /* check whether converged */
	    if (r_norm <= tol) break; 
	}			/* end-for-while(i)-sub-space */

	/* solve H * y = s */
	for (j = i; j > 0; j--) {
	    y[j] = s[j] / H[j][j];
	    for (k = j - 1; k > 0; k--) {
		s[k] -= H[j][k] * y[j];
	    }
	}

	/* Update x */
	/*  z =  V * y   */
	if (i > 0) {
	    phgMatVec(MAT_OP_N, y[1], NULL, flv[1], 0.0, &z);
	    for (j = 2; j <= i; j++) {
		phgMatVec(MAT_OP_N, y[j], NULL, flv[j], 1.0, &z);
	    }
	}

	/*  x = x +  z  */
	phgMatVec(MAT_OP_N, 1.0, NULL, z, 1.0, &x);

	if (r_norm <= tol) {

            /* verification, output should be 1. */
#if 0
            /* r = b - A * x */
            phgVecCopy(b, &r);
            phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
            residual = phgVecNorm2(r, 0, NULL);

            phgPrintf("r_norm / residual %e\n",  r_norm / residual);
#endif

            residual = r_norm;

            if (verb > 0)
                phgPrintf("GMRES iteration %d ||b-Ax||_2 = %e  rel.res %e\n", nits, (double)residual,
                        (double)(residual / r0_norm));

            break;
        }
        else {
            phgVecCopy(b, &z);
            phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &z);
	}
    }				/* main-loop */

    if (res != NULL) *res = residual;

    phgVecDestroy(&z);

    for (i = 0; i < restart + 2; i++) {
	phgVecDestroy(&v_ptr[i]);
	if (i > 0) phgVecDestroy(&flv[i]);
    }

    phgFree(v_ptr);
    phgFree(flv);
    phgFree(H);
    phgFree(HH);
    phgFree(h_data);
    phgFree(hh_data);
    phgFree(y);

    return nits;
}

/* LGMRES, left preconditioned */
static int
solver_LGMRES_l(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc, int restart,
		int augk, int max_its, FLOAT rtol, FLOAT btol, FLOAT atol,
		FLOAT *res, INT verb)
{
    int nits, i, j, k;
    BOOLEAN break_flag = FALSE, false_conv = FALSE;
    INT nlocal;
    FLOAT residual = 0., r0_norm, r_norm, b_norm, tol, v_norm;
    FLOAT tmp, h1, h2;
    FLOAT **H = NULL, *h_data = NULL, **HH = NULL, *hh_data = NULL;
    FLOAT *y = NULL, *s = NULL, *c1 = NULL, *c2 = NULL;
    VEC *r = NULL, *z = NULL, *v1 = NULL;
    VEC **v_ptr = NULL, *x;
    VEC **z_ptr = NULL;
    int zi;
    INT itr_out;
    INT mk;

    assert(b->nvec == 1);

    if (restart <= 0)
	restart = max_its + 1;
    if (augk <= 0)
	augk = 2;

    x = *x_ptr;
    nlocal = x->map->nlocal;

    b_norm = phgVecNorm2(b, 0, NULL);
    if (IsZero(b_norm)) {
	bzero(x->data, sizeof(*x->data) * nlocal);
	return 0;
    }

    /* rewrite restart */
    mk = restart;
    restart += augk;

    /*  allocation of matrix H and Vec v */
    v_ptr = (VEC **)phgAlloc((restart + 2) * sizeof(*v_ptr));
    z_ptr = (VEC **)phgAlloc(augk * sizeof(*z_ptr));
    H = (FLOAT **)phgAlloc((restart + 2) * sizeof(*H));
    HH = (FLOAT **)phgAlloc((restart + 2) * sizeof(*HH));
    h_data =
	(FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));
    hh_data =
	(FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));
    for (i = 0; i < restart + 2; i++) {
	v_ptr[i] = phgMapCreateVec(x->map, 1);
	H[i] = h_data + i * (restart + 2);
	HH[i] = hh_data + i * (restart + 2);
    }

    /* z */
    for (i = 0; i < augk; i++) {
	z_ptr[i] = phgMapCreateVec(x->map, 1);
    }

    y = (FLOAT *)phgAlloc((restart + 2) * 4 * sizeof(*y));
    s = y + restart + 2;
    c1 = y + (restart + 2) * 2;
    c2 = y + (restart + 2) * 3;
    z = phgMapCreateVec(x->map, 1);

    /* loop for 'MaxIter' GMRES cycles */
    /* v_ptr[0]= r = b - A * x */
    r = v_ptr[0];
    phgVecCopy(b, &r);
    if (phgVecNorm1(x, 0, NULL) != 0.)
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);

    r0_norm = r_norm = residual = phgVecNorm2(r, 0, NULL);

    tol = get_tol(r0_norm * rtol, b_norm * btol, atol);

    if (verb > 0) {
	phgPrintf("======================= LGMRES ====================\n");
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
		  (double)atol, (double)rtol, (double)btol);
	phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
	phgPrintf("-----   ------------   -------------  ------------\n");
	phgPrintf("%5d   %12le   %12le   %12le\n", 0,
		  (double)r_norm,
		  (double)(r_norm / r0_norm), (double)(r_norm / b_norm));
    }
    nits = 0;

    /* M*z = b - A*x */
    assert(pc != NULL);
    pc->pc_proc(pc->ctx, r, &z);

    s[0] = 1.;
    itr_out = 0;
    while (nits < max_its) {
	/* tune restart */
	if (itr_out <= augk)
	    restart = mk + itr_out;
	else
	    restart = mk + augk;

	v_norm = s[1] = r_norm = phgVecNorm2(z, 0, NULL);
	phgMatVec(MAT_OP_N, 1.0 / r_norm, NULL, z, 0.0, &v_ptr[1]);

	/* GMRES iteration */
	i = 0;
	bzero(h_data, (restart + 2) * (restart + 2) * sizeof(*h_data));
	bzero(hh_data, (restart + 2) * (restart + 2) * sizeof(*h_data));
	/* v_ptr[1...i]= (M^-1 * A)^(i-1) * (M^-1)(b-A*x)  */
	while (i < restart && (r_norm > tol || false_conv) && nits < max_its) {
	    i++, nits++;
	    v1 = v_ptr[i + 1];

	    if (i <= mk) {
		/* M*v[j+1] = A*v[j] */
		phgMatVec(MAT_OP_N, 1.0, A, v_ptr[i], 0.0, &r);
	    }
	    else {
		/* M*z[j+1] = A*v[j - mk - 1] */
		phgMatVec(MAT_OP_N, 1.0, A, z_ptr[i - mk - 1], 0.0, &r);
	    }

	    assert(pc != NULL);
	    pc->pc_proc(pc->ctx, r, &v1);

	    /* modified Gram-Schmidt orthogonalization */
	    for (k = 1; k <= i; k++) {
		/* H[i][k]=<v[i+1],v[k]> */
		HH[k][i] = H[i][k] = phgVecDot(v1, 0, v_ptr[k], 0, NULL);
		/* v[i+1]=v[i+1]- H[i][k]*v[k] */
		phgMatVec(MAT_OP_N, -H[i][k], NULL, v_ptr[k], 1.0, &v1);
	    }			/* end for k-iter */

	    HH[i + 1][i] = H[i][i + 1] = phgVecNorm2(v1, 0, NULL);
	    if (IsZero(H[i][i + 1])) {
		phgInfo(2, "%le happy break i=%d ", (double)H[i][i + 1], i);
		i--;
		break_flag = TRUE;
		break;
	    }
	    phgMatVec(MAT_OP_N, 0.0, NULL, NULL, 1.0 / H[i][i + 1], &v1);

	    /* Q-R algoritm */
	    for (j = 1; j < i; j++) {
		h1 = c1[j] * H[i][j] + c2[j] * H[i][j + 1];
		h2 = -c2[j] * H[i][j] + c1[j] * H[i][j + 1];
		H[i][j] = h1;
		H[i][j + 1] = h2;
	    }

	    tmp = Sqrt(H[i][i] * H[i][i] + H[i][i + 1] * H[i][i + 1]);
	    if (IsZero(tmp)) {
		tmp = 1.e-16;
	    }
	    c1[i] = H[i][i] / tmp;
	    c2[i] = H[i][i + 1] / tmp;
	    s[i + 1] = -c2[i] * s[i];
	    s[i] = c1[i] * s[i];
	    H[i][i] = c1[i] * H[i][i] + c2[i] * H[i][i + 1];
	    r_norm = Fabs(s[i + 1]);
	    /* check whether real converged */
	    if (false_conv == TRUE && r_norm * s[0] <= tol) {
		false_conv = FALSE;
	    }

	    if (verb > 0) {
		phgPrintf("%5d   %12e   %12e   %12e", nits,
			  (double)r_norm, (double)(r_norm / r0_norm),
			  (double)(r_norm / b_norm));
		if (pc->solver_owned != NULL &&
#if 0
		    pc->solver_owned != solver &&
#endif
		    pc->solver_owned->nits > 0)
		    phgPrintf(" (%d it%s in PC)\n", pc->solver_owned->nits,
			      pc->solver_owned->nits > 1 ? "s" : ".");
		else
		    phgPrintf("\n");
	    }
	}			/* end-for-while(i)-sub-space */

	/* solve H * y = s */
	for (j = i; j > 0; j--) {
	    y[j] = s[j] / H[j][j];
	    for (k = j - 1; k > 0; k--) {
		s[k] -= H[j][k] * y[j];
	    }
	}

	/* Update x */
	/*  r =  V * y   */
	if (i > 0) {
	    if (i <= mk) {
		phgMatVec(MAT_OP_N, y[1], NULL, v_ptr[1], 0.0, &z);
		for (j = 2; j <= i; j++) {
		    phgMatVec(MAT_OP_N, y[j], NULL, v_ptr[j], 1.0, &z);
		}
	    }
	    else {
		int ak;

		phgMatVec(MAT_OP_N, y[1], NULL, v_ptr[1], 0.0, &z);
		for (j = 2; j <= mk; j++) {
		    phgMatVec(MAT_OP_N, y[j], NULL, v_ptr[j], 1.0, &z);
		}

		if (itr_out < augk)
		    ak = itr_out;
		else
		    ak = augk;

		for (j = 1; j <= ak; j++) {
		    phgMatVec(MAT_OP_N, y[j + mk], NULL, z_ptr[j - 1], 1.0,
			      &z);
		}
	    }
	}

	/*  x = x +  r  */
	phgMatVec(MAT_OP_N, 1.0, NULL, z, 1.0, &x);

	/* store z, z[i] = r  */
	zi = itr_out % augk;
	phgVecCopy(z, &z_ptr[zi]);

	false_conv = FALSE;

	if (r_norm <= tol) {
	    /* r = b - A * x */
	    phgVecCopy(b, &r);
	    phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
	    residual = phgVecNorm2(r, 0, NULL);
	    if (residual <= tol) {
		if (verb > 0)
		    phgPrintf
			("GMRES iteration %d ||b-Ax||_2 = %e  rel.res %e\n",
			 nits, (double)residual,
			 (double)(residual / r0_norm));
		break;
	    }
	    else if (break_flag) {
		if (verb > 0) {
		    phgPrintf("You reached the happy break down.\n");
		    phgPrintf("GMRES iteration %d ||b-Ax||_2 = %e  "
			      "rel.res %e\n", nits, (double)residual,
			      (double)(residual / r0_norm));
		}
		break;
	    }
	    else if (i < restart) {
		if (verb > 0) {
		    phgPrintf("false convergence, restart iterations\n");
		}
		false_conv = TRUE;
		/* compute the factor */
		s[0] = residual / r_norm;
	    }

	    assert(pc != NULL);
	    pc->pc_proc(pc->ctx, r, &z);
	    r_norm = phgVecNorm2(z, 0, NULL);
	}
	else {
	    /* computer residual vector and continue loop 
	     * z_{i+1}=M^-1 * r_{i+1} =M^-1 * (b  - A*x_{i+1})
	     *         = M^-1 *(b  - A(x0 + V_i y))
	     *         = z0 - M^-1 *A*V_i* y 
	     *         = z0 - V_{i+1} * {HH_i}^T * y
	     */
#if 0
	    for (j = 1; j <= i + 1; j++) {
		s[j] = 0.;
		for (k = 1; k <= i; k++)
		    s[j] -= HH[j][k] * y[k];
	    }

	    if (i > 0) {
		phgMatVec(MAT_OP_N, s[1] + v_norm, NULL, v_ptr[1], 0.0, &z);
		for (j = 2; j <= i + 1; j++)
		    phgMatVec(MAT_OP_N, s[j], NULL, v_ptr[j], 1.0, &z);
	    }
#else
	    Unused(v_norm);
	    phgVecCopy(b, &r);
	    phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
	    assert(pc != NULL);
	    pc->pc_proc(pc->ctx, r, &z);
#endif
	}

	/*  out itr */
	itr_out++;
    }

    if (res != NULL)
	*res = residual;

    phgVecDestroy(&z);

    /* z */
    for (i = 0; i < augk; i++)
	phgVecDestroy(&z_ptr[i]);
    phgFree(z_ptr);

    /* restart */
    restart = mk + augk;
    for (i = 0; i < restart + 2; i++)
	phgVecDestroy(&v_ptr[i]);
    phgFree(v_ptr);

    phgFree(H);
    phgFree(HH);
    phgFree(h_data);
    phgFree(hh_data);
    phgFree(y);

    return nits;
}

/* LGMRES, right preconditioned */
static int
solver_LGMRES_r(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc, int restart,
		int augk, int max_its, FLOAT rtol, FLOAT btol, FLOAT atol,
		FLOAT *res, INT verb)
{
    int nits, i, j, k;
    INT nlocal;
    FLOAT residual = 0., r0_norm, r_norm, b_norm, tol;
    FLOAT tmp, h1, h2;
    FLOAT **H = NULL, *h_data = NULL, **HH = NULL, *hh_data = NULL;
    FLOAT *y = NULL, *s = NULL, *c1 = NULL, *c2 = NULL;
    VEC *r = NULL, *z = NULL, *v1 = NULL;
    VEC **v_ptr = NULL, *x;
    VEC **z_ptr = NULL;
    int zi;
    INT itr_out;
    INT mk;

    assert(b->nvec == 1);

    if (restart <= 0)
	restart = max_its + 1;
    if (augk <= 0)
	augk = 2;

    x = *x_ptr;
    nlocal = x->map->nlocal;

    b_norm = phgVecNorm2(b, 0, NULL);
    if (IsZero(b_norm)) {
	bzero(x->data, sizeof(*x->data) * nlocal);
	return 0;
    }

    /* rewrite restart */
    mk = restart;
    restart += augk;

    /*  allocation of matrix H and Vec v */
    v_ptr = (VEC **)phgAlloc((restart + 2) * sizeof(*v_ptr));
    z_ptr = (VEC **)phgAlloc(augk * sizeof(*z_ptr));
    H = (FLOAT **)phgAlloc((restart + 2) * sizeof(*H));
    HH = (FLOAT **)phgAlloc((restart + 2) * sizeof(*HH));
    h_data = (FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));
    hh_data = (FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));

    for (i = 0; i < restart + 2; i++) {
	v_ptr[i] = phgMapCreateVec(x->map, 1);
	H[i] = h_data + i * (restart + 2);
	HH[i] = hh_data + i * (restart + 2);
    }

    /* z */
    for (i = 0; i < augk; i++) {
	z_ptr[i] = phgMapCreateVec(x->map, 1);
    }

    y = (FLOAT *)phgAlloc((restart + 2) * 4 * sizeof(*y));
    s = y + restart + 2;
    c1 = y + (restart + 2) * 2;
    c2 = y + (restart + 2) * 3;
    z = phgMapCreateVec(x->map, 1);

    /* loop for 'MaxIter' GMRES cycles */
    /* v_ptr[0]= r = b - A * x */
    r = v_ptr[0];
    phgVecCopy(b, &z);
    if (phgVecNorm1(x, 0, NULL) != 0.)
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &z);

    r0_norm = r_norm = residual = phgVecNorm2(z, 0, NULL);

    tol = get_tol(r0_norm * rtol, b_norm * btol, atol);

    if (verb > 0) {
	phgPrintf("======================= LGMRES ====================\n");
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
		  (double)atol, (double)rtol, (double)btol);
	phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
	phgPrintf("-----   ------------   -------------  ------------\n");
	phgPrintf("%5d   %12le   %12le   %12le\n", 0,
		  (double)r_norm,
		  (double)(r_norm / r0_norm), (double)(r_norm / b_norm));
    }
    nits = 0;

    s[0] = 1.;
    itr_out = 0;
    while (nits < max_its) {
	/* tune restart */
	if (itr_out <= augk)
	    restart = mk + itr_out;
	else
	    restart = mk + augk;

	s[1] = r_norm = phgVecNorm2(z, 0, NULL);
	phgMatVec(MAT_OP_N, 1.0 / r_norm, NULL, z, 0.0, &v_ptr[1]);

	/* GMRES iteration */
	i = 0;
	bzero(h_data, (restart + 2) * (restart + 2) * sizeof(*h_data));
	bzero(hh_data, (restart + 2) * (restart + 2) * sizeof(*h_data));

	/* v_ptr[1...i]= (M^-1 * A)^(i-1) * (M^-1)(b-A*x)  */
	while (i < restart && r_norm > tol && nits < max_its) {
	    i++, nits++;
	    v1 = v_ptr[i + 1];

            /* preconditioning */
	    assert(pc != NULL);

	    if (i <= mk) {
                pc->pc_proc(pc->ctx, v_ptr[i], &r);
	    }
	    else {
                pc->pc_proc(pc->ctx, z_ptr[i - mk - 1], &r);
	    }

            phgMatVec(MAT_OP_N, 1.0, A, r, 0.0, &v1);

	    /* modified Gram-Schmidt orthogonalization */
	    for (k = 1; k <= i; k++) {
		/* H[i][k]=<v[i+1],v[k]> */
		HH[k][i] = H[i][k] = phgVecDot(v1, 0, v_ptr[k], 0, NULL);

		/* v[i+1]=v[i+1]- H[i][k]*v[k] */
		phgMatVec(MAT_OP_N, -H[i][k], NULL, v_ptr[k], 1.0, &v1);
	    }			/* end for k-iter */

	    HH[i + 1][i] = H[i][i + 1] = phgVecNorm2(v1, 0, NULL);
	    if (IsZero(H[i][i + 1])) {
		phgInfo(2, "%le happy break i=%d ", (double)H[i][i + 1], i);
		i--;
		break;
	    }
	    phgMatVec(MAT_OP_N, 0.0, NULL, NULL, 1.0 / H[i][i + 1], &v1);

	    /* Q-R algoritm */
	    for (j = 1; j < i; j++) {
		h1 = c1[j] * H[i][j] + c2[j] * H[i][j + 1];
		h2 = -c2[j] * H[i][j] + c1[j] * H[i][j + 1];
		H[i][j] = h1;
		H[i][j + 1] = h2;
	    }

	    tmp = Sqrt(H[i][i] * H[i][i] + H[i][i + 1] * H[i][i + 1]);
	    if (IsZero(tmp)) {
		tmp = 1.e-16;
	    }
	    c1[i] = H[i][i] / tmp;
	    c2[i] = H[i][i + 1] / tmp;
	    s[i + 1] = -c2[i] * s[i];
	    s[i] = c1[i] * s[i];
	    H[i][i] = c1[i] * H[i][i] + c2[i] * H[i][i + 1];
	    r_norm = Fabs(s[i + 1]);

	    if (verb > 0) {
		phgPrintf("%5d   %12e   %12e   %12e", nits, (double)r_norm, (double)(r_norm / r0_norm),
			  (double)(r_norm / b_norm));
		if (pc->solver_owned != NULL &&
#if 0
		    pc->solver_owned != solver &&
#endif
		    pc->solver_owned->nits > 0)
		    phgPrintf(" (%d it%s in PC)\n", pc->solver_owned->nits,
			      pc->solver_owned->nits > 1 ? "s" : ".");
		else
		    phgPrintf("\n");
	    }

	    /* check whether converged */
	    if (r_norm <= tol) break;
	}			/* end-for-while(i)-sub-space */

	/* solve H * y = s */
	for (j = i; j > 0; j--) {
	    y[j] = s[j] / H[j][j];
	    for (k = j - 1; k > 0; k--) {
		s[k] -= H[j][k] * y[j];
	    }
	}

	/* Update x */
	/*  r =  V * y   */
	if (i > 0) {
	    if (i <= mk) {
		phgMatVec(MAT_OP_N, y[1], NULL, v_ptr[1], 0.0, &z);
		for (j = 2; j <= i; j++) {
		    phgMatVec(MAT_OP_N, y[j], NULL, v_ptr[j], 1.0, &z);
		}
	    }
	    else {
		int ak;

		phgMatVec(MAT_OP_N, y[1], NULL, v_ptr[1], 0.0, &z);
		for (j = 2; j <= mk; j++) {
		    phgMatVec(MAT_OP_N, y[j], NULL, v_ptr[j], 1.0, &z);
		}

		if (itr_out < augk)
		    ak = itr_out;
		else
		    ak = augk;

		for (j = 1; j <= ak; j++) {
		    phgMatVec(MAT_OP_N, y[j + mk], NULL, z_ptr[j - 1], 1.0, &z);
		}
	    }
	}

        /* pc */
        pc->pc_proc(pc->ctx, z, &r);

	/*  x = x +  r  */
	phgMatVec(MAT_OP_N, 1.0, NULL, r, 1.0, &x);

	/* store z, z[i] = z  */
	zi = itr_out % augk;
	phgVecCopy(z, &z_ptr[zi]);

	if (r_norm <= tol) {
	    residual = r_norm;

            if (verb > 0)
                phgPrintf("GMRES iteration %d ||b-Ax||_2 = %e  rel.res %e\n", nits, (double)residual,
                     (double)(residual / r0_norm));
            break;
        }
        else {
            phgVecCopy(b, &z);
            phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &z);
        }

        /*  out itr */
        itr_out++;
    }

    if (res != NULL) *res = residual;

    phgVecDestroy(&z);

    /* z */
    for (i = 0; i < augk; i++)
	phgVecDestroy(&z_ptr[i]);
    phgFree(z_ptr);

    /* restart */
    restart = mk + augk;
    for (i = 0; i < restart + 2; i++)
	phgVecDestroy(&v_ptr[i]);
    phgFree(v_ptr);

    phgFree(H);
    phgFree(HH);
    phgFree(h_data);
    phgFree(hh_data);
    phgFree(y);

    return nits;
}

/* Flexible LGMRES, right preconditioned */
static int
solver_FLGMRES_r(MAT *A, VEC *b, VEC **x_ptr, SOLVER_PC *pc, int restart,
		 int augk, int max_its, FLOAT rtol, FLOAT btol, FLOAT atol,
		 FLOAT *res, INT verb)
{
    int nits, i, j, k;
    INT nlocal;
    FLOAT residual = 0., r0_norm, r_norm, b_norm, tol;
    FLOAT tmp, h1, h2;
    FLOAT **H = NULL, *h_data = NULL, **HH = NULL, *hh_data = NULL;
    FLOAT *y = NULL, *s = NULL, *c1 = NULL, *c2 = NULL;
    VEC *r = NULL, *z = NULL, *v1 = NULL;
    VEC **v_ptr = NULL, *x;
    VEC **z_ptr = NULL;
    VEC **flv;
    int zi;
    INT itr_out;
    INT mk;

    assert(b->nvec == 1);

    if (restart <= 0)
	restart = max_its + 1;
    if (augk <= 0)
	augk = 2;

    x = *x_ptr;
    nlocal = x->map->nlocal;

    b_norm = phgVecNorm2(b, 0, NULL);
    if (IsZero(b_norm)) {
	bzero(x->data, sizeof(*x->data) * nlocal);
	return 0;
    }

    /* rewrite restart */
    mk = restart;
    restart += augk;

    /*  allocation of matrix H and Vec v */
    v_ptr = (VEC **)phgAlloc((restart + 2) * sizeof(*v_ptr));
    z_ptr = (VEC **)phgAlloc(augk * sizeof(*z_ptr));
    flv = (VEC **)phgAlloc((restart + 2) * sizeof(*v_ptr));

    H = (FLOAT **)phgAlloc((restart + 2) * sizeof(*H));
    HH = (FLOAT **)phgAlloc((restart + 2) * sizeof(*HH));
    h_data = (FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));
    hh_data = (FLOAT *)phgAlloc((restart + 2) * (restart + 2) * sizeof(*h_data));

    for (i = 0; i < restart + 2; i++) {
	v_ptr[i] = phgMapCreateVec(x->map, 1);
	H[i] = h_data + i * (restart + 2);
	HH[i] = hh_data + i * (restart + 2);
    }

    /* z */
    for (i = 0; i < augk; i++) {
	z_ptr[i] = phgMapCreateVec(x->map, 1);
    }

    for (i = 1; i < restart + 2; i++) {
	flv[i] = phgMapCreateVec(x->map, 1);
    }

    y = (FLOAT *)phgAlloc((restart + 2) * 4 * sizeof(*y));
    s = y + restart + 2;
    c1 = y + (restart + 2) * 2;
    c2 = y + (restart + 2) * 3;
    z = phgMapCreateVec(x->map, 1);

    /* loop for 'MaxIter' GMRES cycles */
    /* v_ptr[0]= r = b - A * x */
    r = v_ptr[0];
    phgVecCopy(b, &z);
    if (phgVecNorm1(x, 0, NULL) != 0.)
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &z);

    r0_norm = r_norm = residual = phgVecNorm2(z, 0, NULL);

    tol = get_tol(r0_norm * rtol, b_norm * btol, atol);

    if (verb > 0) {
	phgPrintf("======================= FLGMRES ===================\n");
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
		  (double)atol, (double)rtol, (double)btol);
	phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
	phgPrintf("-----   ------------   -------------  ------------\n");
	phgPrintf("%5d   %12le   %12le   %12le\n", 0,
		  (double)r_norm,
		  (double)(r_norm / r0_norm), (double)(r_norm / b_norm));
    }
    nits = 0;

    s[0] = 1.;
    itr_out = 0;
    while (nits < max_its) {
	/* tune restart */
	if (itr_out <= augk)
	    restart = mk + itr_out;
	else
	    restart = mk + augk;

	s[1] = r_norm = phgVecNorm2(z, 0, NULL);
	phgMatVec(MAT_OP_N, 1.0 / r_norm, NULL, z, 0.0, &v_ptr[1]);

	/* GMRES iteration */
	i = 0;
	bzero(h_data, (restart + 2) * (restart + 2) * sizeof(*h_data));
	bzero(hh_data, (restart + 2) * (restart + 2) * sizeof(*h_data));

	/* v_ptr[1...i]= (M^-1 * A)^(i-1) * (M^-1)(b-A*x)  */
	while (i < restart && r_norm > tol && nits < max_its) {
	    i++, nits++;
	    v1 = v_ptr[i + 1];

            /* preconditioning */
	    assert(pc != NULL);

	    if (i <= mk) {
                pc->pc_proc(pc->ctx, v_ptr[i], &r);
	    }
	    else {
                pc->pc_proc(pc->ctx, z_ptr[i - mk - 1], &r);
	    }

            phgMatVec(MAT_OP_N, 1.0, A, r, 0.0, &v1);

            /* backup */
            phgVecCopy(r, &flv[i]);

	    /* modified Gram-Schmidt orthogonalization */
	    for (k = 1; k <= i; k++) {
		/* H[i][k]=<v[i+1],v[k]> */
		HH[k][i] = H[i][k] = phgVecDot(v1, 0, v_ptr[k], 0, NULL);

		/* v[i+1]=v[i+1]- H[i][k]*v[k] */
		phgMatVec(MAT_OP_N, -H[i][k], NULL, v_ptr[k], 1.0, &v1);
	    }			/* end for k-iter */

	    HH[i + 1][i] = H[i][i + 1] = phgVecNorm2(v1, 0, NULL);
	    if (IsZero(H[i][i + 1])) {
		phgInfo(2, "%le happy break i=%d ", (double)H[i][i + 1], i);
		i--;
		break;
	    }
	    phgMatVec(MAT_OP_N, 0.0, NULL, NULL, 1.0 / H[i][i + 1], &v1);

	    /* Q-R algoritm */
	    for (j = 1; j < i; j++) {
		h1 = c1[j] * H[i][j] + c2[j] * H[i][j + 1];
		h2 = -c2[j] * H[i][j] + c1[j] * H[i][j + 1];
		H[i][j] = h1;
		H[i][j + 1] = h2;
	    }

	    tmp = Sqrt(H[i][i] * H[i][i] + H[i][i + 1] * H[i][i + 1]);
	    if (IsZero(tmp)) {
		tmp = 1.e-16;
	    }
	    c1[i] = H[i][i] / tmp;
	    c2[i] = H[i][i + 1] / tmp;
	    s[i + 1] = -c2[i] * s[i];
	    s[i] = c1[i] * s[i];
	    H[i][i] = c1[i] * H[i][i] + c2[i] * H[i][i + 1];
	    r_norm = Fabs(s[i + 1]);

	    if (verb > 0) {
		phgPrintf("%5d   %12e   %12e   %12e", nits, (double)r_norm, (double)(r_norm / r0_norm),
			  (double)(r_norm / b_norm));
		if (pc->solver_owned != NULL &&
#if 0
		    pc->solver_owned != solver &&
#endif
		    pc->solver_owned->nits > 0)
		    phgPrintf(" (%d it%s in PC)\n", pc->solver_owned->nits,
			      pc->solver_owned->nits > 1 ? "s" : ".");
		else
		    phgPrintf("\n");
	    }

	    /* check whether converged */
	    if (r_norm <= tol) break;
	}			/* end-for-while(i)-sub-space */

	/* solve H * y = s */
	for (j = i; j > 0; j--) {
	    y[j] = s[j] / H[j][j];
	    for (k = j - 1; k > 0; k--) {
		s[k] -= H[j][k] * y[j];
	    }
	}

	/* Update x */
	/*  r =  V * y   */
	if (i > 0) {
            phgMatVec(MAT_OP_N, y[1], NULL, flv[1], 0.0, &z);
            for (j = 2; j <= i; j++) {
                phgMatVec(MAT_OP_N, y[j], NULL, flv[j], 1.0, &z);
            }
	}

	/*  x = x +  z  */
	phgMatVec(MAT_OP_N, 1.0, NULL, z, 1.0, &x);

	if (i > 0) {
	    if (i <= mk) {
		phgMatVec(MAT_OP_N, y[1], NULL, v_ptr[1], 0.0, &z);
		for (j = 2; j <= i; j++) {
		    phgMatVec(MAT_OP_N, y[j], NULL, v_ptr[j], 1.0, &z);
		}
	    }
	    else {
		int ak;

		phgMatVec(MAT_OP_N, y[1], NULL, v_ptr[1], 0.0, &z);
		for (j = 2; j <= mk; j++) {
		    phgMatVec(MAT_OP_N, y[j], NULL, v_ptr[j], 1.0, &z);
		}

		if (itr_out < augk)
		    ak = itr_out;
		else
		    ak = augk;

		for (j = 1; j <= ak; j++) {
		    phgMatVec(MAT_OP_N, y[j + mk], NULL, z_ptr[j - 1], 1.0, &z);
		}
	    }
	}

	/* store z, z[i] = z  */
	zi = itr_out % augk;
	phgVecCopy(z, &z_ptr[zi]);

        if (r_norm <= tol) {
            residual = r_norm;

            if (verb > 0)
                phgPrintf("GMRES iteration %d ||b-Ax||_2 = %e  rel.res %e\n", nits, (double)residual,
                        (double)(residual / r0_norm));
            break;
        }
        else {
            phgVecCopy(b, &z);
            phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &z);
        }

        /*  out itr */
        itr_out++;
    }

    if (res != NULL) *res = residual;

    phgVecDestroy(&z);

    /* z */
    for (i = 0; i < augk; i++)
	phgVecDestroy(&z_ptr[i]);
    phgFree(z_ptr);

    /* restart */
    restart = mk + augk;
    for (i = 0; i < restart + 2; i++) {
	phgVecDestroy(&v_ptr[i]);

        if (i > 0) phgVecDestroy(&flv[i]);
    }

    phgFree(v_ptr);
    phgFree(flv);

    phgFree(H);
    phgFree(HH);
    phgFree(h_data);
    phgFree(hh_data);
    phgFree(y);

    return nits;
}

static int
Init(SOLVER *solver)
{
    phgInfo(2, "Create %s solver\n", solver->oem_solver->name);
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));

    if (solver->oem_solver == SOLVER_PCG ||
	solver->oem_solver == SOLVER_BiCGstab) {
	_t->pc_type = pcg_pc_type;
	_t->pc_opts = (pcg_pc_opts == NULL ? NULL : strdup(pcg_pc_opts));
	solver->pc_option_flag = phgOptionsIfUsed("-pcg_pc_type");
    }
    else if (solver->oem_solver == SOLVER_GMRES) {
	_t->lgmres_augk = lgmres_augk;
	_t->gmres_type = gmres_type;
	_t->pc_type = gmres_pc_type;
	_t->pc_opts = (gmres_pc_opts == NULL ? NULL : strdup(gmres_pc_opts));
	_t->gmres_restart = gmres_restart;
	solver->pc_option_flag = phgOptionsIfUsed("-gmres_pc_type");
    }
    else {
	assert(solver->oem_solver == SOLVER_PreOnly);
	_t->pc_type = preonly_pc_type;
	_t->pc_opts = (preonly_pc_opts == NULL ?
		       NULL : strdup(preonly_pc_opts));
	solver->pc_option_flag = phgOptionsIfUsed("-preonly_pc_type");
    }

    /* parameters for block Jacobi preconditioner */
    _t->bjacobi_solver_id = bjacobi_solver_id;
    if (bjacobi_solver_opts != NULL)
	_t->bjacobi_solver_opts = strdup(bjacobi_solver_opts);

    return 0;
}


static int
Create(SOLVER *solver)
{
    /* Set up PC if not already set */

    if (solver->pc == NULL) {
	if (pc_list[_t->pc_type] != NULL) {
	    phgSolverSetPC_(solver, solver, pc_list[_t->pc_type]);
	}
	else {
	    /* Create PC solver */
	    SOLVER *pc_solver;
	    phgOptionsPush();
	    phgSolverSetDefaultSuboptions();
	    phgOptionsSetOptions
		("-solver_rtol 0 -solver_atol 0 -solver_btol 0"
		 " -solver_maxit 1 +solver_warn_maxit");
	    phgOptionsSetOptions(_t->pc_opts);
	    pc_solver = phgSolverMat2Solver(SOLVER_DEFAULT, solver->mat);
	    phgOptionsPop();
	    if (pc_solver->oem_solver == solver->oem_solver) {
		/* avoid regression */
		((OEM_DATA *) pc_solver->oem_data)->pc_type = 0;
	    }
	    phgSolverSetPC_(solver, pc_solver, NULL);
	    solver->pc->solver_owned = pc_solver;	/* destroy afterwards */
	}
    }

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    phgInfo(2, "Destroy Solver\n");
    if (_t->bjacobi_solver != NULL)
	phgSolverDestroy(&_t->bjacobi_solver);
    if (_t->pc_opts != NULL)
	phgFree(_t->pc_opts);
    if (_t->bjacobi_solver_opts != NULL)
	phgFree(_t->bjacobi_solver_opts);
    phgFree(solver->oem_data);

    return 0;
}

static int
PCG_Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    INT nits;
    FLOAT residual = 1.0;
    VEC *B;

    B = phgMapCreateVec(x->map, 0);
    B->data = solver->rhs->data;
    B->offp_data = solver->rhs->offp_data;
    assert(solver->rhs->assembled);
    B->nvec = 1;
    nits = solver_pcg(solver->mat, B, &x, solver->pc, solver->maxit,
		      solver->rtol, solver->btol, solver->atol,
		      &residual, &_t->warned, solver->monitor);
    B->data = NULL;
    B->offp_data = NULL;
    if (destroy) {
	phgMatDestroy(&solver->mat);
	solver->sent_to_oem = TRUE;
    }
    phgVecDestroy(&B);

    solver->residual = residual;
    return solver->nits = nits;
}

static int
BiCGstab_Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    INT nits;
    FLOAT residual = 1.0;
    VEC *B;

    B = phgMapCreateVec(x->map, 0);
    B->data = solver->rhs->data;
    B->offp_data = solver->rhs->offp_data;
    assert(solver->rhs->assembled);
    B->nvec = 1;
    nits = solver_BiCGstab(solver->mat, B, &x, solver->pc, solver->maxit,
			   solver->rtol, solver->btol, solver->atol,
			   &residual, solver->monitor);
    B->data = NULL;
    B->offp_data = NULL;
    if (destroy) {
	phgMatDestroy(&solver->mat);
	solver->sent_to_oem = TRUE;
    }
    phgVecDestroy(&B);

    solver->residual = residual;
    return solver->nits = nits;
}

static int
GMRES_Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    INT nits;
    FLOAT residual = 1.0;
    VEC *B;

    B = phgMapCreateVec(x->map, 0);
    B->data = solver->rhs->data;
    B->offp_data = solver->rhs->offp_data;
    assert(solver->rhs->assembled);
    B->nvec = 1;
    nits = gmres_funcs[_t->gmres_type](
		solver->mat, B, &x, solver->pc, _t->gmres_restart,
		_t->lgmres_augk, solver->maxit, solver->rtol,
		solver->btol, solver->atol, &residual,
		solver->monitor);
    B->data = NULL;
    B->offp_data = NULL;
    if (destroy) {
	phgMatDestroy(&solver->mat);
	solver->sent_to_oem = TRUE;
    }
    phgVecDestroy(&B);

    solver->residual = residual;

    return solver->nits = nits;
}

static int
PreOnly_Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    SOLVER_PC *pc = solver->pc;

    assert(pc != NULL);
    pc->pc_proc(pc->ctx, solver->rhs, &x);

    if (destroy) {
	phgMatDestroy(&solver->mat);
	solver->sent_to_oem = TRUE;
    }

    solver->residual = 0.0;

    return solver->nits = 0;
}

/*--------------------------------------------------------------------*/

#define SetPC			NULL

OEM_SOLVER phgSolverPCG_ = {
    "PCG", PCGRegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, PCG_Solve, NULL, NULL, NULL,
    M_SPD, TRUE, TRUE, TRUE, TRUE
};

OEM_SOLVER phgSolverBiCGstab_ = {
    "BiCGstab", /*PCGRegisterOptions */ NULL,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, BiCGstab_Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, TRUE, TRUE
};

OEM_SOLVER phgSolverGMRES_ = {
    "GMRES", GMRESRegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, GMRES_Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, TRUE, TRUE
};

OEM_SOLVER phgSolverPreOnly_ = {
    "PreOnly", PreOnlyRegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, PreOnly_Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, TRUE, TRUE
};
