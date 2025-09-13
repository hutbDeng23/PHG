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

/* This file implements the ASP solver (Auxiliary Space Preconditioner).
 *
 * Note: this file was originated from the file solver-aps.c which solves the
 * following Poisson equation:
 *	div(alpha grad) u + beta u = f
 * with alpha > 0 and beta >= 0, discretized using H1-conforming elements
 *
 * $Id: solver-asp.c,v 1.17 2022/09/21 02:35:28 zlb Exp $
 */

#include "phg.h"

#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

#define PENALIZE_BDRY 0		/* 0: penalize diagonal rows,
				 * 1: set RHS of bdry rows to 0. */

typedef struct {
    DOF		*alpha, *beta;	/* The coefficients alpha and beta */
    MAT		*P;		/* The transfer matrix */
    MAT		*A;		/* The Poisson matrix */
    BUILD_FUNC	build_P;
    BUILD_FUNC	build_A;

    DOF		*dof_P;		/* DOF for the aux space */
    SOLVER	*auxsolver, *smoother;

    BOOLEAN	assembled;

    char	*cycle_type;
    char	*dof_type;	/* DOF_TYPE for the auxspace */
    int		aux_solver_id;
    char	*aux_solver_opts;
    char	*smoother_opts;

    FLOAT	constant_alpha, constant_beta;
} OEM_DATA;

static char	*cycle_type = "10";
static char	*dof_type = NULL;

/* The  aux solver options */
static int	aux_solver_id = -1;
static char	*aux_solver_opts = NULL;

/* Options for the smoother.
 * If it's NULL then the built-in GS smoother is used, otherwise a solver
 * with the given options is constructed and used as the smoother. */
static char	*smoother_opts = NULL;

static FLOAT	constant_alpha = _F(-1.0);
static FLOAT	constant_beta = _F(-1.0);

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

/*---------------------------- auxiliary functions -------------------------*/

void
phgSolverASPSetCoefficients(SOLVER *solver, DOF *alpha, DOF *beta)
/* defines the coefficients alpha and beta (NULL means 1) */
{
    if (solver->oem_solver != SOLVER_ASP) {
	phgPrintf("*** WARNING from %s: not an ASP solver, do nothing.\n",
			__func__);
	return;
    }

    if (_t->assembled){
        phgWarning("phgSolverASPSetCoefficients() should be \
                 called before phgSolverAssemble()!\n");
        return;
    }

    if (_t->alpha != NULL)
	phgDofFree(&_t->alpha);
    if (_t->beta != NULL)
	phgDofFree(&_t->beta);

    if (alpha != NULL){
	_t->alpha = phgDofCopy(alpha, NULL, NULL, NULL);
    }
    if (beta != NULL){
	_t->beta = phgDofCopy(beta, NULL, NULL, NULL);
    }
   
    if (_t->A != NULL)
        phgMatDestroy(&_t->A);

    return ;
}

void
phgSolverASPSetConstantCoefficients(SOLVER *solver, FLOAT a, FLOAT b)
{
    GRID *g = solver->mat->rmap->dofs[0]->g;

    if (solver->oem_solver != SOLVER_ASP) {
	phgPrintf("*** WARNING from %s: not ASP solver, do nothing.", __func__);
	return;
    }

    if (_t->assembled){
        phgWarning("phgSolverASPSetConstantCoefficients() should be \
                    called before phgSolverAssemble() %d\n", _t->assembled);
        return ;
    }

    if (_t->alpha != NULL)
	phgDofFree(&_t->alpha);
    if (_t->beta != NULL)
	phgDofFree(&_t->beta);

    _t->alpha = phgDofNew(g, DOF_CONSTANT, 1, "alpha", DofNoAction);
    phgDofSetDataByValue(_t->alpha, a);

    if (b != 0.0) {
	_t->beta = phgDofNew(g, DOF_CONSTANT, 1, "beta", DofNoAction);
	phgDofSetDataByValue(_t->beta, b);
    }

    if (_t->A != NULL)
        phgMatDestroy(&_t->A);

    return;
}

static MAT *
build_P0(SOLVER *solver)
/*----------------------------------------------------------------------------
 * Builds the transfer matrix P (H1->H1).
 *
 * We denote by {\phi} and {\psi}, respectively, the finite element bases of
 * the original and aux spaces, then the elements of the transfer matrix P
 * is determined by representing \phi with \psi:
 *		\psi_j = \sum_k P_{k,j} \phi_k
 *--------------------------------------------------------------------------*/
{
    MAP *map = solver->mat->rmap;
    DOF *u;
    DOF_TYPE *type;
    ELEMENT *e;
    GRID *g;

    assert(/*map->ndof == 1 &&*/ map->dofs != NULL);

    u = map->dofs[0];
    g = u->g;
    assert(/*DofFESpace(u) == FE_H1 &&*/ u->dim == 1);

    if (!DofIsHP(u) || u->hp->max_order == u->hp->min_order ||
	_t->dof_type != NULL) {
	if (_t->dof_type == NULL) {
	    int order;
	    order = DofIsHP(u) ? u->hp->min_order : u->type->order;
	    if (!DofIsHP(u) && u->type->name[0] == 'P')
		order -= (order > 1 ? 1 : 0);
	    type = DOF_Pn[order];
	}
	else {
	    phgOptionsPush();
	    phgOptionsSetHandler("-dof_type", _t->dof_type);
	    type = DOF_DEFAULT;
	    phgOptionsPop();
	}
	_t->dof_P = phgDofNew(g, type, 1, "dof_P", DofNoData);
    }
    else {
	HP_TYPE *hp;

	hp = phgHPNew(g, HP_HB);
	ForAllElements(g, e)
	    e->hp_order = u->hp->elem_order[e->index];
	phgHPSetup(hp, FALSE);
	_t->dof_P = phgHPDofNew(g, hp, 1, "u_h", DofNoData);
	phgHPFree(&hp);
    }

    return phgDofGetTransferMatrix(_t->dof_P, u, NULL, NULL);
}

static MAT *
build_poisson(SOLVER *solver)
/* builds the Poisson matrix:
 * 	-div(alpha grad) + beta
 * Note:
 *    1. P^tAP is used if alpha == NULL
 *    2. beta == NULL means beta == 0.0 */
{
    MAT *mat;
    DOF *u, *alpha, *beta;
    int i, j, n, N_P;
    GRID *g;
    ELEMENT *e;
    FLOAT *A, *A0, *A1, d;
    INT *I, *J;
    BTYPE *b, DB_mask;

    if (_t->alpha == NULL) {
	/* set mat = P^t A P (PtAP) */
	MAT *T;
	T = phgMatMat(MAT_OP_T, MAT_OP_N, 1.0, _t->P, solver->mat, 0.0, NULL);
	mat = phgMatMat(MAT_OP_N, MAT_OP_N, 1.0, T, _t->P, 0.0, NULL);
	phgMatDestroy(&T);
	return mat;
    }

    u = solver->mat->rmap->dofs[0];
    g = u->g;
    alpha = _t->alpha;
    beta = _t->beta;

    /* TODO: check DOF_TYPE of solver->dofs[] for consistency */
    assert(beta == NULL || beta->dim == 1);

    DB_mask = solver->mat->rmap->dofs[0]->DB_mask;

    mat = phgMapCreateMat(_t->P->cmap, _t->P->cmap);

    /* number of basis functions in an element */
    if (!DofIsHP(_t->dof_P))
	N_P = _t->dof_P->type->nbas;
    else
	N_P = _t->dof_P->hp->info->types[_t->dof_P->hp->max_order]->nbas;

    A0 = phgAlloc(N_P * N_P * sizeof(*A0));

    A = phgAlloc(N_P * N_P * sizeof(*A));
    A1 = phgAlloc(N_P * sizeof(*A1));
    I = phgAlloc(N_P * sizeof(*I));
    J = phgAlloc(N_P * sizeof(*J));
    b = phgAlloc(N_P * sizeof(*b));

    /* Note: if beta is constant then the element mass matrices only vary
       by element volumes */
    if (!DofIsHP(u) && (beta == NULL || beta->type == DOF_CONSTANT)) {
	assert(g->roots != NULL);
	e = g->roots;
	d = (beta == NULL ? 0.0 : *beta->data) / phgGeomGetVolume(g, e);
	for (i = 0; i < N_P; i++)
	    for (j = 0; j <= i; j++)
		A0[i * N_P + j] = A0[j * N_P + i] =
		    d * phgQuadBasDotBas(e, _t->dof_P, j, _t->dof_P, i, -1);
    }

    ForAllElements(g, e) {

	if (DofIsHP(u))
	    N_P = DofNBas(_t->dof_P, e);

	/* Matrix *A */
	if (DofIsHP(u) || (beta != NULL && beta->type != DOF_CONSTANT)) {
	    /* compute \int beta * phi_j * phi_i */
	    for (i = 0; i < N_P; i++)
		for (j = 0; j <= i; j++)
		    A[i * N_P + j] =
			phgQuadBasABas(e, _t->dof_P, j, beta, _t->dof_P, i, -1);
	}
	else {
	    d = phgGeomGetVolume(g, e);
	    for (i = 0; i < N_P; i++)
		for (j = 0; j <= i; j++)
		    A[i * N_P + j] = A0[i * N_P + j] * d;
	}

	for (i = 0; i < N_P; i++) {
	    I[i] = phgMapE2L(mat->rmap, 0, e, i);
	    b[i] = phgDofGetElementBoundaryType(_t->dof_P, e, i*_t->dof_P->dim);
	    for (j = 0; j <= i; j++) {
		A[j * N_P + i] = (A[i * N_P + j] +=
		    phgQuadGradBasAGradBas(e, _t->dof_P, j, alpha, _t->dof_P,
						i, -1));
	    }
	}

	for (i = 0; i < N_P; i++) {
	    if (b[i] & DB_mask) {
#if 0
		phgMatAddEntry(mat, I[i], I[i], 1.0);
#else
		phgMatAddEntry(mat, I[i], I[i], A[i * N_P + i]);
#endif
	    }
	    else {
		for (j = 0, n = 0; j < N_P; j++) {
		    if (b[j] & DB_mask)
			continue;
		    J[n] = I[j];
		    A1[n] = A[i * N_P + j];
		    n++;
		}
		phgMatAddEntries(mat, 1, I + i, n, J, A1);
	    }
	}
    }	/* end-ForAllElement */

    phgMatAssemble(mat);

    phgFree(A);
    phgFree(A0);
    phgFree(A1);
    phgFree(I);
    phgFree(J);
    phgFree(b);

    return mat;
}

static BUILD_FUNC build_P = build_P0;
static BUILD_FUNC build_A = build_poisson;

DOF *
phgSolverASPGetAuxDof(SOLVER *solver)
{
    if (solver->oem_solver != SOLVER_ASP) {
	phgPrintf("*** WARNING from %s: not an ASP solver.\n", __func__);
	return NULL;
    }

    return _t->dof_P;
}

BUILD_FUNC
phgSolverASPSetBuildP(BUILD_FUNC func)
{
    BUILD_FUNC old_func;

    old_func = build_P;
    build_P = (func == NULL ? build_P0 : func);

    return old_func;
}

BUILD_FUNC
phgSolverASPSetBuildA(BUILD_FUNC func)
{
    BUILD_FUNC old_func;

    old_func = build_A;
    build_A = (func == NULL ? build_poisson : func);

    return old_func;
}

/*----------------------------------------------------------------------*/

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nThe ASP solver options:", "\n", "asp");

    phgOptionsRegisterString("-asp_cycle_type", "ASP cycle type", &cycle_type);
    phgOptionsRegisterString("-asp_dof_type", "DOF type for the aux space",
				&dof_type);

    phgOptionsRegisterString("-asp_smoother_opts",
		"Smoother options (will use an external solver if it's set)",
		&smoother_opts);

    /* options for the aux solver */

#if USE_HYPRE && HYPRE_VERSION_MAJOR >= 2
    /* set default aux solver to Hypre BoomerAMG */
    for (aux_solver_id = 0;
	 phgSolverNames[aux_solver_id] != NULL
	 && strcmp(phgSolverNames[aux_solver_id], "hypre") != 0;
	 aux_solver_id++);
    assert(phgSolverNames[aux_solver_id] != NULL);
#endif	/* USE_HYPRE && HYPRE_VERSION_MAJOR >= 2 */

#if USE_TRILINOS
    if (aux_solver_id < 0 || phgSolverNames[aux_solver_id] == NULL) {
	for (aux_solver_id = 0;
	     phgSolverNames[aux_solver_id] != NULL
	     && strcmp(phgSolverNames[aux_solver_id], "trilinos") != 0;
	     aux_solver_id++);
	assert(phgSolverNames[aux_solver_id] != NULL);
    }
#endif	/* USE_TRILINOS */

    /* Fall back to first available solver */
    if (aux_solver_id < 0)
	aux_solver_id = 0;

    if (phgSolverNames[aux_solver_id] == NULL)
	aux_solver_id = 0;

    phgOptionsRegisterKeyword("-asp_aux_solver", "Auxiliary solver",
			      phgSolverNames, &aux_solver_id);
    phgOptionsRegisterString("-asp_aux_solver_opts", "Auxiliary solver options",
			      &aux_solver_opts);

    /* other options */
    phgOptionsRegisterFloat("-asp_alpha", "Constant alpha value",
			    &constant_alpha);
    phgOptionsRegisterFloat("-asp_beta", "Constant beta value",
			    &constant_beta);

    return 0;
}

static int
Init(SOLVER *solver)
{
    char *s;

    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));

    _t->build_P = build_P;
    _t->build_A = build_A;

    /* copy relevant cmdline arguments to OEM_DATA */

    _t->cycle_type = strdup(cycle_type);

    if (dof_type != NULL) {
	/* remove leading spaces */
	s = dof_type;
	while (*s != '\0' && isspace(*s))
	    s++;
	/* copy the string only when it's not empty */
	if (*s != '\0') {
	    s = _t->dof_type = strdup(s);
	    /* remove trailing spaces */
	    s += strlen(s);
	    while (--s >= _t->dof_type && isspace(*s));
	    *(++s) = '\0';
	}
    }

    _t->aux_solver_id = aux_solver_id;
    if (aux_solver_opts != NULL) {
	/* remove leading spaces */
	s = aux_solver_opts;
	while (*s != '\0' && isspace(*s))
	    s++;
	/* copy the string only when it's not empty */
	if (*s != '\0') {
	    s = _t->aux_solver_opts = strdup(s);
	    /* remove trailing spaces */
	    s += strlen(s);
	    while (--s >= _t->aux_solver_opts && isspace(*s));
	    *(++s) = '\0';
	}
    }

    if (smoother_opts != NULL) {
	/* remove leading spaces */
	s = smoother_opts;
	while (*s != '\0' && isspace(*s))
	    s++;
	/* copy the string only when it's not empty */
	if (*s != '\0') {
	    s = _t->smoother_opts = strdup(s);
	    /* remove trailing spaces */
	    s += strlen(s);
	    while (--s >= _t->smoother_opts && isspace(*s));
	    *(++s) = '\0';
	}
    }

    _t->alpha=NULL;
    _t->beta=NULL;
    _t->A=NULL;
    _t->assembled=FALSE;

    _t->constant_alpha = constant_alpha;
    _t->constant_beta = constant_beta;

    return 0;
}

static int
Create(SOLVER *solver)
{
    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }
    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    if (_t->alpha != NULL)
	phgDofFree(&_t->alpha);

    if (_t->beta != NULL)
	phgDofFree(&_t->beta);

    if (_t->P != NULL)
	phgMatDestroy(&_t->P);
    if (_t->A != NULL)
	phgMatDestroy(&_t->A);
    if (_t->dof_P != NULL)
	phgDofFree(&_t->dof_P);

    if (_t->auxsolver != NULL)
	phgSolverDestroy(&_t->auxsolver);

    if (_t->smoother != NULL)
	phgSolverDestroy(&_t->smoother);

    phgFree(_t->cycle_type);
    phgFree(_t->dof_type);
    phgFree(_t->aux_solver_opts);
    phgFree(_t->smoother_opts);

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    double t = 0.0, t1;

    if (_t->assembled)
	return 0;

    if (_t->alpha == NULL && _t->A == NULL) {
	/* use coefficients constant_alpha and constant_beta if set */

	/* adjust constant values of alpha and beta */
	assert(_t->constant_alpha != 0.0);
	if (_t->constant_alpha > 0.0 && _t->constant_beta < 0.0)
	    _t->constant_beta = 0.0;
	if (_t->constant_alpha < 0.0 && _t->constant_beta >= 0.0)
	    _t->constant_alpha = 1.0;

	if (_t->constant_alpha > 0.0) {
	    if (solver->monitor)
		phgPrintf("ASP: constant coefficients: alpha=%lg, beta=%lg\n",
			(double)_t->constant_alpha, (double)_t->constant_beta);
	    phgSolverASPSetConstantCoefficients(solver,
				_t->constant_alpha, _t->constant_beta);
	}
	else {
	    if (solver->monitor)
		phgPrintf("ASP: using P^tAP as the coarse grid operator.\n");
	}
    }

    _t->assembled = TRUE;

    /* preconditioners */
    if (solver->monitor)
	t = phgGetTime(NULL);

    if (_t->P == NULL) {
	_t->P = (*_t->build_P)(solver);
	if (solver->monitor) {
	    phgPrintf("build_P time: %lg\n", (t1 = phgGetTime(NULL)) - t);
	    t = t1;
	}
    }

    if (_t->A == NULL) {
	_t->A = (*_t->build_A)(solver);
	if (solver->monitor) {
	    phgPrintf("build_A time: %lg\n", (t1 = phgGetTime(NULL)) - t);
	    t = t1;
	}
    }

    /* build solver for the smoother */
    if (_t->smoother_opts != NULL) {
	phgOptionsPush();
	phgSolverSetDefaultSuboptions();
	phgOptionsSetOptions("-solver_rtol 0. -solver_btol 0. -solver_atol 0. "
			     "-solver_maxit 1");
	phgOptionsSetOptions(_t->smoother_opts);
	_t->smoother = phgMat2Solver(SOLVER_DEFAULT, solver->mat);
	_t->smoother->warn_maxit = FALSE;
	phgVecDestroy(&_t->smoother->rhs);
	_t->smoother->rhs_updated = TRUE;
	phgOptionsPop();
    }

    /* build aux solver */
    phgOptionsPush();
    phgSolverSetDefaultSuboptions();
    phgOptionsSetKeyword("-solver", phgSolverNames[_t->aux_solver_id]);
    phgOptionsSetOptions("-solver_rtol 0. "
			 "-solver_btol 0. "
			 "-solver_atol 0. "
			 "-solver_maxit 1");
    /* some default solver options */
#if USE_HYPRE
    phgOptionsSetOptions("-hypre_solver boomeramg -hypre_pc none");
    phgOptionsSetKeyword("-hypre_amg_coarsen_type", "hmis");
    phgOptionsSetInt("-hypre_amg_max_levels", (INT)25);
    /* HYPRE_BoomerAMGSetAggNumLevels 1 (default) */
    phgOptionsSetKeyword("-hypre_amg_relax_type", "gs-h-symmetric");

    /* Gauss-Seidel instead of direct solver on coarsest grid */
    phgOptionsSetKeyword("-hypre_amg_coarsest_relax_type", "gs-h-forward");

    /* HYPRE_BoomerAMGSetNumSweeps 1 (done by solver-hypre.c) */
    phgOptionsSetFloat("-hypre_amg_strength", (FLOAT)0.25);	/* theta */
#endif	/* USE_HYPRE */
#if USE_TRILINOS
    phgOptionsSetOptions("-aztec_solver cg -aztec_pc ml");
    phgOptionsSetOptions("-aztec_disable_errmsg");
# if HAVE_FEENABLEEXCEPT
    phgOptionsSetOptions("-aztec_disable_fpetrap");
# endif	/* HAVE_FEENABLEEXCEPT */
#endif	/* USE_TRILINOS */
    /* append user options */
    phgOptionsSetOptions(_t->aux_solver_opts);
    _t->auxsolver = phgMat2Solver(SOLVER_DEFAULT, _t->A);
    _t->auxsolver->warn_maxit = FALSE;
    phgVecDestroy(&_t->auxsolver->rhs);
    _t->auxsolver->rhs_updated = TRUE;
    phgOptionsPop();

#if PENALIZE_BDRY
{
    INT i;
    MAT_ROW *row;

    /* penalize boundary rows */
    for (i = 0; i < _t->A->rmap->nlocal; i++) {
	row = _t->A->rows + i;
	assert(row->ncols > 0);
	if (row->ncols == 1 && row->cols[0] == i)
	    row->data[0] = DBL_MAX;
    }
}
#endif	/* PENALIZE_BDRY */


    return 0;
}

#if USE_OMP
/* TODO: balance # of non-zeros between threads */
#define ThreadRange(n, tid, start, end)					\
    {									\
	int d = (n) / phgMaxThreads, r = (n) - d * phgMaxThreads;	\
	start = d * (tid) + ((tid) < r ? (tid) : r);			\
 	end = start + d + ((tid) < r ? 1 : 0);				\
    }
#endif	/* USE_OMP */

static VEC *
gauss_seidel(MAT *mat, VEC *rhs, VEC **x_ptr, INT maxit, FLOAT tol,
	     FLOAT *res_ptr, INT verb)
/* scaled symmetric l1-GS relaxation */
{
    VEC *x;
    INT it, i, j, n, *pc = NULL;
    FLOAT *pd = NULL, a, b, res = FLOAT_MAX, *diag1;
#if USE_MPI || USE_OMP
    FLOAT *tmp = NULL;
# if USE_MPI
    INT *pc_offp = NULL;
    FLOAT *pd_offp = NULL, *offp_data = NULL;
# endif	/* USE_MPI */
# if USE_OMP
    FLOAT *xsave, *resp;
# endif	/* USE_OMP */
#endif	/* USE_MPI || USE_OMP */

    x = (x_ptr == NULL ? NULL : *x_ptr);
    if (x == NULL) {
	x = phgMapCreateVec(mat->cmap, rhs->nvec);
	if (x_ptr != NULL)
	    *x_ptr = x;
    }

    if (mat->type == PHG_UNPACKED)
	phgMatPack(mat);

    assert(mat->type == PHG_PACKED);

    if (mat->diag1 == NULL) {
	/* setup off-diagonal L1 norm1 */
	diag1 = mat->diag1 = phgAlloc(mat->rmap->nlocal * sizeof(*mat->diag1));
	if (mat->diag == NULL)
	    phgMatSetupDiagonal(mat);
#if USE_OMP
	if (phgMaxThreads == 1) {
#endif	/* USE_OMP */
#if USE_MPI
	    pd_offp = mat->packed_data + mat->packed_ind[mat->rmap->nlocal];
#endif	/* USE_MPI */
	    for (i = 0; i < mat->rmap->nlocal; i++) {
		a = Fabs(mat->diag[i]);
#if USE_MPI
		/* off-process columns */
		j = mat->rmap->nlocal + i;
		n = (INT)(mat->packed_ind[j + 1] - mat->packed_ind[j]);
		for (j = 0; j < n; j++)
		    a += Fabs(pd_offp[j]);
		pd_offp += n;
#endif	/* USE_MPI */
		assert(a != 0.);
		diag1[i] = 1.0 / a;
	    }
#if USE_OMP
	}
	else {
#if USE_MPI
#pragma omp parallel private(i, j, pd_offp, pc, pd, n, a)
#else	/* USE_MPI */
#pragma omp parallel private(i, j, pc, pd, n, a)
#endif	/* USE_MPI */
{
	    int k;
	    INT startind, endind, l;

	    ThreadRange(x->map->nlocal, phgThreadId, startind, endind)

	    pc = mat->packed_cols + mat->packed_ind[startind];
	    pd = mat->packed_data + mat->packed_ind[startind];
#if USE_MPI
	    j = mat->rmap->nlocal + startind;
	    pd_offp = mat->packed_data + mat->packed_ind[j];
#endif	/* USE_MPI */

# pragma omp for schedule(static)
	    for (k = 0; k < phgMaxThreads; k++) {
		for (i = startind ; i < endind; i++) {
		    a = Fabs(mat->diag[i]);
#if USE_MPI
		    /* off-process columns */
		    j = mat->rmap->nlocal + i;
		    n = (INT)(mat->packed_ind[j + 1] - mat->packed_ind[j]);
		    for (j = 0; j < n; j++)
			a += Fabs(pd_offp[j]);
		    pd_offp += n;
#endif	/* USE_MPI */

		    /* in-process but off-thread columns */
		    n = (INT)(mat->packed_ind[i + 1] - mat->packed_ind[i]);
		    for (j = 0; j < n; j++){
			if ((l = pc[j]) < startind || l >= endind)
			    a += Fabs(pd[j]);
		    }
		    pc += n;
		    pd += n;
    
		    assert(a != 0.);
		    diag1[i] = 1.0 / a;
		} /* i loop */
	    } /* k loop */
} /* omp  parallel*/
	}
#endif	/* USE_OMP */
    } /* mat->diag1 == NULL */
    diag1 = mat->diag1;

#if USE_MPI
    tmp = phgAlloc(x->map->nlocal * sizeof(*tmp));
    if (x->map->nprocs > 1)
	offp_data = phgAlloc(mat->cinfo->rsize * sizeof(*offp_data));
#elif USE_OMP
    tmp = phgAlloc(x->map->nlocal * sizeof(*tmp));
#endif	/* USE_MPI */

#if USE_OMP
    resp = phgAlloc(phgMaxThreads * sizeof(*resp));
    xsave = phgAlloc(x->map->nlocal * sizeof(*xsave));
#endif

#if 0
    double t0 = phgGetTime(NULL);
    maxit = 1000;
#endif

    for (it = 0; it < maxit; it++) {
#if USE_OMP
	if (phgMaxThreads == 1) {
#endif
	/* forward scan */
#if USE_MPI
	    if (x->map->nprocs > 1) {
		phgMapScatterBegin(mat->cinfo, x->nvec, x->data, offp_data);
		phgMapScatterEnd  (mat->cinfo, x->nvec, x->data, offp_data);
	    }
#endif	/* USE_MPI */
	    pc = mat->packed_cols;
	    pd = mat->packed_data;
#if USE_MPI
	    pc_offp = mat->packed_cols + mat->packed_ind[mat->rmap->nlocal];
	    pd_offp = mat->packed_data + mat->packed_ind[mat->rmap->nlocal];
#endif	/* USE_MPI */
	    for (i = 0; i < mat->rmap->nlocal; i++) {
		a = rhs->data[i];
#if USE_MPI
		/* off-process columns */
		j = mat->rmap->nlocal + i;
		n = (INT)(mat->packed_ind[j + 1] - mat->packed_ind[j]);
		for (j = 0; j < n; j++)
		   a -= pd_offp[j] * offp_data[pc_offp[j]];
		tmp[i] = a;
		pc_offp += n;
		pd_offp += n;
#endif	/* USE_MPI */

		/* in-process columns */
		n = (INT)(mat->packed_ind[i + 1] - mat->packed_ind[i]); 
		for (j = 0; j < n; j++)
		    a -= pd[j] * x->data[pc[j]];
		x->data[i] += a * diag1[i];
		pc += n;
		pd += n;
	    }

	    /* backward scan (note: don't exchange and use the new off_p data
	     * here, or the convergence will be slower) */
	    res = 0.0;
	    for (i = mat->rmap->nlocal - 1; i >= 0; i--) {
#if USE_MPI
		a = tmp[i];
#else	/* USE_MPI */
		a = rhs->data[i];
#endif	/* USE_MPI */

		/* in-process columns */
		n = (INT)(mat->packed_ind[i + 1] - mat->packed_ind[i]); 
		pc -= n;
		pd -= n;
		for (j = n - 1; j >= 0; j--)
		    a -= pd[j] * x->data[pc[j]];
		b = x->data[i];
		x->data[i] += a * diag1[i];
		b = Fabs(x->data[i] - b);
		if (res <= b)
		    res = b;
	   }
#if USE_OMP
	}
	else {
	    /* forward scan */
#if USE_MPI
	    if (x->map->nprocs > 1) {
		phgMapScatterBegin(mat->cinfo, x->nvec, x->data, offp_data);
		phgMapScatterEnd  (mat->cinfo, x->nvec, x->data, offp_data);
	    }
#endif	/* USE_MPI */
	    if (x->map->nlocal > 0)
		memcpy(xsave, x->data, x->map->nlocal * sizeof(*xsave));

#if USE_MPI
# pragma omp parallel default(shared) \
	private(a, b, i, j, pc_offp, pd_offp, n, pc, pd)
#else	/* USE_MPI */
# pragma omp parallel default(shared) private(a, b, i, j, n, pc, pd)
#endif	/* USE_MPI */
{ 
	    int k;
	    INT startind, endind, l;

	    ThreadRange(x->map->nlocal, phgThreadId, startind, endind)

	    pc = mat->packed_cols + mat->packed_ind[startind];
	    pd = mat->packed_data + mat->packed_ind[startind];
#if USE_MPI
	    j = mat->rmap->nlocal + startind;
	    pc_offp = mat->packed_cols + mat->packed_ind[j];
	    pd_offp = mat->packed_data + mat->packed_ind[j];
#endif	/* USE_MPI */

#pragma omp for schedule(static)
	    for (k = 0; k < phgMaxThreads; k++) {
		for (i = startind; i < endind; i++) {
		    a = rhs->data[i];
#if USE_MPI
		    /* off-process columns */
		    j = mat->rmap->nlocal + i;
		    n = (INT)(mat->packed_ind[j + 1] - mat->packed_ind[j]);
		    for (j = 0; j < n; j++)
			a -= pd_offp[j] * offp_data[pc_offp[j]];
		    pc_offp += n;
		    pd_offp += n;
#endif	/* USE_MPI */

		    /* in-process columns */
		    n = (INT)(mat->packed_ind[i + 1] - mat->packed_ind[i]);
		    b = 0.0;
		    for (j = 0; j < n; j++) {
			if ((l = pc[j]) < startind || l >= endind) {
			    /* off-thread column */
			    a -= pd[j] * xsave[l];
			}
			else {
			    /* in-thread column */
			    b -= pd[j] * x->data[l];
			}
		    }
		    tmp[i] = a;
		    x->data[i] += (a + b) * diag1[i];
		    pc += n;
		    pd += n;
		}
	    } /* k loop */

	    resp[phgThreadId] = 0.0;

#pragma omp for schedule(static)
	    for (k = 0; k < phgMaxThreads; k++) {
		for (i = endind - 1; i >= startind; i--) {
		    a = tmp[i];
		    /* in-process columns */
		    n = (INT)(mat->packed_ind[i + 1] - mat->packed_ind[i]);
		    pc -= n;
		    pd -= n;
		    for (j = n - 1; j >= 0; j--) {
			if ((l = pc[j]) >= startind && l < endind)
			    a -= pd[j] * x->data[l];
		    }
		    b = x->data[i];
		    x->data[i] += a * diag1[i];
		    b = Fabs(x->data[i] - b);
		    if (resp[phgThreadId] <= b)
			resp[phgThreadId] = b;
		}
	    }
} /* omp parallel */
	    res=resp[0];
	    for (i = 1; i < phgMaxThreads; i++)
		if (res < resp[i]) 
		    res = resp[i]; 
	} /* if phgMaxThreads == 1 */
#endif
	if (verb > 0)
	    phgPrintf("\tG-S sweep %d, res %le\n", it, (double)res);

	if (res <= tol)
	    break;
    } /* it - loop */

#if 0
    printf("time = %lg\n", phgGetTime(NULL) - t0);
    MPI_Finalize();
    exit(0);
#endif

#if USE_MPI
    phgFree(tmp);
    phgFree(offp_data);
#elif USE_OMP
    phgFree(tmp);
#endif

#if USE_OMP
    phgFree(xsave);
    phgFree(resp);
#endif	/* USE_OMP */

    if (res_ptr != NULL)
	*res_ptr = res;

    return x;
}

static void
mult_prec(SOLVER *solver, VEC *r, VEC *x, const char *cycle)
{
    VEC *r0 = NULL, *r1 = NULL, *x1 = NULL;

    while (*cycle != '\0') {
	if (*cycle < '0' && *cycle > '1')
	    phgError(1, "invalid ASP cycle type string: %s\n", _t->cycle_type);
	switch (*(cycle++) - '0') {
	    case 0:	/* smoothing */
		if (phgVerbosity > 1)
		    phgPrintf("\tASP smoothing\n");
		if (_t->smoother == NULL) {
		    /* external smoother undefined, perform GS smoothing */
		    gauss_seidel(solver->mat, r, &x, 1, 0., NULL, 0);
		    break;
		}
		/* call external smoother */
		_t->smoother->rhs = r;
		_t->smoother->rhs->assembled = TRUE;
		phgSolverVecSolve(_t->smoother, FALSE, x);
		_t->smoother->rhs = NULL;
		break;
	    case 1:	/* space P correction */
		if (phgVerbosity > 1)
		    phgPrintf("\tASP aux space correction\n");
		phgVecCopy(r, &r0);
		phgMatVec(MAT_OP_N, -1.0, solver->mat, x, 1.0, &r0);
		if (x1 == NULL || x1->map != _t->P->cmap) {
		    phgVecDestroy(&x1);
		    x1 = phgMapCreateVec(_t->P->cmap, 1);
		}
		else {
		    bzero(x1->data, x1->map->nlocal * sizeof(*x1->data));
		}
		if (r1 != NULL && r1->map != _t->P->cmap)
		    phgVecDestroy(&r1);
		phgMatVec(MAT_OP_T, 1.0, _t->P, r0, 0.0, &r1);

#if !PENALIZE_BDRY
		{
		    /* For boundary DOF: set RHS of bdry rows to 0 */
		    INT i;
		    MAP *map = _t->auxsolver->mat->rmap;

		    for (i = 0; i < map->bdry_nlocal; i++) {
		    	r1->data[map->bdry[i]] = 0;
		    }
		}		
#endif	/* PENALIZE_BDRY */

		_t->auxsolver->rhs = r1;
		_t->auxsolver->rhs->assembled = TRUE;
		phgSolverVecSolve(_t->auxsolver, FALSE, x1);
		_t->auxsolver->rhs = NULL;
		phgMatVec(MAT_OP_N, 1.0, _t->P, x1, 1.0, &x);
		break;
	}	/* case */
    }	/* while */

    if (r0 != NULL)
	phgVecDestroy(&r0);
    if (r1 != NULL)
	phgVecDestroy(&r1);
    if (x1 != NULL)
	phgVecDestroy(&x1);
}

static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    char *p = _t->cycle_type, *q;
    VEC *y0 = NULL, *x0 = NULL;
    FLOAT res = FLOAT_MAX, ib_norm = 0.0, ires0 = 0.0, tol = 0.0;

    Assemble(solver);

    if (solver->maxit > 0 || solver->btol > 0. || solver->rtol > 0. ||
	solver->atol > 0. || solver->monitor) {
	ib_norm = phgVecNorm2(solver->rhs, 0, NULL);
	if (ib_norm == 0.0) {
	    solver->nits = 0;
	    solver->residual = 0.;
	    bzero(x->data, x->nvec * x->map->nlocal * sizeof(*x->data));
	    return 0;
	}
	tol = solver->btol * ib_norm;
	if (tol < solver->atol)
	    tol = solver->atol;
	ib_norm = 1.0 / ib_norm;	/* inversed rhs norm */
    }

    solver->nits = 0;
    if (solver->monitor) {
	DOF_TYPE *type = _t->A->rmap->dofs[0]->type;
        phgPrintf("======================== ASP =====================\n");
	phgPrintf("Auxiliary space: %s\n", type == NULL ? "hp" : type->name);
	phgPrintf("Auxiliary solver: %s\n", _t->auxsolver->oem_solver->name);
	phgPrintf("Smoother: %s\n", _t->smoother == NULL ?
					"built-in point-wise Gauss-Seiddel" :
					_t->smoother->oem_solver->name);
	phgPrintf("Cycle type: %s\n", _t->cycle_type);
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
	    (double)solver->atol, (double)solver->rtol, (double)solver->btol);
        phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
        phgPrintf("-----   ------------   -------------  ------------\n");
    }
    while (solver->nits < solver->maxit) {
	solver->nits++;

	if ((q = strchr(p, '+')) == NULL) {
	    mult_prec(solver, solver->rhs, x, p);
	}
	else {
	    phgVecCopy(solver->rhs, &y0);
	    phgMatVec(MAT_OP_N, -1.0, solver->mat, x, 1.0, &y0);
	    x0 = phgMapCreateVec(y0->map, y0->nvec);
	    while (*p != '\0') {
		if (q != NULL)
		    *q = '\0';
		bzero(x0->data, sizeof(*x0->data) * x0->map->nlocal);
		mult_prec(solver, y0, x0, p);
		phgVecAXPBY(1.0, x0, 1.0, &x);
		if (q == NULL)
		    break;
		*q = '+';
		q = strchr(p = q + 1, '+');
	    }
	}

	if (solver->rtol > 0.0 || tol > 0.0 || solver->monitor) {
	    phgVecCopy(solver->rhs, &x0);
	    phgMatVec(MAT_OP_N, -1.0, solver->mat, x, 1.0, &x0);
	    res = phgVecNorm2(x0, 0, NULL);
	    if (ires0 == 0.0) {
		ires0 = solver->rtol * res;
		if (tol < ires0)
		    tol = ires0;
		ires0 = (res == 0. ? 1.0 : 1.0 / res);	/* inversed res norm */
	    }
	    if (solver->monitor)
		phgPrintf("% 5d   %12le   %12le   %12le\n", solver->nits,
				(double)res, (double)(res * ires0),
				(double)(res * ib_norm));
	    if (res <= tol)
		break;
	}
    }

    phgVecDestroy(&x0);
    phgVecDestroy(&y0);

    if (destroy) {
	Destroy(solver);
	phgMatDestroy(&solver->mat);
    }

    solver->residual = res/*_b*/;

    return solver->nits;
}

OEM_SOLVER phgSolverASP_ = {
    "ASP", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_SPD, TRUE, TRUE, TRUE, FALSE
};
