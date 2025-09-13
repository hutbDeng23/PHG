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

/* AMS, Auxiliary space Maxwell Solver, solves the following time-harmonic
 * Maxwell equation:
 *
 *	curl(alpha curl) E + beta E = J
 *
 * with alpha > 0 and beta >= 0, discretized using edge elements.
 *
 * $Id: solver-ams.c,v 1.129 2022/09/21 02:35:28 zlb Exp $
 */

#define HYPRE_TEST	0	/* 1 ==> test with Hypre AMS preconditioner */

#include "phg.h"

#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdlib.h>

typedef struct {
    DOF		*alpha, *beta;	/* The coefficients alpha and beta */
    MAT		*G, *Pi[3];	/* The G and Pi matrices */
    MAT		*A_G, *A_Pi;	/* The Poisson matrices */
    DOF		*dof_G, *dof_Pi;/* DOFs for the spaces G and Pi */
    SOLVER	*solver_G, *solver_Pi;

#if HYPRE_TEST
    SOLVER	*solver_hypre;
#endif	/* HYPRE_TEST */
    BOOLEAN	assembled;

    char	*cycle_type;

    int		aux_solver_id;
    char	*aux_solver_opts;

    BOOLEAN	warn_aux, full_order_G;
    FLOAT	default_alpha, default_beta;
} OEM_DATA;

static char *cycle_type = "02120";	/* default cycle type */

static int	aux_solver_id	= -1;
static char	*aux_solver_opts	= NULL;

static int	grad_solver_id	= 0;
static char	*grad_solver_opts	= NULL;

static BOOLEAN	warn_aux = TRUE;
static BOOLEAN	full_order_G = FALSE;
static FLOAT	default_alpha = _F(1.0);
static FLOAT	default_beta = _F(1.0);

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

/*---------------------------- auxiliary functions -------------------------*/

static void
gradient_build(SOLVER *solver, DOF *u_h, DOF *f_h, DOF *beta)
{
    int N = u_h->type->nbas;    /* number of basis functions in an element */
    GRID *g = u_h->g;
    ELEMENT *e;
    int i, j;
    FLOAT A[N][N], B[N], buffer[N], *row;
    INT I[N];

    assert(u_h->dim == 1);

    ForAllElements(g, e) {
	/* Note: only works with piecewise constant beta */
	if (beta != NULL && *DofElementData(beta, e->index) != 0.0) {
	    for (i = 0; i < N; i++) {
		j = phgSolverMapE2L(solver, 0, e, i);
		phgSolverAddMatrixEntry(solver, j, j, 1.0);
		phgSolverAddRHSEntry(solver, j, 0.0);
	    }
	    continue;
	}

	/* compute \int \grad\phi_j \cdot \grad\phi_i making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (j = 0; j <= i; j++)
		A[j][i] = A[i][j] =
		    phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	}

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, NULL,
				  buffer, B+ i, DOF_PROJ_NONE)) {
		row = buffer;
	    }
	    else if (phgDofGetElementBoundaryType(u_h, e, i) & NEUMANN) {
		/* TODO */
		row = NULL;
		phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__);
	    }
	    else {
		/* right hand side: \int f * phi_i */
		B[i] = phgQuadDofDotGradBas(e, f_h, u_h, i, QUAD_DEFAULT);
		row = A[i];
	    }
	    phgSolverAddMatrixEntries(solver, 1, I + i, N, I, row); 
	}
	phgSolverAddRHSEntries(solver, N, I, B);
    }
}

void
phgSolverAMSRemoveGradient(DOF **fptr, DOF *beta, int order)
/* Returns a new DOF which is equal to f - \grad\psi, where \psi is 0 in
 * the region beta != 0, and satisfies the following poisson equation:
 * 	<\grad\psi, \grad\phi> = <f, \grad\phi>
 * in the region beta == 0.
 *
 * beta == NULL means beta == 0 everywhere */
{
    SOLVER *solver;
    DOF *f, *psi, *grad;

    assert(fptr != NULL && *fptr != NULL);

    assert(beta == NULL || beta->type == DOF_P0 || beta->type == DOF_DG0
	   || beta->type == DOF_CONSTANT);

    f = *fptr;

    assert(f->type != DOF_ANALYTIC && f->type != DOF_CONSTANT
		&& DofDim(f) == Dim && order >= 1 && order <= 4);

    psi = phgDofNew(f->g, DOF_Pn[order], 1, "psi", DofNoAction);

    /* Set up the Poisson solver */
    phgOptionsPush();
    phgSolverSetDefaultSuboptions();
    phgOptionsSetKeyword("-solver", phgSolverNames[grad_solver_id]);
    phgOptionsSetOptions("-solver_rtol 1e-16 "
			 "-solver_btol 1e-16 "
			 "-solver_maxit 2000 "
			 "-pcg_pc_type jacobi "
			 "-gmres_pc_type jacobi");
    phgOptionsSetOptions(grad_solver_opts);
    solver = phgSolverCreate(SOLVER_DEFAULT, psi, NULL);
    phgOptionsPop();

    gradient_build(solver, psi, f, beta);
    phgSolverSolve(solver, TRUE, psi, NULL);
    if (solver->monitor)
	phgPrintf("Solve gradient equation: nits = %d, residual = %0.2le ",
			solver->nits, (double)solver->residual);
    phgSolverDestroy(&solver);

    grad = phgDofGradient(psi, NULL, f->type, NULL);
    phgDofFree(&psi);
    phgDofAXPY(-1.0, grad, fptr);
    phgDofFree(&grad);

    return;
}

void
phgSolverAMSSetCoefficients(SOLVER *solver, DOF *alpha, DOF *beta)
/* defines the coefficients alpha and beta (NULL means 1) */
{
    assert(solver->oem_solver == SOLVER_AMS);

    if (_t->assembled){
        phgWarning("phgSolverAMSSetCoefficients() should be \
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
   
    if (_t->A_G != NULL)
        phgMatDestroy(&_t->A_G);
    if (_t->A_Pi != NULL)
        phgMatDestroy(&_t->A_Pi);

    return ;
}

void
phgSolverAMSSetConstantCoefficients(SOLVER *solver, FLOAT a, FLOAT b)
{
    GRID *g = solver->mat->rmap->dofs[0]->g;

    assert(solver->oem_solver == SOLVER_AMS);

    if (_t->assembled){
        phgWarning("phgSolverAMSSetConstantCoefficients() should be \
                    called before phgSolverAssemble() %d\n", _t->assembled);
        return ;
    }

    if (_t->alpha != NULL)
	phgDofFree(&_t->alpha);
    if (_t->beta != NULL)
	phgDofFree(&_t->beta);

    _t->alpha = phgDofNew(g, DOF_CONSTANT, 1, "alpha", DofNoAction);
    phgDofSetDataByValue(_t->alpha, a);

    if (b != 1.0) {
	_t->beta = phgDofNew(g, DOF_CONSTANT, 1, "beta", DofNoAction);
	phgDofSetDataByValue(_t->beta, b);
    }

    if (_t->A_G != NULL)
        phgMatDestroy(&_t->A_G);
    if (_t->A_Pi != NULL)
        phgMatDestroy(&_t->A_Pi);

    return;
}

void
phgSolverAMSSetPoisson(SOLVER *solver, MAT *A_G, MAT *A_Pi)
{
   if (_t->assembled){
        phgWarning("phgSolverAMSSetPoisson() should be \
                        called before phgSolverAssemble()!\n");
        return ;
   }

   assert(A_Pi != NULL);

   if (_t->A_G != NULL)
        phgMatDestroy(&_t->A_G);
   if (_t->A_Pi != NULL)
       phgMatDestroy(&_t->A_Pi);

   A_G->refcount++;
   _t->A_G = A_G;
   A_Pi->refcount++;
   _t->A_Pi = A_Pi;

   if (_t->alpha != NULL)
        phgDofFree(&_t->alpha);
   if (_t->beta != NULL)
        phgDofFree(&_t->beta);

   return;
}

static DOF *dof_G, *dof_Pi;	/* scratch variables */

static void
func_G(DOF *u, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* dummy function for computing G.
 * Note: the results are stored as values[u->dim][Dim], here u->dim == nbas */
{
    int i;
    const FLOAT *p, (*J)[Dim + 1];
    DOF_TYPE *type = dof_G->type;

    Unused(bno);
    CheckD_order(u);
    if (type == NULL)
	type = dof_G->hp->info->types[dof_G->hp->elem_order[e->index]];
    /* TODO: cache grad_lambda of the bases of dof_G */
    p = type->BasGrads(dof_G, e, 0, -1, lambda);
    J = (void *)phgGeomGetJacobian(u->g, e);
    for (i = 0; i < u->dim; i++, values += Dim, p += Dim + 1) {
	values[0] = p[0]*J[0][0] + p[1]*J[1][0] + p[2]*J[2][0] + p[3]*J[3][0];
	values[1] = p[0]*J[0][1] + p[1]*J[1][1] + p[2]*J[2][1] + p[3]*J[3][1];
	values[2] = p[0]*J[0][2] + p[1]*J[1][2] + p[2]*J[2][2] + p[3]*J[3][2];
    }
}

static void
func_Pi(DOF *u, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* dummy function for computing Pi.
 * Note: the results are stored as values[u->dim][Dim] */
{
    int i;
    const FLOAT *p;
    DOF_TYPE *type = dof_Pi->type;

    Unused(bno);
    CheckD_order(u);
    if (type == NULL)
	type = dof_Pi->hp->info->types[dof_G->hp->elem_order[e->index]];
    /* TODO: cache bases of dof_Pi */
    p = type->BasFuncs(dof_Pi, e, 0, -1, lambda);
    for (i = 0; i < u->dim / Dim; i++, values += 3 * Dim, p++) {
	values[0] = *p;
	values[1] = 0.;
	values[2] = 0.;

	values[3] = 0.;
	values[4] = *p;
	values[5] = 0.;

	values[6] = 0.;
	values[7] = 0.;
	values[8] = *p;
    }
}

static void
build_G_Pi(SOLVER *solver)
/*----------------------------------------------------------------------------
 * Builds the transfer matrices G (\grad H1->Hcurl) and Pi ((H1)^3->Hcurl).
 *
 * We denote by {\phi} and {\psi}, respectively, the finite element bases of
 * H1 space and Hcurl space.
 *
 *   1. The elements of the gradient matrix G are determined by representing
 *	\grad\phi with \psi:
 *		\grad\phi_j = \sum_k g_{k,j} \psi_k
 *
 *   2. The elements of the transfer matrix Pi_i (i=1,2,3) are determined by
 *	representing e_i\phi with \psi:
 *		e_i\phi_j = \sum_k {Pi_i}_{k,j} \psi_k
 *	where e_1=(1,0,0), e_2=(0,1,0) and e_3=(0,0,1).
 *--------------------------------------------------------------------------*/
{
    MAP *map = solver->mat->rmap, *map_G, *map_Pi;
    DOF *u;
    DOF_TYPE *type;
    GRID *g;
    ELEMENT *e;
    DOF *dof_Hcurl;
    INT *cols, *rows;
    FLOAT *p, *buffer, *buffer0, *buffer1, *buffer2, *pdata[15];
    int i, j, N_Hcurl = 0, N_G = 0, N_Pi = 0;

    assert(map->ndof == 1 && map->dofs != NULL);
    u = map->dofs[0];
    g = u->g;
    type = u->type;
    assert(DofFESpace(u) == FE_Hcurl && u->dim == 1);

    if (!DofIsHP(u)) {
	int order_G, order_Pi;

	order_Pi = u->type->order;
	if (_t->full_order_G) {
	    /* Note: increase order by 1 except for HC0 */
	    order_G = order_Pi + (u->type->np_edge == 1 ? 0 : 1);
	}
	else {
	    /* Note: faster for HC1, though more iterations.
	     * FIXME: can't use GMRES as the aux solver in this case (/0!) */
	    order_G = order_Pi + (order_Pi == 1 ? 0 : 1);
	}

	dof_Pi = phgDofNew(g, DOF_Pn[order_Pi], 1, "dof_Pi", DofNoData);
	dof_G = (order_G == order_Pi ? dof_Pi :
			phgDofNew(g, DOF_Pn[order_G], 1, "dof_G", DofNoData));
	dof_Hcurl = phgDofNew(g, u->type, 1, "dof_Hcurl", DofNoData);

	N_Pi = dof_Pi->type->nbas;
	N_G = dof_G->type->nbas;
	N_Hcurl = u->type->nbas;
    }
    else {
	HP_TYPE *hp;

	hp = phgHPNew(g, HP_HB);
	ForAllElements(g, e)
	    e->hp_order = u->hp->elem_order[e->index];
	phgHPSetup(hp, FALSE);
	dof_Pi = phgHPDofNew(g, hp, 1, "u_h", DofNoData);
	phgHPFree(&hp);

	hp = phgHPNew(g, HP_HB);
	ForAllElements(g, e)
	    e->hp_order = u->hp->elem_order[e->index] + 1;
	phgHPSetup(hp, FALSE);
	dof_G = phgHPDofNew(g, hp, 1, "u_h", DofNoData);
	phgHPFree(&hp);

	dof_Hcurl = phgHPDofNew(g, u->hp, 1, "dof_Hcurl", DofNoData);

	N_Pi = dof_Pi->hp->info->types[dof_Pi->hp->max_order]->nbas;
	N_G = dof_G->hp->info->types[dof_G->hp->max_order]->nbas;
	N_Hcurl = u->hp->info->types[u->hp->max_order]->nbas;
    }

    _t->dof_G = dof_G;
    _t->dof_Pi = dof_Pi;

    map_Pi = phgMapCreate(dof_Pi, NULL);
    map_G = (dof_G == dof_Pi ? map_Pi : phgMapCreate(dof_G, NULL));

    _t->G = phgMapCreateMat(map, map_G);
    phgMatSetMode(_t->G, PHG_REPLACE);

    _t->Pi[0] = phgMapCreateMat(map, map_Pi);
    phgMatSetMode(_t->Pi[0], PHG_REPLACE);
    _t->Pi[1] = phgMapCreateMat(map, map_Pi);
    phgMatSetMode(_t->Pi[1], PHG_REPLACE);
    _t->Pi[2] = phgMapCreateMat(map, map_Pi);
    phgMatSetMode(_t->Pi[2], PHG_REPLACE);

    assert(N_G >= N_Pi);
    buffer = phgAlloc(2 * Dim * N_Hcurl * N_G * sizeof(*buffer));
    buffer0 = buffer  + Dim * N_Hcurl * N_Pi;
    buffer1 = buffer0 + N_Hcurl * N_Pi;
    buffer2 = buffer1 + N_Hcurl * N_Pi;

    rows = phgAlloc((N_Hcurl + N_G) * sizeof(*rows));
    cols = rows  + N_Hcurl;

    ForAllElements(g, e) {

	if (DofIsHP(u)) {
	    N_G = DofNBas(dof_G, e);
	    N_Pi = DofNBas(dof_Pi, e);
	    N_Hcurl = DofNBas(u, e);
	}

	for (j = 0; j < N_Hcurl; j++)
	    rows[j] = phgMapE2L_(map, 0, e, j, FALSE, FALSE);
	for (j = 0; j < N_G; j++)
	    cols[j] = phgMapE2L_(map_G, 0, e, j, TRUE, FALSE);

	/* The gradient matrix */
	dof_Hcurl->dim = N_G;		/* this is hacky! */
	p = buffer;
	if (!DofIsHP(u))
	    j = dof_Hcurl->dim * type->np_edge;
	for (i = 0; i < NEdge; i++) {
	    if (DofIsHP(u)) {
		type = u->hp->info->types[u->hp->edge_order[e->edges[i]]];
		j = dof_Hcurl->dim * type->np_edge;
	    }
	    type->InitFunc(dof_Hcurl, e, EDGE, i, NULL, func_G, NULL, p, pdata);
	    pdata[4 + i] = p;
	    p += j;
	}

	if (!DofIsHP(u))
	    j = dof_Hcurl->dim * type->np_face;
	for (i = 0; i < NFace; i++) {
	    if (DofIsHP(u)) {
		type = u->hp->info->types[u->hp->face_order[e->faces[i]]];
		j = dof_Hcurl->dim * type->np_face;
	    }
	    if (j == 0)
		continue;
	    type->InitFunc(dof_Hcurl, e, FACE, i, NULL, func_G, NULL, p, pdata);
	    pdata[10 + i] = p;
	    p += j;
	}

	if (DofIsHP(u))
	    type = u->hp->info->types[u->hp->elem_order[e->index]];
	if (type->np_elem > 0)
	    type->InitFunc(dof_Hcurl, e, VOLUME, 0, NULL, func_G,
				NULL, p, pdata);

	phgMatSetEntries(_t->G, N_Hcurl, rows, N_G, cols, buffer);

	/* The (H1)^3 matrix */
	dof_Hcurl->dim = Dim * N_Pi;	/* this is hacky! */
	if (dof_Pi != dof_G) {
	    for (j = 0; j < N_Pi; j++)
		cols[j] = phgMapE2L_(map_Pi, 0, e, j, TRUE, FALSE);
	}

	p = buffer;
	if (!DofIsHP(u))
	    j = dof_Hcurl->dim * type->np_edge;
	for (i = 0; i < NEdge; i++) {
	    if (DofIsHP(u)) {
		type = u->hp->info->types[u->hp->edge_order[e->edges[i]]];
		j = dof_Hcurl->dim * type->np_edge;
	    }
	    type->InitFunc(dof_Hcurl, e, EDGE, i, NULL, func_Pi,
				NULL, p, pdata);
	    pdata[4 + i] = p;
	    p += j;
	}

	if (!DofIsHP(u))
	    j = dof_Hcurl->dim * type->np_face;
	for (i = 0; i < NFace; i++) {
	    if (DofIsHP(u)) {
		type = u->hp->info->types[u->hp->face_order[e->faces[i]]];
		j = dof_Hcurl->dim * type->np_face;
	    }
	    if (j == 0)
		continue;
	    type->InitFunc(dof_Hcurl, e, FACE, i, NULL, func_Pi,
				NULL, p, pdata);
	    pdata[10 + i] = p;
	    p += j;
	}

	if (DofIsHP(u))
	    type = u->hp->info->types[u->hp->elem_order[e->index]];
	if (type->np_elem > 0)
	    type->InitFunc(dof_Hcurl, e, VOLUME, 0, NULL, func_Pi,
				NULL, p, pdata);

	dof_Hcurl->dim = 1;		/* restore dof_Hcurl->dim */

	p = buffer;
	for (i = 0; i < N_Hcurl; i++) {
	    for (j = 0; j < N_Pi; j++) {
		buffer0[i * N_Pi + j] = *(p++);
		buffer1[i * N_Pi + j] = *(p++);
		buffer2[i * N_Pi + j] = *(p++);
	    }
	}

	phgMatSetEntries(_t->Pi[0], N_Hcurl, rows, N_Pi, cols, buffer0);
	phgMatSetEntries(_t->Pi[1], N_Hcurl, rows, N_Pi, cols, buffer1);
	phgMatSetEntries(_t->Pi[2], N_Hcurl, rows, N_Pi, cols, buffer2);
    }	/* ForAllElements */

    if (map_G != map_Pi)
	phgMapDestroy(&map_G);
    else
	map_G = NULL;
    phgMapDestroy(&map_Pi);

    phgMatAssemble(_t->G);
    phgMatAssemble(_t->Pi[0]);
    phgMatAssemble(_t->Pi[1]);
    phgMatAssemble(_t->Pi[2]);

    phgFree(rows);
    phgFree(buffer);

    phgDofFree(&dof_Hcurl);
 
    return;
}

static void
build_poisson(SOLVER *solver)
/* builds the Poisson matrices:
 *	_t->A_G <==> -div(beta grad)
 *	_t->A_Pi <==> -div(alpha grad) + beta
 * alpha==NULL means alpha=1, beta==NULL means beta=1 */
{
    DOF *u, *alpha, *beta;
    int i, j, n, N_G, N_Pi;
    GRID *g;
    ELEMENT *e;
    FLOAT *A, *A1, *A0, d;
    INT *I, *J;
    BTYPE *b, DB_mask;

    if (_t->A_Pi != NULL)
	return;

    u = solver->mat->rmap->dofs[0];
    g = u->g;
    alpha = _t->alpha;
    beta = _t->beta;

    /* TODO: check DOF_TYPE of solver->dofs[] for consistency */
    assert(beta == NULL || beta->dim == 1);

    DB_mask = solver->mat->rmap->dofs[0]->DB_mask;

    _t->A_Pi = phgMapCreateMat(_t->Pi[0]->cmap, _t->Pi[0]->cmap);

    if (beta != NULL && beta->type == DOF_CONSTANT && *beta->data == 0.)
	_t->A_G = NULL;
    else
	_t->A_G = phgMapCreateMat(_t->G->cmap, _t->G->cmap);

    /* number of basis functions in an element */
    if (!DofIsHP(u)) {
	N_G = dof_G->type->nbas;
	N_Pi = dof_Pi->type->nbas;
    }
    else {
	N_G = dof_G->hp->info->types[dof_G->hp->max_order]->nbas;
	N_Pi = dof_Pi->hp->info->types[dof_Pi->hp->max_order]->nbas;
    }

    A0 = phgAlloc(N_Pi * N_Pi * sizeof(*A0));

    assert(N_G >= N_Pi);
    A = phgAlloc(N_G * N_G * sizeof(*A));
    A1 = phgAlloc(N_G * sizeof(*A1));
    I = phgAlloc(N_G * sizeof(*I));
    J = phgAlloc(N_G * sizeof(*J));
    b = phgAlloc(N_G * sizeof(*b));

    /* Note: if beta is constant then the element mass matrices only vary
       by element volumes */
    if (g->roots != NULL && !DofIsHP(u) &&
	(beta == NULL || beta->type == DOF_CONSTANT)) {
	e = g->roots;
	d = (beta == NULL ? 1.0 : *beta->data) / phgGeomGetVolume(g, e);
	for (i = 0; i < N_Pi; i++)
	    for (j = 0; j <= i; j++)
		A0[i * N_Pi + j] = A0[j * N_Pi + i] =
		    d * phgQuadBasDotBas(e, dof_Pi, j, dof_Pi, i, -1);
    }

    ForAllElements(g, e) {

	if (DofIsHP(u)) {
	    N_G = DofNBas(dof_G, e);
	    N_Pi = DofNBas(dof_Pi, e);
	}

	/* Matrix _t->A_Pi */
	if (DofIsHP(u) || (beta != NULL && beta->type != DOF_CONSTANT)) {
	    /* compute \int beta * phi_j * phi_i */
	    for (i = 0; i < N_Pi; i++)
		for (j = 0; j <= i; j++)
		    A[i * N_Pi + j] =
			phgQuadBasABas(e, dof_Pi, j, beta, dof_Pi, i, -1);
	}
	else {
	    d = phgGeomGetVolume(g, e);
	    for (i = 0; i < N_Pi; i++)
		for (j = 0; j <= i; j++)
		    A[i * N_Pi + j] = A0[i * N_Pi + j] * d;
	}

	for (i = 0; i < N_Pi; i++) {
	    I[i] = phgMapE2L((_t->A_Pi)->rmap, 0, e, i);
	    b[i] = phgDofGetElementBoundaryType(dof_Pi, e, i * dof_Pi->dim);
	    for (j = 0; j <= i; j++) {
		A[j * N_Pi + i] = (A[i * N_Pi + j] +=
		    phgQuadGradBasAGradBas(e, dof_Pi, j, alpha, dof_Pi, i, -1));
	    }
	}

	for (i = 0; i < N_Pi; i++) {
	    if (b[i] & DB_mask) {
#if 0
		phgMatAddEntry(_t->A_Pi, I[i], I[i], 1.0);
#else
		phgMatAddEntry(_t->A_Pi, I[i], I[i], A[i * N_Pi + i]);
#endif
	    }
	    else {
		for (j = 0, n = 0; j < N_Pi; j++) {
		    if (b[j] & DB_mask)
			continue;
		    J[n] = I[j];
		    A1[n] = A[i * N_Pi + j];
		    n++;
		}
		phgMatAddEntries(_t->A_Pi, 1, I + i, n, J, A1);
	    }
	}

	/* Matrix _t->A_G */
	if (_t->A_G == NULL)
	    continue;

	for (i = 0; i < N_G; i++) {
	    I[i] = phgMapE2L((_t->A_G)->rmap, 0, e, i);
	    b[i] = phgDofGetElementBoundaryType(dof_G, e, i * dof_G->dim);
	    for (j = 0; j <= i; j++) {
		A[j * N_G + i] = A[i * N_G + j] =
		    phgQuadGradBasAGradBas(e, dof_G, j, beta, dof_G, i, -1);
	    }
	}

	for (i = 0; i < N_G; i++) {
	    if (b[i] & DB_mask) {
#if 0
		phgMatAddEntry(_t->A_G, I[i], I[i], 1.0);
#else
		phgMatAddEntry(_t->A_G, I[i], I[i], A[i * N_G + i]);
#endif
	    }
	    else {
		for (j = 0, n = 0; j < N_G; j++) {
		    if (b[j] & DB_mask)
			continue;
		    J[n] = I[j];
		    A1[n] = A[i * N_G + j];
		    n++;
		}
		phgMatAddEntries(_t->A_G, 1, I + i, n, J, A1);
	    }
	}
    }	/* end-ForAllElement */

    phgMatAssemble(_t->A_Pi);
    if (_t->A_G != NULL)
	phgMatAssemble(_t->A_G);

    phgFree(A);
    phgFree(A0);
    phgFree(A1);
    phgFree(I);
    phgFree(J);
    phgFree(b);
}

/*----------------------------------------------------------------------*/

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nThe AMS solver options:", "\n", "ams");

    phgOptionsRegisterString("ams_cycle_type", "AMS cycle type", &cycle_type);

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

    phgOptionsRegisterKeyword("-ams_aux_solver", "Auxiliary solver",
			      phgSolverNames, &aux_solver_id);
    phgOptionsRegisterString("-ams_aux_solver_opts", "Auxiliary solver options",
			      &aux_solver_opts);

    /* options for the Gradient solver */
    phgOptionsRegisterKeyword("-ams_grad_solver", "Gradient solver",
			      phgSolverNames, &grad_solver_id);
    phgOptionsRegisterString("-ams_grad_solver_opts", "Gradient solver options",
			      &grad_solver_opts);

    /* other options */
    phgOptionsRegisterNoArg("-ams_warn_aux", "Warn if auxspaces are undefined",
			    &warn_aux);
    phgOptionsRegisterFloat("-ams_default_alpha", "Default value for alpha",
			    &default_alpha);
    phgOptionsRegisterFloat("-ams_default_beta", "Default value for beta",
			    &default_beta);
    phgOptionsRegisterNoArg("-ams_full_order_G", "Use full order for G",
			    &full_order_G);

    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));

    /* copy relavant cmdline arguments */
    _t->cycle_type = strdup(cycle_type);

    _t->aux_solver_id = aux_solver_id;
    if (aux_solver_opts != NULL)
	_t->aux_solver_opts = strdup(aux_solver_opts);

    _t->alpha = NULL;
    _t->beta = NULL;
    _t->A_G = NULL;
    _t->A_Pi = NULL;
    _t->assembled = FALSE;

    _t->warn_aux = warn_aux;
    _t->full_order_G = full_order_G;
    _t->default_alpha = default_alpha;
    _t->default_beta = default_beta;

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

    if (_t->G != NULL)
	phgMatDestroy(&_t->G);
    if (_t->Pi[0] != NULL)
	phgMatDestroy(&_t->Pi[0]);
    if (_t->Pi[1] != NULL)
	phgMatDestroy(&_t->Pi[1]);
    if (_t->Pi[2] != NULL)
	phgMatDestroy(&_t->Pi[2]);
    if (_t->A_G != NULL)
	phgMatDestroy(&_t->A_G);
    if (_t->A_Pi != NULL)
	phgMatDestroy(&_t->A_Pi);

#if HYPRE_TEST
    if (_t->solver_hypre != NULL)
	phgSolverDestroy(&_t->solver_hypre);
#endif	/* HYPRE_TEST */
    if (_t->solver_G != NULL)
	phgSolverDestroy(&_t->solver_G);
    if (_t->solver_Pi != NULL)
	phgSolverDestroy(&_t->solver_Pi);

    if (_t->dof_G != _t->dof_Pi && _t->dof_G != NULL)
	phgDofFree(&_t->dof_G);
    if (_t->dof_Pi != NULL)
	phgDofFree(&dof_Pi);

    phgFree(_t->cycle_type);
    phgFree(_t->aux_solver_opts);

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    int i;
    MAT_ROW *row;

    if (_t->assembled)
	return 0;

    if (_t->alpha == NULL && _t->A_Pi == NULL) {
	/* use default coefficients alpha=default_alpha and beta=default_beta */
	if (_t->warn_aux)
	    phgPrintf("AMS: auxspaces undefined, "
			"using default values (%lg,%lg).\n",
			(double)_t->default_alpha, (double)_t->default_beta);
	phgSolverAMSSetConstantCoefficients(solver,
				_t->default_alpha, _t->default_beta);
    }

    _t->assembled = TRUE;

    /* preconditioners */
    build_G_Pi(solver);
    build_poisson(solver);

    /* build aux solvers (FIXME: GMRES produces FPE (/0) with HC0 and HC1) */
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
    if (_t->A_G != NULL) {
	_t->solver_G = phgMat2Solver(SOLVER_DEFAULT, _t->A_G);
	_t->solver_G->warn_maxit = FALSE;
	phgVecDestroy(&_t->solver_G->rhs);
	_t->solver_G->rhs_updated = TRUE;
    }
    _t->solver_Pi = phgMat2Solver(SOLVER_DEFAULT, _t->A_Pi);
    _t->solver_Pi->warn_maxit = FALSE;
    phgVecDestroy(&_t->solver_Pi->rhs);
    _t->solver_Pi->rhs_updated = TRUE;
    phgOptionsPop();

#if HYPRE_TEST
    /* create a HYPRE AMS preconditioner for testing */
    phgOptionsPush();
    phgSolverSetDefaultSuboptions();
    phgOptionsSetOptions("-hypre_solver ams -hypre_pc none");
    _t->solver_hypre = phgMat2Solver(SOLVER_HYPRE, solver->mat);
    phgVecDestroy(&_t->solver_hypre->rhs);
    _t->solver_hypre->rhs_updated = TRUE;
    phgSolverHypreAMSSetConstantPoisson(_t->solver_hypre, 1.0, 1.0);
    phgOptionsPop();
#endif	/* HYPRE_TEST */

    /* penalize boundary rows of A_G */
    for (i = 0; _t->A_G != NULL && i < _t->A_G->rmap->nlocal; i++) {
	row = _t->A_G->rows + i;
	if (row->ncols == 0) {
	    /* empty row */
	    row->ncols = row->alloc = 1;
	    row->cols = phgAlloc(sizeof(*row->cols));
	    row->cols[0] = i;
	    row->data = phgAlloc(sizeof(*row->data));
	    row->data[0] = DBL_MAX;
	}
	if (row->ncols == 1 && row->cols[0] == i)
	    row->data[0] = DBL_MAX;
    }

    /* penalize boundary rows of A_Pi */
    for (i = 0; i < _t->A_Pi->rmap->nlocal; i++) {
	row = _t->A_Pi->rows + i;
	assert(row->ncols > 0);
	if (row->ncols == 1 && row->cols[0] == i)
	    row->data[0] = DBL_MAX;
    }

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
	    res = resp[0];
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
    int i;
    VEC *r0 = NULL, *r1 = NULL, *x1 = NULL;

    while (*cycle != '\0') {
	if (*cycle < '0' || *cycle > '2')
	    phgError(1, "invalid AMS cycle type string: %s\n", _t->cycle_type);
	switch (*(cycle++) - '0') {
	    case 0:	/* smoothing */
		if (phgVerbosity > 1)
		    phgPrintf("\tAMS smoothing\n");
		gauss_seidel(solver->mat, r, &x, 1, 0., NULL, phgVerbosity - 1);
		break;
	    case 1:	/* space G correction */
		if (_t->solver_G == NULL)
		    break;
		if (phgVerbosity > 1)
		    phgPrintf("\tAMS space G correction\n");
		phgVecCopy(r, &r0);
		phgMatVec(MAT_OP_N, -1.0, solver->mat, x, 1.0, &r0);
		if (x1 == NULL || x1->map != _t->G->cmap) {
		    phgVecDestroy(&x1);
		    x1 = phgMapCreateVec(_t->G->cmap, 1);
		}
		else {
		    bzero(x1->data, x1->map->nlocal * sizeof(*x1->data));
		}
		if (r1 != NULL && r1->map != _t->G->cmap)
		    phgVecDestroy(&r1);
		phgMatVec(MAT_OP_T, 1.0, _t->G, r0, 0.0, &r1);
		_t->solver_G->rhs = r1;
		_t->solver_G->rhs->assembled = TRUE;
		phgSolverVecSolve(_t->solver_G, FALSE, x1);
		_t->solver_G->rhs = NULL;
		phgMatVec(MAT_OP_N, 1.0, _t->G, x1, 1.0, &x);
		break;
	    case 2:	/* space Pi correction */
		if (phgVerbosity > 1)
		    phgPrintf("\tAMS space Pi correction\n");
		phgVecCopy(r, &r0);
		phgMatVec(MAT_OP_N, -1.0, solver->mat, x, 1.0, &r0);
		if (x1 == NULL || x1->map != _t->Pi[0]->cmap) {
		    phgVecDestroy(&x1);
		    x1 = phgMapCreateVec(_t->Pi[0]->cmap, 1);
		}
		if (r1 != NULL && r1->map != _t->Pi[0]->cmap)
		    phgVecDestroy(&r1);
		for (i = 0; i < Dim; i++) {
#if 0
		    if (i > 0) {
			/* multiplicative corrections of the three components */
			phgVecCopy(r, &r0);
			phgMatVec(MAT_OP_N, -1.0, solver->mat, x, 1.0, &r0);
		    }
#endif
		    phgMatVec(MAT_OP_T, 1.0, _t->Pi[i], r0, 0.0, &r1);
		    _t->solver_Pi->rhs = r1;
		    _t->solver_Pi->rhs->assembled = TRUE;
		    bzero(x1->data, x1->map->nlocal * sizeof(*x1->data));
		    phgSolverVecSolve(_t->solver_Pi, FALSE, x1);
		    _t->solver_Pi->rhs = NULL;
		    phgMatVec(MAT_OP_N, 1.0, _t->Pi[i], x1, 1.0, &x);
		}
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
	DOF_TYPE *type_G = _t->A_G->rmap->dofs[0]->type,
		 *type_Pi = _t->A_Pi->rmap->dofs[0]->type;
        phgPrintf("======================== AMS =====================\n");
	phgPrintf("Auxiliary solver: %s\n", _t->solver_Pi->oem_solver->name);
	phgPrintf("AMS cycle type: %s, Space G: %s, Space Pi: %s\n",
			_t->cycle_type,
			type_G == NULL ? "hp" : type_G->name,
			type_Pi == NULL ? "hp" : type_Pi->name);
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
	    (double)solver->atol, (double)solver->rtol, (double)solver->btol);
        phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
        phgPrintf("-----   ------------   -------------  ------------\n");
    }
    while (solver->nits < solver->maxit) {
	solver->nits++;
#if HYPRE_TEST
	/* test preconditioning using Hypre AMS */
	_t->solver_hypre->rtol = solver->rtol;
	_t->solver_hypre->maxit = solver->maxit;
	_t->solver_hypre->rhs = solver->rhs;
	phgSolverVecSolve(_t->solver_hypre, FALSE, x);
	_t->solver_hypre->rhs = NULL;
	solver->residual = _t->solver_hypre->residual;
	return solver->nits = _t->solver_hypre->nits;
#endif	/* HYPRE_TEST */

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

OEM_SOLVER phgSolverAMS_ = {
    "AMS", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, TRUE, FALSE
};
