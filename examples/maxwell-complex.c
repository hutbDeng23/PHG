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

/* Testing solving time-harmonic Maxwell equations with complex coefficients.
 *
 *	curl mu curl E - k^2 E = J in Omega,
 *	E x n = g x n at boundary
 *
 *      J= (1,1,1) + i (1,1,1),
 *	g = (1,1,1) + i (1,1,1)
 *
 *	Omega = (-1,1)^3,
 *	Omega_1 = (-.5,.5)x(-.5,.5)x(0,.5),
 *	Omega_2 = Omega - Omega_1.
 *
 * Test case 1:
 *      mu=1, k^2=1-i in Omega_1, k^2=1-2i in Omega_2
 *
 * Test case 2:
 *      mu=1 in Omega_1, mu=2 in Omega_2, k^2=1
 *
 * Test case 3:
 *      mu=1 in Omega_1, mu=1000 in Omega_2, k^2=1
 *
 * Note: This is the second example in Wang Long's  Ph.D thesis.
 *
 * The discrete system is written as:
 *
 *	|  A  -B | | E_r |   | J_r |
 *	|	 | |     | = |	   |
 *	| -B  -A | | E_i |   |-J_i |
 *
 * where A <=> curl curl - Re(k^2), B <=> -Im(k^2), and diag(A+B,A+B) is
 * used to precondition the system (which is solved by PCG-AMS).
 *
 * $Id: maxwell-complex.c,v 1.51 2022/09/21 06:00:10 zlb Exp $
 */

#ifndef USE_BLOCK_MATRIX
# define USE_BLOCK_MATRIX 1	/* using block matrix saves about 25% memory */
#endif

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdlib.h>

static const char *test_cases[] = {"1", "2", "3", NULL};
static int test_case = 0;
static BOOLEAN use_ams = TRUE;

static void
func_mu(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    switch (test_case) {
	case 0:
	    *value = 1.;
	    break;
	case 1:
	    *value = (x > -.5 && x < .5 && y > -.5 && y < .5 && z > 0. && z < .5
			? 1. : 2.);
	    break;
	case 2:
	    *value = (x > -.5 && x < .5 && y > -.5 && y < .5 && z > 0. && z < .5
			? 1. : 1000.);
	    break;
    }
}

static void
func_kk_re(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    switch (test_case) {
	case 0:
	    *value = 1.;
	    break;
	case 1:
	case 2:
	    *value = 1.;
	    break;
    }
}

static void
func_kk_im(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    switch (test_case) {
	case 0:
	    *value = (x > -.5 && x < .5 && y > -.5 && y < .5 && z > 0. && z < .5
			? -1.0 : -2.0);
	    break;
	case 1:
	case 2:
	    *value = 0.;
	    break;
    }
}

static double
elapsed_time(GRID *g, BOOLEAN flag, double mflops)
/* returns and prints elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;
 
    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    mem = phgMemoryUsage(g, NULL);

    if (flag) {
	if (mflops <= 0)
	    phgPrintf("[%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
	else
	    phgPrintf("[%0.4lgMB %0.4lfs %0.4lgGF]\n", mem / (1024.0 * 1024.0),
			 et, mflops*1e-3);
    }

    return et;
}

static void
func_g_re(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    *(values++) = 1.0;
    *(values++) = 1.0;
    *(values++) = 1.0;
    return;
}

static void
func_g_im(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    *(values++) = 1.0;
    *(values++) = 1.0;
    *(values++) = 1.0;
    return;
}

static void
#if USE_BLOCK_MATRIX
build_linear_system(MAT *a, MAT *b, VEC *rhs, SOLVER *pc,
		    DOF *f_re, DOF *f_im, DOF *mu, DOF *kk_re, DOF *kk_im)
#else	/* USE_BLOCK_MATRIX */
build_linear_system(SOLVER *solver, SOLVER *pc,
		    DOF *f_re, DOF *f_im, DOF *mu, DOF *kk_re, DOF *kk_im)
#endif	/* USE_BLOCK_MATRIX */
{
    DOF *u_h = pc->mat->rmap->dofs[0];
    int N = u_h->type->nbas;	/* number of basis functions in an element */
    int M = 2 * N;
    GRID *g = u_h->g;
    ELEMENT *e;
    int i, j;
    FLOAT buffer[N], C[M], A0[N][N], kkre, kkim, stiffness, mass;
    INT I[M];
#if USE_BLOCK_MATRIX
    FLOAT A[N][N], B[N][N];
#else	/* USE_BLOCK_MATRIX */
    VEC *rhs = solver->rhs;
    FLOAT A[M][M];
    assert(solver->mat->rmap->nlocal == 2 * pc->mat->rmap->nlocal);
#endif	/* USE_BLOCK_MATRIX */

    assert(u_h->dim == 1);

    ForAllElements(g, e) {
	kkre = *DofElementData(kk_re, e->index);
	kkim = *DofElementData(kk_im, e->index);
	for (i = 0; i < N; i++) {
	    I[i]	= phgMapE2L(rhs->map, 0, e, i);
	    I[N + i]	= phgMapE2L(rhs->map, 1, e, i);
	    for (j = 0; j <= i; j++) {
		stiffness = phgQuadCurlBasACurlBas(e, u_h, j, mu, u_h, i, -1);
		mass = phgQuadBasDotBas(e, u_h, j, u_h, i, -1);
		/* the preconditioner matrix */
		A0[j][i] = A0[i][j] = stiffness + mass*(Fabs(kkre)+Fabs(kkim));
#if USE_BLOCK_MATRIX
		A[i][j] = A[j][i] = stiffness - mass * kkre;
		B[i][j] = B[j][i] = mass * kkim;
#else	/* USE_BLOCK_MATRIX */
		/* the original matrix, upper (real) part */
		A[j][i] = A[i][j] = stiffness - mass * kkre;
		A[N + j][i] = A[i][N + j] = mass * kkim;
		/* the original matrix, lower (imaginary) part */
		A[N + j][N + i] = A[N + i][N + j] = -stiffness + mass * kkre;
		A[j][N + i] = A[N + i][j] = mass * kkim;
#endif	/* USE_BLOCK_MATRIX */
	    }
	}

	/* compute local matrix and RHS */
	for (i = 0; i < N; i++) {	/* loop on basis functions */
	    if (phgDofDirichletBC(u_h, e, i, func_g_re, 
				buffer, C + i, DOF_PROJ_CROSS)) {
		/* Dirichlet boundary */
		if (use_ams)
		    phgSolverAddMatrixEntries(pc, 1, I + i, N, I, buffer); 
#if USE_BLOCK_MATRIX
		phgMatAddEntries(a, 1, I + i, N, I, buffer);
		phgDofDirichletBC(u_h, e, i, func_g_im,
				buffer, C + N + i, DOF_PROJ_CROSS);
		/* Note: the second diagonal block is multiplied by -1,
		 *	 so the RHS must also be negated */
		C[N + i] = -C[N + i];
#else	/* USE_BLOCK_MATRIX */
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer); 
		phgDofDirichletBC(u_h, e, i, func_g_im,
				buffer, C + N + i, DOF_PROJ_CROSS);
#if 0
		/* (Optionally) negate the coefficients and RHS */
		C[N + i] = -C[N + i];
		for (j = 0; j < N; j++)
		    buffer[j] = -buffer[j];
#endif
		phgSolverAddMatrixEntries(solver, 1, I+N+i, N, I+N, buffer); 
#endif	/* USE_BLOCK_MATRIX */
	    }
	    else {
		/* the preconditioning matrix */
		if (use_ams)
		    phgSolverAddMatrixEntries(pc, 1, I+i, N, I, &A0[i][0]); 
		/* right hand side: \int f * phi_i */
		C[i]	 =  phgQuadDofDotBas(e, f_re, u_h, i, QUAD_DEFAULT);
		C[N + i] = -phgQuadDofDotBas(e, f_im, u_h, i, QUAD_DEFAULT);
		/* the original matrix */
#if USE_BLOCK_MATRIX
		phgMatAddEntries(a, 1, I+i, N, I, &A[i][0]);
		phgMatAddEntries(b, 1, I+i, N, I, &B[i][0]);
#else	/* USE_BLOCK_MATRIX */
		phgSolverAddMatrixEntries(solver, 1, I+i, M, I, &A[i][0]);
		phgSolverAddMatrixEntries(solver, 1, I+N+i, M, I, &A[N+i][0]);
#endif	/* USE_BLOCK_MATRIX */
	    }
	}
	phgVecAddEntries(rhs, 0, M, I, C);
    }
}

static FLOAT
estimate_error(DOF *u_re, DOF *u_im, DOF *f_re, DOF *f_im, DOF *mu,
		DOF *kk_re, DOF *kk_im, DOF *error)
{
    GRID *g = u_re->g;
    ELEMENT *e;
    DOF *jmp1_re, *jmp1_im, *jmp2_re, *jmp2_im, *res_re, *res_im;
    DOF *curl_re, *curl_im, *div_re, *div_im, *tmp;

    tmp = phgDofCurl(u_re, NULL, NULL, "curl_u_re");
    curl_re = phgDofMM(MAT_OP_N,MAT_OP_N, 1,3,1, 1.0, mu, 1, tmp, 0.0, NULL);
    phgDofFree(&tmp);

    tmp = phgDofCurl(u_im, NULL, NULL, "curl_u_im");
    curl_im = phgDofMM(MAT_OP_N,MAT_OP_N, 1,3,1, 1.0, mu, 1, tmp, 0.0, NULL);
    phgDofFree(&tmp);

    res_re = phgDofGetSameOrderDG(u_re, -1, NULL);
    phgDofCopy(f_re, &res_re, NULL, NULL);
    phgDofMM(MAT_OP_N, MAT_OP_N, 1, 3, 1, +1.0, kk_re, 1, u_re, 1.0, &res_re);
    phgDofMM(MAT_OP_N, MAT_OP_N, 1, 3, 1, -1.0, kk_im, 1, u_im, 1.0, &res_re);
    jmp1_re = phgQuadFaceJump(curl_re, DOF_PROJ_CROSS, NULL, QUAD_DEFAULT);
    jmp2_re = phgQuadFaceJump(res_re, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    tmp = phgDofCopy(res_re, NULL, NULL, NULL);
    div_re = phgDofCurl(curl_re, NULL, NULL, "residual");
    phgDofAXPBY(-1.0, div_re, 1.0, &res_re);	/* residual */
    phgDofFree(&div_re);
    div_re = phgDofDivergence(tmp, NULL, NULL, "div(f+k^2 u_h)");
    phgDofFree(&tmp);

    res_im = phgDofGetSameOrderDG(u_im, -1, NULL);
    phgDofCopy(f_im, &res_im, NULL, NULL);
    phgDofMM(MAT_OP_N, MAT_OP_N, 1, 3, 1, +1.0, kk_re, 1, u_im, 1.0, &res_im);
    phgDofMM(MAT_OP_N, MAT_OP_N, 1, 3, 1, +1.0, kk_im, 1, u_re, 1.0, &res_im);
    jmp1_im = phgQuadFaceJump(curl_im, DOF_PROJ_CROSS, NULL, QUAD_DEFAULT);
    jmp2_im = phgQuadFaceJump(res_im, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    tmp = phgDofCopy(res_im, NULL, NULL, NULL);
    div_im = phgDofCurl(curl_im, NULL, NULL, "residual");
    phgDofAXPBY(-1.0, div_im, 1.0, &res_im);	/* residual */
    phgDofFree(&div_im);
    div_im = phgDofDivergence(tmp, NULL, NULL, "div(f+k^2 u_h)");
    phgDofFree(&tmp);

    ForAllElements(g, e) {
	int i;
	FLOAT eta, h;
	FLOAT diam = phgGeomGetDiameter(g, e);
	e->mark = 0;		/* clear refinement mmark */
	eta = 0.0;
	for (i = 0; i < NFace; i++) {
	    INT fno;
	    if (e->bound_type[i] == DIRICHLET || e->bound_type[i] == NEUMANN)
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    fno = e->faces[i];
	    eta += (*DofFaceData(jmp1_re, fno) + *DofFaceData(jmp2_re, fno) +
		    *DofFaceData(jmp1_im, fno) + *DofFaceData(jmp2_im, fno)
		   ) * h;
	}
	eta += (phgQuadDofDotDof(e, res_re, res_re, QUAD_DEFAULT) +
		phgQuadDofDotDof(e, res_im, res_im, QUAD_DEFAULT) +
		phgQuadDofDotDof(e, div_re, div_re, QUAD_DEFAULT) +
		phgQuadDofDotDof(e, div_im, div_im, QUAD_DEFAULT)
	       ) * diam * diam;
	*DofElementData(error, e->index) = eta;
    }
    phgDofFree(&curl_re);
    phgDofFree(&jmp1_re);
    phgDofFree(&jmp2_re);
    phgDofFree(&res_re);
    phgDofFree(&div_re);
    phgDofFree(&curl_im);
    phgDofFree(&jmp1_im);
    phgDofFree(&jmp2_im);
    phgDofFree(&res_im);
    phgDofFree(&div_im);
    return Sqrt(phgDofNormInftyVec(error));
}

static void
pc_proc(void *ctx, VEC *b0, VEC **x0)
{
    SOLVER *pc_solver = ctx;
    FLOAT *old_rhs;
    VEC *b, *x;
    INT n = pc_solver->mat->rmap->nlocal;

    assert(2 * n == b0->map->nlocal);

    phgSolverAssemble(pc_solver);
	
    b = phgMapCreateVec(pc_solver->mat->rmap, 1);
    x = phgMapCreateVec(pc_solver->mat->rmap, 1);

    /* setup RHS */  
    old_rhs = pc_solver->rhs->data;
    pc_solver->rhs->data = b->data;
    pc_solver->rhs->assembled=TRUE;
    pc_solver->rhs_updated=TRUE;
 
    /* real part */
    memcpy(b->data, b0->data, sizeof(*b->data) * n);
    bzero(x->data, n * sizeof(*x->data));
    phgSolverVecSolve(pc_solver, FALSE, x);
    memcpy((*x0)->data, x->data, sizeof(*x->data) * n);
 
    /* imaginary part */
    memcpy(b->data, b0->data + n, sizeof(*b->data) * n);
    bzero(x->data, n * sizeof(*x->data));
    phgSolverVecSolve(pc_solver, FALSE, x);
    memcpy((*x0)->data + n, x->data, sizeof(*x->data) * n);
    
    pc_solver->rhs->data = old_rhs;

    phgVecDestroy(&b);
    phgVecDestroy(&x);

    return;  
}

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *u_re, *u_im, *f_re, *f_im, *error, *mu, *kk_re, *kk_im, *tmp;
    SOLVER *solver, *pc = NULL, *pc_ams = NULL;
#if USE_BLOCK_MATRIX
    MAT *D, *O, *pmat[4];
    FLOAT coef[4];
#endif	/* USE_BLOCK_MATRIX */
    FLOAT a, b, max_err, L2_err;
    size_t mem, mem_peak;

    char *fn = "../test/cube3.dat";
    char *vtk = NULL;
    INT pre_refines = 0;
    INT depth = 1;	/* <=0 => uniformly refine -depth or 3 times */
    INT mem_max = 400;	/* max. memory */
    FLOAT tol = 1e-3;	/* convergence tolerance */
    INT submesh_threshold = 1000;
    FLOAT lif_threshold = 1.2;

    phgOptionsRegisterNoArg("-use_ams", "Use AMS preconditioner", &use_ams);
    phgOptionsRegisterKeyword("-test_case", "Test case", test_cases, &test_case);
    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **) &fn);
    phgOptionsRegisterFilename("-vtk_file", "VTK file", (char **) &vtk);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-refine_depth", "Refinement depth (<=0 means "
			  "uniformly refine -depth or 3 times)", &depth);
    phgOptionsRegisterInt("-mem_max", "Max memory (MB)", &mem_max);
    phgOptionsRegisterFloat("-tol", "Convergence tolerance", &tol);
    phgOptionsRegisterFloat("-lif_threshold", "LIF threshold", &lif_threshold);
    phgOptionsRegisterInt("-submesh_threshold", "Submesh threshold",
							&submesh_threshold);

    phgOptionsPreset("-dof_type ND1");
#if USE_MINRES
    phgOptionsPreset("-solver minres -solver_symmetry sym");
#else	/* USE_MINRES */
    phgOptionsPreset("-solver gmres -gmres_restart 5");
#endif	/* USE_MINRES */

    phgInit(&argc, &argv);
    phgOptionsShowUsed();
    if (depth == 0)
	depth = -3;

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    phgRefineAllElements(g, pre_refines);

    u_re = phgDofNew(g, DOF_DEFAULT, 1, "u_re", DofInterpolation);
    u_re->DB_mask = BDRY_MASK;
    u_im = phgDofNew(g, DOF_DEFAULT, 1, "u_im", DofInterpolation);
    u_im->DB_mask = BDRY_MASK;
    f_re = phgDofNew(g, DOF_CONSTANT, 3, "f_re", DofNoAction);
    phgDofSetDataByValuesV(f_re, (FLOAT)1.0, (FLOAT)1.0, (FLOAT)1.0);
    f_im = phgDofNew(g, DOF_CONSTANT, 3, "f_im", DofNoAction);
    phgDofSetDataByValuesV(f_im, (FLOAT)1.0, (FLOAT)1.0, (FLOAT)1.0);

    mu = phgDofNew(g, DOF_P0, 1, "mu", DofInterpolation);
    phgDofSetDataByFunction(mu, func_mu);

    kk_re = phgDofNew(g, DOF_P0, 1, "kk_re", DofInterpolation);
    phgDofSetDataByFunction(kk_re, func_kk_re);
    kk_im = phgDofNew(g, DOF_P0, 1, "kk_im", DofInterpolation);
    phgDofSetDataByFunction(kk_im, func_kk_im);

    error = phgDofNew(g, DOF_P0, 1, "error", DofNoAction);

    while (TRUE) {
	elapsed_time(g, FALSE, 0.);
	phgPrintf("\n------ %"dFMT" DOF, %"dFMT" elements, mesh LIF = %lg\n",
	    2 * DofDataCountGlobal(u_re), g->nleaf_global, (double)g->lif);

	if (phgBalanceGrid(g, lif_threshold, submesh_threshold, NULL, 0.)) {
	    phgPrintf("------ Repartition mesh: nprocs = %d, LIF = %lg ",
			    g->nprocs, (double)g->lif);
	    elapsed_time(g, TRUE, 0.);
	}

	phgPrintf("Set up linear solver: ");
	solver = phgSolverCreate(SOLVER_DEFAULT, u_re, u_im, NULL);
	/*solver->mat->handle_bdry_eqns = FALSE;*/
	phgPrintf("solver LIF = %lg ", (double)solver->rhs->map->lif);

	if (TRUE/*use_ams*/) {
	    phgOptionsPush();
	    phgOptionsSetOptions("-solver pcg -solver_symmetry spd");
	    pc = phgSolverCreate(SOLVER_DEFAULT, u_re, NULL);
	    phgOptionsPop();
	    pc->maxit = 10;
	    pc->rtol = 1e-2;
	    pc->warn_maxit = FALSE;
	    elapsed_time(g, TRUE, 0.);
	}

#if USE_BLOCK_MATRIX
	/* diagonal block */
	D = phgMapCreateMat(pc->mat->rmap, pc->mat->cmap);
	D->handle_bdry_eqns = solver->mat->handle_bdry_eqns;
	/* off-diagonal block */
	O = phgMapCreateMat(pc->mat->rmap, pc->mat->cmap);
	O->handle_bdry_eqns = FALSE;
	pmat[0] = D;	pmat[1] = O;
	coef[0] = 1.0;	coef[1] = 1.0;
	pmat[2] = O;	pmat[3] = D;
	coef[2] = 1.0;	coef[3] = -1.0;
	phgMatDestroy(&solver->mat);
	solver->mat = phgMatCreateBlockMatrix(g->comm, 2, 2, pmat, coef, NULL);
	solver->rhs->mat = solver->mat;
#endif	/* USE_BLOCK_MATRIX */

	phgPrintf("Build linear system ");
#if USE_BLOCK_MATRIX
	build_linear_system(D, O, solver->rhs, pc, f_re,f_im, mu, kk_re,kk_im);
	phgMatDestroy(&D);
	phgMatDestroy(&O);
#else	/* USE_BLOCK_MATRIX */
	build_linear_system(solver, pc, f_re,f_im, mu, kk_re,kk_im);
#endif	/* USE_BLOCK_MATRIX */
	elapsed_time(g, TRUE, 0.);

	if (use_ams) {
	    phgPrintf("Set up preconditioner: ");
	    tmp = phgDofAFXPBY1(1.0, Fabs, kk_re, 0.0, NULL);
	    phgDofAFXPBY1(1.0, Fabs, kk_im, 1.0, &tmp);
	    pc_ams = phgMat2Solver(SOLVER_AMS, pc->mat);
	    phgSolverAMSSetCoefficients(pc_ams, mu, tmp);
	    phgDofFree(&tmp);
	    pc_ams->rtol = 0.;
	    pc_ams->maxit = 1;
	    phgSolverSetPC(pc, pc_ams, NULL);
	    phgSolverSetPC(solver, pc, pc_proc);
	    elapsed_time(g, TRUE, 0.);
	}

	phgPrintf("Solve linear system: ");
	phgSolverSolve(solver, TRUE, u_re, u_im, NULL);
	phgPrintf("nits=%d, resid=%0.4lg ", solver->nits,
			(double)solver->residual);
	elapsed_time(g, TRUE, 0.);
	if (pc != NULL)
	    phgSolverDestroy(&pc);
	if (pc_ams != NULL)
	    phgSolverDestroy(&pc_ams);
	phgSolverDestroy(&solver);

	phgPrintf("Norms: L2(u_re)=%0.8lg, L2(u_im)=%0.8lg\n",
			(double)phgDofNormL2(u_re), (double)phgDofNormL2(u_im));
	max_err = estimate_error(u_re,u_im, f_re,f_im, mu, kk_re,kk_im, error);
	L2_err = Sqrt(phgDofNormL1Vec(error));
	phgPrintf("Estimate: PE_oo=%0.8le, PE_L2=%0.8le ",
			(double)max_err, (double)L2_err);
	elapsed_time(g, TRUE, 0.);
	mem = phgMemoryUsage(g, &mem_peak);
	a = mem / (1024.0 * 1024.0);
	b = mem_peak / (1024.0 * 1024.0);
	phgPrintf("Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		(double)a, (double)b);
	if (L2_err < tol || mem_peak >= 1024 * (size_t)mem_max * 1024) {
	    if (vtk != NULL) {
		const char *vtk_fn;
		phgDofFree(&f_re);
		phgDofFree(&f_im);
		phgRefineAllElements(g, 0);
		phgPrintf("Export mesh to \"%s\" ",
			vtk_fn = phgExportVTK(g, vtk, u_re, u_im, kk_im, NULL));
                if (phgRank == 0) {
                    FILE *f;
                    long long size;
                    f = fopen(vtk_fn, "r");
                    fseek(f, 0L, SEEK_END);
                    size = ftell(f);
                    fclose(f);
                    phgPrintf("size = %lld ", size);
                }
		elapsed_time(g, TRUE, 0.);
	    }
	    break;
	}

	phgPrintf("Refine mesh ");
	if (depth <= 0) {
	    phgRefineAllElements(g, -depth);
	}
	else {
	    phgMarkRefine(MARK_DEFAULT, error, Pow(0.7, 2), NULL, 0.0,
				depth, Pow(tol, 2) / g->nleaf_global);
	    phgRefineMarkedElements(g);
	}
	elapsed_time(g, TRUE, 0.);
    }
    phgDofFree(&u_re);
    phgDofFree(&u_im);
    phgDofFree(&f_re);
    phgDofFree(&f_im);
    phgDofFree(&error);
    phgDofFree(&kk_re);
    phgDofFree(&kk_im);
    phgDofFree(&mu);

    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
