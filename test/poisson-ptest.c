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

/* $Id: poisson-ptest.c,v 1.35 2020/12/06 09:03:08 zlb Exp $ */

#include "phg.h"

#include <string.h>
#include <math.h>

/*
    This test code solves the Poisson equation with Dirichlet BC:

	- \Delta u = f
	u |_{\partial\Omega} = g
*/

/* functions for analytical solution, BC and RHS */
#ifndef TEST_CASE
# define TEST_CASE 3
#endif
#include "../examples/functions.h"

/* constants for evaluating L2 error indicator */
#define C0L2	1.0
#define C1L2	0.5
/* constants for evaluating H1 error indicator */
#define C0H1	1.0
#define C1H1	0.5

static INT dim_u = 1;

#define FUNC_DEF(func_f, f_)				\
static void						\
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)		\
{							\
    int i;						\
    FLOAT d;						\
							\
    if (dim_u == 1) {					\
	*value = f_(x, y, z);				\
	return;						\
    }							\
							\
    d = f_(x, y, z);					\
    for (i = 0; i < dim_u; i++)				\
	*(value++) = d;					\
}

FUNC_DEF(func_f, f_)
FUNC_DEF(func_g, g_)
FUNC_DEF(func_u, u_)

static void
func_grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    int i;

    for (i = 0; i < dim_u; i++, value += Dim)
	grad_u_(value, x, y, z);
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
	    phgPrintf("[%0.4lgMB %0.4lfs]\n",
			(double)mem / (1024.0 * 1024.0), et);
	else
	    phgPrintf("[%0.4lgMB %0.4lfs %0.4lgGF]\n",
			(double)mem / (1024.0 * 1024.0), et, mflops*1e-3);
    }

    return et;
}

static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
    int N = u_h->type->nbas;    /* number of basis functions in an element */
    int dim = u_h->dim;         /* DOF dimension */
    GRID *g = u_h->g;
#if USE_OMP
#pragma omp parallel
  {
#endif	/* USE_OMP */
    ELEMENT *e;
    int i, j, k;
    FLOAT (*A)[N], (*B)[dim], *buffer, *row;
    INT (*I)[dim], (*J)[N];

    A = phgAlloc(N * sizeof(*A));
    B = phgAlloc(N * sizeof(*B));
    buffer = phgAlloc(N * sizeof(*buffer));
    I = phgAlloc(N * sizeof(*I));
    J = phgAlloc(dim * sizeof(*J));
#if USE_OMP
#pragma omp for
#endif	/* USE_OMP */
    ForAllElementsBegin(g, e)
	/* compute \int \grad\phi_j \cdot \grad\phi_i making use of symmetry */
	for (i = 0; i < N; i++) {
	    for (k = 0; k < dim; k++)
		I[i][k] = J[k][i] = phgSolverMapE2L(solver, 0, e, i * dim + k);
	    for (j = 0; j <= i; j++)
		A[j][i] = A[i][j] =
		    phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	}

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, func_g,
				  buffer, B[i], DOF_PROJ_NONE)) {
		row = buffer;
	    }
	    else if (phgDofGetElementBoundaryType(u_h, e, i * dim) & NEUMANN) {
		/* TODO */
		row = NULL;
		phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__);
	    }
	    else {
		/* right hand side: \int f * phi_i */
		phgQuadDofTimesBas(e, f_h, u_h, i, QUAD_DEFAULT, B[i]);
		row = A[i];
	    }
	    for (k = 0; k < dim; k++)
		phgSolverAddMatrixEntries(solver, 1, I[i] + k, N, J[k], row); 
	}
	phgSolverAddRHSEntries(solver, N * dim, I[0], B[0]);
    ForAllElementsEnd
    phgFree(A);
    phgFree(B);
    phgFree(buffer);
    phgFree(I);
    phgFree(J);
#if USE_OMP
  }
#endif	/* USE_OMP */
}

static void
estimate_error(DOF *u_h, DOF *f_h, DOF *grad_u, DOF *e_L2, DOF *e_H1)
/* compute L2 and H1 error indicators */
{
    int i;
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *jump, *residual;
    FLOAT eta0, eta1, d, h;
    FLOAT diam;

    jump = phgQuadFaceJump(grad_u, DOF_PROJ_DOT, "jumps", QUAD_DEFAULT);
    residual = phgDofDivergence(grad_u, NULL, NULL, "div");
    phgDofAXPY(1.0, f_h, &residual);

    ForAllElements(g, e) {
	diam = phgGeomGetDiameter(g, e);
	e->mark = 0;		/* clear refinement mmark */
	eta0 = eta1 = 0.0;
	/* for each face F compute [grad_u \cdot n] */
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] & (DIRICHLET | NEUMANN))
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    d = *DofFaceData(jump, e->faces[i]);
	    eta0 += d * h * h * h;
	    eta1 += d * h;
	}
	d = phgQuadDofDotDof(e, residual, residual, -1);
	eta0 = eta0 * C0L2 + diam * diam * diam * diam * d * C1L2;
	eta1 = eta1 * C0H1 + diam * diam * d * C1H1;

	/* add curved boundary errors (FIXME: how to normalize?) */
	d = phgGetBoundaryError(g, e);
	eta0 += d;
	eta1 += d;

	*DofElementData(e_L2, e->index) = Sqrt(eta0);
	*DofElementData(e_H1, e->index) = Sqrt(eta1);
    }
    phgDofFree(&jump);
    phgDofFree(&residual);

    return;
}

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *u_h, *f_h, *grad_u = NULL, *e_L2, *e_H1;
    SOLVER *solver;
    FLOAT L2_error = 0.0, H1_error = 0.0, EL2_error = 0.0, EH1_error = 0.0;
    INT nits, p0 = 1, p1 = 1000;
    size_t mem, mem_peak;
    int p;

    Unused(analytic);
    pre_refines = 3;
    fn = "cube.dat";
    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-order0", "First order", &p0);
    phgOptionsRegisterInt("-order1", "Last order", &p1);

    phgOptionsPreset("-solver_rtol 1e-30 -solver_btol 0 -solver_atol 0");

    phgOptionsPreset("-solver pcg -pcg_pc_type solver");
    phgOptionsPreset("-pcg_pc_opts \"-solver=mumps\"");
    /*------------------- End command-line options ----------------------*/

    phgVerbosity = 0;
    phgInit(&argc, &argv);
    phgOptionsShowUsed();
    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    /* pre-refine */
    phgRefineAllElements(g, pre_refines);
    phgBalanceGrid(g, 1.0, 1, NULL, 0.);

    f_h = phgDofNew(g, DOF_ANALYTIC, dim_u, "f_h",  func_f);
    e_L2 = phgDofNew(g, DOF_P0, 1, "L2 error indicator", DofNoAction);
    e_H1 = phgDofNew(g, DOF_P0, 1, "H1 error indicator", DofNoAction);
    u_h = phgDofNew(g, DOF_P0, dim_u, "u_h", DofNoAction);
    phgDofSetDataByValue(u_h, 1.0);

    for (p = p0; p <= p1 && DOF_Pn[p] != NULL; p++) {
	double t0, t1;
	phgPrintf("\nTesting uniform p-refinement, order = %d\n", p);
	/* interpolate solution */
	grad_u = phgDofCopy(u_h, NULL, DOF_Pn[p], NULL);
	phgDofFree(&u_h);
	u_h = grad_u;
	if (has_analytic_solution) {
	    DOF *u0 = phgDofNew(g, u_h->type, u_h->dim, "u0", func_u);
	    phgDofAXPY(-1.0, u_h, &u0);
	    phgPrintf("Interp error: %0.3le\n", (double)phgDofNormL2(u0));
	    phgDofFree(&u0);
	}
	elapsed_time(g, FALSE, 0.);	/* reset timer */
	phgPrintf("DOF: %d, T: %d, V: %d, E: %d, F: %d\n",
		  DofDataCountGlobal(u_h),
		  g->nleaf_global, g->nvert_global, g->nedge_global,
		  g->nface_global);
	phgPrintf("Grid distribution: %d submesh%s, LIF = %0.3lg\n",
			g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);

	phgPrintf("Prepare solver: ");
	t0 = phgGetTime(NULL);
	if ((solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL)) == NULL) {
	    phgPrintf("abort.\n");
	    break;
	}
	phgPrintf("LIF %lg ", (double)solver->mat->rmap->lif);
	elapsed_time(g, TRUE, 0.);
    
	phgPrintf("Build linear system: ");
	phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	build_linear_system(solver, u_h, f_h);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	phgPrintf("Assemble linear system: ");
	phgSolverAssemble(solver);
	elapsed_time(g, TRUE, 0.);

	phgPrintf("Solve linear system: ");
	phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	if ((nits = phgSolverSolve(solver, TRUE, u_h, NULL)) < 0) {
	    phgPrintf("abort.\n");
	    break;
	}
	phgPrintf("%d its, %0.4le ", nits, (double)solver->residual);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
	t1 = phgGetTime(NULL);
	phgPrintf("Total solution time: %lg\n", t1 - t0);

	phgSolverDestroy(&solver);

	phgPrintf("Estimate: ");
	grad_u = phgDofGradient(u_h, NULL, NULL, "grad u_h");
	estimate_error(u_h, f_h, grad_u, e_L2, e_H1);
	EL2_error = phgDofNormL2Vec(e_L2);
	EH1_error = phgDofNormL2Vec(e_H1);
	phgPrintf("EL2=%0.2le EH1=%0.2le ", (double)EL2_error,
					    (double)EH1_error);
	if (has_analytic_solution) {
	    DOF *u, *grad;
	    u = phgDofNew(g, u_h->type, u_h->dim, "analytic solution", func_u);
	    grad = phgDofNew(g, u_h->type->grad_type, Dim * u_h->dim,
				"analytic gradient", func_grad);
	    phgDofAXPY(-1.0, u_h, &u);
	    /* Note: the following works only if u is smooth! */
	    phgDofAXPY(-1.0, grad_u, &grad);
	    L2_error = phgDofNormL2(u);
	    H1_error = phgDofNormL2(grad);
	    phgPrintf("L2=%0.2le H1=%0.2le ", (double)L2_error,
			    		      (double)H1_error);
	    phgDofFree(&grad);
	    phgDofFree(&u);
	}
	elapsed_time(g, TRUE, 0.);

	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		mem / (1024.0 * 1024.0), mem_peak / (1024.0 * 1024.0));
	phgDofFree(&grad_u);

	phgQuadClearDofCache(&u_h->type->cache_basfunc, NULL, TRUE);
	phgQuadClearDofCache(&u_h->type->cache_basgrad, NULL, TRUE);
	phgQuadClearDofCache(&u_h->type->cache_gradient, NULL, TRUE);
	phgQuadClearDofCache(&u_h->cache_func, NULL, TRUE);
    }

    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&e_L2);
    phgDofFree(&e_H1);
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
