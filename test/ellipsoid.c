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

/* $Id: ellipsoid.c,v 1.12 2020/12/06 09:03:08 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>

static void
build_linear_system(SOLVER *solver, DOF *u, DOF *f, DOF *a)
{
    int N = u->type->nbas;	/* number of basis functions in an element */
    int i, j;
    GRID *g = u->g;
    ELEMENT *e;
    FLOAT A[N][N], rhs[N], buffer[N];
    INT I[N];

    assert(u->dim == 1);
    ForAllElements(g, e) {
	/* compute \int \grad\phi_j \cdot \grad\phi_i making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (j = 0; j <= i; j++)
		A[j][i] = A[i][j] =
		    phgQuadGradBasAGradBas(e, u, j, a, u, i, QUAD_DEFAULT);
	}

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u, e, i, NULL, buffer, rhs+i,
					DOF_PROJ_NONE)) {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer); 
	    }
	    else {	/* interior node */
		/* right hand side = \int f * phi_i */
		phgQuadDofTimesBas(e, f, u, i, QUAD_DEFAULT, rhs + i);
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
	    }
	}
	phgSolverAddRHSEntries(solver, N, I, rhs);
    }
}

static FLOAT
estimate_error(DOF *u, DOF *f, DOF *grad_u, DOF *a, DOF *error)
/* compute H1 error indicator */
{
    GRID *g = u->g;
    ELEMENT *e;
    DOF *jump, *residual;

    grad_u = phgDofMM(FALSE, FALSE, 1, 3, 1, 1.0, a, 0, grad_u, 0.0, NULL);
    jump = phgQuadFaceJump(grad_u, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    residual = phgDofDivergence(grad_u, NULL, NULL, NULL);
    phgDofAXPY(1.0, f, &residual);
    ForAllElements(g, e) {
	int i;
	FLOAT eta, h;
	FLOAT diam = phgGeomGetDiameter(g, e);
	eta = 0.0;
	/* for each face F compute [grad_u \cdot n] */
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] & (DIRICHLET | NEUMANN))
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    eta += *DofFaceData(jump, e->faces[i]) * h;
	}
	eta = eta * .5 + diam * diam * phgQuadDofDotDof(e, residual, residual,
								QUAD_DEFAULT);
	*DofElementData(error, e->index) = Sqrt(eta);
    }
    phgDofFree(&grad_u);
    phgDofFree(&jump);
    phgDofFree(&residual);
    return phgDofNormInftyVec(error);
}

int
main(int argc, char *argv[])
{
    ELEMENT *e;
    GRID *g;
    const char *fn = "ellipsoid.dat";
    FLOAT vol;
    DOF *u, *a, *grad_u, *f, *error;
    SOLVER *solver;
    FLOAT PE_error = 0.0, thres, tol = 1e-6;
    INT nits, mem_max = 400;
    size_t mem, mem_peak;

    phgOptionsRegisterInt("mem_max", "memory size limit", &mem_max);
    phgOptionsRegisterFloat("tol", "convergence tolerance", &tol);
    phgOptionsPreset("-dof_type P1");
    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    vol = 0.0;
    ForAllElements(g, e)
	vol += phgGeomGetVolume(g, e);
    phgPrintf("Initial mesh: volume = %lg\n", (double)vol);

    f = phgDofNew(g, DOF_P0, 1, "f", DofInterpolation);
    a = phgDofNew(g, DOF_P0, 1, "a", DofInterpolation);
    ForAllElements(g, e) {
	/* f == constant inside/outside the ellipsoid */
	static FLOAT lam[] = {0.25, 0.25, 0.25, 0.25};
	FLOAT x, y, z, r;
	phgGeomLambda2XYZ(g, e, lam, &x, &y, &z);
	x -= 0.5;
	y -= 0.5;
	z -= 0.5;
	r = x * x * 16. + y * y * 25. + z * z * 25.;
	*DofElementData(a, e->index) = (r < 1.0 ? 1.0 : 0.5);
	*DofElementData(f, e->index) = (r < 1.0 ? -10.0 : 5.0);
    }

    u = phgDofNew(g, DOF_DEFAULT, 1, "u", DofInterpolation);
    /*phgDofSetDirichletBoundaryMask(u, BDRY_MASK);*/
    error = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);

    while (TRUE) {
	phgPrintf("%d DOF, %d elements, %d submeshes, load imbalance: %lg\n",
			DofDataCountGlobal(u), g->nleaf_global, g->nprocs,
			(double)g->lif);
	if (phgBalanceGrid(g, 1.2, -1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",
			(double)g->lif);
	solver = phgSolverCreate(SOLVER_DEFAULT, u, NULL);
	build_linear_system(solver, u, f, a);
	phgSolverSolve(solver, TRUE, u, NULL);
	nits = solver->nits;
	phgSolverDestroy(&solver);
	grad_u = phgDofGradient(u, NULL, NULL, NULL);
	PE_error = estimate_error(u, f, grad_u, a, error);
	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("  nit = %d, error indicator = %0.3le, memory = %0.2lfMB\n",
			nits, (double)PE_error, mem_peak / (1024.0 * 1024.0));
	if (PE_error < tol || mem_peak >= 1024 * (size_t)mem_max * 1024) {
	    phgPrintf("Final mesh written to \"%s\".\n",
		phgExportVTK(g, "ellipsoid.vtk", u, grad_u, a, f, error, NULL));
	    phgDofFree(&grad_u);
	    break;
	}
	phgDofFree(&grad_u);
	/*thres = tol / Sqrt((double)g->nleaf_global);*/
	thres = 0.25 * PE_error;
	ForAllElements(g, e)
	    e->mark = (*DofElementData(error, e->index) > thres) ? 1 : 0;
	phgRefineMarkedElements(g);
    }

    phgDofFree(&u);
    phgDofFree(&f);
    phgDofFree(&error);
    phgDofFree(&a);
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
