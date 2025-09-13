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

/* $Id: simplest.c,v 1.117 2020/12/06 09:03:08 zlb Exp $
 *
 * This sample code solves the Helmholtz equation:
 * 	 -\Delta u + a u = f
 * with Dirichlet or periodic boundary conditions */

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <string.h>
#include <math.h>

static FLOAT a = 1.0;

static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = Cos(2. * M_PI * x) * Cos(2. * M_PI * y) * Cos(2. * M_PI * z);
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    func_u(x, y, z, value);
    *value = 12. * M_PI * M_PI * *value + a * *value;
}

static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
    int i, j;
    GRID *g = u_h->g;
    ELEMENT *e;

    assert(u_h->dim == 1);
    ForAllElements(g, e) {
	int N = DofNBas(u_h, e);	/* number of bases in the element */
	FLOAT A[N][N], rhs[N], buffer[N];
	INT I[N];

	/* compute \int \grad\phi_j \cdot \grad\phi_i making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (j = 0; j <= i; j++) {
		A[j][i] = A[i][j] =
		    /* stiffness */
		    phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, QUAD_DEFAULT) +
		    /* mass */
		    a * phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	    }
	}

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, func_u, buffer, rhs+i,
					DOF_PROJ_NONE)) {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer); 
	    }
	    else {	/* interior node */
		/* right hand side = \int f * phi_i */
		phgQuadDofTimesBas(e, f_h, u_h, i, QUAD_DEFAULT, rhs + i);
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
	    }
	}
	phgSolverAddRHSEntries(solver, N, I, rhs);
    }
}

static void
estimate_error(DOF *u_h, DOF *f_h, DOF *grad_u, DOF *error)
/* compute H1 error indicator */
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *jump, *residual, *tmp;

    jump = phgQuadFaceJump(grad_u, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    tmp = phgDofDivergence(grad_u, NULL, NULL, NULL);
    residual = phgDofGetSameOrderDG(u_h, -1, NULL);
    phgDofCopy(f_h, &residual, NULL, NULL); 
    phgDofAXPY(-a, u_h, &residual);
    phgDofAXPY(1.0, tmp, &residual);
    phgDofFree(&tmp);
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
	*DofElementData(error, e->index) = eta;
    }
    phgDofFree(&jump);
    phgDofFree(&residual);

    return;
}

int
main(int argc, char *argv[])
{
    INT periodicity = 0 /* X_MASK | Y_MASK | Z_MASK */;
    INT mem_max = 400;
    FLOAT tol = 5e-1;
    char *fn = "../test/cube4.dat";
    GRID *g;
    ELEMENT *e;
    DOF *u_h, *f_h, *grad_u, *error, *u;
    SOLVER *solver;
    FLOAT L2error, indicator;
    size_t mem_peak;

    phgOptionsRegisterFloat("a", "Coefficient", &a);
    phgOptionsRegisterFloat("tol", "Tolerance", &tol);
    phgOptionsRegisterInt("mem_max", "Maximum memory (MB)", &mem_max);
    phgOptionsRegisterInt("periodicity", "Set periodicity", &periodicity);

    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    phgSetPeriodicity(g, periodicity);

    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    /* The discrete solution */
    if (FALSE) {
	/* u_h is h-p type */
	HP_TYPE *hp = phgHPNew(g, HP_HB);
	ForAllElements(g, e)
	e->hp_order = DOF_DEFAULT->order + GlobalElement(g, e->index) % 3;
	phgHPSetup(hp, FALSE);
	u_h = phgHPDofNew(g, hp, 1, "u_h", DofInterpolation);
	phgHPFree(&hp);
    }
    else {
	/* u_h is non h-p type */
	u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
    }
    phgDofSetDataByValue(u_h, 0.0);

    /* RHS function */
    f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h",  func_f);

    /* DOF for storing a posteriori error estimates */
    error = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);

    /* The analytic solution */
    u = phgDofNew(g, DOF_ANALYTIC, 1, "u", func_u);

    while (TRUE) {
	double t0 = phgGetTime(NULL);
	if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",
			(double)g->lif);
	phgPrintf("%"dFMT" DOF, %"dFMT" elements, %"dFMT
		  " submeshes, load imbalance: %lg\n",
			DofDataCountGlobal(u_h), g->nleaf_global, g->nprocs,
			(double)g->lif);
	solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
	phgPrintf("  DOF: %"dFMT", unknowns: %"dFMT
		  ", Dirichlet bdry: %"dFMT"\n",
		DofDataCountGlobal(u_h), solver->rhs->map->nglobal,
		solver->rhs->map->bdry_nglobal);
	build_linear_system(solver, u_h, f_h);
	phgSolverSolve(solver, TRUE, u_h, NULL);
	phgPrintf("  nits = %d, ", solver->nits);
	phgSolverDestroy(&solver);
	grad_u = phgDofGradient(u_h, NULL, NULL, NULL);
	estimate_error(u_h, f_h, grad_u, error);
	indicator = Sqrt(phgDofNormL1Vec(error));
	phgDofFree(&grad_u);
	grad_u = phgDofCopy(u_h, NULL, NULL, NULL);
	phgDofAXPY(-1.0, u, &grad_u);
	L2error = phgDofNormL2(grad_u);
	phgDofFree(&grad_u);
	phgMemoryUsage(g, &mem_peak);
	phgPrintf("L2 error = %0.3le, indicator = %0.3le, mem = %0.2lfMB\n",
			(double)L2error, (double)indicator,
			(double)mem_peak / (1024.0 * 1024.0));
	phgPrintf("  Wall time: %0.3le\n", phgGetTime(NULL) - t0);
	if (indicator < tol || mem_peak >= 1024 * (size_t)mem_max * 1024)
	    break;
	phgMarkRefine(MARK_DEFAULT, error, Pow(0.8,2), NULL, 0., 1,
			Pow(tol, 2) / g->nleaf_global);
	phgRefineMarkedElements(g);
    }

#if 0
    phgPrintf("Final mesh written to \"%s\".\n",
	phgExportDX(g, "simplest.dx", u_h, error, NULL));
#elif 0
    phgPrintf("Final mesh written to \"%s\".\n",
	phgExportVTK(g, "simplest.vtk", u_h, error, NULL));
#endif

    phgDofFree(&u);
    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&error);

    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
