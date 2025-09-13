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

/* $Id: robin.c,v 1.7 2020/12/06 09:03:08 zlb Exp $ */

/* 
 * This code solves the following Poisson's equation:
 *	- \Delta u(x) = f(x), in \Omega
 * with Robin BC:
 *	(\frac{\partial}{\partial n} + \alpha) u(x) = b(x), on \partial\Omega
 *
 * Note: to extract DOFs and errors from the output:
 *	awk '/DOF/ {n=$2} /Error/ {print n, $4, $7, $10}'
 *
 * Weak form:
 * 	\int_{\Omega} \grad u \grad v + \int_{\partial\Omega} \alpha u v
 * 		= \int_{\Omega} f v + \int_{\partial\Omega} b v
 */

#include "phg.h"

#include <string.h>
#include <math.h>

static FLOAT alpha = 1.0;		/* the constant alpha */

/* unit out normal vector of the current face, must be set before calling the
 * function func_b() (via the DOF b) */
static FLOAT nv[Dim];
#if USE_OMP
#pragma omp threadprivate(nv)
#endif	/* USE_OMP */

static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
/* The analytic solution */
{
    *value = Exp(-10.0 * (x * x + y * y + z * z));
}

static void
func_grad_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
/* Gradient of the analytic solution */
{
    FLOAT u;

    func_u(x, y, z, &u);
    values[0] = -20.0 * x * u;
    values[1] = -20.0 * y * u;
    values[2] = -20.0 * z * u;
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
/* The source (RHS) fucntion f(x) */
{
    FLOAT r2 = x * x + y * y + z * z;

    *value = -(400.0 * r2 - 60.0) * Exp(-10.0 * r2);
}

static void
func_b(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
/* The fucntion b(x) for Robin BC */
{
    FLOAT u, grad_u[Dim];

    func_u(x, y, z, &u);
    func_grad_u(x, y, z, grad_u);
    *value = grad_u[0] * nv[0] + grad_u[1] * nv[1] + grad_u[2] * nv[2] +
	     alpha * u;
}

static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f, DOF *b)
{
    GRID *g = u_h->g;
    ELEMENT *e;

    assert(u_h->dim == 1);
    ForAllElements(g, e) {	/* loop on elements */
	int N = DofNBas(u_h, e);	/* number of bases in the element */
	int ii, jj, face;
	INT i, j;
	FLOAT mat, rhs;

	for (ii = 0; ii < N; ii++) {	/* loop on bases in current element */
	    i = phgSolverMapE2L(solver, 0, e, ii);
	    for (jj = 0; jj < N; jj++) {
		j = phgSolverMapE2L(solver, 0, e, jj);
		mat = phgQuadGradBasDotGradBas(e, u_h, jj, u_h, ii, -1);
		phgSolverAddMatrixEntry(solver, i, j, mat);
	    }
	    phgQuadDofTimesBas(e, f, u_h, ii, -1, &rhs);
	    phgSolverAddRHSEntry(solver, i, rhs);
	}

	/* Check for boundary faces */
	for (face = 0; face < NFace; face++) {
	    if (e->bound_type[face] & BDRY_MASK) {	/* boundary face */
		int n = phgDofGetBasesOnFace(u_h, e, face, NULL);
		SHORT bases[n];
		phgDofGetBasesOnFace(u_h, e, face, bases);
		phgGeomGetFaceOutNormal(g, e, face, nv); /* for func_b() */
		for (ii = 0; ii < n; ii++) {
		    i = phgSolverMapE2L(solver, 0, e, bases[ii]);
		    for (jj = 0; jj < n; jj++) {
			j = phgSolverMapE2L(solver, 0, e, bases[jj]);
			mat = alpha * phgQuadFaceBasDotBas(e, face,
					u_h, bases[jj], u_h, bases[ii], -1);
			phgSolverAddMatrixEntry(solver, i, j, mat);
		    }
		    rhs = phgQuadFaceDofDotBas(e, face,
				b, DOF_PROJ_NONE, u_h, bases[ii], -1);
		    phgSolverAddRHSEntry(solver, i, rhs);
		}
	    }
	}
    }
}

static void
estimate_error(DOF *u_h, DOF *f, DOF *error)
/* compute H1 error indicator */
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *grad, *jump, *residual;

    /* \grad u_h */
    grad = phgDofGradient(u_h, NULL, NULL, NULL);
    /* jump := \int_f [\grad u_h \cdot n], one value on each face */
    jump = phgQuadFaceJump(grad, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    residual = phgDofDivergence(grad, NULL, NULL, NULL);
    phgDofAXPY(1.0, f, &residual);

    ForAllElements(g, e) {
	int i;
	FLOAT eta1, eta2, h;
	/* compute sum of face jumps */
	eta1 = 0.0;
	/* for each face F compute [\grad u_h \cdot n] */
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] & BDRY_MASK)
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);	/* face diameter */
	    eta1 += *DofFaceData(jump, e->faces[i]) * h;
	}
	/* compute residual */
	h = phgGeomGetDiameter(g, e);	/* element diameter */
	eta2 = h * h * phgQuadDofDotDof(e, residual, residual, -1);
	*DofElementData(error, e->index) = eta1 * 0.5 + eta2;
    }
    phgDofFree(&jump);
    phgDofFree(&residual);
    phgDofFree(&grad);

    return;
}

int
main(int argc, char *argv[])
{
    char *fn = "../albert-files/cube5.dat";
				/* input mesh file */
    char *vtk = NULL;		/* VTK file (NULL -> don't write VTK file) */
    FLOAT tol = 0.0;		/* convergence criterion */
    INT mem_max = 200;		/* memory (per process) limit (MB) */

    GRID *g;
    DOF *u_h, *f, *b, *error, *u, *grad_u, *tmp;
    SOLVER *solver;
    double estimator, H1err, L2err;
    size_t mem;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file (input)", (char **)&fn);
    phgOptionsRegisterFilename("-vtk_file", "VTK file (output)", (char **)&vtk);
    phgOptionsRegisterFloat("-tol", "Convergence criterion", &tol);
    phgOptionsRegisterInt("-mem_max", "Maximum memory (MB)", &mem_max);
    phgOptionsRegisterFloat("-alpha", "alpha", &alpha);

    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read mesh file \"%s\".\n", fn);

    /* The numerical solution, default type */
    u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
    u_h->DB_mask = 0;	/* no Dirichlet boundary (important!) */
    phgPrintf("DOF type: %s, mesh file: %s, alpha = %lg\n",
		u_h->type->name, fn, alpha);

    /* RHS function */
    f = phgDofNew(g, DOF_ANALYTIC, 1, "f", func_f);

    /* The function for the Robin boundary condition */
    b = phgDofNew(g, DOF_ANALYTIC, 1, "b",  func_b);

    /* The analytic solution and its gradient */
    u = phgDofNew(g, DOF_ANALYTIC, 1, "u",  func_u);
    grad_u = phgDofNew(g, DOF_ANALYTIC, 3, "grad_u",  func_grad_u);

    /* DOF for storing error indicators, P0 type (1 indicator per element) */
    error = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);

    while (TRUE) {	/* mesh adaptation loop */
	phgPrintf("*** %ld DOF, %ld elements, %d submesh%s\n",
			(long)DofDataCountGlobal(u_h),
			(long)g->nleaf_global,
			g->nprocs, g->nprocs <= 1 ? "" : "es");
	/* This is the only line needed for MPI parallel execution */
	if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("    Repartition mesh, LIF: %lg\n", (double)g->lif);
	/* create solver, with DOF of u_h as its unknowns */
	solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
	/* compute matrix and RHS */
	build_linear_system(solver, u_h, f, b);
	/* solve linear system */
	phgSolverSolve(solver, TRUE, u_h, NULL);
	phgPrintf("    Solve linear system, nits = %d, residual = %le\n",
			solver->nits, (double)solver->residual);
	/* destroy solver */
	phgSolverDestroy(&solver);
	/* compute error indicators */
	estimate_error(u_h, f, error);
	estimator = Sqrt(phgDofNormL1Vec(error));
	phgMemoryUsage(g, &mem);
	/* Compute L2 error */
	tmp = phgDofCopy(u_h, NULL, NULL, NULL);
	phgDofAXPY(-1.0, u, &tmp);
	if (alpha == 0.0) {
	    /* The solution is non-unique up to a constant */
	    DOF *avg = phgDofNew(g, DOF_CONSTANT, 1, "average", DofNoAction);
	    phgDofSetDataByValue(avg, phgDofDotL2(tmp, NULL));
	    phgDofAXPY(-1.0, avg, &tmp);	/* tmp -= \int_\Omega tmp */
	    phgDofFree(&avg);
	}
	L2err = phgDofNormL2(tmp);
	phgDofFree(&tmp);
	/* Compute H1 error */
	tmp = phgDofGradient(u_h, NULL, NULL, NULL);
	phgDofAXPY(-1.0, grad_u, &tmp);
	H1err = phgDofNormL2(tmp);
	H1err = Sqrt(L2err * L2err + H1err * H1err);
	phgDofFree(&tmp);
	phgPrintf("    Errors: est. = %0.3le, H1 = %0.3le, L2 = %0.3le, "
		  "mem = %0.2lfMB\n", estimator, H1err, L2err,
			(double)mem / (1024.0 * 1024.0));
	if (estimator <= tol || mem > 1024 * (size_t)mem_max * 1024)
	    break;
	/* mark elements to be refined */
	phgMarkRefine(MARK_DEFAULT, error, 0.5, NULL, 0., 1,
					Pow(tol, 2) / g->nelem_global);
	/* refine marked elements */
	phgRefineMarkedElements(g);
    }	/* while */

    if (vtk != NULL) {
	/* output a VTK file for postprocessing/visualization */
	phgPrintf("Write final solution to \"%s\".\n", vtk);
	phgExportVTK(g, vtk, u_h, error, NULL);
    }

    phgDofFree(&u_h);
    phgDofFree(&f);
    phgDofFree(&b);
    phgDofFree(&error);
    phgDofFree(&u);
    phgDofFree(&grad_u);

    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
