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

/* $Id: specht.c,v 1.26 2020/12/06 09:03:08 zlb Exp $
 *
 * This sample code tests solving a 2nd order problem:
 * 	 -\Delta u = f, u\in H^1_0(\Omega),
 * or a 4th order problem:
 * 	Delta^2 u = f, u\in H^2_0(\Omega),
 * using the Specht element. */

#include "phg.h"
#include <string.h>
#include <math.h>

/* Note: the Dirichlet BC is enforced with a penalizing parameter \alpha >> 1.
 *    . 2nd order problem - approximate u(x)=b(x) with a Robin type one:
 * 		Du/Dn + \alpha u = b(x)
 *    . 4th order problem - force all boundary DOF to 0 by penalizing the
 *	corresponding equations, i.e., adding the following equation to both
 *	sides of the equations:
 *		\alpha u = 0
 */
static FLOAT alpha = 1e10;

enum {SECOND = 0, FOURTH = 1};
static const char *orders[] = {"2", "4", NULL};
static int order = SECOND;

static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if (order == SECOND)
	*value = Sin(2. * M_PI * x) * Sin(2. * M_PI * y) * Sin(2. * M_PI * z);
    else
	*value = (1. - Cos(2. * M_PI * x)) * (1. - Cos(2. * M_PI * y)) *
		 (1. - Cos(2. * M_PI * z));
}

#if 0
static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if (order == SECOND) {
	func_u(x, y, z, value);
	*value *= 12. * M_PI * M_PI;
    }
    else {
	x *= 2.*M_PI; y *= 2.*M_PI; z *= 2.*M_PI;
	*value = (-(Cos(x)+Cos(y)+Cos(z)) - 9.*Cos(x)*Cos(y)*Cos(z)
		  + 4.*(Cos(x)*Cos(y)+Cos(y)*Cos(z)
			+ Cos(z)*Cos(x))) * Pow(2.*M_PI,4);
    }
}
#endif

static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
    int i, j, face;
    FLOAT d;
    GRID *g = u_h->g;
    ELEMENT *e;

    ForAllElements(g, e) {
	int N = DofNBas(u_h, e);	/* number of bases in the element */
	INT I[N];

	/* Note: I[] maps element DOF to global DOF */
	for (i = 0; i < N; i++)
	    I[i] = phgSolverMapE2L(solver, 0, e, i);

	/* process BC */
	for (face = 0; face < NFace; face++) {
	    if ((e->bound_type[face] & u_h->DB_mask)) {
		/* This is a boundary face, apply BC. */
		int n = phgDofGetBasesOnFace(u_h, e, face, NULL);
		SHORT B[n];
		phgDofGetBasesOnFace(u_h, e, face, B);
		for (i = 0; i < n; i++) {
		    if (order == SECOND) {
			/* 2nd order problem */
			for (j = 0; j < n; j++) {
			    /* A_{ij} += \alpha \int_{face} \phi_j \phi_i */
			    d = alpha * phgQuadFaceBasDotBas(e, face,
                            			u_h, B[j], u_h, B[i], -1);
			    phgSolverAddMatrixEntry(solver, I[B[i]],I[B[j]], d);
			} /* for j */
		    }
		    else {
			/* 4th order problem: A_{ii} += \alpha */
			phgSolverAddMatrixEntry(solver, I[B[i]],I[B[i]], alpha);
		    } /* order == SECOND */
		    /* RHS_i += \int_{face} b(x) \phi_i */
		    /* Nothing to do here since b(x) == 0 */
		} /* for i */
	    } /* if */
	} /* for face */

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    /* interior DOF: A_{ij} += \int_\Omage \phi_j \phi_i */
	    for (j = 0; j < N; j++) {
		if (order == SECOND)
		    phgSolverAddMatrixEntry(solver, I[i], I[j],
			phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, -1)); 
		else
		    phgSolverAddMatrixEntry(solver, I[i], I[j],
			phgQuadLapBasDotLapBas(e, u_h, j, u_h, i, -1)); 
	    }
	    /* RHS_i += \int f * phi_i */
	    phgSolverAddRHSEntry(solver, I[i],
			*phgQuadDofTimesBas(e, f_h, u_h, i, -1, &d));
	}
    }
}

int
main(int argc, char *argv[])
{
    INT mem_max = 1024, refine = 3, nits;
    char *fn = "../albert-files/cube5.dat", kmg;
    GRID *g;
    DOF *u_h, *f_h, *u, *error, *tmp;
    SOLVER *solver;
    size_t mem_peak;
    FLOAT tol = 1e-5, L2err, H1err, mem;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file (input)", (char **)&fn);
    phgOptionsRegisterKeyword("-order", "Problem order", orders, &order);
    phgOptionsRegisterInt("-refine", "Pre-refinement level", &refine);
    phgOptionsRegisterFloat("-tol", "Tolerance", &tol);
    phgOptionsRegisterInt("-mem_max", "Maximum memory (MB)", &mem_max);

    phgOptionsPreset("-dof_type SPT");
    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgRefineAllElements(g, refine);

    phgPrintf("Test solving a %s order problem using the Specht element.\n",
		order == SECOND ? "2nd" : "4th");

    /* The discrete solution */
    u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);

    /* The analytic solution */
    u = phgDofNew(g, order == SECOND ? DOF_P4 : DOF_P6, 1, "u", func_u);

    while (TRUE) {
	double t0 = phgGetTime(NULL);
	phgBalanceGrid(g, 1.2, 1, NULL, 0.);
	solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
	solver->mat->handle_bdry_eqns = FALSE;
	/* create f_h */
#if 0
	/* use analytic expression */
	f_h = phgDofNew(g, /*DOF_ANALYTIC*/DOF_P2, 1, "f_h", func_f);
#else
	/* approximately compute f_h from u */
	tmp = phgDofGradient(u, NULL, DOF_Pn[u->type->grad_type->order], NULL);
	f_h = phgDofDivergence(tmp, NULL,
			DOF_Pn[tmp->type->grad_type->order], NULL);
	phgDofFree(&tmp);
	if (order == SECOND) {
	    /* 2nd order problem: f_h = -\Delta u */
	    phgDofAXPBY(0.0, NULL, -1.0, &f_h);
	}
	else {
	    /* 4th order problem: f_h = \Delta^2 u */
	    tmp = phgDofGradient(f_h, NULL,
			DOF_Pn[f_h->type->grad_type->order], NULL);
	    f_h = phgDofDivergence(tmp, &f_h,
			DOF_Pn[tmp->type->grad_type->order], NULL);
	    phgDofFree(&tmp);
	}
#endif
	/* stiffness matrix and RHS vector */
	build_linear_system(solver, u_h, f_h);
	phgDofFree(&f_h);
	/* solve linear system of equations */
	phgSolverSolve(solver, TRUE, u_h, NULL);
	nits = solver->nits;
	/*FLOAT residual = solver->residual;*/
	phgSolverDestroy(&solver);
	/* compute error = u_h - u */
	error = phgDofCopy(u_h, NULL, u->type, "u_h-u");
	phgDofAXPY(-1.0, u, &error);
	/* L2 error */
	L2err = phgDofNormL2(error);
	/* H1 error */
	tmp = phgDofGradient(error, NULL, NULL, "grad(u_h-u)");
	H1err = Sqrt(Pow(L2err,2) + Pow(phgDofNormL2(tmp),2));
	phgDofFree(&tmp);
	phgDofFree(&error);
	phgMemoryUsage(g, &mem_peak);
	mem = (double)mem_peak;
	if (mem < 1024.0 * 1024.0) {
	    kmg = 'K'; mem /= 1024.0;
	}
	else if (mem < 1024. * 1024. * 1024.) {
	    kmg = 'M'; mem /= 1024.0 * 1024.0;
	}
	else if (mem < 1024. * 1024. * 1024. * 1024.) {
	    kmg = 'G'; mem /= 1024.0 * 1024.0 * 1024.0;
	}
	else {
	    kmg = 'T'; mem /= 1024.0 * 1024.0 * 1024.0 * 1024.0;
	}
	phgPrintf("DOF: %"dFMT"  L2err: %0.2le  H1err: %0.2le  nits: %d  "
		  "mem: %0.2lf%cB  wtime: %0.3le\n",
			DofDataCountGlobal(u_h),
			(double)L2err, (double)H1err, (int)nits,
			(double)mem, kmg, phgGetTime(NULL) - t0);
	if (H1err <= tol || mem_peak >= 1024 * (size_t)mem_max * 1024)
	    break;
	phgRefineAllElements(g, 1);
    }

    phgDofFree(&u);
    phgDofFree(&u_h);

    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
