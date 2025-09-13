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

/* This code tests the H(div) bases by solving the following grad-div type
 * equation:
 *
 *	<\div u,\div v> + tau <u,v> = <f,v>, u,v \in H(\div)
 *
 * $Id: HD_test.c,v 1.14 2020/12/20 10:33:56 zlb Exp $
 */

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <math.h>
#include <stdlib.h>

static FLOAT tau = 1.0;

#define E	_F(2.7182818284590452353602874713526625)
#define Pi	_F(3.1415926535897932384626433832795029)

/* The analytic solution is simply taken from examples/maxwell.c, which has
 * zero trace on the boundary faces of cube1o.dat (the Fichera corner).
 *
 * u1 = x y z (x-1) (y-1) (z-1) (x-1/2) (y-1/2) (z-1/2)
 * u2 = Sin[2 Pi x] Sin[2 Pi y] Sin[2 Pi z]
 * u3 = (1-E^x)(E-E^x)(E-E^(2x))(1-E^y)(E-E^y)(E-E^(2y))(1-E^z)(E-E^z)(E-E^(2z))
 *
 * div(u):
 * 	CForm[Simplify[D[u1,x]+D[u2,y]+D[u3,z]]]
 *
 * grad(div(u)):
 * 	CForm[Simplify[D[D[u1,x]+D[u2,y]+D[u3,z],x]]]
 *	CForm[Simplify[D[D[u1,x]+D[u2,y]+D[u3,z],y]]]
 *	CForm[Simplify[D[D[u1,x]+D[u2,y]+D[u3,z],z]]]
 */

static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *(value++) = x * y * z * (x - 1.0) * (y - 1.0) * (z - 1.0) *
		(x - .5) * (y - .5) * (z - .5);
    *(value++) = Sin(2.0 * Pi * x) * Sin(2.0 * Pi * y) *
		Sin(2.0 * Pi * z);
    *(value++) = (1.0 - Exp(x)) * (E - Exp(x)) * (E - Exp(x + x)) *
		 (1.0 - Exp(y)) * (E - Exp(y)) * (E - Exp(y + y)) *
		 (1.0 - Exp(z)) * (E - Exp(z)) * (E - Exp(z + z));
    return;
}

static void
func_div_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = -2*Exp(2*z)*(-1 + Exp(x))*(-E + Exp(x))*(-E + Exp(2*x))*
	     (-1 + Exp(y))*(-E + Exp(y))*(-E + Exp(2*y))*
	     (-1 + Exp(z))*(-E + Exp(z)) - 
	    Exp(z)*(-1 + Exp(x))*(-E + Exp(x))*(-E + Exp(2*x))*
	     (-1 + Exp(y))*(-E + Exp(y))*(-E + Exp(2*y))*
	     (-1 + Exp(z))*(-E + Exp(2*z)) - 
	    Exp(z)*(-1 + Exp(x))*(-E + Exp(x))*(-E + Exp(2*x))*
	     (-1 + Exp(y))*(-E + Exp(y))*(-E + Exp(2*y))*
	     (-E + Exp(z))*(-E + Exp(2*z)) + 
	    (-1 + x)*(-0.5 + x)*(-1 + y)*(-0.5 + y)*y*(-1 + z)*(-0.5 + z)*z + 
	    (-1 + x)*x*(-1 + y)*(-0.5 + y)*y*(-1 + z)*(-0.5 + z)*z + 
	    (-0.5 + x)*x*(-1 + y)*(-0.5 + y)*y*(-1 + z)*(-0.5 + z)*z + 
	    2*Pi*Cos(2*Pi*y)*Sin(2*Pi*x)*Sin(2*Pi*z);
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
/* Note: f = -\grad\div u + \tau * u */
{
    FLOAT u[3];

    func_u(x, y, z, u);

    *(values++) = tau * u[0] -
	  (-4*Exp(2*(x + z))*(-1 + Exp(x))*(-E + Exp(x))*
	     (-1 + Exp(y))*(-E + Exp(y))*(-E + Exp(2*y))*
	     (-1 + Exp(z))*(-E + Exp(z)) - 
	    2*Exp(x + 2*z)*(-1 + Exp(x))*(-E + Exp(2*x))*
	     (-1 + Exp(y))*(-E + Exp(y))*(-E + Exp(2*y))*
	     (-1 + Exp(z))*(-E + Exp(z)) - 
	    2*Exp(x + 2*z)*(-E + Exp(x))*(-E + Exp(2*x))*
	     (-1 + Exp(y))*(-E + Exp(y))*(-E + Exp(2*y))*
	     (-1 + Exp(z))*(-E + Exp(z)) - 
	    2*Exp(2*x + z)*(-1 + Exp(x))*(-E + Exp(x))*(-1 + Exp(y))*
	     (-E + Exp(y))*(-E + Exp(2*y))*(-1 + Exp(z))*
	     (-E + Exp(2*z)) - Exp(x + z)*(-1 + Exp(x))*
	     (-E + Exp(2*x))*(-1 + Exp(y))*(-E + Exp(y))*
	     (-E + Exp(2*y))*(-1 + Exp(z))*(-E + Exp(2*z)) - 
	    Exp(x + z)*(-E + Exp(x))*(-E + Exp(2*x))*(-1 + Exp(y))*
	     (-E + Exp(y))*(-E + Exp(2*y))*(-1 + Exp(z))*
	     (-E + Exp(2*z)) - 2*Exp(2*x + z)*(-1 + Exp(x))*
	     (-E + Exp(x))*(-1 + Exp(y))*(-E + Exp(y))*
	     (-E + Exp(2*y))*(-E + Exp(z))*(-E + Exp(2*z)) - 
	    Exp(x + z)*(-1 + Exp(x))*(-E + Exp(2*x))*(-1 + Exp(y))*
	     (-E + Exp(y))*(-E + Exp(2*y))*(-E + Exp(z))*
	     (-E + Exp(2*z)) - Exp(x + z)*(-E + Exp(x))*
	     (-E + Exp(2*x))*(-1 + Exp(y))*(-E + Exp(y))*
	     (-E + Exp(2*y))*(-E + Exp(z))*(-E + Exp(2*z)) + 
	    2*(-1 + x)*(-1 + y)*(-0.5 + y)*y*(-1 + z)*(-0.5 + z)*z + 
	    2*(-0.5 + x)*(-1 + y)*(-0.5 + y)*y*(-1 + z)*(-0.5 + z)*z + 
	    2*x*(-1 + y)*(-0.5 + y)*y*(-1 + z)*(-0.5 + z)*z + 
	    4*Pow(Pi,2)*Cos(2*Pi*x)*Cos(2*Pi*y)*Sin(2*Pi*z));

    *(values++) = tau * u[1] -
	    ((-8*Exp(y + z)*(-1 + Exp(x))*
	        (Exp(2) + Exp(3*x) - Exp(1 + x) - Exp(1 + 2*x))*
	        (E + Exp(2) - 3*Exp(2*y) + 4*Exp(3*y) - 
	          3*Exp(1 + 2*y))*(E + Exp(2) - 3*Exp(2*z) + 
	          4*Exp(3*z) - 3*Exp(1 + 2*z)) + 
	       (1 - 6*x + 6*Pow(x,2))*(1 - 6*y + 6*Pow(y,2))*z - 
	       3*(1 - 6*x + 6*Pow(x,2))*(1 - 6*y + 6*Pow(y,2))*Pow(z,2) + 
	       2*(1 - 6*x + 6*Pow(x,2))*(1 - 6*y + 6*Pow(y,2))*Pow(z,3))/8. - 
	    4*Pow(Pi,2)*Sin(2*Pi*x)*Sin(2*Pi*y)*Sin(2*Pi*z));

    *(values++) = tau * u[2] -
	    (-4*Exp(3*z)*(-1 + Exp(x))*(-E + Exp(x))*(-E + Exp(2*x))*
	     (-1 + Exp(y))*(-E + Exp(y))*(-E + Exp(2*y))*(-1 + Exp(z))
	      - 4*Exp(3*z)*(-1 + Exp(x))*(-E + Exp(x))*
	     (-E + Exp(2*x))*(-1 + Exp(y))*(-E + Exp(y))*
	     (-E + Exp(2*y))*(-E + Exp(z)) - 
	    4*Exp(2*z)*(-1 + Exp(x))*(-E + Exp(x))*(-E + Exp(2*x))*
	     (-1 + Exp(y))*(-E + Exp(y))*(-E + Exp(2*y))*
	     (-1 + Exp(z))*(-E + Exp(z)) - 
	    2*Exp(2*z)*(-1 + Exp(x))*(-E + Exp(x))*(-E + Exp(2*x))*
	     (-1 + Exp(y))*(-E + Exp(y))*(-E + Exp(2*y))*
	     (-E + Exp(2*z)) - Exp(z)*(-1 + Exp(x))*(-E + Exp(x))*
	     (-E + Exp(2*x))*(-1 + Exp(y))*(-E + Exp(y))*
	     (-E + Exp(2*y))*(-1 + Exp(z))*(-E + Exp(2*z)) - 
	    Exp(z)*(-1 + Exp(x))*(-E + Exp(x))*(-E + Exp(2*x))*
	     (-1 + Exp(y))*(-E + Exp(y))*(-E + Exp(2*y))*
	     (-E + Exp(z))*(-E + Exp(2*z)) + 
	    (-1 + x)*(-0.5 + x)*(-1 + y)*(-0.5 + y)*y*(-1 + z)*(-0.5 + z) + 
	    (-1 + x)*x*(-1 + y)*(-0.5 + y)*y*(-1 + z)*(-0.5 + z) + 
	    (-0.5 + x)*x*(-1 + y)*(-0.5 + y)*y*(-1 + z)*(-0.5 + z) + 
	    (-1 + x)*(-0.5 + x)*(-1 + y)*(-0.5 + y)*y*(-1 + z)*z + 
	    (-1 + x)*x*(-1 + y)*(-0.5 + y)*y*(-1 + z)*z + 
	    (-0.5 + x)*x*(-1 + y)*(-0.5 + y)*y*(-1 + z)*z + 
	    (-1 + x)*(-0.5 + x)*(-1 + y)*(-0.5 + y)*y*(-0.5 + z)*z + 
	    (-1 + x)*x*(-1 + y)*(-0.5 + y)*y*(-0.5 + z)*z + 
	    (-0.5 + x)*x*(-1 + y)*(-0.5 + y)*y*(-0.5 + z)*z + 
	    4*Pow(Pi,2)*Cos(2*Pi*y)*Cos(2*Pi*z)*Sin(2*Pi*x));
}

#undef E
#undef Pi

static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    assert(u_h->dim == 1);

#if USE_OMP
# pragma omp parallel for private(e)
#endif	/* USE_OMP */
    ForAllElementsBegin(g, e) {	/* \PHGindex{ForAllElements} */
	int N = DofNBas(u_h, e); /* number of basis functions in the element */
	int i, j;
	FLOAT A[N][N], B[N], buffer[N];
	INT I[N];

	/* compute \int \div\phi_j \cdot \div\phi_i, making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (j = 0; j <= i; j++)
		A[j][i] = A[i][j] =
		    phgQuadDivBasDotDivBas(e, u_h, j, u_h, i, QUAD_DEFAULT) +
		    tau * phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	}

	/* compute local matrix and RHS */
	for (i = 0; i < N; i++) {	/* loop on basis functions */
	    if (phgDofDirichletBC(u_h, e, i, func_u, 
				buffer, B + i, DOF_PROJ_DOT)) {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer); 
	    }
	    else {
		/* right hand side: \int f \cdot phi_i */
#if 0
		B[i] = phgQuadFuncDotBas(e, func_f, u_h, i, QUAD_DEFAULT);
#else
		B[i] = phgQuadDofDotBas(e, f_h, u_h, i, QUAD_DEFAULT);
#endif
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
	    }
	}
	phgSolverAddRHSEntries(solver, N, I, B);
    } ForAllElementsEnd
}

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *u_h, *u, *div_u, *f_h, *tmp;
    SOLVER *solver;
    FLOAT error, a, b;

    char *fn = "cube1o.dat";
    INT pre_refines = 0;
    INT depth = 1;			/* <=0 => uniform refine -depth or 3 */
    INT mem_max = 400;			/* max. memory */
    FLOAT tol = 1e-10;			/* convergence tolerance */
    size_t mem, mem_peak;
    INT submesh_threshold = 1;
    FLOAT lif_threshold = 1.2;


    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **) &fn);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-refine_depth", "Refinement depth", &depth);
    phgOptionsRegisterInt("-mem_max", "Max memory (MB)", &mem_max);
    phgOptionsRegisterFloat("-tau", "The value of tau", &tau);
    phgOptionsRegisterFloat("-lif_threshold", "Load imbalance factor threshold",
				&lif_threshold);
    phgOptionsRegisterInt("-submesh_threshold", "Submesh size threshold",
				&submesh_threshold);

    phgOptionsPreset("-dof_type HD1");

    /* By default choose a direct solver if available */
#if USE_MUMPS
    phgOptionsPreset("-solver mumps");
#elif USE_PARDISO
    phgOptionsPreset("-solver pardiso");
#elif USE_SPOOLES
    phgOptionsPreset("-solver spooles");
#elif USE_SUPERLU
    phgOptionsPreset("-solver superlu");
#elif USE_USE_SSPARSE
    phgOptionsPreset("-solver ssparse");
#endif	/* USE_MUMPS */

    phgInit(&argc, &argv);
    phgOptionsShowUsed();

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    phgRefineAllElements(g, pre_refines);

    /* The analytic solution and its divergence */
    u = phgDofNew(g, DOF_Pn[DOF_DEFAULT->order], 3, "u", func_u);
    div_u = phgDofNew(g, u->type, 1, "div_u", func_div_u);

    u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation/*DofNoAction*/);
    u_h->DB_mask = BDRY_MASK;	/* Dirichlet BC only */
#if 0
    f_h = phgDofNew(g, DOF_ANALYTIC, 3, "f_h", func_f);
#else
    f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h", func_f);
#endif

    phgPrintf("DOF types: u_h = %s, f_h = %s, u and div_u = %s\n",
		u_h->type->name, f_h->type->name, u->type->name);

    while (TRUE) {
	double t0;
	phgPrintf("\n------ %"dFMT" DOF, %"dFMT" elements, mesh LIF = %g\n",
		  DofDataCountGlobal(u_h), g->nleaf_global, (double)g->lif);

	if (phgBalanceGrid(g, lif_threshold, submesh_threshold, NULL, 0.)) {
	    phgPrintf("------ Repartition mesh: nprocs = %d, LIF = %lg\n",
			    g->nprocs, (double)g->lif);
	}

#if 0
DOF *div = phgDofDivergence(u, NULL, NULL, NULL);
tmp = phgDofGradient(div, NULL, NULL, NULL);
phgDofFree(&div);
phgDofCopy(u, &div, NULL, NULL);
phgDofAXPBY(-1.0, tmp, tau, &div);
phgDofFree(&tmp);
phgDofAXPBY(-1.0, f_h, 1.0, &div);
phgPrintf("|grad(div(u)) - f| = %lg\n", (double)phgDofNormL2(div));
phgDofFree(&div);
error = 1e10;
goto skip;
#endif

	phgPrintf("Set up linear solver: ");
	solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
	phgPrintf("solver LIF = %lg\n", (double)solver->mat->rmap->lif);
	phgPrintf("Build linear system: ");
	t0 = phgGetTime(NULL);
	build_linear_system(solver, u_h, f_h);
	phgPrintf("time = %lg\n", phgGetTime(NULL) - t0);

	phgPrintf("Solve linear system: ");
	t0 = phgGetTime(NULL);
	phgSolverSolve(solver, TRUE, u_h, NULL);
	phgPrintf("nits = %d, resid = %0.4lg, time = %g\n", solver->nits,
			(double)solver->residual, phgGetTime(NULL) - t0);
	phgSolverDestroy(&solver);

	phgPrintf("Errors: ");
	t0 = phgGetTime(NULL);
	tmp = phgDofCopy(u, NULL, u_h->type, NULL);
	phgDofAXPY(-1.0, u_h, &tmp);
	phgPrintf("L2 = %0.8le, ", (double)(error = phgDofNormL2(tmp)));
	phgDofFree(&tmp);
	tmp = phgDofDivergence(u_h, NULL, DOF_DGn[div_u->type->order], NULL);
	phgDofAXPY(-1.0, div_u, &tmp);
	phgPrintf("H(div) = %0.8le, ", (double)(error += phgDofNormL2(tmp)));
	phgPrintf("time = %lg\n", phgGetTime(NULL) - t0);
	phgDofFree(&tmp);
	goto skip;	/* to make gcc happy */
skip:
	mem = phgMemoryUsage(g, &mem_peak);
	a = mem / (1024.0 * 1024.0);
	b = mem_peak / (1024.0 * 1024.0);
	phgPrintf("Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		(double)a, (double)b);
	if (error < tol || mem_peak >= 1024 * (size_t)mem_max * 1024) {
	    break;
	}

	phgPrintf("Refine mesh ");
	phgRefineAllElements(g, depth);
	phgPrintf("\n");
    }

    //phgExportVTK(g, "HD_test.vtk", u_h, u, div_u, NULL);

    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&u);
    phgDofFree(&div_u);

    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
