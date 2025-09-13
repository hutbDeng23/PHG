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

/* $Id: maxwell-ptest.c,v 1.52 2020/12/06 09:03:08 zlb Exp $ */

/* Note: due to a possible bug of GCC on O3800, MIPSPro cc should be
 * used to compile this file when using long double. The code showing
 * the GCC bug is as follows:
#include <math.h>
#include <stdio.h>

int main()
{
    long double x = 1.0, pi = (long double)PI;
    long double a = 4.0 * cosl(2.0 * PI * x), b = 4.0 * cosl(2.0 * pi * x);

    printf("a = %lg, b = %lg\n", (double)a, (double)b);
}
 */

#include "phg.h"

#define _GNU_SOURCE

#include <string.h>
#include <math.h>

/* Standard constants in math.h/quadmath.h: ME[wlq], PI[wlq] */
#define ME	_F(2.7182818284590452353602874713526625)
#define PI	_F(3.1415926535897932384626433832795029)
#define ZERO	_F(0.0)
#define ONE	_F(1.0)
#define TWO	_F(2.0)
#define FOUR	_F(4.0)
#define EIGHT	_F(8.0)
#define HALF	_F(0.5)

/* This code solves the time-harmonic Maxwell equations with Dirichlet BC:

	curl(1/mu curl) E - k^2 E = J,	interior
	E x n = g,			boundary

*/

#define kk	-ONE	/* kk = k^2 */

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
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *(value++) = x * y * z * (x - ONE) * (y - ONE) * (z - ONE) *
		(x - HALF) * (y - HALF) * (z - HALF);
    *(value++) = Sin(TWO * PI * x) * Sin(TWO * PI * y) *
		Sin(TWO * PI * z);
    *(value++) = (ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * (ME - Exp(y)) * (ME - Exp(TWO * y)) *
		(ONE - Exp(z)) * (ME - Exp(z)) * (ME - Exp(TWO * z));
    return;
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    FLOAT u[3];

    func_u(x, y, z, u);

    *(values++) = (-ONE - kk) * u[0] +
	FOUR * Cos(TWO * PI * x) * PI * PI * Cos(TWO * PI * y) *
		Sin(TWO * PI * z) -
	TWO * x * z * (x - ONE) * (z - ONE) * (x - HALF) *
		(y - HALF) * (z - HALF) - TWO * x * z * (x - ONE) * (y - ONE) *
		(z - ONE) * (x - HALF) * (z - HALF) - TWO * x * y * z *
		(x - ONE) * (z - ONE) * (x - HALF) * (z - HALF) -
	TWO * x * y * (x - ONE) * (y - ONE) * (x - HALF) * (y - HALF) *
		(z - HALF) - TWO * x * y * (x - ONE) * (y - ONE) * (z - ONE) *
		(x - HALF) * (y - HALF) - TWO * x * y * z * (x - ONE) *
		(y - ONE) * (x - HALF) * (y - HALF) +
	Exp(x) * (ME - Exp(x)) * (ME - Exp(TWO * x)) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * Exp(z) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) +
	Exp(x) * (ME - Exp(x)) * (ME - Exp(TWO * x)) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		Exp(z) * (ME - Exp(TWO * z)) +
	TWO * Exp(x) * (ME - Exp(x)) * (ME - Exp(TWO * x)) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		(ME - Exp(z)) * Exp(TWO * z) +
	(ONE - Exp(x)) * Exp(x) * (ME - Exp(TWO * x)) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * Exp(z) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) +
	(ONE - Exp(x)) * Exp(x) * (ME - Exp(TWO * x)) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		Exp(z) * (ME - Exp(TWO * z)) +
	TWO * (ONE - Exp(x)) * Exp(x) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * (ME - Exp(y)) * (ME - Exp(TWO * y)) *
		(ONE - Exp(z)) * (ME - Exp(z)) * Exp(TWO * z) +
	TWO * (ONE - Exp(x)) * (ME - Exp(x)) * Exp(TWO * x) *
		(ONE - Exp(y)) * (ME - Exp(y)) * (ME - Exp(TWO * y)) *
		Exp(z) * (ME - Exp(z)) * (ME - Exp(TWO * z)) +
	TWO * (ONE - Exp(x)) * (ME - Exp(x)) * Exp(TWO * x) *
		(ONE - Exp(y)) * (ME - Exp(y)) * (ME - Exp(TWO * y)) *
		(ONE - Exp(z)) * Exp(z) * (ME - Exp(TWO * z)) +
	FOUR * (ONE - Exp(x)) * (ME - Exp(x)) * Exp(TWO * x) *
		(ONE - Exp(y)) * (ME - Exp(y)) * (ME - Exp(TWO * y)) *
		(ONE - Exp(z)) * (ME - Exp(z)) * Exp(TWO * z) +
	x * y * z * (x - ONE) * (y - ONE) * (z - ONE) * (x - HALF) *
		(y - HALF) * (z - HALF);

    *(values++) = (-ONE - kk) * u[1] +
	(ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) * Exp(y) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * Exp(z) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) +
	(ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) * Exp(y) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		Exp(z) * (ME - Exp(TWO * z)) +
	TWO * (ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		Exp(y) * (ME - Exp(y)) * (ME - Exp(TWO * y)) *
		(ONE - Exp(z)) * (ME - Exp(z)) * Exp(TWO * z) +
	(ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * Exp(y) * (ME - Exp(TWO * y)) * Exp(z) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) +
	(ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * Exp(y) * (ME - Exp(TWO * y)) *
		(ONE - Exp(z)) * Exp(z) * (ME - Exp(TWO * z)) +
	TWO * (ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * Exp(y) * (ME - Exp(TWO * y)) *
		(ONE - Exp(z)) * (ME - Exp(z)) * Exp(TWO * z) +
	TWO * (ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * (ME - Exp(y)) * Exp(TWO * y) * Exp(z) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) +
	TWO * (ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * (ME - Exp(y)) * Exp(TWO * y) *
		(ONE - Exp(z)) * Exp(z) * (ME - Exp(TWO * z)) +
	FOUR * (ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * (ME - Exp(y)) * Exp(TWO * y) *
		(ONE - Exp(z)) * (ME - Exp(z)) * Exp(TWO * z) +
	EIGHT * Sin(TWO * PI * x) * Sin(TWO * PI * y) *
		Sin(TWO * PI * z) * PI * PI +
	z * (x - ONE) * (y - ONE) * (z - ONE) * (x - HALF) * (y - HALF) *
		(z - HALF) + x * z * (y - ONE) * (z - ONE) * (x - HALF) *
		(y - HALF) * (z - HALF) +
	x * z * (x - ONE) * (y - ONE) * (z - ONE) * (y - HALF) * (z - HALF) +
	y * z * (x - ONE) * (z - ONE) * (x - HALF) * (y - HALF) * (z - HALF) +
	x * y * z * (z - ONE) * (x - HALF) * (y - HALF) * (z - HALF) +
	x * y * z * (x - ONE) * (z - ONE) * (y - HALF) * (z - HALF) +
	y * z * (x - ONE) * (y - ONE) * (z - ONE) * (x - HALF) * (z - HALF) +
	x * y * z * (y - ONE) * (z - ONE) * (x - HALF) * (z - HALF) +
	x * y * z * (x - ONE) * (y - ONE) * (z - ONE) * (z - HALF) +
	Sin(TWO * PI * x) * Sin(TWO * PI * y) * Sin(TWO * PI * z);

    *(values++) = (-ONE - kk) * u[2] +
	y * (x - ONE) * (y - ONE) * (z - ONE) * (x - HALF) * (y - HALF) *
		(z - HALF) +
	x * y * (y - ONE) * (z - ONE) * (x - HALF) * (y - HALF) * (z - HALF) +
	x * y * (x - ONE) * (y - ONE) * (z - ONE) * (y - HALF) * (z - HALF) +
	y * z * (x - ONE) * (y - ONE) * (x - HALF) * (y - HALF) * (z - HALF) +
	x * y * z * (y - ONE) * (x - HALF) * (y - HALF) * (z - HALF) +
	x * y * z * (x - ONE) * (y - ONE) * (y - HALF) * (z - HALF) +
	y * z * (x - ONE) * (y - ONE) * (z - ONE) * (x - HALF) * (y - HALF) +
	x * y * z * (y - ONE) * (z - ONE) * (x - HALF) * (y - HALF) +
	x * y * z * (x - ONE) * (y - ONE) * (z - ONE) * (y - HALF) +
	Exp(x) * (ME - Exp(x)) * (ME - Exp(TWO * x)) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) -
	TWO * Exp(x) * Exp(x) * (ME - Exp(TWO * x)) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) -
	FOUR * Exp(x) * (ME - Exp(x)) * Exp(TWO * x) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) +
	(ONE - Exp(x)) * Exp(x) * (ME - Exp(TWO * x)) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) -
	FOUR * (ONE - Exp(x)) * Exp(x) * Exp(TWO * x) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) +
	FOUR * (ONE - Exp(x)) * (ME - Exp(x)) * Exp(TWO * x) * (ONE - Exp(y)) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) + (ONE - Exp(x)) *
		(ME - Exp(x)) * (ME - Exp(TWO * x)) * Exp(y) *
		(ME - Exp(y)) * (ME - Exp(TWO * y)) * (ONE - Exp(z)) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) - TWO * (ONE - Exp(x)) *
		(ME - Exp(x)) * (ME - Exp(TWO * x)) * Exp(y) * Exp(y) *
		(ME - Exp(TWO * y)) * (ONE - Exp(z)) * (ME - Exp(z)) *
		(ME - Exp(TWO * z)) - FOUR * (ONE - Exp(x)) * (ME - Exp(x)) *
		(ME - Exp(TWO * x)) * Exp(y) * (ME - Exp(y)) * Exp(TWO * y) *
		(ONE - Exp(z)) * (ME - Exp(z)) * (ME - Exp(TWO * z)) +
	(ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * Exp(y) * (ME - Exp(TWO * y)) *
		(ONE - Exp(z)) * (ME - Exp(z)) * (ME - Exp(TWO * z)) -
	FOUR * (ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * Exp(y) * Exp(TWO * y) * (ONE - Exp(z)) *
		(ME - Exp(z)) * (ME - Exp(TWO * z)) +
	FOUR * (ONE - Exp(x)) * (ME - Exp(x)) * (ME - Exp(TWO * x)) *
		(ONE - Exp(y)) * (ME - Exp(y)) * Exp(TWO * y) *
		(ONE - Exp(z)) * (ME - Exp(z)) * (ME - Exp(TWO * z)) +
	FOUR * Sin(TWO * PI * x) * Cos(TWO * PI * y) * PI * PI *
		Cos(TWO * PI * z) + (ONE - Exp(x)) * (ME - Exp(x)) *
		(ME - Exp(TWO * x)) * (ONE - Exp(y)) * (ME - Exp(y)) *
		(ME - Exp(TWO * y)) * (ONE - Exp(z)) * (ME - Exp(z)) *
		(ME - Exp(TWO * z));
    return;
}

static void
func_curl_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
/*  u = x y z (x-ONE)(y-ONE)(z-ONE)(x-HALF)(y-HALF)(z-HALF)
    v = Sin[2 Pi x] Sin[2 Pi y] Sin[2 Pi z]
    w = (ONE - Exp[x])(E - Exp[x])(E - Exp[TWO x]) (ONE - Exp[y])(E - Exp[y])(E - Exp[TWO y]) (ONE - Exp[z])(E - Exp[z])(E - Exp[TWO z])
    CForm[Simplify[D[w,y] - D[v,z]]]
    CForm[Simplify[D[u,z] - D[w,x]]]
    CForm[Simplify[D[v,x] - D[u,y]]] */
{
    *(value++) = -TWO*Exp(TWO*y)*(ONE - Exp(x))*(ME - Exp(x))*(ME - Exp(TWO*x))*
     (ONE - Exp(y))*(ME - Exp(y))*(ONE - Exp(z))*(ME - Exp(z))*
     (ME - Exp(TWO*z)) - Exp(y)*(ONE - Exp(x))*(ME - Exp(x))*
     (ME - Exp(TWO*x))*(ONE - Exp(y))*(ME - Exp(TWO*y))*
     (ONE - Exp(z))*(ME - Exp(z))*(ME - Exp(TWO*z)) - 
    Exp(y)*(ONE - Exp(x))*(ME - Exp(x))*(ME - Exp(TWO*x))*
     (ME - Exp(y))*(ME - Exp(TWO*y))*(ONE - Exp(z))*(ME - Exp(z))*
     (ME - Exp(TWO*z)) - 2*PI*Cos(2*PI*z)*Sin(2*PI*x)*Sin(2*PI*y);

    *(value++) = TWO*Exp(TWO*x)*(ONE - Exp(x))*(ME - Exp(x))*(ONE - Exp(y))*
     (ME - Exp(y))*(ME - Exp(TWO*y))*(ONE - Exp(z))*(ME - Exp(z))*
     (ME - Exp(TWO*z)) + Exp(x)*(ONE - Exp(x))*(ME - Exp(TWO*x))*
     (ONE - Exp(y))*(ME - Exp(y))*(ME - Exp(TWO*y))*(ONE - Exp(z))*
     (ME - Exp(z))*(ME - Exp(TWO*z)) + 
    Exp(x)*(ME - Exp(x))*(ME - Exp(TWO*x))*(ONE - Exp(y))*
     (ME - Exp(y))*(ME - Exp(TWO*y))*(ONE - Exp(z))*(ME - Exp(z))*
     (ME - Exp(TWO*z)) + (-ONE + x)*(-HALF + x)*x*(-ONE + y)*(-HALF + y)*y*
     (-ONE + z)*(-HALF + z) + (-ONE + x)*(-HALF + x)*x*(-ONE + y)*(-HALF + y)*y*
     (-ONE + z)*z + (-ONE + x)*(-HALF + x)*x*(-ONE + y)*(-HALF + y)*y*(-HALF + z)*z;

    *(value++) =
   -((-ONE + x)*(-HALF + x)*x*(-ONE + y)*(-HALF + y)*(-ONE + z)*(-HALF + z)*z) - 
    (-ONE + x)*(-HALF + x)*x*(-ONE + y)*y*(-ONE + z)*(-HALF + z)*z - 
    (-ONE + x)*(-HALF + x)*x*(-HALF + y)*y*(-ONE + z)*(-HALF + z)*z + 
    TWO*PI*Cos(TWO*PI*x)*Sin(TWO*PI*y)*Sin(TWO*PI*z);

    return;
}

static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
    GRID *g = u_h->g;
    int N = u_h->type->nbas;	/* number of basis functions in an element */
#if USE_OMP
#pragma omp parallel
  {
#endif  /* USE_OMP */
    ELEMENT *e;
    int i, j;
    FLOAT (*A)[N], *B, *buffer;
    INT *I;

    assert(u_h->dim == 1);

    A = phgAlloc(N * sizeof(*A));
    B = phgAlloc(N * sizeof(*B));
    I = phgAlloc(N * sizeof(*I));
    buffer = phgAlloc(N * sizeof(*buffer));
#if USE_OMP
#pragma omp for
#endif	/* USE_OMP */
    ForAllElementsBegin(g, e)
	/* compute \int \curl\phi_j \cdot \curl\phi_i making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (j = 0; j <= i; j++)
		A[j][i] = A[i][j] =
		    phgQuadCurlBasDotCurlBas(e, u_h, j, u_h, i, QUAD_DEFAULT) -
		    kk * phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	}

	/* compute local matrix and RHS */
	for (i = 0; i < N; i++) {	/* loop on basis functions */
	    if (phgDofDirichletBC(u_h, e, i, func_u, 
				buffer, B + i, DOF_PROJ_CROSS)) {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer); 
	    }
	    else {
		/* right hand side: \int f * phi_i */
#if 0
		B[i] = phgQuadFuncDotBas(e, func_f, u_h, i, QUAD_DEFAULT);
#else
		B[i] = phgQuadDofDotBas(e, f_h, u_h, i, QUAD_DEFAULT);
#endif
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
	    }
	}
	phgSolverAddRHSEntries(solver, N, I, B);
    ForAllElementsEnd

    phgFree(A);
    phgFree(B);
    phgFree(I);
    phgFree(buffer);
#if USE_OMP
  }
#endif  /* USE_OMP */
}

static FLOAT
estimate_error(DOF *u_h, DOF *f_h, DOF *curl_u_h, DOF *error)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *jmp1, *jmp2, *res, *div, *tmp;

    tmp = phgDofGetSameOrderDG(u_h, -1, NULL);
    phgDofCopy(f_h, &tmp, NULL, "f-k^2u_h");
    phgDofAXPY(kk, u_h, &tmp);			/* tmp = f_h - k^2 u_h */
    jmp1 = phgQuadFaceJump(curl_u_h, DOF_PROJ_CROSS, NULL, QUAD_DEFAULT);
    jmp2 = phgQuadFaceJump(tmp, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    res = phgDofCurl(curl_u_h, NULL, tmp->type, "residual");
    phgDofAXPBY(ONE, tmp, -ONE, &res);		/* residual */
    div = phgDofDivergence(tmp, NULL, NULL, "div(f+k^2 u_h)");
    phgDofFree(&tmp);
    ForAllElements(g, e) {
	int i;
	FLOAT eta, h;
	FLOAT diam = phgGeomGetDiameter(g, e);
	e->mark = 0;		/* clear refinement mmark */
	eta = ZERO;
	for (i = 0; i < NFace; i++) {
	    INT fno;
	    if (e->bound_type[i] == DIRICHLET || e->bound_type[i] == NEUMANN)
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    fno = e->faces[i];
	    eta += (*DofFaceData(jmp1, fno) + *DofFaceData(jmp2, fno)) * h;
	}
	eta += (phgQuadDofDotDof(e, res, res, QUAD_DEFAULT) +
		phgQuadDofDotDof(e, div, div, QUAD_DEFAULT)) * diam * diam;
	*DofElementData(error, e->index) = Sqrt(eta);
    }
    phgDofFree(&jmp1);
    phgDofFree(&jmp2);
    phgDofFree(&res);
    phgDofFree(&div);
    return phgDofNormInftyVec(error);
}

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *u_h, *f_h = NULL, *curl_u = NULL, *e_Hcurl;
    SOLVER *solver;
    FLOAT L2_error = ZERO, Hcurl_error = ZERO, indicator = ZERO;
    INT nits;
    INT pre_refines = 0, p0 = 1, p1 = 1000;
    const char *fn;
    size_t mem, mem_peak;
    int p;

    pre_refines = 0;
    fn = "cube1o.dat";
    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-order0", "First order", &p0);
    phgOptionsRegisterInt("-order1", "Last order", &p1);

    phgOptionsPreset("-solver_maxit 1000000");
    phgOptionsPreset("-solver_rtol 1e-30 -solver_btol 0 -solver_atol 0");
    phgOptionsPreset("-solver pcg -pcg_pc_type solver");
    phgOptionsPreset("-pcg_pc_opts {"
#if 1	/* AMS with GMRES+ASM as the aux solver (FIXME: FPE (/0) with HC1!) */
		     "	-solver ams"
		     "  -ams_aux_solver_opts {"
		     "	    -solver_maxit 10"
		     "	    -solver_rtol 0."
		     "	    -solver_atol 0."
		     "	    -solver_btol 1e-2"
	#if 1
		     /* use PHG's GMRES solver */
		     "	    -solver gmres"
		     "	    -gmres_pc_type solver"
		     "	    -gmres_pc_opts {"
	#else
		     /* use PETSc's GMRES solver */
		     "	    -solver petsc"
		     "	    -petsc_pc_opts {"
	#endif
		     "		-solver asm"
		     "		-asm_sub_solver_opts {"
		     "		    -solver mumps"
		     "		    -mumps_symmetry spd"
		     "		    -mumps_precision double"
    		     "		}"
		     "		-asm_restriction global"
		     "  	-asm_prolongation local"
		     "  	-asm_overlap 4"
		     "      }"
		     "	    -verbosity 0"
    		     "	}"
#elif 0	/* AMS with MUMPS as the aux solver, poor parallel scalability */
		     "	-solver ams"
		     "  -ams_aux_solver_opts {"
		     "	    -solver mumps"
		     "	    -verbosity 0"
    		     "	}"
#elif 0	/* AMS with default aux solver (hypre BoomerAMG), crashes with HC10 */
		     "	-solver ams"
#else	/* ASM using MUMPS as sub-solver, slow or non convergence above HC13 */
		     "	-solver asm"
		     "	-asm_sub_solver_opts {"
		     "	    -solver mumps"
		     "	    -mumps_symmetry spd"
		     "	    -mumps_precision double"
    		     "	}"
		     "  -asm_restriction global"
		     "  -asm_prolongation local"
		     "  -asm_overlap 8"
#endif
		     "  +ams_warn_aux"
		     "  -ams_default_alpha 1.0"
		     "  -ams_default_beta 1.0"
		     "  -verbosity 0"
		     "}");

    /*------------------- End command-line options ----------------------*/

    phgVerbosity = 0;
    phgInit(&argc, &argv);
    phgOptionsShowUsed();
    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    /* pre-refine */
    phgRefineAllElements(g, pre_refines);
    phgBalanceGrid(g, ONE, 1, NULL, ZERO);

    e_Hcurl = phgDofNew(g, DOF_P0, 1, "Hcurl error indicator", DofNoAction);
    u_h = phgDofNew(g, DOF_HC0, 1, "u_h", DofNoAction);
    phgDofSetDataByValue(u_h, ZERO);

    for (p = p0; p <= p1 && DOF_HCn[p] != NULL; p++) {
	double t0, t1;
	phgPrintf("\nTesting uniform p-refinement, order = %d (%s)\n",
			DOF_HCn[p]->order, DOF_HCn[p]->name);
	/* interpolate solution */
	curl_u = phgDofCopy(u_h, NULL, DOF_HCn[p], NULL);
	phgDofFree(&u_h);
	u_h = curl_u;
	if (FALSE) {
	    DOF *u0 = phgDofNew(g, u_h->type, u_h->dim, "u0", func_u);
	    phgDofAXPY(-ONE, u_h, &u0);
	    phgPrintf("Interp error: %0.3le\n", (double)phgDofNormL2(u0));
	    phgDofFree(&u0);
	}
	phgDofFree(&f_h);
	f_h = phgDofNew(g, /*DOF_HCn[p], 1*/DOF_ANALYTIC, 3, "f_h",  func_f);
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
	phgPrintf("Total solution time: %lg\n", (double)(t1 - t0));

	phgSolverDestroy(&solver);

	phgPrintf("Errors: ");
	curl_u = phgDofCurl(u_h, NULL, NULL, "curl u_h");
	estimate_error(u_h, f_h, curl_u, e_Hcurl);
	indicator = phgDofNormL2Vec(e_Hcurl);
	phgPrintf("indicator=%0.2le ", (double)indicator);
	if (TRUE) {
	    DOF *u, *curl;
	    u = phgDofNew(g, u_h->type, u_h->dim, "analytic solution", func_u);
	    curl = phgDofNew(g, u_h->type->grad_type, DofDim(u_h),
				"analytic curl", func_curl_u);
	    phgDofAXPY(-ONE, u_h, &u);
	    /* Note: the following works only if u is smooth! */
	    phgDofAXPY(-ONE, curl_u, &curl);
	    L2_error = phgDofNormL2(u);
	    Hcurl_error = phgDofNormL2(curl);
	    phgPrintf("L2=%0.2le Hcurl=%0.2le ", (double)L2_error,
			    		      (double)Hcurl_error);
	    phgDofFree(&curl);
	    phgDofFree(&u);
	}
	elapsed_time(g, TRUE, 0.);

	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		mem / (1024.0 * 1024.0), mem_peak / (1024.0 * 1024.0));
	phgDofFree(&curl_u);

	phgQuadClearDofCache(&u_h->type->cache_basfunc, NULL, TRUE);
	phgQuadClearDofCache(&u_h->type->cache_basgrad, NULL, TRUE);
	phgQuadClearDofCache(&u_h->type->cache_gradient, NULL, TRUE);
	phgQuadClearDofCache(&u_h->type->cache_curl, NULL, TRUE);
	phgQuadClearDofCache(&u_h->cache_func, NULL, TRUE);
    }

    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&e_Hcurl);
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
