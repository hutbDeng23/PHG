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

/* Testing solving time-harmonic Maxwell equations with an analytic solution.
 *
 *	curl(1/mu curl) E - k^2 E = J
 *
 * $Id: maxwell-test.c,v 1.38 2020/12/06 09:03:08 zlb Exp $
 */

#include "phg.h"

#include <string.h>
#include <math.h>
#include <stdlib.h>

#define C 0.1                   //Scale eta
#define EC 5.8e7                /*electric conductivity of copper */
#define MP (4.0e-7*M_PI)        /*magnetic permeability of vacuum */
#define SC 1e-5                 /*scale parameter */
#define HZ (1.0e9*2.0*M_PI)
#define ss (SC*SC*HZ*EC*MP)

static void
func_kk(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = (x < 0.5 /*&& y < 0.5 && z < 0.5*/ ? 0.0 : -1.0);
//*value = -1.0;
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
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *(value++) = x * y * z * (x - 1.0) * (y - 1.0) * (z - 1.0) *
		(x - .5) * (y - .5) * (z - .5);
    *(value++) = Sin(2.0 * M_PI * x) * Sin(2.0 * M_PI * y) *
		Sin(2.0 * M_PI * z);
    *(value++) = (1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * (M_E - Exp(y)) * (M_E - Exp(2.0 * y)) *
		(1.0 - Exp(z)) * (M_E - Exp(z)) * (M_E - Exp(2.0 * z));
    return;
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    FLOAT u[3];
    FLOAT kk;

    /* Note: the original expression equals curlcurl u + u, by adding (-1-kk) u
     * to it we obtain curlcurl u - kk u */

    func_kk(x, y, z, &kk);
/**(values++) = kk;
*(values++) = 0.0;
*(values++) = 0.0;
return;*/
    func_u(x, y, z, u);

    *(values++) = (-1. - kk) * u[0] +
	4.0 * Cos(2.0 * M_PI * x) * M_PI * M_PI * Cos(2.0 * M_PI * y) *
		Sin(2.0 * M_PI * z) -
	2.0 * x * z * (x - 1.0) * (z - 1.0) * (x - 0.5) *
		(y - 0.5) * (z - 0.5) - 2.0 * x * z * (x - 1.0) * (y - 1.0) *
		(z - 1.0) * (x - 0.5) * (z - 0.5) - 2.0 * x * y * z *
		(x - 1.0) * (z - 1.0) * (x - 0.5) * (z - 0.5) -
	2.0 * x * y * (x - 1.0) * (y - 1.0) * (x - 0.5) * (y - 0.5) *
		(z - 0.5) - 2.0 * x * y * (x - 1.0) * (y - 1.0) * (z - 1.0) *
		(x - 0.5) * (y - 0.5) - 2.0 * x * y * z * (x - 1.0) *
		(y - 1.0) * (x - 0.5) * (y - 0.5) +
	Exp(x) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * Exp(z) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	Exp(x) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		Exp(z) * (M_E - Exp(2.0 * z)) +
	2.0 * Exp(x) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		(M_E - Exp(z)) * Exp(2.0 * z) +
	(1.0 - Exp(x)) * Exp(x) * (M_E - Exp(2.0 * x)) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * Exp(z) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	(1.0 - Exp(x)) * Exp(x) * (M_E - Exp(2.0 * x)) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		Exp(z) * (M_E - Exp(2.0 * z)) +
	2.0 * (1.0 - Exp(x)) * Exp(x) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * (M_E - Exp(y)) * (M_E - Exp(2.0 * y)) *
		(1.0 - Exp(z)) * (M_E - Exp(z)) * Exp(2.0 * z) +
	2.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * Exp(2.0 * x) *
		(1.0 - Exp(y)) * (M_E - Exp(y)) * (M_E - Exp(2.0 * y)) *
		Exp(z) * (M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	2.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * Exp(2.0 * x) *
		(1.0 - Exp(y)) * (M_E - Exp(y)) * (M_E - Exp(2.0 * y)) *
		(1.0 - Exp(z)) * Exp(z) * (M_E - Exp(2.0 * z)) +
	4.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * Exp(2.0 * x) *
		(1.0 - Exp(y)) * (M_E - Exp(y)) * (M_E - Exp(2.0 * y)) *
		(1.0 - Exp(z)) * (M_E - Exp(z)) * Exp(2.0 * z) +
	x * y * z * (x - 1.0) * (y - 1.0) * (z - 1.0) * (x - 0.5) *
		(y - 0.5) * (z - 0.5);

    *(values++) = (-1. - kk) * u[1] +
	(1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) * Exp(y) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * Exp(z) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	(1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) * Exp(y) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		Exp(z) * (M_E - Exp(2.0 * z)) +
	2.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		Exp(y) * (M_E - Exp(y)) * (M_E - Exp(2.0 * y)) *
		(1.0 - Exp(z)) * (M_E - Exp(z)) * Exp(2.0 * z) +
	(1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * Exp(y) * (M_E - Exp(2.0 * y)) * Exp(z) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	(1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * Exp(y) * (M_E - Exp(2.0 * y)) *
		(1.0 - Exp(z)) * Exp(z) * (M_E - Exp(2.0 * z)) +
	2.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * Exp(y) * (M_E - Exp(2.0 * y)) *
		(1.0 - Exp(z)) * (M_E - Exp(z)) * Exp(2.0 * z) +
	2.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * (M_E - Exp(y)) * Exp(2.0 * y) * Exp(z) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	2.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * (M_E - Exp(y)) * Exp(2.0 * y) *
		(1.0 - Exp(z)) * Exp(z) * (M_E - Exp(2.0 * z)) +
	4.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * (M_E - Exp(y)) * Exp(2.0 * y) *
		(1.0 - Exp(z)) * (M_E - Exp(z)) * Exp(2.0 * z) +
	8.0 * Sin(2.0 * M_PI * x) * Sin(2.0 * M_PI * y) *
		Sin(2.0 * M_PI * z) * M_PI * M_PI +
	z * (x - 1.0) * (y - 1.0) * (z - 1.0) * (x - 0.5) * (y - 0.5) *
		(z - 0.5) + x * z * (y - 1.0) * (z - 1.0) * (x - 0.5) *
		(y - 0.5) * (z - 0.5) +
	x * z * (x - 1.0) * (y - 1.0) * (z - 1.0) * (y - 0.5) * (z - 0.5) +
	y * z * (x - 1.0) * (z - 1.0) * (x - 0.5) * (y - 0.5) * (z - 0.5) +
	x * y * z * (z - 1.0) * (x - 0.5) * (y - 0.5) * (z - 0.5) +
	x * y * z * (x - 1.0) * (z - 1.0) * (y - 0.5) * (z - 0.5) +
	y * z * (x - 1.0) * (y - 1.0) * (z - 1.0) * (x - 0.5) * (z - 0.5) +
	x * y * z * (y - 1.0) * (z - 1.0) * (x - 0.5) * (z - 0.5) +
	x * y * z * (x - 1.0) * (y - 1.0) * (z - 1.0) * (z - 0.5) +
	Sin(2.0 * M_PI * x) * Sin(2.0 * M_PI * y) * Sin(2.0 * M_PI * z);

    *(values++) = (-1. - kk) * u[2] +
	y * (x - 1.0) * (y - 1.0) * (z - 1.0) * (x - 0.5) * (y - 0.5) *
		(z - 0.5) +
	x * y * (y - 1.0) * (z - 1.0) * (x - 0.5) * (y - 0.5) * (z - 0.5) +
	x * y * (x - 1.0) * (y - 1.0) * (z - 1.0) * (y - 0.5) * (z - 0.5) +
	y * z * (x - 1.0) * (y - 1.0) * (x - 0.5) * (y - 0.5) * (z - 0.5) +
	x * y * z * (y - 1.0) * (x - 0.5) * (y - 0.5) * (z - 0.5) +
	x * y * z * (x - 1.0) * (y - 1.0) * (y - 0.5) * (z - 0.5) +
	y * z * (x - 1.0) * (y - 1.0) * (z - 1.0) * (x - 0.5) * (y - 0.5) +
	x * y * z * (y - 1.0) * (z - 1.0) * (x - 0.5) * (y - 0.5) +
	x * y * z * (x - 1.0) * (y - 1.0) * (z - 1.0) * (y - 0.5) +
	Exp(x) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) -
	2.0 * Exp(x) * Exp(x) * (M_E - Exp(2.0 * x)) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) -
	4.0 * Exp(x) * (M_E - Exp(x)) * Exp(2.0 * x) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	(1.0 - Exp(x)) * Exp(x) * (M_E - Exp(2.0 * x)) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) -
	4.0 * (1.0 - Exp(x)) * Exp(x) * Exp(2.0 * x) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	4.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * Exp(2.0 * x) * (1.0 - Exp(y)) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) + (1.0 - Exp(x)) *
		(M_E - Exp(x)) * (M_E - Exp(2.0 * x)) * Exp(y) *
		(M_E - Exp(y)) * (M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) - 2.0 * (1.0 - Exp(x)) *
		(M_E - Exp(x)) * (M_E - Exp(2.0 * x)) * Exp(y) * Exp(y) *
		(M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) * (M_E - Exp(z)) *
		(M_E - Exp(2.0 * z)) - 4.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) *
		(M_E - Exp(2.0 * x)) * Exp(y) * (M_E - Exp(y)) * Exp(2.0 * y) *
		(1.0 - Exp(z)) * (M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	(1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * Exp(y) * (M_E - Exp(2.0 * y)) *
		(1.0 - Exp(z)) * (M_E - Exp(z)) * (M_E - Exp(2.0 * z)) -
	4.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * Exp(y) * Exp(2.0 * y) * (1.0 - Exp(z)) *
		(M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	4.0 * (1.0 - Exp(x)) * (M_E - Exp(x)) * (M_E - Exp(2.0 * x)) *
		(1.0 - Exp(y)) * (M_E - Exp(y)) * Exp(2.0 * y) *
		(1.0 - Exp(z)) * (M_E - Exp(z)) * (M_E - Exp(2.0 * z)) +
	4.0 * Sin(2.0 * M_PI * x) * Cos(2.0 * M_PI * y) * M_PI * M_PI *
		Cos(2.0 * M_PI * z) + (1.0 - Exp(x)) * (M_E - Exp(x)) *
		(M_E - Exp(2.0 * x)) * (1.0 - Exp(y)) * (M_E - Exp(y)) *
		(M_E - Exp(2.0 * y)) * (1.0 - Exp(z)) * (M_E - Exp(z)) *
		(M_E - Exp(2.0 * z));
    return;
}

static void
func_curl_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
/*  u = x y z (x-1.0)(y-1.0)(z-1.0)(x-0.5)(y-0.5)(z-0.5)
    v = Sin[2 Pi x] Sin[2 Pi y] Sin[2 Pi z]
    w = (1. - Exp[x])(E - Exp[x])(E - Exp[2. x]) (1. - Exp[y])(E - Exp[y])(E - Exp[2. y]) (1. - Exp[z])(E - Exp[z])(E - Exp[2. z])
    CForm[Simplify[D[w,y] - D[v,z]]]
    CForm[Simplify[D[u,z] - D[w,x]]]
    CForm[Simplify[D[v,x] - D[u,y]]] */
{
    (void)func_curl_u;

    *(value++) = -2.*Exp(2.*y)*(1. - Exp(x))*(M_E - Exp(x))*(M_E - Exp(2.*x))*
     (1. - Exp(y))*(M_E - Exp(y))*(1. - Exp(z))*(M_E - Exp(z))*
     (M_E - Exp(2.*z)) - Exp(y)*(1. - Exp(x))*(M_E - Exp(x))*
     (M_E - Exp(2.*x))*(1. - Exp(y))*(M_E - Exp(2.*y))*
     (1. - Exp(z))*(M_E - Exp(z))*(M_E - Exp(2.*z)) - 
    Exp(y)*(1. - Exp(x))*(M_E - Exp(x))*(M_E - Exp(2.*x))*
     (M_E - Exp(y))*(M_E - Exp(2.*y))*(1. - Exp(z))*(M_E - Exp(z))*
     (M_E - Exp(2.*z)) - 2*M_PI*Cos(2*M_PI*z)*Sin(2*M_PI*x)*Sin(2*M_PI*y);

    *(value++) = 2.*Exp(2.*x)*(1. - Exp(x))*(M_E - Exp(x))*(1. - Exp(y))*
     (M_E - Exp(y))*(M_E - Exp(2.*y))*(1. - Exp(z))*(M_E - Exp(z))*
     (M_E - Exp(2.*z)) + Exp(x)*(1. - Exp(x))*(M_E - Exp(2.*x))*
     (1. - Exp(y))*(M_E - Exp(y))*(M_E - Exp(2.*y))*(1. - Exp(z))*
     (M_E - Exp(z))*(M_E - Exp(2.*z)) + 
    Exp(x)*(M_E - Exp(x))*(M_E - Exp(2.*x))*(1. - Exp(y))*
     (M_E - Exp(y))*(M_E - Exp(2.*y))*(1. - Exp(z))*(M_E - Exp(z))*
     (M_E - Exp(2.*z)) + (-1. + x)*(-0.5 + x)*x*(-1. + y)*(-0.5 + y)*y*
     (-1. + z)*(-0.5 + z) + (-1. + x)*(-0.5 + x)*x*(-1. + y)*(-0.5 + y)*y*
     (-1. + z)*z + (-1. + x)*(-0.5 + x)*x*(-1. + y)*(-0.5 + y)*y*(-0.5 + z)*z;

    *(value++) =
   -((-1. + x)*(-0.5 + x)*x*(-1. + y)*(-0.5 + y)*(-1. + z)*(-0.5 + z)*z) - 
    (-1. + x)*(-0.5 + x)*x*(-1. + y)*y*(-1. + z)*(-0.5 + z)*z - 
    (-1. + x)*(-0.5 + x)*x*(-0.5 + y)*y*(-1. + z)*(-0.5 + z)*z + 
    2*M_PI*Cos(2*M_PI*x)*Sin(2*M_PI*y)*Sin(2*M_PI*z);

    return;
}


static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h, DOF *kk)
{
    int N = u_h->type->nbas;	/* number of basis functions in an element */
    GRID *g = u_h->g;
    ELEMENT *e;
    int i, j;
    FLOAT A[N][N], B[N], buffer[N];
    INT I[N];

    assert(u_h->dim == 1);

    ForAllElements(g, e) {	/* \PHGindex{ForAllElements} */
	/* compute \int \curl\phi_j \cdot \curl\phi_i making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (j = 0; j <= i; j++)
		A[j][i] = A[i][j] =
		    phgQuadCurlBasDotCurlBas(e, u_h, j, u_h, i, QUAD_DEFAULT) -
		    phgQuadBasABas(e, u_h, j, kk, u_h, i, QUAD_DEFAULT);
	}

	/* compute local matrix and RHS */
	for (i = 0; i < N; i++) {	/* loop on basis functions */
	    if (phgDofDirichletBC(u_h, e, i, func_u, 
				buffer, B + i, DOF_PROJ_CROSS)) {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer); 
	    }
	    else {
		/* right hand side: \int f * phi_i */
		B[i] = phgQuadDofDotBas(e, f_h, u_h, i, QUAD_DEFAULT);
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
	    }
	}
	phgSolverAddRHSEntries(solver, N, I, B);
    }
}

static FLOAT
estimate_error(DOF *u_h, DOF *f_h, DOF *curl_u_h, DOF *error)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *jmp1, *jmp2, *res, *div, *tmp;

    res = phgDofCopy(f_h, NULL, u_h->type, "f+k^2u_h");
    phgDofAXPY(1.0, u_h, &res);		/* res = f_h + k^2 u_h */
    //phgDofMM(MAT_OP_N, MAT_OP_N, 1, 3, 1, 1.0, kk, 1, u_h, 1.0, &res);
    jmp1 = phgQuadFaceJump(curl_u_h, DOF_PROJ_CROSS, NULL, QUAD_DEFAULT);
    jmp2 = phgQuadFaceJump(res, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    tmp = phgDofCopy(res, NULL, NULL, NULL);	/* save f_h + k^2 u_h */
    div = phgDofCurl(curl_u_h, NULL, NULL, "residual");
    phgDofAXPBY(-1.0, div, 1.0, &res);		/* residual */
    phgDofFree(&div);
    div = phgDofDivergence(tmp, NULL, NULL, "div(f+k^2 u_h)");
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
    ELEMENT *e;
    GRID *g;
    DOF *u_h, *u, *f_h, *error, *curl, *tmp, *kk;
    SOLVER *solver, *pc = NULL;
    FLOAT L2_err, Hcurl_err, a, b, thres, max_error;
    size_t mem, mem_peak;

    static char *fn = "cube1o.dat";
    static char *vtk = NULL;
    static INT pre_refines = 0;
    static INT depth = 1;	/* <=0 => uniformly refine -depth or 3 times */
    static INT mem_max = 400;	/* max. memory */
    static FLOAT tol = 1e-3;	/* convergence tolerance */

    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **) &fn);
    phgOptionsRegisterFilename("vtk_file", "VTK file", (char **) &vtk);
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("refine_depth", "Refinement depth (<=0 means "
			  "uniformly refine -depth or 3 times)", &depth);
    phgOptionsRegisterInt("mem_max", "Max memory (MB)", &mem_max);
    phgOptionsRegisterFloat("tol", "Convergence tolerance", &tol);

    phgOptionsPreset("-dof_type ND1");
    if (SOLVER_HYPRE != NULL)
	phgOptionsPreset("-hypre_pc ams");

    phgInit(&argc, &argv);
    phgOptionsShowUsed();
   if (depth == 0)
	depth = -3;

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    phgRefineAllElements(g, pre_refines);

    u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation/*DofNoAction*/);
    phgDofSetDataByValue(u_h, 0.0);
    u = phgDofNew(g, DOF_DEFAULT, 1, "u", func_u);
    f_h = phgDofNew(g, DOF_ANALYTIC, 3/*DOF_DEFAULT, 1*/, "f_h", func_f);

    kk = phgDofNew(g, DOF_P0, 1, "kk", DofInterpolation);
    phgDofSetDataByFunction(kk, func_kk);

    error = phgDofNew(g, DOF_P0, 1, "error estimates", DofNoAction);

    while (TRUE) {
	elapsed_time(g, FALSE, 0.);
	phgPrintf("\n------ %d DOF, %d elements, mesh LIF = %lg\n",
		  DofDataCountGlobal(u_h), g->nleaf_global, (double)g->lif);

	if (phgBalanceGrid(g, 1.2, 1000, NULL, 0.)) {
	    phgPrintf("------ Repartition mesh: nprocs = %d, LIF = %lg ",
			    g->nprocs, (double)g->lif);
	    elapsed_time(g, TRUE, 0.);
	}

	phgPrintf("Set up linear solver: ");
	solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
	phgPrintf("solver LIF = %lg ", (double)solver->mat->rmap->lif);
	elapsed_time(g, TRUE, 0.);

	phgPrintf("Build linear system ");
	phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
#if 0
	build_linear_system(solver, u_h, f_h, kk);
#else
	tmp = phgDofCopy(f_h, NULL, DOF_DGn[u_h->type->order], NULL);
	phgSolverAMSRemoveGradient(&tmp, kk, u_h->type->order);
	build_linear_system(solver, u_h, tmp, kk);
	phgDofFree(&tmp);
#endif
	elapsed_time(g, TRUE, 0.);

	tmp = phgDofAXPY(-1.0, kk, NULL);
	if (solver->oem_solver == SOLVER_HYPRE) {
	    phgPrintf("Build Poisson-matrices ");
	    /*phgSolverHypreAMSSetPoisson(solver, NULL, NULL);*/
	    phgSolverHypreAMSSetPoisson(solver, NULL, tmp);
	    elapsed_time(g, TRUE, 0.);
	}
	else {
	    pc = phgMat2Solver(SOLVER_AMS, solver->mat);
	    phgSolverAMSSetCoefficients(pc, NULL, tmp);
	    pc->rtol = 0.;
	    pc->maxit = 1;
	    phgSolverSetPC(solver, pc, NULL);
	}
	phgDofFree(&tmp);

	phgPrintf("Solve linear system: ");
	phgSolverSolve(solver, TRUE, u_h, NULL);
	phgPrintf("nits=%d, resid=%0.4lg ", solver->nits,
			(double)solver->residual);
	elapsed_time(g, TRUE, 0.);
	if (pc != NULL)
	    phgSolverDestroy(&pc);
	phgSolverDestroy(&solver);

	phgPrintf("Errors: ");
	tmp = phgDofCopy(u, NULL, NULL, NULL);
	phgDofAXPY(-1.0, u_h, &tmp);
	phgSolverAMSRemoveGradient(&tmp, kk, u_h->type->order);
	phgPrintf("L2=%0.3le, ", (double)(L2_err = phgDofNormL2(tmp)));
	curl = phgDofCurl(tmp, NULL, NULL, NULL);
	phgDofFree(&tmp);
	tmp = curl;
	curl = phgDofCurl(u_h, NULL, NULL, NULL);
	phgPrintf("Hcurl=%0.3le, ", (double)(Hcurl_err = phgDofNormL2(tmp)));
	phgDofFree(&tmp);
	max_error = estimate_error(u_h, f_h, curl, error);
	phgPrintf("Est=%0.3le ", (double)phgDofNormL2Vec(error));
	elapsed_time(g, TRUE, 0.);
	phgPrintf("Norms: L2=%0.3lg, Hcurl=%0.3lg\n",
			(double)phgDofNormL2(u_h), (double)phgDofNormL2(curl));
	mem = phgMemoryUsage(g, &mem_peak);
	a = mem / (1024.0 * 1024.0);
	b = mem_peak / (1024.0 * 1024.0);
	phgPrintf("Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
			(double)a, (double)b);
	phgPrintf("Peformance of the step: %0.4lgMflops\n",
			phgPerfGetMflops(g, NULL, NULL));
	if (Hcurl_err < tol || mem_peak >= 1024 * (size_t)mem_max * 1024) {
	    if (vtk != NULL) {
		const char *vtk_fn;
		phgDofFree(&f_h);
		phgRefineAllElements(g, 0);
		tmp = phgDofCopy(u, NULL, u_h->type, "u-u_h");
		phgDofAXPY(-1.0, u_h, &tmp);
		phgPrintf("Export mesh to \"%s\" ",
			  vtk_fn = phgExportVTK(g, vtk, u_h, tmp, u, kk, NULL));
		phgDofFree(&tmp);
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
	    phgDofFree(&curl);
	    break;
	}

	phgPrintf("Refine mesh ");
	if (depth <= 0) {
	    phgRefineAllElements(g, -depth);
	}
	else {
	    thres = 0.5 * max_error;
	    ForAllElements(g, e)
		if (*DofElementData(error, e->index) >= thres)
		    e->mark = depth;
	    phgRefineMarkedElements(g);
	}
	phgDofFree(&curl);
	elapsed_time(g, TRUE, 0.);
    }
    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&u);
    phgDofFree(&error);
    phgDofFree(&kk);

    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
