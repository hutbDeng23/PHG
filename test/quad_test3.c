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

/*
 * This is a sample program on using phgQuadInterfaceTetrahedron().
 *
 * The element is the standard tetrahedron and the interface is the sphere of
 * radius 1/2 centered at (0,0,0) describe by the level set function:
 *		ls(x,y,z) = x^2 + y^2 + z^2 - 0.5^2
 * and the integrand is defined by:
 * 		u(x,y,z) = (x + y + z)^2
 *
 * $Id: quad_test3.c,v 1.19 2021/02/12 01:29:06 zlb Exp $ */

#include "phg.h"
#include <math.h>

#ifndef COMPUTE_REFERENCE_VALUE
# define COMPUTE_REFERENCE_VALUE	0
#endif

#if COMPUTE_REFERENCE_VALUE
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Compute a reference value on the fly, using the following transformation
 * (r,s,t)->(x,y,z) in the tetrahedron:
 * 	[x, y, z] = [t*r, t*s, t*(1-r-s)]
 * or:
 * 	[r, s, t] = [x/(x+y+z), y/(x+y+z), x+y+z]
 * with r in [0,1], s [0,1-r] and t in [0,1].
 *
 * r, s and 1-r-s correspond to the barycentric coordinates of the triangle
 * (1,0,0), (0,1,0) and (0,0,1).
 *
 * t corresponds to the line segment linking the origin (0,0,0) to the
 * point on the triangle. The intersection point of the line segment with
 * the interface is given by:
 *	t0(r,s) = 1/2/sqrt(r^2+s^2+(1-r-s)^2)
 *
 * The integrand u(x,y,z) = (x+y+z)^2 = t^2.
 *
 * (1) I^- is computed as
 *	  \int_{r=0}^{1} \int_{s=0}^{1-r} \int_{t=0}^{t0(r,s)} det(J) dt dr ds
 *     Mathematica code: N[Re[Integrate[Integrate[Integrate[t^2*t^2, {t,0,1/2/Sqrt[r^2+s^2+(1-r-s)^2]}], {s,0,1-r}], {r,0,1}]], 40]
 *
 * (2) I^0 is computed as
 *	  \int_{r=0}^{1} \int_{s=0}^{1-r} S(r,s) u(r,s) dr ds
 *     Mathematica code: N[Re[Integrate[Integrate[1/4*(2r^2+2(r-1)s+2s^2-2r+1)^(-3/2)*1/4/(r*r+s*s+(1-r-s)^2), {s,0,1-r}], {r,0,1}]], 40]
 *
 ***************************** SageMath code for computing det(J) and S(r,s):
var('r,s,t')
xyz = [t*r, t*s, t*(1-r-s)]
print "det(J) =", factor(det(jacobian(xyz, (r,s,t))))

def cross(u,v):
    return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
def D(u,r):
    return [derivative(u[0],r), derivative(u[1],r), derivative(u[2],r)]
t = 1/2/sqrt(r^2+s^2+(1-r-s)^2)
xyz = [t*r, t*s, t*(1-r-s)]
print "S(r,s) =", factor(norm(vector(cross(D(xyz,r), D(xyz,s)))))
 *****************************
 *	det(J) = t^2
 *	S(r,s) = 1/4*(2*r^2+2*(r-1)*s+2*s^2-2*r+1)^(-3/2)
 */

static QUAD *quad;
static FLOAT r;

static void
func_t(FLOAT t, FLOAT *res, const FLOAT wgt, void *ctx)
{
    *res = (t * t) * (t * t);	/* Note J = integrand = t^2 */
}

static void
func_s(FLOAT s, FLOAT *res, const FLOAT wgt, void *ctx)
{
    INT type = *(INT *)ctx;
    FLOAT t2;

    t2 = 0.25 / (r*r+s*s+Pow(1.-r-s,2.));	/* t^2 */
    if (type != 0)
	phgQuad1D(func_t, 1, 0.0, Sqrt(t2), quad, res, ctx);
    else
	*res = 1./4.*Pow(2.*r*r+2.*(r-1.)*s+2.*s*s-2.*r+1., -1.5) * t2;
}

static void
func_r(FLOAT r0, FLOAT *res, const FLOAT wgt, void *ctx)
{
    r = r0;
    phgQuad1D(func_s, 1, 0.0, 1.0 - r, quad, res, ctx);
}

static FLOAT
get_ref(INT type)
{
    FLOAT ref;
#if SIZEOF_PHG_FLOAT == 16
    quad = phgQuadGetQuad1D(99);
#else	/* SIZEOF_PHG_FLOAT == 16 */
    quad = phgQuadGetQuad1D(59);
#endif	/* SIZEOF_PHG_FLOAT == 16 */
    phgQuad1D(func_r, 1, 0.0, 1.0, quad, &ref, &type);
    if (type > 0)
	ref = _F(0.1) - ref;
    return ref;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#endif	/* COMPUTE_REFERENCE_VALUE */

static void
ls(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
/* the level set function */
{
    *value = x * x + y * y + z * z - .5 * .5;
}

static void
ls_grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *grad)
/* the gradient of the level set function */
{
    grad[0] = x + x;
    grad[1] = y + y;
    grad[2] = z + z;
}

static void
u(FLOAT x, FLOAT y, FLOAT z, FLOAT *f)
/* the integrand */
{
    *f = (x + y + z) * (x + y + z);
}

int
main(int argc, char *argv[])
{
    INT order = 15, type = -1, repeat = 1;
    FLOAT tet[4][3] =  {{1., 0., 0.},	/* vertex 0 */
			{0., 1., 0.},	/* vertex 1 */
			{0., 0., 1.},	/* vertex 2 */
			{0., 0., 0.}};	/* vertex 3 */
    FLOAT res = 0., ref = 1.0, *rule, **rules[3] = {NULL, NULL, NULL};
    int i, n = 0;
    double t;
    char flt1[48], flt2[48], flt3[48];

    phgOptionsRegisterInt("-order", "Gaussian quadrature order", &order);
    phgOptionsRegisterInt("-type",  "Quadrature type (-1, 0, 1)", &type);
    phgOptionsRegisterInt("-repeat", "Times to repeat the test", &repeat);

    /* phgInitSetComm(MPI_COMM_SELF); */
    phgInit(&argc, &argv);

#if COMPUTE_REFERENCE_VALUE
    ref = get_ref(type);
#else	/* COMPUTE_REFERENCE_VALUE */
    if (type == 0) {
	ref = _F(0.2231747704246810387019576057274844651312);
    }
    else {
	type = type < 0 ? -1 : 1;
	ref = _F(0.02231747704246810387019576057274844651312);
	if (type > 0)
	    ref = _F(0.1) - ref;
    }
#endif	/* COMPUTE_REFERENCE_VALUE */

    rules[type + 1] = &rule;

    t = phgGetTime(NULL);
    for (i = 0; i < repeat; i++) {
#if 1
	/* Test code, should compute the integral on the whole tetra,
	 * which is the sum of the results for '-type -1' and 'type 1'.  */
	phgQuadInterfaceTetrahedron(
		NULL,		/* the level set function */
		2,		/* polynomial order of the level set function */
		NULL,		/* the gradient of the level set function */
		tet,		/* coordinates of the vertices of the tetra */
		order,		/* order of the 1D Gaussian quadrature */
		rules[0], rules[1], rules[2]);
	n = phgQuadInterfaceRuleApply(u, 1, PROJ_NONE, rule, &res);
	phgFree(rule);
	printf("Integral on the whole tetrahedron: %lg\n", (double)res);
#endif
	phgQuadInterfaceTetrahedron(
		ls,		/* the level set function */
		2,		/* polynomial order of the level set function */
		ls_grad,	/* the gradient of the level set function */
		tet,		/* coordinates of the vertices of the tetra */
		order,		/* order of the 1D Gaussian quadrature */
		rules[0], rules[1], rules[2]);
	n = phgQuadInterfaceRuleApply(u, 1, PROJ_NONE, rule, &res);
	phgFree(rule);
    }
    t = phgGetTime(NULL) - t;

#if FT_PHG == FT___FLOAT128
	quadmath_snprintf(flt1, sizeof(flt1), "%0.34Qe", ref);
	quadmath_snprintf(flt2, sizeof(flt2), "%0.34Qe", res);
	quadmath_snprintf(flt3, sizeof(flt3), "%0.34Qe", Fabs(1. - res / ref));
#else
	snprintf(flt1, sizeof(flt1), "%0.16le", (double)ref);
	snprintf(flt2, sizeof(flt2), "%0.16le", (double)res);
	snprintf(flt3, sizeof(flt3), "%0.16le", (double)Fabs(1. - res / ref));
#endif
    printf("Domain:     %s\n"
	   "Quad order: %d\n"
	   "Reference:  %s\n"
	   "Computed:   %s\n"
	   "Rel error:  %s\n"
	   "Func evals: %d\n"
	   "Time:       %lg\n",
		type == 0 ? "surface" : (type < 0 ? "Omega-" : "Omega+"),
		order, flt1, flt2, flt3, n, t);

    phgFinalize();

    return 0;
}
