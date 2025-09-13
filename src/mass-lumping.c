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

/* $Id: mass-lumping.c,v 1.100 2022/09/02 14:02:44 zlb Exp $ */

/*-----------------------------------------------------------------------------
 * This file implements the finite element bases with mass lumping (MLn).
 *
 * Note: the FE bases consist of linear combinations of elementary bases.
 * For elementary bases, we use hierarchical homogeneous monomials as follows:
 * 	For vertex nodes:
 * 		\lambda^O1
 * 	For edge nodes:
 * 		b\cdot\lambda^p, \lambda\in\R^2, p\in\N^2, |p| = O1-2
 * 		(equivalent to \lambda^p, p\in{\N^+}^2, |p| = O1)
 * 	For face nodes:
 * 		b\cdot\lambda^p, \lambda\in\R^3, p\in\N^3, |p| = O2-3
 * 		(equivalent to \lambda^p, p\in{\N^+}^3, |p| = O2)
 * 	For interior nodes:
 * 		b\cdot\lambda^p, \lambda\in\R^4, p\in\N^4, |p| = O3-4
 * 		(equivalent to \lambda^p, p\in{\N^+}^4, |p| = O3)
 *---------------------------------------------------------------------------*/

#include "phg.h"
#include "phg/quad-permu.h"
#include <string.h>		/* memset() */
#include <math.h>		/* pow() */
#include <stdlib.h>		/* qsort() */

static GTYPE orbit_type[] =
 {VERTEX, EDGE, EDGE, FACE, FACE, FACE, VOLUME, VOLUME, VOLUME, VOLUME, VOLUME};
typedef enum {S1, S2, S11, S3, S21, S111, S4, S31, S22, S211, S1111} ORBIT_TYPE;
static char orbit_size[] =
	     {1,  1,  2,   1,  3,   6,    1,  4,   6,   12,   24};

static const FLOAT *
get_node(DOF_TYPE *type, int n)
/* returns the barycentric coordinates of the 'n'-th node. */
{
    static FLOAT lambda[4];
    int i;

    if (n < NVert * type->np_vert) {		/* Vertex */
	assert(type->np_vert == 1);
	lambda[0] = lambda[1] = lambda[2] = lambda[3] = 0.0;
	lambda[n] = type->points[0];
	return lambda;
    }
    n -= NVert * type->np_vert;

    if (n < NEdge * type->np_edge) {		/* Edge */
	i = n / type->np_edge;	/* edge no. */
	n %= type->np_edge;
	lambda[0] = lambda[1] = lambda[2] = lambda[3] = 0.0;
	n = 1 * type->np_vert + 2 * n;
	lambda[GetEdgeVertex(i,0)] = type->points[n];
	lambda[GetEdgeVertex(i,1)] = type->points[n + 1];
	return lambda;
    }
    n -= NEdge * type->np_edge;

    if (n < NFace * type->np_face) {		/* Face */
	i = n / type->np_face;	/* face no. */
	n %= type->np_face;
	lambda[i] = 0.0;
	n = 1 * type->np_vert + 2 * type->np_edge + 3 * n;
	lambda[GetFaceVertex(i,0)] = type->points[n];
	lambda[GetFaceVertex(i,1)] = type->points[n + 1];
	lambda[GetFaceVertex(i,2)] = type->points[n + 2];
	return lambda;
    }
    n -= NFace * type->np_face;

    /* interior */
    assert(n >= 0 && n < type->np_elem);
    return type->points + 1 * type->np_vert
			+ 2 * type->np_edge
			+ 3 * type->np_face
			+ 4 * n;
}

static void
get_monomials(DOF_TYPE *type, const FLOAT *lambda, FLOAT *vals, BOOLEAN grad)
/* evaluates the value or the gradient of all the monomial at barycentric
 * coordinates lambda[] and return them in 'vals' */
{
    int i, n, order, v0, v1, v2;
    BYTE *po;
    FLOAT p0, p1, p2, p3, power[type->order + 1][Dim + 1];

    /* pre-compute and store powers of lambda[] */
    for (n = 0; n < Dim + 1; n++) {
	power[0][n] = 1.0;
	p0 = lambda[n];
	for (i = 0; i < type->order; i++)
	    power[i + 1][n] = power[i][n] * p0;
    }

    if (type->mass_lumping->monobas != NULL) {
	const BYTE (*pb)[Dim + 1] = (void *)type->mass_lumping->monobas;
	for (i = 0; i < type->nbas; i++) {
	    p0 = power[pb[i][0]][0];
	    p1 = power[pb[i][1]][1];
	    p2 = power[pb[i][2]][2];
	    p3 = power[pb[i][3]][3];

	    if (!grad) {
		/* compute function values */
		*(vals++) = p0 * p1 * p2 * p3;
		continue;
	    }

	    /* compute gradients */
	    if ((n = pb[i][0]) > 0)
		*(vals++) = n * p1 * p2 * p3 * power[n - 1][0];
	    else
		*(vals++) = 0.0;

	    if ((n = pb[i][1]) > 0)
		*(vals++) = n * p0 * p2 * p3 * power[n - 1][1];
	    else
		*(vals++) = 0.0;

	    if ((n = pb[i][2]) > 0)
		*(vals++) = n * p0 * p1 * p3 * power[n - 1][2];
	    else
		*(vals++) = 0.0;

	    if ((n = pb[i][3]) > 0)
		*(vals++) = n * p0 * p1 * p2 * power[n - 1][3];
	    else
		*(vals++) = 0.0;
	}
	return;
    }

    order = type->mass_lumping->order0;
    assert(type->np_vert == 1 || order == 0);
    for (i = 0; i < NVert && order > 0; i++) {		/* Vertices */
	if (!grad) {
	    *vals = power[order][i];
	    vals++;
	}
	else {
	    vals[0] = vals[1] = vals[2] = vals[3] = 0.;
	    if (order >= 1)
		vals[i] = order * power[order - 1][i];
	    vals += 4;
	}
    }

    order = type->mass_lumping->order1;
    for (i = 0; i < NEdge; i++) {		/* Edges */
	for (n = 1; n <= type->np_edge; n++) {
	    if (!grad) {
		*vals = power[n][GetEdgeVertex(i,0)] *
			power[order - n][GetEdgeVertex(i,1)];
		vals++;
	    }
	    else {
		v0 = GetEdgeVertex(i,0);
		v1 = GetEdgeVertex(i,1);
		vals[0] = vals[1] = vals[2] = vals[3] = 0.;
		vals[v0] = n * power[n - 1][v0] *
			   power[order - n][v1];
		vals[v1] = (order - n) * power[n][v0] *
			   power[order - n - 1][v1];
		vals += 4;
	    }
	}
    }

    order = type->mass_lumping->order2;
    for (i = 0; i < NFace; i++) {		/* Faces */
	for (n = 0; n < type->np_face; n++) {
	    po = type->mass_lumping->face_pow + 2 * n;
	    v0 = GetFaceVertex(i,0);
	    v1 = GetFaceVertex(i,1);
	    v2 = GetFaceVertex(i,2);
	    p0 = power[po[0]][v0];
	    p1 = power[po[1]][v1];
	    p2 = power[order - po[0] - po[1]][v2];
	    if (!grad) {
		*vals = p0 * p1 * p2;
		vals++;
	    }
	    else {
		vals[i] = 0.0;
		vals[v0] = po[0] * p1 * p2 * power[po[0] - 1][v0];
		vals[v1] = po[1] * p0 * p2 * power[po[1] - 1][v1];
		vals[v2] = (order - po[0] - po[1]) * p0 * p1 *
			   power[order - po[0] - po[1] - 1][v2];
		vals += 4;
	    }
	}
    }

    /* interior */
    order = type->mass_lumping->order3;
    for (n = 0; n < type->np_elem; n++) {
	po = type->mass_lumping->elem_pow + 3 * n;
	p0 = power[po[0]][0];
	p1 = power[po[1]][1];
	p2 = power[po[2]][2];
	p3 = power[order - po[0] - po[1] - po[2]][3];
	if (!grad) {
	    *vals = p0 * p1 * p2 * p3;
	    vals++;
	}
	else {
	    vals[0] = po[0] * p1 * p2 * p3 * power[po[0] - 1][0];
	    vals[1] = po[1] * p0 * p2 * p3 * power[po[1] - 1][1];
	    vals[2] = po[2] * p0 * p1 * p3 * power[po[2] - 1][2];
	    vals[3] = (order - po[0] - po[1] - po[2]) * p0 * p1 * p2 *
		      power[order - po[0] - po[1] - po[2] - 1][3];
	    vals += 4;
	}
    }
}

#if TEST_TYPES == 3
/*----------------------------------------------------------------------------*/
/* The warp & bend nodes by:
 * T. Warburton, An explicit construction of interpolation nodes on the simplex,
 * J Eng Math (2006) 56:247-262 */

/* alpha[k] are for order k, with alpha[k][0/1] for triangle/tetrahedron */
static FLOAT alpha[][2] = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {1.4152, 0.0},
	{0.1001, 0.1002}, {0.2751, 1.1332}, {0.9808, 1.5608}, {1.0999, 1.3413},
	{1.2832, 1.2577}, {1.3648, 1.1603}, {1.4773, 1.0153}, {1.4959, 0.6080},
	{1.5743, 0.4523}, {1.5770, 0.8856}, {1.6223, 0.8717}, {1.6258, 0.9655}};

static FLOAT
warp(int order, FLOAT x, const FLOAT GLL[])
/* The 1D warping function. Note x in [-1, 1] */
{
    int i, j;
    FLOAT res = 0., d;

    assert(order > 0);
#define X0(i)	((2 * (i) - order) / (FLOAT)order)
    for (i = 0; i <= order; i++) {
	d = (2.0 * GLL[i] - 1.0) - X0(i);
	for (j = 0; j <= order; j++)
	    if (j != i)
		d *= (x - X0(j)) / (X0(i) - X0(j));
	res += d;
    }

    return res;
}

#define B(L0, L1)	((L0) == 0. ? 1. :  2.* (L1) / (2. * (L1) + (L0)))

static void
warp_bend2(int order, FLOAT lam[3], const FLOAT GLL[], FLOAT alpha)
/* 2D bend function */
{
    FLOAT a2 = alpha * alpha;
    FLOAT g0, g1, g2;

    g0 = 1.0 / (lam[0] + lam[1] + lam[2]);
    lam[0] *= g0;
    lam[1] *= g0;
    lam[2] *= g0;

    assert(order > 0);
#define G(L0, L1, L2) warp(order, (L2)-(L1), GLL) * B(L0, L1) * B(L0, L2)
    g0 = (1. + a2 * lam[0] * lam[0]) * G(lam[0], lam[1], lam[2]) * 0.5;
    g1 = (1. + a2 * lam[1] * lam[1]) * G(lam[1], lam[2], lam[0]) * 0.5;
    g2 = (1. + a2 * lam[2] * lam[2]) * G(lam[2], lam[0], lam[1]) * 0.5;

    lam[1] -= g0;
    lam[2] += g0;

    lam[2] -= g1;
    lam[0] += g1;

    lam[0] -= g2;
    lam[1] += g2;
}

static void
warp_bend3(int order, FLOAT lam[4], const FLOAT GLL[], FLOAT alpha)
/* 3D warp & bend function */
{
    FLOAT a2 = alpha * alpha;
    FLOAT g0, g1, g2, res[4];
    int i, i0, i1, i2;

    assert(order > 0);

    g0 = 1.0 / (lam[0] + lam[1] + lam[2] + lam[3]);
    lam[0] *= g0;
    lam[1] *= g0;
    lam[2] *= g0;
    lam[3] *= g0;

#undef G
#define G(L0,L1,L2) warp(order,(L2)-(L1), GLL) * \
	B(lam[i],L0) * B(lam[i],L1) * B(lam[i],L2)
    memcpy(res, lam, sizeof(res));
    for (i = 0; i < NFace; i++) {
	i0 = GetFaceVertex(i, 0);
	i1 = GetFaceVertex(i, 1);
	i2 = GetFaceVertex(i, 2);
	g0 = (1. + a2 * lam[i] * lam[i]) * G(lam[i0], lam[i1], lam[i2]) * .5;
	g1 = (1. + a2 * lam[i] * lam[i]) * G(lam[i1], lam[i2], lam[i0]) * .5;
	g2 = (1. + a2 * lam[i] * lam[i]) * G(lam[i2], lam[i0], lam[i1]) * .5;

	res[i1] -= g0;
	res[i2] += g0;

	res[i2] -= g1;
	res[i0] += g1;

	res[i0] -= g2;
	res[i1] += g2;
    }

    memcpy(lam, res, sizeof(res));
}
/*----------------------------------------------------------------------------*/
#endif /* TEST_TYPES == 3 */

static inline int
fltcmp(const FLOAT *p0, const FLOAT *p1, int n)
{
    int i;

    for (i = 0; i < n; i++, p0++, p1++)
	if (Fabs(*p0 - *p1) > FLOAT_EPSILON * 10.0 * (Fabs(*p0) + Fabs(*p1)))
	    return 1;
    return 0;
}

static FLOAT *work = NULL;

static int
comp_2float(const void *q1, const void *q2)
/* compares two pairs of FLOATS (for sorting the edge nodes).
 * Note: the nodes are sorted into decreasing order to match DOF_P[1234] */
{
    const FLOAT *p1 = work + 2 * (*(const int *)q1),
		*p2 = work + 2 * (*(const int *)q2);

    return *p1 == *p2 ? 0 : (*p1 > *p2 ? -1 : 1);
}

static void
init(DOF_TYPE *type)
/* this defines the 'Initialize()' member function of DOF_TYPE */
{
    int i, n;
    FLOAT d, *p, *pw, *pw0;

#if USE_OMP
    assert(phgThreadId == 0);
#endif	/* USE_OMP */

    if (type->mass_lumping->weights0[0] < 0.0)
	type->mass_lumping->weights0 = NULL;

    if (type->mass_lumping->monobas[0] < 0)
	type->mass_lumping->monobas = NULL;

    /* compute np_{vert,edge,face,elem}, nbas, order, and grad_type */
    type->np_vert = type->np_edge = type->np_face = type->np_elem = 0;
    for (i = 0; i < type->mass_lumping->norbit; i++) {
	switch (orbit_type[type->mass_lumping->orbits[i]]) {
	    case VERTEX:
		type->np_vert += orbit_size[type->mass_lumping->orbits[i]];
		break;
	    case EDGE:
		type->np_edge += orbit_size[type->mass_lumping->orbits[i]];
		break;
	    case FACE:
		type->np_face += orbit_size[type->mass_lumping->orbits[i]];
		break;
	    case VOLUME:
		type->np_elem += orbit_size[type->mass_lumping->orbits[i]];
		break;
	}
    }

    assert(type->mass_lumping->npoint ==
		type->np_vert + 2 * type->np_edge + 3 * type->np_face
			+ 4 * type->np_elem);

    if (type->mass_lumping->weights0 != NULL)
	assert(type->mass_lumping->nweights0 ==
	       type->np_vert + type->np_edge + type->np_face + type->np_elem);

    type->nbas = 4 * type->np_vert +
		 6 * type->np_edge +
		 4 * type->np_face +
		 1 * type->np_elem;

    if (type->mass_lumping->order3 > 0) {
	type->order = type->mass_lumping->order3;
    }
    else {
	const char (*pb)[4] = (void *)type->mass_lumping->monobas;
	assert(type->mass_lumping->monobas != NULL);
	assert(type->nbas * 4 == type->mass_lumping->nmonobas);
	type->order = type->mass_lumping->order0;
	for (i = 0; i < type->nbas; i++)
	    if ((n = pb[i][0] + pb[i][1] + pb[i][2] + pb[i][3]) > type->order)
		type->order = n;
    }

    if (type->grad_type == NULL)
	type->grad_type = DOF_DGn[type->order - 1];

#if TEST_TYPES
    if (!memcmp(type->name, "TST", 3)) {
	/* Bend the nodes into non-uniform distribution while keeping symmetry.
	 * Change TRUE => FALSE below for no bending (uniform lattice) */
	if (TRUE && type->order > 2) {
	    extern FLOAT *_phg_GLP[];	/* GLL points scaled to [0,1] */
	    FLOAT tmp[type->order + 1];

	    /* sort list of GLL points */
	    memcpy(tmp, _phg_GLP[type->order], sizeof(tmp));
	    qsort(tmp, type->order + 1, sizeof(*tmp), phgCompFLOAT);
	    p = type->points;
#if TEST_TYPES == 1
	    /* Scaled-GLL points */
	    n = type->np_vert + 2 * type->np_edge + 3 * type->np_face +
		4 * type->np_elem;
	    for (i = 0; i < n; i++) {
		assert(p[i] >= 0 && p[i] <= type->order);
		p[i] = tmp[(int)(p[i] + 0.5)];
	    }
#elif TEST_TYPES == 2
	    /* Lobatto grid over the tetrahedron (LTT):
	     * 	H. Luo and C. Pozrikidis,
	     * 	A Lobatto interpolation grid in the tetrahedron,
	     * 	IMA Journal of Applied Mathematics (2006) 71, 298-313 */
	    p += type->np_vert;
	    /* edge DOFs (GLL points) */
	    for (i = 0; i < type->np_edge; i++, p += 2) {
		p[0] = tmp[(int)(p[0] + 0.5)];
		p[1] = tmp[(int)(p[1] + 0.5)];
	    }
	    /* face DOFs */
	    for (i = 0; i < type->np_face; i++, p += 3) {
		FLOAT v0 = tmp[(int)(p[0] + 0.5)];
		FLOAT v1 = tmp[(int)(p[1] + 0.5)];
		FLOAT v2 = tmp[(int)(p[2] + 0.5)];
		p[0] = (1. + 2. * v0 - v1 - v2) / _F(3.);
		p[1] = (1. + 2. * v1 - v2 - v0) / _F(3.);
		p[2] = (1. + 2. * v2 - v0 - v1) / _F(3.);
	    }
	    /* volume DOFs */
	    for (i = 0; i < type->np_elem; i++, p += 4) {
		FLOAT v0 = tmp[(int)(p[0] + 0.5)];
		FLOAT v1 = tmp[(int)(p[1] + 0.5)];
		FLOAT v2 = tmp[(int)(p[2] + 0.5)];
		FLOAT v3 = tmp[(int)(p[3] + 0.5)];
		p[0] = (1. + 3. * v0 - v1 - v2 - v3) * _F(.25);
		p[1] = (1. + 3. * v1 - v2 - v3 - v0) * _F(.25);
		p[2] = (1. + 3. * v2 - v3 - v0 - v1) * _F(.25);
		p[3] = (1. + 3. * v3 - v0 - v1 - v2) * _F(.25);
	    }
#else	/* TEST_TYPES == 3 */
	    /* Warp & bend [T. Warburton, 2006] */
	    p += type->np_vert;
	    /* edge DOFs (GLL points) */
	    for (i = 0; i < type->np_edge; i++, p += 2) {
		p[0] = tmp[(int)(p[0] + 0.5)];
		p[1] = tmp[(int)(p[1] + 0.5)];
	    }
	    /* face DOFs */
	    for (i = 0; i < type->np_face; i++, p += 3)
		warp_bend2(type->order, p, tmp, alpha[type->order][1]);
	    /* volume DOFs */
#warning warp_bend3() is unverified
	    for (i = 0; i < type->np_elem; i++, p += 4)
		warp_bend3(type->order, p, tmp, alpha[type->order][1]);
#endif	/* TEST_TYPES == */
	}

	/* The coordinates in points[] of the "TST*" types may need to be
	 * rescaled such that the sum of the coordinates of each node == 1 */
	p = type->points;
	for (i = 0; i < type->np_vert; i++) {
	    p[i] = 1.0;
	}
	p += type->np_vert;
	for (i = 0; i < type->np_edge; i++) {
	    d = 1.0 / (p[i * 2] + p[i * 2 + 1]);
	    p[i * 2 + 0] *= d;
	    p[i * 2 + 1] *= d;
	}
	p += 2 * type->np_edge;
	for (i = 0; i < type->np_face; i++) {
	    d = 1.0 / (p[i * 3] + p[i * 3 + 1] + p[i * 3 + 2]);
	    p[i * 3 + 0] *= d;
	    p[i * 3 + 1] *= d;
	    p[i * 3 + 2] *= d;
	}
	p += 3 * type->np_face;
	for (i = 0; i < type->np_elem; i++) {
	    d = 1.0 / (p[i * 4] + p[i * 4 + 1] + p[i * 4 + 2] + p[i * 4 + 3]);
	    p[i * 4 + 0] *= d;
	    p[i * 4 + 1] *= d;
	    p[i * 4 + 2] *= d;
	    p[i * 4 + 3] *= d;
	}

#if 0
#warning
	/* print orbits for debugging */
	if (!strcmp(type->name, "TST5")) {
	    p = type->points;
	    for (i = 0; i < type->mass_lumping->norbit; i++) {
		/* Note: the assertions must match the macros in quad-permu.h */
		ORBIT_TYPE o = type->mass_lumping->orbits[i];
		switch (o) {
		    case S1:
			printf("Cons1(%g)\n", (double)p[0]);
			break;
		    case S2:
			assert(!fltcmp(p + 0, p + 1, 1));
			printf("Cons2(%g)\n", (double)p[0]);
			break;
		    case S11:
			printf("Cons11(%g)\n", (double)p[0]);
			break;
		    case S3:
			assert(!fltcmp(p + 0, p + 1, 1) &&
			       !fltcmp(p + 0, p + 2, 1));
			printf("Cons3(%g)\n", (double)p[0]);
			break;
		    case S21:
			assert(!fltcmp(p + 1, p + 2, 1));
			printf("Cons21(%g)\n", (double)p[1]);
			break;
		    case S111:
			printf("Cons111(%g,%g)\n", (double)p[0], (double)p[1]);
			break;
		    case S4:
			assert(!fltcmp(p + 0, p + 1, 1) &&
			       !fltcmp(p + 0, p + 2, 1) &&
			       !fltcmp(p + 0, p + 3, 1));
			printf("Perm4(%g)\n", (double)p[0]);
			break;
		    case S31:
			assert(!fltcmp(p + 1, p + 2, 1) &&
			       !fltcmp(p + 1, p + 3, 1));
			printf("Perm31(%g)\n", (double)p[1]);
			break;
		    case S22:
			assert(!fltcmp(p + 0, p + 1, 1) &&
			       !fltcmp(p + 2, p + 3, 1));
			printf("Perm22(%g)\n", (double)p[0]);
			break;
		    case S211:
			assert(!fltcmp(p + 0, p + 1, 1));
			printf("Perm211(%g,%g)\n", (double)p[0], (double)p[1]);
			break;
		    case S1111:
			printf("Perm1111(%g,%g,%g)\n", (double)p[0],
				(double)p[1], (double)p[2]);
			break;
		}
		p += orbit_size[o - S1] * (o==S1 ? 1 : o<S3 ? 2 : o<S4 ? 3 : 4);
	    }
	}
#endif

    }
#endif	/* TEST_TYPES */

    /* sort the edge nodes (for phgDofMapEdgeData to work correctly) */
    p = type->points + type->np_vert;
    if (type->np_edge > 1) {
	int ind[type->np_edge];
	FLOAT tmp[2 * type->np_edge];

	for (i = 0; i < type->np_edge; i++)
	    ind[i] = i;
	work = p;
	qsort(ind, type->np_edge, sizeof(ind[0]), comp_2float);
	/* reorder the points and rescale the coordinates */
	memcpy(tmp, p, 2 * type->np_edge * sizeof(*tmp));
	for (i = 0; i < type->np_edge; i++) {
	    d = 1.0 / (tmp[2 * ind[i] + 0] + tmp[2 * ind[i] + 1]);
	    p[2 * i + 0] = tmp[2 * ind[i] + 0] * d;
	    p[2 * i + 1] = tmp[2 * ind[i] + 1] * d;
	}
	if (type->mass_lumping->weights0 != NULL) {
	    /* reorder weights0[] accordingly */
	    work = type->mass_lumping->weights0 + type->np_vert;
	    memcpy(tmp, work, type->np_edge * sizeof(*tmp));
	    for (i = 0; i < type->np_edge; i++)
		work[i] = tmp[ind[i]];
	}
	/* reset work */
	work = NULL;
    }

    if (type->mass_lumping->weights0 == NULL)
	return;		/* non mass-lumping (ordinary) element */

    /* set up the weights[] array */
    pw0 = type->mass_lumping->weights0;

    pw = type->mass_lumping->weights = phgAlloc(type->nbas * sizeof(*pw));
    FreeAtExitNoCheck(type->mass_lumping->weights);

    if (type->np_vert > 0) {
	for (i = 0; i < NVert; i++, pw += type->np_vert)
	    memcpy(pw, pw0, type->np_vert * sizeof(*pw));
	pw0 += type->np_vert;
    }

    if (type->np_edge > 0) {
	for (i = 0; i < NEdge; i++, pw += type->np_edge)
	    memcpy(pw, pw0, type->np_edge * sizeof(*pw));
	pw0 += type->np_edge;
    }

    if (type->np_face > 0) {
	for (i = 0; i < NFace; i++, pw += type->np_face)
	    memcpy(pw, pw0, type->np_face * sizeof(*pw));
	pw0 += type->np_face;
    }

    if (type->np_elem > 0)
	memcpy(pw, pw0, type->np_elem * sizeof(*pw));
}

static void
init_buffers(DOF_TYPE *type)
{
    int i, j, k, order;
    BYTE *p;
    FLOAT *A;

    if (type->mass_lumping->coeff != NULL)
	return;

    /* generate orders of monomials on a face */
    if ((order = type->mass_lumping->order2) > 0) {
	type->mass_lumping->face_pow = p =
			phgAlloc(2 * type->np_face * sizeof(*p));
	FreeAtExitNoCheck(type->mass_lumping->face_pow);
	for (i = 1; i < order; i++) {
	    for (j = 1; j < order - i; j++) {
		*(p++) = i;
		*(p++) = j;
	    }
	}
	assert((int)(p - type->mass_lumping->face_pow) == 2 * type->np_face);
    }

    /* generate orders of monomials in the element */
    if ((order = type->mass_lumping->order3) > 0) {
	type->mass_lumping->elem_pow = p =
				phgAlloc(3 * type->np_elem * sizeof(*p));
	FreeAtExitNoCheck(type->mass_lumping->elem_pow);
	if (type->fe_space == FE_H1) {
	    for (i = 1; i < order; i++) {
		for (j = 1; j < order - i; j++) {
		    for (k = 1; k < order - i - j; k++) {
			*(p++) = i;
			*(p++) = j;
			*(p++) = k;
		    }
		}
	    }
	}
	else {
	    for (i = 0; i <= order; i++) {
		for (j = 0; j <= order - i; j++) {
		    for (k = 0; k <= order - i - j; k++) {
			*(p++) = i;
			*(p++) = j;
			*(p++) = k;
		    }
		}
	    }
	}
	assert((int)(p - type->mass_lumping->elem_pow) == 3 * type->np_elem);
    }

    /* buffers for caching monomials */
    type->mass_lumping->bas_lambda =
	phgCalloc(phgMaxThreads * (Dim + 1), sizeof(FLOAT));
    FreeAtExitNoCheck(type->mass_lumping->bas_lambda);
    type->mass_lumping->bas_values =
	phgAlloc(phgMaxThreads * type->nbas * sizeof(FLOAT));
    FreeAtExitNoCheck(type->mass_lumping->bas_values);

    type->mass_lumping->grad_lambda =
	phgCalloc(phgMaxThreads * (Dim + 1), sizeof(FLOAT));
    FreeAtExitNoCheck(type->mass_lumping->grad_lambda);
    type->mass_lumping->grad_values =
	phgAlloc((Dim + 1) * phgMaxThreads * type->nbas * sizeof(FLOAT));
    FreeAtExitNoCheck(type->mass_lumping->grad_values);

    /* compute the coeficient matrix of the FE bases. */
    type->mass_lumping->coeff =
		phgCalloc(type->nbas * type->nbas, sizeof(FLOAT));
    FreeAtExitNoCheck(type->mass_lumping->coeff);

    A = phgAlloc(type->nbas * type->nbas * sizeof(FLOAT));
    for (i = 0; i < type->nbas; i++) {
	type->mass_lumping->coeff[i + i * type->nbas] = 1.0;
	/* Note: a(i,j) = \psi_j(x_i) */
	get_monomials(type, get_node(type, i), A + i * type->nbas, FALSE);
    }
    phgSolverDenseSolver(type->nbas, type->nbas, A, type->mass_lumping->coeff);
    phgFree(A);
}

static const FLOAT *
ML_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT **buffers = NULL;
#if USE_OMP
#pragma omp threadprivate(buffers)
#endif  /* USE_OMP */
    int i;
    FLOAT *p, *p0, *bas_lambda, *bas_values, v;
    DOF_TYPE *type = dof->type;
    
    assert(type != NULL && type->mass_lumping != NULL);

    if (buffers == NULL) {
	buffers = phgCalloc(ORDER_MAX, sizeof(*buffers));
	FreeAtExitNoCheck(buffers);
    }

    i = type->mass_lumping->order0;
    assert(i <= ORDER_MAX);
    if (buffers[i] == NULL) {
	buffers[i] = phgCalloc(type->nbas, sizeof(buffers[i][0]));
	FreeAtExitNoCheck(buffers[i]);
    }
    p0 = p = buffers[i];

    if (no0 < 0)
	no0 = 0;
    if (no1 < 0 || no1 > type->nbas)
	no1 = type->nbas;
    if (no1 <= no0)
	return NULL;

#if 0
if (!memcmp(type->name, "TST", 5)) {
    static int *map = NULL;
    DOF_TYPE *type1 = DOF_Pn[type->name[3] - '0'];
    const FLOAT *b;
    if (map == NULL) {
	const FLOAT *pts0 = type->points, *pts1 = type1->points;
	int k = 0, j;
	map = phgAlloc(type->nbas * sizeof(*map));
	for (i = 0; i < NVert; i++, k++)
	    map[k] = k;
	phgPrintf("0 - %g <-> %g\n", (double)pts0[0], (double)pts1[0]);
	pts0 += 1;
	pts1 += 1;
	for (i = 0; i < type->np_edge; i++) {
	    phgPrintf("%d - %g %g <-> %g %g\n", i,
			(double)pts0[2 * i + 0], (double)pts0[2 * i + 1],
			(double)pts1[2 * i + 0], (double)pts1[2 * i + 1]);
	    for (j = 0; j < type->np_edge; j++)
		if (Fabs(pts0[2 * i + 0] - pts1[2 * j + 0]) <= FLOAT_EPSILON &&
		    Fabs(pts0[2 * i + 1] - pts1[2 * j + 1]) <= FLOAT_EPSILON)
		   break;
		assert(j < type->np_edge);
	    map[k + i] = k + j;
	}
	for (i = 1; i < NEdge; i++)
	    for (j = 0; j < type->np_edge; j++)
		map[k + i * type->np_edge + j] = map[k + j] + i * type->np_edge;
	k += NEdge * type->np_edge;
	pts0 += 2 * type->np_edge;
	pts1 += 2 * type->np_edge;
	for (i = 0; i < type->np_face; i++) {
	    phgPrintf("%d - %g %g %g <-> %g %g %g\n", i,
			(double)pts0[3 * i + 0], (double)pts0[3 * i + 1],
			(double)pts0[3 * i + 2],
			(double)pts1[3 * i + 0], (double)pts1[3 * i + 1],
			(double)pts1[3 * i + 2]);
	    for (j = 0; j < type->np_face; j++)
		if (Fabs(pts0[3 * i + 0] - pts1[3 * j + 0]) <= FLOAT_EPSILON &&
		    Fabs(pts0[3 * i + 1] - pts1[3 * j + 1]) <= FLOAT_EPSILON &&
		    Fabs(pts0[3 * i + 2] - pts1[3 * j + 2]) <= FLOAT_EPSILON)
		   break;
		assert(j < type->np_face);
	    map[k + i] = k + j;
	}
	for (i = 1; i < NFace; i++)
	    for (j = 0; j < type->np_face; j++)
		map[k + i * type->np_face + j] = map[k + j] + i * type->np_face;
	k += NFace * type->np_face;
	pts0 += 3 * type->np_face;
	pts1 += 3 * type->np_face;
	for (i = 0; i < type->np_elem; i++) {
	    phgPrintf("%d - %g %g %g %g <-> %g %g %g %g\n", i,
			(double)pts0[4 * i + 0], (double)pts0[4 * i + 1],
			(double)pts0[4 * i + 2], (double)pts0[4 * i + 3],
			(double)pts1[4 * i + 0], (double)pts1[4 * i + 1],
			(double)pts1[4 * i + 2], (double)pts1[4 * i + 3]);
	    for (j = 0; j < type->np_elem; j++)
		if (Fabs(pts0[4 * i + 0] - pts1[4 * j + 0]) <= FLOAT_EPSILON &&
		    Fabs(pts0[4 * i + 1] - pts1[4 * j + 1]) <= FLOAT_EPSILON &&
		    Fabs(pts0[4 * i + 2] - pts1[4 * j + 2]) <= FLOAT_EPSILON &&
		    Fabs(pts0[4 * i + 3] - pts1[4 * j + 3]) <= FLOAT_EPSILON)
		   break;
		assert(j < type->np_elem);
	    map[k + i] = k + j;
	}
	for (i = 0; i < type->nbas; i++)
	    phgPrintf("map[%d] = %d\n", i, map[i]);
    }
    dof->type = type1;
    b = type1->BasFuncs(dof, e, 0, -1, lambda);
    dof->type = type;
    for (i = no0; i < no1; i++)
	*(p++) = b[map[i]];
    return p0;
}
#endif	/* 0 */

    if (type->mass_lumping->coeff == NULL) {
#if USE_OMP
# pragma omp critical (ML_bas)
#endif	/* USE_OMP */
	init_buffers(type);
    }

    bas_lambda = type->mass_lumping->bas_lambda + (Dim + 1) * phgThreadId;
    bas_values = type->mass_lumping->bas_values + type->nbas * phgThreadId;
    if (memcmp(bas_lambda, lambda, (Dim + 1) * sizeof(*lambda))) {
	memcpy(bas_lambda, lambda, (Dim + 1) * sizeof(*lambda));
	get_monomials(type, lambda, bas_values, FALSE);
    }

    for (; no0 < no1; no0++) {
	v = 0.;
	for (i = 0; i < type->nbas; i++)
	    v += type->mass_lumping->coeff[no0  + i * type->nbas]
			* bas_values[i];
	*(p++) = v;
    }

    return p0;
}

static const FLOAT *
ML_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT **buffers = NULL;
#if USE_OMP
#pragma omp threadprivate(buffers)
#endif  /* USE_OMP */
    int i;
    FLOAT *p, *p0, *grad_lambda, *grad_values, v0, v1, v2, v3;
    DOF_TYPE *type = dof->type;
    
    assert(type != NULL && type->mass_lumping != NULL);

    if (buffers == NULL) {
	buffers = phgCalloc(ORDER_MAX, sizeof(*buffers));
	FreeAtExitNoCheck(buffers);
    }

    i = type->mass_lumping->order0;
    assert(i <= ORDER_MAX);
    if (buffers[i] == NULL) {
	buffers[i] = phgCalloc(type->nbas * (Dim + 1),
					sizeof(buffers[i][0]));
	FreeAtExitNoCheck(buffers[i]);
    }
    p0 = p = buffers[i];

    if (no0 < 0)
	no0 = 0;
    if (no1 < 0 || no1 > type->nbas)
	no1 = type->nbas;
    if (no1 <= no0)
	return NULL;

    if (type->mass_lumping->coeff == NULL) {
#if USE_OMP
# pragma omp critical (ML_grad)
#endif	/* USE_OMP */
	init_buffers(type);
    }

    grad_lambda = type->mass_lumping->grad_lambda + (Dim + 1) * phgThreadId;
    grad_values = type->mass_lumping->grad_values
					+ (Dim + 1) * type->nbas * phgThreadId;
    if (memcmp(grad_lambda, lambda, (Dim + 1) * sizeof(*lambda))) {
	memcpy(grad_lambda, lambda, (Dim + 1) * sizeof(*lambda));
	get_monomials(type, lambda, grad_values, TRUE);
    }

    for (; no0 < no1; no0++) {
	v0 = v1 = v2 = v3 = 0.;
	for (i = 0; i < type->nbas; i++) {
	    v0 += type->mass_lumping->coeff[no0 + i * type->nbas] *
		  grad_values[4 * i];
	    v1 += type->mass_lumping->coeff[no0 + i * type->nbas] *
		  grad_values[4 * i + 1];
	    v2 += type->mass_lumping->coeff[no0 + i * type->nbas] *
		  grad_values[4 * i + 2];
	    v3 += type->mass_lumping->coeff[no0 + i * type->nbas] *
		  grad_values[4 * i + 3];
	}
	*(p++) = v0;
	*(p++) = v1;
	*(p++) = v2;
	*(p++) = v3;
    }

    return p0;
}

/* Note: d = 2 (face) or 3 (element), n = orbit length */
#define MAP_LOC(u0,u1,u2,u3, d, n)	\
	((d == 2 ? u0 * (d + 1) + u1 : (u0 * (d + 1) + u1) * (d + 1) + u2) * n)
#define SET_MAP(u0,u1,u2,u3, d, n, map)				\
    for (i = 0; i < n; i++) {					\
	temp1[0] = temp0[(d + 1) * i + u0];			\
	temp1[1] = temp0[(d + 1) * i + u1];			\
	temp1[2] = temp0[(d + 1) * i + u2];			\
	for (j = 0; j < n; j++)					\
	    if (!fltcmp(temp1, temp0 + (d + 1) * j, d))	\
		break;						\
	assert(j < n);						\
	map[MAP_LOC(u0,u1,u2,u3,d,n) + j] = i;			\
    }

#define SET_MAP3(n, map) 					\
    SET_MAP(0,1,2, -1, 2,n,map) SET_MAP(0,2,1, -1, 2,n,map)	\
    SET_MAP(1,0,2, -1, 2,n,map) SET_MAP(1,2,0, -1, 2,n,map)	\
    SET_MAP(2,0,1, -1, 2,n,map) SET_MAP(2,1,0, -1, 2,n,map)

static void
ML_fmap(const DOF_TYPE *type, const ELEMENT *e, int face_no, int *M)
/* maps ordering of local face data to ordering of data in DOF */
{
    static BYTE *map21 = NULL, *map111 = NULL;
    int i, j, k, n, size, v0, v1, v2, u0, u1, u2;

    assert(type != NULL && type->mass_lumping != NULL);

    if (type->np_face <= 1) {
	for (k = 0; k < type->np_face; k++)
	    (*M++) = k;
	return;
    }

    v0 = MapVertex(e, GetFaceVertex(face_no, 0));
    v1 = MapVertex(e, GetFaceVertex(face_no, 1));
    v2 = MapVertex(e, GetFaceVertex(face_no, 2));

    /* map (v0,v1,v2) to (u0,u1,u2), the latter is a permutation of {0,1,2} */
    u0 = u1 = u2 = 0;
    (v0 > v1) ? u0++ : u1++;
    (v0 > v2) ? u0++ : u2++;
    (v1 > v2) ? u1++ : u2++;

    for (n = k = 0; k < type->mass_lumping->norbit; k++) {
	switch (type->mass_lumping->orbits[k]) {
	    case S1:
	    case S2:
	    case S11:
	    case S4:
	    case S31:
	    case S22:
	    case S211:
	    case S1111:
		break;
	    case S3:
		*(M++) = n++;
		break;
	    case S21:
		size = 3;
#if USE_OMP
#pragma omp critical (ML_fmap21)
#endif  /* USE_OMP */
		if (map21 == NULL) {
		    FLOAT temp0[] = {Perm21(.25)};
		    FLOAT temp1[3];
		    map21 = phgCalloc(3 * 3 * size, sizeof(*map21));
		    FreeAtExitNoCheck(map21);
		    SET_MAP3(size, map21)
		}
		j = MAP_LOC(u0, u1, u2, -1, 2, size);
		for (i = 0; i < size; i++)
		    *(M++) = n + map21[j + i];
		n += size;
		break;
	    case S111:
		size = 6;
#if USE_OMP
#pragma omp critical (ML_fmap111)
#endif  /* USE_OMP */
		if (map111 == NULL) {
		    FLOAT temp0[] = {Perm111(.125,.375)};
		    FLOAT temp1[3];
		    map111 = phgCalloc(3 * 3 * size, sizeof(*map111));
		    FreeAtExitNoCheck(map111);
		    SET_MAP3(size, map111)
		}
		j = MAP_LOC(u0, u1, u2, -1, 2, size);
		for (i = 0; i < size; i++)
		    *(M++) = n + map111[j + i];
		n += size;
		break;
	    default:
		phgError(1, "invalid orbit type (%d), abort.\n",
				type->mass_lumping->orbits[k]);
		break;
	}
    }

    assert(n == type->np_face);
}

#define SET_MAP4(n, map)					\
    SET_MAP(0,1,2,3, 3,n,map) SET_MAP(0,2,1,3, 3,n,map)	\
    SET_MAP(1,0,2,3, 3,n,map) SET_MAP(1,2,0,3, 3,n,map)	\
    SET_MAP(2,0,1,3, 3,n,map) SET_MAP(2,1,0,3, 3,n,map)	\
							\
    SET_MAP(0,1,3,2, 3,n,map) SET_MAP(0,2,3,1, 3,n,map)	\
    SET_MAP(1,0,3,2, 3,n,map) SET_MAP(1,2,3,0, 3,n,map)	\
    SET_MAP(2,0,3,1, 3,n,map) SET_MAP(2,1,3,0, 3,n,map)	\
							\
    SET_MAP(0,3,1,2, 3,n,map) SET_MAP(0,3,2,1, 3,n,map)	\
    SET_MAP(1,3,0,2, 3,n,map) SET_MAP(1,3,2,0, 3,n,map)	\
    SET_MAP(2,3,0,1, 3,n,map) SET_MAP(2,3,1,0, 3,n,map)	\
							\
    SET_MAP(3,0,1,2, 3,n,map) SET_MAP(3,0,2,1, 3,n,map)	\
    SET_MAP(3,1,0,2, 3,n,map) SET_MAP(3,1,2,0, 3,n,map)	\
    SET_MAP(3,2,0,1, 3,n,map) SET_MAP(3,2,1,0, 3,n,map)

static void
ML_elem_map(const DOF_TYPE *type, const ELEMENT *e, int *M)
/* maps ordering of local element data to ordering of data in DOF */
{
    static BYTE *map31 = NULL, *map22 = NULL, *map211 = NULL, *map1111 = NULL;
    int i, j, k, n, size, u0, u1, u2;

    assert(type != NULL && type->mass_lumping != NULL);

    if (type->np_elem <= 1 || ElementIsInOrder(e)) {
	for (k = 0; k < type->np_elem; k++)
	    (*M++) = k;
	return;
    }

    u0 = MapVertex(e, 0);
    u1 = MapVertex(e, 1);
    u2 = MapVertex(e, 2);

    for (n = k = 0; k < type->mass_lumping->norbit; k++) {
	switch (type->mass_lumping->orbits[k]) {
	    case S1:
	    case S2:
	    case S11:
	    case S3:
	    case S21:
	    case S111:
		break;
	    case S4:
		*(M++) = n++;
		break;
	    case S31:
		size = 4;
#if USE_OMP
#pragma omp critical (ML_emap31)
#endif  /* USE_OMP */
		if (map31 == NULL) {
		    FLOAT temp0[] = {Perm31(.125)};
		    FLOAT temp1[3];
		    map31 = phgCalloc(4 * 4 * 4 * size, sizeof(*map31));
		    FreeAtExitNoCheck(map31);
		    SET_MAP4(size, map31)
		}
		j = MAP_LOC(u0, u1, u2, 1 - u0 - u1 - u2, 3, size);
		for (i = 0; i < size; i++)
		    *(M++) = n + map31[j + i];
		n += size;
		break;
	    case S22:
		size = 6;
#if USE_OMP
#pragma omp critical (ML_emap22)
#endif  /* USE_OMP */
		if (map22 == NULL) {
		    FLOAT temp0[] = {Perm22(.125)};
		    FLOAT temp1[3];
		    map22 = phgCalloc(4 * 4 * 4 * size, sizeof(*map22));
		    FreeAtExitNoCheck(map22);
		    SET_MAP4(size, map22)
		}
		j = MAP_LOC(u0, u1, u2, 1 - u0 - u1 - u2, 3, size);
		for (i = 0; i < size; i++)
		    *(M++) = n + map22[j + i];
		n += size;
		break;
	    case S211:
		size = 12;
#if USE_OMP
#pragma omp critical (ML_emap211)
#endif  /* USE_OMP */
		if (map211 == NULL) {
		    FLOAT temp0[] = {Perm211(.125,.25)};
		    FLOAT temp1[3];
		    map211 = phgCalloc(4 * 4 * 4 * size, sizeof(*map211));
		    FreeAtExitNoCheck(map211);
		    SET_MAP4(size, map211)
		}
		j = MAP_LOC(u0, u1, u2, 1 - u0 - u1 - u2, 3, size);
		for (i = 0; i < size; i++)
		    *(M++) = n + map211[j + i];
		n += size;
		break;
	    case S1111:
		size = 24;
#if USE_OMP
#pragma omp critical (ML_emap1111)
#endif  /* USE_OMP */
		if (map1111 == NULL) {
		    FLOAT temp0[] = {Perm1111(.0625,.125,.1875)};
		    FLOAT temp1[3];
		    map1111 = phgCalloc(4 * 4 * 4 * size, sizeof(*map1111));
		    FreeAtExitNoCheck(map1111);
		    SET_MAP4(size, map1111)
		}
		j = MAP_LOC(u0, u1, u2, 1 - u0 - u1 - u2, 3, size);
		for (i = 0; i < size; i++)
		    *(M++) = n + map1111[j + i];
		n += size;
		break;
	    default:
		phgError(1, "invalid orbit type (%d), abort.\n",
				type->mass_lumping->orbits[k]);
		break;
	}
    }

    assert(n == type->np_elem);
}

static void
ML_map(const DOF_TYPE *type, const ELEMENT *e, int no, int *M)
{
    if (no < 0)
	ML_elem_map(type, e, M);
    else
	ML_fmap(type, e, no, M);
}

#undef MAP_LOC
#undef SET_MAP
#undef SET_MAP3
#undef SET_MAP4

/*---------------------------------------------------------------------------
 * Rules for transforming the constrained orbits.
 *
 *   Vertex nodes:
 * 	Cons31(0)	-> Cons1(1)		Vertex
 *   Edge nodes:
 * 	Cons22(x)	-> Cons2(0.5)		Edge center
 * 	Cons211(0,b) 	-> Cons11(b)		Edge (s11)
 *   Face nodes:
 *	Cons31(1/3)	-> Cons3(1./3.)		Face center
 *	Cons211(a,0)	-> Cons21(a)		Face (s21)
 *	Cons1111(0,b,c)	-> Cons111(b,c)		Face (s111)
 *
 * The script 'utils/quad2ml' sorts and converts orbits using the rules above.
 *---------------------------------------------------------------------------*/

/* The macro 'DEFINE_ML_TYPE' defines a DOF_TYPE using macros WEIGHTS, POINTS,
 * and ORBITS */

/* The following are C-preprocessor dependent */
#define Orbits(type)	type ## _orbits
#define Weights0(type)	type ## _weights0
#define MLstruct(type)	type ## _MLstruct
#define Points(type)	type ## _points
#define MonoBas(type)	type ## _monobas
#define Dof(type)	DOF_ ## type ## _
#define Name(type)	#type

#define	DEFINE_ML_TYPE0(O0, O1, O2, O3, type, grad_type, space, cont)	\
    static int		Orbits(type)[] = {ORBITS};			\
    static FLOAT	Weights0(type)[] = {WEIGHTS};			\
    static FLOAT	Points(type)[] = {POINTS};			\
    static char		MonoBas(type)[] = {BASES};			\
    static MASS_LUMPING	MLstruct(type) = {				\
	O0, O1, O2, O3, 						\
	sizeof(Orbits(type)) / sizeof(Orbits(type)[0]),			\
	sizeof(Points(type)) / sizeof(Points(type)[0]),			\
	sizeof(MonoBas(type)) / sizeof(MonoBas(type)[0]),		\
	sizeof(Weights0(type)) / sizeof(Weights0(type)[0]),		\
	Orbits(type),		/* orbits[] */				\
	Weights0(type),		/* weights0[] */			\
	MonoBas(type),		/* monobas[][4] */			\
	/* The pointers below are internal buffers */			\
	NULL,			/* weights[] */				\
	NULL, NULL,		/* face_pow, elem_pow */		\
	NULL,			/* coeff */				\
	NULL, NULL,		/* bas_lambda, bas_values */		\
	NULL, NULL		/* grad_lambda, grad_values */		\
    };									\
    DOF_TYPE Dof(type) = {						\
	DofCache, &MLstruct(type), init, Name(type), Points(type),	\
	NULL,						/* orders[] */	\
	grad_type,					/* grad_type */	\
	NULL,						/* base_type */	\
	NULL,						/* hp_info */	\
	phgDofInterpC2FGeneric,				/* Prolong. */	\
	phgDofInterpF2CGeneric,				/* Interp. */	\
	phgDofInitFuncPoint,				/* Ini. func */	\
	ML_bas, ML_grad, ML_map,	 				\
	space,						/* FE_space */	\
	TRUE,						/* is_nodal */ \
	TRUE,						/* invariant */	\
	FALSE,						/* free_after_use */ \
	-1,						/* id */	\
	-1,						/* nbas */ \
	-1,						/* order */	\
	0,						/* D_order */	\
	cont,						/* continuity */\
	1,						/* dimension */	\
	-1, -1, -1, -1};

#define	DEFINE_ML_TYPE(O0, O1, O2, O3, type, grad_type) \
	DEFINE_ML_TYPE0(O0, O1, O2, O3, type, grad_type, FE_H1, 0)

/* Permutation macros for constructing symmetric set of monomials */
#define Mono4(a)	a,a,a,a
#define Mono31(a,b)	a,a,a,b, a,a,b,a, a,b,a,a, b,a,a,a
#define Mono22(a,b)	a,a,b,b, a,b,a,b, a,b,b,a, b,b,a,a, b,a,b,a, b,a,a,b
#define Mono211(a,b,c)	a,a,b,c, a,a,c,b, a,b,a,c, a,b,c,a, a,c,a,b, a,c,b,a, \
			b,a,a,c, b,a,c,a, b,c,a,a, c,a,a,b, c,a,b,a, c,b,a,a
#define Mono0111(p,a,b,c) p,a,b,c, p,a,c,b, p,b,a,c, p,b,c,a, p,c,a,b, p,c,b,a
#define Mono1111(a,b,c,d) \
    Mono0111(a,b,c,d), Mono0111(b,a,c,d), Mono0111(c,a,b,d), Mono0111(d,a,b,c)

/*-------------------------------- Begin ML1 --------------------------------*/
/* The following macros are generated from the corresponding constrained
 * quadrature rule by the script 'utils/quad2ml' */
#define POINTS	Cons1(1.)
#define WEIGHTS	Dup1(.25)
#define ORBITS	S1
#define BASES	-1
DEFINE_ML_TYPE(1, 1, 1, 1, ML1, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End ML1 ---------------------------------*/

/*-------------------------------- Begin ML2 --------------------------------*/
/* Order 2 14-node basis found with "./quadrule-ml 3 2 2 2 5".
 * The actual solution was the one with largest minimal weight, obtained
 * using symbolic computations ("quadrule dump=...") */
#define POINTS \
	Cons1(1.), \
	Cons2(.5), \
	Perm31(.31685117118244328709490761618832365)
#define WEIGHTS \
	Dup1(.02431334243766409389014280607218878), \
	Dup2(.02431334243766409389014280607218878), \
	Dup31(.18921664390583976527464298481952805)
#define ORBITS	S1,S2,S31
#define BASES \
	Mono31(0,2), \
	Mono22(0,1), \
	Mono31(1,2)
DEFINE_ML_TYPE(2, -1, -1, -1, ML2, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End ML2 ---------------------------------*/

/*-------------------------------- Begin ML2m --------------------------------*/
/* Order 2 15-node basis from the paper by S. Geevers et.al,
 *	New higher-order mass-lumped tetrahedral elements for wave propagation
 *	modelling, SIAM J. Sci. Comput., 2018. */
#define POINTS \
	Cons1(1.), \
	Cons2(.5), \
	Cons3(.3333333333333333333333333333333333), \
	Perm4(.25)
#define WEIGHTS \
	/* 6*17/5040 */ \
	Dup1(.02023809523809523809523809523809524), \
	/* 6*2/315 */ \
	Dup2(.03809523809523809523809523809523810), \
	/* 6*9/560 */ \
	Dup3(.09642857142857142857142857142857143), \
	/* 6*16/315 */ \
	Dup4(.30476190476190476190476190476190476)
#define ORBITS	S1,S2,S3,S4
#define BASES \
	Mono31(0,1), \
	Mono22(1,0), \
	Mono31(1,0), \
	Mono4(1)
DEFINE_ML_TYPE(2, -1, -1, -1, ML2m, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End ML2m ---------------------------------*/

/*-------------------------------- Begin ML3 --------------------------------*/
/* Order 3 32-node basis from the paper by S. Geevers et.al,
 *	New higher-order mass-lumped tetrahedral elements for wave propagation
 *	modelling, SIAM J. Sci. Comput., 2018. */
#define POINTS \
	Cons1 (1.), \
	/* (3-sqrt(3*(sqrt(2)-1)))/6 */ \
	Cons11(.314210342418032893883175061743391), \
	/* (4-sqrt(2))/12 */ \
	Cons21(.215482203135575412599859272982525), \
	/* 1/6 */ \
	Perm31(.166666666666666666666666666666667)
#define WEIGHTS \
	/* 6*(41-9*sqrt(2))/41160 */ \
	Dup1 (.004121294160151916116732478350162204), \
	/* 6*(8+9*sqrt(2))/13720 */ \
	Dup11(.009064688948115680221231136378084818), \
	/* 6*(10-sqrt(2))/1715 */ \
	Dup21(.030037736808024157263667561314718257), \
	/* 6*3/140 */ \
	Dup31(.128571428571428571428571428571428571)
#define ORBITS	S1,S11,S21,S31
#define BASES \
	Mono31(0,1), \
	Mono211(0,2,1), \
	Mono211(1,2,0), \
	Mono31(1,2)
DEFINE_ML_TYPE(3, -1, -1, -1, ML3, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End ML3 ---------------------------------*/

/*-------------------------------- Begin ML4 --------------------------------*/
/* Order 4 65-node basis from the paper by S. Geevers et.al,
 *	New higher-order mass-lumped tetrahedral elements for wave propagation
 *	modelling, SIAM J. Sci. Comput., 2018. */
#define POINTS \
	Cons1(1.), \
	Cons11(.17249194077490857535314475496313593), \
	Cons2 (.5), \
	Cons21(.14741779690136858426180821049471841), \
	Cons21(.45403952722710673616480685584005258), \
	Cons3 (.33333333333333333333333333333333333), \
	Perm31(.12822093162909793115292586197383820), \
	Perm22(.08742182088664352997036484977702327), \
	Perm31(.31240614520708106211729705978204781), \
	Perm4 (.25)
#define WEIGHTS \
	Dup1 (.00072962552706739279394822308370445), \
	Dup11(.00282247451924664655730192913050967), \
	Dup2 (.00106023955505008515322220389159024), \
	Dup21(.01184849151957705922982759713985786), \
	Dup21(.00715479187061820479339285166833011), \
	Dup3 (.00626818558580473610057070895289025), \
	Dup31(.05304855114341457450055043729745327), \
	Dup22(.04134607754640933979637644238744208), \
	Dup31(.04499738112310261786636956547553540), \
	Dup4 (.06347802896038323829038384782310102)
#define ORBITS	S1,S11,S2,S21,S21,S3,S31,S22,S31,S4
#define BASES \
	Mono31(0,1), \
	Mono211(0,2,1), \
	Mono22(2,0), \
	Mono211(1,2,0), \
	Mono211(2,1,0), \
	Mono31(2,0), \
	Mono31(1,2), \
	Mono22(2,1), \
	Mono31(2,1), \
	Mono4(2)
DEFINE_ML_TYPE(4, -1, -1, -1, ML4, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End ML4 ---------------------------------*/

/*-------------------------------- Begin ML4a --------------------------------*/
/* Order 4 60-node basis from the paper by S. Geevers et.al,
 *	New higher-order mass-lumped tetrahedral elements for wave propagation
 *	modelling, SIAM J. Sci. Comput., 2018. */
#define POINTS \
	Cons1 (1.), \
	Cons11(.16148658334966762839187222355198378), \
	Cons2 (.5), \
	Cons21(.14902192884695984835550556230941622), \
	Cons21(.39445919721717828782908721585020459), \
	Perm31(.13020588463725643564034875550580697), \
	Perm22(.06386116838612691411152895029275609), \
	Perm31(.30121792340790865275261948454963855)
#define WEIGHTS \
	Dup1 (.00055914881734603058563689056792032), \
	Dup11(.00289759942588405879073894765409403), \
	Dup2 (.00120330227528155206029421709116674), \
	Dup21(.01201862451504915271393453928966253), \
	Dup21(.00676109620080009857208399429186061), \
	Dup31(.05495546693997779061882907121885302), \
	Dup22(.04035193592435867996125286692982867), \
	Dup31(.06712056651801590053294096847488199)
#define ORBITS	S1,S11,S2,S21,S21,S31,S22,S31
#define BASES \
	Mono31(0,1), \
	Mono211(0,2,1), \
	Mono22(2,0), \
	Mono211(1,2,0), \
	Mono211(2,1,0), \
	Mono31(1,2), \
	Mono22(2,1), \
	Mono31(2,1)
DEFINE_ML_TYPE(4, -1, -1, -1, ML4a, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End ML4a ---------------------------------*/

/*-------------------------------- Begin ML4b --------------------------------*/
/* Order 4 61-node basis from the paper by S. Geevers et.al,
 *	New higher-order mass-lumped tetrahedral elements for wave propagation
 *	modelling, SIAM J. Sci. Comput., 2018. */
#define POINTS \
	Cons1 (1.), \
	Cons11(.20016281047078478593140176926489303), \
	Cons2 (.5), \
	Cons21(.13973509722383662156257565315513392), \
	Cons21(.43194362351776820994051366110052126), \
	Perm31(.12822093162909793115292586197383820), \
	Perm22(.08742182088664352997036484977702327), \
	Perm31(.31240614520708106211729705978204781), \
	Perm4 (.25)
#define WEIGHTS \
	Dup1 (.00095584162254363836517722082676238), \
	Dup11(.00267679510900574363020163923685022), \
	Dup2 (.00222949796742357579240133338574757), \
	Dup21(.01130576978794261389902869248966678), \
	Dup21(.00927255363641630457481638520171274), \
	Dup31(.05304855114341457450055043729745327), \
	Dup22(.04134607754640933979637644238744208), \
	Dup31(.04499738112310261786636956547553540), \
	Dup4 (.06347802896038323829038384782310102)
#define ORBITS	S1,S11,S2,S21,S21,S31,S22,S31,S4
#define BASES \
	Mono31(0,1), \
	Mono211(0,2,1), \
	Mono22(2,0), \
	Mono211(1,2,0), \
	Mono211(2,1,0), \
	Mono31(1,2), \
	Mono22(2,1), \
	Mono31(2,1), \
	Mono4(2)
DEFINE_ML_TYPE(4, -1, -1, -1, ML4b, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End ML4b ---------------------------------*/

/*-------------------------------- Begin ML2o --------------------------------*/
/* The following macros are generated from the corresponding constrained
 * quadrature rule by the script 'utils/quad2ml' */
#define POINTS \
    Cons1(1.),\
    Cons2(0.5),\
    Cons21(.18858048469644503927115437402941689),\
    Perm4(.25)
#define WEIGHTS \
    Dup1(.00129961081762382864323948583189792),\
    Dup2(.00751330903878115632153864252437151),\
    Dup21(.05374664981124274879942069139702288),\
    Dup4(.30476190476190476190476190476190476)
#define ORBITS	S1,S2,S21,S4
#define BASES	-1
DEFINE_ML_TYPE(2, 2, 4, 4, ML2o, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End ML2o ---------------------------------*/

/*-------------------------------- Begin ML3o --------------------------------*/
/* The following macros are generated from the corresponding constrained
 * quadrature rule by the script 'utils/quad2ml' */
#define POINTS \
    Cons1(1.),\
    Cons11(.69474012433043399696452214373799579),\
    Cons21(.14804629800083268127155972966379546),\
    Cons21(.42045997555404366358045983652415104),\
    Perm31(.10486452489170353984304348629272939),\
    Perm22(.12587961966825070744487424837645650)
#define WEIGHTS \
    Dup1(.00139318132340935829885454503458246),\
    Dup11(.00439720814497923320460068363824448),\
    Dup21(.00938677154027050153276493215837929),\
    Dup21(.01517875558886845408624234213923300),\
    Dup31(.04276746867938747299136451172239785),\
    Dup22(.07930076278323240149263804629029959)
#define ORBITS	S1,S11,S21,S21,S31,S22
#define BASES	-1
DEFINE_ML_TYPE(3, 3, 5, 6, ML3o, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End ML3o ---------------------------------*/

/*-------------------------------- Begin ML4o --------------------------------*/
/* WARNING: the current rule has negative weights and should not be used */
/* The following macros are generated from the corresponding constrained
 *  * quadrature rule by the script 'utils/quad2ml' */
#define POINTS \
    Cons1(1.),\
    Cons11(.80854443961101619213427935540700553),\
    Cons11(.30601867253492915323828632954666129),\
    Cons3(1./3.),\
    Cons21(.25867647251242276280547601945611744),\
    Cons111(.18358911895300261675702374239715381,\
	    .00799833687507618645573675055882157),\
    Perm31(.30062604540622798826003985009648294),\
    Perm31(.12819524754074949404694586327517572),\
    Perm211(.07304272679126745377338651048218843,\
	    .35770884123478518379857662808127490)
#define WEIGHTS \
    Dup1(.00036028843803848727831482031703821),\
/**/Dup11(-.04188920908297572294780990987494447),\
    Dup11(.00437666250190729126960626330898023),\
/**/Dup3(-.02677414236318327749238687899671355),\
    Dup21(.02167000187725150565693278672597593),\
    Dup111(.02421578012002043115377259617036411),\
    Dup31(.06494829681676811150761297509556652),\
    Dup31(.04806667921869404922983168397071333),\
    Dup211(.02187727709367027353926813403705858)
#define ORBITS	S1,S11,S11,S3,S21,S111,S31,S31,S211
#define BASES	-1
DEFINE_ML_TYPE(4, 5, 6, 7, ML4o, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End ML4o ---------------------------------*/

/*-------------------------------- Begin P1B --------------------------------*/
/* This implements the $P_1 \oplus \span{bubble} FE space which may be used
 * to implement, e.g., the multiscale stabilization algorithms:
 * 	Lina Song, Yanren Hou, Haibiao Zheng,
 * 	    A variational multiscale method based on bubble functions
 * 		for convection-dominated convection¨Cdiffusion equation,
 * 	Applied Mathematics and Computation, 217 (2010) 2226¨C2237
 *
 * Note: examples/poisson.c does not work when using DG3 as the gradient type.
 */
#define POINTS	Cons1(1.), Perm4(.25)
#define WEIGHTS	-1
#define ORBITS	S1,S4
#define BASES	-1
DEFINE_ML_TYPE(1, 1, 1, 4, P1B, DOF_DG0)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
/*--------------------------------- End P1B ---------------------------------*/

DOF_TYPE *DOF_MLn[] = { DOF_ML1, DOF_ML2, DOF_ML3, DOF_ML4, DOF_ML4a, DOF_ML4b,
			DOF_P1B, NULL};

DOF_TYPE *DOF_MLon[] = {DOF_ML2m, DOF_ML2o, DOF_ML3o, DOF_ML4o, NULL};

/*-------------------------------- Begin MLDG --------------------------------*/
/* MLDG1 (quadrature order: 2) */
#define POINTS	Perm31(.13819660112501051517954131656343619)
#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S31
#define BASES	-1
DEFINE_ML_TYPE0(0, 0, 0, 1, MLDG1, DOF_DG0, FE_L2, -1)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS

/* MLDG2 (quadrature order: 3, multiple solutions) */
#define POINTS	Perm31(.11241946973578680000000000000000000), \
		Perm22(.09301814747511176528289006204768939)
#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S31,S22
#define BASES	-1
DEFINE_ML_TYPE0(0, 0, 0, 2, MLDG2, DOF_DG1, FE_L2, -1)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS

/* MLDG3 (quadrature order: 5, multiple solutions) */
#if 1
#define POINTS	Perm31(.08980373984413880875518672073678476), \
		Perm31(.28545818671466010000000000000000000), \
		Perm211(.42640467172723715351818908060930334,.02013270817172337294732700260455974)
#else
/* FIXME: this results in condition number 8.68e+11! */
#define POINTS	Perm31(.08945436400890264579966883396740960), \
		Perm31(.25000718218609900000000000000000000), \
		Perm211(.42143943120900576379663436549581748,.13258109984515821268875827682570903)
#endif

#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S31,S31,S211
#define BASES	-1
DEFINE_ML_TYPE0(0, 0, 0, 3, MLDG3, DOF_DG2, FE_L2, -1)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS

DOF_TYPE *DOF_MLDGn[] = {&DOF_MLDG1_, &DOF_MLDG2_, &DOF_MLDG3_, NULL};

/*--------------------------------- End MLDG ---------------------------------*/

#if TEST_TYPES
/* Some ordinary, non mass-lumping bases */

#if 0
/*-------------------------------------------- Test basis of orders 2 3 4 5 */
#define POINTS	Cons1  (2./2.),		/* vert: 2 */ \
		Cons11 (1./3.),		/* edge: 12 */ \
		Cons21 (1./4.),		/* face: 112 */ \
		Perm31 (1./5.),		/* elem: 1112 */
#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S1,S11,S21,S31
#define BASES	-1
DEFINE_ML_TYPE(2, 3, 4, 5, TEST, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
#elif 0
/*--------------------- The 3rd order basis by Chen & Babuska (CMAME, 1996) */
#define POINTS	Cons1  (1.), \
		Cons11 (0.7251957), \
		Cons3  (0.3333333)
#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S1,S11,S3
#define BASES	-1
DEFINE_ML_TYPE(3, 3, 3, 3, TEST, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
#elif 0
/*--------------------- The 4th order basis by Chen & Babuska (CMAME, 1996) */
#define POINTS	Cons1  (1.), \
		Cons2  (0.5), \
		Cons11 (0.8306024), \
		Cons21 (0.2208880), \
		Perm4  (0.25)
#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S1,S2,S11,S21,S4
#define BASES	-1
DEFINE_ML_TYPE(4, 4, 4, 4, TEST, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
#elif 0
/*--------------------- The 5th order basis by Chen & Babuska (CMAME, 1996) */
#define POINTS	Cons1  (1.), \
		Cons11 (0.8866427), \
		Cons11 (0.6431761), \
		Cons21 (0.1525171), \
		Cons21 (0.4168658), \
		Perm31 (0.1823054)
#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S1,S11,S11,S21,S21,S31
#define BASES	-1
DEFINE_ML_TYPE(5, 5, 5, 5, TEST, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
#elif 0
/*--------------------- The 6th order basis by Chen & Babuska (CMAME, 1996) */
#define POINTS	Cons1  (1.), \
		Cons11 (0.9194021), \
		Cons11 (0.7349105), \
		Cons2  (0.5), \
		Cons3  (0.3333333), \
		Cons21 (0.1097139), \
		Cons111(0.3157892,0.5586077), \
		Perm31 (0.1357838), \
		Perm22 (0.3559336)
#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S1,S11,S11,S2,S3,S21,S111,S31,S22
#define BASES	-1
DEFINE_ML_TYPE(6, 6, 6, 6, TEST, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
#elif 0
/*--------------------- The 7th order basis by Chen & Babuska (CMAME, 1996) */
#define POINTS	Cons1  (1.), \
		Cons11 (0.9398927), \
		Cons11 (0.7957614), \
		Cons11 (0.6042138), \
		Cons21 (0.0871370), \
		Cons21 (0.4494208), \
		Cons21 (0.2663399), \
		Cons111(0.2447528,0.6584392), \
		Perm31 (0.1046666), \
		Perm31 (0.2936310), \
		Perm211(0.1141973,0.4885725)
#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S1,S11,S11,S11,S21,S21,S21,S111,S31,S31,S211
#define BASES	-1
DEFINE_ML_TYPE(7, 7, 7, 7, TEST, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
#elif 0
/*--------------------- The 8th order basis by Chen & Babuska (CMAME, 1996) */
#define POINTS	Cons1  (1.), \
		Cons2  (0.5), \
		Cons11 (0.9533797), \
		Cons11 (0.8375919), \
		Cons11 (0.6801403), \
		Cons21 (0.0627331), \
		Cons21 (0.2153606), \
		Cons21 (0.3891297), \
		Cons111(0.3657423,0.5524728), \
		Cons111(0.1942206,0.7294168), \
		Perm4  (0.25), \
		Perm31 (0.0834511), \
		Perm22 (0.4060462), \
		Perm211(0.0927818,0.5844475), \
		Perm211(0.2432058,0.4157364)
#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S1,S2,S11,S11,S11,S21,S21,S21,S111,S111,S4,S31,S22,S211,S211
#define BASES	-1
DEFINE_ML_TYPE(8, 8, 8, 8, TEST, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
#elif 0
/*--------------------- The 9th order basis by Chen & Babuska (CMAME, 1996) */
#define POINTS	Cons1  (1.), \
		Cons11 (0.9626819), \
		Cons11 (0.8672666), \
		Cons11 (0.7361751), \
		Cons11 (0.5815151), \
		Cons3  (0.3333333), \
		Cons21 (0.0493729), \
		Cons21 (0.4658361), \
		Cons21 (0.1769439), \
		Cons111(0.3020146,0.6309227), \
		Cons111(0.1575680,0.7808733), \
		Cons111(0.3261032,0.4887991), \
		Perm31 (0.0682526), \
		Perm31 (0.2123036), \
		Perm211(0.0774288,0.6543249), \
		Perm211(0.0781245,0.5016150), \
		Perm211(0.2049128,0.5083123), \
		Perm211(0.3552639,0.2073089)
#define WEIGHTS	-1.0			/* no weights */
#define ORBITS	S1,S11,S11,S11,S11,S3,S21,S21,S21,S111,S111,S111,S31,S31,S211,S211,S211,S211
#define BASES	-1
DEFINE_ML_TYPE(9, 9, 9, 9, TEST, NULL)
#undef BASES
#undef WEIGHTS
#undef POINTS
#undef ORBITS
#else
/*--------------------------------------------------- disable the TEST type */
# undef DOF_TEST
# define DOF_TEST NULL
#endif

/*--------------------------------------------------------------------------*/
/* The TST* types, code generated by ../utils/mass-lumping-gen */
#include "mass-lumping.inc"
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
DOF_TYPE *DOF_TSTn[] =
  {DOF_TST1, DOF_TST2, DOF_TST3, DOF_TST4, DOF_TST5, DOF_TST6, DOF_TST7,
   DOF_TST8, DOF_TST9, DOF_TST10,DOF_TST11,DOF_TST12,DOF_TST13,DOF_TST14,
   DOF_TST15,DOF_TST16,DOF_TEST, NULL};
/*--------------------------------------------------------------------------*/

/* Tables of Gauss-Lobatto points used above */
/* TODO: create a separate file gauss-lobatto.c for permanent use */
#undef Perm2
#undef Perm11
#undef Cons2
#undef Cons11
#define Gtab	_phg_GLP
#include "../p4est/gauss-lobatto.inc"

#endif	/* TEST_TYPES */
