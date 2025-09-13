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

/* $Id: hfeb.c,v 1.47 2022/04/09 01:28:34 zlb Exp $ */

#include "phg.h"

#include <stdlib.h>
#include <string.h>
#include <strings.h>	/* for bzero */
#include <math.h>

static const FLOAT *
HB1_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
  /* evaluation of basis functions */
{
    static FLOAT values[2 * NEdge][Dim];
    int I, J, k, n;
    GRID *g = dof->g;
    FLOAT (*nabla)[Dim + 1] = (void *)(phgGeomGetJacobian(g, e));
    FLOAT *p;

    CheckThread

    if (no1 <= 0)
	no1 = 2 * NEdge;

    n = (no0 / 2) * 2;

    for (k = n, p = (FLOAT *)values; k < no1; k += 2, p += 2 * 3) {
	GetEdgeVertices(e, k / 2, I, J);
	p[0] = lambda[J] * nabla[I][0] - lambda[I] * nabla[J][0];
	p[1] = lambda[J] * nabla[I][1] - lambda[I] * nabla[J][1];
	p[2] = lambda[J] * nabla[I][2] - lambda[I] * nabla[J][2];
    }

    for (k = n + 1, p = (FLOAT *)(values + 1); k < no1; k += 2, p += 2 * 3) {
	GetEdgeVertices(e, k / 2, I, J);
	p[0] = lambda[I] * nabla[J][0] + lambda[J] * nabla[I][0];
	p[1] = lambda[I] * nabla[J][1] + lambda[J] * nabla[I][1];
	p[2] = lambda[I] * nabla[J][2] + lambda[J] * nabla[I][2];
    }

    return (FLOAT *)values + (no0 - n) * Dim;
}

static const FLOAT *
HB1_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT values[2 * NEdge][Dim][Dim + 1];
    int I, J, k, n;
    GRID *g = dof->g;
    FLOAT (*nabla)[Dim + 1] = (void *)(phgGeomGetJacobian(g, e));
    FLOAT (*p)[Dim][Dim + 1] = values;

    CheckThread

    if (no1 <= 0)
	no1 = 2 * NEdge;

    n = (no0 / 2) * 2;

    bzero(values, sizeof(values));

    for (k = n, p = values; k < no1; k += 2, p += 2) {
	GetEdgeVertices(e, k / 2, I, J);
	/* grad (\lambda[J]\nabla[I] - \lambda[I]\nabla[J]) */
	(*p)[0][I] = -nabla[J][0];
	(*p)[0][J] = nabla[I][0];

	(*p)[1][I] = -nabla[J][1];
	(*p)[1][J] = nabla[I][1];

	(*p)[2][I] = -nabla[J][2];
	(*p)[2][J] = nabla[I][2];
    }

    for (k = n + 1, p = values + 1; k < no1; k += 2, p += 2) {
	GetEdgeVertices(e, k / 2, I, J);
	/* grad (\lambda[J]\nabla[I] + \lambda[I]\nabla[J]) */
	(*p)[0][I] = nabla[J][0];
	(*p)[0][J] = nabla[I][0];

	(*p)[1][I] = nabla[J][1];
	(*p)[1][J] = nabla[I][1];

	(*p)[2][I] = nabla[J][2];
	(*p)[2][J] = nabla[I][2];

    }

    return (FLOAT *)values + (no0 - n) * Dim * (Dim + 1);
}

static void
HB1_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	 DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	 const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
{
#if 0 /* old code */
    GRID *g = dof->g;
    int i, k, u0, u1;
    INT v0, v1;
    FLOAT x0, y0, z0, x1, y1, z1, x2, y2, z2;
    FLOAT x, y, z;
    COORD *c;
    QUAD *quad;
    FLOAT *lambda, *w, lam[Dim + 1] = {0., 0., 0., 0.};
    FLOAT *temp = phgAlloc(dof->dim * sizeof(*temp));
    FLOAT *fvals = phgAlloc(Dim * dof->dim * sizeof(*fvals));

    Unused(pdofvalues);

    assert(type == EDGE);
    assert(funcvalues == NULL);
    assert(userfunc_lambda != NULL || userfunc != NULL);

    for (k = 0; k < dof->dim; k++)
	temp[k] = 0.0;

    GetEdgeVertices(e, index, u0, u1);
    v0 = e->verts[u0];
    v1 = e->verts[u1];
    x0 = (*(c = g->verts + v0))[0];
    y0 = (*c)[1];
    z0 = (*c)[2];

    x1 = (*(c = g->verts + v1))[0];
    y1 = (*c)[1];
    z1 = (*c)[2];

    x2 = x1 - x0;
    y2 = y1 - y0;
    z2 = z1 - z0;
    quad = phgQuadGetQuad1D(3);
    lambda = quad->points;
    w = quad->weights;
    if (userfunc_lambda == NULL) {
	userfunc(x0 + x2 * .5, y0 + y2 * .5, z0 + z2 * .5, fvals);
    }
    else {
	lam[u0] = lam[u1] = 0.5;
	userfunc_lambda(dof, e, index * 2, lam, fvals);
    }
    for (k = 0; k < dof->dim; k++) {
	*(dofvalues++) = -(x2 * fvals[0] + y2 * fvals[1] + z2 * fvals[2]);
	fvals += Dim;
    }
    fvals -= Dim * dof->dim;
    for (i = 0; i < quad->npoints; i++) {
	if (userfunc_lambda == NULL) {
	    x = x0 + x2 * lambda[1];
	    y = y0 + y2 * lambda[1];
	    z = z0 + z2 * lambda[1];
	    userfunc(x, y, z, fvals);
	}
	else {
	    lam[u0] = lambda[0];
	    lam[u1] = lambda[1];
	    userfunc_lambda(dof, e, index * 2 + 1, lam, fvals);
	}
	for (k = 0; k < dof->dim; k++) {
	    temp[k] -= 3.0 * (1 - 2.0 * lambda[0])
		* (x2 * fvals[0] + y2 * fvals[1] + z2 * fvals[2]) * w[i];
	    fvals += Dim;
	}
	fvals -= Dim * dof->dim;
	lambda += 2;
    }
    for (k = 0; k < dof->dim; k++)
	*(dofvalues++) = temp[k];
    phgFree(temp);
#else	/* new code based on L2 projection */
    GRID *g = dof->g;
    FLOAT lam[] = {0., 0., 0., 0.};
    int i, j, k, n, u0, u1, n0, n1, m = dof->dim;
    FLOAT x0, y0, z0, x1, y1, z1, x2, y2, z2, x, y, z;
    COORD *c;
    QUAD *quad;
    FLOAT *fvals = phgAlloc(Dim * m * sizeof(*fvals));
    FLOAT *temp0, *temp1, *A, *B = dofvalues;
    FLOAT w0, *r;
    const FLOAT *p, *q, *qw, *qp, *b;

    Unused(pdofvalues);

    assert(type == EDGE);
    assert(funcvalues == NULL);
    assert(userfunc_lambda != NULL || userfunc != NULL);

    n = dof->type->np_edge;
    temp0 = phgAlloc(n * sizeof(*temp0));
    temp1 = phgAlloc(m * sizeof(*temp1));
    A = phgCalloc(n * n, sizeof(*A));

    for (i = 0; i < n * m; i++)
	B[i] = 0.0;

    GetEdgeVertices(e, index, u0, u1);
    x0 = (*(c = g->verts + e->verts[u0]))[0];
    y0 = (*c)[1];
    z0 = (*c)[2];
    x1 = (*(c = g->verts + e->verts[u1]))[0];
    y1 = (*c)[1];
    z1 = (*c)[2];
    x2 = x1 - x0;
    y2 = y1 - y0;
    z2 = z1 - z0;

    quad = phgQuadGetQuad1D(2 * dof->type->order);
    qp = quad->points;
    qw = quad->weights;
    /* range of basis functions on the current edge */
    n0 = n * index;
    n1 = n0 + n;
    for (k = 0; k < quad->npoints; k++) {
	lam[u0] = *(qp++);
	lam[u1] = *(qp++);
	if (userfunc_lambda != NULL) {
	    userfunc_lambda(dof, e, -1, lam, fvals);
	}
	else {
	    phgGeomLambda2XYZ(g, e, lam, &x, &y, &z);
	    userfunc(x, y, z, fvals);
	}
	/* inner product of fvals with the tangent vector */
	for (p = fvals, r = temp1, i = 0; i < m; i++, p += 3)
	    *(r++) = x2 * p[0] + y2 * p[1] + z2 * p[2];
	b = dof->type->BasFuncs(dof, e, n0, n1, lam);
	/* inner product of the basis funcs with the tangent vector */
	for (p = b, r = temp0, i = 0; i < n; i++, p += 3)
	    *(r++) = x2 * p[0] + y2 * p[1] + z2 * p[2];
	/* coefficients and RHS */
	w0 = *(qw++);
	for (i = 0, p = temp0; i < n; i++, p++) {
	    for (j = i, q = p; j < n; j++, q++)
		A[i*n+j] += p[0] * q[0] * w0;
	    for (j = 0, q = temp1; j < m; j++, q++)
		B[i*m+j] += p[0] * q[0] * w0;
	}
    }
    for (i = 1; i < n; i++)
	for (j = 0; j < i; j++)
	    A[i*n+j] = A[j*n+i];
    /* solver the linear system of equations */
    if (!phgSolverDenseSolver(n, m, (void *)A, (void *)B))
	phgError(1, "%s:%d, singular system (shouldn't happen).\n",
			__FILE__, __LINE__);
    phgFree(A);
    phgFree(temp0);
    phgFree(temp1);
#endif
    phgFree(fvals);
}

static const FLOAT *
HB2_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT values[3 * NEdge + 3 * NFace][Dim];
    int I, J, K, L;
    int k;
    GRID *g = dof->g;
    FLOAT (*nabla)[Dim + 1] = (void *)(phgGeomGetJacobian(g, e));

    CheckThread

    if (no1 <= 0)
	no1 = 3 * NEdge + 3 * NFace;

    for (k = 3*(no0/3); k < no1; k++) {
	if (k < 3 * NEdge) {
	    GetEdgeVertices(e, k / 3, I, J);
	    switch (k % 3) {
		case 0:
		    values[k][0] =
			lambda[J] * nabla[I][0] - lambda[I] * nabla[J][0];
		    values[k][1] =
			lambda[J] * nabla[I][1] - lambda[I] * nabla[J][1];
		    values[k][2] =
			lambda[J] * nabla[I][2] - lambda[I] * nabla[J][2];
		    break;
		case 1:
		    values[k][0] =
			lambda[J] * nabla[I][0] + lambda[I] * nabla[J][0];
		    values[k][1] =
			lambda[J] * nabla[I][1] + lambda[I] * nabla[J][1];
		    values[k][2] =
			lambda[J] * nabla[I][2] + lambda[I] * nabla[J][2];
		    break;
		case 2:
		    values[k][0] =
			values[k - 1][0] * (lambda[J] -
					    lambda[I]) * 1.5 -
			values[k - 2][0] * .5;
		    values[k][1] =
			values[k - 1][1] * (lambda[J] -
					    lambda[I]) * 1.5 -
			values[k - 2][1]  * .5;
		    values[k][2] =
			values[k - 1][2] * (lambda[J] -
					    lambda[I]) * 1.5 -
			values[k - 2][2] * .5;
		    break;
	    }
	}
	else {
	    GetFaceVertices(e, (k - 3 * NEdge) / 3, I, J, K, L);
	    switch ((k - 3 * NEdge) % 3) {
		case 0:
		    values[k][0] = lambda[I] * lambda[J] * nabla[K][0];
		    values[k][1] = lambda[I] * lambda[J] * nabla[K][1];
		    values[k][2] = lambda[I] * lambda[J] * nabla[K][2];
		    break;
		case 1:
		    values[k][0] = lambda[I] * lambda[K] * nabla[J][0];
		    values[k][1] = lambda[I] * lambda[K] * nabla[J][1];
		    values[k][2] = lambda[I] * lambda[K] * nabla[J][2];
		    break;
		case 2:
		    values[k][0] = lambda[K] * lambda[J] * nabla[I][0];
		    values[k][1] = lambda[K] * lambda[J] * nabla[I][1];
		    values[k][2] = lambda[K] * lambda[J] * nabla[I][2];
		    break;
	    }
	}
    }

    return (FLOAT *)(values + no0);
}

static const FLOAT *
HB2_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT values[3 * NEdge + 3 * NFace][Dim][Dim + 1];
    int I, J, K, L;
    int k;
    GRID *g = dof->g;
    FLOAT (*nabla)[Dim + 1] = (void *)(phgGeomGetJacobian(g, e));

    CheckThread

    if (no1 <= 0)
	no1 = 3 * NEdge + 3 * NFace;

    bzero(values, sizeof(values));

    for (k = no0; k < no1; k++) {
	if (k < 3 * NEdge) {
	    GetEdgeVertices(e, k / 3, I, J);
	    switch (k % 3) {
		case 0:
		    values[k][0][I] = -nabla[J][0];
		    values[k][0][J] = nabla[I][0];

		    values[k][1][I] = -nabla[J][1];
		    values[k][1][J] = nabla[I][1];

		    values[k][2][I] = -nabla[J][2];
		    values[k][2][J] = nabla[I][2];
		    break;
		case 1:
		    values[k][0][I] = nabla[J][0];
		    values[k][0][J] = nabla[I][0];

		    values[k][1][I] = nabla[J][1];
		    values[k][1][J] = nabla[I][1];

		    values[k][2][I] = nabla[J][2];
		    values[k][2][J] = nabla[I][2];
		    break;
		case 2:
		  /*printf("no:%d  e->index, %d I %d J %d\n",no,e->index,I,J);*/
		    values[k][0][I] =
			-3.0 * (lambda[J] * nabla[I][0] +
				lambda[I] * nabla[J][0]) / 2.0 +
			3.0 * (lambda[J] -
			       lambda[I]) * nabla[J][0] / 2.0 +
			nabla[J][0] / 2.0;
		    values[k][0][J] =
			3.0 * (lambda[J] * nabla[I][0] +
			       lambda[I] * nabla[J][0]) / 2.0 +
			3.0 * (lambda[J] -
			       lambda[I]) * nabla[I][0] / 2.0 -
			nabla[I][0] / 2.0;

		    values[k][1][I] =
			-3.0 * (lambda[J] * nabla[I][1] +
				lambda[I] * nabla[J][1]) / 2.0 +
			3.0 * (lambda[J] -
			       lambda[I]) * nabla[J][1] / 2.0 +
			nabla[J][1] / 2.0;
		    values[k][1][J] =
			3.0 * (lambda[J] * nabla[I][1] +
			       lambda[I] * nabla[J][1]) / 2.0 +
			3.0 * (lambda[J] -
			       lambda[I]) * nabla[I][1] / 2.0 -
			nabla[I][1] / 2.0;

		    values[k][2][I] =
			-3.0 * (lambda[J] * nabla[I][2] +
				lambda[I] * nabla[J][2]) / 2.0 +
			3.0 * (lambda[J] -
			       lambda[I]) * nabla[J][2] / 2.0 +
			nabla[J][2] / 2.0;
		    values[k][2][J] =
			3.0 * (lambda[J] * nabla[I][2] +
			       lambda[I] * nabla[J][2]) / 2.0 +
			3.0 * (lambda[J] -
			       lambda[I]) * nabla[I][2] / 2.0 -
			nabla[I][2] / 2.0;
		    break;
	    }
	}
	else {
	    GetFaceVertices(e, (k - 3 * NEdge) / 3, I, J, K, L);
	    switch ((k - 3 * NEdge) % 3) {
		case 0:
		    values[k][0][I] = lambda[J] * nabla[K][0];
		    values[k][0][J] = lambda[I] * nabla[K][0];

		    values[k][1][I] = lambda[J] * nabla[K][1];
		    values[k][1][J] = lambda[I] * nabla[K][1];

		    values[k][2][I] = lambda[J] * nabla[K][2];
		    values[k][2][J] = lambda[I] * nabla[K][2];
		    break;
		case 1:
		    values[k][0][I] = lambda[K] * nabla[J][0];
		    values[k][0][K] = lambda[I] * nabla[J][0];

		    values[k][1][I] = lambda[K] * nabla[J][1];
		    values[k][1][K] = lambda[I] * nabla[J][1];

		    values[k][2][I] = lambda[K] * nabla[J][2];
		    values[k][2][K] = lambda[I] * nabla[J][2];
		    break;
		case 2:
		    values[k][0][J] = lambda[K] * nabla[I][0];
		    values[k][0][K] = lambda[J] * nabla[I][0];

		    values[k][1][J] = lambda[K] * nabla[I][1];
		    values[k][1][K] = lambda[J] * nabla[I][1];

		    values[k][2][J] = lambda[K] * nabla[I][2];
		    values[k][2][K] = lambda[J] * nabla[I][2];
		    break;
	    }
	}
    }

    return (FLOAT *)(values + no0);
}

static void
HB2_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	 DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	 const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
{
    GRID *g = dof->g;
    FLOAT lam[] = {0., 0., 0., 0.};
    int i, j, k, n, u0, u1, u2, u3, n0, n1, eno, m = dof->dim;
    FLOAT x0, y0, z0, x1, y1, z1, x2, y2, z2, x, y, z;
    COORD *c;
    QUAD *quad;
    FLOAT *fvals = phgAlloc(Dim * m * sizeof(*fvals));
    FLOAT *temp, *A, *B = dofvalues;
    FLOAT w0, *r;
    const FLOAT *p, *q, *qw, *qp, *b, *nv;

    if (type == EDGE) {
	HB1_init(dof, e, type, index, userfunc, userfunc_lambda, funcvalues,
			dofvalues, pdofvalues);
	return;
    }

    assert(type == FACE);
    assert(userfunc_lambda != NULL || userfunc != NULL);
    assert(funcvalues == NULL);

    n = dof->type->np_face;
    temp = phgAlloc(Dim * n * sizeof(*temp));
    A = phgCalloc(n * n, sizeof(*A));

    for (i = 0; i < n * m; i++)
	B[i] = 0.0;

    /* range of basis functions on the current face */
    n0 = dof->type->np_edge * NEdge + n * index;
    n1 = n0 + n;

    GetFaceVertices(e, index, u0, u1, u2, u3);
    if (IsLeaf(e)) {
	nv = phgGeomGetFaceNormal(g, e, index);
    }
    else {
	static FLOAT n[3];
	x0 = (*(c = g->verts + e->verts[u0]))[0];
	y0 = (*c)[1];
	z0 = (*c)[2];
	x1 = (*(c = g->verts + e->verts[u1]))[0] - x0;
	y1 = (*c)[1] - y0;
	z1 = (*c)[2] - z0;
	x2 = (*(c = g->verts + e->verts[u2]))[0] - x0;
	y2 = (*c)[1] - y0;
	z2 = (*c)[2] - z0;
	n[0] = y1 * z2 - y2 * z1;
	n[1] = z1 * x2 - z2 * x1;
	n[2] = x1 * y2 - x2 * y1;
	nv = n;
    }
    quad = phgQuadGetQuad2D(2 * dof->type->order);
    qp = quad->points;
    qw = quad->weights;
    for (k = 0; k < quad->npoints; k++) {
	lam[u0] = *(qp++);
	lam[u1] = *(qp++);
	lam[u2] = *(qp++);
	if (userfunc_lambda != NULL) {
	    userfunc_lambda(dof, e, -1, lam, fvals);
	}
	else {
	    phgGeomLambda2XYZ(g, e, lam, &x, &y, &z);
	    userfunc(x, y, z, fvals);
	}
	/* evaluate basis functions, since we also need values of the
	 * basis functions on all edges, the range is set to
	 * (0, n1) */
	b = dof->type->BasFuncs(dof, e, 0, n1, lam);
	/* fvals -= contributions from edges. Note that at this stage
	 * dof->data already contains valid data on all edges */
	for (eno = 0, p = b; eno < NEdge; eno++) {
	    /* Note: only need to remove contributions of the
	     * basis functions on the three edges of the face */
	    if (eno != GetEdgeNo(u0,u1) && eno != GetEdgeNo(u0,u2) &&
		eno != GetEdgeNo(u1,u2)) {
		p += dof->type->np_edge * 3;
		continue;
	    }
	    if (pdofvalues == NULL)
		q = DofEdgeData(dof, e->edges[eno]);
	    else
		q = pdofvalues[NVert + eno];
	    for (i = 0; i < dof->type->np_edge; i++, p += 3) {
		for (j = 0, r = fvals; j < m; j++) {
		    w0 = *(q++);
		    *(r++) -= w0 * p[0];
		    *(r++) -= w0 * p[1];
		    *(r++) -= w0 * p[2];
		}
	    }
	}
	/* cross product of fvals with the face normal */
	for (i = 0, r = fvals; i < m; i++) {
	    x = r[1] * nv[2] - r[2] * nv[1];
	    y = r[2] * nv[0] - r[0] * nv[2];
	    z = r[0] * nv[1] - r[1] * nv[0];
	    (*r++) = x;
	    (*r++) = y;
	    (*r++) = z;
	}
	b += 3 * n0;
	/* cross product of basis functions with the face normal */
	for (i = 0, r = temp, p = b; i < n; i++, p += 3) {
	    *(r++) = p[1] * nv[2] - p[2] * nv[1];
	    *(r++) = p[2] * nv[0] - p[0] * nv[2];
	    *(r++) = p[0] * nv[1] - p[1] * nv[0];
	}
	/* coefficients and RHS */
	w0 = *(qw++);
	for (i = 0, p = temp; i < n; i++, p += 3) {
	    for (j = i, q = p; j < n; j++, q += 3)
		A[i*n+j] += (p[0]*q[0] + p[1]*q[1] + p[2]*q[2]) * w0;
	    for (j = 0, q = fvals; j < m; j++, q += 3)
		B[i*m+j] += (p[0]*q[0] + p[1]*q[1] + p[2]*q[2]) * w0;
	}
    }

    for (i = 1; i < n; i++)
	for (j = 0; j < i; j++)
	    A[i*n+j] = A[j*n+i];
    /* solver the linear system of equations */
    if (!phgSolverDenseSolver(n, m, (void *)A, (void *)B))
	phgError(1, "%s:%d, singular system (shouldn't happen).\n",
			__FILE__, __LINE__);

    phgFree(A);
    phgFree(temp);
    phgFree(fvals);
}

static BYTE HB1_orders[] = {1, 1};
DOF_TYPE DOF_HFEB1_ = {
    DofReserved,
    "HFEB1",
    NULL,	/* points */
    HB1_orders,
    DOF_DG0,
    NULL,	/* base_type */
    NULL,	/* hp_info */
    phgDofInterpC2FGeneric,
    phgDofInterpF2CGeneric,
    HB1_init,
    HB1_bas,
    HB1_grad,
    NULL,
    FE_Hcurl,
    FALSE,
    FALSE,
    FALSE,
    -1,
    2 * NEdge,
    1,
    0,	/* D_order */
    -1,
    Dim,
    0,
    2,
    0,
    0
};

static BYTE HB2_orders[] = {1, 1, 2, 2, 2, 2};
DOF_TYPE DOF_HFEB2_ = {
    DofReserved,
    "HFEB2",
    NULL,	/* points */
    HB2_orders,
    DOF_DG1,
    NULL,	/* base_type */
    NULL,	/* hp_info */
    phgDofInterpC2FGeneric,
    phgDofInterpF2CGeneric,
    HB2_init,
    HB2_bas,
    HB2_grad,
    NULL,
    FE_Hcurl,
    FALSE,
    FALSE,
    FALSE,
    -1,
    3 * (NEdge + NFace),
    2,
    0,	/* D_order */
    -1,
    Dim,
    0,
    3,
    3,
    0
};

DOF_TYPE *DOF_HFEBn[] = {DOF_HC0, DOF_HFEB1, DOF_HFEB2, NULL};
