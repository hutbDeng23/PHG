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

/* $Id: hierarchical-basis.c,v 1.141 2022/08/19 09:54:51 zlb Exp $
 *
 * This file implements general p-th order H1 conforming and Hcurl conforming
 * hierarchical basis functions */

#include "phg.h"

#include <stdlib.h>
#include <string.h>
#include <strings.h>	/* for bzero */
#include <math.h>

/* Note: use 'make USER_CPPFLAGS=-DSHOW_CONDITION lib' to show the
 * condition nunmbers of local mass matrices used in the projection. */

#if defined(SHOW_CONDITION)

static double
condition_number(int n, const FLOAT a0[])
{
#if !USE_LAPACK
    phgPrintf("*** WARNING: no LAPACK, condition number not computed.\n");
    return 0.0;
#else
    double *a, *work, anorm, rcond;
    int *iwork, i;
    extern void F77_FUNC(dgecon,DGECON)(const char *, int *, double *, int *,
			double *, double *, double *, int *, int *, int);
    extern void F77_FUNC(dgetrf,DGETRF)(int *, int *, double *, int *,
			int *, int *);

    a = phgAlloc(n * (size_t)n * sizeof(*a));
    work = phgAlloc(4 * n * sizeof(*work));
    iwork = phgAlloc(n * sizeof(*iwork));

    anorm = rcond = 0.0;
    for (i = 0; i < n * n; i++) {
	a[i] = (double)(a0[i]);
		/*a[i] = (i % (n+1) == 0 ? 2.:0.);*/
	rcond += fabs(a[i]);
	if ((i + 1) % n == 0) {
	    if (anorm < rcond)
		anorm = rcond;
	    rcond = 0.0;
	}
    }

    F77_FUNC(dgetrf,DGETRF)(&n, &n, a, &n, iwork, &i);
    F77_FUNC(dgecon,DGECON)
	("1", &n, a, &n, &anorm, &rcond, work, iwork, &i, 1);
    if (i != 0)
	printf("Error in DGECON, info = %d\n", i);

    phgFree(a);
    phgFree(work);
    phgFree(iwork);

    return rcond;
#endif
}

static void
show_condition(int n, const FLOAT a[], GTYPE type, int index)
{
    static int edge_flag = 0, face_flag = 0, elem_flag = 0;
    int *flag;
    const char *name;

    if (type == EDGE) {
	flag = &edge_flag;
	name = "Edge";
    }
    else if (type == FACE) {
	flag = &face_flag;
	name = "Face";
    }
    else {
	assert(type == VOLUME);
	flag = &elem_flag;
	name = "Element";
    }
    if (*flag == 0) {
	*flag = 1;
	if (type != VOLUME)
	    phgPrintf("*** rcond(A) on %s %d: %e (n = %d)\n",
				name, index, condition_number(n, a), n);
	else
	    phgPrintf("*** rcond(A) on %s: %e (n = %d)\n",
				name, condition_number(n, a), n);
    }
}

#else	/* defined(SHOW_CONDITION) */
# define show_condition(n, a, type, index)
#endif	/* defined(SHOW_CONDITION) */

/* Evaluating Legendre polynomials and derivatives using the recurrence
 * relation (P_l being the Legendre polynomial of order l):
 *
 *	P_0 = 1,
 *	P_1 = x,
 *	(l+1) * P_{l+1} - (2l+1) * x * P_l + l * P_{l-1} = 0,
 *	l = 1, 2, ...
 */

void
phgLegendreP(int p, FLOAT x, FLOAT *values, FLOAT *values_d)
/* evaluates Legendre polynomials and derivatives up to order p at x in [-1,1]:
 *	 values[k] = L_k(x),
 * and if values_d != NULL:
 *	values_d[k] = L'_k(x),
 * for k = 0, ..., p
 */
{
    int i;
    FLOAT L0, L1, L2, DL0, DL1, DL2;

    assert(p >= 0);

    if (p == 0) {
	*values = 1.0;
	if (values_d != NULL)
	    *values_d = 0.0;
	return;
    }
    else if (p == 1) {
	*(values++) = 1.0;
	*values = x;
	if (values_d != NULL) {
	    *(values_d++) = 0.0;
	    *values_d = 1.0;
	}
	return;
    }

    if (values_d == NULL) {
	*(values++) = L0 = 1.0;
	*(values++) = L1 = x;
	for (i = 1; i < p; i++) {
	    L2 = ((i + i + 1) * x * L1 - i * L0) / (FLOAT)(i + 1);
	    L0 = L1;
	    L1 = *(values++) = L2;
	}
    }
    else {
	*(values++) = L0 = 1.0;
	*(values++) = L1 = x;
	*(values_d++) = DL0 = 0.0;
	*(values_d++) = DL1 = 1.0;
	for (i = 1; i < p; i++) {
	    DL2 = ((i + i + 1) * (L1 + x * DL1) - i * DL0) / (FLOAT)(i + 1);
	    DL0 = DL1;
	    DL1 = *(values_d++) = DL2;
	    L2 = ((i + i + 1) * x * L1 - i * L0) / (FLOAT)(i + 1);
	    L0 = L1;
	    L1 = *(values++) = L2;
	}
    }

    return;
}

#if 0
/* testing evaluation of Legendre polynomials and derivatives */
static FLOAT
L10(FLOAT x)
{
    return -63 / (FLOAT)256.
	   +(3465. * Pow(x, 2)) / 256.
	   -(15015. * Pow(x, 4)) / 128.
	   +(45045. * Pow(x, 6)) / 128.
	   -(109395. * Pow(x, 8)) / 256.
	   +(46189. * Pow(x, 10)) / 256.;
}

static FLOAT
DL10(FLOAT x)
{
    return +(3465. * x) / 128.
	   -(15015. * Pow(x, 3)) / 32. 
	   +(135135. * Pow(x, 5)) / 64.
	   -(109395. * Pow(x, 7)) / 32.
	   +(230945. * Pow(x, 9)) / 128.;
}

int
main(int argc, char **argv[])
{
    int i, n = 10000;
    FLOAT x, d = 2.0 / n, L[11], DL[11];

    for (i = 0; i <= n; i++) {
	x = -1. + i * d;
	phgLegendreP(10, x, L, DL);
	/* Plot 1:2 and 1:3 for function/derivative, 1:4 and 1:5 for errors */
	printf("%lg %lg %lg %lg %lg\n",
		(double)x, (double)L[10], (double)DL[10],
		(double)(L[10] - L10(x)), (double)(DL[10] - DL10(x)));
    }

    return 0;
}
#endif	/* 0 */

/*-------------------------- Supporting functions --------------------------*/

/* buffers for caching evaluated Legendre polynomials. */

/* TODO: make use of symmetry of the basis functions and evaluation
 *	 points (mostly come from quadrature rules) to save more evaluations */

#define CACHE_LP 1	/* !0 ==> cache evaluated Legendre polynomials */

static FLOAT (*L0)[ORDER_MAX] = NULL, (*DL0)[ORDER_MAX] = NULL;
#if USE_OMP
#pragma omp threadprivate(L0, DL0)
#endif	/* USE_OMP */

#if CACHE_LP
static BYTE *cache_order = NULL;
static BOOLEAN cache_reset = TRUE;
#if USE_OMP
#pragma omp threadprivate(cache_order, cache_reset)
#endif	/* USE_OMP */

static inline void
get_L(int order, int v0, int v1, const FLOAT lambda[], FLOAT **L, FLOAT **DL)
/* evaluate and cache Legendre polynomials on a given edge */
{
    int no = GetEdgeNo(v0, v1);

    assert(order >= 0);
    assert(order <= ORDER_MAX - 1);

    if (cache_order[no] < order + 1) {
	cache_order[no] = order + 1;
	phgLegendreP(order, lambda[v1] - lambda[v0],
				L0[no], DL == NULL ? NULL : DL0[no]);
    }

    *L = L0[no];
    if (DL != NULL)
	*DL = DL0[no];
}
#endif	/* CACHE_LP */

/*-------------------------- H1 conforming bases --------------------------*/

static const FLOAT *
HB_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    /* Static buffers for returning values of basis functions.
     * 
     * Note: For uniform order bases, one buffer is used for each order to
     *	     reduce possible buffer conflicts.
     *
     *	     For hp bases, the buffer used is determined by the highest order
     *	     in the element. */
    static FLOAT *buffers[ORDER_MAX + 1];
    static BOOLEAN initialized = FALSE;
#if USE_OMP
#pragma omp threadprivate(buffers, initialized)
#endif	/* USE_OMP */
    FLOAT *L1, *L2, *L3;

    DOF_TYPE *t;
    int order0, order, nbas, bufferlen;
    FLOAT *bas, *p, *p0, beta;
    int i, j, k, i0, v0, v1, v2, v3, v4, v5, v6, v7, m, n;

    if (!initialized) {
	initialized = TRUE;
	bzero(buffers, sizeof(buffers));
    }

    t = dof->type;
    if (t != NULL) {
	order = order0 = t->order;
	bufferlen = nbas = t->nbas;
    }
    else {
	assert(DofIsHP(dof));
	/* FIXME: it is assumed that vert/edge/face orders <= elem order */
	order = order0 = dof->hp->elem_order[e->index];
	nbas = DofNBas(dof, e);
	bufferlen = dof->hp->info->types[order0]->nbas;
    }

    assert(order0 >= 1 && order0 <= ORDER_MAX);

    if (no1 <= 0)
	no1 = nbas;

    if (buffers[order0] == NULL) {
	buffers[order0] = phgAlloc(bufferlen * sizeof(*bas));
	FreeAtExitNoCheck(buffers[order0]);
    }
    bas = buffers[order0];

    p = bas;
    p0 = p + (no1 - no0);

    /*------------------------- vertex functions --------------------------*/
    for (i = no0; i < NVert && i < no1; i++)
	*(p++) = lambda[i];

    if (i >= no1)
	return bas;

    if ((no0 -= NVert) < 0)
	no0 = 0;

#if CACHE_LP
    if (L0 == NULL) {
	L0 = phgAlloc(NEdge * sizeof(*L0));
	DL0 = phgAlloc(NEdge * sizeof(*DL0));
	cache_order = phgAlloc(NEdge * sizeof(*cache_order));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
	FreeAtExitNoCheck(cache_order);
    }
    if (cache_reset)
	memset(cache_order, 0, NEdge * sizeof(*cache_order));
#else	/* CACHE_LP */
    if (L0 == NULL) {
	L0 = phgAlloc(3 * sizeof(*L0));
	DL0 = phgAlloc(3 * sizeof(*DL0));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
    }
#endif	/* CACHE_LP */

    /*------------------------- edge functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NEdge && p < p0; i++) {
	if (DofIsHP(dof)) {
	    order = dof->hp->edge_order[e->edges[i]];
	    assert(order <= order0);
	    t = dof->hp->info->types[order];
	}
	if (t->np_edge == 0)
	    continue;
	n = i0 + t->np_edge;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	GetEdgeVertices(e, i, v0, v1);
	beta = lambda[v0] * lambda[v1];
#if CACHE_LP
	get_L(order0 - 2, v0, v1, lambda, &L1, NULL);
#else	/* CACHE_LP */
	phgLegendreP(order - 2, lambda[v1] - lambda[v0], L1 = L0[0], NULL);
#endif	/* CACHE_LP */
	for (j = 0; i0 < n && p < p0; i0++, j++) {
	    if (i0 < no0)
		continue;
	    *(p++) = beta * L1[j];
	}
    }

    if (p >= p0)
	return bas;

    if ((no0 -= i0) < 0)
	no0 = 0;

    /*------------------------- face functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NFace && p < p0; i++) {
	if (DofIsHP(dof)) {
	    order = dof->hp->face_order[e->faces[i]];
	    assert(order <= order0);
	    t = dof->hp->info->types[order];
	}
	if (t->np_face == 0)
	    continue;
	n = i0 + t->np_face;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	GetFaceVertices(e, i, v0, v1, v2, v3);
	beta = lambda[v0] * lambda[v1] * lambda[v2];
#if CACHE_LP
	get_L(order0 - 3, v0, v1, lambda, &L1, NULL);
	get_L(order0 - 3, v0, v2, lambda, &L2, NULL);
#else	/* CACHE_LP */
	phgLegendreP(order - 3, lambda[v1] - lambda[v0], L1 = L0[0], NULL);
	phgLegendreP(order - 3, lambda[v2] - lambda[v0], L2 = L0[1], NULL);
#endif	/* CACHE_LP */
	for (m = 0; m <= order - 3 && p < p0; m++) {
	    if (i0 + m + 1 <= no0) {
		i0 += m + 1;
		continue;
	    }
	    for (j = 0; j <= m && p < p0; j++, i0++) {
		if (i0 < no0)
		    continue;
		*(p++) = beta * L1[j] * L2[m - j];
	    }
	}
    }

    if (p >= p0)
	return bas;

    if ((no0 -= i0) < 0)
	no0 = 0;

    /*------------------------- element functions --------------------------*/
    if (DofIsHP(dof)) {
	order = order0;
	t = dof->hp->info->types[order];
    }
    n = t->np_elem;
    assert(no0 < n);
    assert(order >= 4);
    GetElementVertices(e, 0, v0, v1, v2, v3, 
		       v4, v5, v6, v7);
    beta = lambda[v0] * lambda[v1] * lambda[v2] * lambda[v3];
#if CACHE_LP
    get_L(order0 - 4, v0, v1, lambda, &L1, NULL);
    get_L(order0 - 4, v0, v2, lambda, &L2, NULL);
    get_L(order0 - 4, v0, v3, lambda, &L3, NULL);
#else	/* CACHE_LP */
    phgLegendreP(order - 4, lambda[v1] - lambda[v0], L1 = L0[0], NULL);
    phgLegendreP(order - 4, lambda[v2] - lambda[v0], L2 = L0[1], NULL);
    phgLegendreP(order - 4, lambda[v3] - lambda[v0], L3 = L0[2], NULL);
#endif	/* CACHE_LP */
    i0 = 0;
    for (m = 0; m <= order - 4 && p < p0; m++) {
	if (i0 + (m + 1) * (m + 2) / 2 <= no0) {
	    i0 += (m + 1) * (m + 2) / 2;
	    continue;
	}
	for (k = 0; k <= m && p < p0; k++) {
	    if (i0 + m - k + 1 <= no0) {
		i0 += m - k + 1;
		continue;
	    }
	    for (j = 0; j <= m - k && p < p0; j++, i0++) {
		if (i0 < no0)
		    continue;
		*(p++) = beta * L1[j] * L2[k] * L3[m - j - k];
	    }
	}
    }

    if (p != p0)
	phgError(1, "(%s:%d) HB_bas: p=%d, p0=%d\n", __FILE__, __LINE__,
			(int)(p - bas), (int)(p0 - bas));

    return bas;
}

static const FLOAT *
HB_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT *buffers[ORDER_MAX + 1];
    static BOOLEAN initialized = FALSE;
#if USE_OMP
#pragma omp threadprivate(buffers, initialized)
#endif	/* USE_OMP */
    FLOAT *L1, *L2, *L3, *DL1, *DL2, *DL3;

    DOF_TYPE *t;
    int order, order0, nbas, bufferlen;
    FLOAT *grad, *p, *p0, beta, beta0, beta1, beta2, beta3, a, b, c, d;
    int i, j, k, ii, i0, v0, v1, v2, v3, v4, v5, v6, v7, m, n;

    if (!initialized) {
	initialized = TRUE;
	bzero(buffers, sizeof(buffers));
    }

    t = dof->type;
    if (t != NULL) {
	order = order0 = t->order;
	bufferlen = (nbas = t->nbas) * (Dim + 1);
    }
    else {
	assert(DofIsHP(dof));
	/* FIXME: it is assumed that vert/edge/face orders <= elem order */
	order = order0 = dof->hp->elem_order[e->index];
	nbas = DofNBas(dof, e);
	bufferlen = dof->hp->info->types[order0]->nbas * (Dim + 1);
    }

    assert(order0 >= 1 && order0 <= ORDER_MAX);

    if (no1 <= 0)
	no1 = nbas;

    if ((grad = buffers[order]) == NULL) {
	buffers[order] = phgAlloc(bufferlen * sizeof(*grad));
	FreeAtExitNoCheck(buffers[order]);
	grad = buffers[order];
    }

    if (order == 0) {
	grad[0] = grad[1] = grad[2] = grad[3] = 0.0;
	return grad;
    }

    p = grad;
    p0 = p + (no1 - no0) * (Dim + 1);

    /*------------------------- vertex functions --------------------------*/
    for (i = no0; i < NVert && i < no1; i++)
	for (j = 0; j < Dim + 1; j++)
	    *(p++) = (i == j ? 1.0 : 0.0);

    if (i >= no1)
	return grad;

    if ((no0 -= NVert) < 0)
	no0 = 0;

#if CACHE_LP
    if (L0 == NULL) {
	L0 = phgAlloc(NEdge * sizeof(*L0));
	DL0 = phgAlloc(NEdge * sizeof(*DL0));
	cache_order = phgAlloc(NEdge * sizeof(*cache_order));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
	FreeAtExitNoCheck(cache_order);
    }
    if (cache_reset)
	memset(cache_order, 0, NEdge * sizeof(*cache_order));
#else	/* CACHE_LP */
    if (L0 == NULL) {
	L0 = phgAlloc(3 * sizeof(*L0));
	DL0 = phgAlloc(3 * sizeof(*DL0));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
    }
#endif	/* CACHE_LP */

    /*------------------------- edge functions --------------------------*/
    i0 = 0;
    for (i = i0; i < NEdge && p < p0; i++) {
	if (DofIsHP(dof)) {
	    order = dof->hp->edge_order[e->edges[i]];
	    assert(order <= order0);
	    t = dof->hp->info->types[order];
	}
	if (t->np_edge == 0)
	    continue;
	n = i0 + t->np_edge;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	GetEdgeVertices(e, i, v0, v1);
	beta = lambda[v0] * lambda[v1];
#if CACHE_LP
	get_L(order0 - 2, v0, v1, lambda, &L1, &DL1);
#else	/* CACHE_LP */
	phgLegendreP(order - 2, lambda[v1] - lambda[v0],
					L1 = L0[0], DL1 = DL0[0]);
#endif	/* CACHE_LP */
	for (j = 0; i0 < n && p < p0; i0++, j++) {
	    if (i0 < no0)
		continue;
	    a = beta * DL1[j];
	    for (ii = 0; ii < Dim + 1; ii++) {
		if (ii == v0)
		    *(p++) = lambda[v1] * L1[j] - a;
		else if (ii == v1)
		    *(p++) = lambda[v0] * L1[j] + a;
		else
		    *(p++) = 0.0;
	    }
	}
    }

    if (p >= p0)
	return grad;

    if ((no0 -= i0) < 0)
	no0 = 0;

    /*------------------------- face functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NFace && p < p0; i++) {
	if (DofIsHP(dof)) {
	    order = dof->hp->face_order[e->faces[i]];
	    assert(order <= order0);
	    t = dof->hp->info->types[order];
	}
	if (t->np_face == 0)
	    continue;
	n = i0 + t->np_face;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	GetFaceVertices(e, i, v0, v1, v2, v3);
	beta0 = lambda[v1] * lambda[v2];
	beta1 = lambda[v0] * lambda[v2];
	beta2 = lambda[v0] * lambda[v1];
	beta = lambda[v0] * beta0;
#if CACHE_LP
	get_L(order0 - 3, v0, v1, lambda, &L1, &DL1);
	get_L(order0 - 3, v0, v2, lambda, &L2, &DL2);
#else	/* CACHE_LP */
	phgLegendreP(order - 3, lambda[v1] - lambda[v0],
					L1 = L0[0], DL1 = DL0[0]);
	phgLegendreP(order - 3, lambda[v2] - lambda[v0],
					L2 = L0[1], DL2 = DL0[1]);
#endif	/* CACHE_LP */
	for (m = 0; m <= order - 3 && p < p0; m++) {
	    if (i0 + m + 1 <= no0) {
		i0 += m + 1;
		continue;
	    }
	    for (j = 0; j <= m && p < p0; j++, i0++) {
		if (i0 < no0)
		    continue;
		a = L1[j] * L2[m - j];
		b = beta * DL1[j] * L2 [m - j];
		c = beta * L1 [j] * DL2[m - j];
		for (ii = 0; ii < Dim + 1; ii++) {
		    if (ii == v0)
			*(p++) = beta0 * a - (b + c);
		    else if (ii == v1)
			*(p++) = beta1 * a + b;
		    else if (ii == v2)
			*(p++) = beta2 * a + c;
		    else
			*(p++) = 0.;
		}
	    }
	}
    }

    if (p >= p0)
	return grad;

    if ((no0 -= i0) < 0)
	no0 = 0;

    /*------------------------- element functions --------------------------*/
    if (DofIsHP(dof)) {
	order = order0;
	t = dof->hp->info->types[order];
    }
    n = t->np_elem;
    assert(no0 < n);
    assert(order >= 4);
    GetElementVertices(e, 0, v0, v1, v2, v3,
		       v4, v5, v6, v7);
    beta0 = lambda[v1] * lambda[v2] * lambda[v3];
    beta1 = lambda[v0] * lambda[v2] * lambda[v3];
    beta2 = lambda[v0] * lambda[v1] * lambda[v3];
    beta3 = lambda[v0] * lambda[v1] * lambda[v2];
    beta = lambda[v0] * beta0;
#if CACHE_LP
    get_L(order0 - 4, v0, v1, lambda, &L1, &DL1);
    get_L(order0 - 4, v0, v2, lambda, &L2, &DL2);
    get_L(order0 - 4, v0, v3, lambda, &L3, &DL3);
#else	/* CACHE_LP */
    phgLegendreP(order - 4, lambda[v1] - lambda[v0], L1 = L0[0], DL1 = DL0[0]);
    phgLegendreP(order - 4, lambda[v2] - lambda[v0], L2 = L0[1], DL2 = DL0[1]);
    phgLegendreP(order - 4, lambda[v3] - lambda[v0], L3 = L0[2], DL3 = DL0[2]);
#endif	/* CACHE_LP */
    i0 = 0;
    for (m = 0; m <= order - 4 && p < p0; m++) {
	if (i0 + (m + 1) * (m + 2) / 2 <= no0) {
	    i0 += (m + 1) * (m + 2) / 2;
	    continue;
	}
	for (k = 0; k <= m && p < p0; k++) {
	    if (i0 + m - k + 1 <= no0) {
		i0 += m - k + 1;
		continue;
	    }
	    for (j = 0; j <= m - k && p < p0; j++, i0++) {
		if (i0 < no0)
		    continue;
		a = L1[j] * L2[k] * L3[m - j - k];
		b = beta * DL1[j] * L2 [k] * L3 [m - j - k];
		c = beta * L1 [j] * DL2[k] * L3 [m - j - k];
		d = beta * L1 [j] * L2 [k] * DL3[m - j - k];
		p[v0] = beta0 * a - (b + c + d);
		p[v1] = beta1 * a + b;
		p[v2] = beta2 * a + c;
		p[v3] = beta3 * a + d;
		p += Dim + 1;
	    }
	}
    }

    assert(p == p0);

    return grad;
}

static DOF *dof0;

static void
HB_init0(DOF *dof, ELEMENT *e, GTYPE type, int index,
	DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
{
    DOF_TYPE *t;
    int order;
    FLOAT a, b, w, x, y, z, *mat, lambda[Dim + 1] = {0., 0., 0., 0.};
    /* pointers to (known) DOF values at lower dimensional entities */
    const FLOAT *f0, *f1, *f2, *f3,		/* vertices */
		*g0, *g1, *g2, *g3, *g4, *g5,	/* edges */
		*h0, *h1, *h2, *h3;		/* faces */
    FLOAT *f = NULL, *qp, *qw;
    const FLOAT *bas, *p;
    int i, j, k, n, m, no0, no1, ii, np2, np3, n0, *pvt = NULL;
    int v0, v1, v2, v3, v4, v5, v6, v7, i0, i1, i2, i3, i4, i5,
	pe[NEdge + 1], pf[NFace + 1];
    QUAD *quad;

    /* cache for factorized matrices */
    static FLOAT *LU_cache[3 * (ORDER_MAX + 1)];
    static int *PV_cache[3 * (ORDER_MAX + 1)];
    static BOOLEAN initialized = FALSE;
#if USE_OMP
#pragma omp threadprivate(LU_cache, PV_cache, initialized)
#endif  /* USE_OMP */

    if (!initialized) {
	initialized = TRUE;
	bzero(LU_cache, sizeof(LU_cache));
	bzero(PV_cache, sizeof(PV_cache));
    }

    assert(userfunc != NULL || userfunc_lambda != NULL);

    switch (type) {
	case VERTEX:
	    /* assign function values to DOF */
	    assert(index >= 0 && index < NVert);
	    t = dof->type;
	    if (t == NULL)
		t = dof->hp->info->types[dof->hp->vert_order[e->verts[index]]];
	    if ((n = t->np_vert) == 0)
		break;
	    assert(t->np_vert == 1);
	    order = t->order;
	    lambda[index] = 1.0;
	    if (userfunc_lambda != NULL) {
		userfunc_lambda(dof0, e, 0, lambda, dofvalues);
	    }
	    else {
		phgGeomLambda2XYZ(dof->g, e, lambda, &x, &y, &z);
		userfunc(x, y, z, dofvalues);
	    }
	    break;
	case EDGE:
	    /* compute DOF by L2-projection on the edge */
	    assert(index >= 0 && index < NEdge);
	    t = dof->type;
	    if (t == NULL)
		t = dof->hp->info->types[dof->hp->edge_order[e->edges[index]]];
	    if ((n = t->np_edge) == 0)
		break;
	    order = t->order;
	    assert(order >= 1 && order <= ORDER_MAX);
	    GetEdgeVertices(e, index, v0, v1)
	    m = dof->dim;
	    f = phgAlloc(m * sizeof(*f));
	    if (pdofvalues == NULL) {
		f0 = DofVertexData(dof, e->verts[v0]);
		f1 = DofVertexData(dof, e->verts[v1]);
	    }
	    else {
		f0 = pdofvalues[v0];
		f1 = pdofvalues[v1];
	    }
	    /* compute the DOF by solving dof->dim systems of linear equations
	     * of size n, dof->data or pdofvalues[] should point to valid
	     * DOF values on the two vertices */
	    if (LU_cache[3 * order] == NULL) {
		LU_cache[3 * order] = phgCalloc(n * n, sizeof(*mat));
		FreeAtExitNoCheck(LU_cache[3 * order]);
		PV_cache[3 * order] = phgAlloc(n * sizeof(*pvt));
		FreeAtExitNoCheck(PV_cache[3 * order]);
		mat = LU_cache[3 * order];
		pvt = PV_cache[3 * order];
	    }
	    else {
		mat = NULL;
	    }
	    bzero(dofvalues, m * n * sizeof(*dofvalues));
	    quad = phgQuadGetQuad1D(2 * order);
	    qp = quad->points;
	    qw = quad->weights;
	    if (!DofIsHP(dof)) {
		no0 = NVert + index * n;
	    }
	    else {
		no0 = NVert;
		for (ii = 0; ii < index; ii++)
		    no0 += dof->hp->info->types
				[dof->hp->edge_order[e->edges[ii]]]->np_edge;
	    }
	    no1 = no0 + n;
#define A(i,j) mat[(i) * n + (j)]
#define B(i,j) dofvalues[(i) * m + (j)]
	    for (k = 0; k < quad->npoints; k++) {
		lambda[v0] = *(qp++);
		lambda[v1] = *(qp++);
		w = *(qw++);
		if (userfunc_lambda != NULL) {
		    userfunc_lambda(dof0, e, 0, lambda, f);
		}
		else {
		    phgGeomLambda2XYZ(dof->g, e, lambda, &x, &y, &z);
		    userfunc(x, y, z, f);
		}
		/* Note: must call BasFuncs after user_func since the latter
		 * might call BasFuncs, such as with phgDofInterpC2FGeneric */
		bas = t->BasFuncs(dof, e, no0, no1, lambda);
		for (i = 0; i < n; i++) {
		    b = bas[i] * w;
		    if (mat != NULL) {
			for (j = 0; j <= i; j++) {
			    A(i,j) += (a = b * bas[j]);
			    if (j < i)
				A(j,i) += a;
			}
		    }
		    for (j = 0; j < m; j++) {
			a = f[j] - (lambda[v0] * f0[j] + lambda[v1] * f1[j]);
			B(i,j) += b * a;
		    }
		}
	    }
#if 0
if (mat != NULL) {
printf("---------------------------------------------\n");
for (i = 0; i < n; i++) {
for (j = 0; j < n; j++) printf(" %0.4lf",
		(double)(Fabs(A(i,j))<1e-10 ? 0. : A(i,j)));
printf(" = %lg\n", (double)B(i,0));
}
printf("---------------------------------------------\n");
}
#endif
#undef A
#undef B
	    /* solve the equation A x = B */
	    if (mat != NULL) {
    		show_condition(n, mat, type, index);
		if (!phgSolverDenseLU(n, (void *)mat, pvt))
		    phgError(1, "(%s:%d) singular system!\n",__FILE__,__LINE__);
	    }
	    else {
		mat = LU_cache[3 * order];
		pvt = PV_cache[3 * order];
	    }
	    phgSolverDenseSV(n, (void *)mat, pvt, m, (void *)dofvalues);
	    break;
	case FACE:
	    /* compute DOF by L2-projection on the face */
	    assert(index >= 0 && index < NFace);
	    t = dof->type;
	    if (t == NULL)
		t = dof->hp->info->types[dof->hp->face_order[e->faces[index]]];
	    if ((n = t->np_face) == 0)
		break;
	    order = t->order;
	    assert(order >= 1 && order <= ORDER_MAX);
	    GetFaceVertices(e, index, v0, v1, v2, v3);
	    m = dof->dim;
	    f = phgAlloc(m * sizeof(*f));
	    /* the next 2 lines make gcc happy */
	    f0 = f1 = f2 = g0 = g1 = g2 = NULL;
	    i0 = i1 = i2 = 0;
	    if (pdofvalues == NULL) {
		f0 = DofVertexData(dof, e->verts[v0]);
		f1 = DofVertexData(dof, e->verts[v1]);
		f2 = DofVertexData(dof, e->verts[v2]);
		g0 = DofEdgeData(dof, e->edges[i0 = GetEdgeNo(v0, v1)]);
		g1 = DofEdgeData(dof, e->edges[i1 = GetEdgeNo(v0, v2)]);
		g2 = DofEdgeData(dof, e->edges[i2 = GetEdgeNo(v1, v2)]);
	    }
	    else {
		f0 = pdofvalues[v0];
		f1 = pdofvalues[v1];
		f2 = pdofvalues[v2];
		g0 = pdofvalues[NVert + (i0 = GetEdgeNo(v0, v1))];
		g1 = pdofvalues[NVert + (i1 = GetEdgeNo(v0, v2))];
		g2 = pdofvalues[NVert + (i2 = GetEdgeNo(v1, v2))];
	    }
	    /* compute the DOF by solving dof->dim systems of linear equations
	     * of size n, dof->data or pdofvalues[] should point to valid
	     * DOF values on the three vertices and three edges */
	    if (LU_cache[3 * order + 1] == NULL) {
		LU_cache[3 * order + 1] = phgCalloc(n * n, sizeof(*mat));
		FreeAtExitNoCheck(PV_cache[3 * order + 1]);
		PV_cache[3 * order + 1] = phgAlloc(n * sizeof(*pvt));
		FreeAtExitNoCheck(LU_cache[3 * order + 1]);
		mat = LU_cache[3 * order + 1];
		pvt = PV_cache[3 * order + 1];
	    }
	    else {
		mat = NULL;
	    }
	    bzero(dofvalues, m * n * sizeof(*dofvalues));
	    quad = phgQuadGetQuad2D(2 * order);
	    qp = quad->points;
	    qw = quad->weights;
	    /* compute	no0 = first basis on the face,
	     *		no1 = last basis + 1 on the face */
	    if (!DofIsHP(dof)) {
		np2 = t->np_edge;
		no0 = NVert + NEdge * np2 + index * n;
	    }
	    else {
		pe[0] = 0;
		for (ii = 0; ii < NEdge; ii++)
		    pe[ii + 1] = pe[ii] + dof->hp->info->types
			[dof->hp->edge_order[e->edges[ii]]]->np_edge;
		no0 = NVert + pe[NEdge];
		for (ii = 0; ii < index; ii++)
		    no0 += dof->hp->info->types
				[dof->hp->face_order[e->faces[ii]]]->np_face;
		np2 = n0 = 0;	/* avoid gcc warnings */
	    }
	    no1 = no0 + n;
#define A(i,j) mat[(i) * n + (j)]
#define B(i,j) dofvalues[(i) * m + (j)]
	    for (k = 0; k < quad->npoints; k++) {
		lambda[v0] = *(qp++);
		lambda[v1] = *(qp++);
		lambda[v2] = *(qp++);
		w = *(qw++);
		if (userfunc_lambda != NULL) {
		    userfunc_lambda(dof0, e, 0, lambda, f);
		}
		else {
		    phgGeomLambda2XYZ(dof->g, e, lambda, &x, &y, &z);
		    userfunc(x, y, z, f);
		}
		/* Note: must call BasFuncs after user_func since the latter
		 * might call BasFuncs, such as with phgDofInterpC2FGeneric */
		bas = t->BasFuncs(dof, e, NVert, no1, lambda);
		p = bas + no0 - NVert;
		for (i = 0; i < n; i++) {
		    b = p[i] * w;
		    if (mat != NULL) {
			for (j = 0; j <= i; j++) {
			    A(i,j) += (a = b * p[j]);
			    if (j < i)
				A(j,i) += a;
			}
		    }
		    for (j = 0; j < m; j++) {
			a = f[j] - (lambda[v0] * f0[j] + lambda[v1] * f1[j] +
				    lambda[v2] * f2[j]);
			if (!DofIsHP(dof)) {
			    for (ii = 0; ii < np2; ii++)
				a -= bas[i0 * np2 + ii] * g0[ii * m + j] +
				     bas[i1 * np2 + ii] * g1[ii * m + j] +
				     bas[i2 * np2 + ii] * g2[ii * m + j];
			}
			else {
			    n0 = pe[i0];
			    np2 = pe[i0 + 1] - n0;
			    for (ii = 0; ii < np2; ii++)
				a -= bas[n0 + ii] * g0[ii * m + j];
			    n0 = pe[i1];
			    np2 = pe[i1 + 1] - n0;
			    for (ii = 0; ii < np2; ii++)
				a -= bas[n0 + ii] * g1[ii * m + j];
			    n0 = pe[i2];
			    np2 = pe[i2 + 1] - n0;
			    for (ii = 0; ii < np2; ii++)
				a -= bas[n0 + ii] * g2[ii * m + j];
			}
			B(i,j) += b * a;
		    }
		}
	    }
#undef A
#undef B
	    /* solve the equation A x = B */
	    if (mat != NULL) {
    		show_condition(n, mat, type, index);
		if (!phgSolverDenseLU(n, (void *)mat, pvt))
		    phgError(1, "(%s:%d) singular system!\n",__FILE__,__LINE__);
	    }
	    else {
		mat = LU_cache[3 * order + 1];
		pvt = PV_cache[3 * order + 1];
	    }
	    phgSolverDenseSV(n, (void *)mat, pvt, m, (void *)dofvalues);
	    break;
	case VOLUME:
	    /* compute DOF by L2-projection on the element */
	    t = dof->type;
	    if (t == NULL)
		t = dof->hp->info->types[dof->hp->elem_order[e->index]];
	    if ((n = t->np_elem) == 0)
		break;
	    order = t->order;
	    assert(order >= 1 && order <= ORDER_MAX);
	    assert(index == 0);
	    GetElementVertices(e, 0, v0, v1, v2, v3,	v4, v5, v6, v7);
	    m = dof->dim;
	    f = phgAlloc(m * sizeof(*f));
	    /* the next 3 lines make gcc happy */
	    f0 = f1 = f2 = f3 = g0 = g1 = g2 = g3 = g4 = g5 = h0 = h1 = h2 = h3
		= NULL;
	    i0 = i1 = i2 = i3 = i4 = i5 = 0;
	    if (pdofvalues == NULL) {
		f0 = DofVertexData(dof, e->verts[v0]);
		f1 = DofVertexData(dof, e->verts[v1]);
		f2 = DofVertexData(dof, e->verts[v2]);
		f3 = DofVertexData(dof, e->verts[v3]);
		g0 = DofEdgeData(dof, e->edges[i0 = GetEdgeNo(v0, v1)]);
		g1 = DofEdgeData(dof, e->edges[i1 = GetEdgeNo(v0, v2)]);
		g2 = DofEdgeData(dof, e->edges[i2 = GetEdgeNo(v0, v3)]);
		g3 = DofEdgeData(dof, e->edges[i3 = GetEdgeNo(v1, v2)]);
		g4 = DofEdgeData(dof, e->edges[i4 = GetEdgeNo(v1, v3)]);
		g5 = DofEdgeData(dof, e->edges[i5 = GetEdgeNo(v2, v3)]);
		h0 = DofFaceData(dof, e->faces[v0]);
		h1 = DofFaceData(dof, e->faces[v1]);
		h2 = DofFaceData(dof, e->faces[v2]);
		h3 = DofFaceData(dof, e->faces[v3]);
	    }
	    else {
		f0 = pdofvalues[v0];
		f1 = pdofvalues[v1];
		f2 = pdofvalues[v2];
		f3 = pdofvalues[v3];
		g0 = pdofvalues[NVert + (i0 = GetEdgeNo(v0, v1))];
		g1 = pdofvalues[NVert + (i1 = GetEdgeNo(v0, v2))];
		g2 = pdofvalues[NVert + (i2 = GetEdgeNo(v0, v3))];
		g3 = pdofvalues[NVert + (i3 = GetEdgeNo(v1, v2))];
		g4 = pdofvalues[NVert + (i4 = GetEdgeNo(v1, v3))];
		g5 = pdofvalues[NVert + (i5 = GetEdgeNo(v2, v3))];
		h0 = pdofvalues[NVert + NEdge + v0];
		h1 = pdofvalues[NVert + NEdge + v1];
		h2 = pdofvalues[NVert + NEdge + v2];
		h3 = pdofvalues[NVert + NEdge + v3];
	    }
	    /* compute the DOF by solving dof->dim systems of linear equations
	     * of size n, dof->data or pdofvalues[] should point to valid
	     * DOF values on the vertices, edges and faces */
	    if (LU_cache[3 * order + 2] == NULL) {
		LU_cache[3 * order + 2] = phgCalloc(n * n, sizeof(*mat));
		FreeAtExitNoCheck(LU_cache[3 * order + 2]);
		PV_cache[3 * order + 2] = phgAlloc(n * sizeof(*pvt));
		FreeAtExitNoCheck(PV_cache[3 * order + 2]);
		mat = LU_cache[3 * order + 2];
		pvt = PV_cache[3 * order + 2];
	    }
	    else {
		mat = NULL;
	    }
	    bzero(dofvalues, m * n * sizeof(*dofvalues));
	    quad = phgQuadGetQuad3D(2 * order);
	    qp = quad->points;
	    qw = quad->weights;
	    /* compute	no0 = first basis on the element,
	     *		no1 = last basis + 1 on the element */
	    if (!DofIsHP(dof)) {
		np2 = t->np_edge;
		np3 = t->np_face;
		no0 = NVert + (n0 = NEdge * np2) + NFace * np3;
	    }
	    else {
		pe[0] = 0;
		for (ii = 0; ii < NEdge; ii++)
		    pe[ii + 1] = pe[ii] + dof->hp->info->types
			[dof->hp->edge_order[e->edges[ii]]]->np_edge;
		pf[0] = 0;
		for (ii = 0; ii < NFace; ii++)
		    pf[ii + 1] = pf[ii] + dof->hp->info->types
			[dof->hp->face_order[e->faces[ii]]]->np_face;
		no0 = NVert + pe[NEdge] + pf[NFace];
		np2 = np3 = n0 = 0;		/* avoid gcc warnings */
	    }
	    no1 = no0 + n;
#define A(i,j) mat[(i) * n + (j)]
#define B(i,j) dofvalues[(i) * m + (j)]
	    for (k = 0; k < quad->npoints; k++) {
		lambda[v0] = *(qp++);
		lambda[v1] = *(qp++);
		lambda[v2] = *(qp++);
		lambda[v3] = *(qp++);
		w = *(qw++);
		if (userfunc_lambda != NULL) {
		    userfunc_lambda(dof0, e, 0, lambda, f);
		}
		else {
		    phgGeomLambda2XYZ(dof->g, e, lambda, &x, &y, &z);
		    userfunc(x, y, z, f);
		}
		/* Note: must call BasFuncs after user_func since the latter
		 * might call BasFuncs, such as with phgDofInterpC2FGeneric */
		bas = t->BasFuncs(dof, e, NVert, no1, lambda);
		p = bas + no0 -NVert;
		for (i = 0; i < n; i++) {
		    b = p[i] * w;
		    if (mat != NULL) {
			for (j = 0; j <= i; j++) {
			    A(i,j) += (a = b * p[j]);
			    if (j < i)
				A(j,i) += a;
			}
		    }
		    for (j = 0; j < m; j++) {
			a = f[j] - (lambda[v0] * f0[j] + lambda[v1] * f1[j] +
				    lambda[v2] * f2[j] + lambda[v3] * f3[j]);
			if (!DofIsHP(dof)) {
			    for (ii = 0; ii < np2; ii++)
				a -= bas[i0 * np2 + ii] * g0[ii * m + j] +
				     bas[i1 * np2 + ii] * g1[ii * m + j] +
				     bas[i2 * np2 + ii] * g2[ii * m + j] +
				     bas[i3 * np2 + ii] * g3[ii * m + j] +
				     bas[i4 * np2 + ii] * g4[ii * m + j] +
				     bas[i5 * np2 + ii] * g5[ii * m + j];
			    for (ii = 0; ii < np3; ii++)
				a -= bas[n0 + v0 * np3 + ii] * h0[ii * m + j] +
				     bas[n0 + v1 * np3 + ii] * h1[ii * m + j] +
				     bas[n0 + v2 * np3 + ii] * h2[ii * m + j] +
				     bas[n0 + v3 * np3 + ii] * h3[ii * m + j];
		        }
			else {
			    n0 = pe[i0];
			    np2 = pe[i0 + 1] - n0;
			    for (ii = 0; ii < np2; ii++)
				a -= bas[n0 + ii] * g0[ii * m + j];

			    n0 = pe[i1];
			    np2 = pe[i1 + 1] - n0;
			    for (ii = 0; ii < np2; ii++)
				a -= bas[n0 + ii] * g1[ii * m + j];

			    n0 = pe[i2];
			    np2 = pe[i2 + 1] - n0;
			    for (ii = 0; ii < np2; ii++)
				a -= bas[n0 + ii] * g2[ii * m + j];

			    n0 = pe[i3];
			    np2 = pe[i3 + 1] - n0;
			    for (ii = 0; ii < np2; ii++)
				a -= bas[n0 + ii] * g3[ii * m + j];

			    n0 = pe[i4];
			    np2 = pe[i4 + 1] - n0;
			    for (ii = 0; ii < np2; ii++)
				a -= bas[n0 + ii] * g4[ii * m + j];

			    n0 = pe[i5];
			    np2 = pe[i5 + 1] - n0;
			    for (ii = 0; ii < np2; ii++)
				a -= bas[n0 + ii] * g5[ii * m + j];

			    n0 = pf[v0];
			    np3 = pf[v0 + 1] - n0;
			    n0 += pe[NEdge];
			    for (ii = 0; ii < np3; ii++)
				a -= bas[n0 + ii] * h0[ii * m + j];

			    n0 = pf[v1];
			    np3 = pf[v1 + 1] - n0;
			    n0 += pe[NEdge];
			    for (ii = 0; ii < np3; ii++)
				a -= bas[n0 + ii] * h1[ii * m + j];

			    n0 = pf[v2];
			    np3 = pf[v2 + 1] - n0;
			    n0 += pe[NEdge];
			    for (ii = 0; ii < np3; ii++)
				a -= bas[n0 + ii] * h2[ii * m + j];

			    n0 = pf[v3];
			    np3 = pf[v3 + 1] - n0;
			    n0 += pe[NEdge];
			    for (ii = 0; ii < np3; ii++)
				a -= bas[n0 + ii] * h3[ii * m + j];
			}
			B(i,j) += b * a;
		    }
		}
	    }
#undef A
#undef B
	    /* solve the equation A x = B */
	    if (mat != NULL) {
    		show_condition(n, mat, type, index);
		if (!phgSolverDenseLU(n, (void *)mat, pvt))
		    phgError(1, "(%s:%d) singular system!\n",__FILE__,__LINE__);
	    }
	    else {
		mat = LU_cache[3 * order + 2];
		pvt = PV_cache[3 * order + 2];
	    }
	    phgSolverDenseSV(n, (void *)mat, pvt, m, (void *)dofvalues);
	    break;
    }

    phgFree(f);
    return;
}

static void
HB_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
{
    DOF_TYPE *t = dof->type;
    FLOAT *pv[15];
    int i, m = dof->dim;
    DOF *dof1;

    dof0 = dof;
    if (t == NULL)
	t = dof->hp->info->types[dof->hp->elem_order[e->index]];
    if (t->np_vert == 1) {
	/* not DG */
	HB_init0(dof, e, type, index, userfunc, userfunc_lambda,
		 funcvalues, dofvalues, pdofvalues);
	return;
    }

    /* map DG basis functions to HB basis functions */
    assert(type == VOLUME);
    dof1 = phgAlloc(sizeof(*dof));
    memcpy(dof1, dof, sizeof(*dof));
    dof1->type = t->base_type;
    dof1->hp = NULL;
    assert(!DofIsHP(dof1));

    for (i = 0; i < NVert; i++) {
	HB_init0(dof1, e, VERTEX, i, userfunc, userfunc_lambda,
		 funcvalues, dofvalues, NULL);
	pv[i] = dofvalues;
	dofvalues += m;
    }

    for (i = 0; i < NEdge; i++) {
	HB_init0(dof1, e, EDGE, i, userfunc, userfunc_lambda,
		 funcvalues, dofvalues, pv);
	pv[NVert + i] = dofvalues;
	dofvalues += m * dof1->type->np_edge;
    }

    for (i = 0; i < NFace; i++) {
	HB_init0(dof1, e, FACE, i, userfunc, userfunc_lambda,
		 funcvalues, dofvalues, pv);
	pv[NVert + NEdge + i] = dofvalues;
	dofvalues += m * dof1->type->np_face;
    }

    HB_init0(dof1, e, VOLUME, 0, userfunc, userfunc_lambda,
	     funcvalues, dofvalues, pv);

    phgFree(dof1);
}

/*-------------------------- DOF types HBn --------------------------*/

#define NBasEdge(p)	(p > 1 ? (p - 1) : 0)
#define NBasFace(p)	(p > 2 ? (p - 1) * (p - 2) / 2 : 0)
#define NBasElem(p)	(p > 3 ? (p - 1) * (p - 2) * (p - 3) / 6 : 0)

#define	D1(a)	a
#define	D2(a)	a, a
#define	D3(a)	a, D2(a)
#define	D4(a)	D2(a), D2(a)
#define	D5(a)	D2(a), D3(a)
#define	D6(a)	D3(a), D3(a)
#define	D7(a)	D3(a), D4(a)
#define	D8(a)	D4(a), D4(a)
#define	D9(a)	D4(a), D5(a)
#define	D10(a)	D5(a), D5(a)
#define	D11(a)	D5(a), D6(a)
#define	D12(a)	D6(a), D6(a)
#define	D13(a)	D6(a), D7(a)
#define	D14(a)	D7(a), D7(a)
#define	D15(a)	D7(a), D8(a)
#define	D16(a)	D8(a), D8(a)

#define D21(a)  D15(a), D6(a)
#define D28(a)  D21(a), D7(a)
#define D36(a)  D28(a), D8(a)
#define D45(a)  D36(a), D9(a)
#define D55(a)  D45(a), D10(a)
#define D66(a)  D55(a), D11(a)
#define D78(a)  D66(a), D12(a)
#define D91(a)  D78(a), D13(a)

#define NBasEdge(p)	(p > 1 ? (p - 1) : 0)
#define NBasFace(p)	(p > 2 ? (p - 1) * (p - 2) / 2 : 0)
#define NBasElem(p)	(p > 3 ? (p - 1) * (p - 2) * (p - 3) / 6 : 0)

HP_INFO hp_info_HB = {
    DOF_HBn,		/* types */
    "Hierarchic H1",	/* name */
    FE_H1,		/* FE space */
    1,			/* dim */
    1,			/* min_order */
    15,			/* max_order */
    VERT_FLAG | EDGE_FLAG | FACE_FLAG | ELEM_FLAG
};

#define DofType(order, type, name, orders, grad, invariant)		\
    DOF_TYPE type = {DofReserved, name, NULL, orders,			\
			grad, NULL, &hp_info_HB,			\
			phgDofInterpC2FGeneric,				\
			phgDofInterpF2CGeneric,				\
			HB_init, HB_bas, HB_grad, NULL, FE_H1,		\
			FALSE, invariant, FALSE, -1,			\
			NVert + NEdge * NBasEdge(order)	+		\
			    NFace * NBasFace(order) + NBasElem(order),	\
			order,		/* order */			\
			0,		/* D_order */			\
			0,		/* continuity */		\
			1,		/* dim */			\
			1, NBasEdge(order), NBasFace(order), NBasElem(order)};

static BYTE HB1_orders[] = {1};
DofType(1, DOF_HB1_, "HB1", HB1_orders, DOF_DG0, TRUE)

static BYTE HB2_orders[] = {
    1,		/* vertex */
    2		/* edge */
};
DofType(2, DOF_HB2_, "HB2", HB2_orders, DOF_DG1, TRUE)

static BYTE HB3_orders[] = {/* vertex: */ 1, /* edge: */ 2, 3, /* face: */ 3};
DofType(3, DOF_HB3_, "HB3", HB3_orders, DOF_DG2, FALSE)

static BYTE HB4_orders[] = {
    1,		/* vertex */
    2, 3, 4,	/* edge */
    3, D2(4),	/* face */
    4		/* element */
};
DofType(4, DOF_HB4_, "HB4", HB4_orders, DOF_DG3, FALSE)

static BYTE HB5_orders[] = {
    1,			/* vertex */
    2, 3, 4, 5,		/* edge */
    3, D2(4), D3(5),	/* face */
    4, D3(5)		/* element (Note: Dn(k), n = (k-2)*(k-3)/2) */
};
DofType(5, DOF_HB5_, "HB5", HB5_orders, DOF_DG4, FALSE)

static BYTE HB6_orders[] = {
    1,				/* vertex */
    2, 3, 4, 5, 6,		/* edge */
    3, D2(4), D3(5), D4(6),	/* face */
    4, D3(5), D6(6) 		/* element */
};
DofType(6, DOF_HB6_, "HB6", HB6_orders, DOF_DG5, FALSE)

static BYTE HB7_orders[] = {
    1,					/* vertex */
    2, 3, 4, 5, 6, 7,			/* edge */
    3, D2(4), D3(5), D4(6), D5(7),	/* face */
    4, D3(5), D6(6), D10(7)		/* element */
};
DofType(7, DOF_HB7_, "HB7", HB7_orders, DOF_DG6, FALSE)

static BYTE HB8_orders[] = {
    1,						/* vertex */
    2, 3, 4, 5, 6, 7, 8,			/* edge */
    3, D2(4), D3(5), D4(6), D5(7), D6(8),	/* face */
    4, D3(5), D6(6), D10(7), D15(8)		/* element */
};
DofType(8, DOF_HB8_, "HB8", HB8_orders, DOF_DG7, FALSE)

static BYTE HB9_orders[] = {
    1,							/* vertex */
    2, 3, 4, 5, 6, 7, 8, 9,				/* edge */
    3, D2(4), D3(5), D4(6), D5(7), D6(8), D7(9),	/* face */
    4, D3(5), D6(6), D10(7), D15(8), D21(9)		/* element */
};
DofType(9, DOF_HB9_, "HB9", HB9_orders, DOF_DG8, FALSE)

static BYTE HB10_orders[] = {
    1,								/* vertex */
    2, 3, 4, 5, 6, 7, 8, 9, 10,					/* edge */
    3, D2(4), D3(5), D4(6), D5(7), D6(8), D7(9), D8(10),	/* face */
    4, D3(5), D6(6), D10(7), D15(8), D21(9), D28(10)		/* element */
};
DofType(10, DOF_HB10_, "HB10", HB10_orders, DOF_DG9, FALSE)

static BYTE HB11_orders[] = {
    1,								/* vertex */
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11,				/* edge */
    3, D2(4), D3(5), D4(6), D5(7), D6(8), D7(9), D8(10), D9(11),/* face */
    4, D3(5), D6(6), D10(7), D15(8), D21(9), D28(10), D36(11)	/* element */
};
DofType(11, DOF_HB11_, "HB11", HB11_orders, DOF_DG10, FALSE)

static BYTE HB12_orders[] = {
    1,								/* vertex */
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,				/* edge */
    3, D2(4), D3(5), D4(6), D5(7), D6(8), D7(9), D8(10), D9(11), D10(12),
								/* face */
    4, D3(5), D6(6), D10(7), D15(8), D21(9), D28(10), D36(11), D45(12)
								/* element */
};
DofType(12, DOF_HB12_, "HB12", HB12_orders, DOF_DG11, FALSE)

static BYTE HB13_orders[] = {
    1,								/* vertex */
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,			/* edge */
    3, D2(4), D3(5), D4(6), D5(7), D6(8), D7(9), D8(10), D9(11), D10(12),
	D11(13),						/* face */
    4, D3(5), D6(6), D10(7), D15(8), D21(9), D28(10), D36(11), D45(12), D55(13)
								/* element */
};
DofType(13, DOF_HB13_, "HB13", HB13_orders, DOF_DG12, FALSE)

static BYTE HB14_orders[] = {
    1,								/* vertex */
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,			/* edge */
    3, D2(4), D3(5), D4(6), D5(7), D6(8), D7(9), D8(10), D9(11), D10(12),
	D11(13), D12(14),					/* face */
    4, D3(5), D6(6), D10(7), D15(8), D21(9), D28(10), D36(11), D45(12),
	D55(13), D66(14)					/* element */
};
DofType(14, DOF_HB14_, "HB14", HB14_orders, DOF_DG13, FALSE)

static BYTE HB15_orders[] = {
    1,								/* vertex */
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,		/* edge */
    3, D2(4), D3(5), D4(6), D5(7), D6(8), D7(9), D8(10), D9(11), D10(12),
	D11(13), D12(14), D13(15),				/* face */
    4, D3(5), D6(6), D10(7), D15(8), D21(9), D28(10), D36(11), D45(12),
	D55(13), D66(14), D78(15)				/* element */
};
DofType(15, DOF_HB15_, "HB15", HB15_orders, DOF_DG14, FALSE)

static BYTE HB16_orders[] = {
    1,								/* vertex */
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,		/* edge */
    3, D2(4), D3(5), D4(6), D5(7), D6(8), D7(9), D8(10), D9(11), D10(12),
	D11(13), D12(14), D13(15), D14(16),			/* face */
    4, D3(5), D6(6), D10(7), D15(8), D21(9), D28(10), D36(11), D45(12),
	D55(13), D66(14), D78(15), D91(16)			/* element */
};
DofType(16, DOF_HB16_, "HB16", HB16_orders, DOF_DG15, FALSE)

DOF_TYPE *DOF_HBn[] = {DOF_DG0,
    DOF_HB1,  DOF_HB2,  DOF_HB3,  DOF_HB4,  DOF_HB5,
    DOF_HB6,  DOF_HB7,  DOF_HB8,  DOF_HB9,  DOF_HB10,
    DOF_HB11, DOF_HB12, DOF_HB13, DOF_HB14, DOF_HB15,
    DOF_HB16, NULL
};

# if 0
/* Note: to test the length of the orders[] array in DOF_HBx:
int i;
for (i = 1; i <= 15; i++) {
    extern int hb_dummy[];
    DOF_TYPE *t = DOF_HBn[i];
	    printf("Testing %s->orders[]: %s\n", t->name,
	hb_dummy[i] == t->np_vert+t->np_edge+t->np_face+t->np_elem ?
	                                "passed" : "failed");
} */
int hb_dummy[] = {0,
    sizeof(HB1_orders)/sizeof(HB1_orders[0]),
    sizeof(HB2_orders)/sizeof(HB2_orders[0]),
    sizeof(HB3_orders)/sizeof(HB3_orders[0]),
    sizeof(HB4_orders)/sizeof(HB4_orders[0]),
    sizeof(HB5_orders)/sizeof(HB5_orders[0]),
    sizeof(HB6_orders)/sizeof(HB6_orders[0]),
    sizeof(HB7_orders)/sizeof(HB7_orders[0]),
    sizeof(HB8_orders)/sizeof(HB8_orders[0]),
    sizeof(HB9_orders)/sizeof(HB9_orders[0]),
    sizeof(HB10_orders)/sizeof(HB10_orders[0]),
    sizeof(HB11_orders)/sizeof(HB11_orders[0]),
    sizeof(HB12_orders)/sizeof(HB12_orders[0]),
    sizeof(HB13_orders)/sizeof(HB13_orders[0]),
    sizeof(HB14_orders)/sizeof(HB14_orders[0]),
    sizeof(HB15_orders)/sizeof(HB15_orders[0]),
    sizeof(HB16_orders)/sizeof(HB16_orders[0])};
#endif

/*------------------------- Hcurl conforming bases -------------------------*/

static const FLOAT *
HC_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT *buffers[ORDER_MAX + 1];
    static FLOAT phi[2 * Dim];
    static BOOLEAN initialized = FALSE;
#if USE_OMP
#pragma omp threadprivate(buffers, phi, initialized)
#endif	/* USE_OMP */

    FLOAT *L;
    const FLOAT (*nabla)[Dim + 1] = (void *)(phgGeomGetJacobian(dof->g, e));

    DOF_TYPE *t;
    int order0, order, nbas, bufferlen;
    FLOAT *phi0, *phi1, *bas, *p, *p0, a1, b1, c1, a2, b2, c2, d;
    const FLOAT *bas_hb;
    DOF hb;
    COORD *coord;
    int i, j, i0, v0, v1, v2, v3, v4, v5, v6, v7, m, n;

    if (!initialized) {
	initialized = TRUE;
	bzero(buffers, sizeof(buffers));
    }

    t = dof->type;
    if (t != NULL) {
	order = t->order;
	bufferlen = nbas = t->nbas;
	if (order == 1 && t == DOF_HC0)
	    order = 0;
	order0 = order;
    }
    else {
	assert(DofIsHP(dof));
	/* FIXME: it is assumed that vert/edge/face orders <= elem order */
	order = order0 = dof->hp->elem_order[e->index];
	nbas = DofNBas(dof, e);
	bufferlen = dof->hp->info->types[order0]->nbas;
    }

    assert(order >= 0 && order <= ORDER_MAX);
    if (no1 <= 0)
	no1 = nbas;
 
    if (no1 <= no0)
	return NULL;

    if ((bas = buffers[order0]) == NULL) {
	buffers[order0] = phgAlloc(Dim * bufferlen * sizeof(*bas));
	FreeAtExitNoCheck(buffers[order0]);
	bas = buffers[order0];
    }

    p = bas;
    p0 = p + Dim * (no1 - no0);

#if CACHE_LP
    if (L0 == NULL) {
	L0 = phgAlloc(NEdge * sizeof(*L0));
	DL0 = phgAlloc(NEdge * sizeof(*DL0));
	cache_order = phgAlloc(NEdge * sizeof(*cache_order));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
	FreeAtExitNoCheck(cache_order);
    }
    cache_reset = FALSE;
    memset(cache_order, 0, NEdge * sizeof(*cache_order));
#else	/* CACHE_LP */
    if (L0 == NULL) {
	L0 = phgAlloc(3 * sizeof(*L0));
	DL0 = phgAlloc(3 * sizeof(*DL0));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
    }
#endif	/* CACHE_LP */

    /*------------------------- edge functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NEdge && p < p0; i++) {
	if (DofIsHP(dof)) {
	    order = dof->hp->edge_order[e->edges[i]];
	    assert(order <= order0);
	    t = dof->hp->info->types[order];
	}
	if (t->np_edge == 0)
	    continue;
	n = i0 + t->np_edge;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	GetEdgeVertices(e, i, v0, v1);
	/* first edge basis */
	if (i0++ >= no0) {
	    phi0 = p;
	    p += 3;
	}
	else {
	    phi0 = phi;
	}
	phi0[0] = lambda[v1] * nabla[v0][0] - lambda[v0] * nabla[v1][0];
	phi0[1] = lambda[v1] * nabla[v0][1] - lambda[v0] * nabla[v1][1];
	phi0[2] = lambda[v1] * nabla[v0][2] - lambda[v0] * nabla[v1][2];
	if (p >= p0)
	    break;
	if (t->np_edge == 1)
	    continue;
	/* second edge basis */
	if (i0++ >= no0) {
	    phi1 = p;
	    p += 3;
	}
	else {
	    phi1 = phi + 3;
	}
	phi1[0] = lambda[v1] * nabla[v0][0] + lambda[v0] * nabla[v1][0];
	phi1[1] = lambda[v1] * nabla[v0][1] + lambda[v0] * nabla[v1][1];
	phi1[2] = lambda[v1] * nabla[v0][2] + lambda[v0] * nabla[v1][2];
	if (p >= p0)
	    break;
	if (t->np_edge == 2)
	    continue;

#if CACHE_LP
	get_L(order0 - 1, v0, v1, lambda, &L, NULL);
#else	/* CACHE_LP */
	phgLegendreP(order - 1, lambda[v1] - lambda[v0], L = L0[0], NULL);
#endif	/* CACHE_LP */
	for (j = 2; i0 < n && p < p0; i0++, j++) {
	    if (i0 < no0)
		continue;
	    /* note: l = j - 1 */
	    b1 = 1. / (FLOAT)j;
	    a1 = (2 * j - 1) * b1;
	    b1 *= j - 1;
	    *(p++) = a1 * L[j - 1] * phi1[0] - b1 * L[j - 2] * phi0[0];
	    *(p++) = a1 * L[j - 1] * phi1[1] - b1 * L[j - 2] * phi0[1];
	    *(p++) = a1 * L[j - 1] * phi1[2] - b1 * L[j - 2] * phi0[2];
	}
    }

    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return bas;
    }

    if ((no0 -= i0) < 0)
	no0 = 0;

    hb.hp = NULL;

    /*------------------------- face functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NFace && p < p0; i++) {
	if (DofIsHP(dof)) {
	    order = dof->hp->face_order[e->faces[i]];
	    assert(order <= order0);
	    t = dof->hp->info->types[order];
	}
	if (t->np_face == 0)
	    continue;
	n = i0 + t->np_face;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	hb.type = DOF_HBn[order];
	GetFaceVertices(e, i, v0, v1, v2, v3);
	/* edge-based */
#define EDGE_BASED(v0, v1, v2)						\
	m = NVert + GetEdgeNo(v0, v1) * hb.type->np_edge;		\
	bas_hb = HB_bas(&hb, e, m, m + hb.type->np_edge, lambda);	\
	for (j = 0; j < hb.type->np_edge && p < p0; j++, bas_hb++, i0++) { \
	    if (i0 < no0)						\
		continue;						\
            *(p++) = *bas_hb * nabla[v2][0];				\
            *(p++) = *bas_hb * nabla[v2][1];				\
            *(p++) = *bas_hb * nabla[v2][2];				\
	}								\
	if (p >= p0)							\
	    break;
	EDGE_BASED(v0, v1, v2);
	EDGE_BASED(v0, v2, v1);
	EDGE_BASED(v1, v2, v0);
#undef EDGE_BASED
	if (t->np_edge == 3 * hb.type->np_edge)
	    continue;
	/* face-based */
	/* get tangent vectors on the face */
	coord = dof->g->verts + e->verts[v1];
	a1 = (*coord)[0];
	b1 = (*coord)[1];
	c1 = (*coord)[2];
	coord = dof->g->verts + e->verts[v2];
	a2 = (*coord)[0];
	b2 = (*coord)[1];
	c2 = (*coord)[2];
	coord = dof->g->verts + e->verts[v0];
	a1 -= (*coord)[0];
	b1 -= (*coord)[1];
	c1 -= (*coord)[2];
	a2 -= (*coord)[0];
	b2 -= (*coord)[1];
	c2 -= (*coord)[2];
	/* normalize tangent vectors (not necessary?) */
	d = 1. / Sqrt(a1 * a1 + b1 * b1 + c1 * c1);
	a1 *= d;
	b1 *= d;
	c1 *= d;
	d = 1. / Sqrt(a2 * a2 + b2 * b2 + c2 * c2);
	a2 *= d;
	b2 *= d;
	c2 *= d;
	m = NVert + NEdge * hb.type->np_edge + i * hb.type->np_face;
	bas_hb = HB_bas(&hb, e, m, m + hb.type->np_face, lambda);
	for (j = 0; j < hb.type->np_face; j++, bas_hb++) {
	    if (i0++ >= no0) {
		*(p++) = *bas_hb * a1;
		*(p++) = *bas_hb * b1;
		*(p++) = *bas_hb * c1;
		if (p >= p0)
		    break;
	    }

	    if (i0++ >= no0) {
		*(p++) = *bas_hb * a2;
		*(p++) = *bas_hb * b2;
		*(p++) = *bas_hb * c2;
		if (p >= p0)
		    break;
	    }
	}
    }

    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return bas;
    }

    if ((no0 -= i0) < 0)
	no0 = 0;

    /*------------------------- element functions --------------------------*/
    if (DofIsHP(dof)) {
	order = order0;
	t = dof->hp->info->types[order];
    }
    n = t->np_elem;
    assert(no0 < n);
    assert(n > 0);
    hb.type = DOF_HBn[order];
    GetElementVertices(e, 0, v0, v1, v2, v3,
		       v4, v5, v6, v7);
    i0 = 0;
    /* face-based */
#define FACE_BASED(v)							\
    m = NVert + NEdge * hb.type->np_edge + v * hb.type->np_face;	\
    bas_hb = HB_bas(&hb, e, m, m + hb.type->np_face, lambda);		\
    for (j = 0; j < hb.type->np_face && p < p0; i0++, j++, bas_hb++) {	\
	if (i0 < no0)							\
	    continue;							\
	*(p++) = *bas_hb * nabla[v][0];					\
	*(p++) = *bas_hb * nabla[v][1];					\
	*(p++) = *bas_hb * nabla[v][2];					\
    }									\
    if (p >= p0) 							\
	goto cont; 
    FACE_BASED(v0);
    FACE_BASED(v1);
    FACE_BASED(v2);
    FACE_BASED(v3);
#undef FACE_BASED
cont:
    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return bas;
    }
    /* element-based */
    m = NVert + NEdge * hb.type->np_edge + NFace * hb.type->np_face;
    bas_hb = HB_bas(&hb, e, m, m + hb.type->np_elem, lambda);
    for (j = 0; j < hb.type->np_elem; j++, bas_hb++) {
	if (i0++ >= no0) {
	    *(p++) = *bas_hb;
	    *(p++) = 0.;
	    *(p++) = 0.;
	    if (p >= p0)
		break;
	}

	if (i0++ >= no0) {
	    *(p++) = 0.;
	    *(p++) = *bas_hb;
	    *(p++) = 0.;
	    if (p >= p0)
		break;
	}

	if (i0++ >= no0) {
	    *(p++) = 0.;
	    *(p++) = 0.;
	    *(p++) = *bas_hb;
	    if (p >= p0)
		break;
	}
    }

    assert(p == p0);

#if CACHE_LP
    cache_reset = TRUE;
#endif	/* CACHE_LP */
    return bas;
}

static const FLOAT *
HC_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT *buffers[ORDER_MAX + 1];
    static FLOAT dphi[2 * Dim * (Dim + 1)];
    static BOOLEAN initialized = FALSE;
#if USE_OMP
#pragma omp threadprivate(buffers, dphi, initialized)
#endif	/* USE_OMP */

    FLOAT *L, *DL;
    const FLOAT (*nabla)[Dim + 1] = (void *)(phgGeomGetJacobian(dof->g, e));

    DOF_TYPE *t;
    int order, order0, nbas, bufferlen;
    FLOAT *grad, *p, *p0, a1, b1, c1, a2, b2, c2, d;
    FLOAT (*dphi0)[Dim + 1], (*dphi1)[Dim + 1];
    const FLOAT *grad_hb;
    FLOAT phi0[3], phi1[3];
    DOF hb;
    COORD *coord;
    int i, j, i0, v0, v1, v2, v3, v4, v5, v6, v7, m, n;

    if (!initialized) {
	initialized = TRUE;
	bzero(buffers, sizeof(buffers));
    }

    t = dof->type;
    if (t != NULL) {
	order = t->order;
	bufferlen = nbas = t->nbas;
	if (order == 1 && t == DOF_HC0)
	    order = 0;
	order0 = order;
    }
    else {
	assert(DofIsHP(dof));
	/* FIXME: it is assumed that vert/edge/face orders <= elem order */
	order = order0 = dof->hp->elem_order[e->index];
	nbas = DofNBas(dof, e);
	bufferlen = dof->hp->info->types[order0]->nbas;
    }

    assert(order >= 0 && order <= ORDER_MAX);
    if (no1 <= 0)
	no1 = nbas;
 
    if (no1 <= no0)
	return NULL;

    if ((grad = buffers[order]) == NULL) {
	buffers[order] = phgAlloc((Dim + 1) * Dim * bufferlen * sizeof(*grad));
	FreeAtExitNoCheck(buffers[order]);
	grad = buffers[order];
    }

    p = grad;
    p0 = p + (Dim + 1) * Dim * (no1 - no0);
    bzero(p, (no1 - no0) * (Dim + 1) * Dim * sizeof(*p));

#if CACHE_LP
    if (L0 == NULL) {
	L0 = phgAlloc(NEdge * sizeof(*L0));
	DL0 = phgAlloc(NEdge * sizeof(*DL0));
	cache_order = phgAlloc(NEdge * sizeof(*cache_order));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
	FreeAtExitNoCheck(cache_order);
    }
    cache_reset = FALSE;
    memset(cache_order, 0, NEdge * sizeof(*cache_order));
#else	/* CACHE_LP */
    if (L0 == NULL) {
	L0 = phgAlloc(3 * sizeof(*L0));
	DL0 = phgAlloc(3 * sizeof(*DL0));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
    }
#endif	/* CACHE_LP */

    /*------------------------- edge functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NEdge && p < p0; i++) {
	if (DofIsHP(dof)) {
	    order = dof->hp->edge_order[e->edges[i]];
	    assert(order <= order0);
	    t = dof->hp->info->types[order];
	}
	if (t->np_edge == 0)
	    continue;
	n = i0 + t->np_edge;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	GetEdgeVertices(e, i, v0, v1);
	/* first edge basis */
	if (i0++ >= no0) {
	    dphi0 = (void *)p;
	    p += Dim * (Dim + 1);
	}
	else {
	    dphi0 = (void *)dphi;
	    memset(dphi0, 0, Dim * (Dim + 1) * sizeof(FLOAT));
	}
	dphi0[0][v0] = -nabla[v1][0];
	dphi0[0][v1] =  nabla[v0][0];
	dphi0[1][v0] = -nabla[v1][1];
	dphi0[1][v1] =  nabla[v0][1];
	dphi0[2][v0] = -nabla[v1][2];
	dphi0[2][v1] =  nabla[v0][2];
	if (p >= p0)
	    break;
	if (t->np_edge == 1)
	    continue;
	/* second edge basis */
	if (i0++ >= no0) {
	    dphi1 = (void *)p;
	    p += Dim * (Dim + 1);
	}
	else {
	    dphi1 = (void *)(dphi + Dim * (Dim + 1));
	    memset(dphi1, 0, Dim * (Dim + 1) * sizeof(FLOAT));
	}
	dphi1[0][v0] = nabla[v1][0];
	dphi1[0][v1] = nabla[v0][0];
	dphi1[1][v0] = nabla[v1][1];
	dphi1[1][v1] = nabla[v0][1];
	dphi1[2][v0] = nabla[v1][2];
	dphi1[2][v1] = nabla[v0][2];
	if (p >= p0)
	    break;
	if (t->np_edge == 2)
	    continue;

	/* need values of the first edge basis functions */
        phi0[0] = lambda[v1] * nabla[v0][0] - lambda[v0] * nabla[v1][0];
        phi0[1] = lambda[v1] * nabla[v0][1] - lambda[v0] * nabla[v1][1];
        phi0[2] = lambda[v1] * nabla[v0][2] - lambda[v0] * nabla[v1][2];

        phi1[0] = lambda[v1] * nabla[v0][0] + lambda[v0] * nabla[v1][0];
        phi1[1] = lambda[v1] * nabla[v0][1] + lambda[v0] * nabla[v1][1];
        phi1[2] = lambda[v1] * nabla[v0][2] + lambda[v0] * nabla[v1][2];

#if CACHE_LP
	get_L(order0 - 1, v0, v1, lambda, &L, &DL);
#else	/* CACHE_LP */
	phgLegendreP(order - 1, lambda[v1] - lambda[v0],
					L = L0[0], DL = DL0[0]);
#endif	/* CACHE_LP */
	for (j = 2; i0 < n && p < p0; i0++, j++) {
	    if (i0 < no0)
		continue;
	   /* note: l = j - 1 */
	    b1 = 1. / (FLOAT)j;
	    a1 = (2 * j - 1) * b1;
	    b1 *= j - 1;
	    p[v0] = a1 * (-DL[j - 1] * phi1[0] + L[j - 1] * dphi1[0][v0]) -
		    b1 * (-DL[j - 2] * phi0[0] + L[j - 2] * dphi0[0][v0]);
	    p[v1] = a1 * ( DL[j - 1] * phi1[0] + L[j - 1] * dphi1[0][v1]) -
		    b1 * ( DL[j - 2] * phi0[0] + L[j - 2] * dphi0[0][v1]);
	    p += Dim + 1;
	    p[v0] = a1 * (-DL[j - 1] * phi1[1] + L[j - 1] * dphi1[1][v0]) -
		    b1 * (-DL[j - 2] * phi0[1] + L[j - 2] * dphi0[1][v0]);
	    p[v1] = a1 * ( DL[j - 1] * phi1[1] + L[j - 1] * dphi1[1][v1]) -
		    b1 * ( DL[j - 2] * phi0[1] + L[j - 2] * dphi0[1][v1]);
	    p += Dim + 1;
	    p[v0] = a1 * (-DL[j - 1] * phi1[2] + L[j - 1] * dphi1[2][v0]) -
		    b1 * (-DL[j - 2] * phi0[2] + L[j - 2] * dphi0[2][v0]);
	    p[v1] = a1 * ( DL[j - 1] * phi1[2] + L[j - 1] * dphi1[2][v1]) -
		    b1 * ( DL[j - 2] * phi0[2] + L[j - 2] * dphi0[2][v1]);
	    p += Dim + 1;
	}
    }

    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return grad;
    }

    if ((no0 -= i0) < 0)
	no0 = 0;

    hb.hp = NULL;

    /*------------------------- face functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NFace && p < p0; i++) {
	if (DofIsHP(dof)) {
	    order = dof->hp->face_order[e->faces[i]];
	    assert(order <= order0);
	    t = dof->hp->info->types[order];
	}
	if (t->np_face == 0)
	    continue;
	n = i0 + t->np_face;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	hb.type = DOF_HBn[order];
	GetFaceVertices(e, i, v0, v1, v2, v3);
	/* edge-based */
#define EDGE_BASED(v0, v1, v2)						\
	m = NVert + GetEdgeNo(v0, v1) * hb.type->np_edge;		\
	grad_hb = HB_grad(&hb, e, m, m + hb.type->np_edge, lambda);	\
	for (j = 0; j < hb.type->np_edge && p < p0;			\
		    i0++, j++, grad_hb += Dim + 1) {			\
	    if (i0 < no0)						\
		continue;						\
            *(p++) = grad_hb[0] * nabla[v2][0];				\
            *(p++) = grad_hb[1] * nabla[v2][0];				\
            *(p++) = grad_hb[2] * nabla[v2][0];				\
            *(p++) = grad_hb[3] * nabla[v2][0];				\
									\
            *(p++) = grad_hb[0] * nabla[v2][1];				\
            *(p++) = grad_hb[1] * nabla[v2][1];				\
            *(p++) = grad_hb[2] * nabla[v2][1];				\
            *(p++) = grad_hb[3] * nabla[v2][1];				\
									\
            *(p++) = grad_hb[0] * nabla[v2][2];				\
            *(p++) = grad_hb[1] * nabla[v2][2];				\
            *(p++) = grad_hb[2] * nabla[v2][2];				\
            *(p++) = grad_hb[3] * nabla[v2][2];				\
	}								\
	if (p >= p0)							\
	    break;
	EDGE_BASED(v0, v1, v2);
	EDGE_BASED(v0, v2, v1);
	EDGE_BASED(v1, v2, v0);
#undef EDGE_BASED
	if (t->np_edge == 3 * hb.type->np_edge)
	    continue;
	/* face-based */
	/* get tangent vectors of the face */
	coord = dof->g->verts + e->verts[v1];
	a1 = (*coord)[0];
	b1 = (*coord)[1];
	c1 = (*coord)[2];
	coord = dof->g->verts + e->verts[v2];
	a2 = (*coord)[0];
	b2 = (*coord)[1];
	c2 = (*coord)[2];
	coord = dof->g->verts + e->verts[v0];
	a1 -= (*coord)[0];
	b1 -= (*coord)[1];
	c1 -= (*coord)[2];
	a2 -= (*coord)[0];
	b2 -= (*coord)[1];
	c2 -= (*coord)[2];
	/* normalize tangent vectors (not necessary?) */
	d = 1. / Sqrt(a1 * a1 + b1 * b1 + c1 * c1);
	a1 *= d;
	b1 *= d;
	c1 *= d;
	d = 1. / Sqrt(a2 * a2 + b2 * b2 + c2 * c2);
	a2 *= d;
	b2 *= d;
	c2 *= d;
	m = NVert + NEdge * hb.type->np_edge + i * hb.type->np_face;
	grad_hb = HB_grad(&hb, e, m, m + hb.type->np_face, lambda);
	for (j = 0; j < hb.type->np_face; j++, grad_hb += Dim + 1) {
	    if (i0++ >= no0) {
		*(p++) = grad_hb[0] * a1;
		*(p++) = grad_hb[1] * a1;
		*(p++) = grad_hb[2] * a1;
		*(p++) = grad_hb[3] * a1;

		*(p++) = grad_hb[0] * b1;
		*(p++) = grad_hb[1] * b1;
		*(p++) = grad_hb[2] * b1;
		*(p++) = grad_hb[3] * b1;

		*(p++) = grad_hb[0] * c1;
		*(p++) = grad_hb[1] * c1;
		*(p++) = grad_hb[2] * c1;
		*(p++) = grad_hb[3] * c1;

		if (p >= p0)
		    break;
	    }

	    if (i0++ >= no0) {
		*(p++) = grad_hb[0] * a2;
		*(p++) = grad_hb[1] * a2;
		*(p++) = grad_hb[2] * a2;
		*(p++) = grad_hb[3] * a2;

		*(p++) = grad_hb[0] * b2;
		*(p++) = grad_hb[1] * b2;
		*(p++) = grad_hb[2] * b2;
		*(p++) = grad_hb[3] * b2;

		*(p++) = grad_hb[0] * c2;
		*(p++) = grad_hb[1] * c2;
		*(p++) = grad_hb[2] * c2;
		*(p++) = grad_hb[3] * c2;

		if (p >= p0)
		    break;
	    }
	}
    }

    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return grad;
    }

    if ((no0 -= i0) < 0)
	no0 = 0;

    /*------------------------- element functions --------------------------*/
    if (DofIsHP(dof)) {
	order = order0;
	t = dof->hp->info->types[order];
    }
    n = t->np_elem;
    assert(n > 0);
    hb.type = DOF_HBn[order];
    GetElementVertices(e, 0, v0, v1, v2, v3, 
		       v4, v5, v6, v7);
    i0 = 0;
    /* face-based */
#define FACE_BASED(v)							\
    m = NVert + NEdge * hb.type->np_edge + v * hb.type->np_face;	\
    grad_hb = HB_grad(&hb, e, m, m + hb.type->np_face, lambda);		\
    for (j = 0; j < hb.type->np_face && p < p0;				\
		i0++, j++, grad_hb += Dim + 1) {			\
	if (i0 < no0)							\
	    continue;							\
	*(p++) = grad_hb[0] * nabla[v][0];				\
	*(p++) = grad_hb[1] * nabla[v][0];				\
	*(p++) = grad_hb[2] * nabla[v][0];				\
	*(p++) = grad_hb[3] * nabla[v][0];				\
									\
	*(p++) = grad_hb[0] * nabla[v][1];				\
	*(p++) = grad_hb[1] * nabla[v][1];				\
	*(p++) = grad_hb[2] * nabla[v][1];				\
	*(p++) = grad_hb[3] * nabla[v][1];				\
									\
	*(p++) = grad_hb[0] * nabla[v][2];				\
	*(p++) = grad_hb[1] * nabla[v][2];				\
	*(p++) = grad_hb[2] * nabla[v][2];				\
	*(p++) = grad_hb[3] * nabla[v][2];				\
    }									\
    if (p >= p0)							\
	goto cont;
    FACE_BASED(v0);
    FACE_BASED(v1);
    FACE_BASED(v2);
    FACE_BASED(v3);
#undef FACE_BASED
cont:
    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return grad;
    }
    /* element-based */
    m = NVert + NEdge * hb.type->np_edge + NFace * hb.type->np_face;
    grad_hb = HB_grad(&hb, e, m, m + hb.type->np_elem, lambda);
    for (j = 0; j < hb.type->np_elem; j++, grad_hb += Dim + 1) {
	if (i0++ >= no0) {
	    memcpy(p + 0 * (Dim + 1), grad_hb, (Dim + 1) * sizeof(*p));
	    p += Dim * (Dim + 1);
	    if (p >= p0)
	    break;
	}

	if (i0++ >= no0) {
	    memcpy(p + 1 * (Dim + 1), grad_hb, (Dim + 1) * sizeof(*p));
	    p += Dim * (Dim + 1);
	    if (p >= p0)
	    break;
	}

	if (i0++ >= no0) {
	    memcpy(p + 2 * (Dim + 1), grad_hb, (Dim + 1) * sizeof(*p));
	    p += Dim * (Dim + 1);
	    if (p >= p0)
	    break;
	}
    }

    assert(p == p0);

#if CACHE_LP
    cache_reset = TRUE;
#endif	/* CACHE_LP */
    return grad;
}

static const FLOAT *
face_normal(GRID *g, ELEMENT *e, int face)
{
    static FLOAT normal[3 * NFace];
#if USE_OMP
# pragma omp threadprivate(normal)
#endif	/* USE_OMP */

    FLOAT *n = normal + 3 * face;
    FLOAT x0, y0, z0, x1, y1, z1, x2, y2, z2;
    COORD *c;
    int u0, u1, u2, u3;

    if (IsLeaf(e)) {
	memcpy(n, phgGeomGetFaceNormal(g, e, face), 3 * sizeof(*n));
    }
    else {
	GetFaceVertices(e, face, u0, u1, u2, u3);
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
	/* FIXME: normalize n? */
    }
    /* FIXME: adjust direction of n? */

    return n;
}

static void
HC1_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	 DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	 const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
/* computes edge DOF by L2 projection on edge */
{
#if 0	/*********************************************************************/
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
	    for (j = i, q = p; j < n; j++)
		A[i * n + j] += *p * *(q++) * w0;
	    for (j = 0, q = temp1; j < m; j++)
		B[i * m + j] += *p * *(q++) * w0;
	}
    }
    for (i = 1; i < n; i++)
	for (j = 0; j < i; j++)
	    A[i * n + j] = A[j * n + i];
#if 0
printf("---------------------------------------------\n");
for (i = 0; i < n; i++) {
for (j = 0; j < n; j++)
if (Fabs(w0=A[i * n + j])<1e-7) printf("       "); else printf(" %0.2le", (double)w0);
printf(" = %lg\n", (double)B[i]);
}
printf("---------------------------------------------\n");
#endif

    /* solver the linear system of equations */
    if (!phgSolverDenseSolver(n, m, (void *)A, (void *)B))
	phgError(1, "(%s:%d) singular system (shouldn't happen).\n",
			__FILE__, __LINE__);
    phgFree(A);
    phgFree(temp0);
    phgFree(temp1);
    phgFree(fvals);
#else	/*********************************************************************/
    /* Note: for edge basis functions, the tangential component are
     * L2 orthogonal on the edge due to the fact that for any j
     *		tau . phi_j = L_j(xi)
     * where tau is the edge vector (v1-v0). The above relation also
     * implies that <tau . phi_j, tau . phi_j>_e is independent of
     * the geometry (e denotes the edge) */
    GRID *g = dof->g;
    FLOAT lam[] = {0., 0., 0., 0.};
    int i, j, k, n, u0, u1, n0, n1, m = dof->dim;
    FLOAT x0, y0, z0, x1, y1, z1, x2, y2, z2, x, y, z;
    COORD *c;
    QUAD *quad;
    FLOAT *fvals = phgAlloc(Dim * m * sizeof(*fvals));
    FLOAT *temp0, *temp1, *B = dofvalues;
    FLOAT w0, *r;
    const FLOAT *p, *q, *qw, *qp, *b;
    DOF_TYPE *t;
    static FLOAT *A = NULL;
    static int size_A = 0;
    static BOOLEAN fa_flag = FALSE;
#if USE_OMP
# pragma omp threadprivate(A, size_A, fa_flag)
#endif	/* USE_OMP */

    if (!fa_flag) {
	fa_flag = TRUE;
	FreeAtExitNoCheck(A);
    }
    Unused(pdofvalues);

    assert(type == EDGE);
    assert(funcvalues == NULL);
    assert(userfunc_lambda != NULL || userfunc != NULL);

    t = dof->type;
    if (t == NULL)
	t = dof->hp->info->types[dof->hp->edge_order[e->edges[index]]];
    n = t->np_edge;
    temp0 = phgAlloc(n * sizeof(*temp0));
    temp1 = phgAlloc(m * sizeof(*temp1));

    if (size_A < n) {
	A = phgRealloc_(A, n * sizeof(*A), size_A * sizeof(*A));
	for (i = size_A; i < n; i++)
	    A[i] = 0.;
    }

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

    quad = phgQuadGetQuad1D(2 * t->order);
    qp = quad->points;
    qw = quad->weights;
    /* range of basis functions on the current edge */
    if (!DofIsHP(dof)) {
	n0 = n * index;
    }
    else {
	n0 = 0;
	for (i = 0; i < index; i++)
	    n0 += dof->hp->info->types
			[dof->hp->edge_order[e->edges[i]]]->np_edge;
    }
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
	b = t->BasFuncs(dof, e, n0, n1, lam);
	/* inner product of the basis funcs with the tangent vector */
	for (p = b, r = temp0, i = 0; i < n; i++, p += 3)
	    *(r++) = x2 * p[0] + y2 * p[1] + z2 * p[2];
	/* coefficients and RHS */
	w0 = *(qw++);
	for (i = 0, p = temp0; i < n; i++, p++) {
	    if (i >= size_A)
		A[i] += *p * *p * w0;
	    for (j = 0, q = temp1; j < m; j++)
		B[i * m + j] += *p * *(q++) * w0;
	}
    }

    for (i = 0, p = A; i < n; i++, p++) {
	for (j = 0; j < m; j++)
	    *(B++) /= *p;
    }

    if (size_A < n)
	size_A = n;

    phgFree(temp0);
    phgFree(temp1);
    phgFree(fvals);
#endif	/*********************************************************************/
}

static void
HC2_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	 DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	 const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
/* computes face DOF by L2 projection on face */
{
    GRID *g = dof->g;
    FLOAT lam[] = {0., 0., 0., 0.};
    int i, j, k, n, u0, u1, u2, u3, n0, n1, eno, m = dof->dim;
    int e0, e1, e2, npe[NEdge];
    FLOAT x, y, z;
    QUAD *quad;
    FLOAT *fvals, *temp, *A, *B = dofvalues;
    FLOAT w0, *r;
    DOF_TYPE *t;
    const FLOAT *p, *q, *qw, *qp, *b, *nv;

    assert(type == FACE);
    assert(userfunc_lambda != NULL || userfunc != NULL);
    assert(funcvalues == NULL);

    t = dof->type;
    if (t == NULL)
	t = dof->hp->info->types[dof->hp->face_order[e->faces[index]]];
    n = t->np_face;
    if (n == 0)
	return;
    fvals = phgAlloc(Dim * m * sizeof(*fvals));
    temp = phgAlloc(Dim * n * sizeof(*temp));
    A = phgCalloc(n * n, sizeof(*A));

    for (i = 0; i < n * m; i++)
	B[i] = 0.0;

    /* range of basis functions on the current face */
    if (!DofIsHP(dof)) {
	for (i = 0; i < NEdge; i++)
	    npe[i] = t->np_edge;
	n0 = NEdge * t->np_edge + n * index;
    }
    else {
	n0 = 0;
	for (i = 0; i < NEdge; i++) {
	    npe[i] = dof->hp->info->types
				   [dof->hp->edge_order[e->edges[i]]]->np_edge;
	    n0 += npe[i];
	}
	for (i = 0; i < index; i++) {
	    n0 += dof->hp->info->types
			[dof->hp->face_order[e->faces[i]]]->np_face;
	}
    }
    n1 = n0 + n;

    /* the three vertices on the face */
    GetFaceVertices(e, index, u0, u1, u2, u3);
    /* the three edges on the face */
    e0 = GetEdgeNo(u0,u1);
    e1 = GetEdgeNo(u0,u2);
    e2 = GetEdgeNo(u1,u2);
    /* the normal vector */
    nv = face_normal(g, e, index);
    quad = phgQuadGetQuad2D(2 * t->order);
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
	b = t->BasFuncs(dof, e, 0, n1, lam);
	/* fvals -= contributions from edges. Note that at this stage
	 * dof->data already contains valid data on all edges */
	for (eno = 0, p = b; eno < NEdge; eno++) {
	    /* Note: only need to remove contributions of the
	     * basis functions on the three edges of the face */
	    if (eno != e0 && eno != e1 && eno != e2) {
		p += npe[eno] * 3;
		continue;
	    }
	    if (pdofvalues == NULL)
		q = DofEdgeData(dof, e->edges[eno]);
	    else
		q = pdofvalues[NVert + eno];
	    for (i = 0; i < npe[eno]; i++, p += 3) {
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
		A[i * n + j] += (p[0]*q[0] + p[1]*q[1] + p[2]*q[2]) * w0;
	    for (j = 0, q = fvals; j < m; j++, q += 3)
		B[i * m + j] += (p[0]*q[0] + p[1]*q[1] + p[2]*q[2]) * w0;
	}
    }

    for (i = 1; i < n; i++)
	for (j = 0; j < i; j++)
	    A[i * n + j] = A[j * n + i];
#if 0
printf("---------------------------------------------\n");
for (i = 0; i < n; i++) {
for (j = 0; j < n; j++) printf(" %lg",
				(double)(Fabs(A[i * n + j]) < 1e-10 ? 0. : 1.));
printf("\n");
}
printf("---------------------------------------------\n");
printf("\n");
#endif

    show_condition(n, A, type, index);

    /* solver the linear system of equations. */
    if (!phgSolverDenseSolver(n, m, (void *)A, (void *)B))
	phgError(1, "(%s:%d) singular system (shouldn't happen).\n",
			__FILE__, __LINE__);

    phgFree(A);
    phgFree(temp);
    phgFree(fvals);
}

static void
HC3_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	 DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	 const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
/* computes element DOF by L2 projection on element */
{
    GRID *g = dof->g;
    int i, j, k, n, n0, n1, no, m = dof->dim;
    int npe[NEdge], npf[NFace];
    QUAD *quad;
    FLOAT x, y, z, *fvals, *A, *B = dofvalues;
    FLOAT w0, *r;
    DOF_TYPE *t;
    const FLOAT *p, *q, *qw, *qp, *b;

    assert(type == VOLUME);
    assert(index == 0);
    assert(userfunc_lambda != NULL || userfunc != NULL);
    assert(funcvalues == NULL);

    t = dof->type;
    if (t == NULL)
	t = dof->hp->info->types[dof->hp->elem_order[e->index]];
    n = t->np_elem;
    if (n == 0)
	return;
    fvals = phgAlloc(Dim * m * sizeof(*fvals));
    A = phgCalloc(n * n, sizeof(*A));

    for (i = 0; i < n * m; i++)
	B[i] = 0.0;

    /* range of basis functions in the current element */
    if (!DofIsHP(dof)) {
	for (i = 0; i < NEdge; i++)
	    npe[i] = t->np_edge;
	for (i = 0; i < NFace; i++)
	    npf[i] = t->np_face;
	n0 = t->np_edge * NEdge + t->np_face * NFace;
    }
    else {
	n0 = 0;
	for (i = 0; i < NEdge; i++) {
	    npe[i] = dof->hp->info->types
				   [dof->hp->edge_order[e->edges[i]]]->np_edge;
	    n0 += npe[i];
	}
	for (i = 0; i < NFace; i++) {
	    npf[i] = dof->hp->info->types
			[dof->hp->face_order[e->faces[i]]]->np_face;
	    n0 += npf[i];
	}
    }
    n1 = n0 + n;

    quad = phgQuadGetQuad3D(2 * t->order);
    qp = quad->points;
    qw = quad->weights;
    for (k = 0; k < quad->npoints; k++, qp += 4) {
	if (userfunc_lambda != NULL) {
	    userfunc_lambda(dof, e, -1, qp, fvals);
	}
	else {
	    phgGeomLambda2XYZ(g, e, qp, &x, &y, &z);
	    userfunc(x, y, z, fvals);
	}
	/* evaluate basis functions, since we also need values of the
	 * basis functions on all edges/faces, the range is set to
	 * (0, n1) */
	b = t->BasFuncs(dof, e, 0, n1, qp);
	/* fvals -= contributions from edges. Note that at this stage
	 * dof->data already contains valid data on all edges */
	for (no = 0, p = b; no < NEdge; no++) {
	    if (pdofvalues == NULL)
		q = DofEdgeData(dof, e->edges[no]);
	    else
		q = pdofvalues[NVert + no];
	    for (i = 0; i < npe[no]; i++, p += 3) {
		for (j = 0, r = fvals; j < m; j++) {
		    w0 = *(q++);
		    *(r++) -= w0 * p[0];
		    *(r++) -= w0 * p[1];
		    *(r++) -= w0 * p[2];
		}
	    }
	}
	/* fvals -= contributions from faces. Note that at this stage
	 * dof->data already contains valid data on all faces */
	for (no = 0; no < NFace; no++) {
	    if (pdofvalues == NULL)
		q = DofFaceData(dof, e->faces[no]);
	    else
		q = pdofvalues[NVert + NEdge + no];
	    for (i = 0; i < npf[no]; i++, p += 3) {
		for (j = 0, r = fvals; j < m; j++) {
		    w0 = *(q++);
		    *(r++) -= w0 * p[0];
		    *(r++) -= w0 * p[1];
		    *(r++) -= w0 * p[2];
		}
	    }
	}
	b += 3 * n0;
	/* coefficients and RHS */
	w0 = *(qw++);
	for (i = 0, p = b; i < n; i++, p += 3) {
	    for (j = i, q = p; j < n; j++, q += 3)
		A[i * n + j] += (p[0] * q[0] + p[1] * q[1] + p[2] * q[2]) * w0;
	    for (j = 0, q = fvals; j < m; j++, q += 3) {
		B[i * m + j] += (p[0] * q[0] + p[1] * q[1] + p[2] * q[2]) * w0;
	    }
	}
    }

    for (i = 1; i < n; i++)
	for (j = 0; j < i; j++)
	    A[i * n + j] = A[j * n + i];
#if 0
printf("---------------------------------------------\n");
for (i = 0; i < n; i++) {
for (j = 0; j < n; j++)
printf(" %012.4le", (double)(Fabs(w0 = A[i * n + j]) < 1e-10 ? 0. : w0));
printf("\n");
}
printf("---------------------------------------------\n");
printf("\n");
#endif

    show_condition(n, A, type, index);

    /* solver the linear system of equations */
    if (!phgSolverDenseSolver(n, m, (void *)A, (void *)B))
	phgError(1, "(%s:%d) singular system (shouldn't happen).\n",
			__FILE__, __LINE__);

    phgFree(A);
    phgFree(fvals);
}

static void
HC_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	 DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	 const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
{
#if 0 * DEBUG
    static BOOLEAN warned = FALSE;
    if (phgRank == 0 && !warned && dof->type->orders == NULL) {
	warned = TRUE;
	phgWarning("the 'orders' member of \"%s\" is undefined!\n",
			dof->type->name);
    }
#endif

    switch (type) {
	case VERTEX:
	    /*phgError(1, "(%s:%d) unexpected error.\n", __FILE__, __LINE__);*/
	    break;
	case EDGE:
	    HC1_init(dof, e, type, index, userfunc, userfunc_lambda,
			funcvalues, dofvalues, pdofvalues);
	    break;
	case FACE:
	    HC2_init(dof, e, type, index, userfunc, userfunc_lambda,
			funcvalues, dofvalues, pdofvalues);
	    break;
	case VOLUME:
	    HC3_init(dof, e, type, index, userfunc, userfunc_lambda,
			funcvalues, dofvalues, pdofvalues);
	    break;
    }
}

/*-------------------------- DOF types HCn --------------------------*/

HP_INFO hp_info_HC = {
    DOF_HCn,		/* types */
    "Hierarchic Hcurl",	/* name */
    FE_Hcurl,		/* FE space */
    3,			/* dim */
    0,			/* min_order */
    15,			/* max_order */
    EDGE_FLAG | FACE_FLAG | ELEM_FLAG
};

#define NBasEdgeHC(p)	(p + 1)
#define NBasFaceHC(p)	(2 * NBasFace(p) + 3 * NBasEdge(p))
#define NBasElemHC(p)	(3 * NBasElem(p) + 4 * NBasFace(p))

#define DofTypeHC(order, type, name, orders, grad)			\
    DOF_TYPE type = {DofReserved, name, NULL, orders,			\
			grad, NULL, &hp_info_HC,			\
			phgDofInterpC2FGeneric,				\
			phgDofInterpF2CGeneric,				\
			HC_init, HC_bas, HC_grad, NULL, FE_Hcurl,	\
			FALSE, FALSE, FALSE, -1,			\
			NEdge * NBasEdgeHC(order) +			\
			NFace * NBasFaceHC(order) + NBasElemHC(order),	\
			(order == 0 ? 1 : order), /* order */		\
			0,		/* D_order */		\
			-1,		/* continuity */		\
			3,		/* dim */			\
			0, NBasEdgeHC(order),				\
			NBasFaceHC(order), NBasElemHC(order)};

static BYTE HC0_orders[] = {1};
DofTypeHC(0, DOF_HC0_, "HC0", HC0_orders, DOF_DG0)

static BYTE HC1_orders[] = {1, 1};
DofTypeHC(1, DOF_HC1_, "HC1", HC1_orders, DOF_DG0)

static BYTE HC2_orders[] = {1, 1, 2, 2, 2, 2};
DofTypeHC(2, DOF_HC2_, "HC2", HC2_orders, DOF_DG1)

static BYTE HC3_orders[] = {
    1, 1, 2, 3,
    2, 3, 2, 3, 2, 3, 3,
    3, 3, 3, 3};
DofTypeHC(3,  DOF_HC3_, "HC3", HC3_orders, DOF_DG2)

DofTypeHC(4, DOF_HC4_, "HC4", NULL, DOF_DG3)
DofTypeHC(5, DOF_HC5_, "HC5", NULL, DOF_DG4)
DofTypeHC(6, DOF_HC6_, "HC6", NULL, DOF_DG5)
DofTypeHC(7, DOF_HC7_, "HC7", NULL, DOF_DG6)
DofTypeHC(8, DOF_HC8_, "HC8", NULL, DOF_DG7)
DofTypeHC(9, DOF_HC9_, "HC9", NULL, DOF_DG8)
DofTypeHC(10, DOF_HC10_, "HC10", NULL, DOF_DG9)
DofTypeHC(11, DOF_HC11_, "HC11", NULL, DOF_DG10)
DofTypeHC(12, DOF_HC12_, "HC12", NULL, DOF_DG11)
DofTypeHC(13, DOF_HC13_, "HC13", NULL, DOF_DG12)
DofTypeHC(14, DOF_HC14_, "HC14", NULL, DOF_DG13)
DofTypeHC(15, DOF_HC15_, "HC15", NULL, DOF_DG14)

DOF_TYPE *DOF_HCn[] = {
    DOF_HC0, DOF_HC1, DOF_HC2, DOF_HC3, DOF_HC4,
    DOF_HC5, DOF_HC6, DOF_HC7, DOF_HC8, DOF_HC9,
    DOF_HC10, DOF_HC11, DOF_HC12, DOF_HC13, DOF_HC14, DOF_HC15,
    NULL
};

/*------------------------- Hdiv conforming bases -------------------------*/

/* Note: the H(div) bases actually used are from the following two papers:
 *
 *   1. Mark Ainsworth and Joe Coyle, Hierarchic finite element bases on
 * 	unstructured tetrahedral meshes, Int. J. Numer. Meth. Engng,
 * 	2003; 58:2103--2130 (DOI: 10.1002/nme.847)
 *
 *   2. Wei Cai, Jian Wu and Jianguo Xin, Divergence-free H(div)-conforming
 *	hierarchical bases for magnetohydrodynamics (MHD),
 *	Commun Math. Stat., 2013; 1:19--35 (DOI 10.1007/s40304-013-0003-9)
 *
 * see doc/PHG-Hdiv-bases.pdf for more details. */

static const FLOAT *
HD_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT *buffers[ORDER_MAX + 1];
    static BOOLEAN initialized = FALSE;
#if USE_OMP
# pragma omp threadprivate(buffers, initialized)
#endif  /* USE_OMP */

    FLOAT *L;
    const FLOAT (*nabla)[Dim + 1] = (void *)(phgGeomGetJacobian(dof->g, e));

    DOF_TYPE *t;
    int order0, order, nbas, bufferlen;
    FLOAT *bas, *p, *p0, a0, b0, c0, a1 = 0., b1 = 0., c1 = 0., a2, b2, c2, d;
    const FLOAT *bas_hb;
    DOF hb;
    COORD *coord;
    int i, j, i0, v0, v1, v2, v3, m, n;

    if (!initialized) {
	initialized = TRUE;
	bzero(buffers, sizeof(buffers));
    }

    t = dof->type;
    if (t != NULL) {
	order = t->order;
	bufferlen = nbas = t->nbas;
	order0 = order;
    }
    else {
	assert(DofIsHP(dof));
	/* FIXME: it is assumed that vert/edge/face orders <= elem order */
	order = order0 = dof->hp->elem_order[e->index];
	nbas = DofNBas(dof, e);
	bufferlen = dof->hp->info->types[order0]->nbas;
    }

    assert(order >= 1 && order <= ORDER_MAX);
    if (no1 <= 0)
	no1 = nbas;
 
    if (no1 <= no0)
	return NULL;

    if ((bas = buffers[order0]) == NULL) {
	buffers[order0] = phgAlloc(Dim * bufferlen * sizeof(*bas));
	FreeAtExitNoCheck(buffers[order0]);
	bas = buffers[order0];
    }

    p = bas;
    p0 = p + Dim * (no1 - no0);

#if CACHE_LP
    if (L0 == NULL) {
	L0 = phgAlloc(NEdge * sizeof(*L0));
	DL0 = phgAlloc(NEdge * sizeof(*DL0));
	cache_order = phgAlloc(NEdge * sizeof(*cache_order));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
	FreeAtExitNoCheck(cache_order);
    }
    cache_reset = FALSE;
    memset(cache_order, 0, NEdge * sizeof(*cache_order));
#else	/* CACHE_LP */
    if (L0 == NULL) {
	L0 = phgAlloc(3 * sizeof(*L0));
	DL0 = phgAlloc(3 * sizeof(*DL0));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
    }
#endif	/* CACHE_LP */

    hb.hp = NULL;

    /*------------------------- face functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NFace && p < p0; i++) {
	if (DofIsHP(dof)) {
	    order = dof->hp->face_order[e->faces[i]];
	    assert(order <= order0);
	    t = dof->hp->info->types[order];
	}
	assert(t->np_face > 0);
	n = i0 + t->np_face;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	GetFaceVertices(e, i, v0, v1, v2, v3);

	/* edge-based face functions */

#if CACHE_LP
# define LEGENDRE(v0, v1) \
	get_L(order0 - 1, v0, v1, lambda, &L, NULL)
#else	/* CACHE_LP */
# define LEGENDRE(v0, v1) \
	phgLegendreP(order - 1, lambda[v1] - lambda[v0], L = L0[0], NULL)
#endif	/* CACHE_LP */

/* --------- Bases by Ainsworth & Coyle (dependent!) --------- */
#define EDGE_BASED1(v0, v1, v2)						\
	/* (a1,b1,c1) = \lambda_v2 \grad\lambda_v0 x \grad\lambda_v1 */	\
	a1 = lambda[v2] * (a2 = nabla[v0][1] * nabla[v1][2] -		\
				nabla[v0][2] * nabla[v1][1]);		\
	b1 = lambda[v2] * (b2 = nabla[v0][2] * nabla[v1][0] -		\
				nabla[v0][0] * nabla[v1][2]);		\
	c1 = lambda[v2] * (c2 = nabla[v0][0] * nabla[v1][1] -		\
				nabla[v0][1] * nabla[v1][0]);		\
	LEGENDRE(v0, v1);						\
	for (j = 0; j < order && p < p0; j++, i0++) {			\
	    if (i0 < no0)						\
		continue;						\
            *(p++) = L[j] * a1;						\
            *(p++) = L[j] * b1;						\
            *(p++) = L[j] * c1;						\
	}								\
	if (p >= p0)							\
	    break;

/* ---------- 2. 1st kind bases by CAI Wei (TEST ONLY!) ---------- */
#define EDGE_BASED2(v0, v1, v2)						\
	/* (a0,b0,c0) = \lambda_v2 \grad\lambda_v0 x \grad\lambda_v1 */	\
	a0 = lambda[v2] * (a2 = nabla[v0][1] * nabla[v1][2] -		\
				nabla[v0][2] * nabla[v1][1]);		\
	b0 = lambda[v2] * (b2 = nabla[v0][2] * nabla[v1][0] -		\
				nabla[v0][0] * nabla[v1][2]);		\
	c0 = lambda[v2] * (c2 = nabla[v0][0] * nabla[v1][1] -		\
				nabla[v0][1] * nabla[v1][0]);		\
  	/* (a0,b0,c0) /= |\grad\lambda_v0 x \grad\lambda_v1| */		\
	d = 1.0 /*/ Sqrt(a2*a2 + b2*b2 + c2*c2)*/;			\
	a0 *= d;							\
	b0 *= d;							\
	c0 *= d;							\
	/* Note: (a1,b1,c1)=(p_{j-2},p_{j-1},p_j),			\
	 * where p_j=P^(3,0)_j(x), x = 2*\lambda_1/(1-\lambda_0) - 1,	\
	 * P^(a,b)_j denotes the j-th order Jacobi polynomial (a,b) */	\
	d = 1.0 - lambda[v0];						\
	x = 2.0 * lambda[v1] / d - 1.0;					\
	for (j = 0; j < order && p < p0; j++, i0++) {			\
	    switch (j) {						\
		case 0:							\
		    c1 = 1.0;						\
		    break;						\
		case 1:							\
		    b1 = c1;						\
		    c1 = (3.0 + 5.0 * x) * 0.5;				\
		    break;						\
		default:						\
		    a1 = b1;						\
		    b1 = c1;						\
		    c1 = ((9*(2*j+2)+(2*j+1)*(2*j+2)*(2*j+3)*x)*b1	\
			  - 2*(j+2)*(j-1)*(2*j+3)*a1)			\
			 / (2*j*(j+3)*(2*j+1));				\
		    break;						\
	    }								\
	    if (i0 >= no0) {						\
		*(p++) = c1 * a0;					\
		*(p++) = c1 * b0;					\
		*(p++) = c1 * c0;					\
	    }								\
	    a0 *= d;							\
	    b0 *= d;							\
	    c0 *= d;							\
	}								\
	if (p >= p0)							\
	    break;

/* ---------- 3. modified 2nd kind bases by ZHANG Linbo ---------- */
#define EDGE_BASED3(v0, v1, v2)						\
	/* (a0,b0,c0) = \lambda_v2 \grad\lambda_v0 x \grad\lambda_v1 */	\
	a0 = lambda[v2] * (a2 = nabla[v0][1] * nabla[v1][2] -		\
				nabla[v0][2] * nabla[v1][1]);		\
	b0 = lambda[v2] * (b2 = nabla[v0][2] * nabla[v1][0] -		\
				nabla[v0][0] * nabla[v1][2]);		\
	c0 = lambda[v2] * (c2 = nabla[v0][0] * nabla[v1][1] -		\
				nabla[v0][1] * nabla[v1][0]);		\
	if (order > 1) {						\
	    a1 = lambda[v2] * a0;					\
	    b1 = lambda[v2] * b0;					\
	    c1 = lambda[v2] * c0;					\
	}								\
	LEGENDRE(v0, v1);						\
	for (j = 0; j < order && p < p0; j++, i0++) {			\
	    if (i0 < no0)						\
		continue;						\
	    switch (j) {						\
		case 0:							\
        	    *(p++) = a0;					\
        	    *(p++) = b0;					\
		    *(p++) = c0;					\
		    break;						\
		case 1:							\
        	    *(p++) = a1;					\
        	    *(p++) = b1;					\
		    *(p++) = c1;					\
		    break;						\
		default:						\
	            *(p++) = L[j-1] * a1 + L[j-2] * a0;			\
	            *(p++) = L[j-1] * b1 + L[j-2] * b0;			\
	            *(p++) = L[j-1] * c1 + L[j-2] * c0;			\
		    break;						\
	    }								\
	}								\
	if (p >= p0)							\
	    break;

/* ---------- 4. 2nd kind bases by CAI Wei ---------- */
#define EDGE_BASED4(v0, v1, v2)						\
	/* (a0,b0,c0) = \lambda_v2 \grad\lambda_v0 x \grad\lambda_v1 */	\
	a0 = lambda[v2] * (a2 = nabla[v0][1] * nabla[v1][2] -		\
				nabla[v0][2] * nabla[v1][1]);		\
	b0 = lambda[v2] * (b2 = nabla[v0][2] * nabla[v1][0] -		\
				nabla[v0][0] * nabla[v1][2]);		\
	c0 = lambda[v2] * (c2 = nabla[v0][0] * nabla[v1][1] -		\
				nabla[v0][1] * nabla[v1][0]);		\
/* (a1,b1,c1)=\lambda_v0\lambda_v1\grad\lambda_v2 x \grad\lambda_v0 */	\
	if (order > 1) {						\
	    FLOAT d = lambda[v0] * lambda[v1];				\
	    a1 = d * (nabla[v2][1] * nabla[v0][2] -			\
		      nabla[v2][2] * nabla[v0][1]);			\
	    b1 = d * (nabla[v2][2] * nabla[v0][0] -			\
		      nabla[v2][0] * nabla[v0][2]);			\
	    c1 = d * (nabla[v2][0] * nabla[v0][1] -			\
		      nabla[v2][1] * nabla[v0][0]);			\
	}								\
	LEGENDRE(v0, v1);						\
	for (j = 0; j < order && p < p0; j++, i0++) {			\
	    if (i0 < no0)						\
		continue;						\
	    switch (j) {						\
		case 0:							\
        	    *(p++) = a0;					\
        	    *(p++) = b0;					\
		    *(p++) = c0;					\
		    break;						\
		case 1:							\
        	    *(p++) = a1;					\
        	    *(p++) = b1;					\
		    *(p++) = c1;					\
		    break;						\
		default:						\
	            *(p++) = L[j-1] * a1 + L[j-2] * a0;			\
	            *(p++) = L[j-1] * b1 + L[j-2] * b0;			\
	            *(p++) = L[j-1] * c1 + L[j-2] * c0;			\
		    break;						\
	    }								\
	}								\
	if (p >= p0)							\
	    break;

#if 0	/* --------------------------- Bases by Ainsworth & Coyle (dependent) */
# define EDGE_BASIS_TYPE	1
# define EDGE_BASED(v0, v1, v2)	EDGE_BASED1(v0, v1, v2)
#elif 0	/* --------------------------- 1st kind bases by CAI Wei (TEST ONLY!) */
	FLOAT x;
	Unused(L);
# define EDGE_BASIS_TYPE	2
# define EDGE_BASED(v0, v1, v2)	EDGE_BASED2(v0, v1, v2)
#elif 0	/* --------------------------- modified 2nd kind bases by ZHANG Linbo */
# define EDGE_BASIS_TYPE	3
# define EDGE_BASED(v0, v1, v2)	EDGE_BASED3(v0, v1, v2)
#else	/* --------------------------- 2nd kind bases by CAI Wei */
# define EDGE_BASIS_TYPE	4
# define EDGE_BASED(v0, v1, v2) EDGE_BASED4(v0, v1, v2)
#endif	/* ------------------------------------------------------------------ */

	EDGE_BASED(v0, v1, v2);
	EDGE_BASED(v0, v2, v1);
	EDGE_BASED(v1, v2, v0);

#undef EDGE_BASED1
#undef EDGE_BASED2
#undef EDGE_BASED3
#undef EDGE_BASED4
#undef EDGE_BASED
#undef LEGENDRE

	/* face bubble functions */
	if (order < 3)
	    continue;
	hb.type = DOF_HBn[order];
	m = NVert + NEdge * hb.type->np_edge + i * hb.type->np_face;
	bas_hb = HB_bas(&hb, e, m, m + hb.type->np_face, lambda);
	for (j = 0; j < hb.type->np_face && p < p0; j++, bas_hb++, i0++) {
	    if (i0 < no0)
		continue;
	    *(p++) = *bas_hb * a2;
	    *(p++) = *bas_hb * b2;
	    *(p++) = *bas_hb * c2;
	}
    }

    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return bas;
    }

    if ((no0 -= i0) < 0)
	no0 = 0;

    /*------------------------- element functions --------------------------*/
    if (DofIsHP(dof)) {
	order = order0;
	t = dof->hp->info->types[order];
    }
    assert(order >= 2);
    hb.type = DOF_HBn[order];
    i0 = 0;

    /* edge-based element functions */
    for (i = 0; i < NEdge && p < p0; i++) {
	n = i0 + hb.type->np_edge;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	GetEdgeVertices(e, i, v0, v1);
	/* get edge vector (v0, v1) */
	coord = dof->g->verts + e->verts[v1];
	a1 = (*coord)[0];
	b1 = (*coord)[1];
	c1 = (*coord)[2];
	coord = dof->g->verts + e->verts[v0];
	a1 -= (*coord)[0];
	b1 -= (*coord)[1];
	c1 -= (*coord)[2];
	/* normalize edge vector (not necessary?) */
	d = 1. / Sqrt(a1 * a1 + b1 * b1 + c1 * c1);
	a1 *= d;
	b1 *= d;
	c1 *= d;
	m = NVert + i * hb.type->np_edge;
	bas_hb = HB_bas(&hb, e, m, m + hb.type->np_edge, lambda);
	for (j = 0; j < hb.type->np_edge && p < p0; j++, bas_hb++, i0++) {
	    if (i0 < no0)
		continue;
	    *(p++) = *bas_hb * a1;
	    *(p++) = *bas_hb * b1;
	    *(p++) = *bas_hb * c1;
	}
    }

    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return bas;
    }

    /* face-based element functions */
    assert(order >= 2);
    for (i = 0; i < NFace && p < p0; i++) {
	GetFaceVertices(e, i, v0, v1, v2, v3);
	/* get tangent vectors on the face */
	coord = dof->g->verts + e->verts[v1];
	a1 = (*coord)[0];
	b1 = (*coord)[1];
	c1 = (*coord)[2];
	coord = dof->g->verts + e->verts[v2];
	a2 = (*coord)[0];
	b2 = (*coord)[1];
	c2 = (*coord)[2];
	coord = dof->g->verts + e->verts[v0];
	a1 -= (*coord)[0];
	b1 -= (*coord)[1];
	c1 -= (*coord)[2];
	a2 -= (*coord)[0];
	b2 -= (*coord)[1];
	c2 -= (*coord)[2];
	/* normalize tangent vector (not necessary?) */
	d = 1. / Sqrt(a1 * a1 + b1 * b1 + c1 * c1);
	a1 *= d;
	b1 *= d;
	c1 *= d;
	d = 1. / Sqrt(a2 * a2 + b2 * b2 + c2 * c2);
	a2 *= d;
	b2 *= d;
	c2 *= d;
	m = NVert + NEdge * hb.type->np_edge + i * hb.type->np_face;
	bas_hb = HB_bas(&hb, e, m, m + hb.type->np_face, lambda);
	for (j = 0; j < hb.type->np_face; j++, bas_hb++) {
	    if (i0++ >= no0) {
		*(p++) = *bas_hb * a1;
		*(p++) = *bas_hb * b1;
		*(p++) = *bas_hb * c1;
		if (p >= p0)
		    break;
	    }

	    if (i0++ >= no0) {
		*(p++) = *bas_hb * a2;
		*(p++) = *bas_hb * b2;
		*(p++) = *bas_hb * c2;
		if (p >= p0)
		    break;
	    }
	}
    }

    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return bas;
    }

    /* element bubble functions */
    m = NVert + NEdge * hb.type->np_edge + NFace * hb.type->np_face;
    bas_hb = HB_bas(&hb, e, m, m + hb.type->np_elem, lambda);
    for (j = 0; j < hb.type->np_elem; j++, bas_hb++) {
	if (i0++ >= no0) {
	    *(p++) = *bas_hb;
	    *(p++) = 0.;
	    *(p++) = 0.;
	    if (p >= p0)
		break;
	}

	if (i0++ >= no0) {
	    *(p++) = 0.;
	    *(p++) = *bas_hb;
	    *(p++) = 0.;
	    if (p >= p0)
		break;
	}

	if (i0++ >= no0) {
	    *(p++) = 0.;
	    *(p++) = 0.;
	    *(p++) = *bas_hb;
	    if (p >= p0)
		break;
	}
    }

    assert(p == p0);

#if CACHE_LP
    cache_reset = TRUE;
#endif	/* CACHE_LP */

    return bas;
}

static const FLOAT *
HD_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT *buffers[ORDER_MAX + 1];
    static BOOLEAN initialized = FALSE;
#if USE_OMP
# pragma omp threadprivate(buffers, initialized)
#endif  /* USE_OMP */

    FLOAT *L, *DL;
    const FLOAT (*nabla)[Dim + 1] = (void *)(phgGeomGetJacobian(dof->g, e));

    DOF_TYPE *t;
    int order, order0, nbas, bufferlen;
    FLOAT *grad, *p, *p0, a0, b0, c0, a1 = 0, b1 = 0, c1 = 0, a2, b2, c2, d;
    const FLOAT *grad_hb;
    DOF hb;
    COORD *coord;
    int i, j, i0, v0, v1, v2, v3, m, n;

    if (!initialized) {
	initialized = TRUE;
	bzero(buffers, sizeof(buffers));
    }

    t = dof->type;
    if (t != NULL) {
	order = t->order;
	bufferlen = nbas = t->nbas;
	order0 = order;
    }
    else {
	assert(DofIsHP(dof));
	/* FIXME: it is assumed that vert/edge/face orders <= elem order */
	order = order0 = dof->hp->elem_order[e->index];
	nbas = DofNBas(dof, e);
	bufferlen = dof->hp->info->types[order0]->nbas;
    }

    assert(order >= 0 && order <= ORDER_MAX);
    if (no1 <= 0)
	no1 = nbas;
 
    if (no1 <= no0)
	return NULL;

    if ((grad = buffers[order]) == NULL) {
	buffers[order] = phgAlloc((Dim + 1) * Dim * bufferlen * sizeof(*grad));
	FreeAtExitNoCheck(buffers[order]);
	grad = buffers[order];
    }

    p = grad;
    p0 = p + (Dim + 1) * Dim * (no1 - no0);
    bzero(p, (no1 - no0) * (Dim + 1) * Dim * sizeof(*p));

#if CACHE_LP
    if (L0 == NULL) {
	L0 = phgAlloc(NEdge * sizeof(*L0));
	DL0 = phgAlloc(NEdge * sizeof(*DL0));
	cache_order = phgAlloc(NEdge * sizeof(*cache_order));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
	FreeAtExitNoCheck(cache_order);
    }
    cache_reset = FALSE;
    memset(cache_order, 0, NEdge * sizeof(*cache_order));
#else	/* CACHE_LP */
    if (L0 == NULL) {
	L0 = phgAlloc(3 * sizeof(*L0));
	DL0 = phgAlloc(3 * sizeof(*DL0));
	FreeAtExitNoCheck(L0);
	FreeAtExitNoCheck(DL0);
    }
#endif	/* CACHE_LP */

    hb.hp = NULL;

    /*------------------------- face functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NFace && p < p0; i++) {
	if (DofIsHP(dof)) {
	    order = dof->hp->face_order[e->faces[i]];
	    assert(order <= order0);
	    t = dof->hp->info->types[order];
	}
	n = i0 + t->np_face;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	GetFaceVertices(e, i, v0, v1, v2, v3);

	/* edge-based */

#if CACHE_LP
# define LEGENDRE(v0, v1) \
	get_L(order0 - 1, v0, v1, lambda, &L, &DL)
#else	/* CACHE_LP */
# define LEGENDRE(v0, v1) \
	phgLegendreP(order - 1, lambda[v1] - lambda[v0], L = L0[0], DL = DL0[0])
#endif	/* CACHE_LP */

/* ZHANG's modification of CAI's 2nd kind bases */
#define EDGE_BASED3(v0, v1, v2)						\
	/* (a0,b0,c0) = \lambda_v2 \grad\lambda_v0 x \grad\lambda_v1 */	\
	a0 = lambda[v2] * (a2 = nabla[v0][1] * nabla[v1][2] -		\
				nabla[v0][2] * nabla[v1][1]);		\
	b0 = lambda[v2] * (b2 = nabla[v0][2] * nabla[v1][0] -		\
				nabla[v0][0] * nabla[v1][2]);		\
	c0 = lambda[v2] * (c2 = nabla[v0][0] * nabla[v1][1] -		\
				nabla[v0][1] * nabla[v1][0]);		\
	if (order > 1) {						\
	    a1 = lambda[v2] * a0;					\
	    b1 = lambda[v2] * b0;					\
	    c1 = lambda[v2] * c0;					\
	}								\
	LEGENDRE(v0, v1);						\
	for (j = 0; j < order && p < p0; j++, i0++) {			\
	    if (i0 < no0)						\
		continue;						\
									\
	    switch (j) {						\
		case 0:							\
        	    p[v0] = 0.0;					\
        	    p[v1] = 0.0;					\
        	    p[v2] = a2;						\
		    p += Dim + 1;					\
									\
        	    p[v0] = 0.0;					\
        	    p[v1] = 0.0;					\
        	    p[v2] = b2;						\
		    p += Dim + 1;					\
									\
        	    p[v0] = 0.0;					\
        	    p[v1] = 0.0;					\
        	    p[v2] = c2;						\
		    p += Dim + 1;					\
		    break;						\
		case 1:							\
        	    p[v0] = 0.0;					\
        	    p[v1] = 0.0;					\
        	    p[v2] = 2.0 * a0;					\
		    p += Dim + 1;					\
									\
        	    p[v0] = 0.0;					\
        	    p[v1] = 0.0;					\
        	    p[v2] = 2.0 * b0;					\
		    p += Dim + 1;					\
									\
        	    p[v0] = 0.0;					\
        	    p[v1] = 0.0;					\
        	    p[v2] = 2.0 * c0;					\
		    p += Dim + 1;					\
		    break;						\
		default:						\
        	    p[v0] = - DL[j-1] * a1 - DL[j-2] * a0;		\
        	    p[v1] = -p[v0];					\
        	    p[v2] = 2.0 * L[j-1] * a0 + L[j-2] * a2;		\
		    p += Dim + 1;					\
									\
        	    p[v0] = - DL[j-1] * b1 - DL[j-2] * b0;		\
        	    p[v1] = -p[v0];					\
        	    p[v2] = 2.0 * L[j-1] * b0 + L[j-2] * b2;		\
		    p += Dim + 1;					\
									\
        	    p[v0] = - DL[j-1] * c1 - DL[j-2] * c0;		\
        	    p[v1] = -p[v0];					\
        	    p[v2] = 2.0 * L[j-1] * c0 + L[j-2] * c2;		\
		    p += Dim + 1;					\
		    break;						\
	    }								\
	}								\
	if (p >= p0)							\
	    break;

/* CAI's 2nd kind bases */
#define EDGE_BASED4(v0, v1, v2)						\
    {									\
	FLOAT a3 = 0., b3 = 0., c3 = 0.;				\
	/* (a0,b0,c0) = \lambda_v2 \grad\lambda_v0 x \grad\lambda_v1 */	\
	a0 = lambda[v2] * (a2 = nabla[v0][1] * nabla[v1][2] -		\
				nabla[v0][2] * nabla[v1][1]);		\
	b0 = lambda[v2] * (b2 = nabla[v0][2] * nabla[v1][0] -		\
				nabla[v0][0] * nabla[v1][2]);		\
	c0 = lambda[v2] * (c2 = nabla[v0][0] * nabla[v1][1] -		\
				nabla[v0][1] * nabla[v1][0]);		\
	if (order > 1) {						\
	    FLOAT d = lambda[v0] * lambda[v1];				\
/* (a1,b1,c1)=\lambda_v0\lambda_v1\grad\lambda_v2 x \grad\lambda_v0 */	\
	    a1 = d * (a3 = nabla[v2][1] * nabla[v0][2] -		\
				nabla[v2][2] * nabla[v0][1]);		\
	    b1 = d * (b3 = nabla[v2][2] * nabla[v0][0] -		\
				nabla[v2][0] * nabla[v0][2]);		\
	    c1 = d * (c3 = nabla[v2][0] * nabla[v0][1] -		\
				nabla[v2][1] * nabla[v0][0]);		\
	}								\
	LEGENDRE(v0, v1);						\
	for (j = 0; j < order && p < p0; j++, i0++) {			\
	    if (i0 < no0)						\
		continue;						\
									\
	    switch (j) {						\
		case 0:							\
        	    p[v0] = 0.0;					\
        	    p[v1] = 0.0;					\
        	    p[v2] = a2;						\
		    p += Dim + 1;					\
									\
        	    p[v0] = 0.0;					\
        	    p[v1] = 0.0;					\
        	    p[v2] = b2;						\
		    p += Dim + 1;					\
									\
        	    p[v0] = 0.0;					\
        	    p[v1] = 0.0;					\
        	    p[v2] = c2;						\
		    p += Dim + 1;					\
		    break;						\
		case 1:							\
        	    p[v0] = lambda[v1] * a3;				\
        	    p[v1] = lambda[v0] * a3;				\
        	    p[v2] = 0.0;					\
		    p += Dim + 1;					\
									\
        	    p[v0] = lambda[v1] * b3;				\
        	    p[v1] = lambda[v0] * b3;				\
        	    p[v2] = 0.0;					\
		    p += Dim + 1;					\
									\
        	    p[v0] = lambda[v1] * c3;				\
        	    p[v1] = lambda[v0] * c3;				\
        	    p[v2] = 0.0;					\
		    p += Dim + 1;					\
		    break;						\
		default:						\
        	    p[v0] = -DL[j-1] * a1 				\
			    +L[j-1] * lambda[v1] * a3			\
			    -DL[j-2] * a0;				\
        	    p[v1] = +DL[j-1] * a1				\
			    +L[j-1] * lambda[v0] * a3			\
			    +DL[j-2] * a0;				\
        	    p[v2] = L[j-2] * a2;				\
		    p += Dim + 1;					\
									\
        	    p[v0] = -DL[j-1] * b1 				\
			    +L[j-1] * lambda[v1] * b3			\
			    -DL[j-2] * b0;				\
        	    p[v1] = +DL[j-1] * b1				\
			    +L[j-1] * lambda[v0] * b3			\
			    +DL[j-2] * b0;				\
        	    p[v2] = L[j-2] * b2;				\
		    p += Dim + 1;					\
									\
        	    p[v0] = -DL[j-1] * c1 				\
			    +L[j-1] * lambda[v1] * c3			\
			    -DL[j-2] * c0;				\
        	    p[v1] = +DL[j-1] * c1				\
			    +L[j-1] * lambda[v0] * c3			\
			    +DL[j-2] * c0;				\
        	    p[v2] = L[j-2] * c2;				\
		    p += Dim + 1;					\
		    break;						\
	    }								\
	}								\
    }									\
	if (p >= p0)							\
	    break;

#if (EDGE_BASIS_TYPE == 3)
# define EDGE_BASED(v0, v1, v2)	EDGE_BASED3(v0, v1, v2)
#elif (EDGE_BASIS_TYPE == 4)
# define EDGE_BASED(v0, v1, v2)	EDGE_BASED4(v0, v1, v2)
#else
# define EDGE_BASED(v0, v1, v2)	\
	phgWarning("WARNING: %s: EDGE_BASIS_TYPE %d unsupported.\n",	\
			__func__, EDGE_BASIS_TYPE);			\
	EDGE_BASED4(v0, v1, v2)
#endif

	EDGE_BASED(v0, v1, v2);
	EDGE_BASED(v0, v2, v1);
	EDGE_BASED(v1, v2, v0);

#undef EDGE_BASED3
#undef EDGE_BASED4
#undef EDGE_BASED
#undef LEGENDRE
#undef EDGE_BASIS_TYPE

	/* face bubble functions */
	if (order < 3)
	    continue;
	hb.type = DOF_HBn[order];
	m = NVert + NEdge * hb.type->np_edge + i * hb.type->np_face;
	grad_hb = HB_grad(&hb, e, m, m + hb.type->np_face, lambda);
	for (j = 0; j < hb.type->np_face && p < p0;
			j++, grad_hb += Dim + 1, i0++) {
	    if (i0 < no0)
		continue;

	    *(p++) = grad_hb[0] * a2;
	    *(p++) = grad_hb[1] * a2;
	    *(p++) = grad_hb[2] * a2;
	    *(p++) = grad_hb[3] * a2;

	    *(p++) = grad_hb[0] * b2;
	    *(p++) = grad_hb[1] * b2;
	    *(p++) = grad_hb[2] * b2;
	    *(p++) = grad_hb[3] * b2;

	    *(p++) = grad_hb[0] * c2;
	    *(p++) = grad_hb[1] * c2;
	    *(p++) = grad_hb[2] * c2;
	    *(p++) = grad_hb[3] * c2;
	}
    }

    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return grad;
    }

    if ((no0 -= i0) < 0)
	no0 = 0;

    /*------------------------- element functions --------------------------*/
    if (DofIsHP(dof)) {
	order = order0;
	t = dof->hp->info->types[order];
    }
    assert(order >= 2);
    hb.type = DOF_HBn[order];
    i0 = 0;

    /* edge-based element functions */
    for (i = 0; i < NEdge && p < p0; i++) {
	n = i0 + hb.type->np_edge;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	GetEdgeVertices(e, i, v0, v1);
	/* get edge vector (v0, v1) */
	coord = dof->g->verts + e->verts[v1];
	a1 = (*coord)[0];
	b1 = (*coord)[1];
	c1 = (*coord)[2];
	coord = dof->g->verts + e->verts[v0];
	a1 -= (*coord)[0];
	b1 -= (*coord)[1];
	c1 -= (*coord)[2];
	/* normalize edge vector (not necessary?) */
	d = 1. / Sqrt(a1 * a1 + b1 * b1 + c1 * c1);
	a1 *= d;
	b1 *= d;
	c1 *= d;
	m = NVert + i * hb.type->np_edge;
	grad_hb = HB_grad(&hb, e, m, m + hb.type->np_edge, lambda);
	for (j = 0; j < hb.type->np_edge && p < p0;
				j++, grad_hb += Dim + 1, i0++) {
	    if (i0 < no0)
		continue;

	    *(p++) = grad_hb[0] * a1;
	    *(p++) = grad_hb[1] * a1;
	    *(p++) = grad_hb[2] * a1;
	    *(p++) = grad_hb[3] * a1;

	    *(p++) = grad_hb[0] * b1;
	    *(p++) = grad_hb[1] * b1;
	    *(p++) = grad_hb[2] * b1;
	    *(p++) = grad_hb[3] * b1;

	    *(p++) = grad_hb[0] * c1;
	    *(p++) = grad_hb[1] * c1;
	    *(p++) = grad_hb[2] * c1;
	    *(p++) = grad_hb[3] * c1;
	}
    }

    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return grad;
    }

    /* face-based element functions */
    assert(order >= 2);
    for (i = 0; i < NFace && p < p0; i++) {
	GetFaceVertices(e, i, v0, v1, v2, v3);
	/* get tangent vectors on the face */
	coord = dof->g->verts + e->verts[v1];
	a1 = (*coord)[0];
	b1 = (*coord)[1];
	c1 = (*coord)[2];
	coord = dof->g->verts + e->verts[v2];
	a2 = (*coord)[0];
	b2 = (*coord)[1];
	c2 = (*coord)[2];
	coord = dof->g->verts + e->verts[v0];
	a1 -= (*coord)[0];
	b1 -= (*coord)[1];
	c1 -= (*coord)[2];
	a2 -= (*coord)[0];
	b2 -= (*coord)[1];
	c2 -= (*coord)[2];
	/* normalize tangent vector (not necessary?) */
	d = 1. / Sqrt(a1 * a1 + b1 * b1 + c1 * c1);
	a1 *= d;
	b1 *= d;
	c1 *= d;
	d = 1. / Sqrt(a2 * a2 + b2 * b2 + c2 * c2);
	a2 *= d;
	b2 *= d;
	c2 *= d;
	m = NVert + NEdge * hb.type->np_edge + i * hb.type->np_face;
	grad_hb = HB_grad(&hb, e, m, m + hb.type->np_face, lambda);
	for (j = 0; j < hb.type->np_face; j++, grad_hb += Dim + 1) {
	    if (i0++ >= no0) {
		*(p++) = grad_hb[0] * a1;
		*(p++) = grad_hb[1] * a1;
		*(p++) = grad_hb[2] * a1;
		*(p++) = grad_hb[3] * a1;

		*(p++) = grad_hb[0] * b1;
		*(p++) = grad_hb[1] * b1;
		*(p++) = grad_hb[2] * b1;
		*(p++) = grad_hb[3] * b1;

		*(p++) = grad_hb[0] * c1;
		*(p++) = grad_hb[1] * c1;
		*(p++) = grad_hb[2] * c1;
		*(p++) = grad_hb[3] * c1;

		if (p >= p0)
		    break;
	    }

	    if (i0++ >= no0) {
		*(p++) = grad_hb[0] * a2;
		*(p++) = grad_hb[1] * a2;
		*(p++) = grad_hb[2] * a2;
		*(p++) = grad_hb[3] * a2;

		*(p++) = grad_hb[0] * b2;
		*(p++) = grad_hb[1] * b2;
		*(p++) = grad_hb[2] * b2;
		*(p++) = grad_hb[3] * b2;

		*(p++) = grad_hb[0] * c2;
		*(p++) = grad_hb[1] * c2;
		*(p++) = grad_hb[2] * c2;
		*(p++) = grad_hb[3] * c2;

		if (p >= p0)
		    break;
	    }
	}
    }

    if (p >= p0) {
#if CACHE_LP
	cache_reset = TRUE;
#endif	/* CACHE_LP */
	return grad;
    }

    /* element bubble functions */
    m = NVert + NEdge * hb.type->np_edge + NFace * hb.type->np_face;
    grad_hb = HB_grad(&hb, e, m, m + hb.type->np_elem, lambda);
    for (j = 0; j < hb.type->np_elem; j++, grad_hb += Dim + 1) {
	if (i0++ >= no0) {
	    memcpy(p + 0 * (Dim + 1), grad_hb, (Dim + 1) * sizeof(*grad_hb));
	    p += Dim * (Dim + 1);
	    if (p >= p0)
		break;
	}

	if (i0++ >= no0) {
	    memcpy(p + 1 * (Dim + 1), grad_hb, (Dim + 1) * sizeof(*grad_hb));
	    p += Dim * (Dim + 1);
	    if (p >= p0)
		break;
	}

	if (i0++ >= no0) {
	    memcpy(p + 2 * (Dim + 1), grad_hb, (Dim + 1) * sizeof(*grad_hb));
	    p += Dim * (Dim + 1);
	    if (p >= p0)
		break;
	}
    }

    assert(p == p0);

#if CACHE_LP
    cache_reset = TRUE;
#endif	/* CACHE_LP */

    return grad;
}

static void
HD2_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	 DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	 const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
/* computes face DOF by L2 projection on face */
{
    GRID *g = dof->g;
    FLOAT lam[] = {0., 0., 0., 0.};
    int i, j, k, n, u0, u1, u2, u3, n0, n1, m = dof->dim;
    FLOAT x, y, z;
    QUAD *quad;
    FLOAT *fvals, *temp, *A, *B = dofvalues;
    FLOAT w0, *r;
    DOF_TYPE *t;
    const FLOAT *p, *q, *qw, *qp, *b, *nv;

    assert(type == FACE);
    assert(userfunc_lambda != NULL || userfunc != NULL);
    assert(funcvalues == NULL);

    t = dof->type;
    if (t == NULL)
	t = dof->hp->info->types[dof->hp->face_order[e->faces[index]]];
    n = t->np_face;
    if (n == 0)
	return;
    fvals = phgAlloc(Dim * m * sizeof(*fvals));
    temp = phgAlloc(Dim * n * sizeof(*temp));
    A = phgCalloc(n * n, sizeof(*A));

    for (i = 0; i < n * m; i++)
	B[i] = 0.0;

    /* range of basis functions on the current face */
    if (!DofIsHP(dof)) {
	n0 = n * index;
    }
    else {
	n0 = 0;
	for (i = 0; i < index; i++) {
	    n0 += dof->hp->info->types
			[dof->hp->face_order[e->faces[i]]]->np_face;
	}
    }
    n1 = n0 + n;

    /* the three vertices on the face */
    GetFaceVertices(e, index, u0, u1, u2, u3);

    /* the normal vector */
    nv = face_normal(g, e, index);
    quad = phgQuadGetQuad2D(2 * t->order);
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
	/* evaluate basis functions */
	b = t->BasFuncs(dof, e, n0, n1, lam);
#if 1	/* project using normal components (the trace on the face) */
	/* inner product of fvals with the face normal */
	for (i = 0, p = r = fvals; i < m; i++, p += 3)
	    (*r++) = p[0] * nv[0] + p[1] * nv[1] + p[2] * nv[2];
	/* inner product of basis functions with the face normal */
	for (i = 0, r = temp, p = b; i < n; i++, p += 3)
	    *(r++) = p[0] * nv[0] + p[1] * nv[1] + p[2] * nv[2];
	/* coefficients and RHS */
	w0 = *(qw++);
	for (i = 0, p = temp; i < n; i++, p++) {
	    for (j = i, q = p; j < n; j++, q++)
		A[i * n + j] += (*p) * (*q) * w0;
	    for (j = 0, q = fvals; j < m; j++, q++)
		B[i * m + j] += (*p) * (*q) * w0;
	}
#else	/* project all components (won't work, for test only!) */
	Unused(nv); Unused(r);
	w0 = *(qw++);
	for (i = 0, p = b; i < n; i++, p += 3) {
	    for (j = i, q = p; j < n; j++, q += 3)
		A[i * n + j] += (p[0] * q[0] + p[1] * q[1] + p[2] * q[2]) * w0;
	    for (j = 0, q = fvals; j < m; j++, q += 3)
		B[i * m + j] += (p[0] * q[0] + p[1] * q[1] + p[2] * q[2]) * w0;
	}
#endif
    }

    for (i = 1; i < n; i++)
	for (j = 0; j < i; j++)
	    A[i * n + j] = A[j * n + i];

#if 0
printf("\nface %d\n", e->faces[index]);
printf("--------------------------------------------- n=%d, m=%d\n", n, m);
printf("A = [");
for (i = 0; i < n; i++) {
  for (j = 0; j < n; j++)
    printf(" %0.16lg", 0.5 * (double)A[i * n + j]);
  printf("\n");
}
printf("]\nB = [");
for (i = 0; i < n; i++)
  printf("%s%lg", i == 0 ? "" : "; ", 0.5 * (double)B[i]);
printf("]\n");
printf("---------------------------------------------\n");
printf("\n");
#endif

    show_condition(n, A, type, index);

    /* solver the linear system of equations. */
    if (!phgSolverDenseSolver(n, m, (void *)A, (void *)B))
	phgError(1, "(%s:%d) singular system (shouldn't happen).\n",
			__FILE__, __LINE__);

    phgFree(A);
    phgFree(temp);
    phgFree(fvals);
}

static void
HD3_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	 DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	 const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
/* computes element DOF by L2 projection on element */
{
    GRID *g = dof->g;
    int i, j, k, n, n0, n1, no, m = dof->dim;
    int npf[NFace];
    QUAD *quad;
    FLOAT x, y, z, *fvals, *A, *B = dofvalues;
    FLOAT w0, *r;
    DOF_TYPE *t;
    const FLOAT *p, *q, *qw, *qp, *b;

    assert(type == VOLUME);
    assert(index == 0);
    assert(userfunc_lambda != NULL || userfunc != NULL);
    assert(funcvalues == NULL);

    t = dof->type;
    if (t == NULL)
	t = dof->hp->info->types[dof->hp->elem_order[e->index]];
    n = t->np_elem;
    if (n == 0)
	return;
    fvals = phgAlloc(Dim * m * sizeof(*fvals));
    A = phgCalloc(n * n, sizeof(*A));

    for (i = 0; i < n * m; i++)
	B[i] = 0.0;

    /* range of basis functions in the current element */
    if (!DofIsHP(dof)) {
	for (i = 0; i < NFace; i++)
	    npf[i] = t->np_face;
	n0 = t->np_face * NFace;
    }
    else {
	n0 = 0;
	for (i = 0; i < NFace; i++) {
	    npf[i] = dof->hp->info->types
			[dof->hp->face_order[e->faces[i]]]->np_face;
	    n0 += npf[i];
	}
    }
    n1 = n0 + n;

    quad = phgQuadGetQuad3D(2 * t->order);
    qp = quad->points;
    qw = quad->weights;
    for (k = 0; k < quad->npoints; k++, qp += 4) {
	if (userfunc_lambda != NULL) {
	    userfunc_lambda(dof, e, -1, qp, fvals);
	}
	else {
	    phgGeomLambda2XYZ(g, e, qp, &x, &y, &z);
	    userfunc(x, y, z, fvals);
	}
	/* evaluate basis functions, since we also need values of the
	 * basis functions on all edges/faces, the range is set to
	 * (0, n1) */
	b = t->BasFuncs(dof, e, 0, n1, qp);
	/* fvals -= contributions from faces. Note that at this stage
	 * dof->data already contains valid data on all faces */
	for (no = 0, p = b; no < NFace; no++) {
	    if (pdofvalues == NULL)
		q = DofFaceData(dof, e->faces[no]);
	    else
		q = pdofvalues[NVert + NEdge + no];
	    for (i = 0; i < npf[no]; i++, p += 3) {
		for (j = 0, r = fvals; j < m; j++) {
		    w0 = *(q++);
		    *(r++) -= w0 * p[0];
		    *(r++) -= w0 * p[1];
		    *(r++) -= w0 * p[2];
		}
	    }
	}
	b += 3 * n0;
	assert(p == b);
	/* coefficients and RHS */
	w0 = *(qw++);
	for (i = 0, p = b; i < n; i++, p += 3) {
	    for (j = i, q = p; j < n; j++, q += 3)
		A[i * n + j] += (p[0] * q[0] + p[1] * q[1] + p[2] * q[2]) * w0;
	    for (j = 0, q = fvals; j < m; j++, q += 3) {
		B[i * m + j] += (p[0] * q[0] + p[1] * q[1] + p[2] * q[2]) * w0;
	    }
	}
    }

    for (i = 1; i < n; i++)
	for (j = 0; j < i; j++)
	    A[i * n + j] = A[j * n + i];

#if 0
printf("---------------------------------------------\n");
for (i = 0; i < n; i++) {
for (j = 0; j < n; j++)
printf(" %012.4le", (double)(Fabs(w0 = A[i * n + j]) < 1e-10 ? 0. : w0));
printf("\n");
}
printf("---------------------------------------------\n");
printf("\n");
#endif

    show_condition(n, A, type, index);

    /* solver the linear system of equations */
    if (!phgSolverDenseSolver(n, m, (void *)A, (void *)B))
	phgError(1, "(%s:%d) singular system (shouldn't happen).\n",
			__FILE__, __LINE__);

    phgFree(A);
    phgFree(fvals);
}

static void
HD_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	 DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	 const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
{
#if 0 * DEBUG
    static BOOLEAN warned = FALSE;
    if (phgRank == 0 && !warned && dof->type->orders == NULL) {
	warned = TRUE;
	phgWarning("the 'orders' member of \"%s\" is undefined!\n",
			dof->type->name);
    }
#endif

    switch (type) {
	case VERTEX:
	case EDGE:
	    /*phgError(1, "(%s:%d) unexpected error.\n", __FILE__, __LINE__);*/
	    break;
	case FACE:
	    HD2_init(dof, e, type, index, userfunc, userfunc_lambda,
			funcvalues, dofvalues, pdofvalues);
	    break;
	case VOLUME:
	    HD3_init(dof, e, type, index, userfunc, userfunc_lambda,
			funcvalues, dofvalues, pdofvalues);
	    break;
    }
}

/*-------------------------- DOF types HDn --------------------------*/

HP_INFO hp_info_HD = {
    DOF_HDn,		/* types */
    "Hierarchic Hdiv",	/* name */
    FE_Hdiv,		/* FE space */
    3,			/* dim */
    1,			/* min_order */
    15,			/* max_order */
    FACE_FLAG | ELEM_FLAG
};

#define NBasFaceHD(p)	(3 * p + NBasFace(p))
#define NBasElemHD(p)	(NEdge * 1 * NBasEdge(p) +	/* edge-based */ \
			 NFace * 2 * NBasFace(p) +	/* face-based */ \
				 3 * NBasElem(p))	/* elem-based */

#define DofTypeHD(order, type, name, orders, grad)			\
    DOF_TYPE type = {DofReserved, name, NULL, orders,			\
			grad, NULL, &hp_info_HD,			\
			phgDofInterpC2FGeneric,				\
			phgDofInterpF2CGeneric,				\
			HD_init, HD_bas, HD_grad, NULL, FE_Hdiv,	\
			FALSE, FALSE, FALSE, -1,			\
			NFace * NBasFaceHD(order) + NBasElemHD(order),	\
			order,		/* order */			\
			0,		/* D_order */			\
			-1,		/* continuity */		\
			3,		/* dim */			\
			0, 0, NBasFaceHD(order), NBasElemHD(order)};

DofTypeHD(1, DOF_HD1_, "HD1", NULL, DOF_DG0)
DofTypeHD(2, DOF_HD2_, "HD2", NULL, DOF_DG1)
DofTypeHD(3, DOF_HD3_, "HD3", NULL, DOF_DG2)
DofTypeHD(4, DOF_HD4_, "HD4", NULL, DOF_DG3)
DofTypeHD(5, DOF_HD5_, "HD5", NULL, DOF_DG4)
DofTypeHD(6, DOF_HD6_, "HD6", NULL, DOF_DG5)
DofTypeHD(7, DOF_HD7_, "HD7", NULL, DOF_DG6)
DofTypeHD(8, DOF_HD8_, "HD8", NULL, DOF_DG7)
DofTypeHD(9, DOF_HD9_, "HD9", NULL, DOF_DG8)
DofTypeHD(10, DOF_HD10_, "HD10", NULL, DOF_DG9)
DofTypeHD(11, DOF_HD11_, "HD11", NULL, DOF_DG10)
DofTypeHD(12, DOF_HD12_, "HD12", NULL, DOF_DG11)
DofTypeHD(13, DOF_HD13_, "HD13", NULL, DOF_DG12)
DofTypeHD(14, DOF_HD14_, "HD14", NULL, DOF_DG13)
DofTypeHD(15, DOF_HD15_, "HD15", NULL, DOF_DG14)

DOF_TYPE *DOF_HDn[] = {
    NULL, DOF_HD1, DOF_HD2, DOF_HD3, DOF_HD4,
    DOF_HD5, DOF_HD6, DOF_HD7, DOF_HD8, DOF_HD9,
    DOF_HD10, DOF_HD11, DOF_HD12, DOF_HD13, DOF_HD14, DOF_HD15,
    NULL
};

/*------------------------------------------------------------------------*/

#ifdef TEST
static void     
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = x * x - y * y + 2. * z * z;
}

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *u;
    const char *fn = "../test/tetra.dat";

    phgInit(&argc, &argv);
    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "(%s:%d) can't read file \"%s\"\n", __FILE__, __LINE__, fn);
    u = phgDofNew(g, &DOF_HB1_, 1, NULL, func_u);
    phgDofFree(&u);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}

#endif	/* TEST */
