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

/* $Id: lagrange1.c,v 1.33 2022/04/09 01:28:34 zlb Exp $
 *
 * This file implements arbitrary order Lagrange elements on uniform
 * lattice nodes.
 *
 * Note: suppose the lattice points are $\{(p_j,q_k,r_l,s_m)\}$,
 * then the basis for the node with barycentric coordinate $(a,b,c,d)$ is:
 $$
 	\prod_{p_j<a}\frac{\lambda_0-p_j}{a-p_j} \cdot
 	\prod_{q_k<b}\frac{\lambda_1-q_k}{b-q_k} \cdot
 	\prod_{r_l<c}\frac{\lambda_2-r_l}{c-r_l} \cdot
 	\prod_{s_m<d}\frac{\lambda_3-s_m}{d-s_m}
 $$
 */

#include "phg.h"
#include <math.h>	/* fabs() */

static FLOAT eps = (16. * FLOAT_EPSILON);

#if defined(GENERATE_STRUCTS)

/* Note: ../utils/lagrange_gen calls "make CFLAGS=-DGENERATE_STRUCTS lagrange1"
 * to compile this section of code for generating the data structures for
 * DOF_Px and DOF_DGx.  */

#include <strings.h>	/* bzero() */

static int
edge_nodes(int p, int (*lambda)[2])
/* returns the number of interior nodes on an edge for an order-'p' basis,
 * the coordinates times 'p' of the nodes are returned in 'lambda' */
{
    int i;

    for (i = p - 1; i > 0; i--, lambda++) {
	(*lambda)[0] = i;
	(*lambda)[1] = p - i;
    }

    return p <= 1 ? 0 : p - 1;
}

static int
face_nodes(int p, int (*lambda)[3])
/* returns the number of interior nodes on a face for an order-'p' basis,
 * the coordinates times 'p' of the nodes are returned in 'lambda' */
{
    int i, j;

    for (i = p - 2; i > 0; i--) {
	for (j = p - i - 1; j > 0; j--, lambda++) {
	    (*lambda)[0] = i;
	    (*lambda)[1] = j;
	    (*lambda)[2] = p - i - j;
	}
    }

    return p <= 2 ? 0 : (p - 1) * (p - 2) / 2;
}

static int
elem_nodes(int p, int (*lambda)[4])
/* returns the number of interior nodes in an element for an order-'p' basis,
 * the coordinates times 'p' of the nodes are returned in 'lambda' */
{
    int i, j, k;

    for (i = p - 3; i > 0; i--) {
	for (j = p - i - 2; j > 0; j--) {
	    for (k = p - i - j - 1; k > 0; k--, lambda++) {
		(*lambda)[0] = i;
		(*lambda)[1] = j;
		(*lambda)[2] = k;
		(*lambda)[3] = p - i - j - k;
	    }
	}
    }

    return p <= 3 ? 0 : (p - 1) * (p - 2) * (p - 3) / 6;
}

int
main(int argc, char *argv[])
{
    int i, j, k, v0, v1, v2, n1, n2, n3, n, p;
    int lambda1[(ORDER_MAX - 1)][2];
    int lambda2[(ORDER_MAX - 1) * (ORDER_MAX - 2) / 2][3];
    int lambda3[(ORDER_MAX - 1) * (ORDER_MAX - 2) * (ORDER_MAX - 3) / 6][4];
    /* lambda is used to store list of all nodes on an element */
    int lambda[(ORDER_MAX + 1) * (ORDER_MAX + 2) * (ORDER_MAX + 3) / 6][4];
    DOF_TYPE *type;
    DOF *dof;
    GRID *g;

    phgInit(&argc, &argv);

    if (phgNProcs > 1) {
	fprintf(stderr, "Error: this program must be run with 1 process.\n");
	phgFinalize();
	return 1;
    }

    if (argc != 2) {
	fprintf(stderr, "Usage: %s order\n", argv[0]);
	phgFinalize();
	return 2;
    }

    p = atoi(argv[1]);
    if (p < 1 || p > ORDER_MAX) {
	fprintf(stderr, "Invalid order %d (1 <= p <= %d).\n", p, 1, ORDER_MAX);
	phgFinalize();
	return 3;
    }

    n1 = edge_nodes(p, lambda1);
    printf("edge nodes: %d\n", n1);
    for (i = 0; i < n1; i++)
	printf("   _F(%2d.)/_F(%d.), _F(%2d.)/_F(%d.),\n",
		lambda1[i][0], p,
		lambda1[i][1], p);

    n2 = face_nodes(p, lambda2);
    printf("face nodes: %d\n", n2);
    for (i = 0; i < n2; i++)
	printf("   _F(%2d.)/_F(%d.), _F(%2d.)/_F(%d.), _F(%2d.)/_F(%d.),\n",
		lambda2[i][0], p,
		lambda2[i][1], p,
		lambda2[i][2], p);

    n3 = elem_nodes(p, lambda3);
    printf("elem nodes: %d\n", n3);
    for (i = 0; i < n3; i++)
	printf("   _F(%2d.)/_F(%d.), _F(%2d.)/_F(%d.), _F(%2d.)/_F(%d.), "
		"_F(%2d.)/_F(%d.)%s\n",
		lambda3[i][0], p,
		lambda3[i][1], p,
		lambda3[i][2], p,
		lambda3[i][3], p,
		i == n3 - 1 ? "" : ",");

    n = (p + 1) * (p + 2) * (p + 3) / 6;

    i = 0;

    /* The 4 vertices */
    lambda[i][0] = p; 
    lambda[i][1] = 0; 
    lambda[i][2] = 0; 
    lambda[i][3] = 0; 
    i++;
    lambda[i][0] = 0; 
    lambda[i][1] = p; 
    lambda[i][2] = 0; 
    lambda[i][3] = 0; 
    i++;
    lambda[i][0] = 0; 
    lambda[i][1] = 0; 
    lambda[i][2] = p; 
    lambda[i][3] = 0; 
    i++;
    lambda[i][0] = 0; 
    lambda[i][1] = 0; 
    lambda[i][2] = 0; 
    lambda[i][3] = p; 
    i++;

    /* The 6 edges */
    for (j = 0; j < NEdge; j++) {
	v0 = GetEdgeVertex(j, 0);
	v1 = GetEdgeVertex(j, 1);
	for (k = 0; k < n1; k++) {
	    bzero(lambda[i], sizeof(lambda[i]));
	    lambda[i][v0] = lambda1[k][0];
	    lambda[i][v1] = lambda1[k][1];
	    i++;
	}
    }

    /* The 4 faces */
    for (j = 0; j < NFace; j++) {
	v0 = GetFaceVertex(j, 0);
	v1 = GetFaceVertex(j, 1);
	v2 = GetFaceVertex(j, 2);
	for (k = 0; k < n2; k++) {
	    lambda[i][v0] = lambda2[k][0];
	    lambda[i][v1] = lambda2[k][1];
	    lambda[i][v2] = lambda2[k][2];
	    lambda[i][j] = 0;
	    i++;
	}
    }

    for (k = 0; k < n3; k++) {
	lambda[i][0] = lambda3[k][0];
	lambda[i][1] = lambda3[k][1];
	lambda[i][2] = lambda3[k][2];
	lambda[i][3] = lambda3[k][3];
	i++;
    }

    assert(i == n);

    printf("DG nodes: %d\n", n);
    for (i = 0; i < n; i++)
	printf("   _F(%2d.)/_F(%d.), _F(%2d.)/_F(%d.), _F(%2d.)/_F(%d.), "
		"_F(%2d.)/_F(%d.)%s\n",
		lambda[i][0], p,
		lambda[i][1], p,
		lambda[i][2], p,
		lambda[i][3], p,
		i == n - 1 ? "" : ",");

    fprintf(stderr, "Testing evaluation of basis functions: ");
    type = DOF_Pn[p];
    if (type == NULL)
	goto end;
    assert(n == type->nbas);
    g = phgNewGrid(-1);
    if (!phgImport(g, "../test/tetra.dat", FALSE))
	phgError(1, "can't read file \"../test/tetra.dat\".\n");
    dof = phgDofNew(g, type, 1, "dof", DofNoAction);
    for (i = 0; i < n; i++) {
	FLOAT lam[] =  {lambda[i][0] / (FLOAT)p, lambda[i][1] / (FLOAT)p,
			lambda[i][2] / (FLOAT)p, lambda[i][3] / (FLOAT)p};
#if 0
	/* Note: this code produces the "TLS definition in libphg.so section
	 * .tbss mismatches non-TLS reference in xxxx.o" error */
	const FLOAT *bas = Pn_bas(dof, g->roots, 0, n, lam);
#else	/* 0 */
	const FLOAT *bas = type->BasFuncs(dof, g->roots, 0, n, lam);
#endif	/* 0 */
	for (j = 0; j < n; j++) {
	    if (Fabs(bas[j] - (i == j ? 1.0 : 0.0)) > eps) {
#if 0
		phgInfo(-1, "lam = (%g,%g,%g,%g)\n",
			lam[0], lam[1], lam[2], lam[3]);
		phgInfo(-1, "lambda = (%d,%d,%d,%d)\n",
			lambda[i][0], lambda[i][1], lambda[i][2], lambda[i][3]);
#endif
		phgError(1, "error: node %d, basis %d, value %g (expect %s)\n",
			 i, j, bas[j], i == j ? "1" : "0");
	    }
	}
    }
    phgDofFree(&dof);
    phgFreeGrid(&g);
    fprintf(stderr, "passed.\n");

end:
    phgFinalize();
    return 0;
}

#else	/* defined(GENERATE_STRUCTS) */

static const FLOAT *
Pn_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    /* Static buffers for returning values of basis functions.
     * Note: one buffer is used for each order to reduce buffer conflicts. */
    static FLOAT **buffers0 = NULL;
    FLOAT **buffers;

    DOF_TYPE *t;
    int order, nbas, bufferlen;
    FLOAT *bas, *p, *p0, a, b, c, onedp;
    FLOAT lam0, lam1, lam2, lam3, L0, L1, L2, L3;
    const FLOAT *pts;
    int i, j, k, r, s, i0, v0, v1, v2, v3, m, n;

    if (buffers0 == NULL) {
#if USE_OMP
#pragma omp critical (Pn_bas)
#endif	/* USE_OMP */
	if (buffers0 == NULL) {
	    buffers0 = phgCalloc((ORDER_MAX + 1) * phgMaxThreads,
				 sizeof(*buffers));
	    FreeAtExitNoCheck_(buffers0, 1);
	}
    }

    buffers = buffers0 + phgThreadId * (ORDER_MAX + 1);

    t = dof->type;
    if (t == NULL) {
	assert(dof->hp != NULL);
	assert(dof->hp->info == HP_DG);
	t = dof->hp->info->types[dof->hp->elem_order[e->index]];
    }
    if (t->base_type != NULL)
	t = t->base_type;
    order = t->order;
    assert(order >= 1 && order <= ORDER_MAX);
    assert(t == DOF_Pn[order]);
    bufferlen = nbas = t->nbas;

    if (no1 <= 0)
	no1 = nbas;

    if (buffers[order] == NULL) {
	buffers[order] = phgAlloc(bufferlen * sizeof(*bas));
	FreeAtExitNoCheck(buffers[order]);
    }
    bas = buffers[order];

    p = bas;
    p0 = p + (no1 - no0);

    onedp = 1. / (FLOAT)order;
    pts = t->points;

    /*------------------------- vertex functions --------------------------*/
    for (i = no0; i < NVert && i < no1; i++) {
	/* Note: phi = prod_{j < order} (lambda - j/order) / (1. - j/order) */
	lam0 = lambda[i];
	b = 1.0;
	for (r = 0; r < order; r++) {
	    c = r * onedp;
	    b *= (lam0 - c) / (1. - c);
	}
	*(p++) = b;
    }

    if (i >= no1)
	return bas;

    if ((no0 -= NVert) < 0)
	no0 = 0;
    pts += 1 * t->np_vert;

    /*------------------------- edge functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NEdge && p < p0; i++) {
	if (t->np_edge == 0)
	    continue;
	n = i0 + t->np_edge;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	v0 = GetEdgeVertex(i, 0);
	v1 = GetEdgeVertex(i, 1);
	lam0 = lambda[v0];
	lam1 = lambda[v1];
	for (j = 0; i0 < n && p < p0; i0++, j++) {
	    if (i0 < no0)
		continue;

	    L0 = pts[2 * j];
	    L1 = pts[2 * j + 1];
	    b = 1.0;

	    a = L0 - eps;
	    for (r = 0; ; r++) {
		c = r * onedp;
		if (c >= a)
		    break;
		b *= (lam0 - c) / (L0 - c);
	    }

	    a = L1 - eps;
	    for (r = 0; ; r++) {
		c = r * onedp;
		if (c >= a)
		    break;
		b *= (lam1 - c) / (L1 - c);
	    }

	    *(p++) = b;
	}
    }

    if (p >= p0)
	return bas;

    if ((no0 -= i0) < 0)
	no0 = 0;
    pts += 2 * t->np_edge;

    /*------------------------- face functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NFace && p < p0; i++) {
	if (t->np_face == 0)
	    continue;
	n = i0 + t->np_face;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	v0 = GetFaceVertex(i, 0);
	v1 = GetFaceVertex(i, 1);
	v2 = GetFaceVertex(i, 2);
	lam0 = lambda[v0];
	lam1 = lambda[v1];
	lam2 = lambda[v2];
	s = 0;
	for (m = 0; m <= order - 3 && p < p0; m++) {
	    if (i0 + m + 1 <= no0) {
		i0 += m + 1;
		s += 3 * (m + 1);
		continue;
	    }
	    for (j = 0; j <= m && p < p0; j++, i0++, s += 3) {
		if (i0 < no0)
		    continue;

		L0 = pts[s];
		L1 = pts[s + 1];
		L2 = pts[s + 2];
		b = 1.0;

		a = L0 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    b *= (lam0 - c) / (L0 - c);
		}

		a = L1 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    b *= (lam1 - c) / (L1 - c);
		}

		a = L2 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    b *= (lam2 - c) / (L2 - c);
		}

		*(p++) = b;
	    }
	}
    }

    if (p >= p0)
	return bas;

    if ((no0 -= i0) < 0)
	no0 = 0;
    pts += 3 * t->np_face;

    /*------------------------- element functions --------------------------*/
    n = t->np_elem;
    assert(no0 < n);
    assert(order >= 4);
    v0 = 0;
    v1 = 1;
    v2 = 2;
    v3 = 3;
    lam0 = lambda[v0];
    lam1 = lambda[v1];
    lam2 = lambda[v2];
    lam3 = lambda[v3];
    i0 = 0;
    s = 0;
    for (m = 0; m <= order - 4 && p < p0; m++) {
	if (i0 + (m + 1) * (m + 2) / 2 <= no0) {
	    i0 += (m + 1) * (m + 2) / 2;
	    s += 4 * (m + 1) * (m + 2) / 2;
	    continue;
	}
	for (k = 0; k <= m && p < p0; k++) {
	    if (i0 + m - k + 1 <= no0) {
		i0 += m - k + 1;
		s += 4 * (m - k + 1);
		continue;
	    }
	    for (j = 0; j <= m - k && p < p0; j++, i0++, s += 4) {
		if (i0 < no0)
		    continue;

		L0 = pts[s];
		L1 = pts[s + 1];
		L2 = pts[s + 2];
		L3 = pts[s + 3];
		b = 1.0;

		a = L0 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    b *= (lam0 - c) / (L0 - c);
		}

		a = L1 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    b *= (lam1 - c) / (L1 - c);
		}

		a = L2 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    b *= (lam2 - c) / (L2 - c);
		}

		a = L3 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    b *= (lam3 - c) / (L3 - c);
		}

		*(p++) = b;
	    }
	}
    }

    if (p != p0)
	phgError(1, "(%s:%d): p=%d, p0=%d\n", __FILE__, __LINE__,
					(int)(p - bas), (int)(p0 - bas));

    return bas;
}

static const FLOAT *
Pn_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT **buffers0 = NULL;
    FLOAT **buffers;

    DOF_TYPE *t;
    int order, nbas, bufferlen, n0, n1, n2, n3;
    FLOAT *bas, *p, *p0, a, b, c, d, onedp;
    FLOAT lam0, lam1, lam2, lam3, L0, L1, L2, L3;
    const FLOAT *pts;
    int i, j, k, r, s, i0, v0, v1, v2, v3, u0, u1, m, n;

    if (buffers0 == NULL) {
#if USE_OMP
#pragma omp critical (Pn_bas)
#endif	/* USE_OMP */
	if (buffers0 == NULL) {
	    buffers0 = phgCalloc((ORDER_MAX + 1) * phgMaxThreads,
				 sizeof(*buffers));
	    FreeAtExitNoCheck_(buffers0, 1);
	}
    }

    buffers = buffers0 + phgThreadId * (ORDER_MAX + 1);

    t = dof->type;
    if (t == NULL) {
	assert(dof->hp != NULL);
	assert(dof->hp->info == HP_DG);
	t = dof->hp->info->types[dof->hp->elem_order[e->index]];
    }
    if (t->base_type != NULL)
	t = t->base_type;
    order = t->order;
    assert(order >= 1 && order <= ORDER_MAX);
    assert(t == DOF_Pn[order]);
    bufferlen = (nbas = t->nbas) * (Dim + 1);

    if (no1 <= 0)
	no1 = nbas;

    if (buffers[order] == NULL) {
	buffers[order] = phgAlloc(bufferlen * sizeof(*bas));
	FreeAtExitNoCheck(buffers[order]);
    }
    bas = buffers[order];

    p = bas;
    p0 = p + (no1 - no0) * (Dim + 1);

    onedp = 1. / (FLOAT)order;
    pts = t->points;

    /*------------------------- vertex functions --------------------------*/
    for (i = no0; i < NVert && i < no1; i++) {
	/* Note: phi = prod_{j < order} (lambda - j/order) / (1. - j/order) */
	v0 = GetFaceVertex(i, 0);
	v1 = GetFaceVertex(i, 1);
	v2 = GetFaceVertex(i, 2);

	/* zero out irrelevant entries */
	p[v0] = p[v1] = p[v2] = 0.;

	n0 = 0;		/* # of 0 factors */
	lam0 = lambda[i];
	b = 1.0;
	for (r = 0; r < order; r++) {
	    d = (lam0 - (c = r * onedp));
	    if (Fabs(d) < eps) {
		d = 1.0;
		n0++;
	    }
	    b *= d / (1. - c);
	}

	/* compute d\lambda[i] */
	if (n0 > 0) {
	    /* one of the factors is zero ==> only one term */
	    p[i] = b;
	    goto cont1;
	}

	for (d = 0., r = 0; r < order; r++) {
	    c = r * onedp;
	    d += b / (lam0 - c);
	}
	p[i] = d;

    cont1:
	p += Dim + 1;
    }

    if (i >= no1)
	return bas;

    if ((no0 -= NVert) < 0)
	no0 = 0;
    pts += 1 * t->np_vert;

    /*------------------------- edge functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NEdge && p < p0; i++) {
	if (t->np_edge == 0)
	    continue;
	n = i0 + t->np_edge;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	v0 = GetEdgeVertex(i, 0);
	v1 = GetEdgeVertex(i, 1);
	lam0 = lambda[v0];
	lam1 = lambda[v1];
	/* The two opposite vertices (note: the opposite edge is 5-i) */
	u0 = GetEdgeVertex(5 - i, 0);
	u1 = GetEdgeVertex(5 - i, 1);
	for (j = 0; i0 < n && p < p0; i0++, j++) {
	    if (i0 < no0)
		continue;

	    /* zero out irrelevant entries */
	    p[u0] = p[u1] = 0.;

	    L0 = pts[2 * j];
	    L1 = pts[2 * j + 1];
	    b = 1.0;

	    n0 = 0;		/* # of zero factor */
	    a = L0 - eps;
	    for (r = 0; ; r++) {
		c = r * onedp;
		if (c >= a)
		    break;
		if (Fabs(d = lam0 - c) < eps) {
		    d = 1.0;
		    n0++;
		}
		b *= d / (L0 - c);
	    }

	    n1 = 0;
	    a = L1 - eps;
	    for (r = 0; ; r++) {
		c = r * onedp;
		if (c >= a)
		    break;
		if (Fabs(d = lam1 - c) < eps) {
		    if (n0 > 0) {
			/* two or more zero factors */
			p[v0] = p[v1] = 0.;
			goto cont2;
		    }
		    d = 1.0;
		    n1++;
		}
		b *= d / (L1 - c);
	    }

	    assert(n0 + n1 <= 1);
	    if (n0 > 0) {
		p[v0] = b;
		p[v1] = 0.;
		goto cont2;
	    }
	    else if (n1 > 0) {
		p[v0] = 0.;
		p[v1] = b;
		goto cont2;
	    }

	    /* compute d\lambda[v0] */
	    a = L0 - eps;
	    for (d = 0., r = 0; ; r++) {
		c = r * onedp;
		if (c >= a)
		    break;
		d += b / (lam0 - c);
	    }
	    p[v0] = d;

	    /* compute d\lambda[v1] */
	    a = L1 - eps;
	    for (d = 0., r = 0; ; r++) {
		c = r * onedp;
		if (c >= a)
		    break;
		d += b / (lam1 - c);
	    }
	    p[v1] = d;

	cont2:
	    p += Dim + 1;
	}
    }

    if (p >= p0)
	return bas;

    if ((no0 -= i0) < 0)
	no0 = 0;
    pts += 2 * t->np_edge;

    /*------------------------- face functions --------------------------*/
    i0 = 0;
    for (i = 0; i < NFace && p < p0; i++) {
	if (t->np_face == 0)
	    continue;
	n = i0 + t->np_face;
	if (n <= no0) {
	    i0 = n;
	    continue;
	}
	v0 = GetFaceVertex(i, 0);
	v1 = GetFaceVertex(i, 1);
	v2 = GetFaceVertex(i, 2);
	lam0 = lambda[v0];
	lam1 = lambda[v1];
	lam2 = lambda[v2];
	s = 0;
	for (m = 0; m <= order - 3 && p < p0; m++) {
	    if (i0 + m + 1 <= no0) {
		i0 += m + 1;
		s += 3 * (m + 1);
		continue;
	    }
	    for (j = 0; j <= m && p < p0; j++, i0++, s += 3) {
		if (i0 < no0)
		    continue;

		/* zero out irrelevant entry */
		p[i] = 0.;

		L0 = pts[s];
		L1 = pts[s + 1];
		L2 = pts[s + 2];
		b = 1.0;

		n0 = 0;		/* # of zero factor */
		a = L0 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    if (Fabs(d = lam0 - c) < eps) {
			d = 1.0;
			n0++;
		    }
		    b *= d / (L0 - c);
		}

		n1 = 0;		/* # of zero factor */
		a = L1 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    if (Fabs(d = lam1 - c) < eps) {
			if (n0 > 0) {
			    /* two or more zero factors */
			    p[v0] = 0.;
			    p[v1] = 0.;
			    p[v2] = 0.;
			    goto cont3;
			}
			d = 1.0;
			n1++;
		    }
		    b *= d / (L1 - c);
		}

		n2 = 0;		/* # of zero factor */
		a = L2 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    if (Fabs(d = lam2 - c) < eps) {
			if (n0 + n1 > 0) {
			    /* two or more zero factors */
			    p[v0] = 0.;
			    p[v1] = 0.;
			    p[v2] = 0.;
			    goto cont3;
			}
			d = 1.0;
			n2++;
		    }
		    b *= d / (L2 - c);
		}

		assert(n0 + n1 + n2 <= 1);
		if (n0 > 0) {
		    p[v0] = b;
		    p[v1] = 0.;
		    p[v2] = 0.;
		    goto cont3;
		}
		else if (n1 > 0) {
		    p[v0] = 0.;
		    p[v1] = b;
		    p[v2] = 0.;
		    goto cont3;
		}
		else if (n2 > 0) {
		    p[v0] = 0.;
		    p[v1] = 0.;
		    p[v2] = b;
		    goto cont3;
		}

		/* compute d\lambda[v0] */
		a = L0 - eps;
		for (d = 0., r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    d += b / (lam0 - c);
		}
		p[v0] = d;

		/* compute d\lambda[v1] */
		a = L1 - eps;
		for (d = 0., r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    d += b / (lam1 - c);
		}
		p[v1] = d;

		/* compute d\lambda[v2] */
		a = L2 - eps;
		for (d = 0., r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    d += b / (lam2 - c);
		}
		p[v2] = d;

	cont3:
		p += Dim + 1;
	    }
	}
    }

    if (p >= p0)
	return bas;

    if ((no0 -= i0) < 0)
	no0 = 0;
    pts += 3 * t->np_face;

    /*------------------------- element functions --------------------------*/
    n = t->np_elem;
    assert(no0 < n);
    assert(order >= 4);
    v0 = 0;
    v1 = 1;
    v2 = 2;
    v3 = 3;
    lam0 = lambda[v0];
    lam1 = lambda[v1];
    lam2 = lambda[v2];
    lam3 = lambda[v3];
    i0 = 0;
    s = 0;
    for (m = 0; m <= order - 4 && p < p0; m++) {
	if (i0 + (m + 1) * (m + 2) / 2 <= no0) {
	    i0 += (m + 1) * (m + 2) / 2;
	    s += 4 * (m + 1) * (m + 2) / 2;
	    continue;
	}
	for (k = 0; k <= m && p < p0; k++) {
	    if (i0 + m - k + 1 <= no0) {
		i0 += m - k + 1;
		s += 4 * (m - k + 1);
		continue;
	    }
	    for (j = 0; j <= m - k && p < p0; j++, i0++, s += 4) {
		if (i0 < no0)
		    continue;

		L0 = pts[s];
		L1 = pts[s + 1];
		L2 = pts[s + 2];
		L3 = pts[s + 3];
		b = 1.0;

		n0 = 0;		/* # of zero factor */
		a = L0 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    if (Fabs(d = lam0 - c) < eps) {
			d = 1.0;
			n0++;
		    }
		    b *= d / (L0 - c);
		}

		n1 = 0;		/* # of zero factor */
		a = L1 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    if (Fabs(d = lam1 - c) < eps) {
			if (n0 > 0) {
			    /* two or more zero factors */
			    p[v0] = 0.;
			    p[v1] = 0.;
			    p[v2] = 0.;
			    p[v3] = 0.;
			    goto cont4;
			}
			d = 1.0;
			n1++;
		    }
		    b *= d / (L1 - c);
		}

		n2 = 0;		/* # of zero factor */
		a = L2 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    if (Fabs(d = lam2 - c) < eps) {
			if (n0 + n1 > 0) {
			    /* two or more zero factors */
			    p[v0] = 0.;
			    p[v1] = 0.;
			    p[v2] = 0.;
			    p[v3] = 0.;
			    goto cont4;
			}
			d = 1.0;
			n2++;
		    }
		    b *= d / (L2 - c);
		}

		n3 = 0;		/* # of zero factor */
		a = L3 - eps;
		for (r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    if (Fabs(d = lam3 - c) < eps) {
			if (n0 + n1 + n2 > 0) {
			    /* two or more zero factors */
			    p[v0] = 0.;
			    p[v1] = 0.;
			    p[v2] = 0.;
			    p[v3] = 0.;
			    goto cont4;
			}
			d = 1.0;
			n3++;
		    }
		    b *= d / (L3 - c);
		}

		assert(n0 + n1 + n2 + n3 <= 1);
		if (n0 > 0) {
		    p[v0] = b;
		    p[v1] = 0.;
		    p[v2] = 0.;
		    p[v3] = 0.;
		    goto cont4;
		}
		else if (n1 > 0) {
		    p[v0] = 0.;
		    p[v1] = b;
		    p[v2] = 0.;
		    p[v3] = 0.;
		    goto cont4;
		}
		else if (n2 > 0) {
		    p[v0] = 0.;
		    p[v1] = 0.;
		    p[v2] = b;
		    p[v3] = 0.;
		    goto cont4;
		}
		else if (n3 > 0) {
		    p[v0] = 0.;
		    p[v1] = 0.;
		    p[v2] = 0.;
		    p[v3] = b;
		    goto cont4;
		}

		/* compute d\lambda[v0] */
		a = L0 - eps;
		for (d = 0., r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    d += b / (lam0 - c);
		}
		p[v0] = d;

		/* compute d\lambda[v1] */
		a = L1 - eps;
		for (d = 0., r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    d += b / (lam1 - c);
		}
		p[v1] = d;

		/* compute d\lambda[v2] */
		a = L2 - eps;
		for (d = 0., r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    d += b / (lam2 - c);
		}
		p[v2] = d;

		/* compute d\lambda[v3] */
		a = L3 - eps;
		for (d = 0., r = 0; ; r++) {
		    c = r * onedp;
		    if (c >= a)
			break;
		    d += b / (lam3 - c);
		}
		p[v3] = d;

	cont4:
		p += Dim + 1;
	    }
	}
    }

    if (p != p0)
	phgError(1, "(%s:%d): p=%d, p0=%d\n", __FILE__, __LINE__,
					(int)(p - bas), (int)(p0 - bas));

    return bas;
}

static void
Pn_fmap(const DOF_TYPE *type, const ELEMENT *e, int face_no, int *M)
/* maps ordering of local face data to ordering of data in DOF */
{
    int p, i, j, k, l, v0, v1, v2, u0, u1, u2, tmp[3];

    if (type->np_face <= 1) {
	for (i = 0; i < type->np_face; i++)
	    M[i] = i;
	return;
    }

    p = type->order;
    u0 = MapVertex(e, GetFaceVertex(face_no, 0));
    u1 = MapVertex(e, GetFaceVertex(face_no, 1));
    u2 = MapVertex(e, GetFaceVertex(face_no, 2));
    v0 = v1 = v2 = 0;
    (u0 > u1) ? v0++ : v1++;
    (u0 > u2) ? v0++ : v2++;
    (u1 > u2) ? v1++ : v2++;

    for (k = 0, i = p - 2; i > 0; i--) {
	for (j = p - i - 1; j > 0; j--, k++) {
	    tmp[v0] = i;
	    tmp[v1] = j;
	    tmp[v2] = p - i - j;
	    u0 = tmp[0];
	    u1 = tmp[1];
#if 0
	    l = 0;
	    for (u2 = p - 2; u2 > u0; u2--)
		l += p - u2 - 1;
	    l += (p - u0 - u1) - 1;
#else
	    /* Note: l = \sum_{i=u0+1}^{p-2}(p-i-1) + p-u0-u1 - 1
	     *	       = \sum_{i=1}^{p-2-u0}i + p-u0-u1 - 1
	     *	       = (p-u0-2)*(p-u0-1)/2 + p-u0-u1 - 1 */
	    l = p - u0;
	    l = ((l - 1) * (l - 2)) / 2 + (l - u1) - 1;
#endif
	    M[k] = l;
        }
    }

    return;
}

static void
Pn_elem_map(const DOF_TYPE *type, const ELEMENT *e, int *M)
/* maps ordering of local element data to ordering of data in DOF.
 * Note: see elem_nodes() above for the ordering of the nodes. */
{
    int p, i, j, k, m, v0, v1, v2, v3, u0, u1, u2, tmp[4];

    if (type->np_elem <= 1 || ElementIsInOrder(e)) {
	for (i = 0; i < type->np_elem; i++)
	    M[i] = i;
	return;
    }

    v0 = MapVertex(e, 0);
    v1 = MapVertex(e, 1);
    v2 = MapVertex(e, 2);
    v3 = MapVertex(e, 3);

    p = type->order;
    m = 0;
    for (i = p - 3; i > 0; i--) {
	for (j = p - i - 2; j > 0; j--) {
	    for (k = p - i - j - 1; k > 0; k--, m++) {
		tmp[v0] = i;
		tmp[v1] = j;
		tmp[v2] = k;
		tmp[v3] = p - i - j - k;
		u0 = tmp[0];
		u1 = tmp[1];
		u2 = tmp[2];
		/* Note: count # points up to the coord (u0,u1,u2)
		 * Part 1:
		 *	for (i = u0+1; i <= p-3; i++)
		 *	  for (j = 1; j <= p-i-2; j++)
		 *	    for (k = 1; k <= p-i-j-1; k++, count++);
		 *	==> s1:=Sum[Sum[p-i-j-1, {j,1,p-i-2}],{i,u0+1,p-3}]
		 * Part 2:
		 *	for (j = u1+1; j <= p-u0-2; j++)
		 *	  for (k = 1; k <= p-u0-j-1; k++, count++);
		 *	==> s2:=Sum[p-u0-j-1,{j,u1+1,p-u0-2}]
		 * Part 3:
		 *	for (k = u2+1; k <= p-u0-u1-1; k++, count++);
		 *	==> s3:=p-u0-u1-u2-1
		 * Total:
		 *	==> CForm[Simplify[s1+s2+s3]]
		 */
#define Power2(p) ((p)*(p))
#define Power3(p) ((p)*Power2(p))
		M[m] = (-6 + Power3(p) - 8*u0 - 3*Power2(u0) - Power3(u0) - 
			3*Power2(p)*(1 + u0) + p*(8 + 6*u0 + 3*Power2(u0) -
			6*u1) + 3*u1 + 6*u0*u1 + 3*Power2(u1) - 6*u2)/6;
	    }
        }
    }

    return;
}

static void
Pn_map(const DOF_TYPE *type, const ELEMENT *e, int no, int *M)
{
    if (no < 0)
	Pn_elem_map(type, e, M);
    else
	Pn_fmap(type, e, no, M);
}

/* Note: insert the file "/tmp/DOF_Pn.inc" generated by the script
 *	 "../utils/lagrange_gen" below */

/*======= The following structs are generated by 'utils/lagrange_gen' ======*/

static FLOAT P5_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
    _F(4.)/_F(5.),  _F(1.)/_F(5.),
    _F(3.)/_F(5.),  _F(2.)/_F(5.),
    _F(2.)/_F(5.),  _F(3.)/_F(5.),
    _F(1.)/_F(5.),  _F(4.)/_F(5.),
   /* face */
    _F(3.)/_F(5.),  _F(1.)/_F(5.),  _F(1.)/_F(5.),
    _F(2.)/_F(5.),  _F(2.)/_F(5.),  _F(1.)/_F(5.),
    _F(2.)/_F(5.),  _F(1.)/_F(5.),  _F(2.)/_F(5.),
    _F(1.)/_F(5.),  _F(3.)/_F(5.),  _F(1.)/_F(5.),
    _F(1.)/_F(5.),  _F(2.)/_F(5.),  _F(2.)/_F(5.),
    _F(1.)/_F(5.),  _F(1.)/_F(5.),  _F(3.)/_F(5.),
   /* elem */
    _F(2.)/_F(5.),  _F(1.)/_F(5.),  _F(1.)/_F(5.),  _F(1.)/_F(5.),
    _F(1.)/_F(5.),  _F(2.)/_F(5.),  _F(1.)/_F(5.),  _F(1.)/_F(5.),
    _F(1.)/_F(5.),  _F(1.)/_F(5.),  _F(2.)/_F(5.),  _F(1.)/_F(5.),
    _F(1.)/_F(5.),  _F(1.)/_F(5.),  _F(1.)/_F(5.),  _F(2.)/_F(5.)
};
DOF_TYPE DOF_P5_ = {
    DofReserved, "P5", P5_points, NULL, DOF_DG4, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    56, 5, 0, 0, 1, 1, 4, 6, 4
};

static FLOAT P6_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
    _F(5.)/_F(6.),  _F(1.)/_F(6.),
    _F(4.)/_F(6.),  _F(2.)/_F(6.),
    _F(3.)/_F(6.),  _F(3.)/_F(6.),
    _F(2.)/_F(6.),  _F(4.)/_F(6.),
    _F(1.)/_F(6.),  _F(5.)/_F(6.),
   /* face */
    _F(4.)/_F(6.),  _F(1.)/_F(6.),  _F(1.)/_F(6.),
    _F(3.)/_F(6.),  _F(2.)/_F(6.),  _F(1.)/_F(6.),
    _F(3.)/_F(6.),  _F(1.)/_F(6.),  _F(2.)/_F(6.),
    _F(2.)/_F(6.),  _F(3.)/_F(6.),  _F(1.)/_F(6.),
    _F(2.)/_F(6.),  _F(2.)/_F(6.),  _F(2.)/_F(6.),
    _F(2.)/_F(6.),  _F(1.)/_F(6.),  _F(3.)/_F(6.),
    _F(1.)/_F(6.),  _F(4.)/_F(6.),  _F(1.)/_F(6.),
    _F(1.)/_F(6.),  _F(3.)/_F(6.),  _F(2.)/_F(6.),
    _F(1.)/_F(6.),  _F(2.)/_F(6.),  _F(3.)/_F(6.),
    _F(1.)/_F(6.),  _F(1.)/_F(6.),  _F(4.)/_F(6.),
   /* elem */
    _F(3.)/_F(6.),  _F(1.)/_F(6.),  _F(1.)/_F(6.),  _F(1.)/_F(6.),
    _F(2.)/_F(6.),  _F(2.)/_F(6.),  _F(1.)/_F(6.),  _F(1.)/_F(6.),
    _F(2.)/_F(6.),  _F(1.)/_F(6.),  _F(2.)/_F(6.),  _F(1.)/_F(6.),
    _F(2.)/_F(6.),  _F(1.)/_F(6.),  _F(1.)/_F(6.),  _F(2.)/_F(6.),
    _F(1.)/_F(6.),  _F(3.)/_F(6.),  _F(1.)/_F(6.),  _F(1.)/_F(6.),
    _F(1.)/_F(6.),  _F(2.)/_F(6.),  _F(2.)/_F(6.),  _F(1.)/_F(6.),
    _F(1.)/_F(6.),  _F(2.)/_F(6.),  _F(1.)/_F(6.),  _F(2.)/_F(6.),
    _F(1.)/_F(6.),  _F(1.)/_F(6.),  _F(3.)/_F(6.),  _F(1.)/_F(6.),
    _F(1.)/_F(6.),  _F(1.)/_F(6.),  _F(2.)/_F(6.),  _F(2.)/_F(6.),
    _F(1.)/_F(6.),  _F(1.)/_F(6.),  _F(1.)/_F(6.),  _F(3.)/_F(6.)
};
DOF_TYPE DOF_P6_ = {
    DofReserved, "P6", P6_points, NULL, DOF_DG5, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    84, 6, 0, 0, 1, 1, 5, 10, 10
};

static FLOAT P7_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
    _F(6.)/_F(7.),  _F(1.)/_F(7.),
    _F(5.)/_F(7.),  _F(2.)/_F(7.),
    _F(4.)/_F(7.),  _F(3.)/_F(7.),
    _F(3.)/_F(7.),  _F(4.)/_F(7.),
    _F(2.)/_F(7.),  _F(5.)/_F(7.),
    _F(1.)/_F(7.),  _F(6.)/_F(7.),
   /* face */
    _F(5.)/_F(7.),  _F(1.)/_F(7.),  _F(1.)/_F(7.),
    _F(4.)/_F(7.),  _F(2.)/_F(7.),  _F(1.)/_F(7.),
    _F(4.)/_F(7.),  _F(1.)/_F(7.),  _F(2.)/_F(7.),
    _F(3.)/_F(7.),  _F(3.)/_F(7.),  _F(1.)/_F(7.),
    _F(3.)/_F(7.),  _F(2.)/_F(7.),  _F(2.)/_F(7.),
    _F(3.)/_F(7.),  _F(1.)/_F(7.),  _F(3.)/_F(7.),
    _F(2.)/_F(7.),  _F(4.)/_F(7.),  _F(1.)/_F(7.),
    _F(2.)/_F(7.),  _F(3.)/_F(7.),  _F(2.)/_F(7.),
    _F(2.)/_F(7.),  _F(2.)/_F(7.),  _F(3.)/_F(7.),
    _F(2.)/_F(7.),  _F(1.)/_F(7.),  _F(4.)/_F(7.),
    _F(1.)/_F(7.),  _F(5.)/_F(7.),  _F(1.)/_F(7.),
    _F(1.)/_F(7.),  _F(4.)/_F(7.),  _F(2.)/_F(7.),
    _F(1.)/_F(7.),  _F(3.)/_F(7.),  _F(3.)/_F(7.),
    _F(1.)/_F(7.),  _F(2.)/_F(7.),  _F(4.)/_F(7.),
    _F(1.)/_F(7.),  _F(1.)/_F(7.),  _F(5.)/_F(7.),
   /* elem */
    _F(4.)/_F(7.),  _F(1.)/_F(7.),  _F(1.)/_F(7.),  _F(1.)/_F(7.),
    _F(3.)/_F(7.),  _F(2.)/_F(7.),  _F(1.)/_F(7.),  _F(1.)/_F(7.),
    _F(3.)/_F(7.),  _F(1.)/_F(7.),  _F(2.)/_F(7.),  _F(1.)/_F(7.),
    _F(3.)/_F(7.),  _F(1.)/_F(7.),  _F(1.)/_F(7.),  _F(2.)/_F(7.),
    _F(2.)/_F(7.),  _F(3.)/_F(7.),  _F(1.)/_F(7.),  _F(1.)/_F(7.),
    _F(2.)/_F(7.),  _F(2.)/_F(7.),  _F(2.)/_F(7.),  _F(1.)/_F(7.),
    _F(2.)/_F(7.),  _F(2.)/_F(7.),  _F(1.)/_F(7.),  _F(2.)/_F(7.),
    _F(2.)/_F(7.),  _F(1.)/_F(7.),  _F(3.)/_F(7.),  _F(1.)/_F(7.),
    _F(2.)/_F(7.),  _F(1.)/_F(7.),  _F(2.)/_F(7.),  _F(2.)/_F(7.),
    _F(2.)/_F(7.),  _F(1.)/_F(7.),  _F(1.)/_F(7.),  _F(3.)/_F(7.),
    _F(1.)/_F(7.),  _F(4.)/_F(7.),  _F(1.)/_F(7.),  _F(1.)/_F(7.),
    _F(1.)/_F(7.),  _F(3.)/_F(7.),  _F(2.)/_F(7.),  _F(1.)/_F(7.),
    _F(1.)/_F(7.),  _F(3.)/_F(7.),  _F(1.)/_F(7.),  _F(2.)/_F(7.),
    _F(1.)/_F(7.),  _F(2.)/_F(7.),  _F(3.)/_F(7.),  _F(1.)/_F(7.),
    _F(1.)/_F(7.),  _F(2.)/_F(7.),  _F(2.)/_F(7.),  _F(2.)/_F(7.),
    _F(1.)/_F(7.),  _F(2.)/_F(7.),  _F(1.)/_F(7.),  _F(3.)/_F(7.),
    _F(1.)/_F(7.),  _F(1.)/_F(7.),  _F(4.)/_F(7.),  _F(1.)/_F(7.),
    _F(1.)/_F(7.),  _F(1.)/_F(7.),  _F(3.)/_F(7.),  _F(2.)/_F(7.),
    _F(1.)/_F(7.),  _F(1.)/_F(7.),  _F(2.)/_F(7.),  _F(3.)/_F(7.),
    _F(1.)/_F(7.),  _F(1.)/_F(7.),  _F(1.)/_F(7.),  _F(4.)/_F(7.)
};
DOF_TYPE DOF_P7_ = {
    DofReserved, "P7", P7_points, NULL, DOF_DG6, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    120, 7, 0, 0, 1, 1, 6, 15, 20
};

static FLOAT P8_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
    _F(7.)/_F(8.),  _F(1.)/_F(8.),
    _F(6.)/_F(8.),  _F(2.)/_F(8.),
    _F(5.)/_F(8.),  _F(3.)/_F(8.),
    _F(4.)/_F(8.),  _F(4.)/_F(8.),
    _F(3.)/_F(8.),  _F(5.)/_F(8.),
    _F(2.)/_F(8.),  _F(6.)/_F(8.),
    _F(1.)/_F(8.),  _F(7.)/_F(8.),
   /* face */
    _F(6.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),
    _F(5.)/_F(8.),  _F(2.)/_F(8.),  _F(1.)/_F(8.),
    _F(5.)/_F(8.),  _F(1.)/_F(8.),  _F(2.)/_F(8.),
    _F(4.)/_F(8.),  _F(3.)/_F(8.),  _F(1.)/_F(8.),
    _F(4.)/_F(8.),  _F(2.)/_F(8.),  _F(2.)/_F(8.),
    _F(4.)/_F(8.),  _F(1.)/_F(8.),  _F(3.)/_F(8.),
    _F(3.)/_F(8.),  _F(4.)/_F(8.),  _F(1.)/_F(8.),
    _F(3.)/_F(8.),  _F(3.)/_F(8.),  _F(2.)/_F(8.),
    _F(3.)/_F(8.),  _F(2.)/_F(8.),  _F(3.)/_F(8.),
    _F(3.)/_F(8.),  _F(1.)/_F(8.),  _F(4.)/_F(8.),
    _F(2.)/_F(8.),  _F(5.)/_F(8.),  _F(1.)/_F(8.),
    _F(2.)/_F(8.),  _F(4.)/_F(8.),  _F(2.)/_F(8.),
    _F(2.)/_F(8.),  _F(3.)/_F(8.),  _F(3.)/_F(8.),
    _F(2.)/_F(8.),  _F(2.)/_F(8.),  _F(4.)/_F(8.),
    _F(2.)/_F(8.),  _F(1.)/_F(8.),  _F(5.)/_F(8.),
    _F(1.)/_F(8.),  _F(6.)/_F(8.),  _F(1.)/_F(8.),
    _F(1.)/_F(8.),  _F(5.)/_F(8.),  _F(2.)/_F(8.),
    _F(1.)/_F(8.),  _F(4.)/_F(8.),  _F(3.)/_F(8.),
    _F(1.)/_F(8.),  _F(3.)/_F(8.),  _F(4.)/_F(8.),
    _F(1.)/_F(8.),  _F(2.)/_F(8.),  _F(5.)/_F(8.),
    _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(6.)/_F(8.),
   /* elem */
    _F(5.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),
    _F(4.)/_F(8.),  _F(2.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),
    _F(4.)/_F(8.),  _F(1.)/_F(8.),  _F(2.)/_F(8.),  _F(1.)/_F(8.),
    _F(4.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(2.)/_F(8.),
    _F(3.)/_F(8.),  _F(3.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),
    _F(3.)/_F(8.),  _F(2.)/_F(8.),  _F(2.)/_F(8.),  _F(1.)/_F(8.),
    _F(3.)/_F(8.),  _F(2.)/_F(8.),  _F(1.)/_F(8.),  _F(2.)/_F(8.),
    _F(3.)/_F(8.),  _F(1.)/_F(8.),  _F(3.)/_F(8.),  _F(1.)/_F(8.),
    _F(3.)/_F(8.),  _F(1.)/_F(8.),  _F(2.)/_F(8.),  _F(2.)/_F(8.),
    _F(3.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(3.)/_F(8.),
    _F(2.)/_F(8.),  _F(4.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),
    _F(2.)/_F(8.),  _F(3.)/_F(8.),  _F(2.)/_F(8.),  _F(1.)/_F(8.),
    _F(2.)/_F(8.),  _F(3.)/_F(8.),  _F(1.)/_F(8.),  _F(2.)/_F(8.),
    _F(2.)/_F(8.),  _F(2.)/_F(8.),  _F(3.)/_F(8.),  _F(1.)/_F(8.),
    _F(2.)/_F(8.),  _F(2.)/_F(8.),  _F(2.)/_F(8.),  _F(2.)/_F(8.),
    _F(2.)/_F(8.),  _F(2.)/_F(8.),  _F(1.)/_F(8.),  _F(3.)/_F(8.),
    _F(2.)/_F(8.),  _F(1.)/_F(8.),  _F(4.)/_F(8.),  _F(1.)/_F(8.),
    _F(2.)/_F(8.),  _F(1.)/_F(8.),  _F(3.)/_F(8.),  _F(2.)/_F(8.),
    _F(2.)/_F(8.),  _F(1.)/_F(8.),  _F(2.)/_F(8.),  _F(3.)/_F(8.),
    _F(2.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(4.)/_F(8.),
    _F(1.)/_F(8.),  _F(5.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),
    _F(1.)/_F(8.),  _F(4.)/_F(8.),  _F(2.)/_F(8.),  _F(1.)/_F(8.),
    _F(1.)/_F(8.),  _F(4.)/_F(8.),  _F(1.)/_F(8.),  _F(2.)/_F(8.),
    _F(1.)/_F(8.),  _F(3.)/_F(8.),  _F(3.)/_F(8.),  _F(1.)/_F(8.),
    _F(1.)/_F(8.),  _F(3.)/_F(8.),  _F(2.)/_F(8.),  _F(2.)/_F(8.),
    _F(1.)/_F(8.),  _F(3.)/_F(8.),  _F(1.)/_F(8.),  _F(3.)/_F(8.),
    _F(1.)/_F(8.),  _F(2.)/_F(8.),  _F(4.)/_F(8.),  _F(1.)/_F(8.),
    _F(1.)/_F(8.),  _F(2.)/_F(8.),  _F(3.)/_F(8.),  _F(2.)/_F(8.),
    _F(1.)/_F(8.),  _F(2.)/_F(8.),  _F(2.)/_F(8.),  _F(3.)/_F(8.),
    _F(1.)/_F(8.),  _F(2.)/_F(8.),  _F(1.)/_F(8.),  _F(4.)/_F(8.),
    _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(5.)/_F(8.),  _F(1.)/_F(8.),
    _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(4.)/_F(8.),  _F(2.)/_F(8.),
    _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(3.)/_F(8.),  _F(3.)/_F(8.),
    _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(2.)/_F(8.),  _F(4.)/_F(8.),
    _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(1.)/_F(8.),  _F(5.)/_F(8.)
};
DOF_TYPE DOF_P8_ = {
    DofReserved, "P8", P8_points, NULL, DOF_DG7, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    165, 8, 0, 0, 1, 1, 7, 21, 35
};

static FLOAT P9_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
    _F(8.)/_F(9.),  _F(1.)/_F(9.),
    _F(7.)/_F(9.),  _F(2.)/_F(9.),
    _F(6.)/_F(9.),  _F(3.)/_F(9.),
    _F(5.)/_F(9.),  _F(4.)/_F(9.),
    _F(4.)/_F(9.),  _F(5.)/_F(9.),
    _F(3.)/_F(9.),  _F(6.)/_F(9.),
    _F(2.)/_F(9.),  _F(7.)/_F(9.),
    _F(1.)/_F(9.),  _F(8.)/_F(9.),
   /* face */
    _F(7.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),
    _F(6.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),
    _F(6.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),
    _F(5.)/_F(9.),  _F(3.)/_F(9.),  _F(1.)/_F(9.),
    _F(5.)/_F(9.),  _F(2.)/_F(9.),  _F(2.)/_F(9.),
    _F(5.)/_F(9.),  _F(1.)/_F(9.),  _F(3.)/_F(9.),
    _F(4.)/_F(9.),  _F(4.)/_F(9.),  _F(1.)/_F(9.),
    _F(4.)/_F(9.),  _F(3.)/_F(9.),  _F(2.)/_F(9.),
    _F(4.)/_F(9.),  _F(2.)/_F(9.),  _F(3.)/_F(9.),
    _F(4.)/_F(9.),  _F(1.)/_F(9.),  _F(4.)/_F(9.),
    _F(3.)/_F(9.),  _F(5.)/_F(9.),  _F(1.)/_F(9.),
    _F(3.)/_F(9.),  _F(4.)/_F(9.),  _F(2.)/_F(9.),
    _F(3.)/_F(9.),  _F(3.)/_F(9.),  _F(3.)/_F(9.),
    _F(3.)/_F(9.),  _F(2.)/_F(9.),  _F(4.)/_F(9.),
    _F(3.)/_F(9.),  _F(1.)/_F(9.),  _F(5.)/_F(9.),
    _F(2.)/_F(9.),  _F(6.)/_F(9.),  _F(1.)/_F(9.),
    _F(2.)/_F(9.),  _F(5.)/_F(9.),  _F(2.)/_F(9.),
    _F(2.)/_F(9.),  _F(4.)/_F(9.),  _F(3.)/_F(9.),
    _F(2.)/_F(9.),  _F(3.)/_F(9.),  _F(4.)/_F(9.),
    _F(2.)/_F(9.),  _F(2.)/_F(9.),  _F(5.)/_F(9.),
    _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(6.)/_F(9.),
    _F(1.)/_F(9.),  _F(7.)/_F(9.),  _F(1.)/_F(9.),
    _F(1.)/_F(9.),  _F(6.)/_F(9.),  _F(2.)/_F(9.),
    _F(1.)/_F(9.),  _F(5.)/_F(9.),  _F(3.)/_F(9.),
    _F(1.)/_F(9.),  _F(4.)/_F(9.),  _F(4.)/_F(9.),
    _F(1.)/_F(9.),  _F(3.)/_F(9.),  _F(5.)/_F(9.),
    _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(6.)/_F(9.),
    _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(7.)/_F(9.),
   /* elem */
    _F(6.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),
    _F(5.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),
    _F(5.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),
    _F(5.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),
    _F(4.)/_F(9.),  _F(3.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),
    _F(4.)/_F(9.),  _F(2.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),
    _F(4.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),
    _F(4.)/_F(9.),  _F(1.)/_F(9.),  _F(3.)/_F(9.),  _F(1.)/_F(9.),
    _F(4.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(2.)/_F(9.),
    _F(4.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(3.)/_F(9.),
    _F(3.)/_F(9.),  _F(4.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),
    _F(3.)/_F(9.),  _F(3.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),
    _F(3.)/_F(9.),  _F(3.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),
    _F(3.)/_F(9.),  _F(2.)/_F(9.),  _F(3.)/_F(9.),  _F(1.)/_F(9.),
    _F(3.)/_F(9.),  _F(2.)/_F(9.),  _F(2.)/_F(9.),  _F(2.)/_F(9.),
    _F(3.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(3.)/_F(9.),
    _F(3.)/_F(9.),  _F(1.)/_F(9.),  _F(4.)/_F(9.),  _F(1.)/_F(9.),
    _F(3.)/_F(9.),  _F(1.)/_F(9.),  _F(3.)/_F(9.),  _F(2.)/_F(9.),
    _F(3.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(3.)/_F(9.),
    _F(3.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(4.)/_F(9.),
    _F(2.)/_F(9.),  _F(5.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),
    _F(2.)/_F(9.),  _F(4.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),
    _F(2.)/_F(9.),  _F(4.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),
    _F(2.)/_F(9.),  _F(3.)/_F(9.),  _F(3.)/_F(9.),  _F(1.)/_F(9.),
    _F(2.)/_F(9.),  _F(3.)/_F(9.),  _F(2.)/_F(9.),  _F(2.)/_F(9.),
    _F(2.)/_F(9.),  _F(3.)/_F(9.),  _F(1.)/_F(9.),  _F(3.)/_F(9.),
    _F(2.)/_F(9.),  _F(2.)/_F(9.),  _F(4.)/_F(9.),  _F(1.)/_F(9.),
    _F(2.)/_F(9.),  _F(2.)/_F(9.),  _F(3.)/_F(9.),  _F(2.)/_F(9.),
    _F(2.)/_F(9.),  _F(2.)/_F(9.),  _F(2.)/_F(9.),  _F(3.)/_F(9.),
    _F(2.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(4.)/_F(9.),
    _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(5.)/_F(9.),  _F(1.)/_F(9.),
    _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(4.)/_F(9.),  _F(2.)/_F(9.),
    _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(3.)/_F(9.),  _F(3.)/_F(9.),
    _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(4.)/_F(9.),
    _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(5.)/_F(9.),
    _F(1.)/_F(9.),  _F(6.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),
    _F(1.)/_F(9.),  _F(5.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),
    _F(1.)/_F(9.),  _F(5.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),
    _F(1.)/_F(9.),  _F(4.)/_F(9.),  _F(3.)/_F(9.),  _F(1.)/_F(9.),
    _F(1.)/_F(9.),  _F(4.)/_F(9.),  _F(2.)/_F(9.),  _F(2.)/_F(9.),
    _F(1.)/_F(9.),  _F(4.)/_F(9.),  _F(1.)/_F(9.),  _F(3.)/_F(9.),
    _F(1.)/_F(9.),  _F(3.)/_F(9.),  _F(4.)/_F(9.),  _F(1.)/_F(9.),
    _F(1.)/_F(9.),  _F(3.)/_F(9.),  _F(3.)/_F(9.),  _F(2.)/_F(9.),
    _F(1.)/_F(9.),  _F(3.)/_F(9.),  _F(2.)/_F(9.),  _F(3.)/_F(9.),
    _F(1.)/_F(9.),  _F(3.)/_F(9.),  _F(1.)/_F(9.),  _F(4.)/_F(9.),
    _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(5.)/_F(9.),  _F(1.)/_F(9.),
    _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(4.)/_F(9.),  _F(2.)/_F(9.),
    _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(3.)/_F(9.),  _F(3.)/_F(9.),
    _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(2.)/_F(9.),  _F(4.)/_F(9.),
    _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(1.)/_F(9.),  _F(5.)/_F(9.),
    _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(6.)/_F(9.),  _F(1.)/_F(9.),
    _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(5.)/_F(9.),  _F(2.)/_F(9.),
    _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(4.)/_F(9.),  _F(3.)/_F(9.),
    _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(3.)/_F(9.),  _F(4.)/_F(9.),
    _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(2.)/_F(9.),  _F(5.)/_F(9.),
    _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(1.)/_F(9.),  _F(6.)/_F(9.)
};
DOF_TYPE DOF_P9_ = {
    DofReserved, "P9", P9_points, NULL, DOF_DG8, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    220, 9, 0, 0, 1, 1, 8, 28, 56
};

static FLOAT P10_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
    _F(9.)/_F(10.),  _F(1.)/_F(10.),
    _F(8.)/_F(10.),  _F(2.)/_F(10.),
    _F(7.)/_F(10.),  _F(3.)/_F(10.),
    _F(6.)/_F(10.),  _F(4.)/_F(10.),
    _F(5.)/_F(10.),  _F(5.)/_F(10.),
    _F(4.)/_F(10.),  _F(6.)/_F(10.),
    _F(3.)/_F(10.),  _F(7.)/_F(10.),
    _F(2.)/_F(10.),  _F(8.)/_F(10.),
    _F(1.)/_F(10.),  _F(9.)/_F(10.),
   /* face */
    _F(8.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),
    _F(7.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),
    _F(7.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),
    _F(6.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),
    _F(6.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),
    _F(6.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),
    _F(5.)/_F(10.),  _F(4.)/_F(10.),  _F(1.)/_F(10.),
    _F(5.)/_F(10.),  _F(3.)/_F(10.),  _F(2.)/_F(10.),
    _F(5.)/_F(10.),  _F(2.)/_F(10.),  _F(3.)/_F(10.),
    _F(5.)/_F(10.),  _F(1.)/_F(10.),  _F(4.)/_F(10.),
    _F(4.)/_F(10.),  _F(5.)/_F(10.),  _F(1.)/_F(10.),
    _F(4.)/_F(10.),  _F(4.)/_F(10.),  _F(2.)/_F(10.),
    _F(4.)/_F(10.),  _F(3.)/_F(10.),  _F(3.)/_F(10.),
    _F(4.)/_F(10.),  _F(2.)/_F(10.),  _F(4.)/_F(10.),
    _F(4.)/_F(10.),  _F(1.)/_F(10.),  _F(5.)/_F(10.),
    _F(3.)/_F(10.),  _F(6.)/_F(10.),  _F(1.)/_F(10.),
    _F(3.)/_F(10.),  _F(5.)/_F(10.),  _F(2.)/_F(10.),
    _F(3.)/_F(10.),  _F(4.)/_F(10.),  _F(3.)/_F(10.),
    _F(3.)/_F(10.),  _F(3.)/_F(10.),  _F(4.)/_F(10.),
    _F(3.)/_F(10.),  _F(2.)/_F(10.),  _F(5.)/_F(10.),
    _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(6.)/_F(10.),
    _F(2.)/_F(10.),  _F(7.)/_F(10.),  _F(1.)/_F(10.),
    _F(2.)/_F(10.),  _F(6.)/_F(10.),  _F(2.)/_F(10.),
    _F(2.)/_F(10.),  _F(5.)/_F(10.),  _F(3.)/_F(10.),
    _F(2.)/_F(10.),  _F(4.)/_F(10.),  _F(4.)/_F(10.),
    _F(2.)/_F(10.),  _F(3.)/_F(10.),  _F(5.)/_F(10.),
    _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(6.)/_F(10.),
    _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(7.)/_F(10.),
    _F(1.)/_F(10.),  _F(8.)/_F(10.),  _F(1.)/_F(10.),
    _F(1.)/_F(10.),  _F(7.)/_F(10.),  _F(2.)/_F(10.),
    _F(1.)/_F(10.),  _F(6.)/_F(10.),  _F(3.)/_F(10.),
    _F(1.)/_F(10.),  _F(5.)/_F(10.),  _F(4.)/_F(10.),
    _F(1.)/_F(10.),  _F(4.)/_F(10.),  _F(5.)/_F(10.),
    _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(6.)/_F(10.),
    _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(7.)/_F(10.),
    _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(8.)/_F(10.),
   /* elem */
    _F(7.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),
    _F(6.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),
    _F(6.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),
    _F(6.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),
    _F(5.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),
    _F(5.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),
    _F(5.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),
    _F(5.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),
    _F(5.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),
    _F(5.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),
    _F(4.)/_F(10.),  _F(4.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),
    _F(4.)/_F(10.),  _F(3.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),
    _F(4.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),
    _F(4.)/_F(10.),  _F(2.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),
    _F(4.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),
    _F(4.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),
    _F(4.)/_F(10.),  _F(1.)/_F(10.),  _F(4.)/_F(10.),  _F(1.)/_F(10.),
    _F(4.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(2.)/_F(10.),
    _F(4.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(3.)/_F(10.),
    _F(4.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(4.)/_F(10.),
    _F(3.)/_F(10.),  _F(5.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),
    _F(3.)/_F(10.),  _F(4.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),
    _F(3.)/_F(10.),  _F(4.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),
    _F(3.)/_F(10.),  _F(3.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),
    _F(3.)/_F(10.),  _F(3.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),
    _F(3.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),
    _F(3.)/_F(10.),  _F(2.)/_F(10.),  _F(4.)/_F(10.),  _F(1.)/_F(10.),
    _F(3.)/_F(10.),  _F(2.)/_F(10.),  _F(3.)/_F(10.),  _F(2.)/_F(10.),
    _F(3.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(3.)/_F(10.),
    _F(3.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(4.)/_F(10.),
    _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(5.)/_F(10.),  _F(1.)/_F(10.),
    _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(4.)/_F(10.),  _F(2.)/_F(10.),
    _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(3.)/_F(10.),
    _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(4.)/_F(10.),
    _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(5.)/_F(10.),
    _F(2.)/_F(10.),  _F(6.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),
    _F(2.)/_F(10.),  _F(5.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),
    _F(2.)/_F(10.),  _F(5.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),
    _F(2.)/_F(10.),  _F(4.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),
    _F(2.)/_F(10.),  _F(4.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),
    _F(2.)/_F(10.),  _F(4.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),
    _F(2.)/_F(10.),  _F(3.)/_F(10.),  _F(4.)/_F(10.),  _F(1.)/_F(10.),
    _F(2.)/_F(10.),  _F(3.)/_F(10.),  _F(3.)/_F(10.),  _F(2.)/_F(10.),
    _F(2.)/_F(10.),  _F(3.)/_F(10.),  _F(2.)/_F(10.),  _F(3.)/_F(10.),
    _F(2.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(4.)/_F(10.),
    _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(5.)/_F(10.),  _F(1.)/_F(10.),
    _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(4.)/_F(10.),  _F(2.)/_F(10.),
    _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(3.)/_F(10.),  _F(3.)/_F(10.),
    _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(4.)/_F(10.),
    _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(5.)/_F(10.),
    _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(6.)/_F(10.),  _F(1.)/_F(10.),
    _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(5.)/_F(10.),  _F(2.)/_F(10.),
    _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(4.)/_F(10.),  _F(3.)/_F(10.),
    _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(4.)/_F(10.),
    _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(5.)/_F(10.),
    _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(6.)/_F(10.),
    _F(1.)/_F(10.),  _F(7.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),
    _F(1.)/_F(10.),  _F(6.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),
    _F(1.)/_F(10.),  _F(6.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),
    _F(1.)/_F(10.),  _F(5.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),
    _F(1.)/_F(10.),  _F(5.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),
    _F(1.)/_F(10.),  _F(5.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),
    _F(1.)/_F(10.),  _F(4.)/_F(10.),  _F(4.)/_F(10.),  _F(1.)/_F(10.),
    _F(1.)/_F(10.),  _F(4.)/_F(10.),  _F(3.)/_F(10.),  _F(2.)/_F(10.),
    _F(1.)/_F(10.),  _F(4.)/_F(10.),  _F(2.)/_F(10.),  _F(3.)/_F(10.),
    _F(1.)/_F(10.),  _F(4.)/_F(10.),  _F(1.)/_F(10.),  _F(4.)/_F(10.),
    _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(5.)/_F(10.),  _F(1.)/_F(10.),
    _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(4.)/_F(10.),  _F(2.)/_F(10.),
    _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(3.)/_F(10.),  _F(3.)/_F(10.),
    _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(2.)/_F(10.),  _F(4.)/_F(10.),
    _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(1.)/_F(10.),  _F(5.)/_F(10.),
    _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(6.)/_F(10.),  _F(1.)/_F(10.),
    _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(5.)/_F(10.),  _F(2.)/_F(10.),
    _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(4.)/_F(10.),  _F(3.)/_F(10.),
    _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(3.)/_F(10.),  _F(4.)/_F(10.),
    _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(2.)/_F(10.),  _F(5.)/_F(10.),
    _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(1.)/_F(10.),  _F(6.)/_F(10.),
    _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(7.)/_F(10.),  _F(1.)/_F(10.),
    _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(6.)/_F(10.),  _F(2.)/_F(10.),
    _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(5.)/_F(10.),  _F(3.)/_F(10.),
    _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(4.)/_F(10.),  _F(4.)/_F(10.),
    _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(3.)/_F(10.),  _F(5.)/_F(10.),
    _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(2.)/_F(10.),  _F(6.)/_F(10.),
    _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(1.)/_F(10.),  _F(7.)/_F(10.)
};
DOF_TYPE DOF_P10_ = {
    DofReserved, "P10", P10_points, NULL, DOF_DG9, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    286, 10, 0, 0, 1, 1, 9, 36, 84
};

static FLOAT P11_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
   _F(10.)/_F(11.),  _F(1.)/_F(11.),
    _F(9.)/_F(11.),  _F(2.)/_F(11.),
    _F(8.)/_F(11.),  _F(3.)/_F(11.),
    _F(7.)/_F(11.),  _F(4.)/_F(11.),
    _F(6.)/_F(11.),  _F(5.)/_F(11.),
    _F(5.)/_F(11.),  _F(6.)/_F(11.),
    _F(4.)/_F(11.),  _F(7.)/_F(11.),
    _F(3.)/_F(11.),  _F(8.)/_F(11.),
    _F(2.)/_F(11.),  _F(9.)/_F(11.),
    _F(1.)/_F(11.), _F(10.)/_F(11.),
   /* face */
    _F(9.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),
    _F(8.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),
    _F(8.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),
    _F(7.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),
    _F(7.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),
    _F(7.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),
    _F(6.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),
    _F(6.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),
    _F(6.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),
    _F(6.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),
    _F(5.)/_F(11.),  _F(5.)/_F(11.),  _F(1.)/_F(11.),
    _F(5.)/_F(11.),  _F(4.)/_F(11.),  _F(2.)/_F(11.),
    _F(5.)/_F(11.),  _F(3.)/_F(11.),  _F(3.)/_F(11.),
    _F(5.)/_F(11.),  _F(2.)/_F(11.),  _F(4.)/_F(11.),
    _F(5.)/_F(11.),  _F(1.)/_F(11.),  _F(5.)/_F(11.),
    _F(4.)/_F(11.),  _F(6.)/_F(11.),  _F(1.)/_F(11.),
    _F(4.)/_F(11.),  _F(5.)/_F(11.),  _F(2.)/_F(11.),
    _F(4.)/_F(11.),  _F(4.)/_F(11.),  _F(3.)/_F(11.),
    _F(4.)/_F(11.),  _F(3.)/_F(11.),  _F(4.)/_F(11.),
    _F(4.)/_F(11.),  _F(2.)/_F(11.),  _F(5.)/_F(11.),
    _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(6.)/_F(11.),
    _F(3.)/_F(11.),  _F(7.)/_F(11.),  _F(1.)/_F(11.),
    _F(3.)/_F(11.),  _F(6.)/_F(11.),  _F(2.)/_F(11.),
    _F(3.)/_F(11.),  _F(5.)/_F(11.),  _F(3.)/_F(11.),
    _F(3.)/_F(11.),  _F(4.)/_F(11.),  _F(4.)/_F(11.),
    _F(3.)/_F(11.),  _F(3.)/_F(11.),  _F(5.)/_F(11.),
    _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(6.)/_F(11.),
    _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(7.)/_F(11.),
    _F(2.)/_F(11.),  _F(8.)/_F(11.),  _F(1.)/_F(11.),
    _F(2.)/_F(11.),  _F(7.)/_F(11.),  _F(2.)/_F(11.),
    _F(2.)/_F(11.),  _F(6.)/_F(11.),  _F(3.)/_F(11.),
    _F(2.)/_F(11.),  _F(5.)/_F(11.),  _F(4.)/_F(11.),
    _F(2.)/_F(11.),  _F(4.)/_F(11.),  _F(5.)/_F(11.),
    _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(6.)/_F(11.),
    _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(7.)/_F(11.),
    _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(8.)/_F(11.),
    _F(1.)/_F(11.),  _F(9.)/_F(11.),  _F(1.)/_F(11.),
    _F(1.)/_F(11.),  _F(8.)/_F(11.),  _F(2.)/_F(11.),
    _F(1.)/_F(11.),  _F(7.)/_F(11.),  _F(3.)/_F(11.),
    _F(1.)/_F(11.),  _F(6.)/_F(11.),  _F(4.)/_F(11.),
    _F(1.)/_F(11.),  _F(5.)/_F(11.),  _F(5.)/_F(11.),
    _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(6.)/_F(11.),
    _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(7.)/_F(11.),
    _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(8.)/_F(11.),
    _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(9.)/_F(11.),
   /* elem */
    _F(8.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),
    _F(7.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),
    _F(7.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),
    _F(7.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),
    _F(6.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),
    _F(6.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),
    _F(6.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),
    _F(6.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),
    _F(6.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),
    _F(6.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),
    _F(5.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),
    _F(5.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),
    _F(5.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),
    _F(5.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),
    _F(5.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),
    _F(5.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),
    _F(5.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),
    _F(5.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),
    _F(5.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),
    _F(5.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),
    _F(4.)/_F(11.),  _F(5.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),
    _F(4.)/_F(11.),  _F(4.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),
    _F(4.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),
    _F(4.)/_F(11.),  _F(3.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),
    _F(4.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),
    _F(4.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),
    _F(4.)/_F(11.),  _F(2.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),
    _F(4.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),
    _F(4.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),
    _F(4.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),
    _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(5.)/_F(11.),  _F(1.)/_F(11.),
    _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(2.)/_F(11.),
    _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(3.)/_F(11.),
    _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(4.)/_F(11.),
    _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(5.)/_F(11.),
    _F(3.)/_F(11.),  _F(6.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),
    _F(3.)/_F(11.),  _F(5.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),
    _F(3.)/_F(11.),  _F(5.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),
    _F(3.)/_F(11.),  _F(4.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),
    _F(3.)/_F(11.),  _F(4.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),
    _F(3.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),
    _F(3.)/_F(11.),  _F(3.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),
    _F(3.)/_F(11.),  _F(3.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),
    _F(3.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),
    _F(3.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),
    _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(5.)/_F(11.),  _F(1.)/_F(11.),
    _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(4.)/_F(11.),  _F(2.)/_F(11.),
    _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(3.)/_F(11.),
    _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(4.)/_F(11.),
    _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(5.)/_F(11.),
    _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(6.)/_F(11.),  _F(1.)/_F(11.),
    _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(5.)/_F(11.),  _F(2.)/_F(11.),
    _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(3.)/_F(11.),
    _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(4.)/_F(11.),
    _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(5.)/_F(11.),
    _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(6.)/_F(11.),
    _F(2.)/_F(11.),  _F(7.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),
    _F(2.)/_F(11.),  _F(6.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),
    _F(2.)/_F(11.),  _F(6.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),
    _F(2.)/_F(11.),  _F(5.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),
    _F(2.)/_F(11.),  _F(5.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),
    _F(2.)/_F(11.),  _F(5.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),
    _F(2.)/_F(11.),  _F(4.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),
    _F(2.)/_F(11.),  _F(4.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),
    _F(2.)/_F(11.),  _F(4.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),
    _F(2.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),
    _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(5.)/_F(11.),  _F(1.)/_F(11.),
    _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(4.)/_F(11.),  _F(2.)/_F(11.),
    _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(3.)/_F(11.),  _F(3.)/_F(11.),
    _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(4.)/_F(11.),
    _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(5.)/_F(11.),
    _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(6.)/_F(11.),  _F(1.)/_F(11.),
    _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(5.)/_F(11.),  _F(2.)/_F(11.),
    _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(4.)/_F(11.),  _F(3.)/_F(11.),
    _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(4.)/_F(11.),
    _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(5.)/_F(11.),
    _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(6.)/_F(11.),
    _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(7.)/_F(11.),  _F(1.)/_F(11.),
    _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(6.)/_F(11.),  _F(2.)/_F(11.),
    _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(5.)/_F(11.),  _F(3.)/_F(11.),
    _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(4.)/_F(11.),
    _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(5.)/_F(11.),
    _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(6.)/_F(11.),
    _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(7.)/_F(11.),
    _F(1.)/_F(11.),  _F(8.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),
    _F(1.)/_F(11.),  _F(7.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),
    _F(1.)/_F(11.),  _F(7.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),
    _F(1.)/_F(11.),  _F(6.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),
    _F(1.)/_F(11.),  _F(6.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),
    _F(1.)/_F(11.),  _F(6.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),
    _F(1.)/_F(11.),  _F(5.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),
    _F(1.)/_F(11.),  _F(5.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),
    _F(1.)/_F(11.),  _F(5.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),
    _F(1.)/_F(11.),  _F(5.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),
    _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(5.)/_F(11.),  _F(1.)/_F(11.),
    _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(4.)/_F(11.),  _F(2.)/_F(11.),
    _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(3.)/_F(11.),  _F(3.)/_F(11.),
    _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(2.)/_F(11.),  _F(4.)/_F(11.),
    _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(1.)/_F(11.),  _F(5.)/_F(11.),
    _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(6.)/_F(11.),  _F(1.)/_F(11.),
    _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(5.)/_F(11.),  _F(2.)/_F(11.),
    _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(4.)/_F(11.),  _F(3.)/_F(11.),
    _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(3.)/_F(11.),  _F(4.)/_F(11.),
    _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(2.)/_F(11.),  _F(5.)/_F(11.),
    _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(1.)/_F(11.),  _F(6.)/_F(11.),
    _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(7.)/_F(11.),  _F(1.)/_F(11.),
    _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(6.)/_F(11.),  _F(2.)/_F(11.),
    _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(5.)/_F(11.),  _F(3.)/_F(11.),
    _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(4.)/_F(11.),  _F(4.)/_F(11.),
    _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(3.)/_F(11.),  _F(5.)/_F(11.),
    _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(2.)/_F(11.),  _F(6.)/_F(11.),
    _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(1.)/_F(11.),  _F(7.)/_F(11.),
    _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(8.)/_F(11.),  _F(1.)/_F(11.),
    _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(7.)/_F(11.),  _F(2.)/_F(11.),
    _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(6.)/_F(11.),  _F(3.)/_F(11.),
    _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(5.)/_F(11.),  _F(4.)/_F(11.),
    _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(4.)/_F(11.),  _F(5.)/_F(11.),
    _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(3.)/_F(11.),  _F(6.)/_F(11.),
    _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(2.)/_F(11.),  _F(7.)/_F(11.),
    _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(1.)/_F(11.),  _F(8.)/_F(11.)
};
DOF_TYPE DOF_P11_ = {
    DofReserved, "P11", P11_points, NULL, DOF_DG10, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    364, 11, 0, 0, 1, 1, 10, 45, 120
};

static FLOAT P12_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
   _F(11.)/_F(12.),  _F(1.)/_F(12.),
   _F(10.)/_F(12.),  _F(2.)/_F(12.),
    _F(9.)/_F(12.),  _F(3.)/_F(12.),
    _F(8.)/_F(12.),  _F(4.)/_F(12.),
    _F(7.)/_F(12.),  _F(5.)/_F(12.),
    _F(6.)/_F(12.),  _F(6.)/_F(12.),
    _F(5.)/_F(12.),  _F(7.)/_F(12.),
    _F(4.)/_F(12.),  _F(8.)/_F(12.),
    _F(3.)/_F(12.),  _F(9.)/_F(12.),
    _F(2.)/_F(12.), _F(10.)/_F(12.),
    _F(1.)/_F(12.), _F(11.)/_F(12.),
   /* face */
   _F(10.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),
    _F(9.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),
    _F(9.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),
    _F(8.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),
    _F(8.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),
    _F(8.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),
    _F(7.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),
    _F(7.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),
    _F(7.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),
    _F(7.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),
    _F(6.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),
    _F(6.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),
    _F(6.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),
    _F(6.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),
    _F(6.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),
    _F(5.)/_F(12.),  _F(6.)/_F(12.),  _F(1.)/_F(12.),
    _F(5.)/_F(12.),  _F(5.)/_F(12.),  _F(2.)/_F(12.),
    _F(5.)/_F(12.),  _F(4.)/_F(12.),  _F(3.)/_F(12.),
    _F(5.)/_F(12.),  _F(3.)/_F(12.),  _F(4.)/_F(12.),
    _F(5.)/_F(12.),  _F(2.)/_F(12.),  _F(5.)/_F(12.),
    _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(6.)/_F(12.),
    _F(4.)/_F(12.),  _F(7.)/_F(12.),  _F(1.)/_F(12.),
    _F(4.)/_F(12.),  _F(6.)/_F(12.),  _F(2.)/_F(12.),
    _F(4.)/_F(12.),  _F(5.)/_F(12.),  _F(3.)/_F(12.),
    _F(4.)/_F(12.),  _F(4.)/_F(12.),  _F(4.)/_F(12.),
    _F(4.)/_F(12.),  _F(3.)/_F(12.),  _F(5.)/_F(12.),
    _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(6.)/_F(12.),
    _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(7.)/_F(12.),
    _F(3.)/_F(12.),  _F(8.)/_F(12.),  _F(1.)/_F(12.),
    _F(3.)/_F(12.),  _F(7.)/_F(12.),  _F(2.)/_F(12.),
    _F(3.)/_F(12.),  _F(6.)/_F(12.),  _F(3.)/_F(12.),
    _F(3.)/_F(12.),  _F(5.)/_F(12.),  _F(4.)/_F(12.),
    _F(3.)/_F(12.),  _F(4.)/_F(12.),  _F(5.)/_F(12.),
    _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(6.)/_F(12.),
    _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(7.)/_F(12.),
    _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(8.)/_F(12.),
    _F(2.)/_F(12.),  _F(9.)/_F(12.),  _F(1.)/_F(12.),
    _F(2.)/_F(12.),  _F(8.)/_F(12.),  _F(2.)/_F(12.),
    _F(2.)/_F(12.),  _F(7.)/_F(12.),  _F(3.)/_F(12.),
    _F(2.)/_F(12.),  _F(6.)/_F(12.),  _F(4.)/_F(12.),
    _F(2.)/_F(12.),  _F(5.)/_F(12.),  _F(5.)/_F(12.),
    _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(6.)/_F(12.),
    _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(7.)/_F(12.),
    _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(8.)/_F(12.),
    _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(9.)/_F(12.),
    _F(1.)/_F(12.), _F(10.)/_F(12.),  _F(1.)/_F(12.),
    _F(1.)/_F(12.),  _F(9.)/_F(12.),  _F(2.)/_F(12.),
    _F(1.)/_F(12.),  _F(8.)/_F(12.),  _F(3.)/_F(12.),
    _F(1.)/_F(12.),  _F(7.)/_F(12.),  _F(4.)/_F(12.),
    _F(1.)/_F(12.),  _F(6.)/_F(12.),  _F(5.)/_F(12.),
    _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(6.)/_F(12.),
    _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(7.)/_F(12.),
    _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(8.)/_F(12.),
    _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(9.)/_F(12.),
    _F(1.)/_F(12.),  _F(1.)/_F(12.), _F(10.)/_F(12.),
   /* elem */
    _F(9.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),
    _F(8.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),
    _F(8.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),
    _F(8.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),
    _F(7.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),
    _F(7.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),
    _F(7.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),
    _F(7.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),
    _F(7.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),
    _F(7.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),
    _F(6.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),
    _F(6.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),
    _F(6.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),
    _F(6.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),
    _F(6.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),
    _F(6.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),
    _F(6.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),
    _F(6.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),
    _F(6.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),
    _F(6.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),
    _F(5.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),
    _F(5.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),
    _F(5.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),
    _F(5.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),
    _F(5.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),
    _F(5.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),
    _F(5.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),
    _F(5.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),
    _F(5.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),
    _F(5.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),
    _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),
    _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),
    _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),
    _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),
    _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),
    _F(4.)/_F(12.),  _F(6.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),
    _F(4.)/_F(12.),  _F(5.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),
    _F(4.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),
    _F(4.)/_F(12.),  _F(4.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),
    _F(4.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),
    _F(4.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),
    _F(4.)/_F(12.),  _F(3.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),
    _F(4.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),
    _F(4.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),
    _F(4.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),
    _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),
    _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),
    _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),
    _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),
    _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),
    _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(6.)/_F(12.),  _F(1.)/_F(12.),
    _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(2.)/_F(12.),
    _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(3.)/_F(12.),
    _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(4.)/_F(12.),
    _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(5.)/_F(12.),
    _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(6.)/_F(12.),
    _F(3.)/_F(12.),  _F(7.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),
    _F(3.)/_F(12.),  _F(6.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),
    _F(3.)/_F(12.),  _F(6.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),
    _F(3.)/_F(12.),  _F(5.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),
    _F(3.)/_F(12.),  _F(5.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),
    _F(3.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),
    _F(3.)/_F(12.),  _F(4.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),
    _F(3.)/_F(12.),  _F(4.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),
    _F(3.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),
    _F(3.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),
    _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),
    _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),
    _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),
    _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),
    _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),
    _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(6.)/_F(12.),  _F(1.)/_F(12.),
    _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(5.)/_F(12.),  _F(2.)/_F(12.),
    _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(3.)/_F(12.),
    _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(4.)/_F(12.),
    _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(5.)/_F(12.),
    _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(6.)/_F(12.),
    _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(7.)/_F(12.),  _F(1.)/_F(12.),
    _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(6.)/_F(12.),  _F(2.)/_F(12.),
    _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(3.)/_F(12.),
    _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(4.)/_F(12.),
    _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(5.)/_F(12.),
    _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(6.)/_F(12.),
    _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(7.)/_F(12.),
    _F(2.)/_F(12.),  _F(8.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),
    _F(2.)/_F(12.),  _F(7.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),
    _F(2.)/_F(12.),  _F(7.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),
    _F(2.)/_F(12.),  _F(6.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),
    _F(2.)/_F(12.),  _F(6.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),
    _F(2.)/_F(12.),  _F(6.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),
    _F(2.)/_F(12.),  _F(5.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),
    _F(2.)/_F(12.),  _F(5.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),
    _F(2.)/_F(12.),  _F(5.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),
    _F(2.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),
    _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),
    _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),
    _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),
    _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),
    _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),
    _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(6.)/_F(12.),  _F(1.)/_F(12.),
    _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(5.)/_F(12.),  _F(2.)/_F(12.),
    _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(4.)/_F(12.),  _F(3.)/_F(12.),
    _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(4.)/_F(12.),
    _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(5.)/_F(12.),
    _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(6.)/_F(12.),
    _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(7.)/_F(12.),  _F(1.)/_F(12.),
    _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(6.)/_F(12.),  _F(2.)/_F(12.),
    _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(5.)/_F(12.),  _F(3.)/_F(12.),
    _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(4.)/_F(12.),
    _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(5.)/_F(12.),
    _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(6.)/_F(12.),
    _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(7.)/_F(12.),
    _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(8.)/_F(12.),  _F(1.)/_F(12.),
    _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(7.)/_F(12.),  _F(2.)/_F(12.),
    _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(6.)/_F(12.),  _F(3.)/_F(12.),
    _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(4.)/_F(12.),
    _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(5.)/_F(12.),
    _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(6.)/_F(12.),
    _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(7.)/_F(12.),
    _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(8.)/_F(12.),
    _F(1.)/_F(12.),  _F(9.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),
    _F(1.)/_F(12.),  _F(8.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),
    _F(1.)/_F(12.),  _F(8.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),
    _F(1.)/_F(12.),  _F(7.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),
    _F(1.)/_F(12.),  _F(7.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),
    _F(1.)/_F(12.),  _F(7.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),
    _F(1.)/_F(12.),  _F(6.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),
    _F(1.)/_F(12.),  _F(6.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),
    _F(1.)/_F(12.),  _F(6.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),
    _F(1.)/_F(12.),  _F(6.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),
    _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),
    _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),
    _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),
    _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),
    _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),
    _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(6.)/_F(12.),  _F(1.)/_F(12.),
    _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(5.)/_F(12.),  _F(2.)/_F(12.),
    _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(4.)/_F(12.),  _F(3.)/_F(12.),
    _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(3.)/_F(12.),  _F(4.)/_F(12.),
    _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(2.)/_F(12.),  _F(5.)/_F(12.),
    _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(1.)/_F(12.),  _F(6.)/_F(12.),
    _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(7.)/_F(12.),  _F(1.)/_F(12.),
    _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(6.)/_F(12.),  _F(2.)/_F(12.),
    _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(5.)/_F(12.),  _F(3.)/_F(12.),
    _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(4.)/_F(12.),  _F(4.)/_F(12.),
    _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(3.)/_F(12.),  _F(5.)/_F(12.),
    _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(2.)/_F(12.),  _F(6.)/_F(12.),
    _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(1.)/_F(12.),  _F(7.)/_F(12.),
    _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(8.)/_F(12.),  _F(1.)/_F(12.),
    _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(7.)/_F(12.),  _F(2.)/_F(12.),
    _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(6.)/_F(12.),  _F(3.)/_F(12.),
    _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(5.)/_F(12.),  _F(4.)/_F(12.),
    _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(4.)/_F(12.),  _F(5.)/_F(12.),
    _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(3.)/_F(12.),  _F(6.)/_F(12.),
    _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(2.)/_F(12.),  _F(7.)/_F(12.),
    _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(1.)/_F(12.),  _F(8.)/_F(12.),
    _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(9.)/_F(12.),  _F(1.)/_F(12.),
    _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(8.)/_F(12.),  _F(2.)/_F(12.),
    _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(7.)/_F(12.),  _F(3.)/_F(12.),
    _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(6.)/_F(12.),  _F(4.)/_F(12.),
    _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(5.)/_F(12.),  _F(5.)/_F(12.),
    _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(4.)/_F(12.),  _F(6.)/_F(12.),
    _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(3.)/_F(12.),  _F(7.)/_F(12.),
    _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(2.)/_F(12.),  _F(8.)/_F(12.),
    _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(1.)/_F(12.),  _F(9.)/_F(12.)
};
DOF_TYPE DOF_P12_ = {
    DofReserved, "P12", P12_points, NULL, DOF_DG11, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    455, 12, 0, 0, 1, 1, 11, 55, 165
};

static FLOAT P13_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
   _F(12.)/_F(13.),  _F(1.)/_F(13.),
   _F(11.)/_F(13.),  _F(2.)/_F(13.),
   _F(10.)/_F(13.),  _F(3.)/_F(13.),
    _F(9.)/_F(13.),  _F(4.)/_F(13.),
    _F(8.)/_F(13.),  _F(5.)/_F(13.),
    _F(7.)/_F(13.),  _F(6.)/_F(13.),
    _F(6.)/_F(13.),  _F(7.)/_F(13.),
    _F(5.)/_F(13.),  _F(8.)/_F(13.),
    _F(4.)/_F(13.),  _F(9.)/_F(13.),
    _F(3.)/_F(13.), _F(10.)/_F(13.),
    _F(2.)/_F(13.), _F(11.)/_F(13.),
    _F(1.)/_F(13.), _F(12.)/_F(13.),
   /* face */
   _F(11.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
   _F(10.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),
   _F(10.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),
    _F(9.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),
    _F(9.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),
    _F(9.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),
    _F(8.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),
    _F(8.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),
    _F(8.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),
    _F(8.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),
    _F(7.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),
    _F(7.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),
    _F(7.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),
    _F(7.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),
    _F(7.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),
    _F(6.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),
    _F(6.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),
    _F(6.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),
    _F(6.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),
    _F(6.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),
    _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),
    _F(5.)/_F(13.),  _F(7.)/_F(13.),  _F(1.)/_F(13.),
    _F(5.)/_F(13.),  _F(6.)/_F(13.),  _F(2.)/_F(13.),
    _F(5.)/_F(13.),  _F(5.)/_F(13.),  _F(3.)/_F(13.),
    _F(5.)/_F(13.),  _F(4.)/_F(13.),  _F(4.)/_F(13.),
    _F(5.)/_F(13.),  _F(3.)/_F(13.),  _F(5.)/_F(13.),
    _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(6.)/_F(13.),
    _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(7.)/_F(13.),
    _F(4.)/_F(13.),  _F(8.)/_F(13.),  _F(1.)/_F(13.),
    _F(4.)/_F(13.),  _F(7.)/_F(13.),  _F(2.)/_F(13.),
    _F(4.)/_F(13.),  _F(6.)/_F(13.),  _F(3.)/_F(13.),
    _F(4.)/_F(13.),  _F(5.)/_F(13.),  _F(4.)/_F(13.),
    _F(4.)/_F(13.),  _F(4.)/_F(13.),  _F(5.)/_F(13.),
    _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(6.)/_F(13.),
    _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(7.)/_F(13.),
    _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(8.)/_F(13.),
    _F(3.)/_F(13.),  _F(9.)/_F(13.),  _F(1.)/_F(13.),
    _F(3.)/_F(13.),  _F(8.)/_F(13.),  _F(2.)/_F(13.),
    _F(3.)/_F(13.),  _F(7.)/_F(13.),  _F(3.)/_F(13.),
    _F(3.)/_F(13.),  _F(6.)/_F(13.),  _F(4.)/_F(13.),
    _F(3.)/_F(13.),  _F(5.)/_F(13.),  _F(5.)/_F(13.),
    _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(6.)/_F(13.),
    _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(7.)/_F(13.),
    _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(8.)/_F(13.),
    _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(9.)/_F(13.),
    _F(2.)/_F(13.), _F(10.)/_F(13.),  _F(1.)/_F(13.),
    _F(2.)/_F(13.),  _F(9.)/_F(13.),  _F(2.)/_F(13.),
    _F(2.)/_F(13.),  _F(8.)/_F(13.),  _F(3.)/_F(13.),
    _F(2.)/_F(13.),  _F(7.)/_F(13.),  _F(4.)/_F(13.),
    _F(2.)/_F(13.),  _F(6.)/_F(13.),  _F(5.)/_F(13.),
    _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(6.)/_F(13.),
    _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(7.)/_F(13.),
    _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(8.)/_F(13.),
    _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(9.)/_F(13.),
    _F(2.)/_F(13.),  _F(1.)/_F(13.), _F(10.)/_F(13.),
    _F(1.)/_F(13.), _F(11.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.), _F(10.)/_F(13.),  _F(2.)/_F(13.),
    _F(1.)/_F(13.),  _F(9.)/_F(13.),  _F(3.)/_F(13.),
    _F(1.)/_F(13.),  _F(8.)/_F(13.),  _F(4.)/_F(13.),
    _F(1.)/_F(13.),  _F(7.)/_F(13.),  _F(5.)/_F(13.),
    _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(6.)/_F(13.),
    _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(7.)/_F(13.),
    _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(8.)/_F(13.),
    _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(9.)/_F(13.),
    _F(1.)/_F(13.),  _F(2.)/_F(13.), _F(10.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.), _F(11.)/_F(13.),
   /* elem */
   _F(10.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
    _F(9.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
    _F(9.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),
    _F(9.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),
    _F(8.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
    _F(8.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),
    _F(8.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),
    _F(8.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),
    _F(8.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),
    _F(8.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),
    _F(7.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
    _F(7.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),
    _F(7.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),
    _F(7.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),
    _F(7.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),
    _F(7.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),
    _F(7.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),
    _F(7.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),
    _F(7.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),
    _F(7.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),
    _F(6.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
    _F(6.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),
    _F(6.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),
    _F(6.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),
    _F(6.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),
    _F(6.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),
    _F(6.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),
    _F(6.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),
    _F(6.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),
    _F(6.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),
    _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),
    _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),
    _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),
    _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),
    _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),
    _F(5.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
    _F(5.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),
    _F(5.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),
    _F(5.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),
    _F(5.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),
    _F(5.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),
    _F(5.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),
    _F(5.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),
    _F(5.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),
    _F(5.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),
    _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),
    _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),
    _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),
    _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),
    _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),
    _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),
    _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),
    _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),
    _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),
    _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),
    _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),
    _F(4.)/_F(13.),  _F(7.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
    _F(4.)/_F(13.),  _F(6.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),
    _F(4.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),
    _F(4.)/_F(13.),  _F(5.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),
    _F(4.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),
    _F(4.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),
    _F(4.)/_F(13.),  _F(4.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),
    _F(4.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),
    _F(4.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),
    _F(4.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),
    _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),
    _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),
    _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),
    _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),
    _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),
    _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),
    _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),
    _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),
    _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),
    _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),
    _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),
    _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(7.)/_F(13.),  _F(1.)/_F(13.),
    _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(2.)/_F(13.),
    _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(3.)/_F(13.),
    _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(4.)/_F(13.),
    _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(5.)/_F(13.),
    _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(6.)/_F(13.),
    _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(7.)/_F(13.),
    _F(3.)/_F(13.),  _F(8.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
    _F(3.)/_F(13.),  _F(7.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),
    _F(3.)/_F(13.),  _F(7.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),
    _F(3.)/_F(13.),  _F(6.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),
    _F(3.)/_F(13.),  _F(6.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),
    _F(3.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),
    _F(3.)/_F(13.),  _F(5.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),
    _F(3.)/_F(13.),  _F(5.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),
    _F(3.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),
    _F(3.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),
    _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),
    _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),
    _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),
    _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),
    _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),
    _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),
    _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),
    _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),
    _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),
    _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),
    _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),
    _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(7.)/_F(13.),  _F(1.)/_F(13.),
    _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(6.)/_F(13.),  _F(2.)/_F(13.),
    _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(3.)/_F(13.),
    _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(4.)/_F(13.),
    _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(5.)/_F(13.),
    _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(6.)/_F(13.),
    _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(7.)/_F(13.),
    _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(8.)/_F(13.),  _F(1.)/_F(13.),
    _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(7.)/_F(13.),  _F(2.)/_F(13.),
    _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(3.)/_F(13.),
    _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(4.)/_F(13.),
    _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(5.)/_F(13.),
    _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(6.)/_F(13.),
    _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(7.)/_F(13.),
    _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(8.)/_F(13.),
    _F(2.)/_F(13.),  _F(9.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
    _F(2.)/_F(13.),  _F(8.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),
    _F(2.)/_F(13.),  _F(8.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),
    _F(2.)/_F(13.),  _F(7.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),
    _F(2.)/_F(13.),  _F(7.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),
    _F(2.)/_F(13.),  _F(7.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),
    _F(2.)/_F(13.),  _F(6.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),
    _F(2.)/_F(13.),  _F(6.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),
    _F(2.)/_F(13.),  _F(6.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),
    _F(2.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),
    _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),
    _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),
    _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),
    _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),
    _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),
    _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),
    _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),
    _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),
    _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),
    _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),
    _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),
    _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(7.)/_F(13.),  _F(1.)/_F(13.),
    _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(6.)/_F(13.),  _F(2.)/_F(13.),
    _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(5.)/_F(13.),  _F(3.)/_F(13.),
    _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(4.)/_F(13.),
    _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(5.)/_F(13.),
    _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(6.)/_F(13.),
    _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(7.)/_F(13.),
    _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(8.)/_F(13.),  _F(1.)/_F(13.),
    _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(7.)/_F(13.),  _F(2.)/_F(13.),
    _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(6.)/_F(13.),  _F(3.)/_F(13.),
    _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(4.)/_F(13.),
    _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(5.)/_F(13.),
    _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(6.)/_F(13.),
    _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(7.)/_F(13.),
    _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(8.)/_F(13.),
    _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(9.)/_F(13.),  _F(1.)/_F(13.),
    _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(8.)/_F(13.),  _F(2.)/_F(13.),
    _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(7.)/_F(13.),  _F(3.)/_F(13.),
    _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(4.)/_F(13.),
    _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(5.)/_F(13.),
    _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(6.)/_F(13.),
    _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(7.)/_F(13.),
    _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(8.)/_F(13.),
    _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(9.)/_F(13.),
    _F(1.)/_F(13.), _F(10.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.),  _F(9.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.),  _F(9.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),
    _F(1.)/_F(13.),  _F(8.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.),  _F(8.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),
    _F(1.)/_F(13.),  _F(8.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),
    _F(1.)/_F(13.),  _F(7.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.),  _F(7.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),
    _F(1.)/_F(13.),  _F(7.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),
    _F(1.)/_F(13.),  _F(7.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),
    _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),
    _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),
    _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),
    _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),
    _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(6.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),
    _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),
    _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),
    _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),
    _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),
    _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(7.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(6.)/_F(13.),  _F(2.)/_F(13.),
    _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(5.)/_F(13.),  _F(3.)/_F(13.),
    _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(4.)/_F(13.),  _F(4.)/_F(13.),
    _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(3.)/_F(13.),  _F(5.)/_F(13.),
    _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(2.)/_F(13.),  _F(6.)/_F(13.),
    _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(1.)/_F(13.),  _F(7.)/_F(13.),
    _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(8.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(7.)/_F(13.),  _F(2.)/_F(13.),
    _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(6.)/_F(13.),  _F(3.)/_F(13.),
    _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(5.)/_F(13.),  _F(4.)/_F(13.),
    _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(4.)/_F(13.),  _F(5.)/_F(13.),
    _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(3.)/_F(13.),  _F(6.)/_F(13.),
    _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(2.)/_F(13.),  _F(7.)/_F(13.),
    _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(1.)/_F(13.),  _F(8.)/_F(13.),
    _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(9.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(8.)/_F(13.),  _F(2.)/_F(13.),
    _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(7.)/_F(13.),  _F(3.)/_F(13.),
    _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(6.)/_F(13.),  _F(4.)/_F(13.),
    _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(5.)/_F(13.),  _F(5.)/_F(13.),
    _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(4.)/_F(13.),  _F(6.)/_F(13.),
    _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(3.)/_F(13.),  _F(7.)/_F(13.),
    _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(2.)/_F(13.),  _F(8.)/_F(13.),
    _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(1.)/_F(13.),  _F(9.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.), _F(10.)/_F(13.),  _F(1.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(9.)/_F(13.),  _F(2.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(8.)/_F(13.),  _F(3.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(7.)/_F(13.),  _F(4.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(6.)/_F(13.),  _F(5.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(5.)/_F(13.),  _F(6.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(4.)/_F(13.),  _F(7.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(3.)/_F(13.),  _F(8.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(2.)/_F(13.),  _F(9.)/_F(13.),
    _F(1.)/_F(13.),  _F(1.)/_F(13.),  _F(1.)/_F(13.), _F(10.)/_F(13.)
};
DOF_TYPE DOF_P13_ = {
    DofReserved, "P13", P13_points, NULL, DOF_DG12, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    560, 13, 0, 0, 1, 1, 12, 66, 220
};

static FLOAT P14_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
   _F(13.)/_F(14.),  _F(1.)/_F(14.),
   _F(12.)/_F(14.),  _F(2.)/_F(14.),
   _F(11.)/_F(14.),  _F(3.)/_F(14.),
   _F(10.)/_F(14.),  _F(4.)/_F(14.),
    _F(9.)/_F(14.),  _F(5.)/_F(14.),
    _F(8.)/_F(14.),  _F(6.)/_F(14.),
    _F(7.)/_F(14.),  _F(7.)/_F(14.),
    _F(6.)/_F(14.),  _F(8.)/_F(14.),
    _F(5.)/_F(14.),  _F(9.)/_F(14.),
    _F(4.)/_F(14.), _F(10.)/_F(14.),
    _F(3.)/_F(14.), _F(11.)/_F(14.),
    _F(2.)/_F(14.), _F(12.)/_F(14.),
    _F(1.)/_F(14.), _F(13.)/_F(14.),
   /* face */
   _F(12.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
   _F(11.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
   _F(11.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
   _F(10.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),
   _F(10.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),
   _F(10.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),
    _F(9.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),
    _F(9.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),
    _F(9.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),
    _F(9.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),
    _F(8.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),
    _F(8.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),
    _F(8.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),
    _F(8.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),
    _F(8.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),
    _F(7.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),
    _F(7.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),
    _F(7.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),
    _F(7.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),
    _F(7.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),
    _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),
    _F(6.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),
    _F(6.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),
    _F(6.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),
    _F(6.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),
    _F(6.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),
    _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),
    _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),
    _F(5.)/_F(14.),  _F(8.)/_F(14.),  _F(1.)/_F(14.),
    _F(5.)/_F(14.),  _F(7.)/_F(14.),  _F(2.)/_F(14.),
    _F(5.)/_F(14.),  _F(6.)/_F(14.),  _F(3.)/_F(14.),
    _F(5.)/_F(14.),  _F(5.)/_F(14.),  _F(4.)/_F(14.),
    _F(5.)/_F(14.),  _F(4.)/_F(14.),  _F(5.)/_F(14.),
    _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(6.)/_F(14.),
    _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(7.)/_F(14.),
    _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(8.)/_F(14.),
    _F(4.)/_F(14.),  _F(9.)/_F(14.),  _F(1.)/_F(14.),
    _F(4.)/_F(14.),  _F(8.)/_F(14.),  _F(2.)/_F(14.),
    _F(4.)/_F(14.),  _F(7.)/_F(14.),  _F(3.)/_F(14.),
    _F(4.)/_F(14.),  _F(6.)/_F(14.),  _F(4.)/_F(14.),
    _F(4.)/_F(14.),  _F(5.)/_F(14.),  _F(5.)/_F(14.),
    _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(6.)/_F(14.),
    _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(7.)/_F(14.),
    _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(8.)/_F(14.),
    _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(9.)/_F(14.),
    _F(3.)/_F(14.), _F(10.)/_F(14.),  _F(1.)/_F(14.),
    _F(3.)/_F(14.),  _F(9.)/_F(14.),  _F(2.)/_F(14.),
    _F(3.)/_F(14.),  _F(8.)/_F(14.),  _F(3.)/_F(14.),
    _F(3.)/_F(14.),  _F(7.)/_F(14.),  _F(4.)/_F(14.),
    _F(3.)/_F(14.),  _F(6.)/_F(14.),  _F(5.)/_F(14.),
    _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(6.)/_F(14.),
    _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(7.)/_F(14.),
    _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(8.)/_F(14.),
    _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(9.)/_F(14.),
    _F(3.)/_F(14.),  _F(1.)/_F(14.), _F(10.)/_F(14.),
    _F(2.)/_F(14.), _F(11.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.), _F(10.)/_F(14.),  _F(2.)/_F(14.),
    _F(2.)/_F(14.),  _F(9.)/_F(14.),  _F(3.)/_F(14.),
    _F(2.)/_F(14.),  _F(8.)/_F(14.),  _F(4.)/_F(14.),
    _F(2.)/_F(14.),  _F(7.)/_F(14.),  _F(5.)/_F(14.),
    _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(6.)/_F(14.),
    _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(7.)/_F(14.),
    _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(8.)/_F(14.),
    _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(9.)/_F(14.),
    _F(2.)/_F(14.),  _F(2.)/_F(14.), _F(10.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.), _F(11.)/_F(14.),
    _F(1.)/_F(14.), _F(12.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.), _F(11.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.), _F(10.)/_F(14.),  _F(3.)/_F(14.),
    _F(1.)/_F(14.),  _F(9.)/_F(14.),  _F(4.)/_F(14.),
    _F(1.)/_F(14.),  _F(8.)/_F(14.),  _F(5.)/_F(14.),
    _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(6.)/_F(14.),
    _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(7.)/_F(14.),
    _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(8.)/_F(14.),
    _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(9.)/_F(14.),
    _F(1.)/_F(14.),  _F(3.)/_F(14.), _F(10.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.), _F(11.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.), _F(12.)/_F(14.),
   /* elem */
   _F(11.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
   _F(10.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
   _F(10.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
   _F(10.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
    _F(9.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
    _F(9.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
    _F(9.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
    _F(9.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),
    _F(9.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),
    _F(9.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),
    _F(8.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
    _F(8.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
    _F(8.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
    _F(8.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),
    _F(8.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),
    _F(8.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),
    _F(8.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),
    _F(8.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),
    _F(8.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),
    _F(8.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),
    _F(7.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
    _F(7.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
    _F(7.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
    _F(7.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),
    _F(7.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),
    _F(7.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),
    _F(7.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),
    _F(7.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),
    _F(7.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),
    _F(7.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),
    _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),
    _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),
    _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),
    _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),
    _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),
    _F(6.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
    _F(6.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
    _F(6.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
    _F(6.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),
    _F(6.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),
    _F(6.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),
    _F(6.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),
    _F(6.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),
    _F(6.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),
    _F(6.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),
    _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),
    _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),
    _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),
    _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),
    _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),
    _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),
    _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),
    _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),
    _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),
    _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),
    _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),
    _F(5.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
    _F(5.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
    _F(5.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
    _F(5.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),
    _F(5.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),
    _F(5.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),
    _F(5.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),
    _F(5.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),
    _F(5.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),
    _F(5.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),
    _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),
    _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),
    _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),
    _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),
    _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),
    _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),
    _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),
    _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),
    _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),
    _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),
    _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),
    _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),
    _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),
    _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),
    _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),
    _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),
    _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),
    _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),
    _F(4.)/_F(14.),  _F(8.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
    _F(4.)/_F(14.),  _F(7.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
    _F(4.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
    _F(4.)/_F(14.),  _F(6.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),
    _F(4.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),
    _F(4.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),
    _F(4.)/_F(14.),  _F(5.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),
    _F(4.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),
    _F(4.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),
    _F(4.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),
    _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),
    _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),
    _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),
    _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),
    _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),
    _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),
    _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),
    _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),
    _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),
    _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),
    _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),
    _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),
    _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),
    _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),
    _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),
    _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),
    _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),
    _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),
    _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(8.)/_F(14.),  _F(1.)/_F(14.),
    _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(2.)/_F(14.),
    _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(3.)/_F(14.),
    _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(4.)/_F(14.),
    _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(5.)/_F(14.),
    _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(6.)/_F(14.),
    _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(7.)/_F(14.),
    _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(8.)/_F(14.),
    _F(3.)/_F(14.),  _F(9.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
    _F(3.)/_F(14.),  _F(8.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
    _F(3.)/_F(14.),  _F(8.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
    _F(3.)/_F(14.),  _F(7.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),
    _F(3.)/_F(14.),  _F(7.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),
    _F(3.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),
    _F(3.)/_F(14.),  _F(6.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),
    _F(3.)/_F(14.),  _F(6.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),
    _F(3.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),
    _F(3.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),
    _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),
    _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),
    _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),
    _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),
    _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),
    _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),
    _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),
    _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),
    _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),
    _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),
    _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),
    _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),
    _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),
    _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),
    _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),
    _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),
    _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),
    _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),
    _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(8.)/_F(14.),  _F(1.)/_F(14.),
    _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(7.)/_F(14.),  _F(2.)/_F(14.),
    _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(3.)/_F(14.),
    _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(4.)/_F(14.),
    _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(5.)/_F(14.),
    _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(6.)/_F(14.),
    _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(7.)/_F(14.),
    _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(8.)/_F(14.),
    _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(9.)/_F(14.),  _F(1.)/_F(14.),
    _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(8.)/_F(14.),  _F(2.)/_F(14.),
    _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(3.)/_F(14.),
    _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(4.)/_F(14.),
    _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(5.)/_F(14.),
    _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(6.)/_F(14.),
    _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(7.)/_F(14.),
    _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(8.)/_F(14.),
    _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(9.)/_F(14.),
    _F(2.)/_F(14.), _F(10.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.),  _F(9.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.),  _F(9.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
    _F(2.)/_F(14.),  _F(8.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.),  _F(8.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),
    _F(2.)/_F(14.),  _F(8.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),
    _F(2.)/_F(14.),  _F(7.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.),  _F(7.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),
    _F(2.)/_F(14.),  _F(7.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),
    _F(2.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),
    _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),
    _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),
    _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),
    _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),
    _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),
    _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),
    _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),
    _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),
    _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),
    _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),
    _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),
    _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),
    _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),
    _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),
    _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),
    _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(8.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(7.)/_F(14.),  _F(2.)/_F(14.),
    _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(6.)/_F(14.),  _F(3.)/_F(14.),
    _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(4.)/_F(14.),
    _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(5.)/_F(14.),
    _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(6.)/_F(14.),
    _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(7.)/_F(14.),
    _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(8.)/_F(14.),
    _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(9.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(8.)/_F(14.),  _F(2.)/_F(14.),
    _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(7.)/_F(14.),  _F(3.)/_F(14.),
    _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(4.)/_F(14.),
    _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(5.)/_F(14.),
    _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(6.)/_F(14.),
    _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(7.)/_F(14.),
    _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(8.)/_F(14.),
    _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(9.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.), _F(10.)/_F(14.),  _F(1.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(9.)/_F(14.),  _F(2.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(8.)/_F(14.),  _F(3.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(4.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(5.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(6.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(7.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(8.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(9.)/_F(14.),
    _F(2.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.), _F(10.)/_F(14.),
    _F(1.)/_F(14.), _F(11.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.), _F(10.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.), _F(10.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.),  _F(9.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.),  _F(9.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.),  _F(9.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),
    _F(1.)/_F(14.),  _F(8.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.),  _F(8.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.),  _F(8.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),
    _F(1.)/_F(14.),  _F(8.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),
    _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),
    _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),
    _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),
    _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),
    _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),
    _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),
    _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),
    _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(7.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(6.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),
    _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),
    _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),
    _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),
    _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),
    _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(8.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(7.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(6.)/_F(14.),  _F(3.)/_F(14.),
    _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(5.)/_F(14.),  _F(4.)/_F(14.),
    _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(4.)/_F(14.),  _F(5.)/_F(14.),
    _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(3.)/_F(14.),  _F(6.)/_F(14.),
    _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(2.)/_F(14.),  _F(7.)/_F(14.),
    _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(1.)/_F(14.),  _F(8.)/_F(14.),
    _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(9.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(8.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(7.)/_F(14.),  _F(3.)/_F(14.),
    _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(6.)/_F(14.),  _F(4.)/_F(14.),
    _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(5.)/_F(14.),  _F(5.)/_F(14.),
    _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(4.)/_F(14.),  _F(6.)/_F(14.),
    _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(3.)/_F(14.),  _F(7.)/_F(14.),
    _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(2.)/_F(14.),  _F(8.)/_F(14.),
    _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(1.)/_F(14.),  _F(9.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.), _F(10.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(9.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(8.)/_F(14.),  _F(3.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(7.)/_F(14.),  _F(4.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(6.)/_F(14.),  _F(5.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(5.)/_F(14.),  _F(6.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(4.)/_F(14.),  _F(7.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(3.)/_F(14.),  _F(8.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(2.)/_F(14.),  _F(9.)/_F(14.),
    _F(1.)/_F(14.),  _F(2.)/_F(14.),  _F(1.)/_F(14.), _F(10.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.), _F(11.)/_F(14.),  _F(1.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.), _F(10.)/_F(14.),  _F(2.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(9.)/_F(14.),  _F(3.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(8.)/_F(14.),  _F(4.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(7.)/_F(14.),  _F(5.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(6.)/_F(14.),  _F(6.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(5.)/_F(14.),  _F(7.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(4.)/_F(14.),  _F(8.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(3.)/_F(14.),  _F(9.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(2.)/_F(14.), _F(10.)/_F(14.),
    _F(1.)/_F(14.),  _F(1.)/_F(14.),  _F(1.)/_F(14.), _F(11.)/_F(14.)
};
DOF_TYPE DOF_P14_ = {
    DofReserved, "P14", P14_points, NULL, DOF_DG13, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    680, 14, 0, 0, 1, 1, 13, 78, 286
};

static FLOAT P15_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
   _F(14.)/_F(15.),  _F(1.)/_F(15.),
   _F(13.)/_F(15.),  _F(2.)/_F(15.),
   _F(12.)/_F(15.),  _F(3.)/_F(15.),
   _F(11.)/_F(15.),  _F(4.)/_F(15.),
   _F(10.)/_F(15.),  _F(5.)/_F(15.),
    _F(9.)/_F(15.),  _F(6.)/_F(15.),
    _F(8.)/_F(15.),  _F(7.)/_F(15.),
    _F(7.)/_F(15.),  _F(8.)/_F(15.),
    _F(6.)/_F(15.),  _F(9.)/_F(15.),
    _F(5.)/_F(15.), _F(10.)/_F(15.),
    _F(4.)/_F(15.), _F(11.)/_F(15.),
    _F(3.)/_F(15.), _F(12.)/_F(15.),
    _F(2.)/_F(15.), _F(13.)/_F(15.),
    _F(1.)/_F(15.), _F(14.)/_F(15.),
   /* face */
   _F(13.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
   _F(12.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
   _F(12.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
   _F(11.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
   _F(11.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
   _F(11.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
   _F(10.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),
   _F(10.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),
   _F(10.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),
   _F(10.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),
    _F(9.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),
    _F(9.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),
    _F(9.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),
    _F(9.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),
    _F(9.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),
    _F(8.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),
    _F(8.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),
    _F(8.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),
    _F(8.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),
    _F(8.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),
    _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),
    _F(7.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),
    _F(7.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),
    _F(7.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),
    _F(7.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),
    _F(7.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),
    _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),
    _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),
    _F(6.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),
    _F(6.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),
    _F(6.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),
    _F(6.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),
    _F(6.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),
    _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),
    _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),
    _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),
    _F(5.)/_F(15.),  _F(9.)/_F(15.),  _F(1.)/_F(15.),
    _F(5.)/_F(15.),  _F(8.)/_F(15.),  _F(2.)/_F(15.),
    _F(5.)/_F(15.),  _F(7.)/_F(15.),  _F(3.)/_F(15.),
    _F(5.)/_F(15.),  _F(6.)/_F(15.),  _F(4.)/_F(15.),
    _F(5.)/_F(15.),  _F(5.)/_F(15.),  _F(5.)/_F(15.),
    _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(6.)/_F(15.),
    _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(7.)/_F(15.),
    _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(8.)/_F(15.),
    _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(9.)/_F(15.),
    _F(4.)/_F(15.), _F(10.)/_F(15.),  _F(1.)/_F(15.),
    _F(4.)/_F(15.),  _F(9.)/_F(15.),  _F(2.)/_F(15.),
    _F(4.)/_F(15.),  _F(8.)/_F(15.),  _F(3.)/_F(15.),
    _F(4.)/_F(15.),  _F(7.)/_F(15.),  _F(4.)/_F(15.),
    _F(4.)/_F(15.),  _F(6.)/_F(15.),  _F(5.)/_F(15.),
    _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(6.)/_F(15.),
    _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(7.)/_F(15.),
    _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(8.)/_F(15.),
    _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(9.)/_F(15.),
    _F(4.)/_F(15.),  _F(1.)/_F(15.), _F(10.)/_F(15.),
    _F(3.)/_F(15.), _F(11.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.), _F(10.)/_F(15.),  _F(2.)/_F(15.),
    _F(3.)/_F(15.),  _F(9.)/_F(15.),  _F(3.)/_F(15.),
    _F(3.)/_F(15.),  _F(8.)/_F(15.),  _F(4.)/_F(15.),
    _F(3.)/_F(15.),  _F(7.)/_F(15.),  _F(5.)/_F(15.),
    _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(6.)/_F(15.),
    _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(7.)/_F(15.),
    _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(8.)/_F(15.),
    _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(9.)/_F(15.),
    _F(3.)/_F(15.),  _F(2.)/_F(15.), _F(10.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.), _F(11.)/_F(15.),
    _F(2.)/_F(15.), _F(12.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.), _F(11.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.), _F(10.)/_F(15.),  _F(3.)/_F(15.),
    _F(2.)/_F(15.),  _F(9.)/_F(15.),  _F(4.)/_F(15.),
    _F(2.)/_F(15.),  _F(8.)/_F(15.),  _F(5.)/_F(15.),
    _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(6.)/_F(15.),
    _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(7.)/_F(15.),
    _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(8.)/_F(15.),
    _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(9.)/_F(15.),
    _F(2.)/_F(15.),  _F(3.)/_F(15.), _F(10.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.), _F(11.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.), _F(12.)/_F(15.),
    _F(1.)/_F(15.), _F(13.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.), _F(12.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.), _F(11.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.), _F(10.)/_F(15.),  _F(4.)/_F(15.),
    _F(1.)/_F(15.),  _F(9.)/_F(15.),  _F(5.)/_F(15.),
    _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(6.)/_F(15.),
    _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(7.)/_F(15.),
    _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(8.)/_F(15.),
    _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(9.)/_F(15.),
    _F(1.)/_F(15.),  _F(4.)/_F(15.), _F(10.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.), _F(11.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.), _F(12.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.), _F(13.)/_F(15.),
   /* elem */
   _F(12.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
   _F(11.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
   _F(11.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
   _F(11.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
   _F(10.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
   _F(10.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
   _F(10.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
   _F(10.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
   _F(10.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
   _F(10.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
    _F(9.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
    _F(9.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
    _F(9.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
    _F(9.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
    _F(9.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
    _F(9.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
    _F(9.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),
    _F(9.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),
    _F(9.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),
    _F(9.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),
    _F(8.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
    _F(8.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
    _F(8.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
    _F(8.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
    _F(8.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
    _F(8.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
    _F(8.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),
    _F(8.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),
    _F(8.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),
    _F(8.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),
    _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),
    _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),
    _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),
    _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),
    _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),
    _F(7.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
    _F(7.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
    _F(7.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
    _F(7.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
    _F(7.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
    _F(7.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
    _F(7.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),
    _F(7.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),
    _F(7.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),
    _F(7.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),
    _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),
    _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),
    _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),
    _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),
    _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),
    _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),
    _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),
    _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),
    _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),
    _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),
    _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),
    _F(6.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
    _F(6.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
    _F(6.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
    _F(6.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
    _F(6.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
    _F(6.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
    _F(6.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),
    _F(6.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),
    _F(6.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),
    _F(6.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),
    _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),
    _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),
    _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),
    _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),
    _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),
    _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),
    _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),
    _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),
    _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),
    _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),
    _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),
    _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),
    _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),
    _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),
    _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),
    _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),
    _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),
    _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),
    _F(5.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
    _F(5.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
    _F(5.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
    _F(5.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
    _F(5.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
    _F(5.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
    _F(5.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),
    _F(5.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),
    _F(5.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),
    _F(5.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),
    _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),
    _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),
    _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),
    _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),
    _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),
    _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),
    _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),
    _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),
    _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),
    _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),
    _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),
    _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),
    _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),
    _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),
    _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),
    _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),
    _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),
    _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),
    _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),
    _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),
    _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),
    _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),
    _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),
    _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),
    _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),
    _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),
    _F(4.)/_F(15.),  _F(9.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
    _F(4.)/_F(15.),  _F(8.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
    _F(4.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
    _F(4.)/_F(15.),  _F(7.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
    _F(4.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
    _F(4.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
    _F(4.)/_F(15.),  _F(6.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),
    _F(4.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),
    _F(4.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),
    _F(4.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),
    _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),
    _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),
    _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),
    _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),
    _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),
    _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),
    _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),
    _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),
    _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),
    _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),
    _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),
    _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),
    _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),
    _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),
    _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),
    _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),
    _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),
    _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),
    _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),
    _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),
    _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),
    _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),
    _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),
    _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),
    _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),
    _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),
    _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(9.)/_F(15.),  _F(1.)/_F(15.),
    _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(2.)/_F(15.),
    _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(3.)/_F(15.),
    _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(4.)/_F(15.),
    _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(5.)/_F(15.),
    _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(6.)/_F(15.),
    _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(7.)/_F(15.),
    _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(8.)/_F(15.),
    _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(9.)/_F(15.),
    _F(3.)/_F(15.), _F(10.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.),  _F(9.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.),  _F(9.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
    _F(3.)/_F(15.),  _F(8.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.),  _F(8.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
    _F(3.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
    _F(3.)/_F(15.),  _F(7.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.),  _F(7.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),
    _F(3.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),
    _F(3.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),
    _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),
    _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),
    _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),
    _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),
    _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),
    _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),
    _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),
    _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),
    _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),
    _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),
    _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),
    _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),
    _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),
    _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),
    _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),
    _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),
    _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),
    _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),
    _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),
    _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),
    _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),
    _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),
    _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(9.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(8.)/_F(15.),  _F(2.)/_F(15.),
    _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(3.)/_F(15.),
    _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(4.)/_F(15.),
    _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(5.)/_F(15.),
    _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(6.)/_F(15.),
    _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(7.)/_F(15.),
    _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(8.)/_F(15.),
    _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(9.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.), _F(10.)/_F(15.),  _F(1.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(9.)/_F(15.),  _F(2.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(3.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(4.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(5.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(6.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(7.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(8.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(9.)/_F(15.),
    _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.), _F(10.)/_F(15.),
    _F(2.)/_F(15.), _F(11.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.), _F(10.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.), _F(10.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.),  _F(9.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.),  _F(9.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.),  _F(9.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
    _F(2.)/_F(15.),  _F(8.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.),  _F(8.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.),  _F(8.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),
    _F(2.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),
    _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),
    _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),
    _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),
    _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),
    _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),
    _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),
    _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),
    _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),
    _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),
    _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),
    _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),
    _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),
    _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),
    _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),
    _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),
    _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),
    _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),
    _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),
    _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(9.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(8.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(7.)/_F(15.),  _F(3.)/_F(15.),
    _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(4.)/_F(15.),
    _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(5.)/_F(15.),
    _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(6.)/_F(15.),
    _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(7.)/_F(15.),
    _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(8.)/_F(15.),
    _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),  _F(9.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.), _F(10.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(9.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(8.)/_F(15.),  _F(3.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(4.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(5.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(6.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(7.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(8.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(9.)/_F(15.),
    _F(2.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.), _F(10.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.), _F(11.)/_F(15.),  _F(1.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.), _F(10.)/_F(15.),  _F(2.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(9.)/_F(15.),  _F(3.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(4.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(5.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(6.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(7.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(8.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(9.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.), _F(10.)/_F(15.),
    _F(2.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.), _F(11.)/_F(15.),
    _F(1.)/_F(15.), _F(12.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.), _F(11.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.), _F(11.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.), _F(10.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.), _F(10.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.), _F(10.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.),  _F(9.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.),  _F(9.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.),  _F(9.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.),  _F(9.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),
    _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),
    _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),
    _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),
    _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),
    _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),
    _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(7.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),
    _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),
    _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),
    _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),
    _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(8.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(7.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(6.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),
    _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),
    _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),
    _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),
    _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),
    _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(9.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(8.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(7.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(6.)/_F(15.),  _F(4.)/_F(15.),
    _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(5.)/_F(15.),  _F(5.)/_F(15.),
    _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(4.)/_F(15.),  _F(6.)/_F(15.),
    _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(3.)/_F(15.),  _F(7.)/_F(15.),
    _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(2.)/_F(15.),  _F(8.)/_F(15.),
    _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(1.)/_F(15.),  _F(9.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.), _F(10.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(9.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(8.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(7.)/_F(15.),  _F(4.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(6.)/_F(15.),  _F(5.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(5.)/_F(15.),  _F(6.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(4.)/_F(15.),  _F(7.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(3.)/_F(15.),  _F(8.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(2.)/_F(15.),  _F(9.)/_F(15.),
    _F(1.)/_F(15.),  _F(3.)/_F(15.),  _F(1.)/_F(15.), _F(10.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.), _F(11.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.), _F(10.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(9.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(8.)/_F(15.),  _F(4.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(7.)/_F(15.),  _F(5.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(6.)/_F(15.),  _F(6.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(5.)/_F(15.),  _F(7.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(4.)/_F(15.),  _F(8.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(3.)/_F(15.),  _F(9.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(2.)/_F(15.), _F(10.)/_F(15.),
    _F(1.)/_F(15.),  _F(2.)/_F(15.),  _F(1.)/_F(15.), _F(11.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.), _F(12.)/_F(15.),  _F(1.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.), _F(11.)/_F(15.),  _F(2.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.), _F(10.)/_F(15.),  _F(3.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(9.)/_F(15.),  _F(4.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(8.)/_F(15.),  _F(5.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(7.)/_F(15.),  _F(6.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(6.)/_F(15.),  _F(7.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(5.)/_F(15.),  _F(8.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(4.)/_F(15.),  _F(9.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(3.)/_F(15.), _F(10.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(2.)/_F(15.), _F(11.)/_F(15.),
    _F(1.)/_F(15.),  _F(1.)/_F(15.),  _F(1.)/_F(15.), _F(12.)/_F(15.)
};
DOF_TYPE DOF_P15_ = {
    DofReserved, "P15", P15_points, NULL, DOF_DG14, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    816, 15, 0, 0, 1, 1, 14, 91, 364
};

static FLOAT P16_points[] = {
   /* vertex */
    _F(1.),
   /* edge */
   _F(15.)/_F(16.),  _F(1.)/_F(16.),
   _F(14.)/_F(16.),  _F(2.)/_F(16.),
   _F(13.)/_F(16.),  _F(3.)/_F(16.),
   _F(12.)/_F(16.),  _F(4.)/_F(16.),
   _F(11.)/_F(16.),  _F(5.)/_F(16.),
   _F(10.)/_F(16.),  _F(6.)/_F(16.),
    _F(9.)/_F(16.),  _F(7.)/_F(16.),
    _F(8.)/_F(16.),  _F(8.)/_F(16.),
    _F(7.)/_F(16.),  _F(9.)/_F(16.),
    _F(6.)/_F(16.), _F(10.)/_F(16.),
    _F(5.)/_F(16.), _F(11.)/_F(16.),
    _F(4.)/_F(16.), _F(12.)/_F(16.),
    _F(3.)/_F(16.), _F(13.)/_F(16.),
    _F(2.)/_F(16.), _F(14.)/_F(16.),
    _F(1.)/_F(16.), _F(15.)/_F(16.),
   /* face */
   _F(14.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
   _F(13.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
   _F(13.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
   _F(12.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
   _F(12.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
   _F(12.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
   _F(11.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
   _F(11.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
   _F(11.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
   _F(11.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
   _F(10.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),
   _F(10.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),
   _F(10.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),
   _F(10.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),
   _F(10.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),
    _F(9.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),
    _F(9.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),
    _F(9.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),
    _F(9.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),
    _F(9.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),
    _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),
    _F(8.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),
    _F(8.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),
    _F(8.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),
    _F(8.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),
    _F(8.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),
    _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),
    _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),
    _F(7.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),
    _F(7.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),
    _F(7.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),
    _F(7.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),
    _F(7.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),
    _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),
    _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),
    _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),
    _F(6.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),
    _F(6.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),
    _F(6.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),
    _F(6.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),
    _F(6.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),
    _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),
    _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),
    _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),
    _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),
    _F(5.)/_F(16.), _F(10.)/_F(16.),  _F(1.)/_F(16.),
    _F(5.)/_F(16.),  _F(9.)/_F(16.),  _F(2.)/_F(16.),
    _F(5.)/_F(16.),  _F(8.)/_F(16.),  _F(3.)/_F(16.),
    _F(5.)/_F(16.),  _F(7.)/_F(16.),  _F(4.)/_F(16.),
    _F(5.)/_F(16.),  _F(6.)/_F(16.),  _F(5.)/_F(16.),
    _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(6.)/_F(16.),
    _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(7.)/_F(16.),
    _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(8.)/_F(16.),
    _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(9.)/_F(16.),
    _F(5.)/_F(16.),  _F(1.)/_F(16.), _F(10.)/_F(16.),
    _F(4.)/_F(16.), _F(11.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.), _F(10.)/_F(16.),  _F(2.)/_F(16.),
    _F(4.)/_F(16.),  _F(9.)/_F(16.),  _F(3.)/_F(16.),
    _F(4.)/_F(16.),  _F(8.)/_F(16.),  _F(4.)/_F(16.),
    _F(4.)/_F(16.),  _F(7.)/_F(16.),  _F(5.)/_F(16.),
    _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(6.)/_F(16.),
    _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(7.)/_F(16.),
    _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(8.)/_F(16.),
    _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(9.)/_F(16.),
    _F(4.)/_F(16.),  _F(2.)/_F(16.), _F(10.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.), _F(11.)/_F(16.),
    _F(3.)/_F(16.), _F(12.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.), _F(11.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.), _F(10.)/_F(16.),  _F(3.)/_F(16.),
    _F(3.)/_F(16.),  _F(9.)/_F(16.),  _F(4.)/_F(16.),
    _F(3.)/_F(16.),  _F(8.)/_F(16.),  _F(5.)/_F(16.),
    _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(6.)/_F(16.),
    _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(7.)/_F(16.),
    _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(8.)/_F(16.),
    _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(9.)/_F(16.),
    _F(3.)/_F(16.),  _F(3.)/_F(16.), _F(10.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.), _F(11.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.), _F(12.)/_F(16.),
    _F(2.)/_F(16.), _F(13.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.), _F(12.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.), _F(11.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.), _F(10.)/_F(16.),  _F(4.)/_F(16.),
    _F(2.)/_F(16.),  _F(9.)/_F(16.),  _F(5.)/_F(16.),
    _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(6.)/_F(16.),
    _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(7.)/_F(16.),
    _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(8.)/_F(16.),
    _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(9.)/_F(16.),
    _F(2.)/_F(16.),  _F(4.)/_F(16.), _F(10.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.), _F(11.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.), _F(12.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.), _F(13.)/_F(16.),
    _F(1.)/_F(16.), _F(14.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.), _F(13.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.), _F(12.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.), _F(11.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.), _F(10.)/_F(16.),  _F(5.)/_F(16.),
    _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(6.)/_F(16.),
    _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(7.)/_F(16.),
    _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(8.)/_F(16.),
    _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(9.)/_F(16.),
    _F(1.)/_F(16.),  _F(5.)/_F(16.), _F(10.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.), _F(11.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.), _F(12.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.), _F(13.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.), _F(14.)/_F(16.),
   /* elem */
   _F(13.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
   _F(12.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
   _F(12.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
   _F(12.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
   _F(11.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
   _F(11.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
   _F(11.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
   _F(11.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
   _F(11.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
   _F(11.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
   _F(10.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
   _F(10.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
   _F(10.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
   _F(10.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
   _F(10.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
   _F(10.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
   _F(10.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
   _F(10.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
   _F(10.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
   _F(10.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
    _F(9.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
    _F(9.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
    _F(9.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
    _F(9.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
    _F(9.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
    _F(9.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
    _F(9.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
    _F(9.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
    _F(9.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
    _F(9.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
    _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),
    _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),
    _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),
    _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),
    _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),
    _F(8.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
    _F(8.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
    _F(8.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
    _F(8.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
    _F(8.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
    _F(8.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
    _F(8.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
    _F(8.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
    _F(8.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
    _F(8.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
    _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),
    _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),
    _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),
    _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),
    _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),
    _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),
    _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),
    _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),
    _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),
    _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),
    _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),
    _F(7.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
    _F(7.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
    _F(7.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
    _F(7.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
    _F(7.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
    _F(7.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
    _F(7.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
    _F(7.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
    _F(7.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
    _F(7.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
    _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),
    _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),
    _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),
    _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),
    _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),
    _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),
    _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),
    _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),
    _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),
    _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),
    _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),
    _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),
    _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),
    _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),
    _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),
    _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),
    _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),
    _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),
    _F(6.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
    _F(6.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
    _F(6.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
    _F(6.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
    _F(6.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
    _F(6.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
    _F(6.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
    _F(6.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
    _F(6.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
    _F(6.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
    _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),
    _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),
    _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),
    _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),
    _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),
    _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),
    _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),
    _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),
    _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),
    _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),
    _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),
    _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),
    _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),
    _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),
    _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),
    _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),
    _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),
    _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),
    _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),
    _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),
    _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),
    _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),
    _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),
    _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),
    _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),
    _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),
    _F(5.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
    _F(5.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
    _F(5.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
    _F(5.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
    _F(5.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
    _F(5.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
    _F(5.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
    _F(5.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
    _F(5.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
    _F(5.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
    _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),
    _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),
    _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),
    _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),
    _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),
    _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),
    _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),
    _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),
    _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),
    _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),
    _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),
    _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),
    _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),
    _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),
    _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),
    _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),
    _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),
    _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),
    _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),
    _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),
    _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),
    _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),
    _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),
    _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),
    _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),
    _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),
    _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),
    _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),
    _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),
    _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),
    _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),
    _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),
    _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),
    _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),
    _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),
    _F(4.)/_F(16.), _F(10.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.),  _F(9.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
    _F(4.)/_F(16.),  _F(8.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
    _F(4.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
    _F(4.)/_F(16.),  _F(7.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
    _F(4.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
    _F(4.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
    _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),
    _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),
    _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),
    _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),
    _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),
    _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),
    _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),
    _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),
    _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),
    _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),
    _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),
    _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),
    _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),
    _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),
    _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),
    _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),
    _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),
    _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),
    _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),
    _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),
    _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),
    _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),
    _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),
    _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),
    _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),
    _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),
    _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),
    _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),
    _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),
    _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.), _F(10.)/_F(16.),  _F(1.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(2.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(3.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(4.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(5.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(6.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(7.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(8.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(9.)/_F(16.),
    _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.), _F(10.)/_F(16.),
    _F(3.)/_F(16.), _F(11.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.), _F(10.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.), _F(10.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.),  _F(9.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.),  _F(9.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
    _F(3.)/_F(16.),  _F(8.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.),  _F(8.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
    _F(3.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
    _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),
    _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),
    _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),
    _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),
    _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),
    _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),
    _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),
    _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),
    _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),
    _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),
    _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),
    _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),
    _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),
    _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),
    _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),
    _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),
    _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),
    _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),
    _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),
    _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),
    _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),
    _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),
    _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),
    _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),
    _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.), _F(10.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(9.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(3.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(4.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(5.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(6.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(7.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(8.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(9.)/_F(16.),
    _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.), _F(10.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.), _F(11.)/_F(16.),  _F(1.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.), _F(10.)/_F(16.),  _F(2.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(3.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(4.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(5.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(6.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(7.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(8.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(9.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.), _F(10.)/_F(16.),
    _F(3.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.), _F(11.)/_F(16.),
    _F(2.)/_F(16.), _F(12.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.), _F(11.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.), _F(11.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.), _F(10.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.), _F(10.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.), _F(10.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.),  _F(9.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.),  _F(9.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.),  _F(9.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
    _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),
    _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),
    _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),
    _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),
    _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),
    _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),
    _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),
    _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),
    _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),
    _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),
    _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),
    _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),
    _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),
    _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),
    _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),
    _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),
    _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),
    _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),
    _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),
    _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.), _F(10.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(9.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(8.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(4.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(5.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(6.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(7.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(8.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),  _F(9.)/_F(16.),
    _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.), _F(10.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.), _F(11.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.), _F(10.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(9.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(4.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(5.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(6.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(7.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(8.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),  _F(9.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.), _F(10.)/_F(16.),
    _F(2.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.), _F(11.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.), _F(12.)/_F(16.),  _F(1.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.), _F(11.)/_F(16.),  _F(2.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.), _F(10.)/_F(16.),  _F(3.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(4.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(5.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(6.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(7.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(8.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(9.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.), _F(10.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.), _F(11.)/_F(16.),
    _F(2.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.), _F(12.)/_F(16.),
    _F(1.)/_F(16.), _F(13.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.), _F(12.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.), _F(12.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.), _F(11.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.), _F(11.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.), _F(11.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.), _F(10.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.), _F(10.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.), _F(10.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.), _F(10.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),
    _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),
    _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),
    _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),
    _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),
    _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),
    _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(8.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(7.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),
    _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),
    _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),
    _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),
    _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(9.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(8.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(7.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(6.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),  _F(5.)/_F(16.),
    _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),
    _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),
    _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),
    _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.), _F(10.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(9.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(8.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(7.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(6.)/_F(16.),  _F(5.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(5.)/_F(16.),  _F(6.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(4.)/_F(16.),  _F(7.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(3.)/_F(16.),  _F(8.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(2.)/_F(16.),  _F(9.)/_F(16.),
    _F(1.)/_F(16.),  _F(4.)/_F(16.),  _F(1.)/_F(16.), _F(10.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.), _F(11.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.), _F(10.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(9.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(8.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(7.)/_F(16.),  _F(5.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(6.)/_F(16.),  _F(6.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(5.)/_F(16.),  _F(7.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(4.)/_F(16.),  _F(8.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(3.)/_F(16.),  _F(9.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(2.)/_F(16.), _F(10.)/_F(16.),
    _F(1.)/_F(16.),  _F(3.)/_F(16.),  _F(1.)/_F(16.), _F(11.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.), _F(12.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.), _F(11.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.), _F(10.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(9.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(8.)/_F(16.),  _F(5.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(7.)/_F(16.),  _F(6.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(6.)/_F(16.),  _F(7.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(5.)/_F(16.),  _F(8.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(4.)/_F(16.),  _F(9.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(3.)/_F(16.), _F(10.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(2.)/_F(16.), _F(11.)/_F(16.),
    _F(1.)/_F(16.),  _F(2.)/_F(16.),  _F(1.)/_F(16.), _F(12.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.), _F(13.)/_F(16.),  _F(1.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.), _F(12.)/_F(16.),  _F(2.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.), _F(11.)/_F(16.),  _F(3.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.), _F(10.)/_F(16.),  _F(4.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(9.)/_F(16.),  _F(5.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(8.)/_F(16.),  _F(6.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(7.)/_F(16.),  _F(7.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(6.)/_F(16.),  _F(8.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(5.)/_F(16.),  _F(9.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(4.)/_F(16.), _F(10.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(3.)/_F(16.), _F(11.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(2.)/_F(16.), _F(12.)/_F(16.),
    _F(1.)/_F(16.),  _F(1.)/_F(16.),  _F(1.)/_F(16.), _F(13.)/_F(16.)
};
DOF_TYPE DOF_P16_ = {
    DofReserved, "P16", P16_points, NULL, DOF_DG15, NULL, NULL,
    phgDofInterpC2FGeneric, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
    Pn_bas, Pn_grad, Pn_map, FE_H1, TRUE, TRUE, FALSE, -1,
    969, 16, 0, 0, 1, 1, 15, 105, 455
};

/*============ End of structs generated by 'utils/lagrange_gen' ============*/

#endif	/* defined(GENERATE_STRUCTS) */
