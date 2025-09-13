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

/* $Id: intergrid.c,v 1.46 2022/07/08 01:01:05 zlb Exp $ */

/* This file implements inter-grid operations (e.g., DOF transfer) */

#include "phg.h"

#include <string.h>
#include <math.h>

typedef struct BLOCKS_ {	/* struct for a chain of buffers */
    void		*buffer;
    struct BLOCKS_	*next;		/* next block */
    size_t		item_size;	/* item (entry) size */
    size_t		block_size;	/* block size */
} BLOCKS;

/*------------------- utility functions for managing boxes ------------------*/

/* Note: ../mesh-files/dolfin-gear.mesh fails when setting EPS to 1e-12 */
#define EPS 1e-9

/* Function for factorizing an integer: n = n0 * n1 * ...
 * such that n0 * weights[0] + n1 * weights[1] + ...  is minimized. */

#define MAXFACTORS	16

static int n_save, ndim, *np, *npwork;
static double *weights, best_sum;
static int numfac;
static int factors[MAXFACTORS], powers[MAXFACTORS];

/* Note: the functions with a "stack" argument are recursive */

static void next_dimension(int dimen);		/* forward declaration */

static void
next_factor(int facno, int prod, int dimen)
{
    int i, p, f;

    if (facno >= numfac) {
	npwork[dimen++] = prod;
	next_dimension(dimen);
	return;
    }

    p = powers[facno];
    f = factors[facno];
    for (i = 0; i <= p; i++, powers[facno]--) {
	next_factor(facno + 1, prod, dimen);
	prod *= f;
    }
    powers[facno] = p;

    return;
}

static int
compare(double sum, int *npwork, double best_sum, int *np)
/* returns true if (sum,npwork) < (best_sum,np) */
{
    int i;

    if (sum < best_sum)
	return 1;
    if (sum > best_sum)
	return 0;
    for (i = 0; i < ndim; i++) {
	if (npwork[i] < np[i])
	    return 1;
	if (npwork[i] > np[i])
	    return 0;
    }
    return 0;
}

static void
next_dimension(int dimen)
{
    static int i, prod;
    static double sum;

    if (dimen >= ndim - 1) {
	prod = 1;
	for (i = 0; i < dimen; i++)
	    prod *= npwork[i];
	npwork[dimen] = n_save / prod;

	/* evaluate this case */
	sum = 0.0;
	for (i = 0; i < ndim; i++)
	    sum += (npwork[i] - 1) * (weights == NULL ? 1.0 : weights[i]);
	if (compare(sum, npwork, best_sum, np)) {
	    best_sum = sum;
	    memcpy(np, npwork, ndim * sizeof(*np));
	}
    }
    else {
	next_factor(0, 1, dimen);
    }

    return;
}

static int
factorize(int n, int *np0, double *weights0, int ndim0)
/* find np0[0], np0[1],..., np0[*ndim] such that:
	np0[0]*np0[1]*...*np0[*ndim] == n
	(np0[0]-1)*weights0[0] +
	(np0[1]-1)*weights0[1] +
	... +
	(np0[ndim]-1)*weights0[ndim] = minimum,
  returns !0 on errors */
{
    int id, iq, p;

    n_save = n;
    np = np0;
    weights = weights0;
    ndim = ndim0;

    if (n <= 1) {
	for (id = 0; id < ndim; id++)
	    np[id] = n;
	return 0;
    }

    np[0] = n;
    best_sum = (np[0] - 1) * (weights == NULL ? 1.0 : weights[0]);
    for (id = 1; id < ndim; id++)
	np[id] = 1;

    /* First, factorize n into prime factors */
    p = numfac = 0;
    id = 2;
    while (id <= n) {
	iq = n / id;
	if (n == iq * id) {
	    if (!p)
		p = 1;
	    else
		p++;
	    n = iq;
	}
	else {
	    if (p) {
		if (numfac == MAXFACTORS)
		    return 2;
		factors[numfac] = id;
		powers[numfac++] = p;
		p = 0;
	    }
	    id = id == 2 ? 3 : id + 2;
	}
    }
    if (p) {
	if (numfac == MAXFACTORS)
	    return 2;
	factors[numfac] = id;
	powers[numfac++] = p;
    }

    /* Loop on all possible combinations to find the minimum */
    if ((npwork = phgAlloc(ndim * sizeof(*npwork))) == NULL)
	return 3;
    next_dimension(0);
    phgFree(npwork);

    return 0;
}

/*--------------------------------------------------------------------------*/

/* struct for recording an array of boxes and lists of elements intersecting
 * with each box.
 *
 * Note: This information may be stored in GRID in the member 'void *box_info'
 *	 (which only needs to be updated when the grid is redistributed or
 *	 flattened, i.e., when nroot is changed).
 *
 *	 The struct BOX_INFO definition can be kept private to this file,
 *	 the function phgInterGridDofEval() can be used to free GRID.box_info
 *	 with 'u == NULL' and 'n == GRID id'.
 */
typedef struct {
    FLOAT	dx, dy, dz;	/* size of the boxes */
    INT		**lists;	/* nx*ny*nz lists of elements list */
    INT		nx, ny, nz;	/* number of boxes in x, y and z direction */
} BOX_INFO;

static void
dump_medit(GRID *g, ELEMENT *e, FLOAT a0, FLOAT b0, FLOAT c0,
			        FLOAT a1, FLOAT b1, FLOAT c1)
{
    static int no = 0;
    char fn[32];
    FILE *f;
    COORD *c;
    int jj;

    sprintf(fn, "tmp%03d.mesh", no++);
    printf("create file \"%s\".\n", fn);
    f = fopen(fn, "w+t");
    fprintf(f, "MeshVersionFormatted 1\nDimension 3\nVertices 12\n");
    for (jj = 0; jj < NVert; jj++) {
	c = g->verts + e->verts[jj];
	fprintf(f, "%0.16le %0.16le %0.16le 0\n", 
			(double)(*c)[0], (double)(*c)[1], (double)(*c)[2]);
    }
    fprintf(f, "%0.16le %0.16le %0.16le 0\n", (double)a0,(double)b0,(double)c0);
    fprintf(f, "%0.16le %0.16le %0.16le 0\n", (double)a0,(double)b0,(double)c1);
    fprintf(f, "%0.16le %0.16le %0.16le 0\n", (double)a0,(double)b1,(double)c1);
    fprintf(f, "%0.16le %0.16le %0.16le 0\n", (double)a0,(double)b1,(double)c0);
    fprintf(f, "%0.16le %0.16le %0.16le 0\n", (double)a1,(double)b0,(double)c0);
    fprintf(f, "%0.16le %0.16le %0.16le 0\n", (double)a1,(double)b0,(double)c1);
    fprintf(f, "%0.16le %0.16le %0.16le 0\n", (double)a1,(double)b1,(double)c1);
    fprintf(f, "%0.16le %0.16le %0.16le 0\n", (double)a1,(double)b1,(double)c0);
    fprintf(f, "Tetrahedra 1\n1 2 3 4 1\nHexahedra 1\n5 6 7 8 9 10 11 12\n");
    fprintf(f, "End\n");
    fclose(f);
}

#if 0	/* unfinished */
static BOOLEAN
check_2D_triangle(FLOAT a0, FLOAT b0, FLOAT a1, FLOAT b1,
		  FLOAT x0, FLOAT y0, FLOAT x1, FLOAT y1, FLOAT x2, FLOAT y2)
/* returns TRUE if the rectangle [a0,a1]x[b0,b1] intersects with the triangle
 * (x0,y0)-(x1,y1)-(x2,y2) */
{
    FLOAT u0, v0, u1, v1;	/* the bounding box */

    u0 = u1 = x0;
    if (u0 > x1)
	u0 = x1;
    else if (u1 < x1)
	u1 = x1;
    if (u0 > x2)
	u0 = x2;
    else if (u1 < x2)
	u1 = x2;

    v0 = v1 = y0;
    if (v0 > y1)
	v0 = y1;
    else if (v1 < y1)
	v1 = y1;
    if (v0 > y2)
	v0 = y2;
    else if (v1 < y2)
	v1 = y2;

    if (a0 > u1 || a1 < u0 || b0 > v1 || b1 < v0)
	return FALSE;

    if (a0 < u0)
	a0 = u0;
    if (a1 > u1)
	a1 = u1;
    if (b0 < v0)
	b0 = v0;
    if (b1 > v1)
	b1 = v1;

    /* compute the 2D interval barycentric coordinates */

    return TRUE;
}

static BOOLEAN
check_2D(FLOAT a0, FLOAT b0, FLOAT a1, FLOAT b1, FLOAT x0, FLOAT y0,
	 FLOAT x1, FLOAT y1, FLOAT x2, FLOAT y2, FLOAT x3, FLOAT y3)
/* returns TRUE if the rectangle [a0,a1]x[b0,b1] intersects with the triangles
 * (x0,y0)-(x1,y1)-(x2,y2) and (x0,y0)-(x1,y1)-(x3,y3) */
{
    if (check_2D_triangle(a0, b0, a1, b1, x0, y0, x1, y1, x2, y2))
	return TRUE;
    if (check_2D_triangle(a0, b0, a1, b1, x0, y0, x1, y1, x3, y3))
	return TRUE;
    if (check_2D_triangle(a0, b0, a1, b1, x0, y0, x2, y2, x3, y3))
	return TRUE;

printf("------------\n");
    return FALSE;
}
#endif	/* 0 */

static BOOLEAN
check_box(GRID *g, ELEMENT *e, BOX_INFO *box_info, int ix, int iy, int iz,
		FLOAT x0, FLOAT y0, FLOAT z0, FLOAT x1, FLOAT y1, FLOAT z1)
/* returns TRUE if the box (ix,iy,iz) intersects with the element 'e' */
{
    FLOAT a0, b0, c0, a1, b1, c1;
    FLOAT lam0[Dim + 1], lam1[Dim + 1];
    const FLOAT *J;
    COORD *c[Dim + 1];
    int ii, jj;

    /* Note: performing the following tests may reduce number of elements/box,
     *	     but the tests are expensive and may require longer time than the
     *	     time saved with fewer elements per box, so the tests may be
     *	     skipped by uncommenting the next line. */
    return TRUE;

    a0 = g->bbox[0][0] + ix * box_info->dx;
    b0 = g->bbox[0][1] + iy * box_info->dy;
    c0 = g->bbox[0][2] + iz * box_info->dz;
    a1 = a0 + box_info->dx;
    b1 = b0 + box_info->dy;
    c1 = c0 + box_info->dz;

    if (a0 < x0)
	a0 = x0;
    if (a1 > x1)
	a1 = x1;

    if (b0 < y0)
	b0 = y0;
    if (b1 > y1)
	b1 = y1;

    if (c0 < z0)
	c0 = z0;
    if (c1 > z1)
	c1 = z1;

    /* check for vertices lying in (a0,a1)x(b0,b1)x(c0,c1) */
    for (ii = 0; ii < NVert; ii++) {
	c[ii] = g->verts + e->verts[ii];
	if ((*c[ii])[0] > a0 - EPS && (*c[ii])[0] < a1 + EPS &&
	    (*c[ii])[1] > b0 - EPS && (*c[ii])[1] < b1 + EPS &&
	    (*c[ii])[2] > c0 - EPS && (*c[ii])[2] < c1 + EPS)
	    return TRUE;
    }

    /* check interval barycentric coordinates of ([a0,a1], [b0,b1], [c0,c1]) */
    J = phgGeomGetJacobian(g, e);
    for (ii = jj = 0; ii < Dim + 1; ii++) {
	lam0[ii] = J[0] * (J[0] > 0. ? a0 : a1) +
		   J[1] * (J[1] > 0. ? b0 : b1) +
		   J[2] * (J[2] > 0. ? c0 : c1) + J[3];
	lam1[ii] = J[0] * (J[0] < 0. ? a0 : a1) +
		   J[1] * (J[1] < 0. ? b0 : b1) +
		   J[2] * (J[2] < 0. ? c0 : c1) + J[3];
	J += Dim + 1;
	if (lam0[ii] <= -EPS) {
	    if (lam1[ii] <= -EPS) {
#if 0
	    if (ix == 6 && iy == 2 && iz == 7)
		dump_medit(g, e, a0, b0, c0, a1, b1, c1);
#else
	    Unused(dump_medit);
#endif
		return FALSE;
	    }
	    continue;
	}
	if (lam1[ii] > -EPS)
	    jj++;
    }
    if (jj == Dim + 1)
	return TRUE;

    /* If the 2D projections of the box and the tetrahedron are disjoint in
     * one of the 3 space directions, then return FALSE, else return TRUE */

#if 0
    if (!check_2D(a0, b0, a1, b1,
			(*c[0])[0], (*c[0])[1], (*c[1])[0], (*c[1])[1],
			(*c[2])[0], (*c[2])[1], (*c[3])[0], (*c[3])[1]))
	return FALSE;

    if (!check_2D(a0, c0, a1, c1,
			(*c[0])[0], (*c[0])[2], (*c[1])[0], (*c[1])[2],
			(*c[2])[0], (*c[2])[2], (*c[3])[0], (*c[3])[2]))
	return FALSE;

    if (!check_2D(b0, c0, b1, c1,
			(*c[0])[1], (*c[0])[2], (*c[1])[1], (*c[1])[2],
			(*c[2])[1], (*c[2])[2], (*c[3])[1], (*c[3])[2]))
	return FALSE;
#endif	/* 0 */

    return TRUE;
}

static BOX_INFO *
init_boxes(GRID *g)
{
    BOX_INFO *box_info;
    ELEMENT *e;
    double time, weights[Dim];
    int np[Dim];
    int i, k, m;
    INT *nalloc, nx, ny, nz;
    FLOAT dx, dy, dz;

    /* number of boxes = g->nroot / 8 */
    m = g->nroot / 8;
    if (m <= 1) {
	return NULL;
    }

    time = phgGetTime(NULL);

    box_info = phgAlloc(sizeof(*box_info));
    weights[0] = fabs(1.0 / (double)(g->bbox[1][0] - g->bbox[0][0]));
    weights[1] = fabs(1.0 / (double)(g->bbox[1][1] - g->bbox[0][1]));
    weights[2] = fabs(1.0 / (double)(g->bbox[1][2] - g->bbox[0][2]));
#if 0
    if (m >= 60)
	m = ((m + 29) / 30) * 30;
    else
	m = ((m + 5) / 6) * 6;
#else	/* 0 */
    /* round m to power of 2 */
    k = 1;
    /*m = (m + 1) / 2;*/
    while (k < m)
	k += k;
    m = k;
#endif	/* 0 */
    if (factorize(m, np, weights, Dim) != 0)
	phgWarning("failure in factorize(), shouldn't happen\n");
    nx = box_info->nx = np[0];
    ny = box_info->ny = np[1];
    nz = box_info->nz = np[2];
    assert(nx > 0 && ny > 0 && nz > 0 && nx * ny * nz == m);

#if DEBUG
    if (phgVerbosity >= 1) {
	phgInfo(-1, "%s: BBox = (%0.3f,%0.3f)x(%0.3f,%0.3f)x(%0.3f,%0.3f)\n",
			__func__,
			(double)g->bbox[0][0], (double)g->bbox[1][0],
			(double)g->bbox[0][1], (double)g->bbox[1][1],
			(double)g->bbox[0][2], (double)g->bbox[1][2]);
	phgInfo(-1, "%s: %"dFMT" root elements, %d boxes, "
		    "nx*ny*nz = %"dFMT"x%"dFMT"x%"dFMT"\n",
			__func__, g->nroot, m, nx, ny, nz);
    }
#endif	/* DEBUG */

    dx = box_info->dx = (g->bbox[1][0] - g->bbox[0][0]) / nx;
    dy = box_info->dy = (g->bbox[1][1] - g->bbox[0][1]) / ny;
    dz = box_info->dz = (g->bbox[1][2] - g->bbox[0][2]) / nz;

    box_info->lists = phgCalloc(m, sizeof(*box_info->lists));
    nalloc = phgAlloc(m * sizeof(*nalloc));

    for (i = 0; i < g->nroot; i++) {
	/* determine lists of boxes with which the tetrahedron intersects.
	 * cf.: Helmut Ratschek and Jon Rokne,
	 * 	Test for intersection between box and tetrahedron
	 * 	International Journal of Computer Mathematics, Volume 65,
	 * 	Issue 3 & 4 1997, pages 191-204
	 * 	File: bibliography/IntersectionBoxTetrahedron.pdf
	 */
	FLOAT x0, y0, z0, x1, y1, z1;
	COORD *c;
	int ii, ix, iy, iz, ix0, iy0, iz0;

	e = g->roots + i;

	/* compute the isothetic box hull of e (x0,x1)x(y0,y1)x(z0,z1) */

	c = g->verts + e->verts[0];
	x0 = x1 = (*c)[0];
	y0 = y1 = (*c)[1];
	z0 = z1 = (*c)[2];
	for (ii = 1; ii < NVert; ii++) {
	    c = g->verts + e->verts[ii];
	    if (x0 > (*c)[0])
		x0 = (*c)[0];
	    else if (x1 < (*c)[0])
		x1 = (*c)[0];
	    if (y0 > (*c)[1])
		y0 = (*c)[1];
	    else if (y1 < (*c)[1])
		y1 = (*c)[1];
	    if (z0 > (*c)[2])
		z0 = (*c)[2];
	    else if (z1 < (*c)[2])
		z1 = (*c)[2];
	}

	/* shift the box */
//printf("element %"dFMT"\n", i);
//printf("   x0 y0 z0 = %g %g %g, x1 y1 z1 = %g %g %g\n", x0, y0, z0, x1, y1, z1);
	x0 -= g->bbox[0][0]; x1 -= g->bbox[0][0];
	y0 -= g->bbox[0][1]; y1 -= g->bbox[0][1];
	z0 -= g->bbox[0][2]; z1 -= g->bbox[0][2];

	/* loop over boxes intersecting the isothetic box hull of e */
	ix0 = (int)(x0 / dx);
	iy0 = (int)(y0 / dy);
	iz0 = (int)(z0 / dz);
//printf("ix0 iy0 iz0 = %d %d %d, dx dy dz = %g %g %g\n", ix0, iy0, iz0, dx, dy, dz);
	for (ix = ix0; ix < nx && ix * dx <= x1; ix++) {
	    for (iy = iy0; iy < ny && iy * dy <= y1; iy++) {
		for (iz = iz0; iz < nz && iz * dz <= z1; iz++) {
//printf("   checking box %d %d %d\n", ix, iy, iz);
		    if (!check_box(g, e, box_info, ix, iy, iz,
					x0 + g->bbox[0][0],
					y0 + g->bbox[0][1],
					z0 + g->bbox[0][2],
					x1 + g->bbox[0][0],
					y1 + g->bbox[0][1],
					z1 + g->bbox[0][2]))
			continue;
//printf("   overlap\n");
		    /* add e to boxes->lists[k] */
		    k = (iz * ny + iy) * nx + ix;
		    if (box_info->lists[k] == NULL) {
			nalloc[k] = 8;
			box_info->lists[k] = phgAlloc(nalloc[k] *
						sizeof(*box_info->lists[k]));
			box_info->lists[k][0] = 1;
		    }
		    else if (++box_info->lists[k][0] >= nalloc[k]) {
			box_info->lists[k] = phgRealloc_(box_info->lists[k],
				(nalloc[k] + 8) * sizeof(*box_info->lists[k]),
				nalloc[k] * sizeof(*box_info->lists[k]));
			nalloc[k] += 8;
		    }
		    box_info->lists[k][box_info->lists[k][0]] = i;
		}
	    }
	}
    }

    phgFree(nalloc);

#if DEBUG
    if (phgVerbosity >= 1) {
	double total = 0.;
	for (k = i = 0; i < m; i++) {
	    if (box_info->lists[i] == NULL)
		continue;
	    if (k < box_info->lists[i][0])
		k = box_info->lists[i][0];
	    total += box_info->lists[i][0] / (double)m;
	}
	phgInfo(-1, "%s: average/max elements per box = %0.2lf/%d, "
		    "time = %0.4lf\n", __func__,
		    total, k, phgGetTime(NULL) - time);
    }
#else	/* DEBUG */
    Unused(time);
#endif	/* DEBUG */

    return box_info;
}

static void
free_boxes(BOX_INFO *box_info)
{
    int m;

    if (box_info == NULL)
	return;

    m = box_info->nx * box_info->ny * box_info->nz;
    while (--m >= 0)
	phgFree(box_info->lists[m]);
    phgFree(box_info->lists);
    phgFree(box_info);
}

/*----------------- end utility functions for managing boxes ----------------*/

static CHAR *flags = NULL;
static size_t dim;
static ELEMENT *e_found;
static FLOAT lambda_found[Dim + 1];

static BOOLEAN
find_recursive(ELEMENT *e, const FLOAT lam[])
{
    static int v[NVert];
    ELEMENT *e0, *e1;
    FLOAT a, lam1[Dim + 1];

    e0 = e->children[0];
    e1 = e->children[1];

    if (e0 == NULL && e1 == NULL) {
	e_found = e;
	memcpy(lambda_found, lam, sizeof(lambda_found));
	return TRUE;
    }

    if (e0 != NULL) {
	phgMapP2C(e, v, NULL, 0);
	if ((a = lam[0] - lam[1]) > -EPS) {
	    lam1[3]     = 2. * lam[1];
	    lam1[v[0]]  = a;
	    lam1[v[2]]  = lam[2];
	    lam1[v[3]]  = lam[3];
	    if (find_recursive(e0, lam1) == TRUE)
		return TRUE;
	}
    }

    if (e1 != NULL) {
	phgMapP2C(e, v, NULL, 1);
	if ((a = lam[1] - lam[0]) > -EPS) {
	    lam1[3]     = 2. * lam[0];
	    lam1[v[1]]  = a;
	    lam1[v[2]]  = lam[2];
	    lam1[v[3]]  = lam[3];
	    if (find_recursive(e1, lam1) == TRUE)
		return TRUE;
	}
    }

    return FALSE;
}

static const FLOAT *
find_element(GRID *g, ELEMENT **ep, BOX_INFO *box_info, FLOAT *jacobian[],
					FLOAT x, FLOAT y, FLOAT z)
/* finds an element containing the given point, returns a pointer to the
 * barycentric coordinates (NULL if the point is not in the submesh),
 * the element found is returned in '*ep'. */
{
    static INT seed = 0;
    INT i, j, k, nelem, *index;
    ELEMENT *e;
    const FLOAT *J;
    FLOAT a, lambda[Dim + 1], *lam;

    *ep = NULL;

    /* first search in the root elements */
    if (box_info == NULL) {
	nelem = g->nroot;
	index = NULL;
    }
    else {
	int ix, iy, iz;
	ix = (int)((x - g->bbox[0][0]) / box_info->dx);
	if (ix >= box_info->nx)
	    ix = box_info->nx - 1;
	iy = (int)((y - g->bbox[0][1]) / box_info->dy);
	if (iy >= box_info->ny)
	    iy = box_info->ny - 1;
	iz = (int)((z - g->bbox[0][2]) / box_info->dz);
	if (iz >= box_info->nz)
	    iz = box_info->nz - 1;
	i = (iz * box_info->ny + iy) * box_info->nx + ix;
	if ((index = box_info->lists[i]) == NULL)
	    nelem = 0;
	else
	    nelem = *(index++);
    }

    if (nelem == 0)
	return NULL;

    seed %= nelem;		/* make sure seed is in the range */
    for (j = 0; j < nelem; j++) {
	/* Note: j even ==> check seed + j / 2, j odd ==> check seed - j / 2 */
	i = (j + 1) / 2;
	i = (seed + ((j & 2) ? nelem - i : i)) % nelem;
	k = (index == NULL ? i : index[i]);
	e = g->roots + k;
	if (IsLeaf(e)) {
	    phgGeomXYZ2Lambda(g, e, x, y, z, lam = lambda);
	}
	else {
	    if (jacobian[k] == NULL) {
		jacobian[k] = phgAlloc((Dim + 1) * (Dim + 1) * sizeof(*J));
		memcpy(jacobian[k], phgGeomGetJacobian(g, e),
				(Dim + 1) * (Dim + 1) * sizeof(*J));
	    }
	    J = jacobian[k];
	    lambda[0] = J[0] * x + J[1] * y + J[2] * z + J[3];
	    J += 4;
	    lambda[1] = J[0] * x + J[1] * y + J[2] * z + J[3];
	    J += 4;
	    lambda[2] = J[0] * x + J[1] * y + J[2] * z + J[3];
	    J += 4;
	    lambda[3] = J[0] * x + J[1] * y + J[2] * z + J[3];
	    lam = lambda;
	}
	a = lam[0];
	if (a > lam[1])
	    a = lam[1];
	if (a > lam[2])
	    a = lam[2];
	if (a > lam[3])
	    a = lam[3];
	if (a > -EPS) {
	    if (find_recursive(e, lam) == TRUE) {	/* found */
		seed = i;
		*ep = e_found;
		return lambda_found;
	    }
	    if (a > EPS) {
		seed = i;
		/* the point can't be within other elements */
		return NULL;
	    }
	}
    }

    return NULL;
}

#if USE_MPI
static void MPIAPI
mpi_op_func(void *invec, void *inoutvec, int *len, MPI_Datatype *type)
{
    int i;
    FLOAT *p = invec, *q = inoutvec;

    Unused(type);

    for (i = 0; i < *len; i++, p += dim + 1, q += dim + 1) {
	if (*(CHAR *)(p + dim) <= *(CHAR *)(q + dim))
	    continue;
	memcpy(q, p, dim * sizeof(*q) + 1);
    }
}
#endif	/* USE_MPI */

static void
dof_eval(DOF *u, INT n, BLOCKS *cblocks, BLOCKS *vblocks, int flag)
/* Note: cblocks and vblocks may overlap to save memory, as long as they
 * have the same item_size */	 
{
    INT i, m;
    GRID *g = u->g;
    ELEMENT *e;
    char *pnts, *pnts_end, *vals, *vals_end;
    BLOCKS *cb, *vb;
    const FLOAT *lam;
    BOX_INFO *box_info;
    FLOAT **jacobian;	/* for saving Jacobian on non-leaf root elements */

    if (n == 0 && g->nprocs == 1)
	return;

    dim = DofDim(u);
    box_info = init_boxes(g);
    jacobian = phgCalloc(g->nroot, sizeof(*jacobian));

    /* evaluate on local submesh */
    flags = phgAlloc(n * sizeof(*flags));
    cb = cblocks;
    pnts = cb->buffer;
    pnts_end = pnts + cb->block_size * cb->item_size;
    vb = vblocks;
    vals = vb->buffer;
    vals_end = vals + vb->block_size * vb->item_size;
    for (m = i = 0; i < n; ) {	/* loop on points */
	lam = find_element(g, &e, box_info, jacobian, (*(COORD *)pnts)[0],
				(*(COORD *)pnts)[1], (*(COORD *)pnts)[2]);
	if (lam != NULL) {
	    flags[i] = 1;
	    phgDofEval(u, e, lam, (void *)vals);
	}
	else {
	    flags[i] = 0;
	    if (g->nprocs == 1)
		phgError(1, "%s (%s:%d): point %"dFMT", (%lg,%lg,%lg) is "
			    "outside the domain.\n",
			    __func__, __FILE__, __LINE__, i,
			    (*(COORD *)pnts)[0], (*(COORD *)pnts)[1],
			    (*(COORD *)pnts)[2]);
	    m++;
	}

	if (++i >= n)
	    break;

	if ((pnts += cb->item_size) >= pnts_end) {
	    cb = cb->next;
	    if (cb == NULL)
		phgError(1, "%s: error in cblocks chain.\n", __func__);
	    pnts = cb->buffer;
	    pnts_end = pnts + cb->block_size * cb->item_size;
	}

	if ((vals += vb->item_size) >= vals_end) {
	    vb = vb->next;
	    if (vb == NULL)
		phgError(1, "%s: error in vblocks chain.\n", __func__);
	    vals = vb->buffer;
	    vals_end = vals + vb->block_size * vb->item_size;
	}
    }

#if USE_MPI
    if (g->nprocs == 1) {
	phgFree(flags);
	for (i = 0; i < g->nroot; i++)
	    phgFree(jacobian[i]);
	phgFree(jacobian);
	free_boxes(box_info);
	return;
    }

    if (flag >= 0) {
	/* same points on all processes, collect results using MPI_Reduce */
	MPI_Datatype type;

	MPI_Op op;
	FLOAT *values0, *values1, *p;

	values0 = phgAlloc(2 * n * (dim + 1) * sizeof(*values0));
	values1 = values0 + n * (dim + 1);

	p = values0;
	vb = vblocks;
	vals = vb->buffer;
	vals_end = vals + vb->block_size * vb->item_size;
	for (i = 0; i < n; i++) {
	    *(CHAR *)(p + dim) = flags[i];
	    memcpy(p, vals, dim * sizeof(*p));
	    p += dim + 1;
	    if ((vals += vb->item_size) >= vals_end) {
		vb = vb->next;
		if (vb == NULL)
		    break;
		vals = vb->buffer;
		vals_end = vals + vb->block_size * vb->item_size;
	    }
	}

	MPI_Op_create(/* func = */mpi_op_func, /* commute = */0, &op);
	{
	    /* MPI_Datatype for {dim*FLOAT, CHAR}, extent = (dim+1)*FLOAT */
#if 1	/* MPI-2 */
	    MPI_Datatype type0;
	    int counts[] = {dim, 1};
	    MPI_Aint indices[] = {0, dim * sizeof(FLOAT)};
	    MPI_Datatype types[] = {PHG_MPI_FLOAT, PHG_MPI_CHAR};
	    MPI_Type_create_struct(2, counts, indices, types, &type0);
	    MPI_Type_create_resized(type0, 0, (dim+1)*sizeof(FLOAT), &type);
	    MPI_Type_free(&type0);
#else	/* MPI-1 */
	    int counts[] = {dim, 1, 1};
	    MPI_Aint indices[] = {0, dim*sizeof(FLOAT), (dim+1)*sizeof(FLOAT)};
	    MPI_Datatype types[] = {PHG_MPI_FLOAT, PHG_MPI_CHAR, MPI_UB};
	    MPI_Type_struct(3, counts, indices, types, &type);
#endif
	}
	MPI_Type_commit(&type);
	/* FIXME: on LSSC3 the next MPI_Reduce crashes if mpi-2level is used. */
	if (flag == 1)
	    MPI_Allreduce(values0, values1, n, type, op, g->comm);
	else
	    MPI_Reduce(values0, values1, n, type, op, 0, g->comm);
	MPI_Type_free(&type);
	MPI_Op_free(&op);

	if (g->rank == 0 || flag == 1) {
	    p = values1;
	    vb = vblocks;
	    vals = vb->buffer;
	    vals_end = vals + vb->block_size * vb->item_size;
	    for (i = 0; i < n; i++) {
		flags[i] = *(CHAR *)(p + dim);
		if (flags[i] == 0)
		    phgError(1, "%s (%s:%d): point %"dFMT" is outside the "
				"domain.\n", __func__, __FILE__, __LINE__, i);
		memcpy(vals, p, dim * sizeof(*p));
		p += dim + 1;
		if ((vals += vb->item_size) >= vals_end) {
		    vb = vb->next;
		    if (vb == NULL)
			break;
		    vals = vb->buffer;
		    vals_end = vals + vb->block_size * vb->item_size;
		}
	    }
	}

	phgFree(values0);
    }
    else {	/* if (flag >= 0) */
	/* different points on different processes: circular algorithm */
	INT *counts;
	void *tmp, *sbuf, *rbuf, *buf0;
	FLOAT *pf;
	CHAR *pc;
	MPI_Datatype type;
	MPI_Status status;
	int ii, jj, kk, prev, next, dim0;
	size_t size;

	assert(flag == -1);
	counts = phgAlloc(g->nprocs * sizeof(*counts));
	MPI_Allgather(&m, 1, PHG_MPI_INT, counts, 1, PHG_MPI_INT, g->comm);
	size = counts[0];
	for (ii = 1; ii < g->nprocs; ii++) {
	    if (size < counts[ii])
		size = counts[ii];
	}

	if (size == 0) {
	    phgFree(counts);
	    phgFree(flags);
	    for (i = 0; i < g->nroot; i++)
		phgFree(jacobian[i]);
	    phgFree(jacobian);
	    free_boxes(box_info);
	    return;
	}

	dim0 = (Dim >= dim ? Dim : dim);
	size *= dim0 * sizeof(*pf) + sizeof(*pc);
	size = (size + sizeof(*pf) - 1) / sizeof(*pf);
	pf = phgAlloc(2 * size * sizeof(*pf));
	buf0 = sbuf = pf;
	rbuf = pf + size;

	if (m > 0) {
	    /* Collect missing values (note m = number of missing values) */
	    pf = sbuf;
	    pc = (CHAR *)(pf + m * dim0);
	    cb = cblocks;
	    pnts = cb->buffer;
	    pnts_end = pnts + cb->block_size * cb->item_size;
	    for (i = 0; i < n; ) {
		if (flags[i] == 0) {
		    memcpy(pf, pnts, Dim * sizeof(FLOAT));
		    pf += dim0;
		    *(pc++) = 0;
		}
		if (++i >= n)
		    break;
		if ((pnts += cb->item_size) >= pnts_end) {
		    cb = cb->next;
		    assert(cb != NULL);
		    pnts = cb->buffer;
		    pnts_end = pnts + cb->block_size * cb->item_size;
		}
	    }
	    assert(pf - (FLOAT *)sbuf == m * dim0);
	}

	/* circulate the lists (TODO: overlap communication/computation?) */
	MPI_Type_contiguous(dim0 * sizeof(*pf) + 1, MPI_BYTE, &type);
	MPI_Type_commit(&type);
	prev = (g->rank + g->nprocs - 1) % g->nprocs;
	next = (g->rank + 1) % g->nprocs;
	jj = g->rank;
	for (ii = 0; ; ) {
	    kk = (jj + g->nprocs - 1) % g->nprocs;
	    MPI_Sendrecv(sbuf, m,          type, next, 222,
			 rbuf, counts[kk], type, prev, 222,
			 g->comm, &status);
	    tmp = sbuf;
	    sbuf = rbuf;
	    rbuf = tmp;
	    m = counts[jj = kk];
	    if (++ii >= g->nprocs)
		break;
	    /* process received list pointed by buf */
	    if (m == 0)
		continue;
	    pf = sbuf;
	    pc = (CHAR *)(pf + m * dim0);
	    for (i = 0; i < m; i++, pc++, pf += dim0) {
		if (*pc != 0)
		    continue;
		lam = find_element(g, &e, box_info, jacobian,
						pf[0], pf[1], pf[2]);
		if (lam != NULL) {
		    *pc = 1;
		    phgDofEval(u, e, lam, pf);
		}
	    }
	}
	MPI_Type_free(&type);

	/* pick up values from the returned list */
	if (m > 0) {
	    pf = sbuf;
	    pc = (CHAR *)(pf + m * dim0);
	    vb = vblocks;
	    vals = vb->buffer;
	    vals_end = vals + vb->block_size * vb->item_size;
	    for (i = 0; i < n; i++) {
		if (flags[i] == 0) {
		    if (*(pc++) == 0)
			phgError(1, "%s (%s:%d): point %"dFMT" is outside the "
				 "domain.\n", __func__, __FILE__, __LINE__, i);
		    memcpy(vals, pf, dim * sizeof(FLOAT));
		    pf += dim0;
		}
		if ((vals += vb->item_size) >= vals_end) {
		    vb = vb->next;
		    if (vb == NULL)
			break;
		    vals = vb->buffer;
		    vals_end = vals + vb->block_size * vb->item_size;
		}
	    }
	}

	phgFree(buf0);
	phgFree(counts);
    }		/* if (flag >= 0) */
#endif	/* USE_MPI */

    phgFree(flags);
    for (i = 0; i < g->nroot; i++)
	phgFree(jacobian[i]);
    phgFree(jacobian);
    free_boxes(box_info);

    return;
}

void
phgInterGridDofEval(DOF *u, INT n, COORD *pts, FLOAT *values, int flag)
/* evaluates 'dof' on the list of 'n' points 'pts' and returns results in
 * 'values', which points to a buffer FLOAT[n][DofDim(dof)].
 *
 * 'flag' == 0:	the list of points is identical on all processes and
 *		the results are only needed by process 0.
 *
 * 'flag' == 1:	the list of points is identical on all processes and
 *		the results are needed by all processes.
 *
 * 'flag' ==-1:	the list of points is different on each process and
 *		the results are distributed in the same way.
 *
 * Note:
 *
 *   1.	It is assumed that all points belong to the domain (an error is
 *	reported if points outside the domain are found)
 *
 *   2.	'values' may overlap with 'pts' when DofDim(u) <= Dim
 */
{
    BLOCKS *cblocks = phgAlloc(2 * sizeof(*cblocks));
    BLOCKS *vblocks = cblocks + 1;

    cblocks->buffer = pts;
    cblocks->next = NULL;
    cblocks->item_size = sizeof(*pts);
    cblocks->block_size = n;

    vblocks->buffer = values;
    vblocks->next = NULL;
    vblocks->item_size = DofDim(u) * sizeof(*values);
    vblocks->block_size = n;

    dof_eval(u, n, cblocks, vblocks, flag);

    phgFree(cblocks);
}

static int pass;
static BLOCKS *blocks;
static char *ptr, *ptr_end;
static INT total_points;
/* func_ptr is the optional coordinate transformation function */
static void (*func_ptr)(FLOAT x, FLOAT y, FLOAT z,
			FLOAT *xout, FLOAT *yout, FLOAT *zout) = NULL;

static void
copy_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
/* dummy function */
{
    if (pass == 0) {
	/* save the evaluation point */
	if ((ptr += blocks->item_size) >= ptr_end) {
	    /* allocate next block */
	    blocks->next = phgAlloc(sizeof(*blocks));
	    blocks->next->item_size = blocks->item_size;
	    blocks->next->block_size = blocks->block_size;
	    blocks = blocks->next;
	    blocks->buffer = phgAlloc(blocks->item_size * blocks->block_size);
	    blocks->next = NULL;
	    ptr = blocks->buffer;
	    ptr_end = ptr + blocks->block_size * blocks->item_size;
	}
	if (func_ptr == NULL) {
	    (*(COORD *)ptr)[0] = x;
	    (*(COORD *)ptr)[1] = y;
	    (*(COORD *)ptr)[2] = z;
	}
	else {
	    func_ptr(x, y, z, &(*(COORD *)ptr)[0], &(*(COORD *)ptr)[1],
			&(*(COORD *)ptr)[2]);
	}
    }
    else {
	/* return evaluated values */
	memcpy(values, ptr, dim * sizeof(*values));
	if ((ptr += blocks->item_size) >= ptr_end) {
	    blocks = blocks->next;
	    if (blocks != NULL) {
		ptr = blocks->buffer;
		ptr_end = ptr + blocks->block_size * blocks->item_size;
	    }
	}
    }
    total_points++;
}

static void
free_blocks(BLOCKS *cb)
{
    if (cb->next != NULL)
	free_blocks(cb->next);
    phgFree(cb->buffer);
    phgFree(cb);
}

void
phgInterGridDofTransfer_(DOF *src, DOF *dest,
			 void (*xyzfunc)(FLOAT x, FLOAT y, FLOAT z,
					 FLOAT *xout, FLOAT *yout, FLOAT *zout))
/* Warning: this function is slow and not scalable */
{
    BLOCKS *blocks_root;

    assert(!SpecialDofType(dest->type));
    assert(src->g != dest->g);	/* shouldn't be used if same grid */
    dim = DofDim(dest);
    assert(dim == DofDim(src));

    func_ptr = xyzfunc;

    blocks = blocks_root = phgAlloc(sizeof(*blocks_root));
    blocks->item_size = (Dim >= dim ? Dim : dim) * sizeof(FLOAT);
    blocks->block_size = 16384;
    blocks->buffer = phgAlloc(blocks->item_size * blocks->block_size);
    blocks->next = NULL;

    /* collect points on which src values are needed */
    pass = 0;
    total_points = 0;
    blocks = blocks_root;
    ptr = blocks->buffer;
    ptr_end = ptr + blocks->block_size * blocks->item_size;
    ptr -= blocks->item_size;
    CheckD_order(dest);
    phgDofSetDataByFunction(dest, copy_func);

    /* compute src values on collected points (overlap values/coordinates) */
    dof_eval(src, total_points, blocks_root, blocks_root, -1);

    /* compute dest using evaluated src values */
    pass = 1;
    total_points = 0;
    blocks = blocks_root;
    ptr = blocks->buffer;
    ptr_end = ptr + blocks->block_size * blocks->item_size;
    phgDofSetDataByFunction(dest, copy_func);

    free_blocks(blocks_root);
}


INT 
phgInterGridPt2EMap(GRID *g, INT n, COORD *pts, ELEMENT **elist)
/* find the elements containing 'n' points 'pts',
 * and returns results in  'elist', elist[i]=NULL 
 * if the i-th point is not in the submesh. 
 * the number of the points in the submesh will be retured
 */
{
    INT i, m;
    ELEMENT *e;
    const FLOAT *lam;
    BOX_INFO *box_info;
    FLOAT **jacobian;	/* for saving Jacobian on non-leaf root elements */

    if (n == 0)
	return 0;
    
    box_info = init_boxes(g);
    jacobian = phgCalloc(g->nroot, sizeof(*jacobian));

    /* evaluate on local submesh */
    for (m = i = 0; i < n; i++ ) {	/* loop on points */
	lam = find_element(g, &e, box_info, jacobian, pts[i][0],
				pts[i][1], pts[i][2]);
	if (lam != NULL) {
	    elist[i] =  e;
	    m++;
	}
	else {
	    elist[i] = NULL;
#if 0
	    if (g->nprocs == 1)
		phgWarning("%s (%s:%d): point %"dFMT", (%lg,%lg,%lg) is "
			    "outside the domain.\n",
			    __func__, __FILE__, __LINE__, i,
			    (*(COORD *)pts)[0], (*(COORD *)pts)[1],
			    (*(COORD *)pts)[2]);
#endif
	}
    }

    for (i = 0; i < g->nroot; i++)
	phgFree(jacobian[i]);
    phgFree(jacobian);
    free_boxes(box_info);

    return m;
}
