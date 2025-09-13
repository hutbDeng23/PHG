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

/* $Id: periodic.c,v 1.43 2020/12/02 01:38:53 zlb Exp $ */

/* Functions for handling periodic meshes */

#include "phg.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* workaround for a possible bug with icc <=11.1 + OpenMP on xeon, or with
 * icc 13.1.3 on mic */
#if defined(__ICC) /*&& (defined(__MIC__) || __ICC <= 1110)*/
INT
_phg_global_vertex_p(GRID *g, INT no)
{
#if 0
    /* With --enable-shared the next two lines are unnecessary.
     * With --disable-shared and without the next two lines, the programs crash
     * in updrn_callback(), at invocation of GlobalVertexP
     * Note: need to use at least two processes to show the problem. */
#if USE_OMP
# pragma omp threadprivate(dummy)
#endif
    static void *dummy;
    dummy = g->period;
#endif

    volatile INT i;	/* Note: without 'volatile' programs segfault on mic */
    i = g->period == NULL ? GlobalVertex(g, no) : g->period->L2Gmap_vert[no];
    return i;
}
#endif	/* defined(__ICC) ... */

void
phgSetPeriodicity(GRID *g, int flags)
{
    if (g->nelem_global != 0)
	phgError(1, "%s must be called before importing the grid.\n", __func__);

    if (flags == 0)
	return;

    if (g->period == NULL)
	g->period = phgCalloc(1, sizeof(*g->period));
    g->period->flags = flags;
}

void
phgSetPeriodicDirections(GRID *g, const FLOAT dirs[])
{
    if (g->nelem_global != 0)
	phgError(1, "%s must be called before importing the grid.\n", __func__);

    if (dirs == NULL)
	return;

    if (g->period == NULL)
	g->period = phgCalloc(1, sizeof(*g->period));
    phgFree(g->period->directions);
    g->period->directions = phgAlloc(Dim * Dim * sizeof(dirs[0]));
    memcpy(g->period->directions, dirs, Dim * Dim * sizeof(dirs[0]));
}

/* The functions below are PHG's internal functions */

void
phgPeriodFree(GRID *g)
{
    if (g->period == NULL)
	return;

    phgFree(g->period->directions);
    phgFree(g->period->L2Gmap_vert);
    phgFree(g->period->ordering);
    phgFree(g->period);
    g->period = NULL;
}

typedef struct {
    ELEMENT	*e;		/* the element */
    INT		c0, c1;		/* scaled coordinates of the barycentre */
    int		face;		/* face no. */
} FACE_T;

#define IsEqual(a, b)	(((a) - (b)) / 2 == 0)

/* variables for building face-face map */
static INT list_n, list_alloc;
static FACE_T *list;

static FLOAT (*tran)[Dim] = NULL, (*itran)[Dim] = NULL, (*bbox)[Dim] = NULL;
static COORD *verts = NULL;	/* coordinates in the reference orthobrick */

static int
face_comp(const void *p0, const void *p1)
{
    INT i;
    FACE_T *f0 = list + *(const INT *)p0, *f1 = list + *(const INT *)p1;

    if (IsEqual(f0->c0, f1->c0))
	return IsEqual(f0->c1, f1->c1) ? 
			0 : (f0->c1 - f1->c1 > (INT)0 ? 1 : -1);

    return (i = f0->c0 - f1->c0) > 0 ? 1 : (i < 0 ? -1 : 0);
}

static INT clist[6];	/* 2D coordinates of the three vertices of a face */

static int
vertex_comp(const void *p0, const void *p1)
{
    INT i;
    INT *a = clist + *(const int *)p0 * 2, *b = clist + *(const int *)p1 * 2;

    if (IsEqual(*a, *b)) {
	++a;
	++b;
	return IsEqual(*a, *b) ? 0 : (*a - *b > 0 ? 1 : -1);
    }

    return (i = *a - *b) > 0 ? 1 : (i < 0 ? -1 : 0);
}

static void
map_vertices(GRID *g, int dir, INT (*vmap)[Dim])
/* builds face-face maps, returns number of faces as return value.
 * 'dir' should be 0 for x, 1 for y, or 2 for z.
 *
 * The size of the array 'vmap', if non nil,  should be g->nvert with entries
 * initialized to -1 before the first call to this function.
 * On return the leading numbers in vmap[i] which are not -1 are vertices
 * which are periodic to vertex i. */
{
    static int scaling = 10000000;	/* for coordinates compare/sort */
    ELEMENT *e;
    COORD *c[NVert];
    int ii, v0, v1, v2, a0, a1;
    FACE_T *pl0, *pl1;
    FLOAT min0, d0, min1, d1, lower, upper;
    const char *dirname;
    FLOAT eps = 2. / (FLOAT)scaling;
    INT i;

    switch (dir) {
	case 0:
	    a0 = 1;
	    a1 = 2;
	    dirname = "x";
	    break;
	case 1:
	    a0 = 0;
	    a1 = 2;
	    dirname = "y";
	    break;
	case 2:
	    a0 = 0;
	    a1 = 1;
	    dirname = "z";
	    break;
	default:	/* make gcc happy */
	    a0 = a1 = 0;
	    dirname = NULL;
    }

    lower = bbox[0][dir];
    upper = bbox[1][dir];
    min0 = 3. * bbox[0][a0];
    min1 = 3. * bbox[0][a1];
    d0   = scaling / (3. * bbox[1][a0] - min0);
    d1   = scaling / (3. * bbox[1][a1] - min1);

    list = NULL;
    list_n = list_alloc = 0;
    for (e = g->roots; e < g->roots + g->nroot; e++) {
	/* loop on faces of 'e' */
	c[0] = verts + e->verts[0];
	c[1] = verts + e->verts[1];
	c[2] = verts + e->verts[2];
	c[3] = verts + e->verts[3];
	for (ii = 0; ii < NFace; ii++) {
	    v0 = GetFaceVertex(ii, 0);
	    v1 = GetFaceVertex(ii, 1);
	    v2 = GetFaceVertex(ii, 2);
	    if (Fabs((*c[v0])[dir] - lower) > eps &&
		Fabs((*c[v0])[dir] - upper) > eps)
		continue;
	    if (Fabs((*c[v0])[dir] - (*c[v1])[dir]) >= eps ||
		Fabs((*c[v0])[dir] - (*c[v2])[dir]) >= eps ||
		Fabs((*c[v1])[dir] - (*c[v2])[dir]) >= eps)
		continue;
	    /* face on a periodic boundary, save it to the list */
	    if (list_n >= list_alloc) {
		list = phgRealloc_(list, (list_alloc + 2048) * sizeof(*list),
					 list_alloc * sizeof(*list));
		list_alloc += 2048;
	    }
	    pl0 = list + list_n++;
	    pl0->e = e;
	    pl0->c0 = RoundINT(((*c[v0])[a0] + (*c[v1])[a0] + (*c[v2])[a0]
				 - min0) * d0);
	    pl0->c1 = RoundINT(((*c[v0])[a1] + (*c[v1])[a1] + (*c[v2])[a1]
				 - min1) * d1);
	    pl0->face = ii;
	}
    }

    if (list_n % 2)
	phgError(1, "%s:%d, the mesh is non conforming on periodic boundaries "
		    "in the %s direction.\n", __FILE__, __LINE__, dirname);

    if (list_n > 0) {
	INT *index = phgAlloc(list_n * sizeof(*index));
	int u[3], v[3], a[3], b[3];

	for (i = 0; i < list_n; i++)
	    index[i] = i;
	qsort(index, list_n, sizeof(*index), face_comp);

	min0 = min0 / 3.;
	min1 = min1 / 3.;
	d0 = d0 * 3.;
	d1 = d1 * 3.;

	for (i = 0; i < list_n; i += 2) {
	    pl0 = list + index[i];
	    pl1 = list + index[i + 1];

	    /* check that the faces are properly paired up */
	    if (!IsEqual(pl0->c0, pl1->c0) || !IsEqual(pl0->c1, pl1->c1))
		phgError(1, "%s:%d, the mesh is non conforming on periodic "
			    "boundaries in the %s direction.\n",
			    __FILE__, __LINE__, dirname);

	    /* find matching vertex pairs */

	    ii = pl0->face;
	    u[0] = GetFaceVertex(ii, 0);
	    u[1] = GetFaceVertex(ii, 1);
	    u[2] = GetFaceVertex(ii, 2);
	    e = pl0->e;
	    e->bound_type[ii] &= ~BDRY_MASK;	/* clear boundary masks */
	    for (ii = 0; ii < 3; ii++) {
		clist[2 * ii + 0] = RoundINT((verts[e->verts[u[ii]]][a0]
					      - min0) * d0);
		clist[2 * ii + 1] = RoundINT((verts[e->verts[u[ii]]][a1]
					      - min1) * d1);
		a[ii] = ii;
	    }
	    /* sort the three vertices */
	    qsort(a, 3, sizeof(a[0]), vertex_comp);

	    ii = pl1->face;
	    v[0] = GetFaceVertex(ii, 0);
	    v[1] = GetFaceVertex(ii, 1);
	    v[2] = GetFaceVertex(ii, 2);
	    e = pl1->e;
	    e->bound_type[ii] &= ~BDRY_MASK;	/* clear boundary masks */
	    for (ii = 0; ii < 3; ii++) {
		clist[2 * ii + 0] = RoundINT((verts[e->verts[v[ii]]][a0]
					      - min0) * d0);
		clist[2 * ii + 1] = RoundINT((verts[e->verts[v[ii]]][a1]
					      - min1) * d1);
		b[ii] = ii;
	    }
	    /* sort the three vertices */
	    qsort(b, 3, sizeof(b[0]), vertex_comp);

	    /* the two faces should be on opposite sides */
	    if (Fabs(verts[pl0->e->verts[u[a[0]]]][dir] -
		     verts[pl1->e->verts[v[b[0]]]][dir]) < 2. * eps)
		phgError(1, "wrong period in %s direction\n", dirname);

	    for (ii = 0; ii < 3; ii++) {
		INT i0, i1;

		/* vertex u[a[ii]] in pl0 matches v[b[ii]] in pl1 */
		i0 = pl0->e->verts[u[a[ii]]];
		i1 = pl1->e->verts[v[b[ii]]];

		if (Fabs(verts[i0][a0] - verts[i1][a0]) * d0 > 0.5 ||
		    Fabs(verts[i0][a1] - verts[i1][a1]) * d1 > 0.5)
		    phgError(1, "%s:%d, unmatched face!\n", __FILE__, __LINE__);

		/* save 'i1' to 'vmap[i0]' */
		if (vmap[i0][0] != -1 && vmap[i0][0] != i1) {
		    if (vmap[i0][1] != -1 && vmap[i0][1] != i1) {
			if (vmap[i0][2] == -1)
			    vmap[i0][2] = i1;
			else if (vmap[i0][2] != i1)
			    phgError(1, "%s:%d, unexpected error.\n",
					__FILE__, __LINE__);
		    }
		    else {
			if (vmap[i0][1] == -1)
			    vmap[i0][1] = i1;
		    }
		}
		else {
		    if (vmap[i0][0] == -1)
			vmap[i0][0] = i1;
		}

		/* save 'i0' to 'vmap[i1]' */
		if (vmap[i1][0] != -1 && vmap[i1][0] != i0) {
		    if (vmap[i1][1] != -1 && vmap[i1][1] != i0) {
			if (vmap[i1][2] == -1)
			    vmap[i1][2] = i0;
			else if (vmap[i1][2] != i0)
			    phgError(1, "%s:%d, unexpected error.\n",
					__FILE__, __LINE__);
		    }
		    else {
			if (vmap[i1][1] == -1)
			    vmap[i1][1] = i0;
		    }
		}
		else {
		    if (vmap[i1][0] == -1)
			vmap[i1][0] = i0;
		}
	    }	/* for ii */
	}	/* for i */

	phgFree(index);
    }

    phgFree(list);

    return;
}

/* vertex-vertex map (one vertex can be mapped to atmost three vertices) */
static INT (*vmap)[Dim];

static void
find_vertices(INT i, INT *list, int *n)
{
    INT j;
    int k;

    for (k = 0; k < 3 && (j = vmap[i][k]) >= 0; k++) {
	vmap[i][k] = -2;
	if (vmap[j][0] == -2)
	    break;
	assert(*n < 8);
	list[(*n)++] = j;
	find_vertices(j, list, n);
    }
}

void
phgPeriodInit(GRID *g)
/* Note: should be called right after importing the grid */
{
    INT i, list[8], nvert;
    int ii, n;

    if (g->period == NULL || g->nvert == 0)
	return;

    assert(g->nprocs == 1);	/* currently only for unpartitioned grid */

    if (g->period->directions == NULL) {
	bbox = g->bbox;
	verts = g->verts;
    }
    else {
	/* Note: itran is not used, it's saved for possible future use. */
	phgFree(itran);
	itran = phgAlloc(Dim * sizeof(*itran));
	memcpy(itran, g->period->directions, Dim * sizeof(*itran));
	phgFree(tran);
	tran = phgCalloc(Dim, sizeof(*tran));
	tran[0][0] = tran[1][1] = tran[2][2] = 1.0;
	if (!phgSolverDenseSolver(Dim, Dim, g->period->directions,
					(FLOAT *)tran))
	    phgError(1, "%s: periodic directions are dependent!\n", __func__);
	phgFree(g->period->directions);
	g->period->directions = NULL;

	/* compute bbox of the reference orthobrick */
	phgFree(bbox);
	bbox = phgAlloc(2 * sizeof(*bbox));
	phgFree(verts);
	verts = phgAlloc(g->nvert * sizeof(*verts));
	for (ii = 0; ii < Dim; ii++) {
	    verts[0][ii] = tran[0][ii] * g->verts[0][0] +
			   tran[1][ii] * g->verts[0][1] +
			   tran[2][ii] * g->verts[0][2];
	    bbox[0][ii] = bbox[1][ii] = verts[0][ii];
	}
	for (i = 1; i < g->nvert; i++) {
	    for (ii = 0; ii < Dim; ii++) {
		FLOAT a;
		verts[i][ii] = a = tran[0][ii] * g->verts[i][0] +
				   tran[1][ii] * g->verts[i][1] +
				   tran[2][ii] * g->verts[i][2];
		if (bbox[0][ii] > a)
		    bbox[0][ii] = a;
		else if (bbox[1][ii] < a)
		    bbox[1][ii] = a;
	    }
	}
    }

    vmap = phgAlloc(g->nvert * sizeof(*vmap));
    for (i = 0; i < g->nvert; i++) {
	for (ii = 0; ii < 3; ii++)
	    vmap[i][ii] = -1;
    }

    if (g->period->flags & X_MASK)
	map_vertices(g, 0, vmap);

    if (g->period->flags & Y_MASK)
	map_vertices(g, 1, vmap);

    if (g->period->flags & Z_MASK)
	map_vertices(g, 2, vmap);

    nvert = 0;
    for (i = 0; i < g->nvert; i++) {
#if 0
if (vmap[i][0] >= 0)
printf("vertex %"dFMT": %"dFMT" %"dFMT" %"dFMT"\n", i, vmap[i][0], vmap[i][1], vmap[i][2]);
#endif
	if (vmap[i][0] == -2)
	    continue;
	list[0] = i;
	n = 1;
	/* recursively find all vertices (at most 7) mapped to vertex i */
	find_vertices(i, list, &n);
#if 0
printf("n=%d:",n); for (ii=0; ii<n; ii++) printf(" %"dFMT, list[ii]); printf("\n");
#endif

	/* map the vertices to 'nvert' */
	for (ii = 0; ii < n; ii++)
	    vmap[list[ii]][1] = nvert;

	nvert++;
    }

#if 0
for (i = 0; i < g->nvert; i++) printf("vertex %"dFMT" ==> %"dFMT"\n", i, vmap[i][1]);
printf("%"dFMT" vertices (originally %"dFMT")\n", nvert, g->nvert);
#endif

    /* build L2Gmap_vert */
    assert(g->period->L2Gmap_vert == NULL);
    g->period->nvert_global = nvert;
    g->period->L2Gmap_vert = phgAlloc(g->nvert *
					sizeof(*g->period->L2Gmap_vert));
    for (i = 0; i < g->nvert; i++)
	g->period->L2Gmap_vert[i] = vmap[i][1];

    phgFree(vmap);

    phgInfo(2, "%s: nvert = %"dFMT", period->nvert_global = %"dFMT"\n",
		__func__, g->nvert, g->period->nvert_global);

    if (tran != NULL) {
	phgFree(tran);
	phgFree(itran);
	phgFree(bbox);
	phgFree(verts);
    }

    tran = NULL;
    itran = NULL;
    bbox = NULL;
    verts = NULL;
}
