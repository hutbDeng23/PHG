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

/* $Id: refine.c,v 1.102 2022/03/24 06:31:27 zlb Exp $ */

/*--------------------------------------------------------------------------
 * Original bisection diagram by Kossaczky:
 *
 *			v0v1v2v3->v4=(v0+v1)/2		<--- type DIAG(3)
 *		       /		      \
 *		v0v2v3v4->v5=(v0+v2)/2	   v1v3v2v4	<--- type FACE(2)
 *		/		      \
 *	v0v3v4v5->v6=(v0+v3)/2	   v2v3v4v5		<--- type EDGE(1)
 *	/		  \
 * v0v4v5v6		v3v4v5v6			<--- type DIAG(3)
 *
 * 2 new types:
 *
 *			v0v1v2v3->v4=(v0+v1)/2		<--- type OPPOSITE(4)
 *		       /		\
 *		  v3v2v0v4	      v3v2v1v4		<--- type FACE(2)
 *
 *			v0v1v2v3->v4=(v0+v1)/2		<--- type MIXED(5)
 *		       /		\
 *		  v0v2v3v4	      v3v2v1v4		<--- type FACE(2)
 *--------------------------------------------------------------------------*/

#include "phg.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>	/* INT_MAX */

#define REFINE_RECURSION_LIMIT	256

/* Note: CompVertex is only used for equality tests */
#define CompVertex(v0, v1) \
	comp_vertex((v0) < nvert_save ? GlobalVertexP(g_, v0) : \
					(v0) + g_->nvert_global, \
		    (v1) < nvert_save ? GlobalVertexP(g_, v1) : \
					(v1) + g_->nvert_global)

static void refine_path(GRID *g, ELEMENT *e);

static GRID *g_;

typedef struct {
    INT	a, b;			/* the two parent vertices */
    INT	index;			/* local index of the vertex */
} VERTEX_t;

/* saved value of old g->nvert */
static size_t nvert_save;

/* The verts[] array is used by comp_vertex() during refinement,
 * and by period_update_xxxx()/update_L2Gmap_vert() after refinement.
 *
 * During refinement, new vertices are inserted to verts[] by new_vertex(),
 *	If g->period == NULL then verts[] is NULL and not used.
 *	If g->period != NULL then verts[] contains all new vertices.
 */
static VERTEX_t *verts = NULL;
static size_t verts_count = 0, verts_allocated = 0;

static int
comp_index(const void *p0, const void *p1)
{
    INT i;
    return (i = ((VERTEX_t *)p0)->index - ((VERTEX_t *)p1)->index) > 0 ?
		1 : (i < 0 ? -1 : 0);
}

static int
comp_vertex(INT v0, INT v1)
/* compares global indices of two vertices. v0 and v1 are global indices
   for old vertices, and local index + g->nvert_global for new vertices.

   Since the global indices of new vertices are not yet assigned,
   we recursively look up the verts[] array. */
{
    static int ii;
    static INT i;
    VERTEX_t *p0, *p1;

    if (v0 == v1)
	return 0;

    if (v0 < g_->nvert_global || v1 < g_->nvert_global)
	return (i = v0 - v1) > 0 ? 1 : (i < 0 ? -1 : 0);

    /* both v0 and v1 are new vertices, search the verts[] table */

    if (verts == NULL) {
	/* simply compare v0 and v1 (only used for equality test) */
	return (i = v0 - v1) > 0 ? 1 : (i < 0 ? -1 : 0);
    }

    v0 -= g_->nvert_global + nvert_save;
    p0 = verts + v0;
    assert(p0->index == v0 + nvert_save);

    v1 -= g_->nvert_global + nvert_save;
    p1 = verts + v1;
    assert(p1->index == v1 + nvert_save);

    if ((ii = comp_vertex(p0->a, p1->a)) != 0)
	return ii;

    return comp_vertex(p0->b, p1->b);
}

static BYTE
opposite_vertex(GRID *g, ELEMENT *e, BYTE v, ELEMENT *e_op)
/* returns the vertex e_op opposite to the vertex v of e.
 *
 * Note: phgOppositeVertex can't be used during refinement since
 *	 it depends on g->period->L2Gmap_vert[] which is not available yet */
{
    INT i, u0 = 0, u1 = 0, u2 = 0;

#if DEBUG
    if (e_op == NULL)
	phgError(1, "unexpected error in %s (possibly caused by inconsistent "
		 "boundary types in the input mesh file).\n", __func__);
#endif	/* DEBUG */

    switch (v) {
	case 0: u0 = e->verts[1]; u1 = e->verts[2]; u2 = e->verts[3]; break;
	case 1: u0 = e->verts[0]; u1 = e->verts[2]; u2 = e->verts[3]; break;
	case 2: u0 = e->verts[0]; u1 = e->verts[1]; u2 = e->verts[3]; break;
	case 3: u0 = e->verts[0]; u1 = e->verts[1]; u2 = e->verts[2]; break;
    }

    i = e_op->verts[0];
    if (CompVertex(i, u0) && CompVertex(i, u1) && CompVertex(i, u2))
	return 0;

    i = e_op->verts[1];
    if (CompVertex(i, u0) && CompVertex(i, u1) && CompVertex(i, u2))
	return 1;

    i = e_op->verts[2];
    if (CompVertex(i, u0) && CompVertex(i, u1) && CompVertex(i, u2))
	return 2;

    i = e_op->verts[3];
    if (CompVertex(i, u0) && CompVertex(i, u1) && CompVertex(i, u2))
	return 3;

    return 4;	/* make gcc happy */
}

static VERTEX_t *verts_gathered = NULL;
static INT *block_start = NULL;

static int
comp_composite_vertex0(const VERTEX_t *g0, const VERTEX_t *g1)
{
    /* Note: since this function is recursive, the following vars are declared
     * static to save stack space */
    static int ii;
    static INT i;
    static VERTEX_t vtx0, vtx1;
    static VERTEX_t *base0 = NULL, *base1 = NULL;
    VERTEX_t *u0, *u1, *v0, *v1;

    /* vertex a */

    if (g0->a < g_->nvert_global || g1->a < g_->nvert_global) {
	/* If one or both of g0->a and g1->a are old vertices,
	   then simply compare their values */
	if ((i = g0->a - g1->a) != 0)
	    return i > 0 ? 1 : -1;
	if (g0->b < g_->nvert_global || g1->b < g_->nvert_global)
	    return (i = g0->b - g1->b) > 0 ? 1 : (i < 0 ? -1 : 0);
	if (block_start != NULL) {
	    base0 = verts_gathered + block_start[g0 - verts_gathered];
	    base1 = verts_gathered + block_start[g1 - verts_gathered];
	}
	else {
	    base0 = verts_gathered;
	    base1 = verts_gathered;
	}
	u0 = u1 = NULL;
    }
    else {
	/* since g0->a is a new vertex, a corresponding entry must exist
	   in the same block as g0 in 'verts_gathered' and before g0 with
	   index == g0->a.

	   Similarly since g1->a is a new vertex, a corresponding entry must
	   exist in the same block as g1 in 'verts_gathered' and before g1 with
	   index == g1->a. */
	if (block_start != NULL) {
	    base0 = verts_gathered + block_start[g0 - verts_gathered];
	    base1 = verts_gathered + block_start[g1 - verts_gathered];
	}
	else {
	    base0 = verts_gathered;
	    base1 = verts_gathered;
	}
	vtx0.index = g0->a - g_->nvert_global;
	vtx1.index = g1->a - g_->nvert_global;
	u0 = bsearch(&vtx0, base0, g0 - base0, sizeof(*g0), comp_index);
	assert(u0 != NULL);
	u1 = bsearch(&vtx1, base1, g1 - base1, sizeof(*g1), comp_index);
	assert(u1 != NULL);
    }
    
    /* since g0->b is a new vertex, a corresponding entry must exist
       in the same block as g0 in 'verts_gathered' and before g0 with
       index == g0->b

       Similarly since g1->b is a new vertex, a corresponding entry must
       exist in the same block as g1 in 'verts_gathered' and before g1
       index == g1->b. */
    vtx0.index = g0->b - g_->nvert_global;
    vtx1.index = g1->b - g_->nvert_global;
    v0 = bsearch(&vtx0, base0, g0 - base0, sizeof(*g0), comp_index);
    assert(v0 != NULL);
    v1 = bsearch(&vtx1, base1, g1 - base1, sizeof(*g1), comp_index);
    assert(v1 != NULL);

    if (u0 == NULL)
	return comp_composite_vertex0(v0, v1);

    if ((ii = comp_composite_vertex0(u0, u1)) != 0)
	return ii;

    return comp_composite_vertex0(v0, v1);
}

static int
comp_composite_vertex(const void *p0, const void *p1)
/* compare two VERTEX_t structs, p0 and p1 are indices to verts_gathered */
{
    return comp_composite_vertex0(verts_gathered + *(const INT *)p0,
				  verts_gathered + *(const INT *)p1);
}

#if 0
static void
eval_vertex(GRID *g, VERTEX_t *v, FLOAT *x, FLOAT *y, FLOAT *z)
/* compute coordinates of new vertex i (for debugging) */
{
    VERTEX_t *u, vtx;
    FLOAT x0, y0, z0, x1, y1, z1;

    if (v->a < g->nvert_global) {
	x0 = g->verts[v->a_local][0];
	y0 = g->verts[v->a_local][1];
	z0 = g->verts[v->a_local][2];
    }
    else {
	vtx.index = v->a - g->nvert_global;
	u = bsearch(&vtx, verts, v - verts, sizeof(*v), comp_index);
	assert(u != NULL);
	eval_vertex(g, u, &x0, &y0, &z0);
    }

    if (v->b < g->nvert_global) {
	x1 = g->verts[v->b_local][0];
	y1 = g->verts[v->b_local][1];
	z1 = g->verts[v->b_local][2];
    }
    else {
	vtx.index = v->b - g->nvert_global;
	u = bsearch(&vtx, verts, v - verts, sizeof(*v), comp_index);
	assert(u != NULL);
	eval_vertex(g, u, &x1, &y1, &z1);
    }

    *x = (x0 + x1) * 0.5;
    *y = (y0 + y1) * 0.5;
    *z = (z0 + z1) * 0.5;
}

static void
check_vertices(GRID *g)
/* this function is for debugging purpose */
{
    INT i;
    FLOAT x, y, z, error, e;

    error = 0;
    for (i = 0; i < verts_count; i++) {
        eval_vertex(g, verts + i, &x, &y, &z);
	x -= g->verts[verts[i].index][0];
	y -= g->verts[verts[i].index][1];
	z -= g->verts[verts[i].index][2];
        e = x*x + y*y + z*z;
	if (e > error)
	    error = e;
    }
    phgInfo(2, "nvert_save=%"dFMT", g->nvert=%"dFMT", verts_count=%"dFMT"\n",
		(INT)nvert_save, g->nvert, (INT)verts_count);
    phgInfo(2, "error check: %lf\n", (double)error);
}
#endif	/* 0 */

static void
period_update_private(GRID *g)
/* updates g->period->L2Gmap_vert for private vertices. */
{
    INT i;
    for (i = nvert_save; i < g->nvert; i++)
	if (g->period->L2Gmap_vert[i] == -1)
	    g->period->L2Gmap_vert[i] = g->period->nvert_global++;
}

#if 0
/* This version of period_update_shared() uses an algorithm similar to the
 * one used in phgUpdateBoundaryTypes(). It currently does not work because
 * the way for determining to which process to send an entry is unknown
 * (the current method based on root_vertex() doesn't work) */
 *
static VERTEX_t *
root_vertex(VERTEX_t *v)
/* this function recursively looks up the smallest root vertex of a
 * composite vertex, using vertices array verts[]. */
{
    INT a = v->a;

    if (a < g_->nvert_global)
	return v;

    v = verts + (a -= g_->nvert_global) - nvert_save;
    assert(v->index == a);

    return root_vertex(v);
}

static void
period_update_shared(GRID *g)
/* updates g->period->L2Gmap_vert and g->period->nvert_global for shared
 * vertices. The list of new vertices are in the buffer pointed by 'verts'
 *
 * TODO: compare performance of the algorithm used here and the one used in
 *	 update_L2Gmap_vert(), and use the faster one for both functions */
{
    INT i, j, m, n, count, offset, *ordering, *isbuf, *irbuf, *counts;
    VERTEX_t *sbuf, *rbuf;
    int rank, *scnts, *sdsps, *rcnts, *rdsps;
    MPI_Datatype type;

    scnts = phgAlloc(4 * g->nprocs * sizeof(*scnts));
    sdsps = scnts + g->nprocs;
    rcnts = sdsps + g->nprocs;
    rdsps = rcnts + g->nprocs;

    /* count # of vertices to send to each process */
    memset(scnts, 0, g->nprocs * sizeof(*scnts));
    for (i = 0; i < verts_count; i++) {
	assert(g->period->L2Gmap_vert[verts[i].index] == -1);
	rank = root_vertex(verts + i)->a % g->nprocs;
	scnts[rank]++;
    }

    /* exchange # of vertices */
    MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, g->comm);

    /* set up sdsps and rdsps */
    m = n = 0;
    for (rank = 0; rank < g->nprocs; rank++) {
	sdsps[rank] = m;
	rdsps[rank] = n;
	m += scnts[rank];
	n += rcnts[rank];
    }

    /* allocate send/recv buffers */
    sbuf = phgAlloc(m * sizeof(*sbuf));
    rbuf = phgAlloc(n * sizeof(*rbuf));
    isbuf = phgAlloc((m + n + g->nprocs) * sizeof(*irbuf));
    irbuf = isbuf + m;
    counts = irbuf + n;

    /* collect vertices to send */
    for (i = 0; i < verts_count; i++) {
	rank = root_vertex(verts + i)->a % g->nprocs;
	sbuf[sdsps[rank]++] = verts[i];
    }

    /* restore sdsps */
    for (rank = 0; rank < g->nprocs; rank++)
	sdsps[rank] -= scnts[rank];

    /* send forward vertices (FIXME: only need to send (a,b) pairs) */
    MPI_Type_contiguous(sizeof(*sbuf), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(sbuf, scnts, sdsps, type, rbuf, rcnts, rdsps, type, g->comm);
    MPI_Type_free(&type);
    phgFree(sbuf);

    if (n > 0) {
	/* set up and sort ordering, set up block_start */
	verts_gathered = rbuf;	/* for comp_composite_vertex() */
	ordering = phgAlloc(2 * n * sizeof(*ordering));
	block_start = ordering + n;
	i = 0;
	for (rank = 0; rank < g->nprocs; rank++) {
	    assert(rdsps[rank] == i);
	    offset = i;
	    for (j = 0; j < rcnts[rank]; j++, i++) {
		ordering[i] = i;
		block_start[i] = offset;
	   }
	}
	qsort(ordering, n, sizeof(*ordering), comp_composite_vertex);

	/* assign indices */
	count = g->period->nvert_global;
	for (i = 0; i < n; ) {
	    for (j = i + 1; j < n; j++) {
		if (comp_composite_vertex(ordering + i, ordering + j) != 0)
		    break;
	    }
	    for (; i < j; i++)
		irbuf[ordering[i]] = count;
	    count++;
	}
	count -= g->period->nvert_global;

	phgFree(ordering);
	block_start = NULL;
	verts_gathered = NULL;
    }
    else {
	count = 0;
    }

    phgFree(rbuf);

    /* gather all counts */
    MPI_Allgather(&count, 1, PHG_MPI_INT, counts, 1, PHG_MPI_INT, g->comm);
    count = offset = 0;
    for (rank = 0; rank < g->nprocs; rank++) {
	if (rank == g->rank)
	    offset = count;
	count += counts[rank];
    }

    /* adjust indices */
    for (i = 0; i < n; i++)
	irbuf[i] += offset;

    /* update g->period->nvert_global */
    g->period->nvert_global += count;

    /* send back indices */
    MPI_Alltoallv(irbuf, rcnts, rdsps, PHG_MPI_INT,
		  isbuf, scnts, sdsps, PHG_MPI_INT, g->comm);

    /* update g->period->L2Gmap_vert */
    for (i = 0; i < verts_count; i++) {
	rank = root_vertex(verts + i)->a % g->nprocs;
	g->period->L2Gmap_vert[verts[i].index] = isbuf[sdsps[rank]++];
    }

    phgFree(isbuf);
    phgFree(scnts);
}
#else	/* 0 */
static void
period_update_shared(GRID *g, VERTEX_t *verts_shared, INT verts_count)
/* updates g->period->L2Gmap_vert and g->period->nvert_global for shared
 * vertices. The list of new vertices are in the buffer 'verts_shared' */
{
    INT i, j, *ordering;

    if (g->nprocs == 1) {
	verts_gathered = verts_shared;	/* for comp_composite_vertex() */
	ordering = phgAlloc(verts_count * sizeof(*ordering));
	for (i = 0; i < verts_count; i++)
	    ordering[i] = i;
	qsort(ordering, verts_count, sizeof(*ordering), comp_composite_vertex);

	/* assign indices */
	for (i = 0; i < verts_count; ) {
	    for (j = i + 1; j < verts_count; j++) {
		if (comp_composite_vertex(ordering + i, ordering + j) != 0)
		    break;
	    }
	    for (; i < j; i++)
		g->period->L2Gmap_vert[verts_shared[ordering[i]].index] =
							g->period->nvert_global;
	    g->period->nvert_global++;
	}
	phgFree(ordering);
	verts_gathered = NULL;

	return;
    }

}
#endif	/* 0 */


void
phgMapP2C(const ELEMENT *e, int Vmap[NVert], int Emap[NEdge], int child)
/* Vmap/Emap maps the vertices of e to those of its children, i.e.:
	e->verts[i] <=> e->children[child]->verts[Vmap[i]], i = 0,1,2,3
	e->egdes[i] <=> e->children[child]->egdes[Emap[i]], i = 0,1,2,3,4,5
   vertices/edges which do not belong to the child are mapped to -1 */
{
    switch (e->type) {
	case DIAGONAL:
	    if (child == 0) {
		if (Vmap != NULL) {
		    Vmap[0] = 0; Vmap[1] = -1; Vmap[2] = 1; Vmap[3] = 2;
		}
		if (Emap != NULL) {
		    Emap[0] = 2; Emap[1] = 0; Emap[2] = 1;
		    Emap[3] = -1; Emap[4] = -1; Emap[5] = 3;
		}
	    }
	    else {
		if (Vmap != NULL) {
		    Vmap[0] = -1; Vmap[1] = 0; Vmap[2] = 2; Vmap[3] = 1;
		}
		if (Emap != NULL) {
		    Emap[0] = 2; Emap[1] = -1; Emap[2] = -1;
		    Emap[3] = 1; Emap[4] = 0; Emap[5] = 3;
		}
	    }
	    break;
	case FACE:
	    if (child == 0) {
		if (Vmap != NULL) {
		    Vmap[0] = 0; Vmap[1] = -1; Vmap[2] = 1; Vmap[3] = 2;
		}
		if (Emap != NULL) {
		    Emap[0] = 2; Emap[1] = 0; Emap[2] = 1;
		    Emap[3] = -1; Emap[4] = -1; Emap[5] = 3;
		}
	    }
	    else {
		if (Vmap != NULL) {
		    Vmap[0] = -1; Vmap[1] = 0; Vmap[2] = 1; Vmap[3] = 2;
		}
		if (Emap != NULL) {
		    Emap[0] = 2; Emap[1] = -1; Emap[2] = -1;
		    Emap[3] = 0; Emap[4] = 1; Emap[5] = 3;
		}
	    }
	    break;
	case EDGE:
	    if (child == 0) {
		if (Vmap != NULL) {
		    Vmap[0] = 0; Vmap[1] = -1; Vmap[2] = 1; Vmap[3] = 2;
		}
		if (Emap != NULL) {
		    Emap[0] = 2; Emap[1] = 0; Emap[2] = 1;
		    Emap[3] = -1; Emap[4] = -1; Emap[5] = 3;
		}
	    }
	    else {
		if (Vmap != NULL) {
		    Vmap[0] = -1; Vmap[1] = 0; Vmap[2] = 1; Vmap[3] = 2;
		}
		if (Emap != NULL) {
		    Emap[0] = 2; Emap[1] = -1; Emap[2] = -1;
		    Emap[3] = 0; Emap[4] = 1; Emap[5] = 3;
		}
	    }
	    break;
	case OPPOSITE:
	    if (child == 0) {
		if (Vmap != NULL) {
		    Vmap[0] = 2; Vmap[1] = -1; Vmap[2] = 1; Vmap[3] = 0;
		}
		if (Emap != NULL) {
		    Emap[0] = 5; Emap[1] = 3; Emap[2] = 1;
		    Emap[3] = -1; Emap[4] = -1; Emap[5] = 0;
		}
	    }
	    else {
		if (Vmap != NULL) {
		    Vmap[0] = -1; Vmap[1] = 2; Vmap[2] = 1; Vmap[3] = 0;
		}
		if (Emap != NULL) {
		    Emap[0] = 5; Emap[1] = -1; Emap[2] = -1;
		    Emap[3] = 3; Emap[4] = 1; Emap[5] = 0;
		}
	    }
	    break;
	case MIXED:
	    if (child == 0) {
		if (Vmap != NULL) {
		    Vmap[0] = 0; Vmap[1] = -1; Vmap[2] = 1; Vmap[3] = 2;
		}
		if (Emap != NULL) {
		    Emap[0] = 2; Emap[1] = 0; Emap[2] = 1;
		    Emap[3] = -1; Emap[4] = -1; Emap[5] = 3;
		}
	    }
	    else {
		if (Vmap != NULL) {
		    Vmap[0] = -1; Vmap[1] = 2; Vmap[2] = 1; Vmap[3] = 0;
		}
		if (Emap != NULL) {
		    Emap[0] = 5; Emap[1] = -1; Emap[2] = -1;
		    Emap[3] = 3; Emap[4] = 1; Emap[5] = 0;
		}
	    }
	    break;
	default:
	    phgError(1, "%s:%d, invalid element type %d (element %d)\n",
			    __FILE__, __LINE__, e->type, e->index);
    }
}

void
phgMapC2P(const ELEMENT *e, int Vmap[NVert], int Emap[NEdge], int child)
/* computes Vmap[NVert]/Emap[NEdge] which maps the vertices/edges of
   children of e to the vertices/edges of e, i.e.:
	e->children[child]->verts[i] <=> e->verts[Vmap[i]], i = 0,1,2,3
	e->children[child]->edges[i] <=> e->edges[Emap[i]], i = 0,1,2,3,4,5
   Note: the last vertex (vertex 3) of the children is the new vertex) */
{
    if (Vmap != NULL)
	Vmap[3] = -1;
    switch (e->type) {
	case DIAGONAL:
	    if (child == 0) {
		if (Vmap != NULL) {
		    Vmap[0] = 0; Vmap[1] = 2; Vmap[2] = 3;
		}
		if (Emap != NULL) {
		    Emap[0] = 1; Emap[1] = 2; Emap[2] = 0;
		    Emap[3] = 5; Emap[4] = -1; Emap[5] = -1;
		}
	    }
	    else {
		if (Vmap != NULL) {
		    Vmap[0] = 1; Vmap[1] = 3; Vmap[2] = 2;
		}
		if (Emap != NULL) {
		    Emap[0] = 4; Emap[1] = 3; Emap[2] = 0;
		    Emap[3] = 5; Emap[4] = -1; Emap[5] = -1;
		}
	    }
	    break;
	case FACE:
	    if (child == 0) {
		if (Vmap != NULL) {
		    Vmap[0] = 0; Vmap[1] = 2; Vmap[2] = 3;
		}
		if (Emap != NULL) {
		    Emap[0] = 1; Emap[1] = 2; Emap[2] = 0;
		    Emap[3] = 5; Emap[4] = -1; Emap[5] = -1;
		}
	    }
	    else {
		if (Vmap != NULL) {
		    Vmap[0] = 1; Vmap[1] = 2; Vmap[2] = 3;
		}
		if (Emap != NULL) {
		    Emap[0] = 3; Emap[1] = 4; Emap[2] = 0;
		    Emap[3] = 5; Emap[4] = -1; Emap[5] = -1;
		}
	    }
	    break;
	case EDGE:
	    if (child == 0) {
		if (Vmap != NULL) {
		    Vmap[0] = 0; Vmap[1] = 2; Vmap[2] = 3;
		}
		if (Emap != NULL) {
		    Emap[0] = 1; Emap[1] = 2; Emap[2] = 0;
		    Emap[3] = 5; Emap[4] = -1; Emap[5] = -1;
		}
	    }
	    else {
		if (Vmap != NULL) {
		    Vmap[0] = 1; Vmap[1] = 2; Vmap[2] = 3;
		}
		if (Emap != NULL) {
		    Emap[0] = 3; Emap[1] = 4; Emap[2] = 0;
		    Emap[3] = 5; Emap[4] = -1; Emap[5] = -1;
		}
	    }
	    break;
	case OPPOSITE:
	    if (child == 0) {
		if (Vmap != NULL) {
		    Vmap[0] = 3; Vmap[1] = 2; Vmap[2] = 0;
		}
		if (Emap != NULL) {
		    Emap[0] = 5; Emap[1] = 2; Emap[2] = -1;
		    Emap[3] = 1; Emap[4] = -1; Emap[5] = 0;
		}
	    }
	    else {
		if (Vmap != NULL) {
		    Vmap[0] = 3; Vmap[1] = 2; Vmap[2] = 1;
		}
		if (Emap != NULL) {
		    Emap[0] = 5; Emap[1] = 4; Emap[2] = -1;
		    Emap[3] = 3; Emap[4] = -1; Emap[5] = 0;
		}
	    }
	    break;
	case MIXED:
	    if (child == 0) {
		if (Vmap != NULL) {
		    Vmap[0] = 0; Vmap[1] = 2; Vmap[2] = 3;
		}
		if (Emap != NULL) {
		    Emap[0] = 1; Emap[1] = 2; Emap[2] = 0;
		    Emap[3] = 5; Emap[4] = -1; Emap[5] = -1;
		}
	    }
	    else {
		if (Vmap != NULL) {
		    Vmap[0] = 3; Vmap[1] = 2; Vmap[2] = 1;
		}
		if (Emap != NULL) {
		    Emap[0] = 5; Emap[1] = 4; Emap[2] = -1;
		    Emap[3] = 3; Emap[4] = -1; Emap[5] = 0;
		}
	    }
	    break;
	default:
	    phgError(1, "%s:%d, invalid element type %d (element %d)\n",
			    __FILE__, __LINE__, e->type, e->index);
    }
}

static void
divide_simplex(GRID *g, ELEMENT *e)
/* this function actually divides (bisects) a simplex. */
{
    ELEMENT *e0, *e1;

    assert(Dim == 3);

    e0 = phgNewElements(1);
    e1 = phgNewElements(1);
    e0->region_mark = e1->region_mark = e->region_mark;
    e0->index = e1->index = g->nelem;	/* updated later by update_indices() */
    e->children[0] = e0;
    e->children[1] = e1;

    g->nleaf++;
    g->nelem++;
    g->ntree += 2;

    phgDebug((4, "divide element %p, children=%p, %p\n", e, e0, e1));

    /* new elements have their 'parent' member point to the parent */
    e0->parent = e1->parent = e;


    e0->generation = e1->generation = e->generation + 1;
    e0->mark  = e1->mark  = e->mark - 1;

    /* note: the child of e containing the vertex of the refinement edge
       with smaller global index inheritates edge no of the refinement edge,
       the other child gets a new number for its cut-edge. */

    switch (e->type) {
	case DIAGONAL:	/* diagonal of the reference parallelepiped */
	    e0->type = FACE;
	    e0->verts[0] = e->verts[0];	e0->bound_type[0] = INTERIOR;
	    e0->verts[1] = e->verts[2]; e0->bound_type[1] = e->bound_type[2];
	    e0->verts[2] = e->verts[3]; e0->bound_type[2] = e->bound_type[3];
	    e0->verts[3] = -1;		e0->bound_type[3] = e->bound_type[1];

	    e1->type = FACE;
	    e1->verts[0] = e->verts[1];	e1->bound_type[0] = INTERIOR;
	    e1->verts[1] = e->verts[3]; e1->bound_type[1] = e->bound_type[3];
	    e1->verts[2] = e->verts[2]; e1->bound_type[2] = e->bound_type[2];
	    e1->verts[3] = -1;		e1->bound_type[3] = e->bound_type[0];
	    break;
	case FACE:	/* diagonal of face of the reference parallelepiped */
	    e0->type = EDGE;
	    e0->verts[0] = e->verts[0];	e0->bound_type[0] = INTERIOR;
	    e0->verts[1] = e->verts[2]; e0->bound_type[1] = e->bound_type[2];
	    e0->verts[2] = e->verts[3]; e0->bound_type[2] = e->bound_type[3];
	    e0->verts[3] = -1;		e0->bound_type[3] = e->bound_type[1];

	    e1->type = EDGE;
	    e1->verts[0] = e->verts[1];	e1->bound_type[0] = INTERIOR;
	    e1->verts[1] = e->verts[2]; e1->bound_type[1] = e->bound_type[2];
	    e1->verts[2] = e->verts[3]; e1->bound_type[2] = e->bound_type[3];
	    e1->verts[3] = -1;		e1->bound_type[3] = e->bound_type[0];
	    break;
	case EDGE:	/* edge of the reference parallelepiped */
	    e0->type = DIAGONAL;
	    e0->verts[0] = e->verts[0];	e0->bound_type[0] = INTERIOR;
	    e0->verts[1] = e->verts[2]; e0->bound_type[1] = e->bound_type[2];
	    e0->verts[2] = e->verts[3]; e0->bound_type[2] = e->bound_type[3];
	    e0->verts[3] = -1;		e0->bound_type[3] = e->bound_type[1];

	    e1->type = DIAGONAL;
	    e1->verts[0] = e->verts[1];	e1->bound_type[0] = INTERIOR;
	    e1->verts[1] = e->verts[2]; e1->bound_type[1] = e->bound_type[2];
	    e1->verts[2] = e->verts[3]; e1->bound_type[2] = e->bound_type[3];
	    e1->verts[3] = -1;		e1->bound_type[3] = e->bound_type[0];
	    break;
	case OPPOSITE:	/* special case: type O (Liu/Arnold) */
	    e0->type = FACE;
	    e0->verts[0] = e->verts[3];	e0->bound_type[0] = e->bound_type[3];
	    e0->verts[1] = e->verts[2]; e0->bound_type[1] = e->bound_type[2];
	    e0->verts[2] = e->verts[0]; e0->bound_type[2] = INTERIOR;
	    e0->verts[3] = -1;		e0->bound_type[3] = e->bound_type[1];

	    e1->type = FACE;
	    e1->verts[0] = e->verts[3];	e1->bound_type[0] = e->bound_type[3];
	    e1->verts[1] = e->verts[2]; e1->bound_type[1] = e->bound_type[2];
	    e1->verts[2] = e->verts[1]; e1->bound_type[2] = INTERIOR;
	    e1->verts[3] = -1;		e1->bound_type[3] = e->bound_type[0];
	    break;
	case MIXED:	/* special case: type M (Liu/Arnold) */
	    e0->type = FACE;
	    e0->verts[0] = e->verts[0];	e0->bound_type[0] = INTERIOR;
	    e0->verts[1] = e->verts[2]; e0->bound_type[1] = e->bound_type[2];
	    e0->verts[2] = e->verts[3]; e0->bound_type[2] = e->bound_type[3];
	    e0->verts[3] = -1;		e0->bound_type[3] = e->bound_type[1];

	    e1->type = FACE;
	    e1->verts[0] = e->verts[3];	e1->bound_type[0] = e->bound_type[3];
	    e1->verts[1] = e->verts[2]; e1->bound_type[1] = e->bound_type[2];
	    e1->verts[2] = e->verts[1]; e1->bound_type[2] = INTERIOR;
	    e1->verts[3] = -1;		e1->bound_type[3] = e->bound_type[0];
	    break;
	default:
	    phgError(1, "%s:%d, invalid element type %d (element %d)\n",
			    __FILE__, __LINE__, e->type, e->index);
    }

#if ALLOW_CURVED_BOUNDARY
    /* copy boundary functions */
    switch (e->type) {
	case DIAGONAL:
	    e0->bound_func[2] = e1->bound_func[2] = e->bound_func[0];
	    e0->bound_func[0] = e->bound_func[1];
	    e0->bound_func[1] = e->bound_func[2];
	    e1->bound_func[1] = e->bound_func[3];
	    e1->bound_func[0] = e->bound_func[4];
	    e0->bound_func[3] = e1->bound_func[3] = e->bound_func[5];
	    if (e->bound_func[0] == e->bound_func[3] ||
		e->bound_func[0] == e->bound_func[1])
		e0->bound_func[4] = e1->bound_func[5] = e->bound_func[0];
	    if (e->bound_func[0] == e->bound_func[2] ||
		e->bound_func[0] == e->bound_func[4])
		e0->bound_func[5] = e1->bound_func[4] = e->bound_func[0];
	    break;
	case FACE:
	case EDGE:
	    e0->bound_func[2] = e1->bound_func[2] = e->bound_func[0];
	    e0->bound_func[0] = e->bound_func[1];
	    e0->bound_func[1] = e->bound_func[2];
	    e1->bound_func[0] = e->bound_func[3];
	    e1->bound_func[1] = e->bound_func[4];
	    e0->bound_func[3] = e1->bound_func[3] = e->bound_func[5];
	    if (e->bound_func[0] == e->bound_func[1] ||
		e->bound_func[0] == e->bound_func[3])
		e0->bound_func[4] = e1->bound_func[4] = e->bound_func[0];
	    if (e->bound_func[0] == e->bound_func[2] ||
		e->bound_func[0] == e->bound_func[4])
		e0->bound_func[5] = e1->bound_func[5] = e->bound_func[0];
	    break;
	case OPPOSITE:
	    e0->bound_func[5] = e1->bound_func[5] = e->bound_func[0];
	    e0->bound_func[3] = e->bound_func[1];
	    e0->bound_func[1] = e->bound_func[2];
	    e1->bound_func[3] = e->bound_func[3];
	    e1->bound_func[1] = e->bound_func[4];
	    e0->bound_func[0] = e1->bound_func[0] = e->bound_func[5];
	    if (e->bound_func[0] == e->bound_func[1] ||
		e->bound_func[0] == e->bound_func[3])
		e0->bound_func[4] = e1->bound_func[4] = e->bound_func[0];
	    if (e->bound_func[0] == e->bound_func[2] ||
		e->bound_func[0] == e->bound_func[4])
		e0->bound_func[2] = e1->bound_func[2] = e->bound_func[0];
	    break;
	case MIXED:
	    e0->bound_func[2] = e1->bound_func[5] = e->bound_func[0];
	    e0->bound_func[0] = e->bound_func[1];
	    e0->bound_func[1] = e->bound_func[2];
	    e1->bound_func[3] = e->bound_func[3];
	    e1->bound_func[1] = e->bound_func[4];
	    e0->bound_func[3] = e1->bound_func[0] = e->bound_func[5];
	    if (e->bound_func[0] == e->bound_func[1] ||
		e->bound_func[0] == e->bound_func[3])
		e0->bound_func[4] = e1->bound_func[4] = e->bound_func[0];
	    if (e->bound_func[0] == e->bound_func[2] ||
		e->bound_func[0] == e->bound_func[4])
		e0->bound_func[5] = e1->bound_func[2] = e->bound_func[0];
	    break;
    }
#if 0
    /* test */
    {
	EXPR **f;
	int i, j;
	COORD *c;
	ELEMENT *ee;
	FLOAT f0, f1;
	assert(FT_PHG == FT_DOUBLE);
	for (i = 0; i < 2; i++) {
	    ee = (i == 0) ? e0 : e1;
	    for (j = 0; j < NEdge; j++) {
		if (ee->bound_func[j] == -1)
		    continue;
		if (GetEdgeVertex(j, 1) == 3)
		    continue;	/* not yet available */
		f = g->bdry_funcs + ee->bound_func[j] * 4;
		c = g->verts + ee->verts[GetEdgeVertex(j, 0)];
		f0 = phgEvaluate3DFunction(*f, (*c)[0], (*c)[1], (*c)[2], 0.0);
		c = g->verts + ee->verts[GetEdgeVertex(j, 1)];
		f1 = phgEvaluate3DFunction(*f, (*c)[0], (*c)[1], (*c)[2], 0.0);
		if (Fabs(f0) >= 1e-10 || Fabs(f1) >= 1e-10) {
		    phgDumpElement(NULL, e);
		    phgInfo(-1, "bound_funcs:");
		    for (j = 0; j < NEdge; j++)
			printf(" %d", e->bound_func[j]);
		    printf("\n");
		    phgDumpElement(NULL, ee);
		    phgInfo(-1, "bound_funcs:");
		    for (j = 0; j < NEdge; j++)
			printf(" %d", ee->bound_func[j]);
		    printf("\n");
		    phgError(1, "edge %d, wrong function: f0=%lg, f1=%lg\n",
			j, (double)f0, (double)f1);
		}
	    }
	}
    }
#endif	/* 0 */
#endif	/* ALLOW_CURVED_BOUNDARY */
}

#if 0
/* the functions match_xxxx are mainly for testing */
static BOOLEAN
match_pair(ELEMENT *e0, ELEMENT *e1)
/* update links between a pair of neighbor elements (slow!) */
{
    int map[NFace] = {-1, -1, -1, -1};
    int i, j, count = 0;

    for (i = 0; i < NFace; i++)
	for (j = 0; j < NFace; j++)
	    if (CompVertex(e0->verts[i], e1->verts[j]) == 0) {
		map[i] = j;
		count++;
		break;
	    }
    if (count != 3)
	return FALSE;

    for (j = 0; j < NFace; j++) {
	if (map[j] == -1) {
	    i = j;
	    break;
	}
    }

    j = opposite_vertex(g_, e0, i, e1);
    e0->neighbours[i] = e1;
    e1->neighbours[j] = e0;

    return TRUE;
}

static void match_child(ELEMENT *e0, ELEMENT *e)
/* find (and match) the child of e0 which is a neighbor of e */
{
    if (match_pair(e0->children[0], e) || match_pair(e0->children[1], e))
	return;

    /* shouldn't occur */
    phgDumpElement(NULL, e0);
    phgDumpElement(NULL, e0->children[0]);
    phgDumpElement(NULL, e0->children[1]);
    phgDumpElement(NULL, e);
    phgError(1, "error matching neighbours, abort.\n");
}
#endif	/* 0 */

static int
update_neighbours(GRID *g, ELEMENT *e)
/* update neighbour information for a newly bisected simplex.
   returns !0 if the children have hanging nodes with the two neighbours
   not sharing the refinement edge (bit 0 for child 0, bit 1 for child 1) */
{
    int a0 = 0, a1 = 0; /* to make gcc happy */
    int b0, b1, c0, c1;
    ELEMENT *e0;
    int i, ret = 0;

    /* neighbour information between the two children */

    switch (e->type) {
	case OPPOSITE:
	    e->children[0]->neighbours[2] = e->children[1];
	    e->children[1]->neighbours[2] = e->children[0];
	    break;
	case MIXED:
	    e->children[0]->neighbours[0] = e->children[1];
	    e->children[1]->neighbours[2] = e->children[0];
	    break;
	default:
	    e->children[0]->neighbours[0] = e->children[1];
	    e->children[1]->neighbours[0] = e->children[0];
    }
    e->children[0]->neighbours[3] = e->neighbours[1];
    e->children[1]->neighbours[3] = e->neighbours[0];

    /* Case 1: the two neighbours not sharing the refinement edge */

    if ((e0 = GetNeighbour(e, 0)) != NULL) {
	if (GetNeighbour(e0, 0) == e)
	    i = 0;
	else if (GetNeighbour(e0, 1) == e)
	    i = 1;
	else if (GetNeighbour(e0, 2) == e)
	    i = 2;
	else if (GetNeighbour(e0, 3) == e)
	    i = 3;
	else
	    i = -1;
	assert(i != -1);
	e0->neighbours[i] = e->children[1];
	if (e0->children[0] != NULL && opposite_vertex(g, e, 0, e0) >= 2)
	    ret |= 2;
    }

    if ((e0 = GetNeighbour(e, 1)) != NULL) {
	if (GetNeighbour(e0, 0) == e)
	    i = 0;
	else if (GetNeighbour(e0, 1) == e)
	    i = 1;
	else if (GetNeighbour(e0, 2) == e)
	    i = 2;
	else if (GetNeighbour(e0, 3) == e)
	    i = 3;
	else
	    i = -1;
	assert(i != -1);
	e0->neighbours[i] = e->children[0];
	if (e0->children[0] != NULL && opposite_vertex(g, e, 1, e0) >= 2)
	    ret |= 1;
    }

    /* Case 2: the two neighbours sharing the refinement edge */

    c0 = 1;
    c1 = (e->type == DIAGONAL) ? 2 : 1;
    if ((e0 = GetNeighbour(e, 2)) != NULL) {
	if (e0->children[0] != NULL) {
	    assert((CompVertex(e0->verts[0], e->verts[0]) == 0 &&
		    CompVertex(e0->verts[1], e->verts[1]) == 0) ||
		   (CompVertex(e0->verts[0], e->verts[1]) == 0 &&
		    CompVertex(e0->verts[1], e->verts[0]) == 0));
	    a0 = (CompVertex(e0->verts[0], e->verts[0]) == 0 ? 0 : 1);
	    a1 = 1 - a0;
	    e->children[0]->neighbours[c0] = e0->children[a0];
	    e->children[1]->neighbours[c1] = e0->children[a1];
	    /* the lines below are only necessary for refine_path */
	    if (e0->neighbours[2] == e/*opposite_vertex(g, e, 2, e0) == 2*/) {
		b0 = 1;
		b1 = (e0->type == DIAGONAL) ? 2 : 1;
	    }
	    else {
		assert(e0->neighbours[3] == e);
		b0 = (e0->type == OPPOSITE) ? 0 : 2;
		b1 = (e0->type == DIAGONAL) ? 1 :
			((e0->type == OPPOSITE || e0->type == MIXED) ? 0 : 2);
	    }
	    e0->children[0]->neighbours[b0] = e->children[a0];
	    if (e0->children[0]->children[0] != NULL && b0 <= 1) {
#if 1
		e0->children[0]->children[1 - b0]->neighbours[3] =
			e->children[a0];
		e->children[a0]->neighbours[a0 == 0 ? c0 : c1] =
			e0->children[0]->children[1 - b0];
#else	/* 1 */
		match_child(e0->children[0], e->children[a0]);
#endif	/* 1 */
	    }
	    e0->children[1]->neighbours[b1] = e->children[a1];
	    if (e0->children[1]->children[0] != NULL && b1 <= 1) {
#if 1
		e0->children[1]->children[1 - b1]->neighbours[3] =
			e->children[a1];
		e->children[a1]->neighbours[a1 == 0 ? c0 : c1] =
			e0->children[1]->children[1 - b1];
#else	/* 1 */
		match_child(e0->children[1], e->children[a1]);
#endif	/* 1 */
	    }
	}
    }

    c0 = (e->type == OPPOSITE) ? 0 : 2;
    c1 = (e->type == DIAGONAL) ? 1 :
		((e->type == OPPOSITE || e->type == MIXED) ? 0 : 2);
    if ((e0 = GetNeighbour(e, 3)) != NULL) {
	if (e0->children[0] != NULL) {
	    assert((CompVertex(e0->verts[0], e->verts[0]) == 0 &&
		    CompVertex(e0->verts[1], e->verts[1]) == 0) ||
		   (CompVertex(e0->verts[0], e->verts[1]) == 0 &&
		    CompVertex(e0->verts[1], e->verts[0]) == 0));
	    a0 = (CompVertex(e0->verts[0], e->verts[0]) == 0 ? 0 : 1);
	    a1 = 1 - a0;
	    e->children[0]->neighbours[c0] = e0->children[a0];
	    e->children[1]->neighbours[c1] = e0->children[a1];
	    /* the lines below are only necessary for refine_path */
	    if (e0->neighbours[2] == e/*opposite_vertex(g, e, 3, e0) == 2*/) {
		b0 = 1;
		b1 = (e0->type == DIAGONAL) ? 2 : 1;
	    }
	    else {
		assert(e0->neighbours[3] == e);
		b0 = (e0->type == OPPOSITE) ? 0 : 2;
		b1 = (e0->type == DIAGONAL) ? 1 :
			((e0->type == OPPOSITE || e0->type == MIXED) ? 0 : 2);
	    }
	    e0->children[0]->neighbours[b0] = e->children[a0];
	    if (e0->children[0]->children[0] != NULL && b0 <= 1) {
#if 1
		e0->children[0]->children[1 - b0]->neighbours[3] =
			e->children[a0];
		e->children[a0]->neighbours[a0 == 0 ? c0 : c1] =
			e0->children[0]->children[1 - b0];
#else	/* 1 */
		match_child(e0->children[0], e->children[a0]);
#endif	/* 1 */
	    }
	    e0->children[1]->neighbours[b1] = e->children[a1];
	    if (e0->children[1]->children[0] != NULL && b1 <= 1) {
#if 1
		e0->children[1]->children[1 - b1]->neighbours[3] =
			e->children[a1];
		e->children[a1]->neighbours[a1 == 0 ? c0 : c1] =
			e0->children[1]->children[1 - b1];
#else	/* 1 */
		match_child(e0->children[1], e->children[a1]);
#endif	/* 1 */
	    }
	}
    }

    return ret;
}

static size_t g_verts_allocated;
static INT too_close;

static void
add_vertex(GRID *g, ELEMENT *e)
/* insert the coordinates of a new vertex into g->verts[] */
{
    COORD *v, *v0, *v1;
    INT a = e->verts[0], b = e->verts[1];

    if (g->nvert >= g_verts_allocated) {
	g->verts = phgRealloc_(g->verts,
		sizeof(*(g->verts)) * (g_verts_allocated + 32768),
		sizeof(*(g->verts)) * g_verts_allocated);
	g_verts_allocated += 32768;
    }
    v = g->verts + (g->nvert++);
    v0 = g->verts + a;
    v1 = g->verts + b;
   (*v)[0] = ((*v0)[0] + (*v1)[0]) * .5;
   (*v)[1] = ((*v0)[1] + (*v1)[1]) * .5;
   (*v)[2] = ((*v0)[2] + (*v1)[2]) * .5;
#if ALLOW_CURVED_BOUNDARY
#if 0	/* testing code */
    {
	EXPR **f = g->bdry_funcs;
	FLOAT *x0 = (*v0) + 0, *y0 = (*v0) + 1, *z0 = (*v0) + 2;
	FLOAT *x1 = (*v1) + 0, *y1 = (*v1) + 1, *z1 = (*v1) + 2;
	e->bound_func[0] = -1;
	while (*f != NULL) {
	    if (fabs(phgEvaluate3DFunction(*f, *x0, *y0, *z0, 0.0)) < 1e-10 &&
		fabs(phgEvaluate3DFunction(*f, *x1, *y1, *z1, 0.0)) < 1e-10) {
		e->bound_func[0] = (f - g->bdry_funcs) / 4;
		break;
	    }
	    f += 4;
	}
    }
#endif	/* 0 */
    if (e->bound_func[0] != -1) {
	EXPR **f = g->bdry_funcs + e->bound_func[0] * 4;
	double x = (*v)[0], y = (*v)[1], z = (*v)[2];
	(*v)[0] = phgEvaluate3DFunction(f[1], x, y, z, 0.0);
	(*v)[1] = phgEvaluate3DFunction(f[2], x, y, z, 0.0);
	(*v)[2] = phgEvaluate3DFunction(f[3], x, y, z, 0.0);
	phgDebug((3, "curved bdry: (%lg, %lg, %lg) => (%lg, %lg, %lg)\n",
		(double)x, (double)y, (double)z,
		(double)(*v)[0], (double)(*v)[1], (double)(*v)[2]));
    }
#endif	/* ALLOW_CURVED_BOUNDARY */
    if (((*v)[0] == (*v0)[0] || (*v)[0] == (*v1)[0]) &&
	((*v)[1] == (*v0)[1] || (*v)[1] == (*v1)[1]) &&
	((*v)[2] == (*v0)[2] || (*v)[2] == (*v1)[2]))
	too_close++;
}

static int max_recur = 0;

#define CompTuple(a0, b0, a1, b1) ((a0 == a1) ? b0 - b1 : a0 - a1)

static void
new_vertex(GRID *g, ELEMENT *e)
/* If e != NULL saves the new vertex created by bisecting 'e' to a table.
 * If e == NULL updates g->L2Gmap_vert / g->period and resets the table.
 *
 * Note: this function makes use of a hash table for vertex searching */
{
    typedef struct {
	/* a < b < m, a and b the old vertices, m the new vertex */
	INT a, b, m;
	/* whether the vertex has been added to the list of shared vertices */
	BOOLEAN shared;
    } TABLE;

    static struct {
	TABLE	*table;
	short	count;
	short	allocated;
    } hash[1499], *h;
    static size_t hash_size = sizeof(hash) / sizeof(hash[0]);
    static BOOLEAN initialized = FALSE;

    static BTYPE *types = NULL;
    static size_t types_allocated = 0;

    TABLE *table, *first, *last, *mid;
    size_t allocated, count;
    INT a, b, m;
    INT i;

    if (!initialized) {
	initialized = TRUE;
	memset(hash, 0, sizeof(hash));
    }

    if (e == NULL) {
	VERTEX_t *verts_shared = NULL;	/* shared vertices */

	/* assign global indices for new vertices and flush the hash tables. */

	count = g->nvert - nvert_save;

	if (g->period != NULL) {
	    /* update g->period->L2Gmap_vert */
	    for (i = 0; i < count; i++)
		types[i] = 1;		/* mark shared vertices */
	    verts_count = 0;		/* # shared vertices */
	    for (h = hash; h < hash + hash_size; h++) {
		if ((table = h->table) == NULL)
		    continue;
		last = table + h->count;
		for (mid = table; mid < last; mid++) {
		    if (mid->shared) {
			verts_count++;
			continue;
		    }
		    i = mid->m - nvert_save;
		    assert(i >= 0 && i < count);
		    types[i] = 0;
		    verts[i].index = -1;	/* for debugging */
		}
	    }

	    /* collect vertices */
	    verts_shared = phgAlloc(verts_count * sizeof(*verts_shared));
	    verts_count = 0;
	    for (i = 0; i < count; i++) {
		if (types[i] == 0)
		    continue;
		assert(verts[i].index == i + nvert_save);
		verts_shared[verts_count++] = verts[i];
	    }

	    phgFree(verts);
	    verts = NULL;

	    g->period->L2Gmap_vert = phgRealloc_(g->period->L2Gmap_vert,
			g->nvert * sizeof(*g->period->L2Gmap_vert),
			nvert_save * sizeof(*g->period->L2Gmap_vert));
	    phgFree(g->period->ordering);
	    g->period->ordering = NULL;
	    for (i = nvert_save; i < g->nvert; i++)
		g->period->L2Gmap_vert[i] = -1;

	    /* update period->L2Gmap_vert for shared vertices */
	    period_update_shared(g, verts_shared, verts_count);
	    phgFree(verts_shared);
	    verts_shared = NULL;

	    /* update period->L2Gmap_vert for private vertices */
	    period_update_private(g);

	    /* Note: verts_count and first is needed by next code block */
	}

	{
	    g->nvert_global = g->nvert;
	    g->nelem_global = g->nelem;
	}

	phgFree(verts_shared);
	verts_shared = NULL;
	phgFree(verts);
	verts = NULL;
	verts_count = verts_allocated = 0;

	phgFree(types);
	types = NULL;
	types_allocated = 0;

	/* reset hash tables */
	for (h = hash; h < hash + hash_size; h++)
	    phgFree(h->table);
	memset(hash, 0, sizeof(hash));

	return;
    }

    if ((a = e->verts[0]) > (b = e->verts[1])) {
	m = a;
	a = b;
	b = m;
    }

    h = hash + (size_t)((b + (a % hash_size) * 1493) % hash_size);
    table = h->table;
    count = h->count;
    allocated = h->allocated;

#if 0
    /* linear search */
    for (mid = table; mid < table + count; mid++)
	if (mid->m == m)
	    break;
    if (mid - table == count) {
	if (count >= allocated) {
	    table = phgRealloc_(table, (allocated + 16) * sizeof(*table),
					allocated * sizeof(*table));
	    allocated += 16;
	}
    }
    mid = table + count;
    mid->a = a;
    mid->b = b;
    mid->m = m = g->nvert;
    assert(g->period != NULL);	/* unimplemented */
    mid->shared = FALSE;
#else	/* 0 */
    /* binary search */
    i = 0;
    if (count > 0) {
	first = table;
	last = table + count - 1;
	if ((i = CompTuple(a, b, first->a, first->b)) == 0) {
	    m = first->m;
	    mid = first;
	    goto found;
	}
	else if (i < 0) {
	    i = 0;
	    goto new;
	}
	if (last > first) {
	    if ((i = CompTuple(a, b, last->a, last->b)) == 0) {
		m = last->m;
		mid = last;
		goto found;
	    }
	    else if (i > 0) {
		i = count;
		goto new;
	    }
	}
	while (last - first > 1) {
	    mid = first + (last - first + 1) / 2;
	    if ((i = CompTuple(a, b, mid->a, mid->b)) == 0) {
		m = mid->m;
		goto found;
	    }
	    if (i < 0)
		last = mid;
	    else
		first = mid;
	}
	i = first - table + 1;
    }

new:
    /* create a new entry and insert it before the i-th entry */
    phgDebug((4, "Binary search with (%d, %d) located entry %d\n", a, b, i));
    if (count >= allocated) {
	table = phgRealloc_(table, (allocated + 16) * sizeof(*table),
				   allocated * sizeof(*table));
	allocated += 16;
    }
    if (count - i > 0) {
#if 0
	memmove(table + i + 1, table + i, (count - i) * sizeof(*table));
#else	/* 0 */
	int j;
	for (j = count; j > i; j--)
	    table[j] = table[j - 1];
#endif	/* 0 */
    }
    mid = table + i;
    mid->a = a;
    mid->b = b;
    mid->m = m = g->nvert;
    add_vertex(g, e);
    if (g->period != NULL) {
	/* The periodic vertices have the REMOTE flag set, we
	 * save the REMOTE flag in types[] array (the REMOTE flag is set
	 * iff both parents have the REMOTE flag set) */
	i = m - nvert_save;
	if (i >= verts_allocated) {
	    assert(i == verts_allocated);
	    types = phgRealloc_(types,
			(types_allocated + 32768) * sizeof(*types),
			types_allocated * sizeof(*types));
	    types_allocated += 32768;
	    verts = phgRealloc_(verts,
			(verts_allocated + 32768) * sizeof(*verts),
			verts_allocated * sizeof(*verts));
	    verts_allocated += 32768;
	    for (; i < verts_allocated; i++) {
		types[i] = 0;
		verts[i].index = -1;
	    }
	    i = m - nvert_save;
	}
	types[i] = (a < nvert_save ? g->types_vert[a] : types[a - nvert_save]) &
		   (b < nvert_save ? g->types_vert[b] : types[b - nvert_save]);
	verts[i].index = m;
	a = (a < nvert_save ? GlobalVertexP(g, a) : a + g->nvert_global);
	b = (b < nvert_save ? GlobalVertexP(g, b) : b + g->nvert_global);
	if (comp_vertex(a, b) < 0) {
	    verts[i].a = a;
	    verts[i].b = b;
	}
	else {
	    verts[i].a = b;
	    verts[i].b = a;
	}
    }
    mid->shared = FALSE;
    count++;

found:
#endif	/* 0 */
    phgDebug((4, "index of the new vertex=%d\n", m));
    e->children[0]->verts[3] = e->children[1]->verts[3] = m;

    /* register shared vertex */
    if (!mid->shared && ((e->bound_type[2] | e->bound_type[3] |
		(types == NULL ? 0 : types[m - nvert_save])) & REMOTE))
	mid->shared = TRUE;

    h->table = table;
    h->count = count;
    h->allocated = allocated;
}

#if 0

/* recursive version of refine_path (does not work with cyclic paths) */

static int recur = 0;

static void
refine_path(GRID *g, ELEMENT *e)
/* refines the element e as well as any other elements to satisfy conformity.
   This function implements the hanging face approach, i.e., it follows
   paths of incompatible faces created by refinement of e. */
{
    int c0, c1;
    ELEMENT *e0;

    phgDebug((4, "%d refining element %p, type=%c, vertices=%d %d %d %d.\n",
		recur, e, GTypeName(e->type)[0],
		e->verts[0], e->verts[1], e->verts[2], e->verts[3]));

    if (recur >= REFINE_RECURSION_LIMIT)
	phgError(1, "recursion limit (%d) exceeded.\n", REFINE_RECURSION_LIMIT);
    if (max_recur < recur)
	max_recur = recur;

    /* divide the element */
    divide_simplex(g, e);
    new_vertex(g, e);

    /* update neighbour information */
    update_neighbours(g, e);

    /* the two neighbours sharing the refinement edge */

    c0 = 1;
    c1 = (e->type == DIAGONAL) ? 2 : 1;
    if ((e0 = GetNeighbour(e, 2)) != NULL) {
	/* recursively refine the 1st neighbour sharing the refinement edge */
	static int a0, a1;
    	phgDebug((4, "%d neighbour 2=%p, type=%c, vertices=%d %d %d %d.\n",
		recur, e0, GTypeName(e0->type)[0],
		e0->verts[0], e0->verts[1], e0->verts[2], e0->verts[3]));
	if (GetNeighbour(e0, 0) == e || GetNeighbour(e0, 1) == e) {
	    /* the neighbour does not share the refinement edge of e */
	    phgDebug((4, "%d neighbour does not share refinement edge.\n",
			recur));
	    if (e0->children[0] == NULL) {
		recur++;
		refine_path(g, e0);
		recur--;
		assert(e->children[0]->children[0] == NULL);
		assert(e->children[1]->children[0] == NULL);
	    }
	    e0 = GetNeighbour(e, 2);
	}
	if (e0->children[0] == NULL) {
	    recur++;
	    refine_path(g, e0);
	    recur--;
	    assert(e->children[0]->children[0] == NULL);
	    assert(e->children[1]->children[0] == NULL);
	}
	assert(e0->children[0]->children[0] == NULL);
	assert(e0->children[1]->children[0] == NULL);
	a0 = (CompVertex(e0->verts[0], e->verts[0]) == 0 ? 0 : 1);
	a1 = 1 - a0;
	e->children[0]->neighbours[c0] = e0->children[a0];
	e->children[1]->neighbours[c1] = e0->children[a1];
    }

    c0 = (e->type == OPPOSITE) ? 0 : 2;
    c1 = (e->type == DIAGONAL) ? 1 :
		((e->type == OPPOSITE || e->type == MIXED) ? 0 : 2);
    if ((e0 = GetNeighbour(e, 3)) != NULL) {
	/* recursively refine the 2nd neighbour sharing the refinement edge */
	static int b0, b1;
    	phgDebug((4, "%d neighbour 3=%p, type=%c, vertices=%d %d %d %d.\n",
		recur, e0, GTypeName(e0->type)[0],
		e0->verts[0], e0->verts[1], e0->verts[2], e0->verts[3]));
	if (GetNeighbour(e0, 0) == e || GetNeighbour(e0, 1) == e) {
	    /* the neighbour does not share the refinement edge of e */
	    phgDebug((4, "%d neighbour does not share refinement edge.\n",
			recur));
	    if (e0->children[0] == NULL) {
		recur++;
		refine_path(g, e0);
		recur--;
		assert(e->children[0]->children[0] == NULL);
		assert(e->children[1]->children[0] == NULL);
	    }
	    e0 = GetNeighbour(e, 3);
	}
	if (e0->children[0] == NULL) {
	    recur++;
	    refine_path(g, e0);
	    recur--;
	    assert(e->children[0]->children[0] == NULL);
	    assert(e->children[1]->children[0] == NULL);
	}
	assert(e0->children[0]->children[0] == NULL);
	assert(e0->children[1]->children[0] == NULL);
	b0 = (CompVertex(e0->verts[0], e->verts[0]) == 0 ? 0 : 1);
	b1 = 1 - b0;
	e->children[0]->neighbours[c0] = e0->children[b0];
	e->children[1]->neighbours[c1] = e0->children[b1];
    }
}

#else	/* 0 */

/* non recursive version of refine_path, more robust but slower(?) */

static BOOLEAN
share_refinement_edge(const ELEMENT *e1, const ELEMENT *e2)
/* returns TRUE if e2 contains refinement edge of e1, FALSE otherwise */
{
    int i;
    INT u;

    u = e1->verts[0];
    for (i = 0; i < NVert; i++)
	if (!CompVertex(u, e2->verts[i]))
	    break;
    if (i >= NVert)
	return FALSE;

    u = e1->verts[1];
    for (i = 0; i < NVert; i++)
	if (!CompVertex(u, e2->verts[i]))
	    break;
    return (i < NVert);
}

static BOOLEAN
match_neighbour(ELEMENT *e0, ELEMENT *e1)
/* finds whether e0 and e1 share a common face, if they do then returns TRUE
 * and updates the neighbours[] entries accordingly */
{
    int i, j, k, m0[] = {-1,-1,-1,-1}, m1[] = {-1,-1,-1,-1};

    /* Here we do by brute force, may be more efficient
     * to track refinement path of both elements */
    for (k = i = 0; i < NVert; i++) {
	for (j = 0; j < NVert; j++) {
	    if (m1[j] >= 0)
		continue;
	    if (!CompVertex(e0->verts[i], e1->verts[j])) {
		k++;
		m0[i] = j;
		m1[j] = i;
		break;
	    }
	}
    }

    if (k == 3) {
	for (k = 0; k < NVert; k++) {
	    if (m0[k] == -1)
		i = k;
	    if (m1[k] == -1)
		j = k;
	}
	e0->neighbours[i] = e1;
	e1->neighbours[j] = e0;
	return TRUE;
    }

    return FALSE;
}

static void
refine_path(GRID *g, ELEMENT *e)
/* refines the element e as well as any other elements to satisfy conformity.
   This function implements the hanging face approach, i.e., it follows
   paths of incompatible faces created by refinement of e. */
{
    int i;
    ELEMENT *e0;
#if 1
    /* use static variables for the stack, re-entry not allowed */
    static ELEMENT **stack = NULL;	/* pointers to the hanging nodes */
    static size_t stack_ptr = 0, stack_max = 0, stack_allocated = 0;
    FreeAtExit(stack);
#else
    /* use dynamic variables for the stack (TO CHECK: allowing re-entry?) */
    ELEMENT **stack = NULL;	/* pointers to the hanging nodes */
    size_t stack_ptr = 0, stack_max = 0, stack_allocated = 0;
#endif

#define ADD_TO_LIST(e) {						\
    if (stack_ptr >= stack_allocated) {					\
	stack = phgRealloc_(stack,					\
			(stack_allocated + 16) * sizeof(*stack),	\
			stack_allocated * sizeof(*stack));		\
	stack_allocated += 16;						\
    }									\
    stack[stack_ptr++] = e;						\
    if (stack_ptr > stack_max)						\
	stack_max = stack_ptr;						\
    phgDebug((4, "\tAdding %p to list, stack=%d\n", e, stack_ptr));	\
}

    ADD_TO_LIST(e)

    while (stack_ptr > 0) {
	do {
#if 1
	    /* first in last out */
	    e = stack[--stack_ptr];
#else	/* 1 */
	    /* first in first out (does not work! why?) */
	    e = *stack;
	    if (--stack_ptr > 0)
		memmove(stack, stack + 1, sizeof(*stack) * stack_ptr);
#endif	/* 1 */
	    if (e->children[0] == NULL)
		break;
	    if (stack_ptr <= 0) {
		e = NULL;
		break;
	    }
	} while (TRUE);

	if (e == NULL)
	    break;

	/* divide the element */
	divide_simplex(g, e);
	new_vertex(g, e);
	phgDebug((4, "Divide: (%c) %p -> (%c) %p %p\n", GTypeName(e->type)[0],
			e, GTypeName(e->children[0]->type)[0], e->children[0],
			e->children[1]));

#if 0
#warning remove me!
for (itmp = 0; g->rank == 1582 && itmp < ntmp_list; itmp++)
    if (!memcmp(e->verts, tmp_list[itmp].verts, sizeof(tmp_list[itmp].verts)))
	break;
if (g->rank==1582 && itmp < ntmp_list) {
  phgInfo(-1, "------------------------------- %s:\n", tmp_list[itmp].desc);
  phgDumpElement(g, e);
  if (itmp == 0) {
    for (e0 = e; e0->generation > 0; e0 = e0->parent);
    phgInfo(-1, "------------------------------- The tree containing E1\n");
    phgDumpBranch(g, e0);
    phgInfo(-1, "----------------------------------------------------\n"); 
  }
}
#endif

	/* update neighbour information */
	switch (update_neighbours(g, e)) {
	    case 1:
		ADD_TO_LIST(e->children[0])
		break;
	    case 2:
		ADD_TO_LIST(e->children[1])
		break;
	    case 3:
		ADD_TO_LIST(e->children[0])
		ADD_TO_LIST(e->children[1])
		break;
	}

	/* check the two neighbours sharing the refinement edge */
	for (i = 2; i < 4; i++) {
	    if ((e0 = GetNeighbour(e, i)) == NULL)
		continue;
	    if (e0->children[0] == NULL) {
		ADD_TO_LIST(e0)
		continue;
	    }
	    /* check whether children of e0 have been refined */
	    int j, k, l;
	    if ((j = opposite_vertex(g, e, i, e0)) < 2) {
		e->neighbours[i] = (e0 = e0->children[1 - j]);
		if (e0->children[0] == NULL) {
		    ADD_TO_LIST(e0)
		    continue;
		}
	    }
	    for (k = 0; k < 2; k++) {	/* loop on e0's children */
		ELEMENT *e1, *ee;
		if ((e1 = e0->children[k])->children[0] == NULL)
		    continue;	/* e1 has no children */
		/* Pairs of neighbouring elements:
		 * 			 e <-> e0
		 *	ee:=e->children[j] <-> e1:=e0->children[k] */
		j = (!CompVertex(e0->verts[0], e->verts[0])) ? k : 1 - k;
		ee = e->children[j];
		if (share_refinement_edge(e1, ee)) {
		    ADD_TO_LIST(ee)
		    continue;
		}
		/* Added on 2020.05.16:
		 * ee does not contain refinement edge of e1, then one of e1's
		 * children must be ee's neighbour, we need to check whether
		 * this neighbour has been refined and ee contains its
		 * refinement edge, and mark ee to be refined if yes.
		 * (see illustration backup/refine-debug-202005.dia). */
		for (l = 0; l < 2; l++) {
		    if ((e1->children[l])->children[0] != NULL &&
			match_neighbour(ee, e1->children[l])) {
			if (share_refinement_edge(e1->children[l], ee))
			    ADD_TO_LIST(ee)
			break;
		    }
		}
	    }
	}
    }

    phgInfo(4, "refine_path: stack_max = %d\n", stack_max);
    if (max_recur < stack_max)
	max_recur = stack_max;

    phgFree(stack);
    stack = NULL;
    stack_ptr = stack_max = stack_allocated = 0;
}

#endif	/* 0 */

static int
comp_edge(const double *length, ELEMENT *e, int i0, int i1)
/* Compare two edges (global indices) by first comparing their length,
   then their indices. */
{
    INT i;
    INT e0 = e->edges[i0], e1 = e->edges[i1];
    double d = length[e0] - length[e1];

    if (d == 0)
	return (i = e0 - e1) > 0 ? 1 : (i < 0 ? -1 : 0);

    return d < 0 ? -1 : 1;
}

static void
permute_vertices(ELEMENT *e, int v0, int v1)
{
    int i, a0 = 0, b0 = 0, a1 = 0, b1 = 0;
    INT v;
    void *p;

    if (v0 == v1)
	return;

    if (v0 > v1) {
	i = v0;
	v0 = v1;
	v1 = i;
    }

    if (g_->dof != NULL)
	phgWarning("TODO: permute DOFs!\n");

    /* verts[] */
    v = e->verts[v0];
    e->verts[v0] = e->verts[v1];
    e->verts[v1] = v;

    /* neighbours[] and bound_type[] */
    p = e->neighbours[v0];
    e->neighbours[v0] = e->neighbours[v1];
    e->neighbours[v1] = p;

    i = e->bound_type[v0];
    e->bound_type[v0] = e->bound_type[v1];
    e->bound_type[v1] = i;

    /* edges[], permute two pairs of edges and boundary function indices */
    switch (10 * v0 + v1) {
	case 01:
	    a0 = 1; b0 = 3;	/* 0-2 and 1-2 */
	    a1 = 2; b1 = 4;	/* 0-3 and 1-3 */
	    break;
	case 02:
	    a0 = 0; b0 = 3;	/* 0-1 and 1-2 */
	    a1 = 2; b1 = 5;	/* 0-3 and 2-3 */
	    break;
	case 03:
	    a0 = 0; b0 = 4;	/* 0-1 and 1-3 */
	    a1 = 1; b1 = 5;	/* 0-2 and 2-3 */
	    break;
	case 12:
	    a0 = 0; b0 = 1;	/* 0-1 and 0-2 */
	    a1 = 4; b1 = 5;	/* 1-3 and 2-3 */
	    break;
	case 13:
	    a0 = 0; b0 = 2;	/* 0-1 and 0-3 */
	    a1 = 3; b1 = 5;	/* 1-2 and 2-3 */
	    break;
	case 23:
	    a0 = 1; b0 = 2;	/* 0-2 and 0-3 */
	    a1 = 3; b1 = 4;	/* 1-2 and 1-3 */
	    break;
    }

    v = e->edges[a0];
    e->edges[a0] = e->edges[b0];
    e->edges[b0] = v;

    v = e->edges[a1];
    e->edges[a1] = e->edges[b1];
    e->edges[b1] = v;

#if ALLOW_CURVED_BOUNDARY
    v = e->bound_func[a0];
    e->bound_func[a0] = e->bound_func[b0];
    e->bound_func[b0] = v;

    v = e->bound_func[a1];
    e->bound_func[a1] = e->bound_func[b1];
    e->bound_func[b1] = v;
#endif	/* ALLOW_CURVED_BOUNDARY */

#if Dim == 3
    /* permute face indices */
    v = e->faces[v0];
    e->faces[v0] = e->faces[v1];
    e->faces[v1] = v;
#endif	/* Dim == 3 */
}

void
phgRefineInit(GRID *g, BOOLEAN do_marking)
/* Compute initial marking (embedding) using Liu/Arnold's algorithm.
   The function only works with flat (unrefined), undistributed meshes */
{
    INT i, bad_edges = 0;
    ELEMENT *e;
    int j, k;
    double *length;	/* (square of) length of all edges */

    FunctionEntry;

    if (g == NULL)
	return;

    g_ = g;

    if (!do_marking)
	goto pre_refine;

    phgInfo(2, "marking edges using Liu/Arnold's algorithm.\n");

    /* compute and store length of edges to save some computations */
    length = phgAlloc(g->nedge * sizeof(*length));
    for (i = 0; i < g->nedge; i++)
	length[i] = -1.0;
    for (i = 0; i < g->nroot; i++) {
	e = g->roots + i;
	for (j = 0; j < NEdge; j++) {
	    COORD *v0, *v1;
	    double dx, dy, dz;
	    k = e->edges[j];
	    if (length[k] >= 0)
		continue;
	    v0 = g->verts + e->verts[GetEdgeVertex(j, 0)];
	    v1 = g->verts + e->verts[GetEdgeVertex(j, 1)];
	    dx = (*v0)[0] - (*v1)[0];
	    dy = (*v0)[1] - (*v1)[1];
	    dz = (*v0)[2] - (*v1)[2];
	    length[k] = dx * dx + dy * dy + dz * dz;
	    if (length[k] <= 1e-30)
		bad_edges++;
	}
    }
    if (bad_edges > 0)
	phgWarning("%d zero length edge%s found.\n", bad_edges,
			bad_edges > 1 ? "s" : "");

    for (i = 0; i < g->nroot; i++) {
	e = g->roots + i;
	/* First, find the longest edge */
	k = 0;
	for (j = 1; j < NEdge; j++) {
	    if (comp_edge(length, e, k, j) < 0)
		k = j;
	}
	/* permute the vertices to make the longest edge edge-0 */
	permute_vertices(e, 0, GetEdgeVertex(k, 0));
	permute_vertices(e, 1, GetEdgeVertex(k, 1));
	/* find marked edge (j) of face 0 */
	j = 3;						/* edge 1-2 */
	if (comp_edge(length, e, j, 4) < 0)
	    j = 4;					/* edge 1-3 */
	if (comp_edge(length, e, j, 5) < 0)
	    j = 5;					/* edge 2-3 */
	/* find marked edge (k) of face 1 */
	k = 1;						/* edge 0-2 */
	if (comp_edge(length, e, k, 2) < 0)
	    k = 2;					/* edge 0-3 */
	if (comp_edge(length, e, k, 5) < 0)
	    k = 5;					/* edge 2-3 */
	/* Note: see kossazcky-arnold.dia for numbering of vertices */
	switch (j * 10 + k) {
	    case 31:		/* 1-2 and 0-2, type P = FACE */
		e->type = FACE;
		break;
	    case 32:		/* 1-2 and 0-3, type A = DIAG */
		e->type = DIAGONAL;
		permute_vertices(e, 2, 3);
		break;
	    case 35:		/* 1-2 and 2-3, type M = MIXED */
		e->type = MIXED;
		permute_vertices(e, 0, 1);
		break;
	    case 41:		/* 1-3 and 0-2, type A = DIAG */
		e->type = DIAGONAL;
		break;
	    case 42:		/* 1-3 and 0-3, type P = FACE */
		e->type = FACE;
		permute_vertices(e, 2, 3);
		break;
	    case 45:		/* 1-3 and 2-3, type M = MIXED */
		e->type = MIXED;
		permute_vertices(e, 2, 3);
		permute_vertices(e, 0, 1);
		break;
	    case 51:		/* 2-3 and 0-2, type M = MIXED */
		e->type = MIXED;
		break;
	    case 52:		/* 2-3 and 0-3, type M = MIXED */
		e->type = MIXED;
		permute_vertices(e, 2, 3);
		break;
	    case 55:		/* 2-3 and 2-3, type O = OPPOSITE */
		e->type = OPPOSITE;
		break;
	}
    }

    phgFree(length);

pre_refine:
    /* count number of elements with non standard types (MIXED and OPPOSITE) */
    j = 0;	/* number of simplices of type OPPOSITE */
    k = 0;	/* number of simplices of type MIXED */
    for (i = 0; i < g->nroot; i++) {
	e = g->roots + i;
	switch (e->type) {
	    case OPPOSITE:
		e->mark = 1;
		j++;
		break;
	    case MIXED:
		e->mark = 1;
		k++;
		break;
	    default:
		e->mark = 0;
	}
    }
    phgInfo(2, "non standard types: OPPOSITE=%d, MIXED=%d\n", j, k);


    Return;
}

static void
refine_marked_traverse(GRID *g, ELEMENT *e)
{
    BOOLEAN flag = FALSE;

    if (e->children[0] != NULL) {
	flag = TRUE;
	refine_marked_traverse(g, e->children[0]);
    }

    if (e->children[1] != NULL) {
	flag = TRUE;
	refine_marked_traverse(g, e->children[1]);
    }

    if (flag || e->mark <= 0) {
#if 1 /*DEBUG*/
	/* Note: the function refine_path() should always generate a
	   conforming mesh. But as the present algorithm may have bugs,
	   further refines are done here as a workaround when needed. */
	int i, face, child;
	ELEMENT *e0, *e1;
	for (i = 0; i < NFace; i++) {
	    /* e might be refined during the loop */
	    if (!IsLeaf(e))
		break;
	    if (HasLocalNeighbour(e, i) && e->neighbours[i] == NULL) {
		e0 = e;
		face = i;
		phgWarning("hanging face(%d), doing further refines.\n", face);
		while (TRUE) {
		    if (phgVerbosity > 2)
			phgDumpElement(NULL, e0);
		    e1 = e0;
		    e0 = e0->parent;
		    /* find the matching face no of parent */
		    switch (e0->type) {
			case DIAGONAL:
			case FACE:
			case EDGE:
			    child = (e1->verts[0] == e0->verts[0]) ? 0 : 1;
			    switch (face) {
				case 0:
				    phgError(1, "unexpected error (%s:%d)\n",
						__FILE__, __LINE__);
				case 1:
				case 2:
				    if (child == 1 && e0->type == DIAGONAL)
					face = 3 - face;
				    face++;
				    break;
				case 3:
				    face = 1 - child;
				    break;
			    }
			    break;
			case OPPOSITE:
			    child = (e1->verts[2] == e0->verts[0]) ? 0 : 1;
			    switch (face) {
				case 0:
				case 1:
				    face = 3 - face;
				    break;
				case 2:
				    phgError(1, "unexpected error (%s:%d)\n",
						__FILE__, __LINE__);
				case 3:
				    face = 1 - child;
				    break;
			    }
			    break;
			case MIXED:
			    child = (e1->verts[0] == e0->verts[0]) ? 0 : 1;
			    switch (face) {
				case 0:
				    if (child == 0) {
				      phgError(1, "unexpected error (%s:%d)\n",
						__FILE__, __LINE__);
				    }
				    face = 3;
				    break;
				case 1:
				    face = 2;
				    break;
				case 2:
				    if (child == 1) {
				      phgError(1, "unexpected error (%s:%d)\n",
						__FILE__, __LINE__);
				    }
				    face = 3;
				    break;
				case 3:
				    face = 1 - child;
				    break;
			    }
			    break;
		    }
		    assert(HasLocalNeighbour(e0, face));
		    phgDebug((3, "===== parent %p, face %d\n", e0, face));
		    if (e0->neighbours[face] != NULL)
			break;
		}
		if (phgVerbosity > 2)
		    phgDumpElement(NULL, e0);
		e0 = e0->neighbours[face];
		if (phgVerbosity > 2)
		    phgDumpElement(NULL, e0);
		if (e0->mark < 1)
		    e0->mark = 1;
		refine_marked_traverse(g, e0);
	    }
	}
#endif	/* 1 */
	return;
    }

    refine_path(g, e);

    /* continue with the newly created children */
    refine_marked_traverse(g, e);

    return;
}

static void
update_indices(GRID *g, ELEMENT *e)
/* update indices on new elements */
{
    int fmap[NVert], emap[NEdge];
    ELEMENT *e0, *e1;
    BOOLEAN flag;

    if (IsLeaf(e))
	return;

    if ((e0 = e->children[0]) == NULL) {
	assert(e->children[1] == NULL);
	return;
    }
    e1 = e->children[1];
    assert(e1 != NULL);

    /* update e ==> e0 */

    flag = (GlobalVertexP(g, e->verts[0]) < GlobalVertexP(g, e->verts[1]));
    if ((g->flags & (EDGE_FLAG | FACE_FLAG)))
	phgMapP2C(e, fmap, emap, 0);

    if ((g->flags & ELEM_FLAG) && flag)
	e0->index = e->index;

    if ((g->flags & EDGE_FLAG)) {
	e0->edges[emap[1]] = e->edges[1];
	e0->edges[emap[2]] = e->edges[2];
	e0->edges[emap[5]] = e->edges[5];
	if (flag)
	    e0->edges[emap[0]] = e->edges[0];
    }

    if ((g->flags & FACE_FLAG)) {
	e0->faces[3] = e->faces[1];
	if (flag) {
	    e0->faces[fmap[2]] = e->faces[2];
	    e0->faces[fmap[3]] = e->faces[3];
	}
    }

    update_indices(g, e0);

    /* update e ==> e1 */

    flag = !flag;
    if ((g->flags & (EDGE_FLAG | FACE_FLAG)))
	phgMapP2C(e, fmap, emap, 1);

    if ((g->flags & ELEM_FLAG) && flag)
	e1->index = e->index;

    if ((g->flags & EDGE_FLAG)) {
	e1->edges[emap[3]] = e->edges[3];
	e1->edges[emap[4]] = e->edges[4];
	e1->edges[emap[5]] = e->edges[5];
	if (flag)
	    e1->edges[emap[0]] = e->edges[0];
    }

    if ((g->flags & FACE_FLAG)) {
	e1->faces[3] = e->faces[0];
	if (flag) {
	    e1->faces[fmap[2]] = e->faces[2];
	    e1->faces[fmap[3]] = e->faces[3];
	}
    }

    update_indices(g, e1);

    return;
}

static void
update_nonleaf_indices(GRID *g, ELEMENT *e)
/* update indices on non leaf elements */
{
    int fmap[NVert], emap[NEdge];
    ELEMENT *e0, *e1;
    BOOLEAN flag;

    if (IsLeaf(e))
	return;

    if ((e0 = e->children[0]) != NULL)
	update_nonleaf_indices(g, e0);
    if ((e1 = e->children[1]) != NULL)
	update_nonleaf_indices(g, e1);

    /* note: we only need to update newly created elements, which must have
     * both children nil or non-nil */
    if (e0 == NULL || e1 == NULL)
	return;

    flag = (GlobalVertexP(g, e->verts[0]) < GlobalVertexP(g, e->verts[1]));

    /* TODO: set e->flag = 0 on old elements (before refinement),
     * set e->flag = 1 for new elements, and only update the latter */

    /*if (e0 != NULL)*/ {
#if 0
	if ((g->flags & ELEM_FLAG) && e->index == -1 && flag)
	    e->index = e0->index;
#endif	/* 0 */

	if ((g->flags & (EDGE_FLAG | FACE_FLAG)))
	    phgMapP2C(e, fmap, emap, 0);

	if ((g->flags & EDGE_FLAG)) {
	    if (e->edges[1] == -1)
		e->edges[1] = e0->edges[emap[1]];
	    if (e->edges[2] == -1)
		e->edges[2] = e0->edges[emap[2]];
	    if (e->edges[5] == -1)
		e->edges[5] = e0->edges[emap[5]];
	    if (e->edges[0] == -1 && flag)
		e->edges[0] = e0->edges[emap[0]];
	}

	if ((g->flags & FACE_FLAG)) {
	    if (e->faces[1] == -1)
		e->faces[1] = e0->faces[3];
	    if (flag) {
		if (e->faces[2] == -1)
		    e->faces[2] = e0->faces[fmap[2]];
		if (e->faces[3] == -1)
		    e->faces[3] = e0->faces[fmap[3]];
	    }
	}
    }

    flag = !flag;

    /*if (e1 != NULL)*/ {
#if 0
	if ((g->flags & ELEM_FLAG) && e->index == -1 && flag)
	    e->index = e1->index;
#endif	/* 0 */

	if ((g->flags & (EDGE_FLAG | FACE_FLAG)))
	    phgMapP2C(e, fmap, emap, 1);

	if ((g->flags & EDGE_FLAG)) {
	    if (e->edges[3] == -1)
		e->edges[3] = e1->edges[emap[3]];
	    if (e->edges[4] == -1)
		e->edges[4] = e1->edges[emap[4]];
	    if (e->edges[5] == -1)
		e->edges[5] = e1->edges[emap[5]];
	    if (e->edges[0] == -1 && flag)
		e->edges[0] = e1->edges[emap[0]];
	}

	if ((g->flags & FACE_FLAG)) {
	    if (e->faces[0] == -1)
		e->faces[0] = e1->faces[3];
	    if (flag) {
		if (e->faces[2] == -1)
		    e->faces[2] = e1->faces[fmap[2]];
		if (e->faces[3] == -1)
		    e->faces[3] = e1->faces[fmap[3]];
	    }
	}
    }

    return;
}

void
phgRefineMarkedElements(GRID *g)
/* refines all elements,  each element being refined at least 'mark' times */
{
    INT i;
    INT nelem, nedge, nelem_global;
#if Dim == 3
    INT nface;
#endif	/* Dim == 3 */
    int halo_type = -1;

    FunctionEntry_head;
    ParallelFunction(phgRefineMarkedElements, g != NULL && g->nprocs > 1);
    FunctionEntry_tail;

    if (g == NULL || g->nprocs <= phgRank)
	Return;

    if (g->nprocs > 1) {
	phgPrintf("\n\n"
"****************************************************************************\n"
"* Sorry, parallel refinement is unsupported by this version of refine.o.   *\n"
"*                                                                          *\n"
"* You may consider replacing the file src/refine.o with one from the obj/  *\n"
"* directory, please read obj/ReadMe.txt for more details.                  *\n"
"*                                                                          *\n"
"* Abort.                                                                   *\n"
"****************************************************************************\n"
	"\n");
#if USE_MPI
	MPI_Abort(phgComm, 1);
#endif	/* USE_MPI */
	exit(1);
    }

    if (g->halo != NULL) {
	halo_type = g->halo->type;
	assert(halo_type >= 0);
	phgRemoveHalo(g);
	g->serial_no -= 1000 * 1000;	/* restore serial_no */
    }
    g->serial_no++;

    g->alien = NULL;
    phgFree(g->alien_map);
    g->alien_map = NULL;

    phgFree(g->Fmap_elems);
    g->Fmap_elems = NULL;
    phgFree(g->Fmap_faces);
    g->Fmap_faces = NULL;

    nelem_global = g->nelem_global;
    nelem = g->nelem;
    nedge = g->nedge;
#if Dim == 3
    nface = g->nface;
#endif	/* Dim == 3 */

    phgDofInterpolate(g, FALSE);	/* save state of current mesh */

    g_verts_allocated = nvert_save = g->nvert;

    g_ = g;


    max_recur = 0;
    too_close = 0;

    for (i = 0; i < g->nroot; i++)
	refine_marked_traverse(g, g->roots + i);

    if (too_close > 0)
	phgWarning("some vertices got too close during refinement.\n");


    /* update L2Gmap_vert and reset tables */
    new_vertex(g, NULL);

    nvert_save = g->nvert;


    if ((g->flags & (EDGE_FLAG | FACE_FLAG | ELEM_FLAG))) {
	/* update edge/face indices of new elements (vertex
	 * are always up to date during refinement) */
	for (i = 0; i < nelem; i++)
	    if (g->elems[i] != NULL)
		update_indices(g, g->elems[i]);
    }

    /* assign indices to new edges */
    if ((g->flags & EDGE_FLAG)) {
	phgUpdateEdges(g);
	phgInfo(2, "+%d (%d) edges\n", g->nedge - nedge, g->nedge);
    }

#if Dim == 3
    /* assign indices to new faces */
    if ((g->flags & FACE_FLAG)) {
	phgUpdateFaces(g);
	phgInfo(2, "+%d (%d) faces\n", g->nface - nface, g->nface);
    }
#endif	/* Dim == 3 */

    /* assign indices to new elements */
    if ((g->flags & ELEM_FLAG)) {
	phgUpdateElementIndices(g, nelem, nelem_global);
    }

    if ((g->flags & (EDGE_FLAG | FACE_FLAG))) {
	/* update edge/face indices of non leaf element (vertex/element
	 * indices are already up to date) */
	for (i = 0; i < g->nroot; i++)
	    update_nonleaf_indices(g, g->roots + i);
    }

    phgInfo(2, "max recursion = %d\n", max_recur);
    phgInfo(2, "+%d (%d) elements.\n", g->nelem - nelem, g->nleaf);

    /* free extra memory after the end of g->verts */
    if (g->nvert < g_verts_allocated)
	g->verts = phgRealloc_(g->verts, sizeof(*(g->verts)) * g->nvert,
				sizeof(*(g->verts)) * g->nvert);
    phgInfo(1, "%d vertices, %d elements, %d edges.\n",
		g->nvert, g->nleaf, g->nedge);

    phgUpdateBoundaryTypes(g);

    phgDofInterpolate(g, TRUE);		/* do interpolation */

#if 0
    if (halo_type >= 0)
	phgSetupHalo(g, halo_type);
#endif

#if ALLOW_CURVED_BOUNDARY
    if (g->bdry_funcs != NULL) {
	/* update volume */
	FLOAT d = 0.;
	ELEMENT *e;
	if ((g->flags & GEOM_FLAG)) {
	    ForAllElements(g, e)
		d += phgGeomGetVolume(g, e);
	}
	else {
	    FLOAT a;
	    ForAllElements(g, e) {
		_phg_compute_elem_data(g, e, &a, NULL);
		d += a;
	    }
	}
#if USE_MPI
	MPI_Allreduce(&d, &g->volume, 1, PHG_MPI_FLOAT, PHG_SUM, g->comm);
#else	/* USE_MPI */
	g->volume = d;
#endif	/* USE_MPI */
    }
#endif	/* ALLOW_CURVED_BOUNDARY */

    PrintTime(1);

    Return;
}

/*--------------- The functions below are for testing purpose -------------*/

int refine_which_ = -1;	/* -1: refine all, >0: refine given submesh */
static int refine_level = 1;	/* level for phgRefineAllElements */

static BOOLEAN
refine_all_callback CB_ARGS(e)
{
    if (g_->nprocs <= 1 || refine_which_ < 0 || g_->rank == refine_which_)
	e->mark = refine_level;
    return TRUE;
}

void
phgRefineAllElements_(GRID *g)
/* refines all elements, each element is at leas refined once */
{
    ParallelFunction(phgRefineAllElements_, g != NULL && g->nprocs > 1);

    if (g == NULL || g->nprocs <= phgRank)
	return;


    g_ = g;
    phgTraverseElements(g, refine_all_callback);

    phgRefineMarkedElements(g);    

    return;
}

void
phgRefineAllElements(GRID *g, int level)
{
    if ((refine_level = level) <= 0)
	return;
    phgRefineAllElements_(g);
}

/* threshold = RAND_MAX * (percentage of marked elements) */
static long int threshold;
static const char *percent_or_number;

static BOOLEAN
random_refine_callback CB_ARGS(e)
{
    if (g_->nprocs > 1 && refine_which_ >= 0 && g_->rank != refine_which_)
	e->mark = 0;
    else
	e->mark = (rand()<=threshold) ? 1 : 0;
    return TRUE;
}

#include <ctype.h>

void
phgRefineRandom_(GRID *g)
/* refines randomly selected elements, the second argument gives
 * approximate number of elements to refine: "30" means 30 elements
 * while "30%" means 30% of the elements.
 *
 * This function is mainly for testing the refinement code */
{
    double d;
    char *end;

    ParallelFunction(phgRefineRandom_, g != NULL && g->nprocs > 1);

    if (g == NULL || g->nprocs <= phgRank)
	return;


    d = strtod(percent_or_number, &end);
    while (isspace(*(BYTE *)end))
	end++;
    d /= (*end == '%') ? 100 : g->nleaf;
    threshold = (long int)(RAND_MAX * d + 0.5);

    g_ = g;
    phgTraverseElements(g, random_refine_callback);
    phgRefineMarkedElements(g);

    return;
}

void
phgRefineRandomElements(GRID *g, const char *s)
{
    percent_or_number = s;
    phgRefineRandom_(g);
}

static double x_center, y_center, z_center, ra, rb, rc;

static double
func(double x, double y, double z)
/* function defining the surface to refine */
{
    /* an ellipsoid */
    x = (x - x_center) / ra;
    y = (y - y_center) / rb;
    z = (z - z_center) / rc;
    return x * x + y * y + z * z - _F(1.2) * _F(1.2);
}

static BOOLEAN
refine_surface_callback CB_ARGS(e)
{
    double f, min, max;
    COORD *c;
    int i;

    if (g_->nprocs > 1 && refine_which_ >= 0 && g_->rank != refine_which_) {
	e->mark = 0;
	return TRUE;
    }

    /* Note: the code below won't work if the entire surface is contained
       in the tetrahedra. Should implement more advanced tests. */
    c = g_->verts + e->verts[0];
    max = min = func(*(c)[0], (*c)[1], (*c)[2]);
    for (i = 1; i < NVert; i++) {
	c = g_->verts + e->verts[i];
	f = func(*(c)[0], (*c)[1], (*c)[2]);
	if (min > f)
	    min = f;
	if (max < f)
	    max = f;
    }
    e->mark = (min <= 0.0 && max >= 0.0) ? 1 : 0;
    return TRUE;
}

void
phgRefineSurface(GRID *g)
/* refines all elements, each element is at leas refined once */
{
    FunctionEntry_head;
    ParallelFunction(phgRefineSurface, g != NULL && g->nprocs > 1);
    FunctionEntry_tail;

    if (g == NULL || g->nprocs <= phgRank)
	Return;

    x_center = (g->bbox[0][0] + g->bbox[1][0]) * 0.5;
    ra = (g->bbox[1][0] - g->bbox[0][0]) * 0.5;
    y_center = (g->bbox[0][1] + g->bbox[1][1]) * 0.5;
    rb = (g->bbox[1][1] - g->bbox[0][1]) * 0.5;
    z_center = (g->bbox[0][2] + g->bbox[1][2]) * 0.5;
    rc = (g->bbox[1][2] - g->bbox[0][2]) * 0.5;

    phgDebug((2, "(x0,y0,z0)=(%lf,%lf,%lf)\n",
			(double)x_center, (double)y_center, (double)z_center));
    phgDebug((2, "(ra,rb,rc)=(%lf,%lf,%lf)\n",
			(double)ra, (double)rb, (double)rc));

    g_ = g;
    phgTraverseElements(g, refine_surface_callback);

    phgRefineMarkedElements(g);    

    Return;
}
