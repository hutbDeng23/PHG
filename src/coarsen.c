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

/* $Id: coarsen.c,v 1.223 2022/03/24 06:31:27 zlb Exp $
 *
 * The serial algorithm:
 *
 * On entry, for any leaf element e, -(e->mark) is the max number of times e
 * is allowed to be coarsened.
 * 
 *   1. The mesh is traversed using function set_mark(), in which
 *	    e->mark = max(e->children[0]->mark, e->children[1]->mark) + 1
 *	is computed for all non leaf elements, e->flag is set to 0 for all
 *	elements, neighbour links and boundary types are updated for non
 *	leaf elements.
 *
 *   2. check_patch() is called to check all patches around the refinement
 *	edges of non leaf elements. If one of the element on a patch has
 *	e->mark > 0, then the patch can't be coarsened, e->flag is recursively
 *	set to CHECK_FAILED for all elements in the patch, as well as those
 *	on the refinement patches of their ancestors.
 *
 *   3. On each branch, prune the children of the top-most element having
 *	e->flag != CHECK_FAILED && e->mark <= 0.
 *
 * The above algorithm is proved to produce a conforming mesh.
 */

#include "phg.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#define ONE_LEVEL_COARSENING	0 /* level by level coarsening if nprocs > 1 */

#define CHECK_FAILED	1
#define CHECK_PASSED	2
#define REMOVED		255

/* Note: GetNeighbour() relies on the REMOTE bit which is not usable for non
 * leaf elements */
#define GetNeighbour0(e, v) (IsLeaf(e) ? GetNeighbour(e, v) : e->neighbours[v])

static GRID *g_;

static void
set_mark(GRID *g, ELEMENT *e)
/* computes e->mark for non leaf elements, sets e->flag=0 for all elements */
{
    ELEMENT *e0, *e1;
    int mark;

    if (e == NULL)
	return;

    e->flag = 0;

    e0 = e->children[0];
    e1 = e->children[1];

    if (e0 == NULL && e1 == NULL)
	return;

    set_mark(g, e0);
    set_mark(g, e1);

    if (e0 != NULL) {
	e->mark = e0->mark + 1;
	if (e1 != NULL) {
	    e->mark = (e->mark >= (mark = e1->mark + 1)) ? e->mark : mark;
	}
    }
    else {	/* e1 != NULL */
	e->mark = e1->mark + 1;
    }

    return;
}

static void
free_children(ELEMENT *e)
{
#if ONE_LEVEL_COARSENING
    int mark;
#endif	/* ONE_LEVEL_COARSENING */
    ELEMENT *e0 = e->children[0], *e1 = e->children[1];

    if (e0 == NULL && e1 == NULL)
	return;

    if (e0 != NULL)
	free_children(e0);

    if (e1 != NULL)
	free_children(e1);

#if ONE_LEVEL_COARSENING
    /* update e->mark */
    if (e0 != NULL) {
	e->mark = e0->mark + 1;
	if (e1 != NULL) {
	    e->mark = (e->mark >= (mark = e1->mark + 1)) ? e->mark : mark;
	}
    }
    else {	/* e1 != NULL */
	e->mark = e1->mark + 1;
    }
#endif	/* ONE_LEVEL_COARSENING */

    if (e0 != NULL) {
	phgFreeElement(&e0);
	e->children[0] = NULL;
    }

    if (e1 != NULL) {
	phgFreeElement(&e1);
	e->children[1] = NULL;
    }
}

#if DEBUG
static int
get_patch(GRID *g,  ELEMENT *e, ELEMENT ***p)
/* returns the patch composed of all elements having edge 0 of 'e' as their
 * edge 0. Note: phgGetPatch() here requires neighbours links of non-leaf
 * elements, which are updated by the function set_mark() */
{
    int i, count;
    INT u0, u1, v0, v1;

    v0 = GlobalVertexP(g, e->verts[0]);
    v1 = GlobalVertexP(g, e->verts[1]);
    count = phgGetPatch(g, e, 0, 1, p);
    for (i = 0; i < count; i++) {
	e = (*p)[i];
	u0 = GlobalVertexP(g, e->verts[0]);
	u1 = GlobalVertexP(g, e->verts[1]);
	assert((v0 == u0 && v1 == u1) || (v0 == u1 && v1 == u0));
    }

    return count;
}
#else	/* DEBUG */
# define get_patch(g, e, p)	phgGetPatch(g, e, 0, 1, p)
#endif	/* DEBUG */

static void
set_uncoarsen(GRID *g, ELEMENT *e)
{
    ELEMENT **p;
    int i, count;

    if (e == NULL || e->flag == CHECK_FAILED)
	return;

    count = get_patch(g, e, &p);
    for (i = 0; i < count; i++)
	p[i]->flag = CHECK_FAILED;
    for (i = 0; i < count; i++)
	set_uncoarsen(g, p[i]->parent);
    phgFree(p);
}

static void
check_patch(GRID *g, ELEMENT *e)
{
    ELEMENT **p;
    int i, count;

    if (e == NULL || IsLeaf(e))
	return;

    check_patch(g, e->children[0]);
    check_patch(g, e->children[1]);

    if (e->flag == CHECK_FAILED || e->flag == CHECK_PASSED)
	return;

    count = get_patch(g, e, &p);
    for (i = 0; i < count; i++) {
	if (p[i]->mark > 0)
	    break;
	e->flag = CHECK_PASSED;
    }

    if (i >= count) {
	/* all elements on the patch passed the test */
	phgFree(p);
	return;
    }

    /* the patch can't be coarsened, set e->flag == CHECK_FAILED for all
     * elements on the patch and all their ancestors */
    for (i = 0; i < count; i++)
	p[i]->flag = CHECK_FAILED;
    for (i = 0; i < count; i++)
	set_uncoarsen(g, p[i]->parent);
    phgFree(p);
}

/* Algorithm for interpolating DOFs:
 *   0) calculate and save local h-p orders on elements in the coarsened mesh,
 *	the h-p order on an element is the maximum order of all its descendants.
 *   1) update saved h-p orders in remove_duplicates().
 *   2) save removed branches if dof_interp_flag == TRUE.
 *   3) remove removed elements, reconstruct new mesh as usual, pack entries
 *	in hp_orders[] in pack_indices().
 *   4) do DOF interpolation using saved branches and element orders. */

static GRID *g_save = NULL;		/* for saving removed branches */
static BOOLEAN dof_interp_flag = FALSE;	/* whether have DOFs to interpolate */
static int nhp;				/* number of HP_TYPE structs */
static CHAR *hp_orders;			/* h-p orders, size = g->nelem * nhp */
static INT nleaf_new, elem_index;	/* index of current element */

static void
save_hp_orders1(ELEMENT *e)
{
    int ii;

    for (ii = 0; ii < nhp; ii++) {
	if (hp_orders[g_->nelem * ii + elem_index] <
		g_->hp[ii]->elem_order[elem_index])
	    hp_orders[g_->nelem * ii + elem_index] =
				g_->hp[ii]->elem_order[elem_index];
    }

    if (e->children[0] != NULL)
	save_hp_orders1(e->children[0]);
    if (e->children[1] != NULL)
	save_hp_orders1(e->children[1]);
}

static void
save_hp_orders0(ELEMENT *e)
{
    if ((e->mark > 0 || e->flag == CHECK_FAILED) && !IsLeaf(e)) {
	if (e->children[0] != NULL)
	    save_hp_orders0(e->children[0]);
	if (e->children[1] != NULL)
	    save_hp_orders0(e->children[1]);
	return;
    }

    nleaf_new++;		/* increment # of elements on coarsened mesh */
    if (nhp > 0) {
	elem_index = e->index;
	save_hp_orders1(e);
    }

    return;
}

static void
save_hp_orders(GRID *g)
{
    INT i;

    nleaf_new = 0;	/* count number of elements on coarsened mesh */
    for (i = 0; i < g->nroot; i++)
	save_hp_orders0(g->roots + i);
}


static void
do_coarsen(GRID *g, ELEMENT *e, INT *ntree, INT *nleaf)
/* coarsens the mesh, and saves pruned branches to g_save for DOF interpolation
 * if dof_interpolate == TRUE */
{
    ELEMENT *ee;

    if (e == NULL)
	return;

    if (e->flag == REMOVED) {
	/* Note: e must be a root element in this case */
	assert(e->parent == NULL);
	free_children(e);
	return;
    }

    (*ntree)++;

    if ((e->mark <= 0 && e->flag != CHECK_FAILED) || IsLeaf(e)) {
	/* e is a leaf on the new (coarsened mesh), prune its children */
	if (!dof_interp_flag) {
	    /* prune children of e */
	    free_children(e);
	    (*nleaf)++;
	}
	else {
	    /* save the branch to g_save */
	    g_save->roots[*nleaf] = *e;
	    /* points parent of the saved element to the original one */
	    g_save->roots[(*nleaf)++].parent = e;
	    /* remove children */
	    e->children[0] = e->children[1] = NULL;
	}
	return;
    }

    /* e is not a leaf on the new mesh */
    if ((ee = e->children[0]) != NULL) {
	if (ee->flag == REMOVED) {
	    free_children(ee);
	    phgFreeElement(&ee);
	    e->children[0] = NULL;
	}
	else {
	    do_coarsen(g, ee, ntree, nleaf);
	}
    }

    if ((ee = e->children[1]) != NULL) {
	if (ee->flag == REMOVED) {
	    free_children(ee);
	    phgFreeElement(&ee);
	    e->children[1] = NULL;
	}
	else {
	    do_coarsen(g, ee, ntree, nleaf);
	}
    }

    return;
}

static INT *map;
static size_t map_local_offset, map_local_size, map_local_max;

static BOOLEAN
map_collect_cb CB_ARGS(e)
{
    INT i;
    INT *local_map = (INT *)(((BYTE *)e) + map_local_offset);

    /* make phgUpdateBoundaryTypes() recompute ordering */
    e->ordering = 0;

    for (i = 0; i < map_local_size; i++) {
	assert(local_map[i] >=0 && local_map[i] < map_local_max);
	map[local_map[i]] = 1;
    }

    return TRUE;
}

static BOOLEAN
map_update_cb CB_ARGS(e)
{
    INT i;
    INT *local_map = (INT *)(((BYTE *)e) + map_local_offset);

    for (i = 0; i < map_local_size; i++) {
	local_map[i] = map[local_map[i]];
	assert(local_map[i] >= 0);
    }

    return TRUE;
}

static BOOLEAN
pack_indices(GRID *g, GTYPE type)
/* updates indices and L2Gmap, `type' may be VERTEX, EDGE, FACE, or VOLUME */
{
    ELEMENT dummy;
    void (*traverse)(GRID *g, BOOLEAN (*cb) CB_ARGS(e)) = NULL;
    INT i, j, newcount, *count = NULL, *count_global = NULL;

    switch (type) {
	case VERTEX:
	    traverse = phgTraverseAllElements;
	    count = &g->nvert;
	    count_global = &g->nvert_global;
	    map_local_offset = (BYTE *)(dummy.verts) - ((BYTE *)&dummy);
	    map_local_size = NVert;
	    map_local_max = g->nvert;
	    break;
	case EDGE:
	    traverse = phgTraverseAllElements;
	    count = &g->nedge;
	    count_global = &g->nedge_global;
	    map_local_offset = (BYTE *)(dummy.edges) - ((BYTE *)&dummy);
	    map_local_size = NEdge;
	    map_local_max = g->nedge;
	    break;
	case FACE:
	    traverse = phgTraverseAllElements;
	    count = &g->nface;
	    count_global = &g->nface_global;
	    map_local_offset = (BYTE *)(dummy.faces) - ((BYTE *)&dummy);
	    map_local_size = NFace;
	    map_local_max = g->nface;
	    break;
	case VOLUME:
	    traverse = phgTraverseAllElements;
	    count = &g->nelem;
	    count_global = &g->nelem_global;
	    map_local_offset = (BYTE *)(&dummy.index) - ((BYTE *)&dummy);
	    map_local_size = 1;
	    map_local_max = g->nelem;
	    break;
    }

    /* collect referenced indices */
    map = phgCalloc(map_local_max, sizeof(*map));
    (*traverse)(g, map_collect_cb);

    newcount = -1;
    for (i = 0; i < *count; i++) {
	if (map[i] != 0)
	    map[i] = ++newcount;
	else
	    map[i] = newcount;
    }
    newcount++;

    if (newcount == *count) {
	/* no change */
	phgFree(map);
	return FALSE;
    }

    /* assign new indices */
    (*traverse)(g, map_update_cb);

    if (type == VERTEX) {
	/* update g->verts[] and g->period->L2Gmap_vert */
	COORD *verts;
	if (g_save != NULL)
	    verts = phgAlloc(newcount * sizeof(*verts));
	else
	    verts = g->verts;
	for (i = j = 0; i < *count; i++) {
	    if (map[i] >= 0 && (i == 0 || map[i] != map[i - 1])) {
		if (i != j) {
		    memcpy(verts + j, g->verts + i, sizeof(*verts));
		    if (g->period != NULL)
			g->period->L2Gmap_vert[j] = g->period->L2Gmap_vert[i];
		}
		else if (verts != g->verts) {
		    memcpy(verts + j, g->verts + i, sizeof(*verts));
		}
		j++;
	    }
	}
	if (g_save != NULL) {
	    g_save->verts = g->verts;
	    g->verts = verts;
	}
	else {
	    assert(verts == g->verts);
	    g->verts = phgRealloc_(verts, newcount * sizeof(*verts),
					  (*count) * sizeof(*verts));
	}
	if (g->period != NULL)
	    g->period->L2Gmap_vert = phgRealloc_(g->period->L2Gmap_vert,
				newcount * sizeof(*g->period->L2Gmap_vert),
				(*count) * sizeof(*g->period->L2Gmap_vert));
    }


    *count = newcount;
    if (g->nprocs == 1)
	*count_global = newcount;

    phgFree(map);
    return TRUE;
}

static struct {
    ELEMENT	*e0, *e1;
    BYTE	v0, v1;
} *faces, *pf;

static BOOLEAN
update_neighbours_callback CB_ARGS(e)
{
    int v;

    for (v = 0; v < NFace; v++) {
	pf = faces + e->faces[v];
	if (pf->e0 == NULL) {
	    pf->e0 = e;
	    pf->v0 = v;
	}
	else {
	    if (pf->e1 != NULL) {
		phgDumpElement(g_, pf->e0);
		phgDumpElement(g_, pf->e1);
		phgDumpElement(g_, e);
		phgError(0, "%s:%d, duplicate element.\n", __FILE__, __LINE__);
	    }
	    pf->e1 = e;
	    pf->v1 = v;
	}
	/* clear neighbour */
	e->neighbours[v] = NULL;
    }

    return TRUE;
}

static void
update_neighbours(GRID *g)
{
    ELEMENT *e0, *e1;
    INT i;
    int v0, v1;

    if (!(g->flags & FACE_FLAG)) {
	phgUpdateNeighbours(g);
	return;
    }

    faces = phgCalloc(g->nface, sizeof(*faces));
    phgTraverseElements(g, update_neighbours_callback);
    for (i = 0, pf = faces; i < g->nface; i++, pf++) {
	if ((e0 = pf->e0) == NULL)
	    continue;
	if ((e1 = pf->e1) == NULL)
	    continue;
	e0->bound_type[v0 = pf->v0] &= ~REMOTE;
	e0->neighbours[v0] = e1;
	e1->bound_type[v1 = pf->v1] &= ~REMOTE;
	e1->neighbours[v1] = e0;
    }
    phgFree(faces);

    return;
}

static void
pack_period_L2Gmap(GRID *g)
{
    INT i, j, *map;

    assert(g->nprocs == 1 && g->period != NULL);

    map = phgCalloc(g->period->nvert_global, sizeof(*map));

    for (i = 0; i < g->nvert; i++)
	map[g->period->L2Gmap_vert[i]] = 1;

    for (i = j = 0; i < g->period->nvert_global; i++)
	if (map[i] != 0)
	    map[i] = j++;

    g->period->nvert_global = j;

    for (i = 0; i < g->nvert; i++)
	g->period->L2Gmap_vert[i] = map[g->period->L2Gmap_vert[i]];

    phgFree(map);
}


/*----------------------- DOF interpolation ------------------------*/

static DOF *dof;
static HP_TYPE *hp_old, *hp_new;
static FLOAT *dof_data_new;
static BYTE *flags;
static int pass;	/* only used when userfunc != DofInterpolation */
static ELEMENT *e_old, *e_new;

static void
dof_interp_recursive(ELEMENT *e, FLOAT **parent_data)
/* Note: e is an element in g_save,
 *	 dof->g and dof->hp match the new mesh,
 *	 dof->data and dof->data_xxxx point to old data
 *
 * Care should be taken when referring to old data locations, e.g.,
 * 	DofEdgeData(dof, e->edges[0])
 * is correct when dof is non h-p, but is wrong when dof is h-p, for the latter
 * case 'dof->data_edge + dof->dim * hp_old->edge_index[e->edges[0]]'
 * should be used.
 */
{
    ELEMENT *e0 = e->children[0], *e1 = e->children[1];
    DOF_USER_FUNC f = dof->userfunc;
    DOF_TYPE *type = (DofIsHP(dof) ? dof->hp->info->types[dof->hp->max_order] :
				     dof->type);
    FLOAT *p, *tmp, *data[15];
    INT i, k;
    int vmap[NVert], emap[NEdge];
    size_t size;
    static int level = 0;

    if (e0 == NULL && e1 == NULL) {
	/* leaf element, simply copy data */
	if (f != DofInterpolation && pass > 0)
	    return;
	if (level > 0) {
	    /* pruned (intermediate) element, copy all data */
	    assert(f == DofInterpolation);
	    assert(!DofIsHP(dof));	/* must be a non h-p DOF */
	    if (dof->count_vert > 0) {
		size = dof->count_vert * sizeof(FLOAT);
		memcpy(parent_data[0], DofVertexData(dof, e->verts[0]), size);
		memcpy(parent_data[1], DofVertexData(dof, e->verts[1]), size);
		memcpy(parent_data[2], DofVertexData(dof, e->verts[2]), size);
		memcpy(parent_data[3], DofVertexData(dof, e->verts[3]), size);
	    }
	    if (dof->count_edge > 0) {
		size = dof->count_edge * sizeof(FLOAT);
		memcpy(parent_data[4], DofEdgeData(dof, e->edges[0]), size);
		memcpy(parent_data[5], DofEdgeData(dof, e->edges[1]), size);
		memcpy(parent_data[6], DofEdgeData(dof, e->edges[2]), size);
		memcpy(parent_data[7], DofEdgeData(dof, e->edges[3]), size);
		memcpy(parent_data[8], DofEdgeData(dof, e->edges[4]), size);
		memcpy(parent_data[9], DofEdgeData(dof, e->edges[5]), size);
	    }
	    if (dof->count_face > 0) {
		size = dof->count_face * sizeof(FLOAT);
		memcpy(parent_data[10], DofFaceData(dof, e->faces[0]), size);
		memcpy(parent_data[11], DofFaceData(dof, e->faces[1]), size);
		memcpy(parent_data[12], DofFaceData(dof, e->faces[2]), size);
		memcpy(parent_data[13], DofFaceData(dof, e->faces[3]), size);
	    }
	}
	else {		/* if (level > 0) */
	    /* leaf element on both old and new meshes, copy old data */
	    if (dof->count_vert > 0) {
		if (!DofIsHP(dof)) {
		    /* non h-p DOF */
		    size = dof->count_vert * sizeof(FLOAT);
		    if (!(flags[i = e->verts[0]] & 1)) {
			flags[i] |= 1;
			memcpy(parent_data[0], DofVertexData(dof, i), size);
		    }
		    if (!(flags[i = e->verts[1]] & 1)) {
			flags[i] |= 1;
			memcpy(parent_data[1], DofVertexData(dof, i), size);
		    }
		    if (!(flags[i = e->verts[2]] & 1)) {
			flags[i] |= 1;
			memcpy(parent_data[2], DofVertexData(dof, i), size);
		    }
		    if (!(flags[i = e->verts[3]] & 1)) {
			flags[i] |= 1;
			memcpy(parent_data[3], DofVertexData(dof, i), size);
		    }
		}
		else {
		    /* h-p DOF */
		    if (!(flags[i = e->verts[0]] & 1)) {
			flags[i] |= 1;
			size = sizeof(FLOAT) * dof->dim * (
					hp_old->vert_index[i + 1] -
					(k = hp_old->vert_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->vert_index[e_new->verts[0] + 1] -
				hp_new->vert_index[e_new->verts[0]]));
			if (size > 0)
			    memcpy(parent_data[0],
				   dof->data_vert + k * dof->dim, size);
		    }
		    if (!(flags[i = e->verts[1]] & 1)) {
			flags[i] |= 1;
			size = sizeof(FLOAT) * dof->dim * (
					hp_old->vert_index[i + 1] -
					(k = hp_old->vert_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->vert_index[e_new->verts[1] + 1] -
				hp_new->vert_index[e_new->verts[1]]));
			if (size > 0)
			    memcpy(parent_data[1],
				   dof->data_vert + k * dof->dim, size);
		    }
		    if (!(flags[i = e->verts[2]] & 1)) {
			flags[i] |= 1;
			size = sizeof(FLOAT) * dof->dim * (
					hp_old->vert_index[i + 1] -
					(k = hp_old->vert_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->vert_index[e_new->verts[2] + 1] -
				hp_new->vert_index[e_new->verts[2]]));
			if (size > 0)
			    memcpy(parent_data[2],
				   dof->data_vert + k * dof->dim, size);
		    }
		    if (!(flags[i = e->verts[3]] & 1)) {
			flags[i] |= 1;
			size = sizeof(FLOAT) * dof->dim * (
					hp_old->vert_index[i + 1] -
					(k = hp_old->vert_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->vert_index[e_new->verts[3] + 1] -
				hp_new->vert_index[e_new->verts[3]]));
			if (size > 0)
			    memcpy(parent_data[3],
				   dof->data_vert + k * dof->dim, size);
		    }
		}
	    }
	    if (dof->count_edge > 0) {
		if (!DofIsHP(dof)) {
		    /* non h-p DOF */
		    size = dof->count_edge * sizeof(FLOAT);
		    if (!(flags[i = e->edges[0]] & 2)) {
			flags[i] |= 2;
			memcpy(parent_data[4], DofEdgeData(dof, i), size);
		    }
		    if (!(flags[i = e->edges[1]] & 2)) {
			flags[i] |= 2;
			memcpy(parent_data[5], DofEdgeData(dof, i), size);
		    }
		    if (!(flags[i = e->edges[2]] & 2)) {
			flags[i] |= 2;
			memcpy(parent_data[6], DofEdgeData(dof, i), size);
		    }
		    if (!(flags[i = e->edges[3]] & 2)) {
			flags[i] |= 2;
			memcpy(parent_data[7], DofEdgeData(dof, i), size);
		    }
		    if (!(flags[i = e->edges[4]] & 2)) {
			flags[i] |= 2;
			memcpy(parent_data[8], DofEdgeData(dof, i), size);
		    }
		    if (!(flags[i = e->edges[5]] & 2)) {
			flags[i] |= 2;
			memcpy(parent_data[9], DofEdgeData(dof, i), size);
		    }
		}
		else {
		    /* h-p DOF */
		    if (!(flags[i = e->edges[0]] & 2)) {
			flags[i] |= 2;
			size = sizeof(FLOAT) * dof->dim * (
					hp_old->edge_index[i + 1] -
					(k = hp_old->edge_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->edge_index[e_new->edges[0] + 1] -
				hp_new->edge_index[e_new->edges[0]]));
			if (size > 0)
			    memcpy(parent_data[4],
				   dof->data_edge + k * dof->dim, size);
		    }
		    if (!(flags[i = e->edges[1]] & 2)) {
			flags[i] |= 2;
			size = sizeof(FLOAT) * dof->dim * (
					hp_old->edge_index[i + 1] -
					(k = hp_old->edge_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->edge_index[e_new->edges[1] + 1] -
				hp_new->edge_index[e_new->edges[1]]));
			if (size > 0)
			    memcpy(parent_data[5],
				   dof->data_edge + k * dof->dim, size);
		    }
		    if (!(flags[i = e->edges[2]] & 2)) {
			flags[i] |= 2;
			size = sizeof(FLOAT) * dof->dim * (
						hp_old->edge_index[i + 1] -
						(k = hp_old->edge_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->edge_index[e_new->edges[2] + 1] -
				hp_new->edge_index[e_new->edges[2]]));
			if (size > 0)
			    memcpy(parent_data[6],
				   dof->data_edge + k * dof->dim, size);
		    }
		    if (!(flags[i = e->edges[3]] & 2)) {
			flags[i] |= 2;
			size = sizeof(FLOAT) * dof->dim * (
						hp_old->edge_index[i + 1] -
						(k = hp_old->edge_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->edge_index[e_new->edges[3] + 1] -
				hp_new->edge_index[e_new->edges[3]]));
			if (size > 0)
			    memcpy(parent_data[7],
				   dof->data_edge + k * dof->dim, size);
		    }
		    if (!(flags[i = e->edges[4]] & 2)) {
			flags[i] |= 2;
			size = sizeof(FLOAT) * dof->dim * (
						hp_old->edge_index[i + 1] -
						(k = hp_old->edge_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->edge_index[e_new->edges[4] + 1] -
				hp_new->edge_index[e_new->edges[4]]));
			if (size > 0)
			    memcpy(parent_data[8],
				   dof->data_edge + k * dof->dim, size);
		    }
		    if (!(flags[i = e->edges[5]] & 2)) {
			flags[i] |= 2;
			size = sizeof(FLOAT) * dof->dim * (
						hp_old->edge_index[i + 1] -
						(k = hp_old->edge_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->edge_index[e_new->edges[5] + 1] -
				hp_new->edge_index[e_new->edges[5]]));
			if (size > 0)
			    memcpy(parent_data[9],
				   dof->data_edge + k * dof->dim, size);
		    }
		}
	    }
	    if (dof->count_face > 0) {
		if (!DofIsHP(dof)) {
		    /* non h-p DOF */
		    size = dof->count_face * sizeof(FLOAT);
		    if (!(flags[i = e->faces[0]] & 4)) {
			flags[i] |= 4;
			memcpy(parent_data[10], DofFaceData(dof, i), size);
		    }
		    if (!(flags[i = e->faces[1]] & 4)) {
			flags[i] |= 4;
			memcpy(parent_data[11], DofFaceData(dof, i), size);
		    }
		    if (!(flags[i = e->faces[2]] & 4)) {
			flags[i] |= 4;
			memcpy(parent_data[12], DofFaceData(dof, i), size);
		    }
		    if (!(flags[i = e->faces[3]] & 4)) {
			flags[i] |= 4;
			memcpy(parent_data[13], DofFaceData(dof, i), size);
		    }
		}
		else {
		    /* h-p DOF */
		    if (!(flags[i = e->faces[0]] & 4)) {
			flags[i] |= 4;
			size = sizeof(FLOAT) * dof->dim * (
						hp_old->face_index[i + 1] -
						(k = hp_old->face_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->face_index[e_new->faces[0] + 1] -
				hp_new->face_index[e_new->faces[0]]));
			if (size > 0)
			    memcpy(parent_data[10],
				   dof->data_face + k * dof->dim, size);
		    }
		    if (!(flags[i = e->faces[1]] & 4)) {
			flags[i] |= 4;
			size = sizeof(FLOAT) * dof->dim * (
						hp_old->face_index[i + 1] -
						(k = hp_old->face_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->face_index[e_new->faces[1] + 1] -
				hp_new->face_index[e_new->faces[1]]));
			if (size > 0)
			    memcpy(parent_data[11],
				   dof->data_face + k * dof->dim, size);
		    }
		    if (!(flags[i = e->faces[2]] & 4)) {
			flags[i] |= 4;
			size = sizeof(FLOAT) * dof->dim * (
						hp_old->face_index[i + 1] -
						(k = hp_old->face_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->face_index[e_new->faces[2] + 1] -
				hp_new->face_index[e_new->faces[2]]));
			if (size > 0)
			    memcpy(parent_data[12],
				   dof->data_face + k * dof->dim, size);
		    }
		    if (!(flags[i = e->faces[3]] & 4)) {
			flags[i] |= 4;
			size = sizeof(FLOAT) * dof->dim * (
						hp_old->face_index[i + 1] -
						(k = hp_old->face_index[i]));
			assert(size / sizeof(FLOAT) == dof->dim * (
				hp_new->face_index[e_new->faces[3] + 1] -
				hp_new->face_index[e_new->faces[3]]));
			if (size > 0)
			    memcpy(parent_data[13],
				   dof->data_face + k * dof->dim, size);
		    }
		}
	    }
	}		/* if (level > 0) */
	if (dof->count_elem > 0) {
	    if (!DofIsHP(dof)) {
		/* non h-p DOF */
		size = dof->count_elem * sizeof(FLOAT);
		memcpy(parent_data[14], DofElementData(dof, e->index), size);
	    }
	    else {
		/* h-p DOF */
		i = e->index;
		size = sizeof(FLOAT) * dof->dim * (
					hp_old->elem_index[i + 1] -
					(k = hp_old->elem_index[i]));
		assert(size / sizeof(FLOAT) == dof->dim * (
			hp_new->elem_index[e_new->index + 1] -
			hp_new->elem_index[e_new->index]));
		if (size > 0) {
		    memcpy(parent_data[14],
			   dof->data_elem + k * dof->dim, size);
		}
	    }
	}
	return;
    }

    if (f != DofInterpolation) {
	/* DOF with an attached user function, compute values on the element */

	if (pass == 0)
	    return;	/* only compute in pass 1 (fewer computations) */

	assert(level == 0);

	if (dof->count_vert > 0) {
	    if (!(flags[i = e->verts[0]] & 1)) {
		flags[i] |= 1;
		type->InitFunc(dof, e_new, VERTEX, 0, f, dof->userfunc_lambda,
					NULL, parent_data[0], parent_data);
	    }
	    if (!(flags[i = e->verts[1]] & 1)) {
		flags[i] |= 1;
		type->InitFunc(dof, e_new, VERTEX, 1, f, dof->userfunc_lambda,
					NULL, parent_data[1], parent_data);
	    }
	    if (!(flags[i = e->verts[2]] & 1)) {
		flags[i] |= 1;
		type->InitFunc(dof, e_new, VERTEX, 2, f, dof->userfunc_lambda,
					NULL, parent_data[2], parent_data);
	    }
	    if (!(flags[i = e->verts[3]] & 1)) {
		flags[i] |= 1;
		type->InitFunc(dof, e_new, VERTEX, 3, f, dof->userfunc_lambda,
					NULL, parent_data[3], parent_data);
	    }
	}
	if (dof->count_edge > 0) {
	    if (!(flags[i = e->edges[0]] & 2)) {
		flags[i] |= 2;
		type->InitFunc(dof, e_new, EDGE, 0, f, dof->userfunc_lambda,
					NULL, parent_data[4 + 0], parent_data);
	    }
	    if (!(flags[i = e->edges[1]] & 2)) {
		flags[i] |= 2;
		type->InitFunc(dof, e_new, EDGE, 1, f, dof->userfunc_lambda,
					NULL, parent_data[4 + 1], parent_data);
	    }
	    if (!(flags[i = e->edges[2]] & 2)) {
		flags[i] |= 2;
		type->InitFunc(dof, e_new, EDGE, 2, f, dof->userfunc_lambda,
					NULL, parent_data[4 + 2], parent_data);
	    }
	    if (!(flags[i = e->edges[3]] & 2)) {
		flags[i] |= 2;
		type->InitFunc(dof, e_new, EDGE, 3, f, dof->userfunc_lambda,
					NULL, parent_data[4 + 3], parent_data);
	    }
	    if (!(flags[i = e->edges[4]] & 2)) {
		flags[i] |= 2;
		type->InitFunc(dof, e_new, EDGE, 4, f, dof->userfunc_lambda,
					NULL, parent_data[4 + 4], parent_data);
	    }
	    if (!(flags[i = e->edges[5]] & 2)) {
		flags[i] |= 2;
		type->InitFunc(dof, e_new, EDGE, 5, f, dof->userfunc_lambda,
					NULL, parent_data[4 + 5], parent_data);
	    }
	}
	if (dof->count_face > 0) {
	    if (!(flags[i = e->faces[0]] & 4)) {
		flags[i] |= 4;
		type->InitFunc(dof, e_new, FACE, 0, f, dof->userfunc_lambda,
					NULL, parent_data[10 + 0], parent_data);
	    }
	    if (!(flags[i = e->faces[1]] & 4)) {
		flags[i] |= 4;
		type->InitFunc(dof, e_new, FACE, 1, f, dof->userfunc_lambda,
					NULL, parent_data[10 + 1], parent_data);
	    }
	    if (!(flags[i = e->faces[2]] & 4)) {
		flags[i] |= 4;
		type->InitFunc(dof, e_new, FACE, 2, f, dof->userfunc_lambda,
					NULL, parent_data[10 + 2], parent_data);
	    }
	    if (!(flags[i = e->faces[3]] & 4)) {
		flags[i] |= 4;
		type->InitFunc(dof, e_new, FACE, 3, f, dof->userfunc_lambda,
					NULL, parent_data[10 + 3], parent_data);
	    }
	}
	if (dof->count_elem > 0) {
	    if (DofIsHP(dof))
		type = hp_new->info->types[hp_new->elem_order[e_new->index]];
	    type->InitFunc(dof, e_new, VOLUME, 0, f, dof->userfunc_lambda,
					NULL, parent_data[14], parent_data);
	}
	return;
    }

    /* recursively use InterpF2C (for geom and optionally other non h-p DOFs) */

    /* allocate buffer for storing data for e0 and e1 in the following order:
     *		1 new vertex (vn),
     *		2 cut edges (v0-vn and v1-vn),
     *		2 new edges (v2-vn and v3-vn),
     *		2 cut faces on f2 (v0-v3-vn and v1-v3-vn),
     *		2 cut faces on f3 (v0-v2-vn and v1-v2-vn),
     *		1 new face,
     *		2 new elements (e0 and e1) */

    assert(!DofIsHP(dof));	/* must be a non h-p DOF */

    size = 1 * dof->count_vert + 4 * dof->count_edge + 5 * dof->count_face +
	   2 * dof->count_elem;
    tmp = phgCalloc(size, sizeof(FLOAT));

    if (e0 != NULL) {
	/* build parent_data[] pointers (pointing to tmp) for e0 */
	phgMapP2C(e, vmap, emap, 0);
	p = tmp;
	if ((size = dof->count_vert) > 0) {
	    /* vertex data */
	    data[vmap[0]] = parent_data[0];
	    data[vmap[2]] = parent_data[2];
	    data[vmap[3]] = parent_data[3];
	    data[3] = p;
	    p += 1 * size;
	}
	if ((size = dof->count_edge) > 0) {
	    /* edge data */
	    data[4 + emap[1]] = parent_data[4 + 1];
	    data[4 + emap[2]] = parent_data[4 + 2];
	    data[4 + emap[5]] = parent_data[4 + 5];
	    data[4 + GetEdgeNo(3, vmap[0])] = p + 0 * size;
	    data[4 + GetEdgeNo(3, vmap[2])] = p + 2 * size;
	    data[4 + GetEdgeNo(3, vmap[3])] = p + 3 * size;
	    p += 4 * size;
	}
	if ((size = dof->count_face) > 0) {
	    /* face data */
	    data[10 + 3] = parent_data[10 + 1];
	    data[10 + vmap[2]] = p + 0 * size;
	    data[10 + vmap[3]] = p + 2 * size;
	    data[10 + vmap[0]] = p + 4 * size;
	    p += 5 * size;
	}
	if ((size = dof->count_elem) > 0) {
	    /* element data */
	    data[14] = p + 0 * size;
	}
	level++;
	dof_interp_recursive(e0, data);
	level--;
    }

    if (e1 != NULL) {
	/* build parent_data[] pointers (pointing to tmp) for e1 */
	phgMapP2C(e, vmap, emap, 1);
	p = tmp;
	if ((size = dof->count_vert) > 0) {
	    /* vertex data */
	    data[vmap[1]] = parent_data[1];
	    data[vmap[2]] = parent_data[2];
	    data[vmap[3]] = parent_data[3];
	    data[3] = p;
	    p += 1 * size;
	}
	if ((size = dof->count_edge) > 0) {
	    /* edge data */
	    data[4 + emap[3]] = parent_data[4 + 3];
	    data[4 + emap[4]] = parent_data[4 + 4];
	    data[4 + emap[5]] = parent_data[4 + 5];
	    data[4 + GetEdgeNo(3, vmap[1])] = p + 1 * size;
	    data[4 + GetEdgeNo(3, vmap[2])] = p + 2 * size;
	    data[4 + GetEdgeNo(3, vmap[3])] = p + 3 * size;
	    p += 4 * size;
	}
	if ((size = dof->count_face) > 0) {
	    /* face data */
	    data[10 + 3] = parent_data[10 + 0];	/* old face */
	    data[10 + vmap[2]] = p + 1 * size;
	    data[10 + vmap[3]] = p + 3 * size;
	    data[10 + vmap[1]] = p + 4 * size;
	    p += 5 * size;
	}
	if ((size = dof->count_elem) > 0) {
	    /* element data */
	    data[14] = p + 1 * size;
	}
	level++;
	dof_interp_recursive(e1, data);
	level--;
    }

    /* build children_data[] pointers */
    p = tmp;
    data[ 0] = p;   p += dof->count_vert;
    data[ 1] = p;   p += dof->count_edge;
    data[ 2] = p;   p += dof->count_edge;
    data[ 3] = p;   p += dof->count_edge;
    data[ 4] = p;   p += dof->count_edge;
    data[ 5] = p;   p += dof->count_face;
    data[ 6] = p;   p += dof->count_face;
    data[ 7] = p;   p += dof->count_face;
    data[ 8] = p;   p += dof->count_face;
    data[ 9] = p;   p += dof->count_face;
    data[10] = p;   p += dof->count_elem;
    data[11] = p;

    /* interpolate e->children ==> e on old mesh */
    dof->g = g_save;
    dof->hp = hp_old;
    type->InterpF2C(dof, e, parent_data, data);
    dof->hp = hp_new;
    dof->g = g_;

    phgFree(tmp);
}

static void
dof_interp0(void)
/* Note: e is an element in g_save */
{
    FLOAT *parent_data[15];
    FLOAT *p;

    e_new = e_old->parent;	/* points to the matching element in g */

    /* interpolate on e_new */
    p = dof_data_new;
    if (dof->count_vert > 0) {
	if (!DofIsHP(dof)) {
	    parent_data[0] = p + e_new->verts[0] * dof->count_vert;
	    parent_data[1] = p + e_new->verts[1] * dof->count_vert;
	    parent_data[2] = p + e_new->verts[2] * dof->count_vert;
	    parent_data[3] = p + e_new->verts[3] * dof->count_vert;
	}
	else {
	    parent_data[0] = p + hp_new->vert_index[e_new->verts[0]]*dof->dim;
	    parent_data[1] = p + hp_new->vert_index[e_new->verts[1]]*dof->dim;
	    parent_data[2] = p + hp_new->vert_index[e_new->verts[2]]*dof->dim;
	    parent_data[3] = p + hp_new->vert_index[e_new->verts[3]]*dof->dim;
	}
	p += DofVertexDataCount(dof);
    }
    if (dof->count_edge > 0) {
	if (!DofIsHP(dof)) {
	    parent_data[4] = p + e_new->edges[0] * dof->count_edge;
	    parent_data[5] = p + e_new->edges[1] * dof->count_edge;
	    parent_data[6] = p + e_new->edges[2] * dof->count_edge;
	    parent_data[7] = p + e_new->edges[3] * dof->count_edge;
	    parent_data[8] = p + e_new->edges[4] * dof->count_edge;
	    parent_data[9] = p + e_new->edges[5] * dof->count_edge;
	}
	else {
	    parent_data[4] = p + hp_new->edge_index[e_new->edges[0]]*dof->dim;
	    parent_data[5] = p + hp_new->edge_index[e_new->edges[1]]*dof->dim;
	    parent_data[6] = p + hp_new->edge_index[e_new->edges[2]]*dof->dim;
	    parent_data[7] = p + hp_new->edge_index[e_new->edges[3]]*dof->dim;
	    parent_data[8] = p + hp_new->edge_index[e_new->edges[4]]*dof->dim;
	    parent_data[9] = p + hp_new->edge_index[e_new->edges[5]]*dof->dim;
	}
	p += DofEdgeDataCount(dof);
    }
    if (dof->count_face > 0) {
	if (!DofIsHP(dof)) {
	    parent_data[10] = p + e_new->faces[0] * dof->count_face;
	    parent_data[11] = p + e_new->faces[1] * dof->count_face;
	    parent_data[12] = p + e_new->faces[2] * dof->count_face;
	    parent_data[13] = p + e_new->faces[3] * dof->count_face;
	}
	else {
	    parent_data[10] = p + hp_new->face_index[e_new->faces[0]]*dof->dim;
	    parent_data[11] = p + hp_new->face_index[e_new->faces[1]]*dof->dim;
	    parent_data[12] = p + hp_new->face_index[e_new->faces[2]]*dof->dim;
	    parent_data[13] = p + hp_new->face_index[e_new->faces[3]]*dof->dim;
	}
	p += DofFaceDataCount(dof);
    }
    if (dof->count_elem > 0) {
	if (!DofIsHP(dof)) {
	    parent_data[14] = p + e_new->index * dof->count_elem;
	}
	else {
	    parent_data[14] = p + hp_new->elem_index[e_new->index] * dof->dim;
	}
    }

    dof_interp_recursive(e_old, parent_data);
}

static void
dof_interp_func(DOF *u, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    int v[NVert], k;
    FLOAT buffers[3][Dim + 1], *lam, *lam0, *lam1, a0 = 0., a1 = 0.;
    ELEMENT *e0, *e1;

    assert(u == dof);
    assert(e == e_old->parent);

    e = e_old;	/* evaluate on old mesh */

    /* trace to the leaf element matching lambda[] */
    k = 0;	/* lam <=> k, lam0 <=> k + 1, lam1 <=> k + 2 */
    lam = (FLOAT *)lambda;
    lam0 = buffers[1];
    lam1 = buffers[2];
    while (TRUE) {
	e0 = e->children[0];
	e1 = e->children[1];

	if (e0 == NULL && e1 == NULL)
	    break;

	if (e0 != NULL) {
	    phgMapP2C(e, v, NULL, 0);
	    lam0[3]	= 2. * lam[1];
	    lam0[v[0]]	= lam[0] - lam[1];
	    lam0[v[2]]	= lam[2];
	    lam0[v[3]]	= lam[3];
	    a0 = lam0[0];
	    if (a0 > lam0[1])
		a0 = lam0[1];
	    if (a0 > lam0[2])
		a0 = lam0[2];
	    if (a0 > lam0[3])
		a0 = lam0[3];
	}

	if (e0 == NULL || (a0 < 0. && e1 != NULL)) {
	    phgMapP2C(e, v, NULL, 1);
	    lam1[3]	= 2. * lam[0];
	    lam1[v[1]]	= lam[1] - lam[0];
	    lam1[v[2]]	= lam[2];
	    lam1[v[3]]	= lam[3];
	    a1 = lam1[0];
	    if (a1 > lam1[1])
		a1 = lam1[1];
	    if (a1 > lam1[2])
		a1 = lam1[2];
	    if (a1 > lam1[3])
		a1 = lam1[3];
	}

	/* use the child on which the smallest lambda has largest value */
	if (e0 != NULL && (a0 >= 0 || e1 == NULL || a0 >= a1)) {
	    /* use child 0 */
	    lam = lam0;
	    e = e0;
	    k = (k + 1) % 3;
	    lam0 = buffers[(k + 1) % 3];
	    lam1 = buffers[(k + 2) % 3];
	}
	else {
	    /* use child 1 */
	    lam = lam1;
	    e = e1;
	    k = (k + 2) % 3;
	    lam0 = buffers[(k + 1) % 3];
	    lam1 = buffers[(k + 2) % 3];
	}
    }

    /* Evaluate on old mesh */
    dof->g = g_save;
    dof->hp = hp_old;
    phgDofEval(u, e, lam, values);
    /* need \grad u if D_order == 1 */
    k = (u->type != NULL ? u->type->D_order : u->hp->info->types[0]->D_order);
    assert(k <= 1);
    if (k == 1)
	phgDofEvalGradient(u, e, lam, NULL, values + DofDim(u));
    dof->hp = hp_new;
    dof->g = g_;
}

static void
dof_interp(void)
{
    INT i;
    DOF_USER_FUNC f = dof->userfunc;
    DOF_USER_FUNC_LAMBDA f_lambda = dof->userfunc_lambda;
    DOF_TYPE *type = (DofIsHP(dof) ? dof->hp->info->types[dof->hp->max_order] :
				     dof->type);

    hp_new = dof->hp;
    dof_data_new = phgCalloc(DofDataCount(dof), sizeof(FLOAT));

    i = 0;
    if (dof->count_face > 0 && i < g_save->nface)
	i = g_save->nface;
    if (dof->count_edge > 0 && i < g_save->nedge)
	i = g_save->nedge;
    if (dof->count_vert > 0 && i < g_save->nvert)
	i = g_save->nvert;
    flags = phgCalloc(i, sizeof(*flags));

    if (dof == g_->geom) {
	dof->data = g_save->geom->data;
	dof->data_vert = g_save->geom->data_vert;
	dof->data_edge = g_save->geom->data_edge;
	dof->data_face = g_save->geom->data_face;
	dof->data_elem = g_save->geom->data_elem;
    }

    /* Recursive interpolation can't be used for h-p DOFs, since the
     * orders on intermediate elements are unavailable. h-p DOFs are
     * computed with InitFunc using dof_interp_func() */

    if (dof->userfunc == DofInterpolation &&
#if 0 /* Note: set to 1 will force recursive interpolation for non hp DOFs */
	DofIsHP(dof) &&
#endif	/* 0|1 */
	type->BasFuncs != NULL && type->InitFunc != NULL) {
	/* change dof->userfunc_lambda to dof_interp_func() */
	dof->userfunc = NULL;
	dof->userfunc_lambda = dof_interp_func;
    }

    /* Note there are two cases:
     *	1. recursive interpolation (geom and optionally non hp DOFs), one pass.
     *	2. direct evaluation (when dof->userfunc != DofInterpolation, including
     *	   hp DOFs, DOFs with an attached user func, and optionally non h-p
     *	   DOFs), two passes:
     *		pass 0 copies data on uncoarsened elements
     *		pass 1 evaluates on coarsened elements
     */

    pass = 0;
    if (dof->userfunc != DofInterpolation) {
	/* first copy all valid data to avoid redundant computation */
	for (i = 0; i < g_save->nroot; i++) {
	    e_old = g_save->roots + i;
	    dof_interp0();
	}
	pass++;
    }

    for (i = 0; i < g_save->nroot; i++) {
	e_old = g_save->roots + i;
	dof_interp0();
    }

    /* restore userfunc */
    dof->userfunc = f;
    dof->userfunc_lambda = f_lambda;

    phgFree(flags);
    if (dof != g_->geom)
	phgFree(dof->data);

    dof->data = dof_data_new;
    dof->data_vert = dof->data;
    dof->data_edge = dof->data_vert + DofVertexDataCount(dof);
    dof->data_face = dof->data_edge + DofEdgeDataCount(dof);
    dof->data_elem = dof->data_face + DofFaceDataCount(dof);
}

/*----------------------- end DOF interpolation ------------------------*/

void
phgCoarsenMarkedElements(GRID *g)
/* coarsen all elements whose 'mark' member is negative */
{
    INT i;
    DOF **pdof;
    ELEMENT *roots, *e;
    double time;
    int halo_type = -1;

    FunctionEntry_head;
    ParallelFunction(phgCoarsenMarkedElements, g != NULL && g->nprocs > 1);
    FunctionEntry_tail;

    if (g == NULL || g->nprocs <= phgRank)
	Return;

    time = phgGetTime(NULL);

    if (g->nprocs > 1) {
	phgPrintf("\n\n"
"****************************************************************************\n"
"* Sorry, parallel coarsening is unsupported by this version of coarsen.o.  *\n"
"*                                                                          *\n"
"* You may consider replacing the file src/coarsen.o with one from the obj/ *\n"
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

    phgDebug((2, "before coarsening: nleaf=%"dFMT", nleaf_global=%"dFMT"\n", 
		g->nleaf, g->nleaf_global));

    g_ = g;

    /* set e->mark, reset e->flag, and set global indices */
    for (i = 0; i < g->nroot; i++)
	set_mark(g, g->roots + i);

    /* recursively set CHECK_FAILED mark on uncoarsenable patches */
    for (i = 0; i < g->nroot; i++)
	check_patch(g, g->roots + i);


    /* set dof_interp_flag, construct list of h-p DOFs */
    dof_interp_flag = FALSE;
    nhp = 0;
    hp_orders = NULL;
    if (g->dof != NULL || g->hp != NULL) {
	HP_TYPE **php;
	for (php = g->hp, i = 0; php != NULL && *php != NULL; php++)
	    (*php)->id = i++;
	nhp = i;
	if (nhp > 0) {
	    dof_interp_flag = TRUE;	/* interpolation needed */
	}
	else {
	    for (pdof = g->dof; (dof = *pdof) != NULL; pdof++) {
		if (SpecialDofType(dof->type) || (!DofIsHP(dof) &&
		    (dof->data == NULL || dof->userfunc == DofNoAction)))
		    continue;
		dof_interp_flag = TRUE;	/* interpolation needed */
		break;
	    }
	}
	if (dof_interp_flag) {
	    hp_orders = phgAlloc(nhp * g->nelem * sizeof(*hp_orders));
	    memset(hp_orders, -1, nhp * g->nelem * sizeof(*hp_orders));
	    save_hp_orders(g);
	}
    }


    /* now hp_orders[] contains correct element h-p orders on the coarsened
     * mesh, and nleaf_new should equal to nleaf of the new mesh */
    if (dof_interp_flag) {
	g_save = phgNewGrid(FACE_FLAG | ELEM_FLAG);	/* needed by geom.c */
	g_save->nroot = nleaf_new;
	g_save->roots = phgAlloc(nleaf_new * sizeof(*g_save->roots));

	g_save->nelem = g->nelem;
	g_save->nface = g->nface;
	g_save->nedge = g->nedge;
	g_save->nvert = g->nvert;

	phgGeomInit_(g_save, FALSE);
	g_save->geom->data = g->geom->data;
	g_save->geom->data_vert = g_save->geom->data;
	g_save->geom->data_edge = g_save->geom->data_vert
					+ DofVertexDataCount(g_save->geom);
	g_save->geom->data_face = g_save->geom->data_edge
					+ DofEdgeDataCount(g_save->geom);
	g_save->geom->data_elem = g_save->geom->data_face
					+ DofFaceDataCount(g_save->geom);
	g->geom->data = NULL;
	g->geom->data_vert = NULL;
	g->geom->data_edge = NULL;
	g->geom->data_face = NULL;
	g->geom->data_elem = NULL;
    }

    g->ntree = g->nleaf = 0;
    for (i = 0; i < g->nroot; i++) {
	do_coarsen(g, g->roots + i, &g->ntree, &g->nleaf);
    }
    assert(g_save == NULL || g->nleaf == g_save->nroot);

    if (g->nprocs == 1) {
	g->nelem_global = g->nleaf;
    }

    if (g->period != NULL) {
	phgFree(g->period->ordering);
	g->period->ordering = NULL;
    }

    pack_indices(g, VERTEX);
    if ((g->flags & EDGE_FLAG))
	pack_indices(g, EDGE);
    if ((g->flags & FACE_FLAG))
	pack_indices(g, FACE);
    if ((g->flags & ELEM_FLAG))
	pack_indices(g, VOLUME);

    update_neighbours(g);

    {
	if (g->period != NULL)
	    pack_period_L2Gmap(g);
	phgUpdateBoundaryTypes(g);
    }

    /* clear QUAD and DOF cache */
    phgDofClearCache(g, NULL, NULL, NULL, FALSE);

    if (dof_interp_flag) {
	HP_TYPE **oldhps;
	int ihp;

	/* recalculate g->geom (the old geom data are saved in g_save->geom),
	 * FIXME: only recompute missing data */
	phgGeomInit(g);

	/* update HP_TYPES, save old HP_TYPEs in oldhps[] */
	oldhps = phgAlloc(nhp * sizeof(*oldhps));
	for (ihp = 0; ihp < nhp; ihp++) {
	    oldhps[ihp] = phgHPDup__(g->hp[ihp], TRUE);
	    g->hp[ihp]->elem_order = phgCalloc(g->nelem,
					     sizeof(*oldhps[ihp]->elem_order));
	    roots = g_save->roots;
	    for (i = 0; i < g_save->nroot; i++, roots++) {
		e = roots->parent;
		g->hp[ihp]->elem_order[e->index] =
				hp_orders[ihp * g_save->nelem + roots->index];
	    }
	    /* FIXME: combine the communications for all h-p DOFs */
	    phgHPSetupOrders(g->hp[ihp]);
	}

	phgFree(hp_orders);

	/* interpolate DOFs with g->geom containing new goem data */
	for (pdof = g->dof; (dof = *pdof) != NULL; pdof++) {
	    if (dof == g->geom ||
		SpecialDofType(dof->type) ||
		dof->data == NULL)
		continue;

	    if (dof->userfunc == DofNoAction) {
		/* allocate new data buffer */
		phgFree(dof->data);
		dof->data = phgCalloc(DofDataCount(dof), sizeof(*dof->data));
		dof->data_vert = dof->data;
		dof->data_edge = dof->data_vert + DofVertexDataCount(dof);
		dof->data_face = dof->data_edge + DofEdgeDataCount(dof);
		dof->data_elem = dof->data_face + DofFaceDataCount(dof);
		continue;
	    }

	    if (DofIsHP(dof))
		hp_old = oldhps[dof->hp->id];
	    else
		hp_old = NULL;

	    dof_interp();
	}

	for (i = 0; i < nhp; i++) {
	    phgHPFreeBuffers__(oldhps[i]);
	    phgFree(oldhps[i]);
	}
	phgFree(oldhps);
    }

    phgFreeGrid(&g_save);

#if 0
    if (halo_type >= 0)
	phgSetupHalo(g, halo_type);
#endif

    if (g->rank == 0 && phgVerbosity > 0)
	phgPrintf("%s: nleaf_global = %"dFMT", time = %0.4lf\n", __func__,
			g->nleaf_global, phgGetTime(NULL) - time);

    phgDebug((2, "after coarsening nleaf=%"dFMT", nleaf_global=%"dFMT"\n", 
		g->nleaf, g->nleaf_global));

#if ALLOW_CURVED_BOUNDARY
    if (g->bdry_funcs != NULL) {
	/* update volume */
	FLOAT d = 0.;
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

#if USE_MPI
    if (g->nprocs > 1) {
	if (phgSameProcessorSpeed()) {
	    MPI_Allreduce(&g->nleaf, &i, 1, PHG_MPI_INT, MPI_MAX, g->comm);
	    g->lif = i * (FLOAT)g->nprocs / g->nelem_global;
	}
	else {
	    /* Update g->lif */
	    double w0[3], w[3];
	    w0[0] = w0[1] = w0[2] = g->nleaf / phgGetProcessorSpeed(NULL);
	    MPI_Allreduce(w0, w, 1, MPI_3DOUBLE, MPI_MSM, g->comm);
	    g->lif = w[2] * g->nprocs / w[1];
	}
    }
#endif	/* USE_MPI */

    PrintTime(2);

    Return;
}
