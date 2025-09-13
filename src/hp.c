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

/* $Id: hp.c,v 1.113 2022/04/09 01:28:34 zlb Exp $
 *
 * Support functions for hp-adaptivity */

#include "phg.h"
#include <string.h>

HP_TYPE *
phgHPNew_(GRID *g, HP_INFO *info, int order, const char *file, int line)
/* allocates and returns a new HP_TYPE object.
 * 'order' is the initial order to use (-1 ==> no initial order) */
{
    HP_TYPE *hp;
    DOF_TYPE *type = info->types[info->max_order];
    size_t len;
    INT i;


    hp = phgCalloc(1, sizeof(*hp));
    hp->info = info;
    hp->srcfile = strdup(file);
    hp->srcline = line;
    hp->g = g;

    if (order >= 0) {
	if (order < info->min_order)
	    order = info->min_order;
	else if (order > info->max_order)
	    order = info->max_order;
    }

#define SetOrder(np_foo, nfoo, foo_order, flag) \
    if (type->np_foo > 0) { \
	    hp->foo_order = phgCalloc(g->nfoo, sizeof(*hp->foo_order)); \
	    for (i = 0; order >= 0 && i < g->nfoo; i++) \
		hp->foo_order[i] = order; \
    }
    SetOrder(np_vert, nvert, vert_order, VERT_FLAG)
    SetOrder(np_edge, nedge, edge_order, EDGE_FLAG)
    SetOrder(np_face, nface, face_order, FACE_FLAG)
    SetOrder(np_elem, nelem, elem_order, ELEM_FLAG)
#undef SetOrder

    if (order >= 0)
	phgHPSetupDataPointers(hp);

    /* attach hp to g->hp */
    len = 0;
    if (g->hp != NULL) {
	HP_TYPE **p;
	for (p = g->hp; *p != NULL; p++, len++);
    }
    g->hp = phgRealloc_(g->hp, (len + 2) * sizeof(*(g->hp)),
				len * sizeof(*(g->hp)));
    g->hp[len] = hp;
    g->hp[len + 1] = NULL;

    return hp;
}

void
phgHPSetupDataPointers_(HP_TYPE *hp, INT nvert, INT nedge, INT nface, INT nelem,
			BOOLEAN use_types, BOOLEAN update_counters)
 /* initializes data pointers, i.e., sets up {edge,face,elem}_index arrays,
  * if 'update_counters' == TRUE than also updates global data counters.
  * if 'use_types' == FALSE than don't use g->types_xxxx[] */
{
    GRID *g = hp->g;
    INT i, k, a[4];

    assert(!update_counters || use_types);

#define SetIndex(nfoo, np_foo, foo_types, foo_index, foo_order, a, flag) \
    if (hp->foo_index != NULL) { \
	phgFree(hp->foo_index); \
	hp->foo_index = NULL; \
    } \
    a = 0; \
    if ((hp->info->flags & flag)) { \
	hp->foo_index = phgAlloc((nfoo + 1) * sizeof(*hp->foo_index)); \
	hp->foo_index[0] = 0; \
	for (i = 0; i < nfoo; i++) { \
	    if (use_types && g->foo_types[i] == UNREFERENCED) { \
		hp->foo_index[i + 1] = hp->foo_index[i]; \
		continue; \
	    } \
	    k = hp->info->types[hp->foo_order[i]]->np_foo; \
	    hp->foo_index[i + 1] = hp->foo_index[i] + k; \
	    if (update_counters && (!use_types || g->foo_types[i] & OWNER)) \
		a += k; \
	} \
    }
    SetIndex(nvert, np_vert, types_vert, vert_index, vert_order, a[0],
	     VERT_FLAG)
    SetIndex(nedge, np_edge, types_edge, edge_index, edge_order, a[1],
	     EDGE_FLAG)
    SetIndex(nface, np_face, types_face, face_index, face_order, a[2],
	     FACE_FLAG)
    SetIndex(nelem, np_elem, types_elem, elem_index, elem_order, a[3],
	     ELEM_FLAG)
#undef SetIndex

    if (update_counters) {
#if USE_MPI
	INT b[4];
	MPI_Allreduce(a, b, 4, PHG_MPI_INT, MPI_SUM, g->comm);
	hp->vert_data_count_global = b[0];
	hp->edge_data_count_global = b[1];
	hp->face_data_count_global = b[2];
	hp->elem_data_count_global = b[3];

	a[0] = -hp->min_order;
	a[1] = hp->max_order;
	MPI_Allreduce(a, b, 2, PHG_MPI_INT, MPI_MAX, g->comm);
	hp->min_order = -b[0];
	hp->max_order = b[1];
#else	/* USE_MPI */
	hp->vert_data_count_global = a[0];
	hp->edge_data_count_global = a[1];
	hp->face_data_count_global = a[2];
	hp->elem_data_count_global = a[3];
#endif	/* USE_MPI */
    }

    return;
}

#if 0
# define CMP_ORDER(a,b)	((a) > (b))	/* order = max(orders) (don't use!) */
#else
# define CMP_ORDER(a,b)	((a) < (b))	/* order = min(orders) */
#endif

void
phgHPSetupOrders(HP_TYPE *hp)
/* assumes hp->elem_order[] are valid and updates other members */
{
    GRID *g = hp->g;
    INT i;
    ELEMENT *e;
    int ii, order;

    assert((hp->info->flags & ELEM_FLAG));

#define SetOrder(nfoo, NFoo, foo_order, foos, flag) \
    if ((hp->info->flags & flag)) { \
	phgFree(hp->foo_order); \
	hp->foo_order = phgCalloc(g->nfoo, sizeof(*hp->foo_order)); \
	for (i = 0; i < g->nfoo; i++) \
	    hp->foo_order[i] = (CMP_ORDER(255, 0) ? 0 : 255); \
	ForAllElements_(g, e, FOR_ALL) { \
	    order = hp->elem_order[e->index]; \
	    for (ii = 0; ii < NFoo; ii++) { \
		i = e->foos[ii]; \
		if (!CMP_ORDER(hp->foo_order[i], order)) \
		    hp->foo_order[i] = order; \
	    } \
	} \
    }
    SetOrder(nvert, NVert, vert_order, verts, VERT_FLAG)
    SetOrder(nedge, NEdge, edge_order, edges, EDGE_FLAG)
    SetOrder(nface, NFace, face_order, faces, FACE_FLAG)
#undef SetOrder

#if USE_MPI
    if (g->nprocs > 1 && (g->halo != NULL ||
	(hp->info->flags & (VERT_FLAG | EDGE_FLAG | FACE_FLAG)))) {
	/* Note: a dummy MAP is used to build communication pattern */
	static DOF_TYPE type = {DofReserved,
	    "dummy", NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	    NULL, NULL, NULL, NULL, FE_None,
	    FALSE, TRUE, FALSE, -1,
	    NVert + NEdge + NFace + 1, 0, 0, -1, -1, 1, 1, 1, 1
	};
	MAP *map;
	BYTE *sbuff, *rbuff;
	INT *V2Dmap;
	DOF *utmp;
	INT k;

	utmp = phgDofNew(g, &type, 1, "dummy", DofNoData);
	map = phgMapCreate(utmp, NULL);

	/* build V2Dmap (reverse of D2Lmap * L2Vmap) */
	V2Dmap = phgAlloc(map->localsize * sizeof(*V2Dmap));
	for (i = 0; i < DofDataCount(utmp); i++) {
	    if ((k = phgMapD2L(map, 0, i)) < 0)
		continue;
	    V2Dmap[map->L2Vmap == NULL ? k : map->L2Vmap[k]] = i;
	}

	map->cinfo = phgMapCreateCommInfo(
				map->comm,
				map->nprocs,
				map->rank,
				map->nlocal,
				map->localsize,
				map->O2Gmap,
				map->ordering,
				map->partition);
	sbuff = phgCalloc((map->cinfo->ssize + map->cinfo->rsize),
			  sizeof(*sbuff));
	rbuff = sbuff + map->cinfo->ssize;

	/* copy orders on unowned verts/edges/faces/elems to rbuff */
	for (i = 0; i < map->cinfo->rsize; i++) {
	    k = map->nlocal + map->cinfo->rind[i];
	    assert(k < map->localsize);
	    k = V2Dmap[k];
	    if (k < g->nvert) {
		/* vert data */
		if (hp->info->flags & VERT_FLAG)
		    rbuff[i] = hp->vert_order[k];
	    }
	    else if ((k -= g->nvert) < g->nedge) {
		/* edge data */
		if (hp->info->flags & EDGE_FLAG)
		    rbuff[i] = hp->edge_order[k];
	    }
	    else if ((k -= g->nedge) < g->nface) {
		/* face data */
		if (hp->info->flags & FACE_FLAG)
		    rbuff[i] = hp->face_order[k];
	    }
	    else {
		/* elem data */
		k -= g->nface;
		assert(k < g->nelem);
		rbuff[i] = hp->elem_order[k];
	    }
	}

	/* send orders on unowned verts/edges/faces to their owners */
	if (map->cinfo->comm == MPI_COMM_NULL)
	    goto skip_comm;
	MPI_Alltoallv(rbuff, map->cinfo->rcnts, map->cinfo->rdsps, PHG_MPI_BYTE,
		      sbuff, map->cinfo->scnts, map->cinfo->sdsps, PHG_MPI_BYTE,
		      map->cinfo->comm);

	/* compare with incoming orders */
	for (i = 0; i < map->cinfo->ssize; i++) {
	    k = map->cinfo->sind[i];
	    assert(k < map->nlocal);
	    k = V2Dmap[k];
	    if (k < g->nvert) {
		/* vert data */
		if ((hp->info->flags & VERT_FLAG) &&
		    !CMP_ORDER(hp->vert_order[k], sbuff[i]))
		    hp->vert_order[k] = sbuff[i];
	    }
	    else if ((k -= g->nvert) < g->nedge) {
		/* edge data */
		if ((hp->info->flags & EDGE_FLAG) &&
		    !CMP_ORDER(hp->edge_order[k], sbuff[i]))
		    hp->edge_order[k] = sbuff[i];
	    }
	    else if ((k -= g->nedge) < g->nface) {
		/* face data */
		if ((hp->info->flags & FACE_FLAG) &&
		    !CMP_ORDER(hp->face_order[k], sbuff[i]))
		    hp->face_order[k] = sbuff[i];
	    }
	    else {
		/* elem data */
		k -= g->nface;
		assert(k < g->nelem);
		if (!CMP_ORDER(hp->elem_order[k], sbuff[i]))
		    hp->elem_order[k] = sbuff[i];
	    }
	}

	/* update incoming orders */
	for (i = 0; i < map->cinfo->ssize; i++) {
	    k = map->cinfo->sind[i];
	    assert(k < map->nlocal);
	    k = V2Dmap[k];
	    if (k < g->nvert) {
		/* vert data */
		if (hp->info->flags & VERT_FLAG)
		    sbuff[i] = hp->vert_order[k];
	    }
	    else if ((k -= g->nvert) < g->nedge) {
		/* edge data */
		if (hp->info->flags & EDGE_FLAG)
		    sbuff[i] = hp->edge_order[k];
	    }
	    else if ((k -= g->nedge) < g->nface) {
		/* face data */
		if (hp->info->flags & FACE_FLAG)
		    sbuff[i] = hp->face_order[k];
	    }
	    else {
		/* elem data */
		k -= g->nface;
		assert(k < g->nelem);
		sbuff[i] = hp->elem_order[k];
	    }
	}

	/* send back incoming orders */
	MPI_Alltoallv(sbuff, map->cinfo->scnts, map->cinfo->sdsps, PHG_MPI_BYTE,
		      rbuff, map->cinfo->rcnts, map->cinfo->rdsps, PHG_MPI_BYTE,
		      map->cinfo->comm);

	/* update orders on unowned vert/edges/faces */
	for (i = 0; i < map->cinfo->rsize; i++) {
	    k = map->nlocal + map->cinfo->rind[i];
	    assert(k < map->localsize);
	    k = V2Dmap[k];
	    if (k < g->nvert) {
		/* vert data */
		if (hp->info->flags & VERT_FLAG)
		    hp->vert_order[k] = rbuff[i];
	    }
	    else if ((k -= g->nvert) < g->nedge) {
		/* edge data */
		if (hp->info->flags & EDGE_FLAG)
		    hp->edge_order[k] = rbuff[i];
	    }
	    else if ((k -= g->nedge) < g->nface) {
		/* face data */
		if (hp->info->flags & FACE_FLAG)
		    hp->face_order[k] = rbuff[i];
	    }
	    else {
		/* elem data */
		k -= g->nface;
		assert(k < g->nelem);
		hp->elem_order[k] = rbuff[i];
	    }
	}

skip_comm:
	phgFree(V2Dmap);
	phgFree(sbuff);
	phgMapDestroy(&map);
	phgDofFree(&utmp);
    }
#endif	/* USE_MPI */

    /* Note: global min/max order will be computed by phgHPSetupDataPointers */
    hp->min_order = hp->info->max_order;
    hp->max_order = hp->info->min_order;
    ForAllElements_(g, e, FOR_ALL) {
	order = hp->elem_order[e->index];
	assert(order >= hp->info->min_order && order <= hp->info->max_order);
	if (hp->min_order > order)
	    hp->min_order = order;
	if (hp->max_order < order)
	    hp->max_order = order;
    }

    phgHPSetupDataPointers(hp);
}

static DOF *u0;

static void
tmp_func(DOF *dof, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    Unused(bno);
    CheckD_order(dof);
    phgDofEval(u0, e, lambda, values);
}

static void
update_dof(DOF *u, HP_TYPE *hp, BOOLEAN interpolate)
/* update 'u' according to new hp settings */
{
    GRID *g = u->g;
    DOF *utmp;

    if (u->hp == hp)
	return;

    utmp = phgHPDofNew(g, hp, u->dim, "tmp", DofNoAction);

    if (interpolate) {
	if (u->userfunc == DofLambdaFunction) {
	    assert(u->userfunc_lambda != NULL);
	    phgDofSetDataByLambdaFunction(utmp, u->userfunc_lambda);
	}
	else if (u->userfunc != DofNoAction &&
		 u->userfunc != DofInterpolation) {
	    phgDofSetDataByFunction(utmp, u->userfunc);
	}
	else {
	    u0 = u;
	    phgDofSetDataByLambdaFunction(utmp, tmp_func);
	}
    }

    phgFree(u->data);
    u->data = u->data_vert = utmp->data_vert;
    u->data_edge = utmp->data_edge;
    u->data_face = utmp->data_face;
    u->data_elem = utmp->data_elem;

    utmp->data = NULL;
    phgDofFree(&utmp);

    phgDofClearCache(NULL, u, NULL, NULL, FALSE);
}

void
phgHPAttachDof(HP_TYPE *hp, DOF *u, BOOLEAN interpolate)
{
    size_t len;

    if (u->hp == hp)
	return;

    if ((!DofIsHP(u) && hp->info != u->type->hp_info) ||
	( DofIsHP(u) && hp->info != u->hp->info))
	phgError(1, "%s: can't use \"%s\" for DOF \"%s\".\n",
		    hp->info->name, u->name);

    update_dof(u, hp, interpolate);
    phgHPFree(&u->hp);
    u->hp = hp;

    /* attach u to hp->dof */
    len = 0;
    if (hp->dof != NULL) {
	DOF **p;
	for (p = hp->dof; *p != NULL; p++, len++);
    }
    hp->dof = phgRealloc_(hp->dof, (len + 2) * sizeof(*(hp->dof)),
				    len * sizeof(*(hp->dof)));
    hp->dof[len] = u;
    hp->dof[len + 1] = NULL;

    hp->refcount++;

    if (!DofIsHP(u))
	phgDofUnrefId(u->type);

    /* reset count_{vert,edge,face,elem} using max order */
    u->type = hp->info->types[hp->max_order];
    u->count_vert = u->type->np_vert * u->dim;
    u->count_edge = u->type->np_edge * u->dim;
    u->count_face = u->type->np_face * u->dim;
    u->count_elem = u->type->np_elem * u->dim;
    u->type = NULL;	/* accessing u->type will produce segfault */

    phgDofClearCache(NULL, u, NULL, NULL, FALSE);
}

void
phgHPSetup(HP_TYPE *hp, BOOLEAN interpolate)
/* Note: requested order on element e is stored in e->hp_order */
{
    HP_TYPE *hp_old;
    GRID *g = hp->g;
    int i;
    ELEMENT *e;

    if (hp == NULL)
	return;

    if (hp->derived)
	phgError(1, "shouldn't change orders of a derived HP_TYPE "
		    "(created at %s:%d)\n", hp->srcfile, hp->srcline);

    phgDebug((2, "setting up hierarchical orders.\n"));
    phgDebug((2, "g->nelem = %d.\n", g->nelem));

    hp_old = phgHPDup__(hp, TRUE);

    /* clean up order of unowned elements */
    ForAllElements_(g, e, FOR_UNOWNED)
	e->hp_order = hp->info->max_order;

    /* set element orders */
    phgFree(hp->elem_order);
    hp->elem_order = phgCalloc(g->nelem, sizeof(*hp->elem_order));
    ForAllElements_(g, e, FOR_ALL) {
	if (e->hp_order < hp->info->min_order)
	    hp->elem_order[e->index] = hp->info->min_order;
	else if (e->hp_order > hp->info->max_order)
	    hp->elem_order[e->index] = hp->info->max_order;
	else
	    hp->elem_order[e->index] = e->hp_order;
    }

    /* compute vertex, edge and face order */
    phgHPSetupOrders(hp);

    /* update DOFs using the HP_TYPE */
    if (hp->dof != NULL) {
	/* note: hp->dof might be reallocated in update_dof() by phgDofNew_ */
	for (i = 0; hp->dof[i] != NULL; i++) {
	    DOF *u = hp->dof[i];
	    u->hp = hp_old;
	    update_dof(u, hp, interpolate);
	    u->hp = hp;
	}
    }

    phgHPFreeBuffers__(hp_old);
    phgFree(hp_old);

    /* This will detach hp->derivative from hp */
    if (hp->derivative != NULL)
	phgHPFree(&hp->derivative);

    /* This will detach hp->dg from hp */
    if (hp->dg != NULL)
	phgHPFree(&hp->dg);

    phgDofClearCache(NULL, NULL, NULL, hp, FALSE);
}

HP_TYPE *
phgHPGetDerivative_(HP_TYPE *old_hp, const char *srcfile, int srcline)
/* returns a new HP_TYPE for the derivatives of 'old_hp' */
{
    GRID *g;
    HP_TYPE *hp;
    INT i;
    int order, min_order;

    if ((hp = old_hp->derivative) != NULL) {
	hp->refcount++;
	return hp;
    }

    hp = phgHPNew_(g = old_hp->g, HP_DG, -1, srcfile, srcline);
    old_hp->derivative = hp;
    hp->refcount++;	/* referenced by old_hp->derivative */

    min_order = hp->info->min_order;

    for (i = 0; i < g->nelem; i++) {
	if (g->types_elem[i] == UNREFERENCED)
	    continue;
	order = ((int)old_hp->info->types[old_hp->elem_order[i]]->order) - 1;
	if (order < min_order)
	    order = min_order;
	hp->elem_order[i] = order;
    }

    order = ((int)old_hp->info->types[old_hp->min_order]->order) - 1;
    if (order < min_order)
	order = min_order;
    hp->min_order = order;

    order = ((int)old_hp->info->types[old_hp->max_order]->order) - 1;
    if (order < min_order)
	order = min_order;
    hp->max_order = order;

    phgHPSetupDataPointers(hp);

    return hp;
}

HP_TYPE *
phgHPGetDG_(HP_TYPE *old_hp, const char *srcfile, int srcline)
/* returns a new DG type HP_TYPE of the same order as 'old_hp' */
{
    GRID *g;
    HP_TYPE *hp;
    INT i;
    int order, min_order, max_order;

    if ((hp = old_hp->dg) != NULL) {
	hp->refcount++;
	return hp;
    }

    hp = phgHPNew_(g = old_hp->g, HP_DG, -1, srcfile, srcline);
    old_hp->dg = hp;
    hp->refcount++;	/* referenced by old_hp->dg */

    min_order = HP_DG->max_order;
    max_order = HP_DG->min_order;
    for (i = 0; i < g->nelem; i++) {
	if (g->types_elem[i] == UNREFERENCED)
	    continue;
	order = old_hp->info->types[old_hp->elem_order[i]]->order;
	if (order < HP_DG->min_order)
	    order = HP_DG->min_order;
	else if (order > HP_DG->max_order)
	    order = HP_DG->max_order;

	hp->elem_order[i] = order;

	if (min_order > order)
	    min_order = order;
	if (max_order < order)
	    max_order = order;
    }

    hp->min_order = min_order;
    hp->max_order = max_order;

    phgHPSetupDataPointers(hp);

    return hp;
}

HP_TYPE *
phgHPDup__(HP_TYPE *hp_old, BOOLEAN free_orig_data)
/* Note: this function creates an 'invisible' HP_TYPE which is not registered
 * to GRID, should only be used as scratch variable, and be freed using
 * phgHPFreeBuffers__().
 *
 * If 'free_orig_data' then the xxxx_order and xxxx_index pointers in the
 * original HP_TYPE will be set to NULL. */
{
    GRID *g = hp_old->g;
    HP_TYPE *hp;

    hp = phgAlloc(sizeof(*hp));
    memcpy(hp, hp_old, sizeof(*hp));
    hp->dof = NULL;
    hp->refcount = 0;
    hp->srcfile = NULL;
    hp->derivative = NULL;
    hp->dg = NULL;

    if (free_orig_data) {
	hp_old->vert_order = NULL;
	hp_old->edge_order = NULL;
	hp_old->face_order = NULL;
	hp_old->elem_order = NULL;
	hp_old->vert_index = NULL;
	hp_old->edge_index = NULL;
	hp_old->face_index = NULL;
	hp_old->elem_index = NULL;
    }
    else {
#define CopyData(nfoo, foo_order, foo_index, flag) \
	if (hp->info->flags & flag) { \
	    hp->foo_order = phgAlloc(g->nfoo * sizeof(*hp->foo_order)); \
	    if (hp_old->foo_order != NULL) \
		memcpy(hp->foo_order, hp_old->foo_order, \
			g->nfoo * sizeof(*hp->foo_order)); \
	    hp->foo_index = phgAlloc((g->nfoo + 1) * sizeof(*hp->foo_index)); \
	    assert(hp_old->foo_index != NULL); \
	    memcpy(hp->foo_index, hp_old->foo_index, \
			(g->nfoo + 1) * sizeof(*hp->foo_index)); \
	}
	CopyData(nvert, vert_order, vert_index, VERT_FLAG)
	CopyData(nedge, edge_order, edge_index, EDGE_FLAG)
	CopyData(nface, face_order, face_index, FACE_FLAG)
	CopyData(nelem, elem_order, elem_index, ELEM_FLAG)
#undef CopyData
    }

    return hp;
}

void
phgHPFreeBuffers__(HP_TYPE *hp)
/* frees all memory blocks allocated by 'hp' */
{
    phgFree(hp->srcfile);
    hp->srcfile = NULL;

    phgFree(hp->vert_order);
    hp->vert_order = NULL;
    phgFree(hp->edge_order);
    hp->edge_order = NULL;
    phgFree(hp->face_order);
    hp->face_order = NULL;
    phgFree(hp->elem_order);
    hp->elem_order = NULL;

    phgFree(hp->vert_index);
    hp->vert_index = NULL;
    phgFree(hp->edge_index);
    hp->edge_index = NULL;
    phgFree(hp->face_index);
    hp->face_index = NULL;
    phgFree(hp->elem_index);
    hp->elem_index = NULL;
}

void
phgHPFree(HP_TYPE **hp_ptr)
{
    HP_TYPE *hp;
    HP_TYPE **p, **p0 = NULL;

    if (hp_ptr == NULL || (hp = *hp_ptr) == NULL)
	return;

    if (hp->refcount > 0) {
	--hp->refcount;
	return;
    }

    *hp_ptr = NULL;

    /* detach hp from g->hp */
    for (p = hp->g->hp; p != NULL && *p != NULL; p++) {
	if (*p == hp)
	    p0 = p;
    }
    if (p == NULL || p0 == NULL) {
	phgError(1, "%s:%d: unexpected error(grid freed before HP_TYPE?)\n",
			__FILE__, __LINE__);
    }
    assert(p > p0);
    memmove(p0, p0 + 1, sizeof(*p0) * (p - p0));

    /* free hp->dof */
    if (hp->dof != NULL) {
	if (*hp->dof != NULL)
	    phgWarning("%s: possibly unfreed DOFs found.\n", __func__);
	phgFree(hp->dof);
    }

    phgDofClearCache(NULL, NULL, NULL, hp, TRUE);

    phgHPFree(&hp->derivative);
    phgHPFree(&hp->dg);
    phgHPFreeBuffers__(hp);
    phgFree(hp);
}
