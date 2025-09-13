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

/* Functions managing maps
 *
 * $Id: map.c,v 1.230 2022/06/24 02:15:56 zlb Exp $ */

#define NEED_GetVariableArgs

#include "phg.h"

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdarg.h>
#include <math.h>
#include <limits.h>	/* INT_MAX */

enum {SFC_AUTO = 0, SFC_YES = 1, SFC_NO = 2};
static const char *sfc_flags[] = {
    "auto",	/* reorder unknowns using SFC if phgMaxThreads > 1 */
    "yes",	/* always reorder unknowns using SFC */
    "no",	/* don't reorder unknowns */
    NULL};
static int sfc_reorder = SFC_AUTO;

enum {SFC_HILBERT = 0, SFC_MORTON = 1, SFC_NONE = 2};
static const char *sfc_types[] = {"hilbert", "morton", "none", NULL};
static int sfc_type = SFC_HILBERT;

/* The new (much faster) algorithm using g->owner_xxxx_xxxx arrays.
 *
 * Note:
 *
 *   1. The unreferenced entries are not handled and map->D2Lmap is used to
 *	remove them from the map (with localsize adjusted accordingly)
 *
 *   2. When there are periodic boundaries, map->D2Lmap is also used to
 *	remove duplicate entries. */

static void
build_periodic_D2Lmap(MAP *map, INT *nonowned_count_vert)
/* builds map->D2Lmap for periodic boundaries, 
 *
 * Each group of periodic vertices is designated a local owner (which
 * matches the global ownership of the vertex). If a vertex is not the
 * local owner of its group then its D2Lmap entry is negated (-1 -.).
 *
 * This function also sets up map->D2Lmap_size, and computes nonowned_count_vert
 * on the fly. */
{
    GRID *g = map->dofs[0]->g;
    DOF *u;
    DOF_TYPE *type;
    INT i, j, k, l, m, dindex, lindex;
    int ii, jj, count;

    *nonowned_count_vert = g->nvert - g->nvert_owned;

    if (g->period == NULL) {
	map->D2Lmap_size = map->localsize;
	return;
    }

#if DEBUG
    if (g->nvert > 0) {
	assert(g->owner_index_vert != NULL && g->owner_rank_vert != NULL);
	assert(g->period != NULL && g->period->ordering != NULL);
    }
#endif	/* DEBUG */

    map->D2Lmap = phgAlloc(map->localsize * sizeof(*map->D2Lmap));

    dindex = 0;
    lindex = 0;
    for (ii = 0; ii < map->ndof; ii++) {
	u = map->dofs[ii];
	if ((type = u->type) == NULL)
	    type = u->hp->info->types[u->hp->max_order];

	/* vertex DOF */
	if (type->np_vert > 0) {
	    i = 0;
	    *nonowned_count_vert = 0;
	    while (i < g->nvert) {
		m = g->period->L2Gmap_vert[g->period->ordering[i]];
		l = i;
		/* choose between identical vertices the one having the
		 * OWNER flag, or the one with the smallest index */
		for (k = i; k < g->nvert; k++) {
		    if (g->period->L2Gmap_vert[g->period->ordering[k]] != m)
			break;
		    if ((g->types_vert[g->period->ordering[k]] & OWNER))
			l = k;
		}
		/* change local owner to a referenced one (if one exists) */
		while (l < k &&
			g->types_vert[g->period->ordering[l]] == UNREFERENCED)
		    l++;
		if (l >= k) {
		    /* all vertices in the group are unreferenced */
		    l = i;
		}
		
		/* set D2Lmap, the non-owner entries are negated. */
		count = -1;
		for (; i < k; i++) {
		    int cnt;
		    cnt = DofDataCountOnVertex(u, g->period->ordering[i]);
		    /* all the vertices must have the same data count */
		    assert(count == -1 || count == cnt);
		    count = cnt;
		    j = DofVertexDataIndex(u, g->period->ordering[i]) + dindex;
		    if (i == l) {
			for (jj = 0; jj < count; jj++)
			    map->D2Lmap[j + jj] = lindex + jj;
			if (!(g->types_vert[g->period->ordering[i]] & OWNER))
			    (*nonowned_count_vert)++;
		    }
		    else {
			for (jj = 0; jj < count; jj++)
			    map->D2Lmap[j + jj] = -1 - lindex - jj;
		    }
		}
		lindex += count;
	    }
	    dindex += DofVertexDataCount(u);
	}

	/* edge DOF */
	if (type->np_edge > 0) {
	    for (i = 0; i < g->nedge; i++) {
		count = DofDataCountOnEdge(u, i);
		for (jj = 0; jj < count; jj++)
		    map->D2Lmap[dindex++] = lindex++;
	    }
	}

	/* face DOF */
	if (type->np_face > 0) {
	    for (i = 0; i < g->nface; i++) {
		count = DofDataCountOnFace(u, i);
		for (jj = 0; jj < count; jj++)
		    map->D2Lmap[dindex++] = lindex++;
	    }
	}

	/* elem DOF */
	if (type->np_elem > 0) {
	    for (i = 0; i < g->nelem; i++) {
		count = DofDataCountOnElement(u, i);
		for (jj = 0; jj < count; jj++)
		    map->D2Lmap[dindex++] = lindex++;
	    }
	}
    }

    assert(map->localsize == dindex);
    map->nlocal = map->localsize = lindex;
    map->D2Lmap_size = dindex;

    return;
}

static void
build_maps(MAP *map)
{
    int ii;
    DOF *u = map->dofs[0];
    GRID *g = u->g;
    INT i, j, nonowned_count_vert;
#if USE_MPI
    INT dindex, lindex, vindex;
    BOOLEAN vert_flag, edge_flag, face_flag, elem_flag;
    int rank = 0, jj, count, n;
    int *scnts = NULL, *sdsps = NULL, *rcnts = NULL, *rdsps = NULL;
    INT oindex, nsend, nrecv, *sbuf, *rbuf, *buf;
    INT *index_vert, *index_edge, *index_face, *index_elem;
    int *rank_vert, *rank_edge, *rank_face, *rank_elem;
    /* list for recording information abbbox non owned entries, used
     * to avoid multiple traversal of vertices/edges/faces/elements */
    struct {
	INT	*lindices;	/* indices in L2Vmap[] */
	INT	index;		/* vertex/edge/face/element index */
	GTYPE	gtype;		/* geometric type */
    } *nonowned = NULL, *no = NULL;
    INT *nonowned_buffer = NULL, nonowned_count = 0;
    BOOLEAN *flags_vert, *flags_edge, *flags_face, *flags_elem;
    BOOLEAN have_unreferenced;
    DOF_TYPE *type;
#endif	/* USE_MPI */

    map->O2Gmap = map->L2Vmap = map->D2Lmap = NULL;
    map->partition = phgAlloc((g->nprocs + 1) * sizeof(*map->partition));
    map->offset = phgAlloc(((int)map->ndof + 1) * sizeof(*map->offset));
    map->offset[0] = 0;
    for (ii = 0; ii < map->ndof; ii++) {
	map->offset[ii + 1] = map->offset[ii] + DofDataCount(map->dofs[ii]);
    }

    build_periodic_D2Lmap(map, &nonowned_count_vert);

    if (g->nprocs == 1) {
	map->partition[0] = 0;
	map->partition[1] = map->nlocal = map->nglobal = map->localsize;
	map->lif = 1.0;
	if (map->D2Lmap != NULL && map->D2Lmap_size > map->localsize) {
	    /* restore negated indices in D2Lmap */
	    for (i = 0; i < map->D2Lmap_size; i++) {
		if ((j = map->D2Lmap[i]) < 0) {
		    map->D2Lmap[i] = -1 - j;
		    assert((j = map->D2Lmap[i]) >= 0 && j < map->localsize);
		}
	    }
	}
	map->L2Vmap_size = map->localsize;
	return;
    }

#if USE_MPI

    /* TODO: use the following (simpler) procedure for building the map:
     * 		for (ii = 0; ii < map->ndof; ii++) {
     * 		    dof = map->dofs[ii];
     * 		    ForAllElements(g, e) {
     * 			nbas = DofNBas(dof, e);
     * 		        for (jj = 0; jj < nbas; jj++) {
     * 		            phgDofGetElementBasisInfo(dof, e, jj, &gtype,
     * 		            				&gindex, &bindex);
     *			    for (kk = 0; kk < dof->dim; kk++) {
     *				if (the gtype is not owned by present process)
     *				    off_proc entry;
     *				else
     *				    local entry;
     *			    }
     * 		        }
     * 		    }
     * 		}
     * and merge the function sfc_ordering() here.
     */

    assert(g->nprocs > 1);
    vert_flag = edge_flag = face_flag = elem_flag = FALSE;
    for (ii = 0; ii < map->ndof; ii++) {
	u = map->dofs[ii];
	if ((type = u->type) == NULL)
	    type = u->hp->info->types[u->hp->max_order];
	if (type->np_vert > 0)
	    vert_flag = TRUE;
	if (type->np_edge > 0)
	    edge_flag = TRUE;
	if (type->np_face > 0)
	    face_flag = TRUE;
	if (type->np_elem > 0)
	    elem_flag = TRUE;
    }

    index_vert = g->owner_index_vert;
    index_edge = g->owner_index_edge;
    index_face = g->owner_index_face;
    index_elem = g->owner_index_elem;

    rank_vert = g->owner_rank_vert;
    rank_edge = g->owner_rank_edge;
    rank_face = g->owner_rank_face;
    rank_elem = g->owner_rank_elem;

    if (!vert_flag)
	nonowned_count_vert = 0;

    /* create owner rank/index mapping if not available */
    if (face_flag) {
	if (index_face == NULL) {
	    assert(rank_face == NULL);
	    index_face = phgAlloc(g->nface * sizeof(*index_face));
	    rank_face = phgAlloc(g->nface * sizeof(*rank_face));
	    for (i = 0; i < g->nface; i++) {
		index_face[i] = -1;
		rank_face[i] = -1;
	    }
	    for (i = 0; i < g->neighbours.count; i++) {
		RNEIGHBOUR *rn = g->neighbours.list + i;
		lindex = rn->local->faces[rn->vertex];
		if ((g->types_face[lindex] & OWNER))
		    continue;
		index_face[lindex] = rn->peer_face;
		rank_face[lindex] = rn->rank;
	    }
	}
	else {
	    assert(rank_face != NULL);
	}
    }

    scnts = phgAlloc(4 * g->nprocs * sizeof(*scnts));
    sdsps = scnts + g->nprocs;
    rcnts = sdsps + g->nprocs;
    rdsps = rcnts + g->nprocs;
    memset(scnts, 0, 2 * g->nprocs * sizeof(*scnts));

    nonowned_count = nonowned_count_vert +
		     (edge_flag ? g->nedge - g->nedge_owned : 0) +
		     (face_flag ? g->nface - g->nface_owned : 0) +
		     (elem_flag ? g->nelem - g->nelem_owned : 0);
    nonowned = phgAlloc(nonowned_count * sizeof(*nonowned));
    nonowned_buffer = phgAlloc(nonowned_count * map->ndof
					* sizeof(*nonowned_buffer));
    for (i = 0; i < nonowned_count; i++)
	nonowned[i].lindices = nonowned_buffer + i * map->ndof;

    map->L2Vmap = phgAlloc(map->localsize * sizeof(*map->L2Vmap));
    map->L2Vmap_size = map->localsize;	/* temporary size of L2Vmap/O2Gmap */

    dindex = lindex = vindex = 0;
    oindex = map->localsize;

    flags_vert = (vert_flag ? phgCalloc(g->nvert, sizeof(*flags_vert)) : NULL);
    flags_edge = (edge_flag ? phgCalloc(g->nedge, sizeof(*flags_edge)) : NULL);
    flags_face = (face_flag ? phgCalloc(g->nface, sizeof(*flags_face)) : NULL);
    flags_elem = (elem_flag ? phgCalloc(g->nelem, sizeof(*flags_elem)) : NULL);

    for (ii = 0; ii < map->ndof; ii++) {
	u = map->dofs[ii];
	phgDebug((4, "processing DOF %s:\n", u->name));
	no = nonowned;

	/* vertex DOF */
	if (vert_flag) {
	    for (i = 0; i < g->nvert; i++, dindex += count) {
		count = DofDataCountOnVertex(u, i);
		if (map->D2Lmap != NULL && map->D2Lmap[dindex] < 0)
		    continue;
		if (g->types_vert[i] & OWNER) {
		    phgDebug((4, "vert %d, lindex=%d, vindex=%d\n",
							i, lindex, vindex));
		    if (g->period == NULL) {
			for (jj = 0; jj < count; jj++)
			    map->L2Vmap[lindex++] = vindex++;
		    }
		    else {
			j = map->D2Lmap[dindex];
			for (jj = 0; jj < count; jj++)
			    map->L2Vmap[j++] = vindex++;
			lindex += count;
		    }
		}
		else {
		    if (flags_vert[i] == 0) {
			no->index = i;
			no->gtype = VERTEX;
			if (rank_vert != NULL && (n = rank_vert[i]) >= 0) {
			    assert(scnts[2 * n] < INT_MAX);
			    scnts[2 * n]++;
			}
			flags_vert[i] = 1;
		    }
		    if (g->period == NULL) {
			no->lindices[ii] = lindex;
			for (jj = 0; jj < count; jj++)
			    map->L2Vmap[lindex++] = oindex++;
		    }
		    else {
			j = map->D2Lmap[dindex];
			no->lindices[ii] = j;
			for (jj = 0; jj < count; jj++)
			    map->L2Vmap[j++] = oindex++;
			lindex += count;
		    }
		    assert(no->gtype == VERTEX);
		    no++;
		}
	    }
	    assert((INT)(no - nonowned) == nonowned_count_vert);
	}
	else {
	    no += nonowned_count_vert;
	}

	/* edge DOF */
	if (edge_flag) {
	    for (i = 0; i < g->nedge; i++) {
		count = DofDataCountOnEdge(u, i);
		if (g->types_edge[i] & OWNER) {
		    phgDebug((4, "edge %d, lindex=%d, vindex=%d\n",
							i, lindex, vindex));
		    for (jj = 0; jj < count; jj++)
			map->L2Vmap[lindex++] = vindex++;
		}
		else {
		    if (flags_edge[i] == 0) {
			no->index = i;
			no->gtype = EDGE;
			if (rank_edge != NULL && (n = rank_edge[i]) >= 0) {
			    assert(scnts[2 * n] < INT_MAX);
			    scnts[2 * n]++;
			}
			flags_edge[i] = 1;
		    }
		    no->lindices[ii] = lindex;
		    assert(no->gtype == EDGE);
		    no++;
		    for (jj = 0; jj < count; jj++)
			map->L2Vmap[lindex++] = oindex++;
		}
	    }
	    dindex += DofEdgeDataCount(u);
	}

	/* face DOF */
	if (face_flag) {
	    for (i = 0; i < g->nface; i++) {
		count = DofDataCountOnFace(u, i);
		if (g->types_face[i] & OWNER) {
		    phgDebug((4, "face %d, lindex=%d, vindex=%d\n",
							i, lindex, vindex));
		    for (jj = 0; jj < count; jj++)
			map->L2Vmap[lindex++] = vindex++;
		}
		else {
		    if (flags_face[i] == 0) {
			no->index = i;
			no->gtype = FACE;
			if (rank_face != NULL && (n = rank_face[i]) >= 0) {
			    assert(scnts[2 * n] < INT_MAX);
			    scnts[2 * n]++;
			}
			flags_face[i] = 1;
		    }
		    no->lindices[ii] = lindex;
		    assert(no->gtype == FACE);
		    no++;
		    for (jj = 0; jj < count; jj++)
			map->L2Vmap[lindex++] = oindex++;
		}
	    }
	    dindex += DofFaceDataCount(u);
	}

	/* elem DOF */
	if (elem_flag) {
	    for (i = 0; i < g->nelem; i++) {
		count = DofDataCountOnElement(u, i);
		if (g->types_elem[i] & OWNER) {
		    phgDebug((4, "elem %d, lindex=%d, vindex=%d\n",
							i, lindex, vindex));
		    for (jj = 0; jj < count; jj++)
			map->L2Vmap[lindex++] = vindex++;
		}
		else {
		    if (flags_elem[i] == 0) {
			no->index = i;
			no->gtype = VOLUME;
			if (rank_elem != NULL && (n = rank_elem[i]) >= 0) {
			    assert(scnts[2 * n] < INT_MAX);
			    scnts[2 * n]++;
			}
			flags_elem[i] = 1;
		    }
		    no->lindices[ii] = lindex;
		    assert(no->gtype == VOLUME);
		    no++;
		    for (jj = 0; jj < count; jj++)
			map->L2Vmap[lindex++] = oindex++;
		}
	    }
	    dindex += DofElementDataCount(u);
	}
    }
    phgFree(flags_vert);
    phgFree(flags_edge);
    phgFree(flags_face);
    phgFree(flags_elem);
    map->nlocal = vindex;
    assert(map->D2Lmap_size == dindex);
    assert(map->localsize == lindex);
    assert(map->nlocal <= map->localsize);
    assert(no - nonowned <= nonowned_count);
    nonowned_count = no - nonowned;

    /* restore negated indices in D2Lmap */
    if (map->D2Lmap != NULL && map->D2Lmap_size > map->localsize) {
	for (i = 0; i < map->D2Lmap_size; i++) {
	    if ((j = map->D2Lmap[i]) < 0) {
		map->D2Lmap[i] = -1 - j;
		assert((j = map->D2Lmap[i]) >= 0 && j < map->localsize);
	    }
	}
    }

    map->O2Gmap = phgAlloc((map->localsize - map->nlocal)
					* sizeof(*map->O2Gmap));

    /* put nlocal in scnts[2k+1] */
    for (ii = 0; ii < g->nprocs; ii++) {
	scnts[2 * ii + 1] = map->nlocal;	assert(map->nlocal <= INT_MAX);
    }
    MPI_Alltoall(scnts, 2, MPI_INT, rcnts, 2, MPI_INT, g->comm);

    /* now scnts[2k]/rcnts[2k] contain the send/recv counts,
     * scnts[2k+1]/rcnts[2k+1] contain nlocal's for all processes */
    map->partition[0] = 0;
    i = 0;
    for (ii = 0; ii < g->nprocs; ii++) {
	if (i < rcnts[2 * ii + 1])
	    i = rcnts[2 * ii + 1];
	map->partition[ii + 1] = map->partition[ii] + rcnts[2 * ii + 1];
	scnts[ii] = scnts[2 * ii];
	rcnts[ii] = rcnts[2 * ii];
    }
    if (g->period != NULL)
	map->nglobal = map->partition[g->nprocs];
    else if (map->nglobal != map->partition[g->nprocs])
	phgError(-1, "%s:%d: map->nglobal (%"dFMT") != "
		     "map->partition[g->nprocs] (%"dFMT") "
		     "(see option -allow_degenerate_neighbours).\n",
		     __FILE__, __LINE__,
		     map->nglobal, map->partition[g->nprocs]);
    if (phgSameProcessorSpeed()) {
	map->lif = (i * (FLOAT)g->nprocs) / map->nglobal;
    }
    else {
	double w0[3], w[3];
	w0[0] = w0[1] = w0[2] = map->nlocal / phgGetProcessorSpeed(NULL);
	MPI_Allreduce(w0, w, 1, MPI_3DOUBLE, MPI_MSM, g->comm);
	map->lif = w[2] * g->nprocs / w[1];
    }

    nsend = nrecv = 0;
    for (ii = 0; ii < g->nprocs; ii++) {
	sdsps[ii] = nsend;
	rdsps[ii] = nrecv;
	nsend += scnts[ii];	assert(nsend <= INT_MAX);
	nrecv += rcnts[ii];	assert(nrecv <= INT_MAX);
    }

    /* allocate max(ndof,2) INTs for each entry in sbuf/rbuf, which contains:
     *     send: INT[0] = the index in the owner process, INT[1] = gtype
     *     recv: the global indices, one for each DOF */
    n = (map->ndof <= 1 ? 2 : map->ndof);
    sbuf = phgAlloc((nsend + nrecv) * n * sizeof(*sbuf));
    rbuf = sbuf + nsend * n;

    /* set up send buffer */
    for (i = 0, no = nonowned; i < nonowned_count; i++, no++) {
	switch (no->gtype) {
	    case VERTEX:
		if (rank_vert == NULL || (rank = rank_vert[no->index]) == -1)
		    break;
		j = n * (sdsps[rank]++);
		sbuf[j + 0] = index_vert[no->index];
		sbuf[j + 1] = no->gtype;
		/*phgDebug((4, "send[%d]: vert %d ==> vert %d on proc %d\n",
			     sdsps[rank] - 1, no->index, sbuf[ii + 0], rank));*/
		break;
	    case EDGE:
		if (rank_edge == NULL || (rank = rank_edge[no->index]) == -1)
		    break;
		j = n * (sdsps[rank]++);
		sbuf[j + 0] = index_edge[no->index];
		sbuf[j + 1] = no->gtype;
		/*phgDebug((4, "send[%d]: edge %d ==> edge %d on proc %d\n",
			     sdsps[rank] - 1, no->index, sbuf[ii + 0], rank));*/
		break;
	    case FACE:
		if (rank_face == NULL || (rank = rank_face[no->index]) == -1)
		    break;
		j = n * (sdsps[rank]++);
		sbuf[j + 0] = index_face[no->index];
		sbuf[j + 1] = no->gtype;
		/*phgDebug((4, "send[%d]: face %d ==> face %d on proc %d\n",
			     sdsps[rank] - 1, no->index, sbuf[ii + 0], rank));*/
		break;
	    case VOLUME:
		if (rank_elem == NULL || (rank = rank_elem[no->index]) == -1)
		    break;
		j = n * (sdsps[rank]++);
		sbuf[j + 0] = index_elem[no->index];
		sbuf[j + 1] = no->gtype;
		/*phgDebug((4, "send[%d]: elem %d ==> elem %d on proc %d\n",
				no->index, sbuf[ii + 0], rank));*/
		break;
	    default:
		phgError(1, "%s:%d, unexpected.\n", __FILE__, __LINE__);
	}
    }
    phgDebug((4, "nsend = %"dFMT", nrecv = %"dFMT"\n", nsend, nrecv));
    assert(sdsps[g->nprocs - 1] == nsend);

    /* restore sdsps */
    for (ii = 0; ii < g->nprocs; ii++)
	sdsps[ii] -= scnts[ii];

    {
	/* construct MPI Datatype. Since only the 1st two bytes need to be sent,
	 * datatype is defined as 2 consequent INTs, with an extent of n INTs */
	MPI_Datatype type;
#if 1	/* MPI-2 */
	MPI_Datatype type0;
	MPI_Type_contiguous(2, PHG_MPI_INT, &type0);
	MPI_Type_create_resized(type0, 0, (MPI_Aint)(n * sizeof(INT)), &type);
	MPI_Type_free(&type0);
#else	/* MPI-1 */
	int blocklens[2];
	MPI_Aint indices[2];
	MPI_Datatype old_types[2];

	blocklens[0]	= 2;
	indices[0]	= 0;
	old_types[0]	= PHG_MPI_INT;

	blocklens[1]	= 1;
	indices[1]	= n * sizeof(INT);
	old_types[1]	= MPI_UB;

	MPI_Type_struct(2, blocklens, indices, old_types, &type);
#endif
	MPI_Type_commit(&type);
	MPI_Alltoallv(sbuf, scnts, sdsps, type, rbuf, rcnts, rdsps, type,
		      g->comm);
	MPI_Type_free(&type);
    }

    /* process received lists */
    rank = 0;			/* source rank, for debugging */
    for (i = 0, buf = rbuf; i < nrecv; i++, buf += n) {
	GTYPE gtype = buf[1];
	assert(buf[1] >= VERTEX && buf[1] <= VOLUME);
	lindex = buf[0];	/* local index of vertex/edge/face/element */
	while (rank < g->nprocs && i >= rdsps[rank] + rcnts[rank])
	    rank++;
	phgDebug((4, "Processing recv[%d] from p%d, gtype=%d, lindex=%d:\n",
			i, rank, gtype, lindex));
	assert(rank < g->nprocs);
	for (ii = 0; ii < map->ndof; ii++) {
	    u = map->dofs[ii];
	    phgDebug((4, "    DOF %s\n", u->name));
	    switch (gtype) {
		case VERTEX:
		    assert((g->types_vert[lindex] & OWNER));
		    count = DofDataCountOnVertex(u, lindex);
		    vindex = map->offset[ii] + DofVertexDataIndex(u, lindex);
		    break;
		case EDGE:
		    assert((g->types_edge[lindex] & OWNER));
		    count = DofDataCountOnEdge(u, lindex);
		    vindex = map->offset[ii] + DofVertexDataCount(u)
					     + DofEdgeDataIndex(u, lindex);
		    break;
		case FACE:
		    assert((g->types_face[lindex] & OWNER));
		    count = DofDataCountOnFace(u, lindex);
		    vindex = map->offset[ii] + DofVertexDataCount(u)
					     + DofEdgeDataCount(u)
					     + DofFaceDataIndex(u, lindex);
		    break;
		case VOLUME:
		    assert((g->types_elem[lindex] & OWNER));
		    count = DofDataCountOnElement(u, lindex);
		    vindex = map->offset[ii] + DofVertexDataCount(u)
					     + DofEdgeDataCount(u)
					     + DofFaceDataCount(u)
					     + DofElementDataIndex(u, lindex);
		    break;
		default:
		    count = 0;		/* avoid gcc warning */
	    }
	    if (count == 0) {
		phgDebug((4, "\tno DOF data.\n"));
		buf[ii] = -1;
		continue;
	    }
	    if (map->D2Lmap != NULL) {
		assert(vindex >= 0 && vindex < map->D2Lmap_size);
		vindex = map->D2Lmap[vindex];
	    }
	    assert(vindex >= 0 && vindex < map->localsize);
	    phgDebug((4, "\tdof index=%d, vindex=%d, count=%d\n",
			  vindex, map->L2Vmap[vindex], count));
	    vindex = map->L2Vmap[vindex];
	    assert(vindex < map->nlocal);	/* must be owned by the proc */
	    buf[ii] = vindex + map->partition[g->rank];
	}
    }

    /* send back lists */
    {
	/* construct MPI Datatype. Since only the first two bytes need to be sent,
	 * datatype is defined as two consequent INTs, with an scaleent of n INTs */
	MPI_Datatype type;

	MPI_Type_contiguous(n, PHG_MPI_INT, &type);
	MPI_Type_commit(&type);
	MPI_Alltoallv(rbuf, rcnts, rdsps, type, sbuf, scnts, sdsps, type,
		      g->comm);
	MPI_Type_free(&type);
    }

    /* process resulting buffer */
    have_unreferenced = FALSE;
    for (i = 0, no = nonowned; i < nonowned_count; i++, no++) {
	for (ii = 0; ii < map->ndof; ii++) {
	    u = map->dofs[ii];
	    switch (no->gtype) {
		case VERTEX:
		    count = DofDataCountOnVertex(u, no->index);
		    rank = (rank_vert == NULL ? -1 : rank_vert[no->index]);
		    break;
		case EDGE:
		    count = DofDataCountOnEdge(u, no->index);
		    rank = (rank_edge == NULL ? -1 : rank_edge[no->index]);
		    break;
		case FACE:
		    count = DofDataCountOnFace(u, no->index);
		    rank = (rank_face == NULL ? -1 : rank_face[no->index]);
		    break;
		case VOLUME:
		    count = DofDataCountOnElement(u, no->index);
		    rank = (rank_elem == NULL ? -1 : rank_elem[no->index]);
		    break;
		default:
		    count = 0;		/* make gcc happy */
		    rank = -1;
		    phgError(-1, "%s:%d, unexpected.\n", __FILE__, __LINE__);
	    }
	    if (rank >= 0) {
		vindex = sbuf[n * (size_t)sdsps[rank] + ii];
		assert(vindex >= 0 || count == 0);
	    }
	    else {
		have_unreferenced = TRUE;
	    }
	    lindex = no->lindices[ii];
	    for (jj = 0; jj < count; jj++) {
		map->L2Vmap[lindex] -= map->localsize - map->nlocal;
		oindex = map->L2Vmap[lindex++] - map->nlocal;
		if (rank >= 0)
		    map->O2Gmap[oindex] = vindex++;
		else
		    map->O2Gmap[oindex] = -1;
	    }
	}
	assert(sdsps[rank] < INT_MAX);
	sdsps[rank]++;
    }

    phgFree(nonowned_buffer);
    phgFree(nonowned);
    phgFree(sbuf);
    phgFree(scnts);
    if (index_face != g->owner_index_face) {
	phgFree(index_face);
	phgFree(rank_face);
    }

    if (!have_unreferenced)
	return;

    /* adjust maps by removing unreferenced off-proc entries */
    nsend = map->localsize - map->nlocal;
    buf = phgAlloc(nsend * sizeof(*buf));
    for (i = 0, oindex = 0; i < nsend; i++) {
	if ((vindex = map->O2Gmap[i]) >= 0) {
	    buf[i] = oindex;
	    map->O2Gmap[oindex++] = vindex;
	}
	else {
	    buf[i] = -1;
	}
    }

    if (g->period == NULL) {
	/* setup map->L2Vmap and map->D2Lmap */
	assert(map->D2Lmap == NULL);
	map->D2Lmap = phgAlloc(map->D2Lmap_size * sizeof(*map->D2Lmap));
	for (i = 0, lindex = 0; i < map->D2Lmap_size; i++) {
	    if ((vindex = map->L2Vmap[i]) < map->nlocal) {
		map->D2Lmap[i] = lindex;
		map->L2Vmap[lindex++] = vindex;
	    }
	    else {
		if ((vindex = buf[vindex - map->nlocal]) >= 0) {
		    map->D2Lmap[i] = lindex;
		    map->L2Vmap[lindex++] = vindex + map->nlocal;
		}
		else {
		    map->D2Lmap[i] = -1;
		}
	    }
	}
    }
    else {
	INT *L2Lmap;
	assert(map->D2Lmap != NULL);
	/* setup map->L2Vmap and save temporary local==>local map in L2Lmap */
	L2Lmap = phgAlloc(map->L2Vmap_size * sizeof(*L2Lmap));
	for (i = 0, lindex = 0; i < map->L2Vmap_size; i++) {
	    if ((vindex = map->L2Vmap[i]) < map->nlocal) {
		L2Lmap[i] = lindex;
		map->L2Vmap[lindex++] = vindex;
	    }
	    else {
		if ((vindex = buf[vindex - map->nlocal]) >= 0) {
		    L2Lmap[i] = lindex;
		    map->L2Vmap[lindex++] = vindex + map->nlocal;
		}
		else {
		    L2Lmap[i] = -1;
		}
	    }
	}

	/* update map->D2Lmap */
	for (i = 0; i < map->D2Lmap_size; i++)
	    map->D2Lmap[i] = L2Lmap[map->D2Lmap[i]];

	phgFree(L2Lmap);
    }

    /* update map->localsize */
    assert(lindex == map->nlocal + oindex);
    map->localsize = lindex;
    if (map->L2Vmap_size != lindex) {
	assert(lindex < map->L2Vmap_size);
	map->O2Gmap = phgRealloc_(map->O2Gmap,
			(lindex - map->nlocal) * sizeof(*map->O2Gmap),
			(lindex - map->nlocal) * sizeof(*map->O2Gmap));
	map->L2Vmap = phgRealloc_(map->L2Vmap,
			lindex * sizeof(*map->O2Gmap),
			lindex * sizeof(*map->O2Gmap));
	map->L2Vmap_size = map->localsize;
    }

    phgFree(buf);
#endif	/* USE_MPI */
}

static void
check_halo(MAP *map, const char *func)
/* checks whether the halo type has changed after map creation */
{
    GRID *g;

    if (map->dofs == NULL)
	return;

    g = map->dofs[0]->g;
    if (g->halo == NULL && map->halo_type < 0)
	return;
    if (g->halo != NULL && g->halo->type == map->halo_type)
	return;
    phgError(1, "%s: halo type changed after creation of the map.\n", func);
}

INT
phgMapD2L_(MAP *map, int dof_no, INT index, BOOLEAN allow_bogus_Lmap)
{
    if (map->bogus_Lmap == TRUE && allow_bogus_Lmap == FALSE)
	phgError(1, "You shouldn't use local indices with this map."
		    " Please use global indices instead.\n");
    assert(map->dofs != NULL);
    check_halo(map, __func__);
    return map->D2Lmap == NULL ?
		map->offset[dof_no] + index : 
		map->D2Lmap[map->offset[dof_no] + index];
}

INT
phgMapL2G(MAP *map, INT index)
{
    DOF *u = map->dofs[0];
    INT i;

    check_halo(map, __func__);
    if (index >= map->L2Vmap_size) {
	phgError(1, "%s:%d: unexpected: index=%d, L2Vmap_size=%d\n",
			__FILE__, __LINE__, index, map->L2Vmap_size);
    }

    if (index < 0 || (u->g->nprocs == 1 && map->L2Vmap == NULL))
	return index;

    i = (map->L2Vmap == NULL ? index : map->L2Vmap[index]);
    if (map->bogus_Lmap == TRUE)
	return i;
    assert(i < map->localsize);

    if (i < 0)
	return i;
    else if (i < map->nlocal)
	return map->partition[map->rank] + i;
    else
	return map->O2Gmap[i - map->nlocal];
}

static void
get_bdry_map(MAP *map)
/* initialize map->bdry_* members */
{
    GRID *g;
    DOF *dof;
    int ii;
    INT n, *counts, count;
    INT i, j, k, index, index0;

    BTYPE *tmp_types[4];
    int tmp_count[4];
    INT tmp_size[4], *tmp_index[4];
    BTYPE DB_mask;

    if (map->dofs == NULL)
	return;

    g = map->dofs[0]->g;

    for (i = 0; i < map->ndof; i++)
	if (map->dofs[i]->DB_masks == NULL) {
	    if (map->dofs[i]->DB_mask != 0)
		break;
	} else {
	    for (j = 0; j < map->dofs[i]->dim; j++)
		if (map->dofs[i]->DB_masks[j] != 0)
		    break;
	    if (j < map->dofs[i]->dim)
		break;
	}
    if (i >= map->ndof)
	return;

    tmp_types[0] = g->types_vert;
    tmp_types[1] = g->types_edge;
    tmp_types[2] = g->types_face;
    tmp_types[3] = g->types_elem;

    tmp_size[0] = g->nvert;
    tmp_size[1] = g->nedge;
    tmp_size[2] = g->nface;
    tmp_size[3] = g->nelem;

    /* count number of boundary DOF owned by the present process.
     * Note: by putting the loop on map->dofs as the inner most one,
     * then entries in map->bdry[] are increasing */
    count = 0;
    phgFree(map->bdry);
    map->bdry = NULL;
    map->bdry_nlocal = 0;
    for (i = 0; i < map->ndof; i++) {
	dof = map->dofs[i];
	DB_mask = dof->DB_mask;
	if (dof->DB_masks == NULL) {
	    if (DB_mask == 0)
		continue;
	} else {
	    for (k = 0; k < dof->dim; k++)
		if (dof->DB_masks[k] != 0)
		    break;
	    if (k >= dof->dim)
		continue;
	}
	if (!DofIsHP(dof)) {
	    tmp_count[0] = dof->count_vert;
	    tmp_count[1] = dof->count_edge;
	    tmp_count[2] = dof->count_face;
	    tmp_count[3] = dof->count_elem;
	    tmp_index[0] = tmp_index[1] = tmp_index[2] = tmp_index[3] = NULL;
	}
	else {
	    DOF_TYPE *type = dof->hp->info->types[dof->hp->max_order];
	    n = dof->dim;
	    tmp_count[0] = type->np_vert * n;
	    tmp_count[1] = type->np_edge * n;
	    tmp_count[2] = type->np_face * n;
	    tmp_count[3] = type->np_elem * n;
	    tmp_index[0] = dof->hp->vert_index;
	    tmp_index[1] = dof->hp->edge_index;
	    tmp_index[2] = dof->hp->face_index;
	    tmp_index[3] = dof->hp->elem_index;
	}

	index0 = 0;
	for (ii = 0; ii < 4; ii++) {
	    if ((n = tmp_count[ii]) == 0)
		continue;
	    for (j = 0; j < tmp_size[ii]; j++) {
		if (DofIsHP(dof) && tmp_index[ii] != NULL)
		    n = tmp_index[ii][j + 1] - tmp_index[ii][j];

		if (dof->DB_masks == NULL) {
		    if (!(tmp_types[ii][j] & DB_mask)) {
			index0 += n;
			continue; /* nscale geo object */
		    }
		    if (count + n > map->bdry_nlocal) {
			map->bdry_nlocal = (count + n + 1024);
			map->bdry = phgRealloc_(map->bdry,
					map->bdry_nlocal * sizeof(*map->bdry),
					count * sizeof(*map->bdry));
		    }
		    for (index = index0; index < index0 + n; index++) {
			k = phgMapD2L(map, i, index);
			assert(k >= 0);
			k = phgMapL2V(map, k);
			assert(k < map->localsize);
#if 0
			/* only mark local bdry entries */
			if (k >= map->nlocal)
			    continue;
#endif
			map->bdry[count++] = k;
		    }
		} else {
		    for (k = 0; k < dof->dim; k++) { 
			int n0, kk; 
			if (!(tmp_types[ii][j] & dof->DB_masks[k])) {
			    continue; /* nscale dof sub component */
			}
			n0 = n / dof->dim;  
			if (count + n0 > map->bdry_nlocal) {
			    map->bdry_nlocal = (count + n0 + 1024);
			    map->bdry = phgRealloc_(map->bdry,
					  map->bdry_nlocal * sizeof(*map->bdry),
					  count * sizeof(*map->bdry));
			}
			for (index = index0 + k; index < index0 + n;
							index += dof->dim) {
			    kk = phgMapD2L(map, i, index);
			    assert(kk >= 0);
			    kk = phgMapL2V(map, kk);
			    assert(kk < map->localsize);
#if 0
			    /* only mark local bdry entries */
			    if (kk >= map->nlocal)
				continue;
#endif
			    map->bdry[count++] = kk;
			} /* end of dof sub index */
		    }	  /* end of dof dim */
		}	  /* end of dof DB_mask[s] */
		index0 += n;
	    }             /* end of geo count */
	}                 /* end of geo type */
    }			  /* end of get ndof */

    if (count > 0) {
	/* sort entries in bdry[] */
	qsort(map->bdry, count, sizeof(*map->bdry), phgCompINT);
	if (g->period != NULL) {
	    /* delete duplicate entries in map->bdry[] */
	    for (i = 0, k = 0; i < count;) {
		for (j = i + 1; j < count && map->bdry[j] == map->bdry[i]; j++);
		map->bdry[k++] = map->bdry[i];
		i = j;
	    }
	    count = k;
	}
    }

    map->bdry_localsize = count;
    for (i = 0; i < count && map->bdry[i] < map->nlocal; i++);
    map->bdry_nlocal = count = i;

    /* gather counts from all processes */
    counts = phgAlloc(g->nprocs * sizeof(*counts));
#if USE_MPI
    MPI_Allgather(&count, 1, PHG_MPI_INT, counts, 1, PHG_MPI_INT, g->comm);
    k = 0;
    for (i = 0; i < g->nprocs; i++) {
	map->bdry_partition[i] = k;
	k += counts[i];
    }
    map->bdry_partition[i] = k;
#else /* USE_MPI */
    k = counts[0] = count;
    map->bdry_partition[0] = 0;
    map->bdry_partition[1] = count;
#endif /* USE_MPI */
    map->bdry_nglobal = k;

    phgInfo(1, "nlocal: %d, boundary DOF: %d (%d total)\n",
	    map->nlocal, count, k);

    phgFree(counts);
}

MAP *
phgMapCreateSimpleMap(MPI_Comm comm, INT n, INT N)
{
    MAP *map = phgCalloc(1, sizeof(*map));
    int i;

    assert(n != PHG_DECIDE || N != PHG_DECIDE);

    map->bogus_Lmap = FALSE;
    map->comm = comm;
#if USE_MPI
    MPI_Comm_size(comm, &map->nprocs);
    MPI_Comm_rank(comm, &map->rank);
#else
    map->nprocs = 1;
    map->rank = 0;
#endif	/* USE_MPI */
    map->halo_type = -1;

    if (n == PHG_DECIDE) {
	n = N / map->nprocs;
	if (map->rank < N % map->nprocs)
	    n++;
    }
    assert(N == PHG_DECIDE || n <= N);
    map->nlocal = n;
    map->partition = phgAlloc((map->nprocs + 1) * sizeof(*map->partition));
    map->partition[0] = 0;
#if USE_MPI
    MPI_Allgather(&n, 1, PHG_MPI_INT, map->partition + 1, 1, PHG_MPI_INT,
		  map->comm);
#else	/* USE_MPI */
    map->partition[1] = n;
#endif	/* USE_MPI */
    for (i = 0; i < map->nprocs; i++)
	map->partition[i+ 1] += map->partition[i];
    map->nglobal = map->partition[i];
    if (N == PHG_DECIDE)
	N = map->nglobal;

    if (map->nglobal > N) {
	/* special case: identical serial map on all processes */
	assert(map->nglobal == map->nprocs * N);
	map->nglobal = N;
#if USE_MPI
	map->comm = MPI_COMM_SELF;
#endif	/* USE_MPI */
	map->nprocs = 1;
	map->rank = 0;
	phgFree(map->partition);
	map->partition = phgAlloc(2 * sizeof(*map->partition));
	map->partition[0] = 0;
	map->partition[1] = N;
    }

    assert(map->nglobal == N);
    map->localsize = map->nlocal;

    map->bdry_nlocal = map->bdry_nglobal = map->bdry_localsize = 0;
    map->bdry = NULL;
    map->bdry_partition =
	phgCalloc(map->nprocs + 1, sizeof(*map->bdry_partition));

    return map;
}

#if USE_MPI
static const INT *O2Gmap;

static int
O2Gmap_comp(const void *p0, const void *p1)
{
    INT i;
    return (i = O2Gmap[*(const INT *)p0] - O2Gmap[*(const INT *)p1]) > 0 ?
		1 : (i < 0 ? -1 : 0);
}


void
_phg_build_ordering(size_t size, const INT map[], INT **ordering)
/* build array ordering[] such that map[ordering[i]] is increasing w.r.t. i */
{
    INT i;

    if (size <= 0) {
	*ordering = NULL;
	return;
    }
    if (*ordering != NULL)
	return;

    O2Gmap = map;
    *ordering = phgAlloc(size * sizeof(**ordering));
    for (i = 0; i < size; i++)
	(*ordering)[i] = i;
    qsort(*ordering, size, sizeof(**ordering), O2Gmap_comp);
}
#endif /* USE_MPI */

#include "phg/partition-sfc.h"

static SFC *sfc;

static int
comp_sfc_key(const void *n1, const void *n2)
{
    SFC *m1 = sfc + *((const INT *)n1), *m2 = sfc + *((const INT *)n2);

    return m1->sfc > m2->sfc ? 1 : (m1->sfc == m2->sfc ? 0 : -1);
}

static void
sfc_ordering(MAP *map)
{
    int ii, jj, id, nbas;
    INT i, j, n, new_index, *ordering, *ind, *sfc_index;
    FLOAT temp;
    FLOAT bbox[6];	/* bounding box */
    FLOAT scale[3];	/* scaling factor */
    ELEMENT *e;
    GRID *g = map->dofs[0]->g;
    DOF *dof;

    if (phgVerbosity > 0)
        phgPrintf("phgMapCreate: reorder local unknowns using SFC.\n");

    /* allocate space for SFC */
    sfc = phgAlloc(g->nleaf * sizeof(*sfc));
    sfc_index = phgAlloc(g->nleaf * sizeof(*sfc_index));

    /* get first valid element */
    ForAllElements(g, e)
	break;
    g->i = -1;

    if (e != NULL) {
	bbox[0] = bbox[3] = g->verts[e->verts[0]][0];
	bbox[1] = bbox[4] = g->verts[e->verts[0]][1];
	bbox[2] = bbox[5] = g->verts[e->verts[0]][2];
    }
    else {
	bbox[0] = bbox[1] = bbox[2] = 0.0;
	bbox[3] = bbox[4] = bbox[5] = 1.0;
    }

    n = 0;
    ForAllElements(g, e) {
	static const FLOAT lam[] = {0.25, 0.25, 0.25, 0.25};
	for (ii = 0; ii < NVert; ii++) {
	    if (bbox[0] > (temp = g->verts[e->verts[ii]][0]))
		bbox[0] = temp;
	    if (bbox[3] < temp)
		bbox[3] = temp;

	    if (bbox[1] > (temp = g->verts[e->verts[ii]][1]))
		bbox[1] = temp;
	    if (bbox[4] < temp)
		bbox[4] = temp;

	    if (bbox[2] > (temp = g->verts[e->verts[ii]][2]))
		bbox[2] = temp;
	    if (bbox[5] < temp)
		bbox[5] = temp;
	}
#if FT_PHG == FT_DOUBLE
	phgGeomLambda2XYZ(g, e, lam, sfc[n].co, sfc[n].co + 1, sfc[n].co + 2);
#else	/* FT_PHG = FT_DOUBLE */
	{
	    FLOAT x, y, z;
	    phgGeomLambda2XYZ(g, e, lam, &x, &y, &z);
	    sfc[n].co[0] = x;
	    sfc[n].co[1] = y;
	    sfc[n].co[2] = z;
	}
#endif	/* FT_PHG = FT_DOUBLE */
        sfc_index[n++] = e->index;
    }  /* scan all leaves */

    /* enlarge bounding box */
    temp = (bbox[3] - bbox[0]) * HSFC_EPSILON;
    bbox[3] += temp;
    bbox[0] -= temp;
    temp = (bbox[4] - bbox[1]) * HSFC_EPSILON;
    bbox[4] += temp;
    bbox[1] -= temp;
    temp = (bbox[5] - bbox[2]) * HSFC_EPSILON;
    bbox[5] += temp;
    bbox[2] -= temp;

    /* scaling coordinates, such that all coordinates belong to (0, 1)^3 */
    scale[0] = bbox[3] - bbox[0]; /* length of x axis of bounding box */
    scale[1] = bbox[4] - bbox[1]; /* length of y axis of bounding box */
    scale[2] = bbox[5] - bbox[2]; /* length of z axis of bounding box */

#if 1
    /* keep original aspect ratio */
    if (scale[0] < scale[1])  scale[0] = scale[1];
    if (scale[0] < scale[2])  scale[0] = scale[2];
    if (scale[0] == 0.) scale[0] = 1.;
    scale[2] = scale[1] = scale[0] = 1. / scale[0];
#else
    /* ignore original aspect ratio */
    scale[0] = 1. / scale[0];
    scale[1] = 1. / scale[1];
    scale[2] = 1. / scale[2];
#endif

    for (i = 0; i < n; i++) {
        sfc[i].co[0] = (sfc[i].co[0] - bbox[0]) * scale[0];
        sfc[i].co[1] = (sfc[i].co[1] - bbox[1]) * scale[1];
        sfc[i].co[2] = (sfc[i].co[2] - bbox[2]) * scale[2];
    }

    /* generate inverse SFC */
    if (sfc_type == SFC_NONE) {
	/* FIXME: no need to use sfc[] in this case */
	for (i = 0; i < n; i++)
	    sfc[i].sfc = i;
    }
    else {
	if (!(sfc_type == SFC_HILBERT ? phgSFCInvHilbert3D :
					phgSFCInvMorton3D)(sfc, n)) {
            phgFree(sfc);
            phgPrintf("phgMapCreate: error in generating inverse SFC\n");
            return;
	}
    }

    /* Set up ordering such that sfc[ordering[]].key is increasing */
    ordering = phgAlloc(n * sizeof(*ordering));
    for (i = 0; i < n; i++)
	ordering[i] = i;
    qsort(ordering, n, sizeof(*ordering), comp_sfc_key);

    phgFree(sfc);
    ind = phgAlloc(map->nlocal * sizeof(*ind));

    /* initialize ind[] to -1 */
    for (i = 0; i < map->nlocal; i++)
	ind[i] = -1;

    /* Compute new numbering.
     * Note: must order DOF by DOF to match the block matrix structure */
    new_index = 0;
    for (id = 0; id < map->ndof; id++) {	/* loop over DOFs */
	dof = map->dofs[id];
	/* loop on elements in space-filling order */
	for (i = 0; i < n; i++) {
	    e = g->elems[sfc_index[ordering[i]]];
	    assert(e != NULL);
	    nbas = DofNBas(dof, e);
	    for (jj = 0; jj < nbas * dof->dim; jj++) {
		j = phgMapL2V(map, phgMapE2L(map, id, e, jj));
		assert(j >= 0);
		if (j >= map->nlocal)
		    continue;		/* off-process entry */
		if (ind[j] >= 0)
		    continue;		/* already mapped */
		ind[j] = new_index++;
	    }
	}
    }
    assert(new_index == map->nlocal);

    phgFree(sfc_index);
    phgFree(ordering);

    if (map->L2Vmap == NULL) {
	map->L2Vmap = phgAlloc(map->localsize * sizeof(*map->L2Vmap));
	for (i = 0; i < map->localsize; i++)
	    map->L2Vmap[i] = i;
    }

#if USE_MPI
    /* Note: use 'while' in place of 'if' to allow using 'break' */
    while (map->nprocs > 1) {
	/* update O2Gmap for the new numbering */
	int mask, *scnts, *sdsps, *rcnts, *rdsps;
	INT *sbuf, *rbuf, ssize, rsize, O2Gsize = map->localsize - map->nlocal;

	scnts = phgAlloc(4 * map->nprocs * sizeof(*scnts));
	sdsps = scnts + map->nprocs;
	rcnts = sdsps + map->nprocs;
	rdsps = rcnts + map->nprocs;

	ii = mask = 0;
	bzero(scnts, map->nprocs * sizeof(*scnts));
	for (i = 0; i < O2Gsize; i++) {
	    j = map->O2Gmap[map->ordering == NULL ? i : map->ordering[i]];
	    while (ii < map->nprocs && map->partition[ii + 1] <= j)
		ii++;
	    assert(ii < map->nprocs && ii != map->rank);
	    scnts[ii]++;
	    mask = 1;
	}

	/* attach mask to scnts[] */
	if (mask != 0) {
	    mask = (mask << (sizeof(mask) * 8 - 1));
	    for (ii = 0; ii < map->nprocs; ii++)
		scnts[ii] |= mask;
	}

	MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, map->comm);

	/* remove mask from scnts and rcnts, and set up sdsps and rdsps */
	mask = (1 << (sizeof(mask) * 8 - 1));
	jj = 0;
	ssize = rsize = 0;
	for (ii = 0; ii < map->nprocs; ii++) {
	    if ((rcnts[ii] & mask))
		jj = 1;
	    sdsps[ii] = ssize;
	    rdsps[ii] = rsize;
	    ssize += (scnts[ii] &= ~mask);
	    rsize += (rcnts[ii] &= ~mask);
	}
	assert(ssize == O2Gsize);

	if (jj == 0) {
	    phgFree(scnts);
	    break;		/* no data to communicate */
	}

	if (map->ordering == NULL) {
	    rbuf = phgAlloc(rsize * sizeof(*sbuf));
	    sbuf = map->O2Gmap;
	}
	else {
	    rbuf = phgAlloc((rsize + ssize) * sizeof(*sbuf));
	    sbuf = rbuf + rsize;
	    for (i = 0; i < O2Gsize; i++)
		sbuf[i] = map->O2Gmap[map->ordering[i]];
	}

	MPI_Alltoallv(sbuf, scnts, sdsps, PHG_MPI_INT,
		      rbuf, rcnts, rdsps, PHG_MPI_INT, map->comm);

	/* update rbuf[] */
	n = map->partition[map->rank];
	for (i = 0; i < rsize; i++) {
	    j = rbuf[i] - n;
	    assert(j >= 0 && j < map->nlocal);
	    rbuf[i] = ind[j] + n;
	}

	/* send back rbuf[] */
	MPI_Alltoallv(rbuf, rcnts, rdsps, PHG_MPI_INT,
		      sbuf, scnts, sdsps, PHG_MPI_INT, map->comm);

	if (map->ordering != NULL) {
	    for (i = 0; i < O2Gsize; i++)
		map->O2Gmap[map->ordering[i]] = sbuf[i];
	}

	phgFree(rbuf);
	phgFree(scnts);

	/* rebuild map->ordering */
	phgFree(map->ordering);
	map->ordering = NULL;
	_phg_build_ordering(map->localsize - map->nlocal, map->O2Gmap,
			    &map->ordering);

	break;
    }
#endif	/* USE_MPI */

    for (i = 0; i < map->localsize; i++) {
	if ((j = map->L2Vmap[i]) >= map->nlocal)
	    continue;
	map->L2Vmap[i] = ind[j];
    }
    phgFree(ind);
}

MAP *
phgMapCreateN(int ndof, DOF **dofs)
{
    MAP *map = phgCalloc(1, sizeof(*map));
    GRID *g;
    DOF *u;
    int j;

    assert(dofs != NULL);

    map->bogus_Lmap = FALSE;
    map->ndof = ndof;
    map->dofs = phgAlloc(ndof * sizeof(*map->dofs));
    memcpy(map->dofs, dofs, ndof * sizeof(*map->dofs));
    u = dofs[0];
    g = u->g;

    /* save halo type */
    map->halo_type = (g->halo == NULL ? -1 : g->halo->type);

#if USE_MPI
    map->comm = g->comm;
#endif	/* USE_MPI */
    map->nprocs = g->nprocs;
    map->rank = g->rank;

    map->nlocal = 0;				/* updated in build_maps() */
    map->nglobal = map->localsize = 0;
    for (j = 0; j < map->ndof; j++) {
	u = map->dofs[j];
	map->nglobal += DofDataCountGlobal(u);
	map->localsize += DofDataCount(u);
    }

    build_maps(map);

#if USE_MPI
    _phg_build_ordering(map->localsize - map->nlocal, map->O2Gmap,
			&map->ordering);
#endif	/* USE_MPI */

    /* reorder unknowns based on SFC (for hybrid/l1 GS smoother efficiency).
     * test case:
     * 	env OMP_NUM_THREADS=4 ./poisson -solver hypre -hypre_pc boomeramg \
     * 		-map_reorder {yes or no} */
    if (sfc_reorder == SFC_YES ||
	(sfc_reorder == SFC_AUTO && phgMaxThreads > 1))
	sfc_ordering(map);

    map->bdry_nlocal = map->bdry_nglobal = map->bdry_localsize = 0;
    map->bdry = NULL;
    map->bdry_partition =
	phgCalloc(g->nprocs + 1, sizeof(*map->bdry_partition));
    get_bdry_map(map);

    return map;
}

MAP *
phgMapCreate(DOF *u, ...)
{
    int ndof;
    DOF **dofs;
    MAP *map;
    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	/* Register options */
	initialized = TRUE;
	phgOptionsRegisterKeyword("-map_reorder",
				  "Whether reorder local entries in MAP",
				  sfc_flags, &sfc_reorder);
	phgOptionsRegisterKeyword("-map_reorder_type",
				  "Method for reordering local entries in MAP",
				  sfc_types, &sfc_type);

	return NULL;
    }

    GetVariableArgs(ndof, dofs, u, DOF *);
    map = phgMapCreateN(ndof, dofs);
    phgFree(dofs);

    return map;
}

void
phgMapDestroy(MAP **map_ptr)
{
    MAP *map = *map_ptr;

    if (map == NULL)
	return;

    *map_ptr = NULL;

    phgMapDestroy(&map->map_nb);

    assert(map->refcount >= 0);

    if (map->refcount > 0) {
	map->refcount--;
	return;
    }

    phgFree(map->D2Lmap);
    phgFree(map->L2Vmap);
    phgFree(map->O2Gmap);
    phgFree(map->ordering);
    phgFree(map->partition);
    phgFree(map->offset);
    phgFree(map->bdry);
    phgFree(map->bdry_partition);

    phgFree(map->dofs);

#if USE_MPI
    phgMapDestroyCommInfo(&map->cinfo);
#endif	/* USE_MPI */

    phgFree(map);

    return;
}

MAP *
phgMapRemoveBoundaryEntries_(MAP *map, BOOLEAN dof_flag)
/* derives a new map from 'map' with all boundary entries removed.
 * if 'dof_flag' == TRUE then preserve or build DOF-to-vector map.
 *
 * See examples/eigen.c for an example usage of this function */
{
    INT n, N, i, k;
    INT *new_index = NULL;
#if USE_MPI
    INT m0, n0, O2Gsize, *sbuff, *rbuff;
    int ii, *scounts, *sdispls, *rcounts, *rdispls;
#endif 

    if (map->bdry_nglobal == 0) {
	map->refcount++;
	return map;
    }

    if (map->map_nb != NULL) {
	map->map_nb->refcount++;
	return map->map_nb;
    }
 
    n = map->nlocal - map->bdry_nlocal;
    N = map->nglobal - map->bdry_nglobal;
    assert(map->dofs != NULL);
    map->map_nb = phgMapCreateSimpleMap(map->comm, n, N);
    map->map_nb->refcount++;	/* note: map_nb is referenced by map */

    if (!dof_flag || map->ndof == 0)
	return map->map_nb;

    map->map_nb->D2Lmap_size = map->D2Lmap_size;
    map->map_nb->L2Vmap_size = map->L2Vmap_size;

    /* Copy list of dofs and offsets */
    map->map_nb->ndof = map->ndof;
    map->map_nb->dofs = phgAlloc(map->ndof * sizeof(*map->dofs));
    memcpy(map->map_nb->dofs, map->dofs, map->ndof * sizeof(*map->dofs));

    map->map_nb->offset = phgAlloc(((int)map->ndof + 1) * sizeof(*map->offset));
    memcpy(map->map_nb->offset, map->offset,
			((int)map->ndof + 1) * sizeof(*map->offset));

    /* Copy D2Lmap */
    if (map->D2Lmap != NULL) {
	map->map_nb->D2Lmap = phgAlloc(map->D2Lmap_size * sizeof(*map->D2Lmap));
	memcpy(map->map_nb->D2Lmap, map->D2Lmap,
		map->D2Lmap_size * sizeof(*map->D2Lmap));
    }

    /* Compute L2Vmap (removing boundary entries), this part of code is similar
     * to the code for building new_index[] in phgMatRemoveBoundaryEntries() */

    new_index = phgAlloc(map->nlocal * sizeof(*new_index));

    /* 1. initialize all entries to -1 */
    for (i = 0; i < map->nlocal; i++) {
	new_index[i] = -1;
    }

    /* 2. set boundary entries to -2 */
    for (i = 0; i < map->bdry_nlocal; i++) {
	assert(map->bdry[i] >= 0 && map->bdry[i] < map->nlocal);
	new_index[map->bdry[i]] = -2;
    }

    /* 3. convert non boundary entries to new local indices */
    for (i = 0, k = 0; i < map->nlocal; i++) {
	if (new_index[i] == -1)
	    new_index[i] = k++;
    }
    assert(k == map->map_nb->nlocal);

    /* Now the meaning of new_index[map->L2Vmap[]]:
     *	. >= 0 means a non boundary entry, and is its new local index;
     *	. -1 means an offprocess non boundary entry;
     *	. -2 means a boundary entry */

    if (map->nprocs <= 1) {
	map->map_nb->localsize = map->map_nb->nlocal;
	map->map_nb->L2Vmap = phgAlloc(map->L2Vmap_size * sizeof(*map->L2Vmap));
	if (map->L2Vmap != NULL) {
	    for (i = 0; i < map->L2Vmap_size; i++)
		map->map_nb->L2Vmap[i] = new_index[map->L2Vmap[i]];
	}
	else {
	    for (i = 0; i < map->L2Vmap_size; i++)
		map->map_nb->L2Vmap[i] = new_index[i];
	}
	phgFree(new_index);

	return map->map_nb;
    }

#if USE_MPI
    O2Gsize = map->localsize - map->nlocal;
    sbuff = phgAlloc(O2Gsize * sizeof(*sbuff));
    for (i = 0; i < O2Gsize; i++) {
	sbuff[i] = map->O2Gmap[map->ordering[i]];
#if DEBUG
	if (i > 0)
	    assert(sbuff[i] > sbuff[i - 1]);
#endif	/* DEBUG */
    }

    scounts = phgAlloc(4 * map->nprocs * sizeof(*scounts));
    sdispls = scounts + map->nprocs * 1;
    rcounts = scounts + map->nprocs * 2;
    rdispls = scounts + map->nprocs * 3;

    /* count number of off-process columns for each process */
    memset(scounts, 0, map->nprocs * sizeof(*scounts));
    for (ii = 1, k = 0; k < O2Gsize; k++) {
	i = sbuff[k];
	while (ii <= map->nprocs && map->partition[ii] <= i)
	    ii++;
	assert(ii <= map->nprocs);
	assert(scounts[ii - 1] < INT_MAX);
	scounts[ii - 1]++;
    }

    /* exchange off-process columns counts */
    MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, map->comm);

    /* compute displacements in the send and receive buffers */
    sdispls[0] = rdispls[0] = 0;
    for (ii = 0; ii < map->nprocs - 1; ii++) {
	sdispls[ii + 1] = sdispls[ii] + scounts[ii];
	assert(sdispls[ii + 1] >= sdispls[ii]);
	rdispls[ii + 1] = rdispls[ii] + rcounts[ii];
	assert(rdispls[ii + 1] >= rdispls[ii]);
    }
    assert(O2Gsize == sdispls[ii] + scounts[ii]);
    n = rdispls[ii] + (INT)rcounts[ii];
			assert(n >= rdispls[ii] && n <= INT_MAX);
    rbuff = phgAlloc(n * sizeof(*rbuff));

    /* send indices to owning process */
    MPI_Alltoallv(sbuff, scounts, sdispls, PHG_MPI_INT,
		  rbuff, rcounts, rdispls, PHG_MPI_INT, map->comm);

    /* convert received indices to new indices for map->map_nb */
    n0 = map->partition[map->rank];
    m0 = map->map_nb->partition[map->rank];
    for (k = 0; k < n; k++) {
	i = rbuff[k] - n0;
	assert(i >= 0 && i < map->nlocal);
	if (new_index[i] >= 0)
	    rbuff[k] = new_index[i] + m0;
	else
	    rbuff[k] = -1;
    }

    /* send back converted indices */
    MPI_Alltoallv(rbuff, rcounts, rdispls, PHG_MPI_INT,
		  sbuff, scounts, sdispls, PHG_MPI_INT, map->comm);

    phgFree(scounts);

    if (O2Gsize > 0) {
	map->map_nb->O2Gmap = phgAlloc(O2Gsize * sizeof(*map->O2Gmap));
	for (i = 0; i < O2Gsize; i++)
	    map->map_nb->O2Gmap[map->ordering[i]] = sbuff[i];
	/* compress map->map_nb->O2Gmap[], and set up sbuff[] as old to new
	 * off-proc index map */ 
	for (i = k = 0; i < O2Gsize; i++) {
	    if (map->map_nb->O2Gmap[i] >= 0) {
		sbuff[i] = map->map_nb->nlocal + k;
		map->map_nb->O2Gmap[k++] = map->map_nb->O2Gmap[i];
	    }
	    else {
		sbuff[i] = -1;
	    }
	}
	map->map_nb->localsize = map->map_nb->nlocal + k;
	map->map_nb->O2Gmap = phgRealloc_(map->map_nb->O2Gmap,
					  k * sizeof(*map->O2Gmap),
					  O2Gsize * sizeof(*map->O2Gmap));
	assert(map->map_nb->ordering == NULL);
	_phg_build_ordering(k, map->map_nb->O2Gmap, &map->map_nb->ordering);
    }

    /* build map->map_nb->L2Vmap[] */
    map->map_nb->L2Vmap = phgAlloc(map->L2Vmap_size * sizeof(*map->L2Vmap));
    for (i = 0; i < map->L2Vmap_size; i++) {
	k = (map->L2Vmap != NULL ? map->L2Vmap[i] : i);
	assert(k >= 0 && k < map->localsize);
	if (k < map->nlocal) {
	    map->map_nb->L2Vmap[i] = new_index[k];
	    assert(map->map_nb->L2Vmap[i] < map->map_nb->nlocal);
	}
	else {
	    k = sbuff[k - map->nlocal];
	    map->map_nb->L2Vmap[i] = k; 
	    assert((k >= map->map_nb->nlocal && k < map->map_nb->localsize)
		   || k < 0);
	}
    }

    phgFree(new_index);
    phgFree(rbuff);
    phgFree(sbuff);
#endif /* USE_MPI */

    return map->map_nb;
}

/*--------- functions for vector scattering/gathering ---------*/

#if USE_MPI
COMM_INFO *
phgMapCreateCommInfo(MPI_Comm comm, int nprocs, int rank, INT nlocal,
		     INT localsize, const INT O2Gmap[], const INT ordering[],
		     const INT partition[])
/* Note: ordering may be NULL which means the natural ordering */
{
    COMM_INFO *cinfo;
    INT i, k, n, *work0 = NULL, O2Gsize = localsize - nlocal;
    const INT *work;
    int ii, mask;

    if (nprocs == 1)
	return NULL;

    cinfo = phgCalloc(1, sizeof(*cinfo));
    cinfo->comm = comm;

    cinfo->scnts = phgAlloc(4 * nprocs * sizeof(*cinfo->scnts));
    cinfo->sdsps = cinfo->scnts + nprocs;
    cinfo->rcnts = cinfo->sdsps + nprocs;
    cinfo->rdsps = cinfo->rcnts + nprocs;

    ii = 0;
    bzero(cinfo->rcnts, nprocs * sizeof(*cinfo->rcnts));
    mask = 0;
    if (ordering != NULL) {
	work = work0 = phgAlloc(O2Gsize * sizeof(*work0));
	n = -1;
	for (i = 0; i < O2Gsize; i++) {
	    k = work0[i] = O2Gmap[ordering[i]];
	    if (k <= n) {
		phgError(1, "%s: bad ordering in O2Gmap[ordering[]] "
			    "(%s:%d, i=%d, k=%d, n=%d).\n",
			    __func__, __FILE__, __LINE__, i, k, n);
	    }
	    n = k;
	    while (ii < nprocs && partition[ii + 1] <= k)
		ii++;
	    assert(ii < nprocs && ii != rank && cinfo->rcnts[ii] < INT_MAX);
	    cinfo->rcnts[ii]++;
	    mask = 1;
	}
    }
    else {
	work = O2Gmap;
	n = -1;
	for (i = 0; i < O2Gsize; i++) {
	    k = O2Gmap[i];
	    if (k <= n) {
		phgError(1, "%s: inconsistent ordering in O2Gmap[] "
			    "(%s:%d, i=%d, k=%d, n=%d).\n",
			    __func__, __FILE__, __LINE__, i, k, n);
	    }
	    n = k;
	    while (ii < nprocs && partition[ii + 1] <= k)
		ii++;
	    assert(ii < nprocs && ii != rank && cinfo->rcnts[ii] < INT_MAX);
	    cinfo->rcnts[ii]++;
	    mask = 1;
	}
    }

    /* attach mask to rcnts[] */
    if (mask != 0) {
	mask = (mask << (sizeof(mask) * 8 - 1));
	for (ii = 0; ii < nprocs; ii++)
	    cinfo->rcnts[ii] |= mask;
    }

    MPI_Alltoall(cinfo->rcnts, 1, MPI_INT, cinfo->scnts, 1, MPI_INT, comm);

    /* check and remove mask */
    mask = (1 << (sizeof(mask) * 8 - 1));
    n = 0;
    for (ii = 0; ii < nprocs; ii++) {
	if ((cinfo->scnts[ii] & mask))
	    n = 1;
	cinfo->scnts[ii] &= ~mask;
	cinfo->rcnts[ii] &= ~mask;
    }
    if (n == 0) {
	/* no communication needed for this map */
	assert(work == NULL);
	phgFree(cinfo->scnts);
	cinfo->scnts = cinfo->sdsps = cinfo->rcnts = cinfo->rdsps = NULL;
	cinfo->comm = MPI_COMM_NULL;
	return cinfo;
    }

    cinfo->size = nlocal;
    cinfo->rind = ordering;

    cinfo->rsize = cinfo->ssize = 0;
    for (ii = 0; ii < nprocs; ii++) {
	cinfo->rdsps[ii] = cinfo->rsize;
	cinfo->sdsps[ii] = cinfo->ssize;
	cinfo->rsize += cinfo->rcnts[ii];	assert(cinfo->rsize <= INT_MAX);
	cinfo->ssize += cinfo->scnts[ii];	assert(cinfo->ssize <= INT_MAX);
    }
    assert(cinfo->rsize == O2Gsize);
    /* prepare indices of data to send */
    cinfo->sind = phgAlloc(cinfo->ssize * sizeof(*cinfo->sind));
    MPI_Alltoallv((INT *)work, cinfo->rcnts, cinfo->rdsps, PHG_MPI_INT,
		  cinfo->sind, cinfo->scnts, cinfo->sdsps, PHG_MPI_INT, comm);
    phgFree(work0);
    /* convert indices in sind[] to local */
    n = partition[rank];
    for (i = 0; i < cinfo->ssize; i++) {
	if (cinfo->sind[i] < n || cinfo->sind[i] >= partition[rank + 1]) {
	    phgError(1, "unexpected error in %s (%s:%d).\n",
			__func__, __FILE__, __LINE__);
	}
	cinfo->sind[i] -= n;
    }

    assert(cinfo->rcnts[rank] == 0 && cinfo->scnts[rank] == 0);

    return cinfo;
}

void
phgMapDestroyCommInfo(COMM_INFO **cinfo)
{
    if (*cinfo == NULL)
	return;

    assert((*cinfo)->nsends == 0 && (*cinfo)->nrecvs == 0
		&& (*cinfo)->buffer == NULL);

    phgFree((*cinfo)->scnts);
    phgFree((*cinfo)->sind);
    phgFree(*cinfo);

    *cinfo = NULL;
}

void
phgMapGatherBegin(COMM_INFO *cinfo, int nvec, FLOAT *data, FLOAT *offp_data)
{
    int ii;
    INT i;
    const INT *rind;
    FLOAT *sbuff, *rbuff, *p, *q;
    MPI_Datatype type;

    if (nvec == 0 || cinfo == NULL || cinfo->comm == MPI_COMM_NULL)
	return;

    assert(cinfo->nsends == 0 && cinfo->nrecvs == 0 && cinfo->buffer == NULL);

#if 0	/* Note: the next return may cause MPI_Alltoallv2 to hang! */
    if (cinfo->ssize == 0 && cinfo->rsize == 0)
	return;
#endif

    if (cinfo->rind != NULL || nvec > 1) {
	/* Note: need rbuff to reorder the data if nvec > 1
	 * (FIXME: use MPI_Datatype with ex=1, size=rsize?) */
	cinfo->buffer = phgAlloc(nvec * (cinfo->ssize + cinfo->rsize)
				      * sizeof(*cinfo->buffer));
	sbuff = cinfo->buffer;
	rbuff = sbuff + nvec * cinfo->ssize;
    }
    else {
	cinfo->buffer = phgAlloc(nvec * cinfo->ssize * sizeof(*cinfo->buffer));
	sbuff = cinfo->buffer;
	rbuff = offp_data;
    }

    if (nvec > 1) {
	MPI_Type_contiguous(nvec, PHG_MPI_FLOAT, &type);
	MPI_Type_commit(&type);
    }
    else {
	type = PHG_MPI_FLOAT;
    }

    /* copy data to send from offp_data to rbuff */
    if (cinfo->rind != NULL) {
	q = rbuff;
	rind = cinfo->rind;
	if (nvec == 1) {
	    /* faster code for nvec == 1 */
	    for (i = 0; i < cinfo->rsize; i++)
		*(q++) = offp_data[*(rind++)];
	}
	else {
	    for (i = 0; i < cinfo->rsize; i++) {
		p = offp_data + *(rind++);
		for (ii = 0; ii < nvec; ii++) {
		    *(q++) = *p;
		    p += cinfo->rsize;
		}
	    }
	}
    }
    else if (nvec > 1) {	/* no need to copy data if nvec == 1 */
	q = rbuff;
	for (i = 0; i < cinfo->rsize; i++) {
	    p = offp_data++;
	    for (ii = 0; ii < nvec; ii++) {
		*(q++) = *p;
		p += cinfo->rsize;
	    }
	}
    }

    MPI_AlltoallvBegin(rbuff, cinfo->rcnts, cinfo->rdsps, type,
		       sbuff, cinfo->scnts, cinfo->sdsps, type, cinfo->comm,
		       &cinfo->requests, &cinfo->nsends, &cinfo->nrecvs);

    if (nvec > 1)
	MPI_Type_free(&type);
}

void
phgMapGatherEnd(COMM_INFO *cinfo, int nvec, FLOAT *data, FLOAT *offp_data)
{
    int ii;
    INT i, *sind;
    FLOAT *p, *q;

    if (nvec == 0 || cinfo == NULL || cinfo->comm == MPI_COMM_NULL)
	return;

#if 0
    if (cinfo->nsends == 0 && cinfo->nrecvs == 0)
	return;
#endif

    MPI_AlltoallvEnd(cinfo->nsends, cinfo->nrecvs, &cinfo->requests);
    cinfo->nsends = cinfo->nrecvs = 0;

    p = cinfo->buffer;
    sind = cinfo->sind;
    if (nvec == 1) {
	/* faster code for nvec == 1 */
	for (i = 0; i < cinfo->ssize; i++) {
	    data[*(sind++)] += *(p++);
	}
    }
    else {
	for (i = 0; i < cinfo->ssize; i++) {
	    q = data + *(sind++);
	    for (ii = 0; ii < nvec; ii++) {
		*q += *(p++);
		q += cinfo->size;
	    }
	}
    }

    phgFree(cinfo->buffer);
    cinfo->buffer = NULL;
}

void
phgMapGather(COMM_INFO *cinfo, int nvec, FLOAT *data, FLOAT *offp_data)
/* data += offp_data */
{
    phgMapGatherBegin(cinfo, nvec, data, offp_data);
    phgMapGatherEnd(cinfo, nvec, data, offp_data);
}

void
phgMapScatterBegin(COMM_INFO *cinfo, int nvec, FLOAT *data, FLOAT *offp_data)
{
    int ii;
    INT i, *sind;
    FLOAT *sbuff, *rbuff, *p, *q;
    MPI_Datatype type;

    if (nvec == 0 || cinfo == NULL || cinfo->comm == MPI_COMM_NULL)
	return;

    assert(cinfo->nsends == 0 && cinfo->nrecvs == 0 && cinfo->buffer == NULL);

#if 0	/* Note: the next return may cause MPI_Alltoallv2 to hang! */
    if (cinfo->ssize == 0 && cinfo->rsize == 0)
	return;
#endif

    if (nvec > 1) {
	static BOOLEAN warned = FALSE;
	if (!warned) {
	    phgWarning("function %s has not been tested!\n", __func__);
	    warned = TRUE;
	}
    }

    if (cinfo->rind != NULL || nvec > 1) {
	/* Note: need rbuff to reorder the data if nvec > 1
	 * (FIXME: use MPI_Datatype with ex=1, size=rsize?) */
	cinfo->buffer = phgAlloc(nvec * (cinfo->ssize + cinfo->rsize)
				      * sizeof(*cinfo->buffer));
	sbuff = cinfo->buffer;
	rbuff = sbuff + nvec * cinfo->ssize;
    }
    else {
	cinfo->buffer = phgAlloc(nvec * cinfo->ssize * sizeof(*cinfo->buffer));
	sbuff = cinfo->buffer;
	rbuff = offp_data;
    }

    if (nvec > 1) {
	MPI_Type_contiguous(nvec, PHG_MPI_FLOAT, &type);
	MPI_Type_commit(&type);
    }
    else {
	type = PHG_MPI_FLOAT;
    }

    /* copy data to send from data to sbuff */
    p = sbuff;
    sind = cinfo->sind;
    if (nvec == 1) {
	/* faster code for nvec == 1 */
	for (i = 0; i < cinfo->ssize; i++) {
	    *(p++) = data[*(sind++)];
	}
    }
    else {
	for (i = 0; i < cinfo->ssize; i++) {
	    q = data + *(sind++);
	    for (ii = 0; ii < nvec; ii++) {
		*(p++) = *q;
		q += cinfo->size;
	    }
	}
    }

    MPI_AlltoallvBegin(sbuff, cinfo->scnts, cinfo->sdsps, type,
		       rbuff, cinfo->rcnts, cinfo->rdsps, type, cinfo->comm,
		       &cinfo->requests, &cinfo->nsends, &cinfo->nrecvs);

    if (nvec > 1)
	MPI_Type_free(&type);
}

void
phgMapScatterEnd(COMM_INFO *cinfo, int nvec, FLOAT *data, FLOAT *offp_data)
{
    int ii;
    INT i;
    const INT *rind;
    FLOAT *p, *q;

    if (nvec == 0 || cinfo == NULL || cinfo->comm == MPI_COMM_NULL)
	return;

#if 0
    if (cinfo->nsends == 0 && cinfo->nrecvs == 0)
	return;
#endif

    MPI_AlltoallvEnd(cinfo->nsends, cinfo->nrecvs, &cinfo->requests);
    cinfo->nsends = cinfo->nrecvs = 0;

    if (cinfo->rind != NULL) {
	q = cinfo->buffer + nvec * cinfo->ssize;
	rind = cinfo->rind;
	if (nvec == 1) {
	    /* faster code for nvec == 1 */
	    for (i = 0; i < cinfo->rsize; i++) {
		offp_data[*(rind++)] = *(q++);
	    }
	}
	else {
	    for (i = 0; i < cinfo->rsize; i++) {
		p = offp_data + *(rind++);
		for (ii = 0; ii < nvec; ii++) {
		    *p = *(q++);
		    p += cinfo->rsize;
		}
	    }
	}
    }
    else if (nvec > 1) {	/* data already in offp_data if nvec == 1 */
	q = cinfo->buffer + nvec * cinfo->ssize;
	for (i = 0; i < cinfo->rsize; i++) {
	    p = offp_data++;
	    for (ii = 0; ii < nvec; ii++) {
		*p = *(q++);
		p += cinfo->rsize;
	    }
	}
    }

    phgFree(cinfo->buffer);
    cinfo->buffer = NULL;
}

void
phgMapScatter(COMM_INFO *cinfo, int nvec, FLOAT *data, FLOAT *offp_data)
/* offp_data := data */
{
    phgMapScatterBegin(cinfo, nvec, data, offp_data);
    phgMapScatterEnd(cinfo, nvec, data, offp_data);
}
#endif	/* USE_MPI */

/*--------- functions for transferring data between DOF list and VEC --------*/

int
phgMapLocalDataToDof(MAP *map, int ndof, DOF **dofs, FLOAT *data)
/* This function copies the solution back to DOFs.
 * 'vec' is the pointer to the local vector data (size = map->nlocal)
 */
{
#if 0	/*==================================================================*/
  /* The old circular algorithm */
  {
    DOF *u;
    INT i, j, k = 0, l, m, nmax;
    INT nlocal = map->nlocal;
    /* struct for storing remote data */
    struct {
	FLOAT value;
	INT gindex;
    } *rdata = NULL;
    INT *offp_counts = NULL;
#if USE_MPI
    MPI_Status status;
    MPI_Datatype type;
    INT first, last;
    int src, dst;

    if (map->nprocs > 1) {
	offp_counts = phgAlloc(map->nprocs * sizeof(*offp_counts));
	m = map->localsize - nlocal;
	MPI_Allgather(&m, 1, PHG_MPI_INT, offp_counts, 1, PHG_MPI_INT,
		      map->comm);
    }
#endif

    nmax = 0;
    if (offp_counts != NULL) {
	for (i = 0; i < map->nprocs; i++) {
	    if (nmax < offp_counts[i])
		nmax = offp_counts[i];
	}
    }
    if (nmax == 0) {
	k = 0;
	goto skip_comm;
    }

#if USE_MPI
    /* TODO: better use Alltoallv */

    /* copy indices of missing data */
    rdata = phgAlloc(2 * nmax * sizeof(*rdata));
    m = map->localsize - map->nlocal;
    for (i = 0; i < m; i++)
	rdata[i].gindex = map->O2Gmap[i];

    /* circulate lists of missing data and collecting missing data */
    j = map->rank;
    k = 0;
    dst = (map->rank == 0) ? map->nprocs - 1 : map->rank - 1;
    src = (map->rank == map->nprocs - 1) ? 0 : map->rank + 1;
    MPI_Type_contiguous(sizeof(*rdata), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    first = map->partition[map->rank];
    last = first + nlocal;
    for (i = 0; i < map->nprocs; i++) {
	if ((l = j + 1) >= map->nprocs)
	    l = 0;
	MPI_Sendrecv(rdata + k, offp_counts[j], type, dst, 111,
		     rdata + (nmax - k), offp_counts[l], type, src,
		     111, map->comm, &status);
	j = l;
	k = nmax - k;
	if (i == map->nprocs - 1)
	    break;
	/* pick up data from local data */
	for (l = 0; l < offp_counts[j]; l++) {
	    if ((m = rdata[k + l].gindex) < first)
		continue;
	    if (m >= last)
		continue;
	    rdata[k + l].value = data[m - first];
	}
    }
    MPI_Type_free(&type);
#endif /* USE_MPI */

  skip_comm:
    /* copy data to DOFs */
    for (l = 0; l < ndof; l++) {
	u = dofs[l];
	m = DofDataCount(u);
	phgDofClearCache(NULL, u, NULL, NULL, FALSE);
	for (i = 0; i < m; i++) {
	    j = phgMapD2L(map, l, i);
	    if (j < 0)
		u->data[i] = 0.0;
	    else if ((j = phgMapL2V(map, j)) < nlocal)
		u->data[i] = data[j];
	    else
		u->data[i] = rdata[k + j - nlocal].value;
	}
    }

    phgFree(offp_counts);
    phgFree(rdata);

    return 0;
  }
#endif	/*==================================================================*/
	
#if 0
  {
    int id;
    VEC *vec = phgMapCreateVec(map, 0);
    DOF ***pdofs = phgAlloc(ndof * sizeof(*pdofs));

    check_halo(map, __func__);
    for (id = 0; id < ndof; id++)
	pdofs[id] = &dofs[id];
    vec->nvec = 1;
    vec->data = data;
    phgMapVecToDofArraysN(map, vec, FALSE, ndof, pdofs);
    vec->nvec = 0;
    vec->data = NULL;
    phgVecDestroy(&vec);
    phgFree(pdofs);
  }
#else	/* 0|1 */
  {
    int id;
    INT i, k, m;
    FLOAT *offp_data;
    DOF *u;

    check_halo(map, __func__);
    assert(map->dofs != NULL);
    assert(ndof == map->ndof);

    offp_data = phgAlloc((map->localsize - map->nlocal) * sizeof(*data));

#if USE_MPI
    /* fetch remote data */
    if (map->nprocs > 1) {
	if (map->cinfo == NULL)
	    map->cinfo = phgMapCreateCommInfo(
				map->comm,
				map->nprocs,
				map->rank,
				map->nlocal,
				map->localsize,
				map->O2Gmap,
				map->ordering,
				map->partition);
	phgMapScatter(map->cinfo, 1, data, offp_data);
    }
#endif	/* USE_MPI */

    for (id = 0; id < ndof; id++) {
	u = dofs[id];
	m = DofDataCount(u);
	phgDofClearCache(NULL, u, NULL, NULL, FALSE);
	for (i = 0; i < m; i++) {
	    k = phgMapD2L(map, id, i);
	    if (k < 0)
		u->data[i] = 0.0;
	    else if ((k = phgMapL2V(map, k)) < map->nlocal)
		u->data[i] = data[k];
	    else
		u->data[i] = offp_data[k - map->nlocal];
	}
    }

    phgFree(offp_data);
  }
#endif	/* 0|1 */

    return 0;
}

int
phgMapDofToLocalData(MAP *map, int ndof, DOF **dofs, FLOAT *data)
/* This function copies the data in the DOFs to local vector 
 * array (for e.g., setting up an initial solution of a solver).
 * 'data' is the pointer to the local data (size = map->nlocal)
 */
{
    DOF *u;
    int k;
    INT i, j, m;
    INT nlocal = map->nlocal;

    check_halo(map, __func__);

    for (k = 0; k < ndof; k++) {
	u = dofs[k];
	m = DofDataCount(u);
	for (i = 0; i < m; i++) {
	    if ((j = phgMapD2L(map, k, i)) >= 0 &&
		(j = phgMapL2V(map, j)) < nlocal && j >= 0)
		data[j] = u->data[i];
	}
    }

    return 0;
}

VEC *
phgMapDofArraysToVecN(MAP *map, int nvec, BOOLEAN remove_bdry,
			 VEC **vecptr, int ndof, DOF ***dofs)
/* copies data from an array of DOFs to a vector */
{
    VEC *vec = (vecptr == NULL ? NULL : *vecptr);
    int id, iv;
    DOF *u, *u0;
    INT i, k, m, *table;
    FLOAT *data;

    assert(map->dofs != NULL);
    assert(ndof == map->ndof);
    MagicCheck(VEC, vec);

    check_halo(map, __func__);

    if (vec != NULL) {
	phgVecDestroy(&vec);
	if (vecptr != NULL)
	    *vecptr = NULL;
    }

    if (nvec <= 0)
	return vec;

    /* expand bdry[] for faster checking */
    if (!remove_bdry) {
	/* no boundary entries */
	vec = phgMapCreateVec(map, nvec);
	table = phgAlloc(map->localsize * sizeof(*table));
	for (i = 0; i < map->localsize; i++)
	    table[i] = i;
    }
    else {
	MAP *map0 = phgMapRemoveBoundaryEntries_(map, FALSE);
	vec = phgMapCreateVec(map0, nvec);
	phgMapDestroy(&map0);
	table = phgCalloc(map->localsize, sizeof(*table));
	for (i = 0; i < map->bdry_localsize; i++)
	    table[map->bdry[i]] = 1;
	for (i = 0, k = 0; i < map->localsize; i++) {
	    if (table[i] == 1)
		table[i] = -1;
	    else
		table[i] = k++;
	}
    }

    if (vecptr != NULL)
	*vecptr = vec;

    for (id = 0; id < ndof; id++, dofs++) {
	u0 = (*dofs)[0];
	m = DofDataCount(u0);
	for (i = 0; i < m; i++) {
	    if ((k = phgMapD2L(map, id, i)) < 0)
		continue;
	    k = table[phgMapL2V(map, k)];
	    if (k == -1 || k >= vec->map->nlocal)
		continue;
	    data = vec->data + k;
	    for (iv = 0; iv < nvec; iv++) {
		u = (*dofs)[iv];
		assert(u->type == u0->type && u->dim == u0->dim);
		*data = u->data[i];
		data += vec->map->nlocal;
	    }
	}
    }
    phgFree(table);

    return vec;
}

VEC *
phgMapDofArraysToVec(MAP *map, int nvec, BOOLEAN remove_bdry,
			VEC **vecptr, DOF **u, ...)
/* copies data from an array of DOFs to a vector */
{
    int ndof;
    VEC *vec;
    DOF ***dofs;

    GetVariableArgs(ndof, dofs, u, DOF **);
    vec = phgMapDofArraysToVecN(map, nvec, remove_bdry, vecptr, ndof, dofs);
    phgFree(dofs);
    return vec;
}

void
phgMapVecToDofArraysN(MAP *map, VEC *vec, BOOLEAN remove_bdry,
			 int ndof, DOF ***dofs)
/* copies data from a vector to an array of DOFs */
{
    int id, iv, nvec;
    INT i, j, k, m, stride;
    FLOAT *data, *offp_data = NULL;
    static FLOAT zero = 0.0;
    DOF *u0, *u;
    INT *table;

    assert(map->dofs != NULL);
    assert(ndof == map->ndof);
    assert(vec != NULL);
    MagicCheck(VEC, vec);

    check_halo(map, __func__);

    if (vec->nvec <= 0)
	return;

    if (!vec->assembled)
	phgVecAssemble(vec);

    nvec = vec->nvec;

    /* expand bdry[] for faster checking */
    if (!remove_bdry || vec->map->nglobal == map->nglobal) {
	/* no boundary entries */
	assert(vec->map->nglobal   == map->nglobal &&
	       vec->map->localsize == map->localsize);
	table = phgAlloc(map->localsize * sizeof(*table));
	for (i = 0; i < map->localsize; i++)
	    table[i] = i;
    }
    else {
	assert(vec->map->nlocal    == map->nlocal    - map->bdry_nlocal &&
	       vec->map->nglobal   == map->nglobal   - map->bdry_nglobal);
	table = phgCalloc(map->localsize, sizeof(*table));
	for (i = 0; i < map->bdry_localsize; i++)
	    table[map->bdry[i]] = 1;
	for (i = 0, k = 0; i < map->localsize; i++) {
	    if (table[i] == 1)
		table[i] = -1;
	    else
		table[i] = k++;
	}
    }

#if USE_MPI
    /* fetch remote data */
    if (map->nprocs > 1 && map->cinfo == NULL)
	map->cinfo = phgMapCreateCommInfo(
				map->comm,
				map->nprocs,
				map->rank,
				map->nlocal,
				map->localsize,
				map->O2Gmap,
				map->ordering,
				map->partition);

    if (map->cinfo != NULL && map->cinfo->comm != MPI_COMM_NULL) {
	/* collect off-process data */
	int ii;
	FLOAT *sbuff, *rbuff, *p, *q;
	MPI_Datatype type;
	COMM_INFO *cinfo = map->cinfo;

	if (nvec > 1) {
	    MPI_Type_contiguous(nvec, PHG_MPI_FLOAT, &type);
	    MPI_Type_commit(&type);
	}
	else {
	    type = PHG_MPI_FLOAT;
	}

	sbuff = phgAlloc(nvec * (cinfo->ssize + cinfo->rsize) * sizeof(*sbuff));
	rbuff = sbuff + nvec * cinfo->ssize;

	/* copy data to send from data to sbuff */
	p = sbuff;
	for (i = 0; i < cinfo->ssize; i++) {
	    k = table[cinfo->sind[i]];
	    assert(k < vec->map->nlocal);
	    if (k == -1) {
		q = &zero;
		stride = 0;
	    }
	    else {
		q = vec->data + k;
		stride = vec->map->nlocal;
	    }
	    for (ii = 0; ii < nvec; ii++) {
		*(p++) = *q;
		q += stride;
	    }
	}

	MPI_Alltoallv(sbuff, cinfo->scnts, cinfo->sdsps, type,
		      rbuff, cinfo->rcnts, cinfo->rdsps, type, cinfo->comm);

	q = rbuff;
	offp_data = phgAlloc(nvec * cinfo->rsize * sizeof(*offp_data));
	for (i = 0; i < cinfo->rsize; i++) {
	    p = offp_data + cinfo->rind[i];
	    for (ii = 0; ii < nvec; ii++) {
		*p = *(q++);
		p += cinfo->rsize;
	    }
	}

	phgFree(sbuff);
	if (nvec > 1)
	    MPI_Type_free(&type);
    }
#endif	/* USE_MPI */

    for (id = 0; id < ndof; id++, dofs++) {
	for (iv = 0; iv < nvec; iv++)
	    phgDofClearCache(NULL, (*dofs)[iv], NULL, NULL, FALSE);
	u0 = (*dofs)[0];
	m = DofDataCount(u0);
	for (i = 0; i < m; i++) {
	    k = phgMapD2L(map, id, i);
	    if (k < 0 || (k = table[j = phgMapL2V(map, k)]) == -1) {
		data = &zero;
		stride = 0;
	    }
	    else if (k < vec->map->nlocal) {
		data = vec->data + k;
		stride = vec->map->nlocal;
	    }
	    else {
		assert(j >= map->nlocal);
		data = offp_data + j - map->nlocal;
		stride = map->localsize - map->nlocal;
	    }
	    for (iv = 0; iv < nvec; iv++) {
		u = (*dofs)[iv];
		assert(u->type == u0->type && u->dim == u0->dim);
		u->data[i] = *data;
		data += stride;
	    }
	}
    }

    phgFree(offp_data);
    phgFree(table);
}

void
phgMapVecToDofArrays(MAP *map, VEC *vec, BOOLEAN remove_bdry, DOF **u, ...)
/* copies data from a vector to an array of DOFs */
{
    int ndof;
    DOF ***dofs;

    GetVariableArgs(ndof, dofs, u, DOF **);
    phgMapVecToDofArraysN(map, vec, remove_bdry, ndof, dofs);
    phgFree(dofs);
}

VEC *
phgMapCreateVec(MAP *map, int nvec)
{
    VEC *v = phgCalloc(1, sizeof(*v));

    map->refcount++;

    v->magic = MAGIC_VEC;
    v->map = map;
    v->mat = NULL;

#if USE_MPI
    v->O2Gmap = NULL;
    v->ordering = NULL;
    v->localsize = 0;
    v->alloc = 0;
#endif	/* USE_MPI */

    v->assembled = TRUE;
    v->nvec = nvec;
    v->data = phgCalloc(nvec * map->nlocal, sizeof(*(v->data)));
    v->offp_data = phgCalloc(nvec * (map->localsize - map->nlocal),
						sizeof(*(v->data)));

    return v;
}

MAT *
phgMapCreateMat(MAP *rmap, MAP *cmap)
{
    MAT *Mp;

    FunctionEntry;

    Mp = phgCalloc(1, sizeof(*Mp));
    Mp->magic = MAGIC_MAT;
    Mp->rmap = rmap;
    Mp->cmap = cmap;
    Mp->type = PHG_UNPACKED;

    /* check for special matrix */
    if (Mp->rmap->nprocs == 1 && Mp->cmap->nprocs > 1) {
	phgPrintf("\n*\n* Warning: Matrices with serial rows and distributed "
		  "columns are inefficient.\n");
	phgPrintf("* Warning: It is better to use distributed rows and serial "
		  "columns instead.\n*\n");
    }

#if USE_MPI
    Mp->O2Gmap = NULL;
    Mp->ordering = NULL;
    Mp->cinfo = NULL;
#endif	/* USE_MPI */
    Mp->localsize = -1;	/* to make convert_row_index() initialize Mp->O2Gmap,
			 * etc., when needed. */

    Mp->rows = phgCalloc(rmap->localsize, sizeof(*Mp->rows));
    Mp->nnz_d = Mp->nnz_o = 0;

    Mp->diag = NULL;
    Mp->packed_data = NULL;
    Mp->packed_cols = NULL;
    Mp->packed_ind = NULL;

    Mp->mv_func = NULL;
    Mp->mv_data = NULL;

    Mp->assembled = FALSE;
    Mp->mode = PHG_UNDEFINED;

    rmap->refcount++;
    cmap->refcount++;

#if USE_OMP
    Mp->nlock = -1;
    Mp->locks = NULL;
    omp_init_lock(&Mp->lock);
#endif	/* USE_OMP */

    Return Mp;

    return Mp;	/* to make MIPS cc happy */
}

MAT *
phgMapCreateMatrixFreeMat(MAP *rmap, MAP *cmap, MV_FUNC mv_func,
			  void *mv_data, ...)
{
    extern void _phg_matrix_free_setup(MAT *, MV_FUNC, int, void **);
    MAT *mat = phgMapCreateMat(rmap, cmap);
    int ndata;
    void **data;

    if (mat == NULL)
	return NULL;

    GetVariableArgs(ndata, data, mv_data, void *);
    _phg_matrix_free_setup(mat, mv_func, ndata, data);
    phgFree(data);

    return mat;
}

BOOLEAN
phgMapCompare(MAP *map1, MAP *map2)
/* compares two maps, returns TRUE if they are compatible, FALSE otherwise */
{
    if (map1 == map2)
	return TRUE;

    if (map1->nprocs != map2->nprocs || map1->rank != map2->rank ||
	map1->nglobal != map2->nglobal || map1->nlocal != map2->nlocal)
	return FALSE;

#if USE_MPI
    if (map1->nprocs > 1) {
	int res;

	MPI_Comm_compare(map1->comm, map2->comm, &res);
	if (res != MPI_IDENT && res != MPI_CONGRUENT)
	    return FALSE;

	if (memcmp(map1->partition, map2->partition,
			(map1->nprocs + 1) * sizeof(*map1->partition)))
	    return FALSE;
    }
#endif	/* USE_MPI */

    return TRUE;
}

/*----------------- functions for exchanging neighbour data ---------------*/

/* Note:
 * 	These functions can be regarded as an extension to the functions
 * 	phgDof*NeighbourData (defined in dof.c). They can be used to fetch
 * 	whole data related to a given MAP from the neighbouring elements.
 *
 * The algorithm:
 *
 * 	1. Collect all remote neighbours, MAP indices, and optionally DOF data.
 *
 * 	2. Construct a temporary GRID for the remote elements, in which the
 * 	   ordering of the vertices preserves the order of global vertices
 * 	   of the original GRID
 *
 * 	3. Build a temporary MAP on the temporary GRID
 *
 * 	4. Set map->L2Vmap[] to the global indices of the original map,
 * 	   and set map->bogus_Lmap = TRUE to prevent users from accidentally
 * 	   using the local indices of the temporary MAP.
 *
 *	5. The function phgMapNeighbourMap returns the temporary MAP and
 *	   elements in the temporary GRID for remote neighbours, and
 *	   the original MAP and elements in the original GRID for local
 *	   neighbours.
 */

#if USE_MPI
static INT *comp_base;

static int
comp_index(const void *p0, const void *p1)
/* compares two vertices. p0 and p1 are pointers to 'comp_base' */
{
    INT i = comp_base[*(INT *)p0] - comp_base[*(INT *)p1];
    return (i < 0) ? -1 : ((i == 0) ? 0 : 1);
}

static int
comp_elem(const void *p0, const void *p1)
/* compares vertices of two elements. p0 and p1 are pointers to 'comp_base' */
{
    return memcmp(comp_base + *(const INT *)p0 * (Dim + 1),
		  comp_base + *(const INT *)p1 * (Dim + 1),
		  (Dim + 1) * sizeof(INT));
}
#endif	/* USE_MPI */

NEIGHBOUR_MAP *
phgMapInitNeighbourMap(MAP *map, GTYPE kind, BOOLEAN dof_data)
/* sets up neighbour data for the given MAP.
 *  'map' is the MAP pointer.
 *  'kind' specifies the kind of neighbours (VERTEX, EDGE or FACE).
 *  'dof_data' specifies whether DOF data are needed.
 */
{
    NEIGHBOUR_MAP *nm = phgCalloc(1, sizeof(*nm));

#if !USE_MPI
    assert(kind == FACE);		/* VERTEX and EDGE unimplemented */
    assert(map->dofs != NULL);
    assert(dof_data == FALSE);		/* unimplemented */

    nm->g_local = map->dofs[0]->g;
    nm->map_local = map;
    nm->kind = kind;
#else	/* !USE_MPI */
    GRID *g, *g1;
    MAP *map1;
    ELEMENT *e, **pp;
    RNEIGHBOUR *rn;
    DOF *dof, **dofs;
    void *sbuf, *rbuf;
    int ii, jj, nbas, fcount, icount, ioffset, pcount, poffset, size;
    int *scnts, *sdsps, *rcnts, *rdsps;
    INT scount, rcount, i, j, *ip;
    MPI_Datatype type;
    /* arrays for phgImportParallelGrid */
    INT nvert, nelem, *elems, *index, *e2rbuf;
    INT *L2Gmap_vert = NULL;	/* for debugging only */
    FLOAT *coord, *fp;
    SIMPLE_TET *tet;

    nm->kind = kind;

    assert(kind == FACE);		/* VERTEX and EDGE unimplemented */
    assert(map->dofs != NULL);
    assert(dof_data == FALSE);		/* unimplemented */

    g = map->dofs[0]->g;

    if (g->nprocs == 1) {
	nm->g_local = map->dofs[0]->g;
	nm->map_local = map;
	nm->kind = kind;
	return nm;
    }

    scnts = rcnts = g->neighbours.counts;
    sdsps = rdsps = g->neighbours.displs;

    /* count number of bases in an element */
    nbas = 0;
    for (jj = 0; jj < map->ndof; jj++) {
	dof = map->dofs[jj];
	assert(g == dof->g);
	if (DofIsHP(dof))
	    phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__);
	nbas += dof->type->nbas;
    }

    /* data to send:
     *	- g_index of the opposite vertex:	1 INT
     *	- e_index of the remote op_vertex:	1 INT
     *	- global MAP indices of the bases:	nbas INTs
     *	- address of the remote neighbour:	1 pointer
     *	- coordinates of the opposite vertex:	Dim FLOATs
     *	- DOF data: TODO */
    icount = nbas + 2;
    pcount = 1;
    fcount = Dim;

#define ALIGN(loc,align) (((loc + align - 1) / align) * align)
    poffset = ALIGN(fcount * sizeof(*fp), sizeof(*pp));
    ioffset = ALIGN(pcount * sizeof(*pp) + poffset, sizeof(*ip));
    size = ioffset + icount * sizeof(*ip);

    scount = sdsps[g->nprocs-1] + scnts[g->nprocs-1];
    rcount = rdsps[g->nprocs-1] + rcnts[g->nprocs-1];
    assert(scount == rcount);	/* only TRUE if no hp-DOFs */

    sbuf = phgAlloc(size * (size_t)scount);
    rbuf = phgAlloc(size * (size_t)rcount);

    /* prepare data to send */
#define OFFSET(pos, offset) (void *)(((BYTE *)(pos)) + (offset))
    for (i = 0; i < g->neighbours.count; i++) {
	fp = OFFSET(sbuf, size * i + 0);
	pp = OFFSET(sbuf, size * i + poffset);
	ip = OFFSET(sbuf, size * i + ioffset);
	rn = g->neighbours.list + i;
	e = rn->local;
	/* g_index of the vertex opposite to the face */
	*(ip++) = GlobalVertex(g, e->verts[rn->vertex]);
	/* e_index of the remote op_vertex */
	*(ip++) = rn->op_vertex;
	/* address (pointer) of the remote neighbour */
	*(pp++) = rn->remote;
	/* coordinates of the vertex opposite to the face */
	for (ii = 0; ii < Dim; ii++)
	    *(fp++) = g->verts[e->verts[rn->vertex]][ii];
	/* Global MAP indices of the bases */
	for (jj = 0; jj < map->ndof; jj++) {
	    dof = map->dofs[jj];
	    if (DofIsHP(dof))
		phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__);
	    for (ii = 0; ii < dof->type->nbas; ii++)
		*(ip++) = phgMapE2G(map, jj, e, ii);
	}
    }

    /* exchange neighbour data */
    MPI_Type_contiguous(size, MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(sbuf, scnts, sdsps, type, rbuf, rcnts, rdsps, type, g->comm);
    MPI_Type_free(&type);
    phgFree(sbuf);
    if (scnts != g->neighbours.counts)
	phgFree(scnts);		/* free dynamically allocated buffer */

    /* create temporary GRID */

    /* collect elements with global vertex numbers */
    elems = phgAlloc(2 * (Dim + 1) * g->neighbours.count * sizeof(*elems));
    index = elems + (Dim + 1) * g->neighbours.count;
    for (i = 0; i < g->neighbours.count; i++) {
	pp = OFFSET(rbuf, size * i + poffset);
	ip = OFFSET(rbuf, size * i + ioffset);
	rn = GetRNeighbour(g, pp[0], ip[1]);
	e = rn->local;
	elems[(Dim + 1) * i + rn->op_vertex] = *ip; /* gindex of op_vertex */
	/* get the 3 vertices on the shared face from the local ELEMENT */
	for (ii = 0; ii < NVertFace; ii++)
	    elems[(Dim + 1) * i + rn->rface[ii]] =
			GlobalVertex(g, e->verts[GetFaceVertex(rn->vertex,ii)]);
	/* initialize index[] for sorting the vertices */
	for (ii = 0; ii < Dim + 1; ii++)
	    index[(Dim + 1) * i + ii] = (Dim + 1) * i + ii;
    }

    /* pack global vertex numbers */
    comp_base = elems;
    if (g->neighbours.count > 0)
	qsort(index, (Dim+1) * g->neighbours.count, sizeof(*index), comp_index);
    nvert = i = 0;
    while (i < (Dim + 1) * g->neighbours.count) {
	j = elems[index[i]];
	elems[index[i]] = nvert;
	while (++i < (Dim + 1) * g->neighbours.count && elems[index[i]] == j)
	    elems[index[i]] = nvert;
#if DEBUG
	/* temporaryly save L2Gmap_vert[] entries in index[] */
	index[nvert] = j;
#endif	/* DEBUG */
	nvert++;
    }

#if DEBUG
    L2Gmap_vert = phgAlloc(nvert * sizeof(*L2Gmap_vert));
    for (i = 0; i < nvert; i++)
	L2Gmap_vert[i] = index[i];
#endif	/* DEBUG */

    /* collect coordinates */
    coord = phgAlloc(Dim * sizeof(*coord) * nvert);
    for (i = 0; i < g->neighbours.count; i++) {
	fp = OFFSET(rbuf, size * i + 0);
	pp = OFFSET(rbuf, size * i + poffset);
	ip = OFFSET(rbuf, size * i + ioffset);
	rn = GetRNeighbour(g, pp[0], ip[1]);
	e = rn->local;
	/* the opposite vertex */
	j = elems[(Dim + 1) * i + rn->op_vertex];
	coord[j * Dim + 0] = fp[0];
	coord[j * Dim + 1] = fp[1];
	coord[j * Dim + 2] = fp[2];
	/* get the 3 vertices on the shared face from the local ELEMENT */
	for (ii = 0; ii < NVertFace; ii++) {
	    e = rn->local;
	    COORD *c = g->verts + e->verts[GetFaceVertex(rn->vertex,ii)];
	    j = elems[(Dim + 1) * i + rn->rface[ii]];
	    coord[j * Dim + 0] = (*c)[0];
	    coord[j * Dim + 1] = (*c)[1];
	    coord[j * Dim + 2] = (*c)[2];
	}
    }

    /* pack elements, removing duplicate ones */
    for (i = 0; i < g->neighbours.count; i++)
	index[i] = i;
    if (g->neighbours.count > 0)
	qsort(index, g->neighbours.count, sizeof(*index), comp_elem);
    nelem = i = 0;
    while (i < g->neighbours.count) {
	j = i;
	while (++i < g->neighbours.count && !comp_elem(index + j, index + i))
	    ;
	nelem++;
    }
    tet = phgAlloc(nelem * sizeof(*tet));
    nm->index = phgAlloc((g->neighbours.count + 1) * sizeof(*nm->index));
    nm->neighbours = phgAlloc(g->neighbours.count * sizeof(*nm->neighbours));
    /* Note: e2rbuf[] stores the map g1->roots[] ==> rbuf[] */
    e2rbuf = phgAlloc(nelem * sizeof(*e2rbuf));
    nelem = i = 0;
    while (i < g->neighbours.count) {
	pp = OFFSET(rbuf, size * index[i] + poffset);
	ip = OFFSET(rbuf, size * index[i] + ioffset);
	rn = GetRNeighbour(g, pp[0], ip[1]);
	/* Note: temporarily use nm->index[] to save the map:
	 *	g->neighbours.list[] ==> g1->roots[] */
	nm->index[rn - g->neighbours.list] = nelem;
	j = i;
	memcpy(tet[nelem].verts,
		elems + (Dim + 1) * index[j], sizeof(tet[0].verts));
	for (ii = 0; ii < Dim + 1; ii++)
	    tet[nelem].bound_type[ii] = NEUMANN; /* may save some processing */
	e2rbuf[nelem] = index[i];	/* g1->roots[] ==> rbuf[]i mapping */
	while (++i < g->neighbours.count && !comp_elem(index + j, index + i)) {
	    pp = OFFSET(rbuf, size * index[i] + poffset);
	    ip = OFFSET(rbuf, size * index[i] + ioffset);
	    rn = GetRNeighbour(g, pp[0], ip[1]);
	    nm->index[rn - g->neighbours.list] = nelem;
	}
	nelem++;
    }

    phgFree(elems);

    g1 = phgImportParallelGrid(nvert, nelem, nvert, nelem,
				L2Gmap_vert, NULL, coord, tet, MPI_COMM_SELF);
    phgFree(L2Gmap_vert);

    phgFree(coord);
    phgFree(tet);

    for (i = 0; i < g->neighbours.count; i++) {
	/* this assumes that phgImportParallelGrid preserves the order
	 * of the input elements */
	nm->neighbours[i] = g1->roots + nm->index[i];
	nm->index[i] = i;
    }
    nm->index[g->neighbours.count] = g->neighbours.count;

    nm->g_remote = g1; 
    nm->g_local = g;

#if 0
    /* dump temporary GRID for visualization */
    {
	char s[1024];
	sprintf(s, "tmp-%d.vtk", g->rank);
	phgInfo(-1, "Writing \"%s\".\n", phgExportVTK(g1, s, NULL));
    }
#endif

    /* Create DOFs with g1 */
    dofs = phgAlloc(map->ndof * sizeof(*dofs));
    for (ii = 0; ii < map->ndof; ii++) {
	dofs[ii] = phgDofNew(g1,
			map->dofs[ii]->type,
			map->dofs[ii]->dim,
			map->dofs[ii]->name,
			DofNoData);
    }

    /* Create the temporary MAP */
    nm->map_local = map;
    nm->map_remote = map1 = phgMapCreateN(map->ndof, dofs);
    phgFree(dofs);
    /* update 'map1'->L2Vmap[] to match the global indices of 'map' */
    assert(map1->L2Vmap == NULL);
    map1->L2Vmap_size = map1->localsize;
    map1->L2Vmap = phgAlloc(map1->L2Vmap_size * sizeof(map1->L2Vmap));
    /* set all entries to -1 to facilitate debugging */
    for (i = 0; i < map1->L2Vmap_size; i++)
	map1->L2Vmap[i] = -1;
    for (i = 0; i < g->neighbours.count; i++) {
	e = nm->neighbours[i];
	ip = OFFSET(rbuf, size * e2rbuf[e - g1->roots] + ioffset);
	ip += 2;	/* skip two entries: 'gindex' and 'op_vertex' */
	for (jj = 0; jj < map1->ndof; jj++) {
	    dof = map1->dofs[jj];
	    if (DofIsHP(dof))
		phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__);
	    for (ii = 0; ii < dof->type->nbas; ii++) {
		map1->L2Vmap[phgMapE2L_(map1, jj, e, ii, TRUE, TRUE)] = *(ip++);
	    }
	}
    }
    phgFree(e2rbuf);
    /* check */
#if DEBUG
    for (i = 0; i < map1->L2Vmap_size; i++) {
	if (map1->L2Vmap[i] < 0 || map1->L2Vmap[i] >= map->nglobal)
	    phgError(1, "invalid entry: map1->L2Vmap[%d] = %d\n",
							i, map1->L2Vmap[i]);
    }
#endif	/* DEBUG */
    map1->bogus_Lmap = TRUE;

    phgFree(rbuf);
#endif	/* !USE_MPI */

    return nm;
}

void
phgMapReleaseNeighbourMap(NEIGHBOUR_MAP **nm_ptr)
{
    NEIGHBOUR_MAP *nm;

    assert(nm_ptr != NULL);
    if ((nm = *nm_ptr) == NULL)
	return;

    if (nm->map_remote != NULL) {
	int ii;
	assert(nm->map_remote->dofs != NULL);
	for (ii = 0; ii < nm->map_remote->ndof; ii++)
	    phgDofFree(nm->map_remote->dofs + ii);
	phgMapDestroy(&nm->map_remote);
    }
    phgFreeGrid(&nm->g_remote);
    phgFree(nm->index);
    phgFree(nm->neighbours);

    phgFree(nm);

    *nm_ptr = NULL;
}

NEIGHBOUR_INFO *
phgMapNeighbourMap(NEIGHBOUR_MAP *nm, ELEMENT *e, int index)
/* 
 * Returns the neighbours of the element 'e' as a {NULL,*,*} terminated list of
 * NEIGHBOUR_INFOs. 'index' specifies the {vertex|edge|face} number.
 */
{
    /* size of ni[] >= n + 1, where n is the max. number of neighbours */
    static NEIGHBOUR_INFO ni[2] = {{NULL, 0, 0}, {NULL, 0, 0}};
#if USE_OMP
#pragma omp threadprivate(ni)
#endif	/* USE_OMP */

    assert(nm->kind == FACE);

#if USE_MPI
    if ((e->bound_type[index] & REMOTE)) {
	/* remote neighbour */
	ni[0].e = nm->neighbours[nm->index[(size_t)(e->neighbours[index])]];
	ni[0].is_remote = TRUE;
	ni[0].index = nm->g_local->neighbours.list
			[(size_t)(e->neighbours[index])].op_vertex;
    }
    else
#endif	/* USE_MPI */
    if ((e->bound_type[index] & INTERIOR)) {
	/* local neighbour */
	ni[0].e = e->neighbours[index];
	ni[0].is_remote = FALSE;
	ni[0].index = phgOppositeFace(nm->g_local, e, index, ni[0].e);
    }
    else {
	/* boundary face, no neighbour */
	ni[0].e = NULL;
    }

    return ni;
}

/*-------------------------------------------------------------------------*/
