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

/* $Id: halo.c,v 1.65 2022/08/13 01:31:42 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>	/* INT_MAX */

ELEMENT *
phgGetNeighbour_(GRID *g, ELEMENT *e, int face, const char *file, int line)
/* returns the neighbour element for interior face and NULL for boundary face.
 * Note: if the neighbour is in another submesh, then g->halo must be set up
 * before calling this function */
{
    if (g->halo != NULL && e->index >= g->halo->nelem_bak)
	phgError(1, "%s:%d: can't access neighbours of halo elements.\n",
			file, line);

    if (!HasRemoteNeighbour(e, face))
	return GetNeighbour(e, face);

    if (g->halo == NULL)
	phgError(1, "%s:%d: cannot access remote neighbours, please call "
		"phgSetupHalo before using this function.\n", file, line);

    return g->halo->elist + g->halo->map[(size_t)(e->neighbours[face])];
}

#if USE_MPI
/* static variables for searching sorted arrays comp_base[order[...]] */
static INT *comp_base;

static int
comp_ordering(const void *p0, const void *p1)
{
    INT i0 = *(const INT *)p0, i1 = *(const INT *)p1, k;

    if ((k = comp_base[i0] - comp_base[i1]) != 0)
	return k < 0 ? -1 : 1;
    return i0 == i1 ? 0 : (i0 < i1 ? -1 : 1);
}

/* static variables for sorting new items */
static struct {
    INT gindex, oindex, pindex;	/* global/owner/periodic index of the item */
    INT rbuf_index;		/* rbuf index: rbuf[rbuf_index].elem */
    int item_index, orank;	/* item index and owner rank, e.g.,
				 * rbuf[rbuf_index].elem.verts[item_index] */
    BTYPE type;			/* boundary type */
} *items, *pi;

/* struct for storing data from the remote neighbour
 */
static struct {
    ELEMENT	elem;		/* Element data, whose index, verts[], edges[]
				 * and faces[] members contain global indices,
				 * neighbours["ii"] stores the pointer of the
				 * remote neighbour in the peer process on face
				 * "ii" (NULL if no remote neighbour in the peer
				 * process), mark stores vertex and op_vertex
				 * of the neighbouring elements */
    COORD	coord[4];	/* coordinates of the vertices */
    INT		pindexes[4];	/* g->period->L2Gmap_vert */
    INT		indexes[15];	/* owner indexes: V+E+F+E = 4+6+4+1 */
    int		ranks[15];	/* owner ranks:	V+E+F+E = 4+6+4+1 */
    BTYPE	types[15];	/* types: V+E+F+E = 4+6+4+1 */
} *sbuf, *rbuf;

static int
comp_item(const void *p0, const void *p1)
{
    INT i0 = items[*(const INT *)p0].gindex,
	i1 = items[*(const INT *)p1].gindex;
    return i0 == i1 ? 0 : (i0 < i1 ? -1 : 1);
}

static int
comp_elem(const void *p0, const void *p1)
{
    INT i0 = rbuf[*(const INT *)p0].elem.index,
	i1 = rbuf[*(const INT *)p1].elem.index;
    return i0 == i1 ? 0 : (i0 < i1 ? -1 : 1);
}

static int
comp_pointer(const void *p0, const void *p1)
{
    const void *q0 = *(const void **)p0, *q1 = *(const void **)p1;
    return q0 == q1 ? 0 : (q0 < q1 ? -1 : 1);
}

void
phgDofUpdateGhost_(DOF *dof, BOOLEAN flag)
/* adjust size of dof->data, and update ghost data if "flag" == TRUE */
{
    GRID *g = dof->g;
    MAP *map;
    FLOAT *data;

    if (SpecialDofType(dof->type) || dof->userfunc == DofNoData)
	return;

    phgInfo(1, "%s: update DOF %s\n", __func__, dof->name);

    /* update halo data */

    if (dof == g->geom) {
	/* Note: use FOR_UNOWNED to only compute halo data */
	phgDofSetDataByFunction__(dof, NULL, NULL, NULL, ~0, FOR_UNOWNED);
	return;
    }

    if (!flag)
	return;

    if (dof->userfunc == DofLambdaFunction ||
	!SpecialDofUserFunction(dof->userfunc)) {
	phgDofSetDataByFunction__(dof,
	    dof->userfunc, dof->userfunc_lambda, NULL, ~0, FOR_UNOWNED);
	return;
    }

    /* fetch data in the halo from remote processes */
    map = phgMapCreate(dof, NULL);
    data = phgAlloc(DofDataCount(dof) * sizeof(*data));
    phgMapDofToLocalData(map, 1, &dof, data);
    phgMapLocalDataToDof(map, 1, &dof, data);
    phgFree(data);
    phgMapDestroy(&map);
}
#else	/* USE_MPI */
void
phgDofUpdateGhost_(DOF *dof, BOOLEAN flag)
{
    /* do nothing */
    Unused(dof);
    Unused(flag);
}
#endif	/* USE_MPI */

#define REALLOC(p, n, m) p = phgRealloc_(p, (n) * sizeof(*p), (m) * sizeof(*p))

void
phgSetupHalo_(GRID *g, int halo_type, BOOLEAN update_dofs)
/* If 'update_dofs' != FALSE, then previously defined DOF objects will be
 * updated with correct data for the halo elements unless userfunc==DofNoAction.
 *
 * The argument "halo_type" specifies the kind and size (width) of the halo,
 * see the macros HaloType, GetHaloSize, GetHaloKind, etc, in halo.h.
 */
{
#if USE_MPI
    INT i, j, k, n, m, *ordering, scount, rcount;
    int ii, jj, rank;
    ELEMENT *e;
    MPI_Datatype type;
    DOF **pdof;
    HP_TYPE **php;
    int halo_kind, halo_size;

    if (g->halo != NULL) {
	if (g->halo->type <= halo_type)
	    return;
	phgRemoveHalo(g);
	assert(g->halo == NULL);
	g->serial_no -= 1000 * 1000;
    }
    g->serial_no += 1000 * 1000;

    if (halo_type < 0)
	return;

    halo_kind = GetHaloKind(halo_type);
    halo_size = GetHaloSize(halo_type);
    assert(halo_type >= HALO_ALL && halo_type <= HALO_FACE);

    if (halo_size > 1)
	phgError(1, "%s (%s:%d): unimplemented feature.\n",
				__func__, __FILE__, __LINE__);

    g->halo = phgCalloc(1, sizeof(*g->halo));
    g->halo->type = halo_type;
    g->halo->nvert_bak = g->nvert;
    g->halo->nedge_bak = g->nedge;
    g->halo->nface_bak = g->nface;
    g->halo->nelem_bak = g->nelem;

    if (g->nprocs == 1)
	return;

    assert(halo_type != HALO_ALL);	/* TODO */

    /* invalidate face->element map */
    phgFree(g->Fmap_elems);
    g->Fmap_elems = NULL;
    phgFree(g->Fmap_faces);
    g->Fmap_faces = NULL;

    if (g->owner_index_face == NULL) {
	/* set up g->owner_{index,rank}_face */
	assert(g->owner_rank_face == NULL);
	g->owner_index_face = phgAlloc(g->nface * sizeof(*g->owner_index_face));
	g->owner_rank_face = phgAlloc(g->nface * sizeof(*g->owner_rank_face));
	for (i = 0; i < g->nface; i++) {
	    g->owner_index_face[i] = -1;
	    g->owner_rank_face[i] = -1;
	}
	for (i = 0; i < g->neighbours.count; i++) {
	    RNEIGHBOUR *rn = g->neighbours.list + i;
	    j = rn->local->faces[rn->vertex];
	    if ((g->types_face[j] & OWNER))
		continue;
	    g->owner_index_face[j] = rn->peer_face;
	    g->owner_rank_face[j] = rn->rank;
	}
    }

    if (g->owner_index_elem == NULL) {
	/* set up g->owner_{index,rank}_elem */
	assert(g->owner_rank_elem == NULL);
	g->owner_index_elem = phgAlloc(g->nface * sizeof(*g->owner_index_elem));
	g->owner_rank_elem = phgAlloc(g->nface * sizeof(*g->owner_rank_elem));
	for (i = 0; i < g->nelem; i++) {
	    if ((g->types_elem[i] & OWNER)) {
		g->owner_index_elem[i] = g->elems[i]->index;
		g->owner_rank_elem[i] = g->rank;
	    }
	    else {
		g->owner_index_elem[i] = -1;
		g->owner_rank_elem[i] = -1;
	    }
	}
    }

    /* Construct list of elements to send to each process in halo->elocal,
     * with the number and location of elements to send/receive to/from
     * proc ii stored in halo->scnt[ii]/rcnt[ii] and halo->sdsp[ii]/rdsp[ii]
     * respectively. */

    g->halo->scnt = phgAlloc(4 * g->nprocs * sizeof(*g->halo->scnt));
    g->halo->rcnt = g->halo->scnt + g->nprocs;
    g->halo->sdsp = g->halo->rcnt + g->nprocs;
    g->halo->rdsp = g->halo->sdsp + g->nprocs;

    if (halo_kind != HALO_FACE) {
	VEF_MAP2 *ve = phgDofSetupVEFMap2_(g, NULL,
				halo_kind == HALO_VERT ? VERT_FLAG : EDGE_FLAG,
				REMOTE, g->types_vert, g->types_edge, NULL,
				FOR_OWNED);
	BTYPE *types;
	ELEMENT ***map;
	int *size, *orank, *isend, *irecv;
	INT n, *send, *recv, *iorder, n0;

	/* The list of elements to send to each process is contructed in the
	 * following steps:
	 *   1. Each process sends the global index of each shared vertex/edge
	 *	to the owner process of the vertex/edge (or better, to
	 *	nvert_global % g->nprocs for load balancing?);
	 *   2. The owner process processes received global indices:
	 *	for each received global index, collects the list of ranks
	 *	having this index and sends the list of ranks back to the
	 *	incoming processes;
	 *   3. Each process can then construct the list of elements to send
	 *	to each other process in elocal[].
	 */

	/* collect list of elements and owner process for each vertex or edge */
	if (halo_kind == HALO_VERT) {
	    n = g->nvert;
	    types = g->types_vert;
	    orank = g->owner_rank_vert;
	    map = ve->Vmap;
	    size = ve->Vsize; 
	}
	else {
	    n = g->nedge;
	    types = g->types_edge;
	    orank = g->owner_rank_edge;
	    map = ve->Emap;
	    size = ve->Esize; 
	}

	for (ii = 0; ii < g->nprocs; ii++)
	    g->halo->scnt[ii] = 0;

	for (i = 0; i < n; i++) {
	    if (!(types[i] & REMOTE))
		continue;
	    assert((size_t)g->halo->scnt[orank[i]] < INT_MAX);
	    g->halo->scnt[orank[i]]++;
	}
	MPI_Alltoall(g->halo->scnt, 1, MPI_INT, g->halo->rcnt, 1, MPI_INT,
		     g->comm);

	scount = rcount = 0;
	for (ii = 0; ii < g->nprocs; ii++) {
	    assert(scount <= INT_MAX && rcount <= INT_MAX);
	    g->halo->sdsp[ii] = scount;
	    g->halo->rdsp[ii] = rcount;
	    scount += g->halo->scnt[ii];
	    rcount += g->halo->rcnt[ii];
	    g->halo->scnt[ii] = 0;	/* reset counter */
	}

	send = phgAlloc(scount * sizeof(*send));
	recv = phgAlloc(rcount * sizeof(*recv));

	/* collect data to send */
	for (i = 0; i < n; i++) {
	    if (!(types[i] & REMOTE))
		continue;
	    ii = orank[i];
	    k = halo_kind == HALO_VERT ? GlobalVertex(g, i) : GlobalEdge(g, i);
	    send[g->halo->sdsp[ii] + (g->halo->scnt[ii]++)] = k;
	}

	MPI_Alltoallv(send, g->halo->scnt, g->halo->sdsp, PHG_MPI_INT,
		      recv, g->halo->rcnt, g->halo->rdsp, PHG_MPI_INT, g->comm);

	/* change entries in send[] to local indices, which will be needed when
	 * processing data received in irecv[] */
	for (ii = 0; ii < g->nprocs; ii++)
	    g->halo->scnt[ii] = 0;
	for (i = 0; i < n; i++) {
	    if (!(types[i] & REMOTE))
		continue;
	    ii = orank[i];
	    send[g->halo->sdsp[ii] + (g->halo->scnt[ii]++)] = i;
	}
	n0 = scount;	/* save number of shared indices */

	/* sort the list of incoming Global indices */
	ordering = phgAlloc(2 * rcount * sizeof(*ordering));
	iorder = ordering + rcount;
	for (i = 0; i < rcount; i++)
	    ordering[i] = i;
	if (rcount > 0) {
	    comp_base = recv;
	    qsort(ordering, rcount, sizeof(*ordering), comp_ordering);
	}
	/* construct the inverse map of ordering[] */
	for (i = 0; i < rcount; i++)
	    iorder[ordering[i]] = i;

	/* construct list of ranks for each received vertex, each list
	 * is preceeded by a count as;
	 *	count, rank_0, ..., rank_{count-1} */

	/* first, count # of ranks for each process */
	for (ii = 0; ii < g->nprocs; ii++)
	    g->halo->scnt[ii] = 0;

	ii = 0;	/* incoming rank */
	for (i = 0; i < rcount; i++) {
	    while (ii < g->nprocs && g->halo->rdsp[ii] + g->halo->rcnt[ii] <= i)
		ii++;
	    assert(ii < g->nprocs);
	    k = recv[i];
	    /* determine the sorted interval [j,m) for the global index 'k' */
	    j = iorder[i];	/* position in the sorted array */
	    for (m = j + 1; m < rcount && recv[ordering[m]] == k; m++);
	    for (; j > 0 && recv[ordering[j - 1]] == k; j--);
	    /* m-j ranks for index k, minus 1 (ii), plus a counter = m-j-1 */
	    assert((size_t)g->halo->scnt[ii] + m - j < INT_MAX);
	    g->halo->scnt[ii] += m - j;
	}

	scount = 0;
	for (ii = 0; ii < g->nprocs; ii++) {
	    assert(scount <= INT_MAX);
	    g->halo->sdsp[ii] = scount;
	    scount += g->halo->scnt[ii];
	    g->halo->scnt[ii] = 0;		/* reset counter */
	}

	/* copy lists to isend[] */
	isend = phgAlloc(scount * sizeof(*isend));
	ii = 0;	/* incoming rank */
	for (i = 0; i < rcount; i++) {
	    while (ii < g->nprocs && g->halo->rdsp[ii] + g->halo->rcnt[ii] <= i)
		ii++;
	    assert(ii < g->nprocs);
	    k = recv[i];
	    /* determine the sorted interval [j,m) for the global index 'k' */
	    j = iorder[i];	/* position in the sorted array */
	    for (m = j + 1; m < rcount && recv[ordering[m]] == k; m++);
	    for (; j > 0 && recv[ordering[j - 1]] == k; j--);
	    isend[g->halo->sdsp[ii] + (g->halo->scnt[ii]++)] = m - j - 1;
	    jj = 0;	/* target rank */
	    for (k = j; k < m; k++) {
		while (jj < g->nprocs &&
		       ordering[k] >= g->halo->rdsp[jj] + g->halo->rcnt[jj])
		    jj++;
		assert(jj < g->nprocs);
		if (jj == ii)
		    continue;
		isend[g->halo->sdsp[ii] + (g->halo->scnt[ii]++)] = jj;
	    }
	}
	phgFree(ordering);
	phgFree(recv);

	/* send back data */
	MPI_Alltoall(g->halo->scnt, 1, MPI_INT, g->halo->rcnt, 1, MPI_INT,
		     g->comm);
	rcount = 0;
	for (ii = 0; ii < g->nprocs; ii++) {
	    g->halo->rdsp[ii] = rcount;
	    rcount += g->halo->rcnt[ii];
	}

	irecv = phgAlloc(rcount * sizeof(*irecv));
	MPI_Alltoallv(isend, g->halo->scnt, g->halo->sdsp, MPI_INT,
		      irecv, g->halo->rcnt, g->halo->rdsp, MPI_INT, g->comm);
	phgFree(isend);

	/* process irecv[], it contains a list of lists, the entries match
	 * the list of global indices with the REMOTE bit set */

	/* count # of elements to send to each process */
	for (ii = 0; ii < g->nprocs; ii++)
	    g->halo->scnt[ii] = 0;
	for (j = i = 0; i < n0; i++) {
	    ii = irecv[j++];	/* # ranks */
	    m = j + ii;
	    assert(m <= rcount);
	    for (; j < m; j++) {
		assert((size_t)g->halo->scnt[irecv[j]] +
					(size_t)size[send[i]] <= INT_MAX);
		g->halo->scnt[irecv[j]] += size[send[i]];	/* # elements */
	    }
	}

	scount = 0;
	for (ii = 0; ii < g->nprocs; ii++) {
	    g->halo->sdsp[ii] = scount;
	    scount += g->halo->scnt[ii];
	    g->halo->scnt[ii] = 0;		/* reset counter */
	}

	/* Now construct elocal[] */
	g->halo->elocal = phgAlloc(scount * sizeof(*g->halo->elocal));
	for (j = i = 0; i < n0; i++) {
	    k = send[i];
	    ii = irecv[j++];	/* # ranks */
	    m = j + ii;
	    for (; j < m; j++)
		for (jj = 0; jj < size[k]; jj++)
		    g->halo->elocal
			[g->halo->sdsp[irecv[j]] + (g->halo->scnt[irecv[j]]++)]
				= map[k][jj];
	}

	phgFree(send);
	phgFree(irecv);
	phgDofFreeVEFMap2(&ve);
    }
    else {
	/* halo_kind == HALO_FACE: construct elocal using RNEIGHBOUR info */
	memcpy(g->halo->scnt, g->neighbours.counts,
				g->nprocs * sizeof(*g->halo->scnt));
	memcpy(g->halo->sdsp, g->neighbours.displs,
				g->nprocs * sizeof(*g->halo->scnt));
	memcpy(g->halo->rcnt, g->neighbours.counts,
				g->nprocs * sizeof(*g->halo->scnt));
	memcpy(g->halo->rdsp, g->neighbours.displs,
				g->nprocs * sizeof(*g->halo->scnt));
	scount = g->neighbours.count;
	g->halo->elocal = phgAlloc(scount * sizeof(*g->halo->elocal));
	for (i = 0; i < g->neighbours.count; i++)
	    g->halo->elocal[i] = g->neighbours.list[i].local;
    }

    /* remove duplicate elements in elocal[] */
    phgInfo(2, "Before removing duplicate elements: scount = %d\n", scount);
    scount = 0;
    for (ii = 0; ii < g->nprocs; ii++) {
	if (g->halo->scnt[ii] <= 0) {
	    g->halo->sdsp[ii] = scount;
	    continue;
	}
	qsort(g->halo->elocal + g->halo->sdsp[ii], g->halo->scnt[ii],
		    sizeof(g->halo->elocal[0]), comp_pointer);
	i = g->halo->sdsp[ii];
	j = 0;
	g->halo->elocal[scount + j] = g->halo->elocal[i++];
	for (; i < g->halo->sdsp[ii] + g->halo->scnt[ii]; i++)
	    if (g->halo->elocal[i] != g->halo->elocal[scount + j])
		g->halo->elocal[scount + (++j)] = g->halo->elocal[i];
	g->halo->sdsp[ii] = scount;
	scount += (g->halo->scnt[ii] = j + 1);
    }
    phgInfo(2, "After removing duplicate elements: scount = %d\n", scount);

    MPI_Alltoall(g->halo->scnt, 1, MPI_INT, g->halo->rcnt, 1, MPI_INT,
		 g->comm);
    rcount = 0;
    for (ii = 0; ii < g->nprocs; ii++) {
	g->halo->rdsp[ii] = rcount;
	rcount += g->halo->rcnt[ii];
    }

#if 0
    phgInfo(-1, "%d elements sent:\n", scount);
    rank = 0;
    for (i = 0; i < scount; i++) {
	while (g->halo->sdsp[rank] + g->halo->scnt[rank] <= i)
	    rank++;
	e = g->halo->elocal[i];
	phgInfo(-1, "    e[%d] at %p to proc %d\n", e->index, e, rank);
    }
#endif

    /* prepare data to send */
    sbuf = phgAlloc(scount * sizeof(*sbuf));
    rank = 0;	/* target rank of the current element */
    for (i = 0; i < scount; i++) {
	int kk;
	BTYPE t;
	while (g->halo->sdsp[rank] + g->halo->scnt[rank] <= i)
	    rank++;
	assert(rank < g->nprocs);
	e = g->halo->elocal[i];
	for (kk = 0; kk < NVert; kk++) {
	    /*----- period->L2Gmap_vert values of the vertices */
	    sbuf[i].pindexes[kk] = GlobalVertexP(g, e->verts[kk]);
	    /*----- coordinates */
	    memcpy(sbuf[i].coord + kk, g->verts + e->verts[kk], sizeof(COORD));
	    /*----- types_vert[] */
	    sbuf[i].types[kk] = (t = g->types_vert[e->verts[kk]]) & ~OWNER;
	    /*----- owner_ranks[] and owner_indexes[] */
	    if ((t & OWNER)) {
		/* Note: corresponding entries in g->owner_* are -1 */
		sbuf[i].indexes[kk] = e->verts[kk];
		sbuf[i].ranks[kk] = g->rank;
	    }
	    else {
		sbuf[i].indexes[kk] = g->owner_index_vert[e->verts[kk]];
		sbuf[i].ranks[kk] = g->owner_rank_vert[e->verts[kk]];
	    }
	}
	/* the edges */
	for (jj = 0; jj < NEdge; jj++) {
	    sbuf[i].types[4 + jj] = (t = g->types_edge[e->edges[jj]]) & ~OWNER;
	    if ((t & OWNER)) {
		/* Note: corresponding entries in g->owner_* are -1 */
		sbuf[i].indexes[4 + jj] = e->edges[jj];
		sbuf[i].ranks[4 + jj] = g->rank;
	    }
	    else {
		sbuf[i].indexes[4 + jj] = g->owner_index_edge[e->edges[jj]];
		sbuf[i].ranks[4 + jj] = g->owner_rank_edge[e->edges[jj]];
	    }
	}
	/* the 4 faces */
	for (jj = 0; jj < NFace; jj++) {
	    sbuf[i].types[10 + jj] = (t = g->types_face[e->faces[jj]]) & ~OWNER;
	    if ((t & OWNER)) {
		/* Note: corresponding entries in g->owner_* are -1 */
		sbuf[i].indexes[10 + jj] = e->faces[jj];
		sbuf[i].ranks[10 + jj] = g->rank;
	    }
	    else {
		sbuf[i].indexes[10 + jj] = g->owner_index_face[e->faces[jj]];
		sbuf[i].ranks[10 + jj] = g->owner_rank_face[e->faces[jj]];
	    }
	}
	/* the element types */
	sbuf[i].types[14] = g->types_elem[e->index] & ~OWNER;
	sbuf[i].indexes[14] = e->index;
	sbuf[i].ranks[14] = g->rank;
	/*----- copy element data */
	memcpy(&sbuf[i].elem, e, sizeof(*e));
	/*----- save information about shared faces: if face "ii" has a remote
	 *	neighbour on process "rank", then the address of the peer
	 *	element is stored in elem.neighbours[ii], and the vertex and
	 *	op_vertex is saved in the 6 bits of elem.mark starting at the
	 *	bit position 6 * "ii". */
	assert(sizeof(sbuf[i].elem.mark) >= 3);	/* needs 4*6 == 24 bits */
	sbuf[i].elem.mark = 0;
	for (ii = 0; ii < NFace; ii++) {
	    RNEIGHBOUR *rn;
	    sbuf[i].elem.neighbours[ii] = NULL;
	    if (!HasRemoteNeighbour(e, ii))
		continue;
	    rn = g->neighbours.list + (size_t)e->neighbours[ii];
	    if (rn->rank != rank)
		continue;
	    sbuf[i].elem.neighbours[ii] = rn->remote;
	    /*----- save vertex + op_vertex in bits ["ii"*6:"ii"*6+5] of mark */
	    sbuf[i].elem.mark |= ((rn->vertex + rn->op_vertex * 8) << (ii * 6));
	}
	/*----- update with global indices */
	e = &sbuf[i].elem;
	for (ii = 0; ii < NVert; ii++)
	    e->verts[ii] = GlobalVertex(g, e->verts[ii]);
	for (ii = 0; ii < NEdge; ii++)
	    e->edges[ii] = GlobalEdge(g, e->edges[ii]);
	for (ii = 0; ii < NFace; ii++)
	    e->faces[ii] = GlobalFace(g, e->faces[ii]);
	e->index = GlobalElement(g, e->index);
    }

    /* exchange data */
    rbuf = phgAlloc(rcount * sizeof(*rbuf));
    MPI_Type_contiguous(sizeof(*sbuf), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(sbuf, g->halo->scnt, g->halo->sdsp, type,
		  rbuf, g->halo->rcnt, g->halo->rdsp, type, g->comm);
    MPI_Type_free(&type);
    phgFree(sbuf);

    /* process received data */

    /* buffer for saving new vertices/edges/faces/elems.
     * note: at most 6 new items per element. */
    items = phgAlloc(rcount * 6 * sizeof(*items));

    /* ordering buffer */
    n = g->nface;
    if (n < g->nedge)
	n = g->nedge;
    if (n < g->nvert)
	n = g->nvert;
    m = rcount * 6;
    ordering = phgAlloc((n > m ? n : m) * sizeof(*ordering));

    /* ------------------- process incoming vertices ------------------- */
    n = g->nvert;
    comp_base = g->L2Gmap_vert;
    for (i = 0; i < n; i++)
	ordering[i] = i;
    qsort(ordering, n, sizeof(*ordering), comp_ordering);
    m = 0;	/* number of (unsorted) new items */
    for (i = 0; i < rcount; i++) {
	e = &rbuf[i].elem;
	for (ii = 0; ii < NVert; ii++) {
	    /* k = global no of the vertex */
	    k = e->verts[ii];
	    /* search for k in comp_base[] */
	    j = phgBinarySearch(n, comp_base, ordering, &k, sizeof(*ordering),
				phgCompINT);
	    if (j < n && k == comp_base[ordering[j]] &&
		g->types_vert[ordering[j]] != UNREFERENCED) {
		/* found in comp_base[] */
		e->verts[ii] = ordering[j];
		continue;
	    }
	    items[m].gindex = k;
	    items[m].rbuf_index = i;
	    items[m].item_index = ii;
	    items[m].pindex = rbuf[i].pindexes[ii];
	    items[m].oindex = rbuf[i].indexes[ii];
	    items[m].orank = rbuf[i].ranks[ii];
	    items[m].type = rbuf[i].types[ii];
	    m++;
	}
    }
    if (m > 0) {
	/* sort list of items according to 'index' */
	for (i = 0; i < m; i++)
	    ordering[i] = i;
	/* sort list of items according to 'index' */
	qsort(ordering, m, sizeof(*ordering), comp_item);
	/* build halo->L2Gmap_vert[] and halo->verts[], and
	 * update g->{verts,L2Gmap_vert}, and e->verts[] entries in the halo */
	REALLOC(g->verts, g->nvert + m, g->nvert);
	REALLOC(g->L2Gmap_vert, g->nvert + m, g->nvert);
	if (g->period != NULL)
	    REALLOC(g->period->L2Gmap_vert, g->nvert + m, g->nvert);
	REALLOC(g->owner_index_vert, g->nvert + m, g->nvert);
	REALLOC(g->owner_rank_vert, g->nvert + m, g->nvert);
	REALLOC(g->types_vert, g->nvert + m, g->nvert);
	pi = items;
	n = g->nvert;
	for (i = 0; i < m; n++) {
	    pi = items + ordering[i];
	    memcpy(g->verts + n, rbuf[pi->rbuf_index].coord + pi->item_index,
		   sizeof(COORD));
	    j = g->L2Gmap_vert[n] = pi->gindex;
	    g->owner_index_vert[n] = pi->oindex;
	    g->owner_rank_vert[n] = pi->orank;
	    if (g->period != NULL) {
		g->period->L2Gmap_vert[n] = pi->pindex;
		if (g->owner_rank_vert[n] == g->rank) {
		    g->owner_index_vert[n] = -1;
		    g->owner_rank_vert[n] = -1;
		}
	    }
	    g->types_vert[n] = (pi->type & ~OWNER) | REMOTE;
	    assert(g->owner_rank_vert[n] != g->rank);
	    do {
		pi = items + ordering[i];
		rbuf[pi->rbuf_index].elem.verts[pi->item_index] = n;
	    } while (++i < m && items[ordering[i]].gindex == j);
	}
	m += g->nvert;
	if (n < m) {
	    REALLOC(g->verts, n, m);
	    REALLOC(g->L2Gmap_vert, n, m);
	    if (g->period != NULL)
		REALLOC(g->period->L2Gmap_vert, n, m);
	    REALLOC(g->owner_index_vert, n, m);
	    REALLOC(g->owner_rank_vert, n, m);
	    REALLOC(g->types_vert, n, m);
	}
	phgInfo(1, "%s: nvert += %d\n", __func__, n - g->nvert);
	g->halo->nvert_bak = g->nvert;
	g->nvert = n;
	if (g->period != NULL) {
	    phgFree(g->period->ordering);
	    g->period->ordering = phgAlloc(n * sizeof(*ordering));
	    comp_base = g->period->L2Gmap_vert;
	    for (i = 0; i < n; i++)
		g->period->ordering[i] = i;
	    qsort(g->period->ordering, n, sizeof(*ordering), comp_ordering);
	}
    }

    /* ------------------- process incoming edges ------------------- */
    n = g->nedge;
    comp_base = g->L2Gmap_edge;
    for (i = 0; i < n; i++)
	ordering[i] = i;
    qsort(ordering, n, sizeof(*ordering), comp_ordering);
    m = 0;	/* number of (unsorted) new items */
    for (i = 0; i < rcount; i++) {
	e = &rbuf[i].elem;
	/* Search L2Gmap_edge[] for the edges */
	for (jj = 0; jj < NEdge; jj++) {
	    /* k = global no of the edge */
	    k = e->edges[jj];
	    /* search for k in comp_base[] */
	    j = phgBinarySearch(n, comp_base, ordering, &k,
				sizeof(*ordering), phgCompINT);
	    if (j < n && k == comp_base[ordering[j]]
		      && g->types_edge[ordering[j]] != UNREFERENCED) {
		/* found in comp_base[] */
		e->edges[jj] = ordering[j];
		continue;
	    }
	    items[m].gindex = k;
	    items[m].rbuf_index = i;
	    items[m].item_index = jj;
	    items[m].oindex = rbuf[i].indexes[4 + jj];
	    items[m].orank = rbuf[i].ranks[4 + jj];
	    items[m].type = rbuf[i].types[4 + jj];
	    m++;
	}
    }
    if (m > 0) {
	/* sort list of items according to 'index' */
	for (i = 0; i < m; i++)
	    ordering[i] = i;
	/* sort list of items according to 'index' */
	qsort(ordering, m, sizeof(*ordering), comp_item);
	/* build halo->L2Gmap_edge[], update e->edges[] entries in the halo */
	REALLOC(g->L2Gmap_edge, g->nedge + m, g->nedge);
	REALLOC(g->owner_index_edge, g->nedge + m, g->nedge);
	REALLOC(g->owner_rank_edge, g->nedge + m, g->nedge);
	REALLOC(g->types_edge, g->nedge + m, g->nedge);
	pi = items;
	n = g->nedge;
	for (i = 0; i < m; n++) {
	    pi = items + ordering[i];
	    j = g->L2Gmap_edge[n] = pi->gindex;
	    g->owner_index_edge[n] = pi->oindex;
	    g->owner_rank_edge[n] = pi->orank;
	    g->types_edge[n] = (pi->type & ~OWNER) | REMOTE;
	    assert(g->owner_rank_edge[n] != g->rank);
	    do {
		pi = items + ordering[i];
		rbuf[pi->rbuf_index].elem.edges[pi->item_index] = n;
	    } while (++i < m && items[ordering[i]].gindex == j);
	}
	m += g->nedge;
	if (n < m) {
	    REALLOC(g->L2Gmap_edge, n, m);
	    REALLOC(g->owner_index_edge, n, m);
	    REALLOC(g->owner_rank_edge, n, m);
	    REALLOC(g->types_edge, n, m);
	}
	phgInfo(1, "%s: nedge += %d\n", __func__, n - g->nedge);
	g->halo->nedge_bak = g->nedge;
	g->nedge = n;
    }

    /* ------------------- process incoming faces ------------------- */
    n = g->nface;
    comp_base = g->L2Gmap_face;
    for (i = 0; i < n; i++)
	ordering[i] = i;
    qsort(ordering, n, sizeof(*ordering), comp_ordering);
    m = 0;	/* number of (unsorted) new items */
    for (i = 0; i < rcount; i++) {
	e = &rbuf[i].elem;
	/* Search L2Gmap_face[] for faces */
	for (jj = 0; jj < NFace; jj++) {
	    /* k = global no of the face */
	    k = e->faces[jj];
	    /* search for k in comp_base[] */
	    j = phgBinarySearch(n, comp_base, ordering, &k,
				sizeof(*ordering), phgCompINT);
	    if (j < n && k == comp_base[ordering[j]]
		      && g->types_face[ordering[j]] != UNREFERENCED) {
		/* found in comp_base[] */
		e->faces[jj] = ordering[j];
		continue;
	    }
	    items[m].gindex = k;
	    items[m].rbuf_index = i;
	    items[m].item_index = jj;
	    items[m].oindex = rbuf[i].indexes[10 + jj];
	    items[m].orank = rbuf[i].ranks[10 + jj];
	    items[m].type = rbuf[i].types[10 + jj];
	    m++;
	}
    }
    if (m > 0) {
	/* sort list of items according to 'index' */
	for (i = 0; i < m; i++)
	    ordering[i] = i;
	/* sort list of items according to 'index' */
	qsort(ordering, m, sizeof(*ordering), comp_item);
	/* build halo->L2Gmap_face[], update e->faces[] entries in the halo */
	REALLOC(g->L2Gmap_face, g->nface + m, g->nface);
	REALLOC(g->owner_index_face, g->nface + m, g->nface);
	REALLOC(g->owner_rank_face, g->nface + m, g->nface);
	REALLOC(g->types_face, g->nface + m, g->nface);
	pi = items;
	n = g->nface;
	for (i = 0; i < m; n++) {
	    pi = items + ordering[i];
	    j = g->L2Gmap_face[n] = pi->gindex;
	    g->owner_index_face[n] = pi->oindex;
	    g->owner_rank_face[n] = pi->orank;
	    g->types_face[n] = (pi->type & ~OWNER) | REMOTE;
	    assert(g->owner_rank_face[n] != g->rank);
	    do {
		pi = items + ordering[i];
		rbuf[pi->rbuf_index].elem.faces[pi->item_index] = n;
	    } while (++i < m && items[ordering[i]].gindex == j);
	}
	m += g->nface;
	if (n < m) {
	    REALLOC(g->L2Gmap_face, n, m);
	    REALLOC(g->owner_index_face, n, m);
	    REALLOC(g->owner_rank_face, n, m);
	    REALLOC(g->types_face, n, m);
	}
	phgInfo(1, "%s: nface += %d\n", __func__, n - g->nface);
	g->halo->nface_bak = g->nface;
	g->nface = n;
    }

    phgFree(items);

    /* ------------------- process incoming elements ------------------- */

    m = rcount;
    if (m > 0) {
	/* sort list of items according to 'index' */
	for (i = 0; i < m; i++)
	    ordering[i] = i;
	/* sort list of items according to 'index' */
	qsort(ordering, m, sizeof(*ordering), comp_elem);
	/* build halo->L2Gmap_elem[], update e->index in the halo */
	g->halo->map = phgAlloc(g->neighbours.count * sizeof(*g->halo->map));
#if DEBUG
	for (j = 0; j < g->neighbours.count; j++)
	    g->halo->map[j] = -1;
#endif
	g->halo->elist = phgAlloc(m * sizeof(*g->halo->elist));
	REALLOC(g->L2Gmap_elem, g->nelem + m, g->nelem);
	REALLOC(g->owner_index_elem, g->nelem + m, g->nelem);
	REALLOC(g->owner_rank_elem, g->nelem + m, g->nelem);
	REALLOC(g->types_elem, g->nelem + m, g->nelem);
	REALLOC(g->elems, g->nelem + m, g->nelem);
	n = g->nelem;
	for (i = 0; i < m; n++) {
	    e = g->halo->elist + n - g->nelem;
	    k = ordering[i];
	    memcpy(e, &rbuf[k].elem, sizeof(ELEMENT));
	    e->index = n;
	    e->mark = 0;
	    e->parent = NULL;
	    j = g->L2Gmap_elem[n] = rbuf[k].elem.index;
	    g->owner_index_elem[n] = rbuf[k].indexes[14];
	    g->owner_rank_elem[n] = rbuf[k].ranks[14];
	    g->types_elem[n] = (rbuf[k].types[14] & ~OWNER) | REMOTE;
	    assert(g->owner_rank_elem[n] != g->rank);
	    g->elems[n] = e;
	    /* Note: the loop below allows for duplicate elements in the
	     * incoming list (e.g., in the old code when
	     *		g->neighbours.list[].local
	     * are sent without compressing) */
	    do {
		int vertex, peer_vertex;
		ELEMENT *peer;
		k = ordering[i];
		/* check for shared faces */
		for (ii = 0; ii < NFace; ii++) {
		    if ((peer = rbuf[k].elem.neighbours[ii]) == NULL)
			continue;
		    vertex = (rbuf[k].elem.mark >> (6 * ii)) & 7;
		    peer_vertex = (rbuf[k].elem.mark >> (6 * ii + 3)) & 7;
		    /* update neighbour links */
		    e->neighbours[vertex] = peer;
		    assert(HasRemoteNeighbour(peer, peer_vertex));
		    g->halo->map[(size_t)peer->neighbours[peer_vertex]]
							= n - g->nelem;
		}
	    } while (++i < m && rbuf[ordering[i]].elem.index == j);
	}
#if DEBUG
	for (j = 0; j < g->neighbours.count; j++) {
	    int v0, v1, v2;
	    if (g->halo->map[j] != -1)
		continue;
	    phgInfo(-1, "unexpected error: missing element %p @ proc %d.\n",
		    g->neighbours.list[j].remote, g->neighbours.list[j].rank);
	    e = g->neighbours.list[j].local;
	    GetFaceVertices(e, g->neighbours.list[j].vertex, v0, v1, v2, ii);
	    phgInfo(-1, "shared face: (%"dFMT", %"dFMT", %"dFMT")\n",
		    GlobalVertex(g, e->verts[v0]),
		    GlobalVertex(g, e->verts[v1]),
		    GlobalVertex(g, e->verts[v2]));
	    phgError(-1, "%s:%d: abort.\n");
	}
#endif
	if (n < g->nelem + m) {
	    REALLOC(g->halo->elist, n - g->nelem, m);
	    m += g->nelem;
	    REALLOC(g->L2Gmap_elem, n, m);
	    REALLOC(g->owner_index_elem, n, m);
	    REALLOC(g->owner_rank_elem, n, m);
	    REALLOC(g->types_elem, n, m);
	    REALLOC(g->elems, n, m);
	    /* update pointers in g->elems[g->nelem:n] */
	    for (i = 0; i < n - g->nelem; i++) {
		e = g->halo->elist + i;
		assert(e->index >= g->nelem);
		g->elems[e->index] = e;
	    }
	}
	phgInfo(1, "%s: nelem += %d\n", __func__, n - g->nelem);
	g->halo->nelem_bak = g->nelem;
	g->nelem = n;
    }

    phgFree(ordering);
    phgFree(rbuf);

    /* update HP_TYPEs */
    for (php = g->hp; php != NULL && *php != NULL; php++) {
	HP_TYPE *hp = *php;
	REALLOC(hp->elem_order, g->nelem, g->halo->nelem_bak);
	ForAllElements_(g, e, FOR_UNOWNED)
	    hp->elem_order[e->index] = hp->info->max_order;
	phgHPSetupOrders(hp);
    }

    /* update DOF data */
    for (pdof = g->dof; pdof != NULL && *pdof != NULL; pdof++) {
	DOF *dof = *pdof;
	FLOAT *data, *p_old, *p_new;
	INT len_old;

	/* adjust data buffer */
	p_new = data = phgCalloc(DofDataCount(dof), sizeof(*data));
	p_old = dof->data;

	len_old = DofVertexDataCount_(dof, g->halo->nvert_bak);
	memcpy(p_new, p_old, len_old * sizeof(*data));
	dof->data_vert = p_new;
	p_old += len_old;
	p_new += DofVertexDataCount(dof);

	len_old = DofEdgeDataCount_(dof, g->halo->nedge_bak);
	memcpy(p_new, p_old, len_old * sizeof(*data));
	dof->data_edge = p_new;
	p_old += len_old;
	p_new += DofEdgeDataCount(dof);

	len_old = DofFaceDataCount_(dof, g->halo->nface_bak);
	memcpy(p_new, p_old, len_old * sizeof(*data));
	dof->data_face = p_new;
	p_old += len_old;
	p_new += DofFaceDataCount(dof);

	len_old = DofElementDataCount_(dof, g->halo->nelem_bak);
	memcpy(p_new, p_old, len_old * sizeof(*data));
	dof->data_elem = p_new;

	phgFree(dof->data);
	dof->data = data;

	/* update data */
	phgDofUpdateGhost_(dof, update_dofs);
    }
#else	/* USE_MPI */
    Unused(g);
    Unused(update_dofs);
#endif	/* USE_MPI */
    Unused(halo_type);

    /* rebuild face-to-element map */
    if (g->nface > 0 && g->Fmap_elems == NULL) {
	VEF_MAP *vef = phgDofSetupVEFMap_(g, NULL, FACE_FLAG, FOR_OWNED);
	phgDofFreeVEFMap(&vef);
    }
}

void
phgRemoveHalo(GRID *g)
{
    DOF **pdof;
    HP_TYPE **php;

    if (g->halo == NULL)
	return;

    g->serial_no += 1000 * 1000;

    if (g->nvert == g->halo->nvert_bak &&
	g->nedge == g->halo->nedge_bak &&
	g->nface == g->halo->nface_bak &&
	g->nelem == g->halo->nelem_bak)
	goto end;

    /* restore DOF data */
    for (pdof = g->dof; pdof!= NULL && *pdof != NULL; pdof++) {
	DOF *dof = *pdof;
	FLOAT *data, *p_old, *p_new;
	INT len, len_vert, len_edge, len_face, len_elem;
	if (dof->data == NULL || SpecialDofType(dof->type))
	    continue;

	phgInfo(1, "%s: update DOF %s\n", __func__, dof->name);

	/* new data sizes */
	len_vert = DofVertexDataCount_(dof, g->halo->nvert_bak);
	len_edge = DofEdgeDataCount_(dof, g->halo->nedge_bak);
	len_face = DofFaceDataCount_(dof, g->halo->nface_bak);
	len_elem = DofElementDataCount_(dof, g->halo->nelem_bak);

	/* allocate new data buffer */
	len = len_vert + len_edge + len_face + len_elem;
	p_new = data = phgAlloc(len * sizeof(*data));
	p_old = dof->data;

	memcpy(p_new, p_old, len_vert * sizeof(*data));
	dof->data_vert = p_new;
	p_old += DofVertexDataCount(dof);
	p_new += len_vert;

	memcpy(p_new, p_old, len_edge * sizeof(*data));
	dof->data_edge = p_new;
	p_old += DofEdgeDataCount(dof);
	p_new += len_edge;

	memcpy(p_new, p_old, len_face * sizeof(*data));
	dof->data_face = p_new;
	p_old += DofFaceDataCount(dof);
	p_new += len_face;

	memcpy(p_new, p_old, len_elem * sizeof(*data));
	dof->data_elem = p_new;

	phgFree(dof->data);
	dof->data = data;
    }

    /* restore HP_TYPEs */
    for (php = g->hp; php != NULL && *php != NULL; php++) {
	HP_TYPE *hp = *php;
	if ((hp->info->flags & VERT_FLAG)) {
	    REALLOC(hp->vert_order, g->halo->nvert_bak, g->nvert);
	    REALLOC(hp->vert_index, g->halo->nvert_bak + 1, g->nvert + 1);
	}
	if ((hp->info->flags & EDGE_FLAG)) {
	    REALLOC(hp->edge_order, g->halo->nedge_bak, g->nedge);
	    REALLOC(hp->edge_index, g->halo->nedge_bak + 1, g->nedge + 1);
	}
	if ((hp->info->flags & FACE_FLAG)) {
	    REALLOC(hp->face_order, g->halo->nface_bak, g->nface);
	    REALLOC(hp->face_index, g->halo->nface_bak + 1, g->nface + 1);
	}
	REALLOC(hp->elem_order, g->halo->nelem_bak, g->nelem);
	REALLOC(hp->elem_index, g->halo->nelem_bak + 1, g->nelem + 1);
    }

    /* invalidate face->element map */
    phgFree(g->Fmap_elems);
    g->Fmap_elems = NULL;
    phgFree(g->Fmap_faces);
    g->Fmap_faces = NULL;

    /* realloc buffers with new sizes */

#undef REALLOC
#define REALLOC(p) \
    p = phgRealloc_(p, g->halo->nvert_bak * sizeof(*p), g->nvert * sizeof(*p))
    if (g->nvert > g->halo->nvert_bak) {
	REALLOC(g->verts);
	REALLOC(g->types_vert);
#if USE_MPI
	REALLOC(g->L2Gmap_vert);
#endif	/* USE_MPI */
	REALLOC(g->owner_index_vert);
	REALLOC(g->owner_rank_vert);
	if (g->period != NULL) {
	    INT i, j;
	    REALLOC(g->period->L2Gmap_vert);
	    for (j = i = 0; i < g->nvert; i++) {
		if (g->period->ordering[i] >= g->halo->nvert_bak)
		    continue;
		g->period->ordering[j++] = g->period->ordering[i];
	    }
	    assert(j == g->halo->nvert_bak);
	    REALLOC(g->period->ordering);
	}
    }

#undef REALLOC
#define REALLOC(p) \
    p = phgRealloc_(p, g->halo->nedge_bak * sizeof(*p), g->nedge * sizeof(*p))
    if (g->nedge > g->halo->nedge_bak) {
	REALLOC(g->types_edge);
#if USE_MPI
	REALLOC(g->L2Gmap_edge);
	REALLOC(g->owner_index_edge);
	REALLOC(g->owner_rank_edge);
#endif	/* USE_MPI */
    }

#undef REALLOC
#define REALLOC(p) \
    p = phgRealloc_(p, g->halo->nface_bak * sizeof(*p), g->nface * sizeof(*p))
    if (g->nface > g->halo->nface_bak) {
	REALLOC(g->types_face);
#if USE_MPI
	REALLOC(g->L2Gmap_face);
#endif	/* USE_MPI */
    }

#undef REALLOC
#define REALLOC(p) \
    p = phgRealloc_(p, g->halo->nelem_bak * sizeof(*p), g->nelem * sizeof(*p))
    if (g->nelem > g->halo->nelem_bak) {
	REALLOC(g->elems);
	REALLOC(g->types_elem);
#if USE_MPI
	REALLOC(g->L2Gmap_elem);
#endif	/* USE_MPI */
    }

#if USE_MPI
    phgFree(g->owner_index_face);
    phgFree(g->owner_index_elem);
    g->owner_index_face = g->owner_index_elem = NULL;
    phgFree(g->owner_rank_face);
    phgFree(g->owner_rank_elem);
    g->owner_rank_face = g->owner_rank_elem = NULL;
#endif	/* USE_MPI */

    g->nvert = g->halo->nvert_bak;
    g->nedge = g->halo->nedge_bak;
    g->nface = g->halo->nface_bak;
    g->nelem = g->halo->nelem_bak;

end:
#if USE_MPI
    if (g->halo->scnt != g->neighbours.counts)
	phgFree(g->halo->scnt);
#endif	/* USE_MPI */
    phgFree(g->halo->elocal);
    phgFree(g->halo->elist);
    phgFree(g->halo->map);
    phgFree(g->halo->verts);
    phgFree(g->halo);
    g->halo = NULL;

    /* rebuild face-to-element map */
    if (g->nface > 0 && g->Fmap_elems == NULL) {
	VEF_MAP *vef = phgDofSetupVEFMap_(g, NULL, FACE_FLAG, FOR_OWNED);
	phgDofFreeVEFMap(&vef);
    }
}
