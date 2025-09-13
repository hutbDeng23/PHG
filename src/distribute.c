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

/* $Id: distribute.c,v 1.221 2022/06/24 05:13:30 zlb Exp $ */

#include "phg.h"
#include "phg/partition-utils.h"

#if USE_MPI	/* whole file */

#include "phg/partition-rtk.h"
#if USE_PARMETIS
#include "phg/partition-metis.h"
#endif	/* USE_PARMETIS */
#if USE_ZOLTAN
#include "phg/partition-zoltan.h"
#endif	/* USE_ZOLTAN */

static const PARTITIONER partitioners_list[] = {
    phgPartitionSFC,
    phgPartitionRTK,
    phgPartitionerRandom,
    phgPartitionerUser,
#if USE_PARMETIS
    phgPartitionParMetis,
#endif	/* USE_PARMETIS */
#if USE_ZOLTAN
    phgPartitionZoltan
#endif	/* USE_ZOLTAN */
};

static const char *partitioners_name[] = {
    "sfc",
    "rtk",
    "random",
    "user",
#if USE_PARMETIS
    "metis",
#endif	/* USE_PARMETIS */
#if USE_ZOLTAN
    "zoltan",
#endif	/* USE_ZOLTAN */
    NULL
};

static int partitioner = 0;
static int partitioners_count = sizeof(partitioners_list) /
				sizeof(partitioners_list[0]);

/* use e->parent as e->peer, which will be restored by phgUpdateBoundaryTypes */
# define peer	parent

static DOF *weights = NULL;
static FLOAT power = 0.;

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>	/* PATH_MAX, INT_MAX */

static GRID *g_;
static size_t pack_count;
static BOOLEAN has_idle_processes = FALSE;

PARTITIONER
phgGetPartitioner(int no)
{
    return no < partitioners_count ? partitioners_list[no] : NULL;
}

int
phgCreateComm_(int nprocs, MPI_Comm comm_world, MPI_Comm *comm)
/* Create a communicator containing the processes [0, nprocs).
 * Note: the idle processes will be looping in this function, until
 * they join a new communicator.
 * This function will also be called by phgFinalizeMPI() with nprocs == -1
 * to tell the idle processes to exit */
{
    int nprocs_max, rank;

#if NO_NEW_COMMUNICATOR
    assert(nprocs == phgNProcs || nprocs == 1 || comm == NULL);
#endif
    phgInfo(3, "phgCreateComm_: nprocs = %d, has_idle_processes = %d\n",
		nprocs, has_idle_processes);
    if (nprocs <= 0 && !has_idle_processes)
	return 0;

    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &nprocs_max);
    do {
	MPI_Bcast(&nprocs, 1, MPI_INT, 0, comm_world);
	if (nprocs <= 0) {
	    has_idle_processes = FALSE;	/* avoid infinite loops */
	    if (comm != NULL)
		*comm = MPI_COMM_NULL;
	    phgFinalize();
	    return 0;
	}
#if NO_NEW_COMMUNICATOR
	*comm = (nprocs == 1 ? MPI_COMM_SELF : comm_world);
	if (rank >= nprocs)
	    *comm = MPI_COMM_NULL;
#else	/* NO_NEW_COMMUNICATOR */
	if (nprocs < nprocs_max) {
	    MPI_Comm_split(comm_world, rank < nprocs ? 1 : MPI_UNDEFINED,
				0, comm);
	    assert(rank >= nprocs || *comm != MPI_COMM_NULL);
	}
	else {
	    MPI_Comm_dup(comm_world, comm);
	}
#endif	/* NO_NEW_COMMUNICATOR */
    } while (!phgMasterSlave && *comm == MPI_COMM_NULL);

    has_idle_processes = (nprocs < nprocs_max);

    return nprocs;
}

static int rebuild_type;
static INT hash_size, map_count;
static struct {
    struct {
	INT local;
	INT global;
    }  *list;
    short count, alloc;
} *hash_table, *hash;

static BOOLEAN
map_callback CB_ARGS(e)
{
    int i, j, n = 0;
    INT index, *loc = NULL;

    switch (rebuild_type) {
	case VERTEX:
	    loc = e->verts;
	    n = NVert;
	    break;
	case EDGE:
	    loc = e->edges;
	    n = NEdge;
	    break;
#if Dim == 3
	case FACE:
	    loc = e->faces;
	    n = NFace;
	    break;
#endif
	case VOLUME:
	    loc = &e->index;
	    n = 1;
	    break;
    }

    for (i = 0; i < n; i++) {
	index = loc[i];
	if (index < 0) {
	    phgInfo(-1, "======================================\n");
	    phgDumpElement(g_, e);
	    phgError(1, "%s:%d, invalid index.\n", __FILE__, __LINE__);
	}
	hash = hash_table + (index % hash_size);
	for (j = 0; j < hash->count; j++)
	    if (hash->list[j].global == index)
		break;
	if (j < hash->count) {
	    loc[i] = hash->list[j].local;
	    continue;
	}
	/* add a new entry to map */
	if (hash->count >= hash->alloc) {
	    assert(hash->alloc < SHRT_MAX);
	    hash->list = phgRealloc_(hash->list,
				     (hash->alloc += 1) * sizeof(*(hash->list)),
				     j * sizeof(*(hash->list)));
	}
	hash->list[j].global = index;
	loc[i] = hash->list[hash->count++].local = map_count++;
    }

    return TRUE;
}

void
phgRebuildL2Gmap_(GRID *g, GTYPE type)
/* rebuilds an L2Gmap array, and replace indices in elements
   with local ones.

   `type' specifies which map to build: VERTEX, EDGE, FACE, VOLUME

   WARNING: this function requires that the elements have global indices */
{
    INT i, nglobal = 0, **map = NULL;
    int j;
    void (*traverse)(GRID *g, BOOLEAN (*cb) CB_ARGS(e)) = NULL;
    double t, t0 = phgGetTime(NULL);

    if (g->nprocs <= 1)
	return;

    g_ = g;

    /* Note: g->nvert/nedge/nface/nelem may not be up to date at this stage,
     *	     so g->nleaf is used to compute the hash size */
    switch (rebuild_type = type) {
	case VERTEX:
	    traverse = phgTraverseAllElements;
	    map = &g->L2Gmap_vert;
	    nglobal = g->nvert_global;
	    hash_size = g->nleaf / 5;
	    break;
	case EDGE:
	    traverse = phgTraverseAllElements;
	    map = &g->L2Gmap_edge;
	    nglobal = g->nedge_global;
	    hash_size = (g->nleaf * 6) / 5;
	    break;
#if Dim == 3
	case FACE:
	    traverse = phgTraverseAllElements;
	    map = &g->L2Gmap_face;
	    nglobal = g->nface_global;
	    hash_size = g->nleaf * 2;
	    break;
#endif
	case VOLUME:
	    traverse = phgTraverseAllElements;
	    map = &g->L2Gmap_elem;
	    nglobal = g->nelem_global;
	    hash_size = g->nleaf;
	    break;
    }

    /* allocate hash table */
    if (nglobal == 0)
	return;

    if (hash_size <= 0)
	hash_size = 1;
    hash_table = phgCalloc(hash_size, sizeof(*hash_table));

    map_count = 0;
    (*traverse)(g, map_callback);

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s (type=%d) - traverse: %0.4lf\n", __func__, type, t0);
	t0 = phgGetTime(NULL);
    }

    /* build L2Gmap and free hash table */
    if (*map != NULL)
	phgFree(*map);
    *map = phgAlloc(map_count * sizeof(**map));
    for (i = 0; i < hash_size; i++) {
	if ((hash = hash_table + i)->list == NULL)
	    continue;
	for (j = 0; j < hash->count; j++) {
#if DEBUG
	    if (hash->list[j].global < 0 || hash->list[j].global >= nglobal)
		phgError(1, "%s:%d, invalid map %d => %d\n", __FILE__, __LINE__,
				hash->list[j].local, hash->list[j].global);
#endif
	    (*map)[hash->list[j].local] = hash->list[j].global;
	}
	phgFree(hash->list);
    }
    phgFree(hash_table);

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s (type=%d) - building: %0.4lf\n", __func__, type, t0);
	t0 = phgGetTime(NULL);
    }

    /* build L2Gmap and free hash table */
    switch (rebuild_type) {
	case VERTEX:
	    g->nvert = map_count;
	    break;
	case EDGE:
	    g->nedge = map_count;
	    break;
#if Dim == 3
	case FACE:
	    g->nface = map_count;
	    break;
#endif
	case VOLUME:
	    g->nelem = map_count;
	    break;
    }
}

typedef struct {
    RNEIGHBOUR rn;		/* the RNEIGHBOUR data */
    INT gindex[NFace - 1];	/* global indices (sorted) */
} RNEIGHBOUR_t;

#define SORT_INDEX(index0, index1)					\
{									\
    INT tmp_;								\
    if (index0[0] > index0[1]) {					\
	tmp_ = index0[0]; index0[0] = index0[1]; index0[1] = tmp_;	\
	tmp_ = index1[0]; index1[0] = index1[1]; index1[1] = tmp_;	\
    }									\
    if (Dim == 3) {							\
	if (index0[0] > index0[2]) {					\
	    tmp_ = index0[0]; index0[0] = index0[2]; index0[2] = tmp_;	\
	    tmp_ = index1[0]; index1[0] = index1[2]; index1[2] = tmp_;	\
	}								\
	if (index0[1] > index0[2]) {					\
	    tmp_ = index0[1]; index0[1] = index0[2]; index0[2] = tmp_;	\
	    tmp_ = index1[1]; index1[1] = index1[2]; index1[2] = tmp_;	\
	}								\
    }									\
}

#if 0	/*----------- old algorithm for updating remote neighbours ----------*/

static RNEIGHBOUR_t *rn_buffer0, *rn_buffer1, *rn_buffer2;

static BOOLEAN
updrn_count_rn CB_ARGS(e)
{
    int i;

    for (i = 0; i < NFace; i++)
	if ((e->bound_type[i] & INTERIOR) && e->neighbours[i] == NULL)
	    pack_count++;

    return TRUE;
}

static BOOLEAN
updrn_collect_rn CB_ARGS(e)
{
    int i, j, k;
    RNEIGHBOUR_t *rn;

    for (i = 0; i < NFace; i++) {
	if (!(e->bound_type[i] & INTERIOR) || e->neighbours[i] != NULL) {
	    e->bound_type[i] &= ~REMOTE;
	    continue;
	}
	e->bound_type[i] |= REMOTE;
	rn = rn_buffer0 + (pack_count++);
	rn->rn.remote = NULL;
	rn->rn.local = e;
	rn->rn.vertex = i;
	for (j = k = 0; j < NVert; j++) {
	    if (j == i)
		continue;
	    rn->rn.rface[k] = j;
	    rn->gindex[k] = GlobalVertexP(g_, e->verts[j]);
	    k++;
	}
	SORT_INDEX(rn->gindex, rn->rn.rface)
    }

    return TRUE;
}

static int
updrn_comp_rn(const void *p0, const void *p1)
{
    const INT *g0 = ((const RNEIGHBOUR_t *)p0)->gindex,
	      *g1 = ((const RNEIGHBOUR_t *)p1)->gindex;
    int i, j;

    for (i = 0; i < NVert - 1; i++)
	if ((j = g0[i] - g1[i]) != 0)
	    return j;

    return 0;
}

static void
updrn_rface_rn(RNEIGHBOUR_t *local, int nlocal, RNEIGHBOUR_t *remote,
	int nremote, int rank)
/* update the local list with entries of the remote list
   (local and remote are sorted lists) */
{
    int i, j, k = 0;
    BYTE tmp[NVert - 1];

    phgDebug((3, "rfaceing faces from proc %d:\n", rank));
    i = j = 0;
    while (i < nlocal && j < nremote) {
	while (i < nlocal && (local->rn.remote != NULL ||
		updrn_comp_rn(local, remote) < 0)) {
	    i++; local++;
	}
	if (i >= nlocal)
	    break;

	while (j < nremote && (/*remote->rn.remote != NULL ||*/
		(k = updrn_comp_rn(remote, local)) < 0)) {
	    j++; remote++;
	}
	if (j >= nremote)
	    break;

	if (k != 0)
	    continue;

	/* update local data */
	phgDebug((4, "rfaceed face: %3d %3d %3d, remote=%p\n",
	    local->gindex[0], local->gindex[1], local->gindex[2],
	    remote->rn.local));
	local->rn.remote	= remote->rn.local;
	remote->rn.remote	= local->rn.local;
	local->rn.op_vertex	= remote->rn.vertex;
	local->rn.rank		= rank;
	memcpy(tmp, local->rn.rface, sizeof(tmp));
	memcpy(local->rn.rface, remote->rn.rface, sizeof(tmp));
	SORT_INDEX(tmp, local->rn.rface)
	i++; local++;
	j++; remote++;
    }
}

void
phgUpdateRNeighbours_(GRID *g)
{
    int i, j;
    int src, dst, i1, i2;
    RNEIGHBOUR_t *rp1, *rp2;
    int *counts;
    MPI_Datatype type;
    MPI_Status  status[2];
    MPI_Request request[2];

    g_ = g;

    if (g->neighbours.counts == NULL)
	 g->neighbours.counts = phgCalloc(g->nprocs,
					  sizeof(*g->neighbours.counts));
    if (g->neighbours.displs == NULL)
	 g->neighbours.displs = phgCalloc(g->nprocs,
					  sizeof(*g->neighbours.displs));

    counts = g->neighbours.counts;

    /* collect REMOTE faces of leaf elements */
    pack_count = 0;
    phgTraverseElements(g, updrn_count_rn);
    rn_buffer0 = phgAlloc(pack_count * sizeof(*rn_buffer0));

    pack_count = 0;
    phgTraverseElements(g, updrn_collect_rn);

    /* sort the list of rneighbours */
    if (pack_count > 0)
	qsort(rn_buffer0, pack_count, sizeof(*rn_buffer0), updrn_comp_rn);

    /* check and remove duplicate entries (may come from coarsen.c)
     * (mpirun -np 2 coarsen_test -refines0=8, with phgPartitionGrid__) */
    for (j = -1, i = 0; i < pack_count; i++) {
	if (i < pack_count - 1 &&
		updrn_comp_rn(rn_buffer0 + i, rn_buffer0 + i + 1) == 0) {
	    /* local neighbours */
	    rp1 = rn_buffer0 + i;
	    rp2 = rn_buffer0 + i + 1;
	    i1 = rp1->rn.vertex;
	    i2 = rp2->rn.vertex;
	    rp1->rn.local->bound_type[i1] &= ~REMOTE;
	    rp2->rn.local->bound_type[i2] &= ~REMOTE;
	    rp1->rn.local->neighbours[i1] = rp2->rn.local;
	    rp2->rn.local->neighbours[i2] = rp1->rn.local;
	    phgDebug((2, "%s:%d, locally matched face: %p (%d)<->%p (%d), "
			 "%d %d %d\n", __FILE__, __LINE__,
			 rp1->rn.local, i1, rp2->rn.local, i2,
			 rp1->gindex[0], rp1->gindex[1], rp1->gindex[2]));
	    i++;
	    continue;
	}
	if (++j != i) {
	    memcpy(rn_buffer0 + j, rn_buffer0 + i, sizeof(*rn_buffer0));
	}
    }
    pack_count = j + 1;

    phgDebug((3, "number of remote neighbours: %d\n", pack_count));

    /* gather number of rneighbours */
    i = pack_count;
    MPI_Allgather(&i, 1, MPI_INT, counts, 1, MPI_INT, g->comm);

    /* compute maxi buffer length */
    j = 0;
    for (i = 0; i < g->nprocs; i++)
	if (j < counts[i])
	    j = counts[i];

    rn_buffer1 = phgAlloc(j * sizeof(*rn_buffer0));
    rn_buffer2 = phgAlloc(j * sizeof(*rn_buffer1));

#if DEBUG
    if (phgVerbosity >= 4)
      for (rp1 = rn_buffer0; rp1 < rn_buffer0 + pack_count; rp1++) {
	ELEMENT *e = rp1->rn.local;
	phgInfo(0, "%3d: %p (g%-2d g%-2d g%-2d g%-2d) %2d %2d %2d, "
		"v=%d, mat=%d%d%d\n",
		rp1 - rn_buffer0, e, GlobalVertexP(g, e->verts[0]),
		GlobalVertexP(g, e->verts[1]), GlobalVertexP(g, e->verts[2]),
		GlobalVertexP(g, e->verts[3]), rp1->gindex[0], rp1->gindex[1],
		rp1->gindex[2], rp1->rn.vertex, rp1->rn.rface[0],
		rp1->rn.rface[1], rp1->rn.rface[2]);
      }
#endif

    /* circulate the lists of rneighbours on the processes ring */
    if ((src = g->rank - 1) < 0)
	src += g->nprocs;
    if ((dst = g->rank + 1) >= g->nprocs)
	dst -= g->nprocs;

    MPI_Type_contiguous(sizeof(*rp1), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    memcpy(rn_buffer1, rn_buffer0, pack_count * sizeof(*rn_buffer0));
    rp1 = rn_buffer1;
    rp2 = rn_buffer2;
    i1 = g->rank;
    i2 = src;
    for (j = 0; j < g->nprocs - 1; j++) {
	MPI_Isend(rp1, counts[i1], type, dst, 1, g->comm, request + 0);
	MPI_Irecv(rp2, counts[i2], type, src, 1, g->comm, request + 1);
	if (j != 0) {
	    /* process list in rp1 */
	    updrn_rface_rn(rn_buffer0, counts[g->rank], rp1, counts[i1], i1);
	}
	MPI_Waitall(2, request, status);
	i1 = i2;
	if (--i2 < 0)
	    i2 += g->nprocs;

	if (rp1 == rn_buffer1) {
	    rp1 = rn_buffer2;
	    rp2 = rn_buffer1;
	}
	else {
	    rp1 = rn_buffer1;
	    rp2 = rn_buffer2;
	}
    }
    /* process last list in rp1 */
    updrn_rface_rn(rn_buffer0, counts[g->rank], rp1, counts[i1], i1);

    MPI_Type_free(&type);

    phgFree(rn_buffer1);
    phgFree(rn_buffer2);

    pack_count = counts[g->rank];
    if (g->neighbours.list != NULL)
	phgFree(g->neighbours.list);
    g->neighbours.list = phgAlloc(pack_count * sizeof(*g->neighbours.list));
    g->neighbours.count = g->neighbours.allocated = pack_count;
    for (i = 0; i < pack_count; i++) {
	g->neighbours.list[i] = rn_buffer0[i].rn;
	if (g->neighbours.list[i].remote == NULL) {
	    phgFree(g->neighbours.list);
	    g->neighbours.list = NULL;
	    g->neighbours.count = 0;
	    phgDumpElement(g, rn_buffer0[i].rn.local);
	    phgError(1, "%s:%d, element %p, face: %d %d %d\n",
			__FILE__, __LINE__, rn_buffer0[i].rn.local,
			rn_buffer0[i].gindex[0], rn_buffer0[i].gindex[1],
			rn_buffer0[i].gindex[2]);
	}
    }

    phgFree(rn_buffer0);

    g->neighbours.count = phgCountRNeighbours(g->nprocs, g->neighbours.list,
				pack_count, counts, g->neighbours.displs);

    for (i = 0; i < pack_count; i++) {
	RNEIGHBOUR *rn = g->neighbours.list + i;
	rn->local->neighbours[rn->vertex] = (void *)((size_t)i);
	phgDebug((4, "rn[%d]: %p => %p@p%d, v=%d, op=%d, rface=%d%d%d\n", i,
		rn->local, rn->remote, rn->rank, rn->vertex,
		rn->op_vertex, rn->rface[0], rn->rface[1], rn->rface[2]));
    }

    return;
}

#else	/*----------- new algorithm for updating remote neighbours ----------*/

static RNEIGHBOUR_t *rn_buff;
static int *rn_counts;

static BOOLEAN
updrn_callback CB_ARGS(e)
{
    int i, j;
    INT k;

    g_->elems[phgTraverseIndex] = e;
    for (i = 0; i < NFace; i++) {
	if (!(e->bound_type[i] & INTERIOR) || e->neighbours[i] != NULL) {
	    e->bound_type[i] &= ~REMOTE;
	    continue;
	}
	for (j = k = 0; j < NVert; j++) {
	    if (j != i)
		k += GlobalVertexP(g_, e->verts[j]);
	}
	assert(rn_counts[k % g_->nprocs] < INT_MAX);
	rn_counts[k % g_->nprocs]++;
    }

    return TRUE;
}

static int
updrn_comp_rn(const void *p0, const void *p1)
{
    INT *g0 = rn_buff[*((const int *)p0)].gindex;
    INT *g1 = rn_buff[*((const int *)p1)].gindex;
    INT j;
    int i;

    for (i = 0; i < NVert - 1; i++)
	if ((j = g0[i] - g1[i]) != 0)
	    return j > 0 ? 1 : -1;

    return 0;
}

void
phgUpdateRNeighbours_(GRID *g)
{
    INT i, k, m, n;
    int ii, jj, nn, *scnts, *sdsps, *rcnts, *rdsps, *order;
    MPI_Datatype type;
    RNEIGHBOUR_t *sbuf, *rbuf, *rp1, *rp2;
    RNEIGHBOUR *sbuf0, *rbuf0;
    ELEMENT *e, **elems_save;
    double time0, time1;

    if (g->nprocs <= 1)
	return;

    g_ = g;

    if (g->neighbours.counts == NULL)
	 g->neighbours.counts = phgAlloc(g->nprocs * sizeof(*scnts));
    if (g->neighbours.displs == NULL)
	 g->neighbours.displs = phgAlloc(g->nprocs * sizeof(*sdsps));

    scnts = g->neighbours.counts;
    sdsps = g->neighbours.displs;

    rcnts = phgAlloc(2 * g->nprocs * sizeof(rcnts));
    rdsps = rcnts + g->nprocs;

    /* collect REMOTE faces of leaf elements */
    time0 = phgGetTime(NULL);
    elems_save = g->elems;
    rn_counts = scnts;
    g->elems = phgAlloc(g->nleaf * sizeof(*g->elems));
    memset(scnts, 0, g->nprocs * sizeof(*scnts));
    phgTraverseElements(g, updrn_callback);

    m = 0;
    for (ii = 0; ii < g->nprocs; ii++) {
	assert(m < INT_MAX);
	sdsps[ii] = m;
	m += scnts[ii];
    }

    /* collect rfaces */
    sbuf = phgAlloc(m * sizeof(*sbuf));
    memset(scnts, 0, g->nprocs * sizeof(*scnts));
    for (k = 0; k < g->nleaf; k++) {
	RNEIGHBOUR_t r;
	e = g->elems[k];
	for (ii = 0; ii < NFace; ii++) {
	    if (!(e->bound_type[ii] & INTERIOR) || e->neighbours[ii] != NULL)
		continue;
	    e->bound_type[ii] |= REMOTE;
	    r.rn.remote = NULL;
	    r.rn.local = e;
	    r.rn.vertex = ii;
	    r.rn.rank = g->rank;
	    for (jj = nn = 0; jj < NVert; jj++) {
		if (jj == ii)
		    continue;
		r.rn.rface[nn] = jj;
		r.gindex[nn] = GlobalVertexP(g_, e->verts[jj]);
		nn++;
	    }
	    SORT_INDEX(r.gindex, r.rn.rface)
	    jj = (r.gindex[0] + r.gindex[1] + r.gindex[2]) % g->nprocs;
	    /*assert(scnts[jj] < INT_MAX - sdsps[jj]);*/
	    sbuf[sdsps[jj] + scnts[jj]++] = r;
	}
    }
    phgFree(g->elems);
    g->elems = elems_save;

    time1 = phgGetTime(NULL);
    if (phgVerbosity > 1) {
	double time = time1 - time0;
	MPI_Reduce(&time, &time0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s: collect = %lg\n", __func__, time0);
	time1 = phgGetTime(NULL);
    }
    time0 = time1;

    MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, g->comm);

    time1 = phgGetTime(NULL);
    if (phgVerbosity > 1)
	phgPrintf("%s: MPI_Alltoall = %lg\n", __func__, time1 - time0);
    time0 = time1;

    n = 0;	/* n is the size of rbuf */
    for (ii = 0; ii < g->nprocs; ii++) {
	rdsps[ii] = n;
	n += rcnts[ii];
	assert(n < INT_MAX);
    }
    rbuf = phgAlloc(n * sizeof(*rbuf));
    order = phgAlloc(n * sizeof(*order));

    phgDebug((2, "%s: send = %d, recv = %d\n", __func__, m, n));

    /* send lists to destination processes */
    MPI_Type_contiguous(sizeof(*sbuf), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(sbuf, scnts, sdsps, type, rbuf, rcnts, rdsps, type, g->comm);
    MPI_Type_free(&type);

    time1 = phgGetTime(NULL);
    if (phgVerbosity > 1)
	phgPrintf("%s: 1st MPI_Alltoallv = %lg\n", __func__, time1 - time0);
    time0 = time1;

    /* sort received lists */
    if (n > 0) {
	for (i = 0; i < n; i++)
	    order[i] = i;
	rn_buff = rbuf;
	qsort(order, n, sizeof(*order), updrn_comp_rn);
    }

    /* process received lists */
    for (i = 0; i < n; i++) {
	BYTE tmp1[NVert - 1], tmp2[NVert - 1];

	rp1 = rbuf + order[i];
	if (i >= n - 1 || updrn_comp_rn(order + i, order + i + 1)) {
	    assert(FALSE);	/* shouldn't happen */
	    rp1->rn.rank = -1;
	    continue;
	}

	/* order[i] and order[i + 1] are matching entries */
	rp2 = rbuf + order[++i];

	assert(i >= n - 1 || updrn_comp_rn(order + i, order + i + 1));

	ii = rp1->rn.rank;
	rp1->rn.rank = rp2->rn.rank;
	rp2->rn.rank = ii;

	rp1->rn.remote = rp2->rn.local;
	rp2->rn.remote = rp1->rn.local;

	if (rp1->rn.rank == rp2->rn.rank) {
	    assert(FALSE);	/* shouldn't happen */
	    continue;
	}

	rp1->rn.op_vertex = rp2->rn.vertex;
	rp2->rn.op_vertex = rp1->rn.vertex;

	memcpy(tmp1, rp1->rn.rface, sizeof(tmp1));
	memcpy(tmp2, rp2->rn.rface, sizeof(tmp2));
	memcpy(rp1->rn.rface, tmp2, sizeof(tmp2));
	memcpy(rp2->rn.rface, tmp1, sizeof(tmp2));
	SORT_INDEX(tmp1, rp1->rn.rface)
	SORT_INDEX(tmp2, rp2->rn.rface)
    }

    time1 = phgGetTime(NULL);
    if (phgVerbosity > 1)
	phgPrintf("%s: processing = %lg\n", __func__, time1 - time0);
    time0 = time1;

    /* send back lists (only send back the rn member) */
    sbuf0 = (RNEIGHBOUR *)sbuf;
    rbuf0 = (RNEIGHBOUR *)rbuf;
    for (i = 0; i < n; i++) {
#if 0
	/* Note: this doesn't work with MacPorts gcc MacOSX */
	rbuf0[i] = rbuf[i].rn;
#else
	memmove(&rbuf0[i], &rbuf[i].rn, sizeof(rbuf[i].rn));
#endif
    }
    MPI_Type_contiguous(sizeof(*sbuf0), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(rbuf0, rcnts, rdsps, type, sbuf0, scnts, sdsps, type,
		  g->comm);
    MPI_Type_free(&type);

    time1 = phgGetTime(NULL);
    if (phgVerbosity > 1)
	phgPrintf("%s: 2nd MPI_Alltoallv = %lg\n", __func__, time1 - time0);
    time0 = time1;

    if (g->neighbours.list != NULL)
	phgFree(g->neighbours.list);
    g->neighbours.list = phgAlloc(m * sizeof(*g->neighbours.list));
    g->neighbours.count = g->neighbours.allocated = m;
    k = 0;
    for (i = n = 0; i < m; i++) {
#if 0
	RNEIGHBOUR *rn = sbuf0 + i;
	if (rn->rank == -1) {
	    /* unmatched rface, assuming to be a boundary face
	     * (may be encountered with coarsen_test1.c?) */
	    phgWarning("%s, unmatched face: elem %d, face %d.\n",
			    __func__, GlobalElement(g, rn->local->index),
			    rn->vertex);
	    rn->local->bound_type[rn->vertex] &= ~(INTERIOR | REMOTE);
	    /*assert((rn->local->bound_type[rn->vertex] & BDRY_MASK));*/
	    rn->local->neighbours[rn->vertex] = NULL;
	    k |= 1;
	    continue;
	}
	else if (rn->rank == g->rank) {
	    /* local neighbour */
	    phgWarning("%s, local neighbour: elem %d, face %d, neighbour %d.\n",
			    __func__, GlobalElement(g, rn->local->index),
			    rn->vertex, GlobalElement(g, rn->remote->index));
	    rn->local->bound_type[rn->vertex] &= ~REMOTE;
	    rn->local->bound_type[rn->vertex] |= INTERIOR;
	    rn->local->neighbours[rn->vertex] = rn->remote;
	    k |= 2;
	    continue;
	}
#endif	/* 0 */
	g->neighbours.list[n++] = sbuf0[i];
    }

#if 0 * DEBUG
#warning comment out me!
    MPI_Allreduce(&k, &i, 1, PHG_MPI_INT, MPI_BOR, g->comm);
    if (i != 0) {
	/* Note: error produced at i=12 by:
	 *    -np 16 coarsen_test1 -partitioner zoltan -zoltan_method block */
	phgPrintf("Dumping submeshes ...\n");
	phgDumpGrid(g);
	phgError(1, "abort.\n");
	phgFinalize();
	exit(1);
    }
#endif	/* DEBUG */

    g->neighbours.count = phgCountRNeighbours(g->nprocs, g->neighbours.list,
				n, scnts, sdsps);

    for (i = 0; i < n; i++) {
	RNEIGHBOUR *rn = g->neighbours.list + i;
	rn->local->neighbours[rn->vertex] = (void *)((size_t)i);
    }

    time1 = phgGetTime(NULL);
    if (phgVerbosity > 1)
	phgPrintf("%s: update g->neighbours = %lg\n", __func__, time1 - time0);
    time0 = time1;

    phgFree(order);
    phgFree(sbuf);
    phgFree(rbuf);
    phgFree(rcnts);

    return;
}

#endif	/*-----------  algorithm for updating remote neighbours ----------*/

static ELEMENT *pack_base;
size_t pack_offset;	/* offset of pack_base with resp. to process 0 */

static INT *iflags;

static BOOLEAN
vertices_callback0 CB_ARGS(e)
{
    int i;
    for (i = 0; i < NVert; i++)
	iflags[e->verts[i]] = 1;
    return TRUE;
}

static BOOLEAN
vertices_callback1 CB_ARGS(e)
{
    int i;
    for (i = 0; i < NVert; i++)
	e->verts[i] = iflags[e->verts[i]];
    return TRUE;
}

static void
pack_vertices(GRID *g)
{
    INT i, n;

    iflags = phgCalloc(g->nvert, sizeof(*(iflags)));

    /* mark referenced vertices */
    phgTraverseAllElements(g, vertices_callback0);

    /* renumber vertices */
    n = 0;
    for (i = 0; i < g->nvert; i++)
	iflags[i] = (iflags[i] != 0) ? n++ : -1;

    /* update vertex indices of all elements */
    phgTraverseAllElements(g, vertices_callback1);

    if (g->L2Gmap_vert == NULL) {
	g->L2Gmap_vert = phgAlloc(g->nvert * sizeof(*(g->L2Gmap_vert)));
	for (i = 0; i < g->nvert; i++)
	    g->L2Gmap_vert[i] = i;
    }

    /* pack tables */
    n = i = 0;
    while (TRUE) {
	while (i < g->nvert && iflags[i] == -1)
	    i++;
	if (i >= g->nvert)
	    break;
	memcpy(g->verts + n, g->verts + i, sizeof(*(g->verts)));
	if (g->period != NULL)
	    g->period->L2Gmap_vert[n] = g->period->L2Gmap_vert[i];
	g->L2Gmap_vert[n++] = g->L2Gmap_vert[i];
	i++;
    }

    phgFree(iflags);

    if (n == 0) {
	phgFree(g->verts);
	g->verts = NULL;
	phgFree(g->L2Gmap_vert);
	g->L2Gmap_vert = NULL;
	if (g->period != NULL) {
	    phgFree(g->period->L2Gmap_vert);
	    g->period->L2Gmap_vert = NULL;
	}
    }
    else {
	g->verts = phgRealloc_(g->verts, n * sizeof(*(g->verts)),
				g->nvert * sizeof(*(g->verts)));
	g->L2Gmap_vert = phgRealloc_(g->L2Gmap_vert,
				n * sizeof(*(g->L2Gmap_vert)),
				g->nvert * sizeof(*(g->L2Gmap_vert)));
	if (g->period != NULL)
	    g->period->L2Gmap_vert = phgRealloc_(g->period->L2Gmap_vert,
				n * sizeof(*(g->period->L2Gmap_vert)),
				g->nvert * sizeof(*(g->period->L2Gmap_vert)));
    }
    g->nvert = n;
}

static void
bcast_bdry_funcs(GRID *g)
{
#if ALLOW_CURVED_BOUNDARY
    EXPR **funcs;
    int n, m, len[2];
    char *buffer = NULL, *p, *q;

    /* broadcast functions for curved boundaries */
    if (g->bdry_funcs != NULL) {
	if (g->rank == 0) {
	    /* functions => strings */
	    n = 0;
	    funcs = g->bdry_funcs;
	    while (*funcs != NULL) {
		m = strlen(q = phgDump3DFunction(*(funcs++))) + 1;
		buffer = phgRealloc_(buffer, n + m, n);
		memcpy(buffer + n, q, m);
		n += m;
	    }
	    len[0] = n;				/* buffer size */
	    len[1] = funcs - g->bdry_funcs;	/* number of functions */
	    assert((len[1] % 4) == 0);
	}
	MPI_Bcast(len, 2, MPI_INT, 0, g->comm);
	n = len[0];
	if (n > 0) {
	    if (g->rank > 0) {
		buffer = phgAlloc(n);
		g->bdry_funcs = phgAlloc((len[1] + 1) * sizeof(*funcs));
	    }
            MPI_Bcast(buffer, n, MPI_BYTE, 0, g->comm);
	    if (g->rank > 0) {
		p = buffer;
		funcs = g->bdry_funcs;
		m = len[1];
	        while (m-- > 0) {
		    *(funcs++) = phgDefine3DFunction(p);
	            p += strlen(p) + 1;
		}
		*funcs = NULL;
	    }
            phgFree(buffer);
	}
    }
#else
    Unused(g);
#endif
}

/*--------------------------------------------------------------------------*/

static int rank_;
static char *cflags;

static BOOLEAN
redistribute_mark_branch_cb CB_ARGS (e)
{
    int i;

    if (rank_ == g_->rank) {
	/* mark elements to be removed */
	if (IsLeaf(e))
	    return TRUE;

	if ((e->children[0] != NULL && e->children[0]->mark == g_->rank) ||
	    (e->children[1] != NULL && e->children[1]->mark == g_->rank))
	    e->mark = g_->rank;
	else
	    e->mark = -1;

	return TRUE;
    }

    e->peer = NULL;		/* needed later */
    if (!IsLeaf(e)) {
	e->mark = ((e->children[0] != NULL && e->children[0]->mark == rank_) ||
		   (e->children[1] != NULL && e->children[1]->mark == rank_)) ?
		  rank_ : -1;
    }
    if (e->mark == rank_) {
	/* mark vertices to be migrated */
	for (i = 0; i < NVert; i++)
	    cflags[e->verts[i]] = 1;
	pack_count++;
    }

    return TRUE;
}

static void
redistribute_copy_branch(ELEMENT *e)
{
    static int i;
    ELEMENT *dest = pack_base++;

    *dest = *e;
    e->peer = dest;

    /* set up global indices */
    for (i = 0; i < NVert; i++)
	dest->verts[i] = GlobalVertex(g_, dest->verts[i]);

    if (e->children[0] != NULL && e->children[0]->mark == rank_)
	redistribute_copy_branch(e->children[0]);

    if (e->children[1] != NULL && e->children[1]->mark == rank_)
	redistribute_copy_branch(e->children[1]);
}

static BOOLEAN
redistribute_remove_neighbours CB_ARGS(e)
{
    int i;
    ELEMENT *e1;

    for (i = 0; i < NFace; i++) {
	if (e->bound_type[i] & INTERIOR) {
	    if ((e1 = e->neighbours[i]) == NULL || e1->mark != g_->rank)
		e->neighbours[i] = NULL;
	}
    }
    return TRUE;
}

static BOOLEAN
redistribute_update_neighbours1 CB_ARGS(e)
{
    int i;

    if (e->mark != g_->rank)
	return TRUE;

    for (i = 0; i < NFace; i++) {
	ELEMENT *e1 = e->neighbours[i];
	if ((e->bound_type[i] & INTERIOR) && e1 != NULL && e1->generation == 0)
	    e->neighbours[i] = e1->peer;
    }

    return TRUE;
}

static BOOLEAN
redistribute_remove_branch(ELEMENT *e)
/* returns TRUE if the branch is removed, FALSE otherwise */
{
    if (e->children[0] != NULL && redistribute_remove_branch(e->children[0]))
	e->children[0] = NULL;
    if (e->children[1] != NULL && redistribute_remove_branch(e->children[1]))
	e->children[1] = NULL;
    if (e->mark == g_->rank)
	return FALSE;
    phgFreeElement(&e);

    return TRUE;
}

static INT nleaf;
static BOOLEAN
redistribute_count_elements CB_ARGS(e)
/* count number of leaf and total elements */
{
    pack_count++;	/* ntree */
    if (IsLeaf(e))
	nleaf++;	/* nleaf */
    return TRUE;
}

typedef struct {
    COORD	coord;
    INT		gindex;
    INT		pindex;		/* periodic L2Gmap_vert[] entry */
    INT		lindex;
} VERT_t;

typedef struct {
    ELEMENT 	*e;		/* pointer to the element */
    INT		verts[NVert];	/* sorted indices of the vertices of the face,
				   the last entry is the opposite vertex */
    int		vertex;		/* the opposite vertex */
} FACE_t;
static FACE_t *faces = NULL;
size_t faces_count = 0, faces_allocated = 0;

typedef struct {
    ELEMENT	*e;		/* pointer to the element */
    INT		verts[NVert];	/* the sorted list of vertices */
} ELEMENT_t;
static ELEMENT_t *simplices = NULL;
static size_t simplices_allocated = 0, simplices_count = 0;

static int
redistribute_comp_simplex(const void *p0, const void *p1)
{
    INT *verts0 = ((ELEMENT_t *)p0)->verts, *verts1 = ((ELEMENT_t *)p1)->verts;
    INT j;
    int i;

    for (i = 0; i < NVert; i++)
	if ((j = verts0[i] - verts1[i]) != 0)
	    return j > 0 ? 1 : -1;

    return 0;
}

static BOOLEAN
is_remote(ELEMENT *e)
/* checks whether 'e' and its descendants have a REMOTE face */
{
    if (e == NULL)
	return FALSE;

    if (IsLeaf(e))
	return  ((e->bound_type[0] & INTERIOR) && e->neighbours[0] == NULL) ||
		((e->bound_type[1] & INTERIOR) && e->neighbours[1] == NULL) ||
		((e->bound_type[2] & INTERIOR) && e->neighbours[2] == NULL) ||
		((e->bound_type[3] & INTERIOR) && e->neighbours[3] == NULL);
    else
	return is_remote(e->children[0]) || is_remote(e->children[1]);
}

static BOOLEAN
redistribute_add_simplex CB_ARGS(e)
/* point 'peer' to the element itself, and add root elements to 'simplices' */
{
    ELEMENT_t *s;

    e->peer = e; /* for checking, see asserts after calling this func. */

    if (e->generation != 0)
	return TRUE;

    if (!is_remote(e))
	return TRUE;

    if (simplices_count >= simplices_allocated) {
	simplices = phgRealloc_(simplices,
			(simplices_allocated + 65536) * sizeof(*simplices),
			simplices_allocated * sizeof(*simplices));
	simplices_allocated += 65536;
    }

    s = simplices + (simplices_count++);
    s->e = e;
    memcpy(s->verts, e->verts, sizeof(e->verts));

    qsort(s->verts, NVert, sizeof(s->verts[0]), phgCompINT);

    return TRUE;
}

static void
redistribute_add_branch(ELEMENT *e0, ELEMENT *e)
{
    ELEMENT *a0, *a1, *b0, *b1;

    e0->peer = e->peer = e0;

    a0 = e0->children[0];
    a1 = e0->children[1];
    b0 = e->children[0];
    b1 = e->children[1];

    if (b0 != NULL) {
	if (a0 == NULL) {
	    a0 = e0->children[0] = phgNewElements(1);
	    *a0 = *b0;
	    a0->children[0] = a0->children[1] = NULL;
	}
	redistribute_add_branch(a0, b0);
    }
    if (b1 != NULL) {
	if (a1 == NULL) {
	    a1 = e0->children[1] = phgNewElements(1);
	    *a1 = *b1;
	    a1->children[0] = a1->children[1] = NULL;
	}
	redistribute_add_branch(a1, b1);
    }

    /* add neighbour links when children are both leaf or non leaf elements */
    if (a0 != NULL && a1 != NULL) {
	if (IsLeaf(a0) == IsLeaf(a1)) {
	    switch (e0->type) {
		case OPPOSITE:
		    a0->bound_type[2] &= ~REMOTE;
		    a1->bound_type[2] &= ~REMOTE;
		    a0->neighbours[2] = a1;
		    a1->neighbours[2] = a0;
		    break;
		case MIXED:
		    a0->bound_type[0] &= ~REMOTE;
		    a1->bound_type[2] &= ~REMOTE;
		    a0->neighbours[0] = a1;
		    a1->neighbours[2] = a0;
		    break;
		default:
		    a0->bound_type[0] &= ~REMOTE;
		    a1->bound_type[0] &= ~REMOTE;
		    a0->neighbours[0] = a1;
		    a1->neighbours[0] = a0;
	    }
	}
    }
}

static inline void
redistribute_add_face(ELEMENT *e, int vertex)
{
    INT i;
    FACE_t *f;

    if (faces_count >= faces_allocated) {
	faces = phgRealloc_(faces, (faces_allocated + 65536) * sizeof(*faces),
				   faces_allocated * sizeof(*faces));
	faces_allocated += 65536;
    }

    f = faces + (faces_count++);
    f->e = e;
    f->vertex = vertex;
    if (g_->period != NULL) {
	f->verts[0] = g_->period->L2Gmap_vert[e->verts[0]];
	f->verts[1] = g_->period->L2Gmap_vert[e->verts[1]];
	f->verts[2] = g_->period->L2Gmap_vert[e->verts[2]];
#if Dim == 3
	f->verts[3] = g_->period->L2Gmap_vert[e->verts[3]];
#endif	/* Dim == 3 */
    }
    else {
	f->verts[0] = e->verts[0];
	f->verts[1] = e->verts[1];
	f->verts[2] = e->verts[2];
#if Dim == 3
	f->verts[3] = e->verts[3];
#endif	/* Dim == 3 */
    }
    if (vertex != NVert - 1) {
	i = f->verts[vertex];
	f->verts[vertex] = f->verts[NVert - 1];
	f->verts[NVert - 1] = i;
    }

    qsort(f->verts, NVert - 1, sizeof(f->verts[0]), phgCompINT);
}

static BOOLEAN
redistribute_collect_faces CB_ARGS(e)
{
    int i;

    for (i = 0; i < NFace; i++) {
	if ((e->bound_type[i] & INTERIOR) && e->neighbours[i] == NULL) {
	    if (faces == NULL)
		faces_count++;
	    else
		redistribute_add_face(e, i);
	}
    }

    return TRUE;
}

static int
redistribute_comp_face0(const void *p0, const void *p1)
{
    INT *verts0 = faces[*(INT *)p0].verts, *verts1 = faces[*(INT *)p1].verts;
    INT j;
    int i;

    for (i = 0; i < NVert - 1; i++)
	if ((j = verts0[i] - verts1[i]) != 0)
	    return j > 0 ? 1 : -1;

    return 0;
}

static int
redistribute_comp_face(const void *p0, const void *p1)
{
    int i;

    if ((i = redistribute_comp_face0(p0, p1)) != 0)
	return i;

    /* this makes simplex with smaller generation come first */
    return faces[*(INT *)p0].e->generation - faces[*(INT *)p1].e->generation;
}

static VERT_t *rverts;

static int
redistribute_comp_gindex0(const void *p0, const void *p1)
{
    INT i;
    return (i = rverts[*(INT *)p0].gindex - rverts[*(INT *)p1].gindex) > 0 ?
		1 : (i < 0 ? -1 : 0);
}

static int
redistribute_comp_gindex(const void *p0, const void *p1)
{
    INT i;
    return (i = ((VERT_t *)p0)->gindex - ((VERT_t *)p1)->gindex) > 0 ?
		1 : (i < 0 ? -1 : 0);
}

static BOOLEAN
redistribute_update_neighbours2 CB_ARGS(e)
/* update neighbours because g->roots has changed (old address in pack_base)
   TODO: can be marged with updating elements from rbranches later */
{
    int i;
    ELEMENT *e1;

    for (i = 0; i < NFace; i++) {
	if (!(e->bound_type[i] & INTERIOR))
	    continue;
	if ((void *)(e1 = e->neighbours[i]) < (void *)pack_base)
	    continue;
	if ((void *)e1 >= (void *)(pack_base + g_->nroot))
	    continue;
	e->neighbours[i] = (void *)(((char *)e1) + pack_offset);
    }

    return TRUE;
}

static BOOLEAN
redistribute_set_global_indices CB_ARGS(e)
{
    int i;

    /* change edge indices to global */
    if ((g_->flags & EDGE_FLAG)) {
	for (i = 0; i < NEdge; i++) {
	    if (e->edges[i] != -1)
		e->edges[i] = GlobalEdge(g_, e->edges[i]);
	}
    }

#if Dim == 3
    /* change face indices to global */
    if ((g_->flags & FACE_FLAG)) {
	for (i = 0; i < NFace; i++) {
	    if (e->faces[i] != -1)
		e->faces[i] = GlobalFace(g_, e->faces[i]);
	}
    }
#endif

    if ((g_->flags & ELEM_FLAG)) {
	e->index = GlobalElement(g_, e->index);
    }

    /* clear REMOTE flags */
    if (e->bound_type[0] & REMOTE) {
	e->bound_type[0] &= ~REMOTE;
	e->neighbours[0] = NULL;
    }
    if (e->bound_type[1] & REMOTE) {
	e->bound_type[1] &= ~REMOTE;
	e->neighbours[1] = NULL;
    }
    if (e->bound_type[2] & REMOTE) {
	e->bound_type[2] &= ~REMOTE;
	e->neighbours[2] = NULL;
    }
    if (e->bound_type[3] & REMOTE) {
	e->bound_type[3] &= ~REMOTE;
	e->neighbours[3] = NULL;
    }

    return TRUE;
}

/* local variables used by phgRedistributeGrid_ and phgBalanceGrid_ */
static int nprocs_requested = 0;
static BOOLEAN grid_changed = FALSE;
static BOOLEAN balancing_grid = FALSE;

static void
hp_delete_derived(HP_TYPE *hp)
/* recursively delete derived HP_TYPEs which are not referenced */
{
    if (hp->derivative != NULL) {
	/* derived HP_TYPE must come after */
	assert(hp->id < hp->derivative->id);
	hp_delete_derived(hp->derivative);
	phgInfo(2, "delete derived HP_TYPE created at %s:%d\n",
		   hp->derivative->srcfile, hp->derivative->srcline);
	if (hp->derivative->refcount == 0)
	    phgHPFree(&hp->derivative);
    }

    if (hp->dg != NULL) {
	/* derived HP_TYPE must come after */
	assert(hp->id < hp->dg->id);
	hp_delete_derived(hp->dg);
	phgInfo(2, "delete derived HP_TYPE created at %s:%d\n",
		   hp->dg->srcfile, hp->dg->srcline);
	if (hp->dg->refcount == 0)
	    phgHPFree(&hp->dg);
    }
}

void
phgRedistributeGrid_(GRID *g, MPI_Comm comm_world, PARTITIONER part_func)
{
    /* k>1 ==> k-ary tree-wise multi-stage mesh redistribution.
     * This option is effective for direct phgRedistributeGrid_ calls only,
     * not for phgBalanceGrid_. */
    static INT partition_k = 8;
    static BOOLEAN partition_k_flag = FALSE;	/* work variable */

    INT i, j, total, total_verts, *order;
    int ii, jj, *counts, *displs, *nverts, *rcounts, *rdispls, *rnverts;
    int nsend, old_nprocs;
    ELEMENT *e, *e0, *branches, *rbranches;
    VERT_t *verts = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    MPI_Datatype type;
    static BOOLEAN initialized = FALSE;
    double t0, t;	/* for timing */
    FLOAT old_lif;
    int nprocs_max, rank;
    int part_id;		/* partitioner id */
    const char *part_name;	/* partitioner name */

    FunctionEntry_head;

    if (!initialized) {
	/* Register options */
	initialized = TRUE;
	if (partitioners_count <= 0)
	    return;
	phgOptionsRegisterTitle("\nMesh partitioning options:", "\n",
								"partition");
	phgOptionsRegisterKeyword("-partitioner", "Mesh partitioner",
		partitioners_name, &partitioner);
	phgOptionsRegisterInt("-partition_k",
		"K for k-ary tree-wise multi-stage initial mesh distribution "
		"in phgImport(.,TRUE)", &partition_k);
	/* register partitioner specific options */
	for (ii = 0; ii < partitioners_count; ii++)
	    partitioners_list[ii](NULL, 0, NULL, 0.);

        /* register specific options for 1D partitioner */
        phgPartition1DOptionsSet();
	return;
    }

    FunctionEntry_tail;

    /* prepare partitioner */
    if (part_func != NULL) {
	/* check whether part_func is a predefined one */
	for (ii = 0; ii < partitioners_count; ii++)
	    if (part_func == partitioners_list[ii])
		break;
	if (ii < partitioners_count) {
	    part_id = ii;
	    part_name = partitioners_name[ii];
	}
	else {
	    part_id = -1;
	    part_name = "custom";
	}
    }
    else {
	assert(partitioner >= 0 && partitioner < partitioners_count);
	part_id = partitioner;
	part_func = partitioners_list[partitioner];
	part_name = partitioners_name[partitioner];
    }

    grid_changed = FALSE;

    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &nprocs_max);

    if (nprocs_max <= 1)
	Return;

    old_nprocs = g->nprocs;
    if (nprocs_requested <= 0 || nprocs_requested > nprocs_max)
	nprocs_requested = nprocs_max;
    if (nprocs_requested < g->nprocs)
	nprocs_requested = g->nprocs;

    /* Use a multi-stage tree-wise mesh distribution algorithm,
     * if not called by phgBalanceGrid_ and partition_k>1 */
    if (!balancing_grid && partition_k > 1 && partition_k_flag == FALSE) {
	int nprocs0 = nprocs_requested;
	partition_k_flag = TRUE;
	while (g->nprocs < nprocs0) {
	    double t0 = phgGetTime(NULL);
	    int nprocs_old = g->nprocs;
	    nprocs_requested = g->nprocs * partition_k;
	    if (nprocs_requested > nprocs0)
		nprocs_requested = nprocs0;
	    phgRedistributeGrid_(g, comm_world, part_func);
	    if (phgVerbosity > 0)
		phgPrintf("***** nprocs %d -> %d: %0.2lfs\n",
				nprocs_old, g->nprocs, phgGetTime(NULL) - t0);
	    if (nprocs_old == g->nprocs)
		break;
	}
	partition_k_flag = FALSE;
	Return;
    }

    g_ = g;

    phgDebug((3, "1: g->nprocs = %d, nprocs_requested = %d\n",
		g->nprocs, nprocs_requested));

    if (!balancing_grid && !partition_k_flag) {
#if NO_NEW_COMMUNICATOR
	comm = comm_world;
#else	/* NO_NEW_COMMUNICATOR */
	MPI_Comm_dup(comm_world, &comm);
#endif	/* NO_NEW_COMMUNICATOR */
    }
    else if (g->nprocs == nprocs_requested) {
	comm = g->comm;
    }
    else {
	nprocs_requested = phgCreateComm_(nprocs_requested, comm_world, &comm);
    }

    phgDebug((3, "2: g->nprocs = %d, nprocs_requested = %d, comm_null = %d\n",
		g->nprocs, nprocs_requested, comm == MPI_COMM_NULL));

    if (comm == MPI_COMM_NULL)
	Return;

    if (nprocs_requested == 1) {
#if !NO_NEW_COMMUNICATOR
	if (comm != g->comm)
	    MPI_Comm_free(&comm);
#endif	/* !NO_NEW_COMMUNICATOR */
	comm = MPI_COMM_NULL;
	Return;
    }
    if (phgVerbosity > 0)
	phgPrintf("Partitioner: %s\n", part_name);
    if (g->last_partitioner >= 0 && part_id != g->last_partitioner
	&& part_func == phgPartitionRTK) {
	static BOOLEAN warned = FALSE;
	if (!warned) {
	    warned = TRUE;
	    phgPrintf("Warning: RTK used after another partitioner.\n");
	}
    }
    old_lif = g->lif;		/* save old LIF */
    g->lif = -1.0;		/* this allows us to detect whether the
				 * partitioner updates the LIF */
    t0 = phgGetTime(NULL);
    grid_changed = part_func(g, comm, weights, power);
    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for partitioning (%s): %0.4lf\n", __func__,
			part_name, t0);
    }

    g->last_partitioner = part_id;

    if (g->nprocs != nprocs_requested) {
	/* bring the grid on newly joined processes up to date */
	int nhp, nhp_min;
	int size = sizeof(grid_changed) + 7 * sizeof(FLOAT) + 7 * sizeof(INT) +
		   sizeof(g->serial_no) + sizeof(g->bdry_funcs) + 
		   sizeof(g->volume) + PATH_MAX;
	char *buf, fn[PATH_MAX];
	HP_TYPE **php;
#define PACK(v) MPI_Pack(&(v), sizeof(v), MPI_BYTE, buf, size, &nsend, comm)
#define UNPACK(v) MPI_Unpack(buf, size, &nsend, &(v), sizeof(v), MPI_BYTE, comm)
	/* delete derived HP_TYPEs which are not currently used,
	 * and count total # of HP_TYPEs */
	for (nhp = 0, php = g->hp; php != NULL && *php != NULL; php++)
	    (*php)->id = nhp++;
	if (nhp > 0) {
	    for (php = g->hp; php != NULL && *php != NULL; php++)
		hp_delete_derived(*php);
	    nhp = (int)(php - g->hp);
	}
#if (1 || DEBUG)
	/* check consistency of HP_TYPEs across processes */
	MPI_Reduce(&nhp, &nhp_min, 1, MPI_INT, MPI_MIN, 0, comm);
	if (rank == 0 && nhp_min != nhp) {
	    phgInfo(-1, "Fatal error: HP_TYPEs mismatch across processes.\n");
	    phgInfo(-1, "The following HP_TYPEs were (possibly) not "
			"properly freed:\n");
	    for (; nhp_min < nhp; nhp_min++)
		phgInfo(-1, "\t%s (created at %s:%d)\n",
			g->hp[nhp_min]->info->name, 
			g->hp[nhp_min]->srcfile, g->hp[nhp_min]->srcline);
	    phgInfo(-1, "Note: an HP_TYPE created after calling "
			"phgRedistributeGrid_ must be freed before the next "
			"call to phgRedistributeGrid_.\n");
	    phgInfo(-1, "Abort.\n");
	}
#endif	/* DEBUG */
	MPI_Bcast(&g->bc_n, 1, MPI_INT, 0, comm);
	size += nhp * (2 + 4 * sizeof(INT)) + 2 * g->bc_n * sizeof(*g->bc_list);
	buf = phgAlloc(size);
	nsend = 0;
	if (rank == 0) {
	    PACK(grid_changed);
	    if (grid_changed || !balancing_grid) {
		PACK(g->lif);
		PACK(g->nvert_global);
		PACK(g->nedge_global);
		PACK(g->nface_global);
		PACK(g->nelem_global);
		PACK(g->serial_no);
		PACK(g->bbox);
		PACK(g->volume);
		if (g->bc_n > 0) {
		    MPI_Pack(g->bc_list, g->bc_n * sizeof(*g->bc_list),
				MPI_BYTE, buf, size, &nsend, comm);
		    MPI_Pack(g->bc_rmap, g->bc_n * sizeof(*g->bc_rmap),
				MPI_BYTE, buf, size, &nsend, comm);
		}
		PACK(g->bdry_funcs);
		if (g->period != NULL) {
		    i = g->period->flags;
		    PACK(i);
		    PACK(g->period->nvert_global);
		}
		else {
		    i = 0;
		    PACK(i);
		}
		/* pack HP_TYPE entries */
		for (php = g->hp; php != NULL && *php != NULL; php++) {
		    PACK((*php)->min_order);
		    PACK((*php)->max_order);
		    PACK((*php)->vert_data_count_global);
		    PACK((*php)->edge_data_count_global);
		    PACK((*php)->face_data_count_global);
		    PACK((*php)->elem_data_count_global);
		}
		strcpy(fn, g->filename == NULL ? "" : g->filename);
		PACK(*fn);
	    }
	    MPI_Bcast(buf, size, MPI_PACKED, 0, comm);
	}
	else {
	    MPI_Bcast(buf, size, MPI_PACKED, 0, comm);
	    UNPACK(grid_changed);
	    if (grid_changed || !balancing_grid) {
		UNPACK(g->lif);
		UNPACK(g->nvert_global);
		UNPACK(g->nedge_global);
		UNPACK(g->nface_global);
		UNPACK(g->nelem_global);
		UNPACK(g->serial_no);
		UNPACK(g->bbox);
		UNPACK(g->volume);

		phgFree(g->bc_list);
		phgFree(g->bc_rmap);
		g->bc_alloc = g->bc_n;
		if (g->bc_n > 0) {
		    g->bc_list = phgAlloc(g->bc_n * sizeof(*g->bc_list));
		    MPI_Unpack(buf, size, &nsend, g->bc_list,
				g->bc_n * sizeof(*g->bc_list), MPI_BYTE, comm);
		    g->bc_rmap = phgAlloc(g->bc_n * sizeof(*g->bc_rmap));
		    MPI_Unpack(buf, size, &nsend, g->bc_rmap,
				g->bc_n * sizeof(*g->bc_rmap), MPI_BYTE, comm);
		}
		else {
		    g->bc_list = NULL;
		    g->bc_rmap = NULL;
		}

		/* Note: bdry_funcs will be updated later */
		UNPACK(g->bdry_funcs);
		UNPACK(i);
		if (i) {
		    if (g->period == NULL)
			g->period = phgCalloc(1, sizeof(*g->period));
		    g->period->flags = i;
		    UNPACK(g->period->nvert_global);
		}
		/* unpack HP_TYPE entries */
		for (php = g->hp; php != NULL && *php != NULL; php++) {
		    UNPACK((*php)->min_order);
		    UNPACK((*php)->max_order);
		    UNPACK((*php)->vert_data_count_global);
		    UNPACK((*php)->edge_data_count_global);
		    UNPACK((*php)->face_data_count_global);
		    UNPACK((*php)->elem_data_count_global);
		}
		UNPACK(*fn);
		phgFree(g->filename);
		g->filename = (fn[0] == '\0' ? NULL : strdup(fn));
	    }
	}
	phgFree(buf);

	if (!grid_changed && balancing_grid) {
#if !NO_NEW_COMMUNICATOR
	    if (comm != g->comm)
		MPI_Comm_free(&comm);
#endif	/* !NO_NEW_COMMUNICATOR */
	    comm = MPI_COMM_NULL;
	    g->lif = old_lif;		/* restore LIF */
	    Return;
	}

#if !NO_NEW_COMMUNICATOR
	if (comm != g->comm)
	    MPI_Comm_free(&g->comm);
#endif	/* !NO_NEW_COMMUNICATOR */
	g->comm = comm;
	comm = MPI_COMM_NULL;
	MPI_Comm_size(g->comm, &g->nprocs);
	MPI_Comm_rank(g->comm, &g->rank);

	if (g->neighbours.counts == NULL) {
	    g->neighbours.counts = phgCalloc(g->nprocs,
					     sizeof(*g->neighbours.counts));
	    g->neighbours.displs = phgCalloc(g->nprocs,
					     sizeof(*g->neighbours.displs));
	}
	else if (old_nprocs < g->nprocs) {
	    g->neighbours.counts = phgRealloc_(g->neighbours.counts,
						g->nprocs *
						sizeof(*g->neighbours.counts),
						old_nprocs *
						sizeof(*g->neighbours.counts));
	    g->neighbours.displs = phgRealloc_(g->neighbours.displs,
						g->nprocs *
						sizeof(*g->neighbours.displs),
						old_nprocs *
						sizeof(*g->neighbours.displs));
	    for (ii = old_nprocs; ii < g->nprocs; ii++) {
		g->neighbours.counts[ii] = 0;
		g->neighbours.displs[ii] = g->neighbours.displs[ii - 1];
	    }
	}

	bcast_bdry_funcs(g);
    }
    else {
#if !NO_NEW_COMMUNICATOR
	if (comm != g->comm)
	    MPI_Comm_free(&comm);
#endif	/* !NO_NEW_COMMUNICATOR */
	comm = MPI_COMM_NULL;
	if (!grid_changed) {
	    g->lif = old_lif;		/* restore LIF */
	    Return;
	}
    }

    g->alien = NULL;
    phgFree(g->alien_map);
    g->alien_map = NULL;
    phgFree(g->Fmap_elems);
    g->Fmap_elems = NULL;
    phgFree(g->Fmap_faces);
    g->Fmap_faces = NULL;

    if (g->halo != NULL) {
	assert(g->halo->type >= 0);
	phgRemoveHalo(g);
	g->serial_no -= 1000 * 1000;	/* restore serial_no */
    }
    g->serial_no += 1000;

    phgDofRedistribute(g, FALSE);	/* save state of current mesh */

    t0 = phgGetTime(NULL);

    /* change edge/face/element indices to global and free L2Gmaps */
    phgTraverseAllElements(g, redistribute_set_global_indices);
    phgFree(g->L2Gmap_edge);
    g->L2Gmap_edge = NULL;
#if Dim == 3
    phgFree(g->L2Gmap_face);
    g->L2Gmap_face = NULL;
#endif
    phgFree(g->L2Gmap_face);
    g->L2Gmap_face = NULL;

    /* clear list of remote neighbours */
    if (g->neighbours.list != NULL) {
	phgFree(g->neighbours.list);
	g->neighbours.list = NULL;
	g->neighbours.allocated = g->neighbours.count = 0;
    }

    counts = phgAlloc(6 * g->nprocs * sizeof(*counts));
    nverts = counts + g->nprocs;
    displs = nverts + g->nprocs;
    rcounts = displs + g->nprocs;
    rnverts = rcounts + g->nprocs;
    rdispls = rnverts + g->nprocs;
    cflags = phgAlloc(g->nvert);

    /* construct individual tree for each processor */
    total_verts = total = 0;
    branches = NULL;
    /* Initialize counts[] such that:
     * 	counts[rank] = number of elements for proc rank
     * It is used to avoid mesh traversals for processes to which no elements
     * are to be sent */
    memset(counts, 0, g->nprocs * sizeof(*counts));
    ForAllElements(g, e) {
	assert(e->mark >= 0 && e->mark < g->nprocs && counts[e->mark]<INT_MAX);
	counts[e->mark]++;
    }
    nsend = 0;		/* nsend = # of procs to which elements are sent */
    for (rank_ = 0; rank_ < g->nprocs; rank_++) {
	pack_count = 0;
	if (rank_ != g->rank && counts[rank_] > 0) {
	    nsend++;
	    if (g->nvert > 0)
		memset(cflags, 0, g->nvert);
	    phgTraverseAllElements(g, redistribute_mark_branch_cb);
	}

	counts[rank_] = pack_count; assert(pack_count <= INT_MAX);
	if (pack_count <= 0) {
	    nverts[rank_] = 0;
	    continue;
	}

	branches = phgRealloc_(branches,
			(total + pack_count) * sizeof(*branches),
			total * sizeof(*branches));
	total += pack_count;

	/* copy elements */
	pack_base = branches + total - pack_count;
#if DEBUG
	j = 0;
#endif
	for (i = 0; i < g->nroot; i++) {
	    if (g->roots[i].mark == rank_) {
		redistribute_copy_branch(g->roots + i);
#if DEBUG
		j++;
#endif
	    }
	}

	phgDebug((3, "# of elements for proc %d: %d (%d roots)\n",
		rank_, pack_count, j));

	/* set up links in the copied branch */
	pack_base -= pack_count;
	for (e = pack_base; e < pack_base + pack_count; e++) {
	    /* pointers to children */
	    if ((e0 = e->children[0]) != NULL && e0->mark == rank_) {
		e->children[0] = (void *)((e0 = e0->peer) - pack_base);
		assert(e0 >= pack_base && e0 < pack_base + pack_count);
	    }
	    else {
		e->children[0] = (void *)-1;
	    }
	    if ((e0 = e->children[1]) != NULL && e0->mark == rank_) {
		e->children[1] = (void *)((e0 = e0->peer) - pack_base);
		assert(e0 >= pack_base && e0 < pack_base + pack_count);
	    }
	    else {
		e->children[1] = (void *)-1;
	    }
	    /* pointers to neighbours */
	    if (!IsLeaf(e)) {
		for (ii = 0; ii < NFace; ii++)
		    e->neighbours[ii] = (void *)-1;
		continue;
	    }
	    for (ii = 0; ii < NFace; ii++) {
		if (e->bound_type[ii] & INTERIOR) {
		    if ((e0 = e->neighbours[ii]) != NULL && e0->mark == rank_) {
			e0 = e0->peer;
#if DEBUG
			if (e0 < pack_base || e0 >= pack_base + pack_count)
			    phgError(1, "%s:%d, unexpected error!\n",
				__FILE__, __LINE__);
#endif
			e->neighbours[ii] = (void *)(e0 - pack_base);
		    }
		    else {
			e->neighbours[ii] = (void *)-1;
		    }
		}
		else {
		    e->neighbours[ii] = (void *)-1;
		}
	    }
	}

	/* copy vertices */
	pack_count = 0;
	for (i = 0; i < g->nvert; i++)
	    pack_count += cflags[i];
	nverts[rank_] = pack_count;	assert(pack_count <= INT_MAX);
	verts = phgRealloc_(verts, (total_verts + pack_count) * sizeof(*verts),
				   total_verts * sizeof(*verts));
	for (i = 0; i < g->nvert; i++) {
	    if (cflags[i] == 0)
		continue;
	    rverts = verts + (total_verts++);
	    memcpy(rverts->coord, g->verts[i], sizeof(COORD));
	    rverts->gindex = GlobalVertex(g, i);
	    if (g->period != NULL)
		rverts->pindex = g->period->L2Gmap_vert[i];
	    rverts->lindex = -1;
	}
    }

    phgFree(cflags);

    /* exchange numbers of elements and vertices to migrate */
#if 0
    {
	MPI_Datatype types[2];
	MPI_Aint indices[2];
	int blocklens[2];

	MPI_Type_vector(2, 1, nverts - counts, MPI_INT, &types[0]);
	blocklens[0] = blocklens[1] = 1;
	indices[0] = 0;
	MPI_Type_extent(MPI_INT, &indices[1]);
	types[1] = MPI_UB;
	MPI_Type_struct(2, blocklens, indices, types, &type);
	MPI_Type_commit(&type);
	MPI_Alltoall(counts, 1, type, rcounts, 1, type, g->comm);
	MPI_Type_free(&type);
	MPI_Type_free(&types[0]);
    }
#else
    /* This is a workaround for a bug of MPICH-1.2.5 */
    {
	int *tmp = phgAlloc(4 * g->nprocs * sizeof(int));
	int *tmp1 = tmp + 2 * g->nprocs;
	for (ii = 0; ii < g->nprocs; ii++) {
	    tmp[2 * ii + 0] = counts[ii];
	    tmp[2 * ii + 1] = nverts[ii];
	}
	MPI_Alltoall(tmp, 2, MPI_INT, tmp1, 2, MPI_INT, g->comm);
	for (ii = 0; ii < g->nprocs; ii++) {
	    rcounts[ii] = tmp1[2 * ii + 0];
	    rnverts[ii] = tmp1[2 * ii + 1];
	}
	phgFree(tmp);
    }
#endif

    /* exchange simplices */
    j = total = 0;
    for (ii = 0; ii < g->nprocs; ii++) {
	displs[ii] = j;
	rdispls[ii] = total;
	j += counts[ii];		assert(j <= INT_MAX);
	total += rcounts[ii];		assert(total <= INT_MAX);
    }
    rbranches = (total > 0) ? phgAlloc(total * sizeof(*rbranches)) : NULL;
    MPI_Type_contiguous(sizeof(*branches), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(branches, counts, displs, type,
		  rbranches, rcounts, rdispls, type, g->comm);
    MPI_Type_free(&type);

    if (phgVerbosity > 0) {
#if DEBUG
	int nsend_max;
	MPI_Reduce(&total, &i, 1, PHG_MPI_INT, MPI_SUM, 0, g->comm);
	MPI_Reduce(&total, &j, 1, PHG_MPI_INT, MPI_MAX, 0, g->comm);
	phgPrintf("%s - maxi number of elements sent: %d, total: %d\n",
			__func__, j, i);
	MPI_Reduce(&nsend, &nsend_max, 1, MPI_INT, MPI_MAX, 0, g->comm);
	phgPrintf("%s - maxi number of target processes: %d\n",
			__func__, nsend_max);
#endif	/* DEBUG */
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for migrating elements: %0.4lf\n", __func__, t0);
	t0 = phgGetTime(NULL);
    }

    /* restore links in the received elements */
    for (ii = 0; ii < g->nprocs; ii++) {
	pack_base = rbranches + rdispls[ii];
	pack_count = rcounts[ii];
	for (e = pack_base; e < pack_base + pack_count; e++) {
	    e->mark = ii;  /* mark where the elem comes from for debugging */
#if 0
	    e->user_ctx = NULL; /* reset user_ctx pointer */
#endif
	    /* restore address of the children */
	    if (e->children[0] != (void *)-1)
		e->children[0] = pack_base + (size_t)(e->children[0]);
	    else
		e->children[0] = NULL;

	    if (e->children[1] != (void *)-1)
		e->children[1] = pack_base + (size_t)(e->children[1]);
	    else
		e->children[1] = NULL;

	    /* restore address of the neighbours */
	    for (jj = 0; jj < NFace; jj++) {
		if (e->neighbours[jj] != (void *)-1)
		    e->neighbours[jj] = pack_base + (size_t)(e->neighbours[jj]);
		else
		    e->neighbours[jj] = NULL;
	    }
	}
    }

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for restoring links: %0.4lf\n", __func__, t0);
	t0 = phgGetTime(NULL);
    }

    phgDebug((3, "old nroot, nleaf, ntree: %d %d %d\n",
		g->nroot, g->nleaf, g->ntree));

    if (branches != NULL) { /* remove migrated branches */
	ELEMENT *p;

	phgFree(branches);

	/* first, mark branches to be removed */
	rank_ = g->rank;
	phgTraverseAllElements(g, redistribute_mark_branch_cb);

	/* remove neighbours which have migrated */
	phgTraverseElements(g, redistribute_remove_neighbours);

	/* prune branches migrated */
	p = g->roots;
	for (i = 0; i < g->nroot; i++) {
	    ELEMENT **c;
	    e = g->roots + i;
	    c = e->children;
	    if (c[0] != NULL && redistribute_remove_branch(c[0]))
		c[0] = NULL;
	    if (c[1] != NULL && redistribute_remove_branch(c[1]))
		c[1] = NULL;
	    if (e->mark == g->rank)
		e->peer = p++;
	}

	if (p - g->roots < g->nroot) {
	    /* update neighbours */
	    phgTraverseElements(g, redistribute_update_neighbours1);

	    /* pack root elements */
	    p = g->roots;
	    for (i = 0; i < g->nroot; i++) {
		e = g->roots + i;
		if (e->mark == g->rank) {
		    if (p != e) *(p++) = *e; else p++;
		}
	    }
	    g->nroot = p - g->roots;
	    if (g->nroot == 0) {
		phgFree(g->roots);
		g->roots = NULL;
	    }
	}
    }

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for removing branches: %0.4lf\n", __func__, t0);
	t0 = phgGetTime(NULL);
    }

    pack_vertices(g);

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for packing vertices: %0.4lf\n", __func__, t0);
	t0 = phgGetTime(NULL);
    }

    /* update g->verts, g->L2Gmap_vert, and g->period->L2Gmap_vert */
    j = total_verts = 0;
    for (ii = 0; ii < g->nprocs; ii++) {
	displs[ii] = j;			assert(j <= INT_MAX);
	rdispls[ii] = total_verts;	assert(total_verts <= INT_MAX);
	j += nverts[ii];
	total_verts += rnverts[ii];
    }
    if (total_verts > 0)
	rverts = phgAlloc(total_verts * sizeof(*rverts));
    else
	rverts = NULL;
    /* Note: don't need to send lindex, and pindex if g->period == NULL */
    MPI_Type_contiguous(sizeof(*verts), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(verts,  nverts,  displs,  type,
		  rverts, rnverts, rdispls, type, g->comm);
    MPI_Type_free(&type);

    if (verts != NULL)
	phgFree(verts);
    phgFree(counts);

    if (total_verts > 0) {
	INT gindex;
	VERT_t v0;

	order = phgAlloc(total_verts * sizeof(*order));
	for (i = 0; i < total_verts; i++)
	    order[i] = i;
	qsort(order, total_verts, sizeof(*order), redistribute_comp_gindex0);

	/* copy sorted entries, removing duplicate vertices */
	j = i = 0;
	verts = phgAlloc(total_verts * sizeof(*verts));
	while (j < total_verts) {
	    gindex = (verts[i++] = rverts[order[j]]).gindex;
	    while (++j < total_verts && rverts[order[j]].gindex == gindex);
	}
	phgFree(order);
	phgFree(rverts);
	rverts = verts;
	j = total_verts = i;
	for (i = 0; i < g->nvert; i++) {
	    v0.gindex = GlobalVertex(g, i);
	    if ((verts = bsearch(&v0, rverts, total_verts, sizeof(*rverts),
			redistribute_comp_gindex)) != NULL) {
		if (verts->lindex != -1) {
		    phgInfo(-1, "duplicate entry: verts[%d]=verts[%d]=%d\n", i,
			verts->lindex, v0.gindex);
		    phgError(1, "%s:%d, abort.\n", __FILE__, __LINE__);
		}
		verts->lindex = i;
		j--;
	    }
	}

	/* Reallocate and update L2Gmap_vert and g->verts */
	if (j > 0) {
	    INT k = g->nvert;
	    g->verts = phgRealloc_(g->verts,
				(j + g->nvert) * sizeof(*(g->verts)),
				g->nvert * sizeof(*(g->verts)));
	    if (g->L2Gmap_vert == NULL) {
		g->L2Gmap_vert =
			phgAlloc((j + g->nvert) * sizeof(*(g->L2Gmap_vert)));
		for (i = 0; i < g->nvert; i++)
		    g->L2Gmap_vert[i] = i;
	    }
	    else {
		g->L2Gmap_vert = phgRealloc_(g->L2Gmap_vert,
				(j + g->nvert) * sizeof(*(g->L2Gmap_vert)),
				g->nvert * sizeof(*(g->L2Gmap_vert)));
	    }
	    if (g->period != NULL) {
		g->period->L2Gmap_vert = phgRealloc_(g->period->L2Gmap_vert,
			(j + g->nvert) * sizeof(*(g->period->L2Gmap_vert)),
			g->nvert * sizeof(*(g->period->L2Gmap_vert)));
		phgFree(g->period->ordering);
		g->period->ordering = NULL;	/* updated in utils.c */
	    }
	    g->nvert += j;
	    for (i = 0; i < total_verts; i++) {
		if (rverts[i].lindex != -1)
		    continue;
		rverts[i].lindex = k;
		memcpy(g->verts + k, rverts[i].coord, sizeof(g->verts[0]));
		if (g->period != NULL)
		    g->period->L2Gmap_vert[k] = rverts[i].pindex;
		g->L2Gmap_vert[k++] = rverts[i].gindex;
	    }
	    if (k != g->nvert)
		phgError(1, "%s:%d, unexpected error: k=%d, g->nvert=%d\n",
				__FILE__, __LINE__, k, g->nvert);
	}

	/* update verts[] in elements */
	for (i = 0; i < total; i++) {
	    for (jj = 0; jj < NVert; jj++) {
		v0.gindex = rbranches[i].verts[jj];
		if ((verts = bsearch(&v0, rverts, total_verts, sizeof(*rverts),
			redistribute_comp_gindex)) == NULL)
		    phgError(1, "%s:%d, error in bsearch, key = %"dFMT", "
				"size = %"dFMT".\n", __FILE__, __LINE__,
				v0.gindex, total_verts);
		rbranches[i].verts[jj] = verts->lindex;
	    }
	}
    }

    if (rverts != NULL)
	phgFree(rverts);

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for updating g->verts/g->L2Gmap: %0.4lf\n",
			__func__, t0);
	t0 = phgGetTime(NULL);
    }

    /* merge incoming branches */
    if (total > 0) {
	INT i0, i1, i2;
	FACE_t *f0, /**f1,*/ *f2;
	ELEMENT_t *s0, *s1, *s2;

	/* collect all root elements */
	for (e = g->roots; e < g->roots + g->nroot; e++)
	    redistribute_add_simplex(e);
	for (e = rbranches; e < rbranches + total; e++)
	    redistribute_add_simplex(e);

	/* sort root elements */
	qsort(simplices, simplices_count, sizeof(*simplices),
		redistribute_comp_simplex);

	/* select the one to keep among each set of duplicate root elements,
	   point the peer member of the elements to it,
	   and count number of new root elements.

	   For a set of identical elements, if one of them is in the local
	   tree, that one is selected. Otherwise the one with the smallest
	   index in 'rbranches' is selected.
	   */
	s0 = simplices;
	while (s0 < simplices + simplices_count) {
	    s1 = s2 = s0;
	    while (++s1 < simplices + simplices_count &&
			!redistribute_comp_simplex(s0, s1)) {
		if ((void *)s2->e < (void *)rbranches
			|| (void *)s2->e >= (void *)(rbranches + total))
		    continue;	/* s2 already points to a local element */
		if ((void *)s1->e < (void *)rbranches
			|| (void *)s1->e >= (void *)(rbranches + total))
		    s2 = s1;	/* select the element in the local tree */
		else if (s2->e > s1->e)
		    s2 = s1;
	    }
	    /* The elements between [s0,s1) are identical */
	    phgDebug((4, "identical root elements: %d\n", s1 - s0));
	    /* deselect elements other than s2 by pointing 'peer' to s2 */
	    while (s0 < s1) {
		s0->e->peer = s2->e;
		s0++;
	    }
	}

	/* count number of new root elements */
	pack_count = 0;
	for (e = rbranches; e < rbranches + total; e++) {
	    if (e->generation == 0 && e->peer == e)
		pack_count++;
	}

	if ((pack_count += g->nroot) > g->nroot) {
	    /* need to reallocate g->roots */
	    pack_base = g->roots;
	    g->roots = phgRealloc_(g->roots, pack_count * sizeof(*(g->roots)),
					g->nroot * sizeof(*(g->roots)));

	    /* update neighbour pointers in the local tree */
	    if (pack_base != g->roots) {
		/* note: have to use (char *) because old and new g->roots
		   may not be aligned */
		pack_offset = (char *)g->roots - (char *)pack_base;
		phgTraverseElements(g, redistribute_update_neighbours2);

		/* update pointers in simplices */
		for (s0 = simplices; s0 < simplices + simplices_count; s0++) {
		    e = s0->e;
		    if ((void *)e >= (void *)pack_base
			&& (void *)e < (void *)(pack_base + g->nroot))
			e = s0->e = (void *)(((char *)e) + pack_offset);
		    e0 = e->peer;
		    if ((void *)e0 >= (void *)pack_base
			&& (void *)e0 < (void *)(pack_base + g->nroot))
			e->peer = (void *)(((char *)e0) + pack_offset);
		}
	    }
	}

	if (simplices != NULL) {
	    phgFree(simplices);
	    simplices = NULL;
	    simplices_count = simplices_allocated = 0;
	}

	/* attach incoming branches to the local tree */
	for (e = rbranches; e < rbranches + total; e++) {
	    if (e->generation != 0)
		continue;
	    if (e->peer == e) {
		/* new root element */
		e0 = g->roots + (g->nroot++);
		*e0 = *e;
		e0->children[0] = e0->children[1] = NULL;
	    }
	    else {
		if ((void *)(e0 = e->peer) < (void *)g->roots
		    || (void *)e0 >= (void *)(g->roots + g->nroot))
		    e0 = e0->peer;
	    }
	    assert((void *)e0 >= (void *)g->roots
			    && (void *)e0 < (void *)(g->roots + g->nroot));
	    redistribute_add_branch(e0, e);
	}

	assert(g->nroot == pack_count);

	/* update neighbours of new elements */
	for (e = rbranches; e < rbranches + total; e++) {
	    e0 = e->peer;
	    assert(e0 != NULL &&
			((void *)e0 < (void *)rbranches
			  || (void *)e0 >= (void *)(rbranches + total)));
	    for (ii = 0; ii < NFace; ii++) {
		ELEMENT *e1 = e0->neighbours[ii];
		if ((void *)e1 >= (void *)rbranches
				&& (void *)e1 < (void *)(rbranches + total))
		    e0->neighbours[ii] = e1->peer;
	    }
	}

	if (rbranches != NULL)
	    phgFree(rbranches);

	/* match new neighbours */
	/* TODO: use update_neighbours() function in coarsen.c */

	/* collect unmatched faces */
	faces_count = 0;
	phgTraverseElements(g, redistribute_collect_faces);
	faces_allocated = faces_count;
	faces = phgAlloc(faces_count * sizeof(*faces));
	faces_count = 0;
	phgTraverseElements(g, redistribute_collect_faces);

	phgDebug((3, "faces_count=%d\n", faces_count));

	order = phgAlloc(faces_count * sizeof(*order));
	for (i = 0; i < faces_count; i++)
	    order[i] = i;
	assert(faces_count <= INT_MAX);
	qsort(order, faces_count, sizeof(*order), redistribute_comp_face);

	/* find neighbours */
	i0 = 0;
	while (i0 < faces_count) {
	    f0 = faces + order[i0];
	    if (f0->e == NULL)
		continue;
	    i1 = i0;
	    while (++i1 < faces_count &&
			!redistribute_comp_face0(order + i0, order + i1));
	    /*f1 = faces + order[i1];*/
	    /* the face is shared by (at most 4) elements between [i0, i1) */
	    while (i0 < i1) {
		i2 = i0;
		f0 = faces + order[i0];
		while (++i2 < i1) {
		    f2 = faces + order[i2];
		    phgDebug((4, "common face: %d %d %d, %p %p\n", f0->verts[0],
				f0->verts[1], f0->verts[2], f0->e, f2->e));
		    f0->e->neighbours[f0->vertex] = f2->e;
		    f0->e->bound_type[f0->vertex] &= ~REMOTE;
		    f2->e->neighbours[f2->vertex] = f0->e;
		    f2->e->bound_type[f2->vertex] &= ~REMOTE;
		}
		++i0;
	    }
	}
	phgFree(order);

	if (faces != NULL) {
	    phgFree(faces);
	    faces = NULL;
	    faces_count = faces_allocated = 0;
	}
    }

    /* adjust counters (FIXME: avoid this traversal?) */
    pack_count = nleaf = 0;
    phgTraverseAllElements(g, redistribute_count_elements);
    g->nleaf = nleaf;
    g->ntree = pack_count;
    phgDebug((3, "new nroot, nleaf, ntree: %d %d %d\n",
		g->nroot, g->nleaf, g->ntree));

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for merging branches: %0.4lf\n", __func__, t0);
	t0 = phgGetTime(NULL);
    }

    /* update RNEIGHBOURS */
    phgUpdateRNeighbours_(g);

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for updating rneighbours: %0.4lf\n", __func__, t0);
	t0 = phgGetTime(NULL);
    }

    phgInfo(2, "%d vertices, %d elements.\n", g->nvert, g->nleaf);

    if ((g->flags & EDGE_FLAG))
	phgRebuildL2Gmap_(g, EDGE);
#if Dim == 3
    if ((g->flags & FACE_FLAG))
	phgRebuildL2Gmap_(g, FACE);
#endif
    if ((g->flags & ELEM_FLAG))
	phgRebuildL2Gmap_(g, VOLUME);

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for rebuilding L2Gmaps: %0.4lf\n", __func__, t0);
	t0 = phgGetTime(NULL);
    }

    phgUpdateBoundaryTypes(g);

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for updating bdry types: %0.4lf\n", __func__, t0);
	t0 = phgGetTime(NULL);
    }

    /* redistribute DOFs or cancel redistribution */
    phgDofRedistribute(g, TRUE);

    if (phgVerbosity > 0) {
	t = phgGetTime(NULL) - t0;
	MPI_Reduce(&t, &t0, 1, MPI_DOUBLE, MPI_MAX, 0, g->comm);
	phgPrintf("%s - time for migrating DOF data: %0.4lf\n", __func__, t0);
	t0 = phgGetTime(NULL);
    }

    /* update g->lif (note: g->lif may have been updated by the partitioner).
     * TODO: merge this communication with other communications */
    if (g->lif < 0) {
	double w0[3], w[3];
	w0[0] = w0[1] = w0[2] = g->nleaf / phgGetProcessorSpeed(NULL);
	MPI_Allreduce(w0, w, 1, MPI_3DOUBLE, MPI_MSM, g->comm);
	g->lif = w[2] * g->nprocs / w[1];
    }
    if (phgVerbosity > 0)
	phgPrintf("Grid LIF: old = %lg, new = %lg\n",
				(double)old_lif, (double)g->lif);

#if DEBUG
    if (phgVerbosity >= 1 && weights != NULL && power != 0.) {
	FLOAT w = 0.;
	ForAllElements(g, e)
	    w += Pow(*DofElementData(weights, e->index), power);
	phgInfo(1, "weight = %lf, nleaf = %d\n", (double)w, g->nleaf);
    }
#endif

    PrintTime(1);

    Return;
}

int
phgBalanceGrid_(GRID *g, FLOAT lif_threshold, INT submesh_threshold,
			DOF *wgts, FLOAT pwr, MPI_Comm comm_world,
			PARTITIONER part_func)
/* returns 0 if grid unchanged, !0 if mesh partitioned or redistributed.
 *
 * 'lif_threshold' is the minimum LIF for redistributing the mesh,
 * used to control when to redistribute the mesh.
 *
 * 'submesh_threshold' gives minimum number of elements in each submesh,
 * used to control the number of submeshes.
 *
 * 'wgts' is a P0-type DOF providing element weights, the real weights
 * used is 'wgts' to the power 'pwr'.
 */
{
    FLOAT lif;
    int nprocs_max, rank;

    MPI_Comm_rank(comm_world, &rank);
    MPI_Comm_size(comm_world, &nprocs_max);

    if (submesh_threshold <= 0)
	submesh_threshold = 1;

    if (g->nelem_global >= submesh_threshold * g->nprocs) {
	nprocs_requested = (g->nelem_global + submesh_threshold - 1)
				/ submesh_threshold;
	if (nprocs_requested < g->nprocs)
	    nprocs_requested = g->nprocs;
    }
    else {
	nprocs_requested = g->nprocs;
    }

#if NO_NEW_COMMUNICATOR
    if (nprocs_requested > 1)
	nprocs_requested = nprocs_max;
#endif

    if (nprocs_requested <= 0)
	nprocs_requested = 1;
    else if (nprocs_requested > nprocs_max)
	nprocs_requested = nprocs_max;

    /* adjust lif according to new number of procs */
    lif = nprocs_requested * g->lif / g->nprocs;
    has_idle_processes = (g->nprocs < nprocs_max);
    if (rank >= g->nprocs || nprocs_requested > g->nprocs
	|| lif >= lif_threshold) {
	do {
	    weights = wgts;
	    power = pwr;
	    balancing_grid = TRUE;
	    phgRedistributeGrid_(g, comm_world, part_func);
	    /* restore flags */
	    balancing_grid = FALSE;
	    nprocs_requested = 0;
	    weights = NULL;
	    power = 0.;
	    if (rank < g->nprocs) {
		if (!grid_changed)
		    has_idle_processes = (g->nprocs < nprocs_max);
		return grid_changed;
	    }
	} while (!grid_changed);
    }

    return 0;
}

#else	/* USE_MPI */

PARTITIONER
phgGetPartitioner(int no)
{
    return NULL;
}

void
phgRedistributeGrid_(GRID *g, MPI_Comm comm_world, PARTITIONER part_func)
{
    Unused(g);
    Unused(comm_world);
    Unused(part_func);

    return;
}

int
phgBalanceGrid_(GRID *g, FLOAT lif_threshold, INT submesh_threshold,
			struct DOF_ *weights, FLOAT power, MPI_Comm comm_world,
			PARTITIONER part_func)
{
    Unused(g);
    Unused(lif_threshold);
    Unused(submesh_threshold);
    Unused(weights);
    Unused(power);
    Unused(comm_world);
    Unused(part_func);

    return 0;
}

#endif	/* USE_MPI */

void
phgRedistributeGrid(GRID *g)
{

    FunctionEntry_head;
    ParallelFunction(phgRedistributeGrid, g != NULL && phgNProcs > 1);
    FunctionEntry_tail;

    phgRedistributeGrid_(g, phgComm, NULL);

    Return;
}
