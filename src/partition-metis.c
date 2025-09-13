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

/* $Id: partition-metis.c,v 1.32 2020/10/21 07:40:42 zlb Exp $ */

#include "phg.h"

#if USE_MPI

#if USE_METIS || USE_PARMETIS
#include <math.h>

static int *nelem;
#endif	/* USE_METIS || USE_PARMETIS */

#if 0 * USE_METIS

#include <metis.h>

#ifndef METIS_VER_MAJOR
# define METIS_VER_MAJOR 4
#endif

#if METIS_VER_MAJOR <= 4
typedef idxtype idx_t;
#endif	/* METIS_VERSION_MAJOR <= 4 */

static idx_t *idx_array;

static BOOLEAN
metis_callback0 CB_ARGS(e)
{
    *(idx_array++) = e->verts[0];
    *(idx_array++) = e->verts[1];
    *(idx_array++) = e->verts[2];
#if Dim == 3
    *(idx_array++) = e->verts[3];
#endif
    return TRUE;
}

static BOOLEAN
metis_callback1 CB_ARGS(e) {nelem[e->mark = *(idx_array++)]++; return TRUE;}

BOOLEAN
phgPartitionMetis(GRID *g, MPI_Comm newcomm, DOF *weights, FLOAT power)
/* calls METIS to partition the grid into phgNProcs subgrids.
   The result of the partitioning is stored in the 'mark' field, i.e.,
   e->mark is the rank of the process the element is assigned to.

   returns FALSE if failure (e.g., some submesh is empty) */
{
    int ne = g->nleaf, nn = g->nvert, etype = Dim - 1;
    idx_t numflag = 0, edgecut;
    idx_t *elmnts = phgAlloc(NVert * sizeof(idx_t) * ne);
    idx_t *epart = phgAlloc(sizeof(idx_t) * ne);
    idx_t *npart = phgAlloc(sizeof(idx_t) * nn);
    int i, nelem_max;
    static BOOLEAN initialized = FALSE;

    Unused(newcomm);
    Unused(weights);
    Unused(power);

    if (!initialized) {
	/* register options (nothing to do, simply return) */
	initialized = TRUE;
	return;
    }

    assert(g->nprocs == 1);
    assert(nprocs == phgNProcs);

    if (!phgSameProcessorSpeed())
	phgPrintf("WARNING: processor speed ignored by METIS partitioner.\n");

    idx_array = elmnts;
    phgTraverseElements(g, metis_callback0);
    METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &phgNProcs,
			&edgecut, epart, npart);
    nelem = phgCalloc(phgNProcs, sizeof(*nelem));
    idx_array = epart;
    phgTraverseElements(g, metis_callback1);
    nelem_max = 0;
#if DEBUG
    for (i = 0; i < phgNProcs; i++)
	phgDebug((2, "leaf elements for p%d: %d\n", i, nelem[i])); 
#endif
    for (i = 0; i < phgNProcs; i++) {
	if (nelem[i] <= 0) break;
	if (nelem_max < nelem[i]) nelem_max = nelem[i];
    }
    phgFree(nelem);
    phgFree(npart);
    phgFree(epart);
    phgFree(elmnts);

    if (i >= phgNProcs) {
	if (phgSameProcessorSpeed())
	    g->lif = ((FLOAT)nelem_max * (FLOAT)phgNProcs) / ((FLOAT)g->nelem);
	else
	    g->lif = -1.0;	/* updated later in distribute.c */
	return TRUE;
    }
    else {
	return FALSE;
    }
}

#endif	/* USE_METIS */

#if USE_PARMETIS

#include <parmetis.h>

#ifdef SCOTCH_PARMETIS_VERSION		/* PTScotch compatibility lib */
# define PARMETIS_MAJOR_VERSION SCOTCH_PARMETIS_VERSION
typedef SCOTCH_Num idxtype;
#endif	/* defined(SCOTCH_PARMETIS_VERSION) */

#ifndef PARMETIS_MAJOR_VERSION
# define PARMETIS_MAJOR_VERSION 3
#endif

#if PARMETIS_MAJOR_VERSION <= 3
typedef idxtype idx_t;
typedef float real_t;
# define REAL_T		MPI_FLOAT
#endif	/* PARMETIS_MAJOR_VERSION <= 3 */

/* variables used by ParMETIS */
static GRID *g_ = NULL;
static idx_t *part;
static idx_t *xadj, *adjncy, *px, *py;
static INT no0 = 0, *map;	/* global numbering of the first element */
static idx_t nedges, nparts;	/* edges in the graph (adjncy size) */
static int *remote_index;	/* global index of remote neighbour */
static int *local_index;	/* global index of local element */

static BOOLEAN
parmetis_callback0 CB_ARGS(e)
{
    int i;

    if (!(g_->flags & ELEM_FLAG))
	e->index = phgTraverseIndex;
    map[e->index] = phgTraverseIndex;
    for (i = 0; i < NFace; i++) {
	if (HasLocalNeighbour(e, i)) {
	    nedges++;
	}
	else if (HasRemoteNeighbour(e, i)) {
	    nedges++;
	    local_index[(size_t)e->neighbours[i]] = no0 + phgTraverseIndex;
	}
    }

    return TRUE;
}

static BOOLEAN
parmetis_callback1 CB_ARGS(e)
{
    int i;

    for (i = 0; i < NFace; i++) {
	if (HasLocalNeighbour(e, i))
	    *(py++) = no0 + map[GetNeighbour(e, i)->index];
	else if (HasRemoteNeighbour(e, i))
	    *(py++) = remote_index[(size_t)e->neighbours[i]];
    }
    *(++px) = py - adjncy;
    return TRUE;
}

static BOOLEAN
parmetis_callback2 CB_ARGS(e)
{
    assert(*px >= 0 && *px < nparts);
    nelem[e->mark = *(px++)]++;

    return TRUE;
}

BOOLEAN
phgPartitionParMetis(GRID *g, MPI_Comm newcomm, DOF *weights, FLOAT power)
/* calls ParMETIS to redistribute the grid into nparts processes.
 * The result of the partitioning is stored in the 'mark' field, i.e.,
 * e->mark is the rank of the process the element is assigned to.
 *
 * Returns FALSE if failure or grid unchanged
 */
{
    int nprocs;
    ELEMENT *e;
    idx_t wgtflag = 0;
    idx_t numflag = 0;
    idx_t ncon = 0;
    idx_t *vtxdist = NULL;
    idx_t *vwgt = NULL;
    idx_t *adjwgt = NULL;

    real_t *tpwgts = NULL;
    real_t ubvec[1] = {1.05};
    /* large value of ipc2redist makes ParMETIS minimize edge cuts */
#ifndef SCOTCH_PARMETIS_VERSION
    real_t ipc2redist = 1000000.0;	/* 0.000001 to 1000000.0 */
#endif
    idx_t options[4] = {1, /*verbosity=*/ 0, /*seed=*/15, /*coupled*/1};
    idx_t edgecut;
    MPI_Comm comm;
    int i, nprocs0, rank0, nleaf_max, *counts, ret;
    static BOOLEAN initialized = FALSE;
#if DEBUG
    int j;
#endif

    if (!initialized) {
	/* register options (nothing to do, simply return) */
	initialized = TRUE;
	return FALSE;
    }

    MPI_Comm_size(newcomm, &nprocs);
    nparts = nprocs;

    g_ = g;

    if (g->nprocs == 1 && nprocs == 1)
	return FALSE;

    if (!phgSameProcessorSpeed()) {
	real_t speed;
	if (g->comm != MPI_COMM_NULL)
	    tpwgts = phgAlloc(nparts * sizeof(*tpwgts));
	speed = phgGetProcessorSpeed(NULL);
	MPI_Gather(&speed, 1, REAL_T, tpwgts, 1, REAL_T, 0, newcomm);
	if (g->comm != MPI_COMM_NULL) {
	    ncon = 1;
	    if (phgRank == 0) {
		speed = 0.;
		for (i = 0; i < nparts; i++)
		    speed += tpwgts[i];
		for (i = 0; i < nparts; i++)
		    tpwgts[i] /= speed;
	    }
	    MPI_Bcast(tpwgts, nparts, REAL_T, 0, g->comm);
	}
    }

    if (g->comm == MPI_COMM_NULL)
	return FALSE;

    /* one extra position for the flag */
    counts = phgAlloc((nparts + 1) * sizeof(*counts));
    nelem = phgCalloc(nparts + 1, sizeof(*nelem));
    part = phgAlloc(g->nleaf * sizeof(*part));

    /* ParMETIS does not work well with a small number of elements.
     * Below is a quick workaround for special cases */
    if (nparts > 1 && g->nleaf_global <= 3 * nparts) {
	int j, bsize, mod, n0, *part0;
	i = g->nleaf;
	MPI_Scan(&i, &n0, 1, MPI_INT, MPI_SUM, g->comm);
	n0 -= g->nleaf;
	bsize = g->nleaf_global / nparts + 1;
	mod = g->nleaf_global % nparts;
	part0 = phgAlloc((nparts + 1) * sizeof(*part0));
	part0[0] = 0;
	for (i = 0, j = 0; i < nparts; i++) {
	    if (n0 >= part0[i]) {
		j = i;
	    }
	    part0[i + 1] = part0[i] + (i < mod ? bsize : bsize - 1);
	}
	for (i = 0; i < g->nleaf; i++, n0++) {
	    while (n0 >= part0[j + 1])
		j++;
	    part[i] = j;
	}
	phgFree(part0);
	xadj = NULL;
	adjncy = NULL;
	map = NULL;
	remote_index = NULL;
	local_index = NULL;
	goto skip_parmetis;
    }

    /* Communicator for ParMETIS consisting of procs with g->nleaf != 0 */
#if NO_NEW_COMMUNICATOR
    comm = g->comm;
    nprocs0 = g->nprocs;
    rank0 = g->rank;
#else	/* NO_NEW_COMMUNICATOR */
    MPI_Comm_split(g->comm, g->nleaf == 0 ? 0 : 1, 0, &comm);
    MPI_Comm_size(comm, &nprocs0);
    MPI_Comm_rank(comm, &rank0);
    if (nprocs0 == g->nprocs) {
	/* Use MPI_Comm_dup to preserve the attributes of phgComm */
	MPI_Comm_free(&comm);
	MPI_Comm_dup(g->comm, &comm);
	nprocs0 = g->nprocs;
	rank0 = g->rank;
    }
#endif	/* NO_NEW_COMMUNICATOR */

    vtxdist = phgAlloc((nprocs0 + 1) * sizeof(*vtxdist));

    if (g->neighbours.count > 0) {
	remote_index = phgAlloc(g->neighbours.count * sizeof(*remote_index));
	local_index = phgAlloc(g->neighbours.count * sizeof(*local_index));
    }
    else {
	local_index = remote_index = NULL;
	if (g->neighbours.counts == NULL)
	    g->neighbours.counts =
			phgAlloc(g->nprocs * sizeof(*g->neighbours.counts));
	if (g->neighbours.displs == NULL)
	    g->neighbours.displs =
			phgAlloc(g->nprocs * sizeof(*g->neighbours.displs));
	for (i = 0; i < g->nprocs; i++)
	    g->neighbours.counts[i] = g->neighbours.displs[i] = 0;
    }

    if ((i = g->nleaf) > 0) {
	MPI_Allgather(&i, 1, MPI_INT, counts, 1, MPI_INT, comm);
	vtxdist[0] = 0;
	for (i = 0; i < nprocs0; i++) {
	    if (i == rank0)
		no0 = vtxdist[i];
	    vtxdist[i + 1] = vtxdist[i] + counts[i];
	}
    }

    /* first, traverse leaf elements, set mark to the global element index,
       and count number of neighbours (i.e., number of edges of the graph) */
    map = phgAlloc(((g->flags & ELEM_FLAG) ? g->nelem : g->nleaf)
			* sizeof(*map));
    nedges = 0;
    phgTraverseElements(g, parmetis_callback0);

    /* exchange element indices of remote neighbours */
    if (g->nprocs > 1)
	MPI_Alltoallv(local_index, g->neighbours.counts, g->neighbours.displs,
		      MPI_INT,
		      remote_index, g->neighbours.counts, g->neighbours.displs,
		      MPI_INT, g->comm);

    xadj = phgAlloc((g->nleaf + 1) * sizeof(*xadj));
    adjncy = phgAlloc(nedges * sizeof(*adjncy));
    if (weights != NULL && power != 0.) {
	FLOAT a, d, w_max = -1e10, w_min = 1e10, tmp[4];
	if (weights->count_elem != 1)
	    phgError(1, "invalid weight!\n");
	/* compute upper/lower bounds of 'weight' */
	ForAllElements(g, e) {
	    a = *DofElementData(weights, e->index);
	    if (w_max < a)
		w_max = a; 
	    if (w_min > a)
		w_min = a; 
	}
	tmp[0] = w_max;
	tmp[1] = -w_min;
	MPI_Allreduce(tmp, tmp + 2, 2, PHG_MPI_FLOAT, PHG_MAX, comm);
	w_max = tmp[2];
	w_min = -tmp[3];
	phgInfo(2, "w_max = %lf, w_min = %lf\n", (double)w_max, (double)w_min);
	if (w_min < 0.) {
	    phgWarning("negative weights, ignoring arguments 'weights' and "
			"'power'!\n");
	    weights = NULL;
	    power = 0.;
	}
	else if (w_max > w_min) {
	    w_max = Pow(w_max, power);
	    w_min = Pow(w_min, power);
	    vwgt = phgAlloc(g->nleaf * sizeof(*vwgt));
	    d = 10000. / (w_max - w_min);
	    ForAllElements(g, e) {
		a = Pow(*DofElementData(weights, e->index), power);
		vwgt[map[e->index]] = Round((a - w_min) * d);
	    }
	    wgtflag = 2;	/* 2: weights on vertices */
	}
    }

    /* next, traverse leaf elements to fill xadj and adjncy */
    *(px = xadj) = 0;
    py = adjncy;
    phgTraverseElements(g, parmetis_callback1);

#if 0
    /*if (nparts != g->nprocs)*/ {
	options[3] = 2;	/* submeshes and processes are decoupled */
	for (i = 0; i < g->nleaf; i++)
	    part[i] = g->rank;
    }
#endif

    /* call ParMETIS to redistribute the elements.
     *
     * Note: it seems ParMETIS_V3_AdaptiveRepart does not like the case
     * g->nprocs != nparts (contrary to what's described in its manual),
     * also it's not available in the PTScotch compatibility lib,
     * and ParMETIS_V3_PartKway is called instead in these cases.
     */
    ret = METIS_OK;
#ifndef SCOTCH_PARMETIS_VERSION
    if ((NO_NEW_COMMUNICATOR || g->nleaf > 0) && nprocs0 == nparts) {
	idx_t *vsize = NULL;
#if PARMETIS_MAJOR_VERSION >= 4 
	/* ncon can't be 0 and tpwgts/ubvec can't be NULL (ParMETIS-4 bug?) */
	real_t d;
	vsize = phgAlloc(g->nleaf * sizeof(*vsize));
	for (i = 0; i < g->nleaf; i++)
	    vsize[i] = 1;
	if (tpwgts == NULL) {
	    d = (ncon = 1) / (real_t)nparts;
	    tpwgts = phgAlloc(nparts * sizeof(*tpwgts));
	    for (i = 0; i < nparts; i++)
		tpwgts[i] = d;
	}
	if (vwgt == NULL) {
	    vwgt = phgAlloc(g->nleaf * sizeof(*vwgt));
	    for (i = 0; i < g->nleaf; i++)
		vwgt[i] = 1;
	}
#endif	/* PARMETIS_MAJOR_VERSION >= 4 */
	ret = ParMETIS_V3_AdaptiveRepart(vtxdist, xadj, adjncy, vwgt, vsize,
		adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec,
		&ipc2redist, options, &edgecut, part, &comm);
#if PARMETIS_MAJOR_VERSION >= 4 
	phgFree(vsize);
#endif	/* PARMETIS_MAJOR_VERSION >= 4 */
    }
    else if (NO_NEW_COMMUNICATOR || g->nleaf > 0) 
#endif	/* !defined(SCOTCH_PARMETIS_VERSION) */
    {
	if (/*vwgt != NULL &&*/ tpwgts == NULL) {
	    /* ParMETIS_V3_PartKway: tpwgts cannot be NULL */
	    real_t d = (ncon = 1) / (real_t)nparts;
	    tpwgts = phgAlloc(nparts * sizeof(*tpwgts));
	    for (i = 0; i < nparts; i++)
		tpwgts[i] = d;
	}
	ret = ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt,
		&wgtflag, &numflag, &ncon, &nparts,
		tpwgts, ubvec, options, &edgecut, part, &comm);
    }

    if (ret != METIS_OK) {
	phgInfo(-1, "Error occurred in ParMETIS (ret=%d), abort.\n", ret);
	phgAbort(1);
    }

#if !NO_NEW_COMMUNICATOR
    MPI_Comm_free(&comm);
#endif	/* !NO_NEW_COMMUNICATOR */

skip_parmetis:
    /* finally, set mark to the new process number the element is assigned to */
    px = part;
    phgTraverseElements(g, parmetis_callback2);

    /* statistics */
    nelem[nparts] = (nelem[g->rank] == g->nleaf) ? 1 : 0;
    MPI_Allreduce(nelem, counts, nparts + 1, MPI_INT, MPI_SUM, g->comm);
    nleaf_max = 0;
    for (i = 0; i < nparts; i++) {
	/*if (counts[i] <= 0) break;*/
	if (nleaf_max < counts[i])
	    nleaf_max = counts[i];
    }

    numflag = (counts[nparts] == g->nprocs);	/* TRUE ==> grid not changed */
#if DEBUG
    phgInfo(3, "number of leaf elements for each process:\n");
    for (j = 0; j < nparts; j++)
	phgInfo(3, "proc %d: %d\n", j, nelem[j]);
#endif

    phgFree(vwgt);
    phgFree(tpwgts);
    phgFree(nelem);
    phgFree(xadj);
    phgFree(adjncy);
    phgFree(part);
    phgFree(counts);
    phgFree(remote_index);
    phgFree(local_index);
    phgFree(vtxdist);
    phgFree(map);

    if (i >= nparts) {
	if (numflag) {
	    if (g->rank == 0)
		phgInfo(1, "grid unchanged.\n");
	    return FALSE;
	}
	if (phgSameProcessorSpeed())
	    g->lif = (nleaf_max * (FLOAT)nparts) / ((FLOAT)g->nelem_global);
	else
	    g->lif = -1.0;	/* updated later in distribute.c */
	return TRUE;
    }
    else {
	if (g->rank == 0 && phgVerbosity > 0)
	    phgWarning("failure in ParMETIS, grid unchanged.\n");
	return FALSE;
    }
}

#endif	/* USE_PARMETIS */

#endif	/* USE_MPI */
