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

/* $Id: partition-rtk.c,v 1.40 2020/05/07 05:37:21 zlb Exp $ */
#include "phg.h"
#include "phg/partition-utils.h"

#if USE_MPI  /* whole file */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "phg/partition-rtk.h"

typedef struct {
    INT global;
    INT nleaf;
}tree;

/* used for ordering local roots */
static int
comp_tree(const void *n1, const void *n2)
{
    tree *m1, *m2;
    int flag = 0;

    m1 = (tree *) n1;
    m2 = (tree *) n2;
    if (m1->global > m2->global)
	flag = 1;
    else if (m1->global == m2->global)
	flag = 0;
    else if (m1->global < m2->global)
	 flag = -1;
    return flag;
}

/* global variable that records the number of leaves in the tree */
/* traverse a tree */
static INT TreeIndex = 0;
static INT TreeDepth = 0;

/* if e is child[0] */
static BOOLEAN IS_LEFT;

/* if e->parent->type is M or O */
static BOOLEAN EPMO;

static BOOLEAN
traverse_tree(ELEMENT *e, BOOLEAN (*cb) (ELEMENT *e))
{
    BOOLEAN is_leaf = TRUE;
    if (TreeDepth == 1) {
	if (((ELEMENT *)e->parent)->type == OPPOSITE ||
	    ((ELEMENT *)e->parent)->type == MIXED)
	    EPMO = TRUE;
    }
    else {
	EPMO = FALSE;
    }

    if (IS_LEFT || EPMO) {   /* access child[0] */
	if (e->children[0] != NULL) {
	    is_leaf = FALSE;
	    IS_LEFT = TRUE;

	    TreeDepth++;
	    if (!traverse_tree(e->children[0], cb))
		return FALSE;
	    TreeDepth--;
	}
	if (e->children[1] != NULL) {
	    IS_LEFT = FALSE;
	    is_leaf = FALSE;

	    TreeDepth++;
	    if (!traverse_tree(e->children[1], cb))
		return FALSE;
	    TreeDepth--;
	}
    }
    else {	    /* access child[1] */
	if (e->children[1] != NULL) {
	    IS_LEFT = FALSE;
	    is_leaf = FALSE;

	    TreeDepth++;
	    if (!traverse_tree(e->children[1], cb))
		return FALSE;
	    TreeDepth--;
	}
	if (e->children[0] != NULL) {
	    IS_LEFT = TRUE;
	    is_leaf = FALSE;

	    TreeDepth++;
	    if (!traverse_tree(e->children[0], cb))
		return FALSE;
	    TreeDepth--;
	}
    }

    if (is_leaf) {
	BOOLEAN ret = cb(e);
	TreeIndex++;
	return ret;
    }
    else {
	return TRUE;
    }
}

/* calculate number of leaves of roots */
static BOOLEAN
tree_leaf(ELEMENT *e)
{
    return TRUE;
}

/* global = npart * mm + rr */
/* 0 <= rr < npart */
static INT npart;
static INT rr;
static INT mm;

/* mark elements */
/* assign elements to partition */
/* no weights */
static BOOLEAN
part_mark(ELEMENT *e)
{
    INT flag = rr * (mm + 1);
    
    if (mm > 0) {
	if (TreeIndex < flag) {
	    e->mark = TreeIndex / (mm + 1);
	}
	else {
	    e->mark = rr +(TreeIndex - flag) / mm;
	}
    }
    else {
	e->mark = TreeIndex;
    }
    return TRUE;
}

typedef struct {
    INT local;
    INT global;
}PID;

/* comparision function for PID */
static int
comp_pid(const void *n1, const void *n2)
{
    PID *m1, *m2;
    int flag = 0;

    m1 = (PID *) n1;
    m2 = (PID *) n2;
    if (m1->global > m2->global)
	flag = 1;
    else if (m1->global == m2->global)
	flag = 0;
    else if (m1->global < m2->global)
	 flag = -1;
    return flag;
}

/* for weighted condition */
static long int WgtsIndex = 0;

/* weights sum of tree */
static BOOLEAN
part_sum(ELEMENT *e)
{
    WgtsIndex += e->mark;
    return TRUE;
}

/* assign e to a partition */
/* use weights */
static BOOLEAN
part_mark2(ELEMENT *e)
{
    long int index;

    if (WgtsIndex > 0) {
	index = WgtsIndex + e->mark - 1;
	WgtsIndex += e->mark;
    }
    else {
	index = WgtsIndex;
    }

    if (mm > 0) {
	e->mark = index / mm;
    }
    else {
	e->mark = index;
    }
    return TRUE;
}


BOOLEAN
phgPartitionRTK(GRID *g, MPI_Comm newcomm, DOF *weights, FLOAT power)
/* k-way refinement tree repartitioner 
 * redistribute grid into nprocs processors
 * submesh X will be move to pid X
 * this procedure allows multi-grid  */
{
    INT i;
    INT sum;
    int nprocs, nprocs0, rank0;
    INT nleaf_max = 0;
    ELEMENT *e;
    BOOLEAN ret;
    BOOLEAN FirstCall;
    BOOLEAN UseWeights = FALSE;
    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	/* Register initial grid order options */
	phgGridInitOrder(NULL, FALSE);
	initialized = TRUE;
	return FALSE;
    }

    if (g == NULL)
	phgError(1, "GRID is NULL! %s  %d\n", __FILE__, __LINE__);
    MPI_Comm_size(newcomm, &nprocs);
    if (nprocs == 1 && g->nprocs == 1)
	return FALSE;

    if (!phgSameProcessorSpeed())
	phgPrintf("WARNING: processor speed ignored by RTK partitioner.\n");

    if (g->comm == MPI_COMM_NULL)
	return FALSE;

    if (nprocs < 0 || nprocs > phgNProcs) {
	npart = phgNProcs;
    }
    else {
	npart = nprocs;
    }

    /* determine if g is initialized or not */
    if (phgRank == 0) {
	if (!strcmp(phgOptionsGetKeyword("-rtk_init_order"), "none"))
	    phgWarning("\"-rtk_init_order none\" used, grid may not be properly initialized for the RTK partitioner.\n");
	if (g->last_partitioner == -1) {
	    FirstCall = TRUE;
	}
	else {
	    FirstCall = FALSE;
	}
    }

    /* synchronize status */
    if (g->nprocs > 1) {
	MPI_Bcast(&FirstCall, 1, PHG_MPI_INT, 0, g->comm);
    }

    nprocs0 = g->nprocs;
    rank0 = g->rank;


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
	MPI_Allreduce(tmp, tmp + 2, 2, PHG_MPI_FLOAT, PHG_MAX, g->comm);
	w_max = tmp[2];
	w_min = -tmp[3];
	phgInfo(2, "w_max = %lf, w_min = %lf\n", (double)w_max, (double)w_min);
	if (w_min < 0.) {
	    phgWarning("negative weights, ignoring arguments 'weights' and "
		    "'power'!\n");
	}
	else if (w_max == w_min) {
	    phgWarning("constant weights, ingoring arguments 'weights' and "
		      "'power'! \n");
	}
	else if (w_max > w_min) {
	    UseWeights = TRUE;
	    w_max = Pow(w_max, power);
	    w_min = Pow(w_min, power);
	    d = 10000. / w_max;
	    ForAllElements(g, e) {
		a = Pow(*DofElementData(weights, e->index), power);
		e->mark = Round(a * d);
	    }
	}
    } /* weights != NULL && power != 0. */

    if (UseWeights) {	/* use weights */
	if (FirstCall) {
	    tree *leaf, *tleaf;
	    INT nroot = g->nroot;
	    INT ntree, kk;
	    MPI_Datatype type;
	    MPI_Datatype types[2];
	    MPI_Aint indices[2];
	    int blocklens[2], *dis, *cnt;
	    INT *recvcount;
	    INT *nelem;
	    FLOAT enlarge;


	    /* total roots in g */
	    MPI_Allreduce(&nroot, &ntree, 1, PHG_MPI_INT, MPI_SUM, g->comm);
	    leaf = (tree *) phgAlloc(nroot * sizeof(tree));
	    tleaf = (tree *) phgAlloc(ntree * sizeof(tree));
	    nelem = (INT *)phgAlloc(npart * 2 * sizeof (INT));

	    /* collect number of leaf of roots */
	    for (i = 0; i < nroot; i++) {
		e = g->roots + i;
		WgtsIndex = 0;
		IS_LEFT = TRUE;
		if (!traverse_tree(e, part_sum))
		    phgError(1, "Elements traversing error %s:%d\n", __FILE__,
			    __LINE__);
		leaf[i].global = GlobalElement(g, e->index);
		leaf[i].nleaf = WgtsIndex;
	    }

	    /* create a new datatype for exchange information */
	    /* used by structure tree */
	    blocklens[0] = blocklens[1] = 1;
	    types[0] = types[1] = PHG_MPI_INT;
#if 1	/* MPI-2 */
	    MPI_Get_address(&tleaf[0].global, &indices[0]);
	    MPI_Get_address(&tleaf[0].nleaf, &indices[1]);
#else	/* MPI-1 */
	    MPI_Address(&tleaf[0].global, &indices[0]);
	    MPI_Address(&tleaf[0].nleaf, &indices[1]);
#endif
	    indices[1] -= indices[0];
	    indices[0] = 0;
#if 1	/* MPI-2 */
	    MPI_Type_create_struct(2, blocklens, indices, types, &type);
#else
	    MPI_Type_struct(2, blocklens, indices, types, &type);
#endif
	    MPI_Type_commit(&type);

	    /* gather number of leaves */
	    recvcount = phgAlloc(nprocs0 * sizeof(*recvcount));
	    dis = phgAlloc(2 * nprocs0 * sizeof(*dis));
	    cnt = dis + nprocs0;

	    MPI_Allgather(&nroot, 1, PHG_MPI_INT, recvcount, 1, PHG_MPI_INT,
		    g->comm);
	    sum = 0;
	    for (i = 0; i < nprocs0; i++) {
		cnt[i] = recvcount[i];
		dis[i] = sum;
		assert(dis[i] == sum);
		sum += recvcount[i];
	    }

	    /* gather leaves infomation: leaf number */
	    MPI_Allgatherv(leaf, nroot, type, tleaf, cnt, dis, type, g->comm);
	    MPI_Type_free(&type);
	    qsort(tleaf, ntree, sizeof(*tleaf), comp_tree);

	    /* determine element index */
	    tleaf[0].global = 0;
	    for (i = 1; i < ntree; i++)
		tleaf[i].global = tleaf[i - 1].global +tleaf[i - 1].nleaf;

	    /* partitioning */
	    enlarge = (npart - 0.5) / (npart * (npart - 1.));
	    if (enlarge > 1.05 / npart) enlarge = 1.05 / npart;
	    
	    mm = (INT)((tleaf[ntree - 1].global + tleaf[ntree - 1].nleaf) 
		    * enlarge);
	    for (i = 0; i < nroot; i++) {
		e = g->roots + i;
		kk = GlobalElement(g, e->index);
		WgtsIndex = tleaf[kk].global;
		IS_LEFT = TRUE;
		traverse_tree(e, part_mark2);
	    } /* finish partitioning */

	    /* determine nleaf_max */
	    for (i = 0; i < npart; i++) {
		nelem[i] = 0;
	    }
	    ForAllElements(g, e) {
		nelem[e->mark] += 1;
	    }
	    MPI_Scan(nelem, nelem + npart, npart, PHG_MPI_INT, MPI_SUM,
			g->comm);
	    nleaf_max = 0;
	    for (i = 0; i < npart; i++) {
		if (nelem[npart + i] > nleaf_max)
		    nleaf_max = nelem[npart + i];
	    }

	    phgFree(leaf);
	    phgFree(tleaf);
	    phgFree(recvcount);
	    phgFree(dis);
	    phgFree(nelem);

	} /* weights && firstcall */
	else {	      /* weighted */
	    long int *leaf;
	    long int  wgts = 0;
	    PID *proc;
	    INT nroot = g->nroot;
	    INT *nelem;
	    FLOAT enlarge;

	    leaf = (long int *)phgAlloc(nprocs0 * sizeof(long int));
	    proc = (PID *)phgAlloc(nroot * sizeof(PID));
	    nelem = (INT *)phgAlloc(npart * 2 * sizeof (INT));

	    /* determine local accessing order */
	    for (i = 0; i < nroot; i++) {
		e = g->roots + i;
		proc[i].local = i;
		proc[i].global = GlobalElement(g, e->index);
	    }
	    qsort(proc, nroot, sizeof(*proc), comp_pid);

	    wgts = 0;
	    ForAllElements(g, e) {
		wgts += e->mark;
	    }
	    MPI_Allgather(&wgts, 1, MPI_LONG, leaf, 1, MPI_LONG, g->comm);
	    for (i = 1; i < nprocs0; i++)
		leaf[i] += leaf[i - 1];

	    /* partitioning */
	    enlarge = (npart - 0.5) / (npart * (npart - 1.));
	    if (enlarge > 1.05 / npart) enlarge = 1.05 / npart;
	    mm = (INT)(leaf[nprocs0 - 1] * enlarge);
	    WgtsIndex = leaf[rank0] - wgts;
	    for (i = 0; i < nroot; i++){
		e = g->roots + proc[i].local;
		IS_LEFT = TRUE;
		traverse_tree(e, part_mark2);
	    } /* finish partitioning */

	    for (i = 0; i < npart; i++) {
		nelem[i] = 0;
	    }
	    ForAllElements(g, e) {
		nelem[e->mark] += 1;
	    }
	    MPI_Scan(nelem, nelem + npart, npart, PHG_MPI_INT, MPI_SUM,
			g->comm);
	    nleaf_max = 0;
	    for (i = 0; i < npart; i++) {
		if (nelem[npart + i] > nleaf_max)
		    nleaf_max = nelem[npart + i];
	    }

	    phgFree(leaf);
	    phgFree(nelem);
	    phgFree(proc);
	} /* weights && not first call */
    }  /* use weights */
    else {     /* no weights */
	if (FirstCall) {		 /* first time call RTK repartitioner */
	    tree *leaf, *tleaf;
	    INT nroot = g->nroot;
	    INT ntree, kk;		  /* sum of g->nroot of all submeshes */
	    MPI_Datatype type;
	    MPI_Datatype types[2];
	    MPI_Aint indices[2];
	    int blocklens[2], *dis, *cnt;
	    INT *recvcount;

	    /* total roots in g */
	    MPI_Allreduce(&nroot, &ntree, 1, PHG_MPI_INT, MPI_SUM, g->comm);
	    leaf = (tree *) phgAlloc(nroot * sizeof(tree));
	    tleaf = (tree *) phgAlloc(ntree * sizeof(tree));

	    /* collect number of leaf of roots */
	    for (i = 0; i < nroot; i++) {
		e = g->roots + i;
		TreeIndex = 0;
		IS_LEFT = TRUE;
		if (!traverse_tree(e, tree_leaf))
		    phgError(1, "Elements traversing error %s:%d\n", __FILE__,
			    __LINE__);
		leaf[i].global = GlobalElement(g, e->index);
		leaf[i].nleaf = TreeIndex;
	    }

	    /* create a new datatype for exchange information */
	    /* used by structure tree */
	    blocklens[0] = blocklens[1] = 1;
	    types[0] = types[1] = PHG_MPI_INT;
#if 1	/* MPI-2 */
	    MPI_Get_address(&tleaf[0].global, &indices[0]);
	    MPI_Get_address(&tleaf[0].nleaf, &indices[1]);
#else	/* MPI-1 */
	    MPI_Address(&tleaf[0].global, &indices[0]);
	    MPI_Address(&tleaf[0].nleaf, &indices[1]);
#endif
	    indices[1] -= indices[0];
	    indices[0] = 0;
#if 1	/* MPI-2 */
	    MPI_Type_create_struct(2, blocklens, indices, types, &type);
#else	/* MPI-1 */
	    MPI_Type_struct(2, blocklens, indices, types, &type);
#endif
	    MPI_Type_commit(&type);

	    /* gather number of leaves */
	    recvcount = phgAlloc(nprocs0 * sizeof(*recvcount));
	    dis = phgAlloc(2 * nprocs0 * sizeof(*dis));
	    cnt = dis + nprocs0;

	    MPI_Allgather(&nroot, 1, PHG_MPI_INT, recvcount, 1, PHG_MPI_INT,
		    g->comm);
	    kk = 0;
	    for (i = 0; i < nprocs0; i++) {
		cnt[i] = recvcount[i];
		dis[i] = kk;
		assert(dis[i] == kk);
		kk += recvcount[i];
	    }

	    /* gather leaves infomation: leaf number */
	    MPI_Allgatherv(leaf, nroot, type, tleaf, cnt, dis, type, g->comm);
	    MPI_Type_free(&type);
	    qsort(tleaf, ntree, sizeof(*tleaf), comp_tree);

	    /* determine element index */
	    tleaf[0].global = 0;
	    for (i = 1; i < ntree; i++)
		tleaf[i].global = tleaf[i - 1].global +tleaf[i - 1].nleaf;

	    /* partitioning */
	    mm = g->nleaf_global / npart;
	    rr = g->nleaf_global % npart;
	    for (i = 0; i < nroot; i++) {
		e = g->roots + i;
		kk = GlobalElement(g, e->index);
		TreeIndex = tleaf[kk].global;
		IS_LEFT = TRUE;
		traverse_tree(e, part_mark);
	    } /* finish partitioning */

	    if (rr > 0) {
		nleaf_max = mm + 1;
	    }
	    else {
		nleaf_max = mm;
	    }

	    phgFree(leaf);
	    phgFree(tleaf);
	    phgFree(recvcount);
	    phgFree(dis);

	}  /* no weights && firstcall */
	else {
	    INT *leaf;
	    PID *proc;
	    INT nroot = g->nroot;

	    leaf = (INT *)phgAlloc(nprocs0 * sizeof(INT));
	    proc = (PID *)phgAlloc(nroot * sizeof(PID));

	    /* determine local accessing order */
	    for (i = 0; i < nroot; i++) {
		e = g->roots + i;
		proc[i].local = i;
		proc[i].global = GlobalElement(g, e->index);
	    }
	    qsort(proc, nroot, sizeof(*proc), comp_pid);

	    /* determine element index */
	    MPI_Allgather(&g->nleaf, 1,PHG_MPI_INT, leaf, 1, PHG_MPI_INT,
				g->comm);
	    for (i = 1; i < nprocs0; i++) {
		leaf[i] += leaf[i - 1];
	    }

	    /* partitioning */
	    mm = g->nleaf_global / npart;
	    rr = g->nleaf_global % npart;
	    TreeIndex = leaf[rank0] - g->nleaf;
	    for (i = 0; i < nroot; i++){
		e = g->roots + proc[i].local;
		IS_LEFT = TRUE;
		traverse_tree(e, part_mark);
	    } /* finish partitioning */

	    if (rr > 0) {
		nleaf_max = mm + 1;
	    }
	    else {
		nleaf_max = mm;
	    }

	    phgFree(leaf);
	    phgFree(proc);

	} /* no weights && not first call */
    }

    /* determine if GRID g is changed or not */
    i = 0;
    ForAllElements(g, e){
	if (e->mark == phgRank)
	    i += 1;
    }
    if (i == g->nleaf) {   /* grid unchanged */
	i = 0;
    }
    else {		   /* grid changed */
	i = 1;
    }
    MPI_Allreduce(&i, &sum, 1, PHG_MPI_INT, MPI_SUM, g->comm);
    if (sum > 0) {	   /* grid changed */
	ret = TRUE;
    }
    else {		   /* grid unchanged */
	ret = FALSE;
    }


    /* update information */
    if (ret) {
	if (phgSameProcessorSpeed())
	    g->lif = nleaf_max * npart / ((FLOAT)g->nelem_global);
	else
	    g->lif = -1;	/* computed later in distribute.c */
    }
    else {
	if (g->rank == 0 && phgVerbosity > 0)
	    phgWarning("grid unchanged.\n");
    }


    return ret;
}

#endif /* USE_MPI, whole file */
