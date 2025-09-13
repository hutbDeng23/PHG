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

/* $Id: partition-utils.c,v 1.129 2022/04/30 07:18:07 zlb Exp $ */

#include "phg.h"
#include "phg/partition-utils.h"
#include "phg/hamilton.h"

#include <stdio.h>
#include <stdlib.h>

/* getpid() */
#include <sys/types.h>
#include <unistd.h>

#include <string.h>	/* memcpy() */
#include <assert.h>

#if USE_MPI
/* load balance efficiency */
static const FLOAT SFC_LBE = 1.004;
#endif	/* USE_MPI */
static const SFC *sfc;

/* ************************************************
 * Hilbert Order 
 * used by phgImport
 **************************************************/

/* ordering:
 * used for ordering Hilbert order and Morton order*/
static int
comp_key(const void *n1, const void *n2)
{
    const SFC *m1, *m2;
    int flag = 0;

    m1 = sfc + *(const INT *)n1;
    m2 = sfc + *(const INT *)n2;
    if (m1->sfc > m2->sfc)
	flag = 1;
    else if (m1->sfc < m2->sfc)
	flag = -1;
    else 
	 flag = 0;
    return flag;
}

/* generate Hilbert Order */
/* 0: success, or fail */
int
phgHilbertOrder(GRID *g)
{
    SFC *x;
    INT i, *x_index, *ordering;
    SFC_FLOAT xmax, xmin;
    SFC_FLOAT ymax, ymin;
    SFC_FLOAT zmax, zmin;
    SFC_FLOAT etx, ety, etz;
    FLOAT x0, y0, z0;
    static const FLOAT lam[] = {0.25, 0.25, 0.25, 0.25};
    ELEMENT *e;

    assert(g != NULL);
    if (g->nroot <= 0) {
	phgWarning("no element!\n");
	return 0;
    }
    if (g->nroot <= 2) {
	return 0;
    }

    sfc = x = phgAlloc(g->nroot * sizeof(SFC));
    x_index = phgAlloc(g->nroot * sizeof(*x_index));
    if (g->nroot > 0 && (x == NULL || x_index == NULL)) {
	phgError(1, "memory allocation error (%s:%d)\n", __FILE__, __LINE__);
	return -1;
    }

    /* scan all elements 
     * calculate barycentric coordinate */
    for (i = 0; i < g->nroot; i++) {
	e = g->roots + i;
	x_index[i] = e->index;
	phgGeomLambda2XYZ(g, e, lam, &x0, &y0, &z0);
	x[i].co[0] = x0;	  /* x */
	x[i].co[1] = y0;	  /* y */
	x[i].co[2] = z0;	  /* z */
    } /* for */
    
    /* boundary box */
    xmin = g->bbox[0][0];
    ymin = g->bbox[0][1];
    zmin = g->bbox[0][2];
    xmax = g->bbox[1][0];
    ymax = g->bbox[1][1];
    zmax = g->bbox[1][2];

    /* scale bounding box  
     * such that x, y, z in [0, 1)^3 */
    etx = xmax - xmin;
    ety = ymax - ymin;
    etz = zmax - zmin;

    if (etx == 0.)
	etx = 1.;
    if (etx < ety)
	etx = ety;
    if (etx < etz)
	etx = etz;
    etx =  1. / etx;
    ety = etx;
    etz = etx;

    for (i = 0; i < g->nroot; i++) {
	x[i].co[0] = (x[i].co[0] - xmin) * etx;
	x[i].co[1] = (x[i].co[1] - ymin) * ety;
	x[i].co[2] = (x[i].co[2] - zmin) * etz;
    }
    /* generate Hilbert SFC */
    if (!phgSFCInvHilbert3D(x, g->nroot)) {
	phgWarning("error in generating Hilbert SFC (%s:%d)\n",
						__FILE__, __LINE__);
	phgFree(x);
	return -1;
    }

    ordering = phgAlloc(g->nroot * sizeof(*ordering));
    for (i = 0; i < g->nroot; i++)
	ordering[i] = i;
    qsort(ordering, g->nroot, sizeof(*ordering), comp_key);
    phgFree(x);
    /* reorder the indices */
    {
	INT *aux;
	INT kk;

	aux = (INT *)phgAlloc(g->nroot * sizeof(INT));
	if (aux == NULL) {
	    phgWarning("memory allocation error (%s:%d)\n", __FILE__, __LINE__);
	    phgFree(x);
	    return -1;
	}

	for (i = 0; i < g->nroot; i++) {
	    aux[x_index[ordering[i]]] = i;
	}

	for (i = 0; i < g->nroot; i++) {
	    e = g->roots + i;
	    kk = e->index;
	    e->index = aux[kk];
	}
	phgFree(aux);
    }

    phgFree(ordering);
    phgFree(x_index);

    return 0;
}

/* generate Morton Order */
/* 0: success, or fail */
int
phgMortonOrder(GRID *g)
{
    SFC *x;
    INT i, *x_index, *ordering;
    SFC_FLOAT xmax, xmin;
    SFC_FLOAT ymax, ymin;
    SFC_FLOAT zmax, zmin;
    SFC_FLOAT etx, ety, etz;
    FLOAT x0, y0, z0;
    static const FLOAT lam[] = {0.25, 0.25, 0.25, 0.25};
    ELEMENT *e;

    assert(g != NULL);
    if (g->nroot <= 0) {
	phgWarning("no element!\n");
	return 0;
    }
    if (g->nroot <= 2) {
	return 0;
    }

    sfc = x = phgAlloc(g->nroot * sizeof(SFC));
    x_index = phgAlloc(g->nroot * sizeof(*x_index));
    if (x == NULL && g->nroot > 0) {
	phgError(1, "memory allocation error (%s:%d)\n", __FILE__, __LINE__);
	return -1;
    }

    /* scan all elements 
     * calculate barycentric coordinate */
    for (i = 0; i < g->nroot; i++) {
	e = g->roots + i;
	x_index[i] = e->index;
	phgGeomLambda2XYZ(g, e, lam, &x0, &y0, &z0);
	x[i].co[0] = x0;	  /* x */
	x[i].co[1] = y0;	  /* y */
	x[i].co[2] = z0;	  /* z */
    } /* for */
    
    /* boundary box */
    xmin = g->bbox[0][0];
    ymin = g->bbox[0][1];
    zmin = g->bbox[0][2];
    xmax = g->bbox[1][0];
    ymax = g->bbox[1][1];
    zmax = g->bbox[1][2];

    /* scale bounding box 
     * such that x, y, z in [0, 1)^3 */
    etx = xmax - xmin;
    ety = ymax - ymin;
    etz = zmax - zmin;

    if (etx == 0.)
	etx = 1.;
    if (etx < ety)
	etx = ety;
    if (etx < etz)
	etx = etz;
    etx =  1. / etx;
    ety = etx;
    etz = etx;

    for (i = 0; i < g->nroot; i++) {
	x[i].co[0] = (x[i].co[0] - xmin) * etx;
	x[i].co[1] = (x[i].co[1] - ymin) * ety;
	x[i].co[2] = (x[i].co[2] - zmin) * etz;
    }
    /* generate Morton SFC */
    if (!phgSFCInvMorton3D(x, g->nroot)) {
	phgWarning("error in generating Morton SFC (%s:%d)\n",
						__FILE__, __LINE__);
	phgFree(x);
	return -1;
    }

    ordering = phgAlloc(g->nroot * sizeof(*ordering));
    for (i = 0; i < g->nroot; i++)
	ordering[i] = i;
    qsort(ordering, g->nroot, sizeof(*ordering), comp_key);
    phgFree(x);

    /* reorder the indices */
    {
	INT *aux;
	INT kk;

	aux = (INT *)phgAlloc(g->nroot * sizeof(INT));
	if (aux == NULL) {
	    phgWarning("memory allocation error (%s:%d)\n", __FILE__, __LINE__);
	    phgFree(x);
	    return -1;
	}

	for (i = 0; i < g->nroot; i++) {
	    aux[x_index[ordering[i]]] = i;
	}

	for (i = 0; i < g->nroot; i++) {
	    e = g->roots + i;
	    kk = e->index;
	    e->index = aux[kk];
	}
	phgFree(aux);
    }

    phgFree(ordering);
    phgFree(x_index);
    return 0;
}

/* construct initial grid order 
 * default: hilbert order */
static const char *orders[] = {"none", "hilbert", "morton", "hamilton", NULL};
static int GridOrder = 0;

int
phgGridInitOrder(GRID *g, BOOLEAN dist)
{
    int ret = 1;
    double t1[3], t2[3];
    static BOOLEAN initialized = FALSE;

    /* register parameters */
    if (!initialized) {
	initialized = TRUE;
	phgOptionsRegisterTitle("\n[RTK partitioner]", "\n", "partition");
	phgOptionsRegisterKeyword("-rtk_init_order", "Initial ordering",
		orders, &GridOrder);
	return TRUE;
    }

    assert(g != NULL);
    if (dist) {  /* parallel version */
	phgError(1, "Not implemented yet!\n");
    }
    else {	 /* sequential version */
	if (phgVerbosity > 0) {
	    phgPrintf("Initial ordering for RTK: ");
	    phgGetTime(t1);
	}
	
	/* Hilbert order */
	if (GridOrder == 1) {
	    if (phgVerbosity > 0) {
		phgPrintf("Hilbert space-filling curve.\n");
	    }
	    ret = phgHilbertOrder(g);
	} /* Morton order */
	else if (GridOrder == 2) {
	    if (phgVerbosity > 0) {
		phgPrintf("Morton space-filling curve.\n");
	    }
	    ret = phgMortonOrder(g);
	}    /* Hamilton Path */
	else if (GridOrder == 3) {
	    if (phgVerbosity > 0) {
		phgPrintf("\tHamilton path.\n");
	    }
	    ret = phgHamiltonPath(g);
	}
	else {
	    if (phgVerbosity > 0) {
		phgPrintf("none.\n");
	    }
	    ret = 0;
	}
	
	/* success */
	if (ret == 0) {
	    if (phgVerbosity > 0) {
		phgGetTime(t2);
		phgPrintf("Finished computing initial ordering of elements: ");
		phgPrintf("%fs\n", t2[2] - t1[2]);
		phgPrintf("Total roots: %"dFMT"\n", g->nroot);
	    }
	}
	else { /* fail */
	    if (phgVerbosity > 0) {
		phgPrintf("Error in computing initial ordering.");
	    }
	}
    }
    return ret;
}

#if USE_MPI
/***************************************************
 *	       1D partitioner 
 **************************************************/
typedef struct {
    SFC_FLOAT low;
    SFC_FLOAT high;
    int   pn;
}PTNS;

/* tolerance for load balance efficiency */
static FLOAT LB_TOL;

/* return the partition id that key belongs to */
static int bsearch_ptn(PTNS *ptn, int p, SFC_FLOAT key);

/* max loops */
static int MAXLOOPS = 10;

static INT CUTS = 16;

void phgPartition1DOptionsSet()
{
    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	initialized = TRUE;
	LB_TOL = SFC_LBE;
	phgOptionsRegisterTitle("\n[1D partitioner]", "\n", 
		"partition");
	phgOptionsRegisterFloat("-lbe_tol", "Tolerance of load balance "
		"efficiency (>= 1.)", &LB_TOL);
	phgOptionsRegisterInt("-1d_cut", "number of cuts of 1D partitioner:"
		"10 <= cut <= 30 ", &CUTS);

	return;
    }
}

/* cut [0, 1] to p parts, such that lbe <= LB_TOL 
 * dts: store key \in [0, 1]
 * lx: length of dts 
 * p: # of partitions */
FLOAT
phgPartition1DV1(DOTS *dts, INT lx, int p, const FLOAT speeds[], MPI_Comm comm)
{
    int npct = (p - 1) * CUTS + 1;  /* total  cuts */
    PTNS  ptn[npct];

    SFC_FLOAT lbe;		    /* load balance efficiency */
    SFC_FLOAT DotsSum;		    /* total sum */
    SFC_FLOAT sum[npct];
    SFC_FLOAT tsum[npct];

    /* arrays only needed on root */
    SFC_FLOAT *psum = NULL;	    /* partial sum, psum[i] = sf * i * sum */
    SFC_FLOAT (*bbox)[2] = NULL;    /* bounding box, 0: low, 1: high*/
    SFC_FLOAT *tmp = NULL;

    INT i,j;
    INT ctn;			    /* if we start next loop */
    int rank;			    /* pid */
    int flag;			    /* iterative times */
    int map[npct];		    /* partition result */

    /* check parameters */
    if (LB_TOL < 1.) {
	phgWarning("LB_TOL too small, changed to %lg\n", (double)SFC_LBE);
	LB_TOL = SFC_LBE;
    }
    else if (LB_TOL > 1.3) {
	phgWarning("LB_TOL too large (expecting poor load balance)\n");
    }

    assert((CUTS >= 10) && (CUTS <= 30));
    assert(lx >= 0 && p >= 1);
    assert(lx == 0 || dts != NULL);

    /* deal with special case */
    if (p == 1) {  /* only one partition */
	for (i = 0; i < lx; i++) {
	    dts[i].pn = 0;
	}
	return 1.;
    }
    
    /* compute sum of lx (bbox[0][0]) and weights (bbox[0][1]) */
    lbe = 0.;
    for (i = 0; i < lx; i++) {
	lbe += dts[i].w;
    }
    {
	SFC_FLOAT tmp[4];
	tmp[0] = lx;
	tmp[1] = lbe;
	MPI_Allreduce(tmp, tmp + 2, 2, MPI_SFC_FLOAT, MPI_SUM, comm);
	if (tmp[2] == 0.)
	    return 1.;	/* no dot, just return */
	DotsSum = tmp[3];
    }

    /* special case */
    if (DotsSum == 0.) {
	for (i = 0; i < lx; i++)
	    dts[i].pn = 0;
	return 1.;
    }
    else if (DotsSum < 0.) {
	phgError(1, "negative weights (%s:%d)\n", __FILE__, __LINE__);
    }

    MPI_Comm_rank(comm, &rank);

    /* initialize partial sum, ideal condition
     * initialize bounding boxs */

    if (rank == 0) {
	psum = phgAlloc(p * sizeof(*psum));
	bbox = phgAlloc(p * sizeof(*bbox));
	tmp = phgAlloc(p * sizeof(*tmp));
	lbe = 0.0;
	for (i = 0; i < p; i++) {
	    lbe += speeds[i];
	    psum[i] = DotsSum * lbe;
	    bbox[i][0] = 0.;
	    bbox[i][1] = 1.;
	}
    }

    /* initial guess partitions */
    lbe = 1. / npct;
    for (i = 0; i < npct; i++) {
	ptn[i].low = i * lbe;
	ptn[i].high = (i + 1) * lbe;
	ptn[i].pn = i;
    }

    /* partitioning 
     * rank 0 generates results then broadcast */
    flag = 1;
    ctn = 1;
    while (flag <= MAXLOOPS && ctn) {
	SFC_FLOAT mxm, s, d;
	INT ne[npct];  /* record the number of elements in each interval*/
	INT ne1[npct]; 
	int ct[npct];  /* record if current interval is overflow interval */
		       /* 0: no; > 0, yes and number of overflow intervals */
	int p1 = p - 1;
	int kk = 0;

	/* sum of current partitions */
	for (i = 0; i < npct; i++) {
	    tsum[i] = 0.;
	    ne1[i] = 0;
	}
	for (i = 0; i < lx; i++) {
	    ctn = bsearch_ptn(ptn, npct, dts[i].key);
	    dts[i].pn = ctn;
	    tsum[ctn] += dts[i].w;
	    ne1[ctn] += 1;
	}

	/* calculate total weights of all intervals */
	MPI_Reduce(tsum, sum, npct, MPI_SFC_FLOAT, MPI_SUM, 0, comm);

	/* calculate the number of elements in each interval */
	MPI_Reduce(ne1, ne, npct, PHG_MPI_INT, MPI_SUM, 0, comm);

	/* rank 0 generate partitions 
	 * then broadcast to all processes */
	if (rank == 0) {
	    /* prefix sum */
	    tsum[0] = sum[0];
	    ct[0] = 0;
	    for (i = 1; i < npct; i++) {
		tsum[i] = sum[i];
		sum[i] += sum[i - 1];
		ct[i] = 0;
	    }

	    /* update bounding box */
	    j = 0;
	    for (i = 0; i < p1; i++) {
		/* look for new bounding box position */
		while (sum[j] < psum[i]) {
		    j++;
		}
		/* update bounding box */
		if (ptn[j].high < bbox[i][1])
		    bbox[i][1] = ptn[j].high;
		if (ptn[j].low	> bbox[i][0])
		    bbox[i][0] = ptn[j].low;
		ct[j] += 1;

		tmp[i] = 0.;
	    }
	    tmp[p - 1] = 0.;

	    /* calculate max sum of current partition
	     * map map[npct] to p parts */
	    j = 0;
	    for (i = 0; i < npct && j < p; i++) {
		if (sum[i] < psum[j]) {
		    map[i] = j;
		}
		else {
		    map[i] = j;
		    j++;

		    /* check if each overflow invterval has only one
		     * element.
		     * kk = 0 if yes. kk = 1 if no */
		    if ((ne[i] > 1) && (j < p1)) {
			kk = 1;
		    }
		}
	    }
	    /* special case that j >= p */
	    for (j = i; j < npct; j++) {
		map[j] = p1;
	    }

	    /* calculate load balance efficiency */
	    for (i = 0; i < npct; i++) {
		j = map[i];
		tmp[j] += tsum[i];
	    }
	    s = mxm = tmp[0] / speeds[0];
	    for (i = 1; i < p; i++) {
		d = tmp[i] / speeds[i];
		if (mxm < d)
		    mxm = d;
		s += d;
	    }
	    lbe = mxm * p / s;
#if 0
	    printf("lbe = %f, kk = %d\n", (double)lbe, kk);
#endif

	    /* terminates if procedure satisfies stop criteria.
	     * two stop criterias:
	     * (#if 1) the first one terminates if lbe is less than or
	     * equals to LB_TOL or each overflow interval has
	     * only one element.
	     * (#else part)the second one just uses lbe */
#if 1
	    if ((lbe <= LB_TOL) || (kk == 0 && p >= 3)) {
		ctn = 0;
	    }
	    else {
		ctn = 1;
	    }
#else	/* 1 */
	    if (lbe <= LB_TOL) {
		ctn = 0;
	    }
	    else {
		ctn = 1;
	    }
#endif	/* 1 */

	    if (phgVerbosity > 2) {
		phgPrintf("phgPartition1DV1 - LIF: %f\n", (double)lbe);
	    }
	}

	/* synchronize status
	 * 1: start next iteration, 0: stop */
	MPI_Bcast(&ctn, 1, PHG_MPI_INT, 0, comm);
       
	if (ctn) {
	    ctn = npct - 1;

	    /* generate new partitions
	     * cut p - 1 bounding boxes to npct parts */
	    if (rank == 0) {
		SFC_FLOAT len;
		INT k, j1, j2;
		INT ctk;
		
		i = 0;
		j = 0;
		k = 0;
		while (i < p1) {
		    while (ct[j] < 1) {
			j++;
		    }

		    /* new partition 
		     * we use array sum[ctn] for temp variable 
		     * sum[ctn] stores internal partitioning points */
		    j2 = ct[j] * CUTS;
#if 0
		    if (j < ctn) {
			len = 1. / j2;

			ctk = i + ct[j] - 1;
			for (j1 = 1; j1 <= j2; j1++) {
			    sum[k] = (bbox[i][0] * (j2 - j1) 
				    + bbox[ctk][1] * j1) * len;
			    k++;
			}
		    }
		    else {
			len = 1. / (j2 + 1.);
			ctk = i + ct[j] - 1;
			for (j1 = 1; j1 <= j2; j1++) {
			    sum[k] = (bbox[i][0] * (j2 + 1 - j1) 
				    + bbox[ctk][1] * j1) * len;
			    k++;
			}

		    }
#else	/* 0 */
		    len = 1. / j2;
		    ctk = i + ct[j] - 1;

		    for (j1 = 1; j1 <= j2; j1++) {
			sum[k] = (bbox[i][0] * (j2 - j1) 
				+ bbox[ctk][1] * j1) * len;
			k++;
		    }
#endif	/* 0 */
		    i += ct[j];
		    j++;
		}
#if 0
		/* check */
		for (i = 1; i < ctn; i++) {
		    assert(sum[i - 1] <= sum[i]);
		}
#endif	/* 0 */
	    }

	    /* rank 0 broadcasts new partitions.
	     * each process initializes new partitions */
	    MPI_Bcast(sum, ctn, MPI_SFC_FLOAT, 0, comm);

	    ptn[0].low = 0.;
	    ptn[0].high = sum[0];
	    for (i = 1; i < ctn; i++) {
		ptn[i].high = sum[i];
		ptn[i].low = sum[i - 1];
	    }
	    ptn[ctn].high = 1.;
	    ptn[ctn].low = sum[ctn - 1];
	    
	    ctn = 1;
	}
	else {
	    break;
	}

	flag++;
    }  /* while */

    if (rank == 0) {
	phgFree(psum);
	psum = NULL;
	phgFree(bbox);
	bbox = NULL;
	phgFree(tmp);
	tmp = NULL;
    }

    /* synchronize results */
    MPI_Bcast(map, npct, PHG_MPI_INT, 0, comm);
    MPI_Bcast(&lbe, 1, MPI_SFC_FLOAT, 0, comm);

    for (i = 0; i < lx; i++) {
	j = dts[i].pn;
	dts[i].pn = map[j];
    }

    if (phgVerbosity > 0) {
	phgPrintf("phgPartition1DV1 - iterations: %d\n", flag);
	phgPrintf("phgPartition1DV1 - LIF: %lf\n", (double)lbe);
    }

    return lbe;
}

FLOAT
phgPartition1DV2(DOTS *dts, INT lx, int p, const FLOAT speeds[], MPI_Comm comm)
{
    int npct = (p - 1) * CUTS + 1;  /* total  cuts */
    PTNS  ptn[npct];

    SFC_FLOAT lbe;		    /* load balance efficiency */
    SFC_FLOAT DotsSum;		    /* total sum */
    SFC_FLOAT sum[npct];
    SFC_FLOAT tsum[npct];

    /* arrays only needed on root */
    SFC_FLOAT *psum = NULL;	    /* partial sum, psum[i] = sf * i * sum */
    SFC_FLOAT (*bbox)[2] = NULL;    /* bounding box, 0: low, 1: high*/
    SFC_FLOAT *tmp = NULL;

    INT i,j;
    INT ctn;			    /* if we start next loop */
    int rank;			    /* pid */
    int flag;			    /* iterative times */
    int map[npct];		    /* partition result */
    int unpt;

    /* check parameters */
    if (LB_TOL < 1.) {
	phgWarning("LB_TOL too small, changed to %lg\n", (double)SFC_LBE);
	LB_TOL = SFC_LBE;
    }
    else if (LB_TOL > 1.3) {
	phgWarning("LB_TOL too large (expecting poor load balance)\n");
    }

    assert((CUTS >= 10) && (CUTS <= 30));
    assert(lx >= 0 && p >= 1);
    assert(lx == 0 || dts != NULL);

    /* deal with special case */
    if (p == 1) {  /* only one partition */
	for (i = 0; i < lx; i++) {
	    dts[i].pn = 0;
	}
	return 1.;
    }
    
    /* compute sum of lx (bbox[0][0]) and weights (bbox[0][1]) */
    lbe = 0.;
    for (i = 0; i < lx; i++) {
	lbe += dts[i].w;
    }
    {
	SFC_FLOAT tmp[4];
	tmp[0] = lx;
	tmp[1] = lbe;
	MPI_Allreduce(tmp, tmp + 2, 2, MPI_SFC_FLOAT, MPI_SUM, comm);
	if (tmp[2] == 0.)
	    return 1.;	/* no dot, just return */
	DotsSum = tmp[3];
    }

    /* special case */
    if (DotsSum == 0.) {
	for (i = 0; i < lx; i++)
	    dts[i].pn = 0;
	return 1.;
    }
    else if (DotsSum < 0.) {
	phgError(1, "negative weights (%s:%d)\n", __FILE__, __LINE__);
    }

    MPI_Comm_rank(comm, &rank);

    /* initialize partial sum, ideal condition
     * initialize bounding boxs */

    if (rank == 0) {
	psum = phgAlloc(p * sizeof(*psum));
	bbox = phgAlloc(p * sizeof(*bbox));
	tmp = phgAlloc(p * sizeof(*tmp));
	lbe = 0.0;
	for (i = 0; i < p; i++) {
	    lbe += speeds[i];
	    psum[i] = DotsSum * lbe;
	    bbox[i][0] = 0.;
	    bbox[i][1] = 1.;
	}
    }

    /* initial guess partitions */
    lbe = 1. / npct;
    for (i = 0; i < npct; i++) {
	ptn[i].low = i * lbe;
	ptn[i].high = (i + 1) * lbe;
	ptn[i].pn = i;
    }

    DotsSum = p / DotsSum;

    /* partitioning 
     * rank 0 generates results then broadcast */
    flag = 1;
    ctn = 1;
    unpt = npct;
    while (flag <= MAXLOOPS && ctn) {
	SFC_FLOAT mxm, s, d;
	INT ne[npct];  /* record the number of elements in each interval*/
	INT ne1[npct]; 
	int ct[npct];  /* record if current interval is overflow interval */
		       /* 0: no; > 0, yes and number of overflow intervals */
	int p1 = p - 1;
	int kk = 0;
	int skc;

	/* sum of current partitions */
	for (i = 0; i < npct; i++) {
	    tsum[i] = 0.;
	    ne1[i] = 0;
	}
	for (i = 0; i < lx; i++) {
	    ctn = bsearch_ptn(ptn, unpt, dts[i].key);
	    dts[i].pn = ctn;
	    tsum[ctn] += dts[i].w;
	    ne1[ctn] += 1;
	}

	/* calculate total weights of all intervals */
	MPI_Reduce(tsum, sum, npct, MPI_SFC_FLOAT, MPI_SUM, 0, comm);

	/* calculate the number of elements in each interval */
	MPI_Reduce(ne1, ne, npct, PHG_MPI_INT, MPI_SUM, 0, comm);

	/* rank 0 generate partitions 
	 * then broadcast to all processes */
	if (rank == 0) {
	    /* prefix sum */
	    tsum[0] = sum[0];
	    ct[0] = 0;
	    for (i = 1; i < npct; i++) {
		tsum[i] = sum[i];
		sum[i] += sum[i - 1];
		ct[i] = 0;
	    }

	    /* update bounding box */
	    j = 0;
	    for (i = 0; i < p1; i++) {
		/* look for new bounding box position */
		while (sum[j] < psum[i]) {
		    j++;
		}
		/* update bounding box */
		if (ptn[j].high < bbox[i][1])
		    bbox[i][1] = ptn[j].high;
		if (ptn[j].low	> bbox[i][0])
		    bbox[i][0] = ptn[j].low;
		ct[j] += 1;

		tmp[i] = 0.;
	    }
	    tmp[p - 1] = 0.;

	    /* calculate max sum of current partition
	     * map map[npct] to p parts */
	    j = 0;
	    for (i = 0; i < npct && j < p; i++) {
		if (sum[i] < psum[j]) {
		    map[i] = j;
		}
		else {
		    map[i] = j;
		    j++;

		    /* check if each overflow invterval has only one
		     * element.
		     * kk = 0 if yes. kk = 1 if no */
		    if ((ne[i] > 1) && (j < p1)) {
			kk = 1;
		    }
		}
	    }
	    /* special case that j >= p */
	    for (j = i; j < npct; j++) {
		map[j] = p1;
	    }

	    /* calculate load balance efficiency */
	    for (i = 0; i < npct; i++) {
		j = map[i];
		tmp[j] += tsum[i];
	    }
	    s = mxm = tmp[0] / speeds[0];
	    for (i = 1; i < p; i++) {
		d = tmp[i] / speeds[i];
		if (mxm < d)
		    mxm = d;
		s += d;
	    }
	    lbe = mxm * p / s;
#if 0
	    printf("lbe = %f, kk = %d\n", (double)lbe, kk);
#endif

	    /* terminates if procedure satisfies stop criteria.
	     * two stop criterias:
	     * (#if 1) the first one terminates if lbe is less than or
	     * equals to LB_TOL or each overflow interval has
	     * only one element.
	     * (#else part)the second one just uses lbe */
#if 1
	    if ((lbe <= LB_TOL) || (kk == 0 && p >= 3)) {
		ctn = 0;
	    }
	    else {
		ctn = 1;
	    }
#else	/* 1 */
	    if (lbe <= LB_TOL) {
		ctn = 0;
	    }
	    else {
		ctn = 1;
	    }
#endif	/* 1 */

	    if (phgVerbosity > 2) {
		phgPrintf("phgPartition1DV2 - LIF: %f\n", (double)lbe);
	    }
	}

	/* synchronize status
	 * 1: start next iteration, 0: stop */
	MPI_Bcast(&ctn, 1, PHG_MPI_INT, 0, comm);
       
	/* calculate the number of elements in each bounding box 
	 * if all the numbers are 1, then terminate */
	if (ctn) {
	    /* generate new partitions
	     * cut p - 1 bounding boxes to npct parts */
	    ctn = npct - 1;
	    if (rank == 0) {
		SFC_FLOAT len;
		SFC_FLOAT md;	 /* mid point of some partition */
		INT k, j1, j2;
		INT ctk;
#if 0
		SFC_FLOAT lbe1, lbe2;

		kk = 0;
		for (i = 1; i < p1; i++) {
		    if (bbox[i - 1][1] > bbox[i][0])
			kk++;
		}

		i = 0;
		j = 0;
		k = 0;
		skc = 0;
		while (i < p1) {
		    while (ct[j] < 1) {
			j++;
		    }

		    if (ct[j] == 1) {
			lbe1 = sum[j] * DotsSum / (i + 1.);
			lbe2 = tmp[i] * DotsSum;
		    }
		    else {
			lbe1 = LB_TOL + 0.1;
			lbe2 = LB_TOL + 0.1;
		    }
		    /* new partition 
		     * we use array tsum[ctn] for temp variable 
		     * tsum[ctn] stores internal partitioning points */
		    md = (bbox[i][0] + bbox[i][1]) * 0.5;

		    if ( ((lbe1 <= LB_TOL) && (lbe2 <= LB_TOL) && (lbe2 >= 1.0))
		      || ((ct[j] == 1) && (kk == 0) && (ne[j] <= 1)) 
		      || (md == bbox[i][0] || bbox[i][1] == md)){
			tsum[k] = bbox[i][0];
			k += 1;
			tsum[k] = md;
			k += 1;
			tsum[k] = bbox[i][1];
			k += 1;

			i += ct[j];
			j++;
			skc++;
		    }
		    else {
			j2 = ct[j] * CUTS;
			len = 1. / j2;
			ctk = i + ct[j] - 1;

			for (j1 = 1; j1 <= j2; j1++) {
			    tsum[k] = (bbox[i][0] * (j2 - j1) 
				    + bbox[ctk][1] * j1) * len;
			    k++;
			}
			i += ct[j];
			j++;
		    }
		}
#else	/* 0 */
		kk = 0;
		for (i = 1; i < p1; i++) {
		    if (bbox[i - 1][1] > bbox[i][0])
			kk++;
		}

		i = 0;
		j = 0;
		k = 0;
		skc = 0;
		while (i < p1) {
		    while (ct[j] < 1) {
			j++;
		    }

		    /* new partition 
		     * we use array tsum[ctn] for temp variable 
		     * tsum[ctn] stores internal partitioning points */
		    md = (bbox[i][0] + bbox[i][1]) * 0.5;

		    if (((ct[j] == 1) && (kk == 0) && (ne[j] <= 1)) 
		      || (md == bbox[i][0] || bbox[i][1] == md)){
			tsum[k] = bbox[i][0];
			k += 1;
			tsum[k] = md;
			k += 1;
			tsum[k] = bbox[i][1];
			k += 1;

			i += ct[j];
			j++;
			skc++;
		    }
		    else {
			j2 = ct[j] * CUTS;
			len = 1. / j2;
			ctk = i + ct[j] - 1;

			for (j1 = 1; j1 <= j2; j1++) {
			    tsum[k] = (bbox[i][0] * (j2 - j1) 
				    + bbox[ctk][1] * j1) * len;
			    k++;
			}
			i += ct[j];
			j++;
		    }
		}
#endif	/* 0 */
		
		unpt = k + 1;
		for (i = k; i <= ctn; i++) {
		    tsum[i] = 1.;
		}
#if 0
		printf("skc: %d k: %"dFMT" npct: %d \n", skc, k, npct);
		/* check */
		for (i = 1; i < ctn; i++) {
		    assert(tsum[i - 1] <= tsum[i]);
		}
#endif	/* 0 */
	    }

	    /* rank 0 broadcasts new partitions.
	     * each process initializes new partitions */
	    MPI_Bcast(tsum, ctn, MPI_SFC_FLOAT, 0, comm);
	    MPI_Bcast(&skc, 1, MPI_INT, 0, comm);
	    MPI_Bcast(&unpt, 1, MPI_INT, 0, comm);

	    if (skc < p1) {
		ptn[0].low = 0.;
		ptn[0].high = tsum[0];
		for (i = 1; i < ctn; i++) {
		    ptn[i].high = tsum[i];
		    ptn[i].low = tsum[i - 1];
		}
		ptn[ctn].high = 1.;
		ptn[ctn].low = tsum[ctn - 1];

		ctn = 1;
	    }
	    else {
		ctn = 0;
	    }
	}
	else {
	    break;
	}

	flag++;
    }  /* while */

    if (rank == 0) {
	phgFree(psum);
	psum = NULL;
	phgFree(bbox);
	bbox = NULL;
	phgFree(tmp);
	tmp = NULL;
    }

    /* synchronize results */
    MPI_Bcast(map, npct, PHG_MPI_INT, 0, comm);
    MPI_Bcast(&lbe, 1, MPI_SFC_FLOAT, 0, comm);

    for (i = 0; i < lx; i++) {
	j = dts[i].pn;
	dts[i].pn = map[j];
    }

    if (phgVerbosity > 0) {
	phgPrintf("phgPartition1DV2 - iterations: %d\n", flag);
	phgPrintf("phgPartition1DV2 - LIF: %lf\n", (double)lbe);
    }

    return lbe;
}

static int
bsearch_ptn(PTNS *ptn, int p, SFC_FLOAT key)
{
    int lo, hi;
    int md;
    int ret;
    
    assert(p >= 1);
    
    lo = 0;
    hi = p - 1;
    md = (lo + hi) >> 1;
    while (1) {
	if (key >= ptn[md].high)
	    ret = 1;
	else if (key - ptn[md].low < 0.)
	    ret = -1;
	else
	    ret = 0;

	if (ret == 1) {
	    lo = md;
	    md = (lo + hi) >> 1;
	    if (md == lo)
		md = hi;
	}
	else if (ret == -1) {
	    hi = md;
	    md = (lo + hi) >> 1;
	}
	else {
	    break;
	}
    }
    return md;
}

/* phgPartitionRemap remaps the partitions such that the
 * the mount of migrated data will be minimized.
 * newcomm: the new communicator (new number of partitions)
 * return value > 0, if partitions are changed;
 * return value = 0, if partitions arn't changed;
 * return value < 0, if subroutines fail and then partitions keep the same.
 *
 * The algorithm comes from PLUM: 
 * parallel load balancing for adaptive unstructured meshed 
 * Leonid Oliker, Rupak Biswas 
 *
 * heuristic algorithm 
 * partition_remap_v1() stores the whole similarity matrix.
 * partition_remap_v2() is an improved version which stores the nonzero entries
 * only.
 * The second version can use memories efficiently,
 * and will be more competitive when the number of processes
 * are large */

/* similiar matrix */
typedef struct {
    double s;	/* number of elements */
    int pid;	/* current pid */
    int ptn;	/* partition id */
} SM;

/* descending order used by qsort */
static int
comp_sml(const void *m1, const void *m2)
{
    SM *n1, *n2;
    int flag = 0;
    
    n1 = (SM *)m1;
    n2 = (SM *)m2;

    if (n1->s > n2->s)
	flag = -1;
    else if (n1->s == n2->s)
	flag = 0;
    else if (n1->s < n2->s)
	flag = 1;
    return flag;
}

static int
partition_remap_v1(GRID *g, MPI_Comm newcomm)
/* Version 1 of the algorithm for remapping partitions */
{
    SM *local, *sml = NULL;
    ELEMENT *e;
    MPI_Datatype type;
    MPI_Datatype types[3];
    MPI_Aint indices[3];
    INT i, j;
    int *result;
    int rank, size, nprocs;
    int blocklens[3];
    int ret = 0;

    MPI_Comm_size(newcomm, &nprocs);
    MPI_Comm_rank(newcomm, &rank);
    if (g == NULL || rank >= g->nprocs)
	return ret;

    i = g->nleaf;
    MPI_Allreduce(&i, &j, 1, PHG_MPI_INT, MPI_SUM, g->comm);
    if ((j == 0) || (nprocs <= 1)) {
	return ret;
    }

    rank = g->rank;
    size = g->nprocs;

    assert(nprocs >= size);

    local = (SM *) phgAlloc(nprocs * sizeof(SM));
    result = (int *) phgAlloc(nprocs * sizeof(int));

    if (local == NULL || result == NULL) {
	phgError(1, "memory allocation error (%s:%d)\n", __FILE__, __LINE__);
    }

    /* create a new datatype for exchange information */
    /* used by structure tree */
    blocklens[0] = blocklens[1] = blocklens[2] = 1;
    types[0] = MPI_DOUBLE;
    types[1] = types[2] = MPI_INT;
#if 1	/* MPI-2 */
    MPI_Get_address(&local[0].s, &indices[0]);
    MPI_Get_address(&local[0].pid, &indices[1]);
    MPI_Get_address(&local[0].ptn, &indices[2]);
#else	/* MPI-1 */
    MPI_Address(&local[0].s, &indices[0]);
    MPI_Address(&local[0].pid, &indices[1]);
    MPI_Address(&local[0].ptn, &indices[2]);
#endif
    indices[2] -= indices[0];
    indices[1] -= indices[0];
    indices[0] = 0;
#if 1	/* MPI-2 */
    MPI_Type_create_struct(3, blocklens, indices, types, &type);
#else	/* MPI-1 */
    MPI_Type_struct(3, blocklens, indices, types, &type);
#endif
    MPI_Type_commit(&type);


    /* gather local information to creat similiar matrix */
    for (i = 0; i < nprocs; i++) {
	local[i].s = 0;
	local[i].pid = rank;
	local[i].ptn = i;
    }

    ForAllElements(g, e) {
	local[e->mark].s += 1.;
    }

    if (rank == 0) {
	/* allocate space for process 0 */
	sml = (SM *) phgAlloc(nprocs * nprocs * sizeof(SM));
	if (sml == NULL) {
	    phgError(1, "memory allocation error (%s:%d)\n", __FILE__,__LINE__);
	}
    }

    /* root: build similiar matrix */
    MPI_Gather(local, nprocs, type, sml, nprocs, type, 0, g->comm);
    MPI_Type_free(&type);
    phgFree(local);

    /* map algorithm by PLUM */
    if (rank == 0){
	int *part_map;
	int *proc_unmap;
	int *t1;
	int count;
	SM *p;

	part_map = (int *) phgAlloc(nprocs * sizeof(int));
	proc_unmap = (int *) phgAlloc(nprocs * sizeof(int));
	t1 = (int *) phgAlloc(nprocs * sizeof(int));

	if ((part_map == NULL) || (proc_unmap == NULL) || (t1 == NULL)) {
	    phgError(1, "memory allocation error (%s:%d)\n", __FILE__,__LINE__);
	}

	/* inlarge similiar matrix */
	/* such that sml is square */
	if (size < nprocs) {
	    for (i = size; i < nprocs; i++) {
		for (j = 0; j < nprocs; j++) {
		    count = i * nprocs + j;
		    sml[count].pid = i;
		    sml[count].ptn = j;
		    sml[count].s = 0;
		}
	    }
	}

	/* initialize */
	for (i = 0; i < nprocs; i++) {
	    proc_unmap[i] = 1;
	    part_map[i] = 0;
	}
	/* order the similiar matrix, descending order */
	qsort(sml, size * nprocs, sizeof(SM), comp_sml);

	/* map algorithm */
	count = 0;
	p = sml;
	while (count < nprocs) {
	    while ((proc_unmap[p->pid] == 0 || part_map[p->ptn] == 1)) {
		p++;
	    }
	    proc_unmap[p->pid] = 0;
	    part_map[p->ptn] = 1;
	    result[p->ptn] = p->pid;
	    count++;
	} /* map algorithm */

	/* verification 
	 * count = 0: valid result 
	 * count > 0: invalid result.*/
	count = 0;
	for (i = 0; i < nprocs; i++) {
	    t1[i] = 0;

	    /* check if final partitions have proper indices.
	     * count = 0 if yes, > 0 if not */
	    if ((result[i] < 0) || (result[i] > nprocs)) {
		count++;
	    }
	}

	/* check if each index was mapped just once 
	 * count = 0 if yes, > 0 if not */
	if (count == 0) {
	    for (i = 0; i < nprocs; i++) {
		t1[result[i]] += 1;
	    }
	    for (i = 0; i < nprocs; i++) {
		if (t1[i] != 1)
		    count++;
	    }
	}

	if (count == 0) { /* succeed */
	    /* default: partitions unchanged. */
	    ret = 0;

	    /* check if partitions are changed */
	    for (i = 0; i < nprocs; i++) {
		if (result[i] != i)
		    ret = 1;
	    }
	}
	else { /* fail */
	    ret = -1;
	}


	phgFree(part_map);
	phgFree(proc_unmap);
	phgFree(t1);
    }

    /* broadcast results */
    MPI_Bcast(&ret, 1, MPI_INT, 0, g->comm);
    if (ret > 0) {
	MPI_Bcast(result, nprocs, MPI_INT, 0, g->comm);
	ForAllElements(g, e) {
	    i = e->mark;
	    e->mark = result[i];
	}
	if (phgVerbosity > 2) {
	    phgPrintf("%s - map results\n", __func__);
	    for (i = 0; i < nprocs; i++) {
		phgPrintf("map: %6"dFMT"  ---> %6d\n", i, result[i]);
	    }
	}
    }
    else {
	if (phgVerbosity > 0)
	    phgPrintf("%s - partitions unchanged\n", __func__);
    }

    phgFree(sml);
    phgFree(result);

    return ret;
}

static int
partition_remap_v2(GRID *g, MPI_Comm newcomm)
/* Version 2 of the algorithm for remapping partitions */
{
#if 1

    ELEMENT *e;
    double *count;
    int *perm;
    int i, rank, nprocs;
    int ret;

    MPI_Comm_size(newcomm, &nprocs);
    MPI_Comm_rank(newcomm, &rank);
    if (nprocs <= 1)
	return 0;

    count = phgCalloc(nprocs, sizeof(*count));
    perm = phgCalloc(nprocs, sizeof(*perm));

    if (count == NULL || perm == NULL) {
	phgError(1, "memory allocation error (%s:%d)\n", __FILE__, __LINE__);
    }

    if (rank >= g->nprocs)
	goto label;

    /* calculating communication data size */
    ForAllElements(g, e) {
	i = e->mark;
	assert(i >= 0 && i < nprocs);
	count[i] += 1.;
    }

    /* remap */
label:
    ret = phgPartitionRemap0(g->comm, count, perm, newcomm);

    if (rank < g->nprocs && ret > 0)
	ForAllElements(g, e)
	    e->mark = perm[e->mark];

    phgFree(count);
    phgFree(perm);

    return ret;
#else	/* 1 */

    SM *local, *sml = NULL;
    ELEMENT *e;
    MPI_Datatype type;
    MPI_Datatype types[3];
    MPI_Aint indices[3];
    INT i;
    INT nt;
    int *result;
    int rank, size;
    int nonz;
    int blocklens[3];
    int *len;
    int *ofst;
    int ret = 0;

    i = g->nleaf;
    MPI_Allreduce(&i, &nt, 1, PHG_MPI_INT, MPI_SUM, g->comm);
    if ((nt == 0) || (nprocs <= 1)) {
	return ret;
    }

    local = (SM *) phgAlloc(nprocs * sizeof(SM));
    result = (int *) phgAlloc(nprocs * sizeof(int));
    len = (int *)phgAlloc(nprocs * sizeof(int));
    ofst = (int *)phgAllooc(nprocs * sizeof(int));

    if (local == NULL || result == NULL || len == NULL || ofst == NULL) {
	phgError(1, "memory allocation error (%s:%d)\n", __FILE__, __LINE__);
    }

    /* create a new datatype for exchange information */
    /* used by structure tree */
    blocklens[0] = blocklens[1] = blocklens[2] = 1;
    types[0] = MPI_DOUBLE;
    types[1] = types[2] = MPI_INT;
#if 1	/* MPI-2 */
    MPI_Get_address(&local[0].s, &indices[0]);
    MPI_Get_address(&local[0].pid, &indices[1]);
    MPI_Get_address(&local[0].ptn, &indices[2]);
#else	/* MPI-1 */
    MPI_Address(&local[0].s, &indices[0]);
    MPI_Address(&local[0].pid, &indices[1]);
    MPI_Address(&local[0].ptn, &indices[2]);
#endif
    indices[2] -= indices[0];
    indices[1] -= indices[0];
    indices[0] = 0;
#if 1	/* MPI-2 */
    MPI_Type_create_struct(3, blocklens, indices, types, &type);
#else	/* MPI-1 */
    MPI_Type_struct(3, blocklens, indices, types, &type);
#endif
    MPI_Type_commit(&type);

    rank = g->rank;
    size = g->nprocs;

    /* gather information to creat similiar matrix */
    for (i = 0; i < nprocs; i++) {
	local[i].s = 0.;
	local[i].pid = rank;
	local[i].ptn = i;
    }
    ForAllElements(g, e) {
	local[e->mark].s += 1.;
    }

    /* calculate non-zero entries of each process.
     * and move non-zero entries to the front part of the array */
    nonz = 0;
    for (i = 0; i < nprocs; i++) {
	if (local[i].s > 0) {
	    local[nonz].s = local[i].s;
	    local[nonz].pid = local[i].pid;
	    local[nonz].ptn = local[i].ptn;
	    nonz++;
	}
    }
    /* process 0 gathers the number of non-zero entries from 
     * all processes and calculate offset  */
    MPI_Gather(&nonz, 1, MPI_INT, len, 1, MPI_INT, 0, g->comm);

    if (rank == 0) {
	ofst[0] = 0;
	for (i = 1; i < size; i++) {
	    ofst[i] = ofst[i - 1] + len[i - 1];
	}

	/* allocate space for process 0 */
	nt = ofst[size - 1] + len[size - 1];
	sml = (SM *) phgAlloc(nt * sizeof(SM));
	if (sml == NULL) {
	    phgError(1, "memory allocation error (%s:%d)\n", __FILE__,__LINE__);
	}
    }

    /* process 0 build similiar matrix */
    MPI_Gatherv(local, nonz, type, sml, len, ofst, type, 0, g->comm);
    
    MPI_Type_free(&type);
    phgFree(local);
    phgFree(len);
    phgFree(ofst);

    /* map algorithm by PLUM */
    if (rank == 0){
	int *part_map;
	int *proc_unmap;
	int *t1;
	int count;
	int kk;
	SM *p;

	part_map = (int *) phgAlloc(nprocs * sizeof(int));
	proc_unmap = (int *) phgAlloc(nprocs * sizeof(int));
	t1 = (int *) phgAlloc(nprocs * sizeof(int));

	if ((part_map == NULL) || (proc_unmap == NULL) || (t1 == NULL)) {
	    phgError(1, "memory allocation error (%s:%d)\n", __FILE__,__LINE__);
	}

	/* initialize */
	for (i = 0; i < nprocs; i++) {
	    proc_unmap[i] = 1;
	    part_map[i] = 0;
	    result[i] = -1;
	}
	/* order the similiar matrix, descending order */
	qsort(sml, nt, sizeof(SM), comp_sml);

	/* map algorithm */
	count = 0;
	kk = 0;
	p = sml;
	nt -= 1;
	while ((count < nprocs) && (kk < nt)) {
	    while ( (proc_unmap[p->pid] == 0 || part_map[p->ptn] == 1)
		 && (kk < nt)) {
		p++;
		kk++;
	    }
	    if ((proc_unmap[p->pid] != 0) && (part_map[p->ptn] != 1)) {
		proc_unmap[p->pid] = 0;
		part_map[p->ptn] = 1;
		result[p->ptn] = p->pid;
		count++;
	    }
	} /* map algorithm */

	count = 0;
	for (i = 0; i < nprocs; i++) {
	    if (result[i] < 0)
		count++;
	    t1[i] = 0;
	}

	/* deal with special case that some arn't mapped */
	if (count > 0) {
	    for (i = 0; i < nprocs; i++) {
		if (result[i] >= 0) {
		    t1[result[i]] = -1;
		}
	    }
	    kk = 0;
	    for (i = 0; i < nprocs; i++) {
		if (result[i] < 0) {
		    while (t1[kk] < 0) {
			kk++;
		    }
		    result[i] = kk;
		    kk++;
		}
	    }
	}

	/* verification */
	count = 0;
	for (i = 0; i < nprocs; i++) {
	    t1[i] = 0;

	    /* check if final partitions have proper indices.
	     * count = 0 if yes, > 0 if not */
	    if ((result[i] < 0) || (result[i] > nprocs)) {
		count++;
	    }
	}

	/* check if each index was mapped just once 
	 * count = 0 if yes, > 0 if not */
	if (count == 0) {
	    for (i = 0; i < nprocs; i++) {
		t1[result[i]] += 1;
	    }
	    for (i = 0; i < nprocs; i++) {
		if (t1[i] != 1)
		    count++;
	    }
	}

	if (count == 0) { /* succeed */
	    /* default: partitions unchanged. */
	    ret = 0;

	    /* check if partitions are changed */
	    for (i = 0; i < nprocs; i++) {
		if (result[i] != i)
		    ret = 1;
	    }
	}
	else { /* fail */
	    ret = -1;
	}

	phgFree(part_map);
	phgFree(proc_unmap);
	phgFree(t1);
    }

    /* broadcast results */
    MPI_Bcast(&ret, 1, MPI_INT, 0, g->comm);
    if (ret > 0) {
	MPI_Bcast(result, nprocs, MPI_INT, 0, g->comm);
	ForAllElements(g, e) {
	    i = e->mark;
	    e->mark = result[i];
	}
	if (phgVerbosity > 2) {
	    phgPrintf("%s - map results\n", __func__);
	    for (i = 0; i < nprocs; i++) {
		phgPrintf("map: %6"dFMT"  ---> %6d\n", i, result[i]);
	    }
	}
    }
    else {
	if (phgVerbosity > 0)
	    phgPrintf("%s - partitions unchanged\n\n", __func__);
    }

    phgFree(sml);
    phgFree(result);

    return ret;

#endif	/* 1 */
}

int
phgPartitionRemap(GRID *g, MPI_Comm newcomm)
{
    int ret;

    /* first try version 2, and fall back to version 1 if the former fails */
    ret = partition_remap_v2(g, newcomm);
    if (ret < 0) {
	phgPrintf("WARNING - %s: Version 2 failed, try Version 1.\n", __func__);
	ret = partition_remap_v1(g, newcomm);
    }

    return ret;
}

/*
 * This function computes a permutation of {0,1,2,...,nprocs-1} in perm[]
 * which tends to minimize:
 *	\max_{rank in comm}
 *	    \sum_{i>=0 && i<nprocs && perm[i]!=rank}
 *		datasize_{rank}[i]
 * where datasize_{rank}[] represents the datasize[] array on process rank.
 *
 * datasize[] (input) and perm[] (output) are arrays of size nprocs,
 * and the resulting perm[] array is identical all across comm.
 */

static int
partition_remap0(MPI_Comm comm, const double datasize[], int perm[],
		 MPI_Comm newcomm)
/* Note: either ('comm' \subset 'newcomm') or ('newcomm' \subset 'comm'),
 * and this function can be called by all processes of either 'comm' or the
 * larger of 'comm' and 'newcomm' */
{
    SM *local, *sml = NULL;
    MPI_Datatype type;
    MPI_Datatype types[3];
    MPI_Aint indices[3];
    INT i;
    INT nt = 0;
    int rank, size, nprocs;
    int nonz;
    int blocklens[3];
    int *len;
    int *ofst;
    int ret = 0;

    MPI_Comm_size(newcomm, &nprocs);
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(newcomm, &rank);
#if DEBUG
    if (rank < size) {
	int rank0;
	MPI_Comm_rank(comm, &rank0);
	assert(rank == rank0);
    }
#endif	/* DEBUG */

    if (nprocs > 1 && !phgSameProcessorSpeed()) {
	/* make sure all processes have the same speed. */
	int flag = 0;
	if (nprocs < phgNProcs) {
	    double s[6];
	    s[0] = s[1] = s[2] = phgGetProcessorSpeed(NULL);
	    MPI_Reduce(s, s + 3, 1, MPI_3DOUBLE, MPI_MSM, 0, newcomm);
	    if (rank == 0)
		flag = ((s[5] - s[3]) / s[5] <= 1e-3 ? 1 : 0);
	    MPI_Bcast(&flag, 1, MPI_INT, 0, newcomm);
	}
	if (!flag) {
	    phgPrintf("SFC - Warning: different proc. speeds, no remapping.\n");
	    return ret;
	}
    }

    if (rank >= nprocs)
	return ret;

    assert((nprocs >= 1) && (datasize != NULL) && (perm != NULL));
    if (nprocs == 1) {
	perm[0] = 0;
	return ret;
    }

#if DEBUG
    /* communication datasize should not less than 0 */
    for (i = 0; i < nprocs; i++) {
	assert(datasize[i] >= 0.);
    }
#endif	/* DEBUG */

    local = (SM *) phgAlloc(nprocs * sizeof(SM));

    if (local == NULL) {
	phgError(1, "memory allocation error (%s:%d)\n", __FILE__, __LINE__);
    }

    /* create a new datatype for exchange information */
    /* used by structure tree */
    blocklens[0] = blocklens[1] = blocklens[2] = 1;
    types[0] = MPI_DOUBLE;
    types[1] = types[2] = MPI_INT;
#if 1	/* MPI-2 */
    MPI_Get_address(&local[0].s, &indices[0]);
    MPI_Get_address(&local[0].pid, &indices[1]);
    MPI_Get_address(&local[0].ptn, &indices[2]);
#else	/* MPI-1 */
    MPI_Address(&local[0].s, &indices[0]);
    MPI_Address(&local[0].pid, &indices[1]);
    MPI_Address(&local[0].ptn, &indices[2]);
#endif
    indices[2] -= indices[0];
    indices[1] -= indices[0];
    indices[0] = 0;
#if 1	/* MPI-2 */
    MPI_Type_create_struct(3, blocklens, indices, types, &type);
#else	/* MPI-1 */
    MPI_Type_struct(3, blocklens, indices, types, &type);
#endif
    MPI_Type_commit(&type);

    /* assemble local similiar matrix */
    for (i = 0; i < nprocs; i++) {
	local[i].s = datasize[i];
	local[i].pid = rank;
	local[i].ptn = i;
    }

    /* calculate non-zero entries, whose rank is less than nprocs.
     * and move non-zero entries to the front part of the array */
    nonz = 0;
    if (rank < nprocs) {
	for (i = 0; i < nprocs; i++) {
	    if (local[i].s > 0) {
		local[nonz].s = local[i].s;
		local[nonz].pid = local[i].pid;
		local[nonz].ptn = local[i].ptn;
		nonz++;
	    }
	}
    }

    /* process 0 gathers the number of non-zero entries from 
     * all processes and calculate offset  */
    len = (int *)phgAlloc(size * sizeof(int));
    ofst = (int *)phgAlloc(size * sizeof(int));

    if ((len == NULL) || (ofst == NULL)) {
	phgError(1, "memory allocation error (%s: %d)\n", __FILE__, __LINE__);
    }
    MPI_Gather(&nonz, 1, MPI_INT, len, 1, MPI_INT, 0, comm);

    if (rank == 0) {
	ofst[0] = 0;
	for (i = 1; i < size; i++) {
	    ofst[i] = ofst[i - 1] + len[i - 1];
	}

	/* allocate space for process 0 */
	nt = ofst[size - 1] + len[size - 1];
	sml = (SM *) phgAlloc(nt * sizeof(SM));
	if (nt > 0 && sml == NULL) {
	    phgError(1, "memory allocation error (%s:%d)\n", __FILE__,__LINE__);
	}
    }

    /* process 0 build similiar matrix */
    MPI_Gatherv(local, nonz, type, sml, len, ofst, type, 0, comm);

    MPI_Type_free(&type);
    phgFree(local);
    phgFree(len);
    phgFree(ofst);

    /* map algorithm by PLUM */
    if (rank == 0){
	int *part_map;
	int *proc_unmap;
	int *t1;
	int count;
	int kk;
	SM *p;

	part_map = (int *) phgAlloc(nprocs * sizeof(int));
	proc_unmap = (int *) phgAlloc(nprocs * sizeof(int));
	t1 = (int *) phgAlloc(nprocs * sizeof(int));

	if ((part_map == NULL) || (proc_unmap == NULL) || (t1 == NULL)) {
	    phgError(1, "memory allocation error (%s:%d)\n", __FILE__,__LINE__);
	}

	/* initialize */
	for (i = 0; i < nprocs; i++) {
	    proc_unmap[i] = 1;
	    part_map[i] = 0;
	    perm[i] = -1;
	}
	/* order the similiar matrix, descending order */
	qsort(sml, nt, sizeof(SM), comp_sml);

	/* map algorithm */
	count = 0;
	kk = 0;
	p = sml;
	nt -= 1;
	while ((count < nprocs) && (kk < nt)) {
	    while ( (proc_unmap[p->pid] == 0 || part_map[p->ptn] == 1)
		 && (kk < nt)) {
		p++;
		kk++;
	    }
	    if ((proc_unmap[p->pid] != 0) && (part_map[p->ptn] != 1)) {
		proc_unmap[p->pid] = 0;
		part_map[p->ptn] = 1;
		perm[p->ptn] = p->pid;
		count++;
	    }
	} /* map algorithm */

	count = 0;
	for (i = 0; i < nprocs; i++) {
	    if (perm[i] < 0)
		count++;
	    t1[i] = 0;
	}

	/* deal with special case that some arn't mapped */
	if (count > 0) {
	    for (i = 0; i < nprocs; i++) {
		if (perm[i] >= 0) {
		    t1[perm[i]] = -1;
		}
	    }
	    kk = 0;
	    for (i = 0; i < nprocs; i++) {
		if (perm[i] < 0) {
		    while (t1[kk] < 0) {
			kk++;
		    }
		    perm[i] = kk;
		    kk++;
		}
	    }
	}

	/* verification */
	/* count = 0: valid result   */
	/* count > 0: invalid result. unmap. */
	count = 0;
	for (i = 0; i < nprocs; i++) {
	    t1[i] = 0;

	    if ((perm[i] < 0) || (perm[i] > nprocs)) {
		count++;
	    }
	}

	/* check if each index is mapped just once 
	 * count = 0 if yes, > 0 if not */
	if (count == 0) {
	    for (i = 0; i < nprocs; i++) {
		t1[perm[i]] += 1;
	    }
	    for (i = 0; i < nprocs; i++) {
		if (t1[i] != 1)
		    count++;
	    }
	}

	if (count == 0) { /* succeed */
	    /* default: partitions unchanged. */
	    ret = 0;

	    /* check if partitions are changed */
	    for (i = 0; i < nprocs; i++) {
		if (perm[i] != i)
		    ret = 1;
	    }
	}
	else { /* fail */
	    ret = -1;
	}

	phgFree(part_map);
	phgFree(proc_unmap);
	phgFree(t1);
    }

    /* broadcast results */
    MPI_Bcast(&ret, 1, MPI_INT, 0, comm);
    if (ret >= 0) {
	MPI_Bcast(perm, nprocs, MPI_INT, 0, comm);
    }
    else if (ret < 0){
	for (i = 0; i < nprocs; i++) {
	    perm[i] = i;;
	}
	phgWarning("Remap failed (%s:%d). Orders "
		"keep unchanged.\n", __FILE__, __LINE__);
    }

    phgFree(sml);

    return ret;
}

static void
print_stats(MPI_Comm comm, int nprocs, const double *datasize, const int *perm)
/* prints total and max amount of data to be moved with and without remapping */
{
    int i, rank;
    double buffer[24], d, old_max, old_total, new_max, new_total;

    old_max = new_max = -1;
    old_total = new_total = 0;
    MPI_Comm_rank(comm, &rank);
    for (i = 0; i < nprocs; i++) {
	d = datasize[i];
	if (rank != i) {
	    old_total += d;
	    if (d > old_max)
		old_max = d;
	}
	if (rank != perm[i]) {
	    new_total += d;
	    if (d > new_max)
		new_max = d;
	}
    }

    buffer[0] = buffer[1] = buffer[2] = old_max;
    buffer[3] = buffer[4] = buffer[5] = old_total;
    buffer[6] = buffer[7] = buffer[8] = new_max;
    buffer[9] = buffer[10]= buffer[11]= new_total;
    MPI_Reduce(buffer, buffer + 12, 4, MPI_3DOUBLE, MPI_MSM, 0, comm);
    if (rank == 0) {
	printf("    Data move without remapping:\n");
	printf("\tmax(max):     %0.0f\n", buffer[12 + 2]);
	printf("\tmax(total):   %0.0f\n", buffer[15 + 2]);
	printf("\ttotal(total): %0.0f\n", buffer[15 + 1]);
	printf("    Data move with remapping:\n");
	printf("\tmax(max):     %0.0f\n", buffer[18 + 2]);
	printf("\tmax(total):   %0.0f\n", buffer[21 + 2]);
	printf("\ttotal(total): %0.0f\n", buffer[21 + 1]);
    }
}

/* Note: change 'FALSE' to 'TRUE' to enable the following macro */
#define check_permutation(n, perm) if (FALSE)			\
{								\
    int i, *tmp = phgAlloc(n * sizeof(*tmp));			\
    memcpy(tmp, perm, n * sizeof(*tmp));			\
    qsort(tmp, n, sizeof(*tmp), phgCompint);			\
    for (i = 0; i < n; i++)					\
	if (tmp[i] != i)					\
	    break;						\
    if (i < n) {						\
	phgInfo(-1, "%s:%d, permutation check failed.\n",	\
			__FILE__, __LINE__);			\
	phgAbort(1);						\
    }								\
    phgFree(tmp);						\
}

int
phgPartitionRemap0(MPI_Comm comm, const double datasize[], int perm[],
		   MPI_Comm newcomm)
/* Note: either ('comm' \subset 'newcomm') or ('newcomm' \subset 'comm'), and
 * this function must be called by all processes of 'comm' \cup 'newcomm' */
{
    int i, j, k, rank, nprocs, nprocs0, ret = 0;
    L2DATA *l2 = NULL;
    double *ds = NULL;
    int *perm0 = NULL;
    MPI_Comm comm_q;
    int size_q, rank_q;

    MPI_Comm_size(comm, &nprocs0);
    MPI_Comm_size(newcomm, &nprocs);
    MPI_Comm_rank(newcomm, &rank);

    if ((l2 = MPI_L2check(newcomm)) == NULL) {
	ret = partition_remap0(comm, datasize, perm, newcomm);
	goto end;
    }

    if (phgVerbosity > 0)
	phgPrintf("%s - remap partitions while preserving 2-level hierarchy.\n",
			__func__);

    if (l2->Lrank == 0) {
	if (!l2->Lsizes_are_equal) {
	    /* split Qcomm into sub-comms with equal Lsize */
	    MPI_Comm_split(l2->Qcomm, l2->Lsize, l2->Qrank, &comm_q);
	    MPI_Comm_size(comm_q, &size_q);
	    MPI_Comm_rank(comm_q, &rank_q);
	}
	else {
	    comm_q = l2->Qcomm;
	    size_q = l2->Qsize;
	    rank_q = l2->Qrank;
	}
    }

    /* Step 1: find permutation to minimize inter-node data move */

    if (l2->Lrank == 0)
	ds = phgCalloc(2 * l2->Qsize, sizeof(*ds));
    else
	ds = phgCalloc(l2->Qsize, sizeof(*ds));
    if (l2->Lrank == 0 || l2->Lsizes_are_equal)
	perm0 = phgAlloc((l2->Qsize + 1) * sizeof(*perm0));

    for (i = 0; rank < nprocs0 && i < nprocs; i++)
	ds[D2QMAP(l2, G2DMAP(l2, i))] += datasize[i];
    MPI_Reduce(ds, ds+l2->Qsize, l2->Qsize, MPI_DOUBLE, MPI_SUM, 0, l2->Lcomm);
    if (l2->Lrank == 0) {
	double *p = ds + l2->Qsize;
	MPI_Group g0, g1;
	if (comm_q != l2->Qcomm) {
	    MPI_Comm_group(l2->Qcomm, &g0);
	    MPI_Comm_group(comm_q,    &g1);
	    for (i = 0; i < size_q; i++) {
		MPI_Group_translate_ranks(g1, 1, &i, g0, &j);
		ds[i] = p[j];
	    }
	    p = ds;
	}
	ret = partition_remap0(comm_q, p, perm0, comm_q);
	if (comm_q != l2->Qcomm) {
	    /* Synchronize return values of different process groups.
	     * Note: we may get different return values from partition_remap0
	     * with different process groups, and we set ret > 0 on all groups
	     * as long as ret > 0 on one of them */
	    MPI_Allreduce(&ret, &i, 1, MPI_INT, MPI_MAX, l2->Qcomm);
	    if ((ret = i) > 0) {
		/* We need to translate perm0[] from comm_q to l2->Qcomm,
		 * and gather results from all comm_q's */
		int *perm1 = phgAlloc((l2->Qsize + 1) * sizeof(*perm1));
		for (i = 0; i < l2->Qsize; i++)
		    perm1[i] = -1;
		for (i = 0; i < size_q; i++) {
		    MPI_Group_translate_ranks(g1, 1, &i, g0, &j);
		    MPI_Group_translate_ranks(g1, 1, perm0 + i, g0, &k);
		    perm1[j] = k;
		}
		MPI_Allreduce(perm1, perm0, l2->Qsize, MPI_INT, MPI_MAX,
				l2->Qcomm);
		phgFree(perm1);
	    }
	    MPI_Group_free(&g1);
	    MPI_Group_free(&g0);
	    MPI_Comm_free(&comm_q);
	}
#if DEBUG
	if (ret > 0)
	    check_permutation(l2->Qsize, perm0);
#endif	/* DEBUG */
	perm0[l2->Qsize] = ret;
    }
    if (l2->Lsizes_are_equal) {
	MPI_Bcast(perm0, l2->Qsize + 1, MPI_INT, 0, l2->Lcomm);
	ret = perm0[l2->Qsize];
    }
    else {
	MPI_Bcast(&ret, 1, MPI_INT, 0, l2->Lcomm);
    }
    if (ret <= 0)
	goto end;

    /* Note: Q2Dmap is only available on local masters, and it is only needed
     * when l2->Lsizes_are_equal is not TRUE */
    if (l2->Lsizes_are_equal || l2->Lrank == 0) {
	int *D2Gmap = NULL;
	if (l2->G2Dmap != NULL) {
	    D2Gmap = phgAlloc(nprocs * sizeof(*D2Gmap));
	    for (i = 0; i < nprocs; i++)
		D2Gmap[G2DMAP(l2, i)] = i;
	}
#define D2GMAP(i) (D2Gmap == NULL ? (i) : D2Gmap[i])
	for (i = 0; i < nprocs; i++) {
	    j = D2QMAP(l2, G2DMAP(l2, i));
	    k = G2DMAP(l2, i) - Q2DMAP(l2, j);	/* the offset within a group */
	    perm[i] = D2GMAP(Q2DMAP(l2, perm0[j]) + k);
	}
#undef D2GMAP
	phgFree(D2Gmap);
    }
    if (!l2->Lsizes_are_equal)
	MPI_Bcast(perm, nprocs, MPI_INT, 0, l2->Lcomm);

    /* [TODO] Step 2: find permutation to minimize intra-node data move */

end:
    phgFree(ds);
    phgFree(perm0);

    /* print info */
    if (ret > 0 && phgVerbosity > 0 && rank < nprocs0)
	print_stats(comm, nprocs, datasize, perm);

    if (ret <= 0) {
	for (i = 0; i < nprocs; i++)
	    perm[i] = i;
    }
#if DEBUG
    else if (rank == 0) {
	check_permutation(nprocs, perm);
    }
#endif	/* DEBUG */

    return ret;
}

/* print grid quality information */
/* this function must be called just after phgBalanceGrid */
/* or e->mark isn't changed.   */
void
phgPartitionQuality(GRID *g)
{
    INT i, k;
    INT cr, ci, cb;
    INT tcr, tci, tcb;
    INT en;
    FLOAT si, tsi;
    ELEMENT *e;
    INT *pid;

    if (g->nroot <= 0)
	return;
    if (phgRank >= g->nprocs)
	return;

    en = ci = cb = cr = 0;
    ForAllElements(g, e) {
	en++;
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] & REMOTE) {
		cr++;
		continue;
	    }
	    if (e->bound_type[i] & INTERIOR) {
		ci++;
		continue;
	    }
	    cb++;
	}
    }

    /* load imbalance efficiency */
    MPI_Allreduce(&en, &k, 1, PHG_MPI_INT, MPI_MAX, g->comm);
    si = (FLOAT)k * g->nprocs / (FLOAT)g->nelem_global;
    phgPrintf("g->lif:	  %f\n", (double)si);

    /* total cut edges */
    MPI_Allreduce(&cr, &tcr, 1, PHG_MPI_INT, MPI_SUM, g->comm);
    tcr >>= 1;
    phgPrintf("Cut edges: %"dFMT"\n", tcr);

    /* surface indics */
    ci >>= 1;
    k = cb + ci + cr;
    /* local surface index */
    if (k > 0) {
	si = (FLOAT) cr / (FLOAT) k;
    }
    else {
	si = 0;
    }
    MPI_Allreduce(&si, &tsi, 1, PHG_MPI_FLOAT, PHG_MAX, g->comm);
    phgPrintf("Maximum local surface index: %f\n", (double)tsi);

    /* average surface index */
    MPI_Allreduce(&si, &tsi, 1, PHG_MPI_FLOAT, PHG_SUM, g->comm);
    si = tsi / (FLOAT)g->nprocs;
    phgPrintf("Average surface index:	    %f\n", (double)si);

    /* global surface index */
    MPI_Allreduce(&cb, &tcb, 1, PHG_MPI_INT, MPI_SUM, g->comm);
    MPI_Allreduce(&ci, &tci, 1, PHG_MPI_INT, MPI_SUM, g->comm);
    k = tci + tcr + tcb;

    si = (FLOAT) tcr / (FLOAT)k;
    phgPrintf("Global surface index:	    %f\n", (double)si);

    /* compute interprocess connectivity */
    pid = phgAlloc((g->nprocs + 1) * sizeof(*pid));
    for (i = 0; i < g->nprocs; i++) {
	if (i != phgRank && g->neighbours.counts[i] != 0) {
	    pid[i] = 1;
	}
	else {
	    pid[i] = 0;
	}
    }
    k = 0;
    for (i = 0; i < g->nprocs; i++)
	k += pid[i];
    phgFree(pid);
    MPI_Allreduce(&k, &i, 1, PHG_MPI_INT, MPI_MAX, g->comm);
    phgPrintf("Max interprocess connectivity: %"dFMT"\n", i);
}

BOOLEAN
phgPartitionerRandom(GRID *g, MPI_Comm newcomm, DOF *weights, FLOAT power)
/* Partitioner random. */
{
    int nprocs, myrank;
    ELEMENT *e;
    static BOOLEAN initialized = FALSE;
    static INT seed = 0;

    if (!initialized) {
	/* register options (nothing to do, simply return) */
	phgOptionsRegisterTitle("\n[Random partitioner]", "\n", "partition");
	phgOptionsRegisterInt("-partition_random_seed",
		"Seed for the random partitioner (0 => use PID)", &seed);
	initialized = TRUE;
	return FALSE;
    }

    MPI_Comm_size(newcomm, &nprocs);
    MPI_Comm_rank(newcomm, &myrank);
    if (seed != 0)
	srandom(myrank + (unsigned int)seed);
    else
	srandom((unsigned int)getpid());
    ForAllElements(g, e) {
	long int i = random();
	e->mark = (i > 0 ? i : -i) % nprocs;
    }
    
    return TRUE;
}

BOOLEAN
phgPartitionerUser(GRID *g, MPI_Comm newcomm, DOF *weights, FLOAT power)
/* Partitioner user. Assume e->mark is already set by user. */
{
    int nprocs, myrank, rank;
    INT *nelems, *nelems_global;
    ELEMENT *e;
    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	/* register options (nothing to do, simply return) */
	phgInfo(0, "register.\n");
	initialized = TRUE;
	return FALSE;
    }

    MPI_Comm_size(newcomm, &nprocs);
    MPI_Comm_rank(newcomm, &myrank);

    phgInfo(3, "nprocs: %d\n", nprocs);

    nelems = phgAlloc(2 * nprocs * sizeof(*nelems));
    nelems_global = nelems + nprocs;
    bzero(nelems, nprocs * sizeof(*nelems));
    ForAllElements(g, e) {
	assert(e->mark >= 0 && e->mark < nprocs);
	nelems[e->mark]++;
    }
    for (rank = 0; rank < nprocs; rank++)
	phgInfo(3, "  Local  elems on rank %3d: %3d\n", rank, nelems[rank]);
    MPI_Allreduce(nelems, nelems_global, nprocs,
		  PHG_MPI_INT, MPI_SUM, g->comm);
    for (rank = 0; rank < nprocs; rank++)
	phgInfo(3, "  Global elems on rank %3d: %3d\n", rank, nelems[rank]);
    
    phgFree(nelems);

    /* Note: user should make sure that the partition succeed
     * and the grid is changed. */
    return TRUE;
}


#endif	/* USE_MPI */
