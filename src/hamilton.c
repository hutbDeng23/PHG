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

/* $Id: hamilton.c,v 1.34 2015/01/07 09:14:19 zlb Exp $ */

/* we have gave two implementations of Hamilton path.
 * They are separated by "#if HEAD_INSERTION...#else....#endif"
 * Default one has better quality.
 * 
 * If users want to use another one, 
 * users have to replace HEAD_INSERTION0 by 1.
*/

#include "phg.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <strings.h>		/* bzero() */

#define HEAD_INSERTION 0

/* ******************************************
 * new struct for phgHamiltonPath 
 * save information for a node
 * *****************************************/
typedef struct hmlt_ {
    struct hmlt_ *next;
    struct hmlt_ *pre;
    ELEMENT *e;
    INT inv;		 /* in vertex */
    INT outv;		 /* out vertex */
    CHAR flag;		  /* -1: head, 0: interior, 1: tail */
} hmlt;

typedef struct ninfo_ {
    struct ninfo_ *next;
    INT index;
}ninfo;       /* neighbour information node */

static hmlt *ham;	     /* spaces for hmlt */
static INT *in_list;	     /* useful neighbours infor */
static ninfo *pinfo;	     /* save node that has useful neighbour */


/* **********************************************
 * locate spaces, read informatin 
 * ham[i] stors information ahout element i 
 * vertices, neighbours, common faces,
 * index
 * *********************************************/
static BOOLEAN
phgHamiltonInitS(GRID *g)
{
    ELEMENT *e;
    int i, n;
    INT kk;
    hmlt *p;

    /* these spaces are used only by process 0 */
    ham = (hmlt *)phgAlloc(g->nroot * sizeof(hmlt));
    in_list = (INT *)phgAlloc(g->nroot * sizeof(INT));
    pinfo = (ninfo *)phgAlloc(g->nroot * sizeof(ninfo));
    
    p = ham;
    kk = g->nroot - 2;
    for (i = 0; i < g->nroot; i++){
	in_list[i] = 1;
	pinfo[i].next = NULL;

	e = g->roots + i;
	n = e->index;
	p[n].e = e;
	if (n <= kk)
	    p[n].next = &p[n + 1];
	else 
	    p[n].next = NULL;
	if (n > 0)
	    p[n].pre = p + n - 1;
	else
	    p[n].pre = NULL;
	
    } /* scan all roots */

    return TRUE;
}

static BOOLEAN
HasBadEdges(GRID *g)
/* test if grid g has bad edges */
{
    INT i, bad_edges = 0;
    ELEMENT *e;
    int j, k;
    double *length;	   /* (square of) length of all edges */
    BOOLEAN ret;

    if (g == NULL)
	phgError(1, "g is NULL!\n");

    phgInfo(2, "Test if grid g has bad edges\n");

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
	    if (length[k] >= 0) continue;
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
    if (bad_edges > 0) {
	ret = TRUE;
    }
    else {
	ret = FALSE;
    }
    
    phgFree(length);
    return ret;
}
/* ****************************************
 * find any neighour that is avaliable 
 * return script if there is. 0, 1, 2, 3.
 * or return -1		  
 * ***************************************/
#define find_neighbour(p)								      \
((((p)->e->neighbours[0] != NULL) && in_list[((ELEMENT *)(p)->e->neighbours[0])->index])? 0 : \
((((p)->e->neighbours[1] != NULL) && in_list[((ELEMENT *)(p)->e->neighbours[1])->index])? 1 : \
((((p)->e->neighbours[2] != NULL) && in_list[((ELEMENT *)(p)->e->neighbours[2])->index])? 2 : \
((((p)->e->neighbours[3] != NULL) && in_list[((ELEMENT *)(p)->e->neighbours[3])->index])? 3 : \
		     -1))))

/* ***************************************
 * construct Hamiltonian Path 
 * serial version 
 * **************************************/
int
phgHamiltonPath(GRID *g)
{
    hmlt *pnew = NULL;	   /* result list */
    hmlt *pold = NULL;	   /* initial list */
    hmlt *pc, *pre, *pt;   /* tmp */
    hmlt *nhead;   /* head */
    hmlt *ntail;   /* tail */
    ELEMENT *ee;
    int ps;
    ninfo *nep;    /* neighbour information */
    ninfo *ptail;  /* neighbour information */
    INT n = 0;
    INT ret = FALSE;

    /* TODO: deal with bad edges */
    if (HasBadEdges(g)) {
	phgWarning("GRID g hase bad edges!\n");
	phgWarning("Hamiltonian Path was ** NOT ** created!\n");
	phgWarning("Hamiltonian Path is used by RTK repartitioner.\n");
	phgWarning("The RTK repartitioner still works, but the result maybe bad\n");
	phgWarning("You'd better use other repartitioners, eg. metis\n");

	return 1;
    }

    /* initial */
    phgHamiltonInitS(g);

    /* just want to make compiler happy */
    ptail = NULL;

    /* construct the Hamiltonian path */
    /* G1 is the result. pnew */
    if (g->nroot > 2){
	/* initial set G1, has two elements */
	/* construct G1*/
	pnew = ham;
	in_list[pnew->e->index] = 0;
	/* pinfo */
#if HEAD_INSERTION    /* head insert */
	if (find_neighbour(pnew) >= 0){
	    nep = &pinfo[pnew->e->index];
	    nep->index = pnew->e->index;
	}
	else {
	    nep = NULL;
	}/* pinfo */
#else
	if (find_neighbour(pnew) >= 0){
	    nep = &pinfo[pnew->e->index];
	    nep->index = pnew->e->index;
	    ptail = nep;
	}
	else {
	    nep = NULL;
	    ptail = NULL;
	}/* pinfo */
#endif
	
	pnew->pre = NULL;
	pold = pnew->next;
	pold ->pre = NULL;
	pnew->next = NULL;

	/* add the second element to pnew */
	ps = find_neighbour(pnew);
	ee = (ELEMENT *) pnew->e->neighbours[ps];
	pc = &ham[ee->index];
	in_list[pc->e->index] = 0;

	/* pinfo  */
#if HEAD_INSERTION    
	if (find_neighbour(pc) >= 0){
	    if (nep != NULL){
		pinfo[pc->e->index].next = nep;
		nep = &pinfo[pc->e->index];
		nep->index = pc->e->index;
	    }
	    else {
		nep = &pinfo[pc->e->index];
		nep->index = pc->e->index;
	    }
	}
#else
	if (find_neighbour(pc) >= 0){
	    if (ptail != NULL){
		ptail->next = &pinfo[pc->e->index];
		ptail->next->index = pc->e->index;
		ptail = &pinfo[pc->e->index];
	    }
	    else {
		nep = &pinfo[pc->e->index];
		nep->index = pc->e->index;
		ptail = nep;
	    }
	}
#endif	      

	pre = pc->pre;

	/* del an elemnt from pold */
	if (pre != NULL){
	    if (pc->next != NULL) {
		pre->next = pc->next; /* normal node */
		pc->next->pre = pre;
	    }
	    else {		      /* tailed node */
		pre->next = NULL;
	    }
	}
	else {
	    pold = pold->next;	  /* head node */
	    pold->pre = NULL;
	}
	/* add an element to pnew, two elements */
	pnew->next = pc;   
	pc->next = NULL;
	pc->pre = pnew;
	
	/* update information for G1 */
	pnew->flag = -1;
	pnew->next->flag = 1;
	nhead = pnew;  /* head of pnew */
	ntail = pc;    /* head of pnew */

	/* common vertex */
	if (ps == 0){
	    nhead->outv = pnew->e->verts[1];
	    ntail->inv = pnew->e->verts[1];
	}
	else if (ps == 1){
	    nhead->outv = pnew->e->verts[2];
	    ntail->inv = pnew->e->verts[2];
	}
	else if (ps == 2){
	    nhead->outv = pnew->e->verts[3];
	    ntail->inv = pnew->e->verts[3];
	}
	else if (ps == 3){
	    nhead->outv = pnew->e->verts[0];
	    ntail->inv = pnew->e->verts[0];
	}
	
	/* be careful when update information about pnew and pold */
	/* flag: -1, 0, 1 */
	/* vertices, -1, 0, 1 */
	/* construct the whole Hamiltonian path */
	n = g->nroot - 2;
	while (pold != NULL && n > 0){
	    /* find a avaliable neighhour in pold 
	     * pc has a neighbour in pold
	     * and pc is in pnew		   */	     
	    if (nep != NULL){
		pc = &ham[nep->index];
		ps = find_neighbour(pc);
		while (ps < 0){
		     nep = nep->next;
		     if (nep == NULL)
			 goto stat;
		     pc = &ham[nep->index];
		     ps = find_neighbour(pc);
		}
	    }
	    else {
		goto stat;
	    }
	    
	    ee = (ELEMENT *) pc->e->neighbours[ps];
	    pt = &ham[ee->index]; /* neighbour in pold */

	    /* del the element from pold */
	    pre = pt->pre;
	    if (pre != NULL){ /* normal node */
		if (pt->next != NULL){
		    pt->next->pre = pre;
		    pre->next = pt->next;
		}
		else {	  /* tailed node */
		   pre->next = NULL;
		}
	    }
	    else if (n > 1){  /* head */
		pold = pold->next;
		pold->pre = NULL;
	    }
	    else {	      /* pold has no element any more */
		pold = NULL;
	    }
	    pt->next = NULL;
	    pt->pre = NULL;
	    in_list[pt->e->index] = 0;
	/* pinfo  */
#if HEAD_INSERTION    
	if (find_neighbour(pt) >= 0){
	    if (nep != NULL){
		pinfo[pt->e->index].next = nep;
		nep = &pinfo[pt->e->index];
		nep->index = pt->e->index;
	    }
	    else {
		nep = &pinfo[pt->e->index];
		nep->index = pt->e->index;
	    }
	}
#else
	if (find_neighbour(pt) >= 0){
	    if (ptail != NULL){
		ptail->next = pinfo + pt->e->index;
		ptail->next->index = pt->e->index;
		ptail = pinfo + pt->e->index;
		ptail->next = NULL;
	    }
	    else {
		nep = &pinfo[pt->e->index];
		nep->index = pt->e->index;
		ptail = nep;
	    }
	}
#endif

	    /* add element to pnew and update information */
	    /* pt is the element will be added to pnew , pt is in pold*/
	    /* pc is the element which is the neighbour of pt, pc in pnew */
	    /* ps is the number of common face of pc */
	    /* there are three conditions: pc is head, tail, or interior node */

	    assert(ps >= 0);
	    {
		INT co[3];	/* vertices in common face */
		INT c1, ct = 0;
		int i;

		/* update common vertices information */
		bzero(co, sizeof(co));	/* to avoid a gcc-4.4.1 warning */
		if (ps == 0){
		    co[0] = pc->e->verts[1];
		    co[1] = pc->e->verts[2];
		    co[2] = pc->e->verts[3];
		}
		else if (ps == 1){
		    co[0] = pc->e->verts[0];
		    co[1] = pc->e->verts[2];
		    co[2] = pc->e->verts[3];
		}
		else if (ps == 2){
		    co[0] = pc->e->verts[0];
		    co[1] = pc->e->verts[1];
		    co[2] = pc->e->verts[3];
		}
		else if (ps == 3){
		    co[0] = pc->e->verts[0];
		    co[1] = pc->e->verts[1];
		    co[2] = pc->e->verts[2];
		}
		
		/* normal interior node */
		if ((pc != nhead) && (pc != ntail)){
		    /* check information about in-vertex or out-vertex */
		    for (i = 0; i < 3; i++){
			if ((co[i] == pc->inv) || (co[i] == pc->outv))
			    break;
		    }
		    if (i >= 3)
			phgError(1, "unexpected.\n");
		    else
			ct = co[i];
		    if (ct == pc->outv){  /* out-vertex */
			for (i = 0; i < 3; i++){
			    if ((co[i] != pc->inv) && (co[i] != pc->outv))
				break;
			}
			c1 = co[i];
			
			pre = pc->next;
			pt->next = pre;
			pre->pre = pt;
			pt->pre = pc;
			pc->next = pt;
			
			pc->outv = c1;
			pt->inv = c1;
			pt->outv = ct;
			pt->flag = 0;
		    }
		    else if (ct == pc->inv){ /* in-vertex */
			for (i = 0; i < 3; i++){
			    if ((co[i] != pc->inv) && (co[i] != pc->outv))
				break;
			}
			c1 = co[i]; /* in of pc, out of pt */
			pre = pc->pre;

			pre->next = pt;
			pt->pre = pre;
			pt->next = pc;
			pc->pre = pt;

			pre->outv = ct;
			pt->inv = ct;
			pt->outv = c1;
			pt->flag = 0;
			pc->inv = c1;

		    }

		}
		else if (pc == nhead){	/* head */
		    for (i = 0; i < 3; i++){
			if (co[i] != pc->outv)
			    break;
		    }
		    c1 = co[i];
		    pt->next = pc;
		    pt->pre = NULL;
		    pc->pre = pt;
		    
		    pc->flag = 0;
		    pc->inv = c1;
		    pt->flag = -1;
		    pt->outv = c1;
		    nhead = pt;
		    pnew = nhead;
		}
		else if (pc == ntail){ /* tail */
		    for (i = 0; i < 3; i++){
			if (co[i] != pc->inv)
			    break;
		    }
		    ct = co[i];

		    pc->next = pt;
		    pt->next = NULL;
		    pt->pre = pc;

		    pc->flag = 0;
		    pc->outv = ct;
		    pt->flag = 1;
		    pt->inv = ct;
		    ntail = pt;
		}
		else {
		    phgError(1, "Hamiltonian path constructing error \n");
		}
	    }
	    n--;
	}
    }
    else if (g->nroot == 2){
	pnew = ham;
	pnew->next = &ham[1];
	ham[1].next = NULL;

    }
    else if (g->nroot == 1){
	pnew = ham;
	pnew->next = NULL;

    }
    else {
	phgError(1, "Error: g->nroot <= 0\n");
    }

stat:
    /* remap the global number of elements */
    if (n == 0) {  /* success */
	INT *augtree, kk, jj, index;
	ELEMENT *e;
	INT i;
	hmlt *ss = pnew;
	INT sc = 0;

	/* verification result */
	/* initial */
	for (i = 0; i < g->nroot; i++){
	    in_list[i] = 1;
	}
	while (ss != NULL && sc <= g->nroot){
	    sc++;
	    in_list[ss->e->index] = 0;
	    ss = ss->next;
	}
	n = g->nroot;
	for (i = 0; i < g->nroot; i++){
	    if (in_list[i] == 0)  n--;
	}

	if (n == 0) {  /* verification success */
	    augtree = (INT *)phgAlloc(g->nroot * sizeof(INT));
	    kk = g->nroot;

	    for (jj = 0; jj < kk; jj++) {
		augtree[pnew->e->index] = jj;
		pnew = pnew->next;
	    }
	    for (jj = 0; jj < kk; jj++) {
		e = g->roots + jj;
		index = e->index;
		e->index = augtree[index];
	    }
	    phgFree(augtree);
	    ret = 0;
	}
	else {	   /* verification fail */
	    ret = 1;
	}
    }
    else {   /* fail */
	ret = 1;
    }
    phgFree(ham);
    phgFree(pinfo);
    phgFree(in_list);

    return ret;
}
