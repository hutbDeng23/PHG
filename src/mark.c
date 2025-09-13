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

/* This file implements marking strategies:
 *   MAX, MNS, GERS and EQDIST
 *
 * Author: Liu Hui (liuhui@lsec.cc.ac.cn)
 *
 * $Id: mark.c,v 1.44 2015/06/19 08:58:42 zlb Exp $ */

#include "phg.h"

#include <math.h>
#include <string.h>

/* # of recursions */
static int phgMarkRecur;  

/****************** determine factor ********/
/* sum = \sum {data < factor * mx */
/* or sum = \sum {data < factor * mx 
 *
 * x: error indicator
 * mx: sentinel (data that smaller than factor * mx will be summed
 * factor: usually 0 < factor < 1., used with mx
 * flag: determine data are stored in element or face 
 * con:  determine element that satisfies the condition to be summed
 *       or not (1 not, 0, yes). 
 * */
static FLOAT
part_sum(DOF *x, FLOAT mx, FLOAT factor, int flag, int con)
{
    GRID *g = NULL;
    ELEMENT *e;
    INT i;
    FLOAT sum, b;
    FLOAT tmp, a;
    FLOAT ret = 0.;

    assert(x != NULL);
    assert(factor >= 0. && mx >= 0. && factor <= 1.);
    assert((flag == 1 && con == 0) || /* refine: face data: mns */
           (flag == 0 && con == 0) || /* refine: element data */
           (flag == 0 && con == 1) || /* refine: element data */
           (flag == -1));             /* coarsen: element data */

    tmp = mx * factor;
    g = x->g;

    /*********** refine ************/
    if (flag >= 0) {
        if (flag == 1) {
            /* flag = 1,  face data, e->mark = 0 */
            /* just interiro faces */
            if (con == 0) {
                sum = 0.0;
                ForAllElements(g, e){
                    for (i = 0; i < NFace; i++){
                        if (e->bound_type[i] & INTERIOR) {
                            b = *DofFaceData(x, e->faces[i]);
                            if (b > tmp){
                                sum += b;
                            }
                        }
                    }
                }

#if USE_MPI
                MPI_Allreduce(&sum, &b, 1, PHG_MPI_FLOAT, PHG_SUM, x->g->comm);
                ret = b * 0.5;
#else
                ret = sum * 0.5;
#endif
            }
            else {
                phgError(1, "wrong parameter: con %s %d\n",__FILE__,__LINE__);
            }
        } /* flag = 1 */
        /* element data */
        else if (flag == 0) {
            if ( con == 0){    /* element data */
                sum = 0.0;
                ForAllElements(g, e){
                    b = *DofElementData(x, e->index);
                    if (b > tmp){
                        sum += b;
                    }
                }
#if USE_MPI
                MPI_Allreduce(&sum, &b, 1, PHG_MPI_FLOAT, PHG_SUM, x->g->comm);
                ret = b;
#else
                ret = sum;
#endif
            }
            else if (con == 1){
                /* con = 1 and flag = 0, data are stored in element */
                /* some marks are marked */
                b = 0.0;
                ForAllElements(g, e){
                    if (e->mark >= 1){
                        b += *DofElementData(x, e->index);
                    }
                }
#if USE_MPI
                MPI_Allreduce(&b, &sum, 1, PHG_MPI_FLOAT, PHG_SUM, x->g->comm);
#else
                sum = b;
#endif

                a = 0.0;
                ForAllElements(g, e){
                    if (e->mark == 0){
                        b = *DofElementData(x, e->index);
                        if (b > tmp){
                            a += b;
                        }
                    }
                }
#if USE_MPI
                MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_SUM, x->g->comm);
                sum += b;
                ret = sum;
#else
                sum += a;
                ret = sum;
#endif
            }
            else {
                phgError(1, "wrong parameter: con %s %d\n",__FILE__,__LINE__);
            }
        } /* flag == 0 */
    }   /* flag >= 0 */
    else if (flag < 0) {
        /*********** coarsen ***********/
        if (flag == -1) { /* element data only */
            sum = 0.0;
            ForAllElements(g, e){
                b = *DofElementData(x, e->index);
                if (b < tmp){
                    sum += b;
                }
            }
#if USE_MPI
            MPI_Allreduce(&sum, &b, 1, PHG_MPI_FLOAT, PHG_SUM, x->g->comm);
            ret = b;
#else
            ret = sum;
#endif
        }
        else {
            phgError(1, "wrong parameter: flag %s %d\n", __FILE__, __LINE__);
        }
    }
    return ret;
}


/* this function returns a factor  such that 
 * \sum{e \in submesh M'} {error_e} > (or <=) eta
 * 
 * error: error indicator
 * mx: sentinel, usually equals to phgDofNormInftyVec(error)
 * eta: sentinel, usually equals to phgDofNormL1Vec(error) * \alpha
 *      0. < alpha < 1.
 * return value belongs to [low, high], usually they belong to
 * [0., 1.].   In general, low = 0.; high = 1.
 * flag: data type, element or face element 
 * con:  determine element that satisfies the condition to be summed
 *       or not (1 not, 0, yes). used by part_sum.
*/
static FLOAT
mark_cal_factor(DOF *error, FLOAT mx, FLOAT eta, 
        FLOAT low, FLOAT up, int flag, int con)
{
    FLOAT sm;
    FLOAT c1, c2, c3;
    static int tmp = 0;
    static FLOAT low_sum, up_sum;
    BOOLEAN ret = FALSE;

    phgMarkRecur++;
    assert(low >= 0.0 && up >= low && up <= 1.0);
    assert(mx >= 0. && eta >= 0.);
    assert(error != NULL);

    /***************  refine *************/
    if (flag >= 0){
        tmp++;
        sm = (low + up) * 0.5;
        if (tmp > 1){
            c2 = part_sum(error, mx, sm, flag, con);
        }
        else {
            c1 = part_sum(error, mx, low, flag, con);
            c2 = part_sum(error, mx, sm, flag, con);
            c3 = part_sum(error, mx, up, flag, con);
            low_sum = c1;
            up_sum = c3;
        }
#if 1
        if (phgRank == 0) {
            if (low_sum == c2 || c2 == up_sum || 
                    c2 == eta || up - low < 0.0001){
                ret = TRUE;
            }
            else {
                ret = FALSE;
            }
        }
#if USE_MPI
        MPI_Bcast(&ret, 1, PHG_MPI_INT, 0, error->g->comm);
#endif

        if (ret) {
            tmp = 0;
            return sm;
        }
#else
        if (phgRank == 0) {
            if (up - low < 0.0001) {
                ret = TRUE;
            }
            else {
                ret = FALSE;
            }
        }
#if USE_MPI
        MPI_Bcast(&ret, 1, PHG_MPI_INT, 0, error->g->comm);
#endif
        if (ret) {
            tmp = 0;
            return sm;
        }
#endif

        if (c2 > eta ){
            low_sum = c2;
            return mark_cal_factor(error, mx, eta, sm, up, flag, con);
        }
        else {
            up_sum = c2;
            return mark_cal_factor(error, mx, eta, low, sm, flag, con);
        }
    }

    /**************  coarsen **************/
    sm = (low + up) * 0.5;

    if (tmp > 1){
        c2 = part_sum(error, mx, sm, flag, con);
    }
    else {
        c1 = part_sum(error, mx, low, flag, con);
        c2 = part_sum(error, mx, sm, flag, con);
        c3 = part_sum(error, mx, up, flag, con);
        low_sum = c1;
        up_sum = c2;
    }

#if 0
        if (phgRank == 0) {
            if ((low_sum == c2 || c2 == up_sum || c2 == eta) 
                     &&  up - low < 0.0001){
                ret = TRUE;
            }
            else {
                ret = FALSE;
            }
        }
#if USE_MPI
        MPI_Bcast(&ret, 1, PHG_MPI_INT, 0, error->g->comm);
#endif
        if (ret) {
            tmp = 0;
            return sm;
        }
#else
        if (phgRank == 0) {
            if (up - low < 0.0001) {
                ret = TRUE;
            }
            else {
                ret = FALSE;
            }
        }
#if USE_MPI
        MPI_Bcast(&ret, 1, PHG_MPI_INT, 0, error->g->comm);
#endif
        if (ret) {
            tmp = 0;
            return sm;
        }
#endif

    if (c2 < eta ){
        low_sum = c2;
        return mark_cal_factor(error, mx, eta, sm, up, flag, con);
    }
    else {
        up_sum = c2;
        return mark_cal_factor(error, mx, eta, low, sm, flag, con);
    }
}

/***************************************************************/
/*                     refine subroutine                       */
/*                                                             */
/***************************************************************/

/* init: e->mark = 0 */
static void
mark_init(GRID *g)
{
    ELEMENT *e;

    ForAllElements(g, e) {
        e->mark = 0;
    }
}

/* GERS: mark elements only uses data in faces */
static void
mark_refine_gers_face(DOF *error, FLOAT fac, int  depth)
{
    GRID *g = error->g;
    ELEMENT *e;
    FLOAT tmpc;
    FLOAT tmp = 0.;
    FLOAT gama;
    FLOAT sum, flag, PE_error;
    int i;
    long count, cc;

    assert(error != NULL);
    assert(fac >= 0. && fac <= 1.);
    assert(error->type->np_face == 1);
    assert(depth > 0);

    if (phgVerbosity > 0) {
        phgPrintf("\n*** GERS, Face Data, Refine\n");
    }

    /* determine factor */
    phgMarkRecur = -1;
    PE_error = phgDofNormInftyVec(error);
    if (PE_error <= 0.0) {
        return;
    }
    flag = phgDofNormL1Vec(error) * fac;

    gama = mark_cal_factor(error, PE_error, flag, 0.0, 1.0, 1, 0);
    if (phgVerbosity > 0){
        sum = part_sum(error, PE_error, gama, 1, 0);
        phgPrintf("recursion times: %d\n", phgMarkRecur);
        phgPrintf("sentinel: %.8lf  marked: %.8lf  theta: %.8lf\n",
                (double)flag, (double)sum, (double)gama);
    }

    /* mark */
    tmpc = PE_error * gama;
    cc = 0;
    ForAllElements(g, e){
        for (i = 0; i < NFace; i++){
            if (e->bound_type[i] & BDRY_MASK)
                continue;
            tmp = *DofFaceData(error, e->faces[i]);
            if (tmp > tmpc){
                e->mark = depth;
                cc++;
            }
        }
    }

    /* number of elements marked */
    if (phgVerbosity > 0) {
#if USE_MPI
        MPI_Allreduce(&cc, &count, 1, MPI_LONG, MPI_SUM, g->comm);
#else
        count = cc;
#endif
        phgPrintf("Elements marked for refining: %ld\n", count);
    }
    return;
}

/* GERS: mark elements only uses element data */
static void
mark_refine_gers_elem(DOF *error, FLOAT fac, int depth)
{
    GRID *g = error->g;
    ELEMENT *e;
    FLOAT gama;
    FLOAT sum, flag, pe;
    long cc, count;
    
    assert(error != NULL);
    assert(error->type->np_elem == 1);
    assert(fac >= 0. && fac <= 1.);
    assert(depth > 0);
    if (phgVerbosity > 0) {
        phgPrintf("\n*** GERS, Element Data, Refine\n");
    }

    /* determine factor */
    phgMarkRecur = -1;
    pe = phgDofNormInftyVec(error);
    if (pe <= 0.0)
        return;

    flag = phgDofNormL1Vec(error) * fac;
    gama = mark_cal_factor(error, pe, flag, 0.0, 1.0, 0, 0);
    if (phgVerbosity > 0){
        sum = part_sum(error, pe, gama, 0, 0);
        phgPrintf("recursion times: %d\n", phgMarkRecur);
        phgPrintf("sentinel: %.8lf  marked: %.8lf  theta: %.8lf\n",
                        (double)flag, (double)sum, (double)gama);
    }
    
    /* mark */
    pe *= gama;
    cc = 0;
    ForAllElements(g, e){
        gama = *DofElementData(error, e->index);
        if (gama > pe){
            e->mark = depth;
            cc++;
        }
    }

    /* number of elements marked */
    if (phgVerbosity > 0) {
#if USE_MPI
        MPI_Allreduce(&cc, &count, 1, MPI_LONG, MPI_SUM, g->comm);
#else
        count = cc;
#endif
        phgPrintf("Elements marked for refining: %ld\n", count);
    }
    return;
}

/* MNS: GERS on face data (error indicator) 
 *      and GERS on element data (data oscillation)
 * GERS: mark elements uses error in faces 
 * and data oscillation in elements */
static void
mark_refine_mns(DOF *error, FLOAT fac, DOF *osc, FLOAT fac2, int depth)
{
    GRID *g = error->g;
    ELEMENT *e;
    FLOAT gama;
    FLOAT sum, flag, PE_error;
    int i;
    long count, cc;
    
    assert(error != NULL);
    assert(error->type->np_face == 1);
    assert(osc != NULL);
    assert(osc->type->np_elem == 1);
    assert(fac >= 0. && fac <= 1.);
    assert(fac2 >= 0. && fac2 <= 1.);
    assert(depth > 0);
    if (phgVerbosity > 0)
        phgPrintf("\n*** MNS, Face Data, Refine\n");
   
    /* determine factor */
    phgMarkRecur = -1;
    PE_error = phgDofNormInftyVec(error);
    if (PE_error <= 0.)
        return;

    flag = phgDofNormL1Vec(error) * fac;
    gama = mark_cal_factor(error, PE_error, flag, 0.0, 1.0, 1, 0);
    if (phgVerbosity > 0){
        sum = part_sum(error, PE_error, gama, 1, 0);
        phgPrintf("recursion times: %d\n", phgMarkRecur);
        phgPrintf("sentinel: %.8lf  marked: %.8lf  theta: %.8lf\n",
                        (double)flag, (double)sum, (double)gama);
    }
    
    /* mark by error indicator */
    sum = PE_error * gama;
    cc = 0;
    ForAllElements(g, e){
        for (i = 0; i < NFace; i++){
            if (e->bound_type[i] & BDRY_MASK)
                continue;
            gama = *DofFaceData(error, e->faces[i]);
            if (gama > sum){
                e->mark = depth;
                cc++;
            }
        }
    }

    /* mark by data oscillation */
    phgMarkRecur = -1;
    PE_error = phgDofNormInftyVec(osc);
    if (PE_error <= 0.0) 
        goto stat;
    flag = phgDofNormL1Vec(osc) * fac2;

    /* determine factor */
    sum = part_sum(osc, PE_error, 1.0, 0, 1);
    if (sum >= flag)
        goto stat;
    gama = mark_cal_factor(osc, PE_error, flag, 0.0, 1.0, 0, 1);
    if (phgVerbosity > 0){
        sum = part_sum(osc, PE_error, gama, 0, 1);
        phgPrintf("\n***  Marked by Data oscillation\n");
        phgPrintf("recursion times: %d\n", phgMarkRecur);
        phgPrintf("sentinel: %.8lf  marked: %.8lf  theta: %.8lf\n",
                        (double)flag, (double)sum, (double)gama);
    }
    
    sum = PE_error * gama;
    ForAllElements(g, e){
        gama = *DofElementData(osc, e->index);
        if (e->mark == 0 && gama > sum){
            e->mark = depth;
            cc++;
        }
    }
    /* number of elements marked */
stat:
    if (phgVerbosity > 0) {
#if USE_MPI
        MPI_Allreduce(&cc, &count, 1, MPI_LONG, MPI_SUM, g->comm);
#else
        count = cc;
#endif
        phgPrintf("Elements marked for refining: %ld\n", count);
     }
    return;
}

/* GERS: mark elements uses error in elements */
/* and data oscillation in elements */
static void
mark_refine_gers_osc(DOF *error, FLOAT fac, DOF *osc, FLOAT fac2, int depth)
{
    GRID *g = osc->g;
    ELEMENT *e;
    FLOAT gama, pe;
    FLOAT sum, flag;
    long cc, count;
    
    assert(error != NULL);
    assert(error->type->np_elem == 1);
    assert(osc != NULL);
    assert(osc->type->np_elem == 1);
    assert(fac >= 0. && fac <= 1.);
    assert(fac2 >= 0. && fac <= 1.);
    assert(depth > 0);
    if (phgVerbosity > 0)
        phgPrintf("\n*** GERS, element data && data oscillation, Refine\n");
   
    /* determine factor */
    phgMarkRecur = -1;
    pe = phgDofNormInftyVec(error);
    if (pe <= 0.)
        return;

    flag = phgDofNormL1Vec(error) * fac;
    gama = mark_cal_factor(error, pe, flag, 0.0, 1.0, 0, 0);
    if (phgVerbosity > 0){
        sum = part_sum(error, pe, gama, 0, 0);
        phgPrintf("recursion times: %d\n", phgMarkRecur);
        phgPrintf("sentinel: %.8lf  marked: %.8lf  theta: %.8lf\n",
                        (double)flag, (double)sum, (double)gama);
    }
    
    /* mark by error indicator */
    pe *= gama;
    cc = 0;
    ForAllElements(g, e){
        gama = *DofElementData(error, e->index);
        if (gama > pe){
            e->mark = depth;
            cc++;
        }
    }

    /* mark by data oscillation */
    /* determine factor */
    phgMarkRecur = -1;
    pe = phgDofNormInftyVec(osc);
    if (pe == 0.)
        goto stat;
    flag = phgDofNormL1Vec(osc) * fac2;
    sum = part_sum(osc, pe, 1.0, 0, 1);
    if (sum >= flag)
        goto stat;
    gama = mark_cal_factor(osc, pe, flag, 0.0, 1.0, 0, 1);
    if (phgVerbosity > 0){
        sum = part_sum(osc, pe, gama, 0, 1);
        phgPrintf("\n*** Marked by Data oscillation\n");
        phgPrintf("recursion times: %d\n", phgMarkRecur);
        phgPrintf("sentinel: %.8lf  marked: %.8lf  theta: %.8lf\n",
                        (double)flag, (double)sum, (double)gama);
    }
    
    pe *= gama;
    ForAllElements(g, e){
        gama = *DofElementData(osc, e->index);
        if (e->mark == 0 && gama > pe){
            e->mark = depth;
            cc++;
        }
    }

    /* number of elements marked */
stat:
    if (phgVerbosity > 0) {
#if USE_MPI
        MPI_Allreduce(&cc, &count, 1, MPI_LONG, MPI_SUM, g->comm);
#else
        count = cc;
#endif
        phgPrintf("Elements marked for refining: %ld\n", count);
    }
    return;
}

/* equaldistribution strategy 
 * element data */
static void
mark_refine_eqdist(DOF *error, FLOAT theta, FLOAT tol, int depth)
{
    GRID *g = NULL;
    ELEMENT *e;
    FLOAT flag;
    long int cc, count;

    assert(error != NULL);
    assert(error->type->np_elem == 1);
    assert(tol > 0.);
    assert(theta >= 0.);
    assert(depth > 0);
    
    if (phgVerbosity > 0)
        phgPrintf("\n*** EquiDistribution, Element data, Refine\n");
    g = error->g;
    flag = theta * tol;

    cc = 0;
    ForAllElements(g, e) {
        if (*DofElementData(error, e->index) > flag) {
            e->mark = depth;
            cc++;
        }
    }
    if (phgVerbosity > 0) {
#if USE_MPI
        MPI_Allreduce(&cc, &count, 1, MPI_LONG, MPI_SUM, g->comm);
#else
        count = cc;
#endif
        phgPrintf("Elements marked for refining: %ld\n", count);
    }
    return;
}

/* MAX: strategy, element data */
static void
mark_refine_max_elem(DOF *error, FLOAT theta, int depth)
{
    GRID *g = error->g;
    ELEMENT *e;
    FLOAT tmp;
    long int cc, count;
    
    assert(error != NULL);
    assert(error->type->np_elem == 1);
    assert(theta >= 0.);
    assert(depth > 0);

    if (phgVerbosity > 0)
        phgPrintf("\n*** Max, Element data, Refine\n");
    /* max strategy */
    tmp = phgDofNormInftyVec(error) * theta;
    if (tmp <= 0.)
        return;

    if (phgVerbosity > 0)
        phgPrintf("max: %.8lf  sentinel: %.8lf  theta: %.8lf\n",
                (double)phgDofNormInftyVec(error), (double)tmp, (double)theta);
    cc = count = 0;
    ForAllElements(g, e){
        if (*DofElementData(error, e->index) > tmp){
            e->mark = depth;
            cc++;
        }
    }
    if (phgVerbosity > 0) {
#if USE_MPI
        MPI_Allreduce(&cc, &count, 1, MPI_LONG, MPI_SUM, g->comm);
#else
        count = cc;
#endif
        phgPrintf("Elements marked for refining: %ld\n", count);
    }
    return;
}

    



/******************************************************************
 * 
 *                   coarsen subroutine 
 *  if the subroutines are used just followinng the refining 
 *     subroutines, flag should be zero.
 *
 *****************************************************************/
/* GERS: mark elements only uses element data */
static void  
mark_coarsen_gers_elem(DOF *error, FLOAT fac, int depth)
{
    GRID *g = error->g;
    ELEMENT *e;
    FLOAT gama;
    FLOAT sum, flag, pe;
    long cc, count;
    
    assert(error != NULL);
    assert(error->type->np_elem == 1);
    assert(depth < 0);
    assert(fac >= 0. && fac <= 1.);

    if (phgVerbosity > 0)
        phgPrintf("\n*** GERS,Element Data, Coarsen\n");
   
    /* determine factor */
    phgMarkRecur = -1;
    pe = phgDofNormInftyVec(error);
    if (pe <= 0.){
        cc = 0;
        ForAllElements(g, e) {
            e->mark = depth;
            cc++;
        }
        goto stat;
    }

    flag = phgDofNormL1Vec(error) * fac;
    gama = mark_cal_factor(error, pe, flag, 0.0, 1.0, -1, 0);
    if (phgVerbosity > 0){
        sum = part_sum(error, pe, gama, -1, 0);
        phgPrintf("recursion times: %d\n", phgMarkRecur);
        phgPrintf("sentinel: %.8lf  marked: %.8lf  theta: %.8lf\n",
                        (double)flag, (double)sum, (double)gama);
    }
    
    /* mark */
    pe *= gama;
    cc = 0;
    ForAllElements(g, e){
        gama = *DofElementData(error, e->index);
        if (gama <= pe){
            e->mark = depth;
            cc++;
        }
    }
stat:
    if (phgVerbosity > 0) {
#if USE_MPI
        MPI_Allreduce(&cc, &count, 1, MPI_LONG, MPI_SUM, g->comm);
#else
        count = cc;
#endif
        phgPrintf("Elements marked for refining: %ld\n", count);
    }
    return;
}

/* max strategy */
static void 
mark_coarsen_max_elem(DOF *error, FLOAT theta, int depth)
{
    GRID *g = error->g;
    ELEMENT *e;
    FLOAT tmp;
    long int cc, count;
    
    assert(error != NULL);
    assert(error->type->np_elem == 1);
    assert(theta >= 0.);
    assert(depth < 0);
    
    if (phgVerbosity > 0)
        phgPrintf("\n*** Max, Element Data, Coarsen\n");

    /* max strategy */
    tmp = phgDofNormInftyVec(error) * theta;
    cc = count = 0;
    ForAllElements(g, e){
        if (*DofElementData(error, e->index) <= tmp){
            e->mark = depth;
            cc++;
        }
    }
    if (phgVerbosity > 0) {
#if USE_MPI
        MPI_Allreduce(&cc, &count, 1, MPI_LONG, MPI_SUM, g->comm);
#else
        count = cc;
#endif
        phgPrintf("Elements marked for coarsening: %ld\n", count);
    }
    return;
}

/* equaldistribution strategy, element data */
static void
mark_coarsen_eqdist(DOF *error, FLOAT theta, FLOAT tol, int depth)
{
    GRID *g = NULL;
    ELEMENT *e;
    FLOAT flag;
    long int cc, count;

    assert(error != NULL);
    assert(error->type->np_elem == 1);
    assert(tol > 0.);
    assert(theta >= 0.);
    
    if (phgVerbosity > 0)
        phgPrintf("\n*** EquiDistribution, Element data, Refine\n");
    g = error->g;
    flag = theta * tol;

    cc = 0;
    ForAllElements(g, e) {
        if (*DofElementData(error, e->index) <= flag) {
            e->mark = depth;
            cc++;
        }
    }
    if (phgVerbosity > 0) {
#if USE_MPI
        MPI_Allreduce(&cc, &count, 1, MPI_LONG, MPI_SUM, g->comm);
#else
        count = cc;
#endif
        phgPrintf("Elements marked for refining: %ld\n", count);
    }
    return;
}

static const STRATEGY strategy_list[] = {
	MARK_MAX, MARK_GERS, MARK_MNS, MARK_EQDIST, MARK_ALL, MARK_NONE};
static const char *strategy_names[] = {
	"max",    "gers",    "mns",    "eqdist",    "all",    "none", NULL};
static int default_refine_strategy = 1;
static int default_coarsen_strategy = 1;

static const char *
default_strategy(OPTION *o, const char *arg, int prefix)
{
    const char **p;

    Unused(prefix);

    if (arg == NULL) {
	/* print help about this option */
#if 0
	phgPrintf("Valid keywords are: ");
	for (p = strategy_names; *p != NULL; p++)
	    phgPrintf("%s\"%s\"", p == strategy_names ? "":", ", *p);
	phgPrintf("\n");
#endif
    }
    else if (arg[0] != '\0' && strcmp(arg, "none")) {
	/* Note: when arg[0]=='\0' old value is kept */
	for (p = strategy_names; *p != NULL; p++) {
	    if (strcmp(*p, arg) != 0)
		continue;
	    default_refine_strategy = default_coarsen_strategy =
						(int)(p - strategy_names);
	    break;
	}
	if (*p == NULL)
	    phgError(1, "invalid argument for the \"-strategy\" option: %s\n",
			arg); 
    }

    return default_refine_strategy == default_coarsen_strategy ?
		strategy_names[default_refine_strategy] : "none";
}

void
phgMarkElements(STRATEGY strategy, DOF *error, FLOAT theta, DOF *osc,
			FLOAT zeta, int depth, 
		STRATEGY coarsen_strategy, DOF *coarsen_error, FLOAT gamma,
			int coarsen_depth, FLOAT tol)
{
    static BOOLEAN initialized = FALSE;

    GRID *g = NULL;
    ELEMENT *e;
    
    if (!initialized){
        initialized = TRUE;
        phgOptionsRegisterTitle("\nMarking strategy options:", "\n", "marking");
        phgOptionsRegisterHandler("-strategy",
				  "Default refinement and coarsening strategy",
				  default_strategy, FALSE);
        phgOptionsRegisterKeyword("-refine_strategy",
				  "Default refinement strategy",
				  strategy_names, &default_refine_strategy);
        phgOptionsRegisterKeyword("-coarsen_strategy",
				  "Default coarsening strategy",
				  strategy_names, &default_coarsen_strategy);
        return;
    }

    assert(tol >= 0.);
    /* assignment to g */
    if (error != NULL) {
        g = error->g;
    }
    else if (osc != NULL){
        g = osc->g;
    }
    else if (coarsen_error != NULL) {
        g = coarsen_error->g;
    }
    else {
        phgError(1, "no valid DOF pass in\n");
    }

    /*  choose default strategy */
    if (strategy == MARK_DEFAULT)
        strategy = strategy_list[default_refine_strategy];

    /* verify error */
    if (!(strategy == MARK_ALL || strategy == MARK_NONE)) {
        if (error == NULL)
            phgError(1, "Error Indicator shouldn't be NULL!\n");
    }

    /* verify refine level */
    if (strategy != MARK_NONE) {
        assert(depth > 0);
    }

    assert(strategy != MARK_NONE || (theta >= 0. && theta <= 1.));
    assert(strategy != MARK_NONE || (zeta >= 0. && zeta <= 1.));

    /* mark elements for refining */
    /* init grid */
    mark_init(g);
    switch (strategy) {
        case MARK_ALL:
            ForAllElements(g, e)
                e->mark = depth;
            return;
            break;
        case MARK_MAX:
            mark_refine_max_elem(error, theta, depth);
            if (osc != NULL)
                 mark_refine_max_elem(osc, zeta, depth);
            break;
        case MARK_MNS:
            mark_refine_mns(error, theta, osc, zeta, depth);
            break;
        case MARK_GERS:
            if (error->type->np_elem == 1) { /* element data */
                if (osc == NULL)
                    mark_refine_gers_elem(error, theta, depth);
                else
                    mark_refine_gers_osc(error, theta, osc, zeta, depth);
            }
            else if (error->type->np_face == 1){  /* face data */
                if (osc == NULL) {
                    mark_refine_gers_face(error, theta, depth);
                }
                else {
                    mark_refine_mns(error, theta, osc, zeta, depth);
                }
            }
            else {
                phgError(1, "Wrong DOF type %s %d \n", __FILE__, __LINE__);
            }
            break;
        case MARK_EQDIST:
            mark_refine_eqdist(error, theta, tol, depth);
            if (osc != NULL) {
                mark_refine_eqdist(osc, zeta, tol, depth);
            }
            break;
        case MARK_NONE:
            break;
        default:
            phgError(1, "phgMarkElements: invalid refine strategy %d!\n",
                                                                strategy);
    }
    
    /* mark elements for coarsening */
    /* choose default strategy */
    if (coarsen_strategy == MARK_DEFAULT)
        coarsen_strategy = strategy_list[default_coarsen_strategy];

    if (coarsen_strategy == MARK_NONE)
	return;

    /* verify coarsen-error */
    if (coarsen_strategy != MARK_ALL) {
        if (coarsen_error == NULL) 
            phgError(1, "Coarsen error indicator shouldn't be NULL!\n");
        else if (coarsen_error->type->np_elem != 1) {
            phgError(1, "Wrong DOF type coarsen_error->type->np_elem != 1\n");
        }
    }
    /* verify parameters */
    assert(coarsen_depth > 0);
    assert(gamma >= 0. && gamma <= 1.);

    switch (coarsen_strategy) {
        case MARK_MAX:
            mark_coarsen_max_elem(coarsen_error, gamma, -coarsen_depth);
            break;
        case MARK_GERS:
            mark_coarsen_gers_elem(coarsen_error, gamma, -coarsen_depth);
            break;
        case MARK_EQDIST:
            mark_coarsen_eqdist(coarsen_error, gamma, tol, -coarsen_depth);
            break;
        case MARK_ALL:
            ForAllElements(g, e) {
                e->mark = -coarsen_depth;
            }
            break;
        case MARK_NONE:
            break;
        default:
            phgError(1,"phgMarkElements: invalid coarsen strategy %d!\n",
                                                        coarsen_strategy);
    }
    return;
}
