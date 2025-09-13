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

/* $Id: periodic.h,v 1.1 2014/04/04 04:30:01 zlb Exp $ */

#ifndef PHG_PERIODIC_H

/**
 @brief The PERIOD struct

  Note: for a grid 'g' with periodic boundaries:

	- g->nvert, g->nvert_global, and g->verts are the same as without
	  periodicity.

	- g->period->L2Gmap_vert[] is an array of size g->nvert which maps
	  a local vertex to its global index, taking into account the
	  periodicity (i.e., vertices which are periodic to each other are
	  mapped to a same index), and g->period->nvert_global is the total
	  number of indices (g->period->nvert_global < g->nvert_global).

	- all other members of g are constructed by regarding periodic
	  vertices/edges/faces as identical (thus g->nedge and g->nface
	  are normally smaller than without periodicity).
 */
typedef struct PERIOD_ {
    FLOAT	*directions;	/**< periodic directions */
    INT		*L2Gmap_vert;	/**< global indices under periodicity map */
    INT		*ordering;	/**< L2Gmap_vert[ordering[]] is increasing */
    INT		nvert_global;	/**< total # vertices */
    int		flags;		/**< combination of [XYZ]_MASK bits */
} PERIOD;

#define X_MASK  1
#define Y_MASK  2
#define Z_MASK  4

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__ICC) /*&& (defined(__MIC__) || __ICC <= 1110)*/
/* workaround for a possible icc <= 11.1 + OpenMP+MPI (xeon) or
 * icc 13.1.3 (mic) bug */
INT _phg_global_vertex_p(struct GRID_ *g, INT no);
# define GlobalVertexP(g, no)	_phg_global_vertex_p(g, no)
#else	/* defined(__ICC) && (defined(__MIC__) || __ICC <= 1110) */
# define GlobalVertexP(g, no) \
    ((g)->period == NULL /*|| (g)->period->L2Gmap_vert == NULL*/ ? \
	GlobalVertex(g, no) : (g)->period->L2Gmap_vert[no])
#endif	/* defined(__ICC) ... */

void phgSetPeriodicity(GRID *g, int flags);
void phgSetPeriodicDirections(GRID *g, const FLOAT dirs[]);

void phgPeriodFree(GRID *g);
void phgPeriodInit(GRID *g);

#ifdef __cplusplus
}
#endif

#define PHG_PERIODIC_H
#endif
