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

/* $Id: grid.h,v 1.4 2017/03/02 08:06:21 zlb Exp $ */

#ifndef PHG_GRID_H	/* whole file */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SIMPLE_TET_ {
    INT verts[4];
    BTYPE bound_type[4];
    SHORT region_mark;
} SIMPLE_TET;


typedef BOOLEAN (*SUBREGION_FUNC)(GRID *g, SIMPLEX *e);

GRID *phgImportParallelGrid(
	    INT nvert, INT nelem, INT nvert_global, INT nelem_global,
	    INT *L2Gmap_vert, INT *L2Gmap_elem, FLOAT *coord, SIMPLE_TET *tet,
	    MPI_Comm comm);

GRID *phgGetSubGrid2(GRID *g, SUBREGION_FUNC subregion_func);
GRID *phgGetSubGrid(GRID *g, BOOLEAN *elem_marker);

#ifdef __cplusplus
}
#endif

#define PHG_GRID_H
#endif	/* !defined(PHG_GRID_H) */
