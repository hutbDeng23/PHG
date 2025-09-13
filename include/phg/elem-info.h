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

/* $Id: elem-info.h,v 1.5 2022/01/13 05:26:41 zlb Exp $ */

#ifndef PHG_ELEM_INFO_H

typedef void (*ADD_TETRA_FUNC)
			(int v0, int v1, int v2, int v3, int bound_type[]);

/* Note: FACE_INFO[0] = number of nodes, FACE_INFO[1:] = node list.
 *
 * !!! IMPORTANT: the vertices must be listed circularly which is required by
 * 		  update_neighbours()
 */
typedef char EDGE_INFO[2];	/* list of vertices */
typedef char FACE_INFO[5];	/* 0: # of vertices, 1..: list of vertices */
typedef void (*FUNC2TET)(int verts[], ADD_TETRA_FUNC func, int bound_type[]);

typedef struct {
    const char	*name;		/* name of the element type */
    const char	*medit_name;	/* name of the Medit element type */
    EDGE_INFO	*edge_info;	/* list of edges */
    FACE_INFO	*face_info;	/* list of faces */
    FUNC2TET    func2tet;	/* function to covert to tet */
    char	dim;		/* space dimension (0, 1, 2, 3) */
    char	nvert;		/* number of vertices */
    char	nedge;		/* number of edges */
    char	nface;		/* number of faces */
} ELEM_INFO;

/* element types (note: must match "phgElemInfo[]" in order!) */
typedef enum {
    ET_POINT	= 0,
    ET_EDGE	= 1,
    ET_TRIANGLE	= 2,
    ET_QUAD	= 3,
    ET_TETRA	= 4,
    ET_PYRAMID	= 5,
    ET_PRISM	= 6,
    ET_HEXA	= 7,
    /* Warning: no corresponding entries in "phgElemInfo[]" for the keywords
     * below, they are only used by the "quad-interface*.c*" files. */
    ET_RECTANGLE = 8,
    ET_CUBOID	= 9,
    ET_MIXED	= 10
} ELEM_TYPE;

extern int phgElemInfoCount;
extern ELEM_INFO phgElemInfo[];

#ifdef __cplusplus
extern "C" {
#endif

void phgPrism2Tetra(int verts[], ADD_TETRA_FUNC func, int bound_type[]);
void phgHexa2Tetra(int verts[], ADD_TETRA_FUNC func, int bound_type[]);
void phgPyramid2Tetra(int verts[], ADD_TETRA_FUNC func, int bound_type[]);

#ifdef __cplusplus
}
#endif

#define PHG_ELEM_INFO_H
#endif
