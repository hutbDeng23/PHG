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

/* $Id: halo.h,v 1.7 2022/07/30 02:33:25 zlb Exp $ */

#ifndef PHG_HALO_H

struct HALO_ {
    /* elocal[] stores list of elements for the halo of remote processes
     * matching scnt[] and sdsps, while elist[] stores list of halo elements
     * matching rcnt[] and rdsp[]. They can be used for exchaging halo data. */
    int		*scnt, *rcnt, *sdsp, *rdsp;
    ELEMENT	**elocal;	/* list of local elements for remote halos */
    ELEMENT	*elist;		/* list of halo elements */
    INT		*map;		/* neighbours.list[i] => halo->elist[map[i]] */
    COORD	*verts;		/* coordinates of remote vertices */
    INT		nvert_bak;	/* saved number of vertices without halo */
    INT 	nedge_bak;	/* saved number of edges without halo */
    INT 	nface_bak;	/* saved number of faces without halo */
    INT 	nelem_bak;	/* saved number of elems without halo */
    int type;			/* halo type */
};

#define HaloType(size, kind)	((size) * 4 + (kind))
#define GetHaloKind(type)	((type) % 4)
#define GetHaloSize(type)	((type) / 4 + 1)
#define HALO_ALL		0
#define HALO_VERT		1
#define HALO_EDGE		2
#define HALO_FACE		3

#ifdef __cplusplus
extern "C" {
#endif

ELEMENT *phgGetNeighbour_(GRID *g, ELEMENT *e, int face,
					const char *file, int line);
#define phgGetNeighbour(g, e, face) \
	phgGetNeighbour_(g, e, face, __FILE__, __LINE__)

void phgDofUpdateGhost_(DOF *dof, BOOLEAN flag);
#define phgDofUpdateGhost(dof)	phgDofUpdateGhost_(dof, TRUE)

void phgSetupHalo_(GRID *g, int halo_type, BOOLEAN update_dofs);
#define phgSetupHalo(g, halo_type) phgSetupHalo_(g, halo_type, TRUE)
void phgRemoveHalo(GRID *g);

#ifdef __cplusplus
}
#endif

#define PHG_HALO_H
#endif
