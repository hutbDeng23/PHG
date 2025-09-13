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

/* $Id: elem-info.c,v 1.10 2021/11/29 01:25:39 zlb Exp $ */

/* This file defines element types and functions for tetrahedralizing them
 *
 * TODO:
 *
 *  1. use edges of tetrahedra to constrain splitting of quadrilaterals
 *
 *  2. make use of geometrical criteria to generate a better mesh, cf.,
 *	Thomas Apel \& Nico D\"uvelmeyer,
 *	Transformation of hexahedral finite element meshes into
 *	tetrahedral meshes according to quality criteria, 
 *	Computing,  Volume 71, Issue 4, 2003, 293--304,
 *	Preprint: SFB393/03-09.
 *	File: bibliography/sfb03-09-HexahedralToTetrahedral.pdf
 */

#include <stdio.h>
#ifdef __PHG__
# include "phg/config.h"
#endif	/* defined(__PHG__) */
#include "phg/elem-info.h"

#ifndef DEBUG
# define DEBUG 1
#endif

/*========================= List of element types  ========================*/

/*********************************************************************
 * Numbering of the vertices of a tetrahedron is irrelevant.
 **********************************************************************/
static EDGE_INFO tetra_edges[] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
static FACE_INFO tetra_faces[] = {
    {3, 1, 2, 3, -1},		/* the face opposite to vertex 0 */
    {3, 0, 2, 3, -1},		/* the face opposite to vertex 1 */
    {3, 0, 1, 3, -1},		/* the face opposite to vertex 2 */
    {3, 0, 1, 2, -1}		/* the face opposite to vertex 3 */
};

/* Array mapping edge no to vertices, used by phgGetEdgeVertex macro */
char phgTetEdgeVertex[6][2] = {{0,1}, {0,2}, {0,3}, {1,2}, {1,3}, {2,3}};
/* Array mapping face no to vertices, used by GetFaceVertex macro */
char phgTetFaceVertex[4][3] = {{1,2,3}, {0,2,3}, {0,1,3}, {0,1,2}};
/* array mapping pairs of vertices to edges */
char phgTetEdgeNo[4][4] = {{-1, 0, 1, 2},
			   { 0,-1, 3, 4},
			   { 1, 3,-1, 5},
			   { 2, 4, 5,-1}};

#define NETet	(sizeof(tetra_edges) / sizeof(tetra_edges[0]))
#define NFTet	(sizeof(tetra_faces) / sizeof(tetra_faces[0]))

/*********************************************************************
 * Numbering of the vertices of a pyramid, view from above:
 *
 *    3	*-------* 2
 *	|	|		4
 *	|	|		*
 *	|	|	    top vertex
 *    0	*-------* 1
 *	bottom face
 *
 **********************************************************************/
static EDGE_INFO pyramid_edges[] = {
    {0,1}, {1,2}, {2,3}, {3,0}, {0,4}, {1,4}, {2,4}, {3,4}
};
static FACE_INFO pyramid_faces[] = {
    /* The faces are ordered as: side faces 0-1, 1-2, 2-3, 3-0, bottom face
     * This ordering facilitates face mapping after rotating the vertices
     * of the bottom face in phgPyramid2Tetra */
    {3, 0, 1, 4, -1},		/* the side face containing vertices 0,1 */
    {3, 1, 2, 4, -1},		/* the side face containing vertices 1,2 */
    {3, 2, 3, 4, -1},		/* the side face containing vertices 2,3 */
    {3, 3, 0, 4, -1},		/* the side face containing vertices 3,0 */
    {4, 0, 1, 2, 3},		/* the bottom face */
};

#define NEPyr	(sizeof(pyramid_edges) / sizeof(pyramid_edges[0]))
#define NFPyr	(sizeof(pyramid_faces) / sizeof(pyramid_faces[0]))

/*********************************************************************
 * Numbering of the vertices of a prism (or a wedge), view from above:
 *
 *    2	*-------* 1	      5	*-------* 4
 *	 \     /		 \     /
 *	  \   /			  \   /
 *	   \ /			   \ /
 *	    *			    *
 *	    0			    3
 *	bottom face		top face
 *
 **********************************************************************/
static EDGE_INFO prism_edges[] = {
    {0,1}, {1,2}, {2,0}, {3,4}, {4,5}, {5,3}, {0,3}, {1,4}, {2,5}
};
static FACE_INFO prism_faces[] = {
    /* The faces are ordered as: side faces 0-1, 1-2, 2-0, then bottom, top
     * This ordering facilitates face mapping after rotating the vertices
     * in phgPrism2Tetra */
    {4, 0, 1, 4, 3},		/* the side face containing vertices 0,1 */
    {4, 1, 2, 5, 4},		/* the side face containing vertices 1,2 */
    {4, 2, 0, 3, 5},		/* the side face containing vertices 2,0 */
    {3, 0, 1, 2, -1},		/* the bottom face */
    {3, 3, 4, 5, -1},		/* the top face */
};
#define NEPri	(sizeof(prism_edges) / sizeof(prism_edges[0]))
#define NFPri	(sizeof(prism_faces) / sizeof(prism_faces[0]))

/**********************************************************************
 * Numbering of the vertices of a hexahedron (or a brick), view from above:
 *
 *    3	*-------* 2	      7	*-------* 6
 *	|       |		|	|
 *	|	|		|	|
 *	|  	|		|	|
 *    0	*-------* 1	      4	*-------* 5
 *	bottom face		top face
 *
 **********************************************************************/
static EDGE_INFO hexa_edges[] = {
    {0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7}, {7,4},
    {0,4}, {1,5}, {2,6}, {3,7}
};
static FACE_INFO hexa_faces[] = {
    /* The faces are ordered as: side faces 0-1, 1-2, 2-3, 3-0, bottom, top
     * This ordering facilitates face mapping after rotating the vertices
     * in phgHexa2Tetra */
    {4, 0, 1, 5, 4},		/* the side face containing vertices 0,1 */
    {4, 1, 2, 6, 5},		/* the side face containing vertices 1,2 */
    {4, 2, 3, 7, 6},		/* the side face containing vertices 2,3 */
    {4, 3, 0, 4, 7},		/* the side face containing vertices 3,0 */
    {4, 0, 1, 2, 3},		/* the bottom face */
    {4, 4, 5, 6, 7},		/* the top face */
};

/* Array mapping edge no to vertices, used by phgGetEdgeVertex macro */
char phgHexEdgeVertex[12][2] = {
    {0,1}, {1,2}, {2,3}, {3,0},
    {4,5}, {5,6}, {6,7}, {7,4},
    {0,4}, {1,5}, {2,6}, {3,7}
};

/* Array mapping face no to vertices, used by GetFaceVertex macro */
char phgHexFaceVertex[6][4] = {
    {0, 1, 5, 4},		/* the side face containing vertices 0,1 */
    {1, 2, 6, 5},		/* the side face containing vertices 1,2 */
    {2, 3, 7, 6},		/* the side face containing vertices 2,3 */
    {3, 0, 4, 7},		/* the side face containing vertices 3,0 */
    {0, 1, 2, 3},		/* the bottom face */
    {4, 5, 6, 7},		/* the top face */
};

/* Array mapping face no to vertices, used by GetFaceVertex macro */
char phgHexFaceEdge[6][4] = {
    {0, 4, 8,   9},		/* the side face containing vertices 0,1 */
    {1, 5, 9,  10},		/* the side face containing vertices 1,2 */
    {2, 6, 10, 11},		/* the side face containing vertices 2,3 */
    {3, 7, 11,  8},		/* the side face containing vertices 3,0 */
    {0, 1, 2, 3},		/* the bottom face */
    {4, 5, 6, 7},		/* the top face */
};

/* array mapping pairs of vertices to edges */
char phgHexEdgeNo[8][8] = {
    {-1,  0, -1,  3,  8, -1, -1, -1},  
    { 0, -1,  1, -1, -1,  9, -1, -1}, 
    {-1,  1, -1,  2, -1, -1, 10, -1}, 
    { 3, -1,  2, -1, -1, -1, -1, 11}, 
    { 8, -1, -1, -1, -1,  4, -1,  7}, 
    {-1,  9, -1, -1,  4, -1,  5, -1}, 
    {-1, -1, 10, -1, -1,  5, -1,  6}, 
    {-1, -1, -1, 11,  7, -1,  6, -1}
};

char phgHexEdgeFace[12][2] = {
    {0, 4}, {1, 4}, {2, 4}, {3, 4}, 
    {0, 5}, {1, 5}, {2, 5}, {3, 5}, 
    {3, 0}, {0, 1}, {1, 2}, {2, 3}
};

#define NEBri	(sizeof(hexa_edges) / sizeof(hexa_edges[0]))
#define NFBri	(sizeof(hexa_faces) / sizeof(hexa_faces[0]))

ELEM_INFO phgElemInfo[] = {
    /* name, medit name, list of edges, list of faces, func2tet,
     * dim, nvert, nedge, nface */
    {"vertex",	"Vertices",	NULL,		NULL,		NULL,
				0, 1, 0, 0},
    {"edge",	NULL/*Edges*/,	NULL,		NULL,		NULL,
				1, 2, 1, 0},
    {"triangle", "Triangles",	NULL,		NULL,		NULL,
				2, 3, 3, 1},
    {"quadrilateral", "Quadrilaterals", NULL,	NULL,		NULL,
				2, 4, 4, 1},
    {"tetragedron", "Tetrahedra", tetra_edges,	tetra_faces,	NULL,
				3, 4, NETet, NFTet},
    {"pyramid",	NULL,		pyramid_edges,	pyramid_faces, phgPyramid2Tetra,
				3, 5, NEPyr, NFPyr},
    {"prism",	NULL,		prism_edges,	prism_faces,	phgPrism2Tetra,
				3, 6, NEPri, NFPri},
    {"hexahedron", "Hexahedra",	hexa_edges,	hexa_faces,	phgHexa2Tetra,
				3, 8, NEBri, NFBri}
};

/* Global variable containing number of element types */
int phgElemInfoCount = (sizeof(phgElemInfo) / sizeof(phgElemInfo[0]));

/*========================= Exported functions ========================*/

void
phgPrism2Tetra(int verts[], ADD_TETRA_FUNC new_tetra, int bound_type[])
/* splits a prism into 3 tetrahedra. */
{
    int i, j, bc[4], *bt = (bound_type == NULL ? NULL : bc);
    short a0, a1, a2, a3, a4, a5;

#if 0
    /* check for duplicate vertices */
    for (i = 0; i < phgElemInfo[ET_PRISM].nedge; i++) {
	if (verts[phgElemInfo[ET_PRISM].edge_info[i][0]] ==
	    verts[phgElemInfo[ET_PRISM].edge_info[i][1]]) {
	    fprintf(stderr, "Warning: removing invalid prism: "
			    "%d %d %d %d %d %d\n", verts[0], verts[1],
			    verts[2], verts[3], verts[4], verts[5]);
	    return;
	}
    }
#endif

/* Face map after rotating the vertices */
#define Bottom	(a0 < 3 ? 3 : 4)
#define Top	(a0 < 3 ? 4 : 3)
#define F(i)	bound_type[i < 3 ? (a0 + i) % 3 : (i == 3 ? Bottom : Top)]

    /* a0 = the vertex with minimum index */
    j = verts[0];
    a0 = 0;
    for (i = 1; i < 6; i++) {
	if (verts[i] < j) {
	    a0 = i;
	    j = verts[i];
	}
    }

    /* renumber the vertices as {a0, a1, a2, a3, a4, a5} */
    if (a0 >= 3) {
	if ((a1 = a0 + 1) > 5)
	    a1 = 3;
	if ((a2 = a1 + 1) > 5)
	    a2 = 3;
	a3 = a0 - 3;
	a4 = a1 - 3;
	a5 = a2 - 3;
    }
    else {
	if ((a1 = a0 + 1) > 2)
	    a1 = 0;
	if ((a2 = a1 + 1) > 2)
	    a2 = 0;
	a3 = a0 + 3;
	a4 = a1 + 3;
	a5 = a2 + 3;
    }

    /* determine the diagonal of the side face not containing {a0,a3} */
    i = verts[a1];
    if (i > verts[a5])
	i = verts[a5];
    j = verts[a4];
    if (j > verts[a2])
	j = verts[a2];

    /* 1st tetrahedron */
    if (bound_type != NULL) {
	bc[0] = F(4);
	bc[1] = -1;
	bc[2] = F(2);
	bc[3] = F(0);
    }
    (*new_tetra)(verts[a0], verts[a3], verts[a4], verts[a5], bt);

    if (i < j) {
	/* 2nd tetrahedron */
	if (bound_type != NULL) {
	    bc[0] = F(1);
	    bc[1] = -1;
	    bc[2] = -1;
	    bc[3] = F(0);
	}
	(*new_tetra)(verts[a0], verts[a1], verts[a4], verts[a5], bt);

	/* 3rd tetrahedron */
	if (bound_type != NULL) {
	    bc[0] = F(1);
	    bc[1] = -1;
	    bc[2] = F(3);
	    bc[3] = F(2);
	}
	(*new_tetra)(verts[a0], verts[a2], verts[a5], verts[a1], bt);
    }
    else {
	/* 2nd tetrahedron */
	if (bound_type != NULL) {
	    bc[0] = F(1);
	    bc[1] = -1;
	    bc[2] = F(3);
	    bc[3] = F(0);
	}
	(*new_tetra)(verts[a0], verts[a1], verts[a4], verts[a2], bt);

	/* 3rd tetrahedron */
	if (bound_type != NULL) {
	    bc[0] = F(1);
	    bc[1] = -1;
	    bc[2] = -1;
	    bc[3] = F(2);
	}
	(*new_tetra)(verts[a0], verts[a2], verts[a5], verts[a4], bt);
    }

#undef Top
#undef Bottom
#undef F
}

void
phgHexa2Tetra(int verts[], ADD_TETRA_FUNC new_tetra, int bound_type[])
/* splits a hexahedron into 5 or 6 tetrahedra.
 * Note: for a given splitting pattern of the 6 faces, the hexahedron can be
 * conformingly split into 5 or 6 tetrahedra iff there exist 3 face diagonals
 * meeting at a same vertex. */
{
    int i, j, bc[5], *bt = (bound_type == NULL ? NULL : bc);
    short a0, a[8];
    int v[6];

#if 0
    /* check for duplicate vertices */
    for (i = 0; i < phgElemInfo[ET_HEXA].nedge; i++) {
	if (verts[phgElemInfo[ET_HEXA].edge_info[i][0]] ==
	    verts[phgElemInfo[ET_HEXA].edge_info[i][1]]) {
	    fprintf(stderr, "Warning: deleting invalid hexahedron: "
			    "%d %d %d %d %d %d %d %d\n", verts[0], verts[1],
                            verts[2], verts[3], verts[4], verts[5], verts[6],
			    verts[7]);
	    return;
	}
    }
#endif

/* Face map after rotating the vertices */
#define Bottom	(a0 < 4 ? 4 : 5)
#define Top	(a0 < 4 ? 5 : 4)
#define F(i)	bound_type[i < 4 ? (a0 + i) % 4 : (i == 4 ? Bottom: Top)]

    /* a0 = a[0] = the vertex with minimum index */
    j = verts[0];
    a0 = 0;
    for (i = 1; i < 8; i++) {
	if (verts[i] < j) {
	    a0 = i;
	    j = verts[i];
	}
    }
    a[0] = a0;

    /* renumber the vertices as {a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]} */
    if (a0 >= 4) {
	a[4] = a[0] - 4;
	for (i = 1; i < 4; i++) {
	    if ((a[i] = a[i - 1] + 1) > 7)
		a[i] = 4;
	    a[4 + i] = a[i] - 4;
	}
    }
    else {
	a[4] = a[0] + 4;
	for (i = 1; i < 4; i++) {
	    if ((a[i] = a[i - 1] + 1) > 3)
		a[i] = 0;
	    a[4 + i] = a[i] + 4;
	}
    }

    /* if the face diagonals of a pair of faces are parallel, then split the
     * hexahedron into two prisms, otherwise split it into 5 tetrahedra */

/* for face-pair {v0,v1,v2,v3}-{u0,u1,u2,u3}, check if min(u0,u2)<min(u1,u3) */
#define CHECK_FACE_PAIR(v0, v1, v2, v3, u0, u1, u2, u3,		\
		  b0, b1, b2, b3, b4, b5)			\
    i = verts[u0];						\
    if (i > verts[u2])						\
	i = verts[u2];						\
    j = verts[u1];						\
    if (j > verts[u3])						\
	j = verts[u3];						\
    if (i < j) {						\
	if (bound_type != NULL) {				\
	    bc[0] = -1;						\
	    bc[1] = b2;						\
	    bc[2] = b3;						\
	    bc[3] = b4;						\
	    bc[4] = b5;						\
	}							\
	v[0] = verts[v0]; v[1] = verts[v2]; v[2] = verts[v3];	\
	v[3] = verts[u0]; v[4] = verts[u2]; v[5] = verts[u3];	\
	phgPrism2Tetra(v, new_tetra, bt);			\
	if (bound_type != NULL) {				\
	    bc[0] = b0;						\
	    bc[1] = b1;						\
	    bc[2] = -1;						\
	    bc[3] = b4;						\
	    bc[4] = b5;						\
	}							\
	v[0] = verts[v0]; v[1] = verts[v1]; v[2] = verts[v2];	\
	v[3] = verts[u0]; v[4] = verts[u1]; v[5] = verts[u2];	\
	phgPrism2Tetra(v, new_tetra, bt);			\
	return;							\
    }

    CHECK_FACE_PAIR(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7],
		    F(0), F(1), F(2), F(3), F(4), F(5))
    CHECK_FACE_PAIR(a[0], a[1], a[5], a[4], a[3], a[2], a[6], a[7],
		    F(4), F(1), F(5), F(3), F(0), F(2))
    CHECK_FACE_PAIR(a[0], a[4], a[7], a[3], a[1], a[5], a[6], a[2],
		    F(0), F(5), F(2), F(4), F(3), F(1))

    /* no parallel face diagonals, split the cube into 5 tetrahedra */
    if (bound_type != NULL) {
	bc[0] = F(5);
	bc[1] = -1;
	bc[2] = F(3);
	bc[3] = F(0);
    }
    (*new_tetra)(verts[a[0]], verts[a[4]], verts[a[5]], verts[a[7]], bt);
    if (bound_type != NULL) {
	bc[0] = F(2);
	bc[1] = F(3);
	bc[2] = -1;
	bc[3] = F(4);
    }
    (*new_tetra)(verts[a[0]], verts[a[2]], verts[a[3]], verts[a[7]], bt);
    if (bound_type != NULL) {
	bc[0] = F(1);
	bc[1] = -1;
	bc[2] = F(0);
	bc[3] = F(4);
    }
    (*new_tetra)(verts[a[0]], verts[a[1]], verts[a[2]], verts[a[5]], bt);
    if (bound_type != NULL) {
	bc[0] = F(5);
	bc[1] = F(2);
	bc[2] = -1;
	bc[3] = F(1);
    }
    (*new_tetra)(verts[a[2]], verts[a[5]], verts[a[6]], verts[a[7]], bt);
    if (bound_type != NULL) {
	bc[0] = -1;
	bc[1] = -1;
	bc[2] = -1;
	bc[3] = -1;
    }
    (*new_tetra)(verts[a[0]], verts[a[2]], verts[a[5]], verts[a[7]], bt);

#undef Top
#undef Bottom
#undef F
}

void
phgPyramid2Tetra(int verts[], ADD_TETRA_FUNC new_tetra, int bound_type[])
/* splits a pyramid into 2 tetrahedra. */
{
    int i, j, bc[4], *bt = (bound_type == NULL ? NULL : bc);
    short a0, a1, a2, a3;

    /* a0 = the vertex on the bottom face with minimum index */
    j = verts[0];
    a0 = 0;
    for (i = 1; i < 4; i++) {
	if (verts[i] < j) {
	    a0 = i;
	    j = verts[i];
	}
    }

    a1 = (a0 + 1) % 4;
    a2 = (a1 + 1) % 4;
    a3 = (a2 + 1) % 4;

    if (bound_type != NULL) {
	bc[0] = bound_type[a1];
	bc[1] = -1;
	bc[2] = bound_type[a0];
	bc[3] = bound_type[4];
    }
    (*new_tetra)(verts[a0], verts[a1], verts[a2], verts[4], bt);
    if (bound_type != NULL) {
	bc[0] = bound_type[a2];
	bc[2] = bound_type[a3];
    }
    (*new_tetra)(verts[a0], verts[a3], verts[a2], verts[4], bt);
}

#if TEST
int
main(int argc, char *argv[])
{
    int i, j, k;
    ELEM_INFO *ei;

    for (i = 0, ei = phgElemInfo; i < phgElemInfoCount; i++, ei++) {
	printf("Name: %-8s  Nvert: %d  Nedge: %-2d  Nface: %d\n",
		ei->name, ei->nvert, ei->nedge, ei->nface);
	printf("%16sEdges:", "");
	for (j = 0; j < ei->nedge; j++)
	    printf(" %d-%d", ei->edge_info[j][0], ei->edge_info[j][1]);
	printf("\n%16sFaces:", "");
	for (j = 0; j < ei->nface; j++) {
	    for (k = 0; k < ei->face_info[j][0]; k++)
		printf("%c%d", k == 0 ? ' ' : '-', ei->face_info[j][k + 1]);
	}
	printf("\n");
    }
}
#endif	/* TEST */
