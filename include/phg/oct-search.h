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

/*
 * Point-location utilities based on oct-tree search.
 *
 * Author: LENG Wei (wleng@lsec.cc.ac.cn)
 * $Id: oct-search.h,v 1.4 2017/03/03 06:34:31 zlb Exp $
 * */

#ifndef PHG_OCT_SEARCH_H	/* whole file */

#define OCT_TREE_EPS	1e-13

/**********************************/
/* Structure & function prototype */
/**********************************/

typedef struct CUBE_ {
    struct CUBE_ *children[2 * 2 * 2];		/**< Pointers to children */
    struct CUBE_ *father;	/* father */
    int level;			/* oct tree level */
    FLOAT dh;			/* element size */
    FLOAT bbox[Dim][2];		/* bounding box */
    FLOAT midx[Dim];		/* mid point */
    FLOAT diam;

    INT ntet, size;
    ELEMENT **liste;
} CUBE;

typedef struct OCT_TREE_ {
    int nx, ny, nz;
    int ncube;
    int nroot;
    int ncube_list;
    FLOAT dx, dy, dz;
    CUBE *roots;
    FLOAT bbox[Dim][2];		/* bounding box */
    FLOAT (*bboxs)[Dim][2];

    GRID *g;			/* working tree */
    /* COORD *tet_center;               /\* extra tets info *\/ */
    /* FLOAT *tet_inner_diam;   /\* extra tets info *\/ */
    CUBE **cube_list;		/* record new cubes alloc */
} OCT_TREE;

/* oct tree search */
OCT_TREE *phgOctTreeBuild(GRID *g);
BOOLEAN phgOctTreeLocatePoint(OCT_TREE *og, const FLOAT *pt, ELEMENT **e_ptr,
		FLOAT *lambda);
void phgOctTreeDestroy(OCT_TREE **og_ptr);

#define PHG_OCT_SEARCH_H
#endif
