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

/* $Id: quad-interface.h,v 1.68 2021/12/24 06:27:46 zlb Exp $
 *
 * Numerical quadrature on an element containing a smooth interface */

#ifndef PHG_QUAD_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct QI_OPTIONS_ {
    FLOAT EPS;
    FLOAT threshold;	/* max a(w) allowed */

    /* whether store in the rule data unit normal vectors of the interface
     * at the quadrature points (only used for surface rules) */
    BOOLEAN nv_flag;

    /* whether use one of the face normals as the direction nu */
    BOOLEAN use_face_normal;

    /* The following flags control using shortcuts (faster but less accurate) */
    BOOLEAN shortcut_t;
    BOOLEAN shortcut_s;

    int adapt;

    int subdiv_type;
    INT subdiv_limit;	/* maximum recursion depth */
    INT subdiv_level0;	/* initial subdivision level */

    int newton;
    INT newton_maxits;	/* maximum Newton iterations */
    INT newton_porder;	/* polynomial order for initial solution */

    BOOLEAN rect2tri;
    BOOLEAN cuboid2tetra;

    INT quad_n1;	/* nbr of intervals for get_quad_1D */
    INT quad_n2;	/* nbr of intervals/direction for get_quad_2D */

    /* debuging options */
    INT graph_t;
    INT graph_s;
    INT graph_w;

    FLOAT dbg_w;	/* <-1: auto, [-1,1]: specified, >1.0: random */
    INT   dbg_wN;	/* nits for 0.618 search, 0: middle point */

    FLOAT dbg_value_t; /* in [0,1] => debug with this t value */
    FLOAT dbg_value_s; /* in [0,1] => debug with this s value */

    BOOLEAN dbg_roots;

    BOOLEAN dbg_t;
    BOOLEAN dbg_s;
    BOOLEAN dbg_elem;
    BOOLEAN dbg_vtk;

    BOOLEAN show_recursions;

    BOOLEAN show_directions;
    BOOLEAN show_intervals_w;
    BOOLEAN show_intervals_s;
    BOOLEAN show_intervals_t;

    BOOLEAN dbg_off;		/* This disables all dbg_ and show_ options */
} QI_OPTIONS;
extern QI_OPTIONS _phg_quad_interface_options;

/*-------------------------- Common functions --------------------------*/

typedef void (*FUNC_2D)(FLOAT x, FLOAT y, FLOAT *u);
#define FUNC_3D	DOF_USER_FUNC

int phgQuadInterfaceMarkElement(void *ls_func, int ls_order, void *ls_grad,
				ELEM_TYPE elem_type, void *elem_data);

typedef void (*QUAD_1D_FUNC)(FLOAT x, FLOAT *res, const FLOAT wgt, void *ctx);

void phgQuad1D(QUAD_1D_FUNC f, int dim, FLOAT a, FLOAT b, QUAD *quad,
		FLOAT *res, void *ctx);

int phgQuadInterfaceGetEdgeCuts(FUNC_3D ls_func, int ls_order, FUNC_3D ls_grad,
		const FLOAT v0[], const FLOAT v1[], FLOAT tol, FLOAT **cuts);

int phgUnivariatePolyRoots(int n, FLOAT a, FLOAT b, FLOAT tol,
			   FLOAT *buffer, int **mults);

/* RULE_HEADER is the header size (in FLOATs). Currently the header is set to:
 *		6 * npoints + (dim - 1) * 2 + nv_flag
 * where:
 *	npoints:	the number of points of the rule
 * 	dim:		1, 2 (triangle, rectangle) or 3 (tetrahedron, coboid),
 *	nv_flag:	0 (FALSE) or 1 (TRUE) (whether contains normal vectors)
 * the list of unit normal vectors are appended after the rule data.
 */

#define RULE_HEADER	1
#define SetRuleHeader(rule, npoints, dim, nv_flag) \
	(rule)[0] = _F(6.0) * (FLOAT)(npoints) + (FLOAT)( \
				2 * ((dim) - 1)) + ((nv_flag) == FALSE ? 0 : 1)

/*---------------------- Common functions ----------------------*/

#define PROJ		DOF_PROJ
#define PROJ_NONE	DOF_PROJ_NONE
#define PROJ_DOT	DOF_PROJ_DOT
#define PROJ_CROSS	DOF_PROJ_CROSS

FLOAT *phgQuadInterfaceRuleCreate(int dim, int np, const FLOAT *pts, int pinc,
		const FLOAT *wgts, int winc, FLOAT scale,
		const FLOAT *nv, int ninc);

int phgQuadInterfaceRuleInfo(FLOAT *rule, int *dim, FLOAT **pw, FLOAT **nv);

int phgQuadInterfaceRuleApply(void *func, int dim, PROJ projection,
			      FLOAT *rule, FLOAT *res);

int phgQuadInterfaceRuleMerge(FLOAT **rule, const FLOAT *rule0);

/*---------------------- PHG functions ----------------------*/

void phgQuadInterface(DOF *ls, DOF *ls_grad, ELEMENT *e, int quad_order,
		      FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);

INT phgQuadInterfaceMarkGrid(DOF *ls, DOF *ls_grad);
/* for backward compatibility */
#define phgQuadInterfaceMarkElements phgQuadInterfaceMarkGrid

/*---------------------- Tetrahedron/Triangle functions ----------------------*/

void phgQuadInterfaceTetrahedron(FUNC_3D ls, int ls_order, FUNC_3D ls_grad,
		FLOAT const tet[Dim + 1][Dim], int quad_order,
		FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);

void phgQuadInterfaceTriangle(FUNC_2D ls, int ls_order, FUNC_2D ls_grad,
		FLOAT const triangle[3][2], int quad_order,
		FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);

void phgQuadInterfaceTetrahedronFace(FUNC_3D ls, int ls_order, FUNC_3D ls_grad,
		FLOAT const tet[Dim + 1][Dim], int face, int quad_order,
		FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);

void phgQuadInterfaceFace(DOF *ls, DOF *ls_grad, ELEMENT *e, int face,
		int quad_order, FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);

/*------------------------ cuboid/rectangle functions -----------------------*/

void phgQuadInterfaceCuboid(FUNC_3D ls, int ls_order,
		FUNC_3D ls_grad, FLOAT const cuboid[2][3], int quad_order,
		FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);

void phgQuadInterfaceRectangle(FUNC_2D ls, int ls_order,
		FUNC_2D ls_grad, FLOAT const rect[2][2], int quad_order,
		FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);

void phgQuadInterfaceCuboidFace(FUNC_3D ls, int ls_order, FUNC_3D ls_grad,
		FLOAT const cuboid[2][Dim], int face,
		int quad_order, FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);

#ifdef __cplusplus
}
#endif

#define PHG_QUAD_INTERFACE_H
#endif
