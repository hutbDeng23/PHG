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

/* $Id: quad.h,v 1.14 2021/05/08 13:00:14 zlb Exp $ */

#ifndef PHG_QUAD_H

typedef struct QUAD_ {
    char *name;			/* name of the quadrature formulae */
    int dim;			/* dimension, 1: edge, 2: face, 3: tetra */
    int order;			/* exact for polynomials of order 'order' */
    int npoints;		/* number of points */
    FLOAT *points;		/* barycentric coordinates of quad. points */
    FLOAT *weights;		/* weights */
    SHORT id;			/* id (for use with reference count) */
} QUAD;

/* cached data for a quadrature rule in an element or on a face */
typedef struct {
    ELEMENT	*e;		/* cache at element 'e' */
    FLOAT	*data;		/* cached values for 'e' */
    /* For DG: need to cache at most 2 (ELEMENT,face) pairs for a face */
    ELEMENT	*e1;		/* the peer element (neighbour) on the face */
    FLOAT	*data1;		/* cached values for "e1" */
    char	face, face1;	/* which face of the element (\in {0,1,2,3}) */
} QUAD_CACHE;

/* list of cahed data, they are typically stored in DOF_TYPE.cache_* members.
 * The list is of length 'n', indexed by QUAD.id, with n - 1 being the largest
 * id of the QUAD referenced by the DOF_TYPE */
typedef struct {
    QUAD_CACHE **caches;
    SHORT n;
} QUAD_CACHE_LIST;

#define QUAD_DEFAULT	-1

#ifdef __cplusplus
extern "C" {
#endif

    void phgQuadFree(QUAD **quad);

    void phgQuadReset(void);
    void phgQuadClearDofCache(void **clist, QUAD *quad, BOOLEAN final);

    QUAD *phgQuadGetQuad1D(int order);
    QUAD *phgQuadGetQuad2D(int order);
    QUAD *phgQuadGetQuad3D(int order);

    FLOAT *phgQuadGetRule_(GRID *g, ELEMENT *e, int which, int order);
#define phgQuadGetRule3D(g, e, order) \
	phgQuadGetRule_(g, e, 3 - 1, order)
#define phgQuadGetRule2D(g, e, face_no, order) \
	phgQuadGetRule_(g, e, (face_no) * 3 + 2 - 1, order)
#define phgQuadGetRule1D(g, e, edge_no, order) \
	phgQuadGetRule_(g, e, (edge_no) * 3 + 1 - 1, order)

    /* functions for caching basis functions and derivativesi or user functions
     * at quadrature points */
    const FLOAT *phgQuadGetFuncValues(GRID *g, ELEMENT *e, int dim,
				      DOF_USER_FUNC userfunc, QUAD *quad);
    const FLOAT *phgQuadGetBasisValues(ELEMENT *e, DOF *u, int n, QUAD *quad);
    const FLOAT *phgQuadGetBasisGradient(ELEMENT *e, DOF *u, int n, QUAD *quad);
    const FLOAT *phgQuadGetBasisCurl(ELEMENT *e, DOF *u, int n, QUAD *quad);
    const FLOAT *phgQuadGetDofValues(ELEMENT *e, DOF *u, QUAD *quad);

    /*-------------------- 2D functions for DG --------------------*/
    const FLOAT *phgQuadFaceGetBasisValues(ELEMENT *e, int face,
					   DOF *u, int n, QUAD *quad);
    const FLOAT *phgQuadFaceGetBasisGradient(ELEMENT *e, int face,
					     DOF *u, int n, QUAD *quad);
    const FLOAT *phgQuadFaceGetDofValues(ELEMENT *e, int face,
					 DOF *u, QUAD *quad);

    FLOAT phgQuadFaceDofABas(ELEMENT * e, int face, DOF *u, DOF *A,
			     DOF *v, int m, int order);
    FLOAT phgQuadFaceDofADnBas(ELEMENT *e, int face, DOF *u, DOF *A,
			       DOF *v, int m, int order);
    FLOAT phgQuadFaceABasABas(DOF *A_u, ELEMENT *e_u, int face_u, DOF *u, int n,
			      DOF *A_v, ELEMENT *e_v, int face_v, DOF *v, int m,
			      int order);
    FLOAT phgQuadFaceABasADnBas(
			DOF *A_u, ELEMENT *e_u, int face_u, DOF *u, int n,
			DOF *A_v, ELEMENT *e_v, int face_v, DOF *v, int m,
			int order);
    FLOAT phgQuadFaceADnBasADnBas(
			DOF *A_u, ELEMENT *e_u, int face_u, DOF *u, int n,
			DOF *A_v, ELEMENT *e_v, int face_v, DOF *v, int m,
			int order);
    /*-------------------- old 2D functions --------------------*/
    FLOAT phgQuadFaceDofDotBas(ELEMENT *e, int face, DOF *u, DOF_PROJ proj,
				DOF *v, int n, int order);
    FLOAT phgQuadFaceBasDotBas(ELEMENT *e, int face, DOF *u, int n,
				DOF *v, int m, int order);
    FLOAT phgQuadFaceBasABas(ELEMENT * e, int face, DOF *u, int n, DOF *A,
				DOF *v, int m, int order);
    FLOAT phgQuadFaceDofDotDof(ELEMENT *e, int face, DOF *u, DOF_PROJ proj,
				DOF *v, int order);
    FLOAT  *phgQuadFaceADofCrossDof(ELEMENT *e, int face, DOF *A, 
				DOF *u, DOF_PROJ u_proj,
				DOF *v, DOF_PROJ v_proj, 
				int order, FLOAT *reval);
    DOF *phgQuadFaceJumpN(DOF *u, DOF_PROJ proj, const char *name, int order,
				DOF *gn);
#define phgQuadFaceJump(u, proj, name, order) \
	phgQuadFaceJumpN(u, proj, name, order, NULL)

    /*-------------------- 3D functions --------------------*/
    FLOAT phgQuadDofNormP(ELEMENT *e, DOF *u, int order, int p);

    FLOAT phgQuadDofDotDof(ELEMENT *e, DOF *u, DOF *v, int order);
/*
FLOAT phgQuadGradBasDotGradBas(ELEMENT *e, DOF *u, int n, DOF *v, int m,
				int order);
				*/
#define phgQuadGradBasDotGradBas(e, u, n, v, m, order) \
	phgQuadGradBasAGradBas(e, u, n, NULL, v, m, order)
    FLOAT phgQuadGradBasAGradBas(ELEMENT *e, DOF *u, int n, DOF *A, DOF *v,
				 int m, int order);
    FLOAT *phgQuadDofTimesBas(ELEMENT *e, DOF *u, DOF *v, int n, int order,
			      FLOAT *res);
/*
 * FLOAT phgQuadBasDotBas(ELEMENT *e, DOF *u, int n, DOF *v, int m, int order); */
#define  phgQuadBasDotBas(e, u, n, v, m, order) \
	phgQuadBasABas(e, u, n, NULL, v, m, order)
    FLOAT phgQuadBasABas(ELEMENT *e, DOF *u, int n, DOF *A, DOF *v, int m,
			 int order);
/*FLOAT phgQuadCurlBasDotCurlBas(ELEMENT *e, DOF *u, int n, DOF *v, int m,
				 int order); */
#define phgQuadCurlBasDotCurlBas(e, u, n, v, m, order)	\
	phgQuadCurlBasACurlBas(e, u, n, NULL, v, m, order)
    FLOAT phgQuadCurlBasACurlBas(ELEMENT *e, DOF *u, int n, DOF *A, DOF *v,
				 int m, int order);
    FLOAT phgQuadBasACurlBas(ELEMENT *e, DOF *u, int n, DOF *A, DOF *v, int m,
		       int order);
    FLOAT phgQuadDofDotBas(ELEMENT *e, DOF *u, DOF *v, int n, int order);
    FLOAT phgQuadDofABas(ELEMENT *e, DOF *u, DOF *A, DOF *v, int m, int order);
    FLOAT phgQuadFuncDotBas(ELEMENT *e, DOF_USER_FUNC userfunc, DOF *u, int n,
			    int order);
#define phgQuadDofDotGradBas(e, u, v, m, order) \
	phgQuadDofAGradBas(e, u, NULL, v, m, order)
    FLOAT phgQuadDofAGradBas(ELEMENT *e, DOF *u, DOF *A, DOF *v, int m,
			     int order);
#define phgQuadDofDotCurlBas(e, u, v, m, order) \
	phgQuadDofACurlBas(e, u, NULL, v, m, order)
    FLOAT phgQuadDofACurlBas(ELEMENT *e, DOF *u, DOF *A, DOF *v, int m,
			     int order);
    FLOAT *phgQuadDofDotCurlBas_(ELEMENT *e, DOF *u, DOF *v, int m,
			     int order, FLOAT *res);
    FLOAT phgQuadGradBasDotBas(ELEMENT *e, DOF *s, int m, DOF *v, int n,
			       int order);
#define phgQuadDivBasDotDivBas(e, u, n, v, m, order) \
	phgQuadDivBasADivBas(e, u, n, NULL, v, m, order)
    FLOAT phgQuadDivBasADivBas(ELEMENT *e, DOF *u, int n, DOF *A, DOF *v,
				 int m, int order);
    FLOAT phgQuadLapBasDotLapBas(ELEMENT *e, DOF *u, int n, DOF *v, int m,
				 int order);

#ifdef __cplusplus
}
#endif
#define PHG_QUAD_H
#endif
