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

/* $Id: geom.h,v 1.6 2021/12/23 02:45:41 zlb Exp $ */

#ifndef PHG_GEOM_H

#ifdef __cplusplus
extern "C" {
#endif

void _phg_compute_elem_data(const GRID *g, const ELEMENT *e,
					FLOAT *vol, FLOAT *diam);

void phgGeomSaveNonLeafData_(GRID *g);
void phgGeomClearNonLeafData_(void);

void phgGeomCheck(GRID *g);
void  phgGeomInit_(GRID *g, BOOLEAN init_data);
#define phgGeomInit(g)	phgGeomInit_(g, TRUE)
FLOAT phgGeomGetVolume(GRID *g, ELEMENT *e);
FLOAT phgGeomGetDiameter(GRID *g, ELEMENT *e);
#define phgGeomGetJacobian(g, e)               \
    phgGeomGetJacobian_(g, e, NULL, NULL)    
const FLOAT *phgGeomGetJacobian_(GRID *g, SIMPLEX *e, const FLOAT *lambda, FLOAT *det);
#define phgGeomGetFaceArea(g, e, face_no)	\
    phgGeomGetFaceArea_(g, e, face_no, NULL, NULL)
FLOAT phgGeomGetFaceArea_(GRID *g, SIMPLEX *e, int face, const FLOAT *lambda, FLOAT *normal);
FLOAT phgGeomGetFaceAreaByIndex(GRID *g, INT face_no);
FLOAT phgGeomGetFaceDiameter(GRID *g, ELEMENT *e, int face);
const FLOAT *phgGeomGetFaceNormal(GRID *g, ELEMENT *e, int face);
void phgGeomGetFaceOutNormal(GRID *g, ELEMENT *e, int face, FLOAT *nv);
FLOAT (*phgGeomGetCorners(GRID *g, ELEMENT *e, FLOAT (*tet)[Dim]))[Dim];
FLOAT *phgGeomXYZ2Lambda(GRID *g, ELEMENT *e, FLOAT x, FLOAT y, FLOAT z,
			 FLOAT lambda[]);
void  phgGeomLambda2XYZ(GRID *g ,ELEMENT *e, const FLOAT lambda[],
				FLOAT *x, FLOAT *y, FLOAT *z);

#ifdef __cplusplus
}
#endif

#define PHG_GEOM_H
#endif
