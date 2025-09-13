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

/* $Id: xfem-utils.h,v 1.4 2022/07/29 08:37:49 zlb Exp $ */

#ifndef PHG_XFEM_UTILS_H	/* whole file */

#ifdef PHG_TO_P4EST
#define phgXFEMDot_			phgXFEMDot_p4est_
#define phgXFEMFaceJump2_		phgXFEMFaceJump2_p4est_
#endif	/* defined(PHG_TO_P4EST) */

FLOAT phgXFEMDot_(XFEM_INFO *xi, DOF *um, DOF *up, DOF *vm, DOF *vp,
		  FLOAT **e_sum, int added_order);
#define phgXFEMDot(xi, um, up, vm, vp) phgXFEMDot_(xi, um, up, vm, vp, NULL, 3)

FLOAT phgXFEMFaceJump2_(XFEM_INFO *xi, DOF *um, DOF *bm, DOF *up, DOF *bp,
			DOF *intf, PROJ proj, FLOAT h_power, int f_mask,
			FLOAT **e_sum, int added_order);
#define phgXFEMFaceJump2(xi, um, bm, up, bp, intf, proj, h_power, f_mask) \
	phgXFEMFaceJump2_(xi, um, bm, up, bp, intf, proj, h_power, f_mask, \
			  NULL, 3)
#define phgDofFaceJump2_(u, b, proj, h_power, f_mask, e_sum, added_order) \
	phgXFEMFaceJump2_(NULL, u, b, u, b, NULL, proj, h_power, f_mask, \
			  e_sum, added_order)
#define phgDofFaceJump2(u, b, proj, h_power, f_mask) \
	phgDofFaceJump2_(u, b, proj, h_power, f_mask, NULL, 0)

void phgXFEMProcessEmptyRows(SOLVER *solver);
#define process_empty_rows	phgXFEMProcessEmptyRows

FLOAT *phgOrthoMatrixGen(int n, int m, FLOAT *a);
void phgOrthoMatrixApply(int n, int bdim, int stride, int dim, FLOAT *buffer,
			 const FLOAT *P);

#define PHG_XFEM_UTILS_H
#endif	/* whole file */
