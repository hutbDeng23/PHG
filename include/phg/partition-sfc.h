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

/* $Id: partition-sfc.h,v 1.1 2014/04/04 04:30:01 zlb Exp $ */
#ifndef PHG_PARTITION_SFC_H

/* float point number used by SFC */
typedef double SFC_FLOAT;
#define MPI_SFC_FLOAT MPI_DOUBLE
#define HSFC_EPSILON 1.6e-7

/* Struct for a 3-D inverse SFC */
typedef struct {
    SFC_FLOAT co[3];	/* coordinates (input, must be in the range [0,1]) */
    SFC_FLOAT sfc;  	/* key (output, in the range [0,1]) */
} SFC;

/* Struct for a 2-D inverse SFC */
typedef struct {
    SFC_FLOAT co[2];	/* coordinates (input, must be in the range [0,1]) */
    SFC_FLOAT sfc;  	/* key (output, in the range [0,1]) */
} SFC2;

#ifdef __cplusplus
extern "C" {
#endif

/* success: TRUE, fail: FALSE            */
/* n <= length of array hsfc             */
/* in general, n == length of array hsfc */
/* inverse Hilbert SFC                   */
BOOLEAN phgSFCInvHilbert3D(SFC *sfc, INT n);
BOOLEAN phgSFCInvHilbert2D(SFC2 *sfc, INT n);

/* inverse Morton SFC                   */
BOOLEAN phgSFCInvMorton3D(SFC *sfc, INT n);
BOOLEAN phgSFCInvMorton2D(SFC2 *sfc, INT n);

#if USE_MPI
/* Hilbert space-filling curve partitioner */
BOOLEAN phgPartitionSFC(GRID *g, MPI_Comm newcomm, DOF *weights, FLOAT power);
#endif

#ifdef __cplusplus
}
#endif

#define PHG_PARTITION_SFC_H
#endif
