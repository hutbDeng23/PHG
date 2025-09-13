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

/* $Id: partition-utils.h,v 1.2 2020/05/11 08:54:41 zlb Exp $ */

#ifndef PHG_PARTITION_UTILS_H

#include "partition-sfc.h"

typedef struct {
    SFC_FLOAT key;/* key \in [0, 1]      */
    SFC_FLOAT w;  /* weight              */
    int pn;       /* partition number    */
} DOTS;


#ifdef __cplusplus
extern "C" {
#endif

/* success 0, or fail, serial version now */
/* used by phgImport                      */
/* Hilbert order */
int phgHilbertOrder(GRID *g);
/* Morton order */
int phgMortonOrder(GRID *g);

/* initial order */
/* 0: success, or fail */
int phgGridInitOrder(GRID *g, BOOLEAN dist);

#if USE_MPI
/* 1d partitioner, partition [0, 1] to np parts */
/* [0, p0), [p1, p2), ..., [p(np - 2), 1.]      */
/* such that weights of each part are equal     */
/* lenx, length of array x                      */
/* np, number of parts, should >= 1             */
/* return load balance efficiency fator         */
/*         ascending  order only                */
/*            version 1                         */
FLOAT phgPartition1DV1(DOTS *x, INT lenx, int np, const FLOAT speeds[],
			MPI_Comm comm);
FLOAT phgPartition1DV2(DOTS *x, INT lenx, int np, const FLOAT speeds[],
			MPI_Comm comm);

/* internal function used to register parameters */
void phgPartition1DOptionsSet();

/* remap partitions such that data migrations are minimized
 * three version. the second one uses less memories.
 * The first two are applied to grid directly, 
 * and results are stored in e->mark. These two can only be called just 
 * after phgPartition* are called...
 * The third one is general one without any limits.
 *
 * return value > 0, if partitions are changed;
 * return value = 0, if partitions arn't changed;
 * return value < 0, if subroutines fail and partitions keep unchanged*/

int phgPartitionRemap(GRID *g, MPI_Comm newcomm);
int phgPartitionRemap0(MPI_Comm comm, const double datasize[], int perm[],
		       MPI_Comm newcomm);

/* print quality information */
void phgPartitionQuality(GRID *g);

/* Partitioner random */
BOOLEAN phgPartitionerRandom(GRID *g, MPI_Comm newcomm, DOF *weights, FLOAT power);

/* Partitioner user */
BOOLEAN phgPartitionerUser(GRID *g, MPI_Comm newcomm, DOF *weights, FLOAT power);
#endif

#ifdef __cplusplus
}
#endif

#define PHG_PARTITION_UTILS_H
#endif	/* PHG_REMAP_H */
