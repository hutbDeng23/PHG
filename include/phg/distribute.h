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

/* $Id: distribute.h,v 1.7 2022/05/07 08:09:31 zlb Exp $ */

#ifndef PHG_DISTRIBUTE_H

typedef BOOLEAN (*PARTITIONER)(GRID *g, MPI_Comm newcomm,
					DOF *weights, FLOAT power);

#ifdef __cplusplus
extern "C" {
#endif

#if USE_MPI
int phgCreateComm_(int nprocs, MPI_Comm comm_world, MPI_Comm *comm);
void phgUpdateRNeighbours_(GRID *g);
void phgRebuildL2Gmap_(GRID *g, GTYPE type);
#endif
PARTITIONER phgGetPartitioner(int no);
void phgRedistributeGrid_(GRID *g, MPI_Comm comm_world, PARTITIONER part_func);
void phgRedistributeGrid(GRID *g);
#define phgPartitionGrid	phgRedistributeGrid
int phgBalanceGrid_(GRID *g, FLOAT lif_threshold, INT submesh_threshold,
	struct DOF_ *weights, FLOAT power, MPI_Comm comm_world,
	PARTITIONER part_func);
#define phgBalanceGrid(g, lif_threshold, submesh_threshold, weights, power) \
	phgBalanceGrid_(g, lif_threshold, submesh_threshold, weights, power, \
			phgComm, NULL)

#ifdef __cplusplus
}
#endif

#define PHG_DISTRIBUTE_H
#endif
