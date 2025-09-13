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

/* $Id: phg-mpi.h,v 1.1 2014/04/04 04:30:01 zlb Exp $ */

#ifndef PHG_MPI_H

#if USE_MPI

#define RegisterParallelFunction(f)	f(NULL)

#define ParallelFunction(f, flag)					\
    {									\
	static int func_no = -1;					\
	if (func_no < 0) {						\
	    func_no = phgRegisterParallelFunction(f, __func__,		\
			__FILE__, __LINE__);				\
	    return;							\
	}								\
	if (phgMasterSlave) {						\
	    if ((flag) && phgNProcs > 1 && phgRank == 0)		\
		phgCallSlaves(func_no);					\
	    if (!(flag) && phgRank > 0)					\
		return;							\
	}								\
    }

/*typedef GRID *(*PFunc_t)(GRID *);*/
typedef void (*PFunc_t)(GRID *);

extern MPI_Op MPI_MSM;
extern MPI_Datatype MPI_3DOUBLE;

#if FT_PHG == FT___FLOAT128
extern MPI_Op MPI_OP_SUM, MPI_OP_MAX;
extern MPI_Datatype MPI_FLOAT128;
#endif	/* FT_PHG == FT___FLOAT128 */

#ifdef __cplusplus
extern "C" {
#endif
int phgRegisterParallelFunction(PFunc_t func, const char *funcname,
		const char *file, int line);
MPI_Comm phgInitSetComm(MPI_Comm comm);
void phgInitMPI(int *argc, char ***argv);
void phgCallSlaves(int cmd);
void phgFinalizeMPI(void);
void phgSlave(void);

#ifdef __cplusplus
}
#endif

#define IS_RANK0		(phgRank == 0)

#else	/* USE_MPI */

typedef int MPI_Comm;
#define IS_RANK0		TRUE
#define MPI_COMM_WORLD		0
#define MPI_COMM_SELF		1
#define MPI_COMM_NULL		-1

#define ParallelFunction(f, flag)
#define RegisterParallelFunction(f)

#endif	/* USE_MPI */

extern int phgRank, phgNProcs, phgCmdSerialno;
extern MPI_Comm phgComm;
extern GRID *_phg_g_sys;

#define PHG_MPI_H
#endif
