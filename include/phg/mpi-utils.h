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

/* $Id: mpi-utils.h,v 1.2 2017/02/28 07:04:01 zlb Exp $ */

#ifndef PHG_MPI_UTILS_H

#ifndef USE_MPI
# define USE_MPI 1
#endif

#ifndef ALLTOALLV_HACK
# define ALLTOALLV_HACK 0
#endif

#if USE_MPI

#ifdef __cplusplus
extern "C" {
#endif

int MPI_Alltoallv___(
    const void *sbuff, const int *scnts, const int *sdsps, MPI_Datatype stype,
    void *rbuff, const int *rcnts, const int *rdsps, MPI_Datatype rtype,
    MPI_Comm comm, const char *file, int line);

int MPI_AlltoallvBegin(
    const void *sbuff, const int *scnts, const int *sdsps, MPI_Datatype stype,
    void *rbuff, const int *rcnts, const int *rdsps, MPI_Datatype rtype,
    MPI_Comm comm, MPI_Request **reqarray, int *nsends, int *nrecvs);
int MPI_AlltoallvEnd(int nsends, int nrecvs, MPI_Request **reqarray);

int MPI_AlltoallvInit(
    const void *sbuff, const int *scnts, const int *sdsps, MPI_Datatype stype,
    void *rbuff, const int *rcnts, const int *rdsps, MPI_Datatype rtype,
    MPI_Comm comm, MPI_Request **reqarray, int *nsends, int *nrecvs);

#if ALLTOALLV_HACK
# if !defined(DONT_RENAME_MPI_FUNCTIONS)
#  ifdef MPI_Alltoallv
#   undef MPI_Alltoallv
#  endif
/*---------------------------------------------------------------------------*/
#  define MPI_Alltoallv(sbuff,scnts,sdsps,stype, rbuff,rcnts,rdsps,rtype, comm)\
    MPI_Alltoallv___(sbuff,scnts,sdsps,stype, rbuff,rcnts,rdsps,rtype, comm, \
		     __FILE__, __LINE__)
/*---------------------------------------------------------------------------*/
# endif	/* !defined(DONT_RENAME_MPI_FUNCTIONS) */
#else	/* ALLTOALLV_HACK */
/*---------------------------------------------------------------------------*/
# define MPI_AlltoallvBegin(sbuff, scnts, sdsps, stype, rbuff, rcnts, rdsps, \
			   rtype, comm, reqs, nsends, nrecvs) \
    MPI_Alltoallv(sbuff, scnts, sdsps, stype, rbuff, rcnts, rdsps, rtype, comm)
# define MPI_AlltoallvEnd(nsends, nrecvs, reqs)
/*---------------------------------------------------------------------------*/
#endif	/* ALLTOALLV_HACK */

#ifdef __cplusplus
}
#endif

#endif	/* USE_MPI */

#define PHG_MPI_UTILS_H
#endif
