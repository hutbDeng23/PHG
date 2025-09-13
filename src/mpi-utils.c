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

/* $Id: mpi-utils.c,v 1.63 2020/05/07 05:37:21 zlb Exp $ */

#ifdef __PHG__
# include "phg.h"
#else	/* defined(__PHG__) */
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>
# include <mpi.h>
# include "phg/mpi-utils.h"
# define USE_MPI 1
# define phgAlloc(size)			malloc(size)
# define phgRealloc_(p,size,size0)	realloc(p,size)
# define phgFree(p)			if (p != NULL) free(p)
# define phgPrintf			printf
# define phgGetTime(tarray)		0.0
static int phgVerbosity = 0;
#endif	/* defined(__PHG__) */

#include <string.h>	/* memcpy() */

#if USE_MPI	/* whole remaining part of the file */

/*---------------------- Workaround for Sparse Alltoallv --------------------*/

#if ALLTOALLV_HACK

/* Limit on the size of a single message (0 means no limit)
 *
 * Note: This is a workaround for some MPI implementations which don't allow
 * very large single messages, such as MVAPICH2 on LSSC3. */
#define MPI_MESSAGE_LIMIT ((size_t)(1<<30))	/* message size limit (bytes) */

#define TAG 314

#if !defined(MPI_VERSION) || MPI_VERSION < 3
  typedef int (MPIAPI isend_t)(void *, int, MPI_Datatype, int, int,
		MPI_Comm, MPI_Request *);
#else	/* MPI_VERSION < 3 */
  typedef int (MPIAPI isend_t)(const void *, int, MPI_Datatype, int, int,
		MPI_Comm, MPI_Request *);
#endif	/* MPI_VERSION < 3 */
typedef int (MPIAPI irecv_t)(void *, int, MPI_Datatype, int, int, MPI_Comm,
		MPI_Request *);

static isend_t *isend = MPI_Isend;
static irecv_t *irecv = MPI_Irecv;

#ifdef MPI_AlltoallvBegin
# undef MPI_AlltoallvBegin
#endif
int
MPI_AlltoallvBegin(
    const void *sendbuf, const int *sendcnts, const int *sdispls,
	MPI_Datatype sendtype,
    void *recvbuf, const int *recvcnts, const int *rdispls,
	MPI_Datatype recvtype,
    MPI_Comm comm, MPI_Request **reqarray, int *nsends, int *nrecvs)
/* Customized Alltoallv function which does not send/recev 0 length blocks.
   This function is derived from intra_Alltoallv() of MPICH 1.2.6 */
{
    int		nprocs, i, j, cnt, alloc, peer, rank;
    MPI_Aint	lb, send_extent, recv_extent;
    int		mpi_errno = MPI_SUCCESS;
    size_t	limit, remain, size;
    char	*p;

    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);
#if 1	/* MPI-2 */
    MPI_Type_get_extent(sendtype, &lb, &send_extent);
    MPI_Type_get_extent(recvtype, &lb, &recv_extent);
#else	/* MPI-1 */
    MPI_Type_extent(sendtype, &send_extent);
    MPI_Type_extent(recvtype, &recv_extent);
#endif

    assert(MPI_MESSAGE_LIMIT == 0 || recv_extent == send_extent);

    *reqarray = phgAlloc((alloc = 2 * nprocs) * sizeof(**reqarray));
    cnt = 0;

    limit = MPI_MESSAGE_LIMIT / recv_extent;
    for (i = 0; i < nprocs; i++) {
	peer = (rank + i) % nprocs;
	if (recvcnts[peer] <= 0)
	    continue;
	p = (char *)recvbuf + rdispls[peer] * recv_extent;
	remain = recvcnts[peer];
	while (remain > 0) {
	    size = (remain > limit && limit > 0 ? limit : remain);
	    remain -= size;
	    if (cnt >= alloc)
		*reqarray = phgRealloc_(*reqarray,
					(alloc += 16) * sizeof(**reqarray),
					cnt * sizeof(**reqarray));
	    if ((mpi_errno = irecv(p, size, recvtype, peer, TAG, comm,
				   &(*reqarray)[cnt])) != MPI_SUCCESS)
		break;
	    p += size * recv_extent;
	    cnt++;
	}
	if (mpi_errno != MPI_SUCCESS)
	   break;
    }
    if (nrecvs != NULL)
	*nrecvs = cnt;

    if (mpi_errno == MPI_SUCCESS) {
	limit = MPI_MESSAGE_LIMIT / send_extent;
	for (i = 0; i < nprocs; i++) {
	    peer = (rank + i) % nprocs;
	    if (sendcnts[peer] <= 0)
		continue;
	    p = (char *)sendbuf + sdispls[peer] * send_extent;
	    remain = sendcnts[peer];
	    while (remain > 0) {
		size = (remain > limit && limit > 0 ? limit : remain);
		remain -= size;
		if (cnt >= alloc)
		    *reqarray = phgRealloc_(*reqarray,
					    (alloc += 16) * sizeof(**reqarray),
					    cnt * sizeof(**reqarray));
		if ((mpi_errno = isend(p, size, sendtype, peer, TAG, comm,
				       &(*reqarray)[cnt])) != MPI_SUCCESS)
		    break;
		cnt++;
		p += size * send_extent;
	    }
	    if (mpi_errno != MPI_SUCCESS)
		break;
	}
    }
    if (nsends != NULL)
	*nsends = cnt - *nrecvs;

    if (mpi_errno) {
	for (j = 0; j < cnt; j++)
	    MPI_Cancel(&(*reqarray)[j]);
	phgFree(*reqarray);
	*reqarray = NULL;
	if (nsends != NULL)
	    *nsends = 0;
	if (nrecvs != NULL)
	    *nrecvs = 0;
    }

    return mpi_errno;
}

#ifdef MPI_AlltoallvEnd
# undef MPI_AlltoallvEnd
#endif
int
MPI_AlltoallvEnd(int nsends, int nrecvs, MPI_Request **reqarray)
/* Customized Alltoallv function which does not send/recev 0 length blocks.
   This function is derived from intra_Alltoallv() of MPICH 1.2.6 */
{
    int		j, rcnt;
    int		mpi_errno = MPI_SUCCESS;
    MPI_Status	*starray;

    if (*reqarray == NULL || (rcnt = nsends + nrecvs) <= 0) {
	phgFree(*reqarray);
	*reqarray = NULL;
	return mpi_errno;
    }

    starray = phgAlloc(rcnt * sizeof(MPI_Status));
    mpi_errno = MPI_Waitall(rcnt, *reqarray, starray);
    if (mpi_errno == MPI_ERR_IN_STATUS) {
	for (j = 0; j < rcnt; j++) {
	    if (starray[j].MPI_ERROR != MPI_SUCCESS)
		    mpi_errno = starray[j].MPI_ERROR;
	}
    }
    phgFree(starray);
    phgFree(*reqarray);
    *reqarray = NULL;

    return mpi_errno;
}

int
MPI_Alltoallv___(
    const void *sendbuf, const int *sendcnts, const int *sdispls,
	MPI_Datatype sendtype,
    void *recvbuf, const int *recvcnts, const int *rdispls,
	MPI_Datatype recvtype,
    MPI_Comm comm, const char *file, int line)
{
    MPI_Request *req;
    int mpi_errno, nprocs, nsends, nrecvs;
#ifdef PHG_PHG_H
    extern void phgInfo(int, const char *, ...);
#endif
    MPI_Aint lb, send_extent, recv_extent;
    size_t limit;
    int *scnts, *sdsps, *rcnts, *rdsps, *scnts0, *rcnts0;
    int i, flag;

    MPI_Comm_size(comm, &nprocs);
#if 1	/* MPI-2 */
    MPI_Type_get_extent(sendtype, &lb, &send_extent);
    MPI_Type_get_extent(recvtype, &lb, &recv_extent);
#else	/* MPI-1 */
    MPI_Type_extent(sendtype, &send_extent);
    MPI_Type_extent(recvtype, &recv_extent);
#endif

    assert(MPI_MESSAGE_LIMIT == 0 || recv_extent == send_extent);
#if 0
    /* enforce limit on the message sizes in this function */
    limit = MPI_MESSAGE_LIMIT / recv_extent;
#else
    /* don't enforce limit on the message sizes in this function
     * (do it in MPI_AlltoallBegin instead) */
    limit = 0;
#endif

    if (limit > 0) {
	scnts = phgAlloc(6 * nprocs * sizeof(*scnts));
	sdsps = scnts + 1 * nprocs;
	rcnts = scnts + 2 * nprocs;
	rdsps = scnts + 3 * nprocs;
	scnts0 = scnts + 4 * nprocs;
	rcnts0 = scnts + 5 * nprocs;
	memcpy(scnts0, sendcnts, nprocs * sizeof(*scnts));
	memcpy(rcnts0, recvcnts, nprocs * sizeof(*rcnts));
	memcpy(sdsps, sdispls, nprocs * sizeof(*sdsps));
	memcpy(rdsps, rdispls, nprocs * sizeof(*rdsps));
    }
    else {
	scnts = (void *)sendcnts;
	sdsps = (void *)sdispls;
	rcnts = (void *)recvcnts;
	rdsps = (void *)rdispls;
    }

    while (TRUE) {
	flag = 0;
	if (limit > 0){
	    for (i = 0; i < nprocs; i++) {
		if ((scnts0[i] -= (scnts[i] =
			(scnts0[i] > limit ? limit : scnts0[i]))) != 0)
		    flag = 1;
		if ((rcnts0[i] -= (rcnts[i] =
			(rcnts0[i] > limit ? limit : rcnts0[i]))) != 0)
		    flag = 1;
	    }
	}
	if ((mpi_errno = MPI_AlltoallvBegin(sendbuf, scnts, sdsps, sendtype,
				recvbuf, rcnts, rdsps, recvtype,
				comm, &req, &nsends, &nrecvs)) != MPI_SUCCESS)
	    break;
#if 0
	phgInfo(-1, "%s: %s:%d, scnt = %d, rcnt = %d, errno = %d\n", __func__,
		file, line, nsends, nrecvs, mpi_errno);
#endif
	if ((mpi_errno = MPI_AlltoallvEnd(nsends, nrecvs, &req)) != MPI_SUCCESS)
	    break;
	if (limit <= 0)
	    break;
	i = flag;
	MPI_Allreduce(&i, &flag, 1, MPI_INT, MPI_MAX, comm);
	if (!flag)
	    break;
	for (i = 0; i < nprocs; i++) {
	    sdsps[i] += scnts[i];
	    rdsps[i] += rcnts[i];
	}
    }

    if (limit > 0)
	phgFree(scnts);

    return mpi_errno;
}

int
MPI_AlltoallvInit(
    const void *sendbuf, const int *sendcnts, const int *sdispls,
	MPI_Datatype sendtype,
    void *recvbuf, const int *recvcnts, const int *rdispls,
	MPI_Datatype recvtype,
    MPI_Comm comm, MPI_Request **reqarray, int *nsends, int *nrecvs)
{
    int mpi_errno;

    irecv = MPI_Recv_init;
    isend = MPI_Send_init;
    mpi_errno = MPI_AlltoallvBegin(sendbuf, sendcnts, sdispls, sendtype,
				   recvbuf, recvcnts, rdispls, recvtype,
				   comm, reqarray, nsends, nrecvs);
    irecv = MPI_Irecv;
    isend = MPI_Isend;

    return mpi_errno;
}

#endif	/* ALLTOALLV_HACK */

#endif	/* USE_MPI */
