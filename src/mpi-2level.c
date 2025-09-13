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

/* $Id: mpi-2level.c,v 1.30 2020/05/07 05:37:21 zlb Exp $ */

/* This file implements internode-intranode 2-level collective communications.
 *
 * Note:
 *
 * . The files mpi-2level.c and mpi-2level.h can be used independently of PHG
 *   as a standalone library, with the exported functions and variables being
 *   defined in mpi-2level.h.
 *
 * . See test/mpi_test.c for an example on using these functions.
 *
 * . Steps to use the 2-level collective communication functions in an
 *   application program:
 *
 *   1. Include "mpi-2level.h" (after "mpi.h"), which redefines the standard
 *	collective functions (e.g., MPI_Alltoallv) as the corresponding 2-level
 *	functions (e.g., MPI_Alltoallv2).
 *
 *   2. Call MPI_L2setup with appropriate parameters after MPI_Init or after
 *      creation of communicators.
 *
 *	The 2-level hierarchy is attached to a communicator as an user-defined
 *	attribute, which is only effective for that communicator, and will be
 *	automatically released when the communicator is destroyed.
 *
 *   Below is an example of doing the step 1 above through 'hijacking' "mpi.h",
 *   without needing to actually modify the application code:
 *
 *   a) Use the following bash commands to create a dummy MPI header file
 *	(which includes the original "mpi.h" followed by "mpi-2level.h"):
 *---------------------------------- cut --------------------------------------
		real_mpi_h="/usr/include/mpich-i386/mpi.h"	# change me!
		mpi_2level_incdir="$HOME/phg/include/phg"	# change me!
		cat <<END >mpi.h
#ifndef MPI_2LEVEL_MPI_H
# define MPI_2LEVEL_MPI_H
# include "$real_mpi_h"
# include "$mpi_2level_incdir/mpi-2level.h"
# pragma message("*** Notice: 2-level collective communications will be used.")
#endif
END
 *---------------------------------- cut --------------------------------------
 *
 *   b) Modify CPPFLAGS, CFLAGS or CXXFLAGS (adding an appropriate "-I" flag)
 *	to make the compiler find the dummy "mpi.h" instead of the original
 *	"mpi.h".
 */

/* The following macro prevents redefining MPI functions in this file */
#define DONT_RENAME_MPI_FUNCTIONS

#ifdef __PHG__
# include "phg.h"
#else	/* defined(__PHG__) */
# include <stdio.h>
# include <assert.h>
# include <mpi.h>
# include "phg/mpi-2level.h"
# define USE_MPI 1
# define phgAlloc			malloc
# define phgCalloc			calloc
# define phgRealloc_(p,size,size0)	realloc(p,size)
# define phgFree(p)			if (p != NULL) free(p)
# define phgPrintf			if (phgRank==0) printf
# define phgGetTime(tarray)		0.0
int phgVerbosity = 0;
int phgRank, phgNProcs;
MPI_Comm phgComm;
#endif	/* defined(__PHG__) */

#include <stdlib.h>	/* qsort() */
#include <string.h>	/* memcpy() */
#include <limits.h>	/* INT_MAX */

#if USE_MPI	/* whole remaining part of the file */

#ifndef MPIAPI
# define MPIAPI
#endif

/*--------------- Two Level Collective Communication Functions --------------*/

/* The functions below implement inter-node/intra-node two level collective
 * communications, whose usage is controlled by the option '-mpi_2level' . */

/* cast 'const xxx *' pointers to 'void *' for MPI <= 2.0 */
#ifdef CNSTV
# undef CNSTV
#endif
#if defined(MPI_VERSION) && MPI_VERSION >= 3
# define CNSTV			/* do nothing */
#else	/* MPI_VERSION >= 3 */
# define CNSTV	(void *)	/* cast pointers to '(void *)' */
#endif	/* MPI_VERSION >= 3 */

/* Note: the following lines may be useful in tracking in which function
 * the program has stopped in case of MPI error. */
#if 0	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
static int
dummy_(
void *sendbuf, int sendcount, MPI_Datatype sendtype,
void *recvbuf, int recvcount, MPI_Datatype recvtype, int root,
MPI_Comm comm,
    const char *file, int line)
{
    int ret;
    phgInfo(-1, "MPI_Gather (%s:%d).\n", file, line);
    ret = MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
    phgInfo(-1, "MPI_Gather (%s:%d) done.\n", file, line);
    return ret;
}
#define MPI_Gather(a,b,c,d,e,f,g,h) dummy_(a,b,c,d,e,f,g,h,__FILE__,__LINE__)
#endif 	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/* global (exported) variables */

/* MPI_L2coll_flag == 0 disables 2-level collective communications */
#ifdef PHG_PHG_H
BOOLEAN MPI_L2coll_flag = 0;
#else	/* PHG_PHG_H */
int MPI_L2coll_flag = 0;
#endif	/* PHG_PHG_H */

/* Internal variables. */

static int l2key = MPI_KEYVAL_INVALID; 	/* attr. key for 2L comm. */

static int MPIAPI
l2destructor(MPI_Comm comm, int keyval, void *data, void *extra)
{
    L2DATA *l2 = data;
    assert(keyval != MPI_KEYVAL_INVALID && keyval == l2key);
    l2->ref_count -= 1;
    if (l2->ref_count == 0) {
	if (phgVerbosity > 0)
	    phgPrintf("##### Destroy attr for 2-level communications. #####\n");
	if (l2->Lcomm != MPI_COMM_NULL)
	    MPI_Comm_free(&l2->Lcomm);
	if (l2->Qcomm != MPI_COMM_NULL)
	    MPI_Comm_free(&l2->Qcomm);
#if 1	/* MPI-2 */
	MPI_Comm_free_keyval(&l2key);
#else	/* MPI-1 */
	MPI_Keyval_free(&l2key);
#endif
	assert(l2key == MPI_KEYVAL_INVALID);
	phgFree(l2->D2Qmap);
	phgFree(l2->L2Gmap);
	phgFree(l2->G2Dmap);
	phgFree(l2->Q2Dmap);
	phgFree(l2);
    }

    return MPI_SUCCESS;
}

static int MPIAPI
l2copier(MPI_Comm comm, int keyval, void *extra, void *in, void *out, int *flag)
{
    L2DATA *l2 = in;

    assert(keyval == l2key && l2 != NULL);
    if (phgVerbosity > 0)
	phgPrintf("##### Copy attr for 2-level communications. #####\n");
    l2->ref_count += 1;
    *(void **)out = in;
    *flag = 1;

    return MPI_SUCCESS;
}

static char (*processor_names)[MPI_MAX_PROCESSOR_NAME] = NULL;

static int
name_comp(const void *p1, const void *p2)
{
    return strcmp(processor_names[*(const int *)p1],
		  processor_names[*(const int *)p2]);
}

void
MPI_L2setup(MPI_Comm *comm, int reorder_type, int ppn)
/* Sets up attribute for 2-level communications. Parameters:
 *   'reorder_type'
 *   	controls whether or how to reorder processes.
 *   'ppn'
 *	>= 2: local comms consist of procs with same rank%ppn value
 *	 < 2: local comms consist of procs with same processor name
 */
{
    int i, rank, nprocs, tmp[2], *color, *index;
    L2DATA *l2;

    MPI_Comm_rank(*comm, &rank);
    MPI_Comm_size(*comm, &nprocs);

#ifndef PHG_PHG_H
    phgRank = rank;
    phgNProcs = nprocs;
    phgComm = *comm;
#endif	/* defined(PHG_PHG_H) */

    if (phgVerbosity > 0)
	phgPrintf("##### %s: setting up two-level topology.\n", __func__);

    color = phgAlloc((2 * nprocs + 1) * sizeof(*color));

    if (ppn > 0) {
	/* group processes by 'ppn' trunks */
	if (rank == 0) {
	    int c = -1, block = 0, size = 0, nblocks = (nprocs + ppn - 1) / ppn;
	    for (i = 0; i < nprocs; i++) {
		if (i >= block + size) {
		    c++;
		    block += size;
		    size = nprocs / nblocks + (c < nprocs % nblocks ? 1 : 0);
		}
#if 0
		color[i] = i / ppn;
#else
		color[i] = c;
#endif
	    }
	    /* set color[nprocs] = max(color) */
	    color[i] = nblocks - 1;
	}
    }
    else {
	/* group processes by MPI_Get_processor_name() */
	char name[MPI_MAX_PROCESSOR_NAME];

	/* Collect processor names */
	MPI_Get_processor_name(name, &i);
	name[i] = '\0';

	if (rank == 0)
	    processor_names = phgAlloc(nprocs * MPI_MAX_PROCESSOR_NAME);
	MPI_Gather(name, MPI_MAX_PROCESSOR_NAME, MPI_BYTE,
	       processor_names, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, 0, *comm);
	if (rank == 0) {
	    index = color + 1 + nprocs;
	    for (i = 0; i < nprocs; i++)
		index[i] = i;
	    qsort(index, nprocs, sizeof(*index), name_comp);
	    color[index[0]] = 0;
	    for (i = 1; i < nprocs; i++) {
		if (!strcmp(processor_names[index[i]],
			    processor_names[index[i - 1]]))
		    color[index[i]] = color[index[i - 1]];
		else
		    color[index[i]] = color[index[i - 1]] + 1;
	    }
	    phgFree(processor_names);
	    processor_names = NULL;
	    color[nprocs] = color[index[nprocs - 1]];	/* max(color) */
	}
    }
    MPI_Bcast(color, nprocs + 1, MPI_INT, 0, *comm);

    if (color[nprocs] == 0 || color[nprocs] == nprocs - 1) {
	/* number of nodes == 1 or number of processes per node == 1,
	 * don't set attributes. */
	if (phgVerbosity > 0)
	    phgPrintf("##### %s: 1 node or 1 proc/node, do nothing. #####\n",
			__func__);
	phgFree(color);
	return;
    }

    if (l2key == MPI_KEYVAL_INVALID &&
#if 1	/* MPI-2 */
	MPI_Comm_create_keyval(l2copier, l2destructor, &l2key, NULL)
#else	/* MPI-1 */
	MPI_Keyval_create(l2copier, l2destructor, &l2key, NULL)
#endif
	!= MPI_SUCCESS)
	MPI_Abort(*comm, 99);

    if (phgVerbosity > 0)
	phgPrintf("##### Create attr for 2-level communications. #####\n");

    l2 = phgCalloc(1, sizeof(*l2));

    /* Create Lcomm using values in color[] */
    MPI_Comm_split(*comm, color[rank], rank, &l2->Lcomm);
    assert(l2->Lcomm != MPI_COMM_NULL);
    MPI_Comm_rank(l2->Lcomm, &l2->Lrank);
    MPI_Comm_size(l2->Lcomm, &l2->Lsize);
    phgFree(color);

    /* Create Qcomm which contains the masters of Lcomm */
    i = (l2->Lrank == 0 ? 0 : MPI_UNDEFINED);
    MPI_Comm_split(*comm, i, rank, &l2->Qcomm);
    if (l2->Qcomm == MPI_COMM_NULL) {
	l2->Qrank = MPI_PROC_NULL;
	l2->Qsize = 0;
    }
    else {
	MPI_Comm_rank(l2->Qcomm, &l2->Qrank);
	MPI_Comm_size(l2->Qcomm, &l2->Qsize);
    }
    tmp[0] = l2->Qsize;
    tmp[1] = l2->Qrank;
    MPI_Bcast(tmp, 2, MPI_INT, 0, l2->Lcomm);
    l2->Qsize = tmp[0];
    l2->Qrank = tmp[1];

   if (reorder_type != L2_NONE) {
	/* reorder process in *comm according to the two level topology */
	MPI_Comm newcomm = MPI_COMM_NULL;
	switch (reorder_type) {
	    case L2_ROUND_ROBIN:
		if (phgVerbosity > 0)
		    phgPrintf("##### Process ordering: round-robin. #####\n");
		MPI_Comm_split(*comm, 1, l2->Lrank, &newcomm);
		break;
	    case L2_BLOCK:
		if (phgVerbosity > 0)
		    phgPrintf("##### Process ordering: block-wise. #####\n");
#if 0
		/* FIXME: on LSSC3 with "-np 6 -n 2" MVAPICH2 gives wrong
		 * ordering in Lcomm as: (0,1,2), (4,3,5) */
		MPI_Comm_split(*comm, 1, l2->Qrank, &newcomm);
#else
		if (l2->Lrank == 0) {
		    MPI_Scan(&l2->Lsize, &i, 1, MPI_INT, MPI_SUM,
				l2->Qcomm);
		    i -= l2->Lsize;
		}
		MPI_Bcast(&i, 1, MPI_INT, 0, l2->Lcomm);
		MPI_Comm_split(*comm, 1, i + l2->Lrank, &newcomm);
#endif
		break;
	}
	MPI_Comm_free(comm);
	*comm = newcomm;
	MPI_Comm_rank(newcomm, &rank);
	MPI_Comm_size(newcomm, &nprocs);
    }

    /* Gather rank in Lcomm */
    l2->L2Gmap = phgAlloc(l2->Lsize * sizeof(*l2->L2Gmap));
    MPI_Allgather(&rank, 1, MPI_INT, l2->L2Gmap, 1, MPI_INT, l2->Lcomm);

    /* Gather Qrank */
    l2->D2Qmap = phgAlloc(nprocs * sizeof(*l2->D2Qmap));
    MPI_Allgather(&l2->Qrank, 1, MPI_INT, l2->D2Qmap, 1, MPI_INT, *comm);

    if (reorder_type == L2_BLOCK) {
	l2->G2Dmap = NULL;
    }
    else {
	l2->G2Dmap = phgAlloc(nprocs * sizeof(*l2->G2Dmap));
	if (l2->Lrank == 0) {
	    MPI_Scan(&l2->Lsize, &i, 1, MPI_INT, MPI_SUM, l2->Qcomm);
	    i -= l2->Lsize;
	}
	MPI_Bcast(&i, 1, MPI_INT, 0, l2->Lcomm);
	i += l2->Lrank;
	MPI_Allgather(&i, 1, MPI_INT, l2->G2Dmap, 1, MPI_INT, *comm);
	/* Update D2Qmap[] */
	index = phgAlloc(nprocs * sizeof(*index));
	for (i = 0; i < nprocs; i++)
	    index[l2->G2Dmap[i]] = l2->D2Qmap[i];
	phgFree(l2->D2Qmap);
	l2->D2Qmap = index;
    }

    /* build Q2Dmap[] array (only needed on local masters) */
    if (l2->Lrank != 0) {
	l2->Q2Dmap = NULL;
    }
    else {
	l2->Q2Dmap = phgCalloc(l2->Qsize + 1, sizeof(*l2->Q2Dmap));
	for (i = 0; i < nprocs; i++)
	    l2->Q2Dmap[l2->D2Qmap[i] + 1]++;
	/* update Lsizes_are_equal flag */
	for (i = 2; i <= l2->Qsize; i++)
	    if (l2->Q2Dmap[i] != l2->Q2Dmap[i - 1])
		break;
	l2->Lsizes_are_equal = (i > l2->Qsize);
	for (i = 1; i <= l2->Qsize; i++) 
	    l2->Q2Dmap[i] += l2->Q2Dmap[i - 1];
    }
    i = l2->Lsizes_are_equal;
    MPI_Bcast(&i, 1, MPI_INT, 0, l2->Lcomm);
    l2->Lsizes_are_equal = i;

    if (l2->Lsizes_are_equal) {
	/* Don't need D2Qmap and Q2Dmap in this case */
	phgFree(l2->D2Qmap);
	l2->D2Qmap = NULL;
	phgFree(l2->Q2Dmap);
	l2->Q2Dmap = NULL;
    }

    if (l2->G2Dmap != NULL) {
	/* check if G2Dmap can be freed */
	for (i = 0; i < nprocs; i++)
	    if (l2->G2Dmap[i] != i)
		break;
	if (i >= nprocs) {
	    /* G2Dmap == I, free the array */
	    phgFree(l2->G2Dmap);
	    l2->G2Dmap = NULL;
	}
    }

    /*phgInfo(-1, "Qrank = %d, Lrank = %d, rank_d = %d\n", l2->Qrank, l2->Lrank, l2->G2Dmap[rank]);*/

    l2->ref_count = 1;
#if 1	/* MPI-2 */
    MPI_Comm_set_attr(*comm, l2key, l2);
#else	/* MPI-1 */
    MPI_Attr_put(*comm, l2key, l2);
#endif

#ifndef PHG_PHG_H
    phgRank = rank;
    phgNProcs = nprocs;
    phgComm = *comm;
#endif	/* defined(PHG_PHG_H) */
}

static void
pack_counters(int n, const int data[], int buffer[], int *position)
/* packs data[] into buffer[] */
{
    int i, j = *position, k, m;

    buffer += *position;
    i = j = 0;
    while (i < n) {
	buffer[j++] = k = data[i++];
	m = i;
	while (i < n && data[i] == k)
	    i++;
	if ((m = i - m) > 0)
	    buffer[j++] = -m;
    }
    *position += j;
}

L2DATA *
MPI_L2check(MPI_Comm comm)
/* returns data pointer if 'comm' has the attribute l2key, NULL otherwise */
{
    int ret;
    void *data;

    if (l2key == MPI_KEYVAL_INVALID)
	return NULL;

#if 1	/* MPI-2 */
    MPI_Comm_get_attr(comm, l2key, &data, &ret);
#else	/* MPI-1 */
    MPI_Attr_get(comm, l2key, &data, &ret);
#endif
    assert(ret == 0 || data != NULL);

    return ret == 0 ? NULL : data;
}

int
MPI_Alltoallv2(
    const void *sbuff, const int *scnts, const int *sdsps, MPI_Datatype stype,
    void *rbuff, const int *rcnts, const int *rdsps, MPI_Datatype rtype,
    MPI_Comm comm)
{
    int i, j, k, l = 0, nprocs = 0;
    /* Note: scl: scnts for local comm, sdl: sdsps for local comm,
     *	     scq: scnts for quotient comm, etc. */
    int *scl = NULL, *sdl = NULL, *rcl = NULL, *rdl = NULL;
    int *scq = NULL, *sdq = NULL, *rcq = NULL, *rdq = NULL;
    int sctn, rctn, *sct = NULL, *sdt = NULL, *rct = NULL, *rdt = NULL;
    int *cnts = NULL, *cntsq, ncntsq, *p1, *p2;
    char *sbuffer = NULL, *rbuffer = NULL, *p;
    size_t ssize, rsize;	/* byte size for stype and rtype */
    size_t count, rcount, scount, size = 0;
    MPI_Datatype type0, type1, type2;
    double time0 = 0, time;
    L2DATA *l2;

    /* workaround with new MPICH's strict buffer alias check */
    if (sbuff == NULL) {
	MPI_Comm_size(comm, &nprocs);
	for (i = 0; i < nprocs; i++)
	    if (scnts[i] != 0)
		break;
	if (i < nprocs) {
	    fprintf(stderr, "MPI_Alltoallv: sbuf==NULL and scnts!=0, abort.\n");
	}
	sbuff = &ssize;	/* any valid location */
    }

    if (rbuff == NULL) {
	if (nprocs <= 0)
	    MPI_Comm_size(comm, &nprocs);
	for (i = 0; i < nprocs; i++)
	    if (rcnts[i] != 0)
		break;
	if (i < nprocs) {
	    fprintf(stderr, "MPI_Alltoallv: rbuf==NULL and rcnts!=0, abort.\n");
	}
	rbuff = &rsize;	/* any valid location */
    }

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Alltoallv.\n", __func__);
	return MPI_Alltoallv(CNSTV sbuff, CNSTV scnts, CNSTV sdsps, stype,
			     rbuff, CNSTV rcnts, CNSTV rdsps, rtype, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    if (nprocs <= 0)
	MPI_Comm_size(comm, &nprocs);

    /* Byte sizes of stype and rtype (WARNING: 2GB limit!).
     * FIXME: use MPI_Type_get_true_extent instead of MPI_Type_size? */
    MPI_Type_size(stype, &i);
    assert(i >= 0);
    ssize = i;
    MPI_Type_size(rtype, &i);
    assert(i >= 0);
    rsize = i;

    /* Step 0: local exchange (TODO: combine it with Step 1?) */
    time = phgGetTime(NULL);
    scl = phgAlloc(4 * l2->Lsize * sizeof(*scl));
    rcl = scl + l2->Lsize;
    sdl = rcl + l2->Lsize;
    rdl = sdl + l2->Lsize;
    for (i = 0; i < l2->Lsize; i++) {
	j = l2->L2Gmap[i];	/* proc 'i' in Lcomm = proc 'j' in comm */
	scl[i] = scnts[j];
	rcl[i] = rcnts[j];
	sdl[i] = sdsps[j];
	rdl[i] = rdsps[j];
    }
    MPI_Alltoallv(CNSTV sbuff, scl, sdl, stype, rbuff, rcl, rdl, rtype,
		  l2->Lcomm);
    (void)time0;	/* to avoid gcc warning */
    time0 = time; time = phgGetTime(NULL);
    /*phgPrintf("%s: local Alltoallv time = %g\n", __func__, time - time0);*/

#if 0	/*-------------------------------------------------------------------*/
    /* short-cut test: simply call MPI_Alltoallv on the remaining processes */
    phgFree(scl);
    sct = phgAlloc(4 * nprocs * sizeof(*sct));
    rct = sct + nprocs;
    sdt = rct + nprocs;
    rdt = sdt + nprocs;
    memcpy(sct, scnts, nprocs * sizeof(*sct));
    memcpy(rct, rcnts, nprocs * sizeof(*rct));
    memcpy(sdt, sdsps, nprocs * sizeof(*sdt));
    memcpy(rdt, rdsps, nprocs * sizeof(*rdt));
    for (i = 0; i < l2->Lsize; i++) {
	j = l2->L2Gmap[i];
	sct[j] = rct[j] = 0;
    }
    i = MPI_Alltoallv(sbuff, sct, sdt, stype, rbuff, rct, rdt, rtype, comm);
    phgFree(sct);
    time0 = time; time = phgGetTime(NULL);
    /*phgPrintf("%s: global Alltoallv time = %g\n", __func__, time - time0);*/
    return i;
#endif	/*-------------------------------------------------------------------*/

    /* allocate buffers, the 2 extra entries are for ssize and rsize in the
     * packed buffer (cntsq[]) */
    scq = phgAlloc(4 * (l2->Qsize + nprocs) * sizeof(*scq) + 2);
    rcq = scq + l2->Qsize;
    /* Note: sdq and rdq are only needed on the master of Lcomm */
    sdq = rcq + l2->Qsize;
    rdq = sdq + l2->Qsize;
    sct = rdq + l2->Qsize;
    rct = sct + nprocs;
    sdt = rct + nprocs;
    rdt = sdt + nprocs;
    /* Note: cntsq[] is first overlaped with sdt and rdt, and is used to hold
     * locally packed counters. Afterwards it will point to gathered scq[] */
    cntsq = sdt;

    /* Step 1: local gather (in sbuffer) */

    /* collect data blocks to send to local master */
    for (i = 0; i < nprocs; i++) {
	j = G2DMAP(l2, i);
	if (D2QMAP(l2, j) != l2->Qrank) {
	    sct[j] = scnts[i];
	    rct[j] = rcnts[i];
	}
	else {
	    sct[j] = rct[j] = 0;
	}
    }

    /* Note: cntsq[] is the buffer containing the counters to send to root_q,
     * which has the following structure:
     *		ssize, rsize, packed rct[], packed scq[]
     */
    ncntsq = 0;
    cntsq[ncntsq++] = (int)ssize;
    cntsq[ncntsq++] = (int)rsize;

    /* pack rct[] (rcnts[] with local entries reset) */
    pack_counters(nprocs, rct, cntsq, &ncntsq);

    /* compute aggregated counters scq[] */
    memset(scq, 0, l2->Qsize * sizeof(*scq));
    for (i = 0; i < nprocs; i++) {
	if (sct[i] != 0)
	    scq[D2QMAP(l2, i)] += sct[i];
    }
    /* pack scq[] (scnts[] with local entries reset, aggregated w.r.t Lcomm) */
    pack_counters(l2->Qsize, scq, cntsq, &ncntsq);

    /* NOTE: with mvapich2-1.6_gcc-4.6 on LSSC3, the the following Gather may
     * yield wrong result. For example, in Pass 3 of the following run:
     *		mpirun -np 16 ./poisson -mpi_2level_ppn 8
     * p10 sends 9 to p8, p12 sends 11 to p8, while p8 gets 11 (should be 9) in
     * rcl[2] and 9 (should be 11) in rcl[4]. */

    /* gather conters (FIXME: combine Gather and Gatherv?) */
    MPI_Gather(&ncntsq, 1, MPI_INT, rcl, 1, MPI_INT, 0, l2->Lcomm);
    if (l2->Lrank == 0) {
	size = rcl[0];
	rdl[0] = 0;
	for (i = 1; i < l2->Lsize; i++) {
	    assert(size <= INT_MAX);
	    rdl[i] = (int)size;
	    size += rcl[i];
	}
	assert(size <= INT_MAX);
	/*phgInfo(-1, "Gathered data size: %d (v.s. %d)\n", (int)(size - 2 * l2->Lsize), l2->Lsize * l2->Qsize + l2->Lsize * nprocs);*/
	cnts = phgAlloc(size * sizeof(*cnts));
    }
    MPI_Gatherv(cntsq, ncntsq, MPI_INT,
		cnts, rcl, rdl, MPI_INT, 0, l2->Lcomm);
    time0 = time; time = phgGetTime(NULL);
    /*phgPrintf("%s: local Gather time = %g\n", __func__, time - time0);*/
    if (l2->Lrank == 0) {
	/* rearrange the received blocks from
	 * 	(rcnts,scq), (rcnts,scq), ...,
	 * to:
	 *	(scq,scq,...), (rcnts,rcnts,...)
	 */
	int *tmp, n = 0, nq = 0;
	/* Compute size of data for rcnts[] (n) and scq[] (nq) respectively */
	p1 = cnts;
	for (i = 0; i < l2->Lsize; i++) {
	    p1 += 2;	/* skip ssize and rsize */
	    /* step over rcnts[] from process i */
	    p2 = p1;
	    j = 0;
	    while (j < nprocs) {
		if ((k = *(p1++)) >= 0)
		    j++;
		else
		    j += -k;
	    }
	    if (j != nprocs) {
		phgError(1, "(%s:%d) unexpected.\n", __FILE__, __LINE__);
	    }
	    n += (int)(p1 - p2);
	    /* step over scq[] from process i */
	    p2 = p1;
	    j = 0;
	    while (j < l2->Qsize) {
		if ((k = *(p1++)) >= 0)
		    j++;
		else
		    j += -k;
	    }
	    if (j != l2->Qsize) {
		phgError(1, "(%s:%d) unexpected.\n", __FILE__, __LINE__);
	    }
	    nq += (int)(p1 - p2);
	}
	assert(nq + n == size - 2 * l2->Lsize);
	tmp = cnts;
	cnts = phgAlloc((n + nq) * sizeof(*tmp));
	cntsq = cnts + n;
	/* copy counters */
	p1 = tmp;
	n = nq = 0;
	for (i = 0; i < l2->Lsize; i++) {
	    size_t s, r;
	    s = *(p1++);		/* ssize on proc i */
	    r = *(p1++);		/* rsize on proc i */
	    /* step over rcnts[] from process i */
	    j = 0;
	    while (j < nprocs) {
		if ((k = *(p1++)) >= 0) {
		    j++;
		    /* transform: k := (k*r)/rsize (FIXME: avoid overflow!) */
		    size = (k / rsize) * r + ((k % rsize) * r) / rsize;
		    assert(((k % rsize) * r) % rsize == 0);
		    assert(size <= INT_MAX);
		    cnts[n++] = (int)size;
		}
		else {
		    j += -k;
		    cnts[n++] = k;
		}
	    }
	    assert(j == nprocs);
	    /* step over scq[] from process i */
	    j = 0;
	    while (j < l2->Qsize) {
		if ((k = *(p1++)) >= 0) {
		    j++;
		    /* transform: k := (k*s)/ssize (FIXME: avoid overflow!) */
		    size = (k / ssize) * s + ((k % ssize) * s) / ssize;
		    assert(((k % ssize) * s) % ssize == 0);
		    assert(size <= INT_MAX);
		    cntsq[nq++] = (int)size;
		}
		else {
		    j += -k;
		    cntsq[nq++] = k;
		}
	    }
	    assert(j == l2->Qsize);
	}
	phgFree(tmp);
    }

    /* compact arrays sct, sdt, rct and rdt */
    for (i = 0; i < nprocs; i++) {
	j = G2DMAP(l2, i);
	sdt[j] = sdsps[i];
	rdt[j] = rdsps[i];
    }
    for (i = j = k = 0; i < nprocs; i++) {
	if (sct[i] != 0) {
	    if (j != i) {
		sct[j] = sct[i];
		sdt[j] = sdt[i];
	    }
	    j++;
	}
	if (rct[i] != 0) {
	    if (k != i) {
		rct[k] = rct[i];
		rdt[k] = rdt[i];
	    }
	    k++;
	}
    }
    sctn = j;
    MPI_Type_indexed(sctn, sct, sdt, stype, &type1);	/* local gather type */
    MPI_Type_commit(&type1);
    rctn = k;
    MPI_Type_indexed(rctn, rct, rdt, rtype, &type2);	/* local scatter type */
    MPI_Type_commit(&type2);

    if (l2->Lrank == 0) {
	/* Compute aggregated sums
	 * 	rct[j] = sum over Lcomm rcnts[][j]
	 * 	scl[i] = sum over comm rcnts[i][]
	 *	rcl[i] = sum over comm scnts[i][]
	 *		0 <= i < Lsize
	 *
	 * 	scq[i] = sum over Lcomm[i] scnts[]
	 *	rcq[i] = sum over Lcomm[i] rcnts[]
	 *		0 <= i < Qsize
	 *
         * Note: scl and rcl are contiguous (so are scq and rcq)
	 */
	memset(rct, 0, nprocs * sizeof(*rct));
	memset(scl, 0, 2 * l2->Lsize * sizeof(*scl));
	memset(scq, 0, 2 * l2->Qsize * sizeof(*scq));
	p1 = cnts;
	p2 = cntsq;
	for (i = 0; i < l2->Lsize; i++) {
	    /* loop on rcnts[] */
	    j = 0;
	    while (j < nprocs) {
		if ((k = *(p1++)) >= 0) {
		    if ((l = k) > 0) {
			scl[i] += l;			/* local Scatterv */
			rcq[D2QMAP(l2, j)] += l;	/* global Alltoallv */
			rct[j] += l;			/* for rearrange data */
		    }
		    j++;
		}
		else {
		    k = -k;
		    if (l > 0) {
			while (k-- > 0) {
			    scl[i] += l;
			    rcq[D2QMAP(l2, j)] += l;
			    rct[j] += l;
			    j++;
			}
		    }
		    else {
			j += k;
		    }
		}
	    }
	    assert(j == nprocs);
	    /* loop on scnts[] */
	    j = 0;
	    while (j < l2->Qsize) {
		if ((k = *(p2++)) >= 0) {
		    if ((l = k) > 0) {
			rcl[i] += l;			/* local Gatherv */
			scq[j] += l;			/* global Alltoallv */
		    }
		    j++;
		}
		else {
		    k = -k;
		    if (l > 0) {
			while (k-- > 0) {
			    rcl[i] += l;
			    scq[j] += l;
			    j++;
			}
		    }
		    else {
			j += k;
		    }
		}
	    }
	    assert(j == l2->Qsize);
	}
	/* compute displacements for local Gatherv and Scatterv */
	rcount = rcl[0];
	scount = scl[0];
	rdl[0] = sdl[0] = 0;
	for (i = 1; i < l2->Lsize; i++) {
	    rdl[i] = (int)rcount;
	    sdl[i] = (int)scount;
	    assert(rcount <= INT_MAX && scount <= INT_MAX);
	    rcount += (size_t)rcl[i];
	    scount += (size_t)scl[i];
	}
	/* compute displacements for global Alltoallv */
	rcount = rcq[0];
	scount = scq[0];
	rdq[0] = sdq[0] = 0;
	for (j = 1; j < l2->Qsize; j++) {
	    rdq[j] = (int)rcount;
	    sdq[j] = (int)scount;
	    assert(rcount <= INT_MAX && scount <= INT_MAX);
	    rcount += (size_t)rcq[j];
	    scount += (size_t)scq[j];
	}
	scount *= ssize;
	rcount *= rsize;
	count = (rcount < scount ? scount : rcount);
	sbuffer = phgAlloc(count + count);
	rbuffer = sbuffer + count;
	/* compute displacements for rearrange received local data */
	size = rct[0];
	rdt[0] = 0;
	for (j = 1; j < nprocs; j++) {
	    assert(size < INT_MAX);
	    rdt[j] = (int)size;
	    size += rct[j];
	}
    }
    MPI_Type_contiguous((int)ssize, MPI_BYTE, &type0);
    MPI_Type_commit(&type0);
    MPI_Gatherv(CNSTV sbuff, sctn == 0 ? 0 : 1, type1,
		rbuffer, rcl, rdl, type0, 0, l2->Lcomm);
    time0 = time; time = phgGetTime(NULL);
    /*phgPrintf("%s: local Gatherv time = %g\n", __func__, time - time0);*/

    /* Step 2: global exchange */
    MPI_Type_free(&type1);
    MPI_Type_contiguous((int)rsize, MPI_BYTE, &type1);
    MPI_Type_commit(&type1);
    if (l2->Lrank == 0) {
	/* rearrange data in rbuffer[] and put them in sbuffer[] */
	p = rbuffer;
	p2 = cntsq;
	for (i = 0; i < l2->Lsize; i++) {
	    j = 0;
	    while (j < l2->Qsize) {
		if ((k = *(p2++)) >= 0) {
		    if ((l = k) > 0) {
			size = l * ssize;
			memcpy(sbuffer + sdq[j] * ssize, p, size);
			p += size;
			sdq[j] += l;
		    }
		    j++;
		}
		else {
		    k = -k;
		    if (l > 0) {
		 	while (k-- > 0) {
			    memcpy(sbuffer + sdq[j] * ssize, p, size);
			    p += size;
			    sdq[j] += l;
			    j++;
		 	}
		    }
		    else {
			j += k;
		    }
		}
	    }
	    assert(j == l2->Qsize);
	}
	/* restore sdq[] */
	for (j = 0; j < l2->Qsize; j++)
	    sdq[j] -= scq[j];

	MPI_Alltoallv(sbuffer, scq, sdq, type0, rbuffer, rcq, rdq, type1,
			l2->Qcomm);
	time0 = time; time = phgGetTime(NULL);
	/*phgPrintf("%s: global Alltoallv time = %g\n", __func__, time-time0);*/

	/* rearrange data in rbuffer[] and put them in sbuffer[] */
	p = sbuffer;
	p1 = cnts;
	for (i = 0; i < l2->Lsize; i++) {
	    j = 0;
	    while (j < nprocs) {
		if ((k = *(p1++)) >= 0) {
		    if ((l = k) > 0) {
			size = l * rsize;
			memcpy(p, rbuffer + rdt[j] * rsize, size);
			p += size;
			rdt[j] += l;
		    }
		    j++;
		}
		else {
		    k = -k;
		    if (l > 0) {
		 	while (k-- > 0) {
			    memcpy(p, rbuffer + rdt[j] * rsize, size);
			    p += size;
			    rdt[j] += l;
			    j++;
		 	}
		    }
		    else {
			j += k;
		    }
		}
	    }
	    assert(j == nprocs);
	}
	phgFree(cnts);
    }

    /* Step 3: local scatter */

    /* Note: MPI_Scatterv on LSSC3 seems to hang if sizeof(recvtype)==0,
     *	     as a work around we set recvcount=0 in this case */
    MPI_Scatterv(sbuffer, scl, sdl, type1,
		 rbuff, rctn == 0 ? 0 : 1, type2, 0, l2->Lcomm);
    time0 = time; time = phgGetTime(NULL);
    /*phgPrintf("%s: local Scatterv time = %g\n", __func__, time - time0);*/

    MPI_Type_free(&type0);
    MPI_Type_free(&type1);
    MPI_Type_free(&type2);

    phgFree(sbuffer);
    phgFree(scl);
    phgFree(scq);

    return MPI_SUCCESS;
}

int
MPI_Alltoall2(const void *sbuff, int scnt, MPI_Datatype stype,
	      void *rbuff, int rcnt, MPI_Datatype rtype, MPI_Comm comm)
{
    int i, j, k, nprocs;
    int *cnts = NULL, *dsps = NULL;
    char *buff = NULL, *buff0 = NULL, *p, *p0;
    size_t size = 0, rsize;
    MPI_Datatype type = MPI_DATATYPE_NULL;
    L2DATA *l2;

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Alltoall.\n", __func__);
	return MPI_Alltoall(CNSTV sbuff, scnt, stype, rbuff, rcnt, rtype, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    MPI_Comm_size(comm, &nprocs);

    MPI_Type_size(rtype, &i);
    assert(i >= 0);
    rsize = i;
    rsize *= rcnt;

    MPI_Type_size(stype, &i);
    assert(i >= 0);
    size = i;
    size *= scnt;
    assert(size == rsize);

    if (rsize == 0)
	return MPI_SUCCESS;

    if (l2->G2Dmap != NULL)
	i = nprocs;
    else if (l2->Lrank == 0 && !l2->Lsizes_are_equal)
	i = l2->Qsize + 2;
    else
	i = 2;
    if (i > 0) {
	cnts = phgAlloc(2 * nprocs * sizeof(*cnts));
	dsps = cnts + nprocs;
    }

    /* local Alltoall */
    if (l2->G2Dmap == NULL) {
	/* procs of comm are contiguous in Lcomm, do Alltoall */
#if 1	/* MPI-2 */
	MPI_Aint lb, sext, rext;
	MPI_Type_get_extent(stype, &lb, &sext);
	MPI_Type_get_extent(rtype, &lb, &rext);
#else	/* MPI-1 */
	MPI_Aint sext, rext;
	MPI_Type_extent(stype, &sext);
	MPI_Type_extent(rtype, &rext);
#endif
	i = l2->L2Gmap[0];
	assert(l2->Lsize == l2->L2Gmap[l2->Lsize - 1] - i + 1);
	MPI_Alltoall(CNSTV (((const char *)sbuff) + scnt * sext * i), scnt,
		     stype, ((char *)rbuff) + rcnt * rext * i, rcnt, rtype,
		     l2->Lcomm);
    }
    else {
	/* do Alltoallv */
	for (i = 0; i < l2->Lsize; i++) {
	    j = l2->L2Gmap[i];
	    cnts[i] = scnt;
	    assert(scnt * (size_t)j <= INT_MAX);
	    dsps[i] = scnt * j;
	}
	MPI_Alltoallv(CNSTV sbuff, cnts, dsps, stype, rbuff, cnts, dsps, rtype,
		      l2->Lcomm);
    }

    /* local Gather */
    if (l2->Lrank == 0) {
	buff = phgAlloc(2 * (nprocs - l2->Lsize) * rsize * l2->Lsize);
	buff0 = buff + (nprocs - l2->Lsize) * rsize * l2->Lsize;
	MPI_Type_contiguous(rsize, MPI_BYTE, &type);
	MPI_Type_commit(&type);
    }
    if (l2->G2Dmap != NULL) {
	for (i = 0; i < nprocs; i++) {
	    j = l2->G2Dmap[i];
	    cnts[j] = scnt;
	    assert(scnt * (size_t)i <= INT_MAX);
	    dsps[j] = scnt * i;
	}
	for (i = 0; i < l2->Lsize; i++)
	    cnts[l2->G2Dmap[l2->L2Gmap[i]]] = 0;
	/* compact cnts[]/dsps[] */
	for (i = k = 0; i < nprocs; i++) {
	    if (cnts[i] == 0)
		continue;
	    if (k != i) {
		cnts[k] = cnts[i];
		dsps[k] = dsps[i];
	    }
	    k++;
	}
	assert(k == nprocs - l2->Lsize);
	MPI_Type_indexed(k, cnts, dsps, stype, &stype);
    }
    else {
	i = l2->L2Gmap[0];
	k = nprocs - l2->Lsize;
	assert(k * (size_t)scnt <= INT_MAX);
	if (i == 0) {
	    cnts[0] = k * scnt;
	    assert(l2->Lsize * (size_t)scnt <= INT_MAX);
	    dsps[0] = l2->Lsize * scnt;
	    MPI_Type_indexed(1, cnts, dsps, stype, &stype);
	}
	else if (i + l2->Lsize == nprocs) {
	    MPI_Type_contiguous(k * scnt, stype, &stype);
	}
	else {
	    cnts[0] = i * scnt;
	    dsps[0] = 0;
	    cnts[1] = (k - i) * scnt;
	    assert((i + l2->Lsize) * (size_t)scnt <= INT_MAX);
	    dsps[1] = (i + l2->Lsize) * scnt;
	    MPI_Type_indexed(2, cnts, dsps, stype, &stype);
	}
    }
    MPI_Type_commit(&stype);
#if defined(OMPI_MAJOR_VERSION) && \
    (OMPI_MAJOR_VERSION < 1 || (OMPI_MAJOR_VERSION == 1 && \
    (OMPI_MINOR_VERSION < 4 || (OMPI_MINOR_VERSION == 4 && \
    OMPI_RELEASE_VERSION <= 1))))
    /* A workaround for a possible bug in OpenMPI-1.4.1 on LSSC3. Without it
     * 'mpirun -np # ./poisson' would fail for '#'>=160 with the error message
     * 'An error occurred in MPI_Gather'. We believe that this is because the
     * MPI_Gather function of OpenMPI-1.4.1 can't handle complex 'sendtype'. */
    {
	MPI_Status status;
	if (l2->Lrank != 0) {
	    buff0 = phgAlloc((nprocs - l2->Lsize) * rsize);
	    MPI_Type_contiguous(rsize, MPI_BYTE, &type);
	    MPI_Type_commit(&type);
	}
	MPI_Sendrecv(CNSTV sbuff, 1, stype, 0, 111,
			   buff0, k,  type, 0, 111, MPI_COMM_SELF, &status);
	MPI_Gather(buff0, k, type, buff, k, type, 0, l2->Lcomm);
	if (l2->Lrank != 0) {
	    phgFree(buff0);
	    MPI_Type_free(&type);
	}
    }
#else
    MPI_Gather(CNSTV sbuff, 1, stype, buff, k, type, 0, l2->Lcomm);
#endif
    MPI_Type_free(&stype);

    if (l2->Lrank == 0) {
	/* rearrange data: buff -> buff0 */
	p = buff;
	p0 = buff0;
	for (j = 0; j < l2->Qsize; j++) {
	    if (j == l2->Qrank)
		continue;
	    k = LSIZES(l2, j);
	    for (i = 0; i < l2->Lsize; i++) {
		memcpy(p0, p + i * (nprocs - l2->Lsize) * rsize, k * rsize);
		p0 += k * rsize;
	    }
	    p += k * rsize;
	}

	size = 0;
	for (i = 0; i < l2->Qsize; i++) {
	    assert(size <= INT_MAX);
	    dsps[i] = (int)size;
	    if (i == l2->Qrank) {
		cnts[i] = 0;
	    }
	    else {
		assert(LSIZES(l2, i) * (size_t)l2->Lsize <= INT_MAX);
		cnts[i] = LSIZES(l2, i) * l2->Lsize;
		size += cnts[i];
	    }
	}
	MPI_Alltoallv(buff0, cnts, dsps, type, buff, cnts, dsps, type,
		      l2->Qcomm);

	/* rearrange data: buff -> buff0 */
	if (l2->G2Dmap != NULL) {
	    for (i = 0; i < nprocs; i++)
		cnts[i] = 1;
	    for (i = 0; i < l2->Lsize; i++)
		cnts[l2->G2Dmap[l2->L2Gmap[i]]] = 0;
	}
	p = buff;
	p0 = buff0;
	for (j = 0; j < l2->Lsize; j++) {
	    for (i = k = 0; i < nprocs; i++) {
		if (l2->G2Dmap != NULL) {
		    if (cnts[i] == 0)
			continue;
		}
		else {
		    if (i >= l2->L2Gmap[0] &&
			i < l2->L2Gmap[0] + l2->Lsize)
			continue;
		}
		memcpy(p0, p + k * l2->Lsize * rsize, rsize);
		p0 += rsize;
		k++;
	    }
	    p += rsize;
	}
    }

    /* local Scatter */
    if (l2->G2Dmap != NULL) {
	for (i = 0; i < nprocs; i++) {
	    j = l2->G2Dmap[i];
	    cnts[j] = rcnt;
	    assert(rcnt * (size_t)i <= INT_MAX);
	    dsps[j] = rcnt * i;
	}
	for (i = 0; i < l2->Lsize; i++)
	    cnts[l2->G2Dmap[l2->L2Gmap[i]]] = 0;
	for (i = k = 0; i < nprocs; i++) {
	    if (cnts[i] == 0)
		continue;
	    if (k != i) {
		cnts[k] = cnts[i];
		dsps[k] = dsps[i];
	    }
	    k++;
	}
	MPI_Type_indexed(k, cnts, dsps, rtype, &rtype);
    }
    else {
	i = l2->L2Gmap[0];
	k = nprocs - l2->Lsize;
	assert(k * (size_t)rcnt <= INT_MAX);
	if (i == 0) {
	    cnts[0] = k * rcnt;
	    assert(l2->Lsize * (size_t)rcnt <= INT_MAX);
	    dsps[0] = l2->Lsize * rcnt;
	    MPI_Type_indexed(1, cnts, dsps, rtype, &rtype);
	}
	else if (i + l2->Lsize == nprocs) {
	    MPI_Type_contiguous(k * rcnt, rtype, &rtype);
	}
	else {
	    cnts[0] = i * rcnt;
	    dsps[0] = 0;
	    cnts[1] = (k - i) * rcnt;
	    assert((i + l2->Lsize) * (size_t)rcnt <= INT_MAX);
	    dsps[1] = (i + l2->Lsize) * rcnt;
	    MPI_Type_indexed(2, cnts, dsps, rtype, &rtype);
	}
    }
    MPI_Type_commit(&rtype);
    MPI_Scatter(buff0, k, type, rbuff, 1, rtype, 0, l2->Lcomm);
    MPI_Type_free(&rtype);

    phgFree(cnts);
    if (l2->Lrank == 0) {
	MPI_Type_free(&type);
	phgFree(buff);
    }

    return MPI_SUCCESS;
}

int
MPI_Allgatherv2(const void *sbuff, int count, MPI_Datatype stype,
		void *rbuff, const int *rcnts, const int *rdsps,
		MPI_Datatype rtype, MPI_Comm comm)
{
    int i, j, nprocs;
    int *cnts = NULL, *dsps = NULL;
    char *buff = NULL;
    size_t size = 0, rsize;
    MPI_Datatype type = MPI_DATATYPE_NULL;
    L2DATA *l2;

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Allgatherv.\n", __func__);
	return MPI_Allgatherv(CNSTV sbuff, count, stype,
			      rbuff, CNSTV rcnts, CNSTV rdsps, rtype, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    MPI_Comm_size(comm, &nprocs);

    MPI_Type_size(rtype, &i);
    assert(i >= 0);
    rsize = i;

    if (l2->Lrank == 0) {
	/* FIXME: use MPI_Type_get_true_extent instead of MPI_Type_size? */
	MPI_Type_contiguous((int)rsize, MPI_BYTE, &type);
	MPI_Type_commit(&type);

	if (l2->G2Dmap == NULL)
	    i = (l2->Lsize > l2->Qsize ? l2->Lsize : l2->Qsize);
	else
	    i = nprocs;
	cnts = phgAlloc(2 * i * sizeof(*cnts));
	dsps = cnts + i;

	size = cnts[0] = rcnts[l2->L2Gmap[0]];
	dsps[0] = 0;
	for (i = 1; i < l2->Lsize; i++) {
	    assert(size <= INT_MAX);
	    dsps[i] = (int)size;
	    size += (cnts[i] = rcnts[l2->L2Gmap[i]]);
	}
	buff = phgAlloc(size * rsize);
    }

    /* Step 0: local Gather */
    MPI_Gatherv(CNSTV sbuff, count, stype, buff, cnts, dsps, type, 0,
		l2->Lcomm);

    /* Step 1: global Allgatherv */
    if (l2->Lrank == 0) {
	char *buff1;
	MPI_Status status;

	assert(size <= INT_MAX);
	count = (int)size;
	memset(cnts, 0, l2->Qsize * sizeof(*cnts));
	for (i = 0; i < nprocs; i++) {
	    j = D2QMAP(l2, G2DMAP(l2, i));
	    cnts[j] += rcnts[i];
	}
	size = cnts[0];
	dsps[0] = 0;
	for (i = 1; i < l2->Qsize; i++) {
	    assert(size <= INT_MAX);
	    dsps[i] = (int)size;
	    size += cnts[i];
	}
	buff1 = phgAlloc(size * rsize);
	MPI_Allgatherv(buff, count, type, buff1, cnts, dsps, type,
		       l2->Qcomm);
	/* Note: rearrange data according to G2Dmap[]i and rcnts/rdsps/rtype */
	/* Note: stype (who's no longer needed) is used as scratch variable */
	assert(size <= INT_MAX);
	count = (int)size;
	stype = type;
	if (l2->G2Dmap != NULL) {
	    for (i = 0; i < nprocs; i++) {
		j = l2->G2Dmap[i];
		cnts[j] = rcnts[i];
		dsps[j] = rdsps[i];
	    }
	    MPI_Type_indexed(nprocs, cnts, dsps, rtype, &type);
	}
	else {
	    MPI_Type_indexed(nprocs, CNSTV rcnts, CNSTV rdsps, rtype, &type);
	}
	MPI_Type_commit(&type);
	/* replace gathered data into rbuff */
	MPI_Sendrecv(buff1, count, stype, 0, 123, rbuff, 1, type, 0, 123,
		     MPI_COMM_SELF, &status);
	MPI_Type_free(&stype);
	if (l2->G2Dmap != NULL)
	    MPI_Type_free(&type);
	phgFree(buff);
	phgFree(buff1);
	phgFree(cnts);
    }

    if (type == MPI_DATATYPE_NULL) {
	MPI_Type_indexed(nprocs, CNSTV rcnts, CNSTV rdsps, rtype, &type);
	MPI_Type_commit(&type);
    }

    /* Step 2: local Bcast */
    MPI_Bcast(rbuff, 1, type, 0, l2->Lcomm);
    MPI_Type_free(&type);

    return MPI_SUCCESS;
}

int
MPI_Allgather2(const void *sbuff, int scnt, MPI_Datatype stype,
		void *rbuff, int rcnt, MPI_Datatype rtype, MPI_Comm comm)
{
    int nprocs;
    size_t rsize;
    char *buff = NULL;
    MPI_Datatype type = MPI_DATATYPE_NULL;
    L2DATA *l2;

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Allgather.\n", __func__);
	return MPI_Allgather(CNSTV sbuff,scnt,stype, rbuff,rcnt,rtype, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    MPI_Type_size(rtype, &nprocs);
    assert(nprocs >= 0);
    rsize = nprocs;

    MPI_Comm_size(comm, &nprocs);

    /* Local gather */
    if (l2->Lrank == 0) {
	buff = phgAlloc(rcnt * rsize * l2->Lsize);
	MPI_Type_contiguous(rsize, MPI_BYTE, &type);
	MPI_Type_commit(&type);
    }
    MPI_Gather(CNSTV sbuff, scnt, stype, buff, rcnt, type, 0, l2->Lcomm);

    /* Global ALlgather */
    if (l2->Lrank == 0) {
	if (l2->G2Dmap == NULL && l2->Lsizes_are_equal) {
	    MPI_Allgather(buff, rcnt * l2->Lsize, type,
			  rbuff, rcnt * l2->Lsize, rtype, l2->Qcomm);
	}
	else {
	    /* procs not in order or different Lsize: use Allgatherv */
	    size_t size;
	    int i, *cnts, *dsps;

	    i = (l2->G2Dmap == NULL ? l2->Qsize : nprocs);
	    cnts = phgAlloc(2 * i * sizeof(*cnts));
	    dsps = cnts + i;
	    size = 0;
	    for (i = 0; i < l2->Qsize; i++) {
		assert(LSIZES(l2, i) * (size_t)rcnt <= INT_MAX);
		cnts[i] = LSIZES(l2, i) * rcnt;
		assert(size <= INT_MAX);
		dsps[i] = (int)size;
		size += cnts[i];
	    }
	    if (l2->G2Dmap == NULL) {
		MPI_Allgatherv(buff, l2->Lsize * rcnt, type,
			       rbuff, cnts, dsps, rtype, l2->Qcomm);
	    }
	    else {
		char *buff1 = phgAlloc(size * rsize);
		MPI_Status status;
		MPI_Allgatherv(buff, l2->Lsize * rcnt, type,
			       buff1, cnts, dsps, type, l2->Qcomm);
		/* reorder data */
		for (i = 0; i < nprocs; i++) {
		    cnts[i] = rcnt;
		    size = l2->G2Dmap[i] * (size_t)rcnt;
		    assert(size <= INT_MAX);
		    dsps[i] = (int)size;
		}
		MPI_Type_indexed(nprocs, cnts, dsps, type, &stype);
		MPI_Type_commit(&stype);
		assert(nprocs * (size_t)rcnt <= INT_MAX);
		MPI_Sendrecv(buff1, 1, stype, 0, 123,
			     rbuff, nprocs * rcnt, rtype, 0, 123,
			     MPI_COMM_SELF, &status);
		MPI_Type_free(&stype);
		phgFree(buff1);
	    }
	    phgFree(cnts);
	}
	MPI_Type_free(&type);
	phgFree(buff);
    }

    /* Local Bcast */
    if (rcnt * (size_t)nprocs <= INT_MAX) {
	MPI_Bcast(rbuff, rcnt * nprocs, rtype, 0, l2->Lcomm);
    }
    else {
	MPI_Type_contiguous(nprocs, rtype, &type);
	MPI_Type_commit(&type);
	MPI_Bcast(rbuff, rcnt, type, 0, l2->Lcomm);
	MPI_Type_free(&type);
    }

    return MPI_SUCCESS;
}

int
MPI_Gatherv2(const void *sbuff, int scnt, MPI_Datatype stype,
	     void *rbuff, const int *rcnts, const int *rdsps,
	     MPI_Datatype rtype, int root, MPI_Comm comm)
{
    int i, j, rank, nprocs, root_q;
    int *cnts = NULL, *dsps = NULL;
    size_t ssize, rsize, size = 0;
    MPI_Datatype type = MPI_DATATYPE_NULL;
    MPI_Status status;
    char *buff = NULL;
    L2DATA *l2;

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Gatherv.\n", __func__);
	return MPI_Gatherv(CNSTV sbuff, scnt, stype,
			   rbuff, CNSTV rcnts, CNSTV rdsps, rtype,
			   root, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (rank != root)
	rtype = stype;

    MPI_Type_size(stype, &i);
    ssize = i;
    MPI_Type_size(rtype, &i);
    rsize = i;
    root_q = D2QMAP(l2, G2DMAP(l2, root));

    /* Step 1: local Gatherv */
    if (l2->Lrank == 0) {
	if (l2->Qrank != root_q)
	    i = l2->Lsize;
	else if (rank != root)
	    i = 2 * nprocs;
	else
	    i = nprocs;
	cnts = phgAlloc(2 * i * sizeof(*cnts));
	dsps = cnts + i;
    }
    assert(scnt * ssize <= INT_MAX);
    i = scnt * ssize;
    MPI_Gather(&i, 1, MPI_INT, cnts, 1, MPI_INT, 0, l2->Lcomm);
    if (l2->Lrank == 0) {
	size = 0;
	for (i = 0; i < l2->Lsize; i++) {
	    assert(cnts[i] % rsize == 0);
	    cnts[i] /= (int)rsize;
	    assert(size <= INT_MAX);
	    dsps[i] = (int)size;
	    size += cnts[i];
	}
	buff = phgAlloc(size * rsize);
	MPI_Type_contiguous(rsize, MPI_BYTE, &type);
	MPI_Type_commit(&type);
    }
    MPI_Gatherv(CNSTV sbuff, scnt, stype, buff, cnts, dsps, type, 0,
		l2->Lcomm);

    /* Step 2: global Gatherv */
    if (l2->Lrank == 0) {
	assert(size <= INT_MAX);
	scnt = (int)size;
	if (l2->Qrank != root_q) {
	    MPI_Gatherv(buff, (int)size, type, NULL, NULL, NULL, rtype,
			root_q, l2->Qcomm);
	    MPI_Type_free(&type);
	}
	else {
	    char *buff1;
	    if (rank != root) {
		/* need to get counts from root */
		rcnts = cnts + nprocs;
		rdsps = dsps + nprocs;
		MPI_Recv((void *)rcnts, nprocs, MPI_INT, root, 123,
			 comm, &status);
	    }
	    memset(cnts, 0, l2->Qsize * sizeof(*cnts));
	    for (i = 0; i < nprocs; i++)
		cnts[D2QMAP(l2, G2DMAP(l2, i))] += rcnts[i];
	    size = 0;
	    for (i = 0; i < l2->Qsize; i++) {
		assert(size <= INT_MAX);
		dsps[i] = size;
		size += cnts[i];
	    }
	    buff1 = phgAlloc(size * rsize);
	    MPI_Gatherv(buff,scnt,type, buff1,cnts,dsps,type, root_q,
			l2->Qcomm);
	    /* Send data to root */
	    assert(size <= INT_MAX);
	    scnt = (int)size;
	    stype = type;
	    if (rank == root) {
		if (l2->G2Dmap == NULL) {
		    MPI_Type_indexed(nprocs, CNSTV rcnts, CNSTV rdsps,
				     rtype, &type);
		}
		else {
		    for (i = 0; i < nprocs; i++) {
			j = l2->G2Dmap[i];
			cnts[j] = rcnts[i];
			dsps[j] = rdsps[i];
		    }
		    MPI_Type_indexed(nprocs, cnts, dsps, rtype, &type);
		}
		MPI_Type_commit(&type);
		MPI_Sendrecv(buff1, scnt, stype, 0, 123, rbuff, 1, type, 0, 123,
			     MPI_COMM_SELF, &status);
		MPI_Type_free(&type);
	    }
	    else {
		MPI_Send(buff1, scnt, stype, root, 123, comm);
	    }
	    MPI_Type_free(&stype);
	    phgFree(buff1);
	}
	phgFree(cnts);
	phgFree(buff);
    }
    else if (rank == root) {
	/* send counts to local master */
	MPI_Send(CNSTV rcnts, nprocs, MPI_INT, l2->L2Gmap[0], 123, comm);
	/* recv data from local master */
	if (l2->G2Dmap == NULL) {
	    MPI_Type_indexed(nprocs, CNSTV rcnts, CNSTV rdsps, rtype, &type);
	}
	else {
	    cnts = phgAlloc(2 * nprocs * sizeof(*cnts));
	    dsps = cnts + nprocs;
	    for (i = 0; i < nprocs; i++) {
		j = l2->G2Dmap[i];
		cnts[j] = rcnts[i];
		dsps[j] = rdsps[i];
	    }
	    MPI_Type_indexed(nprocs, cnts, dsps, rtype, &type);
	    phgFree(cnts);
	    cnts = NULL;
	}
	MPI_Type_commit(&type);
	MPI_Recv(rbuff, 1, type, l2->L2Gmap[0], 123, comm, &status);
	MPI_Type_free(&type);
    }

    return MPI_SUCCESS;
}

int
MPI_Gather2(const void *sbuff, int scnt, MPI_Datatype stype,
	    void *rbuff, int rcnt, MPI_Datatype rtype, int root, MPI_Comm comm)
{
    int i, j, rank, nprocs, root_q;
    int *cnts = NULL, *dsps = NULL;
    size_t rsize, size = 0;
    MPI_Datatype type = MPI_DATATYPE_NULL;
    MPI_Status status;
    char *buff = NULL;
    L2DATA *l2;

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Gather.\n", __func__);
	return MPI_Gather(CNSTV sbuff, scnt, stype,
				rbuff, rcnt, rtype, root, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (rank != root) {
	rtype = stype;
	rcnt = scnt;
    }

    MPI_Type_size(rtype, &i);
    rsize = i;
    root_q = D2QMAP(l2, G2DMAP(l2, root));

    /* Step 1: local Gather */
    if (l2->Lrank == 0) {
	if (l2->Qrank == root_q) {
	    cnts = phgAlloc(2 * nprocs * sizeof(*cnts));
	    dsps = cnts + nprocs;
	}
	buff = phgAlloc(l2->Lsize * rcnt * rsize);
	MPI_Type_contiguous(rsize, MPI_BYTE, &type);
	MPI_Type_commit(&type);
    }
    MPI_Gather(CNSTV sbuff, scnt, stype, buff, rcnt, type, 0, l2->Lcomm);

    /* Step 2: global Gather or Gatherv */
    if (l2->Lrank == 0) {
	size = l2->Lsize * rcnt;
	assert(size <= INT_MAX);
	scnt = (int)size;
	if (l2->Qrank != root_q) {
	    if (l2->Lsizes_are_equal)
		MPI_Gather(buff, scnt, type, NULL, scnt, rtype,
			   root_q, l2->Qcomm);
	    else
		MPI_Gatherv(buff, scnt, type, NULL, NULL, NULL, rtype,
			    root_q, l2->Qcomm);
	    MPI_Type_free(&type);
	}
	else {
	    char *buff1;
	    if (l2->G2Dmap == NULL && rank == root)
		buff1 = NULL;
	    else
		buff1 = phgAlloc(nprocs * rsize * rcnt);
	    if (l2->Lsizes_are_equal) {
		MPI_Gather(buff, scnt, type,
			   l2->G2Dmap == NULL && rank == root ?
				rbuff : buff1,
			   scnt, type, root_q, l2->Qcomm);
	    }
	    else {
		size = 0;
		for (i = 0; i < l2->Qsize; i++) {
		    assert(LSIZES(l2, i) * (size_t)rcnt <= INT_MAX);
		    cnts[i] = LSIZES(l2, i) * rcnt;
		    assert(size <= INT_MAX);
		    dsps[i] = size;
		    size += cnts[i];
		}
		MPI_Gatherv(buff, scnt, type,
			    l2->G2Dmap == NULL && rank == root ?
				rbuff : buff1,
			    cnts, dsps, type, root_q, l2->Qcomm);
	    }
	    /* Send data to root */
	    assert(nprocs * (size_t)rcnt <= INT_MAX);
	    scnt = nprocs * rcnt;
	    stype = type;
	    if (rank == root) {
		if (l2->G2Dmap != NULL) {
		    for (i = 0; i < nprocs; i++) {
			j = l2->G2Dmap[i];
			cnts[j] = rcnt;
			assert(i * (size_t)rcnt <= INT_MAX);
			dsps[j] = i * rcnt;
		    }
		    MPI_Type_indexed(nprocs, cnts, dsps, rtype, &type);
		    MPI_Type_commit(&type);
		    MPI_Sendrecv(buff1, scnt, stype, 0, 123,
				 rbuff, 1, type, 0, 123,
				 MPI_COMM_SELF, &status);
		    MPI_Type_free(&type);
		}
	    }
	    else {
		MPI_Send(buff1, scnt, stype, root, 123, comm);
	    }
	    MPI_Type_free(&stype);
	    phgFree(buff1);
	}
	phgFree(cnts);
	phgFree(buff);
    }
    else if (rank == root) {
	/* recv data from local master */
	if (l2->G2Dmap == NULL) {
	    assert(nprocs * (size_t)rcnt <= INT_MAX);
	    MPI_Recv(rbuff, nprocs * rcnt, rtype, l2->L2Gmap[0], 123,
		     comm, &status);
	}
	else {
	    cnts = phgAlloc(2 * nprocs * sizeof(*cnts));
	    dsps = cnts + nprocs;
	    for (i = 0; i < nprocs; i++) {
		j = l2->G2Dmap[i];
		cnts[j] = rcnt;
		assert(i * (size_t)rcnt <= INT_MAX);
		dsps[j] = i * rcnt;
	    }
	    MPI_Type_indexed(nprocs, cnts, dsps, rtype, &type);
	    phgFree(cnts);
	    cnts = NULL;
	    MPI_Type_commit(&type);
	    MPI_Recv(rbuff, 1, type, l2->L2Gmap[0], 123, comm, &status);
	    MPI_Type_free(&type);
	}
    }

    return MPI_SUCCESS;
}

int
MPI_Scatterv2(
    const void *sbuff, const int *scnts, const int *sdsps, MPI_Datatype stype,
    void *rbuff, int rcnt, MPI_Datatype rtype, int root, MPI_Comm comm)
{
    int i, j, rank, nprocs, root_q;
    int *cnts = NULL, *dsps = NULL, *scnts1 = NULL, *sdsps1 = NULL;
    char *buff0 = NULL, *buff = NULL;
    size_t size, rsize, ssize;
    MPI_Datatype type;
    MPI_Status status;
    L2DATA *l2;

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Scatterv.\n", __func__);
	return MPI_Scatterv(CNSTV sbuff, CNSTV scnts, CNSTV sdsps, stype,
			    rbuff, rcnt, rtype, root, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (rank != root)
	stype = rtype;

    MPI_Type_size(stype, &i);
    ssize = i;
    MPI_Type_size(rtype, &i);
    rsize = i;
    root_q = D2QMAP(l2, G2DMAP(l2, root));

    /* allocate cnts and dsps */
    if (l2->Lrank == 0) {
	i = 0;
	if (l2->Qrank == root_q)
	    i = l2->Qsize;
	if (i < l2->Lsize)
	    i = l2->Lsize;
	cnts = phgAlloc(2 * i * sizeof(*cnts));
	dsps = cnts + i;
    }

    /* Step 1: rearrange data according to Dcomm and send to local master */
    if (l2->Qrank == root_q) {
	/* reorder and send data to local master */
	if (rank == root) {
	    if (l2->G2Dmap == NULL) {
		size = 0;
		for (i = 0; i < nprocs; i++)
		    size += scnts[i];
		scnts1 = (void *)scnts;
		sdsps1 = (void *)sdsps;
	    }
	    else {
		scnts1 = phgAlloc(2 * nprocs * sizeof(*scnts1));
		sdsps1 = scnts1 + nprocs;
		size = 0;
		for (i = 0; i < nprocs; i++) {
		    j = l2->G2Dmap[i];
		    scnts1[j] = scnts[i];
		    sdsps1[j] = sdsps[i];
		    size += scnts[i];
		}
	    }
	    MPI_Type_indexed(nprocs, scnts1, sdsps1, stype, &stype);
	    MPI_Type_commit(&stype);
	    if (l2->Lrank == 0) {	/* root is local master */
		buff0 = phgAlloc(size * ssize);
		MPI_Type_contiguous(ssize, MPI_BYTE, &type);
		MPI_Type_commit(&type);
		assert(size <= INT_MAX);
		MPI_Sendrecv(CNSTV sbuff, 1, stype, 0, 123,
			     buff0, (int)size, type, 0, 123,
			     MPI_COMM_SELF, &status);
	    }
	    else {			/* root is not local master */
		/* send scnts1/sdsps1 to local master */
		int c[] = {nprocs, nprocs, 1};
		MPI_Aint d[] = {0, (const char *)sdsps1 - (const char *)scnts1,
				(const char *)&i - (const char *)scnts1};
		i = (int)ssize;
#if 1	/* MPI-2 */
		MPI_Type_create_hindexed(3, c, d, MPI_INT, &type);
#else	/* MPI-1 */
		MPI_Type_hindexed(3, c, d, MPI_INT, &type);
#endif
		MPI_Type_commit(&type);
		MPI_Send(scnts1, 1, type, l2->L2Gmap[0], 123, comm);
		MPI_Type_free(&type);
		/* send data to local master */
		MPI_Send(CNSTV sbuff, 1, stype, l2->L2Gmap[0], 123, comm);
		if (scnts1 != scnts)
		    phgFree(scnts1);
	    }
	    MPI_Type_free(&stype);
	}
	else if (l2->Lrank == 0) {	/* local master but not root */
	    /* recv scnts1/sdsps1/ssize from root */
	    scnts1 = phgAlloc((2 * nprocs + 1) * sizeof(*scnts1));
	    sdsps1 = scnts1 + nprocs;
	    MPI_Recv(scnts1, 2 * nprocs + 1, MPI_INT, root, 123, comm, &status);
	    ssize = sdsps1[nprocs];
	    /* recv data from root */
	    size = 0;
	    for (i = 0; i < nprocs; i++)
		size += scnts1[i];
	    buff0 = phgAlloc(size * ssize);
	    MPI_Type_contiguous(ssize, MPI_BYTE, &type);
	    MPI_Type_commit(&type);
	    assert(size <= INT_MAX);
	    MPI_Recv(buff0, (int)size, type, root, 123, comm, &status);
	}
    }

    /* Step 2: Scatter data to local masters */
    if (l2->Lrank == 0) {
	if (l2->Qrank == root_q) {
	    memset(cnts, 0, l2->Qsize * sizeof(*cnts));
	    for (i = 0; i < nprocs; i++)
		cnts[D2QMAP(l2, i)] += scnts1[i];
	    size = 0;
	    for (i = 0; i < l2->Qsize; i++) {
		assert(size <= INT_MAX);
		dsps[i] = (int)size;
		size += cnts[i];
	    }
	    size = 0;
	    for (i = 0; i < l2->Lsize; i++)
		size += scnts1[l2->L2Gmap[i]];
	    buff = phgAlloc(size * ssize);
	    MPI_Scatterv(buff0, cnts, dsps, type, buff, (int)size, type,
			 root_q, l2->Qcomm);
	    phgFree(buff0);
	    size = 0;
	    for (i = 0; i < l2->Lsize; i++) {
		cnts[i] = scnts1[l2->L2Gmap[i]];
		assert(size <= INT_MAX);
		dsps[i] = (int)size;
		size += cnts[i];
	    }
	    if (scnts1 != scnts)
		phgFree(scnts1);
	}
	else {
	    /* Gather local rcnt */
	    MPI_Gather(&rcnt, 1, MPI_INT, cnts, 1, MPI_INT, 0, l2->Lcomm);
	    size = 0;
	    for (i = 0; i < l2->Lsize; i++) {
		assert(size <= INT_MAX);
		dsps[i] = (int)size;
		size += cnts[i];
	    }
	    buff = phgAlloc(size * rsize);
	    MPI_Type_contiguous(rsize, MPI_BYTE, &type);
	    MPI_Type_commit(&type);
	    assert(size <= INT_MAX);
	    MPI_Scatterv(NULL, NULL, NULL, stype, buff, (int)size, type,
			 root_q, l2->Qcomm);
	}
    }
    else if (l2->Qrank != root_q) {
	/* Gather rcnt to local master */
	MPI_Gather(&rcnt, 1, MPI_INT, NULL, 1, MPI_INT, 0, l2->Lcomm);
    }

    /* Step 3: local Scatterv */
    MPI_Scatterv(buff, cnts, dsps, type, rbuff, rcnt, rtype, 0, l2->Lcomm);

    if (l2->Lrank == 0) {
	phgFree(buff);
	phgFree(cnts);
	MPI_Type_free(&type);
    }

    return MPI_SUCCESS;
}

int
MPI_Scatter2(const void *sbuff, int scnt, MPI_Datatype stype,
	     void *rbuff, int rcnt, MPI_Datatype rtype, int root, MPI_Comm comm)
{
    int i, j, rank, nprocs, root_q;
    int *cnts = NULL, *dsps = NULL;
    char *buff0 = NULL, *buff = NULL;
    size_t size, rsize, ssize;
    MPI_Datatype type = MPI_DATATYPE_NULL, type0;
    MPI_Status status;
    L2DATA *l2;

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Scatter.\n", __func__);
	return MPI_Scatter(CNSTV sbuff, scnt, stype,
				 rbuff, rcnt, rtype, root, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    if (rank != root) {
	stype = rtype;
	scnt = rcnt;
    }

    MPI_Type_size(stype, &i);
    ssize = i;
    MPI_Type_size(rtype, &i);
    rsize = i;
    root_q = D2QMAP(l2, G2DMAP(l2, root));

    /* allocate cnts and dsps */
    if (l2->Lrank == 0) {
	i = 0;
	if (l2->Qrank == root_q && !l2->Lsizes_are_equal)
	    i = l2->Qsize;
	if (i < l2->Lsize)
	    i = l2->Lsize;
	cnts = phgAlloc(2 * i * sizeof(*cnts));
	dsps = cnts + i;
    }

    /* Step 1: rearrange data according to Dcomm and send to local master */
    if (l2->Qrank == root_q) {
	/* reorder and send data to local master */
	if (rank == root) {
	    if (l2->G2Dmap == NULL) {
		size = nprocs * (size_t)scnt;
		assert(size <= INT_MAX);
		if (l2->Lrank != 0) {
		    MPI_Type_contiguous((int)size, stype, &type0);
		    MPI_Type_commit(&type0);
		}
	    }
	    else {
		int *scnts, *sdsps;
		scnts = phgAlloc(2 * nprocs * sizeof(*scnts));
		sdsps = scnts + nprocs;
		size = 0;
		for (i = 0; i < nprocs; i++) {
		    j = l2->G2Dmap[i];
		    scnts[j] = scnt;
		    assert(size <= INT_MAX);
		    sdsps[j] = (int)size;
		    size += scnt;
		}
		MPI_Type_indexed(nprocs, scnts, sdsps, stype, &type0);
		MPI_Type_commit(&type0);
		phgFree(scnts);
	    }
	    if (l2->Lrank == 0) {	/* root is local master */
		MPI_Type_contiguous(ssize, MPI_BYTE, &type);
		MPI_Type_commit(&type);
		if (l2->G2Dmap != NULL) {
		    buff0 = phgAlloc(size * ssize);
		    assert(size <= INT_MAX);
		    MPI_Sendrecv(CNSTV sbuff, 1, type0, 0, 123,
				 buff0, (int)size, type, 0, 123,
				 MPI_COMM_SELF, &status);
		    MPI_Type_free(&type0);
		    type0 = type;
		}
		else {
		    buff0 = (void *)sbuff;
		    type0 = stype;
		}
	    }
	    else {			/* root is not local master */
		/* send data to local master */
		MPI_Send(CNSTV sbuff, 1, type0, l2->L2Gmap[0], 123, comm);
		MPI_Type_free(&type0);
	    }
	}
	else if (l2->Lrank == 0) {	/* local master but not root */
	    /* recv data from root */
	    size = nprocs * (size_t)rcnt;
	    buff0 = phgAlloc(size * rsize);
	    MPI_Type_contiguous(rsize, MPI_BYTE, &type);
	    MPI_Type_commit(&type);
	    assert(size <= INT_MAX);
	    MPI_Recv(buff0, (int)size, type, root, 123, comm, &status);
	    type0 = type;
	}
    }

    /* Step 2: Scatter data to local masters */
    if (l2->Lrank == 0) {
	if (l2->Lsizes_are_equal) {
	    size = l2->Lsize * rcnt;
	    buff = phgAlloc(size * rsize);
	    assert(size <= INT_MAX);
	    MPI_Scatter(buff0, (int)size, type0,
			buff, (int)size, rtype, root_q, l2->Qcomm);
	    if (buff0 != sbuff)
		phgFree(buff0);
	}
	else {
	    if (l2->Qrank == root_q) {
		size = 0;
		for (i = 0; i < l2->Qsize; i++) {
		    assert(size <= INT_MAX);
		    dsps[i] = (int)size;
		    assert(LSIZES(l2, i) * (size_t)rcnt <= INT_MAX);
		    cnts[i] = LSIZES(l2, i) * rcnt;
		    size += cnts[i];
		}
		size = l2->Lsize * rcnt;
		buff = phgAlloc(size * rsize);
		assert(size <= INT_MAX);
		MPI_Scatterv(buff0, cnts, dsps, type0, buff, (int)size, type,
			     root_q, l2->Qcomm);
		if (buff0 != sbuff)
		    phgFree(buff0);
	    }
	    else {
		size = l2->Lsize * rcnt;
		buff = phgAlloc(size * rsize);
		MPI_Type_contiguous(rsize, MPI_BYTE, &type);
		MPI_Type_commit(&type);
		assert(size <= INT_MAX);
		MPI_Scatterv(NULL, NULL, NULL, stype, buff, (int)size, type,
			 root_q, l2->Qcomm);
	    }
	}
	phgFree(cnts);
	if (type != MPI_DATATYPE_NULL)
	    MPI_Type_free(&type);
    }

    /* Step 3: local Scatter. FIXME: what if real_extent(rtype)!=rsize? */
    MPI_Scatter(buff, rcnt, rtype, rbuff, rcnt, rtype, 0, l2->Lcomm);

    if (l2->Lrank == 0)
	phgFree(buff);

    return MPI_SUCCESS;
}

int
MPI_Allreduce2(const void *sbuff, void *rbuff, int count, MPI_Datatype type,
	       MPI_Op op, MPI_Comm comm)
{
    char *buffer = NULL, *buffer1 = NULL;
    L2DATA *l2;

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Allreduce.\n", __func__);
	return MPI_Allreduce(CNSTV sbuff, rbuff, count, type, op, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    /* local Reduce */
    if (l2->Lrank == 0) {
	MPI_Aint lb, extent;
	MPI_Type_get_extent(type, &lb, &extent);
	buffer = phgAlloc(extent * count);
	buffer1 = buffer - lb;
    }
    MPI_Reduce(CNSTV sbuff, buffer1, count, type, op, 0, l2->Lcomm);

    /* global Allreduce */
    if (l2->Lrank == 0) {
	MPI_Allreduce(buffer1, rbuff, count, type, op, l2->Qcomm);
	phgFree(buffer);
    }

    /* local Bcast */
    MPI_Bcast(rbuff, count, type, 0, l2->Lcomm);

    return MPI_SUCCESS;
}

int
MPI_Reduce2(const void *sbuff, void *rbuff, int count, MPI_Datatype type,
	       MPI_Op op, int root, MPI_Comm comm)
{
    char *buffer = NULL;
    MPI_Aint lb = 0, extent = 0;
    int rank, root_q;
    L2DATA *l2;

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Reduce.\n", __func__);
	return MPI_Reduce(CNSTV sbuff, rbuff, count, type, op, root, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    /* local Reduce */
    if (l2->Lrank == 0) {
	MPI_Type_get_extent(type, &lb, &extent);
	buffer = phgAlloc(extent * count);
    }
    MPI_Reduce(CNSTV sbuff, buffer - lb, count, type, op, 0, l2->Lcomm);

    /* global Reduce */
    MPI_Comm_rank(comm, &rank);
    root_q = D2QMAP(l2, G2DMAP(l2, root));
    if (l2->Lrank == 0) {
	if (l2->Qrank != root_q) {
	    /* Non destination group, simply call Reduce */
	    MPI_Reduce(buffer - lb, rbuff, count, type, op, root_q,
		       l2->Qcomm);
	}
	else {
	    /* Destination group */
	    if (rank == root) {
		/* Proc 0 of Qcomm == 'root' */
	        MPI_Reduce(buffer - lb, rbuff, count, type, op, root_q,
			   l2->Qcomm);
	    }
	    else {
		char *buffer1 = phgAlloc(extent * count);
		/* put result in buffer1 */
	        MPI_Reduce(buffer - lb, buffer1 - lb, count, type, op, root_q,
			   l2->Qcomm);
		phgFree(buffer);
		buffer = buffer1;
		/* send result to 'root' */
		MPI_Send(buffer - lb, count, type, root, 123, comm);
	    }
	}
	phgFree(buffer);
    }

    if (rank == root && l2->Lrank != 0) {
	MPI_Status status;
	MPI_Recv(rbuff, count, type, l2->L2Gmap[0], 123, comm, &status);
    }

    return MPI_SUCCESS;
}

int
MPI_Bcast2(void *buffer, int count, MPI_Datatype type, int root, MPI_Comm comm)
{
    int rank, root_q;
    MPI_Status status;
    L2DATA *l2;

    if (!MPI_L2coll_flag || (l2 = MPI_L2check(comm)) == NULL) {
	if (phgVerbosity > 1)
	    phgPrintf("##### %s: using original MPI_Bcast.\n", __func__);
	return MPI_Bcast(buffer, count, type, root, comm);
    }

    if (phgVerbosity > 1)
	phgPrintf("##### %s: using two-level algorithm.\n", __func__);

    MPI_Comm_rank(comm, &rank);
    root_q = D2QMAP(l2, G2DMAP(l2, root));

    /* send data to local master */
    if (rank == root && l2->Lrank != 0)
	MPI_Send(buffer, count, type, l2->L2Gmap[0], 123, comm);
    else if (l2->Lrank == 0 && l2->Qrank == root_q && rank != root)
	MPI_Recv(buffer, count, type, root, 123, comm, &status);

    /* Bcast over Qcomm */
    if (l2->Lrank == 0)
	MPI_Bcast(buffer, count, type, root_q, l2->Qcomm);

    /* local Bcast */
    MPI_Bcast(buffer, count, type, 0, l2->Lcomm);

    return MPI_SUCCESS;
}

#endif	/* USE_MPI */
