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

/* $Id: mpi-2level.h,v 1.2 2014/09/15 07:31:12 zlb Exp $
 *
 * Header file for the 2-level collective communications functions
 * (for their usage see comments at the top of the file mpi-2level.c)
 *
 * */

#ifndef PHG_MPI_2LEVEL_H

#ifndef USE_MPI
# define USE_MPI 1
#endif

#if USE_MPI

#ifdef __cplusplus
extern "C" {
#endif

/* process reordering types */
enum {
    L2_NONE = 0,        /* no reordering */
    L2_ROUND_ROBIN = 1, /* round-robin */
    L2_BLOCK = 2        /* block-wise */
};

/* The data stored in the l2key attribute.
 *
 * Note: we are concerned with 4 communicators:
 *	comm:	The global communicator with the attribute
 *	Lcomm:	The local communicator (procs on a same node)
 *	Qcomm:	The quotient communicator (local masters, with Lrank==0)
 *	comm_d: The dummy communicator with L2_BLOCK ordering
 */
typedef struct { /* Size	Meaning			Defined on */
    /* tables */ /* ----	-------			---------- */
    int *D2Qmap; /* [nprocs],	comm_d -> Qcomm,	all processes */
    int *L2Gmap; /* [Lsize],	Lcomm -> comm,		all processes */
    int *G2Dmap; /* [nprocs],	comm -> comm_d,		all processes */
    int *Q2Dmap;/* [Qsize + 1],	comm_q -> comm_d,	local masters */
    /* variables */
    MPI_Comm Lcomm, Qcomm;	/* local and quotient comms */
    int Lrank, Lsize, Qrank, Qsize;
    int Lsizes_are_equal;	/* TRUE if all Lsize are equal */
    int ref_count;		/* reference count of the attribute */
} L2DATA;

/* exported variables */
#ifdef PHG_PHG_H
extern BOOLEAN MPI_L2coll_flag;
#else   /* PHG_PHG_H */
extern int MPI_L2coll_flag;
#endif  /* PHG_PHG_H */

#define G2DMAP(l2, i) ((l2)->G2Dmap==NULL ? (i) : (l2)->G2Dmap[i])
#define D2QMAP(l2, i) ((l2)->Lsizes_are_equal ? \
			(i) / (l2)->Lsize : (l2)->D2Qmap[i])
#define Q2DMAP(l2, i) ((l2)->Lsizes_are_equal ? \
			(i)*(l2)->Lsize : (l2)->Q2Dmap[i])
#define LSIZES(l2, i) ((l2)->Lsizes_are_equal ? \
			(l2)->Lsize : (l2)->Q2Dmap[(i) + 1] - (l2)->Q2Dmap[i])

/* user interface functions */
void	MPI_L2setup(MPI_Comm *comm, int reorder_type, int ppn);
L2DATA *MPI_L2check(MPI_Comm comm);

/* 2-level MPI collective communication functions */
int MPI_Alltoallv2(
    const void *sbuff, const int *scnts, const int *sdsps, MPI_Datatype stype,
    void *rbuff, const int *rcnts, const int *rdsps, MPI_Datatype rtype,
    MPI_Comm comm);
int MPI_Alltoall2(const void *sbuff, int scnt, MPI_Datatype stype,
		  void *rbuff, int rcnt, MPI_Datatype rtype, MPI_Comm comm);
int MPI_Allgatherv2(const void *sbuff, int count, MPI_Datatype stype,
		    void *rbuff, const int *rcnts, const int *rdsps,
		    MPI_Datatype rtype, MPI_Comm comm);
int MPI_Allgather2(const void *sbuff, int scnt, MPI_Datatype stype,
		   void *rbuff, int rcnt, MPI_Datatype rtype, MPI_Comm comm);
int MPI_Gatherv2(const void *sbuff, int scnt, MPI_Datatype stype,
		 void *rbuff, const int *rcnts, const int *rdsps,
		 MPI_Datatype rtype, int root, MPI_Comm comm);
int MPI_Gather2(const void *sbuff, int scnt, MPI_Datatype stype,
		void *rbuff, int rcnt, MPI_Datatype rtype,
		int root, MPI_Comm comm);
int MPI_Scatterv2(const void *sbuff, const int *scnts, const int *sdsps,
		  MPI_Datatype stype, void *rbuff, int rcnt, MPI_Datatype rtype,
		  int root, MPI_Comm comm);
int MPI_Scatter2(const void *sbuff, int scnt, MPI_Datatype stype,
		 void *rbuff, int rcnt, MPI_Datatype rtype,
		 int root, MPI_Comm comm);
int MPI_Allreduce2(const void *sbuff, void *rbuff, int count, MPI_Datatype type,
		   MPI_Op op, MPI_Comm comm);
int MPI_Reduce2(const void *sbuff, void *rbuff, int count, MPI_Datatype type,
		MPI_Op op, int root, MPI_Comm comm);
int MPI_Bcast2(void *buffer, int count, MPI_Datatype type, int root,
	       MPI_Comm comm);

#if !defined(DONT_RENAME_MPI_FUNCTIONS)
/*----------------------------------------------------------------------*/
/* Intercept MPI collective functions by redefining them as macros */
/*----------------------------------------------------------------------*/

/*---------- Alltoallv ----------*/
# ifdef MPI_Alltoallv
#  undef MPI_Alltoallv
# endif
# define MPI_Alltoallv	MPI_Alltoallv2

# ifdef MPI_AlltoallvBegin
#  undef MPI_AlltoallvBegin
#  undef MPI_AlltoallvEnd
# endif
# define MPI_AlltoallvBegin(sbuff, scnts, sdsps, stype, rbuff, rcnts, rdsps, \
			   rtype, comm, reqs, nsends, nrecvs) \
    MPI_Alltoallv2(sbuff, scnts, sdsps, stype, rbuff, rcnts, rdsps, rtype, comm)
# define MPI_AlltoallvEnd(nsends, nrecvs, reqs)

/*---------- Alltoall ----------*/
# ifdef MPI_Alltoall
#  undef MPI_Alltoall
# endif
# define MPI_Alltoall	MPI_Alltoall2

/*---------- Allgatherv ----------*/
# ifdef MPI_Allgatherv
#  undef MPI_Allgatherv
# endif
# define MPI_Allgatherv	MPI_Allgatherv2

/*---------- Allgather ----------*/
# ifdef MPI_Allgather
#  undef MPI_Allgather
# endif
# define MPI_Allgather	MPI_Allgather2

/*---------- Gatherv ----------*/
# ifdef MPI_Gatherv
#  undef MPI_Gatherv
# endif
# define MPI_Gatherv	MPI_Gatherv2

/*---------- Gather ----------*/
# ifdef MPI_Gather
#  undef MPI_Gather
# endif
# define MPI_Gather	MPI_Gather2

/*---------- Scatterv ----------*/
# ifdef MPI_Scatterv
#  undef MPI_Scatterv
# endif
# define MPI_Scatterv	MPI_Scatterv2

/*---------- Scatter ----------*/
# ifdef MPI_Scatter
#  undef MPI_Scatter
# endif
# define MPI_Scatter	MPI_Scatter2

/*---------- Allreduce ----------*/
# ifdef MPI_Allreduce
#  undef MPI_Allreduce
# endif
# define MPI_Allreduce	MPI_Allreduce2

/*---------- Reduce ----------*/
# ifdef MPI_Reduce
#  undef MPI_Reduce
# endif
# define MPI_Reduce	MPI_Reduce2

/*---------- Bcast ----------*/
# ifdef MPI_Bcast
#  undef MPI_Bcast
# endif
# define MPI_Bcast	MPI_Bcast2

/*----------------------------------------------------------------------*/
#endif	/* !defined(DONT_RENAME_MPI_FUNCTIONS) */

#ifdef __cplusplus
}
#endif

#endif	/* USE_MPI */

#define PHG_MPI_2LEVEL_H
#endif
