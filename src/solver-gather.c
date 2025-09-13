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

/* $Id: solver-gather.c,v 1.56 2022/09/21 02:35:28 zlb Exp $
 *
 * This solver splits a linear system into independent blocks, these blocks
 * are then scheduled and gathered to a subset of processes and are solved
 * in parallel by calling a sub-solver. */

#include "phg.h"
#include "phg/partition-utils.h"

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
    SOLVER	*solver;	/* solver for the gathered system */
    INT		*index;		/* map: original => gathered */
    INT		*order;		/* ordering (new indices increasing) */
    FLOAT	*idiag;		/* inversed diagonal */
    MAT		*mat;		/* saved columns (diagonal rows) */

    int		*scnts, *sdsps, *rcnts, *rdsps;
    BOOLEAN	factored;

    /* options variables */
    char	*solver_opts;
    INT		max_procs, min_size;
    int		solver_id;
    int		scheduler_id;
    int		remap_id;
} OEM_DATA;

#define	_t	((OEM_DATA *)solver->oem_data)

static int solver_id = 0;
static char *solver_opts = NULL;

static const char *scheduler_names[] = {
	"round-robin",	/* round-robin */
	"one-block",	/* solve on one process as one block */
	NULL};
enum {SCHED_ROUND_ROBIN = 0, SCHED_ROOT_ONLY = 1};
static int scheduler_id = SCHED_ROUND_ROBIN;

static const char *remap_names[] = {
	"none",
	"vector",	/* optimize data move for vector elements */
	"matrix",	/* optimize data move for matrix elements */
	NULL};
enum {REMAP_NONE = 0, REMAP_VECTOR = 1, REMAP_MATRIX = 2};
static int remap_id = REMAP_VECTOR;

static INT max_procs = 1;	/* max # of procs for a single system */
static INT min_size = 1000;	/* min size / proc */

static int
RegisterOptions(void)
{
    /* set default solver to the first direct solver (if available) */
    for (solver_id = 0;
	 phgSolverNames[solver_id] != NULL
	    && (phgSolverList[solver_id] == SOLVER_GATHER ||
		phgSolverList[solver_id] == SOLVER_DIAGONAL ||
		phgSolverList[solver_id]->iterative == TRUE);
	 solver_id++);
    if (phgSolverNames[solver_id] == NULL)
	 solver_id = 0;	/* fall back to the first solver */

    phgOptionsRegisterTitle("\nThe Gather solver options:", "\n", "gather");
    phgOptionsRegisterKeyword("-gather_solver", "The sub-solver used to "
					"solve the gathered blocks",
					phgSolverNames, &solver_id);
    phgOptionsRegisterString("-gather_solver_opts",
					"Options to pass to the sub-solver",
					&solver_opts);
    phgOptionsRegisterKeyword("-gather_solver_scheduler", "Method to schedule "
					"blocks to processes",
					scheduler_names, &scheduler_id);
    phgOptionsRegisterKeyword("-gather_solver_remap", "Controls remapping "
					"of blocks to processes",
					remap_names, &remap_id);
    phgOptionsRegisterInt("-gather_solver_max_procs", "Max number of procs "
					"per block", &max_procs);
    phgOptionsRegisterInt("-gather_solver_min_size", "Min size per proc",
					&min_size);

    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));
    _t->solver_id = solver_id;
    _t->scheduler_id = scheduler_id;
    _t->remap_id = remap_id;
    _t->max_procs = max_procs;
    _t->min_size = min_size;
    if (solver_opts != NULL)
	_t->solver_opts = strdup(solver_opts);

    return 0;
}

static int
Create(SOLVER *solver)
{
    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    phgFree(_t->solver_opts);
    phgFree(_t->order);
    phgFree(_t->index);
    phgFree(_t->idiag);
    phgFree(_t->scnts);
    phgFree(_t->rcnts);
    phgFree(_t->sdsps);
    phgFree(_t->rdsps);
    phgMatDestroy(&_t->mat);
    phgSolverDestroy(&_t->solver);

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    int rank, nprocs, pass0 = 0, pass1= 0;
    INT nlocal, localsize, i, j, k, n, N0;
    MAP *map;
    MAT *B;
    const MAT_ROW *row;
    MAT_ROW *row1;
    FLOAT *data0, *data1;
    INT alloc, ncols0, *cols0, ncols1, *cols1;
    double time0 = 0., time1;

    INT *counts, ntotal, n1;
    struct {
	INT	block;			/* block number */
	INT	index;			/* global index */
    } *list, *list0, *list1;
    size_t *blocksize = NULL;
    int ii, flag;
#if USE_MPI
    MPI_Comm comm;
#endif	/* USE_MPI */

   /* Note: currently only ordinary matrices allowed since it's easier and
    *	    more efficient to work with local column indices. */
    if (solver->mat->type == PHG_PACKED)
	phgMatUnpack(solver->mat);
    assert(solver->mat->type == PHG_UNPACKED);

    if (_t->factored)
	return 0;

    _t->factored = TRUE;

    if (solver->sent_to_oem == TRUE)
	phgError(1, "%s:%d, data already sent to OEM solver!\n",
						__FILE__, __LINE__);
    solver->sent_to_oem = TRUE;

    rank = solver->mat->rmap->rank;
    nprocs = solver->mat->rmap->nprocs;
    nlocal = solver->mat->rmap->nlocal;
    N0 = solver->mat->rmap->partition[rank];

    _t->scnts = phgAlloc(nprocs * sizeof(*_t->scnts));
    _t->sdsps = phgAlloc(nprocs * sizeof(*_t->sdsps));
    _t->rcnts = phgAlloc(nprocs * sizeof(*_t->rcnts));
    _t->rdsps = phgAlloc(nprocs * sizeof(*_t->rdsps));

#if USE_MPI
    comm = solver->mat->rmap->comm;
#endif	/* USE_MPI */

    /*--------------------- determine diagonal blocks --------------------*/

    time0 = phgGetTime(NULL);

    list = phgAlloc((localsize = solver->mat->localsize) * sizeof(*list));

#if USE_MPI
    if (nprocs > 1 && solver->mat->cinfo == NULL)
	    solver->mat->cinfo = phgMapCreateCommInfo(
					solver->mat->cmap->comm,
					nprocs,
					rank,
					nlocal,
					localsize,
					solver->mat->O2Gmap,
					solver->mat->ordering,
					solver->mat->cmap->partition);
#endif	/* USE_MPI */

    if (_t->scheduler_id == SCHED_ROOT_ONLY) {
	/* mark diagonal rows with block = -1 */
	k = 0;
	j = nlocal;
	for (i = 0; i < nlocal; i++) {
	    row = solver->mat->rows + i;
	    if (row->ncols > 1 || (row->ncols == 1 && row->cols[0] != i)) {
		/* non diagonal row */
		j--;
		list[j].block = 0;
		list[j].index = i;
	    }
	    else {
		/* diagonal row */
		list[k].block = -1;
		list[k].index = i;
		k++;
	    }
	}
#if USE_MPI
	MPI_Allreduce(&k, &n1, 1, PHG_MPI_INT, MPI_SUM,
				solver->mat->rmap->comm);
#else	/* USE_MPI */
	n1 = k;
#endif	/* USE_MPI */

	n = ntotal = (solver->mat->rmap->nglobal > n1 ? 1 : 0);
	if (n > 0) {
	    list0 = phgAlloc(sizeof(*list0));
	    list0[0].block = 0;
	    list0[0].index = solver->mat->rmap->nglobal - n1;
	    if (_t->remap_id != REMAP_NONE) {
		blocksize = phgAlloc(sizeof(*blocksize));
		if (_t->remap_id == REMAP_VECTOR)
		    blocksize[0] = solver->mat->rmap->nlocal;
		else
		    blocksize[0] = solver->mat->nnz_d + solver->mat->nnz_o;
	    }
	}
	else {
	    list0 = NULL;
	}

	goto label1;
    }

    /* The algorithm: 
     *	    while changes; do
     *		scan locally
     *		exchange data
     *	    done
     */

    /* initialize list[] */
    for (i = 0; i < localsize; i++) {
	list[i].index = list[i].block =
#if USE_MPI
		(i < nlocal ? i + N0 : solver->mat->O2Gmap[i - nlocal]);
#else	/* USE_MPI */
		i + N0;
#endif	/* USE_MPI */
    }

    cols0 = cols1 = NULL;
#if USE_MPI
    if (nprocs > 1 && solver->mat->cinfo != NULL) {
	/* Note: size of list0 is rsize, size of list1 is ssize */
	cols0 = phgAlloc((solver->mat->cinfo->ssize + solver->mat->cinfo->rsize)
			* sizeof(*cols0));
	cols1 = cols0 + solver->mat->cinfo->rsize;
    }
#endif	/* USE_MPI */

    flag = 1;
    pass0 = 0;
    while (TRUE) {
	pass0++;
	pass1 = 0;
	while (flag) {
	    pass1++;
	    flag = 0;
	    for (i = 0; i < nlocal; i++) {
		/* Not using phgMatGetRow() here since it returns GLOBAL rows */
		row = solver->mat->rows + i;
		/* find smallest index in the row */
		n = list[i].block;
		for (j = 0; j < row->ncols; j++) {
		    if ((k = list[row->cols[j]].block) < n)
			n = k;
		}
		list[i].block = n;
		for (j = 0; j < row->ncols; j++) {
		    if (list[row->cols[j]].block > n) {
			list[row->cols[j]].block = n;
			flag = 1;
		    }
		}
	    }
	}
	phgInfo(2, "\tlocal passes: %d\n", pass1);

#if USE_MPI
	if (nprocs > 1 && solver->mat->cinfo != NULL) {
	    for (i = 0; i < localsize - nlocal; i++) {
		cols0[i] = list[nlocal + solver->mat->cinfo->rind[i]].block;
	    }
	    MPI_Alltoallv(cols0, solver->mat->cinfo->rcnts,
				solver->mat->cinfo->rdsps, PHG_MPI_INT, 
			  cols1, solver->mat->cinfo->scnts,
				solver->mat->cinfo->sdsps, PHG_MPI_INT, comm);
	    flag = 0;
	    for (i = 0; i < solver->mat->cinfo->ssize; i++) {
		j = solver->mat->cinfo->sind[i];
		assert(j >= 0 && j < nlocal);
		if (list[j].block > cols1[i]) {
		    list[j].block = cols1[i];
		    flag = 1;
		}
	    }
	    MPI_Allreduce(&flag, &ii, 1, MPI_INT, MPI_MAX, comm);
	    if (ii == 0)
		break;
	}
	else
#endif	/* USE_MPI */
	    break;
    }
    if (nprocs == 1)
	pass0 = pass1;

    phgFree(cols0);

    if (nlocal > 0) {
	/* mark diagonal rows with block = -1 */
	for (i = 0; i < nlocal; i++) {
	    row = solver->mat->rows + i;
	    if (row->ncols > 1)
		continue;
	    if (row->ncols == 1 && row->cols[0] != i)
		continue;
	    list[i].block = -1;
	}
	qsort(list, nlocal, sizeof(*list), phgCompINT);
    }

    /* collect group indices, skipping diagonal rows (blocks of size 1) */
    i = n = k = 0;
    while (i < nlocal) {
	if (list[i].block == -1) {
	    /* diagonal row */
	    k++;	/* counter for diagonal rows */
	    i++;
	    continue;
	}
	j = i;
	while (++j < nlocal && list[j].block == list[i].block);
	n++;		/* counter for blocks of size > 1 */
	i = j;
    }
    if (k > 0) {
	/* the first entry is used to save number of diagonal rows */
	n++;
    }
    list0 = phgAlloc(n * sizeof(*list0));
    i = n = 0;
    if (k > 0) {
	/* save number of diagonal rows */
	list0[n].block = -1;
	list0[n++].index = k;
    }
    while (i < nlocal) {
	if (list[i].block == -1) {
	    i++;
	    continue;
	}
	j = i;
	while (++j < nlocal && list[j].block == list[i].block);
	list0[n].block = list[i].block;
	list0[n++].index = j - i;
	i = j;
    }

#if USE_MPI
    if (nprocs > 1) {
	MPI_Datatype type;
	MPI_Type_contiguous(sizeof(*list), MPI_BYTE, &type);
	MPI_Type_commit(&type);
	ii = n;
	MPI_Gather(&ii, 1, PHG_MPI_INT, _t->scnts, 1, PHG_MPI_INT, 0, comm);
	list1 = NULL;		/* avoid a gcc warning */
	if (rank == 0) {
	    j = 0;
	    for (ii = 0; ii < nprocs; ii++) {
		_t->sdsps[ii] = j;
		j += _t->scnts[ii];
	    }
	    list1 = phgAlloc(j * sizeof(*list1));
	}
	MPI_Gatherv(list0, n, type, list1, _t->scnts, _t->sdsps, type, 0, comm);
	if (rank == 0) {
	    n = j;
	    phgFree(list0);
	    list0 = list1;
	    if (n > 0)
		qsort(list0, n, sizeof(*list0), phgCompINT);
	    i = k = 0;
	    while (i < n) {
		j = i;
		list0[k].block = list0[i].block;
		list0[k].index = list0[i].index;
		while (++j < n && list0[j].block == list0[i].block)
		    list0[k].index += list0[j].index;
		/* Note: blocks from only one process, i.e. j - i == 1, are
		 * private to that process. */
		k++;
		i = j;
	    }
	    n = k;
	}
	MPI_Bcast(&n, 1, PHG_MPI_INT, 0, comm);
	if (rank != 0) {
	    phgFree(list0);
	    list0 = phgAlloc(n * sizeof(*list0));
	}
	MPI_Bcast(list0, n, type, 0, comm);

	MPI_Type_free(&type);
    }
#endif	/* USE_MPI */

    /* Now list0[i].block contains the block index of i-th block, and
     * list0[i].index contains the global size of the i-th block. */
    ntotal = n1 = 0;
    k = 0;
    for (i = 0; i < n; i++) {
	k += list0[i].index;
	if (list0[i].block == -1)
	    n1 += list0[i].index;
	else
	    ntotal++;
    }

    assert(n1 == 0 || list0[0].block == -1);
    assert(k == solver->mat->rmap->nglobal);
    assert(n == ntotal + (n1 > 0 ? 1 : 0));

    /* discard list0[0] if it corresponds to diagonal rows */
    if (list0[0].block == -1) {
	assert(n == ntotal + 1);
	--n;
	for (i = 0; i < n; i++)
	    list0[i] = list0[i + 1];
    }

    /* Change list[].block to index in list0[] */
    if (_t->remap_id != REMAP_NONE)
	blocksize = phgCalloc(n, sizeof(*blocksize));
    k = -1;
    i = 0;
    while (i < nlocal) {
	INT a;
	if ((j = list[i].block) == -1) {
	    list[i].block = -1;
	    list[i].index -= N0;
	    assert(solver->mat->rows[list[i].index].ncols <= 1);
	    i++;
	    continue;
	}
	while (++k < n) {
	    if (j == list0[k].block)
		break;
	}
	assert(k < n);
	j = i;
	if (_t->remap_id == REMAP_MATRIX) {
	    /* optimize wrt matrix: blocksize[k] = # nonzeros in this block */
	    blocksize[k] = (solver->mat->rows + list[i].index - N0)->ncols;
	    while (++j < nlocal && list[j].block == list[i].block)
		blocksize[k] += (solver->mat->rows + list[j].index - N0)->ncols;
	}
	else {
	    while (++j < nlocal && list[j].block == list[i].block);
	    if (_t->remap_id == REMAP_VECTOR) {
		/* optimize wrt vector: blocksize[k] = vector components */
		blocksize[k] = j - i;
	    }
	}
	a = k;
	while (i < j) {
	    list[i].block = a;
	    list[i].index -= N0;
	    i++;
	}
    }

label1:
    /* At this stage, the block information is described by 'list':
     * 	    list[i].index is the LOCAL row index.
     *	    Row list[i].index belongs to block list[i].block.
     *	    If list[i].block = -1 then it's a diagonal row.
     * 	    If list[i].block != -1 then list[i].block is the index in list0[]
     *
     * 	    list0[k].index is the global size of block k.
     * 	    list0[k].block is not used.
     *
     *	    n1 is the global number of diagonal rows (blocks of size 1),
     *	    ntotal is the global number of blocks of size > 1
     *
     * 	    if (_t->remap != REMAP_NONE) then blocksize[k] is the local size of
     * 	    block k (k=0,...,ntotal-1).
     */

    if (solver->monitor) {
	time1 = phgGetTime(NULL);
	phgPrintf("Gather solver statistics:\n");
	phgPrintf("    blocks of size = 1: %d (incl. %d Dirichlet DOF)\n",
			n1, solver->mat->rmap->bdry_nglobal);
	phgPrintf("    blocks of size > 1: %d\n", ntotal);
	phgPrintf("    global passes for finding independent blocks: %d\n",
			pass0);
	phgPrintf("    time for finding independent blocks: %0.4lf\n",
			time1 - time0);
	time0 = time1;
    }

    /*------------------------- Construct sub-solvers ----------------------*/

    /* Assign blocks to processes: block k ==> process list0[k].block */

    if (ntotal > 0) {
	/* Copy sizes of blocks to list1[] */
	list1 = phgAlloc(ntotal * sizeof(*list1));
	for (k = 0; k < ntotal; k++) {
	    list1[k].block = -list0[k].index;	/* negated block size */
	    list1[k].index = k;			/* index to list0[] */
	}

	qsort(list1, ntotal, sizeof(*list1), phgCompINT);

	for (k = 0; k < ntotal; k++) {
	    /* Note: the blocks are scheduled in reverse order using combined
	     * forward/backward round-robin as follows:
	     * 		0, 1, ..., nprocs - 1,
	     * 		nprocs - 1, ..., 1, 0,
	     * 		0, 1, ..., nprocs - 1,
	     * 		nprocs - 1, ..., 1, 0,
	     * 		... ...
	     */
	    ii = k % nprocs;
	    if (((k / nprocs) % 2) != 0)
		ii = nprocs - 1 - ii;
	    list1[ntotal - 1 - k].block = ii;
	}

	/* Make use of the algorithm in phgPartitionRemap0 to
	 * minimize inter-process data moves.
	 *	remap_id == REMAP_MATRIX: minimize data move for the matrix
	 *	remap_id == REMAP_VECTOR: minimize data move for the vector
	 */
	if (blocksize != NULL) {
	    double *datasize = phgCalloc(nprocs, sizeof(*datasize));
	    int *perm = phgCalloc(nprocs, sizeof(*perm));

	    if (solver->monitor) {
		phgInfo(1, "before calling phgPartitionRemap0:\n");
		for (k = 0; k < ntotal; k++)
		    phgInfo(1, "    block %d, datasize %lld ==> proc %d\n",
				list1[k].index,
				(long long)blocksize[list1[k].index],
				list1[k].block);
	    }

	    for (k = 0; k < ntotal; k++)
		datasize[list1[k].block] += blocksize[list1[k].index];
#if USE_MPI
	    phgPartitionRemap0(comm, datasize, perm, comm);
	    for (k = 0; k < ntotal; k++)
		list1[k].block = perm[list1[k].block];
#endif	/* USE_MPI */

	    if (solver->monitor) {
		phgInfo(1, "after calling phgPartitionRemap0:\n");
		for (k = 0; k < ntotal; k++)
		    phgInfo(1, "    block %d, datasize %lld ==> proc %d\n",
				list1[k].index,
				(long long)blocksize[list1[k].index],
				list1[k].block);
	    }

	    phgFree(blocksize);
	    phgFree(datasize);
	    phgFree(perm);
	}

	for (k = 0; k < ntotal; k++)
	    list0[list1[k].index].block = list1[k].block;

	phgFree(list1);
    }

    counts = phgCalloc(nprocs, sizeof(*counts));
    for (k = 0; k < ntotal; k++) {
	counts[list0[k].block] += list0[k].index;
    }

    _t->order = phgAlloc(nlocal * sizeof(*_t->order));
    _t->index = phgAlloc(localsize * sizeof(*_t->index));
    _t->idiag = phgAlloc(localsize * sizeof(*_t->idiag));

#if DEBUG
    k = 0;
    for (ii = 0; ii < nprocs; ii++)
	k += counts[ii];
    assert(k == solver->mat->rmap->nglobal - n1);
#endif	/* DEBUG */

#if USE_MPI
    map = phgMapCreateSimpleMap(comm, counts[rank],
			solver->mat->rmap->nglobal - n1);
#else	/* USE_MPI */
    map = phgMapCreateSimpleMap(0, counts[rank],
			solver->mat->rmap->nglobal - n1);
#endif	/* USE_MPI */

    /* Count # of knknowns to send to each process. */
    memset(counts, 0, nprocs * sizeof(*counts));
    for (i = 0; i < nlocal; i++) {
	k = list[i].block;			/* block number */
	if (k == -1)
	    continue;
	/* Note: block k is mapped to process 'list0[k].block' */
	counts[list0[k].block]++;
    }

    for (ii = 0; ii < nprocs; ii++)
	_t->scnts[ii] = counts[ii];

    /* compute start unknown index for each derived system:
     * 		counts[ii] is the first index for the derived system 'ii' */
#if USE_MPI
    {
	INT *tmp = phgAlloc(nprocs * sizeof(*tmp));
	MPI_Scan(counts, tmp, nprocs, PHG_MPI_INT, MPI_SUM, comm);
	for (ii = 0; ii < nprocs; ii++)
	    counts[ii] = tmp[ii] - counts[ii] + map->partition[ii];
	phgFree(tmp);
    }
#else
    counts[0] = 0;
#endif	/* USE_MPI */

    for (i = 0; i < nlocal; i++) {
	j = _t->order[i] = list[i].index;	/* row index */
	k = list[i].block;			/* block number */
	if (k == -1) {
	    /* diagonal row */
	    row = solver->mat->rows + j;
	    assert(row->ncols <= 1);
	    if (row->ncols == 0 || row->data[0] == 0.)
		phgError(1, "(%s:%d) row %d is empty, abort.\n",
				i + N0, __FILE__, __LINE__);
	    assert(row->cols[0] == j);
	    /* save diagonal entry and free the matrix row */
	    _t->idiag[j] = 1. / row->data[0];
	    _t->index[j] = -1;
	    if (solver->mat->blocks == NULL && solver->mat->refcount == 0) {
		row1 = solver->mat->rows + j;
		phgFree(row1->cols);
		phgFree(row1->data);
		row1->cols = NULL;
		row1->data = NULL;
		row1->ncols = row1->alloc = 0;
	    }
	}
	else {
	    /* compute the global index in map */
	    /* Note: block k is mapped to process 'list0[k].block' */
	    _t->index[j] = counts[list0[k].block]++;
	    _t->idiag[j] = 0.0;
	}
    }

    phgFree(counts);

#if USE_MPI
    /* get new indices for off-process entries */
    if (nprocs > 1 && solver->mat->cinfo != NULL) {
	cols0 = phgAlloc((solver->mat->cinfo->rsize +
			  solver->mat->cinfo->ssize) * sizeof(*cols0));
	cols1 = cols0 + solver->mat->cinfo->ssize;
	for (i = 0; i < solver->mat->cinfo->ssize; i++)
	    cols0[i] = _t->index[solver->mat->cinfo->sind[i]];
	MPI_Alltoallv(cols0, solver->mat->cinfo->scnts,
				solver->mat->cinfo->sdsps, PHG_MPI_INT, 
		      cols1, solver->mat->cinfo->rcnts,
				solver->mat->cinfo->rdsps, PHG_MPI_INT, comm);
	for (i = 0; i < localsize - nlocal; i++)
	    _t->index[nlocal + solver->mat->ordering[i]] = cols1[i];
	phgFree(cols0);
    }

    MPI_Alltoall(_t->scnts, 1, MPI_INT, _t->rcnts, 1, MPI_INT, comm);
    k = j = 0;
    for (ii = 0; ii < nprocs; ii++) {
	_t->sdsps[ii] = k;
	_t->rdsps[ii] = j;
	k += _t->scnts[ii];
	j += _t->rcnts[ii];
    }
#else	/* USE_MPI */
    _t->rcnts[0] = _t->scnts[0];
    _t->sdsps[0] = _t->rdsps[0] = 0;
#endif	/* USE_MPI */

    phgFree(list0);
    phgFree(list);

    if (map->nglobal == 0) {
	return 0;
    }

    /* Matrix containing saved columns */
    _t->mat = phgMapCreateMat(solver->mat->rmap, solver->mat->cmap);
    _t->mat->handle_bdry_eqns = FALSE;

    /* Temporary matrix for sending rows to destination processes */
    B = phgMapCreateMat(map, map);
    B->handle_bdry_eqns = FALSE;
    phgMapDestroy(&map);

    /* send rows to B and _t->mat */
    data0 = data1 = NULL;
    cols0 = cols1 = NULL;
    alloc = 0;
    for (i = 0; i < nlocal; i++) {
	if (_t->index[i] < 0)
	    continue;
	row = solver->mat->rows + i;
	if (alloc < row->ncols) {
	    alloc = row->ncols;
	    phgFree(data0);
	    data0 = phgAlloc(2 * alloc * sizeof(*data0));
	    data1 = data0 + alloc;
	    phgFree(cols0);
	    cols0 = phgAlloc(2 * alloc * sizeof(*cols0));
	    cols1 = cols0 + alloc;
	}
	ncols0 = ncols1 = 0;
	for (j = 0; j < row->ncols; j++) {
	    k = row->cols[j];
	    if (_t->index[k] >= 0) {
		data0[ncols0] = row->data[j];
		cols0[ncols0] = _t->index[k];
		ncols0++;
	    }
	    else {
		data1[ncols1] = row->data[j];
#if USE_MPI
		if (k >= nlocal)
		    cols1[ncols1] = solver->mat->O2Gmap[k - nlocal];
		else
#endif	/* USE_MPI */
		    cols1[ncols1] = k + N0;
		ncols1++;
	    }
	}
	if (solver->mat->blocks == NULL && solver->mat->refcount == 0) {
	    /* free the matrix row */
	    row1 = solver->mat->rows + i;
	    phgFree(row1->cols);
	    phgFree(row1->data);
	    row1->cols = NULL;
	    row1->data = NULL;
	    row1->ncols = row1->alloc = 0;
	}
	j = _t->index[i];
	phgMatAddGlobalEntries(B, 1, &j, ncols0, cols0, data0);
	phgMatAddLGEntries(_t->mat, 1, &i, ncols1, cols1, data1);
    }
    phgFree(data0);
    phgFree(cols0);
    phgMatAssemble(B);
    phgMatAssemble(_t->mat);

    if (solver->monitor) {
	time1 = phgGetTime(NULL);
	phgPrintf("    time for setting up derived system: %0.4lg\n",
			time1 - time0);
	time0 = time1;
    }

    n = B->rmap->nlocal;
    assert(B->localsize == n);	/* there should be no off-process entries */
    if (n == 0) {
	phgMatDestroy(&B);
	return 0;
    }

    /* Hacky: replace B->rmap and B->cmap with the serial map */
    map = phgMapCreateSimpleMap(MPI_COMM_SELF, n, n);
    phgMapDestroy(&B->rmap);
    phgMapDestroy(&B->cmap);
    B->rmap = map;
    B->cmap = phgMatGetRowMap(B);	/* this increments B->cmap->refcount */
#if USE_MPI
    if (B->cinfo != NULL)
	phgMapDestroyCommInfo(&B->cinfo);
#endif	/* USE_MPI */

#if 0
char fn[124];
sprintf(fn, "B%d.m", rank);
phgMatDumpMATLAB(B, "B", fn);
#endif

    /* Create the serial solver.
     *
     * Note: it's better to create the solver in 'Init()' to honor all solver
     * options for the serial solver, but we couldn't do it there since we
     * need actual matrix data to determine the serial map. */

    phgOptionsPush();
    phgSolverSetDefaultSuboptions();
    /* set options for the serial solver */
    phgOptionsSetKeyword("-solver", phgSolverNames[_t->solver_id]);
    phgOptionsSetOptions(_t->solver_opts);
    assert(_t->solver == NULL);
    _t->solver = phgMat2Solver(SOLVER_DEFAULT, B);
    phgOptionsPop();
    if (_t->solver->oem_solver == solver->oem_solver)
	phgError(1, "(%s:%d) solver regression detected, abort.\n",
							__FILE__, __LINE__);
    phgMatDestroy(&B);

    phgSolverAssemble(_t->solver);

    if (solver->monitor) {
	time1 = phgGetTime(NULL);
	phgInfo(-1, "solver = %s, unknowns = %d, assembly time = %0.4lg\n",
			_t->solver->oem_solver->name,
			_t->solver->mat->rmap->nglobal,
			time1 - time0);
	time0 = time1;
    }

    return 0;
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    INT i, j, k, n;
    VEC *v;
    FLOAT *buffer0, *buffer1;
    double time0 = 0.;

    assert(x->nvec == 1);
    Assemble(solver);

    if (solver->monitor)
	time0 = phgGetTime(NULL);

    n = (_t->solver == NULL ? 0 : _t->solver->mat->rmap->nglobal);

    /* copy initial solution and updated RHS to buffer0 */
    buffer0 = phgAlloc(2 * solver->mat->rmap->nlocal * sizeof(*buffer0));
    for (i = k = 0; i < solver->mat->rmap->nlocal; i++) {
	j = _t->order[i];
	if (_t->index[j] >= 0) {
	    buffer0[2 * k + 1] = x->data[j];
	    k++;
	    continue;
	}
	assert(_t->idiag[j] != 0.0);
	x->data[j] = solver->rhs->data[j] * _t->idiag[j];
    }
    v = phgMatVec(MAT_OP_N, 1.0, _t->mat, x, 0.0, NULL);
    for (i = k = 0; i < solver->mat->rmap->nlocal; i++) {
	j = _t->order[i];
	if (_t->index[j] < 0)
	    continue;
	buffer0[2 * k] = solver->rhs->data[j] - v->data[j];
	k++;
    }
    phgVecDestroy(&v);

#if USE_MPI
    if (solver->rhs->map->nprocs > 1) {
	MPI_Datatype type;
	buffer1 = phgAlloc(2 * n * sizeof(*buffer1));
	MPI_Type_contiguous(2, PHG_MPI_FLOAT, &type);
	MPI_Type_commit(&type);
	MPI_Alltoallv(buffer0, _t->scnts, _t->sdsps, type,
		      buffer1, _t->rcnts, _t->rdsps, type,
		      solver->mat->rmap->comm);
	MPI_Type_free(&type);
    }
    else
#endif	/* use_MPI */
    {
	buffer1 = buffer0;
    }

    if (solver->monitor) {
	double time1 = phgGetTime(NULL);
	phgPrintf("Setting up initial solution and RHS: wtime = %0.4lg\n",
			time1 - time0);
	time0 = time1;
    }

    if (n > 0) {
	VEC *x0 = phgMapCreateVec(_t->solver->rhs->map, 0);

	x0->nvec = 1;
	x0->data = buffer1;
	for (i = 0; i < n; i++) {
	    _t->solver->rhs->data[i] = buffer1[2 * i];
	    buffer1[i] = buffer1[2 * i + 1];
	}
	phgSolverAssemble(_t->solver);
	phgSolverVecSolve(_t->solver, destroy, x0);
	x0->data = NULL;
	phgVecDestroy(&x0);
    }

    if (solver->monitor) {
	double time1 = phgGetTime(NULL);
	phgPrintf("Solving derived system: wtime = %0.4lg\n", time1 - time0);
	time0 = time1;
    }

#if USE_MPI
    if (solver->rhs->map->nprocs > 1) {
	MPI_Alltoallv(buffer1, _t->rcnts, _t->rdsps, PHG_MPI_FLOAT,
		      buffer0, _t->scnts, _t->sdsps, PHG_MPI_FLOAT,
		      solver->mat->rmap->comm);
	phgFree(buffer1);
    }
#endif	/* use_MPI */

    for (i = k = 0; i < solver->mat->rmap->nlocal; i++) {
	j = _t->order[i];
	if (_t->index[j] < 0)
	    continue;
	x->data[j] = buffer0[k++];
    }

    phgFree(buffer0);

    solver->nits = 1;
    if (n > 0)
	solver->residual = _t->solver->residual;
    else
	solver->residual = 0.0;
#if USE_MPI
    {
	FLOAT b, a = Pow(solver->residual, 2);
	MPI_Allreduce(&a,&b,1, PHG_MPI_FLOAT, PHG_SUM, solver->mat->rmap->comm);
	solver->residual = Sqrt(b);
    }
#endif	/* use_MPI */

    if (destroy) {
	phgFree(_t->solver_opts);
	_t->solver_opts = NULL;
	phgFree(_t->order);
	_t->order = NULL;
	phgFree(_t->index);
	_t->index = NULL;
	phgFree(_t->idiag);
	_t->idiag = NULL;
	phgFree(_t->scnts);
	_t->scnts = NULL;
	phgFree(_t->rcnts);
	_t->rcnts = NULL;
	phgFree(_t->sdsps);
	_t->sdsps = NULL;
	phgFree(_t->rdsps);
	_t->rdsps = NULL;
	phgMatDestroy(&_t->mat);
	phgSolverDestroy(&_t->solver);
    }

    if (solver->monitor) {
	double time1 = phgGetTime(NULL);
	phgPrintf("Scattering solution: wtime = %0.4lg\n", time1 - time0);
	time0 = time1;
    }

    return solver->nits;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverGather_ = {
    "Gather", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy,
    AddMatrixEntries, AddRHSEntries,
    Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, FALSE, FALSE
};
