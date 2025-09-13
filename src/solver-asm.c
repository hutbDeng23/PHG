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

/* $Id: solver-asm.c,v 1.56 2022/09/21 02:35:28 zlb Exp $
 *
 * This is an implementation of the (overlaping) additive Schwarz method.
 *
 * Currently only one iteration with the zero initial solution is implemented,
 * so the solver can only be used as a preconditioner. */

#include "phg.h"

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>	/* INT_MAX */

typedef struct {
    SOLVER	*solver;	/* The subdomain solver */
    BOOLEAN	factored;
    INT		*O2Gmap;
    INT		*ordering;
#if USE_MPI
    COMM_INFO	*cinfo;
#endif	/* USE_MPI */

    /* options variables */
    char	*solver_opts;
    int		solver_id;
    int		restriction;	/* restriction (global->local) type */
    int		prolongation;	/* prolongation (local->global) type */
    INT		overlap;
} OEM_DATA;

#define	_t	((OEM_DATA *)solver->oem_data)

static char *solver_opts = NULL;
static int solver_id = 0;

/* Define types of restriction/prolongation.
 *
 * Some commonly used combinations:
 *
 *	RAS (restricted additive schwarz):
 *		restriction = GLOBAL, prolongation = LOCAL (default)
 *
 *	ASH (additive schwarz with harmonic extension):
 *		restriction = LOCAL, prolongation = GLOBAL
 *
 *	AS (standard additive schwarz):
 *		restriction = GLOBAL, prolongation = GLOBAL
 */
static const char *type_names[] = {
    "local", "global", /*"weighted",*/ NULL};
enum {LOCAL = 0, GLOBAL = 1, WEIGHTED = 2};
static int restriction = GLOBAL;
static int prolongation = LOCAL;

static INT overlap = 1;

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

    phgOptionsRegisterTitle("\nThe additive Schwarz preconditioner options:",
				"\n", "asm");
    phgOptionsRegisterKeyword("-asm_sub_solver",
					"The subdomain solver (obsolete)",
					phgSolverNames, &solver_id);
    phgOptionsRegisterString("-asm_sub_solver_opts",
					"Options for the subdomain solver",
					&solver_opts);
    phgOptionsRegisterInt("-asm_overlap", "Number of overlapping layers",
					&overlap);
    phgOptionsRegisterKeyword("-asm_restriction", "Restriction type",
					type_names, &restriction);
    phgOptionsRegisterKeyword("-asm_prolongation", "Prolongation type",
					type_names, &prolongation);

    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));
    _t->solver_id = solver_id;
    _t->restriction = restriction;
    _t->prolongation = prolongation;
    _t->overlap = overlap;
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

    phgFree(_t->O2Gmap);
    phgFree(_t->ordering);
    phgFree(_t->solver_opts);
    if (_t->solver != NULL)
	phgSolverDestroy(&_t->solver);
#if USE_MPI
    phgMapDestroyCommInfo(&_t->cinfo);
#endif	/* USE_MPI */

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    MAP *map;
    MAT *B;
    MAT_ROW *row, *rows;
    INT i, j, n, nlocal, localsize;
    INT map_alloc, map_size;
    int nprocs = solver->mat->rmap->nprocs;
#if USE_MPI
    int o;
    INT k, r, n0;
#endif	/* USE_MPI */
    double time0, time1;

    if (_t->factored)
	return 0;

    _t->factored = TRUE;

    if (solver->sent_to_oem == TRUE)
	phgError(1, "%s:%d, data already sent to OEM solver!\n",
						__FILE__, __LINE__);
    solver->sent_to_oem = TRUE;

    /*------------------------- Set up local matrix ------------------------*/
    time0 = phgGetTime(NULL);

    assert(phgMapCompare(solver->mat->cmap, solver->mat->rmap));

    n = nlocal = solver->mat->cmap->nlocal;
    localsize = solver->mat->localsize;
    if (localsize == 0)
	localsize = nlocal;

#if USE_MPI
    n0 = solver->mat->cmap->partition[solver->mat->cmap->rank];
    if (localsize > nlocal) {
	size_t size = (localsize - n) * sizeof(*_t->O2Gmap);
	_t->O2Gmap = phgAlloc(size);
	memcpy(_t->O2Gmap, solver->mat->O2Gmap, size);
	_t->ordering = phgAlloc(size);
	memcpy(_t->ordering, solver->mat->ordering, size);
    }
#endif	/* USE_MPI */

    map_size = map_alloc = localsize - nlocal;
    rows = phgAlloc(n * sizeof(*rows));
    for (i = 0; i < nlocal; i++) {
	if (solver->mat->type == PHG_UNPACKED) {
	    row = solver->mat->rows + i;
	    j = row->ncols;
	    rows[i].ncols = rows[i].alloc = j;
	    rows[i].cols = phgAlloc(j * sizeof(*row->cols));
	    rows[i].data = phgAlloc(j * sizeof(*row->data));
	    if (j > 0) {
		memcpy(rows[i].cols, row->cols, j * sizeof(*row->cols));
		memcpy(rows[i].data, row->data, j * sizeof(*row->data));
	    }
	    continue;
	}
	else if (solver->mat->type == PHG_PACKED) {
	    INT *pc, o;
	    FLOAT *pd;
	    pc = solver->mat->packed_cols + solver->mat->packed_ind[i];
	    j = (INT)(solver->mat->packed_ind[i + 1]
				- solver->mat->packed_ind[i]);
#if USE_MPI
	    o = (INT)(solver->mat->packed_ind[i + nlocal + 1]
				- solver->mat->packed_ind[i + nlocal]);
#else	/* USE_MPI */
	    o = 0;
#endif	/* USE_MPI */
	    rows[i].ncols = rows[i].alloc = j + o;
	    rows[i].cols = phgAlloc((j + o) * sizeof(*pc));
	    rows[i].data = phgAlloc((j + o) * sizeof(*pd));
	    if (j > 0) {
		pd = solver->mat->packed_data + solver->mat->packed_ind[i];
		memcpy(rows[i].cols, pc, j * sizeof(*pc));
		memcpy(rows[i].data, pd, j * sizeof(*pd));
	    }
#if USE_MPI
	    if (o > 0) {
		pc = solver->mat->packed_cols +
		     solver->mat->packed_ind[i + nlocal];
		pd = solver->mat->packed_data +
		     solver->mat->packed_ind[i + nlocal];
		memcpy(rows[i].data + j, pd, o * sizeof(*pd));
		while (j < rows[i].ncols)
		    rows[i].cols[j++] = solver->mat->ordering[*(pc++)] + nlocal;
	    }
#endif	/* USE_MPI */
	    continue;
	}
	else {
	    row = (MAT_ROW *)phgMatGetRow(solver->mat, i);
	    j = row->ncols;
	    rows[i].ncols = rows[i].alloc = j;
	    rows[i].cols = phgAlloc(j * sizeof(*row->cols));
	    rows[i].data = phgAlloc(j * sizeof(*row->data));
	    if (j > 0) {
		memcpy(rows[i].cols, row->cols, j * sizeof(*row->cols));
		memcpy(rows[i].data, row->data, j * sizeof(*row->data));
	    }
	    if (nprocs <= 1)
		continue;
	}
#if USE_MPI
	/* convert global indices back to local ones */
	row = rows + i;
	for (j = 0; j < row->ncols; j++) {
	    r = row->cols[j];
	    if (r >= n0 && r < n0 + nlocal) {
		row->cols[j] -= n0;
		continue;
	    }
	    if (_t->overlap == 0) {
		row->cols[j] = -1;
		continue;
	    }
	    /* search O2Gmap */
	    k = phgBinarySearchINT(map_size, _t->O2Gmap, _t->ordering, r);
	    if (k < map_size && r == _t->O2Gmap[_t->ordering[k]]) {
		/* found in _t->O2Gmap */
		row->cols[j] = nlocal + _t->ordering[k];
		continue;
	    }
	    /* add a new global index to _t->O2Gmap */
	    if (map_alloc <= map_size) {
		_t->O2Gmap = phgRealloc_(_t->O2Gmap,
				    (map_alloc + 128) * sizeof(*_t->O2Gmap),
				    map_alloc * sizeof(*_t->O2Gmap));
		_t->ordering = phgRealloc_(_t->ordering,
				    (map_alloc + 128) * sizeof(*_t->ordering),
				    map_alloc * sizeof(*_t->ordering));
		map_alloc += 128;
	    }
	    row->cols[j] = nlocal + map_size;
	    _t->O2Gmap[map_size] = r;
	    /* insert map_size into ordering[] */
	    for (r = map_size; r > k; r--)
		_t->ordering[r] = _t->ordering[r - 1];
	    _t->ordering[k] = map_size++;
	}
#endif	/* USE_MPI */
    }
    localsize = map_size + nlocal;

#if USE_MPI
    if (_t->overlap > 0 && nprocs > 1) {
	/* get the overlapping rows */
	INT *partition = solver->mat->rmap->partition;
	MPI_Comm comm = solver->mat->rmap->comm;
	int *scnts, *sdsps, *rcnts, *rdsps, *stemp;
	INT *list, *temp, nlist, ntemp;

	/* Note: see phgMatAssemble() for format of the buffers */
	unsigned char *sbuff, *rbuff, *align;
	INT ncols, nrows = 0;
	FLOAT *pdata = NULL;
	INT r0, *prows = NULL, *pncols = NULL, *pcols = NULL;
	int rank, rank0;
	MPI_Datatype type;
	size_t size, ssize, rsize, *nnz;
	static size_t align_float = 0, align_int = 0, trunc_size = 128;

	scnts = phgAlloc(5 * nprocs * sizeof(*scnts));
	sdsps = scnts + 1 * nprocs;
	rcnts = scnts + 2 * nprocs;
	rdsps = scnts + 3 * nprocs;
	stemp = scnts + 4 * nprocs;

	for (o = 0; o < _t->overlap; o++) {
	    INT *ordering;

	    if (align_float == 0) {
		/* compute alignment of FLOAT and INT */
		struct {BYTE a; FLOAT b;} a;
		struct {BYTE a; INT b;} b;
		align_float = sizeof(a) - sizeof(FLOAT);
		align_int = sizeof(b) - sizeof(INT);
		if (trunc_size % align_int != 0 ||
		    trunc_size % align_float != 0 ||
	    	    (align_float % align_int != 0 &&
		     align_int % align_float != 0))
		    phgError(1, "%s:%d, unexpected alignment.\n",
						__FILE__, __LINE__);
	    }

	    /* make a backup copy of current ordering */
	    if (localsize > nlocal) {
		ordering = phgAlloc((localsize - nlocal) * sizeof(*ordering));
		memcpy(ordering, _t->ordering,
				(localsize - nlocal) * sizeof(*ordering));
	    }
	    else {
		ordering = NULL;
	    }

	    /* count # of entries for each process.
	     * Note: in this step the meaning of rcnts/rdsps and scnts/sdsps
	     *	     is reversed. */
	    memset(rcnts, 0, nprocs * sizeof(*rcnts));
	    rank = 0;
	    for (i = 0; i < localsize - nlocal; i++) {
		if (nlocal + ordering[i] < n)
		    continue;
		j = _t->O2Gmap[ordering[i]];
		while (rank < nprocs && j >= partition[rank + 1])
		    rank++;
		assert(rank < nprocs && rank != solver->mat->rmap->rank);
		rcnts[rank]++;
	    }
	    /* Attach the flag (localsize==n) to the highest bits of rcnts[] */
#define	TEMP_MASK	(1 << (8 * sizeof(rcnts[0]) - 1))
	    if (localsize > n) {
		for (rank = 0; rank < nprocs; rank++)
		    rcnts[rank] |= TEMP_MASK;
	    }
	    MPI_Alltoall(rcnts, 1, MPI_INT, scnts, 1, MPI_INT, comm);
	    k = nlist = ntemp = 0;
	    for (rank = 0; rank < nprocs; rank++) {
		/* check and remove the highest bit in scnts[] and rcnts[] */
		if ((scnts[rank] & TEMP_MASK))
		    k = 1;
		rcnts[rank] &= ~TEMP_MASK;
		scnts[rank] &= ~TEMP_MASK;
		sdsps[rank] = nlist;
		rdsps[rank] = ntemp;
		nlist += scnts[rank];
		ntemp += rcnts[rank];
	    }
	    if (k == 0) {
		/* no more overlapping rows to exchange */
		phgFree(ordering);
		ordering = NULL;
		break;
	    }
	    list = phgAlloc(nlist * sizeof(*list));
	    temp = phgAlloc(ntemp * sizeof(*temp));
	    /* collect row indices */
	    rank = 0;
	    for (i = 0; i < localsize - nlocal; i++) {
		if (nlocal + ordering[i] < n)
		    continue;
		j = _t->O2Gmap[ordering[i]];
		while (rank < nprocs && j >= partition[rank + 1])
		    rank++;
		temp[rdsps[rank]++] = j;
	    }
	    /* restore sdsps */
	    for (rank = 0; rank < nprocs; rank++)
		rdsps[rank] -= rcnts[rank];
	    MPI_Alltoallv(temp, rcnts, rdsps, PHG_MPI_INT,
			  list, scnts, sdsps, PHG_MPI_INT, comm);
	    phgFree(temp);

	    for (i = 0; i < nlist; i++) {
		list[i] -= n0;
		assert(list[i] >= 0 && list[i] < nlocal);
	    }

	    /* Now list[] are the list of rows to send,
	     * scnts[] contains # of rows to send for each process */

	    /* count nnz (nnz[]) for each process */
	    nnz = phgCalloc(nprocs, sizeof(*nnz));
	    for (rank = i = 0; i < nlist; i++) {
		while (rank < nprocs - 1 && sdsps[rank + 1] <= i)
		    rank++;
		if (solver->mat->type == PHG_UNPACKED)
		    ncols = solver->mat->rows[list[i]].ncols;
		else
		    ncols = phgMatGetNnzInRow(solver->mat, list[i], NULL, NULL);
		if (ncols > 0)
		    nnz[rank] += (size_t)ncols;
	    }

	    /* compute in stemp offsets of data for each process */
	    ssize = 0;
	    for (rank = 0; rank < nprocs; rank++) {
		nrows = scnts[rank];

		/* size of header + INT entries */
		size = (4 + nrows * 2 + nnz[rank]) * sizeof(INT);
		/* adjust size to alignment of FLOAT */
		size = ((size + align_float - 1) / align_float) * align_float;
		/* plus size of FLOAT entries */
		size += nnz[rank] * sizeof(FLOAT);
		/* adjust size to multiple of trunc_size */
		size = ((size + trunc_size - 1) / trunc_size) * trunc_size;

		stemp[rank] = (int)(ssize / trunc_size);
		ssize += size;
		assert(ssize / trunc_size <= INT_MAX);
	    }

	    /* allocate send buffer */
	    sbuff = phgAlloc(ssize);

	    /* setup data counts in sbuff */
	    for (rank = 0; rank < nprocs; rank++) {
		size = stemp[rank];
		if (size == (rank == nprocs - 1 ?
				(int)(ssize / trunc_size) : stemp[rank + 1]))
		    continue;
		size *= trunc_size;
		prows = (void *)(sbuff + size);
		prows[0] = scnts[rank];
		prows[1] = (INT)(nnz[rank] % (size_t)(INT_MAX));
		prows[2] = (INT)(nnz[rank] / (size_t)(INT_MAX));
		prows[3] = 0;
	    }
	    phgFree(nnz);

	    if (ssize != 0 && nlist > 0) {
		/* collect rows to sbuff */
		rank0 = -1;
		rank = 0;
		for (i = 0; i < nlist; i++) {
		    while (rank < nprocs - 1 && sdsps[rank + 1] <= i)
			rank++;
		    row = (MAT_ROW *)phgMatGetRow(solver->mat, list[i]);
		    if (rank != rank0) {
			prows = (void *)(sbuff +
					 ((size_t)stemp[rank]) * trunc_size);
			nrows = prows[0];
			size = ((size_t)prows[1]) +
			       ((size_t)prows[2]) * (size_t)(INT_MAX);

			prows += 4;
			pncols = prows + nrows;
			pcols = pncols + nrows;

			align = (void *)(pcols + size);
			align = sbuff + ((align - sbuff + align_float - 1)
					/ align_float) * align_float;
			pdata = (void *)align;
			rank0 = rank;
		    }
		    if ((ncols = row->ncols) > 0) {
			*(prows++) = n0 + list[i];
			*(pncols++) = ncols;
			memcpy(pcols, row->cols, sizeof(INT) * ncols);
			pcols += ncols;
			memcpy(pdata, row->data, sizeof(FLOAT) * ncols);
			pdata += ncols;
		    }
		}
	    }

	    phgFree(list);

	    /* now sbuff contains packed rows to send */
	    for (rank = 0; rank < nprocs; rank++) {
		k = (rank < nprocs - 1 ? stemp[rank + 1] :
					 (int)(ssize / trunc_size));
		scnts[rank] = k - stemp[rank];
		sdsps[rank] = stemp[rank];
	    }

	    /* exchange counts with other processes */
	    MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, comm);

	    /* prepare receive buffer */
	    k = 0;
	    for (rank = 0; rank < nprocs; rank++) {
		rdsps[rank] = k;
		k += rcnts[rank];
	    }
	    rsize = ((size_t)k) * trunc_size;
	    rbuff = phgAlloc(rsize);

	    /* exchange data */
	    MPI_Type_contiguous((int)trunc_size, MPI_BYTE, &type);
	    MPI_Type_commit(&type);
	    MPI_Alltoallv(sbuff, scnts, sdsps, type,
			  rbuff, rcnts, rdsps, type, comm);
	    MPI_Type_free(&type);

	    phgFree(sbuff);

	    /* process received data */

	    rows = phgRealloc_(rows, localsize * sizeof(*rows),
				      n * sizeof(*rows));
	    r0 = -1;
	    for (rank = 0; rank < nprocs; rank++) {
		if (rcnts[rank] == 0)
		    continue;
		prows = (void *)(rbuff + ((size_t)rdsps[rank]) * trunc_size);
		nrows = prows[0];
		size = ((size_t)prows[1]) +
			((size_t)prows[2]) * (size_t)(INT_MAX);

		prows += 4;
		pncols = prows + nrows;
		pcols = pncols + nrows;

		align = (void *)(pcols + size);
		align = rbuff + ((align - rbuff + align_float - 1)
				 / align_float) * align_float;
		pdata = (void *)align;

		for (i = 0; i < nrows; i++) {
		    ncols = *pncols;
		    /* convert global indices in pcols[] to local indices */
		    for (j = 0; j < ncols; j++) {
			r = pcols[j];
			if (r >= n0 && r < n0 + nlocal) {
			    pcols[j] -= n0;
			    continue;
			}
			/* search the _t->O2Gmap */
			k = phgBinarySearchINT(map_size, _t->O2Gmap,
						_t->ordering, r);
			if (k < map_size && r == _t->O2Gmap[_t->ordering[k]]) {
			    /* found in _t->O2Gmap */
			    pcols[j] = nlocal + _t->ordering[k];
			    continue;
			}
			else if (o == _t->overlap - 1) {
			    /* discard the entry */
			    pcols[j] = -1;
			    continue;
			}
			/* add a new global index to _t->O2Gmap */
			if (map_alloc <= map_size) {
			    _t->O2Gmap = phgRealloc_(_t->O2Gmap,
				(map_alloc + 128) * sizeof(*_t->O2Gmap),
				map_alloc * sizeof(*_t->O2Gmap));
			    _t->ordering = phgRealloc_(_t->ordering,
				(map_alloc + 128) * sizeof(*_t->ordering),
				map_alloc * sizeof(*_t->ordering));
			    map_alloc += 128;
			}
			pcols[j] = nlocal + map_size;
			_t->O2Gmap[map_size] = r;
			/* insert map_size into ordering[] */
			for (r = map_size; r > k; r--)
			    _t->ordering[r] = _t->ordering[r - 1];
			_t->ordering[k] = map_size++;
		    }
		    while (++r0 < localsize - nlocal
				&& nlocal + ordering[r0] < n);
		    assert(r0 < localsize - nlocal);
		    row = rows + nlocal + ordering[r0];
#if DEBUG
		    /* compare with the received row index */
		    r = *(prows++);
		    k = phgBinarySearchINT(localsize - nlocal, _t->O2Gmap,
					    ordering, r);
		    assert(k == r0);
#endif	/* DEBUG */
		    row->ncols = row->alloc = ncols;
		    row->cols = phgAlloc(ncols * sizeof(*row->cols));
		    row->data = phgAlloc(ncols * sizeof(*row->data));
		    assert(ncols > 0);
		    memcpy(row->cols, pcols, ncols * sizeof(*row->cols));
		    memcpy(row->data, pdata, ncols * sizeof(*row->data));
		    pncols++;
		    pcols += ncols;
		    pdata += ncols;
		}
	    }

	    phgFree(rbuff);

	    n = localsize;
	    localsize = map_size + nlocal;

	    phgFree(ordering);
	}

	phgFree(scnts);
    }
#endif	/* USE_MPI */

    if (solver->mat->refcount == 0)
	phgMatFreeMatrix(solver->mat);

    /* discard entries in the range [n, localsize) */
    for (i = 0; i < n; i++) {
	INT k;
	row = rows + i;
	for (j = k = 0; j < row->ncols; j++) {
	    if (row->cols[j] >= n || row->cols[j] < 0)
		continue;
	    if (k != j) {
		row->cols[k] = row->cols[j];
		row->data[k] = row->data[j];
	    }
	    k++;
	}
	row->cols = phgRealloc_(row->cols, k * sizeof(*row->cols),
					   k * sizeof(*row->cols));
	row->data = phgRealloc_(row->data, k * sizeof(*row->data),
					   k * sizeof(*row->data));
	row->ncols = row->alloc = k;
    }

#if USE_MPI
    if ((_t->restriction != LOCAL || _t->prolongation != LOCAL) &&
	_t->overlap > 0 && nprocs > 1) {
	_t->cinfo = phgMapCreateCommInfo(
			solver->mat->rmap->comm,
			solver->mat->rmap->nprocs,
			solver->mat->rmap->rank,
			nlocal,
			n,
			_t->O2Gmap,
			_t->ordering,
			solver->mat->rmap->partition);

	/* compute weights (and apply to B) for weighted interpolation */
	if (_t->restriction == WEIGHTED || _t->prolongation == WEIGHTED) {
	    FLOAT *w = phgAlloc(n * sizeof(*w));
	    INT i, j;
	    for (i = 0; i < n; i++)
		w[i] = 1.0;
	    phgMapGather(_t->cinfo, 1, w, w + nlocal);
	    /*for (i = 0; i < n; i++)
		w[i] = 1.0 / w[i];*/
	    phgMapScatter(_t->cinfo, 1, w, w + nlocal);
	    if (_t->restriction == WEIGHTED) {
		/* B := diag(w) * B */
		for (i = 0; i < n; i++)
		    for (j = 0; j < rows[i].ncols; j++)
			rows[i].data[j] *= w[i];
	    }
	    if (_t->prolongation == WEIGHTED) {
		/* B := B * diag(w) */
		for (i = 0; i < n; i++)
		    for (j = 0; j < rows[i].ncols; j++)
			rows[i].data[j] *= w[rows[i].cols[j]];
	    }
	    phgFree(w);
	}
    }
#endif	/* USE_MPI */

    if (n > 0) {
	map = phgMapCreateSimpleMap(MPI_COMM_SELF, n, n);
	B = phgMapCreateMat(map, map);
	phgMapDestroy(&map);
	phgFree(B->rows);
	B->rows = rows;

	/* Create the subdomain solver.
	 *
	 * Note: it's better to create the solver in Init() to honor all solver
	 * options for the serial solver, but we couldn't do it there since we
	 * need actual matrix data to determine the map. */

	phgOptionsPush();
	phgSolverSetDefaultSuboptions();
	/* set options for the subdomain solver */
	phgOptionsSetKeyword("-solver", phgSolverNames[_t->solver_id]);
	phgOptionsSetOptions(_t->solver_opts);
	assert(_t->solver == NULL);
	_t->solver = phgMat2Solver(SOLVER_DEFAULT, B);
	if (_t->solver->oem_solver == solver->oem_solver)
	    phgError(1, "(%s:%d) solver regression detected, abort.\n",
						__FILE__, __LINE__);
	phgOptionsPop();
	phgMatDestroy(&B);
	phgSolverAssemble(_t->solver);
    }

    if (solver->monitor) {
	INT min0, max0, avg0, min1, max1, avg1;
	double mint, maxt, avgt;

	time1 = phgGetTime(NULL);
#if USE_MPI
	if (nprocs > 1) {
	    double a[3][3], b[3][3];

	    a[0][0] = a[0][1] = a[0][2] = nlocal;
	    a[1][0] = a[1][1] = a[1][2] = n;
	    a[2][0] = a[2][1] = a[2][2] = time1 - time0;

	    MPI_Reduce(a, b, 3, MPI_3DOUBLE, MPI_MSM, 0,
				solver->mat->rmap->comm);

	    if (solver->mat->rmap->rank == 0) {
		min0 = RoundINT(b[0][0]);
		avg0 = RoundINT(b[0][1] / nprocs);
		max0 = RoundINT(b[0][2]);

		min1 = RoundINT(b[1][0]);
		avg1 = RoundINT(b[1][1] / nprocs);
		max1 = RoundINT(b[1][2]);

		mint = b[2][0];
		avgt = b[2][1] / nprocs;
		maxt = b[2][2];
	    }
	    else {
		min0 = max0 = avg0 = min1 = max1 = avg1 = 0;
		mint = maxt = avgt = 0.0;
	    }
	}
	else 
#endif	/* USE_MPI */
	{
	    min0 = max0 = avg0 = min1 = max1 = avg1 = n;
	    mint = maxt = avgt = time1 - time0;
	}
	phgPrintf("* ASM: number of subdomains = %d\n",
		  solver->mat->rmap->nprocs);
	phgPrintf("* ASM: restriction = %s, prolongation = %s, "
		  "subdomain solver = %s\n", type_names[_t->restriction],
		  type_names[_t->prolongation], 
		  _t->solver == NULL ? "n.a." : _t->solver->oem_solver->name);
	phgPrintf("* ASM: subdomain assembly time (min/avg/max): "
		  "%0.4lf/%0.4lf/%0.4lf\n", mint, avgt, maxt);
	if (_t->overlap > 0)
	    phgPrintf("* ASM: subdomain size without overlap "
		      "(min/avg/max): %d/%d/%d\n", min0, avg0, max0);
	phgPrintf("* subdomain size with overlap %d (min/avg/max): "
		  "%d/%d/%d\n", _t->overlap, min1, avg1, max1);
	time0 = time1;
    }

    return 0;
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    VEC *x0 = NULL;
    double time0 = 0.;

    assert(x->nvec == 1);
    Assemble(solver);

    if (phgVerbosity > 1)
	time0 = phgGetTime(NULL);

    if (_t->solver != NULL) {
	x0 = phgMapCreateVec(_t->solver->rhs->map, 1);
	if (x->map->nlocal > 0)
	    memcpy(_t->solver->rhs->data, solver->rhs->data,
		   x->map->nlocal * sizeof(*x->data));
    }

#if USE_MPI
    /* Fetch overlapped part of the RHS */
    if (_t->overlap > 0 && solver->mat->cmap->nprocs > 1) {
	size_t size;
	switch (_t->restriction) {
	    case LOCAL:
		if (_t->solver == NULL)
		    break;
		size = (x0->map->nlocal - x->map->nlocal) * sizeof(*x->data);
		if (size > 0)
		    memset(_t->solver->rhs->data + x->map->nlocal, 0, size);
		break;
	    case GLOBAL:
	    case WEIGHTED:
		if (_t->solver == NULL)
		    phgMapScatter(_t->cinfo, 1, NULL, NULL);
		else
		    phgMapScatter(_t->cinfo, x0->nvec, _t->solver->rhs->data,
				  _t->solver->rhs->data + x->map->nlocal);
		break;
	}
    }
#endif	/* USE_MPI */

    if (_t->solver != NULL) {
	_t->solver->rhs->assembled = TRUE;
	_t->solver->rhs_updated = TRUE;
	phgSolverVecSolve(_t->solver, destroy, x0);
    }

    if (x->map->nlocal > 0)
	memcpy(x->data, x0->data, x->map->nlocal * sizeof(*x->data));

#if USE_MPI
    if (_t->overlap > 0 && solver->mat->cmap->nprocs > 1) {
	switch (_t->prolongation) {
	    case LOCAL:
		/* nothing to do */
		break;
	    case GLOBAL:
	    case WEIGHTED:
		if (_t->solver == NULL)
		    phgMapGather(_t->cinfo, 1, x->data, NULL);
		else
		    phgMapGather(_t->cinfo, x0->nvec, x->data,
				 x0->data + x->map->nlocal);
		break;
	}
    }
#endif	/* USE_MPI */

    if (x0 != NULL)
	phgVecDestroy(&x0);

    if (phgVerbosity > 1) {
	double time1 = phgGetTime(NULL);
	phgInfo(0, "ASM - subdomain solver: wtime = %0.4lg\n", time1 - time0);
	time0 = time1;
    }

    if (destroy) {
	phgFree(_t->O2Gmap);
	_t->O2Gmap = NULL;
	phgFree(_t->ordering);
	_t->ordering = NULL;
	phgFree(_t->solver_opts);
	_t->solver_opts = NULL;
	if (_t->solver != NULL)
	    phgSolverDestroy(&_t->solver);
#if USE_MPI
	phgMapDestroyCommInfo(&_t->cinfo);
#endif	/* USE_MPI */
    }

    return solver->nits;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverASM_ = {
    "ASM", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy,
    AddMatrixEntries, AddRHSEntries,
    Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, FALSE, FALSE
};
