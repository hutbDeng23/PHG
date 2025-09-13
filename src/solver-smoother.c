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

/* $Id: solver-smoother.c,v 1.23 2022/09/21 02:35:28 zlb Exp $
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
#define SetPC   		NULL

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>		/* INT_MAX */



#define DGETRF_FOR dgetrf_
#define DGETRI_FOR dgetri_
#define DGEMM_FOR  dgemm_
void DGETRF_FOR(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
void DGETRI_FOR(int *N, double *A, int *LDA, int *IPIV, double *WORK,
		int *LWORK, int *INFO);
void DGEMM_FOR(char *TRANSA, char *TRANSB, int *M, int *N, int *K,
	       double *ALPHA, double *A, int *LDA, double *B, int *LDB,
	       double *BETA, double *C, int *LDC);

#define DENSE_INTERNAL 1
#define DENSE_LAPACK 0
#define DENSE_EXTERNAL 0
#define MAT_MAX_SIZE 2000

typedef struct {
    SOLVER *solver;		/* The subdomain solver */
    MAT_ROW *rows;
    BOOLEAN factored;

    INT nlocal;
    INT nlocal_o;
    INT localsize;

    INT *O2Gmap;
    INT *ordering;
#if USE_MPI
    COMM_INFO *cinfo;
#endif				/* USE_MPI */

    /* options variables */
    char *solver_opts;
    int solver_id;
    int restriction;		/* restriction (global->local) type */
    int prolongation;		/* prolongation (local->global) type */
    INT overlap;

    int nblock;
    int nblock0;		/* local block */
    int nblock1;		/* inter proc block */
    int bdry_offset;
    int *block_ofs;		/* block indices offset */

    INT *block_indices;
    int maxN;

    INT *packed_col;
    FLOAT *packed_dat;

    FLOAT *x0_dat;
    FLOAT *rhs_dat;

    SOLVER **bk_solvers;
    VEC **bk_vec;
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
    "local", "global", /*"weighted", */ NULL
};

static const char *smoother_names[] = {
    "GS", "boxelem", "boxvert", "pc", NULL
};
enum { LOCAL = 0, GLOBAL = 1, WEIGHTED = 2 };
enum { GS = 0, BOXELEM = 1, BOXVERT = 2 };
static int restriction = GLOBAL;
static int prolongation = LOCAL;
static int smoother_type = BOXELEM;

static INT overlap = 1;
static FLOAT omega = 1.;

static int
RegisterOptions(void)
{
    /* set default solver to the first direct solver (if available) */
    for (solver_id = 0;
	 phgSolverNames[solver_id] != NULL
	 && (phgSolverList[solver_id] == SOLVER_GATHER ||
	     phgSolverList[solver_id] == SOLVER_DIAGONAL ||
	     phgSolverList[solver_id]->iterative == TRUE); solver_id++);
    if (phgSolverNames[solver_id] == NULL)
	solver_id = 0;		/* fall back to the first solver */

    phgOptionsRegisterTitle("\nThe Smoother solver options:",
			    "\n", "smoother");
    phgOptionsRegisterKeyword("-smoother_sub_solver",
			      "The subblock solver (obsolete)",
			      phgSolverNames, &solver_id);
    phgOptionsRegisterString("-smoother_sub_solver_opts",
			     "Options for the subblock solver", &solver_opts);
    phgOptionsRegisterInt("-smoother_overlap", "Number of overlapping layers",
			  &overlap);
    phgOptionsRegisterKeyword("-smoother_restriction", "Restriction type",
			      type_names, &restriction);
    phgOptionsRegisterKeyword("-smoother_prolongation", "Prolongation type",
			      type_names, &prolongation);
    phgOptionsRegisterKeyword("-smoother_type", "Smoother type",
			      smoother_names, &smoother_type);
    phgOptionsRegisterFloat("-smoother_damp_factor", "Damp factor", &omega);

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
    INT i;

    if (solver->oem_data == NULL)
	return 0;

    for (i = 0; i < _t->localsize; i++) {
	phgFree(_t->rows[i].cols);
	phgFree(_t->rows[i].data);
    }
    phgFree(_t->rows);

    for (i = 0; i < _t->nblock; i++) {
	phgSolverDestroy(_t->bk_solvers + i);
	phgVecDestroy(_t->bk_vec + i);
    }
    phgFree(_t->block_ofs);
    phgFree(_t->block_indices);
    phgFree(_t->bk_solvers);
    phgFree(_t->bk_vec);
    phgFree(_t->x0_dat);
    phgFree(_t->rhs_dat);
    phgFree(_t->packed_col);
    phgFree(_t->packed_dat);

    phgFree(_t->O2Gmap);
    phgFree(_t->ordering);
    phgFree(_t->solver_opts);
    phgSolverDestroy(&_t->solver);
#if USE_MPI
    phgMapDestroyCommInfo(&_t->cinfo);
#endif /* USE_MPI */

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    /*MAP *map; */
    MAT_ROW *row, *rows;
    INT i, j, n, l, nlocal, localsize;
    INT map_alloc, map_size;
    int nprocs = solver->mat->rmap->nprocs;
#if USE_MPI
    int o;
    INT r;
#endif /* USE_MPI */
    INT k, n0;
    double time0;
    /*BOOLEAN box_extend; */
    static double zero = 0.;
    INT nblock = 0, *block_indices = NULL;
    int nlocal_o;

    INT ib;
    INT npi, npd;
    INT npi_o, npd_o;
    INT *bidx;
    int maxN, minN;

    INT *packed_col, *pc, *pc_o, *pc_b;
    FLOAT *packed_dat, *pd, *pd_o;

    FLOAT *A;

    INT *col, *col_o = NULL, allocated;
    FLOAT *dat, *dat_o = NULL;

    if (_t->factored)
	return 0;

    _t->factored = TRUE;

    time0 = phgGetTime(NULL);

    if (smoother_type == BOXELEM) {
	/* build block indices for box element. */
	DOF **dofs = solver->mat->rmap->dofs;
	int ndof = solver->mat->rmap->ndof;
	GRID *g;
	SIMPLEX *e;
	int i, k;
	int ne_bdry = 0, ne_innr = 0, nidx_bdry = 0, nidx_innr = 0;

	INT *bidx_innr, *bidx_bdry;
	int ib;

	assert(dofs != NULL);
	g = dofs[0]->g;

	/* count */
	ForAllElements(g, e) {
	    int N[ndof], NN;
	    BTYPE btype;

	    NN = 0;
	    for (k = 0; k < ndof; k++) {
		N[k] = DofNBas(dofs[k], e);
		NN += N[k];
	    }

	    /* The block have ALL element local is marked local block,
	     * otherwize, it's interproc block.  */
	    for (k = 0; k < ndof; k++) {
		for (i = 0; i < N[k]; i++) {
		    btype = phgDofGetElementBasisInfo(dofs[k], e, i,
						      NULL, NULL, NULL);
		    if (btype & REMOTE)
			break;
		}
		if (i < N[k])	/* any remote Dof */
		    break;
	    }
	    if (k < ndof) {	/* any remote Dof */
		ne_bdry++;
		nidx_bdry += NN + 1;
	    }
	    else {
		ne_innr++;
		nidx_innr += NN + 1;
	    }
	}


	bidx = phgCalloc(nidx_bdry + nidx_innr, sizeof(*bidx));
	assert(g->nleaf == ne_innr + ne_bdry);
	bidx_innr = bidx;
	bidx_bdry = bidx_innr + nidx_innr;
	_t->block_ofs = phgCalloc(ne_innr + ne_bdry, sizeof(*_t->block_ofs));

	/* record */
	ForAllElements(g, e) {
	    int N[ndof], NN;
	    INT *I0, *I;
	    BTYPE btype;

	    NN = 0;
	    for (k = 0; k < ndof; k++) {
		N[k] = DofNBas(dofs[k], e);
		NN += N[k];
	    }
	    I0 = phgAlloc(NN * sizeof(*I));

	    /* idx */
	    I = I0;
	    for (k = 0; k < ndof; k++)
		for (i = 0; i < N[k]; i++) {
		    *(I++) = phgMapE2G(solver->rhs->map, k, e, i);
		}

	    for (k = 0; k < ndof; k++) {
		for (i = 0; i < N[k]; i++) {
		    btype = phgDofGetElementBasisInfo(dofs[k], e, i,
						      NULL, NULL, NULL);
		    if (btype & REMOTE)
			break;
		}
		if (i < N[k])	/* any remote Dof */
		    break;
	    }
	    if (k < ndof) {	/* any remote Dof */
		*(bidx_bdry++) = NN;
		memcpy(bidx_bdry, I0, NN * sizeof(*I));
		bidx_bdry += NN;
	    }
	    else {
		*(bidx_innr++) = NN;
		memcpy(bidx_innr, I0, NN * sizeof(*I));
		bidx_innr += NN;
	    }

	    phgFree(I0);
	}

	nblock = ne_innr + ne_bdry;
	block_indices = bidx;

	_t->nblock0 = ne_innr;
	_t->nblock1 = ne_bdry;
	_t->nblock = nblock;
	_t->block_indices = bidx;

	/* record block_indices offset */
	bidx = block_indices;
	for (ib = 0; ib < nblock; ib++) {
	    INT N;
	    _t->block_ofs[ib] = (int)(bidx - block_indices);
	    N = *(bidx++);
	    bidx += N;
	}
    }
    else if (smoother_type == BOXVERT) {
	phgError(1, "smoother box vert has not been implemented yet.\n");
    }

    if (solver->sent_to_oem == TRUE)
	phgError(1, "%s:%d, data already sent to OEM solver!\n",
		 __FILE__, __LINE__);
    solver->sent_to_oem = TRUE;

    /*------------------------- Set up local matrix ------------------------*/

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
#else
    n0 = 0;
#endif /* USE_MPI */

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
#else /* USE_MPI */
	    o = 0;
#endif /* USE_MPI */
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
		    rows[i].cols[j++] =
			solver->mat->ordering[*(pc++)] + nlocal;
	    }
#endif /* USE_MPI */
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
	    while (map_alloc <= map_size) {
		_t->O2Gmap = phgRealloc_(_t->O2Gmap,
					 (map_alloc +
					  128) * sizeof(*_t->O2Gmap),
					 map_alloc * sizeof(*_t->O2Gmap));
		_t->ordering =
		    phgRealloc_(_t->ordering,
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
#endif /* USE_MPI */
    }

    phgInfo(3, "n        : %d\n", n);
    phgInfo(3, "nlocal   : %d\n", nlocal);
    phgInfo(3, "map_size : %d\n", map_size);
    phgInfo(3, "localsize: %d\n", localsize);

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
		struct {
		    BYTE a;
		    FLOAT b;
		} a;
		struct {
		    BYTE a;
		    INT b;
		} b;
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
	     *       is reversed. */
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
	    MPI_Alltoallv(temp, rcnts, rdsps, MPI_INT,
			  list, scnts, sdsps, MPI_INT, comm);
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
		    ncols =
			phgMatGetNnzInRow(solver->mat, list[i], NULL, NULL);
		if (ncols > 0)
		    nnz[rank] += (size_t) ncols;
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
		prows[1] = (INT)(nnz[rank] % (size_t) (INT_MAX));
		prows[2] = (INT)(nnz[rank] / (size_t) (INT_MAX));
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
					 ((size_t) stemp[rank]) * trunc_size);
			nrows = prows[0];
			size = ((size_t) prows[1]) +
			    ((size_t) prows[2]) * (size_t) (INT_MAX);

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
			memcpy(pdata, row->data, sizeof(FLOAT) * ncols);

			pcols += ncols;
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
	    rsize = ((size_t) k) * trunc_size;
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
		prows = (void *)(rbuff + ((size_t) rdsps[rank]) * trunc_size);
		nrows = prows[0];
		size = ((size_t) prows[1]) +
		    ((size_t) prows[2]) * (size_t) (INT_MAX);

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
			/* add a new global index to _t->O2Gmap */
			if (map_alloc <= map_size) {
			    _t->O2Gmap = phgRealloc_(_t->O2Gmap,
						     (map_alloc +
						      128) *
						     sizeof(*_t->O2Gmap),
						     map_alloc *
						     sizeof(*_t->O2Gmap));
			    _t->ordering =
				phgRealloc_(_t->ordering,
					    (map_alloc +
					     128) * sizeof(*_t->ordering),
					    map_alloc *
					    sizeof(*_t->ordering));
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
#if 1 || DEBUG
		    /* compare with the received row index */
		    r = *(prows++);
		    k = phgBinarySearchINT(localsize - nlocal, _t->O2Gmap,
					   ordering, r);
		    assert(k == r0);
#endif /* DEBUG */
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

	    phgInfo(3, "n        : %d\n", n);
	    phgInfo(3, "nlocal   : %d\n", nlocal);
	    phgInfo(3, "map_size : %d\n", map_size);
	    phgInfo(3, "localsize: %d %d\n", localsize, map_size + nlocal);

	    n = localsize;
	    localsize = map_size + nlocal;

	    phgFree(ordering);
	}

	phgFree(scnts);
    }
#endif /* USE_MPI */


    if (0 && solver->mat->refcount == 0)
	phgMatFreeMatrix(solver->mat);	/* mat is not destroyed */


    /*
     * Change the ordering of matrow columes index
     *   as is in matvec.c .
     *
     * */




    nlocal_o = n;

#if USE_MPI
    _t->cinfo = phgMapCreateCommInfo(solver->mat->rmap->comm,
				     solver->mat->rmap->nprocs,
				     solver->mat->rmap->rank,
				     nlocal,
				     localsize,
				     _t->O2Gmap,
				     _t->ordering,
				     solver->mat->rmap->partition);
#endif /* USE_MPI */

    _t->nlocal = nlocal;
    _t->nlocal_o = nlocal_o;
    _t->localsize = localsize;

    _t->rows = rows;


    /* ================================================================================
     *
     *
     *               Build matrix solver
     * 
     *
     * ================================================================================
     * */
    if (smoother_type == BOXVERT) {
	/* Extend block_indices from vert to vert box */
	INT *bidx0, *bidx_o;
	INT ncols, ib, size = 0, size_o = 0, alloc = 1000, alloc_o = 1000;

	phgInfo(0, "Extend for vert box.\n");

	bidx0 = _t->block_indices;	/* old boxes */
	bidx = phgAlloc(alloc * sizeof(*bidx));	/* new interior boxes */
	bidx_o = phgAlloc(alloc_o * sizeof(*bidx));	/* new interproc boxes */

	/* For parallel vert block */
	for (ib = 0; ib < nblock; ib++) {
	    int N = *(bidx0++), size0 = 0;
	    INT iG, iL, iG0;

	    phgInfo(3, "\n###vert: %d\n", ib);
	    for (i = 0; i < N; i++) {
		iG = bidx0[i];
		phgInfo(3, "   [%d] %6d\n", i, iG);

		if (iG >= n0 && iG < n0 + nlocal) {
		    iL = iG - n0;
		}
		else {
		    k = phgBinarySearchINT(localsize - nlocal, _t->O2Gmap,
					   _t->ordering, iG);
		    if (!(k < localsize - nlocal
			  && iG == _t->O2Gmap[_t->ordering[k]]))
			phgError(3, "  Block index is out of localsize.\n",
				 iG);
		    iL = nlocal + _t->ordering[k];
		}

		/* get matrix row */
		assert(iL < nlocal_o);
		row = _t->rows + iL;
		ncols = row->ncols;

		phgInfo(3, "   ncols: %d\n", ncols);

		if (ib < _t->nblock0) {	/* Note: indpendent of i */
		    /* add size ncols+1 */
		    while (alloc <= size + ncols) {
			phgInfo(3, "Realloc bidx: %d --> %d\n",
				alloc, size + ncols);
			bidx = phgRealloc_(bidx,
					   (alloc * 2) * sizeof(*bidx),
					   alloc * sizeof(*bidx));
			alloc *= 2;
		    }
		    if (i == 0) {	/* record size of box
					 * at first position */
			size0 = size++;
			bidx[size0] = 0;
		    }
		    /* convert row[i] col[j] to Global index,
		       to add to box */
		    for (j = 0; j < ncols; j++) {
			if ((k = row->cols[j]) < nlocal)
			    iG0 = n0 + k;
			else
			    iG0 = _t->O2Gmap[k - nlocal];

			for (l = size0 + 1; l < size; l++)
			    if (bidx[l] == iG0)
				break;
			if (l < size)	/* found dup */
			    continue;

			bidx[size0]++;	/* found new */
			bidx[size++] = iG0;
		    }
		    phgInfo(3, "   size add %d = %d\n",
			    size - size0 - 1, bidx[size0]);
		}
		else {
		    /* add size ncols+1 */
		    while (alloc_o <= size_o + ncols) {
			phgInfo(3, "Realloc bidx_o: %d --> %d\n",
				alloc_o, size_o + ncols);
			bidx_o = phgRealloc_(bidx_o,
					     (alloc_o * 2) * sizeof(*bidx),
					     alloc_o * sizeof(*bidx));
			alloc_o *= 2;
		    }
		    if (i == 0) {	/* record first position */
			size0 = size_o++;
			bidx_o[size0] = 0;
		    }
		    /* convert to Global index */
		    for (j = 0; j < ncols; j++) {
			if ((k = row->cols[j]) < nlocal)
			    iG0 = n0 + k;
			else
			    iG0 = _t->O2Gmap[k - nlocal];
			for (l = size0 + 1; l < size_o; l++)
			    if (bidx_o[l] == iG0)
				break;
			if (l < size_o)
			    continue;

			bidx_o[size0]++;
			bidx_o[size_o++] = iG0;
		    }
		}
	    }			/* end vert dofs */

	    bidx0 += N;
	}

	bidx = phgRealloc_(bidx,
			   (alloc + alloc_o) * sizeof(*bidx),
			   alloc * sizeof(*bidx));
	memcpy(bidx + size, bidx_o, size_o * sizeof(*bidx));
	block_indices = bidx;
	phgFree(_t->block_indices);
	_t->block_indices = bidx;
    }


    /* Count cache size */
    npi = 0, npd = 0;
    npi_o = 0, npd_o = 0;
    bidx = block_indices;
    maxN = 0, minN = 100000;

    /* Check */
    if (FALSE)
	for (bidx = block_indices, ib = 0; ib < nblock; ib++) {
	    INT N = *(bidx++);
	    phgInfo(0, "block: %d [%d]\n", ib, N);
	    for (j = 0; j < N; j++) {
		phgInfo(0, "%6d, ", bidx[j]);
	    }
	    phgInfo(0, "\n");
	    bidx += N;
	}

    bidx = block_indices;
    for (ib = 0; ib < nblock; ib++) {
	if (ib == _t->nblock0)
	    _t->bdry_offset = bidx - block_indices;

	{
	    INT N = *(bidx++);
	    INT iG[N], iL[N];

	    if (maxN < N)
		maxN = N;
	    if (minN > N)
		minN = N;

	    memcpy(iG, bidx, N * sizeof(*bidx));
	    bidx += N;

	    /* search local row # */
	    for (i = 0; i < N; i++) {
		if (iG[i] >= n0 && iG[i] < n0 + nlocal) {
		    iL[i] = iG[i] - n0;
		}
		else {
		    k = phgBinarySearchINT(localsize - nlocal, _t->O2Gmap,
					   _t->ordering, iG[i]);
		    if (!(k < localsize - nlocal
			  && iG[i] == _t->O2Gmap[_t->ordering[k]]))
			phgError(3, "  Block index is out of localsize.\n",
				 iG[i]);
		    iL[i] = nlocal + _t->ordering[k];
		}

		assert(iL[i] < nlocal_o);
	    }

	    for (i = 0; i < N; i++) {
		INT n1 = 0, n2 = 0;

		row = _t->rows + iL[i];
		for (j = 0; j < row->ncols; j++) {
		    INT jcol = row->cols[j];

		    for (k = 0; k < N; k++)
			if (jcol == iL[k])
			    break;
		    if (k < N) {
			continue;
		    }

		    if (jcol < nlocal) {
			if (!(row->data[j] == zero))
			    n1++;	/* record */
			continue;
		    }

		    if (!(row->data[j] == zero))
			n2++;	/* record */
		}
		assert(n1 < 10000);
		assert(n2 < 10000);

		npi++;		/* ncol */
		npi += n1;
		npd += n1;

		npi_o++;	/* ncol */
		npi_o += n2;
		npd_o += n2;
	    }

	    if (FT_PHG != FT_DOUBLE || DENSE_INTERNAL) {
		/* internal */
		npd += N * N;	/* record matrix */
		npi += N;	/* pivot */
	    }
	    else if (DENSE_LAPACK) {
		/* Lapack */
		npd += N * N;
	    }
	    npi += N;		/* vec index */
	}
    }

    phgInfo(1, "Block size: [%d %d]\n", minN, maxN);
    phgInfo(1, "* setup time: %lg\n", phgGetTime(NULL) - time0);
    _t->maxN = maxN;

    /*
     * Packed data
     *
     * 
     * 0. # of blocks
     * 0'. block entries for idx and data
     * 1. matrix A[N][N]
     * 2. Vec index I[N], to form b[I[.]]
     * 3. other columes in row
     *    col_i[...] to access x[...], like pc_o
     *    dat_i[...] to access a[][], like pd_o
     *
     * */

    time0 = phgGetTime(NULL);

    packed_col = phgCalloc(4 * nblock + npi + npi_o, sizeof(*pc));
    _t->packed_col = packed_col;
    pc_b = packed_col;		/* block entries */
    pc = pc_b + 4 * nblock;	/* diag idx */
    pc_o = pc + npi;		/* offp idx, and possible pivot */

    packed_dat = phgCalloc(npd + npd_o, sizeof(*pd));
    _t->packed_dat = packed_dat;
    pd = packed_dat;
    pd_o = pd + npd;


    _t->bk_solvers = phgCalloc(nblock, sizeof(*_t->bk_solvers));
    _t->bk_vec = phgCalloc(nblock, sizeof(*_t->bk_vec));
    A = phgCalloc(MAT_MAX_SIZE * MAT_MAX_SIZE, sizeof(*A));

    col = NULL, allocated = 0;
    dat = NULL;

    bidx = block_indices;
    for (ib = 0; ib < nblock; ib++) {
	INT N = *(bidx++);
	INT iG[N], iL[N];
	MAP *map = NULL;
	MAT *B = NULL;
	int pivot[N];

	memcpy(iG, bidx, N * sizeof(*bidx));
	bidx += N;

	pc_b[4 * ib] = pc - packed_col;
	pc_b[4 * ib + 1] = pd - packed_dat;
	pc_b[4 * ib + 2] = pc_o - packed_col;
	pc_b[4 * ib + 3] = pd_o - packed_dat;

	/* search local row # */
	for (i = 0; i < N; i++) {
	    if (iG[i] >= n0 && iG[i] < n0 + nlocal) {
		iL[i] = iG[i] - n0;
	    }
	    else {
		k = phgBinarySearchINT(localsize - nlocal, _t->O2Gmap,
				       _t->ordering, iG[i]);
		if (!(k < localsize - nlocal
		      && iG[i] == _t->O2Gmap[_t->ordering[k]]))
		    phgError(3, "  Block index %3d is out of localsize.\n",
			     iG[i]);
		iL[i] = nlocal + _t->ordering[k];
	    }

	    assert(iL[i] < nlocal_o);
	}


	/* Other direct solver */
	if (FT_PHG != FT_DOUBLE || DENSE_INTERNAL) {
	    /* internal solver */
	    bzero(A, N * N * sizeof(*A));
	}
	else if (DENSE_LAPACK) {
	    /* Lapack */
	    bzero(A, N * N * sizeof(*A));
	}
	else {
	    /* external solver */
	    map = phgMapCreateSimpleMap(MPI_COMM_SELF, N, N);
	    B = phgMapCreateMat(map, map);
	    phgMapDestroy(&map);
	}

	for (i = 0; i < N; i++) {
	    INT n1 = 0, n2 = 0;
	    row = _t->rows + iL[i];
	    assert(row->ncols > 0);

	    if (allocated < row->ncols) {
		allocated = row->ncols;
		phgFree(col);
		col = phgAlloc(2 * allocated * sizeof(*col));
		col_o = col + allocated;
		phgFree(dat);
		dat = phgAlloc(2 * allocated * sizeof(*dat));
		dat_o = dat + allocated;
	    }

	    for (j = 0; j < row->ncols; j++) {
		INT jcol = row->cols[j];

		/* matrix */
		/* Case 1)
		 *    In block
		 *     save dat in A,
		 *     !save col 
		 * */
		for (k = 0; k < N; k++)
		    if (jcol == iL[k])
			break;
		if (k < N) {
		    if (FT_PHG != FT_DOUBLE || DENSE_INTERNAL) {
			/* internal solver */
			A[i * N + k] = row->data[j];
		    }
		    else if (DENSE_LAPACK) {
			/* Lapack */
			A[i * N + k] = row->data[j];
		    }
		    else {
			/* external solver */
			if (!(row->data[j] == zero))
			    phgMatAddEntry(B, i, k, row->data[j]);
		    }
		    continue;
		}

		/* Case 2)
		 *   Out block, local
		 *    save dat in pdat
		 *    save col 
		 * */
		if (jcol < nlocal) {
		    if (!(row->data[j] == zero)) {
			col[n1] = jcol;	/* vec data */
			dat[n1] = row->data[j];
			n1++;
		    }
		    continue;
		}

		/* Case 3)
		 *   Out block, offp
		 *    save dat in pdat
		 *    save col 
		 * */
		assert(jcol < localsize);
		if (!(row->data[j] == zero)) {
		    col_o[n2] = jcol;	/* vec data */
		    dat_o[n2] = row->data[j];
		    n2++;
		}
	    }

#if 0
#warning Note ZLB: remove the next two lines?
	    assert(n1 < 10000);
	    assert(n2 < 10000);
#endif

	    *(pc++) = n1;
	    memcpy(pc, col, n1 * sizeof(*pc));
	    memcpy(pd, dat, n1 * sizeof(*pd));
	    pc += n1;
	    pd += n1;

	    *(pc_o++) = n2;
	    memcpy(pc_o, col_o, n2 * sizeof(*pc_o));
	    memcpy(pd_o, dat_o, n2 * sizeof(*pd_o));
	    pc_o += n2;
	    pd_o += n2;
	}

	if (FT_PHG != FT_DOUBLE || DENSE_INTERNAL) {
	    /* phg internal solver */
	    phgSolverDenseLU(N, A, pivot);
	    memcpy(pd, A, N * N * sizeof(*pd));
	    pd += N * N;

	}
	else if (DENSE_LAPACK) {
	    if (FALSE) {
		/* Check matrix */
		phgInfo(0, "elem mat[%d]\n", ib);
		for (i = 0; i < N; i++) {
		    for (j = 0; j < N; j++) {
			phgInfo(0, "%8.3e ", A[i * N + j]);
		    }
		    phgInfo(0, "\n");
		}
	    }

	    /* Lapack */
	    {
		int OK = 0, LDA = N, M = N, n = N, LWORK = N * 3;
		double WORK[LWORK];

		assert(sizeof(double) == sizeof(FLOAT));
		DGETRF_FOR(&M, &n, (double *)A, &LDA, pivot, &OK);
		assert(OK == 0);
		DGETRI_FOR(&n, (double *)A, &LDA, pivot, WORK, &LWORK, &OK);
		assert(OK == 0);
	    }

	    memcpy(pd, A, N * N * sizeof(*pd));
	    pd += N * N;
	}
	else {
	    /* external dense solver */
	    phgOptionsPush();
	    phgOptionsSetKeyword("-solver", "mumps");
	    phgOptionsSetOptions(_t->solver_opts);
	    _t->bk_solvers[ib] = phgMat2Solver(SOLVER_DEFAULT, B);
	    _t->bk_solvers[ib]->monitor = FALSE;
	    phgOptionsPop();
	    phgMatDestroy(&B);

	    phgSolverAssemble(_t->bk_solvers[ib]);
	    _t->bk_vec[ib] = phgMapCreateVec(_t->bk_solvers[ib]->rhs->map, 1);
	}

	memcpy(pc, iL, N * sizeof(*pc));
	pc += N;

	if (FT_PHG != FT_DOUBLE || DENSE_INTERNAL) {
	    /* internal solver, save pivot */
	    for (i = 0; i < N; i++)
		pc[i] = pivot[i];
	    pc += N;
	}
    }

    phgFree(A);
    phgFree(col);
    phgFree(dat);

    phgInfo(1, "* factorize time: %lg\n", phgGetTime(NULL) - time0);

    /* cache */
    _t->x0_dat = phgCalloc(localsize, sizeof(FLOAT));
    _t->rhs_dat = phgCalloc(localsize, sizeof(FLOAT));

    return 0;
}


static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    /*double time0 = 0.; */
    INT i, j, iter, ib, *bidx;
    INT *packed_col = _t->packed_col;
    FLOAT *packed_dat = _t->packed_dat;
    INT nlocal = _t->nlocal;
    INT nlocal_o = _t->nlocal_o;
    INT localsize = _t->localsize;
    /*INT n0 = solver->rhs->map->partition[solver->rhs->map->rank]; */
    FLOAT *x0_dat = _t->x0_dat;
    FLOAT *rhs_dat = _t->rhs_dat;
    INT *pc, *pc_o, *pc_b = packed_col;
    FLOAT *pd, *pd_o;
    int maxN = _t->maxN;
    int nblock = _t->nblock;
    /*int nblock0 = _t->nblock0;
       int nblock1 = _t->nblock1;
       int bdry_offset = _t->bdry_offset; */
    INT *block_indices = _t->block_indices;
    FLOAT B[maxN], X[maxN];
    VEC *res;
    BOOLEAN print_residual = FALSE;

    assert(x->nvec == 1);
    Assemble(solver);

    /*if (phgVerbosity > 1)
       time0 = phgGetTime(NULL); */

    if (nlocal > 0) {
	memcpy(rhs_dat, solver->rhs->data, nlocal * sizeof(*rhs_dat));
	memcpy(x0_dat, x->data, nlocal * sizeof(*x0_dat));
    }

    if (print_residual)
	res = phgVecCopy(solver->rhs, NULL);

#if USE_MPI
    if (_t->overlap > 0 && solver->mat->cmap->nprocs > 1) {
	phgMapScatter(_t->cinfo,
		      solver->rhs->nvec, rhs_dat, rhs_dat + nlocal);
    }
#endif /* USE_MPI */

    /* Check matrix */
    if (FALSE) {
	MAT_ROW *row;

	phgInfo(0, "\n=================\n");
	phgInfo(0, "nlocal  : %d\n", nlocal);
	phgInfo(0, "nlocal_o: %d\n", nlocal_o);
	phgInfo(0, "localsize: %d\n", localsize);


	/* matrix data:
	     rows[0, nlocal_o]
	     colume entry
	*/
	for (i = 0; i < nlocal_o; i++) {
	    INT k;
	    FLOAT res = 0;
	    row = _t->rows + i;
	    phgInfo(0, "row [%5d]\n", i);

	    res = rhs_dat[i];

	    for (j = 0; j < row->ncols; j++) {
		if ((k = row->cols[j]) < nlocal) {
		    res -= x0_dat[k] * row->data[j];
		    continue;
		}
		res -= x0_dat[k]
		    * row->data[j];
	    }
	    phgInfo(0, "   res: %e\n", res);
	}
    }


    for (iter = 0; iter < solver->maxit; iter++) {

	/* ------------------------------------------------------------
	 * Smoother
	 *    Gauss-Seidel type
	 * ------------------------------------------------------------ */


	int k;
	for (k = 0; k < 2; k++) {
	    INT ib0 = 0, ib1 = 0, di = 1;
	    bidx = block_indices;
	    pc_b = packed_col;
	    switch (k) {	/* forward
				 * note: first inner then interproc */
		case (0):
		    ib0 = 0;
		    ib1 = nblock;
		    di = 1;
		    break;
		case (2):	/* backward */
		    ib0 = nblock;
		    ib1 = 0;
		    di = -1;
		    break;
	    }

	    for (ib = ib0; ib < ib1; ib += di) {
		/* int N; */
		/* bidx = block_indices + _t->block_ofs[ib]; */
		/* N = *(bidx++); */
		INT N = *(bidx++);
		INT *iL;
		bidx += N;

		pc = packed_col + pc_b[4 * ib];
		pd = packed_dat + pc_b[4 * ib + 1];
		pc_o = packed_col + pc_b[4 * ib + 2];
		pd_o = packed_dat + pc_b[4 * ib + 3];

		/* Build rhs */
		bzero(B, maxN * sizeof(*B));
		bzero(X, maxN * sizeof(*X));

		for (i = 0; i < N; i++) {
		    INT ncol = *(pc++);
		    INT ncol_o = *(pc_o++);
		    INT jcol;

		    /* local */
		    for (j = 0; j < ncol; j++) {
			jcol = *(pc++);
			assert(jcol >= 0 && jcol < localsize);
			B[i] -= *(pd++) * x0_dat[jcol];
		    }

		    /* offp */
		    for (j = 0; j < ncol_o; j++) {
			jcol = *(pc_o++);
			assert(jcol >= 0 && jcol < localsize);
			B[i] -= *(pd_o++) * x0_dat[jcol];
		    }
		}

		iL = pc;
		pc += N;

		for (i = 0; i < N; i++) {
		    B[i] += rhs_dat[iL[i]];
		}

		/* Solve */
		if (FT_PHG != FT_DOUBLE || DENSE_INTERNAL) {
		    /* phg internal solver */
		    FLOAT *A = pd;
		    int pivot[N];

		    memcpy(X, B, N * sizeof(*X));
		    phgSolverDenseSV(N, A, pivot, 1, X);
		    for (i = 0; i < N; i++)
			pc[i] = pivot[i];
		    pc += N;
		    pd += N * N;

		}
		else if (DENSE_LAPACK) {
		    /* Lapack */
		    char TRANSA = 'N', TRANSB = 'N';
		    int M = 1, K = N, LDA = 1, LDB = N, n = N, LDC = 1;
		    double ALPHA = 1., BETA = 0;

		    FLOAT *A = pd;
		    pd += N * N;

		    DGEMM_FOR(&TRANSA, &TRANSB, &M, &n, &K, &ALPHA,
			      (double *)B, &LDA, (double *)A, &LDB, &BETA,
			      (double *)X, &LDC);
		}
		else {
		    SOLVER *bk_solver = _t->bk_solvers[ib];
		    VEC *vec_x = _t->bk_vec[ib];
		    double time_solve = phgGetTime(NULL);
		    /*size_t mem, mem_peak; */

		    memcpy(bk_solver->rhs->data, B, N * sizeof(*B));
		    bzero(vec_x->data, N * sizeof(vec_x->data));
		    bk_solver->rhs->assembled = TRUE;
		    bk_solver->rhs_updated = TRUE;
		    phgSolverVecSolve(bk_solver, FALSE, vec_x);
		    memcpy(X, vec_x->data, N * sizeof(*X));

		    phgInfo(3, "  %6d[%6d] time: %lg, mem: %lg\n",
			    ib, nblock, phgGetTime(NULL) - time_solve, 0);
		}

		for (i = 0; i < N; i++) {
		    FLOAT x0, dx;

		    x0 = x0_dat[iL[i]];
		    dx = X[i] - x0;
#if 0
		    phgInfo(0, "x0: %12.4e dx[%3d]: %12.4e\n", x0, iL[i], dx);
		    assert(fabs(dx) < 1e-12);
#else
		    x0_dat[iL[i]] = x0 + omega * dx;
#endif
		}
	    }

#if USE_MPI
	    if (_t->overlap > 0 && solver->mat->cmap->nprocs > 1) {
		phgMapScatter(_t->cinfo,
			      solver->rhs->nvec, x0_dat, x0_dat + nlocal);
	    }
#endif /* USE_MPI */
	}

	memcpy(x->data, x0_dat, nlocal * sizeof(*x0_dat));

	if (print_residual) {
	    phgVecCopy(solver->rhs, &res);
	    phgMatVec(MAT_OP_N, -1., solver->mat, x, 1., &res);
	    phgPrintf("   res: %e\n", phgVecNorm2(res, 0, NULL));
	}
    }

    if (print_residual)
	phgVecDestroy(&res);
    /*phgVecDumpMATLAB(x, "x", "x.m"); */
    return solver->nits;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverSmoother_ = {
    "smoother", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy,
    AddMatrixEntries, AddRHSEntries,
    Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, FALSE, FALSE
};
