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

/* $Id: solver-spooles.c,v 1.66 2022/09/21 02:35:28 zlb Exp $ */

/* The SPOOLES interface */

#include "phg.h"

#if USE_SPOOLES	/* whole file */

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

#include <spoolesMPI.h>

typedef struct {
    InpMtx	*mtxA;
    FrontMtx	*frontmtx;
    ETree	*frontETree;
    SolveMap	*solvemap;
    IV		*vtxmapIV, *ownersIV, *oldToNewIV, *newToOldIV;
    IV		*rowmapIV, *ownedColumnsIV;
    IVL         *symbfacIVL;
    SubMtxManager *mtxmanager;

    double	cpus[20];
    int		stats[20];

    BOOLEAN	factored;

    int		lookahead;
    int		pivoting;
} OEM_DATA;

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)

static INT lookahead = 0;	/* look ahead parameter */
static const char *pivoting_list[] = {"default", "true", "false", NULL};
enum {PIVOTING_DEFAULT, PIVOTING_TRUE, PIVOTING_FALSE};
static int pivoting = PIVOTING_DEFAULT;

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nSPOOLES options:", "\n", "spooles");
    phgOptionsRegisterKeyword("spooles_pivoting", "Pivoting",
				pivoting_list, &pivoting);
    phgOptionsRegisterInt("spooles_lookahead", "Look ahead parameter",
				&lookahead);
    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));

    _t->lookahead	= lookahead;
    _t->pivoting	= pivoting;

    return 0;
}

static int
Create(SOLVER *solver)
{
    phgInfo(2, "create SPOOLES solver.\n");

    /* avoid GCC warnings on TZ/TV */
    Unused(TV);
    Unused(TZ);

    IVzero(20, _t->stats);
    DVzero(20, _t->cpus);

    return 0;
}

static void
destroy_objects(SOLVER *solver)
{
    if (_t->mtxA != NULL) {
	InpMtx_free(_t->mtxA);
	_t->mtxA = NULL;
    }
    if (_t->frontmtx != NULL) {
	FrontMtx_free(_t->frontmtx);
	_t->frontmtx = NULL;
    }
    if (_t->mtxmanager != NULL) {
	SubMtxManager_free(_t->mtxmanager);
	_t->mtxmanager = NULL;
    }
    if (_t->frontETree != NULL) {
	ETree_free(_t->frontETree);
	_t->frontETree = NULL;
    }
    if (_t->solvemap != NULL) {
	SolveMap_free(_t->solvemap);
	_t->solvemap = NULL;
    }
    if (_t->vtxmapIV != NULL) {
	IV_free(_t->vtxmapIV);
	_t->vtxmapIV = NULL;
    }
    if (_t->ownersIV != NULL) {
	IV_free(_t->ownersIV);
	_t->ownersIV = NULL;
    }
    if (_t->oldToNewIV != NULL) {
	IV_free(_t->oldToNewIV);
	_t->oldToNewIV = NULL;
    }
    if (_t->newToOldIV != NULL) {
	 IV_free(_t->newToOldIV);
	_t->newToOldIV = NULL;
    }
    if (_t->rowmapIV != NULL) {
	IV_free(_t->rowmapIV);
	_t->rowmapIV = NULL;
    }
    if (_t->ownedColumnsIV != NULL) {
	IV_free(_t->ownedColumnsIV);
	_t->ownedColumnsIV = NULL;
    }
    if (_t->symbfacIVL != NULL) {
	 IVL_free(_t->symbfacIVL);
	_t->symbfacIVL = NULL;
    }
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    phgInfo(2, "destroy SPOOLES objects.\n");
    destroy_objects(solver);
    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    MAP *map = solver->mat->rmap;
    int i, j, k, nz, msglvl = 0, firsttag = 0;
    int Nlocal = map->nlocal, n0 = map->partition[map->rank];
    const MAT_ROW *row;
    double time0, d, tau;
    int symmetryflag = (solver->symmetry >= M_SYM ? SPOOLES_SYMMETRIC :
						    SPOOLES_NONSYMMETRIC);
    /* use pivoting when the matrix is not SPD */
    int pivotingflag = SPOOLES_NO_PIVOTING;
    int seed = 10101;

    /* Note: MPI_TAG_UB is undefined if comm == MPI_COMM_SELF */
    int *iptr, flag;
    MPI_Errhandler errhandler;
    MPI_Comm comm = map->comm;

    if (_t->factored)
	return 0;
    _t->factored = TRUE;

    /* Note: SPOOLES gets MPI_TAG_UB value by calling MPI_Attr_get, but
     *	     on some MPI systems MPI_COMM_SELF does not have this attribute,
     *	     the code below sets the attribute MPI_TAG_UB to MPI_COMM_SELF
     *	     for use by SPOOLES.
     *
     *	     Moreover, some OpenMPI versions (including 1.3) fail on the
     *	     MPI_Attr_put call below, for these MPI systems we need a modified
     *	     SPOOLES version in which the attribute value of MPI_TAG_UB is
     *	     obtained from MPI_COMM_WORLD instead of 'comm', and we temporarily
     *	     change MPI error handler to MPI_ERRORS_RETURN to not abort. */
    MPI_Comm_get_errhandler(comm, &errhandler);
    MPI_Comm_set_errhandler(comm, MPI_ERRORS_RETURN);
    MPI_Comm_get_attr(comm, MPI_TAG_UB, (void *)&iptr, &flag);
    if (flag == 0) {
	static int tag_bound = -1;
	phgInfo(1, "setting MPI_TAG_UB for map->comm\n");
	if (tag_bound == -1) {
	    MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, (void *)&iptr, &flag);
	    tag_bound = *iptr;
	}
	MPI_Comm_set_attr(comm, MPI_TAG_UB, &tag_bound);
    }
    MPI_Comm_set_errhandler(comm, errhandler);

    assert(solver->mat != NULL);
    if (solver->sent_to_oem == TRUE)
	phgError(1, "%s:%d, data already sent to OEM solver!\n",
						__FILE__, __LINE__);
    solver->sent_to_oem = TRUE;

    assert(solver->mat->type != PHG_DESTROYED);
    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }

    switch (_t->pivoting) {
	case PIVOTING_DEFAULT:
	    pivotingflag = (solver->symmetry == M_SPD ? SPOOLES_NO_PIVOTING :
							SPOOLES_PIVOTING);
	    break;
	case PIVOTING_TRUE:
	    pivotingflag = SPOOLES_PIVOTING;
	    break;
	case PIVOTING_FALSE:
	    pivotingflag = SPOOLES_NO_PIVOTING;
	    break;
    }

    /* FIXME: incorrect result if (np == 1 && SYMMETRIC && PIVOTING) */
    if (map->nprocs == 1 && symmetryflag == SPOOLES_SYMMETRIC &&
	pivotingflag == SPOOLES_PIVOTING) {
#if 0
	static BOOLEAN warned = FALSE;
	if (!warned) {
	    warned = TRUE;
	    phgInfo(-1, "***\n");
	    phgWarning("the result may be incorrect if nprocs == 1\n");
	    phgWarning("with -solver_symmetry=sym\n");
	    phgInfo(-1, "***\n");
	}
#else
	/*phgWarning("nprocs == 1: force unsymmetric mode.\n");
	symmetryflag = SPOOLES_NONSYMMETRIC;*/
#endif
    }

    phgInfo(2, "create SPOOLES input matrix.\n");

    time0 = phgGetTime(NULL);
    /* count # of nonzeros */
    nz = 0;
    tau = 1.0;
    if (symmetryflag == SPOOLES_NONSYMMETRIC) {
	for (i = 0; i < Nlocal; i++) {
	    row = phgMatGetRow(solver->mat, i);
	    nz += row->ncols;
	    for (j = 0; j < row->ncols; j++) {
		if ((d = Fabs(row->data[j])) > tau)
		    tau = d;
	    }
	}
    }
    else {
	for (i = 0; i < Nlocal; i++) {
	    row = phgMatGetRow(solver->mat, i);
	    for (j = 0; j < row->ncols; j++) {
		k = row->cols[j];
		if (k < i + n0)
		    continue;
		nz++;
		if ((d = Fabs(row->data[j])) > tau)
		    tau = d;
	    }
	}
    }

    MPI_Allreduce(&tau, &d, 1, MPI_DOUBLE, MPI_MAX, comm);
    tau = d;

    _t->mtxA = InpMtx_new();
    InpMtx_init(_t->mtxA, INPMTX_BY_ROWS, SPOOLES_REAL, nz, 0);

    for (i = 0; i < Nlocal; i++) {
	row = phgMatGetRow(solver->mat, i);
	for (j = 0; j < row->ncols; j++) {
	    k = row->cols[j];
	    if (symmetryflag == SPOOLES_SYMMETRIC && k < i + n0)
		continue;
	    InpMtx_inputRealEntry(_t->mtxA, i + n0, k, row->data[j]);
	}
	/* free matrix row in solver */
	if (solver->mat->refcount == 0 && solver->mat->rows != NULL) {
	    phgFree(solver->mat->rows[i].cols);
	    phgFree(solver->mat->rows[i].data);
	    solver->mat->rows[i].cols = NULL;
	    solver->mat->rows[i].data = NULL;
	    solver->mat->rows[i].ncols = solver->mat->rows[i].alloc = 0;
	}
    }
    InpMtx_sortAndCompress(_t->mtxA);
    InpMtx_changeStorageMode(_t->mtxA, INPMTX_BY_VECTORS);

    if (solver->mat->refcount == 0)
	phgMatFreeMatrix(solver->mat);

    phgInfo(1, "Time for setting up input matrix: %lg\n",
			(double)(phgGetTime(NULL) - time0));

    phgInfo(2, "analyse and factor matrix.\n");

    /* find a low-fill ordering */
    {
	int nedges, root;
	double *opcounts, minops;
	Graph *graph;
	IVL *adjIVL;

	graph = Graph_new();
	adjIVL = InpMtx_MPI_fullAdjacency(_t->mtxA, _t->stats,
					  msglvl, stdout, comm);
	nedges = IVL_tsize(adjIVL);
	Graph_init2(graph, 0, map->nglobal, 0, nedges, map->nglobal,
		    nedges, adjIVL, NULL, NULL);
	_t->frontETree = orderViaMMD(graph, seed + map->rank,
				     msglvl, stdout);
	Graph_free(graph);
	opcounts = DVinit(map->nprocs, 0.0);
	minops = ETree_nFactorOps(_t->frontETree, SPOOLES_REAL, symmetryflag);
	MPI_Allgather((void *)&minops, 1, MPI_DOUBLE,
		      (void *)opcounts, 1, MPI_DOUBLE, comm);
	minops = DVmin(map->nprocs, opcounts, &root);
	DVfree(opcounts);
	_t->frontETree = ETree_MPI_Bcast(_t->frontETree, root,
					 msglvl, stdout, comm);
    }

    /* get the permutations, permute the front tree and matrix */
    {
	_t->oldToNewIV = ETree_oldToNewVtxPerm(_t->frontETree);
	_t->newToOldIV = ETree_newToOldVtxPerm(_t->frontETree);
	ETree_permuteVertices(_t->frontETree, _t->oldToNewIV);
	InpMtx_permute(_t->mtxA, IV_entries(_t->oldToNewIV),
				 IV_entries(_t->oldToNewIV));
	if (symmetryflag == SPOOLES_SYMMETRIC 
	    || symmetryflag == SPOOLES_HERMITIAN) {
	    InpMtx_mapToUpperTriangle(_t->mtxA);
	}
	InpMtx_changeCoordType(_t->mtxA, INPMTX_BY_CHEVRONS);
	InpMtx_changeStorageMode(_t->mtxA, INPMTX_BY_VECTORS);
    }

    /* generate the owners map IV object and the map from vertices to owners */
    {
	double cutoff = 1. / (2. * map->nprocs);
	DV *cumopsDV = DV_new();
	DV_init(cumopsDV, map->nprocs, NULL);
	_t->ownersIV = ETree_ddMap(_t->frontETree, SPOOLES_REAL, symmetryflag,
				   cumopsDV, cutoff);
	DV_free(cumopsDV);
	_t->vtxmapIV = IV_new();
	IV_init(_t->vtxmapIV, map->nglobal, NULL);
	IVgather(map->nglobal, IV_entries(_t->vtxmapIV),
		 IV_entries(_t->ownersIV), ETree_vtxToFront(_t->frontETree));
    }

    /* redistribute the matrix */
    {
	InpMtx *newA;
	newA = InpMtx_MPI_split(_t->mtxA, _t->vtxmapIV, _t->stats, msglvl,
				stdout, firsttag, comm);
	InpMtx_free(_t->mtxA);
	_t->mtxA = newA ;
	InpMtx_changeStorageMode(_t->mtxA, INPMTX_BY_VECTORS);
    }

    /* compute the symbolic factorization */
    _t->symbfacIVL = SymbFac_MPI_initFromInpMtx(_t->frontETree, _t->ownersIV,
						_t->mtxA, _t->stats, msglvl,
						stdout, firsttag, comm);

    /* initialize the front matrix */
    /* Note: mtxmanager is needed by frontmtx, so we save and free it later */
    _t->mtxmanager = SubMtxManager_new();
    SubMtxManager_init(_t->mtxmanager, NO_LOCK, 0);
    _t->frontmtx = FrontMtx_new();
    FrontMtx_init(_t->frontmtx, _t->frontETree, _t->symbfacIVL,
		  SPOOLES_REAL, symmetryflag, FRONTMTX_DENSE_FRONTS,
		  pivotingflag, NO_LOCK, map->rank, _t->ownersIV,
		  _t->mtxmanager, msglvl, stdout);

    /* compute the factorization */
    {
	ChvManager *chvmanager;
	Chv *rootchv;
	int error;
	static double droptol = 0.;
	chvmanager = ChvManager_new();
	ChvManager_init(chvmanager, NO_LOCK, 0);
	rootchv = FrontMtx_MPI_factorInpMtx(_t->frontmtx, _t->mtxA, tau,
					    droptol, chvmanager, _t->ownersIV,
					    _t->lookahead, &error, _t->cpus,
					    _t->stats, msglvl, stdout,
					    firsttag, comm);
	Unused(rootchv);
	ChvManager_free(chvmanager);
	if (error >= 0)
	    phgError(1, "factorization error at front %d", error);
    }

    /* post-process the factorization
     * and split the factor matrices into submatrices */
    FrontMtx_MPI_postProcess(_t->frontmtx, _t->ownersIV, _t->stats, msglvl,
			     stdout, firsttag, comm);

    /* create the solve map object */
    _t->solvemap = SolveMap_new();
    SolveMap_ddMap(_t->solvemap, _t->frontmtx->symmetryflag,
		   FrontMtx_upperBlockIVL(_t->frontmtx),
		   FrontMtx_lowerBlockIVL(_t->frontmtx),
		   map->nprocs, _t->ownersIV,
		   FrontMtx_frontTree(_t->frontmtx),
		   seed, msglvl, stdout);

    /* redistribute the submatrices of the factors */
    FrontMtx_MPI_split(_t->frontmtx, _t->solvemap, _t->stats, msglvl, stdout,
		       firsttag, comm);

    if (FRONTMTX_IS_PIVOTING(_t->frontmtx))
	_t->rowmapIV = FrontMtx_MPI_rowmapIV(_t->frontmtx, _t->ownersIV, msglvl,
					     stdout, comm);

    _t->ownedColumnsIV = FrontMtx_ownedColumnsIV(_t->frontmtx, map->rank,
						 _t->ownersIV, msglvl, stdout);

    return 0;
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    MAP *map = solver->mat->rmap;
    int i, Nlocal = map->nlocal, firsttag = 0, msglvl = 0;
    int nmycol, nrow, ncol, *rowind, *colind;
    double *entries;
    INT n0 = map->partition[map->rank];
    DenseMtx *mtxX, *mtxY, *mtxT;
    SubMtxManager *solvemanager;
    MPI_Comm comm = map->comm;

    Assemble(solver);

    /* create dense matrix for the RHS (Y) */
    mtxY = DenseMtx_new();
    DenseMtx_init(mtxY, SPOOLES_REAL, 0, 0, Nlocal, 1, 1, Nlocal);

    /* set up RHS */
    DenseMtx_rowIndices(mtxY, &nrow, &rowind);
    DenseMtx_columnIndices(mtxY, &ncol, &colind);
    colind[0] = 0;
    entries = DenseMtx_entries(mtxY);
    for (i = 0; i < nrow; i++) {
	rowind[i] = i + n0;
#if 0
	DenseMtx_setRealEntry(mtxY, i, 0, solver->rhs[i]);
#else
	*(entries++) = solver->rhs->data[i];
#endif
    }

    /* permute the RHS */
    DenseMtx_permuteRows(mtxY, _t->oldToNewIV);

    /* redistribute the right hand side */
    mtxT = DenseMtx_MPI_splitByRows(mtxY, _t->vtxmapIV, _t->stats, msglvl,
				    stdout, firsttag, comm);
    DenseMtx_free(mtxY);
    mtxY = mtxT;

    /* permute and redistribute Y if necessary */
    if (FRONTMTX_IS_PIVOTING(_t->frontmtx)) {
	mtxT = DenseMtx_MPI_splitByRows(mtxY, _t->rowmapIV, _t->stats, msglvl,
					stdout, firsttag, comm);
	DenseMtx_free(mtxY);
	mtxY = mtxT;
    }

    /* create a solution DenseMtx object */
    nmycol = IV_size(_t->ownedColumnsIV);
    mtxX = DenseMtx_new();
    if (nmycol > 0) {
	DenseMtx_init(mtxX, SPOOLES_REAL, 0, 0, nmycol, 1, 1, nmycol);
	DenseMtx_rowIndices(mtxX, &nrow, &rowind);
	IVcopy(nmycol, rowind, IV_entries(_t->ownedColumnsIV));
    }

    /* solve the linear system */
    solvemanager = SubMtxManager_new();
    SubMtxManager_init(solvemanager, NO_LOCK, 0);
    FrontMtx_MPI_solve(_t->frontmtx, mtxX, mtxY, solvemanager, _t->solvemap,
		       _t->cpus, _t->stats, msglvl, stdout, firsttag, comm);
    SubMtxManager_free(solvemanager);
    DenseMtx_free(mtxY);

/*static int iii = 0;
if (++iii == 3 && map->rank == 0) {
IV_writeForHumanEye(_t->newToOldIV, stdout);
DenseMtx_writeForHumanEye(mtxX, stdout);
printf("\n");
}*/

    /* permute the solution into the original ordering */
#if 0
    DenseMtx_permuteRows(mtxX, _t->newToOldIV);
#else
    /* FIXME: DenseMtx_permuteRows returns wrong result if
     *    nprocs == 1 && SPOOLES_SYMMETRY && SPOOLES_PIVOTING
     * Below is a workaround */
    if (map->nprocs > 1) {
	DenseMtx_permuteRows(mtxX, _t->newToOldIV);
    }
    else {
	double *buffer = phgAlloc(map->nglobal * sizeof(*buffer));
	memcpy(buffer, mtxX->entries, map->nglobal * sizeof(*buffer));
	colind = IV_entries(_t->newToOldIV);
	DenseMtx_rowIndices(mtxX, &nrow, &rowind);
	for (i = 0; i < map->nglobal; i++) {
	    rowind[i] = i;
	    mtxX->entries[colind[i]] = buffer[i];
	}
	phgFree(buffer);
    }
#endif

/*static int jjj = 0;
if (++jjj == 3 && map->rank == 0) {
DenseMtx_writeForHumanEye(mtxX, stdout);
printf("\n");
phgFinalize();
exit(0);
}*/

    if (destroy)
	destroy_objects(solver);

    /* copy the solution to VEC x (Note: rowind[] is always sorted) */
#if 0
    DenseMtx_rowIndices(mtxX, &nrow, &rowind);
    entries = DenseMtx_entries(mtxX);
    for (i = 0; i < nrow; i++) {
	/* TODO: need to redistribute the solution */
	assert(rowind[i] >= n0 && rowind[i] - n0 < Nlocal);
	x->data[rowind[i] - n0] = entries[i];
    }
#else
    {
	double *buffer;
	int *cnts, *dsps;
	IV *mapIV;
	/* gather the solution to processor 0 */
	mapIV = IV_new();
	IV_init(mapIV, map->nglobal, NULL);
	IV_fill(mapIV, 0);
	mtxT = DenseMtx_MPI_splitByRows(mtxX, mapIV, _t->stats, msglvl,
					stdout, firsttag, comm);
	DenseMtx_free(mtxX);
	IV_free(mapIV);
	mtxX = mtxT;
	DenseMtx_rowIndices(mtxX, &nrow, &rowind);
	assert((map->rank == 0 && nrow == map->nglobal) ||
	       (map->rank != 0 && nrow == 0));
	entries = DenseMtx_entries(mtxX);
	cnts = phgAlloc(2 * map->nprocs * sizeof(*cnts));
	dsps = cnts + map->nprocs;
#if FT_PHG == FT_DOUBLE
	buffer = x->data;
#else	/* FT_PHG == FT_DOUBLE */
	buffer = phgAlloc(Nlocal * sizeof(*buffer));
#endif	/* FT_PHG == FT_DOUBLE */
	for (i = 0; i < map->nprocs; i++) {
	    cnts[i] = map->partition[i + 1] - map->partition[i];
	    dsps[i] = map->partition[i];
	}
	MPI_Scatterv(entries, cnts, dsps, MPI_DOUBLE, buffer, Nlocal,
		     MPI_DOUBLE, 0, comm);
#if FT_PHG != FT_DOUBLE
	for (i = 0; i < Nlocal; i++)
	    x->data[i] = buffer[i];
	phgFree(buffer);
#endif	/* FT_PHG != FT_DOUBLE */
	phgFree(cnts);
    }
#endif
    DenseMtx_free(mtxX);

    solver->residual = 0.0;	/* fake residual norm */
    return solver->nits = 1;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverSPOOLES_ = {
    "SPOOLES", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, FALSE, TRUE, FALSE, FALSE
};

#endif	/* USE_SPOOLES */
