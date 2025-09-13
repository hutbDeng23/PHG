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

/* Two-grid solver/preconditioner.
 * $Id: solver-2grid.c,v 1.38 2022/09/21 02:35:28 zlb Exp $
 */

#include "phg.h"

#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <limits.h>	/* INT_MAX */

/* PENALIZE_BDRY controls BC treatment for the coarse matrix */
#define PENALIZE_BDRY 1		/* !0: penalize diagonal rows,
				 *  0: set RHS of bdry rows to 0. */

typedef struct {
    struct TEMP_DOF_LIST {
	DOF **list;
	int n, alloc;
    }		dofs;
    GRID	*g2;		/* coarse mesh */
    MAT		*P;		/* The transfer matrix */
    MAT		*C;		/* The coarse matrix */
    SOLVER	*fine_solver, *coarse_solver;
    char	*coarse_grid_file;
    char	*cycle_type;
    char	*fine_solver_opts, *coarse_solver_opts;
    int		coarse_grid_refine;
    int		coarse_matrix_type;
    int		coarse_grid_nprocs;
    BOOLEAN	assembled;
    BOOLEAN	g2_owned;
} OEM_DATA;

static char *coarse_grid_file = NULL, *cycle_type = "10";
static char *fine_solver_opts = NULL, *coarse_solver_opts = NULL;
static INT coarse_matrix_type = 0;
static INT coarse_grid_refine = 0;
static INT coarse_grid_nprocs = 0;	/* <=0 -> same as fine grid */

/* user function for building the coarse matrix */
#if 0	/* simple test BUILD_FUNC for the Laplace operator (test only) */
#pragma message("*** WARNING: test-only code.")
    static MAT *
    build_matrix(SOLVER *solver)
    {
	MAT *C;
	MAP *map2 = _t->P->cmap;
	DOF *u;
	ELEMENT *e;

	assert(map2->ndof == 1);
	u = map2->dofs[0];
	assert(u->dim == 1);
	C = phgMapCreateMat(map2, map2);
	ForAllElements(u->g, e) {
	    int i, j, N = DofNBas(u, e);
	    FLOAT A[N][N], buffer[N], *row = NULL;
	    INT I[N];

	    /* compute \int \grad\phi_j \cdot \grad\phi_i using symmetry */
	    for (i = 0; i < N; i++) {
		I[i] = phgMapE2L(C->rmap, 0, e, i);
		for (j = 0; j <= i; j++)
		    A[j][i] = A[i][j] =
			phgQuadGradBasDotGradBas(e, u, j, u, i, QUAD_DEFAULT);
	    }

	    /* loop on basis functions */
	    for (i = 0; i < N; i++) {
		if (phgDofDirichletBC(u,e,i, NULL,buffer,NULL, DOF_PROJ_NONE))
		    row = buffer;
		else if (phgDofGetElementBoundaryType(u, e, i) & NEUMANN)
		    phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__);
		else
		    row = A[i];
		phgMatAddEntries(C, 1, I + i, N, I, row); 
	    }
	}
	phgMatAssemble(C);
	return C;
    }
    static BUILD_FUNC coarse_matrix_func = build_matrix;
#else
    static BUILD_FUNC coarse_matrix_func = NULL;
#endif
static GRID *coarse_grid = NULL;

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

/*---------------------------- exported functions -------------------------*/

BUILD_FUNC
phgSolver2GridSetCoarseMatFunc(BUILD_FUNC f)
{
    BUILD_FUNC old_func = coarse_matrix_func;
    coarse_matrix_func = f;

    return old_func;
}

void
phgSolver2GridSetCoarseGrid(GRID *g)
{
    coarse_grid = g;
    return;
}

/*-------------------------------------------------------------------------*/

static DOF *
get_dof(SOLVER *solver, GRID *g, DOF_TYPE *type, int dim, const char *name,
	const char *file, int line)
/* This function manages the list of temporary DOFs attached to the SOLVER.
 * It returns a stored or dynamically created DOF of required type and dim. */
{
    static char *prefix = "temp_";
    char *temp_name;
    int i;
    DOF *u;

    for (i = 0; i < _t->dofs.n; i++)
	if ((u = _t->dofs.list[i])->g == g && u->type == type && u->dim == dim
		/*&& !strcmp(u->name + strlen(prefix), name)*/)
	    return u;		/* return the DOF found */

    /* create a new DOF and append it to the list */
    temp_name = phgAlloc(strlen(prefix) + strlen(name) + 1);
    sprintf(temp_name, "%s%s", prefix, name);
    u = phgDofNew_(g, type, NULL, dim, temp_name, DofNoData, file, line);
    phgFree(temp_name);
    if (_t->dofs.n >= _t->dofs.alloc)
	_t->dofs.list = phgRealloc_(_t->dofs.list,
			    (_t->dofs.alloc += 4) * sizeof(*_t->dofs.list),
			    _t->dofs.n * sizeof(*_t->dofs.list));
    _t->dofs.list[_t->dofs.n++] = u;

    return u;
}

/* ---------------------------------------------------
 * Build interpolation matrix for one single DOF_TYPE
 * ---------------------------------------------------*/

typedef struct PT_INFO_ {
    int rank;			/* rank */
    INT id;			/* vector index */
    FLOAT x[Dim];		/* xyz */
} PT_INFO;

static void
func_xyz(FLOAT x, FLOAT y, FLOAT z, FLOAT *v)
{
    *(v++) = x;
    *(v++) = y;
    *(v++) = z;
}

#define CHECK_INTERP_MATRIX 0
/* Note CHECK_INTERP_MATRIX != 0 => test only: print |Pu_H - u_h| and exit.
 * 	for ((i=0; i<16; i+=3)); do
 * 	    ./poisson -pre_refines $i -mesh_file ../test/cube.dat \
 * 	    	-solver 2grid -2grid_cycle_type "1" \
 * 	    	-2grid_coarse_grid_file ../albert-files/cube5.dat \
 * 	    	-2grid_coarse_grid_refine $i -dof_type P1 \
 * 	    	| grep "Mat projection error:" */
#if CHECK_INTERP_MATRIX
static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *v)
{
    *v = Exp(x + 2. * y - 3. * z);
}
#endif	/* CHECK_INTERP_MATRIX */

static MAT *
build_P0(SOLVER *solver, DOF *u, OCT_TREE *og)
/* returns the transfer matrix for a single DOF_TYPE.
 *
 * Author: LENG Wei (wleng@lsec.cc.ac.cn)
 * */
{
    GRID *g = u->g;
    GRID *g2 = _t->g2;
    SIMPLEX *e2;
    DOF *u1, *u2, *coord1;;
    MAP *map1, *map2, *map1_coord, *map2_ext;
    VEC *vec1, *vec2, *vec1_coord, *vec2_ext;
    MAT *matP;
    VEC *vecP;
#if USE_MPI
    MPI_Datatype type;
#endif	/* USE_MPI */
    FLOAT *xyz, lambda[Dim + 1];
    INT i, nsend, nrecv;
    int ii, rank, nprocs2 = g2->nprocs;
    int *scnts, *sdsps, *rcnts, *rdsps;
    BOOLEAN *flag;
    FLOAT (*bboxs)[Dim][2];
    PT_INFO *sbuf_pt, *rbuf_pt, *pt_info;

    /* this code only works with non-DG Lagrangian type bases */
    assert(u->type != NULL && u->type->is_nodal);
    assert(u->type->fe_space != FE_L2);	/* to be removed */

    if (u->type->fe_space == FE_L2) {
	/* check availability of remote vertex neighbours */
	if (g->halo == NULL || GetHaloKind(g->halo->type) != HALO_VERT)
	    phgError(1, "%s:%d: vertex halo is required for DG type.\n",
			__FILE__, __LINE__);
    }

    /* bcast bboxs of coarse grid */
    bboxs = phgCalloc(nprocs2, sizeof(*bboxs));
    if (phgRank == 0)
	memcpy(&bboxs[0][0][0], &og->bboxs[0][0][0],
	       nprocs2 * sizeof(*bboxs));
#if USE_MPI
    MPI_Bcast(&bboxs[0][0][0], nprocs2 * Dim * 2, PHG_MPI_FLOAT, 0, g->comm);
#endif	/* USE_MPI */

    if (phgVerbosity >= 1) {
	for (rank = 0; rank < nprocs2; rank++) {
	    phgInfo(1,
		"rank %d: [%12.4f,%12.4f]x[%12.4f,%12.4f]x[%12.4f,%12.4f]\n",
		rank, bboxs[rank][0][0], bboxs[rank][0][1], bboxs[rank][1][0],
		bboxs[rank][1][1], bboxs[rank][2][0], bboxs[rank][2][1]
	    );
	}
    }

    /* DOF1 */
    if (u->dim == 1)
	u1 = u;
    else
	u1 = get_dof(solver, g, u->type, 1, u->name, __FILE__, __LINE__);
    map1 = phgMapCreate(u1, NULL);
    vec1 = phgMapCreateVec(map1, 1);

    /* DOF2 */
    if (g2 != NULL) {
	u2 = get_dof(solver, g2, u->type, 1, u->name, __FILE__, __LINE__);
	map2 = phgMapCreate(u2, NULL);
	vec2 = phgMapCreateVec(map2, 1);
    }
#if 0
    /* extend map2 */
    if (nprocs2 < g->nprocs) {
	/* TODO: remove {map2,vec2}_ext */
	INT nglobal;
#if USE_MPI
	if (g->rank == 0)
	    nglobal = map2->nglobal;
	MPI_Bcast(&nglobal, 1, PHG_MPI_INT, 0, g->comm);
#else
	nglobal = map2->nglobal;
#endif	/* USE_MPI */
	map2_ext = phgMapCreateSimpleMap(g->comm,
				(g2 != NULL) ? map2->nlocal : 0, nglobal);
	vec2_ext = phgMapCreateVec(map2_ext, 1);
    }
    else
#endif
    {
	map2_ext = map2;
	vec2_ext = vec2;
    }

    /* Construct the coarse to fine transfer matrix P.
     *
     * The transfer matrix is constricted by interpolating/projecting the
     * coarse grid basis functions as linear combinations of the fine grid
     * basis functions.
     *
     * Here, both set of basis functions are assumed of Lagrangian type,
     * so we just need to evaluate values of the coarse grid basis functions
     * at the nodes of the fine grid basis functions. */

    matP = phgMapCreateMat(map1, map2_ext);
    vecP = phgMapCreateVec(map1, 1);
    phgVecDisassemble(vecP);
    phgInfo(1, "map1: %d\n", map1->nlocal);

    /* coord */
    coord1 = phgDofNew(g, u1->type, Dim, "coord1", func_xyz);
    map1_coord = phgMapCreate(coord1, NULL);
    vec1_coord = phgMapCreateVec(map1_coord, 1);
    phgMapDofToLocalData(map1_coord, 1, &coord1, vec1_coord->data);

    /* search */
    flag = phgCalloc(map1->nlocal, sizeof(*flag));

    xyz = vec1_coord->data;
    for (i = 0; i < map1->nlocal; i++) {	/* vector index */
	BOOLEAN found;
	INT iG = map1->partition[map1->rank] + i;	/* v2G */
	phgInfo(3, "dof coord: %12.8e, %12.8e, %12.8e\n",
		xyz[0], xyz[1], xyz[2]);

	/* local */
	if (og)
	    found = phgOctTreeLocatePoint(og, xyz, &e2, lambda);
	else
	    found = FALSE;
	if (found) {
	    int N = u2->type->nbas;
	    INT IG[N];		/* global index */
	    FLOAT *bas, one = 1.;

	    for (ii = 0; ii < N; ii++)
		IG[ii] = phgMapE2G(map2, 0, e2, ii);

	    bas = (FLOAT *)u2->type->BasFuncs(u2, e2, 0, -1, lambda);
	    phgMatAddGlobalEntries(matP, 1, &iG, N, IG, bas);
	    phgVecAddGlobalEntry(vecP, 0, iG, one);

	    if (phgVerbosity >= 3) {
		phgInfo(3, "g1: %d\n", iG);
		for (ii = 0; ii < N; ii++)
		    phgInfo(3, "   %d %12.8e\n", IG[ii], bas[ii]);
	    }

	    flag[i] = TRUE;	/* added */
	}

	xyz += Dim;
    }

    /*
     *
     * send search info
     *
     * */
    scnts = phgCalloc(4 * g->nprocs, sizeof(*scnts));
    sdsps = scnts + g->nprocs;
    rcnts = sdsps + g->nprocs;
    rdsps = rcnts + g->nprocs;

    /* count send */
    xyz = vec1_coord->data;
    for (i = 0; i < map1->nlocal; i++) {	/* vector index */
	if (!flag[i]) {

	    for (rank = 0; rank < nprocs2; rank++) {
		if (rank != g->rank
		    && bboxs[rank][0][0] - OCT_TREE_EPS < xyz[0]
		    && xyz[0] < bboxs[rank][0][1] + OCT_TREE_EPS
		    && bboxs[rank][1][0] - OCT_TREE_EPS < xyz[1]
		    && xyz[1] < bboxs[rank][1][1] + OCT_TREE_EPS
		    && bboxs[rank][2][0] - OCT_TREE_EPS < xyz[2]
		    && xyz[2] < bboxs[rank][2][1] + OCT_TREE_EPS) {
		    assert(scnts[rank] < INT_MAX);
		    scnts[rank]++;
		}
	    }
	}
	xyz += Dim;
    }

#if USE_MPI
    MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, g->comm);
#else	/* USE_MPI */
    rcnts[0] = scnts[0];
#endif	/* USE_MPI */
    nsend = nrecv = 0;
    for (rank = 0; rank < g->nprocs; rank++) {
	assert(nsend <= INT_MAX && nrecv <= INT_MAX);
	sdsps[rank] = nsend;
	rdsps[rank] = nrecv;
	nsend += scnts[rank];
	nrecv += rcnts[rank];
    }

    sbuf_pt = phgCalloc((size_t)nsend + (size_t)nrecv, sizeof(*sbuf_pt));
    rbuf_pt = sbuf_pt + nsend;

    /* fill send buf */
    xyz = vec1_coord->data;
    for (i = 0; i < map1->nlocal; i++) {	/* vector index */
	if (!flag[i]) {
	    INT iG = map1->partition[map1->rank] + i;	/* v2G */

	    for (rank = 0; rank < nprocs2; rank++) {
		if (rank != g->rank
		    && bboxs[rank][0][0] - OCT_TREE_EPS < xyz[0]
		    && xyz[0] < bboxs[rank][0][1] + OCT_TREE_EPS
		    && bboxs[rank][1][0] - OCT_TREE_EPS < xyz[1]
		    && xyz[1] < bboxs[rank][1][1] + OCT_TREE_EPS
		    && bboxs[rank][2][0] - OCT_TREE_EPS < xyz[2]
		    && xyz[2] < bboxs[rank][2][1] + OCT_TREE_EPS) {
		    assert(sdsps[rank] < INT_MAX);
		    pt_info = sbuf_pt + sdsps[rank]++;
		    pt_info->rank = g->rank;
		    pt_info->id = iG;
		    memcpy(pt_info->x, xyz, Dim * sizeof(FLOAT));
		}
	    }
	}
	xyz += Dim;
    }

    /* restore face sdsps[] */
    for (rank = 0; rank < g->nprocs; rank++)
	sdsps[rank] -= scnts[rank];

#if USE_MPI
    /* alltoallv */
    MPI_Type_contiguous(sizeof(*sbuf_pt), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(sbuf_pt, scnts, sdsps, type,
		  rbuf_pt, rcnts, rdsps, type, g->comm);
    MPI_Type_free(&type);
#else	/* USE_MPI */
    memcpy(rbuf_pt + rdsps[0], sbuf_pt + sdsps[0], scnts[0] * sizeof(*sbuf_pt));
#endif	/* USE_MPI */

    /*
     *
     * step 3. remote local search
     *
     * */
    if (og == NULL)
	assert(nrecv == 0);

    pt_info = rbuf_pt;
    for (i = 0; i < nrecv; i++, pt_info++) {
	BOOLEAN found;
	INT iG = pt_info->id;

	phgInfo(3, "to find %12.8e %d %d\n",
		pt_info->x[0], pt_info->rank, pt_info->id);
	found = phgOctTreeLocatePoint(og, pt_info->x, &e2, lambda);

	if (found) {
	    int N = u2->type->nbas;
	    INT IG[N];		/* global index */
	    FLOAT *bas;

	    for (ii = 0; ii < N; ii++)
		IG[ii] = phgMapE2G(map2, 0, e2, ii);

	    bas = (FLOAT *)u2->type->BasFuncs(u2, e2, 0, -1, lambda);
	    phgMatAddGlobalEntries(matP, 1, &iG,	/* add to remote row */
				   N, IG, bas);
	    /* phgVecAddGlobalEntry(vecP, 0, iG, one); */
	    /* phgInfo(0, "vecP: %x %d\n", vecP->O2Gmap); */

	    if (phgVerbosity >= 3) {
		phgInfo(3, "g1: %d\n", iG);
		for (ii = 0; ii < N; ii++)
		    phgInfo(3, "   %d %12.8e\n", IG[ii], bas[ii]);
	    }

	    pt_info->id = -1;	/* marked as found */
	}
    }

    /*
     * Send back found info.
     *  */
    /* alltoallv */
#if USE_MPI
    MPI_Type_contiguous(sizeof(*sbuf_pt), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(rbuf_pt, rcnts, rdsps, type,
		  sbuf_pt, scnts, sdsps, type, g->comm);
    MPI_Type_free(&type);
#else	/* USE_MPI */
    memcpy(sbuf_pt + sdsps[0], rbuf_pt + rdsps[0], scnts[0] * sizeof(*sbuf_pt));
#endif	/* USE_MPI */

    xyz = vec1_coord->data;
    for (i = 0; i < map1->nlocal; i++) {	/* vector index */
	if (!flag[i]) {
	    INT iG = map1->partition[map1->rank] + i;	/* v2G */

	    int found = 0;
	    for (rank = 0; rank < nprocs2; rank++) {
		if (rank != g->rank
		    && bboxs[rank][0][0] - OCT_TREE_EPS < xyz[0]
		    && xyz[0] < bboxs[rank][0][1] + OCT_TREE_EPS
		    && bboxs[rank][1][0] - OCT_TREE_EPS < xyz[1]
		    && xyz[1] < bboxs[rank][1][1] + OCT_TREE_EPS
		    && bboxs[rank][2][0] - OCT_TREE_EPS < xyz[2]
		    && xyz[2] < bboxs[rank][2][1] + OCT_TREE_EPS) {
		    assert(sdsps[rank] < INT_MAX);
		    pt_info = sbuf_pt + sdsps[rank]++;
		    if (pt_info->id == -1)
			found++;
		}
	    }
	    if (found == 0) {
		phgError(1, "can't find point (%12.8e, %12.8e, %12.8e)\n",
			 pt_info->x[0], pt_info->x[1], pt_info->x[2]
		    );
	    }

	    phgVecAddGlobalEntry(vecP, 0, iG, (FLOAT)found);
	}
	xyz += Dim;
    }

    /* restore face sdsps[] */
    for (rank = 0; rank < g->nprocs; rank++)
	sdsps[rank] -= scnts[rank];

    /* Assemble */
    phgMatAssemble(matP);	/* Costly !!! */
    phgVecAssemble(vecP);

    /* rescaling (?????) */
    if (TRUE) {
	MAT_ROW *row = matP->rows;
	const FLOAT *v = vecP->data;
	FLOAT sum;

	assert(matP->type == PHG_UNPACKED);
	for (i = 0; i < matP->rmap->nlocal; i++) {
	    INT iG = map1->partition[map1->rank] + i;
	    assert(*v > .5);
	    phgInfo(3, "row: %d, scale: %e\n", iG, *v);

	    sum = 0.;
	    for (ii = 0; ii < row->ncols; ii++) {
		phgInfo(3, "   col: %d, %e\n", row->cols[ii], row->data[ii]);
		row->data[ii] /= *v;
		sum += row->data[ii];
	    }
	    if (Fabs(sum - 1.) > 1e-12)
		phgError(1, "   row[%d] sum %e not one !\n", iG, sum);

	    row++;
	    v++;
	}
    }

    phgFree(flag);
    phgFree(scnts);
    phgFree(sbuf_pt);

#if CHECK_INTERP_MATRIX
    /* Function test */
    phgDofSetDataByFunction(u1, func_u);
    phgMapDofToLocalData(map1, 1, &u1, vec1->data);
    if (g2 != NULL) {
	phgDofSetDataByFunction(u2, func_u);
	phgMapDofToLocalData(map2, 1, &u2, vec2->data);
	memcpy(vec2_ext->data, vec2->data,
	       map2->nlocal * sizeof(*vec2->data));
    }

    //phgVecDumpMATLAB(vec1, "v1", "v1_.m");
    //phgVecDumpMATLAB(vec2, "v2", "v2_.m");
    phgMatVec(MAT_OP_N, -1., matP, vec2_ext, 1., &vec1);
    //phgVecDumpMATLAB(vec1, "v1", "r_.m");
    phgPrintf("\n\n\n---\nMat projection error: %e\n",
	      phgVecNorm2(vec1, 0, NULL));

MPI_Finalize();
exit(0);
#endif	/* CHECK_INTERP_MATRIX */

    phgVecDestroy(&vecP);

    phgVecDestroy(&vec1_coord);
    phgMapDestroy(&map1_coord);
    phgDofFree(&coord1);

    phgVecDestroy(&vec1);
    phgMapDestroy(&map1);

    if (vec2_ext != vec2)
	phgVecDestroy(&vec2_ext);
    if (map2_ext != map2)
	phgMapDestroy(&map2_ext);
    if (g2 != NULL) {
	phgVecDestroy(&vec2);
	phgMapDestroy(&map2);
    }

    phgFree(bboxs);

    return matP;
}

static void
build_P(SOLVER *solver, OCT_TREE *og)
/*----------------------------------------------------------------------------
 * Builds the transfer matrix _t->P(map_f, map_c).
 * 'solver' is the top-level solver (corresponding to the fine grid)
 *
 * Denote by {\phi} and {\psi}, respectively, the finite element bases on
 * the fine and coarse grid, then the elements of the transfer matrix P
 * is determined by representing \phi with \psi:
 *		\psi_j = \sum_k P_{k,j} \phi_k
 *--------------------------------------------------------------------------*/
{
    MAT *P = NULL;
    DOF *u = solver->mat->rmap->dofs[0];
    BOOLEAN free_og = FALSE;

    assert(solver->oem_solver == SOLVER_2Grid);

    if (og == NULL) {
	/* import coarse grid */
	assert(_t->g2 == NULL);
	if (coarse_grid != NULL) {
	    /* use user-provided coarse grid */
	    _t->g2 = coarse_grid;
	    _t->g2_owned = FALSE;
	}
	else {
	    /* import the coarse grid */
	    int nprocs = _t->coarse_grid_nprocs;
	    int refine = _t->coarse_grid_refine, r;
	    MPI_Comm comm;
	    GRID *g = u->g;

	    if (_t->coarse_grid_file == NULL)
		phgError(1, "%s (%s:%d): coarse grid file unspecified.\n",
			__func__, __FILE__, __LINE__);
	    _t->g2 = phgNewGrid(-1);
	    if (nprocs <= 0 || nprocs > g->nprocs)
		nprocs = g->nprocs;
	    if (solver->monitor)
		phgPrintf("* 2-grid: reading coarse mesh file \"%s\"\n.",
				_t->coarse_grid_file);
	    comm = g->comm;
#if USE_MPI
	    if (nprocs < g->nprocs) {
		MPI_Comm_split(g->comm, phgRank < nprocs ? 1 : MPI_UNDEFINED,
			       0, &comm);
		assert(g->rank >= nprocs || comm != MPI_COMM_NULL);
	    }
#endif	/* USE_MPI */
	    if (comm != MPI_COMM_NULL &&
		!phgImport_(_t->g2, _t->coarse_grid_file, comm))
		phgError(1, "%s (%s:%d): can't read grid file \"%s\".\n",
			__func__, __FILE__, __LINE__, _t->coarse_grid_file);
	    while (comm != MPI_COMM_NULL && refine > 0) {
		r = refine > 3 ? 3 : refine;
		phgRefineAllElements(_t->g2, r);
		phgRedistributeGrid_(_t->g2, comm, NULL);
		refine -= r;
		if (solver->monitor)
		    phgPrintf("* 2-grid: refining coarse mesh, "
				"%d procs, %"dFMT" elements\n.",
				_t->g2->nprocs, _t->g2->nelem_global);
	    }
#if USE_MPI
	    if (comm != MPI_COMM_NULL && comm != g->comm)
		MPI_Comm_free(&comm);
#endif	/* USE_MPI */
	    _t->g2_owned = TRUE;
	}
	og = phgOctTreeBuild(_t->g2);
	free_og = TRUE;
    }

    if (solver->mat->rmap->ndof == 0) {
	phgError(1, "%s (%s:%d): unimplemented case.\n",
			__func__, __FILE__, __LINE__);
	return;
    }

    if (solver->mat->rmap->ndof == 1) {
	if (u->dim > 1) {
	    /* TODO: expand each element of matrix P to a dim*dim block */
	    phgError(1, "%s (%s:%d): unimplemented case.\n",
			__func__, __FILE__, __LINE__);
	}
	P = build_P0(solver, u, og);
    }
    else {
	phgError(1, "%s (%s:%d): unimplemented case.\n",
				__func__, __FILE__, __LINE__);
    }

    if (free_og)
	phgOctTreeDestroy(&og);

    _t->P = P;
}

static void
build_C(SOLVER *solver)
/* builds the coarse grid matrix.
 *	_t->coarse_matrix_type == 0: use L_H
 *	_t->coarse_matrix_type == 1: use PtAP
 *	_t->coarse_matrix_type == 2: use PtAP with fill-in of L_H
 */
{
    MAT *T;

    assert(solver->oem_solver == SOLVER_2Grid && _t->P != NULL);

    switch (_t->coarse_matrix_type) {
	case 0:		/* C = L_H */
	    if (coarse_matrix_func == NULL)
		phgError(1, "%s (%s:%d): user function for building coarse "
			"grid matrix is not set!\n",
			__func__, __FILE__, __LINE__);
	    _t->C = (*coarse_matrix_func)(solver);
	    break;
	case 1:		/* C = PtAP */
	    if (_t->P->rmap->nglobal < _t->P->cmap->nglobal)
		phgError(1, "%s:%d: more DOF on coarse grid (%"dFMT") than on "
			 "fine grid (%"dFMT")!\n", __FILE__, __LINE__,
			 _t->P->cmap->nglobal, _t->P->rmap->nglobal);
#if 1
	    /* (PtA)P */
	    T = phgMatMat(MAT_OP_T, MAT_OP_N, 1., _t->P, solver->mat, 0., NULL);
	    _t->C = phgMatMat(MAT_OP_N, MAT_OP_N, 1., T, _t->P, 0., NULL);
#else
	    /* Pt(AP) */
	    T = phgMatMat(MAT_OP_N, MAT_OP_N, 1., solver->mat, _t->P, 0., NULL);
	    _t->C = phgMatMat(MAT_OP_T, MAT_OP_N, 1., _t->P, T, 0., NULL);
#endif
	    phgMatDestroy(&T);
	    break;
	case 2:		/* C = PtAP with fill-in of L_H */
	    phgError(1, "%s (%s:%d): unimplemented coarse matrix type (%d).\n", 
			__func__, __FILE__, __LINE__, _t->coarse_matrix_type);
	    break;
	default:
	    phgError(1, "%s:%d: invalid coarse matrix type (%d).\n",
			__FILE__, __LINE__, _t->coarse_matrix_type);
    }
}

/*----------------------------------------------------------------------*/

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nThe 2-grid solver options:", "\n", "2grid");
    phgOptionsRegisterString("-2grid_cycle_type",
		"Cycle type ('0': fine grid solution, "
			"'1': coarse grid correction)", &cycle_type);
    phgOptionsRegisterString("-2grid_fine_solver_opts",
		"Fine-grid solver options", &fine_solver_opts);
    phgOptionsRegisterString("-2grid_coarse_solver_opts",
		"Coarse-grid solver options", &coarse_solver_opts);
    phgOptionsRegisterInt("-2grid_coarse_matrix_type",
		"Coarse matrix type (0: L_H, 1: PtAP, "
			"2: PtAP with fill-in of L_H)", &coarse_matrix_type);
    phgOptionsRegisterFilename("-2grid_coarse_grid_file",
		"Coarse grid filename", &coarse_grid_file);
    phgOptionsRegisterInt("-2grid_coarse_grid_refine",
		"Coarse grid refinement depth", &coarse_grid_refine);
    phgOptionsRegisterInt("-2grid_coarse_grid_nprocs",
		"Number of processes for the coarse grid", &coarse_grid_nprocs);

    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));

    /* copy relavant cmdline arguments to OEM_DATA */

    _t->cycle_type = strdup(cycle_type);
    _t->fine_solver_opts =
		fine_solver_opts == NULL ? NULL : strdup(fine_solver_opts);
    _t->coarse_solver_opts =
		coarse_solver_opts == NULL ? NULL : strdup(coarse_solver_opts);
    _t->coarse_grid_file =
		coarse_grid_file == NULL ? NULL : strdup(coarse_grid_file);
    _t->coarse_grid_refine = coarse_grid_refine;
    _t->coarse_matrix_type = coarse_matrix_type;
    _t->coarse_grid_nprocs = coarse_grid_nprocs;
    _t->assembled = FALSE;

    return 0;
}

static int
Create(SOLVER *solver)
{
    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }
    return 0;
}

static int
Destroy(SOLVER *solver)
{
    int i;

    if (solver->oem_data == NULL)
	return 0;

    if (_t->fine_solver != NULL)
	phgSolverDestroy(&_t->fine_solver);
    if (_t->coarse_solver != NULL)
	phgSolverDestroy(&_t->coarse_solver);
    if (_t->P != NULL)
	phgMatDestroy(&_t->P);
    if (_t->C != NULL)
	phgMatDestroy(&_t->C);

    for (i = 0; i < _t->dofs.n; i++)
	phgDofFree(_t->dofs.list + i);
    phgFree(_t->dofs.list);

    if (_t->g2 != NULL && _t->g2_owned)
	phgFreeGrid(&_t->g2);

    phgFree(_t->cycle_type);
    phgFree(_t->fine_solver_opts);
    phgFree(_t->coarse_solver_opts);
    phgFree(_t->coarse_grid_file);

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    double t = 0.0, t1;

    if (_t->assembled)
	return 0;
    _t->assembled = TRUE;

    if (solver->monitor)
	t = phgGetTime(NULL);

    if (strchr(_t->cycle_type, '1') != NULL) {
	/* setup coarse grid correction */
	build_P(solver, NULL);
	if (solver->monitor) {
	    phgPrintf("build_P time: %lg\n", (t1 = phgGetTime(NULL)) - t);
	    t = t1;
	}
	build_C(solver);
	if (solver->monitor) {
	    phgPrintf("build coarse matrix time: %lg\n",
				(t1 = phgGetTime(NULL)) - t);
	    t = t1;
	}
    }

    /* build fine solver */
    if (strchr(_t->cycle_type, '0') != NULL) {
	phgOptionsPush();
	phgSolverSetDefaultSuboptions();
	phgOptionsSetOptions("-solver asm");
	phgOptionsSetOptions("-solver_rtol 0. -solver_btol 0. -solver_atol 0. "
			     "-solver_maxit 1");
	phgOptionsSetOptions(_t->fine_solver_opts);
	_t->fine_solver = phgMat2Solver(SOLVER_DEFAULT, solver->mat);
	_t->fine_solver->warn_maxit = FALSE;
	phgVecDestroy(&_t->fine_solver->rhs);
	_t->fine_solver->rhs_updated = TRUE;
	phgOptionsPop();
    }

    /* build coarse solver */
    if (_t->C != NULL) {
	phgOptionsPush();
	phgSolverSetDefaultSuboptions();
	phgOptionsSetOptions("-solver asm");
	phgOptionsSetOptions("-solver_rtol 0. -solver_btol 0. -solver_atol 0. "
			     "-solver_maxit 1");
	phgOptionsSetOptions(_t->coarse_solver_opts);
	_t->coarse_solver = phgMat2Solver(SOLVER_DEFAULT, _t->C);
	_t->coarse_solver->warn_maxit = FALSE;
	phgVecDestroy(&_t->coarse_solver->rhs);
	_t->coarse_solver->rhs_updated = TRUE;
	phgOptionsPop();
    }

#if PENALIZE_BDRY
    if (_t->C != NULL) {
	INT i;
	MAT_ROW *row;

	/* penalize boundary rows */
	for (i = 0; i < _t->C->rmap->nlocal; i++) {
	    row = _t->C->rows + i;
#if 1
	    assert(row->ncols != 0);
#else
	    if (row->ncols == 0) {
		if (row->alloc == 0) {
		    row->alloc = 1;
		    row->cols = phgAlloc(row->alloc * sizeof(*row->cols));
		    row->data = phgAlloc(row->alloc * sizeof(*row->data));
		}
		row->ncols = 1;
		row->cols[0] = i;
		row->data[0] = 0.0;
	    }
#endif
	    if (row->ncols == 1 && row->cols[0] == i)
		row->data[0] = DBL_MAX;
	}
    }
#endif	/* PENALIZE_BDRY */

    return 0;
}

static void
mult_prec(SOLVER *solver, VEC *r, VEC *x, const char *cycle)
{
    VEC *r0 = NULL, *r1 = NULL, *x1 = NULL;

    while (*cycle != '\0') {
	if (*cycle < '0' && *cycle > '1')
	    phgError(1, "invalid 2-grid cycle type string: %s\n",
							_t->cycle_type);
	switch (*(cycle++) - '0') {
	    case 0:	/* solving on fine mesh */
		if (phgVerbosity > 1)
		    phgPrintf("\t2-grid fine grid solution\n");
		_t->fine_solver->rhs = r;
		_t->fine_solver->rhs->assembled = TRUE;
		phgSolverVecSolve(_t->fine_solver, FALSE, x);
		_t->fine_solver->rhs = NULL;
		break;
	    case 1:	/* space P correction */
		if (phgVerbosity > 1)
		    phgPrintf("\t2-Grid coarse grid correction\n");
		phgVecCopy(r, &r0);
		phgMatVec(MAT_OP_N, -1.0, solver->mat, x, 1.0, &r0);
		if (x1 == NULL || x1->map != _t->P->cmap) {
		    phgVecDestroy(&x1);
		    x1 = phgMapCreateVec(_t->P->cmap, 1);
		}
		else {
		    bzero(x1->data, x1->map->nlocal * sizeof(*x1->data));
		}
		if (r1 != NULL && r1->map != _t->P->cmap)
		    phgVecDestroy(&r1);
		phgMatVec(MAT_OP_T, 1.0, _t->P, r0, 0.0, &r1);

#if !PENALIZE_BDRY
		if (_t->C != NULL) {
		    /* For boundary DOF: set RHS of bdry rows to 0 */
		    INT i;
		    MAP *map = /*_t->coarse_solver->mat*/C->rmap;
		    for (i = 0; i < map->bdry_nlocal; i++)
		    	r1->data[map->bdry[i]] = 0.;
		}		
#endif	/* PENALIZE_BDRY */

		_t->coarse_solver->rhs = r1;
		_t->coarse_solver->rhs->assembled = TRUE;
		phgSolverVecSolve(_t->coarse_solver, FALSE, x1);
		_t->coarse_solver->rhs = NULL;
		phgMatVec(MAT_OP_N, 1.0, _t->P, x1, 1.0, &x);
		break;
	}	/* case */
    }	/* while */

    if (r0 != NULL)
	phgVecDestroy(&r0);
    if (r1 != NULL)
	phgVecDestroy(&r1);
    if (x1 != NULL)
	phgVecDestroy(&x1);
}

static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    char *p = _t->cycle_type, *q;
    VEC *y0 = NULL, *x0 = NULL;
    FLOAT res = FLOAT_MAX, ib_norm = 0.0, ires0 = 0.0, tol = 0.0;

    Assemble(solver);

    if (solver->maxit > 0 || solver->btol > 0. || solver->rtol > 0. ||
	solver->atol > 0. || solver->monitor) {
	ib_norm = phgVecNorm2(solver->rhs, 0, NULL);
	if (ib_norm == 0.0) {
	    solver->nits = 0;
	    solver->residual = 0.;
	    bzero(x->data, x->nvec * x->map->nlocal * sizeof(*x->data));
	    return 0;
	}
	tol = solver->btol * ib_norm;
	if (tol < solver->atol)
	    tol = solver->atol;
	ib_norm = 1.0 / ib_norm;	/* inversed rhs norm */
    }

    solver->nits = 0;
    if (solver->monitor) {
        phgPrintf("======================== 2-Grid =====================\n");
	phgPrintf("Cycle type: %s\n", _t->cycle_type);
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
	    (double)solver->atol, (double)solver->rtol, (double)solver->btol);
        phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs\n");
        phgPrintf("-----   ------------   -------------  ------------\n");
    }
    while (solver->nits < solver->maxit) {
	solver->nits++;

	if ((q = strchr(p, '+')) == NULL) {
	    mult_prec(solver, solver->rhs, x, p);
	}
	else {
	    phgVecCopy(solver->rhs, &y0);
	    phgMatVec(MAT_OP_N, -1.0, solver->mat, x, 1.0, &y0);
	    x0 = phgMapCreateVec(y0->map, y0->nvec);
	    while (*p != '\0') {
		if (q != NULL)
		    *q = '\0';
		bzero(x0->data, sizeof(*x0->data) * x0->map->nlocal);
		mult_prec(solver, y0, x0, p);
		phgVecAXPBY(1.0, x0, 1.0, &x);
		if (q == NULL)
		    break;
		*q = '+';
		q = strchr(p = q + 1, '+');
	    }
	}

	if (solver->rtol > 0.0 || tol > 0.0 || solver->monitor) {
	    phgVecCopy(solver->rhs, &x0);
	    phgMatVec(MAT_OP_N, -1.0, solver->mat, x, 1.0, &x0);
	    res = phgVecNorm2(x0, 0, NULL);
	    if (ires0 == 0.0) {
		ires0 = solver->rtol * res;
		if (tol < ires0)
		    tol = ires0;
		ires0 = (res == 0. ? 1.0 : 1.0 / res);	/* inversed res norm */
	    }
	    if (solver->monitor)
		phgPrintf("% 5d   %12le   %12le   %12le\n", solver->nits,
				(double)res, (double)(res * ires0),
				(double)(res * ib_norm));
	    if (res <= tol)
		break;
	}
    }

    phgVecDestroy(&x0);
    phgVecDestroy(&y0);

    if (destroy) {
	Destroy(solver);
	phgMatDestroy(&solver->mat);
    }

    solver->residual = res/*_b*/;

    return solver->nits;
}

OEM_SOLVER phgSolver2Grid_ = {
    "2-Grid", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, TRUE, FALSE
};
