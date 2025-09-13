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

/* $Id: maxwell-eigen.c,v 1.58 2012/11/20 08:35:55 zlb Exp $
 *
 * Compute eigenvalues of time harmonic Maxwell's equations:
 *	\curl \curl u = \lambda u
 * 
 * This code computes eigenvalues of the following modified system:
 * 	(A - C*S*C^t) U = \lambda U
 * where (phi_i \in H(\curl) and \psi_i \in H^1):
 * 	A_{i,j} = <\curl\phi_j, \curl\phi_i>,
 * 	C_{i,j} = <\grad\psi_j, \phi_i>,
 * 	S is a SPD matrix.
 * The above modified system has both negative and positive eigenvalues, with
 * the positive ones being the eigenvalues of the original system.
 * See, e.g., the following paper:
 * 	J. Coyle and P.D. Ledger,
 *	Evidence of Exponential Convergence in the Computation of
 *	Maxwell Eigenvalues,
 *	Computer Methods in Applied Mechanics and Engineering,
 *	2004; 194: 587-604.
 *
 * Here we set S = s*diag(M)^(-1) instead of s*I as suggested in the paper,
 * where s > 0 is a parameter and M is the mass matrix of the H^1 bases.
 * It makes the magnitude of spurious modes independent of the mesh size.
 */

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>

static void
build_matrices(MAT *matA, MAT *matM, MAT *matC, MAT *S, FLOAT s,
		DOF *u_h, DOF *p_h)
/* S is used to store s*diag(M(p_h))^(-1) */
{
    int N = u_h->type->nbas;	/* number of basis functions in an element */
    int M = p_h->type->nbas;
    int i, j;
    GRID *g = u_h->g;
    ELEMENT *e;
    FLOAT A[N][N], B[N][N], C[N][M];
    INT I[N], Ip[M];
    INT k, n0;
    VEC *V = phgMapCreateVec(S->rmap, 1);

    phgVecDisassemble(V);	/* for phgVecAddEntry */
    ForAllElements(g, e) {
	for (i = 0; i < N; i++) {
	    I[i] = phgMapE2L(matA->rmap, 0, e, i);
	    for (k = 0; k < M; k++) {
		/* \int \grad\psi_k\cdot\phi_i */
		C[i][k] = phgQuadGradBasDotBas(e, p_h, k, u_h, i, QUAD_DEFAULT);
	    }
	    for (j = 0; j <= i; j++) {
		/* \int \phi_i\cdot\phi_j */
		B[j][i] = B[i][j] = 
		    phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
		/* \int \curl\phi_i\cdot\curl\phi_j */
		A[j][i] = A[i][j] =
		    phgQuadCurlBasDotCurlBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	    }
	}

	for (i = 0; i < M; i++) {
	    Ip[i] = phgMapE2L(matC->cmap, 0, e, i);
	    if (Ip[i] < 0)	/* boundary entry */
		continue;
	    phgVecAddEntry(V, 0, Ip[i],
			phgQuadBasDotBas(e, p_h, i, p_h, i, QUAD_DEFAULT));
	}

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, NULL, NULL, NULL, DOF_PROJ_CROSS))
		continue;
	    phgMatAddEntries(matA, 1, I + i, N, I, A[i]); 
	    phgMatAddEntries(matM, 1, I + i, N, I, B[i]); 
	    phgMatAddEntries(matC, 1, I + i, M, Ip, C[i]); 
	}
    }

    phgVecAssemble(V);
    n0 = V->map->partition[V->map->rank];
    for (k = 0; k < V->map->nlocal; k++)
	phgMatAddGlobalEntry(S, k + n0, k + n0, s / V->data[k]);
    phgVecDestroy(&V);
    phgMatAssemble(S);
    phgMatSetupDiagonal(S);
    phgMatAssemble(matA);
    phgMatAssemble(matM);
    phgMatAssemble(matC);
}

static FLOAT
estimate_error(FLOAT lambda, DOF *u_h, DOF *error)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *curl_u_h, *jump1, *jump2, *residual, *tmp;

    curl_u_h = phgDofCurl(u_h, NULL, NULL, NULL);
    residual = phgDofGetSameOrderDG(u_h, -1, NULL);
    phgDofCopy(u_h, &residual, NULL, "residual");
    jump2 = phgQuadFaceJump(residual, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    tmp = phgDofCurl(curl_u_h, NULL, NULL, NULL);
    phgDofAXPBY(-1.0, tmp, lambda, &residual);
    phgDofFree(&tmp);
    jump1 = phgQuadFaceJump(curl_u_h, DOF_PROJ_CROSS, NULL, QUAD_DEFAULT);
    ForAllElements(g, e) {
	int i;
	FLOAT eta, h;
	FLOAT diam = phgGeomGetDiameter(g, e);
	e->mark = 0;		/* clear refinement mmark */
	eta = 0.0;
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] == DIRICHLET || e->bound_type[i] == NEUMANN)
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    eta += (*DofFaceData(jump1, e->faces[i]) +
		    *DofFaceData(jump2, e->faces[i])) * h;
	}
	eta = eta + diam * diam * phgQuadDofDotDof(e, residual, residual,
						   QUAD_DEFAULT);
	*DofElementData(error, e->index) = eta;
    }
    phgDofFree(&jump1);
    phgDofFree(&jump2);
    phgDofFree(&curl_u_h);
    phgDofFree(&residual);
    return Sqrt(phgDofNormInftyVec(error));
}

/* callback function which computes (B = A - s*C*C^t)*x */
static int
funcB(MAT_OP op, MAT *B, VEC *x, VEC *y)
/* computes y = op(A) * x */
{
    MAT *A = ((void **)B->mv_data)[0];
    MAT *C = ((void **)B->mv_data)[1];
    MAT *S = ((void **)B->mv_data)[2];
    VEC *z;
    INT i;

    if (op == MAT_OP_D) {
	phgMatSetupDiagonal(A);
	for (i = 0; i < x->map->nlocal; i++) {
	    const MAT_ROW *row = phgMatGetRow(C, i);
	    INT j;
	    FLOAT d;
	    d = 0.0;
	    for (j = 0; j < row->ncols; j++)
		d += row->data[j] * S->diag[j] * row->data[j];
	    y->data[i] = A->diag[i] - d * x->data[i];
	}
	return 0;
    }

    phgMatVec(op, 1.0, A, x, 0.0, &y);			/* y := op(A)*x */
    z = phgMatVec(MAT_OP_T, 1.0, C, x, 0.0, NULL);	/* z := C^t*x */
    for (i = 0; i < z->map->nlocal; i++)		/* z := S*z */
	z->data[i] *= S->diag[i];
    phgMatVec(MAT_OP_N, -1.0, C, z, 1.0, &y);		/* y := y - C*z */
    phgVecDestroy(&z);

    return 0;
}

int
main(int argc, char *argv[])
{
#if 1
    /* \lambda_lnm = \pi^2(l^2+m^2+n^2) with lm+ln+nm > 0 */
    static char *fn = "../test/cube.dat";
    INT nev = 20;
#else
    /* j_m(\sqrt{\lambda}) == 0 (TE) or j'_m(\sqrt{\lambda'}) == 0 (TM),
     * where j_m is the m-th order spherical Bessel function, j'_m its
     * derivative */
    static char *fn = "../test/sphere.dat";
    INT nev = 20;
#endif

    FLOAT s = 10.0;
    BOOLEAN show_spurious = FALSE;
    BOOLEAN matrix_free = FALSE;
    INT pre_refines = 0, mem_max = 400;
    size_t mem, mem_peak;
    int i, j, n, nit;
    GRID *g;
    ELEMENT *e;
    DOF **u_h, **error, *p_h, *div;
    MAP *map, *map1, *m;
    MAT *A, *M, *S, *C, *B;
    FLOAT tol = 1e-10;
    FLOAT PEoo, thres;
    FLOAT *evals;
    double t0, t1, t2;
    long long nnz0, nnz1;

    phgOptionsPreset("-dof_type HC4");

#if USE_MUMPS
    phgOptionsPreset("-mumps_symmetry sym");
#endif	/* USE_MUMPS */
#if USE_PARDISO
    phgOptionsPreset("-pardiso_symmetry sym");
#endif	/* USE_PARDISO */

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", &fn);
    phgOptionsRegisterInt("-nev", "Number of eigen pairs", &nev);
    phgOptionsRegisterInt("-pre_refines", "Pre_refines", &pre_refines);
    phgOptionsRegisterInt("-mem_max", "Maximum memory", &mem_max);
    phgOptionsRegisterFloat("-tol", "Tolerance", &tol);
    phgOptionsRegisterFloat("-s", "Positive parameter", &s);
    phgOptionsRegisterNoArg("-show_spurious", "Show spurious modes",
				&show_spurious);
    phgOptionsRegisterNoArg("-matrix_free", "Use matrix-free matrix",
				&matrix_free);

    phgInit(&argc, &argv);

    assert(DOF_DEFAULT->fe_space == FE_Hcurl);
#if USE_ARPACK
    if (matrix_free)
	phgOptionsSetKeyword("-arpack_mode", "regular-invert");
#endif	/* USE_ARPACK */

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    /* pre-refinement */
    phgRefineAllElements(g, pre_refines);

    evals = phgCalloc(nev, sizeof(*evals));
    u_h = phgAlloc(nev * sizeof(*u_h));
    error = phgAlloc(nev * sizeof(*error));
    for (i = 0; i < nev; i++) {
	u_h[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
	phgDofSetDirichletBoundaryMask(u_h[i], BDRY_MASK);
	/* set random initial values for the eigenvectors */
	phgDofRandomize(u_h[i], i == 0 ? 123 : 0);
	error[i] = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);
    }

    i = DOF_DEFAULT->order;
    assert(i + 1 <= 15);
    if (DOF_DEFAULT->np_edge == 1) {
	/* First kind Nedelec element */
	p_h = phgDofNew(g, DOF_Pn[i], 1, "p_h", DofNoAction);
    }
    else {
	/* Second kind Nedelec element */
	p_h = phgDofNew(g, DOF_Pn[i + 1], 1, "p_h", DofNoAction);
    }

    while (TRUE) {
	phgPrintf("\n");
	if (phgBalanceGrid(g, 1.2, -1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",(double)g->lif);
	m = phgMapCreate(p_h, NULL);
	map1 = phgMapRemoveBoundaryEntries(m);
	phgMapDestroy(&m);
	map = phgMapCreate(u_h[0], NULL);
	m = phgMapRemoveBoundaryEntries(map);
	phgPrintf("*** %"dFMT" elems, %"dFMT" submesh%s, "
		  "DOF(u_h) = %"dFMT", DOF(p_h) = %"dFMT"\n",
			g->nleaf_global, g->nprocs, g->nprocs > 1 ? "es" : "",
			m->nglobal, map1->nglobal);
	A = phgMapCreateMat(m, m);
	M = phgMapCreateMat(m, m);
	C = phgMapCreateMat(m, map1);
	S = phgMapCreateMat(map1, map1);
	phgMapDestroy(&map1);
	phgMapDestroy(&m);
	build_matrices(A, M, C, S, s, u_h[0], p_h);
	if (matrix_free) {
	    /* matrix-free: B = A - C*S*C^t (passing A, C and S in mv_data) */
	    B = phgMapCreateMatrixFreeMat(A->rmap, A->cmap, funcB,
					  /* mv_data[] = */ A, C, S, NULL);
	}
	else {	/* (matrix_free) */
	    nnz0 = A->nnz_d + A->nnz_o;
	    /* B := C*S */
	    t0 = phgGetTime(NULL);
	    B = phgMatMat(MAT_OP_N, MAT_OP_N, 1.0, C, S, 0.0, NULL);
	    /* A = A - B*C^t */
	    t1 = phgGetTime(NULL);
	    phgMatMat(MAT_OP_N, MAT_OP_T, -1.0, B, C, 1.0, &A);
	    t2 = phgGetTime(NULL);
	    nnz1 = A->nnz_d + A->nnz_o;
#if USE_MPI
	    if (g->nprocs > 1) {
		long long d0[2] = {nnz0, nnz1}, d[2];
		MPI_Reduce(d0, d, 2, MPI_LONG_LONG, MPI_SUM, 0, g->comm);
		nnz0 = d[0];
		nnz1 = d[1];
	    }
#endif	/* USE_MPI */
	    phgPrintf("phgMatMat time: C*S^t=%g, C*C^t=%g\n", t1 - t0, t2 - t1);
	    phgPrintf("nnz(A) = %lld, nnz(B) = %lld\n", nnz0, nnz1);
	    phgMatDestroy(&B);
	    B = A;
	    A = NULL;
	}	/* (matrix_free) */
	t0 = phgGetTime(NULL);
	n = phgDofEigenSolve(B, M, nev, EIGEN_CLOSEST, 1.0, &nit, evals,
				map, u_h, NULL);
	phgMapDestroy(&map);
	phgMatDestroy(&B);
	phgMatDestroy(&C);
	phgMatDestroy(&M);
	phgMatDestroy(&A);
	phgMatDestroy(&S);
	for (i = n; i < nev; i++)
	    phgDofRandomize(u_h[i], 0);
	phgPrintf("%d its, converged eigenvalues: %d, wtime: %0.2lfs\n",
			nit, n, phgGetTime(NULL) - t0);
	if (show_spurious) {
	    /* count number of spurious modes */
	    for (i = 0; i < n; i++)
		if (evals[i] > 0.0)
		    break;
	    for (j = 0; j < i; j++)
		phgPrintf("  %4d: tau =%0.6e, |u_h| = %0.6e, spurious mode\n",
		    (j - i), (double)evals[j], (double)phgDofNormL2(u_h[j]));
	}
	PEoo = 0.0;
	for (i = j = 0; i < n; i++) {
	    if (evals[i] <= 0.0)
		continue;
	    div = phgDofDivergence(u_h[i], NULL, NULL, NULL);
	    phgPrintf("  %4d: tau = %0.6e, |u_h| = %0.6e, |div(u_h)| = %0.6e\n",
			j, (double)evals[i], (double)phgDofNormL2(u_h[i]),
			(double)phgDofNormL2(div));
	    phgDofFree(&div);
	    thres = estimate_error(evals[i], u_h[i], error[i]);
	    if (PEoo < thres)
		PEoo = thres;
	    j++;
	}
	if (j == 0)
	    PEoo = 1e10;
	phgPrintf("indicator=%0.3le\n", (double)PEoo);
	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		(double)mem / (1024.0 * 1024.0),
		(double)mem_peak / (1024.0 * 1024.0));
	if (PEoo <= tol || mem_peak >= 1024 * (size_t)1024 * mem_max)
	    break;
	if (j == 0) {
	    phgRefineAllElements(g, 1);
	    continue;
	}
	ForAllElements(g, e)
	    e->mark = 0;
	/*thres = Pow(1e-2,2) / (double)g->nleaf_global;*/
	thres = Pow(PEoo * 0.5, 2);
	for (i = 0; i < n; i++)
	    ForAllElements(g, e)
		if (*DofElementData(error[i],e->index) >= thres)
		    e->mark = 1;
	phgRefineMarkedElements(g);
    }

#if 1
    phgPrintf("Creating \"%s\".\n", phgExportVTKn(g, "output.vtk", n, u_h));
#endif

    phgDofFree(&p_h);
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h + i);
	phgDofFree(error + i);
    }

    phgFree(evals);
    phgFree(u_h);
    phgFree(error);

    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
