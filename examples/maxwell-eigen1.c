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

/* $Id: maxwell-eigen1.c,v 1.19 2012/11/20 08:35:55 zlb Exp $
 *
 * Compute eigenvalues of time harmonic Maxwell's equations:
 *	\curl \curl E = \lambda E,	in \Omega,
 *	\div E = 0,			in \Omega,
 *	E\times n = 0,			on \partial\Omega
 * This code uses the following mixed (saddle point) formulation
 * (with $p\in H_0^1(\Omega)$):
 *	\curl \curl E + \grad p = \lambda E,	in \Omega,
 *	\div E = 0,				in \Omega,
 *	E\times n = 0,				on \partial\Omega
 * Cf., e.g.:
 * 	S. Adam, P. Arbenz, R. Geus,
 * 	Eigenvalue Solvers for Electromagnetic Fields in Cavities,
 * 	Tech. Rep. 275, Institute of Scientific Computing, ETH Z\"urich, 1997 */

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>

static void
build_matrices(MAT *matA, MAT *matM, DOF *u_h, DOF *p_h)
{
    int N = u_h->type->nbas;	/* number of basis functions in an element */
    int M = p_h->type->nbas;
    int i, j;
    GRID *g = u_h->g;
    ELEMENT *e;
    FLOAT A[N][N], B[N][N], C[M][N];
    INT I[N], Ip[N];
    INT k;

    assert(DofTypeDim(u_h) == Dim);

    assert(u_h->dim == 1);
    ForAllElements(g, e) {
	for (i = 0; i < N; i++)
	    I[i] = phgMapE2L(matA->rmap, 0, e, i);
	for (i = 0; i < M; i++)
	    Ip[i] = phgMapE2L(matA->rmap, 1, e, i);
	for (i = 0; i < N; i++) {
	    for (k = 0; k < M; k++) {
		/* \int \grad\psi_k\cdot\phi_i */
		C[k][i] =
		    phgQuadGradBasDotBas(e, p_h, k, u_h, i, QUAD_DEFAULT);
	    }
	    for (j = 0; j <= i; j++) {
		/* \int \phi_j\cdot\phi_i */
		B[j][i] = B[i][j] = 
		    phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
		/* \int \curl\phi_j\cdot\curl\phi_i */
		A[j][i] = A[i][j] =
		    phgQuadCurlBasDotCurlBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	    }
	}

	/* loop on basis functions of u_h */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, NULL, NULL, NULL, DOF_PROJ_CROSS))
		continue;
	    phgMatAddEntries(matA, 1, I + i, N, I, A[i]); 
	    phgMatAddEntries(matM, 1, I + i, N, I, B[i]); 
	    for (j = 0; j < M; j++)
		phgMatAddEntry(matA, I[i], Ip[j], C[j][i]); 
	}

	/* loop on basis functions of p_h */
	for (i = 0; i < M; i++)
	    phgMatAddEntries(matA, 1, Ip + i, N, I, C[i]); 
    }
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

    INT pre_refines = 0, mem_max = 400;;
    size_t mem, mem_peak;
    int i, n, nit;
    GRID *g;
    ELEMENT *e;
    DOF **u_h, **error, **p_h;
    DOF_TYPE *type;
    MAP *map, *m;
    MAT *A, *M;
    FLOAT PEoo, thres;
    FLOAT *evals;
    double wtime;
    FLOAT tol = 1e-10;

    phgOptionsPreset("-dof_type HC4 -mumps_symmetry sym");

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", &fn);
    phgOptionsRegisterInt("-nev", "Number of eigen pairs", &nev);
    phgOptionsRegisterInt("-pre_refines", "Pre_refines", &pre_refines);
    phgOptionsRegisterInt("-mem_max", "Maximum memory", &mem_max);
    phgOptionsRegisterFloat("-tol", "Convergence tolerance", &tol);

    phgInit(&argc, &argv);
    assert(DOF_DEFAULT->fe_space == FE_Hcurl);

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    /* pre-refinement */
    phgRefineAllElements(g, pre_refines);

    evals = phgCalloc(nev, sizeof(*evals));
    u_h = phgAlloc(nev * sizeof(*u_h));
    p_h = phgAlloc(nev * sizeof(*p_h));
    error = phgAlloc(nev * sizeof(*error));
    i = DOF_DEFAULT->order;
    assert(i + 1 <= 15);
    type = DOF_Pn[i + (DOF_DEFAULT->np_edge == 1 ? 0 : 1)];
    for (i = 0; i < nev; i++) {
	u_h[i] = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
	phgDofSetDirichletBoundaryMask(u_h[i], BDRY_MASK);
	p_h[i] = phgDofNew(g, type, 1, "p_h", DofInterpolation);
	phgDofSetDirichletBoundaryMask(p_h[i], BDRY_MASK);
	/* set random initial values for the eigenvectors */
	phgDofRandomize(u_h[i], i == 0 ? 123 : 0);
	phgDofRandomize(p_h[i], 0);
	error[i] = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);
    }

    while (TRUE) {
	phgPrintf("\n");
	if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",(double)g->lif);
	map = phgMapCreate(u_h[0], p_h[0], NULL);
	m = phgMapRemoveBoundaryEntries(map);
	phgPrintf("*** %"dFMT" elems, %"dFMT" submesh%s, "
		  "DOF(u_h) + DOF(p_h) = %"dFMT"\n", g->nleaf_global,
			g->nprocs, g->nprocs > 1 ? "es" : "", m->nglobal);
	A = phgMapCreateMat(m, m);
	M = phgMapCreateMat(m, m);
	phgMapDestroy(&m);
	build_matrices(A, M, u_h[0], p_h[0]);
	wtime = phgGetTime(NULL);
	n = phgDofEigenSolve(A, M, nev, EIGEN_CLOSEST, 0.0, &nit, evals,
				map, u_h, p_h, NULL);
	phgMatDestroy(&M);
	phgMatDestroy(&A);
	phgMapDestroy(&map);
	phgPrintf("%d its, converged eigenvalues: %d, wtime: %0.2lfs\n",
			nit, n, phgGetTime(NULL) - wtime);
	for (i = 0; i < n; i++) {
	    phgPrintf("  %4d: tau = %0.6e, |u_h| = %0.6e, |p_h| = %0.6e\n",
			i, (double)evals[i], (double)phgDofNormL2(u_h[i]),
			   (double)phgDofNormL2(p_h[i]));
	}
	for (; i < nev; i++) {
	    phgDofRandomize(u_h[i], 0);
	    phgDofRandomize(p_h[i], 0);
	}

	PEoo = -1.0;
	for (i = 0; i < n; i++) {
	    thres = estimate_error(evals[i], u_h[i], error[i]);
	    if (PEoo < thres)
		PEoo = thres;
	}
	PEoo = Fabs(PEoo);
	phgPrintf("indicator=%0.3le\n", (double)PEoo);
	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		(double)mem / (1024.0 * 1024.0),
		(double)mem_peak / (1024.0 * 1024.0));
	if (PEoo < tol || mem_peak >= 1024 * (size_t)mem_max * 1024)
	    break;
	if (n == 0) {
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

    phgFree(evals);
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h + i);
	phgDofFree(p_h + i);
	phgDofFree(error + i);
    }
    phgFree(u_h);
    phgFree(p_h);
    phgFree(error);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
