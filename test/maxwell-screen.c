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

/* $Id: maxwell-screen.c,v 1.61 2020/12/06 09:03:08 zlb Exp $
 *
 * see maxwell-screen.tex for a description of the problem */

#define TEST_PRECOND	0	/* code for testing preconditioners */

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* for mkdir() */
#include <sys/stat.h>
#include <sys/types.h>

static FLOAT KK = 1.0;

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    *(values++) = 1.;
    *(values++) = 1.;
    *(values++) = 1.;
    return;
}

static void
func_g(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    *(values++) = 0.;
    *(values++) = 0.;
    *(values++) = 0.;
}

static void
func_mu(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    /* *value = ((x >= -0.5 && x <= 0.5 && y >= -0.5 && y <= 0.5 &&
       z >= -0.5 && z <= 0.5) ? 1.0 : 2.0); */
    *value = 1.0;
}

static double
elapsed_time(GRID *g, BOOLEAN flag, double mflops)
/* returns elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;

    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    mem = phgMemoryUsage(g, NULL);
    if (flag) {
	if (mflops <= 0)
	    phgPrintf("[%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
	else
	    phgPrintf("[%0.4lgMB %0.4lfs %0.4lgGF]\n",
		      mem / (1024.0 * 1024.0), et, mflops * 1e-3);
    }

    return et;
}

static FLOAT KK_pc = 0.0;	/* will be set to -KK after phgInit() */

#if TEST_PRECOND
static MAT *matA;
static SOLVER *m_solver;
#endif	/* TEST_PRECOND */

static void
build_linear_system(SOLVER *solver, DOF *mu, DOF *u_h, DOF *f_h, FLOAT kk)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    int i, j;
    int N = u_h->type->nbas * u_h->dim;	/* size of local matrix */
    FLOAT mass, A[N][N], buffer[N], rhs[N];
#if TEST_PRECOND
    FLOAT M[N][N];
#endif	/* TEST_PRECOND */
    INT I[N];

    ForAllElements(g, e) {
	/* compute \int \curl\phi_j \cdot \curl\phi_i making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (j = 0; j <= i; j++) {
		mass = kk * phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
		A[j][i] = A[i][j] =
		    phgQuadCurlBasACurlBas(e, u_h, j, mu, u_h, i,
					   QUAD_DEFAULT) - mass;
#if TEST_PRECOND
		M[j][i] = M[i][j] = mass;
#endif	/* TEST_PRECOND */
	    }
	}

	/* compute local matrix and RHS */
	for (i = 0; i < N; i++) {	/* loop on basis functions */
	    if (phgDofDirichletBC(u_h, e, i, kk == KK_pc ? NULL : func_g,
				  buffer, rhs + i, DOF_PROJ_CROSS)) {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer);
#if TEST_PRECOND
		phgSolverAddMatrixEntries(m_solver, 1, I + i, N, I, buffer);
#endif	/* TEST_PRECOND */
	    }
	    else {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, &A[i][0]);
		/* right hand side: \int f * phi_i */
		rhs[i] = phgQuadDofDotBas(e, f_h, u_h, i, QUAD_DEFAULT);
#if TEST_PRECOND
		phgSolverAddMatrixEntries(m_solver, 1, I + i, N, I, &A[i][0]);
#endif	/* TEST_PRECOND */
	    }
	}
	phgSolverAddRHSEntries(solver, N, I, rhs);
    }
}

static FLOAT
estimate_error(DOF *mu, DOF *u_h, DOF *f_h, DOF *error)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *jump1, *jump2, *res, *div, *curl_u_h, *tmp;

    tmp = phgDofCurl(u_h, NULL, NULL, NULL);
    curl_u_h =
	phgDofMM(MAT_OP_N, MAT_OP_N, 1, 3, 1, 1.0, mu, 1, tmp, 0.0, NULL);
    phgDofFree(&tmp);

    res = phgDofGetSameOrderDG(u_h, -1, NULL);
    phgDofCopy(f_h, &res, NULL, "residual");
    phgDofAXPY(KK, u_h, &res);	/* res = f + k^2 u_h */
    jump1 = phgQuadFaceJump(curl_u_h, DOF_PROJ_CROSS, NULL, QUAD_DEFAULT);
    jump2 = phgQuadFaceJump(res, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    tmp = phgDofCopy(res, NULL, NULL, NULL);	/* save f_h + k^2 u_h */
    div = phgDofCurl(curl_u_h, NULL, NULL, "tmp");
    phgDofAXPBY(-1.0, div, 1.0, &res);
    phgDofFree(&div);
    div = phgDofDivergence(tmp, NULL, NULL, "div(f+k^2 u_h)");
    phgDofFree(&tmp);
    ForAllElements(g, e) {
	int i;
	FLOAT eta, h;
	FLOAT diam = phgGeomGetDiameter(g, e);
	eta = 0.0;
	for (i = 0; i < NFace; i++) {
	    INT fno;
	    if (e->bound_type[i] & BDRY_MASK)
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    fno = e->faces[i];
	    eta += (*DofFaceData(jump1, fno) + *DofFaceData(jump2, fno)) * h;
	}
	eta =
	    eta + diam * diam * (phgQuadDofDotDof(e, res, res, QUAD_DEFAULT) +
				 phgQuadDofDotDof(e, div, div, QUAD_DEFAULT));
	*DofElementData(error, e->index) = eta;
    }
    phgDofFree(&jump1);
    phgDofFree(&jump2);
    phgDofFree(&res);
    phgDofFree(&div);
    phgDofFree(&curl_u_h);
    return Sqrt(phgDofNormInftyVec(error));
}

/*==================================================================*/

static void
pc_proc(void *ctx, VEC *b, VEC **x)
{
    SOLVER *pc_solver = ctx;
    FLOAT *old_rhs;
 
    phgSolverAssemble(pc_solver);      /* assemble matrix and RHS */

    if (pc_solver->mat->rmap->nlocal != b->map->nlocal) {
	phgPrintf("Invalid PC!\n");
	return;
    }

    /* updating RHS */  
    old_rhs = pc_solver->rhs->data;
    pc_solver->rhs->data = b->data;
    pc_solver->rhs->assembled=TRUE;
    pc_solver->rhs_updated=TRUE;
    bzero((*x)->data, (*x)->map->nlocal * sizeof(*(*x)->data));

    /* solving M x = b */
    phgSolverVecSolve(pc_solver, FALSE, *x);
    pc_solver->rhs->data = old_rhs;

#if TEST_PRECOND
#endif	/* TEST_PRECOND */

    return;  
}

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *mu, *u_h, *f_h, *error;
    SOLVER *solver, *pc, *pc_pc;
    FLOAT PE_l2 = 0.0, PE_infinity = 0.0, Norm_u_h;
    size_t mem, mem_peak;

    char *pc_solver = NULL;
    char *pc_opts = NULL;
    char *solver_opts = NULL;
    char *fn = "screen3.dat";
    char *vtk_fn = "maxwell-screen.vtk";
    char *matrix_path = "mat";
    INT pre_refines = 2;

    INT nmax = 5000000;
    INT nmax_global = 900000000;
    INT vtk_limit = 0;
    FLOAT tol = 1e-8;		/* convergence tolerance */
    BOOLEAN export_mesh = FALSE;
    BOOLEAN dump_matrix = FALSE;
    INT refine_depth = 1;	/* <=0 => uniformly refine -depth or 3 times */
    INT mem_max = 400;		/* max. memory (MB) / process */
    FLOAT lif_threshold = 1.2;
    INT submesh_threshold = 1000;
    INT pc_maxit = 1;		/* -1 means exact solution */
    BOOLEAN enable_pc = TRUE;
    BOOLEAN use_ams = TRUE;
    int refine_counts = 1;	/* current refinement level */
    double t0, t1, tt[3];
    BOOLEAN stop;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterFloat("-lif_threshold",
		"LIF threshold for mesh redistribution", &lif_threshold);
    phgOptionsRegisterInt("-submesh_threshold",
		"Min # of elements in a submesh", &submesh_threshold);
    phgOptionsRegisterNoArg("-dump_matrix", "Save matrix and RHS for MATLAB",
		&dump_matrix);
    phgOptionsRegisterFilename("-matrix_path",
		"Directory for storing matrix and RHS files",
		(char **)&matrix_path);
    phgOptionsRegisterFloat("-kk", "k^2", &KK);
    phgOptionsRegisterFloat("-kk_pc", "k^2 for preconditioning system", &KK_pc);
    phgOptionsRegisterNoArg("-export_mesh", "Export mesh in VTK format",
		&export_mesh);
    phgOptionsRegisterFilename("-vtk_file", "VTK filename", (char **)&vtk_fn);
    phgOptionsRegisterInt("-vtk_limit", "Max # of element exported to VTK file",
		&vtk_limit);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-mem_max", "Max-memory (MB) / process", &mem_max);
    phgOptionsRegisterInt("-max_mem", "Max-memory (MB) / process", &mem_max);
    phgOptionsRegisterInt("-nmax", "Maxi. number of elements/process", &nmax);
    phgOptionsRegisterInt("-nmax_global", "Maxi. number of elements",
				&nmax_global);
    phgOptionsRegisterFloat("-tol", "Convergence tolerance", &tol);
    phgOptionsRegisterInt("-pc_maxit",
		"Maxi. it number in pc (-1 ==> exact PC)", &pc_maxit);
    phgOptionsRegisterInt("-refine_depth",
		"Refinement depth (<=0 ==> uniform refinement)",
		&refine_depth);
    phgOptionsRegisterString("-solver_opts",
		"Options for the main solver", (char **)&solver_opts);
    phgOptionsRegisterNoArg("-enable_pc", "Enable SPD preconditioner",
		&enable_pc);
    phgOptionsRegisterNoArg("-use_ams", "Use AMS preconditioner", &use_ams);
    phgOptionsRegisterString("-pc_solver",
		"Solver for the preconditioning system", (char **)&pc_solver);
    phgOptionsRegisterString("-pc_opts",
		"Options for the preconditioner", (char **)&pc_opts);

    phgOptionsPreset("-refine_strategy gers");
    phgOptionsPreset("-dof_type ND1");

    /* The main solver */
#if USE_MINRES
    phgOptionsPreset("-solver minres");
#else	/* USE_MINRES */
    phgOptionsPreset("-solver pcg");
#endif	/* USE_MINRES */
    phgInit(&argc, &argv);

    /* Set default solver for the preconditioning system */
    if (pc_solver == NULL)
	phgOptionsSetOptions("-pc_solver pcg");

    phgOptionsShowUsed();

    if (KK_pc == 0.0)
	KK_pc = -KK;

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    phgRefineAllElements(g, pre_refines);

    phgGetTime(tt);
    t1 = t0 = tt[2];

    u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
    phgDofSetDataByValue(u_h, 0.0);

    f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h", func_f);
    error = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);

    mu = phgDofNew(g, DOF_P0, 1, "mu", DofInterpolation);
    phgDofSetDataByFunction(mu, func_mu);

    while (TRUE) {
	elapsed_time(g, FALSE, 0);
	phgPrintf("\n%d DOF, %d elements, %d submesh%s, LIF: %lg\n",
		  DofDataCountGlobal(u_h), g->nleaf_global,
		  g->nprocs, g->nprocs > 1 ? "es" : "",
		  (double)g->lif);
	if (phgBalanceGrid(g, lif_threshold, submesh_threshold, NULL, 0)) {
	    phgPrintf("  Repartition mesh, %d submeshes, LIF =  %lg ",
		      g->nprocs, (double)g->lif);
	    elapsed_time(g, TRUE, 0);
#if USE_MPI
	    MPI_Bcast(&refine_counts, 1, MPI_INT, 0, g->comm);
#endif	/* USE_MPI */
	}
	phgPrintf("  Set up linear solver: ");
	phgOptionsPush();
	phgOptionsSetInt("-verbosity", 0);
	phgOptionsSetOptions(solver_opts);
	solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
#if TEST_PRECOND
	matA = solver->mat;
	{
	    MAT *M = phgMapCreateMat(matA->rmap, matA->cmap);
	    m_solver = phgMat2Solver(SOLVER_PCG, M);
	    phgMatDestroy(&M);
	}
#endif	/* TEST_PRECOND */
	phgOptionsPop();
	phgPrintf("solver LIF = %lg ", (double)solver->mat->rmap->lif);
	elapsed_time(g, TRUE, 0);

	phgPrintf("  Build linear system: ");
	phgPerfGetMflops(g, NULL, NULL);	/* reset Mflops count */
	build_linear_system(solver, mu, u_h, f_h, KK);
	elapsed_time(g, TRUE, 0);

	if (dump_matrix) {
	    char s[1024];
	    if (g->rank == 0 && refine_counts == 1)
		mkdir(matrix_path, 0755);
	    sprintf(s, "%s/matrix_%d.m", matrix_path, refine_counts);
	    phgPrintf("  Dumping matrix to file \"%s\".\n", s);
	    phgMatDumpMATLAB(solver->mat, "A", s);
	    sprintf(s, "%s/RHS_%d.m", matrix_path, refine_counts);
	    phgPrintf("  Dumping RHS to file \"%s\".\n", s);
	    phgVecDumpMATLAB(solver->rhs, "b", s);
	}

	if (enable_pc && solver->oem_solver->allow_extern_pc) {
	    phgPrintf("  Build preconditioner: ");
	    phgOptionsPush();
	    phgOptionsSetKeyword("-solver", pc_solver);
	    phgOptionsSetInt("-verbosity", 0);
	    phgOptionsSetOptions(pc_opts);
	    pc = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
	    build_linear_system(pc, mu, u_h, f_h, KK_pc);

	    /* build the aux preconditioner for the preconditioning system */
	    if (use_ams) {
		pc_pc = phgMat2Solver(SOLVER_AMS, pc->mat);
		phgSolverAMSSetConstantCoefficients(pc_pc, 1.0, -KK_pc);
		pc_pc->rtol = 0.;
		pc_pc->maxit = 1;
		phgSolverSetPC(pc, pc_pc, NULL);
	    }
	    else {
		pc_pc = NULL;
	    }
#if 0
	    if (pc->oem_solver == SOLVER_HYPRE && -KK_pc != 1.0)
		phgSolverHypreAMSSetConstantPoisson(pc, 1.0, -KK_pc);
#endif

	    if (pc_maxit > 0) {
		phgSolverSetMaxIt(pc, pc_maxit);
		pc->warn_maxit = FALSE;
	    }

	    phgSolverSetPC(solver, pc, pc_proc);
	    phgOptionsPop();
	    elapsed_time(g, TRUE, 0);
	}
	else {
	    pc = pc_pc = NULL;
	}

	phgPrintf("  Assemble solver and preconditioners: ");
	if (pc_pc != NULL)
	    phgSolverAssemble(pc_pc);
	if (pc != NULL)
	    phgSolverAssemble(pc);
	phgSolverAssemble(solver);
	elapsed_time(g, TRUE, 0);

	phgPrintf("  Solve linear system: ");
phgInfo(0, "=========================== %d DOF\n", DofDataCountGlobal(u_h));
	phgSolverSolve(solver, TRUE, u_h, NULL);
	phgPrintf("nits = %d, resid = %0.4lg ", solver->nits,
		  (double)solver->residual);
	elapsed_time(g, TRUE, 0);
	if (pc_pc != NULL)
	    phgSolverDestroy(&pc_pc);
#if TEST_PRECOND
	phgSolverDestroy(&m_solver);
#endif	/* TEST_PRECOND */
	if (pc != NULL)
	    phgSolverDestroy(&pc);
	phgSolverDestroy(&solver);

	phgPrintf("  Estimate: ");
	Norm_u_h = phgDofNormL2(u_h);
	phgPrintf("|u_h|=%0.3le, ", (double)Norm_u_h);

	PE_infinity = estimate_error(mu, u_h, f_h, error);
	PE_l2 = Sqrt(phgDofNormL1Vec(error));
	phgPrintf("PEoo=%0.3le, PE_L2=%0.3le ",
		  (double)PE_infinity, (double)PE_l2);
	elapsed_time(g, TRUE, 0);
	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("  Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		  mem / (1024.0 * 1024.0), mem_peak / (1024.0 * 1024.0));
	phgPrintf("  Peformance of the step: %0.4lgMflops\n",
			phgPerfGetMflops(g, NULL, NULL));
	phgGetTime(tt);
	phgPrintf("  Adaptive step: %d, elapsed time: %lg, step time: %lg\n",
		  refine_counts, (double)(tt[2] - t0), (double)(tt[2] - t1));
	t1 = tt[2];

	stop = (PE_l2 < tol || 
		(nmax > 0 && g->nleaf_global / phgNProcs > nmax) ||
		(nmax_global > 0 && g->nleaf_global > nmax_global) ||
		mem_peak >= 1024 * (size_t) 1024 * mem_max);

	if (export_mesh &&
	    ((vtk_limit > 0 && g->nleaf_global > vtk_limit) || stop)) {
	    export_mesh = FALSE;
	    phgPrintf("  ===> \"%s\", ",
		      phgExportVTK(g, vtk_fn, mu, u_h, error, NULL));
	    if (phgRank == 0) {
		FILE *f;
		long long size;
		f = fopen(vtk_fn, "r");
		fseek(f, 0L, SEEK_END);
		size = ftell(f);
		fclose(f);
		phgPrintf("size = %lld ", size);
	    }
	    elapsed_time(g, TRUE, 0.);
	}

	if (stop)
	    break;

	phgPrintf("  Refine mesh: ");
	elapsed_time(g, FALSE, 0);
	if (refine_depth <= 0) {
	    /*  uniform mesh refinement */
	    phgRefineAllElements(g, refine_depth == 0 ? 3 : -refine_depth);
	}
	else {
	    phgMarkRefine(MARK_DEFAULT, error, Pow(0.8, 2), NULL, 0.0,
			  refine_depth, Pow(tol, 2) / g->nleaf_global);
	    phgRefineMarkedElements(g);
	}
	elapsed_time(g, TRUE, 0);
	refine_counts++;
    }

    phgDofFree(&mu);
    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&error);

    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
