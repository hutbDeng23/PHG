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

/* $Id: poisson.c,v 1.339 2022/09/19 06:22:50 zlb Exp $ */

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <string.h>
#include <math.h>

/*
    This test code solves the Poisson equation with Dirichlet BC:

	- \Delta u = f
	u |_{\partial\Omega} = g
*/

/* The analytical solution, BC and RHS are defined in the file functions.h */
#ifndef TEST_CASE
# define TEST_CASE 3
#endif
#include "functions.h"

/* constants for evaluating L2 error indicator */
#define C0L2	1.0
#define C1L2	0.5
/* constants for evaluating H1 error indicator */
#define C0H1	1.0
#define C1H1	0.5

#define FUNC_DEF(func_f, f_)				\
static void						\
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)		\
{							\
    *value = f_(x, y, z);				\
    return;						\
}

FUNC_DEF(func_f, f_)	/* function for the right-hand-side */
FUNC_DEF(func_g, g_)	/* function for the boundary condition */
FUNC_DEF(func_u, u_)	/* function for the analytical solution */

static void
func_grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
/* returns the gradient of the analytical solution u */
{
    grad_u_(value, x, y, z);
}

static double
elapsed_time(GRID *g, BOOLEAN flag, double mflops)
/* returns and prints elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;
 
    /* get current wall time */
    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    /* get memory usage */
    mem = phgMemoryUsage(g, NULL);

    if (flag) {
	if (mflops <= 0)
	    phgPrintf("[%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
	else
	    phgPrintf("[%0.4lgMB %0.4lfs %0.4lgGF]\n", mem / (1024.0 * 1024.0),
			 et, mflops*1e-3);
    }

    return et;
}

static BOOLEAN qcache = FALSE;

static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    QCACHE *qc = NULL;
    int Q_f = 0;

    assert(u_h->dim == 1);

    if (qcache) {
	qc = phgQCNew(QD_DEFAULT, u_h);
	Q_f = phgQCAddFEFunction(qc, f_h);
    }

#if USE_OMP
#pragma omp parallel for private(e)
#endif  /* USE_OMP */
    ForAllElementsBegin(g, e) {
	int i, j;
	int N = DofNBas(u_h, e);    /* # basis functions in the element */
	FLOAT A[N][N], B[N], buffer[N], *row;
	INT I[N];

	if (qcache)
	    phgQCSetQuad(qc, phgQuadGetQuad3D(2 * DofTypeOrder(u_h, e)),
							phgGeomGetVolume(g, e));

	/* compute \int \grad\phi_j \cdot \grad\phi_i making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (j = 0; j <= i; j++) {
		FLOAT val;
		if (qcache)
		    val = phgQCIntegrate(qc, e->index, Q_GRAD, j,
					 qc, e->index, Q_GRAD, i);
		else
		    val = phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, -1);
		A[j][i] = A[i][j] = val;
	    }
	}

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, func_g,
				  buffer, B + i, DOF_PROJ_NONE)) {
		row = buffer;
	    }
	    else if (phgDofGetElementBoundaryType(u_h, e, i) & NEUMANN) {
		/* TODO */
		row = NULL;
		phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__);
	    }
	    else {
		/* right hand side: \int f * phi_i */
		if (qcache)
		    B[i] = phgQCIntegrate(qc, e->index, Q_BAS, i,
					  qc, e->index, Q_f, 0);
		else
		    phgQuadDofTimesBas(e, f_h, u_h, i, QUAD_DEFAULT, B + i);
		row = A[i];
	    }
	    phgSolverAddMatrixEntries(solver, 1, I + i, N, I, row); 
	}
	phgSolverAddRHSEntries(solver, N, I, B);
    } ForAllElementsEnd

    if (qcache)
	phgQCFree(&qc);
}

static void
estimate_error(DOF *u_h, DOF *f_h, DOF *grad_u, DOF *e_L2, DOF *e_H1)
/* compute L2 and H1 error indicators */
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *residual, *tmp;

    tmp = phgDofDivergence(grad_u, NULL, NULL, NULL);
    residual = phgDofGetSameOrderDG(u_h, -1, "residual");
    phgDofCopy(f_h, &residual, NULL, NULL);
    phgDofAXPY(1.0, tmp, &residual);
    phgDofFree(&tmp);

    if (qcache) {
	/* New code based on phgXFEMDot_ and phgDofFaceJump2_ */
	FLOAT *data, *data0, *data1;
	/* data := \int_e residual^2 */
	phgXFEMDot_(NULL, residual, NULL, residual, NULL, &data, 0);
	/* data0 := \sum_{F\in e(F)} h_F^3 * \int_F [grad_u.n]^2 */
	phgDofFaceJump2_(grad_u, NULL, PROJ_DOT, 3.0, INTERIOR|NEUMANN,
			 &data0, 0);
	/* data1 := \sum_{F\in e(F)} h_F^1 * \int_F [grad_u.n]^2 */
	phgDofFaceJump2_(grad_u, NULL, PROJ_DOT, 1.0, INTERIOR|NEUMANN,
			 &data1, 0);
#if USE_OMP
#pragma omp parallel for private(e)
#endif  /* USE_OMP */
	ForAllElementsBegin(g, e) {
	    FLOAT diam = phgGeomGetDiameter(g, e);
	    FLOAT d = phgGetBoundaryError(g, e);
	    e_L2->data[e->index] = Pow(diam, 4.0) * data[e->index] * C1L2
				    + data0[e->index] * C0L2 + d;
	    e_H1->data[e->index] = Pow(diam, 2.0) * data[e->index] * C1H1
				    + data1[e->index] * C0H1 + d;
	} ForAllElementsEnd
	phgFree(data);
	phgFree(data0);
	phgFree(data1);
    }
    else {
	/* Old code using phgQuadFaceJump */
	DOF *jump = phgQuadFaceJump(grad_u, DOF_PROJ_DOT, "jump", QUAD_DEFAULT);
    
#if USE_OMP
#pragma omp parallel for private(e)
#endif  /* USE_OMP */
	ForAllElementsBegin(g, e) {
	    int i;
	    FLOAT eta0, eta1, d, h;
	    FLOAT diam = phgGeomGetDiameter(g, e);
	    e->mark = 0;	    /* clear refinement mmark */
	    eta0 = eta1 = 0.0;
	    /* for each face F compute [grad_u \cdot n] */
	    for (i = 0; i < NFace; i++) {
		if (e->bound_type[i] & (DIRICHLET | NEUMANN))
		    continue;	    /* boundary face */
		h = phgGeomGetFaceDiameter(g, e, i);
		d = *DofFaceData(jump, e->faces[i]);
		eta0 += d * h * h * h;
		eta1 += d * h;
	    }
	    d = phgQuadDofDotDof(e, residual, residual, QUAD_DEFAULT);
	    eta0 = eta0 * C0L2 + diam * diam * diam * diam * d * C1L2;
	    eta1 = eta1 * C0H1 + diam * diam * d * C1H1;
    
	    /* add curved boundary errors (FIXME: how to normalize?) */
	    d = phgGetBoundaryError(g, e);
	    eta0 += d;
	    eta1 += d;
    
	    *DofElementData(e_L2, e->index) = eta0;
	    *DofElementData(e_H1, e->index) = eta1;
	} ForAllElementsEnd
	phgDofFree(&jump);
    }

    phgDofFree(&residual);

    return;
}

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *u_h, *f_h, *grad_u = NULL, *e_L2, *e_H1;
    SOLVER *solver;
    FLOAT L2_error = 0.0, H1_error = 0.0, EL2_error = 0.0, EH1_error = 0.0;
    int nits;
    BOOLEAN flag;
    size_t current_min, current_avg, current_max, peak_min, peak_avg, peak_max;

    /*------------------- Command-line options ----------------------*/

    /* List of valid export functions */
    const char *export_names[] = {"vtk", "opendx", "ensight", NULL};
    const char *(*export_list[])(GRID *, const char *, DOF *, ...) =
				{phgExportVTK, phgExportDX, phgExportEnsight};
    int export_index = 0;
    /* file extensions */
    const char *export_ext[] = {"vtk", "dx", "en"};
    FLOAT lif_threshold = 1.2;
    INT submesh_threshold = 1;

    /* Scalar parameters */
    INT		refine_depth = 1;
    FLOAT	tol = 4e-3;
    INT		mem_max = 400;
    INT		nmax = 0;
    BOOLEAN	export_mesh = FALSE;
    BOOLEAN	export_all = FALSE;
    BOOLEAN	dump_matrix = FALSE;
    BOOLEAN	dump_rhs = FALSE;
    BOOLEAN	L2adapt = FALSE;
    BOOLEAN	use_weights = FALSE;
    BOOLEAN	hp_test = FALSE;
    BOOLEAN	flatten = FALSE;
    BOOLEAN	init0 = FALSE;
    BOOLEAN	show_cond = FALSE;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-refine_depth",
			  "Minimum refine depth (<=0 ==> uniform refinement)",
			  &refine_depth);
    phgOptionsRegisterFloat("-tol", "Convergence criterion", &tol);
    phgOptionsRegisterFloat("-lif_threshold", "Load imbalance factor threshold",
			    &lif_threshold);
    phgOptionsRegisterInt("-submesh_threshold", "Submesh size threshold",
			  &submesh_threshold);
    phgOptionsRegisterInt("-mem_max", "Maximum memory (MB)", &mem_max);
    phgOptionsRegisterInt("-nmax", "Maximum number of elements", &nmax);
    phgOptionsRegisterNoArg("-export_mesh", "Export mesh", &export_mesh);
    phgOptionsRegisterNoArg("-export_all", "Export all meshes", &export_all);
    phgOptionsRegisterKeyword("-export_format", "Format of Export files",
			      export_names, &export_index);
    phgOptionsRegisterNoArg("-dump_matrix", "Dump matrix", &dump_matrix);
    phgOptionsRegisterNoArg("-dump_rhs", "Dump RHS", &dump_rhs);
    phgOptionsRegisterNoArg("-L2adapt", "Adaptation using L2 error indicator",
			    &L2adapt);
    phgOptionsRegisterNoArg("-use_weights", "Use error indicator as weights",
			    &use_weights);
    phgOptionsRegisterNoArg("-hp_test", "Test h-p adaptivity functions",
			    &hp_test);
    phgOptionsRegisterNoArg("-flatten_grid",
			    "Flatten the grid after refinement", &flatten);
    phgOptionsRegisterNoArg("-init0",
			    "Set initial solution of the linear "
			    "systems to 0", &init0);
    phgOptionsRegisterNoArg("-show_cond", "Show condition number of the "
			    "stiffness matrix", &show_cond);
    phgOptionsRegisterNoArg("-qcache", "Numerical quadrature interface to use. "
			    "True: QCache, False: legacy phgQuad* functions",
			    &qcache);

    /*------------------- End command-line options ----------------------*/

    phgOptionsPreset("-dof_type P3");
    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    /* pre-refine */
    if (pre_refines > 3) {
	INT k;
	k = g->nelem_global;
#if USE_MPI
	MPI_Bcast(&k, 1, PHG_MPI_INT, 0, phgComm);
#endif	/* USE_MPI */
	k *= Pow(2.0, pre_refines);
	if (0.5 * k > (submesh_threshold * lif_threshold) * phgNProcs) {
	    /* parallel refinement with load rebalancing */
	    while (pre_refines > 0) {
		k = pre_refines > 3 ? 3 : pre_refines;
		phgRefineAllElements(g, k);
#if 0
		/* FIXME: program hangs in the phgBalanceGrid call below
		 * if (g->nprocs > 1 && g->nprocs < phgNProcs) occurs */
		phgBalanceGrid(g, lif_threshold, submesh_threshold, NULL, 0.);
#else
		phgRedistributeGrid(g);
#endif
		pre_refines -= k;
		phgPrintf("Pre-refine: %d procs, %"dFMT" elements.\n",
				g->nprocs, g->nelem_global);
	    }
	    assert(g->nprocs == phgNProcs);
	}
    }
    if (pre_refines > 0)
	phgRefineAllElements(g, pre_refines);

    /* The discrete solution u_h */
    if (hp_test) {
	/* For testing the h-p code */
	ELEMENT *e;
	HP_TYPE *hp = phgHPNew(g, HP_HB);
	ForAllElements(g, e)
	    e->hp_order = DOF_DEFAULT->order + GlobalElement(g, e->index) % 3;
	phgHPSetup(hp, FALSE);
	u_h = phgHPDofNew(g, hp, 1, "u_h",
			init0 ? DofNoAction : DofInterpolation);
	phgHPFree(&hp);
    }
    else {
	u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h",
			init0 ? DofNoAction : DofInterpolation);
    }
    u_h->DB_mask = BDRY_MASK;	/* Dirichlet BC on all boundaries */

    /* The source function (RHS) f_h */
#if 0
    f_h = phgDofNew(g, DOF_ANALYTIC, 1, "f_h",  func_f);
#else
    f_h = phgDofClone(u_h, "f_h", func_f);
#endif

    /* The L2/H1 error indicators */
    e_L2 = phgDofNew(g, DOF_P0, 1, "L2 error indicator", DofNoAction);
    e_H1 = phgDofNew(g, DOF_P0, 1, "H1 error indicator", DofNoAction);

    phgPrintf("Test case %d, mesh file: \"%s\"\n", TEST_CASE, fn);
    phgPrintf("Analytic solution: \"%s\"\n", analytic);
    phgPrintf("Approximation: u_h = %s, f_h = %s\n\n",
		DofTypeName(u_h), DofTypeName(f_h));
    phgOptionsShowUsed();		/* show options */

    while (TRUE) {
	elapsed_time(g, FALSE, 0.);	/* reset timer */
	phgPrintf("\n*** Pass %d\n", g->serial_no % 1000);
#if DEBUG
	phgPrintf("--- Check grid: ");
	phgCheckConformity(g);
	elapsed_time(g, TRUE, 0.);
#endif
	phgPrintf("DOF: %"dFMT", T: %"dFMT", V: %"dFMT", E: %"dFMT
		  ", F: %"dFMT"\n",
		  DofDataCountGlobal(u_h),
		  g->nleaf_global, g->nvert_global, g->nedge_global,
		  g->nface_global);
	if (!use_weights && phgBalanceGrid(g, lif_threshold, submesh_threshold,
						 NULL, 0.)) {
	    phgPrintf("Redistribute grid ");
	    elapsed_time(g, TRUE, 0.);
#if DEBUG
	    phgPrintf("--- Check grid: ");
	    phgCheckConformity(g);
	    elapsed_time(g, TRUE, 0.);
#endif
	}
	phgPrintf("Grid distribution: %d submesh%s, LIF = %0.3lg\n",
			g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);

	phgPrintf("Prepare solver: ");
	if ((solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL)) == NULL) {
	    phgPrintf("abort.\n");
	    break;
	}
	phgPrintf("LIF %lg ", (double)solver->mat->rmap->lif);
	elapsed_time(g, TRUE, 0.);
    
	phgPrintf("Build linear system: ");
	phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	build_linear_system(solver, u_h, f_h);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

#if 0
/* test block matrix */
{
    MAT *a = phgMatCreateBlockMatrix(g->comm, 1, 1, &solver->mat, NULL, NULL);
    phgMatDestroy(&solver->mat);
    solver->rhs->mat = solver->mat = a;
}
#endif

#if 0
/* test solving multiple RHS with a same matrix */
{
    INT nlocal = solver->rhs->map->nlocal;
    INT osize = solver->rhs->map->localsize - nlocal;
    FLOAT *rhs = phgAlloc(nlocal * sizeof(*rhs));
    FLOAT *rhs_offp = phgAlloc(osize * sizeof(*rhs_offp));
    FLOAT *data = phgAlloc(DofDataCount(u_h) * sizeof(*u_h->data));
    /* save (unassembled) RHS and initial solution */
    memcpy(rhs, solver->rhs->data, nlocal * sizeof(*rhs));
    if (osize > 0)
	memcpy(rhs_offp, solver->rhs->offp_data, osize * sizeof(*rhs_offp));
    memcpy(data, u_h->data, DofDataCount(u_h) * sizeof(*u_h->data));
    /*bzero(solver->rhs, localsize * sizeof(*rhs));*/
    phgSolverSolve(solver, FALSE, u_h, NULL);
    /* reset RHS */
    phgSolverResetRHS(solver);
    memcpy(solver->rhs->data, rhs, nlocal * sizeof(*rhs));
    if (osize > 0)
	memcpy(solver->rhs->offp_data, rhs_offp, osize * sizeof(*rhs_offp));
    memcpy(u_h->data, data, DofDataCount(u_h) * sizeof(*u_h->data));
    phgFree(rhs);
    phgFree(rhs_offp);
    phgFree(data);
}
#endif

	phgPrintf("Assemble solver: ");
	phgSolverAssemble(solver);
	phgSolverUpdateRHS(solver);
	elapsed_time(g, TRUE, 0.);

	if (dump_matrix || dump_rhs) {
	    static int no = 0;
	    char fn[128];

#if USE_MPI
	    MPI_Bcast(&no, 1, MPI_INT, 0, g->comm);
#endif	/* USE_MPI */
	    no++;
#if 0
	    if (dump_matrix) {
		sprintf(fn, "matrix%d.dat", no);
		phgMatDumpCSR(solver->mat, fn);
	    }
	    if (dump_rhs) {
		sprintf(fn, "rhs%d.dat", no);
		phgVecDump(solver->rhs, fn);
	    }
#else
	    if (dump_matrix) {
		sprintf(fn, "matrix%d.m", no);
		phgMatDumpMATLAB(solver->mat, "A", fn);
	    }
	    if (dump_rhs) {
		sprintf(fn, "rhs%d.m", no);
		phgVecDumpMATLAB(solver->rhs, "b", fn);
	    }
#endif
	}

#if 0 && USE_PETSC
	if (solver->oem_solver == SOLVER_PETSC) {
	    /* test setting customized PETSc options */
	    void **ksp = solver->oem_solver->GetSolver(solver);
	    PetscOptionsSetValue("-pc_type", "none");
	    KSPSetFromOptions(*ksp);
	}
#endif	/* USE_PETSC */

	if (show_cond) {
	    /* diagonal scaling */
	    INT i, i0;
	    MAT *D, *tmp, *A;
	    phgMatSetupDiagonal(A = solver->mat);
	    D = phgMapCreateMat(A->rmap, A->cmap);
	    i0 = D->rmap->partition[g->rank];
	    for (i = 0; i < D->rmap->nlocal; i++)
		phgMatAddGlobalEntry(D, i + i0, i + i0,
						1.0 / Sqrt(Fabs(A->diag[i])));
	    phgMatAssemble(D);
	    tmp = phgMatMat(MAT_OP_N, MAT_OP_N, 1.0, D, A, 0.0, NULL);
	    A = phgMatMat(MAT_OP_N, MAT_OP_N, 1.0, tmp, D, 0.0, NULL);
	    phgMatDestroy(&D);
	    phgMatDestroy(&tmp);
	    phgPrintf("Condition number of the stiffness matrix "
			"(after diagonal scaling): %0.2le\n",
			  (double)phgMatConditionNumber(A));
	    phgMatDestroy(&A);
	}

	phgPrintf("Solve linear system: ");
	if (init0)
	    phgDofSetDataByValue(u_h, 0.0);
	phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	if ((nits = phgSolverSolve(solver, TRUE, u_h, NULL)) < 0) {
	    phgPrintf("abort.\n");
	    break;
	}
	phgPrintf("%d its, %0.4le ", nits, (double)solver->residual);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	phgSolverDestroy(&solver);

	phgPrintf("Estimate: ");
	grad_u = phgDofGradient(u_h, NULL, NULL, "grad u_h");
	estimate_error(u_h, f_h, grad_u, e_L2, e_H1);
	EL2_error = Sqrt(phgDofNormL1Vec(e_L2));
	EH1_error = Sqrt(phgDofNormL1Vec(e_H1));
	phgPrintf("EL2=%0.2le EH1=%0.2le ", (double)EL2_error,
					    (double)EH1_error);
	phgMemoryUsage1(g, &current_min, &current_avg, &current_max,
				 &peak_min, &peak_avg, &peak_max);
	flag = (EH1_error < tol || peak_max > 1024 * (size_t)mem_max * 1024 ||
		(nmax > 0 && g->nelem_global >= nmax));

	if (has_analytic_solution) {
	    DOF *u, *grad;
	    u = phgDofClone(u_h, "analytic solution", func_u);
	    grad = phgDofGradient(u, NULL, NULL, "analytic gradient");
	    phgDofSetDataByFunction(grad, func_grad);
	    phgDofAXPY(-1.0, u_h, &u);
	    /* Note: the following works only if u is smooth! */
	    phgDofAXPY(-1.0, grad_u, &grad);
	    L2_error = phgDofNormL2(u);
	    H1_error = phgDofNormL2(grad);
	    phgPrintf("L2=%0.2le H1=%0.2le ", (double)L2_error,
			    		      (double)H1_error);
	    phgDofFree(&grad);
	    phgDofFree(&u);
	}
	elapsed_time(g, TRUE, 0.);

	phgPrintf("Memory usage: current %0.4lg/%0.4lg/%0.4lgMB, "
		  "peak %0.4lg/%0.4lg/%0.4lgMB\n",
		(double)current_min / (1024.0 * 1024.0),
		(double)current_avg / (1024.0 * 1024.0),
		(double)current_max / (1024.0 * 1024.0),
		(double)peak_min / (1024.0 * 1024.0),
		(double)peak_avg / (1024.0 * 1024.0),
		(double)peak_max / (1024.0 * 1024.0));

    	if (export_all || (flag && export_mesh)) {
	    static char fn[64];
	    const char *p;
	    if (export_all)
		sprintf(fn, "mesh%03d.%s", g->serial_no % 1000,
				export_ext[export_index]);
	    else
		sprintf(fn, "mesh.%s", export_ext[export_index]);
	    phgPrintf("Export %s file %s ",
			    export_names[export_index],
			    p = export_list[export_index](g, fn, u_h, grad_u,
							e_L2, e_H1, NULL));
            if (phgRank == 0) {
		FILE *f;
		long long size = 0L;
		if ((f = fopen(p, "r")) != NULL) {
		    fseek(f, 0L, SEEK_END);
		    size = ftell(f);
		    fclose(f);
		}
		phgPrintf("size = %lld ", size);
	    }
	    elapsed_time(g, TRUE, 0.);
	    /*phgExportALBERT(g, "mesh.dat");*/
    	}

	phgDofFree(&grad_u);

#define Idle(g)	((g) == NULL || phgRank >= (g)->nprocs)
	if (flag && !Idle(g))
	    break;

	if (use_weights && phgBalanceGrid(g, lif_threshold, submesh_threshold,
				L2adapt ? e_L2 : e_H1, 1.)) {
	    phgPrintf("Redistribute grid, %d submeshes, LIF = %lg ",
			g->nprocs, (double)g->lif);
	    elapsed_time(g, TRUE, 0.);
#if DEBUG
	    phgPrintf("--- Check grid: ");
	    phgCheckConformity(g);
	    elapsed_time(g, TRUE, 0.);
#endif
	}

	phgPrintf("Refine mesh: ");
	if (refine_depth <= 0) {
	    /* uniform refine */
	    phgRefineAllElements(g, refine_depth == 0 ? 3 : -refine_depth);
	}
	else {
	    DOF *ei = L2adapt ? e_L2 : e_H1;
	    /* adaptive refine */
	    if (hp_test) {
		/* test changing orders of u_h */
		ELEMENT *e;
		ForAllElements(g, e)
		    e->hp_order = DOF_DEFAULT->order + 
					GlobalElement(g, e->index) % 5;
		phgHPSetup(u_h->hp, TRUE);
	    }
	    phgMarkRefine(MARK_DEFAULT, ei, Pow(0.8, 2), NULL, 0.0,
				refine_depth, Pow(tol, 2) / g->nleaf_global);
	    phgRefineMarkedElements(g);
	}
	elapsed_time(g, TRUE, 0.);

	if (flatten) {
	    /* flatten the mesh */
	    GRID *g1 = phgDupGrid(g, TRUE);
	    /* TODO: copy u_h data */
	    assert(hp_test == FALSE);		/* TODO */
	    phgDofFree(&u_h);
	    phgDofFree(&f_h);
	    phgDofFree(&e_L2);
	    phgDofFree(&e_H1);
	    u_h = phgDofNew(g1, DOF_DEFAULT, 1, "u_h", DofNoAction);
	    u_h->DB_mask = BDRY_MASK;	/* Dirichlet BC on all boundaries */
	    f_h = phgDofClone(u_h, "f_h", func_f);
	    e_L2 = phgDofNew(g1, DOF_P0, 1, "L2 error indicator", DofNoAction);
	    e_H1 = phgDofNew(g1, DOF_P0, 1, "H1 error indicator", DofNoAction);
	    phgFreeGrid(&g);
	    g = g1;
	}

	if (!init0 && has_analytic_solution) {
	    DOF *u0 = phgDofClone(u_h, "u0", func_u);
	    phgDofAXPY(-1.0, u_h, &u0);
	    phgPrintf("Interpolation error: %0.3le\n",
			    (double)phgDofNormInftyVec(u0));
	    phgDofFree(&u0);
	}
    }

    phgDofFree(&u_h);
    phgDofFree(&f_h);
    phgDofFree(&e_L2);
    phgDofFree(&e_H1);
    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
