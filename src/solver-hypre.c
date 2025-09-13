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

/* $Id: solver-hypre.c,v 1.156 2022/09/21 02:35:28 zlb Exp $ */

/* The HYPRE interface for PHG */

#include "phg.h"

#if USE_HYPRE	/* whole file */

/* HAVE_HYPRE_AMS_SET_PI controls whether explicitly compute the matrix Pi.
 * Note: HAVE_HYPRE_AMS_SET_PI != 0 requires the HYPRE_AMSSetDiscretePi
 * function which does not exist in the standard distribution and is provided
 * by the non oficial patch 'hypre-ams-pi-patch' */
#define HAVE_HYPRE_AMS_SET_PI	0

#include <string.h>
#include <math.h>
#include <HYPRE_utilities.h>
#include <krylov.h>
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>

/* Note: HYPRE_Int is defined since version 2.7.0b */
#if HYPRE_VERSION_MAJOR < 2 || \
    (HYPRE_VERSION_MAJOR == 2 && HYPRE_VERSION_MINOR < 7)
typedef int HYPRE_Int;
#endif

#define g_comm	map->comm

typedef enum {
    NONE, DIAG, PCG, ParaSails, GMRES, BiCGSTAB, Euclid, BoomerAMG,
#if HYPRE_VERSION_MAJOR >= 2
    AMS,
#if HYPRE_VERSION_MINOR >= 4
    LGMRES, FGMRES,
#endif	/* HYPRE_VERSION_MINOR >= 4 */
#endif	/* HYPRE_VERSION_MAJOR >= 2 */
    MLI
} SOLVER_ID;

static const char *solver_names[] = {
    "pcg", "gmres", "bicgstab", "boomeramg",
#if HYPRE_VERSION_MAJOR >= 2
    "ams",
#if HYPRE_VERSION_MINOR >= 4
    "lgmres", "fgmres",
#endif	/* HYPRE_VERSION_MINOR >= 4 */
#endif	/* HYPRE_VERSION_MAJOR >= 2 */
    NULL
};
static const SOLVER_ID solver_list[] = {
    PCG, GMRES, BiCGSTAB, BoomerAMG,
#if HYPRE_VERSION_MAJOR >= 2
    AMS,
#if HYPRE_VERSION_MINOR >= 4
    LGMRES,FGMRES
#endif	/* HYPRE_VERSION_MINOR >= 4 */
#endif	/* HYPRE_VERSION_MAJOR >= 2 */
};

static const char *precon_names[] = {
    "none", "diag", "parasails", "euclid", "boomeramg",
#if HYPRE_VERSION_MAJOR >= 2
    "ams",
#if HYPRE_VERSION_MINOR >= 4
    "lgmres",
#endif	/* HYPRE_VERSION_MINOR >= 4 */
#endif	/* HYPRE_VERSION_MAJOR >= 2 */
    NULL
};
static const SOLVER_ID precon_list[] = {
    NONE, DIAG, ParaSails, Euclid, BoomerAMG,
#if HYPRE_VERSION_MAJOR >= 2
    AMS,
#if HYPRE_VERSION_MINOR >= 4
    LGMRES, FGMRES,
#endif	/* HYPRE_VERSION_MINOR >= 4 */
#endif	/* HYPRE_VERSION_MAJOR >= 2 */
};

#if HYPRE_VERSION_MAJOR <= 1 && HYPRE_VERSION_MINOR < 9
# define HYPRE_IJVectorSetMaxOffProcElmts(a, b)
# define HYPRE_PCGSetPrintLevel(a, b)
# define HYPRE_ParCSRGMRESSetPrintLevel(a, b)
# define HYPRE_ParCSRBiCGSTABSetPrintLevel(a, b)
# define HYPRE_BoomerAMGSetPrintLevel(a, b)
# define HYPRE_BoomerAMGSetRelaxType(a, b)
# define HYPRE_BoomerAMGSetCycleRelaxType(a, b, c)
# define HYPRE_BoomerAMGSetNumSweeps(a, b)
#endif	/* HYPRE_VERSION_MAJOR <= 1 && HYPRE_VERSION_MINOR < 9 */

static int amg_coarsen_types[] = {
    0,	/* CLJP */
    3,	/* RS on each processor followed by boundary coarsening */
    6,	/* Falgout */
    8,	/* PMIS */
    10,	/* HMIS */
};

static const char *amg_coarsen_names[] = {
  "cljp", "rs+bc", "falgout", "pmis", "hmis", NULL
};
static const char *amg_cycle_names[] = {"default", "v-cycle", "w-cycle", NULL};

static const char *sai_sym_list[] = {"nonsym-indefinite", "spd",
					"nonsym-definite", NULL};

static const char *amg_relax_names[] = {
    "jacobi",		/* 0: Jacobi */
    "gs-seq",		/* 1: GS, sequential (very slow!) */
    "gs",		/* 2: GS, interior parallel, boundary seq. (slow!) */
    "gs-h-forward",	/* 3: hybrid GS or SOR, forward solve */
    "gs-h-backward",	/* 4: hybrid Gauss-Seidel or SOR, backward solve */
    "gs-h-chaotic",	/* 5: hybrid chaotic Gauss-Seidel (OpenMP only) */
    "gs-h-symmetric",	/* 6: hybrid symmetric Gauss-Seidel or SSOR */
#if HYPRE_VERSION_MAJOR > 2 || \
    (HYPRE_VERSION_MAJOR == 2 && HYPRE_VERSION_MINOR > 9)
    "l1-gs",		/* 8: l1-scaled hybrid symmetric Gauss-Seidel */
    "l1-jacobi",	/* 18: l1-scaled jacobi */
#endif	/* l1-scaled smoother */
    "direct",		/* 9: Gaussian elimination (only on coarsest level) */
    NULL
};

#if HYPRE_VERSION_MAJOR > 2 || \
    (HYPRE_VERSION_MAJOR == 2 && HYPRE_VERSION_MINOR > 9)
static const int amg_relax_types[] = {
	0, 1, 2, 3, 4, 5, 6, 8, 18, 9};
#else
static const int amg_relax_types[] = {0, 1, 2, 3, 4, 5, 6, 9};
#endif

#if HYPRE_VERSION_MAJOR >= 2
static const char *cycle_names[] = {"01210", "0+1+2", "02120", "010+2",
	"0102010", "1+020", "0201020", "010+020",
#if HYPRE_VERSION_MINOR >= 2
	"013454310", "0+1+3+4+5", "034515430", "01310+01410+01510",
#endif
       	NULL};
#endif	/* HYPRE_VERSION_MAJOR >= 2 */

/* cmdline arguments to be passed to Euclid solver */
static int argc = 0;
static char **argv = NULL;

typedef struct {
    int		solver_id;		/* default solver = PCG */
    int		precon_id;		/* default preconditioner = ParaSails */
    FLOAT	pc_rtol;		/* preconditioner rtol */
    FLOAT	pc_atol;		/* preconditioner atol */
    INT		pc_maxit;		/* preconditioner maxit */

    /* Note: the output matrices/vectors can be used as input for ams_driver.c
     * from the HYPRE distribution. For example,
     *	maxwell -pre_refines=2 -tol=10 +hypre_two_norm -hypre_dump_matvec=aFEM
     * and
     *	ams_driver -coord
     * should produce identical results. */
    char	*dump_matvec;
    BOOLEAN	two_norm;		/* whether use two norm */
    INT		gmres_kdim;		/* <=0 means default */

    /* cmdline arguments for boomerAMG solver */
    INT		amg_num_funcs;		/* size of the system of PDEs */
    INT		amg_max_levels;		/* max MG levels */
    FLOAT	amg_strength;		/* strength threshold */
    FLOAT	amg_max_row_sum;	/* max row sum */

    int		amg_coarsen_type;	/* default coarsening type = Falgout */
    int		amg_cycle_type;		/* MG cycle type */
    int		amg_relax_type;		/* relaxation type */
    int		amg_coarsest_relax_type;/* relax type on the coarsest grid */

#if HYPRE_VERSION_MAJOR >= 2
    /* default parameters for AMS solver */
    int		cycle_type;	/* 3-level cycle type - 1 */
    INT		ams_dim;		/* space dimension */
    INT		rlx_type;		/* relaxation type */
    INT		rlx_sweeps;		/* number of relaxation sweeps */
    INT		ams_coarsen_type;	/* BoomerAMG coarsening type */
    INT		ams_agg_levels;		/* Levels of BoomerAMG agg. coarsening*/
    INT		ams_rlx_type;		/* BoomerAMG relaxation type */
    FLOAT	theta;			/* BoomerAMG threshold */
    BOOLEAN	ams_singular;		/* singular problem (beta=0) */
#endif	/* HYPRE_VERSION_MAJOR >= 2 */

    /* default parameters for Parasails-PCG solver */
    INT		sai_max_levels;
    FLOAT	sai_threshold;
    FLOAT	sai_filter;
    int		sai_sym;
} PARAMETERS;

static PARAMETERS global_params = {
    /* solver_id */		0,
    /* precon_id */		1,
    /* pc_rtol */		0.0,
    /* pc_atol */		0.0,
    /* pc_maxit */		1,
    /* dump_matvec */		NULL,
    /* two_norm */		TRUE,
    /* gmres_kdim */		0/*5*/,

    /* amg_num_funcs */		-1,
    /* amg_max_levels */	-1,
    /* amg_strength */		0.5,
    /* amg_max_row_sum */	0.9,

    /* amg_coarsen_type */	3,	/* pmis */
    /* amg_cycle_type */	0,
    /* amg_relax_type */	6,	/* gs-h-symmetric */
    /* amg_coarsest_relax_type */ 7,

#if HYPRE_VERSION_MAJOR >= 2
    /* cycle_type */		0,
    /* ams_dim */		3,
    /* rlx_type */		2,
    /* rlx_sweeps */		1,
    /* ams_coarsen_type */	10,
    /* ams_agg_levels */	1,
    /* ams_rlx_type */		3,
    /* theta */			0.25,	
    /* ams_singular */		FALSE,
#endif	/* HYPRE_VERSION_MAJOR >= 2 */

    /* sai_max_levels */	1,
    /* sai_threshold */		0.1,
    /* sai_filter */		0.05,
    /* sai_sym */		1
};

typedef int (*SolveFcn)(HYPRE_Solver solver, HYPRE_ParCSRMatrix A,
			HYPRE_ParVector x, HYPRE_ParVector b);
typedef int (*SetPCFcn)(HYPRE_Solver solver, HYPRE_PtrToSolverFcn solve,
			HYPRE_PtrToSolverFcn setup, HYPRE_Solver precon);
typedef int (*GetItFcn)(HYPRE_Solver solver, HYPRE_Int *nits);
typedef int (*GetResidFcn)(HYPRE_Solver solver, double *resid);
typedef int (*DestroyFcn)(HYPRE_Solver solver);
typedef int (*SetTolFcn)(HYPRE_Solver solver, double tol);
typedef int (*SetMaxIterFcn)(HYPRE_Solver solver, HYPRE_Int maxiter);

#if HYPRE_VERSION_MAJOR >= 2
typedef struct{
    HYPRE_IJMatrix G;			/* the discrete gradient matrix */
    HYPRE_IJMatrix Pi;			/* the Hcurl=>(H1)^3 matrix */
    HYPRE_IJMatrix A_alpha;		/* alpha possion matrix */
    HYPRE_IJMatrix A_beta;		/* beta poisson matrix */
    HYPRE_IJVector x, y, z;		/* coordinates of vertices */
} AMS_DATA;
#endif	/* HYPRE_VERSION_MAJOR >= 2 */

typedef struct {
    HYPRE_Int		ilower, iupper;
    HYPRE_IJMatrix	A;		/* The coefficient matrix */
    HYPRE_IJVector	B, U;		/* RHS & solution */
    BOOLEAN		assembled;
    HYPRE_Solver	hsolver, precond;

    /* function pointers to be initialized by setup_hypre_solver() */
    SolveFcn	solver_setup, precon_setup; 
    SolveFcn	solver_solve, precon_solve; 
    SetTolFcn	solver_set_tol;
    SetTolFcn	solver_set_atol;
    SetMaxIterFcn solver_set_maxiter;
    DestroyFcn	solver_destroy, precon_destroy;
    SetPCFcn	solver_setpc;
    GetItFcn	solver_get_nits;
    GetResidFcn solver_get_resid;

    PARAMETERS	params;

#if HYPRE_VERSION_MAJOR >= 2
    /* AMS data */
    AMS_DATA	*ams_data;
#endif	/* HYPRE_VERSION_MAJOR >= 2 */
} OEM_DATA;

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)
#define _ams	(_t->ams_data)
#define _prm	(_t->params)

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nHYPRE options:", "\n", "hypre");
    phgOptionsRegisterKeyword("hypre_solver", "HYPRE solver", solver_names,
				&global_params.solver_id);
    phgOptionsRegisterKeyword("hypre_pc", "HYPRE preconditioner", precon_names,
				&global_params.precon_id);
    phgOptionsRegisterFloat("hypre_pc_rtol", "Preconditioner rtol",
				&global_params.pc_rtol);
    phgOptionsRegisterFloat("hypre_pc_atol", "Preconditioner atol",
				&global_params.pc_atol);
    phgOptionsRegisterInt("hypre_pc_maxit", "Preconditioner maxit",
				&global_params.pc_maxit);
    phgOptionsRegisterNoArg("hypre_two_norm", "Use two-norm (PCG)",
				&global_params.two_norm);
    phgOptionsRegisterFilename("hypre_dump_matvec", "Dump HYPRE mat/vecs",
				&global_params.dump_matvec);

    /* cmdline options for Parasails-PCG */
    phgOptionsRegisterTitle("\n[Options for GMRES|LGMRES|FGMRES solver]",
				"\n", "hypre");
    phgOptionsRegisterInt("hypre_gmres_kdim",
			  "Maxi. dimension of Krylov space (<=0 means default)",
				&global_params.gmres_kdim);

    /* cmdline options for boomerAMG (incomplete) */
    phgOptionsRegisterTitle("\n[Options for BoomerAMG solver]", "\n", "hypre");
    phgOptionsRegisterInt("hypre_amg_num_functions",
				"Size of the system of PDEs",
				&global_params.amg_num_funcs);
    phgOptionsRegisterInt("hypre_amg_max_levels", "Max multigrid levels",
				&global_params.amg_max_levels);
    phgOptionsRegisterFloat("hypre_amg_strength", "Strength threshold",
				&global_params.amg_strength);
    phgOptionsRegisterFloat("hypre_amg_max_row_sum", "Max row sum",
				&global_params.amg_max_row_sum);
    phgOptionsRegisterKeyword("hypre_amg_coarsen_type", "Coarsening type",
				amg_coarsen_names,
				&global_params.amg_coarsen_type);
    phgOptionsRegisterKeyword("hypre_amg_cycle_type", "Multigrid cycle type",
				amg_cycle_names,
				&global_params.amg_cycle_type);
    phgOptionsRegisterKeyword("hypre_amg_relax_type", "Multigrid smoother",
				amg_relax_names,
				&global_params.amg_relax_type);
    phgOptionsRegisterKeyword("hypre_amg_coarsest_relax_type",
				"Coarsest grid solver",
				amg_relax_names,
				&global_params.amg_coarsest_relax_type);

    /* cmdline options for Parasails-PCG */
    phgOptionsRegisterTitle("\n[Options for ParaSails solver]", "\n", "hypre");
    phgOptionsRegisterInt("hypre_sai_max_levels", "Max levels",
				&global_params.sai_max_levels);
    phgOptionsRegisterFloat("hypre_sai_threshold", "Threshold",
				&global_params.sai_threshold);
    phgOptionsRegisterFloat("hypre_sai_filter", "Filter",
				&global_params.sai_filter);
    phgOptionsRegisterKeyword("hypre_sai_sym", "Symmetry", sai_sym_list,
				&global_params.sai_sym);
#if HYPRE_VERSION_MAJOR >= 2
    /* cmdline options for AMS */
    phgOptionsRegisterTitle("\n[Options for AMS solver]", "\n", "hypre");
    phgOptionsRegisterKeyword("hypre_ams_type", "Three-level cycle type",
				cycle_names,
				&global_params.cycle_type);
    phgOptionsRegisterInt("hypre_ams_dim", "Space dimension",
				&global_params.ams_dim);
    phgOptionsRegisterFloat("hypre_ams_theta", "BoomerAMG threshold",
				&global_params.theta);
    phgOptionsRegisterInt("hypre_ams_rlx", "Relaxation type",
				&global_params.rlx_type);
    phgOptionsRegisterInt("hypre_ams_rlxn", "Number of relaxation sweeps",
				&global_params.rlx_sweeps);
    phgOptionsRegisterInt("hypre_ams_ctype", "BoomerAMG coarsening type",
				&global_params.ams_coarsen_type);
    phgOptionsRegisterInt("hypre_ams_agg",
				"Levels of BoomerAMG agg. coarsening",
				&global_params.ams_agg_levels);
    phgOptionsRegisterInt("hypre_ams_amgrlx", "BoomerAMG relaxation type",
				&global_params.ams_rlx_type);
    phgOptionsRegisterNoArg("hypre_ams_singular",
				"Set the singular flag (beta=0)",
				&global_params.ams_singular);
#endif	/* HYPRE_VERSION_MAJOR >= 2 */

    /* TODO: cmdline options for setting parameters for PCG, AMG, etc. */

    return 0;
}

static int
Initialize(int *argcp, char ***argvp)
/* duplcate command-line arguments (for use with Euclid preconditioner) */
{
    int i;

    if (phgVerbosity > 0 && phgRank == 0)
	printf("sizeof(HYPRE_Int) = %d\n", (int)sizeof(HYPRE_Int));

    /* TODO:
     *	1. only copy options handled by Euclid, and remove them from
     *	   the arguments list? 
     *	2. change '-oem_options' to '-petsc_options' and '-hypre_options'?
     */

    argv = phgAlloc(((argc = *argcp) + 1) * sizeof(*argv));
    for (i = 0; i < argc; i++)
	argv[i] = strdup((*argvp)[i]);
    argv[i] = NULL;

    return 0;
}

static int
Finalize(void)
{
    int i;

    for (i = 0; i < argc; i++)
	phgFree(argv[i]);

    phgFree(argv);

    argc = 0;
    argv = NULL;

    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));
    memcpy(&_prm, &global_params, sizeof(PARAMETERS));
    if (global_params.dump_matvec != NULL)
	_prm.dump_matvec = strdup(global_params.dump_matvec);

    return 0;
}

static int
Create(SOLVER *solver)
{
    MAP *map = solver->mat->rmap;
#if HYPRE_VERSION_MAJOR >= 2
	HYPRE_Int *d_nnz, *o_nnz;
	INT Nlocal, i, d, o;
#endif

    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }
    assert(solver->oem_data != NULL);

    _t->ilower = map->partition[map->rank];
    _t->iupper = map->partition[map->rank + 1] - 1;


    _t->A = NULL;
    _t->B = _t->U = NULL;
    _t->assembled = FALSE;

#if HYPRE_VERSION_MAJOR >= 2
    _ams = NULL;

#endif	/* HYPRE_VERSION_MAJOR >= 2 */
    /* create matrix */
    HYPRE_IJMatrixCreate(g_comm, _t->ilower, _t->iupper, _t->ilower,
			 _t->iupper, &_t->A);
    HYPRE_IJMatrixSetObjectType(_t->A, HYPRE_PARCSR);

#if HYPRE_VERSION_MAJOR >= 2
    /* PreAlloc for matrix */
    Nlocal = _t->iupper - _t->ilower + 1;
    d_nnz = phgAlloc(Nlocal * sizeof(*d_nnz));
    o_nnz = phgAlloc(Nlocal * sizeof(*o_nnz));
    for (i = 0; i < Nlocal; i++) {
	phgMatGetNnzInRow(solver->mat, i, &d, &o);
	d_nnz[i] = d;
	o_nnz[i] = o;
    }
    HYPRE_IJMatrixSetDiagOffdSizes(_t->A, d_nnz, o_nnz);
    phgFree(d_nnz);
    phgFree(o_nnz);
#endif

    HYPRE_IJMatrixInitialize(_t->A);

    return 0;
}

static int
AddMatrixEntries(SOLVER *solver, INT nrows, INT *rows, INT ncols, INT *cols,
			FLOAT *values)
{
    INT i;
    HYPRE_Int *pcols, nc = ncols, r;

    if (nrows == 0 || ncols == 0)
	return 0;

    if (sizeof(INT) != sizeof(HYPRE_Int)) {
	pcols = phgAlloc(ncols * sizeof(*pcols));
	for (i = 0; i < ncols; i++)
	    pcols[i] = (HYPRE_Int)cols[i];
    }
    else {
	pcols = (HYPRE_Int *)cols;
    }

#if FT_PHG == FT_DOUBLE
    for (i = 0; i < nrows; i++) {
	r = (HYPRE_Int)*(rows++);
	HYPRE_IJMatrixSetValues(_t->A, 1, &nc, &r, pcols, values);
	values += ncols;
    }
#else	/* FT_PHG == FT_DOUBLE */
    {
	INT j;
	double *v = phgAlloc(ncols * sizeof(*v));
	for (i = 0; i < nrows; i++) {
	    for (j = 0; j < ncols; j++)
		v[j] = *(values++);
	    r = (HYPRE_Int)*(rows++);
	    HYPRE_IJMatrixSetValues(_t->A, 1, &nc, &r, pcols, v);
	    values += ncols;
	}
	phgFree(v);
    }
#endif	/* FT_PHG == FT_DOUBLE */

    if (sizeof(INT) != sizeof(HYPRE_Int))
	phgFree(pcols);

    return 0;
}

static int
AddRHSEntries(SOLVER *solver, INT n, INT *ni, FLOAT *values)
{
    MAP *map = solver->mat->rmap;
    INT i;
    HYPRE_Int *pni;

    if (_t->B == NULL) {
	/* create RHS vector */
	HYPRE_IJVectorCreate(g_comm, _t->ilower, _t->iupper, &_t->B);
	HYPRE_IJVectorSetObjectType(_t->B, HYPRE_PARCSR);
	HYPRE_IJVectorSetMaxOffProcElmts(_t->B, 0);
	HYPRE_IJVectorInitialize(_t->B);
    }

    if (n <= 0)
	return 0;

    if (sizeof(INT) != sizeof(HYPRE_Int)) {
	pni = phgAlloc(n * sizeof(*pni));
	for (i = 0; i < n; i++)
	    pni[i] = (HYPRE_Int)ni[i];
    }
    else {
	pni = (HYPRE_Int *)ni;
    }

#if FT_PHG == FT_DOUBLE
    HYPRE_IJVectorSetValues(_t->B, n, pni, values);
#else	/* FT_PHG == FT_DOUBLE */
    {
	INT j;
	double *v = phgAlloc(n * sizeof(*v));
	for (j = 0; j < n; j++)
	    v[j] = *(values++);
	HYPRE_IJVectorSetValues(_t->B, n, pni, v);
	phgFree(v);
    }
#endif	/* FT_PHG == FT_DOUBLE */

    if (sizeof(INT) != sizeof(HYPRE_Int))
	phgFree(pni);

    return 0;
}

static void
dump_vector(SOLVER *solver, HYPRE_ParVector vec, const char *name)
{
    char fn[128];

    if (_prm.dump_matvec == NULL)
	return;
    sprintf(fn , "%s.%s", _prm.dump_matvec, name);
    phgInfo(1, "Dumping HYPRE ParVector to file \"%s\"\n", fn);
    HYPRE_ParVectorPrint(vec, fn);
}

static void
dump_matrix(SOLVER *solver, HYPRE_ParCSRMatrix mat, const char *name)
{
    char fn[128];

    if (_prm.dump_matvec == NULL)
	return;
    sprintf(fn , "%s.%s", _prm.dump_matvec, name);
    phgInfo(1, "Dumping HYPRE ParCSRMatrix to file \"%s\"\n", fn);
    HYPRE_ParCSRMatrixPrint(mat, fn);
}

#if USE_HYPRE && HYPRE_VERSION_MAJOR >= 2

static int
send_data_to_hypre_ijmatrix(MAT *mat, HYPRE_IJMatrix A, BOOLEAN destroy)
{
    MAP *map = mat->rmap;
    INT i;
    int nprocs = map->nprocs;
    MAT_ROW *row;
    HYPRE_Int j, k, n0, r, n, *pcols = NULL, size = 0;

    phgMatAssemble(mat);

    assert(mat->type == PHG_UNPACKED);
    n0 = mat->rmap->partition[map->rank];
    if (nprocs > 1) {
	row = mat->rows;
	for (i = 0; i < mat->rmap->nlocal; i++, row++) {
	    for (k = 0; k < row->ncols; k++)
		row->cols[k] = ((r = row->cols[k]) < mat->cmap->nlocal ?
				r + n0 : mat->O2Gmap[r - mat->cmap->nlocal]);
	}
    }
    
    j = n0;
    row = mat->rows;
#if FT_PHG == FT_DOUBLE
    for (i = 0; i < mat->rmap->nlocal; i++, row++, j++) {
	if (row->ncols <= 0) {
#if 1
	    row->ncols = 1;
	    phgFree(row->cols);
	    row->cols = phgAlloc(sizeof(*row->cols));
	    row->cols[0] = j;
	    phgFree(row->data);
	    row->data = phgAlloc(sizeof(*row->data));
	    row->data[0] = 0.;
#else	/* 0|1 */
	    phgError(1, "%s: matrix row %"dFMT" is empty!\n", __func__, i);
#endif	/* 0|1 */
	}
	n = row->ncols;
	if (sizeof(INT) != sizeof(HYPRE_Int)) {
	    if (size < n) {
		phgFree(pcols);
		pcols = phgAlloc((size = n) * sizeof(*pcols));
	    }
	    for (k = 0; k < n; k++)
		pcols[k] = (HYPRE_Int)row->cols[k];
	    HYPRE_IJMatrixSetValues(A, 1, &n, &j, pcols, row->data);
	}
	else {
	    HYPRE_IJMatrixSetValues(A, 1, &n, &j, (HYPRE_Int *)row->cols,
					row->data);
	}
	if (destroy) {
	    phgFree(row->cols);
	    phgFree(row->data);
	    row->cols = NULL;
	    row->data = NULL;
	    row->ncols = row->alloc = 0;
	}
    } /* end-i-loop */
    phgFree(pcols);
#else	/* FT_PHG == FT_DOUBLE */
    {
	double *data = NULL;

	for (i = 0; i < mat->rmap->nlocal; i++, row++, j++) {
	    if (row->ncols <= 0) {
#if 1
		row->ncols = 1;
		phgFree(row->cols);
		row->cols = phgAlloc(sizeof(*row->cols));
		row->cols[0] = j;
		phgFree(row->data);
		row->data = phgAlloc(sizeof(*row->data));
		row->data[0] = 0.;
#else	/* 0|1 */
		phgError(1, "%s: matrix row %"dFMT" is empty!\n", __func__, i);
#endif	/* 0|1 */
	    }
	    n = row->ncols;
	    if (size < n) {
		phgFree(data);
		data = phgAlloc((size = n) * sizeof(*data));
		if (sizeof(INT) != sizeof(HYPRE_Int)) {
		    phgFree(pcols);
		    pcols = phgAlloc(size * sizeof(*pcols));
		}
	    }
	    for (k = 0; k < n; k++)
		data[k] = (double)row->data[k];
	    if (sizeof(INT) != sizeof(HYPRE_Int)) {
		for (k = 0; k < n; k++)
		    pcols[k] = (double)row->cols[k];
		HYPRE_IJMatrixSetValues(A, 1, &n, &j, pcols, data);
	    }
	    else {
		HYPRE_IJMatrixSetValues(A, 1, &n, &j, (void *)row->cols, data);
	    }
	    if (destroy) {
		phgFree(row->data);
		row->data = NULL;
		phgFree(row->cols);
		row->cols = NULL;
		row->ncols = row->alloc = 0;
	    }
	} /* end-i-loop */
	phgFree(data);
	phgFree(pcols);
    }
#endif	/* FT_PHG == FT_DOUBLE */
    HYPRE_IJMatrixAssemble(A);

    if (destroy)
	mat->type = PHG_DESTROYED;

    return 0;
}

static void
ams_build_G_xyz(SOLVER *solver, HYPRE_IJMatrix *G, HYPRE_IJVector *x,
		HYPRE_IJVector *y, HYPRE_IJVector *z)
/* builds the gradients matrix G and coordinates vectors x, y, z.
 * (only valid for lowest order elements) */
{
    MAP *map = solver->mat->rmap;
    DOF *u = map->dofs[0], *dof;
    GRID *g = u->g;
    ELEMENT *e;
    SOLVER *solver_vert;
    HYPRE_Int ilower_edge = map->partition[map->rank];
    HYPRE_Int iupper_edge = map->partition[map->rank + 1] - 1;
    HYPRE_Int ilower_vert,iupper_vert;
    int i;
    INT v0, v1;
    HYPRE_Int row, cols[2], ncols = 2;
    double values[2];
    double xyz[Dim];
    COORD *c;

    assert(u->dim == 1 && DofTypeDim(u) == 3 && u->type->order <= 1);
    assert(u->type == DOF_ND1 || u->type == DOF_HC0 || G == NULL);

    dof = phgDofNew(g, DOF_P1, 1, "test", DofNoData);
    solver_vert = phgSolverCreate(SOLVER_DEFAULT, dof, NULL);
    ilower_vert = solver_vert->mat->rmap->partition[map->rank];
    iupper_vert = solver_vert->mat->rmap->partition[map->rank+1] - 1;

    if (G != NULL) {
	if (*G != NULL)
	    HYPRE_IJMatrixDestroy(*G);
	HYPRE_IJMatrixCreate(g_comm, ilower_edge, iupper_edge, 
				      ilower_vert, iupper_vert, G);
	HYPRE_IJMatrixSetObjectType(*G, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(*G);
    }

    if (*x != NULL)
	HYPRE_IJVectorDestroy(*x);
    HYPRE_IJVectorCreate(g_comm, ilower_vert, iupper_vert, x);
    HYPRE_IJVectorSetObjectType(*x, HYPRE_PARCSR);
    HYPRE_IJVectorSetMaxOffProcElmts(*x, 0);
    HYPRE_IJVectorInitialize(*x);

    if (*y != NULL)
	HYPRE_IJVectorDestroy(*y);
    HYPRE_IJVectorCreate(g_comm, ilower_vert, iupper_vert, y);
    HYPRE_IJVectorSetObjectType(*y, HYPRE_PARCSR);
    HYPRE_IJVectorSetMaxOffProcElmts(*y, 0);
    HYPRE_IJVectorInitialize(*y);

    if (*z != NULL)
	HYPRE_IJVectorDestroy(*z);
    HYPRE_IJVectorCreate(g_comm, ilower_vert, iupper_vert, z);
    HYPRE_IJVectorSetObjectType(*z, HYPRE_PARCSR);
    HYPRE_IJVectorSetMaxOffProcElmts(*z, 0);
    HYPRE_IJVectorInitialize(*z);

    values[0]=1.0, values[1]=-1.0;
    ForAllElements(g, e){
	if (G != NULL) {
	    for (i = 0; i < NEdge; i++) {
		row = phgSolverMapE2G(solver, 0, e, i);
		if (row < ilower_edge || row > iupper_edge)
		    continue;
		GetEdgeVertices(e, i, v0, v1);
		cols[0] = phgSolverMapE2G(solver_vert, 0, e, v0);
		cols[1] = phgSolverMapE2G(solver_vert, 0, e, v1);
		HYPRE_IJMatrixSetValues(*G, 1, &ncols, &row, cols, values);
	    }
	}

	for (i = 0; i < NVert; i++) {
	    row = phgSolverMapE2G(solver_vert, 0, e, i);
	    if (row >= ilower_vert && row <= iupper_vert) {
		xyz[0] = (*(c = g->verts + e->verts[i]))[0];
		xyz[1] = (*c)[1];
		xyz[2] = (*c)[2];
		HYPRE_IJVectorSetValues(*x, 1, &row, xyz);
		HYPRE_IJVectorSetValues(*y, 1, &row, xyz + 1);
		HYPRE_IJVectorSetValues(*z, 1, &row, xyz + 2);
	     }
	 }
    }	/* ForAllElements */

    phgDofFree(&dof);
    phgSolverDestroy(&solver_vert);
    if (G != NULL)
	HYPRE_IJMatrixAssemble(*G);
    HYPRE_IJVectorAssemble(*x);
    HYPRE_IJVectorAssemble(*y);
    HYPRE_IJVectorAssemble(*z);
    
    return;
}

static DOF *dof_H1;

static void
ams_func_G(DOF *u, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* dummy function for computing G.
 * Note: the results are stored as values[u->dim][Dim] */
{
    int i;
    const FLOAT *p, (*J)[Dim + 1];

    Unused(bno);
    CheckD_order(u);
    /* TODO: cache grad_lambda of the bases of dof_H1 */
    p = dof_H1->type->BasGrads(dof_H1, e, 0, -1, lambda);
    J = (void *)phgGeomGetJacobian(u->g, e);
    for (i = 0; i < u->dim; i++, values += Dim, p += Dim + 1) {
	values[0] = p[0]*J[0][0] + p[1]*J[1][0] + p[2]*J[2][0] + p[3]*J[3][0];
	values[1] = p[0]*J[0][1] + p[1]*J[1][1] + p[2]*J[2][1] + p[3]*J[3][1];
	values[2] = p[0]*J[0][2] + p[1]*J[1][2] + p[2]*J[2][2] + p[3]*J[3][2];
    }
}

#if HAVE_HYPRE_AMS_SET_PI
static void
ams_func_Pi(DOF *u, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* dummy function for computing Pi.
 * Note: the results are stored as values[u->dim][Dim] */
{
    int i;
    const FLOAT *p;

    Unused(bno);
    CheckD_order(u);
    /* TODO: cache bases of dof_H1 */
    p = dof_H1->type->BasFuncs(dof_H1, e, 0, -1, lambda);
    for (i = 0; i < u->dim / Dim; i++, values += 3 * Dim, p++) {
	values[0] = *p;
	values[1] = 0.;
	values[2] = 0.;

	values[3] = 0.;
	values[4] = *p;
	values[5] = 0.;

	values[6] = 0.;
	values[7] = 0.;
	values[8] = *p;
    }
}
#endif	/* HAVE_HYPRE_AMS_SET_PI */

static void
ams_build_G_Pi(SOLVER *solver, HYPRE_IJMatrix *G, HYPRE_IJMatrix *Pi)
/* build the gradient matrix G whose elements are determined by:
 *	\grad\phi_j = \sum g_{kj} \psi_k
 * where \phi_j is the j-th basis function of the H1 space,
 * and \psi_k is the k-th basis function of the Hcurl space */
{
    MAP *map = solver->mat->rmap;
    DOF *u = map->dofs[0];
    GRID *g = u->g;
    ELEMENT *e;
    SOLVER *solver_H1;
    DOF *dof_Hcurl;
#if HAVE_HYPRE_AMS_SET_PI
    DOF *dof_Hcurl3;
#endif	/* HAVE_HYPRE_AMS_SET_PI */
    DOF_TYPE *h1_type = NULL;
    HYPRE_Int ilower_Hcurl = map->partition[map->rank];
    HYPRE_Int iupper_Hcurl = map->partition[map->rank + 1] - 1;
    HYPRE_Int ilower_H1, iupper_H1;
    HYPRE_Int i, j, N_Hcurl, N_H1, row, ncols, *cols0, *cols;
    FLOAT *p, *floats, *pdata[15];
    double *doubles;

    assert(map->ndof == 1 && u->type->dim == 3);

    h1_type = DOF_Pn[u->type->order];

    N_H1 = h1_type->nbas;
    N_Hcurl = u->type->nbas;

    dof_H1 = phgDofNew(g, h1_type, 1, "tmp", DofNoData);
    solver_H1 = phgSolverCreate(SOLVER_DEFAULT, dof_H1, NULL);
    ilower_H1 = solver_H1->mat->rmap->partition[map->rank];
    iupper_H1 = solver_H1->mat->rmap->partition[map->rank + 1] - 1;

    dof_Hcurl = phgDofNew(g, u->type, N_H1, "tmp", DofNoData);
    if (*G != NULL)
	HYPRE_IJMatrixDestroy(*G);
    HYPRE_IJMatrixCreate(g_comm, ilower_Hcurl, iupper_Hcurl, 
	     ilower_H1, iupper_H1, G);
    HYPRE_IJMatrixSetObjectType(*G, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(*G);

#if HAVE_HYPRE_AMS_SET_PI
    dof_Hcurl3 = phgDofNew(g, u->type, Dim * N_H1, "tmp3", DofNoData);
    if (*Pi != NULL)
	HYPRE_IJMatrixDestroy(*Pi);
    HYPRE_IJMatrixCreate(g_comm, ilower_Hcurl, iupper_Hcurl, 
	     Dim * ilower_H1, Dim * (iupper_H1 + 1) - 1, Pi);
    HYPRE_IJMatrixSetObjectType(*Pi, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(*Pi);
    i = Dim;
#else	/* HAVE_HYPRE_AMS_SET_PI */
    i = 1;
#endif	/* HAVE_HYPRE_AMS_SET_PI */

    floats = phgAlloc(i * N_Hcurl * N_H1 * sizeof(*floats));
    doubles = phgAlloc(i * N_H1 * sizeof(*doubles));
    cols0 = phgAlloc((1 + i) * N_H1 * sizeof(*cols0));
    cols = cols0 + N_H1;

    ForAllElements(g, e) {
	for (j = 0; j < N_H1; j++)
	    cols0[j] = phgSolverMapE2G(solver_H1, 0, e, j);
	/* The gradient matrix */
	p = floats;
	j = dof_Hcurl->dim * dof_Hcurl->type->np_edge;
	for (i = 0; i < NEdge; i++) {
	    u->type->InitFunc(dof_Hcurl, e, EDGE, i, NULL, ams_func_G,
				NULL, p, pdata);
	    pdata[4 + i] = p;
	    p += j;
	}

	j = dof_Hcurl->dim * dof_Hcurl->type->np_face;
	if (j > 0) {
	    for (i = 0; i < NFace; i++) {
		u->type->InitFunc(dof_Hcurl, e, FACE, i, NULL, ams_func_G,
					NULL, p, pdata);
		pdata[10 + i] = p;
		p += j;
	    }
	}

	if (dof_Hcurl->type->np_elem > 0)
	    u->type->InitFunc(dof_Hcurl, e, VOLUME, 0, NULL, ams_func_G,
				NULL, p, pdata);

	p = floats;
	for (i = 0; i < N_Hcurl; i++) {
	    row = phgSolverMapE2G(solver, 0, e, i);
	    for (j = 0, ncols = 0; j < N_H1; j++, p++) {
		if (*p == 0.0)
		    continue;
		cols[ncols] = cols0[j];
		doubles[ncols++] = *p;
	    }
	    HYPRE_IJMatrixSetValues(*G, 1, &ncols, &row, cols, doubles);
	}

#if HAVE_HYPRE_AMS_SET_PI
	/* The (H1)^3 matrix */
	p = floats;
	j = dof_Hcurl3->dim * dof_Hcurl3->type->np_edge;
	for (i = 0; i < NEdge; i++) {
	    u->type->InitFunc(dof_Hcurl3, e, EDGE, i, NULL, ams_func_Pi,
				NULL, p, pdata);
	    pdata[4 + i] = p;
	    p += j;
	}

	j = dof_Hcurl3->dim * dof_Hcurl3->type->np_face;
	if (j > 0) {
	    for (i = 0; i < NFace; i++) {
		u->type->InitFunc(dof_Hcurl3, e, FACE, i, NULL, ams_func_Pi,
					NULL, p, pdata);
		pdata[10 + i] = p;
		p += j;
	    }
	}

	if (dof_Hcurl3->type->np_elem > 0)
	    u->type->InitFunc(dof_Hcurl3, e, VOLUME, 0, NULL, ams_func_Pi,
				NULL, p, pdata);

	p = floats;
	for (i = 0; i < N_Hcurl; i++) {
	    row = phgSolverMapE2G(solver, 0, e, i);
	    for (j = 0, ncols = 0; j < Dim * N_H1; j++, p++) {
		if (*p == 0.0)
		    continue;
		cols[ncols] = Dim * cols0[j / Dim] + (j % Dim);
		doubles[ncols++] = *p;
	    }
	    HYPRE_IJMatrixSetValues(*Pi, 1, &ncols, &row, cols, doubles);
	}
#endif	/* HAVE_HYPRE_AMS_SET_PI */
    }	/* ForAllElements */

    phgFree(cols0);
    phgFree(floats);
    phgFree(doubles);

    phgSolverDestroy(&solver_H1);
    phgDofFree(&dof_H1);
    phgDofFree(&dof_Hcurl);
    HYPRE_IJMatrixAssemble(*G);
#if HAVE_HYPRE_AMS_SET_PI
    phgDofFree(&dof_Hcurl3);
    HYPRE_IJMatrixAssemble(*Pi);
#endif
    
    return;
}

static void
ams_destroy(SOLVER *solver)
{
    if (_ams == NULL)
	return;

    if (_ams->G != NULL) {
	HYPRE_IJMatrixDestroy(_ams->G);
	_ams->G = NULL;
    }
    if (_ams->Pi != NULL) {
	HYPRE_IJMatrixDestroy(_ams->Pi);
	_ams->Pi = NULL;
    }
    if (_ams->x != NULL) {
	HYPRE_IJVectorDestroy(_ams->x);
	_ams->x = NULL;
    }
    if (_ams->y != NULL) {
	HYPRE_IJVectorDestroy(_ams->y);
	_ams->y = NULL;
    }
    if (_ams->z != NULL) {
	HYPRE_IJVectorDestroy(_ams->z);
	_ams->z = NULL;
    }

    if (_ams->A_alpha != NULL) {
	HYPRE_IJMatrixDestroy(_ams->A_alpha);
	_ams->A_alpha = NULL;
    }
    if (_ams->A_beta != NULL) {
	HYPRE_IJMatrixDestroy(_ams->A_beta);
	_ams->A_beta = NULL;
    }

    phgFree(_ams);
    _ams = NULL;
    phgInfo(1, "AMS solver destroyed.\n");

    return;
}

static void
ams_init(SOLVER *solver)
{
    if (!solver->oem_created) {
	solver->oem_created = TRUE;
	solver->oem_solver->Create(solver);
    }

    if (_ams != NULL)
	return;

    phgInfo(1, "Initializing AMS solver.\n");
    _ams = phgAlloc(sizeof(*_ams));
    _ams->G = NULL;
    _ams->Pi = NULL;
    _ams->A_alpha = NULL;
    _ams->A_beta = NULL;
    _ams->x = NULL;
    _ams->y = NULL;
    _ams->z = NULL;

    return;
}

void
phgSolverHypreAMSSetPoisson_(SOLVER *solver, MAT *alpha, MAT *beta,
			     BOOLEAN destroy)
{
    MAP *map = solver->mat->rmap;
    HYPRE_Int ilower_vert, iupper_vert;
    HYPRE_IJMatrix A_alpha=NULL;
    HYPRE_IJMatrix A_beta=NULL;
    
    ams_init(solver);

    if (alpha == NULL)
	return;

    if (_ams->A_alpha != NULL)
	HYPRE_IJMatrixDestroy(_ams->A_alpha);
    if (_ams->A_beta != NULL)
	HYPRE_IJMatrixDestroy(_ams->A_beta);

    ilower_vert = alpha->rmap->partition[map->rank];
    iupper_vert = alpha->rmap->partition[map->rank + 1] - 1;
    HYPRE_IJMatrixCreate(g_comm, ilower_vert, iupper_vert, ilower_vert,
				     iupper_vert, &A_alpha);
    HYPRE_IJMatrixSetObjectType(A_alpha, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A_alpha);
    send_data_to_hypre_ijmatrix(alpha, A_alpha, destroy);
    _ams->A_alpha = A_alpha;

    if (beta != NULL) {
	ilower_vert = beta->rmap->partition[map->rank];
	iupper_vert = beta->rmap->partition[map->rank + 1] - 1;
	HYPRE_IJMatrixCreate(g_comm, ilower_vert, iupper_vert, ilower_vert,
						   iupper_vert, &A_beta);
	HYPRE_IJMatrixSetObjectType(A_beta, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A_beta);
	send_data_to_hypre_ijmatrix(beta, A_beta, destroy);
	_ams->A_beta = A_beta;
    }
    else {
	_ams->A_beta = NULL;
    }

    return;
}

/* user-callable functions */
void
phgSolverHypreAMSSetPoisson(SOLVER *solver, DOF *alpha, DOF *beta)
{
    DOF *tmp, *u;
    MAP *map;
    MAT *ma, *mb;
    int i, j, n, N;
    GRID *g;
    ELEMENT *e;
    FLOAT *A, *B = NULL, *A1, *B1 = NULL, *A0, d;
    INT *I, *J;
    BTYPE *b, DB_mask;

    if (solver->oem_solver != SOLVER_HYPRE)
	return;

    u = solver->mat->rmap->dofs[0];
    g = u->g;

    /* TODO: check DOF_TYPE of solver->dofs[] for consistency */
    assert(solver->mat->rmap->ndof == 1 && DofTypeDim(u) == 3 &&
	   u->dim == 1);
    assert(beta == NULL || beta->dim == 1);
    assert(solver->mat->rmap->dofs != NULL);

    DB_mask = solver->mat->rmap->dofs[0]->DB_mask;

    tmp = phgDofNew(g, DOF_Pn[u->type->order], 1, "H1", DofNoData);
    map = phgMapCreate(tmp, NULL);
    ma = phgMapCreateMat(map, map);
    if (beta != NULL && beta->type == DOF_CONSTANT && *beta->data == 0.) {
	mb = NULL;
    }
    else {
	mb = phgMapCreateMat(map, map);
    }
    phgMapDestroy(&map);

    N = tmp->type->nbas;	/* number of basis functions in an element */
    A = phgAlloc(N * N * sizeof(*A));
    A0 = phgAlloc(N * N * sizeof(*A0));
    A1 = phgAlloc(N * sizeof(*A1));
    if (mb != NULL) {
	B = phgAlloc(N * N * sizeof(*B));
	B1 = phgAlloc(N * sizeof(*B1));
    }
    I = phgAlloc(N * sizeof(*I));
    J = phgAlloc(N * sizeof(*J));
    b = phgAlloc(N * sizeof(*b));

    /* Note: if beta is constant then the element mass matrices only vary
       by element volumes */
    if (beta == NULL || beta->type == DOF_CONSTANT) {
	e = g->roots;
	d = (beta == NULL ? 1.0 : *beta->data) / phgGeomGetVolume(g, e);
	for (i = 0; i < N; i++)
	    for (j = 0; j <= i; j++)
		A0[i * N + j] = A0[j * N + i] =
		    d * phgQuadBasDotBas(e, tmp, j, tmp, i, -1);
    }

    ForAllElements(g, e) {
	if (beta != NULL && beta->type != DOF_CONSTANT) {
	    /* compute \int beta * phi_i * phi_j */
	    for (i = 0; i < N; i++)
		for (j = 0; j <= i; j++)
		    A[i * N + j] = phgQuadBasABas(e, tmp, j, beta, tmp, i, -1);
	}
	else {
	    d = phgGeomGetVolume(g, e);
	    for (i = 0; i < N; i++)
		for (j = 0; j <= i; j++)
		    A[i * N + j] = A0[i * N + j] * d;
	}
	for (i = 0; i < N; i++) {
	    I[i] = phgMapE2L(ma->rmap, 0, e, i);
	    b[i] = phgDofGetElementBoundaryType(tmp, e, i * tmp->dim);
	    for (j = 0; j <= i; j++) {
		/* FIXME: use phgQuadGradBasDotGradBas if both alpha and beta
		 * are constant */
		A[j * N + i] = (A[i * N + j] +=
			phgQuadGradBasAGradBas(e, tmp, j, alpha, tmp, i, -1));
		if (mb != NULL)
		    B[j * N + i] = B[i * N + j] =
			phgQuadGradBasAGradBas(e, tmp, j, beta, tmp, i, -1);
	    }
	}
	/* Stiffness matrix */
	for (i = 0; i < N; i++) {
	    if (b[i] & DB_mask) {
#if 0
		phgMatAddEntry(ma, I[i], I[i], 1.0);
		phgMatAddEntry(mb, I[i], I[i], 1.0);
#else
		phgMatAddEntry(ma, I[i], I[i], A[i * N + i]);
		if (mb != NULL)
		    phgMatAddEntry(mb, I[i], I[i], B[i * N + i]);
#endif
	    }
	    else {
		for (j = 0, n = 0; j < N; j++) {
		    if (b[j] & DB_mask)
			continue;
		    J[n] = I[j];
		    A1[n] = A[i * N + j];
		    if (mb != NULL)
			B1[n] = B[i * N + j];
		    n++;
		}
		phgMatAddEntries(ma, 1, I + i, n, J, A1);
		if (mb != NULL)
		    phgMatAddEntries(mb, 1, I + i, n, J, B1);
	    }
	}
    }	/* end-ForAllElement */

    phgFree(A);
    phgFree(A0);
    phgFree(B);
    phgFree(A1);
    phgFree(B1);
    phgFree(I);
    phgFree(J);
    phgFree(b);

    phgSolverHypreAMSSetPoisson_(solver, ma, mb, TRUE);
    phgMatDestroy(&ma);
    phgMatDestroy(&mb);
    phgDofFree(&tmp);
}

void
phgSolverHypreAMSSetConstantPoisson(SOLVER *solver, FLOAT alpha, FLOAT beta)
{
    GRID *g;
    DOF *da, *db;

    if (solver->oem_solver != SOLVER_HYPRE)
	return;

    g = solver->mat->rmap->dofs[0]->g;

    if (alpha == 1.0) {
	da = NULL;
    }
    else {
	da = phgDofNew(g, DOF_CONSTANT, 1, "alpha", DofNoAction);
	phgDofSetDataByValue(da, alpha);
    }

    if (beta == 1.0) {
	db = NULL;
    }
    else {
	db = phgDofNew(g, DOF_CONSTANT, 1, "beta", DofNoAction);
	phgDofSetDataByValue(db, beta);
    }

    phgSolverHypreAMSSetPoisson(solver, da, db);

    phgDofFree(&da);
    phgDofFree(&db);
}

#endif	/* USE_HYPRE && HYPRE_VERSION_MAJOR >= 2 */

static void
setup_hypre_solver(SOLVER *solver, SOLVER_ID id, HYPRE_Solver *hsolver,
		   FLOAT rtol, FLOAT atol, INT maxit,
		   SolveFcn *setup, SolveFcn *solve, DestroyFcn *destroy,
		   SetTolFcn *set_tol, SetTolFcn *set_atol,
		   SetMaxIterFcn *set_maxiter,
		   SetPCFcn *set_pc, GetItFcn *get_nits,
		   GetResidFcn *get_resid, int level0)
{
    MAP *map = solver->mat->rmap;
#if HYPRE_VERSION_MAJOR >= 2
    HYPRE_ParCSRMatrix  parG, parPi, parA_alpha, parA_beta;
    HYPRE_ParVector	par_x, par_y, par_z;
#endif /* HYPRE_VERSION_MAJOR >= 2 */

    switch (id) {
	case PCG:	/* Preconditioned CG solver */
	    HYPRE_ParCSRPCGCreate(g_comm, hsolver);
	    if (set_maxiter == NULL)
		HYPRE_PCGSetMaxIter(*hsolver, maxit);	/* max iterations */
	    else
		*set_maxiter = HYPRE_PCGSetMaxIter;
	    if (set_tol == NULL)
		HYPRE_PCGSetTol(*hsolver, rtol);	/* conv. tolerance */
	    else
		*set_tol = HYPRE_PCGSetTol;
#if HYPRE_VERSION_MAJOR >= 2 && HYPRE_VERSION_MINOR >= 4
	    if (set_atol == NULL)
		HYPRE_PCGSetAbsoluteTol(*hsolver, atol);/* conv. tolerance */
	    else
		*set_atol = HYPRE_PCGSetAbsoluteTol;
#endif	/* HYPRE_VERSION_MAJOR >= 2 && HYPRE_VERSION_MINOR >= 4 */
	    HYPRE_PCGSetTwoNorm(*hsolver, _prm.two_norm); /* whether 2-norm */
	    HYPRE_PCGSetPrintLevel(*hsolver, solver->monitor + level0);
	    HYPRE_PCGSetLogging(*hsolver, 0);
	    *setup = HYPRE_ParCSRPCGSetup;
	    *solve = HYPRE_ParCSRPCGSolve;
	    *destroy = HYPRE_ParCSRPCGDestroy;
	    if (set_pc != NULL)
		*set_pc = HYPRE_PCGSetPrecond;
	    if (get_nits != NULL)
		*get_nits = HYPRE_PCGGetNumIterations;
	    if (get_resid != NULL)
		*get_resid = HYPRE_PCGGetFinalRelativeResidualNorm;
	    break;
	case GMRES:	/* GMRES */
	    HYPRE_ParCSRGMRESCreate(g_comm, hsolver);
	    if (set_maxiter == NULL)
		HYPRE_ParCSRGMRESSetMaxIter(*hsolver, maxit);
	    else
		*set_maxiter = HYPRE_ParCSRGMRESSetMaxIter;
	    if (set_tol == NULL)
		HYPRE_ParCSRGMRESSetTol(*hsolver, rtol);
	    else
		*set_tol = HYPRE_ParCSRGMRESSetTol;
#if HYPRE_VERSION_MAJOR >= 2 && HYPRE_VERSION_MINOR >= 4
	    if (set_atol == NULL)
		HYPRE_ParCSRGMRESSetAbsoluteTol(*hsolver, atol);
	    else
		*set_atol = HYPRE_ParCSRGMRESSetAbsoluteTol;
#endif	/* HYPRE_VERSION_MAJOR >= 2 && HYPRE_VERSION_MINOR >= 4 */
	    HYPRE_ParCSRGMRESSetPrintLevel(*hsolver, solver->monitor + level0);
	    HYPRE_ParCSRGMRESSetLogging(*hsolver, 0);
	    if (_prm.gmres_kdim > 0)
		HYPRE_ParCSRGMRESSetKDim (*hsolver, _prm.gmres_kdim);
	    *setup = HYPRE_ParCSRGMRESSetup;
	    *solve = HYPRE_ParCSRGMRESSolve;
	    *destroy = HYPRE_ParCSRGMRESDestroy;
	    if (set_pc != NULL)
		*set_pc = HYPRE_GMRESSetPrecond/*HYPRE_ParCSRGMRESSetPrecond*/;
	    if (get_nits != NULL)
		*get_nits = HYPRE_ParCSRGMRESGetNumIterations;
	    if (get_resid != NULL)
		*get_resid = HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm;
	    break;
#if HYPRE_VERSION_MAJOR >= 2 && HYPRE_VERSION_MINOR >= 4
	case LGMRES:	/* LGMRES */
	    HYPRE_ParCSRLGMRESCreate(g_comm, hsolver);
	    if (set_maxiter == NULL)
		HYPRE_ParCSRLGMRESSetMaxIter(*hsolver, maxit);
	    else
		*set_maxiter = HYPRE_ParCSRLGMRESSetMaxIter;
	    if (set_tol == NULL)
		HYPRE_ParCSRLGMRESSetTol(*hsolver, rtol);
	    else
		*set_tol = HYPRE_ParCSRLGMRESSetTol;
	    if (set_atol == NULL)
		HYPRE_ParCSRLGMRESSetAbsoluteTol(*hsolver, atol);
	    else
		*set_atol = HYPRE_ParCSRLGMRESSetAbsoluteTol;
	    HYPRE_ParCSRLGMRESSetPrintLevel(*hsolver, solver->monitor + level0);
	    HYPRE_ParCSRLGMRESSetLogging(*hsolver, 0);
	    if (_prm.gmres_kdim > 0)
		HYPRE_ParCSRLGMRESSetKDim (*hsolver, _prm.gmres_kdim);
	    *setup = HYPRE_ParCSRLGMRESSetup;
	    *solve = HYPRE_ParCSRLGMRESSolve;
	    *destroy = HYPRE_ParCSRLGMRESDestroy;
	    if (set_pc != NULL)
		*set_pc = HYPRE_LGMRESSetPrecond/*HYPRE_ParCSRLGMRESSetPrecond*/;
	    if (get_nits != NULL)
		*get_nits = HYPRE_ParCSRLGMRESGetNumIterations;
	    if (get_resid != NULL)
		*get_resid = HYPRE_ParCSRLGMRESGetFinalRelativeResidualNorm;
	    break;
	case FGMRES:	/* FlexGMRES */
	    HYPRE_ParCSRFlexGMRESCreate(g_comm, hsolver);
	    if (set_maxiter == NULL)
		HYPRE_ParCSRFlexGMRESSetMaxIter(*hsolver, maxit);
	    else
		*set_maxiter = HYPRE_ParCSRFlexGMRESSetMaxIter;
	    if (set_tol == NULL)
		HYPRE_ParCSRFlexGMRESSetTol(*hsolver, rtol);
	    else
		*set_tol = HYPRE_ParCSRFlexGMRESSetTol;
	    if (set_atol == NULL)
		HYPRE_ParCSRFlexGMRESSetAbsoluteTol(*hsolver, atol);
	    else
		*set_atol = HYPRE_ParCSRFlexGMRESSetAbsoluteTol;
	    HYPRE_ParCSRFlexGMRESSetPrintLevel(*hsolver, solver->monitor + level0);
	    HYPRE_ParCSRFlexGMRESSetLogging(*hsolver, 0);
	    if (_prm.gmres_kdim > 0)
		HYPRE_ParCSRFlexGMRESSetKDim (*hsolver, _prm.gmres_kdim);
	    *setup = HYPRE_ParCSRFlexGMRESSetup;
	    *solve = HYPRE_ParCSRFlexGMRESSolve;
	    *destroy = HYPRE_ParCSRFlexGMRESDestroy;
	    if (set_pc != NULL)
		*set_pc = HYPRE_FlexGMRESSetPrecond/*HYPRE_ParCSRFlexGMRESSetPrecond*/;
	    if (get_nits != NULL)
		*get_nits = HYPRE_ParCSRFlexGMRESGetNumIterations;
	    if (get_resid != NULL)
		*get_resid = HYPRE_ParCSRFlexGMRESGetFinalRelativeResidualNorm;
	    break;
#endif	/* HYPRE_VERSION_MAJOR >= 2 && HYPRE_VERSION_MINOR >= 4 */
	case BiCGSTAB:	/* BiCGSTAB */
	    HYPRE_ParCSRBiCGSTABCreate(g_comm, hsolver);
	    if (set_maxiter == NULL)
		HYPRE_ParCSRBiCGSTABSetMaxIter(*hsolver, maxit);
	    else
		*set_maxiter = HYPRE_ParCSRBiCGSTABSetMaxIter;
	    if (set_tol == NULL)
		HYPRE_ParCSRBiCGSTABSetTol(*hsolver, rtol);
	    else
		*set_tol = HYPRE_ParCSRBiCGSTABSetTol;
#if HYPRE_VERSION_MAJOR >= 2 && HYPRE_VERSION_MINOR >= 4
	    if (set_atol == NULL)
		HYPRE_ParCSRBiCGSTABSetAbsoluteTol(*hsolver, atol);
	    else
		*set_atol = HYPRE_ParCSRBiCGSTABSetAbsoluteTol;
#endif	/* HYPRE_VERSION_MAJOR >= 2 && HYPRE_VERSION_MINOR >= 4 */
	    HYPRE_ParCSRBiCGSTABSetPrintLevel(*hsolver, solver->monitor + level0);
	    HYPRE_ParCSRBiCGSTABSetLogging(*hsolver, 0);
	    *setup = HYPRE_ParCSRBiCGSTABSetup;
	    *solve = HYPRE_ParCSRBiCGSTABSolve;
	    *destroy = HYPRE_ParCSRBiCGSTABDestroy;
	    if (set_pc != NULL)
		*set_pc = HYPRE_BiCGSTABSetPrecond;
	    if (get_nits != NULL)
		*get_nits = HYPRE_ParCSRBiCGSTABGetNumIterations;
	    if (get_resid != NULL)
		*get_resid = HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm;
	    break;
	case BoomerAMG:	/* BoomerAMG */
	    HYPRE_BoomerAMGCreate(hsolver);
	    HYPRE_BoomerAMGSetPrintLevel(*hsolver,
		solver->monitor+level0 > 0 ? 2*(solver->monitor+level0-1) : 0);
	    HYPRE_BoomerAMGSetStrongThreshold(*hsolver, _prm.amg_strength);
	    HYPRE_BoomerAMGSetMaxRowSum(*hsolver, _prm.amg_max_row_sum);
	    HYPRE_BoomerAMGSetCoarsenType(*hsolver,	/* coarsening type */
			amg_coarsen_types[_prm.amg_coarsen_type]);
	    HYPRE_BoomerAMGSetRelaxType(*hsolver,
			amg_relax_types[_prm.amg_relax_type]);
	    HYPRE_BoomerAMGSetCycleRelaxType(*hsolver,
			amg_relax_types[_prm.amg_coarsest_relax_type], 3);
	    HYPRE_BoomerAMGSetNumSweeps(*hsolver, 1);	/* swps on each level */
	    if (_prm.amg_num_funcs > 0)			/* size of the system */
		HYPRE_BoomerAMGSetNumFunctions(*hsolver, _prm.amg_num_funcs);
	    if (_prm.amg_max_levels > 0)		/* maxi nbr of levels */
		HYPRE_BoomerAMGSetMaxLevels(*hsolver, _prm.amg_max_levels);
	    if (set_maxiter == NULL)
		HYPRE_BoomerAMGSetMaxIter(*hsolver, maxit);
	    else
		*set_maxiter = HYPRE_BoomerAMGSetMaxIter;
	    if (set_tol == NULL)
		HYPRE_BoomerAMGSetTol(*hsolver, rtol);
	    else
		*set_tol = HYPRE_BoomerAMGSetTol;
#if 0	/* Note: HYPRE_BoomerAMGSetAbsoluteTol() not available */
	    if (set_atol == NULL)
		HYPRE_BoomerAMGSetAbsoluteTol(*hsolver, atol);
	    else
		*set_atol = HYPRE_BoomerAMGSetAbsoluteTol;
#endif
	    *setup = HYPRE_BoomerAMGSetup;
	    *solve = HYPRE_BoomerAMGSolve;
	    *destroy = HYPRE_BoomerAMGDestroy;
	    if (set_pc != NULL)
		*set_pc = NULL;
	    if (get_nits != NULL)
		*get_nits = HYPRE_BoomerAMGGetNumIterations;
	    if (get_resid != NULL)
		*get_resid = HYPRE_BoomerAMGGetFinalRelativeResidualNorm;
	    break;
	case ParaSails:	/* ParaSails solver */
	    HYPRE_ParaSailsCreate(g_comm, hsolver);
	    HYPRE_ParaSailsSetParams(*hsolver,	_prm.sai_threshold,
						_prm.sai_max_levels);
	    HYPRE_ParaSailsSetFilter(*hsolver, _prm.sai_filter);
	    HYPRE_ParaSailsSetSym(*hsolver, _prm.sai_sym);
	    HYPRE_ParaSailsSetLogging(*hsolver, 0);
	    *setup = HYPRE_ParaSailsSetup;
	    *solve = HYPRE_ParaSailsSolve;
	    *destroy = HYPRE_ParaSailsDestroy;
	    if (set_pc != NULL)
		*set_pc = NULL;
	    if (get_nits != NULL)
		*get_nits = NULL;
	    if (get_resid != NULL)
		*get_resid = NULL;
	    break;
	case DIAG:	/* diagonal preconditioning */
	    *setup = HYPRE_ParCSRDiagScaleSetup;
	    *solve = HYPRE_ParCSRDiagScale;
	    *destroy = NULL;
	    *hsolver = NULL;
	    if (set_pc != NULL)
		*set_pc = NULL;
	    if (get_nits != NULL)
		*get_nits = NULL;
	    if (get_resid != NULL)
		*get_resid = NULL;
	    break;
	case Euclid:	/* Euclid solver (ILU) */
	    HYPRE_EuclidCreate(g_comm, hsolver);
	    /* pass options from command line */
	    HYPRE_EuclidSetParams(*hsolver, argc, argv);
	    /* alternatively, pass options from a configuration file
	     * HYPRE_EuclidSetParamsFromFile(*hsolver, "filename"); */
	    *setup = HYPRE_EuclidSetup;
	    *solve = HYPRE_EuclidSolve;
	    *destroy = HYPRE_EuclidDestroy;
	    if (set_pc != NULL)
		*set_pc = NULL;
	    if (get_nits != NULL)
		*get_nits = NULL;
	    if (get_resid != NULL)
		*get_resid = NULL;
	    break;
#if HYPRE_VERSION_MAJOR >= 2
	case AMS:
	    HYPRE_AMSCreate(hsolver);
	    HYPRE_AMSSetDimension(*hsolver, _prm.ams_dim);
	    if (set_maxiter == NULL)
		HYPRE_AMSSetMaxIter(*hsolver, maxit);
	    else
		*set_maxiter = HYPRE_AMSSetMaxIter;
	    if (set_tol == NULL)
		HYPRE_AMSSetTol(*hsolver, rtol);
	    else
		*set_tol = HYPRE_AMSSetTol;
#if 0	/* Note: HYPRE_AMSSetAbsoluteTol() not available */
	    if (set_atol == NULL)
		HYPRE_AMSSetAbsoluteTol(*hsolver, atol);
	    else
		*set_atol = HYPRE_AMSSetAbsoluteTol;
#endif
	    if(_prm.cycle_type > 7)
	    	HYPRE_AMSSetCycleType(*hsolver, _prm.cycle_type + 3);
	    else
	    	HYPRE_AMSSetCycleType(*hsolver, _prm.cycle_type + 1);
	    HYPRE_AMSSetPrintLevel(*hsolver, solver->monitor + level0);
	    ams_init(solver);
	    /* allow G and Pi to be set by user prio to calling this function */
	    if (_ams->G == NULL) {
		if (map->dofs[0]->type==DOF_ND1 ||
		    map->dofs[0]->type==DOF_HC0) {
		    /* faster code for the lowest order elements */
		    ams_build_G_xyz(solver, &_ams->G,
					    &_ams->x, &_ams->y, &_ams->z);
		}
		else {
		    ams_build_G_Pi(solver, &_ams->G, &_ams->Pi);
#if !HAVE_HYPRE_AMS_SET_PI
		    ams_build_G_xyz(solver, NULL, &_ams->x, &_ams->y, &_ams->z);
#endif	/* !HAVE_HYPRE_AMS_SET_PI */
		}
		HYPRE_IJMatrixGetObject(_ams->G, (void **)(void *)&parG);
		dump_matrix(solver, parG, "G");
		HYPRE_AMSSetDiscreteGradient(*hsolver, parG);
		if (_ams->Pi != NULL) {
		    HYPRE_IJMatrixGetObject(_ams->Pi, (void **)(void *)&parPi);
		    dump_matrix(solver, parPi, "Pi");
#if HAVE_HYPRE_AMS_SET_PI
		    HYPRE_AMSSetDiscretePi(*hsolver, parPi);
#endif	/* HAVE_HYPRE_AMS_SET_PI */
		}
		else {
		    assert(_ams->x!=NULL && _ams->y!=NULL && _ams->z!=NULL);
		    HYPRE_IJVectorGetObject(_ams->x, (void **)(void *)&par_x);
		    HYPRE_IJVectorGetObject(_ams->y, (void **)(void *)&par_y);
		    HYPRE_IJVectorGetObject(_ams->z, (void **)(void *)&par_z);
		    dump_vector(solver, par_x, "x");
		    dump_vector(solver, par_y, "y");
		    dump_vector(solver, par_z, "z");
		    HYPRE_AMSSetCoordinateVectors(*hsolver, par_x,par_y,par_z);
		}
	    }
	    else {
		if (phgVerbosity > 1)
		    phgPrintf("HYPRE AMS: using user supplied G/Pi matrices\n");
		HYPRE_IJMatrixGetObject(_ams->G, (void **)(void *)&parG);
		HYPRE_AMSSetDiscreteGradient(*hsolver, parG);
#if HAVE_HYPRE_AMS_SET_PI
		HYPRE_IJMatrixGetObject(_ams->Pi, (void **)(void *)&parPi);
		HYPRE_AMSSetDiscretePi(*hsolver, parPi);
#endif	/* HAVE_HYPRE_AMS_SET_PI */
	    }

	    /* Poisson matrices */
	    if (_ams->A_alpha != NULL) {
		HYPRE_IJMatrixGetObject(_ams->A_alpha,
						(void **)(void *)&parA_alpha);
		dump_matrix(solver, parA_alpha, "Aalpha");
		if (!_prm.ams_singular && _ams->A_beta != NULL) {
		    HYPRE_IJMatrixGetObject(_ams->A_beta,
						(void **)(void *)&parA_beta);
		    dump_matrix(solver, parA_beta, "Abeta");
		}
		else {
		    parA_beta = NULL;
		}
		HYPRE_AMSSetAlphaPoissonMatrix(*hsolver, parA_alpha);
		HYPRE_AMSSetBetaPoissonMatrix(*hsolver, parA_beta);
	    }
	    else if (_prm.ams_singular)
		HYPRE_AMSSetBetaPoissonMatrix(*hsolver, NULL);
	    HYPRE_AMSSetSmoothingOptions(*hsolver,
			_prm.rlx_type, _prm.rlx_sweeps, 1.0, 1.0);
# if (HYPRE_VERSION_MAJOR == 2 && HYPRE_VERSION_MINOR == 0)
	    HYPRE_AMSSetAlphaAMGOptions(*hsolver, _prm.ams_coarsen_type,
			_prm.ams_agg_levels, _prm.ams_rlx_type, _prm.theta);
	    HYPRE_AMSSetBetaAMGOptions(*hsolver, _prm.ams_coarsen_type,
			_prm.ams_agg_levels, _prm.ams_rlx_type, _prm.theta);
# else	/* HYPRE_VERSION_MAJOR > 2 || HYPRE_VERSION_MINOR > 0 */
	    HYPRE_AMSSetAlphaAMGOptions(*hsolver, _prm.ams_coarsen_type,
			_prm.ams_agg_levels, _prm.ams_rlx_type, _prm.theta,
			0, 0);
	    HYPRE_AMSSetBetaAMGOptions(*hsolver, _prm.ams_coarsen_type,
			_prm.ams_agg_levels, _prm.ams_rlx_type, _prm.theta,
			0, 0);
# endif	/* HYPRE_VERSION_MAJOR > 2 || HYPRE_VERSION_MINOR > 0 */
	    *setup = HYPRE_AMSSetup;
	    *solve = HYPRE_AMSSolve;
	    *destroy = HYPRE_AMSDestroy;
	    if (set_pc != NULL)
		*set_pc = NULL;
	    if (get_nits != NULL)
		*get_nits = HYPRE_AMSGetNumIterations;
	    if (get_resid != NULL)
		*get_resid = HYPRE_AMSGetFinalRelativeResidualNorm;
	    break;
#endif	/* HYPRE_VERSION_MAJOR >= 2 */
#if 0
	case MLI_PCG:		/* Another AMG preconditioner */
	    HYPRE_LSI_MLICreate(g_comm, hsolver);
	    HYPRE_LSI_MLISetParams(*hsolver, "MLI strengthThreshold 0.08");
	    *setup = HYPRE_LSI_MLISetup;
	    *solve = HYPRE_LSI_MLISolve;
	    *destroy = HYPRE_LSI_MLIDestroy;
	    if (set_pc != NULL)
		*set_pc = (void *)HYPRE_LSI_MLISetPrecond;
	    if (get_nits != NULL)
		*get_nits = HYPRE_LSI_MLIGetNumIterations;
	    if (get_resid != NULL)
		*get_resid = HYPRE_LSI_MLIGetFinalRelativeResidualNorm;
	    break;
#endif
	default:
	    *hsolver = NULL;
	    break;
    }
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

#if HYPRE_VERSION_MAJOR >= 2
    ams_destroy(solver); 
#endif	/* HYPRE_VERSION_MAJOR >= 2 */

    phgFree(_prm.dump_matvec);

    if (_t->precond != NULL && _t->precon_destroy != NULL)
	_t->precon_destroy(_t->precond);

    if (_t->hsolver != NULL)
	_t->solver_destroy(_t->hsolver);

    if (_t->A != NULL) {
	HYPRE_IJMatrixDestroy(_t->A);
    }
    if (_t->B != NULL) {
	HYPRE_IJVectorDestroy(_t->B);
    }
    if (_t->U != NULL) {
	HYPRE_IJVectorDestroy(_t->U);
    }

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    if (_t->assembled)
	return 0;
    _t->assembled = TRUE;

    HYPRE_IJMatrixAssemble(_t->A);

    /* setup solver */
    phgInfo(1, "setup HYPRE solver \"%s\".\n", solver_names[_prm.solver_id]);
    setup_hypre_solver(solver, solver_list[_prm.solver_id], &_t->hsolver,
		       0., 0., 0, &_t->solver_setup, &_t->solver_solve,
		       &_t->solver_destroy,
		       &_t->solver_set_tol, &_t->solver_set_atol,
		       &_t->solver_set_maxiter, &_t->solver_setpc,
		       &_t->solver_get_nits, &_t->solver_get_resid, 1);
    if (_t->hsolver == NULL)
	phgError(1, "solver \"%s\" is unavailable or unimplemented.\n",
			solver_names[_prm.solver_id]);
    if (_t->solver_setpc != NULL) {
	/* setup preconditioner */
	setup_hypre_solver(solver, precon_list[_prm.precon_id], &_t->precond,
			   _prm.pc_rtol, _prm.pc_atol, _prm.pc_maxit,
			   &_t->precon_setup,
			   &_t->precon_solve, &_t->precon_destroy,
			   NULL, NULL, NULL, NULL, NULL, NULL, 0);
	if (_t->precon_solve != NULL) {
	    phgInfo(1, "setup HYPRE preconditioner \"%s\".\n",
			precon_names[_prm.precon_id]);
	    _t->solver_setpc(_t->hsolver,
			     (HYPRE_PtrToSolverFcn)_t->precon_solve,
			     (HYPRE_PtrToSolverFcn)_t->precon_setup,
			     _t->precond);
	}
    }
    else {
	phgInfo(1, "solver \"%s\" does not support preconditioner.\n",
			solver_names[_prm.solver_id]);
	_t->precond = NULL;
	_t->precon_destroy = NULL;
    }

    return 0;
}

static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    MAP *map = solver->mat->rmap;
    INT i, offset, *ni;
    HYPRE_Int num_iterations = 1, *pni;
    double res_norm = 0.;
#if FT_PHG != FT_DOUBLE
    double *v;
#endif	/* FT_PHG == FT_DOUBLE */

    HYPRE_ParCSRMatrix	par_A;
    HYPRE_ParVector	par_B, par_U;

    ni = phgAlloc(map->nlocal * sizeof(*ni));
    offset = map->partition[map->rank];
    for (i = 0; i < map->nlocal; i++)
	ni[i] = i + offset;
    if (sizeof(INT) != sizeof(HYPRE_Int)) {
	pni = phgAlloc(map->nlocal * sizeof(*pni));
	for (i = 0; i < map->nlocal; i++)
	    pni[i] = (HYPRE_Int)ni[i];
    }
    else {
	pni = (HYPRE_Int *)ni;
    }

#if FT_PHG != FT_DOUBLE
    v = phgAlloc(map->nlocal * sizeof(*v));
#endif	/* FT_PHG == FT_DOUBLE */

    if (_t->B == NULL) {
	/* copy data to RHS vector */
	AddRHSEntries(solver, map->nlocal, ni, solver->rhs->data);
    }

    Assemble(solver);

    /* set up initial solution */
    HYPRE_IJVectorCreate(g_comm, offset, offset + map->nlocal - 1, &_t->U);
    HYPRE_IJVectorSetObjectType(_t->U, HYPRE_PARCSR);
    HYPRE_IJVectorSetMaxOffProcElmts(_t->U, 0);
    HYPRE_IJVectorInitialize(_t->U);
#if FT_PHG == FT_DOUBLE
    HYPRE_IJVectorSetValues(_t->U, map->nlocal, pni, x->data);
#else	/* FT_PHG == FT_DOUBLE */
    for (i = 0; i < map->nlocal; i++)
	v[i] = x->data[i];
    HYPRE_IJVectorSetValues(_t->U, map->nlocal, pni, v);
#endif	/* FT_PHG == FT_DOUBLE */

    HYPRE_IJMatrixGetObject(_t->A, (void **)(void *)&par_A);
    dump_matrix(solver, par_A, "A");
    HYPRE_IJVectorAssemble(_t->B);
    HYPRE_IJVectorGetObject(_t->B, (void **)(void *)&par_B);
    dump_vector(solver, par_B, "b");
    HYPRE_IJVectorAssemble(_t->U);
    HYPRE_IJVectorGetObject(_t->U, (void **)(void *)&par_U);
    dump_vector(solver, par_U, "x0");

    if (solver->rtol > solver->btol)
	_t->solver_set_tol(_t->hsolver, solver->rtol);
    else
	_t->solver_set_tol(_t->hsolver, solver->btol);
    if (_t->solver_set_atol != NULL)
	_t->solver_set_atol(_t->hsolver, solver->atol);
    _t->solver_set_maxiter(_t->hsolver, solver->maxit);
    if (_t->solver_setup != NULL) {
	_t->solver_setup(_t->hsolver, par_A, par_B, par_U);
	/* only call setup once */
	_t->solver_setup = NULL;
    }
    _t->solver_solve(_t->hsolver, par_A, par_B, par_U);

    if (destroy) {
#if HYPRE_VERSION_MAJOR >= 2
	ams_destroy(solver); 
#endif	/* HYPRE_VERSION_MAJOR >= 2 */
	if (_t->precond != NULL && _t->precon_destroy != NULL) {
	    _t->precon_destroy(_t->precond);
	    _t->precond = NULL;
	    _t->precon_destroy = NULL;
	}
    }

    if (_t->solver_get_nits != NULL)
	_t->solver_get_nits(_t->hsolver, &num_iterations);
    if (_t->solver_get_resid != NULL)
	_t->solver_get_resid(_t->hsolver, &res_norm);

    HYPRE_IJVectorDestroy(_t->B);
    _t->B = NULL;

    if (destroy) { 
	_t->solver_destroy(_t->hsolver);
	_t->hsolver = NULL;
    }

    solver->residual = res_norm;
    solver->nits = num_iterations;

    if (destroy) {
	HYPRE_IJMatrixDestroy(_t->A);
	_t->A = NULL;
    }

    /* copy solution to DOFs */
#if FT_PHG == FT_DOUBLE
    HYPRE_IJVectorGetValues(_t->U, map->nlocal, pni, x->data);
#else	/* FT_PHG == FT_DOUBLE */
    HYPRE_IJVectorGetValues(_t->U, map->nlocal, pni, v);
    for (i = 0; i < map->nlocal; i++)
	x->data[i] = v[i];
    phgFree(v);
#endif	/* FT_PHG == FT_DOUBLE */
    if (sizeof(INT) != sizeof(HYPRE_Int))
	phgFree(pni);
    phgFree(ni);

    HYPRE_IJVectorDestroy(_t->U);
    _t->U = NULL;

    return solver->nits;
}

#define SetPC			NULL

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverHYPRE_ = {
    "HYPRE", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, FALSE, FALSE
};

#endif		/* USE_HYPRE */
