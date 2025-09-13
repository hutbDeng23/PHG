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

/* $Id: solver.h,v 1.24 2022/09/21 02:35:28 zlb Exp $ */

#ifndef PHG_SOLVER_H

/* The SOLVER struct */
typedef struct SOLVER_ {
    FLOAT	residual;	/* residual or error norm (filled by Solve() */
    FLOAT       rtol;		/* relative convergence tolerance */
    FLOAT       atol;		/* absolute tolerance */
    FLOAT       btol;		/* b-tolerance */
    double	cond;		/* spectral condition number (-solver_cond) */
    MAT		*mat, *mat_bak;	/* the matrix */
    VEC		*rhs;		/* local right-hand-side vector */
    FLOAT	*D;		/* data for the diagonal scaling matrix */
    struct OEM_SOLVER_	*oem_solver;
    void		*oem_data;
    struct SOLVER_PC_	*pc;	/* PreConditioner for solver */
    struct SOLVER_	*St;	/* the St solver for phgMatConditionNumber */

    int		symmetry;	/* symmetry (unsym, ssym, sym, spd, auto) */
    int		matrix_free;	/* allow, disallow, or enforce */
    int		nits;		/* number of iterations (filled by Solve() */
    int		maxit;		/* Maxi. number of iterations for oem_solver */
    BOOLEAN	monitor;	/* show convergence history (-solver_monitor) */
    BOOLEAN	scaling;	/* perform diagonal scaling (-solver_scaling) */
    BOOLEAN     warn_maxit;	/* TRUE ==> warn when maxit achieved */

    /* status flags */
    BOOLEAN	pc_option_flag;	/* TRUE ==> pc is set by cmdline option */
    BOOLEAN	rhs_updated;	/* TRUE ==> RHS updated with bdry eqns */
    BOOLEAN	sent_to_oem;	/* TRUE ==> data sent to OEM solver */
    BOOLEAN	oem_created;	/* TRUE ==> oem_solver->Create() called */
    BOOLEAN	assembled;	/* TRUE ==> solver already assembled */
} SOLVER;

/* values for the "-solver_matrix_free" option */
enum {ALLOW = 0, DISALLOW = 1, ENFORCE = 2};

/* The OEM_SOLVER struct */
typedef struct OEM_SOLVER_ {
    const char *name;				/* solver name */
    int (*RegisterOptions)(void);
    int (*Initialize)(int *argc, char ***argv); /* by phgSolverInitialize */
    int (*Finalize)(void);

    int (*Init)(SOLVER *solver);	/* called by phgSolverCreate */
    int (*Create)(SOLVER *solver);	/* called by phgSolverSolve */
    int (*Destroy)(SOLVER *solver);
    int (*AddMatrixEntries)(SOLVER *solver, INT nrows, INT *rows,
		INT ncols, INT *cols, FLOAT *values);
    int (*AddRHSEntries)(SOLVER *solver, INT n, INT *ni, FLOAT *values);
    int (*Assemble)(SOLVER *solver);
    int (*SetPC)(SOLVER *solver);
    int (*Solve)(SOLVER *solver, VEC *x, BOOLEAN destroy);

    /* optional user-callable functions (note: double pointers are returned
     * so that user can make changes to (e.g., reallocate) the objects) */
    void **(*GetSolver)(SOLVER *solver); /* returns pointer to OEM solver */
    void **(*GetMatrix)(SOLVER *solver); /* returns pointer to OEM matrix */
    void **(*GetPC)(SOLVER *solver);	 /* returns pointer to OEM precond */

    /* required symmetry */
    SYMMETRY req_symmetry;

    /* flags */
    BOOLEAN iterative;			/* whether an iterative solver */
    BOOLEAN parallel;			/* whether a parallel solver */
    BOOLEAN allow_matrix_free;		/* whether allow shell matrix */
    BOOLEAN allow_extern_pc;		/* whether use phgSolverSetPC() */
} OEM_SOLVER;

/* pre-conditioner */

/* preconditioner interface */
typedef void (*PC_PROC)(void *ctx, VEC *b, VEC **x);

/* SOLVER_PC struct */
typedef struct SOLVER_PC_ {
    char	*name;
    PC_PROC	pc_proc;	/* function performing the PC */
    void	*ctx;		/* context data (passed to PC_PROC) */
    SOLVER	*solver_owned;	/* solver to be freed on destruction */
} SOLVER_PC;

/* export phgSolverNames and phgSolverList to user programs */
extern const char *phgSolverNames_[];
extern OEM_SOLVER *phgSolverList_[];
#define phgSolverNames	(phgSolverNames_ + 1)
#define phgSolverList	(phgSolverList_ + 1)

/* default solver */
#define SOLVER_DEFAULT	((OEM_SOLVER *)NULL)	/* default solver */

/* PETSc solver */
#if USE_PETSC
extern OEM_SOLVER phgSolverPetsc_;
#define SOLVER_PETSC	(&phgSolverPetsc_)
#else
#define SOLVER_PETSC	((OEM_SOLVER *)NULL)
#endif

/* Trilinos solver */
#if USE_TRILINOS
extern OEM_SOLVER phgSolverTrilinos_;
#define SOLVER_TRILINOS	(&phgSolverTrilinos_)
#else
#define SOLVER_TRILINOS	((OEM_SOLVER *)NULL)
#endif

/* SuperLU solver */
#if USE_SUPERLU
extern OEM_SOLVER phgSolverSuperLU_;
#define SOLVER_SUPERLU	(&phgSolverSuperLU_)
#else
#define SOLVER_SUPERLU	((OEM_SOLVER *)NULL)
#endif

/* MUMPS solver */
#if USE_MUMPS
extern OEM_SOLVER phgSolverMUMPS_;
#define SOLVER_MUMPS	(&phgSolverMUMPS_)
#else
#define SOLVER_MUMPS	((OEM_SOLVER *)NULL)
#endif

/* PaStiX solver */
#if USE_PASTIX
extern OEM_SOLVER phgSolverPaStiX_;
#define SOLVER_PASTIX	(&phgSolverPaStiX_)
#else
#define SOLVER_PASTIX	((OEM_SOLVER *)NULL)
#endif

/* HIPS solver */
#if USE_HIPS
extern OEM_SOLVER phgSolverHIPS_;
#define SOLVER_HIPS	(&phgSolverHIPS_)
#else
#define SOLVER_HIPS	((OEM_SOLVER *)NULL)
#endif

/* X9AMG solver */
#if USE_X9AMG
extern OEM_SOLVER phgSolverX9AMG_;
#define SOLVER_X9AMG	(&phgSolverX9AMG_)
#else
#define SOLVER_X9AMG	((OEM_SOLVER *)NULL)
#endif

/* SPOOLES solver */
#if USE_SPOOLES
extern OEM_SOLVER phgSolverSPOOLES_;
#define SOLVER_SPOOLES	(&phgSolverSPOOLES_)
#else
#define SOLVER_SPOOLES	((OEM_SOLVER *)NULL)
#endif

/* PARDISO solver */
#if USE_PARDISO
extern OEM_SOLVER phgSolverPARDISO_;
#define SOLVER_PARDISO	(&phgSolverPARDISO_)
#else
#define SOLVER_PARDISO	((OEM_SOLVER *)NULL)
#endif

/* SuiteSparse solvers */
#if USE_SSPARSE
extern OEM_SOLVER phgSolverSSPARSE_;
#define SOLVER_SSPARSE	(&phgSolverSSPARSE_)
#else
#define SOLVER_SSPARSE	((OEM_SOLVER *)NULL)
#endif

/* MINRES solver */
#if USE_MINRES
extern OEM_SOLVER phgSolverMINRES_;
#define SOLVER_MINRES	(&phgSolverMINRES_)
#else
#define SOLVER_MINRES	((OEM_SOLVER *)NULL)
#endif

/* LASPack solver */
#if USE_LASPACK
extern OEM_SOLVER phgSolverLASPack_;
#define SOLVER_LASPACK	(&phgSolverLASPack_)
#else
#define SOLVER_LASPACK	((OEM_SOLVER *)NULL)
#endif

/* hSParse solver */
#if USE_HSPARSE
extern OEM_SOLVER phgSolverHSparse_;
#define SOLVER_HSPARSE	(&phgSolverHSparse_)
#else
#define SOLVER_HSPARSE	((OEM_SOLVER *)NULL)
#endif

/* HYPRE solver */
#if USE_HYPRE
extern OEM_SOLVER phgSolverHYPRE_;
#define SOLVER_HYPRE	(&phgSolverHYPRE_)
#else
#define SOLVER_HYPRE	((OEM_SOLVER *)NULL)
#endif

#if USE_HYPRE && HYPRE_VERSION_MAJOR >= 2
void phgSolverHypreAMSSetPoisson_(SOLVER *solver, MAT *alpha, MAT *beta,
				  BOOLEAN destroy);
void phgSolverHypreAMSSetPoisson(SOLVER *solver, DOF *alpha, DOF *beta);
void phgSolverHypreAMSSetConstantPoisson(SOLVER *solver,
						FLOAT alpha, FLOAT beta);
#else	/* USE_HYPRE && HYPRE_VERSION_MAJOR >= 2 */
# define phgSolverHypreAMSSetPoisson(solver, a, b)	    /* empty macro */
# define phgSolverHypreAMSSetConstantPoisson(solver, a, b)  /* empty macro */
#endif	/* USE_HYPRE && HYPRE_VERSION_MAJOR >= 2 */

/* Diagonal solver */
extern OEM_SOLVER phgSolverDiagonal_;
#define SOLVER_DIAGONAL	(&phgSolverDiagonal_)

/* Built-in solver */
extern OEM_SOLVER phgSolverPCG_;
#define SOLVER_PCG	(&phgSolverPCG_)

extern OEM_SOLVER phgSolverBiCGstab_;
#define SOLVER_BiCGstab	(&phgSolverBiCGstab_)

extern OEM_SOLVER phgSolverGMRES_;
#define SOLVER_GMRES	(&phgSolverGMRES_)

extern OEM_SOLVER phgSolverPreOnly_;
#define SOLVER_PreOnly	(&phgSolverPreOnly_)

extern OEM_SOLVER phgSolverASP_;
#define SOLVER_ASP	(&phgSolverASP_)

/* function prototype for building the prolongation/auxiliary matrix */
typedef MAT * (*BUILD_FUNC)(SOLVER *solver);

BUILD_FUNC phgSolverASPSetBuildP(BUILD_FUNC func);
BUILD_FUNC phgSolverASPSetBuildA(BUILD_FUNC func);
DOF *phgSolverASPGetAuxDof(SOLVER *solver);
void phgSolverASPSetCoefficients(SOLVER *solver, DOF *alpha, DOF *beta);
void phgSolverASPSetConstantCoefficients(SOLVER *solver, FLOAT a, FLOAT b);

extern OEM_SOLVER phgSolver2Grid_;
#define SOLVER_2Grid	(&phgSolver2Grid_)

BUILD_FUNC phgSolver2GridSetCoarseMatFunc(BUILD_FUNC f);
void phgSolver2GridSetCoarseGrid(GRID *g);


extern OEM_SOLVER phgSolverAMS_;
#define SOLVER_AMS	(&phgSolverAMS_)

void phgSolverAMSRemoveGradient(DOF **fptr, DOF *beta, int order);
void phgSolverAMSSetCoefficients(SOLVER *solver, DOF *alpha, DOF *beta);
void phgSolverAMSSetConstantCoefficients(SOLVER *solver, FLOAT a, FLOAT b);
void phgSolverAMSSetPoisson(SOLVER *solver, MAT *A_G, MAT *A_Pi);

extern OEM_SOLVER phgSolverGather_;
#define SOLVER_GATHER	(&phgSolverGather_)

extern OEM_SOLVER phgSolverASM_;
#define SOLVER_ASM	(&phgSolverASM_)

extern OEM_SOLVER phgSolverSmoother_;
#define SOLVER_SMOOTHER	(&phgSolverSmoother_)

#ifdef __cplusplus
extern "C" {
#endif

#define phgSolverMapD2L(solver, dof_no, index) \
	phgMapD2L((solver)->mat->cmap, dof_no, index)

#define phgSolverMapE2L(solver, dof_no, e, index) \
	phgMapE2L((solver)->mat->cmap, dof_no, e, index)

#define phgSolverMapL2G(solver, index) \
	phgMapL2G((solver)->mat->cmap, index)

#define phgSolverMapE2G(solver, dof_no, e, index) \
	phgMapE2G((solver)->mat->cmap, dof_no, e, index)

/* functions provided by solver interface */
int phgSolverRegisterOptions(void);
int phgSolverInitialize(int *argc, char ***argv);
int phgSolverFinalize(void);
void phgSolverSetDefaultSuboptions(void);

SOLVER *phgSolverMat2Solver(OEM_SOLVER *oem_solver, MAT *mat);
MAT *phgSolverGetMat(SOLVER *solver);
SOLVER *phgSolverCreate(OEM_SOLVER *oem_solver, DOF *u, ...);
int phgSolverDestroy(SOLVER **solver);

#define phgSolverAddMatrixEntry(solver, row, col, value) \
	phgMatAddEntry((solver)->mat, row, col, value)
#define phgSolverAddGlobalMatrixEntry(solver, row, col, value) \
	phgMatAddGlobalEntry((solver)->mat, row, col, value)

#define phgSolverAddMatrixEntries(solver, nrows, rows, ncols, cols, values) \
	phgMatAddEntries((solver)->mat, nrows, rows, ncols, cols, values)
#define phgSolverAddGlobalMatrixEntries(solver,nrows,rows,ncols,cols,values) \
	phgMatAddGlobalEntries((solver)->mat, nrows, rows, ncols, cols, values)

#define phgSolverAddRHSEntry(solver, index, value) \
	phgVecAddEntry((solver)->rhs, 0, index, value)
#define phgSolverAddGlobalRHSEntry(solver, index, value) \
	phgVecAddGlobalEntry((solver)->rhs, 0, index, value)

#define phgSolverAddRHSEntries(solver, n, indices, values) \
	phgVecAddEntries((solver)->rhs, 0, n, indices, values)
#define phgSolverAddGlobalRHSEntries(solver, n, indices, values) \
	phgVecAddGlobalEntries((solver)->rhs, 0, n, indices, values)

BOOLEAN phgSolverDenseLU(int n, FLOAT *a, int pvt[]);
void phgSolverDenseSV(int n, FLOAT *a, int pvt[], int m, FLOAT *b);
BOOLEAN phgSolverDenseSolver(int n, int m, FLOAT *a, FLOAT *b);

int phgSolverAssemble(SOLVER *solver);
#define phgSolverAssembleRHS(solver)	phgVecAssemble((solver)->rhs)

void phgSolverUpdateRHS(SOLVER *solver);
int phgSolverSolve(SOLVER *solver, BOOLEAN destroy, DOF *u, ...);
void phgSolverSetMaxIt(SOLVER *solver, int maxit);
void phgSolverSetTol(SOLVER *solver, FLOAT user_tol);

int phgSolverResetRHS(SOLVER *solver);
int phgSolverVecSolve(SOLVER *solver, BOOLEAN destroy, VEC *x);

void phgSolverSetPC_(SOLVER *solver, void *ctx, PC_PROC pc_func);
void phgSolverSetPC(SOLVER *solver, void *ctx, PC_PROC pc_func);

void phgSolverDumpMATLAB_(SOLVER *solver, const char *mat_name,
			  const char *rhs_name, const char *perm_name,
			  BOOLEAN reorder);
#define phgSolverDumpMATLAB(solver, mat_name, rhs_name) \
	phgSolverDumpMATLAB_(solver, mat_name, rhs_name, NULL, FALSE)
#ifdef __cplusplus
}
#endif

#define PHG_SOLVER_H
#endif
