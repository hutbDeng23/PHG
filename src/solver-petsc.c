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

/* $Id: solver-petsc.c,v 1.104 2022/09/21 02:35:28 zlb Exp $ */

/* The PETSc interface for PHG */

/* Note: "petscksp.h" must be included in "phg.h" before "mpi-utils.h" */
#define NEED_PETSCKSP_H

#include "phg.h"

#if USE_PETSC

#include "phg/petsc-utils.h"

#include <stdlib.h>
#include <string.h>
#include <strings.h>		/* bzero() */

static char *petsc_pc_opts = NULL;
static char *petsc_prefix = NULL;	/* For KSPSetOptionsPrefix() */

typedef struct {
    char	*pc_opts;
    char	*prefix;
    Vec		B;	/* B=RHS */
    Vec		U;	/* U=solution */
    Mat		A;	/* A=matrix */
    KSP		ksp;
    BOOLEAN	assembly_rhs;
    BOOLEAN	assembly_mat;
    VEC		*x, *y;	/* work vectors */
} OEM_DATA;

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)

static PetscErrorCode code;

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nPETSc options:", " \n"
  "Note: '-oem_options \"PETSc options\"' may be used to pass options to\n"
  "PETSc. Please check PETSc's manual or use '-oem_options \"-help\"' to\n"
  "show PETSc options (E.g.: '-oem_options \"-ksp_type cg -pc_type none\"')\n",
  "petsc");
    phgOptionsRegisterString("-petsc_prefix",
		"Prefix for the KSP options (passed to KSPSetOptionsPrefix)",
		&petsc_prefix);
    phgOptionsRegisterString("-petsc_pc_opts",
		"Preconditioner options (if not null then PHG's solver will "
		"be used as the preconditioner with the string as its options)",
		&petsc_pc_opts);

    return 0;
}

static int
Initialize(int *argc, char ***argv)
{
    static char help[] = "\nThe PETSc interface for PHG.\n\n";

    /* Hacky: also initialize SLEPc here */
#if USE_SLEPC
    extern PetscErrorCode SlepcInitialize(int*, char***, char[], const char[]);
    code = SlepcInitialize(argc, argv, NULL, help); CHKERRQ(code);
#else
    code = PetscInitialize(argc, argv, NULL, help); CHKERRQ(code);
#endif

    /* print information about the solver */
    if (phgRank == 0) {
	PetscTruth flg;

	/* check for obsolete PETSc options */
#if PETSC_VERSION_MAJOR<3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<7)
	PetscOptionsHasName(PETSC_NULL,"-ksp_max_it", &flg);
#else	/* PETSC_VERSION < 3.7 */
	PetscOptionsHasName(NULL, PETSC_NULL,"-ksp_max_it", &flg);
#endif	/* PETSC_VERSION < 3.7 */
	if (flg)
	    phgPrintf("WARNING: PETSc option \"-ksp_max_it\" is obsolete, "
		      "use PHG option \"-solver_max_it\" instead.\n");

#if PETSC_VERSION_MAJOR<3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<7)
	PetscOptionsHasName(PETSC_NULL,"-ksp_rtol", &flg);
#else	/* PETSC_VERSION < 3.7 */
	PetscOptionsHasName(NULL, PETSC_NULL,"-ksp_rtol", &flg);
#endif	/* PETSC_VERSION < 3.7 */
	if (flg)
	    phgPrintf("WARNING: PETSc option \"-ksp_rtol\" is obsolete, "
		      "use PHG option \"-solver_rtol\" instead.\n");
	if (phgVerbosity > 0)
	    printf("sizeof(PetscInt) = %d, sizeof(PetscBLASInt) = %d\n",
			(int)sizeof(PetscInt), (int)sizeof(PetscBLASInt));
    }

    return 0;
}

static int
Finalize(void)
{
    if (!PetscInitializeCalled)
	return 0;

#if USE_SLEPC
    extern PetscErrorCode SlepcFinalize(void);
    code = SlepcFinalize(); CHKERRQ(code);
#else
    code = PetscFinalize(); CHKERRQ(code);
#endif

    return 0;
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));
    _t->pc_opts = (petsc_pc_opts == NULL ? NULL : strdup(petsc_pc_opts));
    _t->prefix = (petsc_prefix == NULL ? NULL : strdup(petsc_prefix));
    _t->B = PETSC_NULL;
    _t->U = PETSC_NULL;
    _t->A = PETSC_NULL;
    _t->ksp = PETSC_NULL;
    _t->assembly_rhs = FALSE;
    _t->assembly_mat = FALSE;
    _t->x = _t->y = NULL;

    solver->pc_option_flag = phgOptionsIfUsed("-petsc_pc_opts");

    return 0;
}

static int
to_petsc(VEC *x, Vec v)
/* copy PHG VEC 'x' to Petsc Vec 'v' */
{
    PetscScalar *vec;
    PetscInt size;
    INT i;

    if (v == PETSC_NULL)
	phgError(1, "PETSc CopySolution: invalid PETSc vector.\n");

    code = VecGetLocalSize(v, &size); CHKERRQ(code);
    code = VecGetArray(v, &vec); CHKERRQ(code);
    for (i = 0; i < size; i++)
	vec[i] = x->data[i];
    code = VecRestoreArray(v, &vec); CHKERRQ(code);

    return 0;
}

static int
from_petsc(VEC *x, Vec v)
/* copy Petsc Vec 'v' to PHG VEC 'x' */
{
    PetscScalar *vec;
    PetscInt size;
    INT i;

    if (v == PETSC_NULL)
	phgError(1, "PETSc CopySolution: invalid PETSc vector.\n");

    code = VecGetLocalSize(v, &size); CHKERRQ(code);
    code = VecGetArray(v, &vec); CHKERRQ(code);
    for (i = 0; i < size; i++)
	x->data[i] = vec[i];
    code = VecRestoreArray(v, &vec); CHKERRQ(code);

    return 0;
}

static MAT_OP mat_op = MAT_OP_N;

static int
matop_mult(Mat mat, Vec x, Vec y)
{
    SOLVER *solver;

    /* Note: ctx points to PHG SOLVER struct */
    MatShellGetContext(mat, (void **)(void *)&solver);
    if (_t->x == NULL)
	_t->x = phgMapCreateVec(solver->rhs->map, 1);
    from_petsc(_t->x, x);
    phgMatVec(mat_op, 1.0, solver->mat, _t->x, 0.0, &_t->y);
    to_petsc(_t->y, y);

    return 0;
}

static int
matop_mult_transpose(Mat mat, Vec x, Vec y)
{
    mat_op = MAT_OP_T;
    matop_mult(mat, x, y);
    mat_op = MAT_OP_N;

    return 0;
}

static int
matop_get_diagonal(Mat mat, Vec x)
{
    SOLVER *solver;
    PetscScalar *vec;
    INT i;

    /* Note: ctx points to PHG SOLVER struct */
    MatShellGetContext(mat, (void **)(void *)&solver);
    if (solver->mat->diag == NULL)
	phgMatSetupDiagonal(solver->mat);
    VecGetArray(x, &vec);
    for (i = 0; i < solver->rhs->map->nlocal; i++)
	vec[i] = solver->mat->diag[i];
    VecRestoreArray(x, &vec);

    return 0;
}

static int
Create(SOLVER *solver)
{
    MAP *map = solver->mat->rmap;
    PetscInt N = map->nglobal, Nlocal = map->nlocal;
    PetscInt i, prealloc, *d_nnz, *o_nnz;
    INT d, o;

    assert(solver->mat->assembled);

    if (solver->pc == NULL && solver->pc_option_flag &&
	_t->pc_opts != NULL && strcmp(_t->pc_opts, "none") != 0) {
	/* Create PC solver */
	SOLVER *pc_solver;
	phgOptionsPush();
	phgSolverSetDefaultSuboptions();
	phgOptionsSetOptions("-solver_rtol 0 -solver_atol 0 -solver_btol 0"
				 " -solver_maxit 1 +solver_warn_maxit");
	phgOptionsSetOptions(_t->pc_opts);
	pc_solver = phgSolverMat2Solver(SOLVER_DEFAULT, solver->mat);
	phgOptionsPop();
	if (pc_solver->oem_solver == solver->oem_solver) {
	    /* avoid regression */
	    phgFree(((OEM_DATA *)pc_solver->oem_data)->pc_opts);
	    ((OEM_DATA *)pc_solver->oem_data)->pc_opts = NULL;
	}
	phgSolverSetPC_(solver, pc_solver, NULL);
	solver->pc->solver_owned = pc_solver; /* destroy pc_solver afterwards */
    }

    code = VecCreate(map->comm, &_t->U); CHKERRQ(code);
    code = VecSetSizes(_t->U, Nlocal, N); CHKERRQ(code);
    code = VecSetFromOptions(_t->U); CHKERRQ(code);
    code = VecDuplicate(_t->U, &_t->B); CHKERRQ(code);

    if ((solver->matrix_free == ALLOW && solver->mat->type == PHG_MATRIX_FREE)
	|| solver->matrix_free == ENFORCE) {
	/* matrix-free matrix */
	code = MatCreateShell(map->comm, Nlocal, Nlocal, N, N, solver, &_t->A);
	CHKERRQ(code);
	code = MatShellSetOperation(_t->A, MATOP_MULT,
					(void(*)(void))matop_mult);
	code = MatShellSetOperation(_t->A, MATOP_MULT_TRANSPOSE,
					(void(*)(void))matop_mult_transpose);
	code = MatShellSetOperation(_t->A, MATOP_GET_DIAGONAL,
					(void(*)(void))matop_get_diagonal);
    }
    else {
	/* CSR matrix */
	prealloc = PETSC_DECIDE;
	d_nnz = phgAlloc(Nlocal * sizeof(*d_nnz));
	o_nnz = phgAlloc(Nlocal * sizeof(*o_nnz));
	for (i = 0; i < Nlocal; i++) {
	    phgMatGetNnzInRow(solver->mat, i, &d, &o);
	    d_nnz[i] = d;
	    o_nnz[i] = o;
	}
	code = MatCreateMPIAIJ(map->comm, Nlocal, Nlocal, N, N,
				prealloc, d_nnz, prealloc, o_nnz, &_t->A);
	phgFree(d_nnz);
	phgFree(o_nnz);
    }
    CHKERRQ(code);
    code = MatSetFromOptions(_t->A); CHKERRQ(code);

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    if (_t->pc_opts != NULL)
	phgFree(_t->pc_opts);
    if (_t->ksp != PETSC_NULL)
	phgPetscKSPDestroy(&_t->ksp);
    if (_t->A != PETSC_NULL)
	phgPetscMatDestroy(&_t->A);
    if (_t->B != PETSC_NULL)
	phgPetscVecDestroy(&_t->B);
    if (_t->U != PETSC_NULL)
	phgPetscVecDestroy(&_t->U);

    if (_t->x != NULL)
	phgVecDestroy(&_t->x);
    if (_t->y != NULL)
	phgVecDestroy(&_t->y);

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
AddMatrixEntries(SOLVER *solver, INT nrows, INT *rows, INT ncols, INT *cols,
			FLOAT *values)
{
    PetscInt *prows, *pcols;

    if (nrows == 0 || ncols == 0)
	return 0;

    if (_t->A == PETSC_NULL)
	phgError(1, "PETSc AddMatrixEntries: invalid matrix!\n");

    if (sizeof(INT) != sizeof(PetscInt)) {
	INT i;
	prows = phgAlloc((nrows + ncols) * sizeof(*prows));
	pcols = prows + nrows;
	for (i = 0; i < nrows; i++)
	    prows[i] = (PetscInt)rows[i];
	for (i = 0; i < ncols; i++)
	    pcols[i] = (PetscInt)cols[i];
    }
    else {
	prows = (PetscInt *)rows;
	pcols = (PetscInt *)cols;
    }

    if (sizeof(FLOAT) == sizeof(PetscScalar)) {
	code = MatSetValues(_t->A, nrows, prows, ncols, pcols,
				(PetscScalar *)values, INSERT_VALUES);
    }
    else {
	PetscInt i, n = nrows * ncols;
	PetscScalar *buffer = phgAlloc(n * sizeof(*buffer));
	for (i = 0; i < n; i++)
	    buffer[i] = (PetscScalar)values[i];
	code = MatSetValues(_t->A, nrows, prows, ncols, pcols,
				buffer, INSERT_VALUES);
	phgFree(buffer);
    }

    if (sizeof(INT) != sizeof(PetscInt))
	phgFree(prows);

    CHKERRQ(code);
    _t->assembly_mat = TRUE;

    return 0;
}

static int
AddRHSEntries(SOLVER *solver, INT n, INT *ni, FLOAT *values)
{
    PetscInt *pni;

    if (n == 0)
	return 0;

    if (_t->B == PETSC_NULL)
	phgError(1, "PETSc AddRHSEntries: invalid vector!\n");

    if (sizeof(INT) != sizeof(PetscInt)) {
	INT i;
	pni = phgAlloc(n * sizeof(*pni));
	for (i = 0; i < n; i++)
	    pni[i] = (PetscInt)ni[i];
    }
    else {
	pni = (PetscInt *)ni;
    }

    if (sizeof(FLOAT) == sizeof(PetscScalar)) {
	code = VecSetValues(_t->B, n, pni, (PetscScalar*)values, INSERT_VALUES);
    }
    else {
	PetscInt i;
	PetscScalar *buffer = phgAlloc(n * sizeof(*buffer));
	for (i = 0; i < n; i++)
	    buffer[i] = (PetscScalar)values[i];
	code = VecSetValues(_t->B, n, pni, buffer, INSERT_VALUES);
	phgFree(buffer);
    }

    if (sizeof(INT) != sizeof(PetscInt))
	phgFree(pni);

    CHKERRQ(code);
    _t->assembly_rhs = TRUE;

    return 0;
}
 
static PetscErrorCode
#if PETSC_VERSION_MAJOR<3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<1)
 ShellPreconditioner(void *ctx, Vec x, Vec y)
#else	/* PETSc version <= 3.0 */
 ShellPreconditioner(PC petsc_pc, Vec x, Vec y)
#endif	/* PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 1 */
{
    SOLVER_PC *phg_pc;
    PetscScalar *vec_x, *vec_y;
    VEC *x_phg, *y_phg;
    INT i;
#if PETSC_VERSION_MAJOR>3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>0)
    void *ctx;
    PCShellGetContext(petsc_pc, &ctx);
#endif	/* ! (PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 1) */

    phg_pc = ((SOLVER*)ctx)->pc;

    /* copy petsc-vec to phg-vec */
    code = VecGetArray(x, &vec_x); CHKERRQ(code);
    code = VecGetArray(y, &vec_y); CHKERRQ(code);
    if (sizeof(PetscScalar) == sizeof(FLOAT)) {
	x_phg=phgMapCreateVec(((SOLVER*)ctx)->mat->rmap, 0);
	y_phg=phgMapCreateVec(((SOLVER*)ctx)->mat->rmap, 0);
	x_phg->data = (FLOAT *)vec_x;
	y_phg->data = (FLOAT *)vec_y;
	x_phg->nvec = 1;
	y_phg->nvec = 1;
    }
    else {
	x_phg=phgMapCreateVec(((SOLVER*)ctx)->mat->rmap, 1);
	y_phg=phgMapCreateVec(((SOLVER*)ctx)->mat->rmap, 1);
	for (i = 0; i < x_phg->map->nlocal; i++) {
	    x_phg->data[i] = (FLOAT)vec_x[i];
	    y_phg->data[i] = (FLOAT)vec_y[i];
	}
    }

    /*  Precond Proc*/
    phg_pc->pc_proc(phg_pc->ctx, x_phg, &y_phg);

    /* copy phg-cvec to petsc-vec */
    if (sizeof(PetscScalar) == sizeof(FLOAT)) {
	x_phg->data = NULL;
	y_phg->data = NULL;
    }
    else {
	for (i = 0; i < x_phg->map->nlocal; i++) {
	    vec_x[i] = (PetscScalar)x_phg->data[i];
	    vec_y[i] = (PetscScalar)y_phg->data[i];
	}
    }
    code = VecRestoreArray(x, &vec_x); CHKERRQ(code);
    code = VecRestoreArray(y, &vec_y); CHKERRQ(code);
    phgVecDestroy(&x_phg);
    phgVecDestroy(&y_phg);

    return code;
}

static int
Assemble(SOLVER *solver)
{
    if (_t->assembly_rhs) {
	code = VecAssemblyBegin(_t->B); CHKERRQ(code);
	code = VecAssemblyEnd(_t->B); CHKERRQ(code);
	_t->assembly_rhs = FALSE;
    }

    if (_t->assembly_mat) {
	code = MatAssemblyBegin(_t->A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
	code = MatAssemblyEnd(_t->A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
	_t->assembly_mat = FALSE;
    }

    if (_t->ksp == PETSC_NULL) {
	code = KSPCreate(solver->mat->rmap->comm, &_t->ksp); CHKERRQ(code);
	CHKERRQ(code);
	if (_t->prefix != NULL)
	    KSPSetOptionsPrefix(_t->ksp, _t->prefix); 
	
#if PETSC_VERSION_MAJOR<3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<5)
	code = KSPSetOperators(_t->ksp, _t->A, _t->A,
				DIFFERENT_NONZERO_PATTERN);
#else	/* PETSc < 3.5 */
	code = KSPSetOperators(_t->ksp, _t->A, _t->A);
#endif	/* PETSc < 3.5 */
	CHKERRQ(code);
#if 0	/* The following line is disabled since it breaks the preonly solver */
	code = KSPSetInitialGuessNonzero(_t->ksp, PETSC_TRUE); CHKERRQ(code);
	CHKERRQ(code);
#endif

	/* user-defined prcond */ 
 	if(solver->pc!=NULL){
	    PC pc;
	    KSPGetPC(_t->ksp, &pc);
	    PCSetType(pc, PCSHELL);
#if PETSC_VERSION_MAJOR < 2 || \
    (PETSC_VERSION_MAJOR == 2 && PETSC_VERSION_MINOR < 2)
	    PCShellSetApply(pc, ShellPreconditioner, (void*)solver);
#else	/* PETSC_VERSION_MAJOR < 2 || ... */
	    PCShellSetContext(pc, (void*)solver);
	    PCShellSetApply(pc, ShellPreconditioner);
#endif	/* PETSC_VERSION_MAJOR < 2 || ... */
 	}

	code = KSPSetFromOptions(_t->ksp); CHKERRQ(code);
	CHKERRQ(code);

#if 0
	if (phgRank == 0) {
	    PetscTruth flg, flg1;
	    char name[32], name1[32];
	    PetscOptionsGetString(NULL, "-ksp_type", name, sizeof(name), &flg);
	    PetscOptionsGetString(NULL, "-pc_type", name1, sizeof(name), &flg1);
	    phgPrintf("KSP method = %s, precoditioner = %s\n",
			flg ? name : "default", flg1 ? name1 : "default");
    	}
#endif
    }
    
#if 0
    phgPrintf("The RHS:\n");
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
    VecView(_t->B, PETSC_VIEWER_STDOUT_SELF);
    phgPrintf("The matrix:\n");
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
    MatView(_t->A, PETSC_VIEWER_STDOUT_SELF);
#endif

    return 0;
}

#define set_initial_solution(solver, x)		to_petsc(x, _t->U)
#define copy_solution(solver, x)		from_petsc(x, _t->U)

/* RHS norm, initial and previous residual */
static FLOAT bnorm, residual0, residual1;

static PetscErrorCode
conv_test(KSP ksp, PetscInt it, PetscReal residual, KSPConvergedReason *reason,
	  void *phg_solver)
/* Note: set *reason to 0 if not converged, otherwise set *reason to one of
 *	 KSP_CONVERGED_RTOL and KSP_CONVERGED_ATOL */
{
    double rres, ratio;
    SOLVER *solver = phg_solver;

    if (residual0 == 0.) {
	ratio = rres = 1.0;
    }
    else {
#if 1
	/* work around for a possible icc bug on Altix 4700:
	 *   without the next statement "mpirun -np 1 poisson -solver petsc"
	 *   would produce a "divide by 0" exception */
	if (phgVerbosity > 1000000)
	    phgPrintf("");
#endif
	rres = residual / residual0;
	ratio = residual / residual1;
    }

    if (solver->monitor)
	phgPrintf("% 5d   %12le   %12le   %12le   %0.8lf\n", it,
			(double)residual, rres,
			(double)(residual / bnorm), ratio);

    if (residual <= solver->atol) {
	*reason = KSP_CONVERGED_ATOL;
	return 0;
    }

    if (residual0 <= 0.)
	residual0 = residual;

    if (residual <= solver->rtol * residual0) {
	*reason = KSP_CONVERGED_RTOL;
	return 0;
    }

    if (residual <= solver->btol * bnorm) {
	*reason = KSP_CONVERGED_RTOL;	/* note: no KSP_CONVERGED_BTOL */
	return 0;
    }

    residual1 = residual;
    *reason = (KSPConvergedReason)0;

    return 0;
}

static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    PetscInt nits;
    PetscReal residual;

    if (_t->A == PETSC_NULL ||
	_t->B == PETSC_NULL ||
	_t->U == PETSC_NULL)
	phgError(1, "PETSc Solve: invalid linear system.\n");

    code = Assemble(solver); CHKERRQ(code);

    if (solver->monitor) {
	phgPrintf("\n======================= PETSc ====================\n");
	phgPrintf("Tolerances: atol=%lg, rtol=%lg, btol=%lg\n",
	    (double)solver->atol, (double)solver->rtol, (double)solver->btol);
	phgPrintf("Iters   Resid. norm    Resid./resid0  Resid./rhs     "
		  "Conv. ratio\n");
	phgPrintf("-----   ------------   -------------  ------------   "
		  "-----------\n");
    }

    /* initialize constants used by conv_test().
     * FIXME: should use the preconditioned RHS for bnorm */
    residual0 = residual1 = 0.0;
    bnorm = phgVecNorm2(solver->rhs, 0, NULL);
    if (bnorm == 0.0) {
	/* set the solution to 0 and return */
	bzero(x->data, x->map->nlocal * sizeof(*x->data));
	nits = 0;
	residual = 0.;
	goto skip;
    }

    set_initial_solution(solver, x);

#if 0
    {
	PetscReal norma, normb, normu;
	MatNorm(_t->A, NORM_FROBENIUS, &norma);
	VecNorm(_t->B, NORM_2, &normb);
	VecNorm(_t->U, NORM_2, &normu);
	phgPrintf("\nNorm(A) = %lg, Norm(b) = %lg, Norm(x0) = %lg\n",
			(double)norma, (double)normb, (double)normu);
    }
#endif

    code = KSPSetTolerances(_t->ksp, solver->rtol, solver->atol,
			    PETSC_DEFAULT, solver->maxit);
    CHKERRQ(code);

    code = KSPSetConvergenceTest(_t->ksp,	/* KSP */
				 conv_test,	/* convergence test function */
				 solver		/* userdata pointer */
#if PETSC_VERSION_MAJOR >= 3 
				 , NULL		/* func destroying userdata */
#endif	/* PETSC_VERSION_MAJOR >= 3 */
				);
    CHKERRQ(code);

#if PETSC_VERSION_MAJOR < 2 || \
    (PETSC_VERSION_MAJOR == 2 && (PETSC_VERSION_MINOR < 2 || \
     (PETSC_VERSION_MINOR == 2 && PETSC_VERSION_SUBMINOR < 1)))
    code = KSPSetRhs(_t->ksp, _t->B); CHKERRQ(code);
    code = KSPSetSolution(_t->ksp, _t->U); CHKERRQ(code);
    code = KSPSolve(_t->ksp); CHKERRQ(code);
#else
    code = KSPSolve(_t->ksp, _t->B, _t->U); CHKERRQ(code);
#endif
    code = KSPGetIterationNumber(_t->ksp, &nits); CHKERRQ(code);
    code = KSPGetResidualNorm(_t->ksp, &residual); CHKERRQ(code);

    copy_solution(solver, x);

skip:
    if (destroy) {
	phgPetscKSPDestroy(&_t->ksp);
	_t->ksp = PETSC_NULL;
	phgPetscMatDestroy(&_t->A);
	_t->A = PETSC_NULL;
	phgPetscVecDestroy(&_t->B);
	_t->B = PETSC_NULL;
	phgPetscVecDestroy(&_t->U);
	_t->U = PETSC_NULL;

	if (_t->x != NULL)
	    phgVecDestroy(&_t->x);
	if (_t->y != NULL)
	    phgVecDestroy(&_t->y);
    }

    solver->residual = residual;
    solver->nits = nits;

    return nits;
}

static void **
GetSolver(SOLVER *solver)
/* returns PETSc KSP to the user */
{
    if (solver == NULL || solver->oem_data == NULL)
	return NULL;
    return (void *)&_t->ksp;
}

static void **
GetMatrix(SOLVER *solver)
/* returns PETSc MAT to the user */
{
    if (solver == NULL || solver->oem_data == NULL)
	return NULL;
    return (void *)&_t->A;
}

/*--------------------------------------------------------------------*/

#define SetPC			NULL

OEM_SOLVER phgSolverPetsc_ = {
    "PETSc", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, GetSolver, GetMatrix, NULL,
    M_UNSYM, TRUE, TRUE, TRUE, TRUE
};

#endif		/* USE_PETSC */
