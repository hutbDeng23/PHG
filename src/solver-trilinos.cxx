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

/* $Id: solver-trilinos.cxx,v 1.39 2022/09/21 02:35:28 zlb Exp $ */

/* Trilinos linear solvers interface */

#include "phg.h"

#if USE_TRILINOS		/* whole file */

#include <Epetra_ConfigDefs.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <ml_include.h>
#include <ml_MultiLevelOperator.h>
#include <ml_epetra_utils.h>
#if HAVE_TRILINOS_VERSION_H
# include <Trilinos_version.h>
#endif	/* HAVE_TRILINOS_VERSION_H */

#include <stdlib.h>
#if HAVE_FEENABLEEXCEPT
# include <fenv.h>
#endif	/* HAVE_FEENABLEEXCEPT */

// Note: see TRILINOS_INCLUDE/az_aztec_defs.h for parameters

static const char *solver_name[] = {
    "cg", "gmres", "cgs", "tfqmr", "bicgstab",
    "slu", "symmlq", "gmresr", "fixed_pt", "analyze",
    "lu", "cg_condnum", NULL};
static int solver_list[] = {
    AZ_cg, AZ_gmres, AZ_cgs, AZ_tfqmr, AZ_bicgstab,
    AZ_slu, AZ_symmlq, AZ_GMRESR, AZ_fixed_pt,
    AZ_analyze, AZ_lu, AZ_cg_condnum};
static int solver_index = 0;

static const char *pc_name[] = {
    "none", "jacobi", "sym_gs", "neumann", "ls", "ilu", "bilu", "lu",
    "icc", "ilut", "rilu", "recursive", "smoother", "dom_decomp",
    // Note: borrow the keyword AZ_multilevel for the ML preconditioner
    "ml", NULL};
static int pc_list[] = {
    AZ_none, AZ_Jacobi, AZ_sym_GS, AZ_Neumann, AZ_ls, AZ_ilu, AZ_bilu, AZ_lu,
    AZ_icc, AZ_ilut, AZ_rilu, AZ_recursive, AZ_smoother, AZ_dom_decomp,
    AZ_multilevel};
static int pc_index = 1;	/* default to Jacobi PC */

#if HAVE_FEENABLEEXCEPT
static BOOLEAN disable_fpetrap = TRUE;
#endif	/* HAVE_FEENABLEEXCEPT */
static BOOLEAN disable_errmsg = TRUE;

typedef struct {
    Epetra_MpiComm			*comm;
    Epetra_Map				*map;
    Epetra_Vector			*U, *B;
    Epetra_CrsMatrix			*A;
    Epetra_LinearProblem		*linsys;
    AztecOO				*solver;
    ML_Epetra::MultiLevelOperator	*mlop;
    ML_Aggregate			*agg_object;
    ML					*ml_handle;
    BOOLEAN				assembled;

    int					solver_index;
    int					pc_index;
#if HAVE_FEENABLEEXCEPT
    BOOLEAN				disable_fpetrap;
#endif	/* HAVE_FEENABLEEXCEPT */
    BOOLEAN				disable_errmsg;
} OEM_DATA;

/* convenience macro */
#define _t	((OEM_DATA *)solver->oem_data)

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nTrilinos AztecOO options:", "\n", "trilinos");
    phgOptionsRegisterKeyword("-aztec_solver", "AztecOO solver type",
				solver_name, &solver_index);
    phgOptionsRegisterKeyword("-aztec_pc", "AztecOO PC type",
				pc_name, &pc_index);
#if HAVE_FEENABLEEXCEPT
    phgOptionsRegisterNoArg("-aztec_disable_fpetrap", "Disable fpetrap when "
				"calling AztecOO functions", &disable_fpetrap);
#endif	/* HAVE_FEENABLEEXCEPT */
    phgOptionsRegisterNoArg("-aztec_disable_errmsg", "Disable (most) AztecOO "
				"error messages", &disable_errmsg);

    return 0;
}

#define Initialize	NULL
#define Finalize	NULL
#define SetPC		NULL

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));
    _t->solver_index = solver_index;
    _t->pc_index = pc_index;
#if HAVE_FEENABLEEXCEPT
    _t->disable_fpetrap = disable_fpetrap;
#endif	/* HAVE_FEENABLEEXCEPT */
    _t->disable_errmsg = disable_errmsg;

    _t->assembled = FALSE;

    return 0;
}

static int
Create(SOLVER *solver)
{
    _t->assembled = FALSE;

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    if (_t->mlop != NULL)
	delete _t->mlop;
    if (_t->agg_object != NULL)
	ML_Aggregate_Destroy(&_t->agg_object);
    if (_t->ml_handle != NULL)
	ML_Destroy(&_t->ml_handle);
    if (_t->solver != NULL)
	delete _t->solver;
    if (_t->linsys != NULL)
	delete _t->linsys;
    if (_t->A != NULL)
	delete _t->A;
    if (_t->B != NULL)
	delete _t->B;
    if (_t->U != NULL)
	delete _t->U;
    if (_t->map != NULL)
	delete _t->map;
    if (_t->comm != NULL)
	delete _t->comm;

    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
 
static int
az_printf_null(const char *format, ...)
{
    Unused(format);
    return 0;
}

static void
aztec_errmsg(BOOLEAN flag)
{
    // FIXME: use extern void AZOO_set_stream_err(std::ostream& ostrm);
    extern int (*azoo_printf_err)(const char *format, ...);
    static int (*azoo_printf_err_save)(const char *format, ...) = NULL;

    if (flag == FALSE) {
	/* disable error message */
	azoo_printf_err_save = azoo_printf_err;
	azoo_printf_err = az_printf_null;
    }
    else {
	/* restore error message */
	azoo_printf_err = azoo_printf_err_save;
    }
}

static int
Assemble(SOLVER *solver)
{
    MAP *map0 = solver->mat->rmap;
    int Nglobal = map0->nglobal, Nlocal = map0->nlocal;
    int i, n0;
    const MAT_ROW *row;
#if HAVE_FEENABLEEXCEPT
    int saved_excepts = 0;
#endif	/* HAVE_FEENABLEEXCEPT */

    if (_t->assembled == TRUE)
	return 0;

#if HAVE_FEENABLEEXCEPT
    if (_t->disable_fpetrap)
	saved_excepts = fedisableexcept(FE_UNDERFLOW | FE_DIVBYZERO |
					FE_INVALID | FE_OVERFLOW);
#endif	/* HAVE_FEENABLEEXCEPT */

    if (_t->disable_errmsg)
	aztec_errmsg(FALSE);

    _t->assembled = TRUE;
    assert(solver->mat->assembled);
    assert(solver->mat != NULL);
    assert(solver->sent_to_oem == FALSE);

    _t->comm = new Epetra_MpiComm(map0->comm);

    n0 = map0->partition[map0->rank];

    /* create map */
    int *work = new int[Nlocal];
    for (i = 0; i < Nlocal; i++)
	work[i] = i + n0;
    _t->map = new Epetra_Map(Nglobal, Nlocal, work, 0, *_t->comm);

    /* create vectors */
    _t->U = new Epetra_Vector(*_t->map);
    _t->B = new Epetra_Vector(*_t->map);

    /* create matrix */
    for (i = 0; i < Nlocal; i++)
	work[i] = phgMatGetNnzInRow(solver->mat, i, NULL, NULL);
    _t->A = new Epetra_CrsMatrix(Copy, *_t->map, work);
    delete [] work;
    /* send matrix entries */
    assert(sizeof(int) == sizeof(INT) && sizeof(double) == sizeof(FLOAT));
    for (i = 0; i < Nlocal; i++) {
	row = phgMatGetRow(solver->mat, i);
	if (row->ncols == 0)
	    continue;
	_t->A->InsertGlobalValues(n0 + i,
		row->ncols, (double *)row->data, (int *)row->cols);
	if (solver->mat->refcount == 0 && solver->mat->type == PHG_UNPACKED) {
	    phgFree(solver->mat->rows[i].data);
	    solver->mat->rows[i].data = NULL;
	    phgFree(solver->mat->rows[i].cols);
	    solver->mat->rows[i].cols = NULL;
	    solver->mat->rows[i].ncols = 0;
	}
    }
    if (solver->mat->refcount == 0)
	phgMatFreeMatrix(solver->mat);
    solver->sent_to_oem = TRUE;
    _t->A->FillComplete();

    /* define the linear problem */
    _t->linsys = new Epetra_LinearProblem(_t->A, _t->U, _t->B);
    _t->solver = new AztecOO(*_t->linsys);

    // set solver options
    int *options = new int[AZ_OPTIONS_SIZE];
    double *params = new double[AZ_PARAMS_SIZE];
    // AZ_defaults(options, params);

    /* retrieve internal options and parameters */
    memcpy(options, _t->solver->GetAllAztecOptions(),
			AZ_OPTIONS_SIZE * sizeof(*options));
    memcpy(params, _t->solver->GetAllAztecParams(),
			AZ_PARAMS_SIZE * sizeof(*params));

    if (phgVerbosity <= 0)
	options[AZ_output]	= AZ_none;
    if (_t->disable_errmsg)
	options[AZ_diagnostics] = AZ_none;
    options[AZ_conv]		= AZ_Anorm;
    options[AZ_solver]		= solver_list[_t->solver_index];

    int levels = 10;	/* ML levels */
    switch (pc_list[_t->pc_index]) {
	case AZ_multilevel:	/* set up ML preconditioner */
#if !defined(TRILINOS_MAJOR_VERSION) || TRILINOS_MAJOR_VERSION < 9
	    if (phgNProcs > 1) {
		/* The ML initialization code uses MPI_COMM_WORLD which
		 * may cause the program to hang if communicator changes */
		static int nprocs_save = 0;
		if (nprocs_save > 0 && map0->nprocs != nprocs_save)
		    phgPrintf("Warning! AZTEC/ML hangs if the number "
			  "of processes changes.\n", __func__);
		nprocs_save = map0->nprocs;
	    }
#endif	/* !defined(TRILINOS_MAJOR_VERSION) || TRILINOS_MAJOR_VERSION < 9 */
	    /* use ML_Epetra::MultiLevelOperator class */
	    ML_Set_PrintLevel(phgVerbosity);
	    ML_Create(&_t->ml_handle, levels);

	    ML_Set_Comm_Communicator(_t->ml_handle, map0->comm);
	    ML_Set_Comm_MyRank(_t->ml_handle, map0->rank);
	    ML_Set_Comm_Nprocs(_t->ml_handle, map0->nprocs);

	    ML_Set_Symmetrize(_t->ml_handle, FALSE);
	    ML_Set_Tolerance(_t->ml_handle, 0.0);
	    ML_Set_MaxIterations(_t->ml_handle, 1);

	    EpetraMatrix2MLMatrix(_t->ml_handle, 0, _t->A);
	    ML_Aggregate_Create(&_t->agg_object);
	    levels = ML_Gen_MGHierarchy_UsingAggregation(_t->ml_handle,
			0/*finest level*/, ML_INCREASING, _t->agg_object);
	    ML_Gen_Smoother_SymGaussSeidel(_t->ml_handle, ML_ALL_LEVELS,
			ML_BOTH, 1, ML_DEFAULT);
	    ML_Gen_Solver(_t->ml_handle, ML_MGV, 0, levels - 1);
	    _t->mlop = new ML_Epetra::MultiLevelOperator(_t->ml_handle,
			*_t->comm, *_t->map, *_t->map);
	    _t->solver->SetPrecOperator(_t->mlop);
	    break;
	case AZ_ilu:
	    options[AZ_precond] = AZ_dom_decomp;
	    options[AZ_subdomain_solve] = AZ_ilu;
	case AZ_lu:
	    options[AZ_precond] = AZ_dom_decomp;
	    options[AZ_subdomain_solve] = AZ_lu;
	case AZ_ilut:
	    options[AZ_precond] = AZ_dom_decomp;
	    options[AZ_subdomain_solve] = AZ_ilut;
	    break;
	default:
	    options[AZ_precond] = pc_list[_t->pc_index];
    }
    _t->solver->SetAllAztecOptions(options);
    _t->solver->SetAllAztecParams(params);

    delete [] params;
    delete [] options;

    if (_t->disable_errmsg)
	aztec_errmsg(TRUE);

#if HAVE_FEENABLEEXCEPT
    if (_t->disable_fpetrap) {
	feclearexcept(saved_excepts);
	feenableexcept(saved_excepts);
    }
#endif	/* HAVE_FEENABLEEXCEPT */

    return 0;
}

static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    int i;
    double *values;
#if HAVE_FEENABLEEXCEPT
    int saved_excepts = 0;
#endif	/* HAVE_FEENABLEEXCEPT */

    Assemble(solver);

#if HAVE_FEENABLEEXCEPT
    if (_t->disable_fpetrap)
	saved_excepts = fedisableexcept(FE_UNDERFLOW | FE_DIVBYZERO |
					FE_INVALID | FE_OVERFLOW);
#endif	/* HAVE_FEENABLEEXCEPT */

    if (_t->disable_errmsg)
	aztec_errmsg(FALSE);

    /* set RHS */
    _t->B->ExtractView(&values);
    for (i = 0; i < solver->mat->rmap->nlocal; i++)
	values[i] = (double)solver->rhs->data[i];

    /* set initial solution */
    _t->U->ExtractView(&values);
    for (i = 0; i < solver->mat->rmap->nlocal; i++)
	values[i] = (double)x->data[i];

    /* solve */
    if (solver->rtol > solver->btol)
	_t->solver->/*recursive*/Iterate(solver->maxit, solver->rtol);
    else
	_t->solver->/*recursive*/Iterate(solver->maxit, solver->btol);

    /* get nbr iters and final residual */
    //solver->residual = _t->solver->RecursiveResidual();
    //solver->residual = _t->solver->ScaledResidual();
    solver->residual = _t->solver->TrueResidual();
    solver->nits = _t->solver->NumIters();

    /* copy solution */
    for (i = 0; i < solver->mat->rmap->nlocal; i++)
	x->data[i] = values[i];

    if (_t->disable_errmsg)
	aztec_errmsg(FALSE);

#if HAVE_FEENABLEEXCEPT
    if (_t->disable_fpetrap) {
	feclearexcept(saved_excepts);
	feenableexcept(saved_excepts);
    }
#endif	/* HAVE_FEENABLEEXCEPT */

    if (destroy)
	Destroy(solver);

    return solver->nits;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverTrilinos_ = {
    "Trilinos", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy, AddMatrixEntries,
    AddRHSEntries, Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, TRUE, TRUE, FALSE, FALSE
};

#endif		/* USE_TRILINOS */
