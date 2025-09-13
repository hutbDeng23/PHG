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

/* $Id: solver-mumps.c,v 1.95 2022/09/21 02:35:28 zlb Exp $ */

/* The MUMPS interface (http://mumps.enseeiht.fr/) */

#define _GNU_SOURCE	/* for feenableexcept() */

#include "phg.h"

#if USE_MUMPS	/* whole file */

#define Initialize		NULL
#define Finalize		NULL
#define AddMatrixEntries	NULL
#define AddRHSEntries		NULL
#define SetPC			NULL

#include <smumps_c.h>
#include <dmumps_c.h>

#include <stdlib.h>
#include <string.h>

#if HAVE_FEENABLEEXCEPT
# include <fenv.h>
#endif  /* HAVE_FEENABLEEXCEPT */

/* macros s.t. indices match documentation */
#define ICNTL(I)	icntl[(I)-1]
#define CNTL(I)		cntl[(I)-1]
#define INFO(I)		info[(I)-1]
#define INFOG(I)	infog[(I)-1]

typedef enum {MUMPS_SINGLE = 0, MUMPS_DOUBLE = 1} PRECISION;

typedef struct {
    SMUMPS_STRUC_C	*smumps_id;
    DMUMPS_STRUC_C	*dmumps_id;
    BOOLEAN		factored;

    PRECISION		precision;
    int			ordering;	/* ICNTL(7), ordering_values[] */
    int			analysis_type;		/* ICNTL(28), 0, 1, 2 */
    int			parallel_ordering;	/* ICNTL(29), 0, 1, 2 */
    INT			scalapack_hack;
    INT			refine;
    INT			fillin;
    char		*filename;
#if HAVE_FEENABLEEXCEPT
    BOOLEAN		disable_fpetrap;
#endif  /* HAVE_FEENABLEEXCEPT */
} OEM_DATA;

#define	_t	((OEM_DATA *)solver->oem_data)
#define sid	(_t->smumps_id)
#define did	(_t->dmumps_id)

static const char *precision_keywords[] = {"single", "double", NULL};
int precision = MUMPS_DOUBLE;	/* use double precision by default */
static INT scalapack_hack = 0;	/* <=0: use ScaLAPACK, >0: disable ScaLAPACK */
static const char *ordering_keywords[] =
    {"auto", "amd", "amf", "scotch", "pord", "metis", "qamd", NULL};
static const int ordering_values[] = {7, 0, 2, 3, 4, 5, 6};
static int ordering = 0;	/* ordering scheme */
static const char *analysis_type_keywords[] =
    {"auto", "sequential", "parallel", NULL};
static int analysis_type = 0;
static const char *parallel_ordering_keywords[] =
    {"auto", "ptscotch", "parmetis", NULL};
static int parallel_ordering = 0;
static INT refine = 0;		/* max nbr of iterative refinements */

static INT fillin = 50;		/* percentage of extra fillin. Note: MUMPS's
				 * built-in default is 15 (SPD) and 20 (gen.) */

static char *filename = NULL;	/* filename to write the matrix */

#if HAVE_FEENABLEEXCEPT
static BOOLEAN disable_fpetrap = TRUE;
#endif  /* HAVE_FEENABLEEXCEPT */

static int
RegisterOptions(void)
{
    phgOptionsRegisterTitle("\nMUMPS options:", "\n", "mumps");
    phgOptionsRegisterKeyword("-mumps_precision", "Data type (precision)",
				precision_keywords, &precision);
    phgOptionsRegisterKeyword("-mumps_ordering", "Ordering scheme (only "
			      "effective when not doing analysis and symbolic "
			      "factorization in parallel)",
				ordering_keywords, &ordering);
    phgOptionsRegisterKeyword("-mumps_analysis_type",
			      "Whether perform the analysis and symbolic "
			      "factorization phase in parallel",
				analysis_type_keywords, &analysis_type);
    phgOptionsRegisterKeyword("-mumps_parallel_ordering",
			      "Type of parallel ordering scheme "
			      "(only effective when doing analysis and "
			      "symbolic factorization in parallel)",
				parallel_ordering_keywords,
				&parallel_ordering);

    phgOptionsRegisterInt("-mumps_refine", "Max number of iterative "
				"refinements", &refine);
    phgOptionsRegisterInt("-mumps_fillin", "Extra fill-in (percent, "
				"< 0 means use default setting)", &fillin);
#if HAVE_FEENABLEEXCEPT
    phgOptionsRegisterNoArg("-mumps_disable_fpetrap", "Disable FPE trap",
				&disable_fpetrap);
#endif  /* HAVE_FEENABLEEXCEPT */

    phgOptionsRegisterFilename("-mumps_dump_matrix",
			       "Dump the matrix to given file", &filename);

    phgOptionsRegisterInt("-mumps_scalapack_hack", "Whether use ScaLAPACK "
				"(<=0: use ScaLAPACK, >0: disable ScaLAPACK)",
				&scalapack_hack);
    return 0;
}

static void
call_mumps(SOLVER *solver)
{
#if HAVE_FEENABLEEXCEPT
    int saved_excepts = 0;
    if (_t->disable_fpetrap)
	saved_excepts = fedisableexcept(FE_UNDERFLOW | FE_DIVBYZERO |
					FE_INVALID | FE_OVERFLOW);
#endif  /* HAVE_FEENABLEEXCEPT */

    if (_t->precision == MUMPS_DOUBLE)
	dmumps_c(did);
    else
	smumps_c(sid);

#if HAVE_FEENABLEEXCEPT
    if (_t->disable_fpetrap) {
	feclearexcept(saved_excepts);
	feenableexcept(saved_excepts);
    }
#endif  /* HAVE_FEENABLEEXCEPT */
}

static int
Init(SOLVER *solver)
{
    solver->oem_data = phgCalloc(1, sizeof(OEM_DATA));

    _t->scalapack_hack	= scalapack_hack;
    _t->ordering	= ordering;
    _t->analysis_type	= analysis_type;
    _t->parallel_ordering = parallel_ordering;
    _t->refine		= refine;
    _t->fillin		= fillin;
    _t->precision	= precision;
    if (filename != NULL)
	_t->filename	= strdup(filename);
#if HAVE_FEENABLEEXCEPT
    _t->disable_fpetrap	= disable_fpetrap;
#endif  /* HAVE_FEENABLEEXCEPT */

    return 0;
}

static int
Create(SOLVER *solver)
{
    if (solver->monitor)
	phgPrintf("*** Creating MUMPS solver: symmetry = %s, precision = %s\n",
			solver->symmetry < M_SYM ? "unsym" :
			solver->symmetry == M_SYM ? "sym" : "spd",
			_t->precision == MUMPS_DOUBLE ? "double" : "single");
    if (_t->precision == MUMPS_DOUBLE) {
	did = phgCalloc(1, sizeof(*did));
	did->par = 1;
	did->sym= solver->symmetry < M_SYM  ? 0 :
		  solver->symmetry == M_SYM ? 2 : 1;
	did->comm_fortran = MPI_Comm_c2f(solver->mat->rmap->comm);
	did->job = -1;
    }
    else {
	sid = phgCalloc(1, sizeof(*sid));
	sid->par = 1;
	sid->sym= solver->symmetry < M_SYM  ? 0 :
		  solver->symmetry == M_SYM ? 2 : 1;
	sid->comm_fortran = MPI_Comm_c2f(solver->mat->rmap->comm);
	sid->job = -1;
    }
    call_mumps(solver);
    _t->factored = FALSE;
    phgInfo(2, "create MUMPS solver.\n");

    return 0;
}

static int
Destroy(SOLVER *solver)
{
    if (solver->oem_data == NULL)
	return 0;

    phgInfo(2, "destroy MUMPS solver.\n");

    if (did != NULL) {
	assert(_t->precision == MUMPS_DOUBLE);
	phgFree(did->irn_loc);
	phgFree(did->jcn_loc);
	phgFree(did->a_loc);
	if (solver->mat->rmap->rank == 0)
	    phgFree(did->rhs);
	did->job = -2;
	call_mumps(solver);
	phgFree(did);
    }
    else if (sid != NULL) {
	assert(_t->precision == MUMPS_SINGLE);
	phgFree(sid->irn_loc);
	phgFree(sid->jcn_loc);
	phgFree(sid->a_loc);
	if (solver->mat->rmap->rank == 0)
	    phgFree(sid->rhs);
	sid->job = -2;
	call_mumps(solver);
	phgFree(sid);
    }
    phgFree(_t->filename);
    phgFree(solver->oem_data);
    solver->oem_data = NULL;

    return 0;
}

static int
Assemble(SOLVER *solver)
{
    MAP *map = solver->mat->rmap;
    int i, j, k, nz, *icntl, *job;
    int Nlocal = map->nlocal, n0 = map->partition[map->rank];
    const MAT_ROW *row;
    double time0;
    FLOAT a;

    if (_t->factored)
	return 0;

    if (solver->sent_to_oem == TRUE)
	phgError(1, "%s:%d: data already sent to OEM solver!\n",
						__FILE__, __LINE__);
    solver->sent_to_oem = TRUE;

    assert(solver->mat->type != PHG_DESTROYED);
    if (solver->mat->type == PHG_MATRIX_FREE && solver->mat->blocks == NULL) {
	phgError(1, "%s:%d: only ordinary or block matrix is allowed.\n",
						__FILE__, __LINE__);
    }

    phgInfo(2, "send data to MUMPS solver.\n");

    /* count # of nonzeros */
    nz = 0;
    if (solver->symmetry < M_SYM) {
	for (i = 0; i < Nlocal; i++)
	    nz += phgMatGetNnzInRow(solver->mat, i, NULL, NULL);
    }
    else {
	for (i = 0; i < Nlocal; i++) {
	    row = phgMatGetRow(solver->mat, i);
	    for (j = 0; j < row->ncols; j++) {
		if (row->cols[j] >= i + n0)
		    nz++;
	    }
	}
    }

    if (_t->precision == MUMPS_DOUBLE) {
	did->nz_loc = nz;
	did->irn_loc = phgAlloc(nz * sizeof(*did->irn_loc));
	did->jcn_loc = phgAlloc(nz * sizeof(*did->jcn_loc));
	did->a_loc = phgAlloc(nz * sizeof(*did->a_loc));

	if (map->rank == 0) {
	    did->n = solver->mat->rmap->nglobal;
	    did->nrhs = 1;
	    did->lrhs = solver->mat->rmap->nglobal;
	    did->rhs = phgAlloc(did->n * sizeof(*did->rhs));
	}
    }
    else {
	sid->nz_loc = nz;
	sid->irn_loc = phgAlloc(nz * sizeof(*sid->irn_loc));
	sid->jcn_loc = phgAlloc(nz * sizeof(*sid->jcn_loc));
	sid->a_loc = phgAlloc(nz * sizeof(*sid->a_loc));

	if (map->rank == 0) {
	    sid->n = solver->mat->rmap->nglobal;
	    sid->nrhs = 1;
	    sid->lrhs = solver->mat->rmap->nglobal;
	    sid->rhs = phgAlloc(sid->n * sizeof(*sid->rhs));
	}
    }

    /* setup MUMPS matrix */
    nz = 0;
    for (i = 0; i < Nlocal; i++) {
	row = phgMatGetRow(solver->mat, i);
	if (_t->precision == MUMPS_DOUBLE) {
	    for (j = 0; j < row->ncols; j++) {
		k = row->cols[j];
		if (solver->symmetry >= M_SYM && k < n0 + i)
		    continue;
		a = row->data[j];
		if (a > DBL_MAX)
		    a = DBL_MAX;
		did->a_loc[nz] = a;
		did->irn_loc[nz] = 1 + n0 + i;
		did->jcn_loc[nz++] = 1 + k;
	    }
	}
	else {
	    for (j = 0; j < row->ncols; j++) {
		k = row->cols[j];
		if (solver->symmetry >= M_SYM && k < n0 + i)
		    continue;
		a = row->data[j];
		if (a > FLT_MAX)
		    a = FLT_MAX;
		sid->a_loc[nz] = a;
		sid->irn_loc[nz] = 1 + n0 + i;
		sid->jcn_loc[nz++] = 1 + k;
	    }
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

    if (solver->mat->refcount == 0)
	phgMatFreeMatrix(solver->mat);

    if (_t->filename != NULL) {
	/* dump matrix to fn.rank.mat */
	static int no = 0;
	char *fn = phgAlloc(strlen(_t->filename) + 32);
	FILE *fp;
	if (no == 0 && map->rank == 0) {
	    sprintf(fn, "/bin/rm -f %s[0-9][0-9].[0-9]*.mat", _t->filename);
	    if (system(fn));	/* Note: avoid warn_unused_result warning */
	}
#if USE_MPI
	MPI_Bcast(&no, 1, MPI_INT, 0, map->comm);
#endif	/* USE_MPI */
	sprintf(fn , "%s%02d.%02d.mat", _t->filename, no++, map->rank);
	phgInfo(1, "saving matrix to \"%s\".\n", fn);
	if ((fp = fopen(fn, "w+t")) == NULL)
	    phgError(1, "cannot open output file \"%s\".\n", fn);
	phgFree(fn);
	fprintf(fp, "# global matrix size, symmetry, nproc\n");
	fprintf(fp, "%"dFMT" %d %d\n", map->nglobal, solver->symmetry,
								map->nprocs);
	fprintf(fp, "# number of local nonzero entries\n");
	if (_t->precision == MUMPS_DOUBLE) {
	    fprintf(fp, "%d\n", did->nz_loc);
	    fprintf(fp, "# entries: i j a_ij\n");
	    for (i = 0; i < did->nz_loc; i++)
		fprintf(fp, "%d %d %0.17lg\n", did->irn_loc[i], did->jcn_loc[i],
				(double)did->a_loc[i]);
	}
	else {
	    fprintf(fp, "%d\n", sid->nz_loc);
	    fprintf(fp, "# entries: i j a_ij\n");
	    for (i = 0; i < sid->nz_loc; i++)
		fprintf(fp, "%d %d %0.17g\n", sid->irn_loc[i], sid->jcn_loc[i],
				(float)sid->a_loc[i]);
	}
	fclose(fp);
    }

    _t->factored = TRUE;

    phgInfo(2, "analyse and factor matrix.\n");

    /* set MUMPS parameters */
    if (_t->precision == MUMPS_DOUBLE) {
	icntl = (&did->ICNTL(1)) - 1;
	job = &did->job;
    }
    else {
	icntl = (&sid->ICNTL(1)) - 1;
	job = &sid->job;
    }

    /* No outputs */
    if (!solver->monitor) {
	icntl[1] = -1;
	icntl[2] = -1;
	icntl[3] = -1;
    }
    icntl[4] = solver->monitor ? 1 : 0;
    icntl[18] = 3;			/* matrix in distributed format */

    /* ordering scheme */
    icntl[7] = ordering_values[_t->ordering];

    /* ICNTL(13)>0 disables use of ScaLAPACK for the root frontal matrix */
    icntl[13] = _t->scalapack_hack;

    /* Whether doing analysis and symbolic factorization in parallel */
    icntl[28] = _t->analysis_type;

    /* Ordering scheme for parallel analysis and symbolic factorization */
    icntl[29] = _t->parallel_ordering;

    /* call MUMPS analyser (JOB 1 = analyse) */
    time0 = phgGetTime(NULL);
    *job = 1;
    call_mumps(solver);
    phgInfo(1, "Time for Analysis: %lg\n", (double)(phgGetTime(NULL) - time0));
    if (_t->fillin >= 0)
	icntl[14] = _t->fillin;	/* extra fill-in percentage */

    /* stopping criterion for refinement */
    if (_t->precision == MUMPS_DOUBLE)
	did->CNTL(2) = 0.;
    else
	sid->CNTL(2) = 0.;
    icntl[10] = _t->refine;	/* number of iterative refinements */
    /* FIXME: add an option to use distributed solution vector, it's said to be
     * much faster than using undistributed solution vector on slow networks */
    icntl[21] = 0;		/* undistributed solution vector */

    time0 = phgGetTime(NULL);
    *job = 2;			/* factorize */
    call_mumps(solver);
    phgInfo(1, "Time for Factorize: %lg\n", (double)(phgGetTime(NULL) - time0));

    /* check return code */
    i = (_t->precision == MUMPS_DOUBLE ? did->INFOG(1) : sid->INFOG(1));
    if (i < 0) {
	if (map->rank == 0) {
	    switch (i) {
		case -8:  case -9:  case -11: case -12:
		case -14: case -15: case -18: case -19:
		    phgInfo(-1, "\n\n");
		    phgInfo(-1, "The return code of MUMPS is %d, which means "
			       "an error has occurred in\n", i);
		    phgInfo(-1, "the solution phase (the value of "
				"'-mumps_fillin' is %d, you may try\n",
				icntl[14]);
		    phgInfo(-1, "rerun the program with a larger value).\n");
		    break;
		case -5: case -7: case -13:
		    phgInfo(-1, "\n\n");
		    phgInfo(-1, "The return code of MUMPS is %d, which means "
			       "a failed alloc (retry with\n", i);
		    phgInfo(-1, "a smaller value of -mumps_fillin?)\n");
		    break;
		default:
		    phgInfo(-1, "\n\n");
		    phgInfo(-1, "The return code of MUMPS is %d\n", i);
		    break;
	    }
	}
	phgError(1, "abort.\n");
    }
    else if (i > 0 && map->rank == 0) {
	phgWarning("MUMPS returns a non-zero code: %d\n", i);
    }

    return 0;
}
 
static int
Solve(SOLVER *solver, VEC *x, BOOLEAN destroy)
{
    MAP *map = solver->mat->rmap;
    int i, Nlocal = map->nlocal;
#if USE_MPI
    float  *sbuffer;
    double *dbuffer;
    int *cnts, *dsps;
#endif	/* USE_MPI */
    double time0;
    FLOAT a;

    Assemble(solver);

    /* gather RHS data */
#if USE_MPI
    cnts = phgAlloc(2 * map->nprocs * sizeof(*cnts));
    dsps = cnts + map->nprocs;
    for (i = 0; i < map->nprocs; i++) {
	cnts[i] = map->partition[i + 1] - map->partition[i];
	dsps[i] = map->partition[i];
    }
    if (_t->precision == MUMPS_DOUBLE) {
#if FT_PHG == FT_DOUBLE
	dbuffer = solver->rhs->data;
#else	/* FT_PHG == FT_DOUBLE */
	dbuffer = phgAlloc(Nlocal * sizeof(*dbuffer));
	for (i = 0; i < Nlocal; i++) {
	    a = solver->rhs->data[i];
	    if (a > DBL_MAX)
		a = DBL_MAX;
	    dbuffer[i] = a;
	}
#endif	/* FT_PHG == FT_DOUBLE */
	MPI_Gatherv(dbuffer, Nlocal, MPI_DOUBLE, did->rhs, cnts, dsps,
			MPI_DOUBLE, 0, map->comm);
	if (dbuffer != (void *)solver->rhs->data)
	    phgFree(dbuffer);
    }
    else {
#if FT_PHG == FT_FLOAT
	sbuffer = solver->rhs->data;
#else	/* FT_PHG == FT_FLOAT */
	sbuffer = phgAlloc(Nlocal * sizeof(*sbuffer));
	for (i = 0; i < Nlocal; i++) {
	    a = solver->rhs->data[i];
	    if (a > FLT_MAX)
		a = FLT_MAX;
	    sbuffer[i] = a;
	}
#endif	/* FT_PHG == FT_FLOAT */
	MPI_Gatherv(sbuffer, Nlocal, MPI_FLOAT, sid->rhs, cnts, dsps,
			MPI_FLOAT, 0, map->comm);
	if (sbuffer != (void *)solver->rhs->data)
	    phgFree(sbuffer);
    }
    phgFree(cnts);
#else	/* USE_MPI */
    if (_t->precision == MUMPS_DOUBLE) {
	for (i = 0; i < Nlocal; i++) {
	    a = solver->rhs->data[i];
	    if (a > DBL_MAX)
		a = DBL_MAX;
	    did->rhs[i] = a;
	}
    }
    else {
	for (i = 0; i < Nlocal; i++) {
	    a = solver->rhs->data[i];
	    if (a > FLT_MAX)
		a = FLT_MAX;
	    sid->rhs[i] = a;
	}
    }
#endif	/* USE_MPI */

    time0 = phgGetTime(NULL);
    if (_t->precision == MUMPS_DOUBLE)
	did->job = 3;		/* solve */
    else
	sid->job = 3;		/* solve */
    call_mumps(solver);
    phgInfo(1, "Time for Solve: %lg\n", (double)(phgGetTime(NULL) - time0));

    if (destroy) {
	if (_t->precision == MUMPS_DOUBLE) {
	    phgFree(did->irn_loc);
	    phgFree(did->jcn_loc);
	    did->irn_loc = did->jcn_loc = NULL;
	    phgFree(did->a_loc);
	    did->a_loc = NULL;
	}
	else {
	    phgFree(sid->irn_loc);
	    phgFree(sid->jcn_loc);
	    sid->irn_loc = sid->jcn_loc = NULL;
	    phgFree(sid->a_loc);
	    sid->a_loc = NULL;
	}
    }

    /* Scatter solution */
#if USE_MPI
    cnts = phgAlloc(2 * map->nprocs * sizeof(*cnts));
    dsps = cnts + map->nprocs;
    if (_t->precision == MUMPS_DOUBLE) {
#if FT_PHG == FT_DOUBLE
	dbuffer = x->data;
#else	/* FT_PHG == FT_DOUBLE */
	dbuffer = phgAlloc(Nlocal * sizeof(*dbuffer));
#endif	/* FT_PHG == FT_DOUBLE */
	for (i = 0; i < map->nprocs; i++) {
	    cnts[i] = map->partition[i + 1] - map->partition[i];
	    dsps[i] = map->partition[i];
	}
	MPI_Scatterv(did->rhs, cnts, dsps, MPI_DOUBLE, dbuffer, Nlocal,
			MPI_DOUBLE, 0, map->comm);
#if FT_PHG != FT_DOUBLE
	for (i = 0; i < Nlocal; i++)
	    x->data[i] = dbuffer[i];
	phgFree(dbuffer);
#endif	/* FT_PHG != FT_DOUBLE */
    }
    else {
#if FT_PHG == FT_FLOAT
	sbuffer = x->data;
#else	/* FT_PHG == FT_FLOAT */
	sbuffer = phgAlloc(Nlocal * sizeof(*sbuffer));
#endif	/* FT_PHG == FT_FLOAT */
	for (i = 0; i < map->nprocs; i++) {
	    cnts[i] = map->partition[i + 1] - map->partition[i];
	    dsps[i] = map->partition[i];
	}
	MPI_Scatterv(sid->rhs, cnts, dsps, MPI_FLOAT, sbuffer, Nlocal,
			MPI_FLOAT, 0, map->comm);
#if FT_PHG != FT_FLOAT
	for (i = 0; i < Nlocal; i++)
	    x->data[i] = sbuffer[i];
	phgFree(sbuffer);
#endif	/* FT_PHG != FT_FLOAT */
    }
    phgFree(cnts);
#else	/* USE_MPI */
    if (_t->precision == MUMPS_DOUBLE) {
	for (i = 0; i < Nlocal; i++)
	    x->data[i] = did->rhs[i];
    }
    else {
	for (i = 0; i < Nlocal; i++)
	    x->data[i] = sid->rhs[i];
    }
#endif	/* USE_MPI */

    if (destroy && map->rank == 0) {
	if (_t->precision == MUMPS_DOUBLE) {
	    phgFree(did->rhs);
	    did->rhs = NULL;
	}
	else {
	    phgFree(sid->rhs);
	    sid->rhs = NULL;
	}
    }

    solver->residual = 0.0;	/* fake residual norm */

    return solver->nits = 1;
}

/*--------------------------------------------------------------------*/

OEM_SOLVER phgSolverMUMPS_ = {
    "MUMPS", RegisterOptions,
    Initialize, Finalize, Init, Create, Destroy,
    AddMatrixEntries, AddRHSEntries,
    Assemble, SetPC, Solve, NULL, NULL, NULL,
    M_UNSYM, FALSE, TRUE, FALSE, FALSE
};

#endif	/* USE_MUMPS */
