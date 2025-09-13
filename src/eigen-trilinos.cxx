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

/* $Id: eigen-trilinos.cxx,v 1.21 2018/07/05 06:46:51 zlb Exp $ */

#include "phg.h"

#if USE_TRILINOS_ANASAZI

#ifndef HAVE_CONFIG_H
# define HAVE_CONFIG_H
#endif

#if DEBUG_MPI_COMM
#undef MPI_Comm_dup
#undef MPI_Comm_split
#undef MPI_Comm_create
#undef MPI_Comm_free
#endif	/* DEBUG_MPI_COMM */

#include <Epetra_ConfigDefs.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <AnasaziEpetraAdapter.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <Teuchos_RefCountPtr.hpp>

#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziBlockDavidsonSolMgr.hpp>
#include <AnasaziLOBPCGSolMgr.hpp>

typedef double ST;
typedef Epetra_MultiVector MV;
typedef Epetra_Operator OP;

static Epetra_CrsMatrix *
create_mat(Epetra_Map map, MAT *A, int work[])
{
    int i, n0, nmax, nlocal;
    MAT_ROW *row;
    Epetra_CrsMatrix *a;
    int j, k, *cols;

    if (A == NULL)
	return NULL;

    if (A->type == PHG_PACKED)
	phgMatUnpack(A);

    assert(A->type == PHG_UNPACKED);

    nlocal = A->rmap->nlocal;
    n0 = A->rmap->partition[A->rmap->rank];
    assert(sizeof(ST) == sizeof(FLOAT) && sizeof(int) == sizeof(INT));
    nmax = 0;
    for (i = 0, row = A->rows; i < nlocal; i++, row++) {
	if (nmax < row->ncols)
	    nmax = row->ncols;
        work[i] = row->ncols;
    }
    cols = (int *)phgAlloc(nmax * sizeof(*cols));
    a = new Epetra_CrsMatrix(Copy, map, work);
    for (i = 0, row = A->rows; i < nlocal; i++, row++) {
	if (row->ncols == 0)
	    continue;
	for (j = 0; j < row->ncols; j++) {
	    k = row->cols[j];
	    cols[j] = (k < nlocal ? k + n0 : A->O2Gmap[k - nlocal]);
	}
	a->InsertGlobalValues(n0 + i, row->ncols,
				(double *)row->data, (int *)cols);
#if 0
	if (A->refcount == 0) {
	    phgFree(row->data);
	    row->data = NULL;
	    phgFree(row->cols);
	    row->cols = NULL;
	    row->ncols = 0;
	}
#endif
    }
#if 0
    if (solver->mat->refcount == 0)
	phgMatFreeMatrix(A);
#endif
    a->FillComplete();
    phgFree(cols);

    return a;
}

EigenSolverPrototype(phgEigenSolverTrilinos)
/* (A, B, n, which, tau, evals, evecs, eps, itmax, nit) */
{
    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	/* register options */
	initialized = TRUE;
	/*phgOptionsRegisterTitle("\nTrilinos Anasazi options:", NULL, "trilinos");*/
	return 0;
    }

    if (n > A->rmap->nglobal)
	n = A->rmap->nglobal;

    if (n <= 0 || A->rmap->nglobal <= 0)
	return 0;

    INT i, j, n0;

    Epetra_MpiComm comm(A->rmap->comm);

    n0 = A->rmap->partition[A->rmap->rank];
    VEC *vec = *evecs;
    MagicCheck(VEC, vec);

    // map
    int *work = new int[A->rmap->nlocal];
    for (i = 0, j = n0; i < A->rmap->nlocal; i++)
	work[i] = j++;
    Epetra_Map map(A->rmap->nglobal, A->rmap->nlocal, work, 0, comm);

    // matrices
    Epetra_CrsMatrix *EA = create_mat(map, A, work);
    Epetra_CrsMatrix *EB = create_mat(map, B, work);
    delete [] work;

    int nev = (n <= A->rmap->nglobal ? n : A->rmap->nglobal);
    int blocksize = nev;
    Teuchos::RefCountPtr<Epetra_MultiVector> X = 
	Teuchos::rcp(new Epetra_MultiVector(map, blocksize));
    if (vec == NULL) {
	*evecs = vec = phgMapCreateVec(A->cmap, nev);
	X->Random();
    }
    else if (vec->nvec < nev) {
	vec->nvec = nev;
	phgFree(vec->data);
	vec->data = (FLOAT *)phgAlloc(nev * A->rmap->nlocal *
						sizeof(*vec->data));
	X->Random();
    } else {
	vec->nvec = nev;
	/* TODO: set up X using vec->data values */
	X->Random();
    }

    Teuchos::RefCountPtr<Epetra_CrsMatrix> tEA = Teuchos::rcp(EA, false);
    Teuchos::RefCountPtr<Epetra_CrsMatrix> tEB = Teuchos::rcp(EB, false);
    Teuchos::RefCountPtr< Anasazi::BasicEigenproblem<ST,MV,OP> > problem =
        Teuchos::rcp(new Anasazi::BasicEigenproblem<ST,MV,OP>(tEA, tEB, X));
    problem->setHermitian(true);
    problem->setNEV(n);
    bool ret = problem->setProblem();
    if (ret != true)
	phgError(1, "%s: setProblem() failed.\n", __func__);

    Teuchos::ParameterList mypl;
    mypl.set( "Which", which <= 0 ? "SM" : "LM");
#if 0
    mypl.set( "Verbosity", Anasazi::Errors + Anasazi::Warnings +
                           Anasazi::TimingDetails + Anasazi::FinalSummary +
                           Anasazi::IterationDetails);
#endif
    if (blocksize > A->rmap->nglobal / 2)
	blocksize = A->rmap->nglobal / 2;
    mypl.set( "Block Size", blocksize);
    mypl.set( "Num Blocks", /*A->rmap->nglobal / blocksize*/2);
    mypl.set( "Maximum Restarts", itmax * n);
    mypl.set( "Maximum Iterations", itmax * n);
    mypl.set( "Convergence Tolerance", eps);

    //Anasazi::BlockDavidsonSolMgr<ST,MV,OP> solver(problem, mypl);
    //Anasazi::BlockKrylovSchurSolMgr<ST,MV,OP> solver(problem, mypl);
    Anasazi::LOBPCGSolMgr<double,MV,OP> solver(problem, mypl);

    Anasazi::ReturnType sret = solver.solve();

    delete EA;
    if (EB != NULL)
	delete EB;

//cerr << A->rmap->rank << ": " << X << endl;
//cerr << A->rmap->rank << ": " << problem->getSolution().Evecs << endl;
    Anasazi::Eigensolution<ST, Epetra_MultiVector> sol = problem->getSolution();
    if (phgVerbosity > 0)
	phgPrintf("Anasazi: converged = %d, numVecs = %d\n", sret, sol.numVecs);
    if (sret != Anasazi::Converged)
	nev = sol.numVecs;
    else
	assert(nev == sol.numVecs);

    if (nev > 0) {
	Teuchos::RefCountPtr<Epetra_MultiVector> vecs = sol.Evecs;
	std::vector<Anasazi::Value<ST> > vals = sol.Evals;
	double *a, *b;
	int lda;
	vecs->ExtractView(&a, &lda);
	assert(FT_PHG == FT_DOUBLE);
	b = (double *)vec->data;
	for (i = 0; i < nev; i++) {
	    evals[i] = vals[i].realpart;
	    for (j = 0; j < A->rmap->nlocal; j++)
		*(b++) = *(a++);
	    a += lda - A->rmap->nlocal;
	}
    }
    else {
	phgFree(vec->data);
	vec->data = NULL;
    }

    vec->nvec = nev;
    *nit = 1;

    return nev;
}

#else	/* USE_TRILINOS_ANASAZI */
EigenSolverVoid(phgEigenSolverTrilinos)
#endif	/* USE_TRILINOS_ANASAZI */
