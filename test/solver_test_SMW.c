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

/* This program tests solving 'sparse + low rank' systems using
 * Sherman-Morrison-Woodbury formula:
 *
 *	inv(A + U*V^t) = (I - inv(A)*U*C*V^t) * inv(A),
 *	C = inv(I + V^t*inv(A)*U)
 *
 * $Id: solver_test_SMW.c,v 1.19 2021/03/11 08:10:48 zlb Exp $ */

#include "phg.h"
#include <string.h>
#include <math.h>

/* PC_PROC function which applies the Sherman-Morrison-Woodbury formula */
static void
sherman(void *ctx, VEC *b0, VEC **x0)
{
    SOLVER *solver = ctx;
    VEC *U = ((void **)solver->mat->mv_data)[1];
    SOLVER *solver1 = ((void **)solver->mat->mv_data)[2];
    FLOAT *C = ((void **)solver->mat->mv_data)[3];
    INT *pvt = ((void **)solver->mat->mv_data)[4];

    VEC *x = *x0, *y;
    INT i, j, n = U->map->nlocal, K = U->nvec;
    FLOAT *z;

    phgPrintf("Note: applying the Sherman-Morrison-Woodbury formula.\n");

    /* x := inv(A) * b */
#if 0
    memset(x->data, 0, n * sizeof(*x->data));
#endif
    z = solver1->rhs->data; 	/* save solver1->rhs->data */
    solver1->rhs->data = b0->data;
    solver1->rhs->assembled = TRUE;
    phgSolverVecSolve(solver1, FALSE, x);
    solver1->rhs->data = z; 	/* restore solver1->rhs->data */

    if (K == 0)
	return;

    /* z := U^t * x */
    z = phgCalloc(K, sizeof(*z));
    for (j = 0; j < K; j++)
	for (i = 0; i < n; i++)
	    z[j] += U->data[j * n + i] * x->data[i];
#if USE_MPI
    if (U->map->nprocs > 1) {
	FLOAT *tmp = phgAlloc(K * sizeof(*tmp));
	MPI_Allreduce(z, tmp, K, PHG_MPI_FLOAT, PHG_SUM, U->map->comm);
	phgFree(z);
	z = tmp;
    }
#endif	/* USE_MPI */

    /* z := C*z */
    phgSolverDenseSV(K, C, pvt, 1, z);

    /* y := inv(A) * (U*z) */
    memset(solver1->rhs->data, 0, n * sizeof(*solver1->rhs->data));
    for (i = 0; i < n; i++)
	for (j = 0; j < K; j++)
	    solver1->rhs->data[i] += U->data[j * n + i] * z[j];
    y = phgMapCreateVec(U->map, 1);
    phgSolverVecSolve(solver1, FALSE, y);
    phgFree(z);

    /* x -= y */
    for (i = 0; i < n; i++)
	x->data[i] -= y->data[i];

    phgVecDestroy(&y);
}

/* callback function which computes A+UU^t */
static int
funcB(MAT_OP op, MAT *B, VEC *x, VEC *y)
/* computes y = op(A) * x */
{
    MAT *A = ((void **)B->mv_data)[0];
    VEC *U = ((void **)B->mv_data)[1];
    INT i, j, k, n = U->map->nlocal, K = U->nvec;

    /* compute y = op(A)*x */
    phgMatVec(op, 1.0, A, x, 0.0, &y);

    if (K == 0)
	return 0;

    if (op != MAT_OP_D) {
	/* compute y += U*U^t*x */
	FLOAT *z = phgCalloc(K, sizeof(*z));
	/* z = U^t*x */
	for (j = 0; j < K; j++)
	    for (i = 0; i < n; i++)
		z[j] += U->data[j * n + i] * x->data[i];
#if USE_MPI
	if (U->map->nprocs > 1) {
	    FLOAT *tmp = phgCalloc(K, sizeof(*tmp));
	    MPI_Allreduce(z, tmp, K, PHG_MPI_FLOAT, PHG_SUM, U->map->comm);
	    phgFree(z);
	    z = tmp;
	}
#endif	/* USE_MPI */
	/* y += U*z */
	for (i = 0; i < n; i++)
	    for (j = 0; j < K; j++)
		y->data[i] += U->data[j * n + i] * z[j];
	phgFree(z);
	return 0;
    }

    for (i = 0; i < n; i++) {
	FLOAT d = 0.0;
	for (j = 0; j < K; j++)
	    for (k = 0; k < K; k++)
		d += U->data[j * n + i] * U->data[k * n + i];
	y->data[i] += d * x->data[i];
    }

    return 0;
}

int
main(int argc, char *argv[])
{
    MAT *A, *B;
    VEC *U = NULL, *x;
    FLOAT *C;
    SOLVER *solver, *solver1 = NULL;
    INT i, j, k, n, *pvt, N = 1000, K = 2;
    char *main_opts = NULL, *sub_opts = NULL;

    phgOptionsRegisterInt("-n", "N value", &N);
    phgOptionsRegisterInt("-k", "K value", &K);
    phgOptionsRegisterString("-main_solver_opts", "Options for the main solver",
				&main_opts);
    phgOptionsRegisterString("-sub_solver_opts", "Options for the subsolver",
				&sub_opts);

    /* a direct solver is preferable for the sparse matrix */
    phgOptionsPreset("-solver mumps");
    phgInit(&argc, &argv);

    phgPrintf(
"----------------------------------------------------------------------------\n"
"This code solves (A+UU^t)x=b using the Sherman-Morrison-Woodbury formula.\n"
"Note: may use the following to disable use of the Sherman-Morrison-Woodbury\n"
"algorithm and change to the default solver instead:\n"
"	-preonly_pc_type solver -preonly_pc_opts \"-solver_maxit 2000\"\n"
"----------------------------------------------------------------------------\n"
    );

   phgPrintf("Generating the linear system: N = %"dFMT", K = %"dFMT"\n", N, K);

    /* A is a distributed NxN SPD tridiagonal matrix (A = [-1, 2, -1]) */
    n = N / phgNProcs + (phgRank < (N % phgNProcs) ? 1 : 0);
    A = phgMatCreate(phgComm, n, N);
    phgPrintf("  Generating matrix A.\n");
    for (i = 0; i < n; i++) {
	/* diagonal */
	phgMatAddEntry(A, i, i, 2.0);
	/* diagonal - 1 */
	if (i > 0)
	    phgMatAddEntry(A, i, i - 1, -1.0);
	else if (phgRank > 0)
	    phgMatAddLGEntry(A, i, A->rmap->partition[phgRank] - 1, -1.0);
	/* diagonal + 1 */
	if (i < n - 1)
	    phgMatAddEntry(A, i, i + 1, -1.0);
	else if (phgRank < phgNProcs - 1)
	    phgMatAddLGEntry(A, i, A->rmap->partition[phgRank] + n, -1.0);
    }
    phgMatAssemble(A);

    /* U is a K-component vector */
    U = phgMapCreateVec(A->rmap, K);
    phgVecRandomize(U, 123);

    /* solver1 is the solver for A */
    phgOptionsPush();
    phgOptionsSetOptions(sub_opts);
    solver1 = phgMat2Solver(SOLVER_DEFAULT, A);
    phgOptionsPop();

    /* x is a scratch vector */
    x = phgMapCreateVec(A->rmap, 1);

    /* C is a KxK dense matrix, pvt is an integer array, they store the LU
     * factorization of (I + U^t*inv(A)*U) */
    phgPrintf("  Generating the dense matrix I+U^t*inv(A)*U.\n");
    C = phgCalloc(K * K, sizeof(*C));
    pvt = phgAlloc(K * sizeof(*pvt));
    for (i = 0; i < K; i++) {
	for (j = 0; j < n; j++) {
	    solver1->rhs->data[j] = U->data[i * n + j];
	    x->data[j] = 0.0;
	}
	solver1->rhs->assembled = TRUE;
	phgSolverVecSolve(solver1, FALSE, x);
	for (j = 0; j < K; j++)
	    for (k = 0; k < n; k++)
		C[i * K + j] += U->data[j * n + k] * x->data[k];
    }
#if USE_MPI
    if (U->map->nprocs > 1) {
	FLOAT *tmp = phgAlloc(K * K * sizeof(*tmp));
	MPI_Allreduce(C, tmp, K * K, PHG_MPI_FLOAT, PHG_SUM, U->map->comm);
	phgFree(C);
	C = tmp;
    }
#endif	/* USE_MPI */
    for (i = 0; i < K; i++)
	C[i * K + i] += 1.0;
    phgPrintf("  Factorizing the dense matrix I+U^t*inv(A)*U.\n");
    phgSolverDenseLU(K, C, pvt);

    /* B is a matrix-free matrix representing A + U*U^t, B->mv_data is used
     * to pass A, U, solver1, C and pvt to callback functions */
    B = phgMapCreateMatrixFreeMat(A->rmap, A->cmap, funcB,
				  /* arguments carried over to CB functions */
				  A, U, solver1, C, pvt, NULL);

    /* solver is a PreOnly solver for B whose pc_proc is set to sherman().
     * 
     * Note: can also use pcg, gmres, or petsc for this solver, in this case
     * the solution obtained with the Sherman-Morisson formula is iteratively
     * refined. */
    phgOptionsPush();
    phgOptionsSetOptions("-solver preonly");
    phgOptionsSetOptions(main_opts);
    solver = phgMat2Solver(SOLVER_DEFAULT, B);
    phgSolverSetPC(solver, solver, sherman);
    phgOptionsPop();

    for (i = 0; i < n; i++)
	x->data[i] = 1.0;
    phgMatVec(MAT_OP_N, 1.0, B, x, 0.0, &solver->rhs);

    phgPrintf("Solving the linear system.\n");

    /* reset initial solution to zero */
    memset(x->data, 0, n * sizeof(*x->data));
    phgSolverVecSolve(solver, TRUE, x);
    for (i = 0; i < n; i++)
	solver->rhs->data[i] = 1.0;
    phgVecAXPBY(-1.0, solver->rhs, 1.0, &x);
    phgPrintf("Checking the result: |x - x_exact| / |x_exact| = %lg\n",
			(double)phgVecNorm2(x, 0, NULL) / sqrt((double)N));

    phgSolverDestroy(&solver);
    phgSolverDestroy(&solver1);
    phgMatDestroy(&A);
    phgMatDestroy(&B);
    phgVecDestroy(&U);
    phgVecDestroy(&x);
    phgFree(C);
    phgFree(pvt);

    phgFinalize();

    return 0;
}
