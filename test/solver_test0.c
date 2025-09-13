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

/* $Id: solver_test0.c,v 1.5 2022/08/30 06:03:01 zlb Exp $ */

#include "phg.h"
#include <string.h>
#include <math.h>

int
main(int argc, char *argv[])
{
    VEC *sol;
    MAT *mat;
    SOLVER *solver;
    INT i, N0, N = 3 * 4;

    phgOptionsPreset("-solver gmres -solver_monitor");

    phgVerbosity = 0;
    phgInit(&argc, &argv);

    assert(N % phgNProcs == 0);
    mat = phgMatCreate(phgComm, N / phgNProcs, N);
    sol = phgMapCreateVec(mat->rmap, 1);
    solver = phgMat2Solver(SOLVER_DEFAULT, mat);

    /* construct matrix:
     * 		1 1 1 ... ... 1 1 1
     * 		0 1 0 ... ... 0 0 0
     * 		0 0 1 ... ... 0 0 0
     * 		...................
     * 		0 0 0 ... ... 1 0 0
     * 		0 0 0 ... ... 0 1 0
     * 		0 0 0 ... ... 0 0 1
     */

    /* row 0 */
    if (mat->rmap->rank == 0) {
	for (i = 0; i < mat->rmap->nglobal; i++)
	    phgSolverAddGlobalMatrixEntry(solver, 0, i, 1.0);
	solver->rhs->data[0] = (FLOAT)mat->rmap->nglobal;
    }

    /* rows 1..N-1 */
    for (i = (mat->rmap->rank == 0 ? 1 : 0); i < mat->rmap->nlocal; i++) {
	phgSolverAddMatrixEntry(solver, i, i, 1.0);
	solver->rhs->data[i] = 1.0;
    }

    phgSolverVecSolve(solver, TRUE, sol);

    N0 = mat->rmap->partition[mat->rmap->rank];
    for (i = 0; i < mat->rmap->nlocal; i++)
	phgInfo(-1, "error[%d] = %e\n", i + N0, (double)(sol->data[i] - 1.));

    phgSolverDestroy(&solver);
    phgMatDestroy(&mat);
    phgVecDestroy(&sol);

    phgFinalize();

    return 0;
}
