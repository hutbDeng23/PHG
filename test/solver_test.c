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

/* $Id: solver_test.c,v 1.16 2020/12/06 09:03:08 zlb Exp $ */

#include "phg.h"
#include <string.h>
#include <math.h>

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *u, *v, *w, *u0, *v0, *w0;
    SOLVER *solver;
    INT i, j;
    char *fn = "cube.dat";

    phgVerbosity = 0;
    phgInit(&argc, &argv);
    /*phgPause(0);*/
    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    for (i = 0; i < 4; i++)
	phgRefineAllElements(g, 1);
    phgBalanceGrid(g, 1.0, 0, NULL, 0.);
    for (i = 0; i < 8; i++) {
	phgRefineRandomElements(g, "25%");
	phgBalanceGrid(g, 1.2, -1, NULL, 0.);
    }
    phgCheckConformity(g);

    u = phgDofNew(g, DOF_P1,  1, "u", DofNoAction);
    v = phgDofNew(g, DOF_ND1, 1, "v", DofNoAction);
    w = phgDofNew(g, DOF_P1,  1, "w", DofNoAction);

    phgDofSetDataByValue(u, 0.0);
    phgDofSetDataByValue(v, 0.0);
    phgDofSetDataByValue(w, 0.0);

    u0 = phgDofNew(g, u->type, u->dim, "u0", DofNoAction);
    v0 = phgDofNew(g, v->type, v->dim, "v0", DofNoAction);
    w0 = phgDofNew(g, w->type, w->dim, "w0", DofNoAction);

    solver = phgSolverCreate(SOLVER_DEFAULT, u, v, w, NULL);

    for (i = 0; i < DofDataCount(u); i++) {
	j = phgSolverMapD2L(solver, 0, i);
	u0->data[i] = (FLOAT)(phgSolverMapL2G(solver, j) % 3);
	phgSolverAddMatrixEntry(solver, j, j, 1.0);
	phgSolverAddRHSEntry(solver, j, u0->data[i]);
    }

    for (i = 0; i < DofDataCount(v); i++) {
	j = phgSolverMapD2L(solver, 1, i);
	v0->data[i] = (FLOAT)(phgSolverMapL2G(solver, j) % 5);
	phgSolverAddMatrixEntry(solver, j, j, 1.0);
	phgSolverAddRHSEntry(solver, j, v0->data[i]);
    }

    for (i = 0; i < DofDataCount(w); i++) {
	j = phgSolverMapD2L(solver, 2, i);
	w0->data[i] = (FLOAT)(phgSolverMapL2G(solver, j) % 7);
	phgSolverAddMatrixEntry(solver, j, j, 1.0);
	phgSolverAddRHSEntry(solver, j, w0->data[i]);
    }

    phgSolverSolve(solver, TRUE, u, v, w, NULL);
    phgSolverDestroy(&solver);

    phgDofAXPY(-1.0, u0, &u);
    phgDofAXPY(-1.0, v0, &v);
    phgDofAXPY(-1.0, w0, &w);
    phgPrintf("Error: u = %lg, v = %lg, w = %lg\n",
		(double)phgDofNormInftyVec(u),
		(double)phgDofNormInftyVec(v),
		(double)phgDofNormInftyVec(w));

    phgDofFree(&u);
    phgDofFree(&v);
    phgDofFree(&w);
    phgDofFree(&u0);
    phgDofFree(&v0);
    phgDofFree(&w0);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
