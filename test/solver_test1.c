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

/* $Id: solver_test.c,v 1.14 2007/09/27 07:37:21 zlb Exp $ */

#include "phg.h"
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int
main(int argc, char *argv[])
{
    VEC *sol;
    MAT *mat;
    SOLVER *solver;
    INT i, row, col, n;
    static INT nglobal;
    static INT *nlocals;
    FILE *f;
    char *fn = "solver_test1.dat", *p, *q;
    static char buffer[1024];

    phgOptionsPreset("-solver gather -verbosity 1");

    phgOptionsRegisterFilename("-input_file", "Input filename", &fn);

    phgInit(&argc, &argv);

    if ((f = fopen(fn, "r")) == NULL)
	phgError(1, "Cannot open \"%s\".\n", fn);

    while (TRUE) {
	if (fgets(buffer, sizeof(buffer), f) == NULL)
	    phgError(1, "unexpected error reading input file.\n");
	p = buffer;
	while (isspace((int)*p))
	    p++;
	if (*p != '#')
	   break;
    }
    nglobal = strtol(p, &q, 10);
    assert(p != q);
    phgPrintf("system size = %d, nprocs = %d\n", nglobal, phgNProcs);

    nlocals = phgAlloc(phgNProcs * sizeof(*nlocals));
    n = 0;
    for (i = 0; i < phgNProcs; i++) {
	p = q;
	nlocals[i] = strtol(p, &q, 10);
	n += nlocals[i];
	assert(p != q);
    }
    assert(n == nglobal);

    mat = phgMatCreate(phgComm, nlocals[phgRank], nglobal);
    sol = phgMapCreateVec(mat->rmap, 1);
    solver = phgMat2Solver(SOLVER_DEFAULT, mat);

    if (mat->rmap->rank == 0) {
	for (i = 0; i < mat->rmap->nglobal; i++)
	    phgSolverAddGlobalMatrixEntry(solver, i, i, 100.0);

	while (TRUE) {
	    if (fgets(buffer, sizeof(buffer), f) == NULL)
		break;
	    p = buffer;
	    while (isspace(*p))
		p++;
	    if (*p == '#' || *p == '\0')
		continue;
	    p = buffer;
	    row = strtol(p, &q, 10);
	    if (p == q)
		break;
	    while (TRUE) {
		p = q;
		col = strtol(p, &q, 10);
		if (p == q)
		    break;
		phgSolverAddGlobalMatrixEntry(solver, row, col, -1.0);
	    }
	}
    }

    fclose(f);

    phgSolverVecSolve(solver, TRUE, sol);

    phgSolverDestroy(&solver);
    phgMatDestroy(&mat);
    phgVecDestroy(&sol);

    phgFinalize();

    return 0;
}
