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

/* $Id: mat_test1.c,v 1.8 2012/06/21 13:50:51 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void
do_test(MAP *rmap, MAP *cmap)
{
    MAT *A;
    VEC *x, *y;
    INT i, j, m0, n0;
    INT *cols = phgAlloc(cmap->nglobal * sizeof(*cols));
    FLOAT *data = phgAlloc(cmap->nglobal * sizeof(*data));

    A = phgMapCreateMat(rmap, cmap);
    m0 = A->rmap->partition[A->rmap->rank];
    n0 = A->cmap->partition[A->cmap->rank];

    /* Matrix entries: A(I,J) = 1 + (I-1) + (J-1) */
    for (i = m0; i < m0 + A->rmap->nlocal; i++) {
#if 0
	/* Test MatAddEntry */
	for (j = 0; j < A->cmap->nglobal; j++)
	    phgMatAddGlobalEntry(A, i, j, 1.0 + i + j);
#else
	/* Test MatAddEntries */
	for (j = 0; j < A->cmap->nglobal; j++) {
	    cols[j] = j;
	    data[j] = 1.0 + i + j;
	}
	phgMatAddGlobalEntries(A, 1, &i, A->cmap->nglobal, cols, data);
#endif
    }
    phgFree(cols);
    phgFree(data);
    phgMatAssemble(A);

    phgInfo(-1, "y = A * x\n");
    x = phgMapCreateVec(A->cmap, 1);
    for (i = 0; i < x->map->nlocal; i++)
	x->data[i] = 1.0 + i + n0;
    phgVecAssemble(x);
    y = phgMatVec(MAT_OP_N, 1.0, A, x, 0.0, NULL);
    for (i = 0; i < y->map->nlocal; i++)
	phgInfo(-1, "    y->data[%d] = %lg\n", i + m0, (double)y->data[i]);
    phgVecDestroy(&x);
    phgVecDestroy(&y);

    phgInfo(-1, "y = A' * x\n");
    x = phgMapCreateVec(A->rmap, 1);
    for (i = 0; i < x->map->nlocal; i++)
	x->data[i] = 1.0 + i + m0;
    phgVecAssemble(x);
    y = phgMatVec(MAT_OP_T, 1.0, A, x, 0.0, NULL);
    for (i = 0; i < y->map->nlocal; i++)
	phgInfo(-1, "    y->data[%d] = %lg\n", i + n0, (double)y->data[i]);
    phgVecDestroy(&x);
    phgVecDestroy(&y);

    phgMatDestroy(&A);
}

int
main(int argc, char *argv[])
{
    GRID *g;
    MAP *map, *map0;
    INT n = 4;
    DOF *u;
    static INT pre_refines = 0;
    static char *fn = "../test/cube.dat";
    /*static char *fn = "../albert-files/cube5.dat";*/

    phgOptionsPreset("-dof_type P1");
    phgOptionsRegisterInt("-pre_refines", "Pre-refinements", &pre_refines);
    phgInit(&argc, &argv);
    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgRefineAllElements(g, pre_refines);
    phgPartitionGrid(g);

    u = phgDofNew(g, DOF_DEFAULT, 1, "u", DofNoAction);
    map0 = phgMapCreate(u, NULL);

    phgInfo(-1, "\n");
    phgInfo(-1, "Test special map with size = 0 on all but one process:\n");
    map = phgMapCreateSimpleMap(g->comm, g->rank == g->nprocs - 1 ? n : 0, n);
    phgInfo(-1, "---- DOF rows, special columns:\n");
    do_test(map0, map);
    phgInfo(-1, "---- Special rows, DOF columns:\n");
    do_test(map, map0);
    phgInfo(-1, "---- Special rows, special columns:\n");
    do_test(map, map);
    phgMapDestroy(&map);

    phgInfo(-1, "\n");
    phgInfo(-1, "Test serial map:\n");
    map = phgMapCreateSimpleMap(g->comm, n, n);
    phgInfo(-1, "---- DOF rows, serial columns:\n");
    do_test(map0, map);
    phgInfo(-1, "---- serial rows, DOF columns:\n");
    do_test(map, map0);
    phgInfo(-1, "---- Serial rows, serial columns:\n");
    do_test(map, map);
    phgMapDestroy(&map);

    phgMapDestroy(&map0);
    phgDofFree(&u);
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
