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

/* $Id: dof_test4.c,v 1.14 2012/09/13 08:12:31 zlb Exp $ */

#include "phg.h"
#include <math.h>
#include <stdlib.h>

static void
func_A(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
	    *values = x;
}
static void
func_A3(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    *(values++) = x;
    *(values++) = y;
    *(values++) = z;
}

static void
f(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
	    *values = -10.0 * (x * x + y * y + z * z);
}
static void
f3(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    *(values++) = 1.;
    *(values++) = 1.;
    *(values++) = 1.;
}

int
main(int argc, char **argv)
{
    GRID *g;
    ELEMENT *e;
    DOF *u_h, *A, *tmp;
    char *fn = "cube.dat";
    int i, j;



    FunctionEntry;

    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&fn);

    phgVerbosity = 0;
    phgInit(&argc, &argv);
    /*phgPause(0); */

    g = phgNewGrid(-1);
    phgImport(g, fn, FALSE);

#if 0
    for (i = 0; i < 4; i++)
	phgRefineAllElements(g, 1);
    phgPartitionGrid(g);
    for (i = 0; i < 8; i++) {
	phgRefineRandomElements(g, "25%");
	phgRedistributeGrid(g);
    }
#endif
    /* test for Lagrange element */
    u_h = phgDofNew(g, DOF_P2, 1, "u", DofInterpolation);
    phgDofSetDataByFunction(u_h, f);
    A = phgDofNew(g, DOF_ANALYTIC, 1, "u", func_A);
    {
	int N = u_h->type->nbas * u_h->dim;	/* size of local matrix */
	FLOAT Mat[N][N];
	phgPrintf("Testing phgQuadGradBasAGradBas ...\n");
	ForAllElements(g, e) {
	    phgPrintf("Element: %p\n", e);
	    for (i = 0; i < N; i++) {
		for (j = 0; j <= i; j++) {
		    Mat[j][i] = Mat[i][j] =
			phgQuadGradBasAGradBas(e, u_h, i,A, u_h, j,
				QUAD_DEFAULT);
		}
	    }
	    for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
		    phgPrintf(" %lg",
			(double)(Fabs(Mat[i][j]) < 1e-15 ? 0. : Mat[i][j]));
		phgPrintf("\n");
	    }
	}
    }
    phgDofFree(&u_h);
    phgDofFree(&A);

    u_h = phgDofNew(g, DOF_P2, 1, "u", DofInterpolation);
    phgDofSetDataByFunction(u_h, f);
    A = phgDofNew(g, DOF_ANALYTIC, 1, "u", func_A);
    {
	int N = u_h->type->nbas * u_h->dim;	/* size of local matrix */
	FLOAT Mat[N][N];
	phgPrintf("Testing phgQuadBasABas ...\n");
	ForAllElements(g, e) {
	    phgPrintf("Element: %p\n", e);
	    for (i = 0; i < N; i++) {
		for (j = 0; j <= i; j++) {
		    Mat[j][i] = Mat[i][j] =
			phgQuadBasABas(e, u_h, i,A, u_h, j, QUAD_DEFAULT);
		}
	    }
	    for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
		    phgPrintf(" %lg",
			(double)(Fabs(Mat[i][j]) < 1e-15 ? 0. : Mat[i][j]));
		phgPrintf("\n");
	    }
	}
    }
    phgDofFree(&u_h);
    phgDofFree(&A);


    u_h = phgDofNew(g, DOF_P2, 1, "u", DofInterpolation);
    phgDofSetDataByFunction(u_h, f);
    tmp = phgDofGradient(u_h, NULL, NULL, "tmp");
    A = phgDofNew(g, DOF_ANALYTIC, 1, "u", func_A);
    {
	int N = u_h->type->nbas * u_h->dim;	/* size of local matrix */
	FLOAT Mat[N];
	phgPrintf("Testing phgQuadDofAGradBas ...\n");
	ForAllElements(g, e) {
	    phgPrintf("Element: %p\n", e);
	    for (i = 0; i < N; i++)
		Mat[i] = phgQuadDofAGradBas(e, tmp, A, u_h, i, QUAD_DEFAULT);
	    for (i = 0; i < N; i++)
		phgPrintf(" %lg",
			(double)(Fabs(Mat[i]) < 1e-15 ? 0. : Mat[i]));
	    phgPrintf("\n");
	}
    }

    phgDofFree(&u_h);
    phgDofFree(&tmp);
    phgDofFree(&A);




    u_h = phgDofNew(g, DOF_ND1, 1, "u", DofInterpolation);
    phgDofSetDataByFunction(u_h, f3);
    A = phgDofNew(g, DOF_ANALYTIC, 3, "u", func_A3);
    {
	int N = u_h->type->nbas * u_h->dim;	
	FLOAT Mat[N][N];
	phgPrintf("Testing phgQuadCurlBasACurlBas...\n");
	ForAllElements(g, e) {
	    phgPrintf("Element: %p\n", e);
	    for (i = 0; i < N; i++) {
		for (j = 0; j <= i; j++) {
		    Mat[j][i] = Mat[i][j] =
			phgQuadCurlBasACurlBas(e, u_h, i, A, u_h, j, -1) -
			phgQuadBasABas(e, u_h, i, A, u_h, j, -1);
		}
	    }
	    for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++)
		    phgPrintf(" %lg",
			(double)(Fabs(Mat[i][j]) < 1e-15 ? 0. : Mat[i][j]));
		phgPrintf("\n");
	    }
	}
    }
    phgDofFree(&u_h);
    phgDofFree(&A);


    phgFreeGrid(&g);
    phgFinalize();

    Return 0;
}
