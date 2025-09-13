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

/* $Id: coarsen_test1.c,v 1.33 2021/03/11 08:10:48 zlb Exp $ */

#if 0
for f in *.{dat,bz2,mesh} ../{albert,netgen,mesh}-files/{*.bz2,*.dat,*.mesh}; do
    ! test -r $f || ./coarsen_test -mesh_file $f -refines1=16 || break; done
#endif

#include "phg.h"
#include <math.h>

#define log2(x)	(log(x) / log(2.))

static FLOAT
volume(const GRID *g, const ELEMENT *e)
{
    static FLOAT oned6 = 1.0/(FLOAT)6.0;
    FLOAT x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;

    x0 = g->verts[e->verts[0]][0]; 
    y0 = g->verts[e->verts[0]][1];
    z0 = g->verts[e->verts[0]][2];
    x1 = g->verts[e->verts[1]][0] - x0;
    y1 = g->verts[e->verts[1]][1] - y0;
    z1 = g->verts[e->verts[1]][2] - z0;
    x2 = g->verts[e->verts[2]][0] - x0;
    y2 = g->verts[e->verts[2]][1] - y0;
    z2 = g->verts[e->verts[2]][2] - z0;
    x3 = g->verts[e->verts[3]][0] - x0;
    y3 = g->verts[e->verts[3]][1] - y0;
    z3 = g->verts[e->verts[3]][2] - z0;

    return Fabs (x1*y2*z3 + x2*y3*z1 + y1*z2*x3 -
		(z1*y2*x3 + y1*x2*z3 + z2*y3*x1)) * oned6;
}

int
main(int argc, char *argv[])
{
    GRID *g;
    ELEMENT *e;
    INT i, n = 50, periodicity = 0;
    char *fn = "cube.dat";
    FLOAT theta, v, vol, factor0 = 10000.0, rfactor = 512.0;
    FLOAT x, y;
    double t0, t1;

    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterInt("periodicity", "Periodicity", &periodicity);

    phgVerbosity = 0;
    phgInit(&argc, &argv);
    /*phgPause(0);*/

#if 0*USE_MPI
    phgDebug((1, "_phg_comm_count = %d, _phg_comm_count_max = %d\n",
		_phg_comm_count, _phg_comm_count_max));
#endif

    g = phgNewGrid(-1 /*& ~(FACE_FLAG | EDGE_FLAG | GEOM_FLAG)*/);
    phgSetPeriodicity(g, periodicity);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);

    phgPrintf("Mesh file: %s, factor0: %lg, refinement factor: %lg\n", fn,
		    (double)factor0, (double)rfactor);

    /* compute reference volume of the domain */
    vol = 0.0;
    ForAllElements(g, e)
	vol += volume(g, e);
    vol /= factor0;
#if USE_MPI
    x = vol;
    MPI_Allreduce(&x, &vol, 1, PHG_MPI_FLOAT, PHG_SUM, MPI_COMM_WORLD);
#endif

    theta = 2. * 3.1415926535897932 / n;
    for (i = -2; i < n; i++) {
	FLOAT x0 = 0.5 + 0.4 * cos(i * theta), y0 = 0.5 + 0.4 * sin(i * theta);
	phgPrintf("i = %d, theta = %lgPI\n", i, i * 2.0 / n);
	ForAllElements(g, e) {
#if 0
	    FLOAT lambda[] = {0.25, 0.25, 0.25, 0.25}, z;
	    v = volume(g, e);
	    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
	    x -= x0; y -= y0;
	    if (sqrt(x * x + y * y) < 0.2)
		v *= rfactor;
#else
	    int j, k;
	    v = volume(g, e);
	    k = 0;
	    for (j = 0; j < NVert; j++) {
		x = g->verts[e->verts[j]][0] - x0;
		y = g->verts[e->verts[j]][1] - y0;
		k += sqrt(x * x + y * y) < 0.2 ? -1 : 1;
	    }
	    if (k != NVert && k != -NVert)
		v *= rfactor;
#endif
	    e->mark = (v < vol) ? -(int)(log2(vol / v)) :
				   (int)(log2(v / vol) + 0.999999);
	}
	t0 = phgGetTime(NULL);
	phgCoarsenMarkedElements(g);
	t1 = t0;
	phgPrintf("Coarsening: %d elements, %0.4lgs\n", g->nleaf_global,
			(double)((t0 = phgGetTime(NULL)) - t1));
	phgRefineMarkedElements(g);
	t1 = t0;
	phgPrintf("Refinement: %d elements, %0.4lgs\n", g->nleaf_global,
			(double)((t0 = phgGetTime(NULL)) - t1));
	phgBalanceGrid(g, 1.0, 0, NULL, 0.);
	phgPrintf("Repartitioning: %0.4lgs\n", (double)(phgGetTime(NULL) - t0));
	phgCheckConformity(g);
#if 0
	if (i >= 0) {
	    char s[1024];
	    sprintf(s, "output%03d.vtk", i);
	    phgExportVTK(g, s, NULL);
	}
#endif
    }

    phgFreeGrid(&g);

#if 0*USE_MPI
    phgDebug((1, "_phg_comm_count = %d, _phg_comm_count_max = %d\n",
		_phg_comm_count, _phg_comm_count_max));
#endif

    phgFinalize();
    return 0;
}
