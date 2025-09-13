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

/* $Id: coarsen_test.c,v 1.100 2014/06/12 01:53:12 zlb Exp $ */

#if 0
for f in *.{dat,mesh}* ../{albert,netgen,mesh,TrueGrid}-files/{*.dat,*.mesh}*; do ! test -r $f || ! echo ====== $f || mpirun -np 4 ./coarsen_test -mesh_file $f -refines1=16 -refines0=3 -level=1 || break; done
#endif

#include "phg.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static FLOAT
volume(const GRID *g, const ELEMENT *e)
{
    static double oned6 = 1.0/(double)6.0;
    double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;

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

    return fabs (x1*y2*z3 + x2*y3*z1 + y1*z2*x3 -
		(z1*y2*x3 + y1*x2*z3 + z2*y3*x1)) * oned6;
}

static void
f(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    int i, order = DOF_DEFAULT->order;

    for (i = 0; i < DOF_DEFAULT->dim; i++)
	values[i] = Pow(x, order) - i * Pow(y, order) + 3. * Pow(z, order);
}

static void
f1(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
/* (1 2 3) \cross (x y z)  +  (4 5 6) */
{
    *(values++) = 2. * z - 3. * y + 4.;
    *(values++) = 3. * x - 1. * z + 5.;
    *(values++) = 1. * y - 2. * x + 6.;
    return ;
}

int
main(int argc, char *argv[])
{
    GRID *g;
    ELEMENT *e;
    INT i, nroot, nlast, mem_max = 400;
    size_t mem_peak = 0;
    INT refines0 = 3, refines1 = 5;
    char *fn = "tetra.dat";
    double vol0, vol, d;
    INT level = 3, periodicity = 0;

    /* variables for checking DOF interpolation */
    DOF *u, *u1, *tmp;
    FLOAT error, norm;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterInt("-mem_max", "Maximum memory (MB)", &mem_max);
    phgOptionsRegisterInt("-refines0", "refines0", &refines0);
    phgOptionsRegisterInt("-refines1", "refines1", &refines1);
    phgOptionsRegisterInt("-level", "Coarsening level (0 disables coarsening)",
				&level);
    phgOptionsRegisterInt("-periodicity", "Periodicity", &periodicity);

    phgOptionsPreset("-dof_type P1");

    phgVerbosity = 0;
    phgInit(&argc, &argv);
#if 0
    g = phgNewGrid(-1 & ~(FACE_FLAG | EDGE_FLAG | GEOM_FLAG | ELEM_FLAG));
#else
    g = phgNewGrid(-1);
    phgSetPeriodicity(g, periodicity);
#endif
    phgPrintf("Mesh file: %s\n", fn);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    for (i = 0; i < refines0; i++) {
	GRID *g1;
	phgRefineAllElements(g, 1);
	/* flatten the grid */
	g1 = phgDupGrid(g, TRUE);
	phgFreeGrid(&g);
	g = g1;
	phgPrintf("Pre-refine:  %"dFMT" elements, mem: %dMB\n",
			g->nelem_global,
			(int)(phgMemoryUsage(g, &mem_peak) / (1024 * 1024)));
    }

#if 0
    /* test DOF assignment using userfunc */
    u = phgDofNew(g, DOF_P5, 1, "u", f);
    u1 = phgDofNew(g, DOF_ND1, 1, "u1", f1);
#else
    /* test DOF interpolation */
    u = phgDofNew(g, DOF_DEFAULT, 1, "u", DofInterpolation);
    phgDofSetDataByFunction(u, f);
    u1 = phgDofNew(g, DOF_ND1, 1, "u1", DofInterpolation);
    phgDofSetDataByFunction(u1, f1);
#endif

    nroot = g->nroot;

    i = 0;
    do {
	phgBalanceGrid(g, 1.0, -1, NULL, 0.);
	MPI_Bcast(&i, 1, PHG_MPI_INT, 0, g->comm);
	if (i == 0)
	    phgMemoryUsageReset();
	phgMemoryUsage(g, &mem_peak);
	phgPrintf("Balance: nprocs = %d, lif = %lg, mem = %dMB/proc\n",
			g->nprocs, (double)g->lif,
			(int)((mem_peak + 1024*512) / (1024*1024)));

	if (i >= refines1 || mem_peak >= mem_max*(size_t)1024*1024)
	    break;

	/*phgRefineRandomElements(g, "25%");*/
	phgRefineAllElements(g, 1);
	/*phgRefineSurface(g);*/
	phgPrintf("Refine:  %"dFMT" elements\n", g->nelem_global);

	i++;
    } while (TRUE);

    phgBalanceGrid(g, 1.0, -1, NULL, 0.);
    phgMemoryUsage(g, &mem_peak);
    phgPrintf("Balance: nprocs = %d, lif = %lg, mem = %dMB/proc\n",
			g->nprocs, (double)g->lif,
			(int)((mem_peak + 1024*512) / (1024*1024)));

    assert(g->nprocs == phgNProcs);
    phgCheckConformity(g);

    d = 0.0;
    ForAllElements(g, e) {
	d += volume(g, e);
    }

#if USE_MPI
    MPI_Allreduce(&d, &vol0, 1, MPI_DOUBLE, MPI_SUM, g->comm);
    MPI_Bcast(&nroot, sizeof(nroot), MPI_BYTE, 0, g->comm);
#else
    vol0 = d;
#endif

    phgPrintf("After refinements: %"dFMT" elements, vol = %lg\n",
		g->nelem_global, (double)vol0);

    phgPrintf("Coarsening level = %d\n", level);
    nlast = -1;
    if (level > 0) do {
	if (g->nleaf_global == nlast) {
	    if (++level > 3)
		break;
	    phgPrintf("Coarsening level = %d\n", level);
	}
	d = 0.0;
	ForAllElements(g, e) {
	    d += volume(g, e);
	    e->mark = -level;
	}
#if USE_MPI
	MPI_Allreduce(&d, &vol, 1, MPI_DOUBLE, MPI_SUM, g->comm);
#else
	vol = d;
#endif
#if 0
	if (g->bdry_funcs == NULL && fabs(vol - vol0) > 1e-8 * vol0)
	    phgError(1, "volume check: vol=%lg, vol0=%lg, error = %lg\n",
			(double)vol, (double)vol0, (double)(vol - vol0));
#endif
	if (g->nleaf_global == nroot)
	    break;
	nlast = g->nleaf_global;
	phgCoarsenMarkedElements(g);
	phgPrintf("Coarsen: %"dFMT" elements, vol = %lg\n",
				g->nelem_global, (double)vol);
	phgBalanceGrid(g, 1.0, -1, NULL, 0.);
	phgPrintf("Balance: nprocs = %d, lif = %lg\n", g->nprocs, g->lif);
	phgCheckConformity(g);
    } while (TRUE);

    if (level >= 3 && g->nleaf_global != nroot)
	phgError(1, "Program stops before g->nleaf == g->nroot\n");

    /* check geom data interpolation */
    phgGeomCheck(g);

    /* check DOF interpolation */
    if (g->bdry_funcs == NULL) {
	tmp = phgDofNew(g, u->type, u->dim, "tmp", f);
	norm = phgDofNormInftyVec(tmp);
	phgDofAXPY(-1.0, u, &tmp);
	error = phgDofNormInftyVec(tmp);
	phgPrintf("%s: error = %lg (norm = %lg)\n", tmp->type->name,
				error, norm);
	if (error / norm > 1e-8)
	    /*phgWarning("interpolation error with %s!.\n", tmp->type->name);*/
	    phgError(1, "interpolation error with %s!.\n", tmp->type->name);
	phgDofFree(&tmp);
	tmp = phgDofNew(g, u1->type, u1->dim, "tmp", f1);
	norm = phgDofNormInftyVec(tmp);
	phgDofAXPY(-1.0, u1, &tmp);
	error = phgDofNormInftyVec(tmp);
	phgPrintf("%s: error = %lg (norm = %lg)\n", tmp->type->name,
				(double)error, (double)norm);
	if (error / norm > 1e-8)
	    /*phgWarning("interpolation error with %s!.\n", tmp->type->name);*/
	    phgError(1, "interpolation error with %s!.\n", tmp->type->name);
	phgDofFree(&tmp);
    }

    /*phgExportALBERT(g, "output.dat");*/

    phgDofFree(&u);
    phgDofFree(&u1);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
