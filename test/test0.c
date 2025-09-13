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

/* $Id: test0.c,v 1.46 2012/09/13 08:12:31 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <string.h>

int
main(int argc, char *argv[])
{
    GRID *g;
    char *fn = "cube4.dat";
    INT depth = 1, nelem0, nelem1;
    INT periodicity = X_MASK | Y_MASK | Z_MASK;

    phgOptionsRegisterInt("depth", "Refinement depth", &depth);    
    phgOptionsRegisterInt("periodicity", "Periodicity", &periodicity);    
    phgOptionsRegisterFilename("mesh_file", "Mesh file", &fn);    
    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    phgSetPeriodicity(g, periodicity);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);

    phgPrintf("Initial mesh: vert=%d (%d), edge=%d, face=%d, elem = %d\n",
		g->period == NULL ? g->nvert_global : g->period->nvert_global,
		g->nvert_global, g->nedge_global, g->nface_global,
		g->nelem_global);

    nelem0 = g->nelem_global;

    phgExportVTK(g, "test0.vtk", NULL);
    phgExportDX(g, "test0.dx", NULL);
    phgExportMedit(g, "test0.mesh");

    /* test refinement */
    while (g->nelem_global < 20000) {
	phgRefineAllElements(g, depth);
	phgCheckConformity(g);
	phgBalanceGrid(g, 1.0, 1, NULL, 0.0);
	phgPrintf("np=%d, lif=%4.2f, vert=%d (%d), edge=%d, face=%d, elem=%d\n",
		  g->nprocs, (double)g->lif,
		  g->period == NULL ? g->nvert_global : g->period->nvert_global,
		  g->nvert_global, g->nedge_global, g->nface_global,
		  g->nelem_global);
    }

    /* test coarsening */
    nelem1 = -1;
    while (g->nelem_global != nelem1 && g->nelem_global != nelem0) {
	ELEMENT *e;
	nelem1 = g->nelem_global;
	ForAllElements(g, e)
	    e->mark = -depth;
	phgCoarsenMarkedElements(g);
	phgCheckConformity(g);
	phgBalanceGrid(g, 1.0, 1, NULL, 0.0);
	phgPrintf("np=%d, lif=%4.2f, vert=%d (%d), edge=%d, face=%d, elem=%d\n",
		  g->nprocs, (double)g->lif,
		  g->period == NULL ? g->nvert_global : g->period->nvert_global,
		  g->nvert_global, g->nedge_global, g->nface_global,
		  g->nelem_global);
    }

    //phgDumpGridInfo(g);

    phgFreeGrid(&g);
    phgFinalize();
    
    return 0;
}
