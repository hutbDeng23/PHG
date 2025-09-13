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

/* $Id: phg_test.c,v 1.25 2012/09/13 08:12:31 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <math.h>

static const char *refine_types[] = {"uniform", "random", NULL};
static int refine_type = 0;
static char *random_percent = "25%";

static const char *output_types[] = {"none", "vtk", "opendx", NULL};
static int output_type = 0;
static char *output_file = "output";

int
main(int argc, char *argv[])
{
    int i, n, m;
    GRID *g;
    DOF *u, *v;
    ELEMENT *e;

    FunctionEntry;

    phgOptionsRegisterKeyword("refine_type", "Refine type", refine_types,
				&refine_type);
    phgOptionsRegisterString("random_percent", "Random refine percent",
				&random_percent);
    phgOptionsRegisterFilename("output_file", "Output filename", &output_file);
    phgOptionsRegisterKeyword("output_type", "Output type", output_types,
				&output_type);

    phgVerbosity = 0;
    phgInit(&argc, &argv);
    phgOptionsShowUsed();

    if (argc != 2 && argc != 3 && argc != 4) {
	fprintf(stderr, "Usage:\n  phg_test [options] mesh_file "
			"[serial_refines [parallel_refines]] [options]\n");
	fprintf(stderr, "Example:\n  %s cube.dat 2 3\n", argv[0]);
	fprintf(stderr, "(2 serial refines, 3 parallel refines)\n");
	fprintf(stderr, "Type \"phg_test -help\" to show a list of options.\n");
	exit(1);
    }

#if 0
    phgPause(0);
#endif

    g = phgNewGrid(-1);
    if (!phgImport(g, argv[1], FALSE))
	phgError(1, "can't read file \"%s\".\n", argv[1]);
    phgCheckConformity(g);

#if 0
/* Dump selected element for debugging */
#define GEOM_TYPE_TO_ALBERT(n) \
        ((n) == DIAGONAL ? 0 : ((n) == FACE ? 1 : ((n) == EDGE ? 2 : \
        ((n) == MIXED ? 3 : 4))))
COORD *c;
/* The index can be found with ParaView through
	Mesh quality ==> Cell data to point data ==> Threshold => ancestor */
e = g->roots + 94876;
phgPrintf("ALBERT type: %d\n", GEOM_TYPE_TO_ALBERT(e->type));
for (i = 0; i < 4; i++) {
    c = g->verts + e->verts[i];
    phgPrintf("%lf %lf %lf\n", (double)(*c)[0], (double)(*c)[1], (double)(*c)[2]);
}
exit(0);
#endif

    /* u will carry the index (in g->roots[]) of the ancestor before the
     * mesh is distributed */
    u = phgDofNew(g, DOF_P0, 1, "ancestor", DofInterpolation);
    for (i = 0; i < g->nroot; i++)
	*DofElementData(u, g->roots[i].index) = (FLOAT)i;

    /* print some statistics on the initial mesh */
    phgPrintf("********* Initial mesh ***********\n");
    phgPrintf("Number of elements: %d\n", g->nleaf_global);
    phgPrintf("Number of vertices: %d\n", g->nvert_global);
    phgPrintf("Number of edges:    %d\n", g->nedge_global);
    phgPrintf("Number of faces:    %d\n", g->nface_global);
    phgPrintf("**********************************\n");

    n = (argc >= 3) ? atoi(argv[2]) : 0;
    m = (argc >= 4) ? atoi(argv[3]) : 0;
    phgPrintf("Serial refines: %d, parallel refines: %d\n", n, m);

    for (i = 0; i < n; i++) {
	if (refine_type == 0) {
	    phgRefineAllElements(g, 1);
	}
	else {
	    phgRefineRandomElements(g, random_percent);
	}
    }

    for (i = 0; i < m; i++) {
	phgBalanceGrid(g, 1.0, 0, NULL, 0.);
	phgCheckConformity(g);
#if USE_MPI
	MPI_Bcast(&i, sizeof(i), MPI_BYTE, 0, g->comm);
#endif
	if (refine_type == 0)
	    phgRefineAllElements(g, 1);
	else
	    phgRefineRandomElements(g, random_percent);
	phgCheckConformity(g);
    }

    /*phgDumpGrid(g);*/
    PrintTime(1);

    /* print some statistics on the initial mesh */
    phgPrintf("********** Final mesh ************\n");
    phgPrintf("Number of elements: %d\n", g->nleaf_global);
    phgPrintf("Number of vertices: %d\n", g->nvert_global);
    phgPrintf("Number of edges:    %d\n", g->nedge_global);
    phgPrintf("Number of faces:    %d\n", g->nface_global);
    phgPrintf("**********************************\n");

    /* v holds the generation of the element */
    v = phgDofNew(g, DOF_P0, 1, "generation", DofNoAction);
    ForAllElements(g, e)
	*DofElementData(v, e->index) = (FLOAT)(e->generation);

    if (output_file != NULL) {
	switch (output_type) {
	    case 0:
		break;
	    case 1:
		phgPrintf("Write VTK file \"%s\".\n",
				phgExportVTK(g, output_file, u, v, NULL));
		break;
	    case 2:
		phgPrintf("Write OpenDX file \"%s\".\n", 
				phgExportDX(g, output_file, u, v, NULL));
		break;
	}
    }

    PrintTime(1);

    phgDofFree(&u);
    phgDofFree(&v);
    phgFreeGrid(&g);

    PrintTime(1);

    phgFinalize();
    
    Return 0;
}
