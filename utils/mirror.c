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

/* $Id: mirror.c,v 1.19 2021/02/25 01:47:19 zlb Exp $ 
 *
 * This utility program creates mirrored mesh. The resulting mesh satisfies
 * periodic conformity constraints and can be used with periodic boundary
 * conditions. */

#include "phg.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

static INT depth = 0;
static FLOAT eps0 = 1e-6;
static char *vtk_file = NULL;
static char *dx_file = NULL;

static const char *symmetry_names[] = {"min", "max", NULL};
static int symmetry = 0;

static const char *format_names[] = {"medit", "albert", NULL};
static int format = 0;

/* The name of the file containing the direction vectors */
static char *dir_file = NULL;

/* Transposed matrix transforming the orthobrick into the parallelepiped,
 * whose rows are composed of the three direction vectors */
static FLOAT (*itran)[Dim] = NULL;

/* Transposed matrix transforming the parallelepiped into the orthobrick,
 * which is the inverse of 'tran' */
static FLOAT (*tran)[Dim] = NULL;

#undef min
#undef max
static FLOAT min = 0., max = 0.;

static INT *
vertices_on_symmetry_plane(GRID *g, int dir)
/* 'dir' should be 0 (x), 1 (y), or 2 (z).
 *
 * Returns a dynamically allocated INT array of size g->nvert, in which the
 * entries are -1 if corresponding vertices are on the symmetry plane,
 * and 0 otherwise.
 *
 * This function also sets up 'min', 'max' */
{
    FLOAT x, x0, eps;
    INT i, *vmap;

    if (g->nvert == 0)
	return NULL;

    if (tran == NULL) {
	x0 = g->bbox[symmetry][dir];
	eps  = (g->bbox[1][dir] - g->bbox[0][dir]) * eps0;
    }
    else {
	min = max = tran[0][dir] * g->verts[0][0] +
		    tran[1][dir] * g->verts[0][1] +
		    tran[2][dir] * g->verts[0][2];
	for (i = 1; i < g->nvert; i++) {
	    x = tran[0][dir] * g->verts[i][0] +
		tran[1][dir] * g->verts[i][1] +
		tran[2][dir] * g->verts[i][2];
	    if (min > x)
		min = x;
	    else if (max < x)
		max = x;
	}
	x0 = (symmetry == 0 ? min : max);
	eps  = (max - min) * eps0;
    }

    vmap = phgCalloc(g->nvert, sizeof(*vmap));

    for (i = 0; i < g->nvert; i++) {
	if (tran == NULL)
	    x = g->verts[i][dir];
	else
	    x = tran[0][dir] * g->verts[i][0] +
		tran[1][dir] * g->verts[i][1] +
		tran[2][dir] * g->verts[i][2];
	if (Fabs(x - x0) < eps)
	    vmap[i] = -1;
    }

    return vmap;
}

int
main(int argc, char *argv[])
{
    GRID *g;
    int ii, c, dir;
    const char *in, *out, *dirs = NULL;
    FILE *f;
    INT i, n, *vmap;
    ELEMENT *roots, *e, *e1;

    phgOptionsRegisterInt("-refine_depth", "Pre-refinement depth", &depth);    
    phgOptionsRegisterFloat("-tolerance", "Tolerance in comparing coordinates",
				&eps0);    
    phgOptionsRegisterKeyword("-symmetry_plane", "Symmetry plane",
				symmetry_names, &symmetry);    
    phgOptionsRegisterKeyword("-output_format", "Output format",
				format_names, &format);    
    phgOptionsRegisterFilename("-dir_file", "Filename for the direction "
			       "vectors", &dir_file);    
    phgOptionsRegisterFilename("-vtk_file", "VTK filename", &vtk_file);    
    phgOptionsRegisterFilename("-opendx_file", "OpenDX filename", &dx_file);    

    phgOptionsPreset("+warn_nonopt");
    phgInit(&argc, &argv);

    if (phgNProcs > 1)
	phgError(1, "parallel execution of this program is not allowed.\n");

    if (dir_file != NULL) {
	double d;
	FLOAT (*A)[Dim] = phgAlloc(Dim * sizeof(*A));

	/* read the direction vectors */
	phgPrintf("Loading directions from \"%s\".\n", dir_file);
	itran = phgAlloc(Dim * sizeof(*itran));
	FreeAtExit(itran);
	if ((f = fopen(dir_file, "r")) == NULL)
	    phgError(1, "cannot read file \"%s\".\n", dir_file);
	for (i = 0; i < Dim; i++) {
	    for (ii = 0; ii < Dim; ii++) {
		if (fscanf(f, "%le", &d) != 1)
		    phgError(1, "failure reading file \"%s\".\n", dir_file);
		A[i][ii] = itran[i][ii] = d;
	    }
	}
	fclose(f);

	/* Compute the matrix 'tran' (inverse of itran). */
	tran = phgCalloc(Dim, sizeof(*tran));
	FreeAtExit(tran);
	tran[0][0] = tran[1][1] = tran[2][2] = 1.;
	if (!phgSolverDenseSolver(Dim, Dim, (FLOAT *)A, (FLOAT *)tran))
	    phgError(1, "periodic directions are dependent!\n");
	phgFree(A);
    }

    if (argc != 3 && argc != 4) {
	phgPrintf("This program creates a mirror mesh with respect to given "
		  "space directions.\n");
	phgPrintf("Type \"%s -help user\" for cmdline options.\n", argv[0]);
	phgPrintf("Usage:\n");
	phgPrintf("\t%s input_mesh output_mesh [dirs] [options]\n",
		  argv[0]);
        phgPrintf("\t(dirs is a string composed of 'x', 'y' and 'z')\n");
	phgPrintf("Examples:\n");
	phgPrintf("\t%s cube.dat out.mesh xyz -refine_depth 3\n", argv[0]);
	phgPrintf("\t%s cube.dat out.mesh xxx\n", argv[0]);
	phgFinalize();
	exit(1);
    }

    in = argv[1];
    out = argv[2];
    if (argc == 4)
	dirs = argv[3];

    /* confirm overwriting existing output file */
    if ((f = fopen(out, "r")) != NULL && phgRank == 0) {
	char answer[2];
	fclose(f);
	fprintf(stderr, "Overwrite output file \"%s\" (y/n) ? ", out);
	fgets(answer, sizeof(answer), stdin);
	if (toupper(answer[0]) != 'Y')
	    phgAbort(1);
    }

    /* Note: FACE_FLAG is needed by phgDumpMedit */
    g = phgNewGrid(VERT_FLAG | FACE_FLAG | ELEM_FLAG);
    if (!phgImport(g, in, FALSE))
	phgError(1, "can't read file \"%s\".\n", in);

    phgPrintf("Input mesh: %d elements, %d vertices.\n", g->nelem, g->nvert);

    if (depth > 0) {
	GRID *g1;
	phgRefineAllElements(g, depth);
	/* flatten the mesh */
	g1 = phgDupGrid(g, TRUE);
	phgFreeGrid(&g);
	g = g1;
	phgPrintf("After %d refinements: %d elements, %d vertices.\n",
			depth, g->nelem, g->nvert);
    }

    /* reset neighbours and bound_type (since g->roots will be reallocated),
     * and reset faces[] and edges[] in elements. */
    for (i = 0; i < g->nroot; i++) {
	e = g->roots + i;
	for (ii = 0; ii < NFace; ii++) {
	    e->neighbours[ii] = NULL;
	    e->bound_type[ii] &= BDRY_MASK;
	    e->faces[ii] = -1;
	}
	for (ii = 0; ii < NEdge; ii++) {
	    e->edges[ii] = -1;
	}
    }
    g->nedge = g->nedge_global = g->nface = g->nface_global = 0;

    while ( dirs != NULL && (c = tolower(*(dirs++))) != '\0') {
	if (c != 'x' && c != 'y' && c != 'z')
	    phgError(1, "invalid argument \"%s\" (should only contain 'x', "
			"'y', and 'z').\n", argv[3]);
	dir = c - 'x';
	vmap = vertices_on_symmetry_plane(g, dir);

	/* count # of vertices not on the symmetry plane, and set up vmap[] */
	n = g->nvert;
	for (i = 0; i < g->nvert; i++) {
	    if (vmap[i] == -1)
		continue;
	    vmap[i] = n++;
	}
	/* prepare g->verts */
	g->verts = phgRealloc_(g->verts, n * sizeof(*g->verts),
					 g->nvert * sizeof(*g->verts));
	for (i = 0; i < g->nvert; i++) {
	    COORD *c0, *c1;
	    if (vmap[i] == -1)
		continue;
	    c0 = g->verts + i;
	    c1 = g->verts + vmap[i];
	    if (tran == NULL) {
		memcpy(c1, c0, sizeof(*c1));
		(*c1)[dir] = 2. * g->bbox[symmetry][dir] - (*c1)[dir];
	    }
	    else {
		COORD c2;
		for (ii = 0; ii < Dim; ii++) {
		    c2[ii] = tran[0][ii] * (*c0)[0] +
			     tran[1][ii] * (*c0)[1] +
			     tran[2][ii] * (*c0)[2];
		}
		c2[dir] = 2. * (symmetry == 0 ? min : max) - c2[dir];
		for (ii = 0; ii < Dim; ii++) {
		    (*c1)[ii] = itran[0][ii] * c2[0] +
				itran[1][ii] * c2[1] +
				itran[2][ii] * c2[2];
		}
	    }
	}
	g->nvert_global = g->nvert = n;

	roots = phgNewElements(2 * g->nroot);
	memcpy(roots, g->roots, g->nroot * sizeof (*roots));
	phgFree(g->roots);
	g->roots = roots;
	for (i = 0; i < g->nroot; i++) {
	    e = g->roots + i;
	    e1 = g->roots + g->nroot + i;
	    e1->index = g->nroot + i;
	    e1->region_mark = e->region_mark;
	    memcpy(e1->bound_type, e->bound_type, sizeof(e->bound_type));
	    for (ii = 0; ii < NVert; ii++) {
		if (vmap[e->verts[ii]] >= 0)
		    e1->verts[ii] = vmap[e->verts[ii]];
		else
		    e1->verts[ii] = e->verts[ii];
	    }
	}
	g->nelem_global = g->nleaf = g->nelem = (g->nroot *= 2);

	/* update g->bbox for next mirroring operation */
	if (*dirs != '\0' && tran == NULL)
	    g->bbox[symmetry][dir] = 2. * g->bbox[symmetry][dir] -
					  g->bbox[1 - symmetry][dir];
	else
	    ;	/* do nothing */
	g->volume *= 2.;
	phgPrintf("Mirroring in %c direction: %d elements, %d vertices.\n",
			c, g->nelem, g->nvert);

	phgFree(vmap);
    }

    phgUpdateNeighbours(g);
    phgUpdateEdges(g);
    phgUpdateFaces(g);
    phgUpdateBoundaryTypes(g);
    phgRefineInit(g, TRUE);

    (format == 0 ? phgExportMedit : phgExportALBERT)(g, out);
    phgPrintf("Mesh file \"%s\" created.\n", out);

    if (vtk_file != NULL)
	phgExportVTK(g, vtk_file, NULL);
    if (dx_file != NULL)
	phgExportDX(g, dx_file, NULL);

    phgFreeGrid(&g);
    phgFinalize();
    
    return 0;
}
