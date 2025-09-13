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

/* $Id: hp_test.c,v 1.66 2020/12/06 09:03:08 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static void
func(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    int order = DOF_DEFAULT->order;
    FLOAT px = Pow(x, order), py = Pow(y, order), pz = Pow(z, order);

    if (DOF_DEFAULT == DOF_HC0) {
	/* (1,2,3) x (x,y,z) + (4,5,6) */
	*(values++) = 2. * z - 3. * y + 4.;
	*(values++) = 3. * x - 1. * z + 5.;
	*(values) = 1. * y - 2. * x + 6.;
	return;
    }

    *(values++) = px + py + pz;
    *(values++) = px + py - pz;
    *(values++) = px - py - pz;
}

int
main(int argc, char *argv[])
{
    ELEMENT *e;
    GRID *g;
    DOF *u, *v;
    HP_TYPE *hp;
    INT pre_refines = 0, post_refines = 0;
    const char *fn = "cube.dat";
    INT i, n;

    phgOptionsPreset("-dof_type HB2");
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-post_refines", "Post-refines", &post_refines);
    phgInit(&argc, &argv);

    if (DOF_DEFAULT->hp_info == NULL)
	phgError(1, "\"%s\" is not hierarchical, abort.\n", DOF_DEFAULT->name);

    g = phgNewGrid(-1);
    if (argc > 1)
	fn = argv[1];
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgRefineAllElements(g, pre_refines);
    v = phgDofNew(g, DOF_DEFAULT, 3 / DOF_DEFAULT->dim, "v", DofInterpolation);
    phgDofSetDataByFunction(v, func);
    hp = phgHPNew(g, v->type->hp_info);
    ForAllElements(g, e)
	e->hp_order = DofElementOrder(v, e->index) +
				GlobalElement(g, e->index) % 5;
    phgHPSetup(hp, FALSE);
#if 0
    u = phgDofCopy(v, NULL, NULL, "u");
    phgHPAttachDof(hp, u, TRUE);
#else
    u = phgHPDofNew(g, hp, DofDim(v) / hp->info->dim, "u", DofInterpolation);
    phgDofSetDataByFunction(u, func);
#endif
    phgHPFree(&hp);

#if 0
    /* This will produce an uniform order h-p DOF */
    ForAllElements(g, e)
	e->hp_order = DofElementOrder(v, e->index) + 0;
    phgHPSetup(u->hp, TRUE);
#endif

#if 0
if (phgRank == 0) phgDumpGrid(g);
for (i = 0; i < g->nedge; i++) if (g->types_edge[i] != UNREFERENCED) phgPrintf("edge %d, order %d\n", GlobalEdge(g, i), u->hp->edge_order[i]);
for (i = 0; i < g->nface; i++) if (g->types_face[i] != UNREFERENCED) phgPrintf("face %d, order %d\n", GlobalFace(g, i), u->hp->face_order[i]);
for (i = 0; i < g->nelem; i++) if (g->types_elem[i] != UNREFERENCED) phgPrintf("face %d, order %d\n", GlobalElement(g, i), u->hp->elem_order[i]);
#endif

    //phgDofSetDataByFunction(u, func);
    phgBalanceGrid(g, 1.0, 0, NULL, 0.);

    phgPrintf("Initial mesh: %d elems, %d procs, %d non hp DOF, %d hp DOF\n",
		g->nelem_global, g->nprocs, DofDataCountGlobal(v),
		DofDataCountGlobal(u));

    for (i = 0; i < post_refines; i++) {
	phgRefineSurface(g);
	//phgRefineAllElements(g, 1);
	phgBalanceGrid(g, 1.0, 0, NULL, 0.);
    }

    phgPrintf("Finest mesh:  %d elems, %d procs, %d non hp DOF, %d hp DOF\n",
		g->nelem_global, g->nprocs, DofDataCountGlobal(v),
		DofDataCountGlobal(u));

    /* test phgExportVTK */
    phgExportVTK(g, "hp_test.vtk", u, v, NULL);

    /* test coarsening */
    n = -1;
    while (TRUE) {
	ForAllElements(g, e)
	    e->mark = -3;
	phgCoarsenMarkedElements(g);
	if (g->nelem_global == n)
	    break;
	n = g->nelem_global;
    }

    phgPrintf("Coarsening:   %d elems, %d procs, %d non hp DOF, %d hp DOF\n",
		g->nelem_global, g->nprocs, DofDataCountGlobal(v),
		DofDataCountGlobal(u));

#if 0
    ForAllElements(g, e)
	e->hp_order = DOF_DEFAULT->order + 3;
    phgHPSetup(u, TRUE);
    //phgDofSetDataByFunction(u, func);
#endif

#if 0
    /* Print out orders for checking.
     * Note: the outputs below should be identical for nprocs > 1, but
     * 	     different for nprocs == 1 because elements are reordered
     * 	     according to the Hamilton or Hilbert path for nprocs > 1 */
    ForAllElements(g, e) {
	int i;
	for (i = 0; i < NEdge; i++) {
	    BTYPE b = g->types_edge[e->edges[i]];
	    phgInfo(-1, "Element %d: order %d, edge %d: order %d %c\n",
			GlobalElement(g, e->index),
			u->hp->elem_order[e->index],
			GlobalEdge(g, e->edges[i]),
			u->hp->edge_order[e->edges[i]],
			(b & OWNER) ? '+' : '-');
	}
	for (i = 0; i < NFace; i++) {
	    BTYPE b = g->types_face[e->faces[i]];
	    phgInfo(-1, "Element %d: order %d, face %d: order %d %c\n",
			GlobalElement(g, e->index),
			u->hp->elem_order[e->index],
			GlobalFace(g, e->faces[i]),
			u->hp->face_order[e->faces[i]],
			(b & OWNER) ? '+' : '-');
	}
    }
#endif
    phgPrintf("Vertex  DOF: %d\n", DofVertexDataCountGlobal(u));
    phgPrintf("Edge    DOF: %d\n", DofEdgeDataCountGlobal(u));
    phgPrintf("Face    DOF: %d\n", DofFaceDataCountGlobal(u));
    phgPrintf("Element DOF: %d\n", DofElementDataCountGlobal(u));
    phgPrintf("Total   DOF: %d\n", DofDataCountGlobal(u));

    phgPrintf("Test evaluation: ");
    ForAllElements(g, e) {
	static FLOAT lambda[] = {0.25, 0.25, 0.25, 0.25};
	int dim = DofDim(v);
	FLOAT x, y, z, value0[dim * Dim], value1[dim * Dim], error, norm;
	phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
	phgDofEval(v, e, lambda, value0);
	phgDofEval(u, e, lambda, value1);
	norm = error = 0.0;
	for (i = 0; i < dim; i++) {
	    norm += value0[i] * value0[i];
	    error += Pow(value0[i] - value1[i], 2);
	}
	norm = Sqrt(norm);
	error = Sqrt(error) / (norm < 1e-8 ? 1.0 : norm);
	if (error > 1e-6) {
	    printf("*** proc %d, ERROR: wrong function values (error = %lg).\n",
			g->rank, (double)error);
	    break;
	}
	phgDofEvalGradient(v, e, lambda, NULL, value0);
	phgDofEvalGradient(u, e, lambda, NULL, value1);
	norm = error = 0.0;
	for (i = 0; i < dim * Dim; i++) {
	    norm += value0[i] * value0[i];
	    error += Pow(value0[i] - value1[i], 2);
	}
	norm = Sqrt(norm);
	error = Sqrt(error) / (norm < 1e-8 ? 1.0 : norm);
	if (error > 1e-6) {
	    printf("*** proc %d, ERROR: wrong gradient values (error = %lg).\n",
			g->rank, (double)error);
	    break;
	}
    }
    phgPrintf("done.\n");

    /* test gradient/divergence/curl operation */
    if (TRUE) {
	DOF *du = NULL, *dv = NULL;

	phgPrintf("Test gradient:   ");
	phgDofGradient(u, &du, NULL, NULL);
	phgDofGradient(v, &dv, NULL, NULL);
	phgPrintf("|u|=%14.7le, |v|=%14.7le, ",
			(double)phgDofNormL2(du), (double)phgDofNormL2(dv));
	phgDofAXPY(-1.0, dv, &du);
	phgPrintf("|v-u|=%14.7le\n", (double)phgDofNormL2(du));

	phgPrintf("Test divergence: ");
	phgDofDivergence(u, &du, NULL, NULL);
	phgDofDivergence(v, &dv, NULL, NULL);
	phgPrintf("|u|=%14.7le, |v|=%14.7le, ",
			(double)phgDofNormL2(du), (double)phgDofNormL2(dv));
	phgDofAXPY(-1.0, dv, &du);
	phgPrintf("|v-u|=%14.7le\n", (double)phgDofNormL2(du));

	phgPrintf("Test curl:       ");
	phgDofCurl(u, &du, NULL, NULL);
	phgDofCurl(v, &dv, NULL, NULL);
	phgPrintf("|u|=%14.7le, |v|=%14.7le, ",
			(double)phgDofNormL2(du), (double)phgDofNormL2(dv));
	phgDofAXPY(-1.0, dv, &du);
	phgPrintf("|v-u|=%14.7le\n", (double)phgDofNormL2(du));

	phgDofFree(&du);
	phgDofFree(&dv);
    }

    phgPrintf("Test AXPBY:      ");
    phgPrintf("|u|=%14.7le, |v|=%14.7le, ",
			(double)phgDofNormL2(u), (double)phgDofNormL2(v));
    phgDofAXPY(-1.0, v, &u);
    phgPrintf("|v-u|=%14.7le\n", (double)phgDofNormL2(u));

    phgDofFree(&v);
    phgDofFree(&u);
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
