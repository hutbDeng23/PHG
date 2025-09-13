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

/*
 * This program tests the phgMap*NeighbourMap functions.
 *
 * $Id: map_test1.c,v 1.24 2022/04/02 06:20:33 zlb Exp $ */

#include "phg.h"
#include <math.h>	/* pow() */
#include <string.h>	/* memcpy() */

static FLOAT
compute_jump(NEIGHBOUR_MAP *nm, int dof_id, ELEMENT *e, int face, DOF_PROJ proj)
{
    GRID *g;
    ELEMENT *e1;
    DOF *dof, *dof1;
    MAP *map, *map1;
    NEIGHBOUR_INFO *ni;
    QUAD *q;
    int i, j, k, nbas, face1, u0, u1, u2, v0, v1, v2;
    SHORT *bi, *bi1;
    FLOAT jump, d, lambda[Dim + 1], lambda1[Dim + 1];
    const FLOAT *n = NULL;

    map = nm->map_local;
    dof = map->dofs[dof_id];
    g = nm->g_local;

    if (proj != DOF_PROJ_NONE)
	n = phgGeomGetFaceNormal(g, e, face);

    nbas = phgDofGetBasesOnFace(dof, e, face, NULL);
    bi = phgAlloc(2 * nbas * sizeof(*bi));
    bi1 = bi + nbas;
    phgDofGetBasesOnFace(dof, e, face, bi);

    /* get the neighour */
    ni = phgMapNeighbourMap(nm, e, face);
    assert(ni[0].e != NULL && ni[1].e == NULL); /* only one neighbour */
    e1 = ni[0].e;
    face1 = ni[0].index;
    map1 = ni[0].is_remote ? nm->map_remote : nm->map_local;
    dof1 = map1->dofs[dof_id];
    assert(nbas == phgDofGetBasesOnFace(dof1, e1, face1, NULL));
    phgDofGetBasesOnFace(dof1, e1, face1, bi1);

    GetFaceVertices(e, face, u0, u1, u2, i);
    lambda[face] = 0.0;

    GetFaceVertices(e1, face1, v0, v1, v2, i);
    lambda1[face1] = 0.0;

    jump = 0.0;
    q = phgQuadGetQuad2D(2 * dof->type->order);
    for (k = 0; k < q->npoints; k++) {
	lambda[u0] = lambda1[v0] = q->points[k * 3 + 0];
	lambda[u1] = lambda1[v1] = q->points[k * 3 + 1];
	lambda[u2] = lambda1[v2] = q->points[k * 3 + 2];
#if DEBUG
    {
	FLOAT x, y, z, x1, y1, z1;
	phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
	phgGeomLambda2XYZ(dof1->g, e1, lambda1, &x1, &y1, &z1);
	assert(Fabs(x-x1) + Fabs(y-y1) + Fabs(z-z1) <= 3. * FLOAT_EPSILON);
    }
#endif	/* DEBUG */
	d = 0.0;
	for (j = 0; j < nbas; j++) {
	    int dim = dof->type->dim;
	    FLOAT bv[dim];	/* values of the basis in e */
	    const FLOAT *p;	/* values of the basis in e1 */
#if DEBUG
	    assert(dof->type->fe_space != FE_L2 /* DG */ ||
			phgMapE2G(map, dof_id, e, bi[j]) ==
				phgMapE2G(map1, dof_id, e1, bi1[j]));
#endif	/* DEBUG */
	   /* WARNING: must save the returned values of dof->type->BasFuncs(),
	    * otherwise they would be overwritten by the following call to
	    * dof1->type->BasFuncs(), since dof->type == dof1->type and
	    * BasFuncs() returns the pointer to an internal static buffer. */
	    p = dof->type->BasFuncs(dof, e, bi[j], bi[j] + 1, lambda);
	    memcpy(bv, p, sizeof(bv));
	    p = dof1->type->BasFuncs(dof1, e1, bi1[j], bi1[j] + 1, lambda1);
	    /* compute and store the jump in bv[] */
	    for (i = 0; i < dim; i++)
		bv[i] -= p[i];
	    switch (proj) {
		case DOF_PROJ_NONE:	/* no projection */
		    for (i = 0; i < dim; i++)
			d += bv[i] * bv[i];
		    break;
		case DOF_PROJ_DOT:	/* normal projection */
		    assert(dim == 3);
		    d += Pow(bv[0] * n[0] + bv[1] * n[1] + bv[2] * n[2], 2);
		    break;
		case DOF_PROJ_CROSS:	/* tangential projection */
		    assert(dim == 3);
		    d += Pow(bv[1] * n[2] - bv[2] * n[1], 2) +
			 Pow(bv[2] * n[0] - bv[0] * n[2], 2) +
			 Pow(bv[0] * n[1] - bv[1] * n[0], 2);
		    break;
		default:
		    phgError(1, "Invalid projection type (%d)\n", proj);
	    }
	}
	jump += d * q->weights[k];
    }
    phgFree(bi);

    return jump * phgGeomGetFaceArea(g, e, face);
}

int
main(int argc, char *argv[])
{
    char *fn = "tetra.dat";
    INT pre_refines = 1;
    char *dof_u = "HC2", *dof_v = "HD2";
    int proj = 0;
    const char *proj_names[] = {"none", "dot", "cross", NULL};
    const DOF_PROJ proj_types[] = {DOF_PROJ_NONE, DOF_PROJ_DOT, DOF_PROJ_CROSS};

    GRID *g;
    DOF *u, *v;
    MAP *map;
    NEIGHBOUR_MAP *nm;		/* struct containing remote neighbours */
    VEF_MAP *vef_map;		/* face -> element map */
    FLOAT jumpu, jumpv;
    INT i;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", &fn);
    phgOptionsRegisterInt("-pre_refines", "Pre-refinements", &pre_refines);
    phgOptionsRegisterString("-dof_u", "DOF type for u", &dof_u);
    phgOptionsRegisterString("-dof_v", "DOF type for v", &dof_v);
    phgOptionsRegisterKeyword("-projection", "Projection (w.r.t. face normal)",
				proj_names, &proj);

    phgInit(&argc, &argv);
    g = phgNewGrid(-1);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgRefineAllElements(g, pre_refines);

    phgBalanceGrid(g, 1.1, 0, NULL, 0.);
    phgPrintf("Partition mesh, %d submeshes, load imbalance: %lg\n",
	      g->nprocs, (double)g->lif);

    phgOptionsSetHandler("-dof_type", dof_u);
    u = phgDofNew(g, DOF_DEFAULT, 1, "u", DofInterpolation);
    phgOptionsSetHandler("-dof_type", dof_v);
    v = phgDofNew(g, DOF_DEFAULT, 1, "v", DofInterpolation);
    phgPrintf("u->type=%s, v->type=%s, projection=%s\n",
			u->type->name, v->type->name, proj_names[proj]);

    map = phgMapCreate(u, v, NULL);
    nm = phgMapInitNeighbourMap(map, FACE, /*dof_data=*/FALSE);

    /* compute:
     *    jumpu = \sum_{f\in faces} \int_f \sum_{b\in bases of u} |jump(b)|^2
     *    jumpv = \sum_{f\in faces} \int_f \sum_{b\in bases of v} |jump(b)|^2 */
    jumpu = jumpv = 0.0;
    vef_map = phgDofSetupVEFMap(g, NULL, FACE_FLAG);
    for (i = 0; i < g->nface; i++) {
	if ((g->types_face[i] & OWNER) == FALSE)
	    continue;	/* face owned by another process */
	if ((g->types_face[i] & INTERIOR) == FALSE)
	    continue;	/* boundary face */
	jumpu += compute_jump(nm, 0, vef_map->Fmap[i], vef_map->Find[i],
				proj_types[proj]);
	jumpv += compute_jump(nm, 1, vef_map->Fmap[i], vef_map->Find[i],
				proj_types[proj]);
    }
    phgDofFreeVEFMap(&vef_map);
#if USE_MPI
    if (g->nprocs > 1) {
	FLOAT a[2] = {jumpu, jumpv}, b[2];
	MPI_Reduce(a, b, 2, PHG_MPI_FLOAT, PHG_SUM, 0, g->comm);
	jumpu = b[0];
	jumpv = b[1];
    }
#endif	/* USE_MPI */
    phgPrintf("jumpu = %lg, jumpv = %lg\n", (double)jumpu, (double)jumpv);

    phgMapReleaseNeighbourMap(&nm);
    phgMapDestroy(&map);

    phgDofFree(&u);
    phgDofFree(&v);

    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
