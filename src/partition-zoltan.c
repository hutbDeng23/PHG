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

/* $Id: partition-zoltan.c,v 1.29 2014/03/31 12:10:42 wleng Exp $
 *
 * The Zoltan interface
 */

#include "phg.h"

#if USE_MPI	/* whole file */
#if USE_ZOLTAN	/* whole file */

/* Note: the REFTREE partitioner is not truely scalapable since it requires
 *	 all children of root elements, thus it's not implemented yet. */
#define ENABLE_REFTREE	0 /* remove or change to 1 after REFTREE is done */

#include <zoltan.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

static GRID *g_;
static FLOAT pwr = 0.0;

typedef enum {
    LB_RCB, LB_RIB, LB_HSFC,	/* Geometric */
#if 0 * ENABLE_REFTREE
    LB_REFTREE,			/* Refinement tree */
#endif	/* ENABLE_REFTREE */
    /*LB_GRAPH,*/ LB_HYPERGRAPH,/* Graph */
    /*LB_HIER,*/
    LB_BLOCK, LB_RANDOM		/* Test */
} LB_METHOD;
static const char *partitioners[] = {
    "rcb", "rib", "hsfc",
#if 0 * ENABLE_REFTREE
    "reftree",
#endif	/* ENABLE_REFTREE */
    /*"graph",*/ "hypergraph",
    /*"hier",*/
    "block", "random",
    NULL
};
static int partitioner = 0;

/* Common query functions */

static int
get_nlocal(void *userdata, int *err)
{
    *err = 0;
    return g_->nleaf;
}

static void
get_local_list(void *userdata, int gdim, int ldim, ZOLTAN_ID_PTR gids,
		ZOLTAN_ID_PTR lids, int wgt_dim, float *wgts, int *err)
{
    ELEMENT *e;
    DOF *w = userdata;

    assert(gdim == 1);
    assert((w == NULL) == (wgts == NULL));
    assert(w == NULL || w->count_elem == 1);

    *err = 0;
    ForAllElements(g_, e) {
	*(gids++) = GlobalElement(g_, e->index);
	*(lids++) = e->index;
	if (wgts != NULL)
	    *(wgts++) = (float)(Fabs(Pow(*DofElementData(w, e->index), pwr)));
    }
    
    return;
}

/* Query functions for RCB/RIB/HSFC/REFTREE */

static int
geom_num(void *data, int *err)
{
    *err = 0;
    return 3;
}

static void
geom_multi(void *data, int num_gid_entries, int num_lid_entries,
	   int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
	   int num_dim, double *geom_vec, int *err)
{
    int i;
    FLOAT x, y, z;
    static const FLOAT lam[] = {0.25, 0.25, 0.25, 0.25};

    *err = 0;
    for (i = 0; i < num_obj; i++) {
	/* return the bary-center of elements */
	phgGeomLambda2XYZ(g_, g_->elems[local_ids[i]], lam, &x, &y, &z);
	*(geom_vec++) = x;
	*(geom_vec++) = y;
	*(geom_vec++) = z;
    }

    return;
}

#if 0 * ENABLE_REFTREE
/* Query functions for REFTREE */

static int
num_coarse_obj(void *data, int *err)
{
    *err = 0;
    return g->nroot_global;
}

static void
coarse_obj_list(void *data, int num_gid_entries, int num_lid_entries,
		ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids,
		int *assigned, int *num_vert, ZOLTAN_ID_PTR vertices,
		int *in_order, ZOLTAN_ID_PTR in_vertex,
		ZOLTAN_ID_PTR out_vertex, int *err)
{
    int i;
    ELEMENT *e;

    *err = 0;
    for (i = 0, e = g_->roots; i < g_->nroot; i++, e++) {
	*(global_ids++) = GlobalElement(g_, e->index);	/* index */
	*(global_ids++) = 0;				/* generation */
	*(global_ids++) = GlobalElement(g_, e->index);	/* index */
	*(global_ids++) = 0;				/* generation */
    }
}
#endif	/* ENABLE_REFTREE */

/* Graph/Hypergraph query functions */

static void
HG_size_CS(void *data, int *num_lists, int *num_pins, int *format, int *err)
{
    ELEMENT *e;
    int i, pins;

    *err = 0;
    *format = ZOLTAN_COMPRESSED_VERTEX;	/* use vertex-edges list format */
    *num_lists = g_->nleaf;
    /* count number of pins (neighbours) */
    pins = 0;
    ForAllElements(g_, e) {
	for (i = 0; i < NFace; i++)
	    if (!IS_BDRY(e->bound_type[i]))
		pins++;
    }
    *num_pins = pins;
}

static void
HG_CS(void *data, int num_gid_entries, int nvtxedge, int npins, int format,
      ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr, ZOLTAN_ID_PTR pin_GID,
      int *err)
{
    ELEMENT *e;
    int i, pins;

    *err = 0;
    assert((g_->flags & FACE_FLAG));
    assert(format == ZOLTAN_COMPRESSED_VERTEX);
    assert(num_gid_entries == 1);
    assert(nvtxedge == g_->nleaf);

    pins = 0;
    ForAllElements(g_, e) {
	*(vtxedge_GID++) = GlobalElement(g_, e->index);
	*(vtxedge_ptr++) = pins;
	for (i = 0; i < NFace; i++) {
	    if (IS_BDRY(e->bound_type[i]))
		continue;
	    *(pin_GID++) = GlobalFace(g_, e->faces[i]);
	    pins++;
	}
    }
    assert(pins == npins);
}

/* Parameters for HyperGraph method */
static const char *phg_coarsening_methods[] = {
    "agg",		/* good quality, slow */
    "ipm",		/* good quality, slow, default */
    "l-ipm",		/* poor quality, faster */
    "a-ipm",		/* compromise of IPM and L-IPM */
    NULL
};
static int phg_coarsening_method = 1;

static const char *phg_refinement_methods[] = {
    "fm",	/* approximate Fiduccia-Mattheyses, default */
    "no",	/* no refinement */
    NULL
};
static int phg_refinement_method = 0;

static FLOAT phg_refinement_quality = 1.0;
static BOOLEAN phg_check_hypergraph = FALSE;

/*--------------------------------------------------------------------------*/

void
phgPartitionZoltanInitialize(int argc, char **argv)
{
    float ver;

    if (Zoltan_Initialize(argc, argv, &ver) != ZOLTAN_OK)
	phgError(1, "error in Zoltan_Initialize.\n");
}

BOOLEAN
phgPartitionZoltan(GRID *g, MPI_Comm newcomm, DOF *weights, FLOAT power)
{
    int i, rc, nprocs;
    char str[1024];
    int changed, ngids, nlids, nimport, nexport;
    struct Zoltan_Struct *zs;
    ZOLTAN_ID_PTR import_gids, import_lids, export_gids, export_lids;
    int *import_procs, *import_parts, *export_procs, *export_parts;

    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	/* register Zoltan options */
	initialized = TRUE;
	phgOptionsRegisterTitle("\n[Zoltan options]", "\n", "partition");
	phgOptionsRegisterKeyword("-zoltan_method", "Zoltan method",
		partitioners, &partitioner);
	phgOptionsRegisterTitle("\n[Zoltan HyperGraph options]",
				"\n", "partition");
	phgOptionsRegisterKeyword("-zoltan_coarsening_method",
		"Coarsening method (l-ipm is the fastest)",
		phg_coarsening_methods, &phg_coarsening_method);
	phgOptionsRegisterKeyword("-zoltan_refinement_method",
		"Refinement method",
		phg_refinement_methods, &phg_refinement_method);
	phgOptionsRegisterFloat("-zoltan_refinement_quality",
		"Refinement quality (>0)",
		&phg_refinement_quality);
	phgOptionsRegisterNoArg("-zoltan_check_hypergraph",
		"Check HyperGraph", &phg_check_hypergraph);
	return FALSE;
    }

    MPI_Comm_size(newcomm, &nprocs);

    if (g->comm == MPI_COMM_NULL)
	return FALSE;

    if (g->nprocs == 1 && nprocs == 1)
	return FALSE;

    if (power == 0.0)
	weights = NULL;
    g_ = g;
    pwr = power;

    zs = Zoltan_Create(g->comm);

    Zoltan_Set_Param(zs, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zs, "LB_METHOD", partitioners[partitioner]);
    Zoltan_Set_Param(zs, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param(zs, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zs, "OBJ_WEIGHT_DIM", weights == NULL ? "0" : "1");
    Zoltan_Set_Param(zs, "RETURN_LISTS", "PARTITION ASSIGNMENTS");
    sprintf(str, "%d", nprocs);
    Zoltan_Set_Param(zs, "NUM_GLOBAL_PARTITIONS", str);
    Zoltan_Set_Param(zs, "LB_APPROACH",
			g->nprocs == 1 ? "PARTITION" : "REPARTITION");

    /* Common query functions required by all (except REFTREE) partitioners */
    Zoltan_Set_Num_Obj_Fn(zs, get_nlocal, NULL);
    Zoltan_Set_Obj_List_Fn(zs, get_local_list, weights);

    switch ((LB_METHOD)partitioner) {
	case LB_RCB:
	case LB_RIB:
	case LB_HSFC:
	    /* Mandatory query functions */
	    Zoltan_Set_Num_Geom_Fn(zs, geom_num, NULL);
	    Zoltan_Set_Geom_Multi_Fn(zs, geom_multi, NULL);
	    break;
#if 0 * ENABLE_REFTREE
	case LB_REFTREE:
	    /* Note: using (index,generation) pairs to identify elements */
	    Zoltan_Set_Param(zs, "NUM_GID_ENTRIES", "2");
	    Zoltan_Set_Param(zs, "NUM_LID_ENTRIES", "2");
	    /*Zoltan_Set_Param(zs, "REFTREE_INITPATH", "SIERPINSKI");*/
	    Zoltan_Set_Param(zs, "REFTREE_INITPATH", "CONNECTED");
	    break;
#endif	/* ENABLE_REFTREE */
	/*case LB_GRAPH:*/
	case LB_HYPERGRAPH:
	    /* valid values for coarsening method:
	     *	"AGG" (good),
	     *	"IPM" (good, default),
	     *	"L-IPM" (faster),
	     *	"A-IPM" (compromise) */
	    Zoltan_Set_Param(zs, "PHG_COARSENING_METHOD",
				phg_coarsening_methods[phg_coarsening_method]);
	    Zoltan_Set_Param(zs, "PHG_REFINEMENT_METHOD",
				phg_refinement_methods[phg_refinement_method]);
	    sprintf(str, "%lf", (double)phg_refinement_quality);
	    Zoltan_Set_Param(zs, "PHG_REFINEMENT_QUALITY", str);
	    sprintf(str, "%d", phg_check_hypergraph);
	    Zoltan_Set_Param(zs, "CHECK_HYPERGRAPH", str);
	    /* Mandatory query functions */
	    Zoltan_Set_HG_Size_CS_Fn(zs, HG_size_CS, NULL);
	    Zoltan_Set_HG_CS_Fn(zs, HG_CS, NULL);
	    break;
#if 0
	case LB_HIER:
	    break;
#endif
	case LB_BLOCK:
	case LB_RANDOM:
	    break;
    }

    rc = Zoltan_LB_Partition(zs, &changed, &ngids, &nlids,
			     &nimport, &import_gids, &import_lids,
			     &import_procs, &import_parts,
			     &nexport,
			     &export_gids, &export_lids,
			     &export_procs, &export_parts);

    /* Note: Zoltan may return inconsistent rc/changed flags if a process has
     *	     no element (rc < 0 and changed == 0 on that process). Below is a
     *	     quick workaround.
     * TODO: construct a communicator which only contains procs on which
     *	     g->nleaf > 0 */
    {
	int a[2], b[2];

	if (rc > 0)
	    rc = 0;		/* ignore non fatal error/warning */
	if (g->nleaf == 0)
	    rc = 0;

	a[0] = (rc == 0 ? 0 : 1);
	a[1] = (changed == 0 ? 0 : 1);
	MPI_Allreduce(&a, &b, 2, MPI_INT, MPI_MAX, g->comm);
	if ((rc = -b[0]) < 0)
	    return FALSE;	/* grid unchanged if error */
	changed = b[1];
    }

    if (rc < 0)
	return FALSE;

    assert(nexport == g->nleaf);

    for (i = 0; i < nexport; i++)
	g->elems[export_lids[i]]->mark = export_parts[i];

    Zoltan_LB_Free_Part(&import_gids, &import_lids,
			&import_procs, &import_parts);
    Zoltan_LB_Free_Part(&export_gids, &export_lids,
			&export_procs, &export_parts);
    Zoltan_Destroy(&zs);

    return changed;
}

#endif	/* USE_ZOLTAN */
#endif	/* USE_MPI */
