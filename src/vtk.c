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

/* VTK output functions
 *
 * $Id: vtk.c,v 1.61 2021/12/06 07:13:27 zlb Exp $
 */

#define USE_VECTORS	1

#include "phg.h"
#include "phg/io.h"

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>	/* PATH_MAX */
#include <ctype.h>

static BOOLEAN vtk_nonlinear_flag = TRUE;
static BOOLEAN vtk_boundary_flag = TRUE;

enum {VTK_SINGLE = 0, VTK_DOUBLE = 1};
static const char *vtk_precision_list[] = {"single", "double", NULL};
static int vtk_precision = VTK_SINGLE;

#define VTK_PUT_FLOAT(a) {			\
    if (vtk_precision == VTK_SINGLE) {		\
	float x = (float)(a);			\
	phgIOBigEndianAppend(&x, sizeof(x), 1);	\
    }						\
    else {					\
	double x = (double)(a);			\
	phgIOBigEndianAppend(&x, sizeof(x), 1);	\
    }						\
}

/* buffer for broadcasting filename */
static const char *filename;
/* buffer for broadcasting list of DOF ids */
static unsigned short *dof_ids = NULL;	/* last element for filename length */
static unsigned short dof_ids_len = 0;

static char *
parse_name(GRID *g, int no)
{
    static char name0[128 + 16];
    char name[128];
    char *p = name;
    const char *s;
    int i, n;

    s = phgDofIndex2Pointer(g, dof_ids[no + 1])->name;

    /* check for duplicate names */
    for (i = 0, n = 0; i < no; i++)
	if (!strcmp(s, phgDofIndex2Pointer(g, dof_ids[i + 1])->name))
	    n++;

    while (*s != '\0' && p - name < sizeof(name) - 1) {
	*p = isspace(*(const BYTE *)s) ? '_' : *s;
	p++;
	s++;
    }
    *p = '\0';

    if (n == 0)
	sprintf(name0, "%s", name);
    else
	sprintf(name0, "%s_%d", name, n);

    return name0;
}

void
phgExportVTK_(GRID *g)
/* Output grid in VTK unstructured grid format. */
{
    int k, a, dim;
    BYTE face;
    size_t size, vtk_float_size;
    INT i;
    INT nleaf, nvert, nedge, nbface, nbface_global;
    ELEMENT *e;
    DOF *dof;
    BOOLEAN has_cell_data = FALSE, has_point_data = FALSE;
    FLOAT *values;
    static char fn[PATH_MAX];
    unsigned char c, output_order;
    size_t buffer_size;
    void *buffer;

    ParallelFunction(phgExportVTK_, TRUE);

    if (g == NULL || g->comm == MPI_COMM_NULL)
	return;

    vtk_float_size = (vtk_precision == VTK_SINGLE ?
				sizeof(float) : sizeof(double));

    nleaf = g->nleaf;
    nvert = g->nvert;
    nedge = g->nedge;

    assert((g->flags & VERT_FLAG));
    assert(g->rank != 0 || sizeof(a) == 4);	/* to be fixed later */

#if USE_MPI
    if (phgMasterSlave) {
	if (g->rank == 0) {
	    memcpy(fn, filename, a = strlen(filename) + 1);
	    dof_ids[dof_ids_len - 1] = a ;
	}
	filename = fn;
	MPI_Bcast(&dof_ids_len, sizeof(dof_ids_len), MPI_BYTE, 0, g->comm);
	MPI_Bcast(dof_ids, sizeof(*dof_ids)*dof_ids_len, MPI_BYTE, 0, g->comm);
	a = dof_ids[dof_ids_len - 1];
	MPI_Bcast((void *)filename, a, MPI_BYTE, 0, g->comm);
    }
    else
#else
    {
	strcpy(fn, filename);
	filename = fn;
    }
#endif

    filename = phgIOAddExtension(filename, "vtk");
    if (!phgIOOpen(g->comm, (char *)filename)) {
	phgError(0, "cannot create VTK file \"%s\".\n", filename);
	return;
    }
    if (g->rank == 0)
	phgInfo(1, "creating VTK file \"%s\".\n", filename);

    /* determine output order */
    if (vtk_nonlinear_flag == FALSE) {
	output_order = 1;
    }
    else {
	output_order = 0;
	for (a = 0; a < dof_ids[0]; a++) {
	    dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);
	    if (DofIsHP(dof)) {
		if (output_order < dof->hp->max_order)
		    output_order = dof->hp->max_order;
	    }
	    else {
		if (output_order < dof->type->order)
		    output_order = dof->type->order;
	    }
	}
    }

    if (output_order > 1 && g->period != NULL) {
	static BOOLEAN warned = FALSE;
	if (warned == FALSE && g->rank == 0)
	    phgWarning("nonlinear element for periodic boundaries unimplemented yet, using linear element instead.\n");
	warned = TRUE;
	output_order = 1;
    }

    nbface = nbface_global = 0;
    if (vtk_boundary_flag) {
	/* count number of boundary faces */
	ForAllFaces(g, e, face) {
	    if (!(e->bound_type[face] & BDRY_MASK))
		continue;
	    assert((g->types_face[e->faces[face]] & OWNER));	/* FIXME! */
	    nbface++;
	}
#if USE_MPI
	MPI_Allreduce(&nbface, &nbface_global, 1,PHG_MPI_INT, MPI_SUM, g->comm);
#else	/* USE_MPI */
	nbface_global = nbface;
#endif	/* USE_MPI */
    }

    if (output_order <= 1) {
	/* use linear tetrahedra */
	output_order = 1;
	buffer_size = nvert * vtk_float_size * 3;
	if ((size = nleaf * 5 * sizeof(a)) > buffer_size)
	    buffer_size = size;
	if ((size = nbface * 4 * sizeof(a)) > buffer_size)
	    buffer_size = size;
    }
    else {
	/* use quadratic tetrahedra */
	output_order = 2;
	buffer_size = (nvert > nedge ? nvert : nedge) * vtk_float_size * 3;
	if ((size = nleaf * 11 * sizeof(a)) > buffer_size)
	    buffer_size = size;
	if ((size = nbface * 7 * sizeof(a)) > buffer_size)
	    buffer_size = size;
    }

    /* compute size of DOF data */
    for (a = 0; a < dof_ids[0]; a++) {
	dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);
        dim = DofDim(dof);
	dim = dim > 4 ? 4 : dim;
	if (dof->type == DOF_P0 || dof->type == DOF_DG0) {
	    size = (nleaf + nbface) * dim * vtk_float_size;
	    has_cell_data = TRUE;
	}
	else {
	    if (output_order <= 1)
		size = nvert * dim * vtk_float_size;
	    else
		size = (nvert + nedge) * dim * vtk_float_size;
	    has_point_data = TRUE;
	    if (DofIsHP(dof)) {
		has_cell_data = TRUE;	/* hp orders are saved as cell data */
		if (buffer_size < size)
		    buffer_size = size;
		size = (nleaf + nbface) * sizeof(unsigned char);
	    }
	}
	if (buffer_size < size)
	    buffer_size = size;
    }

    buffer = phgIOAllocBuffer(buffer_size);

    /* File header and vertices */
    phgIOPrintf(g->comm, "# vtk DataFile Version 2.0\n"
		"Tetrahedral grid created by PHG %s.\n"
		"BINARY\nDATASET UNSTRUCTURED_GRID\n"
		"POINTS %"dFMT" %s\n", PHG_VERSION,
		g->nvert_global + (output_order > 1 ? g->nedge_global : 0),
		vtk_precision == VTK_DOUBLE ? "double" : "float");
    phgIOResetBuffer();
    for (i = 0; i < nvert; i++) {
	for (a = 0; a < 3; a++) {
	    VTK_PUT_FLOAT(g->verts[i][a])
	}
    }
#if USE_MPI
    phgIOWrite(g->comm, buffer, vtk_float_size * 3, nvert, g->nvert_global,
		g->types_vert, g->L2Gmap_vert, TRUE);
#else
    phgIOWrite(g->comm, buffer, vtk_float_size * 3, nvert, g->nvert_global,
		NULL, NULL, FALSE);
#endif

    if (output_order > 1) {
	/* save list of nodes at edge centers */
	VEF_MAP *vef = phgDofSetupVEFMap(g, NULL, EDGE_FLAG);
	assert(vef != NULL || nedge == 0);
	phgIOResetBuffer();
	for (i = 0; i < nedge; i++) {
	    COORD *p0, *p1;
	    e = vef->Emap[i];
	    if (e == NULL) {
		/* unreferenced edge */
		VTK_PUT_FLOAT(0.)
		VTK_PUT_FLOAT(0.)
		VTK_PUT_FLOAT(0.)
		continue;
	    }
	    k = vef->Eind[i];
	    p0 = g->verts + e->verts[GetEdgeVertex(k, 0)];
	    p1 = g->verts + e->verts[GetEdgeVertex(k, 1)];
	    for (a = 0; a < 3; a++) {
		VTK_PUT_FLOAT(((*p0)[a] + (*p1)[a]) * 0.5)
	    }
	}
#if USE_MPI
	phgIOWrite(g->comm, buffer, vtk_float_size * 3, nedge, g->nedge_global,
		   g->types_edge, g->L2Gmap_edge, TRUE);
#else
	phgIOWrite(g->comm, buffer, vtk_float_size * 3, nedge, g->nedge_global,
		   NULL, NULL, FALSE);
#endif
	phgDofFreeVEFMap(&vef);
    }

    /* element cells */
    phgIOPrintf(g->comm, "\nCELLS %"dFMT" %"dFMT"\n", g->nelem_global + nbface_global,
		g->nelem_global * (output_order > 1 ? 11 : 5) +
		nbface_global * (output_order > 1 ? 7 : 4));

    /* tetrahdra */
    phgIOResetBuffer();
    ForAllElements(g, e) {
	/* Note: VTK_QUADRATIC_TETRA has type 24 with nodes ordered as:
	 * 	v0, v1, v2, v3, e01, e12, e02, e03, e13, e23
	 * (see /usr/include/vtk/vtkQuadraticTetra.h) */
	static int edge_order[] = {0, 3, 1, 2, 4, 5};
	a = (output_order > 1 ? 10 : 4);
	phgIOBigEndianAppend(&a, sizeof(a), 1);
	/* vertices */
	for (i = 0; i < 4; i++) {
	    a = GlobalVertex(g, e->verts[i]);
	    phgIOBigEndianAppend(&a, sizeof(a), 1);
	}
	if (output_order <= 1)
	    continue;
	/* edge centers */
	for (i = 0; i < 6; i++) {
	    a = GlobalEdge(g, e->edges[edge_order[i]]) + g->nvert_global;
	    phgIOBigEndianAppend(&a, sizeof(a), 1);
	}
    }
    phgIOWrite(g->comm, buffer, sizeof(a) * (output_order > 1 ? 11 : 5), nleaf,
		g->nelem_global, NULL, NULL, FALSE); 

    if (nbface_global > 0) {
	/* boundary faces */
	phgIOResetBuffer();
	ForAllFaces(g, e, face) {
	    /* Note: VTK_QUADRATIC_TRIANGLE has type 22 with nodes ordered as:
	     *		v0, v1, v2, e01, e12, e02
	     * (see /usr/include/vtk/vtkQuadraticTriangle.h) */
	    static int edges[][2] = {{0,1}, {1,2}, {0,2}};
	    int v[3];
	    if (!(e->bound_type[face] & BDRY_MASK))
		continue;
	    assert((g->types_face[e->faces[face]] & OWNER));	/* FIXME! */
	    a = (output_order > 1 ? 6 : 3);
	    phgIOBigEndianAppend(&a, sizeof(a), 1);
	    /* vertices */
	    for (i = 0; i < 3; i++) {
		a = GlobalVertex(g, e->verts[GetFaceVertex(face, i)]);
		phgIOBigEndianAppend(&a, sizeof(a), 1);
	    }
	    if (output_order <= 1)
		continue;
	    /* edge centers */
	    v[0] = GetFaceVertex(face, 0);
	    v[1] = GetFaceVertex(face, 1);
	    v[2] = GetFaceVertex(face, 2);
	    for (i = 0; i < 3; i++) {
		int edge_no = GetEdgeNo(v[edges[i][0]],v[edges[i][1]]);
		a = GlobalEdge(g, e->edges[edge_no]) + g->nvert_global;
		phgIOBigEndianAppend(&a, sizeof(a), 1);
	    }
	}
	phgIOWrite(g->comm, buffer, sizeof(a) * (output_order > 1 ? 7 : 4),
			nbface, nbface_global, NULL, NULL, FALSE); 
    }

    phgIOPrintf(g->comm, "\nCELL_TYPES %"dFMT"\n",
				g->nelem_global + nbface_global);

    /* types for tetrahedra */
    phgIOResetBuffer();
    a = (output_order > 1 ? /* VTK_QUADRATIC_TETRA = */24 : /*VTK_TETRA = */10);
    for (i = 0; i < nleaf; i++)
	phgIOBigEndianAppend(&a, sizeof(a), 1);
    phgIOWrite(g->comm, buffer, sizeof(a), nleaf, g->nelem_global,
							NULL, NULL, FALSE);

    /* types for triangles */
    if (nbface_global > 0) {
	phgIOResetBuffer();
	a = (output_order > 1 ? 22 : 5);
	for (i = 0; i < nbface; i++)
	    phgIOBigEndianAppend(&a, sizeof(a), 1);
	phgIOWrite(g->comm, buffer, sizeof(a), nbface, nbface_global,
						NULL, NULL, FALSE);
    }

    if (has_point_data) {	/* point data */
	phgIOPrintf(g->comm, "\nPOINT_DATA %"dFMT"\n",
		g->nvert_global + (output_order > 1 ? g->nedge_global : 0));
	for (a = 0; a < dof_ids[0]; a++) {
	    DOF *dof0;
	    dof = dof0 = phgDofIndex2Pointer(g, dof_ids[a + 1]);
	    if (dof->type == DOF_P0 || dof->type == DOF_DG0)
		continue;
	    if (dof->type != DOF_Pn[output_order])
		dof = phgDofCopy(dof0, NULL, DOF_Pn[output_order], dof0->name);
	    if (g->rank == 0)
		phgInfo(2, "writing DOF %d: \"%s\"\n", a, dof->name);
	    dim = dof->dim;
	    if (dim > 4) {
		if (g->rank == 0)
		    phgWarning("only saving the first 4 components of \"%s\"\n",
				dof->name);
		dim = 4;
	    }
	    if (USE_VECTORS && dim == Dim)
	        phgIOPrintf(g->comm, "\nVECTORS %s %s\n", parse_name(g, a),
			    vtk_precision == VTK_DOUBLE ? "double" : "float");
	    else
	        phgIOPrintf(g->comm, "\nSCALARS %s %s %d\n"
			    "LOOKUP_TABLE default\n",
			    parse_name(g, a),
			    vtk_precision == VTK_DOUBLE ? "double" : "float",
			    dim);
	    /* collect vertex data */
	    phgIOResetBuffer();
	    values = DofVertexData(dof, 0);
	    for (i = 0; i < dim * nvert; i++) {
		VTK_PUT_FLOAT(*(values++))
	    }
#if USE_MPI
	    phgIOWrite(g->comm, buffer, vtk_float_size * dim, nvert,
			g->nvert_global, g->types_vert, g->L2Gmap_vert, TRUE); 
#else
	    phgIOWrite(g->comm, buffer, vtk_float_size * dim, nvert,
			g->nvert_global, NULL, NULL, FALSE); 
#endif
	    if (output_order > 1) {
		/* collect edge data */
		phgIOResetBuffer();
		values = DofEdgeData(dof, 0);
		for (i = 0; i < dim * nedge; i++) {
		    VTK_PUT_FLOAT(*(values++))
		}
#if USE_MPI
		phgIOWrite(g->comm, buffer, vtk_float_size * dim, nedge,
			g->nedge_global, g->types_edge, g->L2Gmap_edge, TRUE); 
#else
		phgIOWrite(g->comm, buffer, vtk_float_size * dim, nedge,
			g->nedge_global, NULL, NULL, FALSE); 
#endif
	    }
	    if (dof != dof0)
		phgDofFree(&dof);
	}
    }

    phgIOPrintf(g->comm, "\nCELL_DATA %"dFMT"\n",
				g->nelem_global + nbface_global);

    if (has_cell_data) {
	for (a = 0; a < dof_ids[0]; a++) {
	    dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);

	    if (DofIsHP(dof)) {
		/* variable order DOF, save orders as cell data */
		phgIOPrintf(g->comm, "\nSCALARS %s_orders unsigned_char 1\n"
			   "LOOKUP_TABLE default\n", parse_name(g, a));
		/* tetrahedra */
		phgIOResetBuffer();
		ForAllElements(g, e) {
		    c = dof->hp->elem_order[e->index];
		    phgIOBigEndianAppend(&c, sizeof(c), 1);
		}
		phgIOWrite(g->comm, buffer, sizeof(c), nleaf, g->nelem_global,
				NULL, NULL, FALSE); 
		if (nbface_global > 0) {
		    /* boundary faces */
		    phgIOResetBuffer();
		    ForAllFaces(g, e, face) {
			if (!(e->bound_type[face] & BDRY_MASK))
			    continue;
			assert((g->types_face[e->faces[face]] & OWNER));
			c = dof->hp->face_order[e->faces[face]];
			phgIOBigEndianAppend(&c, sizeof(c), 1);
		    }
		    phgIOWrite(g->comm, buffer, sizeof(c),
				nbface, nbface_global, NULL, NULL, FALSE); 
		}
		continue;
	    }

	    if (dof->type != DOF_P0 && dof->type != DOF_DG0)
		continue;
	    if (g->rank == 0)
		phgInfo(2, "writing DOF %d: \"%s\"\n", a, dof->name);
	    dim = dof->dim;
	    if (dim > 4) {
		if (g->rank == 0)
		    phgWarning("only saving the first 4 components of \"%s\"\n",
				dof->name);
		dim = 4;
	    }
	    if (USE_VECTORS && dim == Dim)
	 	phgIOPrintf(g->comm, "\nVECTORS %s %s\n", parse_name(g, a),
			    vtk_precision == VTK_DOUBLE ? "double" : "float");
	    else
	 	phgIOPrintf(g->comm, "\nSCALARS %s %s %d\n"
			    "LOOKUP_TABLE default\n",
			    parse_name(g, a),
			    vtk_precision == VTK_DOUBLE ? "double" : "float",
			    dim);
	    phgIOResetBuffer();
	    ForAllElements(g, e) {
		values = DofElementData(dof, e->index);
		for (k = 0; k < dim; k++) {
		    VTK_PUT_FLOAT(*(values++))
		}
	    }
	    phgIOWrite(g->comm, buffer, vtk_float_size * dim, nleaf,
					g->nelem_global, NULL, NULL, FALSE); 
	    if (nbface_global > 0) {
		phgIOResetBuffer();
		ForAllFaces(g, e, face) {
		    if (!(e->bound_type[face] & BDRY_MASK))
			continue;
		    assert((g->types_face[e->faces[face]] & OWNER));
		    values = DofElementData(dof, e->index);
		    for (k = 0; k < dim; k++) {
			VTK_PUT_FLOAT(*(values++))
		    }
		}
		phgIOWrite(g->comm, buffer, vtk_float_size * dim, nbface,
					nbface_global, NULL, NULL, FALSE); 
	    }
	}
    }

    if (nbface_global > 0) {
	/* save boundary types as cell data */
	phgIOPrintf(g->comm, "\nSCALARS boundary_mask int 1\n"
			     "LOOKUP_TABLE default\n");
	phgIOResetBuffer();
	a = 0;
	ForAllElements(g, e)
	    phgIOBigEndianAppend(&a, sizeof(a), 1);
	phgIOWrite(g->comm, buffer, sizeof(a), nleaf, g->nelem_global,
						NULL, NULL, FALSE);
	phgIOResetBuffer();
	ForAllFaces(g, e, face) {
	    if (!((a = e->bound_type[face]) & BDRY_MASK))
		continue;
	    assert((g->types_face[e->faces[face]] & OWNER));
	    phgIOBigEndianAppend(&a, sizeof(a), 1);
	}
	phgIOWrite(g->comm, buffer, sizeof(a), nbface, nbface_global,
						NULL, NULL, FALSE);
    }

    /* save regional mark as cell data */
    phgIOPrintf(g->comm, "\nSCALARS region_mark int 1\nLOOKUP_TABLE default\n");
    phgIOResetBuffer();
    ForAllElements(g, e) {
	a = e->region_mark;
	phgIOBigEndianAppend(&a, sizeof(a), 1);
    }
    phgIOWrite(g->comm, buffer, sizeof(a), nleaf, g->nelem_global,
							NULL, NULL, FALSE);
    if (nbface_global > 0) {
	phgIOResetBuffer();
	ForAllFaces(g, e, face) {
	    if (!((a = e->bound_type[face]) & BDRY_MASK))
		continue;
	    assert((g->types_face[e->faces[face]] & OWNER));
	    a = e->region_mark;
	    phgIOBigEndianAppend(&a, sizeof(a), 1);
	}
	phgIOWrite(g->comm, buffer, sizeof(a), nbface, nbface_global,
							NULL, NULL, FALSE);
    }

    /* save submesh indices as cell data */
    phgIOPrintf(g->comm, "\nSCALARS submesh_no int 1\nLOOKUP_TABLE default\n");
    a = g->rank;
    phgIOResetBuffer();
    for (i = 0; i < nleaf; i++)
	phgIOBigEndianAppend(&a, sizeof(a), 1);
    phgIOWrite(g->comm, buffer, sizeof(a), nleaf, g->nelem_global,
							NULL, NULL, FALSE);
    if (nbface_global > 0) {
	phgIOResetBuffer();
	for (i = 0; i < nbface; i++)
	    phgIOBigEndianAppend(&a, sizeof(a), 1);
	phgIOWrite(g->comm, buffer, sizeof(a), nbface, nbface_global,
						NULL, NULL, FALSE);
    }

    phgIOClose();

    phgIOFreeBuffer();

    return;
}

const char *
phgExportVTKn(GRID *g, const char *fn, int ndof, DOF **dofs)
/* Outputs mesh in VTK format, returns the output filename.
 * The variable part of the arguments is a NULL terminated list of DOFs
 * to save to the VTK file */
{
    int i, j, k;
    DOF *dof;
    BTYPE *types_vert = NULL;

    FreeAtExit(dof_ids);

    assert(fn != NULL && Dim == 3);
    if (g == NULL)
	return NULL;

    /* output fn */
    filename = fn;
    /* initialize dof_ids */
    dof_ids_len = ndof + 2;
    phgFree(dof_ids);
    dof_ids = phgAlloc(sizeof(*dof_ids) * dof_ids_len);

    /* list of DOFs */
    for (i = 0, k = 0; i < ndof; i++) {
	if ((dof = dofs[i])->type == DOF_ANALYTIC) {
	    if (g->rank == 0) {
		phgWarning("DOF type 'DOF_ANALYTIC' is not supported by "
			   "phgExportVTK(). \n");
		phgWarning("DOF '%s' will not be saved.\n", dof->name);
	    }
	    continue;
	}
	if (dof->g != g)
	    phgError(1, "%s: DOF \"%s\" doesn't belong to the grid.\n",
			__func__, dof->name);
	/* test for missing NULL terminator */
	if (!phgDofIsValid(g, dof)) {
	    if (g->rank == 0)
		phgWarning("phgExportVTK: missing NULL terminator"
			   " in argument list.\n");
	    break;
	}
	if ((j = phgDofPointer2Index(dof)) < 0)
	    break;
	dof_ids[++k] = j;
    }
    dof_ids[0] = k;

    if (g->period != NULL) {
	/* reconstruct types_vert[] without periodicity */
	BYTE	flags;
	PERIOD	*period;
	INT	nvert_owned, *owner_index_vert;
	int	*owner_rank_vert;

	period = g->period;
	g->period = NULL;
	types_vert = g->types_vert;
	g->types_vert = NULL;
	owner_index_vert = g->owner_index_vert;
	g->owner_index_vert = NULL;
	owner_rank_vert = g->owner_rank_vert;
	g->owner_rank_vert = NULL;
	flags = g->flags;
	g->flags = VERT_FLAG;
	nvert_owned = g->nvert_owned;

	phgUpdateBoundaryTypes(g);

	g->period = period;
	phgFree(g->owner_index_vert);
	g->owner_index_vert = owner_index_vert;
	phgFree(g->owner_rank_vert);
	g->owner_rank_vert = owner_rank_vert;
	g->flags = flags;
	g->nvert_owned = nvert_owned;
    }

    phgExportVTK_(g);

    if (g->period != NULL) {
	phgFree(g->types_vert);
	g->types_vert = types_vert;
    }

    return filename;
}

const char *
phgExportVTK(GRID *g, const char *fn, DOF *dof, ...)
{
    int ndof;
    DOF **dofs;
    va_list ap;

    assert(g != NULL || phgInitialized == FALSE);

    if (g == NULL) {
	/* register options for VTK export */
	phgOptionsRegisterTitle("\nVTK export options:", "\n", "vtk");
	phgOptionsRegisterKeyword("-vtk_precision", "VTK data precision",
				  vtk_precision_list, &vtk_precision);
	phgOptionsRegisterNoArg("-vtk_nonlinear_elements",
		"Use VTK's non linear elements for high order elements",
		&vtk_nonlinear_flag);
	phgOptionsRegisterNoArg("-vtk_boundary_faces",
		"Output boundary faces, see BDRY_* in phg.h for the boundary "
		"masks", &vtk_boundary_flag);

	return NULL;
    }

    dofs = phgAlloc(256 * sizeof(*dofs));

    va_start(ap, dof);
    for (ndof = 0; ndof < 256; ndof++) {
	if (dof == NULL)
	    break;
	dofs[ndof] = dof;
	dof = va_arg(ap, DOF *);
    }
    va_end(ap);

    fn = phgExportVTKn(g, fn, ndof, dofs);
    phgFree(dofs);

    return fn;
}
