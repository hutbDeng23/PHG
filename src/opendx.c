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

/* IBM Data Explorer output functions
 *
 * $Id: opendx.c,v 1.41 2021/12/06 06:47:41 zlb Exp $
 */

#define SEPARATE_DATA_FILE	0	/* 1 ==> separate data file */

#include "phg.h"
#include "phg/io.h"

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>	/* PATH_MAX */
#include <ctype.h>

/* buffer for broadcasting filename */
static const char *filename;
/* buffer for broadcasting list of DOF ids */
static unsigned short *dof_ids = NULL;	/* last element for filename length */
static unsigned short dof_ids_len = 0;

static unsigned char *buffer = NULL;
static size_t buffer_pos, buffer_size;

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

#if 0
    while (*s != '\0' && p - name < sizeof(name) - 1) {
	*p = isspace(*(BYTE *)s) ? '_' : *s;
	p++;
	s++;
    }
    *p = '\0';
#else
    strcpy(p, s);
#endif

    if (n == 0)
	sprintf(name0, "%s", name);
    else
	sprintf(name0, "%s_%d", name, n);

    return name0;
}


#define DX_FLOAT float

#define AppendData(buf, size, count) {			\
    size_t size0 = (size) * (count);			\
    assert(buffer_pos + size0 <= buffer_size);	\
    memcpy(buffer + buffer_pos, buf, size0);	\
    buffer_pos += size0;				\
}

void
phgExportDX_(GRID *g)
/* Output grid in VTK unstructured grid format. */
{
    int k, a, dim;
    size_t size;
    INT i;
    INT nleaf, nvert;
    ELEMENT *s;
    DOF *dof;
    FLOAT *values = NULL;
    static char fn[PATH_MAX];
    char *endian;
#if SEPARATE_DATA_FILE
    FILE *fp;
# define PRINT(g, fp)	phgFPrintf(g->comm, fp, 
#else
# define PRINT(g, fp)	phgIOPrintf(g->comm, 
#endif

    ParallelFunction(phgExportDX_, TRUE);

    if (g == NULL || g->comm == MPI_COMM_NULL)
	return;

    nleaf = g->nleaf;
    nvert = g->nvert;

    assert((g->flags & VERT_FLAG));
    assert(g->rank != 0 || sizeof(a) == 4);	/* to be fixed later */

#if USE_MPI
    if (g->rank == 0) {
	memcpy(fn, filename, a = strlen(filename) + 1);
	dof_ids[dof_ids_len - 1] = a ;
    }
    filename = fn;
    MPI_Bcast(&dof_ids_len, sizeof(dof_ids_len), MPI_BYTE, 0, g->comm);
    MPI_Bcast(dof_ids, sizeof(*dof_ids) * dof_ids_len, MPI_BYTE, 0, g->comm);
    a = dof_ids[dof_ids_len - 1];
    MPI_Bcast((void *)filename, a, MPI_BYTE, 0, g->comm);
#else
    strcpy(fn, filename);
    filename = fn;
#endif

#if SEPARATE_DATA_FILE
    /* open header file */
    if ((fp = phgFOpen(g, filename, "w+t")) == NULL)
	    phgError(1, "cannot create  file \"%s\".\n", filename);
    strcat(filename, ".bin");
#endif	/* SEPARATE_DATA_FILE */
    filename = phgIOAddExtension(filename, "dx");
    if (!phgIOOpen(g->comm, (char *)filename)) {
	phgError(0, "cannot create  file \"%s\".\n", filename);
	return;
    }
    if (g->rank == 0)
	phgInfo(1, "creating DX file \"%s\".\n", filename);

    buffer_size = nvert * sizeof(DX_FLOAT) * 3;
    if ((size = nleaf * 4 * sizeof(a)) > buffer_size)
	buffer_size = size;
    /* compute size of DOF data */
    for (a = 0; a < dof_ids[0]; a++) {
	dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);
        dim = DofDim(dof);
	if (dof->type == DOF_P0 || dof->type == DOF_DG0) {
	    size = nleaf * dim * sizeof(DX_FLOAT);
	}
	else {
	    size = nvert * dim * sizeof(DX_FLOAT);
	}
	if (buffer_size < size)
	    buffer_size = size;
    }
    buffer = phgAlloc(buffer_size);

#if SEPARATE_DATA_FILE
    size = 0;
    /* strip path from filename */
    while ((endian = strchr(filename, '/')) != NULL)
	filename = endian + 1;
#endif

#if PHG_BIGENDIAN
    endian = "msb";
#else
    endian = "lsb";
#endif

    /* Comments */
    PRINT(g, fp) "# Data Explorer file created by PHG %s.\n", PHG_VERSION);

    /* vertices */
    PRINT(g, fp) "# vertices\nobject 1 class array type %s rank 1 shape 3 "
		      "items %"dFMT" %s binary\n",
		      sizeof(DX_FLOAT) == 8 ? "double" : "float",
		      g->nvert_global, endian);
#if SEPARATE_DATA_FILE
    PRINT(g, fp) "data file %s,%"dFMT"\n", filename, (INT)size);
#else
    PRINT(g, fp) "data follows\n");
#endif
    buffer_pos = 0;
    for (i = 0; i < nvert; i++) {
	for (a = 0; a < 3; a++) {
	    DX_FLOAT x = g->verts[i][a];
	    AppendData(&x, sizeof(x), 1);
	}
    }
#if USE_MPI
    phgIOWrite(g->comm, buffer, sizeof(DX_FLOAT)*3, nvert, g->nvert_global,
		g->types_vert, g->L2Gmap_vert, TRUE);
#else
    phgIOWrite(g->comm, buffer, sizeof(DX_FLOAT)*3, nvert, g->nvert_global,
		NULL, NULL, FALSE);
#endif
#if SEPARATE_DATA_FILE
    size += sizeof(DX_FLOAT) * 3 * g->nvert_global;
#endif
    PRINT(g, fp) "attribute \"dep\" string \"positions\"\n");
    
    /* tetrahedra */
    PRINT(g, fp) "# elements\nobject 2 class array type int rank 1 shape 4 "
		      "items %"dFMT" %s binary\n", g->nelem_global, endian);
#if SEPARATE_DATA_FILE
    PRINT(g, fp) "data file %s,%"dFMT"\n", filename, (INT)size);
#else
    PRINT(g, fp) "data follows\n");
#endif
    buffer_pos = 0;
    /* elements */
    ForAllElements(g, s) {
	for (i = 0; i < 4; i++) {
	    a = GlobalVertex(g, s->verts[i]);
	    AppendData(&a, sizeof(a), 1);
	}
    }
    phgIOWrite(g->comm, buffer, sizeof(a) * 4, nleaf, g->nelem_global,
							NULL, NULL, FALSE); 
#if SEPARATE_DATA_FILE
    size += sizeof(a) * 4 * g->nelem_global;
#endif
    PRINT(g, fp) "attribute \"element type\" string \"tetrahedra\"\n");
    PRINT(g, fp) "attribute \"ref\" string \"positions\"\n");

    for (a = 0; a < dof_ids[0]; a++) {
	size_t count, nglobal;
	BTYPE *flags;
	INT *map;
	DOF *dof0;
	BOOLEAN map_flag;
	dof = dof0 = phgDofIndex2Pointer(g, dof_ids[a + 1]);
	if (g->rank == 0)
	    phgInfo(2, "writing DOF %d: \"%s\"\n", a, dof->name);
	if (dof->type != DOF_P0 && dof->type != DOF_DG0) {
	    if (dof->type != DOF_P1)
		dof = phgDofCopy(dof0, NULL, DOF_P1, dof0->name);
	    count = nvert;
	    nglobal = g->nvert_global;
#if USE_MPI
	    flags = g->types_vert;
	    map = g->L2Gmap_vert;
	    map_flag = TRUE;
#else
	    flags = NULL;
	    map = NULL;
	    map_flag = FALSE;
#endif
	}
	else {
	    count = nleaf;
	    nglobal = g->nelem_global;
	    flags = NULL;
	    map = NULL;
	    map_flag = FALSE;
	}
	dim = dof->dim;
	if (dim == 1)
	    PRINT(g, fp) "#\nobject %d class array type %s rank 0 "
			      "items %"dFMT" %s binary\n", a + 3,
			      sizeof(DX_FLOAT) == 8 ? "double" : "float",
			      nglobal, endian);
	else
	    PRINT(g, fp) "#\nobject %d class array type %s rank 1 shape %d "
			      "items %"dFMT" %s binary\n", a + 3,
			      sizeof(DX_FLOAT) == 8 ? "double" : "float",
			      dim, nglobal, endian);
#if SEPARATE_DATA_FILE
	PRINT(g, fp) "data file %s,%"dFMT"\n", filename, (INT)size);
#else
	PRINT(g, fp) "data follows\n");
#endif
	/* collect data */
	if (dof->type != DOF_P0 && dof->type != DOF_DG0) {
#if 0
	    ForAllElements(g, s) {
		int j;
		for (k = 0; k < NVert; k++) {
		    values = DofVertexData(dof, s->verts[k]);
		    buffer_pos = s->verts[k] * dim * sizeof(DX_FLOAT);
		    for (j = 0; j < dim; j++) {
			DX_FLOAT v = *(values++);
			AppendData(&v, sizeof(v), 1);
		    }
		}
	    }
#else
	    buffer_pos = 0;
	    values = DofVertexData(dof, 0);
	    for (i = 0; i < dim * nvert; i++) {
		DX_FLOAT v = *(values++);
		AppendData(&v, sizeof(v), 1);
	    }
#endif
	}
	else {
	    buffer_pos = 0;
	    ForAllElements(g, s) {
		values = DofElementData(dof, s->index);
		for (k = 0; k < dim; k++) {
		    DX_FLOAT v = *(values++);
		    AppendData(&v, sizeof(v), 1);
		}
	    }
	}
	phgIOWrite(g->comm, buffer, sizeof(DX_FLOAT)*dim, count, nglobal,
							flags, map, map_flag);
#if SEPARATE_DATA_FILE
	size += sizeof(DX_FLOAT) * dim * g->nvert_global;
#endif
	if (dof->type != DOF_P0 && dof->type != DOF_DG0)
	    PRINT(g, fp) "attribute \"dep\" string \"positions\"\n");
	else
	    PRINT(g, fp) "attribute \"dep\" string \"connections\"\n");

	PRINT(g, fp) "#\nobject \"%s\" class field\n"
		     "component \"positions\" value 1\n"
		     "component \"connections\" value 2\n"
		     "component \"data\" value %d\n"
		     "attribute \"name\" string \"%s\"\n",
		     parse_name(g, a), a + 3, parse_name(g, a));

	if (dof != dof0)
	    phgDofFree(&dof);
    }

    /* save region_mark as cell data */
    PRINT(g, fp) "#\nobject %d class array type int rank 0 "
		 "items %"dFMT" %s binary\n", dof_ids[0] + 3,
			g->nelem_global, endian);
#if SEPARATE_DATA_FILE
    PRINT(g, fp) "data file %s,%"dFMT"\n", filename, (INT)size);
#else
    PRINT(g, fp) "data follows\n");
#endif
    buffer_pos = 0;
    ForAllElements(g, s) {
	a = s->region_mark;
	AppendData(&a, sizeof(a), 1);
    }
    phgIOWrite(g->comm, buffer, sizeof(a), nleaf, g->nelem_global,
							NULL, NULL, FALSE);
    PRINT(g, fp) "attribute \"dep\" string \"connections\"\n");
    PRINT(g, fp) "#\nobject \"region mark\" class field\n"
		     "component \"positions\" value 1\n"
		     "component \"connections\" value 2\n"
		     "component \"data\" value %d\n"
		     "attribute \"name\" string \"region mark\"\n",
			dof_ids[0] + 3);

    /* save submesh indices as cell data */
    PRINT(g, fp) "#\nobject %d class array type int rank 0 "
		 "items %"dFMT" %s binary\n", dof_ids[0] + 4,
			g->nelem_global, endian);
#if SEPARATE_DATA_FILE
    PRINT(g, fp) "data file %s,%"dFMT"\n", filename, (INT)size);
#else
    PRINT(g, fp) "data follows\n");
#endif
    buffer_pos = 0;
    ForAllElements(g, s) {
	a = g->rank;
	AppendData(&a, sizeof(a), 1);
    }
    phgIOWrite(g->comm, buffer, sizeof(a), nleaf, g->nelem_global,
							NULL, NULL, FALSE);
    PRINT(g, fp) "attribute \"dep\" string \"connections\"\n");
    PRINT(g, fp) "#\nobject \"submesh index\" class field\n"
		     "component \"positions\" value 1\n"
		     "component \"connections\" value 2\n"
		     "component \"data\" value %d\n"
		     "attribute \"name\" string \"submesh index\"\n",
			dof_ids[0] + 4);

    /* write group members */
    PRINT(g, fp) "#\nobject \"default\" class group\n");
    if (dof_ids[0] > 0) {
	/* group data */
	for (a = 0; a < dof_ids[0]; a++) {
	    dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);
	    PRINT(g, fp) "member \"%s\" value \"%s\"\n",
		parse_name(g, a), parse_name(g, a));
	}
    }
    PRINT(g, fp) "member \"region mark\" value \"region mark\"\n");
    PRINT(g, fp) "member \"submesh index\" value \"submesh index\"\n");

    PRINT(g, fp) "#\nend\n");

#if SEPARATE_DATA_FILE
    phgFClose(g, fp);
#endif

    phgIOClose();

    phgFree(buffer);
    buffer = NULL;

    return;
}

const char *
phgExportDXn(GRID *g, const char *fn, int ndof, DOF **dofs)
/* Outputs mesh in IBM OpenDX native format, returns the output filename.
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

    /* initialize dof_ids */
    dof_ids_len = ndof + 2;
    phgFree(dof_ids);
    dof_ids = phgAlloc(sizeof(*dof_ids) * dof_ids_len);

    for (i = 0, k = 0; i < ndof; i++) {
	if ((dof = dofs[i])->type == DOF_ANALYTIC) {
	    if (g->rank == 0) {
		phgWarning("DOF type 'DOF_ANALYTIC' is not supported "
			   "by phgExportDX(). \n");
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

    filename = fn;
    phgExportDX_(g);

    if (g->period != NULL) {
	phgFree(g->types_vert);
	g->types_vert = types_vert;
    }

    return filename;
}

const char *
phgExportDX(GRID *g, const char *fn, DOF *dof, ...)
{
    int ndof;
    DOF **dofs = phgAlloc(256 * sizeof(*dofs));
    va_list ap;

    va_start(ap, dof);
    for (ndof = 0; ndof < 256; ndof++) {
	if (dof == NULL)
	    break;
	dofs[ndof] = dof;
	dof = va_arg(ap, DOF *);
    }
    va_end(ap);

    fn = phgExportDXn(g, fn, ndof, dofs);
    phgFree(dofs);

    return fn;
}
