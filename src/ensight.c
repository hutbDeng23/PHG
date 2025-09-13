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

/* Ensight output functions
 *
 * $Id: ensight.c,v 1.13 2021/12/06 06:47:41 zlb Exp $
 */

#include "phg.h"
#include "phg/io.h"

#include <stdlib.h>
#include <string.h>
#include <strings.h>	/* bzero() */
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
static FLOAT time = 0;
static int tstep = 0;

#define ENSIGHT_FLOAT float

#define AppendData(buf, size, count) {			\
	size_t size0 = (size) * (count);		\
	assert(buffer_pos + size0 <= buffer_size);	\
	memcpy(buffer + buffer_pos, buf, size0);	\
	buffer_pos += size0;				\
    }
#define writeString(g, str) {				\
	char buffer[80];				\
	bzero(buffer, 80 * sizeof(char));		\
	strcpy(buffer, str);				\
	phgIOWriteRoot(g->comm, buffer, sizeof(char), 80);	\
    }
#define writeOneInt(g, i) {				\
	int _tmp_i = i;					\
	phgIOWriteRoot(g->comm, &_tmp_i, sizeof(int), 1);	\
    }

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
phgExportEnsight_(GRID *g)
/* Output grid in Ensight unstructured grid format. */
{
    int k, a, dim;
    size_t size;
    INT i;
    INT nleaf, nvert;
    ELEMENT *s;
    DOF *dof;
    FLOAT *values;
    static char fn[PATH_MAX];
    char case_file[PATH_MAX];
    char geo_file[PATH_MAX];
    char var_file[PATH_MAX];

    ParallelFunction(phgExportEnsight_, TRUE);

    if (g == NULL || g->comm == MPI_COMM_NULL)
	return;

    nleaf = g->nleaf;
    nvert = g->nvert;

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

    buffer_size = nvert * sizeof(ENSIGHT_FLOAT) * 3;
    if ((size = nleaf * 5 * sizeof(a)) > buffer_size)
	buffer_size = size;
    /* compute size of DOF data */
    for (a = 0; a < dof_ids[0]; a++) {
	dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);
        dim = DofDim(dof);
	assert(dim == 1 || dim == 3 || dim == 9);
	if (dof->type == DOF_P0 || dof->type == DOF_DG0) {
	    size = nleaf * dim * sizeof(ENSIGHT_FLOAT);
	}
	else {
	    size = nvert * dim * sizeof(ENSIGHT_FLOAT);
	}
	if (buffer_size < size)
	    buffer_size = size;
    }
    buffer = phgAlloc(buffer_size);

    /*****************/
    /* Geometry file */
    /*****************/
    sprintf(geo_file, "%s.%05d.geo", filename, tstep);
    if (g->rank == 0) phgInfo(1, "creating Ensight file \"%s\".\n", geo_file);
    phgIOOpen(g->comm, geo_file);

    /* Geometry file: description */
    writeString(g, "C Binary");
    writeString(g, "1st description");
    writeString(g, "2nd description");

    /* Geometry file: id options */
    writeString(g, "node id off");
    writeString(g, "element id off");

    /* Geometry file: coordinates */
    writeString(g, "coordinates");
    writeOneInt(g, g->nvert_global);	

    buffer_pos = 0;
    for (i = 0; i < nvert; i++) {
	for (a = 0; a < 3; a++) {
	    ENSIGHT_FLOAT x = g->verts[i][a];
	    AppendData(&x, sizeof(x), 1);
	}
    }
#if USE_MPI
    phgIOWrite(g->comm, buffer, sizeof(ENSIGHT_FLOAT)*3, nvert, g->nvert_global,
	       g->types_vert, g->L2Gmap_vert, TRUE);
#else
    phgIOWrite(g->comm, buffer, sizeof(ENSIGHT_FLOAT)*3, nvert, g->nvert_global,
	       NULL, NULL, FALSE);
#endif
    
    /* Geometry file: parts */
    writeString(g, "part 1");
    writeString(g, "description");
    writeString(g, "tetra4");
    writeOneInt(g, g->nelem_global);

    buffer_pos = 0;
    /* element vertices */
    ForAllElements(g, s) {
	for (i = 0; i < 4; i++) {
	    a = GlobalVertex(g, s->verts[i]) + 1;
	    AppendData(&a, sizeof(a), 1);
	}
    }
    phgIOWrite(g->comm, buffer, sizeof(a) * 4, nleaf, g->nelem_global,
							NULL, NULL, FALSE); 
    phgIOClose();


    /*****************/
    /* Variable file */
    /*****************/
    for (a = 0; a < dof_ids[0]; a++) {
	char *dof_name;

	dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);
	dof_name = parse_name(g, a);
	sprintf(var_file, "%s.%05d.%s", filename, tstep, dof_name);
	if (g->rank == 0) phgInfo(1, "creating Ensight file \"%s\".\n", var_file);
	phgIOOpen(g->comm, var_file);

	if (DofIsHP(dof)) {
	    phgError(0, "Unimplement output ensight for HP Dof.\n");
	}

	if (g->rank == 0)
	    phgInfo(2, "writing DOF %d: \"%s\"\n", a, dof->name);
	dim = DofDim(dof);

	if (dof->type == DOF_P0 || dof->type == DOF_DG0) {
	    /* Cell Data */
	    writeString(g, "cell data");
	    writeString(g, "part 1");
	    writeString(g, "tetra4");
	    buffer_pos = 0;
	    ForAllElements(g, s) {
		values = DofElementData(dof, s->index);
		for (k = 0; k < dim; k++) {
		    ENSIGHT_FLOAT v = *(values++);
		    AppendData(&v, sizeof(v), 1);
		}
	    }
	    phgIOWrite(g->comm, buffer, sizeof(ENSIGHT_FLOAT), nleaf,
		       g->nelem_global, NULL, NULL, FALSE); 
	} else {
	    /* Vertex Date */
	    DOF *dof0;

	    writeString(g, "node data");
	    dof = dof0 = phgDofIndex2Pointer(g, dof_ids[a + 1]);
	    if (dof->type != DOF_P1)
		dof = phgDofCopy(dof0, NULL, DOF_P1, dof0->name);

	    buffer_pos = 0;
	    values = DofVertexData(dof0, 0);
	    for (i = 0; i < dim * nvert; i++) {
		ENSIGHT_FLOAT v = *(values++);
		AppendData(&v, sizeof(v), 1);
	    }
#if USE_MPI
	    phgIOWrite(g->comm, buffer, sizeof(ENSIGHT_FLOAT) * dim, nvert,
		       g->nvert_global, g->types_vert, g->L2Gmap_vert, TRUE); 
#else
	    phgIOWrite(g->comm, buffer, sizeof(ENSIGHT_FLOAT) * dim, nvert,
		       g->nvert_global, NULL, NULL, FALSE); 
#endif

	    if (dof != dof0)
		phgDofFree(&dof);
	}
	
	phgIOClose();
    } /* end of variables */


    /* save regional mark as cell data */
    sprintf(var_file, "%s.%05d.region", filename, tstep);
    if (g->rank == 0) phgInfo(1, "creating Ensight file \"%s\".\n", var_file);
    phgIOOpen(g->comm, var_file);

    writeString(g, "region mask");
    writeString(g, "part 1");
    writeString(g, "tetra4");
    buffer_pos = 0;
    ForAllElements(g, s) {
	ENSIGHT_FLOAT v = s->region_mark;
	AppendData(&v, sizeof(v), 1);
    }
    phgIOWrite(g->comm, buffer, sizeof(ENSIGHT_FLOAT), nleaf, g->nelem_global,
							NULL, NULL, FALSE);
    phgIOClose();

    /* save submesh indices as cell data */
    sprintf(var_file, "%s.%05d.rank", filename, tstep);
    if (g->rank == 0) phgInfo(1, "creating Ensight file \"%s\".\n", var_file);
    phgIOOpen(g->comm, var_file);

    writeString(g, "submesh indicies");
    writeString(g, "part 1");
    writeString(g, "tetra4");
    buffer_pos = 0;
    ForAllElements(g, s) {
	ENSIGHT_FLOAT v = g->rank;
	AppendData(&v, sizeof(v), 1);
    }
    phgIOWrite(g->comm, buffer, sizeof(ENSIGHT_FLOAT), nleaf, g->nelem_global,
							NULL, NULL, FALSE);
    phgIOClose();

    /*************/
    /* Case file */
    /*************/
    sprintf(case_file, "%s.%05d.case", filename, tstep);
    if (g->rank == 0) phgInfo(1, "creating Ensight file \"%s\".\n", case_file);
    phgIOOpen(g->comm, case_file);
    phgIOPrintf(g->comm, "# Ensight file, format ensight6, created by PHG\n");

    /* Case file: Format section */
    phgIOPrintf(g->comm, "\nFORMAT\n");
    phgIOPrintf(g->comm, "type:       ensight\n");

    /* Case file: Geometry section */
    phgIOPrintf(g->comm, "\nGEOMETRY\n");
    phgIOPrintf(g->comm, "model:      %s.*****.%s\n", filename, "geo");
    
    /* Case file: Variable section */
    phgIOPrintf(g->comm, "\nVARIABLE\n");
    for (a = 0; a < dof_ids[0]; a++) {
	char *dof_name = parse_name(g, a);
	dof = phgDofIndex2Pointer(g, dof_ids[a + 1]);
	dim = DofDim(dof);
	if (dim == 1) 
	    phgIOPrintf(g->comm, "scalar ");
	else if (dim == 3)
	    phgIOPrintf(g->comm, "vector ");
	else if (dim == 9)
	    phgIOPrintf(g->comm, "tensor ");
	else 
	    phgError(0, "Not allow dim = %d in ensight.", dim);
	    
	if (dof->type == DOF_P0 || dof->type == DOF_DG0)
	    phgIOPrintf(g->comm, "per element:   ");
	else 
	    phgIOPrintf(g->comm, "per node:      ");
	phgIOPrintf(g->comm, "%s %s.*****.%s\n", dof_name, filename, dof_name);
    }
    phgIOPrintf(g->comm, "scalar per element:   %s %s.*****.%s\n",
		"region_mark", filename, "region");
    phgIOPrintf(g->comm, "scalar per element:   %s %s.*****.%s\n",
		"submesh_indice", filename, "rank");

    /* Case file: Time section */
    phgIOPrintf(g->comm, "\nTIME\n");
    phgIOPrintf(g->comm, "time set:              %1d\n", 1);
    phgIOPrintf(g->comm, "number of steps:       %5d\n", 1);
    phgIOPrintf(g->comm, "filename start number: %05d\n", tstep);
    phgIOPrintf(g->comm, "filename increment:    %1d\n", 1);
    phgIOPrintf(g->comm, "time values:\n");
    phgIOPrintf(g->comm, "%14.8e\n", time);

    phgIOPrintf(g->comm, "\n# End of case file");
    phgIOClose();

    phgFree(buffer);
    buffer = NULL;

    return;
}

const char *
phgExportEnsightTn(GRID *g, const char *fn, FLOAT new_time, INT new_tstep, int ndof, DOF **dofs)
/* Outputs mesh in Ensight format, returns the output filename.
 * The variable part of the arguments is a NULL terminated list of DOFs
 * to save to the Ensight file */
{
    int i, j, k;
    DOF *dof;
    BTYPE *types_vert = NULL;

    tstep = new_tstep;
    time = new_time;
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
			   "phgExportEnsight(). \n");
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
		phgWarning("phgExportEnsight: missing NULL terminator"
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

    phgExportEnsight_(g);

    if (g->period != NULL) {
	phgFree(g->types_vert);
	g->types_vert = types_vert;
    }

    return filename;
}


const char *
phgExportEnsightT(GRID *g, const char *fn, FLOAT new_time, INT new_tstep, DOF *dof, ...)
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

    tstep = new_tstep;
    time = new_time;

    fn = phgExportEnsightTn(g, fn, new_tstep, new_time, ndof, dofs);
    phgFree(dofs);

    return fn;
}

const char *
phgExportEnsight(GRID *g, const char *fn, DOF *dof, ...)
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

    fn = phgExportEnsightTn(g, fn, 0, 0, ndof, dofs);
    phgFree(dofs);

    return fn;
}

const char *
phgExportEnsightn(GRID *g, const char *fn, int ndof, DOF **dofs) {
    return phgExportEnsightTn(g, fn, 0, 0, ndof, dofs);
}

