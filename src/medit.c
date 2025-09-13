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

/* Medit (http://www-rocq.inria.fr/gamma/medit/medit.html)
 * "mesh" file format import/export.
 *
 * $Id: medit.c,v 1.86 2022/06/21 06:44:17 zlb Exp $
*/

#include "phg.h"
#include "phg/elem-info.h"
#include "phg/io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

static int get_token_ret = 0;
static BOOLEAN
get_token(FILE *fp, char *token)
{
    int c;
    char *p;

    while (TRUE) {
	if (fscanf(fp, "%s", token) != 1)
	    return get_token_ret = FALSE;
	if (token[0] != '#')
	    break;
	/* skip to newline */
	do
	    if ((c = fgetc(fp)) == EOF)
		return get_token_ret = FALSE;
	while (c != '\n');
    }
    if ((p = strchr(token, '#')) != NULL)
	*p = '\0';

    phgInfo(4, "  got token \"%s\"\n", token);

    return get_token_ret = TRUE;
}

static GRID *g_;

static int
comp_face3(const void *p0, const void *p1)
/* compares two triangular faces */
{
    INT (*f0)[4] = (void *)p0, (*f1)[4] = (void *)p1;
    int i;

    if ((i = (*f0)[0] - (*f1)[0]) != 0)
	return i;
    if ((i = (*f0)[1] - (*f1)[1]) != 0)
	return i;
    return (*f0)[2] - (*f1)[2];
}

static void
process_tetrahedra(GRID *g, INT nlist, INT (*list)[4], INT *vtypes)
/* assign boundary types (e->bound_type[]) using the list of triangles */
{
    /* Note: each face is stored as (v0, v1, v2, elem_no, face_no) */
    INT i;
    INT key[4], (*found)[4];
    ELEMENT *e;
    int j;

    if (nlist > 0)
	qsort(list, nlist, sizeof(*list), comp_face3);

    for (i = 0; i < g->nroot; i++) {
	e = g->roots + i;
	/* element type was saved in mark */
	for (j = 0; j < NFace; j++) {
	    key[0] = e->verts[GetFaceVertex(j, 0)];
	    key[1] = e->verts[GetFaceVertex(j, 1)];
	    key[2] = e->verts[GetFaceVertex(j, 2)];
	    e->bound_type[j] = UNDEFINED;
	    /* adjust bdry type according to types of triangles */
	    if (nlist == 0)
		continue;
	    qsort(key, 3, sizeof(INT), phgCompINT);
	    found = bsearch(key, list, nlist, sizeof(*list), comp_face3);
	    if (found != NULL) {
		e->bound_type[j] = _phg_bcmap(g_, (*found)[3]);
	    }
	}
    }
}

static int
comp_face4(const void *p0, const void *p1)
/* compares two triangular faces */
{
    INT (*f0)[5] = (void *)p0, (*f1)[5] = (void *)p1;
    int i;

    if ((i = (*f0)[0] - (*f1)[0]) != 0)
	return i;
    if ((i = (*f0)[1] - (*f1)[1]) != 0)
	return i;
    if ((i = (*f0)[2] - (*f1)[2]) != 0)
	return i;
    return (*f0)[3] - (*f1)[3];
}

static int region_mark;
static INT ntetra, ntetra_allocated;
static ELEMENT *tetras;

static void
new_tetra(int v0, int v1, int v2, int v3, int bound_type[])
{
    ELEMENT *e;
    int i, b;

    if (v0 == v1 || v0 == v2 || v0 == v3 || v1 == v2 || v1 == v3 || v2 == v3) {
	phgWarning("deleting invalid tetrahedron: %d %d %d %d\n", v0,v1,v2,v3);
	return;
    }
    if (ntetra >= ntetra_allocated) {
	ntetra_allocated += 1024;
	tetras = phgReallocElements(tetras, ntetra, ntetra_allocated);
    }

    e = tetras + (ntetra++);
    e->verts[0] = v0;
    e->verts[1] = v1;
    e->verts[2] = v2;
    e->verts[3] = v3;
    e->region_mark = region_mark;
    if (bound_type != NULL) {
	for (i = 0; i < NVert; i++) {
	    b = bound_type[i];
	    if (b >= 0 && b != UNDEFINED)
		e->bound_type[i] = b;
	}
    }
}

static void
process_hexahedra(GRID *g, int nhexa, int (*hexas)[9],
		  INT nlist, INT (*list)[5], INT *vtypes)
/* assign boundary types (e->bound_type[]) using the list of triangles */
{
    /* Note: each face is stored as (v0, v1, v2, elem_no, face_no) */
    INT i;
    INT key[5], (*found)[5];
    int j, bound_type[6];

    if (nlist > 0)
	qsort(list, nlist, sizeof(*list), comp_face4);

    ntetra = ntetra_allocated = g->nroot;
    tetras = g->roots;
    for (i = 0; i < nhexa; i++) {
	for (j = 0; j < 6; j++) {
	    key[0] = hexas[i][(int)phgElemInfo[ET_HEXA].face_info[j][1]];
	    key[1] = hexas[i][(int)phgElemInfo[ET_HEXA].face_info[j][2]];
	    key[2] = hexas[i][(int)phgElemInfo[ET_HEXA].face_info[j][3]];
	    key[3] = hexas[i][(int)phgElemInfo[ET_HEXA].face_info[j][4]];
	    bound_type[j] = UNDEFINED;
	    /* adjust bdry type according to types of quadrilaterals */
	    if (nlist == 0)
		continue;
	    qsort(key, 4, sizeof(INT), phgCompINT);
	    found = bsearch(key, list, nlist, sizeof(*list), comp_face4);
	    if (found != NULL)
		bound_type[j] = _phg_bcmap(g_, (*found)[4]);
	}
	/* convert to tetrahedra */
	region_mark = hexas[i][8];
	phgHexa2Tetra(hexas[i], new_tetra, bound_type);
    }
    phgInfo(2, "%d hexahedr%s converted to %d tetrahedra\n",
		nhexa, nhexa > 1 ? "a" : "on", ntetra - g->nroot);
    g->ntree = g->nleaf = g->nelem = g->nroot = g->nelem_global = ntetra;
    g->roots = phgRealloc_(tetras, ntetra * sizeof(*g->roots),
					ntetra * sizeof(*g->roots));
}

int
phgImportMedit(GRID *g, const char *filename, BOOLEAN parallel)
{
    FILE *fp;
    char token[1024]; 
    int i, n;
    INT *vtypes = NULL;			/* vertex types */
    INT nlist3 = 0, (*list3)[4] = NULL;	/* list of triangles */
    INT nlist4 = 0, (*list4)[5] = NULL;	/* list of quadrilaterals */
    INT nlistt = 0;			/* number of terahedra */
    int nlisth = 0, (*listh)[9] = NULL;	/* list of hexahedra */

#define READ_NUMBER						\
    if (!get_token(fp, token)) strcpy(token, "End");		\
    if (isalpha((int)(token[0]))) {				\
	phgWarning("fewer entries (%d) than expected.\n", i);	\
	break;							\
    }

    FunctionEntry;

    (void)(parallel);	/* avoid gcc warning */

    if (phgRank >= g->nprocs)
	Return 0;

    if ((fp = phgOpenInputFile_(filename)) == NULL) {
	phgWarning("can't open mesh file \"%s\"!\n", filename);
	phgFreeGrid(&g);
	Return -1;
    }

    g_ = g;

    phgInfo(1, "loading file \"%s\".\n", filename);

    token[0] = '\0';
    if (!get_token(fp, token) || strcasecmp(token, "MeshVersionFormatted") ||
	!get_token(fp, token) || (strcmp(token, "1") && strcmp(token, "2")) ||
	!get_token(fp, token) || strcasecmp(token, "Dimension") ||
	!get_token(fp, token) || strcmp(token, "3")) {
	n = __LINE__;
error:
	phgError(1, "(%s:%d) invalid input file: ret=%s, token=%s\n",
			__FILE__, n, get_token_ret ? "TRUE" : "FALSE", token);
    }

#undef ERROR
#define ERROR	{n = __LINE__; goto error;}

    while (TRUE) {
	if (!get_token(fp, token))
	    break;
next_token:
	if (!strcasecmp(token, "End")) {
	    break;
	}
	else if (!strcasecmp(token, "Vertices")) {
	    BOOLEAN flag = FALSE;
	    if (!get_token(fp, token))
		ERROR
	    n = atoi(token);
	    phgInfo(2, "number of vertices: %d\n", n);
	    g->verts = phgNewVertices(n);
	    vtypes = phgAlloc(n * sizeof(*vtypes)); 
	    for (i = 0; i < n; i++) {
		READ_NUMBER
		g->verts[i][0] = atof(token);
		if (!get_token(fp, token))
		    ERROR
		g->verts[i][1] = atof(token);
		if (!get_token(fp, token))
		    ERROR
		g->verts[i][2] = atof(token);
		/* the type of the vertex */
		if (!get_token(fp, token))
		    ERROR
		if ((vtypes[i] = atoi(token)) != 0)
		    flag = TRUE;
	    }
	    g->nvert_global = g->nvert = i;
	    if (!flag) {
		phgFree(vtypes);
		vtypes = NULL;
	    }
	    if (i < n)
		goto next_token;
	}
	else if (!strcasecmp(token, "Edges")) {
	    if (!get_token(fp, token))
		ERROR
	    n = atoi(token);
	    phgInfo(2, "number of edges: %d (ignored)\n", n);
	    for (i = 0; i < n; i++) {
		READ_NUMBER
		if (!get_token(fp, token))
		    ERROR
		if (!get_token(fp, token))
		    ERROR
	    }
	    if (i < n)
		goto next_token;
	}
	else if (!strcasecmp(token, "Triangles")) {
	    if (!get_token(fp, token))
		ERROR
	    /* extract boundary types from list of triangles */
	    n = atoi(token);
	    phgInfo(2, "number of triangles: %d\n", n);
	    list3 = phgRealloc_(list3, (n + nlist3) * sizeof(*list3),
					nlist3 * sizeof(*list3));
	    for (i = 0; i < n; i++) {
		READ_NUMBER
		list3[nlist3][0] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		list3[nlist3][1] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		list3[nlist3][2] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		list3[nlist3][3] = atoi(token);	/* boundary type */
		/* sort the three vertices */
		qsort(list3[nlist3++], 3, sizeof(INT), phgCompINT);
	    }
	    if (nlist3 == 0) {
		phgFree(list3);
		list3 = NULL;
	    }
	    if (i < n)
		goto next_token;
	}
	else if (!strcasecmp(token, "Quadrilaterals")) {
	    if (!get_token(fp, token))
		ERROR
	    n = atoi(token);
	    phgInfo(2, "number of quadrilaterals: %d\n", n);
	    list4 = phgRealloc_(list4, (n + nlist4) * sizeof(*list4),
					nlist4 * sizeof(*list4));
	    for (i = 0; i < n; i++) {
		READ_NUMBER
		list4[nlist4][0] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		list4[nlist4][1] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		list4[nlist4][2] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		list4[nlist4][3] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		list4[nlist4][4] = atoi(token);	/* boundary type */
		/* sort the four vertices */
		qsort(list4[nlist4++], 4, sizeof(INT), phgCompINT);
	    }
	    if (nlist4 == 0) {
		phgFree(list4);
		list4 = NULL;
	    }
	    if (i < n)
		goto next_token;
	}
	else if (!strcasecmp(token, "Tetrahedra")) {
	    if (!get_token(fp, token))
		    ERROR
	    n = atoi(token);
	    phgInfo(2, "number of tetrahedra: %d\n", n);
	    g->roots = phgReallocElements(g->roots, nlistt, n + nlistt);
	    for (i = 0; i < n; i++) {
		READ_NUMBER
		g->roots[nlistt].verts[0] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		g->roots[nlistt].verts[1] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		g->roots[nlistt].verts[2] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		g->roots[nlistt].verts[3] = atoi(token) - 1;
		/* save type of tetrahedron in 'mark' (currently not used) */
		if (!get_token(fp, token))
		    ERROR
		g->roots[nlistt].region_mark =
			g->roots[nlistt].mark = atoi(token);
		nlistt++;
	    }
	    if (i < n)
		goto next_token;
	}
	else if (!strcasecmp(token, "Hexahedra")) {
	    if (!get_token(fp, token))
		ERROR
	    n = atoi(token);
	    phgInfo(2, "number of hexahedra: %d\n", n);
	    listh = phgRealloc_(listh, (n + nlisth) * sizeof(*listh),
					nlisth * sizeof(*listh));
	    for (i = 0; i < n; i++) {
		READ_NUMBER
		listh[nlisth][0] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		listh[nlisth][1] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		listh[nlisth][2] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		listh[nlisth][3] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		listh[nlisth][4] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		listh[nlisth][5] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		listh[nlisth][6] = atoi(token) - 1;
		if (!get_token(fp, token))
		    ERROR
		listh[nlisth][7] = atoi(token) - 1;
		/* The element type */
		if (!get_token(fp, token))
		    ERROR
		listh[nlisth++][8] = atoi(token);
	    }
	    if (i < n)
		goto next_token;
	}
	else if (!strcasecmp(token, "Corners")) {
	    if (!get_token(fp, token))
		    ERROR
	    n = atoi(token);
	    for (i = 0; i < n; i++) {
		READ_NUMBER
	    }
	    if (i < n)
		goto next_token;
	}
	else if (!strcasecmp(token, "RequiredVertices")) {
	    if (!get_token(fp, token))
		    ERROR
	    n = atoi(token);
	    for (i = 0; i < n; i++) {
		READ_NUMBER
	    }
	    if (i < n)
		goto next_token;
	}
	else if (!strcasecmp(token, "Ridges")) {
	    if (!get_token(fp, token))
		    ERROR
	    n = atoi(token);
	    for (i = 0; i < n; i++) {
		READ_NUMBER
	    }
	    if (i < n)
		goto next_token;
	}
	else if (!strcasecmp(token, "RequiredEdges")) {
	    if (!get_token(fp, token))
		    ERROR
	    n = atoi(token);
	    for (i = 0; i < n; i++) {
		READ_NUMBER
	    }
	    if (i < n)
		goto next_token;
	}
	else {
#if 0
	    int lineno = 1;
	    size_t pos, current_pos = ftell(fp);
	    fseek(fp, 0, SEEK_SET);
	    for (pos = 0; pos < current_pos; pos++)
		if (fgetc(fp) == '\n')
		    lineno++;
	    phgCloseInputFile_(fp);
	    phgError(1, "unkown token \"%s\" near line %d.\n", token, lineno);
#else
	    /* tetgen 1.4.0 seems to generate inconsistent (smaller) number of
	     * triangles, as a workaround extra numbers are silently ignored */
	    static BOOLEAN warned = FALSE;
	    if (!warned) {
		int lineno = 1;
		size_t pos, current_pos = ftell(fp);
		warned = TRUE;
		fseek(fp, 0, SEEK_SET);
		for (pos = 0; pos < current_pos; pos++)
		    if (fgetc(fp) == '\n')
			lineno++;
		phgWarning("extra token \"%s\" near line %d ignored.\n",
				token, lineno);
	    }
	    continue;
#endif
	}
    }

    phgCloseInputFile_(fp);

    g->nelem_global = g->ntree = g->nleaf = g->nelem = g->nroot = nlistt;

    if (g->nroot > 0)
	process_tetrahedra(g, nlist3, list3, vtypes);

    if (list3 != NULL) {
	phgFree(list3);
	list3 = NULL;
	nlist3 = 0;
    }

    if (listh != NULL) {
	process_hexahedra(g, nlisth, listh, nlist4, list4, vtypes);
	phgFree(listh);
	listh = NULL;
	nlisth = 0;
    }

    if (list4 != NULL) {
	phgFree(list4);
	list4 = NULL;
	nlist4 = 0;
    }

    if (vtypes != NULL) {
	phgFree(vtypes);
	vtypes = NULL;
    }

    Return 0;
    return 0;	/* for avoiding MIPS cc warning */
}

BOOLEAN
phgExportMedit(GRID *g, const char *filename)
/* Output Medit mesh file. filename == NULL => output to stdout. */
{
    INT i, j, n;
    FILE *f = NULL;
    ELEMENT *s;
    VEF_MAP *vef;
    /* dof and map are used to reorder vertices */
    DOF *dof;
    MAP *map;
    INT *v2vmap;
    int bc, *found;

    /* variables for backing up members of 'g' if g->period != NULL */
    BYTE	flags;
    PERIOD	*period;
    BTYPE	*types_vert;
    INT		nvert_owned, *owner_index_vert;
    int		*owner_rank_vert;

#if USE_MPI
# define TRUNK_SIZE 65536
    MPI_Status status;
    FLOAT *fbuffer;
    INT *ibuffer;
#endif	/* USE_MPI */

    FunctionEntry;

    if (g == NULL || phgRank >= g->nprocs)
	Return TRUE;

    assert((g->flags & FACE_FLAG));

    if (g->period != NULL) {
	/* FIXME: also make a backup of g->elems? */
	period = g->period;
	g->period = NULL;
	types_vert = g->types_vert;
	g->types_vert = NULL;
	owner_index_vert = g->owner_index_vert;
	g->owner_index_vert = NULL;
	owner_rank_vert = g->owner_rank_vert;
	g->owner_rank_vert = NULL;
	flags = g->flags;
	g->flags = VERT_FLAG;	/* phgUpdateBoundaryTypes will only change
				 * xxxx_vert[] */
	nvert_owned = g->nvert_owned;
	phgUpdateBoundaryTypes(g);
    }
    else {
	period = NULL;
	/* the following lines make gcc happy */
	types_vert = NULL;
	owner_index_vert = NULL;
	owner_rank_vert = NULL;
	flags = 0;
	nvert_owned = 0;
    }

    dof = phgDofNew(g, DOF_P1, 1, "dof", DofNoData);
    dof->DB_mask = 0;
    map = phgMapCreate(dof, NULL);

    /*----------------------- File header and vertices --------------------*/

#define VMap(i)	((phgMapD2V(map, 0, i) < map->nlocal ? \
			phgMapD2V(map, 0, i) + map->partition[map->rank] : \
			map->O2Gmap[phgMapD2V(map, 0, i) - map->nlocal]) + 1) 

    /* construct VEC ==> VERTEX map */
    v2vmap = phgAlloc(map->nlocal * sizeof(*v2vmap));
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert != NULL && !(g->types_vert[i] & OWNER))
	    continue;
	assert(phgMapD2V(map, 0, i) <= map->nlocal);
	v2vmap[phgMapD2V(map, 0, i)] = i;
    }

    if (g->rank == 0) {
	if (filename == NULL) {
	    f = stdout;
	}
	else {
	    if ((f = fopen(filename, "w+t")) == NULL) {
		phgError(0, "cannot open output file \"%s\".\n", filename);
		Return FALSE;
	    }
	    phgInfo(1, "creating Medit mesh file \"%s\".\n", filename);
	}
	/* File header */
	fprintf(f, "MeshVersionFormatted 1\n");
	fprintf(f, "Dimension %d\n", Dim);
	/* Vertices */
	fprintf(f, "\nVertices %"dFMT"\n", g->nvert_global);
	for (i = 0; i < map->nlocal; i++) {
	    j = v2vmap[i];
	    fprintf(f, " %0.16lg %0.16lg %0.16lg 0\n", (double)g->verts[j][0],
			(double)g->verts[j][1], (double)g->verts[j][2]);
	}
#if USE_MPI
	/* receive vertices from other processes */
	fbuffer = phgAlloc(TRUNK_SIZE * Dim * sizeof(*fbuffer));
	for (j = 1; j < g->nprocs; j++) {
	    while (TRUE) {
		MPI_Probe(j, 111, g->comm, &status);
		MPI_Get_count(&status, PHG_MPI_FLOAT, &bc);
		n = bc;
		assert(n <= (Dim * TRUNK_SIZE) && n % Dim == 0);
		MPI_Recv(fbuffer, n, PHG_MPI_FLOAT, j, 111, g->comm, &status);
		/* process received vertices */
		for (i = 0; i < n; i += Dim) {
		    fprintf(f, " %0.16lg %0.16lg %0.16lg 0\n",
				(double)fbuffer[i + 0],
				(double)fbuffer[i + 1],
				(double)fbuffer[i + 2]);
		}
		if (n < Dim * TRUNK_SIZE)
		    break;
	    }
	}
	phgFree(fbuffer);
#endif	/* USE_MPI */
    }
#if USE_MPI
    else {
	/* send vertices to root process */
	fbuffer = phgAlloc(TRUNK_SIZE * Dim * sizeof(*fbuffer));
	for (i = n = 0; i < map->nlocal; i++) {
	    j = v2vmap[i];
	    memcpy(fbuffer + Dim * n, g->verts[j], Dim * sizeof(*fbuffer));
	    if (++n >= TRUNK_SIZE) {
		MPI_Send(fbuffer, Dim * n, PHG_MPI_FLOAT, 0, 111, g->comm);
		n = 0;
	    }
	}
	/* send the last block (may be size 0), which also marks EOD */
	MPI_Send(fbuffer, Dim * n, PHG_MPI_FLOAT, 0, 111, g->comm);
	phgFree(fbuffer);
    }
#endif	/* USE_MPI */
    phgFree(v2vmap);

    /*------------------------------ Tetrahedra ---------------------------*/
    if (g->rank == 0) {
	fprintf(f, "\nTetrahedra %"dFMT"\n", g->nelem_global);
	ForAllElements(g, s)
	    fprintf(f, " %"dFMT" %"dFMT" %"dFMT" %"dFMT" %d\n",
			VMap(s->verts[0]), VMap(s->verts[1]),
			VMap(s->verts[2]), VMap(s->verts[3]),
			s->region_mark);
#if USE_MPI
	/* receive elements from other processes */
	ibuffer = phgAlloc(TRUNK_SIZE * (Dim + 2) * sizeof(*ibuffer));
	for (j = 1; j < g->nprocs; j++) {
	    while (TRUE) {
		MPI_Probe(j, 222, g->comm, &status);
		MPI_Get_count(&status, PHG_MPI_INT, &bc);
		n = bc;
		assert(n <= ((Dim + 2) * TRUNK_SIZE) && n % (Dim + 2) == 0);
		MPI_Recv(ibuffer, n, PHG_MPI_INT, j, 222, g->comm, &status);
		/* process received vertices */
		for (i = 0; i < n; i += Dim + 2) {
		    fprintf(f, " %"dFMT" %"dFMT" %"dFMT" %"dFMT" %"dFMT"\n",
				ibuffer[i + 0], ibuffer[i + 1], ibuffer[i + 2],
				ibuffer[i + 3], ibuffer[i + 4]);
		}
		if (n < (Dim + 2) * TRUNK_SIZE)
		    break;
	    }
	}
	phgFree(ibuffer);
#endif	/* USE_MPI */
    }
#if USE_MPI
    else {
	/* send elements to root process */
	ibuffer = phgAlloc(TRUNK_SIZE * (Dim + 2) * sizeof(*ibuffer));
	n = 0;
	ForAllElements(g, s) {
	    ibuffer[(Dim + 2) * n + 0] = VMap(s->verts[0]);
	    ibuffer[(Dim + 2) * n + 1] = VMap(s->verts[1]);
	    ibuffer[(Dim + 2) * n + 2] = VMap(s->verts[2]);
	    ibuffer[(Dim + 2) * n + 3] = VMap(s->verts[3]);
	    ibuffer[(Dim + 2) * n + 4] = s->region_mark;
	    if (++n >= TRUNK_SIZE) {
		MPI_Send(ibuffer, (Dim + 2) * n, PHG_MPI_INT, 0, 222, g->comm);
		n = 0;
	    }
	}
	/* send the last block (may be size 0), which also marks EOD */
	MPI_Send(ibuffer, (Dim + 2) * n, PHG_MPI_INT, 0, 222, g->comm);
	phgFree(ibuffer);
    }
#endif	/* USE_MPI */

    /*-------------------------------- Faces ------------------------------*/
    vef = phgDofSetupVEFMap(g, NULL, FACE_FLAG);

    /* count # of boundary faces.
     *
     * Note: when g->period != NULL, each set of periodic faces has only one
     *	     owner, this is not a problem because they are all non boundary
     *	     faces */
    for (i = n = 0; i < g->nface; i++) {
	if ((g->types_face[i] & OWNER) == 0 ||
	    (g->types_face[i] & BDRY_MASK) == UNDEFINED)
	    continue;
	assert(vef->Fmap[i] != NULL);
	s = vef->Fmap[i];
	bc = (s->bound_type[vef->Find[i]] & BDRY_MASK) >> BDRY_SHIFT;
	if (bc == 0)
	    continue;
	if (g->bc_list != NULL) {
	    found = bsearch(&bc, g->bc_list, g->bc_n, sizeof(*g->bc_list),
			phgCompint);
	    if (found == NULL)
		continue;
	}
	n++;
    }

#if USE_MPI
    j = n;
    MPI_Allreduce(&j, &n, 1, PHG_MPI_INT, MPI_SUM, g->comm);
#endif	/* USE_MPI */

    if (n == 0) {
	if (g->rank == 0) {
	    fprintf(f, "\nEnd\n");
	    ((f != stdout && f != stderr) ? fclose : fflush)(f);
	}
    }
    else if (g->rank == 0) {
	fprintf(f, "\nTriangles %"dFMT"\n", n);
	for (i = 0; i < g->nface; i++) {
	    if ((g->types_face[i] & OWNER) == 0 ||
		(g->types_face[i] & BDRY_MASK) == UNDEFINED)
		continue;
	    s = vef->Fmap[i];
	    bc = (s->bound_type[vef->Find[i]] & BDRY_MASK) >> BDRY_SHIFT;
	    if (bc == 0)
		continue;
	    if (g->bc_list != NULL) {
		found = bsearch(&bc, g->bc_list, g->bc_n, sizeof(*g->bc_list),
				phgCompint);
		if (found == NULL)
		    continue;
		bc = g->bc_rmap[found - g->bc_list];
	    }
	    fprintf(f, " %"dFMT" %"dFMT" %"dFMT" %d\n",
			VMap(s->verts[GetFaceVertex(vef->Find[i], 0)]),
			VMap(s->verts[GetFaceVertex(vef->Find[i], 1)]),
			VMap(s->verts[GetFaceVertex(vef->Find[i], 2)]),
			bc);
	}

#if USE_MPI
	/* receive faces from other processes */
	ibuffer = phgAlloc(TRUNK_SIZE * (Dim + 1) * sizeof(*ibuffer));
	for (j = 1; j < g->nprocs; j++) {
	    while (TRUE) {
		MPI_Probe(j, 333, g->comm, &status);
		MPI_Get_count(&status, PHG_MPI_INT, &bc);
		n = bc;
		assert(n <= ((Dim + 1) * TRUNK_SIZE) && n % (Dim + 1) == 0);
		MPI_Recv(ibuffer, n, PHG_MPI_INT, j, 333, g->comm, &status);
		/* process received vertices */
		for (i = 0; i < n; i += Dim + 1) {
		    fprintf(f, " %"dFMT" %"dFMT" %"dFMT" %"dFMT"\n",
				ibuffer[i + 0], ibuffer[i + 1], ibuffer[i + 2],
				ibuffer[i + 3]);
		}
		if (n < (Dim + 1) * TRUNK_SIZE)
		    break;
	    }
	}
	phgFree(ibuffer);
#endif	/* USE_MPI */

	fprintf(f, "\nEnd\n");
	((f != stdout && f != stderr) ? fclose : fflush)(f);
    }
#if USE_MPI
    else {
	/* send boundary faces to root process */
	ibuffer = phgAlloc(TRUNK_SIZE * (Dim + 1) * sizeof(*ibuffer));
	for (i = n = 0; i < g->nface; i++) {
	    if (!(g->types_face[i] & OWNER))
		continue;
	    s = vef->Fmap[i];
	    bc = (s->bound_type[vef->Find[i]] & BDRY_MASK) >> BDRY_SHIFT;
	    if (bc == 0)
		continue;
	    if (g->bc_list != NULL) {
		found = bsearch(&bc, g->bc_list, g->bc_n, sizeof(*g->bc_list),
				phgCompint);
		if (found == NULL)
		    continue;
		bc = g->bc_rmap[found - g->bc_list];
	    }
	    ibuffer[(Dim + 1) * n + 0] =
			VMap(s->verts[GetFaceVertex(vef->Find[i], 0)]);
	    ibuffer[(Dim + 1) * n + 1] =
			VMap(s->verts[GetFaceVertex(vef->Find[i], 1)]);
	    ibuffer[(Dim + 1) * n + 2] =
			VMap(s->verts[GetFaceVertex(vef->Find[i], 2)]);
	    ibuffer[(Dim + 1) * n + 3] = bc;
	    if (++n >= TRUNK_SIZE) {
		MPI_Send(ibuffer, (Dim + 1) * n, PHG_MPI_INT, 0, 333, g->comm);
		n = 0;
	    }
	}
	/* send the last block (may be size 0), which also marks EOD */
	MPI_Send(ibuffer, (Dim + 1) * n, PHG_MPI_INT, 0, 333, g->comm);
	phgFree(ibuffer);
    }
#endif	/* USE_MPI */

    phgDofFreeVEFMap(&vef);

    /*------------------------------- Cleanup -----------------------------*/

    phgMapDestroy(&map);
    phgDofFree(&dof);

    if (period != NULL) {
	g->period = period;
	phgFree(g->types_vert);
	g->types_vert = types_vert;
	phgFree(g->owner_index_vert);
	g->owner_index_vert = owner_index_vert;
	phgFree(g->owner_rank_vert);
	g->owner_rank_vert = owner_rank_vert;
	g->flags = flags;
	g->nvert_owned = nvert_owned;
    }

    Return TRUE;
    return TRUE;	/* for avoiding MIPS cc warning */
}
