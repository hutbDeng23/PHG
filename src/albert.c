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

/* $Id: albert.c,v 1.33 2021/12/08 02:48:52 zlb Exp $
 *
 * ALBERTA macro triangulation importer/exporter */

#include "phg.h"
#include "phg/io.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* macro for converting geometric type from/to ALBERT */
#define GEOM_TYPE_FROM_ALBERT(n) \
	((GTYPE)((n) == 0 ? DIAGONAL : ((n) == 1 ? FACE : ((n) == 2 ? EDGE : \
	((n) == 3 ? MIXED : OPPOSITE)))))
#define GEOM_TYPE_TO_ALBERT(n) \
	((n) == DIAGONAL ? 0 : ((n) == FACE ? 1 : ((n) == EDGE ? 2 : \
	((n) == MIXED ? 3 : 4))))

/* macro for converting boundary type from/to ALBERT */
#define BDRY_TYPE_FROM_ALBERT(n) \
	((BTYPE)((n) == 0 ? INTERIOR: \
		((n) < 0  ? NEUMANN : \
		((n) == 1 ? DIRICHLET : \
		((n) == 2 ? BDRY_USER1 : \
		((n) == 3 ? BDRY_USER2 : \
		((n) == 4 ? BDRY_USER3 : \
		((n) == 5 ? BDRY_USER4 : \
		((n) == 6 ? BDRY_USER5 : \
		((n) == 7 ? BDRY_USER6 : \
		((n) == 8 ? BDRY_USER7 : \
		((n) == 9 ? BDRY_USER8 : \
		((n) == 10 ? BDRY_USER9 : \
		((n) == 11 ? BDRY_USER0 : UNDEFINED))))))))))))))
#define BDRY_TYPE_TO_ALBERT(n) \
	(((n) & REMOTE)	    ? 4 : \
	(((n) & DIRICHLET)  ? 1 : \
	(((n) & NEUMANN)    ? -1 : \
	(((n) & BDRY_USER1) ? 2 : \
	(((n) & BDRY_USER2) ? 3 : \
	(((n) & BDRY_USER3) ? 4 : \
	(((n) & BDRY_USER4) ? 5 : \
	(((n) & BDRY_USER5) ? 6 : \
	(((n) & BDRY_USER6) ? 7 : \
	(((n) & BDRY_USER7) ? 8 : \
	(((n) & BDRY_USER8) ? 9 : \
	(((n) & BDRY_USER9) ? 10 : \
	(((n) & BDRY_USER0) ? 11 : \
	(((n) & UNDEFINED)  ? 4 : 0))))))))))))))
	
static GRID grid;
static FILE *fps[32], **fp;	/* max nested #include depth = 32 */
static char line[1024];
static BOOLEAN has_element_type;

static void read_elements(void);
static void read_boundary_types(void);
static void read_coordinates(void);
static void read_neighbours(void);
static void read_bdry_funcs(void);
#if (Dim > 2)
static void read_types(void);
#endif

static struct {
    char	*key;		/* key word */
    char	*fmt;		/* scanf format string */
    union {
        void *loc;		/* addr of the var if fmt!=null */
        void (*func)(void);	/* function if fmt==null */
    } pointer;
    CHAR	mandatory;	/* >0: mandatory, ==0: warn if missing, 
				   <0: don't warn */
    BOOLEAN	flag;		/* FALSE ==> uninitialized */
} *pp, *pp1, phgmp, albertParms[] = {
    {"DIM",			"%" dFMT, {NULL}, -1, FALSE},
    {"DIM_OF_WORLD",		"%" dFMT, {NULL}, -1, FALSE},
    {"number of vertices",	"%" dFMT, {&grid.nvert}, 1, FALSE},
    {"number of elements",	"%" dFMT, {&grid.nroot}, 1, FALSE},
    {"vertex coordinates",	NULL, {read_coordinates}, 1, FALSE},
    {"element vertices",	NULL, {read_elements}, 1, FALSE},
    {"element boundaries",	NULL, {read_boundary_types}, 0, FALSE},
    {"element neighbours",	NULL, {read_neighbours}, -1, FALSE},
    {"curved boundaries",	NULL, {read_bdry_funcs}, -1, FALSE},
#if (Dim > 2)
    {"element type",		NULL, {read_types}, -1, FALSE},
#endif
};
#define albertParms_n	(sizeof(albertParms)/sizeof(albertParms[0]))

static void
read_elements(void)
/* Read vertex indices of all elements */
{
#if (Dim == 2)
    INT i, n1, n2, n3;
#else
    INT i, n1, n2, n3, n4;
#endif

    if (grid.nroot <= 0) {
	phgError(1, "%s:%d, \"number of elements\" should come before "
		   "\"element vertices\".\n", __FILE__, __LINE__);
    }

    if (grid.roots == NULL) grid.roots = phgNewElements(grid.nroot);

    for (i = 0; i < grid.nroot; i++) {
	INT *v = grid.roots[i].verts;
	if (fgets(line, sizeof(line), *fp) == NULL)
	    phgError(1, "%s:%d, unexpected eof.\n", __FILE__, __LINE__);
#if (Dim == 2)
	if (sscanf(line, "%"dFMT" %"dFMT" %"dFMT, &n1, &n2, &n3) != 3)
#else
	if (sscanf(line, "%"dFMT" %"dFMT" %"dFMT" %"dFMT, &n1, &n2, &n3, &n4) != 4)
#endif
	    phgError(1, "%s:%d, invalid input line: %s",
			line, __FILE__, __LINE__);
	v[0] = n1;
	v[1] = n2;
	v[2] = n3;
#if Dim == 2
	if (n1 < 0 || n1 >= grid.nvert || n2 < 0 || n2 >= grid.nvert ||
	    n3 < 0 || n3 >= grid.nvert)
#else
	v[3] = n4;
	if (n1 < 0 || n1 >= grid.nvert || n2 < 0 || n2 >= grid.nvert ||
	    n3 < 0 || n3 >= grid.nvert || n4 < 0 || n4 >= grid.nvert)
#endif
	    phgError(1, "%s:%d, invalid vertex number on input line: %s",
			line, __FILE__, __LINE__);
    }
}

static void
read_boundary_types(void)
/* Read boundary types of all elements */
{
#if (Dim == 2)
    int i, n1, n2, n3;
#else
    int i, n1, n2, n3, n4;
#endif

    if (grid.nroot <= 0) {
	phgError(1, "%s:%d, \"number of elements\" should come before "
		   "\"element boundaries\".\n", __FILE__, __LINE__);
    }

    if (grid.roots == NULL) grid.roots = phgNewElements(grid.nroot);

    for (i = 0; i < grid.nroot; i++) {
	BTYPE *b = grid.roots[i].bound_type;
	if (fgets(line, sizeof(line), *fp) == NULL)
	    phgError(1, "%s:%d, unexpected eof.\n", __FILE__, __LINE__);
#if (Dim == 2)
	if (sscanf(line, "%d %d %d", &n1, &n2, &n3) != 3)
#else
	if (sscanf(line, "%d %d %d %d", &n1, &n2, &n3, &n4) != 4)
#endif
	    phgError(1, "%s:%d, invalid input line: %s",
			line, __FILE__, __LINE__);
	b[0] = BDRY_TYPE_FROM_ALBERT(n1);
	b[1] = BDRY_TYPE_FROM_ALBERT(n2);
	b[2] = BDRY_TYPE_FROM_ALBERT(n3);
#if (Dim == 3)
	b[3] = BDRY_TYPE_FROM_ALBERT(n4);
#endif
    }
}

#if (Dim > 2)
static void
read_types(void)
/* Read element type */
{
    int i, n;

    if (grid.nroot <= 0) {
	phgError(1, "%s:%d, \"number of elements\" should come before "
		   "\"element type\".\n", __FILE__, __LINE__);
    }

    if (grid.roots == NULL) grid.roots = phgNewElements(grid.nroot);

    for (i = 0; i < grid.nroot; i++) {
	if (fgets(line, sizeof(line), *fp) == NULL)
	    phgError(1, "%s:%d, unexpected eof.\n", __FILE__, __LINE__);
	if (sscanf(line, "%d", &n) != 1 || n < 0 || n > 4)
	    phgError(1, "%s:%d, invalid input line: %s",
			line, __FILE__, __LINE__);
	grid.roots[i].type = GEOM_TYPE_FROM_ALBERT(n);
    }

    has_element_type = TRUE;
}
#endif

static void
read_neighbours(void)
/* Read neighbours */
{
#if (Dim == 2)
    int i, n1, n2, n3;
#else
    int i, n1, n2, n3, n4;
#endif

    phgInfo(1, "neighbours information skipped.\n");

    if (grid.nroot <= 0) {
	phgError(1, "%s:%d, \"number of elements\" should come before "
		   "\"element neighbours\".\n", __FILE__, __LINE__);
    }

    if (grid.roots == NULL) grid.roots = phgNewElements(grid.nroot);

    for (i = 0; i < grid.nroot; i++) {
	if (fgets(line, sizeof(line), *fp) == NULL)
	    phgError(1, "%s:%d, unexpected eof.\n", __FILE__, __LINE__);
#if (Dim == 2)
	if (sscanf(line, "%d %d %d", &n1, &n2, &n3) != 3)
#else
	if (sscanf(line, "%d %d %d %d", &n1, &n2, &n3, &n4) != 4)
#endif
	    phgError(1, "%s:%d, invalid input line: %s",
			line, __FILE__, __LINE__);

#if Dim == 2
	if (n1 >= grid.nroot || n2 >= grid.nroot || n3 >= grid.nroot)
#else
	if (n1 >= grid.nroot || n2 >= grid.nroot ||
	    n3 >= grid.nroot || n4 >= grid.nroot)
#endif
	    phgError(1, "%s:%d, invalid element index, input line: %s",
			line, __FILE__, __LINE__);
    }
}

static void
read_coordinates(void)
/* Read coordinates of all vertices */
{
    int i;
#if (Dim == 2)
    double x, y;
#else
    double x, y, z;
#endif

    if (grid.nvert <= 0) {
	phgError(1, "%s:%d, \"number of vertices\" should come before "
		   "\"vertex coordinates\".\n", __FILE__, __LINE__);
    }

    if (grid.verts == NULL) grid.verts = phgNewVertices(grid.nvert);

    for (i = 0; i < grid.nvert; i++) {
	if (fgets(line, sizeof(line), *fp) == NULL)
	    phgError(1, "%s:%d, unexpected eof.\n", __FILE__, __LINE__);
#if (Dim == 2)
	if (sscanf(line, "%lf %lf", &x, &y) != 2)
#else
	if (sscanf(line, "%lf %lf %lf", &x, &y, &z) != 3)
#endif
	    phgError(1, "%s:%d, invalid input line: %s",
				line, __FILE__, __LINE__);
	grid.verts[i][0] = (FLOAT)x;
	grid.verts[i][1] = (FLOAT)y;
#if (Dim == 3)
	grid.verts[i][2] = (FLOAT)z;
#endif
    }
}

static void
read_bdry_funcs(void)
{
    int i, j, n;
    char *p[4], *q;

    if (fgets(line, sizeof(line), *fp) == NULL) {
error0:
	phgError(1, "%s:%d, unexpected eof.\n", __FILE__, __LINE__);
    }

    if ((n = atoi(line)) <= 0 || n > 127)
	phgError(1, "%s:%d, invalid input line: %s", line, __FILE__, __LINE__);

#if ALLOW_CURVED_BOUNDARY
    grid.bdry_funcs = phgAlloc((4 * n + 1) * sizeof(*grid.bdry_funcs));
#endif
    for (i = 0; i < n; i++) {
	line[0] = '\0';
	q = line;
	for (j = 0, q = line; j < 4; j++) {
	    if (*q == '\0') {
		if (fgets(q, sizeof(line) - (q - line), *fp) == NULL)
		    goto error0;
	    }
	    p[j] = q;
	    if ((q = strchr(q, ';')) != NULL) {
		*q = 0;
		do ++q; while (isspace(*(BYTE *)q));
	    }
	    else {
		q = p[j] + strlen(p[j]);
	    }
	}
#if ALLOW_CURVED_BOUNDARY
	for (j = 0; j < 4; j++) {
	    grid.bdry_funcs[i * 4 + j] = phgDefine3DFunction(p[j]);
	}
#endif
    }
#if ALLOW_CURVED_BOUNDARY
    grid.bdry_funcs[4 * n] = NULL;
#endif
}

static int
comp_keys(const void *p1, const void *p2)
{
    pp = (void *)p1;
    pp1 = (void *)p2;

    return strcmp(pp->key, pp1->key);
}

int
phgImportALBERT(GRID *g, const char *filename, BOOLEAN parallel)
/* Reads grid from an ALBERT parameter file. returns NULL on error */
{
    int i;
    char *p, *q;
    static BOOLEAN initialized = FALSE;

    FunctionEntry;

    Unused(parallel);

    if (phgRank > 0)
	Return 0;

    if (!initialized) {
	qsort(albertParms, albertParms_n, sizeof(albertParms[0]), comp_keys);
	initialized = TRUE;
    }

    memcpy(&grid, g, sizeof(GRID));
    grid.nvert = grid.nroot = -1;
    grid.roots = NULL;
    grid.verts = NULL;
#if 0
    grid.edges = NULL;
#endif

    has_element_type = FALSE;
    
    fp = fps;
    if ((*fp = phgOpenInputFile_(filename)) == NULL) {
	phgWarning("can't open parameter file \"%s\"!\n", filename);
	Return -1;
    }

    phgInfo(1, "loading file \"%s\".\n", filename);

    while (TRUE) {
	if (fgets(line, sizeof(line), *fp) == NULL) {
	    phgCloseInputFile_(*fp);
	    if (fp == fps) break;
	    fp--;
	    continue;
	}
	for (p = line + strlen(line) - 1;
	     p >= line && isspace(*(BYTE *)p);
	     p--);
	*(++p) = '\0';
	for (p = line; isspace(*(BYTE *)p); p++);
	if (*p == '\0') continue;
	if (!strncmp(p, "#include", sizeof("#include") - 1)) {
	    /* process include file */
	    p += sizeof("#include") - 1;
	    while (isspace(*(BYTE *)p)) p++;
	    if ((*p != '\'' && *p != '"') || (q = strchr(p + 1, *p)) == NULL) {
		phgWarning("invalid #include statement.\n");
		continue;
	    }
	    p++; *q = '\0';
	    if (fp - fps >= sizeof(fps)/sizeof(fps[0]) - 1) {
		phgWarning("too many nested #include statements.\n");
		continue;
	    }
	    if ((*(++fp) = phgOpenInputFile_(p)) == NULL) {
		phgWarning("can't open input file \"%s\".\n", p);
		fp--;
	    }
	    phgInfo(1, "processing include file \"%s\".\n",  p);
	    continue;
	}
	/* Parse input line */
	if ((q = strchr(p, ':')) != NULL) *q = '\0';
	phgmp.key = p;
	if (q == NULL || (pp = bsearch(&phgmp, albertParms, albertParms_n,
		sizeof(albertParms[0]), comp_keys)) == NULL) {
	    phgWarning("invalid input line \"%s\", ignored.\n", p);
	    continue;
	}
	phgInfo(3, "key = \"%s\".\n",  pp->key);
	if (pp->flag)
	    phgWarning("duplicate section \"%s\".\n", pp->key);
	else
	    pp->flag = TRUE;
	if (pp->fmt != NULL) {
	    if (pp->pointer.loc != NULL)
		sscanf(q + 1, pp->fmt, pp->pointer.loc);
	} else {
	    void (*phgunc)(void) = pp->pointer.func;
	    phgunc();
	}
    }

    for (i = 0; i < albertParms_n; i++) {
	if (!albertParms[i].flag) {
	    char *key = albertParms[i].key;
	    if (albertParms[i].mandatory > 0) 
		phgError(1, "\"%s\" section missing.\n", key);
	    else if (albertParms[i].mandatory == 0)
		phgWarning("\"%s\" section missing.\n", key);
	} else {
	    albertParms[i].flag = FALSE;	/* reset for next call */
	}
    }

    memcpy(g, &grid, sizeof(GRID));

    g->ntree = g->nleaf = g->nelem = g->nroot;
    g->nelem_global = g->nelem;
    g->nvert_global = g->nvert;

    Return has_element_type;
    return has_element_type;	/* for avoiding MIPS cc warning */
}

static FILE *f;
static INT *indices;
static GRID *g_;

static BOOLEAN
albert_vertices_callback CB_ARGS(s)
{
#if Dim == 2
    fprintf(f, " %"dFMT" %"dFMT" %"dFMT"\n", s->verts[0], s->verts[1],
		s->verts[2]);
#else
    fprintf(f, " %"dFMT" %"dFMT" %"dFMT" %"dFMT"\n", s->verts[0],
		s->verts[1], s->verts[2], s->verts[3]);
#endif
    if (!(g_->flags & ELEM_FLAG))
	s->index = phgTraverseIndex;
    indices[s->index] = phgTraverseIndex;
    return TRUE;
}

static BOOLEAN
albert_boundaries_callback CB_ARGS(s)
{
#if Dim == 2
    fprintf(f, " %d %d %d\n", 
		BDRY_TYPE_TO_ALBERT(s->bound_type[0]), 
		BDRY_TYPE_TO_ALBERT(s->bound_type[1]),
		BDRY_TYPE_TO_ALBERT(s->bound_type[2]));
#else
    fprintf(f, " %d %d %d %d\n", 
		BDRY_TYPE_TO_ALBERT(s->bound_type[0]), 
		BDRY_TYPE_TO_ALBERT(s->bound_type[1]),
		BDRY_TYPE_TO_ALBERT(s->bound_type[2]),
		BDRY_TYPE_TO_ALBERT(s->bound_type[3]));
#endif
    return TRUE;
}

#if Dim > 2
static BOOLEAN
albert_type_callback CB_ARGS(s)
{
    fprintf(f, " %d\n", GEOM_TYPE_TO_ALBERT(s->type));
    return TRUE;
}
#endif

static BOOLEAN
albert_neighbours_callback CB_ARGS(s)
{
    int j;

    for (j = 0; j < NFace; j++) {
	INT n;
	if (!(s->bound_type[j] & INTERIOR) || (s->bound_type[j] & REMOTE)) {
	    /* boundary element */
	    n = -1;
	}
	else if (s->neighbours[j] == NULL) {
	    n = -1;
	}
	else {
	    n = indices[((ELEMENT *)s->neighbours[j])->index];
	}
	fprintf(f, " %"dFMT, n);
    }
    fprintf(f, "\n");

    return TRUE;
}

BOOLEAN
phgExportALBERT(GRID *g, const char *filename)
/* Output ALBERT parameter file. filename == NULL => output to stdout. */
{
    INT i;
    INT nelem = g->nleaf;
    INT nvert = g->nvert;
    char *fn = (char *)filename;

    FunctionEntry;

    if (g == NULL || g->nleaf == 0)
	Return TRUE;

    if (g->period != NULL) {
	phgPrintf("%s: periodic boundary unsupported.\n", __func__);
	Return FALSE;
    }

    if (g->nprocs > 1 && fn != NULL) {
	/* append process number to filename */
	fn = phgAlloc(strlen(filename) + 1 + 5);
	sprintf(fn, "%s-%d", filename, g->rank);
    }

    if (fn == NULL) {
	f = stdout;
    }
    else {
	if ((f = fopen(fn, "w+t")) == NULL) {
	    phgError(0, "cannot open output file \"%s\".\n", fn);
	    if (g->nprocs > 1) phgFree(fn);
	    Return FALSE;
	}
	phgInfo(1, "creating ALBERT file \"%s\".\n", fn);
	if (g->nprocs > 1) phgFree(fn);
    }

    /* The array 'indices' is used to save the global element index in traverse
     * order. It is initialized by albert_vertices_callback */
    indices = phgAlloc(g->nelem * sizeof(*indices));
    g_ = g;

    fprintf(f, "DIM: %d\nDIM_OF_WORLD: %d\n", Dim, Dim);
    fprintf(f, "number of vertices: %"dFMT"\n", nvert);
    fprintf(f, "number of elements: %"dFMT"\n", nelem);

    fprintf(f, "\nvertex coordinates:\n");
    for (i = 0; i < nvert; i++) {
#if Dim == 2
	fprintf(f, " %0.16lg %0.16lg\n", (double)g->verts[i][0],
					 (double)g->verts[i][1]);
#else
	fprintf(f, " %0.16lg %0.16lg %0.16lg\n", (double)g->verts[i][0],
			(double)g->verts[i][1], (double)g->verts[i][2]);
#endif
    }
    
    fprintf(f, "\nelement vertices:\n");
    phgTraverseElements(g, albert_vertices_callback);

    fprintf(f, "\nelement boundaries:\n");
    phgTraverseElements(g, albert_boundaries_callback);

#if Dim > 2
    fprintf(f, "\nelement type:\n");
    phgTraverseElements(g, albert_type_callback);
#endif

    fprintf(f, "\nelement neighbours:\n");
    phgTraverseElements(g, albert_neighbours_callback);

#if ALLOW_CURVED_BOUNDARY
    if (g->bdry_funcs != NULL) {
	EXPR **funcs;
	for (i = 0, funcs = g->bdry_funcs; *funcs != NULL; i++, funcs += 4);
	if (i > 0) {
	    fprintf(f, "\ncurved boundaries:\n%"dFMT"\n", i);
	    funcs = g->bdry_funcs;
	    while (*funcs != NULL) {
		fprintf(f, "%s;", phgDump3DFunction(funcs[0]));
		fprintf(f, " %s;", phgDump3DFunction(funcs[1]));
		fprintf(f, " %s;", phgDump3DFunction(funcs[2]));
		fprintf(f, " %s\n", phgDump3DFunction(funcs[3]));
		funcs += 4;
	    }
	}
    }
#endif

    phgFree(indices);
    ((f != stdout && f != stderr) ? fclose : fflush)(f);

    Return TRUE;
    return TRUE;	/* for avoiding MIPS cc warning */
}
