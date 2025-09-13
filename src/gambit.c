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
 * $Id: gambit.c,v 1.16 2021/12/08 02:48:52 zlb Exp $
 */

#include "phg.h"
#include "phg/elem-info.h"
#include "phg/io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* Note: must match "phgElemInfo[]" */
static int Gambit_Elemtype[] = {
    -1,				/* Point */
    -1,				/* Edge */
    -1,				/* Triangle */
    -1,				/* Quadrilateral */
    6, 				/* Tetrahedron */
    7, 				/* Pyramid */
    5, 				/* Wedge */
    4				/* Hexahedron */
};

static int Gambit_Vertmap[][8] = {
    {},				/* Point */
    {},				/* Edge */
    {},				/* Triangle */
    {},				/* Quadrilateral */
    {0, 1, 2, 3},               /* Tetrahedron */
    {0, 1, 3, 2, 4}, 		/* Pyramid */
    {0, 1, 2, 3, 4, 5},		/* Wedge */
    {0, 1, 3, 2, 4, 5, 7, 6}	/* Hexahedron */
};

static int Gambit_Facemap[][8] = {
    {},				/* Point */
    {},				/* Edge */
    {},				/* Triangle */
    {},				/* Quadrilateral */
    {3, 2, 0, 1},               /* Tetrahedron */
    {4, 0, 1, 2, 3}, 		/* Pyramid */
    {0, 1, 2, 3, 4, 5},		/* Wedge */
    {0, 1, 2, 3, 4, 5, 6, 7}	/* Hexahedron */
};

static BOOLEAN
get_token(FILE *fp, char *token)
{
    int c;
    char *p;

    while (TRUE) {
	if (fscanf(fp, "%s", token) != 1)
	    return FALSE;
	if (token[0] != '#')
	    break;
	/* skip to newline */
	do
	    if ((c = fgetc(fp)) == EOF)
		return FALSE;
	while (c != '\n');
    }
    if ((p = strchr(token, '#')) != NULL)
	*p = '\0';
    return TRUE;
}

static GRID *g_;

static int
comp_face2(const void *p0, const void *p1)
/* compares two element faces: (elem index, face no, type)
 * acoording to first to component */
{
    INT (*f0)[3] = (void *)p0, (*f1)[3] = (void *)p1;
    int i;

    if ((i = (*f0)[0] - (*f1)[0]) != 0)
	return i;
    return (*f0)[1] - (*f1)[1];
}

static void
process_tetrahedra(GRID *g, INT nlist, INT (*list)[3], INT *vtypes)
/* assign boundary types (e->bound_type[]) using the list of triangles */
{
    /* Note: each face is stored as (elem, face, type) */
    INT i;
    INT key[3], (*found)[3];
    ELEMENT *e;
    int j;

    for (i = 0; i < g->nroot; i++) {
	e = g->roots + i;
	key[0] = i;
	/* element type was saved in mark */
	for (j = 0; j < NFace; j++) {
	    e->bound_type[j] = UNDEFINED;
	    key[1] = j;
	    /* adjust bdry type according to types of triangles */
	    if (nlist == 0)
		continue;
	    found = bsearch(key, list, nlist, sizeof(*list), comp_face2);
	    if (found != NULL) {
		e->bound_type[j] = _phg_bcmap(g_, (*found)[2]);
		/* printf(" ! set e:[%d,%d], %d:%d\n", e->index, j,  */
		/*        (*found)[3], _phg_bcmap(g_, (*found)[3])); */
	    }
	}
    }
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
process_element(GRID *g, int elem_type, int nelem, void *list_elem,
	     INT nlist, INT (*list)[3], INT *vtypes)
/* assign boundary types (e->bound_type[]) using the list of quadrilaterals */
{
    /* Note: each face is stored as (elem, face_no, type) */
    int nvert_type = phgElemInfo[elem_type].nvert;
    int nface_type = phgElemInfo[elem_type].nface;
    FUNC2TET func2tet = phgElemInfo[elem_type].func2tet;
    int (*liste)[nvert_type + 1];
    INT i, key[3], (*found)[3];
    int j, bound_type[nface_type];

    liste = list_elem;
    ntetra = ntetra_allocated = g->nroot;
    tetras = g->roots;
    for (i = 0; i < nelem; i++) {
	key[0] = i;
	for (j = 0; j < nface_type; j++) {
	    key[1] = j;
	    bound_type[j] = UNDEFINED;
	    /* adjust bdry type according to types of quadrilaterals */
	    if (nlist == 0)
		continue;
	    found = bsearch(key, list, nlist, sizeof(*list), comp_face2);
	    if (found != NULL)
		bound_type[j] = _phg_bcmap(g_, (*found)[2]);
	}
	/* convert to tetrahedra */
	region_mark = liste[i][nvert_type];
	func2tet(liste[i], new_tetra, bound_type);
    }
    phgInfo(2, "%d  %s converted to %d tetrahedra\n",
	    (int)nelem, phgElemInfo[elem_type].name, (int)(ntetra - g->nroot));
    g->ntree = g->nleaf = g->nelem = g->nroot = g->nelem_global = ntetra;
    g->roots = phgRealloc0(tetras, ntetra * sizeof(*g->roots));
}

int
phgImportGambit(GRID *g, const char *filename, BOOLEAN parallel)
{
    FILE *fp;
    char token[1024]; 
    INT i, j, n;
    INT *vtypes = NULL;			/* vertex types */
    INT nlistt = 0;			/* number of terahedra */
    INT nlistE[phgElemInfoCount], size_listE[phgElemInfoCount],
	nlistF[phgElemInfoCount], size_listF[phgElemInfoCount];
    void *listE[phgElemInfoCount];
    INT (*listF[phgElemInfoCount])[3];
    INT (*listEtype)[2] = NULL;	/* list of all element type */
    INT numnp, nelem, ngrps, nbsets;

    FunctionEntry;

    bzero(nlistE, sizeof(nlistE));
    bzero(size_listE, sizeof(size_listE));
    bzero(listE, sizeof(listE));
    bzero(nlistF, sizeof(nlistF));
    bzero(size_listF, sizeof(size_listF));
    bzero(listF, sizeof(listF));

#define READ_NUMBER						\
    if (!get_token(fp, token)) strcpy(token, "End ALL");	\
    if (isalpha((int)(token[0]))) {				\
	phgWarning("fewer entries (%"dFMT") than expected.\n", i);	\
	break;							\
    }

#define GET_TOKEN {				\
	if (!get_token(fp, token))		\
	    goto error;				\
    }

#define UNUSED_LINE { 				\
	if (fgets(token, 1024, fp) == NULL);	\
    }
#define NEXT_LINE UNUSED_LINE

#define UNUSED_LINE_CHECK(str) {					\
	if (fgets(token, 1024, fp) != NULL && !strcmp(token, str))	\
	    phgError(-1, "Unmatched phase when read mesh file,"		\
		     " line:%d!\n", __LINE__);				\
    }

    //	fprintf(stderr, "Unuesd: %s", token);

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

    /***********************/
    /* Control Information */
    /***********************/
    UNUSED_LINE_CHECK("        CONTROL INFO 2.3.16");
    UNUSED_LINE_CHECK("** GAMBIT NEUTRAL FILE");
    UNUSED_LINE; /* gambit session id */
    if (!get_token(fp, token) || strcmp(token, "PROGRAM:") ||
	!get_token(fp, token) || strcmp(token, "Gambit") ||
	!get_token(fp, token) || strcmp(token, "VERSION:") ||
	!get_token(fp, token) || strcmp(token, "2.3.16")) {
    error:
	phgError(1, "invalid or unsupported input file, abort.\n");
    }
    NEXT_LINE; 
    UNUSED_LINE; /* time when gambit file is created */
    UNUSED_LINE_CHECK("     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL");
    GET_TOKEN; numnp = atoi(token);
    GET_TOKEN; nelem = atoi(token);
    GET_TOKEN; ngrps = atoi(token);
    GET_TOKEN; nbsets = atoi(token);
    if (!get_token(fp, token) || strcmp(token, "3") ||
	!get_token(fp, token) || strcmp(token, "3"))
	goto error;
    NEXT_LINE; 
    UNUSED_LINE_CHECK("ENDOFSECTION");


    /*********************/
    /* Nodal Coordinates */
    /*********************/
    UNUSED_LINE_CHECK("   NODAL COORDINATES 2.3.16");
    n = numnp;
    phgInfo(2, "number of vertices: %"dFMT"\n", n);
    g->verts = phgNewVertices(n);
    for (i = 0; i < n; i++) {
	READ_NUMBER;
	GET_TOKEN; g->verts[i][0] = atof(token);
	GET_TOKEN; g->verts[i][1] = atof(token);
	GET_TOKEN; g->verts[i][2] = atof(token);
    }
    g->nvert_global = g->nvert = n;
    NEXT_LINE;
    UNUSED_LINE_CHECK("ENDOFSECTION"); 


    /*****************************/
    /* Element/Cell Connectivity */
    /*****************************/
    UNUSED_LINE_CHECK("      ELEMENTS/CELLS 2.3.16");
    n = nelem;
    phgInfo(2, "number of elements: %"dFMT"\n", n);
    listEtype = phgCalloc(nelem, sizeof(*listEtype));
    g->roots = phgReallocElements(g->roots, nlistt, n + nlistt);
    for (i = 0; i < n; i++) {
	int type, ndp, elem_type, nvert_type, nface_type;
	/* element index */
	READ_NUMBER;	
	/* element type */
	if (!get_token(fp, token)) 
	    goto error;
	type = atoi(token);
	/* element ndp */
	if (!get_token(fp, token)) 
	    goto error;
	ndp = atoi(token);

	for (j = 0; j < phgElemInfoCount; j++) 
	    if (type == Gambit_Elemtype[j]) 
		break;
	assert(j < phgElemInfoCount);
	elem_type = j;
	
	nvert_type = phgElemInfo[elem_type].nvert;
	nface_type = phgElemInfo[elem_type].nface;
	Unused(nface_type);
	Unused(ndp);
	assert(ndp == nvert_type);
	if (elem_type == ET_TETRA) {
	    GET_TOKEN; g->roots[nlistt].verts[0] = atoi(token) - 1;
	    GET_TOKEN; g->roots[nlistt].verts[1] = atoi(token) - 1;
	    GET_TOKEN; g->roots[nlistt].verts[2] = atoi(token) - 1;
	    GET_TOKEN; g->roots[nlistt].verts[3] = atoi(token) - 1;
	    listEtype[i][0] = elem_type;
	    listEtype[i][1] = nlistt;
	    nlistt++;
	} else {
	    INT (*liste)[nvert_type + 1], nliste, size_liste;
	    
	    liste = listE[elem_type];
	    nliste = nlistE[elem_type];
	    size_liste = size_listE[elem_type];
	    if (nliste >= size_liste) 
		listE[elem_type] = liste = 
		    phgRealloc0(liste, (size_listE[elem_type] += 1024) * sizeof(*liste));
	    for (j = 0; j < nvert_type; j++) {
		GET_TOKEN; liste[nliste][Gambit_Vertmap[elem_type][j]] = atoi(token) - 1;
		//printf("elem: %"dFMT"[%"dFMT"]: %d\n", i, j, atoi(token) - 1);
	    }
	    listEtype[i][0] = elem_type;
	    listEtype[i][1] = nliste;
	    nlistE[elem_type]++;
	}
    }
    NEXT_LINE;
    UNUSED_LINE_CHECK("ENDOFSECTION"); 

    /*****************************/
    /* Element Group Information */
    /*****************************/
    for (j = 0; j < ngrps; j++) {
	int nflags;
	UNUSED_LINE_CHECK("       ELEMENT GROUP 2.3.16");
	/* group index */
	if (!get_token(fp, token) || strcmp(token, "GROUP:") ||
	    !get_token(fp, token))
	    goto error;
	assert(j == atoi(token) - 1); 
	/* num of group record */
	if (!get_token(fp, token) || strcmp(token, "ELEMENTS:") ||
	    !get_token(fp, token))
	    goto error;
	n = atoi(token); 	
	/* materials and flags */
	if (!get_token(fp, token) || strcmp(token, "MATERIAL:") ||
	    !get_token(fp, token) || 
	    !get_token(fp, token) || strcmp(token, "NFLAGS:") ||
	    !get_token(fp, token))
	    goto error;
	nflags = atoi(token);
	NEXT_LINE;
	
	/* group name */
	UNUSED_LINE; 		
	/* solver dependent flags */
	for (i = 0; i < nflags; i++)
	    UNUSED_LINE; 	
	
	/* element in group */
	for (i = 0; i < n; i++) {
	    INT elem, elem_type, elem_idx;
   
	    GET_TOKEN;
	    elem = atoi(token) - 1;
	    elem_type = listEtype[elem][0];
	    elem_idx = listEtype[elem][1];
	    
	    if (elem_type == ET_TETRA) {
		g->roots[elem_idx].region_mark = j;
	    } else {
		int nvert_type = phgElemInfo[elem_type].nvert;
		INT (*liste)[nvert_type + 1];
		
		liste = listE[elem_type];
		liste[elem_idx][nvert_type] = j;
	    }
	}
	NEXT_LINE;
	UNUSED_LINE_CHECK("ENDOFSECTION"); 
    }


    /****************************/
    /* Boundary Conditions Sets */
    /****************************/
    for (j = 0; j < nbsets; j++) {
	UNUSED_LINE_CHECK(" BOUNDARY CONDITIONS 2.3.16");
	/* bdry cond name */
	GET_TOKEN;
	phgInfo(2, "boundary condition %s:", token);
	/* bdry cond type */
	GET_TOKEN;
	/* Note: elem:face type is allowed,
	 *       vert type is NOT allowed! */
	assert(1 == atoi(token));
	/* num of bdry cond record */
	GET_TOKEN;
	n = atoi(token);
	phgInfo(2, "%"dFMT"\n", n);
	/* bdry cond: value and code, unused */
	UNUSED_LINE;
	
	for (i = 0; i < n; i++) {
	    int elem, face, elem_type;
	    INT (*listf)[3], nlistf;
	    /* element */
	    READ_NUMBER;
	    elem = atoi(token) - 1; 
	    /* element type */
	    GET_TOKEN;
	    /* bdry face */
	    GET_TOKEN;
	    face = atoi(token) - 1;

	    elem_type = listEtype[elem][0];
	    listf = listF[elem_type];
	    nlistf = nlistF[elem_type];
	    assert(elem_type < phgElemInfoCount);

	    if (nlistf >= size_listF[elem_type])
		listF[elem_type] = listf = 
		    phgRealloc0(listf, (size_listF[elem_type] += 1024) * sizeof(*listf));
	    listf[nlistf][0] = listEtype[elem][1];
	    listf[nlistf][1] = Gambit_Facemap[elem_type][face];
	    listf[nlistf][2] = j;
	    nlistF[elem_type]++;
	}
	NEXT_LINE;
	UNUSED_LINE_CHECK("ENDOFSECTION"); 
    }
	
    phgCloseInputFile_(fp);

    g->nelem_global = g->ntree = g->nleaf = g->nelem = g->nroot = nlistt;

    for (j = 0; j < phgElemInfoCount; j++) {
	INT (*listf)[3], nlistf;
	listf = listF[j];
	nlistf = nlistF[j];
	
	if (nlistf > 0)
	    qsort(listf, nlistf, sizeof(*listf), comp_face2);
    }

    if (g->nroot > 0)
	process_tetrahedra(g, nlistF[ET_TETRA], listF[ET_TETRA], vtypes);

    for (j = 0; j < phgElemInfoCount; j++) {
	if (listE[j] != NULL) {
	    process_element(g, j, nlistE[j], listE[j], nlistF[j], listF[j], vtypes);
	    phgFree(listE[j]);
	    listE[j] = NULL;
	    nlistE[j] = 0;
	}
	if (listF[j] != NULL) {
	    phgFree(listF[j]);
	    listF[j] = NULL;
	    nlistF[j] = 0;
	}
    }

    if (vtypes != NULL) {
	phgFree(vtypes);
	vtypes = NULL;
    }

    phgFree(listEtype);

    Return 0;
    return 0;	/* for avoiding MIPS cc warning */
}
