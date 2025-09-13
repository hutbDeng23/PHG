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

/* A test program for phgQuadInterfaceTriangle.
 *
 * $Id: quad_test3-triangle.c,v 1.49 2021/02/12 01:29:06 zlb Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "phg.h"

#define PI      _F(3.141592653589793238462643383279503)

typedef FLOAT TRIANGLE[3][2];	/* triangle defined by its three vertices */

/* To test a case:
 * 	% rm -f quad_test3-triangle quad_test3-triangle.o
 * 	% make -s USER_CFLAGS="-DTEST_CASE=#" quad_test3-triangle
 * where '#' is the case to test (1-8)
 */
#ifndef TEST_CASE
# define TEST_CASE 1
#endif

/*------------------------------------------------------------------------*/
#if TEST_CASE == 1
#define DESC	""
static TRIANGLE elist0[] = {};		/* empty list ==> read from file */
static FLOAT xc = 0.5, yc = 0.5, r = 0.25,
	     exact[] = {PI*.25*.25, 2.*PI*.25, 1.-PI*.25*.25};
/*------------------------------------------------------------------------*/
#elif TEST_CASE == 2
#define DESC	"Standard triangle, with Omega- entirely inside.\n"
static TRIANGLE elist0[] = { {{0.,0.},{1.,0.},{0.,1.}} };
static FLOAT xc = 0.25, yc = 0.25, r = .2, exact[]={PI*.04,.4*PI,.5-PI*.04};
/*------------------------------------------------------------------------*/
#elif TEST_CASE == 3
#define DESC	"Standard triangle entirely inside Omega-.\n"
static TRIANGLE elist0[] = { {{0.,0.},{1.,0.},{0.,1.}} };
static FLOAT xc = 0., yc = 0., r = 1.25, exact[]={.5,0.,0.};
/*------------------------------------------------------------------------*/
#elif TEST_CASE == 4
#define DESC	"A triangle intersecting with the interface.\n"
static TRIANGLE elist0[] = { {{0.5,1},{1,0.5},{0.5,0.5}} };
static FLOAT xc = 0.5, yc = 0.5, r = .25, exact[]={PI/64.,PI/8.,.125-PI/64.};
/*------------------------------------------------------------------------*/
#elif TEST_CASE == 5
#define DESC	"A triangle inside Omega+ with one vertex on the interface.\n"
static TRIANGLE elist0[] = { {{0.5,0.75},{0.5,0.765625},{0.484375,0.765625}} };
static FLOAT xc = 0.5, yc = 0.5, r = .25, exact[]={0.,0.,1.220703125e-04};
/*------------------------------------------------------------------------*/
#elif TEST_CASE == 6
#define DESC	"The unit square divided into two right triangles.\n"
static TRIANGLE elist0[]={{{0.,0.},{0.,1.},{1.,0.}},{{1.,1.},{1.,0.},{0.,1.}}};
static FLOAT xc = 0.5, yc = 0.5, r = 0.25,
	     exact[] = {PI*.25*.25, 2.*PI*.25, 1.-PI*.25*.25};
/*------------------------------------------------------------------------*/
#elif TEST_CASE == 7
#define DESC	"The standard triangle, and a circle centered at (0,0).\n"
static TRIANGLE elist0[] = { {{0.,0.},{0.,1.},{1.,0.}} };
static FLOAT xc = 0., yc = 0., r = .5, exact[]={PI/16.,PI/4.,.5-PI/16.};
/*------------------------------------------------------------------------*/
#elif TEST_CASE == 8
#define DESC	"A triangle, and a circle centered at (0,0).\n"
static TRIANGLE elist0[] = { {{.5,.5},{0.,.5},{.5,0.}} };
static FLOAT xc=0., yc=0., r=.5, exact[]={PI/16.-.125, PI/4., .25-PI/16.};
/*------------------------------------------------------------------------*/
#else
#error Invalid value for TEST_CASE!
#endif

#define DIM_MAX	4
static INT order = 15, type = -1, dim = 1;
static const char *proj_keys[] = {"none", "dot", "cross", NULL};
static PROJ    proj_vals[] = {PROJ_NONE, PROJ_DOT, PROJ_CROSS};
static int proj = 0;

/* The level set function and its gradient */
static int LS_order = 2;	/* polynomial order */
static void LS_func(FLOAT x, FLOAT y, FLOAT *value) {
    x -= xc; y -= yc; *value = x * x + y * y - r * r; }
static void LS_grad(FLOAT x, FLOAT y, FLOAT *grad) {
    x -= xc; y -= yc; grad[0] = x + x; grad[1] = y + y; }

/* The integrant */
static void
F(FLOAT x, FLOAT y, FLOAT *value)
/* F is defined such that proj(F(x,y)) = (1, 0.5, 0.25, ...) */
{
    int i;
    FLOAT n[2];

    assert(dim > 0 && dim <= DIM_MAX);

    if (proj == PROJ_NONE || type != 0) {
	value[0] = 1.0;
	for (i = 1; i < dim; i++)
	    value[i] = 0.5 * value[i - 1];
	return;
    }

    /* Get normal vector */
    LS_grad(x, y, n);
    x = 1.0 / Sqrt(n[0] * n[0] + n[1] * n[1]);
    if (proj == PROJ_DOT) {
	value[0] = n[0] * x;
	value[1] = n[1] * x;
    }
    else {
	value[0] = -n[1] * x;
	value[1] = n[0] * x;
    }
    for (i = 1; i < dim; i++) {
	value[2 * i] = 0.5 * value[2 * (i - 1)];
	value[2 * i + 1] = 0.5 * value[2 * (i - 1) + 1];
    }
}

static INT
read_mesh(const char *fn, TRIANGLE **pelist)
/* reads triangles into elist[] from the given mesh file and
 * returns # triangles */
{
    int i, nvert, nelem;
    char s[1024];
    FILE *f;
    int i0, i1, i2, i3;
    double a, b;
    FLOAT (*verts)[2];
    TRIANGLE *elems;

    if ((f = fopen(fn, "r")) == NULL)
	phgError(1, "Cannot open file \"%s\".\n", fn);

#define READ_ARG(fmt, var) \
	if (fscanf(f, fmt, var) != 1) \
	    phgError(1, "error reading file \"%s\".\n", fn)

#define FIND_KEY(key) \
    while (TRUE) { \
	READ_ARG("%s", s); \
	if (!strcmp(s, key)) \
	    break; \
    }

    fprintf(stderr, "Reading mesh from file \"%s\" ...\n", fn);

    FIND_KEY("Dimension");
    READ_ARG("%d", &i);
    if (i != 2)
	phgError(1, "Not a 2d mesh, abort.\n");

    FIND_KEY("Vertices");
    READ_ARG("%d", &nvert);
    fprintf(stderr, "number of vertices: %d\n", nvert);
    verts = phgAlloc(nvert * sizeof(*verts));
    for (i = 0; i < nvert; i++) {
	if (fscanf(f, "%lf %lf %d", &a, &b, &i0) != 3)
	    phgError(1, "error reading vertices.\n");
	verts[i][0] = a;
	verts[i][1] = b;
    }

    FIND_KEY("Triangles");
    READ_ARG("%d", &nelem);
    fprintf(stderr, "number of triangles: %d\n", nelem);
    elems = phgAlloc(nelem * sizeof(*elems));
    for (i = 0; i < nelem; i++) {
	if (fscanf(f, "%d %d %d %d", &i0, &i1, &i2, &i3) != 4)
	    phgError(1, "error reading triangles.\n");
	elems[i][0][0] = verts[i0 - 1][0];
	elems[i][0][1] = verts[i0 - 1][1];

	elems[i][1][0] = verts[i1 - 1][0];
	elems[i][1][1] = verts[i1 - 1][1];

	elems[i][2][0] = verts[i2 - 1][0];
	elems[i][2][1] = verts[i2 - 1][1];
    }

    phgFree(verts);
    *pelist = elems;

    fclose(f);

    return nelem;
}

int
main(int argc, char *argv[])
{
    INT i, n;
    double time0;
    TRIANGLE *elist;
    FLOAT res[DIM_MAX], d, *rule, **rules[3] = {NULL, NULL, NULL};
    /* Note: mesh file generated by
     *		~/src/MeshGeneration/triangle/quad_test3-triangle.sh */
    char *mesh_file = "quad_test3-triangle.mesh";

    phgOptionsRegisterInt("-type", "Quadrature type", &type);
    phgOptionsRegisterInt("-order", "Quadrature order", &order);
    phgOptionsRegisterInt("-dim", "Dimension of proj(u)", &dim);
    phgOptionsRegisterKeyword("-proj", "Projection type", proj_keys, &proj);
    phgOptionsRegisterFilename("-mesh_file", "Mesh filename", &mesh_file);

    phgInit(&argc, &argv);
    assert(phgNProcs == 1);

    type = type == 0 ? 0 : (type < 0 ? -1 : 1);
    if (type != 0)
	proj = 0;

    rules[type + 1] = &rule;

    phgPrintf("Test case: %d\n", TEST_CASE);
    if (DESC[0] != '\0')
	strchr(DESC, '%') == NULL ?
		phgPrintf(DESC) : phgPrintf(DESC, mesh_file);

    n = (INT)(sizeof(elist0) / sizeof(elist0[0]));
    if (n == 0)
	n = read_mesh(mesh_file, &elist);
    else
	elist = elist0;

    bzero(res, sizeof(res));

    time0 = phgGetTime(NULL);

#if USE_OMP && defined(_OPENMP) && _OPENMP >= 201511 /* OpenMP 4.5 */
    if (phgMaxThreads > 1)
	phgPrintf("OpenMP version: %d, parallel loop enabled, %d threads.\n",
			_OPENMP, phgMaxThreads);
# pragma omp parallel for reduction(+:res[:dim])
#endif	/* USE_OMP */
    for (i = 0; i < n; i++) {
	int ii;
	FLOAT res0[DIM_MAX];
#if TEST_CASE == 1
	ii = phgQuadInterfaceMarkElement(LS_func, LS_order, LS_grad,
					 ET_TRIANGLE, elist[i]);
	if (ii != 0 && ((ii > 0) != (type > 0)))
	    continue;
#else	/* TEST_CASE == 1 */
	ii = 0;
#endif	/* TEST_CASE == 1 */
	phgQuadInterfaceTriangle(ii == 0 ? LS_func : NULL, LS_order, LS_grad,
				 elist[i], order, rules[0], rules[1], rules[2]);
	phgQuadInterfaceRuleApply(F, dim, proj_vals[proj], rule, res0);
	phgFree(rule);
	for (ii = 0; ii < dim; ii++)
	    res[ii] += res0[ii];
    }

    time0 = phgGetTime(NULL) - time0;

    d = exact[type < 0 ? 0 : (type == 0 ? 1 : 2)]; 
    printf("Result: type = %d, quad_order = %d, dim = %d, proj = %s\n",
		type, order, dim, proj_keys[proj]);
    for (i = 0; i < dim; i++, d *= 0.5)
	printf("    res[%d] = %0.16lg, rel. error = %le\n",
			i, (double)res[i],
			fabs((double)((res[i] - d) / (d == 0. ? 1.0 : d))) );

    printf("Elasped time: %0.4lgs.\n", time0);

    if (elist != elist0)
	phgFree(elist);

    phgFinalize();

    return 0;
}
