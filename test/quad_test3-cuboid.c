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
 * Test program for phgQuadInterfaceCuboid and phgQuadInterfaceRectangle.
 * $Id: quad_test3-cuboid.c,v 1.58 2022/09/14 08:27:19 zlb Exp $ */

#include "phg.h"
#include "phg/sturm.h"
#include <math.h>

#undef HAVE_P4EST
#define HAVE_P4EST 0

#if HAVE_P4EST		/* Note: enabled with "make quad_test3-p4est" */
#include <p4est_vtk.h>
#include <p4est_iterate.h>

#include <p8est_vtk.h>
#include <p8est_iterate.h>
#endif	/* HAVE_P4EST */

/* The domain */
FLOAT xL = -1.0, yL = -2.0, zL = -3.0, xU = 1.0, yU = 2.0, zU = 3.0;

/* interface (level set) parameters */
static FLOAT c = .3;				/* radius or offset */
static FLOAT xc = .5, yc = .5, zc = .5;		/* center or normal vector */

static int ls_order = 2;
static INT dim = 3, type = -1;	/* space dimension and subdomain type */
static int proj = 0;		/* projection type (for type==0 only) */
static INT order = 15, n = 1;
static BOOLEAN marking = TRUE;

static const char *proj_names[] = {"none", "dot", /*"cross",*/ NULL};
static PROJ proj_ids[] = {PROJ_NONE, PROJ_DOT, PROJ_CROSS};

static FLOAT **rules[3] = {NULL, NULL, NULL}, *rule;

/******************************************************************************

		2D			3D
----------------------------------------------------------------------	ls_order
Domain: 	(xL,xU)x(yL,yU)		(xL,xU)x(yL,yU)x(zL,zU)
Levelset:	|(x-xc,y-yc)|^2 - c^2	|(x-xc,y-yc,z-zc)|^2 - c^2	2
		(x,y).(xc,yc)+c		(x,y,z).(xc,yc,zc)+c		1
Integrand:	(x+1)*(y+2)		(x+1)*(y+2)*(z+3)
----------------------------------------------------------------------

Sphere/circle
-------------

var('x,y,z,xc,yc,zc,c,r,theta,phi,xL,yL,zL,xU,yU,zU')

#------------- 2D
xy = [r*cos(theta)+xc, r*sin(theta)+yc]
Jacobian: factor(det(jacobian(xy, (r,theta)))) = r
dS: vector([derivative(xy[0],theta), derivative(xy[1],theta)]).norm() = r

f2(x,y) = (x+1)*(y+2)
whole: integrate(integrate(f2(xL+x*(xU-xL),yL+y*(yU-yL)), x,0,1), y,0,1) * (xU-xL) * (yU-yL) = 1/4*(xL + xU + 2)*(xL - xU)*(yL + yU + 4)*(yL - yU)
inside: integrate(integrate(f2(xy[0],xy[1])*r, r, 0, c), theta, 0, 2*pi) = 2*pi*c^2*xc + 2*pi*c^2 + (pi*c^2*xc + pi*c^2)*yc
curve: integrate(f2(xy[0],xy[1])*r, theta, 0, 2*pi)|_{r=c} = 2*(2*pi + 2*pi*xc + (pi + pi*xc)*yc)*c

#------------- 3D (theta: [0,2*pi], phi: [0,pi])
xyz = [r*cos(theta)*sin(phi)+xc, r*sin(theta)*sin(phi)+yc, r*cos(phi)+zc]

Jacobian: factor(det(jacobian(xyz, (r,theta,phi)))) = r^2*sin(phi)
dS:
def cross(u,v):
    return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
def D(u,r):
    return [derivative(u[0],r), derivative(u[1],r), derivative(u[2],r)]
factor(norm(vector(cross(D(xyz,theta), D(xyz,phi))))) = r^2*sin(phi)

f3(x,y,z) = (x+1)*(y+2)*(z+3)
whole: integrate(integrate(integrate(f3(xL+x*(xU-xL),yL+y*(yU-yL),zL+z*(zU-zL)),x,0,1),y,0,1),z,0,1) * (xU-xL) * (yU - yL) * (zU - zL) = -1/8*(xL + xU + 2)*(xL - xU)*(yL + yU + 4)*(yL - yU)*(zL + zU + 6)*(zL - zU)
inside: integrate(integrate(integrate(r^2*sin(phi)*f3(xyz[0],xyz[1],xyz[2]),r,0,c),theta,0,2*pi),phi,0,pi) = 4/3*pi*((yc*zc + 2*zc)*xc + yc*zc + 2*zc)*c^3 + 4*pi*(xc*yc + yc)*c^3 + 8*pi*c^3*xc + 8*pi*c^3
surface: integrate(integrate(r^2*sin(phi)*f3(xyz[0],xyz[1],xyz[2]),theta,0,2*pi),phi,0,pi)_{r=c} = pi*(xc*(yc + 2) + yc + 2)*c^2*((c^2 + 2*(c + 3)*zc + zc^2 + 6*c + 9)/c - (c^2 - 2*(c - 3)*zc + zc^2 - 6*c + 9)/c)

Plane/line
----------

var('x,y,z')

f2(x,y) = (x+1)*(y+2)
whole:	integrate(integrate(f2(x,y), x, -1, 1), y, -2, 2)	16
inside:	integrate(integrate(f2(x,y), x, -1, 1), y, -2, 0)	4
line:	integrate(f2(x,0), x, -1, 1)				4

f3(x,y,z) = (x+1)*(y+2)*(z+3)
whole:	integrate(integrate(integrate(f3(x,y,z), x, -1, 1), y, -2, 2), z, -3, 3)
	288
inside:	integrate(integrate(integrate(f3(x,y,z), x, -1, 1), y, -2, 2), z, -3, 0)
	72
plane:	integrate(integrate(f3(x,y,0), x, -1, 1), y, -2, 2)
	48

*******************************************************************************/

/* The 2D levelset function  */
static void
ls2(FLOAT x, FLOAT y, FLOAT *value)
{
    switch (ls_order) {
	case 1:
	    *value = x * xc + y * yc + c;
	    break;
	case 2:
	    x -= xc; y -= yc;
	    *value = x * x + y * y - c * c;
	    break;
    }
}

/* The gradient of the 2D levelset function */
static void
ls2_grad(FLOAT x, FLOAT y, FLOAT *grad)
{
    switch (ls_order) {
	case 1:
	    grad[0] = xc; grad[1] = yc;
	    break;
	case 2:
	    x -= xc; y -= yc;
	    grad[0] = x + x; grad[1] = y + y;
	    break;
    }
}

/* The 3D levelset function */
static void
ls3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    switch (ls_order) {
	case 1:
	    *value = x * xc + y * yc + z * zc + c;
	    break;
	case 2:
	    x -= xc; y -= yc; z -= zc;
	    *value = x * x + y * y + z * z - c * c;
	    break;
    }
}

/* The gradient of the 3D levelset function */
static void
ls3_grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *grad)
{
    switch (ls_order) {
	case 1:
	    grad[0] = xc; grad[1] = yc; grad[1] = zc;
	    break;
	case 2:
	    x -= xc; y -= yc; z -= zc;
	    grad[0] = x + x; grad[1] = y + y; grad[2] = z + z;
	    break;
    }
}

/* The 2D integrand */
static void
u2(FLOAT x, FLOAT y, FLOAT *f)
{
    FLOAT d;
    if (type != 0 || proj == 0) {
	*f = (x + 1.) * (y + 2.);
	return;
    }
    ls2_grad(x, y, f);
    d = Sqrt(f[0] * f[0] + f[1] * f[1]);
    assert(d > 1e-15);
    d = (x + 1.) * (y + 2.) / d;
    f[0] *= d;
    f[1] *= d;
}

/* The 3D integrand */
static void
u3(FLOAT x, FLOAT y, FLOAT z, FLOAT *f)
{
    FLOAT d;
    if (type != 0 || proj == 0) {
	*f = (x + 1.) * (y + 2.) * (z + 3.);
	return;
    }
    ls3_grad(x, y, z, f);
    d = Sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    assert(d > 1e-15);
    d = (x + 1.) * (y + 2.) * (z + 3.) / d;
    f[0] *= d;
    f[1] *= d;
    f[2] *= d;
}

/* check_element() checks the relative position of an element with respect to
 * the interface. It returns 0 if the interface might intersect the cuboid,
 * and +1/-1 if the element is in the +/- subdomain. */
static int
check_element(void *hypercube)
{
    int ret;

    if (!marking)
	return 0;

    if (dim == 2)
	ret = phgQuadInterfaceMarkElement(ls2, 2, ls2_grad,
				ET_RECTANGLE, hypercube);
    else
	ret = phgQuadInterfaceMarkElement(ls3, 2, ls3_grad,
				ET_CUBOID, hypercube);

    return ret;
}

static long long nevals, nelems, ntotal;
static FLOAT res, res0, ref;
#if USE_OMP
#pragma omp threadprivate(nevals, nelems, ntotal, res, res0, ref)
#endif  /* USE_OMP */

#if HAVE_P4EST
static int
refine_fn2(p4est_t *p4est, p4est_topidx_t which_tree, p4est_quadrant_t *quad)
/* returns !0 if the quadrant is to be refined */
{
    FLOAT rect[2][2];

#if FT_PHG == FT__DOUBLE
    p4est_qcoord_to_vertex(p4est->connectivity, which_tree,
					quad->x, quad->y, rect[0]);
#else	/* FT_PHG == FT__DOUBLE */
    double rect0[2][2];
    int i, j;
    p4est_qcoord_to_vertex(p4est->connectivity, which_tree,
					quad->x, quad->y, rect0[0]);
    for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
	    rect[i][j] = (FLOAT)rect0[i][j];
#endif	/* FT_PHG == FT__DOUBLE */

    rect[0][0] = xL + (xU - xL) * rect[0][0] / n;
    rect[0][1] = yL + (yU - yL) * rect[0][1] / n;
    rect[1][0] = rect[0][0] + (xU - xL) / (n * (FLOAT)(1 << quad->level));
    rect[1][1] = rect[0][1] + (yU - yL) / (n * (FLOAT)(1 << quad->level));

    return check_element(rect) == 0;
}

static int
refine_fn3(p8est_t *p8est, p4est_topidx_t which_tree, p8est_quadrant_t *quad)
/* returns !0 if the quadrant is to be refined */
{
    FLOAT cuboid[2][3];

#if FT_PHG == FT__DOUBLE
    p8est_qcoord_to_vertex(p8est->connectivity, which_tree,
					quad->x, quad->y, quad->z, cuboid[0]);
#else	/* FT_PHG == FT__DOUBLE */
    double cuboid0[2][3];
    int i, j;
    p8est_qcoord_to_vertex(p8est->connectivity, which_tree,
					quad->x, quad->y, quad->z, cuboid0[0]);
    for (i = 0; i < 2; i++)
	for (j = 0; j < 3; j++)
	    cuboid[i][j] = (FLOAT)cuboid0[i][j];
#endif	/* FT_PHG == FT__DOUBLE */

    cuboid[0][0] = xL + (xU - xL) * cuboid[0][0] / n;
    cuboid[0][1] = yL + (yU - yL) * cuboid[0][1] / n;
    cuboid[0][2] = zL + (zU - zL) * cuboid[0][2] / n;
    cuboid[1][0] = cuboid[0][0] + (xU - xL) / (n * (FLOAT)(1 << quad->level));
    cuboid[1][1] = cuboid[0][1] + (yU - yL) / (n * (FLOAT)(1 << quad->level));
    cuboid[1][2] = cuboid[0][2] + (zU - zL) / (n * (FLOAT)(1 << quad->level));

    return check_element(cuboid) == 0;
}

static void
callback_fn(void *info_ptr, void *user_data)
{
    FLOAT rect[2][2], cuboid[2][3], *elem;
    int np, flag;
    void *ls, *ls0, *u0, *ls0_grad;
    int (*function_ptr)(void *ls, int ls_order, void *ls_grad,
        void *, int quad_order, FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);

    ntotal++;

    if (dim == 2) {
	p4est_iter_volume_info_t *info = info_ptr;
	p4est_quadrant_t *quad = info->quad;
#if FT_PHG == FT__DOUBLE
	p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid,
					quad->x, quad->y, rect[0]);
#else	/* FT_PHG == FT__DOUBLE */
	double rect0[2][2];
	int i, j;
	p4est_qcoord_to_vertex(info->p4est->connectivity, info->treeid,
					quad->x, quad->y, rect0[0]);
	for (i = 0; i < 2; i++)
	    for (j = 0; j < 2; j++)
		rect[i][j] = (FLOAT)rect0[i][j];
#endif	/* FT_PHG == FT__DOUBLE */
	rect[0][0] = xL + (xU - xL) * rect[0][0] / n;
	rect[0][1] = yL + (yU - yL) * rect[0][1] / n;
	rect[1][0] = rect[0][0] + (xU - xL) / (n * (FLOAT)(1 << quad->level));
	rect[1][1] = rect[0][1] + (yU - yL) / (n * (FLOAT)(1 << quad->level));
	flag = check_element(rect);
	ls0 = ls2;
	ls0_grad = ls2_grad;
	u0 = u2;
	function_ptr = (void *)phgQuadInterfaceRectangle;
	elem = (FLOAT *)rect;
    }
    else {
	p8est_iter_volume_info_t *info = info_ptr;
	p8est_quadrant_t *quad = info->quad;
#if FT_PHG == FT__DOUBLE
	p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid,
					quad->x, quad->y, quad->z, cuboid[0]);
#else	/* FT_PHG == FT__DOUBLE */
	double cuboid0[2][3];
	int i, j;
	p8est_qcoord_to_vertex(info->p4est->connectivity, info->treeid,
					quad->x, quad->y, quad->z, cuboid0[0]);
	for (i = 0; i < 2; i++)
	    for (j = 0; j < 3; j++)
		cuboid[i][j] = (FLOAT)cuboid0[i][j];
#endif	/* FT_PHG == FT__DOUBLE */
	cuboid[0][0] = xL + (xU - xL) * cuboid[0][0] / n;
	cuboid[0][1] = yL + (yU - yL) * cuboid[0][1] / n;
	cuboid[0][2] = zL + (zU - zL) * cuboid[0][2] / n;
	cuboid[1][0] = cuboid[0][0] + (xU - xL) / (n * (FLOAT)(1<<quad->level));
	cuboid[1][1] = cuboid[0][1] + (yU - yL) / (n * (FLOAT)(1<<quad->level));
	cuboid[1][2] = cuboid[0][2] + (zU - zL) / (n * (FLOAT)(1<<quad->level));
	flag = check_element(cuboid);
	ls0 = ls3;
	ls0_grad = ls3_grad;
	u0 = u3;
	function_ptr = (void *)phgQuadInterfaceCuboid;
	elem = (FLOAT *)cuboid;
    }

    ls = ls0;
    if (flag != 0) {
	if (flag * type <= 0)
	    return;
	ls = NULL;
    }
    function_ptr(
		ls,		/* the level set function */
		ls_order,	/* poly. order */
		ls0_grad,	/* gradient of levelset */
		elem,		/* the element */
		order,		/* Gaussian quadrature order */
		rules[0], rules[1], rules[2]
    );
    np = phgQuadInterfaceRuleApply(u0, 1, proj_ids[proj], rule, &res0);
    phgFree(rule);

    res += res0;
    if (ls != NULL) {
	/* the element is on the interface */
	nelems++;
	nevals += np;
    }
}

static void
callback_fn2(p4est_iter_volume_info_t *info, void *user_data)
{
    callback_fn(info, user_data);
}

static void
callback_fn3(p8est_iter_volume_info_t *info, void *user_data)
{
    callback_fn(info, user_data);
}
#endif	/* HAVE_P4EST */

int
main(int argc, char *argv[])
{
    FLOAT PI = _F(3.14159265358979323846264338327950288);
    FLOAT error;
    int i;
    double t;
    char flt1[48], flt2[48], flt3[48];
    size_t mem_peak;

#if HAVE_P4EST
    INT refine = 0;
    p4est_connectivity_t *conn4;
    p4est_t *p4est;
    p8est_connectivity_t *conn8;
    p8est_t *p8est;
    BOOLEAN output_vtk = FALSE;
#else
    FLOAT dx, dy, dz;
    int j, k, np, flag;
    void *ls;
#endif	/* HAVE_P4EST */

    phgOptionsRegisterInt("-ls_order", "Polynomial order of levelset function",
				&ls_order);
    phgOptionsRegisterInt("-dim", "Space dimension (2 or 3)", &dim);
    phgOptionsRegisterInt("-order", "Gaussian quadrature order", &order);
    phgOptionsRegisterInt("-type",  "Quadrature type (-1, 0, 1)", &type);
    phgOptionsRegisterKeyword("-proj",  "Projection", proj_names, &proj);
    phgOptionsRegisterInt("-n", "Number of sub-intervals", &n);
    phgOptionsRegisterNoArg("-marking", "Whether mark element", &marking);

    phgOptionsRegisterFloat("-xc", "x coordinate of the center or normal", &xc);
    phgOptionsRegisterFloat("-yc", "y coordinate of the center or normal", &yc);
    phgOptionsRegisterFloat("-zc", "z coordinate of the center or normal", &zc);
    phgOptionsRegisterFloat("-c", "The radius or offset", &c);

    phgOptionsRegisterFloat("-xL", "The x coordinate of the lower corner", &xL);
    phgOptionsRegisterFloat("-yL", "The y coordinate of the lower corner", &yL);
    phgOptionsRegisterFloat("-zL", "The z coordinate of the lower corner", &zL);
    phgOptionsRegisterFloat("-xU", "The x coordinate of the upper corner", &xU);
    phgOptionsRegisterFloat("-yU", "The y coordinate of the upper corner", &yU);
    phgOptionsRegisterFloat("-zU", "The z coordinate of the upper corner", &zU);

#if HAVE_P4EST
    phgOptionsRegisterInt("-refine", "Refinement depth", &refine);
    phgOptionsRegisterNoArg("-output_vtk", "Output final mesh in VTK format",
				&output_vtk);
#endif

    phgInit(&argc, &argv);

    if (n == 1)
	marking = FALSE;

    /* check and normalize parameters */
    if (type < 0)
	type = -1;
    else if (type > 0)
	type = 1;
    if (type != 0)
	proj = 0;

    rules[type + 1] = &rule;

    assert(dim == 2 || dim == 3);
    assert(order >= 0);

    res = 0.0;
    nevals = nelems = ntotal = 0;

    if (dim == 2)
	goto rectangle;

    /*------------------------ cuboid --------------------------*/
#if HAVE_P4EST
    sc_init(phgComm, 0, 0, NULL, /*SC_LP_ESSENTIAL*/SC_LP_SILENT);
    p4est_init(NULL, SC_LP_PRODUCTION);
    conn8 = p8est_connectivity_new_brick(n, n, n, 0, 0, 0);
    p8est = p8est_new(phgComm, conn8, 0, NULL, NULL);
    for (i = 0; i < refine; i++) {
	int recursive = 0, partforcoarsen = 0;
	p8est_refine(p8est, recursive, refine_fn3, NULL);
	p8est_partition(p8est, partforcoarsen, NULL);
    }
    /* perform 1:2 balance */
    /*p8est_balance(p8est, P4EST_CONNECT_FACE, NULL); */
    if (output_vtk)
	p8est_vtk_write_file(p8est, NULL, "quad_test3-p8est");
    t = phgGetTime(NULL);
    p8est_iterate(p8est, NULL, NULL, callback_fn3, NULL, NULL, NULL);
    p8est_destroy(p8est);
    p8est_connectivity_destroy(conn8);
    sc_finalize();
#else	/* HAVE_P4EST */
    ntotal = (long long)n * (long long)n * (long long)n;
    t = phgGetTime(NULL);
    dx = (xU - xL) / n;
    dy = (yU - yL) / n;
    dz = (zU - zL) / n;
    for (i = 0; i < n; i++) {
	FLOAT cuboid[2][3];
	cuboid[1][0] = (cuboid[0][0] = xL + i * dx) + dx;
	for (j = 0; j < n; j++) {
	    cuboid[1][1] = (cuboid[0][1] = yL + j * dy) + dy;
	    for (k = 0; k < n; k++) {
		cuboid[1][2] = (cuboid[0][2] = zL + k * dz) + dz;
		ls = ls3;
		flag = check_element(cuboid);
		if (flag != 0) {
		    if (flag * type <= 0)
			continue;
		    ls = NULL;
		}
		phgQuadInterfaceCuboid(
			ls,		/* the level set function */
			ls_order,	/* poly. order */
			ls3_grad,	/* gradient of the ls func. */
			cuboid,		/* the element */
			order,		/* Gaussian quadrature order */
			rules[0], rules[1], rules[2]
		);
		np = phgQuadInterfaceRuleApply(u3,1,proj_ids[proj], rule,&res0);
		phgFree(rule);
		res += res0;
		if (ls != NULL) {
		    /* the element is on the interface */
		    nelems++;
		    nevals += np;
		}
	    }
	}
    }
#endif	/* HAVE_P4EST */
    goto cont;

rectangle:
    /*------------------------ rectangle --------------------------*/
#if HAVE_P4EST
    sc_init(phgComm, 0, 0, NULL, /*SC_LP_ESSENTIAL*/SC_LP_SILENT);
    p4est_init(NULL, SC_LP_PRODUCTION);
    conn4 = p4est_connectivity_new_brick(n, n, 0, 0);
    p4est = p4est_new(phgComm, conn4, 0, NULL, NULL);
    for (i = 0; i < refine; i++) {
	int recursive = 0, partforcoarsen = 0;
	p4est_refine(p4est, recursive, refine_fn2, NULL);
	p4est_partition(p4est, partforcoarsen, NULL);
    }
    /* perform 1:2 balance */
    /*p4est_balance(p4est, P4EST_CONNECT_FACE, NULL); */
    if (output_vtk)
	p4est_vtk_write_file(p4est, NULL, "quad_test3-p4est");
    t = phgGetTime(NULL);
    p4est_iterate(p4est, NULL, NULL, callback_fn2, NULL, NULL);
    p4est_destroy(p4est);
    p4est_connectivity_destroy(conn4);
    sc_finalize();
#else	/* HAVE_P4EST */
    ntotal = (long long)n * (long long)n;
    t = phgGetTime(NULL);
    dx = (xU - xL) / n;
    dy = (yU - yL) / n;
    for (i = 0; i < n; i++) {
	FLOAT rect[2][2];
	rect[1][0] = (rect[0][0] = xL + i * dx) + dx;
	for (j = 0; j < n; j++) {
	    rect[1][1] = (rect[0][1] = yL + j * dy) + dy;
	    ls = ls2;
	    flag = check_element(rect);
	    if (flag != 0) {
		if (flag * type <= 0)
		    continue;
		ls = NULL;
	    }
	    phgQuadInterfaceRectangle(
			ls,		/* the level set function */
			2,		/* poly. order */
			ls2_grad,	/* gradient of levelset */
			rect,		/* the element */
			order,		/* Gaussian quadrature order */
			rules[0], rules[1], rules[2]
	    );
	    np = phgQuadInterfaceRuleApply(u2, 1, proj_ids[proj], rule, &res0);
	    phgFree(rule);
	    res += res0;
	    if (ls != NULL) {
		/* the element is on the interface */
		nelems++;
		nevals += np;
	    }
	}
    }
#endif	/* HAVE_P4EST */
    /*-------------------------------------------------------------*/

cont:

#if HAVE_P4EST && USE_MPI
    if (phgNProcs > 1) {
	long long tmp0[3] = {ntotal, nelems, nevals}, tmp1[3];
	MPI_Reduce(tmp0, tmp1, 3, MPI_LONG_LONG, MPI_SUM, 0, phgComm);
	ntotal = tmp1[0];
	nelems = tmp1[1];
	nevals = tmp1[2];
	MPI_Reduce(&res, &res0, 1, PHG_MPI_FLOAT, PHG_SUM, 0, phgComm);
	res = res0;
    }
#endif	/* HAVE_P4EST && USE_MPI */

    t = phgGetTime(NULL) - t;

    if (dim == 2) {
	if (ls_order == 1 ||
	    xL > -c + xc || xU < c + xc || yL > -c + yc  || yU < c + yc) {
	    FLOAT rect[2][2] = {{xL, yL}, {xU, yU}};
	    phgOptionsPush();
	    phgOptionsSetOptions("-qi_dbg_off -qi_rect2tri -qi_threshold 0.95");
	    phgQuadInterfaceRectangle(ls2, 2, ls2_grad, rect,
				      (ls_order == 1 && type != 0) ? 3 : 55,
				      rules[0], rules[1], rules[2]);
	    phgQuadInterfaceRuleApply(u2, 1, proj_ids[proj], rule, &ref);
	    phgFree(rule);
	    phgOptionsPop();
	}
	else {
	    ref = type == 0 ? 2.*(2.*PI + 2.*PI*xc + (PI + PI*xc)*yc)*c :
			2.*PI*c*c*xc + 2.*PI*c*c + (PI*c*c*xc + PI*c*c)*yc;
	    if (type > 0)
		ref = .25 * (xL + xU + 2.) * (xL - xU) * (yL + yU + 4.)
			* (yL - yU) - ref;
	}
    }
    else {
	if (ls_order == 1 || xL > -c + xc || xU < c + xc || yL > -c + yc ||
	    yU < c + yc || zL > -c + zc || zU < c + zc) {
	    FLOAT cuboid[2][3] = {{xL, yL, zL}, {xU, yU, zU}};
	    phgOptionsPush();
	    phgOptionsSetOptions("-qi_dbg_off -qi_cuboid2tetra -qi_threshold 0.95");
	    phgQuadInterfaceCuboid(ls3, ls_order, ls3_grad, cuboid,
				   (ls_order == 1 && type != 0) ? 5 : 55,
				   rules[0], rules[1], rules[2]);
	    phgQuadInterfaceRuleApply(u3, 1, proj_ids[proj], rule, &ref);
	    phgFree(rule);
	    phgOptionsPop();
	}
	else {
	    ref = type == 0 ? PI*(xc*(yc + 2.) + yc + 2.)*c*c*((c*c + 2.*(c + 3.)*zc + zc*zc + 6.*c + 9.)/c - (c*c - 2.*(c - 3.)*zc + zc*zc - 6.*c + 9.)/c) :
			_F(4.)/_F(3.)*PI*((yc*zc + 2.*zc)*xc + yc*zc + 2.*zc)*c*c*c + 4.*PI*(xc*yc + yc)*c*c*c + 8.*PI*c*c*c*xc + 8.*PI*c*c*c;
	    if (type > 0)
		ref = -.125 * (xL + xU + 2.) * (xL - xU) * (yL + yU + 4.)
				* (yL - yU) * (zL + zU + 6.) * (zL - zU) - ref;
	}
    }

    error = ref == 0.0 ? Fabs(res) : Fabs(1. - res / ref); 
#if FT_PHG == FT___FLOAT128
	quadmath_snprintf(flt1, sizeof(flt1), "%0.34Qe", ref);
	quadmath_snprintf(flt2, sizeof(flt2), "%0.34Qe", res);
	quadmath_snprintf(flt3, sizeof(flt3), "%0.34Qe", error);
#else
	snprintf(flt1, sizeof(flt1), "%0.16le", (double)ref);
	snprintf(flt2, sizeof(flt2), "%0.16le", (double)res);
	snprintf(flt3, sizeof(flt3), "%0.16le", (double)error);
#endif
    phgMemoryUsage(NULL, &mem_peak);
    phgPrintf("Element type:       %s\n"
	      "Total elements:     %lld\n"
	      "Interface elements: %lld (%lld quadrature points)\n"
	      "Subdomain:          %s\n"
	      "Quadrature order:   %d\n"
	      "Reference value:    %s\n"
	      "Computed value:     %s\n"
	      "Relative error:     %s\n"
	      "Wall-clock time:    %lg\n"
	      "Memory usage:       %0.0lfMB\n",
		dim == 2 ? "Rectangle" : "Cuboid", ntotal, nelems, nevals,
		type == 0 ? "Gamma" : (type < 0 ? "Omega-" : "Omega+"),
		order, flt1, flt2, flt3, t,
		((double)mem_peak) / (1024. * 1024.));

    phgFinalize();

    return 0;
}
