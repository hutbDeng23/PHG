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

/* Numerical quadrature in a part of a triangle bounded by an implicitly
 * defined interface.
 *
 * Initial code by LIU Ziyang.
 *
 * $Id: quad-interface-triangle.c,v 1.115 2022/08/15 01:47:03 zlb Exp $ */

#include "phg.h"
#include "phg/sturm.h"

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <setjmp.h>

#undef Dim
#define Dim 2

#define INTER_MAX 20	/* 最大的交点个数 */
#define Opts	_phg_quad_interface_options
static FLOAT D3 = _F(1.0) / _F(3.0);

#define ONE _F(1.0)

typedef FLOAT TRIANGLE[Dim + 1][Dim];	/* triangle defined by three vertices */

enum {OM = 0, GM = 1, OP = 2};
static struct QuadInfo {
    FUNC_2D	ls;
    FUNC_2D	ls_grad;
    FLOAT	smax, smin, rmax, rmin, s, w_s;
    TRIANGLE	tri;
    QUAD	*quad;
    FLOAT	transform[4], reverse[4];

    BOOLEAN	reverse_flag, skip_rule;
    int		ls_order, quad_type, quad_order, level;

    /* variables for recording the quadrature rule */
    struct {FLOAT *pw, *nv; int n, alloc, trunc;} *rules[3];

    jmp_buf	jmp_env;
} QI;
#if USE_OMP
#pragma omp threadprivate(QI)
#endif /* USE_OMP */

#define LONGJMP if (Opts.subdiv_limit > 0) longjmp(QI.jmp_env, __LINE__)

static FLOAT
max(FLOAT x, FLOAT y)
{
    return x > y ? x : y;
}

static FLOAT
min(FLOAT x, FLOAT y)
{
    return x < y ? x : y;
}

static void
output_vtk(int line)
{
    FILE *fp;
    char fn[128];

    if (Opts.dbg_off || !Opts.dbg_vtk)
	return;
    sprintf(fn, "tmp%d-%d.vtk", QI.level, line);
    fp = fopen(fn, "wt");
    fprintf(fp, "# vtk DataFile Version 2.0\n"
		"vtk output\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n"
		"POINTS 3 float\n"
		"%0.16lg %0.16lg 0\n%0.16lg %0.16lg 0\n%0.16lg %0.16lg 0\n"
		"CELLS 1 4\n"
		"3 0 1 2\n"
		"CELL_TYPES 1\n"
		"5\n",
		(double)QI.tri[0][0], (double)QI.tri[0][1],
		(double)QI.tri[1][0], (double)QI.tri[1][1],
		(double)QI.tri[2][0], (double)QI.tri[2][1]);
    fprintf(stderr, "*** VTK file \"%s\" created.\n", fn);
}

static FLOAT
intersect_seg_seg(const FLOAT *pt1, const FLOAT *pt2,
		  const FLOAT *pt3, const FLOAT *pt4)
/* Returns t in [0,1] such that pt1+t*(pt2-pt1) is the intersection point of
 * the two segments pt1-pt2 and pt3-pt4, or 2.0 if they don't intersect. */
{
    /* 计算线段pt1-pt2与线段pt3-pt4的交点 */
    FLOAT a11, a12, a21, a22, b1, b2;
    FLOAT t, s, d, tol;

    /* We solve the following linear equations:
     *		pt1+t*(pt2-pt1) == pt3+s*(pt4-pt3)
     * for (t,s) in [0,1]^2 */
    a11 = pt2[0] - pt1[0], a12 = pt3[0] - pt4[0],
    a21 = pt2[1] - pt1[1], a22 = pt3[1] - pt4[1],
    b1  = pt3[0] - pt1[0],
    b2  = pt3[1] - pt1[1];

    d = a11 * a22 - a21 * a12;
    if (d == 0.0)
	return 2.0;

    tol = Opts.EPS * d;

    s = a11 * b2 - a21 * b1;
    if ((d > 0. && (s < -tol || s > d + tol)) ||
	(d < 0. && (s > -tol || s < d + tol)))
	return 2.0;

    t = a22 * b1 - a12 * b2;
    if ((d > 0. && (t < -tol || t > d + tol)) ||
	(d < 0. && (t > -tol || t < d + tol)))
	return 2.0;

    t /= d;
    if (t < 0.)
	t = 0.;
    else if (t > 1.)
	t = 1.;

    return t;
}

static int
intersect_seg_triangle(const FLOAT *pt1, const FLOAT *pt2, FLOAT t[3])
/* Returns the number of intersections of the segment pt1-pt2 with the
 * triangle. The positions of the intersections, sorted in ascending order,
 * are returned in t[] (pt1+t*(pt2-pt1)) */
{
    /* 计算线段pt1-pt2与三角形triangled的交点； */
    int i, n;

    t[0] = intersect_seg_seg(pt1, pt2, QI.tri[0], QI.tri[1]);
    t[1] = intersect_seg_seg(pt1, pt2, QI.tri[0], QI.tri[2]);
    t[2] = intersect_seg_seg(pt1, pt2, QI.tri[1], QI.tri[2]);
    qsort(t, 3, sizeof(t[0]), phgCompFLOAT);
    /* count # intersections, remove duplicate intersection points */
    for (n = 0; n < 3 && t[n] < 1. + Opts.EPS; n++) {
	if (n == 0 || Fabs(t[n] - t[n - 1]) > Opts.EPS * 100.0)
	    continue;
	/* delete current entry "n", which is identical to the entry "n-1" */
	for (i = n; i < 3 - 1; i++)
	    t[i] = t[i + 1];
	t[i] = 2.0;
	n--;
    }
    if (n == 3) {
	/* Check the case the segment is almost aligned with an edge,
	 * i.e., 2 of the intersection points are on vertices of the triangle.
	 * Some examples of this case are given by:
	 *	FLOAT tri[][2] = {{0.,0.}, {1.,0.}, {0.,1.}};
	 *
	 *	FLOAT pt1[] = {3.586020369539256e-14, 1.000000000000036};
	 *	FLOAT pt2[] = {0.9999999999997176, -2.824962486158711e-13};
	 *
	 *	FLOAT pt1[] = {-1.448036540806387e-13, 1};
	 *	FLOAT pt2[] = {2.255942614377446e-12, -2.400746268446905e-12};
	 * A test code for the case is available at the end of this file.
	 */
	BOOLEAN flag[3];	/* records whether point i is on a vertex */
	FLOAT pt[Dim];
	int j;
	/* set n := number of intersection points which are a vertex */
	n = 0;
	for (i = 0; i < 3; i++) {
	    pt[0] = (1.0 - t[i]) * pt1[0] + t[i] * pt2[0];
	    pt[1] = (1.0 - t[i]) * pt1[1] + t[i] * pt2[1];
	    flag[i] = FALSE;
	    for (j = 0; j < 3; j++) {
		if (Fabs(pt[0] - QI.tri[j][0]) > Opts.EPS * 100.0 ||
		    Fabs(pt[1] - QI.tri[j][1]) > Opts.EPS * 100.0)
		    continue;
		flag[i] = TRUE;
		n++;
		break;
	    }
	}
	if (n != 2) {
	    output_vtk(__LINE__);
	    phgInfo(-1, "unexpected error:\n"
			"FLOAT pt1[] = {%0.16lg, %0.16lg};\n"
		    	"FLOAT pt2[] = {%0.16lg, %0.16lg};\n"
		    	"          t = {%0.16lg, %0.16lg, %0.16lg};\n"
			"You may try \"-qi_eps\" with a larger or smaller "
			"value (current is %0.16le).\n",
		    	pt1[0], pt1[1], pt2[0], pt2[1], t[0], t[1], t[2],
			(double)Opts.EPS);
	    phgError(1, "%s:%d: abort.\n", __FILE__, __LINE__);
	}
	/* remove the point which is not a vertex (flag[i] == FALSE) */
	for (i = 0; i < 3 && flag[i]; i++);
	for (; i < 3 - 1; i++)
	    t[i] = t[i + 1];
	n = 2;
    }
    return n;
}

static void
unitization(FLOAT *r)
{
    /* 单位化 */
#if 1
    FLOAT d = Sqrt(r[0] * r[0] + r[1] * r[1]);
    assert(d /*>= Opts.EPS*/ != 0.0);
    d = 1.0 / d;
#else
    FLOAT d = 1.0 / Sqrt(r[0] * r[0] + r[1] * r[1]);
#endif
    r[0] *= d;
    r[1] *= d;
}

static void
update_QI(void)
{
    FLOAT x1, y1, x2, y2, x0, y0;

    /* 寻找积分方向 */
    QI.ls_grad( (QI.tri[0][0] + QI.tri[1][0] + QI.tri[2][0]) * D3,
		(QI.tri[0][1] + QI.tri[1][1] + QI.tri[2][1]) * D3,
		QI.transform);	/* s方向坐标 */
    if (QI.transform[0] == 0.0 && QI.transform[1] == 0.0)
	/* shouldn't happen if the interface is well resolved */
#if 1
	QI.transform[0] = QI.transform[1] = Sqrt(_F(0.5));
#else
	LONGJMP;
#endif
    else
	unitization(QI.transform);
    QI.transform[2] = -QI.transform[1];
    QI.transform[3] = QI.transform[0];	/* t方向坐标 */

    QI.reverse_flag = FALSE;	/* force recomputing reverse */

    /* 找矩形把三角形区域框起来 */
    x0 = QI.transform[0] * QI.tri[0][0] +
	 QI.transform[1] * QI.tri[0][1];
    y0 = QI.transform[2] * QI.tri[0][0] +
	 QI.transform[3] * QI.tri[0][1];
    x1 = QI.transform[0] * QI.tri[1][0] +
	 QI.transform[1] * QI.tri[1][1];
    y1 = QI.transform[2] * QI.tri[1][0] +
	 QI.transform[3] * QI.tri[1][1];
    x2 = QI.transform[0] * QI.tri[2][0] +
	 QI.transform[1] * QI.tri[2][1];
    y2 = QI.transform[2] * QI.tri[2][0] +
	 QI.transform[3] * QI.tri[2][1];
    QI.rmax = max(x0, max(x1, x2));
    QI.rmin = min(x0, min(x1, x2));
    QI.smax = max(y0, max(y1, y2));
    QI.smin = min(y0, min(y1, y2));

    if (Opts.dbg_off || !Opts.show_directions)
	return;

    fprintf(stderr, "Triangle: [[%0.16lg,%0.16lg],[%0.16lg,%0.16lg],"
		    "[%0.16lg,%0.16lg]], level %d\n"
		    "Integration direction: [%0.16lg,%0.16lg]\n"
		    "BoundingBox: [%0.16lg,%0.16lg]-[%0.16lg,%0.16lg]\n",
			(double)QI.tri[0][0], (double)QI.tri[0][1],
			(double)QI.tri[1][0], (double)QI.tri[1][1],
			(double)QI.tri[2][0], (double)QI.tri[2][1], QI.level,
			(double)QI.transform[0], (double)QI.transform[1],
			(double)QI.rmin, (double)QI.smin,
			(double)QI.rmax, (double)QI.smax);
}

static void
add_point(FLOAT x, FLOAT y, const FLOAT *n, const FLOAT w)
{
    int i, k, pos;

    /* debugging quad_type */
    if (QI.quad_type < -1 || QI.quad_type > 1) {
	phgError(1, "unexpected!\n");
    }

    i = QI.quad_type  + 1;

    if (QI.skip_rule || QI.rules[i] == NULL)
	return;

    if (QI.rules[i]->n >= QI.rules[i]->alloc) {
	k = QI.rules[i]->alloc + QI.rules[i]->trunc;
	QI.rules[i]->pw = phgRealloc_(QI.rules[i]->pw,
			(RULE_HEADER + k * (size_t)(Dim + 1))
				       * sizeof(*QI.rules[i]->pw),
			(RULE_HEADER + QI.rules[i]->n * (size_t)(Dim + 1))
				       * sizeof(*QI.rules[i]->pw));
	if (QI.quad_type == 0 && Opts.nv_flag)
	    QI.rules[i]->nv = phgRealloc_(QI.rules[i]->nv,
			k * Dim * sizeof(*QI.rules[i]->nv),
			QI.rules[i]->n * Dim * sizeof(*QI.rules[i]->nv));
	QI.rules[i]->alloc = k;
    }

    pos = RULE_HEADER + QI.rules[i]->n * (Dim + 1);

    /* the coordinates */
    QI.rules[i]->pw[pos++] = x;
    QI.rules[i]->pw[pos++] = y;

    /* the weight */
    QI.rules[i]->pw[pos++] = w;

    /* the unit normal vector */
    if (QI.quad_type == 0 && Opts.nv_flag) {
	pos = QI.rules[i]->n * Dim;
	QI.rules[i]->nv[pos++] = n[0];
	QI.rules[i]->nv[pos++] = n[1];
    }

    QI.rules[i]->n++;
}

static void
quad_triangle(FLOAT *res)
{
    /* 计算三角形上u(x,y)的积分 */
    QUAD *quad;
    FLOAT *p, *w, s, d;	/* Note: proj is not needed here */
    int j;
 
    quad = phgQuadGetQuad2D(QI.quad_order);
    p = quad->points;
    w = quad->weights;

    if (res != NULL)
	*res = 0.;

    /* area of the triangle = |det([x0,y0,1; x1,y1,1; x2,y2,1])| / 2! =
	Fabs(x0*y1 + x1*y2 + x2*y0 - y0*x1 - y1*x2 - y2*x0) / 2! */
    s = Fabs(QI.tri[0][0] * QI.tri[1][1] + QI.tri[1][0] * QI.tri[2][1] +
	     QI.tri[2][0] * QI.tri[0][1] - QI.tri[0][1] * QI.tri[1][0] -
	     QI.tri[1][1] * QI.tri[2][0] - QI.tri[2][1] * QI.tri[0][0]) * .5;

    for (j = 0; j < quad->npoints; j++, p += 3, w++) {
	d = s * *w;
	add_point(p[0] * QI.tri[0][0] +
		  p[1] * QI.tri[1][0] +
		  p[2] * QI.tri[2][0],
		  p[0] * QI.tri[0][1] +
		  p[1] * QI.tri[1][1] +
		  p[2] * QI.tri[2][1], NULL, d);
	if (res != NULL)
	    *res += ONE * d;
    }
}

static void
shorten_segment(const FLOAT *v0, const FLOAT *v1, FLOAT *t0, FLOAT *t1)
{
    /* 判断v0-v1哪一部分在三角形内部 */
    int n;
    FLOAT t[3];

    n = intersect_seg_triangle(v0, v1, t);
    assert(n == 1 || n == 2);
    switch (n) {
	case 0: *t0 = 1.0; *t1 = 0.0;	break;
	case 1: *t0 = *t1 = t[0];	break;
	case 2: *t0 = t[0]; *t1 = t[1]; break;
    }
}

typedef struct {
    const FLOAT *v0, *v1;
    FLOAT dv[Dim];		/* = v1 - v0 */
} NS_CONTEXT;

static BOOLEAN
newton_func(const FLOAT *t, FLOAT *f, void *pctx)
{
    NS_CONTEXT *nc = pctx;

    QI.ls(nc->v0[0] + (*t) * nc->dv[0], nc->v0[1] + (*t) * nc->dv[1], f);

    return TRUE;
}

static BOOLEAN
newton_jac(const FLOAT *t, FLOAT *jac, void *pctx)
{
    NS_CONTEXT *nc = pctx;
    FLOAT f[Dim];

    QI.ls_grad(nc->v0[0] + (*t) * nc->dv[0],
		    nc->v0[1] + (*t) * nc->dv[1], f);
    *jac = f[0] * nc->dv[0] + f[1] * nc->dv[1];

    return TRUE;
}

static int
intersect_edge(const FLOAT *v0, const FLOAT *v1, BOOLEAN shorten_flag,
	       FLOAT **roots, int **mults)
{
/* 计算线段v0-v1与ls的交点，shorten_flag=TRUE时，只计算在三角形内部的部分 */
    int i, j, order, n;
    FLOAT t0 = 0., t1 = 1., t, *c;
    FLOAT a, b, fa, fb, f;
    BOOLEAN use_newton;
    NS_CONTEXT ctx;

    if (shorten_flag)
	shorten_segment(v0, v1, &t0, &t1);

    if ((t0 -= Opts.EPS) < 0.0)
	t0 = 0.0;
    if ((t1 += Opts.EPS) > 1.0)
	t1 = 1.0;

    if (t0 >= t1) {
	*roots = NULL;
	*mults = NULL;
	return 0;
    }

    order = QI.ls_order;
    use_newton = (Opts.newton == 1 || order < 0);
    order = Opts.newton_porder > 0 ? Opts.newton_porder :
				     (order < 0 ? -order : order);
    *roots = c = phgAlloc((order + 1) * sizeof(*c));
    for (i = 0; i <= order; i++) {
	t = t0 + i / (FLOAT)order *(t1 - t0);
	QI.ls( v0[0] + t * (v1[0] - v0[0]),
		v0[1] + t * (v1[1] - v0[1]), &c[i] );
    }
    n = phgUnivariatePolyRoots(order, t0, t1, Opts.EPS, c, mults);
    ctx.v0 = v0;
    ctx.v1 = v1;
    ctx.dv[0] = v1[0] - v0[0];
    ctx.dv[1] = v1[1] - v0[1];
    j = -1;
    for (i = 0; i < n; i++) {
	int its;
	FLOAT residual;

	if ((*mults)[i] % 2 == 0)
	    continue;

	t = c[i];
	if (use_newton != TRUE) {
	    j++;
	    c[j] = t;
	    (*mults)[j] = 1;
	    continue;
	}

	if (j >= 0 && t <= c[j])
	    continue;

	its = phgNonlinearSolver(1, newton_func, newton_jac, &t,
				 Opts.newton_maxits, Opts.EPS, Opts.EPS,
				 &residual, &ctx);
	if (its >= 0 && its < Opts.newton_maxits) {
	    /* check for duplicate roots */
	    int k;
	    for (k = 0; k <= j; k++)
		if (Fabs(t - c[k]) < Opts.EPS * 10.0)
		    break;
	    if (k <= j)
		continue;	/* skip duplicate root */
	    j++;
	    c[j] = t;
	    (*mults)[j] = 1;
	}
	else {
	    /* try bisection */
	    a = (j < 0 ? t0 : (c[i] + c[j]) * 0.5);
	    b = (i == n - 1 ? t1 : (c[i] + c[i + 1]) * 0.5);
	    assert(a <= c[i] && c[i] <= b);
	    newton_func(&a, &fa, &ctx);
	    newton_func(&b, &fb, &ctx);
	    if (fa * fb > 0.0)
		continue;	/* no root */
	    if (Fabs(fa) == 0.0) {
		t = a;
	    }
	    else if (Fabs(fb) == 0.0) {
		t = b;
	    }
	    else
		while (TRUE) {
		    t = (a + b) * 0.5;
		    if (t == a || t == b || b - a <= Opts.EPS)
			break;
		    newton_func(&t, &f, &ctx);
		    if (Fabs(f) == 0.0)
			break;
		    if ((f > 0.0) == (fa > 0.0)) {
			a = t;
			fa = f;
		    }
		    else {
			b = t;
			fb = f;
		    }
		}
	    j++;
	    c[j] = t;
	    (*mults)[j] = 1;
	}
    }
    n = j + 1;

    return n;
}

static void
get_reverse_transform(const FLOAT *a, FLOAT *reverse)
{
    /* 矩阵求逆 */
    FLOAT d, det = a[0] * a[3] - a[1] * a[2];
    if (Fabs(det) < Opts.EPS) {
	printf("坐标系转换有误\n");
	exit(0);
    }
    d = 1.0 / det;
    reverse[0] = a[3] * d;
    reverse[1] = -a[1] * d;
    reverse[2] = -a[2] * d;
    reverse[3] = a[0] * d;
}

static void
xy2rs(FLOAT x, FLOAT y, FLOAT *r, FLOAT *s, const FLOAT *transform)
{
    /*（x,y）坐标转化为(r,s)坐标 */
    if (r != NULL)
	*r = x * transform[0] + y * transform[1];
    if (s != NULL)
	*s = x * transform[2] + y * transform[3];
}

static void
rs2xy(FLOAT r, FLOAT s, FLOAT *x, FLOAT *y)
{
    /*（r,s）坐标转化为(x,y)坐标 */
    if (!QI.reverse_flag) {
	get_reverse_transform(QI.transform, QI.reverse);
	QI.reverse_flag = TRUE;
    }
    xy2rs(r, s, x, y, QI.reverse);
}

static void
func_f(FLOAT r, FLOAT *res, const FLOAT wgt, void *ctx)
{
    FLOAT x, y;

    Unused(ctx);
    rs2xy(r, QI.s, &x, &y);
    add_point(x, y, NULL, wgt * QI.w_s);
    if (res != NULL)
	*res = ONE;

    return;
}

static void
func_g(FLOAT s, FLOAT res[3], const FLOAT wgt, void *ctx)
{
    FLOAT res0, inter[INTER_MAX + 2];
    FLOAT rpt1[Dim] = { QI.rmin, s }, rpt2[Dim] = { QI.rmax, s };
    FLOAT pt1[Dim], pt2[Dim], r, rmin, rmax, x, y, *roots;
    FLOAT t[3];
    int i, n, *mults;

    Unused(ctx);

    if (res != NULL)
	res[0] = res[1] = res[2] = 0.0;

    rs2xy(rpt1[0], rpt1[1], &pt1[0], &pt1[1]);
    rs2xy(rpt2[0], rpt2[1], &pt2[0], &pt2[1]);
    i = intersect_seg_triangle(pt1, pt2, t);
    assert(i == 1 || i == 2);
    if (i <= 1)
	return;
    xy2rs(pt1[0] + t[0] * (pt2[0] - pt1[0]),
	   pt1[1] + t[0] * (pt2[1] - pt1[1]), &rmin, NULL, QI.transform);
    xy2rs(pt1[0] + t[1] * (pt2[0] - pt1[0]),
	   pt1[1] + t[1] * (pt2[1] - pt1[1]), &rmax, NULL, QI.transform);
    if (rmin > rmax) {
	r = rmin;
	rmin = rmax;
	rmax = r;
    }
    if (rmax - rmin <= Opts.EPS + Opts.EPS)
	return;

    n = 0;
    inter[n++] = rmin;
    i = intersect_edge(pt1, pt2, TRUE, &roots, &mults);
    if (i >= 2 || (i == 1 && mults[0] > 1)) {
	phgFree(roots);
	phgFree(mults);
	LONGJMP;
    }
    phgFree(mults);
    if (i > 0) {
	xy2rs( pt1[0] + roots[0] * (pt2[0] - pt1[0]),
		pt1[1] + roots[0] * (pt2[1] - pt1[1]), &r, NULL, QI.transform);
	if (r > rmin + Opts.EPS && r < rmax - Opts.EPS)
	    inter[n++] = r;
    }
    inter[n++] = rmax;
    phgFree(roots);

    QI.s = s;
    QI.w_s = wgt;
    for (i = 0; i < n - 1; i++) {
	FLOAT ls;
	/* the interval is entirely either inside or outside the subdomain,
	 * we simply check the middle point of the interval to decide whether
	 * to skip it. */
	rs2xy(0.5 * (inter[i] + inter[i + 1]), s, &x, &y);
	QI.ls(x, y, &ls);
	QI.quad_type = ls <= 0. ? -1 : 1;
	phgQuad1D(func_f, 1, inter[i], inter[i + 1], QI.quad,
		  res == NULL ? NULL : &res0, NULL);
	if (!Opts.dbg_off && Opts.show_intervals_t && !Opts.dbg_t)
	    printf("t_interval: s=%g, i=%d, (%0.16lg,%0.16lg), res0=%0.16lg\n",
			(double)s, i, (double)inter[i], (double)inter[i + 1],
			(double)res0);
	if (res != NULL)
	    res[QI.quad_type + 1] += res0;
    }

    return;
}

static void
func_g_curve(FLOAT s, FLOAT *res, const FLOAT wgt, void *ctx)
{
    int *mults, n;
    FLOAT v0[Dim], v1[Dim], *roots, pt[Dim];

    Unused(ctx);

    rs2xy(QI.rmin, s, &v0[0], &v0[1]);
    rs2xy(QI.rmax, s, &v1[0], &v1[1]);
    n = intersect_edge(v0, v1, TRUE, &roots, &mults);
    if (n <= 0) {
	if (res != NULL)
	    *res = 0.;
	phgFree(roots);
	phgFree(mults);
	return;
    }
    else if (n >= 2 || mults[0] > 1) {
	phgFree(roots);
	phgFree(mults);
	LONGJMP;
    }
    else {
	FLOAT grad[Dim], gradrs[Dim], d;
	pt[0] = v0[0] + (*roots) * (v1[0] - v0[0]);
	pt[1] = v0[1] + (*roots) * (v1[1] - v0[1]);
	QI.ls_grad(pt[0], pt[1], grad);
	if (!QI.reverse_flag) {
	    get_reverse_transform(QI.transform, QI.reverse);
	    QI.reverse_flag = TRUE;
	}
	gradrs[0] = grad[0] * QI.reverse[0] + grad[1] * QI.reverse[2];
	gradrs[1] = grad[0] * QI.reverse[1] + grad[1] * QI.reverse[3];
	d = Sqrt(1. + (gradrs[1] * gradrs[1]) / (gradrs[0] * gradrs[0]));
	unitization(grad);
	QI.quad_type = 0;
	add_point(pt[0], pt[1], grad, d * wgt);
	QI.quad_type = 99;
	*res = ONE * d;
    }
    phgFree(mults);
    phgFree(roots);
    return;
}

/*******************************************************************************************************************************************/

static void quad_interface_triangle(FLOAT *res);

static void
subdivide_triangle(FLOAT res[3])
{
    FLOAT pt1[Dim] = {(QI.tri[0][0] + QI.tri[1][0]) * 0.5,
		      (QI.tri[0][1] + QI.tri[1][1]) * 0.5};
    FLOAT pt2[Dim] = {(QI.tri[0][0] + QI.tri[2][0]) * 0.5,
		      (QI.tri[0][1] + QI.tri[2][1]) * 0.5};
    FLOAT pt3[Dim] = {(QI.tri[2][0] + QI.tri[1][0]) * 0.5,
		      (QI.tri[2][1] + QI.tri[1][1]) * 0.5};
    FLOAT res0[3];
    TRIANGLE triangles[4];
    int i, j;

    QI.level++;

    triangles[0][0][0] = QI.tri[0][0];
    triangles[0][0][1] = QI.tri[0][1];
    triangles[0][1][0] = pt1[0];
    triangles[0][1][1] = pt1[1];
    triangles[0][2][0] = pt2[0];
    triangles[0][2][1] = pt2[1];

    triangles[1][0][0] = QI.tri[1][0];
    triangles[1][0][1] = QI.tri[1][1];
    triangles[1][1][0] = pt1[0];
    triangles[1][1][1] = pt1[1];
    triangles[1][2][0] = pt3[0];
    triangles[1][2][1] = pt3[1];

    triangles[2][0][0] = QI.tri[2][0];
    triangles[2][0][1] = QI.tri[2][1];
    triangles[2][1][0] = pt3[0];
    triangles[2][1][1] = pt3[1];
    triangles[2][2][0] = pt2[0];
    triangles[2][2][1] = pt2[1];

    triangles[3][0][0] = pt3[0];
    triangles[3][0][1] = pt3[1];
    triangles[3][1][0] = pt1[0];
    triangles[3][1][1] = pt1[1];
    triangles[3][2][0] = pt2[0];
    triangles[3][2][1] = pt2[1];

    memcpy(QI.tri, triangles[0], sizeof(QI.tri));
    quad_interface_triangle(res);
    for (j = 1; j < 4; j++) {
	memcpy(QI.tri, triangles[j], sizeof(QI.tri));
	quad_interface_triangle(res == NULL ? NULL : res0);
	if (res != NULL)
	    for (i = 0; i < 3; i++)
		res[i] += res0[i];
    }

    QI.level--;

    return;
}

static void
quad_interface_triangle(FLOAT *res)
{
/* 计算三角形triangle和曲线ls相交区域上的积分 */
    FLOAT s, *roots, thres;
    int i, j, k, n, rule_n[3], *mults;

    /* variables for recording intersections */
    int sinter_n;
    /* use a static buffer to save memory with regard to recursion */
    static FLOAT sinter[INTER_MAX];
#if USE_OMP
#pragma omp threadprivate(sinter)
#endif /* USE_OMP */

    if (QI.level < Opts.subdiv_level0) {
	/* This is the 1st of the 2 entries for subdividing the element */
	subdivide_triangle(res);
	return;
    }

    for (i = 0; i < 3; i++)
	if (QI.rules[i] != NULL)
	    rule_n[i] = QI.rules[i]->n;	/* save position */
    if ((i = setjmp(QI.jmp_env)) != 0) {
	/* This is the 2nd of the 2 entries for subdividing the element */
	if (!Opts.dbg_off && Opts.show_recursions == TRUE
			  && QI.level >= Opts.subdiv_level0) {
	    fprintf(stderr, "Recursion from line %d, level %d\n", i, QI.level);
	    if (Opts.dbg_vtk)
		output_vtk(i);
	}
	for (i = 0; i < 3; i++)
	    if (QI.rules[i] != NULL)
		QI.rules[i]->n = rule_n[i];	/* save position */
	subdivide_triangle(res);
	return;
    }

    if (QI.level > Opts.subdiv_limit + Opts.subdiv_level0)
	phgError(1, "%s: recursion limit %d exceeded, abort.\n",
			__func__, (int)Opts.subdiv_limit);

    /* update 'QI' to match the current triangle */
    update_QI();

    /* Construct integration intervals for s */

    /* Note: cos(a) > Opts.threshold <==> sin(a)^2 < 1 - Opts.threshold^2 */
    thres = 1.0 - Opts.threshold * Opts.threshold;

    /* First, add projections of the intersections on the 3 edges */
    n = 0;
    for (i = 0; i < 3; i++) {	/* loop on the 3 egdes */
	int v0 = i, v1 = (i == 2 ? 0 : i + 1);
	j = intersect_edge(QI.tri[v0], QI.tri[v1], FALSE, &roots, &mults);
	assert(n + j <= INTER_MAX);
	for (k = 0; k < j; k++) {
	    FLOAT x, y, pt[Dim];
	    if (mults[k] > 1) {
		phgFree(roots);
		phgFree(mults);
		LONGJMP;
	    }
	    x = QI.tri[v0][0] + roots[k] * (QI.tri[v1][0] - QI.tri[v0][0]);
	    y = QI.tri[v0][1] + roots[k] * (QI.tri[v1][1] - QI.tri[v0][1]);
	    QI.ls_grad(x, y, pt);
	    /* check angle between tangent vector and r direction */
	    s = QI.transform[0] * pt[0] + QI.transform[1] * pt[1];
	    /* Note: s*s/|pt|^2 = sin(theta)^2 */
	    if (s * s < thres * (pt[0] * pt[0] + pt[1] * pt[1])) {
		phgFree(roots);
		phgFree(mults);
		LONGJMP; /* integr diretion too parallel to the interface */
	    }
	    xy2rs(x, y, NULL, &s, QI.transform);
	    sinter[n++] = s;
	}
	phgFree(roots);
	phgFree(mults);
    }

    /* Next, add projections of the 3 vertices */
    assert(n + 3 <= INTER_MAX);
    for (i = 0; i < 3; i++) {
	xy2rs(QI.tri[i][0], QI.tri[i][1], NULL, &s, QI.transform);
	sinter[n++] = s;
    }

    sinter_n = n;

    /* sort sinter[] */
    qsort(sinter, n, sizeof(sinter[0]), phgCompFLOAT);

    /* remove duplicate entries in sinter[] */
    for (i = 0, j = 1; j < sinter_n; ) {
	while (j < sinter_n && Fabs(sinter[j] - sinter[i]) < Opts.EPS)
	    j++;
	if (j >= sinter_n) 
	    break;
	sinter[++i] = sinter[j];
    }
    sinter_n = ++i;

    if (QI.rules[OM] != NULL || QI.rules[OP] != NULL) {
	/* integration in Omega- and Omega+ */
	FLOAT res0[3];
	if (res != NULL)
	    res[OM] = res[OP] = 0.;
	for (i = 0; i < sinter_n - 1; i = j) {
	    for (j = i + 1; j < sinter_n && sinter[j] - sinter[i] <= Opts.EPS;
									j++);
	    if (j >= sinter_n)
		break;
	    phgQuad1D(func_g, 3, sinter[i], sinter[j],
			  QI.quad, res == NULL ? NULL : res0, NULL);
	    if (!Opts.dbg_off
		&& res != NULL && Opts.show_intervals_s && !Opts.dbg_s)
		printf("s_interval (%0.16g,%0.16g), res0=%0.16lg\n",
			(double)sinter[i], (double)sinter[j], (double)res0[0]);
	    if (!Opts.dbg_off && res != NULL && Opts.dbg_s) {
		static FLOAT error0 = 0.0;
#if USE_OMP
# pragma omp threadprivate(error0)
#endif  /* USE_OMP */
		FLOAT error, res1[3];
		QI.skip_rule = TRUE;
		phgQuad1D(func_g, 3, sinter[i], sinter[j],
			  phgQuadGetQuad1D(QI.quad_order + 2), res1, NULL);
		QI.skip_rule = FALSE;
		error = Fabs((res0[0] - res1[0]) /
				(res1[0] == 0. ? 1. : res1[0]));
		if (error > error0 || Opts.show_intervals_s) {
		    printf("s_interval (%0.16g,%0.16g), res0=%0.16g, err=%e\n",
				(double)sinter[i], (double)sinter[j],
				(double)res0[0], (double)error);
		    if (error > error0)
			error0 = error;
		}
	    }
	    if (res != NULL)
		for (k = 0; k < 3; k++)
		    res[k] += res0[k];
	}
    }

    if (QI.rules[GM] != NULL) {
	/* integration on Gamma */
	FLOAT res0;
	if (res != NULL)
	    res[GM] = 0.;
	for (i = 0; i < sinter_n - 1; i = j) {
	    for (j = i + 1; j < sinter_n && sinter[j] - sinter[i] <= Opts.EPS;
									j++);
	    if (j >= sinter_n)
		break;
	    phgQuad1D(func_g_curve, 1, sinter[i], sinter[j],
			  QI.quad, res == NULL ? NULL : &res0, NULL);
	    if (res != NULL)
		res[GM] += res0;
	}
    }

    if (!Opts.dbg_off && res != NULL && Opts.dbg_elem && res[0] != 0.0) {
	static double error = 0.0;
	FLOAT err, res0[3];
	QUAD *quad_save = QI.quad;
	QI.quad = phgQuadGetQuad1D(QI.quad_order + 2);
	Opts.dbg_elem = FALSE;
	QI.skip_rule = TRUE;
	quad_interface_triangle(res0);
	QI.skip_rule = FALSE;
	Opts.dbg_elem = TRUE;
	err = Fabs(1.0 - res0[0] / res[0]);
	if (error < err) {
	    error = err;
	    fprintf(stderr, "Elem={{%0.16lg,%0.16lg},{%0.16lg,%0.16lg},"
		    "{%0.16lg,%0.16lg}}, res=%lg, res0=%lg, error=%le\n",
		    (double)QI.tri[0][0], (double)QI.tri[0][1],
		    (double)QI.tri[1][0], (double)QI.tri[1][1],
		    (double)QI.tri[2][0], (double)QI.tri[2][1],
		    (double)res[0], (double)res0[0], error);
	}
	QI.quad = quad_save;
    }
}

/* FIXME: potential bug with long double on x86_64 (with gcc, not icc):
 *   ./configure --with-float="long double"
 *  (cd src; rm -f *triangle.o; make USER_CFLAGS=-DLONG_DOUBLE_FIX=0 -s lib -j4)
 *  (cd examples; make interface; ./interface -xfem_ctol=0)
 *  *** ERROR: floating point exception: invalid operation.
 * Tested with gcc version 12.1.1 20220507 (Red Hat 12.1.1-1) (GCC) 
 */
#ifndef LONG_DOUBLE_FIX
# define LONG_DOUBLE_FIX	1	/* !0 => workaround for the bug */
#endif

void
phgQuadInterfaceTriangle(FUNC_2D ls, int ls_order, FUNC_2D ls_grad,
			 FLOAT const triangle[3][2], int quad_order,
			 FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
/* See comments in the function phgQuadInterface (in quad-interface.c),
 * or doc/quad-XFEM.doc for documentations about the arguments. */
{
    int i, j;
    FLOAT **rules[] = {rule_m, rule_0, rule_p};
#if !LONG_DOUBLE_FIX
# warning this produces invalid FP operations with long double
    FLOAT res[3] = {0., 0., 0.};
#else	/* !LONG_DOUBLE_FIX */
    FLOAT *res = NULL;
#endif	/* !LONG_DOUBLE_FIX */

    if (rule_m == NULL && rule_0 == NULL && rule_p == NULL)
	return;

    QI.ls = ls;
    QI.ls_order = ls_order;
    QI.ls_grad = ls_grad;

    for (i = 0; i < Dim + 1; i++) {
	for (j = 0; j < Dim; j++)
	    QI.tri[i][j] = triangle[i][j];
    }
    QI.quad = phgQuadGetQuad1D(quad_order);
    QI.quad_order = quad_order;
    QI.quad_type = 99;	/* to facilitate debugging */
    QI.level = 0;

    for (i = 0; i < 3; i++) {
	int np = (quad_order - 1) / 2 + 1; /* # points of 1D Gaussian rule */
	if (rules[i] == NULL) {
	    QI.rules[i] = NULL;	/* ignore rule i */
	    continue;
	}
	QI.rules[i] = phgCalloc(1, sizeof(*QI.rules[i]));
	QI.rules[i]->n = QI.rules[i]->alloc = 0;
	QI.rules[i]->trunc = 10 * np * (i == GM ? 1 : 2 * np);
	QI.rules[i]->pw = phgCalloc(RULE_HEADER, sizeof(FLOAT));
	QI.rules[i]->nv = NULL;
    }

    QI.skip_rule = FALSE;
    if (ls == NULL) {
	for (i = 0; i < 3; i += 2)
	    if (QI.rules[i] != NULL)
		break;
	if (i < 3) {
	    QI.quad_type = i - 1;
	    quad_triangle(res == NULL ? NULL : &res[i]);
	}
    }
    else {
	quad_interface_triangle(res);
    }

    for (i = 0; i < 3; i++) {
	if (QI.rules[i] == NULL)
	    continue;
	assert(QI.rules[i]->pw != NULL);
#if 0
	/* Set rule = NULL if (n == 0).
	 * Note: this may lead redundant calls since user programs may 
	 * decide to either regenerate the rule based on (rule == NULL) */
	if (QI.rules[i]->n == 0) {
	    phgFree(QI.rules[i]->pw);
	    phgFree(QI.rules[i]->nv);
	    phgFree(QI.rules[i]);
	    *rules[i] = NULL;
	    continue;
	}
#endif
	*rules[i] = phgQuadInterfaceRuleCreate(Dim, QI.rules[i]->n,
			QI.rules[i]->pw + RULE_HEADER, Dim + 1,	/* pts */
			QI.rules[i]->pw + RULE_HEADER + Dim, Dim + 1, /* wgts */
			1.0, QI.rules[i]->nv, Dim /* nv */);
	phgFree(QI.rules[i]->pw);
	phgFree(QI.rules[i]->nv);
	phgFree(QI.rules[i]);
    }
}

#ifndef TEST
# define TEST 0
#endif	/* !defined(TEST) */
#if TEST
int
main(int argc, char *argv[])
{
    FLOAT tri[][2] = {{0.,0.}, {1.,0.}, {0.,1.}}, t[3];
#if 1
    /* Case 1 */
    FLOAT pt1[] = {3.586020369539256e-14, 1.000000000000036};
    FLOAT pt2[] = {0.9999999999997176, -2.824962486158711e-13};
#else
    /* Case 2 */
    FLOAT pt1[] = {-1.448036540806387e-13, 1};
    FLOAT pt2[] = {2.255942614377446e-12, -2.400746268446905e-12};
#endif
    int i, n;

    phgInit(&argc, &argv);

    memcpy(QI.tri, tri, sizeof(tri));
    n = intersect_seg_triangle(pt1, pt2, t);
    printf("Number of cuts: %d\n", n);    
    for (i = 0; i < n; i++)
	printf("\tt[%d] = %0.16lg\n", i, t[i]);

    phgFinalize();
}
#endif	/* TEST */
