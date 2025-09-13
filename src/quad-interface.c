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

/* $Id: quad-interface.c,v 1.617 2022/09/12 09:39:53 zlb Exp $
 *
 * Numerical quadrature on an element containing a smooth interface */

#include "phg.h"
#include "phg/sturm.h"
#include "phg/elem-info.h"	/* ELEM_TYPE */

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <setjmp.h>

#define PI	_F(3.141592653589793238462643383279503)
#define ONE	_F(1.0)

static const char *adapt_args[] = {"0", "1", "2", NULL}; 
static const char *subdiv_type_args[] = {"bisection", "regular", NULL}; 
static const char *newton_args[] = {"auto", "yes", "no", NULL};

/*================= variables for cmdline options ===============*/
QI_OPTIONS _phg_quad_interface_options = {
    FLOAT_EPSILON * 1000.0,	/* EPS */
    0.95,	/* threshold */
    TRUE,	/* nv_flag */
    FALSE,	/* use_face_normal */
    FALSE,	/* shortcut_t */
    TRUE,	/* shortcut_s */
    0,		/* adapt */
    0,		/* subdiv_type */
    128,	/* subdiv_limit */
    0,		/* subdiv_level0 */
    0,		/* newton */
    10,		/* newton_maxits */
    0,		/* newton_porder */
    FALSE,	/* rect2tri */
    FALSE,	/* cuboid2tetra */
    1,		/* quad_n1 */
    1,		/* quad_n2 */
    0,		/* graph_t */
    0,		/* graph_s */
    0,		/* graph_w */
    -2.0,	/* dbg_w */
    8,		/* dbg_wN */
    -1.0,	/* dbg_value_t */
    -1.0,	/* dbg_value_s */
    FALSE,	/* dbg_roots */
    FALSE,	/* dbg_t */
    FALSE,	/* dbg_s */
    FALSE,	/* dbg_elem */
    FALSE,	/* dbg_vtk */
    FALSE,	/* show_recursions */
    FALSE,	/* show_directions */
    FALSE,	/* show_intervals_w */
    FALSE,	/* show_intervals_s */
    FALSE,	/* show_intervals_t */
    FALSE	/* dbg_off */
};
#define Opts	_phg_quad_interface_options

static QUAD *
get_quad_1D(int order)
{
    static QUAD quad, *q0 = NULL;
    static int n0 = -1, order0 = -1;
#if USE_OMP
# pragma omp threadprivate(quad, q0, n0, order0)
#endif  /* USE_OMP */
    QUAD *q = &quad;
    FLOAT *w, *p;
    int i, j;

    if (q0 == NULL) {
	q->weights = NULL;
	q->points = NULL;
    }

    if (q0 == NULL || n0 != Opts.quad_n1 || order0 != order) {
	phgFree(q->weights);
	phgFree(q->points);
	n0 = Opts.quad_n1;
	order0 = order;
	q0 = phgQuadGetQuad1D(order);
	if (n0 <= 1)
	    return q0;
	q->npoints = q0->npoints * n0;
	q->weights = w = phgAlloc(q->npoints * sizeof(*w));
	q->points = p = phgAlloc(2 * q->npoints * sizeof(*p));
	for (i = 0; i < n0; i++) {
	    for (j = 0; j < q0->npoints; j++, w++, p += 2) {
		*w = q0->weights[j] / n0;
		p[0] = (q0->points[2 * j] + i) / n0;
		p[1] = 1.0 - p[0];
	    }
	}
	q0 = q;
    }

    return q0;
}

static QUAD *
get_quad_2D(int order)
{
    static QUAD quad, *q0 = NULL;
    static int n0 = -1, order0 = -1;
#if USE_OMP
# pragma omp threadprivate(quad, q0, n0, order0)
#endif  /* USE_OMP */
    QUAD *q = &quad;
    FLOAT *w, *p, *p0, d;
    int i, j, j0, k, N;

    if (q0 == NULL) {
	q->weights = NULL;
	q->points = NULL;
    }

    if (q0 == NULL || n0 != Opts.quad_n2 || order0 != order) {
	phgFree(q->weights);
	phgFree(q->points);
	n0 = Opts.quad_n2;
	order0 = order;
	q0 = phgQuadGetQuad2D(order);
	if (n0 <= 1)
	    return q0;
	N = n0 * n0;
	q->npoints = q0->npoints * N;
	q->weights = w = phgAlloc(q->npoints * sizeof(*w));
	q->points = p = phgAlloc(3 * q->npoints * sizeof(*p));
	d = 1.0 / n0;
	for (i = 0; i < n0; i++) {
	    for (j = 0; j < n0 - i; j++) {
		for (j0 = 0; j0 < 2; j0++) {
		    FLOAT v0[2], v1[2], v2[2];
		    if (j0 == 0) {
			v0[0] = i * d;		v0[1] = j * d;
			v1[0] = i * d;		v1[1] = (j + 1) * d;
			v2[0] = (i + 1) * d;	v2[1] = j * d;
		    }
		    else {
			if (i + j + 2 > n0)
			    continue;
			v0[0] = (i + 1) * d;	v0[1] = j * d;
			v1[0] = i * d;		v1[1] = (j + 1) * d;
			v2[0] = (i + 1) * d;	v2[1] = (j + 1) * d;
		    }
		    for (k = 0, p0 = q0->points; k < q0->npoints;
				k++, w++, p += 3, p0 += 3) {
			*w = q0->weights[k] / N;
			p[0] = p0[0] * v0[0] + p0[1] * v1[0] + p0[2] * v2[0];
			p[1] = p0[0] * v0[1] + p0[1] * v1[1] + p0[2] * v2[1];
			p[2] = 1.0 - p[0] - p[1];
		    }
		}
	    }
	}
	q0 = q;
    }

    return q0;
}

#if 0
#define phgQuadGetQuad1D	get_quad_1D
#define phgQuadGetQuad2D	get_quad_2D
#endif

/*================= end test code using a composite rule ===============*/

/*--------------------- generic 1D numerical quadrature ----------------------*/

void
phgQuad1D(QUAD_1D_FUNC func, int dim, FLOAT a, FLOAT b, QUAD *quad,
						FLOAT *res, void *ctx)
/* computes and returns the integral \int_a^b f(x) dx using a Gauss quadrature
 * rule of order >= 'order.
 *
 * This function can be called recursively to compute multiple integrals */
{
    int i, k;
    FLOAT d, f[dim], *p, *w;

    assert(quad != NULL);

    d = b - a;
    if (d == 0.0) {
	if (res != NULL)
	    memset(res, 0, sizeof(*res) * dim);
	return;
    }

    p = quad->points;
    w = quad->weights;

    if (res != NULL)
	for (k = 0; k < dim; k++)
	    res[k] = 0.;

    for (i = 0; i < quad->npoints; i++, w++, p += 2) {
	func(a + d * *p, f, *w * d, ctx);
	if (res != NULL)
	    for (k = 0; k < dim; k++)
		res[k] += f[k] * (*w);
    }

    if (res != NULL)
	for (k = 0; k < dim; k++)
	    res[k] *= d;
}

/*--------------------------- utility functions ------------------------------*/

int
phgUnivariatePolyRoots(int n, FLOAT a, FLOAT b, FLOAT tol,
					FLOAT *buffer, int **mults)
/* computes the roots in the interval [a,b] of a univariate polynomial of order
 * 'n' and returns the number of roots found, or -1 if all coefficients of the
 * polynomial are 0.
 *
 * 'buffer'   - Array of length n+1.
 *		On input buffer[i] is the value of the polynomial at the point
 *			a+i*(b-a)/n (i=0,...,n).
 *		On output buffer[k] is the value of the k-th root.
 *
 * 'mults'    -	If not NULL, then it will point to a dynamically allocated
 *		integer array on output (to be freed by the calling function).
 *		with (*mults)[k] equal to the multiplicity of the k-th root.
 *
 * 'tol'      - It controls the precision of the computed roots. If tol<0, then
 *		the roots are not actually computed and the function returns
 *		!0 if there're roots in [a,b] and 0 otherwise.
 */
{
    /* cache the invariant factorized matrix (TODO: cache multiple orders) */
    static BOOLEAN initialized = FALSE;
    static int p = -1;		/* polynomial order */
    static FLOAT *A = NULL;	/* Vandermonde matrix */
    static int *pvt = NULL;	/* pivots of the LU factorization */
#if USE_OMP
# pragma omp threadprivate(p, A, pvt, initialized)
#endif  /* USE_OMP */
    int i, j, k, *iwk, n_save;
    FLOAT t, tp, *rwk, *buffer_save;

    if (initialized != TRUE) {
	initialized = TRUE;
	FreeAtExitNoCheck(A);
	FreeAtExitNoCheck(pvt);
    }

    if (p != n) {
	p = n;
	phgFree(A);
	phgFree(pvt);
	A = phgAlloc((p+1) * (p+1) * sizeof(*A));
	for (i = 0; i <= p; i++) {
	    tp = 1.0;
	    t = i / (FLOAT)p;
	    for (j = 0; j <= p; j++) {
		A[i * (p + 1) + (p - j)] = tp;
		tp *= t;
	    }
	}
	pvt = phgAlloc((p + 1) * sizeof(*pvt));
	i = phgSolverDenseLU(p + 1, A, pvt);
	assert(i == TRUE);
    }

#if 0
# warning test-only code!
    fprintf(stderr, "------------------------------------------------------\n");
    fprintf(stderr, "p = %d, tol = %le, function values:\n", p, tol);
    for (i = 0; i <= p; i++)
	fprintf(stderr, "\t%0.16le\n", buffer[i]);
#endif

    phgSolverDenseSV(p + 1, A, pvt, 1, buffer);

#if 0
# warning debugging code!
    fprintf(stderr, "polynomial values:\n");
    for (i = 0; i <= p; i++) {
	t = i/(FLOAT)p;
	tp = buffer[0];
	for (j = 1; j <=p; j++)
	    tp = tp * t + buffer[j];
	fprintf(stderr, "\t%0.16le\n", tp);
    }
    //exit(0);
#endif

#if 0
# warning debugging code!
    /* print polynomial for debugging */
    fprintf(stderr, "the polynomial:\n");
    for (i = 0; i <= n; i++)
	if (n - i > 1)
	    fprintf(stderr, "%s%0.16lg*x**%d",
			buffer[i] >= 0. && i != 0 ? "+" : "", buffer[i], n - i);
	else if (n - i > 0)
	    fprintf(stderr, "%s%0.16lg*x",
			buffer[i] >= 0. && i != 0 ? "+" : "", buffer[i]);
	else
	    fprintf(stderr, "%s%0.16lg",
			buffer[i] >= 0. && i != 0 ? "+" : "", buffer[i]);
    fprintf(stderr, "\n");
#endif

    /* Find the largest absolute value of the coefficients */
    t = 0.0;
    for (i = 0; i <= p; i++)
	if (Fabs(buffer[i]) > t)
	    t = Fabs(buffer[i]);
    if (t == 0.0) {
	/* the level-set function is constantly zero along the line */
#if 0
	/* return the two end-points as solutions with multiplicity 2 */
	if (mults != NULL) {
	    *mults = phgAlloc(2 * sizeof(**mults));
	    (*mults)[0] = (*mults)[1] = 2;
	}
	buffer[0] = a;
	buffer[1] = b;
	return 2;
#else
	if (mults != NULL)
	    *mults = NULL;
	return -1;
#endif
    }

#if 0
    /* normalize the coefficients (not needed, will be done in phgQuadSturm */
    t = 1.0 / t;
    for (i = 0; i <= p; i++)
	buffer[i] *= t;
    t = 1.0;
#endif

    for (i = 0; i <= p; i++) {
	if (Fabs(buffer[i]) > t * Opts.EPS)
	    break;
    }

    assert(i <= p);
    n_save = (n -= i);
    rwk = phgAlloc(((n+1)*(n+2)+4)/2 * sizeof(*rwk));
    iwk = phgAlloc((n+1) * sizeof(*iwk));
    buffer_save = phgAlloc((n+1) * sizeof(*rwk));
    if (i > 0) {
	for (j = 0; i <= p; i++, j++)
	    buffer[j] = buffer[i];
    }
    memcpy(buffer_save, buffer, (n + 1) * sizeof(*buffer));

    n = phgQuadSturm(buffer, n, 0.0, 1.0, tol, Opts.EPS, rwk, iwk);

#if 0
    /* check the solutions */
    for (i = 0; i < n; i++) {
	FLOAT residual, x = buffer[i];
	residual = buffer_save[0];
	for (j = 1; j <= n_save; j++)
	    residual = x * residual + buffer_save[j];
	if (Fabs(residual) > 1e-6) {
		fprintf(stderr, "p = %d, coeff =", n_save);
		for (i = 0; i <= n_save; i++)
		    fprintf(stderr, " %0.16le", (double)buffer_save[i]);
		fprintf(stderr, "\nwrong solution: %0.16le, residual: %le\n",
			(double)x, (double)residual);
		fprintf(stderr, "abort.\n");
		exit(1);
	}
    }
#else
    Unused(n_save);
#endif
    phgFree(buffer_save);

    if (n == 0 || tol < 0.) {
	phgFree(rwk);
	phgFree(iwk);
	if (mults != NULL)
	    *mults = NULL;
	return n;
    }

    /* rescale the roots */
    if (a != 0.0 || b != 1.0) {
	for (i = 0; i < n; i++)
	    rwk[i] = a + buffer[i] * (b - a);
    }
    else {
	for (i = 0; i < n; i++)
	    rwk[i] = buffer[i];
    }

    /* merge very close roots */
    qsort(rwk, n, sizeof(*rwk), phgCompFLOAT);
    for (k = i = 0; i < n; i++, k++) {
	buffer[k] = rwk[i];
	iwk[k] = iwk[i];
	for (j = i + 1; j < n; j++) {
	    tp = Fabs(rwk[j] - rwk[i]);
	    if (tp * tp > Opts.EPS * 10.0)
		break;
	    iwk[k] += iwk[j];
	}
	i = j - 1;
    }
    n = k;

#if 0
    fprintf(stderr, "Number of roots: %d\n", n);
    for (i = 0; i < n; i++)
	fprintf(stderr, "root %d: multiplicity %d, t = %lg\n", i+1, iwk[i],
							(double)buffer[i]);
#endif

    phgFree(rwk);
    if (mults != NULL)
	*mults = iwk;
    else
	phgFree(iwk);

    return n;
}

typedef struct TET_ {
    FLOAT	volume;			/* volume (optional) */
    FLOAT	(*verts)[Dim];		/* coordinates of the vertices */
    BOOLEAN	is_bdry[NVert];		/* TRUE => no neighbour on the face */
    BOOLEAN	owns_verts;		/* whether it owns the verts[] array */
} TET;

typedef struct {
    const FLOAT	*verts[NVert], *u_st, *nu, *v_t, *nv, *x0, *nw;
    const TET	*pt;
    FLOAT	Jacobian[(Dim + 1) * (Dim + 1)];
    FLOAT	hu, hv, hw;
    const FLOAT	*normal[NFace];		/* unit normal vectors of the faces */
    struct EDGE_CUTS {
	int	n;		/* number of cuts on each edge */
	FLOAT	cuts[2];	/* list of cuts on each edge (at most 2) */
    }		edge_cuts[NEdge];	/* list of edge cuts */
    struct FACE_CUTS {
	int	n;		/* number of cuts on each face */
	int	info[4];	/* pointers to the edge cuts (i*NEdge+eno) */
	COORD	tv[4];		/* tangents of the interface curve */
    }		face_cuts[NFace];	/* list of face cuts */
    int		interv_n;	/* number of inadmissible intervals */
    FLOAT	interv[5 * 2];	/* list of inadmissible intervals */

    /* jmp_buf's for skipping empty integration intervals */
    jmp_buf	jmp_env_s;
    jmp_buf	jmp_env_t;
} QUAD_CONTEXT;

enum {OM = 0, GM = 1, OP = 2}; 
static struct {
    QUAD_CONTEXT	ctx;		/* context data */
    DOF			*ls, *ls_grad;
    ELEMENT		*e;
    FUNC_3D		func;
    FLOAT		H;		/* reference length */
    struct {FLOAT *pw, *nv; int n, alloc, trunc;} *rules[3];
    INT			index;		/* element index */
    int			ls_order;
    int			quad_type;
    int			quad_order;
    BOOLEAN		skip_rule;
} quad_info;
#define QI quad_info

static int depth = 0;

#if 0
/* Note: seems unnecessary to use a stack of jmp_buf's */
# define MAX_DEPTH	32
static jmp_buf		jmp_env[MAX_DEPTH + 1];
# define JMP_ENV		jmp_env[depth]
#else
# ifdef MAX_DEPTH
#  undef		MAX_DEPTH
# endif
static jmp_buf		jmp_env;
# define JMP_ENV		jmp_env
#endif
/* Note: subdiv_limit<=0 ==> no subdivision performed */
#define LONGJMP		if (Opts.subdiv_limit > 0) longjmp(JMP_ENV, __LINE__)

#if USE_OMP
# pragma omp threadprivate(jmp_env, depth, quad_info)
#endif  /* USE_OMP */

#define MAX_EDGE_CUTS \
	((int)(sizeof(qc->edge_cuts->cuts) / sizeof(qc->edge_cuts->cuts[0])))
#define MAX_FACE_CUTS \
	((int)(sizeof(qc->face_cuts->info) / sizeof(qc->face_cuts->info[0])))

#define output_vtk3(fn, v0, v1, v2, v3)	output_vtk3_(fn, v0,v1,v2,v3, __LINE__)
#define output_vtk2(fn, v0, v1, v2)	output_vtk2_(fn, v0, v1, v2, __LINE__)
#define output_vtk(fn, verts) \
    output_vtk3(fn, verts[0], verts[1], verts[2], verts[3]);

static const char *
output_vtk_fn(const char *fn, int line_no)
{
    static char s[32];
#if USE_OMP
# pragma omp threadprivate(s)
#endif  /* USE_OMP */
    char s1[sizeof(s)], *ext;

    assert(strlen(fn) + 10 < sizeof(s));

    ext = strrchr(fn, '.');
    if (ext == NULL) {
	ext = ".vtk";
    }
    else {
	strcpy(s1, fn);
	s1[ext - fn] = '\0';
	fn = s1;
    }

#if defined(__GNUC__) && (__GNUC__ > 7)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wformat-overflow"
#endif

    sprintf(s, "%s-%d%s", fn, line_no, ext);

#if defined(__GNUC__) && (__GNUC__ > 7)
# pragma GCC diagnostic pop
#endif

    return s;
}

static void
output_vtk3_(const char *fn, const FLOAT *v0, const FLOAT *v1, const FLOAT *v2,
	     const FLOAT *v3, int line_no)
{
    FILE *fp;

    if (Opts.dbg_off || !Opts.dbg_vtk)
	return;
    fp = fopen(fn = output_vtk_fn(fn, line_no), "wt");
    fprintf(fp, "# vtk DataFile Version 2.0\n"
		"vtk output\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n"
		"POINTS 4 float\n"
		"%0.16lg %0.16lg %0.16lg\n%0.16lg %0.16lg %0.16lg\n"
		"%0.16lg %0.16lg %0.16lg\n%0.16lg %0.16lg %0.16lg\n"
		"CELLS 1 5\n"
		"4 0 1 2 3\n"
		"CELL_TYPES 1\n"
		"10\n",
		(double)v0[0], (double)v0[1], (double)v0[2],
		(double)v1[0], (double)v1[1], (double)v1[2],
		(double)v2[0], (double)v2[1], (double)v2[2],
		(double)v3[0], (double)v3[1], (double)v3[2]);
    fclose(fp);
    fprintf(stderr, "*** VTK file \"%s\" created.\n", fn);
}

static void
output_vtk2_(const char *fn, const FLOAT *v0, const FLOAT *v1, const FLOAT *v2,
	     int line_no)
{
    FILE *fp;

    if (Opts.dbg_off || !Opts.dbg_vtk)
	return;
    fp = fopen(fn = output_vtk_fn(fn, line_no), "wt");
    fprintf(fp, "# vtk DataFile Version 2.0\n"
		"vtk output\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n"
		"POINTS 3 float\n"
		"%0.16lg %0.16lg %0.16lg\n%0.16lg %0.16lg %0.16lg\n"
		"%0.16lg %0.16lg %0.16lg\n"
		"CELLS 1 4\n"
		"3 0 1 2\n"
		"CELL_TYPES 1\n"
		"5\n",
		(double)v0[0], (double)v0[1], (double)v0[2],
		(double)v1[0], (double)v1[1], (double)v1[2],
		(double)v2[0], (double)v2[1], (double)v2[2]);
    fprintf(stderr, "*** VTK file \"%s\" created.\n", fn);
}

static int
adapt_order(const FLOAT h)
/* computes and returns the quadrature order according to interval size h */
{
    FLOAT fac, tol, d = h / QI.H;
    int i, k, q, p = QI.quad_order;

    if (Opts.adapt == 0)
	return QI.quad_order;

    fac = 1.0;
    for (i = 2; i <= p + 1; i++)
	fac *= i;

    if (Opts.adapt == 1) {
	/* assume H == 1: find q such that 1/(p+1)! >= d^(q+1)/(q+1)! */
	tol = 1.0 / fac;
	k = 1;
    }
    else {
	/* assume h <= H^2:  find q such that d^p/(p+1)! >= d^(2*q+1)/(q+1)! */
	tol = Pow(d, p) / fac;
	k = 2;
    }

    for (q = p - 1; q >= 1; q--) {
	fac /= q + 2;
	if (Pow(d, k * q + 1) / fac > tol)
	    break;
    }
    q++;

    return q;
}

static void
levelset_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if (QI.ls->type == DOF_ANALYTIC && QI.ls->userfunc != NULL) {
	QI.ls->userfunc(x, y, z, value);
    }
    else {
	FLOAT lam[4];
	assert(QI.ls->g != NULL && QI.e != NULL);
	phgGeomXYZ2Lambda(QI.ls->g, QI.e, x, y, z, lam);
	phgDofEval(QI.ls, QI.e, lam, value);
    }
}

static void
levelset_grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *grad, BOOLEAN normalize)
{
    if (QI.ls_grad != NULL && QI.ls_grad->type == DOF_ANALYTIC &&
	QI.ls_grad->userfunc != NULL) {
	QI.ls_grad->userfunc(x, y, z, grad);
    }
    else {
	FLOAT lam[4];
	phgGeomXYZ2Lambda(QI.ls->g, QI.e, x, y, z, lam);
	if (QI.ls_grad == NULL)
	    phgDofEvalGradient(QI.ls, QI.e, lam, NULL, grad);
	else
	    phgDofEval(QI.ls_grad, QI.e, lam, grad);
    }

    if (normalize == TRUE) {
	FLOAT d;
	d = Sqrt(grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]);
	assert(d > Opts.EPS);
	d = 1.0 / d;
	grad[0] *= d;
	grad[1] *= d;
	grad[2] *= d;
    }
}

static void
compute_Jacobian(void)
/* computes the Jabobian
 *	J = D\lambda / Dx = [ D\lambda_0/Dx, D\lambda_0/Dy, D\lambda_0/Dz, c0;
 *			      D\lambda_1/Dx, D\lambda_1/Dy, D\lambda_1/Dz, c1;
 *			      D\lambda_2/Dx, D\lambda_2/Dy, D\lambda_2/Dz, c2;
 *			      D\lambda_3/Dx, D\lambda_3/Dy, D\lambda_3/Dz, c3 ]
 * of the element, by inverting the matrix
 *
 *	[x0, x1, x2, x3; y0, y1, y2, y3; z0, z1, z2, z3; 1 1 1 1].
 *
 * Note: lambda = J[0..Dim][0..Dim-1] * x + J[0..Dim][Dim]
 */
{
    int i, j, k;
    int n = Dim + 1;
    FLOAT d, t, a[Dim + 1][2 * (Dim + 1)];
    FLOAT *J = QI.ctx.Jacobian;

    /* right hand side (I) */
    for (j = 0; j < n; j++) 
	for (i = 0; i < n; i++) 
	    a[j][n + i] = (i == j) ? 1.0 : 0.0;

    /* coefficient matrix */
    for (j = 0; j < n; j++) {
	for (i = 0; i < Dim; i++)
	    a[j][i] = QI.ctx.verts[j][i];
	a[j][Dim] = 1.0;
    }

    /* solve linear systems of equations */
    for (j = 0; j < n; j++) {
	/* find pivot */
	d = Fabs(a[j][j]);
	k = j;
	for (i = j + 1; i < n; i++)
	    if (d < (t = Fabs(a[i][j]))) {
		d = t;
		k = i;
	    }
	if (k != j) /* swap rows j and k */
	    for (i = j; i < n + n; i++) {
		t = a[j][i];
		a[j][i] = a[k][i];
		a[k][i] = t;
	    }
	d = 1.0 / ((t = a[j][j]) == 0. ? 1. : t);
	for (i = j + 1; i < n + n; i++)
	    a[j][i] *= d;
	for (i = j + 1; i < n; i++) {
	    d = a[i][j];
	    for (k = j + 1; k < n + n; k++)
		a[i][k] -= d * a[j][k];
	}
    }

    for (j = n - 2; j >= 0; j--)
	for (i = j + 1; i < n; i++) {
	    d = a[j][i];
	    for (k = n; k < n + n; k++)
		a[j][k] -= d * a[i][k];
	}

    /* Jacobian = trans(a[0..3, 4..7]) */
    for (j = 0; j < n; j++)
	for (i = 0; i < n; i++)
	    *(J++) = a[i][j + n];
}

static FLOAT
levelset_face(int k, FLOAT x, FLOAT y, FLOAT z)
/* returns the k-th bary-coordinate of the point (x,y,z) */
{
    const FLOAT *J = QI.ctx.Jacobian;

    assert(k >= 0 && k <= Dim);
    J += k * 4;
    return J[0] * x + J[1] * y + J[2] * z + J[3];
}

typedef struct {
    const FLOAT	*v0, *v1;
    COORD	dv;		/* = v1 - v0 */
} NS_CONTEXT;

static BOOLEAN
newton_func(const FLOAT *t, FLOAT *f, void *pctx)
{
    NS_CONTEXT *nc = pctx;

    levelset_func(nc->v0[0] + (*t) * nc->dv[0], nc->v0[1] + (*t) * nc->dv[1],
		  nc->v0[2] + (*t) * nc->dv[2], f);

    return TRUE;
}

static BOOLEAN
newton_jac(const FLOAT *t, FLOAT *jac, void *pctx)
{
    NS_CONTEXT *nc = pctx;
    FLOAT f[3];

    levelset_grad(nc->v0[0] + (*t) * nc->dv[0], nc->v0[1] + (*t) * nc->dv[1],
		  nc->v0[2] + (*t) * nc->dv[2], f, FALSE);
    *jac = f[0] * nc->dv[0] + f[1] * nc->dv[1] + f[2] * nc->dv[2];

    return TRUE;
}

static int
intersect_edge(const FLOAT *v0, const FLOAT *v1, FLOAT tol, BOOLEAN in_flag,
			FLOAT **roots, int **mults)
/* computes intersections of the edge 'v0'-'v1' with the interface given by the
 * level-set function 'QI.ls' and returns the number of intersection points.
 *
 * The pointers to the roots and multiplicities are optionally returned in
 * 'roots' and 'mults' respectively and the returned arrays are to be freed by
 * the calling function.
 *
 * If in_flag == TRUE then only the roots inside the tetrahedron are computed.
 *
 */
{
    int i, j, order, n;
    FLOAT t, t0, t1, lam[Dim + 1], lam0[Dim + 1], lam1[Dim + 1], *c;
    FLOAT a, b, fa, fb, f;
    BOOLEAN use_newton;
    NS_CONTEXT ctx;

    t0 = 0.0;
    t1 = 1.0;
#if 1
    /* find the intersection part of {v0+t*(v1-v0)} with the tetrahedron
     * if in_flag == TRUE */
    for (i = 0; in_flag != FALSE && i < NFace; i++) {
	a = levelset_face(i, v0[0], v0[1], v0[2]);
	b = levelset_face(i, v1[0], v1[1], v1[2]);
	/* The inequality:
	 *	a + t * (b - a) >= 0,
	 * implies:
	 *	t >= -a / (b -  a), if (b > a)
	 *	t <= -a / (b -  a), if (b < a)
	 */
	t = b - a;
	if (t * t < Opts.EPS) {
	    /* FIXME: floating-point instability when v0-v1 is (almost)
	     * parallel to the face. Test case (sphere or cylinder):
	     *		-order 35 -mesh_file ellipsoid.dat
	     */
	    if (a > -Opts.EPS || a * a < Opts.EPS)
		continue;
	    /* empty set */
	    t0 = t1 = 0.0;
	    break;
	}
	else if (t > 0.0) {
	    t = -a / t;
	    if (t0 < t)
		t0 = t;
	}
	else {
	    t = -a / t;
	    if (t1 > t)
		t1 = t;
	}
	if (t0 >= t1 + Opts.EPS)
	    break;
    }

    if (t0 >= t1 - FLOAT_EPSILON * 100.0) {
	*roots = NULL;
	*mults = NULL;
	return 0;
    }

    if ((t0 -= Opts.EPS) < 0.0)
	t0 = 0.0;
    if ((t1 += Opts.EPS) > 1.0)
	t1 = 1.0;
#endif

    order = QI.ls_order;
    use_newton = (Opts.newton == 1 || (Opts.newton == 0 && order < 0));
    order = Opts.newton_porder > 0 ? Opts.newton_porder : 
				     (order < 0 ? -order : order);

    /* construct QI.ls(v0+t*(v1-v0)) (an univariate polynomial in t) */
    *roots = c = phgAlloc((order + 1) * sizeof(*c));
    if (QI.ls->type != DOF_ANALYTIC || QI.ls->userfunc == NULL) {
	/* This saves some (x,y,z) to lambda conversions */
	assert(QI.ls->g != NULL && QI.e != NULL);
	phgGeomXYZ2Lambda(QI.ls->g, QI.e, v0[0], v0[1], v0[2], lam0);
	phgGeomXYZ2Lambda(QI.ls->g, QI.e, v1[0], v1[1], v1[2], lam1);
	for (i = 0; i <= order; i++) {
	    t = t0 + i / (FLOAT)order * (t1 - t0);
	    for (j = 0; j <= Dim; j++)
		lam[j] = lam0[j] + t * (lam1[j] - lam0[j]);
	    phgDofEval(QI.ls, QI.e, lam, &c[i]);
	}
    }
    else {
	for (i = 0; i <= order; i++) {
	    t = t0 + i / (FLOAT)order * (t1 - t0);
	    levelset_func(v0[0] + t * (v1[0] - v0[0]),
			  v0[1] + t * (v1[1] - v0[1]),
			  v0[2] + t * (v1[2] - v0[2]), &c[i]);
	}
    }
    n = phgUnivariatePolyRoots(order, t0, t1, tol, c, mults);

    if (n == 0)
	return n;

#if 0
#warning remove me!
	if (random() % 1000 == 0) LONGJMP;
#endif

#if 0
    if (n > 2)
	LONGJMP;
#endif

    /* perform Newton iterations */
    ctx.v0 = v0;
    ctx.v1 = v1;
    ctx.dv[0] = v1[0] - v0[0];
    ctx.dv[1] = v1[1] - v0[1];
    ctx.dv[2] = v1[2] - v0[2];
    j = -1;
    for (i = 0; i < n; i++) {
	int its;
	FLOAT res;

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
			       Opts.newton_maxits, Opts.EPS, Opts.EPS, &res,
			       &ctx);
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
	    /*phgWarning("FIXME: Newton iterations failed.\n");*/
#if 0
# warning debugging code!
	fprintf(stderr, "  root %d/%d: %lg, diff: %le, its: %d\n",
			i, n, c[i], t - c[i], its);
	fprintf(stderr, "  v0=vector([%lg,%lg,%lg])\n"
			"  v1=vector([%lg,%lg,%lg])\n",
			v0[0], v0[1], v0[2], v1[0], v1[1], v1[2]);
	output_vtk("tmp.vtk", QI.ctx.verts);
	/* Sample SageMath code to print the function L(v0+x*(v1-v0)):
		var('t')
		[x,y,z] = v0+t*(v1-v0)
		g(t) = cos(5*x-2.5)*sin(5*y-2.5)+cos(5*y-2.5)*sin(5*z-2.5)+cos(5*z-2.5)*sin(5*x-2.5)
		var('x')
		g(x)
	 */
#endif
	    /* try bisection */
	    a = (j < 0      ? t0 : (c[i] + c[j]) * 0.5);
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
	    else while (TRUE) {
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
#if 0
	    newton_func(&t, &f, &ctx);
	    newton_func(&c[i], &fa, &ctx);
	    printf("before bisection: f(%lg)=%le, after bisection: f(%lg)=%le, diff: %le\n", c[i], fa, t, f, c[i] - t);
#endif
	    j++;
	    c[j] = t;
	    (*mults)[j] = 1;
	}
    }
    n = j + 1;

    return n;
}

/*---------- function for building a tree of subdivided tetrahedra ----------*/

TET *
tet_new(BOOLEAN alloc_verts)
{
    int i;
    TET *pt = phgCalloc(1, sizeof(*pt));

    pt->owns_verts = alloc_verts;
    pt->verts = (alloc_verts == FALSE ? NULL :
					phgAlloc(NVert * sizeof(*pt->verts)));
    pt->volume = -1.0;
    for (i = 0; i < NVert; i++)
	pt->is_bdry[i] = FALSE;		/* TODO */
    return pt;
}

static void
tet_free(TET *pt)
{
    if (pt == NULL)
	return;

    if (pt->owns_verts != FALSE)
	phgFree(pt->verts);
    phgFree(pt);
}

static FLOAT
get_volume(const FLOAT *v0, const FLOAT *v1, const FLOAT *v2, const FLOAT *v3)
{
    FLOAT x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3;

    x = v0[0]; 
    y = v0[1];
    z = v0[2];
    x1 = v1[0] - x;
    y1 = v1[1] - y;
    z1 = v1[2] - z;
    x2 = v2[0] - x;
    y2 = v2[1] - y;
    z2 = v2[2] - z;
    x3 = v3[0] - x;
    y3 = v3[1] - y;
    z3 = v3[2] - z;
    return Fabs(x1*y2*z3 + x2*y3*z1 + y1*z2*x3 -
		z1*y2*x3 - y1*x2*z3 - z2*y3*x1) * (1.0 / (FLOAT)6.0);
}

static FLOAT
tet_volume(TET *pt)
{
    if (pt->volume < 0.0) {
	pt->volume = get_volume(pt->verts[0], pt->verts[1], pt->verts[2],
				pt->verts[3]);
#if 0
	if (pt->volume == 0.0)
	    phgWarning("volume=0!\n");
#endif
    }
    return pt->volume;
}

#define	nChild		(Opts.subdiv_type == 0 ? 2 : 8)

static void
tet_bisect(TET *pt, TET *ppt[])
/* subdivide the tetra 'pt' by longest edge bisection */
{
    int i, k, i0, i1;
    FLOAT a, b, *v0, *v1;

    /* find the longest edge */
    a = -1.0;
    k = -1;
    for (i = 0; i < NEdge; i++) {
	v0 = pt->verts[GetEdgeVertex(i, 0)];
	v1 = pt->verts[GetEdgeVertex(i, 1)];
	b = Pow(v1[0]-v0[0], 2) + Pow(v1[1]-v0[1], 2) + Pow(v1[2]-v0[2], 2);
	if (a < b) {
	    a = b;
	    k = i;
	}
    }

    ppt[0] = tet_new(TRUE);
    ppt[1] = tet_new(TRUE);

    i0 = GetEdgeVertex(k, 0);
    i1 = GetEdgeVertex(k, 1);
#if 0
    i2 = GetEdgeVertex(5 - k, 0);
    i3 = GetEdgeVertex(5 - k, 1);
#endif

    v0 = pt->verts[i0];
    v1 = pt->verts[i1];

    /* t0: v0 v v2 v3 */
    for (i = 0; i < NVert; i++) {
	if (i != i1) {
	    memcpy(ppt[0]->verts[i], pt->verts[i], sizeof(COORD));
	    if (i != i0)
		ppt[0]->is_bdry[i] = pt->is_bdry[i];
	    else
		ppt[0]->is_bdry[i] = FALSE;	/* new face */
	}
	else {
	    ppt[0]->verts[i][0] = (v0[0] + v1[0]) * 0.5;
	    ppt[0]->verts[i][1] = (v0[1] + v1[1]) * 0.5;
	    ppt[0]->verts[i][2] = (v0[2] + v1[2]) * 0.5;
	    ppt[0]->is_bdry[i] = pt->is_bdry[i1];
	}
    }
    ppt[0]->volume = (pt->volume < 0.0 ? -1.0 : pt->volume * 0.5);

    /* t1: v v1 v2 v3 */
    for (i = 0; i < NVert; i++) {
	if (i != i0) {
	    memcpy(ppt[1]->verts[i], pt->verts[i], sizeof(COORD));
	    if (i != i1)
		ppt[1]->is_bdry[i] = pt->is_bdry[i];
	    else
		ppt[1]->is_bdry[i] = FALSE;	/* new face */
	}
	else {
	    memcpy(ppt[1]->verts[i], ppt[0]->verts[i1], sizeof(COORD));
	    ppt[1]->is_bdry[i] = pt->is_bdry[i0];
	}
    }
    ppt[1]->volume = ppt[0]->volume;
}

static void
tet_regular_refine(TET *pt, TET *ppt[])
/* subdivide the tetra 'pt' by regular refinement */
{
    int i, i0, i1;
    COORD c[NEdge];
    TET *pt_new;
    FLOAT volume;

    volume = (pt->volume < 0.0 ? -1.0 : pt->volume * 0.125);
    for (i = 0; i < 8; i++) {
	ppt[i] = tet_new(TRUE);
	ppt[i]->volume = volume;
    }

    for (i = 0; i < NEdge; i++) {
	i0 = GetEdgeVertex(i, 0);
	i1 = GetEdgeVertex(i, 1);
	c[i][0] = (pt->verts[i0][0] + pt->verts[i1][0]) * 0.5;
	c[i][1] = (pt->verts[i0][1] + pt->verts[i1][1]) * 0.5;
	c[i][2] = (pt->verts[i0][2] + pt->verts[i1][2]) * 0.5;
    }

    /* The ordering of vertices below ensures that at most 3 congruent types
     * of tetrahedron will be generated by successive refinements, see, e.g.,
     *
     *	H. Freudentahl. Simplizialzerlegung von Beschr\"ankter Flachheit.
     *	Annals of Mathematics, 3:580-582, 1942.
     *
     * 	J. Bey, Tetrahedral Grid Refinement, Computing, 55:355-378, 1995. */


    /* The four corners */
    pt_new = ppt[0];
    memcpy(pt_new->verts[0], pt->verts[0], sizeof(COORD));
    pt_new->is_bdry[0] = FALSE;	/* inner face */
    memcpy(pt_new->verts[1], c[GetEdgeNo(0,1)], sizeof(COORD));
    pt_new->is_bdry[1] = pt->is_bdry[1];
    memcpy(pt_new->verts[2], c[GetEdgeNo(0,2)], sizeof(COORD));
    pt_new->is_bdry[2] = pt->is_bdry[2];
    memcpy(pt_new->verts[3], c[GetEdgeNo(0,3)], sizeof(COORD));
    pt_new->is_bdry[3] = pt->is_bdry[3];

    pt_new = ppt[1];
    memcpy(pt_new->verts[0], c[GetEdgeNo(1,0)], sizeof(COORD));
    pt_new->is_bdry[0] = pt->is_bdry[0];
    memcpy(pt_new->verts[1], pt->verts[1], sizeof(COORD));
    pt_new->is_bdry[1] = FALSE;	/* inner face */
    memcpy(pt_new->verts[2], c[GetEdgeNo(1,2)], sizeof(COORD));
    pt_new->is_bdry[2] = pt->is_bdry[2];
    memcpy(pt_new->verts[3], c[GetEdgeNo(1,3)], sizeof(COORD));
    pt_new->is_bdry[3] = pt->is_bdry[3];

    pt_new = ppt[2];
    memcpy(pt_new->verts[0], c[GetEdgeNo(2,0)], sizeof(COORD));
    pt_new->is_bdry[0] = pt->is_bdry[0];
    memcpy(pt_new->verts[1], c[GetEdgeNo(2,1)], sizeof(COORD));
    pt_new->is_bdry[1] = pt->is_bdry[1];
    memcpy(pt_new->verts[2], pt->verts[2], sizeof(COORD));
    pt_new->is_bdry[2] = FALSE;	/* inner face */
    memcpy(pt_new->verts[3], c[GetEdgeNo(2,3)], sizeof(COORD));
    pt_new->is_bdry[3] = pt->is_bdry[3];

    pt_new = ppt[3];
    memcpy(pt_new->verts[0], c[GetEdgeNo(3,0)], sizeof(COORD));
    pt_new->is_bdry[0] = pt->is_bdry[0];
    memcpy(pt_new->verts[1], c[GetEdgeNo(3,1)], sizeof(COORD));
    pt_new->is_bdry[1] = pt->is_bdry[1];
    memcpy(pt_new->verts[2], c[GetEdgeNo(3,2)], sizeof(COORD));
    pt_new->is_bdry[2] = pt->is_bdry[2];
    memcpy(pt_new->verts[3], pt->verts[3], sizeof(COORD));
    pt_new->is_bdry[3] = FALSE;	/* inner face */

    /* The four interior ones (the octohedron) */
    pt_new = ppt[4];
    memcpy(pt_new->verts[0], c[GetEdgeNo(0,1)], sizeof(COORD));
    pt_new->is_bdry[0] = pt->is_bdry[1];
    memcpy(pt_new->verts[1], c[GetEdgeNo(0,2)], sizeof(COORD));
    pt_new->is_bdry[1] = pt->is_bdry[2];
    memcpy(pt_new->verts[2], c[GetEdgeNo(0,3)], sizeof(COORD));
    pt_new->is_bdry[2] = pt->is_bdry[3];
    memcpy(pt_new->verts[3], c[GetEdgeNo(1,3)], sizeof(COORD));
    pt_new->is_bdry[3] = FALSE;	/* inner face */

    pt_new = ppt[5];
    memcpy(pt_new->verts[0], c[GetEdgeNo(0,1)], sizeof(COORD));
    pt_new->is_bdry[0] = pt->is_bdry[0];
    memcpy(pt_new->verts[1], c[GetEdgeNo(0,2)], sizeof(COORD));
    pt_new->is_bdry[1] = FALSE;	/* inner face */
    memcpy(pt_new->verts[2], c[GetEdgeNo(1,2)], sizeof(COORD));
    pt_new->is_bdry[2] = pt->is_bdry[2];
    memcpy(pt_new->verts[3], c[GetEdgeNo(1,3)], sizeof(COORD));
    pt_new->is_bdry[3] = pt->is_bdry[3];

    pt_new = ppt[6];
    memcpy(pt_new->verts[0], c[GetEdgeNo(0,2)], sizeof(COORD));
    pt_new->is_bdry[0] = FALSE;	/* inner face */
    memcpy(pt_new->verts[1], c[GetEdgeNo(0,3)], sizeof(COORD));
    pt_new->is_bdry[1] = pt->is_bdry[0];
    memcpy(pt_new->verts[2], c[GetEdgeNo(1,3)], sizeof(COORD));
    pt_new->is_bdry[2] = pt->is_bdry[1];
    memcpy(pt_new->verts[3], c[GetEdgeNo(2,3)], sizeof(COORD));
    pt_new->is_bdry[3] = pt->is_bdry[2];

    pt_new = ppt[7];
    memcpy(pt_new->verts[0], c[GetEdgeNo(0,2)], sizeof(COORD));
    pt_new->is_bdry[0] = pt->is_bdry[0];
    memcpy(pt_new->verts[1], c[GetEdgeNo(1,2)], sizeof(COORD));
    pt_new->is_bdry[1] = pt->is_bdry[1];
    memcpy(pt_new->verts[2], c[GetEdgeNo(1,3)], sizeof(COORD));
    pt_new->is_bdry[2] = FALSE;	/* inner face */
    memcpy(pt_new->verts[3], c[GetEdgeNo(2,3)], sizeof(COORD));
    pt_new->is_bdry[3] = pt->is_bdry[3];
}

static void
tet_subdivide(TET *pt, TET *ppt[])
/* subdivide the tetra 'pt' by regular refinement or bisection */
{
    (Opts.subdiv_type == 0 ? tet_bisect : tet_regular_refine)(pt, ppt);
}

static void
get_face_normal(FLOAT n[], const FLOAT *v0, const FLOAT *v1, const FLOAT *v2,
		const FLOAT *v3)
/* returns in n[] the unit normal vector of the triangle v0, v1, v2.
 * if v3 is not NULL, n is set to the opposite direction of v3-v0 */
{
    FLOAT d;

    n[0] = (v1[1]-v0[1]) * (v2[2]-v0[2]) - (v1[2]-v0[2]) * (v2[1]-v0[1]);
    n[1] = (v1[2]-v0[2]) * (v2[0]-v0[0]) - (v1[0]-v0[0]) * (v2[2]-v0[2]);
    n[2] = (v1[0]-v0[0]) * (v2[1]-v0[1]) - (v1[1]-v0[1]) * (v2[0]-v0[0]);
    d = Sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
#if 0
    /* check for degenerated tetra when using bisection */
    assert(d > Opts.EPS /** Opts.EPS*/);
#endif
    d = 1.0 / d;
    if (v3 != NULL &&
	n[0]*(v3[0]-v0[0]) + n[1]*(v3[1]-v0[1]) + n[2]*(v3[2]-v0[2]) > 0.0)
	d = -d;
    n[0] *= d;
    n[1] *= d;
    n[2] *= d;
}

static void
get_cuts(const FLOAT *v0, const FLOAT *v1, BOOLEAN in_flag,
			int *cuts_n, FLOAT **cuts_r)
{
    int n, *mults;
    FLOAT *roots;

    n = intersect_edge(v0, v1, Opts.EPS, in_flag, &roots, &mults);
    *cuts_n = 0;
    *cuts_r = NULL;
    if (n > 0) {
	int k;
	*cuts_r = phgAlloc(n * sizeof(*roots));
	for (k = 0; k < n; k++) {
	    if (mults[k] % 2 == 0)
		continue;	/* skip root with even multiplicity */
	    if (roots[k] < Opts.EPS || roots[k] > 1.0 - Opts.EPS) {
		(*cuts_r)[(*cuts_n)++] = roots[k];
		continue;		/* skip v1 */
	    }
	    (*cuts_r)[(*cuts_n)++] = roots[k];
	}
    } /* if (n>0) */
    phgFree(roots);
    phgFree(mults);
    if (*cuts_n == 0) {
	phgFree(*cuts_r);
	*cuts_r = NULL;
	return;
    }

#if 0
    /* sorting the roots (unnecessary?) */
    if (*cuts_n > 1) {
	qsort(*cuts_r, *cuts_n, sizeof(**cuts_r), phgCompFLOAT);
    }
#endif

    if (Opts.dbg_off || !Opts.dbg_roots)
	return;

    for (n = 0; n < *cuts_n; n++) {
	FLOAT res;
	static FLOAT res0 = 0.0;
#if USE_OMP
# pragma omp threadprivate(res0)
#endif  /* USE_OMP */
	levelset_func(v0[0] + (*cuts_r)[n] * (v1[0] - v0[0]),
		      v0[1] + (*cuts_r)[n] * (v1[1] - v0[1]),
		      v0[2] + (*cuts_r)[n] * (v1[2] - v0[2]), &res);
	if ((res = Fabs(res)) > res0)
	    fprintf(stderr, "\troots[%d]=%lg, res=%lg\n",
			    n, (double)(*cuts_r)[n], (double)(res0 = res));
    }
}

static void
add_point(FLOAT x, FLOAT y, FLOAT z, const FLOAT *n, const FLOAT w)
/* adds a point to the quadrature rule, with w the weight, n the unit normal.
 * the type of quarrature rule is stored in QI.quad_type */
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
    QI.rules[i]->pw[pos++] = z;

    /* the weights */
    QI.rules[i]->pw[pos++] = w;

    /* the unit normal vector */
    if (QI.quad_type == 0 && Opts.nv_flag) {
	pos = QI.rules[i]->n * Dim;
	QI.rules[i]->nv[pos++] = n[0];
	QI.rules[i]->nv[pos++] = n[1];
	QI.rules[i]->nv[pos++] = n[2];
    }

    QI.rules[i]->n++;
}

static void
quad_tetra(TET *pt, FLOAT *res)
/* adds the volume integral on the whole tetrahedron to *res */
{
    QUAD *quad = phgQuadGetQuad3D(QI.quad_order);
    FLOAT x, y, z, vol, res1, *p, *w;
    int i;

    res1 = 0.0;

    vol = tet_volume(pt);

    p = quad->points;
    w = quad->weights;
    for (i = 0; i < quad->npoints; i++, p += Dim + 1, w++) {
	x = p[0] * pt->verts[0][0] + p[1] * pt->verts[1][0] +
	    p[2] * pt->verts[2][0] + p[3] * pt->verts[3][0];
	y = p[0] * pt->verts[0][1] + p[1] * pt->verts[1][1] +
	    p[2] * pt->verts[2][1] + p[3] * pt->verts[3][1];
	z = p[0] * pt->verts[0][2] + p[1] * pt->verts[1][2] +
	    p[2] * pt->verts[2][2] + p[3] * pt->verts[3][2];
	add_point(x, y, z, NULL, *w * vol);
	res1 += (*w) * ONE;
    }

    /* multiply res1 by the volume of pt and add to res */
    *res += res1 * vol;

    return;
}

static void
quad_triangle(const COORD v0, const COORD v1, const COORD v2, const COORD v3,
	      BOOLEAN ignore_direction, FLOAT *res)
/* adds the surface integral on the triangle (v0,v1,v2) to *res.
 * v3 is the opposite vertex used to determine the outward direction.
 *
 * Note: to avoid double counting, the integral is only computed when the
 * out normal of the face is in the same direction as the interface normal.
 * The argument 'ignore_direction', when is TRUE (e.g. for a boundary face),
 * disables this check so that the integral is alway computed regardless of
 * the direction of the out normal. */
{
    FLOAT n[3], c[3], oned3 = (_F(1.0) / _F(3.0));
    FLOAT *p, *w, S, d;
    QUAD *quad;
    int j;

    /* compute the out normal vector */
    get_face_normal(n, v0, v1, v2, v3);

    if (ignore_direction == FALSE) {
	FLOAT *grad_L = c;
	c[0] = (v0[0] + v1[0] + v2[0]) * oned3;
	c[1] = (v0[1] + v1[1] + v2[1]) * oned3;
	c[2] = (v0[2] + v1[2] + v2[2]) * oned3;
	levelset_grad(c[0], c[1], c[2], grad_L, FALSE);
	if (n[0] * grad_L[0] + n[1] * grad_L[1] + n[2] * grad_L[2] < 0.0)
	    return;
    }

    /* compute 0.5*dS = 0.5 * |(v1-v0)x(v2-v0)| */
    c[0] = (v1[1]-v0[1])*(v2[2]-v0[2]) - (v1[2]-v0[2])*(v2[1]-v0[1]);
    c[1] = (v1[2]-v0[2])*(v2[0]-v0[0]) - (v1[0]-v0[0])*(v2[2]-v0[2]);
    c[2] = (v1[0]-v0[0])*(v2[1]-v0[1]) - (v1[1]-v0[1])*(v2[0]-v0[0]);
    S = 0.5 * Sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    quad = phgQuadGetQuad2D(QI.quad_order);
    p = quad->points;
    w = quad->weights;
    for (j = 0; j < quad->npoints; j++, p += Dim, w++) {
	c[0] = p[0] * v0[0] + p[1] * v1[0] + p[2] * v2[0];
	c[1] = p[0] * v0[1] + p[1] * v1[1] + p[2] * v2[1];
	c[2] = p[0] * v0[2] + p[1] * v1[2] + p[2] * v2[2];
	d = S * (*w);
	add_point(c[0], c[1], c[2], n, d);
	*res += d * ONE;
    }
}

static void
quad_planar_sub(const FLOAT *v0, const FLOAT *v1, const FLOAT *v2,
		const FLOAT *v3, FLOAT Lv3, FLOAT res[3], BOOLEAN face3_flag)
/* computes quadrature in the tetrahedron (v0,v1,v2,v3).
 *
 * The tetrahedron is entirely contained in either Omega- or Omega+ according
 * to sign(Lv3).
 *
 * If face3_flag == TRUE then face 3 is on (part of) the interface */
{
    TET *pt;

    if (face3_flag && QI.rules[GM] != NULL) {
	QI.quad_type = 0;
	quad_triangle(v0, v1, v2, v3, FALSE, &res[GM]);
	QI.quad_type = 99;	/* for debugging */
    }

    if ((Lv3 < 0.0 && QI.rules[OM] == NULL) ||
	(Lv3 > 0.0 && QI.rules[OP] == NULL))
	return;

    pt = tet_new(TRUE);
    memcpy(pt->verts[0], v0, sizeof(FLOAT) * Dim);
    memcpy(pt->verts[1], v1, sizeof(FLOAT) * Dim);
    memcpy(pt->verts[2], v2, sizeof(FLOAT) * Dim);
    memcpy(pt->verts[3], v3, sizeof(FLOAT) * Dim);
    pt->volume = -1.0;
    if (Lv3 <= 0.0 && QI.rules[OM] != NULL) {
	QI.quad_type = -1;
	quad_tetra(pt, &res[OM]);
	QI.quad_type = 99;	/* for debugging */
    }
    if (Lv3 >= 0.0 && QI.rules[OP] != NULL) {
	QI.quad_type = 1;
	quad_tetra(pt, &res[OP]);
	QI.quad_type = 99;	/* for debugging */
    }
    tet_free(pt);
}

static void
quad_planar_interface(TET *pt, FLOAT res[3])
/* computes the integral when the interface is planar */
{
    FLOAT cuts[NEdge], s[NVert], *v0, *v1, *v2, *v3;
    FLOAT v4[3], v5[3], v6[3], v7[3];
    int i, j, k, l, i0, i1, i2, i3, ii;

    res[0] = res[1] = res[2] = 0.;

    for (j = 0; j < NVert; j++) {
	levelset_func(pt->verts[j][0], pt->verts[j][1], pt->verts[j][2], &s[j]);
	if (Fabs(s[j]) <= Opts.EPS)
	    s[j] = 0.0;
    }

    /* count # of "-" and "+" vertices (stored in j and k),
     * i stores the index of an arbitrary vertex not on the interface */
    j = k = 0;
    i = -1;
    for (l = 0; l < NVert; l++) {
	if (s[l] == 0.0)
	    continue;
	i = l;
	if (s[l] < 0.0)
	    j++;
	else
	    k++;
    }
    /*assert(i >= 0);*/

    if (j + k == 1) {
	/* face i is on the interface */
	v0 = pt->verts[GetFaceVertex(i, 0)];
	v1 = pt->verts[GetFaceVertex(i, 1)];
	v2 = pt->verts[GetFaceVertex(i, 2)];
	v3 = pt->verts[i];
	QI.quad_type = s[i] < 0.0 ? -1 : 1;
	quad_tetra(pt, &res[QI.quad_type + 1]);
	QI.quad_type = 0;
	quad_triangle(v0, v1, v2, v3, pt->is_bdry[i], &res[GM]);
	QI.quad_type = 99;	/* for debugging */
	return;
    }

    if (j == 0 || k == 0) {
	/* the element is entirely on one side of the interface, with
	 * at most 2 vertices on the interface. */
	QI.quad_type = s[i] < 0.0 ? -1 : 1;
	quad_tetra(pt, &res[QI.quad_type + 1]);
	QI.quad_type = 99;	/* for debugging */
	return;
    }

    /* construct cut-points in cuts[] */
    j = 0;	/* number of cuts */
    for (i = 0; i < NEdge; i++) {
	FLOAT a, b;
	a = s[GetEdgeVertex(i,0)];
	b = s[GetEdgeVertex(i,1)];
	if (a * b >= 0.0) {
	    cuts[i] = -1.0;
	    continue;
	}
#if 0
	/* bisection */
	{
	    FLOAT t0 = 0.0, t1 = 1.0, t, c;
	    v0 = pt->verts[GetEdgeVertex(i,0)];
	    v1 = pt->verts[GetEdgeVertex(i,1)];
	    while (TRUE) {
		t = (t0 + t1) * 0.5;
		if (t == t0 || t == t1)
		    break;
		levelset_func(v0[0] + t * (v1[0] - v0[0]),
			      v0[1] + t * (v1[1] - v0[1]),
			      v0[2] + t * (v1[2] - v0[2]), &c);
		if (c == 0.0)
		    break;
		if (c * a > 0.0)
		    t0 = t;
		else
		    t1 = t;
	    }
	    cuts[i] = t;
	}
#else
	/* linear approximation */
#ifdef __ICC
/* The following line is a workaround for a potential Intel compiler bug
 * (oneapi-2021.1.1 on LSSC4, with -O2). Test case:
 * 	make -s USER_CFLAGS="-DCASE=1" quad_test2
 * 	./quad_test2 -mesh_file cube.dat
 * 	PHG: *** ERROR: floating point exception: divide by zero.
 * (Note: an independent code showing this bug is saved in "test/intel-bug.c")
 *
 * Alt. solution: "-fp-speculation=safe" (TAN Guangming, tgm@ict.ac.cn)
 */
	assert(a != b);
#endif	/* defined(__ICC) */
	cuts[i] = a / (a-b);
#endif
	j++;
    }

/* The LINK_EDGES macro reorders the vertices such that the two co-planar edges
 * (v0,v1) and (v2,v3) become (v0,v1) and (v0, v2), preserving the order of the
 * two edges */
#define LINK_EDGES(v0, v1, v2, v3)				\
    if (v2 == v1) {						\
	v2 = v0; v0 = v1; v1 = v2; v2 = v3;			\
    }								\
    else if (v3 == v1) {					\
	v3 = v0; v0 = v1; v1 = v3;				\
    }								\
    else if (v2 == v0) {					\
	v2 = v3;						\
    }								\

    switch (j) {
	case 1:		/* 1 cut point, bisect the tetra */
	    /* i = get cut edge */
	    for (i = 0; i < NEdge; i++)
		if (cuts[i] > 0.0)
		    break;
	    j = 5 - i;		/* opposite edge */
	    v0 = pt->verts[i0 = GetEdgeVertex(i,0)];
	    v1 = pt->verts[i1 = GetEdgeVertex(i,1)];
	    v2 = pt->verts[GetEdgeVertex(j,0)];
	    v3 = pt->verts[GetEdgeVertex(j,1)];
	    v4[0] = v0[0] + cuts[i] * (v1[0] - v0[0]);
	    v4[1] = v0[1] + cuts[i] * (v1[1] - v0[1]);
	    v4[2] = v0[2] + cuts[i] * (v1[2] - v0[2]);
	    quad_planar_sub(v2, v3, v4, v0, s[i0], res, TRUE);
	    quad_planar_sub(v2, v3, v4, v1, s[i1], res, TRUE);
	    break;
	case 2:		/* 2 cut points, 1 tetra + 1 pyramid */
	    /* i,j = 2 cut edges */
	    for (i = 0; i < NEdge; i++)
		if (cuts[i] > 0.0)
		    break;
	    for (j = i + 1; j < NEdge; j++)
		if (cuts[j] > 0.0)
		    break;
	    assert(j < NEdge);
	    assert(i + j != 5);		/* the two edges are coplanar */
	    v0 = pt->verts[i0 = GetEdgeVertex(i,0)];
	    v1 = pt->verts[i1 = GetEdgeVertex(i,1)];
	    v4[0] = v0[0] + cuts[i] * (v1[0] - v0[0]);
	    v4[1] = v0[1] + cuts[i] * (v1[1] - v0[1]);
	    v4[2] = v0[2] + cuts[i] * (v1[2] - v0[2]);
	    v0 = pt->verts[i2 = GetEdgeVertex(j,0)];
	    v1 = pt->verts[i3 = GetEdgeVertex(j,1)];
	    v5[0] = v0[0] + cuts[j] * (v1[0] - v0[0]);
	    v5[1] = v0[1] + cuts[j] * (v1[1] - v0[1]);
	    v5[2] = v0[2] + cuts[j] * (v1[2] - v0[2]);
	    LINK_EDGES(i0, i1, i2, i3)
	    i3 = (0 + 1 + 2 + 3) - (i0 + i1 + i2);
	    v0 = pt->verts[i0];
	    v1 = pt->verts[i1];
	    v2 = pt->verts[i2];
	    v3 = pt->verts[i3];
	    /* the tetrahedron (v3,v4,v5,v0) */
	    quad_planar_sub(v3, v4, v5, v0, s[i0], res, TRUE);
	    /* the pyramid (v1,v2,v4,v5,v3) */
	    quad_planar_sub(v3, v4, v5, v1, s[i1], res, TRUE);
	    quad_planar_sub(v3, v1, v5, v2, s[i2], res, FALSE);
	    break;
	case 3:		/* 3 cut points, 1 tetra + 1 prism */
	    for (i = 0; i < NEdge; i++)
		if (cuts[i] > 0.0)
		    break;
	    for (j = i + 1; j < NEdge; j++)
		if (cuts[j] > 0.0)
		    break;
	    for (k = j + 1; k < NEdge; k++)
		if (cuts[k] > 0.0)
		    break;
	    assert(k < NEdge);
	    /* the three edges must share a common vertex */
	    v0 = pt->verts[i0 = GetEdgeVertex(i,0)];
	    v1 = pt->verts[i1 = GetEdgeVertex(i,1)];
	    v4[0] = v0[0] + cuts[i] * (v1[0] - v0[0]);
	    v4[1] = v0[1] + cuts[i] * (v1[1] - v0[1]);
	    v4[2] = v0[2] + cuts[i] * (v1[2] - v0[2]);
	    v0 = pt->verts[i2 = GetEdgeVertex(j,0)];
	    v1 = pt->verts[i3 = GetEdgeVertex(j,1)];
	    v5[0] = v0[0] + cuts[j] * (v1[0] - v0[0]);
	    v5[1] = v0[1] + cuts[j] * (v1[1] - v0[1]);
	    v5[2] = v0[2] + cuts[j] * (v1[2] - v0[2]);
	    v0 = pt->verts[GetEdgeVertex(k,0)];
	    v1 = pt->verts[GetEdgeVertex(k,1)];
	    v6[0] = v0[0] + cuts[k] * (v1[0] - v0[0]);
	    v6[1] = v0[1] + cuts[k] * (v1[1] - v0[1]);
	    v6[2] = v0[2] + cuts[k] * (v1[2] - v0[2]);
	    LINK_EDGES(i0, i1, i2, i3)
	    i3 = (0 + 1 + 2 + 3) - (i0 + i1 + i2);
	    assert(i3 == GetEdgeVertex(k,0) || i3 == GetEdgeVertex(k,1));
	    v0 = pt->verts[i0];
	    v1 = pt->verts[i1];
	    v2 = pt->verts[i2];
	    v3 = pt->verts[i3];
	    /* the tetrahedron (v4,v5,v6,v0) */
	    quad_planar_sub(v4, v5, v6, v0, s[i0], res, TRUE);
	    /* the prism (v4,v5,v6)-(v1,v2,v3) */
	    /* - the tetra (v4,v5,v6,v1) */
	    quad_planar_sub(v4, v5, v6, v1, s[i1], res, TRUE);
	    /* - the pyramid (v2,v3,v5,v6,v1), with no face on Gamma */
	    quad_planar_sub(v2, v5, v6, v1, s[i1], res, FALSE);
	    quad_planar_sub(v2, v3, v6, v1, s[i1], res, FALSE);
	    break;
	case 4:		/* 4 cut points, 1 prism + 1 prism */
	    for (i = 0; i < NEdge; i++)
		if (cuts[i] > 0.0)
		    break;
	    for (j = i + 1; j < NEdge; j++)
		if (cuts[j] > 0.0)
		    break;
	    for (k = j + 1; k < NEdge; k++)
		if (cuts[k] > 0.0)
		    break;
	    for (l = k + 1; l < NEdge; l++)
		if (cuts[l] > 0.0)
		    break;
	    assert(l < NEdge);
	    assert(i + j + k + l == 10);	/* 2 pairs of opposite edges */
	     /* reorder the edges as (v0,v1), (v1,v2), (v2,v3), (v3,v0) */
	    i0 = GetEdgeVertex(i,0);
	    i1 = GetEdgeVertex(i,1);
	    if (j == 5 - i) {
		/* edge j is opposite to edge i, exchange edges j and k */
		i2 = k;
		k = j;
		j = i2;
		assert(j != 5 - i);
	    }
	    i2 = GetEdgeVertex(j, 0);
	    i3 = GetEdgeVertex(j, 1);
	    LINK_EDGES(i1, i0, i2, i3)
	    if (k == 5 - j) {
		/* edge k is opposite to edge j, exchange edges k and l */
		i3 = l;
		l = k;
		k = i3;
		assert(k != 5 - j);
	    }
	    i3 = GetEdgeVertex(k, 0);
	    ii = GetEdgeVertex(k, 1);
	    assert(i3 != i1 && ii != i1);
	    LINK_EDGES(i2, i1, i3, ii)
	    /* now edge l must equal to (i3,i0) or (i0,i3) */
	    assert(((ii=GetEdgeVertex(l, 0)) == i0 || (ii == i3)) &&
                   ((ii=GetEdgeVertex(l, 1)) == i0 || (ii == i3)));
	    /* compute the four new vertices */
	    v0 = pt->verts[GetEdgeVertex(i, 0)];
	    v1 = pt->verts[GetEdgeVertex(i, 1)];
	    v4[0] = v0[0] + cuts[i] * (v1[0] - v0[0]);
	    v4[1] = v0[1] + cuts[i] * (v1[1] - v0[1]);
	    v4[2] = v0[2] + cuts[i] * (v1[2] - v0[2]);
	    v0 = pt->verts[GetEdgeVertex(j, 0)];
	    v1 = pt->verts[GetEdgeVertex(j, 1)];
	    v5[0] = v0[0] + cuts[j] * (v1[0] - v0[0]);
	    v5[1] = v0[1] + cuts[j] * (v1[1] - v0[1]);
	    v5[2] = v0[2] + cuts[j] * (v1[2] - v0[2]);
	    v0 = pt->verts[GetEdgeVertex(k, 0)];
	    v1 = pt->verts[GetEdgeVertex(k, 1)];
	    v6[0] = v0[0] + cuts[k] * (v1[0] - v0[0]);
	    v6[1] = v0[1] + cuts[k] * (v1[1] - v0[1]);
	    v6[2] = v0[2] + cuts[k] * (v1[2] - v0[2]);
	    v0 = pt->verts[GetEdgeVertex(l, 0)];
	    v1 = pt->verts[GetEdgeVertex(l, 1)];
	    v7[0] = v0[0] + cuts[l] * (v1[0] - v0[0]);
	    v7[1] = v0[1] + cuts[l] * (v1[1] - v0[1]);
	    v7[2] = v0[2] + cuts[l] * (v1[2] - v0[2]);
	    v0 = pt->verts[i0];
	    v1 = pt->verts[i1];
	    v2 = pt->verts[i2];
	    v3 = pt->verts[i3];
	    /* Prism v2,v6,v5, v0,v7,v4 */
	    /* Prism v1,v5,v4, v3,v6,v7 */
	    quad_planar_sub(v6, v5, v7, v2, s[i2], res, TRUE);
	    quad_planar_sub(v5, v4, v7, v2, s[i2], res, TRUE);
	    quad_planar_sub(v5, v4, v6, v1, s[i1], res, TRUE);
	    quad_planar_sub(v4, v7, v6, v1, s[i1], res, TRUE);
	    /* note: s[0]*s[2] > 0, s[1]*s[3] > 0 */
	    quad_planar_sub(v4, v7, v2, v0, s[i0], res, FALSE);
	    quad_planar_sub(v7, v6, v1, v3, s[i3], res, FALSE);
	    break;
	default:
	    phgError(1, "%s:%d: unexpected case %d.\n", __FILE__, __LINE__, j);
    }
}

static FLOAT w_s = 1.0, w_t = 1.0;
#if USE_OMP
# pragma omp threadprivate(w_s, w_t)
#endif  /* USE_OMP */

static void
quad_func_r(FLOAT r, FLOAT *res, const FLOAT w, void *pctx)
{
    QUAD_CONTEXT *qc = pctx;

    add_point(qc->u_st[0] + r * qc->nu[0], qc->u_st[1] + r * qc->nu[1],
		qc->u_st[2] + r * qc->nu[2], NULL, w * w_s * w_t);
    *res = ONE;
}

#if 0
/* GCC pragmas and attribute controlling optimization of functions */
#pragma GCC push_options
#pragma GCC optimize ("-O0")
#pragma GCC pop_options

static void quad_func_s() __attribute__((optimize("-O2")));
#endif

static BOOLEAN flag_GM = TRUE;
#if USE_OMP
# pragma omp threadprivate(flag_GM)
#endif	/* USE_OMP */

static void
quad_func_s(FLOAT s, FLOAT res[3], const FLOAT w, void *pctx)
/* quadrature along the direction qc->nu */
{
    QUAD_CONTEXT *qc = pctx;
    COORD u_st, tmp;
    int j, n, flag;
    FLOAT a, b, d, sign, res0[3], *roots, cuts[4 + 1];

    res[0] = res[1] = res[2] = 0.0;

    w_s = w;

    u_st[0] = qc->v_t[0] + s * qc->nv[0];
    u_st[1] = qc->v_t[1] + s * qc->nv[1];
    u_st[2] = qc->v_t[2] + s * qc->nv[2];

    tmp[0] = u_st[0] + qc->nu[0];
    tmp[1] = u_st[1] + qc->nu[1];
    tmp[2] = u_st[2] + qc->nu[2];

    /* intersection with the interface */
    get_cuts(u_st, tmp, TRUE, &n, &roots);
#if 0
    fprintf(stderr, "r integral range = [%lg %lg %lg] - [%lg %lg %lg]\n",
			(double)u_st[0], (double)u_st[1], (double)u_st[2],
			(double)tmp[0], (double)tmp[1], (double)tmp[2]);
#endif
    if (n > 0)
	cuts[0] = roots[0];
    phgFree(roots);
    if (n > 1)
	LONGJMP;

    if (n == 1 && flag_GM && QI.rules[GM] != NULL) {
	/* surface integral */
	COORD G, tmp1;
	d = cuts[0];
	tmp1[0] = u_st[0] + d * qc->nu[0];
	tmp1[1] = u_st[1] + d * qc->nu[1];
	tmp1[2] = u_st[2] + d * qc->nu[2];
	/* Test unnecessary, the interval is already checked in quad_func_t() */
	/*if (levelset_face(0, tmp1[0], tmp1[1], tmp1[2]) > -Opts.EPS &&
	    levelset_face(1, tmp1[0], tmp1[1], tmp1[2]) > -Opts.EPS &&
	    levelset_face(2, tmp1[0], tmp1[1], tmp1[2]) > -Opts.EPS &&
	    (Opts.use_face_normal != FALSE ||
	    levelset_face(3, tmp1[0], tmp1[1], tmp1[2]) > -Opts.EPS))*/ {
	    levelset_grad(tmp1[0], tmp1[1], tmp1[2], G, TRUE);
	    d = (G[0]*qc->nu[0] + G[1]*qc->nu[1] + G[2]*qc->nu[2]) / qc->hu;
	    if (Fabs(d) < 0.01)
		LONGJMP;
	    d = 1.0 / Fabs(d);
	    QI.quad_type = 0;
	    add_point(tmp1[0], tmp1[1], tmp1[2], G, w_s * w_t * d);
	    QI.quad_type = 99;	/* for debugging */
	    res[GM] = d * ONE;
	}
    }

    /* construct list of cut-points for the inner integral (r) */
    if (QI.rules[OM] == NULL && QI.rules[OP] == NULL)
	return;

    /* intersections with the faces:
     *	    the three faces containing v3 if Opts.use_face_normal == TRUE, and
     *	    all faces otherwise. */
    for (j = 0; j < (Opts.use_face_normal ? 3 : NFace); j++) {
	a = levelset_face(j, u_st[0], u_st[1], u_st[2]);
	b = levelset_face(j, tmp[0], tmp[1], tmp[2]);
	if (a * b < -Opts.EPS)
	    cuts[n++] = a / (a - b);
    }

    if (n > 1)
	qsort(cuts, n, sizeof(cuts[0]), phgCompFLOAT);

    flag = 0;
    qc->u_st = u_st;
    res[0] = res[1] = res[2] = 0.0;
    for (j = 0; j <= n; j++) {
	a = (j == 0 ? 0.0 : cuts[j - 1]);
	b = (j == n ? 1.0 : cuts[j]);
	d = (a + b) * 0.5;
	tmp[0] = u_st[0] + d * qc->nu[0];
	tmp[1] = u_st[1] + d * qc->nu[1];
	tmp[2] = u_st[2] + d * qc->nu[2];
	if (levelset_face(0, tmp[0], tmp[1], tmp[2]) < -Opts.EPS ||
	    levelset_face(1, tmp[0], tmp[1], tmp[2]) < -Opts.EPS ||
	    levelset_face(2, tmp[0], tmp[1], tmp[2]) < -Opts.EPS ||
	   (Opts.use_face_normal == FALSE &&
	    levelset_face(3, tmp[0], tmp[1], tmp[2]) < -Opts.EPS))
	    continue;
	levelset_func(tmp[0], tmp[1], tmp[2], &sign);
	flag = 1;
	if (b - a < Opts.EPS)
	    continue;
	QI.quad_type = sign <= 0.0 ? -1 : 1;
	phgQuad1D(quad_func_r, 1, a, b,
		  phgQuadGetQuad1D(adapt_order((b - a) * qc->hu)),
						&res0[QI.quad_type + 1], qc);
	res[QI.quad_type + 1] += res0[QI.quad_type + 1];
	QI.quad_type = 99;	/* for debugging */
    }

    qc->u_st = NULL;

    /* test cases:
     *		-mesh_file ellipsoid.dat -order 15 -refine 3 [-+]qi_shortcut_s
     */
    if (Opts.shortcut_s == TRUE && !flag && Opts.graph_s <= 0)
	longjmp(qc->jmp_env_s, 1);
}

/* Projected length of 'v-v0' on the vector 'n' with |n|='nnorm' */
#ifdef ProjectionOnVector
# undef ProjectionOnVector
#endif
#define ProjectionOnVector(v, v0, n, nnorm_square) \
    	(((v[0]-v0[0]) * n[0] + (v[1]-v0[1]) * n[1] + (v[2]-v0[2]) * n[2]) \
		/ (nnorm_square))

static FLOAT
inter_edge(const FLOAT *v_t, int edge)
{
    QUAD_CONTEXT *qc = &QI.ctx;
#if 0
    /* returns a such that v_t+a*qc->nv is the intersection point of the line
     * segment (v_t,v_t+qc->nv) with the edge 'edge' projected to the base
     * plane along the direction qc->nu, or -1 if no intersection. */
    COORD u0, u1;
    FLOAT a11, a12, a21, a22, b1, b2, d, t, s;
    const FLOAT *v;

    /* project the two vertices of the edge to the base plane:
     *		v -= d * qc->nu
     * where
     *		d = (v - qc->x0, qc->nu) / (qc->hu * qc->hu)
     */
    v = qc->verts[GetEdgeVertex(edge, 0)];
    d = ((v[0] - qc->x0[0]) * qc->nu[0] +
	 (v[1] - qc->x0[1]) * qc->nu[1] +
	 (v[2] - qc->x0[2]) * qc->nu[2]) / (qc->hu * qc->hu);
    u0[0] = v[0] - d * qc->nu[0];
    u0[1] = v[1] - d * qc->nu[1];
    u0[2] = v[2] - d * qc->nu[2];

    v = qc->verts[GetEdgeVertex(edge, 1)];
    d = ((v[0] - qc->x0[0]) * qc->nu[0] +
	 (v[1] - qc->x0[1]) * qc->nu[1] +
	 (v[2] - qc->x0[2]) * qc->nu[2]) / (qc->hu * qc->hu);
    u1[0] = v[0] - d * qc->nu[0];
    u1[1] = v[1] - d * qc->nu[1];
    u1[2] = v[2] - d * qc->nu[2];

    a11 = qc->hv * qc->hv;
    a12 = qc->nv[0] * (u0[0] - u1[0]) +
	  qc->nv[1] * (u0[1] - u1[1]) +
	  qc->nv[2] * (u0[2] - u1[2]);
    b1  = qc->nv[0] * (u0[0] - v_t[0]) +
	  qc->nv[1] * (u0[1] - v_t[1]) +
	  qc->nv[2] * (u0[2] - v_t[2]);
    a21 = a12;
    a22 = (u0[0] - u1[0]) * (u0[0] - u1[0]) +
	  (u0[1] - u1[1]) * (u0[1] - u1[1]) +
	  (u0[2] - u1[2]) * (u0[2] - u1[2]);
    b2  = (u0[0] - u1[0]) * (u0[0] - v_t[0]) +
	  (u0[1] - u1[1]) * (u0[1] - v_t[1]) +
	  (u0[2] - u1[2]) * (u0[2] - v_t[2]);
    d = a11 * a22 - a12 * a21;
    assert(d >= -Opts.EPS);
    if (d <= /*Opts.EPS * (Fabs(a11)+Fabs(a12)+Fabs(a21)+Fabs(a22))*/0.0)
	return -1.0;
    t = (b1 * a22 - b2 * a12);
    if (t <= 0.0 || t >= d)
	return -1.0;
    s = (b1 * a22 - b2 * a12);
    if (s <= 0.0 || s >= d)
	return -1.0;

    return t / d;
#else
    /* returns a such that v_t+a*qc->nv is the projection to (v_t,v_t+nv) of
     * the intersection of the edge 'edge' with the plane:
     * 		{v_t + x\cross qc->nw | x \in \R^3}
     *
     * This code is faster and should be more robust than above.
     */
    const FLOAT *v0, *v1;
    FLOAT a, b, d;
    COORD u;

    v0 = qc->verts[GetEdgeVertex(edge, 0)];
    v1 = qc->verts[GetEdgeVertex(edge, 1)];

    /* a := (v0-v_t, qc->nu) */
    a = (v0[0] - v_t[0]) * qc->nw[0] +
	(v0[1] - v_t[1]) * qc->nw[1] +
	(v0[2] - v_t[2]) * qc->nw[2];

    /* b := (v1-v_t, qc->nu) */
    b = (v1[0] - v_t[0]) * qc->nw[0] +
	(v1[1] - v_t[1]) * qc->nw[1] +
	(v1[2] - v_t[2]) * qc->nw[2];

    if (a * b > 0.0)
	return -1.0;	/* no intersection */

    d = b - a;
    if (d == 0.0)
	return -1.0;
    d = 1.0 / d;

    u[0] = (v0[0] * b - v1[0] * a) * d;
    u[1] = (v0[1] * b - v1[1] * a) * d;
    u[2] = (v0[2] * b - v1[2] * a) * d;
    d = ProjectionOnVector(u, v_t, qc->nv, qc->hv * qc->hv);

    return d;
#endif
}

static BOOLEAN
proj_to_face(FLOAT *v, const FLOAT *n, int face)
/* updates the vector 'v' with its projection on the face 'face' along the
 * direction 'n', returns TRUE if successful and FALSE if 'n' is parallel to
 * the face. */
{
    FLOAT a, b;

    a = levelset_face(face, v[0], v[1], v[2]);
    b = levelset_face(face, v[0] + n[0], v[1] + n[1], v[2] + n[2]);
    if (Fabs(a - b) < Opts.EPS)
	return FALSE;
    a /= a - b;
    v[0] = v[0] + a * n[0];
    v[1] = v[1] + a * n[1];
    v[2] = v[2] + a * n[2];

    return TRUE;
}

static void
quad_func_t(FLOAT t, FLOAT res[3], const FLOAT w, void *pctx)
/* quadrature along the direction qc->nv */
{
    QUAD_CONTEXT *qc = pctx;
    COORD v_t, v1, tmp, tmp1;
    int i, j, n, m, flag, n_save[3];
    FLOAT a, b, c, d, res0[3], cuts[20], *roots;

    res[0] = res[1] = res[2] = 0.0;

    w_t = w;

    v_t[0] = qc->x0[0] + t * qc->nw[0];
    v_t[1] = qc->x0[1] + t * qc->nw[1];
    v_t[2] = qc->x0[2] + t * qc->nw[2];

    v1[0] = v_t[0] + qc->nv[0];
    v1[1] = v_t[1] + qc->nv[1];
    v1[2] = v_t[2] + qc->nv[2];

    qc->v_t = v_t;

    if (!Opts.dbg_off && Opts.dbg_value_s >= 0.0 && Opts.dbg_value_s <= 1.0) {
	/* debugging integral on s for a given value of t */
	QI.skip_rule = TRUE;
	flag_GM = TRUE;
	quad_func_s(Opts.dbg_value_s, res, 1.0, qc);
	QI.skip_rule = FALSE;
	printf("g(%0.16lg,%0.16lg) = %0.16lg\n",
			(double)t, (double)Opts.dbg_value_s, (double)res[0]);
	qc->v_t = NULL;
	return;
    }

    if (Opts.graph_s > 0) {
	assert(Opts.graph_t > 0);
	flag_GM = TRUE;
	QI.skip_rule = TRUE;
	for (i = 0; i <= Opts.graph_s; i++) {
	    a = i / (FLOAT)Opts.graph_s;
	    quad_func_s(a, res0, 1.0, qc);
	    b = res0[0];
	    QI.quad_order += 2;
	    quad_func_s(a, res0, 1.0, qc);
	    QI.quad_order -= 2;
	    printf("%0.16lg %0.16lg %0.16lg %0.16lg\n",
			(double)t, (double)a, (double)b,
			(double)Fabs((b - res0[0]) / (b == 0.0 ? 1.0 : b)));
	}
	QI.skip_rule = FALSE;
	qc->v_t = NULL;
	return;
    }

    /* construct list of cut-points for the middle integral (s) */

    n = 0;

    /* intersections with (v_t,v1) of the 6 edges projected to the base plane,
     * or equivalently the projections to (v_t,v1) of the intersections of
     * the 6 edges with the (nu,nv) plane at v_t */
    for (i = 0; i < NEdge; i++)
	if ((a = inter_edge(v_t, i)) > 0.0)
	    cuts[n++] = a;

    /* Intersections with the projections on the base plane along the direction
     * qc->nu of the interface on the faces.
     *
     * Note: instead of projecting the interface to the base plane and
     * computing its intersections with (v_t,v1), we project v_t and v1 to
     * each face, compute the intersections with the interface in the face,
     * and then project the intersection points back to the base plane.
     */

#ifdef IntersectionsOnFace
# undef IntersectionsOnFace
#endif
#define IntersectionsOnFace(face) \
    memcpy(tmp, v_t, sizeof(tmp)); \
    memcpy(tmp1, v1, sizeof(tmp1)); \
    if (proj_to_face(tmp, qc->nu, face) == TRUE && \
	proj_to_face(tmp1, qc->nu, face) == TRUE) { \
	/* find the segment [a,b] within the triangle */ \
	tmp1[0] -= tmp[0]; \
	tmp1[1] -= tmp[1]; \
	tmp1[2] -= tmp[2]; \
	b = (a = 0.0) + 1.0; \
	for (j = 0; j < 3; j++) { \
	    c = levelset_face(GetFaceVertex(face, j), \
				  tmp[0] + a * tmp1[0], \
				  tmp[1] + a * tmp1[1], \
				  tmp[2] + a * tmp1[2]); \
	    d = levelset_face(GetFaceVertex(face, j), \
				  tmp[0] + b * tmp1[0], \
				  tmp[1] + b * tmp1[1], \
				  tmp[2] + b * tmp1[2]); \
	    if (c < 0.0 && d < 0.0) { a = (b = 0.0) + 1.0; break; } \
	    if (c >= 0.0 && d >= 0.0) continue; \
	    if (c < 0.0) a = (d * a - c * b) / (d - c); \
	    else b = (d * a - c * b) / (d - c); \
	    if (a >= b) break; \
	} \
	if (a < b) { \
	    tmp[0] += a * tmp1[0]; \
	    tmp[1] += a * tmp1[1]; \
	    tmp[2] += a * tmp1[2]; \
	    tmp1[0] = tmp[0] + (b - a) * tmp1[0]; \
	    tmp1[1] = tmp[1] + (b - a) * tmp1[1]; \
	    tmp1[2] = tmp[2] + (b - a) * tmp1[2]; \
	    get_cuts(tmp, tmp1, TRUE, &m, &roots); \
	    if (m > 1) { \
		phgFree(roots); \
		LONGJMP; /* multiple cuts with the trace of the interface */ \
	    } \
	    if (m > 0) { \
		assert(n < sizeof(cuts) / sizeof(cuts[0])); \
		cuts[n++] = a + roots[0] * (b - a); \
	    } \
	    phgFree(roots); \
	} \
    }

    IntersectionsOnFace(0)
    IntersectionsOnFace(1)
    IntersectionsOnFace(2)
    IntersectionsOnFace(3)
#undef IntersectionsOnFace

    if (n > 1)
	qsort(cuts, n, sizeof(cuts[0]), phgCompFLOAT);

    flag = 0;
    res[0] = res[1] = res[2] = 0.0;
    for (j = 0; j <= n; j++) {
	a = (j == 0 ? 0.0 : cuts[j - 1]);
	b = (j == n ? 1.0 : cuts[j]);
	if (!Opts.dbg_off && Opts.show_intervals_s && !Opts.dbg_s)
	    fprintf(stderr, "quad_s: t=0.16%lg, interval (%0.16lg,%0.16lg)\n",
			(double)t, (double)a, (double)b);
	if (b - a < Opts.EPS)
	    continue;
	flag_GM = TRUE;
	/*if (QI.rules[GM] != NULL)*/ {
	    /* for surface integral, check that \Gamma\cap T \ne \emptyset */
	    d = (a + b) * 0.5;
	    v1[0] = v_t[0] + d * qc->nv[0];
	    v1[1] = v_t[1] + d * qc->nv[1];
	    v1[2] = v_t[2] + d * qc->nv[2];
	    tmp[0] = v1[0] + qc->nu[0];
	    tmp[1] = v1[1] + qc->nu[1];
	    tmp[2] = v1[2] + qc->nu[2];
	    get_cuts(v1, tmp, TRUE, &m, &roots);
	    if (m != 1) {
		phgFree(roots);
		if (m > 1)
		    LONGJMP;
		flag_GM = FALSE;
	    }
	    else {
		d = roots[0];
		phgFree(roots);
		v1[0] += d * qc->nu[0];
		v1[1] += d * qc->nu[1];
		v1[2] += d * qc->nu[2];
		if (levelset_face(0, v1[0], v1[1], v1[2]) < -Opts.EPS ||
		    levelset_face(1, v1[0], v1[1], v1[2]) < -Opts.EPS ||
		    levelset_face(2, v1[0], v1[1], v1[2]) < -Opts.EPS ||
		   (Opts.use_face_normal == FALSE &&
		    levelset_face(3, v1[0], v1[1], v1[2]) < -Opts.EPS))
		    flag_GM = FALSE;
	    }
	}
	for (i = 0; i < 3; i++)
	    if (QI.rules[i] != NULL)
		n_save[i] = QI.rules[i]->n;	/* save current position */
	if ((i = setjmp(qc->jmp_env_s)) != 0) {
	    for (i = 0; i < 3; i++)
		if (QI.rules[i] != NULL)
		    QI.rules[i]->n = n_save[i];		/* restore possition */
	    continue;
	}
	phgQuad1D(quad_func_s, 3, a, b,
		  get_quad_1D(adapt_order((b - a) * qc->hv)), res0, qc);
	flag = 1;

	if (!Opts.dbg_off && Opts.dbg_s) {
	    static FLOAT error0 = 0.0;
#if USE_OMP
# pragma omp threadprivate(error0)
#endif  /* USE_OMP */
	    FLOAT error, res1[3];
	    QI.skip_rule = TRUE;
	    phgQuad1D(quad_func_s, 3, a, b,
			get_quad_1D(adapt_order((b - a) * qc->hv) + 2),
			res1, qc);
	    QI.skip_rule = FALSE;
	    error = Fabs((res0[0] - res1[0]) / (res1[0] == 0. ? 1. : res1[0]));
	    if (error > error0 || Opts.show_intervals_s) {
		fprintf(stderr, "quad_s: t=%lg, (%lg,%lg), res=%0.4le, "
				"|error|=%0.4le\n",
				(double)t, (double)a, (double)b,
				(double)res0[0], (double)error);
		if (error > error0)
		    output_vtk3("tmp.vtk", qc->verts[0], qc->verts[1],
				qc->verts[2], qc->verts[3]);
		error0 = error;
	    }
	}

	for (i = 0; i < 3; i++)
	    res[i] += res0[i];
    }

    qc->v_t = NULL;

    if (Opts.shortcut_t == TRUE && !flag && Opts.graph_t <= 0)
	longjmp(qc->jmp_env_t, 1);
}

static void
quad_start(FLOAT res[3])
{
    int i, j, n, m, n_save[3], n_last[3];
    COORD tmp;
    FLOAT a = 0., b, scale, scale0, scale1, cuts[32], res0[3];
    QUAD_CONTEXT *qc = &QI.ctx;
    /* convenience variables */
    const FLOAT *v0, *v1;

    if (!Opts.dbg_off && Opts.dbg_value_t >= 0.0 && Opts.dbg_value_t <= 1.0) {
	/* debugging integral on s for a given value of t */
	QI.skip_rule = TRUE;
	quad_func_t(Opts.dbg_value_t, res, 1.0, qc);
	QI.skip_rule = FALSE;
	printf("h(%0.16lg) = %0.16lg\n",
			(double)Opts.dbg_value_t, (double)res[0]);
	return;
    }

    if (Opts.graph_t > 0) {
	printf("==================== start dumping %s\n",
					Opts.graph_s > 0 ? "h(s,t)" : "g(t)");
	QI.skip_rule = TRUE;
	for (i = 0; i <= Opts.graph_t; i++) {
	    a = i / (FLOAT)Opts.graph_t;
	    quad_func_t(a, res0, 1.0, qc);
	    if (Opts.graph_s <= 0) {
		b = res0[0];
		QI.quad_order += 2;
		quad_func_t(a, res0, 1.0, qc);
		QI.quad_order -= 2;
		printf("%0.16lg %0.16lg %0.16le\n", (double)a, (double)b,
			(double)Fabs((b - res0[0]) / (b == 0.0 ? 1.0 : b)));
	    }
	}
	QI.skip_rule = FALSE;
	printf("==================== end dumping %s\n",
					Opts.graph_s > 0 ? "h(s,t)" : "g(t)");
	return;
    }

    /* find cut-points for the outer integral (t) */

    n = 0;

    /* The projections of vertices */
    cuts[n++] = ProjectionOnVector(qc->verts[2], qc->x0, qc->nw, qc->hw*qc->hw);
    cuts[n++] = ProjectionOnVector(qc->verts[3], qc->x0, qc->nw, qc->hw*qc->hw);
    cuts[n++] = ProjectionOnVector(qc->verts[0], qc->x0, qc->nw, qc->hw*qc->hw);
    cuts[n++] = ProjectionOnVector(qc->verts[1], qc->x0, qc->nw, qc->hw*qc->hw);

    /* projections of all cut-points on all edges */
    for (i = 0; i < NEdge; i++) {
	v0 = qc->verts[GetEdgeVertex(i, 0)];
	v1 = qc->verts[GetEdgeVertex(i, 1)];
	m = qc->edge_cuts[i].n;
	for (j = 0; j < m; j++) {
	    a = qc->edge_cuts[i].cuts[j];
	    tmp[0] = v0[0] + a * (v1[0] - v0[0]);
	    tmp[1] = v0[1] + a * (v1[1] - v0[1]);
	    tmp[2] = v0[2] + a * (v1[2] - v0[2]);
	    if (n >= sizeof(cuts) / sizeof(cuts[0]))
		phgError(1, "%s:%d: cuts[] buffer overflow!\n",
						__FILE__, __LINE__);
	    cuts[n++] = ProjectionOnVector(tmp, qc->x0, qc->nw, qc->hw*qc->hw);
	}
    }

    if (n > 1)
	qsort(cuts, n, sizeof(cuts[0]), phgCompFLOAT);

    scale0 = qc->hv * qc->hw;
    scale1 = qc->hu * qc->hv * qc->hw;

    res[0] = res[1] = res[2] = 0.0;

    for (j = 0; j < 3; j++)
	if (QI.rules[j] != NULL)
	    n_last[j] = QI.rules[j]->n;

    for (j = 0; j <= n; j++) {
	a = (j == 0 ? 0.0 : cuts[j - 1]);
	b = (j == n ? 1.0 : cuts[j]);
	if (!Opts.dbg_off && Opts.show_intervals_t && !Opts.dbg_t)
	    fprintf(stderr, "quad_t: interval (%0.16lg,%0.16lg)\n",
						(double)a, (double)b);
	if (b - a < Opts.EPS)
	    continue;
	for (i = 0; i < 3; i++)
	    if (QI.rules[i] != NULL)
		n_save[i] = QI.rules[i]->n;	/* save current position */
	if (setjmp(qc->jmp_env_t) != 0) {
	    for (i = 0; i < 3; i++)
		if (QI.rules[i] != NULL)
		    QI.rules[i]->n = n_save[i];		/* restore possition */
	    continue;
	}
	/* TODO: use lower order rules for smaller intervals */
	phgQuad1D(quad_func_t, 3, a, b,
		  get_quad_1D(adapt_order((b - a) * qc->hw)), res0, qc);

	if (!Opts.dbg_off && Opts.dbg_t) {
	    static FLOAT error0 = 0.0;
#if USE_OMP
# pragma omp threadprivate(error0)
#endif  /* USE_OMP */
	    FLOAT error, res1[3];
	    QI.skip_rule = TRUE;
	    phgQuad1D(quad_func_t, 3, a, b,
			get_quad_1D(adapt_order((b - a) * qc->hw) + 2),
			res1, qc);
	    QI.skip_rule = FALSE;
	    error = Fabs((res0[0] - res1[0]) / (res1[0] == 0. ? 1. : res1[0]));
	    if (error > error0 || Opts.show_intervals_t) {
		fprintf(stderr, "quad_t: (%lg,%lg), res=%le, |error|=%le\n",
			(double)a, (double)b, (double)(res0[0] * scale1),
			(double)error);
		if (error > error0)
		    output_vtk3("tmp.vtk", qc->verts[0], qc->verts[1],
				qc->verts[2], qc->verts[3]);
		error0 = error;
	    }
	}

	for (i = 0; i < 3; i++)
	    res[i] += res0[i];
    }

    for (j = 0; j < 3; j++) {
	if (QI.rules[j] == NULL)
	    continue;
	scale = j == GM ? scale0 : scale1;
	for (i = n_last[j]; i < QI.rules[j]->n; i++)
	    QI.rules[j]->pw[RULE_HEADER + i * (Dim + 1) + Dim] *= scale;
    }

    for (i = 0; i < 3; i++) {
	if (QI.rules[i] == NULL)
	    continue;
	scale = i == GM ? scale0 : scale1;
	res[i] *= scale;
    }
}

#ifdef GetFaceEdgeNo
# undef GetFaceEdgeNo
#endif
#define GetFaceEdgeNo(face, edge) \
	GetEdgeNo(GetFaceVertex(face, edge == 0 ? 1 : 0), \
		  GetFaceVertex(face, edge == 2 ? 1 : 2))

static void
get_face_edge_cuts(void)
{
    int i, n, m, f, e;
    COORD w;
    FLOAT d, r, s, *t;
    QUAD_CONTEXT *qc = &QI.ctx;
    const FLOAT *v0, *v1, *v2;

    for (f = 0; f < NFace; f++)
	qc->face_cuts[f].n = 0;

    for (e = 0; e < NEdge; e++) {
	FLOAT *roots;
	get_cuts(qc->verts[GetEdgeVertex(e,0)],
		 qc->verts[GetEdgeVertex(e,1)], FALSE, &m, &roots);
	/* Note: 5 - e is the opposite edge */
	if (qc->face_cuts[GetEdgeVertex(5 - e, 0)].n > MAX_FACE_CUTS - m ||
	    qc->face_cuts[GetEdgeVertex(5 - e, 1)].n > MAX_FACE_CUTS - m
#if 0
	    /* Unterminated recursion if an edge is entirely on the interface.
	     * test case: the cylinder, -mesh_file cube.dat -refine 6 */
	    || (m > 1 && (/* one of 3 vertical edges */e == 2 || e >= 4))
#endif
	) {
	    phgFree(roots);
	    LONGJMP;
	}
	qc->edge_cuts[e].n = m;
	for (i = 0; i < m; i++) {
	    qc->edge_cuts[e].cuts[i] = roots[i];
	    f = GetEdgeVertex(5 - e, 0);
	    qc->face_cuts[f].info[qc->face_cuts[f].n++] = i * NEdge + e;
	    f = GetEdgeVertex(5 - e, 1);
	    qc->face_cuts[f].info[qc->face_cuts[f].n++] = i * NEdge + e;
	}
	phgFree(roots);
    }

    /* compute tangent vectors on the faces */
    for (f = 0; f < NFace; f++) {
	for (i = n = 0; i < qc->face_cuts[f].n; i++) {
	    int j;
	    e = (m = qc->face_cuts[f].info[i]) % NEdge;
	    m /= NEdge;
	    v0 = qc->verts[GetEdgeVertex(e,0)];
	    v1 = qc->verts[GetEdgeVertex(e,1)];
	    w[0] = v0[0] + qc->edge_cuts[e].cuts[m] * (v1[0] - v0[0]);
	    w[1] = v0[1] + qc->edge_cuts[e].cuts[m] * (v1[1] - v0[1]);
	    w[2] = v0[2] + qc->edge_cuts[e].cuts[m] * (v1[2] - v0[2]);

	    t = qc->face_cuts[f].tv[n];

	    /* remove duplicate cuts */
	    for (j = i + 1; j < qc->face_cuts[f].n; j++) {
		e = (m = qc->face_cuts[f].info[j]) % NEdge;
		m /= NEdge;
		v0 = qc->verts[GetEdgeVertex(e,0)];
		v1 = qc->verts[GetEdgeVertex(e,1)];
		/* note: use qc->face_cuts[f].tv[n] as work space */
		t[0] = v0[0]+qc->edge_cuts[e].cuts[m]*(v1[0]-v0[0]) - w[0];
		t[1] = v0[1]+qc->edge_cuts[e].cuts[m]*(v1[1]-v0[1]) - w[1];
		t[2] = v0[2]+qc->edge_cuts[e].cuts[m]*(v1[2]-v0[2]) - w[2];
		if (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]
						< Opts.EPS * Opts.EPS)
		    break;
	    }
	    if (j < qc->face_cuts[f].n)
		continue;

	    if (n != i)
		qc->face_cuts[f].info[n] = qc->face_cuts[f].info[i];
	    /* Compute the tangent vector of the interface curve.
	     *
	     * Let \u and \w be two independent vectors in the plane, then
	     * the interface curve in the plane can be parameterized as:
	     * 	\x = \x0 + a(t)*\u + b(t)*\w.
	     * Since dL(\x)/dt == 0, we have:
	     * 	a'(t)*(\grad L\cdot\u) + b'(t)(\grad L\cdot\w)
	     * thus the tangent vector of the curve is given by:
	     *	\x'(t) = a'(t)\u + b'(t)\w,
	     * which is parallel to:
	     * 	-(\grad L\cdot\u)\w + (\grad L\cdot\w)\u,
	     */
	    v0 = qc->verts[GetFaceVertex(f,0)];
	    v1 = qc->verts[GetFaceVertex(f,1)];
	    v2 = qc->verts[GetFaceVertex(f,2)];
	    levelset_grad(w[0], w[1], w[2], t, FALSE);
	    r = t[0]*(v1[0]-v0[0])+t[1]*(v1[1]-v0[1])+t[2]*(v1[2]-v0[2]);
	    s = t[0]*(v2[0]-v0[0])+t[1]*(v2[1]-v0[1])+t[2]*(v2[2]-v0[2]);
	    t[0] = -r * (v2[0]-v0[0]) + s * (v1[0]-v0[0]);
	    t[1] = -r * (v2[1]-v0[1]) + s * (v1[1]-v0[1]);
	    t[2] = -r * (v2[2]-v0[2]) + s * (v1[2]-v0[2]);
	    d = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
	    if (d < Opts.EPS * Opts.EPS) {
		/* The face is tangential to the interface */
		/*LONGJMP;*/
		continue;
	    }
	    d = 1.0 / Sqrt(d);
	    t[0] *= d;
	    t[1] *= d;
	    t[2] *= d;
#if 0
# warning debugging code!
	    fprintf(stderr, "face/edge %d/%d, t=[%lg %lg %lg]-[%lg %lg %lg]\n",
			f, qc->face_cuts[f].info[i] % NEdge, w[0], w[1], w[2],
			w[0] + t[0], w[1] + t[1], w[2] + t[2]);
#endif
	    n++;
	}
	if (/*f < 3 &&*/ n > 2)
	    LONGJMP;		/* at most 2 cuts on a face */
	qc->face_cuts[f].n = n;
    }
}
#undef GetFaceEdgeNo

static void
add_w_interv(FLOAT a, FLOAT b)
{
    QUAD_CONTEXT *qc = &QI.ctx;
    int n_max = (int)(sizeof(qc->interv) / sizeof(qc->interv[0])) / 2;
    FLOAT interv[n_max * 2], a1, b1;
    int i, n;

    if (Fabs(a - b) < Opts.EPS)
#if 0
	return;		/* skip 0-length interval */
#else
	a = b;		/* keep 0-length interval */
#endif

    if (a > b) {
	fprintf(stderr, "%s:%d: unexpected error!\n", __FILE__, __LINE__);
	exit(1);
    }

    for (i = n = 0; i < qc->interv_n; i++) {
	a1 = qc->interv[2 * i + 0];
	b1 = qc->interv[2 * i + 1];
	if (b1 < a - Opts.EPS || a1 > b + Opts.EPS) {
	    /* no overlap, save the interval */
	    assert(n < n_max - 1);
	    interv[2 * n + 0] = a1;
	    interv[2 * n + 1] = b1;
	    n++;
	    continue;
	}
	/* (a,b) := (a,b) \cup (a1,b1) */
	if (a > a1)
	    a = a1;
	if (b < b1)
	    b = b1;
    }
    if (n > 0)
	memcpy(qc->interv + 2, interv, 2 * n * sizeof(interv[0]));
    qc->interv[0] = a;
    qc->interv[1] = b;
    qc->interv_n = n + 1;
}

static BOOLEAN
add_u_interv(const FLOAT *pa, const FLOAT *pb)
/* Converts the interval (*pa, *pb) for u=sqrt(1-w^2)/w to intervals for w,
 * and add them to the list (makes union of them).
 *
 * Returns true if the result = {[-1,1]}, false otherwise.
 *
 * Note that *pb < *pa and pa==NULL => *pa=-infty, *pb==NULL ==> *pb=infty */
{
    QUAD_CONTEXT *qc = &QI.ctx;
    FLOAT a, b;

#ifdef S
# undef S
#endif
#define S(a)	(1.0 / Sqrt(1.0 + (a) * (a)))

    if (pa == NULL && pb == NULL) {
	add_w_interv(0.0, 0.0);
    }
    else if (pa == NULL) {
	b = *pb;
	/* sqrt{1-w^2}/w <= b */
	if (b >= 0.0) {
	    add_w_interv(-1.0, 0.0);
	    add_w_interv(S(b), 1.0);
	}
	else {
	    add_w_interv(-S(b), 0.0);
	}
    }
    else if (pb == NULL) {
	a = *pa;
	/* sqrt{1-w^2}/w >= a */
	if (a >= 0.0) {
	    add_w_interv(0.0, S(a));
	}
	else {
	    add_w_interv(0.0, 1.0);
	    add_w_interv(-1.0, -S(a));
	}
    }
    else {
	a = *pa;
	b = *pb;
	assert(a <= b);
	if (a >= 0.0) {
	    assert(b >= 0.0);
	    add_w_interv(S(b), S(a));
	}
	else if (b <= 0.0) {
	    assert(a < 0.0);
	    add_w_interv(-S(b), -S(a));
	}
	else {
	    assert(a < 0.0 && b > 0.0);
	    add_w_interv(-1.0, -S(a));
	    add_w_interv(S(b), 1.0);
	}
    }

#undef S

    return (qc->interv_n == 1 &&
	qc->interv[0] < -1.0 + Opts.EPS && qc->interv[1] > 1.0 - Opts.EPS);
}

static void
set_candidate_intervals(const FLOAT *v0, const FLOAT *v1, const FLOAT *v2)
/* constructs intervals for w such that w*(v1-v0)+sqrt(1-w^2)*(v2-v0) may be
 * an admissible direction for nw, i.e., the intersection line of the
 * integration plane (nu,nv) with a face is never tangential to the
 * intersection curve of the interface with that face.
 *
 * The intervals are stored in qc->intervals.
 */
{
    QUAD_CONTEXT *qc = &QI.ctx;
    const FLOAT *t1, *t2, *nf;
    COORD tmp;
    FLOAT a1, b1, a2, b2, a, b, c, d, r1, r2;
    int i, n, n_max = (int)(sizeof(qc->interv) / sizeof(qc->interv[0])) / 2;
    FLOAT interv[n_max * 2];

    /* First, find the inadmissible intervals on all faces */

    qc->interv_n = 0;
    for (i = 0; i < NFace; i++) {
	if (qc->face_cuts[i].n <= 1)
	    continue;
	assert(qc->face_cuts[i].n == 2);

	nf = qc->normal[i];	/* unit normal vector of the face f */
	t1 = qc->face_cuts[i].tv[0];
	t2 = qc->face_cuts[i].tv[1];
 
	/* tmp := (v1-v0)\cross qc->normal[i] */
	tmp[0] = (v1[1] - v0[1]) * nf[2] - (v1[2] - v0[2]) * nf[1];
	tmp[1] = (v1[2] - v0[2]) * nf[0] - (v1[0] - v0[0]) * nf[2];
	tmp[2] = (v1[0] - v0[0]) * nf[1] - (v1[1] - v0[1]) * nf[0];
	a1 = (tmp[1] * t1[2] - tmp[2] * t1[1]) * nf[0] +
	     (tmp[2] * t1[0] - tmp[0] * t1[2]) * nf[1] +
	     (tmp[0] * t1[1] - tmp[1] * t1[0]) * nf[2];
	a2 = (tmp[1] * t2[2] - tmp[2] * t2[1]) * nf[0] +
	     (tmp[2] * t2[0] - tmp[0] * t2[2]) * nf[1] +
	     (tmp[0] * t2[1] - tmp[1] * t2[0]) * nf[2];

	/* tmp := (v2-v0)\cross qc->normal[i] */
	tmp[0] = (v2[1] - v0[1]) * nf[2] - (v2[2] - v0[2]) * nf[1];
	tmp[1] = (v2[2] - v0[2]) * nf[0] - (v2[0] - v0[0]) * nf[2];
	tmp[2] = (v2[0] - v0[0]) * nf[1] - (v2[1] - v0[1]) * nf[0];
	b1 = (tmp[1] * t1[2] - tmp[2] * t1[1]) * nf[0] +
	     (tmp[2] * t1[0] - tmp[0] * t1[2]) * nf[1] +
	     (tmp[0] * t1[1] - tmp[1] * t1[0]) * nf[2];
	b2 = (tmp[1] * t2[2] - tmp[2] * t2[1]) * nf[0] +
	     (tmp[2] * t2[0] - tmp[0] * t2[2]) * nf[1] +
	     (tmp[0] * t2[1] - tmp[1] * t2[0]) * nf[2];

	/* transform into: a * u^2 + b * u + c */
	a = b1 * b2;
	b = a1 * b2 + a2 * b1;
	c = a1 * a2;

#if 0
	/* unnecessary, always 2 real roots unless a==0 */
	if (b * b <= 4.0 * a * c + Opts.EPS * Opts.EPS) {
	    if (c >= -Opts.EPS) {
		continue;	/* empty */
	    }
	    else {
		/* [-\infty, \infty] */
		qc->interv_n = 1;
		qc->interv[0] = -1.0;
		qc->interv[1] = 1.0;
		break;
	    }
	}
#endif

	if (a == 0.0) {
	    /* linear */
	    if (b < 0.0) {
		/* [-c/b, infty] */
		a = -c / b;
		if (add_u_interv(&a, NULL) == TRUE)
		    break;
	    }
	    else if (b > 0.0) {
		/* [-infty, -c/b] */
		b = -c / b;
		if (add_u_interv(NULL, &b) == TRUE)
		    break;
	    }
	    else {
		if (Fabs(c) > Opts.EPS) {	/* c != 0.0 */
		    /* Note: this case may occur with a flat interface, as with
		     * the example below:
		     *	make USER_CFLAGS="-DCASE=1 -DLS_ORDER=2" -s quad_test2
		     *	./quad_test2 -mesh_file cube.dat -order 2
		     * In this case we should have c != 0, and the inadmissible
		     * interval [0,0] for w. */
		    if (add_u_interv(NULL, NULL) == TRUE)
			break;
		}
		else if (Fabs(t1[0] * t2[0] + t1[1] * t2[1] + t1[2] * t2[2]
				+ 1.0) <= FLOAT_EPSILON) {
		    /* The angle between the 2 tangent vectors is \pi, as in
		     * the following example (from cube5.dat):
		     *	tet: [1,1,0],[0,1,1],[0,0,0],[0,1,0] (from cube5.dat),
		     *	interface: z-cylinder, radius .3, center (.5,.5,.5).
		     */
		    LONGJMP;
		}
		else {
		    output_vtk3("tmp.vtk", qc->verts[0], qc->verts[1],
					   qc->verts[2], qc->verts[3]);
		    phgError(1, "%s:%d: unexpected case on face %d "
				"(-qi_dbg_vtk to dump the element, which may "
				"be too large to resolve the interface).\n",
				__FILE__, __LINE__, i);
		}
	    }
	    continue;
	}

	/* (r1,r2) = the two roots of the quadratic equation */
#if 0
	d = Sqrt(Fabs(b * b - 4.0 * a * c));
	a += a;
	if (b >= 0.0) {
	    r1 = -(c + c) / (b + d);	/* -2c/(b+Delta) */
	    r2 = -(b + d) / a;		/* -(b+Delta)/(2a) */

	}
	else {
	    r1 = (-b + d) / a;		/* (-b+Delta)/(2a) */
	    r2 = (c + c) / (-b + d);	/* 2c/(-b+Delta) */
	}
#else
	r1 = -a1 / b1;
	r2 = -a2 / b2;
#endif

	/* reorder such that r1 <= r2 */
	if (r1 > r2) {
	    d = r1; r1 = r2; r2 = d;
	}

	if (a > 0.0) {
	    /* [r1, r2] */
	    if (add_u_interv(&r1, &r2) == TRUE)
		break;
	}
	else {
	    /* [-infty, r1] and [r2, infty] */
	    if (add_u_interv(NULL, &r1) == TRUE)
		break;
	    if (add_u_interv(&r2, NULL) == TRUE)
		break;
	}
	
    }

    /* Now construct candidate intervals := [-1,1]\setminus{inadmissible} */

    if (qc->interv_n == 0) {
	/* the entire interval [-1,1] */
	qc->interv_n = 1;
	qc->interv[0] = -1.0;
	qc->interv[1] = 1.0;
	return;
    }

    memcpy(interv, qc->interv, 2 * qc->interv_n * sizeof(*interv));

    /* sort the intervals */
    for (i = 0; i < qc->interv_n - 1; i++) {
	for (n = i + 1; n < qc->interv_n; n++) {
	    if (interv[2 * i + 1] <= interv[2 * n + 0])
		continue;
	    d = interv[2 * i + 0];
	    interv[2 * i + 0] = interv[2 * n + 0];
	    interv[2 * n + 0] = d;

	    d = interv[2 * i + 1];
	    interv[2 * i + 1] = interv[2 * n + 1];
	    interv[2 * n + 1] = d;
	}
    }

    n = qc->interv_n;
    qc->interv_n = 0;

    if (interv[0] > -1.0 + Opts.EPS) {
	qc->interv[0] = -1.0;
	qc->interv[1] = interv[0];
	qc->interv_n++;
    }

    assert(n - 1 + qc->interv_n <= n_max);
    for (i = 0; i < n - 1; i++) {
	qc->interv[2 * qc->interv_n + 0] = interv[2 * i + 1];
	qc->interv[2 * qc->interv_n + 1] = interv[2 * (i + 1) + 0];
	qc->interv_n++;
    }

    if (interv[2 * i + 1] < 1.0 - Opts.EPS) {
	assert(qc->interv_n < n_max);
	qc->interv[2 * qc->interv_n + 0] = interv[2 * i + 1];
	qc->interv[2 * qc->interv_n + 1] = 1.0;
	qc->interv_n++;
    }

    /* merge [-1,a] and [b,1] into [b-2, a] */
    if (qc->interv_n > 1 && qc->interv[0] == -1.0
	&& qc->interv[2 * qc->interv_n - 1] == 1.0) {
	qc->interv[0] = qc->interv[2 * qc->interv_n - 2] - 2.0;
	qc->interv_n--;
    }

    if (!Opts.dbg_off && Opts.show_intervals_w) {
	fprintf(stderr, "candidate intervals:%s\n",
		qc->interv_n == 0 ? " none" : "");
	for (i = 0; i < qc->interv_n; i++)
	    fprintf(stderr, "    [%lg,%lg]\n",
		    (double)qc->interv[2 * i], (double)qc->interv[2 * i + 1]);
    }
}

static FLOAT
set_dir(const FLOAT *v0, const FLOAT *v1, const FLOAT *v2, FLOAT w,
	FLOAT *nw)
/* This function sets the integration direction nw. It computes:
 *	nw := w*(v1-v0) + \sqrt{1-w^2}*(v2-v0),
 * and returns a(w), which is a number in [0,1] indicating how the intersection
 * line of the integration plane (nu,nv) with a face is tangential to the
 * intersection curve of the interface with the face, at the intersection
 * points (the closer to 1 the returned value, the more tangential the
 * intersection line to the curve).
 *
 * More precisely, the function returns 1.0 if w is not in the candidate
 * intervals. Otherwise it returns:
 *	max_{face in {0,1,2,(3?)}, t\in T(face)} \{|(v,t)|\},
 * where:
 * 	v := 'nw'\cross normal(face), normalized
 *	T(face) := the set of unit tangent vectors of the interface curve
 *		   on the face 'face' at the cut points with the edges.
 *
 * Note this function is extended periodically from [-1,1] to R.
 */
{
    QUAD_CONTEXT *qc = &QI.ctx;
    const FLOAT *t;
    COORD v;
    FLOAT ret = 0.0, d;
    int i, f;

    /* translate w into the interval [-1,1] */
    assert(w >= -3.0 && w <= 3.0);
    if (w < -1.0)
	w += 2.0;
    else if (w > 1.0)
	w -= 2.0;

    /* construct vector nw */
    d = Sqrt(1.0 - w * w);
    nw[0] = w * (v1[0] - v0[0]) + d * (v2[0] - v0[0]);
    nw[1] = w * (v1[1] - v0[1]) + d * (v2[1] - v0[1]);
    nw[2] = w * (v1[2] - v0[2]) + d * (v2[2] - v0[2]);

    /* check whether w is in the candidate intervals and return 1 if not */
    for (i = 0; i < qc->interv_n; i++)
	if ((w >= qc->interv[2 * i] && w <= qc->interv[2 * i + 1]) ||
	    (w - 2.0 >= qc->interv[2 * i] && w - 2.0 <= qc->interv[2 * i + 1]))
	    break;
    if (i >= qc->interv_n)
	return 1.0;

    for (f = 0; f < (Opts.use_face_normal == FALSE ? NFace : 3); f++) {
	int n = qc->face_cuts[f].n;

	if (n == 0)
	    continue;

	/* v:= 'nw'\cross normal(f), normalized */
	v[0] = nw[1] * qc->normal[f][2] - nw[2] * qc->normal[f][1];
	v[1] = nw[2] * qc->normal[f][0] - nw[0] * qc->normal[f][2];
	v[2] = nw[0] * qc->normal[f][1] - nw[1] * qc->normal[f][0];
	d = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
	if (d < Opts.EPS * Opts.EPS)
	    continue;		/* 'nw' orthogonal to face(f) */
	d = 1.0 / Sqrt(d);
	v[0] *= d;
	v[1] *= d;
	v[2] *= d;
	for (i = n = 0; i < qc->face_cuts[f].n; i++) {
	    t = qc->face_cuts[f].tv[i];	/* tangent vector in the face */
	    /* d := Fabs(cos(\theta)) = |(t, v)| */
	    d = Fabs(v[0] * t[0] + v[1] * t[1] + v[2] * t[2]);
	    if (ret < d) {
		ret = d;
		if (ret > Opts.threshold)
		    break;
	    }
	}
    }

    return ret;
}

static FLOAT
search_618(const FLOAT *v0, const FLOAT *v1, const FLOAT *v2,
	   FLOAT a, FLOAT b, int n, FLOAT *nw, FLOAT *val)
/* returns (w,nw) := argmin_{w\in(a,b)} set_dir(v0,v1,v2,w,nw) using
 * the 0.618 method. n is the number of iterations. */
{
    FLOAT w1, f1, w2, f2;
    /* 0.618 = (Sqrt(5.0) - 1.0) * 0.5 */
    static const FLOAT d618 = _F(0.61803398874989484820458683436563812);
    w1 = a + (1.0 - d618) * (b - a);
    f1 = set_dir(v0, v1, v2, w1, nw);
    while (n-- > 0) {
	COORD nw2;
	w2 = a + d618 * (b - a);
	f2 = set_dir(v0, v1, v2, w2, nw2);
	if (f1 < f2) {
	    /* [a,w2] */
	    b = a;
	    a = w2;
	}
	else {
	    /* [w1,b] */
	    a = w1;
	    w1 = w2;
	    f1 = f2;
	    memcpy(nw, nw2, sizeof(nw2));
	}
    }

    *val = f1;
    return w1;
}

static void
setup_xwv(FLOAT *x0, FLOAT *nw, FLOAT *nv, FLOAT *nu)
/* sets up the integration coordinates x0 + t * nw + s * nv + r * nu */
{
    QUAD_CONTEXT *qc = &QI.ctx;
    FILE *fvtk;
    FLOAT a, b, d;

    qc->hw = Sqrt(nw[0] * nw[0] + nw[1] * nw[1] + nw[2] * nw[2]);

    /* 'nv' = nu\cross nw */
    nv[0] = nu[1] * nw[2] - nu[2] * nw[1];
    nv[1] = nu[2] * nw[0] - nu[0] * nw[2];
    nv[2] = nu[0] * nw[1] - nu[1] * nw[0];
    qc->hv = Sqrt(nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]);

    /* use v0 as initial x0 */
    x0[0] = qc->verts[0][0];
    x0[1] = qc->verts[0][1];
    x0[2] = qc->verts[0][2];

    qc->nv = nv;
    qc->nw = nw;
    qc->x0 = x0;

    /* adjust the BoundingBox and x0 */
#ifdef BB1
# undef BB1
#endif
#define BB1(v, n, h2) \
    d = ProjectionOnVector(v, x0, n, h2); \
    if (a > d) \
	a = d; \
    if (b < d) \
	b = d; 

#ifdef BB
# undef BB
#endif
#define BB(n, h2, h) \
    a = b = 0.0; \
    BB1(qc->verts[3], n, h2) \
    BB1(qc->verts[2], n, h2) \
    BB1(qc->verts[1], n, h2) \
    if (a != 0.0) { \
	x0[0] += a * n[0]; \
	x0[1] += a * n[1]; \
	x0[2] += a * n[2]; \
    } \
    if (b - a != 1.0) { \
	n[0] *= b - a; \
	n[1] *= b - a; \
	n[2] *= b - a; \
	h *= b - a; \
    }

    BB(nu, qc->hu * qc->hu, qc->hu);
    BB(nv, qc->hv * qc->hv, qc->hv);
    BB(nw, qc->hw * qc->hw, qc->hw);

#undef BB1
#undef BB

    if (Opts.dbg_off || !Opts.show_directions)
	return;

#ifdef C1
# undef C1
#endif
#define C1(a) (double)a[0], (double)a[1], (double)a[2]

#ifdef C2
# undef C2
#endif
#define C2(a,b) (double)(a[0]+b[0]), (double)(a[1]+b[1]), (double)(a[2]+b[2])

#ifdef C3
# undef C3
#endif
#define C3(a,b,c) \
    (double)(a[0]+b[0]+c[0]), (double)(a[1]+b[1]+c[1]), (double)(a[2]+b[2]+c[2])

#ifdef C4
# undef C4
#endif
#define C4(a,b,c,d) (double)(a[0]+b[0]+c[0]+d[0]), \
    (double)(a[1]+b[1]+c[1]+d[1]), (double)(a[2]+b[2]+c[2]+d[2])

    /* print x0, nu, nv, nw */
    fprintf(stderr,
		"x0=[%lg %lg %lg]\nx0+nu=[%lg %lg %lg]\n"
		"x0+nv=[%lg %lg %lg]\nx0+nw=[%lg %lg %lg]\n",
		C1(x0), C2(x0,nu), C2(x0,nv), C2(x0,nw));

    if (Opts.dbg_off || !Opts.dbg_vtk)
	return;

    /* create a VTK file for the integration domain */
#if 0
    fvtk = fopen("box.vtk", "wt");
    /* type 12 for VTK_HEXAHEDRON, the vertices are ordered as:
     *		lower 4 then upper 4, counter-clockwisely */
    fprintf(fvtk,
		"# vtk DataFile Version 2.0\n"
		"phgQuadInterface integration box\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n"
		"POINTS 8 float\n"
		"%0.12lg %0.12lg %0.12lg\n%0.12lg %0.12lg %0.12lg\n"
		"%0.12lg %0.12lg %0.12lg\n%0.12lg %0.12lg %0.12lg\n"
		"%0.12lg %0.12lg %0.12lg\n%0.12lg %0.12lg %0.12lg\n"
		"%0.12lg %0.12lg %0.12lg\n%0.12lg %0.12lg %0.12lg\n"
		"CELLS 1 9\n"
		"8 0 1 2 3 4 5 6 7\n"
		"CELL_TYPES 1\n"
		"12\n",
		C1(x0),    C2(x0,nw),    C3(x0,nw,nv),    C2(x0,nv),
		C2(x0,nu), C3(x0,nw,nu), C4(x0,nw,nv,nu), C3(x0,nv,nu));
    fclose(fvtk);
    fprintf(stderr, "*** VTK file \"box.vtk\" created.\n");
#else
    /* Type 3 for VTK_LINE */
    fvtk = fopen("dirs.vtk", "wt");
    fprintf(fvtk,
		"# vtk DataFile Version 2.0\n"
		"phgQuadInterface integration box\n"
		"ASCII\n"
		"DATASET POLYDATA\n"
		"POINTS 4 float\n"
		"%0.12lg %0.12lg %0.12lg\n%0.12lg %0.12lg %0.12lg\n"
		"%0.12lg %0.12lg %0.12lg\n%0.12lg %0.12lg %0.12lg\n"
		"LINES 3 9\n"
		"2 0 1\n"
		"2 0 2\n"
		"2 0 3\n"
		"CELL_DATA 3\n"
		"SCALARS cell_id unsigned_char\n"
		"LOOKUP_TABLE default\n"
		"0 1 2\n",
		C1(x0), C2(x0,nu), C2(x0,nv), C2(x0,nw));
    fclose(fvtk);
    fprintf(stderr, "*** VTK file \"dirs.vtk\" created"
				" (blue:nu, white:nv, red:nw).\n");
#endif    

#undef C1
#undef C2
#undef C3
#undef C4
}

static void
quad_sub(TET *pt, FLOAT res[3])
{
    int i, j;
    const FLOAT *v0, *v1, *v2, *v3;

    static FLOAT d0, d, a, b, w, w0;
    static COORD g, tmp, x0, nu = {0.,0.,0.}, nv, nw, n0, n1, n2, n3;
    static QUAD_CONTEXT *qc = NULL;
#if USE_OMP
# pragma omp threadprivate(d0, d, a, b, w, w0, g, tmp, x0, nu, nv, nw, n0, n1, n2, n3, qc)
#endif  /* USE_OMP */

#if 0
# warning debugging code!
/* Setup debugging flags for tracking a specific element */
static COORD verts[4] = {{0.5,0.25,0.625},
			 {0.5,0.3125,0.5625},
			 {0.5625,0.3125,0.625},
			 {0.53125,0.28125,0.5625}};
# if 0
j = !memcmp(verts, pt->verts, sizeof(verts));
# else
for (i = 0; i < NVert; i++) if (Fabs(verts[i][0]-pt->verts[i][0]) +
				Fabs(verts[i][1]-pt->verts[i][1]) +
				Fabs(verts[i][2]-pt->verts[i][2]) > Opts.EPS) break;
j = (i >= NVert);
# endif
Opts.show_directions = Opts.show_intervals_w = j;
#endif

    if (qc == NULL)
	qc = &QI.ctx;

    qc->pt = pt;

    v0 = pt->verts[0];
    v1 = pt->verts[1];
    v2 = pt->verts[2];
    v3 = pt->verts[3];

    phgInfo(2, "Tetra: [%lg %lg %lg][%lg %lg %lg][%lg %lg %lg][%lg %lg %lg]\n",
		v0[0],v0[1],v0[2], v1[0],v1[1],v1[2],
		v2[0],v2[1],v2[2], v3[0],v3[1],v3[2]);

    /* rearrange (v0,v1,v2,v3) such that the face (v0,v1,v2) is "the most
     * parallel" to the interface and store its normal direction in nu */

    /* get normal direction of the interface at the barycenter */
    levelset_grad((v0[0] + v1[0] + v2[0] + v3[0]) * 0.25,
		  (v0[1] + v1[1] + v2[1] + v3[1]) * 0.25,
		  (v0[2] + v1[2] + v2[2] + v3[2]) * 0.25, g, FALSE);
    d0 = Sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);	/* d0 = |grad_L| */
    if (d0 < /*Sqrt*/(Opts.EPS)) {
	levelset_func((v0[0] + v1[0] + v2[0] + v3[0]) * 0.25,
		      (v0[1] + v1[1] + v2[1] + v3[1]) * 0.25,
		      (v0[2] + v1[2] + v2[2] + v3[2]) * 0.25, &d0);
	if (Fabs(d0) < Opts.EPS)
	    phgWarning("%s:%d: degenerate interface.\n", __FILE__, __LINE__);
	if (!Opts.dbg_off && Opts.dbg_vtk)
	    output_vtk("tmp.vtk", pt->verts);
	LONGJMP;
    }
    /* normalize grad_L */
    d = 1.0 / d0;
    g[0] *= d;
    g[1] *= d;
    g[2] *= d;

    d = -1.0;
    j = -1;
    for (i = 0; i < NFace; i++) {
	/* compute normal of the face i */
	v0 = pt->verts[GetFaceVertex(i, 0)];
	v1 = pt->verts[GetFaceVertex(i, 1)];
	v2 = pt->verts[GetFaceVertex(i, 2)];
	get_face_normal(tmp, v0, v1, v2, NULL);
	a = Fabs(tmp[0] * g[0] + tmp[1] * g[1] + tmp[2] * g[2]);
	if (d < a) {
	    d = a;
	    j = i;
	    memcpy(nu, tmp, sizeof(nu));
	}
    }

    if (Opts.use_face_normal == FALSE)
	memcpy(nu, g, sizeof(nu));	/* set nu to the interface's normal */
    else
	phgInfo(2, "Using face %d for nu, cos(angle) = %lg\n", j, (double)d);

    /*if (QI.quad_type == 0)*/ {
	/* Surface integral: check that \grad L\cdot nu has constant sign and
	 * |\grad L|/|\grad L\cdot nu| <= 4.0 */
	if (Fabs(g[0] * nu[0] + g[1] * nu[1] + g[2] * nu[2]) < 0.25)
	    LONGJMP;
#if 0
	/* check at the 4 vertices */
	for (i = 0; i < NVert; i++) {
	    levelset_grad(pt->verts[i][0], pt->verts[i][1], pt->verts[i][2],
			  g, FALSE);
	    d = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
	    if (d <= Opts.EPS /** Opts.EPS*/)
		LONGJMP;
	    d = Sqrt(d);
	    a = g[0] * nu[0] + g[1] * nu[1] + g[2] * nu[2];
	    if (a == 0.0 || d > 4.0 * Fabs(a))
		LONGJMP;
	}
#endif
    }

    qc->verts[0] = v0 = pt->verts[GetFaceVertex(j, 0)];
    qc->verts[1] = v1 = pt->verts[GetFaceVertex(j, 1)];
    qc->verts[2] = v2 = pt->verts[GetFaceVertex(j, 2)];
    qc->verts[3] = v3 = pt->verts[j];

    compute_Jacobian();

    /* compute height in the 'nu' direction, which is equal to:
     *		max((v3 - v[012]).nu) - min((v3 - v[012]).nu) == b - a */
    b = a = 0.0;
    for (i = 0; i < 3; i++) {
	d = (v3[0] - qc->verts[i][0]) * nu[0] +
	    (v3[1] - qc->verts[i][1]) * nu[1] +
	    (v3[2] - qc->verts[i][2]) * nu[2];
	if (a > d)
	    a = d;
	else if (b < d)
	    b = d;
	if (Opts.use_face_normal != FALSE)
	    break;
    }
    d = b - a;
#if 0
    if (Opts.use_face_normal != FALSE)	/* for the old behavior */
	d = (v3[0]-v0[0])*nu[0] + (v3[1]-v0[1])*nu[1] + (v3[2]-v0[2])*nu[2];
#endif
    nu[0] *= d;
    nu[1] *= d;
    nu[2] *= d;
    qc->nu = nu;
    qc->hu = Fabs(d);

    /* n0,n1,n2,n3 := the normal vectors of the faces F0, F1, F2, F3 */
    get_face_normal(n0, v1, v2, v3, NULL);
    get_face_normal(n1, v0, v2, v3, NULL);
    get_face_normal(n2, v0, v1, v3, NULL);
    if (Opts.use_face_normal != FALSE) {
	n3[0] = nu[0] / qc->hu;
	n3[1] = nu[1] / qc->hu;
	n3[2] = nu[2] / qc->hu;
    }
    else {
	get_face_normal(n3, v0, v1, v2, NULL);
    }

    qc->normal[0] = n0;
    qc->normal[1] = n1;
    qc->normal[2] = n2;
    qc->normal[3] = n3;

    /* compute and store cut-points on all edges/faces */
    get_face_edge_cuts();

    if (Opts.use_face_normal == FALSE) {
	/* cumpute (w1,w2) such that w[12]-v0 are orthogonal to nu */
	static COORD w1, w2;
#if USE_OMP
# pragma omp threadprivate(w1, w2)
#endif  /* USE_OMP */

	/* w1 is obtained from nu by setting the coordinate with smallest
	 * absolute value to 0, and swapping the other two coordinates with
	 * one of them negated. */
	a = Fabs(nu[0]);
	i = 1;
	j = 2;
	if (a > (b = Fabs(nu[1]))) {
	    a = b;
	    i = 0;
	    j = 2;
	}
	if (a > (b = Fabs(nu[2]))) {
	    i = 0;
	    j = 1;
	}
	w1[3 - i - j] = 0.0;
	w1[i] = -nu[j];
	w1[j] = nu[i];

	/* w2 = nu \cross w1 */
	w2[0] = nu[1] * w1[2] - nu[2] * w1[1];
	w2[1] = nu[2] * w1[0] - nu[0] * w1[2];
	w2[2] = nu[0] * w1[1] - nu[1] * w1[0];

	w1[0] += v0[0];
	w1[1] += v0[1];
	w1[2] += v0[2];
	v1 = w1;

	w2[0] += v0[0];
	w2[1] += v0[1];
	w2[2] += v0[2];
	v2 = w2;
    }

    /* get candidate intervals for w */
    set_candidate_intervals(v0, v1, v2);

    if (!Opts.dbg_off && Opts.dbg_w >= -1.0 && Opts.dbg_w <= 1.0) {
	/* use specified w value */
	a = set_dir(v0, v1, v2, w0 = Opts.dbg_w, nw);
    }
    else if (!Opts.dbg_off && Opts.dbg_w > 1.0) {
	/* use a random w value */
	/*srand48((long int)(Opts.dbg_w));*/
	a = set_dir(v0, v1, v2, w0 = 2.0 * drand48() - 1.0, nw);
    }
    else if (qc->interv_n == 0) {
	/* FIXME: better direction? */
	a = set_dir(v0, v1, v2, w0 = -1.0, nw);
    }
    else {
	/* find vector nw:=nw(w) with w := argmin_{w\in[-1,1]} set_dir(w). */
	w0 = -2.0;
	a = Opts.threshold + 1.0;
	for (j = 0; j < qc->interv_n; j++) {
	    if (!Opts.dbg_off && Opts.dbg_wN <= 0) {
		/* middle point of the interval */
		w = (qc->interv[2 * j] + qc->interv[2 * j + 1]) * 0.5;
		b = set_dir(v0, v1, v2, w, nv);
	    }
	    else {
		/* 0.618 search */
		w = search_618(v0, v1, v2, qc->interv[2*j], qc->interv[2*j+1],
					Opts.dbg_wN, nv, &b);
	    }
	    if (b < a /*&& b <= Opts.threshold*/) {
		w0 = w;
		a = b;
		memcpy(nw, nv, sizeof(nv));
	    }
	}
    }

    /* Note: without the first test, some cases lead to intolerably many
     * recursions when using 128-bit FP. Below are a few of the examples
     * (with the level-set function x^2+y^2+z^2-0.75^2):
     *	1. Domain: (-1,1)x(-2,2)x(-3,3), uniform 5x5x5 grid refined 5 times
     *	   ./quad_test3-p4est -n 5 -order 3 -refine 5 -qi_cuboid2tetra
     *
     *	2. Some tetrahedra (derived from the case above):
     *		{-0.001171875, -0.44921875, 0.6},
     *		{-0.001171875, -0.449609375, 0.5994140625},
     *		{-0.00078125, -0.45, 0.6},
     *		{0.0, -0.45, 0.6}
     *
     *		{-0.0005859375, -0.449609375, 0.6},
     *		{-0.0005859375, -0.4498046875, 0.59970703125},
     *		{-0.00078125, -0.45, 0.6},
     *		{0.0, -0.45, 0.6}	
     */
    if (tet_volume(pt) > Sqrt(Opts.EPS) &&
	a > Opts.threshold && !Opts.dbg_off &&
	(Opts.dbg_w < -1.0/* || Opts.dbg_w > 1.0*/)) {
	if (phgVerbosity > 1)
	    phgInfo(-1, "a(w) (%lg) is greater than threshold (%lg), "
		"subdividing element.\n", (double)a, (double)Opts.threshold);
	LONGJMP;
    }

    if (!Opts.dbg_off && Opts.show_directions)
	fprintf(stderr, "set_dir(%lg) = %lg\n", (double)w0, (double)a);

    setup_xwv(x0, nw, nv, nu);

    if (Opts.graph_w > 0) {
	/* output "w a(w) error" (for plotting) */
	static FILE *fp = NULL;
	static char fn[32];
	static int flag;
	static FLOAT ref;
#if USE_OMP
# pragma omp threadprivate(fp, fn, flag, ref)
#endif  /* USE_OMP */

	for (i = 0; i < NFace; i++)
	    if (qc->face_cuts[i].n == 2)
		break;
	if (i >= NFace)
	    return;	/* no intersection with the interface */

	/* compute a reference value */
	i = QI.quad_order;
	QI.quad_order = 65;
	quad_start(res);
	QI.quad_order = i;
	ref = res[0];

	if (fp != NULL)
	    fclose(fp);
	sprintf(fn, "graph_w%d.out", phgRank);
	fp = fopen(fn, "wt");
	assert(fp != NULL);
	Opts.threshold = 1.0; /* set_dir always return correct value */
	flag = 0;
	QI.skip_rule = TRUE;
	for (i = 0; i <= Opts.graph_w; i++) {
	    w = (i + i - Opts.graph_w) / (FLOAT)Opts.graph_w;
	    a = set_dir(v0, v1, v2, w, nw);
	    setup_xwv(x0, nw, nv, nu);
	    quad_start(res);
	    d = Fabs(res[0] - ref);
	    fprintf(fp, "%0.16lg %0.16lg %0.16lg %le\n", (double)w, (double)a,
		   (double)(1.125 + Log(d + FLOAT_EPSILON) * 0.01), (double)d);
	    if (a != 0.0)
		flag = 1;
	}
	QI.skip_rule = FALSE;
	fclose(fp);
	fprintf(stderr, "Created \"%s\".\n", fn);
	fp = NULL;
	if (flag)
	    exit(0);
	return;
    }

    quad_start(res);
}

#undef ProjectionOnVector

static void
quad_interface0(TET *pt, FLOAT res[3])
{
    int i, rule_last[3];

    if (QI.ls == NULL ||
	(QI.ls->type == DOF_ANALYTIC && QI.ls->userfunc == NULL)) {
	res[0] = res[1] = res[2] = 0.;
	for (i = 0; i < 3; i += 2)
	    if (QI.rules[i] != NULL)
		break;
	if (i >= 3)
	    return;
	QI.quad_type = i - 1;
	quad_tetra(pt, &res[i]);	/* integral on the whole tetra */
	QI.quad_type = 99;		/* for debugging */
	return;
    }

    if ((QI.ls_order <= 1 && QI.ls_order >= 0) ||
	QI.quad_order <= 1 || tet_volume(pt) <= Opts.EPS /** Opts.EPS*/) {
	/* plane interface or very small element */
	quad_planar_interface(pt, res);
	return;
    }

    if (!Opts.dbg_off && Opts.dbg_vtk && Opts.show_recursions) {
	char s[32];
	sprintf(s, "tmp%d.vtk", depth);
	output_vtk(s, pt->verts);
    }

    if (depth > Opts.subdiv_limit) {
        fprintf(stderr, "Recursion limit (%d) exceeded (vol=%lg).\n",
			(int)Opts.subdiv_limit, (double)tet_volume(pt));
	exit(1);
    }

    Unused(output_vtk2_);
    Unused(get_quad_1D);
    Unused(get_quad_2D);

    QI.skip_rule = FALSE;
    for (i = 0; i < 3; i++)
	if (QI.rules[i] != NULL)
	    rule_last[i] = QI.rules[i]->n;	/* save position */

    if ((i = setjmp(JMP_ENV)) != 0) {
	/* This is the only entry for subdividing the element */
	FLOAT res0[3];
	TET *ppt[nChild];
	int j;

	/* restore position. Test case:
	 *	quad_test2.c, CASE 0, -mesh_file cube4.dat -refine >=12 */
	for (j = 0; j < 3; j++)
	    if (QI.rules[j] != NULL)
		QI.rules[j]->n = rule_last[j];
	depth++;

	QI.skip_rule = FALSE;

	if (!Opts.dbg_off && Opts.show_recursions)
	    fprintf(stderr, "elem %"dFMT": recursion from line %d, depth %d\n",
			QI.index, i, depth);
	tet_subdivide(pt, ppt);
	phgInfo(2, "%s: child 0\n", __func__);
	quad_interface0(ppt[0], res);
	tet_free(ppt[0]);
	for (j = 1; j < nChild; j++) {
	    phgInfo(2, "%s: child %d\n", __func__, j);
	    quad_interface0(ppt[j], res0);
	    for (i = 0; i < 3; i++)
		res[i] += res0[i];
	    tet_free(ppt[j]);
	}
	depth--;
	return;
    }

    quad_sub(pt, res);

#if 0
#warning test-only code
    if (res[0] != 0.0) {
	static FLOAT error0 = -1.0;
	FLOAT res0[QI.dim], error;
	int order_bak = QI.quad_order;
	QI.quad_order = 25;
	fprintf(stderr, "***--------------- Check with the order %d result.\n",
		QI.quad_order);
	quad_sub(pt, res0);
	QI.quad_order = order_bak;
	error = Fabs(res[0] - res0[0]);
	if (res0[0] != 0.0)
	    error /= Fabs(res0[0]);
	if (error > error0) {
	    fprintf(stderr, "*** res=%lg, rel. error=%le\n", res0[0], error);
	    error0 = error;
	    output_vtk("tmp.vtk", pt->verts);
	}
    }
#endif
}

static void
quad_recursive(TET *pt, FLOAT res[3])
/* recursively subdivides the tetrahedron to 'Opts.subdiv_level0' */
{
    FLOAT res1[3];
    int i, k;
    static int level = 0;
#if USE_OMP
# pragma omp threadprivate(level)
#endif  /* USE_OMP */

    if (level < Opts.subdiv_level0) {
	TET *ppt[nChild];
	tet_subdivide(pt, ppt);
	level++;
	quad_recursive(ppt[0], res);
	tet_free(ppt[0]);
	for (i = 1; i < nChild; i++) {
	    quad_recursive(ppt[i], res1);
	    for (k = 0; k < 3; k++)
		res[k] += res1[k];
	    tet_free(ppt[i]);
	}
	level--;
	return;
    }

    quad_interface0(pt, res);

    if (!Opts.dbg_off && Opts.dbg_elem && res[0] != 0.0) {
	static FLOAT quad_error = 0.0;
#if USE_OMP
# pragma omp threadprivate(quad_error)
#endif  /* USE_OMP */
	FLOAT res0[3];
	QI.quad_order += 2;
	quad_interface0(pt, res0);
	QI.quad_order -= 2;
	res0[0] = Fabs(1.0 - res0[0] / res[0]);
	if (res0[0] > quad_error) {
	    quad_error = res0[0];
	    fprintf(stderr, "elem %"dFMT", relative error: %lg\n",
					QI.index, (double)quad_error);
	    output_vtk("tmp.vtk", pt->verts);
	}
    }
}

static void
quad_interface(DOF *ls, DOF *ls_grad, TET *pt, int quad_order,
		FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
/* Same as phgQuadInterface, except the argument 'e' is replaced by 'pt'.
 * This function is added for implementing phgQuadInterfaceTetrahedron.
 */
{
    int i;
    FLOAT res[3] = {0., 0., 0.};
    FLOAT **rules[] = {rule_m, rule_0, rule_p};

    if (rule_m == NULL && rule_0 == NULL && rule_p == NULL)
	return;		/* nothing to do */

    QI.ls = ls;
    QI.ls_grad = ls_grad;
    QI.quad_order = quad_order;
    assert(quad_order >= 0);

    for (i = 0; i < 3; i++) {
	int np = (quad_order - 1) / 2 + 1; /* # points of 1D Gaussian rule */
	if (rules[i] == NULL) {
	    QI.rules[i] = NULL;	/* ignore rule i */
	    continue;
	}
	QI.rules[i] = phgCalloc(1, sizeof(*QI.rules[i]));
	QI.rules[i]->n = QI.rules[i]->alloc = 0;
	QI.rules[i]->trunc = 10 * np * np * (i == GM ? 1 : 2 * np);
	QI.rules[i]->pw = phgCalloc(RULE_HEADER, sizeof(FLOAT));
	QI.rules[i]->nv = NULL;
    }

    /* QI.H := reference length = shortest edge */
    QI.H = 0.0;
    for (i = 0; i < NEdge; i++) {
	FLOAT d;
	const FLOAT *v0, *v1;
	v0 = pt->verts[GetEdgeVertex(i, 0)];
	v1 = pt->verts[GetEdgeVertex(i, 1)];
	d = Sqrt(Pow(v1[0]-v0[0],2) + Pow(v1[1]-v0[1],2) + Pow(v1[2]-v0[2],2));
	if (i == 0 || QI.H > d)
	    QI.H = d;
    }

    quad_recursive(pt, res);

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

void
phgQuadInterface(DOF *ls, DOF *ls_grad, ELEMENT *e, int quad_order,
		 FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
/* 
 * This function computes volume or surface integrals in subdomains of
 * a tetrahedral element cut by an interface implicitly defined by the
 * zero level set of the function 'ls'.
 *
 * 'ls' defines the level set function.
 *
 *	As a special convention, and for the convenience of users, if ls==NULL,
 *	then the integral on the whole tetrahedron is returned in *rule_m.
 *	In this case it is required ls_grad!=NULL && ls_grad->!=NULL since
 *	ls_grad->g is needed (the calling function may supply any DOF with
 *	a valid 'g' member in place of ls_grad)
 *
 * 'ls_grad' is optional, it can be set to NULL and in this case grad(ls)
 *	will be evaluated using phgDofEvalGradient().
 *
 * 'e' specifies the element, only the part of the integral within the element
 *	is computed.
 *
 * 'quad_order' specifies the order of 1D Gaussian quadrature used.
 *
 * 'res', when != NULL, points to the buffer for returning the computed result.
 *
 * 'rule_m', 'rule_0', and 'rule_p', return pointers to quadrature rules
 * 	in Omega-, Gamma, and Omega+, respectively, which are dynamically
 * 	allocated, and should be freed by the calling function after use.
 * 	The data format for the rules consists of a header followed by rule
 * 	data, and optionally followed by a list of unit normal vectors at the
 * 	quadrature points (only for rule_0 and when Opts.nv_flag is true.
 *
 * 	See the macro SetRuleHeader and the function phgQuadInterfaceRuleInfo
 * 	for actual data format in the rules.
 *
 * 	One or more of the pointers can be NULL and the corresponding rules
 * 	are not computed.
 */
{
    int i;
    GRID *g;
    TET *pt;

    if (e == NULL) {
	/* register options */
	phgOptionsRegisterTitle("\nOptions for phgQuadInterface():", "\n",
		"quad_interface");
	phgOptionsRegisterFloat("-qi_eps", "Floating point tolerance",
				&Opts.EPS);
	phgOptionsRegisterFloat("-qi_threshold",
		"Threshold on cos(angle) for nearly parallel directions",
		&Opts.threshold);
	phgOptionsRegisterNoArg("-qi_nv_flag",
		"Whether store unit normal vectors in computed rules",
		&Opts.nv_flag);

	phgOptionsRegisterKeyword("-qi_subdiv_type",
		"Subdivision method", subdiv_type_args, &Opts.subdiv_type);
	phgOptionsRegisterInt("-qi_subdiv_limit",
		"Limit on recursion depth (<=0: no subdivision)",
		&Opts.subdiv_limit);

	phgOptionsRegisterNoArg("-qi_use_face_normal",
		"Whether use one of the face normals as the direction nu",
		&Opts.use_face_normal);

	phgOptionsRegisterKeyword("-qi_adapt",
		"Adaptive quadrature order (0: off, 1: H==1, 2: h<=H^2)",
		adapt_args, &Opts.adapt);

	phgOptionsRegisterKeyword("-qi_newton",
		"Use Newton method for the nonlinear equations",
		newton_args, &Opts.newton);
	phgOptionsRegisterInt("-qi_newton_maxits",
		"Maximum number of Newton iterations", &Opts.newton_maxits);
	phgOptionsRegisterInt("-qi_newton_porder",
		"Polynomial order to get initial solutions for Newton method",
		&Opts.newton_porder);

	/* options for debugging only */
	phgOptionsRegisterTitle("\n[Development/debugging options]", "\n",
		"quad_interface");

	phgOptionsRegisterInt("-qi_subdiv_level0",
		"Initial subdivision level", &Opts.subdiv_level0);

	phgOptionsRegisterNoArg("-qi_rect2tri",
		"Split a rectangle into 2 triangles", &Opts.rect2tri);

	phgOptionsRegisterNoArg("-qi_cuboid2tetra",
		"Split a cuboid into 5 tetrahera", &Opts.cuboid2tetra);

	phgOptionsRegisterNoArg("-qi_shortcut_t",
		"Enable shortcuts in t integrals (faster, less accurate)",
		&Opts.shortcut_t);
	phgOptionsRegisterNoArg("-qi_shortcut_s",
		"Enable shortcuts in s integrals (faster, less accurate)",
		&Opts.shortcut_s);

	phgOptionsRegisterInt("-qi_graph_t",
		"Output values of g(t) on an uniform grid for plotting",
		&Opts.graph_t);
	phgOptionsRegisterInt("-qi_graph_s",
		"Output values of h(s,t) on an uniform grid for plotting",
		&Opts.graph_s);
	phgOptionsRegisterInt("-qi_graph_w", "Output set_dir() and exit",
		&Opts.graph_w);

	phgOptionsRegisterNoArg("-qi_dbg_off",
		"Disable all '-qi_dbg_*' and '-qi_show_*' options",
		&Opts.dbg_off);

	phgOptionsRegisterFloat("-qi_dbg_value_t",
		"Debug with the given t value if within [0,1]",
		&Opts.dbg_value_t);
	phgOptionsRegisterFloat("-qi_dbg_value_s",
		"Debug with the given s value if within [0,1]",
		&Opts.dbg_value_s);

	phgOptionsRegisterNoArg("-qi_dbg_roots",
		"Print information about the roots", &Opts.dbg_roots);

	phgOptionsRegisterFloat("-qi_dbg_w",
		"Set w value (<-1: auto, [-1,1]: specified, >1: random)",
		&Opts.dbg_w);
	phgOptionsRegisterInt("-qi_dbg_wN",
		"N for 0.618 search (0: use middle point)", &Opts.dbg_wN);

	phgOptionsRegisterNoArg("-qi_dbg_vtk",
		"Output VTK files \"tmp-###.vtk\" for debugging",
		&Opts.dbg_vtk);

	phgOptionsRegisterNoArg("-qi_dbg_t",
		"Estimate quadrature errors in the t-integrals", &Opts.dbg_t);
	phgOptionsRegisterNoArg("-qi_dbg_s",
		"Estimate quadrature errors in the s-integrals", &Opts.dbg_s);
	phgOptionsRegisterNoArg("-qi_dbg_elem",
		"Estimate quadrature errors in the element", &Opts.dbg_elem);

	phgOptionsRegisterNoArg("-qi_show_recursions", "Print recursions",
		&Opts.show_recursions);

	phgOptionsRegisterNoArg("-qi_show_directions",
		"Print x0, nu, nv and nw", &Opts.show_directions);
	phgOptionsRegisterNoArg("-qi_show_intervals_w",
		"Print candidate intervals for w",
		&Opts.show_intervals_w);
	phgOptionsRegisterNoArg("-qi_show_intervals_t",
		"Print t intervals", &Opts.show_intervals_t);
	phgOptionsRegisterNoArg("-qi_show_intervals_s",
		"Print s intervals", &Opts.show_intervals_s);

	phgOptionsRegisterInt("-qi_n1",
		"n for 1D composite quadrature rule", &Opts.quad_n1);
#if 0
	phgOptionsRegisterInt("-qi_n2",
		"n for 2D composite quadrature rule", &Opts.quad_n2);
#endif
	return;
    }

    assert(ls != NULL || ls_grad != NULL);
    g = (ls != NULL ? ls->g : ls_grad->g);
    assert(g != NULL);

    assert(ls != NULL || ls_grad != NULL);
#define LSOrder(ls, ls_grad, e) ((ls) != NULL ? DofTypeOrder(ls, e) : \
	 DofTypeOrder(ls_grad, e) + (DofTypeOrder(ls_grad, e) < 0 ? 0 : 1))
    QI.e = e;
    QI.index = e->index;
    QI.ls_order = LSOrder(ls, ls_grad, e);

    pt = tet_new(TRUE);
    for (i = 0; i < NVert; i++) {
	memcpy(pt->verts[i], g->verts[e->verts[i]], sizeof(FLOAT) * Dim);
	pt->is_bdry[i] = ((e->bound_type[i] & INTERIOR) == 0);
    }

    quad_interface(ls, ls_grad, pt, quad_order, rule_m, rule_0, rule_p);

    tet_free(pt);
}

void
phgQuadInterfaceTetrahedron(FUNC_3D ls, int ls_order, FUNC_3D ls_grad,
		  FLOAT const tet[Dim + 1][Dim], int quad_order,
		  FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
/* 'ls' and 'ls_grad' are function pointers for computing the level set
 * function and its gradient, 'ls_order' specifies the polynomial order of
 * the level set function (negative => non polynomial).
 *
 * 'tet' specifies the coordinates of the tetrahedron.
 *
 * The other arguments are similar as in the function phgQuadInterface(). */
{
    TET *pt;
    static DOF *dof = NULL, *grad = NULL;
#if USE_OMP
# pragma omp threadprivate(dof, grad)
#endif  /* USE_OMP */

    if (dof == NULL) {
	dof = phgCalloc(1, sizeof(*dof));
	dof->type = DOF_ANALYTIC;
	dof->dim = 1;
	grad = phgCalloc(1, sizeof(*grad));
	grad->type = DOF_ANALYTIC;
	grad->dim = Dim;

	FreeAtExitNoCheck(dof);
	FreeAtExitNoCheck(grad);
    }

    dof->userfunc = ls;
    grad->userfunc = ls_grad;

    QI.e = NULL;
    QI.index = -1;
    QI.ls_order = ls_order;

    pt = tet_new(FALSE);
    pt->verts = (void *)tet;
    quad_interface(dof, grad, pt, quad_order, rule_m, rule_0, rule_p);
    tet_free(pt);
}

/*---------------------------------------------------------------------------*/

static DOF *mark_ls, *mark_grad;
static ELEMENT *mark_e;
#if USE_OMP
#pragma omp threadprivate(mark_e)
#endif	/* USE_OMP */

static void
dummy_ls(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT xref[Dim + 1];
    phgGeomXYZ2Lambda(mark_ls->g, mark_e, x, y, z, xref);
    phgDofEval(mark_ls, mark_e, xref, value);
}

static void
dummy_grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    FLOAT xref[Dim + 1];
    phgGeomXYZ2Lambda(mark_ls->g, mark_e, x, y, z, xref);
    phgDofEval(mark_grad, mark_e, xref, values);
}

INT
phgQuadInterfaceMarkGrid(DOF *ls, DOF *ls_grad)
/* sets e->mark to:
 * 	-1 if e is in Omega- := {x: ls(x) < 0},
 *	 1 if e is in Omega+ := {x: ls(x) > 0}, and
 *	 0 if e might intersect with \Gamma := {x: ls(x) == 0}.
 * This function assumes the interface is reasonably resolved by the elements.
 *
 * Note: same (redundant) code in phg-to-p4est.c and quad-interface.c */
{
    GRID *g = ls->g;
    ELEMENT *e;
    INT n0 = 0;

    mark_ls = ls;
    if (ls_grad != NULL)
	mark_grad = ls_grad;
    else
	mark_grad = phgDofGradient(ls, NULL, NULL, NULL);

#if USE_OMP
#pragma omp parallel for private(e) reduction(+:n0)
#endif	/* USE_OMP */
    ForAllElementsBegin_(g, e, FOR_ALL)
	FLOAT (*corners)[Dim] = phgGeomGetCorners(g, e, NULL);
	mark_e = e;
	e->mark = phgQuadInterfaceMarkElement(dummy_ls, ls->poly_order,
					dummy_grad, e->elem_type, corners);
	if (e->mark == 0)
	    n0++;
    ForAllElementsEnd
 
    if (ls_grad == NULL)
	phgDofFree(&mark_grad);

    return n0;
}

/*---------------------------------------------------------------------------*/

static void
get_rule_info(const FLOAT *rule, int *npoints, int *dim, BOOLEAN *nv_flag)
/* returns the number of point in *npoints, projection type in *proj, and
 * element type (0=triangle/rectangle, 1=tetrahedron/cuboid) in *type,
 * of a quadrature rule.
 *
 * See the macro SetRuleHeader in quad-interface.h for the header format. */
{
    size_t H;
    int h;

    if (rule == NULL) {
	*npoints = *dim = 0;
	*nv_flag = FALSE;
	return;
    }

    /* Note: must match the macro SetRuleHeader */
    H = rule == NULL ? 0 : (size_t)(*rule + 0.5);
    h = (int)(H % 6);
    *npoints = (int)(H / 6);
    *dim = h / 2 + 1;
    *nv_flag = h % 2;
}

FLOAT *
phgQuadInterfaceRuleCreate(int dim, int np, const FLOAT *pts, int pinc,
			const FLOAT *wgts, int winc, FLOAT scale,
			const FLOAT *nv, int ninc)
/* returns a rule of dimension 'dim' with 'np' points. The list of points,
 * weights and normal vectors are pointed by 'pts', 'wgts' and 'nv'
 * respectively. The strides between data for consequetive two points are
 * specified by 'pinc', 'winc' and 'ninc' respectively. */
{
    int i, j;
    FLOAT *rule, *p;

    if (nv == NULL)
	rule = phgAlloc((RULE_HEADER + np * (dim + 1)) * sizeof(FLOAT));
    else
	rule = phgAlloc((RULE_HEADER + np * (2 * dim + 1)) * sizeof(FLOAT));

    p = rule + RULE_HEADER;
    for (i = 0; i < np; i++, pts += pinc, wgts += winc) {
	for (j = 0; j < dim; j++)
	    *(p++) = pts[j];
	*(p++) = *wgts * scale;
    }

    if (nv != NULL)
	for (i = 0; i < np; i++, nv += ninc)
	    for (j = 0; j < dim; j++)
		*(p++) = nv[j];

    SetRuleHeader(rule, np, dim, nv != NULL);

    return rule;
}

int
phgQuadInterfaceRuleApply(void *func_ptr, int dim, PROJ proj,
			  FLOAT *rule, FLOAT *res)
/* Computes the integral using the computed rule 'rule' and returns
 * the results in 'res'.
 *
 * Arguments:
 * 'func_ptr' -	points to the integrand function.
 *
 * 'dim'      - the dimension of the results (length of 'res'), the dimension
 *		of 'func_ptr' is assumed to match 'dim' and 'proj'.
 *
 * 'proj'     -	the projection type w.r.t. to the normal vector of the
 *		interface (only effective for surface integral)
 *
 * The return value of the function is the number of points in the rule.
 */
{
    FUNC_3D func_3D = func_ptr;
    FUNC_2D func_2D = func_ptr;
    int i, j, dimf, n, rule_dim /* dimension of the rule */;
    FLOAT *f, w;
    FLOAT *nv, *pw;

    assert(func_ptr != NULL);
    assert(proj >= 0 && proj <= PROJ_CROSS);

    n = phgQuadInterfaceRuleInfo(rule, &rule_dim, &pw, &nv);
    if (n == 0) {
	if (dim > 0)
	    memset(res, 0, dim * sizeof(*res));
	return 0;
    }
	
    if (dim <= 0 || func_ptr == NULL)
	return n;

    switch (proj) {
	case PROJ_DOT:
	    if (nv == NULL)
		phgError(1, "%s: PROJ_DOT is specified, but normal "
			    "vectors are not available.\n", __func__);
	    dimf = rule_dim * dim;
	    break;
	case PROJ_CROSS:
	    if (nv == NULL)
		phgError(1, "%s: PROJ_CROSS is specified, but normal "
			    "vectors are not available.\n", __func__);
	    if (rule_dim == 1)
		dimf = dim;
	    else if (rule_dim == 2)
		dimf = 2 * dim;
	    else
		dimf = dim;
	    assert(dimf % rule_dim == 0);
	    break;
	default:
	    dimf = dim;
	    break;
    }

    f = phgAlloc(dimf * sizeof(*f));
    memset(res, 0, dim * sizeof(*res));
    for (i = 0; i < n; i++, pw += rule_dim + 1) {
	if (rule_dim == 1)
	    phgError(1, "(%s:%d) unsupported case.\n", __FILE__, __LINE__);
	else if (rule_dim == 2)
	    func_2D(pw[0], pw[1], f);
	else
	    func_3D(pw[0], pw[1], pw[2], f);
	w = pw[rule_dim];
	switch (proj) {
	    case PROJ_NONE:
		for (j = 0; j < dim; j++)
		    res[j] += f[j] * w;
		break;
	    case PROJ_DOT:
		if (rule_dim == 2) {
		    for (j = 0; j < dim; j++)
			res[j] += (f[2 * j + 0] * nv[0] +
				   f[2 * j + 1] * nv[1]) * w;
		    nv += 2;
		}
		else {
		    for (j = 0; j < dim; j++)
			res[j] += (f[3 * j + 0] * nv[0] +
				   f[3 * j + 1] * nv[1] +
				   f[3 * j + 2] * nv[2]) * w;
		    nv += 3;
		}
		break;
	    case PROJ_CROSS:
		if (rule_dim == 2) {
		    for (j = 0; j < dim; j++)
			res[j] += (-f[2 * j + 0] * nv[1] +
				    f[2 * j + 1] * nv[0]) * w;
		    nv += 2;
		}
		else {
		    for (j = 0; j < dim; j += 3) {
			res[j + 0] += (f[j + 1] * nv[2] - f[j + 2] * nv[1]) * w;
			res[j + 1] += (f[j + 2] * nv[0] - f[j + 0] * nv[2]) * w;
			res[j + 2] += (f[j + 0] * nv[1] - f[j + 1] * nv[0]) * w;
		    }
		    nv += 3;
		}
		break;
	}
    }

    phgFree(f);

    return n;
}

int
phgQuadInterfaceRuleInfo(FLOAT *rule, int *dim, FLOAT **pw, FLOAT **nv)
/* Retrieves information about a quadrature rule.
 *
 * 'dim'	- if not NULL, returns the space dimension,
 * 'pw'		- if not NULL, returns a pointer to the list:
 *			2D: {x_0, y_0, w_0, ..., x_n, y_n, w_n}
 *			3D: {x_0, y_0, z_0, w_0, ..., x_n, y_n, z_n, w_n}
 * 'nv'		- if not NULL, returns a pointer to the list of normal vectors.
 *
 * The return value is the number of points of rule.
 */
{
    int n, dim0;
    FLOAT *pw0 = NULL, *nv0 = NULL;
    BOOLEAN nv_flag;

    get_rule_info(rule, &n, &dim0, &nv_flag);
    if (n > 0) {
	assert(dim0 == 2 || dim0 == 3);
	pw0 = rule + RULE_HEADER;
	if (nv_flag)
	    nv0 = pw0 + n * (dim0 + 1);
    }

    if (dim != NULL)
	*dim = dim0;
    if (pw != NULL)
	*pw = pw0;
    if (nv != NULL)
	*nv = nv0;

    return n;
}

int
phgQuadInterfaceRuleMerge(FLOAT **rule, const FLOAT *rule0)
/* merges two rules: rule := rule + rule0 */
{
    int n, n0, dim, dim0;
    BOOLEAN nv_flag, nv_flag0;
    size_t m;
    FLOAT *rule1;

    assert(rule != NULL);

    get_rule_info(*rule, &n, &dim, &nv_flag);
    get_rule_info(rule0, &n0, &dim0, &nv_flag0);

    assert(dim0 == 0 || dim == 0 || dim0 == dim);
    assert(n0 == 0 || n == 0 || nv_flag == nv_flag0);

    if (n0 == 0)
	return n;

    assert(dim0 == 2 || dim0 == 3);
    if (dim == 0)
	dim = dim0;

    nv_flag = nv_flag || nv_flag0;

    m = dim + 1 + (nv_flag ? dim : 0);
    *rule = phgRealloc_(*rule,  (RULE_HEADER + (n + n0) * m) * sizeof(FLOAT),
				(RULE_HEADER + n * m) * sizeof(FLOAT));
    SetRuleHeader(*rule, n + n0, dim, nv_flag);

    rule1 = *rule + RULE_HEADER;
    rule0 += RULE_HEADER;

    if (nv_flag) {
	/* normal vectors */
	if (n > 0)
	    memmove(rule1 + (dim + 1) * (n + n0),
		    rule1 + (dim + 1) * n,
		    dim * n * sizeof(FLOAT));
	if (n0 > 0)
	    memcpy (rule1 + (2 * dim + 1) * n + (dim + 1) * n0,
		    rule0 + (dim + 1) * n0,
		    dim * n0 * sizeof(FLOAT));
    }

    /* rule data */
    if (n0 > 0)
	memcpy(rule1 + (dim + 1) * n, rule0, n0 * (dim + 1) * sizeof(FLOAT));

    return n + n0;
}

/******************************************************************************/

static int
check_element_2d(FUNC_2D ls_func, FUNC_2D ls_grad, int ls_order,
		 int elem_type, void *elem_data)
{
    int dim = 2, nv, i, j, clockwise;
    FLOAT (*v)[dim], a, b, d, f;
    FLOAT tmp[4][dim];

    v = elem_data;
    nv = 3;
    if (elem_type == ET_RECTANGLE) {
	tmp[0][0] = v[0][0];
	tmp[0][1] = v[0][1];
	tmp[1][0] = v[1][0];
	tmp[1][1] = v[0][1];
	tmp[2][0] = v[1][0];
	tmp[2][1] = v[1][1];
	tmp[3][0] = v[0][0];
	tmp[3][1] = v[1][1];
	v = tmp;
	nv = 4;
    }

    /* Step 1: check signs of the vertices */
    a = f = 0.0;
    j = nv - 1;
    for (i = 0, j = nv - 1; i < nv; j = i, i++) {
	/* check sign of vertex i */
	ls_func(v[i][0], v[i][1], &d);
	if (Fabs(d) <= FLOAT_EPSILON * 1000.0 || f * d < 0.)
	    return 0;
	f = d;
	/* determine direction of the edges by computing the (2x) signed area */
	a += v[j][0] * v[i][1] - v[j][1] * v[i][0];
    }
    clockwise = (a < 0.);

    /* Step 2: check signs of the edge centers moved outwards */
    /* loop on edges, checking the sign of the middle point moved outward */
    for (i = 0, j = nv - 1; i < nv; j = i, i++) {
	a = (v[i][0] - v[j][0]) * 0.125;
	b = (v[i][1] - v[j][1]) * 0.125;
	if (clockwise)
	    ls_func((v[i][0] + v[j][0]) * 0.5 - b,
		    (v[i][1] + v[j][1]) * 0.5 + a, &d);
	else
	    ls_func((v[i][0] + v[j][0]) * 0.5 + b,
		    (v[i][1] + v[j][1]) * 0.5 - a, &d);
	if (Fabs(d) <= FLOAT_EPSILON * 1000.0 || f * d < 0.)
	    return 0;
    }

    return d < 0. ? -1 : 1;
}

static BOOLEAN
check_element_face(FUNC_3D ls_func, FLOAT f, const FLOAT v0[],
			const FLOAT v1[], const FLOAT v2[], const FLOAT v3[])
/* check the center of the face (v0,v1,v2), v3 is a vertex not on the face
 * which is used to determine the outward direction */
{
    FLOAT d1[3] = {v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]},
	  d2[3] = {v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]},
	  d [3] = {v3[0] - v0[0], v3[1] - v0[1], v3[2] - v0[2]};
    FLOAT n[3] = {d1[1] * d2[2] - d1[2] * d2[1],
		  d1[0] * d2[2] - d1[2] * d2[0],
		  d1[0] * d2[1] - d1[1] * d2[0]};
    FLOAT a, b;

    a = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2];
    b = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2];
    if (a < b)
	a = b;
    a = 0.125 / Sqrt(a);
    if (d[0] * n[0] + d[1] * n[1] + d[2] * n[2] > 0.)
	a = -a; /* reverse normal direction to point outwards */

    b = _F(1.) / _F(3.);
    ls_func(b * (v0[0] + v1[0] + v2[0]) + a * n[0],
	    b * (v0[1] + v1[1] + v2[1]) + a * n[1],
	    b * (v0[2] + v1[2] + v2[2]) + a * n[2], &a);

    return Fabs(a) <= FLOAT_EPSILON * 1000.0 || f * a < 0. ? TRUE : FALSE;
}

#define check_element_edge(ls, ls_grad, ls_order, v0, v1) \
	phgQuadInterfaceGetEdgeCuts(ls, ls_order, ls_grad, v0, v1, -1.0, NULL)

static int
check_element_3d(FUNC_3D ls_func, FUNC_3D ls_grad, int ls_order,
		 int elem_type, void *elem_data)
{
    int dim = 3, nv, i;
    FLOAT (*v)[dim], d, f;
    FLOAT tmp[8][dim];

    v = elem_data;
    nv = 4;
    if (elem_type == ET_CUBOID) {
	tmp[0][0] = v[0][0];
	tmp[0][1] = v[0][1];
	tmp[0][2] = v[0][2];
	tmp[1][0] = v[1][0];
	tmp[1][1] = v[0][1];
	tmp[1][2] = v[0][2];
	tmp[2][0] = v[1][0];
	tmp[2][1] = v[1][1];
	tmp[2][2] = v[0][2];
	tmp[3][0] = v[0][0];
	tmp[3][1] = v[1][1];
	tmp[3][2] = v[0][2];

	tmp[4][0] = v[0][0];
	tmp[4][1] = v[0][1];
	tmp[4][2] = v[1][2];
	tmp[5][0] = v[1][0];
	tmp[5][1] = v[0][1];
	tmp[5][2] = v[1][2];
	tmp[6][0] = v[1][0];
	tmp[6][1] = v[1][1];
	tmp[6][2] = v[1][2];
	tmp[7][0] = v[0][0];
	tmp[7][1] = v[1][1];
	tmp[7][2] = v[1][2];

	v = tmp;
	nv = 8;
    }

    /* Step 1: check signs of the vertices */
    f = 0.0;
    for (i = 0; i < nv; i++) {
	/* check sign of vertex i */
	ls_func(v[i][0], v[i][1], v[i][2], &d);
	if (Fabs(d) <= FLOAT_EPSILON * 1000.0 || f * d < 0.)
	    return 0;
	f = d;
    }

    /* Step 2: check signs of the face centers moved outwards */
    if (elem_type == ET_TETRA) {
	if (check_element_face(ls_func, f, v[0], v[1], v[2], v[3]) ||
	    check_element_face(ls_func, f, v[0], v[1], v[3], v[2]) ||
	    check_element_face(ls_func, f, v[0], v[2], v[3], v[1]) ||
	    check_element_face(ls_func, f, v[1], v[2], v[3], v[0]))
	    return 0;
    }
    else {
	if (check_element_face(ls_func, f, v[0], v[1], v[3], v[4]) ||
	    check_element_face(ls_func, f, v[0], v[1], v[4], v[3]) ||
	    check_element_face(ls_func, f, v[0], v[3], v[4], v[1]) ||
	    check_element_face(ls_func, f, v[6], v[2], v[5], v[7]) ||
	    check_element_face(ls_func, f, v[6], v[2], v[7], v[5]) ||
	    check_element_face(ls_func, f, v[6], v[5], v[7], v[2]))
	    return 0;
    }

    /* Step 3: check whether the edges intersect the interface */
    if (elem_type == ET_TETRA) {
	if (check_element_edge(ls_func, ls_grad, ls_order, v[0], v[1]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[0], v[2]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[0], v[3]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[1], v[2]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[1], v[3]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[2], v[3]))
	    return 0;
    }
    else {
	if (check_element_edge(ls_func, ls_grad, ls_order, v[0], v[1]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[1], v[2]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[2], v[3]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[3], v[0]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[4], v[5]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[5], v[6]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[6], v[7]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[7], v[4]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[0], v[4]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[1], v[5]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[2], v[6]) ||
	    check_element_edge(ls_func, ls_grad, ls_order, v[3], v[7]))
	    return 0;
    }

    return d < 0. ? -1 : 1;
}

int
phgQuadInterfaceMarkElement(void *ls_func, int ls_order, void *ls_grad,
			    ELEM_TYPE elem_type, void *elem_data)
/* returns relative position of the element w.r.t. the interface:
 *	-1: in \Omega^-, 1: in \Omega^+, 0: (possibly) on the interface.
 *
 * ls_func:	pointer to the levelset function
 * ls_grad:	pointer to the gradient of the levelset function
 * ls_order:	polynomial order of the levelset function
 * elem_type:	ET_TRIANGLE, ET_RECTANGLE, ET_TETRA, ET_CUBOID
 * elem_data:	coordinates of the vertices of the element, in the same format
 *		as the corresponding argument in phgQuadInterface*
 *
 * This function requires the interface to be sufficiently flat in the element.
 */
{
    if (elem_type == ET_TRIANGLE || elem_type == ET_RECTANGLE)
	return check_element_2d(ls_func, ls_grad, ls_order,
						elem_type, elem_data);
    return check_element_3d(ls_func, ls_grad, ls_order, elem_type, elem_data);
}

int
phgQuadInterfaceGetEdgeCuts(FUNC_3D ls_func, int ls_order, FUNC_3D ls_grad,
		const FLOAT v0[], const FLOAT v1[], FLOAT tol, FLOAT **cuts)
/* returns the number of intersection points (cut points) of the edge v0-v1
 * with the interface, ignoring the multiplicity.
 *
 * "cuts", if not NULL, will point to a dynamically allocated FLOAT array
 * containing the cuts with values in [0,1].
 *
 * If "tol" < 0.0, then the function just checks whether the edge and the
 * interface intersect, and returns !0 if yes and 0 if no, and in this case
 * the cuts are not actually computed and "cuts" must be NULL.
 */
{
    int i, ret;
    FLOAT t, *c = NULL;

    assert(tol >= 0.0 || cuts == NULL);

#if 0
ls_order += 8;	/* ./quad_test3-p4est -debug -dim 3 -n {3,5} fails */
#endif
    /* construct ls_func(v0+t*(v1-v0)) (an univariate polynomial in t) */
    c = phgAlloc((ls_order + 1) * sizeof(*c));
    for (i = 0; i <= ls_order; i++) {
	t = i / (FLOAT)ls_order;
	ls_func(v0[0] + t * (v1[0] - v0[0]),
		v0[1] + t * (v1[1] - v0[1]),
		v0[2] + t * (v1[2] - v0[2]), &c[i]);
    }
    ret = phgUnivariatePolyRoots(ls_order, _F(0.), _F(1.), tol, c, NULL);
    if (cuts != NULL)
	*cuts = c;
    else
	phgFree(c);

    return ret;
}

/*---------------------------------------------------------------------------*/

static struct {
    FUNC_3D	ls, ls_grad;
    DOF		*dof_ls, *dof_ls_grad;
    ELEMENT	*e;
    FLOAT const	(*tet)[Dim];
    /* The 3 orthonormal vectors for the local coordinates (w = face normal) */
    FLOAT	u[Dim], v[Dim], w[Dim], w0;
    int		ls_order, face, v0, v1, v2;
} *qf_ctx = NULL;
#if USE_OMP
# pragma omp threadprivate(qf_ctx)
#endif	/* USE_OMP */

static void
face2xyz(FLOAT x1, FLOAT x2, FLOAT *x, FLOAT *y, FLOAT *z)
{
    *x = qf_ctx->u[0] * x1 + qf_ctx->v[0] * x2 + qf_ctx->w[0] * qf_ctx->w0;
    *y = qf_ctx->u[1] * x1 + qf_ctx->v[1] * x2 + qf_ctx->w[1] * qf_ctx->w0;
    *z = qf_ctx->u[2] * x1 + qf_ctx->v[2] * x2 + qf_ctx->w[2] * qf_ctx->w0;
}

static void
ls2(FLOAT x1, FLOAT x2, FLOAT *val)
{
    FLOAT lambda[Dim + 1], x, y, z;

    face2xyz(x1, x2, &x, &y, &z);
    if (qf_ctx->ls != NULL) {
	qf_ctx->ls(x, y, z, val);
    }
    else {
	assert(qf_ctx->dof_ls != NULL && qf_ctx->e != NULL);
	phgGeomXYZ2Lambda(qf_ctx->dof_ls->g, qf_ctx->e, x, y, z, lambda);
	phgDofEval(qf_ctx->dof_ls, qf_ctx->e, lambda, val);
    }
}

static void
ls2_grad(FLOAT x1, FLOAT x2, FLOAT *grad)
{
    FLOAT lambda[Dim + 1], G[Dim], x, y, z;

    face2xyz(x1, x2, &x, &y, &z);
    if (qf_ctx->ls_grad != NULL) {
	qf_ctx->ls_grad(x, y, z, G);
    }
    else {
	assert(qf_ctx->dof_ls_grad != NULL && qf_ctx->e != NULL);
	phgGeomXYZ2Lambda(qf_ctx->dof_ls_grad->g, qf_ctx->e, x, y, z, lambda);
	phgDofEval(qf_ctx->dof_ls_grad, qf_ctx->e, lambda, G);
    }

    /* Project back to face */
    grad[0] = qf_ctx->u[0] * G[0] + qf_ctx->u[1] * G[1] + qf_ctx->u[2] * G[2];
    grad[1] = qf_ctx->v[0] * G[0] + qf_ctx->v[1] * G[1] + qf_ctx->v[2] * G[2];
}

static void
quad_face(FLOAT const tet[Dim + 1][Dim], int face, int quad_order,
	  FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
{
    FLOAT tri[3][2];
    FLOAT d, *pw, *pw2, *pw3, *nv2, *nv3, *nv;
    FLOAT *u = qf_ctx->u, *v = qf_ctx->v, *w = qf_ctx->w;
    int k, i, np;

    qf_ctx->tet = tet;
    qf_ctx->v0 = GetFaceVertex(face, 0);
    qf_ctx->v1 = GetFaceVertex(face, 1);
    qf_ctx->v2 = GetFaceVertex(face, 2);
    qf_ctx->face = face;

    /* local coordinates: (u,v) = the plane, w = face normal
     *		u := v1-v0,
     *		w := face normal=u \cross (v2-v0),
     *		v := w \cross u
     *
     * Note: in previous revisions of this file (<= r1.611) non-orthogonal
     * barycentric coordinates of the triangle were used, which led to
     * inaccurate results because under non-orthogonal coordinates the normal
     * vectors of the projected interface are not preserved, which breaks the
     * algorithm for chosing integration directions to avoid essential
     * singularities.
     */
    for (k = 0; k < Dim; k++) {
	/* u := v1 - v0, v := v2 - v0 */
	u[k] = tet[qf_ctx->v1][k] - tet[qf_ctx->v0][k];
	v[k] = tet[qf_ctx->v2][k] - tet[qf_ctx->v0][k];
    }
    w[0] = u[1] * v[2] - u[2] * v[1];
    w[1] = u[2] * v[0] - u[0] * v[2];
    w[2] = u[0] * v[1] - u[1] * v[0];
    /* normalize u */
    d = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    assert(d != 0.0);
    d = 1.0 / Sqrt(d);
    u[0] *= d;
    u[1] *= d;
    u[2] *= d;
    /* normalize w */
    d = w[0] * w[0] + w[1] * w[1] + w[2] * w[2];
    assert(d != 0.0);
    d = 1.0 / Sqrt(d);
    w[0] *= d;
    w[1] *= d;
    w[2] *= d;
    /* v := w \cross u */
    v[0] = w[1] * u[2] - w[2] * u[1];
    v[1] = w[2] * u[0] - w[0] * u[2];
    v[2] = w[0] * u[1] - w[1] * u[0];

    qf_ctx->w0 = tet[qf_ctx->v0][0] * w[0] +
		 tet[qf_ctx->v0][1] * w[1] +
		 tet[qf_ctx->v0][2] * w[2];

    /* project triangle to new coordinates */
    for (k = 0; k < 3; k++) {
	i = GetFaceVertex(face, k);
	tri[k][0] = tet[i][0] * u[0] + tet[i][1] * u[1] + tet[i][2] * u[2];
	tri[k][1] = tet[i][0] * v[0] + tet[i][1] * v[1] + tet[i][2] * v[2];
    }
#if 0
#warning
phgInfo(-1, "uvw=[[%g,%g,%g];[%g,%g,%g];[%g,%g,%g]], w0=%g\n", u[0], u[1], u[2], v[0], v[1], v[2], w[0], w[1], w[2], qf_ctx->w0);
phgInfo(-1, "projected triangle: [%g,%g], [%g,%g], [%g,%g]\n", tri[0][0], tri[0][1], tri[1][0], tri[1][1], tri[2][0], tri[2][1]);
#endif

    if (qf_ctx->ls != NULL || qf_ctx->dof_ls != NULL) {
	phgQuadInterfaceTriangle(ls2, qf_ctx->ls_order, ls2_grad,
				 (void *)tri, quad_order,
				 rule_m, rule_0, rule_p);
    }
    else {
	phgQuadInterfaceTriangle(NULL, 0, NULL, (void *)tri, quad_order,
				 rule_m, rule_0, rule_p);
    }

    /* convert 2D barycentric coordinates to 3D xyz */
    for (k = -1; k <= 1; k++) {
	FLOAT **rule = k < 0 ? rule_m : (k > 0 ? rule_p : rule_0);
	if (rule == NULL)
	    continue;
	np = phgQuadInterfaceRuleInfo(*rule, NULL, &pw2, &nv2);
	pw3 = phgAlloc((Dim + 1 + (nv2 == NULL ? 0 : Dim)) * np * sizeof(*pw3));
	nv = nv3 = nv2 == NULL ? NULL : pw3 + (Dim + 1) * np;
	for (i = 0, pw = pw3; i < np; i++, pw += Dim + 1, pw2 += 2 + 1) {
	    face2xyz(pw2[0], pw2[1], pw, pw + 1, pw + 2);
	    pw[3] = pw2[2];
	    if (nv2 == NULL)
		continue;
	    /* convert vector on face nv2[2] to 3D coordinates */
	    nv[0] = qf_ctx->u[0] * nv2[0] + qf_ctx->v[0] * nv2[1];
	    nv[1] = qf_ctx->u[1] * nv2[0] + qf_ctx->v[1] * nv2[1];
	    nv[2] = qf_ctx->u[2] * nv2[0] + qf_ctx->v[2] * nv2[1];
	    nv += 3;
	    nv2 += 2;
	}
	phgFree(*rule);
	*rule = phgQuadInterfaceRuleCreate(Dim, np, pw3, Dim + 1,
			pw3 + Dim, Dim + 1, 1.0, nv3, Dim);
	phgFree(pw3);
    }

    phgFree(qf_ctx);
    qf_ctx = NULL;
}

void
phgQuadInterfaceTetrahedronFace(FUNC_3D ls, int ls_order, FUNC_3D ls_grad,
		FLOAT const tet[Dim + 1][Dim], int face,
		int quad_order, FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
/* returns qudrature rules *rule_m, rule_0 and *rule_p, for numerical
 * integration in Face\cap\Omega^-, Face\cap\Gamma and Face\cap\Omega^+,
 * respectively.
 *
 * Note: in the resulting rules, the weights are scaled to match the area of
 * the face. 
 *
 * If ls == NULL, then returns a quadrature rule for the whole face in
 * *rule_m, or *rule_p if rule_m == NULL.
 */
{
    qf_ctx = phgCalloc(1, sizeof(*qf_ctx));
    qf_ctx->ls = ls;
    qf_ctx->ls_grad = ls_grad;
    qf_ctx->ls_order = ls_order;
    quad_face(tet, face, quad_order, rule_m, rule_0, rule_p);
}

void
phgQuadInterfaceFace(DOF *ls, DOF *ls_grad, ELEMENT *e, int face,
		int quad_order, FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
{
    FLOAT tet[Dim + 1][Dim];
    GRID *g = ls_grad->g;
    int i;

    qf_ctx = phgCalloc(1, sizeof(*qf_ctx));
    qf_ctx->dof_ls = ls;
    qf_ctx->dof_ls_grad = ls_grad;
    qf_ctx->e = e;
    assert(ls != NULL || ls_grad != NULL);
    qf_ctx->ls_order = LSOrder(ls, ls_grad, e);
    for (i = 0; i < Dim + 1; i++)
	memcpy(tet[i], g->verts + e->verts[i], Dim * sizeof(FLOAT));
    quad_face((void *)tet, face, quad_order, rule_m, rule_0, rule_p);
}
