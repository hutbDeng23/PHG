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

/* $Id: phg-to-p4est.c,v 1.261 2022/09/18 11:51:35 zlb Exp $ */


#define NEED_GetVariableArgs

#if !defined(DONT_CAST_phgMapE2G)
# define DONT_CAST_phgMapE2G
#endif	/* !defined(DONT_CAST_phgMapE2G) */

#include <phg.h>

#if HAVE_P4EST	/* whole file */

#include "phg/io.h"

#include <p8est_vtk.h>
#include <p8est_iterate.h>
#include <p8est_mesh.h>
#include <p8est_extended.h>
#include <p8est_bits.h>
#include <p8est_communication.h>

#include "phg-to-p4est.h"

#include <sys/mman.h>		/* mprotect() */

static void
xyz_to_ref(GRID *g, INT eno, const FLOAT x[], FLOAT xref[])
/* returns in xref[] the coordinates of point x[] in the ref. element (-1,1)^3.
 * dx[] is an optional array which, if not NULL, returns the relative size
 * of the element w.r.t. the reference element. */
{
    FLOAT (*corners)[3] = g->elems[eno]->corners;

    xref[0] = (x[0] - corners[0][0]) / (-corners[0][0] + corners[1][0]);
    xref[1] = (x[1] - corners[0][1]) / (-corners[0][1] + corners[1][1]);
    xref[2] = (x[2] - corners[0][2]) / (-corners[0][2] + corners[1][2]);
}

static FLOAT *
ref_to_xyz(GRID *g, INT eno, const FLOAT xref[], FLOAT x[])
{
    FLOAT (*cub)[Dim] = phgGeomGetCorners(g, g->elems[eno], NULL);

    x[0] = cub[0][0] + xref[0] * (cub[1][0] - cub[0][0]);
    x[1] = cub[0][1] + xref[1] * (cub[1][1] - cub[0][1]);
    x[2] = cub[0][2] + xref[2] * (cub[1][2] - cub[0][2]);

    return x;
}

static FLOAT
basis1D(int n, int k, FLOAT x, const FLOAT *points, FLOAT *coeff)
/* returns value of the k-th 1D basis function at the point x */
{
    int i;
    FLOAT d, pk;

#if 0	/* not necessary, it's done in DOF_TYPE_init() */
#if USE_OMP
# pragma omp critical
#endif	/* USE_OMP */
#endif
    if (*coeff == 0.) {
	d = 1.0;
	pk = points[k];
	for (i = 0; i < k; i++)
	    d *= (pk - points[i]);
	for (i = k + 1; i <= n; i++)
	    d *= (pk - points[i]);
	*coeff = 1.0 / d;
    }

    d = *coeff;
    for (i = 0; i < k; i++)
	d *= (x - points[i]);
    for (i = k + 1; i <= n; i++)
	d *= (x - points[i]);

    return d;
}

static FLOAT
grad1D(int n, int k, FLOAT x, const FLOAT *points, FLOAT coeff, FLOAT bas)
/* returns value of derivative of the k-th 1D basis function at the point x */
{
    int i, j;
    FLOAT d;

    for (j = 0; j <= n; j++) {
	if (j == k)
	    continue;
	if (x == points[j])
	    break;
    }
    if (j <= n) {
	/* x == points[j] */
	d = coeff;
	for (i = 0; i <= n; i++) {
	    if (i == j || i == k)
		continue;
	    d *= (x - points[i]);
	}
    }
    else {
	/* x != points[i], forall i */
	d = 0.0;
	for (i = 0; i <= n; i++) {
	    if (i == k)
		continue;
	    d += bas / (x - points[i]);
	}
    }

    return d;
}

/*-----------------------------------------------------------------------------
DGPk type: Pk_KIND == 0: monomials, Pk_KIND == 1: Legendre.
To use monomials for DGPk:
  (cd p4est; rm -f phg-to-p4est.o; make USER_CFLAGS=-DPk_KIND=0 phg-to-p4est.o)
------------------------------------------------------------------------------*/
#ifndef Pk_KIND
# define Pk_KIND 1	/* Legendre polynomials by default */
#endif

static FLOAT inline
power(FLOAT x, int p, int p_max, int dim, BOOLEAN dflag)
/* returns either P_p(x) or P_p'(x), caching computed values for given x,
 * where P_p(x) := x^p if Pk_KIND==0, or L_p(2*x-1) if Pk_KIND==1.
 *
 * "p_max" is usually set to type->order in the calling function.
 * "dim" in [0, Dim) denotes index of space dimension (for selecting cache)
 * "dflag" == FALSE: P_p(x), "dflag" == TRUE: P_p'(x) */
{
    /* Note: cache one x value for eacg space dimension */
    static int order = -1;			/* cached order */
    static struct {
	FLOAT	x[Dim];				/* cached x */
	FLOAT	P[Dim][DOF_TYPE_ORDER_MAX + 1];	/* cache for polynomials */
	FLOAT	D[Dim][DOF_TYPE_ORDER_MAX + 1];	/* cache for derivatives */
    } cache;
#if USE_OMP
# pragma omp threadprivate(cache, order)
#endif
    int i;
    FLOAT *P, *D;

    assert(dim >= 0 && dim < Dim);
    assert(p >= 0 && p <= p_max);
    assert(p_max <= DOF_TYPE_ORDER_MAX);
    if (x < -FLOAT_EPSILON * 10. || x > 1. + FLOAT_EPSILON * 10.)
	phgError(1, "%s: unexpected (x=%g)\n", __func__, (double)x);
    if (x < 0.)
	x = 0.;
    else if (x > 1.)
	x = 1.;

    if (p_max > order) {
	/* clear cache */
	order = p_max;
	for (i = 0; i < Dim; i++)
	    cache.x[i] = -2.0;
    }

    P = cache.P[dim];
    D = cache.D[dim];

    if (x != cache.x[dim]) {
#if Pk_KIND == 0
	P[0] = 1.0;
	D[0] = 0.0;
	for (i = 1; i <= order; i++) {
	    D[i] = P[i - 1] * i;
	    P[i] = P[i - 1] * x;
	}
#else	/* Pk_KIND == */
	phgLegendreP(order, x + x - 1., P, D);
	/* scale derivatives with respect to interval [0,1] */
	for (i = 1; i <= order; i++)
	    D[i] *= 2.0;
#endif	/* Pk_KIND == */
	cache.x[dim] = x;
    }

    return (dflag ? D : P)[p];
}

static void
basisFunc(DOF *u, FLOAT x, FLOAT y, FLOAT z, int index, FLOAT *funcValue)
/* 标准单元上的基函数 */
{
    DOF_TYPE *type = u->type;

    if (type->base_type != NULL) {
	DOF dummy;
	int q = type->dim / type->base_type->dim;
	dummy.type = type->base_type;
	bzero(funcValue, type->dim * sizeof(*funcValue));
	basisFunc(&dummy, x, y, z, index / q,
			funcValue + (index % q) * type->base_type->dim);
	return;
    }

    /* the code below is for DGQn and DGPn */
    if (type->is_nodal) {
	assert(type->tensor_product && type->points != NULL);
#define TENSOR_PRODUCT \
	int i, j, k, ox, oy, oz; \
	FLOAT *px, *py, *pz, *cx, *cy, *cz; \
	if (type->orders == NULL) \
	    phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__); \
	ox = type->orders[0]; \
	oy = type->orders[1]; \
	oz = type->orders[2]; \
	/* Note: index := i + (j + k * (oy+1)) * (ox+1) */ \
	i = index % (ox + 1); \
	index /= (ox + 1); \
	j = index % (oy + 1); \
	k = index / (oy + 1); \
	px = type->points; \
	py = px + ox + 1; \
	pz = py + oy + 1; \
	cx = type->work; \
	cy = cx + ox + 1; \
	cz = cy + oy + 1;
	TENSOR_PRODUCT
	*funcValue = basis1D(ox, i, x, px, cx + i) *
		     basis1D(oy, j, y, py, cy + j) *
		     basis1D(oz, k, z, pz, cz + k);
    }
    else {
	/* Note: x^i y^j z^k, with (i,j,k) encoded in orders[index] as:
		i + (j + k * (n + 1)) * (n + 1) */
	int i, j, k, a, n = type->order;
	assert(!type->tensor_product);
	a = type->orders[index];
	i = a % (n + 1);
	a /= n + 1;
	j = a % (n + 1);
	k = a / (n + 1);
	*funcValue = power(x, i, n, 0, FALSE) *
		     power(y, j, n, 1, FALSE) *
		     power(z, k, n, 2, FALSE);
    }
}

static void
gradBasisFunc(DOF *u, FLOAT x, FLOAT y, FLOAT z, int index, FLOAT *grad)
/* 标准单元基函数梯度 */
{
    DOF_TYPE *type = u->type;

    if (type->dim == Dim) {
	DOF dummy;
	int q = type->dim / type->base_type->dim;
	dummy.type = type->base_type;
	bzero(grad, Dim * type->dim * sizeof(*grad));
	gradBasisFunc(&dummy, x, y, z, index / q,
			grad + (index % q) * type->base_type->dim * Dim);
	return;
    }

    /* the code below is for DGQn and DGPn */
    if (type->is_nodal) {
	FLOAT bx, by, bz;
	TENSOR_PRODUCT
#undef TENSOR_PRODUCT
	bx = basis1D(ox, i, x, px, cx + i);
	by = basis1D(oy, j, y, py, cy + j);
	bz = basis1D(oz, k, z, pz, cz + k);
	grad[0] = grad1D(ox, i, x, px, cx[i], bx) * by * bz;
	grad[1] = grad1D(oy, j, y, py, cy[j], by) * bx * bz;
	grad[2] = grad1D(oz, k, z, pz, cz[k], bz) * bx * by;
    }
    else {
	/* Note: the powers (i,j,k) is encoded in orders[index] as:
		i + (j + k * (n + 1)) * (n + 1) */
	int i, j, k, a, n = type->order;
	a = type->orders[index];
	i = a % (n + 1);
	a /= n + 1;
	j = a % (n + 1);
	k = a / (n + 1);
	grad[0] = power(x, i, n, 0, TRUE) *
		  power(y, j, n, 1, FALSE) *
		  power(z, k, n, 2, FALSE);
	grad[1] = power(x, i, n, 0, FALSE) *
		  power(y, j, n, 1, TRUE) *
		  power(z, k, n, 2, FALSE);
	grad[2] = power(x, i, n, 0, FALSE) *
		  power(y, j, n, 1, FALSE) *
		  power(z, k, n, 2, TRUE);
    }
}

/*------------------------ The QC interface for P4EST ------------------------*/

static void
to_ref(QCACHE *qc, INT eno, FLOAT x[], FLOAT xref[])
{
    DOF *u = qc->fe;
    GRID *g = u->g;

    xyz_to_ref(g, eno, x, xref);
}

static void
from_ref(QCACHE *qc, INT eno, FLOAT xref[], FLOAT x[])
{
    DOF *u = qc->fe;
    GRID *g = u->g;

    ref_to_xyz(g, eno, xref, x);
}


static int
b_dim(QCACHE *qc) { return ((DOF *)qc->fe)->type->dim; }

static int
f_dim(QCACHE *qc, void *fe) { return ((DOF *)fe)->dim * b_dim(qc); }

static int
n_bas(QCACHE *qc, INT eno)
{
    return ((DOF *)qc->fe)->type->nbas;
}

static void
bas(QCACHE *qc, INT eno, const FLOAT xref[], FLOAT *buffer, int stride)
{
    int j, nb = n_bas(qc, eno);

    if (stride <= 0)
	stride = ((DOF *)qc->fe)->type->dim;

    for (j = 0; j < nb; j++, buffer += stride)
	basisFunc(qc->fe, xref[0], xref[1], xref[2], j, buffer);
}

static void
grd(QCACHE *qc, INT eno, const FLOAT xref[], FLOAT *buffer, int stride)
{
    int j, k, dim = b_dim(qc), nb = n_bas(qc, eno);
    FLOAT dx, dy, dz;
    FLOAT (*corners)[3] = ((DOF *)(qc->fe))->g->elems[eno]->corners;

    if (stride <= 0)
	stride = ((DOF *)qc->fe)->type->dim * Dim;

    dx = 1.0 / (-corners[0][0] + corners[1][0]);
    dy = 1.0 / (-corners[0][1] + corners[1][1]);
    dz = 1.0 / (-corners[0][2] + corners[1][2]);
    for (j = 0; j < nb; j++, buffer += stride) {
	gradBasisFunc(qc->fe, xref[0], xref[1], xref[2], j, buffer);
	for (k = 0; k < dim; k++) {
	    buffer[k * Dim + 0] *= dx;
	    buffer[k * Dim + 1] *= dy;
	    buffer[k * Dim + 2] *= dz;
	}
    }
    return;
}

static void
fun(QCACHE *qc, void *fe, INT eno, const FLOAT x[], FLOAT *buffer)
{
    GRID *g = ((DOF *)fe)->g;
    ELEMENT *e = g->elems[eno];

    phgDofEval(fe, e, x, buffer);
}

QDESC QD_P4EST_ = {"P4EST", to_ref, from_ref, b_dim, f_dim, n_bas, bas, grd,
			fun, Dim, Dim};

/*-------------------- End of the QC interface for P4EST --------------------*/

/*--------------------------- PHG_to_p8est wrapper ---------------------------*/

INT p8est_i;	/* loop variable used by ForAllElements */

/* DOF types */

/* dummy DOF_TYPEs */
static DOF_TYPE DOF_ANALYTIC_ = {
	"ANALYTIC",	/* name */
	NULL,		/* grad_type */
	NULL,		/* base_type */
	NULL,		/* points */
	NULL,		/* work */
	NULL,		/* orders */
	NULL,		/* pvt */
	1,		/* dim */
	0,		/* nbas */
	-1,		/* order */
	FALSE,		/* tensor_product */
	FALSE,		/* is_nodal */
	FE_L2		/* FE_SPACE */ };
DOF_TYPE *DOF_ANALYTIC = &DOF_ANALYTIC_;

static DOF_TYPE DOF_CONSTANT_ = {
	"CONSTANT",	/* name */
	NULL,		/* grad_type */
	NULL,		/* base_type */
	NULL,		/* points */
	NULL,		/* work */
	NULL,		/* orders */
	NULL,		/* pvt */
	1,		/* dim */
	0,		/* nbas */
	-1,		/* order */
	FALSE,		/* tensor_product */
	FALSE,		/* is_nodal */
	FE_L2		/* FE_SPACE */ };
DOF_TYPE *DOF_CONSTANT = &DOF_CONSTANT_;

DOF_TYPE **DOF_DGQn = NULL, **DOF_DGPn = NULL;
DOF_TYPE *DOF_DEFAULT = NULL;

int p8est_edge_no[NVert][NVert];

static void
DOF_TYPE_init(void)
{
    int i, j, k, l, o, n = DOF_TYPE_ORDER_MAX + 1;
    DOF_TYPE *type;
    char name[128];
    FLOAT *p;

    DOF_DGQn = phgAlloc(n * sizeof(*DOF_DGQn));
    DOF_DGPn = phgAlloc(n * sizeof(*DOF_DGPn));
    for (o = 0; o < n; o++) {
	/* DGQn */
	DOF_DGQn[o] = type = phgCalloc(1, sizeof(*type));
	type->tensor_product = TRUE;
	sprintf(name, "DGQ%d", o);
	type->name = strdup(name);	/* name */
	type->grad_type = type;		/* grad_type */
	type->base_type = NULL;		/* base_type */
	/* points */
	type->is_nodal = TRUE;
	type->points = p = phgCalloc(3 * (o + 1), sizeof(*type->points));
	if (o == 0) {
	    p[0] = p[1] = p[2] = 0.5;
	}
	else {
/*******************************************************************************
 * By default Gauss-Lobatto nodes are used for DGQk elements.
 *
 * To use uniform nodes:
 * (cd p4est; rm -f phg-to-p4est.o; make USER_CFLAGS=-DQk_KIND=0 phg-to-p4est.o)
 *
 * To use Gauss-Legendre nodes:
 * (cd p4est; rm -f phg-to-p4est.o; make USER_CFLAGS=-DQk_KIND=1 phg-to-p4est.o)
 ******************************************************************************/
#ifndef Qk_KIND
# define Qk_KIND 2	/* 0: uniform, 1: Gauss-Legendre, 2: Gauss-Lobatto */
#endif	/* !defined(Qk_KIND) */
#if Qk_KIND == 0	/* uniform points */
	    FLOAT d = 1.0 / o;
	    for (j = 0; j <= o; j++)
		p[j] = p[o + 1 + j] = p[2 * (o + 1) + j] = j * d;
#elif Qk_KIND == 1	/* Gauss-Legendre points */
	    QUAD *q = phgQuadGetQuad1D(2 * o + 1);
	    for (j = 0; j <= o; j++)
		p[j] = p[o + 1 + j] = p[2 * (o + 1) + j] = q->points[2 * j];
#else	/* Gauss-Lobatto points */
# include "gauss-lobatto.inc"
	    memcpy(p, Gtab[o], (o + 1) * sizeof(*p));
	    qsort(p, o + 1, sizeof(*p), phgCompFLOAT);
	    memcpy(p + o + 1, p, (o + 1) * sizeof(*p));
	    memcpy(p + 2 * (o + 1), p, (o + 1) * sizeof(*p));
#endif	/* Qk_KIND == */
	}
	type->work = phgCalloc(3 * (o + 1), sizeof(*type->work));
	/* orders */
	type->orders = phgCalloc(3, sizeof(*type->orders));
	type->orders[0] = type->orders[1] = type->orders[2] = o;
	type->pvt = NULL;		/* unused */
	type->dim = 1;			/* dim */
	type->nbas = (o + 1) * (o + 1) * (o + 1);	/* nbas */
	type->order = o /*+ o + o*/;		/* order */
	type->fe_space = FE_L2;		/* FE_SPACE */
	/* force initialization of type->work */
	for (j = 0; j <= o; j++) {
	    basis1D(o, j, 0.0, p,		type->work + j);
	    basis1D(o, j, 0.0, p + o + 1,	type->work + o + 1 + j);
	    basis1D(o, j, 0.0, p + 2 * (o + 1), type->work + 2 * (o + 1) + j);
	}

	/* DGPn */
	DOF_DGPn[o] = type = phgCalloc(1, sizeof(*type));
	type->tensor_product = FALSE;
	sprintf(name, "DGP%d", o);
	type->name = strdup(name);	/* name */
	type->grad_type = DOF_DGPn[o == 0 ? o : o - 1]; /* grad_type */
	type->base_type = NULL;		/* base_type */
	/* points */
	type->is_nodal = FALSE;
	type->points = NULL;
	type->work = NULL;		/* set in phgDofSetDataByFunction_ */
	type->pvt = NULL;		/* set in phgDofSetDataByFunction_ */
	/* orders */
	type->dim = 1;			/* dim */
	type->nbas = (o + 1) * (o + 2) * (o + 3) / 6;	/* nbas */
	type->order = o;		/* order */
	type->fe_space = FE_L2;		/* FE_SPACE */
	/* The basis functions are monomials x^i * y^j * z^k, i + j + k <= o,
	 * they are order as (here o := type->order):
	 * 	for (k = 0; k <= o; k++)
	 * 	    for (j = 0; j <= o - k; j++)
	 * 		for (i = 0; i <= o - k - j; i++)
	 * We store the powers for the l-th basis function in type->orders[l]
	 * 	i + (j + k * (o + 1)) * (o + 1)
	 */
	type->orders = phgAlloc(type->nbas * sizeof(*type->orders));
	for (l = k = 0; k <= o; k++)
	    for (j = 0; j <= o - k; j++)
	 	for (i = 0; i <= o - k - j; i++)
		    type->orders[l++] = i + (j + k * (o + 1)) * (o + 1);
	assert(l == type->nbas);
    }
    DOF_DEFAULT = DOF_DGQn[1];
}

static void
DOF_TYPE_free(void)
{
    int i, n = DOF_TYPE_ORDER_MAX + 1;

    for (i = 0; i < n; i++) {
	phgDofTypeFree(DOF_DGQn + i);
	phgDofTypeFree(DOF_DGPn + i);
    }
    phgFree(DOF_DGQn);
    phgFree(DOF_DGPn);
    DOF_DGQn = DOF_DGPn = NULL;
    DOF_DEFAULT = NULL;
}

DOF_TYPE *
phgDofTypeVector(DOF_TYPE *base_type, int dim, const char *name)
/* returns a vector DOF_TYPE of dimension dim * "base_type"->dim
 * The optional argument "name" gives the name of the new type */
{
    DOF_TYPE temp = {
	NULL,				/* name */
	base_type->grad_type,		/* grad_type */
	base_type,			/* base_type */
	NULL,				/* points (or nodes, in ref coord) */
	NULL,				/* work (unused) */
	NULL,				/* orders (unused) */
	NULL,				/* pvt (unused) */
	dim * base_type->dim,		/* dim */
	dim * base_type->nbas,		/* nbas */
	base_type->order,		/* polynomial order of the basis */
	base_type->tensor_product,	/* whether tensor product */
	base_type->is_nodal,		/* whether tensor product */
	base_type->fe_space		/* finite element space */
    };
    DOF_TYPE *type = phgAlloc(sizeof(*type));

    memcpy(type, &temp, sizeof(temp));
    if (name != NULL)
	type->name = strdup(name);
    else
	sprintf(type->name = phgAlloc(strlen(base_type->name) + 1 + 1),
		"%sv", base_type->name);
    return type;
}

void
phgDofTypeFree(DOF_TYPE **ptype)
{
    DOF_TYPE *type;

    if (ptype == NULL)
	return;
    type = *ptype;
    *ptype = NULL;
    if (type == NULL)
	return;
    phgFree(type->name);
    phgFree(type->points);
    phgFree(type->work);
    phgFree(type->orders);
    phgFree(type->pvt);
    phgFree(type);
}

#define ElementIsLocal(g, e) ((e)->index < (g)->nelem && (e)->index >= 0)

GRID *
phgNewGrid(int flags)
{
    GRID *g;
    Unused(flags);
    g = phgCalloc(1, sizeof(*g));
    return g;
}

static void
log_handler(FILE *log_stream, const char *filename, int lineno,
            int package, int category, int priority, const char *msg)
{
    /* do nothing, just to get rid of the annoying messages from p4est */
}

static const char *
default_dof_handler(OPTION *o, const char *arg, int prefix)
{
    char s0[128];

    Unused(prefix);

    if (arg == NULL) {
	/* print help message */
	{
	    sprintf(s0, "  Available types: DGQ0...DGQ%d, DGP0...DGP%d",
			DOF_TYPE_ORDER_MAX, DOF_TYPE_ORDER_MAX);
	    phgOptionsPrintfHelp(o, s0);
	}
    }
    else if (arg[0] != '\0') {	/* Note the convention with arg[0]=='\0' */
	const char *s = arg;
	if (s[0] == 'D' && s[1] == 'G' && isdigit(arg[2])) {
	    /* use "DGQ" for "DG" */
	    strcpy(s0, "DGQ");
	    strcpy(s0 + strlen(s0), s + 2);
	    s = s0;
	}
	if (!strncmp(s, "DGQ", 3) && isdigit(s[3])) {
	    int o = atoi(s + 3);
	    if (o < 0 || o > DOF_TYPE_ORDER_MAX)
		return NULL;
	    DOF_DEFAULT = DOF_DGQn[o];
	}
	else if (!strncmp(s, "DGP", 3) && isdigit(s[3])) {
	    int o = atoi(s + 3);
	    if (o < 0 || o > DOF_TYPE_ORDER_MAX)
		return NULL;
	    DOF_DEFAULT = DOF_DGPn[o];
	}
	else {
	    return NULL;
	}
    }

    return DOF_DEFAULT->name;
}

static void
modify_options()
{
    phgOptionsModify("-dof_type", "Default DOF type", NULL, default_dof_handler,
		     VT_HANDLER, FALSE);
}

void
phgInit(int *argc, char ***argv)
{
#undef phgInit
    int i, j;
    extern void (*_phg_modify_options_callback)(void);	/* in option.c */

    DOF_TYPE_init();
    _phg_modify_options_callback = modify_options;
    phgInit(argc, argv);
    phgPrintf("PHG for hexahedral meshes based on PHG-to-p4est.\n");
    if (phgVerbosity > 0)
	phgPrintf("sizeof(ELEMENT) = %d\n", (int)sizeof(ELEMENT));
    sc_init(phgComm, 0, 0, NULL,
	    SC_LP_SILENT/*|SC_LP_PRODUCTION|SC_LP_ESSENTIAL*/);
    sc_set_log_defaults(NULL, log_handler, 0);

    /* initialize p8est_edge_no[][] */
    for (i = 0; i < NVert; i++) {
	p8est_edge_no[i][i] = -1;
	for (j = 0; j < i; j++) {
	    int k;
	    for (k = 0; k < NEdge; k++) {
		if ((GetEdgeVertex(k, 0) == i && GetEdgeVertex(k, 1) == j) ||
		    (GetEdgeVertex(k, 0) == j && GetEdgeVertex(k, 1) == i))
		    break;
	    }
	    p8est_edge_no[i][j] = p8est_edge_no[j][i] = k < NEdge ? k : -1;
	}
    }
}

void
phgFinalize(void)
{
    DOF_TYPE_free();
    sc_finalize();
#undef phgFinalize
    phgFinalize();
}

static FLOAT (*
get_corners(GRID *g, ELEMENT *e, FLOAT (*box)[Dim]))[Dim]
{
    static FLOAT corners[2][Dim];
#if USE_OMP
# pragma omp threadprivate(corners)
#endif
    p8est_quadrant_t node;
    int i;
    double box0[2][Dim];

    p8est_quadrant_corner_node(e->quad, 0, &node);
    p8est_qcoord_to_vertex(g->conn, e->treeid, node.x, node.y, node.z, box0[0]);

    p8est_quadrant_corner_node(e->quad, 7, &node);
    p8est_qcoord_to_vertex(g->conn, e->treeid, node.x, node.y, node.z, box0[1]);

    if (box == NULL)
	box = corners;
    for (i = 0; i < Dim; i++) {
	box[0][i] = box0[0][i];
	box[1][i] = box0[1][i];
    }

    return box;
}

static void
cb_elem(p8est_iter_volume_info_t *info, void *user_data)
{
    GRID *g = user_data;
    p8est_quadrant_t *quad = info->quad;
    ELEMENT *e = quad->p.user_data;

    e->elem_type = ET_CUBOID;
    e->index = g->i;
    e->tree = p8est_tree_array_index(g->p8est->trees, info->treeid);
    e->treeid = info->treeid;
    e->quad = quad;
    get_corners(g, e, e->corners);
    g->elems[g->i++] = e;
}

static void
update_elems(GRID *g)
{
    int i, nmax;

    g->nelem = g->p8est->local_num_quadrants;
    g->nelem_global = g->p8est->global_num_quadrants;

    phgFree(g->elems);
    g->elems = phgCalloc(g->nelem, sizeof(*g->elems));

    g->i = 0;
    p8est_iterate(g->p8est, g->ghost, /*user_data=*/g, cb_elem, NULL,NULL,NULL);
    assert(g->i == g->nelem);

    if (g->partition == NULL)
	g->partition = phgAlloc((g->nprocs + 1) * sizeof(*g->partition));
#if USE_MPI
    MPI_Allgather(&g->nelem, 1, PHG_MPI_INT, g->partition + 1, 1, PHG_MPI_INT,
		  g->comm);
#else	/* USE_MPI */
    g->partition[1] = g->nelem;
#endif	/* USE_MPI */
    g->partition[0] = 0;
    nmax = 0;
    for (i = 0; i < g->nprocs; i++) {
	if (nmax < g->partition[i + 1])
	    nmax = g->partition[i + 1];
	g->partition[i + 1] += g->partition[i];
    }
    g->lif = (FLOAT)nmax * (FLOAT)g->nprocs / (FLOAT)g->nelem_global;
}

# if USE_MPI
static void
setup_ghost_comm(GRID *g)
{
    int ii, rank, *rcnts, *rdsps, *scnts, *sdsps;
    INT i, nsend, nrecv, *sbuf, *rbuf;

    if (g->nprocs == 1 || (g->ghost_comm != NULL &&
			   g->ghost_comm->serial_no == g->serial_no))
	return;

    if (g->ghost_comm != NULL) {
	phgFree(g->ghost_comm->rcnts);
	phgFree(g->ghost_comm->rdsps);
	phgFree(g->ghost_comm->scnts);
	phgFree(g->ghost_comm->sdsps);
	phgFree(g->ghost_comm->slist);
	phgFree(g->ghost_comm);
    }
    g->ghost_comm = phgCalloc(1, sizeof(*g->ghost_comm));
    g->ghost_comm->serial_no = g->serial_no;

    rcnts = phgCalloc(g->nprocs, sizeof(*rcnts));
    rdsps = phgAlloc(g->nprocs * sizeof(*rdsps));
    scnts = phgAlloc(g->nprocs * sizeof(*scnts));
    sdsps = phgAlloc(g->nprocs * sizeof(*sdsps));

    /* count # entries from each process */
    if (g->nghost > 0) {
	assert(g->nghost < INT_MAX);
	rank = 0;
	for (i = 0; i < g->nghost; i++) {
	    INT k = g->O2Gmap[g->ordering[i]];
	    while (rank < g->nprocs && g->partition[rank + 1] <= k)
		rank++;
	    assert(rank < g->nprocs);
	    rcnts[rank]++;
	}
    }
    MPI_Alltoall(rcnts, 1, MPI_INT, scnts, 1, MPI_INT, g->comm);
    nsend = nrecv = 0;
    for (ii = 0; ii < g->nprocs; ii++) {
	assert(nsend < INT_MAX && nrecv < INT_MAX);
	rdsps[ii] = nrecv;
	sdsps[ii] = nsend;
	nrecv += rcnts[ii];
	nsend += scnts[ii];
    }

    MPI_Allreduce(&nsend, &i, 1, PHG_MPI_INT, MPI_MAX, g->comm);
    if (i == 0) {
	phgFree(rcnts);
	phgFree(rdsps);
	phgFree(scnts);
	phgFree(sdsps);
	g->ghost_comm->do_nothing = TRUE;
	return;
    }

    sbuf = phgAlloc(nsend * sizeof(*sbuf));
    rbuf = phgAlloc(nrecv * sizeof(*rbuf));
    for (i = 0; i < g->nghost; i++)
	rbuf[i] = g->O2Gmap[g->ordering[i]];

    MPI_Alltoallv(rbuf, rcnts, rdsps, PHG_MPI_INT,
		  sbuf, scnts, sdsps, PHG_MPI_INT, g->comm);
    phgFree(rbuf);
    for (i = 0; i < nsend; i++) {
	sbuf[i] -= g->partition[g->rank];
	assert(sbuf[i] >= 0 && sbuf[i] < g->nelem);
    }	

    g->ghost_comm->rcnts = rcnts;
    g->ghost_comm->rdsps = rdsps;
    g->ghost_comm->scnts = scnts;
    g->ghost_comm->sdsps = sdsps;
    g->ghost_comm->slist = sbuf;
    g->ghost_comm->nsend = nsend;
}
# endif	/* USE_MPI */

void
phgDofUpdateGhost(DOF *dof)
/* get DOF ghost data.
 * FIXME: work with multiple DOFs for better communication performance */
{
#if 0
    /* Use p4est ghost_exchange, see:
     *	p4est_ghost.c/p4est_ghost_exchange_data_begin */
    p8est_ghost_t *ghost = dof->g->ghost;
    void **mirror_data;
    int data_count = dof->dim * DofNBas(dof, NULL);	/* FLOATs/element */
    size_t data_size = data_count * sizeof(*dof->data);
    p8est_quadrant_t *mirror, *q;
    p8est_topidx_t which_tree;
    p8est_locidx_t which_quad;
    p8est_tree_t *tree;
    ELEMENT *e;
    INT i;

    phgFree(dof->ghost_data);
    dof->ghost_data = NULL;
    if (dof->g->ghost == NULL)
	return;

    dof->ghost_data = phgCalloc(dof->g->nghost, data_size);
    mirror_data = phgAlloc(ghost->mirrors.elem_count * sizeof(void *));

    for (i = 0; i < ghost->mirrors.elem_count; ++i) {
	mirror = p8est_quadrant_array_index(&ghost->mirrors, i);
	which_tree = mirror->p.piggy3.which_tree;
	tree = p8est_tree_array_index (dof->g->p8est->trees, which_tree);
	which_quad = mirror->p.piggy3.local_num - tree->quadrants_offset;
	q = p8est_quadrant_array_index (&tree->quadrants, which_quad);
	e = q->p.user_data;
	mirror_data[i] =dof->data + e->index * data_count;
    }
    p8est_ghost_exchange_custom(dof->g->p8est, ghost, data_size,
					mirror_data, dof->ghost_data);
    phgFree(mirror_data);
#else
#if USE_MPI
    GRID *g = dof->g;
    int data_count;	/* FLOATs/element */
    size_t data_size;
    FLOAT *sbuf, *rbuf;
    INT i, nsend, nrecv;
    MPI_Datatype type;

    if (g->nprocs == 1)
	return;

    if (SpecialDofType(dof->type) || dof->userfunc == DofNoData)
	return;

    data_count = dof->dim * DofNBas(dof, NULL);	/* FLOATs/element */
    if (data_count == 0)
	return;
    data_size = data_count * sizeof(*dof->data);

    if (g->ghost_comm == NULL || g->serial_no != g->ghost_comm->serial_no)
	setup_ghost_comm(g);

    if (g->ghost_comm->do_nothing)
	return;

    nsend = g->ghost_comm->nsend;
    nrecv = g->nghost;

    phgFree(dof->ghost_data);
    dof->ghost_data = phgCalloc(g->nghost, data_size);

    sbuf = phgAlloc((nsend + nrecv) * data_size);
    rbuf = sbuf + nsend * data_count;

    for (i = 0; i < nsend; i++)
	memcpy(sbuf + i * data_count,
	       dof->data + g->ghost_comm->slist[i] * data_count,
	       data_size);

    MPI_Type_contiguous(data_count, PHG_MPI_FLOAT, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(sbuf, g->ghost_comm->scnts, g->ghost_comm->sdsps, type,
		  rbuf, g->ghost_comm->rcnts, g->ghost_comm->rdsps, type,
		  g->comm);
    MPI_Type_free(&type);

    for (i = 0; i < nrecv; i++)
	memcpy(dof->ghost_data + g->ordering[i] * data_count,
	       rbuf  + i * data_count, data_size);

    phgFree(sbuf);
#endif	/* USE_MPI */
#endif
}

static void
remove_ghost(GRID *g)
{
    DOF **pdof;

    if (g->mesh != NULL)
	p8est_mesh_destroy(g->mesh);
    g->mesh = NULL;
    if (g->ghost != NULL)
	p8est_ghost_destroy (g->ghost);
    g->ghost = NULL;
    phgFree(g->ghost_data);
    g->ghost_data = NULL;
    phgFree(g->O2Gmap);
    g->O2Gmap = NULL;
    phgFree(g->ordering);
    g->ordering = NULL;
    g->nghost = 0;

    if (g->ghost_comm != NULL) {
	phgFree(g->ghost_comm->rcnts);
	phgFree(g->ghost_comm->rdsps);
	phgFree(g->ghost_comm->scnts);
	phgFree(g->ghost_comm->sdsps);
	phgFree(g->ghost_comm->slist);
	phgFree(g->ghost_comm);
	g->ghost_comm = NULL;
    }

    for (pdof = g->dofs; pdof != NULL && *pdof != NULL; pdof++) {
	phgFree((*pdof)->ghost_data);
	(*pdof)->ghost_data = NULL;
    }

    g->serial_no += 1000 * 1000;
}

static INT *O2Gmap = NULL;

static int
comp_O2Gmap(const void *p0, const void *p1)
{
    INT i = O2Gmap[*(const INT *)p0] - O2Gmap[*(const INT *)p1];
    return i < 0 ? -1 : i > 0 ? 1 : 0;
}

static void
setup_ghost(GRID *g, int halo_type)
{
    DOF **pdof;
    INT i;

    if (g->ghost != NULL) {
	if (g->halo_type <= halo_type)
	    return;
    }

    if (g->ghost_comm != NULL || g->ghost != NULL)
	remove_ghost(g);

    g->serial_no += 1000 * 1000;

    g->halo_type = halo_type;
    if (halo_type == HALO_FACE)
	halo_type = P8EST_CONNECT_FACE;
    else if (halo_type == HALO_EDGE)
	halo_type = P8EST_CONNECT_EDGE;
    else if (halo_type == HALO_VERT)
	halo_type = P8EST_CONNECT_CORNER;
    else
	phgError(1, "%s: unsupported halo type %d\n", __func__, halo_type);

    g->ghost = p8est_ghost_new(g->p8est, halo_type);
    g->ghost_data = phgAlloc(sizeof(ELEMENT) * g->ghost->ghosts.elem_count);

    /* FIXME: temporarily change e->index to global before ghost exchange */
    for (i = 0; i < g->nelem; i++)
	g->elems[i]->index += g->p8est->global_first_quadrant[g->rank];
    p8est_ghost_exchange_data(g->p8est, g->ghost, g->ghost_data);
    /* restore e->index to local */
    for (i = 0; i < g->nelem; i++)
	g->elems[i]->index -= g->p8est->global_first_quadrant[g->rank];

    g->mesh = p8est_mesh_new/*_ext*/(g->p8est, g->ghost, halo_type);

    g->nghost = g->ghost->ghosts.elem_count;
    g->elems = phgRealloc0(g->elems, (g->nghost+g->nelem) * sizeof(*g->elems));
    g->O2Gmap = phgCalloc(g->nghost, sizeof(*g->O2Gmap));
    /* update user_data members of ghost elements */
    for (i = 0; i < g->nghost; i++) {
	p8est_quadrant_t *q = p8est_quadrant_array_index(&g->ghost->ghosts, i);
	ELEMENT *e = ((ELEMENT *)g->ghost_data) + i;
	e->quad = q;
	e->tree = NULL;
	g->O2Gmap[i] = e->index;
	e->index = g->nelem + i;
	g->elems[e->index] = e;
    }

    if (g->nghost > 0) {
	g->ordering = phgAlloc(g->nghost * sizeof(*g->ordering));
	for (i = 0; i < g->nghost; i++)
	    g->ordering[i] = i;
	assert(g->nghost <= INT_MAX);
	O2Gmap = g->O2Gmap;
	qsort(g->ordering, g->nghost, sizeof(*g->ordering), comp_O2Gmap);
    }

    for (pdof = g->dofs; pdof != NULL && *pdof != NULL; pdof++)
	phgDofUpdateGhost(*pdof);
}

void
phgFreeGrid(GRID **gptr)
{
    GRID *g;
    DOF *dof;

    if (gptr == NULL || *gptr == NULL)
	return;
    g = *gptr;
    *gptr = NULL;

    remove_ghost(g);

    if (g->p8est != NULL)
	p8est_destroy(g->p8est);
    if (g->conn != NULL)
	p8est_connectivity_destroy(g->conn);
    phgFree(g->elems);
    phgFree(g->roots);
    phgFree(g->partition);
    while (g->dofs != NULL && *g->dofs != NULL) {
	dof = *g->dofs;
	if (g->rank == 0)
	    phgPrintf("*** WARNING: unfreed DOF \"%s\" (created at %s:%d)\n",
			dof->name, dof->srcfile, dof->srcline);
	phgDofFree(&dof);
    }
    phgFree(g->dofs);
    phgFree(g);
}

static const char *
get_token(FILE *fp)
/* reads and returns the next token from the file "fp", returns NULL if EOF.
 * "fp" == NULL ==> reset this function and return NULL */
{
    static char *token, line[1024];

    if (fp == NULL) {
	token = NULL;
	return NULL;
    }

    while (TRUE) {
	char *p;
	if (token == NULL) {
	    if (fgets(line, sizeof(line), fp) == NULL)
		return NULL;
	    if ((token = strchr(line, '#')) != NULL)
		*token = '\0';
	    token = line;
	}
	while (isspace(*token))
	    token++;
	if (*token == '\0') {
	    token = NULL;
	    continue;	/* read next line */
        }
	p = token;
	while (!isspace(*token) && *token != '\0')
	    token++;
	if (*token == '\0')
	    token = NULL;
	else
	    *(token++) = '\0';
	return p;
    }

    /*return NULL;*/	/* make gcc happy */
}

static void
read_medit_file(const char *filename, INT *nvert, INT *nquad, INT *nhexa,
		FLOAT (**verts)[Dim + 1], INT (**quads)[5], INT (**hexa)[9])
/* Reads an Medit mesh file, returns TRUE if success or FALSE if failure.
 *
 * The number of vertices, quadrilaterals and hexahedra are respectively
 * returned in "nvert", "nquad" and "nhexa".
 *
 * The list of vertices, quadrilaterals and hexahedra are respectively
 * returned in dynamically allocated arrays "verts[]", "quads[]" and "hexa[]",
 * which are to be freed by the calling function. The last number of each entry
 * is set to the tag of the vertex (casted to FLOAT) or element.
 */
{
    FILE *fp;
    const char *token;
    INT i, n;
    int j;

#define MEDIT_ERROR \
	phgError(1, "%s:%d error reading mesh file \"%s\".\n", \
			__FILE__, __LINE__, filename);

    get_token(NULL);	/* reset get_token() */

    if ((fp = phgOpenInputFile_(filename)) == NULL ||
	(token = get_token(fp)) == NULL ||
	strcasecmp(token, "MeshVersionFormatted") ||
	(token = get_token(fp)) == NULL || strcmp(token, "1") ||
	(token = get_token(fp)) == NULL || strcasecmp(token, "Dimension") ||
	(token = get_token(fp)) == NULL || strcmp(token, "3"))
	MEDIT_ERROR;

    *nvert = *nquad = *nhexa = 0;
    *verts = NULL;
    *quads = NULL;
    *hexa  = NULL;

    while (TRUE) {
	if ((token = get_token(fp)) == NULL || !strcasecmp(token, "End"))
	    break;
	if (!strcasecmp(token, "Vertices")) {
	    if ((token = get_token(fp)) == NULL)
		MEDIT_ERROR;
	    assert(*verts == NULL);
	    *nvert = n = (INT)atoll(token);
	    *verts = phgAlloc(n * sizeof((*verts)[0]));
	    for (i = 0; i < n; i++) {
		for (j = 0; j < Dim + 1; j++) {
		    if ((token = get_token(fp)) == NULL)
			MEDIT_ERROR;
		    (*verts)[i][j] = atof(token);
		}
	    }
	}
	else if (!strcasecmp(token, "Quadrilaterals")) {
	    if ((token = get_token(fp)) == NULL)
		MEDIT_ERROR;
	    assert(*quads == NULL);
	    *nquad = n = (INT)atoll(token);
	    *quads = phgAlloc((n + 1) * sizeof((*quads)[0]));
	    for (i = 0; i < n; i++) {
		for (j = 0; j < 4 + 1; j++) {
		    if ((token = get_token(fp)) == NULL)
			MEDIT_ERROR;
		    (*quads)[i][j] = (INT)atoll(token);
		}
	    }
	}
	else if (!strcasecmp(token, "Hexahedra")) {
	    if ((token = get_token(fp)) == NULL)
		MEDIT_ERROR;
	    assert(*hexa == NULL);
	    *nhexa = n = (INT)atoll(token);
	    *hexa = phgAlloc(n * sizeof((*hexa)[0]));
	    for (i = 0; i < n; i++) {
		for (j = 0; j < 8 + 1; j++) {
		    if ((token = get_token(fp)) == NULL)
			MEDIT_ERROR;
		    (*hexa)[i][j] = (INT)atoll(token);
		}
	    }
	}
    }

    phgCloseInputFile_(fp);
}

/* The following struct is used to align vertices of an element with those of
 * the standard unit cube */
typedef struct {
    FLOAT vertex[Dim];		/* vertex coordinates */
    int	  index;		/* vertex local numbering */
} VERTEX_t;

static int
comp_vertex(const void *p0, const void *p1)
{
    const FLOAT *v0 = ((const VERTEX_t *)p0)->vertex,
		*v1 = ((const VERTEX_t *)p1)->vertex;
    int i;

    for (i = 0; i < 3; i++) {
	FLOAT d;
	if (Fabs(d = v0[i] - v1[i]) > 100.0 * FLOAT_EPSILON)
	    return d < 0.0 ? -1 : 1;
    }

    return 0;
}

/* the following variables are used in comp_quad() */
static INT (*quads)[5] = NULL;

static int
comp_quad(const void *p0, const void *p1)
/* note: the vertex indices are assumed sorted in both "p0" and "p1" */
{
    INT *q0 = quads[*(const INT *)p0], *q1 = quads[*(const INT *)p1];
    int i;

    for (i = 0; i < 3; i++)
	if (q0[i] != q1[i])
	    return q0[i] - q1[i];

    return q0[i] - q1[i];
}

static void
init0_fn(p8est_t *p8est, p8est_topidx_t which_tree, p8est_quadrant_t *quad)
/* init_fn() for p8est_new() */
{
    ELEMENT *e = quad->p.user_data;

    bzero(e, sizeof(*e));
}

BOOLEAN
phgImport_(GRID *g, const char *mesh_fn, MPI_Comm comm)
{
    INT i, nvert = 0, nquad = 0, nhexa = 0;
    INT (*hexa)[9] = NULL, *quad_ordering = NULL;
    BTYPE default_btype = UNDEFINED;
    FLOAT (*verts)[Dim + 1] = NULL;
    int ii;
    ELEMENT *e;

    g->comm = comm;
#if USE_MPI
    MPI_Comm_size(comm, &g->nprocs);
    MPI_Comm_rank(comm, &g->rank);
#else	/* USE_MPI */
    g->nprocs = 1;
    g->rank = 0;
#endif	/* USE_MPI */

    if (mesh_fn == NULL) {
	/* test with a built-in unit cube (used to create a dummy GRID) */
	g->conn = p8est_connectivity_new_brick(1, 1, 1, 0, 0, 0);
	g->p8est = p8est_new_ext(g->comm, g->conn, 0, 0, 1, sizeof(ELEMENT), 
				init0_fn, NULL);
	g->balanced = TRUE;
	g->bbox[0][0] = g->bbox[0][1] = g->bbox[0][2] = 0.0;
	g->bbox[1][0] = g->bbox[1][1] = g->bbox[1][2] = 1.0;
	update_elems(g);
	return TRUE;
    }

    /* Seems all processes need the mesh data. Shall we read the same file
     * in all processes, or just read it in process 0 then broadcast? */
    /*if (g->rank == 0)*/ {
	read_medit_file(mesh_fn, &nvert, &nquad, &nhexa, &verts, &quads, &hexa);
	if (phgVerbosity > 0)
	    phgPrintf("Loading mesh \"%s\":\n  nvert=%d, nquad=%d, nhexa=%d\n",
			mesh_fn, nvert, nquad, nhexa);
	if (nhexa == 0)
	    phgError(1, "the file \"%s\" contains no hexahedra.\n", mesh_fn);
    }
#if 0*USE_MPI
    if (g->nprocs > 1) {
	INT a[] = {nvert, nhexa};
	MPI_Bcast(a, 2, PHG_MPI_INT, 0, g->comm);
	nvert = a[0];
	nhexa = a[1];
    }
#endif	/* USE_MPI */
    g->conn = p8est_connectivity_new(nvert, nhexa, 0, 0, 0, 0);
    for (i = 0; i < nvert; i++)
	for (ii = 0; ii < Dim; ii++)
	    g->conn->vertices[Dim * i + ii] = verts[i][ii];
    for (i = 0; i < nhexa; i++) {
	p4est_topidx_t *v = g->conn->tree_to_vertex + 8 * i;
#if 0	/******************* for general hexhedral mesh ********************/
	for (ii = 0; ii < 8; ii++) {
	    /* right-hand vertex order -> z-order */
	    static int order[] = {0, 1, 3, 2, 4, 5, 7, 6};
	    v[order[ii]] = hexa[i][ii] - 1;
	}
#else	/************************ for cuboid mesh **************************/
	/* Check the element is a cuboid and reorder the vertices so that v0
	 * and v7 have respectively the smallest and largest coordinates. */

	/* First, get the two corners with smallest and largest coordinates and
	 * store them in box[0] and box[1] */
	FLOAT box[2][Dim];
	for (ii = 0; ii < 8; ii++) {
	    FLOAT *t = verts[hexa[i][ii] - 1];
	    int jj;
	    for (jj = 0; jj < Dim; jj++) {
		if (ii == 0 || box[0][jj] > t[jj])
		    box[0][jj] = t[jj];
		if (ii == 0 || box[1][jj] < t[jj])
		    box[1][jj] = t[jj];
	    }
	}

	/* Then match (align) the vertices of the current element with the
	 * reference element (assumed by functions like phgGeomGetCorners) */

	/* "ref[]" stores the vertices and canonical ordering of the unit cube,
	 * sorted according to the coordinates of the vertices. */
	VERTEX_t ref[] = {
	    /* Note: data for "ref" are generated with the bash command below:
-------------------------------------------------------------------------------
cat << END | sort | \
awk '{print "\t    {{box["$1"][0],box["$2"][1],box["$3"][2]},"$4"},"}'
0 0 0 0
1 0 0 1
0 1 0 2
1 1 0 3
0 0 1 4
1 0 1 5
0 1 1 6
1 1 1 7
END
-----------------------------------------------------------------------------*/
	    {{box[0][0],box[0][1],box[0][2]},0},
	    {{box[0][0],box[0][1],box[1][2]},4},
	    {{box[0][0],box[1][1],box[0][2]},2},
	    {{box[0][0],box[1][1],box[1][2]},6},
	    {{box[1][0],box[0][1],box[0][2]},1},
	    {{box[1][0],box[0][1],box[1][2]},5},
	    {{box[1][0],box[1][1],box[0][2]},3},
	    {{box[1][0],box[1][1],box[1][2]},7},
	}, key, *found;
	for (ii = 0; ii < 8; ii++) {
	    memcpy(key.vertex, verts[hexa[i][ii] - 1], sizeof(key.vertex));
	    found = bsearch(&key, ref, 8, sizeof(key),comp_vertex);
	    if (found == NULL)
		phgError(1, "the file \"%s\" contains non cuboid elements.\n",
				 mesh_fn);
	    v[found->index] = hexa[i][ii] - 1;
	}
#endif	/*******************************************************************/

	for (ii = 0; ii < 6; ii++) {
	    g->conn->tree_to_tree[6 * i + ii] = i;
	    g->conn->tree_to_face[6 * i + ii] = ii;
	}
    }

    phgFree(verts);

#if DEBUG
    P4EST_ASSERT(p8est_connectivity_is_valid(g->conn));
#endif	/* DEBUG */
    p8est_connectivity_complete(g->conn);

    g->p8est = p8est_new_ext(g->comm, g->conn, 0, 0, 1, sizeof(ELEMENT), 
				init0_fn, NULL);
    update_elems(g);
    g->balanced = TRUE;
    setup_ghost(g, HALO_VERT);

    /* sort indices in quads[*] and set up quad_ordering[] */
    if (nquad > 0) {
	quad_ordering = phgAlloc(nquad * sizeof(*quad_ordering));
	for (i = 0; i < nquad; i++) {
	    /* sort indices in quads[i] */
	    qsort(quads[i], 4, sizeof(INT), phgCompINT);
	    quad_ordering[i] = i;
	}
	qsort(quad_ordering, nquad, sizeof(*quad_ordering), comp_quad);
	phgLoadBCmap(mesh_fn);
    }

    /* get default boundary type */
    default_btype = phgBTypesValue(phgOptionsGetHandler("-default_bdry_type"));
    assert(default_btype != -1);

    /* set region_mark and boundary types */
    ForAllElements(g, e) {
	int jj;
	INT key = nquad, *found;
	p4est_topidx_t *verts = g->conn->tree_to_vertex + 8 * e->treeid;

	assert(e->treeid < nhexa);
	e->region_mark = hexa[e->treeid][8];

	/* set boundary type of boundary faces to default_btype */
	for (ii = 0; ii < 6; ii++) {
	    if (phgGetNeighbour(g, e, ii) != NULL) {
		e->bound_type[ii] = INTERIOR;		/* interior face */
		continue;
	    }

	    e->bound_type[ii] = default_btype;	/* boundary face */
	    if (nquad == 0)
		continue;
	    /* get the vertices of the 'ii'-th face */
	    for (jj = 0; jj < 4; jj++)
		quads[key][jj] = verts[p8est_face_corners[ii][jj]] + 1;
	    qsort(quads[key], 4, sizeof(INT), phgCompINT);
	    /* search the list of quadrilaterals from the input mesh */
	    found = bsearch(&key, quad_ordering, nquad, sizeof(key), comp_quad);
	    if (found)
		e->bound_type[ii] = phgMapBC(quads[*found][4]);
	}
    }

    phgFree(hexa);
    phgFree(quads);
    phgFree(quad_ordering);

    /* determine BBox(g) */
    i = 0;
    ForAllElements_(g, e, FOR_ALL) {
	FLOAT (*corners)[Dim] = e->corners;
	if (i == 0) {
	    memcpy(g->bbox, e->corners, sizeof(g->bbox));
	    i = 1;
	    continue;
	}
	for (ii = 0; ii < Dim; ii++) {
	    if (g->bbox[0][ii] > corners[0][ii])
		g->bbox[0][ii] = corners[0][ii];
	    if (g->bbox[1][ii] < corners[1][ii])
		g->bbox[1][ii] = corners[1][ii];
	}
    }
#if USE_MPI
    if (g->nprocs > 1) {
	FLOAT bbox0[2][Dim], bbox1[2][Dim];
	for (ii = 0; ii < Dim; ii++) {
	    bbox0[0][ii] = -g->bbox[0][ii];
	    bbox0[1][ii] = g->bbox[1][ii];
	}
	MPI_Allreduce(bbox0, bbox1, 2 * Dim, PHG_MPI_FLOAT, MPI_MAX, g->comm);
	for (ii = 0; ii < Dim; ii++) {
	    g->bbox[0][ii] = -bbox1[0][ii];
	    g->bbox[1][ii] = bbox1[1][ii];
	}
    }
#endif	/* USE_MPI */

    return TRUE;
}

void
phgSetupHalo_(GRID *g, int halo_type, BOOLEAN update_dofs)
{
    if (g->conn == NULL)
	return;

    if (g->nprocs == 1 && halo_type == HALO_ALL)
	halo_type = HALO_VERT;

    setup_ghost(g, halo_type);

#if 0	/* for debugging */
#warning debugging code!
    phgInfo(-1, "%d local elements:\n", g->nelem);
    for (INT i = 0; i < g->nelem; i++) {
	ELEMENT *e = g->elems[i];
	phgInfo(-1, "elem %d: box=(%lg,%lg,%lg)-(%lg,%lg,%lg)\n",
			GlobalElement(g, e->index), 
			e->corners[0][0], e->corners[0][1], e->corners[0][2],
			e->corners[1][0], e->corners[1][1], e->corners[1][2]);
	for (int j = 0; j < NFace; j++) {
	    ELEMENT *e1 = phgGetNeighbour(g, e, j);
	    if (e1 == NULL || e1->index < 0)
		continue;
	    phgInfo(-1, "  neighbor[%d]=%d: box=(%lg,%lg,%lg)-(%lg,%lg,%lg)\n",
			j, GlobalElement(g, e1->index),
			e1->corners[0][0], e1->corners[0][1], e1->corners[0][2],
			e1->corners[1][0], e1->corners[1][1], e1->corners[1][2]);
	}
    }
#endif
}

size_t
phgMemoryUsage(GRID *g, size_t *peak)
{
    /* FIXME: modify phgMemoryUsage1 to use MPI_Comm instead of GRID */
#ifdef phgMemoryUsage
# undef phgMemoryUsage	/* we need the original phgMemoryUsage of PHG here */
#endif
    return phgMemoryUsage(NULL, peak);
}

FLOAT *
phgGeomXYZ2Lambda(GRID *g, ELEMENT *e, FLOAT x, FLOAT y, FLOAT z, FLOAT xref[])
{
    FLOAT xyz[] = {x, y, z};

    xyz_to_ref(g, e->index, xyz, xref);
    return xref;
}

void
phgGeomLambda2XYZ(GRID *g ,ELEMENT *e, const FLOAT xref[],
					FLOAT *x, FLOAT *y, FLOAT *z)
{
    FLOAT xyz[Dim];

    ref_to_xyz(g, e->index, xref, xyz);
    *x = xyz[0];
    *y = xyz[1];
    *z = xyz[2];
}

void
phgGeomGetFaceOutNormal(GRID *g, ELEMENT *e, int face, FLOAT *nv)
/* Note: the faces are ordered as -x, +x, -y, +y, -z, +z */
{
    switch (face) {
	case 0:	nv[0] =-1.; nv[1] = 0.; nv[2] = 0.; break;
	case 1: nv[0] = 1.; nv[1] = 0.; nv[2] = 0.; break;
	case 2: nv[0] = 0.; nv[1] =-1.; nv[2] = 0.; break;
	case 3: nv[0] = 0.; nv[1] = 1.; nv[2] = 0.; break;
	case 4: nv[0] = 0.; nv[1] = 0.; nv[2] =-1.; break;
	case 5: nv[0] = 0.; nv[1] = 0.; nv[2] = 1.; break;
	default: phgError(1, "%s: invalid face no. (%d)\n", __func__, face);
    }
}

FLOAT (*
phgGeomGetCorners(GRID *g, ELEMENT *e, FLOAT (*cub)[Dim]))[Dim]
{
    static COORD cuboid[2];
#if USE_OMP
# pragma omp threadprivate(cuboid)
#endif

    if (cub == NULL)
	cub = cuboid;

    memcpy(cub, e->corners, sizeof(e->corners));

    return cub;
}

ELEMENT *
phgGetNeighbour_(GRID *g, ELEMENT *e, int face, const char *file, int line)
/* Note: the faces are ordered as -x, +x, -y, +y, -z, +z.
 *
 * See p4est_mesh_face_neighbor_next() (src/p4est_mesh.c) in the source files
 * of p4est for information about related structs */
{
    int q2f;
    p8est_locidx_t q2q, quadfacecode;
    ELEMENT *e1;

    /* e must be local! */
    if (!ElementIsLocal(g, e))
	phgError(1, "%s:%d: can't access neighbours of halo elements.\n",
		 file, line);

    assert(g->ghost != NULL);
    assert(face >= 0 && face < P8EST_FACES);
    quadfacecode = P8EST_FACES * (e->tree->quadrants_offset +
		   	sc_array_position(&e->tree->quadrants, e->quad)) + face;
    q2q = g->mesh->quad_to_quad[quadfacecode];
    q2f = (int)g->mesh->quad_to_face[quadfacecode];
    if (q2f < 0) {
	/* neighbour is half size, return a dummy element with max generation
	 * and index -1 */
	static ELEMENT dummy = {{{0.,0.,0.},{0.,0.,0.}},/* corners */
				NULL,			/* tree */
				NULL,			/* quad */
				NULL, NULL,		/* reserved, user_ctx */
				-1,			/* index */
				-1,			/* treeid */
				0, 0,			/* mark, region_mark */
				{0,0,0,0,0,0},		/* bound_type */
				127/*P8EST_MAXLEVEL*/,	/* generation */
				ET_CUBOID		/* elem_type */
			       };
	return &dummy;
    }
    if (q2q < g->mesh->local_num_quadrants) {
	p8est_quadrant_t *q;
	q = p8est_mesh_quadrant_cumulative(g->p8est, q2q, NULL, NULL);
	e1 = q->p.user_data;
    }
    else {
	q2q -= g->mesh->local_num_quadrants;
	assert(q2q < g->mesh->ghost_num_quadrants);
	e1 = ((ELEMENT *)g->ghost_data) + q2q;
    }
    return e1 == e ? NULL : e1;
}

BYTE
phgOppositeFace(GRID *g, ELEMENT *e, BYTE v, ELEMENT *e_op)
{
    int q2f;
    p8est_locidx_t quadfacecode;

    /* The code below only allows e to be local! */
    assert(ElementIsLocal(g, e));
    assert(v >= 0 && v < P8EST_FACES);
    quadfacecode = P8EST_FACES * (e->tree->quadrants_offset +
		   	sc_array_position(&e->tree->quadrants, e->quad)) + v;
    q2f = g->mesh->quad_to_face[quadfacecode];
    assert(q2f >= 0);
    return q2f % 6;
}

DOF *
phgDofNew_(GRID *g, DOF_TYPE *type, HP_TYPE *hp, SHORT dim,
	   const char *name, DOF_USER_FUNC userfunc,
	   const char *srcfile, int srcline)
{
    DOF *dof = phgCalloc(1, sizeof(*dof));
    int i;

    assert(hp == NULL);
    dof->g = g;
    dof->dim = dim;
    dof->type = type;
    dof->name = name == NULL ? NULL : strdup(name);
    dof->userfunc = userfunc;
    dof->userfunc_lambda = NULL;
    dof->srcfile = srcfile;
    dof->srcline = srcline;
    dof->poly_order = -1;

    if (dof->type != DOF_ANALYTIC) {
	int nbas = DofNBas(dof, NULL);
	if (userfunc != DofNoData)
	    dof->data = phgCalloc(g->nelem, nbas * dim * sizeof(FLOAT));
	dof->poly_order = DofTypeOrder(dof, NULL);
	if (!SpecialDofUserFunction(userfunc))
	    phgDofSetDataByFunction(dof, userfunc);
    }

    /* attache "dof" to g->dofs[] */
    for (i = 0; g->dofs != NULL && g->dofs[i] != NULL; i++);
    g->dofs = phgRealloc0(g->dofs, (i + 2) * sizeof(*g->dofs));
    g->dofs[i] = dof;
    g->dofs[i + 1] = NULL;

    return dof;
}

/* context used by the callback functions */
typedef struct {
    int i1, i2;			/* ints */
    void *p1, *p2, *p3, *p4;	/* pointers */
} CTX;
static CTX *ctx = NULL;
#if USE_OMP
# pragma omp threadprivate(ctx)
#endif	/* USE_OMP */

static void
f_wrapper(FLOAT x, FLOAT y, FLOAT z, FLOAT *res)
{
    if (ctx->i1 == 0) {
	/* DOF_USER_FUNC */
	((DOF_USER_FUNC)ctx->p1)(x, y, z, res);
    }
    else {
	/* DOF_USER_FUNC_LAMBDA */
	FLOAT xref[3], xyz[3] = {x, y, z};
	DOF *u = ctx->p2;
	ELEMENT *e = ctx->p3;
	xyz_to_ref(u->g, e->index, xyz, xref);
	((DOF_USER_FUNC_LAMBDA)ctx->p1)(u, e, ctx->i2, xref, res); 
    }
}

void
phgDofSetDataByFunction__(DOF *dof, DOF_USER_FUNC func,
		DOF_USER_FUNC_LAMBDA func_lam, FLOAT *funcvalues,
		INT nelem, ELEMENT **elems)
/* compute DOFs by L2 projection or interpolation in each element */
{
    GRID *g = dof->g;
    ELEMENT *e;
    DOF_TYPE *type = dof->type;
    FLOAT *data;
    INT i;
    QCACHE *qc = NULL;
    int Q_f = 0;

    assert(funcvalues == NULL);		/* unimplemented */

    if (nelem < 0) {
	elems = g->elems;
	nelem = g->nelem;
    }

    if (func == NULL && func_lam == NULL) {
	int n = DofNBas(dof, NULL) * dof->dim;
	size_t size = n * sizeof(FLOAT);
	for (i = 0; i < nelem; i++) {
	    e = elems[i];
	    memset(dof->data + e->index * n, 0, size);
	}
	phgDofUpdateGhost(dof);
	return;
    }

    if (type->base_type != NULL) {
	/* Note: should not create a new DOF here when this function is called
	 * by replace_fn(), because it would alter g->dofs which is used in
	 * replace_fn() */
	dof->type = type->base_type;
	dof->dim *= type->dim / type->base_type->dim,
	phgDofSetDataByFunction__(dof, func,func_lam, funcvalues, nelem, elems);
	dof->dim /= type->dim / type->base_type->dim,
	dof->type = type;
	return;
    }

    /* Note: implement phgQCIntegrate functions returning multiple FLOATs */
#if USE_OMP
# pragma omp parallel
#endif	/* USE_OMP */
    {
	ctx = phgCalloc(1, sizeof(*ctx));
	ctx->i1 = func != NULL ? 0 : 1;	/* 0: FUNC, 1: FUNC_LAMBDA */
	ctx->p1 = func != NULL ? (void *)func : (void *)func_lam;
	ctx->p2 = dof;
    }

    if (!type->is_nodal) {
	qc = phgQCNew(QD_P4EST, dof);
	Q_f = phgQCAddXYZFunction(qc, f_wrapper, dof->dim);
    }

    for (i = 0; i < nelem; i++) {
	int j, n;
	e = elems[i];
	n = DofNBas(dof, e);
	if (type->is_nodal) {	/* Interpolation */
	    FLOAT x[Dim], tmp[Dim], *pt = tmp;
	    assert(type->tensor_product);	/* TODO otherwise */
	    data = dof->data + e->index * n * dof->dim;
	    for (j = 0; j < n; j++, data += dof->dim) {
		if (type->orders == NULL) {
		    pt = type->points + j * Dim;
		}
		else  {	/* j = jx + (jy + jz * (orders[1]+1)) * (orders[0]+1) */
		    int jx, jy, jz;
		    jx = j % (type->orders[0] + 1);
		    jz = j / (type->orders[0] + 1);
		    jy = jz % (type->orders[1] + 1);
		    jz /= (type->orders[1] + 1);
		    tmp[0] = type->points[jx];
		    tmp[1] = type->points[type->orders[0] + 1 + jy];
		    tmp[2] = type->points[type->orders[0] + 1 +
					  type->orders[1] + 1 + jz];
		}
		if (func == NULL) {
		    func_lam(dof, e, j, pt, data);
		    continue;
		}
		ref_to_xyz(dof->g, e->index, pt, x);
		func(x[0], x[1], x[2], data);
	    }
	}
	else {					/* elementwise L2 projection */
	    FLOAT (*corners)[Dim] = phgGeomGetCorners(g, e, NULL);
	    FLOAT d = 1.0 / phgGeomGetVolume(g, e);
	    int k;
	    /* FIXME: use projection on the reference element, save LU of A */
	    FLOAT *rule;
	    assert(!type->tensor_product);	/* TODO otherwise */
	    phgQuadInterfaceCuboid(NULL, 0, NULL, corners, type->order * 2,
				   &rule, NULL, NULL);
	    phgQCSetRule(qc, rule, -1.0);
	    ctx->p3 = e;
	    if (type->work == NULL) {
		FLOAT (*A)[n];
		/* set type->work = LU factorization of the mass matrix
		 * for the reference element (the unit cube) */
		type->work = phgAlloc(n * sizeof(*type->work) * n);
		type->pvt = phgAlloc(n * sizeof(*type->pvt));
		A = (void *)type->work;
		for (j = 0; j < n; j++)
		    for (k = 0; k <= j; k++)
			A[j][k] = A[k][j] = d * phgQCIntegrate
			    (qc, e->index, Q_BAS, j, qc, e->index, Q_BAS, k);
#if 0
#warning
		printf("A:\n");
		for (j = 0; j < n; j++) {
		    for (k = 0; k < n; k++)
			printf(" %lg", A[k][j]);
		    printf("\n");
		}
#endif
		phgSolverDenseLU(n, type->work, type->pvt);
	    }
	    data = dof->data + e->index * n * dof->dim;
	    for (j = 0; j < n; j++) {
		ctx->i2 = j;	/* pass bno to func_lam */
/* FIXME: without one of the two lines below interface-p4est segfaults.
 *			10.3.1 20210422 (Red Hat 10.3.1-1) (GCC)
 *  2020.05.05: it seems the bug doesn't show up if QCINTEGRATE_HACK is defined
 *		(see quad-cache.h) */
#ifndef QCINTEGRATE_HACK
//bzero(data + j * dof->dim, sizeof(*data) * dof->dim);
data[j * dof->dim] = *
#endif
		phgQCIntegrateM(qc, e->index, Q_BAS, j, qc, e->index, Q_f, 0,
				/*N*/1, /*M*/dof->dim, /*K*/1,
				data + j * dof->dim);
		for (k = 0; k < dof->dim; k++)
		    data[j * dof->dim + k] *= d;
	    }
	    phgFree(rule);
#if 0
#warning
if (phgVerbosity) {
    printf("b:\n");
    for (j = 0; j < n; j++) {
	for (k = 0; k < dof->dim; k++)
	    printf(" %lg", data[j * dof->dim + k]);
	printf("\n");
    }
}
#endif
	    phgSolverDenseSV(n, type->work, type->pvt, dof->dim, data);
#if 0
#warning
if (phgVerbosity) {
    printf("solution:\n");
    for (j = 0; j < n; j++) {
	for (k = 0; k < dof->dim; k++)
	    printf(" %lg", data[j * dof->dim + k]);
	printf("\n");
    }
}
#endif
	    data += n * dof->dim;
	}
    }
#if USE_OMP
# pragma omp parallel
#endif	/* USE_OMP */
    {
	phgFree(ctx);
    }

    phgQCFree(&qc);
    phgDofUpdateGhost(dof);
}

void
phgDofSetPolyOrder(DOF *u, SHORT order)
{
    u->poly_order = order;
}

void
phgDofSetFunction(DOF *u, DOF_USER_FUNC func)
{
    assert(u->type != DOF_CONSTANT);
    u->userfunc = func;
    u->userfunc_lambda = NULL;
    if (!SpecialDofUserFunction(func) && u->data != NULL)
	phgDofSetDataByFunction(u, func);
}

void
phgDofSetLambdaFunction(DOF *u, DOF_USER_FUNC_LAMBDA func)
{
    assert(u->type != DOF_CONSTANT);
    u->userfunc = DofLambdaFunction;
    u->userfunc_lambda = func;
    if (!SpecialDofType(u->type) && u->data != NULL)
	phgDofSetDataByLambdaFunction(u, func);
}

void
phgDofFree(DOF **pdof)
{
    DOF *dof, **pd;
    GRID *g;

    if (pdof == NULL || *pdof == NULL)
	return;
    dof = *pdof;
    g = dof->g;
    assert(g->dofs != NULL);

    phgFree(dof->name);
    phgFree(dof->data);
    phgFree(dof->ghost_data);
    phgFree(dof);

    /* remove "dof" from g->dofs[] */
    for (pd = g->dofs; *pd != dof && *pd != NULL; pd++);
    assert(*pd == dof);
    *pdof = NULL;	/* note: allow the case pdof == pd! */
    do {
	pd[0] = pd[1];
	pd++;
    } while (*pd != NULL);
    if (g->dofs[0] == NULL) {
	/* no DOF */
	phgFree(g->dofs);
	g->dofs = NULL;
    }
}

FLOAT *
phgDofEval(DOF *dof, ELEMENT *e, const FLOAT xref[], FLOAT *values)
{
    int i, j, k, nb;
    FLOAT bas[dof->type->dim], *data;

    if (dof->type == DOF_ANALYTIC) {
	if (dof->userfunc == DofLambdaFunction) {
	    dof->userfunc_lambda(dof, e, -1, xref, values);
	}
	else {
	    FLOAT x[3];
	    assert(!SpecialDofUserFunction(dof->userfunc));
	    ref_to_xyz(dof->g, e->index, (FLOAT *)xref, x);
	    dof->userfunc(x[0], x[1], x[2], values);
	}
	return values;
    }

    nb = DofNBas(dof, e);
    if (e->index < dof->g->nelem)
	data = dof->data + dof->dim * (size_t)e->index * nb;
    else
	data = dof->ghost_data +
			dof->dim * (size_t)(e->index - dof->g->nelem) * nb;
    for (i = 0; i < dof->dim * dof->type->dim; i++)
	values[i] = 0.;
    for (j = 0; j < nb; j++) {
	basisFunc(dof, xref[0], xref[1], xref[2], j, bas);
	for (i = 0; i < dof->dim; i++, data++)
	    for (k = 0; k < dof->type->dim; k++)
		values[i * dof->type->dim + k] += *data * bas[k];
    }

    return values;
}

static void
axpby_fn(DOF *y, ELEMENT *e, int bno, const FLOAT *ref, FLOAT *values)
{
    void **ctx = y->ctx;
    int i, n = DofDim(y);
    DOF *x = ctx[0];
    FLOAT a = *(FLOAT *)ctx[1], b = *(FLOAT *)ctx[2];

    if (n == 0)
	return;

    Unused(bno);
    if (a != 0.0) {
	phgDofEval(x, e, ref, values);
	if (a == -1.0)
	    for (i = 0; i < n; i++)
		values[i] = -values[i];
	else if (a != 1.0)
	    for (i = 0; i < n; i++)
		values[i] *= a;
    }
    else {
	bzero(values, n * sizeof(*values));
    }

    if (b != 0.0) {
	FLOAT tmp[n];
	phgDofEval(y, e, ref, tmp);
	if (b == -1.0)
	    for (i = 0; i < n; i++)
		values[i] -= tmp[i];
	else if (b == 1.0)
	    for (i = 0; i < n; i++)
		values[i] += tmp[i];
	else
	    for (i = 0; i < n; i++)
		values[i] += b * tmp[i];
    }
}

DOF *
phgDofAXPBY_(FLOAT a, DOF *x, FLOAT b, DOF **py,
	     const char *srcfile, int srcline, BOOLEAN check)
{
    size_t i;
    DOF *y = py == NULL ? NULL : *py;
    int nbas = DofNBas(x, NULL);

    if (y == NULL) {
	y = phgDofNew_(x->g, x->type, NULL, x->dim, "y", DofNoAction,
							srcfile, srcline);
	b = 0.0;
    }
    if (py != NULL)
	*py = y;

    assert(DofDim(x) == DofDim(y) && x->g == y->g);

    if (!SpecialDofType(x->type)) {
	for (i = 0; i < x->g->nelem * (size_t)nbas * x->dim; i++)
	    y->data[i] = a * x->data[i] + b * y->data[i];
	phgDofUpdateGhost(y);
    }
    else {
	void *ctx[] = {x, &a, &b};
	y->ctx = ctx;
	phgDofSetDataByLambdaFunction(y, axpby_fn);
	y->ctx = NULL;
    }

    return y;
}

DOF *
phgDofCopy_(DOF *src, DOF **pdest, DOF_TYPE *newtype,
			const char *name, const char *srcfile, int srcline)
{
    DOF *dest = NULL;

    assert(src != NULL);

    if (pdest != NULL)
	dest = *pdest;
    if (newtype == NULL)
	newtype = dest == NULL ? src->type : dest->type;
    if (dest != NULL)
	phgDofFree(&dest);
    dest = phgDofNew_(src->g, newtype, NULL, DofDim(src) / newtype->dim,
		      name == NULL ? src->name : name, DofNoAction,
		      srcfile, srcline);
    phgDofAXPBY_(1.0, src, 0.0, &dest, srcfile, srcline, FALSE);
    if (pdest != NULL)
	*pdest = dest;

    return dest;
}

DOF *
dof_derivation(DOF *src, DOF **newdof, DOF_TYPE *newtype,
	 	const char *name, int dim_mul, int dim_div,
	 	DOF_USER_FUNC_LAMBDA cb_func,
		  const char *srcfile, int srcline)
/* This is the generic function for computing derivatives of a DOF object. */
{
    DOF *dest = NULL;
    int dim = DofDim(src) * dim_mul / dim_div;

    assert(DofDim(src) % dim_div == 0);
    assert(src->type != DOF_ANALYTIC);
    if (newtype == NULL)
	newtype = src->type->grad_type;
    assert(dim % newtype->dim == 0);
    assert(name != NULL);

    if (newdof != NULL)
	dest = *newdof;
    if (dest != NULL)
	phgDofFree(&dest);
    dest = phgDofNew_(src->g, newtype, NULL, dim / newtype->dim, name,
		      DofNoAction, srcfile, srcline);
    dest->ctx = src;	/* for use in cb_func */
    phgDofSetDataByLambdaFunction(dest, cb_func);
    if (newdof != NULL)
	*newdof = dest;

    return dest;
}

static void
cb_grad(DOF *u, ELEMENT *e, int bno, const FLOAT xref[], FLOAT *values)
/* callback function for phgDofGradient_.
 * This function evaluates the gradient of "u->ctx" at the given point.
 */
{
    FLOAT dx[Dim], *data, (*corners)[Dim];
    DOF *src = u->ctx;
    int i, j, k, n, dof_dim = src->dim, type_dim = src->type->dim;

    n = DofNBas(src, NULL);
    data = src->data + e->index * n * dof_dim;
    corners = phgGeomGetCorners(src->g, e, NULL);
    dx[0] = 1.0 / (corners[1][0] - corners[0][0]);
    dx[1] = 1.0 / (corners[1][1] - corners[0][1]);
    dx[2] = 1.0 / (corners[1][2] - corners[0][2]);
    /* data[n][dof_dim] -> values[dof_dim][type_dim][Dim] */
    for (j = 0; j < dof_dim * type_dim * Dim; j++)
	values[j] = 0.;
    for (i = 0; i < n; i++, data += dof_dim) {
	FLOAT grad0[type_dim][Dim];
	gradBasisFunc(src, xref[0], xref[1], xref[2], i, (FLOAT *)grad0);
	for (j = 0; j < dof_dim; j++) {
	    for (k = 0; k < type_dim; k++) {
		values[(j * type_dim + k) * Dim + 0]
				+= grad0[k][0] * dx[0] * data[j];
		values[(j * type_dim + k) * Dim + 1]
				+= grad0[k][1] * dx[1] * data[j];
		values[(j * type_dim + k) * Dim + 2]
				+= grad0[k][2] * dx[2] * data[j];
	    }
	}
    }
}

DOF *
phgDofGradient_(DOF *src, DOF **newdof, DOF_TYPE *newtype,
		const char *name, const char *srcfile, int srcline)
{
    return dof_derivation(src, newdof, newtype, name == NULL ? "grad" : name,
			  Dim, 1, cb_grad, srcfile, srcline);
}

static void
cb_curl(DOF *u, ELEMENT *e, int bno, const FLOAT xref[], FLOAT *values)
/* callback function for phgDofCurl_.
 * This function evaluates the curl of "u->ctx" at the given point.
 */
{
    DOF *src = u->ctx;
    FLOAT grad[DofDim(src)][Dim], (*g)[Dim] = grad;
    int i, n = DofDim(src) / Dim;

    /* First, compute grad(src) stored as grad[DofDof(src)][Dim] */
    cb_grad(u, e, bno, xref, (FLOAT *)g);

    /* Then, convert gradient to curl */
    for (i = 0; i < n; i++, values += Dim, g += Dim) {
	values[0] = g[2][1] - g[1][2];
	values[1] = g[0][2] - g[2][0];
	values[2] = g[1][0] - g[0][1];
    }
}

DOF *
phgDofCurl_(DOF *src, DOF **newdof, DOF_TYPE *newtype,
		const char *name, const char *srcfile, int srcline)
{
    return dof_derivation(src, newdof, newtype, name == NULL ? "curl" : name,
			  Dim, Dim, cb_curl, srcfile, srcline);
}

static void
cb_div(DOF *u, ELEMENT *e, int bno, const FLOAT xref[], FLOAT *values)
/* callback function for phgDofCurl_.
 * This function evaluates the curl of "u->ctx" at the given point.
 */
{
    DOF *src = u->ctx;
    FLOAT grad[DofDim(src)][Dim], (*g)[Dim] = grad;
    int i, n = DofDim(src) / Dim;

    /* First, compute grad(src) stored as grad[DofDof(src)][Dim] */
    cb_grad(u, e, bno, xref, (FLOAT *)g);

    /* Then, convert gradient to curl */
    for (i = 0; i < n; i++, values += 1, g += Dim)
	*values = g[0][0] + g[1][1] + g[2][2];
}

DOF *
phgDofDivergence_(DOF *src, DOF **newdof, DOF_TYPE *newtype,
		  const char *name, const char *srcfile, int srcline)
{
    return dof_derivation(src, newdof, newtype, name == NULL ? "div" : name,
			  1, Dim, cb_div, srcfile, srcline);
}

FLOAT
phgDofNormL2_(DOF *u, int quad_order)
{
    ELEMENT *e;
    FLOAT res;
    QCACHE *qc;
    int Q_u;

    qc = phgQCNew(QD_P4EST, u);
    Q_u = phgQCAddFEFunction(qc, u);
    res = 0.;
#if USE_OMP
# pragma omp parallel for private(e) reduction(+:res)
#endif	/* USE_OMP */
    ForAllElementsBegin(u->g, e)
	FLOAT *rule;
	int o = 2 * DofTypeOrder(u, e);
	if (o < 0)
	    o = 2;
	if (quad_order >= 0)
	    o = quad_order;
	phgQuadInterfaceCuboid(NULL, 0, NULL, e->corners, o, &rule, NULL, NULL);
	phgQCSetRule(qc, rule, -1.0);
	res += phgQCIntegrate(qc, e->index, Q_u, 0, qc, e->index, Q_u, 0);
	phgFree(rule);
    ForAllElementsEnd
    phgQCFree(&qc);

#if USE_MPI
    if (u->g->nprocs > 1) {
	FLOAT res0 = res;
	MPI_Allreduce(&res0, &res, 1, PHG_MPI_FLOAT, PHG_SUM, u->g->comm);
    }
#endif	/* USE_MPI */

    return Sqrt(res);
}

/*------------------------- BEGIN list of neighbours -------------------------*/

static INT (*entries)[Dim + 1];	/* integer coordinates, e->index */

static int
comp_entry(const void *p0, const void *p1)
{
    int i;
    INT *c0 = entries[*(const int *)p0], *c1 = entries[*(const int *)p1];

    for (i = Dim - 1; i >= 0; i--) {
	INT a = (c0[i] & ~1), b = (c1[i] & ~1);
	if (a < b - 1)
	    return -1;
	else if (a > b + 1)
	    return 1;
    }

    return 0;
}

static INT
*get_vertex_coord(const FLOAT corners[][Dim], int v, const FLOAT trans[][Dim])
{
    static INT out[Dim];
#if USE_OMP
# pragma omp threadprivate(out)
#endif	/* USE_OMP */
    int bits[] = {v & 1, (v >> 1) & 1, (v >> 2) & 1};

    out[0] = (INT)((corners[bits[0]][0] - trans[0][0]) * trans[1][0] + 0.5);
    out[1] = (INT)((corners[bits[1]][1] - trans[0][1]) * trans[1][1] + 0.5);
    out[2] = (INT)((corners[bits[2]][2] - trans[0][2]) * trans[1][2] + 0.5);

    return out;
}

int
pack_list(int n, INT *list, int nex, INT *exlist)
/* performs "sort -u list". exlist[nex] is the list of excluded indices */ 
{
    int i, j;

    if (n <= 0)
	return n;

    if (nex > 0)
	qsort(exlist, nex, sizeof(*exlist), phgCompint);
    qsort(list, n, sizeof(*list), phgCompint);

    for (j = i = 0; i < n; ) {
	list[j] = list[i++];
	while (i < n && list[i] == list[j])
	    i++;
	if (nex <= 0 ||
	    bsearch(&list[j], exlist, nex, sizeof(*exlist), phgCompint) == NULL)
	    j++;
    }

    return j;
}

NLIST *
phgNeighboursListNew(GRID *g)
/* This function constructs list of edge and vertex neighbours.
 *
 * Note:
 *   1.	The list are actually valid only for elements with e->mark == 0.
 *   2.	A slow algorithm which collects all (vert,elem) pairs with (O(nlog(n))
 *	complexity) is used, instead of tracing face neighbours (O(n)),
 *	because phgGetNeighbour() does not work with ghost elements.
 *   3. The algorithm requires 2:1 element size ratio, since in this case
 *	all neighbours, whether they are smaller, larger or of the same size,
 *	share at least one vertex thus will be found in the (vert, elem) pairs.
 *   4. Only one of the smaller neighbours on an edge or a face will be
 *   	collected, which fullfills our needs.
 *   5. For each element, the face, edge and vertex neighbours are ordered
 *	according to elementwise face, edge and edge numbers, with -1 for
 *	unavailable or non-existent neighbour (outside of mesh or submesh).
 */
{
    NLIST *nl = phgCalloc(1, sizeof(*nl));
    ELEMENT *e;
    FLOAT bbox[2][Dim], dd;
    double t0 = phgGetTime(NULL);
    INT i, n;
    size_t total;
    int ii, kk, nn, *order;

    if (g->nelem == 0)
	return nl;

    phgSetupHalo(g, HALO_VERT);

    /* Note: minus 2 more bits, the lowest bit encodes the vertex no. */
    dd = (FLOAT)(1L << (sizeof(INT) * 8 - 3));
    bbox[1][0] = dd / (g->bbox[1][0] - (bbox[0][0] = g->bbox[0][0]));
    bbox[1][1] = dd / (g->bbox[1][1] - (bbox[0][1] = g->bbox[0][1]));
    bbox[1][2] = dd / (g->bbox[1][2] - (bbox[0][2] = g->bbox[0][2]));

    assert(NEdge >= NVert);		/* can use the entries[] buffer */
    total = g->nelem + g->nghost;
    assert(total * (size_t)NVert <= INT_MAX);
    entries = phgAlloc(total * NVert * sizeof(*entries));
    order = phgAlloc(total * (size_t)NVert * sizeof(*order));

    /* entries := all (nvert,elem) pairs, NVert entries per element, with:
     *		entries[][0:Dim-1] := scaled coordinates of the vertex,
     *		entries[][Dim] := element index */
    n = 0;
    ForAllElements_(g, e, FOR_ALL) {
	/* Note: encode vertex number into scaled coordinates */
	for (ii = 0; ii < NVert; ii++, n++) {
	    INT *c = get_vertex_coord(e->corners, ii, bbox);
	    entries[n][0] = (c[0] << 1) | ((ii >> 0) & 1);
	    entries[n][1] = (c[1] << 1) | ((ii >> 1) & 1);
	    entries[n][2] = (c[2] << 1) | ((ii >> 2) & 1);
	    entries[n][Dim] = e->index;
	}
    }
    /* sort entries */
    for (i = 0; i < n; i++)
	order[i] = i;
    qsort(order, n, sizeof(*order), comp_entry);
#if 0
#warning
FILE *f = fopen("tmp.lst", "w");
for (i = 0; i < n; i++)
fprintf(f, "entries[%d] = {%d, %d, %d, %d}\n", i, entries[order[i]][0], entries[order[i]][1], entries[order[i]][2], entries[order[i]][3]);
fclose(f);
#endif

    nl->vloc = phgAlloc((g->nelem + g->nghost + 1) * sizeof(*nl->vloc));
    nl->eloc = phgAlloc((g->nelem + g->nghost + 1) * sizeof(*nl->eloc));
    nl->floc = phgAlloc((g->nelem + g->nghost + 1) * sizeof(*nl->floc));

    /* total := # of elements with mark == 0 */
    total = 0;
    i = 0;
    ForAllElements_(g, e, FOR_ALL) {
	assert(i == e->index);
	nl->vloc[i] = NVert * total;
	nl->eloc[i] = NEdge * total;
	nl->floc[i] = NFace * total;
	if (e->mark == 0)
	    total++;
	i++;
    }
    nl->vloc[i] = NVert * total;
    nl->eloc[i] = NEdge * total;
    nl->floc[i] = NFace * total;
    nl->vlist = phgCalloc(NVert * (size_t)total, sizeof(*nl->vlist));
    nl->elist = phgCalloc(NEdge * (size_t)total, sizeof(*nl->elist));
    nl->flist = phgCalloc(NFace * (size_t)total, sizeof(*nl->flist));
    for (i = 0; i < NVert * total; i++)
	nl->vlist[i] = -1;
    for (i = 0; i < NEdge * total; i++)
	nl->elist[i] = -1;
    for (i = 0; i < NFace * total; i++)
	nl->flist[i] = -1;

    for (i = 0; i < n; ) {
	INT end, tmp[8];
	/* [i, end) := the list of all elements on the vertex */
	for (end = i + 1;
	     end < n && !comp_entry(order + i, order + end); end++);
	assert(end - i <= 8);
	if (end - i <= 1) {
	    i = end;
	    continue;
	}
	/* arrange elements as a 2x2x2 array in tmp[] */
	for (ii = 0; ii < 8; ii++)
	    tmp[ii] = -1;
	for (; i < end; i++) {
	    ii = ((entries[order[i]][0] & 1) << 0) |
		 ((entries[order[i]][1] & 1) << 1) |
		 ((entries[order[i]][2] & 1) << 2);
	    tmp[7 - ii] = entries[order[i]][Dim];
	}
	for (ii = 0; ii < 8; ii++) {
	    if (tmp[ii] < 0 || g->elems[tmp[ii]]->mark != 0)
		continue;
	    /* face neighbours */
	    for (kk = 0; kk < 3; kk++) {
		int b = 1 & ~(ii >> kk);
		nn = (ii & ~(1 << kk)) | (b << kk);
		if (tmp[nn] < 0)
		    continue;
		nl->flist[nl->floc[tmp[ii]] + 2 * kk + b] = tmp[nn];
	    }
	    /* edge neighbours */
	    for (kk = 0; kk < 3; kk++) {
		int f0, f1, edge;
		if ((f0 = 0) == kk)
		    f0++;
		if ((f1 = f0 + 1) == kk)
		    f1++;
		nn = (ii & ~(1 << f0)) | ((1 & ~(ii >> f0)) << f0);
		nn = (nn & ~(1 << f1)) | ((1 & ~(nn >> f1)) << f1);
		if (tmp[nn] < 0)
		    continue;
		f0 = 2 * f0 + (1 & ~(ii >> f0));
		f1 = 2 * f1 + (1 & ~(ii >> f1));
		for (edge = 0; edge < NEdge; edge++)	/* FIXME: use a table */
		    if (p8est_edge_faces[edge][0] == f0 &&
			p8est_edge_faces[edge][1] == f1)
			break;
		assert(edge < NEdge);
		nl->elist[nl->eloc[tmp[ii]] + edge] = tmp[nn];
	    }
	    /* vertex neighbour */
	    nn = 7 - ii;
	    if (tmp[nn] >= 0)
		nl->vlist[nl->vloc[tmp[ii]] + 7 - ii] = tmp[nn];
	}
    }	/* i */

    phgFree(entries);
    phgFree(order);

    if (phgVerbosity > 0)
	phgPrintf("*** phgNeighboursListNew: %0.2gs.\n", phgGetTime(NULL) - t0);

    return nl;
}

void
phgNeighboursListFree(NLIST **nl)
{
    if (nl == NULL || *nl == NULL)
	return;

    phgFree((*nl)->flist);
    phgFree((*nl)->floc);
    phgFree((*nl)->elist);
    phgFree((*nl)->eloc);
    phgFree((*nl)->vlist);
    phgFree((*nl)->vloc);
    phgFree(*nl);

    *nl = NULL;
}

/*-------------------------- END list of neighbours --------------------------*/

/* The "ref_info" struct saves needed data of the parent during refinement */
static struct {
    GRID	*g;		/* The GRID */
    GRID	*g0, *g1;	/* The dummy GRIDs */
    DOF		*dof0;		/* pass dof0 to interp_fn() */
    INT		data_count;	/* data count per element of all DOFs */
    INT		count;		/* for debugging */
} ref_info;
#if USE_OMP
# pragma omp threadprivate(ref_info)
#endif	/* USE_OMP */

static void
interp_fn(DOF *u, ELEMENT *e, int bno, const FLOAT *ref, FLOAT *values)
{
    ELEMENT *e0 = ref_info.g0->elems[0];
    FLOAT xyz[Dim], ref0[Dim];

    Unused(bno);
    ref_to_xyz(u->g, e->index, ref, xyz);
    xyz_to_ref(ref_info.g0, e0->index, xyz, ref0);
    phgDofEval(ref_info.dof0, e0, ref0, values);
}

static int
refine_fn(p8est_t *p8est, p8est_topidx_t which_tree, p8est_quadrant_t *quad)
/* refine_fn() for p8est_refine() */
{
    ELEMENT *e = quad->p.user_data;

    return e->mark <= 0 ? 0 : 1;
}

static void
replace_fn(p8est_t * p8est, p4est_topidx_t which_tree,
		int num_outgoing, p8est_quadrant_t *outgoing[],
		int num_incoming, p8est_quadrant_t *incoming[])
/* Note: outgoing and incoming are respectively parents and children, during
 * refine normally num_outgoing == 1 and num_incoming == 8 */
{
    GRID *g = ref_info.g;
    int i, j;

    ref_info.count += num_outgoing;
    assert(num_outgoing * 8 == num_incoming);

    for (; num_outgoing > 0; num_outgoing--, outgoing++) {
	ELEMENT *ep = (*outgoing)->p.user_data, *dummy_e0;

	/* parent ELEMENT */
	ep->quad = *outgoing;	/* Note: the quad address may change */
	dummy_e0 = ref_info.g0->elems[0];
	memcpy(dummy_e0->corners, ep->corners, sizeof(ep->corners));

	/* transfer data: parent -> children */

	for (j = 0; j < 8; j++, incoming++) {
	    ELEMENT *e = (*incoming)->p.user_data, *dummy_e1;
	    DOF *dof, **pdof, **pd0, **pd1;
	    FLOAT *p0, *p1;

	    bzero(e, sizeof(*e));
	    e->quad = *incoming;
	    e->generation = ep->generation + 1;
	    e->treeid = ep->treeid;
	    e->tree = ep->tree;
	    get_corners(g, e, e->corners);

#if DEBUG
	    p8est_quadrant_t tmp;
	    assert(p8est_quadrant_is_inside_tree(e->tree, e->quad));
	    p8est_quadrant_parent(e->quad, &tmp);
	    assert(tmp.x == ep->quad->x && tmp.y == ep->quad->y &&
		   tmp.z == ep->quad->z);
#endif	/* DEBUG */

	    /* set up region_mark and bound_type[] using data of the parent */
	    e->mark = ep->mark - 1;
	    e->region_mark = ep->region_mark;
	    for (i = 0; i < 6; i += 2) {
		p4est_qcoord_t c = 0;
		switch (i / 2) {
		    case 0: c = e->quad->x - ep->quad->x; break;
		    case 1: c = e->quad->y - ep->quad->y; break;
		    case 2: c = e->quad->z - ep->quad->z; break;
		}
		if (c == 0) {
		    e->bound_type[i + 0] = ep->bound_type[i + 0];
		    e->bound_type[i + 1] = INTERIOR;	/* interior face */
		}
		else {
		    e->bound_type[i + 0] = INTERIOR;	/* interior face */
		    e->bound_type[i + 1] = ep->bound_type[i + 1];
		}
	    }

	    if (ref_info.data_count == 0)
		continue;

	    /* set up dummy_e1->corners for interpolation */
	    dummy_e1 = ref_info.g1->elems[0];
	    memcpy(dummy_e1->corners, e->corners, sizeof(e->corners));

	    /* interpolate DOFs */
	    pd0 = ref_info.g0->dofs;
	    pd1 = ref_info.g1->dofs;
	    e->reserved = phgAlloc(ref_info.data_count * sizeof(FLOAT));
	    p0 = ep->reserved;
	    p1 = e->reserved;
	    for (pdof = g->dofs; pdof != NULL && *pdof != NULL; pdof++) {
		dof = *pdof;
		if (SpecialDofType(dof->type) ||
		    dof->userfunc != DofInterpolation)
		    continue;
		(*pd0)->data = p0;
		ref_info.dof0 = *pd0;	/* for interp_fn() */
		(*pd1)->data = p1;
		phgDofSetDataByLambdaFunction(*pd1, interp_fn);
		i = DofNBas(dof, NULL) * dof->dim;
		p0 += i;
		p1 += i;
		pd0++;
		pd1++;
	    } /* pdof */
	} /* incoming */

	phgFree(ep->reserved);
	ep->reserved = NULL;
    } /* outgoing */
}

void
phgRefineMarkedElements(GRID *g)
{
    DOF *dof, **pdof;
    ELEMENT *e;
    FLOAT *p;
    double t0 = phgGetTime(NULL), tt;
    INT nrefined;
    int n;

    remove_ghost(g);

    g->serial_no++;

    /* Attach data of DOFs of with the "DofInterpolation" tag to elements,
     * and prepare for performing interpolation */

    /* prepare the dummy grid */
    ref_info.g0 = phgNewGrid(-1);
    phgImport_(ref_info.g0, NULL, MPI_COMM_SELF);	/* the unit cube */
    ref_info.g1 = phgNewGrid(-1);
    phgImport_(ref_info.g1, NULL, MPI_COMM_SELF);	/* the unit cube */

    /* parse list of DOFs */
    ref_info.data_count = 0;
    for (pdof = g->dofs; pdof != NULL && (dof = *pdof) != NULL; pdof++) {
	if (SpecialDofType(dof->type) || dof->userfunc != DofInterpolation)
	    continue;
	ref_info.data_count += DofNBas(dof, NULL) * dof->dim;
	/* create dummy DOFs for interpolation for this DOF */
	phgDofNew(ref_info.g0, dof->type, dof->dim, "src", DofNoData);
	phgDofNew(ref_info.g1, dof->type, dof->dim, "dst", DofNoData);
    }

    if (ref_info.data_count > 0) {
	/* attach DOF data to elements */
	ForAllElements(g, e) {
	    p = e->reserved = phgAlloc(ref_info.data_count * sizeof(FLOAT));
	    for (pdof = g->dofs; pdof != NULL && *pdof != NULL; pdof++) {
		dof = *pdof;
		if (SpecialDofType(dof->type) ||
		    dof->userfunc != DofInterpolation)
		    continue;
		n = DofNBas(dof, NULL) * dof->dim;
		memcpy(p, dof->data + e->index * n, n * sizeof(FLOAT));
		p += n;
	    }
	    assert(p - (FLOAT *)e->reserved == ref_info.data_count);
	}
    }

    ref_info.g = g;

#if 0
#warning
INT count = 0;
ForAllElements(g, e) count += e->mark;
phgInfo(-1, "%s: %"dFMT" elements marked\n", __func__, count);
#endif
    ref_info.count = 0;
    tt = phgGetTime(NULL);
    p8est_refine_ext(g->p8est, /*recursive=*/1, /*maxlevel*/-1,
		     refine_fn, /*init_fn*/NULL, replace_fn);
    nrefined = ref_info.count;
    if (phgVerbosity > 1)
	phgPrintf("*** %s: refine grid, %0.2gs\n", __func__,
						phgGetTime(NULL) - tt);

    /* Call p8est_balance to ensures 2:1 size difference with more refines */
    /* FIXME: if using P8EST_CONNECT_FACE, then for the following case:
     * 		touch p4est/phg-to-p4est.c
     * 		make USER_CFLAGS="-DFORCE_BALANCE=1" -s lib
     * 		cd examples
     * 		make USER_CFLAGS="-DPHG_TO_P4EST" -s clean ipdg
     * 		mpirun -np 6 ./ipdg -refine 2	(or any np>6)
     * the following assertion at p4est_mesh.c:502 (p4est-2.2) fails:
     *		P4EST_ASSERT (side2->is.hanging.quad[h] != NULL);
     *
     * FIXME: p8est_balance_ext() much slower (50x) than p8est_refine_ext()!
     */
    ref_info.count = 0;
    tt = phgGetTime(NULL);
    p8est_balance_ext(g->p8est, P8EST_CONNECT_CORNER, /*init_fn*/NULL,
		      replace_fn);
    if (phgVerbosity > 1)
	phgPrintf("*** %s: balance grid, %0.2gs\n", __func__,
						phgGetTime(NULL) - tt);

    /* Free dummy DOFs and GRID */
    while ((pdof = ref_info.g0->dofs) != NULL && (dof = *pdof) != NULL) {
	dof->data = NULL;
	phgDofFree(&dof);
    }
    phgFreeGrid(&ref_info.g0);

    while ((pdof = ref_info.g1->dofs) != NULL && (dof = *pdof) != NULL) {
	dof->data = NULL;
	phgDofFree(&dof);
    }
    phgFreeGrid(&ref_info.g1);

    /* update ELEMENT data */
    g->balanced = FALSE;
    g->serial_no++;
    update_elems(g);

    tt = phgGetTime(NULL);
    /* restore DOF data */
    for (pdof = g->dofs; pdof != NULL && *pdof != NULL; pdof++) {
	dof = *pdof;
	if (SpecialDofType(dof->type) || dof->userfunc == DofNoData)
	    continue;
	phgFree(dof->data);
	dof->data = phgCalloc(g->nelem,
			      DofNBas(dof, NULL) * dof->dim * sizeof(FLOAT));
	if (!SpecialDofUserFunction(dof->userfunc))
	    phgDofSetDataByFunction(dof, dof->userfunc);
	else if (dof->userfunc == DofLambdaFunction)
	    phgDofSetDataByLambdaFunction(dof, dof->userfunc_lambda);
    }

    if (ref_info.data_count > 0) {
	/* copy DOF data from elements */
	ForAllElements(g, e) {
	    assert(e->reserved != NULL);
	    p = e->reserved;
	    for (pdof = g->dofs; pdof != NULL && *pdof != NULL; pdof++) {
		dof = *pdof;
		if (SpecialDofType(dof->type) ||
		    dof->userfunc != DofInterpolation)
		    continue;
		n = DofNBas(dof, NULL) * dof->dim;
		memcpy(dof->data + e->index * n, p, n * sizeof(FLOAT));
		p += n;
	    }
	    assert(p - (FLOAT *)e->reserved == ref_info.data_count);
	    phgFree(e->reserved);
	    e->reserved = NULL;
	}
    }
    if (phgVerbosity > 1)
	phgPrintf("*** %s: update DOFs, %0.2gs\n", __func__,
						phgGetTime(NULL) - tt);

    if (phgVerbosity > 0) {
#if USE_MPI
	if (g->nprocs > 1) {
	    INT n0[2] = {nrefined, ref_info.count}, n1[2]; 
	    MPI_Reduce(n0, n1, 2, PHG_MPI_INT, MPI_SUM, 0, g->comm);
	    nrefined = n1[0];
	    ref_info.count = n1[1];
	}
#endif	/* USE_MPI */
	phgPrintf("*** phgRefineMarkedElements: %"dFMT"/%"dFMT" refined/"
		  "balanced, %"dFMT" total, %0.2gs.\n",
		  nrefined, ref_info.count, g->nelem_global,
		  phgGetTime(NULL) - t0);
    }
}

void
phgRefineAllElements(GRID *g, int level)
{
    ELEMENT *e;
    int d = 0;

    ForAllElements(g, e) {
#ifndef FORCE_BALANCE
# define FORCE_BALANCE 0
#endif
#if FORCE_BALANCE
#warning debugging code!
#if 0
	d = (g->nelem_global > 1 && e->index == 0) ? 2 : 0;
#else
	d = ((e->corners[0][0] + e->corners[0][1] + e->corners[0][2] == 0. ||
	      e->corners[1][0] + e->corners[1][1] + e->corners[1][2] == 3.) &&
	     g->nelem_global > 1) ? 2 : 0;
#endif
#endif	/* FORCE_BALANCE */
	e->mark = level + d;
    }
    phgRefineMarkedElements(g);
}

int
phgBalanceGrid(GRID *g, FLOAT lif_threshold, INT submesh_threshold,
		DOF *weights, FLOAT power)
{
    p8est_t		*back = NULL;
    p8est_topidx_t	tid;
    size_t		zz;
    p8est_tree_t	*tree;
    ELEMENT *e;
    DOF *dof, **pdof;
    FLOAT *dst_data = NULL, *src_data = NULL, *p;
    int count = 0, n;

    if (g->balanced || g->lif <= lif_threshold)
	return 1;

    remove_ghost(g);

    g->serial_no += 1000;

    /* count DOF data size */
    for (pdof = g->dofs; pdof != NULL && (dof = *pdof) != NULL; pdof++) {
	if (SpecialDofType(dof->type) || dof->userfunc == DofNoData)
	    continue;
	count += dof->dim * DofNBas(dof, NULL);
    }

#define TEST_TRANSF 0

    if (count > 0) {
#if TEST_TRANSF
	count += 6;	
#endif	/* TEST_TRANSF */
	/* copy DOFs data to src_data */
	back = p8est_copy(g->p8est, /*copy_data*/0);
	src_data = phgAlloc(count * sizeof(FLOAT) * g->nelem);
	p = src_data;
	for (tid = g->p8est->first_local_tree;
	     tid <= g->p8est->last_local_tree; ++tid) {
	    tree = p8est_tree_array_index(g->p8est->trees, tid);
	    for (zz = 0; zz < tree->quadrants.elem_count; ++zz) {
		e = p8est_quadrant_array_index(&tree->quadrants, zz)
								->p.user_data;
#if TEST_TRANSF
		memcpy(p, e->corners, sizeof(e->corners));
		p += 6;
#endif	/* TEST_TRANSF */
		for (pdof = g->dofs; pdof != NULL && *pdof != NULL; pdof++) {
		    dof = *pdof;
		    if (SpecialDofType(dof->type) || dof->userfunc == DofNoData)
			continue;
		    n = dof->dim * DofNBas(dof, NULL);
		    memcpy(p, dof->data + e->index * n, n * sizeof(FLOAT));
		    p += n;
		}
	    }
	}
	assert(p - src_data == count * g->nelem);
    }
#if TEST_TRANSF
    count -= 6;	
#endif	/* TEST_TRANSF */

    p8est_partition(g->p8est, /*allow_coarsening*/0, /*weight_fn*/NULL);
    g->serial_no += 65536;
    update_elems(g);

    if (count == 0)
	return 1;

#if TEST_TRANSF
    count += 6;	
#endif	/* TEST_TRANSF */
    dst_data = phgAlloc(count * sizeof(FLOAT) * g->nelem);
    p8est_transfer_fixed(g->p8est->global_first_quadrant,
			 back->global_first_quadrant,
			 g->p8est->mpicomm, /*tag*/0,
			 dst_data, src_data, count * sizeof(FLOAT));

    phgFree(src_data);
    p8est_destroy(back);

    /* copy dst_data to DOFs data */
    for (pdof = g->dofs; pdof != NULL && (dof = *pdof) != NULL; pdof++) {
	if (SpecialDofType(dof->type) || dof->userfunc == DofNoData)
	    continue;
	phgFree(dof->data);
	dof->data = phgAlloc(DofNBas(dof, NULL) * dof->dim
					* g->nelem * sizeof(FLOAT));
    }

    p = dst_data;
    for (tid = g->p8est->first_local_tree;
	 tid <= g->p8est->last_local_tree; ++tid) {
	tree = p8est_tree_array_index(g->p8est->trees, tid);
	for (zz = 0; zz < tree->quadrants.elem_count; ++zz) {
	    e = p8est_quadrant_array_index(&tree->quadrants, zz)->p.user_data;
#if TEST_TRANSF
	    assert(!memcmp(p, e->corners, sizeof(e->corners)));
	    p += 6;
#endif	/* TEST_TRANSF */
	    for (pdof = g->dofs; pdof != NULL && *pdof != NULL; pdof++) {
		dof = *pdof;
		if (SpecialDofType(dof->type) || dof->userfunc == DofNoData)
		    continue;
		n = dof->dim * DofNBas(dof, NULL);
		assert(e->index < g->nelem);
		memcpy(dof->data + e->index * n, p, n * sizeof(FLOAT));
		p += n;
	    }
	}
    }

    phgFree(dst_data);

    g->balanced = TRUE;

    return 1;
}

/*---------------------------------------------------------------------------*/

struct {
    DOF *ls, *grad;
    ELEMENT *e;
} dummy_ctx;
#if USE_OMP
#pragma omp threadprivate(dummy_ctx)
#endif	/* USE_OMP */

static void
dummy_ls(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT xref[Dim + 1];

    if (!SpecialDofUserFunction(dummy_ctx.ls->userfunc)) {
	dummy_ctx.ls->userfunc(x, y, z, value);
	return;
    }
    phgGeomXYZ2Lambda(dummy_ctx.ls->g, dummy_ctx.e, x, y, z, xref);
    phgDofEval(dummy_ctx.ls, dummy_ctx.e, xref, value);
}

static void
dummy_grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    FLOAT xref[Dim + 1];

    if (!SpecialDofUserFunction(dummy_ctx.grad->userfunc)) {
	dummy_ctx.grad->userfunc(x, y, z, values);
	return;
    }
    phgGeomXYZ2Lambda(dummy_ctx.grad->g, dummy_ctx.e, x, y, z, xref);
    phgDofEval(dummy_ctx.grad, dummy_ctx.e, xref, values);
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
    DOF *grad;
    ELEMENT *e;
    INT n0 = 0;

    grad = ls_grad != NULL ? ls_grad : phgDofGradient(ls, NULL, NULL, NULL);

#if USE_OMP
#pragma omp parallel for private(e) reduction(+:n0)
#endif	/* USE_OMP */
    ForAllElementsBegin_(g, e, FOR_ALL)
	FLOAT (*corners)[Dim] = phgGeomGetCorners(g, e, NULL);
	dummy_ctx.ls = ls;
	dummy_ctx.grad = grad;
	dummy_ctx.e = e;
	e->mark = phgQuadInterfaceMarkElement(dummy_ls, ls->poly_order,
					dummy_grad, e->elem_type, corners);
	if (e->mark == 0)
	    n0++;
    ForAllElementsEnd
 
    if (ls_grad == NULL)
	phgDofFree(&grad);

    return n0;
}

/*---------------------------------------------------------------------------*/

FLOAT
phgGeomGetDiameter(GRID *g, ELEMENT *e)
{
    FLOAT (*box)[Dim] = phgGeomGetCorners(g, e, NULL);
    FLOAT dx = box[1][0] - box[0][0],
	  dy = box[1][1] - box[0][1],
	  dz = box[1][2] - box[0][2];

    return Sqrt(dx * dx + dy * dy + dz * dz);
}

FLOAT
phgGeomGetFaceDiameter(GRID *g, ELEMENT *e, int face)
{
    FLOAT (*box)[Dim] = phgGeomGetCorners(g, e, NULL);
    FLOAT dx, dy, dz;

    assert(face >= 0 && face < P8EST_FACES);
    face /= 2;
    dx = face == 0 ? 0. : box[1][0] - box[0][0];
    dy = face == 1 ? 0. : box[1][1] - box[0][1];
    dz = face == 2 ? 0. : box[1][2] - box[0][2];

    return Sqrt(dx * dx + dy * dy + dz * dz);
}

FLOAT
phgGeomGetVolume(GRID *g, ELEMENT *e)
{
    FLOAT (*box)[Dim] = phgGeomGetCorners(g, e, NULL);

    return (box[1][0] - box[0][0]) *
	   (box[1][1] - box[0][1]) *
	   (box[1][2] - box[0][2]);
}

INT
phgMapE2G(MAP *map, int dof_id, ELEMENT *e, int index)
{
    DOF *u;
    GRID *g;
    INT no, nelem;
    int rank;

    assert(dof_id >= 0 && dof_id < map->ndof);
    u = (DOF *)map->dofs[dof_id];
    g = u->g;
    if (ElementIsLocal(u->g, e)) {
	no = u->dim * e->index * DofNBas(u, e) + index;
	rank = g->rank;
	nelem = g->nelem;
    }
    else {
	/* find owner process of e by binary search */
	no = g->O2Gmap[e->index - g->nelem];
	rank = (int)phgBinarySearchINT(g->nprocs, g->partition, NULL, no);
	if (no != g->partition[rank])
	    rank--;
	/* skip processes with no element */
	while (rank < g->nprocs && g->partition[rank + 1] <= no)
	    rank++;
	assert(rank < g->nprocs);
	no = u->dim * (no - g->partition[rank]) * DofNBas(u, e) + index;
	nelem = g->partition[rank + 1] - g->partition[rank];
    }

    while (dof_id > 0) {
	u = (DOF *)map->dofs[--dof_id];
	no += u->dim * nelem * DofNBas(u, e);
    }

    return no + map->partition[rank];
}

#include <stdarg.h>

const char *
phgExportVTK(GRID *g, const char *fn, DOF *dof, ...)
{
    int k, a;
    BYTE face;
    size_t size;
    INT i;
    INT nleaf, nvert, nbface, nbface_global, e_no0;
    ELEMENT *e;
    size_t buffer_size;
    void *buffer;
    size_t float_size = !strcmp(phgOptionsGetKeyword("-vtk_precision"),
				"single") ? sizeof(float) : sizeof(double);
    /* vtk -> p4est vertex ordering in an element (the 1st 4 in a face) */
    int v_map[] = {0, 1, 3, 2, 4, 5, 7, 6};
    const char *vtk_precision = 
		(float_size == sizeof(double)) ? "double" : "float";
    /* variables for exporting DOFs */
    int dofs_n = 0, dofs_alloc = 0;
    DOF **dofs = NULL;
    char **name_list = NULL;
    int name_n = 0, name_alloc = 0;

#define VTK_PUT_FLOAT(a) {			\
    if (float_size == sizeof(float)) {		\
	float x = (float)(a);			\
	phgIOBigEndianAppend(&x, sizeof(x), 1);	\
    }						\
    else {					\
	double x = (double)(a);			\
	phgIOBigEndianAppend(&x, sizeof(x), 1);	\
    }						\
}

#if 0
    p8est_vtk_write_file(g->p8est, NULL, fn);
    return fn;
#endif

    nleaf = g->nelem;
    nvert = g->nelem * 8;

    fn = phgIOAddExtension(fn, "vtk");
    if (!phgIOOpen(g->comm, (char *)fn)) {
	phgError(0, "cannot create VTK file \"%s\".\n", fn);
	return NULL;
    }
    if (g->rank == 0)
	phgInfo(1, "creating VTK file \"%s\".\n", fn);

    nbface = nbface_global = 0;
    if (phgOptionsGetNoArg("-vtk_boundary_faces")) {
	/* count number of boundary faces */
	ForAllElements(g, e)
	    for (face = 0; face < 6; face++)
		if ((e->bound_type[face] & BDRY_MASK))
		    nbface++;
#if USE_MPI
	MPI_Allreduce(&nbface, &nbface_global, 1,PHG_MPI_INT, MPI_SUM, g->comm);
#else	/* USE_MPI */
	nbface_global = nbface;
#endif	/* USE_MPI */
    }

    buffer_size = nvert * float_size * 3;
    if ((size = nleaf * 9 * sizeof(a)) > buffer_size)
	buffer_size = size;
    if ((size = nbface * 5 * sizeof(a)) > buffer_size)
	buffer_size = size;

    buffer = phgIOAllocBuffer(buffer_size);

    /* File header and vertices */
    phgIOPrintf(g->comm, "# vtk DataFile Version 2.0\n"
		"Hexahedral grid created by PHG %s.\n"
		"BINARY\nDATASET UNSTRUCTURED_GRID\n"
		"POINTS %"dFMT" %s\n", PHG_VERSION,
		g->nelem_global * 8,	/* 8 vertices per element */
		vtk_precision);
    phgIOResetBuffer();
    ForAllElements(g, e) {
	FLOAT (*corners)[Dim] = phgGeomGetCorners(g, e, NULL);
	/* TODO: scale the element */
	for (k = 0; k < 8; k++) {
	    int bits[] = {v_map[k]&1, (v_map[k]>>1)&1, (v_map[k]>>2)&1};
	    for (a = 0; a < 3; a++)
		VTK_PUT_FLOAT(corners[bits[a]][a])
	}
    }
    phgIOWrite(g->comm, buffer, float_size * 3, nvert, g->nelem_global * 8,
		NULL, NULL, FALSE);

    /* element cells */
    phgIOPrintf(g->comm, "\nCELLS %"dFMT" %"dFMT"\n",
		g->nelem_global + nbface_global,
		g->nelem_global * 9 + nbface_global * 5);

#if USE_MPI
    MPI_Scan(&g->nleaf, &e_no0, 1, PHG_MPI_INT, MPI_SUM, g->comm);
    e_no0 -= g->nleaf;
#else	/* USE_MPI */
    e_no0 = 0;
#endif	/* USE_MPI */

    /* hexahedra */
    phgIOResetBuffer();
    i = e_no0 * 8;
    ForAllElements(g, e) {
	a = 8;
	phgIOBigEndianAppend(&a, sizeof(a), 1);
	/* vertices */
	for (k = 0; k < 8; k++) {
	    a = i + /*v_map[*/k/*]*/;
	    phgIOBigEndianAppend(&a, sizeof(a), 1);
	}
	i += 8;
    }
    phgIOWrite(g->comm, buffer, sizeof(a) * 9, nleaf,
		g->nelem_global, NULL, NULL, FALSE); 

    if (nbface_global > 0) {
	/* boundary faces */
	phgIOResetBuffer();
	i = e_no0 * 8;
	ForAllElements(g, e) {
	    for (face = 0; face < 6; face++) {
		if (!(e->bound_type[face] & BDRY_MASK))
		    continue;
		a = 4;
		phgIOBigEndianAppend(&a, sizeof(a), 1);
		/* vertices */
		for (k = 0; k < 4; k++) {
#define v_map_inv	v_map	/* v_map * v_map == I */
		    a = i + v_map_inv[p8est_face_corners[face][v_map[k]]];
		    phgIOBigEndianAppend(&a, sizeof(a), 1);
		}
	    }
	    i += 8;
	}
	phgIOWrite(g->comm, buffer, sizeof(a) * 5, nbface,
			nbface_global, NULL, NULL, FALSE); 
    }

    phgIOPrintf(g->comm, "\nCELL_TYPES %"dFMT"\n",
					g->nelem_global + nbface_global);

    /* types for hexahedra */
    phgIOResetBuffer();
    a = /*VTK_HEXAHEDRON=*/12;
    for (i = 0; i < nleaf; i++)
	phgIOBigEndianAppend(&a, sizeof(a), 1);
    phgIOWrite(g->comm, buffer, sizeof(a), nleaf, g->nelem_global,
							NULL, NULL, FALSE);

    /* types for quadrilaterals */
    if (nbface_global > 0) {
	phgIOResetBuffer();
	a = /*VTK_QUAD=*/9;
	for (i = 0; i < nbface; i++)
	    phgIOBigEndianAppend(&a, sizeof(a), 1);
	phgIOWrite(g->comm, buffer, sizeof(a), nbface, nbface_global,
						NULL, NULL, FALSE);
    }

    /* point data */
    if (dof != NULL) {
	va_list ap;
	int dim;
	char *name0, *name, *p;

	phgIOPrintf(g->comm, "\nPOINT_DATA %"dFMT"\n", g->nelem_global * 8);
	va_start(ap, dof);
	while (dof != NULL) {
	    if (dof->type == DOF_DGPn[0] || dof->type == DOF_DGQn[0]) {
		/* save cell dofs */
		if (dofs_n >= dofs_alloc)
		    dofs = phgRealloc_(dofs, (dofs_alloc += 8) * sizeof(*dofs),
					dofs_n * sizeof(*dofs));
		dofs[dofs_n++] = dof;
		dof = va_arg(ap, DOF *);
		continue;
	    }
	    dim = DofDim(dof) > 4 ? 4 : DofDim(dof); /* at most 4 components */
	    name = strdup(dof->name);
	    for (p = name; *p != '\0'; p++)
		if (isspace(*p))
		    *p = '_';
	    /* add "%d" suffix to avoid dulicate name */
	    a = 0;
	    name0 = strdup(name);
	    while (TRUE) {
		/* check for duplicate name */
		for (k = 0; k < name_n; k++)
		    if (!strcasecmp(name, name_list[k]))
			break;
		if (k >= name_n)
		    break;	/* no duplicate found */
		a++;
		p = malloc(strlen(name) + 20);
		free(name);
		name = p;
		sprintf(name, "%s_%d", name0, a);
	    }
	    free(name0);
	    /* add name to name_list */
	    if (name_n >= name_alloc)
		name_list = realloc(name_list,
				    sizeof(*name_list) * (name_alloc += 8));
	    name_list[name_n++] = name;
	    if (dim == Dim)
	        phgIOPrintf(g->comm, "\nVECTORS %s %s\n", name, vtk_precision);
	    else
	        phgIOPrintf(g->comm, "\nSCALARS %s %s %d\n"
			    "LOOKUP_TABLE default\n", name, vtk_precision, dim);
	    phgIOResetBuffer();
	    ForAllElements(g, e) {
		FLOAT values[DofDim(dof)];
		for (k = 0; k < 8; k++) {
		    FLOAT xref[] = {(FLOAT)(v_map[k]&1),
				    (FLOAT)((v_map[k]>>1)&1),
				    (FLOAT)((v_map[k]>>2)&1)};
		    phgDofEval(dof, e, xref, values);
		    for (a = 0; a < dim; a++)
			VTK_PUT_FLOAT(values[a])
		}
	    }
	    phgIOWrite(g->comm, buffer, float_size * dim, nvert,
			g->nelem_global * 8, NULL, NULL, FALSE);
	    dof = va_arg(ap, DOF *);
	}
	va_end(ap);
    }

    /* cell data */
    phgIOPrintf(g->comm, "\nCELL_DATA %"dFMT"\n",
					g->nelem_global + nbface_global);

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
	ForAllElements(g, e) {
	    for (face = 0; face < 6; face++) {
		if (!(e->bound_type[face] & BDRY_MASK))
		    continue;
		a = e->bound_type[face];
		phgIOBigEndianAppend(&a, sizeof(a), 1);
	    }
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
	ForAllElements(g, e) {
	    a = e->region_mark;
	    for (face = 0; face < 6; face++)
		if ((e->bound_type[face] & BDRY_MASK))
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

    for (i = 0; i < dofs_n; i++) {	/* save P0/Q0 DOFs as cell data */
	FLOAT xref[] = {0.5, 0.5, 0.5};
	int dim;
	char *name0, *name, *p;

	dof = dofs[i];

	dim = DofDim(dof) > 4 ? 4 : DofDim(dof); /* at most 4 components */
	name = strdup(dof->name);
	for (p = name; *p != '\0'; p++)
	    if (isspace(*p))
		*p = '_';
	/* add "%d" suffix to avoid dulicate name */
	a = 0;
	name0 = strdup(name);
	while (TRUE) {
	    /* check for duplicate name */
	    for (k = 0; k < name_n; k++)
		if (!strcasecmp(name, name_list[k]))
		    break;
	    if (k >= name_n)
		break;	/* no duplicate found */
	    a++;
	    p = malloc(strlen(name) + 20);
	    free(name);
	    name = p;
	    sprintf(name, "%s_%d", name0, a);
	}
	free(name0);
	/* add name to name_list */
	if (name_n >= name_alloc)
	    name_list = realloc(name_list,
				sizeof(*name_list) * (name_alloc += 8));
	name_list[name_n++] = name;
	if (dim == Dim)
	    phgIOPrintf(g->comm, "\nVECTORS %s %s\n", name, vtk_precision);
	else
	    phgIOPrintf(g->comm, "\nSCALARS %s %s %d\nLOOKUP_TABLE default\n",
				name, vtk_precision, dim);
	phgIOResetBuffer();
	ForAllElements(g, e) {
	    FLOAT values[DofDim(dof)];
	    phgDofEval(dof, e, xref, values);
	    for (a = 0; a < dim; a++)
		VTK_PUT_FLOAT(values[a])
	}
	phgIOWrite(g->comm, buffer, float_size * dim, nleaf,
				g->nelem_global, NULL, NULL, FALSE);
	phgIOResetBuffer();
	ForAllElements(g, e) {
	    for (face = 0; face < 6; face++) {
		if (!(e->bound_type[face] & BDRY_MASK))
		    continue;
		FLOAT values[DofDim(dof)];
		phgDofEval(dof, e, xref, values);
		for (a = 0; a < dim; a++)
		    VTK_PUT_FLOAT(values[a])
	    }
	}
	phgIOWrite(g->comm, buffer, sizeof(a), nbface, nbface_global,
						NULL, NULL, FALSE);
    }
    phgFree(dofs);

    /* free saved names */
    for (k = 0; k < name_n; k++)
	free(name_list[k]);
    free(name_list);

    phgIOClose();

    phgIOFreeBuffer();

    return fn;
}

MAP *
phgMapCreateN(int ndof, DOF **dofs)
{
    GRID *g = dofs[0]->g;
    MAP *map;
    INT n, N;
    int i;

    /* count # local unkowns */
    n = 0;
    for (i = 0; i < ndof; i++) {
	assert(dofs[i]->type != DOF_ANALYTIC);
	n += dofs[i]->dim * g->nelem * DofNBas(dofs[i], NULL);
    }
    N = n;
#if USE_MPI
    if (g->nprocs > 1)
	MPI_Allreduce(&n, &N, 1, PHG_MPI_INT, MPI_SUM, g->comm);
#endif	/* USE_MPI */

    /* create MAP */
    map = phgMapCreateSimpleMap(g->comm, n, N);
    /* save DOFs list to map->dofs */
    phgFree(map->dofs);
    map->ndof = ndof;
    map->dofs = (void *)dofs;

    return map;
}

MAP *
phgMapCreate(DOF *u, ...)
{
    int ndof;
    DOF **dofs;

    GetVariableArgs(ndof, dofs, u, DOF *);

    return phgMapCreateN(ndof, dofs);
}

SOLVER *
phgSolverCreate(OEM_SOLVER *oem_solver, DOF *u, ...)
{
    MAP *map;
    MAT *mat;
    int ndof;
    DOF **dofs;
    SOLVER *solver;

    GetVariableArgs(ndof, dofs, u, DOF *);

    /* create MAP */
    map = phgMapCreateN(ndof, dofs);
    /* create MAT */
    mat = phgMapCreateMat(map, map);
    phgMapDestroy(&map);
    mat->handle_bdry_eqns = FALSE;

    /* create SOLVER */
    solver = phgSolverMat2Solver(oem_solver, mat);
    phgMatDestroy(&mat);
    phgVecDisassemble(solver->rhs);
    solver->rhs_updated = FALSE;

    return solver;
}

int
phgSolverSolve(SOLVER *solver, BOOLEAN destroy, DOF *u, ...)
{
    GRID *g = u->g;
    DOF **dofs;
    VEC *v = phgMapCreateVec(solver->rhs->map, 1);
    FLOAT *data;
    int i, ndof;
    INT n;

    GetVariableArgs(ndof, dofs, u, DOF *);

    /* copy DOF to v as initial solution */
    data = v->data;
    for (i = 0; i < ndof; i++) {
	n = dofs[i]->dim * g->nelem * DofNBas(dofs[i], NULL);
	memcpy(data, dofs[i]->data, n * sizeof(*data));
	data += n;
    }

    phgSolverVecSolve(solver, destroy, v);

    /* copy solution in v to DOF */
    data = v->data;
    for (i = 0; i < ndof; i++) {
	n = dofs[i]->dim * g->nelem * DofNBas(dofs[i], NULL);
	memcpy(dofs[i]->data, data, n * sizeof(*data));
	data += n;
	phgDofUpdateGhost(dofs[i]);
    }

    phgVecDestroy(&v);

    phgFree(dofs);

    return 0;
}

static void
quad_interface_setup(DOF *ls, DOF *ls_grad, ELEMENT *e,
			int *ls_order, FUNC_3D *ls_func, FUNC_3D *grad_func)
{
   assert(ls != NULL || ls_grad != NULL);
   dummy_ctx.e = e;
   dummy_ctx.ls = ls;
   dummy_ctx.grad = ls_grad;
   if (ls != NULL)
      *ls_order = DofTypeOrder(ls, e);
   else
      *ls_order = DofTypeOrder(ls_grad, e) +
			(DofTypeOrder(ls_grad, e) < 0 ? 0 : 1);
   *ls_func = ls == NULL ? NULL : dummy_ls;
   *grad_func = ls_grad == NULL ? NULL : dummy_grad;
}

void
phgQuadInterface(DOF *ls, DOF *ls_grad, ELEMENT *e, int order,
		 FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
{
   FUNC_3D ls_func, grad_func;
   int ls_order;

   quad_interface_setup(ls, ls_grad, e, &ls_order, &ls_func, &grad_func);
   phgQuadInterfaceCuboid(ls_func, ls_order, grad_func,
		phgGeomGetCorners((ls != NULL ? ls : ls_grad)->g, e, NULL),
		order, rule_m, rule_0, rule_p);
}

void
phgQuadInterfaceFace(DOF *ls, DOF *ls_grad, ELEMENT *e, int face, int order,
		     FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
{
   FUNC_3D ls_func, grad_func;
   int ls_order;

   quad_interface_setup(ls, ls_grad, e, &ls_order, &ls_func, &grad_func);
   phgQuadInterfaceCuboidFace(ls_func, ls_order, grad_func,
		phgGeomGetCorners((ls != NULL ? ls : ls_grad)->g, e, NULL),
		face, order, rule_m, rule_0, rule_p);
}

FLOAT *
phgQuadGetRule_(GRID *g, ELEMENT *e, int which, int order)
/* returns a rule of order "order" for integration in the element "e", or
 * one of its face or edge. "which" specifies which part of the element is 
 * the integration domain:
 * 		which = 3 * no + (dim - 1)
 * where "no" is the face ("dim" == 2) or edge ("dim" == 1) number. */
{
    int dim = which % 3 + 1, no = which / 3;
    FLOAT *rule = NULL;

    if (dim == 3)
	phgQuadInterfaceCuboid(NULL, 0, NULL, e->corners, order,
							&rule, NULL, NULL);
    else if (dim == 2)
	phgQuadInterfaceCuboidFace(NULL, 0, NULL, e->corners, no, order,
							&rule, NULL, NULL);
    else
	phgError(1, "%s: unimplemented case (dim == 1)\n");

    return rule;
}

GRID *
phgCreateSimpleGrid(INT nelem, FLOAT (*list_of_corners)[2][Dim])
{
    INT i;
    ELEMENT *e;
    GRID *g = phgNewGrid(-1);

    g->comm = MPI_COMM_SELF;
    g->nprocs = 1;
    g->rank = 0;
    g->balanced = TRUE;

    g->nelem_global = g->nelem = nelem;
    if (g->partition == NULL)
	g->partition = phgAlloc((g->nprocs + 1) * sizeof(*g->partition));
    g->partition[0] = 0;
    g->partition[1] = nelem;
    g->lif = 1.0;

    g->roots = phgCalloc(nelem, sizeof(*g->roots));
    g->nroot = nelem;
    g->elems = phgAlloc(nelem * sizeof(*g->elems));
    for (i = 0; i < nelem; i++) {
	g->elems[i] = e = g->roots + i;
	e->elem_type = ET_CUBOID;
	e->index = i;
	memcpy(e->corners, list_of_corners[i], sizeof(e->corners));
    }

    return g;
}

/*----------------------- End of PHG_to_p8est wrapper -----------------------*/

#include "../src/solver.c"	/* for phgSolverDumpMATLAB */

#endif	/* HAVE_P4EST */
