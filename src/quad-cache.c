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

/* $Id: quad-cache.c,v 1.58 2022/06/18 09:25:54 zlb Exp $ */

/* These functions are used to cache (precompute) values of functions at the
 * quadrature points of a given numerical quadrature rule, in order to avoid
 * recompute the same functions multiple times.
 */

#include "phg.h"
#include "string.h"

/*============================ PHG-specific functions ========================*/

#ifndef PHG_TO_P4EST	/* Note: this allows to include this file inline */

static void
to_ref(QCACHE *qc, INT eno, FLOAT x[], FLOAT xref[])
{
    DOF *u = qc->fe;
    GRID *g = u->g;
    ELEMENT *e = g->elems[eno];

    phgGeomXYZ2Lambda(g, e, x[0], x[1], x[2], xref);
}

static void
from_ref(QCACHE *qc, INT eno, FLOAT xref[], FLOAT x[])
{
    DOF *u = qc->fe;
    GRID *g = u->g;
    ELEMENT *e = g->elems[eno];

    phgGeomLambda2XYZ(g, e, xref, x + 0, x + 1, x + 2);
}

/*------------------------- The QC interface for PHG -------------------------*/

static int
b_dim(QCACHE *qc) { return qc->fe == NULL ? 0 : DofTypeDim((DOF *)qc->fe); }

static int
f_dim(QCACHE *qc, void *fe) { return DofDim((DOF *)fe); }

static int
n_bas(QCACHE *qc, INT eno)
{ DOF *u = qc->fe; return u == NULL ? 0 : DofNBas(u, u->g->elems[eno]); }

static void
bas(QCACHE *qc, INT eno, const FLOAT lam[], FLOAT *buffer, int stride)
{
    DOF *u = qc->fe;
    GRID *g = u->g;
    ELEMENT *e = g->elems[eno];
    DOF_TYPE *t;
    const FLOAT *p;
    int dim = b_dim(qc), nb = n_bas(qc, eno);
    int j;

    t = u->hp == NULL ? u->type : u->hp->info->types[u->hp->max_order];
    p = t->BasFuncs(u, e, 0, -1, lam);
    for (j = 0; j < nb; j++, buffer += stride, p += dim)
	memcpy(buffer, p, dim * sizeof(*p));
}

static void
grd(QCACHE *qc, INT eno, const FLOAT lam[], FLOAT *buffer, int stride)
{
    DOF *u = qc->fe;
    GRID *g = u->g;
    ELEMENT *e = g->elems[eno];
    DOF_TYPE *t;
    const FLOAT *p, *J = phgGeomGetJacobian(g, e);
    int dim = b_dim(qc), nb = n_bas(qc, eno);
    int j, k;

    t = u->hp == NULL ? u->type : u->hp->info->types[u->hp->max_order];
    p = t->BasGrads(u, e, 0, -1, lam);
    /* convert D/Dlambda to D/Dxyz */
    for (j = 0; j < nb; j++, buffer += stride)
	for (k = 0; k < dim; k++, p += Dim + 1) {
	    buffer[k * Dim + 0] = p[0] * J[0 * (Dim + 1) + 0] +
				  p[1] * J[1 * (Dim + 1) + 0] +
				  p[2] * J[2 * (Dim + 1) + 0] +
				  p[3] * J[3 * (Dim + 1) + 0];
	    buffer[k * Dim + 1] = p[0] * J[0 * (Dim + 1) + 1] +
				  p[1] * J[1 * (Dim + 1) + 1] +
				  p[2] * J[2 * (Dim + 1) + 1] +
				  p[3] * J[3 * (Dim + 1) + 1];
	    buffer[k * Dim + 2] = p[0] * J[0 * (Dim + 1) + 2] +
				  p[1] * J[1 * (Dim + 1) + 2] +
				  p[2] * J[2 * (Dim + 1) + 2] +
				  p[3] * J[3 * (Dim + 1) + 2];
	}
}

static void
fun(QCACHE *qc, void *fe, INT eno, const FLOAT lam[], FLOAT *buffer)
{
    DOF *u = fe;
    GRID *g = u->g;
    ELEMENT *e = g->elems[eno];

    phgDofEval(fe, e, lam, buffer);
}

QDESC QD_PHG_ = {"PHG", to_ref, from_ref, b_dim, f_dim, n_bas, bas, grd, fun,
		 /*dim=*/Dim, /*ref_dim=*/Dim + 1};

/*----------------------------------------------------------------------------*/

/*======================= end of PHG-specific functions =====================*/

#endif	/* PHG_TO_P4EST */

int Q_DIV = -1, Q_CURL = -1;	/* assigned in phgQCNew() */

static void
qc_reset_data(QCACHE *qc, int index, int thread)
/* index specifies which element's data (0, 1 or -1, -1 means all elements)
 * thread specifies reset data for which thread, -1 mean all threads */
{
    int i, i0, i1, k, k0, k1, j;

    qc->in_use = TRUE;

    i0 = index < 0 ? 0 : index;
    i1 = index < 0 ? NEC : index + 1;

    k0 = thread < 0 ? 0 : thread;
    k1 = thread < 0 ? phgMaxThreads : thread + 1;

    for (k = k0; k < k1; k++) {
	for (i = i0; i < i1; i++) {
#define FreeData(o) \
	{ if (qc->cache[k].o.owner[i]) phgFree(qc->cache[k].o.data[i]); \
	  phgFree(qc->cache[k].o.dot_n[i]); \
	  phgFree(qc->cache[k].o.cross_n[i]); \
	  qc->cache[k].o.data[i] = qc->cache[k].o.dot_n[i] \
				 = qc->cache[k].o.cross_n[i] = NULL; }
	    FreeData(bas)
	    FreeData(grad)
	    if (qc->cache[k].func != NULL)
		for (j = 0; j < qc->n; j++)
		    FreeData(func[j]);
#undef FreeData
	}
    }
}

#define Thread(qc)	(qc)->cache[phgThreadId]

static int
qc_set_element_face(QCACHE *qc, INT eno, int face)
/* prepares data cache for a given element and face, clear previous cached
 * elements when necessary, and returns the cache index for the element/face.
 *
 * If eno < 0 then reset all caches (rule change or final reset) */
{
    int index;

    qc->in_use = TRUE;

    if (eno < 0) {
	qc_reset_data(qc, -1, phgThreadId);
	for (index = 0; index < NEC; index++)
	    Thread(qc).eno[index] = -1;
	return -1;
    }

    for (index = 0; index < NEC; index++)
	if (Thread(qc).eno[index] == eno && Thread(qc).face[index] == face)
	    return index;

    /* find an unused slot */
    for (index = 0; index < NEC; index++)
	if (Thread(qc).eno[index] < 0)
	    break;
    if (index >= NEC) {
	/* no unused slot, discard one */
	index = Thread(qc).index0;
	if (++Thread(qc).index0 >= NEC)
	    Thread(qc).index0 = 0;
    }

    if (Thread(qc).eno[index] >= 0)
	qc_reset_data(qc, index, phgThreadId);
    Thread(qc).eno[index] = eno;
    Thread(qc).face[index] = face;
    Thread(qc).nb[index] = eno < 0 ? 0 : qc->qd->n_bas(qc, eno);

    return index;
}

static int
Op_dot(int cdim, const FLOAT *coeff, int dim, const FLOAT *in, FLOAT *out,
       QCACHE *qc, int iq, const void *ctx)
/* applies projection to data, returns dimension of the projected data */
{
    int k, out_dim;
    const FLOAT *nv = phgQCGetNormal(qc, iq, NULL);

    if (nv == NULL)
	phgError(1, "%s: normal vector not available.\n", __func__);

    if (dim % Dim != 0)
	phgError(1, "%s: invalid dimension (%d)!.\n", __func__, dim);

    out_dim = dim / Dim;

    if (in == NULL)
	return out_dim;

    for (k = 0; k < dim; k += Dim, in += Dim)
	*(out++) = in[0] * nv[0] + in[1] * nv[1] + in[2] * nv[2];

    return out_dim;
}

static int
Op_mult(int cdim, const FLOAT *coeff, int dim, const FLOAT *in, FLOAT *out,
        QCACHE *qc, int iq, const void *ctx)
/* multiplies coefficient with function. This is the default operator. */
{
    int out_dim = dim, i, j;

    if (cdim == 0) {
	return dim;
    }
    else if (cdim == 1) {
	out_dim = dim;
	if (in == NULL)
	    return out_dim;
	for (i = 0; i < dim; i++)
	    *(out++) = *(in++) * *coeff;
    }
    else if (cdim == dim) {
	out_dim = dim;
	if (in == NULL)
	    return out_dim;
	for (i = 0; i < dim; i++)
	    *(out++) = *(in++) * *(coeff++);
    }
    else if (cdim == dim * dim) {
	out_dim = dim;
	if (in == NULL)
	    return out_dim;
	for (i = 0; i < dim; i++, out++)
	    for (j = 0, *out = 0.0; j < dim; j++)
		*out += in[j] * coeff[i * dim + j];
    }
    else if (dim == 1) {
	out_dim = cdim;
	if (in == NULL)
	    return out_dim;
	for (i = 0; i < cdim; i++)
	    *(out++) = *in * *(coeff++);
    }
    else {
	phgError(1, "%s: inconsistent dimensions: dim = %d, cdim = %d\n",
		    __func__, dim, cdim);
    }

    return out_dim;
}

static int
Op_cross(int cdim, const FLOAT *coeff, int dim, const FLOAT *in, FLOAT *out,
         QCACHE *qc, int iq, const void *ctx)
/* applies projection to data, returns dimension of the projected data */
{
    int k;
    const FLOAT *nv = phgQCGetNormal(qc, iq, NULL);

    if (nv == NULL)
	phgError(1, "%s: normal vector not available.\n", __func__);

    if (dim % Dim != 0)
	phgError(1, "%s: invalid dimension (%d)!.\n", __func__, dim);

    if (in == NULL)
	return dim;

    for (k = 0; k < dim; k += Dim, in += Dim) {
	*(out++) = in[1] * nv[2] - in[2] * nv[1];
	*(out++) = in[2] * nv[0] - in[0] * nv[2];
	*(out++) = in[0] * nv[1] - in[1] * nv[0];
    }

    return dim;
}

static int
Op_div(int cdim, const FLOAT *coeff, int dim, const FLOAT *in, FLOAT *out,
       QCACHE *qc, int iq, const void *ctx)
{
    if (in != NULL) {
	assert(Dim == 3 && dim == 9);
	assert(out != NULL);
	*out = in[0 * Dim + 0] + in[1 * Dim + 1] + in[2 * Dim + 2];
    }
    return 1;
}

static int
Op_curl(int cdim, const FLOAT *coeff, int dim, const FLOAT *in, FLOAT *out,
        QCACHE *qc, int iq, const void *ctx)
{
    if (in != NULL) {
	if (!(Dim == 3 && dim == 9)) {
	    phgError(1, "%s: invalid dimension (%d) of input data.\n",
			__func__, dim);
	}
	assert(out != NULL);
	out[0] = in[2 * Dim + 1] - in[1 * Dim + 2]; /* du2/dx1 - du1/dx2 */
	out[1] = in[0 * Dim + 2] - in[2 * Dim + 0]; /* du0/dx2 - du2/dx0 */
	out[2] = in[1 * Dim + 0] - in[0 * Dim + 1]; /* du1/dx0 - du0/dx1 */
    }
    return 3;
}

int
phgQCAddFunction_(QCACHE *qc,
		  void *f_fe, FLOAT *f_data, FUNC_3D f_xyz, int f_fid, int dim,
		  OP_FUNC op, int base_fid, const void *ctx)
/* Adds (attaches) a function to the QCache object, returns the function no.
 * "dim" specifies the dimension of output data for "op", and the
 * dimension (number of components) of "f_{data,xyz}".
 *
 * If "base_fid" is not equal to Q_NONE, then the new function is obtained
 * by multiplying the function with "base_fid" (a coefficient function).
 * */
{
    int i;

    /* "base_fid" should be a valid FID */
    if (!(base_fid >= Q_NONE && base_fid < qc->n)) {
	phgError(1, "%s:%d: invalid arguments.\n", __FILE__, __LINE__);
    }

    /* "f_fid" should be a valid FID */
    if (!(f_fid >= Q_NONE && f_fid < qc->n)) {
	phgError(1, "%s:%d: invalid arguments.\n", __FILE__, __LINE__);
    }

    /* "f_fid" is only for phgQCAdd*Coefficient ("base_fid" != Q_NONE) */
    if (!(f_fid == Q_NONE || base_fid != Q_NONE)) {
	phgError(1, "%s:%d: invalid arguments.\n", __FILE__, __LINE__);
    }

    /* "op" requires "base_fid" */
    if (!(op == NULL || base_fid != Q_NONE)) {
	phgError(1, "%s:%d: invalid arguments.\n", __FILE__, __LINE__);
    }

    /* exactly one of the function pointers should be non NULL */
    i = (f_fe != NULL) + (f_data != NULL) + (f_xyz != NULL) + (f_fid != Q_NONE);
    if (!(i == 1 || (i == 0 && op != NULL))) {
	phgError(1, "%s:%d: invalid arguments.\n", __FILE__, __LINE__);
    }
 
    if (qc->in_use)
	phgError(1, "phgQCAdd*Function can only be called before any other "
		    "phgQC* functions (except phgQCNew).\n");

    if (qc->n >= qc->alloc) {
	/* reserve and initialize pointers for more functions */
	void *old = qc->f;
	qc->alloc += 4;
	qc->f = phgCalloc(qc->alloc, sizeof(*qc->f));
	if (qc->n > 0)
	    memcpy(qc->f, old, qc->n * sizeof(*qc->f));
	phgFree(old);
    }

    if (f_data != NULL && dim > 0)
	memcpy(qc->f[qc->n].f_data = phgAlloc(dim * sizeof(*f_data)), f_data,
		dim * sizeof(*f_data));
    else
	qc->f[qc->n].f_data = NULL;

    qc->f[qc->n].f_fe = f_fe;
    qc->f[qc->n].f_xyz = f_xyz;
    qc->f[qc->n].f_fid = f_fid;
    qc->f[qc->n].op = op;
    qc->f[qc->n].ctx = ctx;
    if (op != NULL) {
	int base_dim = 0;
	if (base_fid >= 0)
	    base_dim = qc->f[base_fid].f_dim;
	else if (base_fid == Q_BAS)
	    base_dim = qc->qd->b_dim(qc);
	else if (base_fid == Q_GRAD)
	    base_dim = qc->qd->b_dim(qc) * Dim;
	else
	    phgError(1, "%s:%d: invalid base_fid (%d)!\n",
				__FILE__, __LINE__, base_fid);
	qc->f[qc->n].f_dim = op(dim, NULL, base_dim, NULL, NULL, qc, 0, ctx);
    }
    else if (f_fe != NULL)
	qc->f[qc->n].f_dim = qc->qd->f_dim(qc, f_fe);
    else if (f_fid != Q_NONE)
	qc->f[qc->n].f_dim = qc->f[f_fid].f_dim;
    else
	qc->f[qc->n].f_dim = dim;
    qc->f[qc->n].base_fid = base_fid;

    /* set const_flag */
    if (base_fid == Q_NONE)
	qc->f[qc->n].const_flag = (f_data != NULL);
#if 1
    /* the data is always considered non constant after applying an operator,
     * since the latter may depend on (x,y,z) or nv
     *
     * TODO: to save computations and memory, modify OP_FUNC to return an int
     * flag indicating whether the resulting function is constant when both
     * the coefficient and the base_fid are constant, with 3 possibilities:
     *		0: constant, 1: non constant, 2: constant if nv is constant
     */
    else
	qc->f[qc->n].const_flag = FALSE;
#else
    else if (base_fid < 0)
	qc->f[qc->n].const_flag = FALSE;
    else
	qc->f[qc->n].const_flag = qc->f[base_fid].const_flag;
#endif

    /* use Op_mult() for multiplying by a coefficient function */
    if (qc->f[qc->n].op == NULL)
	qc->f[qc->n].op = Op_mult;

    return qc->n++;
}

int
phgQCAddProjection(QCACHE *qc, PROJ proj, int fid)
{
    switch (proj) {
	case PROJ_NONE:
	    return fid;	/* nothing to do! */
	case PROJ_DOT:
	    return phgQCAddOperator(qc, Op_dot, fid);
	case PROJ_CROSS:
	    return phgQCAddOperator(qc, Op_cross, fid);
    }

    return Q_NONE;	/* to avoid gcc warning */
}

void
phgQCSetRule_(QCACHE *qc, FLOAT scale, FLOAT *rule, QUAD *quad)
/* Note: the quadrature rule is in reference coordinates iff scale >= 0.0 */
{
    int dim;
    FLOAT *nv = NULL;

    qc->in_use = TRUE;

    assert(rule == NULL || quad == NULL);

    /* invalidate cached data */
    qc_set_element_face(qc, -1, -1);

    Thread(qc).rule = rule;
    Thread(qc).quad = quad;
    Thread(qc).scale = scale < 0.0 ? 1.0 : scale;
    Thread(qc).ref_flag = (scale >= 0.0);

    if (quad != NULL) {
	Thread(qc).np = quad->npoints;
	if (qc->qd->dim != quad->dim)
	    phgError(1, "Dimensions of QDESC (%d) and quad (%d) mismatch!\n",
			qc->qd->dim, quad->dim);
    }
    else if (rule != NULL) {
	Thread(qc).np = phgQuadInterfaceRuleInfo(rule, &dim, NULL, &nv);
	if (dim != qc->qd->dim)
	    phgError(1, "Dimensions of QDESC (%d) and rule (%d) mismatch!\n",
			qc->qd->dim, dim);
    }
    else {
	Thread(qc).np = 0;
    }

    if (nv != NULL && Thread(qc).nv != NULL) {
	phgFree(Thread(qc).nv);
	Thread(qc).nv = NULL;
    }
}

void
phgQCSetConstantNormal(QCACHE *qc, const FLOAT *nv)
/* If nv != NULL, then the face is flat and nv gives its unit normal vector,
 * otherwise the quadrature rule must provide the list of unit normal vectors
 * at the quadratuee points when needed.
 */
{
    assert(nv != NULL);

    qc->in_use = TRUE;

    if (Thread(qc).nv == NULL ||
	memcmp(nv, Thread(qc).nv, Dim * sizeof(*Thread(qc).nv)) != 0) {
	/* reset projected data */
	int i, j;
#define ResetProjectedData(f) \
	for (j = 0; j < NEC; j++) { \
	    assert(Thread(qc).f.owner[j] || Thread(qc).f.data[j] == NULL); \
	    phgFree(Thread(qc).f.data[j]); \
	    phgFree(Thread(qc).f.dot_n[j]); \
	    phgFree(Thread(qc).f.cross_n[j]); \
	    Thread(qc).f.data[j] = Thread(qc).f.dot_n[j] \
				 = Thread(qc).f.cross_n[j] = NULL; \
	}
	ResetProjectedData(bas)
	ResetProjectedData(grad)
	if (Thread(qc).func != NULL)
	    for (i = 0; i < qc->n; i++)
		if (qc->f[i].op == Op_dot || qc->f[i].op == Op_cross)
		    ResetProjectedData(func[i])
#undef ResetProjectedData
    }

    if (Thread(qc).nv == NULL)
	Thread(qc).nv = phgAlloc(Dim * sizeof(*Thread(qc).nv));
    memcpy(Thread(qc).nv, nv, Dim * sizeof(*Thread(qc).nv));
}

QCACHE *
phgQCNew(QDESC *qd, void *fe)
{
    QCACHE *qc = phgCalloc(1, sizeof(*qc));

    qc->qd = qd;
    qc->fe = fe;
    qc->cache = phgCalloc(phgMaxThreads, sizeof(*qc->cache));

    /* create Q_DIV and Q_CURL operators */
    Q_DIV = phgQCAddOperator(qc, Op_div, Q_GRAD);
    Q_CURL = phgQCAddOperator(qc, Op_curl, Q_GRAD);

    return qc;
}

QCACHE *
phgQCClone(QCACHE *qc0)
/* returns a QCACHE object with the same FE space, functions, coefficients,
 * operators, projections, and constant normal vector */
{
    int i;

    QCACHE *qc= phgQCNew(qc0->qd, qc0->fe);

    /* add functions, coefficients, operators, and projections */
    for (i = qc->n; i < qc0->n; i++)
	phgQCAddFunction_(qc, qc0->f[i].f_fe, qc0->f[i].f_data,
			  qc0->f[i].f_xyz, qc0->f[i].f_fid, qc0->f[i].f_dim,
			  qc0->f[i].op, qc0->f[i].base_fid, qc0->f[i].ctx);

    /* set normal vector */
    if (Thread(qc0).nv != NULL)
	phgQCSetConstantNormal(qc, Thread(qc0).nv);

    return qc;
}

/* prototype for qc_get_data() */
static int qc_get_data(QCACHE *qc, INT eno, int face, int fid, const int bno,
		PROJ proj, FLOAT **data, BOOLEAN *const_flag, BOOLEAN *owner);

static FLOAT *
ref_to_lambda(QCACHE *qc, FLOAT xref[], FLOAT lambda[])
/* (for QD_PHG only) converts reference coord to barycentric coord */
{
    int i, dim = qc->qd->dim;
    FLOAT x = 1.0;

    for (i = 0; i < dim; i++)
	x -= (lambda[i] = xref[i]);
    lambda[dim] = x;

    return lambda;
}

static int
qc_get_nb(QCACHE *qc, int fid, int index)
/* returns nb for "fid", depending on whether "fid" is a function or a basis */
{
    while (fid >= 0)
	fid = qc->f[fid].f_fid;
    return fid == Q_NONE ? 1 : Thread(qc).nb[index];	
}

static FLOAT *
qc_init_values(QCACHE *qc, INT eno, int face, int fid, int index,
		BOOLEAN *const_flag, BOOLEAN *owner_flag)
/* computes the values at the quadrature points, ignoring base_fid.
 *
 * The values are returned in the follwing format:
 * 1. for fid==Q_BAS, where dim==qc->qd->b_dim(qc):
 * 	FLOAT[nbas][npoints][dim]
 * 2. for fid==Q_GRAD, where dim==qc->qd->b_dim(qc):
 * 	FLOAT[nbas][npoints][dim][Dim]
 * 3. for fid==Q_FUNC, where dim==dimension of the function:
 * 	*const_flag == FALSE:	FLOAT[npoints][dim]
 *	*const_flag == TRUE:	FLOAT[dim]
 *
 * "*owner_flag" indicates whether the data is owned by the calling function */
{
#ifndef PHG_TO_P4EST
    DOF *u = NULL;	/* for QD_PHG only */
#endif	/* PHG_TO_P4EST */
    FLOAT *buffer, *p, *q;
    int nb, dim, i, qdim, fi = 0;

    *const_flag = FALSE;

#if 1
    assert(eno >= 0);
#else
    if (/*Thread(qc).np == 0 ||*/ eno < 0)
	return NULL;
#endif

    assert(qc->qd->dim == Dim);	/* TODO: 2D case */

    if (fid >= 0) {
	fi = fid;
	fid = Q_FUNC;
    }

#ifndef PHG_TO_P4EST
    if (qc->qd == QD_PHG) {
	if ((u = qc->fe) == NULL) {
	    assert(fid == Q_FUNC && qc->f[fi].f_fe != NULL);
	    u = qc->f[fi].f_fe;
	}
    }
#endif	/* PHG_TO_P4EST */

#define Me (Thread(qc).me)
    switch (fid) {
	case Q_FUNC:
		Me = &Thread(qc).func[fi];
		*const_flag = qc->f[fi].const_flag;
		break;
	case Q_BAS:
		Me = &Thread(qc).bas;
		break;
	case Q_GRAD:
		Me = &Thread(qc).grad;
		break;
	default:
		phgError(1, "%s:%d: invalid FID %d\n", fid);
    }

    /* Check whether we can use existing data */
    if (Me->data[index] != NULL) {
	*owner_flag = FALSE;
	return Me->data[index];
    }
#undef Me

    i = (qc->f[fi].f_fe != NULL) + (qc->f[fi].f_xyz  != NULL) +
	(qc->f[fi].f_data != NULL) + (qc->f[fi].f_fid != Q_NONE);
    if (fid == Q_FUNC && i == 0) {
	/* Note: this is the case op(void, f2) */
	*owner_flag = FALSE;
	return NULL;
    }
    assert(fid != Q_FUNC || i == 1);

#ifndef PHG_TO_P4EST
#if 0	/* TODO! */
    if (Thread(qc).rule == NULL) {
	ELEMENT *e = u->g->elems[eno];
	assert(Thread(qc).quad != NULL && qc->qd == QD_PHG);
	/* use cached data from phgQuadGet{Basis,Dof}{Values,Gradient}, this
	 * has the advantage of taking into account type->invariant to avoid
	 * unnecessary recomputation of basis functions (for PHG only!). */
	*owner_flag = FALSE;
	if (fid == Q_BAS)
	    return (FLOAT *)(Thread(qc).dim == 3 ?
			phgQuadGetBasisValues(e, u, 0, Thread(qc).quad) :
			phgQuadFaceGetBasisValues(e, face, u, 0,
						  Thread(qc).quad));
	else if (fid == Q_GRAD)
	    return (FLOAT *)(Thread(qc).dim == 3 ?
			phgQuadGetBasisGradient(e, u, 0, Thread(qc).quad) :
			phgQuadFaceGetBasisGradient(e, face, u, 0,
						    Thread(qc).quad));
	else if (qc->f[fi].f_fe != NULL)
	    return (FLOAT *)(Thread(qc).dim == 3 ?
			phgQuadGetDofValues(e, qc->f[fi].f_fe,
					    Thread(qc).quad) :
			phgQuadFaceGetDofValues(e, face, qc->f[fi].f_fe,
						Thread(qc).quad));
#if 0
	else if (qc->f[fi].f_xyz != NULL)
	    return /* TODO: missing phgQuadFaceGetFuncValues function! */
#endif
    }
#endif
#endif	/* PHG_TO_P4EST */

    if (fid != Q_FUNC) {
	nb = qc->qd->n_bas(qc, eno);
	dim = qc->qd->b_dim(qc);
	if (fid == Q_GRAD)
	    dim *= Dim;
    }
    else {
	if (qc->f[fi].f_data != NULL) {
	    *owner_flag = FALSE;
	    *const_flag = TRUE;
	    return qc->f[fi].f_data;
	}
	nb = qc_get_nb(qc, fi, index);
	dim = qc->f[fi].f_dim;

#if 1	/* for PHG only! */
# define CHECK_CONSTANT_DOF 1
	if (qc->qd == QD_PHG && qc->f[fi].f_fe != NULL &&
	    DofIsConstant((DOF *)qc->f[fi].f_fe)) {
	    *owner_flag = FALSE;
	    *const_flag = TRUE;
	    return ((DOF *)qc->f[fi].f_fe)->data;
	}
#else
# define CHECK_CONSTANT_DOF 0
#endif
    }

    buffer = phgAlloc(nb * (size_t)Thread(qc).np * dim * sizeof(*buffer));
    *owner_flag = TRUE;

    if (qc->f[fi].f_fid != Q_NONE) {
	int bno;
	BOOLEAN const_flag0, owned;
	FLOAT *data;
	if (dim == 0 || Thread(qc).np == 0 || nb == 0) {
	    return NULL;
	}
	for (bno = 0; bno < nb; bno++) {
	    qc_get_data(qc, eno, face, qc->f[fi].f_fid, bno,
			PROJ_NONE, &data, &const_flag0, &owned);
	    memcpy(buffer + Thread(qc).np * dim * bno, data,
		   Thread(qc).np * dim * sizeof(*data));
	    if (owned)
		phgFree(data);
	}
	return buffer;
    }

    if (Thread(qc).quad != NULL) {
	q = Thread(qc).quad->points;
	qdim = Dim + 1;
    }
    else {
	phgQuadInterfaceRuleInfo(Thread(qc).rule, &qdim, &q, NULL);
	qdim += 1;
    }

    assert(qc->qd->ref_dim >= qc->qd->dim);	/* x0[] size */
    for (i = 0; i < Thread(qc).np; i++, q += qdim) {
	FLOAT x0[qc->qd->ref_dim], *x;
	p = buffer + dim * i;
	switch (fid) {
	    case Q_BAS:
		if (qc->qd->bas == NULL)
		    phgError(1, "%s: member function bas() undefined.\n",
				qc->qd->name);
		if (Thread(qc).ref_flag)
		    x = qdim == qc->qd->ref_dim ? q : ref_to_lambda(qc, q, x0);
		else
		    qc->qd->to_ref(qc, eno, q, x = x0);
		(qc->bas_hook != NULL ? qc->bas_hook : qc->qd->bas)
					(qc, eno, x, p, Thread(qc).np * dim);
		break;
	    case Q_GRAD:
		if (qc->qd->grd == NULL)
		    phgError(1, "%s: member function grd() undefined.\n",
				qc->qd->name);
		if (Thread(qc).ref_flag)
		    x = qdim == qc->qd->ref_dim ? q : ref_to_lambda(qc, q, x0);
		else
		    qc->qd->to_ref(qc, eno, q, x = x0);
		(qc->grd_hook != NULL ? qc->grd_hook : qc->qd->grd)
					(qc, eno, x, p, Thread(qc).np * dim);
		break;
	    case Q_FUNC:
		if (qc->f[fi].f_fe != NULL) {
		    if (qc->qd->fun == NULL)
			phgError(1, "%s: member function fun() undefined.\n",
				    qc->qd->name);
		    if (Thread(qc).ref_flag)
			x = qdim == qc->qd->ref_dim ? q :
						      ref_to_lambda(qc, q, x0);
		    else
			qc->qd->to_ref(qc, eno, q, x = x0);
		    qc->qd->fun(qc, qc->f[fi].f_fe, eno, x, p);
		}
		else if (qc->f[fi].f_xyz != NULL) {
		    if (Thread(qc).ref_flag) {
			FLOAT x1[qc->qd->ref_dim];
			x = qdim == qc->qd->ref_dim ? q :
						      ref_to_lambda(qc, q, x1);
			qc->qd->from_ref(qc, eno, x, x0);
			x = x0;
		    }
		    else {
			x = q;
		    }
		    qc->f[fi].f_xyz(x[0], x[1], x[2], p);
		}
		else {
		    phgError(1, "%s:%d: unexpected!\n", __FILE__, __LINE__);
		}
		break;
	}
    }

    return buffer;
}

static int
qc_apply_op(QCACHE *qc, OP_FUNC op, int nb, int np,
	    int cdim, const FLOAT *coeff, BOOLEAN const_c,
	    int dim, const FLOAT *in_data, BOOLEAN const_d,
	    FLOAT **out_data, const void *ctx)
/* applies "op" to "data", results returned in *out (dynamically allocated),
 * the return value is the dimension of transformed data */
{
    int i, j, inc_c, out_dim;
    FLOAT *p;
    const FLOAT *q = in_data, *c;

    inc_c = const_c ? 0 : cdim;
    out_dim = op(cdim, NULL, dim, NULL, NULL, qc, 0, ctx);

    p = *out_data = phgAlloc(out_dim * (size_t)nb * np * sizeof(*p));

    for (j = 0; j < nb; j++)
	for (i = 0, c = coeff; i < np; i++, p += out_dim, q += dim, c += inc_c)
	    op(cdim, c, dim, q, p, qc, i, ctx);

    return out_dim;
}

static int
qc_get_data(QCACHE *qc, INT eno, int face, int fid, const int bno,
	    PROJ proj, FLOAT **data, BOOLEAN *const_flag, BOOLEAN *owner_flag)
/* sets "*data" to point to the data for the attached function "fid" 
 * (values of the given function at all quadrature points), and
 * returns data dimension (number of components of the function).
 * "fid" is the FID, "bno" is basis or function no.
 * "*const_flag" indicates whether data is constant (invariant in space).
 * "*owner_flag" indicates whether the data is owned by the calling function.
 */
{
    int dim = 0, index, np, fno;
    OP_FUNC op = NULL;

    /* make sure FID_min is correctly defined */
    assert(fid >= FID_min && fid < qc->n);

    if (fid >= 0) {
	fno = fid;
	fid = Q_FUNC;
    }
    else {
	fno = -1;	/* unused */
    }

    qc->in_use = TRUE;

    index = qc_set_element_face(qc, eno, face);
    np = Thread(qc).np;

#define Me (Thread(qc).me)
    switch (fid) {
	case Q_FUNC:
	    assert(fno >= 0 && fno < qc->n);
	    *const_flag = qc->f[fno].const_flag;
	    if (Thread(qc).func == NULL)
		Thread(qc).func = phgCalloc(qc->n, sizeof(*Thread(qc).func));
	    if (Thread(qc).func[fno].data[index] == NULL && np != 0) {
		int nb = 1, fdim = 0, dim = 0;
		FLOAT *fdata;
		BOOLEAN fowned, fconst, dconst;
		FLOAT *data = NULL;
		fdata = qc_init_values(qc, eno, face, fno, index,
							&fconst, &fowned);
		fdim = qc->f[fno].f_dim;
		/* check coefficient fid */
		fid = qc->f[fno].base_fid;
		if (fid != Q_NONE) {
		    /* this is coefficient function, get base function values */
		    BOOLEAN owned;
		    assert(qc_get_nb(qc, fno, index) == 1); /* coef not basis */
		    dim = qc_get_data(qc, eno, face, fid, 0, PROJ_NONE,
				      &data, &dconst, &owned);
		    assert(owned == FALSE);	/* PROJ_NONE => FALSE */
		    nb = fid < 0 ? Thread(qc).nb[index] :
				   Thread(qc).func[fid].nb[index];
		    assert(qc->f[fno].op != NULL);
                    /* apply op with coeff and data */
		    Thread(qc).func[fno].dim =
				qc_apply_op(qc, qc->f[fno].op, nb, np,
					    fdim, fdata, fconst,
					    dim, data, dconst, &data,
					    qc->f[fno].ctx);
		    Thread(qc).func[fno].data[index] = data;
		    *const_flag = FALSE;
		    if (fowned)
			phgFree(fdata);
		}
		else  {
		    /* ordinary function, return function values */
		    Thread(qc).func[fno].data[index] = fdata;
		    Thread(qc).func[fno].dim = fdim;
		}
		Thread(qc).func[fno].nb[index] = nb;
		Thread(qc).func[fno].owner[index] = TRUE;
	    }
	    Me = &Thread(qc).func[fno];
	    break;
	case Q_BAS:
	    *const_flag = FALSE;
	    if (Thread(qc).bas.data[index] == NULL && np != 0) {
		/* get values of all basis functions */
		Thread(qc).bas.data[index] =
			qc_init_values(qc, eno, face, Q_BAS, index, const_flag,
				       &Thread(qc).bas.owner[index]);
		Thread(qc).bas.dim = qc->qd->b_dim(qc);
		Thread(qc).bas.nb[index] = Thread(qc).nb[index];
		assert(Thread(qc).bas.data[index] != NULL);
	    }
	    Me = &Thread(qc).bas;
	    break;
	case Q_GRAD:
	    *const_flag = FALSE;
	    if (Thread(qc).grad.data[index] == NULL && np != 0) {
		/* get gradients of all basis functions */
		Thread(qc).grad.data[index] =
			qc_init_values(qc, eno, face, Q_GRAD, index, const_flag,
				       &Thread(qc).grad.owner[index]);
		Thread(qc).grad.dim = qc->qd->b_dim(qc) * Dim;
		Thread(qc).grad.nb[index] = Thread(qc).nb[index];
		assert(Thread(qc).grad.data[index] != NULL);
	    }
	    Me = &Thread(qc).grad;
	    break;
    }

    dim = Me->dim;
    *data = Me->data[index] + dim * (*const_flag ? 1 : np) * bno;

    *owner_flag = FALSE;
    if (proj != PROJ_NONE && dim > 0) {
	FLOAT **pdata = NULL;	/* projected data pointer */
	np = *const_flag ? 1 : Thread(qc).np;
	if (proj == PROJ_DOT) {
	    op = Op_dot;
	    pdata = &Me->dot_n[index];
	}
	else if (proj == PROJ_CROSS) {
	    op = Op_cross;
	    pdata = &Me->cross_n[index];
	}
	else {
	    phgError(1, "%s:%d: unexpected error!\n", __FILE__, __LINE__);
	}
#if 0
	/* don't cache projected data, recompute them each time.
	 * (the "dot_n" and "cross_n" members are not used in this case) */
	Unused(pdata);
	dim = qc_apply_op(qc, op, 1, np, 0, NULL, FALSE, dim, *data, data, ctx);
	*owner_flag = TRUE;
#else
	/* compute and cache projected data */
	dim = *pdata != NULL ?  op(0, NULL, dim, NULL, NULL,
				   qc, 0, fno >= 0 ? qc->f[fno].ctx : NULL) :
			qc_apply_op(qc, op, Me->nb[index], np, 0, NULL, FALSE,
				    dim, Me->data[index], *const_flag,
				    pdata, fno >= 0 ? qc->f[fno].ctx : NULL);
	*data = (*pdata) + dim * np * bno;
#endif
    }

#undef Me

    return dim;
}

FLOAT *
phgQCIntegrate_(const char *file, const int line,
	QCACHE *qc1, INT eno1, int face1, int fid1, PROJ proj1, int i1,
	QCACHE *qc2, INT eno2, int face2, int fid2, PROJ proj2, int i2,
	int M, int N, int K, FLOAT *res)
/* This is the general UI for computing integral of the production of two
 * given functions using the quadrature rule attached to "qc1" and "qc2".
 *
 * The first function is specified by "qc1" and "fid1", which is evaluated in
 * the element "eno1", and in the face "face1" of the element with projection
 * "proj1" when applicable. In the case of a basis function, "i1" specifies
 * the basis number.
 *
 * Similarly, the second function is specified by "qc2", "eno2", "face2",
 * "fid2", "proj2" and "i2".
 *
 * The dimensions of the two functions must be consistent.
 * "qc1" and "qc2", if unequal, must have the same number of quadrature points.
 *
 * The arguments "M", "N" and "K" specify matrix dimentions in a GEMM-like way.
 * The integrand is regarded as a matrix function which is the product of two
 * matrix functions. The 1st function is of dimension "MxK" with values stored
 * as "FLOAT[M][K]", the 2nd function is of dimension "KxN" with values stored
 * as "FLOAT[K][N]", and the integrand is of dimension "MxN" with values stored
 * as "FLOAT[M][N]".
 *
 * It is required that "M*K<=dim1" and "K*N<=dim2", where "dim1" and "dim2" are
 * respectively the total number of components in the 1st and the 2nd function.
 * 
 * Zero or negative values of "M", "N" and "K" are regarded as unspecified.
 * If some of the parameters "M", "N" and "K" are unspecified, then they are
 * calculated using "dim1", "dim2" and the specified parameters, in a
 * comprehensive way. For example, if "M" is unspecified and "K" is specified,
 * then "M:=dim1/K".
 *
 * If all "M", "N" and "K" are unspecified, then it is assumed that
 * "max(dim1,dim2)" is a multiple of "min(dim1,dim2)", and we set
 * "K:=min(dim1,dim2)", "M:=dim1/K" and "N:=dim2/K".
 *
 * In the case "M==1" and "N==1", "res" can be NULL and in this case it is set
 * to point to a static and threadprivate internal FLOAT variable.
 *
 * This function returns res, or &res0 if res == NULL.
 */
{
    int i, j, k, np, dim1, dim2, winc;
    FLOAT val = 0.0, d;
    FLOAT *data1, *data2, *p1, *p2, *pw;
    BOOLEAN const_flag1, const_flag2, owner1, owner2;
    static FLOAT res0;
#if USE_OMP
# pragma omp threadprivate(res0)
#endif	/* USE_OMP */

    if (Thread(qc1).quad != NULL) {
	if (Thread(qc1).quad != Thread(qc2).quad)
	    phgError(1, "%s:%d: quads mismatch!\n", file, line);
	pw = Thread(qc1).quad->weights;
	winc = 1;
    }
    else {
	if (Thread(qc1).rule != Thread(qc2).rule)
	    phgError(1, "%s:%d: rules mismatch!\n", file, line);
	phgQuadInterfaceRuleInfo(Thread(qc1).rule, &i, &pw, NULL);
	pw += i;
	winc = i + 1;
    }

    np = Thread(qc1).np;
    dim1 = qc_get_data(qc1, eno1, face1, fid1, i1, proj1,
					&data1, &const_flag1, &owner1);
    p1 = data1;
    dim2 = qc_get_data(qc2, eno2, face2, fid2, i2, proj2,
					&data2, &const_flag2, &owner2);
    p2 = data2;

    if (np == 0 || dim1 == 0 || dim2 == 0) {
	if (owner1)
	    phgFree(data1);
	if (owner2)
	    phgFree(data2);
	if (M * N > 1) {
	    assert(res != NULL);
	    memset(res, 0, M * N * sizeof(*res));
	}
	else {
	    if (M * N == 0) {
		static BOOLEAN warned = FALSE;
		if (!warned)
		    phgWarning("%s:%d: some values may be uninitialized.\n",
				file, line);
#if USE_OMP
# pragma omp critical(uninitialized)
#endif
		warned = TRUE;
	    }
	    if (res == NULL)
		res = &res0;
	    *res = 0.0;
	}
	return res;
    }

    if (M <= 0 && N <= 0 && K <= 0) {
	if (dim1 % dim2 != 0 && dim2 % dim1 != 0)
	    phgError(1, "%s:%d: inconsistent data dimensions (%d <-> %d)!\n",
		 	file, line, dim1, dim2);
	/* set K:=min(dim1,dim2) */
	K = dim1 < dim2 ? dim1 : dim2;
	M = dim1 / K;
	N = dim2 / K;
    }

    if (K <= 0) {
	if (M > 0) {
	    if (dim1 % M != 0)
		phgError(1, "%s:%d: inconsistent K (%d) and dim1 (%d)\n",
				file, line, K, dim1);
	    K = dim1 / M;
	}
	else {
	    assert(N > 0);
	    if (dim2 % N != 0)
		phgError(1, "%s:%d: inconsistent K (%d) and dim2 (%d)\n",
				file, line, K, dim2);
	    K = dim2 / N;
	}
    }

    if (M <= 0) {
	if (dim1 % K != 0)
	    phgError(1, "%s:%d: inconsistent K (%d) and dim1 (%d)\n",
			file, line, K, dim1);
	M = dim1 / K;
    }

    if (N <= 0) {
	if (dim2 % K != 0)
	    phgError(1, "%s:%d: inconsistent K (%d) and dim2 (%d)\n",
			file, line, K, dim2);
	N = dim2 / K;
    }

    if (M * K > dim1)
	phgError(1, "%s:%d: M*K (%d) > dim1 (%d)\n", file, line, M * K, dim1);

    if (K * N > dim2)
	phgError(1, "%s:%d: K*N (%d) > dim2 (%d)\n", file, line, K * N, dim2);

    assert(p1 != NULL && p2 != NULL);
    assert(Thread(qc1).scale == Thread(qc2).scale);

    if (M == 1 && N == 1) {
	/*******************************************************************/
	/* This is the original scalar (inner product) case */
	if (res == NULL)
	    res = &res0;
	if (!const_flag1 && const_flag2) {
	    /* exchange data 1 and 2 */
	    FLOAT *p = p1;
	    p1 = p2;
	    p2 = p;
	    i = dim1;
	    dim1 = dim2;
	    dim2 = i;
	    const_flag1 = TRUE;
	    const_flag2 = FALSE;
	}
	if (!const_flag1) {
	    /* both data non constant */
	    for (; np > 0; np--, pw += winc, p1 += dim1, p2 += dim2) {
		d = 0.0;
		for (i = 0; i < K; i++)
		    d += p1[i] * p2[i];
		val += d * *pw;
	    }
	}
	else if (!const_flag2) {
	    /* data1 constant, data2 non constant */
	    if (K == 1) {
		for (i = 0; i < np; i++, pw += winc, p2 += dim2)
		    val += *p2 * *pw;
		val *= *p1;
	    }
	    else {
		for (; np > 0; np--, pw += winc, p2 += dim2) {
		    d = 0.0;
		    for (i = 0; i < K; i++)
			d += p1[i] * p2[i];
		    val += d * *pw;
		}
	    }
	}
	else {
	    /* both data constant */
	    d = 0.;
	    for (; np > 0; np--, pw += winc)
		d += *pw;
	    for (i = 0; i < K; i++)
		val += *(p1++) * *(p2++);
	    val *= d;
	}
	*res = val * Thread(qc1).scale;
	/*******************************************************************/
    }
    else {
	/*******************************************************************/
	/* compute: sum_{i=0,np} (data1[i][M][K] * data2[i][K][N]) * pw[i] */
	assert(res != NULL);
	if (!const_flag1) {
	    /* both data non constant */
	    memset(res, 0, sizeof(*res) * M * N);
	    for (; np > 0; np--, pw += winc, p1 += dim1, p2 += dim2) {
		FLOAT tmp[M * N];
		memset(tmp, 0, sizeof(tmp));
		for (i = 0; i < M; i++)
		    for (j = 0; j < N; j++)
			for (k = 0; k < K; k++)
			    tmp[i * N + j] += p1[i * K + k] * p2[k * N + j];
		for (i = 0; i < M * N; i++)
		    res[i] += tmp[i] * *pw;
	    }
	}
	else if (!const_flag2) {
	    /* data1 constant, data2 non constant */
	    FLOAT tmp[M * N];
	    memset(tmp, 0, sizeof(tmp));
	    for (; np > 0; np--, pw += winc, p2 += dim2) {
		for (j = 0; j < N; j++)
		    for (k = 0; k < K; k++)
			tmp[k * N + j] += p2[k * N + j] * *pw;
	    }
	    memset(res, 0, sizeof(*res) * M * N);
	    for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
		    for (k = 0; k < K; k++)
			res[i * N + j] += p1[i * K + k] * tmp[k * N + j];
	}
	else {
	    /* both data constant */
	    d = 0.;
	    for (; np > 0; np--, pw += winc)
		d += *pw;
	    memset(res, 0, sizeof(*res) * M * N);
	    for (i = 0; i < M; i++)
		for (j = 0; j < N; j++)
		    for (k = 0; k < K; k++)
			res[i * N + j] += p1[i * K + k] * p2[k * N + j];
	    for (i = 0; i < M * N; i++)
		res[i] *= d;
	}
	/*******************************************************************/
	if (Thread(qc1).scale != 1.0) {
	    for (i = 0; i < M * N; i++)
		res[i] *= Thread(qc1).scale;
	}
    }

    if (owner1)
	phgFree(data1);
    if (owner2)
	phgFree(data2);

    return res;
}

void
phgQCFree(QCACHE **qc)
{
    int i;

    if (*qc == NULL)
	return;

    /* free cached data */
    qc_reset_data(*qc, -1, -1);
    for (i = 0; i < phgMaxThreads; i++) {
	phgFree((*qc)->cache[i].func);
	phgFree((*qc)->cache[i].nv);
    }
    phgFree((*qc)->cache);

    if ((*qc)->f != NULL) {
	/* free qc->f[i].f_data and qc->f */
	for (i = 0; i < (*qc)->n; i++)
	    phgFree((*qc)->f[i].f_data);
	phgFree((*qc)->f);
    }

    if ((*qc)->owns_user_ctx && (*qc)->user_ctx != NULL)
	phgFree((*qc)->user_ctx);

    phgFree(*qc);
    *qc = NULL;
}

int
phgQCGetNP(QCACHE *qc)
/* returns the number of points of the current quadrature rule in "qc". */
{
    return qc == NULL ? 0 : Thread(qc).np;
}

const FLOAT *
phgQCGetNormal(QCACHE *qc, int iq, int *inc)
/* returns a pointer to the unit normal vectors at the "iq"-th quadrature point,
 * or NULL if not available,
 *
 * If "inc" != NULL, "*inc" is set to the increment of normal vectors
 * (0 for constant normal vectors).
 */
{
    int dim = 0;
    FLOAT *nv;

    if ((nv = Thread(qc).nv) == NULL && Thread(qc).rule != NULL)
	phgQuadInterfaceRuleInfo(Thread(qc).rule, &dim, NULL, &nv);

    if (inc != NULL)
	*inc = dim;

    return nv + iq * dim;
}

const FLOAT *
phgQCGetPoint(QCACHE *qc, int iq, int *inc, BOOLEAN *ref_flag)
/* returns a pointer to the coordinates of the "iq"-th quadrature point.
 *
 * If "inc" != NULL, "*inc" is set to the increment of the coordinates.
 *
 * If "ref_flag" != NULL, "*ref_flag" indicates whether the returned coordinates
 * are reference coordinates (TRUE) or physical coordinates (FALSE).
 */
{
    int dim = Thread(qc).ref_flag ? qc->qd->ref_dim : qc->qd->dim;
    FLOAT *pw = NULL;

    if (Thread(qc).rule != NULL)
	phgQuadInterfaceRuleInfo(Thread(qc).rule, NULL, &pw, NULL);

    if (inc != NULL)
	*inc = dim + 1;

    if (ref_flag != NULL)
	*ref_flag = Thread(qc).ref_flag;

    return pw + iq * (dim + 1);
}
