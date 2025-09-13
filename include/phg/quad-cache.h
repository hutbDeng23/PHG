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

/* $Id: quad-cache.h,v 1.41 2022/05/13 02:02:17 zlb Exp $ */

#ifndef PHG_QUAD_CACHE_H

#ifdef __cplusplus
extern "C" {
#endif

/* FID is used to identify functions attached to the QCACHE struct,
 *	< 0	basis functions and derivatives
 *	>= 0	user functions
 */
enum {
    Q_FUNC	=  0,	/* user function */
    Q_BAS	= -1,	/* basis functions */
    Q_GRAD	= -2,	/* gradient of basis functions */
    Q_NONE	= -3	/* used by base_fid, the smallest FID value */
#define FID_min	(Q_NONE + 1)
};

extern int Q_DIV, Q_CURL;	/* define in qcache.c */

/* OP_FUNC defines the prototype for the operator function in
 * phgAddCallBackCoefficient().
 *
 * The function performs the "coeff" op "in" -> "out" operation at a given
 * quadrature point.
 * 	"cdim" and "dim" specify resp. the dimension of "coeff" and "data"
 * 		("coeff"  may be NULL, meaning no coefficient)
 *
 * The other arguments provides optional information which may be needed
 * by the function:
 *     "ctx"	points to user-defined optional data,
 *     "qc"	points to the QCACHE object,
 *     "x"	points to the coordinates of the current point,
 *     "nv"	points to the unit normal vector of the interface at "x".
 *
 * The return value is the dimension of the transformed data. If "data" == NULL
 * then the function simply returns the dimension of "out".
 *
 * As examples, see Op_div() and Op_curl() below.
 */
typedef struct QCACHE_ QCACHE;
typedef int (*OP_FUNC) (int dim1, const FLOAT *data1,
			int dim2, const FLOAT *data2, FLOAT *out,
			QCACHE *qc, int iq, const void *ctx);

#define NEC	2	/* # element/face pairs simultaneously cached */

/* QCACHE stores values of basis functions, gradient of basis functions,
 * and user functions at the quadrature points, for use by phgQCIntegrate
 *
 * The "owner" members indicates whether the data is owned and should be freed
 * by this QCACHE object. */
struct QCACHE_ {
    /* pointer to the QDESC object (the UI) */
    struct QDESC_ *qd;
    /* pointer to the finite element space */
    void *fe;	/* for PHG it's a (DOF *) pointer */

    /* list of added functions, dynamically allocated */
    struct {
	/* pointers for the user function, at most one of them !=NULL */
	void	*f_fe;		/* FE function */
	FLOAT	*f_data;	/* data for constant user function */
	FUNC_3D f_xyz;		/* xyz user function */

	OP_FUNC op;		/* operator (base_fid!=Q_NONE) */
	const void *ctx;	/* user context */

	int	f_fid;		/* a pre-registered function */
	int	f_dim;		/* data dimension for f_{fe,data,xyz} */

	int	base_fid;	/* If base_fid==Q_NONE, then the function
				 * values are stored. Otherwise the values
				 * are regarded as the coefficient of the
				 * function base_fid, and the resulting function
				 * values are given by:
				 * 	op(func, base_fid). */
	BOOLEAN	const_flag;
    } *f;

    /* the following struct stores information about actual cached data.
     * For OpenMP support, "cache" is an allocated array of phgMaxThreads
     * entries, with "cache[t]" for thread "t" */
    struct {
	/* the quadrature rule, exactly one of them is not NULL */
	QUAD *quad;
	FLOAT *rule;
	BOOLEAN ref_flag;	/* whether "rule" uses reference coordinates */

	/* scaling factor (depends on rule or quad) */
	FLOAT scale;

	FLOAT *nv;

	INT eno[NEC];	/* numbers of cached elements */
	/* nbr quad. points, current cache index, face no., nbr bases */
	int np, index0, face[NEC], nb[NEC];

	/* cache[i] stores cachaed data for thread i */
	struct {
	    FLOAT *data[NEC];	/* cached data */
	    FLOAT *dot_n[NEC];	/* cached projected data */
	    FLOAT *cross_n[NEC];/* cached projected data */
	    int	dim;		/* data dimension */
	    int nb[NEC];	/* nb */
	    BOOLEAN owner[NEC];	/* whether data owner by this object */
	} bas, grad, *func, *me; /* Note: me is a work variable */
    } *cache;

    /* If "bas_hook" != NULL, it will be called inplace of qd->bas() */
    void (*bas_hook)(struct QCACHE_ *qc, INT eno, const FLOAT *xref,
		     FLOAT *buffer, int stride);
    /* If "grd_hook" != NULL, it will be called inplace of qd->grd() */
    void (*grd_hook)(struct QCACHE_ *qc, INT eno, const FLOAT *xref,
		     FLOAT *buffer, int stride);

    /* "user_ctx" is reserved for use by user program */
    void *user_ctx;
    BOOLEAN owns_user_ctx;	/* TRUE ==> phgQCFree() should free user_ctx */

    int	n, alloc;		/* # functions and allocated size of f[] */
    BOOLEAN in_use;		/* preventing adding functions when in use */
};

/* QDESC is used to interact with specific finite element packages.
 * Please see the implementation of QD_PHG for an example. */
typedef struct QDESC_ {
    const char *name;	/* interface name */

    /* Note: any of the following function pointers may be set to NULL
     * if not needed. */

    /* to_ref() transforms coordinates "x" to reference coordinates "xref" */
    void (*to_ref)(QCACHE *qc, INT eno, FLOAT x[], FLOAT xref[]);

    /* from_ref() transforms reference coordinates "xref" to coordinates "x" */
    void (*from_ref)(QCACHE *qc, INT eno, FLOAT xref[], FLOAT x[]);

    /* b_dim() returns the space dimension of the basis functions (1 or 3) */
    int (*b_dim)(QCACHE *qc);

    /* f_dim() returns the space dimension of the FE function "fe" */
    int (*f_dim)(QCACHE *qc, void *fe);

    /* n_bas() returns number of basis functions in the given element */
    int (*n_bas)(QCACHE *qc, INT eno);

    /* bas() computes the values of all basis functions at the given point
     * (specified by "xref", in reference coordinates). The computed values
     * are stored as FLOAT buffer[nb][stride][dim], where:
     * 		nb is the number of basis functions in the element,
     * 		dim is the space dimension of the basis functions,
     * 		stride is an argument provided by the calling function. */
    void (*bas)(QCACHE *qc, INT eno, const FLOAT *xref, FLOAT *buffer,
		int stride);

    /* grad() computes the gradients of all basis functions at the given point
     * (specified by "xref", in reference coordinates).
     * The computed values are stored as FLOAT buffer[nb][stride][dim][Dim],
     * where Dim==3 for 3D, and nb and dim are the same as in bas(). */
    void (*grd)(QCACHE *qc, INT eno, const FLOAT *xref, FLOAT *buffer,
		int stride);

    /* fun() computes the values of the FE function "fe" at the given point
     * (specified by "xref", in reference coordinates).
     * The computed values are stored as FLOAT buffer[dim], where dim is given
     * by the "dim" argument of the function phgQCAddFEFunction(). */
    void (*fun)(QCACHE *qc, void *fe, INT eno, const FLOAT *xref,
		FLOAT *buffer);

    /* Note:
     *	 1. Only for QD_PHG "ref_dim" == "dim" + 1 (barycentric coordinates),
     *	    otherwise we always have "ref_dim" == "dim".
     *   2. If reference coordinates are stored in a rule, only the first
     *	    "dim" coordinates of each point are stored. */
    int dim, ref_dim;	/* space dimension, reference coordinates dimension */
} QDESC;

extern QDESC QD_PHG_;
#define QD_PHG		(&QD_PHG_)
#define QD_DEFAULT	QD_PHG

#define phgQCAddFECoefficient(qc, f_fe, base_fid) \
	phgQCAddFunction_(qc, f_fe, NULL, NULL, Q_NONE, 0, \
			  NULL, base_fid, NULL)
#define phgQCAddConstantCoefficient(qc, f_data, dim, base_fid) \
	phgQCAddFunction_(qc, NULL, f_data, NULL, Q_NONE, dim, \
			  NULL, base_fid, NULL)
#define phgQCAddXYZCoefficient(qc, f_func, dim, base_fid) \
	phgQCAddFunction_(qc, NULL, NULL, f_func, Q_NONE, dim, \
			  NULL, base_fid, NULL)
#define phgQCAddFIDCoefficient(qc, f_fid, base_fid) \
	phgQCAddFunction_(qc, NULL, NULL, NULL, f_fid, 0, \
			  NULL, base_fid, NULL)

#define phgQCAddOperator(qc, op, base_fid) \
	phgQCAddFunction_(qc, NULL, NULL, NULL, Q_NONE, 0, op, base_fid, NULL)

#define phgQCAddFEFunction(qc, fe) \
	phgQCAddFECoefficient(qc, fe, Q_NONE)
#define phgQCAddConstantFunction(qc, data, dim) \
	phgQCAddConstantCoefficient(qc, data, dim, Q_NONE)
#define phgQCAddXYZFunction(qc, func, dim) \
	phgQCAddXYZCoefficient(qc, func, dim, Q_NONE)
int phgQCAddFunction_(QCACHE *qc,
		void *f_fe, FLOAT *f_data, FUNC_3D f_xyz, int f_fid, int dim,
		OP_FUNC op, int base_fid, const void *ctx);

int phgQCAddProjection(QCACHE *qc, PROJ proj, int fid);

#define phgQCSetRule(qc, rule, scale) \
	phgQCSetRule_(qc, scale, rule, NULL)
#define phgQCSetQuad(qc, quad, scale) \
	phgQCSetRule_(qc, scale, NULL, quad)
void phgQCSetRule_(QCACHE *qc, FLOAT scale, FLOAT *rule, QUAD *quad);

void phgQCSetConstantNormal(QCACHE *qc, const FLOAT *nv);
QCACHE *phgQCNew(QDESC *qd, void *fe);
QCACHE *phgQCClone(QCACHE *qc0);
void phgQCFree(QCACHE **qc);

int phgQCGetNP(QCACHE *qc);
const FLOAT *phgQCGetNormal(QCACHE *qc, int iq, int *inc);
const FLOAT *phgQCGetPoint(QCACHE *qc, int iq, int *inc, BOOLEAN *ref_flag);

FLOAT *phgQCIntegrate_(const char *file, const int line,
	QCACHE *qc1, INT eno1, int face1, int fid1, PROJ proj1, int i1,
	QCACHE *qc2, INT eno2, int face2, int fid2, PROJ proj2, int i2,
	int M, int N, int K, FLOAT *res);

#define phgQCIntegrateFaceM(qc1, e1, face1, fid1, proj1, i1, \
			    qc2, e2, face2, fid2, proj2, i2, M, N, K, res) \
	phgQCIntegrate_(__FILE__, __LINE__, \
			qc1, e1, face1, fid1, proj1, i1, \
			qc2, e2, face2, fid2, proj2, i2, M, N, K, res)

#define phgQCIntegrateM(qc1, e1, fid1, i1, qc2, e2, fid2, i2, M, N, K, res) \
	phgQCIntegrateFaceM(qc1, e1, -1, fid1, PROJ_NONE, i1, \
			    qc2, e2, -1, fid2, PROJ_NONE, i2, M, N, K, res) \

#ifndef __INTEL_COMPILER
# define phgQCIntegrateFace(qc1, e1, face1, fid1, proj1, i1, \
			   qc2, e2, face2, fid2, proj2, i2) \
	(*phgQCIntegrateFaceM(qc1, e1, face1, fid1, proj1, i1, \
			      qc2, e2, face2, fid2, proj2, i2, 1, 1, 0, NULL))
#else	/* !defined(__INTEL_COMPILER) */
# warning enable ICC-2021.1 workaround on LSSC4
/*  2021.05.05: this is a workaround for icc version 2021.1 (gcc 4.8.5 compat.)
 *		on LSSC4, without this, "p4est/interface -refine k" produces
 *		wrong results for k>0. The bug seems to be related to the
 *			"phgQCIntegrateFace(...) - phgQCIntegrateFace(...)"
 *		expressions in build_linear_system() of interface-p4est.c */
static FLOAT
phgQCIntegrateFace(
	QCACHE *qc1, INT eno1, int face1, int fid1, PROJ proj1, int i1,
	QCACHE *qc2, INT eno2, int face2, int fid2, PROJ proj2, int i2)
{
    return *phgQCIntegrate_(__FILE__, __LINE__,
	qc1, eno1, face1, fid1, proj1, i1,
	qc2, eno2, face2, fid2, proj2, i2, 1, 1, 0, NULL);
}
#endif	/* !defined(__INTEL_COMPILER) */

#define phgQCIntegrate(qc1, e1, t1, i1, qc2, e2, t2, i2) \
	phgQCIntegrateFace(qc1, e1, -1, t1, PROJ_NONE, i1, \
			   qc2, e2, -1, t2, PROJ_NONE, i2)

#ifdef __cplusplus
}
#endif

#define PHG_QUAD_CACHE_H
#endif
