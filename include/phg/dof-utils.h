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

/* $Id: dof-utils.h,v 1.6 2022/09/21 05:40:59 zlb Exp $ */

#ifndef PHG_DOF_UTILS_H

/* function prototypes */

#ifdef __cplusplus
extern "C" {
#endif

DOF_TYPE *phgDofCreateSimpleType(const char *name,
                        int np_vert, int np_edge, int np_face, int np_elem);

FLOAT *phgDofEval_(DOF *dof, ELEMENT *e, const FLOAT lambda[],
			const FLOAT *basvalues, int bas_stride,
			FLOAT *values, FLOAT **data);
#define phgDofEval(dof, e, lambda, values) \
	phgDofEval_(dof, e, lambda, NULL, 0, values, NULL)

FLOAT *phgDofEvalGradient_(DOF *dof, ELEMENT *e, const FLOAT lambda[],
			const FLOAT *gradbas, FLOAT *values, FLOAT **data);
#define phgDofEvalGradient(dof, e, lambda, gradbas, values) \
	phgDofEvalGradient_(dof, e, lambda, gradbas, values, NULL)
DOF *phgDofGradient_(DOF *src, DOF **newdof, DOF_TYPE *newtype,
			const char *name, const char *srcfile, int srcline);
#define phgDofGradient(src, dest, newtype, name) \
	phgDofGradient_(src, dest, newtype, name, __FILE__, __LINE__)

FLOAT *phgDofEvalDivergence(DOF *dof, ELEMENT *e, const FLOAT lambda[],
			const FLOAT *gradbas, FLOAT *values);
#define phgDofDivergence(src, dest, newtype, name) \
	phgDofDivergence_(src, dest, newtype, name, __FILE__, __LINE__)
DOF *phgDofDivergence_(DOF *src, DOF **newdof, DOF_TYPE *newtype,
			const char *name, const char *srcfile, int srcline);
FLOAT *phgDofEvalCurl(DOF *dof, ELEMENT *e, const FLOAT lambda[],
			const FLOAT *gradbas, FLOAT *values);
#define phgDofCurl(src, dest, newtype, name) \
	phgDofCurl_(src, dest, newtype, name, __FILE__, __LINE__)
DOF *phgDofCurl_(DOF *src, DOF **newdof, DOF_TYPE *newtype,
			const char *name, const char *srcfile, int srcline);

void phgDofEvalBasisHessian(DOF *dof, ELEMENT *e, const FLOAT *lambda,
			const char *which, FLOAT *values);

/* Utility functions */

/* Note: The use of 'dof_index' instead of 'bas_no' in phgDofDirichletBC__ is
 * for supporting multiple DB masks through the array DOF->DB_masks[] */
#define phgDofDirichletBC(u, e, bas_no, func, mat, rhs, proj)	\
    phgDofDirichletBC_(u, e, bas_no * u->dim, func, mat, rhs, proj)
#define phgDofDirichletBC_(u, e, dof_index, func, mat, rhs, proj)	\
    phgDofDirichletBC__(u, e, dof_index, func, NULL, NULL, mat, rhs, proj)
BOOLEAN phgDofDirichletBC__(DOF *u, ELEMENT *e, int dof_index,
	DOF_USER_FUNC func, DOF_USER_FUNC_LAMBDA func_lambda, DOF *func_dof,
	FLOAT mat[], FLOAT rhs[], DOF_PROJ proj);

/* BLAS-like functions */
#define phgDofCopy(src, dest, newtype, name) \
	phgDofCopy_(src, dest, newtype, name, __FILE__, __LINE__)
DOF *phgDofCopy_(DOF *src, DOF **dest, DOF_TYPE *newtype,
		 const char *name, const char *srcfile, int srcline);

FLOAT phgDofDotL2Vec(DOF *x, DOF *y);
#define phgDofNormL2Vec(x)	Sqrt(phgDofDotL2Vec(x, x))
FLOAT phgDofNormLpVec(DOF *x, FLOAT p);
FLOAT phgDofNormL1Vec(DOF *x);
FLOAT phgDofNormInftyVec(DOF *x);
FLOAT phgDofMinValVec(DOF *x);
FLOAT phgDofMaxValVec(DOF *x);

FLOAT phgDofNormL1_(DOF *u, int quad_order);
FLOAT phgDofNormL2_(DOF *u, int quad_order);
FLOAT phgDofDotL2_(DOF *u, DOF *v, int quad_order);
FLOAT phgDofNormH1_(DOF *u, int quad_order);

#define phgDofNormL1(u)		phgDofNormL1_(u, -1)
#define phgDofNormL2(u)		phgDofNormL2_(u, -1)
#define phgDofDotL2(u, v)	phgDofDotL2_(u, v, -1)
#define phgDofNormH1(u)		phgDofNormH1_(u, -1)

void phgDofRandomize(DOF *u, long int seed);

#define phgDofMM(transa, transb, M, N, K, alpha, A, blka, B, beta, C) \
	phgDofMM_(transa, transb, M, N, K, alpha, A, blka, B, beta, C, \
		  __FILE__, __LINE__, TRUE)
DOF *phgDofMM_(MAT_OP transa, MAT_OP transb, int M, int N, int K,
        FLOAT alpha, DOF *A, int blka, DOF *B, FLOAT beta, DOF **Cptr,
        const char *srcfile, int srcline, BOOLEAN check);

DOF *phgDofAXPBY_(FLOAT a, DOF *x, FLOAT b, DOF **y,
	const char *srcfile, int srcline, BOOLEAN check);
#define phgDofAXPBY(a, x, b, y) \
	phgDofAXPBY_(a, x, b, y, __FILE__, __LINE__, TRUE)
#define phgDofAXPY(a, x, y)	phgDofAXPBY(a, x, (FLOAT)1.L, y)

DOF *phgDofMatVec_(FLOAT alpha, DOF *A, DOF *x, FLOAT beta, DOF **y_ptr,
		const char *srcfile, int srcline, BOOLEAN check);
#define phgDofMatVec(alpha, A, x, beta, y) \
	phgDofMatVec_(alpha, A, x, beta, y, __FILE__, __LINE__, TRUE)

DOF *phgDofAFXPBY_(FLOAT a, void (*f)(FLOAT *in, FLOAT *out), DOF *x, FLOAT b,
		   DOF **yptr, const char *srcfile, int srcline);
#define phgDofAFXPBY(a, f, x, b, yptr) \
	phgDofAFXPBY_(a, f, x, b, yptr, __FILE__, __LINE__)

DOF *phgDofAFXPBY1_(FLOAT a, FLOAT (*f)(FLOAT xvalue), DOF *x,
		    FLOAT b, DOF **yptr, const char *srcfile, int srcline);
#define phgDofAFXPBY1(a, f, x, b, yptr) \
	phgDofAFXPBY1_(a, f, x, b, yptr, __FILE__, __LINE__)

int phgDofEigenSolve_(MAT *A, MAT *B, int n, int which, FLOAT tau, int *nit,
		      FLOAT *evals, void *ssd, struct MAP_ *map, DOF **u, ...);
int phgDofEigenSolve(MAT *A, MAT *B, int n, int which, FLOAT tau, int *nit,
		     FLOAT *evals, struct MAP_ *map, DOF **u, ...);

DOF *phgDofGetSameOrderDG_(DOF *dof, int dim, const char *name,
			const char *file, int line);
#define phgDofGetSameOrderDG(dof, dim, name) \
	phgDofGetSameOrderDG_(dof, dim, name, __FILE__, __LINE__)

/*void phgDofMapDataInElement(DOF_TYPE *type, const ELEMENT *e, int dim,
			  FLOAT *values);*/
MAT *phgDofGetTransferMatrix(DOF *u, DOF *v, MAP *map_u, MAP *map_v);

#ifdef __cplusplus
}
#endif

#define PHG_DOF_UTILS_H
#endif
