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

/* utility functions for eXtended FEM
 * $Id: xfem-utils.c,v 1.16 2022/08/12 06:42:25 zlb Exp $ */

#include "phg.h"

#include <math.h>
#include <strings.h>	/* bzero() */
#include <string.h>	/* memcpy() */

FLOAT
phgXFEMDot_(XFEM_INFO *xi, DOF *um, DOF *up, DOF *vm, DOF *vp,
	    FLOAT **e_sum, int added_order)
/* computes L2 dot-product (u,v), where
 *	u:=um and v:=vm in \Omega^-,
 * and
 *	u:=up and v:=vp in \Omega^+.
 *
 * See phgXFEM_FaceJump2_ for the arguments "e_sum" and "added_order".
 */
{
    FLOAT dot = 0.0, val, *e_data = NULL;
    GRID *g = um->g;
    ELEMENT *e;
    QCACHE *qc_m, *qc_p = NULL;
    DOF *ls = NULL, *ls_grad = NULL;
    int Q_um, Q_vm, Q_up = 0, Q_vp = 0;

    assert(xi != NULL || up == NULL || up == um);
    assert(xi != NULL || vp == NULL || vp == vm);

    qc_m = phgQCNew(QD_DEFAULT, um);
    Q_um = phgQCAddFEFunction(qc_m, um);
    Q_vm = um == vm ? Q_um : phgQCAddFEFunction(qc_m, vm);

    if (xi != NULL) {
	ls = xi->ls;
	ls_grad = xi->ls_grad;
	qc_p = phgQCNew(QD_DEFAULT, up);
	Q_up = phgQCAddFEFunction(qc_p, up);
	Q_vp = up == vp ? Q_up : phgQCAddFEFunction(qc_p, vp);
    }

    assert(DofDim(um) == DofDim(vm));
    assert(xi == NULL || DofDim(um) == DofDim(up));
    assert(xi == NULL || DofDim(um) == DofDim(vp));

    if (e_sum != NULL)
	*e_sum = e_data = phgCalloc(g->nelem, sizeof(*e_data));

#if USE_OMP
#pragma omp parallel for private(e) reduction(+:dot)
#endif	/* USE_OMP */
    ForAllElementsBegin(g, e) {
	FLOAT *rule_m = NULL, *rule_p = NULL;
	FLOAT corners[NVert][Dim];
	int i, j = 0, order;	/* quadrature order */
	int mark = xi == NULL ? -1 : xi->info[e->index].mark;

	i = DofTypeOrder(um, e);
	if (xi != NULL && i < DofTypeOrder(up, e))
	    i = DofTypeOrder(up, e);
	j = DofTypeOrder(vm, e);
	if (xi != NULL && j < DofTypeOrder(vp, e))
	    j = DofTypeOrder(vp, e);
	order = i + j + (mark == 0 ? added_order : 0);

	if (mark == 0)
	    phgQuadInterface(ls, ls_grad, e, order,
			 mark <= 0 ? &rule_m : NULL, NULL,
			 mark >= 0 ? &rule_p : NULL);
	else
#ifdef PHG_TO_P4EST
	    phgQuadInterfaceCuboid
#else	/* PHG_TO_P4EST */
	    phgQuadInterfaceTetrahedron
#endif	/* PHG_TO_P4EST */
			(NULL, 0, NULL, phgGeomGetCorners(g, e, corners), order,
			 mark <= 0 ? &rule_m : NULL, NULL,
			 mark >= 0 ? &rule_p : NULL);

	phgQCSetRule(qc_m, rule_m, -1.);
	val = phgQCIntegrate(qc_m, e->index, Q_um, 0, qc_m, e->index, Q_vm, 0);
	dot += val;
	if (e_data != NULL)
	    e_data[e->index] += val;
	phgFree(rule_m);

	if (xi == NULL)
	    continue;
	phgQCSetRule(qc_p, rule_p, -1.);
	val = phgQCIntegrate(qc_p, e->index, Q_up, 0, qc_p, e->index, Q_vp, 0);
	dot += val;
	if (e_data != NULL)
	    e_data[e->index] += val;
	phgFree(rule_p);
    } ForAllElementsEnd

    phgQCFree(&qc_m);
    if (xi != NULL)
	phgQCFree(&qc_p);

#if USE_MPI
    if (g->nprocs > 1) {
	FLOAT d = dot;
	MPI_Allreduce(&d, &dot, 1, PHG_MPI_FLOAT, PHG_SUM, g->comm);
    }
#endif	/* USE_MPI */

    return dot;
}

/* Information for jmp_func. */
static struct {
    DOF *u, *u1, *intf;
    ELEMENT *e, *e1;
} jmp_info;
#if USE_OMP
# pragma omp threadprivate(jmp_info)
#endif	/* USE_OMP */

static void
jmp_func(DOF *unused_dof, ELEMENT *e, int bno, const FLOAT *xref,
	 FLOAT *values)
/* callback function for computing jumps of functions */
{
    int i, dim = DofDim(jmp_info.u1);
    FLOAT values1[dim];

    assert(e == jmp_info.e);
    phgDofEval(jmp_info.u, e, xref, values);
    if (jmp_info.intf != NULL) {
	phgDofEval(jmp_info.intf, e, xref, values1);
	for (i = 0; i < dim; i++)
	    values[i] -= values1[i];
    }
    if (jmp_info.e1 != NULL && jmp_info.u1 != NULL) {
	/* interior face or interface, compute jump */
	FLOAT xref1[Dim + 1];
	assert(dim == DofDim(jmp_info.u));
	if (jmp_info.e1 != e) {
	    FLOAT x, y, z;
	    phgGeomLambda2XYZ(jmp_info.u1->g, jmp_info.e, xref, &x, &y, &z);
	    phgGeomXYZ2Lambda(jmp_info.u1->g, jmp_info.e1, x, y, z, xref1);
	    xref = xref1;
	}
	phgDofEval(jmp_info.u1, jmp_info.e1, xref, values1);
	for (i = 0; i < dim; i++)
	    values[i] -= values1[i];
    }
}

FLOAT
phgXFEMFaceJump2_(XFEM_INFO *xi, DOF *um, DOF *bm, DOF *up, DOF *bp,
		  DOF *intf, PROJ proj, FLOAT h_power, int f_mask,
		  FLOAT **e_sum, int added_order)
/* This function computes the following integrals:
 *	\sum_{F \in subset(f_mask)} (
 *		\int_{F\cap\Omega^-} proj([[um]])^2 +
 *		\int_{F\cap\Omega^+} proj([[up]])^2
 *	) * h^{h_power}
 *	+ \int_{\Gamma} proj(um - up - intf)^2 * h^{h_power}
 * where:
 *	"\Omega^-" and "\Omega^+" are the subdomains,
 *	"\Gamma" is the interface,
 *	"F" denotes a face,
 *	"proj" denotes some projection w.r.t. the normal direction,
 *	"h" denotes the diameter of the face or interface.
 *
 * Given a face F, "[[.]]" denotes the jump across F, defined as follows:
 *	- if F is the common face of element e and element e1, then
 *		[[um]] := um|_e - um|_e1, [[up]] := up|_e - up|_e1,
 *	- if F is a boundary face, then
 *		[[um]] := um - bm, [[up]] := up - bp.
 * For "bm", "bp" and "intf", NULL pointer means zero.
 *
 *
 * The parameter "h_power" specifies the power, if it's not 0 then the result
 * in each face will be multiplied by the diameter of the face to the power
 * of "h_power".
 *
 * The parameter "f_mask" defines the subset of element faces to include, e.g.:
 *	-1			- all faces,
 *	INTERIOR		- all interior faces,
 *	BDRY_MASK		- all boundary faces,
 *	DIRICHLET|NEUMANN	- Dirichlet and Neumann faces.
 *
 * The parameter "added_order" defines the added numerical quadrature order
 * for interface elements, the actual numerical quadrature order used is:
 *	2 * DofTypeOrder("um", e) in non interface elements, and
 *	2 * DofTypeOrder("um", e) + "added_order" in interface elements
 *
 * The argument "e_sum", if not NULL, will be set to point to a dynamically
 * allocated buffer of size g->nelem, with "(*e_sum)[i]" set to the sum of
 * the above integrals on all faces and interface piece of the element with
 * local index "i".
 *
 * Note: this function can also be used to compute face jumps for non-interface
 * problems with "xi" == NULL, in this case it is required that either
 * "up" == "um" or "up == NULL".
 */
{
    GRID *g = um->g;
    ELEMENT *e;
    DOF *dummy, *J = NULL;	/* for computing face or interface jump */
    VEC *vec = NULL;
    QCACHE *q;
    FLOAT val = 0., a;
    int Qj;

    assert(bm == NULL || DofDim(bm) == DofDim(um));
    assert(up == NULL || DofDim(up) == DofDim(um));
    assert(bp == NULL || DofDim(bp) == DofDim(um));
    assert(intf == NULL || DofDim(intf) == DofDim(um));
    assert(xi != NULL || up == NULL || up == um || bp == NULL || bp == bm);

    phgSetupHalo(g, HALO_FACE);

    dummy = phgDofNew(g, DOF_ANALYTIC, DofDim(um), "dummy", DofLambdaFunction);
    phgDofSetLambdaFunction(dummy, jmp_func);
    q = phgQCNew(QD_DEFAULT, /*um*/dummy);
    Qj = phgQCAddFEFunction(q, dummy);

    if (e_sum != NULL) {
	MAP *map;
	J = phgDofNew(g, DOF_DGn[0], 1, "sums", DofNoAction);
	map = phgMapCreate(J, NULL);
	vec = phgMapCreateVec(map, 1);
	phgVecDisassemble(vec);
	phgMapDestroy(&map);
    }

#ifdef Mp
# undef Mp
#endif
#ifdef PHG_TO_P4EST
#define Mp	PHG_to_p8est_map_e2g 
#else	/* PHG_TO_P4EST */
#define Mp	phgMapE2G
#endif	/* PHG_TO_P4EST */

#if USE_OMP
#pragma omp parallel for private(e) reduction(+:val)
#endif  /* USE_OMP */
    ForAllElementsBegin(g, e) {
	ELEMENT *e1;
	FLOAT *rule_m = NULL, *rule = NULL, *rule_p = NULL;
	INT eno = e->index;
	int k, face, order = DofTypeOrder(um, e) * 2;

	/*==================== the interface part ====================*/
	if (xi != NULL && xi->info[eno].mark == 0) {
	    phgQuadInterface(xi->ls, xi->ls_grad, e, order + added_order,
				NULL, &rule, NULL);
	    phgQCSetRule(q, rule, -1.);
	    jmp_info.u = um;
	    jmp_info.u1 = up;
	    jmp_info.intf = intf;
	    jmp_info.e = jmp_info.e1 = e;
	    a = phgQCIntegrateFace(q, eno, -1, Qj, proj, 0,
				   q, eno, -1, Qj, proj, 0);
	    jmp_info.intf = NULL;
	    if (h_power != 0.)
		/* FIXME: compute diameter of the interface piece, and handle
		 * merged elements */
		a *= Pow(phgGeomGetDiameter(g, e), h_power);
	    val += a;
	    if (e_sum != NULL)
		phgVecAddGlobalEntry(vec, 0, Mp(vec->map, 0, e, 0), a);
	    phgFree(rule);
	    rule = NULL;
	}

	/* the faces of the element */
	for (face = 0; face < NFace; face++) {
	    FLOAT nv[Dim];
	    e1 = phgGetNeighbour(g, e, face);
	    if (e1 != NULL) {
		/* a face is only processed by the smaller of the two
		 * neighbouring elements, and by the one with smaller
		 * global index if the two elements are of the same size,
		 * to avoid double counting or redundant computation */
		if (e->generation < e1->generation)
		    continue;
		if (e->generation == e1->generation &&
		    GlobalElement(g, e->index) > GlobalElement(g, e1->index))
		    continue; /* process each interior face just once */
	    }
	    jmp_info.e = e;
	    jmp_info.e1 = e1;
	    phgGeomGetFaceOutNormal(g, e, face, nv);
	    phgQCSetConstantNormal(q, nv);
	    /* Note: e1 is either NULL or valid (since gen(e) >= gen(e1)) */
	    if (e1 != NULL && xi != NULL &&
		xi->info[e->index].mark * xi->info[e1->index].mark < 0) {
		/* The face is flat and is part of the interface
		 * (e and e1 in different subdomains) */
		if (xi->info[e->index].mark < 0) {
		    jmp_info.u = um;
		    jmp_info.u1 = up;
		}
		else {
		    jmp_info.u = up;
		    jmp_info.u1 = um;
		}
		phgQuadInterfaceFace(NULL, xi->ls_grad, e, face,
						order, &rule, NULL, NULL);
		phgQCSetRule(q, rule, -1.);
		a = phgQCIntegrateFace(q, eno, face, Qj, proj, 0,
				       q, eno, face, Qj, proj, 0);
		if (h_power != 0.)
		    /* FIXME: handle merged elements */
		    a *= Pow(phgGeomGetFaceDiameter(g, e, face), h_power);
		val += a;
		if (e_sum != NULL) {
		    phgVecAddGlobalEntry(vec, 0, Mp(vec->map, 0, e, 0), a);
		    phgVecAddGlobalEntry(vec, 0, Mp(vec->map, 0, e1, 0), a);
		}
		phgFree(rule);
		rule = NULL;
		continue;
	    }
	    /* check against f_mask */
	    if (!(GetFaceBTYPE(g, e, face) & f_mask))
		continue;
	    if (xi == NULL || xi->info[eno].mark < 0) {
#ifndef PHG_TO_P4EST
		FLOAT corners[Dim + 1][Dim];
		phgGeomGetCorners(g, e, corners);
		phgQuadInterfaceTetrahedronFace(NULL, 0, NULL,
			corners, face, order, &rule_m, NULL, NULL);
#else
		phgQuadInterfaceCuboidFace(NULL, 0, NULL,
			e->corners, face, order, &rule_m, NULL, NULL);
#endif
	    }
	    else if (xi->info[eno].mark > 0)
		phgQuadInterfaceFace(NULL, xi->ls_grad, e, face,
				     order, NULL, NULL, &rule_p);
	    else
		phgQuadInterfaceFace(xi->ls, xi->ls_grad, e, face,
				order + added_order, &rule_m, NULL, &rule_p);
	    if (e1 == NULL)
		jmp_info.e1 = e;
	    for (k = 0; k < 2; k++) {
		phgQCSetRule(q, k == 0 ? rule_m : rule_p, -1.);
		if (phgQCGetNP(q) == 0)
		    continue;
		jmp_info.u = jmp_info.u1 = k == 0 ? um : up;
		if (e1 == NULL)
		    jmp_info.u1 = k == 0 ? bm : bp;
		a = phgQCIntegrateFace(q, eno, face, Qj, proj, 0,
				       q, eno, face, Qj, proj, 0);
		if (h_power != 0.)
		    /* FIXME: handle merged elements */
		    a *= Pow(phgGeomGetFaceDiameter(g, e, face), h_power);
		val += a;
		if (e_sum != NULL) {
		    phgVecAddGlobalEntry(vec, 0, Mp(vec->map, 0, e, 0), a);
		    if (e1 != NULL)
			phgVecAddGlobalEntry(vec, 0, Mp(vec->map, 0, e1, 0), a);
		}
	    }
	    phgFree(rule_m);
	    phgFree(rule_p);
	    rule_m = rule_p = NULL;
	}
    } ForAllElementsEnd
#undef Mp

    phgQCFree(&q);
    phgDofFree(&dummy);

    if (e_sum != NULL) {
	phgVecAssemble(vec);
#ifdef PHG_TO_P4EST
	phgFree(J->data);
	J->data = vec->data;
	vec->data = NULL;
	phgDofUpdateGhost(J);
#else	/* PHG_TO_P4EST */
	phgMapLocalDataToDof(vec->map, 1, &J, vec->data);
#endif	/* PHG_TO_P4EST */
	phgVecDestroy(&vec);
	*e_sum = J->data;
	J->data = NULL;
	phgDofFree(&J);
    }

#if USE_MPI
    if (g->nprocs > 1) {
	a = val;
	MPI_Allreduce(&a, &val, 1, PHG_MPI_FLOAT, PHG_SUM, g->comm);
    }
#endif

    return val;
}

#ifndef PHG_TO_P4EST	/* avoid doubly defined functions */

void
phgXFEMProcessEmptyRows(SOLVER *solver)
/* process empty row */
{
    int i;
    MAT *A = phgSolverGetMat(solver);
    MAT_ROW *row;

    phgMatAssemble(A);
    phgMatDisassemble(A);
    assert(A->type == PHG_UNPACKED);
    for (i = 0; i < A->rmap->nlocal; i++) {
	row = A->rows + i;
	if (row->ncols == 0) {
	    /* empty row */
	    row->ncols = row->alloc = 1;
	    phgFree(row->cols);
	    row->cols = phgAlloc(sizeof(*row->cols));
	    row->cols[0] = i;
	    phgFree(row->data);
	    row->data = phgAlloc(sizeof(*row->data));
	    row->data[0] = 1.0;
	}
    }
    phgMatAssemble(A);
    phgMatDestroy(&A);
}

FLOAT *
phgOrthoMatrixGen(int n, int m, FLOAT *a)
/* This function is used by xfem.c and dg.c for computing orthonormaled basis
 * functions. It performs QR decomposition of the 'm'x'n' matrix 'a',
 * returns the pointer to a dynamically allocated buffer of 'n'*('n'+1)/2
 * FLOATs containing the triangular matrix P = R^(-1).
 *
 * It is required that 'a' is full rank and 'n' <= 'm'.
 *
 * See the macros "mA" and "mP" for ordering of matrix elements in the buffers.
 *
 * Here, 'n' is the number of basis functions in an element,
 * 'm' is the number of quadrature points used to compute the mass matrix,
 * and 'a'(i,j) = \phi_j(x_i), 0 <= i < m, 0 <= j < n, where \phi_j is the j-th
 * basis function and x_i is the i-th quadrature point.
 *
 * The elements of 'a', 'R' and 'P' are in column-major (i.e., FORTRAN) order.
 */
{
    int i, j, k;
    FLOAT *p = phgAlloc(n * (n + 1) / 2 * sizeof(*p));	/* allocate matrix P */
    FLOAT d;

    assert(n <= m);

#define mP(i,j)	p[(j) * ((j) + 1) / 2 + (i)]
#define mA(i,j)	a[(j) * m + (i)]

#if 0
#warning
/* print matrices for debugging: A = QR */
printf("A = [");
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
	printf("%0.15g%s", (double)mA(i,j), j == n - 1 ? "" : ",");
    printf(i == m - 1 ? "];\n" : "; ");
}
#endif

    /* compute QR decomposition of A */

#if USE_LAPACK
    if (sizeof(FLOAT) == sizeof(double)) {
	/* use LAPACK DGEQRF */
	int lwork, info;
	FLOAT *tau = NULL, *work;
	/* Note: may use DGEQP3 (with pivoting) for better precision */
	extern void F77_FUNC(dgeqrf,DGEQRF)(int *m, int *n, FLOAT *a, int *lda,
			     FLOAT *tau, FLOAT *work, int *lwork, int *info);

	/* get lwork size */
	lwork = -1;
	work = &d;
	F77_FUNC(dgeqrf,DGEQRF)(&m, &n, a, &m, tau, work, &lwork, &info);
	assert(info == 0);
	lwork = (int)(d + 0.5);
	work = phgAlloc((lwork + n) * sizeof(*work));
	tau = work + lwork;
	F77_FUNC(dgeqrf,DGEQRF)(&m, &n, a, &m, tau, work, &lwork, &info);
	phgFree(work);
	assert(info == 0);
    }
    else
#endif	/* USE_LAPACK */
    {
	/* use modified Gram-Schmidt scheme */
	for (i = 0; i < n; i++) {
	    /* normalize a(:,i) */
	    d = 0.0;
	    for (k = 0; k < m; k++)	
		d += mA(k,i) * mA(k,i);
	    assert(d != 0.0);
	    d = Sqrt(d);
	    mP(i,i) = d;
	    d = 1.0 / d;
	    for (k = 0; k < m; k++)
		mA(k,i) *= d;
	    for (j = i + 1; j < n; j++) {
		d = 0.0;
		for (k = 0; k < m; k++)
		    d += mA(k,i) * mA(k,j);
		mP(i,j) = d;
		for (k = 0; k < m; k++)
		    mA(k,j) -= d * mA(k,i);
	    }
	}
#if 0
#warning
printf("Q = [");
for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++)
	printf("%0.15g%s", (double)mA(i,j), j == n - 1 ? "" : ",");
    printf(i == m - 1 ? "];\n" : "; ");
}
#endif
	/* copy R to the upper-triangular part of A */
	for (i = 0; i < n; i++)
	    for (j = i; j < n; j++)
		mA(i,j) = mP(i,j);
    }

    /* Compute "p" := R^(-1), R is stored in the upper triangular part of mA */
    /* first, set P := I */
    for (i = 0; i < n; i++)
	for (j = i; j < n; j++)
	    mP(i,j) = i == j ? 1. : 0.;
    /* then, compute P = R^(-1) */
    for (i = n - 1; i >= 0; i--) {
	d = 1.0 / mA(i, i);	/* diagonal entry */
	for (j = i; j < n; j++) {
	    mP(i,j) *= d;
	    for (k = 0; k < i; k++)
		mP(k,j) -= mA(k,i) * mP(i,j);
	}
    }

#if 0
#warning
printf("P = [");
for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
	printf("%0.15g%s", i > j ? 0. : (double)mP(i,j), j == n - 1 ? "" : ",");
    printf(i == n - 1 ? "];\n" : "; ");
}
#endif

#undef mA
#undef mP

    return p;
}

void
phgOrthoMatrixApply(int n, int bdim, int stride, int dim, FLOAT *buffer,
		    const FLOAT *P)
/* This function is used by xfem.c and dg.c for computing orthonormaled basis
 * functions. It applies the orthonormalization matrix P to values of all basis
 * functions (or their gradients) at a given point.
 *
 *	n	- number of basis functions,
 *	bdim	- number of components in the basis functions (1 for scalar
 *		  basis, Dim for vector basis, etc.), 
 *	dim	- 1 for basis functions, Dim for gradient(basis functions),
 *	stride	- increment in "buffer" between basis functions,
 *	buffer	- contains values for the original basis functions (input), and
 *		  for the orthonormalized basis functions (output), at a given
 *		  point, stored as buffer[n][stride]
 */
{
    int i, j, k, l;
    FLOAT tmp[n][bdim][bdim * dim];

    if (P == NULL)
	return;

    dim *= bdim;

    for (i = 0; i < n; i++)
	for (l = 0; l < bdim; l++)
	    for (k = 0; k < dim; k++) {
		tmp[i][l][k] = buffer[(i * bdim + l) * stride + k];
		buffer[(i * bdim + l) * stride + k] = 0.;
	    }

    for (i = 0; i < n; i++)
	for (j = 0; j <= i; j++, P++)
	    for (l = 0; l < bdim; l++)
		for (k = 0; k < dim; k++)
		    buffer[(i * bdim + l) * stride + k] += *P * tmp[j][l][k];
}

#endif	/* !defined(PHG_TO_P4EST) */
