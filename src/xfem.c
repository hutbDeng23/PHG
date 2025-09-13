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
 * $Id: xfem.c,v 1.154 2022/09/23 02:47:13 zlb Exp $ */

#include "phg.h"

#include <math.h>
#include <strings.h>	/* bzero() */
#include <string.h>	/* memcpy() */

/* parameters with cmdline opptions */
typedef struct {
/* 
 * "ctol" is the tolerance on the local curvature: h * kapa <= ctol,
 * where kappa is the max absolute curvature of the interface in the element,
 * and h the diameter of the element. If ctol <= 0 then the curvatures are
 * not checked
 *
 * "dtol" is the interface deviation tolerance: interface deviations <= dtol.
 * If "dtol" <= 0.0 then interface deviations are not checked.
 *
 * "etol" is the edge fraction on a cut edge below which the element
 * is considered "small". If "etol" <= 0 then edge cuts are not checked.
 *
 * "vtol" is the volume fraction in a cut element below which the element
 * is considered "small". If "vtol" <= 0 then volumes are not checked.
 *
 * "horder" is the polynomial order used to approximate the Hessian of the
 * levelset function when estimating curvatures of the interface (only if the
 * levelset function is non-polynomial, otherwise it's unused).
 *
 * If "dbg_vtk" is not FALSE, then each call to phgXFEMInit outputs a file
 * "merge.vtk" which can be used to debug or visualize the merged elements.
 *
 * "ortho" selects the set of elements in which finite element basis functions
 * are orthogonalized (by multiplying the original basis functions functions
 * with an upper triangular matrix in each element belonging to the set.
 */
    FLOAT ctol, dtol, etol, vtol;
    INT refine_limit, horder;
    INT dbg_level;		/* verbosity (-xfem_dbg_evel) */
    int ortho;
    BOOLEAN dbg_vtk;
} XFEM_OPTIONS;

/* ORTHO_AUTO <=>
 * 	ORTHO_CUT for Qk (nodal),
 * 	ORTHO_ALL for DGPk (modal),
 * 	ORTHO_NONE for Pk (H1 conforming) */
enum {ORTHO_NONE = 0, ORTHO_CUT = 1, ORTHO_ALL = 2, ORTHO_AUTO = 3};

#ifdef PHG_TO_P4EST
# if 0
/* for Pk_KIND == 0 */
# define OrthoType(u)	(Opts.ortho != ORTHO_AUTO ? Opts.ortho : \
			 ((DOF *)(u))->type->is_nodal ? ORTHO_CUT : ORTHO_ALL)
# else
/* for Pk_KIND == 1 */
# define OrthoType(u)	(Opts.ortho != ORTHO_AUTO ? Opts.ortho : ORTHO_CUT)
# endif
#else	/* defined(PHG_TO_P4EST) */
#define OrthoType(u)	(((DOF *)(u))->type->fe_space != FE_L2 ? ORTHO_NONE : \
				      Opts.ortho == ORTHO_AUTO ? ORTHO_CUT : \
								 Opts.ortho)
#endif	/* defined(PHG_TO_P4EST) */

#ifndef PHG_TO_P4EST
static const char *ortho_keywords[] = {"none", "cut", "all", "auto", NULL};
XFEM_OPTIONS phg_xfem_options = {
	/* ctol */	0.5,
	/* dtol */	0.05,
	/* etol */	0.25,
	/* vtol */	0.0,
	/* refine_limit */ 20,
	/* horder */	1,
	/* dbg_level */	0,
	/* ortho */	ORTHO_AUTO,
	/* dbg_vtk */	FALSE
};
#else	/* !defined(PHG_TO_P4EST) */
extern XFEM_OPTIONS phg_xfem_options;
#endif	/* !defined(PHG_TO_P4EST) */
#define Opts	phg_xfem_options

static FLOAT
get_curvature(XFEM_INFO *xi, ELEMENT *e, int np, const FLOAT pts[][Dim])
/* Returns approx(max(|k1|,|k2|)), k1 and k2 are the principal curvatures.
 * Note: the formulas for computing curvatures were taken from:
 *	Eric Albin, Ronnie Knikker, Shihe Xin, Christian Oliver Paschereit,
 *	Yves d¡¯Angelo.
 *
 *	Computational assessment of curvatures and principal
 *	directions of implicit surfaces from 3D scalar data.
 *
 *	Lecture Notes in Computer Science, Springer, 2017,
 *	Mathematical Methods for Curves and Surfaces, 10521, pp.1-22.
 *	10.1007/978-3-319-67885-6_1. hal-01486547
 * https://www.archives-ouvertes.fr/hal-01486547/document
 *
 * The array "pts[np]" contains intersections of \Gamma with the edges of K.
 */
{
    int i, j, k;
#ifdef PHG_TO_P4EST
    FLOAT xref[] = {0.5, 0.5, 0.5};
#else
    FLOAT xref[] = {0.25, 0.25, 0.25, 0.25};
#endif
    FLOAT n[Dim], u[Dim], v[Dim], d, kh, kk, kmax;
    FLOAT hess[Dim][Dim], Luu, Luv, Lvv;

    if (np == 0) {
	pts = NULL;
	np = 1;
    }

    kmax = 0.0;
    for (i = 0; i < np; i++) {
	if (pts != NULL)
	    phgGeomXYZ2Lambda(xi->ls->g, e,
			pts[i][0], pts[i][1], pts[i][2], xref);
	phgDofEval(xi->ls_grad, e, xref, n);
	d = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
	/* The case below may occur when the surface is degenerated,
	 * or the barycenter is used (when pts == NULL) */
	if (d < FLOAT_EPSILON) {
#if 0
# warning
:wq
/* ./xfem_test -xfem_ctol 1000.0 -refine 4 -xfem_dbg_vtk -xfem_etol 0 */
	    phgInfo(-1, "*** %s: WARNING: interface possibly degenerated, "
			"element marked with -9999\n", __func__);
	    e->region_mark = -9999;
#else
	    static FLOAT d_min = FLOAT_EPSILON;
#if USE_OMP
# pragma omp threadprivate(d_min)
#endif	/* USE_OMP */
	    if (d < d_min) {
		/* note: pts may be NULL */
		FLOAT x[3];
		phgGeomLambda2XYZ(xi->ls->g, e, xref, x, x + 1, x + 2);
		phgInfo(-1, "*** %s: WARNING: the interface is possibly "
			    "degenerated at [%0.16g,%0.16g,%0.16g], "
			    "ignoring ctol for the element.\n",
			    __func__, x[0], x[1], x[2]);
		d_min = d;
	    }
#endif
	    return /*1e10*/0.0;
	}
	d = 1.0 / Sqrt(d);	/* d := 1/|\grad L| = 1/|Ln| */
	n[0] *= d;
	n[1] *= d;
	n[2] *= d;

	u[0] = n[0];
	u[1] = -n[2];
	u[2] = n[1];

	v[0] = n[1] * u[2] - n[2] * u[1];
	v[1] = n[2] * u[0] - n[0] * u[2];
	v[2] = n[0] * u[1] - n[1] * u[0];

	/* Note:
	 * 	L_n = (\grad L)^T n = (\grad L)^T\grad L/|\grad L| = |\grad L|
	 *
	 * 	L_uu = u^T \hess(L) u
	 * 	L_uv = u^T \hess(L) v
	 * 	L_vv = v^T \hess(L) v */
	phgDofEval(xi->ls_hess, e, xref, hess[0]);
	Luu = Luv = Lvv = 0.0;
	for (j = 0; j < Dim; j++) {
	    for (k = 0; k < Dim; k++) {
		Luu += u[j] * hess[j][k] * u[k];
		Luv += u[j] * hess[j][k] * v[k];
		Lvv += v[j] * hess[j][k] * v[k];
	    }
	}
	/* The mean cavature:
	 * kh := (Luu + Lvv) / (2 * |L_n|) = (Luu + Lvv) * d / 2 */
	kh = (Luu + Lvv) * d * 0.5;
	/* The Gaussian cavature:
	 * kk := (Luu * Lvv - Luv^2) / |L_n|^2 = (Luu * Lvv - Luv^2) * d^2 */
	kk = (Luu * Lvv - Luv * Luv) * d * d;
	d = kh * kh - kk;
	if (d < 0.0)
	    d = 0.0;
	if (kh > 0.0)
	    d = kh + Sqrt(d);
	else
	    d = -kh - Sqrt(d);
	if (kmax < d)
	    kmax = d;
	if (xi->ls_hess->poly_order == 0)
	    break;	/* constant Hessian ==> constant curvatures */
    }

    return kmax;
}

static struct {
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

static void
ortho_aux(FLOAT x, FLOAT y, FLOAT z, DOF *dof, FLOAT *values)
/* auxiliary function for ortho_ls() and ortho_ls_grad() */
{
    GRID *g = dummy_ctx.ls->g;
    INT *mlist = (void *)dummy_ctx.e;
    ELEMENT *e = g->elems[mlist[0]];
    FLOAT xref[Dim + 1];


    phgGeomXYZ2Lambda(g, e, x, y, z, xref);
    phgDofEval(dof, e, xref, values);
}

static void
ortho_ls(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    ortho_aux(x, y, z, dummy_ctx.ls, values);
}

static void
ortho_ls_grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    ortho_aux(x, y, z, dummy_ctx.grad, values);
}


static int horder = 0;

static void
resolve_interface(XFEM_INFO *xi, int quad_order)
/* refines the mesh such that constraints on the intersection types and
 * local curvature are satisfied.
 */
{
    GRID *g = xi->ls->g;
    ELEMENT *e;
    DOF *tmp;
    double t0;
    int i;

    if (Opts.ctol <= 0.0)
	return;			/* do nothing */

    horder = Opts.horder;
    assert(xi->ls_grad->poly_order >= 0 || horder >= 0);
    if (xi->ls_grad->poly_order >= 0)
	horder = xi->ls_grad->poly_order - 1;
    tmp = phgDofCopy(xi->ls_grad, NULL, DOF_DGn[horder + 1], "G");
    xi->ls_hess = phgDofGradient(tmp, NULL, NULL, "Hessian");
    phgDofFree(&tmp);
    xi->ls_hess->userfunc = DofInterpolation;
    if (Opts.dbg_level > 0)
	phgPrintf("*** %s: ctol = %lg\n", __func__, (double)Opts.ctol);
    t0 = phgGetTime(NULL);
    while (TRUE) {
	double t = phgGetTime(NULL);
	INT n0, count = 0;	/* # of elements marked for refinement */
	ForAllElements(g, e) {
	    int np;		/* # of recorded intersection points */
	    FLOAT pts[2 * NEdge][Dim];
	    FLOAT corners[NVert][Dim];

	    dummy_ctx.ls = xi->ls;
	    dummy_ctx.grad = xi->ls_grad;
	    dummy_ctx.e = e;
	    phgGeomGetCorners(g, e, corners);

	    e->mark = 0;

	    /* First, a quick check on the element */
	    if (phgQuadInterfaceMarkElement(dummy_ls, xi->ls->poly_order,
			dummy_grad, e->elem_type, corners) != 0)
		continue;	/* no intersection is possible */

	    /* Check intersection with the edges. TODO: maintain a list of
	     * processed edges to save and reuse the cuts. */
	    np = 0;
	    for (i = 0; i < NEdge; i++) {
		int i0 = GetEdgeVertex(i, 0), i1 = GetEdgeVertex(i, 1);
#ifndef PHG_TO_P4EST
		FLOAT *v0 = corners[i0], *v1 = corners[i1];
#else	/* PHG_TO_P4EST */
		FLOAT v0[Dim] = {corners[(i0 >> 0) & 1][0],
				 corners[(i0 >> 1) & 1][1],
				 corners[(i0 >> 2) & 1][2]};
		FLOAT v1[Dim] = {corners[(i1 >> 0) & 1][0],
				 corners[(i1 >> 1) & 1][1],
				 corners[(i1 >> 2) & 1][2]};
#endif	/* PHG_TO_P4EST */
		FLOAT *cuts;
		int n;
		n = phgQuadInterfaceGetEdgeCuts(dummy_ls, xi->ls->poly_order,
						dummy_grad, v0, v1, 0.0, &cuts);
		if (n == 0) {
		    phgFree(cuts);
		    continue;
		}
		assert(np + n <= sizeof(pts) / sizeof(pts[0]));
		while (--n >= 0) {
		    pts[np][0] = v0[0] + cuts[n] * (v1[0] - v0[0]);
		    pts[np][1] = v0[1] + cuts[n] * (v1[1] - v0[1]);
		    pts[np][2] = v0[2] + cuts[n] * (v1[2] - v0[2]);
		    np++;
		}
		phgFree(cuts);
	    }

	    /* Check curvature */
	    if (phgGeomGetDiameter(g, e) * get_curvature(xi, e, np, pts)
								<= Opts.ctol)
		continue;
#if 0
#warning
FLOAT (*corners)[Dim] = phgGeomGetCorners(g, e, NULL);
phgInfo(-1, "**** element [[%g;%g;%g],[%g;%g;%g]], K = %g, h = %g\n", corners[0][0], corners[0][1], corners[0][2], corners[1][0], corners[1][1], corners[1][2], get_curvature(xi, e, np, pts), phgGeomGetDiameter(g, e));
#endif
	    e->mark = 1;
	    count += e->mark;
	}
#if USE_MPI
	if (g->nprocs > 1) {
	    INT count0 = count;
	    MPI_Allreduce(&count0, &count, 1, PHG_MPI_INT, MPI_SUM, g->comm);
	}
#endif	/* USE_MPI */
	if (Opts.dbg_level > 0)
	    phgPrintf("*** %s: %"dFMT" elements to be refined, %0.2lgs.\n",
			__func__, count, phgGetTime(NULL) - t);
	if (count == 0)
	    break;
	n0 = g->nelem_global;
	phgRefineMarkedElements(g);
	phgBalanceGrid(g, 1.1, 0, NULL, 0.);
	if (n0 == g->nelem_global) {
	    /* marked elements not refined (p4est bug?). Test case:
	     * 	  ./interface-lv -xfem_ctol 4.0 -refine0 2 -xfem_dbg_level 1
	     * (unit cube, LS = (4 * (x-xc)^2 + 2 * (y-yc)^2 + (z-zc)^2 - 1/50)^3 - (3 * (x-xc)^2 + 2 * (y-yc)^2) * (z-zc)^3 */
	    phgPrintf("*** %s: [FIXME] marked elements not refined, "
			"stop resolving interface.\n", __func__);
	    break;
	}
    }
    if (Opts.dbg_level > 0)
	phgPrintf("*** %s: %"dFMT" elements, %0.2lgs.\n",
			__func__, g->nelem_global, phgGetTime(NULL) - t0);
    phgDofFree(&xi->ls_hess);
}

static void
mark_elements(XFEM_INFO *xi, int quad_order, DOF *old_info)
/* computes and stores element mark and volume fraction in xi->info[].mark
 * and xi->info[].frac.
 *
 * This function also initializes xi->info[*].macro[*] to -1, sets up
 * xi->work[*] and xi->info[*], and sets e->mark = 0.
 *
 * The argument "old_info", if not NULL, provides index, rank and generation
 * of the element or its parent in the mesh before refinement.
 */
{
    GRID *g = xi->ls->g;
    ELEMENT *e;
    XFEM_INFO *xi0;
    /* statistics on the minimum fractions of the subdomains */
    double tmp[2][phgMaxThreads], t0;
    INT copy_count = 0;
    INT etmp[3][phgMaxThreads], nlocal;       /* ecnts: 0=- 1=0 2=+ */
    int ii;

    /* avoid recomputing if the mesh is unchanged */
    if (xi->marked && xi->serial_no == g->serial_no)
	return;

    /* Note: using HALO_FACE below the following cases fail (p4est bug?):
     *	mpirun -np {10|20} ./interface -refine0=2 -xfem_ctol=0.5
     *	src/p4est_mesh.c:502: mesh_iter_face: Assertion `side2->is.hanging.quad[h] != ((void *)0)' failed. */
    phgSetupHalo(g, HALO_VERT);

    xi->serial_no = g->serial_no;
    xi->marked = TRUE;

    bzero(xi->minfrac0, sizeof(xi->minfrac0));
    bzero(xi->minfrac, sizeof(xi->minfrac));
    bzero(xi->ecnts, sizeof(xi->ecnts));

    /* allocate and initialize xi->info */
    if (xi->work == NULL && xi->info != NULL) {
	phgFree(xi->info);
	xi->info = NULL;
    }
    xi0 = phgCalloc(1, sizeof(*xi0));
    xi0->info = xi->info;
    xi0->work = xi->work;
#ifdef PHG_TO_P4EST
    xi->info = phgCalloc(g->nelem + g->nghost, sizeof(*xi->info));
    xi->work = phgCalloc(g->nelem + g->nghost, sizeof(*xi->work));
#else	/* defined(PHG_TO_P4EST) */
    xi->info = phgCalloc(g->nelem, sizeof(*xi->info));
    xi->work = phgCalloc(g->nelem, sizeof(*xi->work));
#endif	/* defined(PHG_TO_P4EST) */

    /* mark all elements as uninitialized */
    ForAllElements_(g, e, FOR_ALL)
	xi->work[e->index].frac = 2.0;

    /* copy available data from old_info to save some recomputations */
    ForAllElements(g, e) {
	int rank, generation;
	INT index;
#ifdef PHG_TO_P4EST
	xi->work[e->index].is_small[0] = xi->work[e->index].is_small[1] = 0;
	assert(sizeof(xi->info->vsign) * 8 >= 2 * NVert);
#endif	/* defined(PHG_TO_P4EST) */
	if (xi0->info == NULL || old_info == NULL)
	    continue;
	/* copy xi->info[] entries for elements which have not been refined */
	index = (INT)(old_info->data[e->index * 2 + 1] + 0.5);
	rank = (int)(old_info->data[e->index * 2] + 0.5);
	generation = rank % 256;
	rank /= 256;
	if (rank != g->rank || generation != e->generation)
	    continue;
	memcpy(&xi->info[e->index], &xi0->info[index], sizeof(xi0->info[0]));
	memcpy(&xi->work[e->index], &xi0->work[index], sizeof(xi0->work[0]));
	assert(xi->work[e->index].frac >= 0 && xi->work[e->index].frac <= 1.0);
	copy_count++;
    }
    phgFree(xi0->info);
    phgFree(xi0->work);
    phgFree(xi0);

    phgQuadInterfaceMarkGrid(xi->ls, xi->ls_grad);

    for (ii = 0; ii < phgMaxThreads; ii++) {
	tmp[0][ii] = tmp[1][ii] = 2.0;
	etmp[0][ii] = etmp[1][ii] = etmp[2][ii] = 0;
    }

#ifdef PHG_TO_P4EST
    nlocal = g->nelem;
#else	/* defined(PHG_TO_P4EST) */
    nlocal = g->halo == NULL ? g->nelem : g->halo->nelem_bak;
#endif	/* defined(PHG_TO_P4EST) */

    t0 = phgGetTime(NULL);
#if USE_OMP
#pragma omp parallel for private(e)
#endif  /* USE_OMP */
    ForAllElementsBegin_(g, e, FOR_ALL) {
	FLOAT *p, val;
	INT eno = e->index;
	int i, n;
	BOOLEAN initialized = (xi->work[eno].frac < 1.5);
#if defined(PHG_TO_P4EST)
	/* initalize "macro" members */
	xi->info[eno].macro[0] = xi->info[eno].macro[1] = -1;
#endif	/* defined(PHG_TO_P4EST) */

	/* set up "frac", "nv" and update "mark", "vsign". */

	if (!initialized) {
	    FLOAT *rule_m = NULL, *rule_p = NULL;
	    if (e->mark == 0) {
		phgQuadInterface(xi->ls, xi->ls_grad, e, quad_order,
					&rule_m, NULL, &rule_p);
		if (phgQuadInterfaceRuleInfo(rule_p, NULL,NULL,NULL) == 0)
		    e->mark = -1;
		else if (phgQuadInterfaceRuleInfo(rule_m, NULL,NULL,NULL) == 0)
		    e->mark = 1;
	    }
	    xi->info[eno].mark = e->mark < 0 ? -1 : (e->mark > 0 ? 1 : 0);
	    if (xi->info[eno].mark == 0) {
		n = phgQuadInterfaceRuleInfo(rule_m, NULL, &p, NULL);
		for (i = 0, val = 0.; i < n; i++, p += Dim + 1)
		    val += p[Dim];
		xi->work[eno].frac = (val /= phgGeomGetVolume(g, e));
		if (xi->work[eno].frac < Sqrt(FLOAT_EPSILON))
		    xi->work[eno].frac = 0.0;
		else if (xi->work[eno].frac > 1.0 - Sqrt(FLOAT_EPSILON))
		    xi->work[eno].frac = 1.0;
	    }
	    else {
		xi->work[eno].frac = xi->info[eno].mark < 0 ? 1.0 : 0.0;
	    }
	    phgFree(rule_m);
	    phgFree(rule_p);
	}

	e->mark = 0;

	if (e->index < nlocal) {
	    /* statistics about elements in different subdomains */
	    if (xi->info[eno].mark < 0) {
		if (tmp[0][phgThreadId] > 1.0)
		    tmp[0][phgThreadId] = 1.0;
		etmp[0][phgThreadId]++;
	    }
	    else if (xi->info[eno].mark > 0) {
		if (tmp[1][phgThreadId] > 1.0)
		    tmp[1][phgThreadId] = 1.0;
		etmp[2][phgThreadId]++;
	    }
	    else {
		val = xi->work[eno].frac;
		if (val > 0. && tmp[0][phgThreadId] > val)
		    tmp[0][phgThreadId] = val;
		val = 1.0 - val;
		if (val > 0. && tmp[1][phgThreadId] > val)
		    tmp[1][phgThreadId] = val;
		etmp[1][phgThreadId]++;
	    }
	}

    } ForAllElementsEnd

    xi->minfrac0[0] = xi->minfrac0[1] = 1.0;
    xi->ecnts[0] = xi->ecnts[1] = xi->ecnts[2] = 0;
    for (ii = 0; ii < phgMaxThreads; ii++) {
	if (xi->minfrac0[0] > tmp[0][ii])
	    xi->minfrac0[0] = tmp[0][ii];
	if (xi->minfrac0[1] > tmp[1][ii])
	    xi->minfrac0[1] = tmp[1][ii];
	xi->ecnts[0] += etmp[0][ii];
	xi->ecnts[1] += etmp[1][ii];
	xi->ecnts[2] += etmp[2][ii];
    }
#if USE_MPI
    if (g->nprocs > 1) {
	double tmp[2];
	INT etmp0[4], etmp1[4];
	tmp[0] = xi->minfrac0[0];
	tmp[1] = xi->minfrac0[1];
	MPI_Allreduce(tmp, xi->minfrac0, 2, MPI_DOUBLE, MPI_MIN, g->comm);
	etmp0[0] = xi->ecnts[0];
	etmp0[1] = xi->ecnts[1];
	etmp0[2] = xi->ecnts[2];
	etmp0[3] = copy_count;
	MPI_Allreduce(etmp0, etmp1, 4, PHG_MPI_INT, MPI_SUM, g->comm);
	xi->ecnts[0] = etmp1[0];
	xi->ecnts[1] = etmp1[1];
	xi->ecnts[2] = etmp1[2];
	copy_count = etmp1[3];
    }
#endif	/* USE_MPI */

    xi->minfrac[0] = xi->minfrac0[0];
    xi->minfrac[1] = xi->minfrac0[1];

    if (Opts.dbg_level > 0)
	phgPrintf("*** %s: %"dFMT" new entries, %0.2lgs.\n", __func__,
			g->nelem_global - copy_count, phgGetTime(NULL) - t0);
}

#ifdef PHG_TO_P4EST
# define Anchor(xi, g, e, sign) \
    (xi->info[e->index].macro[sign] < 0 ? e : \
		g->elems[xi->mlist[xi->info[e->index].macro[sign]][0]])
#else	/* defined(PHG_TO_P4EST) */
# define Anchor(xi, g, e, sign) (e)
#endif	/* defined(PHG_TO_P4EST) */


/*---------------------------------------------------------------------------*/

XFEM_INFO *
phgXFEMInit(DOF *ls, DOF *ls_grad, int ls_order, int quad_order)
/* This function creats and returns an XFEM_INFO object, and prepares the mesh
 * for XFEM computation by performing the following steps:
 *   1.	resolve the interface,
 *   2.	merge small slements.
 *
 * "ls"		is the levelset function
 * "ls_grad"	is the gradient of the levelset function (may be NULL)
 * "ls_order"	is the polynomial order of of the levelset function (only used
 *		if "ls"->type == DOF_ANALYTIC && "ls_order" > 0)
 * "quad_order"	is the numerical quadrature order for computing volume ratios.
 */
{
    XFEM_INFO *xi;
    DOF *old_info = NULL;
    double t0;
#ifdef PHG_TO_P4EST
#else	/* defined(PHG_TO_P4EST) */
    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	if (ls != NULL)
	    phgError(1, "%s must be called with ls=NULL once before phgInit.\n",
		     __func__);
	initialized = TRUE;
	phgOptionsRegisterTitle("\nXFEM options:", "\n", "xfem");
	phgOptionsRegisterFloat("-xfem_ctol",
				"Tolerance on the curvature (h*kappa)",
				&Opts.ctol);
	phgOptionsRegisterFloat("-xfem_dtol",
				"Tolerance on the interface deviation "
				"(unenforced)", &Opts.dtol);
	phgOptionsRegisterFloat("-xfem_etol",
				"Edge fraction for small elements",
				&Opts.etol);
	phgOptionsRegisterFloat("-xfem_vtol",
				"Volume fraction for small elements",
				&Opts.vtol);
	phgOptionsRegisterInt("-xfem_horder",
				"Polynomial order for approximating Hessian of "
				"the levelset function", &Opts.horder);
	phgOptionsRegisterInt("-xfem_refine_limit",
				"Limit on the number of refinement passes when "
				"merging elements", &Opts.refine_limit);
	phgOptionsRegisterKeyword("-xfem_ortho",
				"Orthogonalize FE basis functions",
				ortho_keywords, &Opts.ortho);
	phgOptionsRegisterInt("-xfem_dbg_level",
				"Level of debugging messages", &Opts.dbg_level);
	phgOptionsRegisterNoArg("-xfem_dbg_vtk",
				"Output VTK file for debugging", &Opts.dbg_vtk);
	return NULL;
    }
#endif	/* defined(PHG_TO_P4EST) */

    assert(ls != NULL);

    if (ls_order > 0) {
	if (ls->type == DOF_ANALYTIC)
	    phgDofSetPolyOrder(ls, ls_order);;
	if (ls_grad != NULL && ls_grad->type == DOF_ANALYTIC)
	    phgDofSetPolyOrder(ls_grad, ls_order - 1);;
    }

    xi = phgCalloc(1, sizeof(*xi));
    xi->ls = ls;
    if ((xi->ls_grad = ls_grad) == NULL) {
	xi->ls_grad = phgDofGradient(ls, NULL, NULL, NULL);
	xi->ls_grad->userfunc = DofInterpolation;
	xi->grad_owned = TRUE;
    }
    xi->quad_order = quad_order;
    xi->serial_no = -1;

    if (Opts.dbg_level > 0)
	phgPrintf("*** xfem_init: ctol = %g, dtol = %g, etol = %g, vtol = %g, "
		  "quad_order = %d\n", (double)Opts.ctol, (double)Opts.dtol,
		  (double)Opts.etol, (double)Opts.vtol, quad_order);

    resolve_interface(xi, quad_order);

    t0 = phgGetTime(NULL);

	mark_elements(xi, quad_order, old_info);

    if (Opts.dbg_level > 0)
	phgPrintf("*** mark/merge elements: %0.2lgs\n", phgGetTime(NULL) - t0);

    return xi;
}

void
phgXFEMFree(XFEM_INFO **xi)
{
#if defined(PHG_TO_P4EST)
    ELEMENT *e;
#endif	/* defined(PHG_TO_P4EST) */
    int i;

    if (xi == NULL || *xi == NULL)
	return;

#if defined(PHG_TO_P4EST)
    if ((*xi)->ls != NULL)
	ForAllElements((*xi)->ls->g, e)
	    e->reserved = NULL;	/* reset e->reserved */

    for (i = 0; i < (*xi)->qc_n; i++) {
	DOF *u = (*xi)->qc_list[i]->fe;
	phgDofFree(&u);
	phgQCFree(&(*xi)->qc_list[i]);
    }
    phgFree((*xi)->qc_list);

    phgFree((*xi)->mlist);
    phgFreeGrid(&(*xi)->g_mac);
#endif	/* defined(PHG_TO_P4EST) */

    for (i = 0; i < (*xi)->northo; i++) {
	int sign;
	INT nelem = (*xi)->ls->g->nelem;
#if defined(PHG_TO_P4EST)
	nelem += (*xi)->ls->g->nghost;
#endif	/* defined(PHG_TO_P4EST) */
	for (sign = 0; sign < 2; sign++) {
	    INT k;
	    if ((*xi)->ortho[i].P[sign] == NULL)
		continue;
	    for (k = 0; k < nelem; k++)
		phgFree((*xi)->ortho[i].P[sign][k]);
	    phgFree((*xi)->ortho[i].P[sign]);
	}
	phgFree((*xi)->ortho[i].P0);
    }
    phgFree((*xi)->ortho);

    phgFree((*xi)->info);
    phgFree((*xi)->work);

    if ((*xi)->grad_owned)
	phgDofFree(&(*xi)->ls_grad);

    phgFree(*xi);
    *xi = NULL;
}

static const FLOAT *
get_ortho_matrix(XFEM_INFO *xi, QCACHE *qc, int sign, INT eno)
{
    GRID *g, *g0;
    ELEMENT *e, *e0, *anchor;
    DOF *u = qc->fe;
    INT mno, *mlist, mlist0[2];	/* mno: macro element index, -1 if unmerged */
    DOF_TYPE *type_bak, *type;
    int i, j, o, m, n, nelem;
    BOOLEAN cut;
    FLOAT **P, *rule, *ptwt, *a;
#ifndef PHG_TO_P4EST
    FLOAT corners[NVert][Dim];
#endif	/* !defined(PHG_TO_P4EST) */

    if (OrthoType(u) == ORTHO_NONE)
	return NULL;

    g0 = g = xi->ls->g;
    e0 = e = g->elems[eno];
    mlist0[1] = -1;
    {
	if (OrthoType(u) != ORTHO_ALL && xi->info[eno].mark != 0)
	    return NULL;
	mno = -1;
	mlist0[0] = e->index;
	mlist = mlist0;
	nelem = 1;
    }

    /* compute orthogonalization matrix */
    type = type_bak = u->type;
    while (type->base_type != NULL)
	type = type->base_type;
    assert(type->dim == 1);
    n = type->nbas;
 
    /* check previously computed P matrix */
    for (o = 0; o < xi->northo; o++)
	if (xi->ortho[o].type == type)
	    break;
    if (o >= xi->northo) {
	INT nelem0 = g0->nelem;
#ifdef PHG_TO_P4EST
	nelem0 += g0->nghost;
#endif	/* defined(PHG_TO_P4EST) */
	xi->ortho = phgRealloc_(xi->ortho, (o + 1) * sizeof(*xi->ortho),
				    o * sizeof(*xi->ortho));
	xi->ortho[o].type = type;
	xi->ortho[o].P[0] = phgCalloc(nelem0, sizeof(*xi->ortho[o].P[0]));
	xi->ortho[o].P[1] = phgCalloc(nelem0, sizeof(*xi->ortho[o].P[1]));
	xi->ortho[o].P0 = NULL;
	xi->northo = o + 1;
    }

    anchor = Anchor(xi, g0, e0, sign);
    if (mno < 0 && xi->info[e0->index].mark != 0) {
	/* same matrix for all non-interface elements */
	assert(e0 == anchor);
	P = &xi->ortho[o].P0;
	cut = FALSE;
    }
    else {
	/* interface or macro element */
	P = &xi->ortho[o].P[sign][anchor->index];
	cut = TRUE;
    }
    if (*P != NULL)
	return *P;	/* matrix already computed */

    dummy_ctx.ls = xi->ls;
    dummy_ctx.grad = xi->ls_grad;
#if 1	/*--------------------------------------------------------------------*/
    /* build the quadrature rule using the original elements */
    dummy_ctx.e = (void *)mlist0;
    rule = NULL;
    for (i = 0; i < nelem; i++) {
	FLOAT *rule0;
	e0 = g0->elems[mlist0[0] = mlist[i]];
# ifdef PHG_TO_P4EST
	phgQuadInterfaceCuboid(cut ? ortho_ls : NULL, DofTypeOrder(xi->ls, e0),
			   ortho_ls_grad, e0->corners,
			   cut ? xi->quad_order : 2 * type->order,
			   sign == 0 ? &rule0 : NULL, NULL,
			   sign == 1 ? &rule0 : NULL);
# else	/* defined(PHG_TO_P4EST) */
	phgGeomGetCorners(g0, e0, corners);
	phgQuadInterfaceTetrahedron(cut ? ortho_ls : NULL,
				DofTypeOrder(xi->ls,e0), ortho_ls_grad, corners,
				cut ? xi->quad_order : 2 * type->order,
				sign == 0 ? &rule0 : NULL, NULL,
				sign == 1 ? &rule0 : NULL);
# endif	/* defined(PHG_TO_P4EST) */
	phgQuadInterfaceRuleMerge(&rule, rule0);
	phgFree(rule0);
    }
#else	/*--------------------------------------------------------------------*/
    /* build the quadrature rule using the macro element */
    Unused(nelem);
    dummy_ctx.e = (void *)mlist;
# ifdef PHG_TO_P4EST
    phgQuadInterfaceCuboid(cut ? ortho_ls : NULL, DofTypeOrder(xi->ls,e0),
			   ortho_ls_grad, e->corners,
			   cut ? xi->quad_order : 2 * type->order,
			   sign == 0 ? &rule : NULL, NULL,
			   sign == 1 ? &rule : NULL);
# else	/* defined(PHG_TO_P4EST) */
    phgGeomGetCorners(g, e, corners);
    phgQuadInterfaceTetrahedron(cut ? ortho_ls : NULL,
				DofTypeOrder(xi->ls, e), ortho_ls_grad, corners,
				cut ? xi->quad_order : 2 * type->order,
				sign == 0 ? &rule : NULL, NULL,
				sign == 1 ? &rule : NULL);
# endif	/* defined(PHG_TO_P4EST) */
#endif	/*--------------------------------------------------------------------*/

    m = phgQuadInterfaceRuleInfo(rule, NULL, &ptwt, NULL);
    if (m == 0) {
	phgFree(rule);
	return NULL;
    }
    assert(m >= n);

#define mA(i,j) a[(j) * m + (i)]

    if (Opts.dbg_level > 2)
	phgInfo(-1, "*** %s: sign %d, macro %d: m,n = %d,%d\n",
		    __func__, sign, mno, m, n);
    a = phgAlloc(m * n * sizeof(*a));
    /* compute a = [\phi_j(x_i)] */
    for (i = 0; i < m; i++, ptwt += Dim + 1) {
	FLOAT sqrtw = Sqrt(ptwt[Dim]), xref[Dim + 1];
	phgGeomXYZ2Lambda(g, e, ptwt[0], ptwt[1], ptwt[2], xref);
	/* WARNING: inconsistency between e and qc if e is a macro element! */
	u->type = type;     /* hacky! */
	qc->qd->bas(qc, e->index, xref, &mA(i,0), m);
	u->type = type_bak;
	for (j = 0; j < n; j++)
	    mA(i,j) *= sqrtw;
    }
    phgFree(rule);
#undef mA
    *P = phgOrthoMatrixGen(n, m, a);
    phgFree(a);

    return *P;
}

typedef struct {
    XFEM_INFO	*xi;
    QCACHE	*qc;	/* QC for the macro elements */
    int		subdomain_no;
} QCACHE_CTX;

static const FLOAT *
get_macro_element(QCACHE **qc, INT *eno, const FLOAT xref[], FLOAT buffer[])
/* sets "*eno" to the macro element number if it's part of a macro element,
 * sets "*qc" to the QCACHE object for xi->g_mac (stored in qc->user_ctx),
 * stores the coordinates in the macro element in buffer[], and returns buffer,
 *
 * otherwise "*eno" and "*qc" are unchanged and xref[] is returned.
 *
 * "buffer" is of size qc->qd->ref_dim, provided by the caller function.
 */
{
#ifndef PHG_TO_P4EST
    return xref;
#else	/* !defined(PHG_TO_P4EST) */
    FLOAT x, y, z;
    QCACHE_CTX *qctx = (*qc)->user_ctx;
    XFEM_INFO *xi = qctx->xi;
    GRID *g;
    ELEMENT *e, *e1;
    INT j;
    int sign = qctx->subdomain_no;

    if ((j = xi->info[*eno].macro[sign]) < 0)
	return xref;

    g = xi->ls->g;
    e = g->elems[*eno];
    e1 = xi->g_mac->elems[j];
    phgGeomLambda2XYZ(g, e, xref, &x, &y, &z);
    phgGeomXYZ2Lambda(xi->g_mac, e1, x, y, z, buffer);

    *qc = qctx->qc;
    *eno = j;

    return buffer;
#endif	/* !defined(PHG_TO_P4EST) */
}

#define BaseTypeNBAS(type)	((type)->nbas / (type)->dim)

static void
hooked_bas(QCACHE *qc, INT eno, const FLOAT xref0[], FLOAT *buffer, int stride)
{
    QCACHE_CTX *qctx = qc->user_ctx;
    XFEM_INFO *xi = qctx->xi;
    int sign = qctx->subdomain_no;
    int n = BaseTypeNBAS(((DOF *)qc->fe)->type);
    int m = ((DOF *)qc->fe)->type->dim;
    FLOAT tmp[qc->qd->ref_dim];
    const FLOAT *xref;
    const FLOAT *P;	/* the orthogonalization matrix */

    assert(qc->user_ctx != NULL);
    P = get_ortho_matrix(xi, qc, sign, eno);
    xref = get_macro_element(&qc, &eno, xref0, tmp);
    qc->qd->bas(qc, eno, xref, buffer, stride);
    phgOrthoMatrixApply(n, m, stride, 1, buffer, P);
}

static void
hooked_grd(QCACHE *qc, INT eno, const FLOAT xref0[], FLOAT *buffer, int stride)
{
    QCACHE_CTX *qctx = qc->user_ctx;
    XFEM_INFO *xi = qctx->xi;
    int sign = qctx->subdomain_no;
    int n = BaseTypeNBAS(((DOF *)qc->fe)->type);
    int m = ((DOF *)qc->fe)->type->dim;
    FLOAT tmp[qc->qd->ref_dim];
    const FLOAT *xref;
    const FLOAT *P;	/* the orthogonalization matrix */

    assert(qc->user_ctx != NULL);
    P = get_ortho_matrix(xi, qc, sign, eno);
    xref = get_macro_element(&qc, &eno, xref0, tmp);
    qc->qd->grd(qc, eno, xref, buffer, stride);
    phgOrthoMatrixApply(n, m, stride, Dim, buffer, P);
}

#define CheckGrid(xi) if ((xi)->serial_no != (xi)->ls->g->serial_no) \
    phgError(1, "%s:%d: grid changed after phgXFEMInit.\n", __FILE__, __LINE__);

QCACHE *
phgXFEMNewQC(XFEM_INFO *xi, void *fe, int subdomain_no)
/* Note: subdomain_no == 0: \Omega^-, subdomain_no == 1: \Omega^+ */
{
    QCACHE *qc = phgQCNew(QD_DEFAULT, fe);
    QCACHE_CTX *qctx = phgCalloc(1, sizeof(*qctx));
    DOF *u = fe;

    qctx->xi = xi;
    qctx->subdomain_no = subdomain_no;

    qc->user_ctx = qctx;
    qc->owns_user_ctx = TRUE;

    qc->bas_hook = hooked_bas;
    qc->grd_hook = hooked_grd;


    CheckGrid(xi)

    if (Opts.ortho != ORTHO_NONE && Opts.ortho != ORTHO_AUTO &&
	u->type->fe_space != FE_L2 && u->g->rank == 0) {
	static BOOLEAN warned = FALSE;
	if (!warned)
	    phgWarning("Cannot normalize \"%s\" basis, option "
			"\"-xfem_ortho=%s\" ignored.\n", u->type->name,
			Opts.ortho == ORTHO_CUT ? "cut" : "all");
	warned = TRUE;
    }

    return qc;
}

#if defined(PHG_TO_P4EST)
INT
phgXFEMMapE2G(MAP *map, int dof_no, ELEMENT *e, int index)
{
    XFEM_INFO *xi = e->reserved;

    assert(xi == NULL || map->ndof % 2 == 0);

    if (xi != NULL) {
	CheckGrid(xi)
	/* DOFs must be ordered as (minus,plus) pairs */
	assert(((DOF *)map->dofs[2 * (dof_no / 2)])->type == 
	       ((DOF *)map->dofs[2 * (dof_no / 2) + 1])->type);
	e = Anchor(xi, xi->ls->g, e, dof_no % 2);
    }

    return PHG_to_p8est_map_e2g(map, dof_no, e, index);
}
#endif	/* defined(PHG_TO_P4EST) */


void
phgXFEMSetInitialSolution(XFEM_INFO *xi, DOF *um, DOF *up)
/* prepares DOFs on the macro elements (original elements => macro elements) */
{
    GRID *g = um->g;
    ELEMENT *e;
    DOF *u;
    INT i, j;
    int sign;


    CheckGrid(xi)

    if (OrthoType(um) == ORTHO_NONE)
	return;

#define mP(i,j) p[(j) * ((j) + 1) / 2 + (i)]
    /* convert DOF w.r.t. original basis to orthogonalized basis */
    for (sign = 0; sign < 2; sign++) {
	QCACHE *qc = phgQCNew(QD_DEFAULT, u = sign == 0 ? um : up);
	ForAllElements(g, e) {
	    int k, n = BaseTypeNBAS(u->type), dim = DofDim(u);
	    FLOAT *dof;
	    const FLOAT *p; /* the orthogonalization matrix */
	    if ((p = get_ortho_matrix(xi, qc, sign, e->index)) == NULL)
		continue;
	    dof = u->data + e->index * dim * n;
	    /* Compute P_{nxn}^{-1} * dof_[nxdim], where
	     * P is an upper triangular matrix in column-major order,
	     * while dof are column vectors in row-major order */
	    for (i = n - 1; i >= 0; i--)
		for (j = 0; j < dim; j++) {
		    for (k = i + 1; k < n; k++)
			dof[i * dim + j] -= mP(i,k) * dof[k * dim + j];
		    dof[i * dim + j] /= mP(i,i);
		}
	}
	phgQCFree(&qc);
    }
#undef mP
}


void
phgXFEMUpdateSolution(XFEM_INFO *xi, DOF *um, DOF *up)
/* updates DOFs on the merged elements (macro elements => original elements) */
{
    GRID *g = um->g;
    ELEMENT **elems, *e;
    DOF *u;
    DOF_TYPE *type;
    INT i, o;
    int sign;

    CheckGrid(xi)

    elems = phgAlloc(g->nelem * sizeof(*elems));
    for (sign = 0; sign < 2; sign++) {
	u = sign == 0 ? um : up;

#if 0
#warning
#define CHECK_UPDATE
Opts.dbg_level = 1;
printf("**************** Sign = %d, before %s:\n",sign,  __func__);
size = u->dim * DofNBas(u, e);
for (i = 0; i < g->nelem; i++) {
FLOAT *p = u->data + i * size;
/*if ((index = xi->info[i].macro[sign]) < 0) continue;*/
e = g->elems[i];
printf("macro %d, (%g %g %g)-(%g %g %g), %s[%d] =", xi->info[i].macro[sign], e->corners[0][0], e->corners[0][1], e->corners[0][2], e->corners[1][0], e->corners[1][1], e->corners[1][2], u->name, e->index);
for (int j = 0; j < size; j++) printf(" %g", p[j]);
printf("\n");
}
#endif

	/* convert DOF w.r.t. orthogonalized basis back to original basis */
	type = u->type;
	while (type->base_type != NULL)
	    type = type->base_type;
	for (o = 0; o < xi->northo; o++)
	    if (xi->ortho[o].type == type)
		break;
	if (o < xi->northo && xi->ortho[o].P[sign] != NULL) {
	    int j, k, n = type->nbas, dim = DofDim(u);
	    FLOAT dof0[dim * n], *dof;
	    const FLOAT *p;	/* the orthogonalization matrix */
	    QCACHE *qc = phgQCNew(QD_DEFAULT, u);
	    ForAllElements(g, e) {
		if ((p = get_ortho_matrix(xi, qc, sign, e->index)) == NULL)
		    continue;
		dof = u->data + e->index * dim * n;
		memcpy(dof0, dof, sizeof(dof0));
		/* Compute dof_{nxdim} := P_{nxn} * dof0_[nxdim], where
		 * P is an upper triangular matrix in column-major,
		 * while dof/dof0 are column vectors in row-major:
		 * 	dof[i][j] = sum_k P(i,k) * dof0[k][j], 
		 * 		0 <= i < n, 0 <= j < dim. */
		memset(dof, 0, sizeof(dof0));
		for (k = 0; k < n; k++)
		    for (i = 0; i <= k; i++, p++)
			for (j = 0; j < dim; j++)
			    dof[i * dim + j] += *p * dof0[k * dim + j];
	    }
	    phgQCFree(&qc);
	}
	if (OrthoType(u) != ORTHO_NONE)
	    phgDofUpdateGhost(u);

    }
    phgFree(elems);
}
