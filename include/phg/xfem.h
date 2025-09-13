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

/* $Id: xfem.h,v 1.296 2022/09/03 01:06:48 zlb Exp $ */

#ifndef PHG_XFEM_H	/* whole file */

#ifdef PHG_TO_P4EST

# ifndef MERGE_MAX
/* note: only macro types with size <= MERGE_MAX will be used */
#  define MERGE_MAX	27	/* 1, 2, 3, 4, 6, 8, 9, 12, 18 or 27 */
# endif
#define XFEM_INFO	XFEM_INFO_p4est

#endif	/* defined(PHG_TO_P4EST) */

typedef struct {
    /* the levelset function, defined either by DOFs or function pointers */
    DOF		*ls, *ls_grad;	/* the level set function and its gradient */
    DOF		*ls_hess;	/* scratch variable used by get_curvature() */
    struct {
	/* orthogonalization matrices, which are triangular matrices T such that
	 * (phi_1,...,phi_n)*T is orthogonal.
	 *
	 * P[sign][nelem+nghost] stores the matrices for cut elements
	 * (for merged elements only those of the anchor elements are stored).
	 *
	 * P0 stores the (same) matrix for all non-interface elements. */
	DOF_TYPE *type;		/* DOF_TYPE */
	FLOAT	*P0, **P[2];	/* orthogonalization matrices */
    } *ortho;
    int		northo;		/* size of the ortho[] */
    int		quad_order;	/* quadrature order in interface elements */

#if defined(PHG_TO_P4EST)
    GRID	*g_mac;		/* temporary GRID for the macro elements */
    QCACHE	**qc_list;	/* list of QCs to be freed by phgXFEMFree */
    void	*ctx;		/* pointer for use by callback function */

    /* mlist[i] contains the list of native elements for the macro element i,
     * filled by trailing -1's if fewer than MERGE_MAX native elements */
    INT		m, (*mlist)[MERGE_MAX];

    INT		mcnts[MERGE_MAX + 1];	/* # macro elements pairs by types */
    INT		fcnts[2];		/* # elements failed to merge */
    int		qc_n;			/* size of qc_list[] */
#endif	/* defined(PHG_TO_P4EST) */

    /* The following 2-componet arrays stores information about merges,
     * with component 0 for Omega^- and component 1 for Omega^+ */

    /* "info" is an array of size g->nelem, dynamically allocated */
    struct {
#if defined(PHG_TO_P4EST)
	INT	macro[2];	/* macro element number (-1 if not merged) */
	SHORT	vsign;		/* vertex sign, 2 bits/vert: 0, 1(-), 2(+) */
#endif	/* defined(PHG_TO_P4EST) */
	CHAR	mark;		/* -1: Omega^-, 1: Omega^+, 0: interface */
    } *info;

    /* work is only needed during merging elements */
    struct {
	FLOAT	frac;		/* volume fraction Omega^- in the elem */
#if defined(PHG_TO_P4EST)
	FLOAT	nv[Dim];	/* normal direction of the interface */
	FLOAT	ecut[NEdge];	/* edge cut, <0: no cut, >1: multiple cuts */
	CHAR	vmask[2];	/* short vertex masks, 1 bit per vertex */
	CHAR	is_small[2];	/* !0 <=> element is small and to be merged */
#endif	/* defined(PHG_TO_P4EST) */
    } *work;

    /* Statistics on the minimum fractions of elements in the subdomains:
     *		[0]:		Omega^-,
     *		[1]:		Omega^+,
     *		minfrac0:	min volume fraction before merging,
     *		minfrac:	min volume fraction after merging.
     *		maxdev:		largest interface deviation. */
    double minfrac0[2], minfrac[2], maxdev[2];
    INT ecnts[3];	/* # elems in \Omega-(0), \Gamma(1) and \Omega+(2) */

    int		serial_no;		/* serial_no of the current mesh */
    BOOLEAN	marked;			/* TRUE => grid is marked */
    BOOLEAN	grad_owned;		/* TRUE => owns ls_grad */
} XFEM_INFO;

#ifdef PHG_TO_P4EST

INT phgXFEMMapE2G(MAP *map, int dof_no, ELEMENT *e, int index);

#if !defined(DONT_CAST_phgMapE2G)
# ifdef phgMapE2G
#  undef phgMapE2G
# endif
# define phgMapE2G			phgXFEMMapE2G
#endif 	/* !defined(DONT_CAST_phgMapE2G) */

/* Append _p4est for the p4est version of functions */
#define phgXFEMInit			phgXFEMInit_p4est
#define phgXFEMFree			phgXFEMFree_p4est
#define phgXFEMSetInitialSolution	phgXFEMSetInitialSolution_p4est
#define phgXFEMUpdateSolution		phgXFEMUpdateSolution_p4est
#define phgXFEMNewQC			phgXFEMNewQC_p4est

#endif	/* defined(PHG_TO_P4EST) */

XFEM_INFO *phgXFEMInit(DOF *ls, DOF *ls_grad, int ls_order, int quad_order);
void phgXFEMFree(XFEM_INFO **xi);
QCACHE *phgXFEMNewQC(XFEM_INFO *xi, void *fe, int subdomain_no);
void phgXFEMSetInitialSolution(XFEM_INFO *xi, DOF *um, DOF *up);
void phgXFEMUpdateSolution(XFEM_INFO *xi, DOF *um, DOF *up);

#define PHG_XFEM_H
#endif	/* whole file */
