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

/* $Id: phg-to-p4est.h,v 1.57 2022/09/05 06:50:14 zlb Exp $ */

#ifndef PHG_TO_P4EST_H	/* whole file */
#define PHG_TO_P4EST_H

#if HAVE_P4EST		/* the remaining part of the file */

#ifndef PHG_TO_P4EST
# define PHG_TO_P4EST 1
# endif

#include <p8est_vtk.h>
#include <p8est_iterate.h>
#include <p8est_mesh.h>

#ifdef __cplusplus
extern "C" {
#endif

extern QDESC QD_P4EST_;
#define QD_P4EST (&QD_P4EST_)

#define p8est_topidx_t	p4est_topidx_t
#define p8est_locidx_t	p4est_locidx_t

#undef QD_DEFAULT
#define QD_DEFAULT	QD_P4EST

#undef QD_PHG
#define QD_PHG	QD_P4EST	/* for backward compatibility */

#undef NVert
#define NVert	8
#undef NEdge
#define NEdge	12
#undef NFace
#define NFace	6

#undef GetEdgeVertex
#define GetEdgeVertex(e, v)	p8est_edge_corners[e][v]

#undef GetFaceVertex
#define GetFaceVertex(f, v)	p8est_face_corners[f][v]

#undef GetEdgeNo
#define GetEdgeNo(v0, v1)	(p8est_edge_no[v0][v1])

extern INT p8est_i;
#undef ForAllElements_
#undef ForAllElementsBegin_
#undef ForAllElementsEnd
#define ForAllElements_(g, e, for_what) \
	for (p8est_i = ((for_what) == FOR_OWNED || (for_what) == FOR_ALL) ? \
				0 : g->nelem; \
	     p8est_i < ((for_what) == FOR_NONE ? -1 : g->nelem) + \
		       (((for_what) == FOR_OWNED || (for_what) == FOR_NONE) ? \
				0 : g->nghost) \
		&& (e = g->elems[p8est_i]) != NULL; \
	     p8est_i++)
#define ForAllElementsBegin_(g, e, for_what) \
	for (p8est_i = ((for_what) == FOR_OWNED || (for_what) == FOR_ALL) ? \
				0 : g->nelem; \
	     p8est_i < ((for_what) == FOR_NONE ? -1 : g->nelem) + \
		       (((for_what) == FOR_OWNED || (for_what) == FOR_NONE) ? \
				0 : g->nghost); p8est_i++) { \
	    e = g->elems[p8est_i];
#define ForAllElementsEnd }

#undef GlobalElement
#define GlobalElement(g, no) \
	((no) < 0 ? (no) : ((no) < (g)->nelem ? \
		(no) + (INT)(g)->p8est->global_first_quadrant[(g)->rank] : \
			(g)->O2Gmap[(no) - (g)->nelem]))

/* redefine DOF and DOF_TYPE */
#define DOF		P4EST_DOF
#define DOF_		P4EST_DOF_
#define DOF_TYPE_	P4EST_DOF_TYPE_
#define DOF_TYPE	P4EST_DOF_TYPE

/* redefine ELEMENT, stored as *p8est_quadrant_t.p.userdata */
#define ELEMENT P4EST_ELEMENT
typedef struct {
    /* Note: 'corners' stores the coordinates of the corners (to be removed
     * if we can figure out how to get coordinates of ghost elements) */
    FLOAT	corners[2][Dim];
    p8est_tree_t *tree;		/* pointer to the corresponding tree */
    p8est_quadrant_t *quad;	/* pointer to the corresponding quadrant */
    void	*reserved;	/* pointer reserved for internal use */
    void	*user_ctx;	/* pointer for use by user */
    INT		index, treeid;	/* local element index */
    int		mark, region_mark;
    BTYPE	bound_type[6];	/* boundary types */
    CHAR	generation;	/* level of refinement */
    ELEM_TYPE	elem_type;	/* element type (elem-info.h) */
} ELEMENT;

/* redefine GRID */
#define GRID P4EST_GRID
typedef struct {
    p8est_connectivity_t *conn;
    p8est_t *p8est;
    p8est_mesh_t *mesh;
    p8est_ghost_t *ghost;

    ELEMENT **elems;
    ELEMENT *roots;	/* allocated array of elements (optional) */
    ELEMENT *ghost_data;
    INT *O2Gmap;	/* off-process (ghost) element's global indices */
    INT *ordering;	/* O2Gmap[ordering[]] strictly increasing */
    INT *partition;	/* size nprocs+1, # elements of all processes */

    /* for ghost data exchange */
    struct {
	INT	*slist, nsend;
	int	*rcnts, *rdsps, *scnts, *sdsps;
	int	serial_no;
	BOOLEAN	do_nothing;
    } *ghost_comm;

    /* list of DOF objects (NULL terminated) */
    struct DOF_ **dofs;

    FLOAT lif, bbox[2][Dim];

#define nleaf_global	nelem_global
#define nleaf		nelem
    INT i, nelem, nelem_global, nghost, nroot;	/* iterator and # elements */
    BOOLEAN balanced; /* for only calling p8est_partition() once */

    int nprocs, rank, halo_type, serial_no;
    MPI_Comm comm;
} GRID;

typedef struct DOF_TYPE_ {
    char *name;
    struct DOF_TYPE_ *grad_type;	/* gradient type */
    struct DOF_TYPE_ *base_type;	/* base type (for derived DOF_TYPEs) */

    /* 1. If "tensor_product" == TRUE (for DGQn), then "orders" stores degrees
     *	  in each space direction as follows:
     * 		"orders"[0] the degree in x,
     * 		"orders"[1] the degree in y,
     * 		"orders"[2] the degree in z.
     *    and "points", if not NULL, stores the nodes in the following format:
     * 		"points"[("orders"[0]+1) + ("orders"[1]+1) + ("orders"[2]+1)],
     *	  the actual nodes (interpolation points) are given by:
     * 		{(px[i], py[j], pz[k])},
     * 			i in [0, "orders"[0]],
     * 			j in [0, "orders"[1]],
     * 			k in [0, "orders"[2]],
     *	  where:
     *		px = "points", py = px + "orders"[0]+1, pz = py + "orders"[1]+1,
     *	  and "work" points to an array of the same structure as "points",
     *	  which stores coefficients of the basis functions.
     *
     * 2. If "tensor_product" == FALSE (for DGPn), then "orders[i][0:2]" stores
     *	  the orders in x, y and z of the basis function i, and "points", if
     *	  not NULL, stores the full list of nodes in the format:
     *		"points"[nbas][Dim],
     *	  and "work" will store information about the basis functions.
     *
     * If "points" == NULL, then element-wise L2 projection will be used,
     * and "work"/"pvt" will store respectively the LU factorization and pivots
     * of the element mass matrix of the reference element (the unit cube).
     */
    FLOAT *points;	/* interpolation points (in ref. coord.) */
    FLOAT *work;	/* work area, depends on the basis type */
    int *orders, *pvt;

    /* dim:	space dimension of the basis (e.g. 1 for scalar, 3 for vector)
     * nbas:	# basis functions in an element
     * order:	highest polynomial order of the basis functions */
    int dim, nbas, order;

    BOOLEAN tensor_product;	/* whether is a tensor product basis */
    BOOLEAN is_nodal;		/* whether is a nodal basis */

    /* finite element space */
    FE_SPACE fe_space;
} DOF_TYPE;

#define DOF_USER_FUNC_LAMBDA	PHG_to_p8est_DOF_USER_FUNC_LAMBDA
typedef void (*DOF_USER_FUNC_LAMBDA)(struct DOF_ *u, ELEMENT *e, int bno,
                const FLOAT *xref, FLOAT *values);

typedef struct DOF_ {
    char *name;
    GRID *g;
    DOF_TYPE *type;
    /* data and ghost_data point respectively to dynamically allocated arrays
     * 		FLOAT[nelem][n_bas][dim],
     * 		FLOAT[nghost][n_bas][dim]
     */
    FLOAT *data, *ghost_data;
    FUNC_3D userfunc;	/* function pointer for DOF_ANALYTIC */
    DOF_USER_FUNC_LAMBDA userfunc_lambda;
    void *ctx;		/* for users to pass data to callback functions */
    const char *srcfile;
    int dim, poly_order, srcline;
} DOF;

#undef DOF_ANALYTIC
#define DOF_ANALYTIC	p8est_dof_analytic
#define DOF_ANALYTIC_	p8est_dof_analytic_
extern DOF_TYPE *DOF_ANALYTIC;

#undef DOF_CONSTANT
#define DOF_CONSTANT	p8est_dof_constant
#define DOF_CONSTANT_	p8est_dof_constant_
extern DOF_TYPE *DOF_CONSTANT;

#define DOF_DEFAULT	p8est_dof_default
extern DOF_TYPE *DOF_DEFAULT;

#define DOF_TYPE_ORDER_MAX	16
#define DOF_DGn	DOF_DGPn
extern DOF_TYPE **DOF_DGQn, **DOF_DGPn;

extern int p8est_edge_no[NVert][NVert];

#undef DofDataCountGlobal
#define DofDataCountGlobal(u)	\
	((u)->g->nelem_global * (INT)(u)->dim * (INT)(u)->type->nbas)
#define phgDofSetDirichletBoundaryMask(u, mask)

#undef DofTypeOrder
#define DofTypeOrder(dof, e) (SpecialDofType((dof)->type) ? \
				(dof)->poly_order : (dof)->type->order)
#undef DofDim
#define DofDim(u) ((u)->dim * (u)->type->dim)
#undef DofTypeDim
#define DofTypeDim(u) ((u)->type->dim)
#undef DofNBas
#define DofNBas(dof, e)	(dof)->type->nbas
#undef DofFESpace
#define DofFESpace(dof) (dof)->type->fe_space

/* NLIST (neighbour list) stores lists of neighbours of all or some elements */
typedef struct {
    /* flist := list of face neighbours,
     * elist := list of edge neighbours.
     * vlist := list of vertex neighbours.
     * Note:	the list of edge neighbours exclude face neighbours,
     *		the list of vertex neighbours exclude face and edge neighbours.
     */
    INT	*flist, *elist, *vlist;
    /* floc/eloc/vloc have size g->nelem + g->nghost + 1 and
     * for i, 0 <= i < g->nelem + g->nghost,
     *		floc[i] is the start position in flist[],
     *		eloc[i] is the start position in elist[],
     *		vloc[i] is the start position in vlist[],
     * of the list of edge/vertex neighbours of element i */
    INT	*floc, *eloc, *vloc;
} NLIST;

NLIST *phgNeighboursListNew(GRID *g);
void phgNeighboursListFree(NLIST **nl);

/* PHG functions are encapsulated below, with exactly the same prototypes */

#define phgQuadInterface	PHG_to_p8est_quad_interface
#define phgQuadInterfaceFace	PHG_to_p8est_quad_interface_face
#define phgQuadInterfaceMarkGrid PHG_to_p8est_mark_grid
#define phgQuadGetRule_		PHG_to_p8est_quad_get_rule_
#define phgNewGrid		PHG_to_p8est_new_grid
#define phgInit			PHG_to_p8est_init
#define phgFinalize		PHG_to_p8est_finalize
#define phgFreeGrid		PHG_to_p8est_free_grid
#define phgImport_		PHG_to_p8est_import_
#define phgRefineAllElements	PHG_to_p8est_refine_all_elements
#define phgRefineMarkedElements	PHG_to_p8est_refine_marked_elements
#ifdef phgBalanceGrid
# undef phgBalanceGrid
#endif
#define phgBalanceGrid		PHG_to_p8est_balance_grid
#define phgSetupHalo_		PHG_to_p8est_setup_halo_
#define phgSetPeriodicity(g, flags) \
	if (flags != 0) phgWarning("phgSetPeriodicity ignored.\n");
#define phgMemoryUsage		PHG_to_p8est_memory_usage

#define phgGeomGetFaceOutNormal	PHG_to_p8est_get_out_normal
#define phgGeomGetCorners	PHG_to_p8est_get_corners
#define phgGeomLambda2XYZ	PHG_to_p8est_lamdba2xyz
#define phgGeomXYZ2Lambda	PHG_to_p8est_xyz2lamdba
#define phgGeomGetDiameter	PHG_to_p8est_get_diameter
#define phgGeomGetFaceDiameter	PHG_to_p8est_get_face_diameter
#define phgGeomGetVolume	PHG_to_p8est_get_volume

#define phgGetNeighbour_	PHG_to_p8est_get_neighbour_
#define phgOppositeFace		PHG_to_p8est_opposite_face
#define phgDofNew_		PHG_to_p8est_dof_new_
#define phgDofTypeVector	PHG_to_p8est_dof_type_vector
#define phgDofTypeFree		PHG_to_p8est_dof_type_free

#define phgDofSetDataByFunction__ PHG_to_p8est_dof_set_data_by_function__
#undef phgDofSetDataByFunction_
#define phgDofSetDataByFunction_(dof, userfunc, userfunc_lambda, funcvalues) \
	phgDofSetDataByFunction__(dof, userfunc, userfunc_lambda, funcvalues, \
				  -1, NULL)

#define phgDofSetPolyOrder	PHG_to_p8est_dof_poly_order
#define phgDofSetFunction	PHG_to_p8est_dof_set_function
#define phgDofSetLambdaFunction	PHG_to_p8est_dof_set_lambda_function
#define phgDofFree		PHG_to_p8est_dof_free
#undef phgDofEval
#define phgDofEval		PHG_to_p8est_dof_eval
#define phgDofAXPBY_		PHG_to_p8est_dof_axpby_
#define phgDofCopy_		PHG_to_p8est_dof_copy_
#define phgDofGradient_		PHG_to_p8est_dof_gradient_
#define phgDofCurl_		PHG_to_p8est_dof_curl_
#define phgDofDivergence_	PHG_to_p8est_dof_divergence_
#define phgDofNormL2_		PHG_to_p8est_dof_norm_L2_
#undef phgDofUpdateGhost
#define phgDofUpdateGhost	PHG_to_p8est_dof_update_ghost
#undef phgMapE2G
#define phgMapE2G		PHG_to_p8est_map_e2g
#define phgExportVTK		PHG_to_p8est_export_vtk
#define phgMapCreate		PHG_to_p8est_map_create
#define phgMapCreateN		PHG_to_p8est_map_create_n
#define phgSolverCreate		PHG_to_p8est_solver_create
#define phgSolverSolve		PHG_to_p8est_solver_solve
#define phgSolverDumpMATLAB_	PHG_to_p8est_solver_dump_matlab

GRID *phgNewGrid(int flags);
void phgInit(int *argc, char ***argv);
void phgFinalize(void);
void phgFreeGrid(GRID **gptr);
BOOLEAN phgImport_(GRID *g, const char *filename, MPI_Comm comm);
void phgRefineAllElements(GRID *g, int level);
void phgRefineMarkedElements(GRID *g);
int phgBalanceGrid(GRID *g, FLOAT lif_threshold, INT submesh_threshold,
					DOF *weights, FLOAT power);
void phgSetupHalo_(GRID *g, int halo_type, BOOLEAN update_dofs);
size_t phgMemoryUsage(GRID *g, size_t *peak);
void phgGeomGetFaceOutNormal(GRID *g, ELEMENT *e, int face, FLOAT *nv);
FLOAT (*phgGeomGetCorners(GRID *g, ELEMENT *e, FLOAT (*cub)[Dim]))[Dim];
FLOAT phgGeomGetDiameter(GRID *g, ELEMENT *e);
FLOAT phgGeomGetFaceDiameter(GRID *g, ELEMENT *e, int face);
FLOAT phgGeomGetVolume(GRID *g, ELEMENT *e);
void  phgGeomLambda2XYZ(GRID *g ,ELEMENT *e, const FLOAT lambda[],
			FLOAT *x, FLOAT *y, FLOAT *z);
FLOAT *phgGeomXYZ2Lambda(GRID *g, ELEMENT *e, FLOAT x, FLOAT y, FLOAT z,
                         FLOAT lambda[]);
ELEMENT *phgGetNeighbour_(GRID *g, ELEMENT *e, int face,
					const char *file, int line);
BYTE phgOppositeFace(GRID *g, ELEMENT *e, BYTE v, ELEMENT *e_op);
DOF_TYPE *phgDofTypeVector(DOF_TYPE *base_type, int dim, const char *name);
void phgDofTypeFree(DOF_TYPE **ptype);
DOF *phgDofNew_(GRID *g, DOF_TYPE *type, HP_TYPE *hp, SHORT dim,
		const char *name, DOF_USER_FUNC userfunc,
		const char *srcfile, int srcline);
void phgDofSetDataByFunction__(DOF *dof, DOF_USER_FUNC func,
			  DOF_USER_FUNC_LAMBDA func_lam, FLOAT *funcvalues,
			  INT nelem, ELEMENT **elems);
void phgDofSetPolyOrder(DOF *u, SHORT order);
void phgDofSetFunction(DOF *u, DOF_USER_FUNC func);
void phgDofSetLambdaFunction(DOF *u, DOF_USER_FUNC_LAMBDA func);
void phgDofFree(DOF **pdof);
void phgDofUpdateGhost(DOF *dof);
FLOAT *phgDofEval(DOF *dof, ELEMENT *e, const FLOAT x[], FLOAT *values);
DOF *phgDofAXPBY_(FLOAT a, DOF *x, FLOAT b, DOF **py,
		  const char *srcfile, int srcline, BOOLEAN check);
DOF *phgDofCopy_(DOF *src, DOF **dest, DOF_TYPE *newtype,
			const char *name, const char *srcfile, int srcline);
DOF *phgDofGradient_(DOF *src, DOF **newdof, DOF_TYPE *newtype,
		     const char *name, const char *srcfile, int srcline);
DOF *phgDofCurl_(DOF *src, DOF **newdof, DOF_TYPE *newtype,
		 const char *name, const char *srcfile, int srcline);
DOF *phgDofDivergence_(DOF *src, DOF **newdof, DOF_TYPE *newtype,
		       const char *name, const char *srcfile, int srcline);
FLOAT phgDofNormL2_(DOF *u, int quad_order);
INT phgQuadInterfaceMarkGrid(DOF *ls, DOF *ls_grad);
INT phgMapE2G(MAP *map, int dof_id, ELEMENT *e, int index);
const char *phgExportVTK(GRID *g, const char *fn, DOF *dof, ...);
MAP *phgMapCreateN(int ndof, DOF **dofs);
MAP *phgMapCreate(DOF *u, ...);
SOLVER *phgSolverCreate(OEM_SOLVER *oem_solver, DOF *u, ...);
int phgSolverSolve(SOLVER *solver, BOOLEAN destroy, DOF *u, ...);
void phgSolverDumpMATLAB_(SOLVER *solver, const char *mat_name,
		const char *rhs_name, const char *perm_name, BOOLEAN reorder);
void phgQuadInterface(DOF *ls, DOF *ls_grad, ELEMENT *e, int order,
		      FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);
void phgQuadInterfaceFace(DOF *ls, DOF *grad, ELEMENT *e, int face, int order,
			  FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p);
FLOAT *phgQuadGetRule_(GRID *g, ELEMENT *e, int which, int order);
GRID *phgCreateSimpleGrid(INT nelem, FLOAT (*list_of_corners)[2][Dim]);

#ifdef __cplusplus
}
#endif

#else	/* HAVE_P4EST */
# ifdef PHG_TO_P4EST
#  undef PHG_TO_P4EST
# endif

#endif	/* HAVE_P4EST */

#endif	/* defined(PHG_TO_P4EST_H) */
