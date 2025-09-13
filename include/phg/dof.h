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

/* $Id: dof.h,v 1.31 2022/04/09 01:28:34 zlb Exp $ */

#ifndef PHG_DOF_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Data structures and functions for administration of DOFs.
 *
 * The DOFs are divided into four kinds:
 * 	those uniquely owned by a vertex (vertex DOF)
 *	those uniquely owned by an edge (edge DOF)
 *	those uniquely owned by a face (face DOF)
 *	those uniquely owned by an element (volume DOF)
 *
 * A DOF is described by np_vert, np_edge, np_face, np_elem which describe
 * the number of points per vertex, edge, face, and element respectively,
 * only 0 or 1 is allow for np_vert.
 *
 * Note: when np_edge or np_face > 1, the order of points is determined
 * by the global indices of the vertices of the edge or face!
*/

/* 'DOF_USER_FUNC' and 'DOF_USER_FUNC_LAMBDA' are respectively the prototypes
 * for 'user_func' and 'user_func_lambda' in 'DOF_TYPE.InitFunc'. They are used
 * to interpolate a given function (say 'F') to the DOF object 'u', and should
 * return in 'values' the values of 'F', and its derivatives if D_order>0, as
 * listed in the following table.
 *
 * Note:
 *    The `bno' argument is provided for saving redundant computations
 *    (compute values of basis functions only once when type->invariant==TRUE).
 *    It is for internal use only and should not used by users.
 *
 * D_order	Returned values		Ordering
 * ---------------------------------------------------------------------------
 *   0		F 			values[u->dim][u->type->dim]
 *   1		F,Fx,Fy,Fz		values[u->dim][u->type->dim][Dim + 1]
 *  >1		unimplemented
 * ---------------------------------------------------------------------------
 * (in the table (Fx,Fy,Fz) denotes the gradient of F)
 */
typedef void (*DOF_USER_FUNC)(FLOAT x, FLOAT y, FLOAT z, FLOAT *values);

/* similar to DOF_USER_FUNC, but uses lambda[] instead of xyz */
typedef void (*DOF_USER_FUNC_LAMBDA)(struct DOF_ *u, ELEMENT *e, int bno,
		const FLOAT *lambda, FLOAT *values);

/* prototype of functions assigning DOF values.
 * Note: the function should assign all DOF values at the given vertex, edge,
 * face, or element (denoted by 'type'), using either 'userfunc' or
 * 'funcvalues' (exactly one of them is non nil, 'funcvalues' should contain
 * dof->dim * dof->type->dim * dof->type->np_xxxx numbers, with 'xxxx'=
 * 'vert', 'edge', 'face', or 'elem'
 *
 * Note: if the points[] array is not NULL, then the DOF_INIT_FUNC should
 * only call the user function for the points in the points[] array (and
 * in order). Some functions (like DofCopy) rely on this assumption to
 * make use of cached basis function values. For example, schematically,
 * user_func_lambda is called like (assuming type == EDGE):
 *
 *	for (i = 0; i < np_edge; i++) {
 *	    lambda := coordinates of the i-th points of the edge;
 *	    user_func_lambda(dof, e, lambda, dof_buffer);
 *	    dof_buffer += dof->dim;
 *	}
 */
typedef void (*DOF_INIT_FUNC)(struct DOF_ *dof, ELEMENT *e,
			      GTYPE type, int index, DOF_USER_FUNC userfunc,
			      DOF_USER_FUNC_LAMBDA userfunc_lambda,
			      const FLOAT *funcvalues, FLOAT *dofvalues,
			      FLOAT **pdofvalues);

/* Temporary macro for checking unimplemented case */
#define WARN_DOF_PERMUTATION(u) \
    if (((u)->type->np_edge > 1 || (u)->type->np_face > 1) \
	&& (u)->type->points != NULL) \
	phgError(1, "%s:%d: unimplemented case for DOF \"%s\".\n", \
			__FILE__, __LINE__, (u)->type->name)

/* Macro checking whether an element is in global order,
 * i.e., its ordering of vertices matches the global numbering of vertices */
#define ElementIsInOrder(e) \
	((e)->ordering == (0 | (1 << 3) | (2 << 6) | (3 << 9)))

#define MapVertex(e, no)	(((e)->ordering >> (3 * (no))) & 7)

#define GetEdgeVertices(e, index, v0, v1)	\
    v0 = GetEdgeVertex(index, 0);		\
    v1 = GetEdgeVertex(index, 1);		\
    if (MapVertex(e, v0) > MapVertex(e, v1)) {	\
	int _tmp = v0; v0 = v1; v1 = _tmp;	\
    }

#define GetFaceVertices(e, index, v0, v1, v2, v3)			\
    phgGetFaceVertices_(e, index, &(v0), &(v1), &(v2), &(v3), NULL)
void phgGetFaceVertices_(const ELEMENT *e, int index,
			 int *v0, int *v1, int *v2, int *v3, int *map);

#define GetElementVertices(e, index, v0, v1, v2, v3,		\
			   v4, v5, v6, v7)			\
    phgGetElementVertices_(e, index, &(v0), &(v1), &(v2), &(v3),\
				     &(v4), &(v5), &(v6), &(v7))
void phgGetElementVertices_(const ELEMENT *e, int index,
			    int *v0, int *v1, int *v2, int *v3,
			    int *v4, int *v5, int *v6, int *v7);

void _phg_get_lambda_on_vert(int v0, const FLOAT *qp, FLOAT *lambda);
void _phg_get_lambda_on_edge(int v0, int v1, const FLOAT *qp, FLOAT *lambda);
void _phg_get_lambda_on_face(int v0, int v1, int v2, int v3,
			     const FLOAT *qp, FLOAT *lambda);

/*------------------- DOF function types ---------------------*/

/* DOF_INTERP_FUNC interpolates the DOF values between parent/children.

   'parent_data' is a list of pointers to DOF data for the element 'e':
	old_data[0-3]	v0, v1, v2, v3
	old_data[4-9]	e0, e1, e2, e3, e4, e5
	old_data[10-13]	f0, f1, f2, f3
	old_data[14]	volume data

   'children_data' is a list of pointers to the new locations
	new_data[0]	the new vertex (at center of e0)
	new_data[1]	the cut edge owning v0 of e
	new_data[2]	the cut edge owning v1 of e
	new_data[3]	the new edge owning v2 of e
	new_data[4]	the new edge owning v3 of e
	new_data[5]	the cut face on f2 owning v0
	new_data[6]	the cut face on f2 owning v1
	new_data[7]	the cut face on f3 owning v0
	new_data[8]	the cut face on f3 owning v1
	new_data[9]	the new face
	new_data[10]	volume data of child 0
	new_data[11]	volume data of child 1 */
typedef void (*DOF_INTERP_FUNC)(struct DOF_ *dof, ELEMENT *e,
				FLOAT **parent_data, FLOAT **children_data);

/* A DOF is described by 4 numbers: np_vert, np_edge, np_face, and np_elem,
   which represent # of points per vertex (only 0 or 1 is allowed),
   # of points per edge, # of points per face, and # of points per volume,
   respectively. The 'points' field gives barycentric coordinates of the
   distribution of points, which must have
	np_vert * 1 + np_edge * 2 + np_face * 3 + np_elem * 4
   entries.

   Data for a given DOF are stored in an array as follows,
   suppose n0 = g->nvert, n1 = n0 + g->nedge, n2 = n1 + g->nface,
   and n3 = n2 + g->nelem, D represents 'dim' FLOATs:

   Vertex data (g->nvert * np_vert * dim FLOATs):
	    D[   0][0]	... ...	D[   0][np_vert-1]
	    D[   1][0]	... ...	D[   1][np_vert-1]
	    ... ...	... ...	... ...
	    D[n0-1][0]	... ...	D[n0-1][np_vert-1]
   Edge data (g->nedge * np_edge * dim FLOATs):
	    D[n0+0][0]	... ...	D[n0+0][np_edge-1]
	    D[n0+1][0]	... ...	D[n0+1][np_edge-1]
	    ... ...	... ...	... ...
	    D[n1-1][0]	... ...	D[n1-1][np_edge-1]
   Face data (g->nface * np_face * dim FLOATs):
	    D[n1+0][0]	... ...	D[n1+0][np_face-1]
	    D[n1+1][0]	... ...	D[n1+1][np_face-1]
	    ... ...	... ...	... ...
	    D[n2-1][0]	... ...	D[n2-1][np_face-1]
   Element data (g->nelem * np_elem * dim FLOATs):
	    D[n2+0][0]	... ...	D[n2+0][np_elem-1]
	    D[n2+1][0]	... ...	D[n2+1][np_elem-1]
	    ... ...	... ...	... ...
	    D[n3-1][0]	... ...	D[n3-1][np_elem-1] */

/* function returning a pointer to the values of the basis functions
   'no0'--'no1 - 1' at a given point, in the order FLOAT[][type->dim],
   'no1 <= 0' means 'no1' == nbas */
typedef const FLOAT *(*DOF_BASIS_FUNC)(struct DOF_ *dof, ELEMENT *e, int no0,
						int no1, const FLOAT *lambda);

/* function returning a pointer to the values of the gradient with respect
   to lambda of the basis functions 'no0'--'no1 - 1' at a given point,
   in the order FLOAT[][type->dim][Dim + 1], 'no1 <= 0' means 'no1' == nbas */
typedef const FLOAT *(*DOF_BASIS_GRAD)(struct DOF_ *dof, ELEMENT *e, int no0,
						int no1, const FLOAT *lambda);

typedef struct DOF_TYPE_ DOF_TYPE;

/* Function for handling (face and element) DOF data ordering, it should
 * return in M[] the mapping (local ordering)->(DOF ordering), which maps
 * [0,np_{face,elem})->[0,np_{face,elem}).
 *
 * The function maps DOF on face `no' if no>=0 and element DOF if no<0
 */
typedef void (*DOF_MAP_FUNC)(const DOF_TYPE *type, const ELEMENT *e,
			      int no, int *M);

/* macros for use in DOF_TYPE initialization */
#define	DofCache NULL, NULL, NULL, NULL, NULL
#define	DofReserved	DofCache, \
			NULL,	/* mass_lumping = */ \
			NULL	/* initialize() */

/* finite element spaces */
typedef enum {FE_None, FE_L2, FE_Hcurl, FE_Hdiv, FE_H1, FE_Cinfty} FE_SPACE;

typedef struct HP_INFO_ HP_INFO;	/* defined later in hp.h */

struct DOF_TYPE_ {
    /*----------------------- pointers in DofCache -----------------------*/
    /* caches for 3D quadrature */
    void	*cache_basfunc, *cache_basgrad; /* cached bases, Dbas/Dlambda */
    void	*cache_gradient, *cache_curl, *cache_hessian;

    struct MASS_LUMPING_ *mass_lumping;		/* mass lumping struct */
    void (*Initialize)(struct DOF_TYPE_ *);	/* Initialization function */
    /*-----------------------------------------------------------------------*/

    const char	*name;		/* name (or description) */
    FLOAT	*points;	/* barycentric coordinates of the points,
				   (np_vert + np_edge*2 + np_face*3 + np_elem*4)
				   entries */
    BYTE	*orders;	/* polynomial order of individual basis
				   functions (may be NULL), should contain
				   np_vert + np_edge + np_face + np_elem
				   entries */

    DOF_TYPE	*grad_type;	/* DOF type of the gradients */
    DOF_TYPE	*base_type;	/* for phgDofTypeVector: the base type
				 * for DGn: the corresponding Pn/HBn type */
    HP_INFO	*hp_info;	/* NULL for non-hierarchic bases */

    /* pointers to DOF functions */
    DOF_INTERP_FUNC	InterpC2F; /* coarse-fine interpolation function */
    DOF_INTERP_FUNC	InterpF2C; /* fine-coarse interpolation function */
    DOF_INIT_FUNC	InitFunc;  /* function for assigning DOF values */
    DOF_BASIS_FUNC	BasFuncs;  /* function for evaluating basis funcs */
    DOF_BASIS_GRAD	BasGrads;  /* function for evaluating gradients */
    DOF_MAP_FUNC	DofMap;	   /* function for ordering face data or
				    * element data ("no"<0 for the latter) */

    FE_SPACE	fe_space;	/* the FE space the functions belong to */

    BOOLEAN	is_nodal;	/* TRUE ==> have nodal points */
    BOOLEAN	invariant;	/* TRUE ==> same bas functions for all elems */
    BOOLEAN	free_after_use;	/* TRUE ==> DOF_TYPE is dynamically allocated
					    and should be freed after last
					    DOF using it is freed */

    SHORT	id;		/* id (assigned when referenced) */

    SHORT	nbas;		/* number of basis functions per element:
					NVert * np_vert + NEdge * np_edge +
				 	NFace * np_face + np_elem */
    BYTE	order;		/* polynomials order */
    BYTE	D_order;	/* Order of derivatives required by InitFunc */
    CHAR	continuity;	/* continuity across element boundaries:
				   < 0 means discontinuous, 0 means C^0,
				   1 means C^1, etc. */
    SHORT	dim;		/* space dimension (# of values) */

    SHORT	np_vert;	/* # of DOFs per vertex (must be 0 or 1) */
    SHORT	np_edge;	/* # of DOFs per edge */
    SHORT	np_face;	/* # of DOFs per face */
    SHORT	np_elem;	/* # of DOFs per element */
};

typedef struct HP_TYPE_ HP_TYPE;	/* defined later in hp.h */
typedef struct DOF_ {
    int		magic;		/* magic number */
    void        *cache_func;    /* cached function values */
    char	*name;		/* name of the dof */
    const char	*srcfile;	/* src file name in which DOF is created */
    GRID	*g;		/* the mesh */

    /* Note: exactly one of 'type' and 'hp' members is not NULL,
     * and 'hp' != NULL means a variable order DOF */
    DOF_TYPE	*type;		/* type of the DOF */
    HP_TYPE	*hp;		/* variable order DOF (h-p adaptivity) */

    /* data pointers */
    FLOAT	*data;		/* length = dim for DOF_CONSTANT,
				   0 for DOF_ANALYTIC, and
				     dim * (
					nvert * np_vert + nedge * np_edge +
					nface * np_face + nelem * np_elem ) */
    FLOAT	*data_vert;	/* pointer to vertex data */
    FLOAT	*data_edge;	/* pointer to edge data */
    FLOAT	*data_face;	/* pointer to face data */
    FLOAT	*data_elem;	/* pointer to volume data */

    BTYPE	*DB_masks;	/* Dirichlet boundary masks for each dof dim */

    /* (x,y,z) function attached to DOF */
    DOF_USER_FUNC userfunc;
    DOF_USER_FUNC_LAMBDA userfunc_lambda;

    /* some constants */
    size_t	count_vert;	/* data count per vertex */
    size_t	count_edge;	/* data count per edge */
    size_t	count_face;	/* data count per face */
    size_t	count_elem;	/* data count per volume */
    int		srcline;	/* line no in src file */

    SHORT	dim;		/* # of variables */
    SHORT	poly_order;	/* polynomial order (<0 ==> non polynomial) */
    BTYPE	DB_mask;	/* Dirichlet boundary mask */

    BOOLEAN     invariant;      /* TRUE ==> same values on all elems,
				 * default FALSE */
} DOF;

extern DOF_TYPE *DOF_DEFAULT;
extern DOF_TYPE DOF_CONSTANT_, DOF_ANALYTIC_;

/* special DOF_TYPEs */
#define DOF_TYPE_VOID		((DOF_TYPE *)NULL)
#define DOF_CONSTANT		(&DOF_CONSTANT_)
#define DOF_ANALYTIC		(&DOF_ANALYTIC_)
#define SpecialDofType(type)	((type) == DOF_ANALYTIC || \
				 (type) == DOF_CONSTANT)

#define DofTypeDim(dof)		(DofIsHP(dof) ? \
				    ((dof)->hp->info->dim) : \
				    (SpecialDofType((dof)->type) ? \
					1 : (dof)->type->dim))
#define DofTypeName(dof)	(DofIsHP(dof) ?			\
				 (dof)->hp->info->name : (dof)->type->name)
#define DofTypeOrder(dof, e)	(SpecialDofType((dof)->type) ? \
				 (dof)->poly_order : \
				 ((!DofIsHP(dof) ? (dof)->type->order : \
	(dof)->hp->info->types[(dof)->hp->elem_order[(e)->index]]->order)))

#define DofIsConstant(dof)	((dof)->type == DOF_CONSTANT)
#define DofIsAnalytic(dof)	((dof)->type == DOF_ANALYTIC)

#define DofVertexOrder(u, no)	(!DofIsHP(u) ? (u)->type->order :	\
	((u)->hp->info->flags & VERT_FLAG) ? \
		(u)->hp->info->types[(u)->hp->edge_order[no]]->order : 0)
#define DofEdgeOrder(u, no)	(!DofIsHP(u) ? (u)->type->order :	\
	((u)->hp->info->flags & EDGE_FLAG) ? \
		(u)->hp->info->types[(u)->hp->edge_order[no]]->order : 0
#define DofFaceOrder(u, no)	(!DofIsHP(u) ? (u)->type->order :	\
	((u)->hp->info->flags & FACE_FLAG) ? \
		(u)->hp->info->types[(u)->hp->face_order[no]]->order : 0)
#define DofElementOrder(u, no)	(!DofIsHP(u) ? (u)->type->order :	\
	(u)->hp->info->types[(u)->hp->elem_order[no]]->order)

/* special DOF user functions */
void phgDofUserFuncNoData_(FLOAT x, FLOAT y, FLOAT z, FLOAT *values);
void phgDofUserFuncNoAction_(FLOAT x, FLOAT y, FLOAT z, FLOAT *values);
void phgDofUserFuncInterpolation_(FLOAT x, FLOAT y, FLOAT z, FLOAT *values);
void phgDofUserFuncLambda_(FLOAT x, FLOAT y, FLOAT z, FLOAT *values);

void *phgDofUserFuncGetPointer_(int which);
#define DofUserFuncVoid		((DOF_USER_FUNC)NULL)
#define DofNoData		((DOF_USER_FUNC)phgDofUserFuncGetPointer_(0))
#define DofNoAction		((DOF_USER_FUNC)phgDofUserFuncGetPointer_(1))
#define DofInterpolation	((DOF_USER_FUNC)phgDofUserFuncGetPointer_(2))
#define DofLambdaFunction	((DOF_USER_FUNC)phgDofUserFuncGetPointer_(3))
#define SpecialDofUserFunction(func) \
	((func) == DofUserFuncVoid || (func) == DofNoData || \
	 (func) == DofNoAction || (func) == DofInterpolation || \
	 (func) == DofLambdaFunction)

#ifdef __cplusplus
}
#endif

/* Macro for testing an h-p DOF */
#if 0
# define DofIsHP(dof)	((dof)->type == NULL)
#else	/* 0|1 */
# define DofIsHP(dof)	((dof)->hp != NULL)
#endif	/* 0|1 */

/* Macro for getting FE_SPACE of a DOF */
#define DofFESpace(dof) \
	(DofIsHP(dof) ? (dof)->hp->info->fe_space : (dof)->type->fe_space)

/* Macro returning DOF dimension */
#define DofDim(dof)	(!DofIsHP(dof) ? (dof)->dim * (dof)->type->dim : \
				(dof)->dim * (dof)->hp->info->types[0]->dim)

/* Macros for accessing DOF data */
#define DofVertexData(dof, vert_no) \
	((dof)->data_vert + DofVertexDataIndex(dof, vert_no))
#define DofEdgeData(dof, edge_no) \
	((dof)->data_edge + DofEdgeDataIndex(dof, edge_no))
#define DofFaceData(dof, face_no) \
	((dof)->data_face + DofFaceDataIndex(dof, face_no))
#define DofElementData(dof, elem_no) \
	((dof)->data_elem + DofElementDataIndex(dof, elem_no))

#define DofVertexDataIndex(dof, vert_no) ( \
	(!DofIsHP(dof) ? (vert_no) * (dof)->count_vert : \
		((dof)->hp->info->flags & VERT_FLAG) ? \
			(dof)->dim * (dof)->hp->vert_index[vert_no] : 0))
#define DofEdgeDataIndex(dof, edge_no) ( \
	(!DofIsHP(dof) ? (edge_no) * (dof)->count_edge : \
		((dof)->hp->info->flags & EDGE_FLAG) ? \
			(dof)->dim * (dof)->hp->edge_index[edge_no] : 0))
#define DofFaceDataIndex(dof, face_no) ( \
	(!DofIsHP(dof) ? (face_no) * (dof)->count_face :  \
		((dof)->hp->info->flags & FACE_FLAG) ? \
			(dof)->dim * (dof)->hp->face_index[face_no] : 0))
#define DofElementDataIndex(dof, elem_no) ( \
	(!DofIsHP(dof) ? (elem_no) * (dof)->count_elem :  \
		(dof)->dim * (dof)->hp->elem_index[elem_no]))
#define DofData(dof) ((dof)->data)

/* Macros for getting data counts for the subgrid */
#define DofVertexDataCount_(dof, nvert) \
	(!DofIsHP(dof) ? (dof)->count_vert * (nvert) : \
		((dof)->hp->info->flags & VERT_FLAG) == 0 ? 0 : \
		 (dof)->dim * (dof)->hp->vert_index[nvert])
#define DofVertexDataCount(dof) \
	DofVertexDataCount_(dof, (dof)->g->nvert)

#define DofEdgeDataCount_(dof, nedge) \
	(!DofIsHP(dof) ? (dof)->count_edge * (nedge) : \
		((dof)->hp->info->flags & EDGE_FLAG) == 0 ? 0 : \
		 (dof)->dim * (dof)->hp->edge_index[nedge])
#define DofEdgeDataCount(dof) \
	DofEdgeDataCount_(dof, (dof)->g->nedge)

#define DofFaceDataCount_(dof, nface) \
	(!DofIsHP(dof) ? (dof)->count_face * (nface) : \
		((dof)->hp->info->flags & FACE_FLAG) == 0 ? 0 : \
		 (dof)->dim * (dof)->hp->face_index[nface])
#define DofFaceDataCount(dof) \
	DofFaceDataCount_(dof, (dof)->g->nface)

#define DofElementDataCount_(dof, nelem) \
	(!DofIsHP(dof) ? (dof)->count_elem * (nelem) : \
		((dof)->hp->info->flags & ELEM_FLAG) == 0 ? 0 : \
		 (dof)->dim * (dof)->hp->elem_index[nelem])
#define DofElementDataCount(dof) \
	DofElementDataCount_(dof, (dof)->g->nelem)

#define DofDataCount(dof)	\
	(DofVertexDataCount(dof) + DofEdgeDataCount(dof) + \
	 DofFaceDataCount(dof)   + DofElementDataCount(dof))
#define DofGetDataCount(dof)	DofDataCount(dof)

/* Macros for getting data counts for the whole grid */
#define DofVertexDataCountGlobal(dof) \
	(DofIsHP(dof) ? (dof)->dim * (dof)->hp->vert_data_count_global : \
	 ((dof)->count_vert * ((dof)->g->period == NULL ? \
					(dof)->g->nvert_global : \
					(dof)->g->period->nvert_global)))
#define DofEdgeDataCountGlobal(dof) \
	(!DofIsHP(dof) ? (dof)->count_edge * (dof)->g->nedge_global : \
			 (dof)->dim * (dof)->hp->edge_data_count_global)
#define DofFaceDataCountGlobal(dof) \
	(!DofIsHP(dof) ? (dof)->count_face * (dof)->g->nface_global : \
			 (dof)->dim * (dof)->hp->face_data_count_global)
#define DofElementDataCountGlobal(dof) \
	(!DofIsHP(dof) ? (dof)->count_elem * (dof)->g->nelem_global : \
			 (dof)->dim * (dof)->hp->elem_data_count_global)
#define DofDataCountGlobal(dof)	\
	(DofVertexDataCountGlobal(dof) + DofEdgeDataCountGlobal(dof) + \
	 DofFaceDataCountGlobal(dof)   + DofElementDataCountGlobal(dof))
#define DofGetDataCountGlobal(dof)	DofDataCountGlobal(dof)

/* Macros for getting basis/data counts on a given vertex/edge/face/element */
#define DofDataCountOnVertex(dof, vert_no) \
	(!DofIsHP(dof) ? (dof)->count_vert : \
		((dof)->hp->info->flags & VERT_FLAG) == 0 ? 0 : (dof)->dim * \
			((dof)->hp->vert_index[vert_no + 1] - \
			 (dof)->hp->vert_index[vert_no]))
#define DofDataCountOnEdge(dof, edge_no) \
	(!DofIsHP(dof) ? (dof)->count_edge : \
		((dof)->hp->info->flags & EDGE_FLAG) == 0 ? 0 : (dof)->dim * \
			((dof)->hp->edge_index[edge_no + 1] - \
			 (dof)->hp->edge_index[edge_no]))
#define DofDataCountOnFace(dof, face_no) \
	(!DofIsHP(dof) ? (dof)->count_face : \
		((dof)->hp->info->flags & FACE_FLAG) == 0 ? 0 : (dof)->dim * \
			((dof)->hp->face_index[face_no + 1] - \
			 (dof)->hp->face_index[face_no]))
#define DofDataCountOnElement(dof, elem_no) \
	(!DofIsHP(dof) ? (dof)->count_elem : (dof)->dim * \
		((dof)->hp->elem_index[elem_no + 1] - \
		 (dof)->hp->elem_index[elem_no]))

/* Macro returning type and # basis functions at a vert/edge/face/elem */
#define DofNBasOnVertex(dof, vert_no) (!DofIsHP(dof) ? (dof)->type->np_vert : \
	HPNBasOnVertex((dof)->hp, vert_no))
#define DofNBasOnEdge(dof, edge_no) (!DofIsHP(dof) ? (dof)->type->np_edge : \
	HPNBasOnEdge((dof)->hp, edge_no))
#define DofNBasOnFace(dof, face_no) (!DofIsHP(dof) ? (dof)->type->np_face : \
	HPNBasOnFace((dof)->hp, face_no))
#define DofNBasOnElement(dof, elem_no) (!DofIsHP(dof) ? (dof)->type->np_elem : \
	HPNBasOnElement((dof)->hp, elem_no))
#define DofNBas(dof, e) \
	(!DofIsHP(dof) ? (dof)->type->nbas : HPNBas((dof)->hp, e))
#define DofGetNBas(dof, e) DofNBas(dof, e)	/* for legacy programs */

#define DofTypeOnVertex(dof, vert_no) \
	(!DofIsHP(dof) ? (dof)->type : HPTypeOnVertex((dof)->hp, vert_no))
#define DofTypeOnEdge(dof, edge_no) \
	(!DofIsHP(dof) ? (dof)->type : HPTypeOnEdge((dof)->hp, edge_no))
#define DofTypeOnFace(dof, face_no) \
	(!DofIsHP(dof) ? (dof)->type : HPTypeOnFace((dof)->hp, face_no))
#define DofTypeOnElement(dof, elem_no) \
	(!DofIsHP(dof) ? (dof)->type : HPTypeOnElement((dof)->hp, elem_no))

/* structure for exchanging neighbour data */
typedef struct {
    DOF		*dof;
    struct MAP_	*map;
    FLOAT	*data;	/* remote data, stored in the order as return by
			 * phgDofGetBasesOnFace for each face as:
			 * 	data[m][n][dof->dim]
			 * where:
			 *	m is total number of shared faces
			 *	n = 3 * (np_vert + np_edge) + np_face
			 *	    (or 1 for DOF_DG0 or DOF_P0).
			 *
			 * Note: for non hp DOF 'n' == nbas, for hp DOF
			 *	 'n' varies and 'n'[i] == loc[i + 1] - loc[i]
			 *	 for shared face i.
			 */
    INT		*index;	/* global vector indices stored as: ind[m][n] */
    INT		*loc;	/* loc[i] and loc[i + 1] are respectively start and end
			 * position in 'data' and 'index' for shared face i. */
    BYTE	*order;	/* order[i] is the order on remote neighbour i */

    /* The following variables are for caching data for local neighbours */
    struct {
	ELEMENT	*peer;	/* local neighbour */
	SHORT	*plist;	/* buffer for storing list of bases on the face
			   for the local neighbour */
	int	face;	/* the face no (with respect to 'peer'). */
    } *cache;

    /* common constants */
    int		nbas;	/* == length of plist */
    int		dof_no;	/* map[dof_no] == dof */
} NEIGHBOUR_DATA;

typedef struct {
    ELEMENT	**Vmap, **Emap, **Fmap;	/* vertex/edge/face to element maps */
    BYTE	*Vind, *Eind, *Find;	/* local indices in the element */
} VEF_MAP;

typedef struct {
    ELEMENT     **elist;                        /* all element in vef map*/
    ELEMENT	***Vmap, ***Emap, ***Fmap;	/* vertex/edge/face to element patch start */
    BYTE	**Vind, **Eind, **Find;	        /* local indices in the element */
    int         *Vsize, *Esize, *Fsize;         /* size of vertex/edge/face to element patch */
    BYTE        *idxlist;           		/* all local indices in the element in vef map */
} VEF_MAP2;

/* projection modes onto a face */
typedef enum {DOF_PROJ_NONE, DOF_PROJ_DOT, DOF_PROJ_CROSS} DOF_PROJ;
extern const char *_phg_dof_proj_name[];

typedef struct {
    DOF_TYPE	**types_array;	/* DOF_TYPEs array (Pn, DGn, etc.) */
    const char	*desc;		/* a short description of the DOF_TYPEs */
} DOF_TYPES_INFO;

extern DOF_TYPES_INFO phgDofTypesList[];

/* function prototypes */

#ifdef __cplusplus
extern "C" {
#endif

void phgDofTypesInit(void);
DOF_TYPE *phgDofTypeVector(DOF_TYPE *base_type, int dim, const char *name);
void phgDofTypeFree(DOF_TYPE **ptype);

const char *phgDofFESpaceName(FE_SPACE fe_space);

void phgDofRefId(DOF_TYPE *type);
void phgDofUnrefId(DOF_TYPE *type);
size_t phgDofRefCount(SHORT id);

DOF *phgDofNew_(GRID *g, DOF_TYPE *type, HP_TYPE *hp, SHORT dim,
		const char *name, DOF_USER_FUNC userfunc,
		const char *srcfile, int srcline);
#define phgDofNew(g, type, dim, name, userfunc) \
	phgDofNew_(g, type, NULL, dim, name, userfunc, __FILE__, __LINE__)
#define phgHPDofNew(g, hp, dim, name, userfunc) \
	phgDofNew_(g, NULL, hp,   dim, name, userfunc, __FILE__, __LINE__)
#define phgDofClone(dof, name, userfunc) \
	phgDofNew_((dof)->g, (dof)->type, (dof)->hp, (dof)->dim, \
		   name, userfunc, __FILE__, __LINE__)
BTYPE phgDofSetDirichletBoundaryMask(DOF *u, BTYPE mask);
void phgDofSetDirichletBoundaryMasks(DOF *u, const BTYPE masks[]);
BOOLEAN phgDofIsValid(GRID *g, DOF *dof);
DOF *phgDofIndex2Pointer(GRID *g, int id);
int  phgDofPointer2Index(DOF *dof);

void phgDofSetPolyOrder(DOF *u, SHORT order);
void phgDofSetFunction(DOF *u, DOF_USER_FUNC func);
void phgDofSetLambdaFunction(DOF *u, DOF_USER_FUNC_LAMBDA func);

void phgDofClearCache(GRID *g, DOF *dof, DOF_TYPE *type, HP_TYPE *hp,
		      BOOLEAN final);
void phgDofFree(DOF **dof);
void phgDofMapEdgeData(const DOF_TYPE *type, const ELEMENT *e, int edge_no,
			int *M);
void phgDofMapFaceData(const DOF_TYPE *type, const ELEMENT *e, int face_no,
			int *M);
void phgDofMapElementData(const DOF_TYPE *type, const ELEMENT *e, int *M);
int  phgDofGetBasesOnFace(DOF *dof, ELEMENT *e, int face_no, SHORT bases[]);
void phgDofDump(DOF *dof);
void phgDofInterpC2FGeneric(DOF *dof, ELEMENT *e, FLOAT **parent_data,
		FLOAT **children_data);
void phgDofInterpF2CGeneric(DOF *dof, ELEMENT *e, FLOAT **parent_data,
		FLOAT **children_data);
void phgDofInitFuncPoint(DOF *dof, ELEMENT *e, GTYPE type, int index,
		DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
		const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues);

VEF_MAP *phgDofSetupVEFMap_(GRID *g, DOF *dof, int flags, FORALL_TYPE for_what);
#define phgDofSetupVEFMap(g, dof, flags) \
	phgDofSetupVEFMap_(g, dof, flags, FOR_ALL)
void phgDofFreeVEFMap(VEF_MAP **vef_ptr);

VEF_MAP2 *phgDofSetupVEFMap2_(GRID *g, DOF *dof, int flags, int mask,
		const BTYPE *types_vert, const BTYPE *types_edge,
		const BTYPE *types_face, FORALL_TYPE for_what);
#define phgDofSetupVEFMap2(g, dof, flags) \
	phgDofSetupVEFMap2_(g, dof, flags, -1, NULL, NULL, NULL, FOR_OWNED)
void phgDofFreeVEFMap2(VEF_MAP2 **vef_ptr);

void phgDofSetDataByValue(DOF *dof, FLOAT value);
void phgDofSetDataByValues(DOF *dof, const FLOAT array_of_values[]);
void phgDofSetDataByValuesV(DOF *dof, FLOAT v0, ...);
void phgDofSetDataByFunction__(DOF *dof, DOF_USER_FUNC userfunc,
	DOF_USER_FUNC_LAMBDA userfunc_lambda, FLOAT *funcvalues,
	BTYPE bdry_mask, FORALL_TYPE for_what);
#define phgDofSetDataByFunction_(dof, userfunc, userfunc_lambda, funcvalues) \
	phgDofSetDataByFunction__(dof, userfunc, userfunc_lambda, funcvalues, \
	~0, FOR_ALL)
#define phgDofSetDataByFunction(dof, f) \
	phgDofSetDataByFunction_(dof, f, NULL, NULL)
#define phgDofSetDataByLambdaFunction(dof, f) \
	phgDofSetDataByFunction_(dof, NULL, f, NULL)
#define phgDofSetBdryDataByFunction(dof, f, bdry_mask) \
	phgDofSetDataByFunction__(dof, f, NULL, NULL, bdry_mask, FOR_OWNED)
#define phgDofSetBdryDataByLambdaFunction(dof, f, bdry_mask) \
	phgDofSetDataByFunction__(dof, NULL, f, NULL, bdry_mask, FOR_OWNED)


/* functions called by distribute.c */
void phgDofRedistribute(GRID *g, BOOLEAN flag);

/* functions called by refine.c */
void phgDofInterpolate(GRID *g, BOOLEAN flag);

/* Utility functions */
INT phgDofMapE2D_(DOF *u, ELEMENT *e, int index, BOOLEAN global_ordering);
#define phgDofMapE2D(dof, e, index) phgDofMapE2D_(dof, e, index, TRUE)

#define BasisOrder(u, e, i) (!DofIsHP(u) ? (u)->type->order : \
	(u)->hp->info->types[(u)->hp->elem_order[e->index]]->order)

BTYPE phgDofGetElementBasisInfo_(DOF *u, ELEMENT *e, int bas_no,
				 GTYPE *gtype, int *gindex, int *bindex,
				 FLOAT xyz[], FLOAT lambda[]);
#define phgDofGetElementBasisInfo(u, e, bas_no, gtype, gindex, bindex) \
	phgDofGetElementBasisInfo_(u, e, bas_no, gtype, gindex, bindex, \
				   NULL, NULL)
#define phgDofGetElementBoundaryType(u, e, index) \
	phgDofGetElementBasisInfo(u, e, (index) / (u)->dim, NULL, NULL, NULL)
FLOAT *phgDofGetElementCoordinates(DOF *u, ELEMENT *e, int index);

void phgDofGetElementDatas(DOF *dof, SIMPLEX *e, FLOAT *value);
void phgDofSetElementDatas(DOF *dof, SIMPLEX *e, FLOAT *value);

NEIGHBOUR_DATA *phgDofInitNeighbourData(DOF *dof, struct MAP_ *map);
void phgDofReleaseNeighbourData(NEIGHBOUR_DATA **nd_ptr);
FLOAT *phgDofNeighbourData(NEIGHBOUR_DATA *nd, ELEMENT *e, int face_no,
			int bas_no, INT *gindex);
int phgDofNeighbourNBas(NEIGHBOUR_DATA *nd, ELEMENT *e, int face_no,
			DOF_TYPE **peer_type);

#ifdef __cplusplus
}
#endif

#define PHG_DOF_H
#endif
