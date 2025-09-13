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

/* $Id: dof.c,v 1.393 2022/09/02 01:40:16 zlb Exp $ */

#include "phg.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <limits.h>	/* INT_MAX */

DOF_TYPE *DOF_DEFAULT = DOF_P2;

/* Special DOFs */

DOF_TYPE DOF_CONSTANT_ = {DofReserved, "Constant", NULL, NULL, &DOF_CONSTANT_,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, FE_Cinfty,
    FALSE,			/* is_nodal */
    TRUE, FALSE, -1,		/* invariant, free_after_use, id */
    0, 0, 0, 0, 1,		/* nbas, order, D_order, cont, dim */
    0, 0, 0, 0
};				/* np_vert, np_edge, np_face, np_elem */

DOF_TYPE DOF_ANALYTIC_ = {DofReserved, "Analytic", NULL, NULL, NULL,
    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, FE_None,
    FALSE,			/* is_nodal */
    FALSE, FALSE, -1,		/* invariant, free_after_use, id */
    0, -1, 0, -1, 1,		/* nbas, order, D_order, cont, dim */
    0, 0, 0, 0
};				/* np_vert, np_edge, np_face, np_elem */

DOF_TYPES_INFO phgDofTypesList[] = {
    {DOF_Pn,	"H^1 nodal (uniform):             "},
    {DOF_MLn,	"H^1 nodal mass-lumped (new):     "},
    {DOF_MLon,	"H^1 nodal mass-lumped (old):     "},
    {DOF_MLDGn,	"L2 nodal mass-lumped:            "},
#if TEST_TYPES
    {DOF_TSTn,	"H^1 nodal (non uniform):         "},
#endif	/* TEST_TYPES */
    {DOF_NDn,	"1st kind H(curl) (by Cui Tao):   "},
    {DOF_HFEBn,	"Hierarchical H(curl) (by CJQ):   "},
    {DOF_HBn,	"Hierarchical H^1:                "},
    {DOF_HCn,	"Hierarchical H(curl):            "},
    {DOF_HDn,	"Hierarchical H(div):             "},
    {DOF_DGn,	"Discontinuous (L^2):             "},
    {DOF_SPTn,	"Specht (H^2 nonconforming):      "},
    {NULL,	NULL}
};


void
phgGetFaceVertices_(const ELEMENT *e, int face_no,
		    int *v0, int *v1, int *v2, int *v3,
		    int *map)
/* Get face vertices in the global vertices order,
 *   if map is not NULL, record the inverse map.
 * Tetra: v3 not uesd.
 * Hexa : 4 verts maped in the order that starting
 *   with min, then adjacent smaller, then goes along
 *   the edges.
 * 
 * */
{
    if (e->elem_type == ET_TETRA) {
	int _u0, _u1, _u2, _tmp;			
	*v0 = GetFaceVertex(face_no, 0);		
	*v1 = GetFaceVertex(face_no, 1);		
	*v2 = GetFaceVertex(face_no, 2);		
	_u0 = MapVertex(e, *v0);			
	_u1 = MapVertex(e, *v1);			
	_u2 = MapVertex(e, *v2);			
	if (_u0 > _u1) {				
	    _tmp = *v0, *v0 = *v1; *v1 = _tmp;		
	    _tmp = _u0, _u0 = _u1; _u1 = _tmp;	
	}						
	if (_u0 > _u2) {				
	    _tmp = *v0, *v0 = *v2; *v2 = _tmp;		
	    _u2 = _u0;				
	}						
	if (_u1 > _u2) {				
	    _tmp = *v1, *v1 = *v2; *v2 = _tmp;		
	}
	*v3 = -1;		/* v3 not used for triangle */
	if (map != NULL)
	    map = NULL; 	/* To be added for triangle */
    }
#if PHG_GRID_TYPE != GRID_TYPE_TET
    else if (e->elem_type == ET_HEXA) {
	int i, v[4], u[4], p[4], u0;
	const int before[] = {3, 0, 1, 2};
	const int after[] = {1, 2, 3, 0};

	for (i = 0; i < 4; i++) {
	    v[i] = GetFaceVertex(face_no, i);
	    u[i] = MapVertex(e, v[i]);
	} 
    
	/* min */
	u0 = u[0];
	p[0] = 0;
	for (i = 0; i < 4; i++) 
	    if (u0 > u[i]) {
		u0 = u[i];
		p[0] = i;
	    }

	/* direction */
	if (u[before[p[0]]] < u[after[p[0]]]) {
	    p[1] = before[p[0]];
	    p[2] = before[p[1]];
	    p[3] = before[p[2]];
	} else {
	    p[1] = after[p[0]];
	    p[2] = after[p[1]];
	    p[3] = after[p[2]];
	}
	
	*v0 = v[p[0]]; 
	*v1 = v[p[1]]; 
	*v2 = v[p[2]]; 
	*v3 = v[p[3]]; 

	phgInfo(3, "  u %3d, %3d, %3d, %3d\n", u[0], u[1], u[2], u[3]);
	phgInfo(3, "  v %3d, %3d, %3d, %3d\n", v[0], v[1], v[2], v[3]);
	phgInfo(3, "  p %3d, %3d, %3d, %3d\n", p[0], p[1], p[2], p[3]);

	if (map != NULL) {
	    for (i = 0; i < 4; i++)
		map[i] = p[i];
	}
    }
#endif	/* PHG_GRID_TYPE != GRID_TYPE_TET */
    else {
	phgError(1, "%s:%d: %s unimplemented.\n", __FILE__, __LINE__, __func__);
    }

    return;
}


void
phgGetElementVertices_(const ELEMENT *e, int index,
		       int *v0, int *v1, int *v2, int *v3,
		       int *v4, int *v5, int *v6, int *v7)
/* Get element vertices in the global vertices order.
 * Tetra: v[4-7] not uesd. */
{
    if (e->elem_type == ET_TETRA) {
	int tmp[4];					
	assert((index) == 0);			
	tmp[MapVertex(e, 0)] = 0;			
	tmp[MapVertex(e, 1)] = 1;			
	tmp[MapVertex(e, 2)] = 2;			
	tmp[MapVertex(e, 3)] = 3;			
	*v0 = tmp[0];				
	*v1 = tmp[1];				
	*v2 = tmp[2];				
	*v3 = tmp[3];	
	*v4 = 0;
	*v5 = 0;
	*v6 = 0;
	*v7 = 0;
    } else if (e->elem_type == ET_HEXA) {
	phgError(1, "%s: %s unimplemented.\n", __FILE__, __FUNCTION__);
    } else {
	phgError(1, "%s: %s unimplemented.\n", __FILE__, __FUNCTION__);
    }
}

void
phgDofTypesInit(void)
/* This function calls Initialize() function of all pre-defined DOF_TYPEs */
{
    int i;
    DOF_TYPE **t;

    for (i = 0; phgDofTypesList[i].types_array != NULL; i++) {
	t = phgDofTypesList[i].types_array;
	for (; *t == NULL; t++);	/* skip leading NULLs */
	for (; *t != NULL; t++) {
	    if ((*t)->Initialize == NULL)
		continue;
	    (*t)->Initialize(*t);
	    (*t)->Initialize = NULL;
	}
    }
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

    phgFree((void *)type->name);
    phgFree(type->points);
    phgFree(type->orders);
    phgFree(type);
}

static void
vtype_init(DOF *dof, ELEMENT *e, GTYPE gtype, int index,
	    DOF_USER_FUNC userfunc,
            DOF_USER_FUNC_LAMBDA userfunc_lambda,
            const FLOAT *funcvalues, FLOAT *dofvalues,
            FLOAT **pdofvalues)
{
    DOF_TYPE *type = dof->type;
    int d;
    void *hp = dof->hp;		/* backup dof->hp */

    assert(type != NULL && type->base_type != NULL &&
	   type->base_type->InitFunc != NULL);
    d = type->dim / type->base_type->dim;

    dof->dim *= d;
    dof->hp = NULL;
    dof->type = type->base_type;
    dof->type->InitFunc(dof, e, gtype, index, userfunc, userfunc_lambda,
			funcvalues, dofvalues, pdofvalues);
    dof->type = type;		/* restore dof->type */
    dof->hp = hp;			/* restore dof->hp */
    dof->dim /= d;
}

static const FLOAT *
vtype_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT *buffer = NULL;
    static int buffer_size = 0;
#if USE_OMP
#pragma omp threadprivate(buffer, buffer_size)
#endif  /* USE_OMP */
    const FLOAT *p0;
    FLOAT *p;
    DOF_TYPE *type = dof->type;
    int i, d, dim0;

    if (no1 <= 0)
	no1 = type->nbas;

    if (no1 <= no0)
	return NULL;
    assert(type != NULL && type->base_type != NULL);
    dim0 = type->base_type->dim;
    d =  type->dim / dim0;
    dof->type = type->base_type;
    p0 = type->base_type->BasFuncs(dof, e, no0 / d, (no1 + d - 1) / d,
					lambda);
    dof->type = type;
    if (d == 1)
	return p0;

    i = type->nbas * type->dim;
    if (buffer_size < i) {
	BOOLEAN flag = (buffer == NULL);
	phgFree(buffer);
	buffer = phgAlloc((buffer_size = i) * sizeof(*buffer));
	if (flag)
	    FreeAtExitNoCheck(buffer);
    }

    bzero(buffer, (no1 - no0) * type->dim * sizeof(*buffer));
    for (i = no0, p = buffer; i < no1; i++, p += type->dim)
	memcpy(p + (i % d) * dim0, p0 + (i / d) * dim0, dim0 * sizeof(*p));

    return buffer;
}

static const FLOAT *
vtype_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    static FLOAT *buffer = NULL;
    static int buffer_size = 0;
#if USE_OMP
#pragma omp threadprivate(buffer, buffer_size)
#endif  /* USE_OMP */
    const FLOAT *p0;
    FLOAT *p;
    DOF_TYPE *type = dof->type;
    int i, d, dim0;

    if (no1 <= 0)
	no1 = type->nbas;

    if (no1 <= no0)
	return NULL;
    assert(type != NULL && type->base_type != NULL);
    dim0 = type->base_type->dim;
    d =  type->dim / dim0;
    dof->type = type->base_type;
    p0 = type->base_type->BasGrads(dof, e, no0 / d, (no1 + d - 1) / d,
					lambda);
    dof->type = type;

    if (d == 1)
	return p0;

    i = type->nbas * type->dim * (Dim + 1);
    if (buffer_size < i) {
	BOOLEAN flag = (buffer == NULL);
	phgFree(buffer);
	buffer = phgAlloc((buffer_size = i) * sizeof(*buffer));
	if (flag)
	    FreeAtExitNoCheck(buffer);
    }

    bzero(buffer, (Dim + 1) * (no1 - no0) * type->dim * sizeof(*buffer));
    for (i = no0, p = buffer; i < no1; i++, p += type->dim * (Dim + 1))
	memcpy(p + (i % d) * dim0 * (Dim + 1),
	       p0 + (i / d) * dim0 * (Dim + 1), dim0 * (Dim + 1) * sizeof(*p));

    return buffer;
}

static void
vtype_map(const DOF_TYPE *type, const ELEMENT *e, int no, int *M)
/* returns ordering of local face or element data */
{
    assert(type != NULL && type->base_type != NULL &&
	   type->base_type->DofMap != NULL);
    if (type->dim == type->base_type->dim) {
	type->base_type->DofMap(type->base_type, e, no, M);
    }
    else {
	int n = no < 0 ? type->base_type->np_elem : type->base_type->np_face;
	int M0[n];
	int i, j, d = type->dim / type->base_type->dim;
	type->base_type->DofMap(type->base_type, e, no, M0);
	for (i = 0; i < n; i++)
	    for (j = 0; j < d; j++)
		M[i * d + j] = d * M0[i] + j;
    }
}

DOF_TYPE *
phgDofTypeVector(DOF_TYPE *base_type, int dim, const char *name)
/* returns a vector DOF_TYPE of dimension dim * "base_type"->dim
 * The optional argument "name" gives the name of the new type */
{
    DOF_TYPE temp = {
	DofReserved,			/* reserved pointers */
	NULL,				/* name */
	NULL,				/* points */
	NULL,				/* orders */
	base_type->grad_type,		/* grad_type */
	base_type,			/* base_type */
	NULL,				/* hp_info */
	phgDofInterpC2FGeneric,		/* InterpC2F */
	phgDofInterpF2CGeneric,		/* InterpF2C */
	vtype_init,			/* InitFunc */
	vtype_bas,			/* BasFuncs */
	vtype_grad,			/* BasGrads */
	base_type->DofMap == NULL ? NULL : vtype_map, /* DofMap */
	base_type->fe_space,		/* fe_space */
	base_type->is_nodal,		/* is_nodal */
	base_type->invariant,		/* invariant */
	FALSE,				/* free_after_use */
	-1,				/* id */
	base_type->nbas * dim,		/* nbas */
	base_type->order,		/* order */
	base_type->D_order,		/* D_order */
	base_type->continuity,		/* continuity */
	base_type->dim * dim,		/* dim */
	base_type->np_vert * dim,	/* np_vert */
	base_type->np_edge * dim,	/* np_edge */
	base_type->np_face * dim,	/* np_face */
	base_type->np_elem * dim	/* np_elem */
    };
    DOF_TYPE *type = phgAlloc(sizeof(*type));
    char *s;

    assert(dim >= 1);
    memcpy(type, &temp, sizeof(temp));

    if (name != NULL)
	s = strdup(name);
    else
	sprintf(s = phgAlloc(strlen(base_type->name) + 1 + 1),
		"%sv", base_type->name);
    type->name = s;

    if (base_type->orders != NULL) {
	int i, n;
	n = type->np_vert + type->np_edge + type->np_face + type->np_elem;
	type->orders = phgAlloc(n * sizeof(*type->orders));
	for (i = 0; i < n; i++)
	    type->orders[i] = base_type->orders[i / dim];
    }

    return type;
}

const char *
phgDofFESpaceName(FE_SPACE fe_space)
{
    switch (fe_space) {
	case FE_None:
	    return "none";
	case FE_L2:
	    return "L_2";
	case FE_Hcurl:
	    return "H(curl)";
	case FE_Hdiv:
	    return "H(div)";
	case FE_H1:
	    return "H^1";
	case FE_Cinfty:
	    return "C(infty)";
    }

    return "unknown";
}

/* Dummy USER_FUNCs */

void
phgDofUserFuncNoData_(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    Unused(x);
    Unused(y);
    Unused(z);
    Unused(values);
    phgError(1, "%s:%d, calling dummy DOF function.\n", __FILE__, __LINE__);
    return;
}

void
phgDofUserFuncNoAction_(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    Unused(x);
    Unused(y);
    Unused(z);
    Unused(values);
    phgError(1, "%s:%d, calling dummy DOF function.\n", __FILE__, __LINE__);
    return;
}

void
phgDofUserFuncInterpolation_(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    Unused(x);
    Unused(y);
    Unused(z);
    Unused(values);
    phgError(1, "%s:%d, calling dummy DOF function.\n", __FILE__, __LINE__);
    return;
}

void
phgDofUserFuncLambda_(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    Unused(x);
    Unused(y);
    Unused(z);
    Unused(values);
    phgError(1, "%s:%d, calling dummy DOF function.\n", __FILE__, __LINE__);
    return;
}

void *
phgDofUserFuncGetPointer_(int which)
/* returns a pointer to one of phgDofUserFuncNoData_, phgDofUserFuncNoAction_,
 * phgDofUserFuncInterpolation_ and phgDofUserFuncLambda_, according to which.
 *
 * Note: this is to make the tests like 'userfunc == phgDofUserFuncNoAction_'
 * work when using libphg.dll under MinGW. */
{
    switch (which) {
	case 0:	return phgDofUserFuncNoData_;
	case 1: return phgDofUserFuncNoAction_;
	case 2: return phgDofUserFuncInterpolation_;
	case 3: return phgDofUserFuncLambda_;
	default: phgError(1, "%s: invalid argument %d\n", which);
    }

    return NULL;
}

/*--------------- list of referenced DOF_TYPEs and HP_INFOs ----------------*/

/* list of referenced DOF_TYPES */
static struct {
    /* DOF_TYPE pointer */
    DOF_TYPE	*type;
    size_t	refcount;	/* Note: move to DOF_TYPE? */
} *type_list = NULL;
static size_t type_list_count = 0, type_list_allocated = 0;

void
phgDofRefId(DOF_TYPE *type)
/* updates id of 'type' */
{
    SHORT id;

    if ((id = type->id) >= 0) {
	assert(id < type_list_count);
	type_list[id].refcount++;
	return;
    }
    
    /* add a new entry */
    if (type_list_count >= type_list_allocated) {
	size_t len = (type_list_allocated += 8) * sizeof(*type_list);
	type_list = phgRealloc_(type_list, len,
				(type_list_allocated - 8) * sizeof(*type_list));
    }

    type_list[type_list_count].type = type;
    type->id = type_list_count;

    type_list[type_list_count++].refcount = 1;

    return;
}

void
phgDofUnrefId(DOF_TYPE *type)
/* unreferences the id of 'type' */
{
    SHORT id = type->id;
    int i;

    if (id < 0)
	return;

    assert(id < type_list_count);

    if (--type_list[id].refcount > 0)
	return;

    phgDofClearCache(NULL, NULL, type, NULL, TRUE);

    if (id < --type_list_count) {
	/* remove the type from type_list[] */
	memmove(type_list + id, type_list + id + 1,
				(type_list_count - id) * sizeof(*type_list));
	for (i = id; i < type_list_count; i++) {
	    type_list[i].type->id = i;
	}
    }

    if (type->free_after_use)
	phgDofTypeFree(&type);
    else
	type->id = -1;

    if (type_list_count == 0) {
	phgFree(type_list);
	type_list = NULL;
	type_list_allocated = 0;
    }

    return;
}

size_t
phgDofRefCount(SHORT id)
{
    return id < type_list_count ? type_list[id].refcount : 0;
}

/*--------------------------------------------------------------------------*/

DOF *
phgDofNew_(GRID *g, DOF_TYPE *type, HP_TYPE *hp, SHORT dim, const char *name,
	   DOF_USER_FUNC userfunc, const char *srcfile, int srcline)
{
    DOF *dof;
    size_t len;

    FunctionEntry;

    if ((type == NULL) == (hp == NULL))
	phgError(1, "%s: one of 'type' and 'hp' should be NULL (%s:%d).\n",
			__func__, srcfile, srcline);

    if (phgMasterSlave && phgNProcs > 1) {
	/* Master/slave mode is not yet implemented */
	phgError(1, "phgDofNew does not allow master/slave mode\n");
    }

    if (type != NULL)
	phgDofRefId(type);
    else
	type = hp->info->types[hp->max_order];

    if (!SpecialDofType(type) && !phgGridIsValid(g))
	phgError(1, "%s:%d, invalid grid.\n",
			srcfile == NULL ? "" : (void *)srcfile, srcline);

    dof = phgCalloc(1, sizeof(*dof));
    dof->magic = MAGIC_DOF;
    if (name == NULL)
	name = "noname";
    phgInfo(2, "creating DOF \"%s\" (call from %s:%d).\n",
		    			name, srcfile, srcline);
    dof->name = phgAlloc(strlen(name) + 1);
    strcpy(dof->name, name);
    if (hp != NULL) {
	dof->type = NULL;
	dof->hp = hp;
	hp->refcount++;
    }
    else {
	dof->hp = NULL;
	dof->type = type;
    }
    dof->userfunc = userfunc;
    dof->userfunc_lambda = NULL;
    dof->invariant = FALSE;
    dof->DB_mask = DIRICHLET;
    dof->DB_masks = NULL;
    if (type == DOF_CONSTANT && userfunc != DofNoAction
	&& userfunc != DofNoData) {
	phgPrintf("changing userfunc of constant DOF \"%s\" to DofNoAction\n",
		  dof->name);
	dof->userfunc = DofNoAction;
    }
    if (userfunc == DofInterpolation &&
	type != NULL && type->InterpC2F == NULL) {
	phgPrintf("DOF \"%s\": no interpolation function for type \"%s\", "
		  "changing userfunc to DofNoAction\n", dof->name,
		  type->name);
	dof->userfunc = DofNoAction;
    }
    dof->g = g;
    if ((dof->dim = dim) <= 0)
	phgError(1, "invalid value for dim (%d)", dim);
    dof->srcfile = srcfile;
    dof->srcline = srcline;

    if (type == DOF_CONSTANT)
	dof->poly_order = 0;
    else if (SpecialDofType(type) || type == NULL)
	dof->poly_order = -1;
    else
	dof->poly_order = type->order;

    if (!SpecialDofType(type)) {
	if (type->np_vert > 0 && !(g->flags & VERT_FLAG))
	    phgError(1, "%s (%s:%d) elements has vertex DOFs "
		     "but VERT_FLAG is FALSE.\n", dof->type->name,
			srcfile, srcline);
	if (type->np_edge > 0 && !(g->flags & EDGE_FLAG))
	    phgError(1, "%s (%s:%d) elements require edge indices "
		     "but EDGE_FLAG is FALSE.\n", dof->type->name,
			srcfile, srcline);
	if (type->np_face > 0 && !(g->flags & FACE_FLAG))
	    phgError(1, "%s (%s:%d) elements require face indices "
		     "but FACE_FLAG is FALSE.\n", dof->type->name,
			srcfile, srcline);
	if (type->np_elem > 0 && !(g->flags & ELEM_FLAG))
	    phgError(1, "%s (%s:%d) elements require element indices "
		     "but ELEM_FLAG is FALSE.\n", dof->type->name,
			srcfile, srcline);
    }

    /* attach the new DOF to g */
    if (g != NULL) {
	len = 0;
	if (g->dof != NULL) {
	    DOF **p;
	    for (p = g->dof; *p != NULL; p++, len++);
	}
	g->dof = phgRealloc_(g->dof, (len + 2) * sizeof(*(g->dof)),
				 len * sizeof(*(g->dof)));
	g->dof[len] = dof;
	g->dof[len + 1] = NULL;
    }

    /* attach dof to hp->dof */
    if (hp != NULL) {
	len = 0;
	if (hp->dof != NULL) {
	   DOF **p;
	   for (p = hp->dof; *p != NULL; p++, len++);
	}
	hp->dof = phgRealloc_(hp->dof, (len + 2) * sizeof(*(hp->dof)),
					len * sizeof(*(hp->dof)));
	hp->dof[len] = dof;
	hp->dof[len + 1] = NULL;
    }

    /* init dof's cache */
    dof->cache_func = NULL;

    dof->count_vert = dim * (SpecialDofType(type) ? 0 : type->np_vert);
    dof->count_edge = dim * (SpecialDofType(type) ? 0 : type->np_edge);
    dof->count_face = dim * (SpecialDofType(type) ? 0 : type->np_face);
    dof->count_elem = dim * (SpecialDofType(type) ? 0 : type->np_elem);

    /* init DOF values */
    if (dof->userfunc == DofNoData) {
	dof->data = NULL;
	dof->userfunc = DofNoAction;
    }
    else if (userfunc != DofLambdaFunction) {
	phgDofSetDataByFunction(dof, userfunc);
    }

    Return dof;
    return dof;		/* to make MIPS cc happy */
}

void
phgDofFree(DOF **dof_ptr)
{
    DOF *dof = (dof_ptr == NULL ? NULL : *dof_ptr);

    MagicCheck(DOF, dof)

    if (dof == NULL)
	return;

    if (dof->g != NULL && !phgDofIsValid(dof->g, dof))
	phgError(1, "phgDofFree: invalid DOF object \"%s\".\n", dof->name);

    phgInfo(2, "freeing DOF \"%s\" (created at %s:%d).\n", 
		dof->name == NULL ? "nonname" : dof->name,
		dof->srcfile, dof->srcline);

    /* detach the dof from the mesh */
    if (phgGridIsValid(dof->g)) {
	DOF **p, **p0 = NULL;
	for (p = dof->g->dof; p != NULL && *p != NULL; p++) {
	    if (*p == dof)
		p0 = p;
	}
	if (p == NULL || p0 == NULL)
	    phgError(1, "%s:%d: unexpected error (grid freed before DOF?)\n",
		     __FILE__, __LINE__);
	assert(p > p0);
	memmove(p0, p0 + 1, sizeof(*p0) * (p - p0));
    }

    /* Note: don't use the DofIsHP macro because some functions, e.g.,
     * phgQuadFaceJumpN, may explicitly set dof->type to NULL (in this
     * case both dof->type and dof->hp are NULL pointers) before
     * calling phgDofFree */

    if (dof->hp != NULL) {
	/* detach the dof from dof->hp */
	DOF **p, **p0 = NULL;
	for (p = dof->hp->dof; p != NULL && *p != NULL; p++) {
	    if (*p == dof) {
		assert(p0 != dof_ptr);	/* trap phgDofFree(g->dofs + i) */
		p0 = p;
	    }
	}
	if (p == NULL || p0 == NULL)
	    phgError(1, "%s:%d: unexpected error (grid freed before DOF?)\n",
		     __FILE__, __LINE__);
	assert(p > p0);
	memmove(p0, p0 + 1, sizeof(*p0) * (p - p0));
	phgHPFree(&dof->hp);
    }

    if (dof->type != NULL)
	phgDofUnrefId(dof->type);

    if (dof->DB_masks != NULL)
	phgFree(dof->DB_masks);
    phgFree(dof->name);
    phgFree(dof->data);
    phgDofClearCache(NULL, dof, NULL, NULL, TRUE);

    phgFree(dof);

    *dof_ptr = NULL;

    return;
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
    if (!SpecialDofUserFunction(func))
	phgDofSetDataByFunction(u, func);
}

void
phgDofSetLambdaFunction(DOF *u, DOF_USER_FUNC_LAMBDA func)
{
    assert(u->type != DOF_CONSTANT);
    u->userfunc = DofLambdaFunction;
    u->userfunc_lambda = func;
    phgDofSetDataByLambdaFunction(u, func);
}

BTYPE
phgDofSetDirichletBoundaryMask(DOF *u, BTYPE mask)
{
    BTYPE old_mask = u->DB_mask;

    u->DB_mask = mask;
    if (u->DB_masks != NULL) {
	phgFree(u->DB_masks);
	u->DB_masks = NULL;
    }

    return old_mask;
}

void
phgDofSetDirichletBoundaryMasks(DOF *u, const BTYPE masks[])
{
    if (u->dim <= 1) {
	assert(u->DB_masks == NULL);
	u->DB_mask = masks[0];
	return;
    }

    if (u->DB_masks == NULL)
	u->DB_masks = phgAlloc(u->dim * sizeof(*u->DB_masks));

    memcpy(u->DB_masks, masks, u->dim * sizeof(*u->DB_masks));
    return;
}

BOOLEAN
phgDofIsValid(GRID *g, DOF *dof)
{
    DOF **p;

    if (g == NULL || (p = g->dof) == NULL)
	return FALSE;

    while (*p != dof && *p != NULL)
	p++;

    return *p == dof;
}

DOF *
phgDofIndex2Pointer(GRID *g, int id)
{
    DOF **p;

    if ((p = g->dof) == NULL)
	return NULL;

    while (id-- > 0 && *p != NULL) {
	p++;
    }

    return *p;
}

int
phgDofPointer2Index(DOF *dof)
{
    DOF **p;

    if (dof == NULL)
	return -1;

    p = dof->g->dof;
    while (*p != NULL && *p != dof) {
	p++;
    }

    return (*p == dof) ? p - dof->g->dof : -1;
}

static void
clear_cache0(DOF_TYPE *type, HP_TYPE *hp, BOOLEAN final)
{
    if (hp != NULL) {
	assert(type == NULL);
	phgQuadClearDofCache(&hp->cache_basfunc, NULL, final);
	phgQuadClearDofCache(&hp->cache_basgrad, NULL, final);
	phgQuadClearDofCache(&hp->cache_gradient, NULL, final);
	phgQuadClearDofCache(&hp->cache_curl, NULL, final);
	phgQuadClearDofCache(&hp->cache_hessian, NULL, final);
	return;
    }

    assert(type != NULL);
    if (final || !type->invariant) {
	phgQuadClearDofCache(&type->cache_basfunc, NULL, final);
	phgQuadClearDofCache(&type->cache_basgrad, NULL, final);
    }
    phgQuadClearDofCache(&type->cache_gradient, NULL, final);
    phgQuadClearDofCache(&type->cache_curl, NULL, final);
    phgQuadClearDofCache(&type->cache_hessian, NULL, final);

    return;
}

void
phgDofClearCache(GRID *g, DOF *dof, DOF_TYPE *type, HP_TYPE *hp, BOOLEAN final)
/*
 * 1. If 'dof' != NULL, then clear cached function values in 'dof' and return.
 *
 * 2. Otherwise if 'g' != NULL, then clear cached function/basis values in
 *    'g'->dof[] and return.
 *
 * 3. Otherwise if 'dof' != NULL, then clear cached function values in
 *    'dof' and return.
 *
 * 4. Otherwise if 'type' != NULL, then clear caches in 'type' and return.
 *
 * 5. Otherwise if 'hp' != NULL, then clear caches in 'hp' and return.
 *
 * 6. Otherwise clear all caches in type_list[] (unused).
 *
 * 'cache_basfunc' and 'cache_basgrad' are cleared only if:
 * 	'final' == TRUE || type == NULL || type->invariant == FALSE
 */
{
    DOF **pdof;
    size_t i;

    if (dof != NULL) {
	phgQuadClearDofCache(&(dof->cache_func), NULL, final);
	return;
    }

    if (g != NULL) {
	if (g->dof == NULL)
	    return;
	for (pdof = g->dof; *pdof != NULL; pdof++) {
	    clear_cache0((*pdof)->type, (*pdof)->hp, final);
	    phgQuadClearDofCache(&((*pdof)->cache_func), NULL, final);
	}
	return;
    }

    if (type != NULL || hp != NULL) {
	clear_cache0(type, hp, final);
	return;
    }

    /* clear all quad caches */
    assert(FALSE);		/* unused */
    for (i = 0; i < type_list_count; i++) {
	clear_cache0(type_list[i].type, NULL, final);
    }

    return;
}

void
phgDofMapEdgeData(const DOF_TYPE *type, const ELEMENT *e, int edge_no, int *M)
/* returns in M[] ordering of data according to edge direction */
{
    int i;

    if (type->np_edge == 0)
	return;

    if (type->np_edge <= 1 || type->BasFuncs == NULL || !type->is_nodal
	|| MapVertex(e, GetEdgeVertex(edge_no, 0))
	< MapVertex(e, GetEdgeVertex(edge_no, 1))) {
	for (i = 0; i < type->np_edge; i++)
	    M[i] = i;
	return;
    }

    if (type->base_type != NULL) {
	int M0[type->base_type->np_edge];
	int d = type->np_edge / type->base_type->np_edge;
	assert(d == type->dim / type->base_type->dim);
	phgDofMapEdgeData(type->base_type, e, edge_no, M0);
	for (i = 0; i < type->np_edge; i++)
	    M[i] = d * M0[i / d] + (i % d);
	return;
    }

    for (i = 0; i < type->np_edge; i++)
	M[i] = type->np_edge - 1 - i;
}

void
phgDofMapFaceData(const DOF_TYPE *type, const ELEMENT *e, int face_no, int *M)
/* compute ordering of data according to global indices of face vertices */
{
    int i;

    if (type->np_face <= 1 || type->DofMap == NULL) {
	/* Note: for non-point DOF the bases are already ordered according to
	 * global face orientation, return an identical map in this case */
	for (i = 0; i < type->np_face; i++)
	    M[i] = i;
	return;
    }

    type->DofMap(type, e, face_no, M);

    return;
}

void
phgDofMapElementData(const DOF_TYPE *type, const ELEMENT *e, int *M)
/* compute ordering of data according to global indices of vertices */
{
    int i;

    if (type->np_elem <= 1 || type->DofMap == NULL || ElementIsInOrder(e)) {
	for (i = 0; i < type->np_elem; i++)
	    M[i] = i;
	return;
    }

    type->DofMap(type, e, -1, M);

    return;
}

int
nonhp_bases_on_face(DOF_TYPE *type, ELEMENT *e, int face_no, SHORT bases[])
{
    int i, j, np, m0, m, v0, v1, v[4], *M;
    SHORT *bases0 = bases;

    if (type->base_type != NULL) {
	int dim = type->dim / type->base_type->dim, n;
	if (dim == 1 || bases == NULL) {
	    n = nonhp_bases_on_face(type->base_type, e, face_no, bases);
	}
	else {
	    SHORT tmp[type->base_type->nbas];
	    int i, j;
	    n  = nonhp_bases_on_face(type->base_type, e, face_no, tmp);
	    for (i = 0; i < n; i++)
		for (j = 0; j < dim; j++)
		    bases[i * dim + j] = tmp[i] * dim + j;
	}
	return n;
    }

    /* Note: the following lines are to be removed in the future, after
     * adding an appropriate flag in DOF_TYPE indicating which bases
     * functions are nonzero on a given face */
    if (type == DOF_P0) {
	if (bases != NULL)
	    bases[0] = 0;
	return 1;
    }

    if (bases == NULL)
	return 3 * (type->np_vert + type->np_edge) + type->np_face;

    GetFaceVertices(e, face_no, v[0], v[1], v[2], v[3]);

    M = phgAlloc(sizeof(*M) *
	(type->np_edge > type->np_face ? type->np_edge : type->np_face));

    /* The VERTEX bases */
    np = type->np_vert;
    for (i = 0; i < 3; i++) {
	m = v[i] * np;
	for (j = 0; j < np; j++)
	    *(bases++) = m++;
    }

    /* The EDGE bases */
    m0 = NVert * np;
    np = type->np_edge;
    for (i = 0; i < 3; i++) {
	switch (i) {
	    case 0: 	v0 = v[0]; v1 = v[1]; break;
	    case 1: 	v0 = v[0]; v1 = v[2]; break;
	    case 2: 	v0 = v[1]; v1 = v[2]; break;
	    default:	v0 = v1 = 0;	/* to make gcc happy */
	}
	j = GetEdgeNo(v0,v1);
	m = m0 + np * j;
	phgDofMapEdgeData(type, e, j, M);
	for (j = 0; j < np; j++)
	    bases[M[j]] = m++;
	bases += np;
    }

    /* The FACE bases */
    m = m0 + NEdge * np;
    np = type->np_face;
    m += face_no * np;
    phgDofMapFaceData(type, e, face_no, M);
    for (j = 0; j < np; j++)
	bases[M[j]] = m++;
    bases += np;

    phgFree(M);

    return (int)(bases - bases0);
}

int
phgDofGetBasesOnFace(DOF *dof, ELEMENT *e, int face_no, SHORT bases[])
/* computes the list of element bases of DOF on the element 'e' which are not
 * identically zero on the face 'face_no' (including VERTEX, EDGE, FACE and
 * VOLUME bases), returns number of bases.
 *
 * The bases are ordered according to the global orientation of the FACE,
 * i.e., global indices of the three vertices of the face, such that the
 * ordering is identical for the two elements sharing the face.
 *
 * Note:
 *   1.	For non hp DOF, the length of the array bases[] is equal to
 * 	3 * (type->np_vert + type->np_edge) + type->np_face (or 1 if type
 * 	is DOF_P0), where type is the base type of dof->type.
 *
 *   2. If bases == NULL the function simply returns number of
 * 	bases on the face.
 *
 *   3. The current implementation assumes only VERTEX/EDGE/FACE bases
 *      on the face are not zero (except for DG0 or P0).
 *
 *      TODO: add a flag in DOF_TYPE to indicate whether it has this property.
 */
{
    DOF_TYPE *type;
    int i, j, k, m, np, v0, v1, v[4], *M, index[NEdge + 1];
    SHORT *bases0 = bases;

#if 0
    FE_SPACE fe_space = (!DofIsHP(dof) ? dof->type->fe_space :
					 dof->hp->info->fe_space);

    /* The next test is commented out by Leng Wei, since this function is also
     * needed in handling Neumann BC with H1 elements. */
    if (fe_space != FE_L2 && fe_space != FE_None)
	phgError(1, "%s: don't know how to handle DOF \"%s\" "
		    "(only DG types are allowed).\n", __func__, dof->name);
#endif

    if (!DofIsHP(dof))
	return nonhp_bases_on_face(dof->type, e, face_no, bases);

    assert((dof->hp->info->flags & ELEM_FLAG));
    type = dof->hp->info->types[dof->hp->elem_order[e->index]];

    if (type->base_type != NULL) {
	assert(type->base_type->dim == type->dim);	/* TODO otherwise */
	assert(!(dof->hp->info->flags & (VERT_FLAG | EDGE_FLAG | FACE_FLAG)));
	return nonhp_bases_on_face(type->base_type, e, face_no, bases);
    }

    GetFaceVertices(e, face_no, v[0], v[1], v[2], v[3]);

    if (bases == NULL) {
	np = 0;
	if (dof->hp->info->flags & VERT_FLAG)
	    np += dof->hp->info->types[dof->hp->vert_order
                                [e->verts[v[0]]]]->np_vert +
		  dof->hp->info->types[dof->hp->vert_order
				[e->verts[v[1]]]]->np_vert +
		  dof->hp->info->types[dof->hp->vert_order
				[e->verts[v[2]]]]->np_vert;
	if (dof->hp->info->flags & EDGE_FLAG)
	    np += dof->hp->info->types[dof->hp->edge_order
                                [e->edges[GetEdgeNo(v[0],v[1])]]]->np_edge +
		  dof->hp->info->types[dof->hp->edge_order
                                [e->edges[GetEdgeNo(v[0],v[2])]]]->np_edge +
		  dof->hp->info->types[dof->hp->edge_order
                                [e->edges[GetEdgeNo(v[1],v[2])]]]->np_edge;
	if (dof->hp->info->flags & FACE_FLAG)
	    np += dof->hp->info->types[dof->hp->face_order
				[e->faces[face_no]]]->np_face;
	return np;
    }

    M = phgAlloc(sizeof(*M) *
	(type->np_edge > type->np_face ? type->np_edge : type->np_face));

    /* The VERTEX bases */
    m = 0;
    if (dof->hp->info->flags & VERT_FLAG) {
	index[0] = j = 0;
	for (i = 0; i < NVert; i++) {
	    j += dof->hp->info->types
			[dof->hp->vert_order[e->verts[i]]]->np_vert;
	    index[i + 1] = j;
	}
	for (i = 0; i < 3; i++)
	    for (k = index[v[i]]; k < index[v[i] + 1]; k++)
		*(bases++) = k;
	m = index[NVert];
    }

    /* The EDGE bases */
    if (dof->hp->info->flags & EDGE_FLAG) {
	index[0] = j = 0;
	for (i = 0; i < NEdge; i++) {
	    j += dof->hp->info->types
			[dof->hp->edge_order[e->edges[i]]]->np_edge;
	    index[i + 1] = j;
	}
	for (i = 0; i < 3; i++) {
	    switch (i) {
		case 0: 	v0 = v[0]; v1 = v[1]; break;
		case 1: 	v0 = v[0]; v1 = v[2]; break;
		case 2: 	v0 = v[1]; v1 = v[2]; break;
		default:	v0 = v1 = 0;	/* to make gcc happy */
	    }
	    j = GetEdgeNo(v0,v1);
	    type = dof->hp->info->types[dof->hp->edge_order[e->edges[j]]];
	    phgDofMapEdgeData(type, e, j, M);
	    for (k = index[j]; k < index[j + 1]; k++)
		bases[M[k - index[j]]] = k + m;
	    bases += index[j + 1] - index[j];
	}
	m += index[NEdge];
    }

    /* The FACE bases */
    if (dof->hp->info->flags & FACE_FLAG) {
	for (i = 0; i < face_no; i++)
	    m += dof->hp->info->types
			[dof->hp->face_order[e->faces[i]]]->np_face;
	type = dof->hp->info->types[dof->hp->face_order[e->faces[i]]];
	np = type->np_face;
	phgDofMapFaceData(type, e, i, M);
	for (k = 0; k < np; k++)
	    bases[M[k]] = m + k;
	bases += np;
    }

    phgFree(M);

    return (int)(bases - bases0);
}

VEF_MAP *
phgDofSetupVEFMap_(GRID *g, DOF *dof, int flags, FORALL_TYPE for_what)
/* creates vertex/edge/face to element maps.
 *
 * 'for_what' indicates to build the map for which subset of the elements.
 *
 * 'flags' is a bitwise combination of VERT_FLAG, EDGE_FLAG and FACE_FLAG,
 * indicating which maps are needed.
 *
 * 'dof' is optional and can be set to NULL, it's also used to indicate which
 * maps are needed.
 *
 * Only the maps required by both 'flags' and 'dof' are constructed. */
{
    size_t size, size0, size1, size2;
    DOF_TYPE *type;
    ELEMENT *e;
    VEF_MAP *vef;
    int i;
    INT j;
#if USE_OMP
    omp_lock_t *locks0, *locks1, *locks2;
#endif	/* USE_OMP */

    if (dof != NULL) {
	type = dof->type;
	if (type == NULL)
	    type = dof->hp->info->types[dof->hp->max_order];
	while (type->base_type != NULL)
	    type = type->base_type;
    }
    else {
	type = DOF_P4;		/* has VERTEX, EDGE, FACE, and VOLUME data */
    }

    if (type->np_vert == 0)
	flags &= ~VERT_FLAG;
    if (type->np_edge == 0)
	flags &= ~EDGE_FLAG;
    if (type->np_face == 0)
	flags &= ~FACE_FLAG;

    size0 = ((flags & VERT_FLAG) ? g->nvert : 0);
    size1 = ((flags & EDGE_FLAG) ? g->nedge : 0);
    size2 = ((flags & FACE_FLAG) ? g->nface : 0);
    size = size0 + size1 + size2;

    vef = phgCalloc(1, sizeof(*vef));

    if (size == 0)
	return vef;

    vef->Vmap = phgCalloc(size, sizeof(*vef->Vmap));
    vef->Emap = vef->Vmap + size0;
    vef->Fmap = vef->Emap + size1;

    vef->Vind = phgCalloc(size, sizeof(*vef->Vind));
    vef->Eind = vef->Vind + size0;
    vef->Find = vef->Eind + size1;

#if USE_OMP
    locks0 = phgAlloc(size * sizeof(*locks0));
    locks1 = locks0 + size0;
    locks2 = locks1 + size1;
    for (j = 0; j < size; j++)
	omp_init_lock(&locks0[j]);
# pragma omp parallel for schedule(static) private(e, i, j)
#endif	/* USE_OMP */
    ForAllElementsBegin_(g, e, for_what) {
	for (i = 0; (flags & VERT_FLAG) && i < NVert; i++) {
	    if (size0 == 0)
		break;
	    j = e->verts[i];
	    if (vef->Vmap[j] != NULL)
		continue;
#if USE_OMP
	    omp_set_lock(&locks0[j]);
#endif	/* USE_OMP */
	    vef->Vmap[j] = e;
	    vef->Vind[j] = i;
#if USE_OMP
	    omp_unset_lock(&locks0[j]);
#endif	/* USE_OMP */
	}

	for (i = 0; (flags & EDGE_FLAG) && i < NEdge; i++) {
	    if (size1 == 0)
		break;
	    j = e->edges[i];
	    if (vef->Emap[j] != NULL)
		continue;
#if USE_OMP
	    omp_set_lock(&locks1[j]);
#endif	/* USE_OMP */
	    vef->Emap[j] = e;
	    vef->Eind[j] = i;
#if USE_OMP
	    omp_unset_lock(&locks1[j]);
#endif	/* USE_OMP */
	}

	for (i = 0; (flags & FACE_FLAG) && i < NFace; i++) {
	    if (size2 == 0)
		break;
	    j = e->faces[i];
	    if (vef->Fmap[j] != NULL)
		continue;
#if USE_OMP
	    omp_set_lock(&locks2[j]);
#endif	/* USE_OMP */
	    vef->Fmap[j] = e;
	    vef->Find[j] = i;
#if USE_OMP
	    omp_unset_lock(&locks2[j]);
#endif	/* USE_OMP */
	}
    } ForAllElementsEnd
#if USE_OMP
    for (j = 0; j < size; j++)
	omp_destroy_lock(&locks0[j]);
    phgFree(locks0);
#endif	/* USE_OMP */

    if ((flags & FACE_FLAG) && g->Fmap_elems == NULL && g->nface != 0 &&
	(for_what == FOR_OWNED || for_what == FOR_ALL)) {
	/* save face-to-element map in g->Fmap_{elems,faces} */
	assert(g->Fmap_faces == NULL);
	g->Fmap_elems = phgAlloc(g->nface * sizeof(*vef->Fmap));
	memcpy(g->Fmap_elems, vef->Fmap, g->nface * sizeof(*vef->Fmap));
	g->Fmap_faces = phgAlloc(g->nface * sizeof(*vef->Find));
	memcpy(g->Fmap_faces, vef->Find, g->nface * sizeof(*vef->Find));
    }

    return vef;
}

void
phgDofFreeVEFMap(VEF_MAP **vef_ptr)
{
    VEF_MAP *vef = *vef_ptr;

    if (vef == NULL)
	return;

    phgFree(vef->Vmap);
    phgFree(vef->Vind);
    phgFree(vef);

    *vef_ptr = NULL;
}

VEF_MAP2 *
phgDofSetupVEFMap2_(GRID *g, DOF *dof, int flags, int mask,
			const BTYPE *types_vert, const BTYPE *types_edge,
			const BTYPE *types_face, FORALL_TYPE for_what)
/* Creates local vertex/edge/face to element links.
 *
 * 'flags' is a bitwise combination of VERT_FLAG, EDGE_FLAG and FACE_FLAG,
 * indicating which maps are needed.
 *
 * 'dof' is optional and can be set to NULL, it's also used to indicate which
 * maps are needed.
 *
 * The 'mask' and 'types_*' arrays are used to select subsets of vertices,
 * edges and faces. The links are generated for 'xxxx' i if and only if:
 * 	types_xxxx == NULL || (types_xxxx[i] & mask) != 0
 * where 'xxxx' denotes 'vert', 'edge' or 'face'.
 *
 * Only the maps required by both 'flags' and 'dof' are constructed.
 *
 * Difference between VEFMap and VEFMap2:
 * VEFMap links a vertex/edge/face to ONE OF the local containing elements,
 * while VEFMap2 links a vertex/edge/face to ALL local containing elements. 
 *
 */
{
    size_t size, size0, size1, size2;
    DOF_TYPE *type;
    ELEMENT *e, **elist;
    BYTE *Vn, *En, *Fn;
    VEF_MAP2 *vef;
    BYTE *idxlist;
    INT i, ne, list_size, pos;
    int ii, jj;

    if (dof != NULL) {
	type = dof->type;
	if (type == NULL)
	    type = dof->hp->info->types[dof->hp->max_order];
	while (type->base_type != NULL)
	    type = type->base_type;
    }
    else {
	type = DOF_P4;		/* has VERTEX, EDGE, FACE, and VOLUME data */
    }

    size0 = (((flags & VERT_FLAG) && type->np_vert > 0) ? g->nvert : 0);
    size1 = (((flags & EDGE_FLAG) && type->np_edge > 0) ? g->nedge : 0);
    size2 = (((flags & FACE_FLAG) && type->np_face > 0) ? g->nface : 0);
    size = size0 + size1 + size2;

    vef = phgCalloc(1, sizeof(*vef));

    if (size == 0)
	return vef;

    vef->Vmap = phgCalloc(size, sizeof(*vef->Vmap));
    vef->Emap = vef->Vmap + size0;
    vef->Fmap = vef->Emap + size1;

    vef->Vind = phgCalloc(size, sizeof(*vef->Vind));
    vef->Eind = vef->Vind + size0;
    vef->Find = vef->Eind + size1;

    vef->Vsize = phgCalloc(size, sizeof(*vef->Vsize));
    vef->Esize = vef->Vsize + size0;
    vef->Fsize = vef->Esize + size1;

    /* record local pos */
    Vn = phgCalloc(size, sizeof(*Vn));
    En = Vn + size0;
    Fn = En + size1;

    ne = 0;
    ForAllElements_(g, e, for_what) {
	for (ii = 0; (flags & VERT_FLAG) && ii < NVert; ii++) {
	    if (size0 == 0)
		break;
	    i = e->verts[ii];
	    if (types_vert != NULL && !(types_vert[i] & mask))
		continue;
	    vef->Vsize[i]++;
	    ne++;
	}

	for (ii = 0; (flags & EDGE_FLAG) && ii < NEdge; ii++) {
	    if (size1 == 0)
		break;
	    i = e->edges[ii];
	    if (types_edge != NULL && !(types_edge[i] & mask))
		continue;
	    vef->Esize[i]++;
	    ne++;
	}

	for (ii = 0; (flags & FACE_FLAG) && ii < NFace; ii++) {
	    if (size2 == 0)
		break;
	    i = e->faces[ii];
	    if (types_face != NULL && !(types_face[i] & mask))
		continue;
	    vef->Fsize[i]++;
	    ne++;
	}
    }

    list_size = ne;
    elist = vef->elist = phgCalloc(list_size, sizeof(*elist));
    idxlist = vef->idxlist = phgCalloc(list_size, sizeof(*idxlist));

    pos = 0;
    if (size0 > 0)
	for (i = 0; (flags & VERT_FLAG) && i < g->nvert; i++) {
	    if (types_vert != NULL && !(types_vert[i] & mask))
		continue;
	    vef->Vmap[i] = elist + pos;
	    vef->Vind[i] = idxlist + pos;
	    pos += vef->Vsize[i];
	}

    if (size1 > 0)
	for (i = 0; (flags & EDGE_FLAG) && i < g->nedge; i++) {
	    if (types_edge != NULL && !(types_edge[i] & mask))
		continue;
	    vef->Emap[i] = elist + pos;
	    vef->Eind[i] = idxlist + pos;
	    pos += vef->Esize[i];
	}

    if (size2 > 0)
	for (i = 0; (flags & FACE_FLAG) && i < g->nface; i++) {
	    if (types_face != NULL && !(types_face[i] & mask))
		continue;
	    vef->Fmap[i] = elist + pos;
	    vef->Find[i] = idxlist + pos;
	    pos += vef->Fsize[i];
	}

    ForAllElements_(g, e, for_what) {
	for (ii = 0; (flags & VERT_FLAG) && ii < NVert; ii++) {
	    if (size0 == 0)
		break;
	    i = e->verts[ii];
	    if (types_vert != NULL && !(types_vert[i] & mask))
		continue;
	    jj = Vn[i];
	    vef->Vmap[i][jj] = e;
	    vef->Vind[i][jj] = ii;
	    Vn[i]++;
	}

	for (ii = 0; (flags & EDGE_FLAG) && ii < NEdge; ii++) {
	    if (size1 == 0)
		break;
	    i = e->edges[ii];
	    if (types_edge != NULL && !(types_edge[i] & mask))
		continue;
	    jj = En[i];
	    vef->Emap[i][jj] = e;
	    vef->Eind[i][jj] = ii;
	    En[i]++;
	}

	for (ii = 0; (flags & FACE_FLAG) && ii < NFace; ii++) {
	    if (size2 == 0)
		break;
	    i = e->faces[ii];
	    if (types_face != NULL && !(types_face[i] & mask))
		continue;
	    jj = Fn[i];		
	    vef->Fmap[i][jj] = e;
	    vef->Find[i][jj] = ii;
	    Fn[i]++;
	}
    }

    phgFree(Vn);

    return vef;
}

void
phgDofFreeVEFMap2(VEF_MAP2 **vef_ptr)
{
    VEF_MAP2 *vef = *vef_ptr;

    if (vef == NULL)
	return;

    phgFree(vef->elist);
    phgFree(vef->idxlist);
    phgFree(vef->Vmap);
    phgFree(vef->Vind);
    phgFree(vef->Vsize);
    phgFree(vef);

    *vef_ptr = NULL;
}

void
phgDofDump(DOF *dof)
{
    char s[1024], *fmt;
    int i, j, k, index = 0, d;
    GRID *g;
    FLOAT *data, *pt = NULL, x, y, z;
    ELEMENT *e = NULL;
    COORD *p0 = NULL, *p1 = NULL, *p2 = NULL, *p3 = NULL;
    DOF_TYPE *type, *t, *t0;
    int *M = NULL, M_size = 0;
    VEF_MAP *vef = NULL;

    if (dof == NULL)
	return;

    if (dof->type == DOF_ANALYTIC) {
	if (dof->userfunc_lambda != NULL)
	    phgInfo(-1, "Analytic DOF \"%s\", dim = %d, lambda function: %p\n",
		    dof->name, dof->dim, dof->userfunc_lambda);
	else
	    phgInfo(-1, "Analytic DOF \"%s\", dim = %d, xyz function: %p\n",
		    dof->name, dof->dim, dof->userfunc);
	return;
    }
    else if (dof->type == DOF_CONSTANT) {
	if (dof->data == NULL) {
	    phgInfo(-1, "Constant DOF \"%s\", dim = %d, no data\n",
		    dof->name, dof->dim);
	    return;
	}
	phgInfo(-1, "Constant DOF \"%s\", dim = %d, values:\n",
		dof->name, dof->dim);
	for (i = 0, j = 0; i < dof->dim / 4; i++, j += 4) {
	    phgInfo(-1, "    %10.4le %10.4le %10.4le %10.4le\n",
		    (double)dof->data[j], (double)dof->data[j + 1],
		    (double)dof->data[j + 2], (double)dof->data[j + 3]);
	}
	switch (dof->dim % 4) {
	    case 1:
		phgInfo(-1, "    %10.4le\n", (double)dof->data[j]);
		break;
	    case 2:
		phgInfo(-1, "    %10.4le %10.4le\n", (double)dof->data[j],
			(double)dof->data[j + 1]);
		break;
	    case 3:
		phgInfo(-1, "    %10.4le %10.4le %10.4le\n",
			(double)dof->data[j], (double)dof->data[j + 1],
			(double)dof->data[j + 2]);
		break;
	}
	return;
    }

    if (!DofIsHP(dof))
	phgInfo(-1, "DOF \"%s\": dim=%dx%d, "
		    "np_vert=%d, np_edge=%d, np_face=%d, np_elem=%d\n",
		    dof->name, dof->dim, dof->type->dim,
		    dof->type->np_vert, dof->type->np_edge,
		    dof->type->np_face, dof->type->np_elem);
    else
	phgInfo(-1, "DOF \"%s\": dim=%dx%d, hierachical.\n",
		    dof->name, dof->dim, DofTypeDim(dof));

    if ((g = dof->g) == NULL || (data = dof->data) == NULL) {
	phgInfo(-1, "no data.\n");
	return;
    }

    type = dof->type;
    if (type == NULL)
	type = dof->hp->info->types[dof->hp->max_order];

    vef = phgDofSetupVEFMap(g, dof, EDGE_FLAG | FACE_FLAG);

    if ((k = type->np_vert) > 0) {
	fmt = (k <= 1) ? "%5d  " : ((k < 10) ? "%5d[%d]  " : "%5d[%02d]  ");
	phgInfo(-1, "Vertices:\n");
	for (i = 0; i < g->nvert; i++) {
	    t = DofTypeOnVertex(dof, i);
	    if (g->types_vert[i] == UNREFERENCED) {
		    data += dof->dim * t->np_vert;
		continue;
	    }
	    t0 = t;
	    while (t0->points == NULL && t0->base_type != NULL)
		t0 = t0->base_type;
	    for (j = 0; j < t->np_vert; j++) {
		sprintf(s, fmt, GlobalVertex(g, i), j);
		if (t0->points != NULL)
		    sprintf(s + strlen(s), "%s(%6.3lf,%6.3lf,%6.3lf) = ",
			    dof->name, (double)g->verts[i][0],
			    (double)g->verts[i][1], (double)g->verts[i][2]);
		for (k = 0; k < dof->dim; k++)
		    sprintf(s + strlen(s), "%10.4le ", (double)*(data++));
		phgInfo(-1, "%s\n", s);
	    }
	}
    }

    if ((k = type->np_edge) > 0) {
	if (M_size < k) {
	    phgFree(M);
	    M = phgAlloc((M_size = k) * sizeof(*M));
	}
	fmt = (k <= 1) ? "%5d  " : ((k < 10) ? "%5d[%d]  " : "%5d[%02d]  ");
	phgInfo(-1, "Edges:\n");
	for (i = 0; i < g->nedge; i++) {
	    if (!DofIsHP(dof))
		t = type;
	    else
	 	t = dof->hp->info->types[dof->hp->edge_order[i]];
	    if (g->types_edge[i] == UNREFERENCED) {
		if (!DofIsHP(dof))
		    data += dof->count_edge;
		else
		    data += dof->dim *
			(dof->hp->edge_index[i + 1] - dof->hp->edge_index[i]);
		continue;
	    }
	    t0 = t;
	    while (t0->points == NULL && t0->base_type != NULL)
		t0 = t0->base_type;
	    d = t->dim / t0->dim;
	    e = vef->Emap[i];
	    index = vef->Eind[i];
	    if (t0->points != NULL) {
		assert(vef != NULL && vef->Emap[i] != NULL);
		p0 = g->verts + e->verts[GetEdgeVertex(index, 0)];
		p1 = g->verts + e->verts[GetEdgeVertex(index, 1)];
	    }
	    phgDofMapEdgeData(t, e, index, M);
	    for (j = 0; j < t->np_edge; j++) {
		sprintf(s, fmt, GlobalEdge(g, i), j);
		if (t0->points != NULL) {
		    pt = t0->points + t0->np_vert + 2 * (j / d);
		    x = (*p0)[0] * pt[0] + (*p1)[0] * pt[1];
		    y = (*p0)[1] * pt[0] + (*p1)[1] * pt[1];
		    z = (*p0)[2] * pt[0] + (*p1)[2] * pt[1];
		    sprintf(s + strlen(s), "%s(%6.3lf,%6.3lf,%6.3lf) = ",
			    dof->name, (double)x, (double)y, (double)z);
		}
		for (k = 0; k < dof->dim; k++)
		    sprintf(s + strlen(s), "%10.4le ",
			    (double)*(data + dof->dim * M[j] + k));
		phgInfo(-1, "%s\n", s);
	    }
	    data += dof->dim * t->np_edge;
	}
    }

    if ((k = type->np_face) > 0) {
	if (M_size < k) {
	    phgFree(M);
	    M = phgAlloc((M_size = k) * sizeof(*M));
	}
	fmt = (k <= 1) ? "%5d  " : ((k < 10) ? "%5d[%d]  " : "%5d[%02d]  ");
	phgInfo(-1, "Faces:\n");
	for (i = 0; i < g->nface; i++) {
	    if (!DofIsHP(dof))
		t = type;
	    else
	 	t = dof->hp->info->types[dof->hp->face_order[i]];
	    if (g->types_face[i] == UNREFERENCED) {
		if (!DofIsHP(dof))
		    data += dof->count_face;
		else
		    data += dof->dim *
			(dof->hp->face_index[i + 1] - dof->hp->face_index[i]);
		continue;
	    }
	    t0 = t;
	    while (t0->points == NULL && t0->base_type != NULL)
		t0 = t0->base_type;
	    d = t->dim / t0->dim;
	    e = vef->Fmap[i];
	    index = vef->Find[i];
	    if (t0->points != NULL) {
		assert(vef != NULL && vef->Fmap[i] != NULL);
		p0 = g->verts + e->verts[GetFaceVertex(index, 0)];
		p1 = g->verts + e->verts[GetFaceVertex(index, 1)];
		p2 = g->verts + e->verts[GetFaceVertex(index, 2)];
	    }
	    phgDofMapFaceData(t, e, index, M);
	    for (j = 0; j < t->np_face; j++) {
		sprintf(s, fmt, GlobalFace(g, i), j);
		if (t0->points != NULL) {
		    pt = t0->points + t0->np_vert + 2 * t0->np_edge +
				      3 * (j / d);
		    x = (*p0)[0] * pt[0] + (*p1)[0] * pt[1] + (*p2)[0] * pt[2];
		    y = (*p0)[1] * pt[0] + (*p1)[1] * pt[1] + (*p2)[1] * pt[2];
		    z = (*p0)[2] * pt[0] + (*p1)[2] * pt[1] + (*p2)[2] * pt[2];
		    sprintf(s + strlen(s), "%s(%6.3lf,%6.3lf,%6.3lf) = ",
			    dof->name, (double)x, (double)y, (double)z);
		}
		for (k = 0; k < dof->dim; k++)
		    sprintf(s + strlen(s), "%10.4le ",
			    (double)*(data + dof->dim * M[j] + k));
		phgInfo(-1, "%s\n", s);
	    }
	    data += dof->dim * t->np_face;
	}
    }

    if ((k = type->np_elem) > 0) {
	if (M_size < k) {
	    phgFree(M);
	    M = phgAlloc((M_size = k) * sizeof(*M));
	}
	fmt = (k <= 1) ? "%5d  " : ((k < 10) ? "%5d[%d]  " : "%5d[%02d]  ");
	phgInfo(-1, "Elements:\n");
	for (i = 0; i < g->nelem; i++) {
	    if (!DofIsHP(dof))
		t = type;
	    else
	 	t = dof->hp->info->types[dof->hp->elem_order[i]];
	    if (g->types_elem[i] == UNREFERENCED) {
		if (!DofIsHP(dof))
		    data += dof->count_elem;
		else
		    data += dof->dim *
			(dof->hp->elem_index[i + 1] - dof->hp->elem_index[i]);
		continue;
	    }
	    t0 = t;
	    while (t0->points == NULL && t0->base_type != NULL)
		t0 = t0->base_type;
	    d = t->dim / t0->dim;
	    e = g->elems[i];
	    if (t0->points != NULL) {
		assert(e != NULL && e->index == i);
		p0 = g->verts + e->verts[0];
		p1 = g->verts + e->verts[1];
		p2 = g->verts + e->verts[2];
		p3 = g->verts + e->verts[3];
	    }
	    phgDofMapElementData(t, e, M);
	    for (j = 0; j < t->np_elem; j++) {
		sprintf(s, fmt, GlobalElement(g, i), j);
		if (t0->points != NULL) {
		    pt = t0->points + t0->np_vert + 2 * t0->np_edge +
				      3 * t0->np_face + 4 * (j / d);
		    x = (*p0)[0] * pt[0] + (*p1)[0] * pt[1] +
			(*p2)[0] * pt[2] + (*p3)[0] * pt[3];
		    y = (*p0)[1] * pt[0] + (*p1)[1] * pt[1] +
			(*p2)[1] * pt[2] + (*p3)[1] * pt[3];
		    z = (*p0)[2] * pt[0] + (*p1)[2] * pt[1] +
			(*p2)[2] * pt[2] + (*p3)[2] * pt[3];
		    sprintf(s + strlen(s), "%s(%6.3lf,%6.3lf,%6.3lf) = ",
			    dof->name, (double)x, (double)y, (double)z);
		}
		for (k = 0; k < dof->dim; k++)
		    sprintf(s + strlen(s), "%10.4le ",
				(double)*(data + dof->dim * M[j] + k));
		phgInfo(-1, "%s\n", s);
	    }
	    data += dof->dim * t->np_elem;
	}
    }

    phgDofFreeVEFMap(&vef);
    phgFree(M);
}

static ELEMENT *elem;
static FLOAT **pdata = NULL, **cdata = NULL, *cdata0[15], *cdata1[15];
static BOOLEAN cdata0_flag = FALSE, cdata1_flag = FALSE;
static int vmap0[NVert], emap0[NEdge], vmap1[NVert], emap1[NEdge];

static void
c2f_func0(DOF *dof, ELEMENT *e0, int bno, const FLOAT *lambda, FLOAT *values)
{
    FLOAT lam[Dim + 1];

    Unused(bno);

    lam[0] = lambda[vmap0[0]] + (lam[1] = lambda[3] * 0.5);
    lam[2] = lambda[vmap0[2]];
    lam[3] = lambda[vmap0[3]];

    CheckD_order(dof);
    phgDofEval_(dof, elem, lam, NULL, 0, values, pdata);
}

static void
c2f_func1(DOF *dof, ELEMENT *e1, int bno, const FLOAT *lambda, FLOAT *values)
{
    FLOAT lam[Dim + 1];

    Unused(bno);

    lam[1] = lambda[vmap1[1]] + (lam[0] = lambda[3] * 0.5);
    lam[2] = lambda[vmap1[2]];
    lam[3] = lambda[vmap1[3]];

    CheckD_order(dof);
    phgDofEval_(dof, elem, lam, NULL, 0, values, pdata);
}

static void
f2c_func(DOF *dof, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    int *map;
    FLOAT lam[Dim + 1], **data;

    Unused(bno);

    if ((e->children[1] == NULL) ||
	(lambda[0] >= lambda[1] && e->children[0] != NULL)) {
	/* use child 0 */
	map = vmap0;
	lam[map[0]] = lambda[0] - lambda[1];
	lam[3] = 2. * lambda[1];
	lam[map[2]] = lambda[2];
	lam[map[3]] = lambda[3];
	e = e->children[0];
	data = cdata0;
	if (cdata0_flag == FALSE) {
	    cdata0_flag = TRUE;
	    /* VERTEX data */
	    data[3] = cdata[0];
	    data[map[0]] = pdata[0];
	    data[map[2]] = pdata[2];
	    data[map[3]] = pdata[3];
	    /* EDGE data */
	    data[NVert + GetEdgeNo(map[0], 3)] = cdata[1];
	    data[NVert + GetEdgeNo(map[2], 3)] = cdata[3];
	    data[NVert + GetEdgeNo(map[3], 3)] = cdata[4];
	    data[NVert + GetEdgeNo(map[0], map[2])] = pdata[NVert + 1];
	    data[NVert + GetEdgeNo(map[0], map[3])] = pdata[NVert + 2];
	    data[NVert + GetEdgeNo(map[2], map[3])] = pdata[NVert + 5];
	    /* FACE data */
	    data[NVert + NEdge + map[2]] = cdata[5];
	    data[NVert + NEdge + map[3]] = cdata[7];
	    data[NVert + NEdge + map[0]] = cdata[9];
	    data[NVert + NEdge + 3] = pdata[NVert + NEdge + 1];
	    /* VOLUME data */
	    data[NVert + NEdge + NFace] = cdata[10];
	}
    }
    else {
	/* use child 1 */
	map = vmap1;
	lam[map[1]] = lambda[1] - lambda[0];
	lam[3] = 2. * lambda[0];
	lam[map[2]] = lambda[2];
	lam[map[3]] = lambda[3];
	e = e->children[1];
	data = cdata1;
	if (cdata1_flag == FALSE) {
	    cdata1_flag = TRUE;
	    /* VERTEX data */
	    data[3] = cdata[0];
	    data[map[1]] = pdata[1];
	    data[map[2]] = pdata[2];
	    data[map[3]] = pdata[3];
	    /* EDGE data */
	    data[NVert + GetEdgeNo(map[1], 3)] = cdata[2];
	    data[NVert + GetEdgeNo(map[2], 3)] = cdata[3];
	    data[NVert + GetEdgeNo(map[3], 3)] = cdata[4];
	    data[NVert + GetEdgeNo(map[1], map[2])] = pdata[NVert + 3];
	    data[NVert + GetEdgeNo(map[1], map[3])] = pdata[NVert + 4];
	    data[NVert + GetEdgeNo(map[2], map[3])] = pdata[NVert + 5];
	    /* FACE data */
	    data[NVert + NEdge + map[2]] = cdata[6];
	    data[NVert + NEdge + map[3]] = cdata[8];
	    data[NVert + NEdge + map[1]] = cdata[9];
	    data[NVert + NEdge + 3] = pdata[NVert + NEdge + 0];
	    /* VOLUME data */
	    data[NVert + NEdge + NFace] = cdata[11];
	}
    }

    CheckD_order(dof);
    phgDofEval_(dof, e, lam, NULL, 0, values, data);
}

void
phgDofInterpC2FGeneric(DOF *dof, ELEMENT *e, FLOAT **parent_data,
		       FLOAT **children_data)
/* Default coarse-to-fine interpolation function. It can be used
 * as the 'InterpC2F' member in any DOF_TYPE which has valid 'InitFunc'
 * and 'BasFuncs' members.
 *
 * Note:
 *	For better efficiency specialized functions should be implemented
 *	for specific DOF_TYPEs */
{
    ELEMENT *e0 = e->children[0], *e1 = e->children[1];
    FLOAT **child0_data, **child1_data;

    phgMapP2C(e, vmap0, emap0, 0);
    phgMapP2C(e, vmap1, emap1, 1);

    /* Note: only need child[01]_data for HBn, HCn and HFEBn with n >= 2 */
    /* FIXME: do we need to add a flag in DOF_TYPE for this test? */
    if (dof->type->order > 1 /* other tests */) {
	/* setup complete list of pointers for children data */
	child0_data = phgAlloc(2 * 15 * sizeof(*child0_data));
	child1_data = child0_data + 15;

	/* the vertices of child0 */
	child0_data[vmap0[0]] = parent_data[0];
	child0_data[vmap0[2]] = parent_data[2];
	child0_data[vmap0[3]] = parent_data[3];
	child0_data[3] = children_data[0];

	/* the vertices of child1 */
	child1_data[vmap1[1]] = parent_data[1];
	child1_data[vmap1[2]] = parent_data[2];
	child1_data[vmap1[3]] = parent_data[3];
	child1_data[3] = children_data[0];

	/* the edges of child0 */
	child0_data[4 + emap0[1]] = parent_data[5];	/* parent edge 0-2 */
	child0_data[4 + emap0[2]] = parent_data[6];	/* parent edge 0-3 */
	child0_data[4 + emap0[5]] = parent_data[9];	/* parent edge 2-3 */
	/* cut edge (v0 mapped to 2 or 0) */
	if (vmap0[0] == 2)
	    child0_data[4 + 5] = children_data[1];	/* cut edge (2-3) */
	else
	    child0_data[4 + 2] = children_data[1];	/* cut edge (0-3) */
	/* new edge on v2 (v2 mapped to 1) */
	child0_data[4 + 4] = children_data[3];		/* new edge (1-3) */
	/* new edge on v3 (v3 mapped to 0 or 2) */
	if (vmap0[3] == 0)
	    child0_data[4 + 2] = children_data[4];	/* new edge (0-3) */
	else
	    child0_data[4 + 5] = children_data[4];	/* new edge (2-3) */

	/* the edges of child1 */
	child1_data[4 + emap1[3]] = parent_data[7];	/* parent edge 1-2 */
	child1_data[4 + emap1[4]] = parent_data[8];	/* parent edge 1-3 */
	child1_data[4 + emap1[5]] = parent_data[9];	/* parent edge 2-3 */
	/* cut edge  (v1 mapped to 2 or 0) */
	if (vmap1[1] == 2)
	    child1_data[4 + 5] = children_data[2];	/* cut edge (2-3) */
	else
	    child1_data[4 + 2] = children_data[2];	/* cut edge (0-3) */
	/* new edge on v2 (v2 mapped to 2 or 1) */
	if (vmap1[2] == 2)
	    child1_data[4 + 5] = children_data[3];	/* new edge (2-3) */
	else
	    child1_data[4 + 4] = children_data[3];	/* new edge (1-3) */
	/* new edge on v3 (v3 mapped to 0, 1 or 2) */
	if (vmap1[3] == 0)
	    child1_data[4 + 2] = children_data[4];	/* new edge (0-3) */
	else if (vmap1[3] == 1)
	    child1_data[4 + 4] = children_data[4];	/* new edge (1-3) */
	else
	    child1_data[4 + 5] = children_data[4];	/* new edge (2-3) */

	/* the faces of child0 */
	child0_data[10 + 3] = parent_data[11];		/* parent face 1 */
	child0_data[10 + vmap0[0]] = children_data[9];	/* new face */
	child0_data[10 + vmap0[2]] = children_data[5];	/* cut face */
	child0_data[10 + vmap0[3]] = children_data[7];	/* cut face */

	/* the faces of child1 */
	child1_data[10 + 3] = parent_data[10];		/* parent face 0 */
	child1_data[10 + vmap1[1]] = children_data[9];	/* new face */
	child1_data[10 + vmap1[2]] = children_data[6];	/* cut face */
	child1_data[10 + vmap1[3]] = children_data[8];	/* cut face */

	/* the volume data */
	child0_data[14] = children_data[10];
	child1_data[14] = children_data[11];
    }
    else {
	child0_data = child1_data = NULL;
    }

    elem = e;
    pdata = parent_data;

    if (dof->type->np_vert > 0) {
	/* The new vertex */
	dof->type->InitFunc(dof, e0, VERTEX, 3, NULL,
			    c2f_func0, NULL, children_data[0], child0_data);
    }

    if (dof->type->np_edge > 0) {
	/* The cut-edges */
	dof->type->InitFunc(dof, e0, EDGE, emap0[0], NULL,
			    c2f_func0, NULL, children_data[1], child0_data);
	dof->type->InitFunc(dof, e1, EDGE, emap1[0], NULL,
			    c2f_func1, NULL, children_data[2], child1_data);
	/* The new edges */
	dof->type->InitFunc(dof, e0, EDGE, GetEdgeNo(vmap0[2], 3), NULL,
			    c2f_func0, NULL, children_data[3], child0_data);
	dof->type->InitFunc(dof, e0, EDGE, GetEdgeNo(vmap0[3], 3), NULL,
			    c2f_func0, NULL, children_data[4], child0_data);
    }

    if (dof->type->np_face > 0) {
	/* the cut faces of e0 */
	dof->type->InitFunc(dof, e0, FACE, vmap0[2], NULL,
			    c2f_func0, NULL, children_data[5], child0_data);
	dof->type->InitFunc(dof, e0, FACE, vmap0[3], NULL,
			    c2f_func0, NULL, children_data[7], child0_data);
	/* the cut faces of e1 */
	dof->type->InitFunc(dof, e1, FACE, vmap1[2], NULL,
			    c2f_func1, NULL, children_data[6], child1_data);
	dof->type->InitFunc(dof, e1, FACE, vmap1[3], NULL,
			    c2f_func1, NULL, children_data[8], child1_data);
	/* the new face */
	dof->type->InitFunc(dof, e0, FACE, vmap0[0], NULL,
			    c2f_func0, NULL, children_data[9], child0_data);
    }

    if (dof->type->np_elem > 0) {
	/* volume data on e0 */
	dof->type->InitFunc(dof, e0, VOLUME, 0, NULL,
			    c2f_func0, NULL, children_data[10], child0_data);
	/* volume data on e1 */
	dof->type->InitFunc(dof, e1, VOLUME, 0, NULL,
			    c2f_func1, NULL, children_data[11], child1_data);
    }

    phgFree(child0_data);
}

void
phgDofInterpF2CGeneric(DOF *dof, ELEMENT *e, FLOAT **parent_data,
		       FLOAT **children_data)
/* Default fine-to-coarse interpolation function. It can be used
 * as the 'InterpF2C' member in any DOF_TYPE which has valid 'InitFunc'
 * and 'BasFuncs' members.
 *
 * Note:
 *   1. It is assumed that 'parent_data' already contains valid values at
 *	the locations common to parent and children, only values at edge 0,
 *	faces 2 and 3, and element data need to be computed.
 *
 *   2. For better efficiency specialized functions should be implemented
 *	for specific DOF_TYPEs */
{
#if USE_MPI
    ELEMENT *e0 = e->children[0], *e1 = e->children[1];
#endif

    pdata = parent_data;
    cdata = children_data;
    cdata0_flag = FALSE;
    cdata1_flag = FALSE;

    phgMapP2C(e, vmap0, emap0, 0);
    phgMapP2C(e, vmap1, emap1, 1);

#if USE_MPI
    if (dof->type->np_vert > 0) {
	if (e0 == NULL)
	    dof->type->InitFunc(dof, e, VERTEX, 0, NULL, f2c_func, NULL,
				parent_data[0], parent_data);
	else if (e1 == NULL)
	    dof->type->InitFunc(dof, e, VERTEX, 1, NULL, f2c_func, NULL,
				parent_data[1], parent_data);
    }
#endif

    if (dof->type->np_edge > 0) {
	dof->type->InitFunc(dof, e, EDGE, 0, NULL, f2c_func, NULL,
			    parent_data[NVert + 0], parent_data);
#if USE_MPI
	if (e0 == NULL) {
	    dof->type->InitFunc(dof, e, EDGE, 1, NULL, f2c_func, NULL,
				parent_data[NVert + 1], parent_data);
	    dof->type->InitFunc(dof, e, EDGE, 2, NULL, f2c_func, NULL,
				parent_data[NVert + 2], parent_data);
	}
	else if (e1 == NULL) {
	    dof->type->InitFunc(dof, e, EDGE, 3, NULL, f2c_func, NULL,
				parent_data[NVert + 3], parent_data);
	    dof->type->InitFunc(dof, e, EDGE, 4, NULL, f2c_func, NULL,
				parent_data[NVert + 4], parent_data);
	}
#endif
    }

    if (dof->type->np_face > 0) {
	dof->type->InitFunc(dof, e, FACE, 2, NULL, f2c_func, NULL,
			    parent_data[NVert + NEdge + 2], parent_data);
	dof->type->InitFunc(dof, e, FACE, 3, NULL, f2c_func, NULL,
			    parent_data[NVert + NEdge + 3], parent_data);
#if USE_MPI
	if (e0 == NULL)
	    dof->type->InitFunc(dof, e, FACE, 1, NULL, f2c_func, NULL,
				parent_data[NVert + NEdge + 1], parent_data);
	else if (e1 == NULL)
	    dof->type->InitFunc(dof, e, FACE, 0, NULL, f2c_func, NULL,
				parent_data[NVert + NEdge + 0], parent_data);
#endif
    }

    if (dof->type->np_elem > 0) {
	dof->type->InitFunc(dof, e, VOLUME, 0, NULL, f2c_func, NULL,
			    parent_data[NVert + NEdge + NFace], parent_data);
    }
}

void
phgDofInitFuncPoint(DOF *dof, ELEMENT *e, GTYPE type, int index,
		    DOF_USER_FUNC userfunc,
		    DOF_USER_FUNC_LAMBDA userfunc_lambda,
		    const FLOAT *funcvalues, FLOAT *dofvalues,
		    FLOAT **pdofvalues)
{
    const GRID *g = dof->g;
    COORD *p0 = NULL, *p1 = NULL, *p2 = NULL, *p3;
    FLOAT *points;
    FLOAT x, y, z, lambda[Dim + 1];
    int n, j, v0, v1, v2, bno = 0;
    DOF_TYPE *t;

    Unused(pdofvalues);
    CheckD_order(dof);

    switch (type) {
	case VERTEX:
	    t = DofTypeOnVertex(dof, e->verts[index]);
	    n = t->np_vert * t->dim;
	    break;
	case EDGE:
	    t = DofTypeOnEdge(dof, e->edges[index]);
	    n = t->np_edge * t->dim;
	    break;
	case FACE:
	    t = DofTypeOnFace(dof, e->faces[index]);
	    n = t->np_face * t->dim;
	    break;
	case VOLUME:
	    t = DofTypeOnElement(dof, e->index);
	    n = t->np_elem * t->dim;
	    break;
	default:
	    t = NULL;
	    n = 0;		/* make gcc happy */
    }

    if (n == 0)
	return;

    if (funcvalues != NULL) {
	memcpy(dofvalues, funcvalues, n * dof->dim * sizeof(*dofvalues));
	return;
    }

    if (t->points == NULL)
	phgError(1, "%s:%d, not a point DOF.\n", __FILE__, __LINE__);

    if (userfunc == NULL && userfunc_lambda == NULL)
	phgError(1, "%s:%d, invalid user function!\n", __FILE__, __LINE__);

    switch (type) {
	case VERTEX:
	    if (userfunc_lambda != NULL) {
		lambda[0] = lambda[1] = lambda[2] = lambda[3] = 0.;
		lambda[index] = 1.0;
		bno = index * n;
		for (j = 0; j < n; j++) {
		    userfunc_lambda(dof, e, bno + j, lambda, dofvalues);
		    dofvalues += dof->dim;
		}
		break;
	    }
	    x = g->verts[e->verts[index]][0];
	    y = g->verts[e->verts[index]][1];
	    z = g->verts[e->verts[index]][2];
	    for (j = 0; j < n; j++) {
		userfunc(x, y, z, dofvalues);
		dofvalues += dof->dim;
	    }
	    break;
	case EDGE:
	    v0 = GetEdgeVertex(index, 0);
	    v1 = GetEdgeVertex(index, 1);
	    if (userfunc_lambda != NULL) {
		lambda[0] = lambda[1] = lambda[2] = lambda[3] = 0.;
		bno = NVert * t->np_vert + index * n;
	    }
	    else {
		p0 = g->verts + e->verts[v0];
		p1 = g->verts + e->verts[v1];
	    }
	    points = t->points + t->np_vert;
	    if (n > 1) {
		int M[n];
		phgDofMapEdgeData(t, e, index, M);
		if (userfunc_lambda == NULL) {
		    for (j = 0; j < n; j++) {
			x = (*p0)[0] * points[0] + (*p1)[0] * points[1];
			y = (*p0)[1] * points[0] + (*p1)[1] * points[1];
			z = (*p0)[2] * points[0] + (*p1)[2] * points[1];
			userfunc(x, y, z, dofvalues + M[j] * dof->dim);
			points += 2;
		    }
		}
		else {
		    for (j = 0; j < n; j++) {
			lambda[v0] = points[0];
			lambda[v1] = points[1];
			userfunc_lambda(dof, e, bno + j, lambda,
					dofvalues + M[j] * dof->dim);
			points += 2;
		    }
		}
	    }
	    else {
		if (userfunc_lambda == NULL) {
		    for (j = 0; j < n; j++) {
			x = (*p0)[0] * points[0] + (*p1)[0] * points[1];
			y = (*p0)[1] * points[0] + (*p1)[1] * points[1];
			z = (*p0)[2] * points[0] + (*p1)[2] * points[1];
			userfunc(x, y, z, dofvalues);
			points += 2;
			dofvalues += dof->dim;
		    }
		}
		else {
		    for (j = 0; j < n; j++) {
			lambda[v0] = points[0];
			lambda[v1] = points[1];
			userfunc_lambda(dof, e, bno + j, lambda, dofvalues);
			points += 2;
			dofvalues += dof->dim;
		    }
		}
	    }
	    break;
	case FACE:
	    v0 = GetFaceVertex(index, 0);
	    v1 = GetFaceVertex(index, 1);
	    v2 = GetFaceVertex(index, 2);
	    if (userfunc_lambda != NULL) {
		lambda[0] = lambda[1] = lambda[2] = lambda[3] = 0.;
		bno = NVert * t->np_vert + NEdge * t->np_edge + index * n;
	    }
	    else {
		p0 = g->verts + e->verts[v0];
		p1 = g->verts + e->verts[v1];
		p2 = g->verts + e->verts[v2];
	    }
	    points = t->points + t->np_vert + t->np_edge * 2;
	    if (n > 1) {
		int M[n];
		phgDofMapFaceData(t, e, index, M);
		if (userfunc_lambda == NULL) {
		    for (j = 0; j < n; j++) {
			x = (*p0)[0] * points[0] + (*p1)[0] * points[1] +
			    (*p2)[0] * points[2];
			y = (*p0)[1] * points[0] + (*p1)[1] * points[1] +
			    (*p2)[1] * points[2];
			z = (*p0)[2] * points[0] + (*p1)[2] * points[1] +
			    (*p2)[2] * points[2];
			userfunc(x, y, z, dofvalues + M[j] * dof->dim);
			points += 3;
		    }
		}
		else {
		    for (j = 0; j < n; j++) {
			lambda[v0] = points[0];
			lambda[v1] = points[1];
			lambda[v2] = points[2];
			userfunc_lambda(dof, e, bno + j, lambda,
					dofvalues + M[j] * dof->dim);
			points += 3;
		    }
		}
	    }
	    else {
		if (userfunc_lambda == NULL) {
		    for (j = 0; j < n; j++) {
			x = (*p0)[0] * points[0] + (*p1)[0] * points[1] +
			    (*p2)[0] * points[2];
			y = (*p0)[1] * points[0] + (*p1)[1] * points[1] +
			    (*p2)[1] * points[2];
			z = (*p0)[2] * points[0] + (*p1)[2] * points[1] +
			    (*p2)[2] * points[2];
			userfunc(x, y, z, dofvalues);
			points += 3;
			dofvalues += dof->dim;
		    }
		}
		else {
		    for (j = 0; j < n; j++) {
			lambda[v0] = points[0];
			lambda[v1] = points[1];
			lambda[v2] = points[2];
			userfunc_lambda(dof, e, bno + j, lambda, dofvalues);
			points += 3;
			dofvalues += dof->dim;
		    }
		}
	    }
	    break;
	case VOLUME:
	    points = t->points + t->np_vert + t->np_edge * 2 + t->np_face * 3;
	    if (n > 1 && !ElementIsInOrder(e)) {
		int M[n];
		phgDofMapElementData(t, e, M);
		if (userfunc_lambda == NULL) {
		    p0 = g->verts + e->verts[0];
		    p1 = g->verts + e->verts[1];
		    p2 = g->verts + e->verts[2];
		    p3 = g->verts + e->verts[3];
		    for (j = 0; j < n; j++) {
			x = (*p0)[0] * points[0] + (*p1)[0] * points[1] +
			    (*p2)[0] * points[2] + (*p3)[0] * points[3];
			y = (*p0)[1] * points[0] + (*p1)[1] * points[1] +
			    (*p2)[1] * points[2] + (*p3)[1] * points[3];
			z = (*p0)[2] * points[0] + (*p1)[2] * points[1] +
			    (*p2)[2] * points[2] + (*p3)[2] * points[3];
			userfunc(x, y, z, dofvalues + M[j] * dof->dim);
			points += 4;
		    }
		}
		else {
		    bno = NVert * t->np_vert + NEdge * t->np_edge +
				NFace * t->np_face + index * n;
		    for (j = 0; j < n; j++) {
			userfunc_lambda(dof, e, bno + j, points,
						dofvalues + M[j] * dof->dim);
			points += 4;
		    }
		}
	    }
	    else {
		if (userfunc_lambda == NULL) {
		    p0 = g->verts + e->verts[0];
		    p1 = g->verts + e->verts[1];
		    p2 = g->verts + e->verts[2];
		    p3 = g->verts + e->verts[3];
		    for (j = 0; j < n; j++) {
			x = (*p0)[0] * points[0] + (*p1)[0] * points[1] +
			    (*p2)[0] * points[2] + (*p3)[0] * points[3];
			y = (*p0)[1] * points[0] + (*p1)[1] * points[1] +
			    (*p2)[1] * points[2] + (*p3)[1] * points[3];
			z = (*p0)[2] * points[0] + (*p1)[2] * points[1] +
			    (*p2)[2] * points[2] + (*p3)[2] * points[3];
			userfunc(x, y, z, dofvalues);
			points += 4;
			dofvalues += dof->dim;
		    }
		}
		else {
		    bno = NVert * t->np_vert + NEdge * t->np_edge +
				NFace * t->np_face + index * n;
		    for (j = 0; j < n; j++) {
			userfunc_lambda(dof, e, bno + j, points, dofvalues);
			points += 4;
			dofvalues += dof->dim;
		    }
		}
	    }
	    break;
    }
}

static BOOLEAN flagged_only = FALSE;

void
phgDofSetDataByFunction__(DOF *dof, DOF_USER_FUNC userfunc,
			  DOF_USER_FUNC_LAMBDA userfunc_lambda,
			  FLOAT *funcvalues, BTYPE bdry_mask,
			  FORALL_TYPE for_what)
/* Note:
 *  . 'bdry_mask' is used to only set boundary DOFs.
 *  . 'for_what' controls setting DOFs in which part of the submesh. */
{
    GRID *g = dof->g;
    INT i;
    size_t size;
    DOF_TYPE *type, *t;
    VEF_MAP *vef;

    if (dof->type == DOF_ANALYTIC)
	return;

    phgDofClearCache(NULL, dof, NULL, NULL, FALSE);

    size = (dof->type == DOF_CONSTANT) ? dof->dim : DofDataCount(dof);
    if (dof->data == NULL) {
	dof->data = phgCalloc(size, sizeof(*dof->data));
	dof->data_vert = dof->data;
	dof->data_edge = dof->data_vert + DofVertexDataCount(dof);
	dof->data_face = dof->data_edge + DofEdgeDataCount(dof);
	dof->data_elem = dof->data_face + DofFaceDataCount(dof);
    }

    if (userfunc == DofNoAction ||
	(userfunc == DofInterpolation && userfunc_lambda == NULL))
	return;

    if (dof->type == DOF_CONSTANT) {
	/* constant DOF, using funcvalues */
	if (funcvalues != NULL) {
	    memcpy(dof->data, funcvalues, dof->dim * sizeof(*dof->data));
	}
	else if (userfunc_lambda != NULL) {
	    FLOAT lambda[] = { 0., 0., 0., 1. };
	    userfunc_lambda(dof, dof->g->roots, -1, lambda, dof->data);
	}
	else {
	    assert(userfunc != NULL);
	    userfunc(0., 0., 0., dof->data);
	}
	return;
    }

    type = dof->type;
    if (type == NULL)
	type = dof->hp->info->types[dof->hp->max_order];
    if ((userfunc != NULL || userfunc_lambda != NULL)
	&& type->InitFunc == NULL)
	phgError(1, "%s: no initialization function!\n", __func__);

    vef = phgDofSetupVEFMap_(g, dof, -1, for_what);

    if (type->np_vert > 0) {
#if USE_OMP
# pragma omp parallel for schedule(static) private(i, t)
#endif	/* USE_OMP */
	for (i = 0; i < g->nvert; i++) {
	    if (!(vef->Vmap[i] != NULL
		  && (!flagged_only || vef->Vmap[i]->flag < 1)
		  && (g->types_vert[i] & bdry_mask)))
		continue;
	    t = DofTypeOnVertex(dof, i);
	    t->InitFunc(dof, vef->Vmap[i], VERTEX, vef->Vind[i],
				userfunc, userfunc_lambda, funcvalues,
				DofVertexData(dof, i), NULL);
	}
    }

    if (type->np_edge > 0) {
#if USE_OMP
# pragma omp parallel for schedule(static) private(i, t)
#endif	/* USE_OMP */
	for (i = 0; i < g->nedge; i++) {
	    if (!(vef->Emap[i] != NULL
		  && (!flagged_only || vef->Emap[i]->flag < 1)
		  && (g->types_edge[i] & bdry_mask)))
		continue;
	    t = DofTypeOnEdge(dof, i);
	    t->InitFunc(dof, vef->Emap[i], EDGE, vef->Eind[i],
			userfunc, userfunc_lambda, funcvalues,
			DofEdgeData(dof, i), NULL);
	}
    }

    if (type->np_face > 0) {
#if USE_OMP
# pragma omp parallel for schedule(static) private(i, t)
#endif	/* USE_OMP */
	for (i = 0; i < g->nface; i++) {
	    if (!(vef->Fmap[i] != NULL
		  && (!flagged_only || vef->Fmap[i]->flag < 1)
		  && (g->types_face[i] & bdry_mask)))
		continue;
	    t = DofTypeOnFace(dof, i);
	    t->InitFunc(dof, vef->Fmap[i], FACE, vef->Find[i],
			userfunc, userfunc_lambda, funcvalues,
			DofFaceData(dof, i), NULL);
	}
    }

    if (type->np_elem > 0) {
	ELEMENT *e;
#if USE_OMP
# pragma omp parallel for schedule(static) private(e, i, t)
#endif	/* USE_OMP */
	ForAllElementsBegin_(g, e, for_what) {
	    i = e->index;
	    if (!((!flagged_only || e->flag < 1) &&
		  (g->types_elem[i] & bdry_mask)))
		continue;
	    t = DofTypeOnElement(dof, i);
	    t->InitFunc(dof, e, VOLUME, 0,
			userfunc, userfunc_lambda, funcvalues,
			DofElementData(dof, i), NULL);
	} ForAllElementsEnd
    }

    phgDofFreeVEFMap(&vef);
}

static int dummy_n;
static const FLOAT *dummy_values;
static void
dummy_userfunc_lambda(DOF *dof, ELEMENT *e, int bno, const FLOAT *lambda,
			FLOAT *values)
{
    int i, n;
    const FLOAT *p;

    Unused(e);
    Unused(lambda);
    Unused(bno);

    /* Note: if dof->type->D_order > 0 then the derivatives are also required */
    CheckD_order(dof);
    n = DofDim(dof);
    for (i = 0, p = dummy_values; i < n; i++) {
	*(values++) = *(p++);
	if (p >= dummy_values + dummy_n)
	    p = dummy_values;
    }
}

void
phgDofSetDataByValue(DOF *dof, FLOAT value)
{
    size_t i, np, nvalues;
    FLOAT *values;

    assert(dof->type != DOF_ANALYTIC);

    if (dof->type != DOF_CONSTANT) {
	DOF_TYPE *type = dof->type;
	if (type == NULL)
	    type = dof->hp->info->types[dof->hp->max_order];
	nvalues = DofDim(dof);
	if (type->InitFunc != phgDofInitFuncPoint) {
	    if (type->InitFunc == NULL) {
		/* no init function, simply copy the data */
		for (i = 0; i < DofDataCount(dof); i++)
		    dof->data[i] = value;
		return;
	    }
	    dummy_n = 1;
	    dummy_values = &value;
	    phgDofSetDataByFunction_(dof, NULL, dummy_userfunc_lambda, NULL);
	    return;
	}
	np = type->np_vert;
	if (np < (i = type->np_edge))
	    np = i;
	if (np < (i = type->np_face))
	    np = i;
	if (np < (i = type->np_elem + 1))
	    np = i;
    }
    else {
	nvalues = dof->dim;
	np = 1;
    }
    values = phgAlloc(nvalues * np * sizeof(*values));
    /* make 'np' * 'nvalues' duplicate copies of the value */
    for (i = 0; i < np * nvalues; i++)
	values[i] = value;
    phgDofSetDataByFunction_(dof, NULL, NULL, values);
    phgFree(values);
}

void
phgDofSetDataByValues(DOF *dof, const FLOAT array_of_values[])
{
    size_t i, np, nvalues;
    FLOAT *values;

    assert(dof->type != DOF_ANALYTIC);

    if (dof->type != DOF_CONSTANT) {
	nvalues = DofDim(dof);
	if (dof->type->InitFunc != phgDofInitFuncPoint) {
	    dummy_n = nvalues;
	    dummy_values = array_of_values;
	    phgDofSetDataByFunction_(dof, NULL, dummy_userfunc_lambda, NULL);
	    return;
	}
	np = dof->type->np_vert;
	if (np < (i = dof->type->np_edge))
	    np = i;
	if (np < (i = dof->type->np_face))
	    np = i;
	if (np < (i = dof->type->np_elem))
	    np = i;
    }
    else {
	nvalues = dof->dim;
	np = 1;
    }
    values = phgAlloc(nvalues * np * sizeof(*values));
    /* make 'np' duplicate copies of the values */
    for (i = 0; i < nvalues; i++)
	values[i] = (i < nvalues) ? array_of_values[i] : values[i - nvalues];
    phgDofSetDataByFunction_(dof, NULL, NULL, values);
    phgFree(values);
}

void
phgDofSetDataByValuesV(DOF *dof, FLOAT v0, ...)
/* set DOF data using variable arguments */
{
    int i, n = DofDim(dof);
    FLOAT *values;
    va_list ap;                                                         

    values = phgAlloc(n * sizeof(*values));
    values[0] = v0;
    va_start(ap, v0);
    for (i = 1; i < n; i++) {
#if FT_PHG == FT_FLOAT
	/* Note: 'float' in variable list automatically promoted to 'double' */
	values[i] = (FLOAT)va_arg(ap, double);
#else	/* FT_PHG == FT_FLOAT */
	values[i] = va_arg(ap, FLOAT);
#endif	/* FT_PHG == FT_FLOAT */
    }
    va_end(ap);
    phgDofSetDataByValues(dof, values);
    phgFree(values);
}

static int level = 0;		/* recursion level */
/* common values and pointers used by the interpolation routines */
static int n0, n1, n2, n3;
static FLOAT *data0, *data1, *data2, *data3;

void
update_dof_data(DOF *dof, ELEMENT *e, FLOAT **data)
/* places DOF data for element 'e' pointed by 'data' to dof->data */
{
    int i, j;
    int size[4];
    int count[] = { NVert, NEdge, NFace, 1 };
    INT *index[4];
    FLOAT *loc[4];
    FLOAT *p;

    size[0] = n0;
    size[1] = n1;
    size[2] = n2;
    size[3] = n3;
    index[0] = e->verts;
    index[1] = e->edges;
    index[2] = e->faces;
    index[3] = &e->index;
    loc[0] = data0;
    loc[1] = data1;
    loc[2] = data2;
    loc[3] = data3;

    for (j = 0; j < 4; j++) {
	for (i = 0; i < count[j]; i++) {
	    p = loc[j] + index[j][i] * size[j];
	    if (size[j] > 0 && p != *data) {
		memcpy(p, *data, size[j] * sizeof(FLOAT));
	    }
	    data++;
	}
    }
}

static void
build_child_pointers(DOF *dof, ELEMENT *e, FLOAT **old_data, FLOAT **new_data,
		     FLOAT *buffer, int child)
/* construct list of pointers for children */
{
    int i, Vmap[NVert], Emap[NEdge];
    FLOAT *p1, *p2, *p3;

    phgMapC2P(e, Vmap, Emap, child);

    p1 = buffer;
    p2 = p1 + n1 * 4;		/* 4 new edges */
    p3 = p2 + n2 * 5;		/* 5 new faces */

    /* vertices */
    if (n0 > 0) {
	for (i = 0; i < NVert - 1; i++) {
	    *(new_data++) = old_data[Vmap[i]];
	}
	*(new_data++) = data0 + e->children[0]->verts[NVert - 1] * n0;
    }
    else {
	new_data += NVert;
    }

    /* edges */
    if (n1 > 0) {
	for (i = 0; i < NEdge; i++) {
	    switch (Emap[i]) {
		case -1:	/* new edge */
		    (*new_data++) = p1 + (Vmap[GetEdgeVertex(i, 0)]) * n1;
		    break;
		case 0:	/* cut edge */
		    (*new_data++) = p1 + child * n1;
		    break;
		default:	/* old edge */
		    (*new_data++) = old_data[NVert + Emap[i]];
		    break;
	    }
	}
    }
    else {
	new_data += NEdge;
    }

    /* faces */
    if (n2 > 0) {
	for (i = 0; i < NFace; i++) {
	    switch (Vmap[i]) {
		case 0:
		case 1:	/* new face */
		    (*new_data++) = p2 + 4 * n2;
		    break;
		case 2:	/* cut face on F2 */
		    (*new_data++) = p2 + child * n2;
		    break;
		case 3:	/* cut face on F3 */
		    (*new_data++) = p2 + (child + 2) * n2;
		    break;
		case -1:	/* old faces */
		    (*new_data++) = old_data[NVert + NEdge + 1 - child];
		    break;
	    }
	}
    }
    else {
	new_data += NFace;
    }

    /* volume */
    if (n3 > 0) {
	(*new_data) = p3 + child * n3;
    }
}

static void
interp_recursive(DOF *dof, ELEMENT *e, FLOAT **old_data)
/* recursively interpolates DOFs down the binary tree starting from e */
{
    static FLOAT *p1, *p2, *p3;
    static ELEMENT *e0, *e1;
    FLOAT *buffer, *new_data[15];

    if ((e0 = e->children[0]) == NULL || (e1 = e->children[1]) == NULL) {
	/* copy data to dof->data */
	update_dof_data(dof, e, old_data);
	return;
    }

    level++;

    /* allocate buffer for new DOFs (only allocate for edge/face/vol data) */
    if (n1 > 0 || n2 > 0 || n3 > 0)
	buffer = phgAlloc((n1 * 4 + n2 * 5 + n3 * 2) * sizeof(*buffer));
    else
	buffer = NULL;

    p1 = buffer;
    p2 = p1 + n1 * 4;
    p3 = p2 + n2 * 5;

    memset(new_data, 0, sizeof(new_data));

    if (dof->type->InterpC2F != NULL && dof->userfunc == DofInterpolation) {
	/* interpolation */
	if (n0 > 0) {
	    /* the new vertex */
	    new_data[0] = data0 + e0->verts[3] * n0;
	}

	if (n1 > 0) {
	    /* the cut edge */
	    new_data[1] = p1;
	    new_data[2] = p1 + n1;
	    /* the new edge on V2 */
	    new_data[3] = p1 + n1 * 2;
	    /* the new edge on V3 */
	    new_data[4] = p1 + n1 * 3;
	}

	if (n2 > 0) {
	    /* the cut faces on F2 */
	    new_data[5] = p2;
	    new_data[6] = p2 + n2;
	    /* the cut faces on F3 */
	    new_data[7] = p2 + n2 * 2;
	    new_data[8] = p2 + n2 * 3;
	    /* the new face */
	    new_data[9] = p2 + n2 * 4;
	}

	if (n3 > 0) {
	    /* volume data for child 0 */
	    new_data[10] = p3;
	    /* volume data for child 1 */
	    new_data[11] = p3 + n3;
	}

	dof->type->InterpC2F(dof, e, old_data, new_data);
    }

    if (e0 != NULL) {
	e0->flag = 0;		/* clear flag on new elements */
	build_child_pointers(dof, e, old_data, new_data, buffer, 0);
	interp_recursive(dof, e0, new_data);
	e1 = e->children[1];	/* note: e1 changed by recursion */
    }

    if (e1 != NULL) {
	e1->flag = 0;		/* clear flag on new elements */
	build_child_pointers(dof, e, old_data, new_data, buffer, 1);
	interp_recursive(dof, e1, new_data);
    }

    phgFree(buffer);

    level--;
}

static INT nelem = 0, nvert = 0, nedge = 0, nface = 0;
static BYTE elem_order, *interp_flags;
static DOF *u;
static HP_TYPE *old_hp, *new_hp;

static void
interp_hp_orders(ELEMENT *e,
		 BYTE vert_order[], BYTE edge_order[], BYTE face_order[])
/* compute orders on new edges/faces/elements */
{
    static int Vmap[NVert];
    static ELEMENT *ee;
    BYTE vert_order1[NVert], edge_order1[NEdge], face_order1[NFace];

    /* child 0 */
    phgMapP2C(e, Vmap, NULL, 0);

    vert_order1[Vmap[0]]	= vert_order[0];	/* vert 0 */
    vert_order1[Vmap[2]]	= vert_order[2];	/* vert 2 */
    vert_order1[Vmap[3]]	= vert_order[3];	/* vert 3 */
    vert_order1[3]		= edge_order[0];	/* new vertex */

    edge_order1[GetEdgeNo(3, Vmap[0])]		= edge_order[0]; /* edge 0-1 */
    edge_order1[GetEdgeNo(Vmap[0], Vmap[2])]	= edge_order[1]; /* edge 0-2 */
    edge_order1[GetEdgeNo(Vmap[0], Vmap[3])]	= edge_order[2]; /* edge 0-3 */
    edge_order1[GetEdgeNo(Vmap[2], Vmap[3])]	= edge_order[5]; /* edge 2-3 */
    edge_order1[GetEdgeNo(3, Vmap[2])]		= face_order[3]; /* face 3 */
    edge_order1[GetEdgeNo(3, Vmap[3])]		= face_order[2]; /* face 2 */

    face_order1[Vmap[0]]	= elem_order;		/* new face */
    face_order1[3]		= face_order[1];	/* face 1 */
    face_order1[Vmap[2]]	= face_order[2];	/* face 2 */
    face_order1[Vmap[3]]	= face_order[3];	/* face 3 */

    ee = e->children[0];
    if (!IsLeaf(ee)) {
	interp_hp_orders(ee, vert_order1, edge_order1, face_order1);
    }
    else {
	if (new_hp->info->flags & VERT_FLAG) {
	    new_hp->vert_order[ee->verts[0]] = vert_order1[0];
	    new_hp->vert_order[ee->verts[1]] = vert_order1[1];
	    new_hp->vert_order[ee->verts[2]] = vert_order1[2];
	    new_hp->vert_order[ee->verts[3]] = vert_order1[3];
	}
	if (new_hp->info->flags & EDGE_FLAG) {
	    new_hp->edge_order[ee->edges[0]] = edge_order1[0];
	    new_hp->edge_order[ee->edges[1]] = edge_order1[1];
	    new_hp->edge_order[ee->edges[2]] = edge_order1[2];
	    new_hp->edge_order[ee->edges[3]] = edge_order1[3];
	    new_hp->edge_order[ee->edges[4]] = edge_order1[4];
	    new_hp->edge_order[ee->edges[5]] = edge_order1[5];
	}
	if (new_hp->info->flags & FACE_FLAG) {
	    new_hp->face_order[ee->faces[0]] = face_order1[0];
	    new_hp->face_order[ee->faces[1]] = face_order1[1];
	    new_hp->face_order[ee->faces[2]] = face_order1[2];
	    new_hp->face_order[ee->faces[3]] = face_order1[3];
	}
	if (new_hp->info->flags & ELEM_FLAG)
	    new_hp->elem_order[ee->index] = elem_order;
    }

    /* child 1 */
    phgMapP2C(e, Vmap, NULL, 1);

    vert_order1[Vmap[1]]	= vert_order[1];	/* vert 1 */
    vert_order1[Vmap[2]]	= vert_order[2];	/* vert 2 */
    vert_order1[Vmap[3]]	= vert_order[3];	/* vert 3 */
    vert_order1[3]		= edge_order[0];	/* new vertex */

    edge_order1[GetEdgeNo(3, Vmap[1])]		= edge_order[0]; /* edge 0-1 */
    edge_order1[GetEdgeNo(Vmap[1], Vmap[2])]	= edge_order[3]; /* edge 1-2 */
    edge_order1[GetEdgeNo(Vmap[1], Vmap[3])]	= edge_order[4]; /* edge 1-3 */
    edge_order1[GetEdgeNo(Vmap[2], Vmap[3])]	= edge_order[5]; /* edge 2-3 */
    edge_order1[GetEdgeNo(3, Vmap[2])]		= face_order[3]; /* face 3 */
    edge_order1[GetEdgeNo(3, Vmap[3])]		= face_order[2]; /* face 2 */

    face_order1[Vmap[1]]	= elem_order;		/* new face */
    face_order1[3]		= face_order[0];	/* face 0 */
    face_order1[Vmap[2]]	= face_order[2];	/* face 2 */
    face_order1[Vmap[3]]	= face_order[3];	/* face 3 */

    ee = e->children[1];
    if (!IsLeaf(ee)) {
	interp_hp_orders(ee, vert_order1, edge_order1, face_order1);
    }
    else {
	if (new_hp->info->flags & VERT_FLAG) {
	    new_hp->vert_order[ee->verts[0]] = vert_order1[0];
	    new_hp->vert_order[ee->verts[1]] = vert_order1[1];
	    new_hp->vert_order[ee->verts[2]] = vert_order1[2];
	    new_hp->vert_order[ee->verts[3]] = vert_order1[3];
	}
	if (new_hp->info->flags & EDGE_FLAG) {
	    new_hp->edge_order[ee->edges[0]] = edge_order1[0];
	    new_hp->edge_order[ee->edges[1]] = edge_order1[1];
	    new_hp->edge_order[ee->edges[2]] = edge_order1[2];
	    new_hp->edge_order[ee->edges[3]] = edge_order1[3];
	    new_hp->edge_order[ee->edges[4]] = edge_order1[4];
	    new_hp->edge_order[ee->edges[5]] = edge_order1[5];
	}
	if (new_hp->info->flags & FACE_FLAG) {
	    new_hp->face_order[ee->faces[0]] = face_order1[0];
	    new_hp->face_order[ee->faces[1]] = face_order1[1];
	    new_hp->face_order[ee->faces[2]] = face_order1[2];
	    new_hp->face_order[ee->faces[3]] = face_order1[3];
	}
	if (new_hp->info->flags & ELEM_FLAG)
	    new_hp->elem_order[ee->index] = elem_order;
    }
}

static void
dof_interp_func(DOF *u, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    HP_TYPE *hp;
#if 0
    /* lambda on e ==> x,y,z ==> lambda on elem */
    FLOAT lam[Dim + 1];
    FLOAT x, y, z;

    Unused(bno);
    phgGeomLambda2XYZ(u->g, e, lambda, &x, &y, &z);
    phgGeomXYZ2Lambda(u->g, elem, x, y, z, lam);
#else	/* 0 */
    /* child ==> parent loop */
    int v[NVert], k;
    FLOAT lambda0[Dim + 1], lambda1[Dim + 1], *lam, *lam1;
    ELEMENT *ep;

    lam = (FLOAT *)lambda;
    lam1 = lambda1;
    while (e != elem) {
	ep = e->parent;
	if (e == ep->children[0]) {
	    phgMapP2C(ep, v, NULL, 0);
	    lam1[1] = 0.5 * lam[3];
	    lam1[0] = lam[v[0]] + lam1[1];
	    lam1[2] = lam[v[2]];
	    lam1[3] = lam[v[3]];
	}
	else {
	    phgMapP2C(ep, v, NULL, 1);
	    lam1[0] = 0.5 * lam[3];
	    lam1[1] = lam[v[1]] + lam1[0];
	    lam1[2] = lam[v[2]];
	    lam1[3] = lam[v[3]];
	}
	lam = lam1;
	lam1 = (lam == lambda1 ? lambda0 : lambda1);
	e = ep;
    }
#endif	/* 0 */

    hp = u->hp;		/* backup u->hp */
    u->hp = old_hp;	/* use old_hp */
    phgDofEval_(u, elem, lam, NULL, 0, values, cdata0);
    /* need \grad u if D_order > 0 */
    k = (u->type != NULL ? u->type->D_order : u->hp->info->types[0]->D_order);
    assert(k <= 1);
    if (k == 1)
	phgDofEvalGradient_(u, elem, lam, NULL, values + DofDim(u), cdata0);
    u->hp = hp;		/* restore u->hp */
}

static void
dof_interp0(ELEMENT *e)
/* compute DOF on element e.
 *	- 'elem' points to its ancestor in the original mesh
 *	- 'cdata0' points to data locations on 'elem' in the original data
 *	- 'old_hp' points to the old HP_TYPE struct
 *	- 'u' is the DOF with new data and hp pointers. */
{
    INT i;
    int ii;
    DOF_TYPE *t;

    if (!DofIsHP(u))
	t = u->type;
    else
	t = u->hp->info->types[u->hp->elem_order[e->index]];

    if (!DofIsHP(u)) {	/* non h-p DOF */
	if (u->count_vert > 0) {
	    for (ii = 0; ii < NVert; ii++) {
		if ((i = e->verts[ii]) < nvert || (interp_flags[i] & 1))
		    continue;
		interp_flags[i] |= 1;
		t->InitFunc(u, e, VERTEX, ii, NULL, dof_interp_func,
				  NULL, u->data_vert + i * u->count_vert, NULL);
	    }
	}

	if (u->count_edge > 0) {
	    for (ii = 0; ii < NEdge; ii++) {
		if ((interp_flags[i = e->edges[ii]] & 2))
		    continue;
		interp_flags[i] |= 2;
		t->InitFunc(u, e, EDGE, ii, NULL, dof_interp_func,
				  NULL, u->data_edge + i * u->count_edge, NULL);
	    }
	}

	if (u->count_face > 0) {
	    for (ii = 0; ii < NFace; ii++) {
		if ((interp_flags[i = e->faces[ii]] & 4))
		    continue;
		interp_flags[i] |= 4;
		t->InitFunc(u, e, FACE, ii, NULL, dof_interp_func,
				  NULL, u->data_face + i * u->count_face, NULL);
	    }
	}

	if (u->count_elem > 0) {
	    i = e->index;
	    t->InitFunc(u, e, VOLUME, 0, NULL, dof_interp_func,
			      NULL, u->data_elem + i * u->count_elem, NULL);
	}
    }
    else {			/* h-p DOF */
	if (u->count_vert > 0) {
	    for (ii = 0; ii < NVert; ii++) {
		if ((i = e->verts[ii]) < nvert || (interp_flags[i] & 1))
		    continue;
		interp_flags[i] |= 1;
		DofTypeOnVertex(u, e->verts[ii])->InitFunc(
				u, e, VERTEX, ii, NULL, dof_interp_func, NULL,
				u->data_vert + u->hp->vert_index[i] * u->dim,
				NULL);
	    }
	}

	if (t->np_edge > 0) {
	    for (ii = 0; ii < NEdge; ii++) {
		if ((interp_flags[i = e->edges[ii]] & 2))
		    continue;
		interp_flags[i] |= 2;
		t->InitFunc(u, e, EDGE, ii, NULL, dof_interp_func, NULL,
			    u->data_edge + u->hp->edge_index[i] * u->dim,
			    NULL);
	    }
	}

	if (t->np_face > 0) {
	    for (ii = 0; ii < NFace; ii++) {
		if ((interp_flags[i = e->faces[ii]] & 4))
		    continue;
		interp_flags[i] |= 4;
		t->InitFunc(u, e, FACE, ii, NULL, dof_interp_func, NULL,
			    u->data_face + u->hp->face_index[i] * u->dim,
			    NULL);
	    }
	}

	if (t->np_elem > 0) {
	    i = e->index;
	    t->InitFunc(u, e, VOLUME, 0, NULL, dof_interp_func, NULL,
			u->data_elem + u->hp->elem_index[i] * u->dim,
			NULL);
	}
    }
}

static void
dof_interp(ELEMENT *e)
{
    if (e->children[0] == NULL) {
	/* leaf element */
	dof_interp0(e);
    }
    else {
	/* non leaf element */
	dof_interp(e->children[0]);
	dof_interp(e->children[1]);
    }
}

static inline void
interp_clear_flag(ELEMENT *e)
{
    if (e->children[0] == NULL)
	return;

    e->children[0]->flag = 0;
    interp_clear_flag(e->children[0]);

    e->children[1]->flag = 0;
    interp_clear_flag(e->children[1]);
}

/* For non h-p FE DOFs, there're two options for interpolation:
 * 	1. recursive interpolation using the InterpF2C function.
 * 	2. direct evaluation of finite element functions.
 * For h-p DOFs, since the recursive version is not implemented, direct
 * evaluation is always used.
 *
 * INTERP_DIRECT_EVAL controls whether use recursion (0) or evaluation (1)
 * when interpolating non h-p DOFs (1 is generally faster than 0) */
#define INTERP_DIRECT_EVAL	1

void
phgDofInterpolate(GRID *g, BOOLEAN flag)
/* interpolates DOFs, should be called with flag==FALSE
   before refinement, then with flag==TRUE after refinement. */
{
    static ELEMENT **elems = NULL;
    static BOOLEAN elems_flag = FALSE;
    static INT nleaf;
    DOF **dof;
    double time0;
    ELEMENT **pe;
    size_t size, flag_size = 0;
    FLOAT *old_data, *p0, *p1, *p2, *p3;
    int ndof;	/* number of non special DOFs which have data */
    int nhp;	/* number of hp DOFs */
    int ii;
    HP_TYPE **oldhps;

    /* check whether interpolation is really needed */
    if (g->dof == NULL && g->hp == NULL)
	return;

    ndof = 0;
    for (dof = g->dof; *dof != NULL; dof++) {
	if (SpecialDofType((*dof)->type))
	    continue;
	ndof++;
    }

    nhp = 0;
    for (oldhps = g->hp; oldhps != NULL && *oldhps != NULL; oldhps++, nhp++);

    if (ndof == 0 && nhp == 0)
	return;

    if ((!flag && elems_flag) || (flag && !elems_flag))
	phgError(1, "phgDofInterpolate: incorrect calling sequence.\n");

    if (!flag) {
	ELEMENT *e;

	/* save current DOF counts */
	nvert = g->nvert;
	nedge = g->nedge;
	nface = g->nface;
	nelem = g->nelem;
	nleaf = g->nleaf;

	/* store current leaf elements */
	pe = elems = phgAlloc(nleaf * sizeof(*elems));
	ForAllElements(g, e)
	    *(pe++) = e;
	assert(pe - elems == nleaf);

	elems_flag = TRUE;
	return;
    }

    /* reset cache, map, etc. (not necessary since elements not changed) */
    /* phgDofClearCache(g, NULL, NULL, NULL, FALSE); */

    time0 = phgGetTime(NULL);

    if (g->nleaf == nleaf)
	goto end;		/* no new elements created */

    /* save old geom data */
    phgGeomSaveNonLeafData_(g);

    /* Note:
     * 
     * During interpolation of the geometric data which is the first DOF,
     * e->flag will be set to 1 for old elements and 0 for new elements,
     * which is needed for both setting DOF values with phgDofSetDataByFunction
     * (flagged_only == TRUE) and the phgGeomXXXX functions.
     *
     * When g->geom == NULL, e->flag needs to be explicitly set
     */
    if (g->geom == NULL) {
	for (pe = elems; pe < elems + nleaf; pe++) {
	    interp_clear_flag(*pe);
	    (*pe)->flag = 1;
	}
    }

    /* set up hp-orders for new elements, save old orders in oldhps */
    oldhps = phgAlloc(nhp * sizeof(*oldhps));
    for (ii = 0; ii < nhp; ii++) {
	BYTE vert_order[NVert], edge_order[NEdge], face_order[NFace];
	new_hp = g->hp[ii];
	new_hp->id = ii;
	old_hp = oldhps[ii] = phgHPDup__(g->hp[ii], TRUE);
	/* Note: vert/edge/face/elem numbers are preserved during refinement */
	if (old_hp->info->flags & VERT_FLAG) {
	    phgFree(new_hp->vert_order);
	    new_hp->vert_order = phgCalloc(g->nvert,
					sizeof(*new_hp->vert_order));
	    memcpy(new_hp->vert_order, old_hp->vert_order,
					nvert * sizeof(*old_hp->vert_order));
	}
	if (old_hp->info->flags & EDGE_FLAG) {
	    phgFree(new_hp->edge_order);
	    new_hp->edge_order = phgCalloc(g->nedge,
					sizeof(*new_hp->edge_order));
	    memcpy(new_hp->edge_order, old_hp->edge_order,
					nedge * sizeof(*old_hp->edge_order));
	}
	if (old_hp->info->flags & FACE_FLAG) {
	    phgFree(new_hp->face_order);
	    new_hp->face_order = phgCalloc(g->nface,
					sizeof(*new_hp->face_order));
	    memcpy(new_hp->face_order, old_hp->face_order,
					nface * sizeof(*old_hp->face_order));
	}
	if (old_hp->info->flags & ELEM_FLAG) {
	    phgFree(new_hp->elem_order);
	    new_hp->elem_order = phgCalloc(g->nelem,
					sizeof(*new_hp->elem_order));
	    memcpy(new_hp->elem_order, old_hp->elem_order,
					nelem * sizeof(*old_hp->elem_order));
	}
	/* copy orders on new verts/edges/faces/elements */
	for (pe = elems; pe < elems + nleaf; pe++) {
	    if ((*pe)->children[0] == NULL)
		continue;
	    if (old_hp->info->flags & VERT_FLAG) {
		vert_order[0] = old_hp->vert_order[(*pe)->verts[0]];
		vert_order[1] = old_hp->vert_order[(*pe)->verts[1]];
		vert_order[2] = old_hp->vert_order[(*pe)->verts[2]];
		vert_order[3] = old_hp->vert_order[(*pe)->verts[3]];
	    }
	    else {
		memset(vert_order, 0, sizeof(vert_order));
	    }
	    if (old_hp->info->flags & EDGE_FLAG) {
		edge_order[0] = old_hp->edge_order[(*pe)->edges[0]];
		edge_order[1] = old_hp->edge_order[(*pe)->edges[1]];
		edge_order[2] = old_hp->edge_order[(*pe)->edges[2]];
		edge_order[3] = old_hp->edge_order[(*pe)->edges[3]];
		edge_order[4] = old_hp->edge_order[(*pe)->edges[4]];
		edge_order[5] = old_hp->edge_order[(*pe)->edges[5]];
	    }
	    else {
		memset(edge_order, 0, sizeof(edge_order));
	    }
	    if (old_hp->info->flags & FACE_FLAG) {
		face_order[0] = old_hp->face_order[(*pe)->faces[0]];
		face_order[1] = old_hp->face_order[(*pe)->faces[1]];
		face_order[2] = old_hp->face_order[(*pe)->faces[2]];
		face_order[3] = old_hp->face_order[(*pe)->faces[3]];
	    }
	    else {
		memset(face_order, 0, sizeof(face_order));
	    }
	    if (old_hp->info->flags & ELEM_FLAG)
		elem_order = old_hp->elem_order[(*pe)->index];
	    interp_hp_orders(*pe, vert_order, edge_order, face_order);
	}
	/* Note: the global counters will be updated together */
	phgHPSetupDataPointers_(new_hp, g->nvert, g->nedge, g->nface, g->nelem,
				TRUE, FALSE);
    }

    /* loop on DOFs */
    interp_flags = NULL;
    for (dof = g->dof; (u = *dof) != NULL; dof++) {
	if (SpecialDofType(u->type) || u->data == NULL)
	    continue;

	/* allocate new buffer */
	size = DofDataCount(u);

	old_data = u->data;
	u->data = phgCalloc(size, sizeof(*u->data));
	u->data_vert = u->data;
	u->data_edge = u->data_vert + DofVertexDataCount(u);
	u->data_face = u->data_edge + DofEdgeDataCount(u);
	u->data_elem = u->data_face + DofFaceDataCount(u);

	if (u->userfunc == DofNoAction) {
	    assert(u != g->geom);
	    phgFree(old_data);
	    continue;
	}

	/* copy old data to new buffer */
	p0 = old_data;

	if (!DofIsHP(u)) {
	    old_hp = NULL;
	    if ((size = u->count_vert * nvert) > 0) {
		memcpy(u->data_vert, p0, size * sizeof(*old_data));
		p0 += size;
	    }

	    if ((size = u->count_edge * nedge) > 0) {
		memcpy(u->data_edge, p0, size * sizeof(*old_data));
		p0 += size;
	    }

	    if ((size = u->count_face * nface) > 0) {
		memcpy(u->data_face, p0, size * sizeof(*old_data));
		p0 += size;
	    }

	    if ((size = u->count_elem * nelem) > 0)
		memcpy(u->data_elem, p0, size * sizeof(*old_data));
	}
	else {
	    old_hp = oldhps[u->hp->id];
	    if (old_hp->info->flags & VERT_FLAG) {
		assert(old_hp->vert_index[nvert] == u->hp->vert_index[nvert]);
		if ((size = u->dim * old_hp->vert_index[nvert]) > 0) {
		    memcpy(u->data_vert, p0, size * sizeof(*old_data));
		    p0 += size;
		}
	    }
	    if (old_hp->info->flags & EDGE_FLAG) {
		assert(old_hp->edge_index[nedge] == u->hp->edge_index[nedge]);
		if ((size = u->dim * old_hp->edge_index[nedge]) > 0) {
		    memcpy(u->data_edge, p0, size * sizeof(*old_data));
		    p0 += size;
		}
	    }
	    if (old_hp->info->flags & FACE_FLAG) {
		assert(old_hp->face_index[nface] == u->hp->face_index[nface]);
		if ((size = u->dim * old_hp->face_index[nface]) > 0) {
		    memcpy(u->data_face, p0, size * sizeof(*old_data));
		    p0 += size;
		}
	    }
	    if (old_hp->info->flags & ELEM_FLAG) {
		assert(old_hp->elem_index[nelem] == u->hp->elem_index[nelem]);
		if ((size = u->dim * old_hp->elem_index[nelem]) > 0) {
		    memcpy(u->data_elem, p0, size * sizeof(*old_data));
		}
	    }
	}

	if (u->userfunc != DofInterpolation) {
	    flagged_only = TRUE;
	    phgDofSetDataByFunction_(u, u->userfunc, u->userfunc_lambda, NULL);
	    flagged_only = FALSE;
	    assert(u != g->geom);
	    phgFree(old_data);
	    continue;
	}

	if (DofIsHP(u) /* h-p DOF */
#if INTERP_DIRECT_EVAL
	    || (u->type->BasFuncs != NULL && u->type->InitFunc != NULL)
#endif	/* INTERP_DIRECT_EVAL */
	   ) {
	    /* method 1: direct evaluation */
	    if (interp_flags == NULL) {
		flag_size = g->nface;
		if (flag_size < g->nedge)
		    flag_size = g->nedge;
		if (flag_size < g->nvert)
		    flag_size = g->nvert;
		interp_flags = phgAlloc(flag_size * sizeof(*interp_flags));
	    }
	    memset(interp_flags, 0, flag_size * sizeof(*interp_flags));
	    for (pe = elems; pe < elems + nleaf; pe++) {
		/* set flags on elements which are not refined */
		if ((elem = *pe)->children[0] != NULL)
		    continue;
		interp_flags[elem->verts[0]] |= 1;
		interp_flags[elem->verts[1]] |= 1;
		interp_flags[elem->verts[2]] |= 1;
		interp_flags[elem->verts[3]] |= 1;

		interp_flags[elem->edges[0]] |= 2;
		interp_flags[elem->edges[1]] |= 2;
		interp_flags[elem->edges[2]] |= 2;
		interp_flags[elem->edges[3]] |= 2;
		interp_flags[elem->edges[4]] |= 2;
		interp_flags[elem->edges[5]] |= 2;

		interp_flags[elem->faces[0]] |= 4;
		interp_flags[elem->faces[1]] |= 4;
		interp_flags[elem->faces[2]] |= 4;
		interp_flags[elem->faces[3]] |= 4;
	    }
	    for (pe = elems; pe < elems + nleaf; pe++) {
		/* interpolate on elements which are refined */
		if ((elem = *pe)->children[0] == NULL)
		    continue;
		data0 = old_data;
		if (!DofIsHP(u)) {
		    data1 = data0 + u->count_vert * nvert;
		    data2 = data1 + u->count_edge * nedge;
		    data3 = data2 + u->count_face * nface;

		    cdata0[0] = data0 + elem->verts[0] * u->count_vert;
		    cdata0[1] = data0 + elem->verts[1] * u->count_vert;
		    cdata0[2] = data0 + elem->verts[2] * u->count_vert;
		    cdata0[3] = data0 + elem->verts[3] * u->count_vert;

		    cdata0[4] = data1 + elem->edges[0] * u->count_edge;
		    cdata0[5] = data1 + elem->edges[1] * u->count_edge;
		    cdata0[6] = data1 + elem->edges[2] * u->count_edge;
		    cdata0[7] = data1 + elem->edges[3] * u->count_edge;
		    cdata0[8] = data1 + elem->edges[4] * u->count_edge;
		    cdata0[9] = data1 + elem->edges[5] * u->count_edge;

		    cdata0[10] = data2 + elem->faces[0] * u->count_face;
		    cdata0[11] = data2 + elem->faces[1] * u->count_face;
		    cdata0[12] = data2 + elem->faces[2] * u->count_face;
		    cdata0[13] = data2 + elem->faces[3] * u->count_face;

		    cdata0[14] = data3 + elem->index * u->count_elem;
		}
		else {
		    data1 = data0;
		    if (old_hp->info->flags & VERT_FLAG)
			data1 += old_hp->vert_index[nvert] * u->dim;
		    data2 = data1;
		    if (old_hp->info->flags & EDGE_FLAG)
			data2 += old_hp->edge_index[nedge] * u->dim;
		    data3 = data2;
		    if (old_hp->info->flags & FACE_FLAG)
			data3 += old_hp->face_index[nface] * u->dim;

		    cdata0[0] = cdata0[1] = cdata0[2] = cdata0[3] = data0;
		    if (old_hp->info->flags & VERT_FLAG) {
			cdata0[0] += old_hp->vert_index[elem->verts[0]]*u->dim;
			cdata0[1] += old_hp->vert_index[elem->verts[1]]*u->dim;
			cdata0[2] += old_hp->vert_index[elem->verts[2]]*u->dim;
			cdata0[3] += old_hp->vert_index[elem->verts[3]]*u->dim;
		    }

		    cdata0[4] = cdata0[5] = cdata0[6] = cdata0[7] =
				cdata0[8] = cdata0[9] = data1;
		    if (old_hp->info->flags & EDGE_FLAG) {
			cdata0[4] += old_hp->edge_index[elem->edges[0]]*u->dim;
			cdata0[5] += old_hp->edge_index[elem->edges[1]]*u->dim;
			cdata0[6] += old_hp->edge_index[elem->edges[2]]*u->dim;
			cdata0[7] += old_hp->edge_index[elem->edges[3]]*u->dim;
			cdata0[8] += old_hp->edge_index[elem->edges[4]]*u->dim;
			cdata0[9] += old_hp->edge_index[elem->edges[5]]*u->dim;
		    }

		    cdata0[10] = cdata0[11] = cdata0[12] = cdata0[13] = data2;
		    if (old_hp->info->flags & FACE_FLAG) {
			cdata0[10] += old_hp->face_index[elem->faces[0]]*u->dim;
			cdata0[11] += old_hp->face_index[elem->faces[1]]*u->dim;
			cdata0[12] += old_hp->face_index[elem->faces[2]]*u->dim;
			cdata0[13] += old_hp->face_index[elem->faces[3]]*u->dim;
		    }

		    cdata0[14] = data3;
		    if (old_hp->info->flags & ELEM_FLAG)
			cdata0[14] += old_hp->elem_index[elem->index] * u->dim;
		}
		dof_interp(elem);
	    }
	    assert(u != g->geom);
	    phgFree(old_data);
	    continue;
	}

	/* method 2: recursive interpolation using type->InterpC2F */
	assert(!DofIsHP(u));	/* h-p support unimplemented yet */
	n0 = u->count_vert;
	n1 = u->count_edge;
	n2 = u->count_face;
	n3 = u->count_elem;

	data0 = u->data;
	data1 = data0 + n0 * g->nvert;
	data2 = data1 + n1 * g->nedge;
	data3 = data2 + n2 * g->nface;

	p0 = old_data;
	p1 = p0 + n0 * nvert;
	p2 = p1 + n1 * nedge;
	p3 = p2 + n2 * nface;

	for (pe = elems; pe < elems + nleaf; pe++) {
	    FLOAT *old[15];
	    if (n0 > 0) {
		old[0] = p0 + (*pe)->verts[0] * n0;
		old[1] = p0 + (*pe)->verts[1] * n0;
		old[2] = p0 + (*pe)->verts[2] * n0;
		old[3] = p0 + (*pe)->verts[3] * n0;
	    }
	    if (n1 > 0) {
		old[4] = p1 + (*pe)->edges[0] * n1;
		old[5] = p1 + (*pe)->edges[1] * n1;
		old[6] = p1 + (*pe)->edges[2] * n1;
		old[7] = p1 + (*pe)->edges[3] * n1;
		old[8] = p1 + (*pe)->edges[4] * n1;
		old[9] = p1 + (*pe)->edges[5] * n1;
	    }
	    if (n2 > 0) {
		old[10] = p2 + (*pe)->faces[0] * n2;
		old[11] = p2 + (*pe)->faces[1] * n2;
		old[12] = p2 + (*pe)->faces[2] * n2;
		old[13] = p2 + (*pe)->faces[3] * n2;
	    }
	    if (n3 > 0) {
		old[14] = p3 + (*pe)->index * n3;
	    }

	    interp_recursive(u, *pe, old);
	    (*pe)->flag = 1;
	}
	if (u != g->geom)
	    phgFree(old_data);
    }
    phgFree(interp_flags);

    /* clear (free) old geom data */
    phgGeomClearNonLeafData_();

    /* free saved old HP_TYPEs */
    for (ii = 0; ii < nhp; ii++) {
	phgHPFreeBuffers__(oldhps[ii]);
	phgFree(oldhps[ii]);
    }
    phgFree(oldhps);

end:

    /* update global counters hp->xxxx_data_count_global */
    if (nhp > 0) {
	INT i, (*a)[4];
	a = phgAlloc(nhp * sizeof(*a));
	for (ii = 0; ii < nhp; ii++) {
	    new_hp = g->hp[ii];

	    nvert = 0;
	    for (i = 0; (new_hp->info->flags & VERT_FLAG) && i < g->nvert; i++)
		if ((g->types_vert[i] & OWNER))
		    nvert += new_hp->vert_index[i + 1] - new_hp->vert_index[i];

	    nedge = 0;
	    for (i = 0; (new_hp->info->flags & EDGE_FLAG) && i < g->nedge; i++)
		if ((g->types_edge[i] & OWNER))
		    nedge += new_hp->edge_index[i + 1] - new_hp->edge_index[i];

	    nface = 0;
	    for (i = 0; (new_hp->info->flags & FACE_FLAG) && i < g->nface; i++)
		if ((g->types_face[i] & OWNER))
		    nface += new_hp->face_index[i + 1] - new_hp->face_index[i];

	    nelem = 0;
	    for (i = 0; (new_hp->info->flags & ELEM_FLAG) && i < g->nelem; i++)
		if ((g->types_elem[i] & OWNER))
		    nelem += new_hp->elem_index[i + 1] - new_hp->elem_index[i];

	    a[ii][0] = nvert;
	    a[ii][1] = nedge;
	    a[ii][2] = nface;
	    a[ii][3] = nelem;
	}

#if USE_MPI
	if (g->nprocs > 1) {
	    INT (*b)[4];
	    b = phgAlloc(nhp * sizeof(*b));
	    MPI_Allreduce(a, b, 4 * nhp, PHG_MPI_INT, MPI_SUM, g->comm);
	    phgFree(a);
	    a = b;
	}
#endif	/* USE_MPI */

	for (ii = 0; ii < nhp; ii++) {
	    new_hp = g->hp[ii];
	    new_hp->vert_data_count_global = a[ii][0];
	    new_hp->edge_data_count_global = a[ii][1];
	    new_hp->face_data_count_global = a[ii][2];
	    new_hp->elem_data_count_global = a[ii][3];
	}

	phgFree(a);
    }

    /* Free memory */
    phgFree(elems);
    elems = NULL;
    elems_flag = FALSE;

    if (phgVerbosity >= 1)
	phgPrintf("* %s: time = %0.4lf\n", __func__, phgGetTime(NULL) - time0);
}

#if USE_MPI
static INT *current_map;
static int
sort_map_indices_comp(const void *p0, const void *p1)
{
    INT i = current_map[*((const INT *)p0)] - current_map[*((const INT *)p1)];
    return i > 0 ? 1 : (i < 0 ? -1 : 0);
}

static INT *
sort_map_indices(INT count, INT *map)
{
    INT i;

    if (count == 0)
	return NULL;

    current_map = map;
    map = phgAlloc(count * sizeof(*map));
    for (i = 0; i < count; i++)
	map[i] = i;
    if (current_map != NULL)
	qsort(map, count, sizeof(*map), sort_map_indices_comp);

    return map;
}

static void
copy_dof_data(DOF *dof,
	FLOAT *d0, INT n0[4], INT *index0[4], INT *map0[4], HP_TYPE *hp0,
	BTYPE *type0[4],
	FLOAT *d1, INT n1[4], INT *index1[4], INT *map1[4], HP_TYPE *hp1,
	BTYPE *type1[4])
/* copy data at 'd0' indexed with 'index0' and 'map0' to
	data at 'd1' indexed with 'index1' and 'map1'.

   Note: index0 or index1 is allowed to be NULL meaning identity map */
{
    INT *M0, *M1, *I0, *I1, *L0, *L1, N0, N1;
    BTYPE *T0, *T1;
    INT dof_size[4], *dindex0[4], *dindex1[4];
    int i0, i1, j0 = 0, j1, ii;
    size_t size;

    if (hp0 != NULL) {
	assert(hp1 != NULL);
	dindex0[0] = hp0->vert_index;
	dindex0[1] = hp0->edge_index;
	dindex0[2] = hp0->face_index;
	dindex0[3] = hp0->elem_index;
	dindex1[0] = hp1->vert_index;
	dindex1[1] = hp1->edge_index;
	dindex1[2] = hp1->face_index;
	dindex1[3] = hp1->elem_index;
    }
    else {
	dindex0[0] = dindex0[1] = dindex0[2] = dindex0[3] = NULL;
	dindex1[0] = dindex1[1] = dindex1[2] = dindex1[3] = NULL;
    }

    dof_size[0] = dof->count_vert;
    dof_size[1] = dof->count_edge;
    dof_size[2] = dof->count_face;
    dof_size[3] = dof->count_elem;

    for (ii = 0; ii < 4; ii++) {
	M0 = *(map0++);
	M1 = *(map1++);
	I0 = index0 == NULL ? NULL : *(index0++);
	I1 = index1 == NULL ? NULL : *(index1++);
	N0 = *(n0++);
	N1 = *(n1++);
	L0 = dindex0[ii];
	L1 = dindex1[ii];
	T0 = (type0 == NULL) ? NULL : *(type0++);
	T1 = (type1 == NULL) ? NULL : *(type1++);
	size = dof_size[ii];
	/* copy data */
	i1 = 0;
	i0 = 0;
	while (i1 < N1) {
	    j1 = (I1 == NULL ? i1 : I1[i1]);
	    if (T1 != NULL && T1[j1] == UNREFERENCED) {
		i1++;
		continue;
	    }
	    while (i0 < N0 && M0[j0 = (I0 == NULL ? i0 : I0[i0])] < M1[j1])
		i0++;
	    if (i0 >= N0)
		break;
	    if (M0[j0] > M1[j1] || (T0 != NULL && T0[j0] == UNREFERENCED)) {
		i1++;
		continue;
	    }
	    /* copy data */
	    if (L0 != NULL) {
		assert(L1 != NULL);
		size = dof->dim * (L1[j1 + 1] - L1[j1]);
		if (size == 0 || L0[j0 + 1] == L0[j0]) {
		    /* src or target refers to unreferenced location, skip */
		    i1++;
		    continue;
		}
		assert(size == dof->dim * (L0[j0 + 1] - L0[j0]));
		memcpy(d1 + L1[j1] * dof->dim, d0 + L0[j0] * dof->dim,
			size * sizeof(*d0));
	    }
	    else {
		assert(L1 == NULL);
		memcpy(d1 + j1 * size, d0 + j0 * size, size * sizeof(*d0));
	    }
	    i1++;
	}
	d0 += (L0 == NULL ? size * N0 : L0[N0] * dof->dim);
	d1 += (L1 == NULL ? size * N1 : L1[N1] * dof->dim);
    }
}

static void
copy_hp_orders(BYTE *orders0[4], INT n0[4], INT *index0[4], INT *map0[4],
	       BTYPE *type0[4],
	       BYTE *orders1[4], INT n1[4], INT *index1[4], INT *map1[4])
/* copy orders0 indexed with 'index0' and 'map0' to
	orders1 indexed with 'index1' and 'map1'.
   Note: index0 or index1 is allowed to be NULL meaning identity map */
{
    INT *M0, *M1, *I0, *I1, N0, N1;
    int i0, i1, j0 = 0, j1, ii;
    BTYPE *T0;

    for (ii = 0; ii < 4; ii++, orders0++, orders1++) {
	M0 = *(map0++);
	M1 = *(map1++);
	I0 = index0 == NULL ? NULL : *(index0++);
	I1 = index1 == NULL ? NULL : *(index1++);
	N0 = *(n0++);
	N1 = *(n1++);
	T0 = (type0 == NULL) ? NULL : *(type0++);
	if (*orders0 == NULL || *orders1 == NULL)
	    continue;	/* hp->xxxx_order == NULL */
	i1 = 0;
	i0 = 0;
	while (i1 < N1) {
	    j1 = (I1 == NULL ? i1 : I1[i1]);
	    while (i0 < N0 && M0[j0 = (I0 == NULL ? i0 : I0[i0])] < M1[j1])
		i0++;
	    if (i0 >= N0)
		break;
	    if (M0[j0] > M1[j1] || (T0 != NULL && T0[j0] == UNREFERENCED)) {
		i1++;
		continue;
	    }
	    /* copy order */
	    (*orders1)[j1] = (*orders0)[j0];
	    i1++;
	}
    }
}
#endif

#define IGNORE_NOACTION	0

void
phgDofRedistribute(GRID *g, BOOLEAN flag)
/* redistributes DOFs, should be called with flag==FALSE
 * before the grid is changed, then with flag==TRUE after.
 *
 * Note: some communications may be merged with those in phgRedistributeGrid
 * 	 for better efficiency.
 */
{
#if USE_MPI
    /* old maps */
    static INT	*L2Gmap_elem = NULL, *L2Gmap_vert = NULL,
		*L2Gmap_edge = NULL, *L2Gmap_face = NULL;
    /* old boundary types */
    static BTYPE *types_elem = NULL, *types_vert = NULL,
		 *types_edge = NULL, *types_face = NULL;

    /* lists recording vertices/edges/faces/elements to send,
     * The length of lists is g->nprocs. lists[i] describes data to be sent
     * to proc i, which is an array of INT organized as follows:
     *		lists[i][0] = number of vertices
     *		lists[i][1] = number of edges
     *		lists[i][2] = number of faces
     *		lists[i][3] = number of elements
     *
     *		lists[i][4] = h-p vertex data size
     *		lists[i][5] = h-p edge data size
     *		lists[i][6] = h-p face data size
     *		lists[i][7] = h-p element data size
     *
     * followed by:
     *		list of indices of vertices,
     * 		nhp_vert lists of vert orders, the length of
     * 			each list is lists[i][0] bytes, the last list is
     * 			padded to sizeof(INT).
     *
     *		list of indices of edges,
     * 		nhp_edge lists of edge orders, the length of
     * 			each list is lists[i][1] bytes, the last list is
     * 			padded to sizeof(INT).
     *
     *		list of indices of faces,
     * 		nhp_face lists of face orders, the length of
     * 			each list is lists[i][2] bytes, the last list is
     * 			padded to sizeof(INT).
     *
     *		list of indices of elements,
     * 		nhp_elem lists of elem orders, the length of
     * 			each list is lists[i][3] bytes, the last list is
     * 			padded to sizeof(INT).
     *
     * where nhp* are the count of HP_TYPEs.
     *
     * The total length of list[i] is:
     *	    8 + lists[i][0] +
     *	    lists[i][0] + (nhp_vert * lists[i][0] + sizeof(INT)-1)/sizeof(INT) +
     *	    lists[i][1] + (nhp_edge * lists[i][1] + sizeof(INT)-1)/sizeof(INT) +
     *	    lists[i][2] + (nhp_face * lists[i][2] + sizeof(INT)-1)/sizeof(INT) +
     *	    lists[i][3] + (nhp_elem * lists[i][3] + sizeof(INT)-1)/sizeof(INT)
     *		
     * The data block in sbuff and rbuff for process i consists of list[i][8:]
     * followed by DOF data in the order vertex data edge data, face data,
     * and element data, one DOF after another.
     */
    static INT **lists = NULL, *list;
    /* sorted old local indices */
    static INT *index0[4];

    /* blksize should be multiple of sizeof(INT) and sizeof(FLOAT),
     * and the maximum message size is limited by 2GB*blksize */
#if 0
    static size_t blksize = (sizeof(INT) > sizeof(FLOAT) ?
						sizeof(INT) : sizeof(FLOAT));
#else
    static size_t blksize = 128;
#endif
    static size_t isize = sizeof(INT);

    BTYPE *type0[4], *type1[4];

    /* pointers to old counts and maps */
    INT n0[4], *map0[4];
    /* pointers to new counts, maps and sorted local indices */
    INT n1[4], *map1[4], *index1[4];
    /* pointers to temporary maps */
    INT *map2[4];

    DOF **dof, *u;
    INT i, j;
    int ii, jj, kk, cnt, cnt1, cnt_vert, cnt_edge, cnt_face, cnt_elem;
    /* data count per vertex/edge/face/element for non h-p DOFs */
    int data_count[4] = {0, 0, 0, 0};

    int ndof_total;	/* for debugging */
    int ndof;		/* number of non h-p DOFs which have data */
    int ndof_hp;	/* number of h-p DOFs */
    int nhp;		/* number of HP_TYPEs */
    /* number of HP_TYPEs having vert/edge/face/elem data */
    int nhp_vert, nhp_edge, nhp_face, nhp_elem;
    DOF **hp_dofs;	/* list of h-p DOFs */
    HP_TYPE **tmphps;	/* list of temporary HP_TYPEs */
    HP_TYPE *hp = NULL;	/* scratch variable */
    /* data counts to/from each process */
    INT (*sdata_counts)[10], (*rdata_counts)[10];
    int *scounts, *sdispls, *rcounts, *rdispls;
    BYTE *sbuff, *rbuff, *buff;
    size_t size, size0;

    MPI_Datatype type;

    assert((blksize % sizeof(INT)) == 0 && (blksize % sizeof(FLOAT)) == 0);

    phgDofClearCache(g, NULL, NULL, NULL, FALSE);

    /* count number of non h-p and h-p DOFs */
    ndof_total = ndof = ndof_hp = 0;
    for (dof = g->dof; dof != NULL && (u = *dof) != NULL; dof++) {
	ndof_total++;
#if IGNORE_NOACTION
	if (!DofIsHP(u) && u->userfunc == DofNoAction)
	    continue;
#endif	/* IGNORE_NOACTION */
	if (SpecialDofType(u->type))
	    continue;
	!DofIsHP(u) ? ndof++ : ndof_hp++;
    }

    /* count number of HP_TYPEs */
    nhp = nhp_vert = nhp_edge = nhp_face = nhp_elem = 0;
    for (tmphps = g->hp; tmphps != NULL && *tmphps != NULL; tmphps++) {
	(*tmphps)->id = nhp++;
	if ((*tmphps)->info->flags & VERT_FLAG)
	    nhp_vert++;
	if ((*tmphps)->info->flags & EDGE_FLAG)
	    nhp_edge++;
	if ((*tmphps)->info->flags & FACE_FLAG)
	    nhp_face++;
	if ((*tmphps)->info->flags & ELEM_FLAG)
	    nhp_elem++;
    }

    if (ndof == 0 && nhp == 0 && ndof_hp == 0)
	return;

    hp_dofs = phgAlloc(ndof_hp * sizeof(*hp_dofs));

    /* count total data size on vertex, edge, face and element respectively for
     * non h-p DOFs, and copy h-p DOFs to hp_dofs[] */
    for (dof = g->dof, ii = 0; dof != NULL && (u = *dof) != NULL; dof++) {
#if IGNORE_NOACTION
	if (!DofIsHP(u) && u->userfunc == DofNoAction)
	    continue;
#endif	/* IGNORE_NOACTION */
	if (SpecialDofType(u->type))
	    continue;
	if (!DofIsHP(u)) {
	    /* fixed order DOF */
	    data_count[0] += u->count_vert;
	    data_count[1] += u->count_edge;
	    data_count[2] += u->count_face;
	    data_count[3] += u->count_elem;
	    continue;
	}

	hp_dofs[ii++] = u;
    }
    assert(ii == ndof_hp);

    /* initialize arrays */
    n0[0] = nvert;
    n0[1] = nedge;
    n0[2] = nface;
    n0[3] = nelem;
    n1[0] = g->nvert;
    n1[1] = g->nedge;
    n1[2] = g->nface;
    n1[3] = g->nelem;
    map0[0] = L2Gmap_vert;
    map0[1] = L2Gmap_edge;
    map0[2] = L2Gmap_face;
    map0[3] = L2Gmap_elem;
    map1[0] = g->L2Gmap_vert;
    map1[1] = g->L2Gmap_edge;
    map1[2] = g->L2Gmap_face;
    map1[3] = g->L2Gmap_elem;
    type0[0] = types_vert;
    type0[1] = types_edge;
    type0[2] = types_face;
    type0[3] = types_elem;
    type1[0] = g->types_vert;
    type1[1] = g->types_edge;
    type1[2] = g->types_face;
    type1[3] = g->types_elem;

    if (flag == FALSE) {
	BYTE *flags, *flags0;	/* used to mark vertices/edges/faces/elements */
	BYTE v0, v1, v2;
	ELEMENT *e;
	INT k;

	/* sort old maps */
	index0[0] = sort_map_indices(g->nvert, g->L2Gmap_vert);
	index0[1] = sort_map_indices((g->flags & EDGE_FLAG) ? g->nedge : 0,
					g->L2Gmap_edge);
	index0[2] = sort_map_indices((g->flags & FACE_FLAG) ? g->nface : 0,
					g->L2Gmap_face);
	index0[3] = sort_map_indices((g->flags & ELEM_FLAG) ? g->nelem : 0,
					g->L2Gmap_elem);

	/* construct lists of vertices/edges/faces/elements to send */
	lists = phgCalloc(g->nprocs, sizeof(*lists));
	size = g->nface;
	if (size < g->nvert)
	    size = g->nvert;
	if (size < g->nedge)
	    size = g->nedge;
	if (size < g->nelem)
	    size = g->nelem;
	
	flags = phgAlloc((size + g->nprocs) * sizeof(*flags));
	/* determine to which processes have data to send */
	flags0 = flags + size;
	memset(flags0, 0, g->nprocs * sizeof(*flags0));
	ForAllElements(g, e) {
	    if (e->mark == g->rank)
		continue;
	    flags0[e->mark] = 1;
	}
	/* construct individual lists */
	for (ii = 0; ii < g->nprocs; ii++) {
	    RNEIGHBOUR *rn, *rn1;
	    if (flags0[ii] == 0) {
		lists[ii] = phgCalloc(8, sizeof(*lists[ii]));
		continue;
	    }
	    /* construct list for process ii:
	     *		bit 0/1/2/3 = vert/edge/face/elem */
	    memset(flags, 0, size * sizeof(*flags));
	    ForAllElements(g, e) {
		if (e->mark != ii)
		    continue;
		if (ndof_hp > 0 || data_count[0] > 0)
		    for (jj = 0; jj < NVert; jj++)
			flags[e->verts[jj]] |= VERT_FLAG;
		if (ndof_hp > 0 || data_count[1] > 0)
		    for (jj = 0; jj < NEdge; jj++)
			flags[e->edges[jj]] |= EDGE_FLAG;
		if (ndof_hp > 0 || data_count[2] > 0)
		    for (jj = 0; jj < NFace; jj++)
			flags[e->faces[jj]] |= FACE_FLAG;
		if (ndof_hp > 0 || data_count[3] > 0)
		    flags[e->index] |= ELEM_FLAG;
	    }
	    /* remove vertices/edges/faces on shared faces (remote process
	     * already has data) */
	    rn = g->neighbours.list + g->neighbours.displs[ii];
	    rn1 = rn + g->neighbours.counts[ii];
	    for (; rn < rn1; rn++) {
		if (rn->rank != ii)
		    continue;
		jj = rn->vertex;
		e = rn->local;
		v0 = GetFaceVertex(jj, 0);
		v1 = GetFaceVertex(jj, 1);
		v2 = GetFaceVertex(jj, 2);
		/* remove vertices on the face from the list */
		flags[e->verts[v0]] &= ~VERT_FLAG;
		flags[e->verts[v1]] &= ~VERT_FLAG;
		flags[e->verts[v2]] &= ~VERT_FLAG;
		/* remove edges on the face from the list */
		flags[e->edges[GetEdgeNo(v0, v1)]] &= ~EDGE_FLAG;
		flags[e->edges[GetEdgeNo(v0, v2)]] &= ~EDGE_FLAG;
		flags[e->edges[GetEdgeNo(v1, v2)]] &= ~EDGE_FLAG;
		/* remove the face */
		flags[e->faces[jj]] &= ~FACE_FLAG;
	    }
	    /* count # of vertices/edges/faces/elements */
	    nvert = 0;
	    if (ndof_hp > 0 || data_count[0] > 0) {
		for (i = 0; i < g->nvert; i++) {
		    if (flags[i] & VERT_FLAG)
			nvert++;
		}
	    }
	    nedge = 0;
	    if (ndof_hp > 0 || data_count[1] > 0) {
		for (i = 0; i < g->nedge; i++) {
		    if (flags[i] & EDGE_FLAG)
			nedge++;
		}
	    }
	    nface = 0;
	    if (ndof_hp > 0 || data_count[2] > 0) {
		for (i = 0; i < g->nface; i++) {
		    if (flags[i] & FACE_FLAG)
			nface++;
		}
	    }
	    nelem = 0;
	    if (ndof_hp > 0 || data_count[3] > 0) {
		for (i = 0; i < g->nelem; i++) {
		    if (flags[i] & ELEM_FLAG)
			nelem++;
		}
	    }
	    /* allocate list */
	    size0 = 8 + nvert + nedge + nface + nelem +
		   (nhp_vert * nvert + isize - 1) / isize +
		   (nhp_edge * nedge + isize - 1) / isize +
		   (nhp_face * nface + isize - 1) / isize +
		   (nhp_elem * nelem + isize - 1) / isize;
	    list = lists[ii] = phgAlloc(size0 * sizeof(*list));
	    list[0] = nvert;	/* number of vertices to send */
	    list[1] = nedge;	/* number of edges to send */
	    list[2] = nface;	/* number of faces to send */
	    list[3] = nelem;	/* number of elements to send */

	    /* construct list of vertices/edges/faces/elements,
	     * and count data size for h-p DOFs */
	    list[4] = list[5] = list[6] = list[7] = 0;	/* h-p data sizes */
	    j = 8;
#define SaveOrder(nfoo, np_foo, mask, GlobalFoo, nhp_foo, \
		  index0_foo, list_foo, foo_order) \
	    if (nfoo > 0) { \
		buff = (BYTE *)(list + j + nfoo); \
		for (i = 0; i < g->nfoo; i++) { \
		    k = index0_foo[i]; \
		    if (flags[k] & mask) { \
			list[j++] = GlobalFoo(g, k); \
			/* save orders */ \
			for (jj = kk = 0; jj < nhp; jj++) \
			    if ((g->hp[jj]->info->flags & mask)) \
				buff[kk++ * nfoo] = g->hp[jj]->foo_order[k]; \
			assert(kk == nhp_foo); \
			/* count h-p data */ \
			for (jj = 0; jj < ndof_hp; jj++) { \
			    u = hp_dofs[jj]; \
			    if (!(u->hp->info->flags & mask)) \
				continue; \
			    cnt = u->hp->foo_order[k]; \
			    cnt = u->hp->info->types[cnt]->np_foo; \
			    list_foo += u->dim * cnt; \
			} \
			buff++; \
		    } \
		} \
		j += (nhp_foo * nfoo + sizeof(*list) - 1) / sizeof(*list); \
	    }
	    SaveOrder(nvert, np_vert, VERT_FLAG, GlobalVertex,
		      nhp_vert, index0[0], list[4], vert_order)
	    SaveOrder(nedge, np_edge, EDGE_FLAG, GlobalEdge,
		      nhp_edge, index0[1], list[5], edge_order)
	    SaveOrder(nface, np_face, FACE_FLAG, GlobalFace,
		      nhp_face, index0[2], list[6], face_order)
	    SaveOrder(nelem, np_elem, ELEM_FLAG, GlobalElement,
		      nhp_elem, index0[3], list[7], elem_order)
#undef SaveOrder
	    assert(j == size0);
	}
	phgFree(flags);

	/* save maps and types arrays */
	if ((g->flags & ELEM_FLAG)) {
	    /* save L2Gmap_elem and types_elem */
	    size = (nelem = g->nelem) * sizeof(*L2Gmap_elem);
	    phgFree(L2Gmap_elem);
	    L2Gmap_elem = phgAlloc(size);
	    if (g->L2Gmap_elem != NULL) {
		memcpy(L2Gmap_elem, g->L2Gmap_elem, size);
	    }
	    else {
		for (i = 0; i < nelem; i++)
		    L2Gmap_elem[i] = i;
	    }
	    types_elem = phgAlloc(nelem * sizeof(*types_elem));
	    memcpy(types_elem, g->types_elem, nelem * sizeof(*types_elem));
	}

	/* save L2Gmap_vert and types_vert */
	size = (nvert = g->nvert) * sizeof(*L2Gmap_vert);
	phgFree(L2Gmap_vert);
	L2Gmap_vert = phgAlloc(size);
	if (g->L2Gmap_vert != NULL) {
	    memcpy(L2Gmap_vert, g->L2Gmap_vert, size);
	}
	else {
	    for (i = 0; i < nvert; i++)
		L2Gmap_vert[i] = i;
	}
	types_vert = phgAlloc(nvert * sizeof(*types_vert));
	memcpy(types_vert, g->types_vert, nvert * sizeof(*types_vert));

	if ((g->flags & EDGE_FLAG)) {
	    /* save L2Gmap_edge and types_edge */
	    size = (nedge = g->nedge) * sizeof(*L2Gmap_edge);
	    phgFree(L2Gmap_edge);
	    L2Gmap_edge = phgAlloc(size);
	    if (g->L2Gmap_edge != NULL) {
		memcpy(L2Gmap_edge, g->L2Gmap_edge, size);
	    }
	    else {
		for (i = 0; i < nedge; i++)
		    L2Gmap_edge[i] = i;
	    }
	    types_edge = phgAlloc(nedge * sizeof(*types_edge));
	    memcpy(types_edge, g->types_edge, nedge * sizeof(*types_edge));
	}

	if ((g->flags & FACE_FLAG)) {
	    /* save L2Gmap_face and types_face */
	    size = (nface = g->nface) * sizeof(*L2Gmap_face);
	    phgFree(L2Gmap_face);
	    L2Gmap_face = phgAlloc(size);
	    if (g->L2Gmap_face != NULL) {
		memcpy(L2Gmap_face, g->L2Gmap_face, size);
	    }
	    else {
		for (i = 0; i < nface; i++)
		    L2Gmap_face[i] = i;
	    }
	    types_face = phgAlloc(nface * sizeof(*types_face));
	    memcpy(types_face, g->types_face, nface * sizeof(*types_face));
	}

	phgFree(hp_dofs);
	return;
    }

    phgDebug((3, "old: nvert=%d, nedge=%d, nface=%d, nelem=%d\n",
	      n0[0], n0[1], n0[2], n0[3]));
    phgDebug((3, "new: nvert=%d, nedge=%d, nface=%d, nelem=%d\n",
	      n1[0], n1[1], n1[2], n1[3]));

    if (!(g->flags & EDGE_FLAG))
	n0[1] = n1[1] = 0;
    if (!(g->flags & FACE_FLAG))
	n0[2] = n1[2] = 0;
    if (!(g->flags & ELEM_FLAG))
	n0[3] = n1[3] = 0;

    /* sort new maps */
    for (ii = 0; ii < 4; ii++) {
	index1[ii] = sort_map_indices(n1[ii], map1[ii]);
    }

    /* exchange data counts */
    sdata_counts = phgAlloc(2 * g->nprocs * sizeof(*sdata_counts)); 
    rdata_counts = sdata_counts + g->nprocs;
    for (ii = 0; ii < g->nprocs; ii++) {
	memcpy(sdata_counts[ii], lists[ii], 8 * sizeof(*(lists[ii])));
	sdata_counts[ii][8] = ndof_total;
	sdata_counts[ii][9] = ndof_hp;
    }
    MPI_Alltoall(sdata_counts, 10, PHG_MPI_INT,
		 rdata_counts, 10, PHG_MPI_INT, g->comm);

    /* check that total number of DOFs match */
    cnt = cnt1 = rdata_counts[0][8];
    for (ii = 1; ii < g->nprocs; ii++) {
	if (cnt > rdata_counts[ii][8])
	    cnt = rdata_counts[ii][8];	/* cnt = min ndof_total */
	else if (cnt1 < rdata_counts[ii][8])
	    cnt1 = rdata_counts[ii][8];	/* cnt1 = max ndof_total */
    }
    if (cnt != cnt1) {
	if (g->rank == 0) {
	    assert(cnt1 == ndof_total);
	    phgInfo(-1, "Fatal error: DOFs mismatch between submeshes.\n");
	    phgInfo(-1, "The following DOFs were (possibly) not "
			"properly freed:\n");
	    for (ii = cnt; ii < cnt1; ii++)
		phgInfo(-1, "\t\"%s\" (created at %s:%d)\n", g->dof[ii]->name,
				g->dof[ii]->srcfile, g->dof[ii]->srcline);
	    phgInfo(-1, "Note: a DOF created after phgRedistributeGrid must be "
			"freed before the next call to phgRedistributeGrid.\n");
	    phgInfo(-1, "Abort.\n");
	}
	MPI_Barrier(g->comm);
	MPI_Abort(phgComm, 1);
	exit(1);
    }

    /* check that total number of h-p DOFs match */
    cnt = cnt1 = rdata_counts[0][9];
    for (ii = 1; ii < g->nprocs; ii++) {
	if (cnt > rdata_counts[ii][9])
	    cnt = rdata_counts[ii][9];	/* cnt = min ndof_hp */
	else if (cnt1 < rdata_counts[ii][9])
	    cnt1 = rdata_counts[ii][9];	/* cnt1 = max ndof_hp */
    }
    if (cnt != cnt1) {
	if (g->rank == 0) {
	    assert(cnt1 == ndof_hp);
	    phgInfo(-1, "Fatal error: h-p DOFs mismatch between submeshes.\n");
	    phgInfo(-1, "Abort.\n");
	}
	MPI_Barrier(g->comm);
	MPI_Abort(phgComm, 1);
	exit(1);
    }

    assert(sizeof(FLOAT) % sizeof(INT) == 0 ||
	   sizeof(INT) % sizeof(FLOAT) == 0);

    /* prepare lists and data to send */
    scounts = phgAlloc(4 * g->nprocs * sizeof(*scounts));
    sdispls = scounts + g->nprocs * 1;
    rcounts = scounts + g->nprocs * 2;
    rdispls = scounts + g->nprocs * 3;
    for (ii = 0; ii < g->nprocs; ii++) {
	/* send size */
	size =	/* size for the list */
	       (sdata_counts[ii][0] + sdata_counts[ii][1] +
		sdata_counts[ii][2] + sdata_counts[ii][3]) * sizeof(INT) + 
		/* size for non h-p data */
	       (sdata_counts[ii][0] * data_count[0] +
		sdata_counts[ii][1] * data_count[1] +
		sdata_counts[ii][2] * data_count[2] +
		sdata_counts[ii][3] * data_count[3]) * sizeof(FLOAT) +
		/* size for h-p orders */
	       ((nhp_vert * sdata_counts[ii][0] + isize - 1) / isize +
	        (nhp_edge * sdata_counts[ii][1] + isize - 1) / isize +
		(nhp_face * sdata_counts[ii][2] + isize - 1) / isize +
		(nhp_elem * sdata_counts[ii][3] + isize - 1) / isize) * isize +
		/* size for h-p data */
	       (sdata_counts[ii][4] + sdata_counts[ii][5] +
		sdata_counts[ii][6] + sdata_counts[ii][7]) * sizeof(FLOAT);
	if ((jj = size % blksize) != 0)
	    size += blksize - jj;
	size /= blksize;
	assert(size <= INT_MAX);
	scounts[ii] = size;
	sdispls[ii] = (ii == 0 ? 0 : sdispls[ii - 1] + scounts[ii - 1]);
	assert(sdispls[ii] + (INT)scounts[ii] <= INT_MAX);

	/* receive size */
	size =	/* size for the list */
	       (rdata_counts[ii][0] + rdata_counts[ii][1] +
		rdata_counts[ii][2] + rdata_counts[ii][3]) * sizeof(INT) + 
		/* size for non h-p data */
	       (rdata_counts[ii][0] * data_count[0] +
		rdata_counts[ii][1] * data_count[1] +
		rdata_counts[ii][2] * data_count[2] +
		rdata_counts[ii][3] * data_count[3]) * sizeof(FLOAT) +
		/* size for h-p orders */
	       ((nhp_vert * rdata_counts[ii][0] + isize - 1) / isize +
	        (nhp_edge * rdata_counts[ii][1] + isize - 1) / isize +
		(nhp_face * rdata_counts[ii][2] + isize - 1) / isize +
		(nhp_elem * rdata_counts[ii][3] + isize - 1) / isize) * isize +
		/* size for h-p data */
	       (rdata_counts[ii][4] + rdata_counts[ii][5] +
		rdata_counts[ii][6] + rdata_counts[ii][7]) * sizeof(FLOAT);
	if ((jj = size % blksize) != 0)
	    size += blksize - jj;
	size /= blksize;
	assert(size <= INT_MAX);
	rcounts[ii] = size;
	rdispls[ii] = (ii == 0 ? 0 : rdispls[ii - 1] + rcounts[ii - 1]);
	assert(rdispls[ii] + (INT)rcounts[ii] <= INT_MAX);
    }

    tmphps = phgAlloc(nhp * sizeof(*tmphps));

    /* copy lists and data to send */
    size = sdispls[g->nprocs - 1] + scounts[g->nprocs - 1];
    sbuff = phgAlloc(size * blksize);
    for (cnt = 0; cnt < nhp; cnt++) {
	hp = tmphps[cnt] = phgCalloc(1, sizeof(*hp));
	hp->info = g->hp[cnt]->info;
    }
    for (ii = 0; ii < g->nprocs; ii++) {
	if (scounts[ii] == 0)
	    continue;
	/* copy list */
	buff = sbuff + sdispls[ii] * blksize;
	size = (sdata_counts[ii][0] + sdata_counts[ii][1] +
                sdata_counts[ii][2] + sdata_counts[ii][3] +
	        (nhp_vert * sdata_counts[ii][0] + isize - 1) / isize +
	        (nhp_edge * sdata_counts[ii][1] + isize - 1) / isize +
		(nhp_face * sdata_counts[ii][2] + isize - 1) / isize +
		(nhp_elem * sdata_counts[ii][3] + isize - 1) / isize) * isize;
	memcpy(buff, lists[ii] + 8, size);
	if ((jj = size % sizeof(FLOAT)) != 0)
	    size += sizeof(FLOAT) - jj;
	/* copy DOF data for process ii */
	map2[0] = (INT *)buff;
	map2[1] = map2[0] + sdata_counts[ii][0] +
		  (nhp_vert * sdata_counts[ii][0] + isize - 1) / isize;
	map2[2] = map2[1] + sdata_counts[ii][1] +
		  (nhp_edge * sdata_counts[ii][1] + isize - 1) / isize;
	map2[3] = map2[2] + sdata_counts[ii][2] +
		  (nhp_face * sdata_counts[ii][2] + isize - 1) / isize;
	/* setup temporary HP_TYPEs */
	cnt_vert = cnt_edge = cnt_face = cnt_elem = 0;
	for (cnt = 0; cnt < nhp; cnt++) {
	    hp = tmphps[cnt];

	    if (g->hp[cnt]->info->flags & VERT_FLAG)
		hp->vert_order = ((BYTE *)(map2[0] + sdata_counts[ii][0])) +
				 cnt_vert++ * sdata_counts[ii][0];
	    else
		hp->vert_order = NULL;

	    if (g->hp[cnt]->info->flags & EDGE_FLAG)
		hp->edge_order = ((BYTE *)(map2[1] + sdata_counts[ii][1])) +
				 cnt_edge++ * sdata_counts[ii][1];
	    else
		hp->edge_order = NULL;

	    if (g->hp[cnt]->info->flags & FACE_FLAG)
		hp->face_order = ((BYTE *)(map2[2] + sdata_counts[ii][2])) +
				 cnt_face++ * sdata_counts[ii][2];
	    else
		hp->face_order = NULL;

	    if (g->hp[cnt]->info->flags & ELEM_FLAG)
		hp->elem_order = ((BYTE *)(map2[3] + sdata_counts[ii][3])) +
				 cnt_elem++ * sdata_counts[ii][3];
	    else
		hp->elem_order = NULL;

	    phgHPSetupDataPointers_(hp, sdata_counts[ii][0],
					sdata_counts[ii][1],
					sdata_counts[ii][2],
					sdata_counts[ii][3], FALSE, FALSE);
	}
	assert(cnt_vert == nhp_vert && cnt_edge == nhp_edge &&
	       cnt_face == nhp_face && cnt_elem == nhp_elem);
	/* copy DOF data */
	for (dof = g->dof; dof != NULL && (u = *dof) != NULL; dof++) {
#if IGNORE_NOACTION
	    if (u->userfunc == DofNoAction)
		continue;
#endif	/* IGNORE_NOACTION */
	    if (DofIsHP(u))
		hp = tmphps[u->hp->id];
	    else
		hp = NULL;
	    copy_dof_data(u, u->data, n0, index0, map0, u->hp, type0,
		(FLOAT *)(buff + size), sdata_counts[ii], NULL, map2, hp, NULL);
	    if (hp == NULL) {
		size += (sdata_counts[ii][0] * u->count_vert +
			 sdata_counts[ii][1] * u->count_edge +
			 sdata_counts[ii][2] * u->count_face +
			 sdata_counts[ii][3] * u->count_elem) * sizeof(FLOAT);
	    }
	    else {
		if (u->hp->info->flags & VERT_FLAG)
		    size += sizeof(FLOAT) * u->dim *
					hp->vert_index[sdata_counts[ii][0]];
		if (u->hp->info->flags & EDGE_FLAG)
		    size += sizeof(FLOAT) * u->dim *
					hp->edge_index[sdata_counts[ii][1]];
		if (u->hp->info->flags & FACE_FLAG)
		    size += sizeof(FLOAT) * u->dim *
					hp->face_index[sdata_counts[ii][2]];
		if (u->hp->info->flags & ELEM_FLAG)
		    size += sizeof(FLOAT) * u->dim *
					hp->elem_index[sdata_counts[ii][3]];
	    }
	}
	assert((size + blksize - 1) / blksize == scounts[ii]);
    }

    /* free arrays in tmphps[] */
    for (cnt = 0; cnt < nhp; cnt++) {
	hp = tmphps[cnt];
	hp->vert_order = NULL;
	hp->edge_order = NULL;
	hp->face_order = NULL;
	hp->elem_order = NULL;
	phgHPFreeBuffers__(hp);
    }

    /* free lists */
    for (ii = 0; ii < g->nprocs; ii++)
	phgFree(lists[ii]);
    phgFree(lists);
    lists = NULL;

    /* send lists and data */
    size = rdispls[g->nprocs - 1] + rcounts[g->nprocs - 1];
    rbuff = phgAlloc(size * blksize);
    MPI_Type_contiguous(blksize, MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(sbuff, scounts, sdispls, type, rbuff, rcounts, rdispls, type,
		  g->comm);
    MPI_Type_free(&type);

    /* free send buffer */
    phgFree(sbuff);

    /* FIXME: use flags[] to only copy data once, and check no data missing? */

    /* backup old HP_TYPEs */
    cnt_vert = cnt_edge = cnt_face = cnt_elem = 0;
    for (cnt = 0; cnt < nhp; cnt++) {
	BYTE *orders0[4], *orders1[4];

	hp = g->hp[cnt];
	phgFree(tmphps[cnt]);
	tmphps[cnt] = phgHPDup__(hp, TRUE);	/* save old HP_TYPE */

	/* set up new orders and pointers */
	if (hp->info->flags & VERT_FLAG) {
	    hp->vert_order = phgCalloc(g->nvert, sizeof(*hp->vert_order));
	    hp->vert_index = phgAlloc((g->nvert + 1) * sizeof(*hp->vert_index));
	}

	if (hp->info->flags & EDGE_FLAG) {
	    hp->edge_order = phgCalloc(g->nedge, sizeof(*hp->edge_order));
	    hp->edge_index = phgAlloc((g->nedge + 1) * sizeof(*hp->edge_index));
	}

	if (hp->info->flags & FACE_FLAG) {
	    hp->face_order = phgCalloc(g->nface, sizeof(*hp->face_order));
	    hp->face_index = phgAlloc((g->nface + 1) * sizeof(*hp->face_index));
	}

	if (hp->info->flags & ELEM_FLAG) {
	    hp->elem_order = phgCalloc(g->nelem, sizeof(*hp->elem_order));
	    hp->elem_index = phgAlloc((g->nelem + 1) * sizeof(*hp->elem_index));
	}

	orders1[0] = hp->vert_order;
	orders1[1] = hp->edge_order;
	orders1[2] = hp->face_order;
	orders1[3] = hp->elem_order;

	/* pick-up orders from remote entries */
	for (ii = 0; ii < g->nprocs; ii++) {
	    if (rdata_counts[ii][0] + rdata_counts[ii][1] +
		rdata_counts[ii][2] + rdata_counts[ii][3] == 0)
		continue;
	    map2[0] = (INT *)(rbuff + rdispls[ii] * blksize);
	    map2[1] = map2[0] + rdata_counts[ii][0] +
			(nhp_vert * rdata_counts[ii][0] + isize - 1) / isize;
	    map2[2] = map2[1] + rdata_counts[ii][1] +
			(nhp_edge * rdata_counts[ii][1] + isize - 1) / isize;
	    map2[3] = map2[2] + rdata_counts[ii][2] +
			(nhp_face * rdata_counts[ii][2] + isize - 1) / isize;

	    orders0[0] = (BYTE *)(map2[0] + rdata_counts[ii][0])
					  + cnt_vert * rdata_counts[ii][0];
	    orders0[1] = (BYTE *)(map2[1] + rdata_counts[ii][1])
					  + cnt_edge * rdata_counts[ii][1];
	    orders0[2] = (BYTE *)(map2[2] + rdata_counts[ii][2])
					  + cnt_face * rdata_counts[ii][2];
	    orders0[3] = (BYTE *)(map2[3] + rdata_counts[ii][3])
					  + cnt_elem * rdata_counts[ii][3];

	    copy_hp_orders(orders0, rdata_counts[ii], NULL, map2, NULL,
			   orders1, n1, index1, map1);
	}

	/* pick-up orders from local entries (u->hp[]) */
	orders0[0] = tmphps[cnt]->vert_order;
	orders0[1] = tmphps[cnt]->edge_order;
	orders0[2] = tmphps[cnt]->face_order;
	orders0[3] = tmphps[cnt]->elem_order;

	copy_hp_orders(orders0, n0, index0, map0, type0,
		       orders1, n1, index1, map1);

	/* set up data pointers */
	phgHPSetupDataPointers_(hp, g->nvert, g->nedge, g->nface, g->nelem,
				TRUE, FALSE);

	if (hp->info->flags & VERT_FLAG)
	    cnt_vert++;
	if (hp->info->flags & EDGE_FLAG)
	    cnt_edge++;
	if (hp->info->flags & FACE_FLAG)
	    cnt_face++;
	if (hp->info->flags & ELEM_FLAG)
	    cnt_elem++;
    }
    assert(cnt_vert == nhp_vert && cnt_edge == nhp_edge &&
	   cnt_face == nhp_face && cnt_elem == nhp_elem);

    /* copy local DOF data */
    for (dof = g->dof; dof != NULL && (u = *dof) != NULL; dof++) {
	FLOAT *old_data, *new_data;

#if IGNORE_NOACTION
	if (u->userfunc == DofNoAction) {
	    if (!SpecialDofType(u->type)) {
		phgFree(u->data);
		u->data = phgCalloc(DofDataCount(u), sizeof(*old_data));
		u->data_vert = u->data;
		u->data_edge = u->data_vert + DofVertexDataCount(u);
		u->data_face = u->data_edge + DofEdgeDataCount(u);
		u->data_elem = u->data_face + DofFaceDataCount(u);
	    }
	    continue;
	}
#endif	/* IGNORE_NOACTION */
	if (SpecialDofType((*dof)->type))
	    continue;

	/* hp = old HP_TYPE */
	if (DofIsHP(u))
	    hp = tmphps[u->hp->id];
	else
	    hp = NULL;

	old_data = u->data;
	u->data = phgCalloc(DofDataCount(u), sizeof(*old_data));
	u->data_vert = new_data = u->data;
	u->data_edge = u->data_vert + DofVertexDataCount(u);
	u->data_face = u->data_edge + DofEdgeDataCount(u);
	u->data_elem = u->data_face + DofFaceDataCount(u);
	copy_dof_data(u, old_data, n0, index0, map0, hp, type0,
			 new_data, n1, index1, map1, u->hp, type1);
	phgFree(old_data);
    }

    /* copy received data */
    for (ii = 0; ii < g->nprocs; ii++) {
	if (rcounts[ii] == 0)
	    continue;
	/* copy list */
	buff = rbuff + rdispls[ii] * blksize;
	size = (rdata_counts[ii][0] + rdata_counts[ii][1] +
                rdata_counts[ii][2] + rdata_counts[ii][3] +
	        (nhp_vert * rdata_counts[ii][0] + isize - 1) / isize +
	        (nhp_edge * rdata_counts[ii][1] + isize - 1) / isize +
		(nhp_face * rdata_counts[ii][2] + isize - 1) / isize +
		(nhp_elem * rdata_counts[ii][3] + isize - 1) / isize) * isize;
	if ((jj = size % sizeof(FLOAT)) != 0)
	    size += sizeof(FLOAT) - jj;

	map2[0] = (INT *)buff;
	map2[1] = map2[0] + rdata_counts[ii][0] +
		  (nhp_vert * rdata_counts[ii][0] + isize - 1) / isize;
	map2[2] = map2[1] + rdata_counts[ii][1] +
		  (nhp_edge * rdata_counts[ii][1] + isize - 1) / isize;
	map2[3] = map2[2] + rdata_counts[ii][2] +
		  (nhp_face * rdata_counts[ii][2] + isize - 1) / isize;

	/* set up tmphps[] for received data */
	cnt_vert = cnt_edge = cnt_face = cnt_elem = 0;
	for (cnt = 0; cnt < nhp; cnt++) {
	    hp = tmphps[cnt];
	    phgHPFreeBuffers__(hp);
	    if (g->hp[cnt]->info->flags & VERT_FLAG)
		hp->vert_order = (BYTE *)(map2[0] + rdata_counts[ii][0])
					  + cnt_vert++ * rdata_counts[ii][0];
	    if (g->hp[cnt]->info->flags & EDGE_FLAG)
		hp->edge_order = (BYTE *)(map2[1] + rdata_counts[ii][1])
					  + cnt_edge++ * rdata_counts[ii][1];
	    if (g->hp[cnt]->info->flags & FACE_FLAG)
		hp->face_order = (BYTE *)(map2[2] + rdata_counts[ii][2])
					  + cnt_face++ * rdata_counts[ii][2];
	    if (g->hp[cnt]->info->flags & ELEM_FLAG)
		hp->elem_order = (BYTE *)(map2[3] + rdata_counts[ii][3])
					  + cnt_elem++ * rdata_counts[ii][3];
	    phgHPSetupDataPointers_(hp, rdata_counts[ii][0],
					rdata_counts[ii][1],
					rdata_counts[ii][2],
					rdata_counts[ii][3], FALSE, FALSE);
	}
	assert(cnt_vert == nhp_vert && cnt_edge == nhp_edge &&
	       cnt_face == nhp_face && cnt_elem == nhp_elem);

	for (dof = g->dof, cnt = 0; dof != NULL && (u = *dof) != NULL; dof++) {
#if IGNORE_NOACTION
	    if (u->userfunc == DofNoAction)
		continue;
#endif	/* IGNORE_NOACTION */

	    if (SpecialDofType(u->type))
		continue;

	    if (DofIsHP(u))
		hp = tmphps[u->hp->id];
	    else
		hp = NULL;

	    copy_dof_data(u, (FLOAT *)(buff + size), rdata_counts[ii], NULL,
		map2, hp, NULL, u->data, n1, index1, map1, u->hp, type1);
	    if (hp == NULL) {
		size += (rdata_counts[ii][0] * u->count_vert +
			 rdata_counts[ii][1] * u->count_edge +
			 rdata_counts[ii][2] * u->count_face +
			 rdata_counts[ii][3] * u->count_elem) * sizeof(FLOAT);
	    }
	    else {
		if (u->hp->info->flags & VERT_FLAG)
		    size += sizeof(FLOAT) * u->dim *
					hp->vert_index[rdata_counts[ii][0]];
		if (u->hp->info->flags & EDGE_FLAG)
		    size += sizeof(FLOAT) * u->dim *
					hp->edge_index[rdata_counts[ii][1]];
		if (u->hp->info->flags & FACE_FLAG)
		    size += sizeof(FLOAT) * u->dim *
					hp->face_index[rdata_counts[ii][2]];
		if (u->hp->info->flags & ELEM_FLAG)
		    size += sizeof(FLOAT) * u->dim *
					hp->elem_index[rdata_counts[ii][3]];
	    }
	}
	assert((size + blksize - 1) / blksize == rcounts[ii]);
	for (cnt = 0; cnt < nhp; cnt++) {
	    hp = tmphps[cnt];
	    hp->vert_order = NULL;
	    hp->edge_order = NULL;
	    hp->face_order = NULL;
	    hp->elem_order = NULL;
	}
    }

    /* free temporary HP_TYPEs */
    for (cnt = 0; cnt < nhp; cnt++) {
	hp = tmphps[cnt];
	phgHPFreeBuffers__(hp);
	phgFree(hp);
    }
    phgFree(tmphps);

    /* free memory blocks and reset variables */
    phgFree(rbuff);
    phgFree(scounts);
    phgFree(sdata_counts);

    for (ii = 0; ii < 4; ii++) {
	phgFree(map0[ii]);
	phgFree(index0[ii]);
	phgFree(type0[ii]);
	phgFree(index1[ii]);
    }

    L2Gmap_elem = L2Gmap_vert = L2Gmap_edge = L2Gmap_face = NULL;
    types_elem = types_vert = types_edge = types_face = NULL;
    nelem = nvert = nedge = nface = 0;

    phgFree(hp_dofs);
#else	/* USE_MPI */
    Unused(g);
    Unused(flag);
#endif	/* USE_MPI */

    return;
}

/*----------------- utilities ---------------*/

INT
phgDofMapE2D_(DOF *u, ELEMENT *e, int index, BOOLEAN global_ordering)
/* returns DOF index of the unknown with element index 'index' on 'e'.
 *
 * if 'global_ordering' == TRUE then the order of edge and face DOFs are
 * adjusted according to global orientation of the edge/face */
{
    GRID *g = u->g;
    DOF_TYPE *type = u->type;
    HP_TYPE *hp;
    int m, n;
    INT no, index0;

    if (type == NULL)
	goto hp_dof;

    if (index < (n = (m = u->count_vert) * NVert)) {
	return m * e->verts[index / m] + (index % m);
    }
    index -= n;
    index0 = m * g->nvert;

    if (index < (n = (m = u->count_edge) * NEdge)) {
	if (global_ordering && type->np_edge > 1) {
	    int k = index / m, M[type->np_edge];
	    phgDofMapEdgeData(type, e, k, M);
	    index %= m;
	    return index0 + m * e->edges[k] + M[index / u->dim] * u->dim +
			index % u->dim;
	}
	else {
	    return index0 + m * e->edges[index / m] + (index % m);
	}
    }
    index -= n;
    index0 += m * g->nedge;

    if (index < (n = (m = u->count_face) * NFace)) {
	if (global_ordering && type->np_face > 1) {
	    int k = index / m, M[type->np_face];
	    phgDofMapFaceData(type, e, k, M);
	    index %= m;
	    return index0 + m * e->faces[k] + M[index / u->dim] * u->dim +
			index % u->dim;
	}
	else {
	    return index0 + m * e->faces[index / m] + (index % m);
	}
    }
    index -= n;
    index0 += m * g->nface;

    if (global_ordering && type->np_elem > 1 && type->DofMap != NULL &&
		!ElementIsInOrder(e)) {
	int M[type->np_elem];
	phgDofMapElementData(type, e, M);
	return index0 + u->count_elem * e->index +
			M[index / u->dim] * u->dim + index % u->dim;
    }
    else {
	return index0 + u->count_elem * e->index + index;
    }

hp_dof:
    index0 = index % u->dim;
    index /= u->dim;
    hp = u->hp;

    if (hp->info->flags & VERT_FLAG) {
	for (m = 0; m < NVert; m++) {
	    n = hp->info->types[hp->vert_order[no = e->verts[m]]]->np_vert;
	    if (index < n)
		return index0 + (hp->vert_index[no] + index) * u->dim;
	    index -= n;
	}
	index0 += hp->vert_index[g->nvert] * u->dim;
    }

    if (hp->info->flags & EDGE_FLAG) {
	for (m = 0; m < NEdge; m++) {
	    n = hp->info->types[hp->edge_order[no = e->edges[m]]]->np_edge;
	    if (index < n)
		return index0 + (hp->edge_index[no] + index) * u->dim;
	    index -= n;
	}
	index0 += hp->edge_index[g->nedge] * u->dim;
    }

    if (hp->info->flags & FACE_FLAG) {
	for (m = 0; m < NFace; m++) {
	    n = hp->info->types[hp->face_order[no = e->faces[m]]]->np_face;
	    if (index < n)
		return index0 + (hp->face_index[no] + index) * u->dim;
	    index -= n;
	}
	index0 += hp->face_index[g->nface] * u->dim;
    }

    assert(hp->info->flags & ELEM_FLAG);
    n = hp->info->types[hp->elem_order[no = e->index]]->np_elem;
    assert(index < n);
    return index0 + (hp->elem_index[no] + index) * u->dim;
}

BTYPE
phgDofGetElementBasisInfo_(DOF *u, ELEMENT *e, int bas_no,
			   GTYPE *gtype, int *gindex, int *bindex,
			   FLOAT xyz[], FLOAT lambda[])
/* returns:
 *	boundary type flags as return value,
 *	geometric type in '*gtype',
 *	vertex, edge or face index in '*gindex',
 *	and basis index on the edge, face or element in 'bindex',
 * of the basis 'bas_no' of 'u' on element 'e' */
{
    int l, m, n;
    GRID *g = u->g;
    DOF_TYPE *type = u->type;
    HP_TYPE *hp;
    FLOAT lam0[Dim + 1], *lam = (lambda == NULL ? lam0 : lambda);


    if (type == NULL)
	type = u->hp->info->types[u->hp->max_order];

    if (type->base_type != NULL) {
	int dim = type->dim / type->base_type->dim;
	BTYPE t;
	assert(!DofIsHP(u));	/* TODO: HP DOF */
	u->type = type->base_type;
	/*u->dim *= dim;*/
	t = phgDofGetElementBasisInfo_(u, e, bas_no / dim,
			   gtype, gindex, bindex, xyz, lambda);
	/*u->dim /= dim;*/
	u->type = type;
	if (bindex != NULL)
	    *bindex = *bindex * dim + bas_no % dim;
	return t;
    }

    if (!type->is_nodal && (xyz != NULL || lambda != NULL))
	phgError(1, "%s (%s:%d), DOF type \"%s\" has no nodal points.\n",
		 __func__, __FILE__, __LINE__, type->name);

    if (DofIsHP(u))
	goto hp_dof;

    if (bas_no < (n = (m = type->np_vert) * NVert)) {
	l = bas_no % m;	/* index on the vertex */
	m = bas_no / m;	/* vertex # */
	if (bindex != NULL)
	    *bindex = l;
	if (gtype != NULL)
	    *gtype = VERTEX;
	if (gindex != NULL)
	    *gindex = m;
	if (xyz != NULL || lambda != NULL) {
	    memset(lam, 0, (Dim + 1) * sizeof(*lam));
	    lam[m] = type->points[l];
	    if (xyz != NULL)
		phgGeomLambda2XYZ(g, e, lam, xyz, xyz + 1, xyz + 2);
	}
	return g->types_vert[e->verts[m]];
    }
    bas_no -= n;

    if (bas_no < (n = (m = type->np_edge) * NEdge)) {
	l = bas_no % m;	/* index on the edge */
	m = bas_no / m;	/* edge # */
	if (bindex != NULL)
	    *bindex = l;
	if (gtype != NULL)
	    *gtype = EDGE;
	if (gindex != NULL)
	    *gindex = m;
	if (xyz != NULL || lambda != NULL) {
	    memset(lam, 0, (Dim + 1) * sizeof(*lam));
	    lam[GetEdgeVertex(m, 0)] = type->points[type->np_vert + l * 2 + 0];
	    lam[GetEdgeVertex(m, 1)] = type->points[type->np_vert + l * 2 + 1];
	    if (xyz != NULL)
		phgGeomLambda2XYZ(g, e, lam, xyz, xyz + 1, xyz + 2);
	}
	return g->types_edge[e->edges[m]];
    }
    bas_no -= n;

    if (bas_no < (n = (m = type->np_face) * NFace)) {
	l = bas_no % m;	/* index on the face */
	m = bas_no / m;	/* face # */
	if (bindex != NULL)
	    *bindex = l;
	if (gtype != NULL)
	    *gtype = FACE;
	if (gindex != NULL)
	    *gindex = m;
	if (xyz != NULL || lambda != NULL) {
	    lam[m] = 0.;
	    lam[GetFaceVertex(m, 0)] = type->points[type->np_vert +
						    type->np_edge * 2 +
						    l * 3 + 0];
	    lam[GetFaceVertex(m, 1)] = type->points[type->np_vert +
						    type->np_edge * 2 +
						    l * 3 + 1];
	    lam[GetFaceVertex(m, 2)] = type->points[type->np_vert +
						    type->np_edge * 2 +
						    l * 3 + 2];
	    if (xyz != NULL)
		phgGeomLambda2XYZ(g, e, lam, xyz, xyz + 1, xyz + 2);
	}
	return g->types_face[e->faces[m]];
    }
    bas_no -= n;

    assert(bas_no < type->np_elem);
    if (bindex != NULL)
	*bindex = bas_no;
    if (gtype != NULL)
	*gtype = VOLUME;
    if (gindex != NULL)
	*gindex = 0;
    if (xyz != NULL || lambda != NULL) {
	lam = type->points + type->np_vert + type->np_edge * 2 +
		type->np_face * 3 + bas_no * 4;
	if (lambda != NULL)
	    memcpy(lambda, lam, 4 * sizeof(*lam));
	if (xyz != NULL)
	    phgGeomLambda2XYZ(g, e, lam, xyz, xyz + 1, xyz + 2);
    }
    return g->types_elem[e->index];

hp_dof:
    hp = u->hp;

    assert(xyz == NULL && lambda == NULL);	/* TODO */

    for (m = 0; (hp->info->flags & VERT_FLAG) && m < NVert; m++) {
	n = hp->info->types[hp->vert_order[e->verts[m]]]->np_vert;
	if (bas_no < n) {
	    if (bindex != NULL)
		*bindex = bas_no;
	    if (gtype != NULL)
		*gtype = VERTEX;
	    if (gindex != NULL)
		*gindex = m;
	    return g->types_vert[e->verts[m]];
	}
	bas_no -= n;
    }

    for (m = 0; (hp->info->flags & EDGE_FLAG) && m < NEdge; m++) {
	n = hp->info->types[hp->edge_order[e->edges[m]]]->np_edge;
	if (bas_no < n) {
	    if (bindex != NULL)
		*bindex = bas_no;
	    if (gtype != NULL)
		*gtype = EDGE;
	    if (gindex != NULL)
		*gindex = m;
	    return g->types_edge[e->edges[m]];
	}
	bas_no -= n;
    }

    for (m = 0; (hp->info->flags & FACE_FLAG) && m < NFace; m++) {
	n = hp->info->types[hp->face_order[e->faces[m]]]->np_face;
	if (bas_no < n) {
	    if (bindex != NULL)
		*bindex = bas_no;
	    if (gtype != NULL)
		*gtype = FACE;
	    if (gindex != NULL)
		*gindex = m;
	    return g->types_face[e->faces[m]];
	}
	bas_no -= n;
    }

    assert(hp->info->flags & ELEM_FLAG);
    assert(bas_no < hp->info->types[hp->elem_order[e->index]]->np_elem);
    if (bindex != NULL)
	*bindex = bas_no;
    if (gtype != NULL)
	*gtype = VOLUME;
    if (gindex != NULL)
	*gindex = 0;
    return g->types_elem[e->index];
}

FLOAT *
phgDofGetElementCoordinates(DOF *u, ELEMENT *e, int index)
/* returns world coordinates of the unknown with local DOF index 'index' */
{
    static FLOAT c[Dim];
#if USE_OMP
#pragma omp threadprivate(c)
#endif	/* USE_OMP */

    phgDofGetElementBasisInfo_(u, e, index / u->dim, NULL, NULL, NULL, c, NULL);

    return c;
}

void 
phgDofGetElementDatas(DOF *dof, SIMPLEX *e, FLOAT *value) 
/* Get all Dof element datas, stored in values[nbas][dof->dim] in order of bas point,
 * value buffer is provided by user. */
{
    int i, j, k, n;
    FLOAT *v;
    int dim = dof->dim;

    i = 0;

    if ((n = dof->type->np_vert) > 0) {
	for (k = 0; k < NVert; k++) {
	    v = DofVertexData(dof, e->verts[k]);
	    for (j = 0; j < n; j++, i++)
		memcpy(value + i*dim, v + j*dim, dim * sizeof(FLOAT));
	}
    }

    if ((n = dof->type->np_edge) > 1) {
	static int *M = NULL, M_size = 0;
	if (M_size < n) {
	    phgFree(M);
	    M = phgAlloc((M_size = n) * sizeof(*M));
	    FreeAtExit(M);
	}
	for (k = 0; k < NEdge; k++) {
	    v = DofEdgeData(dof, e->edges[k]);
	    phgDofMapEdgeData(dof->type, e, k, M);
	    for (j = 0; j < n; j++, i++)
		memcpy(value + i*dim, v + M[j]*dim, dim * sizeof(FLOAT));
	}
    }
    else if (n > 0) {
	for (k = 0; k < NEdge; k++) {
	    v = DofEdgeData(dof, e->edges[k]);
	    for (j = 0; j < n; j++, i++)
		memcpy(value + i*dim, v + j*dim, dim * sizeof(FLOAT));
	}
    }

    if ((n = dof->type->np_face) > 1) {
	static int *M = NULL, M_size = 0;
	if (M_size < n) {
	    phgFree(M);
	    M = phgAlloc((M_size = n) * sizeof(*M));
	    FreeAtExit(M);
	}
	for (k = 0; k < NFace; k++) {
	    v = DofFaceData(dof, e->faces[k]);
	    phgDofMapFaceData(dof->type, e, k, M);
	    for (j = 0; j < n; j++, i++)
		memcpy(value + i*dim, v + M[j]*dim, dim * sizeof(FLOAT));
	}
    }
    else if (n > 0) {
	for (k = 0; k < NFace; k++) {
	    v = DofFaceData(dof, e->faces[k]);
	    for (j = 0; j < n; j++, i++)
		memcpy(value + i*dim, v + j*dim, dim * sizeof(FLOAT));
	}
    }

    if ((n = dof->type->np_elem) > 1 && dof->type->DofMap != NULL &&
	!ElementIsInOrder(e)) {
	static int *M = NULL, M_size = 0;
	if (M_size < n) {
	    phgFree(M);
	    M = phgAlloc((M_size = n) * sizeof(*M));
	    FreeAtExit(M);
	}
	v = DofElementData(dof, e->index);
	phgDofMapElementData(dof->type, e, M);
	for (j = 0; j < n; j++, i++)
	    memcpy(value + i*dim, v + M[j]*dim, dim * sizeof(FLOAT));
    }
    else if (n > 0) {
	v = DofElementData(dof, e->index);
	for (j = 0; j < n; j++, i++)
	    memcpy(value + i*dim, v + j*dim, dim * sizeof(FLOAT));
    }

    return;
}

void 
phgDofSetElementDatas(DOF *dof, SIMPLEX *e, FLOAT *value) 
/* Set all Dof element datas, stored in values[nbas][dof->dim] in order of bas point,
 * value is provided by user.
 * Note: this function only change local Dof data, Dof data on other procs is not
 * affected. User should be manage the consistence of global DOF data.
 * */

{
    int i, j, k, n;
    FLOAT *v;
    int dim = dof->dim;

    i = 0;
    
    if ((n = dof->type->np_vert) > 0) {
	for (k = 0; k < NVert; k++) {
	    v = DofVertexData(dof, e->verts[k]);
	    for (j = 0; j < n; j++, i++)
		memcpy(v + j*dim, value + i*dim, dim * sizeof(FLOAT));
	}
    }

    if ((n = dof->type->np_edge) > 1) {
	static int *M = NULL, M_size = 0;
	if (M_size < n) {
	    phgFree(M);
	    M = phgAlloc((M_size = n) * sizeof(*M));
	    FreeAtExit(M);
	}
	for (k = 0; k < NEdge; k++) {
	    v = DofEdgeData(dof, e->edges[k]);
	    phgDofMapEdgeData(dof->type, e, k, M);
	    for (j = 0; j < n; j++, i++)
		memcpy(v + M[j]*dim, value + i*dim, dim * sizeof(FLOAT));
	}
    }
    else if (n > 0) {
	for (k = 0; k < NEdge; k++) {
	    v = DofEdgeData(dof, e->edges[k]);
	    for (j = 0; j < n; j++, i++)
		memcpy(v + j*dim, value + i*dim, dim * sizeof(FLOAT));
	}
    }

    if ((n = dof->type->np_face) > 1) {
	static int *M = NULL, M_size = 0;
	if (M_size < n) {
	    phgFree(M);
	    M = phgAlloc((M_size = n) * sizeof(*M));
	    FreeAtExit(M);
	}
	for (k = 0; k < NFace; k++) {
	    v = DofFaceData(dof, e->faces[k]);
	    phgDofMapFaceData(dof->type, e, k, M);
	    for (j = 0; j < n; j++, i++)
		memcpy(v + M[j]*dim, value + i*dim, dim * sizeof(FLOAT));
	}
    }
    else if (n > 0) {
	for (k = 0; k < NFace; k++) {
	    v = DofFaceData(dof, e->faces[k]);
	    for (j = 0; j < n; j++, i++)
		memcpy(v + j*dim, value + i*dim, dim * sizeof(FLOAT));
	}
    }

    if ((n = dof->type->np_elem) > 1 && dof->type->DofMap != NULL &&
	!ElementIsInOrder(e)) {
	static int *M = NULL, M_size = 0;
	if (M_size < n) {
	    phgFree(M);
	    M = phgAlloc((M_size = n) * sizeof(*M));
	    FreeAtExit(M);
	}
	v = DofElementData(dof, e->index);
	phgDofMapElementData(dof->type, e, M);
	for (j = 0; j < n; j++, i++)
	    memcpy(v + M[j]*dim, value + i*dim, dim * sizeof(FLOAT));
    }
    else if (n > 0) {
	v = DofElementData(dof, e->index);
	for (j = 0; j < n; j++, i++)
	    memcpy(v + j*dim, value + i*dim, dim * sizeof(FLOAT));
    }

    return;
}

/*----------------- functions for exchanging neighbour data ---------------*/

/* Note:
 * 	These functions can only fetch face data from the neighbours,
 * 	and only work with certain DOF_TYPEs such as the DOF_DG types.
 *
 *	For more general cases the functions phgMap*NeighbourData (defined in
 *	map.c) can be used.
 */

NEIGHBOUR_DATA *
phgDofInitNeighbourData(DOF *dof, MAP *map)
/* prepares neighbour data on shared faces.
 *
 * Note:
 *	The data are stored in the order as return by phgDofGetBasesOnFace
 *	such that the ordering is identical on all elements.
 *
 *	If 'map' != NULL then the corresponding global vector indices of the
 *	data are stored in the array 'index'
 */
{
    NEIGHBOUR_DATA *nd;
    DOF_TYPE *t;
#if USE_MPI
    int nbas;
    GRID *g = dof->g;
    void *send, *recv;
    char *p = NULL;
    ELEMENT *e;
    size_t size;
    INT i, j, k, *stmp = NULL, *rtmp = NULL;
    MPI_Datatype type;
    int *scnts, *sdsps, *rcnts, *rdsps;
    RNEIGHBOUR *rn;
    SHORT *bases;
#endif	/* USE_MPI */
    FE_SPACE fe_space = (!DofIsHP(dof) ? dof->type->fe_space :
					 dof->hp->info->fe_space);

    /* Note: only DG type is allowed (to be fixed when needed) */
    if (fe_space != FE_L2 && fe_space != FE_None)
	phgError(1, "%s: don't know how to handle DOF \"%s\" "
		    "(only DG types are allowed).\n", __func__, dof->name);

    if ((t = dof->type) == NULL)
	t = dof->hp->info->types[dof->hp->max_order];

    if (t->base_type != NULL) {
	assert(t->dim == 1);
	assert(t == DOF_DGn[t->order]);
	t = t->base_type;
    }

    nd = phgCalloc(1, sizeof(*nd));
    CheckThread
    nd->cache = phgCalloc(phgMaxThreads, sizeof(*nd->cache));
    nd->dof = dof;

    if (map != NULL) {
	for (nd->dof_no = 0; nd->dof_no < map->ndof; nd->dof_no++)
	    if (map->dofs[nd->dof_no] == dof)
		break;
	assert(nd->dof_no < map->ndof);
	nd->map = map;
    }

    /* Note: the test below is to be removed in the future, after adding an
     * appropriate flag in DOF_TYPE indicating which bases functions are
     * nonzero on a given face */
    if (t != DOF_P0)
	nd->nbas = 3 * (t->np_vert + t->np_edge) + t->np_face;
    else
	nd->nbas = 1;		/* special case */

#if USE_MPI
    nbas = nd->nbas;		/* kept unchanged if !DofIsHP(dof) */

    if (g->nprocs <= 1)
	return nd;

/* Note: each data block in send/recv buffers contains:
 * 	nbas * FLOAT + (map == NULL ? 0 : nbas) * INT + 1 * ELEMENT + 2 * BYTE
 * and is padded according to sizeof(FLOAT). The following macros compute
 * offsets for each kind of data and depend on the current value of 'nbas' */
#define P_data(buf)	(buf)
#define P_index(buf)	(P_data(buf)	+ dof->dim * (nbas) * sizeof(FLOAT))
#define P_e(buf)	(P_index(buf)	+ (map==NULL ? 0 : nbas) * sizeof(INT))
#define P_op_vert(buf)	(P_e(buf)	+ sizeof(ELEMENT))
#define P_order(buf)	(P_op_vert(buf)	+ sizeof(BYTE))
#define P_end(buf)	(P_order(buf)	+ sizeof(BYTE))

    /* buffer for storing bases no. */
    bases = phgAlloc(nd->nbas * sizeof(*bases));

    if (!DofIsHP(dof)) {
	/* non hp DOF, nbas is constant */
	size = ((P_end(p) - p + sizeof(FLOAT) - 1) / sizeof (FLOAT))
		* sizeof(FLOAT);	/* Note: size will be kept unchanged */
	scnts = rcnts = g->neighbours.counts;
	sdsps = rdsps = g->neighbours.displs;
	send = phgAlloc(size * g->neighbours.count);
	recv = phgAlloc(size * g->neighbours.count);
	MPI_Type_contiguous(size, MPI_BYTE, &type);
	MPI_Type_commit(&type);
    }
    else {
	/* hp DOF, nbas varies */
	scnts = phgAlloc(4 * g->nprocs * sizeof(*scnts));
	sdsps = scnts + 1 * g->nprocs;
	rcnts = scnts + 2 * g->nprocs;
	rdsps = scnts + 3 * g->nprocs;
	memset(scnts, 0, g->nprocs * sizeof(*scnts));
	/* note: stmp/rtmp contains local/remote nbas values for each face */
	stmp = phgAlloc(2 * g->neighbours.count * sizeof(*stmp));
	rtmp = stmp + g->neighbours.count;
	/* store local nbas values in stmp[] */
	for (i = 0; i < g->neighbours.count; i++) {
	    rn = g->neighbours.list + i;
	    nbas = phgDofGetBasesOnFace(dof, rn->local, rn->vertex, NULL);
	    stmp[i] = nbas;
	    size = (P_end(p) - p + sizeof(FLOAT) - 1) / sizeof (FLOAT);
	    scnts[rn->rank] += size;
	}
	size = 0;
	for (i = 0; i < g->nprocs; i++) {
	    sdsps[i] = size;
	    size += scnts[i];
	}
	/* allocate send buffer */
	send = phgAlloc(size * sizeof(FLOAT));

	/* communicate nbas values */
	type = PHG_MPI_INT;
	MPI_Alltoallv(stmp, g->neighbours.counts, g->neighbours.displs, type,
		      rtmp, g->neighbours.counts, g->neighbours.displs, type,
		      g->comm);

	/* compute rcnts and rdsps */
	memset(rcnts, 0, g->nprocs * sizeof(*rcnts));
	for (i = 0; i < g->neighbours.count; i++) {
	    rn = g->neighbours.list + i;
	    nbas = rtmp[i];
	    size = (P_end(p) - p + sizeof(FLOAT) - 1) / sizeof (FLOAT);
	    rcnts[rn->rank] += size;
	}
	size = 0;
	for (i = 0; i < g->nprocs; i++) {
	    rdsps[i] = size;
	    size += rcnts[i];
	}
	/* allocate recv buffer */
	recv = phgAlloc(size * sizeof(FLOAT));

	type = PHG_MPI_FLOAT;
    }

    /* prepare send buffer */
    for (i = 0, p = send; i < g->neighbours.count; i++) {
	if (DofIsHP(dof)) {
	    nbas = stmp[i];
	    assert(nbas <= nd->nbas);
	    size = ((P_end(p) - p + sizeof(FLOAT) - 1) / sizeof (FLOAT))
		    * sizeof(FLOAT);
	}
	rn = g->neighbours.list + i;
	phgDofGetBasesOnFace(dof, rn->local, rn->vertex, bases);
	/* copy data to current send buffer block pointed by 'p' */
	for (j = 0; j < nbas; j++) {
	    INT k = phgDofMapE2D(dof, rn->local, bases[j] * dof->dim);
	    /* copy DOF data j */
	    memcpy(P_data(p) + j * dof->dim * sizeof(*dof->data),
				dof->data + k, dof->dim * sizeof(*dof->data));
	    if (map != NULL) {
		/* copy DOF index j */
		k = phgMapD2G(map, nd->dof_no, k);
		memcpy(P_index(p) + j * sizeof(k), &k, sizeof(k));
	    }
	}
	/* copy remote element pointer */
	memcpy(P_e(p), &rn->remote, sizeof(rn->remote));
	/* copy op_vertex value */
	*(BYTE *)P_op_vert(p) = rn->op_vertex;
	if (DofIsHP(dof))
	    *(BYTE *)P_order(p) = dof->hp->elem_order[rn->local->index];
	p += size;
    }
    phgFree(bases);

    MPI_Alltoallv(send, scnts, sdsps, type, recv, rcnts, rdsps, type, g->comm);

    phgFree(send);
    if (type != PHG_MPI_FLOAT)
	MPI_Type_free(&type);
    if (scnts != g->neighbours.counts)
	phgFree(scnts);

    /* copy and reorder received data */
    if (DofIsHP(dof)) {
	/* set up nd->loc[] (which is then needed when copying data) */
	nd->loc = phgAlloc((g->neighbours.count + 1) * sizeof(*nd->loc));
	nd->order = phgAlloc(g->neighbours.count * sizeof(*nd->order));
	for (i = 0, p = recv; i < g->neighbours.count; i++) {
	    nbas = rtmp[i];
	    size = ((P_end(p) - p + sizeof(FLOAT) - 1) / sizeof (FLOAT))
		    * sizeof(FLOAT);
	    /* use memcpy to set 'e' to avoid possible alignment problems */
	    memcpy(&e, P_e(p), sizeof(e));

	    j = GetRNeighbour(g, e, *(BYTE *)P_op_vert(p)) - g->neighbours.list;
	    nd->loc[j] = nbas;
	    p += size;
	}
	size = 0;
	for (i = 0; i < g->neighbours.count; i++) {
	    nbas = nd->loc[i];
	    nd->loc[i] = size;
	    size += nbas;
	}
	nd->loc[i] = size;
    }

    if (!DofIsHP(dof))
	j = g->neighbours.count * nd->nbas;
    else
	j = nd->loc[g->neighbours.count];
    nd->data = phgAlloc(sizeof(*nd->data) * dof->dim * j);
    if (map != NULL)
	nd->index = phgAlloc(sizeof(*nd->index) * j);

    /* pick up data in recv buffer */
    for (i = 0, p = recv; i < g->neighbours.count; i++) {
	if (DofIsHP(dof))
	    nbas = rtmp[i];
	/* use memcpy to set 'e' to avoid possible alignment problems */
	memcpy(&e, P_e(p), sizeof(e));
	j = GetRNeighbour(g, e, *(BYTE *)P_op_vert(p)) - g->neighbours.list;
	if (!DofIsHP(dof)) {
	    k = nbas * j;
	}
	else {
	    k = nd->loc[j];
	    size = ((P_end(p) - p + sizeof(FLOAT) - 1) / sizeof (FLOAT))
		    * sizeof(FLOAT);
	}
	memcpy(nd->data + k * dof->dim, P_data(p),
				nbas * dof->dim * sizeof(FLOAT));
	if (map != NULL)
	    memcpy(nd->index + k, P_index(p), nbas * sizeof(INT));
	if (DofIsHP(dof))
	    nd->order[j] = *(BYTE *)P_order(p);
	p += size;
    }
    phgFree(recv);
    phgFree(stmp);
#endif	/* USE_MPI */

    return nd;
}

void
phgDofReleaseNeighbourData(NEIGHBOUR_DATA **nd_ptr)
{
    NEIGHBOUR_DATA *nd;

    assert(nd_ptr != NULL);
    if ((nd = *nd_ptr) == NULL)
	return;

    phgFree(nd->data);
    phgFree(nd->index);
    phgFree(nd->loc);
    phgFree(nd->order);
    if (nd->cache != NULL) {
	int i;
	for (i = 0; i < phgMaxThreads; i++)
	    phgFree(nd->cache[i].plist);
	phgFree(nd->cache);
    }

    phgFree(nd);

    *nd_ptr = NULL;
}

FLOAT *
phgDofNeighbourData(NEIGHBOUR_DATA *nd, ELEMENT *e, int face_no, int bas_no,
		    INT *gindex)
/* returns the address of the DOF data for a given basis function on the face
 * 'face_no' of 'e' (note: bas_no is the basis number with respect to the face)
 *
 * If 'gindex' != NULL then it contains the global vector index of the DOF. */
{
    DOF *dof = nd->dof;
    ELEMENT *e1;
    int face_no1;

#if USE_MPI
    if (e->bound_type[face_no] & REMOTE) {
	/* Remote neighbour */
	INT j, k;

	j = (size_t)(e->neighbours[face_no]);
	if (!DofIsHP(dof)) {
	    assert(bas_no < nd->nbas);
	    k = j * nd->nbas + bas_no;
	}
	else {
	    k = nd->loc[j] + bas_no;
	    assert(k < nd->loc[j + 1]);
	}

	if (gindex != NULL) {
	    assert(nd->index != NULL);
	    *gindex = nd->index[k];
	}

	return nd->data + k * dof->dim;
    }
#endif

    e1 = e->neighbours[face_no];
    if (e1 == NULL)
	return NULL;
    face_no1 = phgOppositeFace(dof->g, e, face_no, e1);
    if (nd->cache[phgThreadId].peer != e1 ||
	nd->cache[phgThreadId].face != face_no1) {
	if (nd->cache[phgThreadId].plist == NULL) {
	    nd->cache[phgThreadId].plist = phgAlloc(nd->nbas
				* sizeof(*nd->cache[phgThreadId].plist));
	}
	nd->cache[phgThreadId].peer = e1;
	nd->cache[phgThreadId].face = face_no1;
	/* compute and cache list of bases on the face */
	phgDofGetBasesOnFace(dof, e1, face_no1, nd->cache[phgThreadId].plist);
    }

    if (gindex != NULL) {
	assert(nd->map != NULL);
	*gindex = phgMapE2G(nd->map, nd->dof_no, e1,
			nd->cache[phgThreadId].plist[bas_no] * dof->dim);
    }

    return dof->data + phgDofMapE2D(dof, e1,
			nd->cache[phgThreadId].plist[bas_no] * dof->dim);
}

int
phgDofNeighbourNBas(NEIGHBOUR_DATA *nd, ELEMENT *e, int face_no,
						DOF_TYPE **peer_type)
/* returns the number of bases on the given face in the peer element (the
 * value which would be returned by calling phgDofGetBasesOnFace with the
 * peer element), and optionally the DOF_TYPE on the peer element. */
{
    DOF *dof = nd->dof;
    ELEMENT *e1;
    int face_no1;

    if (!DofIsHP(dof)) {
	if (peer_type != NULL)
	    *peer_type = dof->type;
	return nd->nbas;
    }

#if USE_MPI
    if (e->bound_type[face_no] & REMOTE) {
	/* Remote neighbour */
	INT j = (size_t)(e->neighbours[face_no]);

	if (peer_type != NULL)
	    *peer_type = dof->hp->info->types[nd->order[j]];

	return nd->loc[j + 1] - nd->loc[j];
    }
#endif

    e1 = e->neighbours[face_no];
    assert(e1 != NULL);
    face_no1 = phgOppositeFace(dof->g, e, face_no, e1);

    if (peer_type != NULL)
	*peer_type = dof->hp->info->types[dof->hp->elem_order[e1->index]];

    return phgDofGetBasesOnFace(dof, e1, face_no1, NULL);
}
