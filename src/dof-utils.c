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

/* $Id: dof-utils.c,v 1.204 2022/09/05 06:50:14 zlb Exp $ */

#define NEED_GetVariableArgs

#include "phg.h"

#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <ctype.h>	/* toupper() */

DOF_TYPE *
phgDofCreateSimpleType(const char *name,
			int np_vert, int np_edge, int np_face, int np_elem)
/* Creates a simple DOF_TYPE with given np_{vert,edge,face,elem} and returns
 * a pointer to the DOF_TYPE created. */
{
    static DOF_TYPE template = {DofReserved, NULL,
	NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	NULL, NULL, NULL, NULL, FE_None,
	FALSE, FALSE, FALSE, -1, 0, 0, 0, -1, 1, 0, 0, 0, 0
    };
    DOF_TYPE *type = phgAlloc(sizeof(template));
    memcpy(type, &template, sizeof(template));
    type->name = strdup(name);
    type->np_vert = np_vert;
    type->np_edge = np_edge;
    type->np_face = np_face;
    type->np_elem = np_elem;
    type->nbas = type->np_vert * NVert + type->np_edge * NEdge +
		 type->np_face * NFace + type->np_elem * 1;

    return type;
}

FLOAT *
phgDofEval_(DOF *dof, ELEMENT *e, const FLOAT lambda[], const FLOAT *basvalues,
	    int bas_stride, FLOAT *values, FLOAT **data)
/* evaluates DOF value at the barycentric coordinates 'lambda'.
 * 'values' is an array of DofDim(dof) entries provided
 * by the calling function containing the evaluated values.
 *
 * returns 'values'
 *
 * 'basvalues', if not NULL, points to the array of values of the basis
 * functions at 'lambda' evaluated previously ('lambda' not used in this case)
 *
 * Values of the basis functions are stored as basvalues[nbas][bas_stride],
 * if bas_stride == 0 then it is set to type->dim.
 *
 * 'data', if not NULL, points to an array of pointers pointing to the DOF
 * data of the element (NVert + NEdge + NFace + 1 entries). This is mainly
 * for use by the interpolation functions.
 */
{
    int i, j, k, l, n, nvalues;
    FLOAT *u;
    const FLOAT *v, *a, *b;
    int dim = dof->dim, dim0 = DofTypeDim(dof);

    if (dof->type == DOF_CONSTANT) {
	memcpy(values, dof->data, dim * sizeof(*values));
	return values;
    }
    else if (dof->type == DOF_ANALYTIC) {
	if (lambda == NULL ||
	    (dof->userfunc == NULL && dof->userfunc_lambda == NULL)) {
	    phgError(1, "%s:%d, unexpected error.\n", __FILE__, __LINE__);
	}
	if (dof->userfunc_lambda == NULL) {
	    FLOAT x, y, z;
	    phgGeomLambda2XYZ(dof->g, e, lambda, &x, &y, &z);
	    dof->userfunc(x, y, z, values);
	}
	else {
	    dof->userfunc_lambda(dof, e, 0, lambda, values);
	}
	return values;
    }

    nvalues = dim * dim0;
    for (i = 0; i < nvalues; i++)
	values[i] = 0.0;

    if (DofIsHP(dof)) {
	/* hierarchical */
	int order;

	assert(dof->hp->info->types[dof->hp->elem_order[e->index]]->DofMap
		== NULL);

	/* compute values of basis functions */

	order = dof->hp->info->min_order;
	a = b = dof->hp->info->types[order]->BasFuncs(dof, e, 0, -1, lambda);

	/* vertices */
	for (k = 0; (dof->hp->info->flags & VERT_FLAG) && k < NVert; k++) {
	    order = dof->hp->vert_order[e->verts[k]];
	    n = dof->hp->info->types[order]->np_vert;
	    v = (data == NULL) ? DofVertexData(dof, e->verts[k]) : *(data++);
	    for (j = 0; j < n; j++) {
		u = values;
		for (i = 0; i < dim; i++, v++) {
		    a = b;
		    for (l = 0; l < dim0; l++) {
			*(u++) += (*(a++)) * (*v);
		    }
		}
		b = a;
	    }
	}

	/* edges */
	/* Note: no need to call phgDofMapEdgeData for hierarchical bases. */
	for (k = 0; (dof->hp->info->flags & EDGE_FLAG) && k < NEdge; k++) {
	    order = dof->hp->edge_order[e->edges[k]];
	    n = dof->hp->info->types[order]->np_edge;
	    v = (data == NULL) ? DofEdgeData(dof, e->edges[k]) : *(data++);
	    for (j = 0; j < n; j++) {
		u = values;
		for (i = 0; i < dim; i++, v++) {
		    a = b;
		    for (l = 0; l < dim0; l++) {
			*(u++) += (*(a++)) * (*v);
		    }
		}
		b = a;
	    }
	}

	/* faces */
	/* Note: no need to call phgDofMapFaceData for hierarchical bases. */
	for (k = 0; (dof->hp->info->flags & FACE_FLAG) && k < NFace; k++) {
	    order = dof->hp->face_order[e->faces[k]];
	    n = dof->hp->info->types[order]->np_face;
	    v = (data == NULL) ? DofFaceData(dof, e->faces[k]) : *(data++);
	    for (j = 0; j < n; j++) {
		u = values;
		for (i = 0; i < dim; i++, v++) {
		    a = b;
		    for (l = 0; l < dim0; l++) {
			*(u++) += (*(a++)) * (*v);
		    }
		}
		b = a;
	    }
	}

	/* element */
	/* Note: no need to call phgDofMapElementData for hierarchical bases. */
	order = dof->hp->elem_order[e->index];
	n = dof->hp->info->types[order]->np_elem;
	v = (data == NULL) ? DofElementData(dof, e->index) : *(data++);
	for (j = 0; j < n; j++) {
	    u = values;
	    for (i = 0; i < dim; i++, v++) {
		a = b;
		for (l = 0; l < dim0; l++) {
		    *(u++) += (*(a++)) * (*v);
		}
	    }
	    b = a;
	}

	return values;
    }

    assert(!DofIsHP(dof));

    if (dof->type->BasFuncs == NULL)
	phgError(1, "%s:%d, function 'BasFuncs' unvailable in \"%s\"!\n",
		 __FILE__, __LINE__, dof->type->name);

    /* get basis function values, adjust bas_stride -= type->dim */
    if (basvalues == NULL) {
	a = b = dof->type->BasFuncs(dof, e, 0, -1, lambda);
	bas_stride = 0;
    }
    else {
	a = b = basvalues;
	if (bas_stride <= 0)
	    bas_stride = 0;
	else
	    bas_stride -= dim0;
    }

    if ((n = dof->type->np_vert) > 0) {
	for (k = 0; k < NVert; k++) {
	    v = (data == NULL) ? DofVertexData(dof, e->verts[k]) : data[k];
	    for (j = 0; j < n; j++) {
		u = values;
		for (i = 0; i < dim; i++, v++) {
		    a = b;
		    for (l = 0; l < dim0; l++) {
			*(u++) += (*(a++)) * (*v);
		    }
		}
		b = a + bas_stride;
	    }
	}
    }
    if (data != NULL)
	data += NVert;

    if ((n = dof->type->np_edge) > 1) {
	int M[n];
	for (k = 0; k < NEdge; k++) {
	    v = (data == NULL) ? DofEdgeData(dof, e->edges[k]) : data[k];
	    phgDofMapEdgeData(dof->type, e, k, M);
	    for (j = 0; j < n; j++) {
		u = values;
		for (i = 0; i < dim; i++) {
		    a = b;
		    for (l = 0; l < dim0; l++) {
			*(u++) += (*(a++)) * v[M[j] * dof->dim + i];
		    }
		}
		b = a + bas_stride;
	    }
	}
    }
    else if (n > 0) {
	for (k = 0; k < NEdge; k++) {
	    v = (data == NULL) ? DofEdgeData(dof, e->edges[k]) : data[k];
	    for (j = 0; j < n; j++) {
		u = values;
		for (i = 0; i < dim; i++, v++) {
		    a = b;
		    for (l = 0; l < dim0; l++) {
			*(u++) += (*(a++)) * (*v);
		    }
		}
		b = a + bas_stride;
	    }
	}
    }
    if (data != NULL)
	data += NEdge;

    if ((n = dof->type->np_face) > 1) {
	int M[n];
	for (k = 0; k < NFace; k++) {
	    v = (data == NULL) ? DofFaceData(dof, e->faces[k]) : data[k];
	    phgDofMapFaceData(dof->type, e, k, M);
	    for (j = 0; j < n; j++) {
		u = values;
		for (i = 0; i < dim; i++) {
		    a = b;
		    for (l = 0; l < dim0; l++) {
			*(u++) += (*(a++)) * v[M[j] * dof->dim + i];
		    }
		}
		b = a + bas_stride;
	    }
	}
    }
    else if (n > 0) {
	for (k = 0; k < NFace; k++) {
	    v = (data == NULL) ? DofFaceData(dof, e->faces[k]) : data[k];
	    for (j = 0; j < n; j++) {
		u = values;
		for (i = 0; i < dim; i++, v++) {
		    a = b;
		    for (l = 0; l < dim0; l++) {
			*(u++) += (*(a++)) * (*v);
		    }
		}
		b = a + bas_stride;
	    }
	}
    }
    if (data != NULL)
	data += NFace;

    if ((n = dof->type->np_elem) > 1 && dof->type->DofMap != NULL &&
	!ElementIsInOrder(e)) {
	int M[n];
	v = (data == NULL) ? DofElementData(dof, e->index) : *data;
	phgDofMapElementData(dof->type, e, M);
	for (j = 0; j < n; j++) {
	    u = values;
	    for (i = 0; i < dim; i++) {
		a = b;
		for (l = 0; l < dim0; l++) {
		    *(u++) += (*(a++)) * v[M[j] * dof->dim + i];
		}
	    }
	    b = a + bas_stride;
	}
    }
    else if (n > 0) {
	v = (data == NULL) ? DofElementData(dof, e->index) : *data;
	for (j = 0; j < n; j++) {
	    u = values;
	    for (i = 0; i < dim; i++, v++) {
		a = b;
		for (l = 0; l < dim0; l++) {
		    *(u++) += (*(a++)) * (*v);
		}
	    }
	    b = a + bas_stride;
	}
    }

    return values;
}

static void
get_bas_funcs(DOF *src, DOF *dof, ELEMENT *e, FLOAT *values)
/* computes values of basis funcs of 'src' at nodal points of 'dof'.
 * values points to an array of
	dof->type->nbas * src->type->nbas * src->type->dim
 * entries to hold the results */
{
    int count, j, k, i0, i1, i2;
    FLOAT *points, *v, *p;
    FLOAT lambda[Dim + 1] = { 0., 0., 0., 0. };
    DOF_TYPE *t = dof->type;

    if (e == NULL)
	return;

    while (t->points == NULL) {
	assert(t->base_type != NULL);
	t = t->base_type;
    }
    points = t->points;

    count = src->type->nbas * src->type->dim;
    k = t->nbas * count;
    v = values;

    if (t->np_vert > 0) {
	for (k = 0; k < NVert; k++) {
	    p = points;
	    for (j = 0; j < t->np_vert; j++) {
		lambda[k] = *p;
		memcpy(v, src->type->BasFuncs(src, e, 0, -1, lambda),
		       count * sizeof(*values));
		p += 1;
		v += count;
	    }
	    lambda[k] = 0.;
	}
	points = p;
    }

    if (t->np_edge > 0) {
	for (k = 0; k < NEdge; k++) {
	    p = points;
	    i0 = GetEdgeVertex(k, 0);
	    i1 = GetEdgeVertex(k, 1);
	    for (j = 0; j < t->np_edge; j++) {
		lambda[i0] = p[0];
		lambda[i1] = p[1];
		memcpy(v, src->type->BasFuncs(src, e, 0, -1, lambda),
		       count * sizeof(*values));
		p += 2;
		v += count;
	    }
	    lambda[i0] = lambda[i1] = 0.;
	}
	points = p;
    }

    if (t->np_face > 0) {
	for (k = 0; k < NFace; k++) {
	    p = points;
	    i0 = GetFaceVertex(k, 0);
	    i1 = GetFaceVertex(k, 1);
	    i2 = GetFaceVertex(k, 2);
	    for (j = 0; j < t->np_face; j++) {
		lambda[i0] = p[0];
		lambda[i1] = p[1];
		lambda[i2] = p[2];
		memcpy(v, src->type->BasFuncs(src, e, 0, -1, lambda),
		       count * sizeof(*values));
		p += 3;
		v += count;
	    }
	    lambda[i0] = lambda[i1] = lambda[i2] = 0.;
	}
	points = p;
    }

    if (t->np_elem > 0) {
	for (j = 0; j < t->np_elem; j++) {
	    memcpy(v, src->type->BasFuncs(src, e, 0, -1, points),
		   count * sizeof(*values));
	    points += 4;
	    v += count;
	}
    }
}

static FLOAT *
get_bas_gradient(DOF *src, DOF *dof, ELEMENT *e)
/* computes gradient of basis funcs of 'src' (with respect to lambda)
 * at nodal points of 'dof'.
 * returns a pointer to a buffer containing the computed values,
 * which should be freed by the calling functions */
{
    int grad_count, j, k, i0, i1, i2;
    FLOAT *points, *grad_values, *g, *p;
    FLOAT lambda[Dim + 1] = { 0., 0., 0., 0. };
    DOF_TYPE *t = dof->type;

    if (e == NULL)
	return NULL;

    while (t->points == NULL) {
	assert(t->base_type != NULL);
	t = t->base_type;
    }
    points = t->points;


    assert(src->type->invariant == TRUE);

    grad_count = src->type->nbas * (Dim + 1) * src->type->dim;
    k = t->nbas * grad_count;
    g = grad_values = phgAlloc(k * sizeof(*grad_values));

    if (t->np_vert > 0) {
	for (k = 0; k < NVert; k++) {
	    p = points;
	    for (j = 0; j < t->np_vert; j++) {
		lambda[k] = *p;
		memcpy(g, src->type->BasGrads(src, e, 0, -1, lambda),
		       grad_count * sizeof(*grad_values));
		p += 1;
		g += grad_count;
	    }
	    lambda[k] = 0.;
	}
	points = p;
    }

    if (t->np_edge > 0) {
	for (k = 0; k < NEdge; k++) {
	    p = points;
	    i0 = GetEdgeVertex(k, 0);
	    i1 = GetEdgeVertex(k, 1);
	    for (j = 0; j < t->np_edge; j++) {
		lambda[i0] = p[0];
		lambda[i1] = p[1];
		memcpy(g, src->type->BasGrads(src, e, 0, -1, lambda),
		       grad_count * sizeof(*grad_values));
		p += 2;
		g += grad_count;
	    }
	    lambda[i0] = lambda[i1] = 0.;
	}
	points = p;
    }

    if (t->np_face > 0) {
	for (k = 0; k < NFace; k++) {
	    p = points;
	    i0 = GetFaceVertex(k, 0);
	    i1 = GetFaceVertex(k, 1);
	    i2 = GetFaceVertex(k, 2);
	    for (j = 0; j < t->np_face; j++) {
		lambda[i0] = p[0];
		lambda[i1] = p[1];
		lambda[i2] = p[2];
		memcpy(g, src->type->BasGrads(src, e, 0, -1, lambda),
		       grad_count * sizeof(*grad_values));
		p += 3;
		g += grad_count;
	    }
	    lambda[i0] = lambda[i1] = lambda[i2] = 0.;
	}
	points = p;
    }

    if (t->np_elem > 0) {
	for (j = 0; j < t->np_elem; j++) {
	    memcpy(g, src->type->BasGrads(src, e, 0, -1, points),
		   grad_count * sizeof(*grad_values));
	    points += 4;
	    g += grad_count;
	}
    }

    return grad_values;
}

static FLOAT *
get_grad_lambda(DOF *dof, ELEMENT *e, const FLOAT lambda[],
		const FLOAT *gradbas, FLOAT **data)
/* evaluates the gradient with respect to lambda at the barycentric
 * coordinates 'lambda'. returns a pointer to a buffer containing the
 *		[dof->dim][dof->type->dim][Dim+1]
 * evaluated values, which should be freed by the calling function.
 *
 * 'gradbas', if not NULL, points to the array of gradients of the basis
 * functions at 'lambda' evaluated previously ('lambda' not used in this case)
 *
 * 'data', if not NULL, points to an array of pointers pointing to the DOF
 * data of the element (NVert + NEdge + NFace + 1 entries). This is mainly
 * for use by the interpolation functions.
 */
{
    int i, j, k, l, m, dim;
    const FLOAT *b;
    FLOAT *v, *grad_lambda, *grad, *g;
    DOF_TYPE *type0 = dof->type, *type = type0;

    if (e == NULL)
	return NULL;

    if (type0 == NULL)
	type0 = dof->hp->info->types[dof->hp->max_order];

    dim = DofDim(dof);
    grad = grad_lambda = phgAlloc(dim * (Dim + 1) * sizeof(*grad_lambda));
    if (gradbas == NULL) {
	assert(lambda != NULL);
	gradbas = type0->BasGrads(dof, e, 0, -1, lambda);
    }

    for (i = 0; i < dim * (Dim + 1); i++)
	grad_lambda[i] = 0.0;

    for (l = 0; l < Dim + 1; l++) {

	b = gradbas++;

	if (type0->np_vert > 0) {
	    type = type0;
	    for (k = 0; k < NVert; k++) {
		type = DofTypeOnVertex(dof, e->verts[k]);
		v = (data == NULL) ? DofVertexData(dof, e->verts[k]) : data[k];
		for (j = 0; j < type->np_vert; j++) {
		    g = grad;
		    for (i = 0; i < dof->dim; i++) {
			for (m = 0; m < type->dim; m++) {
			    *g += b[m * (Dim + 1)] * (*v);
			    g += Dim + 1;
			}
			v++;
		    }
		    b += type->dim * (Dim + 1);
		}
	    }
	}
	if (data != NULL)
	    data += NVert;

	if (type0->np_edge > 1) {
	    for (k = 0; k < NEdge; k++) {
		INT eno = e->edges[k];
		int M[type0->np_edge];
		if (DofIsHP(dof))
		    type = dof->hp->info->types[dof->hp->edge_order[eno]];
		v = (data == NULL) ? DofEdgeData(dof, eno) : data[k];
		phgDofMapEdgeData(type, e, k, M);
		for (j = 0; j < type->np_edge; j++) {
		    g = grad;
		    for (i = 0; i < dof->dim; i++) {
			for (m = 0; m < type->dim; m++) {
			    *g += b[m * (Dim + 1)] * v[M[j] * dof->dim + i];
			    g += Dim + 1;
			}
		    }
		    b += type->dim * (Dim + 1);
		}
	    }
	}
	else if (type0->np_edge > 0) {
	    for (k = 0; k < NEdge; k++) {
		INT eno = e->edges[k];
		if (DofIsHP(dof))
		    type = dof->hp->info->types[dof->hp->edge_order[eno]];
		v = (data == NULL) ? DofEdgeData(dof, eno) : data[k];
		for (j = 0; j < type->np_edge; j++) {
		    g = grad;
		    for (i = 0; i < dof->dim; i++) {
			for (m = 0; m < type->dim; m++) {
			    *g += b[m * (Dim + 1)] * (*v);
			    g += Dim + 1;
			}
			v++;
		    }
		    b += type->dim * (Dim + 1);
		}
	    }
	}
	if (data != NULL)
	    data += NEdge;

	if (type0->np_face > 1) {
	    for (k = 0; k < NFace; k++) {
		INT fno = e->faces[k];
		int M[type0->np_face];
		if (DofIsHP(dof))
		    type = dof->hp->info->types[dof->hp->face_order[fno]];
		v = (data == NULL) ? DofFaceData(dof, fno) : data[k];
		phgDofMapFaceData(type, e, k, M);
		for (j = 0; j < type->np_face; j++) {
		    g = grad;
		    for (i = 0; i < dof->dim; i++) {
			for (m = 0; m < type->dim; m++) {
			    *g += b[m * (Dim + 1)] * v[M[j] * dof->dim + i];
			    g += Dim + 1;
			}
		    }
		    b += type->dim * (Dim + 1);
		}
	    }
	}
	else if (type0->np_face > 0) {
	    for (k = 0; k < NFace; k++) {
		INT fno = e->faces[k];
		if (DofIsHP(dof))
		    type = dof->hp->info->types[dof->hp->face_order[fno]];
		v = (data == NULL) ? DofFaceData(dof, fno) : data[k];
		for (j = 0; j < type->np_face; j++) {
		    g = grad;
		    for (i = 0; i < dof->dim; i++) {
			for (m = 0; m < type->dim; m++) {
			    *g += b[m * (Dim + 1)] * (*v);
			    g += Dim + 1;
			}
			v++;
		    }
		    b += type->dim * (Dim + 1);
		}
	    }
	}
	if (data != NULL)
	    data += NFace;

	if (type0->np_elem > 1 && type0->DofMap != NULL &&
	    !ElementIsInOrder(e)) {
	    int M[type0->np_elem];
	    if (DofIsHP(dof))
		type = dof->hp->info->types[dof->hp->elem_order[e->index]];
	    v = (data == NULL) ? DofElementData(dof, e->index) : *data;
	    phgDofMapElementData(type, e, M);
	    for (j = 0; j < type->np_elem; j++) {
		g = grad;
		for (i = 0; i < dof->dim; i++) {
		    for (m = 0; m < type->dim; m++) {
			*g += b[m * (Dim + 1)] * v[M[j] * dof->dim + i];
			g += Dim + 1;
		    }
		}
		b += type->dim * (Dim + 1);
	    }
	}
	else if (type0->np_elem > 0) {
	    if (DofIsHP(dof))
		type = dof->hp->info->types[dof->hp->elem_order[e->index]];
	    v = (data == NULL) ? DofElementData(dof, e->index) : *data;
	    for (j = 0; j < type->np_elem; j++) {
		g = grad;
		for (i = 0; i < dof->dim; i++) {
		    for (m = 0; m < type->dim; m++) {
			*g += b[m * (Dim + 1)] * (*v);
			g += Dim + 1;
		    }
		    v++;
		}
		b += type->dim * (Dim + 1);
	    }
	}
	if (data != NULL)
	    data -= NVert + NEdge + NFace;	/* restore data */

	grad++;
    }

    return grad_lambda;
}

FLOAT *
phgDofEvalGradient_(DOF *dof, ELEMENT *e, const FLOAT lambda[],
		   const FLOAT *gradbas, FLOAT *values, FLOAT **data)
/* evaluates the gradient at the barycentric coordinates 'lambda'.
 * 'values' is an array, provided by the calling function, of
 * [DofDim(dof)][Dim] entries containing the evaluated values.
 *
 * returns 'values'
 *
 * 'gradbas', if not NULL, points to the array of gradients of the basis
 * functions at 'lambda' evaluated previously ('lambda' not used in this case)
 *
 * 'data', if not NULL, points to an array of pointers pointing to the DOF
 * data of the element (NVert + NEdge + NFace + 1 entries). This is mainly
 * for use by the interpolation functions.
 */
{
    int i, j, k, nvalues;
    FLOAT *grad_lambda, *g;
    const FLOAT (*J)[Dim + 1];

    if (SpecialDofType(dof->type))
	phgError(1, "%s:%d, cannot evaluate gradient of special DOF types.\n");

    if (!DofIsHP(dof) && dof->type->BasGrads == NULL)
	phgError(1, "%s:%d, function 'BasGrads' unvailable in \"%s\"!\n",
		 __FILE__, __LINE__, dof->type->name);

    nvalues = DofDim(dof);
    J = (void *)(phgGeomGetJacobian(dof->g, e));

    /* get gradient with respect to lambda */
    grad_lambda = get_grad_lambda(dof, e, lambda, gradbas, data);

    /* multiply grad_lambda by the Jacobian to get grad_xyz */
    for (i = 0; i < nvalues * Dim; i++)
	values[i] = 0.0;

    g = grad_lambda;
    for (k = 0; k < Dim + 1; k++) {
	for (j = 0; j < Dim; j++) {
	    for (i = 0; i < nvalues; i++) {
		values[i * Dim + j] += g[i * (Dim + 1) + k] * J[k][j];
	    }
	}
    }
    phgFree(grad_lambda);

    return values;
}

FLOAT *
phgDofEvalDivergence(DOF *dof, ELEMENT *e, const FLOAT lambda[],
		     const FLOAT *gradbas, FLOAT *values)
/* evaluates the divergence at the barycentric coordinates 'lambda'.
 * 'values' is an array, provided by the calling function, of
 * [DofDim(dof) / Dim] entries containing the evaluated values.
 *
 * returns 'values'
 *
 * 'gradbas', if not NULL, points to the array of gradients of the basis
 * functions 'lambda' evaluated previously ('lambda' not used in this case)
 */
{
    int i, j, k, dim;
    FLOAT *grad_lambda, *g, d;
    FLOAT (*J)[Dim + 1];

    if (SpecialDofType(dof->type))
	phgError(1,
		 "%s:%d, can't evaluate divergence of special DOF types.\n");

    assert(DofDim(dof) % Dim == 0);

    if (!DofIsHP(dof) && dof->type->BasGrads == NULL)
	phgError(1, "%s:%d, function 'BasGrads' unvailable in \"%s\"!\n",
		 __FILE__, __LINE__, dof->type->name);

    J = (void *)(phgGeomGetJacobian(dof->g, e));

    /* get gradient with respect to lambda */
    grad_lambda = get_grad_lambda(dof, e, lambda, gradbas, NULL);

    /* multiply grad_lambda with columns of the Jacobian to compute div_xyz */
    dim = DofDim(dof) / Dim;
    g = grad_lambda;
    for (i = 0; i < dim; i++) {
	d = 0.0;
	for (j = 0; j < Dim; j++)
	    for (k = 0; k < Dim + 1; k++)
		d += g[(i * Dim + j) * (Dim + 1) + k] * J[k][j];
	values[i] = d;
    }

    phgFree(grad_lambda);

    return values;
}

FLOAT *
phgDofEvalCurl(DOF *dof, ELEMENT *e, const FLOAT lambda[],
	       const FLOAT *gradbas, FLOAT *values)
/* evaluates the Curl at the barycentric coordinates 'lambda'.
 * 'values' is an array, provided by the calling function, of
 * [DofDim(dof) / Dim] entries containing the evaluated values.
 *
 * returns 'values'
 *
 * 'gradbas', if not NULL, points to the array of gradients of the basis
 * functions 'lambda' evaluated previously ('lambda' not used in this case)
 */
{
    int i, j, k, m, dim;
    FLOAT *grad_lambda, *grad_xyz, *g, d;
    FLOAT (*J)[Dim + 1];

    if (SpecialDofType(dof->type))
	phgError(1, "%s:%d, can't evaluate curl of special DOF types.\n");

    assert(DofDim(dof) % Dim == 0);

    if (!DofIsHP(dof) && dof->type->BasGrads == NULL)
	phgError(1, "%s:%d, function 'BasGrads' unvailable in \"%s\"!\n",
		 __FILE__, __LINE__, dof->type->name);

    J = (void *)(phgGeomGetJacobian(dof->g, e));

    /* get gradient with respect to lambda */
    grad_lambda = get_grad_lambda(dof, e, lambda, gradbas, NULL);

    /* multiply grad_lambda with columns of the Jacobian to compute grad_xyz */
    grad_xyz = phgAlloc(DofDim(dof) * Dim * sizeof(*grad_xyz));

    dim = DofDim(dof) / Dim;
    g = grad_lambda;
    for (i = 0; i < dim; i++) {
	for (m = 0; m < Dim; m++)
	    for (j = 0; j < Dim; j++) {
		d = 0.0;
		for (k = 0; k < Dim + 1; k++)
		    d += g[m * (Dim + 1) + k] * J[k][j];
		grad_xyz[m * Dim + j] = d;
	    }
	g += Dim * (Dim + 1);
    }
    /* compute curl_xyz using grad_xyz */
    g = grad_xyz;
    for (i = 0; i < dim; i++) {
	*(values++) = g[7] - g[5];
	*(values++) = g[2] - g[6];
	*(values++) = g[3] - g[1];
	g += (Dim * Dim);
    }

    phgFree(grad_lambda);
    phgFree(grad_xyz);

    return values;
}

static DOF *cache_dof;
static FLOAT *cache_buf = NULL;
static int cache_dim = 0;

static void
grad_func(DOF *dof, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* temporary Init function used by phgDofGradient */
{
    CheckD_order(dof);
    phgDofEvalGradient(cache_dof, e, lambda, cache_buf + bno * cache_dim,
				values);
}

DOF *
phgDofGradient_(DOF *src, DOF **dest_ptr, DOF_TYPE *newtype,
		const char *name, const char *srcfile, int srcline)
/* returns a new DOF which is the gradient of 'src' */
{
    char *s = NULL;
    DOF *dest = (dest_ptr == NULL ? NULL : *dest_ptr);
    HP_TYPE *hp = NULL;
    int dim;

    assert(!DofIsHP(src) || newtype == NULL);

    MagicCheck(DOF, dest)

    if (dest != NULL)
	phgDofFree(&dest);

    if (DofIsHP(src)) {
	assert(newtype == NULL);
	hp = phgHPGetDerivative_(src->hp, srcfile, srcline);
    }
    else if (newtype == NULL) {
	assert(!DofIsHP(src));
	newtype = src->type->grad_type;
	if (newtype == NULL)
	    phgError(1, "%s: gradient type for DOF \"%s\" undefined.\n",
		 __func__, src->type->name);
    }

    /* create the new DOF */
    if (name != NULL) {
	s = phgAlloc(strlen(name) + 1);
	strcpy(s, name);
    }
    else {
	s = phgAlloc(strlen(src->name) + 1 + 6);
	sprintf(s, "grad(%s)", src->name);
    }

    dim = (newtype != NULL ? newtype->dim : hp->info->dim);
    dest = phgDofNew_(src->g, newtype, hp, (Dim * DofDim(src)) / dim,
		      s, DofNoAction, srcfile, srcline);
    dest->DB_mask = 0;
    if (dest->DB_masks != NULL)
	memset(dest->DB_masks, 0, dest->dim * sizeof(*dest->DB_masks));

    phgHPFree(&hp);
    phgFree(s);

    if (dest_ptr != NULL)
	*dest_ptr = dest;

    if (src->data == NULL)
	return dest;

    if (!DofIsHP(src) && src->type->invariant == TRUE && dest->type->is_nodal) {
	cache_buf = get_bas_gradient(src, dest, src->g->roots);
	cache_dim = src->type->nbas * (Dim + 1) * src->type->dim;
    }
    cache_dof = src;

    phgDofSetDataByFunction_(dest, NULL, grad_func, NULL);

    if (cache_buf != NULL) {
	phgFree(cache_buf);
	cache_buf = NULL;
    }
    cache_dim = 0;

    return dest;
}

static void
div_func(DOF *dof, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* temporary Init function used by phgDofDivergence */
{
    CheckD_order(dof);
    phgDofEvalDivergence(cache_dof, e, lambda, cache_buf + bno * cache_dim,
			 values);
}

DOF *
phgDofDivergence_(DOF *src, DOF **dest_ptr, DOF_TYPE *newtype,
			const char *name, const char *srcfile, int srcline)
/* returns a new DOF which is the divergence of 'src' */
{
    char *s = NULL;
    DOF *dest = (dest_ptr == NULL ? NULL : *dest_ptr);
    HP_TYPE *hp = NULL;
    int dim;

    assert(!DofIsHP(src) || newtype == NULL);
    assert(DofDim(src) % Dim == 0);

    MagicCheck(DOF, dest)

    if (dest != NULL)
	phgDofFree(&dest);

    if (DofIsHP(src)) {
	assert(newtype == NULL);
	hp = phgHPGetDerivative_(src->hp, srcfile, srcline);
    }
    else if (newtype == NULL) {
	assert(!DofIsHP(src));
	newtype = src->type->grad_type;
	if (newtype == NULL)
	    phgError(1, "%s: gradient type for DOF \"%s\" undefined.\n",
		 __func__, src->type->name);
    }

    /* create the new DOF */
    if (name != NULL) {
	s = phgAlloc(strlen(name) + 1);
	strcpy(s, name);
    }
    else {
	s = phgAlloc(strlen(src->name) + 1 + 6);
	sprintf(s, "div(%s)", src->name);
    }

    /* assert(DofDim(src) % (Dim * newtype->dim) == 0);*/
    dim = (newtype != NULL ? newtype->dim : hp->info->dim);
    dest = phgDofNew_(src->g, newtype, hp, DofDim(src) / (Dim * dim),
		      s, DofNoAction, srcfile, srcline);
    dest->DB_mask = 0;
    if (dest->DB_masks != NULL)
	memset(dest->DB_masks, 0, dest->dim * sizeof(*dest->DB_masks));
    phgHPFree(&hp);
    phgFree(s);

    if (dest_ptr != NULL)
	*dest_ptr = dest;

    if (src->data == NULL)
	return dest;

    if (!DofIsHP(src) && src->type->invariant == TRUE && dest->type->is_nodal) {
	cache_buf = get_bas_gradient(src, dest, src->g->roots);
	cache_dim = src->type->nbas * (Dim + 1) * src->type->dim;
    }
    cache_dof = src;

    phgDofSetDataByFunction_(dest, NULL, div_func, NULL);

    if (cache_buf != NULL) {
	phgFree(cache_buf);
	cache_buf = NULL;
    }
    cache_dim = 0;

    return dest;
}

static void
curl_func(DOF *dof, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* temporary Init function used by phgDofCurl */
{
    CheckD_order(dof);
    phgDofEvalCurl(cache_dof, e, lambda, cache_buf + bno * cache_dim, values);
}

DOF *
phgDofCurl_(DOF *src, DOF **dest_ptr, DOF_TYPE *newtype,
		const char *name, const char *srcfile, int srcline)
{
    char *s;
    DOF *dest = (dest_ptr == NULL ? NULL : *dest_ptr);
    HP_TYPE *hp = NULL;
    int dim;

    assert(!DofIsHP(src) || newtype == NULL);
    assert(DofDim(src) % Dim == 0);

    MagicCheck(DOF, dest)

    if (dest != NULL)
	phgDofFree(&dest);

    if (DofIsHP(src)) {
	assert(newtype == NULL);
	hp = phgHPGetDerivative_(src->hp, srcfile, srcline);
    }
    else if (newtype == NULL) {
	assert(!DofIsHP(src));
	newtype = src->type->grad_type;
	if (newtype == NULL)
	    phgError(1, "%s: gradient type for DOF \"%s\" undefined.\n",
		 __func__, src->type->name);
    }

    /* create the new DOF */
    if (name != NULL) {
	s = phgAlloc(strlen(name) + 1);
	strcpy(s, name);
    }
    else {
	s = phgAlloc(strlen(src->name) + 1 + 6);
	sprintf(s, "curl(%s)", src->name);
    }

    dim = (newtype != NULL ? newtype->dim : hp->info->dim);
    assert(DofDim(src) % dim == 0);
    dest = phgDofNew_(src->g, newtype, hp, DofDim(src) / dim,
		      s, DofNoAction, srcfile, srcline);
    dest->DB_mask = 0;
    if (dest->DB_masks != NULL)
	memset(dest->DB_masks, 0, dest->dim * sizeof(*dest->DB_masks));
    phgHPFree(&hp);
    phgFree(s);

    if (dest_ptr != NULL)
	*dest_ptr = dest;

    if (src->data == NULL)
	return dest;

#if 0
    /* FIXME: incorrect (see test/maxwell-ptest.c) */
    if (newtype->InitFunc == phgDofInitFuncPoint &&
	src->type->InitFunc == phgDofInitFuncPoint) {
	/* simpler code when grad_type is a point function */
	FLOAT *g, *a;
	INT k, n;
	DOF *grad;

	grad = phgDofGradient(src, NULL, NULL, NULL);
	g = DofData(grad);
	a = DofData(dest);
	n = DofDataCount(dest);
	for (k = 0; k < n; k += Dim) {
	    *(a++) = g[7] - g[5];
	    *(a++) = g[2] - g[6];
	    *(a++) = g[3] - g[1];
	    g += (Dim * Dim);
	}
	phgDofFree(&grad);
	return dest;
    }
#endif

    if (!DofIsHP(src) && src->type->invariant == TRUE &&
	dest->type->InitFunc == phgDofInitFuncPoint) {
	cache_buf = get_bas_gradient(src, dest, src->g->roots);
	cache_dim = src->type->nbas * (Dim + 1) * src->type->dim;
    }
    cache_dof = src;

    phgDofSetDataByFunction_(dest, NULL, curl_func, NULL);

    if (cache_buf != NULL) {
	phgFree(cache_buf);
	cache_buf = NULL;
    }
    cache_dim = 0;
    cache_dof = NULL;

    return dest;
}

/* static variables for the Hessian matrix */
static DOF *hessian_dof = NULL;
#if USE_OMP
#pragma omp threadprivate(hessian_dof)
#endif  /* USE_OMP */

static void
hessian_func(DOF *u, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* returns the gradient of all the basis functions of 'hessian_dof',
 * in the order values[type->nbas][type->dim][Dim] */
{
    int i, j;
    const FLOAT *p = hessian_dof->type->BasGrads(hessian_dof, e, 0, -1, lambda);
    const FLOAT (*J)[Dim + 1] = (void *)phgGeomGetJacobian(u->g, e);

    for (i=0; i<hessian_dof->type->nbas*hessian_dof->type->dim; i++, p+=Dim+1)
	for (j = 0; j < Dim; j++) 
	    *(values++) = p[0]*J[0][j]+p[1]*J[1][j]+p[2]*J[2][j]+p[3]*J[3][j];
}

void
phgDofEvalBasisHessian(DOF *dof, ELEMENT *e, const FLOAT *lambda,
			const char *which, FLOAT *values)
/* returns in 'values' the second derivatives, i.e., the Hessian matrix, of
 * all the basis functions of 'dof' on the element 'e' at the point with
 * barycoordinates 'lambda[]', in the order:
 * 		values[dof->type->nbas][dof->type->dim][N]
 *
 * The argument 'which' specifies the part or form of the matrix returned:
 *    -	toupper(*which) == 'D': returns only the diagonal entries:
 *   		\phi_xx, \phi_yy, \phi_zz (N = 3)
 *    -	toupper(*which) == 'L': returns the lower triangular matrix:
 *		\phi_xx, \phi_yx, \phi_yy, \phi_zx, \phi_zy, \phi_zz (N = 6)
 *   - 	toupper(*which) == 'U': returns the upper triangular matrix:
 *		\phi_xx, \phi_xy, \phi_xz, \phi_yy, \phi_yz, \phi_zz (N = 6)
 *   - toupper(*which) == 'F': returns the full Hessian matrix:
 *		\phi_xx, \phi_xy, \phi_xz,
 *		\phi_yx, \phi_yy, \phi_yz,
 *		\phi_zx, \phi_zy, \phi_zz
 *   		(N = 9)
 *
 * Note:
 *   1. A code for testing this function is available in src/specht.c.
 *   2. To save some redundant ops, calls to this function with the same
 *	'index' value should be grouped together.
 */
{
    /* Note: data[] is recomputed only when J, index or dof->type is changed */
    static FLOAT *data = NULL;
    static FLOAT J0[(Dim+1)*(Dim+1)] =
			{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    static DOF_TYPE *type0 = NULL;
#if USE_OMP
#pragma omp threadprivate(data, J0, type0)
#endif  /* USE_OMP */

    DOF *dummy;
    DOF_TYPE *type, *grad_type;
    const FLOAT (*J)[Dim + 1] = (void *)phgGeomGetJacobian(dof->g, e);
    const FLOAT *p;
    int i, j, k, l;

    /* TODO: hp DOF */
    assert(dof->hp == NULL && dof->type != NULL);

    hessian_dof = dof;

    type = dof->type;
    grad_type = type->grad_type;	/* usually a DGn type */
    assert(grad_type != NULL && grad_type->dim == 1);
    assert(grad_type->np_vert + grad_type->np_edge + grad_type->np_face == 0);

    dummy = phgCalloc(1, sizeof(DOF));	/* dummy DOF */
    dummy->type = grad_type;
    dummy->g = dof->g;
    dummy->dim = Dim * type->dim * type->nbas;
    dummy->count_elem = grad_type->np_elem * dummy->dim;

    if (memcmp(J0, J, sizeof(J0)) || type0 != type) {
	/* compute DOFs (in data) for \grad_lambda\phi w.r.t. grad_type */
	if (type0 != type) {
	    phgFree(data);
	    /* Note: data stores temporary DOF w.r.t. grad_type in the order:
	     *		data[grad_type->nbas][type->nbas][type->dim][Dim] */
	    data = phgAlloc(type->nbas * type->dim * Dim
				* grad_type->nbas * sizeof(*data));
	    FreeAtExit(data);
	}
	memcpy(J0, J, sizeof(J0));
	type0 = type;
	grad_type->InitFunc(dummy, e, VOLUME, 0, NULL, hessian_func,
							NULL, data, NULL);
    }

    p = grad_type->BasGrads(dummy, e, 0, -1, lambda);

    for (l = 0; l < type->nbas * type->dim; l++) {
	for (i = 0; i < Dim; i++) {
	    int j0 = 0, j1 = Dim;
	    switch (toupper(*which)) {
		case 'D':	/* Diagonal */
		    j0 = i; j1 = i + 1;
		    break;
		case 'L':	/* Lower triangular matrix */
		    j0 = 0; j1 = i + 1;
		    break;
		case 'U':	/* Upper triangular matrix */
		    j0 = i; j1 = Dim;
		    break;
		case 'F':	/* Full matrix */
		case 'A':
		    j0 = 0; j1 = Dim;
		    break;
		default:
		    phgError(1, "(%s) invalid value for 'which': %s\n", which);
	    }
	    for (j = j0; j < j1; j++, values++) {
		*values = 0.0;
		for (k = 0; k < grad_type->nbas; k++)
		    *values += (p[(Dim + 1) * k + 0] * J[0][j] +
				p[(Dim + 1) * k + 1] * J[1][j] +
				p[(Dim + 1) * k + 2] * J[2][j] +
				p[(Dim + 1) * k + 3] * J[3][j])
			* data[(type->nbas * type->dim * k + l) * Dim + i];
	    }
	}
    }

    phgFree(dummy);
}

const char *_phg_dof_proj_name[] =
	{"DOF_PROJ_NONE", "DOF_PROJ_DOT", "DOF_PROJ_CROSS", NULL};

BOOLEAN
phgDofDirichletBC__(DOF *u, ELEMENT *e, int dof_index, DOF_USER_FUNC func,
		    DOF_USER_FUNC_LAMBDA func_lambda, DOF *func_dof,
		    FLOAT mat[], FLOAT rhs[], DOF_PROJ proj)
/* computes BC using either interpolation or L2-projection on the boundary.
 *
 * Note:
 *
 *  -1)	The BC is determined by the following arguments:
 *		func_dof, func_lambda, func
 *	(with the first non-nil pointer being effective)
 *
 *   0) dof_index = (element DOF index);
 *      bas_index = (element basis index) = (element DOF index / u->dim);
 *
 *   1) This function turns the Dirichlet BC at a given DOF into an equation.
 *	If u->dim == 1, the equation is given by:
 *		\sum_{0 <= j < u->type->nbas} mat[j] u[j] = *rhs
 *	where the coefficients mat[j] = <\phi_j, \phi_{bas_index}>, and
 *	the right hand side *rhs = <func, \phi_{bas_index}>, in which <.,.>
 *	denotes the inner product on an edge or on a face, depending on the
 *	geometric type of the DOF.
 *
 *	Note: if enforce_L2 == FALSE and u->type->is_nodal == TRUE, direct
 *	evaluation (interpolation) is used, i.e.:
 *		mat[bas_index] = 1.0, mat[j] = 0 for j != bas_index
 *		*rhs = the function value at DOF location.
 *	If enforce_L2 == TRUE, then L2 projection is enforced.
 *
 *   2) u->type->nbas coefficients are returned in mat[], which are independent
 *	of u->dim, and u->dim values are returned in rhs[] (if rhs != NULL).
 *
 *	The user program is resposible for duplicating the coefficients in
 *	mat[] according to u->dim when assembling the global matrix.
 *
 *   3) 'proj' must be either DOF_PROJ_NONE or DOF_PROJ_CROSS, the latter case
 *	requires u->type->dim == Dim and is a special case for the Dirichlet
 *	BC of the Maxwell equations: u \cros n = g \cross n
 *
 *   4) (func==NULL && func_lambda==NULL && func_dof==NULL) means func == 0
 */
{
    int bas_index = dof_index / u->dim;
    int sub_index = dof_index % u->dim;
    BOOLEAN enforce_L2 = FALSE;
    GRID *g = u->g;
    DOF_TYPE *type = u->type;
    int i, j, k, n0, n1, gind, bind, v0, v1, v2, v3;
    int N = DofNBas(u, e);
    const FLOAT *b, *qp, *qw, *p;
    FLOAT *q, lambda[Dim + 1], w0, x, y, z;
    FLOAT f[DofDim(u)], w[DofTypeDim(u)];
    QUAD *quad;
    BTYPE btype;
    GTYPE gtype;
    /* flags indicating which basis functions to use in the L2 projection */
    unsigned char flags[N];
    static BOOLEAN warned = FALSE;	/* warn unimplemented proj. only once */

    assert(func_dof == NULL || DofDim(func_dof) == DofDim(u));

    /* Warning: if 'func' calls phgDofEval with the same DOF_TYPE, then the
	values from BasFuncs are destroyed after calling 'func' */
    CheckD_order(u);
    assert(proj == DOF_PROJ_NONE || DofTypeDim(u) == Dim);

    btype = phgDofGetElementBasisInfo(u, e, bas_index, &gtype, &gind, &bind);
    if (u->DB_masks == NULL) {
	if (!(btype & u->DB_mask))
	    return FALSE;
    }
    else {
	if (!(btype & u->DB_masks[sub_index])) 
	    return FALSE;
    }

    if (mat == NULL) {
	assert(rhs == NULL);
	return TRUE;
    }

    if (type == NULL)
	type = u->hp->info->types[u->hp->elem_order[e->index]];

    if (type->base_type != NULL) {
	assert(type->dim == 1);	/* for DG type, TODO otherwise */
	type = type->base_type;
    }

    for (i = 0; i < N; i++)
	mat[i] = 0.0;

#if 0
    /* only evaluate BC when the submesh is the owner */
    if (rhs != NULL && !(btype & OWNER)) {
	for (i = 0; i < DofDim(u); i++)
	    rhs[i] = 0.;
	return TRUE;
    }
#endif

    if (!enforce_L2 && proj == DOF_PROJ_NONE && type->is_nodal) {
	/* Lagrange basis, return diagonal equations for BC */
	mat[bas_index] = 1.0;
	if (rhs == NULL)
	    return TRUE;
	if (func_dof != NULL) {
	    phgDofGetElementBasisInfo_(u, e, bas_index,
				       NULL, NULL, NULL, NULL, lambda);
	    phgDofEval(func_dof, e, lambda, rhs);
	}
	else if (func_lambda != NULL) {
	    phgDofGetElementBasisInfo_(u, e, bas_index,
				       NULL, NULL, NULL, NULL, lambda);
	    func_lambda(u, e, bas_index, lambda, rhs);
	}
	else if (func != NULL) {
	    p = phgDofGetElementCoordinates(u, e, bas_index * u->dim);
	    func(p[0], p[1], p[2], rhs);
	}
	else {
	    bzero(rhs, u->dim * type->dim * sizeof(*rhs));
	}
	return TRUE;
    }

    /* TODO: up to now with all basis functions, bdry DOFs can be
     * evaluated by calling InitFunc which only needs func values
     * at the boundary, moreover phgDofDirichletBC is called
     * in proper order. Thus, direct evaluation of bdry DOFs by
     * calling InitFunc is possible and is faster (when np_xxxx > 1,
     * the values can be cached for the next few calls) */

    /* Note: use different quadrature orders for different basis
     * functions may save some (marginal) computations */

#define FUNC(values)						\
	if (func_dof != NULL) {					\
	    phgDofEval(func_dof, e, lambda, values);		\
	}							\
	else if (func_lambda != NULL) {				\
	    func_lambda(u, e, bas_index, lambda, values);	\
	}							\
	else if (func != NULL) {				\
	    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);	\
	    func(x, y, z, values);				\
	}							\
	else {							\
	    bzero(values, u->dim * type->dim * sizeof(*rhs));	\
	}

    switch (gtype) {
	case VERTEX:
	    assert(proj == DOF_PROJ_NONE);
	    lambda[0] = lambda[1] = lambda[2] = lambda[3] = 0.0;
	    lambda[gind] = 1.0;
	    if (DofIsHP(u)) {
		assert(u->hp->info->flags & VERT_FLAG);
		n0 = 0;
		for (i = 0; i < gind; i++)
		    n0 += u->hp->info->types
				[u->hp->vert_order[e->verts[i]]]->np_vert;
		type = u->hp->info->types[u->hp->vert_order[e->verts[gind]]];
	    }
	    else {
		n0 = gind * type->np_vert;
	    }
	    n1 = n0 + type->np_vert;
	    b = type->BasFuncs(u, e, n0, n1, lambda);
	    if (type->dim == 1) {
		/* faster code for type->dim == 1 */
		w0 = b[bas_index - n0];
		for (i = n0; i < n1; i++)
		    mat[i] = w0 * *(b++);
		if (rhs == NULL)
		    break;
		FUNC(rhs)
		for (i = 0; i < u->dim; i++)
		    rhs[i] *= w0;
	    }
	    else {
		p = b + (bas_index - n0) * type->dim;
		for (i = n0; i < n1; i++) {
		    w0 = *(b++) * p[0];
		    for (j = 1; j < type->dim; j++)
			w0 += *(b++) * p[j];
		    mat[i] = w0;
		}
		if (rhs == NULL)
		    break;
		FUNC(q = f)
		for (i = 0; i < u->dim; i++) {
		    w0 = *(q++) * p[0];
		    for (j = 1; j < type->dim; j++)
			w0 += *(q++) * p[j];
		    *(rhs++) = w0;
		}
	    }
	    break;
	case EDGE:
	    /* L2 projection on the edge */
	    if (rhs != NULL)
		for (i = 0; i < u->dim; i++)
		    rhs[i] = 0.0;
	    GetEdgeVertices(e, gind, v0, v1);
	    memset(flags, 0, N);
	    if (!DofIsHP(u)) {
		if (type->np_vert > 0) {
		    /* basis functions on the vertices of the current edge */
		    memset(flags + v0 * type->np_vert, 1, type->np_vert);
		    memset(flags + v1 * type->np_vert, 1, type->np_vert);
		}
		n0 = type->np_vert * NVert + type->np_edge * gind;
	    }
	    else {
		n0 = 0;
		/* basis functions on the vertices of the current edge */
		for (k = 0;
		    (u->hp->info->flags & VERT_FLAG) && k < NVert; k++) {
		    type = u->hp->info->types[u->hp->vert_order[e->verts[k]]];
		    if (k == v0 || k == v1)
			memset(flags + n0, 1, type->np_vert);
		    n0 += type->np_vert;
		}
		assert(u->hp->info->flags & EDGE_FLAG);
		for (k = 0; k < gind; k++)
		    n0 += u->hp->info->types
				[u->hp->edge_order[e->edges[k]]]->np_edge;
		type = u->hp->info->types[u->hp->edge_order[e->edges[gind]]];
	    }
	    /* basis functions on the current edge */
	    memset(flags + n0, 1, type->np_edge);
	    quad = phgQuadGetQuad1D(2 * type->order);
	    lambda[0] = lambda[1] = lambda[2] = lambda[3] = 0.0;
	    qp = quad->points;
	    qw = quad->weights;
	    for (k = 0; k < quad->npoints; k++) {
		lambda[v0] = *(qp++);
		lambda[v1] = *(qp++);
		b = type->BasFuncs(u, e, 0, -1, lambda);
		w0 = *(qw++);
		if (type->dim == 1) {	/*-----------------------------------*/
		    /* faster code for type->dim == 1 */
		    w0 *= b[bas_index];
		    for (i = 0; i < N; i++) {
			if (flags[i] != 0) {
			    mat[i] += b[i] * w0;
			}
		    }
		    if (rhs != NULL) {
			FUNC(f)
			for (i = 0; i < u->dim; i++)
			    rhs[i] += f[i] * w0;
		    }
		}
		else if (proj == DOF_PROJ_CROSS) {	/*-------------------*/
		    /* Dirichlet BC for Maxwell equations: n x u = n x g */
		    /* Bdry eqn:  <t . u, t . phi_j> = <t . g, t . phg_j>,
		     * where t is the edge vector */
		    COORD *c0, *c1;
		    FLOAT t0, t1, t2;
		    c0 = g->verts + e->verts[v0];
		    c1 = g->verts + e->verts[v1];
		    t0 = (*c1)[0] - (*c0)[0];
		    t1 = (*c1)[1] - (*c0)[1];
		    t2 = (*c1)[2] - (*c0)[2];
		    q = (void *)(b + bas_index * Dim);
		    w0 *= t0 * q[0] + t1 * q[1] + t2 * q[2];
		    for (i = 0; i < N; i++, b += Dim) {
			if (flags[i] != 0) {
			    mat[i] += w0 * (t0 * b[0] + t1 * b[1] + t2 * b[2]);
			}
		    }
		    if (rhs != NULL) {
			FUNC(f)
			for (i = 0, q = f; i < u->dim; i++, q += Dim)
			    rhs[i] += w0 * (t0 * q[0] + t1 * q[1] + t2 * q[2]);
		    }
		}
		else {	/*---------------------------------------------------*/
		    if (proj != DOF_PROJ_NONE && g->rank == 0 && !warned) {
			phgWarning("%s(%s:%d): %s unimplemented, using %s.\n",
					__func__, __FILE__, __LINE__,
					_phg_dof_proj_name[proj],
					_phg_dof_proj_name[DOF_PROJ_NONE]);
			warned = TRUE;
		    }
		    p = b + bas_index * type->dim;
		    for (j = 0; j < type->dim; j++)
			w[j] = *(p++) * w0;
		    for (i = 0; i < N; i++) {
			if (flags[i] == 0) {
			    b += type->dim;
			}
			else {
			    w0 = mat[i];
			    for (j = 0; j < type->dim; j++)
				w0 += *(b++) * w[j];
			    mat[i] = w0;
			}
		    }	/* for i */
		    if (rhs != NULL) {
			FUNC(q = f)
			for (i = 0; i < u->dim; i++) {
			    w0 = rhs[i];
			    for (j = 0; j < type->dim; j++)
				w0 += *(q++) * w[j];
			    rhs[i] = w0;
			}
		    }
		}	/*---------------------------------------------------*/
	    }	/* for k */
#if 0
	    /* multiply by edge length */
	    {
		COORD *c0, *c1;
		c0 = g->verts + e->verts[v0];
		c1 = g->verts + e->verts[v1];
		x = (*c1)[0] - (*c0)[0];
		y = (*c1)[1] - (*c0)[1];
		z = (*c1)[2] - (*c0)[2];
		w = Sqrt(x * x + y * y + z * z);
	    }
	    for (i = 0; i < N; i++)
		mat[i] *= w;
	    if (rhs == NULL)
		break;
	    for (i = 0; i < u->dim; i++)
		rhs[i] *= w;
#endif
	    break;
	case FACE:
	    /* L2 projection on the face */
	    if (rhs != NULL)
		for (i = 0; i < u->dim; i++)
		    rhs[i] = 0.0;
	    GetFaceVertices(e, gind, v0, v1, v2, v3);
	    memset(flags, 0, N);
	    if (!DofIsHP(u)) {
		/* basis functions on the vertices of the current face */
		if (type->np_vert > 0) {
		    memset(flags + v0 * type->np_vert, 1, type->np_vert);
		    memset(flags + v1 * type->np_vert, 1, type->np_vert);
		    memset(flags + v2 * type->np_vert, 1, type->np_vert);
		}
		n0 = type->np_vert * NVert;
		if (type->np_edge > 0) {
		    i = GetEdgeNo(v0, v1);
		    memset(flags + n0 + i * type->np_edge, 1, type->np_edge);
		    i = GetEdgeNo(v0, v2);
		    memset(flags + n0 + i * type->np_edge, 1, type->np_edge);
		    i = GetEdgeNo(v1, v2);
		    memset(flags + n0 + i * type->np_edge, 1, type->np_edge);
	        }
		/* basis functions on the current face */
		n0 += type->np_edge * NEdge + type->np_face * gind;
	    }
	    else {
		int e0, e1, e2;
		n0 = 0;
		for (k = 0; (u->hp->info->flags & VERT_FLAG) && k < NVert;
		     k++) {
		    type = u->hp->info->types[u->hp->vert_order[e->verts[k]]];
		    if (k == v0 || k == v1 || k == v2)
			memset(flags + n0, 1, type->np_vert);
		    n0 += type->np_vert;
		}
		e0 = GetEdgeNo(v0, v1);
		e1 = GetEdgeNo(v0, v2);
		e2 = GetEdgeNo(v1, v2);
		for (k = 0; (u->hp->info->flags & EDGE_FLAG) && k < NEdge;
		     k++) {
		    type = u->hp->info->types[u->hp->edge_order[e->edges[k]]];
		    if (k == e0 || k == e1 || k == e2)
			memset(flags + n0, 1, type->np_edge);
		    n0 += type->np_edge;
		}
		/* basis functions on the current face */
		assert(u->hp->info->flags & FACE_FLAG);
		for (k = 0; k < gind; k++)
		    n0 += u->hp->info->types
				[u->hp->face_order[e->faces[k]]]->np_face;
		type = u->hp->info->types[u->hp->face_order[e->faces[gind]]];
	    }
	    memset(flags + n0, 1, type->np_face);
	    quad = phgQuadGetQuad2D(2 * type->order);
	    lambda[gind] = 0.0;
	    qp = quad->points;
	    qw = quad->weights;
	    for (k = 0; k < quad->npoints; k++) {
		lambda[v0] = *(qp++);
		lambda[v1] = *(qp++);
		lambda[v2] = *(qp++);
		b = type->BasFuncs(u, e, 0, -1, lambda);
		w0 = *(qw++);
		if (type->dim == 1) {
		    /* faster code for type->dim == 1 */
		    w0 *= b[bas_index];
		    for (i = 0; i < N; i++) {
			if (flags[i] != 0) {
			    mat[i] += b[i] * w0;
			}
		    }
		    if (rhs != NULL) {
			FUNC(f)
			for (i = 0; i < u->dim; i++)
			    rhs[i] += f[i] * w0;
		    }
		}
		else {
		    FLOAT n[Dim];
		    phgGeomGetFaceOutNormal(g, e, gind, n);
		    int m = type->dim;
		    if (proj == DOF_PROJ_CROSS) {
			/* Dirichlet BC for Maxwell equations: n x u = n x g */
			/* Bdry eqn:  <u x n, phi_j x n> = <g x n, phg_j x n> */

			/* compute basis x n */
			/* Warning: overriding static buffer in BasFuncs */
			assert(m == Dim);
			for (i = 0, q = (void *)b; i < N; i++, q += Dim) {
			    x = q[1] * n[2] - q[2] * n[1];
			    y = q[2] * n[0] - q[0] * n[2];
			    z = q[0] * n[1] - q[1] * n[0];
			    q[0] = x;
			    q[1] = y;
			    q[2] = z;
			}
		    }
		    else if (proj == DOF_PROJ_DOT) {
			/* compute basis . n */
			/* Warning: overriding static buffer in BasFuncs */
			FLOAT *r;
			assert(m == Dim);
			m /= Dim;
			for (i = 0, r = q = (void *)b; i < N; i++, q += Dim)
			    *(r++) = q[0] * n[0] + q[1] * n[1] + q[2] * n[2];
		    }
		    else if (proj != DOF_PROJ_NONE && g->rank == 0 && !warned) {
			phgWarning("%s(%s:%d): %s unimplemented, using %s.\n",
					__func__, __FILE__, __LINE__,
					_phg_dof_proj_name[proj],
					_phg_dof_proj_name[DOF_PROJ_NONE]);
			warned = TRUE;
		    }
		    p = b + bas_index * m;
		    for (j = 0; j < m; j++)
			w[j] = *(p++) * w0;
		    for (i = 0; i < N; i++) {
			if (flags[i] == 0) {
			    b += m;
			}
			else {
			    w0 = mat[i];
			    for (j = 0; j < m; j++)
				w0 += *(b++) * w[j];
			    mat[i] = w0;
			}
		    }
		    if (rhs != NULL) {
			int m = type->dim;
			FUNC(q = f)
			if (proj == DOF_PROJ_CROSS) {
			    /* compute g x n */
			    for (i = 0; i < u->dim; i++, q += Dim) {
				x = q[1] * n[2] - q[2] * n[1];
				y = q[2] * n[0] - q[0] * n[2];
				z = q[0] * n[1] - q[1] * n[0];
				q[0] = x;
				q[1] = y;
				q[2] = z;
			    }
			    q = f;
			}
			else if (proj == DOF_PROJ_DOT) {
			    FLOAT *r = q;
			    m /= Dim;
			    /* compute g . n */
			    for (i = 0; i < u->dim; i++, q += Dim)
				*(r++) = q[0]*n[0] + q[1]*n[1] + q[2]*n[2];
			    q = f;
			}
			else if (proj != DOF_PROJ_NONE && g->rank == 0
					&& !warned) {
			    phgWarning("%s(%s:%d): %s unimplemented.\n",
					__func__, __FILE__, __LINE__,
					_phg_dof_proj_name[proj]);
			    warned = TRUE;
			}
			for (i = 0; i < u->dim; i++) {
			    w0 = rhs[i];
			    for (j = 0; j < m; j++)
				w0 += *(q++) * w[j];
			    rhs[i] = w0;
			}
		    }
		}
	    }
#if 0
	    /* multiply by face area */
	    w = phgGeomGetFaceArea(g, e, gind);
	    for (i = 0; i < N; i++)
		mat[i] *= w;
	    if (rhs == NULL)
		break;
	    for (i = 0; i < u->dim; i++)
		rhs[i] *= w;
#endif
	    break;
	case VOLUME:
	    assert(proj == DOF_PROJ_NONE);
	    phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__);
	    break;
    }

#undef FUNC

    return TRUE;	/* to make gcc happy */
}

/*--------------------------- BLAS-like functiones ----------------------*/

static FLOAT *bas;
static size_t bas_count;
#if USE_OMP
# pragma omp threadprivate(bas)
#endif	/* USE_OMP */

static void
func(DOF *dof, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
    Unused(bno);
    CheckD_order(dof);
    phgDofEval_(cache_dof, e, lambda, bas, 0, values, NULL);
    bas += bas_count;
}

DOF *
phgDofCopy_(DOF *src, DOF **dest_ptr, DOF_TYPE *newtype,
	    const char *name, const char *srcfile, int srcline)
/* returns a new DOF whose type is 'newtype', and whose values are evaluated
 * using 'src' (a copy of src). */
{
    GRID *g = src->g;
    ELEMENT *e;
    FLOAT w = 1.0, *basvalues = NULL, *a, *d, *buffer = NULL;
    int i, k, dim, dim_new, count = 0, nvalues, nd;
    INT n;
    size_t size;
    DOF *wgts = NULL;
    DOF *dest = (dest_ptr == NULL ? NULL : *dest_ptr);
    BYTE *flags0, *flags;
#if USE_OMP
    static FLOAT *bas0 = NULL;
#pragma omp threadprivate(bas0)
    omp_lock_t *locks0, *locks;
#endif	/* USE_OMP */
    char *auto_name;
    DOF_TYPE *oldtype;
    BOOLEAN basflag = FALSE;

    MagicCheck(DOF, dest)

    if (dest != NULL && newtype != NULL && dest->type != newtype) {
	phgDofFree(&dest);
	dest = NULL;
    }

    dim = DofDim(src);
    /* the name of dest */
    if ((auto_name = (void *)name) == NULL) {
#if 0
	auto_name = phgAlloc(strlen(src->name) + 8 + 1);
	sprintf(auto_name, "copy of %s", src->name);
#else
	auto_name = strdup(src->name);
#endif
    }

    if (dest == NULL) {
	if (newtype == NULL)
	    newtype = src->type;
	dim_new = (newtype == NULL ? DofTypeDim(src) : newtype->dim);
	assert(dim % dim_new == 0);
	dest = phgDofNew_(g, newtype, newtype == NULL ? src->hp : NULL,
			  dim / dim_new, auto_name, DofNoAction,
			  srcfile, srcline);
	if (dest_ptr != NULL)
	    *dest_ptr = dest;
    }
    else {
	assert(dim == DofDim(dest));
	phgFree(dest->name); dest->name = strdup(auto_name);
	dest->srcfile = srcfile;
	dest->srcline = srcline;
	phgFree(dest->cache_func); dest->cache_func = NULL;
	newtype = dest->type;
	if (!SpecialDofType(newtype))
	    memset(dest->data, 0, DofDataCount(dest) * sizeof(*dest->data));
    }
    if (auto_name != name)
	phgFree(auto_name);

    phgDofClearCache(NULL, dest, NULL, NULL, FALSE);

    dest->DB_mask = src->DB_mask;
    if (src->DB_masks != NULL) {
	if (dest->DB_masks == NULL)
	    dest->DB_masks = phgAlloc(dest->dim * sizeof(*dest->DB_masks));
	memcpy(dest->DB_masks, src->DB_masks, dest->dim * sizeof(*dest->DB_masks));
    }
    dest->invariant = src->invariant;

    if (SpecialDofType(newtype)) {
	assert(newtype == src->type);
	dest->userfunc = src->userfunc;
	dest->userfunc_lambda = src->userfunc_lambda;
	if (newtype == DOF_CONSTANT)
	    memcpy(dest->data, src->data, dim * sizeof(*dest->data));
	return dest;
    }

    if ((newtype != NULL && newtype == src->type) ||
	(newtype == NULL && src->hp == dest->hp)) {
	/* simply duplicate the data */
	size_t size = DofDataCount(dest);
	if (size > 0)
	    memcpy(dest->data, src->data, sizeof(FLOAT) * size);
	dest->userfunc = src->userfunc;
	dest->userfunc_lambda = src->userfunc_lambda;
	return dest;
    }

    if (src->type == DOF_ANALYTIC) {
	CheckD_order(dest);
	if (src->userfunc_lambda != NULL)
	    phgDofSetDataByLambdaFunction(dest, src->userfunc_lambda);
	else
	    phgDofSetDataByFunction(dest, src->userfunc);
	return dest;
    }

    if (newtype != NULL && newtype->BasFuncs == NULL)
	phgError(1, "phgDofCopy: basis funcs for DOF type \"%s\" undefined.\n",
		 newtype->name);

    dest->userfunc = src->userfunc;
    dest->userfunc_lambda = src->userfunc_lambda;

    oldtype = src->type;
    if (oldtype == NULL)
	oldtype = src->hp->info->types[src->hp->info->min_order];

    if (!SpecialDofType(oldtype) && newtype != NULL && newtype->is_nodal
	&& !DofIsHP(src)) {
	basflag = TRUE;
	count = oldtype->nbas * oldtype->dim;
	basvalues = phgAlloc(newtype->nbas * count * sizeof(*basvalues));

	if (oldtype->invariant == TRUE)
	    get_bas_funcs(src, dest, src->g->roots, basvalues);
    }

    if (newtype == NULL)
	newtype = dest->hp->info->types[dest->hp->max_order];

    size = (newtype->np_vert > 0 ? g->nvert : 0) +
	   (newtype->np_edge > 0 ? g->nedge : 0) +
	   (newtype->np_face > 0 ? g->nface : 0);
    flags0 = phgCalloc(size, sizeof(*flags0));

    if (!SpecialDofType(oldtype) && oldtype->continuity < 0) {
	static DOF_TYPE DOF_WGTS = {DofReserved,
	    "Weights", NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	    phgDofInitFuncPoint, NULL, NULL, NULL, FE_None,
	    FALSE, FALSE, FALSE, -1, 0, 0, 0, -1, 1, 0, 0, 0, 0
	};
	DOF_WGTS.np_vert = (newtype->np_vert > 0) ? 1 : 0;
	DOF_WGTS.np_edge = (newtype->np_edge > 0) ? 1 : 0;
	DOF_WGTS.np_face = (newtype->np_face > 0) ? 1 : 0;
	DOF_WGTS.nbas = DOF_WGTS.np_vert * NVert + DOF_WGTS.np_edge * NEdge +
	    DOF_WGTS.np_face * NFace;
	if (DOF_WGTS.nbas > 0) {
	    /* Other cases will be implemented later when really needed */
	    wgts = phgDofNew(g, &DOF_WGTS, 1, "weights", DofNoAction);
	    phgDofSetDataByValue(wgts, 0.0);
	    phgDofSetDataByValue(dest, 0.0);
	    /* allocate buffer for storing weighted vertex/edge/face data */
	    if ((n = DofDataCount(dest) - DofElementDataCount(dest)) > 0)
		buffer = phgCalloc(n, sizeof(*buffer));
	}
    }

    nvalues = dest->dim;
    cache_dof = src;
    bas_count = count;

#if USE_OMP
    locks0 = phgAlloc(size * sizeof(*locks0));
    for (n = 0; n < size; n++)
	omp_init_lock(&locks0[n]);

    /* initialize bas0 to NULL */
# pragma omp parallel for schedule(static)
    for (i = 0; i < phgMaxThreads; i++)
	bas0 = NULL;

    /* parallel loop */
# pragma omp parallel for private(e, w, nd, i, k, n, a, d, flags, locks)
#endif  /* USE_OMP */
    ForAllElementsBegin_(g, e, FOR_ALL) {
	if (wgts != NULL) {
#if 0
	    /* use element volume as weight */
	    w = phgGeomGetVolume(g, e);
#else
	    /* use 1 as weight */
	    w = 1.0;
#endif
	}

	if (basflag && oldtype->invariant == FALSE) {
#if USE_OMP
	    if (bas0 == NULL)
		bas0 = phgAlloc(newtype->nbas * count * sizeof(*bas0));
	    get_bas_funcs(src, dest, e, bas0);
	    bas = bas0;
#else	/* USE_OMP */
	    get_bas_funcs(src, dest, e, basvalues);
	    bas = basvalues;
#endif	/* USE_OMP */
	}
	else {
	    bas = basvalues;
	}
	flags = flags0;
#if USE_OMP
	locks = locks0;
#endif  /* USE_OMP */

	if (DofIsHP(dest))
	    newtype = dest->hp->info->types[dest->hp->elem_order[e->index]];

	if (newtype->np_vert > 0) {
	    nd = nvalues * newtype->np_vert;
	    for (k = 0; k < NVert; k++) {
		if (flags[n = e->verts[k]] && wgts == NULL) {
		    /* note: count==0 and bas==NULL for variable order DOF */
		    bas += count * newtype->np_vert;
		    continue;
		}
#if USE_OMP
		omp_set_lock(&locks[n]);
#endif  /* USE_OMP */
		flags[n] = TRUE;
		a = DofVertexData(dest, n);
		newtype->InitFunc(dest, e, VERTEX, k, NULL, func, NULL, a,
					NULL);
		if (wgts != NULL) {
		    d = buffer + (a - DofData(dest));
		    if (DofIsHP(dest))
			nd = dest->dim * (dest->hp->vert_index[n + 1] -
					  dest->hp->vert_index[n]);
		    for (i = 0; i < nd; i++)
			*(d++) += *(a++) * w;
		    *DofVertexData(wgts, n) += w;
		}
#if USE_OMP
		omp_unset_lock(&locks[n]);
#endif  /* USE_OMP */
	    }
	    flags += g->nvert;
#if USE_OMP
	    locks += g->nvert;
#endif  /* USE_OMP */
	}

	if (newtype->np_edge > 0) {
	    nd = nvalues * newtype->np_edge;
	    for (k = 0; k < NEdge; k++) {
		if (flags[n = e->edges[k]] && wgts == NULL) {
		    /* note: count==0 and bas==NULL for variable order DOF */
		    bas += count * newtype->np_edge;
		    continue;
		}
#if USE_OMP
		omp_set_lock(&locks[n]);
#endif  /* USE_OMP */
		flags[n] = TRUE;
		a = DofEdgeData(dest, n);
		newtype->InitFunc(dest, e, EDGE, k, NULL, func, NULL, a, NULL);
		if (wgts != NULL) {
		    d = buffer + (a - DofData(dest));
		    if (DofIsHP(dest))
			nd = dest->dim * (dest->hp->edge_index[n + 1] -
					  dest->hp->edge_index[n]);
		    for (i = 0; i < nd; i++)
			*(d++) += *(a++) * w;
		    *DofEdgeData(wgts, n) += w;
		}
#if USE_OMP
		omp_unset_lock(&locks[n]);
#endif  /* USE_OMP */
	    }
	    flags += g->nedge;
#if USE_OMP
	    locks += g->nedge;
#endif  /* USE_OMP */
	}

	if (newtype->np_face > 0) {
	    nd = nvalues * newtype->np_face;
	    for (k = 0; k < NFace; k++) {
		if (flags[n = e->faces[k]] && wgts == NULL) {
		    /* note: count==0 and bas==NULL for variable order DOF */
		    bas += count * newtype->np_face;
		    continue;
		}
#if USE_OMP
		omp_set_lock(&locks[n]);
#endif  /* USE_OMP */
		flags[n] = TRUE;
		a = DofFaceData(dest, n);
		newtype->InitFunc(dest, e, FACE, k, NULL, func, NULL, a, NULL);
		if (wgts != NULL) {
		    d = buffer + (a - DofData(dest));
		    if (DofIsHP(dest))
			nd = dest->dim * (dest->hp->face_index[n + 1] -
						 dest->hp->face_index[n]);
		    for (i = 0; i < nd; i++)
			*(d++) += *(a++) * w;
		    *DofFaceData(wgts, n) += w;
		}
#if USE_OMP
		omp_unset_lock(&locks[n]);
#endif  /* USE_OMP */
	    }
	}

	if (newtype->np_elem > 0) {
	    a = DofElementData(dest, e->index);
	    newtype->InitFunc(dest, e, VOLUME, 0, NULL, func, NULL, a, NULL);
	}

    } ForAllElementsEnd

#if USE_OMP
    for (n = 0; n < size; n++)
	omp_init_lock(&locks0[n]);
    phgFree(locks0);
# pragma omp parallel for schedule(static)
    for (i = 0; i < phgMaxThreads; i++) {
	phgFree(bas0);
	bas0 = NULL;
    }
#endif	/* USE_OMP */

    phgFree(basvalues);
    phgFree(flags0);

    if (wgts == NULL)
	return dest;

    if ((n = DofDataCount(dest) - DofElementDataCount(dest)) > 0) {
	memcpy(DofData(dest), buffer, n * sizeof(*buffer));
	phgFree(buffer);
    }

    if (DofIsHP(dest))
	newtype = dest->hp->info->types[dest->hp->max_order];

#if USE_MPI
    /* FIXME: directly interacting with the vector assembly code in solver.c
       is more efficient */
    if (g->nprocs > 1) {
	SOLVER *solver = phgSolverCreate(SOLVER_DIAGONAL, dest, NULL);
	INT K = 0, I;
	int j;

	if (newtype->np_vert > 0) {
#if USE_OMP	/* seems slowers with parallel for!!! */
# pragma omp parallel for private(nd, a, d, K, j, I)
#endif	/* USE_OMP */
	    for (n = 0; n < g->nvert; n++) {
		if (g->types_vert[n] == UNREFERENCED)
		    continue;
		nd = DofDataCountOnVertex(dest, n);
		a = DofVertexData(dest, n);
		d = DofVertexData(wgts, n);
		K = (INT)(a - DofData(dest));
		for (j = 0; j < nd; j++, K++, a++) {
		    I = phgSolverMapD2L(solver, 0, K);
		    assert(I >= 0 && I < solver->mat->rmap->localsize);
		    phgSolverAddRHSEntry(solver, I, *a);
		    phgSolverAddMatrixEntry(solver, I, I, *d);
		}
	    }
	}

	if (newtype->np_edge > 0) {
#if USE_OMP
# pragma omp parallel for private(nd, a, d, K, j, I)
#endif	/* USE_OMP */
	    for (n = 0; n < g->nedge; n++) {
		if (g->types_edge[n] == UNREFERENCED)
		    continue;
		nd = DofDataCountOnEdge(dest, n);
		a = DofEdgeData(dest, n);
		d = DofEdgeData(wgts, n);
		K = (INT)(a - DofData(dest));
		for (j = 0; j < nd; j++, K++, a++) {
		    I = phgSolverMapD2L(solver, 0, K);
		    assert(I >= 0 && I < solver->mat->rmap->localsize);
		    phgSolverAddRHSEntry(solver, I, *a);
		    phgSolverAddMatrixEntry(solver, I, I, *d);
		}
	    }
	}

	if (newtype->np_face > 0) {
#if USE_OMP
# pragma omp parallel for private(nd, a, d, K, j, I)
#endif	/* USE_OMP */
	    for (n = 0; n < g->nface; n++) {
		if (g->types_face[n] == UNREFERENCED)
		    continue;
		nd = DofDataCountOnFace(dest, n);
		a = DofFaceData(dest, n);
		d = DofFaceData(wgts, n);
		K = (INT)(a - DofData(dest));
		for (j = 0; j < nd; j++, K++, a++) {
		    I = phgSolverMapD2L(solver, 0, K);
		    assert(I >= 0 && I < solver->mat->rmap->localsize);
		    phgSolverAddRHSEntry(solver, I, *a);
		    phgSolverAddMatrixEntry(solver, I, I, *d);
		}
	    }
	}

	if (newtype->np_elem > 0) {
#if USE_OMP
# pragma omp parallel for private(nd, a, d, K, j, I)
#endif	/* USE_OMP */
	    for (n = 0; n < g->nelem; n++) {
		if (g->types_elem[n] == UNREFERENCED)
		    continue;
		nd = DofDataCountOnElement(dest, n);
		a = DofElementData(dest, n);
		d = DofElementData(wgts, n);
		K = (INT)(a - DofData(dest));
		for (j = 0; j < nd; j++, K++, a++) {
		    I = phgSolverMapD2L(solver, 0, K);
		    assert(I >= 0 && I < solver->mat->rmap->localsize);
		    phgSolverAddRHSEntry(solver, I, *a);
		    phgSolverAddMatrixEntry(solver, I, I, 1.0);
		}
	    }
	}

	phgSolverSolve(solver, TRUE, dest, NULL);
	phgSolverDestroy(&solver);

	phgDofFree(&wgts);

	return dest;
    }
#endif

    a = DofData(dest);
    d = DofData(wgts);

    if (newtype->np_vert > 0) {
	k = nvalues * newtype->np_vert;
	for (n = 0; n < g->nvert; n++) {
	    if (DofIsHP(dest))
		k = dest->dim * (dest->hp->vert_index[n + 1] -
					dest->hp->vert_index[n]);
	    if ((w = *(d++)) == 0.0) {
		a += k;
	    }
	    else if (k == 1) {
		*(a++) /= w;
	    }
	    else {
		w = 1.0 / w;
		for (i = 0; i < k; i++)
		    *(a++) *= w;
	    }
	}
    }

    if (newtype->np_edge > 0) {
	k = nvalues * newtype->np_edge;
	for (n = 0; n < g->nedge; n++) {
	    if (DofIsHP(dest))
		k = dest->dim * (dest->hp->edge_index[n + 1] -
					dest->hp->edge_index[n]);
	    if ((w = *(d++)) == 0.0) {
		a += k;
	    }
	    else if (k == 1) {
		*(a++) /= w;
	    }
	    else {
		w = 1.0 / w;
		for (i = 0; i < k; i++)
		    *(a++) *= w;
	    }
	}
    }

    if (newtype->np_face > 0) {
	k = nvalues * newtype->np_face;
	for (n = 0; n < g->nface; n++) {
	    if (DofIsHP(dest))
		k = dest->dim * (dest->hp->face_index[n + 1] -
					dest->hp->face_index[n]);
	    if ((w = *(d++)) == 0.0) {
		a += k;
	    }
	    else if (k == 1) {
		*(a++) /= w;
	    }
	    else {
		w = 1.0 / w;
		for (i = 0; i < k; i++)
		    *(a++) *= w;
	    }
	}
    }

    phgDofFree(&wgts);

    return dest;
}

FLOAT
phgDofDotL2Vec(DOF *x, DOF *y)
/* vector L2 norm */
{
    GRID *g = x->g;
    INT i;
    FLOAT *dx = x->data, *dy = y->data;
    FLOAT a;
    int ii, n;

    assert(x->type == y->type && x->hp == y->hp);	/* this is handy */
    assert(x->dim == y->dim);
    assert(!SpecialDofType(x->type));

    a = 0.0;

    if (DofVertexDataCount(x) > 0) {
	n = x->count_vert;
	for (i = 0; i < g->nvert; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->vert_index[i + 1] - x->hp->vert_index[i]);
	    if (!(g->types_vert[i] & OWNER)) {
		dx += n;
		dy += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++)
		a += *(dx++) * *(dy++);
	}
    }

    if (DofEdgeDataCount(x) > 0) {
	n = x->count_edge;
	for (i = 0; i < g->nedge; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->edge_index[i + 1] - x->hp->edge_index[i]);
	    if (!(g->types_edge[i] & OWNER)) {
		dx += n;
		dy += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++)
		a += *(dx++) * *(dy++);
	}
    }

    if (DofFaceDataCount(x) > 0) {
	n = x->count_face;
	for (i = 0; i < g->nface; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->face_index[i + 1] - x->hp->face_index[i]);
	    if (!(g->types_face[i] & OWNER)) {
		dx += n;
		dy += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++)
		a += *(dx++) * *(dy++);
	}
    }

    if (DofElementDataCount(x) > 0) {
	n = x->count_elem;
	for (i = 0; i < g->nelem; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->elem_index[i + 1] - x->hp->elem_index[i]);
	    if (!(g->types_elem[i] & OWNER)) {
		dx += n;
		dy += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++)
		a += *(dx++) * *(dy++);
	}
    }

#if USE_MPI
    if (g->nprocs > 1) {
	FLOAT b;
	MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_SUM, x->g->comm);
	a = b;
    }
#endif	/* USE_MPI */

    return a;
}

FLOAT
phgDofNormLpVec(DOF *x, FLOAT p)
/* vector Lp norm */
{
    GRID *g = x->g;
    INT i;
    FLOAT *d = x->data;
    FLOAT a;
    int ii, n;

    assert(p >= 0.);
    assert(!SpecialDofType(x->type));

    a = 0.0;

    if (DofVertexDataCount(x) > 0) {
	n = x->count_vert;
	for (i = 0; i < g->nvert; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->vert_index[i + 1] - x->hp->vert_index[i]);
	    if (!(g->types_vert[i] & OWNER)) {
		d += n;
		continue;
	    }
	    if (p == 0.)
		for (ii = 0; ii < n; ii++, d++)
		    a += (*d == 0. ? 0. : 1.);
	    else
		for (ii = 0; ii < n; ii++, d++)
		    a += Pow(Fabs(*d), p);
	}
    }

    if (DofEdgeDataCount(x) > 0) {
	n = x->count_edge;
	for (i = 0; i < g->nedge; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->edge_index[i + 1] - x->hp->edge_index[i]);
	    if (!(g->types_edge[i] & OWNER)) {
		d += n;
		continue;
	    }
	    if (p == 0.)
		for (ii = 0; ii < n; ii++, d++)
		    a += (*d == 0. ? 0. : 1.);
	    else
		for (ii = 0; ii < n; ii++, d++)
		    a += Pow(Fabs(*d), p);
	}
    }

    if (DofFaceDataCount(x) > 0) {
	n = x->count_face;
	for (i = 0; i < g->nface; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->face_index[i + 1] - x->hp->face_index[i]);
	    if (!(g->types_face[i] & OWNER)) {
		d += n;
		continue;
	    }
	    if (p == 0.)
		for (ii = 0; ii < n; ii++, d++)
		    a += (*d == 0. ? 0. : 1.);
	    else
		for (ii = 0; ii < n; ii++, d++)
		    a += Pow(Fabs(*d), p);
	}
    }

    if (DofElementDataCount(x) > 0) {
	n = x->count_elem;
	for (i = 0; i < g->nelem; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->elem_index[i + 1] - x->hp->elem_index[i]);
	    if (!(g->types_elem[i] & OWNER)) {
		d += n;
		continue;
	    }
	    if (p == 0.)
		for (ii = 0; ii < n; ii++, d++)
		    a += (*d == 0. ? 0. : 1.);
	    else
		for (ii = 0; ii < n; ii++, d++)
		    a += Pow(Fabs(*d), p);
	}
    }

#if USE_MPI
    if (g->nprocs > 1) {
	FLOAT b;
	MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_SUM, x->g->comm);
	a = b;
    }
#endif	/* USE_MPI */

    return Pow(a, p == 0.0 ? 1.0 : 1.0 / p);
}

FLOAT
phgDofNormL1Vec(DOF *x)
/* vector L1 norm */
{
    GRID *g = x->g;
    INT i;
    FLOAT *d = x->data;
    FLOAT a;
    int ii, n;
#if USE_MPI
    FLOAT b;
#endif

    a = 0.0;

    if (DofVertexDataCount(x) > 0) {
	n = x->count_vert;
	for (i = 0; i < g->nvert; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->vert_index[i + 1] - x->hp->vert_index[i]);
	    if (!(g->types_vert[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		a += Fabs(*(d++));
	    }
	}
    }

    if (DofEdgeDataCount(x) > 0) {
	n = x->count_edge;
	for (i = 0; i < g->nedge; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->edge_index[i + 1] - x->hp->edge_index[i]);
	    if (!(g->types_edge[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		a += Fabs(*(d++));
	    }
	}
    }

    if (DofFaceDataCount(x) > 0) {
	n = x->count_face;
	for (i = 0; i < g->nface; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->face_index[i + 1] - x->hp->face_index[i]);
	    if (!(g->types_face[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		a += Fabs(*(d++));
	    }
	}
    }

    if (DofElementDataCount(x) > 0) {
	n = x->count_elem;
	for (i = 0; i < g->nelem; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->elem_index[i + 1] - x->hp->elem_index[i]);
	    if (!(g->types_elem[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		a += Fabs(*(d++));
	    }
	}
    }

#if USE_MPI
    MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_SUM, x->g->comm);
    return b;
#else
    return a;
#endif
}

FLOAT
phgDofNormInftyVec(DOF *x)
{
    GRID *g = x->g;
    INT i;
    FLOAT *d = x->data;
    FLOAT a, b;
    int ii, n;

    a = 0.0;

    if (DofVertexDataCount(x) > 0) {
	n = x->count_vert;
	for (i = 0; i < g->nvert; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->vert_index[i + 1] - x->hp->vert_index[i]);
	    if (!(g->types_vert[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a < (b = Fabs(*(d++))))
		    a = b;
	    }
	}
    }

    if (DofEdgeDataCount(x) > 0) {
	n = x->count_edge;
	for (i = 0; i < g->nedge; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->edge_index[i + 1] - x->hp->edge_index[i]);
	    if (!(g->types_edge[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a < (b = Fabs(*(d++))))
		    a = b;
	    }
	}
    }

    if (DofFaceDataCount(x) > 0) {
	n = x->count_face;
	for (i = 0; i < g->nface; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->face_index[i + 1] - x->hp->face_index[i]);
	    if (!(g->types_face[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a < (b = Fabs(*(d++))))
		    a = b;
	    }
	}
    }

    if (DofElementDataCount(x) > 0) {
	n = x->count_elem;
	for (i = 0; i < g->nelem; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->elem_index[i + 1] - x->hp->elem_index[i]);
	    if (!(g->types_elem[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a < (b = Fabs(*(d++))))
		    a = b;
	    }
	}
    }

#if USE_MPI
    MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_MAX, x->g->comm);
    return b;
#else
    return a;
#endif
}

FLOAT
phgDofMinValVec(DOF *x)
{
    GRID *g = x->g;
    INT i;
    FLOAT *d = x->data;
    FLOAT a, b;
    int ii, n;

    a = FLOAT_MAX;

    if (DofVertexDataCount(x) > 0) {
	n = x->count_vert;
	for (i = 0; i < g->nvert; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->vert_index[i + 1] -
			      x->hp->vert_index[i]);
	    if (!(g->types_vert[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a > (b = *(d++)))
		    a = b;
	    }
	}
    }

    if (DofEdgeDataCount(x) > 0) {
	n = x->count_edge;
	for (i = 0; i < g->nedge; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->edge_index[i + 1] -
			      x->hp->edge_index[i]);
	    if (!(g->types_edge[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a > (b = *(d++)))
		    a = b;
	    }
	}
    }

    if (DofFaceDataCount(x) > 0) {
	n = x->count_face;
	for (i = 0; i < g->nface; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->face_index[i + 1] -
			      x->hp->face_index[i]);
	    if (!(g->types_face[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a > (b = *(d++)))
		    a = b;
	    }
	}
    }

    if (DofElementDataCount(x) > 0) {
	n = x->count_elem;
	for (i = 0; i < g->nelem; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->elem_index[i + 1] -
			      x->hp->elem_index[i]);
	    if (!(g->types_elem[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a > (b = *(d++)))
		    a = b;
	    }
	}
    }

#if USE_MPI
    MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, MPI_MIN, x->g->comm);
    return b;
#else	/* USE_MPI */
    return a;
#endif	/* USE_MPI */
}

FLOAT
phgDofMaxValVec(DOF *x)
{
    GRID *g = x->g;
    INT i;
    FLOAT *d = x->data;
    FLOAT a, b;
    int ii, n;

    a = -FLOAT_MAX;

    if (DofVertexDataCount(x) > 0) {
	n = x->count_vert;
	for (i = 0; i < g->nvert; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->vert_index[i + 1] -
			      x->hp->vert_index[i]);
	    if (!(g->types_vert[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a < (b = *(d++)))
		    a = b;
	    }
	}
    }

    if (DofEdgeDataCount(x) > 0) {
	n = x->count_edge;
	for (i = 0; i < g->nedge; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->edge_index[i + 1] -
			      x->hp->edge_index[i]);
	    if (!(g->types_edge[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a < (b = *(d++)))
		    a = b;
	    }
	}
    }

    if (DofFaceDataCount(x) > 0) {
	n = x->count_face;
	for (i = 0; i < g->nface; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->face_index[i + 1] -
			      x->hp->face_index[i]);
	    if (!(g->types_face[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a < (b = *(d++)))
		    a = b;
	    }
	}
    }

    if (DofElementDataCount(x) > 0) {
	n = x->count_elem;
	for (i = 0; i < g->nelem; i++) {
	    if (DofIsHP(x))
		n = x->dim * (x->hp->elem_index[i + 1] -
			      x->hp->elem_index[i]);
	    if (!(g->types_elem[i] & OWNER)) {
		d += n;
		continue;
	    }
	    for (ii = 0; ii < n; ii++) {
		if (a < (b = *(d++)))
		    a = b;
	    }
	}
    }

#if USE_MPI
    MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_MAX, x->g->comm);
    return b;
#else	/* USE_MPI */
    return a;
#endif	/* USE_MPI */
}

FLOAT
phgDofNormL1_(DOF *u, int quad_order)
{
    int i, j, nvalues;
    FLOAT d, d0, tmp, sum;
    const FLOAT *v1, *w;
    QUAD *quad;
    ELEMENT *e;

    nvalues = DofDim(u);
    assert(u != NULL);
    
    sum = 0.;
    ForAllElements(u->g, e){
	if (quad_order >= 0)
	    i = quad_order;
	else if ((i = DofTypeOrder(u, e)) < 0)
	    i = 1;
	quad = phgQuadGetQuad3D(i);
        d = 0.;
        w = quad->weights;
        v1 = phgQuadGetDofValues(e, u, quad);
        for (i = 0; i < quad->npoints; i++){
            tmp = 0.;
            for (j = 0; j < nvalues; j++){
                d0 = *(v1++);
	        tmp += Fabs(d0);
	    }
            d += tmp * *(w++);
	}
	sum += d * phgGeomGetVolume(u->g, e);
    }

#if USE_MPI
    MPI_Allreduce(&sum, &d0, 1, PHG_MPI_FLOAT, PHG_SUM, u->g->comm);
#else
    d0 = sum;
#endif

    return d0;
}

FLOAT
phgDofNormL2_(DOF *u, int quad_order)
{
    ELEMENT *e;
    FLOAT norm;
#if USE_MPI
    FLOAT a, b;
#endif

    norm = 0.0;

#if USE_OMP
#pragma omp parallel for schedule(static) reduction(+:norm) private(e)
#endif	/* USE_OMP */
    ForAllElementsBegin(u->g, e)
	norm += phgQuadDofDotDof(e, u, u, quad_order);
    ForAllElementsEnd

#if USE_MPI
    a = norm;
    MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_SUM, u->g->comm);
    norm = Sqrt(b);
#else
    norm = Sqrt(norm);
#endif

    return norm;
}

FLOAT
phgDofDotL2_(DOF *u, DOF *v, int quad_order)
/* returns L2 inner product of u and v, v == NULL means v == 1 */
{
    ELEMENT *e;
    FLOAT udotv;
#if USE_MPI
    FLOAT a;
#endif
    BOOLEAN v_is_nil = FALSE;

    if (v == NULL) {
	v_is_nil = TRUE;
	v = phgDofNew(u->g, DOF_CONSTANT, DofDim(u), "one", DofNoAction);
	phgDofSetDataByValue(v, 1.0);
    }

    udotv = 0.0;
    ForAllElements(u->g, e)
	udotv += phgQuadDofDotDof(e, u, v, quad_order);

#if USE_MPI
    a = udotv;
    MPI_Allreduce(&a, &udotv, 1, PHG_MPI_FLOAT, PHG_SUM, u->g->comm);
#endif

    if (v_is_nil)
	phgDofFree(&v);

    return udotv;
}

FLOAT
phgDofNormH1_(DOF *u, int quad_order)
{
    ELEMENT *e;
    FLOAT norm;
    DOF *grad_u = phgDofGradient(u, NULL, NULL, NULL);
#if USE_MPI
    FLOAT a, b;
#endif

    norm = 0.0;
    ForAllElements(u->g, e)
	norm += phgQuadDofDotDof(e, u, u, quad_order) +
		phgQuadDofDotDof(e, grad_u, grad_u, quad_order);

#if USE_MPI
    a = norm;
    MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_SUM, u->g->comm);
    norm = Sqrt(b);
#else
    norm = Sqrt(norm);
#endif

    phgDofFree(&grad_u);

    return norm;
}

void
phgDofRandomize(DOF *u, long int seed)
/* Note: seed == 0 ==> no seed */
{
    GRID *g = u->g;
    INT j;
    FLOAT *data;
    int i, k, n;
    BTYPE btype;

    MagicCheck(DOF, u);

    phgDofClearCache(NULL, u, NULL, NULL, FALSE);

    if (seed != 0)
	srand48(seed + u->g->rank);
    data = DofData(u);
    if (DofVertexDataCount(u) > 0) {
	n = u->type->np_vert;
	for (j = 0; j < g->nvert; j++) {
	    if (DofIsHP(u))
		n = u->hp->vert_index[j + 1] - u->hp->vert_index[j];
	    btype = g->types_vert[j];
	    for (k = 0; k < u->dim; k++) {
		if (btype == UNREFERENCED || (btype & (u->DB_masks == NULL ?
						u->DB_mask : u->DB_masks[k]))) {
		    for (i = 0; i < n; i++)
			*(data++) = 0.0;
		}
		else {
		    for (i = 0; i < n; i++)
			*(data++) = 2.0 * drand48() - 1.0;
		}
	    }
	}
    }

    if (DofEdgeDataCount(u) > 0) {
	n = u->type->np_edge;
	for (j = 0; j < g->nedge; j++) {
	    if (DofIsHP(u))
		n = u->hp->edge_index[j + 1] - u->hp->edge_index[j];
	    btype = g->types_edge[j];
	    for (k = 0; k < u->dim; k++) {
		if (btype == UNREFERENCED || (btype & (u->DB_masks == NULL ?
						u->DB_mask : u->DB_masks[k]))) {
		    for (i = 0; i < n; i++)
		        *(data++) = 0.0;
	        }
	        else {
		    for (i = 0; i < n; i++)
		        *(data++) = 2.0 * drand48() - 1.0;
	        }
	    }
	}
    }

    if (DofFaceDataCount(u) > 0) {
	n = u->type->np_face;
	for (j = 0; j < g->nface; j++) {
	    if (DofIsHP(u))
		n = u->hp->face_index[j + 1] - u->hp->face_index[j];
	    btype = g->types_face[j];
	    for (k = 0; k < u->dim; k++) {
		if (btype == UNREFERENCED || (btype & (u->DB_masks == NULL ?
						u->DB_mask : u->DB_masks[k]))) {
		    for (i = 0; i < n; i++)
		        *(data++) = 0.0;
	        }
	        else {
		    for (i = 0; i < n; i++)
		        *(data++) = 2.0 * drand48() - 1.0;
	        }
	    }
	}
    }

    if (DofElementDataCount(u) > 0) {
	n = u->type->np_elem;
	for (j = 0; j < g->nelem; j++) {
	    if (DofIsHP(u))
		n = u->hp->elem_index[j + 1] - u->hp->elem_index[j];
	    btype = g->types_elem[j];
	    for (k = 0; k < u->dim; k++) {
		if (btype == UNREFERENCED || (btype & (u->DB_masks == NULL ?
						u->DB_mask : u->DB_masks[k]))) {
		    for (i = 0; i < n; i++)
		        *(data++) = 0.0;
	        }
	        else {
		    for (i = 0; i < n; i++)
		        *(data++) = 2.0 * drand48() - 1.0;
	        }
	    }
	}
    }
}

static void
check_FE_spaces(DOF *src0, DOF *src1, DOF *dest, const char *file, int line)
/* check whether the function 'src0*src1' belongs to the FE space of 'dest',
 * and issue a warning if not */
{
    DOF_TYPE *type0, *type1, *type2;
    int sorder, dorder;
    char msg[128];

    assert(dest != NULL);
    if ((type2 = dest->type) == NULL)
	type2 = dest->hp->info->types[dest->hp->max_order];

    if (src0 == NULL)
	type0 = NULL;
    else if ((type0 = src0->type) == NULL)
	type0 = src0->hp->info->types[src0->hp->max_order];

    if (src1 == NULL)
	type1 = NULL;
    else if ((type1 = src1->type) == NULL)
	type1 = src1->hp->info->types[src1->hp->max_order];

    msg[0] = '\0';

    /* first, check finite element spaces. */
    assert(type2->fe_space != FE_None && type2->fe_space != FE_Cinfty);
    switch (type2->fe_space) {
	case FE_H1:
	    if (type0 != NULL &&
		(type0->fe_space == FE_L2 || 
		 type0->fe_space == FE_Hcurl ||
		 type0->fe_space == FE_Hdiv)) {
		sprintf(msg, "Warning: %s:%d, transfering function from %s "
			  "to H1.", file, line,
			  phgDofFESpaceName(type0->fe_space));
		break;
	    }
	    if (type1 != NULL &&
		(type1->fe_space == FE_L2 || 
		 type1->fe_space == FE_Hcurl ||
		 type1->fe_space == FE_Hdiv)) {
		sprintf(msg, "Warning: %s:%d, transfering function from %s "
			  "to H1.", file, line,
			  phgDofFESpaceName(type1->fe_space));
		break;
	    }
	    break;
	case FE_Hcurl:
	    if (type0 != NULL &&
		(type0->fe_space == FE_L2 || type0->fe_space == FE_Hdiv)) {
		sprintf(msg, "Warning: %s:%d, transfering function from %s "
			  "to Hcurl.", file, line,
			  phgDofFESpaceName(type0->fe_space));
		break;
	    }
	    if (type1 != NULL &&
		(type1->fe_space == FE_L2 || type1->fe_space == FE_Hdiv)) {
		sprintf(msg, "Warning: %s:%d, transfering function from %s "
			  "to Hcurl.", file, line,
			  phgDofFESpaceName(type1->fe_space));
		break;
	    }
	    break;
	case FE_Hdiv:
	    if (type0 != NULL &&
		(type0->fe_space == FE_L2 || type0->fe_space == FE_Hcurl)) {
		sprintf(msg, "Warning: %s:%d, transfering function from %s "
			  "to Hdiv.", file, line,
			  phgDofFESpaceName(type0->fe_space));
		break;
	    }
	    if (type1 != NULL &&
		(type1->fe_space == FE_L2 || type1->fe_space == FE_Hcurl)) {
		sprintf(msg, "Warning: %s:%d, transfering function from %s "
			  "to Hdiv.", file, line,
			  phgDofFESpaceName(type1->fe_space));
		break;
	    }
	    break;
	default:
	    break;
    }

    /* next, check polynomial orders */
    if (msg[0] == '\0') {
	sorder = (type0 == NULL || type0 == DOF_ANALYTIC ? 0 : type0->order) +
		 (type1 == NULL || type1 == DOF_ANALYTIC ? 0 : type1->order);
	dorder = (type2 == NULL ? 0 : type2->order);
	if (sorder > dorder)
	    sprintf(msg, "Warning: %s:%d, order of the target DOF is too low.",
		    file, line);
    }

    if (msg[0] != '\0')
	phgPrintf("\n****\n**** %s\n****\n", msg);
}

static DOF *dofa, *dofb, *dofc;
static FLOAT *cache_bufb = NULL, *cache_bufc = NULL;
static int cache_dimb = 0, cache_dimc = 0;
static int M0, N0, K0, Kstart, Ksize, incA;
static FLOAT *bufa, *bufb, alpha0;
static MAT_OP ta, tb;

static void
mm_func(DOF *dof, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* temporary Init function used by phgDofMM */
{
    int i, j, k, m, k0;
    FLOAT tmp, *A, *B, *C;

    CheckThread
    CheckD_order(dof);
    if (dofa != NULL)
	phgDofEval_(dofa, e, lambda, cache_buf + bno * cache_dim, 0,
		    bufa, NULL);
    phgDofEval_(dofb, e, lambda, cache_bufb + bno * cache_dimb, 0, bufb, NULL);
    if (DofIsHP(dofc) || dofc->type->InitFunc != phgDofInitFuncPoint) {
	phgDofEval_(dofc, e, lambda, cache_bufc + bno * cache_dimc, 0,
			values, NULL);
    }
    else for (; TRUE;) {
	INT loc;
	static int *M = NULL, M_size = 0;
	FreeAtExit(M);
	/* faster code: don't need to go through phgDofEval */
	if (bno < (m = dofc->type->np_vert * NVert)) {
	    loc = (e->verts[bno / dofc->type->np_vert] * dofc->type->np_vert +
		   bno % dofc->type->np_vert) * dofc->dim;
	    memcpy(values, dofc->data_vert + loc, dofc->dim * sizeof(*values));
	    break;
	}
	bno -= m;
	if (bno < (m = dofc->type->np_edge * NEdge)) {
	    i = bno / dofc->type->np_edge;
	    if (M_size < dofc->type->np_edge) {
		phgFree(M);
		M = phgAlloc((M_size = dofc->type->np_edge) * sizeof(*M));
	    }
	    phgDofMapEdgeData(dofc->type, e, i, M);
	    loc = (e->edges[i] * dofc->type->np_edge +
		   M[bno % dofc->type->np_edge]) * dofc->dim;
	    memcpy(values, dofc->data_edge + loc, dofc->dim * sizeof(*values));
	    break;
	}
	bno -= m;
	if (bno < (m = dofc->type->np_face * NFace)) {
	    i = bno / dofc->type->np_face;
	    if (M_size < dofc->type->np_face) {
		phgFree(M);
		M = phgAlloc((M_size = dofc->type->np_face) * sizeof(*M));
	    }
	    phgDofMapFaceData(dofc->type, e, i, M);
	    loc = (e->faces[i] * dofc->type->np_face +
		   M[bno % dofc->type->np_face]) * dofc->dim;
	    memcpy(values, dofc->data_face + loc, dofc->dim * sizeof(*values));
	    break;
	}
	bno -= m;
	assert(bno < dofc->type->np_elem);
	if (M_size < dofc->type->np_elem) {
	    phgFree(M);
	    M = phgAlloc((M_size = dofc->type->np_elem) * sizeof(*M));
	}
	phgDofMapElementData(dofc->type, e, M);
	loc = (e->index * dofc->type->np_elem + M[bno]) * dofc->dim;
	memcpy(values, dofc->data_elem + loc, dofc->dim * sizeof(*values));
	break;
    }

    /* compute c += alpha * A * b */
    A = bufa;
    B = bufb;
    C = values;
    k0 = Kstart;
    if (dofa == NULL) {
	m = DofDim(dofc);
	assert(m == DofDim(dofb));
	if (tb == MAT_OP_N || M0 == 1) {
	    if (alpha0 == 1.) {
		for (i = 0; i < m; i++)
		    *(C++) += *(B++);
	    }
	    else if (alpha0 == -1.) {
		for (i = 0; i < m; i++)
		    *(C++) -= *(B++);
	    }
	    else {
		for (i = 0; i < m; i++)
		    *(C++) += *(B++) * alpha0;
	    }
	}
	else {	/* tb == MAT_OP_T */
	    WARN_UNTESTED
	    if (alpha0 == 1.) {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++)
			*(C++) += B[j * M0];
		    B++;
		}
	    }
	    else if (alpha0 == -1.) {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++)
			*(C++) -= B[j * M0];
		    B++;
		}
	    }
	    else {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++)
			*(C++) += B[j * M0] * alpha0;
		    B++;
		}
	    }
	}
    }
    else if (ta == MAT_OP_N) {
	if (tb == MAT_OP_N || M0 == 1) {
	    if (alpha0 == 1.) {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j + k0 * N0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k] * B[k * N0];
			*(C++) += tmp;
		    }
		    A += incA;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	    else if (alpha0 == -1.) {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j + k0 * N0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k] * B[k * N0];
			*(C++) -= tmp;
		    }
		    A += incA;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	    else {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j + k0 * N0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k] * B[k * N0];
			*(C++) += tmp * alpha0;
		    }
		    A += incA;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	}
	else {	/* tb == MAT_OP_T */
	    if (alpha0 == 1.) {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j * K0 + k0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k] * B[k];
			*(C++) += tmp;
		    }
		    A += incA;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	    else if (alpha0 == -1.) {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j * K0 + k0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k] * B[k];
			*(C++) -= tmp;
		    }
		    A += incA;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	    else {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j * K0 + k0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k] * B[k];
			*(C++) += tmp * alpha0;
		    }
		    A += incA;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	}
    }
    else {	/* ta == MAT_OP_T */
	if (tb == MAT_OP_N || M0 == 1) {
	    if (alpha0 == 1.) {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j + k0 * N0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k * incA] * B[k * N0];
			*(C++) += tmp;
		    }
		    A++;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	    else if (alpha0 == -1.) {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j + k0 * N0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k * incA] * B[k * N0];
			*(C++) -= tmp;
		    }
		    A++;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	    else {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j + k0 * N0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k * incA] * B[k * N0];
			*(C++) += tmp * alpha0;
		    }
		    A++;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	}
	else {	/* tb == MAT_OP_T */
	    if (alpha0 == 1.) {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j * K0 + k0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k * incA] * B[k];
			*(C++) += tmp;
		    }
		    A++;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	    else if (alpha0 == -1.) {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j * K0 + k0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k * incA] * B[k];
			*(C++) -= tmp;
		    }
		    A++;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	    else {
		for (i = 0; i < M0; i++) {
		    for (j = 0; j < N0; j++) {
			B = bufb + j * K0 + k0;
			tmp = A[0] * B[0];
			for (k = 1; k < Ksize; k++)
			    tmp += A[k * incA] * B[k];
			*(C++) += tmp * alpha0;
		    }
		    A++;
		    if (M0 == K0 && i % Ksize == Ksize - 1)
			k0 += Ksize;
		}
	    }
	}
    }
}

DOF *
phgDofMM_(MAT_OP transa, MAT_OP transb, int M, int N, int K,
	  FLOAT alpha, DOF *A, int blka, DOF *B, FLOAT beta, DOF **Cptr,
	  const char *srcfile, int srcline, BOOLEAN check)
/******************************************************************************
 * DGEMM-like matrix-matrix operation between individual DOF:
 *
 *	C := alpha * op(A) * op(B) + beta * C
 *
 * (array dimensions: M K x K N = M N)
 *
 * Parameters:
 *
 *  transa - if transa == MAT_OP_T then op(A) = A', otherwise op(A) = A
 *  transb - if transb == MAT_OP_T then op(B) = B', otherwise op(B) = B
 *  M -	specifies  the number of rows of the matrix op(A) and of the matrix C
 *  N -	specifies the number of columns of the matrix op(B) and the number of
 *	columns of the matrix C
 *  K -	specifies the number of columns of the matrix op(A) and the number of
 *	rows of the matrix op(B)
 *  A -	array of dimension [M][K] when transa == MAT_OP_N, and of dimension
 *	[K][M] otherwise
 *  blka - blocking size for A. The following conditions must be satisfied:
 *	(1) If blka <= 0 then (A == NULL and M == K) or DofDim(A) == M * K.
 *	(2) If blka > 0 then M == K and M % blka == 0 and
 *	    DofDim(A) % (blka * blka) == 0 and DofDim(A) / blka == M,
 *	    in this case A is interpreted as a block diagonal matrix with
 *	    M/blka blocks and block size blka,
 *  B -	array of dimension [K][N] when transb == MAT_OP_N, and of dimension
 *	[N][K] otherwise
 *  C -	array of dimension [M][N]
 *
 * Note: all arrays are stored in C (i.e., row-major) order
 *****************************************************************************/
{
    DOF *C = (Cptr == NULL ? NULL : *Cptr);
    INT i, m;
    int ii, jj;
    FLOAT *p, *q;

    MagicCheck(DOF, A);
    MagicCheck(DOF, B);
    MagicCheck(DOF, C);

    assert(C == NULL || (C != A && C != B));
    assert(A != NULL || M == K);
    assert(A == NULL || (blka <= 0 && DofDim(A) == M * K) ||
	   (blka > 0 && M == K && DofDim(A) == M * blka && M % blka == 0));
    assert(B == NULL || DofDim(B) == K * N);
    assert(C == NULL || DofDim(C) == M * N);

    if (C == NULL) {
	assert(B != NULL && (M * N) % DofTypeDim(B) == 0);
	beta = 0.;
	ii = M * N / DofTypeDim(B);
	C = phgDofNew_(B->g, B->type, B->hp, ii, "C", DofNoAction,
				srcfile, srcline);
	C->DB_mask = B->DB_mask;
	if (Cptr != NULL)
	    *Cptr = C;
    }
    else if (beta == 0.) {
	phgDofSetDataByValue(C, (FLOAT)0.);
    }
    else if (beta == -1.) {
	m = DofDataCount(C);
	p = DofData(C);
	for (i = 0; i < m; i++, p++)
	    *p = - *p;
    }
    else if (beta != 1.) {
	m = DofDataCount(C);
	p = DofData(C);
	for (i = 0; i < m; i++)
	    *(p++) *= beta;
    }

    CheckD_order(C);
    if (C->DB_masks != NULL)
	memset(C->DB_masks, 0, C->dim * sizeof(*C->DB_masks));
    if (C->userfunc != DofInterpolation)
	C->userfunc = DofNoAction;
    C->userfunc_lambda = NULL;

    /* warn user if transfering to a smaller FE space */
    if (check)
	check_FE_spaces(A, B, C, srcfile, srcline);

    phgDofClearCache(NULL, C, NULL, NULL, FALSE);

    if (alpha == 0.)
	return C;

    assert(B != NULL);

    if (A == NULL &&
	((!DofIsHP(B) && B->type == C->type) /*|| B->hp === C->hp */)) {
	m = DofDataCount(C);
	p = DofData(B);
	q = DofData(C);
	if (transb == MAT_OP_N) {
	    if (alpha == 1.) {
		for (i = 0; i < m; i++)
		    *(q++) += *(p++);
	    }
	    else if (alpha == -1.) {
		for (i = 0; i < m; i++)
		    *(q++) -= *(p++);
	    }
	    else {
		for (i = 0; i < m; i++)
		    *(q++) += *(p++) * alpha;
	    }
	}
	else {
	    WARN_UNTESTED
	    m /= N * M;
	    if (alpha == 1.) {
		for (i = 0; i < m; i++) {
		    for (ii = 0; ii < M; ii++) {
			for (jj = 0; jj < N; jj++)
			    *(q++) += p[jj * M];
			p++;
		    }
		}
	    }
	    else if (alpha == -1.) {
		for (i = 0; i < m; i++) {
		    for (ii = 0; ii < M; ii++) {
			for (jj = 0; jj < N; jj++)
			    *(q++) -= p[jj * M];
			p++;
		    }
		}
	    }
	    else {
		for (i = 0; i < m; i++) {
		    for (ii = 0; ii < M; ii++) {
			for (jj = 0; jj < N; jj++)
			    *(q++) += p[jj * M] * alpha;
			p++;
		    }
		}
	    }
	}
	return C;
    }

    /* pre-evaluation of basis functions at nodal points of C */
    if (!DofIsHP(C) && C->type->is_nodal) {
	ii = C->type->nbas;
	if (A != NULL && !DofIsHP(A) && A->type->BasFuncs != NULL &&
	    A->type->invariant) {
	    cache_dim = A->type->nbas * A->type->dim;
	    cache_buf = phgAlloc(cache_dim * ii * sizeof(*cache_buf));
	    get_bas_funcs(A, C, C->g->roots, cache_buf);
	}
	if (!DofIsHP(B) && B->type->BasFuncs != NULL &&
	    B->type->invariant) {
	    cache_dimb = B->type->nbas * B->type->dim;
	    cache_bufb = phgAlloc(cache_dimb * ii * sizeof(*cache_bufb));
	    get_bas_funcs(B, C, C->g->roots, cache_bufb);
	}
	if (!DofIsHP(C) && C->type->BasFuncs != NULL &&
	    C->type->invariant && C->type->InitFunc != phgDofInitFuncPoint) {
	    cache_dimc = C->type->nbas * C->type->dim;
	    cache_bufc = phgAlloc(cache_dimc * ii * sizeof(*cache_bufc));
	    get_bas_funcs(C, C, C->g->roots, cache_bufc);
	}
    }

    M0 = M;
    N0 = N;
    K0 = K;
    if (blka == M && M == K)
	blka = 0;
    if (blka <= 0) {
	Kstart = 0;
	Ksize = K;
	incA = (transa == MAT_OP_N ? K : M);
    }
    else {
	Kstart = 0;
	Ksize = blka;
	incA = blka;
    }
    alpha0 = alpha;
    ta = transa;
    tb = transb;
    bufa = (A == NULL ? NULL : phgAlloc(DofDim(A) * sizeof(*bufa)));
    bufb = phgAlloc(DofDim(B) * sizeof(*bufb));
    dofa = A;
    dofb = B;
    dofc = phgDofCopy(C, NULL, NULL, NULL);
    phgDofSetDataByFunction_(C, NULL, mm_func, NULL);
    phgDofFree(&dofc);
    phgFree(bufa);
    phgFree(bufb);

    if (cache_buf != NULL) {
	phgFree(cache_buf);
	cache_buf = NULL;
	cache_dim = 0;
    }
    if (cache_bufb != NULL) {
	phgFree(cache_bufb);
	cache_bufb = NULL;
	cache_dimb = 0;
    }
    if (cache_bufc != NULL) {
	phgFree(cache_bufc);
	cache_bufc = NULL;
	cache_dimc = 0;
    }

    return C;
}

DOF *
phgDofAXPBY_(FLOAT a, DOF *x0, FLOAT b, DOF **y_ptr,
			const char *srcfile, int srcline, BOOLEAN check)
/* y := a * x + b * y.
 * Note: only works when dofvalues are linear function of funcvalues
 * otherwise need to evaluate function values using dof->type->InitFunc
 */
{
#if 0
    int m = DofDim(x0);
    return phgDofMM_(MAT_OP_N, MAT_OP_N, m, 1, m, a, NULL, 0,
					x0, b, y_ptr, srcfile, srcline, check);
#else
    DOF *x = x0, *y = (y_ptr == NULL ? NULL : *y_ptr);
    INT i, n = 0;
    FLOAT *p;

    MagicCheck(DOF, y);

    assert(y == NULL || y->g == x->g);
    assert(y == NULL || y->type != DOF_ANALYTIC);
    assert(y == NULL || y != x);

    if (y == NULL) {
	assert(x0 != NULL && x0->type != DOF_ANALYTIC);
	y = phgDofNew_(x0->g, x0->type, x0->hp, x0->dim, NULL, DofNoAction,
				srcfile, srcline);
	y->DB_mask = x0->DB_mask;
	if (y_ptr != NULL)
	    *y_ptr = y;
    }
    else if (b == 0.) {
	phgDofSetDataByValue(y, (FLOAT)0.);
    }
    else if (b == -1.) {
	n = DofDataCount(y);
	p = DofData(y);
	for (i = 0; i < n; i++, p++)
	    *p = - *p;
    }
    else if (b != 1.) {
	n = DofDataCount(y);
	p = DofData(y);
	for (i = 0; i < n; i++)
	    *(p++) *= b;
    }

    if (y->DB_masks != NULL)
	memset(y->DB_masks, 0, y->dim * sizeof(*y->DB_masks));
    if (y->userfunc != DofInterpolation)
	y->userfunc = DofNoAction;
    y->userfunc_lambda = NULL;

    if (x != NULL && DofDim(y) != DofDim(x))
	phgError(1, "%s:%d: phgDofAXPBY: incompatible DOFs \"%s\" <-> "
		    "\"%s\".\n", srcfile, srcline, x->name, y->name);

    phgDofClearCache(NULL, y, NULL, NULL, FALSE);

    if (a == 0.0)
	return y;			/* nothing to do */

    /* warn user if transfering to a smaller FE space */
    if (check)
	check_FE_spaces(x, NULL, y, srcfile, srcline);

    if (y->type != x->type || y->hp != x->hp) {
	x = phgDofNew_(y->g, y->type, y->hp, y->dim, "AXPBY_tmp", DofNoAction,
			__FILE__, __LINE__);
	phgDofCopy(x0, &x, NULL, NULL);
    }

    if ((n = DofDataCount(y)) == 0) {
	if (x != x0)
	    phgDofFree(&x);
	return y;
    }

    if (a == -1.0) {		/* y := y - x */
	for (i = 0; i < n; i++) {
	    y->data[i] -= x->data[i];
	}
    }
    else if (a == 1.0) {	/* y := y + x */
	for (i = 0; i < n; i++) {
	    y->data[i] += x->data[i];
	}
    }
    else {			/* y := y + a * x */
	for (i = 0; i < n; i++) {
	    y->data[i] += a * x->data[i];
	}
    }

    if (x != x0)
	phgDofFree(&x);

    return y;
#endif
}

DOF *
phgDofMatVec_(FLOAT alpha, DOF *A, DOF *x, FLOAT beta, DOF **y_ptr,
		const char *srcfile, int srcline, BOOLEAN check)
/* BLAS-like function computing y := alpha A x + beta y */
{
    DOF *y = (y_ptr == NULL ? NULL : *y_ptr);
    INT M;

    Unused(y);

    if (A == NULL || alpha == 0.)
	return phgDofAXPBY_(alpha, x, beta, y_ptr, srcfile, srcline, check);

    assert(x != NULL);
    M = DofDim(x);
    assert(y == NULL || M == DofDim(y));
    if (DofDim(A) == 1)
	return phgDofMM_(MAT_OP_N, MAT_OP_N, 1, M, 1,
			alpha, A, 0, x, beta, y_ptr, srcfile, srcline, check);
    if (DofDim(A) == M)
	return phgDofMM_(MAT_OP_N, MAT_OP_N, M, 1, M,
			alpha, A, 1, x, beta, y_ptr, srcfile, srcline, check);
    assert(DofDim(A) == M * M); 
    return phgDofMM_(MAT_OP_N, MAT_OP_N, M, M, M,
			alpha, A, 0, x, beta, y_ptr, srcfile, srcline, check);
}

static FLOAT *AFXPBY_buffer = NULL;
static void (*AFXPBY_f)(FLOAT *in, FLOAT *out) = NULL;

static void
AFXPBY_func(DOF *dof, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* temporary Init function used by phgDofAFXPBY_ */
{
    CheckD_order(dof);
    phgDofEval_(cache_dof, e, lambda, cache_buf + bno * cache_dim,
		0, AFXPBY_buffer, NULL);
    AFXPBY_f(AFXPBY_buffer, values);
}

DOF *
phgDofAFXPBY_(FLOAT a, void (*f)(FLOAT *in, FLOAT *out), DOF *x, FLOAT b,
	      DOF **yptr, const char *srcfile, int srcline)
/* computes and returns y := a f(x) + b y */
{
    DOF *tmp, *y;

    if (f == NULL || a == 0.0 || x == NULL)
	return phgDofAXPBY_(a, x, b, yptr, srcfile, srcline, TRUE);

    /* compute tmp := f(x) */
    y = (yptr == NULL || *yptr == NULL ? x : *yptr);
    assert(!SpecialDofType(y->type));
    tmp = phgDofNew_(y->g, y->type, y->hp, y->dim, "AFXPBY_tmp", DofNoAction,
			 srcfile, srcline);

    if (!DofIsHP(x) && !SpecialDofType(x->type) &&
	x->type->invariant == TRUE && !DofIsHP(tmp) && tmp->type->is_nodal) {
	cache_dim = x->type->nbas * x->type->dim;
	cache_buf = phgAlloc(tmp->type->nbas * cache_dim * sizeof(*cache_buf));
	get_bas_funcs(x, tmp, x->g->roots, cache_buf);
    }
    else {
	cache_buf = NULL;
	cache_dim = 0;
    }
    cache_dof = x;

    AFXPBY_buffer = phgAlloc(DofDim(x) * sizeof(*AFXPBY_buffer));
    AFXPBY_f = f;
    phgDofSetDataByFunction_(tmp, NULL, AFXPBY_func, NULL);
    AFXPBY_f = NULL;
    phgFree(AFXPBY_buffer);
    AFXPBY_buffer = NULL;

    if (cache_buf != NULL) {
	phgFree(cache_buf);
	cache_buf = NULL;
    }
    cache_dim = 0;

    y = phgDofAXPBY_(a, tmp, b, yptr, srcfile, srcline, TRUE);
    phgDofFree(&tmp);

    return y;
}

static FLOAT (*f_tmp)(FLOAT xvalue) = NULL;
static int nin = 0;

static void
f_dummy(FLOAT *in, FLOAT *out)
{
    int i;

    for (i = 0; i < nin; i++)
	*(out++) = f_tmp(*(in++));
}

DOF *
phgDofAFXPBY1_(FLOAT a, FLOAT (*f)(FLOAT xvalue), DOF *x,
	       FLOAT b, DOF **yptr, const char *srcfile, int srcline)
{
    f_tmp = f;
    nin = (x == NULL ? 0 : DofDim(x));
    x = phgDofAFXPBY_(a, f_dummy, x, b, yptr, srcfile, srcline);
    f_tmp = NULL;
    nin = 0;
    return x;
}

static int
eigen_solver(MAT *A, MAT *B, int n, int which, FLOAT tau, int *nit,
		  FLOAT *evals, void *ssd, MAP *map, DOF **u, va_list ap)
{
    VEC *vec = NULL;
    int ndof;
    INT i, nglobal;
    DOF ***dofs, **dofs0;
    BOOLEAN remove_bdry;

    GetVariableArgs0(ndof, dofs, u, DOF **, ap);
    assert(ndof > 0);
    if (map == NULL) {
	dofs0 = phgAlloc(ndof * sizeof(*dofs0));
	for (i = 0; i < ndof; i++)
	    dofs0[i] = dofs[i][0];
	map = phgMapCreateN(ndof, dofs0);
	phgFree(dofs0);	/* note: map will be freed if dofs0 != NULL */
    }
    else {
	dofs0 = NULL;
	nglobal = 0;
        for (i = 0; i < ndof; i++)
            nglobal += DofDataCountGlobal(dofs[i][0]);
	if (nglobal != map->nglobal)
	    phgError(1, "%s: inconsistent map argument!\n", __func__);
    }
    remove_bdry = (map->nglobal - map->bdry_nglobal == A->cmap->nglobal);
    phgMapDofArraysToVecN(map, n, remove_bdry, &vec, ndof, dofs);
    n = phgEigenSolve_(A, B, n, which, tau, evals, &vec, nit, ssd);
    phgMapVecToDofArraysN(map, vec, remove_bdry, ndof, dofs);
    phgVecDestroy(&vec);
    phgFree(dofs);
    if (dofs0 != NULL)
	phgMapDestroy(&map);

    return n;
}

/* new interface */
int
phgDofEigenSolve_(MAT *A, MAT *B, int n, int which, FLOAT tau, int *nit,
		  FLOAT *evals, void *ssd, MAP *map, DOF **u, ...)
{
    va_list ap;
    int ret;

    va_start(ap, u);
    ret = eigen_solver(A, B, n, which, tau, nit, evals, ssd, map, u, ap);
    va_end(ap);

    return ret;
}

/* old interface */
int
phgDofEigenSolve(MAT *A, MAT *B, int n, int which, FLOAT tau, int *nit,
		 FLOAT *evals, MAP *map, DOF **u, ...)
{
    va_list ap;
    int ret;

    va_start(ap, u);
    ret = eigen_solver(A, B, n, which, tau, nit, evals, NULL, map, u, ap);
    va_end(ap);

    return ret;
}

DOF *
phgDofGetSameOrderDG_(DOF *dof, int dim, const char *name,
				const char *srcfile, int srcline)
/* returns a DG DOF of the same order (hp or non-hp) as 'dof'.
 * If 'dim' < 0 then 'dim' = DofDim('dof') */
{
    HP_TYPE *hp;

    if (dim < 0)
	dim = DofDim(dof);

    if (!DofIsHP(dof))
	return phgDofNew_(dof->g, DOF_DGn[dof->type->order], NULL, dim, name,
			  DofNoAction, srcfile, srcline); 

    if (dof->hp->info == HP_DG)
	return phgDofNew_(dof->g, NULL, dof->hp, dim, name, DofNoAction,
			  srcfile, srcline);

    /* the next two lines will validate dof->hp->dg */
    hp = phgHPGetDG_(dof->hp, srcfile, srcline);
    phgHPFree(&hp);

    return phgDofNew_(dof->g, NULL, dof->hp->dg, dim, name, DofNoAction,
		      srcfile, srcline);
}

/*---------------------------------------------------------------------------*/

#if 0
void
phgDofMapDataInElement(DOF_TYPE *type, const ELEMENT *e, int dim, FLOAT *values)
/* The buffer values[] contains data (e.g., values of all bases functions) in
 * the order vertices/edges/faces/element. This function rearranges data on
 * according to their global orientation. */
{
    FLOAT *buffer, *p;
    int i, k, nbas, *M;

    if (type == NULL || type->points == NULL ||
	(type->np_edge <= 1 && type->np_face <= 1))
	return;

    assert(type->BasFuncs != NULL);

    nbas = type->nbas;

    i = (type->np_edge < type->np_face ? type->np_face : type->np_edge);
    i = (i < type->np_elem ? type->np_elem : i);
    M = phgAlloc(i * sizeof(*M));
    buffer = phgAlloc(dim * i * sizeof(*values));

    p = values + NVert * dim * type->np_vert;

    if (type->np_edge <= 1) {
	p += NEdge * dim * type->np_edge;
    }
    else {
	for (i = 0; i < NEdge; i++) {
	    phgDofMapEdgeData(type, e, i, M);
	    memcpy(buffer, p, dim * type->np_edge * sizeof(*p));
	    for (k = 0; k < type->np_edge; k++)
		memcpy(p + dim * M[k], buffer + dim * k, dim * sizeof(*p));
	    p += dim * type->np_edge;
	}
    }

    if (type->np_face > 1) {
	for (i = 0; i < NFace; i++) {
	    phgDofMapFaceData(type, e, i, M);
	    memcpy(buffer, p, dim * type->np_face * sizeof(*p));
	    for (k = 0; k < type->np_face; k++)
		memcpy(p + dim * M[k], buffer + dim * k, dim * sizeof(*p));
	    p += dim * type->np_face;
	}
    }

    /* TODO: reorder element data */

    phgFree(M);
    phgFree(buffer);

    return;
}
#endif	/* 0 */

static void
func_tr(DOF *v, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* dummy function for computing P.
 * Note: the results stored as values[v->dim][type->dim] with v->dim == nbas */
{
    const FLOAT *p;
    DOF_TYPE *type = dofa->type;

    Unused(bno);
    CheckD_order(v);

    if (type == NULL)
	type = dofa->hp->info->types[dofa->hp->elem_order[e->index]];

    p = type->BasFuncs(dofa, e, 0, -1, lambda);
    memcpy(values, p, type->dim * v->dim * sizeof(*p));
}

MAT *
phgDofGetTransferMatrix(DOF *u, DOF *v, MAP *map_u, MAP *map_v)
/* Builds the transfer matrix P bwtween two DOFs, whose elements are given by:
 *	\phi_j = \sum m_{k,j} \psi_k
 * where \phi_j is the j-th basis of u, \psi_k the k-th basis of v.
 *
 * map_u and map_v are maps for u and v respectively. They are optional and
 * both or either of them can be NULL. If one or both of them are NULL then
 * new maps are created accordingly.
 *
 * Note: this function should normally be called with the FE space of u being
 * a subspace of the FE space of v.
 *
 * Note: the present code was imported from solver-asp.c r1.10,
 *	 old code is available in dof-utils.c r1.180
 */
{
    MAT *P;
    DOF_TYPE *type;
    ELEMENT e0, *e, *e1;
    GRID *g = u->g;
    DOF *dof_tmp;
    INT *cols, *rows;
    FLOAT *buffer, *p;
    int i, j, M, N, flag;
    BOOLEAN free_map_u, free_map_v;

    dofa = u;

    if (map_u == NULL) {
	free_map_u = TRUE;
	map_u = phgMapCreate(u, NULL);
    }
    else {
	free_map_u = FALSE;
	assert(map_u->dofs != NULL &&
		map_u->dofs[0] == u && map_u->dofs[1] == NULL);
    }

    if (map_v == NULL) {
	free_map_v = TRUE;
	map_v = phgMapCreate(v, NULL);
    }
    else {
	free_map_v = FALSE;
	assert(map_v->dofs != NULL &&
		map_v->dofs[0] == v && map_v->dofs[1] == NULL);
    }

    if (DofIsHP(u)) {
	M = u->hp->info->types[u->hp->max_order]->nbas;
    }
    else {
	M = u->type->nbas;
    }

    if (DofIsHP(v)) {
	dof_tmp = phgHPDofNew(g, v->hp, 1, "dof_tmp", DofNoData);
	N = v->hp->info->types[v->hp->max_order]->nbas;
    }
    else {
	dof_tmp = phgDofNew(g, v->type, 1, "dof_tmp", DofNoData);
	N = v->type->nbas;
    }

    P = phgMapCreateMat(map_v, map_u);
    phgMatSetMode(P, PHG_REPLACE);

    buffer = phgAlloc(N * M * sizeof(*buffer));
    rows = phgAlloc((N + M) * sizeof(*rows));
    cols = rows  + N;

    flag = 0;
    memset(&e0, 0, sizeof(e0));

    ForAllElements(g, e1) {
	e = e1;
	if (DofIsHP(v)) {
	    int k, row = 0, col = 0;

	    M = DofNBas(u, e);
	    N = DofNBas(v, e);

#if 0
	    /* force recompute the matrix for each element (for debugging) */
	    Unused(k);
	    Unused(row);
	    Unused(col);
	    flag = 0;
#else
	    assert(DofIsHP(u) || !strcmp(u->type->name, "P1") ||
		   !memcmp(u->type->name, "HB", 2));

	    for (j = 0; j < N; j++)
		rows[j] = phgMapE2L_(map_v, 0, e, j, FALSE, FALSE);
	    for (j = 0; j < M; j++)
		cols[j] = phgMapE2L_(map_u, 0, e, j, TRUE, FALSE);

	    for (i = 0; i < NVert; i++) {
		phgMatSetEntry(P, rows[row], cols[col], (FLOAT)1.0);
		row++;
		col++;
	    }

	    for (i = 0; i < NEdge; i++) {
		j = v->hp->info->types[v->hp->edge_order[e->edges[i]]]->np_edge;
		k = (!DofIsHP(u) ? u->type :
		    u->hp->info->types[v->hp->edge_order[e->edges[i]]])
			->np_edge;
		assert(k <= j);
		while (j-- > 0) {
		    if (k-- > 0) {
			phgMatSetEntry(P, rows[row], cols[col], _F(1.));
			col++;
		    }
		    row++;
		}
	    }

	    for (i = 0; i < NFace; i++) {
		j = v->hp->info->types[v->hp->face_order[e->faces[i]]]->np_face;
		k = (!DofIsHP(u) ? u->type :
		    u->hp->info->types[v->hp->face_order[e->faces[i]]])
			->np_face;
		assert(k <= j);
		while (j-- > 0) {
		    if (k-- > 0) {
			phgMatSetEntry(P, rows[row], cols[col], _F(1.0));
			col++;
		    }
		    row++;
		}
	    }

	    j = v->hp->info->types[v->hp->elem_order[e->index]]->np_elem;
	    k = (!DofIsHP(u) ? u->type :
		u->hp->info->types[v->hp->elem_order[e->index]])
		    ->np_elem;
	    assert(k <= j);
	    while (j-- > 0) {
		if (k-- > 0) {
		    phgMatSetEntry(P, rows[row], cols[col], (FLOAT)1.0);
		    col++;
		}
		row++;
	    }

	    continue;
#endif
	}
	else if (v->type->invariant && !ElementIsInOrder(e1)) {
	    /* reorder the vertices of e1 in increasing order */
	    int k, m;
	    e = &e0;
	    e->elem_type = e1->elem_type;
	    e->ordering = 0 | (1 << 3) | (2 << 6) | (3 << 9);
	    e->index = e1->index;
	    m = 0;
	    for (j = 0; j < NVert; j++) {
		i = MapVertex(e1,j);
		e->verts[i] = e1->verts[j];
		e->faces[i] = e1->faces[j];
		for (k = j + 1; k < NVert; k++)
		    e->edges[GetEdgeNo(i,MapVertex(e1,k))] = e1->edges[m++];
	    }
	}

	for (j = 0; j < N; j++)
	    rows[j] = phgMapE2L_(map_v, 0, e, j, FALSE, FALSE);
	for (j = 0; j < M; j++)
	    cols[j] = phgMapE2L_(map_u, 0, e, j, TRUE, FALSE);

	type = v->type;

	if (flag == 0) {
	    /* compute the element transfer matrix */
	    FLOAT *pdata[15];

	    if (v->type->invariant)
		flag = 1;		/* only compute element matrix once */
	    dof_tmp->dim = M;		/* this is hacky! */
	    p = buffer;

	    if (DofIsHP(v))
		type = v->hp->info->types[v->hp->min_order];
	    j = M * type->np_vert;
	    for (i = 0; i < NVert; i++) {
		type->InitFunc(dof_tmp, e, VERTEX, i, NULL, func_tr,
							NULL, p, pdata);
		pdata[i] = p;
		p += j;
	    }

	    if (!DofIsHP(v))
		j = M * type->np_edge;
	    for (i = 0; i < NEdge; i++) {
		if (DofIsHP(v)) {
		    type = v->hp->info->types[v->hp->edge_order[e->edges[i]]];
		    j = M * type->np_edge;
		}
		type->InitFunc(dof_tmp, e, EDGE, i, NULL, func_tr,
							NULL, p, pdata);
		pdata[4 + i] = p;
		p += j;
	    }

	    if (!DofIsHP(v))
		j = M * type->np_face;
	    for (i = 0; i < NFace; i++) {
		if (DofIsHP(v)) {
		    type = v->hp->info->types[v->hp->face_order[e->faces[i]]];
		    j = M * type->np_face;
		}
		if (j == 0)
		    continue;
		type->InitFunc(dof_tmp, e, FACE, i, NULL, func_tr,
							NULL, p, pdata);
		pdata[10 + i] = p;
		p += j;
	    }

	    if (DofIsHP(v))
		type = v->hp->info->types[v->hp->elem_order[e->index]];
	    if (type->np_elem > 0)
		type->InitFunc(dof_tmp, e, VOLUME, 0, NULL, func_tr,
							NULL, p, pdata);
	    dof_tmp->dim = 1;		/* restore dof_tmp->dim */
	}

	phgMatSetEntries(P, N, rows, M, cols, buffer);

#if 0
	/* (for debug only) */
	if (e != e1 && DofIsHP(v)) {
	    /* compare data obtained from the original and permuted element */
	    FLOAT *buffer1 = buffer;
	    INT *cols1 = cols, *rows1 = rows;
	    buffer = phgAlloc(N * M * sizeof(*buffer));
	    rows = phgAlloc((N + M) * sizeof(*rows));
	    cols = rows  + N;
	    e = e1;
	    flag = 0;
/*-------------- These lines are copied from above ------------------*/
	for (j = 0; j < N; j++)
	    rows[j] = phgMapE2L_(map_v, 0, e, j, FALSE);
	for (j = 0; j < M; j++)
	    cols[j] = phgMapE2L_(map_u, 0, e, j, TRUE);

	if (flag == 0) {
	    FLOAT *pdata[15];

	    flag = 1;			/* only recompute once */
	    type = v->type;
	    dof_tmp->dim = M;	/* this is hacky! */
	    p = buffer;

	    if (DofIsHP(v))
		type = v->hp->info->types[v->hp->min_order];
	    j = M * type->np_vert;
	    for (i = 0; i < NVert; i++) {
		type->InitFunc(dof_tmp, e, VERTEX, i, NULL, func_tr,
							NULL, p, pdata);
		pdata[i] = p;
		p += j;
	    }

	    if (!DofIsHP(v))
		j = M * type->np_edge;
	    for (i = 0; i < NEdge; i++) {
		if (DofIsHP(v)) {
		    type = v->hp->info->types[v->hp->edge_order[e->edges[i]]];
		    j = M * type->np_edge;
		}
		type->InitFunc(dof_tmp, e, EDGE, i, NULL, func_tr,
							NULL, p, pdata);
		pdata[4 + i] = p;
		p += j;
	    }

	    if (!DofIsHP(v))
		j = M * type->np_face;
	    for (i = 0; i < NFace; i++) {
		if (DofIsHP(v)) {
		    type = v->hp->info->types[v->hp->face_order[e->faces[i]]];
		    j = M * type->np_face;
		}
		if (j == 0)
		    continue;
		type->InitFunc(dof_tmp, e, FACE, i, NULL, func_tr,
							NULL, p, pdata);
		pdata[10 + i] = p;
		p += j;
	    }

	    if (DofIsHP(v))
		type = v->hp->info->types[v->hp->elem_order[e->index]];
	    if (type->np_elem > 0)
		type->InitFunc(dof_tmp, e, VOLUME, 0, NULL, func_tr,
							NULL, p, pdata);
	    dof_tmp->dim = 1;		/* restore dof_tmp->dim */
	}
/*-------------------------------------------------------------------*/

	    //i = 4 * (v->type->np_vert + v->type->np_face) +
		//6 * v->type->np_edge;
	    i = 0;
	    for (; i < N; i++) {
		int i1;
		for (i1 = 0; i1 < N; i1++)
		    if (rows[i] == rows1[i1])
			break;
		assert(i1 < N);
		//j = 4 * (u->type->np_vert + u->type->np_face) +
		    //6 * u->type->np_edge;
		j = 0;
		for (; j < M; j++) {
		    int j1;
		    double v, v1;
		    for (j1 = 0; j1 < M; j1++)
			if (cols[j] == cols1[j1])
			    break;
		    assert(j1 < M);
		    v = buffer[i * M + j];
		    v1 = buffer1[i1 * M + j1];
		    if (Fabs(v - v1) / (Fabs(v) + Fabs(v1) + 1.) > 1e-5)
			phgError(1, "P(%d,%d): (%g) - (%g) = %g\n",
						i, j, v, v1, v - v1);
		}
	    }
	    phgFree(buffer);
	    phgFree(rows);
	    buffer = buffer1;
	    rows = rows1;
	    cols = cols1;
	}
#endif
    }	/* ForAllElements */

    phgMatAssemble(P);

    phgFree(rows);
    phgFree(buffer);
    phgDofFree(&dof_tmp);
    if (free_map_u)
	phgMapDestroy(&map_u);
    if (free_map_v)
	phgMapDestroy(&map_v);
 
    return P;
}

/*---------------------------------------------------------------------------*/
