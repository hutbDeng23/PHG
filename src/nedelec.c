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

/* Edge elements.
 *
 * $Id: nedelec.c,v 1.52 2022/04/09 01:28:34 zlb Exp $
 */

#include "phg.h"

#include <math.h>
#include <string.h>
#include <strings.h>


/*----------------------------------*/
#if 0
static FLOAT *
lambda_var(GRID *g, ELEMENT *e, int no, FLOAT *var)
/*  Suppose:
 *   lambda = a * x + b * y + c * z + d;
 *   Return value: var
 *    var[0]=a,var[1]=b,var[2]=c;
 */
{
  int i, v0;
  COORD *c;
  FLOAT x[NVert], y[NVert], z[NVert];
  FLOAT *v = var;
  for (i = 0; i < NVert; i++)
    {
      v0 = e->verts[i];
      x[i] = (*(c = g->verts + v0))[0];
      y[i] = (*c)[1];
      z[i] = (*c)[2];
    }
/*	
	*(v++)=-(y[2]-y[1])*(z[3] - z[1]) + (y[3]-y[2])*(z[2] - z[1]);
	*(v++)=(x[2]-x[1])*(z[3] - z[1]) - (x[3]-x[1])*(z[2] - z[1]);
	*(v++)=-(x[2]-x[1])*(y[3] - y[1]) + (x[3]-x[1])*(y[2] - y[1]);
	*(v++)=x[1]*y[2]*z[3] +x[2]*y[3]*z[1] + x[3]*y[1]*z[2] -x[3]*y[2]*z[1] - x[2]*y[1]*z[3] +x[1]* y[3]*z[2]; 
*/

  get_lambda_var(1, 2, 3);
  get_lambda_var(0, 3, 2);
  get_lambda_var(0, 1, 3);
  get_lambda_var(0, 2, 1);

  return var;
}
#endif

static const FLOAT *
ND1_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of basis functions */
{
  static FLOAT values[NEdge][Dim], *v;
#if USE_OMP
#pragma omp threadprivate(values, v)
#endif  /* USE_OMP */
  int I, J, k;
  GRID *g = dof->g;
  FLOAT (*nabla)[Dim + 1] = (void *) (phgGeomGetJacobian(g, e));

  if (no1 <= 0)
    no1 = NEdge;

      for (k = no0; k < no1; k++)
	{
	  v = values[k];
	  GetEdgeVertices (e, k, I, J);
	  v[0] = lambda[I] * nabla[J][0] - lambda[J] * nabla[I][0];
	  v[1] = lambda[I] * nabla[J][1] - lambda[J] * nabla[I][1];
	  v[2] = lambda[I] * nabla[J][2] - lambda[J] * nabla[I][2];
	}

  return (FLOAT *)(values + no0);
}


static const FLOAT *
ND1_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of gradients of basis functions 
 * | \nabla\lambda1	-\nabla\lambda0		0		0      |
 * | \nabla\lambda2		0	-\nabla\lambda0		0      |
 * | \nabla\lambda3		0		0	-\nabla\lambda0|
 * |       0		\nabla\lambda2  -\nabla\lambda1		0      |	
 * |	   0		\nabla\lambda3     	0	-\nabla\lambda1|
 * |	   0			0	\nabla\lambda3  -\nabla\lambda2|
 * */
{
  static FLOAT values[NEdge][Dim][Dim+1], *v[Dim];
#if USE_OMP
#pragma omp threadprivate(values, v)
#endif  /* USE_OMP */
  int I, J, k;
  GRID *g = dof->g;
  FLOAT (*nabla)[Dim + 1] = (void *)(phgGeomGetJacobian(g, e));

  if (no1 <= 0)
    no1 = NEdge;

  bzero(values, sizeof(values));

      for (k = no0; k < no1; k++)
	{
	  v[0] = values[k][0];
	  v[1] = values[k][1];
	  v[2] = values[k][2];
	  GetEdgeVertices (e, k, I, J);
	  v[0][I] = nabla[J][0];
	  v[0][J] -= nabla[I][0];

	  v[1][I] = nabla[J][1];
	  v[1][J] -= nabla[I][1];

	  v[2][I] = nabla[J][2];
	  v[2][J] -= nabla[I][2];
	}

  return (FLOAT *)(values + no0);
}

#if 0
static FLOAT *
ND1_get_u(DOF *dof, ELEMENT *e, FLOAT *lambda, FLOAT *old_data)
{
  int i;
  FLOAT *basv;
  SHORT nbas = dof->type->nbas;
  static FLOAT values[] = { 0, 0, 0 };
#if USE_OMP
#pragma omp threadprivate(values)
#endif  /* USE_OMP */

  basv = Nd1_bas(dof->g, e, -1, lambda);
  for (i = 0; i < nbas; i++)
    {
      values[0] += old_data[i] * (*basv++);
      values[1] += old_data[i] * (*basv++);
      values[2] += old_data[i] * (*basv++);
    }
  return values;
}
#endif

static void
ND1_interp(DOF *dof, ELEMENT *e, FLOAT **old_data, FLOAT **new_data)
/* interpolation  
 * suppose the direction of the new edge is : oldvertex -> newvertex
 * just suit for dof->dim=1 
 * */
{
  INT  V0, V1,V2,V3,V4;
  int i,dim=dof->dim;
  FLOAT u0, u1;
  GRID *g = dof->g;
   
  V0 = GlobalVertexP(g, e->verts[0]);
  V1 = GlobalVertexP(g, e->verts[1]);
  V2 = GlobalVertexP(g, e->verts[2]);
  V3 = GlobalVertexP(g, e->verts[3]);
  V4 = GlobalVertexP(g, e->children[0]->verts[3]);
  for (i = 0; i < dim; i++)
   {
    u0 = 0.5 * old_data[4][i];
    new_data[1][i] = ((V1 - V0 > 0) == (V4 - V0 > 0)) ? u0 : -u0;
    new_data[2][i] = ((V1 - V0 > 0) == (V1 - V4 > 0)) ? u0 : -u0;

    u0 = ((V4 - V2 > 0) == (V0 - V2 > 0)) ? old_data[5][i] : - old_data[5][i];
    u1 = ((V4 - V2 > 0) == (V1 - V2 > 0)) ? old_data[7][i] : - old_data[7][i];
    new_data[3][i] = 0.5 * (u0 + u1);

    u0 = ((V4 - V3 > 0) == (V0 - V3 > 0)) ? old_data[6][i] : - old_data[6][i];
    u1 = ((V4 - V3 > 0) == (V1 - V3 > 0)) ? old_data[8][i] : - old_data[8][i];
    new_data[4][i] = 0.5 * (u0 + u1);
   }
}

static void
ND1_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	  DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	  const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
{
  GRID *g = dof->g;
  int i, v0, v1, u0, u1;
  FLOAT x0, y0, z0, x1, y1, z1;
  COORD *c;
  FLOAT *tmp = NULL;

  Unused(pdofvalues);

  assert (type == EDGE);
  assert (funcvalues == NULL || (userfunc == NULL && userfunc_lambda == NULL));
  v0 = e->verts[u0 = GetEdgeVertex(index, 0)];
  v1 = e->verts[u1 = GetEdgeVertex(index, 1)];
  if (GlobalVertexP(g, v0) > GlobalVertexP(g, v1))
    {
      i = v0;
      v0 = v1;
      v1 = i;
    }

  x0 = (*(c = g->verts + v0))[0];
  y0 = (*c)[1];
  z0 = (*c)[2];

  x1 = (*(c = g->verts + v1))[0];
  y1 = (*c)[1];
  z1 = (*c)[2];

  if (funcvalues == NULL)
    {
      funcvalues = tmp = phgAlloc (dof->dim * Dim * sizeof (*funcvalues));
      if (userfunc_lambda != NULL)
        {
	  FLOAT lambda[Dim + 1] = {0., 0., 0., 0.};
	  lambda[u0] = lambda[u1] = 0.5;
	  userfunc_lambda(dof, e, 0, lambda, tmp);
        }
      else
        {
          userfunc((x0 + x1) * .5, (y0 + y1) * .5, (z0 + z1) * .5, tmp);
        }
    }

  x1 -= x0;
  y1 -= y0;
  z1 -= z0;
  for (i = 0; i < dof->dim; i++)
    {
      *(dofvalues++) = x1 * funcvalues[0] + y1 * funcvalues[1]
	+ z1 * funcvalues[2];
      funcvalues += Dim;
    }

  phgFree(tmp);
}

#if 0
# define ND1_interp	phgDofInterpC2FGeneric
#endif

DOF_TYPE DOF_ND1_ = {
  DofReserved,
  "ND1",			/* name */
  NULL,				/* points */
  NULL,				/* orders */
  DOF_DG0,			/* gradient type */
  NULL,				/* base type */
  NULL,				/* hp_info */
  ND1_interp,			/* InterpC2F */
  phgDofInterpF2CGeneric,	/* InterpF2C */
  ND1_init,			/* initialization function */
  ND1_bas,			/* basis functions */
  ND1_grad,			/* gradient of basis functions */
  NULL,
  FE_Hcurl,
  FALSE,			/* is_nodal */
  FALSE,			/* invariant */
  FALSE,			/* free_after_use */
  -1,				/* id */
  NEdge,			/* nbas */
  1,				/* polynomial order */
  0,				/* D_order */
  -1,				/* continuity (discontinuous) */
  Dim,				/* dim */
  0,				/* np_vert */
  1,				/* np_edge */
  0,			/* np_face */
  0			/* np_elem */
};

DOF_TYPE *DOF_NDn[] = {DOF_ND1, NULL};
