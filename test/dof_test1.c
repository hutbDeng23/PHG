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

/* This code is for testing the ND1 basis.
 * $Id: dof_test1.c,v 1.65 2021/03/11 08:10:48 zlb Exp $ */

#include "phg.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void
print_basf(GRID *g,ELEMENT *e,DOF *u,FLOAT* lambda)
{
    const FLOAT *v;
    FLOAT (*J)[Dim + 1] = (void *)(phgGeomGetJacobian(u->g, e));
    int i;
    v=u->type->BasFuncs (u, e, 0, -1, lambda);
   phgPrintf("----------------------------------------\n");
   phgPrintf("----------------------------------------\n");
   phgPrintf("basis function values:\n");
    phgPrintf("e->index %d : Volume : %lg\n",e->index,(double)phgGeomGetVolume (g, e));
   for(i=0;i<u->type->nbas;i++,v+=3)
     phgPrintf("basf[%d]= %9lf , %9lf , %9lf\n",i,(double)v[0],(double)v[1],(double)v[2]);  
   phgPrintf("Jacobi:\n");
   for(i=0;i<Dim+1;i++)
     phgPrintf("%6.3lf %6.3lf %6.3lf %6.3lf\n",(double)J[i][0],(double)J[i][1],(double)J[i][2],(double)J[i][3]);
   phgPrintf("----------------------------------------\n");
  return; 
}

#if 0
#define ff(v,i,j)\
     (v[i][0]*v[j][0]+v[i][1]*v[j][1]+v[i][2]*v[j][2])
#define ee(e,v,i,j,val) { \
  int vi0,vi1,vj0,vj1;\
  GetEdgeVertices(e, i, vi0, vi1);  \
  GetEdgeVertices(e, j, vj0, vj1);   \
  val=                                                \
  (v[vi0][1]*v[vi1][2]-v[vi0][2]*v[vi1][1])*(v[vj0][1]*v[vj1][2]-v[vj0][2]*v[vj1][1]) +\
  (v[vi0][2]*v[vi1][0]-v[vi0][0]*v[vi1][2])*(v[vj0][2]*v[vj1][0]-v[vj0][0]*v[vj1][2]) +\
  (v[vi0][0]*v[vi1][1]-v[vi0][1]*v[vi1][0])*(v[vj0][0]*v[vj1][1]-v[vj0][1]*v[vj1][0]); \
} 

static void 
check_quad(DOF* u)
{
  GRID *g=u->g;
  ELEMENT *e;
  int i=0,j=0;
  int N=u->type->nbas* u->dim;
  FLOAT FQ[N][N],EQ[N][N],F[N][N],E[N][N];
  FLOAT vol,tmp;
  FLOAT (*J)[Dim + 1];

  ForAllElements(g,e){
    vol=phgGeomGetVolume(g, e);	   
    J=(void *) (phgGeomGetJacobian (g, e));
    F[0][0] = (ff(J,1,1)-ff(J,0,1)+ff(J,0,0))*vol/10.;
    F[0][1] = F[1][0] =(2*ff(J,1,2)-ff(J,1,0)-ff(J,0,2)+ff(J,0,0))*vol/20.; 
    F[0][2] = F[2][0] =(2*ff(J,1,3)-ff(J,1,0)-ff(J,0,3)+ff(J,0,0))*vol/20.;
    F[0][3] = F[3][0] =(ff(J,1,2)+ff(J,0,1)-ff(J,1,1)-2*ff(J,0,2))*vol/20.;
    F[0][4] = F[4][0] =(ff(J,1,3)-ff(J,1,1)+ff(J,0,1)-2*ff(J,0,3))*vol/20.;
    F[0][5] = F[5][0] =(ff(J,1,3)-ff(J,1,2)-ff(J,0,3)+ff(J,0,2))*vol/20.;

    F[1][1] = (ff(J,2,2)-ff(J,0,2)+ff(J,0,0))*vol/10.;
    F[1][2] = F[2][1] =(2*ff(J,2,3)-ff(J,0,2)-ff(J,0,3)+ff(J,0,0))*vol/20.;
    F[1][3] = F[3][1] =(ff(J,2,2)-ff(J,1,2)-ff(J,0,2)+2*ff(J,0,1))*vol/20.;
    F[1][4] = F[4][1] =(ff(J,2,3)-ff(J,1,2)+ff(J,0,1)-ff(J,0,3))*vol/20.;
    F[1][5] = F[5][1] =(ff(J,0,2)-ff(J,2,2)-2*ff(J,0,3)+ff(J,2,3))*vol/20.;

    F[2][2] = (ff(J,3,3)-ff(J,0,3)+ff(J,0,0))*vol/10.;
    F[2][3] = F[3][2] =(ff(J,2,3)-ff(J,1,3)-ff(J,0,2)+ff(J,0,1))*vol/20.;
    F[2][4] = F[4][2] =(ff(J,3,3)-ff(J,1,3)+2*ff(J,0,1)-ff(J,0,3))*vol/20.;
    F[2][5] = F[5][2] =(ff(J,3,3)-ff(J,2,3)-ff(J,0,3)+2*ff(J,0,2))*vol/20.;

    F[3][3] = (ff(J,2,2)-ff(J,1,2)+ff(J,1,1))*vol/10.;
    F[3][4] = F[4][3] =(2*ff(J,2,3)-ff(J,1,2)+ff(J,1,1)-ff(J,1,3))*vol/20.;
    F[3][5] = F[5][3] =(ff(J,2,3)-ff(J,2,2)-2*ff(J,1,3)+ff(J,1,2))*vol/20.;

    F[4][4] = (ff(J,1,1)-ff(J,1,3)+ff(J,3,3))*vol/10.;
    F[4][5] = F[5][4] =(2*ff(J,1,2)-ff(J,1,3)+ff(J,3,3)-ff(J,2,3))*vol/20.;

    F[5][5] =(ff(J,3,3)-ff(J,2,3)+ff(J,2,2))*vol/10.;

    for (i = 0; i < N; i++)
      {
	for (j = 0; j < i; j++)
	 {
	   EQ[i][j]=EQ[j][i]=phgQuadCurlBasDotCurlBas (e,u,i,u,j,-1);
	   FQ[i][j]=FQ[j][i]=phgQuadBasDotBas (e, u, i, u, j,  -1);
	   ee(e,J,i,j,tmp)
	   E[j][i]=E[i][j]=tmp*vol*4.0;
        if(Fabs(FQ[i][j]-F[i][j])>1e-10) 
            phgPrintf("***e->index %d *(%d,%d) FQ %lg F %lg \n",e->index,i,j,(double)FQ[i][j],(double)F[i][j]);
        if(Fabs(EQ[i][j]-E[i][j])>1e-10) 
            phgPrintf("***e->index %d *(%d,%d)  EQ %lg E %lg\n",e->index,i,j,(double)EQ[i][j],(double)E[i][j]);
        }
     }
  }//end of ForAllElements
}
#endif

/* test code for DOF_ND1 */
#define TESTFUN 2
static void
f (FLOAT x, FLOAT y, FLOAT z, FLOAT * values)
{
	switch(TESTFUN)
	{
         case 1:
             *(values++) = x * y * z * (x - 1.0) * (y - 1.0) * (1.0 - z) * 
		                       (x - .5) * (y - .5) * (z - .5);
             *(values++) = sin (2.0 * M_PI * x) * sin (2.0 * M_PI * y) * 
		           sin (2.0 * M_PI * z);
             *(values++) =
	        (1.0 - exp(x))*(M_E - exp(x ))*(M_E - exp (2.0 * x )) * 
                (1.0 - exp(y))*(M_E - exp(y ))*(M_E - exp (2.0 * y )) * 
                (1.0 - exp(z))*(M_E - exp(z ))*(M_E - exp (2.0 * z )) ;
	     break;
        case 2:
/* (1 2 3) \cross (x y z)  +  (4 5 6) */
            *(values++) = 2. * z - 3. * y + 4.;
            *(values++) = 3. * x - 1. * z + 5.;
            *(values++) = 1. * y - 2. * x + 6.;
	    break;
	case 3:
	     *(values++) = 2*exp(x+y+z);
	     *(values++) = -exp(x+y+z);
	     *(values++) = -exp(x+y+z);
	     break;
	case 4:
	     *(values++) = x*(y-1)*z;
	     *(values++) = sin(x)*sin(z)*y;
	     *(values++) = x*x+exp(y)+z;
	     break;
	case 5:
	    *(values++) = x+y+1;
	    *(values++) = z-2*y+2;
	    *(values++) = 2*x+z+3;
	     break;
	default:
	     phgError(1,"Invalid TEST CASE no!\n");
            break;
	}
  return ;
}

static void
grad_f (FLOAT x, FLOAT y, FLOAT z, FLOAT * values)
{
  switch(TESTFUN){
    case 2:
  	 *(values++) = 0.; *(values++) = -3.;  *(values++) = 2.;
         *(values++) = 3.; *(values++) = 0.;  *(values++) = -1.;
         *(values++) = -2.;  *(values++) = 1.;  *(values++) = 0.;
	 break;
    case 3:
   *(values++)=2*exp(x+y+z);*(values++)=2*exp(x+y+z);*(values++)=2*exp(x+y+z);
   *(values++)=- exp(x+y+z);*(values++)=- exp(x+y+z);*(values++)=- exp(x+y+z);
   *(values++)=- exp(x+y+z);*(values++)=- exp(x+y+z);*(values++)=- exp(x+y+z);
	 break;
    case 4:
  	 *(values++) = (y-1)*z; *(values++) = x*z;  *(values++) =x*(y-1);
         *(values++)=cos(x)*y*sin(z);
	 *(values++)=sin(x)*sin(z);
	 *(values++)=sin(x)*y*cos(z);
         *(values++) = 2*x;  *(values++) = exp(y);  *(values++) = 1.;
	 break;
    case 5:
  	 *(values++) = 1.; *(values++) = 1.;  *(values++) = 0.;
         *(values++) = 0.; *(values++) = -2.;  *(values++) = 1.;
         *(values++) = 2.;  *(values++) = 0.;  *(values++) = 1.;
	 break;
    default:
	     phgError(1,"Invalid TEST CASE no!\n");
	 break;
  }
}

static void
div_f (FLOAT x, FLOAT y, FLOAT z, FLOAT * values)
{
  switch(TESTFUN)
  {
    case 2:
         *values = 0.;
	 break;
    case 3:
         *values = 0;
	 break;
    case 4:
         *values = (y-1)*z+sin(x)*sin(z);
	 break;
    case 5:
  	 *(values++) = 0.; 
	 break;
    default:
	 phgError(1,"Invalid TEST CASE no!\n");
	 break;
  }
  return;
}

static void
curl_f (FLOAT x, FLOAT y, FLOAT z, FLOAT * values)
{
  switch(TESTFUN)
  {
    case 2:
         *(values++) = 2.;  *(values++) = 4.;  *(values++) = 6.;
	 break;
    case 3:
         *(values++)=0.;*(values++)=3.*exp(x+y+z);*(values++)=-3.*exp(x+y+z);
	 break;
    case 4:
         *(values++) = exp(y)-sin(x)*cos(z)*y;
       	 *(values++) = x*(y-1.)-2*x;
       	 *(values++) = cos(x)*sin(z)*y-x*z;
	 break;
    case 5:
  	 *(values++) = -1.; *(values++) = -2.;  *(values++) = -1.;
	 break;
    default:
	 phgError(1,"Invalid TEST CASE no!\n");
	 break;
  }
  return;
}

void 
get_u(DOF * dof, ELEMENT * e, FLOAT * lambda, FLOAT* values)
{
	int i;
	FLOAT (*basv)[Dim]=(void *)dof->type->BasFuncs(dof, e, 0, -1, lambda);
	SHORT nbas = dof->type->nbas;
	FLOAT *data;
FLOAT x[NVert], y[NVert], z[NVert];
FLOAT nabla[NVert][Dim];
FLOAT vol=phgGeomGetVolume(dof->g,e)*6;
//FLOAT (*nabla0)[Dim + 1] =
//	        (void *) (phgGeomGetJacobian (dof->g, e));
int j;
int I,J;
FLOAT basv0[NEdge][Dim],*v;
COORD *c;
int v0;
for (i = 0; i < NVert; i++)
{
  v0 = e->verts[i];
  x[i] = (*(c = dof->g->verts + v0))[0];
  y[i] = (*c)[1];
  z[i] = (*c)[2];
}
 
nabla[0][0]=(-(y[2]-y[1])*(z[3] - z[1]) + (y[3]-y[1])*(z[2] - z[1]))/(vol);
nabla[0][1]=((x[2]-x[1])*(z[3] - z[1]) - (x[3]-x[1])*(z[2] - z[1]))/(vol);
nabla[0][2]=(-(x[2]-x[1])*(y[3] - y[1]) + (x[3]-x[1])*(y[2] - y[1]))/(vol);
	
nabla[1][0]=((y[2]-y[0])*(z[3] - z[0]) - (y[3]-y[0])*(z[2] - z[0]))/(vol);
nabla[1][1]=(-(x[2]-x[0])*(z[3] - z[0]) + (x[3]-x[0])*(z[2] - z[0]))/(vol);
nabla[1][2]=((x[2]-x[0])*(y[3] - y[0]) - (x[3]-x[0])*(y[2] - y[0]))/(vol);

nabla[2][0]=(-(y[1]-y[0])*(z[3] - z[0]) + (y[3]-y[0])*(z[1] - z[0]))/(vol);
nabla[2][1]=((x[1]-x[0])*(z[3] - z[0]) - (x[3]-x[0])*(z[1] - z[0]))/(vol);
nabla[2][2]=(-(x[1]-x[0])*(y[3] - y[0]) + (x[3]-x[0])*(y[1] - y[0]))/(vol);

nabla[3][0]=((y[1]-y[0])*(z[2] - z[0]) - (y[2]-y[0])*(z[1] - z[0]))/(vol);
nabla[3][1]=(-(x[1]-x[0])*(z[2] - z[0]) + (x[2]-x[0])*(z[1] - z[0]))/(vol);
nabla[3][2]=((x[1]-x[0])*(y[2] - y[0]) - (x[2]-x[0])*(y[1] - y[0]))/(vol);

for(j = 0; j < NEdge; j++)
{ 
v = basv0[j];
GetEdgeVertices (e, j, I, J);
if (GlobalVertex (dof->g, e->verts[I]) > GlobalVertex (dof->g, e->verts[J]))
{
i = I;
I = J;
J = i;
}
v[0] = lambda[I] * nabla[J][0] - lambda[J] * nabla[I][0];
v[1] = lambda[I] * nabla[J][1] - lambda[J] * nabla[I][1];
v[2] = lambda[I] * nabla[J][2] - lambda[J] * nabla[I][2];
}
	
        values[0]=0;
        values[1]=0;
        values[2]=0;
	for (i = 0; i < nbas; i++)
         {
	    data=DofEdgeData(dof, e->edges[i]);
            values[0] += (*data) * basv[i][0];
            values[1] += (*data) * basv[i][1];
            values[2] += (*data) * basv[i][2];
         }
	return;	 
}

void 
check_basfunc(GRID *g,DOF *u,ELEMENT *e)
       /* Check bas function, s.t.:
        * \int_e_ij bas_ij * \tau ds = 1 
        * */
{
        FLOAT x0, y0, z0, x1, y1, z1,x2, y2, z2;
	const FLOAT *basv;
        FLOAT quad_pt[][4]={{0.5,.5,0.,0.},
	               {.5,0.,0.5,0.},
		       {0.5,0.,0.,0.5},
		       {0.,0.5,0.5,0.},
		       {0.,0.5,0.,0.5},
		       {0.,0.,0.5,0.5}
                       };
        int v0,v1;
        int i,j;
        COORD *c;
        for(i=0;i<NEdge;i++)
         {
           v0 = GetEdgeVertex (i, 0);
           v1 = GetEdgeVertex (i, 1);
           v0 = e->verts[v0];
           v1 = e->verts[v1];
           if (GlobalVertex (g, v0) > GlobalVertex (g, v1))
            {
             j = v0;
             v0 = v1;
             v1 = j;
             }

           x0 = (*(c = g->verts + v0))[0];
           y0 = (*c)[1];
           z0 = (*c)[2];

           x1 = (*(c = g->verts + v1))[0];
           y1 = (*c)[1];
           z1 = (*c)[2];

           x2 = x1-x0;
           y2 = y1-y0;
           z2 = z1-z0;
      
           basv=u->type->BasFuncs(u,e,i, i+1, quad_pt[i]);
           if(Fabs(basv[0]*x2+basv[1]*y2+basv[2]*z2 -1.0)>1e-10)
             phgError(1,"%s:%d:  bas[%d(%d %d)] * e = %lg |  %lg %lg %lg |  %lg %lg %lg|  %lg %lg %lg|  %lg %lg %lg\n", __FILE__, __LINE__, i, v0, v1,
			(double)(basv[0] * x2 + basv[1] * y2 + basv[2] * z2),
			(double)basv[0], (double)basv[1], (double)basv[2],
			(double)x0, (double)y0, (double)z0,
			(double)x1, (double)y1, (double)z1,
			(double)x2, (double)y2, (double)z2);
        }
} 


FLOAT 
compute_HCurlerror(DOF *u,QUAD *quad,DOF_USER_FUNC userfunc,DOF_USER_FUNC usercurl)
{
      int i,j,nvalues=DofDim(u);
      FLOAT *lambda,*lambda0,*w;
      FLOAT *v1,*v2,*v3,*v4,d,d0,err_L2=0.;
      FLOAT values1[nvalues],values2[nvalues],xyz[nvalues];
      FLOAT values3[nvalues],values4[nvalues];
      GRID *g=u->g;
      ELEMENT *e;
#if USE_MPI
      FLOAT a, b;
#endif     
      if (quad == NULL) {
        i = u->type->order;
        quad = phgQuadGetQuad3D(i);
      }
      lambda0 = quad->points;
      ForAllElements (g, e)
      {
        lambda=lambda0; 
        w = quad->weights;
        d=0.;
        for(i=0;i<quad->npoints; i++) {
          phgDofEval(u, e, lambda, v1 = values1);
          phgDofEvalCurl(u, e, lambda, NULL,v3 = values3);
          phgGeomLambda2XYZ(g,e,lambda,xyz,xyz+1,xyz+2);
          userfunc ( xyz[0], xyz[1], xyz[2], v2=values2);
          usercurl ( xyz[0], xyz[1], xyz[2], v4=values4);
          d0 = 0.;
          for(j=0;j<nvalues;j++,v1++,v2++,v3++,v4++)
            d0+=((*v1)-(*v2))* ((*v1)-(*v2))+((*v3)-(*v4))* ((*v3)-(*v4));	
          d+=d0 * (*(w++));
          lambda += Dim + 1;	  
	}
//      err_infinity=(err_infinity>d?err_infinity:d);
      err_L2+=d*phgGeomGetVolume (g, e);
      }
#if USE_MPI
      a = err_L2;
      MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_SUM, g->comm);
      err_L2 = Sqrt(b); 
      /*
      a=err_infinity;
      MPI_Allreduce(&a, &b, 1, PHG_MPI_FLOAT, PHG_MAX, g->comm);
      err_infinity = Sqrt(b);
      */
#else 
      err_L2 = Sqrt(err_L2);
#endif
      return err_L2 ; 
}

int
main (int argc, char **argv)
{
  GRID *g;
  DOF *u, *tmp, *div_u,*curl_u,*grad_u;
  int i,j;
  static char *fn = "cube.dat";
  SOLVER solver;

  FunctionEntry;
  phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&fn);
  phgVerbosity = 0;
  phgInit(&argc, &argv);
  phgPrintf ("Testing DOF operations of \"%s\"\n", DOF_ND1->name);
  g = phgNewGrid (-1);
  phgImport (g, fn, FALSE);
  for (i = 0; i < 2; i++)
    phgRefineAllElements (g, 1);
 /* 
  for (i = 0; i < 3; i++)
    {
      phgRefineRandomElements (g, "25%");
      phgRedistributeGrid (g);
    }
    */
  phgPartitionGrid (g);
  phgCheckConformity (g);

/* Check DofEval ,phgLambda2XYZ and Interpolation*/
while(TRUE)
{
      FLOAT err_L2,err_interp,err_grad,err_curl;
      FLOAT err;
      ELEMENT *e;

      PrintTime (0);
      u = phgDofNew (g, DOF_ND1, 1, "u", DofInterpolation);
      phgDofSetDataByFunction (u, f);
      phgPrintf ("%d DOF, %d elements, load imbalance: %lg\n",
			DofDataCountGlobal(u), g->nleaf_global,
			(double)g->lif);
      
      err_L2=0.;
	
      /* checking gradient */
      grad_u = phgDofGradient(u, NULL, NULL, "grad_u");
      tmp = phgDofNew (g, DOF_P0,  Dim * DofDim(u), "tmp", grad_f);
      phgDofAXPY(-1.0, grad_u, &tmp);
      err_grad = phgDofNormL2 (tmp);
      phgPrintf ("Gradient (%s) error    = %0.3le.\n", grad_u->type->name,
			(double)err_grad);
      phgDofFree (&tmp);

      /* checking divergence */
      div_u = phgDofDivergence(u, NULL, NULL, "div_u");
      tmp = phgDofNew(g, u->type->grad_type,
	DofDim(u) / (Dim * u->type->grad_type->dim), "tmp", div_f);
      phgDofAXPY(-1.0, div_u, &tmp);
      err_grad = phgDofNormL2 (tmp);
      phgDofFree (&tmp);
      phgPrintf ("Divergence (%s) error    = %0.3le.\n", div_u->type->name,
			(double)err_grad);

      /* checking curl */
      curl_u = phgDofCurl (u, NULL, NULL, "curl_u");
      tmp = phgDofNew (g, u->type->grad_type, 
	  DofDim(u) / u->type->grad_type->dim, "tmp", curl_f);
      phgDofAXPY(-1.0, curl_u, &tmp);
      err_curl = phgDofNormL2 (tmp);
      phgDofFree (&tmp);
      phgPrintf("Curl (%s) error = %0.3le.\n", u->type->name, (double)err_curl);
      
      /* checking FaceJump */
      tmp = phgQuadFaceJump(u, DOF_PROJ_DOT, "u . n", QUAD_DEFAULT);
      phgPrintf("Face jump (%s)      = %0.3le.\n",
			tmp->name, (double)phgDofNormInftyVec(tmp));
      phgDofFree(&tmp);

      tmp = phgQuadFaceJump(u, DOF_PROJ_CROSS, "u x n", QUAD_DEFAULT);
      phgPrintf("Face jump (%s)      = %0.3le.\n", tmp->name,
			(double)phgDofNormInftyVec(tmp));
      phgDofFree(&tmp);
      
      tmp = phgQuadFaceJump (curl_u, DOF_PROJ_CROSS, "curl_u x n", QUAD_DEFAULT);
      phgPrintf ("Face jump (%s) = %0.3le.\n", tmp->name,
			(double)phgDofNormInftyVec(tmp));
      phgDofFree (&tmp);

      
      /*Test BasFunc */
      ForAllElements (g, e)
      {
	check_basfunc(g,u,e);
      }//end of ForAllElement

      /*Test QUAD*/
      // check_quad(u);
      
      err_L2=compute_HCurlerror(u,phgQuadGetQuad3D(1),f,curl_f);
      phgPrintf("DofEvaluation L2 error\t=%0.3le.\n", (double)err_L2);
   // phgPrintf("DofEvaluation Infinity error\t=%0.3le.\n", (double)err_infinity);
      
      if (err_L2 < 1e-3 || g->nleaf_global / phgNProcs > 50000) 
      {
          phgDofFree(&div_u);
          phgDofFree(&grad_u);
          phgDofFree(&curl_u);
	  phgDofFree(&u);
	  break;
      }else {
         phgRefineAllElements (g, 3);
      }
      /* checking interpolation */
      tmp = phgDofNew (g, u->type, 1, "tmp", f);
      phgDofAXPY(-1.0, u, &tmp);
      err_interp = phgDofNormL2 (tmp);
      phgPrintf("Interpolation (%s) error\t=%0.3le.\n", 
		      u->type->name, (double)err_interp);
      phgDofFree (&tmp);
      
      /* testing DofCopy */
      tmp = phgDofCopy(u, NULL, DOF_P1, "tmp");
      phgDofFree(&u);
      u = phgDofNew(g, DOF_P1, Dim, "u0", f);
      phgDofAXPY(-1.0, u, &tmp);
      err = phgDofNormInftyVec(tmp);
      phgDofFree(&tmp);
      phgPrintf("Test phgDofCopy: %1.3le.\n",  (double)err);

      phgDofFree(&div_u);
      phgDofFree(&grad_u);
      phgDofFree(&curl_u);
      phgDofFree (&u);
      phgPrintf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
}//end of while 

  phgPrintf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
  phgPrintf("%lf %lf \n", DBL_MIN, DBL_EPSILON);
  phgPrintf("%d doube %d  int %d  struct %d  BOOLEAN %d BYTE %d  ",
sizeof(solver),sizeof(FLOAT),sizeof(INT),sizeof(solver.mat->rmap),sizeof(BOOLEAN),sizeof(BYTE));
  phgFreeGrid(&g);

  
if(FALSE)
{ 
  FLOAT b_dot_b,c_dot_c,d_dot_b,f_dot_b;	
  ELEMENT *e;
  g = phgNewGrid (-1);
  phgImport (g, fn, FALSE);

  for(i=0;i<1;i++)
  {	
    u = phgDofNew (g, DOF_ND1, 1, "u", DofInterpolation);
    phgDofSetDataByFunction (u, f);
    ForAllElements(g,e){ 
      b_dot_b=phgQuadBasDotBas (e, u, 1, u,1,-1);
      c_dot_c=phgQuadCurlBasDotCurlBas (e, u, 1, u,1,-1);
      d_dot_b=phgQuadDofDotBas (e, u, u,1,-1);
      f_dot_b=phgQuadFuncDotBas (e, f, u, 1, -1);
      phgPrintf("%d :  bas_dot_bas:   %lg\n     curl_dot_crul: %lg\n     dof_dot_bas:   %lg\n     fun_dot_bas:   %lg \n",e->index, (double)b_dot_b, (double)c_dot_c, (double)d_dot_b, (double)f_dot_b);

  {
     FLOAT lambda[]={0.25,0.25,0.25,0.25};   
     print_basf(g,e,u,lambda);
     phgPrintf("---------QuadBasDotBas----------------------\n");
     for(i=0;i<6;i++)
	for(j=i;j<6;j++)
	{
          b_dot_b=phgQuadBasDotBas (e, u, i, u,j,-1);
          phgPrintf("%d %d: %lg\n",i,j,(double)b_dot_b); 
	}
  }

    }//end for ForAllElements
     phgDofFree (&u);
     phgRefineAllElements (g,1);
  }//end of For-Block
  phgFreeGrid(&g);
}
  phgFinalize ();
  Return 0;
}
