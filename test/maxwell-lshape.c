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

/* $Id: maxwell-lshape.c,v 1.26 2020/12/06 09:03:08 zlb Exp $
 *
 * L-shaped domain: (-1,1)^3 - (0,1)x(0,-1)x(-1,1),
 * Dirichlet BC, analytic solution:
 *	\grad (r^{1/2}\sin(\phi/2)) = \grad \sqrt{(r-x)/2}
 *		( -\sqrt{r-x}/(2\sqrt{2}r), \sqrt{r+x}/(2\sqrt{2}r), 0)
 * where r = \sqrt{x^2+y^2}
 */

#include "phg.h"
#include "io.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define KK 	1.0
#define KK_pc	(-Fabs(KK))
#define SHIFT	/*1.1*/0.0;

static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
#if 1
    FLOAT r, rsqrt8;

    y += SHIFT;
    if ((r = Sqrt(x * x + y * y)) < 1e-15)
	r = 1e-15;
    rsqrt8 = 2. * M_SQRT2 * r;
    *(values++) = -Sqrt(r - x) / rsqrt8;
    *(values++) = Sqrt(r + x) / rsqrt8;
    *(values++) = 0.;
#else
    /* \grad (x^3 + y^3 + z^3) */
    /**(values++) = 3. * x * x;
    *(values++) = 3. * y * y;
    *(values++) = 3. * z * z;*/

    /* \grad sin(x)sin(y)sin(z) */
    *(values++) = cos(x) * sin(y) * sin(z);
    *(values++) = sin(x) * cos(y) * sin(z);
    *(values++) = sin(x) * sin(y) * cos(z);

    /* (1,2,3) x (x,y,z) */
    /**(values++) = 2. * z - 3. * y;
    *(values++) = 3. * x - 1. * z;
    *(values++) = 1. * y - 2. * x;*/
#endif
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
  func_u(x, y, z, values);
  values[0] *= -KK;
  values[1] *= -KK;
  values[2] *= -KK;
  return;
}

static void
func_mu(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
/* returns 1/\mu */
{
    *value = 1.0;
}

#if 0
static void 
custom_init(DOF *dof, FLOAT a)
{
  GRID *g = dof->g;
  ELEMENT *e;
  BYTE *flag;
  FLOAT x0, y0, r0, x1, y1, r1;
  COORD *c;
  INT no;
  int i, v0, v1;

  assert(dof->type == DOF_ND1 || dof->type == DOF_HC0);

  flag = phgCalloc(g->nedge, sizeof(*dof->data));
  bzero(flag, sizeof(flag));

  ForAllElements (g,e)
  {
	for(i=0;i<NEdge;i++)
	{
	  no=e->edges[i];	
	  if(flag[no]) continue;
	  flag[no]=TRUE;
	  GetEdgeVertices(e, i, v0, v1);
	  x0 = (*(c = g->verts + e->verts[v0]))[0];
	  y0 = (*c)[1] + SHIFT;
	  x1 = (*(c = g->verts + e->verts[v1]))[0];
	  y1 = (*c)[1] + SHIFT;
	  
	  r0 = Sqrt(x0 * x0 + y0 * y0);
	  r1 = Sqrt(x1 * x1 + y1 * y1);
	  *DofEdgeData(dof, no) = a * (Sqrt((r1-x1)*.5) - Sqrt((r0-x0)*.5));
	}
  }
  phgFree(flag);

  return;
}
#endif	/* 0 */

static double 
elapsed_time(GRID *g, BOOLEAN flag, double mflops)
/* returns elapsed time since last call to this function */
{
  static double t0=0.0;
  double et, tt[3];
  size_t mem;

  phgGetTime(tt);
  et=tt[2] - t0;
  t0=tt[2];

   mem = phgMemoryUsage(g, NULL);
      if (flag) {
       if (mflops <= 0)
           phgPrintf("[%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
       else
           phgPrintf("[%0.4lgMB %0.4lfs %0.4lgGF]\n", mem / (1024.0 * 1024.0),
                      et, mflops*1e-3);
      }

   return et;		
}

static void
build_linear_system(SOLVER *solver, DOF *mu, DOF *u_h, DOF *f_h, FLOAT kk)
{
  GRID *g = u_h->g;
  ELEMENT *e;
  int i, j;
  int N = u_h->type->nbas * u_h->dim;	/* size of local matrix */
  FLOAT A[N][N], buffer[N], rhs[N];
  INT I[N];

  ForAllElements(g, e) {	/* \PHGindex{ForAllElements} */
    /* compute \int \curl\phi_j \cdot \curl\phi_i making use of symmetry */
    for (i = 0; i < N; i++) {
	I[i] = phgSolverMapE2L(solver, 0, e, i);
	for (j = 0; j <= i; j++) {
	        A[j][i] = A[i][j] =
		  phgQuadCurlBasACurlBas(e, u_h, j, mu, u_h, i, QUAD_DEFAULT) -
		  kk * phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	}
    }

    /* compute local matrix and RHS */
    for (i = 0; i < N; i++) {	/* loop on basis functions */
	if (phgDofDirichletBC(u_h, e, i, kk == KK_pc ? NULL : func_u,
			      buffer, rhs + i, DOF_PROJ_CROSS)) {
	    if (u_h->type == DOF_ND1 || u_h->type == DOF_HC0) {
		int v0, v1;
		FLOAT x0, y0, r0, x1, y1, r1;
		COORD *c;
		GetEdgeVertices(e, i, v0, v1);
		x0 = (*(c = g->verts + e->verts[v0]))[0];
		y0 = (*c)[1] + SHIFT;
		x1 = (*(c = g->verts + e->verts[v1]))[0];
		y1 = (*c)[1] + SHIFT;
		r0 = Sqrt(x0 * x0 + y0 * y0);
		r1 = Sqrt(x1 * x1 + y1 * y1);
		rhs[i] = Sqrt((r1-x1)*.5) - Sqrt((r0-x0)*.5);
		phgSolverAddMatrixEntry(solver, I[i], I[i], 1.0);
	    }
	    else {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer);
	    }
	}
	else {
	    phgSolverAddMatrixEntries(solver, 1, I + i, N, I, &A[i][0]);
	    /* right hand side: \int f * phi_i */
	    rhs[i] = phgQuadDofDotBas(e, f_h, u_h, i, QUAD_DEFAULT);
	}
    }
    phgSolverAddRHSEntries(solver, N, I, rhs);
  }
}

static FLOAT
estimate_error(DOF *mu, DOF *u_h, DOF *f_h, DOF *error, FLOAT *curl)
{
  GRID *g = u_h->g;
  ELEMENT *e;
  DOF *jump1,*jump2, *res,*div,*curl_u_h, *tmp;

  tmp = phgDofCurl(u_h, NULL, NULL, NULL);
  *curl = phgDofNormL2(u_h);
  curl_u_h = phgDofMM(MAT_OP_N, MAT_OP_N, 1, 3, 1, 1.0, mu, 1, tmp, 0.0, NULL);
  phgDofFree(&tmp);
  tmp = phgDofCopy(f_h, NULL, u_h->type, "f+k^2u_h");
  phgDofAXPY(KK, u_h, &tmp);
  jump1 = phgQuadFaceJump(curl_u_h, DOF_PROJ_CROSS, NULL, QUAD_DEFAULT);
  jump2 = phgQuadFaceJump(tmp, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
  res = phgDofCurl(curl_u_h, NULL, NULL, "residual");
  phgDofAXPBY(1.0, tmp, -1.0, &res);
  div = phgDofDivergence(tmp, NULL, NULL, "div(f+k^2 u_h)");
  phgDofFree(&tmp);
  ForAllElements(g, e)
  {
    int i;
    FLOAT eta, h;
    FLOAT diam = phgGeomGetDiameter(g, e);
    eta = 0.0;
    for (i = 0; i < NFace; i++)
      {
	INT fno;
	if (e->bound_type[i] & BDRY_MASK)
	  continue;		/* boundary face */
	h = phgGeomGetFaceDiameter(g, e,i);
	fno=e->faces[i];
	eta += (*DofFaceData(jump1, fno) + *DofFaceData(jump2, fno)) * h ;
      }
    eta += diam * diam * (phgQuadDofDotDof(e, res, res, QUAD_DEFAULT) +
			  phgQuadDofDotDof(e, div, div, QUAD_DEFAULT));
   *DofElementData(error, e->index) = Sqrt(eta);
  }
  phgDofFree(&jump1);
  phgDofFree(&jump2);
  phgDofFree(&res);
  phgDofFree(&div);
  phgDofFree(&curl_u_h);
  return phgDofNormInftyVec(error);
}

static void
mark_refine(DOF *error, FLOAT thres, INT depth)
{
  GRID *g = error->g;
  ELEMENT *e;

  ForAllElements(g, e)
    e->mark = (*DofElementData(error, e->index) > thres ? depth : 0);
}

/*==================================================================*/

int
main(int argc, char *argv[])
{
  GRID *g;
  DOF *mu, *u_h, *f_h, *error, *tmp;
  SOLVER *solver;
  SOLVER *pc;
  FLOAT PE_l2=0.0, PE_infinity = 0.0, Norm_u_h, curl;
  size_t mem, mem_peak; 
  
  /* refine strategies */
  static const char *strategies[] = {"maximum", "equidist", NULL};
  typedef enum {MAXIMUM, EQUIDIST} STRATEGY;
  static int strategy = MAXIMUM;

  static char *fn = "lshape.dat"; 
  static char *vtk_fn = "maxwell-lshape.vtk";
  static INT pre_refines=2;


  static INT nmax = 500000;
  static INT vtk_limit = 0; 
  static FLOAT tol = 1e-3;    /* convergence tolerance */
  static BOOLEAN export_mesh = TRUE;
  static INT refine_depth = 1; /* <=0 => uniformly refine -depth or 3 times */
  static INT mem_max = 400;   /* max. memory */
  static INT pc_maxit = -1;
  static BOOLEAN enable_pc = TRUE;
  /* total time */
  int refine_counts=1;  
  double t0, tt[3];     
  BOOLEAN stop;

  phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
  phgOptionsRegisterFilename("-vtk_file", "VTK file", (char **)&vtk_fn);
  phgOptionsRegisterInt("-vtk_limit", "max # of element exported to VTK file", 
			&vtk_limit);
  phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
  phgOptionsRegisterInt("-mem_max", "Max-memory(MB)", &mem_max);
  phgOptionsRegisterInt("-nmax", "Maxi. number of unknowns/process", &nmax);
  phgOptionsRegisterInt("-pc_maxit", "Maxi. it number in pc (-1 ==> exact PC)", 
				&pc_maxit);
  phgOptionsRegisterInt("-refine_depth",
                        "Minimum refine depth (0 ==> uniform refines)",
                         &refine_depth);
  phgOptionsRegisterNoArg("-export_mesh", "Export mesh", &export_mesh);
  phgOptionsRegisterKeyword("-refine_strategy", "Adative strategy",
		                            strategies, &strategy);
  phgOptionsRegisterNoArg("-enable_pc", "Enable preconditioner", &enable_pc);

  if (SOLVER_HYPRE == NULL)
    phgError(1, "This program requires 2.0.0 or higher version of HYPRE.\n");

  phgOptionsPreset("-dof_type ND1");
  /* This is for the main solver */
  phgOptionsPreset("-solver pcg");
  /* This is for the preconditioner using HYPRE PCG/AMS */
  phgOptionsPreset("-hypre_solver pcg -hypre_pc ams");
  phgInit(&argc, &argv);	/* 3rd arg = SPMD flag\PHGindex{phgInit} */
  g = phgNewGrid(-1);
  if (!phgImport(g, fn, FALSE))
    phgError(1, "can't read file \"%s\".\n", fn);

  phgRefineAllElements(g, pre_refines);

  phgGetTime(tt);
  t0=tt[2];
      
  u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
  u_h->DB_mask = BDRY_MASK;
  phgDofSetDataByValue(u_h, 0.0);

  f_h = phgDofNew(g, DOF_ANALYTIC, 3, "f_h", func_f);

  error = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);

  mu = phgDofNew(g, DOF_P0, 1, "mu", DofInterpolation);
  phgDofSetDataByFunction(mu, func_mu);

  while (TRUE)
    {
      elapsed_time(g,FALSE,0);
      phgPrintf("\n%d DOF, %d elements, LIF: %lg\n",
		 DofDataCountGlobal(u_h), g->nleaf_global, (double)g->lif);
      if (phgBalanceGrid(g, 1.2, 1000, NULL, 0)){
	phgPrintf("  Repartition mesh, %d proc(s), LIF =  %lg ", 
			 g->nprocs, (double)g->lif);
	elapsed_time(g,TRUE,0); 
      }
      phgPrintf("  Set up linear solver: ");
      solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
      phgPrintf("solver LIF = %lg ",  (double)solver->mat->rmap->lif);
      elapsed_time(g,TRUE,0);

      phgPrintf("  Build linear system: ");
      build_linear_system(solver, mu, u_h, f_h, KK);
      elapsed_time(g,TRUE,0);

      if (enable_pc) {
	phgPrintf("  Build preconditioner: ");
	pc = phgSolverCreate(SOLVER_HYPRE, u_h, NULL);
	build_linear_system(pc, mu, u_h, f_h, KK_pc);
	/*phgSolverHypreAMSSetConstantPoisson(pc, 1.0, -KK_pc);*/
     
	if (pc_maxit > 0) {
          phgSolverSetMaxIt(pc, pc_maxit);
          pc->warn_maxit = FALSE;
        }

        phgSolverSetPC(solver, pc, NULL);
        elapsed_time(g,TRUE,0);
      }

      phgPrintf("  Solve linear system: ");
#if 1
      phgSolverSolve(solver, TRUE, u_h, NULL);
      phgPrintf("nits = %d, resid = %0.4lg ", solver->nits,
			(double)solver->residual);
      elapsed_time(g,TRUE,0);
#else
    {
      VEC *x = phgMapCreateVec(solver->mat->rmap, 1), *b;
      MAT *A = phgSolverGetMat(solver);
      DOF *dofs[1] = {u_h};
      phgSolverSolve(solver, FALSE, u_h, NULL);
      b = phgVecCopy(solver->rhs, NULL);
      phgMapDofToLocalData(solver->mat->rmap, 1, dofs, x->data);
      phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &b);
      phgPrintf("nits = %d, resid = %0.4lg ", solver->nits,
			(double)solver->residual);
      elapsed_time(g,TRUE,0);
      phgPrintf("Residual = %lg\n", (double)phgVecNorm2(b, 0, NULL));
      break;
    }
#endif
      phgSolverDestroy(&solver);
      if (enable_pc)
	phgSolverDestroy(&pc);
      
      phgPrintf("  Estimate: "); 
      Norm_u_h=phgDofNormL2(u_h);
      phgPrintf("|u_h|=%0.3le, ", (double)Norm_u_h);

      PE_infinity = estimate_error(mu, u_h, f_h, error, &curl);
      PE_l2 = phgDofNormL2Vec(error);
      phgPrintf("PEoo=%0.3le, PE_L2=%0.3le ",
			(double)PE_infinity, (double)PE_l2);
      elapsed_time(g, TRUE, 0);
      tmp = phgDofCopy(u_h, NULL, NULL, "error");
      phgDofAXPY(1. / KK, f_h, &tmp);
      phgPrintf("  L2 error: %le, curl u_h: %le\n", (double)phgDofNormL2(tmp),
		(double)curl);
      mem = phgMemoryUsage(g, &mem_peak);
      phgPrintf("  Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		mem / (1024.0 * 1024.0), mem_peak / (1024.0 * 1024.0));
      phgGetTime(tt);
      phgPrintf("  Adaptive step: %d, elapsed time: %lg\n",
			refine_counts, (double)(tt[2] - t0));

      stop = (PE_l2 < tol || g->nleaf_global / phgNProcs > nmax ||
		mem_peak >= 1024 * (size_t)1024 * mem_max);

      if(export_mesh &&
	 ((vtk_limit > 0 && g->nleaf_global > vtk_limit) || stop)) {
	 export_mesh = FALSE;
	 phgPrintf("  ===> \"%s\", ",
			phgExportVTK(g, vtk_fn, u_h, mu, tmp, error, NULL));
	 if (phgRank == 0) {
	    FILE *f;
	    long long size;
	    f = fopen(vtk_fn, "r");
	    fseek(f, 0L, SEEK_END);
	    size = ftell(f);
	    fclose(f);
	    phgPrintf("size = %lld ", size);
	 }
	 elapsed_time(g, TRUE, 0.);
      }

      phgDofFree(&tmp);

      if (stop)
         break;

      phgPrintf("  Refine mesh: ");
      elapsed_time(g,FALSE,0);
      if (refine_depth <= 0) {
	    /*  uniform refine */
        phgRefineAllElements(g, refine_depth == 0 ? 3 : -refine_depth);
      }else{
	 /* adaptive refine */
        switch((STRATEGY)strategy){
	   case  MAXIMUM:
                 mark_refine(error, 0.5*PE_infinity, refine_depth);
                 break;  
           case EQUIDIST:
                 mark_refine(error, tol / Sqrt((double) g->nleaf_global),
				 refine_depth);
		 break;
           default:
		 phgError(1, "Undefined refine strategy\n");
	}
      phgRefineMarkedElements(g);
      }
      elapsed_time(g,TRUE,0);
      refine_counts++;
    }

  phgDofFree(&mu);
  phgDofFree(&u_h);
  phgDofFree(&f_h);
  phgDofFree(&error);		

  phgFreeGrid(&g);
  phgFinalize();
  return 0;
}
