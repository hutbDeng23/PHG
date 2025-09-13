/* Testing solving the following time-harmonic Maxwell equations using hp-FEM.
 *
 *	curl curl E + E = J
 *
 * Based on a code by Jiang XUE (jxue@lsec.cc.ac.cn)
 *
 * $Id: maxwell-hptest.c,v 1.3 2021/03/11 08:10:48 zlb Exp $
 */

#include "phg.h"
#include <math.h>
#include <stdlib.h>

static double
elapsed_time(GRID *g, BOOLEAN flag, double mflops)
/* returns and prints elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;
 
    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

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
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
  FLOAT r = Sqrt(x*x + y*y + z*z + 0.00000005);

  *(value++) = x / r;
  *(value++) = y / r;
  *(value++) = z / r;

    return;
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
  func_u(x, y, z, values);
  return;
}

static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
    GRID *g = u_h->g;
    SIMPLEX *e;
    int i, j, order;

    assert(u_h->dim == 1);

    ForAllElements(g, e) {
      int N = DofNBas(u_h, e); 
      FLOAT A[N][N], B[N], buffer[N];
      INT I[N];

	/* compute \int \curl\phi_i \cdot \curl\phi_j making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (j = 0; j <= i; j++)
		A[j][i] = A[i][j] =
		    phgQuadCurlBasDotCurlBas(e, u_h, i, u_h, j, QUAD_DEFAULT) +
		  phgQuadBasDotBas(e, u_h, i, u_h, j, QUAD_DEFAULT);
	}

	/* compute local matrix and RHS */
	for (i = 0; i < N; i++) {	/* loop on basis functions */
	    if (phgDofDirichletBC(u_h, e, i, func_u, 
				buffer, B + i, DOF_PROJ_CROSS)) {
	        phgPrintf("Error in Dirichlet Boundary Condition\n");
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer); 
	    }
	    else {
		/* right hand side: \int f * phi_i */
		order = DofElementOrder(u_h, e->index) * 2;
		B[i] = phgQuadDofDotBas(e, f_h, u_h, i, order);

		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
	    }
	}
	phgSolverAddRHSEntries(solver, N, I, B);
    }
}

static FLOAT
estimate_error(DOF *u_h, DOF *f_h, DOF *curl_u_h, DOF *error)
{
    GRID *g = u_h->g;
    SIMPLEX *e;
    DOF *jmp1, *jmp2, *res, *div, *tmp;
    FLOAT order, ak;
    NEIGHBOUR_DATA *aa;

    res = phgDofGetSameOrderDG(u_h, -1, NULL);
    phgDofCopy(f_h, &res, NULL, "f+k^2u_h");
    phgDofAXPY(-1., u_h, &res);			/* res = f_h + k^2 u_h */
    jmp1 = phgQuadFaceJump(curl_u_h, DOF_PROJ_CROSS, NULL, QUAD_DEFAULT);
    jmp2 = phgQuadFaceJump(res, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    tmp = phgDofCopy(res, NULL, NULL, NULL);	/* save f_h + k^2 u_h */
    div = phgDofCurl(curl_u_h, NULL, NULL, "residual");	/* div = curlcurl u_h */
    phgDofAXPBY(-1.0, div, 1.0, &res);		/* residual */
    phgDofFree(&div);
    div = phgDofDivergence(tmp, NULL, NULL, "div(f+k^2 u_h)");
    phgDofFree(&tmp);

    tmp = phgDofNew(g, DOF_P0, 1, "order of u_h", DofNoAction);

    ForAllElements(g, e) {
        *DofElementData(tmp, e->index) = DofElementOrder(u_h, e->index);
    }
    aa = phgDofInitNeighbourData(tmp, NULL);  

    ForAllElements(g, e) {
	int i;
	FLOAT eta, h;
	FLOAT diam = phgGeomGetDiameter(g, e);
	e->mark = 0;		/* clear refinement mmark */
	eta = 0.0;
	for (i = 0; i < NFace; i++) {
	    INT fno;
	    if (IS_BDRY(e->bound_type[i]))
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    fno = e->faces[i];

	    order = DofElementOrder(u_h, e->index);
            ak = *phgDofNeighbourData(aa, e, i, 0, NULL);
            order = (order > ak ? order : ak);

	    eta += (*DofFaceData(jmp1, fno) + *DofFaceData(jmp2, fno)) * h / order;
	}
	order = DofElementOrder(u_h, e->index);
	eta = eta * 0.5 + (phgQuadDofDotDof(e, res, res, order * 2) +
		phgQuadDofDotDof(e, div, div, order * 2)) * diam * diam / Pow(order, 2.);
	*DofElementData(error, e->index) = eta;
    }
    phgDofFree(&jmp1);
    phgDofFree(&jmp2);
    phgDofFree(&res);
    phgDofFree(&div);

    phgDofFree(&tmp);
    phgDofReleaseNeighbourData(&aa);
    return Sqrt(phgDofNormInftyVec(error));
}


/* calculate L2 norm of (u_h - u) */
static FLOAT
L2_error(DOF *u_h, DOF *u)
{
    SIMPLEX *e;
    FLOAT tmp, result;
    INT order;

    assert((u != NULL) && (u_h != NULL));
    assert(!SpecialDofType(u_h->type));

    tmp = 0.0;
    ForAllElements(u->g, e){
        order = DofElementOrder(u_h, e->index) * 2 + 3;
        tmp += phgQuadDofDotDof(e, u_h, u_h, order);
        tmp -= 2.0 * phgQuadDofDotDof(e, u_h, u, order);
        tmp += phgQuadDofDotDof(e, u, u, order);
    }
#if USE_MPI
    MPI_Allreduce(&tmp, &result, 1, PHG_MPI_FLOAT, PHG_SUM, u->g->comm);
#else
    result = tmp;
#endif

    return Sqrt(Fabs(result));
}

static void
hp_refine(DOF *u_h, DOF *error, DOF *errorp, DOF *vol, FLOAT theta,
		INT depth, INT max_order)
{
  GRID    *g = u_h->g;
  SIMPLEX *e;
  INT count, count1, count2;
  FLOAT max_err = Sqrt(phgDofNormInftyVec(error));
  DOF *inc_ord;
  INT p_max = 1, p_max2 = 1;
  INT p_min = 15, p_min2 = 15;

  inc_ord = phgDofNew(g, DOF_P0, 1, "increase order?", DofInterpolation);

  phgPrintf("determine p- and h-refinement:");
  phgMarkRefine(MARK_MAX, error, Pow(theta, 2.), NULL, 0.0,
		depth, max_err / g->nleaf_global);
  
  /* mark h- or p-refinement */
  count2 = count1 = 0;
  ForAllElements(g, e) {
    if (e->mark > 0) {
	/* if error <= predicted error */
	if (*DofElementData(error, e->index) 
	    <= *DofElementData(errorp, e->index)) {
	  /* h-refinement */
	  *DofElementData(inc_ord, e->index) = -1.;
	  count2++;
	}
	else {
	  if (DofElementOrder(u_h, e->index) + 1 > max_order) {
	    /* keep order unchanged and h-refinement */
	    *DofElementData(inc_ord, e->index) = -1.;
	    count2++;
	  }
	  else {
	    /* increase order p to p+1 and 
	     * remove h-refinement mark */
	    *DofElementData(inc_ord, e->index) = 1.;
	    e->mark = 0;
	    count1++;
	  }
	}
    } /* if mark > 0 */
    else {
      /* keep order unchaged for mark == 0*/
      *DofElementData(inc_ord, e->index) = -1.;
    }
  } /* ForAllElements */
  
#if USE_MPI
  MPI_Allreduce(&count2, &count, 1, MPI_INT, MPI_SUM, g->comm);
#else
  count = count2;
#endif
  phgPrintf(" [M,h: %d]  ", count);
  
#if USE_MPI
  MPI_Allreduce(&count1, &count, 1, MPI_INT, MPI_SUM, g->comm);
#else
  count = count1;
#endif
  phgPrintf(" [M,p: %d]  ", count);
  
  
  /* h-refinement */
  phgRefineMarkedElements(g);
  
  /* mark elements for p-refinement */
  count2 = 0;
  ForAllElements(g, e) {
    /* increase order or not */
    if (*DofElementData(inc_ord, e->index) > 0.) {
      e->hp_order = DofElementOrder(u_h, e->index) + 1;
      count2++;
    }
    else {
      e->hp_order = DofElementOrder(u_h, e->index);
    }
  }
  
#if USE_MPI
  MPI_Allreduce(&count2, &count, 1, MPI_INT, MPI_SUM, g->comm);
#else
  count = count2;
#endif
  phgPrintf("[R,p: %d]\n", count);
  
  
  /* p-refinement */
  phgHPSetup(u_h->hp, TRUE);

  ForAllElements(g, e) {
    p_max2 = (p_max2 > e->hp_order ? p_max2 : e->hp_order);
    p_min2 = (p_min2 < e->hp_order ? p_min2 : e->hp_order);
  }

#if USE_MPI
  MPI_Allreduce(&p_max2, &p_max, 1, MPI_INT, MPI_MAX, g->comm);
  MPI_Allreduce(&p_min2, &p_min, 1, MPI_INT, MPI_MAX, g->comm);
#else
  p_max = p_max2;
  p_min = p_min2;
#endif

  phgPrintf("\n---  Maximal order in NEXT step: %d \n",  p_max);
  phgPrintf("---  Minimal order in NEXT step: %d \n",  p_min);

  {
    FLOAT prt, vlm;
    INT kk;
    
    ForAllElements(g, e) {
      
      kk =  DofElementOrder(u_h, e->index);
      vlm = phgGeomGetVolume(g, e);
      prt = vlm / *DofElementData(vol, e->index);

      if (*DofElementData(inc_ord, e->index) > 0.) 
	prt = Pow(prt, (kk-1) * 2. / 3.) * Pow((kk-1) / kk, (kk-1));
      else
	prt = Pow(prt, kk * 2. / 3.);

      prt *= 1.0;
      //phgPrintf("---***  prt = %0.6le\n", prt);

      *DofElementData(vol, e->index) = vlm;
      *DofElementData(errorp, e->index) = prt * 
	*DofElementData(error, e->index);
    }
  }
  phgDofFree(&inc_ord);

  return;
}
 
int
main(int argc, char *argv[])
{
    GRID *g;
    SIMPLEX *e;
    HP_TYPE *hp;
    DOF *u_h, *u, *f_h, *error, *curl_u_h;
    SOLVER *solver, *pc = NULL;
    FLOAT L2_err, a, b, max_error;
    FLOAT e_L2, n_L2, n_Hc;
    FLOAT tmp1, tmp2;
    size_t mem, mem_peak;

    DOF *errorp;           /* predicted error */
    DOF *vol;             /* diameter of element */
    DOF *inc_ord;          /* if increase order or not */

    DOF *geneo;            /* generation in last step */
    DOF *hpord;            /* element order in last step */
    INT p_max = 1, p_max2 = 1;
    INT p_min = 15, p_min2 = 15;

    /*char *fn = "cube1o.dat";*/
    char *fn = "NeumannFichera.dat";
    char *vtk = NULL;
    char *pc_opts = NULL;
    INT pre_refines = 3;
    INT depth = 1;			/* <=0 => uniform refine -depth or 3 */
    INT mem_max = 1024;			/* max. memory */
    FLOAT tol = 1e-10;			/* convergence tolerance */
    int max_order = 8;
    INT submesh_threshold = 1000;
    FLOAT lif_threshold = 1.2;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", &fn);
    phgOptionsRegisterFilename("-vtk_file", "VTK file", &vtk);
    phgOptionsRegisterString("-pc_opts", "Preconditioner options", &pc_opts);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-refine_depth", "Refinement depth (<=0 means "
			  "uniformly refine -depth or 3 times)", &depth);
    phgOptionsRegisterInt("-mem_max", "Max memory (MB)", &mem_max);
    phgOptionsRegisterFloat("-tol", "Convergence tolerance", &tol);
    phgOptionsRegisterInt("-submesh_threshold", "Submesh threshold",
				&submesh_threshold);
    phgOptionsRegisterFloat("-lif_threshold", "LIF threshold", &lif_threshold);

    phgOptionsPreset("-dof_type HC1");
    phgOptionsPreset("-strategy max");

    phgInit(&argc, &argv);
 
    assert((max_order >= 1) && (max_order <= 15));
    assert(DOF_DEFAULT->order >= 1);

    phgOptionsShowUsed();

    if (depth == 0)
	depth = -3;

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    phgRefineAllElements(g, pre_refines);
    /* init variables */
    hp = phgHPNew(g, HP_HC);
    ForAllElements(g, e) {
      e->hp_order = DOF_DEFAULT->order;
    }
    phgHPSetup(hp, TRUE);
    u_h = phgHPDofNew(g, hp, 1, "u_h", DofInterpolation);
    phgHPFree(&hp);
 
    //u_h->DB_mask = BDRY_MASK;

    f_h = phgDofNew(g, DOF_ANALYTIC, 3, "f_h", func_f);
    u   = phgDofNew(g, DOF_ANALYTIC, 3, "u", func_u);

    error = phgDofNew(g, DOF_P0, 1, "error estimates", DofNoAction);

    /* useful information for hp-adaptivity */
    geneo = phgDofNew(g, DOF_P0, 1, "generation of element2", DofInterpolation);
    vol = phgDofNew(g, DOF_P0, 1, "volume of element", DofInterpolation);
    inc_ord = phgDofNew(g, DOF_P0, 1, "increase order?", DofInterpolation);
    hpord = phgDofNew(g, DOF_P0, 1, "element order", DofInterpolation);
    errorp = phgDofNew(g, DOF_P0, 1, "errorp indicator", DofInterpolation);

    ForAllElements(g, e) {
        *DofElementData(geneo, e->index) = e->generation + 0.5;

        /* set errorp -1., force h-refinment at first step */
        *DofElementData(errorp, e->index) = -1.;
        *DofElementData(vol, e->index) = phgGeomGetVolume(g, e);
        *DofElementData(hpord, e->index) = DofElementOrder(u_h, e->index) + 0.5;
    }
 

    int flag = 0;
    while (TRUE) {
	elapsed_time(g, FALSE, 0.);
        flag++;
        phgPrintf("\n******** Pass %d \n", flag);
        phgPrintf("-------- %d DOF, %d elements, %d submeshes,load imbalance: %lg\n",
		  DofDataCountGlobal(u_h), g->nleaf_global, g->nprocs,
                (double)g->lif);
    
	if (phgBalanceGrid(g, lif_threshold, submesh_threshold, NULL, 0.)) {
	    phgPrintf("------ Repartition mesh: nprocs = %d, LIF = %lg ",
			    g->nprocs, (double)g->lif);
	    elapsed_time(g, TRUE, 0.);
	}
	
	phgPrintf("Set up linear solver: ");
	solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
	phgPrintf("solver LIF = %lg ", (double)solver->mat->rmap->lif);
	elapsed_time(g, TRUE, 0.);

	phgPrintf("Build linear system ");
	build_linear_system(solver, u_h, f_h);
	elapsed_time(g, TRUE, 0.);

	phgOptionsPush();
	phgOptionsSetOptions("-solver ams -solver_maxit 1 -verbosity 0");
	phgOptionsSetOptions("-solver_rtol 0 -solver_btol 0 -solver_atol 0");
	phgOptionsSetOptions("+solver_warn_maxit");
	phgOptionsSetOptions(pc_opts);
	pc = phgMat2Solver(SOLVER_DEFAULT, solver->mat);
	phgOptionsPop();
	if (pc->oem_solver == SOLVER_AMS)
	    phgSolverAMSSetConstantCoefficients(pc, 1.0, 1.0);
	phgSolverSetPC(solver, pc, NULL);

	phgPrintf("Solve linear system: ");
	phgSolverSolve(solver, TRUE, u_h, NULL);
	phgPrintf("nits=%d, resid=%0.4lg ", solver->nits,
			(double)solver->residual);
	elapsed_time(g, TRUE, 0.);
	phgSolverDestroy(&pc);
	phgSolverDestroy(&solver);


	e_L2 = L2_error(u_h, u);
	n_L2 = phgDofNormL2(u_h);
	phgPrintf(" L2(%d)= %.8le ;  RL2(%d)= %.8le ;\n", flag, e_L2,	\
		  flag, e_L2 / n_L2);

	curl_u_h = phgDofCurl(u_h, NULL, NULL, NULL);
	n_Hc = phgDofNormL2(curl_u_h);
	tmp2 = Sqrt(Pow(e_L2, 2.) + Pow(n_Hc, 2.));
	tmp1 = Sqrt(Pow(n_L2, 2.) + Pow(n_Hc, 2.));
	phgPrintf(" absHcurlSemi(%d)= %.8le ; absHcurl(%d)= %.8le ; \n", flag, n_Hc, flag, tmp2);

	max_error = estimate_error(u_h, f_h, curl_u_h, error);
	L2_err = Sqrt(phgDofNormL1Vec(error));

	phgPrintf(" DOF(%d)= %d ; RHcurl(%d)= %.8le ; Est(%d) = %0.8le ;\n", 
		  flag, DofDataCountGlobal(u_h), flag, tmp2 / tmp1, flag, (double)L2_err);
	phgDofFree(&curl_u_h);


	ForAllElements(g, e) {
	  p_max2 = (p_max2 > e->hp_order ? p_max2 : e->hp_order);
	  p_min2 = (p_min2 < e->hp_order ? p_min2 : e->hp_order);
	}	
#if USE_MPI
	MPI_Allreduce(&p_max2, &p_max, 1, MPI_INT, MPI_MAX, g->comm);
	MPI_Allreduce(&p_min2, &p_min, 1, MPI_INT, MPI_MIN, g->comm);
#else
	p_max = p_max2;
	p_min = p_min2;
#endif	
	phgPrintf("\n---  Maximal order of element this step: %d \n",  p_max);
	phgPrintf("---  Minimal order of element this step: %d \n",  p_min);
	elapsed_time(g, TRUE, 0.);

	mem = phgMemoryUsage(g, &mem_peak);
	a = mem / (1024.0 * 1024.0);
	b = mem_peak / (1024.0 * 1024.0);
	phgPrintf("Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		(double)a, (double)b);
	if (L2_err < tol || mem_peak >= 1024 * (size_t)mem_max * 1024) {
	    break;
	}

	phgPrintf("Refine mesh ");
#if 1
	//max_order = 5;
	phgPrintf("\n The depth of h refinement: %d\n", depth);
	phgPrintf(" The max order of p refinement: %d\n", max_order);
	hp_refine(u_h, error, errorp, vol, 0.1, depth, max_order);
	Unused(max_error);
#else
	int count, count2, count1;
	FLOAT theta = 0.1;
        /* p- and h-refinement */
        phgPrintf("determine p- and h-refinement:");
        phgMarkRefine(MARK_DEFAULT, error, theta, NULL,
                Pow(0.2, 2.), 1, max_error / g->nleaf_global);

        /* mark h- or p-refinement */
        count2 = count1 = 0;
        ForAllElements(g, e) {
            if (e->mark > 0) {
                /* h-refinement or p-refinement last step */
                if ((e->generation > *DofElementData(geneo, e->index))
                        || (DofElementOrder(u_h, e->index)
                            > *DofElementData(hpord, e->index))) {
                    /* if error <= predicted error */
                    if (*DofElementData(error, e->index) 
                            <= *DofElementData(errorp, e->index)) {
                        /* p-refinement */
                        /* hp_order can't greater than max_order, or
                         * use h-refinement instead of p-refinement */
                        if (DofElementOrder(u_h, e->index) + 1 > max_order) {
                            /* keep order unchanged and h-refinement */
                            *DofElementData(inc_ord, e->index) = -1.;
                            count2++;
                        }
                        else {
                            /* increase order p to p+1 and 
                             * remove h-refinement mark */
                            *DofElementData(inc_ord, e->index) = 1.;
                            e->mark = 0;
                            count1++;
                        }
                    }
                    else {
                        /* h-refinement */
                        *DofElementData(inc_ord, e->index) = -1.;
                        count2++;
                    }
                }
                else {
#if 0
                    /* p-refinement */
                    /* hp_order can't greater than max_order */
                    if (DofElementOrder(u_h, e->index) + 1 > max_order) {
                        /* keep order unchanged and h-refinement */
                        *DofElementData(inc_ord, e->index) = -1.;
                        count2++;
                    }
                    else {
                        /* increase order p to p+1 and 
                         * remove h-refinement mark */
                        *DofElementData(inc_ord, e->index) = 1.;
                        e->mark = 0;
                        count1++;
                    }
#else
                    if (*DofElementData(error, e->index) <=
                            *DofElementData(errorp, e->index)) {
                        if (DofElementOrder(u_h, e->index) + 1 > max_order) {
                            /* keep order unchanged and h-refinement */
                            *DofElementData(inc_ord, e->index) = -1.;
                            count2++;
                        }
                        else {
                            /* increase order p to p+1 and 
                             * remove h-refinement mark */
                            *DofElementData(inc_ord, e->index) = 1.;
                            e->mark = 0;
                            count1++;
                        }
                    }
                    else {
                        *DofElementData(inc_ord, e->index) = -1.;
                        count2++;
                    }
#endif
                }
            } /* if mark > 0 */
            else {
                /* keep order unchaged for mark == 0*/
                *DofElementData(inc_ord, e->index) = -1.;
            }
        } /* ForAllElements */

#if USE_MPI
        MPI_Allreduce(&count2, &count, 1, MPI_INT, MPI_SUM, g->comm);
#else
        count = count2;
#endif
        phgPrintf(" [M,h: %d]  ", count);

#if USE_MPI
        MPI_Allreduce(&count1, &count, 1, MPI_INT, MPI_SUM, g->comm);
#else
        count = count1;
#endif
        phgPrintf(" [M,p: %d]  ", count);
        /* save generation information for geneo 
         * We can do this in marking step 
         * but for the sake of clearness we do here */
        ForAllElements(g, e) {
            *DofElementData(geneo, e->index) = e->generation + 0.5;
        }

        /* h-refinement */
        phgRefineMarkedElements(g);

        /* mark elements for p-refinement */
        count2 = 0;
        ForAllElements(g, e) {
            /* increase order or not */
            if (*DofElementData(inc_ord, e->index) > 0.) {
                e->hp_order = DofElementOrder(u_h, e->index) + 1;
                count2++;
            }
            else {
                e->hp_order = DofElementOrder(u_h, e->index);
            }
        }

#if USE_MPI
        MPI_Allreduce(&count2, &count, 1, MPI_INT, MPI_SUM, g->comm);
#else
        count = count2;
#endif
        phgPrintf("[R,p: %d]\n", count);

        /* save element order information for hprod
         * We can do this in marking step 
         * but for the sake of clearness we do here */
        ForAllElements(g, e) {
            *DofElementData(hpord, e->index) = 
                DofElementOrder(u_h, e->index) + 0.5;
        }

        /* p-refinement */
        phgHPSetup(u_h->hp, TRUE);

        /* update predicted error information 
         * update generation information 
         * update diameter information */
        {
            FLOAT prt, vlm;
            INT kk;

            ForAllElements(g, e) {
                if (e->generation > *DofElementData(geneo, e->index)) {
                    /* h-refinement*/
                    if (*DofElementData(inc_ord, e->index) < 0.) {
                        kk =  2 * DofElementOrder(u_h, e->index);
                        vlm = phgGeomGetVolume(g, e);
                        prt = vlm / *DofElementData(vol, e->index);
                        prt = Pow(prt, kk / 3.);

                        *DofElementData(vol, e->index) = vlm;
                        *DofElementData(errorp, e->index) = prt * 
                            *DofElementData(error, e->index);
                    }
                    else { /* h-refinement && p-refinement */
                        kk =  DofElementOrder(u_h, e->index) - 1;
                        vlm = phgGeomGetVolume(g, e);
                        prt = vlm / *DofElementData(vol, e->index);
                        prt = Pow(prt, kk * 2. / 3.) * \
			      Pow(kk / (kk+1), kk);

                        *DofElementData(vol, e->index) = vlm;
                        *DofElementData(errorp, e->index) = prt * 
                            *DofElementData(error, e->index);
                    }
                }
                else { /* p-refinement */
                    if (*DofElementData(inc_ord, e->index) > 0.) {
                        kk =  DofElementOrder(u_h, e->index) - 1;
                        vlm = phgGeomGetVolume(g, e);
                        prt = Pow(kk / (kk+1), kk);

                        *DofElementData(vol, e->index) = vlm;
                        *DofElementData(errorp, e->index) = prt * 
                            *DofElementData(error, e->index);
                    } /* keep unchanged */
                    else {
                        vlm = phgGeomGetVolume(g, e);

                        *DofElementData(vol, e->index) = vlm;
                        *DofElementData(errorp, e->index) =
                            *DofElementData(error, e->index);
                    }
                }/* if */
            }  /* ForAllElements */
        }
#endif

	elapsed_time(g, TRUE, 0.);
    }

    if (vtk != NULL)
	phgPrintf("Export mesh to \"%s\" ",
			phgExportVTK(g, vtk, u_h, error, NULL));

    phgDofFree(&u_h);
    phgDofFree(&error);
    phgDofFree(&f_h);
    phgDofFree(&u);

    phgDofFree(&errorp);
    phgDofFree(&vol);
    phgDofFree(&inc_ord);
    phgDofFree(&geneo);
    phgDofFree(&hpord);

    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
