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

/* ************************************************************
 *   problem: *   u_{t} - \Delta{u} = func_f  (x, t) \in \Omega X (0, T)
 *   u = func_g (x, t) \in \partial\Omega X [0, T]
 *   u = func_u0       x \in \Omega  t = 0.0
 *
 *   The space-time adaptation algorithm is based on:
 *	Zhiming CHEN et. al,
 *	An adaptive finite element method with reliable and efficient error 
 *	control for linear parabolic problems, Math. Comp. 73 (2004), 1163-1197.
 *   
 * $Id: heat.c,v 1.59 2021/03/11 08:09:53 zlb Exp $
 **************************************************************/

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <math.h>
#include <stdlib.h>


/* these parameters can be modified */
static FLOAT htolerance = 0.001;
static FLOAT T = 1.0;           /* time domain: [0, T],T > 0 */
#define FLAGA   25              /* parameter used by func_f, func_g, > 0. */
#define FLAGB   0.7             /* parameter used by func_f, func_g, > 0. */

/* default time step */
/* can be modified */
static FLOAT stime = 0.1;

/* thest can't be modified */
static FLOAT tol_space;
static FLOAT tol_time;
static FLOAT crtime = 0.0;    /* current time. Modified by program */
static FLOAT ptime = 0.0;     /* previous time. */
static FLOAT SpaceTOL;        /* TOL / T, */
static FLOAT TimeTOL;         /* TOL / (2 * T) */
static FLOAT FTOL;            /* Sqrt(TOL) / (2 * T) */

static void     /* rhs: f */
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT tmp, etmp, ct, st;
    ct = Cos(8.0 * M_PI * crtime);
    st = Sin(8.0 * M_PI * crtime);
    tmp = (x - 0.5 - 0.4 * st) * (x - 0.5 - 0.4 * st)
        + (y - 0.5 - 0.4 * ct) * (y - 0.5 - 0.4 * ct)
        + (z - 1.0) * (z - 1.0);
    etmp = 1.0 / (FLAGA * tmp + FLAGB + 0.2);
    *value = ((6.4 * FLAGA * M_PI * ct * (x - 0.5 - 0.4 * st) 
             - 6.4 * FLAGA * M_PI * st * (y - 0.5 - 0.4 * ct)) * etmp * etmp
             - 4.0 * FLAGA * FLAGA * tmp * etmp * etmp * etmp * etmp
             - 8.0 * FLAGA * FLAGA * tmp * etmp * etmp * etmp
             + 6.0 * FLAGA * etmp * etmp) * Exp(etmp) * Exp(-2.5);
}

static void    /* boundary function: Dirichlet */
func_g(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT tmp;
    tmp = 1.0 / (FLAGA * ((x - 0.5 - 0.4 * Sin(8.0 * M_PI * crtime)) 
                            * (x - 0.5 - 0.4 * Sin(8.0 * M_PI * crtime))
                        + (y - 0.5 - 0.4 * Cos(8.0 * M_PI * crtime)) 
                            * (y - 0.5 - 0.4 * Cos(8.0 * M_PI * crtime))
                        + (z - 1.0) * (z - 1.0)) + FLAGB + 0.2);
    *value = exp(tmp) * Exp(-2.5);
}

static void /* initial solution u0: t = 0*/
func_u0(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT tmp;
    tmp = 1.0 / (FLAGA * ((x - 0.5) * (x - 0.5)
                        + (y - 0.9) * (y - 0.9)
                        + (z - 1.0) * (z - 1.0)) + FLAGB + 0.2);
    *value = exp(tmp) * Exp(-2.5);
}

static double
cal_time(BOOLEAN flag)
{
    static double t0 = 0.0;
    double et, tt[3];

    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    if (flag){
        phgPrintf("  %0.4lfs",et);
    }

    return et;
}

/* f_h: \bar{f}^n */
/* mean value of f in [t_(n - 1), t_n] */
/* use Simpson formula */
static void
update_fh(DOF *f_h, DOF *ftmp)
{
    phgDofSetDataByValue(f_h, 0.0);
    phgDofSetDataByFunction(ftmp, func_f);
    phgDofAXPY(1.0 / 6.0, ftmp, &f_h);
    crtime -= stime * 0.5;
    phgDofSetDataByFunction(ftmp, func_f);
    phgDofAXPY(2.0 / 3.0, ftmp, &f_h);
    crtime -= stime * 0.5;
    phgDofSetDataByFunction(ftmp, func_f);
    phgDofAXPY(1.0 / 6.0, ftmp, &f_h);
    crtime += stime;
}


/* build linear systme */
static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f, DOF *u_p)
{
    int N = u_h->type->nbas;    /* number of basis functions in an element */
    int dim = u_h->dim;         /* DOF dimension */
    GRID *g = u_h->g;
    ELEMENT *e;
    int i, j, k, n;
    FLOAT a[N][N], A[N][dim], B[N], buffer[N], *row;
    INT I[N][dim], J[N][dim];
    FLOAT tmp;
   
    ForAllElements(g, e) {
        for (i = 0; i < N; i++) {
            for (k = 0; k < dim; k++)
                I[i][k] = phgSolverMapE2L(solver, 0, e, i * dim + k);
            for (j = 0; j < i; j++) {
                a[j][i] = a[i][j] =
                    phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, QUAD_DEFAULT) * stime
                    + phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
            }
            a[i][i] = phgQuadGradBasDotGradBas(e, u_h, i, u_h, i, QUAD_DEFAULT) * stime
                    + phgQuadBasDotBas(e, u_h, i, u_h, i, QUAD_DEFAULT);
        }

        /* loop on basis functions */
        for (i = 0; i < N; i++) {
            if (phgDofDirichletBC(u_h, e, i, func_g, buffer, &B[i],
                                DOF_PROJ_NONE)) {
                row = buffer;
            }
            else if (phgDofGetElementBoundaryType(u_h, e, i * dim) & NEUMANN) {
                /* TODO */
                row = NULL;
                phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__);
            }
            else {
                /* right hand side: \int f * phi_i */
                phgQuadDofTimesBas(e, f, u_h, i, QUAD_DEFAULT, &B[i]);
                B[i] = stime * B[i];
                phgQuadDofTimesBas(e, u_p, u_h, i, QUAD_DEFAULT, &tmp);
                B[i] += tmp;
                row = a[i];
            }
            for (j = 0, n = 0; j < N; j++) {
                for (k = 0; k < dim; k++) {
                    A[n][k] = row[j];
                    J[n++][k] = I[j][k];
                }
            }
            phgSolverAddMatrixEntries(solver, 1, I[i], n, J[0], A[0]); 
        }
        phgSolverAddRHSEntries(solver, N * dim, I[0], &B[0]);
    }
}

static void  /* solve Au = b */
SolveLinearSystem(OEM_SOLVER *sol, DOF *u_h, DOF *f_h, DOF *u_p)
{
    GRID *g = u_h->g;
    SOLVER *solver;
    int nits;
    phgPrintf("Number of DOF:           %d\n", DofDataCountGlobal(u_h));
    phgPrintf("Number of elements:      %d\n", g->nleaf_global);
    phgPrintf("Load imbalance:          %lg\n", (double)g->lif);
                       
    phgPrintf("create solver:         ");
    cal_time(FALSE);
    solver = phgSolverCreate(sol, u_h, NULL);
    cal_time(TRUE);
    phgPrintf("\n");
                
    phgPrintf("build linear system:   ");
    cal_time(FALSE);
    build_linear_system(solver, u_h, f_h, u_p);
    cal_time(TRUE);
    phgPrintf("\n");
              
    phgPrintf("solve solver:          ");
    cal_time(FALSE);
    nits = phgSolverSolve(solver, TRUE, u_h, NULL);
    cal_time(TRUE);
    phgPrintf("     [%d]    \n", nits);
    phgSolverDestroy(&solver);

}

/* time residual, global */
static FLOAT
estimate_time_error(DOF *u_h, DOF *u_p)
{
    DOF *tmp, *grad;
    FLOAT ss;

    tmp = phgDofCopy(u_h, NULL, DOF_DEFAULT, "tmp2");
    phgDofAXPY( -1.0, u_p, &tmp);
    grad = phgDofGradient(tmp, NULL, NULL, "grad");
    ss = phgDofNormL2(grad);

    phgDofFree(&tmp);
    phgDofFree(&grad);

    return ss * ss * (FLOAT)(1.L/3.L);
}


/* global source error, use Simpson's rule 
 * divide [ptime, crtime] into [ptime, ptime + stime / 3.0]
 * [ptime + stime / 3.0, ptime + 2.0 * stime / 3.0]
 * [ptime + 2.0 * stime / 3.0, crtime]
 * use Simpson's rule separately
 */
static FLOAT
estimate_source_error(DOF *f_h)
{
    GRID *g = f_h->g;
    DOF *tmp;
    FLOAT result = 0.0, time, h;

    time = crtime;
    h = stime / 6.0;
   
    /* P:1/3 */
    tmp = phgDofNew(g, DOF_DEFAULT, 1, "tmp", DofNoAction);
    phgDofSetDataByFunction(tmp, func_f);
    phgDofAXPY(-1.0, f_h, &tmp);
    result += phgDofNormL2(tmp);

    crtime -= h;
    phgDofSetDataByFunction(tmp, func_f);
    phgDofAXPY(-1.0, f_h, &tmp);
    result += 4.0 * phgDofNormL2(tmp);
    
    crtime -= h;
    phgDofSetDataByFunction(tmp, func_f);
    phgDofAXPY(-1.0, f_h, &tmp);
    result +=  2.0 * phgDofNormL2(tmp);
    
    /* M: 1/ 3 */
    crtime -= h;
    phgDofSetDataByFunction(tmp, func_f);
    phgDofAXPY(-1.0, f_h, &tmp);
    result +=  4.0 * phgDofNormL2(tmp);
    
    crtime -= h;
    phgDofSetDataByFunction(tmp, func_f);
    phgDofAXPY(-1.0, f_h, &tmp);
    result += 2.0 *  phgDofNormL2(tmp);
   
    /* L:1/3 */
    crtime -= h;
    phgDofSetDataByFunction(tmp, func_f);
    phgDofAXPY(-1.0, f_h, &tmp);
    result += 4.0 * phgDofNormL2(tmp);
    
    crtime -= h;
    phgDofSetDataByFunction(tmp, func_f);
    phgDofAXPY(-1.0, f_h, &tmp);
    result += phgDofNormL2(tmp);
    phgDofFree(&tmp);
    
    crtime = time;
    return result / 18.0;
}

/* space error indicator*/
static FLOAT
estimate_space_error(DOF *u_h, DOF *u_p, DOF *f_h, DOF *error)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *tmp, *jump;
    FLOAT eta, h;
    int i;
  
    tmp = phgDofGradient(u_h, NULL, NULL, "tmp");
    jump = phgQuadFaceJump(tmp, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    phgDofFree(&tmp);

    eta = 1.0 / stime;
    tmp = phgDofCopy(f_h, NULL, f_h->type, NULL);
    phgDofAXPY(eta, u_p, &tmp);
    eta *= -1.0;
    phgDofAXPY(eta, u_h, &tmp);

    ForAllElements(g, e){
        eta = 0.0;
        for (i = 0; i < NFace; i++) {
            if (e->bound_type[i] & (DIRICHLET | NEUMANN))
                continue;    /* boundary face */
            h = phgGeomGetFaceDiameter(g, e, i);
            eta +=  (*DofFaceData(jump, e->faces[i])) * h;
        }
        h = phgGeomGetDiameter(g, e);
        eta += eta * 0.5  + h * h * phgQuadDofDotDof(e, tmp, tmp, QUAD_DEFAULT);

        *DofElementData(error, e->index) = eta;
    }
    phgDofFree(&tmp);
    phgDofFree(&jump);
    return phgDofNormInftyVec(error);
}

/* data oscillation */
static FLOAT
estimate_osc(DOF *f_h, DOF *u_h, DOF *u_p, DOF *osc)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *res, *cons;
    double pk, tmp;
    double h, v;
    
    tmp = 1 / stime;
    res = phgDofCopy(f_h, NULL, DOF_DEFAULT, "res");
    phgDofAXPY(tmp, u_p, &res);
    tmp *= -1.0;
    phgDofAXPY(tmp, u_h, &res);

    cons = phgDofNew(g, DOF_DEFAULT, 1, "cons", DofNoAction);
    phgDofSetDataByValue(cons, 1.0);
    
    ForAllElements(g, e){
        h = phgGeomGetDiameter(g, e);
        v = phgGeomGetVolume(g, e);
        pk = phgQuadDofDotDof(e, res, cons, -1) / v;
        tmp = (phgQuadDofDotDof(e, res, res, -1) - pk * pk * v) * h * h;
        *DofElementData(osc, e->index) = tmp;
    }
    phgDofFree(&res);
    phgDofFree(&cons);
    return phgDofNormInftyVec(osc);
}

/* calculate inital error */
static FLOAT
Init_error(DOF *u_h, DOF *u_t, DOF *error)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    FLOAT tmp;
    
    ForAllElements(g, e){
        tmp = 0.0;
        tmp += phgQuadDofDotDof(e, u_h, u_h, 4);
        tmp -= 2.0 * phgQuadDofDotDof(e, u_h, u_t, 4);
        tmp += phgQuadDofDotDof(e, u_t, u_t, 4);
        *DofElementData(error, e->index) = tmp;
    }

    return phgDofNormInftyVec(error);
}

/* calculate L2 norm of (u_t - u_h) */
/* u_t's type is DOF_ANALYTIC */
static FLOAT
L2_error(DOF *u, DOF *u_h)
{
    ELEMENT *e;
    FLOAT tmp, result;

    tmp = 0.0;
    ForAllElements(u->g, e){
        tmp += phgQuadDofDotDof(e, u_h, u_h, 5);
        tmp -= 2.0 * phgQuadDofDotDof(e, u_h, u, 5);
        tmp += phgQuadDofDotDof(e, u, u, 5);
    }
#if USE_MPI
    MPI_Allreduce(&tmp, &result, 1, PHG_MPI_FLOAT, PHG_SUM, u->g->comm);
#else
    result = tmp;
#endif
    return Sqrt(result);
}

/* cal L2 norm for DOF_ANALYTIC */
static FLOAT
L2_norm(DOF *u)  
{
    ELEMENT *e;
    FLOAT tmp, result;

    tmp = 0.0;
    ForAllElements(u->g, e){
        tmp += phgQuadDofDotDof(e, u, u, 5);
    }
#if USE_MPI
    MPI_Allreduce(&tmp, &result, 1, PHG_MPI_FLOAT, PHG_SUM, u->g->comm);
#else
    result = tmp;
#endif
    return Sqrt(result);
}

static void
result_print(DOF *u_t, DOF *u_h, DOF *error, 
          FLOAT PE_error, FLOAT etaspace, FLOAT etaf, FLOAT etatime)
{
    size_t peakmem;
    FLOAT tmp, tmp1;
    
    phgPrintf("time step :     %lg\n", (double)stime);
    phgPrintf("source error:  %*lf         tol:%lg\n",
                               9, (double)Sqrt(etaf), (double)Sqrt(FTOL));
    phgPrintf("time residual: %*lf         tol:%lg\n",
                               9, (double)Sqrt(etatime), (double)Sqrt(TimeTOL));
    
    phgPrintf("max space error:         %lg\n", (double)Sqrt(PE_error));
    phgPrintf("space error indicator:  %*lf  tol:%lg\n",
                       9,(double)Sqrt(etaspace), (double)Sqrt(SpaceTOL));
        
    tmp = L2_norm(u_t);
    tmp1 = L2_error(u_t, u_h);
    phgPrintf("true solution norm:      %0.10lf\n", (double)tmp);
    phgPrintf("numerical solution norm: %0.10lf\n", (double)L2_norm(u_h));
    phgPrintf("True L2 Error:           %0.10lf\n", (double)tmp1);
    phgPrintf("relative Error:          %0.10lf\n", (double)(tmp1/tmp));
    phgMemoryUsage(u_h->g, &peakmem);
    phgPrintf("memory usage: %lf Mb\n", peakmem / (1024. * 1024.));
}

int
main(int argc, char *argv[])
{

    GRID *g;
    ELEMENT *e;
    DOF *u_h;  /* numerical solution at time t_n */
    DOF *u_t;  /* true solution at time t_n */
    DOF *u_p;  /* numerical solution at time t_{n-1} */
    DOF *error, *osc; /* error indicator and data oscillation */
    DOF *ftmp, *f_h;
    FLOAT etaf, etatime, etaspace;  /* sum of error indicator */
    FLOAT PE_error;                 /* max space indicator */
    FLOAT tmp, tmp1;
    FLOAT delta1 = 0.95;    /* time-step size control parameter */
    FLOAT delta2 = 1.1;     /* time-step size control parameter */
    FLOAT thetatime = 0.5;  /* time-step size control parameter */
    double tt[3], tt1[3], tt2[3];
    char *fn = "../test/cube.dat";
    int flag = 0;
    long long Nelem = 0; /* total elements */
    long long Ndof = 0;  /* total DOF */
    char sn[40];         /* output file name */
    static BOOLEAN export_mesh = FALSE;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterNoArg("-export_mesh", "Export mesh", &export_mesh);
    phgOptionsRegisterFloat("-tol", "Convergence criterion", &htolerance);
    phgOptionsRegisterFloat("-T", "computaional domain: [0, T]", &T);
    phgOptionsRegisterFloat("-step", "initial time step: 0. < step < T", &stime);

    phgInit(&argc, &argv); 
    g = phgNewGrid(-1);
    if (!phgImport(g, fn, TRUE))
        phgError(1, "can't read file \"%s\".\n", fn);

    u_t = phgDofNew(g, DOF_ANALYTIC, 1, "u_t", func_g);
    u_p = phgDofNew(g, DOF_DEFAULT, 1, "u_p", func_u0);
    u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
    f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h",  DofNoAction);
    ftmp = phgDofNew(g, DOF_DEFAULT, 1, "ftmp", DofNoAction);
    error = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);
    osc = phgDofNew(g, DOF_P0, 1, "oscillaiton", DofNoAction);

    /* verify parameters */
    if (htolerance <= 0.) {
        phgWarning("*** tol <= 0.!\n");
        phgWarning("*** change tol to 0.1 now\n");
        htolerance = 0.1;
    }
    if (T <= 0.) {
        phgWarning("*** T <= 0.!\n");
        phgWarning("*** change T to 1. now\n");
        T = 1.;
    }
    if (stime <= 0. || stime >= T) {
        phgWarning("*** Wrong time step!\n");
        phgWarning("*** change time step to %f now\n", T / 10.);
        stime = T / 10;
    }

    /* init parameters */
    tol_space = tol_time = htolerance;
    SpaceTOL = (tol_space) / T;         /* TOL / T, */
    TimeTOL = (tol_space) / (2.0 * T) ;  /* TOL / (2 * T) */
    FTOL = Sqrt(tol_time) / (2.0 * T);


    phgGetTime(tt);
    phgPrintf("Number of processes:%d\n", phgNProcs);

    /* init U0 such that L2 norm (u0 - U0) < Sqrt(tol_space) */
    etaspace = 2 * tol_space;
    while(etaspace > tol_space * 0.25){
        PE_error = Init_error(u_p, u_t, error);
        etaspace = phgDofNormL1Vec(error);
        phgPrintf("Init grid. Error:%lg    tol:%lg\n",
                                (double)Sqrt(etaspace), (double)Sqrt(tol_space));
        phgPrintf("Number of DOF:           %d\n", DofDataCountGlobal(u_h));
        phgPrintf("Number of elements:      %d\n\n", g->nleaf_global);
                
        tmp = 0.5 * PE_error;
        ForAllElements(g, e){
            e->mark = 0;
            if (*DofElementData(error, e->index) > tmp)
                e->mark = 1;
        }
        phgRefineMarkedElements(g);
        if (phgBalanceGrid(g, 1.2, -1, NULL, 0.))
            phgPrintf("Repartition mesh, lif:   %lg\n\n", (double)g->lif);
    }
    phgPrintf("Init grid. Error:%lg    tol:%lg\n",
                                (double)Sqrt(etaspace), (double)Sqrt(tol_space));

    while(crtime < T - 1e-8){
        phgGetTime(tt1);
        ptime = crtime;
        crtime += stime;
        phgPrintf("\n/********* start new time layer *************/\n");
        phgPrintf("current time layer: [%lf, %lf]\n",
                                (double)ptime, (double)crtime);

        if (crtime > T){
            crtime = T;
            stime = T - ptime;
            phgPrintf("current time layer: [%lf, %lf]\n",
                                (double)ptime, (double)crtime);
        }

        flag++;
        if (flag > 1){   /* update u_p */
            phgDofFree(&u_p);
            u_p = phgDofCopy(u_h, NULL, DOF_DEFAULT, "u_p2");
            u_p->userfunc = DofInterpolation;
        }

        /* update f_h */
        update_fh(f_h, ftmp);
        
        SolveLinearSystem(SOLVER_DEFAULT, u_h, f_h, u_p);
        
        /* calculate error indicators */
        etaf = estimate_source_error(f_h);
        etatime = estimate_time_error(u_h, u_p);
        PE_error = estimate_space_error(u_h, u_p, f_h, error);
        etaspace = phgDofNormL1Vec(error);
        
        /* output result */
        result_print(u_t, u_h, error, PE_error, etaspace, etaf, etatime);
        
        /* time step control */
        while(etaf > FTOL || etatime > TimeTOL){
            phgPrintf("\n&&&&&&&& adjust time step &&&&&&&&&\n");
            stime *= delta1;
            crtime = stime + ptime;

            phgPrintf("time step: %lg\n", (double)stime);
            phgPrintf("current time layer: [%lf, %lf]\n",
                                (double)ptime, (double)crtime);

            /* update f_h */
            update_fh(f_h, ftmp);
            
            SolveLinearSystem(SOLVER_DEFAULT, u_h, f_h, u_p);
            
            /* calculate errors */
            etaf = estimate_source_error(f_h);
            etatime = estimate_time_error(u_h, u_p);
            PE_error = estimate_space_error(u_h, u_p, f_h, error);
            etaspace = phgDofNormL1Vec(error);
        
            /* output result */
            result_print(u_t, u_h, error, PE_error, etaspace, etaf, etatime);
        }
        
        phgPrintf("\ntime step: %lg\n\n", (double)stime);

        /* space adativity */
        while(etaspace > SpaceTOL){
            /* refine and coarsen elements */
            phgPrintf("\n=============== adjust grid =================\n");
            estimate_osc(f_h, u_h, u_p, osc);
            phgMarkElements(MARK_DEFAULT, error, Pow(0.3,2), osc, Pow(0.3,2), 1,
			    MARK_DEFAULT, error, Pow(0.1,2), 1,
			    Pow(tol_space,2) / g->nleaf_global);
            phgPrintf("refine elements: ");
            cal_time(FALSE);
            phgRefineMarkedElements(g);
            cal_time(TRUE);
            phgPrintf("\ncoarsen elements:");
            cal_time(FALSE);
            phgCoarsenMarkedElements(g);
            cal_time(TRUE);

            /*phgVerbosity = 0;*/
            cal_time(FALSE);
            if (phgBalanceGrid(g, 1.2, -1, NULL, 0.)) {
		phgPrintf("\nredistribute grid:");
		cal_time(TRUE);
            }
            phgPrintf("\n\n");
            /*phgVerbosity = 0;*/
            /* updat f_h */
            update_fh(f_h, ftmp);

            SolveLinearSystem(SOLVER_DEFAULT, u_h, f_h, u_p);

            /* calculate errors */
            etaf = estimate_source_error(f_h);
            etatime = estimate_time_error(u_h, u_p);
            PE_error = estimate_space_error(u_h, u_p, f_h, error);
            etaspace = phgDofNormL1Vec(error);
        
            /* output result */
            result_print(u_t, u_h, error, PE_error, etaspace, etaf, etatime);

            /* time step control */
            while(etaf > FTOL || etatime > TimeTOL){
                phgPrintf("\n&&&&&&&&& adjust time step &&&&&&&&&\n");
                stime *= delta1;
                crtime = stime + ptime;
                phgPrintf("time step: %lg\n", (double)stime);
                phgPrintf("current time: [%le, %le]\n",
                                (double)ptime, (double)crtime);

                /* update f_h */
                update_fh(f_h, ftmp);

                SolveLinearSystem(SOLVER_DEFAULT, u_h, f_h, u_p);
           
                /* cal error indicators */
                etaf = estimate_source_error(f_h);
                etatime = estimate_time_error(u_h, u_p);
                PE_error = estimate_space_error(u_h, u_p, f_h, error);
                etaspace = phgDofNormL1Vec(error);
        
                /* output result */
                result_print(u_t, u_h, error, PE_error, etaspace, etaf, etatime);
            }
            phgPrintf("max error:%lg   etaspace:%lg  tol:%lg\n",
                        (double)Sqrt(PE_error), (double)Sqrt(etaspace), (double)Sqrt(SpaceTOL));
        }
        

        phgDofSetDataByFunction(ftmp, func_g);
        phgDofAXPY(-1.0, u_h, &ftmp);
        if (export_mesh) {
            *sn = '\0';
            if (flag % 40 == 0){
                sprintf(sn, "heat-4-%5.5d.vtk", flag / 40);
                phgPrintf("output mesh to %s:\n", sn);
                phgExportVTK(g, sn, u_h, NULL);
            }
        }

        tmp = thetatime * tol_time;
        tmp1 = Sqrt(tmp) / (2.0 * T);
        tmp = thetatime * TimeTOL;

        /* increase time step if time error indicator is too small */
        if (etatime <= tmp && etaf <= tmp1){
            stime *= delta2;
            phgPrintf("adjust time-step to: %lf\n", (double)stime);
        }

        Nelem += g->nleaf_global;
        Ndof += DofDataCountGlobal(u_h);
        phgGetTime(tt2);
        phgPrintf("\ntime usage of this step: %lfs\n",
                                (double)(tt2[2] - tt1[2]));
    }

    
    /* cal total errors */
    phgPrintf("\nTotal time steps:%d\n", flag);
    phgPrintf("\nTotal elements:%ld\n", Nelem);
    Nelem /= flag;
    phgPrintf("Average elements:%ld\n", Nelem);
    phgPrintf("\nTotal DOF:%ld\n", Ndof);
    Ndof /= flag;
    phgPrintf("Average DOF:%ld\n", Ndof);
    phgGetTime(tt2);
    phgPrintf("\nTotal time:%lfs\n", (double)(tt2[2] - tt[2]));

    phgDofFree(&ftmp);
    phgDofFree(&osc);
    phgDofFree(&u_h);
    phgDofFree(&u_p);
    phgDofFree(&f_h);
    phgDofFree(&u_t);
    phgDofFree(&error);
    phgFreeGrid(&g);
    phgFinalize();
    return 0;
}
