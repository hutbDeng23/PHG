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

/*******************************************************
* $Id: non-smooth.c,v 1.64 2021/03/11 08:09:53 zlb Exp $
*
* non-smooth coefficients elliptic PDE
* 
* -div(a(x) grad u) = f  x \in \Omega
* u = g  x \in \partial \Omega
*
* a(x): piesewise constant  a(x) > \alpha > 0
* a(x) = R, (0, 1)^3
* a(x) = 1, (0, 1)^3 \ (0.5, 1)^3
*
* ture solution is hard to find, thus we use a high-order 
* one to approximate.
* 
*******************************************************/

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <math.h>

/* discontinuous coefficient */
static FLOAT R = 1000.;

#if 1
static void   /* rhs */
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT r2;
    if ((x >= 0.5) && (y >= 0.5) && (z >= 0.5)) {
        r2 =  ((x - 0.8) * (x - 0.8)  + (y - 0.8) * (y - 0.8) 
                + (z - 0.8) * (z - 0.8)) * 40 + 0.2;
        *value = (1./r2);
    }
    else {
        *value = 0.;
    }

    *value *= 0.01;
}

static void  /* boundary function */
func_g(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT r2;

    r2 = ((x - 0.8) * (x - 0.8) + (y - 0.8) * (y - 0.8)
        + (z - 0.8) * (z - 0.8)) * 40 + 1.;
    *value = Exp(1 / r2);
    *value *= 0.01;
}

static void /* coefficent a(x)*/
func_ax(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if ((x >= 0.5) && (y >= 0.5) && (z >= 0.5)) {
        *value = R;
    }
    else {
        *value = 1.0;
    }
}

static void /* coefficent a(x), 3 dim */
func_ax3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if ((x >= 0.5) && (y >= 0.5) && (z >= 0.5)) {
        *(value++) = R;
        *(value++) = R;
        *(value++) = R;
    }
    else {
        *(value++) = 1.0;
        *(value++) = 1.0;
        *(value++) = 1.0;
    }
}
#else
static void   /* rhs */
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT r2;
    if (x <= 0.5) {
        r2 =  ((x - 0.4) * (x - 0.4)  + (y - 0.5) * (y - 0.5) 
                + (z - 0.5) * (z - 0.5)) * 40 + 0.2;
        *value = (1./r2);
    }
    else {
        *value = 0.;
    }

    *value *= 0.01;
}

static void  /* boundary function */
func_g(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT r2;

    r2 = ((x - 0.) * (x - 0.) + (y - 0.5) * (y - 0.5)
        + (z - 0.5) * (z - 0.5)) * 40 + 1.;
    *value = Exp(1 / r2);
    *value *= 0.01;
}

static void /* coefficent a(x)*/
func_ax(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if (x <= 0.5) {
        *value = R;
    }
    else {
        *value = 1.0;
    }
}

static void /* coefficent a(x), 3 dim */
func_ax3(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    if (x <= 0.5) {
        *(value++) = R;
        *(value++) = R;
        *(value++) = R;
    }
    else {
        *(value++) = 1.0;
        *(value++) = 1.0;
        *(value++) = 1.0;
    }
}
#endif

static FLOAT  /* calculate time */
cal_time(BOOLEAN flag)
{
    static FLOAT t0 = 0.0;
    double et, tt[3];

    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    if (flag){
        phgPrintf("%0.6lfs", (double)et);
    }

    return et;
}


/* build linear systme */
static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h, DOF *a)
{
    int N = u_h->type->nbas;    /* number of basis functions in an element */
    GRID *g = u_h->g;
    ELEMENT *e;
    int i, j;
    FLOAT A[N][N], B[N], buffer[N], *row, tmp;
    INT I[N];

    assert(u_h->dim == 1);

    ForAllElements(g, e) {
    /* compute \int a(x) \grad\phi_j \cdot \grad\phi_i making use of symmetry */
        tmp = *DofElementData(a, e->index);
        for (i = 0; i < N; i++) {
            I[i] = phgSolverMapE2L(solver, 0, e, i);
            for (j = 0; j <= i; j++)
                A[j][i] = A[i][j] =
                tmp * phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
        }

        /* loop on basis functions */
        for (i = 0; i < N; i++) {
            if (phgDofDirichletBC(u_h, e, i, func_g,
                  buffer, B + i, DOF_PROJ_NONE)) {
                row = buffer;
            }
            else {
                /* right hand side: \int f * phi_i */
                phgQuadDofTimesBas(e, f_h, u_h, i, QUAD_DEFAULT, B + i);
                row = A[i];
            }
            phgSolverAddMatrixEntries(solver, 1, I + i, N, I, row); 
        }
        phgSolverAddRHSEntries(solver, N, I, B);
    }
}

static void  /* solve Au = b */
SolveEqua(OEM_SOLVER *sol, DOF *u_h, DOF *f_h, DOF *ax)
{
    GRID *g = u_h->g;
    SOLVER *solver;
    int nits;
    phgPrintf("DOF: %"dFMT" ", DofDataCountGlobal(u_h));
    phgPrintf("Elements: %"dFMT"  ", g->nleaf_global);
    phgPrintf("Load imbalance: %lg\n",  (double)g->lif);
                       
    solver = phgSolverCreate(sol, u_h, NULL);
    phgPrintf("build linear system:   ");
    cal_time(FALSE);
    build_linear_system(solver, u_h, f_h, ax);
    cal_time(TRUE);
              
    phgPrintf(" solve equation:  ");
    cal_time(FALSE);
    nits = phgSolverSolve(solver, TRUE, u_h, NULL);
    cal_time(TRUE);
    phgPrintf(" [%d]\n", nits);
    phgSolverDestroy(&solver);

}

/* calulate error indicator */
static FLOAT
estimate_space(DOF *u_h, DOF *f_h, DOF *a, DOF *ax3, DOF *error)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    NEIGHBOUR_DATA *aa;
    DOF *jump, *grad, *agrad, *ddiv, *td;
    FLOAT sum, tmp;
    FLOAT h, ak, ae;
    int i;

    grad = phgDofGradient(u_h, NULL, NULL, "grad");
    agrad = phgDofMM(MAT_OP_N, MAT_OP_N, 3, 1, 3, 1.0, ax3, 1, grad, 0.0, NULL);
    ddiv = phgDofGetSameOrderDG(agrad, DofDim(f_h), "div");

    phgDofCopy(f_h, &ddiv, NULL, NULL);
    td = phgDofDivergence(agrad, NULL, NULL, "div");
    phgDofAXPY(1.0, td, &ddiv);
    phgDofFree(&td);

    jump = phgQuadFaceJump(agrad, DOF_PROJ_DOT, "jump", -1);

    phgDofSetDataByValue(error, 0.0);
    aa = phgDofInitNeighbourData(a, NULL);  
    ForAllElements(g, e){ /* neighbours are local */
        sum = 0.0;
        tmp = *DofElementData(a, e->index);
        h = phgGeomGetDiameter(g, e);
        sum = h * h / tmp * phgQuadDofDotDof(e, ddiv, ddiv, -1);
        for (i = 0; i < NFace; i++) {
            if (e->bound_type[i] & (DIRICHLET | NEUMANN))
                continue;    /* boundary face */
            ak = *phgDofNeighbourData(aa, e, i, 0, NULL);
            ae = (tmp > ak) ? tmp : ak;
            h = phgGeomGetFaceDiameter(g, e, i);
            sum +=  (*DofFaceData(jump, e->faces[i])) * h * 0.5 / ae;
        }
        *DofElementData(error, e->index) = sum;
    }
    
    phgDofFree(&ddiv);
    phgDofFree(&jump);
    phgDofFree(&grad);
    phgDofFree(&agrad);
    phgDofReleaseNeighbourData(&aa);

    return phgDofNormInftyVec(error);

}

/* data oscillation */
static FLOAT
estimate_osc(DOF *f_h, DOF *f_t, DOF *a, DOF *osc)
{
    GRID *g = f_h->g;
    ELEMENT *e;
    FLOAT tmp, h, aa;
    INT order;
    
    assert(!SpecialDofType(f_h->type));
    ForAllElements(f_h->g, e){
        tmp = 0.0;
        h = phgGeomGetDiameter(g, e);
        aa = *DofElementData(a, e->index);
        order = DofElementOrder(f_h, e->index) * 2;
        tmp += phgQuadDofDotDof(e, f_h, f_h, order);
        tmp -= 2.0 * phgQuadDofDotDof(e, f_h, f_t, order);
        tmp += phgQuadDofDotDof(e, f_t, f_t, order);
        *DofElementData(osc, e->index) = Fabs(tmp) * h * h / aa;
    }
    return phgDofNormInftyVec(osc);
}
/* calcute L2 norm 
 * u1 != NULL
 * if (u2 != NULL), calculate L2(u1 - u2)
*/
FLOAT
Norm_L2(DOF *u1, DOF *u2)
{
    ELEMENT *e;
    FLOAT result;
    FLOAT tmp;

    assert(u1 != NULL);

    if (u2 == NULL) {
        result = phgDofNormL2(u1);
    }
    else {
        tmp = 0.;
        ForAllElements(u1->g, e) {
            tmp += phgQuadDofDotDof(e, u1, u1, -1);
            tmp -= 2. * phgQuadDofDotDof(e, u1, u2, -1);
            tmp += phgQuadDofDotDof(e, u2, u2, -1);
        }
#if USE_MPI
        MPI_Allreduce(&tmp, &result, 1, PHG_MPI_FLOAT, PHG_SUM, u1->g->comm);
#else
        result = tmp;
#endif
        result = Sqrt(Fabs(result));
    }

    return result;
}

/* cal energy norm */
/* (\int a(x) \grad u dot \grad u dx */
static FLOAT 
Norm_Energy(DOF *u, DOF *u2, DOF *a)
{
    GRID *g;
    ELEMENT *e;
    DOF *grad = NULL, *grad2 = NULL;
    FLOAT sum = 0.0;
    FLOAT tmp = 0.0;

    assert(u != NULL);
    g = u->g;

    if (u2 == NULL) {
        grad = phgDofGradient(u, NULL, NULL, "u grad");

        ForAllElements(g, e){
            tmp = *DofElementData(a, e->index);
            sum += tmp * phgQuadDofDotDof(e, grad, grad, -1);
        }
    }
    else {
        grad = phgDofGradient(u, NULL, NULL, "u grad");
        grad2 = phgDofGradient(u2, NULL, NULL, "u2 grad");

        ForAllElements(g, e){
            tmp = *DofElementData(a, e->index);
            sum += tmp * phgQuadDofDotDof(e, grad, grad, -1);
            sum += tmp * phgQuadDofDotDof(e, grad2, grad2, -1);
            sum -= 2. * tmp * phgQuadDofDotDof(e, grad, grad2, -1);
        }
    }
#if USE_MPI
    MPI_Allreduce(&sum, &tmp, 1, PHG_MPI_FLOAT, PHG_SUM, g->comm);
    sum = tmp;
#endif

    phgDofFree(&grad);
    phgDofFree(&grad2);

    return Sqrt(Fabs(sum));
}

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *u_h, *u_t;
    DOF *f_h, *f_t;
    DOF *ax, *ax3;
    DOF *error, *osc;
    FLOAT etaspace, etamax, tmp, tp;
    double tt1[3], tt2[3], tt3[3];
    size_t peakmem;
    static char *fn = "../test/non-smooth.dat";
    static FLOAT tol = 3e-3;
    static FLOAT theta = 0.3;
    static BOOLEAN export_mesh = TRUE;
    static BOOLEAN U_T = TRUE;
    int flag = 0;
    INT order = 2;
    int pstop = 0;
    FLOAT  err_save = 1.;
    
    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterFloat("-tol", "Convergence criterion", &tol);
    phgOptionsRegisterFloat("-theta", "GERS parameter", &theta);
    phgOptionsRegisterFloat("-ns", "non-smooth parameter( > 0)", &R);
    phgOptionsRegisterInt("-order", "order of FEM 1<= order <= 3", &order);
    phgOptionsRegisterNoArg("-export_mesh", "Export mesh", &export_mesh);
    phgOptionsRegisterNoArg("-u_t", "solve u_t", &U_T);


    phgInit(&argc, &argv);
    g = phgNewGrid(-1);

    if (!phgImport(g, fn, TRUE))
        phgError(1, "can't read file \"%s\".\n", fn);

    if (R <= 0.)
        phgError(1, "non-smooth parameter ns <= 0.\n");
    if (order < 1) {
        phgWarning("Order is too small. It's set as 2.\n");
        order = 2;
    }
    else if (order > 3) {
        phgWarning("Order is too big. It's set as 2.\n");
        order = 2;
    }

    u_h = phgDofNew(g, DOF_Pn[order], 1, "u_h", DofInterpolation);
    u_t = phgDofNew(g, DOF_Pn[order + 1], 1, "u_t", DofInterpolation);
    f_h = phgDofNew(g, DOF_Pn[order], 1, "f_h", func_f);
    f_t = phgDofNew(g, DOF_ANALYTIC, 1, "f_t", func_f);

    u_h->DB_mask = BDRY_MASK;
    u_t->DB_mask = BDRY_MASK;
    phgDofSetDataByValue(u_t, 3.);

    ax = phgDofNew(g, DOF_P0, 1, "a", func_ax);
    ax3 = phgDofNew(g, DOF_P0, 3, "(a, a, a)", func_ax3);
    error = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);
    osc = phgDofNew(g, DOF_P0, 1, "osc", DofNoAction);

    phgGetTime(tt1);
    while(TRUE){
        phgGetTime(tt2);
        flag++;
        phgPrintf("\n\n****************   Pass %d    ****************\n", flag);

        if (phgBalanceGrid(g, 1.2, -1, NULL, 0.))
            phgPrintf("Repartition mesh, load imbalance: %lg\n",
                (double)g->lif);

        /* solve the linear system */
        phgPrintf("solve u_h(numerical solution):\n");
        SolveEqua(SOLVER_DEFAULT, u_h, f_h, ax);
        if (U_T) {
            phgPrintf("\nsolve u_t(approximately analytic solution):\n");
            SolveEqua(SOLVER_DEFAULT, u_t, f_h, ax);
        }


        /* calculate error */
        etamax = estimate_space(u_h, f_h, ax, ax3, error);
        estimate_osc(f_h, f_t, ax, osc);
        etaspace = phgDofNormL1Vec(error);
        phgPrintf("\nMax error indicator: %.9f  Total error:         %.9f\n",
                        (double)Sqrt(etamax), (double)Sqrt(etaspace));

        /* wether we solve u_t and estimate error or not*/
        if (U_T) {
            /* calculate L2 norm */
            tp = Norm_L2(u_h, NULL);
            tmp = Norm_L2(u_t, NULL);
            phgPrintf("L2 norm(u_h) :       %.9f  ", (double)tp);
            phgPrintf("L2 norm(u_t) :       %.9f\n", (double)tmp);
            tp = Norm_L2(u_t, u_h);
            phgPrintf("True L2 error:       %.9f  ", (double)tp);
            phgPrintf("relative error:      %.9f\n",
                    (double)(tp / tmp));

            /* calculate energy norm */
            tp = Norm_Energy(u_h, NULL, ax);
            tmp = Norm_Energy(u_t, NULL, ax);
            phgPrintf("Energy norm(u_h):    %.9f  ", (double)tp);
            phgPrintf("EnergyNorm(u_t):     %.9f\n", (double)tmp);
            tp = Norm_Energy(u_t, u_h, ax);
            phgPrintf("True energy error:   %.9f  ", (double)tp);
            phgPrintf("relative error:      %.9f\n",
                    (double)(tp / tmp));
        }

        /* stop or not */
        if (flag == 1) err_save = Sqrt(etaspace);

        if (phgRank == 0) {
            if (Sqrt(etaspace) / err_save <= tol){
                pstop = 1;
            }
            else {
                pstop = 0;
            }
        }
#if USE_MPI
        MPI_Bcast(&pstop, 1, MPI_INT, 0, g->comm);
#endif
        if (pstop){
            phgPrintf("\nTolorance: %.9f Error reduction(Relative): %lf\n", 
                    tol, Sqrt(etaspace) / err_save);
            phgGetTime(tt3);
            phgPrintf("\nTime usage:       %lfs\n", (double)(tt3[2] - tt2[2]));
            phgPrintf("Memory usage:     %lfMb\n\n",
                    (double)peakmem / (1024. * 1024.));
            break;
        }
        else {
            phgPrintf("\nTolorance: %*lf Error reduction(Relative): %lf\n", 
                    15, tol, Sqrt(etaspace) / err_save);
        }


        /* mark elements for refining */
        phgMarkRefine(MARK_DEFAULT, error, Pow(theta,2), osc, Pow(theta,2.),
			1, Pow(tol, 2) / g->nleaf_global);

        phgPrintf("Refine elements:  ");
        cal_time(FALSE);

        /* refine elements */
        phgRefineMarkedElements(g);
        cal_time(TRUE);

        phgGetTime(tt3);
        phgMemoryUsage(g, &peakmem);
        phgPrintf("\nTime usage:       %.6fs\n", (double)(tt3[2] - tt2[2]));
        phgPrintf("Memory usage:     %.6fMb\n",
                (double)peakmem / (1024. * 1024.));
    }

    phgPrintf("Total time:         %fs\n\n", (double)(tt3[2] - tt1[2]));

    /* export mesh */
    if (export_mesh){
        DOF *grad;

        grad = phgDofGradient(u_h, NULL, NULL, "grad");
        phgPrintf("Final mesh written to \"%s\".\n\n",
                phgExportVTK(g, "non-smooth.vtk",\
                    u_h, grad, error, ax, NULL));
        phgDofFree(&grad);
    }

    phgDofFree(&u_t);
    phgDofFree(&u_h);
    phgDofFree(&ax);
    phgDofFree(&ax3);
    phgDofFree(&f_h);
    phgDofFree(&error);
    phgDofFree(&f_t);
    phgDofFree(&osc);

    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
