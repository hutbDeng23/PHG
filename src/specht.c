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

/* 3D Specht element.
 *  [1] Bernhard Specht, Modified shape functions for the three-node plate
 * 	bending element passing the patch test,
 * 	Int. J. for Numerical Methods in Engineering, Vol. 26, 705-715 (1988)
 *  [2] 石钟慈，陈绍春，Specht九参数板元的分析，计算数学，1989，第三期，312-318
 *  [3] Ming Wang, Zhong-ci Shi, Jinchao Xu, A new class of Zienkiewicz-type
 *	non-conforming element in any dimensions,
 *	Numer. Math. (2007) 106:335–347
 *  [4] Pingbing Ming, ?????????????????
 *
 * $Id: specht.c,v 1.16 2022/04/09 01:28:34 zlb Exp $
 */

#include "phg.h"

#include <ctype.h>	/* toupper */
#include <string.h>
#include <strings.h>
#include <math.h>

static FLOAT
compute_pij(int i, int j, const FLOAT lambda[], const FLOAT (*J)[Dim + 1])
{
    int n, u, v;
    FLOAT val, D0, D1, D2;	/* \grad (\lambda_i - \lambda_j) */

    val = lambda[i] * lambda[j] * (1.0 + lambda[i] - lambda[j]);

    if (i == j)
	return val;

    /* Note: \grad\lambda_i = J[i][0:2] */
    D0 = J[i][0] - J[j][0];
    D1 = J[i][1] - J[j][1];
    D2 = J[i][2] - J[j][2];
    /* Note: sum of the two opposite edge numbers equals to 5 */
    n = 5 - GetEdgeNo(i,j);
    u = GetEdgeVertex(n, 0);
    v = GetEdgeVertex(n, 1);
    /*assert(i + j + u + v == 6);*/
    val += (30.0 * (lambda[i] - lambda[j])
	    +
	    20.0 * ((D0*J[u][0] + D1*J[u][1] + D2*J[u][2]) *
		    (3.*lambda[u]-1.0) /
		    (J[u][0]*J[u][0] + J[u][1]*J[u][1] + J[u][2]*J[u][2])
		    +
		    (D0*J[v][0] + D1*J[v][1] + D2*J[v][2]) *
		    (3.*lambda[v]-1.0) /
		    (J[v][0]*J[v][0] + J[v][1]*J[v][1] + J[v][2]*J[v][2]))
	    ) * lambda[0]*lambda[1]*lambda[2]*lambda[3];

    return val;
}

static const FLOAT *
Specht_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* computes basis functions */
{
    static FLOAT values[NVert*(Dim+1)];
    /* Note: pij is computed only when lambda or J is changed */
    static FLOAT pij[Dim+1][Dim+1], lambda0[Dim+1] = {0.,0.,0.,0.},
	J0[(Dim+1)*(Dim+1)]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0,};
#if USE_OMP
#pragma omp threadprivate(values, pij, lambda0, J0)
#endif  /* USE_OMP */
    int i, j, k, n;
    GRID *g = dof->g;
    const FLOAT (*J)[Dim + 1] = (void *)phgGeomGetJacobian(g, e);

    if (no1 <= 0)
	no1 = NVert * (Dim + 1);

    if (no0 >= no1)
	return NULL;

    if (memcmp(lambda0, lambda, sizeof(lambda0)) || memcmp(J0, J, sizeof(J0))) {
	memcpy(lambda0, lambda, sizeof(lambda0));
	memcpy(J0, J, sizeof(J0));
	for (i = 0; i <= Dim; i++)
	    for (j = 0; j <= Dim; j++)
		pij[i][j] = compute_pij(i, j, lambda, J);
    }

    for (k = no0; k < no1; k++) {
	i = k / (Dim + 1);	/* vertex no. */
	n = k % (Dim + 1);	/* basis no./vertex */
	if (n == 0) {	/* p_i */
	    values[k] = 0.0;
	    for (j = 0; j <= Dim; j++)
		values[k] += pij[i][j];
	}
	else {		/* p_{ij} */
	    COORD *c = g->verts + e->verts[i], *c1;
	    values[k] = 0.0;
	    --n;
	    for (j = 0; j < Dim; j++) {
		int v = GetFaceVertex(i, j);
		c1 = g->verts + e->verts[v];
		values[k] += pij[i][v] * ((*c1)[n] - (*c)[n]) * 0.5;
	    }
	}
    }

    return values + no0;
}

static void
compute_grad_pij(int i, int j, const FLOAT lambda[], const FLOAT (*J)[Dim + 1],
		 FLOAT vals[])
{
    int u, v, n, ii;
    FLOAT D0, D1, D2;	/* \grad (\lambda_i - \lambda_j) */

    for (ii = 0; ii <= Dim; ii++) {
	if (i == j && ii == i)
	    vals[ii] = 2.*lambda[i];
	else if (ii == i)
	    vals[ii] = lambda[j]*(1.+2.*lambda[i]-lambda[j]);
	else if (ii == j)
	    vals[ii] = lambda[i]*(1.+lambda[i]-2.*lambda[j]);
	else
	    vals[ii] = 0.0;
    }

    if (i == j)
	return;

    /* Note: \grad\lambda_i = J[i][0:2] */
    D0 = J[i][0] - J[j][0];
    D1 = J[i][1] - J[j][1];
    D2 = J[i][2] - J[j][2];
    /* Note: sum of the two opposite edge numbers equals to 5 */
    n = 5 - GetEdgeNo(i,j);
    u = GetEdgeVertex(n, 0);
    v = GetEdgeVertex(n, 1);
    for (ii = 0; ii <= Dim; ii++) {
	FLOAT dq0, d1, d2;
	dq0 = 1.0;
	for (n = 0; n <= Dim; n++)
	    if (ii != n)
		dq0 *= lambda[n];
	if (ii == i)
	    d1 = 30.0;
	else if (ii == j)
	    d1 = -30.0;
	else
	    d1 = 0.0;
	if (ii == u)
	    d2 = 60.0 * (D0*J[u][0] + D1*J[u][1] + D2*J[u][2]) /
		(J[u][0]*J[u][0] + J[u][1]*J[u][1] + J[u][2]*J[u][2]);
	else if (ii == v)
	    d2 = 60.0 * (D0*J[v][0] + D1*J[v][1] + D2*J[v][2]) /
		(J[v][0]*J[v][0] + J[v][1]*J[v][1] + J[v][2]*J[v][2]);
	else
	    d2 = 0.0;
	vals[ii] +=
	       (30.0 * (lambda[i] - lambda[j])
		+
		20.0 * ((D0*J[u][0] + D1*J[u][1] + D2*J[u][2]) *
			(3.*lambda[u]-1.0) /
			(J[u][0]*J[u][0] + J[u][1]*J[u][1] + J[u][2]*J[u][2])
			+
			(D0*J[v][0] + D1*J[v][1] + D2*J[v][2]) *
			(3.*lambda[v]-1.0) /
			(J[v][0]*J[v][0] + J[v][1]*J[v][1] + J[v][2]*J[v][2]))
		+
		lambda[ii] * (d1 + d2) ) * dq0;
    }
}

static const FLOAT *
Specht_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* computes gradient w.r.t lambda of basis functions */
{
    static FLOAT values[NVert*(Dim+1)][Dim+1];
    /* Note: grad_pij is computed only when lambda or J is changed */
    static FLOAT grad_pij[Dim+1][Dim+1][Dim+1], lambda0[Dim+1] = {0.,0.,0.,0.},
	J0[(Dim+1)*(Dim+1)]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0,};
#if USE_OMP
#pragma omp threadprivate(values, grad_pij, lambda0, J0)
#endif  /* USE_OMP */
    int i, j, k, n, ii;
    GRID *g = dof->g;
    const FLOAT (*J)[Dim + 1] = (void *)phgGeomGetJacobian(g, e);

    if (no1 <= 0)
	no1 = NVert * (Dim + 1);

    if (no0 >= no1)
	return NULL;

    if (memcmp(lambda0, lambda, sizeof(lambda0)) || memcmp(J0, J, sizeof(J0))) {
	memcpy(lambda0, lambda, sizeof(lambda0));
	memcpy(J0, J, sizeof(J0));
	for (i = 0; i <= Dim; i++)
	    for (j = 0; j <= Dim; j++)
		compute_grad_pij(i, j, lambda, J, grad_pij[i][j]);
    }

    for (k = no0; k < no1; k++) {
	i = k / (Dim + 1);	/* vertex no. */
	n = k % (Dim + 1);	/* basis no./vertex */
	if (n == 0) {	/* p_i */
	    for (ii = 0; ii <= Dim; ii++) {
		values[k][ii] = 0.0;
		for (j = 0; j <= Dim; j++)
		    values[k][ii] += grad_pij[i][j][ii];
	    }
	}
	else {		/* p_{ij} */
	    COORD *c = g->verts + e->verts[i], *c1;
	    --n;
	    for (ii = 0; ii <= Dim; ii++)
		values[k][ii] = 0.0;
	    for (j = 0; j < Dim; j++) {
		int v = GetFaceVertex(i, j);
		c1 = g->verts + e->verts[v];
		for (ii = 0; ii <= Dim; ii++)
		    values[k][ii] += grad_pij[i][v][ii]*((*c1)[n]-(*c)[n])*0.5;
	    }
	}
    }

    return (FLOAT *)(values + no0);
}

static void
Specht_interpC2F(DOF *dof, ELEMENT *e, FLOAT **old_data, FLOAT **new_data)
{
    phgError(1, "unimplemented function %s.\n", __func__);
}

static void
Specht_interpF2C(DOF *dof, ELEMENT *e, FLOAT **old_data, FLOAT **new_data)
{
    phgError(1, "unimplemented function %s.\n", __func__);
}

static void
Specht_init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	  DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	  const FLOAT *funcvalues, FLOAT *dofvalues, FLOAT **pdofvalues)
{
    GRID *g = dof->g;

    Unused(pdofvalues);

    assert(type == VERTEX);
    assert(funcvalues == NULL || (userfunc == NULL && userfunc_lambda == NULL));

    if (funcvalues == NULL) {
	if (userfunc_lambda != NULL) {
	    FLOAT lambda[Dim + 1] = {0., 0., 0., 0.};
	    lambda[index] = 1.0;
	    userfunc_lambda(dof, e, 0, lambda, dofvalues);
	}
	else {
	    COORD *c = g->verts + e->verts[index];
	    userfunc((*c)[0], (*c)[1], (*c)[2], dofvalues);
        }
    }
    else {
	memcpy(dofvalues, funcvalues, dof->dim * (Dim+1) * sizeof(*funcvalues));
    }
}

DOF_TYPE DOF_SPT_ = {
  DofReserved,
  "SPT",		/* name */
  NULL,			/* points */
  NULL,			/* orders */
  DOF_DG4,		/* gradient type */
  NULL,			/* base type */
  NULL,			/* hp_info */
  Specht_interpC2F,	/* InterpC2F */
  Specht_interpF2C,	/* InterpF2C */
  Specht_init,		/* initialization function */
  Specht_bas,		/* basis functions */
  Specht_grad,		/* gradient of basis functions */
  NULL,
  FE_H1,
  FALSE,		/* is_nodal */
  FALSE,		/* invariant */
  FALSE,		/* free_after_use */
  -1,			/* id */
  NVert*4,		/* nbas */
  5,			/* polynomial order */
  1,			/* derivatives order */
  0,			/* continuity (discontinuous) */
  1,			/* dim */
  4,			/* np_vert */
  0,			/* np_edge */
  0,			/* np_face */
  0			/* np_elem */
};

#ifndef DOF_SPT
extern DOF_TYPE DOF_SPT_;
# define DOF_SPT (&DOF_SPT_)
#endif

DOF_TYPE *DOF_SPTn[] = {DOF_SPT, NULL};

/*-------------------------- The test code ----------------------------*/

/* To compile and run the test code:
 * 	cd phg/test
 * 	ln -s ../src/specht.c .
 *	make USER_CFLAGS=-DSPECHT_TEST=1 -s specht
 *	./specht -mesh_file cubde.dat -refine 3
 */

#ifndef SPECHT_TEST
# define SPECHT_TEST 0
#endif
#if SPECHT_TEST

#include <math.h>

#if 1	/*-----------------------------------------------*/
/* Transcendental function, interpolation order is 3 (L2) or 2 (H1). */
#define F	/*Exp*/Sin
#define DF	/*Exp*/Cos
static void
func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = F(x + 2.*y + 3.*z);
    return;
}

static void
grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    values[0] = DF(x + 2.*y + 3.*z);
    values[1] = 2. * values[0];
    values[2] = 3. * values[0];
    return;
}
#else	/*-----------------------------------------------*/
/* Quadratic function, interpolation error should be zero. */
static void
func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = x*x + 2.*y*y + 3.*z*z;
    return;
}

static void
grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    values[0] = 2. * x;
    values[1] = 4. * y;
    values[2] = 6. * z;
    return;
}
#endif	/*-----------------------------------------------*/

static void
func_specht(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
/* initialization function for the Specht element, it returns in `values[]`
 * the values of the function and its three derivatives at the given point
 * in the order:
 * 	f, f_x, f_y, f_z
 */
{
    func(x, y, z, values);
    grad(x, y, z, values + 1);
}

int
main(int argc, char *argv[])
{
    DOF *u, *u0, *grad_u, *grad_u0;
    GRID *g;
    char *fn = "tetra.dat";
    int refine = 0;
    int flag = 1;
    const char *flags[] = {"diag", "lower", "upper", "full", NULL};

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", &fn);
    phgOptionsRegisterInt("-refine", "Refinement count", &refine);
    phgOptionsRegisterKeyword("-hessian", "form of the Hessian matrix",
				flags, &flag);

    phgOptionsPreset("-dof_type SPT");

    phgInit(&argc, &argv);
    g = phgNewGrid(-1);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgRefineAllElements(g, refine);
    phgBalanceGrid(g, 1.0, 0, NULL, 0.);

    u0 = phgDofNew(g, DOF_P5, 1, NULL, func);
    u = phgDofNew(g, DOF_DEFAULT, 1, NULL, DofNoAction);
#if 1
    /* Initialize u by interpolation */
    phgDofSetDataByFunction(u, u->type == DOF_SPT ? func_specht : func);
#else
    /* Initialize u by L2 projection */
    Unused(func_specht);
    u->DB_mask = 0;
    {
	SOLVER *solver = phgSolverCreate(SOLVER_DEFAULT, u, NULL);
	ELEMENT *e;
	int N = u->type->nbas, i, j;
	FLOAT A[N][N], rhs[N];
	INT I[N];
	ForAllElements(g, e) {
	    for (i = 0; i < N; i++) {
		I[i] = phgSolverMapE2L(solver, 0, e, i);
		for (j = 0; j <= i; j++)
		    A[j][i] = A[i][j] =
			phgQuadBasDotBas(e, u, j, u, i, QUAD_DEFAULT);
	    }
	    for (i = 0; i < N; i++) {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
		/* right hand side = \int f * phi_i */
		phgQuadDofTimesBas(e, u0, u, i, QUAD_DEFAULT, rhs + i);
	    }
	    phgSolverAddRHSEntries(solver, N, I, rhs);
	}
	phgSolverSolve(solver, TRUE, u, NULL);
	phgPrintf("L2 projection: nits = %d, residual = %lg\n",
			solver->nits, (double)solver->residual);
	phgSolverDestroy(&solver);
    }
#endif

#if 1
    /* Check orthogonality of bas funcs, should print a 16x16 identity matrix */
    if (g->rank == 0) {
      ELEMENT *e;
      ForAllElements(g, e) {
	int i, j, k;
	FLOAT lambda[] = {0.,0.,0.,0.}, xyz[Dim];
        FLOAT buffer[u->type->nbas*Dim*Dim];
	const FLOAT *v;
	FLOAT (*J)[Dim + 1] = (void *)phgGeomGetJacobian(g, e);
	if (u->type == DOF_SPT) {
	    printf("\nChecking orthogonality of basis functions:\n");
	    for (i = 0; i < 4; i++) {
		lambda[i] = 1.0;
		/* Print \phi on the vertex */
		v = u->type->BasFuncs(u, e, 0, -1, lambda);
		printf("    p%d  =", i);
		for (j = 0; j < u->type->nbas; j++)
		    printf(" %lg", (double)v[j]);
		printf("\n");
		/* Print \grad\phi = D\phi/D\lambda*J on the vertex */
		v = u->type->BasGrads(u, e, 0, -1, lambda);
		for (k = 0; k < Dim; k++) {
		    printf("    p%d%c =", i, k==0 ? 'x' : (k==1 ? 'y':'z'));
		    for (j = 0; j < u->type->nbas; j++) {
			double d;
			d = v[(Dim + 1) * j + 0] * J[0][k] +
			    v[(Dim + 1) * j + 1] * J[1][k] +
			    v[(Dim + 1) * j + 2] * J[2][k] +
			    v[(Dim + 1) * j + 3] * J[3][k];
			printf(" %lg", fabs(d) < 1e-16 ? 0.0 : d);
		    }	/* for (j ...) */
		    printf("\n");
		}	/* for (k ...) */
		lambda[i] = 0.0;
	    }	/* for (i ...) */
	}

	/* Test phgDofEvalBasisHessian() (checking with finite difference) */
	lambda[0] = 0.1; lambda[1] = 0.2; lambda[2] = 0.3; lambda[3] = 0.4;
	phgGeomLambda2XYZ(g, e, lambda, xyz, xyz + 1, xyz + 2);
	phgDofEvalBasisHessian(u, e, lambda, flags[flag], buffer);
	for (i = 0, v = buffer; i < u->type->nbas; i++) {
	    FLOAT h[Dim][Dim];
	    FLOAT eps = 1e-8;	/* finite difference step */
	    printf("\nChecking the Hessian matrix of basis function %d:\n", i);
	    for (j = 0; j < Dim; j++) {
		const FLOAT *v0;
		xyz[j] += 0.5 * eps;
		phgGeomXYZ2Lambda(g, e, xyz[0], xyz[1], xyz[2], lambda);
		v0 = u->type->BasGrads(u, e, i, i+1, lambda);
		for (k = 0; k < Dim; k++)
		    h[j][k] =
			v0[0]*J[0][k]+v0[1]*J[1][k]+v0[2]*J[2][k]+v0[3]*J[3][k];
		xyz[j] -= eps;
		phgGeomXYZ2Lambda(g, e, xyz[0], xyz[1], xyz[2], lambda);
		v0 = u->type->BasGrads(u, e, i, i+1, lambda);
		for (k = 0; k < Dim; k++)
		    h[j][k] = (h[j][k] - (v0[0]*J[0][k]+v0[1]*J[1][k]
					+v0[2]*J[2][k]+v0[3]*J[3][k])) / eps;
		xyz[j] += 0.5 * eps;
	    }
	    switch (toupper(flags[flag][0])) {
		case 'D':
		    printf("    Dxx=%le, Dyy=%le, Dzz=%le\n",
				v[0] - h[0][0], v[1] - h[1][1], v[2] - h[2][2]);
		    v += 3;
		    break;
		case 'L':
		    printf("    Dxx=%lg\n", v[0] - h[0][0]);
		    printf("    Dyx=%lg, Dyy=%lg\n",
				v[1] - h[1][0], v[2] - h[1][1]);
		    printf("    Dzx=%lg, Dzy=%lg, Dzz=%lg\n",
				v[3] - h[2][0], v[4] - h[2][1],
				v[5] - h[2][2]);
		    v += 6;
		    break;
		case 'U':
		    printf("    Dxx=%lg, Dxy=%lg, Dxz=%lg\n",
				v[0] - h[0][0], v[1] - h[0][1],
				v[2] - h[0][2]);
		    printf("    Dyy=%lg, Dyz=%lg\n",
				v[3] - h[1][1], v[4] - h[1][2]);
		    printf("    Dzz=%lg\n", v[5] - h[2][2]);
		    v += 6;
		    break;
		case 'F':
		case 'A':
		    for (j = 0; j < Dim; j++) {
			for (k = 0; k < Dim; k++)
			    printf("  %le", v[j*3+k] - h[j][k]);
			printf("\n");
		    }
		    v += 9;
		    break;
	    }
	    //break;		/* only check basis 0 */
	}

	break;
      }	/* ForAllElements */
    }	/* g->nprocs == 1 */
#endif

    grad_u = phgDofGradient(u, NULL, NULL, "grad_u");
    grad_u0 = phgDofGradient(u0, NULL, NULL, "grad_u0");

    phgDofAXPY(-1.0, u, &u0);
    phgDofAXPY(-1.0, grad_u, &grad_u0);
    phgPrintf("\nApproximation errors:\n");
    phgPrintf("    # elements: %d, L2: %le, H1: %le\n", (int)g->nelem_global,
		(double)phgDofNormL2(u0), (double)phgDofNormL2(grad_u0));

    phgDofFree(&u);
    phgDofFree(&u0);
    phgDofFree(&grad_u);
    phgDofFree(&grad_u0);
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}

#endif	/* SPECHT_TEST */
