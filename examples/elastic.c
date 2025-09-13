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

/* $Id: elastic.c,v 1.54 2020/12/06 09:03:07 zlb Exp $ */

#include "phg.h"

#if (PHG_VERSION_MAJOR <= 0 && PHG_VERSION_MINOR < 9)
# undef ELEMENT
typedef SIMPLEX ELEMENT;
#endif

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef TEST
# define TEST	1	/* !0 => static (BVP), 0 => modal analysis (EVP) */
#endif

static FLOAT E  = 2.0/*e11*/;		/* elasticity modulus */
static FLOAT nu = 0.25;			/* Poisson ratio */
static FLOAT rho = 7850.0/*e-6*/;	/* density (kg/cm^3) */

/**********************************************************************

e[u_,v_,w_] := {D[u,x],D[v,y],D[w,z],D[v,z]+D[w,y],D[u,z]+D[w,x],D[u,y]+D[v,x]};

(* s[u_,v_,w_] := c{{1,a,a,0,0,0},{a,1,a,0,0,0},{a,a,1,0,0,0},
		    {0,0,0,b,0,0},{0,0,0,0,b,0},{0,0,0,0,0,b}} . e[u,v,w]; *)
sxx[u_,v_,w_]  := c {1,a,a,0,0,0} . e[u,v,w];
syy[u_,v_,w_]  := c {a,1,a,0,0,0} . e[u,v,w];
szz[u_,v_,w_]  := c {a,a,1,0,0,0} . e[u,v,w];
syz[u_,v_,w_] := c {0,0,0,b,0,0} . e[u,v,w];
szx[u_,v_,w_] := c {0,0,0,0,b,0} . e[u,v,w];
sxy[u_,v_,w_] := c {0,0,0,0,0,b} . e[u,v,w];

f[u_,v_,w_] := {D[sxx[u,v,w],x] + D[sxy[u,v,w],y] + D[szx[u,v,w],z],
		D[sxy[u,v,w],x] + D[syy[u,v,w],y] + D[syz[u,v,w],z],
		D[szx[u,v,w],x] + D[syz[u,v,w],y] + D[szz[u,v,w],z]};

CForm[f[1.*x*x*y*z, 2.*x*y*y*z, 3.*x*y*z*z]]

***********************************************************************/

#if TEST
#define func_g	func_u
static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *(value++) = 1. * x * x * y * z;
    *(value++) = 2. * x * y * y * z;
    *(value++) = 3. * x * y * z * z;
}
static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT a = nu / (1.0 - nu);
    FLOAT b = (1. - 2. * nu) / (2. * (1. - nu));
    FLOAT c = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu));

    *(value++) = 10.*b*c*y*z + c*(2.*y*z + 10.*a*y*z);
    *(value++) = 8.*b*c*x*z + c*(4.*x*z + 8.*a*x*z);
    *(value++) = 6.*b*c*x*y + c*(6.*x*y + 6.*a*x*y);
}
#else	/* TEST */
# define func_f	NULL
# define func_g	NULL
#endif	/* TEST */

static void
stiffness_matrix(ELEMENT *e, DOF *u, int m, int n, FLOAT E, FLOAT nu,
		 int lda, FLOAT *pA, QUAD *quad)
/* returns in 'A' the contribution of the m-th and n-th basis functions
 * to the element stiffness matrix. Elements of the 3x3 small matrix are
 * stored at A[j * lda + k] (j,k = 0,1,2) */
{
    int i;
    FLOAT (*A)[lda] = (void *)pA;
    const FLOAT *pm, *pn, *w;
    FLOAT a00, a01, a02, a10, a11, a12, a20, a21, a22;
    FLOAT a = nu / (1.0 - nu);
    FLOAT b = (1. - 2. * nu) / (2. * (1. - nu));
    FLOAT c = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu)) *
		phgGeomGetVolume(u->g, e);
 
    pm = phgQuadGetBasisGradient(e, u, m, quad);
    pn = phgQuadGetBasisGradient(e, u, n, quad);
    w = quad->weights;

    a00 = a01 = a02 = a10 = a11 = a12 = a20 = a21 = a22 = 0.;
    for (i = 0; i < quad->npoints; i++) {
	a00 += (pm[0] * pn[0] + b * (pm[1] * pn[1] + pm[2] * pn[2])) * *w;
	a01 += (a * pm[0] * pn[1] + b * pm[1] * pn[0]) * *w;
	a02 += (a * pm[0] * pn[2] + b * pm[2] * pn[0]) * *w;

	a10 += (a * pm[1] * pn[0] + b * pm[0] * pn[1]) * *w;
	a11 += (pm[1] * pn[1] + b * (pm[0] * pn[0] + pm[2] * pn[2])) * *w;
	a12 += (a * pm[1] * pn[2] + b * pm[2] * pn[1]) * *w;

	a20 += (a * pm[2] * pn[0] + b * pm[0] * pn[2]) * *w;
	a21 += (a * pm[2] * pn[1] + b * pm[1] * pn[2]) * *w;
	a22 += (pm[2] * pn[2] + b * (pm[0] * pn[0] + pm[1] * pn[1])) * *w;
	pn += 3;
	pm += 3;
	w++;
    }

    A[0][0] = a00 * c;
    A[0][1] = a01 * c;
    A[0][2] = a02 * c;

    A[1][0] = a10 * c;
    A[1][1] = a11 * c;
    A[1][2] = a12 * c;

    A[2][0] = a20 * c;
    A[2][1] = a21 * c;
    A[2][2] = a22 * c;

/*printf("a%d%d=[%lg %lg %lg; %lg %lg %lg; %lg %lg %lg];\n", m, n, a00*c, a01*c, a02*c, a10*c, a11*c, a12*c, a20*c, a21*c, a22*c);*/

    return;
}

static void
build_matrices(MAT *ma, MAT *mb, VEC *rhs_vec, DOF *u_h, DOF *f)
/* build stiffness (ma) and mass (mb) matrices */
{
    int N = u_h->type->nbas;	/* number of element basis functions */
    int i, j, k, l;
    GRID *g = u_h->g;
    ELEMENT *e;
    FLOAT A[N][3][N][3], buffer[N], rhs[N][3], d;
    static FLOAT *B0 = NULL;
    FLOAT B[N][N];
    INT I[N][3], J[3][N];
    QUAD *quad = phgQuadGetQuad3D((u_h->type->order - 1) * 2);

    assert(u_h->type->dim == 1 && u_h->dim == 3);
    assert(mb == NULL || u_h->type->invariant);

    if (mb != NULL && B0 == NULL && g->nroot > 0) {
	/* (\int \phi_j\cdot\phi_i)/vol is independent of element */
	FreeAtExit(B0);
	B0 = phgAlloc(N * N * sizeof(*B0));
	e = g->roots;
	d = 1. / phgGeomGetVolume(g, e);
	for (i = 0; i < N; i++)
	    for (j = 0; j <= i; j++)
		B0[i * N + j] = B0[i + j * N] = d *
		    phgQuadBasDotBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
    }

    ForAllElements(g, e) {
	d = phgGeomGetVolume(g, e) * rho;
	for (i = 0; i < N; i++) {
	    for (k = 0; k < 3; k++)
		J[k][i] = I[i][k] = phgMapE2L(ma->rmap, 0, e, i * 3 + k);
	    for (j = 0; j <= i; j++) {
		/* the stiffness matrix */
		stiffness_matrix(e, u_h, i, j, E, nu, 3 * N, &A[i][0][j][0],
				 quad);
		if (j != i) {
		    for (k = 0; k < 3; k++)
			for (l = 0; l < 3; l++)
			    A[j][k][i][l] = A[i][l][j][k];
		}
		/* the mass matrix: \int \phi_i\cdot\phi_j */
		if (mb != NULL)
		    B[j][i] = B[i][j] = B0[i * N + j] * d;
	    }
	}
#if 0
/* print element stiffness matrix */
FLOAT *p = (void *)A;
printf("vol = %0.16lg;\n", phgGeomGetVolume(g, e));
printf("a = [\n");
for (i = 0; i < 3 * N; i++) {
 for (j = 0; j < 3 * N; j++) printf("%0.16lg ", p[i * 3 * N + j]); printf("\n");
}
printf("];\n");
printf("b = [\n");
for (i = 0; i < N; i++) for (k = 0; k < 3; k++) {
 for (j = 0; j < N; j++) for (l = 0; l < 3; l++)
  printf("%0.16lg ", l == k ? B[i][j] : 0.0);
 printf("\n");
}
printf("];\n");
exit(0);
#endif

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, func_g, buffer, rhs[i],
					DOF_PROJ_NONE)) {
		/* Dirichlet boundary */
		for (k = 0; k < 3; k++) {
		    phgMatAddEntries(ma, 1, I[i] + k, N, J[k], buffer);
		    if (mb != NULL)
			phgMatAddEntries(mb, 1, I[i] + k, N, J[k], buffer);
		}
	    }
	    else {
		/* interior node */
		phgMatAddEntries(ma, 3, I[i], 3 * N, I[0], &A[i][0][0][0]);
		if (mb == NULL) {
		    phgQuadDofTimesBas(e, f, u_h, i, QUAD_DEFAULT, rhs[i]);
		    for (k = 0; k < 3; k++)
			rhs[i][k] = -rhs[i][k];
		}
		else {
		    for (k = 0; k < 3; k++)
			phgMatAddEntries(mb, 1, I[i] + k, N, J[k], B[i]);
		}
	    }
	}
	if (rhs_vec != NULL)
	    phgVecAddEntries(rhs_vec, 0, N * 3, I[0], rhs[0]);
    }
}

/* Note: the following functions/macros are to be moved to libphg.a */

/* Note: phgDofDisplacementToStress() works with variable E and nu,
 *	 phgDofDisplacementToStress0() works with constant E and nu */
#define phgDofDisplacementToStress(u, E, nu, stress_ptr, name) \
    phgDofDisplacementToStress_(u, E, nu, stress_ptr, name, __FILE__, __LINE__)

#define phgDofDisplacementToStress0(u, E, nu, stress_ptr, name) \
    phgDofDisplacementToStress0_(u, E, nu, stress_ptr, name, __FILE__, __LINE__)

/* static pointers for stress_func() */
static DOF *grad_u, *elasticity, *poisson;

static void
stress_func(DOF *str, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
/* TODO: cache bases of grad_u, elasticity, and poisson if they are invariant */
{
    static FLOAT E0 = -1., nu0 = -1., a, b, c;
    FLOAT E, nu, g[Dim][Dim], (*s)[Dim];

    phgDofEval(grad_u,     e, lambda, (FLOAT *)g);
    phgDofEval(elasticity, e, lambda, &E);
    phgDofEval(poisson,    e, lambda, &nu);

    /* Note: in most cases E and nu are piecewise constant, the following
     * test recomputes a, b, and c only when necessary */
    if (E != E0 || nu != nu0) {
	E0 = E;
	nu0 = nu;
	c = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu));
	a = nu / (1. - nu);
	b = c * (1. - 2. * nu) * .5 / (1. - nu);
    }

    s = (void *)values;
    s[0][0] = (g[0][0] + a * (g[1][1] + g[2][2])) * c;
    s[1][1] = (g[1][1] + a * (g[0][0] + g[2][2])) * c;
    s[2][2] = (g[2][2] + a * (g[0][0] + g[1][1])) * c;
    /* \tau_yz */
    s[1][2] = s[2][1] = b * (g[1][2] + g[2][1]);
    /* \tau_zx */
    s[0][2] = s[2][0] = b * (g[0][2] + g[2][0]);
    /* \tau_xy */
    s[0][1] = s[1][0] = b * (g[0][1] + g[1][0]);

    return;
}

DOF *
phgDofDisplacementToStress_(DOF *u, DOF *E, DOF *nu, DOF **stress_ptr,
		const char *name0, const char *srcfile, int srcline)
/* computes the stress tensor using constitution law for isotropic materials.
 * (see bibliography/Elasticity/constitution.pdf) */
{
    char *name;
    DOF *stress = (stress_ptr == NULL ? NULL : *stress_ptr);

    assert(DofDim(u) == 3);
    assert(DofDim(E) == 1 && DofDim(nu) == 1);

    MagicCheck(DOF, stress)

    if (stress != NULL)
	phgDofFree(&stress);

    /* set up DOF pointers for stress_func() */
    elasticity = E;
    poisson = nu;
    grad_u = phgDofGradient(u, NULL, NULL, NULL);

    /* create the new DOF */
    if (name0 != NULL) {
	name = phgAlloc(strlen(name0) + 1);
	strcpy(name, name0);
    }
    else {
	name = phgAlloc(strlen(u->name) + 1 + 8);
	sprintf(name, "stress(%s)", u->name);
    }
    stress = phgDofGetSameOrderDG_(grad_u, Dim * Dim, name, srcfile, srcline);
    phgFree(name);
    if (stress_ptr != NULL)
	*stress_ptr = stress;

    /* Note: may also use phgDofMM (less efficient?) */
    phgDofSetDataByLambdaFunction(stress, stress_func);

    phgDofFree(&grad_u);

    return stress;
}

DOF *
phgDofDisplacementToStress0_(DOF *u, FLOAT E, FLOAT nu, DOF **stress_ptr,
		const char *name0, const char *srcfile, int srcline)
/* code for constant E/nu */
{
    DOF *stress = (stress_ptr == NULL ? NULL : *stress_ptr);

#if 0
    /* Note: this code only works for constant E/nu and Lagrangian bases */
    char *name;
    FLOAT (*g)[3], (*s)[3], a, b, c;
    INT k, n;
    DOF *grad;

    assert(DofDim(u) == 3);

    MagicCheck(DOF, stress)

    if (stress != NULL)
	phgDofFree(&stress);

    if (u->type->grad_type->InitFunc != phgDofInitFuncPoint)
	phgError(1, "%s: not yet implemented for DOF type \"%s\"\n", __func__,
		 u->type->name);

    /* create the new DOF */
    if (name0 != NULL) {
	name = phgAlloc(strlen(name0) + 1);
	strcpy(name, name0);
    }
    else {
	name = phgAlloc(strlen(u->name) + 1 + 8);
	sprintf(name, "stress(%s)", u->name);
    }
    stress = phgDofNew_(u->g, u->type->grad_type, NULL,
			3 * 3, name, DofNoAction, srcfile, srcline);
    phgFree(name);

    if (stress_ptr != NULL)
	*stress_ptr = stress;

    if (u->data == NULL)
	return stress;

    grad = phgDofGradient(u, NULL, NULL, NULL);
    g = (void *)DofData(grad);
    s = (void *)DofData(stress);
    n = DofDataCount(stress);
    c = E * (1. - nu) / ((1. + nu) * (1. - 2. * nu));
    a = nu / (1. - nu);
    b = c * (1. - 2. * nu) * .5 / (1. - nu);
    for (k = 0; k < n; k += 3 * 3) {
	s[0][0] = (g[0][0] + a * (g[1][1] + g[2][2])) * c;
	s[1][1] = (g[1][1] + a * (g[0][0] + g[2][2])) * c;
	s[2][2] = (g[2][2] + a * (g[0][0] + g[1][1])) * c;
	/* \tau_yz */
	s[1][2] = s[2][1] = b * (g[1][2] + g[2][1]);
	/* \tau_zx */
	s[0][2] = s[2][0] = b * (g[0][2] + g[2][0]);
	/* \tau_xy */
	s[0][1] = s[1][0] = b * (g[0][1] + g[1][0]);
	g += 3;
	s += 3;
    }

    phgDofFree(&grad);
#else
    DOF *E1, *nu1;
    E1 = phgDofNew(u->g, DOF_CONSTANT, 1, "E", DofNoAction);
    phgDofSetDataByValue(E1, E);
    nu1 = phgDofNew(u->g, DOF_CONSTANT, 1, "nu", DofNoAction);
    phgDofSetDataByValue(nu1, nu);
    stress = phgDofDisplacementToStress_(u, E1, nu1, stress_ptr, name0,
					 srcfile, srcline);
    phgDofFree(&E1);
    phgDofFree(&nu1);
#endif

    return stress;
}

static FLOAT
#if TEST
estimate_error(DOF *f, DOF *u_h, DOF *error)
#else
estimate_error(FLOAT lambda, DOF *u_h, DOF *error)
#endif
/* compute H1 error indicator */
{
    GRID *g = u_h->g;
    ELEMENT *e;
    DOF *stress, *jump, *residual, *tmp;

    stress = phgDofDisplacementToStress0(u_h, E, nu, NULL, NULL);
    tmp = phgDofDivergence(stress, NULL, NULL, NULL);
    residual = phgDofGetSameOrderDG(u_h, -1, NULL);
#if TEST
    phgDofCopy(f, &residual, NULL, NULL);
    phgDofAXPBY(1.0, tmp, -1.0, &residual);
#else	/* TEST */
    phgDofCopy(u_h, &residual, NULL, NULL);
    phgDofAXPBY(1.0, tmp, lambda, &residual);
#endif	/* TEST */
    phgDofFree(&tmp);
    jump = phgQuadFaceJump(stress, DOF_PROJ_DOT, NULL, QUAD_DEFAULT);
    ForAllElements(g, e) {
	int i;
	FLOAT eta, h;
	FLOAT diam = phgGeomGetDiameter(g, e);
	eta = 0.0;
	/* for each face F compute [stress \cdot n] */
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] & (DIRICHLET | NEUMANN))
		continue;	/* boundary face */
	    h = phgGeomGetFaceDiameter(g, e, i);
	    eta += *DofFaceData(jump, e->faces[i]) * h;
	}
	eta = eta*.5 + diam*diam * phgQuadDofDotDof(e, residual, residual, -1);
	*DofElementData(error, e->index) = eta;
    }
    phgDofFree(&jump);
    phgDofFree(&residual);
    phgDofFree(&stress);

    return Sqrt(phgDofNormInftyVec(error));
}

static int
bc_map(int bctype)
{
    return DIRICHLET;	/* set Dirichlet BC on all boundaries */
}

int
main(int argc, char *argv[])
{
    //static char *fn = "../test/model.mesh";
    static char *fn = "../test/cube.dat";
    size_t mem, mem_peak;
    GRID *g;
    ELEMENT *e;
#if TEST
    FLOAT tol = 1e-2;
    DOF *u_h, *error, *f, *tmp;
    SOLVER *solver;
#else	/* TEST */
    FLOAT tol = 1e-1;
    int i, n, nit;
    DOF **u_h, **error;
    MAT *B;
    FLOAT *evals;
#endif	/* TEST */
    INT nev = 20, pre_refines = 0;
    MAP *map;
    MAT *A;
    FLOAT PEoo, thres;
    double wtime;

    phgOptionsRegisterFloat("tol", "Convergence criterion", &tol);
    phgOptionsRegisterInt("nev", "Number of eigenvalues", &nev);
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);

    phgInit(&argc, &argv);
    if (argc == 2)
	fn = argv[1];
    g = phgNewGrid(-1);
    phgImportSetBdryMapFunc(bc_map);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    /* pre-refinement */
    phgRefineAllElements(g, pre_refines);

#if TEST
    u_h = phgDofNew(g, DOF_DEFAULT, 3, "u_h", DofInterpolation);
    error = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);
    f = phgDofNew(g, DOF_ANALYTIC, 3, "f", func_f);
#else	/* TEST */
    u_h = phgAlloc(nev * sizeof(*u_h));
    error = phgAlloc(nev * sizeof(*error));
    evals = phgCalloc(nev, sizeof(*evals));
    for (i = 0; i < nev; i++) {
	u_h[i] = phgDofNew(g, DOF_DEFAULT, 3, "u_h", DofInterpolation);
	phgDofSetDirichletBoundaryMask(u_h[i], BDRY_MASK);
	/* set random initial values for the eigenvectors */
	phgDofRandomize(u_h[i], i == 0 ? 123 : 0);
	error[i] = phgDofNew(g, DOF_P0, 1, "error indicator", DofNoAction);
    }
#endif	/* TEST */

    while (TRUE) {
	phgPrintf("\n");
	if (phgBalanceGrid(g, 1.2, -1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",(double)g->lif);
#if TEST
	phgPrintf("%"dFMT" DOF, %"dFMT
		  " elements, %d submesh%s load imbalance: %lg\n",
			DofDataCountGlobal(u_h), g->nleaf_global,
			g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);
	map = phgMapCreate(u_h, NULL);
	A = phgMapCreateMat(map, map);
	A->handle_bdry_eqns = TRUE;
	phgMapDestroy(&map);
	solver = phgMat2Solver(SOLVER_DEFAULT, A);
	phgVecDisassemble(solver->rhs);
	phgMatDestroy(&A);
	build_matrices(solver->mat, NULL, solver->rhs, u_h, f);
	wtime = phgGetTime(NULL);
	phgSolverSolve(solver, TRUE, u_h, NULL);
	phgPrintf("  %d iterations, residual: %0.6le, wtime: %0.4lfs\n",
			solver->nits, (double)solver->residual,
			phgGetTime(NULL) - wtime);
	phgSolverDestroy(&solver);
	PEoo = estimate_error(f, u_h, error);
	tmp = phgDofNew(g, u_h->type, 3, "u", func_u);
	phgDofAXPY(-1.0, u_h, &tmp);
	phgPrintf("  L2(u-u_h) = %0.6le, Loo(u-u_h) = %0.6le\n",
		  (double)phgDofNormL2(tmp), (double)phgDofNormInftyVec(tmp));
	phgDofFree(&tmp);
#else	/* TEST */
	phgPrintf("%"dFMT" DOF, %"dFMT
		  " elements, %d submesh%s load imbalance: %lg\n",
			DofDataCountGlobal(u_h[0]), g->nleaf_global,
			g->nprocs, g->nprocs > 1 ? "es" : "", (double)g->lif);
	map = phgMapCreate(*u_h, NULL);
	{
	    MAP *m = phgMapRemoveBoundaryEntries(map);
	    A = phgMapCreateMat(m, m);
	    B = phgMapCreateMat(m, m);
	    phgMapDestroy(&m);
	}
	build_matrices(A, B, NULL, u_h[0], NULL);
	wtime = phgGetTime(NULL);
	n = phgDofEigenSolve(A, B, nev, EIGEN_SMALLEST, 0.0, &nit, evals,
				map, u_h, NULL);
	phgMapDestroy(&map);
	phgPrintf("  %d its, converged eigenvalues: %d, wtime: %0.4lfs\n",
			nit, n, phgGetTime(NULL) - wtime);
	for (i = 0; i < n; i++)
	    phgPrintf("    tau[%d]=%0.6le, L2(u_h[%d])=%0.6le\n",
		i, (double)evals[i], i, (double)phgDofNormL2(u_h[i]));
	for (; i < nev; i++)
	    phgDofRandomize(u_h[i], 0);
	phgMatDestroy(&A);
	phgMatDestroy(&B);

	PEoo = -1.0;
	for (i = 0; i < n; i++) {
	    thres = estimate_error(evals[i], u_h[i], error[i]);
	    if (PEoo < thres)
		PEoo = thres;
	}
	PEoo = Fabs(PEoo);
#endif	/* TEST */
	phgPrintf("  indicator=%0.6le\n", (double)PEoo);
	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("  Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		(double)mem / (1024.0 * 1024.0),
		(double)mem_peak / (1024.0 * 1024.0));
	if (PEoo < tol || mem_peak >= 1024 * (size_t)1024 * 400)
	    break;
	/*thres = tol * tol / (double)g->nleaf_global;*/
	thres = Pow(PEoo * 0.5, 2);
	ForAllElements(g, e)
	    e->mark = 0;
#if TEST
	ForAllElements(g, e)
	    if (*DofElementData(error, e->index) > thres)
		e->mark = 1;
#else	/* TEST */
	if (n == 0) {
	    ForAllElements(g, e)
		e->mark = 1;
	}
	else {
	    for (i = 0; i < nev; i++)
		ForAllElements(g, e)
		    if (*DofElementData(error[i],e->index) > thres)
			e->mark = 1;
	}
#endif	/* TEST */
	phgRefineMarkedElements(g);
    }

#if 1
# if TEST
    phgPrintf("\nCreating \"%s\".\n", phgExportVTK(g, "output.vtk", u_h,NULL));
# else	/* TEST */
    phgPrintf("\nCreating \"%s\".\n", phgExportVTKn(g, "output.vtk",nev,u_h));
# endif	/* TEST */
#endif	/* 0|1 */

#if TEST
    phgDofFree(&f);
    phgDofFree(&u_h);
    phgDofFree(&error);
#else	/* TEST */
    phgFree(evals);
    for (i = 0; i < nev; i++) {
	phgDofFree(u_h + i);
	phgDofFree(error + i);
    }
    phgFree(u_h);
    phgFree(error);
#endif	/* TEST */
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
