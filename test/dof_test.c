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

/* This code tests basic DOF operations (including basis functions).
 * $Id: dof_test.c,v 1.94 2021/11/23 04:55:39 zlb Exp $ */

#include "phg.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef enum {D0, D1, D2} WHAT;

static INT func_order = -1;

static void
func0(FLOAT x, FLOAT y, FLOAT z, FLOAT *values, WHAT what)
{
    static int order = -1;
    int i;
    int dim = DOF_DEFAULT->dim;
    FLOAT a, b, c;

    if (order < 0) {
	if ((order = func_order) < 0)
	    order = (DOF_DEFAULT->mass_lumping != NULL ?
			DOF_DEFAULT->mass_lumping->order0 :
			DOF_DEFAULT->order);
	phgPrintf("Order of test function = %d\n", order);
    }

    if (order == 0 || DOF_DEFAULT == DOF_HC0 || DOF_DEFAULT == DOF_ND1) {
	/* constant */
	switch (what) {
	    case D0:
		if (dim == 1) {
		    *values = Asin((FLOAT)1.0) * 4.0;
		}
		else if (DOF_DEFAULT == DOF_HC0 || DOF_DEFAULT == DOF_ND1) {
		    /* (1,2,3) x (x,y,z) + (4,5,6) */
		    *(values++) = 2. * z - 3. * y + 4.;
		    *(values++) = 3. * x - 1. * z + 5.;
		    *(values) = 1. * y - 2. * x + 6.;
		}
		else {
		    *(values++) = 1.0;
		    *(values++) = 2.0;
		    *values = 3.0;
		}
		break;
	    case D1:
		if (dim == 1) {
		    values[0] = values[1] = values[2] = 0.;
		}
		else if (DOF_DEFAULT == DOF_HC0 || DOF_DEFAULT == DOF_ND1) {
		    values[0] = 2.;
		    values[1] = 4.;
		    values[2] = 6.;
		}
		else {
		    for (i = 0; i < dim * 3; i++)
			values[i] = 0.;
		}
		break;
	    case D2:
		if (dim == 1) {
		    *values = 0.;
		}
		else {
		    *values = 0.;
		}
		break;
	}
	return;
    }

    switch (what) {
	case D0:
	    a = b = c = 1.;
	    for (i = 0; i < order; i++) {
		a *= x;
		b *= y;
		c *= z;
	    }
	    for (i = 0; i < dim; i++)
    		*(values++) = (i + 1) * (a - b - b + c);
	    break;
	case D1:
	    a = b = c = order;
	    for (i = 0; i < order - 1; i++) {
		a *= x;
		b *= y;
		c *= z;
	    }
	    for (i = 0; i < dim; i++, values += 3) {
		values[0] = (i + 1) * a;
		values[1] = (i + 1) * (- b - b);
		values[2] = (i + 1) * c;
	    }
	    break;
	case D2:
	    a = b = c = order * (order - 1);
	    for (i = 0; i < order - 2; i++) {
		a *= x;
		b *= y;
		c *= z;
	    }
	    for (i = 0; i < dim; i++)
		*(values++) = (i + 1) * (a - b - b + c);
	    break;
    }

    return;
}

static void
func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    func0(x, y, z, value, D0);
}

static void
d1_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    func0(x, y, z, values, D1);
}

static void
d2_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    func0(x, y, z, value, D2);
}

static DOF *dummy_dof;
static void
dummy_func(DOF *dof, ELEMENT *e, int bno, const FLOAT *lambda, FLOAT *values)
{
     phgDofEval_(dummy_dof, e, lambda, NULL, 0, values, NULL);
}

int
main(int argc, char **argv)
{
    GRID *g;
    DOF *u, *d1_u, *d2_u, *tmp;
    INT i, pre_refines = 1, post_refines = 1;
    BOOLEAN dump_dof = FALSE, test_bas = TRUE;
    FLOAT err, nrm;
    char *fn = "cube.dat";
    const char *dname;

    FunctionEntry;

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterInt("-pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("-post_refines", "Post-refines", &post_refines);
    phgOptionsRegisterInt("-order", "Order of the test function", &func_order);
    phgOptionsRegisterNoArg("-dump", "Dump DOF values", &dump_dof);
    phgOptionsRegisterNoArg("-test_bas", "Test bases evaluation", &test_bas);

    phgOptionsPreset("+warn_nonopt");

    phgInit(&argc, &argv);
    /*phgPause(0);*/

    if (argc != 2) {
	phgPrintf("Usage: %s dof_type\n", argv[0]);
	phgFinalize();
	exit(1);
    }

    phgOptionsSetHandler("-dof_type", argv[1]);

    g = phgNewGrid(-1);
    phgImport(g, fn, FALSE);

#if 0
    if (FALSE && DOF_DEFAULT->dim == 3) {
	/* Test with a single basis function */
	FLOAT v[3];
	ELEMENT *e = g->roots;
	FLOAT lambda[4] = {0.15, 0.2, 0.3, 0.35};

	tmp = phgDofNew(g, DOF_ANALYTIC, 1, "u", /*basis_*/func);
	u = phgDofNew(g, DOF_DEFAULT, 1, "u", /*basis_*/func);
	phgDofDump(u);
	phgDofEval(tmp, e, lambda, v);
	printf("analytic = %lg %lg %lg\n", (double)v[0], (double)v[1], (double)v[2]);
/*----------------------------------------------------------------------------*/
if (FALSE) {	/* test with L2 projection on the whole element */
    FLOAT w0, x, y, z, *fvals = v;
    int j, k, n = u->type->nbas;
    FLOAT A[n * n], B[n];
    const FLOAT *qp, *qw, *p, *q, *b;
    QUAD *quad = phgQuadGetQuad3D(6);
    qp = quad->points;
    qw = quad->weights;
    e = g->roots;
    bzero(A, sizeof(A));
    bzero(B, sizeof(B));
    for (k = 0; k < quad->npoints; k++, qp += 4) {
	phgGeomLambda2XYZ(g, e, qp, &x, &y, &z);
	func(x, y, z, fvals);
	b = DOF_DEFAULT->BasFuncs(u, e, 0, -1, qp);
	/* coefficients and RHS */
	w0 = *(qw++);
	for (i = 0, p = b; i < n; i++, p += 3) {
	    for (j = i, q = p; j < n; j++, q += 3)
		A[i * n + j] += (p[0]*q[0] + p[1]*q[1] + p[2]*q[2]) * w0;
	    q = fvals;
	    B[i] += (p[0]*q[0] + p[1]*q[1] + p[2]*q[2]) * w0;
	}
    }
    for (i = 1; i < n; i++)
	for (j = 0; j < i; j++)
	    A[i * n + j] = A[j * n + i];
    /* solver the linear system of equations */
    if (!phgSolverDenseSolver(n, 1, (void *)A, (void *)B))
	phgError(1, "%s:%d, singular system (shouldn't happen).\n",
			__FILE__, __LINE__);
    memcpy(u->data, B, sizeof(B));
    phgDofDump(u);
}
/*----------------------------------------------------------------------------*/
	phgDofEval(u, e, lambda, v);
	printf("interpolated = %lg %lg %lg\n", (double)v[0], (double)v[1], (double)v[2]);
	phgDofFree(&u);
	phgDofFree(&tmp);
	exit(0);
    }
#endif

    /* Testing BasFuncs and Init */
    if (!test_bas)
	goto skip1;
    dummy_dof = tmp = phgDofNew(g, DOF_DEFAULT, 1, NULL, DofNoAction);
    u = phgDofNew(g, tmp->type, 1, NULL, DofNoAction);
    if (DofDataCountGlobal(tmp) >= 500)
	goto skip;
    phgPrintf("========================================================\n");
    phgPrintf("Testing \"%s\": BasFuncs/InitFunc.\n", tmp->type->name);
    for (i = 0; i < DofDataCount(tmp); i++) {
	INT j;
	phgPrintf("Testing basis no. %d\r", i);
	tmp->data[i] = 1.0;
	phgDofSetDataByFunction_(u, NULL, dummy_func, NULL);
	for (j = 0; j < DofDataCount(u); j++) {
	    if (Fabs(u->data[j] - tmp->data[j]) > 1000.0 * FLOAT_EPSILON) {
		phgDofDump(tmp); phgDofDump(u);
		phgPrintf("\nerror: wrong data at i=%d, j=%d: %lg ==> %lg\n",
				i, j, (double)tmp->data[j], (double)u->data[j]);
		phgError(1, "abort.\n");
	    }
	}
	if (j < DofDataCount(u))
	    break;
	tmp->data[i] = 0.0;
    }
    if (i >= DofDataCount(u))
	phgPrintf("\npassed\n");
skip:
    phgDofFree(&u);

    /* Testing calling individual BasFuncs and BasGrads */
    phgPrintf("========================================================\n");
    phgPrintf("Testing \"%s\": evaluating individual bases.\n",
		tmp->type->name);
    {
	const FLOAT lambda[] = {1./29., 5./29., 9./29., 14./29.}, *p;
	FLOAT values[(Dim + 1) * tmp->type->dim];
	ELEMENT *e = g->roots;
	int ii, n;
	for (i = 0; i < tmp->type->nbas; i++) {
	    phgPrintf("Testing basis no. %d\r", i);
	    /* BasFuncs */
	    n = tmp->type->dim;
	    memcpy(values, tmp->type->BasFuncs(tmp, e, i, i + 1, lambda),
		   n * sizeof(*values));
	    p = tmp->type->BasFuncs(tmp, e, 0, -1, lambda) + i * n;
	    for (ii = 0; ii < n; ii++) {
		if (Fabs(p[ii] - values[ii]) > 100. * FLOAT_EPSILON)
		    phgError(1, "(%s:%d) test failed, error = %le.\n", __FILE__, __LINE__, (double)(p[ii] - values[ii]));
	    }
	    /* BasGrads */
	    n = (Dim + 1) * tmp->type->dim;
	    memcpy(values, tmp->type->BasGrads(tmp, e, i, i + 1, lambda),
		   n * sizeof(*values));
	    p = tmp->type->BasGrads(tmp, e, 0, -1, lambda) + i * n;
	    for (ii = 0; ii < n; ii++) {
		if (Fabs(p[ii] - values[ii]) > 100. * FLOAT_EPSILON)
		    phgError(1, "(%s:%d) test failed.\n", __FILE__, __LINE__);
	    }
	}
	phgPrintf("\npassed\n");
    }

    phgDofFree(&tmp);

skip1:
    phgPrintf("Refining grid ...\n");
    phgRefineAllElements(g, pre_refines);
    phgPartitionGrid(g);

    phgPrintf("========================================================\n");
    phgPrintf("Testing \"%s\": DOF operations.\n", DOF_DEFAULT->name);
    u = phgDofNew(g, DOF_DEFAULT, 1, "u", DofInterpolation);
    phgDofSetDataByFunction(u, func);
    if (dump_dof)
	phgDofDump(u);

    PrintTime(0);

    for (i = 0; i < post_refines
			&& DofDataCountGlobal(u) / g->nprocs < 500000; i++) {
	/*phgRefineRandomElements(g, "25%");
	phgRefineSurface(g);*/
	phgRefineAllElements(g, 1);
	phgRedistributeGrid(g);
    }
    /*phgCheckConformity(g);*/

    PrintTime(0);

    /* check results */
#if 0
    /* skip checking derivatives */
    nrm = phgDofNormInftyVec(u);
    tmp = phgDofNew(g, DOF_ANALYTIC, DofDim(u), "analytic", func);
    phgDofAXPY(-1.0, tmp, &u);
    err = phgDofNormInftyVec(u);
    phgPrintf("Interpolation error  = %0.3le (norm = %0.3le).\n", err, nrm);
    phgDofFree(&tmp);
    phgDofFree(&u);
    phgFinalize();
    exit(0);
#else
    if (DOF_DEFAULT != DOF_HC0 && DOF_DEFAULT != DOF_ND1) {
	d1_u = phgDofGradient(u, NULL, NULL, "grad_u");
	d2_u = phgDofDivergence(d1_u, NULL, NULL, "div_grad_u");
	dname = "Gradient";
    }
    else {
	d1_u = phgDofCurl(u, NULL, NULL, "curl_u");
	d2_u = phgDofCurl(d1_u, NULL, NULL, "curl_curl_u");
	dname = "Curl";
    }
    nrm = phgDofNormInftyVec(u);
    tmp = phgDofNew(g, DOF_ANALYTIC, DofDim(u), "analytic", func);
    phgDofAXPY(-1.0, tmp, &u);
    err = phgDofNormInftyVec(u);
    phgPrintf("Interpolation error  = %0.3le (norm = %0.3le).\n",
		(double)err, (double)nrm);
    phgDofFree(&tmp);
    phgDofFree(&u);
#endif

    /* checking gradient */
    nrm = phgDofNormInftyVec(d1_u);
    tmp = phgDofNew(g, DOF_ANALYTIC, DofDim(d1_u), "d1_u0", d1_func);
    phgDofAXPY(-1.0, tmp, &d1_u);
    err = phgDofNormInftyVec(d1_u);
    phgPrintf("%s (%s) error  = %0.3le (norm = %0.3le).\n", dname,
			d1_u->type->name, (double)err, (double)nrm);
    phgDofFree(&d1_u);
    phgDofFree(&tmp);

    /* checking divergence */
    nrm = phgDofNormInftyVec(d2_u);
    tmp = phgDofNew(g, DOF_ANALYTIC, DofDim(d2_u), "d2_u0", d2_func);
    phgDofAXPY(-1.0, tmp, &d2_u);
    err = phgDofNormInftyVec(d2_u);
    phgPrintf("%s (%s) error = %0.3le (norm = %0.3le).\n",
			DOF_DEFAULT == DOF_HC0 ? "CurlCurl" : "Laplacian",
			d2_u->type->name, (double)err, (double)nrm);
    phgDofFree(&d2_u);
    phgDofFree(&tmp);

#if 0
    /* Testing DofAXPY / DofCopy */
    DOF_DEFAULT = DOF_HC2;
    phgPrintf("==============================================\n");
    u = phgDofNew(g, DOF_P2, 3, NULL, DofNoAction);
    tmp = phgDofNew(g, DOF_DEFAULT, 1, NULL, DofNoAction);
    phgDofSetDataByFunction(u, func);
    phgDofSetDataByFunction(tmp, func);
    phgDofAXPY(-1.0, u, &tmp);
    phgPrintf("%s ==> %s: error = %lf\n", u->type->name, tmp->type->name,
		(double)phgDofNormInftyVec(tmp));
    phgDofFree(&u);
    phgDofFree(&tmp);

    tmp = phgDofNew(g, DOF_P2, 3, NULL, DofNoAction);
    u = phgDofNew(g, DOF_DEFAULT, 1, NULL, DofNoAction);
    phgDofSetDataByFunction(u, func);
    phgDofSetDataByFunction(tmp, func);
    phgDofAXPBY_(-1.0, u, 1.0, &tmp, __FILE__, __LINE__, FALSE);
    phgPrintf("%s ==> %s: error = %lf\n", u->type->name, tmp->type->name,
		(double)phgDofNormInftyVec(tmp));
    phgDofFree(&u);
    phgDofFree(&tmp);

    phgPrintf("==============================================\n");
    u = phgDofNew(g, DOF_DG2, 3, NULL, DofNoAction);
    tmp = phgDofNew(g, DOF_DEFAULT, 1, NULL, DofNoAction);
    phgDofSetDataByFunction(u, func);
    phgDofSetDataByFunction(tmp, func);
    phgDofAXPBY_(-1.0, u, 1.0, &tmp, __FILE__, __LINE__, FALSE);
    phgPrintf("%s ==> %s: error = %lf\n", u->type->name, tmp->type->name,
		(double)phgDofNormInftyVec(tmp));
    phgDofFree(&u);
    phgDofFree(&tmp);

    tmp = phgDofNew(g, DOF_DG2, 3, NULL, DofNoAction);
    u = phgDofNew(g, DOF_DEFAULT, 1, NULL, DofNoAction);
    phgDofSetDataByFunction(u, func);
    phgDofSetDataByFunction(tmp, func);
    phgDofAXPY(-1.0, u, &tmp);
    phgPrintf("%s ==> %s: error = %lf\n", u->type->name, tmp->type->name,
		(double)phgDofNormInftyVec(tmp));
    phgDofFree(&u);
    phgDofFree(&tmp);

    phgPrintf("==============================================\n");

    DOF_DEFAULT = DOF_P4;

    tmp = phgDofNew(g, DOF_DG4, 1, NULL, DofNoAction);
    u = phgDofNew(g, DOF_P4, 1, NULL, DofNoAction);
    phgDofSetDataByFunction(u, func);
    phgDofSetDataByFunction(tmp, func);
    phgDofAXPY(-1.0, u, &tmp);
    phgPrintf("%s ==> %s: error = %lf\n", u->type->name, tmp->type->name,
		(double)phgDofNormInftyVec(tmp));
    phgDofFree(&u);
    phgDofFree(&tmp);

    u = phgDofNew(g, DOF_DG4, 1, NULL, DofNoAction);
    tmp = phgDofNew(g, DOF_P4, 1, NULL, DofNoAction);
    phgDofSetDataByFunction(u, func);
    phgDofSetDataByFunction(tmp, func);
    phgDofAXPBY_(-1.0, u, 1.0, &tmp, __FILE__, __LINE__, FALSE);
    phgPrintf("%s ==> %s: error = %lf\n", u->type->name, tmp->type->name,
		(double)phgDofNormInftyVec(tmp));
    phgDofFree(&u);
    phgDofFree(&tmp);

    u = phgDofNew(g, DOF_DG4, 1, NULL, DofNoAction);
    tmp = phgDofNew(g, DOF_DG2, 1, NULL, DofNoAction);
    phgDofSetDataByFunction(u, func);
    phgDofSetDataByFunction(tmp, func);
    phgDofAXPBY_(-1.0, u, 1.0, &tmp, __FILE__, __LINE__, FALSE);
    phgPrintf("%s ==> %s: error = %lf\n", u->type->name, tmp->type->name,
		(double)phgDofNormInftyVec(tmp));
    phgDofFree(&u);
    phgDofFree(&tmp);

    tmp = phgDofNew(g, DOF_DG4, 1, NULL, DofNoAction);
    u = phgDofNew(g, DOF_DG2, 1, NULL, DofNoAction);
    phgDofSetDataByFunction(u, func);
    phgDofSetDataByFunction(tmp, func);
    phgDofAXPY(-1.0, u, &tmp);
    phgPrintf("Note: the next result may vary with different nprocs due to "
	      "random refinement\n");
    phgPrintf("%s ==> %s: error = %lf\n", u->type->name, tmp->type->name,
		(double)phgDofNormInftyVec(tmp));
    phgDofFree(&u);
    phgDofFree(&tmp);
#endif

    phgFreeGrid(&g);
    phgFinalize();

    Return 0;
}
