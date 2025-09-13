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

/* $Id: quad_test2.c,v 1.208 2022/01/12 13:08:22 zlb Exp $
 *
 * Test program for phgQuadInterface() */

#include "phg.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifndef LS_ORDER
# define LS_ORDER	LS_ORDER0
#endif

#ifndef CASE
/* The test cases (make USER_CFLAGS="-DCASE=0")
 *	0: sphere
 *	1: plane
 *	2: cylinder
 *	3: polynomial, non convex
 *	4: two spheres (only works with some meshes, e.g., cube.dat with 3+
 *	   refinements)
 *	5: non polynomial, non convex
 */
# define CASE	0
#endif

#if !(CASE >= 0 && CASE <= 5)
# warning invalid value for CASE, using CASE=0
# undef CASE
# define CASE 0
#endif

#define PI	(4.0 * Atan(1.0))

static void
levelset_func(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
#if CASE == 0	/*-------------------------------------------------------------------*/
# define DESC		"a sphere of radius 1/4 centered at (1/2,1/2,1/2)"
# define LS_ORDER0	2
# define VOL		_F(4.)/_F(3.) * PI * Pow(0.25,3)
# define ARE		4. * PI * Pow(0.25,2)
# define NORMAL		v[0]=2.*(x-.5); v[1]=2.*(y-.5); v[2]=2.*(z-.5);
    x -= 0.5; y -= 0.5; z -= 0.5;
    *value = x * x + y * y + z * z - 0.25 * 0.25;
#elif CASE == 1	/*-------------------------------------------------------------------*/
# define DESC		"a plane at z=0.5"
# define LS_ORDER0	1
# define VOL		0.5
# define ARE		1.0
# define NORMAL		v[0]=0.; v[1]=0.; v[2]=1.;
    *value = z - 0.5;
#elif CASE == 2	/*-------------------------------------------------------------------*/
# define DESC		"a cylinder of radius 0.25 centered at x=z=0.5"
# define LS_ORDER0	2
# define VOL		PI * Pow(0.25,2)
# define ARE		0.5 * PI
# define NORMAL		v[0]=2.*(x-.5); v[1]=0.0; v[2]=2.*(z-.5);
    x -= 0.5; z -= 0.5;
    *value = x * x + z * z - 0.25 * 0.25;
#elif CASE == 3	/*-------------------------------------------------------------------*/
# define DESC		"non-convex, higher order polynomial"
# define LS_ORDER0	4
# define VOL		_F(.51981836309785254023043979577190664)
# define ARE		_F(2.2228899780473178633613776771549132)
# define NORMAL	\
    v[0] = 4.*Pow(x-.5,1); v[1] = -24.*Pow(y-.5,2); v[2] = -64.*Pow(z-.5,3);
    *value = 2.0*Pow(x-.5,2)-8.0*Pow(y-.5,3)-16.0*Pow(z-.5, 4) - 0.02;
#elif CASE == 4	/*-------------------------------------------------------------------*/
/* WARNING: won't work with cube5.dat! (violating alignment requirement) */
# define DESC		"two spheres, r=1/4, c=(1/2,1/2,1/4),(1/2,1/2,3/4)"
# define LS_ORDER0	2
# define VOL		2.0 * 4./3. * 4.*Atan(1.0) * Pow(0.25,3)
# define ARE		2.0 * 4.0 * PI * Pow(0.25,2)
# define NORMAL \
	if (z <= .5) {v[0]=2.*(x-.5); v[1]=2.*(y-.5); v[2]=2.*(z-.25);} \
	else {v[0]=2.*(x-.5); v[1]=2.*(y-.5); v[2]=2.*(z-.75);}
    FLOAT d1, d2;
    x -= 0.5; y -= 0.5; z -= 0.25;
    d1 = x * x + y * y + z * z - 0.25 * 0.25;
    z -= 0.5;
    d2 = x * x + y * y + z * z - 0.25 * 0.25;
    *value = (d1 < d2 ? d1 : d2);
#elif CASE == 5	/*-------------------------------------------------------------------*/
/* Note: analytic results not all known, results below are for LS_ORDER=5 */
# define DESC		"non-polynomial level set function"
# define LS_ORDER0	-2	/* non-polynomial */
# define VOL		_F(.49999941858294764204999750286803778)
# define ARE		_F(4.0234421230249520358971520380217822)
# define NORMAL \
	x = 9.5 * (x-0.5); y = 9.5 * (y-0.5); z = 4.75 * (z-0.5); \
	v[0]=(-Sin(x)*Sin(y) + Cos(z)*Cos(x)) * 9.5; \
	v[1]=( Cos(x)*Cos(y) - Sin(y)*Sin(z)) * 9.5; \
	v[2]=( Cos(y)*Cos(z) - Sin(z)*Sin(x)) * 4.75;
    x = 9.5 * (x-0.5); y = 9.5 * (y-0.5); z = 4.75 * (z-0.5);
    *value = Cos(x)*Sin(y) + Cos(y)*Sin(z) + Cos(z)*Sin(x);
#endif
}

static void
levelset_grad(FLOAT x, FLOAT y, FLOAT z, FLOAT *v)
{
    NORMAL
}

static void
one(FLOAT x, FLOAT y, FLOAT z, FLOAT *v)
{
    *v = 1.0;
}

static void normal(FLOAT x, FLOAT y, FLOAT z, FLOAT *v)
{
    FLOAT d;
    levelset_grad(x, y, z, v);
    d = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if (d > 0.0)
	d = 1.0 / Sqrt(d);
    v[0] *= d; v[1] *= d; v[2] *= d; 
}

static void
refine_interface(DOF *ls, DOF *ls_grad, int depth)
/* refine elements on the interface (needn't be very thorough) */
{
    GRID *g = ls->g;
    ELEMENT *e;
    INT nelem = g->nelem_global;
    double t0;
    size_t mem;
    static int level = 0;

    t0 = phgGetTime(NULL);
#if 0
    ForAllElements(g, e) {
	FLOAT d, lam[] = {0., 0., 0., 0.};
	int i, sign = 0;
	for (i = 0; i < NVert; i++) {
	    lam[i] = 1.0;
	    phgDofEval(ls, e, lam, &d);
	    lam[i] = 0.0;
	    if (Fabs(d) <= FLOAT_EPSILON)
		break;
	    if (i == 0) {
		sign = d > 0.0 ? 1 : -1;
		continue;
	    }
	    if (sign * d < 0.0)
		break;
	}
	e->mark = i < NVert ? depth : 0;
    }
#else
    phgQuadInterfaceMarkGrid(ls, ls_grad);
    ForAllElements(g, e)
	e->mark = (e->mark == 0 ? depth : 0);
#endif

    phgRefineMarkedElements(g);

    if (nelem == g->nelem_global)
	phgRefineAllElements(g, depth);

    phgBalanceGrid(g, 1.0, 0, NULL, 0.);

    t0 = phgGetTime(NULL) - t0;
    phgMemoryUsage(g, &mem);
    phgPrintf("--- refine level %d, elem %"dFMT", mem %0.2lgMB, time %0.4lfs\n",
				level += depth, g->nelem_global,
				mem / (1024.0*1024.0), t0);
}

#ifndef REFINE
# define REFINE 0
#endif

static int refine = REFINE, refine1 = REFINE, refine_step = 3;
static int order = 1, order1 = 1, order_step = 2;

static const char *
opts_func(OPTION *o, const char *arg, int prefix)
{
    char *end;
    int *p, *p1, *p2;

    Unused(prefix);

    if (arg == NULL || arg[0] == '\0')
	return NULL;

    if (!strcmp(o->name, "refine")) {
	p  = &refine;
	p1 = &refine1;
	p2 = &refine_step;
	*p2 = 3;
    }
    else {
	assert(!strcmp(o->name, "order"));
	p  = &order;
	p1 = &order1;
	p2 = &order_step;
	*p2 = 2;
    }
    *p1 = *p = strtol(arg, &end, 10);
    if (end == arg || end == NULL ||
	(*end != '\0' && *end != ':' && *end != ','))
	return NULL;	/* signal error */
    if (*end != '\0') {
	arg = end + 1;
	*p1 = strtol(arg, &end, 10);
	if (end == arg || end == NULL ||
	    (*end != '\0' && *end != ':' && *end != ','))
	    return NULL;
	if (*end != '\0') {
	    arg = end + 1;
	    *p2 = strtol(arg, &end, 10);
	    if (arg == end || end == NULL || *end != '\0')
		return NULL;
	}
    }

    if (*p1 < *p)
	*p1 = *p;

    return arg;
}

int
main(int argc, char *argv[])
{
    ELEMENT *e;
    GRID *g;
    DOF *ls, *ls_grad;
    char *fn = "quad_test2.mesh";
    INT r, m, marked = 0;
    FLOAT a, res, vol, evals;
    size_t mem;
    double t0;
    BOOLEAN use_marks = TRUE, debug_marks = FALSE;
    const char *types[] = {"-", "0", "+", NULL};
    int type = 0, order0;
    char flt1[48], flt2[48], flt3[48];

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterHandler("-refine", "Refinement levels"
				" (in the form \"first[:last[:step]]\")",
				opts_func, FALSE);
    phgOptionsRegisterHandler("-order", "Quadrature orders"
				" (in the form \"first[:last[:step]]\")",
				opts_func, FALSE);
    phgOptionsRegisterNoArg("-use_marks",
				"Use phgQuadInterfaceMarkGrid",
				&use_marks);
    phgOptionsRegisterNoArg("-debug_marks",
				"Debugging phgQuadInterfaceMarkGrid",
				&debug_marks);
    phgOptionsRegisterKeyword("-type", "Quadrature type", types, &type);

    phgInit(&argc, &argv);

    phgPrintf("Test case: %s, LS_ORDER = %d\n", DESC, LS_ORDER);

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgCheckConformity(g);

    /* create the level set function */
    /* don't use DofInterpolation since there may be curved boundaries */
    ls = phgDofNew(g, LS_ORDER > 0 ? DOF_Pn[LS_ORDER] : DOF_ANALYTIC, 1,
			"level-set", levelset_func);

    /* refine elements on the interface */
    r = refine;
    while (r > 0) {
	int depth = r > 3 ? 3 : r;
	refine_interface(ls, ls_grad, depth);
	r -= depth;
    }

    /* create the gradient of the level set function */
    ls_grad = (LS_ORDER > 0 ? NULL :
		phgDofNew(g, DOF_ANALYTIC, 3, "ls_grad", levelset_grad));

    order0 = order;
    for (; refine <= refine1; refine += refine_step) {
	for (order = order0; order <= order1; order += order_step) {
	    phgPrintf("%"dFMT" elements, order=%d, refine=%d, "
		      "subdiv_level0=%d, subdiv_type=%s.\n",
			g->nelem_global, order, refine,
			phgOptionsGetInt("-qi_subdiv_level0"),
			phgOptionsGetKeyword("-qi_subdiv_type"));

	    if (use_marks)
		marked = phgQuadInterfaceMarkGrid(ls, ls_grad);	/* set marks */

	    /* volume (type-1 == -1,1) or surface (type-1 == 0) integral */
	    t0 = phgGetTime(NULL);
	    res = vol = 0.0;
	    evals = 0.0;		/* total # of function evaluations */
	    m = 0;
#if USE_OMP
#pragma omp parallel for private(e, a) reduction(+:res,vol,m,evals)
#endif	/* USE_OMP */
	    ForAllElementsBegin(g, e)
	    {
		FLOAT *rule, **rules[] = {NULL, NULL, NULL};
		int ret = 0;

		rules[type] = &rule;
		if (use_marks == FALSE)
		    e->mark = 0;
		if (e->mark == 0)
		    m++;
		if (type - 1 != 0) {
		    /* volume integral */
		    if (e->mark == 0) /* element possibly on the interface */
			phgQuadInterface(ls, ls_grad, e, order,
					 rules[0], rules[1], rules[2]);
		    else if ((type - 1) * e->mark > 0)
			/* element entirely in the region */
			phgQuadInterface(NULL, ls, e, order,
					 rules[0], rules[1], rules[2]);
		    else
			/* element entirely not in the region */
			rule = NULL;
		    ret = phgQuadInterfaceRuleApply(one, 1, PROJ_NONE,
						    rule, &a);
		    phgFree(rule);
		    /* for debugging */
		    if (e->mark != 0 && debug_marks) {
			FLOAT b;
			phgQuadInterface(ls, ls_grad, e, order,
					 rules[0], rules[1], rules[2]);
		    	phgQuadInterfaceRuleApply(one, 1, PROJ_NONE,
						  rule, &b);
			phgFree(rule);
			if (Fabs(a - b) > 10.0 * FLOAT_EPSILON) {
			    fprintf(stderr, "Possible wrong mark.\n");
			    fprintf(stderr, "a = %lg, b = %lg\n",
							(double)a, (double)b);
			    phgDumpElement(g, e);
			    //MPI_Finalize();
			    //exit(0);
			}
		    }
		    vol += phgGeomGetVolume(g, e);
		}
		else {
		    /* surface integral */
		    if (e->mark != 0)
			rule = NULL;
		    else
			phgQuadInterface(ls, ls_grad, e, 
				order, rules[0], rules[1], rules[2]);
#if 1
		    /* with normal projection */
		    ret = phgQuadInterfaceRuleApply(normal, 1, PROJ_DOT,
						    rule, &a);
#else
		    /* without projection */
		    Unused(normal);
		    ret = phgQuadInterfaceRuleApply(one, 1, PROJ_NONE,
						    rule, &a);
#endif
		    phgFree(rule);
		}
		res += a;
		evals += ret;
	    }
	    ForAllElementsEnd

	    assert(use_marks == FALSE || m == marked);

#if USE_MPI
	    if (g->nprocs > 1) {
		FLOAT src[4] = {res, vol, (FLOAT)marked, evals}, dst[4];
		MPI_Allreduce(src, dst, 4, PHG_MPI_FLOAT, PHG_SUM, g->comm);
		res = dst[0];
		vol = dst[1];
		marked = (INT)(dst[2] + 0.5);
		evals = dst[3];
	    }
#endif
	    phgMemoryUsage(g, &mem);
	    t0 = phgGetTime(NULL) - t0;
	    a = (type - 1 == 0 ? ARE : VOL);
	    if (type - 1 == 1)
		a = vol - a;
#if FT_PHG == FT___FLOAT128
	    quadmath_snprintf(flt1, sizeof(flt1), "%0.34Qe", res);
	    quadmath_snprintf(flt2, sizeof(flt2), "%0.34Qe", a);
	    quadmath_snprintf(flt3, sizeof(flt3), "%0.34Qe",
					Fabs(1. - res / (a == 0. ? 1. : a)));
#else
	    snprintf(flt1, sizeof(flt1), "%0.16le", (double)res);
	    snprintf(flt2, sizeof(flt2), "%0.16le", (double)a);
	    snprintf(flt3, sizeof(flt3), "%0.16le",
				(double)Fabs(1. - res / (a == 0. ? 1. : a)));
#endif
	    phgPrintf("*** Testing integration %s.\n"
		      "***    Computed:       %s\n"
		      "***    Exact:          %s\n"
		      "***    Relative error: %s\n",
			type - 1 == 0 ? "on the interface" :
				  (type - 1 == -1 ? "in Omega-" : "in Omega+"),
			flt1, flt2, flt3);
	    phgPrintf("*** %"dFMT" elems, evals/marked: %0.0lf/%"dFMT", "
			"%0.3lgMB, %0.2les\n", g->nelem_global, (double)evals,
			marked, mem / (1024.0*1024.0), t0);
	}

	if (refine + refine_step <= refine1)
	    refine_interface(ls, ls_grad, refine_step);
    }

    phgDofFree(&ls);
    if (ls_grad != NULL)
	phgDofFree(&ls_grad);

    /* output the level set function for debugging */
    if (phgOptionsGetNoArg("-qi_dbg_vtk") /*&& LS_ORDER < 0*/) {
	/* perform some refinement for better resolution */
	ls = phgDofNew(g, DOF_P2, 1, "level-set", levelset_func);
	while ((LS_ORDER < 0 || LS_ORDER > 1) && g->nelem_global < 150000)
	    refine_interface(ls, ls_grad, 1);
	if (use_marks) {
	    phgQuadInterfaceMarkGrid(ls, ls_grad);
	    ForAllElements(g, e)
		e->region_mark = e->mark;
	}
	phgExportVTK(g, "level-set.vtk", ls, NULL);
	phgDofFree(&ls);
    }

    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
