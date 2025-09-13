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

/* $Id: quad-cache_test.c,v 1.3 2022/04/24 01:38:55 zlb Exp $ */

#include "phg.h"
#include <string.h>
#include <math.h>

static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = Exp(x * y * z);
}

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    func_u(x, y, z, value);
    *value *= -(x * x * y * y + x * x * z * z + y * y * z * z);
}

static void
quad_cache_test(DOF *u_h, DOF *f_h)
{
    GRID *g = u_h->g;
    ELEMENT *e;
    QCACHE *qc;
    DOF *grad_f_h = phgDofGradient(f_h, NULL, NULL, "grad_f_h");
    int Q_fb, Q_fnb, Q_tmp;
    
    assert(u_h->dim == 1);

    qc = phgQCNew(QD_DEFAULT, u_h);

    /* Test 1: compute \int_F (\grad \phi_j).n ((\grad f_h) \phi_i) . n */
    Q_fb = phgQCAddFECoefficient(qc, grad_f_h, Q_BAS);	/* \grad f_h * BAS */

    /* Test 2: compute \int_F (\grad \phi_j).n ((\grad f_h).n) \phi_i */
    Q_tmp = phgQCAddFEFunction(qc, grad_f_h);
    Q_tmp = phgQCAddProjection(qc, PROJ_DOT, Q_tmp);
    Q_fnb = phgQCAddFIDCoefficient(qc, Q_tmp, Q_BAS); /* (\grad f_h.n) * BAS */

#if USE_OMP
#pragma omp parallel for private(e)
#endif	/* USE_OMP */
    ForAllElementsBegin(g, e) {
	int N = DofNBas(u_h, e);
	int i, j, face;
	FLOAT val, val1, *rule, nv[Dim];

	/* loop on the faces of the elements */
	for (face = 0; face < NFace; face++) {
	    rule = phgQuadGetRule2D(g, e, face, 2 * u_h->type->order - 2);
	    phgQCSetRule(qc, rule, -1.);
	    phgGeomGetFaceOutNormal(g, e, face, nv);
	    phgQCSetConstantNormal(qc, nv);
	    for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
		    val = phgQCIntegrateFace(
				qc, e->index, face, Q_GRAD, PROJ_DOT, j,
				qc, e->index, face, Q_fb,   PROJ_DOT, i);
		    val1 = phgQCIntegrateFace(
				qc, e->index, face, Q_GRAD, PROJ_DOT, j,
				qc, e->index, face, Q_fnb,  PROJ_NONE, i);
		    printf("i = %d, j = %d, val = %lg, error = %lg\n",
				i, j, (double)val, (double)Fabs(val1 - val));
		    if (Fabs(val1 - val) > 1e-10 * (1.0 + Fabs(val)))
			phgError(1, "test failed.\n");
		}
	    }
	    phgFree(rule);
	}
    } ForAllElementsEnd

    phgDofFree(&grad_f_h);
    phgQCFree(&qc);
}

int
main(int argc, char *argv[])
{
    char *fn = "../test/cube.dat";
    GRID *g;
    DOF *u_h, *f_h;
    INT refine = 0;

    phgOptionsRegisterInt("-refine", "Refinement depth", &refine);
    phgOptionsRegisterFilename("-mesh_file", "Mesh file", &fn);

    phgOptionsPreset("-dof_type DG1 -solver gmres");

    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);

    while (refine > 0) {
	int level = refine > 3 ? 3 : refine;
	phgRefineAllElements(g, level);
	phgBalanceGrid(g, 1.2, 0, NULL, 1.0);
	refine -= level;
    }

    phgPrintf("%"dFMT" elements, %"dFMT" submeshes, load imbalance: %lg\n",
			g->nleaf_global, g->nprocs, (double)g->lif);

    phgSetupHalo(g, HALO_FACE);

    u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
    f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h",  func_f);

    quad_cache_test(u_h, f_h);

    phgDofFree(&u_h);
    phgDofFree(&f_h);

    phgFreeGrid(&g);

    phgFinalize();

    return 0;
}
