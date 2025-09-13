/* $Id: poisson.c,v 1.286 2007/09/27 07:37:20 zlb Exp $
 *
 * This program verifies the 2D quadrature functions with following formular:
 *
 *	\int_e \div(u)\phi = -\int_e u\cdot\grad(\phi) + \int_f (u\cdot n)\phi
 *
 **/

#include "phg.h"

#include <string.h>
#include <math.h>

static void
func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    FLOAT fx, fy, fz;

    fx = (1. - x * x),
    fy = (1. - y * y),
    fz = (1. - z * z);

    *(values++) = 2. * x * fy * fz;
    *(values++) = 2. * fx * y * fz;
    *(values++) = 2. * fx * fy * z;
}

static void
func_divu(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    FLOAT fx, fy, fz;

    fx = (1. - x * x),
    fy = (1. - y * y),
    fz = (1. - z * z);

    *value = 2. * (fy * fz + fx * fz + fx * fy);
}

int
main(int argc, char *argv[])
{
    ELEMENT *e;
    GRID *g;
    DOF *u, *divu;
    FLOAT a, b;
    int i;

    /*------------------- Command-line options ----------------------*/
    static char *fn = "tetra.dat";

    phgOptionsPreset("-dof_type HB5");
    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    u = phgDofNew(g, DOF_DEFAULT, 3, "u", func_u);
    //divu = phgDofNew(g, u->type->grad_type, 1, "divu", func_divu);
    divu = phgDofDivergence(u, NULL, NULL, NULL);

    ForAllElements(g, e) {
	phgPrintf("Element %d:\n", e->index);
	for (i = 0; i < u->type->nbas; i++) {
	    a = - phgQuadDofDotGradBas(e, u, u, i, -1)
		+ phgQuadFaceDofDotBas(e, 0, u, DOF_PROJ_DOT, u, i, -1)
		+ phgQuadFaceDofDotBas(e, 1, u, DOF_PROJ_DOT, u, i, -1)
		+ phgQuadFaceDofDotBas(e, 2, u, DOF_PROJ_DOT, u, i, -1)
		+ phgQuadFaceDofDotBas(e, 3, u, DOF_PROJ_DOT, u, i, -1);
	    b = phgQuadDofDotBas(e, divu, u, i, -1);
	    phgPrintf("\tBasis %3d:   a=%12.5le, b=%12.5le, error=%12.5le\n",
			i, a, b, a - b);
	}
    }

    phgDofFree(&u);
    phgDofFree(&divu);
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
