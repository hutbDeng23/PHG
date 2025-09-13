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

/* $Id: quad_test.c,v 1.16 2020/10/20 23:16:53 zlb Exp $
 *
 * Numerical quadrature test. */

#include "phg.h"
#include <stdlib.h>
#include <math.h>

#if FT_PHG == FT___FLOAT128 && (!HAVE_LIBQUADMATH || !HAVE_QUADMATH_H)
#undef Pow
static FLOAT
Pow0(FLOAT base, int p)
{
    if (p == 0) {
	return 1.0Q;
    }
    else if (p == 1) {
	return base;
    }
    else if (p & 1) {
	return Pow0(base, p - 1) * base;
    }
    else {
	FLOAT a = Pow0(base, p / 2);
	return a * a; 
    }
}

static FLOAT
Pow(FLOAT base, FLOAT power)
{
    int p = (int)power;
    assert((FLOAT)p == power && p >= 0);
    return Pow0(base, p);
}
#endif

int
main(int argc, char *argv[])
{
    QUAD *quad;
    int i, p, np;
#if FT_PHG == FT___FLOAT128
    int pmax = 150;
#elif FT_PHG == FT_FLOAT
    int pmax = 30;
#else
    int pmax = 75;
#endif
    FLOAT a, b, error;
    const FLOAT *pw, *pp;
    BOOLEAN dump_flag = FALSE;

    phgOptionsRegisterNoArg("-dump_rules", "Dump quadrature rules to stdout",
			    &dump_flag);

    phgInit(&argc, &argv);

    /* 1D quadrature lambda0^p + lambda1^p */
    fprintf(stderr, "\n========================== 1D quadrature rules.\n");
    for (p = 0; p <= pmax; p++) {
	quad = phgQuadGetQuad1D(p);
	np = quad->npoints;
	pp = quad->points;
	pw = quad->weights;
	b = 0.;
	if (dump_flag)
	    printf("dim = %d, order = %d, npoints = %d\n", 1, p, np);
	for (i = 0; i < np; i++, pp += 2) {
	    if (dump_flag)
		printf("%0.16le %0.16le %0.16le\n",
			(double)pp[0], (double)pp[1], (double)*pw);
	    b += (Pow(pp[0],p) + Pow(pp[1],p)) * *(pw++);
	}
	/* analytic: \int lambda^p = 1 / (p + 1) */
	a = 2. / (FLOAT)(p + 1) * 1.;
	fprintf(stderr,
		"    order %2d, points %4d, result = %le, rel. error = %le\n",
		p, np, (double)b, (double)Fabs(error=(a-b)/b));
	if (Fabs(error) > 10000.0 * FLOAT_EPSILON)
	    fprintf(stderr, "*** WARNING: possibly incorrect result\n");
    }

    /* 2D quadrature lambda0^p + lambda1^p + lambda2^p */
    fprintf(stderr, "\n========================== 2D quadrature rules.\n");
    for (p = 0; p <= pmax; p++) {
	quad = phgQuadGetQuad2D(p);
	np = quad->npoints;
	pp = quad->points;
	pw = quad->weights;
	b = 0.;
	if (dump_flag)
	    printf("dim = %d, order = %d, npoints = %d\n", 2, p, np);
	for (i = 0; i < np; i++, pp += 3) {
	    if (dump_flag)
		printf("%0.16le %0.16le %0.16le %0.16le\n", (double)pp[0],
			(double)pp[1], (double)pp[2], (double)*pw);
	    b += (Pow(pp[0],p) + Pow(pp[1],p) + Pow(pp[2],p)) * *(pw++);
	}
	/* analytic: \int lambda^p = 1 / ((p + 1) * (p + 2)) */
	a = 3. / (FLOAT)((p + 1) * (p + 2)) * 2.;
	fprintf(stderr,
		"    order %2d, points %4d, result = %le, rel. error = %le\n",
		p, np, (double)b, (double)Fabs(error=(a - b) / b));
	if (Fabs(error) > 10000.0 * FLOAT_EPSILON)
	    fprintf(stderr, "*** WARNING: possibly incorrect result\n");
    }

    /* 3D quadrature lambda0^p + lambda1^p + lambda2^p + lambda3^p */
    fprintf(stderr, "\n========================== 3D quadrature rules.\n");
    for (p = 0; p <= pmax; p++) {
	quad = phgQuadGetQuad3D(p);
	np = quad->npoints;
	pp = quad->points;
	pw = quad->weights;
	b = 0.;
	if (dump_flag)
	    printf("dim = %d, order = %d, npoints = %d\n", 3, p, np);
	for (i = 0; i < np; i++, pp += 4) {
	    if (dump_flag)
		printf("%0.16le %0.16le %0.16le %0.16le %0.16le\n",
			(double)pp[0], (double)pp[1], (double)pp[2],
			(double)pp[3], (double)*pw);
	    b += (Pow(pp[0],p) + Pow(pp[1],p) + Pow(pp[2],p) + Pow(pp[3],p)) *
		 *(pw++);
	}
	/* analytic: \int lambda^p = 1 / ((p + 1) * (p + 2) * (p + 3)) */
	a = 4. / (FLOAT)((p + 1) * (p + 2) * (p + 3)) * 6.;
	fprintf(stderr,
		"    order %2d, points %4d, result = %le, rel. error = %le\n",
		p, np, (double)b, (double)Fabs(error=(a - b) / b));
	if (Fabs(error) > 10000.0 * FLOAT_EPSILON)
	    fprintf(stderr, "*** WARNING: possibly incorrect result\n");
    }

    phgFinalize();

    return 0;
}
