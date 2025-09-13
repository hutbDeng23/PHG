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

#include <math.h>
#ifndef TEST_CASE
# define TEST_CASE -1
#endif


/***************/
/* Check utils */
/***************/

#define NbasFace(u) (3 * (u->type->np_vert + u->type->np_edge)	\
		     + u->type->np_face)

#define FUNCV_DEF(func_f, F_)                       \
    static void                                     \
    func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)	\
    {                                               \
        F_(value, x, y, z);                         \
    }                                               \

/*
 *  Problems defination
 *
 *  */
#define DRIVEN_CAVITY 100
#define CYLINDER 101

/* *****************************************************************
 * test case 1:
 * Stokes equations, analytic solution: u:P3, p:P1
 *
 * ****************************************************************/
#if   TEST_CASE == 1
static const char *fn = "../../geo/cube.D6.medit";
static INT pre_refines = 0;
static FLOAT Re = 1.0;
#define u_(u, x, y, z) {			\
    u[0] = pow(x, 0.3e1);			\
    u[1] = -0.3e1 * y * (x * x + z * z);	\
    u[2] = pow(z, 0.3e1);			\
  };
#define p_(p, x, y, z) {			\
    p[0] = 1 - 2 * x;				\
  };
#define gradp_(gradp, x, y, z) {		\
    gradp[0] = -2;				\
    gradp[1] = 0;				\
    gradp[2] = 0;				\
  };
#define gradu_(gradu, x, y, z) {		\
    gradu[0] = 3 * x * x;			\
    gradu[1] = 0;				\
    gradu[2] = 0;				\
    gradu[3] = -6 * y * x;			\
    gradu[4] = -3 * x * x - 3 * z * z;		\
    gradu[5] = -6 * y * z;			\
    gradu[6] = 0;				\
    gradu[7] = 0;				\
    gradu[8] = 3 * z * z;			\
  };
#define f_(f, x, y, z) {			\
    f[0] = - nu * 6 * x - 2;			\
    f[1] = nu * 12 * y;				\
    f[2] = - nu * 6 * z;			\
  };
#define g1_(g, x, y, z) {			\
    g[0] = 3 * x * x - 1 + 2 * x;		\
    g[1] = 0;					\
    g[2] = 0;					\
  };
#define g2_(g, x, y, z) {			\
    g[0] = -6 * y * x;				\
    g[1] = -3 * x * x - 3 * z * z - 1 + 2 * x;	\
    g[2] = -6 * y * z;				\
  };
#define g3_(g, x, y, z) {			\
    g[0] = 0;					\
    g[1] = 0;					\
    g[2] = 3 * z * z - 1 + 2 * x;		\
  };


/* ************************************************************
 * test case 3:
 * Navier-Stokes analytic
 *
 * ************************************************************/
#elif   TEST_CASE == 3
static const char *fn = "../test/cube.dat";
static INT pre_refines = 0;
static FLOAT Re = 1.0;
#define u_(u, x, y, z) {			\
    u[0] = pow(x, 0.3e1);			\
    u[1] = -0.3e1 * y * (x * x + z * z);	\
    u[2] = pow(z, 0.3e1);			\
  };
#define p_(p, x, y, z) {			\
    p[0] = 3 - 2 * x;				\
  };
#define gradp_(gradp, x, y, z) {		\
    gradp[0] = -2;				\
    gradp[1] = 0;				\
    gradp[2] = 0;				\
  };
#define gradu_(gradu, x, y, z) {		\
    gradu[0] = 3 * x * x;			\
    gradu[1] = 0;				\
    gradu[2] = 0;				\
    gradu[3] = -6 * y * x;			\
    gradu[4] = -3 * x * x - 3 * z * z;		\
    gradu[5] = -6 * y * z;			\
    gradu[6] = 0;				\
    gradu[7] = 0;				\
    gradu[8] = 3 * z * z;			\
  };
#define f_(f, x, y, z) {						\
    f[0] = -6 * nu * x - 2 + 3 * pow((double) x, (double) 5);		\
    f[1] = 12 * nu * y - 6 * y * pow((double) x, (double) 4) - 3 * (-3 * x * x - 3 * z * z) * y * (x * x + z * z) - 6 * y * pow((double) z, (double) 4); \
    f[2] = -6 * nu * z + 3 * pow((double) z, (double) 5);		\
  };
#define g1_(g, x, y, z) {			\
    g[0] = nu * (3 * x * x) - (3 - 2 * x);	\
    g[1] = 0;					\
    g[2] = 0;					\
  };
#define g2_(g, x, y, z) {				\
    g[0] = nu * (-6 * y * x);				\
    g[1] = nu * (-3 * x * x - 3 * z * z) - (3 - 2 * x);	\
    g[2] = nu * (-6 * y * z);				\
  };
#define g3_(g, x, y, z) {			\
    g[0] = 0;					\
    g[1] = 0;					\
    g[2] = nu * (3 * z * z) - (3 - 2 * x);	\
  };

/* ************************************************************
 * test case 4:
 * Navier-Stokes analytic
 *
 * ************************************************************/
#elif   TEST_CASE == 4
static const char *fn = "../test/cube.dat";
static INT pre_refines = 0;
static FLOAT Re = 1.0;
static FLOAT time;
#define t time
#define u_(u, x, y, z) {					\
	u[0] = pow(x, 0.3e1) * (double) (3 - t);		\
	u[1] = -0.3e1 * y * (x * x + z * z) * (double) (3 - t);	\
	u[2] = pow(z, 0.3e1) * (double) (3 - t);		\
    };
#define p_(p, x, y, z) {			\
	p[0] = 2 - t * (2 - t) - 3 * y * x;	\
    };
#define gradp_(gradp, x, y, z) {		\
	gradp[0] = -3 * y;			\
	gradp[1] = -3 * x;			\
	gradp[2] = 0;				\
    };
#define gradu_(gradu, x, y, z) {			\
	gradu[0] = 3 * x * x * (3 - t);			\
	gradu[1] = 0;					\
	gradu[2] = 0;					\
	gradu[3] = -6 * y * x * (3 - t);		\
	gradu[4] = -3 * (x * x + z * z) * (3 - t);	\
	gradu[5] = -6 * y * z * (3 - t);		\
	gradu[6] = 0;					\
	gradu[7] = 0;					\
	gradu[8] = 3 * z * z * (3 - t);			\
    };
#define f_(f, x, y, z) {						\
	f[0] = -pow(x, 0.3e1) - 0.6e1 * nu * x * (double) (3 - t) - (double) (3 * y) + 0.3e1 * pow(x, 0.5e1) * (double)  pow((double) (3 - t), (double) 2); \
	f[1] = 0.3e1 * (double) y * (x * x + z * z) + 0.12e2 * nu * (double) y * (double) (3 - t) - 0.3e1 * x - 0.6e1 * (double) y * pow(x, 0.4e1) * (double)  pow((double) (3 - t), (double) 2) + 0.9e1 * pow(x * x + z * z, 0.2e1) * (double)  pow((double) (3 - t), (double) 2) * (double) y - 0.6e1 * (double) y * pow(z, 0.4e1) * (double)  pow((double) (3 - t), (double) 2); \
	f[2] = -pow(z, 0.3e1) - 0.6e1 * nu * z * (double) (3 - t) + 0.3e1 * pow(z, 0.5e1) * (double)  pow((double) (3 - t), (double) 2); \
    };
#define conv_(conv, x, y, z) {						\
	conv[0] = 0.3e1 * pow(x, 0.5e1) * (double)  pow((double) (3 - t), (double) 2); \
	conv[1] = -0.6e1 * y * pow(x, 0.4e1) * (double)  pow((double) (3 - t), (double) 2) + 0.9e1 * pow(x * x + z * z, 0.2e1) * (double)  pow((double) (3 - t), (double) 2) * y - 0.6e1 * y * pow(z, 0.4e1) * (double)  pow((double) (3 - t), (double) 2); \
	conv[2] = 0.3e1 * pow(z, 0.5e1) * (double)  pow((double) (3 - t), (double) 2); \
    };
#define g1_(g1, x, y, z) {						\
	g1[0] = 3 * nu * x * x * (3 - t) - 2 + t * (2 - t) + 3 * y * x;	\
	g1[1] = 0;							\
	g1[2] = 0;							\
    };
#define g2_(g2, x, y, z) {						\
	g2[0] = -6 * nu * y * x * (3 - t);				\
	g2[1] = -3 * nu * (x * x + z * z) * (3 - t) - 2 + t * (2 - t) + 3 * y * x; \
	g2[2] = -6 * nu * y * z * (3 - t);				\
    };
#define g3_(g3, x, y, z) {						\
	g3[0] = 0;							\
	g3[1] = 0;							\
	g3[2] = 3 * nu * z * z * (3 - t) - 2 + t * (2 - t) + 3 * y * x;	\
    };

/* ************************************************************
 *
 * Kovasznay flow
 * 
 * ************************************************************/
#elif   TEST_CASE == KOVASZNAY
static const char *fn = "../examples/kovasznay.mesh";
static INT pre_refines = 0;
static FLOAT Re = 40.0;
#define LAMBDA (1./(2*nu) - sqrt( 1./(4*nu*nu) + 4*M_PI*M_PI))

#define u_(u, x, y, z) {                                                \
        u[0] = 0.1e1 - exp(LAMBDA * x) * cos((double) (2 * M_PI * y));  \
        u[1] = LAMBDA * exp(LAMBDA * x) * sin((double) (2 * M_PI * y)) / (double) M_PI / 0.2e1; \
        u[2] = 0;                                                       \
    };
#define p_(p, x, y, z) {                                                \
        p[0] = 0.1e1 / 0.2e1 - exp((double) (2 * LAMBDA * x)) / 0.2e1;  \
    };
#define gradp_(gradp, x, y, z) {                                        \
        gradp[0] = -(double) LAMBDA * exp((double) (2 * LAMBDA * x));   \
        gradp[1] = 0;                                                   \
        gradp[2] = 0;                                                   \
    };
#define gradu_(gradu, x, y, z) {                                        \
        gradu[0] = -LAMBDA * exp(LAMBDA * x) * cos((double) (2 * M_PI * y)); \
        gradu[1] = 0.2e1 * exp(LAMBDA * x) * sin((double) (2 * M_PI * y)) * (double) M_PI; \
        gradu[2] = 0;                                                   \
        gradu[3] = LAMBDA * LAMBDA * exp(LAMBDA * x) * sin((double) (2 * M_PI * y)) / (double) M_PI / 0.2e1; \
        gradu[4] = LAMBDA * exp(LAMBDA * x) * cos((double) (2 * M_PI * y)); \
        gradu[5] = 0;                                                   \
        gradu[6] = 0;                                                   \
        gradu[7] = 0;                                                   \
        gradu[8] = 0;                                                   \
    };
#define f_(f, x, y, z) {                                                \
        f[0] = 0;                                                       \
        f[1] = 0;                                                       \
        f[2] = 0;                                                       \
    };
#define conv_(conv, x, y, z) {                                          \
        conv[0] = -LAMBDA * exp(LAMBDA * x) * cos((double) (2 * M_PI * y)) * (0.1e1 - exp(LAMBDA * x) * cos((double) (2 * M_PI * y))) + pow(exp(LAMBDA * x), 0.2e1) * pow(sin((double) (2 * M_PI * y)), 0.2e1) * LAMBDA; \
        conv[1] = LAMBDA * LAMBDA * exp(LAMBDA * x) * sin((double) (2 * M_PI * y)) / (double) M_PI * (0.1e1 - exp(LAMBDA * x) * cos((double) (2 * M_PI * y))) / 0.2e1 + LAMBDA * LAMBDA * pow(exp(LAMBDA * x), 0.2e1) * cos((double) (2 * M_PI * y)) * sin((double) (2 * M_PI * y)) / (double) M_PI / 0.2e1; \
        conv[2] = 0;                                                    \
    };
#define g1_(gx, x, y, z) {						\
	gx[0] = -nu * LAMBDA * exp(LAMBDA * x) * cos((double) (2 * M_PI * y)) - 0.1e1 / 0.2e1 + exp(0.2e1 * LAMBDA * x) / 0.2e1; \
	gx[1] = 0.2e1 * nu * exp(LAMBDA * x) * sin((double) (2 * M_PI * y)) * (double) M_PI; \
	gx[2] = 0;							\
    };
#define g2_(gy, x, y, z) {						\
	gy[0] = nu * LAMBDA * LAMBDA * exp(LAMBDA * x) * sin((double) (2 * M_PI * y)) / (double) M_PI / 0.2e1; \
	gy[1] = nu * LAMBDA * exp(LAMBDA * x) * cos((double) (2 * M_PI * y)) - 0.1e1 / 0.2e1 + exp(0.2e1 * LAMBDA * x) / 0.2e1; \
	gy[2] = 0;							\
    };
#define g3_(gz, x, y, z) {						\
	gz[0] = 0;							\
	gz[1] = 0;							\
	gz[2] = -0.1e1 / 0.2e1 + exp((double) (2 * LAMBDA * x)) / 0.2e1; \
    };


/* ************************************************************
 * driven cavity: 
 *   u = (1, 0, 0) on top of cube
 *       (0, 0, 0) otherwise;
 *
 *
 * ************************************************************/
#elif   TEST_CASE == DRIVEN_CAVITY
static const char *fn = "../test/cube.dat";
static INT pre_refines = 0;
static FLOAT Re = 1.0;
#define ENCLOSE_FLOW
#define u_(u, x, y, z) {			\
    u[0] = (z>0.99999999)?1:0;			\
    u[1] = 0;					\
    u[2] = 0;					\
  };
#define p_(p, x, y, z) {			\
    p[0] = 0;					\
  };
#define gradp_(gradp, x, y, z) {		\
    gradp[0] = 0;				\
    gradp[1] = 0;				\
    gradp[2] = 0;				\
  };
#define gradu_(gradu, x, y, z) {		\
    gradu[0] =0;				\
    gradu[1] =0;				\
    gradu[2] =0;				\
    gradu[3] =0;				\
    gradu[4] =0;				\
    gradu[5] =0;				\
    gradu[6] =0;				\
    gradu[7] =0;				\
    gradu[8] =0;				\
  };
#define f_(f, x, y, z) {			\
    f[0] =0;					\
    f[1] =0;					\
    f[2] =0;					\
  };
#define g1_(g, x, y, z) {			\
    g[0] =0;					\
    g[1] =0;					\
    g[2] =0;					\
  };
#define g2_(g, x, y, z) {			\
    g[0] =0;					\
    g[1] =-10;					\
    g[2] =0;					\
  };
#define g3_(g, x, y, z) {			\
    g[0] =0;					\
    g[1] =0;					\
    g[2] =0;					\
  };


/* ************************************************************
 *
 * flow past cylinder
 *
 * 
 * ************************************************************/
#elif   TEST_CASE == CYLINDER
static FLOAT nu = 0.1;
static const char *fn = "../../geo/channal1.mesh";
static FLOAT Re = 1.0;
#define u_(u, x, y, z) {						\
    u[0] = (x<0.00001) ? 1.0:0;						\
    u[1] = 0.0;								\
    u[2] = 0.0;								\
  };
#define p_(p, x, y, z) {			\
    p[0] = 0.;					\
  };
#define gradp_(gradp, x, y, z) {		\
    gradp[0] = 0;				\
    gradp[1] = 0;				\
    gradp[2] = 0;				\
  };
#define gradu_(gradu, x, y, z) {		\
    gradu[0] =0;				\
    gradu[1] =0;				\
    gradu[2] =0;				\
    gradu[3] =0;				\
    gradu[4] =0;				\
    gradu[5] =0;				\
    gradu[6] =0;				\
    gradu[7] =0;				\
    gradu[8] =0;				\
  };
#define f_(f, x, y, z) {			\
    f[0] =0;					\
    f[1] =0;					\
    f[2] =0;					\
  };
#define g1_(g, x, y, z) {			\
    g[0] =0;					\
    g[1] =0;					\
    g[2] =0;					\
  };
#define g2_(g, x, y, z) {			\
    g[0] =0;					\
    g[1] =0;					\
    g[2] =0;					\
  };
#define g3_(g, x, y, z) {			\
    g[0] =0;					\
    g[1] =0;					\
    g[2] =0;					\
  };

/* ************************************************************
 * end of test case
 * ************************************************************/
#else
# error invalid TEST_CASE!
#endif
