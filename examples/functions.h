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

#ifndef TEST_CASE
# define TEST_CASE 3
#endif

#if   TEST_CASE == 0	/*-------------------------------------------------*/
    static const char *fn = "../albert-files/lshape.dat.bz2";
    static INT pre_refines = 0;
    static BOOLEAN has_analytic_solution = FALSE;
    static const char *analytic = "unavailable";
#   define u_(x, y, z)		0.0
#   define grad_u_(g, x, y, z)	g[0] = g[1] = g[2] = 0.0
#   define f_(x, y, z)		1.0
#   define g_(x, y, z)		1.0
#elif TEST_CASE == 1	/*-------------------------------------------------*/
    static const char *fn = "../test/cube.dat";
    static INT pre_refines = 0;
    static BOOLEAN has_analytic_solution = TRUE;
    static const char *analytic = "u=x+y+z";
#   define u_(x, y, z)		x + y + z
#   define grad_u_(g, x, y, z)	g[0] = g[1] = g[2] = 1.0
#   define f_(x, y, z)		0.0
#   define g_(x, y, z)		u_(x, y, z)
#elif TEST_CASE == 2	/*-------------------------------------------------*/
    static const char *fn = "../test/cube.dat";
    static INT pre_refines = 0;
    static BOOLEAN has_analytic_solution = TRUE;
    static const char *analytic = "u=x^4-2y^4+3z^4";
#   define u_(x, y, z)		x*x*x*x - 2.0*y*y*y*y + 3.*z*z*z*z
#   define grad_u_(g, x, y, z)	{g[0]=4.*x*x*x; g[1]=-8.*y*y*y; g[2]=12.*z*z*z;}
#   define f_(x, y, z)		-12.*(x*x - 2.*y*y + 3.*z*z)
#   define g_(x, y, z)		u_(x, y, z)
#elif TEST_CASE == 3	/*-------------------------------------------------*/
    static const char *fn = "../albert-files/cube5.dat";
    static INT pre_refines = 0;
    static BOOLEAN has_analytic_solution = TRUE;
    static const char *analytic = "u=e^{-10(x^2+y^2+z^2)}";
#   define u_(x, y, z)		Exp(-10.0*(x*x+y*y+z*z))
#   define grad_u_(g, x, y, z)	{g[0] = -20.0*x*Exp(-10.0*(x*x+y*y+z*z)); \
				g[1] = -20.0*y*Exp(-10.0*(x*x+y*y+z*z)); \
				g[2] = -20.0*z*Exp(-10.0*(x*x+y*y+z*z));}
#   define f_(x, y, z)		-(400.0*(x*x+y*y+z*z) - 60.0) * \
				Exp(-10.0*(x*x+y*y+z*z))
#   define g_(x, y, z)		u_(x, y, z)
#elif TEST_CASE == 4	/*-------------------------------------------------*/
    static const char *fn = "../test/cube.dat";
    static INT pre_refines = 0;
    static BOOLEAN has_analytic_solution = TRUE;
    static const char *analytic = "u=e^{x+y+z}";
#   define u_(x, y, z)		Exp(x+y+z)
#   define grad_u_(g, x, y, z)	{g[0] = Exp(x+y+z); g[1] = Exp(x+y+z); \
				g[2] = Exp(x+y+z);}
#   define f_(x, y, z)		-3.0*Exp(x+y+z)
#   define g_(x, y, z)		u_(x, y, z)
#elif TEST_CASE == 5	/*-------------------------------------------------*/
    static const char *fn = "../albert-files/cub1o.dat.bz2";
    static INT pre_refines = 0;
    static BOOLEAN has_analytic_solution = FALSE;
    static const char *analytic = "unavailable";
#   define u_(x, y, z)		0.0
#   define grad_u_(g, x, y, z)	g[0] = g[1] = g[2] = 0.0
#   define f_(x, y, z)		1.0
#   define g_(x, y, z)		1.0
#elif TEST_CASE == 6	/*-------------------------------------------------*/
    static const char *fn = "../test/cube.dat";
    static INT pre_refines = 0;
    static BOOLEAN has_analytic_solution = FALSE;
    static const char *analytic = "unavailable";
#   define u_(x, y, z)		0.0
#   define grad_u_(g, x, y, z)	g[0] = g[1] = g[2] = 0.0
#   define f_(x, y, z)		1.0
#   define g_(x, y, z)		(z <= 0.5) ? 1.0 : -1.0
#elif TEST_CASE == 7	/*-------------------------------------------------*/
    static const char *fn = "../test/cube.dat";
    static INT pre_refines = 0;
    static BOOLEAN has_analytic_solution = TRUE;
    static const char *analytic = "u=\\Sqrt{xyz}";
#   define u_(x, y, z)		Sqrt(x*y*z)
#   define grad_u_(g, x, y, z)	{g[0] = (x == 0.0) ? 0.0 : 0.5*Sqrt(y*z/x); \
				g[1] = (y == 0.0) ? 0.0 : 0.5*Sqrt(x*z/y); \
				g[2] = (z == 0.0) ? 0.0 : 0.5*Sqrt(x*y/z);}
#   define f_(x, y, z)		(x == 0.0 || y == 0.0 || z == 0.0) ? 0.0 : \
				0.25*(Sqrt(y*z/x)/x+Sqrt(x*z/y)/y+Sqrt(x*y/z)/z)
#   define g_(x, y, z)		u_(x, y, z)
#elif TEST_CASE == 8	/*-------------------------------------------------*/
    /* XIN Jiping's test case */
    static const char *fn = "../test/cube.dat";
    static INT pre_refines = 0;
    static BOOLEAN has_analytic_solution = TRUE;
    static const char *analytic = "u=e^x e^y \\cos(\\sqrt{2}z)";
#   define u_(x, y, z)		Exp(x)*Exp(y)*Cos(Sqrt(2)*z)
#   define grad_u_(g, x, y, z)	{g[0]=u_(x,y,z); g[1]=u_(x,y,z); \
				 g[2]=-Sqrt(2)*Exp(x)*Exp(y)*Sin(Sqrt(2)*z);}
#   define f_(x, y, z)		0.0
#   define g_(x, y, z)		u_(x, y, z)
#else
# error invalid TEST_CASE!
#endif
