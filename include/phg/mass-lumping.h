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

/* $Id: mass-lumping.h,v 1.11 2022/09/02 07:33:24 zlb Exp $ */

#ifndef PHG_MASS_LUMPING_H

/*-----------------------------------------------------------------------*/
extern DOF_TYPE *DOF_MLn[];

extern DOF_TYPE DOF_ML1_;
#define DOF_ML1 (&DOF_ML1_)

extern DOF_TYPE DOF_ML2_;
#define DOF_ML2 (&DOF_ML2_)

extern DOF_TYPE DOF_ML2m_;
#define DOF_ML2m (&DOF_ML2m_)

extern DOF_TYPE DOF_ML3_;
#define DOF_ML3 (&DOF_ML3_)

extern DOF_TYPE DOF_ML4_;
#define DOF_ML4 (&DOF_ML4_)

extern DOF_TYPE DOF_ML4a_;
#define DOF_ML4a (&DOF_ML4a_)

extern DOF_TYPE DOF_ML4b_;
#define DOF_ML4b (&DOF_ML4b_)

extern DOF_TYPE DOF_P1B_;
#define DOF_P1B (&DOF_P1B_)

/*-----------------------------------------------------------------------*/
extern DOF_TYPE *DOF_MLon[];

extern DOF_TYPE DOF_ML2o_;
#define DOF_ML2o (&DOF_ML2o_)

extern DOF_TYPE DOF_ML3o_;
#define DOF_ML3o (&DOF_ML3o_)

extern DOF_TYPE DOF_ML4o_;
#define DOF_ML4o (&DOF_ML4o_)

/*-----------------------------------------------------------------------*/
extern DOF_TYPE *DOF_MLDGn[];

extern DOF_TYPE DOF_MLDG1_;
#define DOF_MLDG1 (&DOF_MLDG1_)

/*-----------------------------------------------------------------------
 * TEST_TYPES defines the type of 'TST*' types:
 * 	0: disabled, 1: scaled GLL, 2: LTT, 3: Warp & bend
 * (note: scaled GLL nodes are used in MFEM)
 *
 * WARNING: to recompile the library with a different value of TEST_TYPES,
 * it does not suffice to only recompile mass-lumping.c, use instead the
 * following commands to force recompiling files depending on mass-lumping.h:
 *	touch include/phg/mass-lumping.h
 *	cd src
 *	make -s -j4 USER_CFLAGS=-DTEST_TYPES=1 lib
 *-----------------------------------------------------------------------*/
#ifndef TEST_TYPES
# define TEST_TYPES	2	/* !0 ==> enable DOF_TYPES DOF_TSTn[] */
#endif

#if TEST_TYPES

extern DOF_TYPE *DOF_TSTn[];

extern DOF_TYPE DOF_TST1_;
#define DOF_TST1 (&DOF_TST1_)

extern DOF_TYPE DOF_TST2_;
#define DOF_TST2 (&DOF_TST2_)

extern DOF_TYPE DOF_TST3_;
#define DOF_TST3 (&DOF_TST3_)

extern DOF_TYPE DOF_TST4_;
#define DOF_TST4 (&DOF_TST4_)

extern DOF_TYPE DOF_TST5_;
#define DOF_TST5 (&DOF_TST5_)

extern DOF_TYPE DOF_TST6_;
#define DOF_TST6 (&DOF_TST6_)

extern DOF_TYPE DOF_TST7_;
#define DOF_TST7 (&DOF_TST7_)

extern DOF_TYPE DOF_TST8_;
#define DOF_TST8 (&DOF_TST8_)

extern DOF_TYPE DOF_TST9_;
#define DOF_TST9 (&DOF_TST9_)

extern DOF_TYPE DOF_TST10_;
#define DOF_TST10 (&DOF_TST10_)

extern DOF_TYPE DOF_TST11_;
#define DOF_TST11 (&DOF_TST11_)

extern DOF_TYPE DOF_TST12_;
#define DOF_TST12 (&DOF_TST12_)

extern DOF_TYPE DOF_TST13_;
#define DOF_TST13 (&DOF_TST13_)

extern DOF_TYPE DOF_TST14_;
#define DOF_TST14 (&DOF_TST14_)

extern DOF_TYPE DOF_TST15_;
#define DOF_TST15 (&DOF_TST15_)

extern DOF_TYPE DOF_TST16_;
#define DOF_TST16 (&DOF_TST16_)

extern DOF_TYPE DOF_TEST_;
#define DOF_TEST (&DOF_TEST_)

#endif	/*------------------------------------------------------- TEST_TYPES */

typedef struct MASS_LUMPING_ {
    int order0;			/* vertex order */
    int order1;			/* edge order */
    int order2;			/* face order */
    int order3;			/* internal order */
    int norbit;			/* # orbits */
    int npoint;			/* size of points[] (for debugging only) */
    int nmonobas;		/* size of monobas[] (for debugging only) */
    int nweights0;		/* size of weights0[] (for debugging only) */
    int *orbits;		/* array of orbits */
    FLOAT *weights0;		/* weights0[np_vert+np_edge+np_face+np_elem].
				 *
				 * Note: weights0==NULL or weights0[0]<0.0
				 * means a non mass-lumping element, which can
				 * be used to construct ordinary Lagrangian
				 * elements (such as the TST* elements) */
    char *monobas;		/* optional list of monomial bases,
				 * of size [nbas][Dim + 1]. It is defined and
				 * used for more general cases which don't fit
				 * into the order[0123] form. */
    /* buffers initialized in init() */
    FLOAT *weights;		/* weights[nbas] */
    /* buffers initialized in init_buffers() */
    BYTE *face_pow;		/* list of powers for monomials on a face */
    BYTE *elem_pow;		/* list of powers for interior monomials */
    FLOAT *coeff;		/* coefficients of the FE bases */
    /* caches */
    FLOAT *bas_lambda;		/* lambda for cached monomials */
    FLOAT *bas_values;		/* cached values of monomials */
    FLOAT *grad_lambda;		/* lambda for cached grad(monomials) */
    FLOAT *grad_values;		/* cached values of grad(monomials) */
} MASS_LUMPING;

#define PHG_MASS_LUMPING_H
#endif
