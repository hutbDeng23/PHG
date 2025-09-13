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

/* $Id: lagrange.h,v 1.1 2014/04/04 04:30:01 zlb Exp $ */

#ifndef PHG_LAGRANGE_H

extern DOF_TYPE *DOF_Pn[];

extern DOF_TYPE DOF_P0_;
#define DOF_P0 (&DOF_P0_)

extern DOF_TYPE DOF_P1_;
#define DOF_P1 (&DOF_P1_)

extern DOF_TYPE DOF_P2_;
#define DOF_P2 (&DOF_P2_)

extern DOF_TYPE DOF_P3_;
#define DOF_P3 (&DOF_P3_)

extern DOF_TYPE DOF_P4_;
#define DOF_P4 (&DOF_P4_)

/* DOF_P5-DOF_P16 are provided by lagrange1.c */

extern DOF_TYPE DOF_P5_;
#define DOF_P5 (&DOF_P5_)

extern DOF_TYPE DOF_P6_;
#define DOF_P6 (&DOF_P6_)

extern DOF_TYPE DOF_P7_;
#define DOF_P7 (&DOF_P7_)

extern DOF_TYPE DOF_P8_;
#define DOF_P8 (&DOF_P8_)

extern DOF_TYPE DOF_P9_;
#define DOF_P9 (&DOF_P9_)

extern DOF_TYPE DOF_P10_;
#define DOF_P10 (&DOF_P10_)

extern DOF_TYPE DOF_P11_;
#define DOF_P11 (&DOF_P11_)

extern DOF_TYPE DOF_P12_;
#define DOF_P12 (&DOF_P12_)

extern DOF_TYPE DOF_P13_;
#define DOF_P13 (&DOF_P13_)

extern DOF_TYPE DOF_P14_;
#define DOF_P14 (&DOF_P14_)

extern DOF_TYPE DOF_P15_;
#define DOF_P15 (&DOF_P15_)

extern DOF_TYPE DOF_P16_;
#define DOF_P16 (&DOF_P16_)

#define PHG_LAGRANGE_H
#endif
