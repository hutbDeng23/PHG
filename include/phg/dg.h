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

/* $Id: dg.h,v 1.1 2014/04/04 04:30:00 zlb Exp $ */

#ifndef PHG_DG_H

extern DOF_TYPE *DOF_DGn[];

extern DOF_TYPE DOF_DG0_;
#define DOF_DG0	(&DOF_DG0_)

extern DOF_TYPE DOF_DG1_;
#define DOF_DG1 (&DOF_DG1_)

extern DOF_TYPE DOF_DG2_;
#define DOF_DG2 (&DOF_DG2_)

extern DOF_TYPE DOF_DG3_;
#define DOF_DG3 (&DOF_DG3_)

extern DOF_TYPE DOF_DG4_;
#define DOF_DG4 (&DOF_DG4_)

extern DOF_TYPE DOF_DG5_;
#define DOF_DG5 (&DOF_DG5_)

extern DOF_TYPE DOF_DG6_;
#define DOF_DG6 (&DOF_DG6_)

extern DOF_TYPE DOF_DG7_;
#define DOF_DG7 (&DOF_DG7_)

extern DOF_TYPE DOF_DG8_;
#define DOF_DG8 (&DOF_DG8_)

extern DOF_TYPE DOF_DG9_;
#define DOF_DG9 (&DOF_DG9_)

extern DOF_TYPE DOF_DG10_;
#define DOF_DG10 (&DOF_DG10_)

extern DOF_TYPE DOF_DG11_;
#define DOF_DG11 (&DOF_DG11_)

extern DOF_TYPE DOF_DG12_;
#define DOF_DG12 (&DOF_DG12_)

extern DOF_TYPE DOF_DG13_;
#define DOF_DG13 (&DOF_DG13_)

extern DOF_TYPE DOF_DG14_;
#define DOF_DG14 (&DOF_DG14_)

extern DOF_TYPE DOF_DG15_;
#define DOF_DG15 (&DOF_DG15_)

extern DOF_TYPE DOF_DG16_;
#define DOF_DG16 (&DOF_DG16_)

extern HP_INFO hp_info_DG;
#define HP_DG (&hp_info_DG)

#define PHG_DG_H
#endif
