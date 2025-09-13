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

/* $Id: specht.h,v 1.4 2017/06/02 06:10:40 zlb Exp $ */

#ifndef PHG_SPECHT_H

extern DOF_TYPE *DOF_SPTn[];

extern DOF_TYPE DOF_SPT_;
#define DOF_SPT (&DOF_SPT_)

/* Reminder macro cheking D_order for DOF_TYPE.InitFunc */
#define CheckD_order(u) \
    if (((u)->type != NULL ? (u)->type : \
		(u)->hp->info->types[(u)->hp->max_order])->D_order > 0) \
	phgError(1, "unimplemented case at %s:%d (DOF '%s').\n", \
		__FILE__, __LINE__, (u)->name);

#define PHG_SPECHT_H
#endif
