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

/* $Id: nonlinear-solver.h,v 1.2 2018/03/22 02:42:19 zlb Exp $ */

#ifndef PHG_NONLINEAR_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif

/* User function prototype for the nonlinear solver.
 * Note:
 *   1. The NS_USER_FUNC can force the solver to stop by returning FALSE.
 *   2. 'ctx' (if !=NULL) carries optional information for the user function.
 */
typedef BOOLEAN (*NS_USER_FUNC)(const FLOAT x[], FLOAT f[], void *ctx);

int phgNonlinearSolver(int n, NS_USER_FUNC f, NS_USER_FUNC fjac, FLOAT x[],
		       int maxits, FLOAT tolf, FLOAT tolx,
		       FLOAT *residual, void *ctx);

#ifdef __cplusplus
}
#endif

#define PHG_NONLINEAR_SOLVER_H
#endif
