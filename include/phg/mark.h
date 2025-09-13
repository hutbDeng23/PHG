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

/* $Id: mark.h,v 1.1 2014/04/04 04:30:01 zlb Exp $ */

#ifndef PHG_MARK_H

/******************************************************************************
 * DOF *error and DOF *osc must have one value on each element or each face
 * which is the (L^p norm)^p of the error on the element or face
 * (normally p = 2)
 *                                                                 
 * The type of error indicator defined over element 
 * must be DOF_P0                                                 
 * the type of error indicator defined over face                   
 * must be the same as the type of Face Jump                       
 * the type of data oscillation must be DOF_P0 and                 
 * defined over element                                            
 * ***************************************************************************/

/* Strategies */
typedef enum {
     MARK_DEFAULT,                /* use default strategy */
     MARK_NONE,                        /* don't mark any elements */
     MARK_ALL,                        /* mark all elements (uniformly refine) */
     MARK_MAX,                        /* MAX strategy */
     MARK_GERS,                        /* GERS strategy */
     MARK_MNS,                        /* MNS strategy */
     MARK_EQDIST        /* equidistribution strategy */
} STRATEGY;

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************
 *                  
 *                  Marking function
 *
 *
 *   'theta' >= 0.0, 'zeta' >= 0., 'gamma' >= 0.
 *   'depth' >= 1, int
 *   'caorsen_depth' >= 1, int
 *
 *   'error', 'osc' and 'coarsen_error' should contain elementwise or facewise
 *   errors (to an appropriate power), 'error' cannot be NULL except for the
 *   first call (initialization) or the strategy is MARK_ALL or MARK_NONE.
 *   'osc' may be NULL.
 *   'coarsen_error' can't be NULL if coarsen_strategy is not MARK_NONE.
 *
 *   It's not necessary to specify the types of DOFs,
 *        the subroutine can probe them, but the strategy and
 *        the types of DOFs must match
 *   DOF of coarsening error should contain element data.
 *                  
 *****************************************************************/

void phgMarkElements(
  STRATEGY strategy, DOF *error, FLOAT theta, DOF *osc, FLOAT zeta, int depth,
  STRATEGY coarsen_strategy, DOF *coarsen_error, FLOAT gamma, int coarsen_depth,
  FLOAT tol);

/* refine only */
#define phgMarkRefine(strategy, error, theta, osc, zeta, depth, tol) \
    phgMarkElements(strategy, error, theta, osc, zeta, depth, \
                     MARK_NONE, NULL, 0.5, 1, tol)

/* coarsen only */
#define phgMarkCoarsen(strategy, error, gama, depth, tol) \
    phgMarkElements(MARK_NONE, NULL, 0.5, NULL, 0.5, 1, \
             strategy, error, gama, depth, tol)

#ifdef __cplusplus
}
#endif

#define PHG_MARK_H
#endif
