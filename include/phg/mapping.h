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

/* $Id: mapping.h,v 1.1 2014/04/04 04:30:01 zlb Exp $ */

#ifndef PHG_MAPPING_H

/*
 *
 * Mapping defines the map form element reference coordinates to real
 * world coordinates,
 *    X(\Xi) = \sum_i X_i \phi(\Xi).  
 * A mapping could be expressed as some special kind of DOFs,
 * e.g., Pn, DGn, Qn or MONn.
 *
 * 
 * */

typedef DOF_USER_FUNC MAPPING_USER_FUNC;
typedef DOF_USER_FUNC_LAMBDA MAPPING_USER_FUNC_LAMBDA;
typedef DOF_INIT_FUNC MAPPING_INIT_FUNC;
typedef DOF_BASIS_FUNC MAPPING_BASIS_FUNC;
typedef DOF_BASIS_GRAD MAPPING_BASIS_GRAD;
typedef DOF_MAP_FUNC MAPPING_MAP_FUNC;
typedef DOF_TYPE MAPPING_TYPE;
typedef DOF MAPPING;

extern MAPPING_TYPE *MAPPING_DEFAULT;
#define MappingNoData DofNoData 
#define MappingNoAction DofNoAction 
#define MappingInterpolation DofInterpolation 
#define MappingLambdaFunction DofLambdaFunction 


#ifdef __cplusplus
extern "C" {
#endif

#define phgMappingNew_ phgDofNew_ 
#define phgHPMappingNew phgHPDofNew 

#ifdef __cplusplus
}
#endif

#define PHG_MAPPING_H
#endif
