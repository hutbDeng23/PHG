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

/* $Id: vtk.h,v 1.3 2021/12/04 13:11:28 zlb Exp $ */

#ifndef PHG_VTK_H

#ifdef __cplusplus
extern "C" {
#endif

void phgExportVTK_(GRID *g);
const char *phgExportVTK(GRID *g, const char *fn, DOF *dof, ...);
const char *phgExportVTKn(GRID *g, const char *fn, int ndof, DOF **dofs);

#ifdef __cplusplus
}
#endif

#define PHG_VTK_H
#endif
