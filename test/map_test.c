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

/* $Id: map_test.c,v 1.5 2012/09/13 08:12:31 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

int
main(int argc, char *argv[])
{
  ELEMENT *e;
  GRID *g;
  DOF *u, *v, *u_hp, *v_hp;
  HP_TYPE *hp;
  MAP *map;
  char *fn = "cube.dat";
  char *dof_u = "P2", *dof_v = "P1";
  INT step = 0, pre_refines = 0;

  phgOptionsRegisterFilename("-mesh_file", "Mesh file", &fn);
  phgOptionsRegisterInt("-pre_refines", "Pre-refinements", &pre_refines);
  phgOptionsRegisterString("-dof_u", "DOF type for u", &dof_u);
  phgOptionsRegisterString("-dof_v", "DOF type for v", &dof_v);

  phgInit(&argc, &argv);
  g = phgNewGrid(-1);
  if (!phgImport(g, fn, FALSE))
    phgError(1, "can't read file \"%s\".\n", fn);
  phgRefineAllElements(g, pre_refines);

  phgOptionsSetHandler("-dof_type", dof_u);
  u = phgDofNew(g, DOF_DEFAULT, 1, "u", DofInterpolation);
  phgOptionsSetHandler("-dof_type", dof_v);
  v = phgDofNew(g, DOF_DEFAULT, 1, "v", DofInterpolation);
  phgPrintf("u->type = %s, v->type = %s\n", u->type->name, v->type->name);

  hp = phgHPNew(g, HP_HB);
  u_hp = phgHPDofNew(g, hp, 1, "u_hp", DofInterpolation);
  phgHPFree(&hp);
  hp = phgHPNew(g, HP_HC);
  v_hp = phgHPDofNew(g, hp, 1, "v_hp", DofInterpolation);
  phgHPFree(&hp);

  while (TRUE) {
    if (phgBalanceGrid(g, 1.1, 1, NULL, 0.))
      phgPrintf("Repartition mesh, %d submeshes, load imbalance: %lg\n",
		g->nprocs, (double)g->lif);

    phgPrintf("Testing map with non HP DOFs:\n");
    map = phgMapCreate(u, v, NULL);
    phgPrintf("    nlocal = %d, nglobal = %d\n", map->nlocal, map->nglobal);
    phgMapDestroy(&map);

    phgPrintf("Testing map with HP DOFs:\n");
    ForAllElements(g, e)
	e->hp_order = 1 + GlobalElement(g, e->index) % 4;
    phgHPSetup(u_hp->hp, FALSE);

    ForAllElements(g, e)
	e->hp_order = 1 + (3 - GlobalElement(g, e->index) % 4);
    phgHPSetup(v_hp->hp, FALSE);

    map = phgMapCreate(u_hp, v_hp, NULL);
    phgPrintf("    nlocal = %d, nglobal = %d\n", map->nlocal, map->nglobal);
    phgMapDestroy(&map);

    phgPrintf("Testing map with HP and non HP DOFs:\n");
    map = phgMapCreate(u, u_hp, v, v_hp, NULL);
    phgPrintf("    nlocal = %d, nglobal = %d\n", map->nlocal, map->nglobal);
    phgMapDestroy(&map);

    if (++step >= 1)
	break;
    phgRefineAllElements(g, 1);
  }

  phgDofFree(&u_hp);
  phgDofFree(&v_hp);
  phgDofFree(&u);
  phgDofFree(&v);

  phgFreeGrid(&g);
  phgFinalize();

  return 0;
}
