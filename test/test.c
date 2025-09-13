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

/* $Id: test.c,v 1.33 2012/12/04 09:51:13 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>

int
main(int argc, char *argv[])
{
    const char *fn = "cube.dat";
    ELEMENT *e;
    GRID *g;
    DOF *u;
    int i, n;

    /* Test initialize MPI before phgInit */
    MPI_Init(&argc, &argv);

    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (argc > 1)
	fn = argv[1];
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);
    phgRefineAllElements(g, 0);

    u = phgDofNew(g, DOF_DEFAULT, 1, "u", DofNoAction);
    n = u->type->nbas;
    ForAllElements(g, e) {
	FLOAT xyz[3];
	const FLOAT *c;
	for (i = 0; i < n; i++) {
	    c = phgDofGetElementCoordinates(u, e, i);
	    phgDofGetElementBasisInfo_(u, e, i, NULL, NULL, NULL, xyz, NULL);
	    assert(!memcmp(c, xyz, sizeof(xyz)));
	}
    }
    phgDofFree(&u);

    phgFreeGrid(&g);
    phgFinalize();

    MPI_Finalize();

    return 0;
}
