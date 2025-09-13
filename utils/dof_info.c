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

/* $Id: dof_info.c,v 1.20 2021/02/24 08:50:24 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <string.h>

/* Note: must match Sxxx macros in mass-lumping.c */
static const char *orbits_name[] =
  {"S1", "S2", "S11", "S3", "S21", "S111", "S4", "S31", "S22", "S211", "S1111"};

int
main(int argc, char *argv[])
{
    int i, n0, n1;

    phgOptionsPreset("+warn_nonopt");
    phgInit(&argc, &argv);

    if (argc != 2) {
	phgPrintf("Print information about a DOF_TYPE.\n");
	phgPrintf("Usage: %s [-verbosity 1] DOF_TYPE\n", argv[0]);
	phgPrintf("(\"%s -help generic\" to get a list of valid DOF_TYPEs)\n",
			argv[0]);
	phgFinalize();
	exit(1);
    }

    phgOptionsSetHandler("-dof_type", argv[1]);

    phgPrintf("DOF type \"%s\":\n", DOF_DEFAULT->name);

    phgPrintf("  F. E. space:       %s\n",
		phgDofFESpaceName(DOF_DEFAULT->fe_space));
    phgPrintf("  Order:             %d\n", DOF_DEFAULT->order);
    phgPrintf("  Dimension:         %d\n", DOF_DEFAULT->dim);
    phgPrintf("  Invariant:         %s\n",
		DOF_DEFAULT->invariant ? "yes" : "no");
    phgPrintf("  Continuity:        C(%d)\n", DOF_DEFAULT->continuity);
    phgPrintf("  Gradient type:     \"%s\"\n", DOF_DEFAULT->grad_type->name);

    if (DOF_DEFAULT->points == NULL || phgVerbosity <= 0) {
	phgPrintf("  Nodal points:      %s\n",
			DOF_DEFAULT->points == NULL ? "no" : "yes");
    }
    else {
	const FLOAT *p = DOF_DEFAULT->points;
	const char *s = "  Nodal points:      ", *s1 = "\t\t     ";
	for (i = 0; i < DOF_DEFAULT->np_vert; i++, p += 1, s = s1)
	    phgPrintf("%svert[%d] = %g\n", s, i, (double)p[0]);
	for (i = 0; i < DOF_DEFAULT->np_edge; i++, p += 2, s = s1)
	    phgPrintf("%sedge[%d] = %g %g\n", s, i,
			(double)p[0], (double)p[1]);
	for (i = 0; i < DOF_DEFAULT->np_face; i++, p += 3, s = s1)
	    phgPrintf("%sface[%d] = %g %g %g\n", s, i,
			(double)p[0], (double)p[1], (double)p[2]);
	for (i = 0; i < DOF_DEFAULT->np_elem; i++, p += 4, s = s1)
	    phgPrintf("%selem[%d] = %g %g %g %g\n", s, i,
			(double)p[0], (double)p[1], (double)p[2], (double)p[3]);
    }


    phgPrintf("  Basis functions:   %d "
		"(vertex %d, edge %d, face %d, element %d)\n",
		DOF_DEFAULT->nbas, DOF_DEFAULT->np_vert, DOF_DEFAULT->np_edge,
		DOF_DEFAULT->np_face, DOF_DEFAULT->np_elem);
    phgPrintf("                     ");
    n0 = n1 = 0;
    if (DOF_DEFAULT->np_vert > 0) {
	n1 = n0 + DOF_DEFAULT->np_vert * NVert;
	phgPrintf("%s%d-%d: VERTEX", n0 == 0 ? "" : ", ", n0, n1 - 1);
	n0 = n1;
    }
    if (DOF_DEFAULT->np_edge > 0) {
	n1 = n0 + DOF_DEFAULT->np_edge * NEdge;
	phgPrintf("%s%d-%d: EDGE", n0 == 0 ? "" : ", ", n0, n1 - 1);
	n0 = n1;
    }
    if (DOF_DEFAULT->np_face > 0) {
	n1 = n0 + DOF_DEFAULT->np_face * NFace;
	phgPrintf("%s%d-%d: FACE", n0 == 0 ? "" : ", ", n0, n1 - 1);
	n0 = n1;
    }
    if (DOF_DEFAULT->np_elem > 0) {
	n1 = n0 + DOF_DEFAULT->np_elem;
	phgPrintf("%s%d-%d: VOLUME", n0 == 0 ? "" : ", ", n0, n1 - 1);
	n0 = n1;
    }
    phgPrintf("\n");

    phgPrintf("  Mass lumping:      ");
    if (DOF_DEFAULT->mass_lumping != NULL) {
	const FLOAT *w = DOF_DEFAULT->mass_lumping->weights0;
	const char *m = DOF_DEFAULT->mass_lumping->monobas;
	const char *indent = "                     ";

	phgPrintf("orders: %d %d %d %d\n",
			DOF_DEFAULT->mass_lumping->order0,
			DOF_DEFAULT->mass_lumping->order1,
			DOF_DEFAULT->mass_lumping->order2,
			DOF_DEFAULT->mass_lumping->order3);

	phgPrintf("%sorbits(%d): ", indent, DOF_DEFAULT->mass_lumping->norbit);

	for (i = 0; i < DOF_DEFAULT->mass_lumping->norbit - 1; i++)
	    phgPrintf("%s,", orbits_name[DOF_DEFAULT->mass_lumping->orbits[i]]);
	phgPrintf("%s\n",orbits_name[DOF_DEFAULT->mass_lumping->orbits[i]]);

	phgPrintf("%smonomials: ", indent);
	if (m != NULL) {
	    INT pb[4], pb0[4];
	    for (i = 0; i < DOF_DEFAULT->nbas; i++) {
		pb[0] = -m[i * 4 + 0];
		pb[1] = -m[i * 4 + 1];
		pb[2] = -m[i * 4 + 2];
		pb[3] = -m[i * 4 + 3];
		qsort(pb, 4, sizeof(pb[0]), phgCompINT);
		if (i > 0 && !memcmp(pb, pb0, sizeof(pb)))
		    continue;
		if (i != 0)
		    phgPrintf(",");
		if (pb[0] != 0) phgPrintf("(%d", -pb[0]);
		if (pb[1] != 0) phgPrintf(",%d", -pb[1]);
		if (pb[2] != 0) phgPrintf(",%d", -pb[2]);
		if (pb[3] != 0) phgPrintf(",%d", -pb[3]);
		phgPrintf(")");
		memcpy(pb0, pb, sizeof(pb));
	    }
	    phgPrintf("\n");
	}
	else {
	    phgPrintf("auto generated\n");
	}

	if (w == NULL || phgVerbosity <= 0) {
	    phgPrintf("%sweights: %s\n", indent, w == NULL ? "no" : "yes");
	}
	else {
	    const char *s = "weights: ", *s1 = "\t\t ";
	    for (i = 0; i < DOF_DEFAULT->np_vert; i++, s = s1)
		phgPrintf("%s%svert[%d] = %g\n", indent, s, i, (double)*(w++));
	    for (i = 0; i < DOF_DEFAULT->np_edge; i++, s = s1)
		phgPrintf("%s%sedge[%d] = %g\n", indent, s, i, (double)*(w++));
	    for (i = 0; i < DOF_DEFAULT->np_face; i++, s = s1)
		phgPrintf("%s%sface[%d] = %g\n", indent, s, i, (double)*(w++));
	    for (i = 0; i < DOF_DEFAULT->np_elem; i++, s = s1)
		phgPrintf("%s%selem[%d] = %g\n", indent, s, i, (double)*(w++));
	}
    }
    else {
	phgPrintf("no\n");
    }

    phgFinalize();

    return 0;
}
