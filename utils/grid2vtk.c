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

/* This program reads a mesh and creates a VTK file for it.
 * $Id: grid2vtk.c,v 1.5 2021/02/24 08:50:24 zlb Exp $ */

#include "phg.h"

#include <stdlib.h>

int
main(int argc, char *argv[])
{
    const char *in, *out;
    int flag;
    GRID *g;

    phgOptionsPreset("+warn_nonopt");
    phgInit(&argc, &argv);

    if (argc != 3) {
	phgPrintf("Usage: %s mesh_file vtk_file\n", argv[0]);
	phgFinalize();
	exit(1);
    }
    in = argv[1];
    out = argv[2];

    g = phgNewGrid(-1);
    if (!phgImport(g, in, TRUE))
	phgError(1, "can't read file \"%s\".\n", in);

    flag = 0;
    if (phgRank == 0) {
	FILE *f;
	phgDumpGridInfo(g);
	if ((f = fopen(out, "rb")) != NULL) {
	    char s[1024];
	    fclose(f);
	    phgPrintf("File \"%s\" exists, overwrite it (y/n)? ", out);
	    fgets(s, sizeof(s), stdin);
	    if (s[0] != 'y' && s[0] != 'Y')
		flag = 1;
	}
    }
#if USE_MPI
    MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif	/* USE_MPI */
    if (!flag)
	phgPrintf("VTK file \"%s\" created.\n", phgExportVTK(g, out, NULL));
    else
	phgPrintf("No VTK file created.\n");

    phgFreeGrid(&g);
    phgFinalize();
    
    return 0;
}
