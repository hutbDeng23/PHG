/* $Id: cut_test2.c,v 1.9 2017/12/26 07:52:06 zlb Exp $ */

#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define PLANE_CUT	0	/* 0 ==> sphere cut */

static void
cut_plane(FLOAT x, FLOAT y, FLOAT z, FLOAT *val)
/* The equation for the cut plane */
{
#if PLANE_CUT	/*----------------------------------------------------------*/
    *val = 1.0 * (x - 0.5) + 2.0 * (y - 0.5) + 3.0 * (z - 0.5);
# define AREA	0.0
# define VOL	0.0
#else		/*----------------------------------------------------------*/
    *val = Sqrt(Pow(x - 0.5, 2) + Pow(y - 0.5, 2) + Pow(z - 0.5, 2)) - 0.25;
# define AREA	(4.0 * M_PI * 0.25 * 0.25)
# define VOL	(4.0 * M_PI * 0.25 * 0.25 * 0.25 / 3.0)
#endif		/*----------------------------------------------------------*/
}

int
main(int argc, char *argv[])
{
    GRID *g, *g1;
    ELEMENT *e;
    const char *fn = "cube.dat";
    INT mem_max = 400, refines = 15;
    FLOAT eps = 1e-5, tol = 1e-4;
    BOOLEAN vtk = FALSE;
    double time;
    int pass;

    phgOptionsRegisterInt("-mem_max", "Maximum memory", &mem_max);
    phgOptionsRegisterInt("-refines", "Number of refinements", &refines);
    phgOptionsRegisterFloat("-eps", "Precision for cut points", &eps);
    phgOptionsRegisterFloat("-tol", "Tolerance for cross-surface elements"
				" (0 => no check)", &tol);
    phgOptionsRegisterNoArg("-vtk", "Output VTK file", &vtk);
    phgInit(&argc, &argv);

    g = phgNewGrid(-1);
    if (argc > 1)
	fn = argv[1];
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);

    g1 = NULL;
    for (pass = 1; pass <= refines; pass++) {
	VEF_MAP *m;
	double area, vol, *fvals;
	size_t mem_peak;
	INT j;
	BTYPE mask = OWNER | BDRY_USER1;

	phgRefineAllElements(g, 1);
	/*phgRefineRandomElements(g, "50%");*/
	phgBalanceGrid(g, 1.0, 1, NULL, 1.0);
	time = phgGetTime(NULL);
	phgFreeGrid(&g1);
	g1 = phgSurfaceCut2(g, eps, BDRY_USER1, cut_plane);
	time = phgGetTime(NULL) - time;
	assert(g1 != NULL);
	phgCheckConformity(g1);

	if (tol > 0.0) {
	    /* check cross-surface elements */
	    fvals = phgAlloc(g1->nvert * sizeof(*fvals));
	    for (j = 0; j < g1->nvert; j++) {
		COORD *c = g1->verts + j;
		cut_plane((*c)[0], (*c)[1], (*c)[2], &fvals[j]);
	    }
	    ForAllElements(g1, e) {
		int v, flag;
		flag = 0;
		for (v = 0; v < NVert; v++) {
		    FLOAT f = fvals[e->verts[v]];
		    if (f <= tol && f >= -tol) {
			continue;
		    }
		    if (flag == 0) {
			flag = (f <= -tol ? -1 : 1);
			continue;
		    }
		    if ((flag == -1 && f >= tol) || (flag == 1 && f <= -tol)) {
			phgError(1, "Cross-surface element found!\n");
		    }
		}
	    }
	    phgFree(fvals);
	}

	/* compute area of the triangles on the surface */
	m = phgDofSetupVEFMap(g1, NULL, FACE_FLAG);
	area = 0.0;
	for (j = 0; j < g1->nface; j++) {
	    if ((g1->types_face[j] & mask) != mask)
		continue;
	    area += phgGeomGetFaceArea(g1, m->Fmap[j], m->Find[j]);
	}
	phgDofFreeVEFMap(&m);
	vol = 0.0;
	ForAllElements(g1, e) {
	    assert(e->mark == -1.0 || e->mark == 1.0);
	    if (e->mark == -1.0)
		vol += phgGeomGetVolume(g1, e);
	}
#if USE_MPI
	if (g1->nprocs > 1) {
	    double d[4] = {area, vol, 0.0, 0.0};
	    MPI_Reduce(d, d + 2, 2, MPI_DOUBLE, MPI_SUM, 0, g1->comm);
	    area = d[2];
	    vol = d[3];
	}
#endif	/* USE_MPI */
	phgMemoryUsage(g, &mem_peak);
	phgPrintf("Pass %d, %d submeshes, %"dFMT"+%"dFMT" vertices, %"dFMT
		  "+%"dFMT" elements\n", pass, g->nprocs,
		g->nvert_global, g1->nvert_global - g->nvert_global,
		g->nelem_global, g1->nelem_global - g->nelem_global);
	phgPrintf("    Errors: area = %0.4e, vol = %0.4e, mem: %0.4gMB, "
		  "wtime: %0.4g\n", fabs(area - AREA) / AREA,
		fabs(vol - VOL) / VOL, mem_peak / ((double)1024 * 1024), time);
	if (phgVerbosity > 0) {
	    phgInfo(-1, "g1->nvert=%d, g1->nvert_global=%d\n",
			g1->nvert, g1->nvert_global);
	    phgInfo(-1, "g1->nedge=%d, g1->nedge_global=%d\n",
			g1->nedge, g1->nedge_global);
	    phgInfo(-1, "g1->nface=%d, g1->nface_global=%d\n",
			g1->nface, g1->nface_global);
	    phgInfo(-1, "g1->nelem=%d, g1->nelem_global=%d\n",
			g1->nelem, g1->nelem_global);
	    phgInfo(-1, "%d+%d vertices, %d+%d elements.\n", 
			g->nvert_global, g1->nvert_global - g->nvert_global,
			g->nelem_global, g1->nelem_global - g->nelem_global);
	}
#if USE_MPI
	MPI_Bcast(&pass, 1, MPI_INT, 0, g->comm);
#endif	/* USE_MPI */
	if ((double)mem_peak >= mem_max * (double)1024 * (double)1024)
	    break;
    }

    /* save the grid (with a marker function) */
    if (vtk && g1 != NULL) {
	ForAllElements(g1, e)
	    e->region_mark = e->mark;
	phgPrintf("File \"%s\" created.\n",
			phgExportVTK(g1, "cut_test2.vtk", NULL));
    }

    phgFreeGrid(&g1);
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
