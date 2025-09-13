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

/* $Id: phg_tcl.c,v 1.21 2014/12/12 14:11:30 zlb Exp $ */

#include "phg.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>	/* for PATH_MAX */

#if USE_TCL
# include <tcl.h>
# if USE_TK
#  include <tk.h>
# endif
#endif
/* FIXME: without using vtktcl.c, phg.tcl segfaults on exit under LINUX */
#if USE_VTK
#ifdef VTKTCL_C
/*# include "/usr/lib/vtk/tcl/vtktcl.c" */
# include VTKTCL_C
#endif
#endif

GRID *grid_ = NULL;	/* made global for vtkPhgData.cxx */
static char *FileName = NULL;

/* The following two variables are defined in refine.c */
static int display_which_ = 0;
static double wall_time = 0;	/* time spent in last operation */
static double tarray[3], time0;

#define StartTimer	phgGetTime(tarray); time0 = tarray[2];
#define StopTimer	phgGetTime(tarray); wall_time = tarray[2] - time0;

static void
tcl_set_result(Tcl_Interp *interp, const char *fmt, ...)
{
    va_list ap;

#if defined(TCL_MAJOR_VERSION) && defined(TCL_MINOR_VERSION) && \
    ((TCL_MAJOR_VERSION==8 && TCL_MINOR_VERSION>=5) || TCL_MAJOR_VERSION>8)
    /* Tcl >= 8.5 */
    static char result[PATH_MAX + 128];

    va_start(ap, fmt);
    vsprintf(result, fmt, ap);
    Tcl_SetResult(interp, result, TCL_STATIC);
#else
    /* Tcl < 8.5 */
    va_start(ap, fmt);
    vsprintf(interp->result, fmt, ap);
#endif
    va_end(ap);
}

static int
Init(ClientData clientData, Tcl_Interp *interp, int argc,
		                const char *argv[])
{
    Unused(&clientData);
    Unused(interp);
    Unused(&argc);
    Unused(argv);

    return TCL_OK;
}

static int
ShowStatus(ClientData clientData, Tcl_Interp *interp, int argc,
		                const char *argv[])
{
    size_t mem;

    Unused(&clientData);
    Unused(argv);
    Unused(&argc);

    if (grid_ == NULL) {
	tcl_set_result(interp, "no mesh");
	return TCL_OK;
    }

    phgMemoryUsage(grid_, &mem);
    tcl_set_result(interp,
		" %"dFMT" vertices, %"dFMT" element%s, load imbalance: %0.2lf, "
		"volume: %lg (%0.2lfs, %0.2lfMB)",
		grid_->nvert_global, grid_->nelem_global,
		grid_->nelem_global == 1 ? "" : "s",
		(double)grid_->lif, (double)grid_->volume,
		wall_time, mem / 1024.0 / 1024.0);

    return TCL_OK;
}

static int
GetNProcs(ClientData clientData, Tcl_Interp *interp, int argc,
		                const char *argv[])
{
    tcl_set_result(interp, "%d", phgNProcs);

    return TCL_OK;
}

static int
GetNSubmeshes(ClientData clientData, Tcl_Interp *interp, int argc,
		                const char *argv[])
{
    tcl_set_result(interp, "%d", (grid_ == NULL) ? 0 : grid_->nprocs);

    return TCL_OK;
}

#if USE_MPI
static void
bcast_display_which(GRID *g)
{
    ParallelFunction(bcast_display_which, TRUE);

    if (phgNProcs > 1)
	MPI_Bcast(&display_which_, 1, MPI_INT, 0, phgComm);

    return;
}

static ELEMENT *pe;
static BOOLEAN get_submesh_callback CB_ARGS(e) {*(pe++) = *e; return TRUE;}
static GRID *g_;

static void
get_submesh_(GRID *g)
{
    int buffer[2];

    ParallelFunction(get_submesh_, g != NULL && phgNProcs > 1);

    if (display_which_ == 0)
	return;
    if (phgRank != 0 && phgRank != display_which_)
	return;

    if (phgRank == 0) {
	MPI_Status status;
	g_ = g = phgNewGrid(0);
	MPI_Recv(buffer, 2, MPI_INT,
			display_which_, 1, phgComm, &status); 
	g->nvert = buffer[0];
	g->nroot = g->nleaf = buffer[1];
	g->roots = phgAlloc(g->nroot * sizeof(*(g->roots)));
	MPI_Recv(g->roots, sizeof(*(g->roots)) * g->nroot, MPI_BYTE,
			display_which_, 2, phgComm, &status); 
	g->verts = phgAlloc(g->nvert * sizeof(*(g->verts)));
	MPI_Recv(g->verts, sizeof(*(g->verts)) * g->nvert, MPI_BYTE,
			display_which_, 3, phgComm, &status); 
    }
    else {
	static ELEMENT *roots;
	buffer[0] = g->nvert;
	buffer[1] = g->nleaf;
	MPI_Send(buffer, 2, MPI_INT, 0, 1, phgComm); 
	pe = roots = phgAlloc(g->nleaf * sizeof(*roots));
	phgTraverseElements(g, get_submesh_callback);
	MPI_Send(roots, sizeof(*roots) * g->nleaf, MPI_BYTE,
			0, 2, phgComm); 
	phgFree(roots);
	MPI_Send(g->verts, sizeof(*(g->verts)) * g->nvert, MPI_BYTE,
			0, 3, phgComm); 
    }

    return;
}

GRID *
get_submesh(GRID *g)
{
    get_submesh_(g_ = g);

    return g_;
}
#endif	/* USE_MPI */

static int
SelectDisplay(ClientData clientData, Tcl_Interp *interp, int argc,
		                const char *argv[])
{
#if USE_MPI
    int i = display_which_;

    if (grid_ == NULL || grid_->nprocs <= 1) {
	i = 0;
    }
    else {
	if (argc > 1) i = atoi(argv[1]);
	if (i < 0 || i >= grid_->nprocs) {
	    /*phgWarning("invalid mesh no: %d\n", i);*/
	    i = display_which_;
	}
    }
    if (i != display_which_) {
	display_which_ = i;
	bcast_display_which(grid_);
#if USE_VTK
	Tcl_Eval(interp, "phgGrid Modified");
#endif
    }
#else	/* USE_MPI */
    Unused(&clientData);
    Unused(interp);
    Unused(argv);
    Unused(&argc);
#endif	/* USE_MPI */

    tcl_set_result(interp, "%d", display_which_);

    return TCL_OK;
}

#if USE_MPI
static void
bcast_refine_which(GRID *g)
{
    Unused(g);
    ParallelFunction(bcast_refine_which, TRUE);

    if (phgNProcs > 1)
	MPI_Bcast(&refine_which_, 1, MPI_INT, 0, phgComm);

    return;
}
#endif	/* USE_MPI */

static int
SelectRefine(ClientData clientData, Tcl_Interp *interp, int argc,
		                const char *argv[])
{
    Unused(&clientData);
    Unused(interp);

#if USE_MPI
    if (argc > 1) {
	refine_which_ = atoi(argv[1]);
	bcast_refine_which(grid_);
    }
#else	/* USE_MPI */
    Unused(argv);
    Unused(&argc);
#endif	/* USE_MPI */
    return TCL_OK;
}

static int
Finalize(ClientData clientData, Tcl_Interp *interp, int argc,
		                const char *argv[])
{
    Unused(&clientData);
    Unused(interp);
    Unused(&argc);
    Unused(argv);

    phgFinalize();
    return TCL_OK;
}

static char *
ConvertUTF8(const char *utf8)
/* converts an UTF-8 string to system encoding */
{
    static char s[PATH_MAX + 1];
    Tcl_DString dstr;

    Tcl_UtfToExternalDString(NULL, utf8, strlen(utf8), &dstr);
    strcpy(s, Tcl_DStringValue(&dstr));
    Tcl_DStringFree(&dstr);
    return s;
}

static int
ConvertFileName(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    Unused(clientData);

    if (argc != 2) {
	phgInfo(0, "%s: error in arguments.\nUsage: %s filename\n", 
		argv[0], argv[0]);
	return TCL_ERROR;
    }

    tcl_set_result(interp, "%s", ConvertUTF8(argv[1]));
    return TCL_OK;
}

static int
Import(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    const char *fname;

    if (argc != 2) {
	tcl_set_result(interp,
		"phg_tcl: error in arguments.\nUsage: %s filename\n", argv[0]);
	return TCL_ERROR;
    }
    phgFreeGrid(&grid_);

    if (FileName != NULL) {
	free(FileName);
	FileName = NULL;
    }

    /* reset mesh no. */
    display_which_ = 0;
    refine_which_ = -1;
#if USE_MPI
    bcast_display_which(grid_);
    bcast_refine_which(grid_);
#endif

    fname = ConvertUTF8(argv[1]);
    grid_ = phgNewGrid(0);
    phgMemoryUsageReset();
    StartTimer
    phgImport(grid_, fname, FALSE);
    StopTimer
    if (grid_->nroot > 0) {
	FileName = strdup(argv[1]);
#if USE_VTK
	Tcl_Eval(interp, "phgGrid Modified");
#endif
    }
    return TCL_OK;
}

static int
DumpALBERT(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    const char *fname;

    if (argc != 2) {
	phgInfo(0, "%s: error in arguments.\nUsage: %s filename\n", 
		argv[0], argv[0]);
	return TCL_ERROR;
    }
    fname = ConvertUTF8(argv[1]);
    phgInfo(1, "Writing ALBERT file \"%s\" ...\n", fname);
    phgMemoryUsageReset();
    StartTimer
    if (phgExportALBERT(grid_, fname)) {
	if (FileName != NULL) free(FileName);
	FileName = strdup(argv[1]);
	return TCL_OK;
    }
    else {
	return TCL_ERROR;
    }
    StopTimer
}

static int
DumpMedit(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    const char *fname;

    if (argc != 2) {
	phgInfo(0, "%s: error in arguments.\nUsage: %s filename\n", 
		argv[0], argv[0]);
	return TCL_ERROR;
    }
    fname = ConvertUTF8(argv[1]);
    phgInfo(1, "Writing Medit mesh file \"%s\" ...\n", fname);
    phgMemoryUsageReset();
    StartTimer
    if (phgExportMedit(grid_, fname)) {
	if (FileName != NULL) free(FileName);
	FileName = strdup(argv[1]);
	return TCL_OK;
    }
    else {
	return TCL_ERROR;
    }
    StopTimer
}

static int
GetFileName(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    Unused(clientData);
    Unused(argc);
    Unused(*argv);

    tcl_set_result(interp, "%s", FileName == NULL ? "" : FileName);
    return TCL_OK;
}

static int
FreeGrid(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    if (grid_ == NULL)
	return TCL_OK;
    phgFreeGrid(&grid_);
#if USE_VTK
    Tcl_Eval(interp, "phgGrid Modified");
#endif
    return TCL_OK;
}

static int
RefineAll(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    int level;
    if (grid_ == NULL)
	return TCL_OK;
    level = atoi(argv[1]);
    phgMemoryUsageReset();
    StartTimer
    phgRefineAllElements(grid_, level);
    StopTimer
#if USE_VTK
    Tcl_Eval(interp, "phgGrid Modified");
#endif
    return TCL_OK;
}

static int
RefineRandom(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    if (grid_ == NULL)
	return TCL_OK;
    if (argc != 2) {
	phgInfo(0, "Wrong # of arguments for phgRefineRandomElements.\n");
	return TCL_OK;
    }
    phgInfo(1, "Refining %s elements.\n", argv[1]);
    phgMemoryUsageReset();
    StartTimer
    phgRefineRandomElements(grid_, argv[1]);
    StopTimer
#if USE_VTK
    Tcl_Eval(interp, "phgGrid Modified");
#endif
    return TCL_OK;
}

static int
RefineSurface(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    if (grid_ == NULL)
	return TCL_OK;
    phgMemoryUsageReset();
    StartTimer
    phgRefineSurface(grid_);
    StopTimer
#if USE_VTK
    Tcl_Eval(interp, "phgGrid Modified");
#endif
    return TCL_OK;
}

static SHORT level;
static BOOLEAN coarsen_callback CB_ARGS(e) {e->mark = level; return TRUE;}

static void
set_coarsen_mark(GRID *g)
{
#if USE_MPI
    ParallelFunction(set_coarsen_mark, TRUE);
    MPI_Bcast(&level, sizeof(level), MPI_BYTE, 0, g->comm);
#endif
    phgTraverseElements(g, coarsen_callback);
}

static int
Coarsen(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    if (grid_ == NULL)
	return TCL_OK;
    if (argc != 2) {
	phgInfo(0, "Wrong # of arguments for phgCoarsenGrid.\n");
	return TCL_OK;
    }
    phgInfo(1, "max coarsen level: %d\n", -(level = atoi(argv[1])));
    if (level >= 0)
	return TCL_OK;
    set_coarsen_mark(grid_);
    phgMemoryUsageReset();
    StartTimer
    phgCoarsenMarkedElements(grid_);
    StopTimer
#if USE_VTK
    Tcl_Eval(interp, "phgGrid Modified");
#endif
    return TCL_OK;
}

static int
PartitionGrid(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    if (grid_ == NULL)
	return TCL_OK;

    if (argc != 1) {
	phgInfo(0, "Wrong # of arguments for phgPartitionGrid.\n");
	return TCL_OK;
    }

    phgMemoryUsageReset();
    if (grid_->nprocs <= 1) {
	StartTimer
	/*phgPartitionGrid(grid_);*/
	phgRedistributeGrid(grid_);
	StopTimer
	if (grid_->nprocs > 1) {
	    tcl_set_result(interp, "1");
# if USE_VTK
	    Tcl_Eval(interp, "phgGrid Modified");
# endif
	}
    }
    else {
	StartTimer
	phgRedistributeGrid(grid_);
	StopTimer
	tcl_set_result(interp, "1");
# if USE_VTK
	Tcl_Eval(interp, "phgGrid Modified");
# endif
    }

    return TCL_OK;
}

static int
CheckConformity(ClientData clientData, Tcl_Interp *interp, int argc,
		const char *argv[])
{
    if (grid_ == NULL)
	return TCL_OK;
    phgCheckConformity(grid_);
    return TCL_OK;
}

int
Tcl_AppInit(Tcl_Interp *interp)
{
    extern int Vtktcl_Init(Tcl_Interp *interp);	/* in vtktcl.c */
    extern int Vtkphgdatatcl_Init(Tcl_Interp *interp); 

#if 1
    /* Hack to set correct system encoding.
       It is not needed with the zh_CN.patch of Tcl */
    const char *s;
    if (!phgMemcmp(s = getenv("LANG"), "zh_CN.GB", 8) ||
	!phgMemcmp(s = getenv("LANG"), "zh_CN.gb", 8) ||
	!strcmp(s, "zh_CN"))
	Tcl_SetSystemEncoding(interp, "euc-cn");
#endif

    if (Tcl_Init(interp) == TCL_ERROR)
	return TCL_ERROR;
#if USE_TK
    if (Tk_Init(interp) == TCL_ERROR)
	return TCL_ERROR;
#endif
#if USE_VTK
#ifdef VTKTCL_C
    if (Vtktcl_Init(interp) == TCL_ERROR)
	return TCL_ERROR;
#endif
    if (Vtkphgdatatcl_Init(interp) == TCL_ERROR)
	return TCL_ERROR;
#endif
    Tcl_CreateCommand(interp, "phgInit", Init, NULL, NULL);
    Tcl_CreateCommand(interp, "phgShowStatus", ShowStatus, NULL, NULL);
    Tcl_CreateCommand(interp, "phgGetNProcs", GetNProcs, NULL, NULL);
    Tcl_CreateCommand(interp, "phgGetNSubmeshes", GetNSubmeshes, NULL, NULL);
    Tcl_CreateCommand(interp, "phgSelectDisplay", SelectDisplay, NULL, NULL);
    Tcl_CreateCommand(interp, "phgSelectRefine", SelectRefine, NULL, NULL);
    Tcl_CreateCommand(interp, "phgFinalize", Finalize, NULL, NULL);
    Tcl_CreateCommand(interp, "phgConvertFileName", ConvertFileName, NULL,NULL);
    Tcl_CreateCommand(interp, "phgImport", Import, NULL, NULL);
    Tcl_CreateCommand(interp, "phgExportALBERT", DumpALBERT, NULL, NULL);
    Tcl_CreateCommand(interp, "phgExportMedit", DumpMedit, NULL, NULL);
    Tcl_CreateCommand(interp, "phgGetFileName", GetFileName, NULL, NULL);
    Tcl_CreateCommand(interp, "phgFreeGrid", FreeGrid, NULL, NULL);
    Tcl_CreateCommand(interp, "phgRefineAllElements", RefineAll, NULL, NULL);
    Tcl_CreateCommand(interp, "phgRefineRandomElements", RefineRandom,
			NULL, NULL);
    Tcl_CreateCommand(interp, "phgRefineSurface", RefineSurface, NULL, NULL);
    Tcl_CreateCommand(interp, "phgCoarsenGrid", Coarsen, NULL, NULL);
    Tcl_CreateCommand(interp, "phgPartitionGrid", PartitionGrid, NULL, NULL);
    Tcl_CreateCommand(interp, "phgCheckConformity", CheckConformity,
			NULL, NULL);

    return TCL_OK;
}

int
main(int argc, char *argv[])
{
#if 0
    if (argc <= 1) {
	static char *argv0[] = {NULL, "phg.tcl", NULL};
	argv0[0] = argv[0];
	argv = argv0;
	argc = 2;
    }
#endif
#if USE_VTK
    /*setenv("TCLLIBPATH", VTK_DIR "/tcl", 1);*/
    putenv("TCLLIBPATH=" VTK_DIR "/tcl");
#endif
    srand(0);
#if USE_MPI
    RegisterParallelFunction(bcast_display_which);
    RegisterParallelFunction(bcast_refine_which);
    RegisterParallelFunction(get_submesh);
    RegisterParallelFunction(set_coarsen_mark);
#endif
    phgMasterSlave = TRUE;
    phgInit(&argc, &argv);
#if 0
    {
	Tcl_Interp *interp;
	int code;
	if (argc < 2)
	    return 1;
	interp = Tcl_CreateInterp();
	Tcl_CreateCommand(interp, "phgImport", phgImport,
			NULL, NULL);
	code = Tcl_EvalFile(interp, argv[1]);
	Tcl_DeleteInterp(interp);
	return code;
    }
#else
# if USE_TK
    Tk_Main(argc, argv, Tcl_AppInit);
# else
    Tcl_Main(argc, argv, Tcl_AppInit);
# endif
    return 0;	/* to make gcc happy */
#endif
}
