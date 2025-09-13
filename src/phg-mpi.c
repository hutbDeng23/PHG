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

/* $Id: phg-mpi.c,v 1.77 2022/08/25 03:00:47 zlb Exp $ */

#include "phg.h"

#if USE_ZOLTAN
#include "phg/partition-zoltan.h"
#endif	/* USE_ZOLTAN */

GRID *_phg_g_sys = NULL;

int phgRank = 0, phgNProcs = 1, phgCmdSerialno = 0;
MPI_Comm phgComm = MPI_COMM_WORLD;

#if USE_MPI

#include <stdlib.h>

/* for getpid() */
#include <sys/types.h>
#include <unistd.h>
#include <limits.h>	/* PATH_MAX */

static int functions_count = 0;
static PFunc_t *functions = NULL;
static int mpi_initialized = FALSE;

static MPI_Errhandler errhandler = MPI_ERRHANDLER_NULL;

MPI_Comm
phgInitSetComm(MPI_Comm comm)
/* sets up the world communicator (must be called before phgInit()).
 * returns old setting. */
{
    MPI_Comm old_comm;

    if (phgInitialized)
	phgError(1, "%s must be called before phgInit().\n", __func__);
    old_comm = phgComm;
    phgComm = comm;

    return old_comm;
}

int
phgRegisterParallelFunction(PFunc_t func, const char *funcname,
				const char *file, int line)
{
    if (phgInitialized) {
	phgInfo(-1, "%s:%d, unregistered parallel function called.\n",
		file, line);
	phgInfo(-1, "Please add the following line:\n");
	phgInfo(-1, "    RegisterParallelFunction(%s);\n", funcname);
	phgInfo(-1, "before phgInit(...).\n");
	phgError(1, "abort.\n");
    }

    phgInfo(3, "registering parallel function %d: %s\n",
	functions_count, funcname);
    functions = phgRealloc_(functions, sizeof(PFunc_t) * (functions_count + 1),
					sizeof(PFunc_t) * functions_count);
    functions[functions_count] = func;

    return functions_count++;
}

MPI_Op		MPI_MSM;
MPI_Datatype	MPI_3DOUBLE;

/* User function for MPI_Op_create which computes min, sum and max */
static void MPIAPI
MSM(void *in, void *inout, int *len, MPI_Datatype *type)
/* computes b := a op b */
{
    double *a = in, *b = inout;
    int i;

    assert(*type == MPI_3DOUBLE);
    for (i = 0; i < *len; i++) {
	b[0] = (b[0] < a[0] ? b[0] : a[0]); 
	b[1] = a[1] + b[1];
	b[2] = (b[2] > a[2] ? b[2] : a[2]); 
	b += 3;
	a += 3;
    }
}

#if FT_PHG == FT___FLOAT128
MPI_Op		MPI_OP_SUM, MPI_OP_MAX;
MPI_Datatype	MPI_FLOAT128;

/* User function for MPI_OP_SUM */
static void MPIAPI
__sum(void *in, void *inout, int *len, MPI_Datatype *type)
/* computes b := a op b */
{
    FLOAT *a = in, *b = inout;
    int i;

    assert(*type == MPI_FLOAT128);
    for (i = 0; i < *len; i++)
	*(b++) += *(a++);
}

/* User function for MPI_OP_MAX */
static void MPIAPI
__max(void *in, void *inout, int *len, MPI_Datatype *type)
/* computes b := a op b */
{
    FLOAT *a = in, *b = inout;
    int i;

    assert(*type == MPI_FLOAT128);
    for (i = 0; i < *len; i++, a++, b++)
	if (*b < *a)
	    *b = *a;
}
#endif	/* FT_PHG == FT___FLOAT128 */

static void
finalize_mpi(GRID *g)
{
    Unused(g);

    ParallelFunction(finalize_mpi, TRUE);

    if (functions != NULL) {
	phgFree(functions);
	functions = NULL;
	functions_count = 0;
    }

    if (phgNProcs <= 0)
	return;

    MPI_Op_free(&MPI_MSM);
    MPI_Type_free(&MPI_3DOUBLE);

#if FT_PHG == FT___FLOAT128
    MPI_Op_free(&MPI_OP_SUM);
    MPI_Op_free(&MPI_OP_MAX);
    MPI_Type_free(&MPI_FLOAT128);
#endif	/* FT_PHG == FT___FLOAT128 */

#if !NO_NEW_COMMUNICATOR
    MPI_Comm_free(&phgComm);
#endif	/* !NO_NEW_COMMUNICATOR */

    if (errhandler != MPI_ERRHANDLER_NULL)
	MPI_Errhandler_free(&errhandler);

    if (!mpi_initialized) {
	MPI_Finalize();
    }
    else {
	/* to prevent slave processes from exiting */
	phgRank = 0;
    }

    if (phgMasterSlave && phgRank > 0)
	exit(0);

    phgNProcs = -1;

    return;
}

/* Variables for command-line options */
static BOOLEAN l2flag = FALSE;
static INT l2ppn = 0;		/* default to 'use processor names' */
static const char *l2ordering_types[] =
    {"none", "round-robin", "block", NULL};
static int l2ordering_type = L2_BLOCK;	/* default to 'block-wise ordering' */
static BOOLEAN mpi_trap_error = FALSE;

static void MPIAPI
errhandler_function(MPI_Comm *comm, int *code, ...)
{
    int len;
    char errmsg[MPI_MAX_ERROR_STRING + 1];

    MPI_Error_string(*code, errmsg, &len);
    phgError(1, "MPI error: %s\n", errmsg);
}

void
phgInitMPI(int *argc, char **argv[])
{
    char cwd[PATH_MAX];
    static BOOLEAN first_call = TRUE;

    /* The first call is just for registering MPI-related options */
    if (first_call) {
	first_call = FALSE;
	phgOptionsRegisterTitle("\nMPI options:", "\n", "mpi");
	phgOptionsRegisterNoArg("-mpi_trap_error",
				"Use custom MPI error handler",
				&mpi_trap_error);
	phgOptionsRegisterNoArg("-mpi_2level",
				"Build intra-node/inter-node process topology",
				&l2flag);
	phgOptionsRegisterInt("-mpi_2level_ppn",
				"Procs per node (used for testing/debugging, "
				"<=0 => use processor names)", &l2ppn);
	phgOptionsRegisterKeyword("-mpi_2level_ordering",
				"Scheme for 2-level process ordering",
				l2ordering_types, &l2ordering_type);
	phgOptionsRegisterNoArg("-mpi_2level_coll",
				"Enable 2-level collective communications",
				&MPI_L2coll_flag);
	return;
    }

    MPI_Initialized(&mpi_initialized);
    /* save cwd */
    if (getcwd(cwd, sizeof(cwd)) == NULL)
	phgError(1, "can't determine cwd, abort.\n");
    if (!mpi_initialized) {
	assert(phgComm == MPI_COMM_WORLD);
	if (MPI_Init(argc, argv))
	    phgError(1, "MPI_Init failed, abort.\n");
    }

    /* make a duplicate of phgComm */
#if !NO_NEW_COMMUNICATOR
    MPI_Comm_dup(phgComm, &phgComm);
#endif	/* !NO_NEW_COMMUNICATOR */

    /* broadcast and restore cwd */
    MPI_Bcast(cwd, sizeof(cwd), MPI_BYTE, 0, phgComm);
    if (chdir(cwd) == -1)
	phgError(1, "can't change to directory \"%s\", abort.\n", cwd);

    MPI_Comm_rank(phgComm, &phgRank);
    MPI_Comm_size(phgComm, &phgNProcs);

    /* Initialize parallel functions */
    RegisterParallelFunction(finalize_mpi);
    RegisterParallelFunction(phgSetVerbosity_);
    RegisterParallelFunction(phgNewGrid_);
    RegisterParallelFunction(phgRedistributeGrid);
    RegisterParallelFunction(phgRefineMarkedElements);
    RegisterParallelFunction(phgCoarsenMarkedElements);
    RegisterParallelFunction(phgRefineAllElements_);
    RegisterParallelFunction(phgRefineRandom_);
    RegisterParallelFunction(phgRefineSurface);
    RegisterParallelFunction(phgCheckConformity);
    RegisterParallelFunction(phgFreeGrid_);
    RegisterParallelFunction(phgDumpGridInfo);
    RegisterParallelFunction(phgDumpGrid);
    RegisterParallelFunction(phgDumpTree);
    RegisterParallelFunction(phgExportVTK_);
    RegisterParallelFunction(phgExportDX_);
    RegisterParallelFunction(phgExportEnsight_);

    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_3DOUBLE);
    MPI_Type_commit(&MPI_3DOUBLE);
    MPI_Op_create(MSM, 1, &MPI_MSM);

#if FT_PHG == FT___FLOAT128
    MPI_Type_contiguous(SIZEOF_PHG_FLOAT, MPI_BYTE, &MPI_FLOAT128);
    MPI_Type_commit(&MPI_FLOAT128);
    MPI_Op_create(__sum, 1, &MPI_OP_SUM);
    MPI_Op_create(__max, 1, &MPI_OP_MAX);
#endif	/* FT_PHG == FT___FLOAT128 */

    phgOptionsParseCmdline(argc, argv);

#if !NO_NEW_COMMUNICATOR
    /* Create attribute for two-level collective communications */
    if (l2flag)
	MPI_L2setup(&phgComm, l2ordering_type, l2ppn);
    MPI_Comm_size(phgComm, &phgNProcs);
    MPI_Comm_rank(phgComm, &phgRank);
#endif	/* !NO_NEW_COMMUNICATOR */

    if (mpi_trap_error) {
	/* Create custom error handler to facilitate MPI error debugging */
#if 1	/* MPI-2 */
	MPI_Comm_create_errhandler(errhandler_function, &errhandler);
	MPI_Comm_set_errhandler(phgComm, errhandler);
#else
	MPI_Errhandler_create(errhandler_function, &errhandler);
	MPI_Errhandler_set(phgComm, errhandler);
#endif
    }

    /* Initialize pre-defined DOF_TYPEs */
    phgDofTypesInit();

    /* Initialize modules after parsing all options.
     * Note: also add the calls below in phgInit() in utils.c */
    phgSolverInitialize(argc, argv);	/* init phgSolver */
#if USE_ZOLTAN
    phgPartitionZoltanInitialize(*argc, *argv);
#endif	/* USE_ZOLTAN */

    phgInitialized = TRUE;
    phgOptionsHelp();

    if (phgRank == 0) {
	fprintf(stderr, "Parallel Hierarchical Grid (PHG " PHG_VERSION
			", build " PHG_BUILD ").\n");
	fflush(stderr);
	if (phgVerbosity >= 2) {
	    phgPrintf("sizeof(INT)        = %d\n", sizeof(INT));
	    phgPrintf("sizeof(FLOAT)      = %d\n", sizeof(FLOAT));
	    phgPrintf("sizeof(ELEMENT)    = %d\n", sizeof(ELEMENT));
	    phgPrintf("sizeof(GRID)       = %d\n", sizeof(GRID));
	    phgPrintf("sizeof(RNEIGHBOUR) = %d\n", sizeof(RNEIGHBOUR));
	}
	if (phgVerbosity >= 1) {
	    phgPrintf("number of process(es): %d\n", phgNProcs);
	    phgPrintf("%s mode.\n", phgMasterSlave ? "Master-slave" : "SPMD");
	}
	fflush(stdout);
    }

    if (phgVerbosity > 0) {
	char name[MPI_MAX_PROCESSOR_NAME];
	int len;
	MPI_Get_processor_name(name, &len);
	name[len] = '\0';
	phgInfo(1, ">>> processor=%s, pid=%d <<<\n", name, getpid());
    }

    phgTrapSignals_(FALSE);
    phgPerfInit();

    if (phgRank > 0 && phgMasterSlave) {
	if (_phg_pause)
	    phgPause(0);
	phgSlave();
    }

    return;
}

void
phgFinalizeMPI(void)
{
    finalize_mpi(_phg_g_sys);
}

void
phgCallSlaves(int cmd)
{
    if (phgMasterSlave) {
	phgInfo(3, "bringing slave to function %d\n", cmd);
	MPI_Bcast(&cmd, 1, MPI_INT, 0, phgComm);
	phgCmdSerialno++;
    }
}

void
phgSlave(void)
{
    int cmd;

    while (TRUE) {
	/* wait for next command */
	phgInfo(3, "phgSlave: waiting for command.\n");
	MPI_Bcast(&cmd, 1, MPI_INT, 0, phgComm);
	if (cmd >= functions_count) 
	    phgError(1, "phgSlave: invalid command %d\n", cmd);
	phgCmdSerialno++;
	phgInfo(3, "entering function %d, g_sys=%p\n", cmd, _phg_g_sys);
	functions[cmd](_phg_g_sys);
	phgInfo(3, "returning from function %d, g_sys=%p\n", cmd, _phg_g_sys);
    }
}

/* This is for old .o files which refer to phgAlltoallv___() */
int
phgAlltoallv___(
    const void *sendbuf, const int *sendcnts, const int *sdispls,
	MPI_Datatype sendtype,
    void *recvbuf, const int *recvcnts, const int *rdispls,
	MPI_Datatype recvtype,
    MPI_Comm comm, const char *file, int line)
{
    Unused(file);
    Unused(line);
    return MPI_Alltoallv(sendbuf, sendcnts, sdispls, sendtype,
			recvbuf, recvcnts, rdispls, recvtype, comm);
}

#endif /* USE_MPI */
