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

/* Linear solver test with a MATLAB matrix (e.g., created by phgMatDumpMATLAB)
 * $Id: solver_test_matlab.c,v 1.17 2022/08/28 06:59:34 zlb Exp $ */

#include "phg.h"
#include <string.h>
#include <math.h>

int
main(int argc, char *argv[])
{
    VEC *x, *x0;
    MAT *A;
    SOLVER *solver;
    INT i, n, m, nnz;
    double buf[3], time0;
    char *cmd = NULL;
    const char *magic = "solver_test_matlab: start matrix data";
    FILE *f_in = NULL;
    BOOLEAN use_octave = FALSE;

    phgOptionsRegisterNoArg("-use_octave", "Read input matrix using octave",
				&use_octave);
    phgOptionsPreset("+warn_nonopt");
    phgInit(&argc, &argv);

    phgPrintf("Linear solver test using a MATLAB matrix.\n");

    if (argc != (use_octave ? 3 : 2)) {
	phgPrintf("Usage:\n");
	phgPrintf("   %s [PHG_options] Fname.m%s [PHG_options]\n",
			argv[0], use_octave ? " Vname" : "");
	phgPrintf("where:\n");
	phgPrintf("  'Fname.m' is the MATLAB file defining the matrix.\n");
	if (use_octave)
	    phgPrintf("  'Vname' is the MATLAB variable for the matrix.\n");
	phgFinalize();
	exit(1);
    }

    time0 = phgGetTime(NULL);

    if (phgRank == 0) {
	char *line = NULL;
	size_t len = 0;
	if (use_octave) {
	    /* construct MATLAB commands for piping the matrix */
	    cmd = phgAlloc(strlen(argv[1]) + strlen(argv[2]) + 1024);
	    sprintf(cmd, "echo \""
			    "source %s; [n,m]=size(%s); [r,c,v]=find(%s); "
				    /*"stdout=fopen('tmp.dat', 'w'); "*/
			    "fprintf(stdout, '\\n%s\\n'); "
			    "fwrite(stdout, [n m length(c)], 'real*8'); "
			    "fwrite(stdout, [r'; c'; v'], 'real*8'); "
				    /*"fclose(stdout);"*/
		     "\" | /usr/bin/octave --no-gui",
		    argv[1], argv[2], argv[2], magic);
    
	    phgPrintf("Processing file \"%s\" with octave.\n", argv[1]);
	    if ((f_in = popen(cmd, "r")) == NULL) {
		fprintf(stderr, "Execution of the following command failed.\n");
		fprintf(stderr, "%s\n", cmd);
	read_err:
		phgAbort(1);
		exit(1);
	    }
    
	    /* skip data before 'magic' */
	    while (TRUE) {
		if (getline(&line, &len, f_in) == -1L)
		    goto read_err;
		assert(line != NULL && line[strlen(line) - 1] == '\n');
		line[strlen(line) - 1] = '\0';
		if (!strcmp(line, magic))
		    break;
	    }
    
	    /* matrix size */
	    assert(sizeof(double) == 8);
	    if (fread(buf, sizeof(buf[0]), 3, f_in) != 3) {
		fprintf(stderr, "Error reading matrix size.\n");
		phgAbort(1);
		exit(1);
	    }
	}
	else {	/* use_octave */
	    /* read matrix header */
	    double nnz_d;
	   phgPrintf("Reading file \"%s\".\n", argv[1]);
	    if ((f_in = fopen(argv[1], "r")) == NULL)
		phgError(1, "Cannot open file \"%s\", abort.\n", argv[1]);
	    if (fscanf(f_in, "%% matrix size: %"dFMT"x%"dFMT", nnz: %lf\n",
						&n, &m, &nnz_d) != 3)
		phgError(1, "Invalid file header, abort.\n");
	    nnz = (INT)(nnz_d + 0.5);
	    while (TRUE) {
		/* look for "xxx = spconvert([" or "xxx_fid = fopen('" */
		if (getline(&line, &len, f_in) == -1L)
		    phgError(1, "Cannot locate sparse matrix data.\n");
		if (strstr(line, "= spconvert([") != NULL)
		    break;	/* data follows */
		if (strstr(line, "_fid = fopen('") != NULL) {
		    /* xxx_fid = fopen('A.m.dat', 'r'); */
		    char *p, *q;
		    p = strchr(line, '\'') + 1;
		    assert(p != NULL);
		    q = strchr(p, '\'');
		    assert(q != NULL);
		    *q = '\0';
		    fclose(f_in);
		    if ((f_in = fopen(p, "r")) == NULL)
			phgError(1, "Cannot open file \"%s\".\n", p);
		    phgPrintf("Reading matrix data from \"%s\".\n", p);
		    break;
		}
	    }
	    buf[0] = n;
	    buf[1] = m;
	    buf[2] = nnz;
	}	/* use_octave */
	phgFree(line);
    }

    /* Create matrix */
#if USE_MPI
    MPI_Bcast(buf, 3, MPI_INT, 0, phgComm);
#endif	/* USE_MPI */
    n = (INT)(buf[0] + 0.5);
    m = (INT)(buf[1] + 0.5);
    nnz = (INT)(buf[2] + 0.5);
    phgPrintf("Matrix size: %"dFMT"x%"dFMT", nnz: %"dFMT"\n", n, m, nnz);
    if (n != m) {
	fprintf(stderr, "Non square matrix not allowed, abort.\n");
	phgAbort(1);
	exit(1);
    }
    m = n / phgNProcs + (phgRank < n % phgNProcs ? 1 : 0);	/* local size */
    A = phgMatCreate(phgComm, m, n);

    if (phgRank == 0 && nnz > 0) {
	/* Read matrix data */
	for (i = 0; i < nnz; i++) {
	    if (use_octave) {
		if (fread(buf, sizeof(buf[0]), 3, f_in) != 3) {
		    fprintf(stderr, "Error reading matrix data.\n");
		    phgAbort(1);
		    exit(1);
		}
		assert(buf[0] >= 0.99 && buf[1] >= 0.99);
	    }
	    else {
		INT j, k;
		if (fscanf(f_in, "%"dFMT" %"dFMT" %lf\n", &j, &k, &buf[2]) != 3)
		    break;
		buf[0] = j;
		buf[1] = k;
	    }
	    if (buf[2] == 0.0)
		continue;
	    phgMatAddGlobalEntry(A, (INT)(buf[0] + 0.5) - 1,
				    (INT)(buf[1] + 0.5) - 1,
				    (FLOAT)buf[2]);
	}
	(use_octave ? pclose : fclose)(f_in); f_in = NULL;
	free(cmd); cmd = NULL;
	phgPrintf("Time reading matrix: %lgs.\n", phgGetTime(NULL) - time0);
    }

    phgMatAssemble(A);
#if USE_MPI
    buf[1] = (double)A->nnz_d + (double)A->nnz_o;
    MPI_Reduce(&buf[1], &buf[0], 1, MPI_DOUBLE, MPI_SUM, 0, phgComm);
#else	/* USE_MPI */
    buf[0] = (double)A->nnz_d + (double)A->nnz_o;
#endif	/* USE_MPI */

    solver = phgMat2Solver(SOLVER_DEFAULT, A);
    phgMatDestroy(&A);

    /* setup test solution and RHS */
    x = phgMapCreateVec(solver->rhs->map, 1);
    x0 = phgMapCreateVec(solver->rhs->map, 1);
    for (i = 0; i < m; i++)
	x0->data[i] = 1.0;
    phgMatVec(MAT_OP_N, 1.0, solver->mat, x0, 0.0, &solver->rhs);
    phgPrintf("Solving linear system: ");
    phgSolverVecSolve(solver, FALSE, x);
    phgVecAXPBY(1.0, x0, -1.0, &x);
    phgPrintf("nits = %d, residual = %le, error = %le\n",
		solver->nits, (double)solver->residual,
		(double)phgVecNorm2(x, 0, NULL));

    phgVecDestroy(&x);
    phgVecDestroy(&x0);
    phgSolverDestroy(&solver);
    phgFinalize();

    return 0;
}
