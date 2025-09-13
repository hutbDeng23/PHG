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

/* $Id: mat_test.c,v 1.58 2021/06/28 10:25:40 zlb Exp $ */

#define NEED_PETSCKSP_H
#include "phg.h"
#include "phg/petsc-utils.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

static INT loop_count = 100;

#if USE_HYPRE
#include "HYPRE_utilities.h"
#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "_hypre_parcsr_mv.h"

#ifdef I
#undef I
#endif

static void
setup_hypre_mat(MAT *A, HYPRE_IJMatrix *pa)
{
    MAT_ROW *row;
    INT i, j, d, I, n0, n, *cols;

    assert(sizeof(INT) == sizeof(int) && sizeof(FLOAT)==sizeof(double));
    phgMatUnpack(A);

    HYPRE_IJMatrixCreate(A->rmap->comm,
			 A->rmap->partition[A->rmap->rank],
			 A->rmap->partition[A->rmap->rank + 1] - 1,
			 A->cmap->partition[A->cmap->rank],
			 A->cmap->partition[A->cmap->rank + 1] - 1,
			 pa);
    HYPRE_IJMatrixSetObjectType(*pa, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(*pa);
    
    I = A->rmap->partition[A->rmap->rank];
    n0 = A->cmap->partition[A->cmap->rank];
    n = 0;
    cols = NULL;
    for (i = 0, row = A->rows; i < A->rmap->nlocal; i++, row++, I++) {
	if (row->ncols > 0) {
	    if (n < row->ncols) {
		phgFree(cols);
		cols = phgAlloc((n = row->ncols) * sizeof(*cols));
	    }
	    for (j = 0; j < row->ncols; j++) {
		if ((d = row->cols[j]) < A->cmap->nlocal)
		    cols[j] = d + n0;
		else
		    cols[j] = A->O2Gmap[d - A->cmap->nlocal];
	    }
	    HYPRE_IJMatrixSetValues(*pa, 1, &row->ncols, &I, cols, row->data);
	}
    }
    phgFree(cols);
    HYPRE_IJMatrixAssemble(*pa);
}
#endif	/* USE_HYPRE */

static void
build_matrix(MAT *mat, DOF *u_h)
{
    int N = u_h->type->nbas;	/* number of basis functions in an element */
    int i, j;
    GRID *g = u_h->g;
    ELEMENT *e;
    FLOAT A[N][N], buf[N];
    INT I[N];

    assert(u_h->dim == 1);
    ForAllElements(g, e) {
	/* compute \int \grad\phi_j \cdot \grad\phi_i making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgMapE2L(mat->cmap, 0, e, i);
	    for (j = 0; j <= i; j++)
		A[j][i] = A[i][j] =
		    phgQuadGradBasDotGradBas(e, u_h, j, u_h, i, QUAD_DEFAULT);
	}

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, NULL, buf, NULL, DOF_PROJ_NONE)) {
		phgMatAddEntries(mat, 1, I + i, N, I, buf); 
	    }
	    else {	/* interior node */
		phgMatAddEntries(mat, 1, I + i, N, I, A[i]); 
	    }
	}
    }
}

int
main(int argc, char *argv[])
{
    GRID *g;
    DOF *u_h;
    MAT *A, *A0, *B;
    MAP *map;
    INT i;
    size_t nnz, mem, mem_peak;
    VEC *x, *y0, *y1, *y2;
    double t0, t1, dnz, dnz1, mflops, mop;
    char *fn = "../test/cube.dat";
    FLOAT mem_max = 300;
    INT refine = 0;
    BOOLEAN petsc = TRUE, hypre = TRUE, matmat = TRUE; 

    phgOptionsRegisterFilename("-mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterInt("-loop_count", "Loop count", &loop_count);
    phgOptionsRegisterInt("-refine", "Refinement level", &refine);
    phgOptionsRegisterFloat("-mem_max", "Maximum memory", &mem_max);
    phgOptionsRegisterNoArg("-petsc", "Test PETSc if available", &petsc);
    phgOptionsRegisterNoArg("-hypre", "Test Hypre if available", &hypre);
    phgOptionsRegisterNoArg("-matmat", "Test MxM", &matmat);

    phgInit(&argc, &argv);
    g = phgNewGrid(-1);
    if (!phgImport(g, fn, TRUE))
	phgError(1, "can't read file \"%s\".\n", fn);
    while (refine > 0) {
	i = refine > 3 ? 3 : refine;
	phgBalanceGrid(g, 1.2, 1, NULL, 0.);
	phgRefineAllElements(g, i);
	phgPrintf("Refining mesh: %"dFMT" elements, %d submeshes, LIF: %lg\n",
			g->nelem_global, g->nprocs, (double)g->lif);
	refine -= i;
    }
    u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofNoAction);

    while (TRUE) {
	phgPrintf("\n");
	if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, %d submeshes, load imbalance: %lg\n",
			g->nprocs, (double)g->lif);
	map = phgMapCreate(u_h, NULL);
	A = phgMapCreateMat(map, map);
	A->handle_bdry_eqns = TRUE;
	build_matrix(A, u_h);
	phgMatAssemble(A);

	/* Note: A is unsymmetric (A' != A) if boundary entries not removed */
	phgMatRemoveBoundaryEntries(A);

#if 0
	/* test block matrix operation */
	A0 = phgMatCreateBlockMatrix(g->comm, 1, 1, &A, NULL);
#else
	A0 = A;
#endif

	phgPrintf("%d DOF, %d elems, %d submeshes, matrix size: %d, LIF: %lg\n",
			DofDataCountGlobal(u_h), g->nleaf_global,
			g->nprocs, A->rmap->nglobal, (double)g->lif);

	/* test PHG mat-vec multiply */
	x = phgMapCreateVec(A->cmap, 1);
	y1 = phgMapCreateVec(A->rmap, 1);
	phgVecRandomize(x, 123);
	phgMatVec(MAT_OP_N, 1.0, A0, x, 0.0, &y1);

	phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	t0 = phgGetTime(NULL);
	for (i = 0; i < loop_count; i++) {
	    phgMatVec(MAT_OP_N, 1.0, A0, x, 0.0, &y1);
	}
	t1 = phgGetTime(NULL);
	mflops = phgPerfGetMflops(g, NULL, NULL);
	y0 = phgVecCopy(y1, NULL);
	nnz = A->nnz_d + A->nnz_o;
#if USE_MPI
	dnz1 = nnz;
	MPI_Reduce(&dnz1, &dnz, 1, MPI_DOUBLE, MPI_SUM, 0, g->comm);
#else
	dnz = nnz;
#endif
	mop = loop_count * (dnz + dnz - A->rmap->nlocal) * 1e-6;

	phgPrintf("\n");
	t1 -= t0;
	phgPrintf("   PHG:  time %0.4lf, nnz %0.16lg, %0.2lfMF (%0.2lfMF)\n",
			t1, dnz, mop / (t1 == 0 ? 1. : t1), mflops);

	/* test trans(A)*x */
	phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	t0 = phgGetTime(NULL);
	for (i = 0; i < loop_count; i++) {
	    phgMatVec(MAT_OP_T, 1.0, A0, x, 0.0, &y1);
	}
	t1 = phgGetTime(NULL);
	mflops = phgPerfGetMflops(g, NULL, NULL);
	t1 -= t0;
	phgPrintf("  A'*x:  time %0.4lf, nnz %0.16lg, %0.2lfMF (%0.2lfMF), "
		  "err: %le\n", t1, dnz, mop / (t1 == 0 ? 1. : t1), mflops,
		 (double)phgVecNorm2(phgVecAXPBY(-1.0, y0, 1.0, &y1), 0, NULL));

if (matmat) {
	/* time A * trans(A) */
	phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	t0 = phgGetTime(NULL);
	B = phgMatMat(MAT_OP_N, MAT_OP_N, 1.0, A, A, 0.0, NULL);
	t1 = phgGetTime(NULL);
	mflops = phgPerfGetMflops(g, NULL, NULL);
	nnz = B->nnz_d + B->nnz_o;
#if USE_MPI
	dnz1 = nnz;
	MPI_Reduce(&dnz1, &dnz, 1, MPI_DOUBLE, MPI_SUM, 0, g->comm);
#else
	dnz = nnz;
#endif
	/* compare B*x <--> A*A*x */
	y2 = phgMatVec(MAT_OP_N, 1.0, B, x, 0.0, NULL);
	phgMatVec(MAT_OP_N, 1.0, A0, y0, 0.0, &y1);
	phgMatDestroy(&B);
	t1 -= t0;
	phgPrintf("   A*A:  time %0.4lf, nnz %0.16lg, %0.2lfMF, err: %le\n",
		  t1, dnz, mflops,
		 (double)phgVecNorm2(phgVecAXPBY(-1.0, y1, 1.0, &y2), 0, NULL));
}

#if USE_PETSC
	MPI_Allreduce(&A->rmap->nlocal, &i, 1, PHG_MPI_INT, MPI_MIN, g->comm);
	if (i <= 0)
	    phgPrintf("\nSkipping PETSc test since A has local size 0.\n");
	if (i > 0 && petsc == TRUE && phgMaxThreads <= 1) {
	    Mat ma, mb;
	    MatInfo info;
	    Vec va, vb, vc;
	    PetscScalar *vec;

	    ma = phgPetscCreateMatAIJ(A);
	    MatCreateVecs(ma, PETSC_NULL, &va);
	    VecDuplicate(va, &vb);
	    VecGetArray(va, &vec);
	    memcpy(vec, x->data, x->map->nlocal * sizeof(*vec));
	    VecRestoreArray(va, &vec);
	    MatMult(ma, va, vb);
	    phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	    t0 = phgGetTime(NULL);
	    for (i = 0; i < loop_count; i++) {
		MatMult(ma, va, vb);
	    }
	    t1 = phgGetTime(NULL);
	    mflops = phgPerfGetMflops(g, NULL, NULL);
	    VecGetArray(vb, &vec);
	    memcpy(y1->data, vec, x->map->nlocal * sizeof(*vec));
	    VecRestoreArray(vb, &vec);

	    MatGetInfo(ma, MAT_GLOBAL_SUM, &info);
	    /*phgPrintf("    --------------------------------------------"
		      "-------------------------\n");*/
	    phgPrintf("\n");
	    t1 -= t0;
	    dnz = info.nz_used;
	    phgPrintf(" PETSc:  time %0.4lf, nnz %0.16lg, %0.2lfMF (%0.2lfMF), "
		      "err: %le\n", t1, dnz, mop / (t1==0 ? 1.:t1), mflops,
		 (double)phgVecNorm2(phgVecAXPBY(-1.0, y0, 1.0, &y1), 0, NULL));

	    phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	    t0 = phgGetTime(NULL);
	    for (i = 0; i < loop_count; i++) {
		MatMultTranspose(ma, va, vb);
	    }
	    t1 = phgGetTime(NULL);
	    mflops = phgPerfGetMflops(g, NULL, NULL);
	    VecGetArray(vb, &vec);
	    memcpy(y1->data, vec, x->map->nlocal * sizeof(*vec));
	    VecRestoreArray(vb, &vec);
	    t1 -= t0;
	    phgPrintf("  A'*x:  time %0.4lf, nnz %0.16lg, %0.2lfMF (%0.2lfMF), "
		      "err: %le\n", t1, dnz, mop / (t1==0 ? 1.:t1), mflops,
		(double)phgVecNorm2(phgVecAXPBY(-1.0, y0, 1.0, &y1), 0, NULL));

if (matmat) {
	    phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	    t0 = phgGetTime(NULL);
	    MatMatMult(ma, ma, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &mb);
	    t1 = phgGetTime(NULL);
	    mflops = phgPerfGetMflops(g, NULL, NULL);
	    t1 -= t0;
	    MatGetInfo(mb, MAT_GLOBAL_SUM, &info);
	    dnz = info.nz_used;
	    VecDuplicate(va, &vc);
	    /* compare B*x <--> A*A*x */
	    MatMult(ma, vb, vc);
	    MatMult(mb, va, vb);
	    VecGetArray(vb, &vec);
	    memcpy(y1->data, vec, x->map->nlocal * sizeof(*vec));
	    VecRestoreArray(vb, &vec);
	    VecGetArray(vc, &vec);
	    memcpy(y2->data, vec, x->map->nlocal * sizeof(*vec));
	    VecRestoreArray(vc, &vec);
	    phgPrintf("   A*A:  time %0.4lf, nnz %0.16lg, %0.2lfMF, err: %le\n",
		  t1, dnz, mflops,
		 (double)phgVecNorm2(phgVecAXPBY(-1.0, y1, 1.0, &y2), 0, NULL));
	    phgPetscVecDestroy(&vc);
}

	    phgPetscMatDestroy(&mb);
	    phgPetscMatDestroy(&ma);
	    phgPetscVecDestroy(&va);
	    phgPetscVecDestroy(&vb);
	}
#endif	/* USE_PETSC */

#if USE_HYPRE
	if (hypre == TRUE && phgMaxThreads <= 1) {
	    HYPRE_IJMatrix ma;
	    HYPRE_IJVector va, vb, vc;
	    HYPRE_ParCSRMatrix  par_ma;
	    hypre_ParCSRMatrix  *par_mb;
	    HYPRE_ParVector	par_va, par_vb, par_vc;
	    HYPRE_Int offset, *ni, start, end;
	    assert(sizeof(INT)==sizeof(int) && sizeof(FLOAT)==sizeof(double));
	    setup_hypre_mat(A, &ma);
	    ni = phgAlloc(2 * A->rmap->nlocal * sizeof(*ni));
	    offset = A->cmap->partition[A->cmap->rank];
	    for (i = 0; i < A->rmap->nlocal; i++)
		ni[i] = i + offset;
	    HYPRE_IJVectorCreate(g->comm, offset, offset + A->rmap->nlocal - 1,
				 &va);
	    HYPRE_IJVectorCreate(g->comm, offset, offset + A->rmap->nlocal - 1,
				 &vb);
	    HYPRE_IJVectorCreate(g->comm, offset, offset + A->rmap->nlocal - 1,
				 &vc);
	    HYPRE_IJVectorSetObjectType(va, HYPRE_PARCSR);
	    HYPRE_IJVectorSetObjectType(vb, HYPRE_PARCSR);
	    HYPRE_IJVectorSetObjectType(vc, HYPRE_PARCSR);
	    HYPRE_IJVectorSetMaxOffProcElmts(va, 0);
	    HYPRE_IJVectorSetMaxOffProcElmts(vb, 0);
	    HYPRE_IJVectorSetMaxOffProcElmts(vc, 0);
	    HYPRE_IJVectorInitialize(va);
	    HYPRE_IJVectorInitialize(vb);
	    HYPRE_IJVectorInitialize(vc);
	    HYPRE_IJMatrixGetObject(ma, (void **)(void *)&par_ma);
	    HYPRE_IJVectorGetObject(va, (void **)(void *)&par_va);
	    HYPRE_IJVectorGetObject(vb, (void **)(void *)&par_vb);
	    HYPRE_IJVectorGetObject(vc, (void **)(void *)&par_vc);
	    HYPRE_IJVectorSetValues(va, A->cmap->nlocal, ni, (double *)x->data);
	    HYPRE_IJVectorAssemble(va);
	    HYPRE_IJVectorAssemble(vb);
	    HYPRE_IJVectorAssemble(vc);

	    HYPRE_IJMatrixGetRowCounts(ma, A->cmap->nlocal,
					ni, ni + A->rmap->nlocal);
	    for (i = 0, nnz = 0; i < A->rmap->nlocal; i++)
		nnz += ni[A->rmap->nlocal + i];
#if USE_MPI
	    dnz1 = nnz;
	    MPI_Reduce(&dnz1, &dnz, 1, MPI_DOUBLE, MPI_SUM, 0, g->comm);
#else
	    dnz = nnz;
#endif

	    HYPRE_ParCSRMatrixMatvec(1.0, par_ma, par_va, 0.0, par_vb);
	    phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	    t0 = phgGetTime(NULL);
	    for (i = 0; i < loop_count; i++) {
		HYPRE_ParCSRMatrixMatvec(1.0, par_ma, par_va, 0.0, par_vb);
	    }
	    t1 = phgGetTime(NULL);
	    mflops = phgPerfGetMflops(g, NULL, NULL);
	    HYPRE_IJVectorGetValues(vb, A->rmap->nlocal, ni, (double*)y1->data);
	    /*phgPrintf("    --------------------------------------------"
		      "-------------------------\n");*/
	    phgPrintf("\n");
	    t1 -= t0;
	    phgPrintf(" HYPRE:  time %0.4lf, nnz %0.16lg, %0.2lfMF (%0.2lfMF), "
		      "err: %le\n", t1, dnz, mop / (t1==0 ? 1.:t1), mflops,
		(double)phgVecNorm2(phgVecAXPBY(-1.0, y0, 1.0, &y1), 0, NULL));

	    phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	    t0 = phgGetTime(NULL);
	    for (i = 0; i < loop_count; i++) {
		HYPRE_ParCSRMatrixMatvecT(1.0, par_ma, par_va, 0.0, par_vb);
	    }
	    t1 = phgGetTime(NULL);
	    mflops = phgPerfGetMflops(g, NULL, NULL);
	    HYPRE_IJVectorGetValues(vb, A->rmap->nlocal, ni, (double*)y1->data);
	    t1 -= t0;
	    phgPrintf("  A'*x:  time %0.4lf, nnz %0.16lg, %0.2lfMF (%0.2lfMF), "
		      "err: %le\n", t1, dnz, mop / (t1==0 ? 1.:t1), mflops,
		(double)phgVecNorm2(phgVecAXPBY(-1.0, y0, 1.0, &y1), 0, NULL));

if (matmat) {
	    phgPerfGetMflops(g, NULL, NULL);	/* reset flops counter */
	    t0 = phgGetTime(NULL);
	    /* Note: 'HYPRE_ParCSRMatrix' is currently typedef'ed to
	     *	     'hypre_ParCSRMatrix *' */
	    par_mb = hypre_ParMatmul((hypre_ParCSRMatrix *)par_ma,
					(hypre_ParCSRMatrix *)par_ma);
	    t1 = phgGetTime(NULL);
	    mflops = phgPerfGetMflops(g, NULL, NULL);
	    start = hypre_ParCSRMatrixFirstRowIndex(par_mb);
	    end = hypre_ParCSRMatrixLastRowIndex(par_mb) + 1;
	    for (i = start, nnz = 0; i < end; i++) {
		HYPRE_Int ncols;
		hypre_ParCSRMatrixGetRow(par_mb, i, &ncols, NULL, NULL);
		hypre_ParCSRMatrixRestoreRow(par_mb, i, &ncols, NULL, NULL);
		nnz += ncols;
	    }
#if USE_MPI
	    dnz1 = nnz;
	    MPI_Reduce(&dnz1, &dnz, 1, MPI_DOUBLE, MPI_SUM, 0, g->comm);
#else
	    dnz = nnz;
#endif
	    /* compare B*x <--> A*A*x */
	    HYPRE_ParCSRMatrixMatvec(1.0, par_ma, par_vb, 0.0, par_vc);
	    HYPRE_ParCSRMatrixMatvec(1.0, (void *)par_mb, par_va, 0.0, par_vb);
	    HYPRE_IJVectorGetValues(vb, A->rmap->nlocal, ni, (double*)y1->data);
	    HYPRE_IJVectorGetValues(vc, A->rmap->nlocal, ni, (double*)y2->data);
	    hypre_ParCSRMatrixDestroy((par_mb));
	    t1 -= t0;
	    phgPrintf("   A*A:  time %0.4lf, nnz %0.16lg, %0.2lfMF, err: %le\n",
		  t1, dnz, mflops,
		 (double)phgVecNorm2(phgVecAXPBY(-1.0, y1, 1.0, &y2), 0, NULL));
}

	    phgFree(ni);
	    HYPRE_IJMatrixDestroy(ma);
	    HYPRE_IJVectorDestroy(va);
	    HYPRE_IJVectorDestroy(vb);
	    HYPRE_IJVectorDestroy(vc);
	}
#endif	/* USE_HYPRE */

	if (A0 != A)
	    phgMatDestroy(&A0);
#if 0
if (A->rmap->nglobal > 1000) {
    VEC *v = phgMapCreateVec(A->rmap, 3);
    for (i = 0; i < v->map->nlocal; i++) {
	v->data[i + 0 * v->map->nlocal] = 1 * (i + v->map->partition[g->rank]);
	v->data[i + 1 * v->map->nlocal] = 2 * (i + v->map->partition[g->rank]);
	v->data[i + 2 * v->map->nlocal] = 3 * (i + v->map->partition[g->rank]);
    }
    phgMatDumpMATLAB(A, "A", "A.m");
    phgVecDumpMATLAB(v, "v", "v.m");
    phgFinalize();
    exit(0);
}
#endif
	phgMatDestroy(&A);
	phgVecDestroy(&x);
	phgVecDestroy(&y0);
	phgVecDestroy(&y1);
	if (matmat)
	    phgVecDestroy(&y2);
	phgMapDestroy(&map);
	mem = phgMemoryUsage(g, &mem_peak);
	dnz = mem / (1024.0 * 1024.0);
	dnz1 = mem_peak / (1024.0 * 1024.0);
	/*phgPrintf("    --------------------------------------------"
		  "-------------------------\n");*/
	phgPrintf("\n");
	phgPrintf("  Memory: current %0.4lgMB, peak %0.4lgMB\n", dnz, dnz1);
#if 0
{
    static int loop_count = 0;
    if (++loop_count == 4)
	break;
}
#endif
	if (sizeof(INT) == sizeof(int) && (g->nelem_global >= INT_MAX / 2 ||
	    DofDataCountGlobal(u_h) >= INT_MAX / 2))
	    break;
	if (mem_peak > 1024 * (size_t)1024 * mem_max)
	    break;
	phgRefineAllElements(g, 1);
    }
    phgDofFree(&u_h);
    phgFreeGrid(&g);
    phgFinalize();

    return 0;
}
