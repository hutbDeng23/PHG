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

/* $Id: matvec.c,v 1.441 2022/09/23 06:18:40 zlb Exp $ */

/* Functions for basic mat/vec operations. */

#define NEED_GetVariableArgs

/* USE_INDP controls whether use *(p++) (1) or p[i] (0) to access arrays
 * in MatVec (may it to 1 on i386, 0 on x86_64) */
#ifndef USE_INDP
# define USE_INDP	0
#endif	/* !defined(USE_INDP) */

/* Note: "petscksp.h" must be included in "phg.h" before "mpi-utils.h" */
#define NEED_PETSCKSP_H

#include "phg.h"
#include "phg/petsc-utils.h"

#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <stdarg.h>
#include <limits.h>	/* INT_MAX */
#include <unistd.h>	/* getpid() */

#define vec_compare(v1, v2) \
	assert((v1)->nvec == (v2)->nvec && phgMapCompare((v1)->map, (v2)->map))

FLOAT
phgVecNorm2(VEC *src, int which, FLOAT *res)
{
    INT i;
    FLOAT v, *p;
#if USE_MPI
    FLOAT v1;
#endif	/* USE_MPI */

    assert(which < src->nvec);

    if (!src->assembled)
	phgVecAssemble(src);

    v = 0.0;
    p = src->data + which * (size_t)src->map->nlocal;
#if USE_OMP
    if (phgMaxThreads == 1) {
#endif	/* USE_OMP */
      for (i = 0; i < src->map->nlocal; i++)
	 v += p[i] * p[i];
#if USE_OMP
    }
    else {
#pragma omp parallel for schedule(static) reduction(+:v)
      for (i = 0; i < src->map->nlocal; i++)
	 v += p[i] * p[i];
    }
#endif	/* USE_OMP */

#if USE_MPI
    MPI_Allreduce(&v, &v1, 1, PHG_MPI_FLOAT, PHG_SUM, src->map->comm);
    v = v1;
#endif	/* USE_MPI */

    v = Sqrt(v);
    if (res != NULL)
	*res = v;

    return v;
}

FLOAT
phgVecNorm1(VEC *src, int which, FLOAT *res)
{
    INT i;
    FLOAT v, *p;
#if USE_MPI
    FLOAT v1;
#endif	/* USE_MPI */

    if (!src->assembled)
	phgVecAssemble(src);

    assert(which < src->nvec);

    v = 0.0;
    p = src->data + which * (size_t)src->map->nlocal;
    for (i = 0; i < src->map->nlocal; i++)
	v += Fabs(p[i]);

#if USE_MPI
    MPI_Allreduce(&v, &v1, 1, PHG_MPI_FLOAT, PHG_SUM, src->map->comm);
    v = v1;
#endif	/* USE_MPI */

    if (res != NULL)
	*res = v;

    return v;
}

FLOAT
phgVecNormInfty(VEC *src, int which, FLOAT *res)
{
    INT i;
    FLOAT v, *p;
    FLOAT cmp;
#if USE_MPI
    FLOAT v1;
#endif	/* USE_MPI */

    assert(which < src->nvec);

    if (!src->assembled)
	phgVecAssemble(src);

    v = 0.0;
    p = src->data + which * (size_t)src->map->nlocal;
    for (i = 0; i < src->map->nlocal; i++) {
	cmp = Fabs(p[i]);
	if (cmp > v)
	    v = cmp;
    }

#if USE_MPI
    MPI_Allreduce(&v, &v1, 1, PHG_MPI_FLOAT, PHG_MAX, src->map->comm);
    v = v1;
#endif	/* USE_MPI */

    if (res != NULL)
	*res = v;

    return v;
}

VEC *
phgVecCopy(VEC *src, VEC **dest_ptr)
{
    VEC *dest = (dest_ptr == NULL ? NULL : *dest_ptr);

    if (!src->assembled)
	phgVecAssemble(src);

    if (dest == NULL) {
	dest = phgMapCreateVec(src->map, src->nvec);
	if (dest_ptr != NULL)
	    *dest_ptr = dest;
	dest->mat = src->mat;	     
    }
    else {
	MagicCheck(VEC, dest);
	assert(phgMapCompare(dest->map, src->map));
	if (dest->nvec != src->nvec) {
	    phgFree(dest->data);
	    dest->data = phgAlloc(src->nvec * src->map->nlocal * sizeof(FLOAT));
	    phgFree(dest->offp_data);
	    dest->offp_data = phgAlloc(src->nvec * (src->map->localsize
					- src->map->nlocal) * sizeof(FLOAT));
	}
	dest->nvec = src->nvec;
    }

    dest->assembled = TRUE;

    memcpy(dest->data, src->data,
		src->nvec * sizeof(*src->data) * src->map->nlocal);

    return dest;
}

#if 0 * USE_MPI
/* Note: this function is not needed for present. It is saved here for
 * possible future extension */
static void
localize_column_indices(MAT *mat)
/* The columns indices of input matrix are global and this functions
 * rebuilds local indices and O2Gmap.
 *
 * This function is not needed at present and is saved here for future use. */
{
    INT i, j, k, l, n0, n, size = 0, O2Gsize = 0, *found;
    MAT_ROW *row;

#warning FIXME: should use cmap->O2Gmap as the initial map
    phgFree(mat->O2Gmap);
    mat->O2Gmap = NULL;
    phgFree(mat->ordering);
    mat->ordering = NULL;

    if (mat->cmap->nprocs <= 1)
	return;

    n0 = mat->cmap->partition[mat->cmap->rank];
    /* build O2Gmap */
    for (i = 0, row = mat->rows; i < mat->rmap->nlocal; i++, row++) {
	n = row->ncols;
	for (j = 0; j < n; j++) {
	    if ((k = row->cols[j]) >= n0 && k < n0 + mat->cmap->nlocal)
		/* in-process entry */
		continue;
	    /* FIXME: use a hash table instead of binary search? */
	    l = phgBinarySearchINT(O2Gsize, mat->O2Gmap, NULL, k);
	    if (l < O2Gsize && mat->O2Gmap[l] == k)
		continue;
	    /* insert a new entry at location l */
	    if (O2Gsize >= size) {
		mat->O2Gmap = phgRealloc_(mat->O2Gmap,
					  (size + 1024) * sizeof(*mat->O2Gmap),
					  size * sizeof(*mat->O2Gmap));
		size += 1024;
	    }
	    if (l < O2Gsize) {
		memmove(mat->O2Gmap + l + 1, mat->O2Gmap + l,
			(O2Gsize - l) * sizeof(*mat->O2Gmap));
	    }
	    mat->O2Gmap[l] = k;
	    O2Gsize++;
	}
    }

    /* convert column indices to local ones */
    for (i = 0, row = mat->rows; i < mat->rmap->nlocal; i++, row++) {
	n = row->ncols;
	for (j = 0; j < n; j++) {
	    if ((k = row->cols[j]) >= n0 && k < n0 + mat->cmap->nlocal) {
		row->cols[j] -= n0;
		continue;
	    }
	    found = bsearch(&k, mat->O2Gmap, O2Gsize, sizeof(k), phgCompINT);
	    assert(found != NULL);
	    row->cols[j] = mat->cmap->nlocal + (INT)(found - mat->O2Gmap);
	}
    }

    if (O2Gsize == 0)
	return;

    mat->ordering = phgAlloc(O2Gsize * sizeof(*mat->ordering));
    for (i = 0; i < O2Gsize; i++)
	mat->ordering[i] = i;

    mat->localsize = mat->cmap->nlocal + O2Gsize;
}
#endif	/* USE_MPI */

static INT *col_comp_base;
#if USE_OMP
# pragma omp threadprivate(col_comp_base)
#endif	/* USE_OMP */

static int
col_comp(const void *p1, const void *p2)
{
    INT i = col_comp_base[*(const INT *)p1] - col_comp_base[*(const INT *)p2];
    return (i < 0 ? -1 : (i == 0 ? 0 : 1));
}

void
phgMatPack(MAT *mat)
/* packs the matrix and sets up communications for mat-vec operations */
{
    INT i, j, n, *pc;
    FLOAT *pd;
    MAT_ROW *row;
    size_t nnz;
#if USE_MPI
    INT k, m, O2Gsize, *pc_o, *work;
    FLOAT *pd_o;
#endif	/* USE_MPI */

    if (!mat->assembled)
	phgMatAssemble(mat);

    assert(mat->type != PHG_DESTROYED);

    if (mat->type != PHG_UNPACKED)
	return;

    phgInfo(1, "packing matrix rows.\n");
    mat->type = PHG_PACKED;

#if USE_MPI
    j = 2 * mat->rmap->nlocal + 1;
#else
    j = mat->rmap->nlocal + 1;
#endif
    nnz = mat->nnz_d + mat->nnz_o;
    mat->packed_cols = phgAlloc(nnz * sizeof(*mat->packed_cols));
    mat->packed_data = phgAlloc(nnz * sizeof(*mat->packed_data));
    mat->packed_ind = phgAlloc(j * sizeof(*mat->packed_ind));

    /* convert indices and pack data */
#if USE_MPI
    if (mat->cmap->nprocs > 1)
	O2Gsize = mat->localsize - mat->cmap->nlocal;
    else
	O2Gsize = 0;

    assert(O2Gsize <= 0 || (mat->O2Gmap != NULL && mat->ordering != NULL));

    /* work[] contains the inverse map of ordering[] */
    work = phgAlloc(O2Gsize * sizeof(*work));
    for (i = 0; i < O2Gsize; i++)
	work[mat->ordering[i]] = i;
#endif	/* USE_MPI */
    pc = mat->packed_cols;
    pd = mat->packed_data;
#if USE_MPI
    pc_o = pc + mat->nnz_d;
    pd_o = pd + mat->nnz_d;
#endif	/* USE_MPI */
    for (i = 0, row = mat->rows; i < mat->rmap->nlocal; i++, row++) {
	n = row->ncols;
	mat->packed_ind[i] = pc - mat->packed_cols;
#if !USE_MPI
	if (n > 0) {
	    memcpy(pc, row->cols, n * sizeof(*row->cols));
	    memcpy(pd, row->data, n * sizeof(*row->data));
	}
	pc += n;
	pd += n;
#else	/* !USE_MPI */
	mat->packed_ind[mat->rmap->nlocal + i] = pc_o - mat->packed_cols;
	m = 0;
	for (j = 0; j < n; j++) {
	    if ((k = row->cols[j]) < mat->cmap->nlocal) {
		/* in-process entry */
		*(pd++) = row->data[j];
		*(pc++) = k;
		m++;
		continue;
	    }
	    /* off-process entry */
	    assert(k < mat->localsize);
	    *(pd_o++) = row->data[j];
	    /* Note: Global index is computed by O2Gmap[ordering[*pc_o]] */
	    *(pc_o++) = work[k - mat->cmap->nlocal];
	}
	n = m;
#endif	/* !USE_MPI */

	/* Sort in-process entries (better performance?) */
	if (FALSE && n > 0) {
	    for (j = 0; j < n; j++)
		row->cols[j] = j;
	    col_comp_base = pc - n;
	    qsort(row->cols, n, sizeof(*row->cols), col_comp);
	    for (j = 0; j < n; j++) {
		row->data[j] = *(pd - n + row->cols[j]);
		row->cols[j] = *(pc - n + row->cols[j]);
	    }
	    memcpy(pd - n, row->data, n * sizeof(*pd));
	    memcpy(pc - n, row->cols, n * sizeof(*pc));
	}

	phgFree(row->cols);
	phgFree(row->data);
    }
    assert(pc == mat->packed_cols + mat->nnz_d);
    assert(pd == mat->packed_data + mat->nnz_d);
    mat->packed_ind[i] = pc - mat->packed_cols;
#if USE_MPI
    assert(pc_o == pc + mat->nnz_o);
    assert(pd_o == pd + mat->nnz_o);
    mat->packed_ind[mat->rmap->nlocal + i] = pc_o - mat->packed_cols;
#endif	/* USE_MPI */

    phgFree(mat->rows);
    mat->rows = NULL;

#if USE_MPI
    /* Note: the code can be replaced by phgMapCreateCommInfo() in map.c */

    if (mat->cmap->nprocs == 1) {
	assert(work == NULL);
	return;
    }

    for (i = 0; i < O2Gsize; i++)
	work[i] = mat->O2Gmap[mat->ordering[i]];

    if (mat->cinfo != NULL)
	phgMapDestroyCommInfo(&mat->cinfo);

    mat->cinfo = phgMapCreateCommInfo(mat->cmap->comm,
				      mat->cmap->nprocs,
				      mat->cmap->rank,
				      mat->cmap->nlocal,
				      mat->localsize,
				      work,	/* sorted O2Gmap */
				      NULL,	/* natural ordering */
				      mat->cmap->partition);

    phgFree(work);
#endif	/* USE_MPI */
}

void
phgMatUnpack(MAT *mat)
/* unpacks matrix */
{
    INT i, n, *pc = NULL;
    MAT_ROW *row;
    FLOAT *pd = NULL;
#if USE_MPI
    INT j, m, *pc_o = NULL;
    FLOAT *pd_o = NULL;
#endif	/* USE_MPI */

    if (mat->type != PHG_PACKED)
	return;

    phgInfo(1, "unpacking matrix rows.\n");
    mat->type = PHG_UNPACKED;

    if (mat->rows == NULL)
	mat->rows = phgAlloc(mat->rmap->localsize * sizeof(*mat->rows));
    if ((n = mat->rmap->localsize - mat->rmap->nlocal) > 0)
	bzero(mat->rows + mat->rmap->nlocal, n * sizeof(*mat->rows));
    for (i = 0, row = mat->rows; i < mat->rmap->nlocal; i++, row++) {
	pc = mat->packed_cols + mat->packed_ind[i];
	pd = mat->packed_data + mat->packed_ind[i];
	n = (INT)(mat->packed_ind[i + 1] - mat->packed_ind[i]);
#if !USE_MPI
	row->alloc = row->ncols = n;
	row->cols = phgAlloc(n * sizeof(*row->cols));
	row->data = phgAlloc(n * sizeof(*row->data));
	if (n > 0) {
	    memcpy(row->data, pd, n * sizeof(*row->data));
	    memcpy(row->cols, pc, n * sizeof(*row->cols));
	    pd += n;
	    pc += n;
	}
#else	/* !USE_MPI */
	j = mat->rmap->nlocal + i;
	pc_o = mat->packed_cols + mat->packed_ind[j];
	pd_o = mat->packed_data + mat->packed_ind[j];
	m = (INT)(mat->packed_ind[j + 1] - mat->packed_ind[j]);
	row->alloc = row->ncols = n + m;
	row->cols = phgAlloc((n + m) * sizeof(*row->cols));
	row->data = phgAlloc((n + m) * sizeof(*row->data));
	if (n > 0) {
	    memcpy(row->data, pd, n * sizeof(*row->data));
	    pd += n;
	    for (j = 0; j < n; j++) {
		assert(*pc < mat->cmap->nlocal);
		row->cols[j] = *(pc++);
	    }
	}
	if (m > 0) {
	    memcpy(row->data + n, pd_o, m * sizeof(*row->data));
	    pd_o += m;
	    for (j = 0; j < m; j++)
		row->cols[n + j] = mat->ordering[*(pc_o++)] + mat->cmap->nlocal;
	}
#endif	/* !USE_MPI */
    }
    assert(pc == mat->packed_cols + mat->nnz_d);
    assert(pd == mat->packed_data + mat->nnz_d);
#if USE_MPI
    assert(pc_o == pc + mat->nnz_o);
    assert(pd_o == pd + mat->nnz_o);
#endif	/* USE_MPI */

    phgFree(mat->packed_cols);
    mat->packed_cols = NULL;
    phgFree(mat->packed_data);
    mat->packed_data = NULL;
    phgFree(mat->packed_ind);
    mat->packed_ind = NULL;

#if USE_MPI
    /* Invalidate mat->cinfo since the ordering of offp-columns is different */
    if (mat->cinfo != NULL)
	phgMapDestroyCommInfo(&mat->cinfo);
#endif	/* USE_MPI */
}

void
phgMatSetupDiagonal(MAT *mat)
{
    INT i, k, n, *pc;
    FLOAT *pd;
    MAT_ROW *row;

    phgMatAssemble(mat);

    if (mat->diag != NULL)
	return;

    assert(mat->type != PHG_DESTROYED);

    mat->diag = phgAlloc(mat->rmap->nlocal * sizeof(*mat->diag));

    if (mat->type == PHG_MATRIX_FREE) {
	/* get the diagonal by multiplying mat with the vector [1 ... 1]^T */
	VEC *x, *y;
	phgFree(mat->diag);
	x = phgMapCreateVec(mat->rmap, 1);
	for (i = 0; i < x->map->nlocal; i++)
	    x->data[i] = 1.0;
	y = phgMatVec(MAT_OP_D, 1.0, mat, x, 0.0, NULL);
	mat->diag = y->data;
	y->data = NULL;
	phgVecDestroy(&x);
	phgVecDestroy(&y);
    }
    else if (mat->type != PHG_PACKED) {
	for (i = 0, row = mat->rows; i < mat->rmap->nlocal; i++, row++) {
	    for (k = 0; k < row->ncols; k++) {
		if (row->cols[k] == i)
		    break;
	    }
	    mat->diag[i] = (k < row->ncols ? row->data[k] : 0.0);
	}
    }
    else {
	for (i = 0; i < mat->rmap->nlocal; i++) {
	    pc = mat->packed_cols + mat->packed_ind[i];
	    pd = mat->packed_data + mat->packed_ind[i];
	    n = (INT)(mat->packed_ind[i + 1] - mat->packed_ind[i]);
	    for (k = 0; k < n; k++) {
		if (pc[k] == i)
		    break;
	    }
	    mat->diag[i] = (k < n ? pd[k] : 0.0);
	}
    }
}

VEC *
phgMatVec(MAT_OP op, FLOAT alpha, MAT *A, VEC *x, FLOAT beta, VEC **y_ptr)
/* BLAS gemv like function computing y := alpha op(A) x + beta y.
 * Special values of alpha and beta are taken into account to reduce
 * the number of computations. if A == NULL, then A is considered
 * to be the identity matrix I of size x->nglobal.
 *
 * op should be one of:
 *	MAT_OP_N (multiply x by A),
 *	MAT_OP_T (multiply x by trans(A)),
 *	MAT_OP_D (multiply x by the diagonal part of A)
 *
 * This function can do many kinds of matrix-vector or vector-vector
 * operations, for example:
 *
 *  phgMatVec(MAT_OP_N, 1.0, NULL, x, 0.0, &y)       . y := x  (vector copy)
 *  phgMatVec(MAT_OP_N, -1.0, NULL, x, 0.0, &y)	     . y := - x
 *  phgMatVec(MAT_OP_N, alpha, NULL, x, 0.0, &y)     . y := alpha x
 *  phgMatVec(MAT_OP_N, 0.0, NULL, NULL, -1.0, &y)   . y := - y
 *  phgMatVec(MAT_OP_N, 0.0, NULL, NULL, beta, &y)   . y := beta y
 *  phgMatVec(MAT_OP_N, alpha, NULL, x, beta, &y)    . y := alpha x + beta y
 *  phgMatVec(MAT_OP_N, 1.0, A, x, 0.0, &y)	     . y := A x
 *  phgMatVec(MAT_OP_T, 1.0, A, x, 0.0, &y)          . y := A' x
 *  phgMatVec(MAT_OP_N, -1.0, A, x, 0.0, &y)	     . y := - A x
 */
{
    INT i, j, n, *pc;
    FLOAT *data, *pd, *v, *v0;
    VEC *y = (y_ptr == NULL ? NULL : *y_ptr);
#if USE_MPI
    FLOAT *offp_data = NULL;
#endif	/* USE_MPI */

    MagicCheck(VEC, y)

    if (alpha != 0.0) {
	if (x == NULL)
	    phgError(1,
		     "phgMatVec, %s:%d, x must be non NULL when alpha!=0\n",
		     __FILE__, __LINE__);
    }
    else {
	if (y == NULL)
	    phgError(1,
		     "phgMatVec, %s:%d, y must be non NULL when alpha==0\n",
		     __FILE__, __LINE__);
    }

    assert(y == NULL || y->nvec == 1);
    assert(x == NULL || x->nvec == 1);

    if (x != NULL && !x->assembled)
	phgVecAssemble(x);

    if (y != NULL && !y->assembled)
	phgVecAssemble(y);

    /* multiply y by beta */
    if (alpha == 0.0 || A != NULL) {
	if (y == NULL) {
	    /* x must be non NULL */
	    if (A != NULL) {
		if (op != MAT_OP_T)
		    y = phgMapCreateVec(A->rmap, x->nvec);
		else
		    y = phgMapCreateVec(A->cmap, x->nvec);
	    }
	    else
		y = phgMapCreateVec(x->map, x->nvec);

	    if (y_ptr != NULL)
		*y_ptr = y;
	    beta = 0.0;
	}
	else {
	    /* compute y := beta y */
	    if (beta == 0.0) {
		bzero(y->data, y->map->nlocal * sizeof(*y->data));
	    }
	    else if (beta == (FLOAT)(-1.0)) {
		for (i = 0, v = y->data; i < y->map->nlocal; i++, v++)
		    *v = -(*v);
	    }
	    else if (beta != (FLOAT)1.0) {
		for (i = 0, v = y->data; i < y->map->nlocal; i++)
		    *(v++) *= beta;
	    }
	}

	if (alpha == 0.0)
	    return y;
    }
    else {			/* alpha != 0 && A == NULL */
	if (y == NULL) {
	    y = phgMapCreateVec(x->map, x->nvec);
	    if (y_ptr != NULL)
		*y_ptr = y;
	    beta = 0.0;
	}
	else {
	    vec_compare(x, y);
	}
	if (alpha == 1.0) {
	    if (beta == 0.0) {	/* y := x */
		memcpy(y->data, x->data, y->map->nlocal * sizeof(*y->data));
	    }
	    else if (beta == (FLOAT)(-1.0)) {	/* y := x - y */
		v0 = x->data;
		for (i = 0, v = y->data; i < y->map->nlocal; i++, v++)
		    *v = *(v0++) - (*v);
	    }
	    else if (beta == (FLOAT)1.0) {	/* y := x + y */
		v = y->data; v0 = x->data;
#if USE_OMP
                if (phgMaxThreads == 1) {
#endif  /* USE_OMP */
		for (i = 0; i < y->map->nlocal; i++)
		    *(v++) += (*(v0++));
#if USE_OMP 
                }else {
#pragma omp parallel for schedule(static)  default(shared)
		  for (i = 0; i < y->map->nlocal; i++)
		    v[i] += v0[i];
                }
#endif /* USE_OMP */
	    }
	    else {		/* y := x + beta y */
#if USE_OMP
                if (phgMaxThreads == 1) {
#endif  /* USE_OMP */
		    for (i = 0, v = y->data, v0 = x->data; i < y->map->nlocal;
		         i++, v++)
		        *v = *(v0++) + beta * (*v);
#if USE_OMP
                }else {
	            v = y->data; v0 = x->data;
#pragma omp parallel for schedule(static)  default(shared)
	 	    for (i = 0; i < y->map->nlocal; i++)
                        v[i] = v0[i] + beta * v[i];
                }
#endif  /* USE_OMP */
	    }
	}
	else if (alpha == -1.0) {
	    if (beta == 0.0) {	/* y := - x */
		for (i = 0, v = y->data, v0 = x->data; i < y->map->nlocal; i++)
		    *(v++) = -*(v0++);
	    }
	    else if (beta == (FLOAT)(-1.0)) {	/* y := - x - y */
		v0 = x->data;
		for (i = 0, v = y->data; i < y->map->nlocal; i++, v++)
		    *v = -*(v0++) - (*v);
	    }
	    else if (beta == (FLOAT)1.0) {	/* y := - x + y */
		for (i = 0, v = y->data, v0 = x->data; i < y->map->nlocal; i++)
		    *(v++) += -*(v0++);
	    }
	    else {		/* y := - x + beta y */
#if USE_OMP
                if (phgMaxThreads == 1) {
#endif  /* USE_OMP */
		    for (i = 0, v = y->data, v0 = x->data; i < y->map->nlocal;
		         i++, v++)
		        *v = -*(v0++) + beta * (*v);
#if USE_OMP
                    }else {
		        v = y->data; v0 = x->data;
#pragma omp parallel for schedule(static) default(shared)
		        for (i = 0; i < y->map->nlocal; i++)
		            v[i] = -v0[i] + beta * v[i];
                    }
#endif  /* USE_OMP */
	    }
	}
	else {
	    if (beta == 0.0) {	/* y := alpha x */
		for (i = 0, v = y->data, v0 = x->data; i < y->map->nlocal; i++)
		    *(v++) = alpha * *(v0++);
	    }
	    else if (beta == (FLOAT)(-1.0)) {	/* y := alpha x - y */
		for (i = 0, v = y->data, v0 = x->data; i < y->map->nlocal;
		     i++, v++)
		    *v = alpha * *(v0++) - (*v);
	    }
	    else if (beta == (FLOAT)1.0) {	/* y := alpha x + y */
		for (i = 0, v = y->data, v0 = x->data; i < y->map->nlocal; i++)
		    *(v++) += alpha * (*(v0++));
	    }
	    else {		/* y := alpha x + beta y */
#if USE_OMP
                if (phgMaxThreads == 1) {
#endif  /* USE_OMP */
		    for (i = 0, v = y->data, v0 = x->data; i < y->map->nlocal;
		        i++, v++)
		       *v = alpha * *(v0++) + beta * (*v);
#if USE_OMP
                }else {
                    v = y->data; v0 = x->data;
#pragma omp parallel for schedule(static) default(shared)
		    for (i = 0; i < y->map->nlocal;i++)
		       v[i] = alpha * v0[i] + beta * v[i];
		}
#endif  /* USE_OMP */
	    }
	}

	return y;
    }

    /* now we must have alpha != 0 and A != NULL */
    assert(A->type != PHG_DESTROYED);

    if (!A->assembled)
	phgMatAssemble(A);

    if (A->type == PHG_MATRIX_FREE) {
	if (A->mv_b != 0.0) {
	    VEC *b;
	    int ret;
	    b = phgMapCreateVec(y->map, y->nvec);
	    ret = A->mv_func(op, A, x, b);
	    if (ret != 0)
		phgError(1, "phgMatVec: nonzero return (%d) from mv_func.\n",
				ret);
	    phgVecAXPBY(alpha * A->mv_b, b, 1.0, &y);
	    phgVecDestroy(&b);
	}

	if (A->mv_a != 0.0 && A->mv_x != NULL)
	    phgMatVec(op, alpha * A->mv_a, A->mv_x, x, 1.0, &y);

	return y;
    }

    if (op == MAT_OP_D) {
	assert(x->map->nlocal == y->map->nlocal);
	if (A->diag == NULL)
	    phgMatSetupDiagonal(A);
	pd = A->diag;
	v0 = x->data;
	v = y->data;
	if (alpha == 1.0) {
#if USE_OMP
	    if (phgMaxThreads == 1) {
#endif  /* USE_OMP */
	    for (i = 0; i < y->map->nlocal; i++)
		*(v++) += *(v0++) * *(pd++);
#if USE_OMP
	    }
	    else {
#pragma omp parallel for schedule(static) default(shared)
	    for (i = 0; i < y->map->nlocal; i++)
		v[i] += v0[i] * pd[i];
	    }
#endif  /* USE_OMP */
	}
	else if (alpha == -1.0) {
#if USE_OMP
            if (phgMaxThreads == 1) {
#endif  /* USE_OMP */
	    for (i = 0; i < y->map->nlocal; i++)
		*(v++) -= *(v0++) * *(pd++);
#if USE_OMP
            }
            else {
#pragma omp parallel for schedule(static) default(shared)
	    for (i = 0; i < y->map->nlocal; i++)
		v[i] -= v0[i] * pd[i];
            }
#endif  /* USE_OMP */
	}
	else {
#if USE_OMP
            if (phgMaxThreads == 1) {
#endif  /* USE_OMP */
	    for (i = 0; i < y->map->nlocal; i++)
		*(v++) += alpha * *(v0++) * *(pd++);
#if USE_OMP
            }
            else {
#pragma omp parallel for schedule(static) default(shared)
	    for (i = 0; i < y->map->nlocal; i++)
		v[i] += alpha*v0[i] * pd[i];
            }
#endif  /* USE_OMP */
	}
	return y;
    }

    phgMatPack(A);

    if (alpha == 1.0) {
	data = x->data;
    }
    else if (alpha == -1.0) {
	v = data = phgAlloc(x->map->nlocal * sizeof(*data));
	v0 = x->data;
	for (i = 0; i < x->map->nlocal; i++)
	    *(v++) = -*(v0++);
    }
    else {
	v = data = phgAlloc(x->map->nlocal * sizeof(*data));
	v0 = x->data;
#if USE_OMP
        if (phgMaxThreads == 1) {
#endif  /* USE_OMP */
           for (i = 0; i < x->map->nlocal; i++)
	       *(v++) = alpha * *(v0++);
#if USE_OMP
        }
        else {
#pragma omp parallel for schedule(static) default(shared)
           for (i = 0; i < x->map->nlocal; i++)
	       v[i] = alpha * v0[i];
            }
#endif  /* USE_OMP */
    }

#if USE_MPI
    if (A->cmap->nprocs > 1) {
	offp_data = phgAlloc(A->cinfo->rsize * sizeof(*offp_data));
    }
#endif	/* USE_MPI */

    if (op == MAT_OP_N) {	/* y += alpha * A * x */
	if (A->cmap->nlocal != x->map->nlocal ||
	    A->rmap->nlocal != y->map->nlocal)
	    phgError(1, "%s:%d: inconsistent matrix-vector.", __FILE__,
		     __LINE__);

#if USE_MPI
	if (A->cmap->nprocs > 1) {
	    phgMapScatterBegin(A->cinfo, y->nvec, data, offp_data);
	}
#endif	/* USE_MPI */

	/* multiply with local data */
#if USE_OMP
    if (phgMaxThreads == 1) {
#endif	/* USE_OMP */
	for (i = 0, v = y->data; i < A->rmap->nlocal; i++) {
	    pc = A->packed_cols + A->packed_ind[i];
	    pd = A->packed_data + A->packed_ind[i];
	    n = (INT)(A->packed_ind[i + 1] - A->packed_ind[i]);
	    if (n == 0) {
		v++;
		continue;
	    }
#if 0
	    /* manual loop unrolling */
	    beta = 0;
	    switch (j = (n%4)) {
		case 3:
		    beta += *(pd++) * data[*(pc++)];
		case 2:
		    beta += *(pd++) * data[*(pc++)];
		case 1:
		    beta += *(pd++) * data[*(pc++)];
	    }
	    for (; j < n; j += 4) {
		beta += *(pd++) * data[*(pc++)];
		beta += *(pd++) * data[*(pc++)];
		beta += *(pd++) * data[*(pc++)];
		beta += *(pd++) * data[*(pc++)];
	    }
#else	/* 0 */
#if USE_INDP
	    beta = *(pd++) * data[*(pc++)];
	    for (j = 1; j < n; j++)
		beta += *(pd++) * data[*(pc++)];
#else	/* USE_INDP */
	    beta = pd[0] * data[pc[0]];
	    for (j = 1; j < n; j++)
		beta += pd[j] * data[pc[j]];
#endif	/* USE_INDP */
#endif	/* 0 */
	    *(v++) += beta;
	}
#if USE_OMP
    }
    else {
	v = y->data;
#pragma omp parallel for private(i, j, beta, pd, pc, n) default(shared)
	for (i = 0; i < A->rmap->nlocal; i++) {
	    /*printf("thread %d/%d: i = %"dFMT"\n", phgThreadId, phgMaxThreads, i);*/
	    pc = A->packed_cols + A->packed_ind[i];
	    pd = A->packed_data + A->packed_ind[i];
	    n = (INT)(A->packed_ind[i + 1] - A->packed_ind[i]);
	    if (n == 0)
		continue;
#if USE_INDP
	    beta = *(pd++) * data[*(pc++)];
	    for (j = 1; j < n; j++)
		beta += *(pd++) * data[*(pc++)];
#else	/* USE_INDP */
	    beta = pd[0] * data[pc[0]];
	    for (j = 1; j < n; j++)
		beta += pd[j] * data[pc[j]];
#endif	/* USE_INDP */
	    v[i] += beta;
	}
    }
#endif	/* USE_OMP */

#if USE_MPI
	if (A->cmap->nprocs > 1) {
	    phgMapScatterEnd(A->cinfo, x->nvec, data, offp_data);
	    /* multiply with remote data */
#if USE_OMP
    if (phgMaxThreads == 1) {
#endif	/* USE_OMP */
	    for (i = 0, v = y->data; i < A->rmap->nlocal; i++) {
		j = A->rmap->nlocal + i;
		pc = A->packed_cols + A->packed_ind[j];
		pd = A->packed_data + A->packed_ind[j];
		n = (INT)(A->packed_ind[j + 1] - A->packed_ind[j]);
		if (n == 0) {
		    v++;
		    continue;
		}
#if USE_INDP
		beta = *(pd++) * offp_data[*(pc++)];
		for (j = 1; j < n; j++)
		    beta += *(pd++) * offp_data[*(pc++)];
#else	/* USE_INDP */
		beta = pd[0] * offp_data[pc[0]];
		for (j = 1; j < n; j++)
		    beta += pd[j] * offp_data[pc[j]];
		*(v++) += beta;
#endif	/* USE_INDP */
	    }
#if USE_OMP
    } else {
	    v = y->data;
#pragma omp parallel for private(i, j, beta, pd, pc, n) default(shared)
	    for (i = 0; i < A->rmap->nlocal; i++) {
		j = A->rmap->nlocal + i;
		pc = A->packed_cols + A->packed_ind[j];
		pd = A->packed_data + A->packed_ind[j];
		n = (INT)(A->packed_ind[j + 1] - A->packed_ind[j]);
		if (n == 0)
		    continue;
#if USE_INDP
		beta = *(pd++) * offp_data[*(pc++)];
		for (j = 1; j < n; j++)
		    beta += *(pd++) * offp_data[*(pc++)];
#else	/* USE_INDP */
		beta = pd[0] * offp_data[pc[0]];
		for (j = 1; j < n; j++)
		    beta += pd[j] * offp_data[pc[j]];
#endif	/* USE_INDP */
		v[i] += beta;
	    }
    }
#endif	/* USE_OMP */
	}
#endif	/* USE_MPI */
    }
    else if (op == MAT_OP_T) {		/* y += alpha * trans(A) * x */
#if USE_OMP
	/* Note: initialize buffer/bufferchunk to NULL/0 to make gcc6 happy */
	FLOAT *buffer = NULL;
	INT bufferchunk = 0;
        if (phgMaxThreads > 1) {
	    bufferchunk = A->cmap->nlocal;
#if USE_MPI
	    if (A->cmap->nprocs > 1 && bufferchunk < A->cinfo->rsize)
		bufferchunk = A->cinfo->rsize;
#endif	/* USE_MPI */
	    buffer = phgAlloc((phgMaxThreads - 1) * bufferchunk
						  * sizeof(*buffer));
	}
#endif  /* USE_OMP */

	if (A->cmap->nlocal != y->map->nlocal ||
	    A->cmap->nglobal != y->map->nglobal ||
	    A->rmap->nlocal != x->map->nlocal ||
	    A->rmap->nglobal != x->map->nglobal)
	    phgError(1, "%s:%d: inconsistent mat-vec.", __FILE__, __LINE__);

#if USE_MPI
	if (A->cmap->nprocs > 1 && A->rmap->nprocs == A->cmap->nprocs) {
	/* off-process entries */
#if USE_OMP
	    if (phgMaxThreads == 1) {
#endif	/* USE_OMP */
		bzero(offp_data, A->cinfo->rsize * sizeof(*offp_data));
		v = data;
		for (i = 0; i < A->rmap->nlocal; i++) {
		    j = A->rmap->nlocal + i;
		    pc = A->packed_cols + A->packed_ind[j];
		    pd = A->packed_data + A->packed_ind[j];
		    n = (INT)(A->packed_ind[j + 1] - A->packed_ind[j]);
		    if (n == 0) {
			v++;
			continue;
		    }
		    beta = *(v++);
		    for (j = 0; j < n; j++) {
			offp_data[pc[j]] += pd[j] * beta;
		    }
		}
#if USE_OMP
	    }
	    else {
		bufferchunk = A->cinfo->rsize;
		bzero(offp_data, bufferchunk * sizeof(*offp_data));
		bzero(buffer, (phgMaxThreads - 1) * bufferchunk
						  * sizeof(*buffer));
		v = data;
#pragma omp parallel default(shared) private(v0)
{
		if (phgThreadId == 0)
		    v0 = offp_data;
		else
		    v0 = buffer + (phgThreadId - 1) * bufferchunk;
/*		phgInfo(-1,"thread %d/%d: i = %"dFMT" rsize %"dFMT" %"dFMT"\n", phgThreadId, phgMaxThreads, i, A->cinfo->rsize, (INT)((v0-buffer)/bufferchunk)); */
#pragma omp for schedule(static) private(i, j, beta, pd, pc, n)
		for (i = 0; i < A->rmap->nlocal; i++) {
		    j = A->rmap->nlocal + i;
		    pc = A->packed_cols + A->packed_ind[j];
		    pd = A->packed_data + A->packed_ind[j];
		    n = (INT)(A->packed_ind[j + 1] - A->packed_ind[j]);
		    if (n == 0) {
			continue;
		    }
		    beta = v[i];
		    for (j = 0; j < n; j++) 
#if USE_INDP
			v0[*(pc++)] += *(pd++) * beta;
#else	/* USE_INDP */
			v0[pc[j]] += pd[j] * beta;
#endif	/* USE_INDP */
		}
#pragma omp for schedule(static) private(i, j)
		/* FIXME: more efficient to swap i and j loops within thread? */
		for (i = 0; i < bufferchunk; i++) {
		    for (j = 0; j < phgMaxThreads - 1; j++) {
		       offp_data[i] += buffer[i + j * bufferchunk];
		   }
		}
} /* omp section*/
	    }
#endif	/* USE_OMP */
	    /* exchange offp_data */
	    phgMapGatherBegin(A->cinfo, y->nvec, y->data, offp_data);
	}
#endif	/* USE_MPI */

	/* in-process entries */
#if USE_OMP
	if (phgMaxThreads == 1) {
#endif	/* USE_OMP */
	    v = data;
	    for (i = 0; i < A->rmap->nlocal; i++) {
		pc = A->packed_cols + A->packed_ind[i];
		pd = A->packed_data + A->packed_ind[i];
		n = (INT)(A->packed_ind[i + 1] - A->packed_ind[i]);
		if (n == 0) {
		    v++;
		    continue;
		}
		beta = *(v++);
		for (j = 0; j < n; j++) {
		    y->data[pc[j]] += pd[j] * beta;
		}
	    }
#if USE_OMP
	}
	else {
	    v = data;
	    bufferchunk = y->map->nlocal;
	    bzero(buffer, (phgMaxThreads - 1) * bufferchunk * sizeof(*buffer));
#pragma omp parallel default(shared) private(v0)
{
	    if (phgThreadId == 0)
		v0 = y->data;
	    else
		v0 = buffer + (phgThreadId - 1) * bufferchunk;
#pragma omp for schedule(static) private(i, j, beta, pd, pc, n)
	    for (i = 0; i < A->rmap->nlocal; i++) {
		pc = A->packed_cols + A->packed_ind[i];
		pd = A->packed_data + A->packed_ind[i];
		n = (INT)(A->packed_ind[i + 1] - A->packed_ind[i]);
		if (n == 0)
		    continue;
		beta = v[i];
		for (j = 0; j < n; j++) 
#if USE_INDP
		    v0[*(pc++)] += *(pd++) * beta;
#else	/* USE_INDP */
		    v0[pc[j]] += pd[j] * beta;
#endif	/* USE_INDP */
	    }
#pragma omp for schedule(static) private(i, j)
	    /* FIXME: more efficient to swap i and j loops within thread? */
	    for (i = 0; i < bufferchunk; i++){
		for (j = 0; j < phgMaxThreads - 1; j++)
		    y->data[i] += buffer[i + j * bufferchunk];
	    }
} /* omp section */
	}
#endif  /* USE_OMP */

#if USE_OMP
	if (phgMaxThreads > 1) {
	    phgFree(buffer);
	}
#endif

#if USE_MPI
	if (A->cmap->nprocs > 1 && A->rmap->nprocs == A->cmap->nprocs) {
	    phgMapGatherEnd(A->cinfo, y->nvec, y->data, offp_data);
	}
	if (A->cmap->nprocs == 1 && A->rmap->nprocs > 1) {
	    /* special case: distributed rows and serial columns */
	    v = y->data;
	    y->data = phgAlloc(y->nvec * y->map->nlocal * sizeof(*y->data));
	    MPI_Allreduce(v, y->data, y->nvec * y->map->nlocal,
			  PHG_MPI_FLOAT, PHG_SUM, A->rmap->comm);
	    phgFree(v);
	}
#endif	/* USE_MPI */

    }
    else {
	phgInfo(-1, "%s: op must be one of MAT_OP_N (%d), MAT_OP_T (%d) and "
		    "MAT_OP_D (%d)", __func__, MAT_OP_N, MAT_OP_T, MAT_OP_D);
	phgError(1, "%s: invalid MAT operation %d!\n", __func__, op);
    }

    if (data != x->data)
	phgFree(data);

#if USE_MPI
    phgFree(offp_data);
#endif	/* USE_MPI */

    return y;
}

FLOAT
phgVecDot(VEC *x, int which_x, VEC *y, int which_y, FLOAT *res)
{
    INT i;
    FLOAT v;
#if USE_MPI
    FLOAT v1;
#endif	/* USE_MPI */

    assert(x->nvec == 1 && y->nvec == 1);	/* to be removed */
    assert(which_x < x->nvec && which_y < y->nvec
		&& phgMapCompare(x->map, y->map));

    if (!x->assembled)
	phgVecAssemble(x);

    if (!y->assembled)
	phgVecAssemble(y);

    v = 0.0;
#if USE_OMP
    if (phgMaxThreads == 1) {
#endif  /* USE_OMP */
	for (i = 0; i < x->map->nlocal; i++)
	    v += x->data[i] * y->data[i];
#if USE_OMP
    }
    else {
#pragma omp parallel for schedule(static) default(shared) reduction(+:v)
	for (i = 0; i < x->map->nlocal; i++)
 	    v += x->data[i] * y->data[i];
    }
#endif  /* USE_OMP */

#if USE_MPI
    MPI_Allreduce(&v, &v1, 1, PHG_MPI_FLOAT, PHG_SUM, x->map->comm);
    v = v1;
#endif	/* USE_MPI */
    if (res != NULL)
	*res = v;

    return v;
}

VEC *
phgVecAXPBY(FLOAT a, VEC *x, FLOAT b, VEC **y_ptr)
{
    INT i, n;
    VEC *y = (y_ptr == NULL ? NULL : *y_ptr);

    MagicCheck(VEC, y)

    assert(x != NULL || y != NULL);
    assert((a == 0.0 || x != NULL) && (b == 0.0 || y != NULL));

    if (a != 0. && !x->assembled)
	phgVecAssemble(x);

    if (y == NULL) {
	y = phgMapCreateVec(x->map, x->nvec);
	if (y_ptr != NULL)
	    *y_ptr = y;
    }

    if (b != 0. && !y->assembled)
	phgVecAssemble(y);

    if (a == 0.0 && b == 0.0) {
	bzero(y->data, y->map->nlocal * y->nvec * sizeof(*y->data));
	return y;
    }

#if USE_OMP
    if (phgMaxThreads == 1) {
#endif  /* USE_OMP */
        if (b == -1.0) {
	    n = y->map->nlocal * y->nvec;
	    for (i = 0; i < n; i++)
	        y->data[i] = -y->data[i];
        }
        else if (b != 1.0) {
	    n = y->map->nlocal * y->nvec;
	    for (i = 0; i < n; i++)
	        y->data[i] = b * y->data[i];
        }

        if (a == 0.0)
	    return y;

        vec_compare(x, y);

        n = y->map->nlocal * y->nvec;
        if (a == 1.0) {
	    for (i = 0; i < n; i++)
	        y->data[i] += x->data[i];
        }
        else if (a == -1.0) {
	    for (i = 0; i < n; i++)
	        y->data[i] -= x->data[i];
        }
        else {
	    for (i = 0; i < n; i++)
	    y->data[i] += a * x->data[i];
        }
#if USE_OMP
    }
    else {
        if (b == -1.0) {
	    n = y->map->nlocal * y->nvec;
	    for (i = 0; i < n; i++)
	        y->data[i] = -y->data[i];
        }
        else if (b != 1.0) {
	    n = y->map->nlocal * y->nvec;
#pragma omp parallel for schedule(static) default(shared)
	    for (i = 0; i < n; i++)
	        y->data[i] = b * y->data[i];
        }

        if (a == 0.0)
	    return y;

        vec_compare(x, y);

        n = y->map->nlocal * y->nvec;
        if (a == 1.0) {
#pragma omp parallel for schedule(static) default(shared)
	    for (i = 0; i < n; i++)
	        y->data[i] += x->data[i];
        }
        else if (a == -1.0) {
#pragma omp parallel for schedule(static) default(shared)
	    for (i = 0; i < n; i++)
	        y->data[i] -= x->data[i];
        }
        else {
#pragma omp parallel for schedule(static) default(shared)
	    for (i = 0; i < n; i++)
	    y->data[i] += a * x->data[i];
        }

    }
#endif  /* USE_OMP */

    return y;
}

VEC *
phgVecCreate(MPI_Comm comm, INT m, INT M, int nvec)
{
    MAP *map = phgMapCreateSimpleMap(comm, m, M);
    VEC *v = phgMapCreateVec(map, nvec);
    phgMapDestroy(&map);

    return v;
}

void
phgVecDestroy(VEC **vec_ptr)
{
    VEC *vec = *vec_ptr;

    MagicCheck(VEC, vec)

    if (vec == NULL)
	return;

    phgFree(vec->data);
    phgFree(vec->offp_data);
    phgMapDestroy(&vec->map);
#if USE_MPI
    phgFree(vec->O2Gmap);
    phgFree(vec->ordering);
#endif	/* USE_MPI */
    phgFree(vec);

    *vec_ptr = NULL;
    return;
}

MAT *
phgMatCreate(MPI_Comm comm, INT m, INT M)
{
    MAP *map = phgMapCreateSimpleMap(comm, m, M);
    MAT *mat = phgMapCreateMat(map, map);
    phgMapDestroy(&map);

    return mat;
}

MAT *
phgMatCreateNonSquare(MPI_Comm comm, INT m, INT M, INT n, INT N)
{
    MAP *rmap = phgMapCreateSimpleMap(comm, m, M);
    MAP *cmap = phgMapCreateSimpleMap(comm, n, N);
    MAT *mat = phgMapCreateMat(rmap, cmap);

    phgMapDestroy(&rmap);
    phgMapDestroy(&cmap);

    return mat;
}

void
_phg_matrix_free_setup(MAT *mat, MV_FUNC mv_func, int ndata, void **data)
{
    if (mat->mv_data != NULL) {
	phgFree(mat->mv_data);
	mat->mv_data = NULL;
    }

    if (ndata > 0) {
	mat->mv_data = phgAlloc((ndata + 1) * sizeof(*mat->mv_data));
	memcpy(mat->mv_data, data, ndata * sizeof(*data));
	mat->mv_data[ndata] = NULL;
    }

    phgFree(mat->rows);
    mat->rows = NULL;
    mat->mv_func = mv_func;
    mat->mv_a = 0.0;
    mat->mv_b = 1.0;
    mat->mv_x = NULL;
    mat->type = PHG_MATRIX_FREE;
}

MAT *
phgMatCreateMatrixFree(MPI_Comm comm, INT m, INT M, MV_FUNC mv_func,
		       void *mv_data, ...)
{
    MAT *mat = phgMatCreate(comm, m, M);
    int ndata;
    void **data;

    if (mat == NULL)
	return NULL;

    GetVariableArgs(ndata, data, mv_data, void *);
    _phg_matrix_free_setup(mat, mv_func, ndata, data);
    phgFree(data);

    return mat;
}

MAT *
phgMatCreateMatrixFreeNonSquare(MPI_Comm comm, INT m, INT M, INT n, INT N,
			MV_FUNC mv_func, void *mv_data, ...)
{
    MAT *mat = phgMatCreateNonSquare(comm, m, M, n, N);
    int ndata;
    void **data;

    if (mat == NULL)
	return NULL;

    GetVariableArgs(ndata, data, mv_data, void *);
    _phg_matrix_free_setup(mat, mv_func, ndata, data);
    phgFree(data);

    return mat;
}

static int
block_mv_func(MAT_OP op, MAT *mat, VEC *x, VEC *y)
/* performs block matrix - vector multiply */
{
    int i, j, p, q;
    MAT *a;
    VEC *xtmp, *ytmp;
    FLOAT *xdata, *ydata, coeff;
    INT m, n;

    assert(mat->blocks != NULL);

    /* TODO: communicate off-process data together instead of block by block */

    p = mat->blocks->p;
    q = mat->blocks->q;

    switch (op) {
	case MAT_OP_D:
	    assert(p == q);
	    ydata = y->data;
	    xdata = x->data;
	    for (i = 0; i < p; i++) {
		a = mat->blocks->pmat[i * q + i];
		m = mat->blocks->nrows[i];
		n = mat->blocks->ncols[i];
		if (a == NULL) {
		    bzero(ydata, m * y->nvec * sizeof(*ydata));
		    xdata += n * x->nvec;
		    ydata += m * y->nvec;
		    continue;
		}
		xtmp = phgMapCreateVec(a->cmap, x->nvec);
		phgFree(xtmp->data);
		xtmp->data = xdata;
		ytmp = phgMapCreateVec(a->rmap, y->nvec);
		phgFree(ytmp->data);
		ytmp->data = ydata;
		coeff = mat->blocks->coeff[i * q + i];
		phgMatVec(op, coeff, a, xtmp, 0.0, &ytmp);
		xtmp->data = NULL;
		ytmp->data = NULL;
		phgVecDestroy(&xtmp);
		phgVecDestroy(&ytmp);
		xdata += n * x->nvec;
		ydata += m * y->nvec;
	    }
	    break;
	case MAT_OP_N:
	    ydata = y->data;
	    bzero(ydata, mat->rmap->nlocal * y->nvec * sizeof(*ydata));
	    for (i = 0; i < p; i++) {
		m = mat->blocks->nrows[i];
		ytmp = NULL;
		xdata = x->data;
		for (j = 0; j < q; j++) {
		    n = mat->blocks->ncols[j];
		    a = mat->blocks->pmat[i * q + j];
		    coeff = mat->blocks->coeff[i * q + j];
		    if (a == NULL || coeff == 0.0) {
			xdata += n * x->nvec;
			continue;
		    }
		    xtmp = phgMapCreateVec(a->cmap, x->nvec);
		    phgFree(xtmp->data);
		    xtmp->data = xdata;
		    if (ytmp == NULL) {
			ytmp = phgMapCreateVec(a->rmap, y->nvec);
			phgFree(ytmp->data);
			ytmp->data = ydata;
			phgMatVec(op, coeff, a, xtmp, 0.0, &ytmp);
		    }
		    else {
			phgMatVec(op, coeff, a, xtmp, 1.0, &ytmp);
		    }
		    xtmp->data = NULL;
		    phgVecDestroy(&xtmp);
		    xdata += n * x->nvec;
		}
		if (ytmp != NULL) {
		    ytmp->data = NULL;
		    phgVecDestroy(&ytmp);
		}
		ydata += m * y->nvec;
	    }
	    break;
	case MAT_OP_T:
	    ydata = y->data;
	    bzero(ydata, mat->cmap->nlocal * y->nvec * sizeof(*ydata));
	    for (j = 0; j < q; j++) {
		n = mat->blocks->ncols[j];
		ytmp = NULL;
		xdata = x->data;
		for (i = 0; i < p; i++) {
		    m = mat->blocks->nrows[i];
		    a = mat->blocks->pmat[i * q + j];
		    if (a == NULL) {
			xdata += m * x->nvec;
			continue;
		    }
		    coeff = mat->blocks->coeff[i * q + j];
		    xtmp = phgMapCreateVec(a->rmap, x->nvec);
		    phgFree(xtmp->data);
		    xtmp->data = xdata;
		    if (ytmp == NULL) {
			ytmp = phgMapCreateVec(a->cmap, y->nvec);
			phgFree(ytmp->data);
			ytmp->data = ydata;
			phgMatVec(op, coeff, a, xtmp, 0.0, &ytmp);
		    }
		    else {
			phgMatVec(op, coeff, a, xtmp, 1.0, &ytmp);
		    }
		    xtmp->data = NULL;
		    phgVecDestroy(&xtmp);
		    xdata += m * x->nvec;
		}
		if (ytmp != NULL) {
		    ytmp->data = NULL;
		    phgVecDestroy(&ytmp);
		}
		ydata += n * y->nvec;
	    }
	    break;
    }

    return 0;
}

MAT *
phgMatCreateBlockMatrix(MPI_Comm comm, int p, int q, MAT *pmat[], FLOAT coeff[],
			MAT_OP trans[])
/*
 * Creates a PxQ block matrix with block (j,k) = coeff[j*q+k] * pmat[j*q+k].
 *
 * Notes:
 *	1) block-matrix is handled as special matrix-free matrix,
 *	2) a NULL entry in pmat[] corresponds to a zero matrix,
 *	3) coeff == NULL means all coefficients are equal to 1,
 *	4) trans == NULL means all blocks are non-transposed.
 */
{
    MAT *a;
    INT k, K, m, M, n, N, o;
    INT *nrows, *ncols, *noffp, *nrows_global, *ncols_global, **partitions;
    int i, j;
    int assembled = TRUE, handle_bdry_eqns = -1;
#if USE_MPI
    MPI_Comm c;
    int ret;
#endif	/* USE_MPI */

    if (trans != NULL)
	phgError(1, "%s (%s:%d): unimplemented case.\n", __func__,
			__FILE__, __LINE__);

    nrows = phgAlloc(p * sizeof(*nrows));
    nrows_global = phgAlloc(p * sizeof(*nrows_global));
    ncols = phgAlloc(q * sizeof(*ncols));
    ncols_global = phgAlloc(q * sizeof(*ncols_global));
    noffp = phgAlloc(q * sizeof(*noffp));
    partitions = phgAlloc(q * sizeof(*partitions));

    /* count number of rows */
    m = M = 0;
    for (i = 0; i < p; i++) {
	k = K = 0;
#if USE_MPI
	c = MPI_COMM_NULL;
#endif	/* USE_MPI */
	for (j = 0; j < q; j++) {
	    if ((a = pmat[i * q + j]) == NULL)
		continue;
	    if (!a->assembled)
		assembled = FALSE;
	    if (j == i) {
		/* diagonal block */
		if (handle_bdry_eqns < 0)
		    handle_bdry_eqns = a->handle_bdry_eqns;
		else if (handle_bdry_eqns != a->handle_bdry_eqns)
		    phgError(1, "%s: handle_bdry_eqns flags mismatch on "
				"diagonal blocks.\n", __func__);
	    }
	    else {
		if (a->handle_bdry_eqns)
		    phgError(1, "%s: handle_bdry_eqns is TRUE on non "
				"diagonal blocks.\n", __func__);
	    }
#if USE_MPI
	    if (c == MPI_COMM_NULL) {
		c = a->rmap->comm;
	    }
	    else {
		MPI_Comm_compare(c, a->rmap->comm, &ret);
		if (ret != MPI_IDENT)
		    phgError(1, "%s: inconsistent rmaps.\n", __func__);
	    }
#endif	/* USE_MPI */
	    a->refcount++;
	    if (K == 0) {
		k = a->rmap->nlocal;
		K = a->rmap->nglobal;
	    }
	    else {
		if (a->rmap->nlocal != k || a->rmap->nglobal != K)
		    phgError(1, "%s: inconsistent row blocks.\n", __func__);
	    }
	}
	assert(K > 0);
	nrows[i] = k;
	nrows_global[i] = K;
	m += k;
	M += K;
    }

    /* count number of columns */
    n = N = 0;
    for (j = 0; j < q; j++) {
	o = k = K = 0;
#if USE_MPI
	c = MPI_COMM_NULL;
#endif	/* USE_MPI */
	for (i = 0; i < p; i++) {
	    if ((a = pmat[i * q + j]) == NULL)
		continue;
#if USE_MPI
	    if (c == MPI_COMM_NULL) {
		c = a->cmap->comm;
	    }
	    else {
		MPI_Comm_compare(c, a->cmap->comm, &ret);
		if (ret != MPI_IDENT)
		    phgError(1, "%s: inconsistent cmaps.\n", __func__);
	    }
#endif	/* USE_MPI */
	    if (K == 0) {
		k = a->cmap->nlocal;
		K = a->cmap->nglobal;
		o = a->cmap->localsize - a->cmap->nlocal;
		partitions[j] = a->cmap->partition;
	    }
	    else {
		if (a->cmap->nlocal != k || a->cmap->nglobal != K ||
		    a->cmap->localsize - a->cmap->nlocal != o)
		    phgError(1, "%s: inconsistent row blocks.\n", __func__);
	    }
	}
	assert(K > 0);
	ncols[j] = k;
	ncols_global[j] = K;
	noffp[j] = o;
	n += k;
	N += K;
    }

    if (m == n && M == N)
	a = phgMatCreateMatrixFree(comm, m, M, block_mv_func, NULL);
    else
	a = phgMatCreateMatrixFreeNonSquare(comm, m, M, n, N, block_mv_func,
					    NULL);

    a->handle_bdry_eqns = handle_bdry_eqns;
    if (!assembled || handle_bdry_eqns)
	a->assembled = FALSE;

    a->blocks = phgAlloc(sizeof(*a->blocks));
    a->blocks->p = p;
    a->blocks->q = q;
    a->blocks->nrows = nrows;
    a->blocks->nrows_global = nrows_global;
    a->blocks->ncols = ncols;
    a->blocks->ncols_global = ncols_global;
    a->blocks->noffp = noffp;
    a->blocks->partitions = partitions;

    a->blocks->pmat = phgAlloc(p * q * sizeof(*a->blocks->pmat));
    memcpy(a->blocks->pmat, pmat, p * q * sizeof(*a->blocks->pmat));

    a->blocks->coeff = phgAlloc(p * q * sizeof(*a->blocks->coeff));
    if (coeff != NULL) {
	memcpy(a->blocks->coeff, coeff, p * q * sizeof(*a->blocks->coeff));
    }
    else {
	for (i = 0; i < p * q; i++)
	    a->blocks->coeff[i] = 1.0;
    }

    return a;
}

static INT
get_matrix_row(MAT *mat, INT irow, const INT **pcols, const FLOAT **pdata)
/* returns global column indices and data of 'irow'-th row.
 * Note: this function uses a stack of buffers for handling nested block
 *	 matrices, which are freed at exit. So the calling function need not
 *	 and should not free *pcols and *pdata */
{
    MAT **m;
    MAT_ROW *r;
    FLOAT *c;
    int ii, jj, rank;
    INT i, j, n, ncols, offset, n0, n1, *cols;
    FLOAT *data;
    const INT *cols0;
    const FLOAT *data0;

#ifndef ROW_STACK_DEPTH
# define ROW_STACK_DEPTH 8
#endif	/* !defined(ROW_STACK_DEPTH) */
    static int stack = -1;
    static FLOAT *d_buffer[ROW_STACK_DEPTH];
    static INT *c_buffer[ROW_STACK_DEPTH];
    static INT c_alloc[ROW_STACK_DEPTH], d_alloc[ROW_STACK_DEPTH];
#if USE_OMP
# pragma omp threadprivate(stack,d_buffer,c_buffer,c_alloc,d_alloc)
#endif	/* USE_OMP */

    assert(irow >= 0 && irow < mat->rmap->nlocal);

    if (stack < 0) {
	/* initialize the stack */
	stack = 0;
	bzero(d_buffer, sizeof(d_buffer));
	bzero(c_buffer, sizeof(c_buffer));
	bzero(c_alloc, sizeof(c_alloc));
	bzero(d_alloc, sizeof(d_alloc));
	for (ii = 0; ii < ROW_STACK_DEPTH; ii++) {
	    FreeAtExitNoCheck(c_buffer[ii]);
	    FreeAtExitNoCheck(d_buffer[ii]);
	}
    }

    if (stack > 1) {
	/* Note: nested block matrices */
	phgWarning("(%s:%d) untested case!\n", __FILE__, __LINE__);
    }

    if (mat->type != PHG_MATRIX_FREE) {
	/* ordinary matrix */
	if (!(mat->type == PHG_PACKED || mat->type == PHG_UNPACKED)) {
	    if (mat->type == PHG_DESTROYED)
		phgError(-1, "%s: the matrix has been destroyed.\n", __func__);
	    if (mat->type == PHG_DENSE)
		phgError(-1, "%s:%d: dense matrix unsupported.\n", __func__);
	}
	if (mat->type == PHG_UNPACKED) {
	    r = mat->rows + irow;
	    *pdata = r->data;
	    if (r->ncols == 0)
		*pcols = NULL;
#if USE_MPI
	    else if (mat->cmap->nprocs > 1) {
		/* convert to global column indices */
		if (r->ncols > c_alloc[stack]) {
		    c_alloc[stack] = r->ncols;
		    phgFree(c_buffer[stack]);
		    c_buffer[stack] = phgAlloc(r->ncols * sizeof(*r->cols));
		}
		*pcols = cols = c_buffer[stack];
		offset = mat->cmap->partition[mat->cmap->rank];
		for (i = 0; i < r->ncols; i++) {
		    if (r->cols[i] < mat->cmap->nlocal)
			cols[i] = offset + r->cols[i];
		    else
			cols[i] = mat->O2Gmap[r->cols[i] - mat->cmap->nlocal];
		}
	    }
	    else
#endif	/* USE_MPI */
		*pcols = r->cols;
	    return r->ncols;
	}
	else {		/* PHG_PACKED */
	    *pcols = mat->packed_cols + mat->packed_ind[irow];
	    *pdata = mat->packed_data + mat->packed_ind[irow];
	    ncols = (INT)(mat->packed_ind[irow + 1] - mat->packed_ind[irow]);
#if USE_MPI
	    n = irow + mat->rmap->nlocal;
	    n = (INT)(mat->packed_ind[n + 1] - mat->packed_ind[n]);
#else	/* USE_MPI */
	    n = 0;
#endif	/* USE_MPI */
	    if (mat->cmap->nprocs == 1) {
		return ncols;
	    }
#if USE_MPI
	    /* convert to global column indices */
	    if (ncols + n > c_alloc[stack]) {
		c_alloc[stack] = ncols + n;
		phgFree(c_buffer[stack]);
		c_buffer[stack] = phgAlloc((ncols + n) * sizeof(*cols));
	    }
	    cols = c_buffer[stack];
	    offset = mat->cmap->partition[mat->cmap->rank];
	    for (i = 0; i < ncols; i++)
		cols[i] = offset + (*pcols)[i];
	    if (n == 0) {
		*pcols = cols;
		return ncols;
	    }
	    *pcols = mat->packed_cols
			+ mat->packed_ind[irow + mat->rmap->nlocal];
	    for (i = 0; i < n; i++)
		cols[ncols + i] = mat->O2Gmap[mat->ordering[(*pcols)[i]]];
	    *pcols = cols;
	    /* copy data to buffer */
	    if (ncols + n > d_alloc[stack]) {
		d_alloc[stack] = ncols + n;
		phgFree(d_buffer[stack]);
		d_buffer[stack] = phgAlloc((ncols + n) * sizeof(*data));
	    }
	    if (ncols > 0)
		memcpy(d_buffer[stack], *pdata, ncols * sizeof(*data));
	    *pdata = mat->packed_data
			+ mat->packed_ind[irow + mat->rmap->nlocal];
	    memcpy(d_buffer[stack] + ncols, *pdata, n * sizeof(*data));
	    *pdata = d_buffer[stack];
#endif	/* USE_MPI */
	    return ncols + n;
	}
    }

    /* TODO: handle mv_a, mv_b and mv_x! */
    assert(mat->mv_a == 0.0 && mat->mv_x == NULL && mat->mv_b == 1.0);
    assert(mat->blocks != NULL);

    /* find the row block containing 'irow'->th row */
    m = mat->blocks->pmat;
    c = mat->blocks->coeff;
    for (ii = 0; ii < mat->blocks->p; ii++) {
	if (irow < mat->blocks->nrows[ii])
	    break;
	irow -= mat->blocks->nrows[ii];
	m += mat->blocks->q;
	c += mat->blocks->q;
    }
    assert(ii < mat->blocks->p);

    /* loop on the row of blocks */
    ncols = 0;
    cols = c_buffer[stack];
    data = d_buffer[stack];
    offset = mat->cmap->partition[mat->cmap->rank];
    for (ii = 0; ii < mat->blocks->q; ii++, m++, c++) {
	if (*m == NULL || *c == 0.0) {
	    offset += mat->blocks->ncols[ii];
	    continue;
	}
	stack++;
	assert(stack < ROW_STACK_DEPTH);
	n = get_matrix_row(*m, irow, &cols0, &data0);
	stack--;
	if (n > 0) {
	    if (ncols + n > c_alloc[stack]) {
		cols = phgRealloc_(cols, (ncols + n) * sizeof(*cols),
					 c_alloc[stack] * sizeof(*cols));
		c_alloc[stack] = ncols + n;
		c_buffer[stack] = cols;
	    }
	    if (ncols + n > d_alloc[stack]) {
		data = phgRealloc_(data, (ncols + n) * sizeof(*data),
					 d_alloc[stack] * sizeof(*data));
		d_alloc[stack] = ncols + n;
		d_buffer[stack] = data;
	    }
	    /* compute column indices */
	    n0 = (*m)->cmap->partition[(*m)->cmap->rank];
	    n1 = n0 + (*m)->cmap->nlocal;
	    for (i = 0; i < n; i++) {
		j = cols0[i];
		/* map global index j of *m to global index of mat */
		if (j >= n0 && j < n1) {
		    cols[ncols + i] = offset + j - n0;
		    continue;
		}
		/* TODO: replace the following map with a binary search */
		for (rank = 0; rank < (*m)->cmap->nprocs; rank++)
		    if (j < (*m)->cmap->partition[rank + 1])
			break;
		assert(rank < (*m)->cmap->nprocs);
		/* compute local index (with resp. to *m) in proc rank */
		j -= (*m)->cmap->partition[rank];
		/* TODO: pre-compute and store partial sums below */
		for (jj = 0; jj < ii; jj++)
		    j += mat->blocks->partitions[jj][rank + 1] -
			 mat->blocks->partitions[jj][rank];
		cols[ncols + i] = j + mat->cmap->partition[rank];
	    }
	    /* compute data */
	    if (*c == 1.0) {
		memcpy(data + ncols, data0, n * sizeof(*data));
	    }
	    else if (*c == -1.0) {
		for (i = 0; i < n; i++)
		    data[ncols + i] = -data0[i];
	    }
	    else {
		for (i = 0; i < n; i++)
		    data[ncols + i] = (*c) * data0[i];
	    }
	    ncols += n;
	}
	offset += mat->blocks->ncols[ii];
    }

    *pcols = cols;
    *pdata = data;

    return ncols;
}

const MAT_ROW *
phgMatGetRow_(MAT *mat, INT irow, MAT_ROW *row)
/* returns the 'irow'-th row (where 'irow' is the local vector index) of 'mat'.
 *
 * Note: let 'row' be the returned pointer, then
 *
 *   1. row->cols[] contains GLOBAL column indices.
 *
 *   2. If row == NULL, then row, row->cols and row->data point to internal
 *	buffers which need/should not be freed by user.
 *
 *   3. If row != NULL then the results are returned in row.
 *	At input row->cols and row->data should be valid buffers of size
 *	row->alloc, which may be reallocated when needed with row->alloc
 *	updated accordingly. The calling function is resposible to free
 *	row->cols and row->data after usage to avoid memory leaks.
 *
 * WARNING: when called with row == NULL, the data of the row obtained in the
 *	    previous call will be destroyed and replaced with the data of the
 *	    new call.
 */
{
    static MAT_ROW row_buffer = {NULL, NULL, 0, 0};
#if USE_OMP
# pragma omp threadprivate(row_buffer)
#endif	/* USE_OMP */
    MAT_ROW *r = &row_buffer;

    assert(irow >= 0 && irow < mat->rmap->nlocal);

    if (mat != NULL && !mat->assembled)
	phgMatAssemble(mat);

    /* WARNING: the function get_matrix_row() will modify cols/data to point to
     * some internal buffers. The calling function should not free them. */
    r->ncols = get_matrix_row(mat, irow, (void *)&r->cols, (void *)&r->data);

    if (row == NULL) {
	row = r;
    }
    else {
	/* copy data to row */
	if (row->alloc < r->ncols) {
	    row->alloc = r->ncols;
	    phgFree(row->cols);
	    row->cols = phgAlloc(row->alloc * sizeof(*row->cols));
	    phgFree(row->data);
	    row->data = phgAlloc(row->alloc * sizeof(*row->data));
	}
	row->ncols = r->ncols;
	if (row->ncols > 0) {
	    memcpy(row->cols, r->cols, r->ncols * sizeof(*row->cols));
	    memcpy(row->data, r->data, r->ncols * sizeof(*row->data));
	}
    }

    return row;
}

INT
phgMatGetNnzInRow(MAT *mat, INT irow, INT *nd, INT *no)
/* returns number of nonzeros in a given row */
{
    MAT **m;
    MAT_ROW *row;
    FLOAT *c;
    int ii;
    INT j, d, o, d1, o1;

    if (mat != NULL && !mat->assembled)
	phgMatAssemble(mat);

    d = o = 0;
    if (mat == NULL) {
	/* do nothing */
    }
    else if (mat->blocks == NULL) {
	assert(mat->type == PHG_PACKED || mat->type == PHG_UNPACKED);
	if (mat->type == PHG_UNPACKED) {
	    row = mat->rows + irow;
	    if (mat->cmap->nprocs == 1 || (nd == NULL && no == NULL)) {
		d = row->ncols;
	    }
	    else {
		for (j = 0; j < row->ncols; j++) {
		    row->cols[j] < mat->cmap->nlocal ? d++ : o++;
		}
	    }
	}
	else {
	    d = (INT)(mat->packed_ind[irow + 1] - mat->packed_ind[irow]);
#if USE_MPI
	    j = irow + mat->rmap->nlocal;
	    o = (INT)(mat->packed_ind[j + 1] - mat->packed_ind[j]);
#endif	/* USE_MPI */
	}
    }
    else {
	assert(mat->type == PHG_MATRIX_FREE);
	m = mat->blocks->pmat;
	c = mat->blocks->coeff;
	for (ii = 0; ii < mat->blocks->p; ii++) {
	    if (irow < mat->blocks->nrows[ii])
		break;
	    irow -= mat->blocks->nrows[ii];
	    m += mat->blocks->q;
	    c += mat->blocks->q;
	}
	assert(ii < mat->blocks->p);
	/* loop on the row of blocks */
	for (ii = 0; ii < mat->blocks->q; ii++, m++, c++) {
	    if (*c == 0.0)
		continue;
	    if (mat->cmap->nprocs == 1 || (nd == NULL && no == NULL)) {
		d += phgMatGetNnzInRow(*m, irow, NULL, NULL);
	    }
	    else {
		phgMatGetNnzInRow(*m, irow, &d1, &o1);
		d += d1;
		o += o1;
	    }
	}
    }

    if (nd != NULL)
	*nd = d;
    if (no != NULL)
	*no = o;

    return d + o;
}

FLOAT
phgMatGetEntry_(MAT *A, INT row, INT col, BOOLEAN global_col)
/* returns the entry of the matrix "A" with local row index "row" and
 * local (if "global_col" == FALSE) or global (if "global_col" == TRUE)
 * column index "col".
 */
{
    const MAT_ROW *R = phgMatGetRow(A, row);
    int i;

#if USE_MPI
    if (!global_col && A->cmap->nprocs > 1) {
	/* convert to global column index */
	assert(col >= 0 && col < A->cmap->localsize);
	if (col < A->cmap->nlocal)
	    col += A->cmap->partition[A->cmap->rank];
	else
	    col = A->cmap->O2Gmap[col - A->cmap->nlocal];
    }
#else
    Unused(global_col);
#endif	/* USE_MPI */

    for (i = 0; i < R->ncols; i++)
	if (col == R->cols[i])
	    return R->data[i];

    return 0.0;	/* not found */
} 

static void
free_matrix_rows(INT size, MAT_ROW **rows)
{
    MAT_ROW *row;

    if (*rows == NULL)
	return;

    /* free PHG's matrix */
    for (row = *rows; row < *rows + size; row++) {
	/* free PHG matrix row */
	phgFree(row->cols);
	phgFree(row->data);
	row->cols = NULL;
	row->data = NULL;
	row->ncols = row->alloc = 0;
    }
    phgFree(*rows);
    *rows = NULL;
}

static void
free_bdry_eqns(MAT *mat)
{
    INT i;

    if (mat->bdry_eqns == NULL)
	return;

    for (i = 0; i < mat->bdry_eqns_count; i++) {
	phgFree(mat->bdry_eqns[i].lu);
	phgFree(mat->bdry_eqns[i].pvt);
	phgFree(mat->bdry_eqns[i].map);
    }
    phgFree(mat->bdry_eqns);
    mat->bdry_eqns = NULL;
    mat->bdry_eqns_count = mat->bdry_eqns_allocated = 0;
}

#if USE_OMP
static void
free_locks(MAT *mat)
{
    INT i;

    CheckThread
    if (mat->locks != NULL) {
	for (i = 0; i < mat->nlock; i++) {
	    omp_destroy_lock(&mat->locks[i]);
	}
	phgFree(mat->locks);
	mat->nlock = -1;
	mat->locks = NULL;
    }
}
#else	/* USE_OMP */
#define free_locks(mat)
#endif	/* USE_OMP */

void
phgMatFreeMatrix(MAT *mat)
{
    int i;

    MagicCheck(MAT, mat)

    if (mat == NULL || mat->type == PHG_DESTROYED || mat->refcount > 0)
	return;

    if (mat->blocks != NULL) {
	/* block matrix */
	for (i = 0; i < mat->blocks->p * mat->blocks->q; i++)
	    phgMatFreeMatrix(mat->blocks->pmat[i]);
    }

    phgFree(mat->diag);
    mat->diag = NULL;
    phgFree(mat->diag1);
    mat->diag1 = NULL;
    phgFree(mat->packed_data);
    mat->packed_data = NULL;
    phgFree(mat->packed_cols);
    mat->packed_cols = NULL;
    phgFree(mat->packed_ind);
    mat->packed_ind = NULL;

    free_locks(mat);

    /*if (mat->type == PHG_MATRIX_FREE)
	return;*/

    phgInfo(1, "freeing matrix data.\n");
    mat->type = PHG_DESTROYED;

#if USE_MPI
    phgMapDestroyCommInfo(&mat->cinfo);
#endif	/* USE_MPI */

    free_matrix_rows(mat->rmap->localsize, &mat->rows);

    return;
}

static void
mat_reset(MAT *mat)
{
    int i;

    MagicCheck(MAT, mat)

    if (mat == NULL)
	return;

    if (mat->blocks != NULL) {
	/* block matrix */
	for (i = 0; i < mat->blocks->p * mat->blocks->q; i++)
	    phgMatDestroy(mat->blocks->pmat + i);
    }

    phgMatFreeMatrix(mat);

#if USE_MPI
    phgFree(mat->O2Gmap);
    mat->O2Gmap = NULL;
    phgFree(mat->ordering);
    mat->ordering = NULL;

    if (mat->cinfo != NULL)
	phgMapDestroyCommInfo(&mat->cinfo);
#endif	/* USE_MPI */

    free_matrix_rows(mat->rmap->nlocal, &mat->bdry_rows);
    free_bdry_eqns(mat);

    phgFree(mat->mv_data);
    mat->mv_data = NULL;
    phgMatDestroy(&mat->mv_x);

    phgFree(mat->bdry);
    mat->bdry = NULL;

    if (mat->blocks != NULL) {
	phgFree(mat->blocks->pmat);
	phgFree(mat->blocks->nrows);
	phgFree(mat->blocks->nrows_global);
	phgFree(mat->blocks->ncols);
	phgFree(mat->blocks->ncols_global);
	phgFree(mat->blocks->noffp);
	phgFree(mat->blocks->partitions);
	phgFree(mat->blocks->coeff);
	phgFree(mat->blocks);
	mat->blocks = NULL;
    }

    mat->nnz_d = mat->nnz_o = 0;
    mat->bdry_localsize = 0;

    mat->assembled = FALSE;
    mat->localsize = -1;
    mat->type = PHG_UNPACKED;

    return;
}

void
phgMatDestroy(MAT **mat_ptr)
{
    MAT *mat;

    if (mat_ptr == NULL)
	return;

    mat = *mat_ptr;
    MagicCheck(MAT, mat)

    if (mat == NULL)
	return;

    *mat_ptr = NULL;

    assert(mat->refcount >= 0);
    if (mat->refcount > 0) {
	mat->refcount--;
	return;
    }

    mat_reset(mat);

    phgMapDestroy(&mat->rmap);
    phgMapDestroy(&mat->cmap);

#if USE_OMP
    omp_destroy_lock(&mat->lock);
#endif	/* USE_OMP */

    phgFree(mat->row_buffer.cols);
    phgFree(mat->row_buffer.data);

    phgFree(mat);

    return;
}

MAP *
phgVecGetMap(VEC *vec)
{
    vec->map->refcount++;
    return vec->map;
}

MAP *
phgMatGetRowMap(MAT *mat)
{
    mat->rmap->refcount++;
    return mat->rmap;
}

MAP *
phgMatGetColumnMap(MAT *mat)
{
    mat->cmap->refcount++;
    return mat->cmap;
}

static void
dup_matrix_rows(INT nrows, MAT_ROW *row0, MAT_ROW **row_ptr)
{
    MAT_ROW *row = *row_ptr;
    INT i;

    if (row0 == NULL || nrows == 0)
	return;

    if (row == NULL)
	*row_ptr = row = phgCalloc(nrows, sizeof(*row));

    for (i = 0; i < nrows; i++, row++, row0++) {
	row->ncols = row->alloc = row0->ncols;
	phgFree(row->cols);
	phgFree(row->data);
	row->cols = phgAlloc(row->ncols * sizeof(*row->cols));
	row->data = phgAlloc(row->ncols * sizeof(*row->data));
	if (row->ncols > 0) {
	    memcpy(row->cols, row0->cols, row->ncols * sizeof(*row->cols));
	    memcpy(row->data, row0->data, row->ncols * sizeof(*row->data));
	}
    }
}

MAT *
phgMatDup(MAT *src)
{
    MAT *mat;
    INT i, n;
    size_t size;

    assert(src->type != PHG_DESTROYED);

    if (src->type == PHG_MATRIX_FREE) {
	src->refcount++;
	return src;
    }

    if (!src->assembled)
	phgMatAssemble(src);

#define DupData(ptr, len)						\
    if (src->ptr != NULL && (size = (len) * sizeof(*src->ptr)) != 0) {	\
	mat->ptr = phgAlloc(size);					\
	memcpy(mat->ptr, src->ptr, size);				\
    }

    mat = phgMapCreateMat(src->rmap, src->cmap);
    mat->type = src->type;
    dup_matrix_rows(src->rmap->nlocal, src->rows, &mat->rows);
    DupData(diag, mat->rmap->nlocal);
    mat->nnz_d = src->nnz_d;
    mat->nnz_o = src->nnz_o;
    DupData(packed_cols, src->nnz_d + src->nnz_o);
    DupData(packed_data, src->nnz_d + src->nnz_o);
    mat->localsize = src->localsize;
#if !USE_MPI
    DupData(packed_ind, src->rmap->nlocal + 1);
#else	/* !USE_MPI */
    DupData(packed_ind, 2 * src->rmap->nlocal + 1);
    n = mat->localsize - mat->cmap->nlocal;
    DupData(O2Gmap, n);
    DupData(ordering, n);
    if (src->cinfo != NULL) {
	mat->cinfo = phgCalloc(1, sizeof(*mat->cinfo));
	mat->cinfo->size = src->cinfo->size;
	mat->cinfo->ssize = src->cinfo->ssize;
	mat->cinfo->rsize = src->cinfo->rsize;
	mat->cinfo->comm = src->cinfo->comm;
	mat->cinfo->rind = src->cinfo->rind;
	DupData(cinfo->sind,  mat->cinfo->ssize);
	DupData(cinfo->scnts, mat->cmap->nprocs);
	DupData(cinfo->rcnts, mat->cmap->nprocs);
	DupData(cinfo->sdsps, mat->cmap->nprocs);
	DupData(cinfo->rdsps, mat->cmap->nprocs);
    }
#endif	/* !USE_MPI */
    mat->handle_bdry_eqns = src->handle_bdry_eqns;
    mat->bdry_localsize = src->bdry_localsize;
    DupData(bdry, mat->bdry_localsize);

    mat->bdry_eqns_allocated = mat->bdry_eqns_count = src->bdry_eqns_count;
    if (mat->bdry_eqns_count > 0) {
	mat->bdry_eqns = phgAlloc(mat->bdry_eqns_count *
					sizeof(*mat->bdry_eqns));
	for (i = 0; i < mat->bdry_eqns_count; i++) {
	    n = src->bdry_eqns[i].n;
	    DupData(bdry_eqns[i].lu, n * n);
	    DupData(bdry_eqns[i].pvt, n);
	    DupData(bdry_eqns[i].map, n);
	}
    }

    dup_matrix_rows(mat->rmap->nlocal, src->bdry_rows, &mat->bdry_rows);

    mat->assembled = TRUE;

#undef DupData

    return mat;
}

static INT digits = 16;		/* digits in MATLAB files */

static char *
add_dot_m(const char *var_name, const char *file_name, const char *func_name)
/* returns a dynamically allocated string for "file_name" with the ".m" ext. */
{
    char *fn;
    int n;

    if (file_name == NULL)
	file_name = var_name;

    n = strlen(file_name);
    if (n >= 2 && !strcmp(file_name + n - 2, ".m")) {
	fn = strdup(file_name);
    }
    else {
	fn = phgAlloc(n + 2 + 1);
	strcpy(fn, file_name);
	strcat(fn, ".m");
    }

    /*if (phgVerbosity > 0)*/
	phgPrintf("*** %s: creating \"%s\".\n", func_name, fn);

    return fn;
}

void
phgMatDumpMATLAB(MAT *A, const char *var_name, const char *file_name)
/* Outputs matrix A to the following MATLAB format:
 *	var_name = spconvert([1 1 5.0; 1 3 4.4; 2 4 5.5])
 * (a seperate data file is used for faster loading when matrix is large)
 *
 * Note: the smallest 8 generalized eigenvalues A x = lambda B x with sparse
 * 	 matrices A and B can be computed by eigs(A,B,8,'SM') */
{
    FILE *fp;
    INT i, j, k, rmax, cmax;
    const MAT_ROW *row;
    double nnz;
    INT n0;
    char format[16], *fn;
#if USE_MPI
    MPI_Status status;
    MPI_Datatype type;
# ifndef MAT_TRUNK_SIZE
#  define MAT_TRUNK_SIZE	65536
# endif	/* !defined(MAT_TRUNK_SIZE) */
    struct {
	FLOAT	a;
	INT	row, col;
    } *buffer;
    int n, tag = 111;
#endif	/* USE_MPI */

    if (!phgInitialized) {
	phgOptionsRegisterTitle("\nMATLAB options:", "\n", "matlab");
	phgOptionsRegisterInt("-matlab_digits", "Digits in output MATLAB files",
				&digits);
	return;
    }

    if (!A->assembled)
	phgMatAssemble(A);

    fn = add_dot_m(var_name, file_name, __func__);

    n0 = A->rmap->partition[A->rmap->rank];
    nnz = A->nnz_d + A->nnz_o;

    if (digits > 40)
	digits = 40;

    if (digits > 0)
	sprintf(format, "%%d %%d %%0.%dlg\n", (int)digits);
    else
	sprintf(format, "%%d %%d %%lg\n");

#if USE_MPI
    MPI_Type_contiguous(sizeof(*buffer), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    buffer = phgAlloc(MAT_TRUNK_SIZE * sizeof(*buffer));

    if (A->rmap->nprocs > 1) {
	double a;
	MPI_Reduce(&nnz, &a, 1, MPI_DOUBLE, MPI_SUM, 0, A->rmap->comm);
	nnz = a;
    }
#endif	/* USE_MPI */

    /*-----------------------------------------------------------------------*/
    if (A->rmap->rank == 0) {
	/* the root process (TODO: overlap communications with file I/O) */

	if ((fp = fopen(fn, "w")) == NULL)
	    phgError(1, "%s: can't create file \"%s\".\n", __func__, fn);

	fprintf(fp, "%% matrix size: %"dFMT"x%"dFMT", nnz: %lg\n",
		A->rmap->nglobal, A->cmap->nglobal, nnz);

	if (nnz <= 10000) {
	    fprintf(fp, "%s = spconvert([\n", var_name);
	}
	else {
	    char fndat[strlen(fn) + 1 + 4];
	    sprintf(fndat, "%s.dat", fn);
	    fprintf(fp, "%s_fid = fopen('%s', 'r');\n", var_name, fndat);
#if 0
    /* Note: may put commands and data in a same file by putting a 'return'
     * followed by data after the MATLAB commands, and skip the commands when
     * reading data (see code below).
     *
     * The problem is MATLAB still parses input after the 'return' command */
    fprintf(fp, "phg_tmp = '';\n"
		"while strcmp(phg_tmp, 'return') == 0\n"
		"    phg_tmp = fgetl(phg_fid);\n"
		"end\n");
#endif	/* 0 */
	    fprintf(fp, "%s = fscanf(%s_fid,'%%lf');\n", var_name, var_name);
	    fprintf(fp, "fclose(%s_fid);\n", var_name);
	    fprintf(fp, "%s = reshape(%s, 3, round(length(%s)/3));\n",
			var_name, var_name, var_name);
	    fprintf(fp, "%s = spconvert(%s');\n", var_name, var_name);
	    fclose(fp);
	    if ((fp = fopen(fndat, "w")) == NULL)
		phgError(1, "%s: can't create file \"%s\".\n", __func__, fndat);
	}

	rmax = cmax = 0;
	for (i = 0; i < A->rmap->nlocal; i++) {
	    row = phgMatGetRow(A, i);
	    for (j = 0; j < row->ncols; j++) {
		k = row->cols[j];
		fprintf(fp, format, i + n0 + 1, k + 1, (double)row->data[j]);
		if (rmax < i + n0 + 1)
		    rmax = i + n0 + 1;
		if (cmax < k + 1)
		    cmax = k + 1;
	    }
	}

#if USE_MPI
	/* receive data from other processes */
	for (k = 1; k < A->rmap->nprocs; k++) {
	    while (TRUE) {
		MPI_Probe(k, tag, A->rmap->comm, &status);
		MPI_Get_count(&status, type, &n);
		assert(n <= MAT_TRUNK_SIZE);
		MPI_Recv(buffer, n, type, k, tag, A->rmap->comm, &status);
		/* process received list */
		for (i = 0; i < n; i++) {
		    fprintf(fp, format, buffer[i].row, buffer[i].col,
						 (double)buffer[i].a);
		    if (rmax < buffer[i].row)
			rmax = buffer[i].row;
		    if (cmax < buffer[i].col)
			cmax = buffer[i].col;
		}
		if (n < MAT_TRUNK_SIZE)
		    break;
	    }
	}
#endif	/* USE_MPI */

	/* add a dummy entry to force matrix size */
	if (rmax < A->rmap->nglobal || cmax < A->cmap->nglobal)
	    fprintf(fp, "%"dFMT" %"dFMT" 0\n",
				A->rmap->nglobal, A->cmap->nglobal);

	if (nnz <= 10000)
	    fprintf(fp, "]);\n");

	fclose(fp);
    }
    /*-----------------------------------------------------------------------*/
#if USE_MPI
    else {
	/* non root processes */
	n = 0;
	for (i = 0; i < A->rmap->nlocal; i++) {
	    row = phgMatGetRow(A, i);
	    for (j = 0; j < row->ncols; j++) {
		k = row->cols[j];
		buffer[n].a = row->data[j];
		buffer[n].row = i + n0 + 1;
		buffer[n].col = k + 1;
		if (++n >= MAT_TRUNK_SIZE) {
		    MPI_Send(buffer, n, type, 0, tag, A->rmap->comm);
		    n = 0;
		}
	    }
	}
	/* send the last block (may be size 0), which also marks EOD */
	MPI_Send(buffer, n, type, 0, tag, A->rmap->comm);
    }

    phgFree(buffer);
    MPI_Type_free(&type);
#endif	/* USE_MPI */
    /*-----------------------------------------------------------------------*/
    
#if 0
    /* output dense format */
    FLOAT *a = phgAlloc(A->cmap->nglobal * sizeof(*a));
    if ((fp = fopen(fn, "w+t")) == NULL)
	phgError(1, "%s: cannot open \"%s\" for write.\n", __func__, fn);
    fprintf(fp, "%s = [\n", var_name);
    for (i = 0; i < A->rmap->nlocal; i++) {
	row = phgMatGetRow(A, i);
	bzero(a, A->cmap->nglobal * sizeof(*a));
	for (k = 0; k < row->ncols; k++)
	    a[row->cols[k]] = row->data[k];
	for (k = 0; k < A->cmap->nglobal; k++) {
	    fprintf(fp, "%0.16lg%s", (double)a[k],
			k < A->cmap->nglobal - 1 ? " " : "\n");
	}
    }
    fprintf(fp, "];\n");
    fclose(fp);
    phgFree(a);
#endif

    phgFree(fn);
}

void
phgVecDumpMATLAB(VEC *v, const char *var_name, const char *file_name)
/* outputs the vector components as columns of a MATLAB (dense) matrix */
{
    FILE *fp;
    INT i, j;
    char format[20], *fn;
#if USE_MPI
    int k, n, vec_trunk_size, tag = 222;
    MPI_Status status;
    FLOAT *buffer;
#endif	/* USE_MPI */

    if (v->nvec <= 0)
	return;
    if (!v->assembled)
	phgVecAssemble(v);

    fn = add_dot_m(var_name, file_name, __func__);

    if (digits > 40)
	digits = 40;
    
    if (digits > 0)
	sprintf(format, " %%0.%dlg", (int)digits);
    else
	sprintf(format, " %%lg");

#if USE_MPI
    /* Note: use a 8MB buffer */
    vec_trunk_size = 8 * 1024 * 1024 / sizeof(FLOAT) / v->nvec ;
#endif	/* USE_MPI */

    /*-----------------------------------------------------------------------*/
    if (v->map->rank == 0) {
	if ((fp = fopen(fn, v->map->rank == 0 ? "w" : "a")) == NULL)
	    phgError(1, "%s: can't create file \"%s\".\n", __func__, fn);
	fprintf(fp, "%% vector size: %"dFMT", vector components: %d\n",
			v->map->nglobal, v->nvec);
	fprintf(fp, "%s = [\n", var_name);
	for (i = 0; i < v->map->nlocal; i++) {
	    for (j = 0; j < v->nvec; j++)
		fprintf(fp, format, (double)v->data[j * v->map->nlocal + i]);
	    fprintf(fp, "\n");
	}

#if USE_MPI
	/* receive data from other processes */
	buffer = phgAlloc(vec_trunk_size * v->nvec * sizeof(FLOAT));
	for (k = 1; k < v->map->nprocs; k++) {
	    while (TRUE) {
		MPI_Probe(k, tag, v->map->comm, &status);
		MPI_Get_count(&status, PHG_MPI_FLOAT, &n);
		assert(n % v->nvec == 0);
		assert(n <= vec_trunk_size * v->nvec);
		MPI_Recv(buffer, n, PHG_MPI_FLOAT, k, tag, v->map->comm,
			 &status);
		/* process received list */
		n /= v->nvec;
		for (i = 0; i < n; i++) {
		    for (j = 0; j < v->nvec; j++)
			fprintf(fp, format, (double)buffer[j * n + i]);
		    fprintf(fp, "\n");
		}
		if (n < vec_trunk_size)
		    break;
	    }
	}
	phgFree(buffer);
#endif	/* USE_MPI */

	fprintf(fp, "];\n");
	fclose(fp);
    }
#if USE_MPI
    /*-----------------------------------------------------------------------*/
    else {
	i = 0;
	while (TRUE) {
	    n = v->map->nlocal - i;
	    if (n > vec_trunk_size)
		n = vec_trunk_size;
	    if (n == 0) {
		MPI_Send(v->data + i, n, PHG_MPI_FLOAT, 0, tag, v->map->comm);
	    }
	    else {
		MPI_Datatype type;
		/* Note: data are sent in the order [nvec][n] */
		MPI_Type_vector(v->nvec, n, v->map->nlocal, PHG_MPI_FLOAT,
				&type);
		MPI_Type_commit(&type);
		MPI_Send(v->data + i, 1, type, 0, tag, v->map->comm);
		MPI_Type_free(&type);
	    }
	    i += n;
	    if (n < vec_trunk_size)
		break;
	}
    }
    /*-----------------------------------------------------------------------*/
#endif	/* USE_MPI */

    phgFree(fn);
}

VEC *
phgVecRandomize(VEC *v, long int seed)
{
    INT i, n;
    FLOAT *p;

    if (v->nvec == 0 || v->map->nlocal == 0)
	return v;
    if (seed != 0)
	srand48(seed + v->map->rank);
    n = v->nvec * v->map->nlocal;
    if (v->data == NULL)
	v->data = phgAlloc(n * sizeof(*v->data));
    for (i = 0, p = v->data; i < n; i++)
	*(p++) = 2.0 * drand48() - 1.0;

    v->assembled = TRUE;

    return v;
}

/*-------------------------------------------------------------------------*/

#ifndef TRUNK_SIZE
# define	TRUNK_SIZE	8	/* size of allocation trunks */
#endif	/* !defined(TRUNK_SIZE) */

/* FIXME: SORT_COLUMNS==1 has not been tested when cols[] contains negated
 *	 global indices */
#ifndef SORT_COLUMNS
# define SORT_COLUMNS	0	/* whether sort columns */
#endif	/* !defined(SORT_COLUMNS) */

#if SORT_COLUMNS
static INT *cols_tmp;
#if USE_OMP
# pragma omp threadprivate(cols_tmp)
#endif	/* USE_OMP */

static int
comp_columns(const void *c0, const void *c1)
{
    INT i;

    return (i = cols_tmp[*((INT *)c0)] - cols_tmp[*((INT *)c1)]) > 0 ?
		1 : (i < 0 ? -1 : 0);
}
#endif /* SORT_COLUMNS */

#if USE_OMP
static void
init_locks(MAT *mat)
{
    if (phgMaxThreads > 1) {
	if (mat->locks == NULL) {
#pragma omp critical (locks)
	    if (mat->locks == NULL) {
		omp_lock_t *locks;
		INT i;
		mat->nlock = mat->localsize;
		if (mat->nlock < mat->rmap->localsize)
		    mat->nlock = mat->rmap->localsize;
		locks = phgAlloc(mat->nlock * sizeof(*mat->locks));
		for (i = 0; i < mat->nlock; i++) {
		    omp_init_lock(&locks[i]);
		}
		mat->locks = locks;
	    }
	}
	assert(mat->locks != NULL || mat->nlock == 0);
    }
}
#endif	/* USE_OMP */

static inline void
add_columns_hashed(INT *ncols1, INT **cols1, FLOAT **data1, INT *size1,
		   INT ncols2, INT *cols2, FLOAT *data2, BOOLEAN repl, FLOAT d,
		   INT hash_size, INT hash_trunk)
/* merge columns using a hash table:
 *	{*ncols1, (*cols1)[], (*data1)[]} += d * {ncols2, cols2, data2}
 *
 * Note:
 * 0. If repl==TRUE then later values overwrite previous ones (PHG_REPLACE),
 *    otherwise they are added up (PHG_ADD).
 * 1. If ncols1 == NULL, then free the hash table and return.
 * 2. If *ncols1 > 0, then the hash table is already in use and matches
 *    {*ncols1, (*cols1)[], (*data1)[]}.
 * 3. cols2[] may contain duplicate entries and (*cols1)[] must contain
 *    distinct entries.
 * 4. If size1!=NULL, then *size1 is the current allocated size of (*cols1)[]
 *    and (*data1)[] and will be updated accordingly, otherwise no realloc of
 *    (*cols1)[] and (*data1)[] is performed.
 * 5. It is allowed that *cols1==cols2 and *data1==data2, in this case we must
 *    have size1==NULL and *ncols1==0 on input.
 * 6. Each entry h = hash[k], for k in [0:hash_size), in the hash table is
 *    organized as follows:
 * 	h[0] 	      = current number of hashed indices
 * 	h[1]	      = current allocated size of h[]
 *	h[2:(h[0]-1)] = list of hashed indices, they are indices to the
 *		 	array (*cols1)[0:(*ncols1)-1].
 *    Let 'col' be a column index in (*cols1)[0:(*ncols1)-1], and let
 *    h := hash[col % hash_size], then we have
 *		col == (*cols1)[h[j]] for some j, j >= 2 && j < h[1]
 */
{
    static INT **hash = NULL;
    static INT hash_size0 = 0;
#if USE_OMP
# pragma omp threadprivate(hash, hash_size0)
#endif	/* USE_OMP */
    INT j, m, n, k, index, ncols, *cols, *h;
    FLOAT *data;

    if (ncols1 == NULL) {
	/* free hash table */
	if (hash != NULL) {
	    if (hash_size <= 0)
		hash_size = hash_size0;
	    for (j = 0; j < hash_size; j++)
		if (hash[j] != NULL)
		    phgFree(hash[j]);
	    phgFree(hash);
	    hash = NULL;
	    hash_size0 = 0;
	}
	return;
    }

    if (ncols2 == 0 || (!repl && d == 0.0))
	return;

    ncols = *ncols1;
    cols = *cols1;
    data = *data1;

    if (ncols == 0) {
	/* reset the hash table */
	if (hash == NULL) {
	    hash = phgCalloc(hash_size0 = hash_size, sizeof(*hash));
	}
	else {
	    for (j = 0; j < hash_size; j++)
		if (hash[j] != NULL)
		    hash[j][0] = 0;
	}
    }

    assert(hash_size == hash_size0);

    for (n = 0; n < ncols2; n++) {
	index = cols2[n];
	if ((h = hash[k = index % hash_size]) == NULL) {
	    h = hash[k] = phgAlloc(hash_trunk * sizeof(*h));
	    h[0] = 0;
	    h[1] = hash_trunk;
	}
	m = h[0];
	/* check if the column index 'index' is already hashed and saved in
	 * cols[0:ncols-1] */
#if 0	/*--------------------------------------------------------------------*/
	/* use binary search, in this case cols[h[]] is sorted and increasing */
    #if 0	/*----------------------------------------------------*/
	k = phgBinarySearchINT(m, cols, h + 2, index);
    #else	/*----------------------------------------------------*/
	if (m <= 1) {
	    k = ((m == 1 && index <= cols[h[2]]) ? 0 : m);
	}
	else {
	    INT a, b, c, l;
	    if ((k = index) <= cols[h[a = 2]]) {
		k = 0;
	    }
	    else if (k >= (l = cols[h[b = 2 + m - 1]])) {
		k = (k == l ? m - 1 : m);
	    }
	    else {
		while (b > a + 1) {
		    if (k == (l = cols[h[c = (a + b) / 2]])) {
			b = c;
			break;
		    }
		    else if (k < l)
			b = c;
		    else
			a = c;
		}
		k = b - 2;
	    }
	}
    #endif	/*----------------------------------------------------*/
	if (k < m && cols[h[2 + k]] == index) {
	    /* column index already hashed, add or replace the value */
	    if (repl)
		data[h[2 + k]] = d * data2[n];
	    else
		data[h[2 + k]] += d * data2[n];
	    continue;
	}
	/* new column index, insert it into the hash table */
	if (m + 2 >= h[1]) {
	    /* enlarge current hash table */
	    h[1] += hash_trunk;
	    h = hash[index % hash_size] = phgRealloc_(h, h[1] * sizeof(*h),
							 (m + 2) * sizeof(*h));
	}
	for (; m > k; m--)
	    h[2 + m] = h[2 + m - 1];
#else	/*--------------------------------------------------------------------*/
	/* use linear search */
	for (k = 0; k < m; k++) {
	    if (index == cols[h[2 + k]]) {
		/* column index already hashed, add to or replace the value */
		if (repl)
		    data[h[2 + k]] = d * data2[n];
		else
		    data[h[2 + k]] += d * data2[n];
		break;
	    }
	}
	if (k < m)
	    continue;
	/* new column index, append it to the hash table */
	if (m + 2 >= h[1]) {
	    /* enlarge current hash table */
	    h[1] += hash_trunk;
	    h = hash[index % hash_size] = phgRealloc_(h, h[1] * sizeof(*h),
							 (m + 2) * sizeof(*h));
	}
#endif	/*--------------------------------------------------------------------*/
	/* append the new column index to cols[] and data[] */
	if (size1 != NULL && ncols >= *size1) {
	    /* enlarge cols[] and data[] */
	    *size1 += 4096;
	    *cols1 = cols = phgRealloc_(cols, *size1 * sizeof(*cols),
					      ncols * sizeof(*cols));
	    *data1 = data = phgRealloc_(data, *size1 * sizeof(*data),
					      ncols * sizeof(*data));
	}
	cols[ncols] = index;
	data[ncols] = d * data2[n];
	h[m + 2] = ncols++;
	h[0]++;
    }
    *ncols1 = ncols;
}

/* TODO: new algorithm (supposedly faster) for phgMatAddEntries functions:
 *	Simply appending newly added entries to existing ones.
 *	Duplicate columns will be handled in phgMatAssemble row by row using
 *	a hash table (see the phgMatMat code).
 *	Drawback: more memory required before matrix assembly.
 *
 * If the above new algorithm is used, then the macro SORT_COLUMNS and the
 * argument 'sorted' will become obsolete */

static int
#if 0
#warning debugging code!
__attribute__ ((noinline))
#endif
add_columns_to_row(MAT *mat, INT irow, INT ncols, INT *cols, FLOAT *data,
		   BOOLEAN sorted, BOOLEAN replace)
/* note: 'sorted' indicates whether the input array 'cols' is sorted,
 * 	 'replace' indicates whether replace or add the values */
{
    MAT_ROW *row;
    INT index;
    INT i;
#if SORT_COLUMNS
    INT index0, *old_cols;
    FLOAT *old_data;
    INT j, k, *ordering;	/* sorted indices of the array 'cols' */
#else /* SORT_COLUMNS */
    Unused(sorted);
#endif /* SORT_COLUMNS */

    /*for (i = 0; i < ncols; i++)
       if (Fabs(data[i]) < 1e-10)
       data[i] = 0.0; */

    if (ncols <= 0)
	return 0;

#if USE_OMP
    if (phgMaxThreads > 1) {
	if (mat->locks == NULL)
	    init_locks(mat);
#if 0
#define omp_unset_lock(l)
#define omp_set_lock(l)
#endif
//phgInfo(-1, "tid %d, setting lock %"dFMT" (nlock=%"dFMT")\n", phgThreadId, irow, mat->nlock);
	omp_set_lock(&mat->locks[irow]);
    }
#endif	/* USE_OMP */

    row = mat->rows + irow;

#if SORT_COLUMNS
    ordering = phgAlloc(ncols * sizeof(*ordering));
    for (i = 0; i < ncols; i++)
	ordering[i] = i;
    if (!sorted) {
	CheckThread
	cols_tmp = cols;
	qsort(ordering, ncols, sizeof(*ordering), comp_columns);
    }
#endif /* SORT_COLUMNS */

    if (row->alloc == 0) {
	/* faster code for a new row */
	row->alloc = ((ncols + TRUNK_SIZE - 1) / TRUNK_SIZE) * TRUNK_SIZE;
	row->cols = phgAlloc(row->alloc * sizeof(*row->cols));
	row->data = phgAlloc(row->alloc * sizeof(*row->data));
	row->ncols = ncols;
#if SORT_COLUMNS
	for (i = 0; i < ncols; i++) {
	    row->cols[i] = cols[ordering[i]];
	    row->data[i] = data[ordering[i]];
	}
	phgFree(ordering);
#else /* SORT_COLUMNS */
	memcpy(row->cols, cols, ncols * sizeof(*cols));
	memcpy(row->data, data, ncols * sizeof(*data));
#endif /* SORT_COLUMNS */

#if USE_OMP
	if (phgMaxThreads > 1) {
//printf("tid %d, unsetting lock %"dFMT"\n", phgThreadId, irow);
	    omp_unset_lock(&mat->locks[irow]);
	}
#endif	/* USE_OMP */

	return 0;
    }

    /* merge new cols with existing cols */

#if SORT_COLUMNS
#if 1
    /* allocate new arrays */
    old_cols = row->cols;
    old_data = row->data;
    row->alloc = (row->ncols <= ncols) ? ncols : row->ncols;
    row->alloc = ((row->alloc + TRUNK_SIZE - 1) / TRUNK_SIZE) * TRUNK_SIZE;
    row->cols = phgAlloc(row->alloc * sizeof(INT));
    row->data = phgAlloc(row->alloc * sizeof(FLOAT));
    i = j = k = 0;
    /* merge arrays */
    while (i < ncols || j < row->ncols) {
	if (k >= row->alloc) {
	    /* need to enlarge the arrays */
	    row->cols = phgRealloc_(row->cols, 
				(row->alloc + TRUNK_SIZE) * sizeof(*row->cols),
				row->alloc * sizeof(*row->cols));
	    row->data = phgRealloc_(row->data,
				(row->alloc + TRUNK_SIZE) * sizeof(*row->data),
				row->alloc * sizeof(*row->data));
	    row->alloc += TRUNK_SIZE;
	}

	if (j >= row->ncols) {
	    row->cols[k] = cols[index = ordering[i++]];
	    row->data[k++] = data[index];
	}
	else if (i >= ncols) {
	    row->cols[k] = old_cols[j];
	    row->data[k++] = old_data[j++];
	}
	else if ((index = cols[ordering[i]]) == (index0 = old_cols[j])) {
	    row->cols[k] = index;
	    if (replace) {
#if 0 * DEBUG
		static BOOLEAN warned = FALSE;
		if (!warned && Fabs(data[ordering[i]] - old_data[j])
					> 1e-6 * (1. + Fabs(old_data[j]))) {
		    warned = TRUE;
		    phgWarning("inconsistent matrix entries: %lf <-> %lf\n",
			       (double)data[ordering[i]], (double)old_data[j]);
		}
#endif	/* DEBUG */
		row->data[k++] = data[ordering[i++]];
		j++;
	    }
	    else {
		row->data[k++] = data[ordering[i++]] + old_data[j++];
	    }
	}
	else if (index > index0) {
	    row->cols[k] = index0;
	    row->data[k++] = old_data[j++];
	}
	else {			/* index < index0 */
	    row->cols[k] = index;
	    row->data[k++] = data[ordering[i++]];
	}
    }
    phgFree(old_data);
    phgFree(old_cols);
    phgFree(ordering);
    row->ncols = k;
#else /* 0|1 */
    /* first, count number of distinct cols */
    k = row->ncols + ncols;
    index = cols[ordering[i = 0]];
    index0 = row->cols[j = 0];
    while (i < ncols && j < row->ncols) {
	if (index == index0) {
	    k--;
	    if (++i >= ncols || ++j >= row->ncols)
		break;
	    index = cols[ordering[i]];
	    index0 = row->cols[j];
	}
	else if (index < index0) {
	    if (++i >= ncols)
		break;
	    index = cols[ordering[i]];
	}
	else {
	    if (++j >= row->ncols)
		break;
	    index0 = row->cols[j];
	}
    }

    if (k > row->alloc) {
	/* need to enlarge the arrays */
	size_t size = ((k + TRUNK_SIZE - 1) / TRUNK_SIZE) * TRUNK_SIZE;
	row->cols = phgRealloc_(row->cols, size * sizeof(*row->cols),
				row->alloc * sizeof(*row->cols));
	row->data = phgRealloc_(row->data, size * sizeof(*row->data),
				row->alloc * sizeof(*row->data));
	row->alloc = size;
    }

    index = cols[ordering[i = ncols - 1]];
    index0 = row->cols[j = row->ncols - 1];
    row->ncols = k;
    while (--k >= 0) {
	assert(i >= 0 || j >= 0);
	if (i < 0) {
	    if (k != j) {
		row->cols[k] = row->cols[j];
		row->data[k] = row->data[j--];
	    }
	    else {
		j--;
	    }
	    if (j >= 0)
		index0 = row->cols[j];
	}
	else if (j < 0) {
	    row->cols[k] = cols[index = ordering[i--]];
	    row->data[k] = data[index];
	    if (i >= 0)
		index = cols[ordering[i]];
	}
	else if (index == index0) {
	    row->cols[k] = index;
	    if (replace) {
#if 0 * DEBUG
		static BOOLEAN warned = FALSE;
		if (!warned && Fabs(data[ordering[i]] - row->data[j])
					> 1e-6 * (1. + Fabs(row->data[j]))) {
		    warned = TRUE;
		    phgWarning("inconsistent matrix entries: %lf <-> %lf\n",
			       (double)data[ordering[i]], (double)row->data[j]);
		}
#endif	/* DEBUG */
		row->data[k] = data[ordering[i--]];
		j--;
	    }
	    else {
		row->data[k] = data[ordering[i--]] + row->data[j--];
	    }
	    if (j >= 0)
		index0 = row->cols[j];
	    if (i >= 0)
		index = cols[ordering[i]];
	}
	else if (index > index0) {
	    row->cols[k] = index;
	    row->data[k] = data[ordering[i--]];
	    if (i >= 0)
		index = cols[ordering[i]];
	}
	else {
	    if (k != j) {
		row->cols[k] = index0;
		row->data[k] = row->data[j--];;
	    }
	    else {
		j--;
	    }
	    if (j >= 0)
		index0 = row->cols[j];
	}
    }
#endif /* 0|1 */
#else /* SORT_COLUMNS */
    while (ncols-- > 0) {
	index = *(cols++);
	for (i = 0; i < row->ncols; i++)
	    if (row->cols[i] == index)
		break;
	if (i < row->ncols) {
	    /* add to existing column */
	    if (replace) {
#if 0 * DEBUG
		static BOOLEAN warned = FALSE;
		if (!warned && Fabs(*data - row->data[i])
					> 1e-6 * (1. + Fabs(row->data[i]))) {
		    warned = TRUE;
		    phgWarning("inconsistent matrix entries: %lf <-> %lf\n",
			       (double)*data, (double)row->data[i]);
		}
#endif	/* DEBUG */
		row->data[i] = *(data++);
	    }
	    else {
		row->data[i] += *(data++);
	    }
	    continue;
	}
	if (row->ncols >= row->alloc) {
	    row->cols = phgRealloc_(row->cols,
				(row->alloc + TRUNK_SIZE) * sizeof(*row->cols),
				row->alloc * sizeof(*row->cols));
	    row->data = phgRealloc_(row->data,
				(row->alloc + TRUNK_SIZE) * sizeof(*row->data),
				row->alloc * sizeof(*row->data));
	    row->alloc += TRUNK_SIZE;
	}
	row->cols[row->ncols] = index;
	row->data[row->ncols++] = *(data++);
    }
#endif /* SORT_COLUMNS */

#if USE_OMP
    if (phgMaxThreads > 1) {
//printf("tid %d, unsetting lock %"dFMT"\n", phgThreadId, irow);
	omp_unset_lock(&mat->locks[irow]);
    }
#endif	/* USE_OMP */

    return 0;
}

void
phgVecDisassemble(VEC *vec)
{
    if (vec->assembled == FALSE)
	return;

    /* Note: set vec->localsize to -1 to make VecAddEntries initialize
     * vec->O2Gmap, etc., when needed. */
#if USE_MPI
    vec->localsize = -1;
#endif	/* USE_MPI */
    vec->assembled = FALSE;
}

void
phgMatDisassemble(MAT *mat)
{
    INT i;
    MAT_ROW *row_bdry;
#if USE_MPI
    MAT_ROW *row;
    INT j, k, ncols;
#endif	/* USE_MPI */

    if (mat->assembled == FALSE)
	return;

    if (mat->type == PHG_MATRIX_FREE) {
	phgError(1, "in %s (%s:%d): unimplemented feature.\n",
					__func__, __FILE__, __LINE__);
    }

    if (mat->type == PHG_PACKED)
	phgMatUnpack(mat);

    /* restore entries in mat->bdry_rows[] */
    if (mat->bdry_rows != NULL) {
	row_bdry = mat->bdry_rows;
	for (i = 0; i < mat->rmap->nlocal; i++, row_bdry++) {
	    add_columns_to_row(mat, i, row_bdry->ncols, row_bdry->cols,
					row_bdry->data, SORT_COLUMNS, TRUE);
	}
	free_matrix_rows(mat->rmap->nlocal, &mat->bdry_rows);
    }

    /* discard mat->bdry[] */
    if (mat->bdry != NULL) {
	phgFree(mat->bdry);
	mat->bdry = NULL;
    }

    mat->bdry_localsize = 0;

    /* TODO: discard bdry eqns if handle_bdry_eqns == TRUE and
     * mat->bdry_localsize > 0 */
    assert(mat->bdry_eqns == NULL);

#if USE_MPI
    /* convert indices above mat->localsize to global ones */
    for (i = 0, row = mat->rows; i < mat->cmap->nlocal; i++, row++) {
	ncols = row->ncols;
	for (j = 0; j < ncols; j++) {
	    /* convert indices above mat->cmap->localsize to global ones */
	    if ((k = row->cols[j]) >= mat->cmap->localsize)
		row->cols[j] = -1 - mat->O2Gmap[k - mat->cmap->nlocal];
	}
    }

    /* Note: set mat->localsize to -1 to tell convert_row_index() that
     *	     mat->O2Gmap, etc., are uninitialized. */
    mat->localsize = -1;
#endif	/* USE_MPI */

    mat->mode = PHG_UNDEFINED;
    mat->assembled = FALSE;
}

MAT *
phgMatAXPBY_(FLOAT a, MAT *x, FLOAT b, MAT **py, BOOLEAN force_matrix_free)
/* computes y := a * x + b * y.
 *
 * If force_matrix_free == TRUE then the resulting matrix will always be
 * matrix-free */
{
    MAT *y;
    INT i, j, ncols, nbuffer = 0, *cols = NULL;
#if USE_MPI
    INT k, maxcols = 0;
#endif	/* USE_MPI */
    FLOAT *buffer = NULL, *data;
    MAT_ROW *row, *row0;

    assert(a == 0.0 || x != NULL);
    assert(x != NULL || (py != NULL && *py != NULL));

    if (x != NULL && !x->assembled) {
	assert(x->type != PHG_DESTROYED);
	phgMatAssemble(x);
    }

    /* If force_matrix_free, or one of the matrices is matrix-free matrix,
     * then the resulting matrix is a matrix-free matrix. */
    if (force_matrix_free
	||
	(b != 0. && py != NULL && *py != NULL && (*py)->type == PHG_MATRIX_FREE)
	||
	(a != 0. && x != NULL && x->type == PHG_MATRIX_FREE)) {
	if (py == NULL || *py == NULL) {
	    y = phgMapCreateMatrixFreeMat(x->rmap, x->cmap, NULL, NULL);
	    y->mv_a = a;
	    y->mv_b = 0.0;
	    y->mv_x = x;
	    if (x != NULL)
		x->refcount++;
	    if (py != NULL)
		*py = y;
	    return y;
	}
	y = *py;
	if (y->type == PHG_MATRIX_FREE) {
	    if (a == 0.0 || x == NULL) {
		y->mv_b *= b;
		return y;
	    }
	    if (y->mv_a == 0. || y->mv_x == NULL) {
		if (y->mv_x != NULL)
		    phgMatDestroy(&y->mv_x);
		y->mv_b *= b;
		y->mv_a = a;
		y->mv_x = x;
		x->refcount++;
		return y;
	    }
	}
	/* Hacky: use an 1x1 block matrix to handle this case. */
	y = phgMatCreateBlockMatrix(y->cmap->comm, 1, 1, &y, NULL, NULL);
	y->mv_b = b;
	y->mv_a = a;
	y->mv_x = x;
	if (x != NULL)
	    x->refcount++;
	*py = y;
	return y;
    }

    if (py == NULL || *py == NULL) {
	b = 0.0;
	y = phgMapCreateMat(x->rmap, x->cmap);
	if (py != NULL)
	    *py = y;
    }
    else {
	y = *py;
	assert(y->type != PHG_DESTROYED);
	phgMatUnpack(y);
	if (x != NULL) {
	    assert(x->rmap->nglobal == y->rmap->nglobal &&
		   x->rmap->nlocal == y->rmap->nlocal &&
		   x->cmap->nglobal == y->cmap->nglobal &&
		   x->cmap->nlocal == y->cmap->nlocal &&
		   x->cmap->localsize == y->cmap->localsize);
	    assert((x->cmap->L2Vmap == NULL) == (y->cmap->L2Vmap == NULL));
#if 1
	    /* TODO: if maps of x and y don't match, then implement the
	     * function using phMatAddGlobalEntries */
	    if (x->cmap->L2Vmap != NULL)
		assert(!memcmp(x->cmap->L2Vmap, y->cmap->L2Vmap,
				x->cmap->localsize * sizeof(*x->cmap->L2Vmap)));
	    if ((i = x->cmap->localsize - x->cmap->nlocal) > 0) {
		assert(x->cmap->O2Gmap != NULL && y->cmap->O2Gmap != NULL);
		assert(!memcmp(x->cmap->O2Gmap, y->cmap->O2Gmap,
				i * sizeof(*x->cmap->O2Gmap)));
	    }
#endif
	}
    }

    if (x != NULL) {
	phgMatUnpack(x);
    }

    if (a == 0.0 && b == 1.0)
	return y;

    /* TODO: faster code for b == 0.0 */

    if (b != 0.0) {
	/* disassemble y to allow adding new entries. */
	phgMatDisassemble(y);
    }

    for (i = 0, row = y->rows; i < y->cmap->nlocal; i++, row++) {
	/* compute b * y */
	if (b == 0.0) {
	    row->ncols = row->alloc = 0;
	    phgFree(row->cols);
	    row->cols = NULL;
	    phgFree(row->data);
	    row->data = NULL;
	}
	else if (b == -1.0) {
	    ncols = row->ncols;
	    for (j = 0; j < ncols; j++) {
		row->data[j] = -row->data[j];
	    }
	}
	else if (b != 1.0) {
	    ncols = row->ncols;
	    for (j = 0; j < ncols; j++) {
		row->data[j] *= b;
	    }
	}

	/* add a * x to y */
	if (a == 0.0)
	    continue;
	row0 = x->rows + i;
	ncols = row0->ncols;
	if (ncols == 0)
	    continue;
	/* copy a*x data to 'buffer' (pointed to by 'data') */
	data = row0->data;
	if (a != 1.0) {
	    if (nbuffer < ncols) {
		phgFree(buffer);
		buffer = phgAlloc((nbuffer = ncols + 100) * sizeof(*buffer));
	    }
	    if (a == -1.0) {
		for (j = 0; j < ncols; j++)
		    buffer[j] = - *(data++);
	    }
	    else {
		for (j = 0; j < ncols; j++)
		    buffer[j] = a * *(data++);
	    }
	    data = buffer;
	}

#if USE_MPI
	/* convert x cols to y cols */
	if (maxcols < ncols) {
	    phgFree(cols);
	    cols = phgAlloc((maxcols = ncols + 100) * sizeof(*cols));
	}
	for (j = 0; j < ncols; j++) {
	    if ((k = row0->cols[j]) < y->cmap->localsize ) {
		cols[j] = k;
		continue;
	    }
	    /* Note: negative value in cols before assembly means
	     * -1 - (global index), see phgMatAssemble */
	    cols[j] = -x->O2Gmap[k - x->cmap->nlocal] - 1;
	}
	add_columns_to_row(y, i, ncols, cols, data, FALSE, FALSE);
#else	/* USE_MPI */
	add_columns_to_row(y, i, ncols, row0->cols, data, FALSE, FALSE);
#endif	/* USE_MPI */

	/* TODO: eleminate zero entries */
    }

    phgFree(buffer);
    phgFree(cols);
    free_matrix_rows(y->rmap->nlocal, &y->bdry_rows);
    free_bdry_eqns(y);
    phgMatAssemble(y);

    return y;
}

MAT *
phgMatGetOffprocRows(MAT *B, INT n, const INT list[],
		     BOOLEAN list_is_ordered, const INT ordering[])
/* returns a MAT containing off-process rows of B listed in the array list[]
 * (whose size is 'n')
 *
 * Note:
 *   1. The function may return NULL if 'B' is not partitioned or 'n' == 0 on
 *	all processes.
 *   2. If 'list_is_ordered' == TRUE, then 'list' should be strictly increasing
 *   	and 'ordering' is not used.
 *   3. If 'list_is_ordered' == FALSE, then if 'ordering' != NULL, then
 *	list[ordering[]] should be strictly increasing, otherwise a temporary
 *	ordering array will be created.
 */
{
    MAT *O = NULL;

#if USE_MPI
    if (B->rmap->nprocs > 1) {
	const MAT_ROW *row;
	MAP *rmap;
	rmap = phgMapCreateSimpleMap(B->rmap->comm, n, PHG_DECIDE);
	if (rmap->nglobal > 0) {
	    MAP *cmap;
	    COMM_INFO *cinfo;
	    INT i, j, *sind, *offset, *ordering1 = NULL;
	    int ii;
	    /* Note: use a serial cmap for O to avoid building O2Gmap and
	     * converting indices */
	    cmap = phgMapCreateSimpleMap(MPI_COMM_SELF, B->cmap->nglobal,
							B->cmap->nglobal);
	    O = phgMapCreateMat(rmap, cmap);
	    /* collect off-process rows of B into O */
	    if (list_is_ordered == TRUE) {
		ordering = NULL;
	    }
	    else if (ordering == NULL && n > 0) {
		_phg_build_ordering(n, list, &ordering1);
		ordering = ordering1;
	    }
	    cinfo = phgMapCreateCommInfo(
				B->rmap->comm, B->rmap->nprocs, B->rmap->rank,
				B->rmap->nlocal, B->rmap->nlocal + n,
				list, ordering, B->rmap->partition);
	    offset = phgAlloc(2 * rmap->nprocs * sizeof(*offset));
	    for (ii = 0; ii < rmap->nprocs; ii++)
		offset[rmap->nprocs + ii] = cinfo->scnts[ii];
	    MPI_Scan(offset + rmap->nprocs, offset, rmap->nprocs,
		     PHG_MPI_INT, MPI_SUM, rmap->comm);
	    for (ii = 0; ii < rmap->nprocs; ii++)
		offset[ii] -= cinfo->scnts[ii];
	    sind = cinfo->sind;
	    for (ii = 0; ii < rmap->nprocs; ii++) {
		for (i = 0; i < cinfo->scnts[ii]; i++, sind++) {
		    row = phgMatGetRow(B, *sind);
		    if (row->ncols <= 0)
			continue;
		    j = rmap->partition[ii] + offset[ii] + i;
		    phgMatAddGlobalEntries(O, 1, &j,
				row->ncols, row->cols, row->data);
		}
	    }
	    phgFree(ordering1);
	    phgFree(offset);
	    phgMapDestroyCommInfo(&cinfo);
	    phgMapDestroy(&cmap);
	    phgMatAssemble(O);
	}
	phgMapDestroy(&rmap);
    }
#endif	/* USE_MPI */

    return O;
}

MAT *
phgMatTranspose(MAT *A)
{
    INT i, gi;
    const MAT_ROW *row;
    MAT *At = phgMapCreateMat(A->cmap, A->rmap);

    if (!A->assembled)
	phgMatAssemble(A);

    gi = A->rmap->partition[A->rmap->rank];
    for (i = 0; i < A->rmap->nlocal; i++, gi++) {
	row = phgMatGetRow(A, i);
	phgMatAddGlobalEntries(At, row->ncols, row->cols, 1, &gi, row->data);
    }

    phgMatAssemble(At);

    return At;
}

/* The following macros controls whether to use HYPRE or PETSc's MatMat
 * functions in phgMatMat. May change them when making the library as, e.g.:
 * make -s USER_CFLAGS="-DUSE_HYPRE_ParMatmul=0 -DUSE_PETSC_MatMatMult=0" lib
 */

/* whether use HYPRE hypre_ParMatmul */
#ifndef USE_HYPRE_ParMatmul
# define USE_HYPRE_ParMatmul	(TRUE && USE_HYPRE && !USE_PETSC)
#endif	/* !undefined(USE_HYPRE_ParMatmul) */

/* whether use PETSc MatMat */
#ifndef USE_PETSC_MatMatMult
# define USE_PETSC_MatMatMult	(TRUE && USE_PETSC)
#endif	/* !undefined(USE_PETSC_MatMatMult) */

#if USE_HYPRE_ParMatmul
#include <HYPRE.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <_hypre_parcsr_mv.h>
#if HYPRE_VERSION_MAJOR < 2 || \
    (HYPRE_VERSION_MAJOR == 2 && HYPRE_VERSION_MINOR < 7)
typedef int HYPRE_Int;
#endif
static HYPRE_IJMatrix
setup_hypre_mat(MAT *A)
{
    HYPRE_IJMatrix a;
    const MAT_ROW *row;
    INT i, j;
    HYPRE_Int gi, ncols, size = 0, *cols = NULL;
    double *data = NULL;

    HYPRE_IJMatrixCreate(A->rmap->comm,
			 A->rmap->partition[A->rmap->rank],
			 A->rmap->partition[A->rmap->rank + 1] - 1,
			 A->cmap->partition[A->cmap->rank],
			 A->cmap->partition[A->cmap->rank + 1] - 1,
			 &a);
    HYPRE_IJMatrixSetObjectType(a, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(a);
    
    gi = A->rmap->partition[A->rmap->rank];
    for (i = 0; i < A->rmap->nlocal; i++, gi++) {
	row = phgMatGetRow(A, i);
	if (row->ncols > 0) {
	    if (size < row->ncols) {
		size = row->ncols;
		phgFree(cols);
		cols = phgAlloc(size * sizeof(*cols));
		phgFree(data);
		data = phgAlloc(size * sizeof(*data));
	    }
	    for (j = 0; j < row->ncols; j++) {
		cols[j] = row->cols[j];
		data[j] = row->data[j];
	    }
	    ncols = row->ncols;
	    HYPRE_IJMatrixSetValues(a, 1, &ncols, &gi, cols, data);
	}
    }
    phgFree(cols);
    phgFree(data);
    HYPRE_IJMatrixAssemble(a);

    /* HYPRE_ParCSRMatrix is in fact typedef'ed to 'hypre_ParCSRMatrix *' */
    return a;
}
#endif	/* USE_HYPRE_ParMatmul */

MAT *
phgMatMat(MAT_OP transa, MAT_OP transb, FLOAT alpha, MAT *A, MAT *B,
	  FLOAT beta, MAT **Cptr)
/* BLAS gemm like function computing C := alpha * op(A) * op(B) + beta * C */
{
    MAT *C, *T = NULL;

    assert((transa == MAT_OP_N || transa == MAT_OP_T) &&
	   (transb == MAT_OP_N || transb == MAT_OP_T));

    if (alpha == 0.0)
	return phgMatAXPBY(alpha, NULL, beta, Cptr);

    assert(A != NULL && B != NULL);
    assert(phgMapCompare(transa == MAT_OP_N ? A->cmap : A->rmap,
			 transb == MAT_OP_N ? B->rmap : B->cmap));

    if (!A->assembled)
	phgMatAssemble(A);
    if (!B->assembled)
	phgMatAssemble(B);

    if (beta == 0.0 && Cptr != NULL && *Cptr != NULL)
	phgMatDestroy(Cptr);

    if (Cptr == NULL || *Cptr == NULL) {
	C = phgMapCreateMat(transa == MAT_OP_N ? A->rmap : A->cmap,
			    transb == MAT_OP_N ? B->cmap : B->rmap);
	if (Cptr != NULL)
	    *Cptr = C;
    }
    else {
	C = *Cptr;
	assert(phgMapCompare(C->rmap, transa == MAT_OP_N ? A->rmap : A->cmap));
	assert(phgMapCompare(C->cmap, transb == MAT_OP_N ? B->cmap : B->rmap));
    }

    /* C := beta * C */
    if (beta != 1.0 && beta != 0.0) {
	phgMatAXPBY(0.0, NULL, beta, &C);
	beta = (beta != 0.0 ? 1.0 : 0.0);
    }

    if (C->diag != NULL) {
	/* invalidate C->diag */
	phgFree(C->diag);
	C->diag = NULL;
    }

    if (C->diag1 != NULL) {
	/* invalidate C->diag1 */
	phgFree(C->diag1);
	C->diag1 = NULL;
    }

#if USE_HYPRE_ParMatmul
/* Note: HYPRE_IJMatrix (if it exists) should be destroyed instead of
 * hypre_ParCSRMatrix, to avoid memory leaks. */
#define DestroyPar(a)
#define DestroyIJ(a)
    if (TRUE) {		/* use HYPRE hypre_ParMatmul */
	INT i, *cols = NULL, size = 0;
	FLOAT *data = NULL;
	HYPRE_Int k, start, end;
	HYPRE_IJMatrix aij = NULL, bij = NULL;
	hypre_ParCSRMatrix *a, *b, *c = NULL, *t = NULL;

	if (phgVerbosity > 0)
	    phgPrintf("Note: using HYPRE ParMatmul for matrix multiply.\n");
	aij = setup_hypre_mat(A);
	HYPRE_IJMatrixGetObject(aij, (void **)(void *)&a);
	if (A == B) {
	    bij = NULL;
	    b = a;
	}
	else {
	    bij = setup_hypre_mat(B);
	    HYPRE_IJMatrixGetObject(bij, (void **)(void *)&b);
	}
	if (transa == MAT_OP_N && transb == MAT_OP_N) {
	    c = hypre_ParMatmul(a, b);
	    DestroyIJ(bij);
	    DestroyIJ(aij);
	}
	else if (transa == MAT_OP_N && transb == MAT_OP_T) {
	    hypre_ParCSRMatrixTranspose(b, &t, 1);
	    DestroyIJ(bij);
	    c = hypre_ParMatmul(a, t);
	    DestroyPar(t);
	    DestroyIJ(aij);
	}
	else if (transa == MAT_OP_T && transb == MAT_OP_N) {
	    hypre_ParCSRMatrixTranspose(a, &t, 1);
	    c = hypre_ParMatmul(t, b);
	    DestroyPar(t);
	    DestroyIJ(bij);
	    DestroyIJ(aij);
	}
	else {
	    c = hypre_ParMatmul(b, a);
	    DestroyIJ(bij);
	    DestroyIJ(aij);
	    hypre_ParCSRMatrixTranspose(c, &t, 1);
	    hypre_ParCSRMatrixDestroy(c);
	    c = t;
	    t = NULL;
	}
#undef DestroyPar
#undef DestroyIJ
	phgMatDisassemble(C);
	start = hypre_ParCSRMatrixFirstRowIndex(c);
	end = hypre_ParCSRMatrixLastRowIndex(c) + 1;
	for (i = start; i < end; i++) {
	    HYPRE_Int ncols;
	    HYPRE_Int *cols1;
	    double *data1;
	    k = hypre_ParCSRMatrixGetRow(c, i, &ncols, &cols1, &data1);
	    if (k != 0)
		phgError(1, "%s: got error %d in hypre_ParCSRMatrixGetRow "
			    "(row %d).\n", __func__, k, i);
	    if (ncols > size) {
		size = ncols;
		phgFree(cols);
		cols = phgAlloc(size * sizeof(*cols));
		phgFree(data);
		data = phgAlloc(size * sizeof(*data));
	    }
	    for (k = 0; k < ncols; k++) {
		cols[k] = cols1[k];
		data[k] = alpha * data1[k];
	    }
	    hypre_ParCSRMatrixRestoreRow(c, i, &ncols, &cols1, &data1);
	    phgMatAddGlobalEntries(C, 1, &i, ncols, cols, data);
	}
	hypre_ParCSRMatrixDestroy(c);
/*----------------------------------------------------------------------------
 * Note: as a feature in HYPRE <= 2.8.0b, in the hypre_ParMatmul function of
 * hypre <= 2.8.0b, the product matrix shares the partitionings with the
 * multiplicands, thus the multiplicands should not be destroyed before the
 * product. Thus the macros DestroyPar and DestroyIJ are defined as empty
 * macros above and redefined below to delay destruction of the multiplicands
 * to after the destruction of the product.
 *
 * If, in future versions of HYPRE, the multiplicands can be safely destroyed
 * before the product, then move the '#define DestroyXXX(a) ...' macros below
 * to right after the '#if USE_HYPRE' above, and delete the remaining lines in
 * the following section. This would potentially reduce total memory usage. */
#define DestroyPar(a) {if (a != NULL) {hypre_ParCSRMatrixDestroy(a); a = NULL;}}
#define DestroyIJ(a) {if (a != NULL) {HYPRE_IJMatrixDestroy(a); a = NULL;}}
DestroyPar(t);
DestroyIJ(bij);
DestroyIJ(aij);
#undef DestroyPar
#undef DestroyIJ
/*---------------------------------------------------------------------------*/
	phgFree(cols);
	phgFree(data);
	phgMatAssemble(C);
	return C;
    }
#endif	/* USE_HYPRE_ParMatmul */

#if USE_PETSC_MatMatMult
    if (TRUE) {		/* use PETSc MatMatMult */
	INT i, k;
	INT *cols = NULL, size = 0;
	FLOAT *data = NULL;
	Mat a, b, c;
	PetscInt start, end, ncols1;
	if (phgVerbosity > 0)
	    phgPrintf("Note: using PETSc MatMatMult for matrix multiply.\n");
	a = phgPetscCreateMatAIJ(A);
	b = (A != B ? phgPetscCreateMatAIJ(B) : a);
	if (transa == MAT_OP_N && transb == MAT_OP_N) {
	    MatMatMult(a, b, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &c);
	    if (b != a)
		phgPetscMatDestroy(&b);
	    phgPetscMatDestroy(&a);
	}
	else if (transa == MAT_OP_N && transb == MAT_OP_T) {
	    Mat bt;
	    MatTranspose(b, MAT_INITIAL_MATRIX, &bt);
	    if (b != a)
		phgPetscMatDestroy(&b);
	    MatMatMult(a, bt, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &c);
	    phgPetscMatDestroy(&bt);
	    phgPetscMatDestroy(&a);
	}
	else if (transa == MAT_OP_T && transb == MAT_OP_N) {
	    Mat at;
	    MatTranspose(a, MAT_INITIAL_MATRIX, &at);
	    MatMatMult(at, b, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &c);
	    phgPetscMatDestroy(&at);
	    if (b != a)
		phgPetscMatDestroy(&b);
	    phgPetscMatDestroy(&a);
	}
	else {
	    /* PETSc doesn't provide MatTransposeMatTransposeMult? */
	    MatMatMult(b, a, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &c);
	    if (b != a)
		phgPetscMatDestroy(&b);
	    phgPetscMatDestroy(&a);
	    a = c;
	    MatTranspose(a, MAT_INITIAL_MATRIX, &c);
	    phgPetscMatDestroy(&a);
	}
	/* add results to C */
	phgMatDisassemble(C);
	MatGetOwnershipRange(c, &start, &end);
	for (i = start; i < end; i++) {
	    PetscInt ncols;
	    const PetscInt *cols1;
	    const PetscScalar *data1;
	    MatGetRow(c, i, &ncols1, &cols1, &data1);
	    ncols = ncols1;
	    if (ncols > size) {
		size = ncols;
		phgFree(cols);
		cols = phgAlloc(size * sizeof(*cols));
		phgFree(data);
		data = phgAlloc(size * sizeof(*data));
	    }
	    for (k = 0; k < ncols; k++) {
		cols[k] = cols1[k];
		data[k] = alpha * data1[k];
	    }
	    MatRestoreRow(c, i, &ncols1, &cols1, &data1);
	    phgMatAddGlobalEntries(C, 1, &i, ncols, cols, data);
	}
	phgPetscMatDestroy(&c);
	phgFree(cols);
	phgFree(data);
	phgMatAssemble(C);
	return C;
    }
#endif	/* USE_PETSC_MatMatMult */

    if (transa == MAT_OP_N && transb == MAT_OP_N) {
	INT *cols = NULL, size = 0;
	FLOAT *data = NULL;
	INT i, n0, n1;
	/* The following arrays store CSR data for B for faster access */
	size_t *Bind = NULL;
	INT *Bcols = NULL;
	FLOAT *Bdata = NULL;

	if (phgMaxThreads > 1)
	    phgWarning("multi-threaded matrix multiply code is broken!\n");

	if (phgVerbosity > 0)
	    phgPrintf("Note: using PHG's own matrix multiply function.\n");
	phgMatDisassemble(C);
#if USE_MPI
	MAT *O;		/* MAT object for collecting off-proc rows of B */
	INT O2Gsize = A->localsize - A->cmap->nlocal;
	/* fetch off-process rows of B matching off-process cols of A.
	 * 
	 * Note: a matrix O consisting of off-process rows of B is
	 * constructed. This method is sub-optimal because it involves
	 * unnecessary costs of converting between local and global 
	 * indices, and phgMatAddEntries. This may be optimized later. */
	if (A->rmap->nprocs > 1 || B->rmap->nprocs > 1) {
	    if (O2Gsize > 0 && A->O2Gmap == NULL) {
		/* build O2Gmap for A */
		phgError(1, "(%s:%d) unimplemented.\n", __FILE__, __LINE__);
	    }
	    if (O2Gsize < 0)
		O2Gsize = 0;
	    O = phgMatGetOffprocRows(B, O2Gsize, A->O2Gmap, FALSE, A->ordering);
	}
	else {
	    O = NULL;
	}
#endif	/* USE_MPI */
#if 1
	/* retrieve CSR data of B in Bind, Bcols, and Bdata for faster access */
	Bind  = phgAlloc((B->rmap->nlocal + 1) * sizeof(*Bind));
	Bcols = phgAlloc((B->nnz_d + B->nnz_o) * sizeof(*Bcols));
	Bdata = phgAlloc((B->nnz_d + B->nnz_o) * sizeof(*Bdata));
	Bind[0] = 0;
	for (i = 0; i < B->rmap->nlocal; i++) {
	    const MAT_ROW *row;
	    row = phgMatGetRow_(B, i, NULL);
	    Bind[i + 1] = Bind[i] + row->ncols;
	    if (row->ncols == 0)
		continue;
	    memcpy(Bcols + Bind[i], row->cols, row->ncols * sizeof(*Bcols));
	    memcpy(Bdata + Bind[i], row->data, row->ncols * sizeof(*Bdata));
	}
#endif
	/* loop on rows of A */
	n0 = B->rmap->partition[B->rmap->rank];
	n1 = B->rmap->partition[B->rmap->rank + 1];
	for (i = 0; i < A->rmap->nlocal; i++) {
	    INT j, k, ncols, ncols1;
	    FLOAT d;
	    MAT_ROW row_buffer1 = {NULL, NULL, 0, 0}, row_buffer2;
	    const MAT_ROW *row1 = NULL, *row2 = NULL;

	    row1 = phgMatGetRow_(A, i, &row_buffer1);
	    if ((ncols1 = row1->ncols) == 0)
		continue;
	    ncols = 0;
	    /* loop on rows of B matching the current row of A */
	    for (j = 0; j < ncols1; j++) {
		d = alpha * row1->data[j];
		k = row1->cols[j];
		if (k >= n0 && k < n1) {
		    if (Bind != NULL) {
			row_buffer2.cols = Bcols + Bind[k -= n0];
			row_buffer2.data = Bdata + Bind[k];
			row_buffer2.ncols = (INT)(Bind[k + 1] - Bind[k]);
			row2 = &row_buffer2;
		    }
		    else {
			row2 = phgMatGetRow(B, k - n0);
		    }
		}
#if USE_MPI
		else {
		    k = phgBinarySearchINT(O2Gsize, A->O2Gmap, A->ordering, k);
		    assert(k >= 0 && k < O2Gsize &&
			   row1->cols[j] == A->O2Gmap[A->ordering[k]]);
		    row2 = phgMatGetRow(O, k);
		}
#endif	/* USE_MPI */
#ifndef HASH_SIZE
# define HASH_SIZE	257
#endif
#ifndef HASH_TRUNK
# define HASH_TRUNK	32	/* must > 2 */
#endif
		add_columns_hashed(&ncols, &cols, &data, &size,
				   row2->ncols, row2->cols, row2->data,
				   FALSE, d, HASH_SIZE, HASH_TRUNK);
	    }
	    k = A->rmap->partition[A->rmap->rank] + i;
	    phgMatAddGlobalEntries(C, 1, &k, ncols, cols, data);
	    phgFree(row_buffer1.cols);
	    row_buffer1.cols = NULL;
	    phgFree(row_buffer1.data);
	    row_buffer1.data = NULL;
	}
#if USE_MPI
	phgMatDestroy(&O);
#endif	/* USE_MPI */
	phgFree(Bind);
	phgFree(Bcols);
	phgFree(Bdata);
	/* free hash table in add_columns_hashed */
	add_columns_hashed(NULL,NULL,NULL,NULL, 0,NULL,NULL, 0,0, FALSE, 0.);
	phgFree(data);
	data = NULL;
	phgFree(cols);
	cols = NULL;
	size = 0;
	phgMatAssemble(C);
    }
    else if (transa == MAT_OP_N && transb == MAT_OP_T) {
	T = phgMatTranspose(B);
	phgMatMat(MAT_OP_N, MAT_OP_N, alpha, A, T, beta, &C);
	phgMatDestroy(&T);
    }
    else if (transa == MAT_OP_T && transb == MAT_OP_N) {
	T = phgMatTranspose(A);
	phgMatMat(MAT_OP_N, MAT_OP_N, alpha, T, B, beta, &C);
	phgMatDestroy(&T);
    }
    else {
	assert(transa == MAT_OP_T && transb == MAT_OP_T);
	phgMatMat(MAT_OP_N, MAT_OP_N, alpha, B, A, beta, &C);
	T = C;
	C = phgMatTranspose(T);
	phgMatDestroy(&T);
	if (Cptr != NULL)
	    *Cptr = C;
    }

    return C;
}

static void
build_bdry(MAT *mat)
/* builds mat->bdry and computes mat->bdry_localsize */
{
    INT k;
#if USE_MPI
    int rank, *scnts, *sdsps, *rcnts, *rdsps;
    INT i, j, n0, ssize, rsize, *sbuf, *rbuf, *rind;
#endif	/* USE_MPI */

    if (mat->cmap->bdry_nglobal == 0)
	return;

    k = mat->bdry_localsize = mat->cmap->bdry_localsize;
    phgFree(mat->bdry);
    mat->bdry = NULL;
    if (k > 0) {
	mat->bdry = phgAlloc(k * sizeof(*mat->bdry));
	memcpy(mat->bdry, mat->cmap->bdry, k * sizeof(*mat->bdry));
    }

#if USE_MPI
    rsize = mat->localsize - mat->cmap->localsize;

    /* count off-proc columns not covered by mat->cmap */
    scnts = phgAlloc(4 * mat->cmap->nprocs * sizeof(*scnts));
    sdsps = scnts + mat->cmap->nprocs;
    rcnts = sdsps + mat->cmap->nprocs;
    rdsps = rcnts + mat->cmap->nprocs;
    rbuf = phgAlloc(2 * rsize * sizeof(*rbuf));
    rind = rbuf + rsize;
    rank = 0;
    bzero(rcnts, mat->cmap->nprocs * sizeof(*rcnts));
    /* count and collect entries in the range
     *	[mat->cmap->locslsize, mat->localsize) */
    for (j = 0, i = 0; i < mat->localsize - mat->cmap->nlocal; i++) {
	if (mat->ordering[i] < mat->cmap->localsize - mat->cmap->nlocal)
	    continue;
	rind[j] = mat->ordering[i];
	rbuf[j++] = k = mat->O2Gmap[mat->ordering[i]];
	while (rank < mat->cmap->nprocs &&
			mat->cmap->partition[rank + 1] <= k)
	    rank++;
	assert(rank < mat->cmap->nprocs && rank != mat->cmap->rank);
	assert(rcnts[rank] < INT_MAX);
	rcnts[rank]++;
    }
    assert(j == rsize);
    MPI_Alltoall(rcnts, 1, MPI_INT, scnts, 1, MPI_INT, mat->cmap->comm);
    i = j = 0;
    for (rank = 0; rank < mat->cmap->nprocs; rank++) {
	rdsps[rank] = i;
	sdsps[rank] = j;
	i += rcnts[rank];	assert(i <= INT_MAX);
	j += scnts[rank];	assert(j <= INT_MAX);
    }
    sbuf = phgAlloc((ssize = j) * sizeof(*sbuf));
    MPI_Alltoallv(rbuf, rcnts, rdsps, PHG_MPI_INT,
		  sbuf, scnts, sdsps, PHG_MPI_INT, mat->cmap->comm);
    n0 = mat->cmap->partition[mat->cmap->rank];
    if (ssize > 0 && mat->cmap->bdry_nlocal > 0) {
	for (i = 0; i < ssize; i++) {
	    k = sbuf[i] - n0;
	    assert(k >= 0 && k < mat->cmap->nlocal);
	    if (bsearch(&k, mat->cmap->bdry, mat->cmap->bdry_nlocal,
			sizeof(k), phgCompINT) == NULL)
		sbuf[i] = -1;	/* non boundary entry */
	}
    }
    else {
	for (i = 0; i < ssize; i++) {
	    sbuf[i] = -1;	/* all are non boundary entries */
	}
    }
    MPI_Alltoallv(sbuf, scnts, sdsps, PHG_MPI_INT,
		  rbuf, rcnts, rdsps, PHG_MPI_INT, mat->cmap->comm);
    for (i = 0, k = 0; i < rsize; i++) {
	if (rbuf[i] >= 0)
	    k++;
    }
    if (k > 0) {
	k += mat->bdry_localsize;
	mat->bdry = phgRealloc_(mat->bdry, k * sizeof(*mat->bdry),
				mat->bdry_localsize * sizeof(*mat->bdry));
	for (i = 0, k = 0; i < rsize; i++) {
	    if (rbuf[i] < 0)
		continue;
	    mat->bdry[mat->bdry_localsize + k] = mat->cmap->nlocal + rind[i];
	    k++;
	}
	qsort(mat->bdry + mat->bdry_localsize, k, sizeof(k), phgCompINT);
	mat->bdry_localsize += k;
    }
    phgFree(scnts);
    phgFree(rbuf);
    phgFree(sbuf);
#endif	/* USE_MPI */
}

static BOOLEAN
collect_local_bdry_equations(MAT *mat, INT i, int map0[],
			     INT **map1, int *n, int *alloc)
/* Recursively collects coupled local rows starting with the i-th row to form
 * a small closed linear system of local bdry equations.
 *
 * Note: since this function is called with increasing value of i, some rows
 * collected with indices < i may have already been processed (e.g., if the
 * small system is a lower triangular one). Thus to get correct result the
 * function solve_local_bdry_equations() must scan the local bdry equations in
 * the same order. */
{
    static int level = 0;
#if USE_OMP
# pragma omp threadprivate(level)
#endif	/* USE_OMP */
    INT j, k;
    MAT_ROW *row = mat->rows + i;

    if (*n >= *alloc) {
	*map1 = phgRealloc_(*map1, (*alloc + 128) * sizeof(**map1),
				   (*alloc) * sizeof(**map1));
	*alloc += 128;
    }
    map0[i] = *n;
    (*map1)[(*n)++] = i;

    for (j = 0; j < row->ncols; j++) {
	if (map0[k = row->cols[j]] >= 0)
	    continue;
	if (map0[k] == -2) {	/* non boundary entry */
	    static BOOLEAN warned = FALSE;
	    if (!warned) {
		phgWarning("(%s:%d) bdry equations involve non bdry entries "
			   "(level=%d, i=%"dFMT", ncols=%"dFMT", "
			   "j=%"dFMT", k=%"dFMT").\n",
	 		   __FILE__, __LINE__, level, i, row->ncols, j, k);
		warned = TRUE;
	    }
	    /* reset map0[] */
	    for (j = 0; j < *n; j++)
		map0[(*map1)[j]] = -1;
	    *n = 0;
	    return FALSE;
	}
	level++;
	if (!collect_local_bdry_equations(mat, k, map0, map1, n, alloc)) {
	    level--;
	    return FALSE;
	}
	level--;
    }

    return TRUE;
}

static void
save_local_bdry_equations(MAT *mat)
/* do LU factorization and save local bdry equations before matrix assembly */
{
    int *map0;
    INT *map1 = NULL;
    int ii, n = 0, alloc = 0;
    MAT_ROW *row, *row1;
    INT i, j;
    double t0 = phgGetTime(NULL);
    int count = 0, nmax = 0;
    BDRY_EQNS *be;

    if (!mat->handle_bdry_eqns)
	return;

    assert(!mat->assembled);

    if (mat->type == PHG_MATRIX_FREE)
	phgError(1, "calling %s with a matrix-free matrix.\n", __func__);

    phgInfo(1, "trying to diagonalize local boundary equations.\n");

    for (ii = 0; ii < mat->rmap->ndof; ii++) {
	DOF *dof = mat->rmap->dofs[ii];
	if (dof->DB_masks == NULL) {
	    if (dof->DB_mask != 0)
		break;
	} else {
	    for (i = 0; i < dof->dim; i++)
		if (dof->DB_masks[i] != 0)
		    break;
	    if (i < dof->dim)
		break;
	}
    }
    if (ii >= mat->rmap->ndof)	/* no Dirichlet boundary DOF */
	return;

    /* Entries in map0[] (TODO: may use -3 to mark processed bdry rows):
     *	-2: non-boundary entries
     *	-1: boundary entries
     */
    map0 = phgAlloc(mat->cmap->localsize * sizeof(*map0));
    for (i = 0; i < mat->cmap->localsize; i++)
	map0[i] = -2;
    for (i = 0; i < mat->cmap->bdry_localsize; i++)
	map0[mat->cmap->bdry[i]] = -1;

    for (i = 0, row = mat->rows; i < mat->rmap->localsize; i++, row++) {
	if (map0[i] == -2)
	    continue;
	/* Skip empty or non negative diagonal rows (for positiveness?).
	 * Note:
	 *   1. depending on phgDofDirichletBC, some bdry rows may be empty
	 *   2. processed rows are skipped since they have row->data[0]==1 */
	if (row->ncols <= 0 ||
	    (row->ncols == 1 && row->cols[0] == i && row->data[0] >= 0.0))
	    continue;
	/* recursively collect connected boundary equations */
	if (!collect_local_bdry_equations(mat, i, map0, &map1, &n, &alloc))
	    continue;
	assert(n > 1 ||
	       (row->ncols == 1 && row->cols[0] == i && row->data[0] < 0.0));
	if (n > 2048) {
	    phgWarning("%s: bdry eqns too large (%d)\n", __func__, n);
	    /* reset map0[] */
	    for (ii = 0; ii < n; ii++)
		map0[map1[ii]] = -1;
	    n = 0;
	    continue;
	}
	if (n > nmax)
	    nmax = n;
	/* save the boundary equations */
	if (count >= mat->bdry_eqns_allocated) {
	    mat->bdry_eqns = phgRealloc_(mat->bdry_eqns,
			(count + 128) * sizeof(*mat->bdry_eqns),
			mat->bdry_eqns_allocated * sizeof(*mat->bdry_eqns));
	    mat->bdry_eqns_allocated = count + 128;
	}
	be = mat->bdry_eqns + (count++);
	be->n = n;
	be->lu = phgCalloc(n * n, sizeof(*be->lu));
	be->pvt = phgAlloc(n * sizeof(*be->pvt));
	be->map = phgAlloc(n * sizeof(*be->map));
	for (ii = 0; ii < n; ii++) {
	    row1 = mat->rows + map1[ii];
	    for (j = 0; j < row1->ncols; j++)
		be->lu[map0[row1->cols[j]] + ii * n] = row1->data[j];
	    be->map[ii] = map1[ii];
	}
	if (!phgSolverDenseLU(n, (void *)be->lu, be->pvt))
	    phgError(1, "%s: singular bdry equations (n=%d).\n", __func__, n);

	for (ii = 0; ii < n; ii++) {
	    row1 = mat->rows + map1[ii];
	    if (row1->ncols != 1) {
		phgFree(row1->cols);
		phgFree(row1->data);
		row1->cols = phgAlloc(sizeof(*row1->cols));
		row1->data = phgAlloc(sizeof(*row1->data));
		row1->ncols = row1->alloc = 1;
	    }
	    row1->cols[0] = map1[ii];
	    row1->data[0] = 1.0;
	}

	/* reset map0[] */
	for (ii = 0; ii < n; ii++)
	    map0[map1[ii]] = -1;
	n = 0;
    }
    mat->bdry_eqns_count = count;
    phgInfo(1, "diagonalize: %d eqns, max size = %d, time = %0.4lf\n",
	    count, nmax, (double)(phgGetTime(NULL) - t0));

    phgFree(map0);
    phgFree(map1);
}

static void
save_bdry_columns(MAT *mat, BOOLEAN diag)
/* moves all boundary columns to bdry_rows[] */
{
    MAT_ROW *row, *row_bdry;
    INT i, j, k, l;
    BYTE *flags;

    if (!mat->handle_bdry_eqns)
	return;

    if (mat->type == PHG_MATRIX_FREE) {
	if (mat->blocks != NULL)
	   phgError(1, "%s unimplemented for block matrix.\n", __func__);
	return;
    }

    flags = phgCalloc(mat->localsize, sizeof(*flags));
    free_matrix_rows(mat->rmap->nlocal, &mat->bdry_rows);
    mat->bdry_rows = phgCalloc(mat->rmap->nlocal, sizeof(*row));

    for (i = 0; i < mat->bdry_localsize; i++)
	flags[mat->bdry[i]] = 1;

    row_bdry = mat->bdry_rows;
    row = mat->rows;
    for (i = 0; i < mat->rmap->nlocal; i++, row++, row_bdry++) {
	/* count # of bdry columns */
	for (k = 0, j = 0; j < row->ncols; j++) {
	    if (flags[row->cols[j]] != 0)
		k++;
	}
	row_bdry->ncols = row->alloc = k;
	row_bdry->cols = phgAlloc(k * sizeof(*row_bdry->cols));
	row_bdry->data = phgAlloc(k * sizeof(*row_bdry->data));
	for (k = 0, l = 0, j = 0; j < row->ncols; j++) {
	    /* keep non boundary entries, and diagonal if diag==TRUE
		in mat->rows */
	    if ((diag && row->cols[j] == i) || flags[row->cols[j]] == 0) {
		if (l != j) {
		    row->cols[l] = row->cols[j];
		    row->data[l] = row->data[j];
		}
		l++;
	    }
	    /* save boundary entries to mat->bdry_rows */
	    if (flags[row->cols[j]] != 0) {
		row_bdry->cols[k] = row->cols[j];
		row_bdry->data[k] = row->data[j];
		k++;
	    }
	}
	row->cols = phgRealloc_(row->cols, l * sizeof(*row->cols),
				row->ncols * sizeof(*row->cols));
	row->data = phgRealloc_(row->data, l * sizeof(*row->data),
				row->ncols * sizeof(*row->data));
	row->ncols = row->alloc = l;
    }
    phgFree(flags);
}

static void
update_nnz(MAT *mat)
{
    INT i;
    MAT_ROW *row;
    size_t nnz_d = 0, nnz_o = 0;
#if USE_MPI
    INT j, n;
#endif	/* USE_MPI */

    if (mat == NULL)
	return;

    if (mat->blocks != NULL) {
	MAT **pmat = mat->blocks->pmat;
	for (i = 0; i < mat->blocks->p * mat->blocks->q; i++) {
	    if (pmat[i] == 0)
		continue;
	    nnz_d += pmat[i]->nnz_d;
	    nnz_o += pmat[i]->nnz_o;
	}
	mat->nnz_d = nnz_d;
	mat->nnz_o = nnz_o;
	return;
    }

    if (mat->rows == NULL) {
	/* matrix-free matrix or nlocal == 0 */
	mat->nnz_d = mat->nnz_o = 0;
	return;
    }

    for (i = 0, row = mat->rows; i < mat->rmap->nlocal; i++, row++) {
#if USE_MPI
	n = row->ncols;
	for (j = 0; j < n; j++) {
	    if (row->cols[j] < mat->cmap->nlocal)
		nnz_d++;
	    else
		nnz_o++;
	}
#else	/* USE_MPI */
	nnz_d += row->ncols;
#endif
    }
    mat->nnz_d = nnz_d;
    mat->nnz_o = nnz_o;
}

static void
assemble_block_matrix(MAT *mat)
{
    int i, j, p, q;
    MAT **pmat, *d, *o;

    assert(mat->type == PHG_MATRIX_FREE);

    if (mat->assembled)
	return; 
    mat->assembled = TRUE;

    p = mat->blocks->p;
    q = mat->blocks->q;
    pmat = mat->blocks->pmat;

    if (!mat->handle_bdry_eqns) {
	for (i = 0; i < p * q; i++)
	    if(pmat[i] != NULL)  phgMatAssemble(pmat[i]);
	return;
    }

    assert(p == q);

    for (i = 0; i < p; i++) {
	d = pmat[i * q + i];
	if (d != NULL) {
	    if (!d->handle_bdry_eqns)
		phgError(1, "diagonal block, handle_bdry_eqns = FALSE\n");
	    phgMatAssemble(d);
	}
	for (j = 0; j < mat->blocks->q; j++) {
	    if (j == i)
		continue;
	    o = pmat[i * q + j];
	    if (o == NULL)
		continue;
	    if (o->handle_bdry_eqns)
		phgError(1, "off-diagonal block, handle_bdry_eqns = TRUE\n");
	    phgMatAssemble(o);
	}
    }

    for (i = 0; i < p * q; i++)
	if ((o = pmat[i]) != NULL)
	    o->handle_bdry_eqns = TRUE;	/* to fool save_bdry_columns */

    for (i = 0; i < p; i++) {
	for (j = 0; j < mat->blocks->q; j++) {
	    if (j == i || (o = pmat[i * q + j]) == NULL)
		continue;
	    save_bdry_columns(o, FALSE);
	    /* reset handle_bdry_eqns to FALSE which also neutralizes future
	     * possible save_bdry_columns() calls with the same matrix */
	    o->handle_bdry_eqns = FALSE;
	    update_nnz(o);
	}
    }

    return;
}

void
phgMatAssemble(MAT *mat)
/* assembles matrix, i.e., move off-process rows to corresponding
 * processes (Alltoallv). */
{
#if USE_MPI
    /* The buffer is organized as follows:
       INT   nrows, nnz_low, nnz_high, pad,
       INT   rows[0], ..., rows[nrows-1],
       INT   ncols[0], ..., ncols[nrows-1],
       INT   cols[0][0], ..., cols[0][ncols[0] - 1],
       INT   cols[1][0], ..., cols[1][ncols[1] - 1],
       (padding to alignment of FLOAT)
       FLOAT data[0][0], ..., data[0][ncols[0] - 1],
       FLOAT data[1][0], ..., data[1][ncols[1] - 1],
       ... ...
       FLOAT data[nrows-1][0], ..., data[nrows-1][ncols[nrows-1] - 1],
       (padding to multiple of trunc_size)

       Note: the above may not work if alignment of FLOAT is not a multiple of
	     alignment of INT, producing BUS error.
    */
    unsigned char *sbuff, *rbuff, *align;
    INT ncols, nrows = 0, *O2Gmap, *ordering;
    MAT_ROW *row;
    FLOAT *pdata = NULL;
    INT *prows = NULL, *pncols = NULL, *pcols = NULL;
    INT i, j, k, r, rindex0, cindex0, O2Gmap_alloc, O2Gmap_size;
    int rank, rank0;
    int *scounts, *rcounts, *sdispls, *rdispls;
    MPI_Datatype type;
    size_t size, ssize, rsize, *nnz;
    /* Note: the maxi. size of buffer is limited by '2 * trunc_size' GB and
     * the overhead is bounded by '(nprocs-1) * trunc_size' */ 
    static size_t align_float = 0, align_int = 0, trunc_size = 128;
#if USE_OMP
# pragma omp threadprivate(align_float, align_int, trunc_size)
#endif	/* USE_OMP */

    if (mat->assembled)
	return;

    /* reset mat->bdry_localsize which was used by convert_row_index() */
    mat->bdry_localsize = 0;

    /* mat->localsize > 0 ==> mat->rows contains off-proc rows */
    if (mat->localsize < 0)
	mat->localsize = 0;

    if (mat->blocks != NULL) {
	assemble_block_matrix(mat);
	update_nnz(mat);
	free_locks(mat);
	return;
    }

    save_local_bdry_equations(mat);

    if (mat->rmap->nprocs <= 1 && mat->cmap->nprocs <= 1) {
	phgFree(mat->O2Gmap);
	mat->O2Gmap = NULL;
	phgFree(mat->ordering);
	mat->ordering = NULL;
	mat->localsize = mat->cmap->localsize;
	build_bdry(mat);
	save_bdry_columns(mat, TRUE);
	mat->assembled = TRUE;
	update_nnz(mat);
	free_locks(mat);
	return;
    }

    if (mat->type == PHG_MATRIX_FREE) {
	/* matrix-free matrix */
	assert(mat->rows == NULL);
	mat->assembled = TRUE;
	free_locks(mat);
	return;
    }

    /* Note: the off-proc rows are saved in mat->rows, mat->O2Gmap
     *	     and mat->ordering, and mat->localsize is total number of rows */

    if (mat->mode == PHG_REPLACE) {
	/* Hacky: simply discard all off-proc rows (this mode is currently only
	 * used by solver-ams.c) */
	row = mat->rows + mat->rmap->nlocal;
	for (i = 0; i < mat->localsize - mat->rmap->nlocal; i++, row++) {
	    phgFree(row->cols);
	    row->cols = NULL;
	    phgFree(row->data);
	    row->data = NULL;
	    row->ncols = row->alloc = 0;
	}
	if (mat->localsize > mat->rmap->localsize) {
	    mat->rows = phgRealloc_(mat->rows,
				mat->rmap->localsize * sizeof(*mat->rows),
				mat->rmap->localsize * sizeof(*mat->rows));
	}
	phgFree(mat->O2Gmap);
	mat->O2Gmap = NULL;
	phgFree(mat->ordering);
	mat->ordering = NULL;
	mat->localsize = mat->cmap->localsize;
	O2Gmap_size = mat->cmap->localsize - mat->cmap->nlocal;
	if (O2Gmap_size > 0) {
	    /* copy cmap->O2Gmap/cmap->ordering */
	    mat->O2Gmap = phgAlloc(O2Gmap_size * sizeof(*mat->O2Gmap));
	    mat->ordering = phgAlloc(O2Gmap_size * sizeof(*mat->ordering));
	    memcpy(mat->O2Gmap, mat->cmap->O2Gmap,
			O2Gmap_size * sizeof(*mat->O2Gmap));
	    memcpy(mat->ordering, mat->cmap->ordering,
			O2Gmap_size * sizeof(*mat->ordering));
	}
	mat->assembled = TRUE;
	update_nnz(mat);
	free_locks(mat);
	return;
    }

    if (align_float == 0) {
	/* compute alignment of FLOAT and INT (TODO: add configure check?) */
	struct {
	    BYTE a;
	    FLOAT b;
	} a;
	struct {
	    BYTE a;
	    INT b;
	} b;
	align_float = sizeof(a) - sizeof(FLOAT);
	align_int = sizeof(b) - sizeof(INT);
	if (trunc_size % align_int != 0 || trunc_size % align_float != 0 ||
	    (align_float % align_int != 0 && align_int % align_float != 0))
	    phgError(1, "%s:%d, unexpected alignment.\n", __FILE__, __LINE__);
    }

    scounts = phgAlloc(4 * mat->rmap->nprocs * sizeof(*scounts));
    rcounts = scounts + 1 * mat->rmap->nprocs;
    sdispls = scounts + 2 * mat->rmap->nprocs;
    rdispls = scounts + 3 * mat->rmap->nprocs;

    /* count # of rows (scounts[]) and nnz (nnz[]) for each process */

    /* reset counts */
    nnz = phgAlloc(mat->rmap->nprocs * sizeof(*nnz));
    for (rank = 0; rank < mat->rmap->nprocs; rank++) {
	scounts[rank] = 0;
	nnz[rank] = 0;
    }

    rank = 0;
    for (i = 0; i < mat->localsize - mat->rmap->nlocal; i++) {
	row = mat->rows + mat->rmap->nlocal + mat->ordering[i];
	if ((ncols = row->ncols) <= 0)
	    continue;
	r = mat->O2Gmap[mat->ordering[i]];
	for (; r >= mat->rmap->partition[rank + 1]; rank++);
	assert(rank < mat->rmap->nprocs);
	if (ncols > 0) {
	    assert(scounts[rank] < INT_MAX);
	    scounts[rank]++;
	    nnz[rank] += (size_t)ncols;
	}
    }

    /* compute in sdispls offsets of data for each process */
    ssize = 0;
    for (rank = 0; rank < mat->rmap->nprocs; rank++) {
	nrows = scounts[rank];

	/* size of header + INT entries */
	size = (4 + nrows * 2 + nnz[rank]) * sizeof(INT);
	/* adjust size to alignment of FLOAT */
	size = ((size + align_float - 1) / align_float) * align_float;
	/* plus size of FLOAT entries */
	size += nnz[rank] * sizeof(FLOAT);
	/* adjust size to multiple of trunc_size */
	size = ((size + trunc_size - 1) / trunc_size) * trunc_size;

	sdispls[rank] = (int)(ssize / trunc_size);
	ssize += size;
	assert(ssize / trunc_size <= INT_MAX);
    }

    /* allocate send buffer */
    sbuff = phgAlloc(ssize);

    /* setup data counts in sbuff */
    for (rank = 0; rank < mat->rmap->nprocs; rank++) {
	size = sdispls[rank];
	if (size == (rank == mat->rmap->nprocs - 1 ?
			(int)(ssize / trunc_size) : sdispls[rank + 1]))
	    continue;
	size *= trunc_size;
	prows = (void *)(sbuff + size);
	prows[0] = scounts[rank];
	prows[1] = (INT)(nnz[rank] % (size_t)(INT_MAX));
	prows[2] = (INT)(nnz[rank] / (size_t)(INT_MAX));
	prows[3] = 0;
    }
    phgFree(nnz);

    rindex0 = mat->rmap->partition[mat->rmap->rank];
    cindex0 = mat->cmap->partition[mat->cmap->rank];

    O2Gmap_size = O2Gmap_alloc = mat->cmap->localsize - mat->cmap->nlocal;
    O2Gmap = phgAlloc(O2Gmap_alloc * sizeof(*O2Gmap));
    ordering = phgAlloc(O2Gmap_alloc * sizeof(*ordering));
    if (O2Gmap_size > 0) {
	memcpy(O2Gmap, mat->cmap->O2Gmap, O2Gmap_size * sizeof(*O2Gmap));
	memcpy(ordering, mat->cmap->ordering, O2Gmap_size * sizeof(*ordering));
    }

    if (ssize != 0 && mat->localsize > mat->rmap->nlocal) {
	/* collect off-process rows to sbuff */
	rank0 = -1;
	rank = 0;
	for (i = 0; i < mat->localsize - mat->rmap->nlocal; i++) {
	    row = mat->rows + mat->rmap->nlocal + mat->ordering[i];
	    r = mat->O2Gmap[mat->ordering[i]];
	    for (; mat->rmap->partition[rank + 1] <= r; rank++);
	    assert(rank < mat->rmap->nprocs);
	    if (rank != rank0) {
		prows = (void *)(sbuff + ((size_t)sdispls[rank]) * trunc_size);
		nrows = prows[0];
		size = ((size_t)prows[1]) +
		       ((size_t)prows[2]) * (size_t)(INT_MAX);

		prows += 4;
		pncols = prows + nrows;
		pcols = pncols + nrows;

		align = (void *)(pcols + size);
		align = sbuff + ((align - sbuff + align_float - 1)
				/ align_float) * align_float;
		pdata = (void *)align;
		rank0 = rank;
	    }
	    if ((ncols = row->ncols) > 0) {
		*(prows++) = r;
		*(pncols++) = ncols;
		for (j = 0; j < ncols; j++) {
		    k = row->cols[j];
		    if (k < 0) {
			pcols[j] = -1 - k; /* restore negated global index */
			continue;
		    }
		    /* convert local index to global index */
		    pcols[j] = (k < mat->cmap->nlocal || O2Gmap == NULL ?
				k + cindex0 : O2Gmap[k - mat->cmap->nlocal]);
		}
		pcols += ncols;
		memcpy(pdata, row->data, sizeof(FLOAT) * ncols);
		pdata += ncols;
	    }
	    if (row->alloc > 0) {
		phgFree(row->cols);
		row->cols = NULL;
		phgFree(row->data);
		row->data = NULL;
		row->ncols = row->alloc = 0;
	    }
	}
    }

    /* now sbuff contains packed off-process data */
    for (rank = 0; rank < mat->rmap->nprocs - 1; rank++)
	scounts[rank] = sdispls[rank + 1] - sdispls[rank];
    scounts[rank] = (int)(ssize / trunc_size) - sdispls[rank];

    /* exchange counts with other processes */
    MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, mat->rmap->comm);

    /* prepare receive buffer */
    rdispls[0] = 0;
    for (rank = 0; rank < mat->rmap->nprocs - 1; rank++)
	rdispls[rank + 1] = rdispls[rank] + rcounts[rank];
    rsize = (rdispls[rank] + rcounts[rank]) * trunc_size;
    rbuff = phgAlloc(rsize);

    /* exchange data */
    MPI_Type_contiguous((int)trunc_size, MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(sbuff, scounts, sdispls, type,
		  rbuff, rcounts, rdispls, type, mat->rmap->comm);
    MPI_Type_free(&type);

    phgFree(sbuff);

    phgFree(mat->O2Gmap);
    phgFree(mat->ordering);
    if (mat->localsize > mat->rmap->localsize)
	mat->rows = phgRealloc_(mat->rows,
				mat->rmap->localsize * sizeof(*mat->rows),
				mat->rmap->localsize * sizeof(*mat->rows));
    mat->localsize = mat->cmap->localsize;

    /* convert negated global indices in local rows to local indices */
    for (i = 0, row = mat->rows; i < mat->rmap->nlocal; i++, row++) {
	if ((ncols = row->ncols) == 0)
	    continue;
	for (j = 0; j < ncols; j++) {
	    r = row->cols[j];
	    if (r >= 0) {
#if DEBUG
		if (r >= mat->cmap->localsize) {
		    phgError(1, "%s:%d, unexpected!\n", __FILE__, __LINE__);
		}
#endif	/* DEBUG */
		continue;
	    }
	    /* search the O2Gmap */
	    r = -1 - r;
	    assert(r < cindex0 || r >= cindex0 + mat->cmap->nlocal);
	    k = phgBinarySearchINT(O2Gmap_size, O2Gmap, ordering, r);
	    if (k < O2Gmap_size && r == O2Gmap[ordering[k]]) {
		/* found in O2Gmap */
		row->cols[j] = mat->cmap->nlocal + ordering[k];
		continue;
	    }
	    /* add a new global index to O2Gmap */
	    if (O2Gmap_alloc <= O2Gmap_size) {
		O2Gmap = phgRealloc_(O2Gmap, 
				(O2Gmap_alloc + 128) * sizeof(*O2Gmap),
				O2Gmap_alloc * sizeof(*O2Gmap));
		ordering = phgRealloc_(ordering,
				(O2Gmap_alloc + 128) * sizeof(*ordering),
				O2Gmap_alloc * sizeof(*ordering));
		O2Gmap_alloc += 128;
	    }
	    row->cols[j] = mat->cmap->nlocal + O2Gmap_size;
	    O2Gmap[O2Gmap_size] = r;
	    /* insert O2Gmap_size into ordering[] */
	    for (r = O2Gmap_size; r > k; r--)
		ordering[r] = ordering[r - 1];
	    ordering[k] = O2Gmap_size++;
	}
    }

    /* process received data */
    for (rank = 0; rank < mat->rmap->nprocs; rank++) {
	if (rcounts[rank] == 0)
	    continue;
	prows = (void *)(rbuff + ((size_t)rdispls[rank]) * trunc_size);
	nrows = prows[0];
	size = ((size_t)prows[1]) + ((size_t)prows[2]) * (size_t)(INT_MAX);

	prows += 4;
	pncols = prows + nrows;
	pcols = pncols + nrows;

	align = (void *)(pcols + size);
	align = rbuff + ((align - rbuff + align_float - 1) / align_float)
	    * align_float;
	pdata = (void *)align;

	for (i = 0; i < nrows; i++) {
	    ncols = *pncols;
	    /* convert global indices in pcols[] to local indices */
	    for (j = 0; j < ncols; j++) {
		r = pcols[j];
		if (r >= cindex0 && r < cindex0 + mat->cmap->nlocal) {
		    pcols[j] -= cindex0;
		    continue;
		}
		/* search the O2Gmap */
		k = phgBinarySearchINT(O2Gmap_size, O2Gmap, ordering, r);
		if (k < O2Gmap_size && r == O2Gmap[ordering[k]]) {
		    /* found in O2Gmap */
		    pcols[j] = mat->cmap->nlocal + ordering[k];
		    continue;
		}
		/* add a new global index to O2Gmap */
		if (O2Gmap_alloc <= O2Gmap_size) {
		    O2Gmap = phgRealloc_(O2Gmap,
				(O2Gmap_alloc + 128) * sizeof(*O2Gmap),
				O2Gmap_alloc * sizeof(*O2Gmap));
		    ordering = phgRealloc_(ordering,
				(O2Gmap_alloc + 128) * sizeof(*ordering),
				O2Gmap_alloc * sizeof(*ordering));
		    O2Gmap_alloc += 128;
		}
		pcols[j] = mat->cmap->nlocal + O2Gmap_size;
		O2Gmap[O2Gmap_size] = r;
		/* insert O2Gmap_size into ordering[] */
		for (r = O2Gmap_size; r > k; r--)
		    ordering[r] = ordering[r - 1];
		ordering[k] = O2Gmap_size++;
	    }
	    r = *(prows++) - rindex0;
	    assert(r >= 0 && r < mat->rmap->nlocal);
	    add_columns_to_row(mat, r, ncols, pcols, pdata, FALSE, FALSE);
	    pncols++;
	    pcols += ncols;
	    pdata += ncols;
	}
    }
    mat->O2Gmap = O2Gmap;
    mat->ordering = ordering;
    mat->localsize = mat->cmap->nlocal + O2Gmap_size;

    phgFree(scounts);
    phgFree(rbuff);
#else	/* USE_MPI */
    if (mat->assembled)
	return;

    if (mat->blocks != NULL) {
	assemble_block_matrix(mat);
	update_nnz(mat);
	free_locks(mat);
	return;
    }

    save_local_bdry_equations(mat);
    mat->localsize = mat->cmap->localsize;
#endif /* USE_MPI */

    build_bdry(mat);
    save_bdry_columns(mat, TRUE);
    mat->assembled = TRUE;
    update_nnz(mat);

    free_locks(mat);
    return;
}

#define dofindex2vecindex(map, dof_ind, vec_ind) {			\
    INT i_ = (dof_ind);							\
    if (map->L2Vmap != NULL && i_ >= map->L2Vmap_size) {		\
	phgInfo(-1, "AddEntries: invalid index: %"dFMT" (%s:%d)\n"	\
		"(note: AddEntries uses local indices)\n",		\
		i_, __FILE__, __LINE__);				\
	phgError(1, "AddEntries: abort.\n");				\
    }									\
    (vec_ind) = (map->L2Vmap == NULL ? i_ : map->L2Vmap[i_]);		\
}

static inline INT
convert_col_index(MAT *mat, INT col, BOOLEAN col_is_global)
/* Converts a local or global col index to vector index.
 *
 * When col_is_global == FALSE, a negetive return value means a removed entry
 * (e.g., a boundary entry in a map with boundary entries removed), which
 * should be ignored (skipped).
 *
 * When col_is_global == TRUE, a negetive return value means a negated global
 * index (-col-1), which means a non local column index not in the array
 * O2Gmap[] thus cannot be converted to vector index. The negated global
 * indices should be added to the matrix as is and will be handled in
 * phgMatAssemble.
 */
{
    INT col1;

#if USE_MPI
    INT i, j, *O2Gmap, *ordering;

    if (col < 0)
	return col;

    if (mat->cmap->nprocs > 1) {
	if (col_is_global) {
	    O2Gmap = mat->cmap->O2Gmap;
	    ordering = mat->cmap->ordering;
	    if (col < 0 || col >= mat->cmap->nglobal)
		phgError(1, "%s (%s:%d): col %"dFMT" out of range "
			    "[0,%"dFMT").\n", __func__, __FILE__, __LINE__,
			    col, mat->cmap->nglobal);
	    if (col >= (j = mat->cmap->partition[mat->cmap->rank]) &&
		col < mat->cmap->partition[mat->cmap->rank + 1]) {
		col1 = col - j;
	    }
	    else {
		j = mat->cmap->localsize - mat->cmap->nlocal;
		i =  phgBinarySearchINT(j, O2Gmap, ordering, col);
		if (i < j && col == O2Gmap[ordering[i]])
		    col1 = mat->cmap->nlocal + ordering[i];
		else
		    col1 = -1 - col;
	    }
	}
	else {
	    dofindex2vecindex(mat->cmap, col, col1);
	}
    }
    else
#endif	/* USE_MPI */
    {
	if (col < 0)
	    return col;

	if (col_is_global)
	    col1 = col;
	else
	    dofindex2vecindex(mat->cmap, col, col1);
	if (col1 >= mat->cmap->localsize) {
	   phgError(1, "%s: unexpected!\n", __func__);
	}
    }

    return col1;
}

static inline void
check_mat_localsize(MAT *mat)
{
    if (mat->localsize >= 0)
	return;

#if USE_OMP
//# pragma omp flush //(mat->localsize, mat->O2Gmap, mat->ordering)
# pragma omp critical (mat_localsize)
#endif  /* USE_OMP */
    if (mat->localsize < 0) {
	/* initialize mat->O2Gmap, mat->ordering, etc. */
	assert(mat->bdry_localsize == 0);
	/* Note: mat->bdry_localsize temporarily stores allocated size of
	 * O2Gmap and ordering */
	mat->bdry_localsize = mat->rmap->localsize - mat->rmap->nlocal;
#if USE_MPI
	phgFree(mat->O2Gmap);
	phgFree(mat->ordering);
	if (mat->rmap->O2Gmap != NULL) {
	    size_t size;
	    size = (mat->rmap->localsize - mat->rmap->nlocal)
			* sizeof(*mat->O2Gmap);
	    mat->O2Gmap = phgAlloc(size);
	    memcpy(mat->O2Gmap, mat->rmap->O2Gmap, size);
	    mat->ordering = phgAlloc(size);
	    memcpy(mat->ordering, mat->rmap->ordering, size);
	}
	else {
	    mat->O2Gmap = NULL;
	    mat->ordering = NULL;
	}
#endif	/* USE_MPI */
	mat->localsize = mat->rmap->localsize;
#if USE_OMP
//# pragma omp flush //(mat->localsize, mat->O2Gmap, mat->ordering)
#endif  /* USE_OMP */
    }
}

static inline INT
convert_row_index(MAT *mat, INT row, BOOLEAN row_is_global)
/* converts a local or global row index to vector index.
 * Note: before assembly, mat->O2Gmap, mat->ordering, and mat->localsize match
 * rmap entries, which are used to store off-proc rows. mat->O2Gmap,
 * mat->ordering, and mat->rows are dynamically enlarged for new off-proc rows
 * and mat->bdry_localsize stores the allocated size of O2Gmap, ordering, and
 * rows. These entries will be cleared in phgMatAssemble before assembly.
 */
{
    INT row1;

#if USE_MPI
    INT j, O2Gmap_size;
    MAT_ROW *mat_row;

    if (row < 0)
	return row;

#if USE_OMP
    /* FIXME: with this lock,
     * 		env OMP_NUM_THREADS=2 ./ipdg -refine 12
     * 		env OMP_NUM_THREADS=4 ./ipdg -refine 12
     * is much slower (0.47s->0.54s, 0.25s->0.34s), but without it,
     *		env OMP_NUM_THREADS=2 mpirun -np 2 ./interface-p4est -refine 4
     * randomly (mostly) crashes. */
    omp_set_lock(&mat->lock);
    check_mat_localsize(mat);
    omp_unset_lock(&mat->lock);
#else	/* USE_OMP */
    check_mat_localsize(mat);
#endif	/* USE_OMP */
	
    if (mat->rmap->nprocs > 1) {
	if (row_is_global) {
	    if (row >= (j = mat->rmap->partition[mat->rmap->rank]) &&
		row < mat->rmap->partition[mat->rmap->rank + 1]) {
		row1 = row - j;
	    }
	    else {
		check_mat_localsize(mat);
		O2Gmap_size = mat->localsize - mat->rmap->nlocal;
		/* search mat->O2Gmap */
		row1 = phgBinarySearchINT(O2Gmap_size, mat->O2Gmap,
					   mat->ordering, row);
		if (row1 >= O2Gmap_size ||
		    row != mat->O2Gmap[mat->ordering[row1]]) {
		    /* New entry ==> enlarge mat->rows if necessary.
		     * Hacky: allocated size stored in mat->bdry_localsize */
#if USE_OMP
		    omp_set_lock(&mat->lock);
#endif	/* USE_OMP */
		    O2Gmap_size = mat->localsize - mat->rmap->nlocal;
		    if (O2Gmap_size >= mat->bdry_localsize) {
			mat->rows = phgRealloc_(mat->rows,
				sizeof(*mat->rows) *
				    (mat->rmap->nlocal + O2Gmap_size + 1024),
				sizeof(*mat->rows) *
				    (mat->rmap->nlocal + mat->bdry_localsize));
#if USE_OMP
			if (phgMaxThreads > 1)
			    mat->locks = phgRealloc_(mat->locks,
				sizeof(*mat->locks) *
				    (mat->rmap->nlocal + O2Gmap_size + 1024),
				sizeof(*mat->locks) *
				    (mat->rmap->nlocal + mat->bdry_localsize));
#endif	/* USE_OMP */
			mat->O2Gmap = phgRealloc_(mat->O2Gmap,
				sizeof(*mat->O2Gmap) * (O2Gmap_size + 1024),
				sizeof(*mat->O2Gmap) * mat->bdry_localsize);
			mat->ordering = phgRealloc_(mat->ordering,
				sizeof(*mat->ordering) * (O2Gmap_size + 1024),
				sizeof(*mat->ordering) * mat->bdry_localsize);
			mat->bdry_localsize = O2Gmap_size + 1024;
		    }
#if USE_OMP
		    if (mat->locks != NULL)
			omp_init_lock(&mat->locks[mat->nlock++]);
#endif	/* USE_OMP */
		    mat->localsize++;
		    /* insert O2Gmap_size to mat->ordering[row1] */
		    for (j = O2Gmap_size; j > row1; j--)
			mat->ordering[j] = mat->ordering[j - 1];
		    mat->ordering[row1] = O2Gmap_size;
		    row1 = O2Gmap_size;
		    /* set entry mat->O2Gmap[row1] */
		    mat->O2Gmap[row1] = row;
		    /* init mat->rows[row1 + mat->rmap->nlocal] */
		    mat_row = mat->rows + row1 + mat->rmap->nlocal;
		    mat_row->ncols = 0;
		    mat_row->alloc = 0;
		    mat_row->cols = NULL;
		    mat_row->data = NULL;
		    row1 += mat->rmap->nlocal;
#if USE_OMP
		    omp_unset_lock(&mat->lock);
#endif	/* USE_OMP */
		}
		else {
		    row1 = mat->rmap->nlocal + mat->ordering[row1];
		}
	    }
	    assert(row1 < mat->rmap->nlocal || row1 < mat->localsize);
	}
	else {
	    check_mat_localsize(mat);
	    dofindex2vecindex(mat->rmap, row, row1);
	    assert(row1 < mat->rmap->localsize);
	}
    }
    else
#endif	/* USE_MPI */
    {
	check_mat_localsize(mat);
	if (row_is_global)
	    row1 = row;
	else
	    dofindex2vecindex(mat->rmap, row, row1);
	assert(row1 < mat->rmap->localsize);
    }

    return row1;
}

void
phgMatSetMode(MAT *mat, MAT_MODE mode)
{
    assert(mode == PHG_ADD || mode == PHG_REPLACE);
    if (mat->mode == mode)
	return;
    if (mat->mode != PHG_UNDEFINED)
	phgError(1, "%s should only be called once and before setting/adding "
		    "any entries to the matrix.\n", __func__);
    mat->mode = mode;

    return;
}

void
phgMatAddEntry_(MAT *mat, INT row, INT col, FLOAT value,
		BOOLEAN row_is_global, BOOLEAN col_is_global,
		BOOLEAN replace)
/* if replace == TRUE then replace the matrix entry with the new value,
 * otherwise add the new value to the existing entry.
 *
 * WARNING: replace == TRUE is a limited implementation (only the local rows
 *	    are set), use with care!
 *
 * Note: before assembly, mat->O2Gmap/mat->ordering/mat->localsize refer to
 *	 off-proc rows. After assembly, they refer to off-proc columns.
 */
{
    static FLOAT zero = 0.0;
    INT row1, col1;

    if (row < 0 || col < 0 || (/*!replace &&*/ value == zero))
	return;

    if (mat->type == PHG_MATRIX_FREE)
	phgError(1, "%s: can't add entries to matrix-free matrix.\n", __func__);

    if (mat->assembled)
	phgError(1, "%s: can't add entries to assembled matrix.\n", __func__);

    if (mat->mode == PHG_UNDEFINED)
	mat->mode = PHG_ADD;

    if (replace && mat->mode == PHG_ADD) {
	phgInfo(-1, "Please call phgMatSetReplaceMode() before calling "
		    "phgMatSetEntries().\n");
	phgError(1, "abort.\n");
    }
    else if (!replace && mat->mode == PHG_REPLACE) {
	phgError(1, "%s: phgMatAddEntries/phgMatSetEntries can't be used "
		    "at the same time.\n", __func__);
    }

    col1 = convert_col_index(mat, col, col_is_global);
    /* Note: if col_is_global == TRUE and col1 < 0, then it's a negated global
     * index, otherwise it's a removed (boundary) entry */
    if (col1 < 0 && !col_is_global)
	return;
    row1 = convert_row_index(mat, row, row_is_global);
    if (row1 < 0)
	return;
    add_columns_to_row(mat, row1, 1, &col1, &value, FALSE, replace);

#if 0
    assert(!row_is_global && !col_is_global);
    if (col != row && (solver->oem_solver == SOLVER_SUPERLU
			   /*|| solver->oem_solver == SOLVER_SPC */ )) {
	/* enforce symmetry of the sparsity */
	dofindex2vecindex(mat->cmap, col, row1);
	dofindex2vecindex(mat->rmap, row, col1);
	add_columns_to_row(mat, row1, 1, &col1, &zero, FALSE);
    }
#endif /* 0|1 */

    return;
}

void
phgMatAddEntries_(MAT *mat,
			INT nrows, const INT *rows,
			INT ncols, const INT *cols,
			const FLOAT *values,
			BOOLEAN rows_are_global,
			BOOLEAN cols_are_global,
			BOOLEAN replace)
/* if replace == TRUE then replace the matrix entries with the new values,
 * otherwise add the new values to existing entries.
 *
 * WARNING: replace == TRUE is a limited implementation (only the local rows
 *	    are set), use with care! */
{
    static INT *allocated = NULL;
    static INT **cols1_ptr = NULL;
    static FLOAT **data1_ptr = NULL;

    static FLOAT zero = 0.0;
    INT i, j, n, row1, *cols1;
    const INT *cols0;
    FLOAT a, *data1;
    const FLOAT *v;

    if (allocated == NULL) {
#if USE_OMP
#pragma omp critical (mat_add_entry)
#endif	/* USE_OMP */
	if (allocated == NULL) {
	    cols1_ptr = phgCalloc(phgMaxThreads, sizeof(*cols1_ptr));
	    data1_ptr = phgCalloc(phgMaxThreads, sizeof(*data1_ptr));
	    FreeAtExitNoCheck_(cols1_ptr, 1);
	    FreeAtExitNoCheck_(data1_ptr, 1);
	    for (i = 0; i < phgMaxThreads; i++) {
		FreeAtExitNoCheck(cols1_ptr[i]);
		FreeAtExitNoCheck(data1_ptr[i]);
	    }
	    allocated = phgCalloc(phgMaxThreads, sizeof(*allocated));
	    FreeAtExitNoCheck(allocated);
	}
    }

    if (ncols <= 0 || nrows <= 0)
	return;

    if (ncols == 1 && nrows == 1) {
	phgMatAddEntry_(mat, *rows, *cols, *values, rows_are_global,
			cols_are_global, replace);
	return;
    }

    if (mat->mode == PHG_UNDEFINED)
	mat->mode = PHG_ADD;

    if (replace && mat->mode == PHG_ADD) {
	phgInfo(-1, "Please call phgMatSetReplaceMode() before calling "
		    "phgMatSetEntries().\n");
	phgError(1, "abort.\n");
    }
    else if (!replace && mat->mode == PHG_REPLACE) {
	phgError(1, "%s: phgMatAddEntries/phgMatSetEntries can't be used "
		    "at the same time.\n", __func__);
    }

    if (mat->type == PHG_MATRIX_FREE)
	phgError(1, "%s: can't add entries to matrix-free matrix.\n", __func__);

    if (mat->assembled)
	phgError(1, "%s: can't add entries to assembled matrix.\n", __func__);

    if (allocated[phgThreadId] < ncols) {
	phgFree(data1_ptr[phgThreadId]);
	phgFree(cols1_ptr[phgThreadId]);
	allocated[phgThreadId] = ncols + 100;
	data1_ptr[phgThreadId] = phgAlloc(allocated[phgThreadId]
					* sizeof(*data1_ptr[phgThreadId]));
	cols1_ptr[phgThreadId] = phgAlloc(2 * allocated[phgThreadId]
					* sizeof(*cols1_ptr[phgThreadId]));
    }
    cols1 = cols1_ptr[phgThreadId];
    data1 = data1_ptr[phgThreadId];

    /* convert column indices to local ones */
    if (mat->cmap->nprocs > 1 || mat->cmap->L2Vmap != NULL) {
	INT *p = cols1 + ncols;
	for (i = 0; i < ncols; i++) {
	    if (cols[i] >= 0)
		p[i] = convert_col_index(mat, cols[i], cols_are_global);
	    else
		p[i] = cols[i];
	}
	cols0 = p;
    }
    else {
	cols0 = cols;
#if DEBUG
	for (i = 0; i < ncols; i++)
	    if (cols0[i] >= mat->cmap->localsize) {
		phgError(1, "%s:%d (%s): column index out of range!\n",
				__FILE__, __LINE__, __func__);
	    }
#endif	/* DEBUG */
    }

    v = values;
    for (i = 0; i < nrows; i++) {
	if (rows[i] < 0)
	    continue;
	for (j = 0, n = 0; j < ncols; j++) {
	    a = *(v++);
	    /* Note: if col_is_global == TRUE and col0[j] < 0, then it's a
	     * negated global index, otherwise it's a remove entry. */
	    if ((cols0[j] < 0 && !cols_are_global) ||
		(/*!replace &&*/ a == zero))
		continue;
	    cols1[n] = cols0[j];
	    data1[n++] = a;
	}
	if (n == 0)
	    continue;
	row1 = convert_row_index(mat, rows[i], rows_are_global);
	if (row1 < 0)
	    continue;
	add_columns_to_row(mat, row1, n, cols1, data1, FALSE, replace);
    }

#if 0
    /* FIXME: skip boundary DOF */
    if (solver->oem_solver == SOLVER_SUPERLU
	/*|| solver->oem_solver == SOLVER_SPC */ ) {
	/* enforce symmetry of the sparsity */
	INT col;
	assert(!rows_are_global && !cols_are_global);
	for (i = 0; i < nrows; i++) {
	    dofindex2vecindex(mat->rmap, rows[i], col);
	    for (j = 0; j < ncols; j++) {
		if (*(values++) == zero || cols1[j] == col)
		    continue;
		dofindex2vecindex(mat->cmap, cols[j], row1);
		add_columns_to_row(mat, row1, 1, &col, &zero, FALSE);
	    }
	}
    }
#endif /* 0|1 */

    return;
}

static void
solve_local_bdry_equations(MAT *mat, VEC *rhs)
{
    int i, n = 0;
    INT j, k;
    FLOAT *x = NULL;
    BDRY_EQNS *be;

    if (!mat->handle_bdry_eqns)
	return;

    if (mat == NULL || rhs == NULL)
	return;

    phgInfo(1, "solving local boundary equations.\n");

    assert(!rhs->assembled);

    if (mat->blocks != NULL) {
	/* Important: the user has to ensure rhs matches the diagonal blocks */
	int p = mat->blocks->p, q = mat->blocks->q;
	INT offpsize = 0, nlocal = 0;

	assert(p == q);
	for (i = 0; i < p; i++) {
	    MAT *diag = mat->blocks->pmat[i * q + i];
	    VEC *rhs0;

	    if (diag == NULL) {
		nlocal += mat->blocks->nrows[i];
		offpsize += mat->blocks->noffp[i];
		continue;
	    }

	    if (diag->blocks != NULL)
		phgError(1, "(%s:%d) unimplemented.\n", __FILE__, __LINE__);

	    /* construct RHS for the diagonal block */
	    rhs0 = phgMapCreateVec(diag->cmap, rhs->nvec);
	    assert(rhs->nvec == 1);
	    for (j = 0; j < diag->cmap->nlocal; j++)
		rhs0->data[j] = rhs->data[nlocal + j];
	    for (j = 0; j < diag->cmap->localsize - diag->cmap->nlocal; j++)
		rhs0->offp_data[j] = rhs->offp_data[j + offpsize];
	    rhs0->assembled = FALSE;
	    solve_local_bdry_equations(diag, rhs0);
	    /* update global RHS */
	    for (j = 0; j < diag->cmap->nlocal; j++)
		rhs->data[nlocal + j] = rhs0->data[j];
	    for (j = 0; j < diag->cmap->localsize - diag->cmap->nlocal; j++)
		rhs->offp_data[offpsize + j] = rhs0->offp_data[j];
	    offpsize += diag->cmap->localsize - diag->cmap->nlocal;
	    nlocal += diag->cmap->nlocal;
	    phgVecDestroy(&rhs0);
	}
	assert(nlocal == rhs->map->nlocal);
	assert(offpsize == rhs->map->localsize - rhs->map->nlocal);

	return;
    }

    if (mat->bdry_eqns_count == 0)
	return;
    be = mat->bdry_eqns;
    for (k = 0; k < mat->bdry_eqns_count; k++, be++) {
	if (n < be->n) {
	    n = be->n;
	    phgFree(x);
	    x = phgAlloc(n * sizeof(*x));
	}
	for (i = 0; i < be->n; i++) {
	    if ((j = be->map[i]) < rhs->map->nlocal)
		x[i] = rhs->data[j];
	    else
		x[i] = rhs->offp_data[j - rhs->map->nlocal];
	}
	phgSolverDenseSV(be->n, (void *)be->lu, be->pvt, 1, (void *)x);
	for (i = 0; i < be->n; i++) {
	    if ((j = be->map[i]) < rhs->map->nlocal)
		rhs->data[j] = x[i];
	    else
		rhs->offp_data[j - rhs->map->nlocal] = x[i];
	}
    }

    phgFree(x);
}

void
phgVecAssemble(VEC *vec)
{
#if USE_MPI
    /* flag indicating whether the user has added extra offproc entries */
    int flag0, flag;
    COMM_INFO *cinfo;
#endif	/* USE_MPI */

    if (vec->assembled)
	return;

    if (vec->mat != NULL)
        solve_local_bdry_equations(vec->mat, vec);
    vec->assembled = TRUE;

#if USE_MPI
    if (vec->map->nprocs <= 1 || vec->nvec <= 0) {
	phgFree(vec->O2Gmap);
	vec->O2Gmap = NULL;
	phgFree(vec->ordering);
	vec->ordering = NULL;
	vec->alloc = 0;
	vec->localsize = 0;
	return;
    }

    assert(vec->localsize < 0 || vec->localsize >= vec->map->localsize);
    flag0 = (vec->localsize > vec->map->localsize ? 1 : 0);
    MPI_Allreduce(&flag0, &flag, 1, MPI_INT, MPI_MAX, vec->map->comm);

    if (!flag) {
	/* no extra offproc entries in vec, use vec->map->cinfo */
	if (vec->map->cinfo == NULL)
	    vec->map->cinfo = phgMapCreateCommInfo(
				vec->map->comm,
				vec->map->nprocs,
				vec->map->rank,
				vec->map->nlocal,
				vec->map->localsize,
				vec->map->O2Gmap,
				vec->map->ordering,
				vec->map->partition);
	cinfo = vec->map->cinfo;
    }
    else {
	/* some procs have extra offproc entries, create a custom cinfo */
	cinfo = phgMapCreateCommInfo(
				vec->map->comm,
				vec->map->nprocs,
				vec->map->rank,
				vec->map->nlocal,
				flag0 ? vec->localsize : vec->map->localsize,
				vec->O2Gmap,
				vec->ordering,
				vec->map->partition);
    }
    phgMapGather(cinfo, vec->nvec, vec->data, vec->offp_data);
    if (cinfo != vec->map->cinfo)
	phgMapDestroyCommInfo(&cinfo);

    phgFree(vec->O2Gmap);
    vec->O2Gmap = NULL;
    phgFree(vec->ordering);
    vec->ordering = NULL;
    vec->alloc = 0;
    if (vec->localsize > vec->map->localsize) {
	phgFree(vec->offp_data);
	vec->offp_data = phgCalloc(vec->map->localsize,
					sizeof(*vec->offp_data));
    }
    vec->localsize = 0;
#endif /* USE_MPI */

    return;
}

void
phgVecAddEntries_(VEC *vec, int which,
		  INT n, const INT *indices, const FLOAT *values,
		  BOOLEAN global)
{
    INT i, j, row, offp_size;
    FLOAT *data, *offp_data;
    INT index;

    assert(which < vec->nvec && vec->assembled == FALSE);
    if (n == 0)
	return;

#if USE_MPI
    if (vec->localsize < 0) {
	/* Init vec->O2Gmap and vec->ordering */
#if USE_OMP
#pragma omp critical (vec_add_entry)
	{
#endif	/* USE_OMP */
	offp_size = vec->map->localsize - vec->map->nlocal;
	if (offp_size > 0) {
	    assert(vec->map->O2Gmap != NULL && vec->map->ordering != NULL);
	    phgFree(vec->O2Gmap);
	    phgFree(vec->ordering);
	    vec->alloc = 0;
	    vec->O2Gmap = phgAlloc(offp_size * sizeof(*vec->O2Gmap));
	    vec->ordering = phgAlloc(offp_size * sizeof(*vec->ordering));
	    memcpy(vec->O2Gmap, vec->map->O2Gmap,
				offp_size * sizeof(*vec->O2Gmap));
	    memcpy(vec->ordering, vec->map->ordering,
				offp_size * sizeof(*vec->ordering));
	}
	vec->localsize = vec->map->localsize;
#if USE_OMP
	}
#endif	/* USE_OMP */
    }
    else {
	offp_size = vec->localsize - vec->map->nlocal;
    }
#else	/* USE_MPI */
    offp_size = vec->map->localsize - vec->map->nlocal;
#endif	/* USE_MPI */

    data = vec->data + which * vec->map->nlocal;
    offp_data = vec->offp_data + which * offp_size;
    for (j = 0; j < n; j++, values++, indices++) {
	if (*values == 0.0 || *indices < 0)
	    continue;
	if (global) {
	    index = *indices;
	    assert(index >= 0 && index < vec->map->nglobal);
	    if (index >= (i = vec->map->partition[vec->map->rank]) &&
		index < vec->map->partition[vec->map->rank + 1]) {
		row = index - i;
	    }
#if USE_MPI
	    else {
		/* search vec->O2Gmap */
		row = phgBinarySearchINT(offp_size, vec->O2Gmap, vec->ordering,
				      index);
		if (row >= offp_size ||
		    index != vec->O2Gmap[vec->ordering[row]]) {
		    CheckThread
		    /* new entry */
		    if (offp_size >= vec->alloc) {
			/* Enlarge offp_size */
			vec->offp_data = phgRealloc_(vec->offp_data,
				sizeof(*vec->offp_data) * (offp_size + 1024),
				sizeof(*vec->offp_data) * vec->alloc);
			vec->O2Gmap = phgRealloc_(vec->O2Gmap,
				sizeof(*vec->O2Gmap) * (offp_size + 1024),
				sizeof(*vec->O2Gmap) * vec->alloc);
			vec->ordering = phgRealloc_(vec->ordering,
				sizeof(*vec->ordering) * (offp_size + 1024),
				sizeof(*vec->ordering) * vec->alloc);
			vec->alloc = offp_size + 1024;
		    }
		    vec->localsize++;
		    /* insert O2Gmap_size to mat->ordering[row] */
		    for (i = offp_size; i > row; i--)
			    vec->ordering[i] = vec->ordering[i - 1];
		    vec->ordering[row] = offp_size;
		    row = offp_size++;
		    /* set entry vec->O2Gmap[row] */
		    vec->O2Gmap[row] = index;
		    /* FIXME: store offp_data as [vec->alloc][] instead of
		     *	      [offp_size][]? */
		    assert(vec->nvec == 1);
		    offp_data = vec->offp_data + which * offp_size;
		    /* reset value in the new entry */
		    offp_data[row] = 0.0;
		    row += vec->map->nlocal;
		}
		else {
		    row = vec->map->nlocal + vec->ordering[row];
		}
	    }
	    assert(row >= 0 && row < vec->localsize);
#else	/* USE_MPI */
	    else row = 0;	/* make gcc happy */
	    assert(row >= 0 && row < vec->map->localsize);
#endif	/* USE_MPI */
	}
	else {
	    dofindex2vecindex(vec->map, *indices, row);
	    assert(row < vec->map->localsize);
	}

	if (row < 0)
	    continue;

#if USE_OMP
#pragma omp critical (vec_add_entry1)
#endif	/* USE_OMP */
	if (row < vec->map->nlocal)
	    data[row] += *values;
	else
	    offp_data[row - vec->map->nlocal] += *values;
    }

    return;
}

void
phgVecAddEntry(VEC *vec, int which, INT index, FLOAT value)
{
    phgVecAddEntries(vec, which, 1, &index, &value);
}

void
phgVecAddGlobalEntry(VEC *vec, int which, INT index, FLOAT value)
{
    phgVecAddGlobalEntries(vec, which, 1, &index, &value);
}

/*--------------------------------------------------------------------------*/

MAT *
phgMatRemoveBoundaryEntries(MAT *mat)
/* removes boundary entries from mat */
{
    INT i, j, k, m;
    MAT_ROW *row, *row1;
    INT *new_index = NULL;
#if USE_MPI
    int ii;
    INT n, n0, m0, *sbuff, *rbuff, *O2Gmap, *ordering;
    int *scounts, *sdispls, *rcounts, *rdispls;
#endif
    MAP *map;

    FunctionEntry;

    if (!mat->assembled)
	phgMatAssemble(mat);

    assert(mat->type != PHG_DESTROYED);

    if (mat->rmap->bdry_nglobal == 0 && mat->cmap->bdry_nglobal == 0)
	Return mat;

    if (mat->type == PHG_PACKED)
	phgMatUnpack(mat);

    if (mat->refcount > 0 && mat->rmap->rank == 0 && mat->cmap->rank == 0)
	phgWarning("matrix referenced by other objects.\n");

    if (mat->rmap->bdry_nglobal != 0) {
	map = phgMapRemoveBoundaryEntries(mat->rmap);
	if (mat->type != PHG_MATRIX_FREE) {
	    /* remove boundary rows */
	    for (i = 0; i < mat->rmap->bdry_nlocal; i++) {
		row = mat->rows + mat->rmap->bdry[i];
		phgFree(row->cols);
		row->cols = NULL;
		phgFree(row->data);
		row->data = NULL;
		row->alloc = 0;
		row->ncols = -1;
	    }
	    row = row1 = mat->rows;
	    for (i = 0; i < mat->rmap->nlocal; i++, row++) {
		if (row->ncols < 0) {
		    /* deleted row */
		    row->ncols = 0;
		    continue;
		}
		if (row1 != row) {
		    /* copy row to row1 */
		    row1->cols = row->cols;
		    row1->data = row->data;
		    row1->alloc = row->alloc;
		    row1->ncols = row->ncols;
		    row->ncols = row->alloc = 0;
		    row->cols = NULL;
		    row->data = NULL;
		}
	 	row1++;
	    }
	    assert(row1 - mat->rows == map->nlocal);
	    mat->rows = phgRealloc_(mat->rows, map->localsize * sizeof(*row),
					       map->localsize * sizeof(*row));
	}

	phgFree(mat->diag);
	mat->diag = NULL;
	phgFree(mat->diag1);
	mat->diag1 = NULL;
	free_matrix_rows(mat->rmap->nlocal, &mat->bdry_rows);
	phgMapDestroy(&mat->rmap);
	mat->rmap = map;
    }

    mat->handle_bdry_eqns = FALSE;
    free_bdry_eqns(mat);

    if (mat->cmap->bdry_nglobal == 0) {
	update_nnz(mat);
	Return mat;
    }

    map = phgMapRemoveBoundaryEntries(mat->cmap);

    if (mat->type == PHG_MATRIX_FREE) {
	phgMapDestroy(&mat->cmap);
	mat->cmap = map;
	goto end;
    }

    /* remove boundary columns */

    /* expand mat->bdry array for faster checking */
    new_index = phgAlloc(mat->localsize * sizeof(*new_index));

    /* Note: the interval [nlocal, localsize) may include
     * unreferenced off-process entries (some of them may have been deleted by
     * save_bdry_columns()), which should be excluded from the communications
     * for getting new indices of off-process entries. So we proceed with the
     * following steps:
     *
     *  Step 1: initialize all entries in new_index to -1.
     *
     *  Step 2: set all boundary entries to -2.
     *
     *  Step 3: set all non boundary entries of new_index in [0, nlocal)
     *          to new local indices for 'mat'.
     *
     *  Step 4: scan rows, copy non boundary entries to 'mat'.
     *          The column indices for entries in [0, nlocal) are
     *          converted to new local indices for 'mat', and the entries
     *          in [nlocal, localsize) are unchanged
     *          (they only exist and require communications when nprocs > 1).
     */

    /* initialize all new_index entries to -1 */
    for (i = 0; i < mat->localsize; i++) {
	new_index[i] = -1;
    }
    /* set boundary entries in new_index to -2 */
    for (i = 0; i < mat->bdry_localsize; i++) {
	assert(mat->bdry[i] >= 0 && mat->bdry[i] < mat->localsize);
	new_index[mat->bdry[i]] = -2;
    }
    /* convert non boundary new_index entries to new local indices for mat */
    for (i = 0, k = 0; i < mat->cmap->nlocal; i++) {
	if (new_index[i] == -1)
	    new_index[i] = k++;
    }
    assert(k == map->nlocal);

    /* new_index >= 0 means a non boundary entry present in mat, and is its new
     * local index; new_index == -1 means an offprocess non boundary entry;
     * new_index == -2 means a boundary entry */

    for (i = 0, row = mat->rows; i < mat->rmap->nlocal; i++, row++) {
	/* count non boundary columns in current row */
	for (j = 0, k = 0; j < row->ncols; j++) {
	    m = row->cols[j];
	    assert(m >= 0 && m < mat->localsize);
	    if (new_index[m] > -2) {
		if (m >= mat->cmap->nlocal) {
		    new_index[m] = 0;
		    row->cols[k] = m;
		}
		else {
		    row->cols[k] = new_index[m];;
		    assert(row->cols[k] < map->nlocal);
		}
		row->data[k++] = row->data[j];
	    }
	}
	row->ncols = k;
	if (k == 0) {
	    phgFree(row->cols);
	    phgFree(row->data);
	    row->cols = NULL;
	    row->data = NULL;
	    row->alloc = 0;
	}
    }

    if (map->nprocs <= 1 || map->nglobal == 0) {
	mat->localsize = map->nlocal;
	goto end;
    }

#if USE_MPI
    m0 = map->partition[map->rank]; /* first new global column index of mat */

    /* count total number of non-boundary off-process columns */
    for (m = 0, i = mat->cmap->nlocal; i < mat->localsize; i++) {
	if (new_index[i] >= 0) {
	    new_index[i] = m++;	/* new_index[i] = index of column i */
	}
    }
    sbuff = phgAlloc(m * sizeof(*sbuff));
    ordering = NULL;
    O2Gmap = phgAlloc(m * sizeof(*O2Gmap));
    if (m > 0) {
	/* collect global indices of non-boundary off-process columns */
	for (k = 0, i = mat->cmap->nlocal; i < mat->localsize; i++) {
	    if (new_index[i] >= 0) {
		O2Gmap[k++] = mat->O2Gmap[i - mat->cmap->nlocal];
	    }
	}
	_phg_build_ordering(m, O2Gmap, &ordering);
	for (i = 0; i < m; i++)
	    sbuff[i] = O2Gmap[ordering[i]];
    }

    scounts = phgAlloc(4 * map->nprocs * sizeof(*scounts));
    sdispls = scounts + map->nprocs * 1;
    rcounts = scounts + map->nprocs * 2;
    rdispls = scounts + map->nprocs * 3;

    /* count number of off-process columns for each process */
    memset(scounts, 0, map->nprocs * sizeof(*scounts));
    for (ii = 1, k = 0; k < m; k++) {
	j = sbuff[k];
	while (mat->cmap->partition[ii] <= j) {
	    ii++;
	    assert(ii <= map->nprocs);
	}
	assert(scounts[ii - 1] < INT_MAX);
	scounts[ii - 1]++;
    }

    /* exchange off-process columns counts */
    MPI_Alltoall(scounts, 1, MPI_INT, rcounts, 1, MPI_INT, map->comm);

    /* compute displacements in the send and receive buffers */
    sdispls[0] = rdispls[0] = 0;
    for (ii = 0; ii < map->nprocs - 1; ii++) {
	sdispls[ii + 1] = sdispls[ii] + scounts[ii];
	rdispls[ii + 1] = rdispls[ii] + rcounts[ii];
	assert(sdispls[ii + 1] <= INT_MAX && rdispls[ii + 1] <= INT_MAX);
    }
    assert(m == sdispls[ii] + scounts[ii]);
    n = rdispls[ii] + rcounts[ii];
    rbuff = phgAlloc(n * sizeof(*rbuff));

    /* send indices to owning process */
    MPI_Alltoallv(sbuff, scounts, sdispls, PHG_MPI_INT,
		  rbuff, rcounts, rdispls, PHG_MPI_INT, map->comm);

    /* convert received indices to new indices for mat */
    n0 = mat->cmap->partition[mat->cmap->rank];
    for (k = 0; k < n; k++) {
	assert(rbuff[k] >= n0 && rbuff[k] < n0 + mat->cmap->nlocal);
	/* the next assertion would fail if the entry rbuff[k] is not
	 * referenced by mat, see note on new_index above. */
	assert(new_index[rbuff[k] - n0] >= 0);
	rbuff[k] = new_index[rbuff[k] - n0] + m0;
	assert(rbuff[k] < m0 + map->nlocal);
    }

    /* send back converted indices */
    MPI_Alltoallv(rbuff, rcounts, rdispls, PHG_MPI_INT,
		  sbuff, scounts, sdispls, PHG_MPI_INT, map->comm);

    for (i = 0; i < m; i++)
	O2Gmap[ordering[i]] = sbuff[i];

    /* Update O2Gmap[] such that its lower part is identical to map->O2Gmap[].
     * The array sbuf[] is temporarily used to store the mapping from
     * the present O2Gmap to the updated O2Gmap (-1 means entries not in
     * map->O2Gmap). */
    for (i = 0; i < m; i++)
	sbuff[i] = -1;    
    for (i = 0; i < map->localsize - map->nlocal; i++) {
	j = phgBinarySearchINT(m, O2Gmap, ordering, map->O2Gmap[i]);
	if (j < m && O2Gmap[ordering[j]] == map->O2Gmap[i]) {
	    /* entry found in O2Gmap[], map it to i */
	    sbuff[ordering[j]] = i;
	}
    }
    /* count # of entries in O2Gmap[] which are not in map->O2Gmap[] */
    for (i = 0, k = map->localsize; i < m; i++)
	if (sbuff[i] == -1)
	    k++;
    mat->localsize = k;
    /* setup updated O2Gmap and ordering */
    phgFree(mat->O2Gmap);
    phgFree(mat->ordering);
    mat->O2Gmap = phgAlloc((k - map->nlocal) * sizeof(*mat->O2Gmap));
    for (i = 0; i < map->localsize - map->nlocal; i++)
	mat->O2Gmap[i] = map->O2Gmap[i];
    for (j = 0, k = map->localsize - map->nlocal; j < m; j++) {
	if (sbuff[j] >= 0) {
	    assert(O2Gmap[j] == map->O2Gmap[sbuff[j]]);
	    continue;
	}
	mat->O2Gmap[i++] = O2Gmap[j];
	sbuff[j] = k++;
    }
    mat->ordering = NULL;
    _phg_build_ordering(i, mat->O2Gmap, &mat->ordering);

    /* convert column indices in mat to new local indices */
    for (i = 0, row = mat->rows; i < mat->rmap->nlocal; i++, row++) {
	for (j = 0; j < row->ncols; j++) {
	    if ((k = row->cols[j]) < map->nlocal)
		continue;
	    assert(new_index[k] < m);
	    row->cols[j] = map->nlocal + sbuff[new_index[k]];
	}
    }

    phgFree(O2Gmap);
    phgFree(ordering);
    phgFree(rbuff);
    phgFree(sbuff);
    phgFree(scounts);
#endif /* USE_MPI */

end:
    phgFree(new_index);
    phgMapDestroy(&mat->cmap);
    mat->cmap = map;

    phgFree(mat->bdry);
    mat->bdry = NULL;
    mat->bdry_localsize = 0;

    update_nnz(mat);

    Return mat;
    return mat;	/* make MIPSPro cc happy */
}

int
phgMatEigenSolve(MAT *A, MAT *B, int n, int which, FLOAT tau, int *nit,
		 FLOAT *evals, MAP *map, DOF **u, ...)
{
    va_list ap;

    /* va_list is not really needed, it avoids a pgcc warning */
    va_start(ap, u);
    va_end(ap);
    (void)ap;
    phgPrintf(
	"****************************************************************\n"
	"* phgMatEigenSolve: the name of the function has been changed. *\n"
	"* Please use phgDofEigenSolve() instead.                       *\n"
	"****************************************************************\n");
    phgFinalize();
    exit(1);

    return 0;	/* make MIPSPro cc happy */
}

int
phgMatDumpCSR(MAT *mat, const char *fn)
/*****************************************************************************
 * This function dumps the matrix in CSR format.
 *
 * [From Li Zhongze]
 *
 * Assume a matrix of size 4 has the following form:
 *
 *	a_{00}   a_{01}   0.0      a_{03}
 *	0.0      a_{11}   a_{12}    0.0
 *	a_{20}   0.0      a_{22}   a_{23}
 *	0.0      0.0      a_{32}   a_{33}
 *
 * The sample CSR file format is as follows (0-based C format):
 *
 *	4
 *	0
 *	3
 *	5
 *	8
 *	10
 *	0
 *	1
 *	3
 *	1
 *	2
 *	0
 *	2
 *	3
 *	2
 *	3
 *	a_{00}  
 *	a_{01}
 *	a_{03}
 *	a_{11}
 *	a_{12}
 *	a_{20} 
 *	a_{22}
 *	a_{23}
 *	a_{32}
 *	a_{33}
 *
 * The first line is the size of the matrix, followed by integer arrays ia
 * of size n+1, and ja of size nnz(the total number of nonzero entries in
 * the matrix, nnz=ia[n]-ia[0]). A real array a of size nnz is stored last
 * in the file.  ia is an array of pointer to the start of each row in
 * array ja and a. Array ja stores the column indices (0-based) of entries
 * on each row. Array a stores the corresponding entry values. 
 *****************************************************************************/
{
    FILE *fp;
    INT i, j;
    MAT_ROW *row;
    char *fn1;

    FunctionEntry;

    if (fn == NULL)
	Return 0;		/* do nothing */

    if (phgNProcs > 1) {
	fn1 = phgAlloc(strlen(fn) + 1 + 5);
	sprintf(fn1, "%s.%04d", fn, phgRank);
	fn = fn1;
    }
    else {
	fn1 = NULL;
    }

    if (!mat->assembled)
	phgMatAssemble(mat);

    assert(mat->type != PHG_DESTROYED);

    if (mat->type == PHG_PACKED)
	phgMatUnpack(mat);

    if ((fp = fopen(fn, "w+t")) == NULL) {
	phgError(0, "cannot open file \"%s\"\n", fn);
	Return 1;
    }
    phgFree(fn1);

    /* pass 1: dump matrix header */
    fprintf(fp, "%"dFMT"\n", mat->rmap->nlocal);
    row = mat->rows;
    for (i = j = 0; i < mat->rmap->nlocal; i++, row++) {
	if (fp != NULL) {
	    fprintf(fp, "%"dFMT"\n", j);
	    j += row->ncols;
	}
    }
    fprintf(fp, "%"dFMT"\n", j);

    /* pass 2: dump column indices */
    row = mat->rows;
    for (i = 0; i < mat->rmap->nlocal; i++, row++)
	for (j = 0; j < row->ncols; j++)
	    fprintf(fp, "%"dFMT"\n", row->cols[j]);

    /* pass 3: dump nonzero matrix entries */
    row = mat->rows;
    for (i = 0; i < mat->rmap->nlocal; i++, row++)
	for (j = 0; j < row->ncols; j++)
	    fprintf(fp, "%0.17lg\n", (double)row->data[j]);

    fclose(fp);
    phgInfo(2, "matrix dumped to file \"%s\".\n", fn);

    Return 0;
    return 0;
}

int
phgVecDump(VEC *vec, const char *fn)
/*****************************************************************************
 * This function dumps a vector to a file.
 *
 * The vector is written as follows (assuming vec size = 4):
 *	4
 *	v_{0}
 *	v_{1}
 *	v_{2}
 *	v_{3}
 *****************************************************************************/
{
    FILE *fp;
    INT i;
    char *fn1;

    FunctionEntry;

    if (fn == NULL)
	Return 0;		/* do nothing */

    if (phgNProcs > 1) {
	fn1 = phgAlloc(strlen(fn) + 1 + 5);
	sprintf(fn1, "%s.%04d", fn, phgRank);
	fn = fn1;
    }
    else {
	fn1 = NULL;
    }

    if (!vec->assembled)
	phgVecAssemble(vec);

    if ((fp = fopen(fn, "w+t")) == NULL) {
	phgError(0, "cannot open file \"%s\"\n", fn);
	Return 1;
    }
    phgFree(fn1);

    fprintf(fp, "%"dFMT"\n", vec->map->nlocal);
    for (i = 0; i < vec->map->nlocal; i++)
	fprintf(fp, "%0.17lg\n", (double)vec->data[i]);

    fclose(fp);
    phgInfo(2, "RHS dumped to file \"%s\".\n", fn);

    Return 0;
    return 0;
}

FLOAT
phgMatConditionNumber_(MAT *A, SYMMETRY symmetry, SOLVER *S, SOLVER *St)
/* returns an estimated condition number of the matrix A.
 *
 * "symmetry" specifies symmetry of the matrix (0 unsym, 1 ssym, 2 sym, 3 spd).
 *
 * "S" and "St" are optional solver for solving A and A^t respectively.
 */
{
    static char *solver_opts =
			"-solver_atol=0 -solver_rtol=0 -solver_btol=1e-5"
			 /*" -solver_maxit=100 "*/;
    static FLOAT tol = 1e-5;	/* tolerance for (inverse) power iterations */
    static INT maxit = 1000;
    static BOOLEAN initialized = FALSE, monitor = FALSE;

    BOOLEAN A_owned = FALSE, At_owned = FALSE;
    MAT *At = NULL;
    VEC *v, *v0, *vt0;
    FLOAT umin, umax;
    double cond, d, t;
    int it;

    /* variables for backup/restore some solver parameters */
    SOLVER *S_bak = NULL, *St_bak = NULL;

    if (!initialized) {
	initialized = TRUE;
	assert(A == NULL);
	phgOptionsRegisterTitle("\nOptions for computing condition number:",
				"\n", "cond");
	phgOptionsRegisterFloat("-cond_tol", "Relative error tolerance", &tol);
	phgOptionsRegisterInt("-cond_maxit",
			"Max number of its for the power method", &maxit);
	phgOptionsRegisterNoArg("-cond_monitor",
			"Show itertions of the power method", &monitor);
	/* Note: for convenience use "append" mode for -cond_solver_opts */
	phgOptionsRegister("-cond_solver_opts",
			"Custom solver options for the inverse power method",
			NULL, &solver_opts, VT_STRING, TRUE);
	return 0.;
    }

    assert(A->rmap->nglobal == A->cmap->nglobal ||
	  (S == NULL && St == NULL && symmetry == M_UNSYM));

    assert(A->type != PHG_DESTROYED);
    phgMatAssemble(A);

    t = phgGetTime(NULL);

    /* set condition solver options */
    phgOptionsPush();
    phgSolverSetDefaultSuboptions();
    phgOptionsSetKeyword("-solver_symmetry", symmetry == M_UNSYM ? "unsym" :
					     symmetry == M_SSYM ? "ssym" :
					     symmetry == M_SYM ? "sym" : "spd");
    phgOptionsSetOptions(solver_opts);  /* user solver options */

#define CopyParams(S_bak, S)			\
	S_bak->atol = S->atol;			\
	S_bak->rtol = S->rtol;			\
	S_bak->btol = S->btol;			\
	S_bak->maxit = S->maxit;		\
	S_bak->monitor = S->monitor;		\
	S_bak->rhs = S->rhs;			\
	S_bak->mat_bak = S->mat_bak;		\
	S_bak->D = S->D;			\
	S_bak->scaling = S->scaling;		\
	S_bak->rhs_updated = S->rhs_updated;

#define SetParams(S)					\
	S->atol = phgOptionsGetFloat("-solver_atol");	\
	S->btol = phgOptionsGetFloat("-solver_btol");	\
	S->rtol = phgOptionsGetFloat("-solver_rtol");	\
	S->maxit = phgOptionsGetInt("-solver_maxit");	\
	S->monitor = phgOptionsGetNoArg("-solver_monitor");	\
	S->mat_bak = NULL;				\
	S->D = NULL;					\
	S->scaling = 0;

    if (S != NULL) {
	S_bak = phgCalloc(1, sizeof(*S_bak));
	CopyParams(S_bak, S);
	SetParams(S);
    }
    if (St != NULL) {
	St_bak = phgCalloc(1, sizeof(*St_bak));
	CopyParams(St_bak, St);
	SetParams(St);
    }

    /* Prepare matrices and solvers for power and inverse power methods */
    if (S == NULL || (symmetry < M_SYM && St == NULL)) {
	if (A->rmap->nglobal != A->cmap->nglobal) {
	    /* non square matrix */
	    if (monitor)
		phgPrintf("*** %s: non square matrix.\n", __func__);
	    if (S == NULL) {
		if (monitor)
		    phgPrintf("*** %s: creating solver S.\n", __func__);
		if (A->rmap->nglobal < A->cmap->nglobal)
		    /* factorize A * A^t */
		    A = phgMatMat(MAT_OP_N, MAT_OP_T, 1., A, A, 0., NULL);
		else
		    /* factorize A^t * A */
		    A = phgMatMat(MAT_OP_T, MAT_OP_N, 1., A, A, 0., NULL);
		A_owned = TRUE;
		phgOptionsSetKeyword("-solver_symmetry", "spd");
		S = phgMat2Solver(SOLVER_DEFAULT, A);
	    }
	}
	else {
	    /* square matrix: factorrize A and A^t */
	    if (monitor)
		phgPrintf("*** %s: square matrix.\n", __func__);
	    if (S == NULL) {
		if (monitor)
		    phgPrintf("*** %s: creating solver S.\n", __func__);
		S = phgMat2Solver(SOLVER_DEFAULT, A);
	    }
	    if (symmetry < M_SYM && St == NULL) {
		if (monitor)
		    phgPrintf("*** %s: creating solver St.\n", __func__);
		At = phgMatTranspose(A);
		St = phgMat2Solver(SOLVER_DEFAULT, At);
		At_owned = TRUE;
	    }
	}
    }

    phgOptionsPop();

    if (symmetry >= M_SYM) {
	if (St == NULL)
	    St = S;
	At = A;
    }

    if (St != NULL && At == NULL) {
	At = phgSolverGetMat(St);
	At_owned = TRUE;
    }

    if (monitor)
	phgPrintf("*** %s: Solver = %s\n", __func__, S->oem_solver->name);

    v = phgMapCreateVec(S->rhs->map, 1);

    if (monitor)
	phgPrintf("*** %s: Computing largest singular value umax:\n", __func__);
    umax = 0.0;
    phgVecRandomize(v, getpid());
    phgVecAXPBY(0.0, NULL, 1.0 / phgVecNorm2(v, 0, NULL), &v);
    for (it = 1; it <= maxit; it++) {
	/* v := A A^t v */
	if (At != NULL)
	    v0 = phgMatVec(MAT_OP_N, 1.0, At, v, 0.0, NULL);
	else
	    v0 = phgVecCopy(v, NULL);
	phgMatVec(MAT_OP_N, 1.0, A, v0, 0.0, &v);
	phgVecDestroy(&v0);
	d = umax;	/* save old umax */
	umax = phgVecNorm2(v, 0, NULL);
	if (umax != 0.0)
	    phgVecAXPBY(0.0, NULL, 1.0 / umax, &v);
	umax = Sqrt(umax);
	if (monitor)
	    phgPrintf("*** %s:   it %d, umax: %0.4le\n", __func__,
						it, (double)umax);
	if (it  == maxit || Fabs(d - umax) <= tol * umax)
	    break;
    }

    if (monitor)
	phgPrintf("*** %s: Computing smallest singular value umin:\n",__func__);
    cond = 0.0;
    v0 = vt0 = NULL;
    phgVecRandomize(v, getpid());
    phgVecAXPBY(0.0, NULL, 1.0 / phgVecNorm2(v, 0, NULL), &v);
    umin = 1.0;
    for (it = 1; it <= maxit; it++) {
	BOOLEAN handle_bdry_eqns_bak;
	if (St != NULL) {
	    /* v := (A A^t)^(-1) v */
	    St->rhs_updated = TRUE;
	    St->rhs = phgVecCopy(v, NULL);
	    St->rhs->mat = NULL;
	    if (vt0 != NULL) {
		phgVecDestroy(&v);
		v = vt0;
		vt0 = NULL;
	    }
	    else {
		phgVecAXPBY(0.0, NULL, 1. / umin, &v);
	    }
	    handle_bdry_eqns_bak = St->mat->handle_bdry_eqns;
	    St->mat->handle_bdry_eqns = FALSE;
	    phgSolverVecSolve(St, FALSE, v);
	    St->mat->handle_bdry_eqns = handle_bdry_eqns_bak;
	    phgVecDestroy(&St->rhs);
	    if (symmetry < M_SYM)
		vt0 = phgVecCopy(v, NULL);	/* the next initial solution */
	}
	S->rhs_updated = TRUE;
	S->rhs = phgVecCopy(v, NULL);
	S->rhs->mat = NULL;
	if (v0 != NULL) {
	    phgVecDestroy(&v);
	    v = v0;
	    v0 = NULL;
	}
	else {
	    phgVecAXPBY(0.0, NULL, 1. / (St == NULL ? umin * umin : umin), &v);
	}
	handle_bdry_eqns_bak = S->mat->handle_bdry_eqns;
	S->mat->handle_bdry_eqns = FALSE;
	phgSolverVecSolve(S, FALSE, v);
	S->mat->handle_bdry_eqns = handle_bdry_eqns_bak;
	phgVecDestroy(&S->rhs);
	if (St != NULL && symmetry < M_SYM)
	    v0 = phgVecCopy(v, NULL);	/* the next inital solution */
	umin = phgVecNorm2(v, 0, NULL);
	umin = umin == 0. ? 1.0 : 1.0 / umin;
	phgVecAXPBY(0.0, NULL, umin, &v);
	umin = Sqrt(umin);
	d = cond;	/* save old cond */
	cond = umax / umin;	/* separate Sqrt's to avoid FPE */
	if (monitor)
	    phgPrintf("*** %s:   it %d, umin: %0.4le, cond: %0.4le\n",
			__func__, it, (double)umin, (double)cond);
	if (it  == maxit || Fabs(d - cond) <= tol * cond)
	    break;
    }

    phgVecDestroy(&v0);
    phgVecDestroy(&vt0);

    if (St != S && St_bak == NULL)
	phgSolverDestroy(&St);
    if (S_bak == NULL)
	phgSolverDestroy(&S);
    if (At_owned)
	phgMatDestroy(&At);
    if (A_owned)
	phgMatDestroy(&A);
    phgVecDestroy(&v);

    if (S_bak != NULL) {
	CopyParams(S, S_bak);
	phgFree(S_bak);
    }
    if (St_bak != NULL) {
	CopyParams(St, St_bak);
	phgFree(St_bak);
    }

    if (monitor)
	phgPrintf("*** %s: Condition number = %0.4e, wall time = %0.4g\n",
			__func__, (double)cond, phgGetTime(NULL) - t);

    return cond;
}

FLOAT
phgMatNormInfty(MAT *A)
/* returns the infty norm (max(sum(abs(row))) of the matrix */
{
    INT j;
    FLOAT d = 0.0, a;

    for (j = 0; j < A->rmap->nlocal; j++) {
	int ii;
	const MAT_ROW *row = phgMatGetRow(A, j);
	for (a = 0., ii = 0; ii < row->ncols; ii++)
	    a += Fabs(row->data[ii]);
	if (d < a)
	    d = a;
    }

#if USE_MPI
    if (A->rmap->nprocs > 1) {
	MPI_Allreduce(&d, &a, 1, PHG_MPI_FLOAT, PHG_MAX, A->rmap->comm);
	d = a;
    }
#endif	/* USE_MPI */

    return d;
}

MAT *
phgMatDiagonalBlockInverse(MAT *A)
/* computes A^{-1}, assuming A is a block diagonal matrix whose diagonal
 * blocks are all local.
 *
 * returns A^{-1} if successful, NULL if failure.
 * A is assembled and otherwise unchanged. */
{
#ifdef I
# undef I	/* defined in complex.h */
#endif
#ifdef MAXn
# undef MAXn
#endif
#define MAXn 4096	/* maximum allowed block size */
    MAT *Ainv, *ret, *At;
    FLOAT *D, *Dinv;
    MAT_ROW *rows;
    const MAT_ROW *row;
    INT I, J, *list, N, N0, n0;
    int *indx;	/* -2: unprocessed, -1: processed, >=0: current */
    int i, j, k, n, *pvt;


    if (!A->assembled)
	phgMatAssemble(A);

    N = A->rmap->nlocal;
    assert(N == A->cmap->nlocal);
    N0 = A->rmap->partition[A->rmap->rank];

    /* construct A^t (NOT NECESSARY IF A IS SYMMETRIC!!!) */
    At = phgMatCreate(MPI_COMM_SELF, N, N);
    for (I = 0; I < N; I++) {
	row = phgMatGetRow(A, I);
	for (j = 0; j < row->ncols; j++) {
	    J = row->cols[j] - N0;
	    if (J < 0 || J >= N) {
		phgWarning("%s:%d: not locally block diagonal, abort.\n",
				    __FILE__, __LINE__);
		phgMatDestroy(&At);
		return NULL;
	    }
	    phgMatAddEntry(At, J, I, row->data[j]);
	}
    }
    phgMatAssemble(At);

    ret = Ainv = phgMapCreateMat(A->rmap, A->cmap);
    phgMatSetReplaceMode(Ainv);
    rows = phgCalloc(MAXn, sizeof(*rows));
    list = phgAlloc(MAXn * sizeof(*list));

    indx = phgAlloc((N + MAXn) * sizeof(*indx));
    pvt = indx + N;

    for (I = 0; I < N; I++)
	indx[I] = -2;

    for (I = 0; I < N; I++) {
	assert(indx[I] < 0);
	if (indx[I] == -1)
	    continue;	/* row already processed */
	/* find the diagonal block starting at I */
	n = 0;
	phgMatGetRow_(A, I, &rows[n]);
	list[n] = I;
	indx[I] = n++;
	for (i = 0; i < n && ret != NULL; i++) {
	    J = list[i];
	    for (k = 0; k < 2; k++) {	    /* k=0: row, k=1: column */
		row = k == 0 ? &rows[indx[list[i]]] : phgMatGetRow(At, I);
		n0 = k == 0 ? N0 : 0;
		if (row->ncols == 0) {
		    phgWarning("%s: empty row found, abort.\n", __func__);
		    ret = NULL;
		    break;
		}
		/* scan columns in the row */
		for (j = 0; j < row->ncols; j++) {
		    J = row->cols[j] - n0;
		    if (J < 0 || J >= N) {
			phgWarning("%s:%d: not locally block diagonal, abort.\n",
					__FILE__, __LINE__);
			ret = NULL;
			break;
		    }
		    if (indx[J] >= 0)
			continue;
		    assert(indx[J] == -2);	/* unprocessed */
		    if (n >= MAXn) {
			phgWarning("%s: block size exceeds %d.\n", __func__,
									MAXn);
			ret = NULL;
			break;
		    }
		    phgMatGetRow_(A, J, &rows[n]);
		    list[n] = J;
		    indx[J] = n++;
		} /* j */
	    } /* k */
	} /* i */

	/* allocate dense matrix block */
	D = phgCalloc(2 * n * n, sizeof(*D));
	Dinv = D + n * n;
	/* Copy matrix data */
	for (i = 0; i < n; i++) {
	    row = &rows[i];
	    for (j = 0; j < row->ncols; j++)
		D[i * n + indx[row->cols[j] - N0]] = row->data[j];
	    Dinv[i * n + i] = 1.0;
	}
	/* compute inverse matrix */
#if 0
#warning
printf("%s: n = %d\n", __func__, n);
for (i = 0; i < n; i++) {
printf(i == 0 ? "D  " : "   ");
for (j = 0; j < n; j++) printf(" %0.16g", (double)D[i * n + j]);
printf("\n");
}
#endif
	if (phgSolverDenseLU(n, D, pvt) == FALSE) {
	    phgWarning("%s: inreversible diagonal block, abort.\n", __func__);
	    phgFree(D);
	    ret = NULL;
	    break;
	}
	phgSolverDenseSV(n, D, pvt, n, Dinv);
#if 0
#warning
for (i = 0; i < n; i++) {
printf(i == 0 ? "Inv" : "   ");
for (j = 0; j < n; j++) printf(" %0.16g", (double)Dinv[i * n + j]);
printf("\n");
}
#endif
	for (i = 0; i < n; i++)
	    for (j = 0; j < n; j++)
		phgMatSetEntry(Ainv, list[i], list[j], Dinv[i * n + j]);
	phgFree(D);

	for (i = 0; i < n; i++)
	    indx[list[i]] = -1;	/* mark processed entries */
    }

    if (ret == NULL)
	phgMatDestroy(&Ainv);
    else
	phgMatAssemble(Ainv);

    phgFree(list);
    phgFree(indx);
    for (i = 0; i < MAXn; i++) {
	phgFree(rows[i].cols);
	phgFree(rows[i].data);
    }
    phgFree(rows);

#undef MAXn
#undef ADD_ENTRY

    if (At != A)
	phgMatDestroy(&At);

    return ret;
}
