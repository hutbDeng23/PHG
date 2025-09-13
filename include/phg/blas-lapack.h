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

/* $Id: blas-lapack.h,v 1.1 2014/04/04 04:30:00 zlb Exp $ */

#ifndef PHG_BLAS_LAPACK_H

#ifdef __cplusplus
extern "C" {
#endif

#if USE_BLAS
/*---------------------------------- BLAS -----------------------------------*/
void F77_FUNC(dgemm,DGEMM)(const char *, const char *,
			   BLAS_INT *, BLAS_INT *, BLAS_INT *,
			   double *, double *, BLAS_INT *, double *,
			   BLAS_INT *, double *, double *, BLAS_INT *,
			   int, int);
void F77_FUNC(sgemm,SGEMM)(const char *, const char *,
			   BLAS_INT *, BLAS_INT *, BLAS_INT *,
			   float *, float *, BLAS_INT *, float *,
			   BLAS_INT *, float *, float *, BLAS_INT *,
			   int, int);

void F77_FUNC(daxpy,DAXPY)(BLAS_INT *, double *, double *, BLAS_INT *,
			   double *, BLAS_INT *);
void F77_FUNC(saxpy,SAXPY)(BLAS_INT *, float *, float *, BLAS_INT *,
			   float *, BLAS_INT *);

void F77_FUNC(dscal,DSCAL)(BLAS_INT *, double *, double *, BLAS_INT *);
void F77_FUNC(sscal,SSCAL)(BLAS_INT *, float *, float *, BLAS_INT *);
#endif	/* USE_BLAS */

#if USE_LAPACK
/*--------------------------------- LAPACK ----------------------------------*/
void F77_FUNC(dgetrf,DGETRF)(BLAS_INT *, BLAS_INT *, double *,
			     BLAS_INT *, BLAS_INT *, BLAS_INT *);
void F77_FUNC(sgetrf,SGETRF)(BLAS_INT *, BLAS_INT *, float *,
			     BLAS_INT *, BLAS_INT *, BLAS_INT *);

void F77_FUNC(dgetrs,DGETRS)(char *, BLAS_INT *, BLAS_INT *,
			     double *, BLAS_INT *, BLAS_INT *,
			     double *, BLAS_INT *, BLAS_INT *, int);
void F77_FUNC(sgetrs,SGETRS)(char *, BLAS_INT *, BLAS_INT *,
			     float *, BLAS_INT *, BLAS_INT *,
			     float *, BLAS_INT *, BLAS_INT *, int);
#endif	/* USE_LAPACK */

#ifdef __cplusplus
}
#endif

#define PHG_BLAS_LAPACK_H
#endif
