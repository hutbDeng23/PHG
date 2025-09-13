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

/* $Id: eigen.h,v 1.1 2014/04/04 04:30:00 zlb Exp $ */

#ifndef PHG_EIGEN_H
#define PHG_EIGEN_H

#ifdef __cplusplus
extern "C" {
#endif

/* phgEigenSolve tries to compute 'n' generalized eigenvalues 'lambda' and
 * corresponding eigenvectors 'x' satisfying A*x = lambda*B*x.
 *
 * The return value is the number of eigenvalues successfully found.
 *
 * 'which' < 0 means compute the 'n' smallest eigenvalues,
 * 'which' > 0 means compute the 'n' largest eigenvalues,
 * 'which' == 0 means compute the 'n' eigenvalues closest to 'tau'
 * 'B' == NULL means 'B' == I (standard eigenvalues)
 *
 * If at input 'eigenvectors' != NULL, then it contains initial guesses
 * for the eigenvectors.
 */

/* Note: the last argument, ssd, points to solver specific data, which can be
 * set to NULL and whose meaning is interpreted by specific solvers. */
#define EigenSolverPrototype(func) int \
    func(MAT *A, MAT *B, int n, int which, FLOAT tau, FLOAT *evals, \
	 VEC **evecs, FLOAT eps, int itmax, int *nit, void *ssd)

#define EigenSolverVoid(func) EigenSolverPrototype(func) {	\
    Unused(A);							\
    Unused(B);							\
    Unused(n);							\
    Unused(which);						\
    Unused(tau);						\
    Unused(evals);						\
    Unused(evecs);						\
    Unused(eps);						\
    Unused(itmax);						\
    Unused(nit);						\
    Unused(ssd);						\
    return 0;							\
}

typedef EigenSolverPrototype((*EIGEN_SOLVER));

enum {EIGEN_SMALLEST = -1, EIGEN_CLOSEST = 0, EIGEN_LARGEST = 1};

extern EigenSolverPrototype(phgEigenSolverJDBSYM);
extern EigenSolverPrototype(phgEigenSolverBLOPEX);
extern EigenSolverPrototype(phgEigenSolverPRIMME);
extern EigenSolverPrototype(phgEigenSolverARPACK);
extern EigenSolverPrototype(phgEigenSolverSLEPc);
extern EigenSolverPrototype(phgEigenSolverTrilinos);
extern EigenSolverPrototype(phgEigenSolverGWLOBPCG);

int phgEigenSolve_(MAT *A, MAT *B, int n, int which, FLOAT tau,
		   FLOAT *evals, VEC **evecs, int *nit, void *ssd);
#define phgEigenSolve(A, B, n, which, tau, evals, evecs, nit) \
	phgEigenSolve_(A, B, n, which, tau, evals, evecs, nit, NULL)

#ifdef __cplusplus
}
#endif

#endif	/* !defined(PHG_EIGEN_H) */
