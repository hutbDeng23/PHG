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

/* $Id: eigen-gwg.c,v 1.19 2014/04/04 04:30:02 zlb Exp $
 *
 * This is an interface to Gao Weiguo's LOBPCG package. */

/* Questions:
 *   1. What are the 'overlapped' and 'non overlapped' forms of vectors?
 *   2. What to do with ijob == 4?
 *   3. What are the preconditioning matrices:
 *		PA \approx A^{-1} and PB \approx B^{-1}?
 *		Why PB is not used in the sample code?
 *   4. What is the matrix C in vGERV()?
 *   5. Can we work with partitioned vectors for better scalability?
 *   6. In vDOT: is the MPI_Allreduce() necessary? It seems that the norms of
 *	the eigenvectors are scaled by sqrt(nprocs) (because of it).
 *
 * Suggestions:
 *   1. Provide a header file 'lobpcg_revcomc.h' with prototype for the
 *	function lobpcg_revcomc()
 *   2. Either move the functions orthog_blockab(), orthog_inneraa() and
 *	orthog_innerab() to liblobpcg.a to avoid dependency of lobpcg_revcomc()
 *	on util.o, or add util.o to liblobpcg.a
 *   3. Fix the annoying messages '... parameter XXX had an illegal value'
 */

#include "phg.h"
#include "phg/blas-lapack.h"
#include <string.h>
#include <math.h>

#if USE_GWLOBPCG

#if 0
#include "phg/lobpcg_revcomc.h"
#else
extern int lobpcg_revcomc(int myid,MPI_Comm comm,int nprocs,int n,int nev,
    double *X,double *AX,double *BX,
    double *W,double *AW,double *BW,
    double *P,double *AP,double *BP,
    double *R,double *SA,double *SB,
    double *lambda,double *nrm,
    double *work,double *Xini,
    int lwork,int maxit,int ng,int *idx,
    int *iter,int *conv,int *nnc,
    int *ijob,int *inxt,int *iidx);
#endif

/*------------------------ Interface functions ------------------------*/

/* 1: LOBPCG I, 0: LOBPCG II, >1: LOBPCG with gap detection */
static int		lobpcg_type = 0;
static const char	*lobpcg_type_keywords[] = {"I", "II", "gap", NULL};
static int		lobpcg_type_values[] = {1, 0, 2};

static MAT *mA, *mB;
static VEC *vx, *vy;
static int rank, nprocs;
static MPI_Comm comm;
#if USE_MPI
static int *dsps = NULL, *cnts = NULL;
#endif	/* USE_MPI */

/* function in lobpcg_mv.c */

static void
vGEMV(char *mat, double *x, double *y)
{
    int i;

    if (!strcmp(mat, "A")) {
	if (nprocs == 1) {
	    vx->data = x;
	    vy->data = y;
	    phgMatVec(MAT_OP_N, 1.0, mA, vx, 0.0, &vy);
	    vx->data = vy->data = NULL;
	    return;
	}
#if USE_MPI
	vx->data = x + dsps[rank];
	phgMatVec(MAT_OP_N, 1.0, mA, vx, 0.0, &vy);
	MPI_Allgatherv(vy->data, mA->cmap->nlocal, PHG_MPI_FLOAT,
		       y, cnts, dsps, PHG_MPI_FLOAT, comm);
	vx->data = NULL;
	return;
#endif	/* USE_MPI */
    }
    else if (!strcmp(mat, "B")) {
	if (nprocs == 1) {
	    vx->data = x;
	    vy->data = y;
	    phgMatVec(MAT_OP_N, 1.0, mB, vx, 0.0, &vy);
	    vx->data = vy->data = NULL;
	    return;
	}
#if USE_MPI
	vx->data = x + dsps[rank];
	phgMatVec(MAT_OP_N, 1.0, mB, vx, 0.0, &vy);
	MPI_Allgatherv(vy->data, mB->cmap->nlocal, PHG_MPI_FLOAT,
		       y, cnts, dsps, PHG_MPI_FLOAT, comm);
	vx->data = NULL;
	return;
#endif	/* USE_MPI */
    }
    else if (!strcmp(mat, "PA")) {
	if (mA->diag == NULL)
	    phgMatSetupDiagonal(mA);
	if (nprocs == 1) {
	    for (i = 0; i < mA->cmap->nlocal; i++)
		y[i] = x[i] / mA->diag[i];
	    return;
	}
#if USE_MPI
	for (i = 0; i < mA->cmap->nlocal; i++)
	    vy->data[i] = x[dsps[rank] + i] / mA->diag[i];
	MPI_Allgatherv(vy->data, mA->cmap->nlocal, PHG_MPI_FLOAT,
		       y, cnts, dsps, PHG_MPI_FLOAT, comm);
	return;
#endif	/* USE_MPI */
    }

    phgError(1, "unimplemented.\n");
}

static void
vGERV(int rank, MPI_Comm comm, double *x, double *y, double *work)
{
    int i;

    for (i = 0; i < mA->cmap->nglobal; i++)
	y[i] += x[i] * 0.;
}

/* function in lobpcg_util.c */

static double
vDOT(int rank, MPI_Comm comm, int n, const double *x1, const double *x2)
{
    double res = 0.;
    int i;

    for (i = 0; i < n; i++)
	res += x1[i] * x2[i];

#if USE_MPI
    /* FIXME: the norms of the eigenvectors are divided by sqrt(nprocs) */
    if (nprocs > 1) {
	double tmp;
	MPI_Allreduce(&res, &tmp, 1, MPI_DOUBLE, MPI_SUM, comm);
	res = tmp;
    }
#endif	/* USE_MPI */

    return res;
}

#define vMix(rank, comm, x)	/* do nothing */
#define dcopy_(n, src, incs, dst, incd) {	\
    assert(*(incd) == 1 && *(incs) == 1);	\
    memcpy(dst, src, *(n) * sizeof(double));	\
}

/*---------------------------------------------------------------------------*/
/* The following functions are copied from lobpcg_util.c.
 * They're required by lobpcg_revcomc(), are currently missing from liblobpcg.a
 * and are to be removed later */

/* the next line disables ScaLAPACK's error messages */
void F77_FUNC(pxerbla,PXERBLA)(void) {return;}

int
orthog_blockab(int myid,MPI_Comm comm,
    int n,int a,double *X,double *AX,double *BX,int b,
    double *P,double *AP,double *BP,double *work)
/* make the block orthog for P about X */
{
    static double one=1.,zero=0.,mone=-1.;
    int lm=a,ll=a*b;
	/* work = X'*BP */
    F77_FUNC(dgemm,DGEMM)("T", "N", &a, &b, &n, &one, X, &n, BP, &n, &zero,
			  work+ll, &lm, 1, 1);
    MPI_Allreduce(work+ll, work, ll, MPI_DOUBLE, MPI_SUM, comm);

    /* P = P - X * work */
    F77_FUNC(dgemm,DGEMM)("N", "N", &n, &b, &a, &mone, X, &n, work, &lm, &one,
			  P, &n, 1, 1);
    /* AP = AP - AX * work */
    F77_FUNC(dgemm,DGEMM)("N", "N", &n, &b, &a, &mone, AX, &n, work, &lm, &one,
			  AP, &n, 1, 1);
	/* BP = BP - BX * work */
    F77_FUNC(dgemm,DGEMM)("N", "N", &n, &b, &a, &mone, BX, &n, work, &lm, &one,
			  BP, &n, 1, 1);

    return 0;
}

int
orthog_innerab(int myid,MPI_Comm comm,int n,
    int a,double *X,double *AX,double *BX,int b,
    double *P,double *AP,double *BP)
/* make MGS orthog for P about X */
{
    int i,j;
    static int ione=1;
    double dot,r;
    double *X_j,*AX_j,*BX_j,*P_i,*AP_i,*BP_i;

    for (i=0;i<b;i++) {
	P_i=P+i*n;
	AP_i=AP+i*n;
	BP_i=BP+i*n;
	dot=vDOT(myid,comm,n,P_i,BP_i);
	for (j=0;j<a;j++) {
		X_j=X+j*n;
		AX_j=AX+j*n;
		BX_j=BX+j*n;
		r=-vDOT(myid,comm,n,P_i,BX_j);
		F77_FUNC(daxpy,DAXPY)(&n,&r,X_j,&ione,P_i,&ione);
		F77_FUNC(daxpy,DAXPY)(&n,&r,AX_j,&ione,AP_i,&ione);
		F77_FUNC(daxpy,DAXPY)(&n,&r,BX_j,&ione,BP_i,&ione);
	}
    r=vDOT(myid,comm,n,P_i,BP_i);

	if (r<0.5*dot) {
		dot=r;
		for (j=0;j<a;j++) {
			X_j=X+j*n;
			AX_j=AX+j*n;
			BX_j=BX+j*n;
			r=-vDOT(myid,comm,n,P_i,BX_j);
			F77_FUNC(daxpy,DAXPY)(&n,&r,X_j,&ione,P_i,&ione);
			F77_FUNC(daxpy,DAXPY)(&n,&r,AX_j,&ione,AP_i,&ione);
			F77_FUNC(daxpy,DAXPY)(&n,&r,BX_j,&ione,BP_i,&ione);
		}
		r=vDOT(myid,comm,n,P_i,BP_i);
	}
	r=1.0/sqrt(r);
	F77_FUNC(dscal,DSCAL)(&n,&r,P_i,&ione);
	F77_FUNC(dscal,DSCAL)(&n,&r,AP_i,&ione);
	F77_FUNC(dscal,DSCAL)(&n,&r,BP_i,&ione);
    }
    return 0;
}

/* make MGS orthog inner P*/
int
orthog_inneraa(int myid,MPI_Comm comm,int n,
    int a,double *X,double *AX,double *BX)
{
    int i,j;
    static int ione=1;
    double dot,r;
    double *X_i,*AX_i,*BX_i,*X_j,*AX_j,*BX_j;

    for (i=0;i<a;i++) {
	X_i=X+i*n;
	AX_i=AX+i*n;
	BX_i=BX+i*n;
	dot=vDOT(myid,comm,n,X_i,BX_i);
	for (j=0;j<i;j++) {
	    X_j=X+j*n;
	    AX_j=AX+j*n;
	    BX_j=BX+j*n;
	    r=-vDOT(myid,comm,n,X_i,BX_j);
	    F77_FUNC(daxpy,DAXPY)(&n,&r,X_j,&ione,X_i,&ione);
	    F77_FUNC(daxpy,DAXPY)(&n,&r,AX_j,&ione,AX_i,&ione);
	    F77_FUNC(daxpy,DAXPY)(&n,&r,BX_j,&ione,BX_i,&ione);
	}
	r=vDOT(myid,comm,n,X_i,BX_i);
	if (r<0.5*dot) {
	    dot=r;
	    for (j=0;j<i;j++) {
		X_j=X+j*n;
		AX_j=AX+j*n;
		BX_j=BX+j*n;
        r=-vDOT(myid,comm,n,X_i,BX_j);
		F77_FUNC(daxpy,DAXPY)(&n,&r,X_j,&ione,X_i,&ione);
		F77_FUNC(daxpy,DAXPY)(&n,&r,AX_j,&ione,AX_i,&ione);
		F77_FUNC(daxpy,DAXPY)(&n,&r,BX_j,&ione,BX_i,&ione);
	    }
        r=vDOT(myid,comm,n,X_i,BX_i);
	}
	r=1.0/sqrt(r);
	F77_FUNC(dscal,DSCAL)(&n,&r,X_i,&ione);
	F77_FUNC(dscal,DSCAL)(&n,&r,AX_i,&ione);
	F77_FUNC(dscal,DSCAL)(&n,&r,BX_i,&ione);
    }

    return 0;
}

/*---------------------------------------------------------------------------*/

EigenSolverPrototype(phgEigenSolverGWLOBPCG)
/* (A, B, n, which, tau, evals, evecs, eps, itmax, nit) */
{
    static BOOLEAN initialized = FALSE;
    VEC *vec;
 
    int i,j,nnc;  
    int lwork, conv, iter;
    int ione=1,info;
    int ijob,inxt,iidx,*idx;
    double *X,*AX,*BX,*W,*AW,*BW,*P,*AP,*BP,*R,*SA,*SB,*lambda,*nrm,*work,*Xini;
    double maxr,minr;
    double dot,dotax,dotbx;
#if 0
    int ln,lc,lnc,mm,len,k;
    int *gidx,*oidx,*midx, flag, start;    
    double gap,tol2=1.e-5,*sa,*sb;
    double one=1.,zero=0.;
    static double maxeig=0.;
#endif

    Unused(ione);

    if (!initialized) {
	/* register options */
	initialized = TRUE;
	phgOptionsRegisterTitle("\nGWLOBPCG options:", "\n", "gwlobpcg");
	phgOptionsRegisterKeyword("-gwlobpcg_type",
			"LOBPCG type ('gap' means LOBPCG with gap detection)",
			lobpcg_type_keywords, &lobpcg_type);
	return 0;
    }

    assert(FT_PHG == FT_DOUBLE);
    assert(sizeof(BLAS_INT) == sizeof(int));
    assert(which < 0 || (which == 0 && tau == 0.0));

    *nit = 0;
    vec = *evecs;
    MagicCheck(VEC, vec);

    /* Note: n is number of eigenvalues to compute (nev) */
    if (n <= 0 || A->cmap->nglobal <= 0)
	return 0;

    if (n >= A->cmap->nglobal)
	n = A->cmap->nglobal;

    rank = A->cmap->rank;
    nprocs = A->cmap->nprocs;
    comm = A->cmap->comm;

    mA = A;
    mB = B;
    vx = phgMapCreateVec(A->cmap, 1);
    vy = phgMapCreateVec(A->cmap, 1);
    phgFree(vx->data);
    vx->data = NULL;
    if (nprocs <= 1) {
	phgFree(vy->data);
	vy->data = NULL;
    }
#if USE_MPI
    else {
	dsps = phgAlloc(2 * nprocs * sizeof(*dsps));
	cnts = dsps + nprocs;
	for (i = 0; i < nprocs; i++) {
	    dsps[i] = A->cmap->partition[i];
	    cnts[i] = A->cmap->partition[i + 1] - dsps[i];
	}
    }
#endif	/* USE_MPI */

#define LA	24
#define N	(A->cmap->nglobal)
#define nev	n
    /* allocate memeory */
    lwork = (18*nev*nev+nev > 2*LA ? 18*nev*nev+nev : 2*LA);
    lwork = (lwork>N ? lwork : N);

    X = calloc(10*N*nev+18*nev*nev+10*nev+lwork + N*nev, sizeof(double));
    Xini = X + 10*N*nev+18*nev*nev+10*nev+lwork;
    idx = calloc(3*nev+9, sizeof(int));

    if (vec == NULL) {
	*evecs = vec = phgMapCreateVec(A->cmap, nev);
    }
    else {
	assert(vec->nvec == n);
	if (nprocs <= 1) {
	    memcpy(Xini, vec->data, nev * A->cmap->nlocal * sizeof(*Xini));
	}
#if USE_MPI
	else {
	    MPI_Datatype type;
	    double *p, *q;
	    MPI_Type_vector(1, nev, nev, PHG_MPI_FLOAT, &type);
	    MPI_Type_commit(&type);
	    MPI_Allgatherv(vec->data, A->cmap->nlocal, type,
			   X, cnts, dsps, type, comm);
	    MPI_Type_free(&type);
	    for (i = 0, p = Xini; i < nev; i++) {
		q = X + i * cnts[0];
		for (j= 0; j < nprocs; j++) {
		    memcpy(p, q, cnts[j] * sizeof(*X));
		    p += cnts[j];
		    if (j < nprocs - 1)
			q += (nev - i) * cnts[j] + i * cnts[j + 1];
		}
	    }
	    memset(X, 0, nev * N * sizeof(*X));
	}
#endif	/* USE_MPI */
    }

    AX=X+3*N*nev;  
    BX=X+6*N*nev;  
    R=X+9*N*nev;
    SA=X+10*N*nev; 
    SB=SA+9*nev*nev+3*nev;
    lambda=SB+9*nev*nev+3*nev;
    nrm=lambda+3*nev;
    work=nrm+nev;

    W=X+N*nev;     P=X+2*N*nev;
    AW=AX+N*nev;   AP=AX+2*N*nev;
    BW=BX+N*nev;   BP=BX+2*N*nev;

    iidx=0;
    inxt=0;
    conv=0;
    nnc=nev;
    while (1) {
        info = lobpcg_revcomc(rank, comm, nprocs,N,nev,
                X,AX,BX,
                W,AW,BW,
                P,AP,BP,
                R,SA,SB,
                lambda,nrm,
                work,Xini,
                lwork,itmax,lobpcg_type_values[lobpcg_type],idx,
                &iter,&conv,&nnc,
                &ijob,&inxt,&iidx);
#if 0
	phgPrintf("revcomc: info=%d, ijob=%d, iter=%d, conv=%d, nnc=%d\n",
			info, ijob, iter, conv, nnc);
#endif
                
        if (info==1) break;

        switch (ijob) {
            case 0:
                /* AX=A*X;    BX=B*X */
                for (j=conv;j<nev;j++) {
                    vGEMV("A",&X[j*N],&AX[j*N]);
                    vGERV(rank,comm,&X[j*N],&AX[j*N],work);
                    vGEMV("B",&X[j*N],&BX[j*N]);
                }
                break;
            case 1:
                /* AW=A*W;    BW=B*W */
                for (j=conv;j<nev;j++) {
                    vGEMV("A",&W[j*N],&AW[j*N]);
                    vGERV(rank,comm,&W[j*N],&AW[j*N],work);
                    vGEMV("B",&W[j*N],&BW[j*N]);
                }
                break;
            case 2:
                /* do precondition */
                for (j=conv;j<nev;j++) {
                    vGEMV("PA",&R[j*N],&W[j*N]);
                }
                break;
            case 3:
                /* get residual norm */
                maxr=0.; minr=1.;
                for (j=conv;j<nev;j++) {
                    dcopy_(&N,&R[j*N],&ione,work,&ione);
                    vMix(rank,comm,&R[j*N]);		
                    dotax=sqrt(vDOT(rank,comm,N,&AX[j*N],&AX[j*N]));
                    dotbx=sqrt(vDOT(rank,comm,N,&BX[j*N],&BX[j*N]));
                    dot=sqrt(vDOT(rank,comm,N,&R[j*N],work));
                    nrm[j]=dot/(dotax+fabs(lambda[j])*dotbx);
                    if (rank==0) {
                        if (minr>nrm[j]) { /* minr: min error */
                            minr=nrm[j];
            	        } 
            	        if (maxr<nrm[j]) { /* maxr: max error */
                            maxr=nrm[j];
            	        }
            	        if (nrm[j]<eps && conv==j) {
                            conv++;
            	        }
                    }
                }
                
                MPI_Bcast(&conv,1,MPI_INT,0,comm);
                if (phgVerbosity > 0 && rank == 0 && conv < nev)
                    printf("iter=%d, conv=%d, nnc=%d, nrmc=%e, minr=%e, "
			   "maxr=%e\n", iter, conv, nnc, nrm[conv], minr, maxr);
                break;
	    case 4:
		/* FIXME: what to do? */
                break;
            default:
		phgError(1, "unexpected 'ijob=%d'.\n", ijob);
		break;
        }
    }

    *nit = iter;
    n = conv;

    Xini = X;
    for (i = 0; i < rank; i++)
	Xini += cnts[i];
    for (i = 0; i < nev; i++) {
	evals[i] = lambda[i];
	for (j = 0; j < A->cmap->nlocal; j++)
	    vec->data[i * A->cmap->nlocal + j] = Xini[i * N + j];
    }

    free(X);
    free(idx);

#if USE_MPI
    free(dsps);
    dsps = cnts = NULL;
#endif	/* USE_MPI */
    phgVecDestroy(&vx);
    phgVecDestroy(&vy);

    return n;
}

#else	/* USE_GWLOBPCG */
EigenSolverVoid(phgEigenSolverGWLOBPCG)
#endif	/* USE_GWLOBPCG */
