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

/* $Id: sturm.c,v 1.51 2022/04/10 02:47:46 zlb Exp $ */

/* Find real roots of a polynamial in a given interval using Sturm method */

/* Define the following to use BRENT's algorithm to compute isoloated roots.
   (it's only faster when high precision is requested (<1e-12)) */
#undef BRENT

/* Don't enable NEWTON before it's complete!!! */
#ifndef USE_NEWTON
# define USE_NEWTON	1	/* use Newton iterations for accel./improve */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "phg.h"

#include "phg/sturm.h"

#ifdef BRENT
#  include "phg/brent.h"
#endif /*  */

/* Lower the precision by 8 bits to avoid some troubles caused by
   round-off errors, especially when deciding if a and b are roots */

static FLOAT *sturm_rwk;	/* the polynomials in the Sturm series */
static int sturm_np;		/* number of polynomials in the Sturm series */
static int *sturm_iwk;		/* degrees of the polynomials */
static FLOAT *sturm_r;		/* contains the roots found */
static int sturm_nr;		/* contains the number of roots found */
static FLOAT sturm_tol;	/* tolerance */
static int  sturm1(FLOAT *p, int n, FLOAT a, FLOAT b, FLOAT tol,
		       FLOAT *rwk, int *iwk);
static void  sturm0(int ia, FLOAT a, FLOAT pa, int ib, FLOAT b, FLOAT pb);
static void  sturm_factor(FLOAT *p, int n, FLOAT a, FLOAT *rwk);
static int  sturm_init(FLOAT *p, int n, FLOAT tol, FLOAT *rwk, int *iwk);
static int  sturm_remainder(int n0, FLOAT *p0, int n1, FLOAT *p1, FLOAT *p2);
static int  sturm_ncis(FLOAT x, FLOAT *px);

#if USE_NEWTON
#include <string.h>
static int sturm_newton(int n, FLOAT *p, FLOAT *p1, FLOAT *x,
					FLOAT a, FLOAT b, int maxits);
#endif

static FLOAT EPS;	/* the 'machine' precision */
static BOOLEAN FLAG;	/* TRUE ==> only check existence of roots (tol < 0) */

#if USE_OMP
# pragma omp threadprivate(sturm_rwk, sturm_np, sturm_iwk, sturm_r, sturm_nr, \
			   sturm_tol, EPS, FLAG)
#endif  /* USE_OMP */

/*static*/ void
poly_print_symbolic(FILE *fout, const char *prefix, int n, FLOAT *p)
/* prints symbolic form of the polynomial to fout */
{
    int i, flag = 0;
    char s[48];

    fprintf(fout, "%s", prefix);
    for (i = 0; i <= n; i++) {
	if (p[i] == 0.)
	    continue;
	fprintf(fout, "%s", p[i] > 0. && flag ? "+" :
			    p[i] == -1 && i != n ? "-" : "");
#if FT_PHG == FT___FLOAT128
	quadmath_snprintf(s, sizeof(s), "%0.34Qg", p[i]);
#else
	sprintf(s, "%0.16lg", (double)p[i]);
#endif
	if (i == n)
	    fprintf(fout, "%s", s);
	else if (Fabs(p[i]) != 1.0)
	    fprintf(fout, "%s*", s);
	if (n - i > 0) {
	    fprintf(fout, "x");
	    /* Note: '**' works with both gnuplot and sagemath */
	    if (n - i > 1)
		fprintf(fout, "**%d", n - i);
	}
	flag = 1;
    }
    fprintf(fout, "%s\n", flag ? "" : "0");
}

static void
poly_print(FILE *fout, const char *prefix, int n, FLOAT *p)
/* prints coefficients (in decreasing power order) of the polynomial */
{
    int i;

    fprintf(fout, "%s%d:", prefix, n);
    for (i = 0; i <= n; i++) {
#if FT_PHG == FT___FLOAT128
	char s[48];
	quadmath_snprintf(s, sizeof(s), "%0.34Qg", p[i]);
	fprintf(fout, " %s", s);
#else	/* FT_PHG == FT___FLOAT128 */
	fprintf(fout, " %0.16lg", (double)p[i]);
#endif	/* FT_PHG == FT___FLOAT128 */
    }
    fprintf(fout, "\n");
}

int
phgQuadSturm(FLOAT *p, int n, FLOAT a, FLOAT b, FLOAT tol, FLOAT eps,
	     FLOAT *rwk, int *iwk)
/*-------------------------------------------------------------------------
 * p  	: contains coefficients of the polynomial as input and the roots
 *	  found as output, p(x) = p[0] * x^n + p[1] * x^(n-1) + ... + p[n]
 * n	: order of the polynomial.
 * a,b	: search for roots in the interval [a,b].
 * tol  : stopping criterion for the bisections and iterations.
 *	  If tol < 0 then this function does not actually computes the roots,
 *	  it returns !0 if the polynomial has roots in [a,b], and 0 otherwise.
 * eps  : the 'machine' precision (can be set to 0 to use FLOAT_EPSILON)
 * rwk	: work vector of length >= ((n+1)*(n+2)+4)/2.
 * 	  It can be NULL, in which case it will be dynamically allocated.
 * iwk	: work vector of length >= n+1, also returns multiplicity of the roots.
 *	  It can be NULL, in which case it will be dynamically allocated.
 *
 * Returns number of different roots found (-1 means a null polynomial).
 *-------------------------------------------------------------------------*/
{
    FLOAT pa, pb, *p0;
    int ia, ib, k, nr = 0, *iwk0 = NULL;
    FLOAT *rwk0 = NULL;

    Unused(poly_print);
    Unused(poly_print_symbolic);

#if 0
#warning
    printf("%s: a = %g, b = %g, ", __func__, (double)a, (double)b);
    poly_print(stdout, "order = ", n, p);
#endif /*  */

    if (tol < 0.) {
	tol = FLOAT_EPSILON;
	FLAG = TRUE;
    }
    else {
	FLAG = FALSE;
    }

    EPS = eps;
    if (EPS < FLOAT_EPSILON)
	EPS = FLOAT_EPSILON * 1000.0;

    p0 = p;
    if (a > b) {
	pa = a;
	a = b;
	b = pa;
    }

    /* normalize the coefficients */
    pa = Fabs(p[0]);
    for (ia = 1; ia <= n; ia++) {
	pb = Fabs(p[ia]);
	if (pa < pb)
	    pa = pb;
    }
    if (pa == 0.0)
	return (-1);
    if (pa != 1.0) {
	pa = 1.0 / pa;
	for (ia = 0; ia <= n; ia++)
	    p[ia] *= pa;
    }

    if (rwk == NULL)
	rwk0 = rwk = malloc(((n + 1) * (n + 2) + 4) / 2 * sizeof(*rwk));    

    if (iwk == NULL)
	iwk0 = iwk = malloc((n + 1) * sizeof(*iwk));    

    /* remove leading zeros */
    ia = 0;
    while (ia <= n && Fabs(p[ia]) <= EPS)
	ia++;
#if 0
    if (ia > n)
	goto end;
#else
    assert(ia <= n);	/* p[] are normalized such that max{|p[]|}==1.0 */
#endif
    if (ia)
	for (ib = ia; ib <= n; ib++)
	    p[ib - ia] = p[ib];
    n -= ia;

    /* count multiplicity of '0' as a root */
    ia = n;
    while (ia >= 0 && Fabs(p[ia]) <= EPS)
	ia--;
    if (ia < n && a <= 0. && b >= 0.) {
	if (FLAG)
	    return 1;
	iwk[0] = n - ia;
	for (ib = ia; ib >= 0; ib--)
	    p[ib + 1] = p[ib];
	p[0] = 0.0;
	p++;
	iwk++;
	nr++;
    }
    n = ia;

    /* check if a is root of p[x] */
    k = 0;
    while (n) {
	sturm_factor(p, n, a, rwk);
	pa = 0.0;
	for (ia = 0; ia < n; ia++) {
	    if (pa < Fabs(rwk[ia]))
		pa = Fabs(rwk[ia]);
	}
	if (Fabs(rwk[n] / pa) > EPS)
	    break;
	if (FLAG)
	    return 1;
	for (ia = 0; ia < n; ia++)
	    p[ia] = rwk[ia];
	k++;
	n--;
    }
    if (k) {
	for (ia = n + 1; ia > 0; ia--)
	    p[ia] = p[ia - 1];
	p[0] = a;
	iwk[0] = k;
	p++;
	iwk++;
	nr++;
    }

    /* check if b is root of p[x] */
    k = 0;
    while (n && a != b) {
	sturm_factor(p, n, b, rwk);
	pa = 0.0;
	for (ia = 0; ia < n; ia++) {
	    if (pa < Fabs(rwk[ia]))
		pa = Fabs(rwk[ia]);
	}
	if (Fabs(rwk[n] / pa) > EPS)
	    break;
	if (FLAG)
	    return 1;
	for (ia = 0; ia < n; ia++)
	    p[ia] = rwk[ia];
	k++;
	n--;
    }
    if (k) {
	for (ia = n + 1; ia > 0; ia--)
	    p[ia] = p[ia - 1];
	p[0] = b;
	iwk[0] = k;
	p++;
	iwk++;
	nr++;
    }

    if (n && a != b) {
#if 0
#warning
	poly_print(stdout, "after pre-precessing: n=", n, p);
#endif
#if USE_NEWTON
	/* save coefficients for post-processing with Newton iterations */
	FLOAT *p_save = NULL;
	if (n > 2) {
	    p_save = malloc((n + 1) * sizeof(*p_save));
	    memcpy(p_save, p, (n + 1) * sizeof(*p_save));
	}
#endif	/* USE_NEWTON */
	ia = sturm1(p, n, a - EPS, b + EPS, tol, rwk, iwk);
	assert(ia >=0 );
#if USE_NEWTON
	if (n > 2) {
	    /* perform Newton iterations to improve the precision */
	    int *ploc;	/* list of locations for the derivatives of p(x) */
	    memcpy(rwk, p_save, (n + 1) * sizeof(*p_save));
	    free(p_save);
	    ploc = malloc((n + 1) * sizeof(*ploc));
	    k = 0;	/* k is the highest order of computed derivatives */
	    ploc[0] = 0;
	    for (ib = 0; ib < ia; ib++) {
		int i, m;
		m = iwk[ib];	/* the multiplicity */
		while (k < m) {
		    FLOAT *p, *p1;
		    /* compute derivatives of p(x) up to order 'm':
		     * 	rwk[ploc[k]] = p^{k}(x) */
		    ploc[k + 1] = ploc[k] + (n - k + 1);
		    p = rwk + ploc[k];
		    p1 = rwk + ploc[k + 1];
		    for (i = 0; i < n - k; i++)
			p1[i] = p[i] * (n - k - i);
		    k++;
		}
		/* perform Newton iterations with p^{m-1}(x) */
		sturm_newton(n - (m - 1),	/* order of p^{m-1} */
			     rwk + ploc[m - 1],	/* p^{m-1} */
			     rwk + ploc[m],	/* p^{m} */
			     p + ib, a, b, 3);
	    }
	    free(ploc);
	}
#endif	/* USE_NEWTON */
	nr += ia;
    }

    if (!FLAG) {
	for (ia = 0; ia < nr; ia++) {
	    if (p0[ia] < a && p0[ia] >= a - EPS)
		p0[ia] = a;
	    else if (p0[ia] > b && p0[ia] <= b + EPS)
		p0[ia] = b;
	}
#if 0
#warning
	poly_print(stdout, "nr=", nr, p0);
#endif /*  */
    }

    if (rwk0 != NULL)
	free(rwk0);

    if (iwk0 != NULL)
	free(iwk0);

    return nr;
}

static int
check_range(int n, FLOAT a, FLOAT b, FLOAT p[], int iwk[])
/* checks the range of the roots in p[]/iwk[] (sorted) and remove those
 * outside of [a,b], returns the number of roots in [a,b] */
{
    int i, i0;

    while (n > 0 && p[n - 1] > b)
	n--;
    if (n == 0)
	return 0;

    for (i0 = 0; i0 < n; i0++)
	if (p[i0] >= a)
	    break;
    if (i0 >= n)
	return 0;
    else if (i0 == 0)
	return n;

    for (i = 0; i0 < n; i0++, i++) {
	p[i] = p[i0];
	iwk[i] = iwk[i0];
    }
    return i;
}

static int 
sturm1(FLOAT *p, int n, FLOAT a, FLOAT b, FLOAT tol, FLOAT *rwk, int *iwk)
/* computes the roots in [a,b] of the given polynomial and returns # roots. */
{
    FLOAT pa, pb;
    int ia, ib, nr;
    static int recursion = 0;
#if USE_OMP
# pragma omp threadprivate(recursion)
#endif  /* USE_OMP */

#if 0
#warning
    printf("%s: recursion = %d, n = %d\n", __func__, recursion, n);
#endif

    nr = 0;
    if (!n)
	return (0);
    else if (n == 1) {		/* linear equation */
	p[0] = -p[1] / p[0];
	iwk[0] = 1;
      check:
	return (((p[0] >= a && p[0] <= b) ? 1 : 0));
    }
    else if (n == 2) {		/* quadratic equation */
	pa = p[1] * p[1] - 4.0 * p[0] * p[2];		/* pa = Delta^2 */
	if (pa < -EPS)		/* no real roots */
	    return (0);
	if (pa <= EPS * EPS) {	/* two equal roots */
	    p[0] = -p[1] / (p[0] + p[0]);
	    iwk[0] = 2;
	    goto check;
	}
	else {			/* two different roots */
	    pa = Sqrt(pa);				/* pa = Delta */
	    if (p[1] >= 0.0) {
		pb = -(p[2] + p[2]) / (p[1] + pa);	/* -2c/(b+Delta) */
		pa = -(p[1] + pa) / (p[0] + p[0]);	/* -(b+Delta)/(2a) */
	    }
	    else {
		pb = (-p[1] + pa) / (p[0] + p[0]);	/* (-b+Delta)/(2a) */
		pa = (p[2] + p[2]) / (-p[1] + pa);	/* 2c/(-b+Delta) */
	    }
	    if (pa <= pb) {
		p[0] = pa;
		p[1] = pb;
	    }
	    else {
		p[0] = pb;
		p[1] = pa;
	    }
	    iwk[0] = iwk[1] = 1;
	    ia = 0;
	    if (p[0] >= a && p[0] <= b)
		ia++;
	    else
		p[0] = p[1];
	    if (p[ia] >= a && p[ia] <= b)
		ia++;
	    return (ia);
	}
    }
    else if (FALSE && n == 3) {		/* cubic equation */
	/* Cardano's formula:
	 *	http://mathworld.wolfram.com/CubicFormula.html
	 *
	 * FIXME: robustness w.r.t. roundoff errors, see test cases in:
	 *		test/sturm_test.c
	 */
	FLOAT Q, R, D, S, T;
	static FLOAT d3=_F(1.)/_F(3.), d9=_F(1.)/_F(9.), d27=_F(1.)/_F(27.);
	static FLOAT twopi = _F(6.28318530717958647692528676655900577);
	pa = 1.0 / p[0];
	p[1] *= pa;
	p[2] *= pa;
	p[3] *= pa;
	Q = (3.0 * p[2] - p[1] * p[1]) * d9;
	R = (p[1] * (9.0 * p[2] - 2.0 * p[1] * p[1]) * d27 - p[3]) * 0.5;
	D = Q * Q * Q + R * R;
	p[0] = - p[1] * d3;	/* pre-compute -a2/3 and store in p[0] */
	if (D > EPS * EPS) {		/*------------------------ D > 0 */
	    /* one real root */
	    D = Sqrt(D);
	    S = R + D;
	    S = S >= 0.0 ? Pow(S, d3) : -Pow(-S, d3);
	    T = R - D;
	    T = T >= 0.0 ? Pow(T, d3) : -Pow(-T, d3);
	    p[0] += S + T;
	    iwk[0] = 1;
	    return check_range(1, a, b, p, iwk);
	}
	else if (D > - EPS * EPS) {	/*------------------------ D == 0 */
	    /* D == 0: 3 real roots, at least two of them are equal */
	    S = T = R >= 0.0 ? Pow(R, d3) : -Pow(-R, d3);
	    if (Fabs(S) > EPS) {
		/* one single root, one double root */
		if (S < 0.0) {
		    p[1] = p[0] - S;
		    iwk[1] = 2;
		    p[0] += S + T;
		    iwk[0] = 1;
		}
		else {
		    p[1] = p[0]  + S + T;
		    iwk[1] = 1;
		    p[0] -= S;
		    iwk[0] = 2;
		}
		return check_range(2, a, b, p, iwk);
	    }
	    else {
		/* one root with multiplicity 3 */
		p[0] += S;
		iwk[0] = 3;
		return check_range(1, a, b, p, iwk);
	    }
	}
	else {				/*------------------------ D < 0 */
	    /* three distinct real roots */
	    FLOAT theta;
	    iwk[0] = iwk[1] = iwk[2] = 1;
	    D = Sqrt(-Q);
	    theta = Acos(R / (-Q * D));
	    D += D;
	    p[2] = p[0] + D * Cos(theta * d3);
	    p[1] = p[0] + D * Cos((theta + twopi) * d3);
	    p[0] = p[0] + D * Cos((theta + twopi + twopi) * d3);
	    qsort(p, 3, sizeof(p[0]), phgCompFLOAT);
	    return check_range(3, a, b, p, iwk);
	}
    }

    if (sturm_init(p, n, tol, rwk, iwk))
	return (0);

    /* sturm_np must be greater than or equal to 2 */
    if (sturm_np == 2) {
	/* p'(x) is factor of p(x) ==> p(x) = p[0]*(x-r)^n
	 *	==> -n*r*p[0]==p[1] ==> r=-p[1]/(n*p[0]) */
	p[0] = pa = -p[1] / (n * p[0]);
	iwk[0] = n;
	return (pa >= a && pa <= b ? 1 : 0);
    }
    else {
	/* find ALL roots of the last polynomial of the Sturm sequence */
	int n0, n1;
	ib = 0;
	for (ia = 0; ia < sturm_np - 1; ia++)
	    ib += iwk[ia] + 1;
	n0 = iwk[ia];
	if (n0) {
	    int j, k;
	    FLOAT d;
	    for (ia = 0; ia <= n0; ia++) {
		rwk[ia] = rwk[ib + ia];
	    }

	    /* find lower and upper bounds for rwk(x) */
	    pa = 0.0;
	    pb = 0.0;
	    for (k = 1; k <= n0; k++) {
		d = Fabs(rwk[k]);
		if (pa < d)
		    pa = d;
		pb += d;
	    }
	    d = Fabs(rwk[0]);
	    pa = 1.0 + pa / d;
	    pb /= d;
	    if (pb < 1.0)
		pb = 1.0;
	    pa = 0.0001 + (pa < pb ? pa : pb);
	    recursion++;
	    n1 = sturm1(rwk, n0, -pa, pa, EPS, rwk + n0 + 1, iwk);
	    recursion--;
	    assert(n1 >= 0);

	    /* store roots found with multiplicity increased by 1,
	       and eliminate corresponding factors from p(x) */
	    for (ib = 0; ib < n1; ib++) {
		pa = rwk[ib];
		j = iwk[ib];
		for (ia = 0; ia <= j; ia++) {
		    sturm_factor(p, n, pa, rwk + n1);
		    for (k = 0; k < n; k++)
			p[k] = rwk[n1 + k];
		    n--;
		}
		if (pa >= a - EPS && pa <= b + EPS) {
		    if (FLAG)
			return 1;
		    for (ia = n + 1; ia > 0; ia--)
			p[ia] = p[ia - 1];
		    *(p++) = pa;
		    iwk[nr++] = j + 1;
		}
	    }			/* end of "for (ib=0;..." */
	    iwk += nr;

	    /* reinit the arrays (this will also restore sturm_rwk, etc.) */
	    if (sturm_init(p, n, tol, rwk, iwk))
		return (nr);
	}			/* end of "if (n0) ..." */
    }

#if 0
#warning
    printf("recursion = %d\n", recursion);
    poly_print(stdout, "reduced p=", n, p);
#endif

    /* at this point, p(x) can only have simple roots */
    ia = sturm_ncis(a, &pa);
    ib = sturm_ncis(b, &pb);
    if (FLAG)
	return ia > ib ? 1 : 0;
    sturm_nr = 0;
    sturm0(ia, a, pa, ib, b, pb);
    for (ia = 0; ia < sturm_nr; ia++)
	iwk[ia] = 1;
    return (nr + sturm_nr);
}

static void 
sturm_factor(FLOAT *p, int n, FLOAT a, FLOAT *rwk)
/* divide the polynomial p[x] by (x-a), return result in rwk,
   rwk[n]=remainder. */
{
    int i;
    rwk[0] = p[0];
    for (i = 1; i <= n; i++)
	rwk[i] = p[i] + a * rwk[i - 1];
}

static int 
sturm_init(FLOAT *p, int n, FLOAT tol, FLOAT *rwk, int *iwk)
/* returns !0 if p(x) is a constant polynomial */
{
    int i;
    int p0, p1, p2, n0, n1, n2;
    FLOAT d, eps;

    p0 = 0;
    /* p0(x)=p(x) */
    if (FALSE && p[0] != 1.0) {
	d = 1.0 / p[0];		/* divide by leading coefficient */
	for (i = 0; i <= n; i++)
	    rwk[i] = (p[i] *= d);
    }
    else {
	for (i = 0; i <= n; i++)
	    rwk[i] = p[i];
    }

    iwk[0] = n0 = n;
    p1 = n0 + 1;

    /* p1(x)=p'(x) */
    d = (n > 0 ? 1.0 / (n * Fabs(rwk[0])) : 0.0);
    for (i = 0; i < n; i++)
	rwk[p1 + i] = rwk[p0 + i] * (n - i) * d;
    iwk[1] = n1 = n - 1;
    p2 = p1 + n1 + 1;
    sturm_rwk = rwk;
    sturm_iwk = iwk;
    sturm_r = p;
    sturm_tol = tol;
    sturm_np = 2;

#if 0
#warning
    poly_print(stdout, "p0=", n0, rwk + p0);
    poly_print(stdout, "p1=", n1, rwk + p1);
#endif

    /* compute p_{k+2}(x): p_k(x)=q_k(x)*p_{k+1}(x)-p_{k+2}(x), k>=0 */
    eps = EPS * 100.0;
    do {
	n2 = sturm_remainder(n0, rwk + p0, n1, rwk + p1, rwk + p2);
	/* remove leading 0's in the remainder */
	while (n2 >= 0 && Fabs(rwk[p2]) <= eps) {
	    p2++;
	    n2--;
	}
#if 0
#warning
	poly_print(stdout, "p2=", n2, rwk + p2);
#endif
	if (n2 < 0)
	    break;
	d = 1.0 / Fabs(rwk[p2]);
	/*eps *= d;*/
#if 1
	/* normalize the leading coefficient of p2 to 1
	 * FIXME: (x-0.1)^3*(x-499999/5000000), or
	 * 	4 0. 1. 1 -0.3999998 0.05999994 -0.003999994 .0000999998
	 * fails when using __float128 */
	for (i = 0; i <= n2; i++)
	    rwk[p2+i] *= d;
# if 0
# warning
	poly_print(stdout, "normalized p2=", n2, rwk + p2);
# endif
#endif
	p0 = p1;
	p1 = p2;
	p2 += n2 + 1;
	n0 = n1;
	n1 = n2;
	iwk[sturm_np++] = n2;
    } while (1);

    return (0);
}

static int 
sturm_remainder(int n0, FLOAT *p0, int n1, FLOAT *p1, FLOAT *p2)
/* compute p2(x) such that p0(x)=q(x)*p1(x)-p2(x). returns the degree of
   p2(x) (returns -1 if p2(x)==0). The highest-order coefficients of p0
   and p1 are supposed nonzero */
{
    int i, j, n2;
    FLOAT d;
    if (!n0)
	return (-1);
    for (i = 0; i <= n0; i++)
	p2[i] = p0[i];
    n2 = n0;
    j = 0;
    while (n2 >= n1) {
	d = p2[j] / p1[0];
	for (i = 1; i <= n1; i++)
	    p2[j + i] -= p1[i] * d;
	j++;
	n2--;
	while (n2 >= 0 && Fabs(p2[j]) <= EPS) {
	    n2--;
	    j++;
	}
    }
    if (n2 >= 0)
	for (i = 0; i <= n2; i++)
	    p2[i] = -p2[j + i];

    return (n2);
}

static FLOAT
poly_eval(int n, FLOAT *p, FLOAT x)
/* evaluates the value of p(x) */
{
    int i;
    FLOAT d;
    d = p[0];
    for (i = 1; i <= n; i++)
	d = d * x + p[i];
    return (d);
}

static int 
sturm_ncis(FLOAT x, FLOAT *px)
/* returns the number of changes in sign. *px contains p0(x) */
{
    int i, s;
    FLOAT d0, d1, *p;
    *px = d0 = poly_eval(sturm_iwk[0], p = sturm_rwk, x);
    p += sturm_iwk[0] + 1;
    i = 1;
    while (d0 == 0.0 && i < sturm_np) {
	d0 = poly_eval(sturm_iwk[i], p, x);
	p += sturm_iwk[i] + 1;
	i++;
    }
    s = 0;
    while (i < sturm_np) {
	d1 = poly_eval(sturm_iwk[i], p, x);
	p += sturm_iwk[i] + 1;
	i++;
	if (d1 == 0.0)
	    continue;
	if ((d0 > 0.0) != (d1 > 0.0))
	    s++;
	d0 = d1;
    }
    return (s);
}


#ifndef BRENT	/* ========================================================= */

#if USE_NEWTON
static int
sturm_newton(int n, FLOAT *p, FLOAT *pprime, FLOAT *px,
					FLOAT a, FLOAT b, int maxit)
/* iterates using Newton's method,
 * if successful returns 0 with *px updated, otherwise *px is unchanged. */
{
    int i;
    FLOAT x0, x = *px, pp, p0, p0_save;

    p0 = p0_save = poly_eval(n, p, x);
    for (i = 0; i < maxit; i++) {
#if 0
#warning
	fprintf(stderr, "it=%d, x=%0.16lg, p(x)=%le\n",
				i, (double)x, (double)p0);
#endif
	pp = poly_eval(n - 1, pprime, x0 = x);
	if (Fabs(pp) /*<= EPS*/== 0.) {
#if 0
#warning
	    phgWarning("%s: derivative too small!\n", __func__);
	    poly_print_symbolic(stdout, "p(x) = ", n, p);
	    poly_print_symbolic(stdout, "p'(x) = ", n - 1, pprime);
	    printf("p'(%0.16g) = %g\n", (double)x, (double)pp);
#endif
	    return -1;
	}
	x = x0 - p0 / pp;
	if (Fabs(x - x0) <= 0.5 * sturm_tol * (Fabs(x) + Fabs(x0)))
	    break;
	p0 = poly_eval(n, p, x);
    }
    if (Fabs(p0) <= Fabs(p0_save) && x >= a - EPS && x <= b + EPS) {
	*px = x;
	return 0;	/* success */
    }
    return -1;
}
#endif /* USE_NEWTON */

static void 
sturm_bisec(FLOAT a, FLOAT pa, FLOAT b, FLOAT pb)
/* compute the root in the interval (a,b) to desired precision using
   a combination of bisection and Newton's method. */
{
    FLOAT c;
    int n;

#undef USE_NEWTON
#define USE_NEWTON 0	/* set to !0 to replace bisection with Newton iters */

#if USE_NEWTON
    FLOAT *pprime;		/* first derivative of the polynomial */
    FLOAT x;
    int i;

#endif /* USE_NEWTON */
    n = sturm_iwk[0];

#if USE_NEWTON
    pprime = malloc(n * sizeof(FLOAT));
    if (pprime != NULL) {
	for (i = 0; i < n; i++)
	    pprime[i] = (i + 1) * sturm_rwk[i + 1];
    }

#endif /* USE_NEWTON */

    do {
	c = (a + b) * 0.5;
	if (c == a || c == b || (b - a) <= sturm_tol * (Fabs(a) + Fabs(b))) {
	  end:
#if 0
	    for (n = 0; n < sturm_nr; n++)
		if (Fabs(c - sturm_r[n]) <= EPS ? ? ?)
		    return;

#endif /* 0 */
	    sturm_r[sturm_nr++] = c;

#if USE_NEWTON
	    if (pprime != NULL)
		free(pprime);
#endif /* USE_NEWTON */
	    return;
	}
	pb = poly_eval(n, sturm_rwk, c);
	if (pb == 0.0)
	    goto end;

#if USE_NEWTON
	/* try Newton's method */
	if (sturm_newton(n, sturm_rwk, pprime, &c, a, b, 10) == 0)
	    goto end;
#endif /* USE_NEWTON */

	if ((pb > 0.0) == (pa > 0.0)) {
	    a = c;
	    pa = pb;
	}
	else
	    b = c;
    } while (1);
}

#else /* undefined(BRENT) */

static FLOAT
sturm_func(FLOAT x)
{
    return poly_eval(sturm_iwk[0], sturm_rwk, x);
}

static void 
sturm_bisec(FLOAT a, FLOAT pa, FLOAT b, FLOAT pb)
{
    sturm_r[sturm_nr++] = zeroin(a, b, pa, pb, sturm_func, sturm_tol);
}
#endif	/* undefined(BRENT) */

static void 
sturm0(int ia, FLOAT a, FLOAT pa, int ib, FLOAT b, FLOAT pb)
/* locate recursively real roots of p0(x) */
{
    FLOAT c, pc;
    int ic;
    static short level = 0;
#if USE_OMP
# pragma omp threadprivate(level)
#endif  /* USE_OMP */

    if (level > 130) {
	fprintf(stderr, "%s:%d: recursion limit exceeded "
			"(a=%lg, ia=%d, b=%lg, ib=%d)\n",
			__FILE__, __LINE__, (double)a, ia, (double)b, ib);
	exit(1);
    }

    if (ia <= ib) {
	return;
    }
    if (ia == ib + 1 && (pa > 0.0) != (pb > 0.0)) {
	sturm_bisec(a, pa, b, pb);
    }
    else {
	c = (a + b) * 0.5;
	if (c == a || c == b || (b - a) <= sturm_tol * (pc = Fabs(a) + Fabs(b))
		|| pc < sturm_tol) {

#if 0
	    for (ic = 0; ic < sturm_nr; i++)
		if (Fabs(c - sturm_r[ic]) <= 1e-8 ? ? ?)
		    return;

#endif /* 0 */
	    sturm_r[sturm_nr] = c;
	    return;
	}
	ic = sturm_ncis(c, &pc);
	level++;
	sturm0(ia, a, pa, ic, c, pc);
	sturm0(ic, c, pc, ib, b, pb);
	level--;
    }
}
