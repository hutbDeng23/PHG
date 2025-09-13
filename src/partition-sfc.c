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

/* $Id: partition-sfc.c,v 1.88 2014/04/04 04:30:02 zlb Exp $ */

/* SFC Space filling curve partitioner */
/* Hilbert and Morton methods */
#include "phg.h"
#include "phg/partition-utils.h"

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* **************************************************************
 * Inverse Hilbert Space-Filling curve (3d)
 * converte coordinates in [0,1]x[0,1]x[0,1] to [0, 1]
 * n <= the length of array x 
 * 
 * Here we give two versions, the first one is general method
 * that based on table, which is based on Zoltan's.
 * We extented original implement, and we fixed some bugs.
 * The second one is based on an improved algorithm, which we
 * modified the first one. The time complexity is O(r), 
 * r = log2(max(x, y, z)) + 1. In general, r < n;
 * Both are fast. User can choose any one.
 * If grid has big coordinates much more than small coordinates, the 
 *   first maybe better.
 * If grid has small coordinats much more than bigs ones, the second
 *   maybe better.
 *
 * assume sizeof(INT) >= 4
 * **************************************************************/

/* data: idata2d, idata3d, istat2d, istat3d
 * are borrowed from zoltan.
 * They are used for generating inverse Hilbert SFC */

/* 2 dimension to nkey conversion */
static unsigned const int idata2d[] =  
 {0, 3, 1, 2,
  0, 1, 3, 2,
  2, 3, 1, 0,
  2, 1, 3, 0};

/* 2 dimension to nkey state transitions */
static unsigned const int istate2d[] = 
 {1, 2, 0, 0,
  0, 1, 3, 1,
  2, 0, 2, 3,
  3, 3, 1, 2};

/* 3 dimension to nkey conversion */
static unsigned const int idata3d [] = {
 0,  7,  3,  4,  1,  6,  2,  5,
 0,  1,  3,  2,  7,  6,  4,  5,
 0,  3,  7,  4,  1,  2,  6,  5,
 2,  3,  5,  4,  1,  0,  6,  7,
 4,  5,  3,  2,  7,  6,  0,  1,
 4,  7,  3,  0,  5,  6,  2,  1,
 6,  7,  5,  4,  1,  0,  2,  3,
 0,  1,  7,  6,  3,  2,  4,  5,
 2,  1,  5,  6,  3,  0,  4,  7,
 6,  1,  5,  2,  7,  0,  4,  3,
 0,  7,  1,  6,  3,  4,  2,  5,
 2,  1,  3,  0,  5,  6,  4,  7,
 4,  7,  5,  6,  3,  0,  2,  1,
 4,  5,  7,  6,  3,  2,  0,  1,
 6,  1,  7,  0,  5,  2,  4,  3,
 0,  3,  1,  2,  7,  4,  6,  5,
 2,  3,  1,  0,  5,  4,  6,  7,
 6,  7,  1,  0,  5,  4,  2,  3,
 2,  5,  1,  6,  3,  4,  0,  7,
 4,  3,  7,  0,  5,  2,  6,  1,
 4,  3,  5,  2,  7,  0,  6,  1,
 6,  5,  1,  2,  7,  4,  0,  3,
 2,  5,  3,  4,  1,  6,  0,  7,
 6,  5,  7,  4,  1,  2,  0,  3};


/* 3 dimension to nkey state transitions */
static unsigned const int istate3d [] ={
 1,  6,  3,  4,  2,  5,  0,  0,
 0,  7,  8,  1,  9,  4,  5,  1,
15, 22, 23, 20,  0,  2, 19,  2,
 3, 23,  3, 15,  6, 20, 16, 22,
11,  4, 12,  4, 20,  1, 22, 13,
22, 12, 20, 11,  5,  0,  5, 19,
17,  0,  6, 21,  3,  9,  6,  2,
10,  1, 14, 13, 11,  7, 12,  7,
 8,  9,  8, 18, 14, 12, 10, 11,
21,  8,  9,  9,  1,  6, 17,  7,
 7, 17, 15, 12, 16, 13, 10, 10,
11, 14,  9,  5, 11, 22,  0,  8,
18,  5, 12, 10, 19,  8, 12, 20,
 8, 13, 19,  7,  5, 13, 18,  4,
23, 11,  7, 17, 14, 14,  6,  1,
 2, 18, 10, 15, 21, 19, 20, 15,
16, 21, 17, 19, 16,  2,  3, 18,
 6, 10, 16, 14, 17, 23, 17, 15,
18, 18, 21,  8, 17,  7, 13, 16,
 3,  4, 13, 16, 19, 19,  2,  5,
16, 13, 20, 20,  4,  3, 15, 12,
 9, 21, 18, 21, 15, 14, 23, 10,
22, 22,  6,  1, 23, 11,  4,  3,
14, 23,  2,  9, 22, 23, 21,  0};

/*  maximum level */
/*  in this file, we force that MAXLEVEL >= 20 && MAXLEVEL <= 30 */
#define  MAXLEVEL  25

/* first algorithm */
#if 1

BOOLEAN
phgSFCInvHilbert3D(SFC *x, INT n)
{
    static unsigned const int *d[] =
	{idata3d,	idata3d + 8,   idata3d + 16,  idata3d + 24,
	 idata3d + 32,	idata3d + 40,  idata3d + 48,  idata3d + 56,
	 idata3d + 64,	idata3d + 72,  idata3d + 80,  idata3d + 88,
	 idata3d + 96,	idata3d + 104, idata3d + 112, idata3d + 120,
	 idata3d + 128, idata3d + 136, idata3d + 144, idata3d + 152,
	 idata3d + 160, idata3d + 168, idata3d + 176, idata3d + 184
	};

    static unsigned const int *s[] =
	{istate3d,	 istate3d + 8,	 istate3d + 16,  istate3d + 24,
	 istate3d + 32,  istate3d + 40,  istate3d + 48,  istate3d + 56,
	 istate3d + 64,  istate3d + 72,  istate3d + 80,  istate3d + 88,
	 istate3d + 96,  istate3d + 104, istate3d + 112, istate3d + 120,
	 istate3d + 128, istate3d + 136, istate3d + 144, istate3d + 152,
	 istate3d + 160, istate3d + 168, istate3d + 176, istate3d + 184
	};

    int level, EL;
    unsigned int key[3], c[3], temp, state;
    INT i;

    static unsigned  IMAX;
    static unsigned EfBit;
    static BOOLEAN initialized = FALSE;
    static int k0 = 0, k1 = 0, k2 = 0;

    if (!initialized) {
	initialized = TRUE;

	assert(sizeof(int) >= 4);
	assert(sizeof(unsigned int) >= 4);

	/* Init IMAX such that all effective bits are 1 */
	/* IMAX = 2^32 - 1 */
	IMAX = 4294967295U;

	/* 30 effective bits, all 1 */
	/* EfBit = 2^30 - 1 */
	EfBit = IMAX >> 2;


	k0 = 60 - MAXLEVEL * 3;
	k1 = 30 - MAXLEVEL * 3;
	k2 = - MAXLEVEL * 3;
    }

    for (i = 0; i < n; i++) {
	/* verification coordinats */
#if 0
	assert(x[i].co[0] >= 0. && x[i].co[0] <= 1.);
	assert(x[i].co[1] >= 0. && x[i].co[1] <= 1.);
	assert(x[i].co[2] >= 0. && x[i].co[2] <= 1.);
#endif
	/* convert x,y,z coordinates to integers in range [0,IMAX] */
	c[0] = (unsigned int)(x[i].co[0] * (double)IMAX);	 /* x */
	c[1] = (unsigned int)(x[i].co[1] * (double)IMAX);	 /* y */
	c[2] = (unsigned int)(x[i].co[2] * (double)IMAX);	 /* z */
	c[1] >>= 1;
	c[2] >>= 2;

	/* use state tables to convert nested quadrant's 
	 * coordinates level by level */
	key[0] = key[1] = key[2] = 0;
	state = 0;
	EL = 30;
	for (level = 0; level < MAXLEVEL; level++) {
	    /* extract 3 bits at current level */
	    EL--;
	    temp = ((c[0] >> EL) & 4)
		  |((c[1] >> EL) & 2)
		  |((c[2] >> EL) & 1);

	    /* treat key[] as long shift register */
	    /* shift in converted coordinate */
	    /* every key[] has thirty effective bit */
	    /* the last 30 are effective */
	    key[0] = (key[0] << 3) | ((key[1] >> 27) & 7);
	    key[1] = (key[1] << 3) | ((key[2] >> 27) & 7);
	    key[2] = (key[2] << 3) | *(d[state] + temp);

	    state = *(s[state] + temp);
	}
	
	key[0] = key[0] & EfBit;
	key[1] = key[1] & EfBit;
	key[2] = key[2] & EfBit;

	/* convert 3 part Hilbert key to double */
	x[i].sfc  = ldexp((double)key[2], k2);
	x[i].sfc += ldexp((double)key[1], k1);
	x[i].sfc += ldexp((double)key[0], k0);
    }

    return (i == n);
}
#else
/* log2(n)	       */
/* n >= 0	       */
/* return 0, if n <= 1 */
static unsigned int
LOG2(unsigned int n)
{
    unsigned int lo, hi, md;

    if (n <= 1) {
	return 0;
    }

    lo = 1;
    hi = MAXLEVEL;
    md = (lo + hi) >> 1;

    while(1) {
	if ((n >> md) > 0) {
	    lo = md;
	}
	else {
	    hi = md;
	}

	md = (lo + hi) >> 1;
	if (md == lo) {
	    break;
	}
    }

    return lo;
}

BOOLEAN
phgSFCInvHilbert3D(SFC * x, INT n)
{
    /* ordering data */
    static unsigned const int *d[] =
	{idata3d,	idata3d + 8,   idata3d + 16,  idata3d + 24,
	 idata3d + 32,	idata3d + 40,  idata3d + 48,  idata3d + 56,
	 idata3d + 64,	idata3d + 72,  idata3d + 80,  idata3d + 88,
	 idata3d + 96,	idata3d + 104, idata3d + 112, idata3d + 120,
	 idata3d + 128, idata3d + 136, idata3d + 144, idata3d + 152,
	 idata3d + 160, idata3d + 168, idata3d + 176, idata3d + 184
	};
    /* orientation(status) data */
    static unsigned const int *s[] =
	{istate3d,	 istate3d + 8,	 istate3d + 16,  istate3d + 24,
	 istate3d + 32,  istate3d + 40,  istate3d + 48,  istate3d + 56,
	 istate3d + 64,  istate3d + 72,  istate3d + 80,  istate3d + 88,
	 istate3d + 96,  istate3d + 104, istate3d + 112, istate3d + 120,
	 istate3d + 128, istate3d + 136, istate3d + 144, istate3d + 152,
	 istate3d + 160, istate3d + 168, istate3d + 176, istate3d + 184
	};
    /* auxiliary parameters, constants */
    static unsigned int  IMAX;
    static unsigned int  EfBit;
    static BOOLEAN initialized = FALSE;
    static int k0, k1, k2;

    int level, EL;
    unsigned int key[3], c[3], temp, state;
    unsigned int r;
    INT i;



    /* Initialize Inverse SFC parameters */
    if (!initialized) {
	if (phgVerbosity > 0) {
	    phgPrintf("Initialize Inverse SFC parameters\n");
	}

	initialized = TRUE;

	assert(sizeof(int) >= 4);
	assert(sizeof(unsigned int) >= 4);


	/* Init IMAX such that all effective bits are 1 */
	/* IMAX = 2^MAXLEVEL - 1 */
	/* Will be reduced later */
	IMAX = 1;
	IMAX <<= MAXLEVEL;
	IMAX -= 1;

	/* 30 effective bits, all 1 */
	/* EfBit = 2^30 - 1 */
	EfBit = 1;
	EfBit <<= 30;
	EfBit -= 1;


	/* index used by ldexp */
	/* effective number */
	k0 = 60 - MAXLEVEL * 3;
	k1 = 30 - MAXLEVEL * 3;
	k2 = - MAXLEVEL * 3;

    } /* end: initialize */

    for (i = 0; i < n; i++) {
	/* verificate coordinats */
#if 0
	assert(x[i].co[0] >= 0. && x[i].co[0] <= 1.);
	assert(x[i].co[1] >= 0. && x[i].co[1] <= 1.);
	assert(x[i].co[2] >= 0. && x[i].co[2] <= 1.);
#endif
	/* convert x,y,z coordinates to integers in range [0,IMAX] */
	c[0] = (unsigned int)(x[i].co[0] * (double)IMAX);	 /* x */
	c[1] = (unsigned int)(x[i].co[1] * (double)IMAX);	 /* y */
	c[2] = (unsigned int)(x[i].co[2] * (double)IMAX);	 /* z */
	
	/* r = log2(max(x, y, z)) */
	r = c[0];
	if (r < c[1])
	    r = c[1];
	if (r < c[2])
	    r = c[2];

	c[0] <<= 2;
	c[1] <<= 1;
	r = LOG2(r);

	/* use state tables to convert nested quadrant's 
	 * coordinates level by level */
	EL = r + 1;
	r = MAXLEVEL - EL;
	key[0] = key[1] = key[2] = 0;
	state = r % 2;

	for (level = r; level < MAXLEVEL; level++) {
	    /* extract 3 bits at current level */
	    EL--;
	    temp = ((c[0] >> (EL)) & 4)
		  |((c[1] >> (EL)) & 2)
		  |((c[2] >> (EL)) & 1);

	    /* treat key[] as long shift register */
	    /* shift in converted coordinate */
	    /* every key[] has thirty effective bit */
	    /* the last 30 are effective */
	    /* the last tree bits are effective */
	    key[0] = (key[0] << 3) | ((key[1] >> 27) & 7);
	    key[1] = (key[1] << 3) | ((key[2] >> 27) & 7);
	    key[2] = (key[2] << 3) | *(d[state] + temp);

	    state = *(s[state] + temp);
	}
	
	key[0] = key[0] & EfBit;
	key[1] = key[1] & EfBit;
	key[2] = key[2] & EfBit;

	/* convert 3 part Hilbert key to double */
	x[i].sfc  = ldexp((double)key[2], k2);
	x[i].sfc += ldexp((double)key[1], k1);
	x[i].sfc += ldexp((double)key[0], k0);
    }

    return (i == n);
}
#endif

/* inverse Hilbert SFC (2d) */
BOOLEAN
phgSFCInvHilbert2D(SFC2 *x, INT n)
{
    static unsigned const int *d[] =
	   {idata2d,  idata2d  +4, idata2d  +8, idata2d  +12};
    static unsigned const int *s[] =
	   {istate2d, istate2d +4, istate2d +8, istate2d +12};
    int level;
    unsigned int key[2], c[2], temp, state;
    INT i;

    static unsigned  IMAX;
    static BOOLEAN initialized = FALSE;
    static int k0 = 0, k1 = 0;

    if (!initialized) {
	initialized = TRUE;

	assert(sizeof(int) >= 4);
	assert(sizeof(unsigned int) >= 4);


	/* Init IMAX such that all effective bits are 1 */
	/* IMAX = 2^32 - 1 */
	IMAX = 4294967295U;

	k1 =  - MAXLEVEL * 2;
	k0 = 32 + k1;
    }

    for (i = 0; i < n; i++) {
	/* verification coordinats */
#if 0
	assert(x[i].co[0] >= 0. && x[i].co[0] <= 1.);
	assert(x[i].co[1] >= 0. && x[i].co[1] <= 1.);
#endif
	/* convert x,y coordinates to integers in range [0,IMAX] */
	c[0] = (unsigned int)(x[i].co[0] * (double)IMAX);	 /* x */
	c[1] = (unsigned int)(x[i].co[1] * (double)IMAX);	 /* y */

	/* use state tables to convert nested quadrant's 
	 * coordinates level by level */
	key[0] = key[1] = 0;
	state = 0;
	for (level = 0; level < MAXLEVEL; level++) {
	    /* extract 2 bits at current level */
	    temp = ((c[0] >> (30 - level)) & 2)
		  |((c[1] >> (31 - level)) & 1);

	    /* treat key[] as long shift register */
	    /* shift in converted coordinate */
	    key[0] = (key[0] << 2) | ((key[1] >> 30) & 3);
	    key[1] = (key[1] << 2) | *(d[state] + temp);
	    state = *(s[state] + temp);
	}
	
	key[0] = key[0] & IMAX;
	key[1] = key[1] & IMAX;

	/* convert 2 part Hilbert key to double */
	x[i].sfc =  ldexp((double)key[0], k0) 
		  + ldexp((double)key[1], k1);
    }

    return (i == n);
}

/* inverse Morton SFC (3d) */
BOOLEAN
phgSFCInvMorton3D(SFC *x, INT n)
{
    int level, EL;
    unsigned int key[3], c[3], temp;
    INT i;

    static unsigned int  IMAX;
    static unsigned int  EfBit;
    static BOOLEAN initialized = FALSE;
    static int k0 = 0, k1 = 0, k2 = 0;

    if (!initialized) {
	initialized = TRUE;

	assert(sizeof(int) >= 4);
	assert(sizeof(unsigned int) >= 4);

	/* Init IMAX such that all effective bits are 1 */
	/* IMAX = 2^32 - 1 */
	IMAX = 4294967295U;

	/* 30 effective bits, all 1 */
	/* EfBit = 2^30 - 1 */
	EfBit = IMAX >> 2;


	k0 = 60 - MAXLEVEL * 3;
	k1 = 30 - MAXLEVEL * 3;
	k2 = - MAXLEVEL * 3;
    }

    for (i = 0; i < n; i++) {
	/* verification coordinats */
#if 0
	assert(x[i].co[0] >= 0. && x[i].co[0] <= 1.);
	assert(x[i].co[1] >= 0. && x[i].co[1] <= 1.);
	assert(x[i].co[2] >= 0. && x[i].co[2] <= 1.);
#endif
	/* convert x,y,z coordinates to integers in range [0,IMAX] */
	c[0] = (unsigned int)(x[i].co[0] * (double)IMAX);	 /* x */
	c[1] = (unsigned int)(x[i].co[1] * (double)IMAX);	 /* y */
	c[2] = (unsigned int)(x[i].co[2] * (double)IMAX);	 /* z */
	c[0] >>= 2;
	c[1] >>= 1;

	/* use state tables to convert nested quadrant's 
	 * coordinates level by level */
	key[0] = key[1] = key[2] = 0;
	EL = 30;
	for (level = 0; level < MAXLEVEL; level++) {
	    /* extract 3 bits at current level */
	    EL--;
	    temp = ((c[0] >> (EL)) & 1)
		  |((c[1] >> (EL)) & 2)
		  |((c[2] >> (EL)) & 4);

	    /* treat key[] as long shift register */
	    /* shift in converted coordinate */
	    /* every key[] has thirty effective bit */
	    /* the last 30 are effective */
	    key[0] = (key[0] << 3) | ((key[1] >> 27) & 7);
	    key[1] = (key[1] << 3) | ((key[2] >> 27) & 7);
	    key[2] = (key[2] << 3) | (temp);

	}
	
	key[0] = key[0] & EfBit;
	key[1] = key[1] & EfBit;
	key[2] = key[2] & EfBit;

	/* convert 3 part Morton key to double */
	x[i].sfc  = ldexp((double)key[2], k2);
	x[i].sfc += ldexp((double)key[1], k1);
	x[i].sfc += ldexp((double)key[0], k0);
    }

    return (i == n);
}

/* inverse Morton SFC (2d) */
BOOLEAN
phgSFCInvMorton2D(SFC2 *x, INT n)
{
    int level;
    unsigned int key[2], c[2], temp;
    INT i;

    static unsigned int  IMAX;
    static BOOLEAN initialized = FALSE;
    static int k0 = 0, k1 = 0;

    if (!initialized) {
	initialized = TRUE;

	assert(sizeof(int) >= 4);
	assert(sizeof(unsigned int) >= 4);

	/* Init IMAX such that all effective bits are 1 */
	/* IMAX = 2^32 - 1 */
	IMAX = 4294967295U;

	k1 =  - MAXLEVEL * 2;
	k0 = 32 + k1;
    }

    for (i = 0; i < n; i++) {
	/* verification coordinats */
#if 0
	assert(x[i].co[0] >= 0. && x[i].co[0] <= 1.);
	assert(x[i].co[1] >= 0. && x[i].co[1] <= 1.);
#endif
	/* convert x,y coordinates to integers in range [0,IMAX] */
	c[0] = (unsigned int)(x[i].co[0] * (double)IMAX);	 /* x */
	c[1] = (unsigned int)(x[i].co[1] * (double)IMAX);	 /* y */

	/* use state tables to convert nested quadrant's 
	 * coordinates level by level */
	key[0] = key[1] = 0;
	for (level = 0; level < MAXLEVEL; level++) {
	    /* extract 2 bits at current level */
	    temp = ((c[0] >> (31 - level)) & 1)
		  |((c[1] >> (30 - level)) & 2);

	    /* treat key[] as long shift register */
	    /* shift in converted coordinate */
	    key[0] = (key[0] << 2) | ((key[1] >> 30) & 3);
	    key[1] = temp;
	}
	
	key[0] = key[0] & IMAX;
	key[1] = key[1] & IMAX;

	/* convert 2 part Hilbert key to double */
	x[i].sfc =  ldexp((double)key[0], k0) 
		  + ldexp((double)key[1], k1);
    }

    return (i == n);
}

#if USE_MPI


static const char *sfc_types[] = {"hilbert", "morton", NULL};
enum {HILBERT = 0, MORTON = 1};
static int sfc_type = 0;

static int sfc_remap = TRUE;

/* If ALLOW_EMPTY_PARTITIONS is set to 0 then phgPartitionSFC will try to avoid
 * generating empty partitions. */
#define ALLOW_EMPTY_PARTITIONS	1

/* Hilbert space-filling curve partitioner */
/* for 3d space, can be extended to 2d and 3d with slight change */
BOOLEAN
phgPartitionSFC(GRID *g, MPI_Comm newcomm, DOF *weights, FLOAT power)
{
    int nprocs, rank;
    DOTS *dots;
    ELEMENT *e;
    BOOLEAN ret;
    SFC *x = NULL;
    FLOAT out[6];	/* min and max coordinats */
    FLOAT ext[3];	/* scaling factor */
    FLOAT temp;
    FLOAT lif = 1.0;
    FLOAT *speeds = NULL;
    INT i, nleaf;
    FLOAT x0, y0, z0;
    BOOLEAN fail;
    static const FLOAT lam[] = {0.25, 0.25, 0.25, 0.25};
    static BOOLEAN initialized = FALSE;

    if (!initialized) {
	initialized = TRUE;
	phgOptionsRegisterTitle("\n[SFC options]", "\n", "partition");
	phgOptionsRegisterKeyword("-sfc_method", "SFC method",
			sfc_types, &sfc_type);
	phgOptionsRegisterNoArg("-sfc_remap",
			"Remap partitions to minimize data move",
			&sfc_remap);
	return FALSE;
    }

    MPI_Comm_size(newcomm, &nprocs);
    MPI_Comm_rank(newcomm, &rank);

    if (nprocs == 1 && g->nprocs == 1)
	return FALSE;

#if 0
    /* For testing only: assign different speeds to processes */
    temp = (SFC_FLOAT)/*(1 << rank)*/1.0 / (rank + 1.0);
    ret = FALSE;
#else
    temp = phgGetProcessorSpeed(NULL);
    ret = phgSameProcessorSpeed();
#endif
    if (rank == 0) {
	speeds = phgAlloc(nprocs * sizeof(*speeds));
	if (ret) {
	    for (i = 0; i < nprocs; i++)
		speeds[i] = 1.0;
	}
	else {
	    MPI_Gather(&temp, 1, PHG_MPI_FLOAT, speeds, 1, PHG_MPI_FLOAT,
				0, newcomm);
	}
	/* normalize speeds */
	temp = speeds[0];
	for (i = 1; i < nprocs; i++)
	    temp += speeds[i];
	for (i = 0; i < nprocs; i++)
	    speeds[i] /= temp;
    }
    else if (!ret) {
	MPI_Gather(&temp, 1, PHG_MPI_FLOAT, speeds, 1, PHG_MPI_FLOAT,
			0, newcomm);
    }

    ret = FALSE;
    if (g->nprocs <= rank)
	goto label;

    /* allocate space for SFC */
    nleaf = g->nleaf;
    x = (SFC *)phgAlloc(nleaf * sizeof(SFC));
    if (nleaf > 0 && x == NULL) {
	phgError(1, "memory allocation error (%s:%d)\n", __FILE__, __LINE__);
	return ret;
    }

    /* calculate coordinates */
    if (phgVerbosity > 0) {
	phgPrintf("SFC - init objects coordinates ... ");
    }
    i = 0;
    ForAllElements(g, e) {
	phgGeomLambda2XYZ(g, e, lam, &x0, &y0, &z0);
	x[i].co[0] = x0;	  /* x */
	x[i].co[1] = y0;	  /* y */
	x[i].co[2] = z0;	  /* z */

	i++;
    }  /* scan all leaves */

    /* bounding box */
    out[0] = g->bbox[0][0];
    out[1] = g->bbox[0][1];
    out[2] = g->bbox[0][2];
    out[3] = g->bbox[1][0];
    out[4] = g->bbox[1][1];
    out[5] = g->bbox[1][2];

    /* enlarge bounding box */
    temp = (out[3] - out[0]) * HSFC_EPSILON;
    out[3] += temp;
    out[0] -= temp;
    temp = (out[4] - out[1]) * HSFC_EPSILON;
    out[4] += temp;
    out[1] -= temp;
    temp = (out[5] - out[2]) * HSFC_EPSILON;
    out[5] += temp;
    out[2] -= temp;

    /* scaling coordinates, such that all coordinates belong to (0, 1)^3 */
    ext[0] = out[3] - out[0]; /* length of x axis of bounding box */
    ext[1] = out[4] - out[1]; /* length of y axis of bounding box */
    ext[2] = out[5] - out[2]; /* length of z axis of bounding box */

    /* preserve domain ratio */
#if 1
    if (ext[0] < ext[1])
	ext[0] = ext[1];
    if (ext[0] < ext[2])
	ext[0] = ext[2];

    if (ext[0] == 0.)
	ext[0] = 1.;
    ext[0] = 1. / ext[0];
    ext[1] =  ext[0];
    ext[2] =  ext[0];
#else
    /* NOT preserve domain ratio */
    ext[0] = 1. / ext[0];
    ext[1] = 1. / ext[1];
    ext[2] = 1. / ext[2];
#endif

    for (i = 0; i < nleaf; i++) {
	x[i].co[0] = (x[i].co[0] - out[0]) * ext[0];
	x[i].co[1] = (x[i].co[1] - out[1]) * ext[1];
	x[i].co[2] = (x[i].co[2] - out[2]) * ext[2];
    }  /* coordinates \in (0, 1)^3 */

    if (phgVerbosity > 0) {
	phgPrintf("ok\n");
    }

    /* generate inverse SFC */
    temp = phgGetTime(NULL);
    if (sfc_type == HILBERT) {
	if (phgVerbosity > 0) {
	    phgPrintf("SFC - Hilbert SFC\n");
	}
	ret = phgSFCInvHilbert3D(x, nleaf);
    }
    else if (sfc_type == MORTON) {
	if (phgVerbosity > 0) {
	    phgPrintf("SFC - Morton SFC\n");
	}
	ret = phgSFCInvMorton3D(x, nleaf);
    }
    temp = phgGetTime(NULL) - temp;
    if (phgVerbosity > 0) {
	phgPrintf("SFC - number of elements: %"dFMT" (in process: %d)\n",
			g->nleaf, rank);
	phgPrintf("SFC - time of generating SFC: %.5fs\n", temp);
    }


    if (!ret) {
	phgFree(x);
	phgError(1, "SFC - Error in generating inverse SFC\n");
	return ret;
    }
    if (phgVerbosity > 0) {
	phgPrintf("SFC - generate inverse SFC ... ok\n");
    }
    
    /* initialize dots */
    dots = (DOTS *)phgAlloc(nleaf * sizeof(DOTS));
    if (nleaf > 0 && dots == NULL) {
	phgFree(x);
	phgError(1, "memory allocation error (%s:%d)\n", __FILE__, __LINE__);
	return ret;
    }
    i = 0;
    ForAllElements(g, e) {
	dots[i].key = x[i].sfc;
	i++;
    }
    phgFree(x);

    /* weighted condition */
    if (weights != NULL && power != 0.) {
	SFC_FLOAT wtp = -1e10, wmax;

	i = 0;
	ForAllElements(g, e) {
	    dots[i].w = Pow(*DofElementData(weights, e->index), power);
	    if (dots[i].w > wtp)
		wtp = dots[i].w;
	    i++;
	}

	/* standardize weights */
	MPI_Allreduce(&wtp, &wmax, 1, MPI_SFC_FLOAT, MPI_MAX, g->comm);
	if (wmax > 0.) {
	    wmax = 1. / wmax * 100.;
	    for (i = 0; i < nleaf; i++) {
		dots[i].w *= wmax;
	    }
	}
	else {
	    for (i = 0; i < nleaf; i++) {
		dots[i].w = 1.;
	    }
	}
    }
    else { /* no weights DOF, all are 1. */
	for (i = 0; i < nleaf; i++) {
	    dots[i].w = 1.;
	}
    }
   
    if (phgVerbosity > 0) {
	phgPrintf("SFC - 1D partitioning:\n");
    }
    /* 1d  partition */
    if (g->nleaf_global <= nprocs) {
	ForAllElements(g, e) {
	    e->mark = GlobalElement(g, e->index) % nprocs;
	}
	/* if lif == -1 then it will be computed in distribute.c */
	lif = phgSameProcessorSpeed() ? nprocs * 1. / g->nleaf_global : -1.0;
	fail = FALSE;
    }
    else {
#if !ALLOW_EMPTY_PARTITIONS
	INT cl[nprocs];
	INT cg[nprocs];
	int kk;
	INT nmax;
#endif	/* !ALLOW_EMPTY_PARTITIONS */

	lif = phgPartition1DV1(dots, nleaf, nprocs, speeds, g->comm);

	/* assign */
	i = 0;
	ForAllElements(g, e) {
	    e->mark = dots[i].pn;
	    assert(e->mark >= 0 && e->mark < nprocs);
	    i++;
	}

#if ALLOW_EMPTY_PARTITIONS
	fail = FALSE;
#else	/* ALLOW_EMPTY_PARTITIONS */
	/* check if some partition has no element 
	 * if yes, then fails. or succeeds */
	for (kk = 0; kk < nprocs; kk++) {
	    cl[kk] = 0;
	}
	ForAllElements(g, e) {
	    cl[e->mark] += 1;
	}

	MPI_Allreduce(cl, cg, nprocs, PHG_MPI_INT, MPI_SUM, g->comm);
	nmax = cg[0];
	for (kk = 1; kk < nprocs; kk++) {
	    if (cg[kk] < cg[0])
		cg[0] = cg[kk];
	    if (cg[kk] > nmax)
		nmax = cg[kk];
	}

	if ((cg[0] == 0) && (nmax > 0)) {
	    fail = TRUE;
	}
	else {
	    fail = FALSE;
	}
#endif	/* ALLOW_EMPTY_PARTITIONS */

	/* if V1 fails, switch to V2 */
	if (fail) {
	    if (rank == 0)
		phgWarning("SFC - phgPartition1DV1 failed, retry with phgPartition1DV2.\n");
	    lif = phgPartition1DV2(dots, nleaf, nprocs, speeds, g->comm);

	    /* assign */
	    i = 0;
	    ForAllElements(g, e) {
		e->mark = dots[i].pn;
		assert(e->mark >= 0 && e->mark < nprocs);
		i++;
	    }

#if ALLOW_EMPTY_PARTITIONS
	    fail = FALSE;
#else	/* ALLOW_EMPTY_PARTITIONS */
	    /* check if some partition has no element 
	     * if yes, then fail. or succeed */
	    for (kk = 0; kk < nprocs; kk++) {
		cl[kk] = 0;
	    }
	    ForAllElements(g, e) {
		cl[e->mark] += 1;
	    }

	    MPI_Allreduce(cl, cg, nprocs, PHG_MPI_INT, MPI_SUM, g->comm);
	    nmax = cg[0];
	    for (kk = 1; kk < nprocs; kk++) {
		if (cg[kk] < cg[0])
		    cg[0] = cg[kk];
		if (cg[kk] > nmax)
		    nmax = cg[kk];
	    }

	    if ((cg[0] == 0) && (nmax > 0)) {
		fail = TRUE;
	    }
	    else {
		fail = FALSE;
	    }
#endif	/* ALLOW_EMPTY_PARTITIONS */
	}
    }
    phgFree(dots);

    /* if 1d partitioners fail,
     * partition the grid according to its index */
    if (fail) {
	INT kk = g->nleaf_global / nprocs;
	INT rr = g->nleaf_global % nprocs;
	INT j;

	if (rank == 0) {
	    phgWarning("SFC - 1D partitioner failed.\n");
	    phgWarning("SFC - grid partitioned according to element index.\n");
	}

	j = (kk + 1) * rr;
	ForAllElements(g, e) {
	    i = GlobalElement(g, e->index);
	    if (i <= j) {
		e->mark = i / (kk + 1);
	    }
	    else {
		e->mark = rr + (i - j) / kk;
	    }
	    assert(e->mark >= 0 && e->mark < nprocs);
	}

	if (!phgSameProcessorSpeed()) {
	    lif = -1.0;		/* will be computed in distribute.c */
	}
	else if (rr == 0) {
	    lif = 1.;
	}
	else {
	    lif = (kk + 1) * nprocs * 1. / g->nleaf_global;
	}
    }

    if (phgVerbosity > 2) {
	INT c[nprocs];
	INT ct[nprocs];
	int kk;

	for (kk = 0; kk < nprocs; kk++) {
	    c[kk] = 0;
	}
	ForAllElements(g, e) {
	    assert(e->mark >= 0);
	    assert(e->mark < nprocs);

	    c[e->mark] += 1;
	}

	MPI_Allreduce(c, ct, nprocs, PHG_MPI_INT, MPI_SUM, g->comm);
	for (kk = 0; kk < nprocs; kk++) {
	    phgPrintf("Rank: %6d     Elements:	 %"dFMT"\n", kk, ct[kk]);
	}
    }
    if (phgVerbosity > 0)
	phgPrintf("SFC - 1D partitioning ... ok\n");

    /* remap partitions */
label:
    if (sfc_remap) {
	if (phgVerbosity > 0)
	    phgPrintf("SFC - Remapping partitions\n");
	i = phgPartitionRemap(g, newcomm);
    }
    else {
	i = 0;
    }

    if (rank == 0) {
	phgFree(speeds);
	speeds = NULL;
    }

    if (g->nprocs <= rank)
	return 0;

    if (i > 0) {
	if (phgVerbosity > 2) {
	    INT c[nprocs];
	    INT ct[nprocs];
	    int kk;

	    for (kk = 0; kk < nprocs; kk++) {
		c[kk] = 0;
	    }
	    ForAllElements(g, e) {
		assert(e->mark >= 0);
		assert(e->mark < nprocs);

		c[e->mark] += 1;
	    }

	    MPI_Allreduce(c, ct, nprocs, PHG_MPI_INT, MPI_SUM, g->comm);
	    for (kk = 0; kk < nprocs; kk++) {
		phgPrintf("Rank: %6d	 Elements:   %"dFMT"\n", kk, ct[kk]);
	    }
	}
    }

    if (phgVerbosity > 0) {
	if (i > 0) {
	    phgPrintf("SFC - Remapping partitions ... ok\n");
	}
    }

    /* update information */
    i = 0;
    ForAllElements(g, e) {
	if (e->mark == rank)
	    i++;
    }
    if (i == g->nleaf) {
	nleaf = 0;
    }
    else {
	nleaf = 1;
    }

    MPI_Allreduce(&nleaf, &i, 1, PHG_MPI_INT, MPI_SUM, g->comm);

    if (i == 0) { /* grid unchanged */
	ret = FALSE;
	if (phgVerbosity > 0)
	    phgPrintf("SFC - %s: ** Grid Unchanged **\n", __func__);
    }
    else {	  /* grid changed */
	ret = TRUE;
	g->lif = lif;
    }


    return ret;
}

#endif
