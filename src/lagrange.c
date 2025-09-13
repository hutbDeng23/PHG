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

/* $Id: lagrange.c,v 1.43 2022/04/09 01:28:34 zlb Exp $ */

/* Lagrange finite elements: DOF_P0, DOF_P1, DOF_P2, DOF_P3, DOF_P4, ...
 *
 * Number of terms of Pn = (n+1)(n+2)(n+3) / 3!
 */

#include "phg.h"

#include <string.h>

/* DEBUG_LAGRANGE1!=0 ==> use code in lagrange1.c (for debugging) */
#define DEBUG_LAGRANGE1		0

/*--------------------------- P0 element ----------------------------*/

static const FLOAT *
P0_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of basis functions */
{
    static FLOAT value = 1.0;

    Unused(lambda);
    Unused(dof);
    Unused(e);
    assert(no0 == 0 && no1 <= 1);

    return &value;
}

static const FLOAT *
P0_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of gradients of basis functions */
{
    static FLOAT values[] = {0., 0., 0., 0.};
    
    Unused(dof);
    Unused(e);
    assert(no0 == 0 && no1 <= 1);

    return values;
}

static void
P0_interp(DOF *dof, ELEMENT *e, FLOAT **old_data, FLOAT **new_data)
/* interpolation */
{
    int i;
    FLOAT *vnew0 = new_data[10], *vnew1 = new_data[11], *vold = old_data[14];

    for (i = 0; i < dof->dim; i++)
	*(vnew0++) = *(vnew1++) = *(vold++);
}

#if !DEBUG_LAGRANGE1

/*--------------------------- P1 element ----------------------------*/

static const FLOAT *
P1_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of basis functions */
{
    static FLOAT *buffer = NULL;
    FLOAT *values;

    Unused(dof);
    Unused(e);

    if (buffer == NULL) {
#if USE_OMP
#pragma omp critical (P1_bas)
#endif	/* USE_OMP */
	if (buffer == NULL) {
	    buffer = phgAlloc(phgMaxThreads * NVert * sizeof(*buffer));
	    FreeAtExitNoCheck(buffer);
	}
    }

    values = buffer + phgThreadId * NVert;

    if (no1 <= 0)
	no1 = NVert;
    assert(no0 < no1);

    memcpy(values, lambda + no0, (no1 - no0) * sizeof(values[0]));

    return values;
}

static const FLOAT *
P1_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of gradients of basis functions */
{
    static FLOAT (*buffer)[Dim + 1] = NULL;
    FLOAT (*values)[Dim + 1];
    int i, j;
 
    Unused(dof);
    Unused(e);

    if (buffer == NULL) {
#if USE_OMP
#pragma omp critical (P1_bas)
#endif	/* USE_OMP */
	if (buffer == NULL) {
	    buffer = phgAlloc(phgMaxThreads * NVert * sizeof(*buffer));
	    FreeAtExitNoCheck(buffer);
	}
    }

    values = buffer + phgThreadId * NVert;

    if (no1 <= 0)
	no1 = NVert;
    assert(no0 < no1);

    for (j = no0; j < no1; j++)
	for (i = 0; i < Dim + 1; i++)
	    values[j][i] = (i == j) ? 1.0 : 0.0;

    return (FLOAT *)(values + no0);
}

static void
P1_interp(DOF *dof, ELEMENT *e, FLOAT **old_data, FLOAT **new_data)
/* interpolation */
{
    int i;
    FLOAT *vnew = new_data[0], *vold0 = old_data[0], *vold1 = old_data[1];

    Unused(e);

    for (i = 0; i < dof->dim; i++)
	*(vnew++) = 0.5 * (*(vold0++) + *(vold1++));
}

/*--------------------------- P2 element ----------------------------*/

static const FLOAT *
P2_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of basis functions */
{
    static FLOAT *buffer = NULL;
    FLOAT *p, *values;
    FLOAT L0 = lambda[0], L1 = lambda[1], L2 = lambda[2], L3 = lambda[3];

    Unused(dof);
    Unused(e);

    if (buffer == NULL) {
#if USE_OMP
#pragma omp critical (P2_bas)
#endif	/* USE_OMP */
	if (buffer == NULL) {
	    buffer = phgAlloc(phgMaxThreads * (NVert + NEdge)
					    * sizeof(*buffer));
	    FreeAtExitNoCheck(buffer);
	}
    }

    p = values = buffer + phgThreadId * (NVert + NEdge);

    if (no1 <= 0)
	no1 = 10;
    assert(no0 < no1);

    if (no0 == 0 && no1 == 10) {
	FLOAT L0 = lambda[0], L1 = lambda[1], L2 = lambda[2], L3 = lambda[3];
	*(p++) = L0 * (L0 + L0 - 1.);
	*(p++) = L1 * (L1 + L1 - 1.);
	*(p++) = L2 * (L2 + L2 - 1.);
	*(p++) = L3 * (L3 + L3 - 1.);
	*(p++) = 4. * L0 * L1;
	*(p++) = 4. * L0 * L2;
	*(p++) = 4. * L0 * L3;
	*(p++) = 4. * L1 * L2;
	*(p++) = 4. * L1 * L3;
	*p     = 4. * L2 * L3;
	return values;
    }

    switch (no0) {
	case 0:
	    *(p++) = L0 * (L0 + L0 - 1.);
	    if (no1 == 1)
		return values;
	case 1:
	    *(p++) = L1 * (L1 + L1 - 1.);
	    if (no1 == 2)
		return values;
	case 2:
	    *(p++) = L2 * (L2 + L2 - 1.);
	    if (no1 == 3)
		return values;
	case 3:
	    *(p++) = L3 * (L3 + L3 - 1.);
	    if (no1 == 4)
		return values;
	case 4:
	    *(p++) = 4. * L0 * L1;
	    if (no1 == 5)
		return values;
	case 5:
	    *(p++) = 4. * L0 * L2;
	    if (no1 == 6)
		return values;
	case 6:
	    *(p++) = 4. * L0 * L3;
	    if (no1 == 7)
		return values;
	case 7:
	    *(p++) = 4. * L1 * L2;
	    if (no1 == 8)
		return values;
	case 8:
	    *(p++) = 4. * L1 * L3;
	    if (no1 == 9)
		return values;
	case 9:
	    *(p++) = 4. * L2 * L3;
		return values;
    }

    return values;	/* to make gcc happy */
}

static const FLOAT *
P2_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of gradients of basis functions */
{
    static FLOAT (*buffer)[Dim + 1] = NULL;
    FLOAT (*values)[Dim + 1];
    FLOAT (*p)[Dim + 1];

    Unused(dof);
    Unused(e);

    if (buffer == NULL) {
#if USE_OMP
#pragma omp critical (P2_bas)
#endif	/* USE_OMP */
	if (buffer == NULL) {
	    buffer = phgAlloc(phgMaxThreads * (NVert + NEdge)
					    * sizeof(*buffer));
	    FreeAtExitNoCheck(buffer);
	}
    }

    p = values = buffer + phgThreadId * (NVert + NEdge);

    if (no1 <= 0)
	no1 = 10;
    assert(no0 < no1);

    switch (no0) {
	case 0:	/* vertex 0 */
	    (*p)[0] = 4. * lambda[0] - 1.;
	    (*p)[1] = 0.;
	    (*p)[2] = 0.;
	    (*p)[3] = 0.;
	    if (no1 == 1)
		return (FLOAT *)values;
	    p++;
	case 1:	/* vertex 1 */
	    (*p)[0] = 0.;
	    (*p)[1] = 4. * lambda[1] - 1.;
	    (*p)[2] = 0.;
	    (*p)[3] = 0.;
	    if (no1 == 2)
		return (FLOAT *)values;
	    p++;
	case 2:	/* vertex 2 */
	    (*p)[0] = 0.;
	    (*p)[1] = 0.;
	    (*p)[2] = 4. * lambda[2] - 1.;
	    (*p)[3] = 0.;
	    if (no1 == 3)
		return (FLOAT *)values;
	    p++;
	case 3:	/* vertex 3 */
	    (*p)[0] = 0.;
	    (*p)[1] = 0.;
	    (*p)[2] = 0.;
	    (*p)[3] = 4. * lambda[3] - 1.;
	    if (no1 == 4)
		return (FLOAT *)values;
	    p++;
	case 4:	/* edge 0-1 */
	    (*p)[0] = 4. * lambda[1];
	    (*p)[1] = 4. * lambda[0];
	    (*p)[2] = 0.;
	    (*p)[3] = 0.;
	    if (no1 == 5)
		return (FLOAT *)values;
	    p++;
	case 5:	/* edge 0-2 */
	    (*p)[0] = 4. * lambda[2];
	    (*p)[1] = 0.;
	    (*p)[2] = 4. * lambda[0];
	    (*p)[3] = 0.;
	    if (no1 == 6)
		return (FLOAT *)values;
	    p++;
	case 6:	/* edge 0-3 */
	    (*p)[0] = 4. * lambda[3];
	    (*p)[1] = 0.;
	    (*p)[2] = 0.;
	    (*p)[3] = 4. * lambda[0];
	    if (no1 == 7)
		return (FLOAT *)values;
	    p++;
	case 7:	/* edge 1-2 */
	    (*p)[0] = 0;
	    (*p)[1] = 4. * lambda[2];
	    (*p)[2] = 4. * lambda[1];
	    (*p)[3] = 0.;
	    if (no1 == 8)
		return (FLOAT *)values;
	    p++;
	case 8:	/* edge 1-3 */
	    (*p)[0] = 0;
	    (*p)[1] = 4. * lambda[3];
	    (*p)[2] = 0.;
	    (*p)[3] = 4. * lambda[1];
	    if (no1 == 9)
		return (FLOAT *)values;
	    p++;
	case 9:	/* edge 2-3 */
	    (*p)[0] = 0;
	    (*p)[1] = 0.;
	    (*p)[2] = 4. * lambda[3];
	    (*p)[3] = 4. * lambda[2];
	    return (FLOAT *)values;
    }

    return (FLOAT *)values;
}

static void
P2_interp(DOF *dof, ELEMENT *e, FLOAT **old_data, FLOAT **new_data)
/* interpolation */
{
    FLOAT *vold0 = old_data[0], *vold1 = old_data[1],
	  *eold0 = old_data[NVert + 0], *eold1 = old_data[NVert + 1],
	  *eold2 = old_data[NVert + 2], *eold3 = old_data[NVert + 3],
	  *eold4 = old_data[NVert + 4],
	  *vnew = new_data[0], *enew0 = new_data[1], *enew1 = new_data[2],
	  *enew2 = new_data[3], *enew3 = new_data[4],
	  v0, v1, e0, tmp;
    int i;

    for (i = 0; i < dof->dim; i++) {
	v0 = .125 * (*(vold0++));
	v1 = .125 * (*(vold1++));
	/* the new vertex */
	*vnew = *eold0;
	/* the two cut edges */
	e0 = .75 * (*(vnew++));
	*(enew0++) = 3. * v0 + e0 - v1;
	*(enew1++) = 3. * v1 + e0 - v0;
	/* the new edge on e->verts[2] */
	tmp = .25 * (*(eold0++)) - v0 - v1;
	*(enew2++) = tmp + .5 * ((*(eold1++)) + (*(eold3++)));
	/* the new edge on e->verts[3] */
	*(enew3++) = tmp + .5 * ((*(eold2++)) + (*(eold4++)));
    }
}

/*--------------------------- P3 element ----------------------------*/

static const FLOAT *
P3_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of basis functions */
{
    static FLOAT *buffer = NULL;
    FLOAT L0, L1, L2, L3, L03, L13, L23, L33, d, *p, *values;

    Unused(dof);
    Unused(e);

    if (buffer == NULL) {
#if USE_OMP
#pragma omp critical (P3_bas)
#endif	/* USE_OMP */
	if (buffer == NULL) {
	    buffer = phgAlloc(phgMaxThreads * (NVert + 2 * NEdge + NFace)
					    * sizeof(*buffer));
	    FreeAtExitNoCheck(buffer);
	}
    }

    p = values = buffer + phgThreadId * (NVert + 2 * NEdge + NFace);

    if (no1 <= 0)
	no1 = 20;
    assert(no0 < no1);

    if (no0 == 0 && no1 == 20) {
	L0 = lambda[0], L1 = lambda[1], L2 = lambda[2], L3 = lambda[3];
	L03 = 3. * L0, L13 = 3. * L1, L23 = 3. * L2, L33 = 3. * L3;
	*(p++) = .5 * (L03 - 1.) * (L03 - 2.) * L0;	/* vertex 0 */
	*(p++) = .5 * (L13 - 1.) * (L13 - 2.) * L1;	/* vertex 1 */
	*(p++) = .5 * (L23 - 1.) * (L23 - 2.) * L2;	/* vertex 2 */
	*(p++) = .5 * (L33 - 1.) * (L33 - 2.) * L3;	/* vertex 3 */
	*(p++) = 4.5 * (L03 - 1.) * (d = L0 * L1);	/* edge 0-1 */
	*(p++) = 4.5 * (L13 - 1.) * d;		/* edge 0-1 */
	*(p++) = 4.5 * (L03 - 1.) * (d = L0 * L2);	/* edge 0-2 */
	*(p++) = 4.5 * (L23 - 1.) * d;		/* edge 0-2 */
	*(p++) = 4.5 * (L03 - 1.) * (d = L0 * L3);	/* edge 0-3 */
	*(p++) = 4.5 * (L33 - 1.) * d;		/* edge 0-3 */
	*(p++) = 4.5 * (L13 - 1.) * (d = L1 * L2);	/* edge 1-2 */
	*(p++) = 4.5 * (L23 - 1.) * d;		/* edge 1-2 */
	*(p++) = 4.5 * (L13 - 1.) * (d = L1 * L3);	/* edge 1-3 */
	*(p++) = 4.5 * (L33 - 1.) * d;		/* edge 1-3 */
	*(p++) = 4.5 * (L23 - 1.) * (d = L2 * L3);	/* edge 2-3 */
	*(p++) = 4.5 * (L33 - 1.) * d;		/* edge 2-3 */
	*(p++) = 27. * L1 * L2 * L3;		/* face 0 */
	*(p++) = 27. * L0 * L2 * L3;		/* face 1 */
	*(p++) = 27. * L0 * L1 * L3;		/* face 2 */
	*p     = 27. * L0 * L1 * L2;		/* face 3 */
	return values;
    }

    switch (no0) {
	case 0:
	    L0 = lambda[0];
	    L03 = 3. * L0;
	    *(p++) = .5 * (L03 - 1.) * (L03 - 2.) * L0;
	    if (no1 == 1)
		return values;
	case 1:
	    L1 = lambda[1];
	    L13 = 3. * L1;
	    *(p++) = .5 * (L13 - 1.) * (L13 - 2.) * L1;
	    if (no1 == 2)
		return values;
	case 2:
	    L2 = lambda[2];
	    L23 = 3. * L2;
	    *(p++) = .5 * (L23 - 1.) * (L23 - 2.) * L2;
	    if (no1 == 3)
		return values;
	case 3:
	    L3 = lambda[3];
	    L33 = 3. * L3;
	    *(p++) = .5 * (L33 - 1.) * (L33 - 2.) * L3;
	    if (no1 == 4)
		return values;
	case 4:
	    *(p++) = (13.5 * lambda[0] - 4.5) * lambda[0] * lambda[1];
	    if (no1 == 5)
		return values;
	case 5:
	    *(p++) = (13.5 * lambda[1] - 4.5) * lambda[0] * lambda[1];
	    if (no1 == 6)
		return values;
	case 6:
	    *(p++) = (13.5 * lambda[0] - 4.5) * lambda[0] * lambda[2];
	    if (no1 == 7)
		return values;
	case 7:
	    *(p++) = (13.5 * lambda[2] - 4.5) * lambda[0] * lambda[2];
	    if (no1 == 8)
		return values;
	case 8:
	    *(p++) = (13.5 * lambda[0] - 4.5) * lambda[0] * lambda[3];
	    if (no1 == 9)
		return values;
	case 9:
	    *(p++) = (13.5 * lambda[3] - 4.5) * lambda[0] * lambda[3];
	    if (no1 == 10)
		return values;
	case 10:
	    *(p++) = (13.5 * lambda[1] - 4.5) * lambda[1] * lambda[2];
	    if (no1 == 11)
		return values;
	case 11:
	    *(p++) = (13.5 * lambda[2] - 4.5) * lambda[1] * lambda[2];
	    if (no1 == 12)
		return values;
	case 12:
	    *(p++) = (13.5 * lambda[1] - 4.5) * lambda[1] * lambda[3];
	    if (no1 == 13)
		return values;
	case 13:
	    *(p++) = (13.5 * lambda[3] - 4.5) * lambda[1] * lambda[3];
	    if (no1 == 14)
		return values;
	case 14:
	    *(p++) = (13.5 * lambda[2] - 4.5) * lambda[2] * lambda[3];
	    if (no1 == 15)
		return values;
	case 15:
	    *(p++) = (13.5 * lambda[3] - 4.5) * lambda[2] * lambda[3];
	    if (no1 == 16)
		return values;
	case 16:
	    *(p++) = 27. * lambda[1] * lambda[2] * lambda[3];
	    if (no1 == 17)
		return values;
	case 17:
	    *(p++) = 27. * lambda[0] * lambda[2] * lambda[3];
	    if (no1 == 18)
		return values;
	case 18:
	    *(p++) = 27. * lambda[0] * lambda[1] * lambda[3];
	    if (no1 == 19)
		return values;
	case 19:
	    *p = 27. * lambda[0] * lambda[1] * lambda[2];
	    return values;
    }

    return values;
}

static const FLOAT *
P3_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of gradients of basis functions */
{
    static FLOAT (*buffer)[Dim + 1] = NULL;
    FLOAT (*values)[Dim + 1];
    FLOAT (*p)[Dim + 1];

    Unused(dof);
    Unused(e);

    if (buffer == NULL) {
#if USE_OMP
#pragma omp critical (P3_bas)
#endif	/* USE_OMP */
	if (buffer == NULL) {
	    buffer = phgAlloc(phgMaxThreads * (NVert + 2 * NEdge + NFace)
					    * sizeof(*buffer));
	    FreeAtExitNoCheck(buffer);
	}
    }

    p = values = buffer + phgThreadId * (NVert + 2 * NEdge + NFace);

    if (no1 <= 0)
	no1 = 20;
    assert(no0 < no1);

    switch (no0) {
	case 0:	/* vertex 0 */
	    (*p)[0] = (13.5 * lambda[0] - 9.) * lambda[0] + 1.;
	    (*p)[1] = 0.;
	    (*p)[2] = 0.;
	    (*p)[3] = 0.;
	    if (no1 == 1)
		return (FLOAT *)values;
	    p++;
	case 1:	/* vertex 1 */
	    (*p)[0] = 0.;
	    (*p)[1] = (13.5 * lambda[1] - 9.) * lambda[1] + 1.;
	    (*p)[2] = 0.;
	    (*p)[3] = 0.;
	    if (no1 == 2)
		return (FLOAT *)values;
	    p++;
	case 2:	/* vertex 2 */
	    (*p)[0] = 0.;
	    (*p)[1] = 0.;
	    (*p)[2] = (13.5 * lambda[2] - 9.) * lambda[2] + 1.;
	    (*p)[3] = 0.;
	    if (no1 == 3)
		return (FLOAT *)values;
	    p++;
	case 3:	/* vertex 3 */
	    (*p)[0] = 0.;
	    (*p)[1] = 0.;
	    (*p)[2] = 0.;
	    (*p)[3] = (13.5 * lambda[3] - 9.) * lambda[3] + 1.;
	    if (no1 == 4)
		return (FLOAT *)values;
	    p++;
	case 4:	/* edge 0-1 */
	    (*p)[0] = (27. * lambda[0] - 4.5) * lambda[1];
	    (*p)[1] = 4.5 * lambda[0] * (3. * lambda[0] - 1.);
	    (*p)[2] = 0.;
	    (*p)[3] = 0.;
	    if (no1 == 5)
		return (FLOAT *)values;
	    p++;
	case 5:	/* edge 0-1 */
	    (*p)[0] = 4.5 * lambda[1] * (3. * lambda[1] - 1.);
	    (*p)[1] = (27. * lambda[1] - 4.5) * lambda[0];
	    (*p)[2] = 0.;
	    (*p)[3] = 0.;
	    if (no1 == 6)
		return (FLOAT *)values;
	    p++;
	case 6:	/* edge 0-2 */
	    (*p)[0] = (27. * lambda[0] - 4.5) * lambda[2];
	    (*p)[1] = 0.;
	    (*p)[2] = 4.5 * lambda[0] * (3. * lambda[0] - 1.);
	    (*p)[3] = 0.;
	    if (no1 == 7)
		return (FLOAT *)values;
	    p++;
	case 7:	/* edge 0-2 */
	    (*p)[0] = 4.5 * lambda[2] * (3. * lambda[2] - 1.);
	    (*p)[1] = 0.;
	    (*p)[2] = (27. * lambda[2] - 4.5) * lambda[0];
	    (*p)[3] = 0.;
	    if (no1 == 8)
		return (FLOAT *)values;
	    p++;
	case 8:	/* edge 0-3 */
	    (*p)[0] = (27. * lambda[0] - 4.5) * lambda[3];
	    (*p)[1] = 0.;
	    (*p)[2] = 0.;
	    (*p)[3] = 4.5 * lambda[0] * (3. * lambda[0] - 1.);
	    if (no1 == 9)
		return (FLOAT *)values;
	    p++;
	case 9:	/* edge 0-3 */
	    (*p)[0] = 4.5 * lambda[3] * (3. * lambda[3] - 1.);
	    (*p)[1] = 0.;
	    (*p)[2] = 0.;
	    (*p)[3] = (27. * lambda[3] - 4.5) * lambda[0];
	    if (no1 == 10)
		return (FLOAT *)values;
	    p++;
	case 10:	/* edge 1-2 */
	    (*p)[0] = 0.;
	    (*p)[1] = (27. * lambda[1] - 4.5) * lambda[2];
	    (*p)[2] = 4.5 * lambda[1] * (3. * lambda[1] - 1.);
	    (*p)[3] = 0.;
	    if (no1 == 11)
		return (FLOAT *)values;
	    p++;
	case 11:	/* edge 1-2 */
	    (*p)[0] = 0.;
	    (*p)[1] = 4.5 * lambda[2] * (3. * lambda[2] - 1.);
	    (*p)[2] = (27. * lambda[2] - 4.5) * lambda[1];
	    (*p)[3] = 0.;
	    if (no1 == 12)
		return (FLOAT *)values;
	    p++;
	case 12:	/* edge 1-3 */
	    (*p)[0] = 0.;
	    (*p)[1] = (27. * lambda[1] - 4.5) * lambda[3];
	    (*p)[2] = 0.;
	    (*p)[3] = 4.5 * lambda[1] * (3. * lambda[1] - 1.);
	    if (no1 == 13)
		return (FLOAT *)values;
	    p++;
	case 13:	/* edge 2-3 */
	    (*p)[0] = 0.;
	    (*p)[1] = 4.5 * lambda[3] * (3. * lambda[3] - 1.);
	    (*p)[2] = 0.;
	    (*p)[3] = (27. * lambda[3] - 4.5) * lambda[1];
	    if (no1 == 14)
		return (FLOAT *)values;
	    p++;
	case 14:	/* edge 2-3 */
	    (*p)[0] = 0.;
	    (*p)[1] = 0.;
	    (*p)[2] = (27. * lambda[2] - 4.5) * lambda[3];
	    (*p)[3] = 4.5 * lambda[2] * (3. * lambda[2] - 1.);
	    if (no1 == 15)
		return (FLOAT *)values;
	    p++;
	case 15:	/* edge 2-3 */
	    (*p)[0] = 0.;
	    (*p)[1] = 0.;
	    (*p)[2] = 4.5 * lambda[3] * (3. * lambda[3] - 1.);
	    (*p)[3] = (27. * lambda[3] - 4.5) * lambda[2];
	    if (no1 == 16)
		return (FLOAT *)values;
	    p++;
	case 16:	/* face 0 */
	    (*p)[0] = 0.;
	    (*p)[1] = 27. * lambda[2] * lambda[3];
	    (*p)[2] = 27. * lambda[1] * lambda[3];
	    (*p)[3] = 27. * lambda[1] * lambda[2];
	    if (no1 == 17)
		return (FLOAT *)values;
	    p++;
	case 17:	/* face 1 */
	    (*p)[0] = 27. * lambda[2] * lambda[3];
	    (*p)[1] = 0.;
	    (*p)[2] = 27. * lambda[0] * lambda[3];
	    (*p)[3] = 27. * lambda[0] * lambda[2];
	    if (no1 == 18)
		return (FLOAT *)values;
	    p++;
	case 18:	/* face 2 */
	    (*p)[0] = 27. * lambda[1] * lambda[3];
	    (*p)[1] = 27. * lambda[0] * lambda[3];
	    (*p)[2] = 0.;
	    (*p)[3] = 27. * lambda[0] * lambda[1];
	    if (no1 == 19)
		return (FLOAT *)values;
	    p++;
	case 19:	/* face 3 */
	    (*p)[0] = 27. * lambda[1] * lambda[2];
	    (*p)[1] = 27. * lambda[0] * lambda[2];
	    (*p)[2] = 27. * lambda[0] * lambda[1];
	    (*p)[3] = 0.;
	    return (FLOAT *)values;
    }

    return (FLOAT *)values;
}

/*--------------------------- P4 element ----------------------------*/

static const FLOAT *
P4_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of basis functions */
{
    static FLOAT *buffer = NULL;
    FLOAT L0, L1, L2, L3, L0T4M1, L1T4M1, L2T4M1, L3T4M1;
    FLOAT t0, t1, t2, t3, d01, d02, d03, d12, d13, d23;
    FLOAT d, d123, d023, d013, d012;
    FLOAT *p, *values;

    Unused(dof);
    Unused(e);

    if (buffer == NULL) {
#if USE_OMP
#pragma omp critical (P4_bas)
#endif	/* USE_OMP */
	if (buffer == NULL) {
	    buffer = phgAlloc(phgMaxThreads * (NVert + 3 * (NEdge + NFace) + 1)
					    * sizeof(*buffer));
	    FreeAtExitNoCheck(buffer);
	}
    }

    p = values = buffer + phgThreadId * (NVert + 3 * (NEdge + NFace) + 1);

    if (no1 <= 0)
	no1 = 35;
    assert(no0 < no1);

    if (no0 == 0 && no1 == 35) {
	L0 = lambda[0], L1 = lambda[1], L2 = lambda[2], L3 = lambda[3];
	L0T4M1 = 4. * L0 - 1.;
	L1T4M1 = 4. * L1 - 1.;
	L2T4M1 = 4. * L2 - 1.;
	L3T4M1 = 4. * L3 - 1.;

	/* vertices */
	*(p++) = L0T4M1*((2./(FLOAT)3.)*L0 - (1./(FLOAT)3.))*(L0T4M1 - 2.)*L0;
	*(p++) = L1T4M1*((2./(FLOAT)3.)*L1 - (1./(FLOAT)3.))*(L1T4M1 - 2.)*L1;
	*(p++) = L2T4M1*((2./(FLOAT)3.)*L2 - (1./(FLOAT)3.))*(L2T4M1 - 2.)*L2;
	*(p++) = L3T4M1*((2./(FLOAT)3.)*L3 - (1./(FLOAT)3.))*(L3T4M1 - 2.)*L3;

	t0 = L0T4M1 * ((32./(FLOAT)3.) * L0 - (16./(FLOAT)3.));
	t1 = L1T4M1 * ((32./(FLOAT)3.) * L1 - (16./(FLOAT)3.));
	t2 = L2T4M1 * ((32./(FLOAT)3.) * L2 - (16./(FLOAT)3.));
	t3 = L3T4M1 * ((32./(FLOAT)3.) * L3 - (16./(FLOAT)3.));

	d01 = L0 * L1;
	d02 = L0 * L2;
	d03 = L0 * L3;
	d12 = L1 * L2;
	d13 = L1 * L3;
	d23 = L2 * L3;

	/* edge 0-1 */
	*(p++) = t0 * d01;
	*(p++) = 4. * L0T4M1 * L1T4M1 * d01;
	*(p++) = t1 * d01;

	/* edge 0-2 */
	*(p++) = t0 * d02;
	*(p++) = 4. * L0T4M1 * L2T4M1 * d02;
	*(p++) = t2 * d02;

	/* edge 0-3 */
	*(p++) = t0 * d03;
	*(p++) = 4. * L0T4M1 * L3T4M1 * d03;
	*(p++) = t3 * d03;

	/* edge 1-2 */
	*(p++) = t1 * d12;
	*(p++) = 4. * L1T4M1 * L2T4M1 * d12;
	*(p++) = t2 * d12;

	/* edge 1-3 */
	*(p++) = t1 * d13;
	*(p++) = 4. * L1T4M1 * L3T4M1 * d13;
	*(p++) = t3 * d13;

	/* edge 2-3 */
	*(p++) = t2 * d23;
	*(p++) = 4. * L2T4M1 * L3T4M1 * d23;
	*(p++) = t3 * d23;

	L0T4M1 *= 32.;
	L1T4M1 *= 32.;
	L2T4M1 *= 32.;
	L3T4M1 *= 32.;

	d123 = L1 * d23;
	d023 = L0 * d23;
	d013 = L0 * d13;
	d012 = L0 * d12;

	/* face 0 */
	*(p++) = L1T4M1 * d123;
	*(p++) = L2T4M1 * d123;
	*(p++) = L3T4M1 * d123;

	/* face 1 */
	*(p++) = L0T4M1 * d023;
	*(p++) = L2T4M1 * d023;
	*(p++) = L3T4M1 * d023;

	/* face 2 */
	*(p++) = L0T4M1 * d013;
	*(p++) = L1T4M1 * d013;
	*(p++) = L3T4M1 * d013;

	/* face 3 */
	*(p++) = L0T4M1 * d012;
	*(p++) = L1T4M1 * d012;
	*(p++) = L2T4M1 * d012;

	/* element */
	*p     = 256. * L0 * d123;

	return values;
    }

    switch (no0) {
	case 0:
	    d = 4. * (L0 = lambda[0]);
	    *(p++) = (d-1.)*((2./(FLOAT)3.)*L0 - (1./(FLOAT)3.))*(d-3.)*L0;
	    if (no1 == 1)
		return values;
	case 1:
	    d = 4. * (L1 = lambda[1]);
	    *(p++) = (d-1.)*((2./(FLOAT)3.)*L1 - (1./(FLOAT)3.))*(d-3.)*L1;
	    if (no1 == 2)
		return values;
	case 2:
	    d = 4. * (L2 = lambda[2]);
	    *(p++) = (d-1.)*((2./(FLOAT)3.)*L2 - (1./(FLOAT)3.))*(d-3.)*L2;
	    if (no1 == 3)
		return values;
	case 3:
	    d = 4. * (L3 = lambda[3]);
	    *(p++) = (d-1.)*((2./(FLOAT)3.)*L3 - (1./(FLOAT)3.))*(d-3.)*L3;
	    if (no1 == 4)
		return values;

	case 4:
	    L0 = lambda[0];
	    L1 = lambda[1];
	    *(p++) = ((64./(FLOAT)3.)*L0 - (16./(FLOAT)3.))*(L0+L0-1.)*L0*L1;
	    if (no1 == 5)
		return values;
	case 5:
	    L0 = lambda[0];
	    L1 = lambda[1];
	    *(p++) = (16. * L0 - 4.) * (4. * L1 - 1.) * L0 * L1;
	    if (no1 == 6)
		return values;
	case 6:
	    L0 = lambda[0];
	    L1 = lambda[1];
	    *(p++) = ((64./(FLOAT)3.)*L1 - (16./(FLOAT)3.))*(L1+L1-1.)*L0*L1;
	    if (no1 == 7)
		return values;

	case 7:
	    L0 = lambda[0];
	    L2 = lambda[2];
	    *(p++) = ((64./(FLOAT)3.)*L0 - (16./(FLOAT)3.))*(L0+L0-1.)*L0*L2;
	    if (no1 == 8)
		return values;
	case 8:
	    L0 = lambda[0];
	    L2 = lambda[2];
	    *(p++) = (16. * L0 - 4.) * (4. * L2 - 1.) * L0 * L2;
	    if (no1 == 9)
		return values;
	case 9:
	    L0 = lambda[0];
	    L2 = lambda[2];
	    *(p++) = ((64./(FLOAT)3.)*L2 - (16./(FLOAT)3.))*(L2+L2-1.)*L0*L2;
	    if (no1 == 10)
		return values;

	case 10:
	    L0 = lambda[0];
	    L3 = lambda[3];
	    *(p++) = ((64./(FLOAT)3.)*L0 - (16./(FLOAT)3.))*(L0+L0-1.)*L0*L3;
	    if (no1 == 11)
		return values;
	case 11:
	    L0 = lambda[0];
	    L3 = lambda[3];
	    *(p++) = (16. * L0 - 4.) * (4. * L3 - 1.) * L0 * L3;
	    if (no1 == 12)
		return values;
	case 12:
	    L0 = lambda[0];
	    L3 = lambda[3];
	    *(p++) = ((64./(FLOAT)3.)*L3 - (16./(FLOAT)3.))*(L3+L3-1.)*L0*L3;
	    if (no1 == 13)
		return values;

	case 13:
	    L1 = lambda[1];
	    L2 = lambda[2];
	    *(p++) = ((64./(FLOAT)3.)*L1 - (16./(FLOAT)3.))*(L1+L1-1.)*L1*L2;
	    if (no1 == 14)
		return values;
	case 14:
	    L1 = lambda[1];
	    L2 = lambda[2];
	    *(p++) = (16. * L1 - 4.) * (4. * L2 - 1.) * L1 * L2;
	    if (no1 == 15)
		return values;
	case 15:
	    L1 = lambda[1];
	    L2 = lambda[2];
	    *(p++) = ((64./(FLOAT)3.)*L2 - (16./(FLOAT)3.))*(L2+L2-1.)*L1*L2;
	    if (no1 == 16)
		return values;

	case 16:
	    L1 = lambda[1];
	    L3 = lambda[3];
	    *(p++) = ((64./(FLOAT)3.)*L1 - (16./(FLOAT)3.))*(L1+L1-1.)*L1*L3;
	    if (no1 == 17)
		return values;
	case 17:
	    L1 = lambda[1];
	    L3 = lambda[3];
	    *(p++) = (16. * L1 - 4.) * (4. * L3 - 1.) * L1 * L3;
	    if (no1 == 18)
		return values;
	case 18:
	    L1 = lambda[1];
	    L3 = lambda[3];
	    *(p++) = ((64./(FLOAT)3.)*L3 - (16./(FLOAT)3.))*(L3+L3-1.)*L1*L3;
	    if (no1 == 19)
		return values;

	case 19:
	    L2 = lambda[2];
	    L3 = lambda[3];
	    *(p++) = ((64./(FLOAT)3.)*L2 - (16./(FLOAT)3.))*(L2+L2-1.)*L2*L3;
	    if (no1 == 20)
		return values;
	case 20:
	    L2 = lambda[2];
	    L3 = lambda[3];
	    *(p++) = (16. * L2 - 4.) * (4. * L3 - 1.) * L2 * L3;
	    if (no1 == 21)
		return values;
	case 21:
	    L2 = lambda[2];
	    L3 = lambda[3];
	    *(p++) = ((64./(FLOAT)3.)*L3 - (16./(FLOAT)3.))*(L3+L3-1.)*L2*L3;
	    if (no1 == 22)
		return values;

	case 22:
	    L1 = lambda[1];
	    L2 = lambda[2];
	    L3 = lambda[3];
	    *(p++) = (128. * L1 - 32.) * L1 * L2 * L3;
	    if (no1 == 23)
		return values;
	case 23:
	    L1 = lambda[1];
	    L2 = lambda[2];
	    L3 = lambda[3];
	    *(p++) = (128. * L2 - 32.) * L1 * L2 * L3;
	    if (no1 == 24)
		return values;
	case 24:
	    L1 = lambda[1];
	    L2 = lambda[2];
	    L3 = lambda[3];
	    *(p++) = (128. * L3 - 32.) * L1 * L2 * L3;
	    if (no1 == 25)
		return values;

	case 25:
	    L0 = lambda[0];
	    L2 = lambda[2];
	    L3 = lambda[3];
	    *(p++) = (128. * L0 - 32.) * L0 * L2 * L3;
	    if (no1 == 26)
		return values;
	case 26:
	    L0 = lambda[0];
	    L2 = lambda[2];
	    L3 = lambda[3];
	    *(p++) = (128. * L2 - 32.) * L0 * L2 * L3;
	    if (no1 == 27)
		return values;
	case 27:
	    L0 = lambda[0];
	    L2 = lambda[2];
	    L3 = lambda[3];
	    *(p++) = (128. * L3 - 32.) * L0 * L2 * L3;
	    if (no1 == 28)
		return values;

	case 28:
	    L0 = lambda[0];
	    L1 = lambda[1];
	    L3 = lambda[3];
	    *(p++) = (128. * L0 - 32.) * L0 * L1 * L3;
	    if (no1 == 29)
		return values;
	case 29:
	    L0 = lambda[0];
	    L1 = lambda[1];
	    L3 = lambda[3];
	    *(p++) = (128. * L1 - 32.) * L0 * L1 * L3;
	    if (no1 == 30)
		return values;
	case 30:
	    L0 = lambda[0];
	    L1 = lambda[1];
	    L3 = lambda[3];
	    *(p++) = (128. * L3 - 32.) * L0 * L1 * L3;
	    if (no1 == 31)
		return values;

	case 31:
	    L0 = lambda[0];
	    L1 = lambda[1];
	    L2 = lambda[2];
	    *(p++) = (128. * L0 - 32.) * L0 * L1 * L2;
	    if (no1 == 32)
		return values;
	case 32:
	    L0 = lambda[0];
	    L1 = lambda[1];
	    L2 = lambda[2];
	    *(p++) = (128. * L1 - 32.) * L0 * L1 * L2;
	    if (no1 == 33)
		return values;
	case 33:
	    L0 = lambda[0];
	    L1 = lambda[1];
	    L2 = lambda[2];
	    *(p++) = (128. * L2 - 32.) * L0 * L1 * L2;
	    if (no1 == 34)
		return values;

	case 34:
	    *p     = 256. * lambda[0] * lambda[1] * lambda[2] * lambda[3];
	    return values;
    }

    return values;
}

static const FLOAT *
P4_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
/* evaluation of gradients of basis functions */
{
    static FLOAT (*buffer)[Dim + 1] = NULL;
    FLOAT L0, L1, L2, L3, (*p)[Dim + 1], (*values)[Dim + 1];

    Unused(dof);
    Unused(e);

    if (buffer == NULL) {
#if USE_OMP
#pragma omp critical (P4_bas)
#endif	/* USE_OMP */
	if (buffer == NULL) {
	    buffer = phgAlloc(phgMaxThreads * (NVert + 3 * (NEdge + NFace) + 1)
					    * sizeof(*buffer));
	    FreeAtExitNoCheck(buffer);
	}
    }

    p = values = buffer + phgThreadId * (NVert + 3 * (NEdge + NFace) + 1);

    if (no1 <= 0)
	no1 = 35;
    assert(no0 < no1);

#define TEMP(k, v, u1, u2, u3) \
   if (no0 <= k) { \
	L0 = lambda[v]; \
	(*p)[v] = ((8./(FLOAT)3.)*L0 - 1.) * (L0*(16.*L0 - 12.) + 1.); \
	(*p)[u1] = 0.; \
	(*p)[u2] = 0.; \
	(*p)[u3] = 0.; \
	if (no1 == k + 1)	\
	    return (FLOAT *)values; \
	p++;	\
    }

    TEMP(0, 0, 1, 2, 3)		/* vertex 0 */
    TEMP(1, 1, 0, 2, 3)		/* vertex 1 */
    TEMP(2, 2, 0, 1, 3)		/* vertex 2 */
    TEMP(3, 3, 0, 1, 2)		/* vertex 3 */

#undef TEMP

#define TEMP0(k, v0, v1, u0, u1) \
    if (no0 <= k) { \
	L0 = lambda[v0]; \
	L1 = lambda[v1]; \
	(*p)[v0] = (L0*((384./(FLOAT)3.)*L0 - 64.) + (16./(FLOAT)3.))*L1; \
	(*p)[v1] = L0*((64./(FLOAT)3.)*L0 - (16./(FLOAT)3.))*(2. * L0 - 1.); \
	(*p)[u0] = 0.; \
	(*p)[u1] = 0.; \
	if (no1 == k + 1) \
	    return (FLOAT *)values; \
	p++;	\
    }
#define TEMP1(k, v0, v1, u0, u1) \
    if (no0 <= k) { \
	L0 = lambda[v0]; \
	L1 = lambda[v1]; \
	(*p)[v0] = (32. * L0 - 4.) * (4. * L1 - 1.) * L1; \
	(*p)[v1] = L0 * (4. * L0 - 1.) * (32. * L1 - 4.); \
	(*p)[u0] = 0.; \
	(*p)[u1] = 0.; \
	if (no1 == k + 1) \
	    return (FLOAT *)values; \
	p++;	\
    }

    TEMP0(4, 0, 1, 2, 3)	/* edge 0-1 point 0 */
    TEMP1(5, 0, 1, 2, 3)	/* edge 0-1 point 1 */
    TEMP0(6, 1, 0, 2, 3)	/* edge 0-1 point 2 */

    TEMP0(7, 0, 2, 1, 3)	/* edge 0-2 point 0 */
    TEMP1(8, 0, 2, 1, 3)	/* edge 0-2 point 1 */
    TEMP0(9, 2, 0, 1, 3)	/* edge 0-2 point 2 */

    TEMP0(10, 0, 3, 1, 2)	/* edge 0-3 point 0 */
    TEMP1(11, 0, 3, 1, 2)	/* edge 0-3 point 1 */
    TEMP0(12, 3, 0, 1, 2)	/* edge 0-3 point 2 */

    TEMP0(13, 1, 2, 0, 3)	/* edge 1-2 point 0 */
    TEMP1(14, 1, 2, 0, 3)	/* edge 1-2 point 1 */
    TEMP0(15, 2, 1, 0, 3)	/* edge 1-2 point 2 */

    TEMP0(16, 1, 3, 0, 2)	/* edge 1-3 point 0 */
    TEMP1(17, 1, 3, 0, 2)	/* edge 1-3 point 1 */
    TEMP0(18, 3, 1, 0, 2)	/* edge 1-3 point 2 */

    TEMP0(19, 2, 3, 0, 1)	/* edge 2-3 point 0 */
    TEMP1(20, 2, 3, 0, 1)	/* edge 2-3 point 1 */
    TEMP0(21, 3, 2, 0, 1)	/* edge 2-3 point 2 */

#undef TEMP0
#undef TEMP1

#define TEMP(k, v1, v2, v3, u) \
    if (no0 <= k) { \
	L1 = lambda[v1]; \
	L2 = lambda[v2]; \
	L3 = lambda[v3]; \
	(*p)[u] = 0.; \
	(*p)[v1] = (256. * L1 - 32.) * L2 * L3; \
	(*p)[v2] = (128. * L1 - 32.) * L1 * L3; \
	(*p)[v3] = (128. * L1 - 32.) * L1 * L2; \
	if (no1 == k + 1) \
	    return (FLOAT *)values; \
	p++;	\
    }

    TEMP(22, 1, 2, 3, 0)	/* face 0 point 0 */
    TEMP(23, 2, 3, 1, 0)	/* face 0 point 1 */
    TEMP(24, 3, 1, 2, 0)	/* face 0 point 2 */

    TEMP(25, 0, 2, 3, 1)	/* face 1 point 0 */
    TEMP(26, 2, 3, 0, 1)	/* face 1 point 1 */
    TEMP(27, 3, 0, 2, 1)	/* face 1 point 2 */

    TEMP(28, 0, 1, 3, 2)	/* face 2 point 0 */
    TEMP(29, 1, 3, 0, 2)	/* face 2 point 1 */
    TEMP(30, 3, 0, 1, 2)	/* face 2 point 2 */

    TEMP(31, 0, 1, 2, 3)	/* face 3 point 0 */
    TEMP(32, 1, 2, 0, 3)	/* face 3 point 1 */
    TEMP(33, 2, 0, 1, 3)	/* face 3 point 2 */

#undef TEMP

    (*p)[0] = 256. * lambda[1] * lambda[2] * lambda[3];
    (*p)[1] = 256. * lambda[0] * lambda[2] * lambda[3];
    (*p)[2] = 256. * lambda[0] * lambda[1] * lambda[3];
    (*p)[3] = 256. * lambda[0] * lambda[1] * lambda[2];

    return (FLOAT *)values;
}

/*------------------- General interpolation function --------------------*/
#if 0
#include <math.h>
static void
check_old_data(DOF *dof, ELEMENT *e, FLOAT **data)
{
    int i, j;
    FLOAT v, *p;
    COORD *p0, *p1, *p2, *p3;
    extern void f(FLOAT, FLOAT, FLOAT, FLOAT *); /* defined in dof_test.c */

    data += NVert;

    for (i = 0; i < NEdge; i++) {
	int M[dof->type->np_edge];
	p0 = dof->g->verts + e->verts[GetEdgeVertex(i, 0)];
	p1 = dof->g->verts + e->verts[GetEdgeVertex(i, 1)];
	phgDofMapEdgeData(dof->type, e, i, M);
	p = dof->type->points + dof->type->np_vert;
	for (j = 0; j < dof->type->np_edge; j++) {
	    f(  (*p0)[0] * p[0] + (*p1)[0] * p[1],
		(*p0)[1] * p[0] + (*p1)[1] * p[1],
		(*p0)[2] * p[0] + (*p1)[2] * p[1], &v);
	    if (Fabs(v - (*data)[M[j]]) > 1e-10)
		phgInfo(-1, "EDGE[%d][%d][%d] %lg <==> %lg!\n", e->index, i, M[j], (double)v, (double)(*data)[M[j]]);
	    p += 2;
	}
	data++;
    }

    for (i = 0; i < NFace; i++) {
	int M[dof->type->np_face];
	p0 = dof->g->verts + e->verts[GetFaceVertex(i, 0)];
	p1 = dof->g->verts + e->verts[GetFaceVertex(i, 1)];
	p2 = dof->g->verts + e->verts[GetFaceVertex(i, 2)];
	phgDofMapFaceData(dof->type, e, i, M);
	p = dof->type->points + dof->type->np_vert + dof->type->np_edge * 2;
	for (j = 0; j < dof->type->np_face; j++) {
	    f(  (*p0)[0] * p[0] + (*p1)[0] * p[1] + (*p2)[0] * p[2],
		(*p0)[1] * p[0] + (*p1)[1] * p[1] + (*p2)[1] * p[2],
		(*p0)[2] * p[0] + (*p1)[2] * p[1] + (*p2)[2] * p[2], &v);
	    if (Fabs(v - (*data)[M[j]]) > 1e-10)
		phgInfo(-1, "FACE[%d][%d][%d] %lg <==> %lg!\n", e->index, i, M[j], (double)v, (double)(*data)[M[j]]);
	    p += 3;
	}
	data++;
    }

    p0 = dof->g->verts + e->verts[0];
    p1 = dof->g->verts + e->verts[1];
    p2 = dof->g->verts + e->verts[2];
    p3 = dof->g->verts + e->verts[3];
    p = dof->type->points + dof->type->np_vert + dof->type->np_edge * 2
	+ dof->type->np_face * 3;
    for (j = 0; j < dof->type->np_elem; j++) {
	f((*p0)[0] * p[0] + (*p1)[0] * p[1] + (*p2)[0] * p[2] + (*p3)[0] * p[3],
	  (*p0)[1] * p[0] + (*p1)[1] * p[1] + (*p2)[1] * p[2] + (*p3)[1] * p[3],
	  (*p0)[2] * p[0] + (*p1)[2] * p[1] + (*p2)[2] * p[2] + (*p3)[2] * p[3],
	  &v);
	if (Fabs(v - (*data)[j]) > 1e-10)
	    phgInfo(-1, "ELEM[%d][%d] %lg <==> %lg!\n", e->index, j, (double)v, (double)(*data)[j]);
	p += 4;
    }
}
#endif

static void
Pn_interp(DOF *dof, ELEMENT *e, FLOAT **old_data, FLOAT **new_data)
/* universal interpolation function for scalar DOFs (slow).
 * TODO: cache values of basis functions to save computations */
{
    FLOAT lambda[] = {0.5, 0.5, 0., 0.};
    FLOAT *points = dof->type->points, *p, *d;
    int i, k, dim = dof->dim;
    INT V0 = 0, V1 = 0, V2 = 0, V3 = 0, V = 0;

#if 0
/*phgInfo(-1, "e = %p, type = %s, gen = %d\n", e,GTypeName(e->type),e->generation);*/
check_old_data(dof, e, old_data);
#endif

    if (dof->type->np_face > 1 || dof->type->np_edge > 0) {
	V = GlobalVertexP(dof->g, e->children[0]->verts[3]);
	V0 = GlobalVertexP(dof->g, e->verts[0]);
	V1 = GlobalVertexP(dof->g, e->verts[1]);
	V2 = GlobalVertexP(dof->g, e->verts[2]);
	V3 = GlobalVertexP(dof->g, e->verts[3]);
/*phgInfo(-1, "V V0 V1 V2 V3 = %d %d %d %d %d\n", V, V0, V1, V2, V3);*/
    }

    if (points == NULL)
	phgError(1, "%s:%d: unsupported DOF type \"%s\"!\n",
			__FILE__, __LINE__, dof->type->name);

    /* The new vertex */
    if (dof->type->np_vert > 0) {
	/* np_vert is only allowed to be 0 or 1 */
	phgDofEval_(dof, e, lambda, NULL, 0, new_data[0], old_data);
	points += dof->type->np_vert;
    }

    if (dof->type->np_edge > 0) {

	/* the cut edge 0 */
	p = points;
	d = new_data[1];
	k = dim;
	if (V0 > V) {
	    d += (dof->type->np_edge - 1) * k;
	    k = -k;
	}
	for (i = 0; i < dof->type->np_edge; i++) {
	    lambda[0] = p[0] + (lambda[1] = 0.5 * p[1]);
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 2;
	    d += k;
	}

	/* the cut edge 1 */
	p = points;
	d = new_data[2];
	k = dim;
	if (V > V1) {
	    d += (dof->type->np_edge - 1) * k;
	    k = -k;
	}
	for (i = 0; i < dof->type->np_edge; i++) {
	    lambda[1] = p[1] + (lambda[0] = 0.5 * p[0]);
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 2;
	    d += k;
	}

	/* the new edge on v2 */
	p = points;
	d = new_data[3];
	k = dim;
	if (V > V2) {
	    d += (dof->type->np_edge - 1) * k;
	    k = -k;
	}
	for (i = 0; i < dof->type->np_edge; i++) {
	    lambda[0] = lambda[1] = 0.5 * p[0];
	    lambda[2] = p[1];
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 2;
	    d += k;
	}

	/* the new edge on v3 */
	lambda[2] = 0.;
	p = points;
	d = new_data[4];
	k = dim;
	if (V > V3) {
	    d += (dof->type->np_edge - 1) * k;
	    k = -k;
	}
	for (i = 0; i < dof->type->np_edge; i++) {
	    lambda[0] = lambda[1] = 0.5 * p[0];
	    lambda[3] = p[1];
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 2;
	    d += k;
	}

	points = p;
    }

    if (dof->type->np_face == 1) {
	/* cut face F2 owning v0 */
	lambda[2] = 0.0;
	d = new_data[5];
	p = points;
	for (i = 0; i < dof->type->np_face; i++) {
	    lambda[0] = p[0] + (lambda[1] = 0.5 * p[1]);
	    lambda[3] = p[2];
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 3;
	    d += dim;
	}

	/* cut face F2 owning v1 */
	d = new_data[6];
	p = points;
	for (i = 0; i < dof->type->np_face; i++) {
	    lambda[1] = p[1] + (lambda[0] = 0.5 * p[0]);
	    lambda[3] = p[2];
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 3;
	    d += dim;
	}

	/* cut face F3 owning v0 */
	lambda[3] = 0.0;
	d = new_data[7];
	p = points;
	for (i = 0; i < dof->type->np_face; i++) {
	    lambda[0] = p[0] + (lambda[1] = 0.5 * p[1]);
	    lambda[2] = p[2];
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 3;
	    d += dim;
	}

	/* cut face F3 owning v1 */
	d = new_data[8];
	p = points;
	for (i = 0; i < dof->type->np_face; i++) {
	    lambda[1] = p[0] + (lambda[0] = 0.5 * p[1]);
	    lambda[2] = p[2];
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 3;
	    d += dim;
	}

	/* the new face */
	d = new_data[9];
	p = points;
	for (i = 0; i < dof->type->np_face; i++) {
	    lambda[0] = lambda[1] = 0.5 * p[0];
	    lambda[2] = p[1];
	    lambda[3] = p[2];
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 3;
	    d += dim;
	}

	points = p;
     }
     else if (dof->type->np_face == 3) {
	int M[3];

	/* cut face F2 owning v0 */
	lambda[2] = 0.0;
	d = new_data[5];
	p = points;
	M[0] = M[1] = M[2] = 0;
	(V0 > V3) ? M[0]++ : M[1]++;
	(V0 > V)  ? M[0]++ : M[2]++;
	(V3 > V)  ? M[1]++ : M[2]++;
	for (i = 0; i < 3; i++) {
	    lambda[0] = p[0] + (lambda[1] = 0.5 * p[2]);
	    lambda[3] = p[1];
	    phgDofEval_(dof, e, lambda, NULL, 0, d + M[i] * dim, old_data);
	    p += 3;
	}

	/* cut face F2 owning v1 */
	d = new_data[6];
	p = points;
	M[0] = M[1] = M[2] = 0;
	(V1 > V3) ? M[0]++ : M[1]++;
	(V1 > V)  ? M[0]++ : M[2]++;
	(V3 > V)  ? M[1]++ : M[2]++;
	for (i = 0; i < 3; i++) {
	    lambda[1] = p[0] + (lambda[0] = 0.5 * p[2]);
	    lambda[3] = p[1];
	    phgDofEval_(dof, e, lambda, NULL, 0, d + M[i] * dim, old_data);
	    p += 3;
	}

	/* cut face F3 owning v0 */
	lambda[3] = 0.0;
	d = new_data[7];
	p = points;
	M[0] = M[1] = M[2] = 0;
	(V0 > V2) ? M[0]++ : M[1]++;
	(V0 > V)  ? M[0]++ : M[2]++;
	(V2 > V)  ? M[1]++ : M[2]++;
	for (i = 0; i < 3; i++) {
	    lambda[0] = p[0] + (lambda[1] = 0.5 * p[2]);
	    lambda[2] = p[1];
	    phgDofEval_(dof, e, lambda, NULL, 0, d + M[i] * dim, old_data);
	    p += 3;
	}

	/* cut face F3 owning v1 */
	d = new_data[8];
	p = points;
	M[0] = M[1] = M[2] = 0;
	(V1 > V2) ? M[0]++ : M[1]++;
	(V1 > V)  ? M[0]++ : M[2]++;
	(V2 > V)  ? M[1]++ : M[2]++;
	for (i = 0; i < 3; i++) {
	    lambda[1] = p[0] + (lambda[0] = 0.5 * p[2]);
	    lambda[2] = p[1];
	    phgDofEval_(dof, e, lambda, NULL, 0, d + M[i] * dim, old_data);
	    p += 3;
	}

	/* the new face */
	d = new_data[9];
	p = points;
	M[0] = M[1] = M[2] = 0;
	(V2 > V3) ? M[0]++ : M[1]++;
	(V2 > V)  ? M[0]++ : M[2]++;
	(V3 > V)  ? M[1]++ : M[2]++;
	for (i = 0; i < 3; i++) {
	    lambda[0] = lambda[1] = 0.5 * p[2];
	    lambda[2] = p[0];
	    lambda[3] = p[1];
	    phgDofEval_(dof, e, lambda, NULL, 0, d + M[i] * dim, old_data);
	    p += 3;
	}

	points = p;
    }
    else if (dof->type->np_face != 0) {
	phgError(1, "%s:%d: unimplemented.\n", __FILE__, __LINE__);
    }

    /* volume data */
    if (dof->type->np_elem == 1) {
	int Vmap0[NVert], Vmap1[NVert];

	phgMapP2C(e, Vmap0, NULL, 0);	
	phgMapP2C(e, Vmap1, NULL, 1);	

	/* child 0 */
	d = new_data[10];
	p = points;
	for (i = 0; i < dof->type->np_elem; i++) {
	    lambda[0] = p[Vmap0[0]] + (lambda[1] = p[3] * 0.5);
	    lambda[2] = p[Vmap0[2]];
	    lambda[3] = p[Vmap0[3]];
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 4;
	    d += dim;
	}

	/* child 1 */
	d = new_data[11];
	p = points;
	for (i = 0; i < dof->type->np_elem; i++) {
	    lambda[1] = p[Vmap1[1]] + (lambda[0] = p[3] * 0.5);
	    lambda[2] = p[Vmap1[2]];
	    lambda[3] = p[Vmap1[3]];
	    phgDofEval_(dof, e, lambda, NULL, 0, d, old_data);
	    p += 4;
	    d += dim;
	}
    }
    else if (dof->type->np_elem != 0) {
	phgError(0, "%s:%d: TODO: permute lambda!\n", __FILE__, __LINE__);
    }
}

static void
P4_map(const DOF_TYPE *type, const ELEMENT *e, int face_no, int *M)
/* maps face data to global DOF for DOF_P4 */
{
    int v0, v1, v2;

    if (face_no < 0) {
	/* map element DOF */
	M[0] = 0;
	return;
    }

    M[0] = 0;
    M[1] = 0;
    M[2] = 0;
    v0 = MapVertex(e, GetFaceVertex(face_no, 0));
    v1 = MapVertex(e, GetFaceVertex(face_no, 1));
    v2 = MapVertex(e, GetFaceVertex(face_no, 2));
    (v0 > v1) ? M[0]++ : M[1]++;
    (v0 > v2) ? M[0]++ : M[2]++;
    (v1 > v2) ? M[1]++ : M[2]++;

    return;
}

/*-------------------------------------------------------------------*/

#else	/* !DEBUG_LAGRANGE1 */

static const FLOAT *
Pn_bas(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    return DOF_P5_.BasFuncs(dof, e, no0, no1, lambda);
}

static const FLOAT *
Pn_grad(DOF *dof, ELEMENT *e, int no0, int no1, const FLOAT *lambda)
{
    return DOF_P5_.BasGrads(dof, e, no0, no1, lambda);
}

static void
Pn_map(const DOF_TYPE *type, const ELEMENT *e, int no, int *M)
/* maps ordering of local face data to ordering of data in DOF */
{
    DOF_P5_.DofMap(type, e, no, M);
}

# define P1_interp	phgDofInterpC2FGeneric
# define P1_bas		Pn_bas
# define P1_grad	Pn_grad

# define P2_interp	phgDofInterpC2FGeneric
# define P2_bas		Pn_bas
# define P2_grad	Pn_grad

# define P3_bas		Pn_bas
# define P3_grad	Pn_grad

# define P4_bas		Pn_bas
# define P4_grad	Pn_grad
# define P4_map		Pn_map

# define Pn_interp	phgDofInterpC2FGeneric

#endif	/* !DEBUG_LAGRANGE1 */

static FLOAT P0_points[] = {_F(0.25), _F(0.25), _F(0.25), _F(0.25)};
DOF_TYPE DOF_P0_ = {DofReserved, "P0", P0_points, NULL, DOF_DG0, NULL, NULL,
			P0_interp, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
			P0_bas, P0_grad, NULL, FE_L2,
			TRUE,	/* is_nodal */
			TRUE, FALSE, -1,/* invariant, free_after_use, id */
			1, 0, 0, -1, 1,	/* nbas, order, D_order, cont, dim */
			0, 0, 0, 1};	/* np_vert, np_edge, np_face, np_elem */

static FLOAT P1_points[] = {_F(1.0)};
DOF_TYPE DOF_P1_ = {DofReserved, "P1", P1_points, NULL, DOF_DG0, NULL, NULL,
			P1_interp, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
			P1_bas, P1_grad, NULL, FE_H1,
			TRUE,	/* is_nodal */
			TRUE, FALSE, -1, NVert, 1, 0, 0, 1,
			1, 0, 0, 0};

static FLOAT P2_points[] = {_F(1.0), _F(0.5),_F(0.5)};
DOF_TYPE DOF_P2_ = {DofReserved, "P2", P2_points, NULL, DOF_DG1, NULL, NULL,
			P2_interp, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
			P2_bas, P2_grad, NULL, FE_H1,
			TRUE,	/* is_nodal */
			TRUE, FALSE, -1, NVert + NEdge, 2, 0, 0, 1,
			1, 1, 0, 0};

static FLOAT P3_points[] = {_F(1.0), _F(2.)/_F(3.),_F(1.)/_F(3.), _F(1.)/_F(3.),
	_F(2.)/_F(3.), _F(1.)/_F(3.),_F(1.)/_F(3.),_F(1.)/_F(3.)};
DOF_TYPE DOF_P3_ = {DofReserved, "P3", P3_points, NULL, DOF_DG2, NULL, NULL,
			Pn_interp, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
			P3_bas, P3_grad, NULL, FE_H1,
			TRUE,	/* is_nodal */
			TRUE, FALSE, -1, NVert + 2*NEdge + NFace, 3, 0, 0, 1,
			1, 2, 1, 0};

static FLOAT P4_points[] = {
    _F(1.),
    _F(.75),_F(.25),	_F(.5),_F(.5),  _F(.25),_F(.75),
    _F(.5),_F(.25),_F(.25),  _F(.25),_F(.5),_F(.25),  _F(.25),_F(.25),_F(.5),
    _F(.25),_F(.25),_F(.25),_F(.25)
};
DOF_TYPE DOF_P4_ = {DofReserved, "P4", P4_points, NULL, DOF_DG3, NULL, NULL,
			Pn_interp, phgDofInterpF2CGeneric, phgDofInitFuncPoint,
			P4_bas, P4_grad, P4_map, FE_H1,
			TRUE,	/* is_nodal */
			TRUE, FALSE, -1, NVert + 3*(NEdge + NFace) + 1,
			4, 0, 0, 1, 1, 3, 3, 1};

DOF_TYPE *DOF_Pn[] = {DOF_P0,  DOF_P1,  DOF_P2,  DOF_P3,  DOF_P4,
		      /* The following are defined in lagrange1.c */
		      DOF_P5,  DOF_P6,  DOF_P7,  DOF_P8,  DOF_P9,
		      DOF_P10, DOF_P11, DOF_P12, DOF_P13, DOF_P14,
		      DOF_P15, DOF_P16, NULL};
