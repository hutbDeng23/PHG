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

/* functions for geometrical computations.
 *
 * These functions defines a special DOF type, DOF_GEOM, which contains
 * geometrical data:
 *	Face data:	area, diameter, normal[Dim]
 *	Element data:	volume, diameter, jacobian[Dim + 1][Dim + 1]
 *
 * $Id: geom.c,v 1.81 2022/04/09 01:28:34 zlb Exp $
 */

#include "phg.h"

#include <math.h>
#include <string.h>


#if USE_LAPACK
/* Dense matrix operation. */
#  define DGETRF_FOR dgetrf_
#  define DGETRI_FOR dgetri_
void DGETRF_FOR(int *M, int *N, double *A, int *LDA, int *IPIV, int *INFO);
void DGETRI_FOR(int *N, double *A, int *LDA, int *IPIV, double *WORK, int *LWORK, int *INFO);
#endif

void
_phg_compute_elem_data(const GRID *g, const ELEMENT *e, FLOAT *vol, FLOAT *diam)
/* computes volume & diameter
 *-----------------------------------------------------------------------
 * General formula for the volume of a k-simplex in n-dimensions:
 *	V = \sqrt(|\det(W * W^T)|)/k!
 * where:
 *	W = [v1-v0,v2-v0,..,vk-v0]
 * For k=n=2:
 *		| x1-x0 y1-y0 |
 *	V = det | x2-x0 y2-y0 | / 2
 * For k=n=3:
 *		| x1-x0 y1-y0 z1-z0 |
 *	V = det | x2-x0 y2-y0 z2-z0 | / 6
 *		| x3-x0 y3-y0 z3-z0 |
 *
 * 3D: V=S*h/3 where S is the area of the base, h the height 
 *----------------------------------------------------------------------
 */
{
    static FLOAT oned6 = 1.0/(FLOAT)6.0;
    FLOAT x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;

    if (vol == NULL && diam == NULL)
	return;

    x0 = g->verts[e->verts[0]][0]; 
    y0 = g->verts[e->verts[0]][1];
    z0 = g->verts[e->verts[0]][2];
    x1 = g->verts[e->verts[1]][0] - x0;
    y1 = g->verts[e->verts[1]][1] - y0;
    z1 = g->verts[e->verts[1]][2] - z0;
    x2 = g->verts[e->verts[2]][0] - x0;
    y2 = g->verts[e->verts[2]][1] - y0;
    z2 = g->verts[e->verts[2]][2] - z0;
    x3 = g->verts[e->verts[3]][0] - x0;
    y3 = g->verts[e->verts[3]][1] - y0;
    z3 = g->verts[e->verts[3]][2] - z0;

    if (vol != NULL)
	*vol = Fabs(x1*y2*z3 + x2*y3*z1 + y1*z2*x3 -
		    (z1*y2*x3 + y1*x2*z3 + z2*y3*x1)) * oned6;

    if (diam != NULL) {
	/* cf. http://mathworld.wolfram.com/Circumsphere.html */
	FLOAT a, bx, by, bz, d1, d2, d3;
	a  = x1*y2*z3 + x2*y3*z1 + x3*y1*z2 - (x3*y2*z1 + x2*y1*z3 + x1*y3*z2);
	d1 = x1*x1 + y1*y1 + z1*z1;
	d2 = x2*x2 + y2*y2 + z2*z2;
	d3 = x3*x3 + y3*y3 + z3*z3;
	bx = d1*y2*z3 + d2*y3*z1 + d3*y1*z2 - (d3*y2*z1 + d2*y1*z3 + d1*y3*z2);
	by = d1*x2*z3 + d2*x3*z1 + d3*x1*z2 - (d3*x2*z1 + d2*x1*z3 + d1*x3*z2);
	bz = d1*x2*y3 + d2*x3*y1 + d3*x1*y2 - (d3*x2*y1 + d2*x1*y3 + d1*x3*y2);
	*diam = Sqrt(bx*bx + by*by + bz*bz) / (a == 0. ? 1. : Fabs(a));
    }
}

static void
compute_jacobian(const GRID *g, const ELEMENT *e, FLOAT *J)
/* computes the Jabobian
 *	J = D\lambda / Dx = [ D\lambda_0/Dx, D\lambda_0/Dy, D\lambda_0/Dz, c0;
 *			      D\lambda_1/Dx, D\lambda_1/Dy, D\lambda_1/Dz, c1;
 *			      D\lambda_2/Dx, D\lambda_2/Dy, D\lambda_2/Dz, c2;
 *			      D\lambda_3/Dx, D\lambda_3/Dy, D\lambda_3/Dz, c3 ]
 * of the element, by inverting the matrix
 *
 *	[x0, x1, x2, x3; y0, y1, y2, y3; z0, z1, z2, z3; 1 1 1 1].
 *
 * Note: lambda = J[0..Dim][0..Dim-1] * x + J[0..Dim][Dim]
 */
{
    int i, j, k;
    int n = Dim + 1;
    FLOAT d, t, a[Dim + 1][2 * (Dim + 1)];

    /* right hand side (I) */
    for (j = 0; j < n; j++) 
	for (i = 0; i < n; i++) 
	    a[j][n + i] = (i == j) ? 1.0 : 0.0;

    /* coefficient matrix */
    for (j = 0; j < n; j++) {
	for (i = 0; i < Dim; i++)
	    a[j][i] = g->verts[e->verts[j]][i];
	a[j][Dim] = 1.0;
    }

    /* solve linear systems of equations */
    for (j = 0; j < n; j++) {
	/* find pivot */
	d = Fabs(a[j][j]);
	k = j;
	for (i = j + 1; i < n; i++)
	    if (d < (t = Fabs(a[i][j]))) {
		d = t;
		k = i;
	    }
	if (k != j) /* swap rows j and k */
	    for (i = j; i < n + n; i++) {
		t = a[j][i];
		a[j][i] = a[k][i];
		a[k][i] = t;
	    }
	d = 1.0 / ((t = a[j][j]) == 0. ? 1. : t);
	for (i = j + 1; i < n + n; i++)
	    a[j][i] *= d;
	for (i = j + 1; i < n; i++) {
	    d = a[i][j];
	    for (k = j + 1; k < n + n; k++)
		a[i][k] -= d * a[j][k];
	}
    }

    for (j = n - 2; j >= 0; j--)
	for (i = j + 1; i < n; i++) {
	    d = a[j][i];
	    for (k = n; k < n + n; k++)
		a[j][k] -= d * a[i][k];
	}

    /* Jacobian = trans(a[0..3, 4..7]) */
    for (j = 0; j < n; j++)
	for (i = 0; i < n; i++)
	    *(J++) = a[i][j + n];
}

static void
jacobian_c2f(ELEMENT *e, const FLOAT J[][Dim + 1], FLOAT J0[][Dim + 1],
		FLOAT J1[][Dim + 1])
/* Evaluate (interpolate) the Jacobian 'J0' for child 0, 'J1' for child 1
   using the Jacobian 'J' of the parent.

\documentclass{article}
\pagestyle{empty}
\usepackage{amsmath}
\renewcommand\arraystretch{1.2}
\begin{document}
Let $\lambda_0$, $\lambda_1$, $\lambda_2$, and $\lambda_3$ be the
barycentric coordinates of the parent element, $J$ its Jacobian,
$\mu_{m_0}$, $\mu_{m_2}$, $\mu_{m_3}$, $\mu_3$ the barycentric
coordinates of child 0, then:
\[
 \begin{aligned}
  \begin{pmatrix}
    \mu_{m_0} \\ \mu_3 \\ \mu_{m_2} \\ \mu_{m_3} \\
  \end{pmatrix}
  &=
  \begin{pmatrix}
    1 & -1 & 0 & 0 \\
    0 & 2  & 0 & 0 \\
    0 & 0  & 1 & 0 \\
    0 & 0  & 0 & 1
  \end{pmatrix}
  \begin{pmatrix}
    \lambda_0 \\ \lambda_1 \\ \lambda_2 \\ \lambda_3
  \end{pmatrix}
  \\
  &=
  \begin{pmatrix}
    1 & -1 & 0 & 0 \\
    0 & 2  & 0 & 0 \\
    0 & 0  & 1 & 0 \\
    0 & 0  & 0 & 1
  \end{pmatrix}
  \begin{pmatrix}
    J_{00} & J_{01} & J_{02} & J_{03} \\
    J_{10} & J_{11} & J_{12} & J_{13} \\
    J_{20} & J_{21} & J_{22} & J_{23} \\
    J_{30} & J_{31} & J_{32} & J_{33}
  \end{pmatrix}
  \begin{pmatrix}
    x \\ y \\ z \\ 1
  \end{pmatrix}
  \\
  &=
  \begin{pmatrix}
    J_{00} - J_{10} & J_{01} - J_{11} & J_{02} - J_{12} & J_{03} - J_{13} \\
    2J_{10} & 2J_{11} & 2J_{12} & 2J_{13} \\
    J_{20} & J_{21} & J_{22} & J_{23} \\
    J_{30} & J_{31} & J_{32} & J_{33}
  \end{pmatrix}
  \begin{pmatrix}
    x \\ y \\ z \\ 1
  \end{pmatrix}
 \end{aligned}
\]
\end{document}
*/
{
    int i, map[NVert];

    if (J0 != NULL) {
	phgMapP2C(e, map, NULL, 0);
	for (i = 0; i < Dim + 1; i++) {
	    J0[map[0]][i] = J[0][i] - J[1][i];
	    J0[3][i]      = 2. * J[1][i];
	    J0[map[2]][i] = J[2][i];
	    J0[map[3]][i] = J[3][i];
	}
    }

    if (J1 != NULL) {
	phgMapP2C(e, map, NULL, 1);
	for (i = 0; i < Dim + 1; i++) {
	    J1[map[1]][i] = J[1][i] - J[0][i];
	    J1[3][i]      = 2. * J[0][i];
	    J1[map[2]][i] = J[2][i];
	    J1[map[3]][i] = J[3][i];
	}
    }
}

static void
jacobian_f2c(ELEMENT *e, FLOAT J[][Dim + 1], const FLOAT J0[][Dim + 1],
		const FLOAT J1[][Dim + 1])
/* Evaluate (interpolate) the Jacobian 'J' using either the Jacobian 'J0'
   of child 0 or 'J1' of child 1

\documentclass{article}
\pagestyle{empty}
\usepackage{amsmath}
\renewcommand\arraystretch{1.2}
\begin{document}
Let $\lambda_0$, $\lambda_1$, $\lambda_2$, and $\lambda_3$ be the
barycentric coordinates of the parent element,
$\mu_{m_0}$, $\mu_{m_2}$, $\mu_{m_3}$, $\mu_3$ the barycentric
coordinates of child 0, $J^0$ the Jacobian of child 0, then:
\[
 \begin{aligned}
  \begin{pmatrix}
    \lambda_0 \\ \lambda_1 \\ \lambda_2 \\ \lambda_3 \\
  \end{pmatrix}
  &=
  \begin{pmatrix}
    1 & \frac12 & 0 & 0 \\
    0 & \frac12  & 0 & 0 \\
    0 & 0  & 1 & 0 \\
    0 & 0  & 0 & 1
  \end{pmatrix}
  \begin{pmatrix}
    \mu_{m_0} \\ \mu_3 \\ \mu_{m_2} \\ \mu_{m_3}
  \end{pmatrix}
  \\
  &=
  \begin{pmatrix}
    1 & \frac12 & 0 & 0 \\
    0 & \frac12 & 0 & 0 \\
    0 & 0  & 1 & 0 \\
    0 & 0  & 0 & 1
  \end{pmatrix}
  \begin{pmatrix}
    J^0_{m_00} & J^0_{m_01} & J^0_{m_02} & J^0_{m_03} \\
    J^0_{30} & J^0_{31} & J^0_{32} & J^0_{33} \\
    J^0_{m_20} & J^0_{m_21} & J^0_{m_22} & J^0_{m_23} \\
    J^0_{m_30} & J^0_{m_31} & J^0_{m_32} & J^0_{m_33}
  \end{pmatrix}
  \begin{pmatrix}
    x \\ y \\ z \\ 1
  \end{pmatrix}
  \\
  &=
  \begin{pmatrix}
    J^0_{m_00} + \frac12 J^0_{30} & J^0_{m_01} + \frac12 J^0_{31}
	& J^0_{m_02} + \frac12 J^0_{32} & J^0_{m_03} + \frac12 J^0_{33} \\
    \frac12 J^0_{30} & \frac12 J^0_{31} & \frac12 J^0_{32} & \frac12 J^0_{33} \\
    J^0_{m_20} & J^0_{m_21} & J^0_{m_22} & J^0_{m_23} \\
    J^0_{m_30} & J^0_{m_31} & J^0_{m_32} & J^0_{m_33}
  \end{pmatrix}
  \begin{pmatrix}
    x \\ y \\ z \\ 1
  \end{pmatrix}
 \end{aligned}
\]
\end{document}
*/
{
    int i, map[NVert];

    assert(J0 != NULL || J1 != NULL);

    if (J0 != NULL) {
	phgMapP2C(e, map, NULL, 0);
	for (i = 0; i < Dim + 1; i++) {
	    J[0][i] = J0[map[0]][i] + (J[1][i] = .5 * J0[3][i]);
	    J[2][i] = J0[map[2]][i];
	    J[3][i] = J0[map[3]][i];
	}
	return;
    }

    /*if (J1 != NULL)*/ {
	phgMapP2C(e, map, NULL, 1);
	for (i = 0; i < Dim + 1; i++) {
	    J[1][i] = J1[map[1]][i] + (J[0][i] = .5 * J1[3][i]);
	    J[2][i] = J1[map[2]][i];
	    J[3][i] = J1[map[3]][i];
	}
	return;
    }
}

static void
compute_face_data(const GRID *g, const ELEMENT *e, int face,
			FLOAT *area, FLOAT *diam, FLOAT *n)
/* face normal ('n'), area and diameter */
{
    FLOAT d, a[3], b[3];
    int v0, v1, v2, v3;
    const FLOAT *p0, *p1, *p2;

    /* indices of the three vertices on the face */
    GetFaceVertices(e, face, v0, v1, v2, v3);
    p0 = g->verts[e->verts[v0]];
    p1 = g->verts[e->verts[v1]];
    p2 = g->verts[e->verts[v2]];

    a[0] = p1[0] - p0[0]; a[1] = p1[1] - p0[1]; a[2] = p1[2] - p0[2];
    b[0] = p2[0] - p0[0]; b[1] = p2[1] - p0[1]; b[2] = p2[2] - p0[2];

    if (n != NULL) {
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = a[2] * b[0] - a[0] * b[2];
	n[2] = a[0] * b[1] - a[1] * b[0];
	d = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
	d = 1.0 / (d == 0.0 ? 1.0 : Sqrt(d));
	n[0] *= d; n[1] *= d; n[2] *= d;
    }

    if (area == NULL && diam == NULL)
	return;

    /* compute the area (Heron's formula) and diameter
	(http://mathworld.wolfram.com/Circumradius.html)
	FIXME: better formula? */
    a[0] = Sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    a[1] = Sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    b[0] = p2[0] - p1[0]; b[1] = p2[1] - p1[1]; b[2] = p2[2] - p1[2];
    a[2] = Sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    d = (a[0] + a[1] + a[2]) * 0.5;

    if (area != NULL) {
	FLOAT c = d * (d - a[0]) * (d - a[1]) * (d - a[2]);
	if (c <= 0.) {
	    static BOOLEAN warned = FALSE;
	    if (!warned) {
		warned = TRUE;
		phgWarning("bad mesh: degenerated face found.\n");
	    }
	    *area = 0.;
	}
	else {
	    *area = Sqrt(c);
	}
    }

    if (diam != NULL) {
	d = d * (a[0]+a[1] - d) * (a[0]+a[2] - d) * (a[1]+a[2] - d);
	*diam = 0.5 * a[0] * a[1] * a[2] / (d <= 0. ? 1. : Sqrt(d));
    }
}

/*---------------- DOF for geometric data ------------------*/

static void
init(DOF *dof, ELEMENT *e, GTYPE type, int index,
	DOF_USER_FUNC userfunc, DOF_USER_FUNC_LAMBDA userfunc_lambda,
	const FLOAT *userdata, FLOAT *values, FLOAT **pvalues)
/* DOF function for geometric data */
{
    Unused(userfunc);
    Unused(userfunc_lambda);
    Unused(userdata);
    Unused(pvalues);

    switch (type) {
	case FACE:
	    compute_face_data(dof->g, e, index, values, values + 1, values + 2);
	    break;
	case VOLUME:
	    _phg_compute_elem_data(dof->g, e, values, values + 1);
	    compute_jacobian(dof->g, e, values + 2);
	    break;
    }
}

static void
adjust_direction(GRID *g, FLOAT *n, INT u0, INT u1, INT u2,
			INT v0, INT v1, INT v2)
/* adjusts the normal vector 'n' on face(u0,u1,u2) such that it matches
 * the direction of the face (v0,v1,v2), here u's and v's are local indices */
{
    INT v;
    int s0, s1;

    u0 = GlobalVertexP(g, u0);
    u1 = GlobalVertexP(g, u1);
    u2 = GlobalVertexP(g, u2);
    s0 = 1;
    if (u0 > u1) {s0 = -s0; v = u0; u0 = u1; u1 = v;}
    if (u0 > u2) {s0 = -s0; v = u0; u0 = u2; u2 = v;}
    if (u1 > u2) s0 = -s0;

    v0 = GlobalVertexP(g, v0);
    v1 = GlobalVertexP(g, v1);
    v2 = GlobalVertexP(g, v2);
    s1 = 1;
    if (v0 > v1) {s1 = -s1; v = v0; v0 = v1; v1 = v;}
    if (v0 > v2) {s1 = -s1; v = v0; v0 = v2; v2 = v;}
    if (v1 > v2) s1 = -s1;

    if (s0 != s1) {
	n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2];
    }
}

static void
interp(DOF *dof, ELEMENT *e, FLOAT **old_data, FLOAT **new_data)
/* interpolates some geometric data to save computations */
{
    int Vmap0[NVert], Vmap1[NVert];
    FLOAT *p0, *p1, *q;
    ELEMENT *e0 = e->children[0], *e1 = e->children[1];

    /*----- the face data, area, diameter, and normal vector -----*/
    phgMapP2C(e, Vmap0, NULL, 0);
    phgMapP2C(e, Vmap1, NULL, 1);

#if ALLOW_CURVED_BOUNDARY
    if (e->bound_func[0] != -1) {	/* curved boundary, recompute */
	/* Face 2 of e */
	p0 = new_data[5];
	p1 = new_data[6];
	compute_face_data(dof->g, e0, Vmap0[2], p0, p0 + 1, p0 + 2);
	compute_face_data(dof->g, e1, Vmap1[2], p1, p1 + 1, p1 + 2);

	/* Face 3 of e */
	p0 = new_data[7];
	p1 = new_data[8];
	compute_face_data(dof->g, e0, Vmap0[3], p0, p0 + 1, p0 + 2);
	compute_face_data(dof->g, e1, Vmap1[3], p1, p1 + 1, p1 + 2);

	/* the new face */
	p0 = new_data[9];
	compute_face_data(dof->g, e0, Vmap0[0], p0, p0 + 1, p0 + 2);

	/*------------ The volume data: volume, diameter ------------*/
	p0 = new_data[10];
	p1 = new_data[11];
	_phg_compute_elem_data(dof->g, e0, p0, p0 + 1);
	_phg_compute_elem_data(dof->g, e1, p1, p1 + 1);
	compute_jacobian(dof->g, e0, p0 + 2);
	compute_jacobian(dof->g, e1, p1 + 2);

	return;
    }
#endif

    /* Face 2 of e */
    p0 = new_data[5];
    p1 = new_data[6];
    q  = old_data[12];
    /* the areas */
    *(p0++) = *(p1++) = 0.5 * *(q++);
#if 1
    /* the diameters */
    compute_face_data(dof->g, e0, Vmap0[2], NULL, p0++, NULL);
    compute_face_data(dof->g, e1, Vmap1[2], NULL, p1++, NULL);
    /* the normals */
    p0[0] = p1[0] = *(++q); p0[1] = p1[1] = *(++q); p0[2] = p1[2] = *(++q);
    adjust_direction(dof->g, p0, e->verts[0], e->verts[1], e->verts[3],
		     e0->verts[Vmap0[0]], e0->verts[3], e0->verts[Vmap0[3]]);
    adjust_direction(dof->g, p1, e->verts[0], e->verts[1], e->verts[3],
		     e1->verts[3], e1->verts[Vmap1[1]], e1->verts[Vmap1[3]]);
#else
    compute_face_data(dof->g, e0, Vmap0[2], NULL, p0, p0 + 1);
    compute_face_data(dof->g, e1, Vmap1[2], NULL, p1, p1 + 1);
#endif

    /* Face 3 of e */
    p0 = new_data[7];
    p1 = new_data[8];
    q  = old_data[13];
    /* the areas */
    *(p0++) = *(p1++) = 0.5 * *(q++);
#if 1
    /* the diameters */
    compute_face_data(dof->g, e0, Vmap0[3], NULL, p0++, NULL);
    compute_face_data(dof->g, e1, Vmap1[3], NULL, p1++, NULL);
    /* the normals */
    p0[0] = p1[0] = *(++q); p0[1] = p1[1] = *(++q); p0[2] = p1[2] = *(++q);
    adjust_direction(dof->g, p0, e->verts[0], e->verts[1], e->verts[2],
		     e0->verts[Vmap0[0]], e0->verts[3], e0->verts[Vmap0[2]]);
    adjust_direction(dof->g, p1, e->verts[0], e->verts[1], e->verts[2],
		     e1->verts[3], e1->verts[Vmap1[1]], e1->verts[Vmap1[2]]);
#else
    compute_face_data(dof->g, e0, Vmap0[3], NULL, p0, p0 + 1);
    compute_face_data(dof->g, e1, Vmap1[3], NULL, p1, p1 + 1);
#endif

    /* the new face */
    p0 = new_data[9];
    compute_face_data(dof->g, e0, Vmap0[0], p0, p0 + 1, p0 + 2);

    /*------------ The volume data: volume, diameter ------------*/
    p0 = new_data[10];
    p1 = new_data[11];
    /* the volumes */
    *(p0++) = *(p1++) = 0.5 * *(old_data[14]);
    /* the diameters */
    _phg_compute_elem_data(dof->g, e0, NULL, p0++);
    _phg_compute_elem_data(dof->g, e1, NULL, p1++);
    /* the Jacobians. */
#if 0
    compute_jacobian(dof->g, e0, p0);
    compute_jacobian(dof->g, e1, p1);
#else
    jacobian_c2f(e, (void *)(old_data[14] + 2), (void *)p0, (void *)p1);
#endif
}

static void
interp_f2c(DOF *dof, ELEMENT *e, FLOAT **parent_data, FLOAT **children_data)
{
    ELEMENT *ee;
    int map[NVert];

    parent_data += NVert + NEdge;

#if ALLOW_CURVED_BOUNDARY
    if (e->bound_func[0] != -1) {	/* curved boundary, recompute */
	if (e->children[0] == NULL)
	    init(dof, e, FACE, 1, NULL, NULL, NULL, parent_data[1], NULL);
	if (e->children[1] == NULL)
	    init(dof, e, FACE, 0, NULL, NULL, NULL, parent_data[0], NULL);
	init(dof, e, FACE, 2, NULL, NULL, NULL, parent_data[2], NULL);
	init(dof, e, FACE, 3, NULL, NULL, NULL, parent_data[3], NULL);
	init(dof, e, VOLUME, 0, NULL, NULL, NULL, parent_data[NFace], NULL);
	return;
    }
#endif

    /* face diameter */
    if (e->children[0] == NULL)
	compute_face_data(dof->g, e, 1, parent_data[1], parent_data[1] + 1,
				parent_data[1] + 2);
    if (e->children[1] == NULL)
	compute_face_data(dof->g, e, 0, parent_data[0], parent_data[0] + 1,
				parent_data[0] + 2);
    compute_face_data(dof->g, e, 2, NULL, parent_data[2] + 1, NULL);
    compute_face_data(dof->g, e, 3, NULL, parent_data[3] + 1, NULL);
    if ((ee = e->children[0]) != NULL) {	/* interpolate using child 0 */
	phgMapP2C(e, map, NULL, 0);
	/* face areas */
	parent_data[2][0] = children_data[5][0] * 2.0;
	parent_data[3][0] = children_data[7][0] * 2.0;
	/* face normals */
	parent_data[2][2] = children_data[5][2];
	parent_data[2][3] = children_data[5][3];
	parent_data[2][4] = children_data[5][4];
	adjust_direction(dof->g, parent_data[2] + 2,
			 e->verts[0], e->verts[1], e->verts[3],
			 ee->verts[map[0]], ee->verts[3], ee->verts[map[3]]);
	parent_data[3][2] = children_data[7][2];
	parent_data[3][3] = children_data[7][3];
	parent_data[3][4] = children_data[7][4];
	adjust_direction(dof->g, parent_data[3] + 2,
			 e->verts[0], e->verts[1], e->verts[2],
			 ee->verts[map[0]], ee->verts[3], ee->verts[map[2]]);
	/* volume */
	parent_data[NFace][0] = children_data[10][0] * 2.0;
    }
    else {					/* interpolate using child 1 */
	ee = e->children[1];
	phgMapP2C(e, map, NULL, 1);
	/* face areas */
	parent_data[2][0] = children_data[6][0] * 2.0;
	parent_data[3][0] = children_data[8][0] * 2.0;
	/* face normals */
	parent_data[2][2] = children_data[6][2];
	parent_data[2][3] = children_data[6][3];
	parent_data[2][4] = children_data[6][4];
	adjust_direction(dof->g, parent_data[2] + 2,
			 e->verts[0], e->verts[1], e->verts[3],
			 ee->verts[3], ee->verts[map[1]], ee->verts[map[3]]);
	parent_data[3][2] = children_data[8][2];
	parent_data[3][3] = children_data[8][3];
	parent_data[3][4] = children_data[8][4];
	adjust_direction(dof->g, parent_data[3] + 2,
			 e->verts[0], e->verts[1], e->verts[2],
			 ee->verts[3], ee->verts[map[1]], ee->verts[map[2]]);
	/* volume */
	parent_data[NFace][0] = children_data[11][0] * 2.0;
    }
    /* diameter */
    _phg_compute_elem_data(dof->g, e, NULL, parent_data[NFace] + 1);
    /* Jacobian */
    if (e->children[0] != NULL)
	jacobian_f2c(e, (void *)(parent_data[NFace] + 2),
			(void *)(children_data[10] + 2), NULL);
    else
	jacobian_f2c(e, (void *)(parent_data[NFace] + 2),
			NULL, (void *)(children_data[11] + 2));
}

static BOOLEAN
comp_floats(int n, const FLOAT *f, const FLOAT *g, BOOLEAN signless)
{
    static FLOAT epsilon = -1.0;
    int i;
    FLOAT e = 0.;

    CheckThread
    if (epsilon < 0.)
	epsilon = Sqrt(FLOAT_EPSILON);

    for (i = 0; i < n; i++, f++, g++) {
	if ((e = Fabs((*f - *g) / (Fabs(*g) < 1.0 ? 1.0: *g))) > epsilon)
	    break;
    }

    if (i >= n)
	return TRUE;

    if (!signless) {
	phgInfo(-1, "value 0: %lg, value 1: %lg, error: %lg\n",
			(double)*f, (double)*g, (double)e);
	return FALSE;
    }

    f -= i;
    g -= i;
    for (i = 0; i < n; i++, f++, g++) {
	if ((e = Fabs((*f + *g) / (Fabs(*g) < 1.0 ? 1.0: *g))) > 1e-8) {
	    phgInfo(-1, "value 0: %lg, value 1: %lg, error: %lg\n",
			(double)*f, (double)*g, (double)e);
	    return FALSE;
	}
    }

    return TRUE;
}

void
phgGeomCheck(GRID *g)
{
    int i;
    ELEMENT *e;
    FLOAT vol, dia, d, buffer[(Dim + 1) * (Dim + 1)];

    if (g->geom == NULL)
	return;

    ForAllElements(g, e) {
	_phg_compute_elem_data(g, e, &vol, &dia);
	d = phgGeomGetVolume(g, e);
	if (!comp_floats(1, &vol, &d, FALSE))
	    phgError(1, "phgGeomCheck: element %d, volume, %lg != %lg\n",
			e->index, (double)vol, (double)d);
	d = phgGeomGetDiameter(g, e);
	if (!comp_floats(1, &dia, &d, FALSE))
	    phgError(1, "phgGeomCheck: element %d, diameter, %lg != %lg\n",
			e->index, (double)dia, (double)d);
	compute_jacobian(g, e, buffer);
	if (!comp_floats((Dim + 1) * (Dim + 1), buffer,
				phgGeomGetJacobian(g, e), FALSE))
	    phgError(1, "phgGeomCheck: element %d, Jacobian\n");
	for (i = 0; i < NFace; i++) {
	    const FLOAT *p0, *p1;
	    compute_face_data(g, e, i, &vol, &dia, buffer);
	    d = phgGeomGetFaceArea(g, e, i);
	    if (!comp_floats(1, &vol, &d, FALSE))
		phgError(1, "phgGeomCheck: element %d, face %d(%d), area, "
			    "%lg != %lg\n", e->index, i, e->faces[i],
			    (double)vol, (double)d);
	    d = phgGeomGetFaceDiameter(g, e, i);
	    if (!comp_floats(1, &dia, &d, FALSE))
		phgError(1, "phgGeomCheck: element %d, face %d(%d), diameter, "
			    "%lg != %lg\n", e->index, i, e->faces[i],
			    (double)dia, (double)d);
	    p0 = buffer;
	    p1 = phgGeomGetFaceNormal(g, e, i);
	    if (!comp_floats(Dim, p0, p1, TRUE))
		phgError(1, "phgGeomCheck: element %d, face %d(%d), normal, "
			    "(%lg,%lg,%lg) != (%lg,%lg,%lg)\n",
			    e->index, i, e->faces[i],
			    (double)p0[0], (double)p0[1], (double)p0[2],
			    (double)p1[0], (double)p1[1], (double)p1[2]);
	}
    }

    if (phgVerbosity > 0)
	phgPrintf("%s: passed.\n", __func__);
}

/*-------------------------- User interface --------------------------*/

/* Temporary geometric data pointers during mesh refinement (and coarsening?)
 *
 * Note: when interpolating DOFs after mesh refinement, the old geom data
 *	 pointers are saved, and used when geom data on non leaf elements
 *	 are requested.
 *
 *	 For mesh refinement the nonleaf elements on which old geom data
 *	 are valid have e->flag == 1. */
static FLOAT *temp_data_vert = NULL;
/*static FLOAT *temp_data_edge = NULL;*/
static FLOAT *temp_data_face = NULL;
static FLOAT *temp_data_elem = NULL;

void
phgGeomSaveNonLeafData_(GRID *g)
{
    assert(g->geom != NULL);
    temp_data_vert = g->geom->data_vert;
    /*temp_data_edge = g->geom->data_edge;*/
    temp_data_face = g->geom->data_face;
    temp_data_elem = g->geom->data_elem;
}

void
phgGeomClearNonLeafData_(void)
{
    phgFree(temp_data_vert);
    temp_data_vert =/*temp_data_edge =*/temp_data_face = temp_data_elem = NULL;
}

void
phgGeomInit_(GRID *g, BOOLEAN init_data)
/* this function must be called before distributing the grid in order to use
   geometric data */
{
    static DOF_TYPE DOF_GEOM = {DofReserved,
	"Geom", NULL, NULL, NULL, NULL, NULL, interp, interp_f2c, init,
	NULL, NULL, NULL, FE_None, FALSE, FALSE, FALSE,
	-1, 0, 0, 0, -1, 1, 0, 0, 2 + Dim, 2 + (Dim + 1) * (Dim + 1)};

    if (g == NULL)
	return;

    if (init_data) {
	if (g->geom == NULL) {
	    g->geom = phgDofNew(g, &DOF_GEOM, 1, "geom", DofInterpolation);
	}
	else if (g->geom->data == NULL) {
	    g->geom->data = phgCalloc(DofDataCount(g->geom),
						sizeof(*g->geom->data));
	    g->geom->data_vert = g->geom->data;
	    g->geom->data_edge = g->geom->data_vert +
						DofVertexDataCount(g->geom);
	    g->geom->data_face = g->geom->data_edge +
						DofEdgeDataCount(g->geom);
	    g->geom->data_elem = g->geom->data_face +
						DofFaceDataCount(g->geom);
	}
	/* Compute geometric data.
	 * Note: geom's init function doesn't use the 'f' (2nd) argument */
	phgDofSetDataByFunction(g->geom, NULL);
    }
    else if (g->geom == NULL) {
	g->geom = phgDofNew(g, &DOF_GEOM, 1, "geom", DofNoData);
    }
}

FLOAT
phgGeomGetVolume(GRID *g, ELEMENT *e)
{
    if (g->geom == NULL)
	phgError(1, "phgGeomGetVolume: geometric data unavailable.\n");

    /* Note: geom data are assumed always available on leaf elements,
     *	     during grid refinement/coarsening */
    if (!IsLeaf(e)) {
	FLOAT vol;
	if (e->flag == 1 && temp_data_elem != NULL)
	    return temp_data_elem[g->geom->count_elem * e->index];
	_phg_compute_elem_data(g, e, &vol, NULL);
	return vol;
    }
    else {
	return *(DofElementData(g->geom, e->index));
    }
}

FLOAT
phgGeomGetDiameter(GRID *g, ELEMENT *e)
{
    if (g->geom == NULL)
	phgError(1, "phgGeomGetDiameter: geometric data unavailable.\n");

    if (IsLeaf(e))
	return *(DofElementData(g->geom, e->index) + 1);

    if (e->flag == 1 && temp_data_elem != NULL)
	return temp_data_elem[g->geom->count_elem * e->index + 1];

    phgError(1, "phgGeomGetDiameter: no data for non leaf element.\n");

    return 0.;
}

#if PHG_GRID_TYPE == GRID_TYPE_HEX
/*
 * Element pointer cache and element coord cahce,
 * used by phgGeomGetJacobian and phgGeomLambda2XYZ.
 */
static SIMPLEX *e_cache = NULL;
static FLOAT *coord_cache = NULL;
#endif  /* PHG_GRID_TYPE */

const FLOAT*
phgGeomGetJacobian_(GRID *g, SIMPLEX *e, const FLOAT *lambda, FLOAT *det)
/* Compute jacobian matrix at lambda in element e, return 1/det(J) and
 * the determinate(volume). */
{
#if PHG_GRID_TYPE == GRID_TYPE_TET
    /* For linear tet element, jacobian and vol do not depend on
     * reference coordnates, hence lambda is not used, and jacobian
     * and vol are cached if possible.  */
    static FLOAT J[(Dim + 1) * (Dim + 1)];
#if USE_OMP
# pragma omp threadprivate(J)
#endif  /* USE_OMP */

    if (g->geom == NULL)
	phgError(1, "phgGeomGetJacobian: geometric data unavailable.\n");

    if (IsLeaf(e))
	return DofElementData(g->geom, e->index) + 2;

    if (e->flag == 1 && temp_data_elem != NULL)
	return temp_data_elem + g->geom->count_elem * e->index + 2;

    /* FIXME: make use jacobian of children */
    CheckThread
    compute_jacobian(g, e, J);

    if (det != NULL) {
	*det = phgGeomGetVolume(g, e);
    }
    return J;

#elif PHG_GRID_TYPE == GRID_TYPE_HEX
    /* Compute the jacobian.
     * 
     * J = dlambda / dX = (dX / dlambda)^-1, 
     * dX / dlambda = \sum_i X_i * dphi_i / dlambda.
     *
     * Cahce jacobian values only at one point lambda_p, if lambda_p
     * changes, cache is recomputed. Note: the point is distinguished
     * by its value.
     *
     * */
    CheckThread
    static FLOAT lam_cache[Dim+1];
    static FLOAT J[Dim+1][Dim+1];
    static FLOAT det0;
    //DOF_TYPE type = mapping->type;
    MAPPING *mapping;
    int i, j, k, N;
    FLOAT *v;
    const FLOAT *a;
    
    mapping = g->mapping;
    N = mapping->type->nbas;

    /* In case lambda is not present,
     * use center of element. */
    static FLOAT lambda0[Dim+1] = {0, 0, 0, 0};
    if (lambda == NULL) {
	lambda = lambda0;
	phgInfo(5, "Compute jacobian with lambda center.\n");
    }

    phgInfo(5, "EC[%x]", e_cache);
    if (e == e_cache 
	&& lambda[0] != 0
	&& lambda[1] != 0
	&& lambda[2] != 0
	&& lambda[3] != 0
	&& fabs(lambda[0] - lam_cache[0]) < 1e-12
	&& fabs(lambda[1] - lam_cache[1]) < 1e-12
	&& fabs(lambda[2] - lam_cache[2]) < 1e-12
	&& fabs(lambda[3] - lam_cache[3]) < 1e-12
	&& coord_cache != NULL) {
	if (det != NULL)
	    *det = det0;
	phgInfo(5, "#JC");
	return J[0];
    }

    /* cache coord in element */
    if (e != e_cache || coord_cache == NULL) {
	phgFree(coord_cache);
	coord_cache = phgCalloc(N * Dim, sizeof(FLOAT));
	phgDofGetElementDatas(mapping, e, coord_cache);
	/* SHOW_M(coord_cache, N, Dim); */
	FreeAtExit(coord_cache);
	e_cache = e;
	phgInfo(5, "\n#JE");
    }
    v = coord_cache;
    memcpy(lam_cache, lambda, (Dim+1)*sizeof(FLOAT));
    phgInfo(5, "#JL");
    
    assert(mapping != NULL);
    assert(mapping->dim == Dim);
    assert(mapping->type->dim == 1);
    a = mapping->type->BasGrads(mapping, e, 0, -1, lambda); /* grad basis [Nbas][Dim+1] */
    /* SHOW_M(a, N, Dim+1); */

    bzero(J, sizeof(J));
    for (i = 0; i < N; i++) {
	for (k = 0; k < Dim; k++) /* xyz */
	    for (j = 0; j < Dim; j++) /* lambda */
		J[k][j] += v[k] * a[j];
	v += Dim;
	a += Dim + 1;
    }

    /* inverse */
#  if USE_LAPACK
    {
	int OK = 0, LDA = Dim+1, M = Dim, 
	    LWORK = M*3, pvt[Dim];
	double WORK[LWORK];

	assert(sizeof(double) == sizeof(FLOAT));
	DGETRF_FOR(&M, &M, &J[0][0], &LDA, pvt, &OK); assert(OK == 0);

	det0 = 1;
	for (i = 0; i < M; i++)
	    if (pvt[i] != i) 
		det0 *= -J[i][i];
	    else
		det0 *= J[i][i];

	DGETRI_FOR(&M, &J[0][0], &LDA, pvt, WORK, &LWORK, &OK); assert(OK == 0);
    }
#  else
    /* FIXME: write PHG dense matrix operation. */
    phgError(1, "Unimplemented!\n");
#  endif

    if (det != NULL)
	*det = det0;
    return J[0];
#elif PHG_GRID_TYPE == GRID_TYPE_MIX
    phgError(1, "Unimplemented!\n");
#endif  /* PHG_GRID_TYPE */
}

#define SQUARE(x) ((x)*(x))
#define DISTANCE(x1, y1, z1, x2, y2, z2)                        \
    Sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2))
#define INNER_PRODUCT(p, q)                     \
    (*(p  ) * *(q  ) +                          \
     *(p+1) * *(q+1) +                          \
     *(p+2) * *(q+2))
#define CROSS_PRODUCT(a, b, n) {                      \
        n[0] =  (a[1] * b[2] - b[1] * a[2]);        \
        n[1] = -(a[0] * b[2] - b[0] * a[2]);        \
        n[2] =  (a[0] * b[1] - b[0] * a[1]);        \
    }

FLOAT
phgGeomGetFaceArea_(GRID *g, SIMPLEX *e, int face_no,
		    const FLOAT *lambda,
		    FLOAT *normal) 
/* Reture the area and out normal of face at reference coord lambda.
 * There are two cases: 1. linear tet element, 2. nonliner tet or
 * hex elements. TODO: choose code by these two case.  */
{
#if PHG_GRID_TYPE == GRID_TYPE_TET
    /* For linear tet element, area and normal do not depend on
     * reference coordnates, hence lambda is not used, and area and
     * normal are cached if possible.  */
    if (g->geom == NULL)
	phgError(1, "phgGeomGetFaceArea: geometric data unavailable.\n");

    if (IsLeaf(e))
	return *(DofFaceData(g->geom, e->faces[face_no]));

    if (e->flag == 1 && temp_data_face != NULL)
	return temp_data_face[g->geom->count_face * e->faces[face_no]];

    phgError(1, "phgGeomGetFaceArea: no data for non leaf element.\n");

    if (normal != NULL)
	phgGeomGetFaceOutNormal(g, e, face_no, normal);

    return 0.;
#elif PHG_GRID_TYPE == GRID_TYPE_HEX
    /* For nonliner tet or hex elements */
    MAPPING *mapping = g->mapping;
    MAPPING_TYPE *type = mapping->type;
    int i, ii, j, k, 
	N = mapping->type->nbas,
	nbas = type->nbas;
    FLOAT DXL[3][Dim], area,
	crd[N*Dim], n[Dim]; 
    const FLOAT *x, *ggi;
    const int *v, *s;
    SHORT bases[nbas];

    const FLOAT lambda0[NFace][Dim+1] = {
	{0, -1, 0, 0},
	{1, 0, 0, 0},
	{0, 1, 0, 0},
	{-1, 0, 0, 0},
	{0, 0, -1, 0},
	{0, 0, 1, 0},
    };
    if (lambda == NULL) {
	lambda = lambda0[face_no];
	phgWarning("Compute face area with lambda center.\n");
    }

    assert(mapping->type == DOF_P1 
	   || mapping->type == DOF_P2
	   || mapping->type == DOF_DG1
	   || mapping->type == DOF_DG2);
    /* GetFaceVertices(e, face, v[0], v[1], v[2], v[3]); */
    phgDofGetElementDatas(mapping, e, crd);
    bzero(DXL, sizeof(DXL));

    /*
     * Compute area & normal:
     * normal = [dx / dxi' , dy / dxi' , dz / dxi' ] \cross
     *          [dx / deta', dy / deta', dz / deta']
     *
     * How (xi', eta') maps (xi, eta, zeta) is defined in
     *    get_lambda_on_face,
     * and here we assume v_i = GetFaceVertex(face_no, i)
     *
     * dx / dxi = \sum_i x_i * d\phi_i / dxi (xi) ;
     *  */
    const int face_xi[NFace][3] = {
	{0, 2, 1},			/* +x, +z */
	{1, 2, 0},			/* +y, +z */
	{0, 2, 1},			/* -x, +z */
	{1, 2, 0},			/* -y, +z */
	{0, 1, 2},			/* +x, +y */
	{0, 1, 2},			/* +x, +y */
    };
    const int sign_xi[NFace][3] = {
	{ 1,  1,  -1},			/* +x, +z   -y*/
	{ 1,  1,  1},			/* +y, +z   +x*/
	{-1,  1,  1},			/* -x, +z   +y*/
	{-1,  1,  -1},			/* -y, +z   -x*/
	{ 1,  1,  -1},			/* +x, +y   -z*/
	{ 1,  1,  1},			/* +x, +y   +z*/
    };
    v = face_xi[face_no];
    s = sign_xi[face_no];
    x = crd;

    for (i = 0; i < nbas; i++) {
	ggi = type->BasGrads(mapping, e, i, i+1, lambda);
	for (j = 0; j < 3; j++) /* xi' */
	    for (k = 0; k < Dim; k++) /* dim */
		DXL[j][k] += x[i*Dim + k]
		    * ggi[v[j]] * s[j];
    }

    CROSS_PRODUCT(DXL[0], DXL[1], n);
    area = INNER_PRODUCT(n, n);
    area = Sqrt(area);
    assert(area > 1e-12);
    for (k = 0; k < Dim; k++)
	n[k] /= area;

    if (normal != NULL) {
	COORD *c0, *c1, *c2;
	FLOAT c[Dim];

	/* Out normal.
	 * This is a safer way than judge by element center.
	 * */
	if (INNER_PRODUCT(DXL[2], n) < 0) {
	    n[0] = -n[0];
	    n[1] = -n[1];
	    n[2] = -n[2];
	}
	memcpy(normal, n, Dim * sizeof(FLOAT));
    }

    return area;
#elif PHG_GRID_TYPE == GRID_TYPE_MIX
    phgError(1, "Unimplemented!\n");
#endif	/* PHG_GRID_TYPE */
}

FLOAT
phgGeomGetFaceAreaByIndex(GRID *g, INT face_no)
{
    if (g->geom == NULL)
	phgError(1, "phgGeomGetFaceArea: geometric data unavailable.\n");

    return *(DofFaceData(g->geom, face_no));
}

FLOAT
phgGeomGetFaceDiameter(GRID *g, ELEMENT *e, int face)
{
    if (g->geom == NULL)
	phgError(1, "phgGeomGetFaceDiameter: geometric data unavailable.\n");

    if (IsLeaf(e))
	return *(DofFaceData(g->geom, e->faces[face]) + 1);

    if (e->flag == 1 && temp_data_face != NULL)
	return temp_data_face[g->geom->count_face * e->faces[face] + 1];

    phgError(1, "phgGeomGetFaceDiameter: no data for non leaf element.\n");

    return 0.;
}

const FLOAT *
phgGeomGetFaceNormal(GRID *g, ELEMENT *e, int face)
{
    if (g->geom == NULL)
	phgError(1, "phgGeomGetFaceNormal: geometric data unavailable.\n");

    if (IsLeaf(e))
	return DofFaceData(g->geom, e->faces[face]) + 2;

    if (e->flag == 1 && temp_data_face != NULL)
	return temp_data_face + g->geom->count_face * e->faces[face] + 2;

    phgError(1, "phgGeomGetFaceNormal: no data for non leaf element.\n");

    return NULL;
}

void
phgGeomGetFaceOutNormal(GRID *g, ELEMENT *e, int face, FLOAT nv[])
{
    COORD *c0, *c1, *c2, *c;

    if (g->geom == NULL)
	phgError(1, "phgGeomGetFaceNormal: geometric data unavailable.\n");

    if (IsLeaf(e))
	memcpy(nv, DofFaceData(g->geom, e->faces[face]) + 2, Dim * sizeof(*nv));
    else if (e->flag == 1 && temp_data_face != NULL)
	memcpy(nv, temp_data_face + g->geom->count_face * e->faces[face] + 2,
			Dim * sizeof(*nv));
    else
	phgError(1, "phgGeomGetFaceOutNormal: no data for non leaf element.\n");

    c = g->verts + e->verts[face];
    c0 = g->verts + e->verts[GetFaceVertex(face, 0)];
    c1 = g->verts + e->verts[GetFaceVertex(face, 1)];
    c2 = g->verts + e->verts[GetFaceVertex(face, 2)];
    if (nv[0] * (3. * (*c)[0] - ((*c0)[0] + (*c1)[0] + (*c2)[0])) +
	nv[1] * (3. * (*c)[1] - ((*c0)[1] + (*c1)[1] + (*c2)[1])) +
	nv[2] * (3. * (*c)[2] - ((*c0)[2] + (*c1)[2] + (*c2)[2])) > 0.0) {
	nv[0] = -nv[0];
	nv[1] = -nv[1];
	nv[2] = -nv[2];
    }
}

FLOAT (*
phgGeomGetCorners(GRID *g, ELEMENT *e, FLOAT (*tetra)[Dim]))[Dim]
/* returns the coordinates of the four vertices of 'e' as FLOAT[4][3].
 * The result is stored in the user-provided buffer 'tetra' if tetra!=NULL,
 * or in a static internal thread private buffer otherwise. */
{
    static COORD tet[NVert];
#if USE_OMP
# pragma omp threadprivate(tet)
#endif
    int j;

    if (tetra == NULL)
	tetra = tet;

    for (j = 0; j < NVert; j++)
	memcpy(tetra[j], g->verts[e->verts[j]], sizeof(tetra[0]));

    return tetra;
}

FLOAT *
phgGeomXYZ2Lambda(GRID *g, ELEMENT *e, FLOAT x, FLOAT y, FLOAT z,
		  FLOAT lambda[])
{
#if PHG_GRID_TYPE == GRID_TYPE_TET
    const FLOAT *J = phgGeomGetJacobian(g, e);

    lambda[0] = J[0] * x + J[1] * y + J[2] * z + J[3];
    J += 4;
    lambda[1] = J[0] * x + J[1] * y + J[2] * z + J[3];
    J += 4;
    lambda[2] = J[0] * x + J[1] * y + J[2] * z + J[3];
    J += 4;
    lambda[3] = J[0] * x + J[1] * y + J[2] * z + J[3];

    return lambda;
#elif PHG_GRID_TYPE == GRID_TYPE_HEX
    /*
     * The mapping from x to \xi is nonliner, e.g.
     *    x(\xi) = \sum_i x_i * phi_i(\xi). 
     * Nonliner iterations are needed to compute \xi.
     *  */
    phgError(1, "Unimplemented!\n");
#elif PHG_GRID_TYPE == GRID_TYPE_MIX
    phgError(1, "Unimplemented!\n");
#endif  /* PHG_GRID_TYPE */
}

void
phgGeomLambda2XYZ(GRID *g ,SIMPLEX *e, const FLOAT lambda[],
				FLOAT *x, FLOAT *y, FLOAT *z)
{
#if PHG_GRID_TYPE == GRID_TYPE_TET
    int i;
    COORD *c;

    *x = *y = *z = 0.;
    for(i = 0; i < NVert; i++, lambda++) {
	c = g->verts + e->verts[i];
	*x += (*c)[0] * (*lambda);
	*y += (*c)[1] * (*lambda);
	*z += (*c)[2] * (*lambda);
    }

    return;

#elif PHG_GRID_TYPE == GRID_TYPE_HEX
    /* Note: use phgDofEval alternatively,
     * in that way, manage the dof data cache.
     * */
    DOF *mapping = g->mapping;
    int i, N = mapping->type->nbas;
    FLOAT *v;
    const FLOAT *a;

    /* cache coord in element */
    if (e != e_cache || coord_cache == NULL) {
	phgFree(coord_cache);
	coord_cache = phgCalloc(N * Dim, sizeof(FLOAT));
	phgDofGetElementDatas(mapping, e, coord_cache);
	FreeAtExit(coord_cache);
	e_cache = e;
	phgInfo(5, "\n#CEx");
    }
    v = coord_cache;
    a = mapping->type->BasFuncs(mapping, e, 0, -1, lambda); /* basis [Nbas] */

    *x = *y = *z = 0.;
    for(i = 0; i < N; i++) {
	*x += v[0] * (*a);
	*y += v[1] * (*a);
	*z += v[2] * (*a);
	v += Dim;
	a++;
    }

    return;
#elif PHG_GRID_TYPE == GRID_TYPE_MIX
    phgError(1, "Unimplemented!\n");
#endif  /* PHG_GRID_TYPE */
}
