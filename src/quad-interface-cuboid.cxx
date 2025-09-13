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

/* This file implements the C interface functions phgQuadInterfaceCuboid and
 * phgQuadInterfaceRectangle for quad-cuboid.
 *
 * $Id: quad-interface-cuboid.cxx,v 1.84 2022/09/15 02:16:25 zlb Exp $ */

#include <fstream>
#include <iostream>
#include <vector>

#include <phg.h>

#define Opts    _phg_quad_interface_options

static void
rect2tri(FUNC_2D ls, int ls_order, FUNC_2D ls_grad, FLOAT const rect[2][2],
	 int quad_order, FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
/* computes the integrals by splitting the rectangle into 2 triangles */
{
    int i;
    FLOAT triangle[3][2], *rule_m0, *rule_00, *rule_p0;
    FLOAT **rules[] = {rule_m, rule_0, rule_p};
    FLOAT **rules0[] = {&rule_m0, &rule_00, &rule_p0};

    for (i = 0; i < 3; i++)
	if (rules[i] == NULL)
	    rules0[i] = NULL;

    triangle[0][0] = rect[0][0];
    triangle[0][1] = rect[0][1];

    triangle[1][0] = rect[1][0];
    triangle[1][1] = rect[1][1];

    triangle[2][0] = rect[0][0];
    triangle[2][1] = rect[1][1];
    phgQuadInterfaceTriangle(ls, ls_order, ls_grad, triangle, quad_order,
				rule_m, rule_0, rule_p);

    triangle[2][0] = rect[1][0];
    triangle[2][1] = rect[0][1];
    phgQuadInterfaceTriangle(ls, ls_order, ls_grad, triangle, quad_order,
				rules0[0], rules0[1], rules0[2]);

    for (i = 0; i < 3; i++) {
	if (rules[i] == NULL)
	    continue;
	phgQuadInterfaceRuleMerge(rules[i], *rules0[i]);
	phgFree(*rules0[i]);
    }
}

#if 0
#include <math.h>
static PROJ proj = PROJ_NONE;
static void
u3(FLOAT x, FLOAT y, FLOAT z, FLOAT *f)
{
    FLOAT d;
    if (proj == PROJ_NONE) {
	*f = (x + 1.) * (y + 2.) * (z + 3.);
	return;
    }
    f[0] = f[1] = f[2] = 0.5;
    d = Sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    assert(d > 1e-15);
    d = (x + 1.) * (y + 2.) * (z + 3.) / d;
    f[0] *= d;
    f[1] *= d;
    f[2] *= d;
}
#endif

static void
cuboid2tetra(FUNC_3D ls, int ls_order, FUNC_3D ls_grad,
	     FLOAT const cuboid[2][3], int quad_order,
	     FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
/* computes the integrals by splitting the cuboid into 5 tetrahedra */
{
    int i, j, k;
    FLOAT tetra[4][3], *rule_m0, *rule_00, *rule_p0;
    FLOAT **rules[] = {rule_m, rule_0, rule_p};
    FLOAT **rules0[] = {&rule_m0, &rule_00, &rule_p0};

    /* Note: 1. The following data are from examples/cube5.dat,
     *	     2. The last tetrahedron has all interior faces. */
    int verts[8][3] =
	{{0,0,0},{1,0,0},{0,0,1},{1,0,1},{1,1,0},{1,1,1},{0,1,0},{0,1,1}};
    int elems[5][4] = {{5,6,2,7}, {5,1,2,3}, {1,6,2,0}, {1,6,5,4}, {1,6,5,2}};

    if (ls_order <= 1 && phgRank == 0) {
	static BOOLEAN warned = FALSE;
	if (!warned) {
	    warned = TRUE;
	    phgWarning("%s may not work with planar interfaces.\n", __func__);
	}
    }

    for (i = 0; i < 3; i++) {
	if (rules[i] == NULL)
	    rules0[i] = NULL;
	else
	    *rules[i] = NULL;
    }

    for (i = 0; i < 5; i++) {
	for (j = 0; j < 4; j++) {
	    for (k = 0; k < 3; k++)
		tetra[j][k] = cuboid[verts[elems[i][j]][k]][k];
	}
	phgQuadInterfaceTetrahedron(ls, ls_order, ls_grad, tetra, quad_order,
					rules0[0], rules0[1], rules0[2]);
	for (j = 0; j < 3; j++) {
	    if (rules[j] == NULL)
		continue;
#if 0
#warning
FLOAT res0;
proj = j == 1 ? PROJ_DOT : PROJ_NONE;
phgQuadInterfaceRuleApply((void *)u3, 1, proj, *rules0[j], &res0);
printf("*** %s: [%g,%g,%g],[%g,%g,%g],[%g,%g,%g],[%g,%g,%g]: %g\n", __func__, (double)tetra[0][0], (double)tetra[0][1], (double)tetra[0][2], (double)tetra[1][0], (double)tetra[1][1], (double)tetra[1][2], (double)tetra[2][0], (double)tetra[2][1], (double)tetra[2][2], (double)tetra[3][0], (double)tetra[3][1], (double)tetra[3][2], (double)res0);
#endif
	    phgQuadInterfaceRuleMerge(rules[j], *rules0[j]);
	    phgFree(*rules0[j]);
	}
    }
}

#if USE_QUAD_CUBOID

static void
output_vtk(int dim, const void *elem_ptr, const char *fn,
				int dir = -1, FLOAT a = 0., FLOAT b = 0.)
{
    FILE *fp;
    typedef FLOAT const (*elem_t)[dim];
    elem_t c = (elem_t)elem_ptr;

    if (Opts.dbg_off || !Opts.dbg_vtk)
	return;
    fp = fopen(fn, "wt");
    fprintf(fp, "# vtk DataFile Version 2.0\n"
		"vtk output\n"
		"ASCII\n"
		"DATASET UNSTRUCTURED_GRID\n"
		"POINTS %d float\n", dim == 2 && dir < 0 ? 4 : 8);

    if (dim == 2) {
	FLOAT rect[][2] = {{c[0][0], c[0][1]}, {c[1][0], c[0][1]},
			     {c[1][0], c[1][1]}, {c[0][0], c[1][1]}};
	if (dir < 0) {
	    for (int i = 0; i < 4; i++)
		fprintf(fp, "%0.16lg %0.16lg 0.0\n", (double)rect[i][0],
						     (double)rect[i][1]);
	    fprintf(fp, "CELLS 1 5\n"
			"4 0 1 2 3\n"
			"CELL_TYPES 1\n"
			"9\n");
	}
	else {
	    for (int i = 0; i < 4; i++) {
		double t[3];
		FLOAT *p = rect[i];
		for (int j = 0; j < 3; j++)
		    if (j != dir)
			t[j] = *(p++);
		t[dir] = a;
		fprintf(fp, "%0.16lg %0.16lg %0.16lg\n", t[0], t[1], t[2]);
		t[dir] = b;
		fprintf(fp, "%0.16lg %0.16lg %0.16lg\n", t[0], t[1], t[2]);
	    }
	    fprintf(fp, "CELLS 2 10\n"
			"4 0 2 4 6\n"
			"4 1 3 5 7\n"
			"CELL_TYPES 2\n"
			"9 9\n");
	}
    }
    else {
	FLOAT cuboid[][3] = {
		{c[0][0], c[0][1], c[0][2]}, {c[1][0], c[0][1], c[0][2]},
		{c[1][0], c[1][1], c[0][2]}, {c[0][0], c[1][1], c[0][2]},
		{c[0][0], c[0][1], c[1][2]}, {c[1][0], c[0][1], c[1][2]},
		{c[1][0], c[1][1], c[1][2]}, {c[0][0], c[1][1], c[1][2]}};
	for (int i = 0; i < 8; i++)
	    fprintf(fp, "%0.16lg %0.16lg %0.16lg\n", (double)cuboid[i][0],
				(double)cuboid[i][1], (double)cuboid[i][2]);
	fprintf(fp, "CELLS 1 9\n"
		    "8 0 1 2 3 4 5 6 7\n"
		    "CELL_TYPES 1\n"
		    "12\n");
    }

    fclose(fp);
    fprintf(stderr, "*** VTK file \"%s\" created.\n", fn);
}

#include "phg/quad-cuboid.hpp"

quad_cuboid::quad_cuboid(unsigned dim,
			 std::vector <FLOAT>lower,
			 std::vector <FLOAT>upper)
{
    if (dim < 1)
	throw std::invalid_argument("The dimension must be at least 1.");
    mDim = dim;
    if (lower.size() != dim || upper.size() != dim)
	throw std::invalid_argument
	    ("The dimension of the lower and upper bounds must match.");
    mLower = lower;
    mUpper = upper;
    mLevel = mLine = 0;

    //if (phgOptionsIfUsed("-qi_threshold"))
	threshold = Sqrt(1.0 - Opts.threshold * Opts.threshold);
}

void
quad_cuboid::initialize(FUNC L, FUNC gradL, int deg, unsigned order)
{
    initialize(FUNCLIST {L}, FUNCLIST {gradL}, deg, order);
}

void
quad_cuboid::initialize(const std::vector <FUNC> &L,
			const std::vector <FUNC> &gradL,
			int deg, unsigned order)
{
    std::vector <quad_struct> mQuad0;		// lower level quadrature rule

    if (mInitialized)
	return;

    mDeg = deg;
    mOrder = order;

    if (mDim > 1) {	// dimension reduction: mDim -> mDim - 1
	// mDirection := the direction the most orthogonal to the interface

	mDirection = mDim;	// invalidate mDirection

	// inadm[j] == true => direction j is inadmissible
	bool inadm[mDim];
	// grad[j] records min(cos(angle)) at the edge cuts in direction j
	std::vector<FLOAT> grad = mLower;
	for (int j = 0; j < mDim; j++) {
	    inadm[j] = false;
	    grad[j] = 1.1;
	}

	// loop on all edges, processing the cut-points on each edge
	for (int j = 0; j < mDim; j++) {
	    // loop on all edges in direction j
	    for (unsigned k = 0; k < (1 << (mDim - 1)); k++) {
		std::vector<FLOAT> x0 = mLower;
		unsigned mask = k;
		for (int i = 0; i < mDim; i++) {
		    if (i == j)
			continue;
		    if ((k & mask))
			x0[i] = mUpper[i];
		    mask >>= 1;
		}
		// loop on level-set functions
		for (int i = 0; i < L.size(); i++) {
		    auto r = solver([&L,i,j,x0, this](std::vector<FLOAT> x)
							-> std::vector<FLOAT> {
					std::vector<FLOAT> x1 = x0;
					x1[j] = x[0];
					return L[i](x1);},
				    [&gradL,i,j,x0, this](std::vector<FLOAT> x)
						-> std::vector<FLOAT> {
					std::vector<FLOAT> x1 = x0;
					x1[j] = x[0];
					return gradL[i](x1);},
				    mLower[j], mUpper[j]);
		    if (r.size() == 0)
			continue;
		    if (r.size() > 1)
			grad[j] = 0.0;
		    for (int k = 0; k < r.size(); k++) {
			x0[j] = r[k];
			auto g = gradL[i](x0);
			if (normalize(g) < Opts.EPS)
			    continue;
			for (int k = 0; k < mDim; k++) {
			    FLOAT d = Fabs(g[k]);
			    if (d <= FLOAT_EPSILON)
				inadm[k] = true;
			    if (grad[k] > d)
				grad[k] = d;
			}
		    }
		}
	    }
	}

	// First, try set mDirection := argmax_{i, grad[i]>=threshold} grad[i]
	FLOAT g = -1.0;
	for (int i = 0; i < mDim; i++) {
	    // Note: (grad[i] == 1.1) ==> no cut-points on edges in direction i
	    if (inadm[i] || grad[i] < threshold || grad[i] > 1.01
		|| grad[i] <= g)
		continue;
	    g = grad[i];
	    mDirection = i;
	}

	if (!Opts.dbg_off && Opts.show_directions) {
	    std::vector <FLOAT> len = mUpper;
	    for (int i = 0; i < len.size(); i++)
		len[i] -= mLower[i];
	    std::cerr << "*** Setting mDirection: mDim=" << mDim
		      << ", mLevel=" << mLevel << ", Elem=" << mLower
		      << "+" << len << "\n";
	    if (mDirection < mDim)
		std::cerr << "***    done in step 1: mDirection=" << mDirection
			  << ", g=" << g << "\n";
	}

	// Fallback to the old code if the previous step failed:
	//	use the direction with largest component of the gradient at the
	//	center, if there are multiple levelset functions, choose the
	//	largest of the smallest component of the normalized gradients
	//	of all the levelset functions
	if (mDirection >= mDim) {
	    mDirection = mDim;
	    // compute grad[j] := min_{i=0,...,gradL.size()} grad_L[i](midp)[j],
	    // where midp := the center, j = 0, ..., mDim
	    std::vector <FLOAT> midp;
	    for (int j = 0; j < mDim; j++)
		midp.emplace_back(0.5 * (mLower[j] + mUpper[j]));
	    std::vector <FLOAT> grad;
	    for (int i = 0; i < gradL.size(); i++) {
		auto g = gradL[i](midp);
		normalize(g);
		for (int j = 0; j < mDim; j++)
		    if (i == 0)
			grad.emplace_back(Fabs(g[j]));
		    else if (grad[j] > Fabs(g[j]))
			grad[j] = Fabs(g[j]);
	    }
	    assert(grad.size() == mDim);
	    FLOAT d = -1.0;
	    for (int j = 0; j < mDim; j++) {
		if (!inadm[j] && Fabs(grad[j]) > d) {
		    mDirection = j;
		    d = Fabs(grad[j]);
		}
	    }

	    if (mDirection >= mDim) {
		// cannot find a suitable integration direction
		SPLIT(L, gradL);
		return;
	    }

	    if (!Opts.dbg_off && Opts.show_directions)
		std::cerr << "***    done in step 2: mDirection=" << mDirection
			  << ", g=" << Fabs(grad[mDirection]) << "\n";
	}

	mInfo[mLevel].dir = mDirection;
	mInfo[mLevel].lower = mLower[mDirection];
	mInfo[mLevel].upper = mUpper[mDirection];

	std::vector <FLOAT> lower = mLower, upper = mUpper;
	lower.erase(lower.begin() + mDirection);
	upper.erase(upper.begin() + mDirection);
	std::vector <FUNC> projL, projgradL;
	projL.reserve(2 * L.size());
	projgradL.reserve(2 * gradL.size());
	for (int i = 0; i < L.size(); i++) {
	    projL.emplace_back([&L, i, this](std::vector <FLOAT>x)
						->std::vector <FLOAT> {
		x.insert(x.begin() + mDirection, mLower[mDirection]);
		return L[i](x);
	    });
	    projL.emplace_back([&L, i, this](std::vector <FLOAT>x)
						->std::vector <FLOAT> {
		x.insert(x.begin() + mDirection, mUpper[mDirection]);
		return L[i](x);
	    });
	    projgradL.emplace_back([&gradL, i, this](std::vector <FLOAT>x)
						->std::vector <FLOAT> {
		x.insert(x.begin() + mDirection, mLower[mDirection]);
		auto y = gradL[i](x);
		y.erase(y.begin() + mDirection);
		return y;
	    });
	    projgradL.emplace_back([&gradL, i, this](std::vector <FLOAT>x)
						->std::vector <FLOAT> {
		x.insert(x.begin() + mDirection, mUpper[mDirection]);
		auto y = gradL[i](x);
		y.erase(y.begin() + mDirection);
		return y;
	    });
	}
	quad_cuboid quadObj(mDim - 1, lower, upper);
	quadObj.mLine = __LINE__;
	quadObj.mLevel = mLevel + 1;
	memcpy(quadObj.mInfo, mInfo, sizeof(mInfo));
	quadObj.initialize(projL, projgradL, deg, order);
	mQuad0 = quadObj.mQuad;
    }

    // Construct tensor product rule mQuad0 x (Gauss-Legendre in mDirection)

    int size = (L.size() + 1) * (mDim == 1 ? 1 : mQuad0.size());
    if (mLevel == 0) {
	mPosQuad.reserve(((size + 1) / 2) * ((mOrder + 1) / 2));
	mNegQuad.reserve(((size + 1) / 2) * ((mOrder + 1) / 2));
	mSurfQuad.reserve(size);
	mSurfNormVecs.reserve(size);
    }
    else {
	mQuad.reserve(size * ((mOrder + 1) / 2));
    }

    // loop on quadrature points from the lower level (mQuad0)
    for (int j = 0; j < (mDim > 1 ? mQuad0.size() : 1); j++) {
	FLOAT wgt = mDim > 1 ? mQuad0[j].mW : 1.0;
	std::vector <FLOAT> x0;
	std::set <FLOAT> roots;
	if (mDim > 1)
	    x0 = mQuad0[j].mAbs;
	x0.insert(x0.begin() + mDirection, 0.0);

	roots.insert(mLower[mDirection]);
	roots.insert(mUpper[mDirection]);

	// add cut-points with the interfaces
	for (unsigned i = 0; i < L.size(); i++) {
	    //**************** intersection of Gamma with the integration line
	    auto r = solver([&L, i, x0, this](std::vector<FLOAT> x)
						-> std::vector<FLOAT> {
				std::vector <FLOAT> x1 = x0;
				x1[mDirection] = x[0];
				return L[i](x1);},
			    [&gradL, i, x0, this](std::vector<FLOAT> x)
						-> std::vector<FLOAT> {
				std::vector<FLOAT> x1 = x0;
				x1[mDirection] = x[0];
				return gradL[i](x1);},
			    mLower[mDirection], mUpper[mDirection]);
	    if (r.size() > 0) {
		if (mDim > 1 && r.size() > 1) {
		    SPLIT(L, gradL);
		    return;
		}
		// check cos(angle) w.r.t. mDirection and split if small!!!
		x0[mDirection] = r[0];
		auto g = gradL[i](x0);
		FLOAT d = normalize(g);
		// Split if theta := angle(normal,mDirection) is close to +-\pi
		// Note: grad[mDirection] = cos(theta), while Opts.threshold in
		// 	 quad-interface.c is for sin(theta), so:
		// 		threshold := Sqrt(1. - Pow(Opts.threshold, 2));
		if (d >= Opts.EPS && mDim > 1
			&& Fabs(g[mDirection]) < threshold) {
		    SPLIT(L, gradL);
		    return;
		}
		if (mLevel == 0) {
		    // save surface integral quadrature point and normal vector
		    if (Fabs(g[mDirection]) <= Opts.EPS) {
			SPLIT(L, gradL);
			return;
		    }
	            FLOAT t = 1.0 / Fabs(g[mDirection]);
		    mSurfQuad.emplace_back(wgt * t, x0);
		    mSurfNormVecs.emplace_back(g);
		}
		FLOAT last = mLower[mDirection];
		FLOAT eps = Fabs(mUpper[mDirection] - last) * Opts.EPS;
		for (int i = 0; i < r.size(); i++) {
		    if (Fabs(last - r[i]) < eps)
			continue;
		    roots.insert(last = r[i]);
		}
	    }
	}

	std::vector <FLOAT>intervals(roots.begin(), roots.end());
	// apply Gauss-Legendre quadrature to each interval
	for (unsigned i = 0; i < intervals.size() - 1; i++) {
	    QUAD *quad = phgQuadGetQuad1D(mOrder);
	    // a := left-end of the interval, d := length of the interval
	    FLOAT a = intervals[i], d = intervals[i + 1] - a, w = wgt * d;
	    if (Fabs(d) < Opts.EPS)
		continue;
	    std::vector <quad_struct> *q;
	    if (mLevel > 0) {
		q = &mQuad;	// ignore sign at lower levels
	    }
	    else {
		assert(L.size() == 1);
		// get sign at the middle point
		x0[mDirection] = 0.5 * d + a;
		FLOAT sign = L[0](x0)[0];
		//assert(Fabs(sign) > Opts.EPS);
		q = sign < 0.0 ? &mNegQuad : &mPosQuad;
	    }
	    for (int k = 0; w != 0 && k < quad->npoints; k++) {
		x0[mDirection] = quad->points[2 * k] * d + a;
		q->emplace_back(w * quad->weights[k], x0);
	    }
	}
    }

#if 0
#warning This is for debugging reserved sizes!
std::cerr << "*** Dim=" << mDim << ", actual/reserved size: "
	  << "mQuad=" << mQuad.size() << "/" << mQuad.capacity() << ", "
	  << "Pos=" << mPosQuad.size() << "/" << mPosQuad.capacity() << ", "
	  << "Neg=" << mNegQuad.size() << "/" << mNegQuad.capacity() << ", "
	  << "Surf=" << mSurfQuad.size() << "/" << mSurfQuad.capacity() << "\n";
#endif

    if (mLevel == 0) {
	mPosQuad.shrink_to_fit();
	mNegQuad.shrink_to_fit();
	mSurfQuad.shrink_to_fit();
	mSurfNormVecs.shrink_to_fit();
	mInitialized = 1;
    }

    return;
}

FLOAT
quad_cuboid::normalize(std::vector<FLOAT> &g)
// normalizes a vector, returns its original norm
{
    FLOAT d = g[0];
    if (mDim == 1) {
	g[0] = 1.0;
	return Fabs(d);
    }
    d *= d;
    for (int j = 1; j < g.size(); j++)
	d += g[j] * g[j];
    if (d == 0.0)
	return d;
    FLOAT f = 1.0 / (d = Sqrt(d));
    for (int j = 0; j < g.size(); j++)
	g[j] *= f;
    return d;
}

void quad_cuboid::split(const std::vector<FUNC> &L,
                        const std::vector<FUNC> &gradL,
			const char *file, int line)
// recursion: subdivide an element into 4 (2D) or 8 (3D) sub-elements.
{
    static int level[4] = {0, 0, 0, 0};
#if USE_OMP
# pragma omp threadprivate(level)
#endif	/* USE_OMP */

    if (!Opts.dbg_off && Opts.show_recursions) {
	const char *fn = strrchr(file, '/');
	fn = fn == NULL ? file : fn + 1;
	phgInfo(-1, "split triggered (line %d): mDim=%d, mDirection=%d, "
		    "level=%d, mLine=%d\n", line, mDim, mDirection,
		    level[mDim], mLine);
	if (Opts.dbg_vtk) {
	    // This saves the VTK files from top down, and always the newest
	    // sequence is saved.
	    char fn[128], prefix[16];
	    if (phgNProcs > 1)
		sprintf(prefix, "tmp%04d_", phgRank);
	    else
		sprintf(prefix, "tmp");
	    if (level[mDim] == 0) {
		// remove previously saved files
		sprintf(fn, "/bin/rm -f %s%d-*.vtk", prefix, mDim);
		system(fn);
	    }
	    sprintf(fn, "%s%d-%d.vtk", prefix, mDim, level[mDim]);
	    if (mDim == 2) {
		FLOAT e[][2] = {{mLower[0], mLower[1]}, {mUpper[0], mUpper[1]}};
		if (mLevel == 1)
		    output_vtk(2, e, fn,
				mInfo[0].dir, mInfo[0].lower, mInfo[0].upper);
		else
		    output_vtk(2, e, fn);
	    }
	    else if (mDim == 3) {
		FLOAT e[][3] = {{mLower[0], mLower[1], mLower[2]},
				{mUpper[0], mUpper[1], mUpper[2]}};
		output_vtk(3, e, fn);			
	    }
	}
    }

    mQuad.clear();
    mPosQuad.clear();
    mNegQuad.clear();
    mSurfQuad.clear();
    mSurfNormVecs.clear();

    std::vector<FLOAT> h = mLower;
    FLOAT measure = 1.0;
    for (unsigned i = 0; i < mDim; i++) {
        h[i] = 0.5 * (mUpper[i] - mLower[i]);
	FLOAT tmp = h[i] + mLower[i];
	if (tmp == mLower[i] || tmp == mUpper[i])
	    return;	// edge length <= FLOAT precision
	measure *= h[i];
    }
    if (Fabs(measure) <= FLOAT_EPSILON)
	return;

    level[mDim]++;

    if (level[mDim] > Opts.subdiv_limit)
	phgError(1, "%s:%d: mDim = %d, subdivision limit exceeded, abort.",
		 __FILE__, __LINE__, mDim);

    // reserve spaces for quadrature points and normal vectors
#if 0	// unnecessary (reallocated 4 or 8 times only, and hard to estimate)
    if (mLevel > 0) {
	mQuad.reserve(6 * (size_t)(pow(mOrder + 1, mDim) + 0.5));
    }
    else {
	mPosQuad.reserve(3 * (size_t)(pow(mOrder + 1, mDim) + 0.5));
	mNegQuad.reserve(3 * (size_t)(pow(mOrder + 1, mDim) + 0.5));
	mSurfQuad.reserve(3 * (size_t)(pow(mOrder + 1, mDim - 1) + 0.5));
	mSurfNormVecs.reserve(3 * (size_t)(pow(mOrder + 1, mDim - 1) + 0.5));
    }
#endif

    for (int i = 0; i < (1 << mDim); i++) {
        std::vector<FLOAT> lower = mLower;
        for (unsigned j = 0; j < mDim; j++) {
            lower[j] += int((i & (1 << j)) != 0) * h[j];
        }
        std::vector<FLOAT> upper = lower;
        for (unsigned j = 0; j < mDim; j++) {
            upper[j] += h[j];
        }
        quad_cuboid quadObj(mDim, lower, upper);
	quadObj.mLine = __LINE__;
	quadObj.mLevel = mLevel;
	memcpy(quadObj.mInfo, mInfo, sizeof(mInfo));
        quadObj.initialize(L, gradL, mDeg, mOrder);
	if (mLevel > 0) {
            mQuad.insert(mQuad.end(), quadObj.mQuad.begin(),
					quadObj.mQuad.end());
	}
	else {
            mPosQuad.insert(mPosQuad.end(), quadObj.mPosQuad.begin(),
					quadObj.mPosQuad.end());
            mNegQuad.insert(mNegQuad.end(), quadObj.mNegQuad.begin(),
					quadObj.mNegQuad.end());
            mSurfQuad.insert(mSurfQuad.end(), quadObj.mSurfQuad.begin(),
					quadObj.mSurfQuad.end());
            mSurfNormVecs.insert(mSurfNormVecs.end(),
					quadObj.mSurfNormVecs.begin(),
					quadObj.mSurfNormVecs.end());
	}
    }

    level[mDim]--;

#if 0
#warning
    if (!Opts.dbg_off && Opts.show_recursions && Opts.dbg_vtk) {
	static int level_max[4] = {-1, -1, -1, -1};
	static bool flag[4] = {false, false, false, false};
	// This saves the VTK files from bottom up, and only when the current
	// sequence is deeper than all previous ones.
#if USE_OMP
# pragma omp threadprivate(level_max, flag)
#endif	/* USE_OMP */
	if (level[mDim] > level_max[mDim]) {
	    level_max[mDim] = level[mDim];
	    flag[mDim] = true;
	}
	if (flag[mDim]) {
	    // output the current sequence of sub-elements
	    char fn[128];
	    if (phgNProcs > 1)
		sprintf(fn, "tmp%04d_%d-%d.vtk", phgRank, mDim, level[mDim]);
	    else
		sprintf(fn, "tmp%d-%d.vtk", mDim, level[mDim]);
	    if (mDim == 2) {
		FLOAT e[][2] = {{mLower[0], mLower[1]}, {mUpper[0], mUpper[1]}};
		if (mLevel == 1)
		    output_vtk(2, e, fn,
				mInfo[0].dir, mInfo[0].lower, mInfo[0].upper);
		else
		    output_vtk(2, e, fn);
	    }
	    else if (mDim == 3) {
		FLOAT e[][3] = {{mLower[0], mLower[1], mLower[2]},
				{mUpper[0], mUpper[1], mUpper[2]}};
		output_vtk(3, e, fn);			
	    }
	    if (level[mDim] == 0)
		flag[mDim] = false;
	}
    }
#endif
}

bool
quad_cuboid::refinement(const FUNC &F, const FUNC &gradF, FLOAT x1, FLOAT x2,
			FLOAT *root)
{
    const int MAXIT = 100;
    FLOAT xh, xl;
    FLOAT fl = F({x1})[0];
    FLOAT fh = F({x2})[0];
    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
	return false;
    }
    if (fl == 0.0) {
	*root = x1;
	return true;
    }
    if (fh == 0.0) {
	*root = x2;
	return true;
    }
    if (fl < 0.0) {
	xl = x1;
	xh = x2;
    }
    else {
	xh = x1;
	xl = x2;
    }
    FLOAT rts = 0.5 * (x1 + x2);
    FLOAT dxold = Fabs(x2 - x1);
    FLOAT dx = dxold;
    FLOAT f = F({rts})[0];
    FLOAT df = gradF({rts})[0];
    for (int j = 0; j < MAXIT; j++) {
	if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) ||
	    (Fabs(2.0 * f) > Fabs(dxold * df))) {
	    dxold = dx;
	    dx = 0.5 * (xh - xl);
	    rts = xl + dx;
	    if (xl == rts) {
		*root = rts;
		return true;
	    }
	}
	else {
	    dxold = dx;
	    dx = f / df;
	    FLOAT temp = rts;
	    rts -= dx;
	    if (temp == rts) {
		*root = rts;
		return true;
	    }
	}
	f = F({rts})[0];
	if (Fabs(dx) < FLOAT_EPSILON) {
	    *root = rts;
	    return true;
	}
	df = gradF({rts})[0];
	if (f < 0.0)
	    xl = rts;
	else
	    xh = rts;
    }
    return false;
}

std::vector<FLOAT>
quad_cuboid::solver(const FUNC &F, const FUNC &gradF, FLOAT a, FLOAT b)
/* solves F==0 in [a,b], returns number of roots found, and set *root to one
 * of the roots found. */
{
    int deg = abs(mDeg), *mults;
    std::vector<FLOAT> roots;

    if (deg == 0) {
	if (Fabs(F({a})[0]) < Opts.EPS)
	    roots.emplace_back(a);
	if (Fabs(F({b})[0]) < Opts.EPS)
	    roots.emplace_back(b);
	return roots;
    }

    FLOAT *p = new FLOAT[deg + 1], h = (b - a) / deg;
    for (int i = 0; i <= deg; i++)
	p[i] = F({a + i * h})[0];
    int n = phgUnivariatePolyRoots(deg, a, b, 0.0, p, &mults);

    if (n > 0) {
	if (mDeg >= 0) {
	    for (int i = 0; i < n; i++)
		for (int j = 0; j < mults[i]; j++)
		    roots.emplace_back(p[i]);
	    std::sort(roots.begin(), roots.end());
	}
	else {
	    FLOAT r;
	    std::sort(p, p + n);
	    std::vector <FLOAT>intervals(n + 1);
	    intervals[0] = a;
	    intervals[n] = b;
	    for (int i = 1; i < n; i++)
		intervals[i] = 0.5 * (p[i - 1] + p[i]);
	    for (int i = 0; i < n; i++)
		if (refinement(F, gradF, intervals[i], intervals[i + 1], &r))
		    roots.emplace_back(r);
	}
    }

    phgFree(mults);
    delete[]p;

    return roots;
}

/******************************************************************************
 * Interface functions
 ******************************************************************************/

using namespace std;

namespace {	// local functions and variables
///////////////////////////////////////////////////////////////////////////////

void *L_func, *L_grad;
int dim, type;	// space dimension (2 or 3) and quadrature type (-1, 0 or 1)
void *elem;	// pointer to the current element (FLOAT[2][dim])
#if USE_OMP
# pragma omp threadprivate(L_func, L_grad, dim, type, elem)
#endif	/* USE_OMP */


// integrand (1) for testing (-qi_dbg_elem)

void
one2(FLOAT x, FLOAT y, FLOAT *res)
{
    if (type != 0) {
	*res = 1.0;
    }
    else {
	((FUNC_2D)L_grad)(x, y, res);
	FLOAT d = Sqrt(res[0] * res[0] + res[1] * res[1]);
	assert(d > 1e-15);
	d = 1.0 / d;
	res[0] *= d;
	res[1] *= d;
    }
}

void
one3(FLOAT x, FLOAT y, FLOAT z, FLOAT *res)
{
    if (type != 0) {
	*res = 1.0;
    }
    else {
	((FUNC_3D)L_grad)(x, y, z, res);
	FLOAT d = Sqrt(res[0] * res[0] + res[1] * res[1] + res[2] * res[2]);
	assert(d > 1e-15);
	d = 1.0 / d;
	res[0] *= d;
	res[1] *= d;
	res[2] *= d;
    }
}

vector<FLOAT>
levelsetfunc(vector<FLOAT> x) {
    FLOAT val;
    if (dim == 2)
	((FUNC_2D)L_func)(x[0], x[1], &val);
    else
	((FUNC_3D)L_func)(x[0], x[1], x[2], &val);
    return {val};
}

vector<FLOAT>
gradlevelsetfunc(vector<FLOAT> x) {
    FLOAT grad[3];
    if (dim == 2) {
	((FUNC_2D)L_grad)(x[0], x[1], grad);
	return {grad[0], grad[1]};
    }
    else {
	((FUNC_3D)L_grad)(x[0], x[1], x[2], grad);
	return {grad[0], grad[1], grad[2]};
    }
}

void
get_rules(void *ls, int ls_order, void *ls_grad,
	  vector<FLOAT> lower, vector<FLOAT> upper,
	  int order, FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
{
    int n;
    FLOAT **rule, *ptr, d;

    if (rule_m != NULL)
	*rule_m = NULL;
    if (rule_0 != NULL)
	*rule_0 = NULL;
    if (rule_p != NULL)
	*rule_p = NULL;

    if (ls == NULL) {
	// integrate on the whole element and store it in rule_m
	rule = rule_m != NULL ? rule_m : rule_p;
	if (rule == NULL) {
	    /* return an empty rule in rule_0 */
	    rule = rule_0;
	    *rule = (FLOAT *)phgAlloc(sizeof(FLOAT) * RULE_HEADER);
	    SetRuleHeader(*rule, 0, dim, FALSE);
	    return;
	}
	QUAD *quad = phgQuadGetQuad1D(order);
	int m = quad->npoints;	    // # points per dimension
	n = m;
	for (int i = 1; i < dim; i++)
	    n *= m;
	int size = dim + 1;
	*rule = (FLOAT *)phgAlloc(sizeof(FLOAT) * (RULE_HEADER + n * size));
	ptr = *rule + RULE_HEADER;
	SetRuleHeader(*rule, n, dim, FALSE);
	if (dim == 2) {
	    FLOAT W = upper[0] - lower[0],
		  H = upper[1] - lower[1];
	    d = Fabs(W * H);
	    for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++, ptr += size) {
		    ptr[0] = lower[0] + quad->points[2 * i] * W;
		    ptr[1] = lower[1] + quad->points[2 * j] * H;
		    ptr[2] = quad->weights[i] * quad->weights[j] * d;
		}
	    }
	}
	else {
	    FLOAT W = upper[0] - lower[0],
		  H = upper[1] - lower[1],
		  D = upper[2] - lower[2];
	    d = Fabs(W * H * D);
	    for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
		    for (int k = 0; k < m; k++, ptr += size) {
			ptr[0] = lower[0] + quad->points[2 * i] * W;
			ptr[1] = lower[1] + quad->points[2 * j] * H;
			ptr[2] = lower[2] + quad->points[2 * k] * D;
			ptr[3] = quad->weights[i] * quad->weights[j] *
				 quad->weights[k] * d;
		    }
		}
	    }
	}
    }
    else {	/* ls == NULL */
	/* integration in a cut element */
	quad_cuboid quadObj(dim, lower, upper);
	L_func = ls;
	L_grad = ls_grad;
	quadObj.initialize(levelsetfunc, gradlevelsetfunc, ls_order, order);
	for (type = -1; type <= 1; type++) {
	    vector<quad_struct> *qRule;
	    if (type < 0) {
		qRule = &quadObj.mNegQuad;
		rule = rule_m;
	    }
	    else if (type > 0) {
		qRule = &quadObj.mPosQuad;
		rule = rule_p;
	    }
	    else {
		qRule = &quadObj.mSurfQuad;
		rule = rule_0;
	    }

	    if (rule == NULL)
		continue;

	    //cout << "Number of quadrature points: " << qRule->size() << '\n';

	    n = qRule->size();
	    BOOLEAN nv_flag = (type == 0 && Opts.nv_flag);
	    int size = dim + 1 + (nv_flag ? dim : 0);
	    *rule = (FLOAT *)phgAlloc(sizeof(FLOAT) * (RULE_HEADER + n * size));
	    ptr = *rule + RULE_HEADER;
	    FLOAT *nv = ptr + n * (dim + 1);
	    SetRuleHeader(*rule, n, dim, nv_flag);
	    // copy quadrature rule
	    for (auto &x : *qRule) {
		for (int i = 0; i < x.mAbs.size(); i++)
		    *(ptr++) = x.mAbs[i];
		*(ptr++) = x.mW;
	    }
	    // copy normal vectors of the interface
	    if (nv_flag)
		for (auto &x : quadObj.mSurfNormVecs)
		    for (int i = 0; i < x.size(); i++)
			*(nv++) = x[i];
	}    /* for (type ...) */
    }	/* ls == NULL */

    if (!Opts.dbg_off && Opts.dbg_elem) {
	/* compare with results obtained by cuboid2tetra/rect2triangle */
	static FLOAT max_error[] = {0., 0., 0.};
#if USE_OMP
# pragma omp threadprivate(max_error)
#endif  /* USE_OMP */
	FLOAT error, res, res0;
	FLOAT *rule0_m, *rule0_0, *rule0_p;
	FLOAT **rules [] = {rule_m, rule_0, rule_p};
	FLOAT **rules0[] = {&rule0_m, &rule0_0, &rule0_p};
	void *func = dim == 2 ? (void *)one2 : (void *)one3;
	for (int i = 0; i < 3; i++)
	    if (rules[i] == NULL)
		rules0[i] = NULL;
	/* compute reference values with rect2tri or cuboid2tetra */
	Opts.dbg_off = TRUE;
	if (dim == 2)
	    rect2tri((FUNC_2D)ls, ls_order, (FUNC_2D)ls_grad,
		     (const FLOAT (*)[2])elem, order,
		     rules0[0], rules0[1], rules0[2]);
	else
	    cuboid2tetra((FUNC_3D)ls, ls_order, (FUNC_3D)ls_grad,
			 (const FLOAT (*)[3])elem, order,
			 rules0[0], rules0[1], rules0[2]);
	Opts.dbg_off = FALSE;
	for (type = -1; type <= 1; type++) {
	    PROJ proj = type == 0 ? PROJ_DOT : PROJ_NONE;
	    if (rules[1 + type] == NULL)
		continue;
	    phgQuadInterfaceRuleApply(func, 1, proj, *rules[type + 1], &res);
	    phgQuadInterfaceRuleApply(func, 1, proj, *rules0[type + 1], &res0);
#if 0
#warning
	    /* This has same effect as "-qi_cuboid2tetra -qi_rect2triangle" */
	    phgFree(*rules[type + 1]);
	    *rules[type + 1] = *rules0[type + 1];
	    *rules0[type + 1] = NULL;
#endif
	    phgFree(*rules0[type + 1]);
	    error = Fabs(res0) > Opts.EPS ? Fabs(1.0 - res / res0) : 0.0;
	    if (error > max_error[1 + type]) {
		max_error[1 + type] = error;
		if (dim == 2) {
		    const FLOAT (*re)[2] = (const FLOAT (*)[2])elem;
		    phgInfo(-1, "E([%0.16lg,%0.16lg]-[%0.16lg,%0.16lg]), "
				"order %d, type %d, "
				"res: %lg, ref: %lg, error: %lg\n",
			(double)re[0][0], (double)re[0][1],
			(double)re[1][0], (double)re[1][1],
			order, type, (double)res, (double)res0, (double)error);
		}
		else {
		    const FLOAT (*cu)[3] = (const FLOAT (*)[3])elem;
		    phgInfo(-1, "E([%0.16lg,%0.16lg,%0.16lg]-"
				"[%0.16lg,%0.16lg,%0.16lg]), order %d, "
				"type %d, res: %lg, ref: %lg, error: %lg\n",
			(double)cu[0][0], (double)cu[0][1], (double)cu[0][2],
			(double)cu[1][0], (double)cu[1][1], (double)cu[1][2],
			order, type, (double)res, (double)res0, (double)error);
		}
		char fn[128];
		if (phgNProcs > 1)
		    sprintf(fn, "tmp%04d_elem-%d.vtk", phgRank, dim);
		else
		    sprintf(fn, "tmp_elem-%d.vtk", dim);
		output_vtk(dim, elem, fn);
	    }
	}
    }
}

///////////////////////////////////////////////////////////////////////////////
}	/* namespace */

#endif	/* USE_QUAD_CUBOID */

extern "C" void
phgQuadInterfaceRectangle(FUNC_2D ls, int ls_order, FUNC_2D ls_grad,
	FLOAT const rect[2][2], int order,
	FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
{
    if (rule_m == NULL && rule_0 == NULL && rule_p == NULL)
	return;
#if USE_QUAD_CUBOID
    if (Opts.rect2tri == FALSE) {
	dim = 2;
	elem = (void *)rect;
	get_rules((void *)ls, ls_order, (void *)ls_grad,
		  {rect[0][0], rect[0][1]}, {rect[1][0], rect[1][1]}, 
		  order, rule_m, rule_0, rule_p);
    } else
#endif	/* USE_QUAD_CUBOID */
    rect2tri(ls, ls_order, ls_grad, rect, order, rule_m, rule_0, rule_p);
}

extern "C" void
phgQuadInterfaceCuboid(FUNC_3D ls, int ls_order, FUNC_3D ls_grad,
	FLOAT const cuboid[2][3], int order,
	FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
{
    if (rule_m == NULL && rule_0 == NULL && rule_p == NULL)
	return;

#if USE_QUAD_CUBOID
    if (Opts.cuboid2tetra == FALSE) {
	dim = 3;
	elem = (void *)cuboid;
	get_rules((void *)ls, ls_order, (void *)ls_grad,
		  {cuboid[0][0], cuboid[0][1], cuboid[0][2]}, 
		  {cuboid[1][0], cuboid[1][1], cuboid[1][2]}, 
		  order, rule_m, rule_0, rule_p);
    } else
#endif	/* USE_QUAD_CUBOID */
    cuboid2tetra(ls, ls_order, ls_grad, cuboid, order, rule_m, rule_0, rule_p);
}

/*---------------------------------------------------------------------------*/

static struct QF_CTX_ {
    FUNC_3D	ls, ls_grad;
    FLOAT const	(*cuboid)[Dim];
    int		ls_order, face;
} *qf_ctx = NULL;
#if USE_OMP
# pragma omp threadprivate(qf_ctx)
#endif	/* USE_OMP */

static void
rect2xyz(FLOAT x1, FLOAT x2, FLOAT *x, FLOAT *y, FLOAT *z)
{
    switch (qf_ctx->face / 2) {
	case 0:
	    *x = qf_ctx->cuboid[qf_ctx->face&1][0];
	    *y = x1;
	    *z = x2;
	    break;
	case 1:
	    *x = x1;
	    *y = qf_ctx->cuboid[qf_ctx->face&1][1];
	    *z = x2;
	    break;
	case 2:
	    *x = x1;
	    *y = x2;
	    *z = qf_ctx->cuboid[qf_ctx->face&1][2];
	    break;
    }
}

static void
ls2(FLOAT x1, FLOAT x2, FLOAT *val)
{
    FLOAT x, y, z;

    rect2xyz(x1, x2, &x, &y, &z);
    qf_ctx->ls(x, y, z, val);
}

static void
ls2_grad(FLOAT x1, FLOAT x2, FLOAT *grad)
{
    FLOAT grad3[Dim], x, y, z;

    rect2xyz(x1, x2, &x, &y, &z);
    qf_ctx->ls_grad(x, y, z, grad3);

    switch (qf_ctx->face / 2) {
	case 0: grad[0] = grad3[1]; grad[1] = grad3[2]; break;
	case 1: grad[0] = grad3[0]; grad[1] = grad3[2]; break;
	case 2: grad[0] = grad3[0]; grad[1] = grad3[1]; break;
    }
}

extern "C" void
phgQuadInterfaceCuboidFace(FUNC_3D ls, int ls_order, FUNC_3D ls_grad,
		FLOAT const cuboid[2][Dim], int face,
		int quad_order, FLOAT **rule_m, FLOAT **rule_0, FLOAT **rule_p)
/* returns qudrature rules *rule_m, rule_0 and *rule_p, for numerical
 * integration in Face\cap\Omega^-, Face\cap\Gamma and Face\cap\Omega^+,
 * respectively.
 *
 * Note: in the resulting rules, the weights are scaled to match the area of
 * the face. 
 *
 * If ls == NULL, then returns a quadrature rule for the whole face in
 * *rule_m, or *rule_p if rule_m == NULL.
 *
 * Face numbering:
 * 	0: x==xmin, 1: x==xmax, 2: y==ymin, 3: y==ymax, 4: z==zmin, 5: z==zmax.
 */
{
    FLOAT rect[2][2];
    int k, i, np, i0, i1, i2;

    qf_ctx = (struct QF_CTX_ *)phgCalloc(1, sizeof(*qf_ctx));
    qf_ctx->ls = ls;
    qf_ctx->ls_grad = ls_grad;
    qf_ctx->ls_order = ls_order;
    qf_ctx->cuboid = cuboid;
    qf_ctx->face = face;

    i0 = face / 2;
    switch (i0) {
	case 0:	i1 = 1; i2 = 2; break;
	case 1: i1 = 0; i2 = 2; break;
	case 2: i1 = 0; i2 = 1; break;
	default: phgError(1, "%s: invalid face no.: %d\n", face);
    }

    rect[0][0] = cuboid[0][i1];
    rect[0][1] = cuboid[0][i2];
    rect[1][0] = cuboid[1][i1];
    rect[1][1] = cuboid[1][i2];

    phgQuadInterfaceRectangle(ls != NULL ? ls2 : NULL, ls_order, ls2_grad,
			      rect, quad_order, rule_m, rule_0, rule_p);

    /* convert 2D barycentric coordinates to 3D xyz */
    for (k = -1; k <= 1; k++) {
	FLOAT **rule = k < 0 ? rule_m : (k > 0 ? rule_p : rule_0);
	FLOAT *pw, *pw2, *pw3, *nv2, *nv3, *nv;
	if (rule == NULL)
	    continue;
	np = phgQuadInterfaceRuleInfo(*rule, NULL, &pw2, &nv2);
	pw3 = (FLOAT *)phgAlloc((Dim + 1 + (nv2 == NULL ? 0 : Dim))
							* np * sizeof(*pw3));
	nv = nv3 = nv2 == NULL ? NULL : pw3 + (Dim + 1) * np;
	for (i = 0, pw = pw3; i < np; i++, pw += Dim + 1, pw2 += 2 + 1) {
	    rect2xyz(pw2[0], pw2[1], pw, pw + 1, pw + 2);
	    pw[3] = pw2[2];
	    if (nv2 == NULL)
		continue;
	    /* convert to 3D */
	    nv[i0] = 0.0;
	    nv[i1] = nv2[0];
	    nv[i2] = nv2[1];
	    nv2 += 2;
	    nv += 3;
	}
	phgFree(*rule);
	*rule = phgQuadInterfaceRuleCreate(Dim, np, pw3, Dim + 1,
			pw3 + Dim, Dim + 1, 1.0, nv3, Dim);
	phgFree(pw3);
    }

    phgFree(qf_ctx);
    qf_ctx = NULL;
}
