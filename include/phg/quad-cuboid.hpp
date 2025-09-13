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

/* High order numerical integration in subdomains of a hyperrectangle cut by
 * an implicit interface.
 *
 * Original Author: WEN Yundi (wenyd@lsec.cc.ac.cn)
 *
 * $Id: quad-cuboid.hpp,v 1.39 2021/02/16 03:40:30 zlb Exp $ */

#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <set>
#include <vector>
#include <cstring>

#include "phg.h"
#include "phg/sturm.h"

#if 0
using namespace
    std::placeholders;
const int
    MAX_ORDER = 10;
#endif

typedef std::function <std::vector <FLOAT>(std::vector <FLOAT>)> FUNC;
typedef std::initializer_list <FUNC> FUNCLIST;

#if FT_PHG == FT___FLOAT128	// libquadmath may not be available
std::ostream &operator<<(std::ostream &os, __float128 x)
{
    char buffer[128];
    quadmath_snprintf(buffer, 128, "%Qe", x);
    os << buffer;
    return os;
}
#endif

std::ostream &operator<<(std::ostream &os, std::vector <FLOAT> v)
{
    os << "[";
    for (int i = 0; i < v.size(); i++) {
	if (i > 0) os << ",";
	os << v[i];
    }
    os << "]";
    return os;
}

struct quad_struct {
    FLOAT mW;
    std::vector <FLOAT> mAbs;
    quad_struct(FLOAT w, std::vector <FLOAT>abs) {
	mW = w;
	mAbs = abs;
    } 
    quad_struct(const quad_struct &src) = default;
    quad_struct(quad_struct &&src) {
	mW = src.mW;
	mAbs = std::move(src.mAbs);
    }
    quad_struct &operator=(const quad_struct &src) = default;
};

std::ostream &operator<<(std::ostream &os, const quad_struct &q)
{
    os << q.mW << ": (" << q.mAbs[0];
    for (unsigned i = 1; i < q.mAbs.size(); i++) {
	os << ", " << q.mAbs[i];
    }
    os << ")\n";
    return os;
}

class quad_cuboid {
  public:
    // user interface
    // constructors
    quad_cuboid(unsigned dim,
		std::vector <FLOAT>lower, std::vector <FLOAT>upper);
    quad_cuboid(unsigned dim,
		std::vector <FLOAT>lower, std::vector <FLOAT>upper,
		int lineno);
    ~quad_cuboid() = default;

    // initializations

    /* Single levelset function initialization.
     * Preset integrals will be enabled after this call.
     * parameters:
     * L:       levelset function
     * gradL:   the gradient of L
     * deg:     >= 0 => L is a polynomial of degree 'deg',
     *		< 0  => L is non polynomial, interpolations with polynomials
     *			of degree '-deg' will be used in root finding.
     * order:   order of the Gauss-Legendre qudrature
     */
    void initialize(FUNC L, FUNC gradL, int deg, unsigned order);

    // quaddrature strutures
    std::vector <quad_struct> mPosQuad;		// quad points in Omega+
    std::vector <quad_struct> mNegQuad;		// quad points in Omega-
    std::vector <quad_struct> mSurfQuad;	// quad points on the interface
    std::vector <std::vector<FLOAT>>
			      mSurfNormVecs;	// list of unit normal vectors

  private:
    std::vector <quad_struct> mQuad;		// quad points in all regions

    /* Multi levelset function initialization.
     * Won't enable preset integrals.
     * May cause unexpected errors when called from outside the quad_cuboid
     * class.
     */
    void initialize(const std::vector <FUNC> &L,
		    const std::vector <FUNC> &gradL,
		    int deg, unsigned order);
    void split(const std::vector<FUNC> &L, const std::vector<FUNC> &gradL,
	       const char *file, int line);
    #define SPLIT(L, gradL) split(L, gradL, __FILE__, __LINE__)
    FLOAT normalize(std::vector<FLOAT> &g);	// return value = original norm

    // data members
 
    FLOAT threshold = 0.312245;	// \approx Sqrt(1.0 - 0.95*0.95)

    int mLine;		// source file line no. (for debugging)
    unsigned mDim, mLevel = 0;
    std::vector <FLOAT>mLower, mUpper;
    int mDeg;
    // order of 1D quadrature points.
    unsigned mOrder = 0;
    // first integration direction
    unsigned mDirection = 0;
    int mInitialized = 0;
#define MAX_DIM 3
    struct {
	FLOAT lower, upper;
	int dir;
    } mInfo[MAX_DIM];	// information about the upper levels

    bool refinement(const FUNC &F, const FUNC &gradF, FLOAT x1, FLOAT x2,
		    FLOAT *root);
    std::vector<FLOAT> solver(const FUNC &F, const FUNC &gradF,
						FLOAT a, FLOAT b);
};
