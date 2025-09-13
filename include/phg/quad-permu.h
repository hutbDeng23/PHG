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

/* $Id: quad-permu.h,v 1.4 2022/08/25 06:59:12 zlb Exp $ */

#ifndef PHG_QUAD_PERMU_H

/*-------------------------- Permutation macros ----------------------------*/

#define PermC(a)	(a)		/* the casting macro */
#define PermOne		_F(1.0)		/* The "1" */
#define PermHalf	_F(0.5)		/* The half */
#define PermThree	_F(3.0)		/* The "3" */

/* 0D */
#define Perm1(a)	_F(a)
#define Dup1(w)		_F(w)

#define Cons1(a)	Perm1(a)

/* unsymmetric, single point orbit */
#define Dup0		Dup1
#define Dup10		Dup0
#define Perm10		Perm1

/* 1D */
#define Perm2(a)	_F(a),_F(a)
#define Dup2(w)		_F(w)
#define Perm11(a)	PermC(_F(a)),PermC(PermOne-(_F(a))), \
			PermC(PermOne-(_F(a))),PermC(_F(a))
#define Dup11(w)	_F(w),_F(w)

#define Cons2(a)	Perm2(a)
#define Cons11(a)	Perm11(a)

/* unsymmetric, single point orbit */
#define Dup20		Dup0
#define Perm20(a)	PermC(_F(a)),PermC(PermOne-(_F(a)))

/* 2D */
#define Perm3(a)	_F(a),_F(a),_F(a)
#define Dup3(w)		_F(w)
#define Perm21(a) \
	PermC(PermOne-(_F(a))-(_F(a))),PermC(_F(a)),PermC(_F(a)), \
	PermC(_F(a)),PermC(PermOne-(_F(a))-(_F(a))),PermC(_F(a)), \
	PermC(_F(a)),PermC(_F(a)),PermC(PermOne-(_F(a))-(_F(a)))
#define Dup21(w)	_F(w),_F(w),_F(w)
#define Perm111(a,b) \
	PermC(_F(a)),PermC(_F(b)),PermC(PermOne-(_F(a))-(_F(b))), \
	PermC(_F(a)),PermC(PermOne-(_F(a))-(_F(b))),PermC(_F(b)), \
	PermC(_F(b)),PermC(_F(a)),PermC(PermOne-(_F(a))-(_F(b))), \
	PermC(_F(b)),PermC(PermOne-(_F(a))-(_F(b))),PermC(_F(a)), \
	PermC(PermOne-(_F(a))-(_F(b))),PermC(_F(a)),PermC(_F(b)), \
	PermC(PermOne-(_F(a))-(_F(b))),PermC(_F(b)),PermC(_F(a))
#define Dup111(w)	_F(w),_F(w),_F(w),_F(w),_F(w),_F(w)

#define Cons3(a)	Perm3(a)
#define Cons21(a)	Perm21(a)
#define Cons111(a,b)	Perm111(a,b)

/* unsymmetric, single point orbit */
#define Dup30		Dup0
#define Perm30(a,b)	PermC(_F(a)),PermC(_F(b)),PermC(PermOne-(_F(a))-(_F(b)))

/* 3D */
#define Perm4(a)	_F(a),_F(a),_F(a),_F(a)
#define Dup4(w)		_F(w)
#define Perm31(a) \
    PermC(PermOne-PermThree*(_F(a))),PermC(_F(a)),PermC(_F(a)),PermC(_F(a)), \
    PermC(_F(a)),PermC(PermOne-PermThree*(_F(a))),PermC(_F(a)),PermC(_F(a)), \
    PermC(_F(a)),PermC(_F(a)),PermC(PermOne-PermThree*(_F(a))),PermC(_F(a)), \
    PermC(_F(a)),PermC(_F(a)),PermC(_F(a)),PermC(PermOne-PermThree*(_F(a)))
#define Dup31(w)	_F(w),_F(w),_F(w),_F(w)
#define Perm22(a) \
    PermC(_F(a)),PermC(_F(a)),PermC(PermHalf-(_F(a))),PermC(PermHalf-(_F(a))), \
    PermC(_F(a)),PermC(PermHalf-(_F(a))),PermC(_F(a)),PermC(PermHalf-(_F(a))), \
    PermC(_F(a)),PermC(PermHalf-(_F(a))),PermC(PermHalf-(_F(a))),PermC(_F(a)), \
    PermC(PermHalf-(_F(a))),PermC(_F(a)),PermC(PermHalf-(_F(a))),PermC(_F(a)), \
    PermC(PermHalf-(_F(a))),PermC(_F(a)),PermC(_F(a)),PermC(PermHalf-(_F(a))), \
    PermC(PermHalf-(_F(a))),PermC(PermHalf-(_F(a))),PermC(_F(a)),PermC(_F(a))
#define Dup22(w)	_F(w),_F(w),_F(w),_F(w),_F(w),_F(w)
#define Perm211(a,b) \
	PermC(_F(a)),PermC(_F(a)),PermC(_F(b)),\
		PermC(PermOne-(_F(a))-(_F(a))-(_F(b))), \
	PermC(_F(a)),PermC(_F(a)),PermC(PermOne-(_F(a))-(_F(a))-(_F(b))),\
		PermC(_F(b)), \
	PermC(_F(a)),PermC(_F(b)),PermC(_F(a)),\
		PermC(PermOne-(_F(a))-(_F(a))-(_F(b))), \
	PermC(_F(a)),PermC(_F(b)),PermC(PermOne-(_F(a))-(_F(a))-(_F(b))),\
		PermC(_F(a)), \
	PermC(_F(a)),PermC(PermOne-(_F(a))-(_F(a))-(_F(b))),PermC(_F(a)),\
		PermC(_F(b)), \
	PermC(_F(a)),PermC(PermOne-(_F(a))-(_F(a))-(_F(b))),PermC(_F(b)),\
		PermC(_F(a)), \
	PermC(_F(b)),PermC(_F(a)),PermC(_F(a)),\
		PermC(PermOne-(_F(a))-(_F(a))-(_F(b))), \
	PermC(_F(b)),PermC(_F(a)),PermC(PermOne-(_F(a))-(_F(a))-(_F(b))),\
		PermC(_F(a)), \
	PermC(_F(b)),PermC(PermOne-(_F(a))-(_F(a))-(_F(b))),PermC(_F(a)),\
		PermC(_F(a)), \
	PermC(PermOne-(_F(a))-(_F(a))-(_F(b))),PermC(_F(a)),PermC(_F(a)),\
		PermC(_F(b)), \
	PermC(PermOne-(_F(a))-(_F(a))-(_F(b))),PermC(_F(a)),PermC(_F(b)),\
		PermC(_F(a)), \
	PermC(PermOne-(_F(a))-(_F(a))-(_F(b))),PermC(_F(b)),PermC(_F(a)),\
		PermC(_F(a))
#define Dup211(w)	_F(w),_F(w),_F(w),_F(w),_F(w),_F(w),\
			_F(w),_F(w),_F(w),_F(w),_F(w),_F(w)
#define Perm0111(p,a,b,c) \
    PermC(p),PermC(a),PermC(b),PermC(c),PermC(p),PermC(a),PermC(c),PermC(b), \
    PermC(p),PermC(b),PermC(a),PermC(c),PermC(p),PermC(b),PermC(c),PermC(a), \
    PermC(p),PermC(c),PermC(a),PermC(b),PermC(p),PermC(c),PermC(b),PermC(a)
#define Perm1111(a,b,c) \
	Perm0111(_F(a),_F(b),_F(c),PermOne-(_F(a))-(_F(b))-(_F(c))), \
	Perm0111(_F(b),_F(a),_F(c),PermOne-(_F(a))-(_F(b))-(_F(c))), \
	Perm0111(_F(c),_F(a),_F(b),PermOne-(_F(a))-(_F(b))-(_F(c))), \
	Perm0111(PermOne-(_F(a))-(_F(b))-(_F(c)),_F(a),_F(b),_F(c))
#define Dup1111(w)	Dup111(w), Dup111(w), Dup111(w), Dup111(w)

#define Cons4(a)	Perm4(a)
#define Cons31(a)	Perm31(a)
#define Cons22(a)	Perm22(a)
#define Cons211(a,b)	Perm211(a,b)
#define Cons1111(a,b,c)	Perm111(a,b,c)

/* unsymmetric, single point orbit */
#define Dup40		Dup0
#define Perm40(a,b,c) \
  PermC(_F(a)),PermC(_F(b)),PermC(_F(c)),PermC(PermOne-(_F(a))-(_F(b))-(_F(c)))

#define PHG_QUAD_PERMU_H
#endif
