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

/*
 * Runtime evaluation of analytical functions (expressions).
 *
 * $Id: expr.c,v 1.10 2014/04/04 06:47:12 zlb Exp $ */

#include "phg/config.h"	/* get the value of HAVE_LIBMATHEVAL */

#include <stdlib.h>
#if !HAVE_LIBMATHEVAL
# include "expression.c"	/* use expression.y as a fallback */
#else	/* !HAVE_LIBMATHEVAL */

#include "phg/expression.h"

/* use libmatheval */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <matheval.h>

EXPR *
phgDefine3DFunction(const char *function_expression)
{
    EXPR *e;
    char *str, *p;
    const char *q;

    /* delete while spaces in the expression */
    str = malloc(strlen(function_expression) + 1);
    for (p = str, q = function_expression; *q != '\0'; q++) {
	if (isspace((int)(*q)))
	    continue;
	*(p++) = *q;
    }
    *p = '\0';
    e = evaluator_create(str);
    free(str);

    return e;
}

double
phgEvaluate3DFunction(EXPR *tree, double x, double y, double z, double t)
{
    char *names[] = {"x", "y", "z", "t"};
    double values[] = {x, y, z, t};

    if (tree == NULL) {
	fprintf(stderr, "phgEvaluate3DFunction: function undefined\n");
	return 0.0;
    }

    return evaluator_evaluate(tree, 4, names, values);
}

EXPR *
phgDup3DFunction(EXPR *tree)
{
    return evaluator_create(evaluator_get_string(tree));
}

void
phgFree3DFunction(EXPR *tree)
{
    evaluator_destroy(tree);
}

char *
phgDump3DFunction(EXPR *tree)
{
    return evaluator_get_string(tree);
}

/*---------------------------------------------------------------------------*/
#endif	/* !HAVE_LIBMATHEVAL */
