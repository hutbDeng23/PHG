/* $Id: float_test.c,v 1.1 2011/04/14 05:27:09 zlb Exp $ */

#include "phg.h"

#define show_type(type, type_name, size) {				\
    type eps, one = 1.0L;						\
    eps = one;								\
    while (one + eps != one)						\
	eps *= 0.5L;							\
    printf("%s: size=%d, eps=%le\n", type_name, size, (double)eps);	\
}

int
main(int argc, char *argv[])
{

    show_type(float, "float", SIZEOF_FLOAT)
    show_type(double, "double", SIZEOF_DOUBLE)
    show_type(long double, "long double", SIZEOF_LONG_DOUBLE)
#if SIZEOF___FLOAT128 > 0
    show_type(__float128, "__float128", SIZEOF___FLOAT128)
#endif

    return 0;
}
