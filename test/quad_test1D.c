/* Test phgQuad1D with:
 * 	\int_0^1 \int_0^{\sin(y)} xy dx dy = -cos(2)/16 - sin(2)/8 + 3/16 */

#include "phg.h"
#include <math.h>

/* Note ctx = &y in f0() */
static void f0(FLOAT x, FLOAT *res, void *y) { *res = x * *((FLOAT *)y); }
static void f(FLOAT y, FLOAT *res, void *ctx)
    { phgQuad1D(f0, 1, 0., Sin(y), 1, res, &y); }

int
main(int argc, char *argv[])
{
    int i, n, order;
    FLOAT y0, y1, res, r, exact;

    if (argc != 3) {
	fprintf(stderr, "Usage: %s n order\n", argv[0]);
	fprintf(stderr, "(n: number of intervals, order: quadrature order)\n");
	return 1;
    }

    n = atoi(argv[1]);
    order = atoi(argv[2]);
    res = 0.0;
    y0 = 0.0; 
    for (i = 0; i < n; i++) {
	y1 = ((FLOAT)(i + 1)) / (FLOAT)n;
	phgQuad1D(f, 1, y0, y1, order, &r, NULL);
	res += r;
	y0 = y1;
    }

    exact = ((3.0 - Cos(2.0)) * 0.5 - Sin(2.0)) * 0.125;
    printf("result = %lg, error = %lg\n", (double)res,
					(double)(1.0 - res / exact));

    return 0;
}
