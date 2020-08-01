/* Andre Augusto Giannotti Scota (https://sites.google.com/view/a2gs/) */

#include <stdio.h>

/* #undef DEBUG */
/* #define FLOAT_t float */

#define LAGRANGE_PI_OK    (0)
#define LAGRANGE_PI_ERROR (1)

int PolynomialInterpolation_Lagrange(FLOAT_t XY[][2], unsigned int n, FLOAT_t x, FLOAT_t *y, char *poly, size_t polySz)
{
	register unsigned int ip = 0;
	register unsigned int is = 0;
	register FLOAT_t P = 0.0;
	size_t pSleft = 0;

	for(is = 0, *y = 0.0; is < n; is++){

		for(ip = 0, P = 1.0; ip < n; ip++){
			if(ip == is) continue;

			P *= (x - XY[ip][0]) / (XY[is][0] - XY[ip][0]);

			if(poly != NULL) pSleft += snprintf(&poly[pSleft], polySz - pSleft, "((x - %E) / (%E - %E)) * ", XY[ip][0], XY[is][0], XY[ip][0]);
		}

		if(poly != NULL) pSleft += snprintf(&poly[pSleft], polySz - pSleft, "%E%s", XY[0][is], (is + 1 < n ? " + " : ""));

		*y += XY[is][1] * P;
	}

	return(LAGRANGE_PI_OK);
}

int main(int argc, char *argv[])
{
	FLOAT_t m[7][2] = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}};
	unsigned int n = 7;
	FLOAT_t x = 3.5;
	FLOAT_t y = 0.0;
#ifdef DEBUG
	#define POLYSTR_SZ (4000)
	char p[POLYSTR_SZ] = {0};
#endif

#ifdef DEBUG
	if(PolynomialInterpolation_Lagrange(m, n, x, &y, p, POLYSTR_SZ) == LAGRANGE_PI_ERROR){
#else
	if(PolynomialInterpolation_Lagrange(m, n, x, &y, NULL, 0) == LAGRANGE_PI_ERROR){
#endif
		return(-1);
	}

	printf("SAMPLE:\n");
	printf("   X  |  Y\n");
	printf("------+------\n");
	printf(" %02.02f | %02.02f\n",   m[0][0], m[0][1]);
	printf(" %02.02f | %02.02f\n",   m[1][0], m[1][1]);
	printf(" %02.02f | %02.02f\n",   m[2][0], m[2][1]);
	printf(" %02.02f | %02.02f\n",   m[3][0], m[3][1]);
	printf(" %02.02f | %02.02f\n",   m[4][0], m[4][1]);
	printf(" %02.02f | %02.02f\n",   m[5][0], m[5][1]);
	printf(" %02.02f | %02.02f\n\n", m[6][0], m[6][1]);

	printf("f(%02.02f) = [%02.03f]\n", x, y);

#ifdef DEBUG
	printf("f(x) = [%s]\n", p);
#endif

	return(0);
}
