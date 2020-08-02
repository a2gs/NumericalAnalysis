/* Andre Augusto Giannotti Scota (https://sites.google.com/view/a2gs/) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* #undef DEBUG */
/* #define FLOAT_t float */
/* #define FLOAT_t double */

#define NEWTON_PI_OK    (0)
#define NEWTON_PI_ERROR (1)

static inline unsigned int offset(unsigned int l, unsigned int c, unsigned int cTot)
{
	return((l * cTot) + c);
}

int DividedDifferenceTableCalculation(FLOAT_t *dd, unsigned int n, unsigned int c)
{
	register unsigned int i  = 0;
	register unsigned int j  = 0;
	register unsigned int nn = 0;

	for(j = 2, nn = n - 1; j < c; j++, nn--){

		for(i = 0; i < nn; i++){

			dd[offset(i, j, c)] = ((dd[offset(i + 1, j - 1, c)] - dd[offset(i, j - 1, c)])
			                                                    /
			                       (dd[offset(i + j - 1, 0, c)] - dd[offset(i, 0, c)])    );

#ifdef DEBUG
			{
				unsigned int ib;
				unsigned int jb;

				printf("[%02d ; %02d] = (%02.02f - %02.02f) / (%02.02f - %02.02f) = %02.02f\n", i, j,
				       dd[offset(i+1, j-1, c)], dd[offset(i, j-1, c)], dd[offset(i+j-1, 0, c)], dd[offset(i, 0, c)], dd[offset(i, j, c)]);

				for(jb = 0; jb < c; jb++)
					printf("\t%02d\t", jb);

				printf("\n");

				for(ib = 0; ib < n; ib++){
					printf("%02d  ", ib);
					for(jb = 0; jb < c; jb++)
						printf("%02.02f\t\t", dd[offset(ib, jb, c)]);
					printf("\n");
				}
				printf("--------------------------\n");
			}
#endif

		}
	}

	return(NEWTON_PI_OK);
}

int DividedDifferenceTable(FLOAT_t XY[][2], unsigned int n, FLOAT_t **dd, unsigned int *columns)
{
	register unsigned int i = 0;
	register FLOAT_t *pdd = NULL;
	size_t sdd = 0;

	*columns = n + 1;
	sdd = sizeof(FLOAT_t) * n * *columns;

	pdd = malloc(sdd);
	if(pdd == NULL)
		return(NEWTON_PI_ERROR);

	memset(pdd, 0, sdd);

	for(i = 0; i < n; i++){
		pdd[offset(i, 0, *columns)] = XY[i][0];
		pdd[offset(i, 1, *columns)] = XY[i][1];
	}

	if(DividedDifferenceTableCalculation(pdd, n, *columns) == NEWTON_PI_ERROR){
		free(pdd);
		return(NEWTON_PI_ERROR);
	}

	*dd = pdd;

	return(NEWTON_PI_OK);
}

int PolynomialInterpolation_Newton(FLOAT_t XY[][2], unsigned int n, FLOAT_t x, FLOAT_t *y, char *poly, size_t polySz)
{
	register unsigned int ip = 0;
	register unsigned int is = 0;
	register FLOAT_t P = 0.0;
	unsigned int c = 0;
	size_t pSleft = 0;
	FLOAT_t *dd = NULL;

	if(DividedDifferenceTable(XY, n, &dd, &c) == NEWTON_PI_ERROR)
		return(NEWTON_PI_ERROR);

	if(poly != NULL) pSleft += snprintf(&poly[pSleft], polySz - pSleft, "%E + ", dd[offset(0, 1, c)]);

	for(is = 2, *y = dd[offset(0, 1, c)]; is < c; is++){

		for(ip = 0, P = 1.0; ip < is - 1; ip++){
			P *= (x - dd[offset(ip, 0, c)]);

			if(poly != NULL) pSleft += snprintf(&poly[pSleft], polySz - pSleft, "(x - %E) * ", dd[offset(ip, 0, c)]);
		}

		if(poly != NULL) pSleft += snprintf(&poly[pSleft], polySz - pSleft, "%E%s", dd[offset(0, is, c)], (is + 1 < c ? " + " : ""));

		*y += dd[offset(0, is, c)] * P;
	}

	free(dd);

	return(NEWTON_PI_OK);
}

int main(int argc, char *argv[])
{
	FLOAT_t m[7][2] = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}};
	unsigned int n = 7;
	FLOAT_t x = 3.5;
	FLOAT_t y = 0.0;
#ifdef DEBUG
	#define POLYSTR_SZ (2000)
	char p[POLYSTR_SZ] = {0};
#endif

#ifdef DEBUG
	if(PolynomialInterpolation_Newton(m, n, x, &y, p, POLYSTR_SZ) == NEWTON_PI_ERROR){
#else
	if(PolynomialInterpolation_Newton(m, n, x, &y, NULL, 0) == NEWTON_PI_ERROR){
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
