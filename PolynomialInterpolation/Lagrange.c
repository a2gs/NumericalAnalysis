/* Andre Augusto Giannotti Scota (https://sites.google.com/view/a2gs/) */

#include <stdio.h>

float PolynomialInterpolation_Lagrange(float XY[][2], unsigned int n, float x)
{
	unsigned int ip = 0;
	unsigned int is = 0;
	float P = 0.0;
	float S = 0.0;

	for(is = 0, S = 0.0; is < n; is++){

		for(ip = 0, P = 1.0; ip < n; ip++){
			if(ip == is) continue;

			P *= (x - XY[ip][0]) / (XY[is][0] - XY[ip][0]);
		}

		S += XY[is][1] * P;
	}

	return(S);
}

int main(int argc, char *argv[])
{
	float m[7][2] = {{1.0, 2.0}, {2.0, 3.0}, {3.0, 4.0}, {4.0, 5.0}, {5.0, 6.0}, {6.0, 7.0}, {7.0, 8.0}};
	unsigned int n = 7;
	float x = 3.5;
	float y = 0.0;

	y = PolynomialInterpolation_Lagrange(m, n, x);

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

	return(0);
}
