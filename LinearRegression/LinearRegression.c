/* Andre Augusto Giannotti Scota (https://sites.google.com/view/a2gs/) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
/*
#include <unistd.h>
#include <signal.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
*/

#define LR_OK   (0)
#define LR_ERRO (1)

float angularCoefficientLR(unsigned int n, float sumXY, float sumX, float sumY, float sumXsquare)
{
	return(((n * sumXY) - (sumX * sumY)) / ((n * sumXsquare) - pow(sumX, 2.0)));
}

float interceptLR(float sumY, float angulCoef, float sumX, unsigned int n)
{
	return(((sumY - (angulCoef * sumX)) / n));
}

void summarization_LinearRegression(float *sumX, float *sumY, float *sumXY, float *sumXsquare, float *sumYsquare, float points[][2], unsigned int n)
{
	unsigned int i = 0;

	for(i = 0, *sumX  = 0.0, *sumY  = 0.0, *sumXY = 0.0, *sumXsquare = 0.0, *sumYsquare = 0.0; i < n; i++){
		*sumX += points[i][0];
		*sumY += points[i][1];

		*sumXsquare += powf(points[i][0], 2.0);
		*sumYsquare += powf(points[i][1], 2.0);

		*sumXY += points[i][0] * points[i][1];
	}

	return;
}

int LinearRegression(float points[][2], const unsigned int n, float *intercept, float *angulCoef)
{
	float sumX  = 0.0;
	float sumY  = 0.0;
	float sumXY = 0.0;
	float sumXsquare = 0.0;
	float sumYsquare = 0.0;

	summarization_LinearRegression(&sumX, &sumY, &sumXY, &sumXsquare, &sumYsquare, points, n);

	*angulCoef = angularCoefficientLR(n, sumXY, sumX, sumY, sumXsquare);
	*intercept = interceptLR(sumY, *angulCoef, sumX, n);

	return(LR_OK);
}

float CorrelationCoefficient_LinearRegression(float points[][2], const unsigned int n)
{
	float sumX  = 0.0;
	float sumY  = 0.0;
	float sumXY = 0.0;
	float sumXsquare = 0.0;
	float sumYsquare = 0.0;

	summarization_LinearRegression(&sumX, &sumY, &sumXY, &sumXsquare, &sumYsquare, points, n);

	return(((n * sumXY) - (sumX * sumY)) / (sqrtf(((n * sumXsquare) - (pow(sumX, 2.0)))) * sqrtf(((n * sumYsquare) - (pow(sumY, 2.0))))));
}

int main(int argc, char *argv[])
{

	float m[4][2] = {{3.0, 7.0}, {2.0, 5.0}, {-1.0, -1.0}, {4.0, 9.0}};
	int n = 4;
	float a = 0.0, b = 0.0;
	float cc = 0.0;

	printf("Sample:\n");
	printf("   X  |  Y\n");
	printf("------+------\n");
	printf("%02.02f  | %02.02f\n", m[0][0], m[0][1]);
	printf("%02.02f  | %02.02f\n", m[1][0], m[1][1]);
	printf("%02.02f | %02.02f\n", m[2][0], m[2][1]);
	printf("%02.02f  | %02.02f\n\n", m[3][0], m[3][1]);

	cc = CorrelationCoefficient_LinearRegression(m, n);

	if(LinearRegression(m, n, &a, &b) == LR_ERRO){
		printf("Linear Regression error.\n");
		return(-1);
	}

	printf("f(x) = [%02.02f]x + [%02.02f]\n", a, b);
	printf("Correlation coefficient: [%02.02f]\n", cc);

	return(0);
}
