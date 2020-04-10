/* Andre Augusto Giannotti Scota (https://sites.google.com/view/a2gs/) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double determinant3x3(const double m[3][3])
{
	return((m[0][0] * m[1][1] * m[2][2]) +
	       (m[0][1] * m[1][2] * m[2][0]) +
	       (m[0][2] * m[1][0] * m[2][1]) -
	       (m[0][2] * m[1][1] * m[2][0]) -
	       (m[0][1] * m[1][0] * m[2][2]) -
	       (m[0][0] * m[1][2] * m[2][1]));
}

double determinant2x2(const double m[2][2])
{
	return((m[0][0] * m[1][1]) -
	       (m[0][1] * m[1][0]));
}

unsigned int offset(unsigned int l, unsigned int c, unsigned int cTot)
{
	return((l * cTot) + c);
}

void printLinearSystem(const unsigned int dim, const double *q)
{
	register unsigned int i = 0, j = 0;
	char var = 'a';

	for(i = 0; i < dim; i++){

		for(j = 0, var = 'a'; j < dim; j++, var++){
			printf("%f%c", q[offset(i, j, dim+1)], var);

			if(j != dim-1) printf(" + ");
			else           printf(" = ");
		}

		printf("%f\n", q[offset(i, j, dim+1)]);
	}

	return;
}

#define GE_RET_UNDETERMINED (1)
#define GE_RET_IMPOSSIBLE   (2)
#define GE_RET_POSSIBLE     (3)
#define GE_RET_ERROR        (-1)

int GaussElimination(unsigned int dim, double *q, double *result)
{
	register unsigned int round = 0;
	register unsigned int line  = 0;
	register unsigned int row   = 0;
	register unsigned int i = 0;
	register unsigned int mtrxRows = 0;
	double det[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
	double *newLine = NULL;
	size_t SzNewLine = 0;

	mtrxRows = dim + 1;

	SzNewLine = mtrxRows * sizeof(double); /* +1 for result */

	newLine = (double *)malloc(SzNewLine);
	if(newLine == NULL)
		return(GE_RET_ERROR);

	printf("--------------------------------------------------\nTriangulation:\n");

	for(round = 0; round < dim; round++){

		for(line = round + 1; line < mtrxRows; line++){

			memset(newLine, 0, SzNewLine);
			det[0][0] = det[0][1] = det[1][0] = det[1][1] = 0.0;

			for(row = round; row < mtrxRows; row++){

				det[0][0] = q[offset(round, round, mtrxRows)];
				det[0][1] = q[offset(round, row+1, mtrxRows)];
				det[1][0] = q[offset(line , round, mtrxRows)];
				det[1][1] = q[offset(line , row+1, mtrxRows)];

				newLine[row+1] = determinant2x2(det);
			}

			for(i = 0; i < mtrxRows; i++)
				q[offset(line, i, mtrxRows)] = newLine[i];
		}

		printf("\n");
		printLinearSystem(dim, q);
	}

	free(newLine);

	printf("--------------------------------------------------\nSolving:\n");








	result[0] = 1.0;
	result[1] = 2.0;
	result[2] = 3.0;



	return(GE_RET_POSSIBLE);
}

int main(int argc, char *argv[])
{
	int ret = 0;
	double *q = NULL, *result = NULL;
	unsigned int dim = 0, i = 0;
	size_t len = 0;

	/* --- sample data START ---------------------- */
	/*

  /-
 | 2a + 3b - 1c = 5
-| 1a - 1b + 2c = 5
 | 1a + 4b - 1c = 6
  \-

  a = 1
  b = 2
  c = 3

	 */

	dim = 3; /* 3 variables with 3 equations */

	len = ((dim * dim) + dim ) * sizeof(double);
	q = (double *)malloc(len);
	if(q == NULL){
		printf("Erro malloc q.\n");
		return(-1);
	}
	memset(q, 0, len);

	len = dim * sizeof(double);
	result = (double *)malloc(len);
	if(result == NULL){
		printf("Erro malloc result.\n");
		return(-1);
	}
	memset(result, 0, len);

	q[0] = 2.0; q[1] =  3.0; q[2]  = -1.0; q[3]  = 5.0;
	q[4] = 1.0; q[5] = -1.0; q[6]  =  2.0; q[7]  = 5.0;
	q[8] = 1.0; q[9] =  4.0; q[10] = -1.0; q[11] = 6.0;
	/* --- sample data END ------------------------ */

	printLinearSystem(dim, q);

	ret = GaussElimination(dim, q, result);

	switch(ret){
		case GE_RET_UNDETERMINED:
			printf("Undetermined linear system.\n");
			break;

		case GE_RET_IMPOSSIBLE:
			printf("Impossible linear system.\n");
			break;

		case GE_RET_POSSIBLE:
			for(i = 0; i < dim; i++)
				printf("%c = %f\n", i + 'a', result[i]);

			break;

		case GE_RET_ERROR:
			break;

		default:
			break;
	}

	free(q);
	free(result);

	return(0);
}
