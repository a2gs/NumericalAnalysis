/* Andre Augusto Giannotti Scota (https://sites.google.com/view/a2gs/) */

#include <stdio.h>
#include <string.h>

#define GE_STATIC (1)
//#define GE_DYNAMIC (1)

static inline double determinant3x3(const double m[3][3])
{
	return((m[0][0] * m[1][1] * m[2][2]) +
	       (m[0][1] * m[1][2] * m[2][0]) +
	       (m[0][2] * m[1][0] * m[2][1]) -
	       (m[0][2] * m[1][1] * m[2][0]) -
	       (m[0][1] * m[1][0] * m[2][2]) -
	       (m[0][0] * m[1][2] * m[2][1]));
}

static inline double determinant2x2(const double m[2][2])
{
	return((m[0][0] * m[1][1]) -
	       (m[0][1] * m[1][0]));
}

static inline unsigned int offset(unsigned int l, unsigned int c, unsigned int cTot)
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

unsigned int countZeros(unsigned int totElem, double *q)
{
	unsigned int i = 0;

	for(i = 0; (q[i] == 0.0) && (i < totElem); i++);

	return(i);
}

/* return:
 * 0 - there was no exchange
 * 1 - there was exchange
 */
int reorder(unsigned int dim, double *q)
{
	int flag = 0;
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int k = 0;
	unsigned int qtd0Current = 0;
	unsigned int qtd0Next    = 0;
	double aux = 0.0;

	flag = 0;
	k    = dim + 1;

	for(i = 0; i < dim-1; i++){
		qtd0Current = countZeros(dim, &q[offset(i  , 0, k)]);
		qtd0Next    = countZeros(dim, &q[offset(i+1, 0, k)]);

		if(qtd0Current > qtd0Next){

			for(j=0; j <= dim; j++){
				aux = q[offset(i, j, k)];
				q[offset(i  , j, k)] = q[offset(i+1, j, k)];
				q[offset(i+1, j, k)] = aux;
			}

			flag = 1;
		}
	}

	return(flag);
}

#define GE_RET_UNDETERMINED (1)
#define GE_RET_IMPOSSIBLE   (2)
#define GE_RET_POSSIBLE     (3)
#define GE_RET_ERROR        (-1)
#define GE_RET_OK           (0)

int GaussElimination_Triangulation(unsigned int dim, double *q)
{
	register unsigned int round    = 0;
	register unsigned int line     = 0;
	register unsigned int row      = 0;
	register unsigned int i        = 0;
	register unsigned int mtrxRows = 0;
	size_t SzNewLine               = 0;
	double det[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
#ifdef GE_STATIC
	double newLine[dim+1];
#elif GE_DYNAMIC
	double *newLine = NULL;
#endif

	mtrxRows = dim + 1; /* +1 for result: q[dim][mtrxRows] */
	SzNewLine = mtrxRows * sizeof(double);

#ifdef GE_DYNAMIC
	newLine = (double *)malloc(SzNewLine);
	if(newLine == NULL)
		return(GE_RET_ERROR);
#endif

	for(round = 0; round < dim-1; round++){

		for(line = round+1; line < dim; line++){

			memset(newLine, 0, SzNewLine);
			det[0][0] = det[0][1] = det[1][0] = det[1][1] = 0.0;

			for(row = round+1; row < mtrxRows; row++){

				det[0][0] = q[offset(round, round, mtrxRows)];
				det[0][1] = q[offset(round, row  , mtrxRows)];
				det[1][0] = q[offset(line , round, mtrxRows)];
				det[1][1] = q[offset(line , row  , mtrxRows)];

				newLine[row] = determinant2x2(det);
			}

			for(i = 0; i < mtrxRows; i++)
				q[offset(line, i, mtrxRows)] = newLine[i];
		}

		reorder(dim, q);

		printf("\n");
		printLinearSystem(dim, q);
	}

#ifdef GE_DYNAMIC
	free(newLine);
#endif

	return(GE_RET_OK);
}

#define GE_ZERO_CLOSED_POSITIVE_CLOSED_LIMIT(__GE_point_, __GE_upper_limit_) ((__GE_point_ >= 0) && (__GE_point_ <= __GE_upper_limit_))
int GaussElimination_Solve(unsigned int dim, double *q, double *result)
{
	register unsigned int mtrxRows = 0;
	register unsigned int line     = 0;
	register unsigned int row      = 0;
	register unsigned int limit    = 0;

	mtrxRows = dim + 1; /* +1 for result: q[dim][mtrxRows] */
	limit    = dim - 1;

	if((q[offset(limit, limit, mtrxRows)] == 0.0) && (q[offset(limit, dim, mtrxRows)] == 0.0))
		return(GE_RET_UNDETERMINED);

	if(q[offset(limit, limit, mtrxRows)] == 0.0)
		return(GE_RET_IMPOSSIBLE);

	for(line = limit; GE_ZERO_CLOSED_POSITIVE_CLOSED_LIMIT(line, limit); line--){

		for(row = limit; GE_ZERO_CLOSED_POSITIVE_CLOSED_LIMIT(row, limit); row--){

			if(row == line){
				result[line] = q[offset(line, dim, mtrxRows)] / q[offset(line, row, mtrxRows)];
				break; /* end work at this line */
			}else{
				q[offset(line, dim, mtrxRows)] +=- (q[offset(line, row, mtrxRows)] * result[row]);
			}
		}
	}

	return(GE_RET_POSSIBLE);
}

int GaussElimination(unsigned int dim, double *q, double *result)
{
	printf("--------------------------------------------------\nTriangulation:\n");

	if(GaussElimination_Triangulation(dim, q) == GE_RET_ERROR){
		return(GE_RET_ERROR);
	}

	printf("--------------------------------------------------\nSolving:\n");

	return(GaussElimination_Solve(dim, q, result));
}

int main(int argc, char *argv[])
{
	int ret = 0;
	unsigned int dim = 0, i = 0;
	double *q = NULL;

#ifdef GE_STATIC
	double qq[4][5] = {{0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0, 0.0}}, result[4] = {0.0, 0.0, 0.0, 0.0};
#elif GE_DYNAMIC
	double *result = NULL;
	size_t  len = 0;
#endif

	dim = 3; /* 3 variables with 3 equations: samples 1, 2, 4 and 5 */

#ifdef GE_DYNAMIC

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

#elif GE_STATIC

	q = *qq;

#endif

#define GE_SAMPLE (1)
#if GE_SAMPLE == 1
/* SAMPLE1
  /-
 | 2a + 3b - 1c = 5
-| 1a - 1b + 2c = 5
 | 1a + 4b - 1c = 6
  \-

	a = 1; b = 2; c = 3
*/
	q[offset(0, 0, 4)] = 2.0; q[offset(0, 1, 4)] =  3.0; q[offset(0, 2, 4)] = -1.0; q[offset(0, 3, 4)] = 5.0;
	q[offset(1, 0, 4)] = 1.0; q[offset(1, 1, 4)] = -1.0; q[offset(1, 2, 4)] =  2.0; q[offset(1, 3, 4)] = 5.0;
	q[offset(2, 0, 4)] = 1.0; q[offset(2, 1, 4)] =  4.0; q[offset(2, 2, 4)] = -1.0; q[offset(2, 3, 4)] = 6.0;
#elif GE_SAMPLE == 2
	q[offset(0, 0, 4)] = 3.0; q[offset(0, 1, 4)] = -1.0; q[offset(0, 2, 4)] =  2.0; q[offset(0, 3, 4)] =  6.5;
	q[offset(1, 0, 4)] = 3.0; q[offset(1, 1, 4)] = -1.0; q[offset(1, 2, 4)] = -1.0; q[offset(1, 3, 4)] = -1.0;
	q[offset(2, 0, 4)] = 1.0; q[offset(2, 1, 4)] =  1.0; q[offset(2, 2, 4)] = -1.0; q[offset(2, 3, 4)] =  0.0;
	/* a = 1.0; b = 1.5; c = 2.5 */
#elif GE_SAMPLE == 3
	q[offset(0, 0, 5)] =  1.0; q[offset(0, 1, 5)] = -1.0; q[offset(0, 2, 5)] =  1.0; q[offset(0, 3, 5)] =  1.0; q[offset(0, 4, 5)] = 0.0;
	q[offset(1, 0, 5)] =  1.0; q[offset(1, 1, 5)] =  1.0; q[offset(1, 2, 5)] = -1.0; q[offset(1, 3, 5)] =  1.0; q[offset(1, 4, 5)] = 1.0;
	q[offset(2, 0, 5)] = -1.0; q[offset(2, 1, 5)] =  1.0; q[offset(2, 2, 5)] =  1.0; q[offset(2, 3, 5)] = -1.0; q[offset(2, 4, 5)] = 0.0;
	q[offset(3, 0, 5)] =  2.0; q[offset(3, 1, 5)] = -1.0; q[offset(3, 2, 5)] = -1.0; q[offset(3, 3, 5)] =  3.0; q[offset(3, 4, 5)] = 1.0;
	/* a = 0.0; b = 0.5; c = 0.0; d = 0.5 */
	dim = 4; /* 4 variables with 4 equations */
#elif GE_SAMPLE == 4 /* GE_RET_UNDETERMINED */
	q[offset(0, 0, 4)] = 3.0; q[offset(0, 1, 4)] = -1.0; q[offset(0, 2, 4)] =  1.0; q[offset(0, 3, 4)] = 8.0;
	q[offset(1, 0, 4)] = 1.0; q[offset(1, 1, 4)] =  2.0; q[offset(1, 2, 4)] = -1.0; q[offset(1, 3, 4)] = 4.0;
	q[offset(2, 0, 4)] = 2.0; q[offset(2, 1, 4)] = -3.0; q[offset(2, 2, 4)] =  2.0; q[offset(2, 3, 4)] = 4.0;
#elif GE_SAMPLE == 5 /* GE_RET_IMPOSSIBLE */
	q[offset(0, 0, 4)] = 2.0; q[offset(0, 1, 4)] =  1.0; q[offset(0, 2, 4)] = -1.0; q[offset(0, 3, 4)] = 4.0;
	q[offset(1, 0, 4)] = 1.0; q[offset(1, 1, 4)] = -1.0; q[offset(1, 2, 4)] =  1.0; q[offset(1, 3, 4)] = 2.0;
	q[offset(2, 0, 4)] = 1.0; q[offset(2, 1, 4)] =  2.0; q[offset(2, 2, 4)] = -2.0; q[offset(2, 3, 4)] = 1.0;
#endif

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

#ifdef GE_DYNAMIC
	free(q);
	free(result);
#endif

	return(0);
}
