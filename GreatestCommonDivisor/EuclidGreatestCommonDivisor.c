/* Andre Augusto Giannotti Scota
 * andre.scota@gmail.com
 * MIT License
 */

#include <stdio.h>
#include <stdlib.h>

/* GCD by Euclid Algorithm */
unsigned long gcd(unsigned long x, unsigned long y)
{
	unsigned long aux = 0;
 
	if((x <= 0) || (y <= 0))
		return(0);

	do{
		if(x < y){
			aux = x;
			x = y;
			y = aux;
		}

		x -= y;

	}while(x > 0);

	return(y);
}
 
int main(int argc, char **argv)
{
	if(argc != 3){
		printf("Usage: %s <NUM1> <NUM2>  (NUM limited by MAX_ULONG)\n", argv[0]);
		return(-1);
	}

	printf("Common Divisor: [%ld]\n", gcd(strtoul(argv[1], NULL, 10), strtoul(argv[2], NULL, 10)));

	return(0);
}
