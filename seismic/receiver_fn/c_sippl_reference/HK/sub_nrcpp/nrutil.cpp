#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil.h"

#define NR_END 1
#define FREE_ARG char*


/*----------------------------------------------------------------------------*/
/* Standard error handler in C                                                */
/*----------------------------------------------------------------------------*/
void nrerror(char error_text[])
{
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "%s\n", error_text);
	fprintf(stderr, "...now exiting to system...\n\n");
	exit(1);
}

/*----------------------------------------------------------------------------*/
/* Standard error handler in C++                                                */
/*----------------------------------------------------------------------------*/
void NRcatch(NRerror err) 
{
	fprintf(stderr, "Numerical Recipes run-time error...\n");
	fprintf(stderr, "ERROR: %s\n, in file %s at line %d\n", 
		err.message, err.file, err.line);
	fprintf(stderr, "...now exiting to system...\n\n");
	exit(1);
}