#ifndef _NR_H_
#define _NR_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cstdio>
#include <cmath>


#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned int *ilob, *iupb, *ncumfq, jdif, nc, minint, nch, ncum, nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif
	
#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned int *icod, *ncod, *left, *right, nch, nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif


/*----------------------------------------------------------------------------*/
/* In ffts.c                                                                  */
/*----------------------------------------------------------------------------*/
void four1(float data[], unsigned int nn, int isign);
void twofft(float data1[], float data2[], float fft1[], float fft2[], unsigned int n);
void realft(float data[], unsigned int n, int isign);
void realft2(float data[], unsigned int n, int isign);
void correl(float data1[], float data2[], unsigned int n, float ans[]);
void convlv(float data[], unsigned int n, float respns[], unsigned int m, 
				int isign, float ans[]);
				
/*----------------------------------------------------------------------------*/
/* In svd.c                                                                   */
/*----------------------------------------------------------------------------*/
int ludcmp(float **a, int n, int *indx);
void lubksb(float **a, int n, int *indx, float b[]);
int invmatrix(float **a, int n);
void mprove(float **a, float **alud, int n, int indx[], float b[], float x[]);
int svdcmp(float **a, int m, int n, float w[], float **v);
float pythag(float a, float b);
void svbksb(float **u, float w[], float **v, int m, int n, float b[], float x[]);
				
/*----------------------------------------------------------------------------*/
/* In eigen.c                                                                 */
/*----------------------------------------------------------------------------*/
void jacobi(float **a, int n, float d[], float **v, int *nrot);
void eigsrt(float d[], float **v, int n);

/*----------------------------------------------------------------------------*/
/* In sort.c                                                                  */
/*----------------------------------------------------------------------------*/
void fflip(float *a, int n);
void iflip(int *a, int n);
int bubble(float *dat, int n, int option);
int bubble_indx(float *dat, int *index, int n, int option);
int shell(float *a, int n, int option);
int quick(float *S, int *IDX, int p, int r, int option);
int quick_indx(float arr[], unsigned int indx[], unsigned int n, int option);
void quick_rank(unsigned int indx[], unsigned int irank[], unsigned int n);


#endif /* _NR_H_ */

 		
