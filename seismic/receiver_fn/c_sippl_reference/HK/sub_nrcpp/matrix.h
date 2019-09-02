#ifndef _NR_MATRIX_H_
#define _NR_MATRIX_H_

#include "complex.h"

float *vector(int nl, int nh);
int *ivector(int nl, int nh);
unsigned char *cvector(int nl, int nh);
unsigned long *lvector(int nl, int nh);
double *dvector(int nl, int nh);
fcomplex *fcvector(int nl, int nh);
dcomplex *dcvector(int nl, int nh);

float **matrix(int nrl, int nrh, int ncl, int nch);
double **dmatrix(int nrl, int nrh, int ncl, int nch);
int **imatrix(int nrl, int nrh, int ncl, int nch);
fcomplex **fcmatrix(int nrl, int nrh, int ncl, int nch);
dcomplex **dcmatrix(int nrl, int nrh, int ncl, int nch);
float **submatrix(float **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl);
float **convert_matrix(float *a, int nrl, int nrh, int ncl, int nch);
float ***f3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh);

void free_vector(float *v, int nl, int nh);
void free_ivector(int *v, int nl, int nh);
void free_cvector(unsigned char *v, int nl, int nh);
void free_lvector(unsigned long *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh);
void free_fcvector(fcomplex *v, int nl, int nh);
void free_dcvector(dcomplex *v, int nl, int nh);

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_fcmatrix(fcomplex **m, int nrl, int nrh, int ncl, int nch);
void free_dcmatrix(dcomplex **m, int nrl, int nrh, int ncl, int nch);
void free_submatrix(float **b, int nrl, int nrh, int ncl, int nch);
void free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch);
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh);


#endif /* _NR_MATRIX_H_ */
