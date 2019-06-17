#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "matrix.h"
#include "nrutil.h"

#define NR_END 1
#define FREE_ARG char*


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Allocate a float vector with subscript range v[nl,...,nh]                  */
/*----------------------------------------------------------------------------*/
float *vector(int nl, int nh)
{
	float *v = NULL;
	
	v = new float[nh-nl+1+NR_END]; 
	if (!v) nrerror("allocation failure in float vector()");
	return v-nl+NR_END;

}		

/*----------------------------------------------------------------------------*/
/* Allocate an int vector with subscript range v[nl,...,nh]                   */
/*----------------------------------------------------------------------------*/
int *ivector(int nl, int nh)
{
	int *v = NULL;
	
	v = new int[nh-nl+1+NR_END];
	if (!v) nrerror("allocation failure in int ivector()");
	return v-nl+NR_END;

}		

/*----------------------------------------------------------------------------*/
/* Allocate an unsigned char vector with subscript range v[nl,...,nh]         */
/*----------------------------------------------------------------------------*/
unsigned char *cvector(int nl, int nh)
{
	unsigned char *v = NULL;
	
	v = new unsigned char[nh-nl+1+NR_END];
	if (!v) nrerror("allocation failure in unsigned char cvector()");
	return v-nl+NR_END;

}		

/*----------------------------------------------------------------------------*/
/* Allocate an unsigned long vector with subscript range v[nl,...,nh]         */
/*----------------------------------------------------------------------------*/
unsigned long *lvector(int nl, int nh)
{
	unsigned long *v = NULL;
	
	v = new unsigned long[nh-nl+1+NR_END];
	if (!v) nrerror("allocation failure in unsigned long vector()");
	return v-nl+NR_END;

}		

/*----------------------------------------------------------------------------*/
/* Allocate a double vector with subscript range v[nl,...,nh]                 */
/*----------------------------------------------------------------------------*/
double *dvector(int nl, int nh)
{
	double *v = NULL;
	
	v = new double[nh-nl+1+NR_END];
	if (!v) nrerror("allocation failure in double dvector()");
	return v-nl+NR_END;

}		

/*----------------------------------------------------------------------------*/
/* Allocate a float complex vector with subscript range v[nl,...,nh]          */
/*----------------------------------------------------------------------------*/
fcomplex *fcvector(int nl, int nh)
{
	fcomplex *v = NULL;
	
	v = new fcomplex[nh-nl+1+NR_END];
	if (!v) nrerror("allocation failure in float complex cvector()");
	return v-nl+NR_END;

}		

/*----------------------------------------------------------------------------*/
/* Allocate a double complex vector with subscript range v[nl,...,nh]         */
/*----------------------------------------------------------------------------*/
dcomplex *dcvector(int nl, int nh)
{
	dcomplex *v = NULL;
	
	v = new dcomplex[nh-nl+1+NR_END];
	if (!v) nrerror("allocation failure in float complex cvector()");
	return v-nl+NR_END;

}		


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Allocate a float matrix with subscript range m[nrl...nrh][ncl...nch]       */
/*----------------------------------------------------------------------------*/
float **matrix(int nrl, int nrh, int ncl, int nch)
{
	int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	float **m;
	
/*----------------------------------------------------------------------------*/
/* Alocate pointers to rows                                                   */
/*----------------------------------------------------------------------------*/
	m = new float*[nrow+NR_END]; 
	if (!m) nrerror("allocation failure 1 in float matrix()");
	m += NR_END;
	m -= nrl;
	
/*----------------------------------------------------------------------------*/
/* Allocate rows and set pointers to them                                     */
/*----------------------------------------------------------------------------*/
	m[nrl] = new float[nrow*ncol+NR_END];
	if (!m[nrl]) nrerror("allocation failure 2 in float matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i-1] + ncol;
	
/*----------------------------------------------------------------------------*/
/* Return pointer to array of pointers to rows                                */
/*----------------------------------------------------------------------------*/
	return m;

}		

/*----------------------------------------------------------------------------*/
/* Allocate a double matrix with subscript range m[nrl...nrh][ncl...nch]      */
/*----------------------------------------------------------------------------*/
double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	double **m;
	
/*----------------------------------------------------------------------------*/
/* Alocate pointers to rows                                                   */
/*----------------------------------------------------------------------------*/
	m = new double*[nrow+NR_END];
	if (!m) nrerror("allocation failure 1 in double dmatrix()");
	m += NR_END;
	m -= nrl;
	
/*----------------------------------------------------------------------------*/
/* Allocate rows and set pointers to them                                     */
/*----------------------------------------------------------------------------*/
	m[nrl] = new double[nrow*ncol+NR_END];
	if (!m[nrl]) nrerror("allocation failure 2 in double dmatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i-1] + ncol;
	
/*----------------------------------------------------------------------------*/
/* Return pointer to array of pointers to rows                                */
/*----------------------------------------------------------------------------*/
	return m;

}		

/*----------------------------------------------------------------------------*/
/* Allocate an int matrix with subscript range m[nrl...nrh][ncl...nch]        */
/*----------------------------------------------------------------------------*/
int **imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	int **m;
	
/*----------------------------------------------------------------------------*/
/* Alocate pointers to rows                                                   */
/*----------------------------------------------------------------------------*/
	m = new int*[nrow+NR_END];
	if (!m) nrerror("allocation failure 1 in int imatrix()");
	m += NR_END;
	m -= nrl;
	
/*----------------------------------------------------------------------------*/
/* Allocate rows and set pointers to them                                     */
/*----------------------------------------------------------------------------*/
	m[nrl] = new int[nrow*ncol+NR_END];
	if (!m[nrl]) nrerror("allocation failure 2 in int imatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i-1] + ncol;
	
/*----------------------------------------------------------------------------*/
/* Return pointer to array of pointers to rows                                */
/*----------------------------------------------------------------------------*/
	return m;

}		

/*----------------------------------------------------------------------------*/
/* allocate a fcomplex matrix with subscript range m[nrl...nrh][ncl...nch]    */
/*----------------------------------------------------------------------------*/
fcomplex **fcmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	fcomplex **m;
	
/*----------------------------------------------------------------------------*/
/* Alocate pointers to rows                                                   */
/*----------------------------------------------------------------------------*/
	m = new fcomplex*[nrow+NR_END];
	if (!m) nrerror("allocation failure 1 in fcomplex matrix()");
	m += NR_END;
	m -= nrl;
	
/*----------------------------------------------------------------------------*/
/* Allocate rows and set pointers to them                                     */
/*----------------------------------------------------------------------------*/
	m[nrl] = new fcomplex[nrow*ncol+NR_END];
	if (!m[nrl]) nrerror("allocation failure 2 in fcomplex matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i-1] + ncol;
	
/*----------------------------------------------------------------------------*/
/* Return pointer to array of pointers to rows                                */
/*----------------------------------------------------------------------------*/
	return m;

}		

/*----------------------------------------------------------------------------*/
/* Allocate a dcomplex matrix with subscript range m[nrl...nrh][ncl...nch]    */
/*----------------------------------------------------------------------------*/
dcomplex **dcmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	dcomplex **m;
	
/*----------------------------------------------------------------------------*/
/* Alocate pointers to rows                                                   */
/*----------------------------------------------------------------------------*/
	m = new dcomplex*[nrow+NR_END];
	if (!m) nrerror("allocation failure 1 in dcomplex matrix()");
	m += NR_END;
	m -= nrl;
	
/*----------------------------------------------------------------------------*/
/* Allocate rows and set pointers to them                                     */
/*----------------------------------------------------------------------------*/
	m[nrl] = new dcomplex[nrow*ncol+NR_END];
	if (!m[nrl]) nrerror("allocation failure 2 in dcomplex matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	
	for (i = nrl + 1; i <= nrh; i++) m[i] = m[i-1] + ncol;
	
/*----------------------------------------------------------------------------*/
/* Return pointer to array of pointers to rows                                */
/*----------------------------------------------------------------------------*/
	return m;

}		


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Point a submatrix [newrl...][newcl...] to a[oldrl...oldrh][oldcl...oldch]. 
	We do not need upper row and column indices of the new matrix, 
	since they are implied by the quantities already given.                    */
/*----------------------------------------------------------------------------*/
float **submatrix(float **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl)
{
	int i, j, nrow=oldrh-oldrl+1, ncol=oldcl-newcl;
	float **m;
	
/*----------------------------------------------------------------------------*/
/* Alocate pointers to rows                                                   */
/*----------------------------------------------------------------------------*/
	m = new float*[nrow+NR_END];
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;
	
	for (i = oldrl, j = newrl; i <= oldrh; i++, j++) m[j] = a[i] + ncol;
	
/*----------------------------------------------------------------------------*/
/* Return pointer to array of pointers to rows                                */
/*----------------------------------------------------------------------------*/
	return m;

}		


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Allocate a float amatrix m[nrl...nrh][ncl...nch]] that point to the matrix
	declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1 
	and ncol=nch-ncl+1. 
	The routine should be called with the address &a[0][0] as the first argument.
	i.e., 'a' is zero-shift.                                                   */
/*----------------------------------------------------------------------------*/
float **convert_matrix(float *a, int nrl, int nrh, int ncl, int nch)
{
	int i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1;
	float **m;
	
/*----------------------------------------------------------------------------*/
/* Alocate pointers to rows                                                   */
/*----------------------------------------------------------------------------*/
	m = new float*[nrow+NR_END];
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;
	
	m[nrl] = a - ncl;
	for (i = 1, j = nrl + 1; i <= nrow; i++, j++) m[j] = m[j-1] + ncol;
	
/*----------------------------------------------------------------------------*/
/* Return pointer to array of pointers to rows                                */
/*----------------------------------------------------------------------------*/
	return m;

}		


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Allocate a float 3D tensor with range t[nrl...nrh][ncl...nch][ndl...ndh]   */
/*----------------------------------------------------------------------------*/
float ***f3tensor(int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
	int i, j, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ndh-ndl+1;
	float ***t;
	
/*----------------------------------------------------------------------------*/
/* Alocate pointers to rows                                                   */
/*----------------------------------------------------------------------------*/
	t = new float**[nrow+NR_END];
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;
	
/*----------------------------------------------------------------------------*/
/* Allocate pointers to rows and set pointers to them                         */
/*----------------------------------------------------------------------------*/
	t[nrl] = new float*[nrow*ncol+NR_END];
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	
/*----------------------------------------------------------------------------*/
/* Allocate rows and set pointers to them                                     */
/*----------------------------------------------------------------------------*/
	t[nrl][ncl] = new float[nrow*ncol*ndep+NR_END];
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	
	for (j = ncl + 1; j <= nch; j++) t[nrl][j] = t[nrl][j-1] + ndep;
	for (i = nrl + 1; i <= nrh; i++) {
		t[i] = t[i-1] + ncol;
		t[i][ncl] = t[i-1][ncl] + ncol * ndep;
		for (j = ncl + 1; j < nch; j++) t[i][j] = t[i][j-1] + ndep;
	}
			
/*----------------------------------------------------------------------------*/
/* Return pointer to array of pointers to rows                                */
/*----------------------------------------------------------------------------*/
	return t;

}		


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Free a float vector allocated with vector()                                */
/*----------------------------------------------------------------------------*/
void free_vector(float *v, int nl, int nh)
{
	float *vtmp;
	
	vtmp = v+nl-NR_END;
	delete [] vtmp;
	vtmp = NULL;
	v = NULL;
}

/*----------------------------------------------------------------------------*/
/* Free an int vector allocated with ivector()                                */
/*----------------------------------------------------------------------------*/
void free_ivector(int *v, int nl, int nh)
{
	int *vtmp;
	
	vtmp = v+nl-NR_END;
	delete [] vtmp;
	vtmp = NULL;
	v = NULL;
}

/*----------------------------------------------------------------------------*/
/* Free an unsinged char vector allocated with cvector()                      */
/*----------------------------------------------------------------------------*/
void free_cvector(unsigned char *v, int nl, int nh)
{
	unsigned char *vtmp;
	
	vtmp = v+nl-NR_END;
	delete [] vtmp;
	vtmp = NULL;
	v = NULL;
}

/*----------------------------------------------------------------------------*/
/* Free an unsinged long vector allocated with lvector()                      */
/*----------------------------------------------------------------------------*/
void free_lvector(unsigned long *v, int nl, int nh)
{
	unsigned long *vtmp;
	
	vtmp = v+nl-NR_END;
	delete [] vtmp;
	vtmp = NULL;
	v = NULL;
}

/*----------------------------------------------------------------------------*/
/* Free a double vector allocated with dvector()                              */
/*----------------------------------------------------------------------------*/
void free_dvector(double *v, int nl, int nh)
{
	double *vtmp;
	
	vtmp = v+nl-NR_END;
	delete [] vtmp;
	vtmp = NULL;
	v = NULL;
}

/*----------------------------------------------------------------------------*/
/* Free a float complex vector allocated with fcvector()                      */
/*----------------------------------------------------------------------------*/
void free_fcvector(fcomplex *v, int nl, int nh)
{
	fcomplex *vtmp;
	
	vtmp = v+nl-NR_END;
	delete [] vtmp;
	vtmp = NULL;
	v = NULL;
}

/*----------------------------------------------------------------------------*/
/* Free a double complex vector allocated with dcvector()                     */
/*----------------------------------------------------------------------------*/
void free_dcvector(dcomplex *v, int nl, int nh)
{
	dcomplex *vtmp;
	
	vtmp = v+nl-NR_END;
	delete [] vtmp;
	vtmp = NULL;
	v = NULL;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Free a float matrix allocated with matrix()                                */
/*----------------------------------------------------------------------------*/
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
{
	float *mrow;
	float **mtmp;
	
	mrow = m[nrl]+ncl-NR_END;
	delete [] mrow;
	mrow = NULL;
	m[nrl] = NULL;
	mtmp = m+nrl-NR_END;
	delete [] mtmp;
	mtmp = NULL;
	m = NULL;
}

/*----------------------------------------------------------------------------*/
/* Free a double matrix allocated with dmatrix()                              */
/*----------------------------------------------------------------------------*/
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	double *mrow;
	double **mtmp;

	mrow = m[nrl]+ncl-NR_END;
	delete [] mrow;
	mrow = NULL;
	m[nrl] = NULL;
	mtmp = m+nrl-NR_END;
	delete [] mtmp;
	mtmp = NULL;
	m = NULL;
}

/*----------------------------------------------------------------------------*/
/* Free an int matrix allocated with imatrix()                                */
/*----------------------------------------------------------------------------*/
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int *mrow;
	int **mtmp;

	mrow = m[nrl]+ncl-NR_END;
	delete [] mrow;
	mrow = NULL;
	m[nrl] = NULL;
	mtmp = m+nrl-NR_END;
	delete [] mtmp;
	mtmp = NULL;
	m = NULL;
}

/*----------------------------------------------------------------------------*/
/* Free a fcomplex matrix allocated with dmatrix()                            */
/*----------------------------------------------------------------------------*/
void free_fcmatrix(fcomplex **m, int nrl, int nrh, int ncl, int nch)
{
	fcomplex *mrow;
	fcomplex **mtmp;

	mrow = m[nrl]+ncl-NR_END;
	delete [] mrow;
	mrow = NULL;
	m[nrl] = NULL;
	mtmp = m+nrl-NR_END;
	delete [] mtmp;
	mtmp = NULL;
	m = NULL;
}

/*----------------------------------------------------------------------------*/
/* Free a dcomplex matrix allocated with dmatrix()                            */
/*----------------------------------------------------------------------------*/
void free_dcmatrix(dcomplex **m, int nrl, int nrh, int ncl, int nch)
{
	dcomplex *mrow;
	dcomplex **mtmp;

	mrow = m[nrl]+ncl-NR_END;
	delete [] mrow;
	mrow = NULL;
	m[nrl] = NULL;
	mtmp = m+nrl-NR_END;
	delete [] mtmp;
	mtmp = NULL;
	m = NULL;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Free a float submatrix allocated with submatrix()                          */
/*----------------------------------------------------------------------------*/
void free_submatrix(float **b, int nrl, int nrh, int ncl, int nch)
{
	float **btmp;
	
	btmp = b+nrl-NR_END;
	delete [] btmp;
	btmp = NULL;
	b = NULL;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Free a float submatrix allocated with convert_matrix()                     */
/*----------------------------------------------------------------------------*/
void free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch)
{
	float **btmp;

	btmp = b+nrl-NR_END;
	delete [] btmp;
	btmp = NULL;
	b = NULL;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Free a float 3D vector allocated with f3tensor()                           */
/*----------------------------------------------------------------------------*/
void free_f3tensor(float ***t, int nrl, int nrh, int ncl, int nch, int ndl, int ndh)
{
	float *trowcol;
	float **trow;
	float ***ttmp;
	
	trowcol = t[nrl][ncl]+ndl-NR_END;
	delete [] trowcol;
	trowcol = NULL;
	t[nrl][ncl] = NULL;
	trow = t[nrl]+ncl-NR_END;
	delete [] trow;
	trow = NULL;
	t[nrl] = NULL;
	ttmp = t+nrl-NR_END;
	delete [] ttmp;
	ttmp = NULL;
	t = NULL;
}


