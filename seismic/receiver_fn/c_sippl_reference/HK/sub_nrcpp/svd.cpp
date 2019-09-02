#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cstdio>
#include <cmath>
#include "nrutil.h"
#include "matrix.h"

/*----------------------------------------------------------------------------*/
/* Define a small number                                                      */
/*----------------------------------------------------------------------------*/
#ifndef TINY
#define TINY 1.0e-20;
#endif

#ifndef MAX_SVD_ITERS
#define MAX_SVD_ITERS 50
#endif


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Given a matrix a[1..n][1..n] , this routine replaces it by 
	the LU decomposition of a rowwise permutation of itself. a and n are input. 
	a = |a11, a12, a13, a14|
	    |a21, a22, a23, a24|
	    |a31, a32, a33, a34|
	    |a41, a42, a43, a44| for an example
	a is output, arranged as 
	a = |u11, u12, u13, u14|
	    |l21, u22, u23, u24|
	    |l31, l32, u33, u34|
	    |l41, l42, l43, u44| 
	Note that all l11, l22, l33, l44 = 1, so do not need to store.
	The original a is destroyed.
	indx[1..n] is an output vector that records the row permutation effected by
	the partial pivoting.
	d is return as +/- 1 depending on whether the number of row interchanges 
	was even or odd, respectively. 
	This routine is used in combination with lubksb to solve linear equations 
	or invert a square matrix.                                                 */
/*----------------------------------------------------------------------------*/
int ludcmp(float **a, int n, int *indx)
{
	int i, j, k;
	int imax = 0;
	int d;
	float big;
	float dum;
	float sum;
	float temp;
/*----------------------------------------------------------------------------*/
/* vv stores the implicit scaling of each row.                                */
/*----------------------------------------------------------------------------*/
	float *vv = NULL; 

	vv = vector(1, n);

/*----------------------------------------------------------------------------*/
/* Initial value of return value d 
	No row interchanges yet.                                                   */
/*----------------------------------------------------------------------------*/
	d = 1;

/*----------------------------------------------------------------------------*/
/* Loop over rows to get the implicit scaling information.                    */
/*----------------------------------------------------------------------------*/
	for (i = 1; i <= n; i++) {

		big = 0.0;	
		for (j = 1; j <= n; j++) {
			if ( (temp = fabs(a[i][j])) > big ) big = temp;
		}
		
/*----------------------------------------------------------------------------*/
/* No nonzero largest element.                                                */
/*----------------------------------------------------------------------------*/
		if (big == 0.0) {
			nrerror("Singular matrix in routine ludcmp"); 
			d = 0;
		}
			
/*----------------------------------------------------------------------------*/
/* Save the scaling.                                                          */
/*----------------------------------------------------------------------------*/
		vv[i] = 1.0 / big;
	}

/*----------------------------------------------------------------------------*/
/* This is the loop over columns of Crout¡¯s method.                           */
/*----------------------------------------------------------------------------*/
	for (j = 1; j <= n; j++) { 
		
/*----------------------------------------------------------------------------*/
/* This is 
                       i-1   
	u[i][j] = a[i][j] - sum(l[i][k] * u[k][j] except for i = j.         
	                    k=l                                      
	in Crout's method                                                          */
/*----------------------------------------------------------------------------*/
		for (i = 1; i < j; i++) {
			sum = a[i][j];	
			for (k = 1; k < i; k++) {
				sum -= a[i][k] * a[k][j];
			}
			a[i][j] = sum;
		}

/*----------------------------------------------------------------------------*/
/* Initialize for the search for largest pivot element.                       */
/*----------------------------------------------------------------------------*/
		big = 0.0;
/*----------------------------------------------------------------------------*/
/* This is i = j of above equation, and
                                    j-1   
	l[i][j] = 1/u[j][i] * (a[i][j] - sum(l[i][k] * u[k][j] 
	                                 k=1  
	for i = j+1 . . . N                                                        */
/*----------------------------------------------------------------------------*/
		for (i = j; i <= n; i++) {
			sum = a[i][j];
			for (k = 1; k<j; k++) {
				sum -= a[i][k] * a[k][j];
			}
			a[i][j] = sum;

/*----------------------------------------------------------------------------*/
/* Is the figure of merit for the pivot better than the best so far?          */
/*----------------------------------------------------------------------------*/
			if ( (dum = vv[i] * fabs(sum)) >= big ) {
				big = dum; 
				imax = i;
			}
		}

/*----------------------------------------------------------------------------*/
/* Do we need to interchange rows? Yes, do so...                              */
/*----------------------------------------------------------------------------*/
		if (j != imax) {
			for (k = 1; k <= n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k]; 
				a[j][k] = dum;
			}

/*----------------------------------------------------------------------------*/
/* ...and change the parity of d.                                            */
/*----------------------------------------------------------------------------*/
			d = -d;
		
/*----------------------------------------------------------------------------*/
/* Also interchange the scale factor.                                         */
/*----------------------------------------------------------------------------*/
			vv[imax] = vv[j];
		}

		indx[j] = imax;
/*----------------------------------------------------------------------------*/
/* If the pivot element is zero, the matrix is singular 
	(at least to the precision of the algorithm). 
	For some applications on singular matrices, it is desirable to 
	substitute TINY for zero.                                                  */
/*----------------------------------------------------------------------------*/
		if (a[j][j] == 0.0) a[j][j] = TINY;

/*----------------------------------------------------------------------------*/
/* Now, finally, divide by the pivot element.                                 */
/*----------------------------------------------------------------------------*/
		if (j != n) {
			dum = 1.0 / (a[j][j]);
			for (i = j+1; i <= n; i++) {
				a[i][j] *= dum;
			}
		}

/*----------------------------------------------------------------------------*/
/* Go back for the next column in the reduction.                              */
/*----------------------------------------------------------------------------*/
	}
	
	free_vector(vv, 1, n);
	vv = NULL;
	return d;	
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Solves the set of n linear equations A ¡¤ X = B using LU solution. 
	Here a[1..n][1..n] is input, not as the matrix A but rather as 
	its LU decomposition, determined by the routine ludcmp. 
	indx[1..n] is input as the permutation vector returned by ludcmp. 
	b[1..n] is input as the right-hand side vector B, and returns with 
	the solution vector X. 
	Note the original b is destroyed.
	a, n, and indx are not modified by this routine and can be left 
	in place for successive calls with different right-hand sides b. 
	This routine takes into account the possibility that b will begin with 
	many zero elements, so it is efficient for use in matrix inversion.        */
/*----------------------------------------------------------------------------*/
void lubksb(float **a, int n, int *indx, float b[])
{
	int i, j;
	int ii = 0;
	int ip;
	float sum;

/*----------------------------------------------------------------------------*/
/* When ii is set to a positive value, it will become the index of the first 
	nonvanishing element of b. 
	We now do the forward substitution, 
	y[1] = b[1]/l[1][1]
	                           i-1   
	y[i] = 1/l[i][i] * (b[i] - sum(l[i][j] * y[j])
	                           j=1    
	for i = 2, 3, ... N                           
	The only new wrinkle is to unscramble the permutation as we go.            */
/*----------------------------------------------------------------------------*/
	for (i = 1; i <= n; i++) {
		ip = indx[i]; 
		sum = b[ip]; 
		b[ip] = b[i]; 
		
		if (ii) {
			for (j = ii; j <= i-1; j++) {
				sum -= a[i][j] * b[j];
			}
		}	
/*----------------------------------------------------------------------------*/
/* A nonzero element was encountered, so from now on we will have to do 
	the sums in the loop above.                                                */
/*----------------------------------------------------------------------------*/
		else if (sum) {
			ii=i; 
		}
		
		b[i] = sum;
	} 

/*----------------------------------------------------------------------------*/
/* Now we do the backsubstitution, 
	x[N] = y[N]/u[N][N]
	                           N   
	x[i] = 1/u[i][i] * (y[i] - sum(u[i][j] * x[j])
	                           j=i+1    
	for i = N-1, N-2, N-3, ... 1                                               */
/*----------------------------------------------------------------------------*/
	for (i = n; i >= 1; i--) {
		sum = b[i];

		for (j = i+1; j <= n; j++) {
			sum -= a[i][j] * b[j];
		}

/*----------------------------------------------------------------------------*/
/* Store a component of the solution vector X . All done!.                    */
/*----------------------------------------------------------------------------*/
		b[i] = sum / a[i][i];
	}
}	

	 
/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Using the above LU decomposition and lubksb() routines, 
	it is completely straightforward to find the inverse of a square matrix 
	column by column.
	a[1..n][1..n] is input matrix, returned as inverse matrix                  */	
/*----------------------------------------------------------------------------*/
int invmatrix(float **a, int n)
{
	int i, j;
	int *indx = NULL;	
	int d;
	float *col = NULL;
	float **y = NULL;
	
	int ludcmp(float **a, int n, int *indx);
	void lubksb(float **a, int n, int *indx, float b[]);
	
	indx = ivector(1, n);
	col = vector(1, n);
	y = matrix(1, n, 1, n); 
	
/*----------------------------------------------------------------------------*/
/* Decompose the matrix just once.                                            */
/*----------------------------------------------------------------------------*/
	if (ludcmp(a, n, indx) == 0) {
		nrerror("Warning in invmatrix: Error in ludcmp, the result is unstable"); 
		d = -1;
	}	 
	else {
		d = 0;
	}
		
/*----------------------------------------------------------------------------*/
/* Find inverse by columns.                                                   */ 
/*----------------------------------------------------------------------------*/
	for (j = 1; j <= n; j++) {

		for (i = 1; i <= n; i++) {
			col[i]=0.0; 
		}
		col[j] = 1.0;
		lubksb(a, n, indx, col);

		for (i = 1; i <= n; i++) {
			y[i][j] = col[i];
		}
	}
	
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) {
			a[i][j] = y[i][j];
		}
	}
				
	free_ivector(indx, 1, n);
	indx = NULL;
	free_vector(col, 1, n);
	col = NULL;
	free_matrix(y, 1, n, 1, n);
	y = NULL;
	
	return d;	
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Improves a solution vector x[1..n] of the linear set of equations A¡¤X = B. 
	The matrix a[1..n][1..n], and the vectors b[1..n] and x[1..n] are input, 
	as is the dimension n. 
	Also input is alud[1..n][1..n], the LU decomposition of a[][] 
	as returned by ludcmp, the vector indx[1..n] also returned by that routine,
	and x[1...n] returned by lubksb. 
	On output, only x[1..n] is modified, to an improved set of values. 
	Note that since the original a[][] is destroyed by calling ludcmp, and 
	the original b[] is destroyed by calling lubksb, make copy of them before 
	using this subroutine.                                                     */
/*----------------------------------------------------------------------------*/
void mprove(float **a, float **alud, int n, int indx[], float b[], float x[])
{
	int i, j;	
	double sdp; 
	float *r = NULL;
	
	void lubksb(float **a, int n, int *indx, float b[]); 

	r = vector(1, n);

/*----------------------------------------------------------------------------*/
/* Calculate the right-hand side, accumulating the residual 
	in double precision.                                                       */
/*----------------------------------------------------------------------------*/
	for (i = 1; i <= n; i++) {
		sdp = -b[i];
		for (j = 1; j <= n; j++) {
			sdp += a[i][j] * x[j]; 
		}	
		r[i] = sdp;
	}
	
/*----------------------------------------------------------------------------*/
/* Solve for the error term, and subtract it from the old solution.           */
/*----------------------------------------------------------------------------*/
	lubksb(alud, n, indx, r);
	for (i = 1; i <= n; i++) {
		x[i] -= r[i];
	}

	free_vector(r, 1, n);
	r = NULL;	
}

		
/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Given a matrix a[1...m][1...n], this routine computes its singular value 
	decomposion, A = U * W * V^T. 
	The matrix U replace matrix a on output. 
	The diagonal matrix of singular values W is output as a vector W[1...n]. 
	The matrix V (not the transpose V^T) is output as V[1...n][1...n].         */
/*----------------------------------------------------------------------------*/
int svdcmp(float **a, int m, int n, float w[], float **v)
{
	float pythag(float a, float b);
	
	int flag, its, nm;
	int i, j, k, l, jj;
	float anorm, c, f, g, h, s, scale, x, y, z;
	float *rv1 = NULL;
	int return_value;

	return_value = 0;
	rv1 = vector(1, n);
/*----------------------------------------------------------------------------*/
/* Householder reduction to bidiagonal form                                   */	
/*----------------------------------------------------------------------------*/
	g = scale = anorm = 0.0;    
	
	for (i = 1; i <= n; i++) {
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (k = i; k <= m; k++)  scale += fabs(a[k][i]);
			if (scale) {
				for (k = i; k <= m; k++) {
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][i] = f - g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = i; k <= m; k++)  s += a[k][i] * a[k][j];
					f = s / h;
					for (k = i; k <= m; k++)  a[k][j] += f * a[k][i];
				}
				for (k = i; k <= m; k++)  a[k][i] *= scale;
			}
		}
		
		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m && i != n) {
			for (k = l; k <= n; k++)  scale += fabs(a[i][k]);
			if (scale) {
				for (k = l; k <= n; k++) {
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][l] = f - g;
				for (k = l; k <= n; k++)  rv1[k] = a[i][k] / h;
				for (j = l; j <= m; j++) {
					for (s = 0.0, k = l; k <= n; k++)  s += a[j][k] * a[i][k];
					for (k = l; k <= n; k++)  a[j][k] += s * rv1[k];
				}
				for (k = l; k <= n; k++)  a[i][k] *= scale;
			}
		}
		anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	
/*----------------------------------------------------------------------------*/
/* Accumulation of right-hand transformations                                 */
/*----------------------------------------------------------------------------*/
	for (i = n; i >= 1; i--) {            
		if (i < n) {
			if (g) {
				for (j = l; j <= n; j++)     
/*----------------------------------------------------------------------------*/
/* Double division to avoid possible underflow                                */
/*----------------------------------------------------------------------------*/
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= n; k++)  s += a[i][k] * v[k][j];
					for (k = l; k <= n; k++)  v[k][j] += s * v[k][i];
				}
			}
			for (j = l; j <= n; j++)  v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	
/*----------------------------------------------------------------------------*/
/* Accumulation of left-hand transformations                                  */
/*----------------------------------------------------------------------------*/
	for (i = IMIN(m, n); i >= 1; i--) {   
		l = i + 1;
		g = w[i];
		for (j = l; j <= n; j++)  a[i][j] = 0.0;
		if (g) {
			g = 1.0 / g;
			for (j = l; j <= n; j++) {
				for (s = 0.0, k = l; k <= m; k++)  s += a[k][i] * a[k][j];
				f = (s / a[i][i]) * g;
				for (k = i; k <= m; k++)  a[k][j] += f * a[k][i];
			}
			for (j = i; j <= m; j++)  a[j][i] *= g;
		}
		else {
			for (j = i; j <= m; j++)  a[j][i] = 0.0;
		}
		++a[i][i];
	}
	
/*----------------------------------------------------------------------------*/
 /* Diagonalization of the bidiagonal form:                                   */
/*----------------------------------------------------------------------------*/
	for (k = n; k >= 1; k--) {           
/*----------------------------------------------------------------------------*/
/* Loop over singular values, and over allowed iterations                     */
/*----------------------------------------------------------------------------*/
		for (its = 1; its <= MAX_SVD_ITERS; its++) {       										
			flag = 1;
/*----------------------------------------------------------------------------*/
/* Test for splitting                                                         */
/*----------------------------------------------------------------------------*/
			for (l = k; l >= 1; l--) {      
				nm = l - 1;                  
/*----------------------------------------------------------------------------*/
/* Note that rv1[1] is always zero                                            */
/*----------------------------------------------------------------------------*/
				if ((float)(fabs(rv1[l]) + anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((float)(fabs(w[nm]) + anorm) == anorm)  break;
			}
			if (flag) {
/*----------------------------------------------------------------------------*/
/* Cancellation of rv1[l], if l > 1                                           */
/*----------------------------------------------------------------------------*/
				c = 0.0;                     
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ((float)(fabs(f) + anorm) == anorm)  break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 1; j <= m; j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y * c + z * s;
						a[j][i] = z * c - y * s;
					}
				}
			}	
			
			z = w[k];
/*----------------------------------------------------------------------------*/
/* Convergence                                                                */
/*----------------------------------------------------------------------------*/
			if (l == k) {                   
/*----------------------------------------------------------------------------*/
/* Singular value is made nonnegative                                         */
/*----------------------------------------------------------------------------*/
				if (z < 0.0) {               
					w[k] = -z;
					for (j = 1; j <= n; j++)  v[j][k] = -v[j][k];
				}
				break;
			}
			
			if (its == MAX_SVD_ITERS) {
				nrerror("No convergence in MAX svdcmp iterations.");
				return_value = 1;
			}
/*----------------------------------------------------------------------------*/
/* Shift form bottom 2-by-2 minor                                             */
/*----------------------------------------------------------------------------*/
			x = w[l];                      
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;                   
/*----------------------------------------------------------------------------*/
/* Next QR transformation                                                     */
/*----------------------------------------------------------------------------*/
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 1; jj <= n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j] = z;                    
/*----------------------------------------------------------------------------*/
/* Rotation can be arbitrary if z = 0                                         */
/*----------------------------------------------------------------------------*/
				if (z) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 1; jj <= m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free_vector(rv1, 1, n);
	rv1 = NULL;						   	
	return(return_value);					 		   				 							
}

	
/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Computes sqrt(a^2 + b^2) without destructive underflow or overflow         */
/*----------------------------------------------------------------------------*/
float pythag(float a, float b)
{
	float absa, absb;
	
	absa = fabs(a);
	absb = fabs(b);
	
	if (absa > absb)  return absa * sqrt(1.0 + SQR(absb / absa));
	else  return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb)));
}	


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Soves A * X = B for a vector X, using svd solution. 
	where A is specified by the arrays u[1...m][1...n],
	W[1...n], v[1...n][1...n] as returned by svdcmp. 
	m and n are the dimensions of a, and will be equal for square matrices. 
	b[1...m] is the input right-hand side. 
	x[1...n] is the output solution vector.
	No input quantites are destroyed, so the routine may be called sequentially 
	with different b's.                                                        */
/*----------------------------------------------------------------------------*/
void svbksb(float **u, float w[], float **v, int m, int n, float b[], float x[])
{
	int jj, j, i;
	float s;
	float *tmp = NULL;
	
	tmp = vector(1, n);
	
/*----------------------------------------------------------------------------*/
/* Calculate U^T * B                                                          */
/*----------------------------------------------------------------------------*/
	for (j = 1; j <= n; j++) {                
		s = 0.0;
/*----------------------------------------------------------------------------*/
/* Nonzero result only if w[j] is nonzero                                     */
/*----------------------------------------------------------------------------*/
		if (w[j]) {                            
			for (i = 1; i <= m; i++)  s += u[i][j] * b[i];
/*----------------------------------------------------------------------------*/
/* This is the divide by w[j]                                                 */
/*----------------------------------------------------------------------------*/
			s /= w[j];                          
		}
		tmp[j] = s;
	}
	
/*----------------------------------------------------------------------------*/
/* Matrix multiply by V to get answer                                         */
/*----------------------------------------------------------------------------*/
	for (j = 1; j <= n; j++) {                
		s = 0.0;
		for (jj = 1; jj <= n; jj++) s += v[j][jj] * tmp[jj];
		x[j] = s;
	}
	
	free_vector(tmp, 1, n);
	tmp = NULL;
}

				         
