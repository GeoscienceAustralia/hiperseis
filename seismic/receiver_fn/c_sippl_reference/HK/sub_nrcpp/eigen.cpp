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

#ifndef ROTATE
#define ROTATE(a, i, j, k, l) g = a[i][j]; h = a[k][l]; a[i][j] = g - s * (h + g * tau); a[k][l] = h + s * (g - h * tau);
#endif

/*----------------------------------------------------------------------------*/
/* Compute all eigenvalues and eigenvectors of a real symmetric matrix 
	a[1...n][1...n].
	On output, elements of 'a' above the diagonal are destroyed.
	d[1...n] returns the eigenvalues of 'a'.
	v[1...n][1...n] is a matrix whose columns contain the normalized 
	eigenvectors of 'a'.
	nrot returns the number of Jacobi rotations that were required.            */
/*----------------------------------------------------------------------------*/
void jacobi(float **a, int n, float d[], float **v, int *nrot)
{
	int j, iq, ip, i;
	float tresh, theta, tau, t, sm, s, h, g, c;
	float *b, *z;

	b = vector(1, n);
	z = vector(1, n);
	
/*----------------------------------------------------------------------------*/
/* Initialize to the identity matrix                                          */
/*----------------------------------------------------------------------------*/
	for (ip = 1; ip <= n; ip++) {
		for (iq = 1; iq <= n; iq++) {
			v[ip][iq] = 0.0;
		}	
		v[ip][ip] = 1.0;
	}
		
/*----------------------------------------------------------------------------*/
/* Initialize b and d to the diagonal of a
	This vector will accumulate terms of the form ta_pq as in 
	equation (11.1.14)                                                         */
/*----------------------------------------------------------------------------*/
	for (ip = 1; ip <= n; ip++) {
		b[ip] = d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}

	*nrot = 0;
	for (i = 1; i <= 50; i++) {
		sm = 0.0;
/*----------------------------------------------------------------------------*/
/* Sum off-diagonal elements                                                  */
/*----------------------------------------------------------------------------*/
		for (ip = 1; ip <= n - 1; ip++) {
			for (iq = ip + 1; iq <= n; iq++) {
				sm += fabs(a[ip][iq]);
			}
		}		
/*----------------------------------------------------------------------------*/
/* The normal return, which relies on quadratic convergence to 
	machine underflow                                                          */
/*----------------------------------------------------------------------------*/
		if (sm == 0.0) {
			free_vector(z, 1, n);
			free_vector(b, 1, n);
			z = NULL;
			b = NULL;
			return;
		}				
/*----------------------------------------------------------------------------*/
/* ... On the first three sweeps                                              */
/*----------------------------------------------------------------------------*/
		if (i < 4) {
			tresh = 0.2 * sm / (n * n);
		}
/*----------------------------------------------------------------------------*/
/* ... thereafter                                                             */
/*----------------------------------------------------------------------------*/
		else {
			tresh = 0.0;
		}			
		
		for (ip = 1; ip <= n - 1; ip++) {
			for (iq = ip + 1; iq <= n; iq++) {
				g = 100.0 * fabs(a[ip][iq]);
/*----------------------------------------------------------------------------*/
/* After 4 sweeps, skip the rotation if the off-diagonal element is small.    */
/*----------------------------------------------------------------------------*/
				if ( i > 4 && (float)(fabs(d[ip]) + g) == (float)fabs(d[ip])
						&& (float)(fabs(d[iq]) + g) == (float)fabs(d[iq]) ) {
					a[ip][iq] = 0.0;
				}
				else if (fabs(a[ip][iq]) > tresh) {
					h = d[iq] - d[ip];
					if ( (float)(fabs(h) + g) == (float)fabs(h) ) {
						t = a[ip][iq] / h;
					}
					else {
						theta = 0.5 * h / a[ip][iq];
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
						if (theta < 0.0) t = -t;
					}
					c = 1.0 / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
/*----------------------------------------------------------------------------*/
/* Case of rotations 1 <= j < p                                               */
/*----------------------------------------------------------------------------*/
					for (j = 1; j <= ip - 1; j++) {
						ROTATE(a, j, ip, j, iq)
					}						
/*----------------------------------------------------------------------------*/
/* Case of rotations p < j < q                                                */
/*----------------------------------------------------------------------------*/
					for (j = ip+1; j <= iq - 1; j++) {
						ROTATE(a, ip, j, j, iq)
					}						
/*----------------------------------------------------------------------------*/
/* Case of rotations q < j <= n                                               */
/*----------------------------------------------------------------------------*/
					for (j = iq+1; j <= n; j++) {
						ROTATE(a, ip, j, iq, j)
					}						
					for (j = 1; j <= n; j++) {
						ROTATE(v, j, ip, j, iq)
					}						
					++(*nrot);
				}
			}
		}			
			
/*----------------------------------------------------------------------------*/
/* Update d with the sum of ta_pq, and reinitialize z                         */
/*----------------------------------------------------------------------------*/
		for (ip = 1; ip <= n; ip++) {
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}
	nrerror("Too many iterations in routine jacobi");
}


/*----------------------------------------------------------------------------*/
/* Given the eigenvalues d[1...n] and eigenvectors v[1...n][1...n] as output 
	from jacobi, this routine sorts the eigenvalues into descending order, 
	and rearranges the columns of v corresponding. 
	The method is straight insertion.                                          */
/*----------------------------------------------------------------------------*/
void eigsrt(float d[], float **v, int n)
{
	int k, j, i;
	float p;
		
	for (i = 1; i <= n; i++) {  
		p = d[k = i];
		for (j = i + 1; j <= n; j++) {
			if (d[j] >= p) p = d[k = j];
		}	
		if (k != i) { 
			d[k] = d[i];
			d[i] = p;
			for (j = 1; j <= n; j++) {
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}			
				

				         
