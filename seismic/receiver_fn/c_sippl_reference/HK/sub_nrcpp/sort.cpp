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
/* This file contains a number of sorting subroutines.
	Besides directly sorting the given array, they may also come with 
	index and/rank return. In that case, the original array are not changed.

	org array:      a[1]=14     a[2]=8      a[3]=32     a[4]=7      a[5]=3      a[6]=15

	index table: indx[1]=[5] indx[2]=[4] indx[3]=[2] indx[4]=[1] indx[5]=[6] indx[6]=[3]
	it means: indx[1] is the min, which is a[5], indx[2] is the 2nd min, which is a[4]...

	rank table:  rank[1]=[4] rank[2]=[3] rank[3]=[6] rank[4]=[2] rank[5]=[1] rank[6]=[5]
	it means: a[1] is on location [4], a[2] is on location [3]... after sorting

	sorted array:   a[1]=3      a[2]=7      a[3]=8       a[4]=14    a[5]=15     a[6]=32   
*/
/*----------------------------------------------------------------------------*/


#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp; 
#define M 7
#define NSTACK 100


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Reverse float array a[1...n].
   When return a[n] = a[1], a[n-1] = a[2] .... a[2] = a[n-1], a[1] = a[n] 
   If n is even, swap a[n/2] and a[n/2+1].
   If n is odd, the center one a[n/2+1] is not changed.  
   Note the original a[] is destroyed.                                        */
/*----------------------------------------------------------------------------*/
void fflip(float *a, int n)
{
	int i;
	int k;
	float tmp;
	
	k = n / 2;
	for (i = 1; i <= k; i++) {
		tmp = a[i];
		a[i] = a[n - i + 1];
		a[n - i + 1] = tmp;
	}
}		


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Reverse int array a[1...n].
   When return a[n] = a[1], a[n-1] = a[2] .... a[2] = a[n-1], a[1] = a[n] 
   If n is even, swap a[n/2] and a[n/2+1].
   If n is odd, the center one a[n/2+1] is not changed.  
   Note the original a[] is destroyed.                                        */
/*----------------------------------------------------------------------------*/
void iflip(int *a, int n)
{
	int i;
	int k;
	int tmp;
	
	k = n / 2;
	for (i = 1; i <= k; i++) {
		tmp = a[i];
		a[i] = a[n - i + 1];
		a[n - i + 1] = tmp;
	}
}		


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Bubble sort, which is most inefficient method, OK for small number n < 20.
   The old "sort" function 
   If option > 0, sort in ascending order 
   If option < 0, sort in decending order 	
   a[0...n-1] is input and replaced by sorted array as return.
   Note the original a[] is destroyed                                         */
/*----------------------------------------------------------------------------*/
int bubble(float *dat, int n, int option)
{
	int i, j;
	float tmp;
	
	void fflip(float *a, int n);
	
	if (option == 0) {
		fprintf(stderr, "bubble: Option == 0, not sort.\n");
		return -1;
	}
	
	for (i = 0; i < n-1; i++) {
		for (j = 0; j < n-1-i; j++) {
			if (dat[j] > dat[j+1]) {
				tmp = dat[j];
				dat[j] = dat[j+1];
				dat[j+1] = tmp;
			}
		}
	}
		
	if (option < 0) {
		fflip(dat-1, n);
	}
	
	return 0;		
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Bubble sort, which is most inefficient method, OK for small number n < 20.
   a[0...n-1] is input array.
   If option > 0, sort in ascending order, return index and sorted a[].
   If option < 0, sort in decending order, return index and sorted a[].
   Note the original a[] is destroyed                                         */
/*----------------------------------------------------------------------------*/
int bubble_indx(float *dat, int *index, int n, int option)
{
	int i, j;
	float tmp;
	int itmp;

	void fflip(float *a, int n);
	void iflip(int *a, int n);
		
	if (option == 0) {
		fprintf(stderr, "bubble_indx: Option == 0, not sort.\n");
		return -1;
	}

	for (i = 0; i < n; i++) {
		index[i] = i;
	}
				
	for (i = 0; i < n-1; i++) {
		for (j = 0; j < n-1-i; j++) {
			if (dat[j] > dat[j+1]) {
				tmp = dat[j];
				dat[j] = dat[j+1];
				dat[j+1] = tmp;
				itmp = index[j];
				index[j] = index[j+1];
				index[j+1] = itmp;
			}
		}
	}
	
	if (option < 0) {
		fflip(dat-1, n);
		iflip(index-1, n);
	}		
	
	return 0;		
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Sort a[1...n] using Shell's method (diminishing increment sort), 
	which is good for small number n < 50.
   a[1...n] is input array.
   If option > 0, sort in ascending order, return sorted a[].
   If option < 0, sort in decending order, return sorted a[].
   Note the original a[] is destroyed                                         */
/*----------------------------------------------------------------------------*/
int shell(float *a, int n, int option)
{
	int i, j, inc;
	float v;
	
	void fflip(float *a, int n);
	
	if (option == 0) {
		fprintf(stderr, "shell: Option == 0, not sort.\n");
		return -1;
	}
	
/*----------------------------------------------------------------------------*/
/* Determine the starting increment                                           */
/*----------------------------------------------------------------------------*/
	inc = 1;
	do {
		inc *= 3;
		inc++;
	} while (inc <= n);
	
/*----------------------------------------------------------------------------*/
/* Loop over the partial sorts                                                */
/*----------------------------------------------------------------------------*/
	do {
		inc /= 3;

/*----------------------------------------------------------------------------*/
/* Outer loop of straight insertion                                           */
/*----------------------------------------------------------------------------*/
		for (i = inc+1; i <= n; i++) {
			v = a[i];
			j = i;

/*----------------------------------------------------------------------------*/
/* Inner loop of straight insertion                                           */
/*----------------------------------------------------------------------------*/
			while (a[j - inc] > v) {
				a[j] = a[j - inc];
				j -= inc;
				if (j <= inc) break;
			}
			a[j] = v;
		}
	} while (inc > 1);
	
	if (option < 0) {
		fflip(a, n);
	}
	
	return 0;		
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Quick sorting,	which is the fastest sorting method for large number n.
   S[0...n-1] is input array.
   p usually is the first point index of S, p = 0.
   r usually is the last point index of S, r = n-1. 
   If option > 0, sort in ascending order, return sorted S[] and index.
   If option < 0, sort in decending order, return sorted S[] and index.
   Note: 
   (1) The original S[] is destroyed
   (2) Because it is an iterative algorithm, IDX must be assigned values 
   outside the quick subroutine. In practice, IDX[0...n-1] or IDX[1...n]      */
/*----------------------------------------------------------------------------*/
int quick(float *S, int *IDX, int p, int r, int option)
{
	int q;
	
	int partition(float *S, int *IDX, int p, int r, int option);
	
	if (option == 0) {
		fprintf(stderr, "quick: Option == 0, not sort.\n");
		return -1;
	}
		
	if (p < r) {
		q = partition(S, IDX, p, r, option);
		quick(S, IDX, p, q-1, option);
		quick(S, IDX, q+1, r, option);
	}			
	return 0;
	
}



/******************************************************************************/
/*----------------------------------------------------------------------------*/
/*	Partition in Quick Sorting                                                 */
/*----------------------------------------------------------------------------*/
int partition(float *S, int *IDX, int p, int r, int option)
{
	float X, tmp;
	int tmpid;
	int i, j;
	
	X = S[r];
	i = p-1;
	
	if (option > 0) { 
		for (j = p; j < r; j++) {
			if (S[j] < X) {
				i++;
				tmp = S[i];
				S[i] = S[j];
				S[j] = tmp;

				tmpid = IDX[i];
				IDX[i] = IDX[j];
				IDX[j] = tmpid;
			}
		}
	}	
	else if (option < 0) {
		for (j = p; j < r; j++) {
			if (S[j] > X) {
				i++;
				tmp = S[i];
				S[i] = S[j];
				S[j] = tmp;

				tmpid = IDX[i];
				IDX[i] = IDX[j];
				IDX[j] = tmpid;
			}
		}
	}			
	
	tmp = S[i+1];
	S[i+1] = S[r];
	S[r] = tmp;

	tmpid = IDX[i+1];
	IDX[i+1] = IDX[r];
	IDX[r] = tmpid;
	
	return i+1;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/*	Indexes an array arr[1..n] , i.e., outputs the array indx[1..n] 
	such that arr[indx[j]] is in ascending/decending order for j = 1,2,..,N. 
	The input quantities n and arr are not changed.
	This subroutine use quick sort algorithm, but different code as before.
	Note: if arr and indx are zero-offset, it works as calling with arr-1 and 
	indx-1, however, the value of indx[0...n-1] is still [1...n]               */
/*----------------------------------------------------------------------------*/
int quick_indx(float arr[], unsigned int indx[], unsigned int n, int option)
{
	unsigned int i, j, k, l = 1;
	unsigned int indxt;
	unsigned int ir = n;
	unsigned int itemp;
	int return_value = 0;
	int jstack = 0;
	int *istack = NULL;
	float a;

	if (option == 0) {
		fprintf(stderr, "quick_indx: Option == 0, not sort.\n");
		return -1;
	}

	istack = ivector(1, NSTACK);

	for (j = 1; j <= n; j++) indx[j] = j; 

	for (; ;) {
		if (ir - l < M) {
			for (j = l+1; j <= ir; j++) {
				indxt = indx[j]; 
				a = arr[indxt];
				for (i = j-1; i >= l; i--) {
					if (arr[indx[i]] <= a) break; 
					indx[i+1] = indx[i];
				}
				indx[i+1] = indxt;
			}

			if (jstack == 0) break; 
			ir = istack[jstack--]; 
			l = istack[jstack--];
		} 
		else {
			k = (l+ir) >> 1;
			SWAP(indx[k], indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l], indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1], indx[ir])
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l], indx[l+1])
			}
			i = l + 1; 
			j = ir;
			indxt = indx[l+1];
			a = arr[indxt]; 
			
			for (; ;) {
				do i++; while (arr[indx[i]] < a); 
				do j--; while (arr[indx[j]] > a); 
				if (j < i) break;
				SWAP(indx[i], indx[j])
			}

			indx[l+1] = indx[j]; 
			indx[j] = indxt; 
			jstack += 2;

			if (jstack > NSTACK) {
				nrerror("quick_indx: NSTACK too small in indexx."); 
				return_value = -1;
			}
			else {
				return_value = 0;
			}
			
			if (ir-i+1 >= j-l) { 
				istack[jstack] = ir; 
				istack[jstack-1] = i; 
				ir = j-1;	
			} 
			else {
				istack[jstack] = j-1; 
				istack[jstack-1] = l; 
				l=i;
			}
		}	
	}

	if (option < 0) {
		iflip((int*)indx, (int)n);
	}
	
	free_ivector(istack, 1, NSTACK);
	istack = NULL;
	return return_value;	
}
	

/******************************************************************************/
/*----------------------------------------------------------------------------*/
/*	Given indx[1..n] as output from the routine quick_indx(), 
	returns an array irank[1..n], the corresponding table of ranks.
	Note: if indx and irank are zero-offset, it works as calling with 
	indx-1 and irank-1, however, the value of irank[0...n-1] is still [1...n]  */
/*----------------------------------------------------------------------------*/
void quick_rank(unsigned int indx[], unsigned int irank[], unsigned int n)
{
	unsigned int j;

	for (j = 1; j <= n; j++) {
		irank[indx[j]] = j;
	}
}
											
