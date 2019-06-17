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

#define PI 3.141592653589793
#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr

/*----------------------------------------------------------------------------*/
/* If 'isign' is input as 1, it replaces data[1...2*nn] by its discrete 
	Fourier transform. 
	If 'isign' is input as -1, it replaces data[1...2*nn] by 'nn' times 
	its inverse transform.
	data is a complex array of length nn or, equivalently, a real array of length 2*nn.
	nn MUST be an integer power of 2 (this is NOT checked for!).               */
/*----------------------------------------------------------------------------*/
void four1(float data[], unsigned int nn, int isign)
{
	unsigned int n, mmax, m, j, istep, i;
/*----------------------------------------------------------------------------*/
/* Double for the trigonometric recurrences                                   */
/*----------------------------------------------------------------------------*/
	double wtemp, wr, wpr, wpi, wi, theta; 
	float tempr, tempi;

	n = nn << 1;
	j = 1;

/*----------------------------------------------------------------------------*/
/* This is the bit-reversal section of the routine.
	Exchange the two complex numbers.                                          */
/*----------------------------------------------------------------------------*/
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			SWAP(data[j], data[i]);
			SWAP(data[j+1], data[i+1]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	
/*----------------------------------------------------------------------------*/
/* Here begins the Danielson-Lanczos section of the routine                   */
/*----------------------------------------------------------------------------*/
	mmax = 2;
/*----------------------------------------------------------------------------*/
/* Outer loop executed log2(nn) times                                         */
/*----------------------------------------------------------------------------*/
	while (n > mmax) {                            
		istep = mmax << 1;
/*----------------------------------------------------------------------------*/
/* Initialize the trigonometric recurrence                                    */
/*----------------------------------------------------------------------------*/
		theta = isign * (2.0 * PI / mmax); 
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
/*----------------------------------------------------------------------------*/
/* Here are the two nested inner loops                                        */
/*----------------------------------------------------------------------------*/
		for (m = 1; m < mmax; m += 2) {            
			for (i = m; i <= n; i += istep) {
/*----------------------------------------------------------------------------*/
/* Follow is the Danielson-Lanczos formula                                    */
/*----------------------------------------------------------------------------*/
				j = i + mmax;    
				tempr = wr * data[j] - wi * data[j+1];
				tempi = wr * data[j+1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
/*----------------------------------------------------------------------------*/
/* Trigonometric recurrence                                                   */
/*----------------------------------------------------------------------------*/
			wr = (wtemp=wr)*wpr - wi*wpi + wr;    
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}
	
}

	
/*----------------------------------------------------------------------------*/
/* Given two real input arrays data[1...n] and data2[1...n].
	It calls four1 and returns two complex output arrays, 
	fft1[1...2n] and fft2[1..2n].
	Each of complex length n (i.e., real length 2*n) contains the discrete 
	Fourier transforms of the respective data arrays.
	n MUST be an integer power of 2.                                           */
/*----------------------------------------------------------------------------*/
void twofft(float data1[], float data2[], float fft1[], float fft2[], unsigned int n)
{
	void four1(float data[], unsigned int nn, int isign);
	
	unsigned int nn3, nn2, jj, j;
	float rep, rem, aip, aim;
	
	nn3 = 1 + (nn2 = 2 + n + n);
	
/*----------------------------------------------------------------------------*/
/* Pack the two real arrays into one complex array.                           */
/*----------------------------------------------------------------------------*/
	for (j = 1, jj = 2; j <= n; j++, jj += 2) {
		fft1[jj-1] = data1[j];
		fft1[jj] = data2[j];
	}
	
/*----------------------------------------------------------------------------*/
/* Transform the complex array.                                               */
/*----------------------------------------------------------------------------*/
	four1(fft1, n, 1);
	fft2[1] = fft1[2];
	fft1[2] = fft2[2] = 0.0;
	for (j = 3; j <= n + 1; j += 2) {
/*----------------------------------------------------------------------------*/
/* Use symmetries to separate the two transforms                              */
/*----------------------------------------------------------------------------*/
		rep = 0.5 * (fft1[j] + fft1[nn2-j]);  
		rem = 0.5 * (fft1[j] - fft1[nn2-j]);					
		aip = 0.5 * (fft1[j+1] + fft1[nn3-j]);
		aim = 0.5 * (fft1[j+1] - fft1[nn3-j]);
		
/*----------------------------------------------------------------------------*/
/* Ship them out in two complex arrays                                        */ 
/*----------------------------------------------------------------------------*/
		fft1[j] = rep;
		fft1[j+1] = aim;
		fft1[nn2-j] = rep;
		fft1[nn3-j] = -aim;
		fft2[j] = aip;
		fft2[j+1] = -rem;
		fft2[nn2-j] = aip;
		fft2[nn3-j] = rem;
	}
}


/*----------------------------------------------------------------------------*/
/* Calculate the Fourier transform of a set of n real data points stored in 
	data[1...n].
	Replace this data by the positive frequency half of its complex Fourier transform.
	The real-valued first and last components of the complex transform are 
	returned as elements data[1] and data[2], respectively.
	n MUST be a power of 2.
	This routine also calculates the inverse transform of a complex data array 
	if it is the transform of real data. 
	(Result in this case must be multiplied by 2/n.)                           */
/*----------------------------------------------------------------------------*/
void realft(float data[], unsigned int n, int isign)
{
	void four1(float data[], unsigned int nn, int isign);
	
	unsigned int i, i1, i2, i3, i4, np3;
	float c1=0.5, c2, h1r, h1i, h2r, h2i;
/*----------------------------------------------------------------------------*/
/* Double for the trigonometric recurrences                                   */
/*----------------------------------------------------------------------------*/
	double wr, wi, wpr, wpi, wtemp, theta; 
	
/*----------------------------------------------------------------------------*/
/* Initialize the recurrence                                                  */
/*----------------------------------------------------------------------------*/
	theta = PI / (double)(n >> 1);  
/*----------------------------------------------------------------------------*/
/* The forward transform                                                      */
/*----------------------------------------------------------------------------*/
	if (isign == 1) {                              
		c2 = -0.5;
		four1(data, n>>1, 1);
	}
/*----------------------------------------------------------------------------*/
/* The inverse transform                                                      */
/*----------------------------------------------------------------------------*/
	else {                                        
		c2 = 0.5;
		theta = -theta;
	}		
	
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0 + wpr;
	wi = wpi;
	np3 = n + 3;
/*----------------------------------------------------------------------------*/
/* Case i = 1 will be done separately below                                   */
/*----------------------------------------------------------------------------*/
	for (i = 2; i <= (n>>2); i++) {               
		i1 = i + i - 1;
		i2 = i1 + 1;
		i3 = np3 - i2;
		i4 = i3 + 1;
		
/*----------------------------------------------------------------------------*/
/* The two separate transforms are separated out of data.                     */
/*----------------------------------------------------------------------------*/
		h1r = c1 * (data[i1] + data[i3]);
		h1i = c1 * (data[i2] - data[i4]);
		h2r = -c2 * (data[i2] + data[i4]);
		h2i = c2 * (data[i1] - data[i3]);
		
/*----------------------------------------------------------------------------*/
/* They are recombined to form the true transform of the original real data.  */
/*----------------------------------------------------------------------------*/
		data[i1] = h1r + wr * h2r - wi * h2i;
		data[i2] = h1i + wr * h2i + wi * h2r;
		data[i3] = h1r - wr * h2r + wi * h2i;
		data[i4] = -h1i + wr * h2i + wi * h2r;
		
/*----------------------------------------------------------------------------*/
/* The recurrence.                                                            */
/*----------------------------------------------------------------------------*/
		wr = (wtemp = wr) * wpr - wi * wpi + wr;
		wi = wi * wpr + wtemp * wpi + wi;
	}	
	
/*----------------------------------------------------------------------------*/
/* Squeeze the first and last data together to get them all within the 
	original array.                                                            */
/*----------------------------------------------------------------------------*/
	if (isign == 1) {
		data[1] = (h1r = data[1]) + data[2];
		data[2] = h1r - data[2];
	}
	else {
		data[1] = c1 * ((h1r = data[1]) + data[2]);
		data[2] = c1 * (h1r - data[2]);
		four1(data, n >> 1, -1);
	}
	
}			
	
	
/*----------------------------------------------------------------------------*/
/* This version is exactly the same as Fortran version in processing data numbers.
	However, due to 'n' complex number == '2n' real number, 
	it may not be used properly.                                               */
/*----------------------------------------------------------------------------*/
void realft2(float data[], unsigned int n, int isign)
{
	void four1(float data[], unsigned int nn, int isign);
	
	unsigned int i, i1, i2, i3, i4, np3;
	float c1=0.5, c2, h1r, h1i, h2r, h2i;
/*----------------------------------------------------------------------------*/
/* Double for the trigonometric recurrences                                   */
/*----------------------------------------------------------------------------*/
	double wr, wi, wpr, wpi, wtemp, theta; 
	
/*----------------------------------------------------------------------------*/
/* Initialize the recurrence                                                  */
/*----------------------------------------------------------------------------*/
	theta = PI / (double)(n);  
/*----------------------------------------------------------------------------*/
/* The forward transform                                                      */
/*----------------------------------------------------------------------------*/
	if (isign == 1) {                              
		c2 = -0.5;
		four1(data, n, 1);
	}
/*----------------------------------------------------------------------------*/
/* The inverse transform                                                      */
/*----------------------------------------------------------------------------*/
	else {                                         
		c2 = 0.5;
		theta = -theta;
	}		
	
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0 + wpr;
	wi = wpi;
	np3 = 2*n + 3;
/*----------------------------------------------------------------------------*/
/* Case i = 1 will be done separately below                                   */
/*----------------------------------------------------------------------------*/
	for (i = 2; i <= (n>>2)+1; i++) {               
		i1 = i + i - 1;
		i2 = i1 + 1;
		i3 = np3 - i2;
		i4 = i3 + 1;
		
/*----------------------------------------------------------------------------*/
/* The two separate transforms are separated out of data.                     */
/*----------------------------------------------------------------------------*/
		h1r = c1 * (data[i1] + data[i3]);
		h1i = c1 * (data[i2] - data[i4]);
		h2r = -c2 * (data[i2] + data[i4]);
		h2i = c2 * (data[i1] - data[i3]);
		
/*----------------------------------------------------------------------------*/
/* They are recombined to form the true transform of the original real data.  */
/*----------------------------------------------------------------------------*/
		data[i1] = h1r + wr * h2r - wi * h2i;
		data[i2] = h1i + wr * h2i + wi * h2r;
		data[i3] = h1r - wr * h2r + wi * h2i;
		data[i4] = -h1i + wr * h2i + wi * h2r;
		
/*----------------------------------------------------------------------------*/
/* The recurrence.                                                            */
/*----------------------------------------------------------------------------*/
		wr = (wtemp = wr) * wpr - wi * wpi + wr;
		wi = wi * wpr + wtemp * wpi + wi;
	}	
	
/*----------------------------------------------------------------------------*/
/* Squeeze the first and last data together to get them all within 
	the original array.                                                        */
/*----------------------------------------------------------------------------*/
	if (isign == 1) {
		data[1] = (h1r = data[1]) + data[2];
		data[2] = h1r - data[2];
	}
	else {
		data[1] = c1 * ((h1r = data[1]) + data[2]);
		data[2] = c1 * (h1r - data[2]);
		four1(data, n, -1);
	}
	
}			
	
	
/*----------------------------------------------------------------------------*/
/* Compute the correlation of two real data sets data1[1...n] and data2[1...n].
	n MUST be an integer power of two, so the calculation includes zero padding 
	if necessary
	The answer is returned as the first n points in ans[1...2*n] stored in 
	wrap-around order, i.e., correlations at increasingly negative lags are in 
	ans[n] on down to ans[n/2+1], while correlations at increasingly positive 
	lags are in ans[1] (zero lag) on up to ans[n/2].
	Note that ans must be supplied in the calling program with length at least 2*n, 
	since it is also used as working space. 
	Sign convention of this routine: if data1 lags data2, i.e., is shifted to 
	the right of it, then ans will show a peak at positive lags.               */
/*----------------------------------------------------------------------------*/
void correl(float data1[], float data2[], unsigned int n, float ans[])
{
	void realft(float data[], unsigned int n, int isign);
	void realft2(float data[], unsigned int n, int isign);
	void twofft(float data1[], float data2[], float fft1[], float fft2[], unsigned int n);
	
	unsigned int no2, i;
	float dum, *fft;
		
	fft = vector(1, n << 1);
	
/*----------------------------------------------------------------------------*/
/* Transform both data vectors at once.                                       */
/*----------------------------------------------------------------------------*/
	twofft(data1, data2, fft, ans, n);
	
/*----------------------------------------------------------------------------*/
/* Normalization for inverse FFT.                                             */
/*----------------------------------------------------------------------------*/
	no2 = n >> 1;
/*----------------------------------------------------------------------------*/
/* Multiply to find FFT of their correlation.                                 */
/*----------------------------------------------------------------------------*/
	for (i = 2; i <= n + 2; i += 2) {
		ans[i-1] = (fft[i-1] * (dum = ans[i-1]) + fft[i] * ans[i]) / no2;
		ans[i] = (fft[i] * dum - fft[i-1] * ans[i]) / no2;
	}
	
/*----------------------------------------------------------------------------*/
/* Pack first and last into one element.                                      */
/*----------------------------------------------------------------------------*/
	ans[2] = ans[n+1];            
	
/*----------------------------------------------------------------------------*/
/* Inverse transform gives correlation.                                       */
/*----------------------------------------------------------------------------*/
	realft(ans, n, -1);           
	
	free_vector(fft, 1, n << 1);
}		

	
/*----------------------------------------------------------------------------*/
/* Convolves or deconvolves a real data set data[1,...,n] 
	(including any user-supplied zero padding) with a response function respns[1,...,n]. 
	The response function must be stored in wrap-around order in the fist m elements 
	of respns, where m is an odd integer <= n. 
	Wrap-around order means that the first half of the array respns contains 
	the impulse response function at positive times, while the second half of 
	the array contains the impluse response function at negative times,
	counting down from the highest element respns[m]. 
	On input isign is +1 for convolution, -1 for deconvolution. 
	The answer is returned in the first n components of ans. 
	However,	ans must be supplied in the calling program with dimensions [1,...,2*n], 
	for consistency with twofft. 
	n MUST be an integer power of two.                                         */
/*----------------------------------------------------------------------------*/
void convlv(float data[], unsigned int n, float respns[], unsigned int m, int isign, float ans[])
{
	void realft(float data[], unsigned int n, int isign);
	void twofft(float data1[], float data2[], float fft1[], float fft2[], unsigned int n);
	unsigned int i, no2;
	float dum, mag2, *fft;		
	
	fft = vector(1, n<<1);
/*----------------------------------------------------------------------------*/
/* Put respns in array of length n.                                           */
/*----------------------------------------------------------------------------*/
	for (i = 1; i <= (m-1)/2; i++)         
		respns[n+1-i] = respns[m+1-i];
/*----------------------------------------------------------------------------*/
/* Pad with zeros.                                                            */
/*----------------------------------------------------------------------------*/
	for (i = (m+3)/2; i <= n-(m-1)/2; i++) 
		respns[i] = 0.0;
		
/*----------------------------------------------------------------------------*/
/* FFT both at once                                                           */
/*----------------------------------------------------------------------------*/
	twofft(data, respns, fft, ans, n);     
	
	no2 = n >> 1;
	for (i = 2; i <= n+2; i+=2) {
		if (isign == 1) {
/*----------------------------------------------------------------------------*/
/* Multiple FFTs to convolve.                                                 */
/*----------------------------------------------------------------------------*/
			ans[i-1] = (fft[i-1] * (dum=ans[i-1]) - fft[i] * ans[i]) / no2;  
			ans[i] = (fft[i] * dum + fft[i-1] * ans[i]) / no2;
		}
		else if (isign == -1) {
			if ((mag2=sqrt(ans[i-1]) + sqrt(ans[i])) == 0.0)
				nrerror("Deconvolving at response zero in convlv");
/*----------------------------------------------------------------------------*/
/* Divide FFTs to deconvolve.                                                 */
/*----------------------------------------------------------------------------*/
			ans[i-1] = (fft[i-1] * (dum=ans[i-1]) + fft[i] * ans[i]) / mag2 / no2;  
			ans[i] = (fft[i] * dum - fft[i-1] * ans[i]) / mag2 / no2;
		}
		else
			nrerror("No meaning for isign in convlv");
	}
	
/*----------------------------------------------------------------------------*/
/* Pack last element with first for realft.                                   */
/*----------------------------------------------------------------------------*/
	ans[2] = ans[n+1];   

/*----------------------------------------------------------------------------*/
/* Inverse transform back to time domain.                                     */
/*----------------------------------------------------------------------------*/
	realft(ans, n, -1); 	

	free_vector(fft, 1, n<<1);
}			
				 			
