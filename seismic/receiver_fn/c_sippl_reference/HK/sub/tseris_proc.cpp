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

int xcorr(float *x, float *y, const int& n, const int& maxdelay, 
			float *r, int *offshift, const int& opt)

{
	int i,j, delay_time;
	double mx, my, sx, sy, sxy, denom;
   
	/* Calculate the mean of the two series x[], y[] */
	mx = 0;
	my = 0;
	if (opt == 1) {   
		for (i = 0; i < n; i++) {
			mx += x[i];
			my += y[i];
		}
		mx /= n;
		my /= n;
	}

	/* Calculate the denominator */
	sx = 0;
	sy = 0;
	for (i = 0; i < n; i++) {
		sx += (x[i] - mx) * (x[i] - mx);
		sy += (y[i] - my) * (y[i] - my);
	}
	denom = sqrt(sx*sy);

	/* Calculate the correlation series */
	for (delay_time = -maxdelay; delay_time < maxdelay; delay_time++) {
		sxy = 0;
		for (i = 0; i < n; i++) {
			j = i + delay_time;
			if (j < 0 || j >= n)
				sxy += (x[i] - mx) * (-my);
			else
				sxy += (x[i] - mx) * (y[j] - my);
		}
		offshift[delay_time + maxdelay] = delay_time;
		r[delay_time + maxdelay] = sxy / denom;      
	}

	return 0;
}

//---------------------------------------------------------------------------
float rms(float *x, float *y, const int& n)
{
	int i;
	double r;
	
	r = 0.0;
	for (i = 0; i < n; i++) {
		r += (x[i]-y[i]) * (x[i]-y[i]);
	}
	
	r = sqrt(r/n);
	
	return (float)r;
}
	

//---------------------------------------------------------------------------
float vrd(float *x, float *y, const int& n)
{
	int i;
	double r, power, power1, power2;
	
	r = 0.0;
	power1 = 0.0;
	power2 = 0.0;
	for (i = 0; i < n; i++) {
		r += (x[i]-y[i]) * (x[i]-y[i]);
		power1 += x[i] * x[i];
		power2 += y[i] * y[i];
	}
	
	power = DMAX(power1, power2);
	if (power == 0.0) {
		r = 0.0;
	}	
	else {
		r = sqrt(r / power);
	}
	
	return (float)r;
}
	

