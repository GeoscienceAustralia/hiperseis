#include <math.h>
#include "complex.h"


/******************************************************************************/
/*----------------------------------------------------------------------------*/
fcomplex Complexf(const float& re, const float& im)
{	
	fcomplex c;
	
	c.r = re;
	c.i = im;
	return c;
}	
	
/*----------------------------------------------------------------------------*/
fcomplex Complexf(const double& re, const float& im)
{	
	fcomplex c;
	
	c.r = (float)re;
	c.i = im;
	return c;
}	
	
/*----------------------------------------------------------------------------*/
fcomplex Complexf(const float& re, const double& im)
{	
	fcomplex c;
	
	c.r = re;
	c.i = (float)im;
	return c;
}	
	
/*----------------------------------------------------------------------------*/
fcomplex Complexf(const double& re, const double& im)
{	
	fcomplex c;
	
	c.r = (float)re;
	c.i = (float)im;
	return c;
}	

	
/******************************************************************************/
/*----------------------------------------------------------------------------*/
dcomplex Complexd(const double& re, const double& im)
{	
	dcomplex c;
	
	c.r = re;
	c.i = im;
	return c;
}	

/*----------------------------------------------------------------------------*/
dcomplex Complexd(const double& re, const float& im)
{	
	dcomplex c;
	
	c.r = re;
	c.i = (double)im;
	return c;
}	

/*----------------------------------------------------------------------------*/
dcomplex Complexd(const float& re, const double& im)
{	
	dcomplex c;
	
	c.r = (double)re;
	c.i = im;
	return c;
}	

/*----------------------------------------------------------------------------*/
dcomplex Complexd(const float& re, const float& im)
{	
	dcomplex c;
	
	c.r = (double)re;
	c.i = (double)im;
	return c;
}	


/******************************************************************************/
/*---------------------------------------------------------------------------*/
fcomplex Complexd2f(const dcomplex& a)
{
	fcomplex c;
	
	c.r = (float)a.r;
	c.i = (float)a.i;
	return c;
}
	
/*----------------------------------------------------------------------------*/
dcomplex Complexf2d(const fcomplex& a)
{
	dcomplex c;
	
	c.r = (double)a.r;
	c.i = (double)a.i;
	return c;
}	


/******************************************************************************/
/*----------------------------------------------------------------------------*/
fcomplex Cadd(fcomplex a, fcomplex b)
{	
	fcomplex c;
	
	c.r = a.r + b.r;
	c.i = a.i + b.i;
	return c;
}

/*----------------------------------------------------------------------------*/
/* Added by Youlin Chen                                                       */									
/*----------------------------------------------------------------------------*/
fcomplex RCadd(float x, fcomplex a)
{	
	fcomplex c;
	
	c.r = a.r + x;
	c.i = a.i;
	return c;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
fcomplex Csub(fcomplex a, fcomplex b)
{	
	fcomplex c;
	
	c.r = a.r - b.r;
	c.i = a.i - b.i;
	return c;
}

/******************************************************************************/
/*----------------------------------------------------------------------------*/
/*fcomplex Cmul(fcomplex a, fcomplex b)
{	
	fcomplex c;
	
	c.r = a.r * b.r - a.i * b.i;
	c.i = a.i * b.r + a.r * b.i;
	return c;
}*/

/*----------------------------------------------------------------------------*/
fcomplex Cmul(fcomplex a, fcomplex b)
{	
	fcomplex c;
	float tmp1, tmp2;
	
	tmp1 = a.r * b.r;
	tmp2 = a.i * b.i;
	
	c.r = tmp1 - tmp2;
	c.i = (a.r + a.i) * (b.r + b.i) - tmp1 -tmp2;
	
	return c;
}

/*----------------------------------------------------------------------------*/
fcomplex RCmul(float x, fcomplex a)
{	
	fcomplex c;
	
	c.r = x * a.r;
	c.i = x * a.i;
	return c;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
fcomplex Cdiv(fcomplex a, fcomplex b)
{	
	fcomplex c;
	float r, den;
	
	if (fabs(b.r) >= fabs(b.i)) {
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = (a.r + r * a.i) / den;
		c.i = (a.i - r * a.r) / den;
	}
	else {
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = (a.r * r + a.i) / den;
		c.i = (a.i * r - a.r) / den;
	}
	return c;		
}

/*----------------------------------------------------------------------------*/
/* Added by Youlin Chen 									
	complex / float                                                            */
/*----------------------------------------------------------------------------*/
fcomplex CRdiv(fcomplex a, float b)
{	
	fcomplex c;
	
	c.r = a.r / b;
	c.i = a.i / b;
	return c;		
}

/*----------------------------------------------------------------------------*/
/* Added by Youlin Chen 									
	float / complex                                                            */
/*----------------------------------------------------------------------------*/
fcomplex RCdiv(float a, fcomplex b)
{	
	fcomplex c;
	float r, den;
	
	if (fabs(b.r) >= fabs(b.i)) {
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = a / den;
		c.i = (-r * a) / den;
	}
	else {
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = (a * r) / den;
		c.i = -a / den;
	}
	return c;		
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
fcomplex Conjg(const fcomplex& z)
{	
	fcomplex c;
	
	c.r = z.r;
	c.i = -z.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex Conjg(const dcomplex& z)
{	
	dcomplex c;
	
	c.r = z.r;
	c.i = -z.i;
	return c;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
float Cabs(const fcomplex& z)
{	
	float x, y, ans, temp;
	
	x = fabs(z.r);
	y = fabs(z.i);
	if (x == 0.0) ans = y;
	else if (y == 0.0) ans = x;
	else if (x > y) {
		temp = y / x;
		ans = x * sqrt(1.0 + temp * temp);
	}
	else {
		temp = x / y;
		ans = y * sqrt(1.0 + temp * temp);
	}
	return ans;
}
		
/*----------------------------------------------------------------------------*/
double Cabs(const dcomplex& z)
{	
	double x, y, ans, temp;
	
	x = fabs(z.r);
	y = fabs(z.i);
	if (x == 0.0) ans = y;
	else if (y == 0.0) ans = x;
	else if (x > y) {
		temp = y / x;
		ans = x * sqrt(1.0 + temp * temp);
	}
	else {
		temp = x / y;
		ans = y * sqrt(1.0 + temp * temp);
	}
	return ans;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
fcomplex Csqrt(const fcomplex& z)
{	
	fcomplex c;
	float x, y, w, r;
	
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r = 0.0;
		c.i = 0.0;
		return c;
	}
	else {
		x = fabs(z.r);
		y = fabs(z.i);
		if (x >= y) {
			r = y / x;
			w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r*r)));
		}
		else {
			r = x / y;
			w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r*r)));
		}
		if (z.r >= 0.0) {
			c.r = w;
			c.i = z.i / (2.0 * w);
		}
		else {
			c.i = (z.i >= 0.0) ? w : -w;
			c.r = z.i / (2.0 * c.i);
		}
		return c;
	}
}			

/*----------------------------------------------------------------------------*/
dcomplex Csqrt(const dcomplex& z)
{	
	dcomplex c;
	double x, y, w, r;
	
	if ((z.r == 0.0) && (z.i == 0.0)) {
		c.r = 0.0;
		c.i = 0.0;
		return c;
	}
	else {
		x = fabs(z.r);
		y = fabs(z.i);
		if (x >= y) {
			r = y / x;
			w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r*r)));
		}
		else {
			r = x / y;
			w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r*r)));
		}
		if (z.r >= 0.0) {
			c.r = w;
			c.i = z.i / (2.0 * w);
		}
		else {
			c.i = (z.i >= 0.0) ? w : -w;
			c.r = z.i / (2.0 * c.i);
		}
		return c;
	}
}			


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Added by Youlin Chen                                                       */									
/*----------------------------------------------------------------------------*/
fcomplex Cexp(const fcomplex& x)
{
	fcomplex z; 
	float r;
	
	r = exp(x.r);
	z.r = r * cos(x.i);
	z.i = r * sin(x.i);
	return z;
}

/*----------------------------------------------------------------------------*/
dcomplex Cexp(const dcomplex& x)
{
	dcomplex z; 
	double r;
	
	r = exp(x.r);
	z.r = r * cos(x.i);
	z.i = r * sin(x.i);
	return z;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Add by Youlin Chen
	operator overload                                                          */
/*----------------------------------------------------------------------------*/
fcomplex operator +(const fcomplex& a, const fcomplex& b) 
{
	fcomplex c;
	
	c.r = a.r + b.r;
	c.i = a.i + b.i;
	return c;
}

/*----------------------------------------------------------------------------*/
fcomplex operator +(const float& x, const fcomplex& a) 
{
	fcomplex c;
	
	c.r = a.r + x;
	c.i = a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
fcomplex operator +(const fcomplex& a, const float& x) 
{
	fcomplex c;

	c.r = a.r + x;
	c.i = a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator +(const dcomplex& a, const dcomplex& b) 
{
	dcomplex c;
	
	c.r = a.r + b.r;
	c.i = a.i + b.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator +(const dcomplex& a, const fcomplex& b) 
{
	dcomplex c;
	
	c.r = a.r + (double)b.r;
	c.i = a.i + (double)b.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator +(const fcomplex& a, const dcomplex& b) 
{
	dcomplex c;
	
	c.r = (double)a.r + b.r;
	c.i = (double)a.i + b.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator +(const double& x, const dcomplex& a) 
{
	dcomplex c;
	
	c.r = a.r + x;
	c.i = a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator +(const float& x, const dcomplex& a)
{
	dcomplex c;
	
	c.r = a.r + (double)x;
	c.i = a.i;
	return c;
}
 
/*----------------------------------------------------------------------------*/
dcomplex operator +(const double& x, const fcomplex& a) 
{
	dcomplex c;
	
	c.r = (double)a.r + x;
	c.i = (double)a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator +(const dcomplex& a, const double& x) 
{
	dcomplex c;

	c.r = a.r + x;
	c.i = a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator +(const dcomplex& a, const float& x)
{
	dcomplex c;

	c.r = a.r + (double)x;
	c.i = a.i;
	return c;
}
 
/*----------------------------------------------------------------------------*/
dcomplex operator +(const fcomplex& a, const double& x) 
{
	dcomplex c;
	
	c.r = (double)a.r + x;
	c.i = (double)a.i;
	return c;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
fcomplex operator -(const fcomplex& a) 
{
	fcomplex c;

	c.r = -a.r;
	c.i = -a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
fcomplex operator -(const fcomplex& a, const fcomplex& b)
{	
	fcomplex c;
	
	c.r = a.r - b.r;
	c.i = a.i - b.i;
	return c;
}

/*----------------------------------------------------------------------------*/
fcomplex operator -(const fcomplex& a, const float& x)
{	
	fcomplex c;
	
	c.r = a.r - x;
	c.i = a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
fcomplex operator -(const float& x, const fcomplex& a)
{	
	fcomplex c;
	
	c.r = x - a.r;
	c.i = -a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator -(const dcomplex& a) 
{
	dcomplex c;

	c.r = -a.r;
	c.i = -a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator -(const dcomplex& a, const dcomplex& b)
{	
	dcomplex c;
	
	c.r = a.r - b.r;
	c.i = a.i - b.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator -(const fcomplex& a, const dcomplex& b)
{	
	dcomplex c;
	
	c.r = (double)a.r - b.r;
	c.i = (double)a.i - b.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator -(const dcomplex& a, const fcomplex& b)
{	
	dcomplex c;
	
	c.r = a.r - (double)b.r;
	c.i = a.i - (double)b.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator -(const dcomplex& a, const double& x)
{	
	dcomplex c;
	
	c.r = a.r - x;
	c.i = a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator -(const dcomplex& a, const float& x)
{	
	dcomplex c;
	
	c.r = a.r - (double)x;
	c.i = a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator -(const fcomplex& a, const double& x)
{	
	dcomplex c;
	
	c.r = (double)a.r - x;
	c.i = (double)a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator -(const double& x, const dcomplex& a)
{	
	dcomplex c;
	
	c.r = x - a.r;
	c.i = -a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator -(const float& x, const dcomplex& a)
{	
	dcomplex c;
	
	c.r = (double)x - a.r;
	c.i = -a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator -(const double& x, const fcomplex& a)
{	
	dcomplex c;
	
	c.r = x - (double)a.r;
	c.i = -((double)(a.i));
	return c;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
fcomplex operator *(const fcomplex& a, const fcomplex& b)
{	
	fcomplex c;
	float tmp1, tmp2;
	
	tmp1 = a.r * b.r;
	tmp2 = a.i * b.i;
	
	c.r = tmp1 - tmp2;
	c.i = (a.r + a.i) * (b.r + b.i) - tmp1 -tmp2;
	
	return c;
}

/*----------------------------------------------------------------------------*/
fcomplex operator *(const float& x, const fcomplex& a)
{	
	fcomplex c;
	
	c.r = x * a.r;
	c.i = x * a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
fcomplex operator *(const fcomplex& a, const float& x)
{	
	fcomplex c;
	
	c.r = x * a.r;
	c.i = x * a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator *(const dcomplex& a, const dcomplex& b)
{	
	dcomplex c;
	double tmp1, tmp2;
	
	tmp1 = a.r * b.r;
	tmp2 = a.i * b.i;
	
	c.r = tmp1 - tmp2;
	c.i = (a.r + a.i) * (b.r + b.i) - tmp1 -tmp2;
	
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator *(const fcomplex& a, const dcomplex& b)
{	
	dcomplex c;
	double tmp1, tmp2;
	
	tmp1 = (double)a.r * b.r;
	tmp2 = (double)a.i * b.i;
	
	c.r = tmp1 - tmp2;
	c.i = ((double)a.r + (double)a.i) * (b.r + b.i) - tmp1 -tmp2;
	
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator *(const dcomplex& a, const fcomplex& b)
{	
	dcomplex c;
	double tmp1, tmp2;
	
	tmp1 = a.r * (double)b.r;
	tmp2 = a.i * (double)b.i;
	
	c.r = tmp1 - tmp2;
	c.i = (a.r + a.i) * ((double)b.r + (double)b.i) - tmp1 -tmp2;
	
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator *(const double& x, const dcomplex& a)
{	
	dcomplex c;
	
	c.r = x * a.r;
	c.i = x * a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator *(const float& x, const dcomplex& a)
{	
	dcomplex c;
	
	c.r = (double)x * a.r;
	c.i = (double)x * a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator *(const double& x, const fcomplex& a)
{	
	dcomplex c;
	
	c.r = x * (double)a.r;
	c.i = x * (double)a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator *(const dcomplex& a, const double& x)
{	
	dcomplex c;
	
	c.r = x * a.r;
	c.i = x * a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator *(const dcomplex& a, const float& x)
{	
	dcomplex c;
	
	c.r = (double)x * a.r;
	c.i = (double)x * a.i;
	return c;
}

/*----------------------------------------------------------------------------*/
dcomplex operator *(const fcomplex& a, const double& x)
{	
	dcomplex c;
	
	c.r = x * (double)a.r;
	c.i = x * (double)a.i;
	return c;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
fcomplex operator /(const fcomplex& a, const fcomplex& b)
{	
	fcomplex c;
	float r, den;
	
	if (fabs(b.r) >= fabs(b.i)) {
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = (a.r + r * a.i) / den;
		c.i = (a.i - r * a.r) / den;
	}
	else {
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = (a.r * r + a.i) / den;
		c.i = (a.i * r - a.r) / den;
	}
	return c;		
}

/*----------------------------------------------------------------------------*/
fcomplex operator /(const fcomplex& a, const float& b)
{	
	fcomplex c;
	
	c.r = a.r / b;
	c.i = a.i / b;
	return c;		
}

/*----------------------------------------------------------------------------*/
fcomplex operator /(const float& a, const fcomplex& b)
{	
	fcomplex c;
	float r, den;
	
	if (fabs(b.r) >= fabs(b.i)) {
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = a / den;
		c.i = (-r * a) / den;
	}
	else {
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = (a * r) / den;
		c.i = -a / den;
	}
	return c;		
}

/*----------------------------------------------------------------------------*/
dcomplex operator /(const dcomplex& a, const dcomplex& b)
{	
	dcomplex c;
	double r, den;
	
	if (fabs(b.r) >= fabs(b.i)) {
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = (a.r + r * a.i) / den;
		c.i = (a.i - r * a.r) / den;
	}
	else {
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = (a.r * r + a.i) / den;
		c.i = (a.i * r - a.r) / den;
	}
	return c;		
}

/*----------------------------------------------------------------------------*/
dcomplex operator /(const dcomplex& a, const fcomplex& b)
{	
	dcomplex c;
	double r, den;
	
	if (fabs((double)b.r) >= fabs((double)b.i)) {
		r = (double)b.i / (double)b.r;
		den = (double)b.r + r * (double)b.i;
		c.r = (a.r + r * a.i) / den;
		c.i = (a.i - r * a.r) / den;
	}
	else {
		r = (double)b.r / (double)b.i;
		den = (double)b.i + r * (double)b.r;
		c.r = (a.r * r + a.i) / den;
		c.i = (a.i * r - a.r) / den;
	}
	return c;		
}

/*----------------------------------------------------------------------------*/
dcomplex operator /(const fcomplex& a, const dcomplex& b)
{	
	dcomplex c;
	double r, den;
	
	if (fabs(b.r) >= fabs(b.i)) {
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = ((double)a.r + r * (double)a.i) / den;
		c.i = ((double)a.i - r * (double)a.r) / den;
	}
	else {
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = ((double)a.r * r + (double)a.i) / den;
		c.i = ((double)a.i * r - (double)a.r) / den;
	}
	return c;		
}

/*----------------------------------------------------------------------------*/
dcomplex operator /(const dcomplex& a, const double& x)
{	
	dcomplex c;
	
	c.r = a.r / x;
	c.i = a.i / x;
	return c;		
}

/*----------------------------------------------------------------------------*/
dcomplex operator /(const dcomplex& a, const float& x)
{	
	dcomplex c;
	
	c.r = a.r / (double)x;
	c.i = a.i / (double)x;
	return c;		
}

/*----------------------------------------------------------------------------*/
dcomplex operator /(const fcomplex& a, const double& x)
{	
	dcomplex c;
	
	c.r = (double)a.r / x;
	c.i = (double)a.i / x;
	return c;		
}

/*----------------------------------------------------------------------------*/
dcomplex operator /(const double& x, const dcomplex& b)
{	
	dcomplex c;
	double r, den;
	
	if (fabs(b.r) >= fabs(b.i)) {
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = x / den;
		c.i = (-r * x) / den;
	}
	else {
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = (x * r) / den;
		c.i = -x / den;
	}
	return c;		
}

/*----------------------------------------------------------------------------*/
dcomplex operator /(const float& x, const dcomplex& b)
{	
	dcomplex c;
	double r, den;
	
	if (fabs(b.r) >= fabs(b.i)) {
		r = b.i / b.r;
		den = b.r + r * b.i;
		c.r = (double)x / den;
		c.i = (-r * (double)x) / den;
	}
	else {
		r = b.r / b.i;
		den = b.i + r * b.r;
		c.r = ((double)x * r) / den;
		c.i = -(double)x / den;
	}
	return c;		
}

/*----------------------------------------------------------------------------*/
dcomplex operator /(const double& x, const fcomplex& b)
{	
	dcomplex c;
	double r, den;
	
	if (fabs((double)b.r) >= fabs((double)b.i)) {
		r = (double)b.i / (double)b.r;
		den = (double)b.r + r * (double)b.i;
		c.r = x / den;
		c.i = (-r * x) / den;
	}
	else {
		r = (double)b.r / (double)b.i;
		den = (double)b.i + r * (double)b.r;
		c.r = (x * r) / den;
		c.i = -x / den;
	}
	return c;		
}
