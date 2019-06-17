#ifndef _COMPLEX_H_
#define _COMPLEX_H_

/*----------------------------------------------------------------------------*/
/*	Note: 
	Use C++ compiler because function overloading and operator overloading
	are used in function bodys.
	Asignment = cannot be done because it must be non-static member function.
	fcomplex = fcomplex and dcomplex = dcomplex are done implicitly.
	fcomplex = dcomplex and dcomplex = fcomplex must be done explicitly from their
	real and imagine parts using Complexf2d() and Complexd2f() functions.      */
/*----------------------------------------------------------------------------*/
			 
typedef struct {	
	float r;
	float i;
} 
fcomplex;

typedef struct {	
	double r;
	double i;
} 
dcomplex;

/*----------------------------------------------------------------------------*/
/* C functions                                                                */
/*----------------------------------------------------------------------------*/
fcomplex Cadd(fcomplex a, fcomplex b);
fcomplex RCadd(float x, fcomplex a);
fcomplex Csub(fcomplex a, fcomplex b);
fcomplex Cmul(fcomplex a, fcomplex b);
fcomplex RCmul(float x, fcomplex a);
fcomplex Cdiv(fcomplex a, fcomplex b);
fcomplex CRdiv(fcomplex a, float b);
fcomplex RCdiv(float a, fcomplex b);

/*----------------------------------------------------------------------------*/
/* Function overloading for C++                                               */
/*----------------------------------------------------------------------------*/
fcomplex Complexf(const float& re, const float& im);
fcomplex Complexf(const double& re, const float& im);
fcomplex Complexf(const float& re, const double& im);
fcomplex Complexf(const double& re, const double& im);
dcomplex Complexd(const double& re, const double& im);
dcomplex Complexd(const double& re, const float& im);
dcomplex Complexd(const float& re, const double& im);
dcomplex Complexd(const float& re, const float& im);
fcomplex Complexd2f(const dcomplex& a);
dcomplex Complexf2d(const fcomplex& a);

fcomplex Conjg(const fcomplex& z);
dcomplex Conjg(const dcomplex& z);
float Cabs(const fcomplex& z);
double Cabs(const dcomplex& z);
fcomplex Csqrt(const fcomplex& z);
dcomplex Csqrt(const dcomplex& z);
fcomplex Cexp(const fcomplex& x);
dcomplex Cexp(const dcomplex& x);

/*----------------------------------------------------------------------------*/
/* Operator overloading for C++. They are call-by-reference.                  */
/*----------------------------------------------------------------------------*/
fcomplex operator +(const fcomplex& a, const fcomplex& b); 
fcomplex operator +(const float& x, const fcomplex& a); 
fcomplex operator +(const fcomplex& a, const float& x); 
dcomplex operator +(const dcomplex& a, const dcomplex& b); 
dcomplex operator +(const dcomplex& a, const fcomplex& b);
dcomplex operator +(const fcomplex& a, const dcomplex& b);
dcomplex operator +(const double& x, const dcomplex& a); 
dcomplex operator +(const float& x, const dcomplex& a); 
dcomplex operator +(const double& x, const fcomplex& a); 
dcomplex operator +(const dcomplex& a, const double& x); 
dcomplex operator +(const dcomplex& a, const float& x); 
dcomplex operator +(const fcomplex& a, const double& x); 

fcomplex operator -(const fcomplex& a); 
fcomplex operator -(const fcomplex& a, const fcomplex& b);
fcomplex operator -(const fcomplex& a, const float& x);
fcomplex operator -(const float& x, const fcomplex& a);
dcomplex operator -(const dcomplex& a); 
dcomplex operator -(const dcomplex& a, const dcomplex& b);
dcomplex operator -(const fcomplex& a, const dcomplex& b);
dcomplex operator -(const dcomplex& a, const fcomplex& b);
dcomplex operator -(const dcomplex& a, const double& x);
dcomplex operator -(const dcomplex& a, const float& x);
dcomplex operator -(const fcomplex& a, const double& x);
dcomplex operator -(const double& x, const dcomplex& a);
dcomplex operator -(const float& x, const dcomplex& a);
dcomplex operator -(const double& x, const fcomplex& a);

fcomplex operator *(const fcomplex& a, const fcomplex& b);
fcomplex operator *(const float& x, const fcomplex& a);
fcomplex operator *(const fcomplex& a, const float& x);
dcomplex operator *(const dcomplex& a, const dcomplex& b);
dcomplex operator *(const fcomplex& a, const dcomplex& b);
dcomplex operator *(const dcomplex& a, const fcomplex& b);
dcomplex operator *(const double& x, const dcomplex& a);
dcomplex operator *(const float& x, const dcomplex& a);
dcomplex operator *(const double& x, const fcomplex& a);
dcomplex operator *(const dcomplex& a, const double& x);
dcomplex operator *(const dcomplex& a, const float& x);
dcomplex operator *(const fcomplex& a, const double& x);

fcomplex operator /(const fcomplex& a, const fcomplex& b);
fcomplex operator /(const fcomplex& a, const float& b);
fcomplex operator /(const float& a, const fcomplex& b);
dcomplex operator /(const dcomplex& a, const dcomplex& b);
dcomplex operator /(const dcomplex& a, const fcomplex& b);
dcomplex operator /(const fcomplex& a, const dcomplex& b);
dcomplex operator /(const dcomplex& a, const double& x);
dcomplex operator /(const dcomplex& a, const float& x);
dcomplex operator /(const fcomplex& a, const double& x);
dcomplex operator /(const double& x, const dcomplex& b);
dcomplex operator /(const float& x, const dcomplex& b);
dcomplex operator /(const double& x, const fcomplex& b);

#endif /* _COMPLEX_H_ */
