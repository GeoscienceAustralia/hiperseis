#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#ifndef SQR
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#endif

#ifndef DSQR
static double dsqrarg;
#define DSQR(a)  ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
#endif

#ifndef DMAX
static double dmaxarg1, dmaxarg2;
#define DMAX(a, b) (dmaxarg1=(a), dmaxarg2=(b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
#endif

#ifndef DMIN
static double dminarg1, dminarg2;
#define DMIN(a, b) (dminarg1=(a), dminarg2=(b), (dminarg1) < (dminarg2) ? (dminarg1) : (dminarg2))
#endif

#ifndef FMAX
static float maxarg1, maxarg2;
#define FMAX(a, b) (maxarg1=(a), maxarg2=(b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#endif

#ifndef FMIN
static float minarg1, minarg2;
#define FMIN(a, b) (minarg1=(a), minarg2=(b), (minarg1) < (minarg2) ? (minarg1) : (minarg2))
#endif

#ifndef LMAX
static long lmaxarg1, lmaxarg2;
#define LMAX(a, b) (lmaxarg1=(a), lmaxarg2=(b), (lmaxarg1) > (lmaxarg2) ? (lmaxarg1) : (lmaxarg2))
#endif

#ifndef LMIN
static long lminarg1, lminarg2;
#define LMIN(a, b) (lminarg1=(a), lminarg2=(b), (lminarg1) < (lminarg2) ? (lminarg1) : (lminarg2))
#endif

#ifndef IMAX
static int imaxarg1, imaxarg2;
#define IMAX(a, b) (imaxarg1=(a), imaxarg2=(b), (imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))
#endif

#ifndef IMIN
static int iminarg1, iminarg2;
#define IMIN(a, b) (iminarg1=(a), iminarg2=(b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
#endif

#ifndef SIGN
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#endif

#ifndef SIGN1
#define SIGN1(x) ((x) >= 0.0 ? 1.0 : -1.0)
#endif


/*----------------------------------------------------------------------------*/
struct NRerror {
	char *message;
	char *file;
	int line;
	NRerror(char *m, char *f, int l) : message(m), file(f), line(l) {}
};		


/*----------------------------------------------------------------------------*/
void nrerror(char error_text[]);
void NRcatch(NRerror err);


#endif /* _NR_UTILS_H_ */
