#include "recordrf.h"
#include "sac.h"
#include "nrutil.h"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>


#ifndef MAX_MSFUNC
#define MAX_MSFUNC 300
#endif

#ifndef	_ZCC_L_
#define	_ZCC_L_		5.     	/* For Moho depth in km	*/
#endif

#ifndef	_CD_FACTOR_
#define	_CD_FACTOR_	4	
#endif

#ifndef	_RCC_L_
#define	_RCC_L_		0.2     	/* For Vp/Vs	*/
#endif


int segchar(char* inchar, char* charp1, char* charp2, char token, int seq);

/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Calculate travel time, pick up amplitude 	
	Stacking all stations                                                      */
/*----------------------------------------------------------------------------*/
void nth_stk(RECORD *rfrec, HKS_INFO *hks, const int& bootstrap)
{	
	int izr;
	int irec;
	int i;
	double tmpx1, tmpx2, tmpx3;
	double sa_0p1s, sa_2p1s, sa_1p2s;
		
/*----------------------------------------------------------------------------*/
/* From here, fork bootstrap and non-bootstrap
	non-bootstrap = !btstrp.flag  
	This is most time consuming part                                           */
/*----------------------------------------------------------------------------*/
	if (!bootstrap) {
		
/*----------------------------------------------------------------------------*/
/* linear stacking                                                            */
/*----------------------------------------------------------------------------*/
		if (hks->nroot == 1) {  
			for (irec = 0; irec < hks->nrfrec; irec++) {
				rfrec[irec].LoadData();
				printf("Processing rfrec[%d/%d]...\n", irec, hks->nrfrec);

				for (izr = 0; izr < hks->nzr; izr++) {
					rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array
					
					rfrec[irec].PeakAmp(hks);

					//if (rfrec[irec].a_0p1s == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					//if (rfrec[irec].a_2p1s == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					//if (rfrec[irec].a_1p2s == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					
					hks->a_0p1s[izr] += rfrec[irec].a_0p1s;
					hks->a_2p1s[izr] += rfrec[irec].a_2p1s;
					hks->a_1p2s[izr] += rfrec[irec].a_1p2s;
				}
				rfrec[irec].FreeData();
			}
			
			for (izr = 0; izr < hks->nzr; izr++) {					
				hks->a_0p1s[izr] /= hks->nrfrec;
				hks->a_2p1s[izr] /= hks->nrfrec;
				hks->a_1p2s[izr] /= hks->nrfrec;	
			}
		}
/*----------------------------------------------------------------------------*/
/* nth-root stacking n == 2                                                   */
/*----------------------------------------------------------------------------*/
		else if (hks->nroot == 2) {  
			for (irec = 0; irec < hks->nrfrec; irec++) {
				rfrec[irec].LoadData();
				printf("Processing rfrec[%d/%d]...\n", irec, hks->nrfrec);

				for (izr = 0; izr < hks->nzr; izr++) {
					rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array

					rfrec[irec].PeakAmp(hks);

					//if (rfrec[irec].a_0p1s == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					//if (rfrec[irec].a_2p1s == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					//if (rfrec[irec].a_1p2s == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					
					tmpx1 = SIGN1(rfrec[irec].a_0p1s); 
					tmpx2 = SIGN1(rfrec[irec].a_2p1s); 
					tmpx3 = SIGN1(rfrec[irec].a_1p2s); 
	
					hks->a_0p1s[izr] += tmpx1 * sqrt(fabs(rfrec[irec].a_0p1s));
					hks->a_2p1s[izr] += tmpx2 * sqrt(fabs(rfrec[irec].a_2p1s));
					hks->a_1p2s[izr] += tmpx3 * sqrt(fabs(rfrec[irec].a_1p2s));
				}
				rfrec[irec].FreeData();
			}
			
			for (izr = 0; izr < hks->nzr; izr++) {
					
				hks->a_0p1s[izr] /= hks->nrfrec;
				hks->a_2p1s[izr] /= hks->nrfrec;
				hks->a_1p2s[izr] /= hks->nrfrec;
	
				tmpx1 = 1.0;
				tmpx2 = 1.0;
				tmpx3 = 1.0;
				for (i = 0; i < hks->nroot; i++)	{
					tmpx1 *= fabs(hks->a_0p1s[izr]);
					tmpx2 *= fabs(hks->a_2p1s[izr]);
					tmpx3 *= fabs(hks->a_1p2s[izr]);
				}	
				hks->a_0p1s[izr] =  SIGN1(hks->a_0p1s[izr]) * tmpx1;
				hks->a_2p1s[izr] =  SIGN1(hks->a_2p1s[izr]) * tmpx2;
				hks->a_1p2s[izr] = -SIGN1(hks->a_1p2s[izr]) * tmpx3;
			}
		}
/*----------------------------------------------------------------------------*/
/* nth-root stacking n == 4                                                   */
/*----------------------------------------------------------------------------*/
		else if (hks->nroot == 4) {  
			for (irec = 0; irec < hks->nrfrec; irec++) {
				rfrec[irec].LoadData();
				printf("Processing rfrec[%d/%d]...\n", irec, hks->nrfrec);

				for (izr = 0; izr < hks->nzr; izr++) {
					//if (izr % (hks->nzr/10) == 0)	{
					//	printf("%d-root Stacking: %d, iz = %d, ir = %d\n", 
					//		hks->nroot, izr, izr / hks->nrv, izr % hks->nrv);
					//}
							
					rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array
					rfrec[irec].PeakAmp(hks);

					//if (rfrec[irec].a_0p1s == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					//if (rfrec[irec].a_2p1s == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					//if (rfrec[irec].a_1p2s == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					
					tmpx1 = SIGN1(rfrec[irec].a_0p1s); 
					tmpx2 = SIGN1(rfrec[irec].a_2p1s); 
					tmpx3 = SIGN1(rfrec[irec].a_1p2s); 
	
					hks->a_0p1s[izr] += tmpx1 * sqrt(sqrt(fabs(rfrec[irec].a_0p1s)));
					hks->a_2p1s[izr] += tmpx2 * sqrt(sqrt(fabs(rfrec[irec].a_2p1s)));
					hks->a_1p2s[izr] += tmpx3 * sqrt(sqrt(fabs(rfrec[irec].a_1p2s)));
				}
				rfrec[irec].FreeData();
			}
			
			for (izr = 0; izr < hks->nzr; izr++) {
					
				hks->a_0p1s[izr] /= hks->nrfrec;
				hks->a_2p1s[izr] /= hks->nrfrec;
				hks->a_1p2s[izr] /= hks->nrfrec;
	
				tmpx1 = 1.0;
				tmpx2 = 1.0;
				tmpx3 = 1.0;
				for (i = 0; i < hks->nroot; i++)	{
					tmpx1 *= fabs(hks->a_0p1s[izr]);
					tmpx2 *= fabs(hks->a_2p1s[izr]);
					tmpx3 *= fabs(hks->a_1p2s[izr]);
				}	
				hks->a_0p1s[izr] =  SIGN1(hks->a_0p1s[izr]) * tmpx1;
				hks->a_2p1s[izr] =  SIGN1(hks->a_2p1s[izr]) * tmpx2;
				hks->a_1p2s[izr] = -SIGN1(hks->a_1p2s[izr]) * tmpx3;
			}
		}
/*----------------------------------------------------------------------------*/
/* nth-root stacking                                                          */
/*----------------------------------------------------------------------------*/
		else {  
			for (irec = 0; irec < hks->nrfrec; irec++) {
				rfrec[irec].LoadData();
				printf("Processing rfrec[%d/%d]...\n", irec, hks->nrfrec);

				for (izr = 0; izr < hks->nzr; izr++) {
					rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array

					rfrec[irec].PeakAmp(hks);

					//if (rfrec[irec].a_0p1s == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					//if (rfrec[irec].a_2p1s == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					//if (rfrec[irec].a_1p2s == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					
					tmpx1 = SIGN1(rfrec[irec].a_0p1s); 
					tmpx2 = SIGN1(rfrec[irec].a_2p1s); 
					tmpx3 = SIGN1(rfrec[irec].a_1p2s); 
	
					hks->a_0p1s[izr] += tmpx1 * pow(fabs(rfrec[irec].a_0p1s),  hks->nrootpw);
					hks->a_2p1s[izr] += tmpx2 * pow(fabs(rfrec[irec].a_2p1s),  hks->nrootpw);
					hks->a_1p2s[izr] += tmpx3 * pow(fabs(rfrec[irec].a_1p2s),  hks->nrootpw);
				}
				rfrec[irec].FreeData();
			}
			
			for (izr = 0; izr < hks->nzr; izr++) {
					
				hks->a_0p1s[izr] /= hks->nrfrec;
				hks->a_2p1s[izr] /= hks->nrfrec;
				hks->a_1p2s[izr] /= hks->nrfrec;
	
				tmpx1 = 1.0;
				tmpx2 = 1.0;
				tmpx3 = 1.0;
				for (i = 0; i < hks->nroot; i++)	{
					tmpx1 *= fabs(hks->a_0p1s[izr]);
					tmpx2 *= fabs(hks->a_2p1s[izr]);
					tmpx3 *= fabs(hks->a_1p2s[izr]);
				}	
				hks->a_0p1s[izr] =  SIGN1(hks->a_0p1s[izr]) * tmpx1;
				hks->a_2p1s[izr] =  SIGN1(hks->a_2p1s[izr]) * tmpx2;
				hks->a_1p2s[izr] = -SIGN1(hks->a_1p2s[izr]) * tmpx3;
			}
		}
	}
		
/*----------------------------------------------------------------------------*/
/* From here, fork bootstrap and non-bootstrap
	bootstrap = btstrp.flag                                                    */
/*----------------------------------------------------------------------------*/
	else {

/*----------------------------------------------------------------------------*/
/* linear stacking                                                            */
/*----------------------------------------------------------------------------*/
		if (hks->nroot == 1) {  
			for (izr = 0; izr < hks->nzr; izr++) {
				if (izr % (hks->nzr/10) == 0)	{
					printf("%d-root Stacking: %d, iz = %d, ir = %d\n", 
						hks->nroot, izr, izr / hks->nrv, izr % hks->nrv);
				}		
				sa_0p1s = 0.0;
				sa_2p1s = 0.0;
				sa_1p2s = 0.0;
				for (irec = 0; irec < hks->nrfrec; irec++) {

					if (rfrec[irec].n_resmpl != 0) {

						rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array
						//rfrec[irec].LoadData(rf_data, hks);
						rfrec[irec].PeakAmp(hks);

						tmpx1 = rfrec[irec].a_0p1s * rfrec[irec].n_resmpl; 
						tmpx2 = rfrec[irec].a_2p1s * rfrec[irec].n_resmpl; 
						tmpx3 = rfrec[irec].a_1p2s * rfrec[irec].n_resmpl; 

						if (tmpx1 == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
						if (tmpx2 == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
						if (tmpx3 == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);

						sa_0p1s += tmpx1;
						sa_2p1s += tmpx2;
						sa_1p2s += tmpx3;
					}	
				}	
				sa_0p1s /= hks->nrfrec;
				sa_2p1s /= hks->nrfrec;
				sa_1p2s /= hks->nrfrec;
		
				hks->a_0p1s[izr] = sa_0p1s;
				hks->a_2p1s[izr] = sa_2p1s;
				hks->a_1p2s[izr] = -sa_1p2s;
			}	
		}
/*----------------------------------------------------------------------------*/
/* nth-root stacking n == 2                                                   */
/*----------------------------------------------------------------------------*/
		else if (hks->nroot == 2) {  
			for (izr = 0; izr < hks->nzr; izr++) {
				if (izr % (hks->nzr/10) == 0)	{
					printf("%d-root Stacking: %d, iz = %d, ir = %d\n", 
						hks->nroot, izr, izr / hks->nrv, izr % hks->nrv);
				}		
				sa_0p1s = 0.0;
				sa_2p1s = 0.0;
				sa_1p2s = 0.0;
				for (irec = 0; irec < hks->nrfrec; irec++) {

					if (rfrec[irec].n_resmpl != 0) {

						rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array
						//rfrec[irec].LoadData(rf_data, hks);
						rfrec[irec].PeakAmp(hks);

						tmpx1 = rfrec[irec].a_0p1s * rfrec[irec].n_resmpl; 
						tmpx2 = rfrec[irec].a_2p1s * rfrec[irec].n_resmpl; 
						tmpx3 = rfrec[irec].a_1p2s * rfrec[irec].n_resmpl; 

						if (tmpx1 == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
						if (tmpx2 == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
						if (tmpx3 == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					
						tmpx1 = SIGN1(tmpx1) * sqrt(fabs(tmpx1));
						tmpx2 = SIGN1(tmpx2) * sqrt(fabs(tmpx2));
						tmpx3 = SIGN1(tmpx3) * sqrt(fabs(tmpx3));

						sa_0p1s += tmpx1;
						sa_2p1s += tmpx2;
						sa_1p2s += tmpx3;
					}	
				}	
				sa_0p1s /= hks->nrfrec;
				sa_2p1s /= hks->nrfrec;
				sa_1p2s /= hks->nrfrec;
	
				tmpx1 = 1.0;
				tmpx2 = 1.0;
				tmpx3 = 1.0;
				for (i = 0; i < hks->nroot; i++)	{
					tmpx1 *= fabs(sa_0p1s);
					tmpx2 *= fabs(sa_2p1s);
					tmpx3 *= fabs(sa_1p2s);
				}	
				hks->a_0p1s[izr] =  SIGN1(sa_0p1s) * tmpx1;
				hks->a_2p1s[izr] =  SIGN1(sa_2p1s) * tmpx2;
				hks->a_1p2s[izr] = -SIGN1(sa_1p2s) * tmpx3;
			}
		}
/*----------------------------------------------------------------------------*/
/* nth-root stacking n == 4                                                   */
/*----------------------------------------------------------------------------*/
		else if (hks->nroot == 4) {  
			for (izr = 0; izr < hks->nzr; izr++) {
				if (izr % (hks->nzr/10) == 0)	{
					printf("%d-root Stacking: %d, iz = %d, ir = %d\n", 
						hks->nroot, izr, izr / hks->nrv, izr % hks->nrv);
				}		
				sa_0p1s = 0.0;
				sa_2p1s = 0.0;
				sa_1p2s = 0.0;
				for (irec = 0; irec < hks->nrfrec; irec++) {

					if (rfrec[irec].n_resmpl != 0) {

						rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array
						//rfrec[irec].LoadData(rf_data, hks);
						rfrec[irec].PeakAmp(hks);

						tmpx1 = rfrec[irec].a_0p1s * rfrec[irec].n_resmpl; 
						tmpx2 = rfrec[irec].a_2p1s * rfrec[irec].n_resmpl; 
						tmpx3 = rfrec[irec].a_1p2s * rfrec[irec].n_resmpl; 

						if (tmpx1 == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
						if (tmpx2 == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
						if (tmpx3 == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					
						tmpx1 = SIGN1(tmpx1) * sqrt(sqrt(fabs(tmpx1)));
						tmpx2 = SIGN1(tmpx2) * sqrt(sqrt(fabs(tmpx2)));
						tmpx3 = SIGN1(tmpx3) * sqrt(sqrt(fabs(tmpx3)));

						sa_0p1s += tmpx1;
						sa_2p1s += tmpx2;
						sa_1p2s += tmpx3;
					}	
				}	
				sa_0p1s /= hks->nrfrec;
				sa_2p1s /= hks->nrfrec;
				sa_1p2s /= hks->nrfrec;
	
				tmpx1 = 1.0;
				tmpx2 = 1.0;
				tmpx3 = 1.0;
				for (i = 0; i < hks->nroot; i++)	{
					tmpx1 *= fabs(sa_0p1s);
					tmpx2 *= fabs(sa_2p1s);
					tmpx3 *= fabs(sa_1p2s);
				}	
				hks->a_0p1s[izr] =  SIGN1(sa_0p1s) * tmpx1;
				hks->a_2p1s[izr] =  SIGN1(sa_2p1s) * tmpx2;
				hks->a_1p2s[izr] = -SIGN1(sa_1p2s) * tmpx3;
			}
		}
/*----------------------------------------------------------------------------*/
/* nth-root stacking                                                          */
/*----------------------------------------------------------------------------*/
		else {  
			for (izr = 0; izr < hks->nzr; izr++) {
				if (izr % (hks->nzr/10) == 0)	{
					printf("%d-root Stacking: %d, iz = %d, ir = %d\n", 
						hks->nroot, izr, izr / hks->nrv, izr % hks->nrv);
				}		
				sa_0p1s = 0.0;
				sa_2p1s = 0.0;
				sa_1p2s = 0.0;
				for (irec = 0; irec < hks->nrfrec; irec++) {

					if (rfrec[irec].n_resmpl != 0) {

						rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array
						//rfrec[irec].LoadData(rf_data, hks);
						rfrec[irec].PeakAmp(hks);

						tmpx1 = rfrec[irec].a_0p1s * rfrec[irec].n_resmpl; 
						tmpx2 = rfrec[irec].a_2p1s * rfrec[irec].n_resmpl; 
						tmpx3 = rfrec[irec].a_1p2s * rfrec[irec].n_resmpl; 

						if (tmpx1 == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
						if (tmpx2 == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
						if (tmpx3 == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					
						tmpx1 = SIGN1(tmpx1) * pow(fabs(tmpx1),  hks->nrootpw);
						tmpx2 = SIGN1(tmpx2) * pow(fabs(tmpx2), hks->nrootpw);
						tmpx3 = SIGN1(tmpx3) * pow(fabs(tmpx3), hks->nrootpw);

						sa_0p1s += tmpx1;
						sa_2p1s += tmpx2;
						sa_1p2s += tmpx3;
					}	
				}	
				sa_0p1s /= hks->nrfrec;
				sa_2p1s /= hks->nrfrec;
				sa_1p2s /= hks->nrfrec;
	
				tmpx1  = 1.0;
				tmpx2 = 1.0;
				tmpx3 = 1.0;
				for (i = 0; i < hks->nroot; i++)	{
					tmpx1 *= fabs(sa_0p1s);
					tmpx2 *= fabs(sa_2p1s);
					tmpx3 *= fabs(sa_1p2s);
				}	
				hks->a_0p1s[izr] =  SIGN1(sa_0p1s) * tmpx1;
				hks->a_2p1s[izr] =  SIGN1(sa_2p1s) * tmpx2;
				hks->a_1p2s[izr] = -SIGN1(sa_1p2s) * tmpx3;
			}
		}
			
	}			
	return;
}	


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* _acc_()                                                                    */
/*----------------------------------------------------------------------------*/
double _acc_(int n, float* data1, float* data2)
{
	int  i;
	double  sum;

	sum = 0.;
	for (i = 0; i < n; i++) {
		sum += (data1[i] * data2[i]);
	}
	return(sum);
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* wgt_xc() combination of wgtrv_xc and wgtmh_xc                              */
/*----------------------------------------------------------------------------*/
int wgt_xc(HKS_INFO *hks)
{
	int iz, ir, it, itk;
	int i, k; 
	int *imh_max0, imh_max;
	int lfmh, rhmh;
	int nmh2, nmh;
	int *irv_max0, irv_max;
	int bmrv, tprv;
	int nrv2, nrv;
	float	mh_max, rv_max;
	float	*ar_0p1ss, *ar_0p1s, *ar_2p1s, *ar_1p2s;
	float *xc_rv, *xc_mh;
	double cc, c11, c12, c13, c22, c23, c33, rc12, rc13, wgt;
	
/*----------------------------------------------------------------------------*/
/* This part is for xcorr as a function of rv                                 */
/*----------------------------------------------------------------------------*/
	nmh2 = (int)rint(_ZCC_L_ * 0.5 / hks->mhdt);
	if (nmh2 < 1 ) {
		fprintf(stderr, "Error in wgt_xc: _cc_wgt_: nmh2 = %3d, %8.5f\n", nmh2, hks->mhdt);
		return -1;
	}	
   
	imh_max0 = new int[hks->nrv]; 
	 
/*----------------------------------------------------------------------------*/
/* for each vp/vs ratio, find the peak in the PS mode                         */
/*----------------------------------------------------------------------------*/
	for (ir = 0; ir < hks->nrv; ir++) {
		mh_max = -MAXFLOAT;
		for (iz = 0; iz < hks->nmh; iz++) {
			it = iz * hks->nrv + ir;
			if (mh_max < hks->a_0p1s[it]) {
				mh_max = hks->a_0p1s[it];
				imh_max0[ir] = iz;
			}
		}
	}

	xc_rv = new float[hks->nrv];	
	ar_0p1ss = new float[hks->nmh];
	ar_0p1s = new float[hks->nmh];
	ar_2p1s = new float[hks->nmh];
	ar_1p2s = new float[hks->nmh];
	
//	rv = hks->rvlb; 
//	hks->xc_max = 0.0;
	for (ir = 0; ir < hks->nrv; ir++) {

		for (iz = 0; iz < hks->nmh; iz++) {
			it = iz * hks->nrv + ir;
   	   ar_0p1s[iz] = hks->a_0p1s[it];
      	ar_2p1s[iz] = hks->a_2p1s[it];
			ar_1p2s[iz] = hks->a_1p2s[it];
			ar_0p1ss[iz] = 0.0;
		}
		
/*----------------------------------------------------------------------------*/
/* imh_max average index over max_mh at nearby rv                             */
/*----------------------------------------------------------------------------*/
		if (ir == 0 || ir == hks->nrv - 1) {
			imh_max = imh_max0[ir];
		}
		else if (ir == 1 || ir == hks->nrv - 2) {
			imh_max = (int)(0.5 * imh_max0[ir] + 0.25 * (imh_max0[ir-1] + imh_max0[ir+1]));
		}
		else{
			imh_max = (int)(0.4 * imh_max0[ir] + 0.2 * (imh_max0[ir-1] + imh_max0[ir+1]) + 
                      0.1 * (imh_max0[ir-2] + imh_max0[ir+2]));
		}

/*----------------------------------------------------------------------------*/
/* set up smoothing range for mh                                              */
/*----------------------------------------------------------------------------*/
		lfmh = MAX(imh_max - nmh2, 0);
		rhmh = MIN(imh_max + nmh2, hks->nmh-1);
		nmh = rhmh - lfmh + 1;
		wgt = (double) nmh / (double)(nmh2 * 2 + 1);
		
		for (i = lfmh; i <= rhmh; i++) {
			//printf("_cc_wgt_: ir= %3d %7.4f i_max= %3d %6.1f i1= %3d i2= %3d zN= %3d ", 
			//		  	 ir, r, i_max, i_max * ds->z_inc + ds->z_min, i1, i2, zN);
			it = (i - imh_max) * _CD_FACTOR_ + imh_max;
			if (it >= 0 && it < hks->nmh) {
				for (k = 0; k < _CD_FACTOR_; k++) {
					itk = it + k;
					if (itk >= hks->nmh) break;	
					ar_0p1ss[i] += ar_0p1s[itk];
				}
//			a_0p1ss[ii] /= (float)k; 
			}
		//printf("ii = %3d it = %3d %10.5f\n", ii, it, a_0p1ss[ii]);
		}
		c11 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_0p1ss[lfmh]);
		c22 = _acc_(nmh, &ar_2p1s[lfmh],  &ar_2p1s[lfmh]);
		c33 = _acc_(nmh, &ar_1p2s[lfmh],  &ar_1p2s[lfmh]);
		c12 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_2p1s[lfmh]);
		c13 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_1p2s[lfmh]);
		c23 = _acc_(nmh, &ar_2p1s[lfmh],  &ar_1p2s[lfmh]);	

		rc12 = c12/sqrt(c11 * c22); 
		rc13 = c13/sqrt(c11 * c33);
		switch(hks->xc_flag) {
			case 7:
		    	xc_rv[ir] = (rc12 + rc13) * 0.5 * wgt;
		    	break;
			case 8:
	   	 	xc_rv[ir] = rc12 * wgt;
	    		break;
			case 9:
	    		xc_rv[ir] = rc13 * wgt;
	    		break;
		}
		
//		if (hks->xc[ir] > hks->xc_max) {
//			hks->xc_max = hks->xc[ir];
//			hks->rv_xcmax = rv;
//		}
//		rv += hks->rvdt;
	
	} /* ratio loop */

/*----------------------------------------------------------------------------*/
/* This part is for xcorr as a function of mh                                 */
/*----------------------------------------------------------------------------*/
	nrv2 = (int)rint(_RCC_L_ * 0.5 / hks->rvdt);
	if (nrv2 < 1 ) {
		fprintf(stderr, "Error in wgt_xc: _cc_wgt_: nrv2 = %3d, %8.5f\n", nrv2, hks->rvdt);
		return -1;
	}	
   
	irv_max0 = new int[hks->nmh];  

/*----------------------------------------------------------------------------*/
/* for each moho depth, find the peak in the PS mode                          */
/*----------------------------------------------------------------------------*/
	for (iz = 0; iz < hks->nmh; iz++) {
		rv_max = -MAXFLOAT;
		for (ir = 0; ir < hks->nrv; ir++) {
			it = iz * hks->nrv + ir;
			if (rv_max < hks->a_0p1s[it]) {
				rv_max = hks->a_0p1s[it];
				irv_max0[iz] = ir;
			}
		}
	}
  
	xc_mh = new float[hks->nmh];
	delete [] ar_0p1ss;
	ar_0p1ss = new float[hks->nrv];
	delete [] ar_0p1s;
	ar_0p1s = new float[hks->nrv];
	delete [] ar_2p1s;
	ar_2p1s = new float[hks->nrv];
	delete [] ar_1p2s;
	ar_1p2s = new float[hks->nrv];
		
//	mh = hks->mhlb; 
//	hks->xc_max = 0.0;
	for (iz = 0; iz < hks->nmh; iz++) {

		for (ir = 0; ir < hks->nrv; ir++) {
			it = iz * hks->nrv + ir;
   	   ar_0p1s[ir] = hks->a_0p1s[it];
      	ar_2p1s[ir] = hks->a_2p1s[it];
			ar_1p2s[ir] = hks->a_1p2s[it];
			ar_0p1ss[ir] = 0.0;
		}
		
/*----------------------------------------------------------------------------*/
/* irv_max average index over max_rv at nearby mh                             */
/*----------------------------------------------------------------------------*/
		if (iz == 0 || iz == hks->nmh - 1) {
			irv_max = irv_max0[iz];
		}
		else if (iz == 1 || iz == hks->nmh - 2) {
			irv_max = (int)(0.5 * irv_max0[iz] + 0.25 * (irv_max0[iz-1] + irv_max0[iz+1]));
		}
		else{
			irv_max = (int)(0.4 * irv_max0[iz] + 0.2 * (irv_max0[iz-1] + irv_max0[iz+1]) + 
                      0.1 * (irv_max0[iz-2] + irv_max0[iz+2]));
		}

/*----------------------------------------------------------------------------*/
/* set up smoothing range for rv                                              */
/*----------------------------------------------------------------------------*/
		bmrv = MAX(irv_max - nrv2, 0);
		tprv = MIN(irv_max + nrv2, hks->nrv-1);
		nrv = tprv - bmrv + 1;
		wgt = (double) nrv / (double)(nrv2 * 2 + 1);
		
		for (i = bmrv; i <= tprv; i++) {
			//printf("_cc_wgt_: ir= %3d %7.4f i_max= %3d %6.1f i1= %3d i2= %3d zN= %3d ", 
			//  	 ir, r, i_max, i_max * ds->z_inc + ds->z_min, i1, i2, zN);
			it = (i - irv_max) * _CD_FACTOR_ + irv_max;
			if (it >= 0 && it < hks->nrv) {
				for (k = 0; k < _CD_FACTOR_; k++) {
					itk = it + k;
					if (itk >= hks->nrv) break;	
					ar_0p1ss[i] += ar_0p1s[itk];
				}
//			a_0p1ss[ii] /= (float)k; 
			}
		//printf("ii = %3d it = %3d %10.5f\n", ii, it, a_0p1ss[ii]);
		}
		c11 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_0p1ss[bmrv]);
		c22 = _acc_(nrv, &ar_2p1s[bmrv],  &ar_2p1s[bmrv]);
		c33 = _acc_(nrv, &ar_1p2s[bmrv],  &ar_1p2s[bmrv]);
		c12 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_2p1s[bmrv]);
		c13 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_1p2s[bmrv]);
		c23 = _acc_(nrv, &ar_2p1s[bmrv],  &ar_1p2s[bmrv]);	

		rc12 = c12/sqrt(c11 * c22); 
		rc13 = c13/sqrt(c11 * c33);
		switch(hks->xc_flag) {
			case 7:
		    	xc_mh[iz] = (rc12 + rc13) * 0.5 * wgt;
		    	break;
			case 8:
	   	 	xc_mh[iz] = rc12 * wgt;
	    		break;
			case 9:
				xc_mh[iz] = rc13 * wgt;
				break;
		}
		
//		if (hks->xc[iz] > hks->xc_max) {
//			hks->xc_max = hks->xc[iz];
//			hks->rv_xcmax = mh;
//		}
//		mh += hks->mhdt;
	
	} /* moho loop */
		
		
/*----------------------------------------------------------------------------*/
/* Combination of Xcorr(rv) * Xcorr(mh)                                       */
/*----------------------------------------------------------------------------*/
	for (ir = 0; ir < hks->nrv; ir++) {
		if (ir == 0 || ir == hks->nrv - 1) {
			cc = xc_rv[ir] + 1.;
		}
		else if (ir == 1 || ir == hks->nrv - 2) {
			cc = 0.50 * xc_rv[ir] + 0.25 * (xc_rv[ir-1] + xc_rv[ir+1]) + 1.;
		}
		else{
      	cc = 0.40 * xc_rv[ir] + 0.20 * (xc_rv[ir-1] + xc_rv[ir+1]) + 
				  0.10 * (xc_rv[ir-2] + xc_rv[ir+2]) + 1.;
		}
//		ccw[ir] = (float) cc;
		for (iz = 0; iz < hks->nmh; iz++) {
			it = iz * hks->nrv + ir;
			hks->a_0p1s[it] *= cc;
			hks->a_2p1s[it] *= cc;
			hks->a_1p2s[it] *= cc;
		}
	}
	
	for (iz = 0; iz < hks->nmh; iz++) {
		if (iz == 0 || iz == hks->nmh - 1) {
			cc = xc_mh[iz] + 1.;
		}
		else if (iz == 1 || iz == hks->nmh - 2) {
			cc = 0.50 * xc_mh[iz] + 0.25 * (xc_mh[iz-1] + xc_mh[iz+1]) + 1.;
		}
		else{
      	cc = 0.40 * xc_mh[iz] + 0.20 * (xc_mh[iz-1] + xc_mh[iz+1]) + 
				  0.10 *	(xc_mh[iz-2] + xc_mh[iz+2]) + 1.;
		}
//		ccw[ir] = (float) cc;
		for (ir = 0; ir < hks->nrv; ir++) {
			it = iz * hks->nrv + ir;
			hks->a_0p1s[it] *= cc;
			hks->a_2p1s[it] *= cc;
			hks->a_1p2s[it] *= cc;
		}
	}
	
	delete [] ar_0p1ss;
	ar_0p1ss = NULL;
	delete [] ar_0p1s;
	ar_0p1s = NULL;
	delete [] ar_2p1s;
	ar_2p1s = NULL;
	delete [] ar_1p2s;
	ar_1p2s = NULL;
	delete [] xc_rv;
	xc_rv = NULL;
	delete [] xc_mh;
	xc_mh = NULL;
	delete [] imh_max0;
	imh_max0 = NULL;
	delete [] irv_max0;  
	irv_max0 = NULL;
	
	return 0;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* wgtrv_xc() xcorr is a function of rv                                       */
/*----------------------------------------------------------------------------*/
int wgtrv_xc(HKS_INFO *hks)
{
	int iz, ir, it, itk;
	int i, k; 
	int *imh_max0, imh_max;
	int lfmh, rhmh;
	int nmh2, nmh;
	float	mh_max;
	float rv;
	float	*ar_0p1ss, *ar_0p1s, *ar_2p1s, *ar_1p2s;
	double cc, c11, c12, c13, c22, c23, c33, rc12, rc13, wgt;
	
	nmh2 = (int)rint(_ZCC_L_ * 0.5 / hks->mhdt);
	if (nmh2 < 1 ) {
		fprintf(stderr, "Error in wgtrv_xc: _cc_wgt_: nmh2 = %3d %8.5f\n", nmh2, hks->mhdt);
		return -1;
	}	
   
	imh_max0 = new int[hks->nrv];  

/*----------------------------------------------------------------------------*/
/* for each vp/vs ratio, find the peak in the PS mode                         */
/*----------------------------------------------------------------------------*/
	for (ir = 0; ir < hks->nrv; ir++) {
		mh_max = -MAXFLOAT;
		for (iz = 0; iz < hks->nmh; iz++) {
			it = iz * hks->nrv + ir;
			if (mh_max < hks->a_0p1s[it]) {
				mh_max = hks->a_0p1s[it];
				imh_max0[ir] = iz;
			}
		}
	}
  
	if (hks->xc != NULL) { 
		delete [] hks->xc;
		hks->xc = NULL;
	}	
	hks->xc = new float[hks->nrv];
	
	ar_0p1ss = new float[hks->nmh];
	ar_0p1s = new float[hks->nmh];
	ar_2p1s = new float[hks->nmh];
	ar_1p2s = new float[hks->nmh];
	
	rv = hks->rvlb; 
	hks->xc_max = 0.0;
	for (ir = 0; ir < hks->nrv; ir++) {

		for (iz = 0; iz < hks->nmh; iz++) {
			it = iz * hks->nrv + ir;
   	   ar_0p1s[iz] = hks->a_0p1s[it];
      	ar_2p1s[iz] = hks->a_2p1s[it];
			ar_1p2s[iz] = hks->a_1p2s[it];
			ar_0p1ss[iz] = 0.0;
		}
		
/*----------------------------------------------------------------------------*/
/* imh_max average index over max_mh at nearby rv                             */
/*----------------------------------------------------------------------------*/
		if (ir == 0 || ir == hks->nrv - 1) {
			imh_max = imh_max0[ir];
		}
		else if (ir == 1 || ir == hks->nrv - 2) {
			imh_max = (int)(0.5 * imh_max0[ir] + 0.25 * (imh_max0[ir-1] + imh_max0[ir+1]));
		}
		else{
			imh_max = (int)(0.4 * imh_max0[ir] + 0.2 * (imh_max0[ir-1] + imh_max0[ir+1]) + 
                      0.1 * (imh_max0[ir-2] + imh_max0[ir+2]));
		}

/*----------------------------------------------------------------------------*/
/* set up smoothing range for mh                                              */
/*----------------------------------------------------------------------------*/
		lfmh = MAX(imh_max - nmh2, 0);
		rhmh = MIN(imh_max + nmh2, hks->nmh-1);
		nmh = rhmh - lfmh + 1;
		wgt = (double) nmh / (double)(nmh2 * 2 + 1);
		
		for (i = lfmh; i <= rhmh; i++) {
			//printf("_cc_wgt_: ir= %3d %7.4f i_max= %3d %6.1f i1= %3d i2= %3d zN= %3d ", 
			//		  	 ir, r, i_max, i_max * ds->z_inc + ds->z_min, i1, i2, zN);
			it = (i - imh_max) * _CD_FACTOR_ + imh_max;
			if (it >= 0 && it < hks->nmh) {
				for (k = 0; k < _CD_FACTOR_; k++) {
					itk = it + k;
					if (itk >= hks->nmh) break;	
					ar_0p1ss[i] += ar_0p1s[itk];
				}
//			a_0p1ss[ii] /= (float)k; 
			}
		//printf("ii = %3d it = %3d %10.5f\n", ii, it, a_0p1ss[ii]);
		}
		c11 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_0p1ss[lfmh]);
		c22 = _acc_(nmh, &ar_2p1s[lfmh],  &ar_2p1s[lfmh]);
		c33 = _acc_(nmh, &ar_1p2s[lfmh],  &ar_1p2s[lfmh]);
		c12 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_2p1s[lfmh]);
		c13 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_1p2s[lfmh]);
		c23 = _acc_(nmh, &ar_2p1s[lfmh],  &ar_1p2s[lfmh]);	

		rc12 = c12/sqrt(c11 * c22); 
		rc13 = c13/sqrt(c11 * c33);
		switch(hks->xc_flag) {
			case 1:
				hks->xc[ir] = (rc12 + rc13) * 0.5 * wgt;
				break;
			case 2:
				hks->xc[ir] = rc12 * wgt;
				break;
			case 3:
				hks->xc[ir] = rc13 * wgt;
				break;
		}
		
		if (hks->xc[ir] > hks->xc_max) {
			hks->xc_max = hks->xc[ir];
			hks->rv_xcmax = rv;
		}
		rv += hks->rvdt;
	
	} /* ratio loop */
		
	for (ir = 0; ir < hks->nrv; ir++) {
		if (ir == 0 || ir == hks->nrv - 1) {
			cc = hks->xc[ir] + 1.;
		}
		else if (ir == 1 || ir == hks->nrv - 2) {
			cc = 0.50 * hks->xc[ir] + 0.25 * (hks->xc[ir-1] + hks->xc[ir+1]) + 1.;
		}
		else{
      	cc = 0.40 * hks->xc[ir] + 0.20 * (hks->xc[ir-1] + hks->xc[ir+1]) + 
				  0.10 * (hks->xc[ir-2] + hks->xc[ir+2]) + 1.;
		}
//    ccw[ir] = (float) cc;
		for (iz = 0; iz < hks->nmh; iz++) {
			it = iz * hks->nrv + ir;
			hks->a_0p1s[it] *= cc;
			hks->a_2p1s[it] *= cc;
			hks->a_1p2s[it] *= cc;
		}
	}
	
	delete [] ar_0p1ss;
	ar_0p1ss = NULL;
	delete [] ar_0p1s;
	ar_0p1s = NULL;
	delete [] ar_2p1s;
	ar_2p1s = NULL;
	delete [] ar_1p2s;
	ar_1p2s = NULL;
	delete [] imh_max0;  
	imh_max0 = NULL;
	
	return 0;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* wgtmh_xc() xcorr is a function of mh                                       */
/*----------------------------------------------------------------------------*/
int wgtmh_xc(HKS_INFO *hks)
{
	int iz, ir, it, itk;
	int i, k; 
	int *irv_max0, irv_max;
	int bmrv, tprv;
	int nrv2, nrv;
	float	rv_max;
	float mh;
	float	*ar_0p1ss, *ar_0p1s, *ar_2p1s, *ar_1p2s;
	double	cc, c11, c12, c13, c22, c23, c33, rc12, rc13, wgt;
	
	nrv2 = (int)rint(_RCC_L_ * 0.5 / hks->rvdt);
	if (nrv2 < 1 ) {
		fprintf(stderr, "Error in wgtmh_xc: _cc_wgt_: nrv2 = %3d, %8.5f\n", nrv2, hks->rvdt);
		return -1;
	}	
   
	irv_max0 = new int[hks->nmh];  

/*----------------------------------------------------------------------------*/
/* for each moho depth, find the peak in the PS mode                          */
/*----------------------------------------------------------------------------*/
	for (iz = 0; iz < hks->nmh; iz++) {
		rv_max = -MAXFLOAT;
		for (ir = 0; ir < hks->nrv; ir++) {
			it = iz * hks->nrv + ir;
			if (rv_max < hks->a_0p1s[it]) {
				rv_max = hks->a_0p1s[it];
				irv_max0[iz] = ir;
			}
		}
	}

  
	if (hks->xc != NULL) { 
		delete [] hks->xc;
		hks->xc = NULL;
	}	
	hks->xc = new float[hks->nmh];
	
	ar_0p1ss = new float[hks->nrv];
	ar_0p1s = new float[hks->nrv];
	ar_2p1s = new float[hks->nrv];
	ar_1p2s = new float[hks->nrv];
	
	
	mh = hks->mhlb; 
	hks->xc_max = 0.0;
	for (iz = 0; iz < hks->nmh; iz++) {

		for (ir = 0; ir < hks->nrv; ir++) {
			it = iz * hks->nrv + ir;
   	   ar_0p1s[ir] = hks->a_0p1s[it];
      	ar_2p1s[ir] = hks->a_2p1s[it];
			ar_1p2s[ir] = hks->a_1p2s[it];
			ar_0p1ss[ir] = 0.0;
		}
		
/*----------------------------------------------------------------------------*/
/* irv_max average index over max_rv at nearby mh                             */
/*----------------------------------------------------------------------------*/
		if (iz == 0 || iz == hks->nmh - 1) {
			irv_max = irv_max0[iz];
		}
		else if (iz == 1 || iz == hks->nmh - 2) {
			irv_max = (int)(0.5 * irv_max0[iz] + 0.25 * (irv_max0[iz-1] + irv_max0[iz+1]));
		}
		else{
			irv_max = (int)(0.4 * irv_max0[iz] + 0.2 * (irv_max0[iz-1] + irv_max0[iz+1]) + 
                      0.1 * (irv_max0[iz-2] + irv_max0[iz+2]));
		}

/*----------------------------------------------------------------------------*/
/* set up smoothing range for rv                                              */
/*----------------------------------------------------------------------------*/
		bmrv = MAX(irv_max - nrv2, 0);
		tprv = MIN(irv_max + nrv2, hks->nrv-1);
		nrv = tprv - bmrv + 1;
		wgt = (double) nrv / (double)(nrv2 * 2 + 1);
		
		for (i = bmrv; i <= tprv; i++) {
			//printf("_cc_wgt_: ir= %3d %7.4f i_max= %3d %6.1f i1= %3d i2= %3d zN= %3d ", 
			//		  	 ir, r, i_max, i_max * ds->z_inc + ds->z_min, i1, i2, zN);
			it = (i - irv_max) * _CD_FACTOR_ + irv_max;
			if (it >= 0 && it < hks->nrv) {
				for (k = 0; k < _CD_FACTOR_; k++) {
					itk = it + k;
					if (itk >= hks->nrv) break;	
					ar_0p1ss[i] += ar_0p1s[itk];
				}
//			a_0p1ss[ii] /= (float)k; 
			}
		//printf("ii = %3d it = %3d %10.5f\n", ii, it, a_0p1ss[ii]);
		}
		c11 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_0p1ss[bmrv]);
		c22 = _acc_(nrv, &ar_2p1s[bmrv],  &ar_2p1s[bmrv]);
		c33 = _acc_(nrv, &ar_1p2s[bmrv],  &ar_1p2s[bmrv]);
		c12 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_2p1s[bmrv]);
		c13 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_1p2s[bmrv]);
		c23 = _acc_(nrv, &ar_2p1s[bmrv],  &ar_1p2s[bmrv]);	

		rc12 = c12/sqrt(c11 * c22); 
		rc13 = c13/sqrt(c11 * c33);
		switch(hks->xc_flag) {
			case 4:
				//ds->cc[ir] = (c12/sqrt(c11 * c22) + c13/sqrt(c11 * c33) + c23/sqrt(c22 * c33))/3. * wgt; 
				hks->xc[iz] = (rc12 + rc13) * 0.5 * wgt;
				break;
			case 5:
				hks->xc[iz] = rc12 * wgt;
				break;
			case 6:
				hks->xc[iz] = rc13 * wgt;
				break;
		}
		
		if (hks->xc[iz] > hks->xc_max) {
			hks->xc_max = hks->xc[iz];
			hks->rv_xcmax = mh;
		}
		mh += hks->mhdt;
	
	} /* moho loop */
		
	for (iz = 0; iz < hks->nmh; iz++) {
		if (iz == 0 || iz == hks->nmh - 1) {
			cc = hks->xc[iz] + 1.;
		}
		else if (iz == 1 || iz == hks->nmh - 2) {
			cc = 0.50 * hks->xc[iz] + 0.25 * (hks->xc[iz-1] + hks->xc[iz+1]) + 1.;
		}
		else{
      	cc = 0.40 * hks->xc[iz] + 0.20 * (hks->xc[iz-1] + hks->xc[iz+1]) + 
				  0.10 *	(hks->xc[iz-2] + hks->xc[iz+2]) + 1.;
		}
//    ccw[ir] = (float) cc;
		for (ir = 0; ir < hks->nrv; ir++) {
			it = iz * hks->nrv + ir;
			hks->a_0p1s[it] *= cc;
			hks->a_2p1s[it] *= cc;
			hks->a_1p2s[it] *= cc;
		}
	}
	
	delete [] ar_0p1ss;
	ar_0p1ss = NULL;
	delete [] ar_0p1s;
	ar_0p1s = NULL;
	delete [] ar_2p1s;
	ar_2p1s = NULL;
	delete [] ar_1p2s;
	ar_1p2s = NULL;
	delete [] irv_max0;  
	irv_max0 = NULL;
	
	return 0;
}


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* statistics                                                                 
	!bootstrap for normal, bootstrap for bootstrap                             */
/*----------------------------------------------------------------------------*/
int statis(HKS_INFO *hks, const int& bootstrap)
{
	int i, j, iz, ir, izr; 
	int irv_max, imh_max, kzr;
	float mh_lb1, mh_ub1, mh_lb2, mh_ub2; 
	float rv_lb1, rv_ub1, rv_lb2, rv_ub2;
	float tmpx1, tmpx2, tmpx;
	double sum_mh, sum_rv, sum_mhrv, sum_mh2, sum_rv2;
	double mhtmp, rvtmp;
	int dist2, dist2_crt;
	int isOK;
	//float *tmpsfunc;
	//int *tmpkzr, *tmpindex;
	MAX_SFUNCTIONS exsf; 
	
	hks->sfunc_max = -MAXFLOAT;
	kzr = -MAXINT;
	irv_max = -MAXINT;
	imh_max = -MAXINT;
	for (izr = 0; izr < hks->nzr; izr++) {

		hks->a3ps[izr] = hks->w_0p1s * hks->a_0p1s[izr] 
							+ hks->w_2p1s * hks->a_2p1s[izr] 
							+ hks->w_1p2s * hks->a_1p2s[izr];
		
		if (hks->sfunc_max < hks->a3ps[izr]) {
			hks->sfunc_max = hks->a3ps[izr];			
			kzr = izr;
		}
	}
	if (kzr == -MAXINT) {
		fprintf(stderr, "Error in statis: NOT found the Max sfunc = %e\n", hks->sfunc_max);
		return -1;
	}
	else {	 		
		printf("The Max sfunc = %e\n", hks->sfunc_max);
	}

/*----------------------------------------------------------------------------*/
/* Normalize      
	Mean and variance (sqrt(sigma)) of S-function                              
	Find out the S-functions > 0.6 (normalized values)                         */
/*----------------------------------------------------------------------------*/
	hks->sfunc_mean = 0.0;
	for (izr = 0; izr < hks->nzr; izr++) {			
		hks->a3ps[izr] /= hks->sfunc_max; 
		hks->sfunc_mean += hks->a3ps[izr];
	}	
	hks->sfunc_max = 1.0;
	hks->sfunc_mean /= hks->nzr;

	hks->sfunc_sigma = 0.0;		
	hks->ntags = 0;
	for (izr = 0; izr < hks->nzr; izr++) {			
		hks->sfunc_sigma += (hks->a3ps[izr] - hks->sfunc_mean)*(hks->a3ps[izr] - hks->sfunc_mean);	
		//sfunc_store[ix].sigma_s += (hks->a3ps[i*hks->nmh + j] - hks->a3ps[irv_max*hks->nmh + imh_max]) 
		//				 * (hks->a3ps[i*hks->nmh + j] - hks->a3ps[irv_max*hks->nmh + imh_max]);
		hks->tag[izr] = 0;
		if (hks->a3ps[izr] > 0.6) {
			hks->tag[izr] = 1;
			hks->ntags++;
		}	
	}
	hks->sfunc_sigma /= (hks->nzr - 1);
	hks->sfunc_sigma = sqrt(hks->sfunc_sigma);

/*----------------------------------------------------------------------------*/
/* Compute correlation coefficient between H and k                            */
/*----------------------------------------------------------------------------*/
	sum_mh = 0.0;
	sum_rv = 0.0;
	sum_mhrv = 0.0;
	sum_mh2 = 0.0;
	sum_rv2 = 0.0;
	for (izr = 0; izr < hks->nzr; izr++) {			
		if (!hks->tag[izr]) {
			continue;
		}	 
		
		ir = izr % hks->nrv;
		iz = izr / hks->nrv;

		rvtmp = hks->rvlb + ir * hks->rvdt;	
		mhtmp = hks->mhlb + iz * hks->mhdt;

		sum_mh += mhtmp;
		sum_rv += rvtmp;
		sum_mhrv += rvtmp * mhtmp;
		sum_mh2 += mhtmp * mhtmp;
		sum_rv2 += rvtmp * rvtmp;			
	}
	hks->mhrv_cc = (float) ((sum_mhrv - sum_mh * sum_rv / hks->ntags) / 
			sqrt((sum_mh2 - sum_mh * sum_mh / hks->ntags) * (sum_rv2 - sum_rv * sum_rv / hks->ntags)));

	
/*----------------------------------------------------------------------------*/
/* !bootstrap for normal                                                      */
/*----------------------------------------------------------------------------*/
	if (!bootstrap) {

		if (hks->msf != NULL)  {
			delete [] hks->msf;
			hks->msf = NULL;
		}
		hks->msf = new MAX_SFUNCTIONS[MAX_MSFUNC];
	
/*----------------------------------------------------------------------------*/
/* Search for all sfunction values close to the maximum 
	Do not count points within pre-set Euclid distance                         */ 	
/*----------------------------------------------------------------------------*/
		hks->n_maxsfunc = 1;
		hks->msf[0].sfunc = hks->sfunc_max;
		hks->msf[0].kzr = kzr;
		hks->msf[0].isOK = 1;

		dist2_crt = (int)((hks->sfcrt_mh/hks->mhdt/2) * (hks->sfcrt_mh/hks->mhdt/2) + 
								(hks->sfcrt_rv/hks->rvdt/2) * (hks->sfcrt_rv/hks->rvdt/2));
		for (izr = 0; izr < hks->nzr; izr++) {
			if ((hks->sfunc_max - hks->a3ps[izr]) / hks->sfunc_max < hks->sfcrt_peak) {
				ir = izr % hks->nrv;
				iz = izr / hks->nrv;
				isOK = 1;
				for (i = 0; i < hks->n_maxsfunc; i++) {
					irv_max = hks->msf[i].kzr % hks->nrv;
					imh_max = hks->msf[i].kzr / hks->nrv;										
					dist2 = (ir - irv_max) * (ir - irv_max) + (iz - imh_max) * (iz - imh_max);
					if (i == 0) {
						if (dist2 <= dist2_crt) {
							isOK = 0;
							break;
						}
					}	
					else {
						if (dist2 <= dist2_crt) {
							if (hks->a3ps[izr] > hks->msf[i].sfunc) {
								isOK = 0;
								hks->msf[i].sfunc = hks->a3ps[izr];
								hks->msf[i].kzr = izr;
								hks->msf[i].isOK = 1;
								break;
							}
							else {
								isOK = 0;
								break;
							}	
						}
					} 	
				}
				if (isOK) {		
					hks->msf[hks->n_maxsfunc].sfunc = hks->a3ps[izr];
					hks->msf[hks->n_maxsfunc].kzr = izr;
					hks->msf[hks->n_maxsfunc].isOK = 1;
					hks->n_maxsfunc++;
				}	
				if (hks->n_maxsfunc >= MAX_MSFUNC) {
					fprintf(stderr, "\nNo. of Max_sfunc = %d > Max, increase buffer.\n\n",hks->n_maxsfunc);
					break;
				}
			}	
		}		
	
/*----------------------------------------------------------------------------*/
/* Sort the maxs of sfunction                                                 */ 	
/*----------------------------------------------------------------------------*/
		printf("No. of Max_sfunc = %d\n",hks->n_maxsfunc);
		for (i = 0; i < hks->n_maxsfunc-1; i++) {
			for (j = 0; j < hks->n_maxsfunc-1-i; j++) {
				if (hks->msf[j].sfunc < hks->msf[j+1].sfunc) {
					exsf.sfunc = hks->msf[j].sfunc;
					hks->msf[j].sfunc = hks->msf[j+1].sfunc;
					hks->msf[j+1].sfunc = exsf.sfunc;
				
					exsf.kzr = hks->msf[j].kzr;
					hks->msf[j].kzr = hks->msf[j+1].kzr;
					hks->msf[j+1].kzr = exsf.kzr;				
				}
			}
		}	

/*
		tmpsfunc = new float[hks->n_maxsfunc];
		tmpkzr = new int[hks->n_maxsfunc];
		tmpindex = new int[hks->n_maxsfunc];
		for (i = 0; i < hks->n_maxsfunc; i++) {
			tmpsfunc[i] = hks->msf[i].sfunc;
			tmpkzr[i] = hks->msf[i].kzr;
			tmpindex[i] = i;
		}	 
		quicksort(tmpsfunc, tmpindex, 0, hks->n_maxsfunc-1);
		for (i = 0; i < hks->n_maxsfunc; i++) {
			hks->msf[i].sfunc = tmpsfunc[i];
			hks->msf[i].kzr = tmpkzr[tmpindex[i]];
		}	 
		delete [] tmpsfunc;
		delete [] tmpkzr;
		delete [] tmpindex;
*/
	
/*----------------------------------------------------------------------------*/
/* Some statistics                                                            */	
/*----------------------------------------------------------------------------*/
		for (i = 0; i < hks->n_maxsfunc; i++) {
		
			hks->msf[i].krv = hks->msf[i].kzr % hks->nrv;
			hks->msf[i].kmh = hks->msf[i].kzr / hks->nrv;
		
			hks->msf[i].frv = hks->rvlb + hks->msf[i].krv * hks->rvdt;	
			hks->msf[i].fmh = hks->mhlb + hks->msf[i].kmh * hks->mhdt;
				
			tmpx1 = (hks->a3ps[hks->msf[i].kmh * hks->nrv + hks->msf[i].krv] - 
						 hks->a3ps[hks->msf[i].kmh * hks->nrv + (hks->msf[i].krv-1)]) / hks->rvdt;
			tmpx2 = (hks->a3ps[hks->msf[i].kmh * hks->nrv + (hks->msf[i].krv+1)] - 
						 hks->a3ps[hks->msf[i].kmh * hks->nrv + hks->msf[i].krv]) / hks->rvdt;
			tmpx = (tmpx1 - tmpx2) / hks->rvdt;
			hks->msf[i].sigma_rv = sqrt(2 * hks->sfunc_sigma / fabs(tmpx));

			tmpx1 = (hks->a3ps[hks->msf[i].kmh * hks->nrv + hks->msf[i].krv] - 
						 hks->a3ps[(hks->msf[i].kmh-1) * hks->nrv + hks->msf[i].krv]) / hks->mhdt;
			tmpx2 = (hks->a3ps[(hks->msf[i].kmh+1) * hks->nrv + hks->msf[i].krv] - 
						 hks->a3ps[hks->msf[i].kmh * hks->nrv + hks->msf[i].krv]) / hks->mhdt;
			tmpx = (tmpx1 - tmpx2) / hks->mhdt;
			hks->msf[i].sigma_mh = sqrt(2 * hks->sfunc_sigma / fabs(tmpx));

		}
			

/*----------------------------------------------------------------------------*/
/*	Select the maxs of sfunction 
	Rule out some Mohos have overlaps                                          */
/*----------------------------------------------------------------------------*/
		for (i = 0; i < hks->n_maxsfunc-1; i++) {
			if (hks->msf[i].isOK == 1) {
				rv_lb1 = hks->msf[i].frv - hks->msf[i].sigma_rv;
				rv_ub1 = hks->msf[i].frv + hks->msf[i].sigma_rv;
				
				mh_lb1 = hks->msf[i].fmh - hks->msf[i].sigma_mh;
				mh_ub1 = hks->msf[i].fmh + hks->msf[i].sigma_mh;
						
				for (j = i+1; j < hks->n_maxsfunc; j++) {
					if (hks->msf[j].isOK == 1) {
						rv_lb2 = hks->msf[j].frv - hks->msf[j].sigma_rv;
						rv_ub2 = hks->msf[j].frv + hks->msf[j].sigma_rv;
					
						mh_lb2 = hks->msf[j].fmh - hks->msf[j].sigma_mh;
						mh_ub2 = hks->msf[j].fmh + hks->msf[j].sigma_mh;

//						if ((mh_ub1 < mh_lb2 || mh_ub2 < mh_lb1) && (rv_ub1 < rv_lb2 || rv_ub2 < rv_lb1)) {
						if (!(mh_ub1 < mh_lb2 || mh_ub2 < mh_lb1)) {
							hks->msf[j].isOK = 0;
							//printf("Max sfunc %d neglect\n", j);
						}
					}
				}
			}
		}					
		int c_maxsfunc = 0;
		for (i = 0; i < hks->n_maxsfunc; i++) {
			if (hks->msf[i].isOK == 1) c_maxsfunc++;
		}
		printf("Total %d/%d Max sfunctions are kept: statis_option = %d\n", 
			c_maxsfunc, hks->n_maxsfunc, bootstrap);
	
	}    

/*----------------------------------------------------------------------------*/
/* for bootstrap                                                              */
/*----------------------------------------------------------------------------*/
	else {

		if (hks->msf != NULL) {
			delete [] hks->msf;
			hks->msf = NULL;
		}
		hks->msf = new MAX_SFUNCTIONS[1];

		hks->n_maxsfunc = 1;
		hks->msf[0].sfunc = hks->sfunc_max;
		hks->msf[0].kzr = kzr;
		hks->msf[0].isOK = 1;

		hks->msf[0].krv = hks->msf[0].kzr % hks->nrv;
		hks->msf[0].kmh = hks->msf[0].kzr / hks->nrv;
		
		hks->msf[0].frv = hks->rvlb + hks->msf[0].krv * hks->rvdt;	
		hks->msf[0].fmh = hks->mhlb + hks->msf[0].kmh * hks->mhdt;

	}

	return 0;
}		


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Generate H-k and S-function file name with path                            
	opt = 1: Use old name format: hk2dx7.AXX, hk2dx7.AXX.00..., sfr2dx7.AXX
	opt = 2: Use new name format: hk2d.HL.AXX.x7, hk2d.HL.AXX.x7.00..., 
				sfr2d.HL.AXX.x7                                                   */
/*----------------------------------------------------------------------------*/
int create_hkfname(char *fpth_file, char *fodir, char *hkfile, 
						const int& xc_flag, const int& opt)
{ 
	char seg1[16];
	char seg2[8];

	if (opt == 1) {
		if (segchar(hkfile, seg1, seg2, '.', 1) != 1) return -3;
		seg1[strlen(seg1)-1] = '\0';
		sprintf(fpth_file, "%s%sx%d.%s", fodir, seg1, xc_flag, seg2);	
	}
	else if (opt == 2) {
		sprintf(fpth_file, "%s%s.x%d", fodir, hkfile, xc_flag);			
	}

	return 0;
}	


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Write solution of H-k search
	Write S-function 
	hks->wsf_flag = 0: binary file
					  = 1: ascii file                                                        */
/*----------------------------------------------------------------------------*/
int wsfunc(char* fpth, HKS_INFO *hks)
{
	FILE* wfpt;
	int ir, iz;
	
	if (!hks->wsf_flag) {
	
		wfpt = fopen(fpth, "wb");
 		if (wfpt == NULL) {
   		fprintf(stderr, "\nError in wsfunc: Cannot open output file %s!!\n\n", fpth);
   		return -1;
 		}

		for (ir=0; ir < hks->nrv; ir++) {	
			for (iz=0; iz < hks->nmh; iz++) {
/*----------------------------------------------------------------------------*/
/* The sequence is to consistent old version                                  */			
/*----------------------------------------------------------------------------*/
				if (fwrite(&(hks->a3ps[iz*hks->nrv + ir]), sizeof(float), 1, wfpt) != 1) {
					fprintf(stderr, "\nError in wsfunc: Writing sfunc %s!! \n\n", fpth);
					fclose(wfpt);
					return -2;
				}				
			}		
		}

		fclose(wfpt);
	}
	
	else {

		wfpt = fopen(fpth, "w");
 		if (wfpt == NULL) {
   		fprintf(stderr, "\nError in wsfunc: Cannot open output file %s!!\n\n", fpth);
   		return -1;
 		}
	
		fprintf(wfpt, "%%K_n = %d, Moho_N = %d\n\n", hks->nrv, hks->nmh); 	
		fprintf(wfpt, "%%%5s  %6s  %9s\n", "K", "Moho", "sfunc");

		float rv, mh;
		for (rv=hks->rvlb, ir=0; ir<hks->nrv; ir++, rv+=hks->rvdt) {	
			for (mh=hks->mhlb, iz=0; iz<hks->nmh; iz++,mh+=hks->mhdt) {
/*----------------------------------------------------------------------------*/
/* The sequence is to consistent old version                                  */			
/*----------------------------------------------------------------------------*/
				fprintf(wfpt, "%6.3f  %6.2f  %11.4e\n", rv, mh, hks->a3ps[iz*hks->nrv + ir]);
			}
		
		}
	
		fclose(wfpt);		
	}
		
	return 0;
}	


/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Write solution of H-k search
	Write summary of H-k search hkr2dx7.xxx, and hkr2dx7.xxx.00, 01, ... 
	This version is for non-bootstrap                                          */
/*----------------------------------------------------------------------------*/
int wmh_rv(char* fpth, HKS_INFO *hks, RECORD *rfrec, BAZ_INFO *bazinfo)
{
	FILE* wfpt;
	FILE* wfpt2;
	char fpth_hk1[128];
	char tmpch[64], dummy[64];
	//double rho_rvmh;
	int i, irec;
	double thetr;
	float **covmat, **eigvct;
	float *eigval;
	int nrot;
	float aspect_ratio;

	extern void jacobi(float **a, int n, float d[], float **v, int *nrot);
		
/*----------------------------------------------------------------------------*/
/* Write summary of H-k search to hkr2dx7.xxx                                 */
/*----------------------------------------------------------------------------*/
	wfpt2 = fopen(fpth, "w");
	if (wfpt2 == NULL) {
 		fprintf(stderr, "\nError in wmh_rv: Cannot open output file %s!!\n\n", fpth);
 		return -1;
	}
	fprintf(wfpt2, "Selected average crust Vp = %6.4f km/s\n", hks->vp);
	fprintf(wfpt2, "Moho depth from %6.3f to %6.3f at step %6.3f with %d samples\n", 
					hks->mhlb, hks->mhub, hks->mhdt, hks->nmh);
	fprintf(wfpt2, "Vp/Vs ratio from %6.3f to %6.3f at step %6.3f with %d samples\n", 
					hks->rvlb, hks->rvub, hks->rvdt, hks->nrv);
	fprintf(wfpt2, "Weighting: t_Ps = %4.2f, t_PpPs = %4.2f, t_PpSs+PsPs = %4.2f\n", 
					hks->w_0p1s, hks->w_2p1s, hks->w_1p2s);
	fprintf(wfpt2, "Smoothing time window = %f\n", hks->twl);
	fprintf(wfpt2, "N-th root = %d\n", hks->nroot);
	fprintf(wfpt2, "Cross-correlation weighting flag = %d\n", hks->xc_flag);
	fprintf(wfpt2, "BAZ range from %5.1f to % 5.1f\n", bazinfo->lb, bazinfo->ub);				
	fprintf(wfpt2, "Search H = %3.1f K = %4.2f sfunction peaks < %4.2f\n", 
					hks->sfcrt_mh, hks->sfcrt_rv, hks->sfcrt_peak);
	if (hks->tsht_flag) fprintf(wfpt2, "Set time shift\n");
	else fprintf(wfpt2, "Unset time shift\n");
	if (hks->wsf_flag) fprintf(wfpt2, "Write sfunction in ASCII\n\n\n");
	else fprintf(wfpt2, "Write sfunction in BINARY\n\n\n");
	
/*----------------------------------------------------------------------------*/
/* Compute error ellipses for each solution                                   */
/*----------------------------------------------------------------------------*/
	covmat = matrix(1, 2, 1, 2);
	eigvct = matrix(1, 2, 1, 2);
	eigval = vector(1, 2);
	aspect_ratio = (hks->rvdt / (hks->rvub - hks->rvlb)) / (hks->mhdt / (hks->mhub - hks->mhlb));

	int ix_count = 0;
	for (i = 0; i < hks->n_maxsfunc; i++) {
		
		if (hks->msf[i].isOK == 0) {
			//printf("%d: sfunc = %e, Neglect.\n", i, hks->msf[i].sfunc);			
			continue;
		}
//		else {
//			printf("%d: sfunc = %e.\n", i, hks->msf[i].sfunc);
//		}
					
		/*rho_rvmh = 0.0;
		for (mh=hks->mhlb, iz=0; iz<hks->nmh; iz++,mh+=hks->mhdt) {
			for (rv=hks->rvlb, ir=0; ir<hks->nrv; ir++, rv+=hks->rvdt) {	
				rho_rvmh += (rv - hks->msf[i].frv)*(mh - hks->msf[i].fmh);
			}
		}	
		rho_rvmh /= (hks->nzr);
	
		thetr = atan2(2 * rho_rvmh * hks->msf[i].sigma_rv * hks->msf[i].sigma_mh, 
						hks->msf[i].sigma_mh * hks->msf[i].sigma_mh - hks->msf[i].sigma_rv * hks->msf[i].sigma_rv) / 2;
		a2 = hks->msf[i].sigma_mh * hks->msf[i].sigma_mh * hks->msf[i].sigma_rv * hks->msf[i].sigma_rv
				 * (sin(thetr) * sin(thetr) - cos(thetr) * cos(thetr)) 
				 / (hks->msf[i].sigma_mh * hks->msf[i].sigma_mh * sin(thetr) * sin(thetr) 
				 - hks->msf[i].sigma_rv * hks->msf[i].sigma_rv * cos(thetr) * cos(thetr));
		hks->msf[i].a = sqrt(fabs(a2));
		b2 = hks->msf[i].sigma_mh * hks->msf[i].sigma_mh * hks->msf[i].sigma_rv * hks->msf[i].sigma_rv
	  		 * (cos(thetr) * cos(thetr) - sin(thetr) * sin(thetr)) 
				 / (hks->msf[i].sigma_mh * hks->msf[i].sigma_mh * cos(thetr) * cos(thetr) 
				- hks->msf[i].sigma_rv * hks->msf[i].sigma_rv * sin(thetr) * sin(thetr));
		hks->msf[i].b = sqrt(fabs(b2));
		hks->msf[i].thetd = thetr * 180 / M_PI;  */
		
		hks->msf[i].sigma_mhrv = hks->mhrv_cc * hks->msf[i].sigma_rv * hks->msf[i].sigma_mh;
		covmat[1][1] = hks->msf[i].sigma_mh * hks->msf[i].sigma_mh;
		covmat[2][2] = hks->msf[i].sigma_rv * hks->msf[i].sigma_rv;
		covmat[1][2] = hks->msf[i].sigma_mhrv;
		covmat[2][1] = hks->msf[i].sigma_mhrv;
		
		jacobi(covmat, 2, eigval, eigvct, &nrot);
		
		hks->msf[i].a = sqrt(FMAX(eigval[1], eigval[2]));
		hks->msf[i].b = sqrt(FMIN(eigval[1], eigval[2]));
		thetr = atan2((1.0 / aspect_ratio) * 2.0 * hks->msf[i].sigma_mhrv,
				hks->msf[i].sigma_mh * hks->msf[i].sigma_mh - hks->msf[i].sigma_rv * hks->msf[i].sigma_rv); 
		hks->msf[i].thetd = (float)(thetr * 180 / M_PI);

/*----------------------------------------------------------------------------*/
/* Write summary of H-k search to hkr2dx7.xxx                                 */
/*----------------------------------------------------------------------------*/
		fprintf(wfpt2, "%d\n", ix_count);
		fprintf(wfpt2, "Max sfunction = %11.4e +/- %11.4e\n", 
						hks->msf[i].sfunc, hks->sfunc_sigma);
		fprintf(wfpt2, "Moho = %6.3f +/- %6.3f km, Vp/Vs = %6.3f +/- %6.3f\n", 
						hks->msf[i].fmh, hks->msf[i].sigma_mh, hks->msf[i].frv, hks->msf[i].sigma_rv);
		fprintf(wfpt2, "a = %6.3f, b = %6.3f, alpha = %6.3f\n", 
						hks->msf[i].a, hks->msf[i].b, hks->msf[i].thetd);
		fprintf(wfpt2, "\n");


/*----------------------------------------------------------------------------*/
/* Write summary of H-k search to each hkr2dx7.xxx.00, 01, ...                */
/*----------------------------------------------------------------------------*/
		sprintf(fpth_hk1, "%s.%02d", fpth, ix_count);
		 
 		wfpt = fopen(fpth_hk1, "w");
 		if (wfpt == NULL) {
   		fprintf(stderr, "\nError in wmh_rv: Cannot open output file %s!!\n\n", fpth_hk1);
   		return -1;
 		}
 	
		fprintf(wfpt, "%%Selected average crust Vp = %6.4f km/s\n", hks->vp);
		fprintf(wfpt, "%%Moho depth from %6.3f to %6.3f at step %6.3f with %d samples\n", 
						hks->mhlb, hks->mhub, hks->mhdt, hks->nmh);
		fprintf(wfpt, "%%Vp/Vs ratio from %6.3f to %6.3f at step %6.3f with %d samples\n", 
						hks->rvlb, hks->rvub, hks->rvdt, hks->nrv);
		fprintf(wfpt, "%%Weighting: t_Ps = %4.2f, t_PpPs = %4.2f, t_PpSs+PsPs = %4.2f\n", 
						hks->w_0p1s, hks->w_2p1s, hks->w_1p2s);
		fprintf(wfpt, "%%Smoothing time window = %f\n", hks->twl);
		fprintf(wfpt, "%%N-th root = %d\n", hks->nroot);
		fprintf(wfpt, "%%Cross-correlation weighting flag = %d\n", hks->xc_flag);
		fprintf(wfpt, "%%BAZ range from %5.1f to % 5.1f\n", bazinfo->lb, bazinfo->ub);				
		fprintf(wfpt, "%%Search H = %3.1f K = %4.2f sfunction peaks < %4.2f\n", 
						hks->sfcrt_mh, hks->sfcrt_rv, hks->sfcrt_peak);
		if (hks->tsht_flag) fprintf(wfpt, "%%Set time shift\n");
		else fprintf(wfpt, "%%Unset time shift\n");
		if (hks->wsf_flag) fprintf(wfpt, "%%Write sfunction in ASCII\n\n\n");
		else fprintf(wfpt, "%%Write sfunction in BINARY\n\n");
				
		fprintf(wfpt, "%%Max sfunction = %11.4e\n", hks->msf[i].sfunc);
		fprintf(wfpt, "%%sigma_sfunction = %11.4e\n", hks->sfunc_sigma);
		fprintf(wfpt, "%%Moho depth = %6.3f +/- %6.3f km\n", hks->msf[i].fmh, hks->msf[i].sigma_mh);
		fprintf(wfpt, "%%Vp/Vs ratio = %6.3f +/- %6.3f\n", hks->msf[i].frv, hks->msf[i].sigma_rv);
		fprintf(wfpt, "%%a = %6.3f, b = %6.3f, alpha = %6.3f\n\n\n", 
						hks->msf[i].a, hks->msf[i].b, hks->msf[i].thetd);
		for (irec = 0; irec < hks->nrfrec; irec++) {
			if (segchar(rfrec[irec].GetName(), dummy, tmpch, '/', 1) != 1) {
				fclose(wfpt);
				return -2;
			}	
			rfrec[irec].Ttps(hks, hks->msf[i].kzr);  //Note how the mh and rv sort into 1-D array
			rfrec[irec].LoadData();
			rfrec[irec].PeakAmp(hks);
			rfrec[irec].FreeData();
			
			fprintf(wfpt, "%s   %9.5f   %5.2f %7.4f   %5.2f %7.4f   %5.2f %7.4f\n",  
						tmpch, rfrec[irec].Ray(),
						rfrec[irec].t_0p1s, rfrec[irec].a_0p1s,
						rfrec[irec].t_2p1s, rfrec[irec].a_2p1s,
						rfrec[irec].t_1p2s, rfrec[irec].a_1p2s);
		}	

		fclose(wfpt);
		ix_count++;
	}
			
	fclose(wfpt2);

	free_matrix(covmat, 1, 2, 1, 2);
	free_matrix(eigvct, 1, 2, 1, 2);
	free_vector(eigval, 1, 2);
	covmat = NULL;
	eigvct = NULL;
	eigval = NULL;

	return 0;
	
}

/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Write solution of H-k search
	Write summary of H-k search hkr2dx7.xxx, and hkr2dx7.xxx.00, 01, ... 
	This version is for bootstrap                                              */
/*----------------------------------------------------------------------------*/
int wmh_rv2(char* fpth, HKS_INFO *hks, RECORD *rfrec, BAZ_INFO *bazinfo, 
				BOOTSTRAP &btstrp)
{
	FILE* wfpt;
	FILE* wfpt2;
	char fpth_hk1[128];
	char tmpch[64], dummy[64];
	int i, irec;
	double a2, b2, thetr;
	double a, b, thetd;
		
	wfpt2 = fopen(fpth, "w");
	if (wfpt2 == NULL) {
 		fprintf(stderr, "\nError in wmh_rv2: Cannot open output file %s!!\n\n", fpth);
 		return -1;
	}
	fprintf(wfpt2, "Neglect Vp loop\n");
	fprintf(wfpt2, "Vp/Vs ratio from %6.3f to %6.3f at step %6.3f with %d samples\n", 
					hks->rvlb, hks->rvub, hks->rvdt, hks->nrv);
	fprintf(wfpt2, "Moho depth from %6.3f to %6.3f at step %6.3f with %d samples\n", 
					hks->mhlb, hks->mhub, hks->mhdt, hks->nmh);
	fprintf(wfpt2, "BAZ range from %5.1f to % 5.1f\n", bazinfo->lb, bazinfo->ub);				
	fprintf(wfpt2, "Selected average crust Vp = %6.4f km/s\n", hks->vp);
	fprintf(wfpt2, "Weighting: t_Ps = %4.2f, t_PpPs = %4.2f, t_PpSs+PsPs = %4.2f\n", hks->w_0p1s, hks->w_2p1s, hks->w_1p2s);
	fprintf(wfpt2, "Smoothing time window = %f\n", hks->twl);
	fprintf(wfpt2, "N-th root = %d\n", hks->nroot);
	fprintf(wfpt2, "Cross-correlation weighting flag = %d\n\n\n", hks->xc_flag);


	btstrp.mh_bm = 0.0;
	btstrp.rv_bm = 0.0;
	for (i = 0; i < btstrp.n; i++) {					
		btstrp.mh_bm += btstrp.fmh[i];
		btstrp.rv_bm += btstrp.frv[i];

		fprintf(wfpt2, "%6.3f  %6.3f\n",	btstrp.fmh[i], btstrp.frv[i]);				
	}	
	btstrp.mh_bm /= btstrp.n;
	btstrp.rv_bm /= btstrp.n;
	
	btstrp.sigma_mh = 0.0;
	btstrp.sigma_rv = 0.0;
	btstrp.rho_rvmh = 0.0;
	for (i = 0; i < btstrp.n; i++) {							
		btstrp.sigma_mh += (btstrp.fmh[i] - btstrp.mh_bm) * (btstrp.fmh[i] - btstrp.mh_bm);
		btstrp.sigma_rv += (btstrp.frv[i] - btstrp.rv_bm) * (btstrp.frv[i] - btstrp.rv_bm);
		btstrp.rho_rvmh += (btstrp.frv[i] - btstrp.rv_bm) * (btstrp.fmh[i] - btstrp.mh_bm);
	}
	btstrp.sigma_mh = sqrt(btstrp.sigma_mh / (btstrp.n-1));
	btstrp.sigma_rv = sqrt(btstrp.sigma_rv / (btstrp.n-1));
	btstrp.rho_rvmh = btstrp.rho_rvmh / (btstrp.n-1);
	
			
	thetr = atan2(2 * btstrp.rho_rvmh * btstrp.sigma_rv * btstrp.sigma_mh, 
					btstrp.sigma_mh * btstrp.sigma_mh - btstrp.sigma_rv * btstrp.sigma_rv) / 2;
	a2 = btstrp.sigma_mh * btstrp.sigma_mh * btstrp.sigma_rv * btstrp.sigma_rv
			 * (sin(thetr) * sin(thetr) - cos(thetr) * cos(thetr)) 
			 / (btstrp.sigma_mh * btstrp.sigma_mh * sin(thetr) * sin(thetr) 
			 - btstrp.sigma_rv * btstrp.sigma_rv * cos(thetr) * cos(thetr));
	a = sqrt(fabs(a2));
	b2 = btstrp.sigma_mh * btstrp.sigma_mh * btstrp.sigma_rv * btstrp.sigma_rv
	  		 * (cos(thetr) * cos(thetr) - sin(thetr) * sin(thetr)) 
			 / (btstrp.sigma_mh * btstrp.sigma_mh * cos(thetr) * cos(thetr) 
			- btstrp.sigma_rv * btstrp.sigma_rv * sin(thetr) * sin(thetr));
	b = sqrt(fabs(b2));
	thetd = thetr * 180 / M_PI;

	fprintf(wfpt2, "\nMoho = %6.3f +/- %6.3f km, Vp/Vs = %6.3f +/- %6.3f\n", 
				btstrp.mh_bm, btstrp.sigma_mh, btstrp.rv_bm, btstrp.sigma_rv);
	fprintf(wfpt2, "a = %6.3f, b = %6.3f, alpha = %6.3f\n", a, b, thetd);
	fprintf(wfpt2, "\n");
	fclose(wfpt2);


	sprintf(fpth_hk1, "%s.00", fpth);
	wfpt = fopen(fpth_hk1, "w");
	if (wfpt == NULL) {
 		fprintf(stderr, "\nError wmh_rv2: Cannot open output file %s!!\n\n", fpth_hk1);
  		return -1;
	}
 	
	fprintf(wfpt, "%%Neglect Vp loop\n");
	fprintf(wfpt, "%%Vp/Vs ratio from %6.3f to %6.3f at step %6.3f with %d samples\n", 
					hks->rvlb, hks->rvub, hks->rvdt, hks->nrv);
	fprintf(wfpt, "%%Moho depth from %6.3f to %6.3f at step %6.3f with %d samples\n", 
					hks->mhlb, hks->mhub, hks->mhdt, hks->nmh);
	fprintf(wfpt, "%%BAZ range from %5.1f to % 5.1f\n", bazinfo->lb, bazinfo->ub);				
	fprintf(wfpt, "%%Selected average crust Vp = %6.4f km/s\n", hks->vp);
	fprintf(wfpt, "%%Weighting: t_Ps = %4.2f, t_PpPs = %4.2f, t_PpSs+PsPs = %4.2f\n", hks->w_0p1s, hks->w_2p1s, hks->w_1p2s);
	fprintf(wfpt, "%%Smoothing time window = %f\n", hks->twl);
	fprintf(wfpt, "%%N-th root = %d\n", hks->nroot);
	fprintf(wfpt, "%%Cross-correlation weighting flag = %d\n", hks->xc_flag);
	fprintf(wfpt, "%%Max sfunction = %11.4e\n", hks->msf[0].sfunc);
	//fprintf(wfpt, "%%sigma_sfunction = %11.4e\n", hks->sfunc_sigma);
	fprintf(wfpt, "%%Vp/Vs ratio = %6.3f +/- %6.3f\n", btstrp.rv_bm, btstrp.sigma_rv);
	fprintf(wfpt, "%%Moho depth = %6.3f +/- %6.3f km\n", btstrp.mh_bm, btstrp.sigma_mh);
	fprintf(wfpt, "%%a = %6.3f, b = %6.3f, alpha = %6.3f\n", a, b, thetd);
	fprintf(wfpt, "\n");
	for (irec = 0; irec < hks->nrfrec; irec++) {
		if (segchar(rfrec[irec].GetName(), dummy, tmpch, '/', 1) != 1) {
			fclose(wfpt);
			return -2;
		}	
		rfrec[irec].Ttps_rvmh(hks->vp, btstrp.mh_bm, btstrp.rv_bm);  //Note how the mh and rv sort into 1-D array
		rfrec[irec].PeakAmp(hks);

		fprintf(wfpt, "%s   %9.5f   %5.2f %7.4f   %5.2f %7.4f   %5.2f %7.4f\n",  
					tmpch, rfrec[irec].Ray(),
					rfrec[irec].t_0p1s, rfrec[irec].a_0p1s,
					rfrec[irec].t_2p1s, rfrec[irec].a_2p1s,
					rfrec[irec].t_1p2s, rfrec[irec].a_1p2s);
	}	

	fclose(wfpt);
	return 0;
	
}

