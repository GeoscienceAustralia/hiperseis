#include "recordrf.h"
#include "nrutil.h"

#ifndef	_ZCC_L_
#define	_ZCC_L_		5.0     	/* For Moho depth in km	*/
#endif

#ifndef	_RCC_L_
#define	_RCC_L_		0.2     	/* For Vp/Vs	*/
#endif

#ifndef	_CD_FACTOR_
#define	_CD_FACTOR_	 4	
#endif

static double _acc_(int n, float* data1, float* data2);


/******************************************************************************/
/* wgt_xc() combination of wgtrv_xc and wgtmh_xc                              */
/******************************************************************************/
int wgt_xc(HKS_INFO *hks)
{
	int iz, ir, izr;
	int it, itk;
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
//	double cc, c11, c12, c13, c22, c23, c33, rc12, rc13, wgt;
	double cc, c11, c12, c13, c22, c33, rc12, rc13, wgt;
	
	imh_max0 = NULL;
	irv_max0 = NULL;
	xc_rv = NULL;
	xc_mh = NULL;
	ar_0p1ss = NULL;
	ar_0p1s = NULL;
	ar_2p1s = NULL;
	ar_1p2s = NULL;
		
/*----------------------------------------------------------------------------*/
/* Part I is for xcorr as a function of rv: wgtrv_xc()                        */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* nmh2 is half range in points for search and summation for sfunctions
	along 'H' of a given 'K'                                                   */
/*----------------------------------------------------------------------------*/
	nmh2 = (int)rint(_ZCC_L_ * 0.5 / hks->mhdt);
	if (nmh2 < 1 ) {
		fprintf(stderr, "Error in wgt_xc - Part I: nmh2 = %3d, %8.5f\n", nmh2, hks->mhdt);
		return -1;
	}	
   
/*----------------------------------------------------------------------------*/
/* Assign memory for the locations of sfunction peaks for each 'K'            */
/*----------------------------------------------------------------------------*/
	imh_max0 = new int[hks->nrv]; 
	if (imh_max0 == NULL) {
		fprintf(stderr, "wgt_xc - Part I: Error in new int imh_max0!\n\n");
		return -1;
	}	
	 
/*----------------------------------------------------------------------------*/
/* For each vp/vs ratio, find the sfunction peak in the 0p1s mode
	Store their location in imh_max0                                           */
/*----------------------------------------------------------------------------*/
	for (ir = 0; ir < hks->nrv; ir++) {
		mh_max = -MAXFLOAT;
		for (iz = 0; iz < hks->nmh; iz++) {
			izr = iz * hks->nrv + ir;
			if (mh_max < hks->a_0p1s[izr]) {
				mh_max = hks->a_0p1s[izr];
				imh_max0[ir] = iz;
			}
		}
	}

/*----------------------------------------------------------------------------*/
/* Assign memory for cross-correlation of all 'K's (xc_rv),
	sfunction amplitude of all 'H' (ar_0p1ss, ar_0p1s, ar_2p1s, ar_1p2s)
	for each 'K'                                                               */
/*----------------------------------------------------------------------------*/
	if (xc_rv != NULL) { 
		delete [] xc_rv;
		xc_rv = NULL;
	}	
	xc_rv = new float[hks->nrv];
	if (xc_rv == NULL) {
		fprintf(stderr, "wgt_xc - Part I: Error in new float xc_rv!\n\n");
		return -1;
	}	
		
	if (ar_0p1ss != NULL) {
		delete [] ar_0p1ss;
		ar_0p1ss = NULL;
	}	
	ar_0p1ss = new float[hks->nmh];
	if (ar_0p1ss == NULL) {
		fprintf(stderr, "wgt_xc - Part I: Error in new float ar_0p1ss!\n\n");
		return -1;
	}
		
	if (ar_0p1s != NULL) {
		delete [] ar_0p1s;
		ar_0p1s = NULL;
	}	
	ar_0p1s = new float[hks->nmh];
	if (ar_0p1s == NULL) {
		fprintf(stderr, "wgt_xc - Part I: Error in new float ar_0p1s!\n\n");
		return -1;
	}
		
	if (ar_2p1s != NULL) {
		delete [] ar_2p1s;
		ar_2p1s = NULL;
	}	
	ar_2p1s = new float[hks->nmh];
	if (ar_2p1s == NULL) {
		fprintf(stderr, "wgt_xc - Part I: Error in new float ar_2p1s!\n\n");
		return -1;
	}
		
	if (ar_1p2s != NULL) {
		delete [] ar_1p2s;
		ar_1p2s = NULL;
	}	
	ar_1p2s = new float[hks->nmh];
	if (ar_1p2s == NULL) {
		fprintf(stderr, "wgt_xc - Part I: Error in new float ar_1p2s!\n\n");
		return -1;
	}	
	
/*----------------------------------------------------------------------------*/
/* Start the loop to compute cross correlation for each 'K'                   */
/*----------------------------------------------------------------------------*/
	for (ir = 0; ir < hks->nrv; ir++) {

/*----------------------------------------------------------------------------*/
/* For a given 'K', copy sfunctions of all 'H' to ar_0p1s, ar_2p1s, 
	and ar_1p2s. 
	Note in each loop of 'K', they take different sfunction values             */
/*----------------------------------------------------------------------------*/
		for (iz = 0; iz < hks->nmh; iz++) {
			izr = iz * hks->nrv + ir;
   	   ar_0p1s[iz] = hks->a_0p1s[izr];
      	ar_2p1s[iz] = hks->a_2p1s[izr];
			ar_1p2s[iz] = hks->a_1p2s[izr];
			ar_0p1ss[iz] = 0.0;
		}
		
/*----------------------------------------------------------------------------*/
/* For each 'K', the location of sfunction peak in the 0p1s mode is stored 
	in imh_max0                                           
	For a given 'K' loop, do a 5-points smoothing for imh_max0
	The average is stored in imh_max, which is different in each 'K' loop.     */
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
/* For each 'K', generate ar_0p1ss, which has non-zero values centered 
	at imh_max with the range of +/-nmh2. 
	The value of each ar_0p1ss is summation over _CD_FACTOR_ points.
	For example:
	imh_max = 9
	nmh2 = 3
	_CD_FACTOR = 4
	lfmh = imh_max - nmh2 = 9 - 3 = 6
	rhmh = imh_max + nmh2 = 9 + 3 = 12
	for (i = lfmh; i <= rhmh; i++) 
		it = (i - imh_max) * _CD_FACTOR_ + imh_max
		i = 6, it = (6 - 9) * 4 + 9  = -3
		i = 7, it = (7 - 9) * 4 + 9  =  1
		i = 8, it = (8 - 9) * 4 + 9  =  5
		i = 9, it = (9 - 9) * 4 + 9  =  9
		i =10, it = (10 - 9) * 4 + 9 = 13
		i =11, it = (11 - 9) * 4 + 9 = 17
		i =12, it = (12 - 9) * 4 + 9 = 21
	The it(i) = it(i-1) + _CD_FACTOR, define the starting point of summation
	ar_0p1ss[i = 6] = 0, because it < 0
	ar_0p1ss[i = 7] = sum(ar_0p1s[1] ... ar_0p1s[4])
	ar_0p1ss[i = 8] = sum(ar_0p1s[5] ... ar_0p1s[8])
	ar_0p1ss[i = 9] = sum(ar_0p1s[9] ... ar_0p1s[12])
	ar_0p1ss[i =10] = sum(ar_0p1s[13] ... ar_0p1s[16])
	ar_0p1ss[i =11] = sum(ar_0p1s[17] ... ar_0p1s[20])
	ar_0p1ss[i =12] = sum(ar_0p1s[21] ... ar_0p1s[24])
	Such that, the non-zero values of ar_0p1ss occupy 
	[imh_max-nmh2 : imh_max+nmh2], but each value is the sum of _CD_FACTOR_ 
	points all over the entire sfunction along 'H'. See ar_0p1ss[imh_max = 9]
	to understand the layout of ar_0p1ss
	Normally wgt = 1                                                           */
/*----------------------------------------------------------------------------*/
		lfmh = IMAX(imh_max - nmh2, 0);
		rhmh = IMIN(imh_max + nmh2, hks->nmh-1);
		nmh = rhmh - lfmh + 1;
		wgt = (double) nmh / (double)(nmh2 * 2 + 1);
		
		for (i = lfmh; i <= rhmh; i++) {
			it = (i - imh_max) * _CD_FACTOR_ + imh_max;
			if (it >= 0 && it < hks->nmh) {
				for (k = 0; k < _CD_FACTOR_; k++) {
					itk = it + k;
					if (itk >= hks->nmh) break;	
					ar_0p1ss[i] += ar_0p1s[itk];
				}
			//a_0p1ss[ii] /= (float)k; 
			}
		//printf("ii = %3d it = %3d %10.5f\n", ii, it, a_0p1ss[ii]);
		}

/*----------------------------------------------------------------------------*/
/* For each 'K', calculate all cross correlations 
	Note the non-value range of ar_0p1ss   
	Obtain a cross correlation coefficient xc_rv[ir] for each 'K'            */
/*----------------------------------------------------------------------------*/
		c11 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_0p1ss[lfmh]);
		c22 = _acc_(nmh, &ar_2p1s[lfmh],  &ar_2p1s[lfmh]);
		c33 = _acc_(nmh, &ar_1p2s[lfmh],  &ar_1p2s[lfmh]);
		c12 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_2p1s[lfmh]);
		c13 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_1p2s[lfmh]);
//		c23 = _acc_(nmh, &ar_2p1s[lfmh],  &ar_1p2s[lfmh]);	

		rc12 = c12/sqrt(c11 * c22); 
		rc13 = c13/sqrt(c11 * c33);
		switch(hks->xc_flag) {
			case 7:
//				xc_rv[ir] = (c12/sqrt(c11 * c22) + c13/sqrt(c11 * c33) + c23/sqrt(c22 * c33))/3. * wgt; 
		    	xc_rv[ir] = (rc12 + rc13) * 0.5 * wgt;
		    	break;
			case 8:
	   	 	xc_rv[ir] = rc12 * wgt;
	    		break;
			case 9:
	    		xc_rv[ir] = rc13 * wgt;
	    		break;
		}
			
	} /* ratio loop */


/*----------------------------------------------------------------------------*/
/* Part II is for xcorr as a function of mh: wgtmh_xc()                       */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* nrv2 is half range in points for search and summation for sfunctions
	along 'K' of a given 'H'                                                   */
/*----------------------------------------------------------------------------*/
	nrv2 = (int)rint(_RCC_L_ * 0.5 / hks->rvdt);
	if (nrv2 < 1 ) {
		fprintf(stderr, "Error in wgt_xc - Part II: nrv2 = %3d, %8.5f\n", nrv2, hks->rvdt);
		return -1;
	}	
   
/*----------------------------------------------------------------------------*/
/* Assign memory for the locations of sfunction peaks for each 'H'            */
/*----------------------------------------------------------------------------*/
	irv_max0 = new int[hks->nmh];  
	if (irv_max0 == NULL) {
		fprintf(stderr, "wgt_xc - Part II: Error in new int irv_max0!\n\n");
		return -1;
	}	

/*----------------------------------------------------------------------------*/
/* For each 'H', find the sfunction peak in the 0p1s mode
	Store their location in irv_max0                                           */
/*----------------------------------------------------------------------------*/
	for (iz = 0; iz < hks->nmh; iz++) {
		rv_max = -MAXFLOAT;
		for (ir = 0; ir < hks->nrv; ir++) {
			izr = iz * hks->nrv + ir;
			if (rv_max < hks->a_0p1s[izr]) {
				rv_max = hks->a_0p1s[izr];
				irv_max0[iz] = ir;
			}
		}
	}
  
/*----------------------------------------------------------------------------*/
/* Assign memory for cross-correlation of all 'H's (xc_mh),
	sfunction amplitude of all 'K's (ar_0p1ss, ar_0p1s, ar_2p1s, ar_1p2s)
	for each 'H'                                                               */
/*----------------------------------------------------------------------------*/  
	if (xc_mh != NULL) { 
		delete [] xc_mh;
		xc_mh = NULL;
	}	
	xc_mh = new float[hks->nmh];
	if (xc_mh == NULL) {
		fprintf(stderr, "wgt_xc - Part II: Error in new float xc_mh!\n\n");
		return -1;
	}	
		
	if (ar_0p1ss != NULL) {
		delete [] ar_0p1ss;
		ar_0p1ss = NULL;
	}	
	ar_0p1ss = new float[hks->nrv];
	if (ar_0p1ss == NULL) {
		fprintf(stderr, "wgt_xc - Part II: Error in new float ar_0p1ss!\n\n");
		return -1;
	}
		
	if (ar_0p1s != NULL) {
		delete [] ar_0p1s;
		ar_0p1s = NULL;
	}	
	ar_0p1s = new float[hks->nrv];
	if (ar_0p1s == NULL) {
		fprintf(stderr, "wgt_xc - Part II: Error in new float ar_0p1s!\n\n");
		return -1;
	}
		
	if (ar_2p1s != NULL) {
		delete [] ar_2p1s;
		ar_2p1s = NULL;
	}	
	ar_2p1s = new float[hks->nrv];
	if (ar_2p1s == NULL) {
		fprintf(stderr, "wgt_xc - Part II: Error in new float ar_2p1s!\n\n");
		return -1;
	}
		
	if (ar_1p2s != NULL) {
		delete [] ar_1p2s;
		ar_1p2s = NULL;
	}	
	ar_1p2s = new float[hks->nrv];
	if (ar_1p2s == NULL) {
		fprintf(stderr, "wgt_xc - Part II: Error in new float ar_1p2s!\n\n");
		return -1;
	}	
		
/*----------------------------------------------------------------------------*/
/* Start the loop to compute cross correlation for each 'H'                   */
/*----------------------------------------------------------------------------*/	
	for (iz = 0; iz < hks->nmh; iz++) {

/*----------------------------------------------------------------------------*/
/* For a given 'H', copy sfunctions of all 'K' to ar_0p1s, ar_2p1s, 
	and ar_1p2s. 
	Note in each loop of 'H', they take different sfunction values             */
/*----------------------------------------------------------------------------*/
		for (ir = 0; ir < hks->nrv; ir++) {
			izr = iz * hks->nrv + ir;
   	   ar_0p1s[ir] = hks->a_0p1s[izr];
      	ar_2p1s[ir] = hks->a_2p1s[izr];
			ar_1p2s[ir] = hks->a_1p2s[izr];
			ar_0p1ss[ir] = 0.0;
		}
		
/*----------------------------------------------------------------------------*/
/* For each 'H', the location of sfunction peak in the 0p1s mode is stored 
	in irv_max0                                           
	For a given 'H' loop, do a 5-points smoothing for irv_max0
	The average is stored in irv_max, which is different in each 'H' loop.     */
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
/* For each 'H', generate ar_0p1ss, which has non-zero values centered 
	at irv_max with the range of +/-nrv2. 
	The value of each ar_0p1ss is summation over _CD_FACTOR_ points.
	For example:
	irv_max = 9
	nrv2 = 3
	_CD_FACTOR = 4
	bmrv = irv_max - nrv2 = 9 - 3 = 6
	tprv = irv_max + nrv2 = 9 + 3 = 12
	for (i = bmrv; i <= tprv; i++) 
		it = (i - irv_max) * _CD_FACTOR_ + irv_max
		i = 6, it = (6 - 9) * 4 + 9  = -3
		i = 7, it = (7 - 9) * 4 + 9  =  1
		i = 8, it = (8 - 9) * 4 + 9  =  5
		i = 9, it = (9 - 9) * 4 + 9  =  9
		i =10, it = (10 - 9) * 4 + 9 = 13
		i =11, it = (11 - 9) * 4 + 9 = 17
		i =12, it = (12 - 9) * 4 + 9 = 21
	The it(i) = it(i-1) + _CD_FACTOR, define the starting point of summation
	ar_0p1ss[i = 6] = 0, because it < 0
	ar_0p1ss[i = 7] = sum(ar_0p1s[1] ... ar_0p1s[4])
	ar_0p1ss[i = 8] = sum(ar_0p1s[5] ... ar_0p1s[8])
	ar_0p1ss[i = 9] = sum(ar_0p1s[9] ... ar_0p1s[12])
	ar_0p1ss[i =10] = sum(ar_0p1s[13] ... ar_0p1s[16])
	ar_0p1ss[i =11] = sum(ar_0p1s[17] ... ar_0p1s[20])
	ar_0p1ss[i =12] = sum(ar_0p1s[21] ... ar_0p1s[24])
	Such that, the non-zero values of ar_0p1ss occupy 
	[irv_max-nrv2 : irv_max+nrv2], but each value is the sum of _CD_FACTOR_ 
	points all over the entire sfunction along 'K'. See ar_0p1ss[irv_max = 9]
	to understand the layout of ar_0p1ss
	Normally wgt = 1                                                           */
/*----------------------------------------------------------------------------*/
		bmrv = IMAX(irv_max - nrv2, 0);
		tprv = IMIN(irv_max + nrv2, hks->nrv-1);
		nrv = tprv - bmrv + 1;
		wgt = (double) nrv / (double)(nrv2 * 2 + 1);
		
		for (i = bmrv; i <= tprv; i++) {
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

/*----------------------------------------------------------------------------*/
/* For each 'H', calculate all cross correlations 
	Note the non-value range of ar_0p1ss   
	Obtain a cross correlation coefficient xc_mh[iz] for each 'H'            */
/*----------------------------------------------------------------------------*/
		c11 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_0p1ss[bmrv]);
		c22 = _acc_(nrv, &ar_2p1s[bmrv],  &ar_2p1s[bmrv]);
		c33 = _acc_(nrv, &ar_1p2s[bmrv],  &ar_1p2s[bmrv]);
		c12 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_2p1s[bmrv]);
		c13 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_1p2s[bmrv]);
//		c23 = _acc_(nrv, &ar_2p1s[bmrv],  &ar_1p2s[bmrv]);	

		rc12 = c12/sqrt(c11 * c22); 
		rc13 = c13/sqrt(c11 * c33);
		switch(hks->xc_flag) {
			case 7:
//				xc_mh[iz] = (c12/sqrt(c11 * c22) + c13/sqrt(c11 * c33) + c23/sqrt(c22 * c33))/3. * wgt; 
		    	xc_mh[iz] = (rc12 + rc13) * 0.5 * wgt;
		    	break;
			case 8:
	   	 	xc_mh[iz] = rc12 * wgt;
	    		break;
			case 9:
				xc_mh[iz] = rc13 * wgt;
				break;
		}
		
	} /* moho loop */
		
		
/*----------------------------------------------------------------------------*/
/* After obtain a cross correlation coefficient xc_rv[ir] for each 'K', and
	xc_mh[iz] for each 'H', conduct a new 'K' loop to smooth 
	xc_rv over 5 xc_rv along 'K', and a new 'H' loop to smooth 
	xc_mh over 5 xc_mh along 'H'.
	And then combine xc_rv and xc_mh                                           */
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
/*----------------------------------------------------------------------------*/
/* For a given 'K', multiple xc_rv[ir] on sfunction of the three phases 
	along 'H' for the given 'K'.
	For different 'K', xc_rv[ir] is different                                */
/*----------------------------------------------------------------------------*/
		for (iz = 0; iz < hks->nmh; iz++) {
			izr = iz * hks->nrv + ir;
			hks->a_0p1s[izr] *= cc;
			hks->a_2p1s[izr] *= cc;
			hks->a_1p2s[izr] *= cc;
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
/*----------------------------------------------------------------------------*/
/* For a given 'H', multiple xc_mh[iz] on sfunction of the three phases 
	along 'K' for the given 'H'.
	For different 'H', xc_mh[iz] is different                                */
/*----------------------------------------------------------------------------*/
		for (ir = 0; ir < hks->nrv; ir++) {
			izr = iz * hks->nrv + ir;
			hks->a_0p1s[izr] *= cc;
			hks->a_2p1s[izr] *= cc;
			hks->a_1p2s[izr] *= cc;
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
/* wgtrv_xc() xcorr is a function of rv                                       */
/******************************************************************************/
int wgtrv_xc(HKS_INFO *hks)
{
	int iz, ir, izr;
	int it, itk;
	int i, k; 
	int *imh_max0, imh_max;
	int lfmh, rhmh;
	int nmh2, nmh;
	float	mh_max;
//	float rv;
	float *xc_rv;
	float	*ar_0p1ss, *ar_0p1s, *ar_2p1s, *ar_1p2s;
//	double cc, c11, c12, c13, c22, c23, c33, rc12, rc13, wgt;
	double cc, c11, c12, c13, c22, c33, rc12, rc13, wgt;
	char fa0p1ss[16];
	FILE* wfpt;
	
	extern int wsfunc_3phs(char* fpth, HKS_INFO *hks, int opt);
	
	imh_max0 = NULL;
	xc_rv = NULL;
	ar_0p1ss = NULL;
	ar_0p1s = NULL;
	ar_2p1s = NULL;
	ar_1p2s = NULL;
		
/*----------------------------------------------------------------------------*/
/* nmh2 is half range in points for search and summation for sfunctions
	along 'H' of a given 'K'                                                   */
/*----------------------------------------------------------------------------*/
	nmh2 = (int)rint(_ZCC_L_ * 0.5 / hks->mhdt);
	if (nmh2 < 1 ) {
		fprintf(stderr, "Error in wgtrv_xc: nmh2 = %3d %8.5f\n", nmh2, hks->mhdt);
		return -1;
	}	
   
/*----------------------------------------------------------------------------*/
/* Assign memory for the locations of sfunction peaks for each 'K'            */
/*----------------------------------------------------------------------------*/
	imh_max0 = new int[hks->nrv];  
	if (imh_max0 == NULL) {
		fprintf(stderr, "wgtrv_xc: Error in new int imh_max0!\n\n");
		return -1;
	}	

/*----------------------------------------------------------------------------*/
/* For each vp/vs ratio, find the sfunction peak in the 0p1s mode
	Store their location in imh_max0                                           */
/*----------------------------------------------------------------------------*/
	for (ir = 0; ir < hks->nrv; ir++) {
		mh_max = -MAXFLOAT;
		for (iz = 0; iz < hks->nmh; iz++) {
			izr = iz * hks->nrv + ir;
			if (mh_max < hks->a_0p1s[izr]) {
				mh_max = hks->a_0p1s[izr];
				imh_max0[ir] = iz;
			}
		}
	}
  
/*----------------------------------------------------------------------------*/
/* Assign memory for cross-correlation of all 'K's (xc_rv),
	sfunction amplitude of all 'H' (ar_0p1ss, ar_0p1s, ar_2p1s, ar_1p2s)
	for each 'K'                                                               */
/*----------------------------------------------------------------------------*/
	if (xc_rv != NULL) { 
		delete [] xc_rv;
		xc_rv = NULL;
	}	
	xc_rv = new float[hks->nrv];
	if (xc_rv == NULL) {
		fprintf(stderr, "wgtrv_xc: Error in new float xc_rv!\n\n");
		return -1;
	}	
	
	ar_0p1ss = new float[hks->nmh];
	if (ar_0p1ss == NULL) {
		fprintf(stderr, "wgtrv_xc: Error in new float ar_0p1ss!\n\n");
		return -1;
	}	
	ar_0p1s = new float[hks->nmh];
	if (ar_0p1s == NULL) {
		fprintf(stderr, "wgtrv_xc: Error in new float ar_0p1s!\n\n");
		return -1;
	}	
	ar_2p1s = new float[hks->nmh];
	if (ar_2p1s == NULL) {
		fprintf(stderr, "wgtrv_xc: Error in new float ar_2p1s!\n\n");
		return -1;
	}	
	ar_1p2s = new float[hks->nmh];
	if (ar_1p2s == NULL) {
		fprintf(stderr, "wgtrv_xc: Error in new float ar_1p2s!\n\n");
		return -1;
	}	
	
/*----------------------------------------------------------------------------*/
/* If hks->verbose == 1, open a file to write ar_0p1ss for illustration.
	First test if file = 'wgtrv_a0p1ss' exists. If yes, clean the content      */
/*----------------------------------------------------------------------------*/
	if (hks->verbose) {
		strcpy(fa0p1ss, "wgtrv_a0p1ss");
		wfpt = fopen(fa0p1ss, "wb");
 		if (wfpt == NULL) {
   		fprintf(stderr, "\nError in wgtrv_xc: Cannot open output file 'wgtrv_a0p1ss!!\n\n");
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
			delete [] imh_max0;  
			imh_max0 = NULL;  		
   		return -2;
 		}
 		fclose(wfpt);
	}

/*----------------------------------------------------------------------------*/
/* Start the loop to compute cross correlation for each 'K'                   */
/*----------------------------------------------------------------------------*/
//	rv = hks->rvlb; 
//	hks->xc_max = 0.0;
	for (ir = 0; ir < hks->nrv; ir++) {

/*----------------------------------------------------------------------------*/
/* For a given 'K', copy sfunctions of all 'H' to ar_0p1s, ar_2p1s, 
	and ar_1p2s. 
	Note in each loop of 'K', they take different sfunction values             */
/*----------------------------------------------------------------------------*/
		for (iz = 0; iz < hks->nmh; iz++) {
			izr = iz * hks->nrv + ir;
   	   ar_0p1s[iz] = hks->a_0p1s[izr];
      	ar_2p1s[iz] = hks->a_2p1s[izr];
			ar_1p2s[iz] = hks->a_1p2s[izr];
			ar_0p1ss[iz] = 0.0;
		}
		
/*----------------------------------------------------------------------------*/
/* For each 'K', the location of sfunction peak in the 0p1s mode is stored 
	in imh_max0                                           
	For a given 'K' loop, do a 5-points smoothing for imh_max0
	The average is stored in imh_max, which is different in each 'K' loop.     */
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
/* For each 'K', generate ar_0p1ss, which has non-zero values centered 
	at imh_max with the range of +/-nmh2. 
	The value of each ar_0p1ss is summation over _CD_FACTOR_ points.
	For example:
	imh_max = 9
	nmh2 = 3
	_CD_FACTOR = 4
	lfmh = imh_max - nmh2 = 9 - 3 = 6
	rhmh = imh_max + nmh2 = 9 + 3 = 12
	for (i = lfmh; i <= rhmh; i++) 
		it = (i - imh_max) * _CD_FACTOR_ + imh_max
		i = 6, it = (6 - 9) * 4 + 9  = -3
		i = 7, it = (7 - 9) * 4 + 9  =  1
		i = 8, it = (8 - 9) * 4 + 9  =  5
		i = 9, it = (9 - 9) * 4 + 9  =  9
		i =10, it = (10 - 9) * 4 + 9 = 13
		i =11, it = (11 - 9) * 4 + 9 = 17
		i =12, it = (12 - 9) * 4 + 9 = 21
	The it(i) = it(i-1) + _CD_FACTOR, define the starting point of summation
	ar_0p1ss[i = 6] = 0, because it < 0
	ar_0p1ss[i = 7] = sum(ar_0p1s[1] ... ar_0p1s[4])
	ar_0p1ss[i = 8] = sum(ar_0p1s[5] ... ar_0p1s[8])
	ar_0p1ss[i = 9] = sum(ar_0p1s[9] ... ar_0p1s[12])
	ar_0p1ss[i =10] = sum(ar_0p1s[13] ... ar_0p1s[16])
	ar_0p1ss[i =11] = sum(ar_0p1s[17] ... ar_0p1s[20])
	ar_0p1ss[i =12] = sum(ar_0p1s[21] ... ar_0p1s[24])
	Such that, the non-zero values of ar_0p1ss occupy 
	[imh_max-nmh2 : imh_max+nmh2], but each value is the sum of _CD_FACTOR_ 
	points all over the entire sfunction along 'H'. See ar_0p1ss[imh_max = 9]
	to understand the layout of ar_0p1ss
	Normally wgt = 1                                                           */
/*----------------------------------------------------------------------------*/
		lfmh = IMAX(imh_max - nmh2, 0);
		rhmh = IMIN(imh_max + nmh2, hks->nmh-1);
		nmh = rhmh - lfmh + 1;
		wgt = (double) nmh / (double)(nmh2 * 2 + 1);
		
		for (i = lfmh; i <= rhmh; i++) {
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

/*----------------------------------------------------------------------------*/
/* If hks->verbose == 1, open a file to write ar_0p1ss for illustration.
	First test if file = 'wgtrv_a0p1ss' exists. If yes, clean the content      */
/*----------------------------------------------------------------------------*/
		if (hks->verbose) {
			wfpt = fopen(fa0p1ss, "ab");
 			if (wfpt == NULL) {
   			fprintf(stderr, "\nError in wgtrv_xc: Cannot append output file 'wgtrv_a0p1ss!!\n\n");
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
				delete [] imh_max0;  
				imh_max0 = NULL;  		
   			return -2;
	 		}
			for (iz = 0; iz < hks->nmh; iz++) {
				if (fwrite(&(ar_0p1ss[iz]), sizeof(float), 1, wfpt) != 1) {
					fprintf(stderr, "\nError in wgtrv_xc: Writing 'wgtrv_a0p1ss!! \n\n");
					fclose(wfpt);
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
					delete [] imh_max0;  
					imh_max0 = NULL;  							
					return -2;
				}				
			}
 			fclose(wfpt);
		}
		
/*----------------------------------------------------------------------------*/
/* For each 'K', calculate all cross correlations 
	Note the non-value range of ar_0p1ss   
	Obtain a cross correlation coefficient xc_rv[ir] for each 'K'            */
/*----------------------------------------------------------------------------*/
		c11 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_0p1ss[lfmh]);
		c22 = _acc_(nmh, &ar_2p1s[lfmh],  &ar_2p1s[lfmh]);
		c33 = _acc_(nmh, &ar_1p2s[lfmh],  &ar_1p2s[lfmh]);
		c12 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_2p1s[lfmh]);
		c13 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_1p2s[lfmh]);
//		c23 = _acc_(nmh, &ar_2p1s[lfmh],  &ar_1p2s[lfmh]);	

		rc12 = c12/sqrt(c11 * c22); 
		rc13 = c13/sqrt(c11 * c33);
		switch(hks->xc_flag) {
			case 1:
//				xc_rv[ir] = (c12/sqrt(c11 * c22) + c13/sqrt(c11 * c33) + c23/sqrt(c22 * c33))/3. * wgt; 
				xc_rv[ir] = (rc12 + rc13) * 0.5 * wgt;
				break;
			case 2:
				xc_rv[ir] = rc12 * wgt;
				break;
			case 3:
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
/* After obtain a cross correlation coefficient xc_rv[ir] for each 'K',
	conduct a new 'K' loop to smooth xc_rv over 5 xc_rv along 'K'          */
/*----------------------------------------------------------------------------*/
	for (ir = 0; ir < hks->nrv; ir++) {
		if (ir == 0 || ir == hks->nrv - 1) {
			cc = xc_rv[ir] + 1.0;
		}
		else if (ir == 1 || ir == hks->nrv - 2) {
			cc = 0.50 * xc_rv[ir] + 0.25 * (xc_rv[ir-1] + xc_rv[ir+1]) + 1.0;
		}
		else{
      	cc = 0.40 * xc_rv[ir] + 0.20 * (xc_rv[ir-1] + xc_rv[ir+1]) + 
				  0.10 * (xc_rv[ir-2] + xc_rv[ir+2]) + 1.0;
		}
/*----------------------------------------------------------------------------*/
/* For a given 'K', multiple xc_rv[ir] on sfunction of the three phases 
	along 'H' for the given 'K'.
	For different 'K', xc_rv[ir] is different                                */
/*----------------------------------------------------------------------------*/
		for (iz = 0; iz < hks->nmh; iz++) {
			izr = iz * hks->nrv + ir;
			hks->a_0p1s[izr] *= cc;
			hks->a_2p1s[izr] *= cc;
			hks->a_1p2s[izr] *= cc;
		}
	}
	
	if (hks->verbose) {
		strcpy(fa0p1ss, "wgtrv_a0p1s");
		if (wsfunc_3phs(fa0p1ss, hks, 1)) {
			fprintf(stderr, "\nError in wgtrc_xc: writing a0p1s!!\n\n");
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
			delete [] imh_max0;  
			imh_max0 = NULL;
 			return -2;
		}
		
		strcpy(fa0p1ss, "wgtrv_a2p1s");
		if (wsfunc_3phs(fa0p1ss, hks, 2)) {
			fprintf(stderr, "\nError in wgtrc_xc: writing a2p1s!!\n\n");
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
			delete [] imh_max0;  
			imh_max0 = NULL;
 			return -2;
		}
		
		strcpy(fa0p1ss, "wgtrv_a1p2s");
		if (wsfunc_3phs(fa0p1ss, hks, 3)) {
			fprintf(stderr, "\nError in wgtrc_xc: writing a1p2s!!\n\n");
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
			delete [] imh_max0;  
			imh_max0 = NULL;
 			return -2;
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
	delete [] imh_max0;  
	imh_max0 = NULL;
	
	return 0;
}


/******************************************************************************/
/* wgtmh_xc() xcorr is a function of mh                                       */
/******************************************************************************/
int wgtmh_xc(HKS_INFO *hks)
{
	int iz, ir, izr;
	int it, itk;
	int i, k; 
	int *irv_max0, irv_max;
	int bmrv, tprv;
	int nrv2, nrv;
	float	rv_max;
//	float mh;
	float *xc_mh;
	float	*ar_0p1ss, *ar_0p1s, *ar_2p1s, *ar_1p2s;
//	double cc, c11, c12, c13, c22, c23, c33, rc12, rc13, wgt;
	double cc, c11, c12, c13, c22, c33, rc12, rc13, wgt;
	
	irv_max0 = NULL;
	xc_mh = NULL;
	ar_0p1ss = NULL;
	ar_0p1s = NULL;
	ar_2p1s = NULL;
	ar_1p2s = NULL;
	
/*----------------------------------------------------------------------------*/
/* nrv2 is half range in points for search and summation for sfunctions
	along 'K' of a given 'H'                                                   */
/*----------------------------------------------------------------------------*/
	nrv2 = (int)rint(_RCC_L_ * 0.5 / hks->rvdt);
	if (nrv2 < 1 ) {
		fprintf(stderr, "Error in wgtmh_xc: nrv2 = %3d, %8.5f\n", nrv2, hks->rvdt);
		return -1;
	}	
   
/*----------------------------------------------------------------------------*/
/* Assign memory for the locations of sfunction peaks for each 'H'            */
/*----------------------------------------------------------------------------*/
	irv_max0 = new int[hks->nmh];  
	if (irv_max0 == NULL) {
		fprintf(stderr, "wgtmh_xc: Error in new int irv_max0!\n\n");
		return -1;
	}	

/*----------------------------------------------------------------------------*/
/* For each 'H', find the sfunction peak in the 0p1s mode
	Store their location in irv_max0                                           */
/*----------------------------------------------------------------------------*/
	for (iz = 0; iz < hks->nmh; iz++) {
		rv_max = -MAXFLOAT;
		for (ir = 0; ir < hks->nrv; ir++) {
			izr = iz * hks->nrv + ir;
			if (rv_max < hks->a_0p1s[izr]) {
				rv_max = hks->a_0p1s[izr];
				irv_max0[iz] = ir;
			}
		}
	}

/*----------------------------------------------------------------------------*/
/* Assign memory for cross-correlation of all 'H's (xc_mh),
	sfunction amplitude of all 'K's (ar_0p1ss, ar_0p1s, ar_2p1s, ar_1p2s)
	for each 'H'                                                               */
/*----------------------------------------------------------------------------*/  
	if (xc_mh != NULL) { 
		delete [] xc_mh;
		xc_mh = NULL;
	}	
	xc_mh = new float[hks->nmh];
	if (xc_mh == NULL) {
		fprintf(stderr, "wgtmh_xc: Error in new float xc_mh!\n\n");
		return -1;
	}	
	
	ar_0p1ss = new float[hks->nrv];
	if (ar_0p1ss == NULL) {
		fprintf(stderr, "wgtmh_xc: Error in new float ar_0p1ss!\n\n");
		return -1;
	}	
	ar_0p1s = new float[hks->nrv];
	if (ar_0p1s == NULL) {
		fprintf(stderr, "wgtmh_xc: Error in new float ar_0p1s!\n\n");
		return -1;
	}	
	ar_2p1s = new float[hks->nrv];
	if (ar_2p1s == NULL) {
		fprintf(stderr, "wgtmh_xc: Error in new float ar_2p1s!\n\n");
		return -1;
	}	
	ar_1p2s = new float[hks->nrv];
	if (ar_1p2s == NULL) {
		fprintf(stderr, "wgtmh_xc: Error in new float ar_1p2s!\n\n");
		return -1;
	}	
	
/*----------------------------------------------------------------------------*/
/* Start the loop to compute cross correlation for each 'H'                   */
/*----------------------------------------------------------------------------*/	
//	mh = hks->mhlb; 
//	hks->xc_max = 0.0;
	for (iz = 0; iz < hks->nmh; iz++) {

/*----------------------------------------------------------------------------*/
/* For a given 'H', copy sfunctions of all 'K' to ar_0p1s, ar_2p1s, 
	and ar_1p2s. 
	Note in each loop of 'H', they take different sfunction values             */
/*----------------------------------------------------------------------------*/
		for (ir = 0; ir < hks->nrv; ir++) {
			izr = iz * hks->nrv + ir;
   	   ar_0p1s[ir] = hks->a_0p1s[izr];
      	ar_2p1s[ir] = hks->a_2p1s[izr];
			ar_1p2s[ir] = hks->a_1p2s[izr];
			ar_0p1ss[ir] = 0.0;
		}
		
/*----------------------------------------------------------------------------*/
/* For each 'H', the location of sfunction peak in the 0p1s mode is stored 
	in irv_max0                                           
	For a given 'H' loop, do a 5-points smoothing for irv_max0
	The average is stored in irv_max, which is different in each 'H' loop.     */
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
/* For each 'H', generate ar_0p1ss, which has non-zero values centered 
	at irv_max with the range of +/-nrv2. 
	The value of each ar_0p1ss is summation over _CD_FACTOR_ points.
	For example:
	irv_max = 9
	nrv2 = 3
	_CD_FACTOR = 4
	bmrv = irv_max - nrv2 = 9 - 3 = 6
	tprv = irv_max + nrv2 = 9 + 3 = 12
	for (i = bmrv; i <= tprv; i++) 
		it = (i - irv_max) * _CD_FACTOR_ + irv_max
		i = 6, it = (6 - 9) * 4 + 9  = -3
		i = 7, it = (7 - 9) * 4 + 9  =  1
		i = 8, it = (8 - 9) * 4 + 9  =  5
		i = 9, it = (9 - 9) * 4 + 9  =  9
		i =10, it = (10 - 9) * 4 + 9 = 13
		i =11, it = (11 - 9) * 4 + 9 = 17
		i =12, it = (12 - 9) * 4 + 9 = 21
	The it(i) = it(i-1) + _CD_FACTOR, define the starting point of summation
	ar_0p1ss[i = 6] = 0, because it < 0
	ar_0p1ss[i = 7] = sum(ar_0p1s[1] ... ar_0p1s[4])
	ar_0p1ss[i = 8] = sum(ar_0p1s[5] ... ar_0p1s[8])
	ar_0p1ss[i = 9] = sum(ar_0p1s[9] ... ar_0p1s[12])
	ar_0p1ss[i =10] = sum(ar_0p1s[13] ... ar_0p1s[16])
	ar_0p1ss[i =11] = sum(ar_0p1s[17] ... ar_0p1s[20])
	ar_0p1ss[i =12] = sum(ar_0p1s[21] ... ar_0p1s[24])
	Such that, the non-zero values of ar_0p1ss occupy 
	[irv_max-nrv2 : irv_max+nrv2], but each value is the sum of _CD_FACTOR_ 
	points all over the entire sfunction along 'K'. See ar_0p1ss[irv_max = 9]
	to understand the layout of ar_0p1ss
	Normally wgt = 1                                                           */
/*----------------------------------------------------------------------------*/
		bmrv = IMAX(irv_max - nrv2, 0);
		tprv = IMIN(irv_max + nrv2, hks->nrv-1);
		nrv = tprv - bmrv + 1;
		wgt = (double) nrv / (double)(nrv2 * 2 + 1);
		
		for (i = bmrv; i <= tprv; i++) {
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

/*----------------------------------------------------------------------------*/
/* For each 'H', calculate all cross correlations 
	Note the non-value range of ar_0p1ss   
	Obtain a cross correlation coefficient xc_mh[iz] for each 'H'            */
/*----------------------------------------------------------------------------*/
		c11 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_0p1ss[bmrv]);
		c22 = _acc_(nrv, &ar_2p1s[bmrv],  &ar_2p1s[bmrv]);
		c33 = _acc_(nrv, &ar_1p2s[bmrv],  &ar_1p2s[bmrv]);
		c12 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_2p1s[bmrv]);
		c13 = _acc_(nrv, &ar_0p1ss[bmrv], &ar_1p2s[bmrv]);
//		c23 = _acc_(nrv, &ar_2p1s[bmrv],  &ar_1p2s[bmrv]);	

		rc12 = c12/sqrt(c11 * c22); 
		rc13 = c13/sqrt(c11 * c33);
		switch(hks->xc_flag) {
			case 4:
//				xc_mh[ir] = (c12/sqrt(c11 * c22) + c13/sqrt(c11 * c33) + c23/sqrt(c22 * c33))/3. * wgt; 
				xc_mh[iz] = (rc12 + rc13) * 0.5 * wgt;
				break;
			case 5:
				xc_mh[iz] = rc12 * wgt;
				break;
			case 6:
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
/* After obtain a cross correlation coefficient xc_mh[iz] for each 'H',
	conduct a new 'H' loop to smooth xc_mh over 5 xc_mh along 'H'          */
/*----------------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------------*/
/* For a given 'H', multiple xc_mh[iz] on sfunction of the three phases 
	along 'K' for the given 'H'.
	For different 'H', xc_mh[iz] is different                                */
/*----------------------------------------------------------------------------*/
		for (ir = 0; ir < hks->nrv; ir++) {
			izr = iz * hks->nrv + ir;
			hks->a_0p1s[izr] *= cc;
			hks->a_2p1s[izr] *= cc;
			hks->a_1p2s[izr] *= cc;
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
	delete [] xc_mh;
	xc_mh = NULL;
	delete [] irv_max0;  
	irv_max0 = NULL;
	
	return 0;
}


/******************************************************************************/
/* _acc_()                                                                    */
/******************************************************************************/
static double _acc_(int n, float* data1, float* data2)
{
	int  i;
	double  sum;

	sum = 0.;
	for (i = 0; i < n; i++) {
		sum += (data1[i] * data2[i]);
	}
	return(sum);
}




