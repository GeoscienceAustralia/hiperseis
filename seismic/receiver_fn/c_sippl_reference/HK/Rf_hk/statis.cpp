#include "recordrf.h"


#ifndef MAX_MSFUNC
#define MAX_MSFUNC 100
#endif


/******************************************************************************/
/* statistics                                                                 */
/******************************************************************************/
int statis(HKS_INFO *hks)
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
	int c_maxsfunc;
	MAX_SFUNCTIONS exsf; 
	
/*----------------------------------------------------------------------------*/
/* Calculate a3ps = w1 * a_0p1s + w2 * a_2p1s + w3 * a_1p2s
	Note a_1p2s is -a_1p2s.
	Search for the max amplitude of sfunction, along with its location in 
	[zr] (1D)                                                                  */
/*----------------------------------------------------------------------------*/	
	hks->sfunc_max = -MAXFLOAT;
	kzr = -INT_MAX;
	for (izr = 0; izr < hks->nzr; izr++) {

		hks->a3ps[izr] = hks->w_0p1s * hks->a_0p1s[izr] 
							+ hks->w_2p1s * hks->a_2p1s[izr] 
							+ hks->w_1p2s * hks->a_1p2s[izr];
		
		if (hks->sfunc_max < hks->a3ps[izr]) {
			hks->sfunc_max = hks->a3ps[izr];			
			kzr = izr;
		}
	}
	if (kzr == -INT_MAX) {
		fprintf(stderr, "Error in statis: NOT found the Max sfunc = %e\n", hks->sfunc_max);
		return -1;
	}
	else {	 		
		printf("The Max sfunc = %e before normalization.\n", hks->sfunc_max);
	}

/*----------------------------------------------------------------------------*/
/* Normalize sfunction      
	Calculate the Mean (sfunc_mean) and variance (sqrt(sfunc_sigma)) of sfunction                              
	Search for the amplitude of sfunctions > 0.6 (normalized values)
	Calcualte the number of these sfunctions > 0.6, and
	tag their position in [zr] (1D):
	hks->tag[izr > 0.6] = 1, hks->[izr else] = 0 (initial values)              */
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
		hks->sfunc_sigma += (hks->a3ps[izr] - hks->sfunc_mean) * (hks->a3ps[izr] - hks->sfunc_mean);	
//		sfunc_store[ix].sigma_s += (hks->a3ps[i*hks->nmh + j] - hks->a3ps[irv_max*hks->nmh + imh_max]) 
//					 * (hks->a3ps[i*hks->nmh + j] - hks->a3ps[irv_max*hks->nmh + imh_max]);
		hks->tag[izr] = 0;
		if (hks->a3ps[izr] > 0.6) {
			hks->tag[izr] = 1;
			hks->ntags++;
		}	
	}
	hks->sfunc_sigma /= (hks->nzr - 1);
	hks->sfunc_sigma = sqrt(hks->sfunc_sigma);

/*----------------------------------------------------------------------------*/
/* Compute correlation coefficient between H and k (not sfunction amplitude)
	for those tagged sfunction.
	hks->mhrv_cc will be used to compute the H and k error ellipse             */
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
		
		iz = izr / hks->nrv;
		ir = izr % hks->nrv;

		mhtmp = hks->mhlb + iz * hks->mhdt;
		rvtmp = hks->rvlb + ir * hks->rvdt;	

		sum_mh += mhtmp;
		sum_rv += rvtmp;
		sum_mhrv += rvtmp * mhtmp;
		sum_mh2 += mhtmp * mhtmp;
		sum_rv2 += rvtmp * rvtmp;			
	}
	hks->mhrv_cc = (float) ((sum_mhrv - sum_mh * sum_rv / hks->ntags) / 
		sqrt((sum_mh2 - sum_mh * sum_mh / hks->ntags) * (sum_rv2 - sum_rv * sum_rv / hks->ntags)));

	
/*----------------------------------------------------------------------------*/
/* Assign memory for hks->msf recording the information of all selected 
	sfunctions close to the global peak                                        */
/*----------------------------------------------------------------------------*/
	if (hks->msf != NULL)  {
		delete [] hks->msf;
		hks->msf = NULL;
	}
	hks->msf = new MAX_SFUNCTIONS[MAX_MSFUNC];
	if (hks->msf == NULL) {
		fprintf(stderr, "\nstatis: Error in assigning memory for hks->msf !!!\n\n");
		return -1;
	}	
	for (i = 0; i < MAX_MSFUNC; i++) {
		hks->msf[i].isOK = 0;
	}
		
/*----------------------------------------------------------------------------*/
/* The global peak is stored in hks->msf[0] 
	Search for all sfunction values close to the global peak
	Do not count points within pre-set Euclid distance
	The Euclid distance is defined by hks->sfcrt_mh and hks->sfcrt_rv,
	which are input parameters                                                 */ 	
/*----------------------------------------------------------------------------*/
	hks->n_maxsfunc = 1;
	hks->msf[0].sfunc = hks->sfunc_max;
	hks->msf[0].kzr = kzr;
	hks->msf[0].isOK = 1;

/*----------------------------------------------------------------------------*/
/* Compute the pre-set Euclid distance                                        */  
/*----------------------------------------------------------------------------*/
	dist2_crt = (int)((hks->sfcrt_mh/hks->mhdt/2) * (hks->sfcrt_mh/hks->mhdt/2) + 
							(hks->sfcrt_rv/hks->rvdt/2) * (hks->sfcrt_rv/hks->rvdt/2));

/*----------------------------------------------------------------------------*/
/* Search for all sfunctions                                                  */  
/*----------------------------------------------------------------------------*/
	for (izr = 0; izr < hks->nzr; izr++) {

/*----------------------------------------------------------------------------*/
/* hks->sfcrt_peak defines a percent of difference between all sfunction 
	amplitudes and the global peak
	Search for those sfunctions < the percent, and
	preliminary select them (isOK = 1), but require further check              */  
/*----------------------------------------------------------------------------*/
		if ((hks->sfunc_max - hks->a3ps[izr]) / hks->sfunc_max < hks->sfcrt_peak) {
			iz = izr / hks->nrv;
			ir = izr % hks->nrv;
			isOK = 1;

/*----------------------------------------------------------------------------*/
/* Seach existing hks->msf database recording the information of all selected 
	sfunctions close to the global peak
	if 'izr' is close to the given hks->msf peak
		if the given peak is the global peak hks->msf[i = 0], ignore it
		if not, but 'izr' has higher amplitude than the given hks->msf peak
			replace the one in hks->msf database with the new one                */
/*----------------------------------------------------------------------------*/
			for (i = 0; i < hks->n_maxsfunc; i++) {
				irv_max = hks->msf[i].kzr % hks->nrv;
				imh_max = hks->msf[i].kzr / hks->nrv;										
				dist2 = (ir - irv_max) * (ir - irv_max) + (iz - imh_max) * (iz - imh_max);
				
				if (dist2 <= dist2_crt) {
					if (i == 0) {
						isOK = 0;
						break;
					}	
					else {
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
			
/*----------------------------------------------------------------------------*/
/* After search for all sfunctions
	if the sfunction is still selected, put it to hks->msf database            */  
/*----------------------------------------------------------------------------*/
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
/* Sort the selected peaks in hks->msf in decending order 
	Because hks->msf has no empty slot before n_maxsfunc, 
	all msf.isOK == 1 for [0, n_maxsfunc-1].
	Only switch sfunction amplitude (msf.sfunc) and location (msf.kzr)         */ 	
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
	
/*----------------------------------------------------------------------------*/
/* Some statistics
	Calcuate variance of sfunction peaks in hks->msf                           */	
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
/*	Since we already know the range of H and K 
	Do the last round of selection for hks->msf bu ruling out those overlaping 
	in H and K
	Note the hks->msf has blanks [0, n_maxsfunc-1]                             */
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

					if (!(mh_ub1 < mh_lb2 || mh_ub2 < mh_lb1) && !(rv_ub1 < rv_lb2 || rv_ub2 < rv_lb1)) {
						hks->msf[j].isOK = 0;
					}
				}
			}
		}
	}					
	
	c_maxsfunc = 0;
	for (i = 0; i < hks->n_maxsfunc; i++) {
		if (hks->msf[i].isOK == 1) c_maxsfunc++;
	}
	printf("Total %d/%d Max sfunctions are kept.\n", c_maxsfunc, hks->n_maxsfunc);
	
	return 0;
}		

