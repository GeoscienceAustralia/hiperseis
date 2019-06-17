#include "recordrf.h"
#include "nr.h"
#include "matrix.h"
#include "nrutil.h"

extern void jacobi(float **a, int n, float d[], float **v, int *nrot);
extern int segchar(char* inchar, char* charp1, char* charp2, char token, int seq);


/******************************************************************************/
/* Write solution of H-k search
	Write summary of H-k search hkr2dx7.xxx, and hkr2dx7.xxx.00, 01, ... 
	This version is for non-bootstrap                                          */
/******************************************************************************/
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
	fprintf(wfpt2, "Process %d RF traces\n", hks->nrfrec);			
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
			printf("%d: Neglect, H = %f+/-%f, k = %f+%f, sfunc = %e\n", 
				i, hks->msf[i].fmh, hks->msf[i].sigma_mh, hks->msf[i].frv, hks->msf[i].sigma_rv, hks->msf[i].sfunc);			
			continue;
		}
		else {
			printf("%d:    keep, H = %f+/-%f, k = %f+%f, sfunc = %e\n", 
				i, hks->msf[i].fmh, hks->msf[i].sigma_mh, hks->msf[i].frv, hks->msf[i].sigma_rv, hks->msf[i].sfunc);			
		}
					
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
		fprintf(wfpt, "%%Process %d RF traces\n", hks->nrfrec);			
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
			
			fprintf(wfpt, "%s   %9.5f   %5.2f %7.4f   %5.2f %7.4f   %5.2f %7.4f   (%5.2f - %6.2f)\n",  
						tmpch, rfrec[irec].Ray(),
						rfrec[irec].t_0p1s, rfrec[irec].a_0p1s,
						rfrec[irec].t_2p1s, rfrec[irec].a_2p1s,
						rfrec[irec].t_1p2s, rfrec[irec].a_1p2s,
						rfrec[irec].Rfstart(), rfrec[irec].Rfend());
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
/* Write solution of H-k search
	Write summary of H-k hk2d.HL.AXX.x7, 
	and hk2d.HL.AXX.x7.b010.00, hk2d.HL.AXX.x7.b010.01, ...
	hk2d.HL.AXX.x7.020.00, hk2d.HL.AXX.x7.b020.01, ... 
	This version is moving BAZ                                                 */
/******************************************************************************/
int wmh_rvb(char* fpth, HKS_INFO *hks, RECORD *rfrec, BAZ_INFO *bazinfo)
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
	
/*----------------------------------------------------------------------------*/
/* Write summary of H-k search to hk2d.HL.AXX.x7                                 */
/*----------------------------------------------------------------------------*/
	wfpt2 = fopen(fpth, "a");
	if (wfpt2 == NULL) {
 		fprintf(stderr, "\nError in wmh_rvb: Cannot open output file %s!!\n\n", fpth);
 		return -1;
	}
	
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
			printf("%d: Neglect, H = %f+/-%f, k = %f+%f, sfunc = %e\n", 
				i, hks->msf[i].fmh, hks->msf[i].sigma_mh, hks->msf[i].frv, hks->msf[i].sigma_rv, hks->msf[i].sfunc);			
			continue;
		}
		else {
			printf("%d:    keep, H = %f+/-%f, k = %f+%f, sfunc = %e\n", 
				i, hks->msf[i].fmh, hks->msf[i].sigma_mh, hks->msf[i].frv, hks->msf[i].sigma_rv, hks->msf[i].sfunc);			
		}
					
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
/* Write summary of H-k search to each hk2d.HL.AXX.x7.b010.00, 01, ...        */
/*----------------------------------------------------------------------------*/
		sprintf(fpth_hk1, "%s.b%03d.%02d", fpth, (int)bazinfo->lb, ix_count);
		 
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
		fprintf(wfpt, "%%Process %d RF traces\n", hks->nrfrec);			
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
			
			fprintf(wfpt, "%s   %9.5f   %5.2f %7.4f   %5.2f %7.4f   %5.2f %7.4f   (%5.2f - %6.2f)\n",  
						tmpch, rfrec[irec].Ray(),
						rfrec[irec].t_0p1s, rfrec[irec].a_0p1s,
						rfrec[irec].t_2p1s, rfrec[irec].a_2p1s,
						rfrec[irec].t_1p2s, rfrec[irec].a_1p2s,
						rfrec[irec].Rfstart(), rfrec[irec].Rfend());
		}	

		fclose(wfpt);
		
/*----------------------------------------------------------------------------*/
/* Write summary of H-k search to hk2d.HL.AXX.x7                                 */
/*----------------------------------------------------------------------------*/
		if (segchar(fpth_hk1, dummy, tmpch, '/', 1) != 1) return -3;

		fprintf(wfpt2, "%-19s  %5.1f  %5.1f  %3d  %d  %6.3f +/- %6.3f  %6.3f +/- %6.3f\n", 
				  tmpch, bazinfo->lb, bazinfo->ub, hks->nrfrec, ix_count, 
				  hks->msf[i].fmh, hks->msf[i].sigma_mh, hks->msf[i].frv, hks->msf[i].sigma_rv);
		
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
/* Write solution of H-k search
	Write summary of H-k search hkr2dx7.xxx, and hkr2dx7.xxx.00, 01, ... 
	This version is for bootstrap                                              */
/******************************************************************************/
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


/******************************************************************************/
/* Write solution of H-k search
	Write S-function 
	hks->wsf_flag = 0: binary file
					  = 1: ascii file                                                        */
/******************************************************************************/
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

/*----------------------------------------------------------------------------*/
/* The sequence is first loop H and next loop K                               */			
/*----------------------------------------------------------------------------*/
		for (ir=0; ir < hks->nrv; ir++) {	
			for (iz=0; iz < hks->nmh; iz++) {
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

/*----------------------------------------------------------------------------*/
/* The sequence is first loop H and next loop K   
	The sequence is to consistent old version                                  */			
/*----------------------------------------------------------------------------*/
		float rv, mh;
		for (rv=hks->rvlb, ir=0; ir<hks->nrv; ir++, rv+=hks->rvdt) {	
			for (mh=hks->mhlb, iz=0; iz<hks->nmh; iz++,mh+=hks->mhdt) {
				fprintf(wfpt, "%6.3f  %6.2f  %11.4e\n", rv, mh, hks->a3ps[iz*hks->nrv + ir]);
			}
		
		}
	
		fclose(wfpt);		
	}
		
	return 0;
}	


/******************************************************************************/
/* Write S-functions of any phase 0p1s, 2p1s or 1p2s in terms of the option
	1, 2, or 3
	Write in binary file in the working dir
	The purpose is to illustrate and debug                                     */
/******************************************************************************/
int wsfunc_3phs(char* fpth, HKS_INFO *hks, int opt)
{
	FILE* wfpt;
	int ir, iz;
		
	wfpt = fopen(fpth, "wb");
	if (wfpt == NULL) {
  		fprintf(stderr, "\nError in wsfunc_3phs: Cannot open output file %s!!\n\n", fpth);
  		return -1;
	}

/*----------------------------------------------------------------------------*/
/* The sequence is first loop H and next loop K                               */			
/*----------------------------------------------------------------------------*/
	if (opt == 1) {
		for (ir=0; ir < hks->nrv; ir++) {	
			for (iz=0; iz < hks->nmh; iz++) {
				if (fwrite(&(hks->a_0p1s[iz*hks->nrv + ir]), sizeof(float), 1, wfpt) != 1) {
					fprintf(stderr, "\nError in wsfunc_3phs: Writing 0p1s to %s!! \n\n", fpth);
					fclose(wfpt);
					return -2;
				}				
			}		
		}
	}
	
	else if (opt == 2) {
		for (ir=0; ir < hks->nrv; ir++) {	
			for (iz=0; iz < hks->nmh; iz++) {
				if (fwrite(&(hks->a_2p1s[iz*hks->nrv + ir]), sizeof(float), 1, wfpt) != 1) {
					fprintf(stderr, "\nError in wsfunc_3phs: Writing 2p1s to %s!! \n\n", fpth);
					fclose(wfpt);
					return -2;
				}				
			}		
		}		
	}
	
	else if (opt == 3) {
		for (ir=0; ir < hks->nrv; ir++) {	
			for (iz=0; iz < hks->nmh; iz++) {
				if (fwrite(&(hks->a_1p2s[iz*hks->nrv + ir]), sizeof(float), 1, wfpt) != 1) {
					fprintf(stderr, "\nError in wsfunc_3phs: Writing 1p2s to %s!! \n\n", fpth);
					fclose(wfpt);
					return -2;
				}				
			}		
		}		
	}
	
	else {
		fprintf(stderr, "\nError in wsfunc_3phs: Unknown option = %d!! \n\n", opt);
		fclose(wfpt);
	}					

	fclose(wfpt);
	return 0;
}
	
