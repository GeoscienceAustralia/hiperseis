#include "recordrf.h"
#include "sac.h"

void help();
int readline(FILE* fp, char* oneline);
int segchar(char* inchar, char* charp1, char* charp2, char token, int seq);
int arrivalps(float &vp, float &rv, float &mh, float &w_para, float &t1, float &t2, float &t3);
//int rsac(char *kname, float *sacd, int& npts, float& B, float& delta, int npts_max);
//int rsac0(char *kname, float* data, int& npts, float& delta, float& b);

/************************* Main *******************************/
int main(int argc, char*argv[])
{
	FILE* contrpt;
	FILE* wfpt;
	float vp, rv, mh; 
	float mhlb, mhub, mhdt;
	float t1, t2, t3;
	double amp1, amp2, amp3;
	double sumamp1, sumamp2, sumamp3, sumamp;
	float w_para[MAX_REC];
	float B;
	float delta;
	float tmpx1, tmpx2, tmpx;
	float sacData[MAX_SACLEN];
	//float* sacData;
	float sfunc_max;
	float sfunc_mean;
	float vp_max;
	float* mhfunc;
	float* sfunc;
	float sigma_s, sigma_rv, sigma_mh;
	float max_rv, max_mh;
	float rho_rv, rho_mh, rho_rvmh;
	float magl, magld, a, a2, b, b2;
	char rf_dir[128];
	char rf_file[128];
	char rf_name[MAX_REC][64];
	char fileout[128];
	char filesum[128];
	char oneline[256];
	int i, j, k;
	int irec;
	int rec_count;
	int npts;
	int rl;
	int ntmpx;
	int nvp, nrv, nmh;
	int irv_max, imh_max;
	int nroot = 4;
	RECORD* rfrec;	
	SAC_HEADER* sachd;

	if (argc == 2) {
	  	contrpt = fopen(argv[1], "r");
  		if (contrpt == NULL) {
    		printf("\nError, can not open control file %s!!\n\n",argv[1]);
    		return 0;
  		}
		
		fscanf(contrpt, "%f", &vp);
		fscanf(contrpt, "%f", &rv);
		fscanf(contrpt, "%f %f %f", &mhlb, &mhub, &mhdt);
		fscanf(contrpt, "%s", fileout);
		fscanf(contrpt, "%s", rf_dir);		
			
		rec_count = 0;
		while(1) {
			rl = readline(contrpt, oneline);
			//printf("rl = %d, %s\n", rl, oneline);
			if  (rl == 0) { 
				break;
			}
			else if (rl == 1) {
				sscanf(oneline, "%s %f", rf_name[rec_count], &w_para[rec_count]);
				rec_count++;
				if (rec_count > MAX_REC) {
					printf("\nNo. of receiver functions = %d > Max_Rec(%d), increase buffer.\n\n",
						 rec_count, MAX_REC);
					fclose(contrpt);
					return 1;
				}
			}    
		}
	}
	else {
		help();
		return 1;
	}	
	
/*------------------------- Set up rf records ---------------------------------*/	
	printf("Reading %d receiver functions\n", rec_count);

	rfrec = new RECORD[rec_count];
	for (irec = 0; irec < rec_count; irec++) {
		rfrec[irec].InitRecord(rf_name[irec]);
		
		strcpy(rf_file, rf_dir);
		if (rf_dir[strlen(rf_dir)-1] != '/') strcat(rf_file, "/");		
		strcat(rf_file, rf_name[irec]);

		/*--------- openning input sac file 1 ------------*/
		if (!rsac0(rf_file, sacData, npts, delta, B)) {
   		printf("\nError in reading sac data!!\n\n");
   		return(0);
 		}
		
		printf("%s, B=%f, DT=%f, N=%d, sac[0]=%f\n", rfrec[irec].GetRecName(),B,delta,npts,
		sacData[0]);
				printf("\n\nHere...........\n");

		rfrec[irec].InstallData(w_para[irec], B, delta, npts, sacData);						
		printf("%s, B=%f, DT=%f, N=%d\n", rfrec[irec].GetRecName(),B,delta,npts);
	
		if (irec == 0) {
				sachd = new SAC_HEADER;
				sachd->rsacheader(rf_file);
				sachd->set_unusedf(13, 4, vp);
				sachd->set_unusedf(13, 5, rv);						
		}
	}
	
/*--------------------- prepare output file ----------------------------*/
// 	wfpt = fopen(fileout, "w");
// 	if (wfpt == NULL) {
//   	printf("\nError, can not open output file %s!!\n\n", fileout);
//   	return 0;
// 	}
//	fprintf(wfpt, "%% Vp/Vs ratio = %4.2f\n", rv);
//	fprintf(wfpt, "%% Average crust Vp = %4.2f km/s\n", vp);
//	fprintf(wfpt, "%% Moho depth km from %4.2f to %4.2f at step %4.2f\n", mhlb, mhub, mhdt);
//	fprintf(wfpt, "\n");

/*------------------------- Processing ---------------------------------*/	
	tmpx = (mhub-mhlb)/mhdt;
	printf("(mhub-mhlb)/mhdt = %f\n", tmpx);
	ntmpx = (int) tmpx;
	if (tmpx - ntmpx >= 0.1) nmh = ntmpx+2;
	else nmh = ntmpx+1;

	mhfunc = new float[nmh];
	sfunc = new float[nmh];
	mh = mhlb;
	for (j = 0; j < nmh; j++) {			
		printf("moho = %f, vp/vs = %f\n", mh, rv);
			
		sumamp1 = 0.0;
		for (irec = 0; irec < rec_count; irec++) {
			arrivalps(vp, rv, mh, w_para[irec], t1, t2, t3);
			amp1 = rfrec[irec].WaveScan(t1);
			if (amp1 != 0.0) amp1 = amp1/fabs(amp1) * pow(fabs(amp1), 1.0/nroot);

			sumamp1 += amp1;
		}	
		sumamp1 /= rec_count;
				
		if (sumamp1 != 0) sumamp1 = sumamp1/fabs(sumamp1) * pow(fabs(sumamp1), (double)nroot);
		
		mhfunc[j] = mh;	
		//sumamp = wgt1*sumamp1 + wgt2*sumamp2 - wgt3*sumamp3;
		sfunc[j] = sumamp1;
		
		//fprintf(wfpt, "%5.1f  %9.2e\n", mhfunc[j], sfunc[j]);
		
		mh += mhdt;
	}
	
	sachd->set_b(mhlb);
	sachd->set_delta(mhdt);
	sachd->set_npts(nmh);
	sachd->set_e((nmh-1)*mhdt+mhlb);
	printf("New header: B = %f, delta = %f, E = %f, npts = %d\n",
		sachd->_b, sachd->_delta, sachd->_e, sachd->_npts); 
	wsac0(fileout, sfunc, sachd);
	
				
	return (1);				
}

/*------------------------------------------------------------------*/
void help()
{
	printf("Usage:\n");
	printf("moho para_file.\n");
	printf("1. Average crust P-wave velocity\n");
	printf("2. Vp/Vs ratio: low bound, up bound, and step\n");
	printf("3. moho depth: low bound, up bound, and step\n");
	printf("4. weighting factor for t_Ps, t_PpPs, and t_PpSs+PsPs\n");
	printf("5. Out put file name for s-function\n");
	printf("6. Out put file name for results\n");
	printf("7. Dir of Input receiver function files\n");
	printf("8-end. receiver function file names and ray parameters\n");	
	printf("\n");
}

/*------------------------------------------------------------------*/
int arrivalps(float &vp, float &rv, float &mh, float &p, float &t1, float &t2, float &t3)
{
	float vs;
	float slow_s, slow_p;
	
	vs = vp/rv;
	slow_s = sqrt(1/vs/vs - p*p);
	slow_p = sqrt(1/vp/vp - p*p);
	
	t1 = mh * (slow_s - slow_p);
	t2 = mh * (slow_s + slow_p);
	t3 = mh * slow_s * 2;
	
	return 1;
}	

/*------------------------------------------------------------------*/
		 	
