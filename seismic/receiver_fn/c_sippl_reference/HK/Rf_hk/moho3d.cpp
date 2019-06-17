#include "recordrf.h"
#include "sac.h"

#ifndef	EPSON
#define	EPSON		0.01
#endif

#ifndef MAX_REC
#define MAX_REC 1024
#endif

#ifndef MAX_MSFUNC
#define MAX_MSFUNC 300
#endif

#ifndef DEFAULT_W0P1S
#define DEFAULT_W0P1S 0.7
#endif

#ifndef DEFAULT_W2P1S
#define DEFAULT_W2P1S 0.2
#endif

#ifndef DEFAULT_W1P2S
#define DEFAULT_W1P2S 0.1
#endif

#ifndef DEFAULT_TWL
#define DEFAULT_TWL 0.1   /*sec*/
#endif

#ifndef DEFAULT_NSTACK
#define DEFAULT_NSTACK 1;
#endif

#ifndef	_ZCC_L_
#define	_ZCC_L_		5.     	/* in km	*/
#endif

#ifndef	_CD_FACTOR_
#define	_CD_FACTOR_	4	
#endif

#ifndef	SIGN
#define	SIGN(x)	((x) > 0 ? 1.: -1.)
#endif


void help();
static int set_n(float minVal, float &maxVal, float step);
int readline(FILE* fp, char* oneline);
int segchar(char* inchar, char* charp1, char* charp2, char token, int seq);
void _errmsg_(char *msg);
int _ss2token_(char *str, char **av, char token);
static double _acc_(int n, float* data1, float* data2);
static int wgt_xc(HKS_INFO *hks);
static int statis(HKS_INFO *hks, int sopt);


/************************* Main *******************************/
int main(int argc, char*argv[])
{
	FILE* contrpt;
	FILE* wfpt;
	float vp, vp_max, sfunc_gmax; 
	float w_para;
	float t_sht;
	double tmpx1, tmpx2, tmpx;
	double sa_0p1s, sa_2p1s, sa_1p2s;
	char *flist = NULL;
	char rf_dir[128];
	char rf_file[128];
	char rf_name[64];
	char *fodir = NULL;
//	char *fsfunc = NULL;
	char *fhk = NULL;
//	char fpth_sfunc[128];
	char fpth_hk[128];
	char seg1[16];
	char seg2[8];
	char oneline[256];
  char *av[64];  
	int i, j, k;
	int ivp;
	int iz, ir, izr;
	int irec;
	int rl;
	int c;
	int ac;
	int sline;
	int tsht_flag;
	RECORD *rfrec;
	HKS_INFO *hks;
	
	hks = new HKS_INFO;
	hks->vplb = 0.0;
	hks->vpub = 0.0;
	hks->vpdt = 0.0;
	hks->rvlb = 0.0;
	hks->rvub = 0.0;
	hks->rvdt = 0.0;
	hks->mhlb = 0.0;
	hks->mhub = 0.0;
	hks->mhdt = 0.0;	
	hks->w_0p1s = DEFAULT_W0P1S;
	hks->w_2p1s = DEFAULT_W2P1S; 
	hks->w_1p2s = DEFAULT_W1P2S;
	hks->twl = DEFAULT_TWL;
	hks->nroot = DEFAULT_NSTACK;
	hks->xc_flag = 0;
	hks->msf = NULL;
	hks->a_0p1s = NULL;
	hks->a_2p1s = NULL;
	hks->a_1p2s = NULL;
	hks->a3ps = NULL;
	hks->xc = NULL;
	tsht_flag = 1;

	while ((c=getopt(argc, argv, "P:Z:K:W:T:N:X:L:H:D:O:s")) != -1) {
		switch (c) {
			case 'P':
	    	ac = _ss2token_(optarg, av, '/');
	    	if (ac != 3) _errmsg_(optarg);
				hks->vplb = (float)atof(av[0]);
				hks->vpub = (float)atof(av[1]);
				hks->vpdt = (float)atof(av[2]);
				break;
			case 'Z':
	    	ac = _ss2token_(optarg, av, '/');
	    	if (ac != 3) _errmsg_(optarg);
				hks->mhlb = (float)atof(av[0]);
				hks->mhub = (float)atof(av[1]);
				hks->mhdt = (float)atof(av[2]);
				break;
			case 'K':
	    	ac = _ss2token_(optarg, av, '/');
	    	if (ac != 3) _errmsg_(optarg);
				hks->rvlb = (float)atof(av[0]);
				hks->rvub = (float)atof(av[1]);
				hks->rvdt = (float)atof(av[2]);
				break;
			case 'W':
	    	ac = _ss2token_(optarg, av, '/');
	    	if (ac != 3) _errmsg_(optarg);
				hks->w_0p1s = (float)atof(av[0]);
				hks->w_2p1s = (float)atof(av[1]);
				hks->w_1p2s = (float)atof(av[2]);
				break;
			case 'T':
	    	hks->twl = (float)atof(optarg);
				break;
			case 'N':
	    	hks->nroot = atoi(optarg);
				break;
			case 'X':
	    	hks->xc_flag = atoi(optarg);
				if (hks->xc_flag < 0 || hks->xc_flag >= 4) {
					fprintf(stderr, "Possible argument (== %d?) for option -%c is [0, 1, 2, 3].\n", hks->xc_flag, optopt);
					help();
					return -1;
				}				 				
				break;
			case 'L':
	    	flist = optarg;
	    	break;
			case 'H':
	    	sline = atoi(optarg);
				if (sline < 1) {
					fprintf(stderr, "Possible argument (== %d?) for option -%c is [1, 2, 3].\n", sline, optopt);
					help();
					return -1;
				}				 
	    	break;
			case 'D':
	    	fodir = optarg;
	    	break;
			case 'O':
				fhk = optarg;
	    	break;
      case '?':
        if (optopt == 'W') {
					hks->w_0p1s = DEFAULT_W0P1S;
					hks->w_2p1s = DEFAULT_W2P1S; 
					hks->w_1p2s = DEFAULT_W1P2S;
				}
				else if (optopt == 'N') {
					hks->nroot = DEFAULT_NSTACK;
				}
				else if (optopt == 'D') {
		    	fodir = new char[8];
					strcpy(fodir, "./");
				}
				else if (optopt == 'P' || optopt == 'Z' || optopt == 'K' || optopt == 'L' || optopt == 'H' || optopt == 'O') { 	
					fprintf(stderr, "Option -%c requires an argument.\n", optopt);
					help();
					return -1;
				}
        else if (isprint(optopt)) {
					fprintf(stderr, "Unknown option `-%c'.\n", optopt);
					help();
					return -1;
				}	
        else {
					fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
					help();
					return -1;
				}	
      default:
				help();
				return -1;
		}
	}
   
	if (sline == 0) {
		fprintf(stderr, "Possible argument (== %d?) for option -%c is [1, 2, 3].\n", sline, optopt);
		help();
		return -1;
	}
		 
	printf("\nVp_search = %f : %f :%f\n",	hks->vplb, hks->vpdt, hks->vpub);
	printf("H_search = %f : %f : %f\n",	hks->mhlb, hks->mhdt, hks->mhub);
	if (hks->mhlb == 0.0 || hks->mhdt == 0.0 || hks->mhub == 0.0) return -1;
	printf("K_search = %f : %f : %f\n",	hks->rvlb, hks->rvdt, hks->rvub);
	if (hks->rvlb == 0.0 || hks->rvdt == 0.0 || hks->rvub == 0.0) return -1;	
	printf("Weighting = %f, %f, %f\n",	hks->w_0p1s, hks->w_2p1s, hks->w_1p2s);
	printf("Smoothing time window = %f\n", hks->twl);
	printf("N-th root = %d\n", hks->nroot);
	printf("Cross-correlation weighting flag = %d\n", hks->xc_flag);
	
	if (flist == NULL) {
		fprintf(stderr, "Time trace list file = NULL\n");
		return -1;
	}
	else {	
		printf("Time trace list file = '%s'\n", flist);
	}
		
	if (fodir == NULL) {
	 	fodir = new char[8];
		strcpy(fodir, "./");
	}		
	printf("Output dir = '%s'\n", fodir);
	
//	if (fsfunc == NULL) {
//		fprintf(stderr, "Output file: sfunc = NULL\n");
//		return -1;
//	}
//	else {	
//		printf("Output file: sfunc = '%s'\n", fsfunc);
//	}
		
	if (fhk == NULL) {
		fprintf(stderr, "Output file: mh_rv = NULL\n");
		return -1;
	}
	else {	
		printf("Output file: mh_rv = '%s'\n", fhk);
	}


/*------------------------- Reading rf list ---------------------------------*/	
/*------------------------- Only find the number of rf records ---------------------------------*/	
 	contrpt = fopen(flist, "r");
	if (contrpt == NULL) {
  	fprintf(stderr, "\nError, can not open control file %s!!\n\n",flist);
   	return -1;
 	}

	for (i = 1; i < sline; i++)  {
		rl = readline(contrpt, oneline);
		if (rl == 0) { 
			fprintf(stderr, "Imcomplete file: %s.\n", flist);
			return -1;
		}
	}	
		
	rl = readline(contrpt, oneline);
	if (rl == 0) { 
		fprintf(stderr, "Imcomplete file: %s.\n", flist);
		return -1;
	}
		
	hks->nrfrec = 0;
	while(1) {
		rl = readline(contrpt, oneline);
		if (rl == 0) break;
		else if (rl == 1) hks->nrfrec++;  
	}
	
	printf("Reading %d receiver functions\n", hks->nrfrec);
	rfrec = new RECORD[hks->nrfrec];

/*------------------------- Processing ---------------------------------*/	
/*---------- Set search loop number -------------------*/	
	hks->nvp = set_n(hks->vplb, hks->vpub, hks->vpdt);
	hks->nrv = set_n(hks->rvlb, hks->rvub, hks->rvdt);
	hks->nmh = set_n(hks->mhlb, hks->mhub, hks->mhdt);
	printf("Vp_search = %f : %f : %f(reset), n = %d\n",	hks->vplb, hks->vpdt, hks->vpub, hks->nvp);
	printf("H_search = %f : %f : %f(reset), n = %d\n",	hks->mhlb, hks->mhdt, hks->mhub, hks->nmh);
	printf("K_search = %f : %f : %f(reset), n = %d\n",	hks->rvlb, hks->rvdt, hks->rvub, hks->nrv);
	hks->nzr = hks->nrv * hks->nmh;

	tmpx = hks->w_0p1s + hks->w_2p1s + hks->w_1p2s;
	hks->w_0p1s /= tmpx;
	hks->w_2p1s /= tmpx;
	hks->w_1p2s /= tmpx;
	printf("Weighting = %f, %f, %f(reset)\n",	hks->w_0p1s, hks->w_2p1s, hks->w_1p2s);

	
/*------------------------- Reading rf list ---------------------------------*/	
/*------------------------- Set up rf records ---------------------------------*/	
 	rewind(contrpt);
	for (i = 1; i < sline; i++) rl = readline(contrpt, oneline);
		
	rl = readline(contrpt, oneline);
	sscanf(oneline, "%s", rf_dir);
	printf("Dir of sacfiles = '%s'\n", rf_dir);
	
	irec = 0;	
	while(1) {
		rl = readline(contrpt, oneline);
		if (rl == 0) break;
		else if (rl == 1) {
			sscanf(oneline, "%s%f%f%*[^\n]", rf_name, &w_para, &t_sht);
		
			strcpy(rf_file, rf_dir);
			if (rf_dir[strlen(rf_dir)-1] != '/') strcat(rf_file, "/");		
			strcat(rf_file, rf_name);

			if (tsht_flag == 0) {
				//printf("Unset time shift.\n");
				t_sht = 0.0;
			}	
			rfrec[irec].Init(rf_file, w_para, t_sht);
			rfrec[irec].LoadData();					
			//printf("%d: '%s'\n", irec, rfrec[irec].GetName());
			//printf("%d: %f %f %f\n", irec, rfrec[irec].ray, rfrec[irec].Delta(), rfrec[irec].Bt());
			 
			if (irec > 0 && rfrec[irec].Delta() != rfrec[irec-1].Delta()) {
				fprintf(stderr, "\nRec %d: Inconsistent delta of '%s'(%f) to (%f)\n\n", irec, rfrec[irec].GetName(), 
					rfrec[irec].Delta(), rfrec[irec-1].Delta());
				fclose(contrpt);					
				return -1;
			}
			irec++;

		}    
	}
	fclose(contrpt);
	
/*---------- Set points of smoothing time window -------------------*/	
  float	rbuf;
  hks->delta = rfrec[0].Delta();
  rbuf = hks->twl / hks->delta;
  if ((rbuf - (int)rbuf) > EPSON * hks->delta) {
		hks->ntwl = (int)ceil(rbuf) + 1;
	}
	else{
		hks->ntwl = (int)rbuf + 1;
	}


/*---------- Prepare output file ---------------------------------*/
/*---------- Out put into mh_rv2d.XXX ----------------------------*/	
	if (segchar(fhk, seg1, seg2, '.', 1) != 1) return 0;
	seg1[strlen(seg1)-1] = '\0';
	if (fodir[strlen(fodir)-1] == '/') sprintf(fpth_hk, "%s%sx%d.%s", fodir, seg1, hks->xc_flag, seg2);	
	else sprintf(fpth_hk, "%s/%sx%d.%s", fodir, seg1, hks->xc_flag, seg2);	

 	wfpt = fopen(fpth_hk, "w");
 	if (wfpt == NULL) {
   	printf("\nError, can not open output file %s!!\n\n", fpth_hk);
   	return 0;
 	}
	fprintf(wfpt, "%%%4s  %15s  %15s  %9s\n", "Vp", "K", "Moho", "sfunc");


	hks->a3ps = new float[hks->nzr];
	hks->a_0p1s = new float[hks->nzr];
	hks->a_2p1s = new float[hks->nzr];
	hks->a_1p2s = new float[hks->nzr];
	sfunc_gmax = -MAXFLOAT;
	vp = hks->vplb;
/*------------- Search Vp ----------------------------------------*/	
	for (ivp = 0; ivp < hks->nvp; ivp++) { 

/*---------- Calculate travel time: ts, tp -------------------*/	
		for (irec = 0; irec < hks->nrfrec; irec++) {
			rfrec[irec].TpTs(vp, hks);   //Note how the mh and rv sort into 1-D array
		}

/*---------- Calculate travel time, pick up amplitude -------------------*/	
/*---------- Stacking all stations -------------------*/	
		if (hks->nroot == 1) {  /* linear stacking */
			for (izr = 0; izr < hks->nzr; izr++) {
				if (izr % (hks->nzr/10) == 0)	printf("%d-root Stacking: %d, iz = %d, ir = %d\n", hks->nroot, izr, izr / hks->nrv, izr % hks->nrv);
				sa_0p1s = 0.0;
				sa_2p1s = 0.0;
				sa_1p2s = 0.0;
				for (irec = 0; irec < hks->nrfrec; irec++) {

					rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array
					rfrec[irec].PeakAmp(hks);

					tmpx  = rfrec[irec].a_0p1s; 
					tmpx1 = rfrec[irec].a_2p1s; 
					tmpx2 = rfrec[irec].a_1p2s; 

					if (tmpx == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					if (tmpx1 == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					if (tmpx2 == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);

					sa_0p1s += tmpx;
					sa_2p1s += tmpx1;
					sa_1p2s += tmpx2;
				}	
				sa_0p1s /= hks->nrfrec;
				sa_2p1s /= hks->nrfrec;
				sa_1p2s /= hks->nrfrec;
	
				hks->a_0p1s[izr] = sa_0p1s;
				hks->a_2p1s[izr] = sa_2p1s;
				hks->a_1p2s[izr] = -sa_1p2s;
			}	
		}
		else {  /* nth-root stacking */
			for (izr = 0; izr < hks->nzr; izr++) {
				if (izr % (hks->nzr/10) == 0)	printf("%d-root Stacking: %d, iz = %d, ir = %d\n", hks->nroot, izr, izr / hks->nrv, izr % hks->nrv);
				sa_0p1s = 0.0;
				sa_2p1s = 0.0;
				sa_1p2s = 0.0;
				for (irec = 0; irec < hks->nrfrec; irec++) {

					rfrec[irec].Ttps(hks, izr);  //Note how the mh and rv sort into 1-D array
					rfrec[irec].PeakAmp(hks);

					tmpx  = rfrec[irec].a_0p1s; 
					tmpx1 = rfrec[irec].a_2p1s; 
					tmpx2 = rfrec[irec].a_1p2s; 

					if (tmpx == 0.0) printf("0p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					if (tmpx1 == 0.0) printf("2p1s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					if (tmpx2 == 0.0) printf("1p2s: irec = %d, izr = %d, == 0.0\n", irec, izr);
					tmpx  = SIGN(tmpx)  * pow(fabs(tmpx),  1.0/(double)hks->nroot);
					tmpx1 = SIGN(tmpx1) * pow(fabs(tmpx1), 1.0/(double)hks->nroot);
					tmpx2 = SIGN(tmpx2) * pow(fabs(tmpx2), 1.0/(double)hks->nroot);

					sa_0p1s += tmpx;
					sa_2p1s += tmpx1;
					sa_1p2s += tmpx2;
				}	
				sa_0p1s /= hks->nrfrec;
				sa_2p1s /= hks->nrfrec;
				sa_1p2s /= hks->nrfrec;
	
				tmpx  = 1.0;
				tmpx1 = 1.0;
				tmpx2 = 1.0;
				for (i = 0; i < hks->nroot; i++)	{
					tmpx  *= fabs(sa_0p1s);
					tmpx1 *= fabs(sa_2p1s);
					tmpx2 *= fabs(sa_1p2s);
				}	
				hks->a_0p1s[izr] = SIGN(sa_0p1s) * tmpx;
				hks->a_2p1s[izr] = SIGN(sa_2p1s) * tmpx1;
				hks->a_1p2s[izr] = -SIGN(sa_1p2s) * tmpx2;
			}
		}		
		 
/*---------- Cross-corelation between 0p1s, 2p1s, 1p2s -------------------------*/
		if (hks->xc_flag > 0) {
			if (wgt_xc(hks) != 1) return -1;
		}	
			
/*--------- Search for the maximum, do statistics ---------------*/	
		if (statis(hks, 1) != 1) return -1;

		for (i = 0; i < hks->n_maxsfunc; i++) {
			if (sfunc_gmax < hks->msf[i].sfunc) {
				vp_max = vp;
				sfunc_gmax = hks->msf[i].sfunc;
			} 
		}

/*--------- Out put into mh_rv2d.XXX.00 for xc wgt--------------------------------*/	
		printf("Vp = %f, Current Max Vp = %f\n", vp, vp_max);
		for (i = 0; i < hks->n_maxsfunc; i++) {
			if (hks->msf[i].isOK == 0) {
				//printf("%d: sfunc = %e, Neglect.\n", i, hks->msf[i].sfunc);			
				continue;
			}			
			fprintf(wfpt, "%5.2f  %2d  %6.3f +/- %6.3f  %6.3f +/- %6.3f  %11.4e  \n", 
							vp, i, hks->msf[i].fmh, hks->msf[i].sigma_mh, hks->msf[i].frv, hks->msf[i].sigma_rv,
							hks->msf[i].sfunc);
		}
		vp += hks->vpdt;				
	
	}
		
	fprintf(wfpt, "\n\n%%Vp loop from %5.2f to %5.2f at step %5.2f\n", hks->vplb, hks->vpub, hks->vpdt);
	fprintf(wfpt, "%%Vp/Vs ratio from %6.3f to %6.3f at step %6.3f\n", hks->rvlb, hks->rvub, hks->rvdt);
	fprintf(wfpt, "%%Moho depth from %6.3f to %6.3f at step %6.3f\n", hks->mhlb, hks->mhub, hks->mhdt);
	fprintf(wfpt, "%%Selected average crust Vp = %6.4f km/s\n", vp);
	fprintf(wfpt, "%%Weighting: t_Ps = %4.2f, t_PpPs = %4.2f, t_PpSs+PsPs = %4.2f\n", hks->w_0p1s, hks->w_2p1s, hks->w_1p2s);
	fprintf(wfpt, "%%Smoothing time window = %f\n", hks->twl);
	fprintf(wfpt, "%%N-th root = %d\n", hks->nroot);
	fprintf(wfpt, "%%Cross-correlation weighting flag = %d\n\n\n", hks->xc_flag);
	fprintf(wfpt, "%%Selected Vp = %5.2f\n", vp_max);
	
	fclose(wfpt);
	delete hks;
	delete [] rfrec;
	return (0);				
}


/*------------------------------------------------------------------*/
void help()
{
	printf("Usage:\n");
	printf("moho2d -P 5.5/6.5/0.1 -Z 10.0/100.0/0.1 -K 1.50/2.00/0.01 -W 0.6/0.2/0.2 -T 0.1 -N n \n"); 
	printf("       -L list -H 3 -s\n");
	printf("       -D dirout -O mh_rv3d \n");
	printf("Option = -P, Argument = 5.5/6.5/0.1: Average crustal velocity start/end/increment\n");
	printf("Option = -Z, Argument = 10/100/0.1: Depth search start/end/increment\n");
	printf("Option = -K, Argument = 1.5/2.0/0.01: Vp/Vs ratio search start/end/increment\n");
	printf("Option = -W, Argument = 0.6/0.2/0.2: Weighing factor for 0p1s/2p1s/1p2s. Default = 0.7/0.2/0.1\n");
	printf("Option = -T, Argument = 0.1: Smoothing time window in sec. Default = 0.1\n");
	printf("Option = -N, Argument = n: n-th stacking. Default = 1\n");
	printf("Option = -X, Argument = 0-3: 0-No XC weight, 1-All, 2-0p1s+2p1s, 3-0p1s+1p2s. Default = 1\n");
	printf("Option = -L, Argument = list: The file name of list of time traces.\n");
	printf("Option = -H, Argument = 1: Read the list from the n-th line. n = 1 for the 1st line. \n");
	printf("                           The first line that read indicates the dir of list. \n");
	printf("Option = -D, Argument = dirout: The dir for output file.\n");
	printf("Option = -O, Argument = output files: mh_rv3d.\n");
	printf("Option = -s, Unset Time shift.\n");
	printf("\n");
}


/*----------------- _set_n_ -----------------------------------*/
static int set_n(float minVal, float &maxVal, float step)
{
	int n;
  float rbuf;

  rbuf = (maxVal - minVal) / step;
  if ((rbuf - (int)rbuf) > EPSON * step) {
		n = (int)ceil(rbuf) + 1;
  }
	else{
    n = (int)rbuf  + 1;
  }
  maxVal = minVal + (n - 1) * step;
  return(n);
}
		 	
/*------------- _acc_() -------------------------*/
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

/*------------- _cc_wgt_() -------------------------*/
static int wgt_xc(HKS_INFO *hks)
{
  int	iz, ir, it, itk;
	int i, k; 
	int *imh_max0, imh_max;
	int lfmh, rhmh;
	int nmh2, nmh;
  float	mh_max;
	float rv, mh;
	float xc_max;
  float	*ar_0p1ss, *ar_0p1s, *ar_2p1s, *ar_1p2s;
	double	cc, c11, c12, c13, c22, c23, c33, rc12, rc13, wgt;
	
  nmh2 = (int)rint(_ZCC_L_ * 0.5 / hks->mhdt);
  if (nmh2 < 1 ) {
	  printf("_cc_wgt_: nmh2 = %3d %8.5f\n", nmh2, hks->mhdt);
		return -1;
	}	
   
	imh_max0 = new int[hks->nzr];  
	/* for each vp/vs ratio, find the peak in the PS mode */
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

  
	if (hks->xc == NULL) {
		//printf("hks->xc == NULL, new\n");
		hks->xc = new float[hks->nrv];
	}	
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
		
		/*- imh_max average index over max_mh at nearby rv -*/
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

		/*- set up smoothing range for mh -*/
		lfmh = MAX(imh_max - nmh2, 0);
		rhmh = MIN(imh_max + nmh2, hks->nmh-1);
		nmh = rhmh - lfmh + 1;
		wgt = (double) nmh / (double)(nmh2 * 2 + 1);
		
		for (i = lfmh; i <= rhmh; i++) {
//			printf("_cc_wgt_: ir= %3d %7.4f i_max= %3d %6.1f i1= %3d i2= %3d zN= %3d ", 
//		  	 ir, r, i_max, i_max * ds->z_inc + ds->z_min, i1, i2, zN);
	    it = (i - imh_max) * _CD_FACTOR_ + imh_max;
	    if (it >= 0 && it < hks->nmh) {
				for (k = 0; k < _CD_FACTOR_; k++) {
					itk = it + k;
			    if (itk >= hks->nmh) break;	
					ar_0p1ss[i] += ar_0p1s[itk];
				}
	/*		a_0p1ss[ii] /= (float)k; */
			}
//	    printf("ii = %3d it = %3d %10.5f\n", ii, it, a_0p1ss[ii]);
		}
		c11 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_0p1ss[lfmh]);
		c22 = _acc_(nmh, &ar_2p1s[lfmh],  &ar_2p1s[lfmh]);
		c33 = _acc_(nmh, &ar_1p2s[lfmh],  &ar_1p2s[lfmh]);
		c12 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_2p1s[lfmh]);
		c13 = _acc_(nmh, &ar_0p1ss[lfmh], &ar_1p2s[lfmh]);
		c23 = _acc_(nmh, &ar_2p1s[lfmh],  &ar_1p2s[lfmh]);	

/*#ifdef	DEBUG_DETAIL
	{
	  char	tmpCC[256];
		float	b;

	  b = i1 * ds->z_inc + ds->z_min;
	  sprintf(tmpCC, "%06.3f.0p1s", r);
	  wt_sac1(tmpCC, &a_0p1ss[i1], zN, b, ds->z_inc);
	  sprintf(tmpCC, "%06.3f.2p1s", r);
	  wt_sac1(tmpCC, &a_2p1s[i1], zN, b, ds->z_inc);
    sprintf(tmpCC, "%06.3f.1p2s", r);
    wt_sac1(tmpCC, &a_1p2s[i1], zN, b, ds->z_inc);

		sprintf(tmpCC, "%06.3f.0P1S", r);
    wt_sac1(tmpCC, &a_0p1s[0], ds->zN, ds->z_min, ds->z_inc);
		sprintf(tmpCC, "%06.3f.0P1SS", r);
    wt_sac1(tmpCC, &a_0p1ss[0], ds->zN, ds->z_min, ds->z_inc);
    sprintf(tmpCC, "%06.3f.2P1S", r);
    wt_sac1(tmpCC, &a_2p1s[0], ds->zN, ds->z_min, ds->z_inc);
    sprintf(tmpCC, "%06.3f.1P2S", r);
    wt_sac1(tmpCC, &a_1p2s[0], ds->zN, ds->z_min, ds->z_inc);
	}
#endif
*/
		rc12 = c12/sqrt(c11 * c22); 
		rc13 = c13/sqrt(c11 * c33);
		switch(hks->xc_flag) {
			case 1:
	    	// ds->cc[ir] = (c12/sqrt(c11 * c22) + 
		    //  c13/sqrt(c11 * c33) + 
		    //  c23/sqrt(c22 * c33))/3. * wgt; 
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
      cc = 0.40 * hks->xc[ir] + 0.20 * (hks->xc[ir-1] + hks->xc[ir+1]) + 0.10 * (hks->xc[ir-2] + hks->xc[ir+2]) + 1.;
    }
//    ccw[ir] = (float) cc;
    for (iz = 0; iz < hks->nmh; iz++) {
	    it = iz * hks->nrv + ir;
	    hks->a_0p1s[it] *= cc;
	    hks->a_2p1s[it] *= cc;
	    hks->a_1p2s[it] *= cc;
		}
	}
	
/*	if (savesac) {
		sprintf(ccwFile, "ccw_%01d.sac", ds->cc_flag);
				wt_sac1(ccwFile, &ccw[0], ds->rN, ds->r_min, ds->r_inc);
	}
*/
	delete [] ar_0p1ss;
	delete [] ar_0p1s;
	delete [] ar_2p1s;
	delete [] ar_1p2s;
	delete [] imh_max0;  
	
	return 1;
}


/*--------- statistics -------------------------*/
static int statis(HKS_INFO *hks, int sopt)
{
	int i, j, izr, kzr;
	float mh_lb1, mh_ub1, mh_lb2, mh_ub2; 
	float rv_lb1, rv_ub1, rv_lb2, rv_ub2;
	float tmpx1, tmpx2, tmpx;
	MAX_SFUNCTIONS exsf;
	
	hks->sfunc_max = -MAXFLOAT;
	hks->sfunc_mean = 0.0;
	for (izr = 0; izr < hks->nzr; izr++) {

		hks->a3ps[izr] = hks->w_0p1s*hks->a_0p1s[izr] + hks->w_2p1s*hks->a_2p1s[izr] + hks->w_1p2s*hks->a_1p2s[izr];
		/*if (hks->xc_flag == 0) 
			hks->a3ps[izr] = hks->w_0p1s*hks->a_0p1s[izr] + hks->w_2p1s*hks->a_2p1s[izr] + hks->w_1p2s*hks->a_1p2s[izr];
		else 
			hks->a3ps[izr] = hks->a_0p1s[izr] + hks->a_2p1s[izr] + hks->a_1p2s[izr];
		*/
		
		if (hks->sfunc_max < hks->a3ps[izr]) {
			hks->sfunc_max = hks->a3ps[izr];			
			kzr = izr;
			//irv_max = izr % hks->nrv;
			//imh_max = izr / hks->nrv;
		}
		hks->sfunc_mean += hks->a3ps[izr];
	}	 		
	hks->sfunc_mean /= hks->nzr;
	printf("The Max sfunc = %e\n", hks->sfunc_max);


//*--------- Search for all sfunction values close to the maximum -------*/ 	
	if (hks->msf == NULL) hks->msf = new MAX_SFUNCTIONS[MAX_MSFUNC];
	
	if (sopt) {
		hks->n_maxsfunc = 0;
		for (izr = 0; izr < hks->nzr; izr++) {			
			if ((hks->sfunc_max - hks->a3ps[izr]) / hks->sfunc_max < 0.009) {
				hks->msf[hks->n_maxsfunc].sfunc = hks->a3ps[izr];
				hks->msf[hks->n_maxsfunc].kzr = izr;
				hks->msf[hks->n_maxsfunc].isOK = 1;
				hks->n_maxsfunc++;
				if (hks->n_maxsfunc >= MAX_MSFUNC) {
					fprintf(stderr, "\nNo. of Max_sfunc = %d > Max, increase buffer.\n\n",hks->n_maxsfunc);
					break;
				}
			}	
		}
	}
	else {
		hks->n_maxsfunc = 1;
		hks->msf[0].sfunc = hks->sfunc_max;
		hks->msf[0].kzr = kzr;
		hks->msf[0].isOK = 1;
	}
	
/*-------------- Sort the maxs of sfunction ------------------*/ 	
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
	
/*------------------------- Some statistics ---------------------------------*/	
	for (i = 0; i < hks->n_maxsfunc; i++) {
		
		hks->msf[i].krv = hks->msf[i].kzr % hks->nrv;
		hks->msf[i].kmh = hks->msf[i].kzr / hks->nrv;
		
		hks->msf[i].frv = hks->rvlb + hks->msf[i].krv * hks->rvdt;	
		hks->msf[i].fmh = hks->mhlb + hks->msf[i].kmh * hks->mhdt;
			
		hks->sfunc_sigma = 0.0;		
		for (izr = 0; izr < hks->nzr; izr++) {			
			hks->sfunc_sigma += (hks->a3ps[izr] - hks->sfunc_mean)*(hks->a3ps[izr] - hks->sfunc_mean);	
			//sfunc_store[ix].sigma_s += (hks->a3ps[i*hks->nmh + j] - hks->a3ps[irv_max*hks->nmh + imh_max]) 
			//				 * (hks->a3ps[i*hks->nmh + j] - hks->a3ps[irv_max*hks->nmh + imh_max]);
		}

		hks->sfunc_sigma /= (hks->nzr - 1);
		hks->sfunc_sigma = sqrt(hks->sfunc_sigma);
	
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
			

/*-------------- Select the maxs of sfunction ------------------
 --------------- Rule out some Mohos have overlaps -------------*/
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

					//if ((mh_ub1 < mh_lb2 || mh_ub2 < mh_lb1) && (rv_ub1 < rv_lb2 || rv_ub2 < rv_lb1)) {
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
	printf("Total %d/%d Max sfunctions are kept: statis_option = %d\n", c_maxsfunc, hks->n_maxsfunc, sopt);
	
	return 1;
}		

