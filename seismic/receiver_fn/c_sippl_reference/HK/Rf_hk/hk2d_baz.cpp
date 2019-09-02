#include "recordrf.h"
#include "sac.h"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <unistd.h>


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

#ifndef DEFAULT_PEAK
#define DEFAULT_PEAK 0.2  //0.009
#endif
        
#ifndef DEFAULT_SFMH
#define DEFAULT_SFMH 3.0  //2.0
#endif

#ifndef DEFAULT_SFRV
#define DEFAULT_SFRV 0.15 //0.15
#endif


void help();
static int quicksort_baz(RF_INFO *rfinfo, int p, int r);
static int partition_baz(RF_INFO *rfinfo, int p, int r);
int create_hkfname(char *fpth_file, char *fodir, char *hkfile, const int& xc_flag, const int& opt);
int create_hkfnameb(char *fpth_file, char *fodir, char *hkfile, const float& baz_lb, const int& xc_flag, const int& opt);
int wsfunc(char* fpth, HKS_INFO *hks);
int wmh_rvb(char* fpth, HKS_INFO *hks, RECORD *rfrec, BAZ_INFO *bazinfo);
int readline(FILE* fp, char* oneline);
void _errmsg_(char *msg);
int _ss2token_(char *str, char **av, char token);
int set_n(float minVal, float &maxVal, float step);
int mkpath(char *str, char token);
double _acc_(int n, float* data1, float* data2);
int wgt_xc(HKS_INFO *hks);
int wgtrv_xc(HKS_INFO *hks);
int wgtmh_xc(HKS_INFO *hks);
int statis(HKS_INFO *hks);
void nth_stk(RECORD *rfrec, HKS_INFO *hks);


/******************************************************************************/
int main(int argc, char*argv[])
{
	FILE* contrpt;
	FILE* wfpt2;
	float w_para;
	float t_sht;
	float baz;
	float tmpx;
	char *flist = NULL;
	char rf_dir[128];
	char rf_file[128];
	char rf_name[64];
	char fodir1[128];
	char *fodir = NULL;
	char *fsfunc = NULL;
	char *fhk = NULL;
	char fpth_sfunc[128];
	char fpth_hk[128];
	char oneline[256];
	char *av[64];  
	int nrfrec_all; 
	int i;
	int irec;
	int rl;
	int c;
	int ac;
	int sline, tline;
	RECORD *rfrec = NULL;
	HKS_INFO *hks = NULL;
	BAZ_INFO *bazinfo = NULL;
	RF_INFO *rfinfo = NULL;
	

/*----------------------------------------------------------------------------*/
/* First set memory for hks and bazinfo                                       */
/*----------------------------------------------------------------------------*/
	hks = new HKS_INFO;
	if (hks == NULL) {
		fprintf(stderr, "Error in new HKS_INFO hks!\n\n");
		return -1;
	}	
	bazinfo = new BAZ_INFO;
	if (bazinfo == NULL) {
		fprintf(stderr, "Error in new BAZ_INFO bazinfo!\n\n");
		return -1;
	}	
	
/*----------------------------------------------------------------------------*/
/* Default values                                                             */
/*----------------------------------------------------------------------------*/
	hks->vp = 6.2;                         // -P: Vp
	hks->mhlb = 0.0;                       // -Z: mhlb/mhub/mhdt 
	hks->mhub = 0.0;
	hks->mhdt = 0.0;	
	hks->rvlb = 0.0;                       // -K: rvlb/rvub/rvdt  
	hks->rvub = 0.0;
	hks->rvdt = 0.0;
	hks->w_0p1s = DEFAULT_W0P1S;           // -W: 0p1s/2p1s/1p2s
	hks->w_2p1s = DEFAULT_W2P1S; 
	hks->w_1p2s = DEFAULT_W1P2S;
	hks->twl = DEFAULT_TWL;                // -T: Time window for smoothing RF data     
	hks->nroot = DEFAULT_NSTACK;           // -N: n-th root
	hks->xc_flag = 0;                      // -X: 0 - 9  
	bazinfo->lb = 0.0;                     // automatic set in this version
	bazinfo->ub = 360.0;
	bazinfo->adap_flag = 0;                // -R: default not use adaptive number of RFs 
	hks->sfcrt_peak = DEFAULT_PEAK;        // -C: peak/mh/rv 
	hks->sfcrt_mh = DEFAULT_SFMH;          // peak is a percent of (global_max - amp)/global_max 
	hks->sfcrt_rv = DEFAULT_SFRV;          // mh and rv define search range   
	hks->wsf_flag = 0;                     // -o: ASCII sfunction  Default in binary  
	hks->tsht_flag = 1;                    // -s: Unset time shift  Default set time shift
	hks->verbose = 0;
	sline = -INT_MAX;                       // -H: sline/tline
	tline = -INT_MAX;
	
	hks->msf = NULL;
	hks->a_0p1s = NULL;
	hks->a_2p1s = NULL;
	hks->a_1p2s = NULL;
	hks->a3ps = NULL;
//	hks->xc = NULL;

/*----------------------------------------------------------------------------*/
/* Input keyboard options                                                     */
/*----------------------------------------------------------------------------*/
	while ((c=getopt(argc, argv, "P:Z:K:W:T:N:X:R:C:L:H:D:O:osvh")) != -1) {
		switch (c) {
			case 'P':
		    	hks->vp = (float)atof(optarg);
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
				break;
			case 'R':
    			ac = _ss2token_(optarg, av, '/');
    			if (ac == 2) {
					bazinfo->adap_flag = 0;
	    			bazinfo->window = (float)atof(av[0]);
					bazinfo->step = (float)atof(av[1]);
				}
				else if (ac == 3) {	
					bazinfo->adap_flag = 1;
   	 			bazinfo->window = (float)atof(av[0]);
					bazinfo->step = (float)atof(av[1]);
					bazinfo->adap_n = (int)atof(av[2]);
				}
				else {	
					_errmsg_(optarg);
				}	
				break;
			case 'C':
    			ac = _ss2token_(optarg, av, '/');
    			if (ac != 3) _errmsg_(optarg);
	    		hks->sfcrt_peak = (float)atof(av[0]);
	    		hks->sfcrt_mh = (float)atof(av[1]);
	    		hks->sfcrt_rv = (float)atof(av[2]);
				break;
			case 'L':
	    		flist = optarg;
	    		break;
			case 'H':
    			ac = _ss2token_(optarg, av, '/');
    			if (ac != 2) _errmsg_(optarg);
	    		sline = atoi(av[0]);
	    		tline = atoi(av[1]);
	    		break;
			case 'D':
	    		fodir = optarg;
				if (mkpath(fodir, '/')) _errmsg_(fodir);
	    		break;
			case 'O':
	    		ac = _ss2token_(optarg, av, '/');
	    		if (ac != 2) _errmsg_(optarg);
				fsfunc = av[0];
				fhk = av[1];
	    		break;
			case 's':
				hks->tsht_flag = 0;
				break;
			case 'o':
				hks->wsf_flag = 1;
				break;
			case 'v':
				hks->verbose = 1;
				break;
			case 'h':
				help();
				return 0;
			case '?':
				if (optopt == 'W') {
					hks->w_0p1s = DEFAULT_W0P1S;
					hks->w_2p1s = DEFAULT_W2P1S; 
					hks->w_1p2s = DEFAULT_W1P2S;
				}
				else if (optopt == 'T') {
					hks->twl = DEFAULT_TWL;
				}	
				else if (optopt == 'N') {
					hks->nroot = DEFAULT_NSTACK;
				}
				else if (optopt == 'X') {
					hks->xc_flag = 0;
				}	
				else if (optopt == 'R') {
					bazinfo->lb = 0.0;
					bazinfo->ub = 360.0;
				}	
				else if (optopt == 'C') {
					hks->sfcrt_peak = DEFAULT_PEAK;
					hks->sfcrt_mh = DEFAULT_SFMH;
					hks->sfcrt_rv = DEFAULT_SFRV;
				}	
				else if (optopt == 'D') {
					fodir = new char[8];
					strcpy(fodir, "./");
				}
				else if (optopt == 'P' || optopt == 'Z' || optopt == 'K' || 
							optopt == 'L' || optopt == 'H' || optopt == 'O' ||
							optopt == 'R') { 	
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
				return 0;
		}
	}

   
/*----------------------------------------------------------------------------*/
/* Print, check and reset keyboard options                                    */
/*----------------------------------------------------------------------------*/
	printf("\n");	 

/*----------------------------------------------------------------------------*/
/* Print -P: crustal Vp                                                       */
/*----------------------------------------------------------------------------*/
	printf("Vp_crust = %f\n",	hks->vp);

/*----------------------------------------------------------------------------*/
/* Print, check and reset -Z and -K options
	Set the number of H and k search 
	The upbound of H and k may change according to their lowbound and increment*/	
/*----------------------------------------------------------------------------*/
	if (hks->mhlb == 0.0 || hks->mhdt == 0.0 || hks->mhub == 0.0) {
		fprintf(stderr, "H_search = %f : %f : %f\n",	hks->mhlb, hks->mhdt, hks->mhub);
		fprintf(stderr, "Incorrect search range for H.\n\n");
		return -1;
	}	
	else {
		hks->nmh = set_n(hks->mhlb, hks->mhub, hks->mhdt);
		printf("H_search = %f : %f : %f(may reset), n = %d\n", 
				hks->mhlb, hks->mhdt, hks->mhub, hks->nmh);		
	}
		
	if (hks->rvlb == 0.0 || hks->rvdt == 0.0 || hks->rvub == 0.0) {
		fprintf(stderr, "K_search = %f : %f : %f\n",	hks->rvlb, hks->rvdt, hks->rvub);
		fprintf(stderr, "Incorrect search range for K.\n\n");
		return -1;	
	}
	else {
		hks->nrv = set_n(hks->rvlb, hks->rvub, hks->rvdt);
		printf("K_search = %f : %f : %f(may reset), n = %d\n", 
				hks->rvlb, hks->rvdt, hks->rvub, hks->nrv);
	}	
		
/*----------------------------------------------------------------------------*/
/* Print, check and reset -W option                                           
	Reset the weight factor of 3 phases if their summation is not equal to 1   */
/*----------------------------------------------------------------------------*/
	tmpx = hks->w_0p1s + hks->w_2p1s + hks->w_1p2s;
	if (tmpx == 1.0) {
		printf("Weighting = %f, %f, %f\n",	hks->w_0p1s, hks->w_2p1s, hks->w_1p2s);
	}
	else {	
		hks->w_0p1s /= tmpx;
		hks->w_2p1s /= tmpx;
		hks->w_1p2s /= tmpx;
		printf("Weighting = %f, %f, %f(reset)\n",	hks->w_0p1s, hks->w_2p1s, hks->w_1p2s);
	}
	
/*----------------------------------------------------------------------------*/
/* Print -T: Time window for smoothing RF data                                */
/*----------------------------------------------------------------------------*/	
	printf("Smoothing time window = %f\n", hks->twl);
	
/*----------------------------------------------------------------------------*/
/* Print -N: n-th root                                                        */
/*----------------------------------------------------------------------------*/		
	printf("N-th root = %d\n", hks->nroot);
	hks->nrootpw = 1.0/(double)hks->nroot;
	
/*----------------------------------------------------------------------------*/
/* Print and check -X: cross-correlation option                               */
/*----------------------------------------------------------------------------*/		
	if (hks->xc_flag < 0 || hks->xc_flag >= 10) {
		fprintf(stderr, "Possible argument (= %d?) for option -X is [0 - 9].\n\n", 
				hks->xc_flag);
		help();
		return -1;
	} 
	else {				 				
		printf("Cross-correlation weighting flag = %d\n", hks->xc_flag);
	}
	
/*----------------------------------------------------------------------------*/
/* Print -R: moving BAZ window info                                           */
/*----------------------------------------------------------------------------*/			
	if (bazinfo->adap_flag) 
		printf("Adaptive BAZ window: Window = %f, Step = %f, adap_N = %d\n", 
				bazinfo->window, bazinfo->step, bazinfo->adap_n);	
	else
		printf("Not use adaptive BAZ window: Window = %f, Step = %f\n", 
				bazinfo->window, bazinfo->step);	
	if (bazinfo->window == 0.0 || bazinfo->step == 0.0) {
		fprintf(stderr, "Error in BAZ moving window and/or step.\n\n");
		return -1;
	}		
	
/*----------------------------------------------------------------------------*/
/* Print and check -L: RF data list file                                      */
/*----------------------------------------------------------------------------*/		
	if (flist == NULL) {
		fprintf(stderr, "RF data list file = NULL for option -L.\n\n");
		return -1;
	}
	else {	
		printf("RF data list file = '%s'\n", flist);
	}
		
/*----------------------------------------------------------------------------*/
/* Print and check -H: sline/tline                                            */
/*----------------------------------------------------------------------------*/		
	if (tline < 1) {
		fprintf(stderr, "Possible total lines (tline = %d?) for option -H must be >= 1.\n\n", 
				tline);
		help();
		return -1;
	}				 
	if (sline == 0) {
		fprintf(stderr, "sline = %d. Not yet implement to use RFs with different Gaussian factors.\n\n",
				sline);
		return -1;
	}
	printf("sline/tline = %d/%d, read line %d from total %d lines\n", 
		sline, tline, sline, tline);
		
/*----------------------------------------------------------------------------*/
/* Print and check -D: output dir
	If no -D, output in the current dir                                        */
/*----------------------------------------------------------------------------*/		
	if (fodir == NULL) {
	 	fodir = new char[8];
		strcpy(fodir, "./");
	}
	if (fodir[strlen(fodir)-1] == '/') sprintf(fodir1, "%s", fodir);	
	else sprintf(fodir1, "%s/", fodir);	
	printf("Output dir = '%s'\n", fodir1);
	
/*----------------------------------------------------------------------------*/
/* Print and check -O and -o: output file name and format                     */
/*----------------------------------------------------------------------------*/		
	if (fsfunc == NULL) {
		fprintf(stderr, "Output file: sfunc = NULL for option -O.\n\n");
		return -1;
	}
	if (fhk == NULL) {
		fprintf(stderr, "Output file: mh_rv = NULL for option -O.\n\n");
		return -1;
	}
	if (hks->wsf_flag) {
		printf("Output file: sfunc = '%s'(in ASCII); hk = '%s'\n", fsfunc, fhk);
	}
	else {
		printf("Output file: sfunc = '%s'(in Binary); hk = '%s'\n", fsfunc, fhk);
	}
	
/*----------------------------------------------------------------------------*/
/* Print -s: Unset time shift                                                 */
/*----------------------------------------------------------------------------*/
	if (hks->tsht_flag) printf("Set time shift\n");
	else printf("Unset time shift\n");	

/*----------------------------------------------------------------------------*/
/* Print -v: set verbose                                                      */
/*----------------------------------------------------------------------------*/
	if (hks->verbose) {
		printf("Output sfunctions in middle way to illustrate. (Only works for x1)\n");
	}
	

/*----------------------------------------------------------------------------*/
/* The total number of grids in H and K                                       */
/*----------------------------------------------------------------------------*/
	hks->nzr = hks->nrv * hks->nmh;

/*----------------------------------------------------------------------------*/
/* The amplitude array turns to 1D, the total number is nzr                   */
/*----------------------------------------------------------------------------*/
	hks->a3ps = new float[hks->nzr];
	if (hks->a3ps == NULL) {
		fprintf(stderr, "Error in new float hks->a3ps!\n\n");
		delete hks;
		delete bazinfo;
		return -1;
	}	
	hks->a_0p1s = new float[hks->nzr];
	if (hks->a_0p1s == NULL) {
		fprintf(stderr, "Error in new float hks->a_0p1s!\n\n");
		delete hks;
		delete bazinfo;
		return -1;
	}	
	hks->a_2p1s = new float[hks->nzr];
	if (hks->a_2p1s == NULL) {
		fprintf(stderr, "Error in new float hks->a_2p1s!\n\n");
		delete hks;
		delete bazinfo;
		return -1;
	}	
	hks->a_1p2s = new float[hks->nzr];
	if (hks->a_1p2s == NULL) {
		fprintf(stderr, "Error in new float hks->a_1p2s!\n\n");
		delete hks;
		delete bazinfo;
		return -1;
	}	
	hks->tag = new int[hks->nzr];
	if (hks->tag == NULL) {
		fprintf(stderr, "Error in new int hks->tag!\n\n");
		delete hks;
		delete bazinfo;
		return -1;
	}	
	for (i = 0; i < hks->nzr; i++) {
		hks->a3ps[i] = 0.0;
		hks->a_0p1s[i] = 0.0;
		hks->a_2p1s[i] = 0.0;
		hks->a_1p2s[i] = 0.0;
		hks->tag[i] = 0;
	}	

/*----------------------------------------------------------------------------*/
/* Read Dir of sac files according to sline
	First scan RF lines and calculate the number of RF records                 */
/*----------------------------------------------------------------------------*/
 	contrpt = fopen(flist, "r");
	if (contrpt == NULL) {
		fprintf(stderr, "\nError, can not open control file %s!!\n\n",flist);
   	return -2;
 	}

/*----------------------------------------------------------------------------*/
/* Locate line "sline" from lines "tline"                                     */
/*----------------------------------------------------------------------------*/
	i = 1;
	while(i <= tline)  {
		rl = readline(contrpt, oneline);
		if (rl == 0) { 
			fprintf(stderr, "Imcomplete file: %s.\n", flist);
			return -1;
		}
		else if (rl == 1) {
			if (i == sline) {
				sscanf(oneline, "%s", rf_dir);
			}
			i++;	
		}
	}		
	printf("Dir of RF files = '%s'\n", rf_dir);
		
/*----------------------------------------------------------------------------*/
/* Calculate all of the number of RF records in the list                      */
/*----------------------------------------------------------------------------*/
	nrfrec_all = 0;
	while(1) {
		rl = readline(contrpt, oneline);
		if (rl == 0) break;
		else if (rl == 1) nrfrec_all++;  						
	}
	
/*----------------------------------------------------------------------------*/
/* Set memory for all RFs                                                     */
/*----------------------------------------------------------------------------*/
	printf("Locate total %d receiver functions for processing\n", nrfrec_all);
	if (nrfrec_all == 0) {
		printf("No RFs are read.\n");
		delete [] hks->a3ps;
		hks->a3ps = NULL;
		delete [] hks->a_0p1s;
		hks->a_0p1s = NULL;
		delete [] hks->a_2p1s;
		hks->a_2p1s = NULL;
		delete [] hks->a_1p2s;
		hks->a_1p2s = NULL;
		delete [] hks->tag;
		hks->tag = NULL;
		delete hks;
		hks = NULL;
		delete bazinfo;
		bazinfo = NULL;
		return 0;
	}
	
	rfrec = new RECORD[nrfrec_all];
	if (rfrec == NULL) {
		fprintf(stderr, "Error in new RECORD rfrec!\n\n");
		delete [] hks->a3ps;
		hks->a3ps = NULL;
		delete [] hks->a_0p1s;
		hks->a_0p1s = NULL;
		delete [] hks->a_2p1s;
		hks->a_2p1s = NULL;
		delete [] hks->a_1p2s;
		hks->a_1p2s = NULL;
		delete [] hks->tag;
		hks->tag = NULL;
		delete hks;
		hks = NULL;
		delete bazinfo;
		bazinfo = NULL;
		return -1;
	}	

	rfinfo = new RF_INFO[nrfrec_all];
	if (rfinfo == NULL) {
		fprintf(stderr, "Error in new RF_INFO rfinfo!\n\n");
		delete [] hks->a3ps;
		hks->a3ps = NULL;
		delete [] hks->a_0p1s;
		hks->a_0p1s = NULL;
		delete [] hks->a_2p1s;
		hks->a_2p1s = NULL;
		delete [] hks->a_1p2s;
		hks->a_1p2s = NULL;
		delete [] hks->tag;
		hks->tag = NULL;
		delete hks;
		hks = NULL;
		delete bazinfo;
		bazinfo = NULL;
		return -1;
	}	

/*----------------------------------------------------------------------------*/
/* Prepare output file name like hk2d.HL.AXX.x7                               */	
/*----------------------------------------------------------------------------*/
	if (create_hkfname(fpth_hk, fodir1, fhk, hks->xc_flag, 2)) return -3;

/*----------------------------------------------------------------------------*/
/* Write summary of H-k search to hk2d.HL.AXX.x7 header                       */
/*----------------------------------------------------------------------------*/
	wfpt2 = fopen(fpth_hk, "w");
	if (wfpt2 == NULL) {
 		fprintf(stderr, "\nError: Cannot open output file %s!!\n\n", fpth_hk);
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
	if (bazinfo->adap_flag) 
		fprintf(wfpt2, "Adaptive BAZ window: Window = %f, Step = %f, adap_N = %d\n", 
				bazinfo->window, bazinfo->step, bazinfo->adap_n);	
	else
		fprintf(wfpt2, "Not use adaptive BAZ window: Window = %f, Step = %f\n", 
				bazinfo->window, bazinfo->step);	
	fprintf(wfpt2, "Process %d RF traces\n", hks->nrfrec);			
	fprintf(wfpt2, "Search H = %3.1f K = %4.2f sfunction peaks < %4.2f\n", 
					hks->sfcrt_mh, hks->sfcrt_rv, hks->sfcrt_peak);
	if (hks->tsht_flag) fprintf(wfpt2, "Set time shift\n");
	else fprintf(wfpt2, "Unset time shift\n");
	if (hks->wsf_flag) fprintf(wfpt2, "Write sfunction in ASCII\n\n\n");
	else fprintf(wfpt2, "Write sfunction in BINARY\n\n\n");
	fprintf(wfpt2, "%%hk_file               BAZLB  BAZUB No.RFs XC  Moho                Vp/Vs\n\n");
	fclose(wfpt2);
	

/*----------------------------------------------------------------------------*/
/* Start BAZ loop: process RFs in a window in one time and increase by step
	The low bound of 1st window bazinfo->lb = 0.0                              */	
/*----------------------------------------------------------------------------*/
	bazinfo->lb = 0.0;	
	
/*----------------------------------------------------------------------------*/
/* BAZ loop continues if bazinfo->lb < 360.0
	bazinfo->lb will increase by bazinfo->step                                 */	
/*----------------------------------------------------------------------------*/
	while (bazinfo->lb < 360.0) {

/*----------------------------------------------------------------------------*/
/* The up bound of BAZ is defined by bazinfo->lb and bazinfo->window          */	
/*----------------------------------------------------------------------------*/
		bazinfo->ub = bazinfo->lb + bazinfo->window;

/*----------------------------------------------------------------------------*/
/* Search RF list and choose RFs                                              */	
/*----------------------------------------------------------------------------*/
	 	rewind(contrpt);
		i = 1;
		while (i <= tline) {
			rl = readline(contrpt, oneline);
			if (rl == 1) i++;
		}
		
		irec = 0;	
		while(1) {
			rl = readline(contrpt, oneline);
			if (rl == 0) break;
			else if (rl == 1) {
				sscanf(oneline, "%s%f%f%f%*[^\n]", rf_name, &w_para, &t_sht, &baz);

/*----------------------------------------------------------------------------*/
/* This case is for bazinfo->ub not beyond 360 
	like using [lb = 30.0 - ub = 120.0]                                        */
/*----------------------------------------------------------------------------*/
				if (bazinfo->ub <= 360.0) { 
					if (baz >= bazinfo->lb && baz < bazinfo->ub) {					
						strcpy(rfinfo[irec].name, rf_name);
						rfinfo[irec].w_para = w_para;
						rfinfo[irec].t_sht = t_sht;
						rfinfo[irec].baz = baz;
						irec++;
					}
				}	
/*----------------------------------------------------------------------------*/
/* This case is for bazinfo->ub beyond 360 
	like using [lb = 350.0 - ub = 10.0]                                        */
/*----------------------------------------------------------------------------*/
				else {
					if (baz >= bazinfo->lb || baz < bazinfo->ub-360.0) {
						
						strcpy(rfinfo[irec].name, rf_name);
						rfinfo[irec].w_para = w_para;
						rfinfo[irec].t_sht = t_sht;
						rfinfo[irec].baz = baz;
						irec++;
					}					
				}
			}
		}		    
		hks->nrfrec = irec;
		
/*----------------------------------------------------------------------------*/
/* If the selected number of RFs in this window is less than 5 
	Neglect them and go to next window by increasing bazinfo->step or
	bazinfo->window                                                            */
/*----------------------------------------------------------------------------*/
		if (hks->nrfrec < 5) {
			printf("Neglect baz %f - %f\n", bazinfo->lb, bazinfo->ub);
			if (bazinfo->ub - bazinfo->lb < bazinfo->step) bazinfo->lb = bazinfo->ub;
			else bazinfo->lb += bazinfo->step;
			continue;
		}

/*----------------------------------------------------------------------------*/
/* This part is for bazinfo->adap_flag is set, and selected number of RFs in 
	this window > the max RFs determined by bazinfo->adap_n.                   */
/*----------------------------------------------------------------------------*/
		if (bazinfo->adap_flag && hks->nrfrec > bazinfo->adap_n) {
		
/*----------------------------------------------------------------------------*/
/* If bazinfo->ub is beyond 360.0, adjust baz for sorting.                    */
/*----------------------------------------------------------------------------*/
			if (bazinfo->ub > 360.0) { 
				for (irec = 0; irec < hks->nrfrec; irec++) {
					if (rfinfo[irec].baz < bazinfo->ub-360.0)
						rfinfo[irec].baz += 360.0;
				}						
			}
					
/*----------------------------------------------------------------------------*/
/* Sort struct rfinfo in terms of baz in ascend order.                        */
/*----------------------------------------------------------------------------*/
			quicksort_baz(rfinfo, 0, hks->nrfrec-1);
			
/*----------------------------------------------------------------------------*/
/* Reset hks->nrfrec                                                          */
/*----------------------------------------------------------------------------*/
			hks->nrfrec = bazinfo->adap_n;
			
/*----------------------------------------------------------------------------*/
/* Reset bazinfo->ub
	Because bazinfo->lb and ub are integer numbers, it is possible that
	bazinfo->lb == bazinfo->ub only use adap_n RFs. So that automatically 
	increase used RFs number until bazinfo->ub > bazinfo->lb                                                           */
/*----------------------------------------------------------------------------*/
			while ((int)bazinfo->lb == (int)rfinfo[hks->nrfrec-1].baz) {
				hks->nrfrec++;
			}	 
			bazinfo->ub = 1.0 * (int)rfinfo[hks->nrfrec-1].baz;
								
		}
		
/*----------------------------------------------------------------------------*/
/* Set up selected RF records. The real RF sac data is not read at this stage 
	in order to handle large number of RF records                              */	
/*----------------------------------------------------------------------------*/
		for (irec = 0; irec < hks->nrfrec; irec++) {

			strcpy(rf_file, rf_dir);	
			if (rf_dir[strlen(rf_dir)-1] != '/') strcat(rf_file, "/");		
			strcat(rf_file, rfinfo[irec].name);

/*----------------------------------------------------------------------------*/
/* Unset time shift if tsht_flag == 0
	Assign RF file name, ray parameter, and shift time.                        */
/*----------------------------------------------------------------------------*/
			if (!hks->tsht_flag) t_sht = 0.0;	
			rfrec[irec].Init(rf_file, w_para, t_sht);
					
/*----------------------------------------------------------------------------*/
/* Calcuate and assign the start time and end time of reading RF data
	Not load RF data at this stage, instead loading it when stacking.          */
/*----------------------------------------------------------------------------*/
			rfrec[irec].Ttps_range(hks);					
						
		}    

		printf("\nReading %d receiver functions\n", hks->nrfrec);
		printf("Processing baz %f - %f\n", bazinfo->lb, bazinfo->ub);

/*----------------------------------------------------------------------------*/
/* Calculate travel times of 0p1s, 2p1s, and 1p2s starting from P, 
	pick up amplitude for each RF trace 	
	stack all three phase for all RF traces
	Actually read RF data at this stage                                        */	
/*----------------------------------------------------------------------------*/
		nth_stk(rfrec, hks);

/*----------------------------------------------------------------------------*/
/* Cross-corelation between 0p1s, 2p1s, 1p2s                                  */
/*----------------------------------------------------------------------------*/
		if (hks->xc_flag == 1 || hks->xc_flag == 2 || hks->xc_flag == 3) {
			if (wgtrv_xc(hks)) return -1;
		}
		else if (hks->xc_flag == 4 || hks->xc_flag == 5 || hks->xc_flag == 6) {
			if (wgtmh_xc(hks)) return -1;
		}
		else if (hks->xc_flag == 7 || hks->xc_flag == 8 || hks->xc_flag == 9) {
			if (wgt_xc(hks)) return -1;
		}	

/*----------------------------------------------------------------------------*/
/* Search for the maximum, do statistics                                      */	
/*----------------------------------------------------------------------------*/
		if (statis(hks)) return -1;

/*----------------------------------------------------------------------------*/
/* Output S-function at each grid for each BAZ                                */
/*----------------------------------------------------------------------------*/
		if (create_hkfnameb(fpth_sfunc, fodir1, fsfunc, bazinfo->lb, hks->xc_flag, 2)) return -3;

		if (wsfunc(fpth_sfunc, hks)) {
			fprintf(stderr, "\nError in writing sfunc in '%s'!!\n\n", fpth_sfunc);
 			return -2;
		}
		printf("Output %s.\n", fpth_sfunc);

/*----------------------------------------------------------------------------*/
/* Output H-k search results into hk2d.HL.AXX.x7 
	and hk2d.HL.AXX.x7.b010.00, hk2d.HL.AXX.x7.b010.01 ...                     */	
/*----------------------------------------------------------------------------*/
		if (wmh_rvb(fpth_hk, hks, rfrec, bazinfo)) {
			fprintf(stderr, "\nError in writing mh_rv in '%s'!!\n\n", fpth_hk);
			return -2;
		}	
		printf("Output %s.\n", fpth_hk);
	
		if (bazinfo->ub - bazinfo->lb < bazinfo->step) bazinfo->lb = bazinfo->ub;
		else bazinfo->lb += bazinfo->step;

/*----------------------------------------------------------------------------*/
/* BAZ loop end                                                               */
/*----------------------------------------------------------------------------*/
	}	
	fclose(contrpt);
	
	delete [] hks->msf;
	hks->msf = NULL;
	delete hks;
	hks = NULL;
	delete [] rfrec;
	rfrec = NULL;
	delete bazinfo;
	bazinfo = NULL;
	delete rfinfo;
	rfinfo = NULL;
	
	return 0;				
}


/******************************************************************************/
/* Help                                                                       */	
/******************************************************************************/
void help()
{
	printf("Usage:\n");
	printf("hd2d_baz -P 6.5254 -Z 10.0/100.0/0.1 -K 1.50/2.00/0.01 -W 0.6/0.2/0.2 -T 0.1 -N n \n"); 
	printf("         -R 10.0/5.0/20 -L list -H 1/3 -s\n");
	printf("         -D dirout -O sfr2d.AH.ANQ/hkr2d.AH.ANQ \n");
	printf("Option = -P, Argument = Average crustal velocity\n");
	printf("Option = -Z, Argument = Moho Depth search. start/end/increment\n");
	printf("Option = -K, Argument = Vp/Vs ratio search. start/end/increment\n");
	printf("Option = -W, Argument = Weighing factor for 0p1s/2p1s/1p2s. Default = 0.7/0.2/0.1\n");
	printf("Option = -T, Argument = Smoothing time window in sec. Default = 0.1\n");
	printf("Option = -N, Argument = n-th stacking. Default = 1\n");
	printf("Option = -X, Argument = XC weight. Default = 0\n"); 
	printf("             0: No XC weight; \n"); 
	printf("             1: All phase; 2: 0p1s+2p1s; 3: 0p1s+1p2s. As a function of Vp/Vs.\n");
	printf("             4: All phase; 5: 0p1s+2p1s; 6: 0p1s+1p2s. As a function of Moho depth.\n");
	printf("             7: All phase; 8: 0p1s+2p1s; 9: 0p1s+1p2s. As a function of Vp/Vs and Moho.\n");
	printf("Option = -R, Argument = BAZ window/moving step.\n"); 
	printf("             or Argument = BAZ window/moving step/adaptive number.\n");
	printf("Option = -C, Argument = Dpeak/mh/rv. The values control peak search. \n");
	printf("             Search for peaks p[i] whose (p_max-p[i])/p_max < Dpeak, \n");
	printf("             but the location of p[i] must be outside of an area defined by mh and rv. \n");
	printf("             Usually use default = 0.2/3.0/0.15.\n"); 
	printf("Option = -L, Argument = list: The file name of list of time traces.\n");
	printf("Option = -H, Argument = Indicate the dir of RF list. sline/tline\n");
	printf("             tline = n: There are n dirs (different Gaussian factors). RFs names are the same\n");
	printf("             sline = i (i <= n) Read and use the i-th dir's RFs.\n");
	printf("             sline = 0 May use all dirs' RFs.\n");
	printf("Option = -D, Argument = dirout: The dir for output file.\n");
	printf("Option = -O, Argument = output files: sfunc and mh_rv2d.\n");
	printf("Option = -o, No argument. To use ASCII output for sfunctions.\n");	
	printf("Option = -s, No argument. To to unset Time shift.\n");
	printf("\n");
}



/******************************************************************************/
/* Quick Sorting in ascend order                                              */	
/******************************************************************************/
static int quicksort_baz(RF_INFO *rfinfo, int p, int r)
{
	int q;
		
	if (p < r) {
		q = partition_baz(rfinfo, p, r);
		quicksort_baz(rfinfo, p, q-1);
		quicksort_baz(rfinfo, q+1, r);
//		printf("(");
//		for(int j = 0; j <= q-1; j++) printf("%4d", S[j]);
//		printf(") %4d (", S[q]);
//		for(int j = q+1; j <= r; j++) printf("%4d", S[j]);
//		printf(")\n");
	}
	
	return 0;
}


/******************************************************************************/
/* Partition in Quick Sorting
	subroutine for quicksort_baz                                               */	
/******************************************************************************/
static int partition_baz(RF_INFO *rfinfo, int p, int r)
{
	float baz, tmp;
	char tmpname[64];
	int i, j;
	
	baz = rfinfo[r].baz;
	i = p-1;
	for (j = p; j < r; j++) {
		if (rfinfo[j].baz <= baz) {
			i++;
			tmp = rfinfo[i].baz;
			rfinfo[i].baz = rfinfo[j].baz;
			rfinfo[j].baz = tmp;
			
			tmp = rfinfo[i].w_para;
			rfinfo[i].w_para = rfinfo[j].w_para;
			rfinfo[j].w_para = tmp;
			
			tmp = rfinfo[i].t_sht;
			rfinfo[i].t_sht = rfinfo[j].t_sht;
			rfinfo[j].t_sht = tmp;
			
			strcpy(tmpname, rfinfo[i].name);
			strcpy(rfinfo[i].name, rfinfo[j].name);
			strcpy(rfinfo[j].name, tmpname);
		}
	}
	
	tmp = rfinfo[i+1].baz;
	rfinfo[i+1].baz = rfinfo[r].baz;
	rfinfo[r].baz = tmp;
	
	tmp = rfinfo[i+1].w_para;
	rfinfo[i+1].w_para = rfinfo[r].w_para;
	rfinfo[r].w_para = tmp;
	
	tmp = rfinfo[i+1].t_sht;
	rfinfo[i+1].t_sht = rfinfo[r].t_sht;
	rfinfo[r].t_sht = tmp;
	
	strcpy(tmpname, rfinfo[i+1].name);
	strcpy(rfinfo[i+1].name, rfinfo[r].name);
	strcpy(rfinfo[r].name, tmpname);
	

	return i+1;
}

