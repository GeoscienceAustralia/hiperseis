#ifndef RECORDRF_H
#define RECORDRF_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cstdio>
#include <cmath>
#include <time.h>
#include "sac.h"

#ifdef LINUX
#include <values.h>
#endif


#ifndef	EARTH_RADIUS
#define	EARTH_RADIUS	6371.227
#endif


/*---------------------- Structure Max sfunctions ------------------------*/
typedef struct { 
	float sfunc;
	int kzr;
	int krv; 
	float frv;
	float sigma_rv;
	int kmh;
	float fmh;  
	float sigma_mh;
	float sigma_mhrv;
	float thetd, a, b;
	int isOK;	
} 
MAX_SFUNCTIONS;


/*---------------------- Structure Bootstrap ------------------------*/
typedef struct { 
	int flag;
	int n;   /* number of bootstraps */
	float *frv;
	float *fmh;  
	float rv_bm;
	float mh_bm;
	float sigma_rv;
	float sigma_mh;
	float rho_rvmh;	  
} 
BOOTSTRAP;


/*---------------------- Structure BAZ ------------------------*/
typedef struct { 
	float window;      /* Used for moving BAZs, window of BAZs, like 10.0 deg */
	float step;        /* Used for moving BAZs, step of moving BAZ window, like 5.0 deg (overlap is OK) */  
	int n;             /* Used for moving BAZs, number of BAZ windows */
	int adap_flag;     /* Used for moving BAZs, for the case number of RFs in a window is too many.
								 if set, use adap_n number of RFs by decreasing window size */
	int adap_n;        /* Used for moving BAZs, if ada_flag is set, adap_n is max number of RFs in a window */
	float lb;          /* For both fixed and moving BAZs, low bound of BAZ. It will change for the later case */
	float ub;          /* For both fixed and moving BAZs, up bound of BAZ. it will change for the later case */
} 
BAZ_INFO;


/*---------------------- Structure RF parameter ------------------------*/
typedef struct { 
	char name[64];
	float w_para;   /* Ray parameter */
	float t_sht;
	float baz;
} 
RF_INFO;


/*---------------------- Structure information of HK search ------------------------*/
typedef struct {
	float vp;								       /* Crustal Vp */	 
	float vplb, vpub, vpdt;                 /* Range of crustal Vp and increment */
	float rvlb, rvub, rvdt;                 /* Range of Vp/Vs and increment */   
	float mhlb, mhub, mhdt;                 /* Range of Moho depth and increment */   
	int nvp;                                /* Number of Vp  */    
	int nrv;                                /* Number of Vp/Vs  */   
	int nmh;                                /* Number of Moho depth */    
	int nzr;                                /* = nrv * nmh */   
	float w_0p1s, w_2p1s, w_1p2s;           /* weight factor of 3 phases */
	float *a_0p1s, *a_2p1s, *a_1p2s;        /* Amplitude of 3 phases, with xc weighting */
	float *amp_0p1s, *amp_2p1s, *amp_1p2s;  /* Amplitude of 3 phases, backing up stacked values */
	float *a3ps;                            /* Amplitude of sum of 3 phases with their weight */   
	int *tag;                               /* Index for certain selection */  
	int ntags;                              /* Number of selections */
	int nrfrec;                             /* Number of RFs */
	float delta;
	float twl;                              /* Smoothing time window in sec */
	int ntwl;                               /* Smoothing time window in points */  
	int nroot;                              /* N-th root */
	double nrootpw;                         /* 1 / (N-th root) */
	int xc_flag;                            /* Cross correlation flag */
	int n_maxsfunc;                         /* Number of peaks in S-function */
	float sfunc_max;                        /* The max amplitude in S-function */
	float sfunc_mean;                       /* The mean amplitude of S-function */
	float sfunc_sigma;                      /* The variance of S-function */ 
	float mhrv_cc;                          /* Correlation coefficient between H and k */
	float sfcrt_mh;                         /* Search ranges in sfcrt_mh and sfcrt_rv */ 
	float sfcrt_rv;
	float sfcrt_peak;                       /* for amplitude < sfcrt_peak */  
	MAX_SFUNCTIONS *msf;                    /* struct for peaks in S-function */ 
	/*float *xc;      */                  
	/*float xc_max;   */
	/*float	rv_xcmax; */
	int tsht_flag;                          /* Set time shift flag. set == 1, unset == 0 */
	int wsf_flag;                           /* write sfunction flag. ASCII == 1, Binary == 0 */ 
	int verbose;							       /* Output sfunctions in wgtrv_xc() if set */
} 
HKS_INFO;


/*------------------------------------------------------------*/
/*---------------------- Class RECORD ------------------------*/
/*------------------------------------------------------------*/
class RECORD {
	public:
		RECORD();
		~RECORD();
		
		void Init(char* recName, const float& rayp, const float& tsht);
		void LoadData();
		void FreeData();
		void PeakAmp(HKS_INFO *hks);
		void Ttps_range(HKS_INFO *hks);
		/*void TpTs(HKS_INFO *hks); */
		void Ttps(HKS_INFO *hks, const int& izr);
		void Ttps_rvmh(const float& vp, const float& mh, const float& rv);
		char* GetName() const {return _recName;}
		float Delta() const {return _delta;}
		float Bt() const {return _Bt;}
		float Ray() const {return _ray;}      
		float Rfstart() const {return _rfstart;}  
		float Rfend() const {return _rfend;}
		float t_0p1s;
		float t_2p1s;
		float t_1p2s;
		float a_0p1s;
		float a_2p1s;
		float a_1p2s;
		int n_resmpl;    /* number of resamples */
	
	private:	
		char *_recName;      /* RF sac file name with path */
		float *_data;        /* array to store sac data */
		int _dsize;          /* size of array _data, different from _npts and _nptsr */ 
		float _delta;        /* delta of sac file */
		int _npts;	         /* npts of sac file */ 
		float _Bt;           /* starting time of sac file */
		float _ray;          /* ray parameter in sec/deg */
		float _tsht;         /* time shift in sec */  
		float _rfstart;      /* start time of read */  
		float _rfend;        /* end time of read */
		int _spt;            /* start point of read, calculated from _Bt and _rfstart */
		int _nptsr;	         /* npts of actually read calcultaed from _rfstart and _rfend*/
		float _tp;          /* Travel time of P for H */
		float _ts;          /* Travel time of S for H and K */
		float WaveScan(const float& t, const float& twl); 	
};


#endif
