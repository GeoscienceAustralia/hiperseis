#include "sac.h"
#include "recordrf.h"


#ifndef EPSON
#define EPSON       0.01
#endif

/******************************************************************************/
/* ----------------------- class Record --------------------------------------*/
/******************************************************************************/
RECORD::RECORD()
{
/*----------------------------------------------------------------------------*/
/* Private                                                                    */
/*----------------------------------------------------------------------------*/
    _recName = NULL;
    _data = NULL;
    _delta = 0.0;        
    _npts = 0;           
    _Bt = 0;           
    _ray = 0.0;          
    _tsht = 0.0;         
    _rfstart = 0.0;      
    _rfend = 0.0;
    _spt = 0;
    _nptsr = 0;                     
    _tp = 0.0;
    _ts = 0.0;

/*----------------------------------------------------------------------------*/
/* Public                                                                     */
/*----------------------------------------------------------------------------*/
    t_0p1s = 0.0;
    t_2p1s = 0.0;
    t_1p2s = 0.0;   
    a_0p1s = 0.0;
    a_2p1s = 0.0;
    a_1p2s = 0.0;   
    n_resmpl = 0;
}


/******************************************************************************/
/******************************************************************************/
RECORD::~RECORD()
{
    if (_recName != NULL) {
        delete[] _recName;
        _recName = NULL;
    }   
    if (_data != NULL) {
        delete[] _data;
        _data = NULL;
    }   
    /*
    if (_tp != NULL) {
        delete[] _tp;
        _tp = NULL;
    }
    if (_ts != NULL) {  
        delete[] _ts;
        _ts = NULL;
    }
    */  
}   
            

/******************************************************************************/
/* Set unchanged and initial values to members in record                      */
/******************************************************************************/
void RECORD::Init(char* recName, const float& rayp, const float& tsht)  
{
    if (_recName != NULL) delete[] _recName;    
    _recName = new char[strlen(recName) + 1];
    strcpy(_recName, recName);
    
    _ray = rayp;
    _tsht = tsht;
}
                                

/******************************************************************************/
/* Load RF sac file into '_data', whose size is determined by rfend and Bt.                       
    When operation for '_data' is finished, free memory using FreeData.
    If the class is claimed as an array, try to make only one sac file exists
    at one time, but some sac header information of all members in the array 
    are always kept.      
     0    1    2    3   4   5   6   7   8   9   10  11  12  13  14  15  16  ...
    -0.4 -0.3 -0.2 -0.1 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 ...
    Bt(adjusted)................rfstart.................rfend..................
    |--------------------------------------------------------------| data size 
    _data: only between rfstart and rfend has sac file values
     0    0    0    0   0   0   d   d   d   d    d   d   d   0   0             */
/******************************************************************************/
void RECORD::LoadData()
{
    SAC_HEADER* sachd;
    float et;
    
    extern int rsac4(char *kname, float* _data, const int& spt, const int& nptsr);
    extern int zero(float*, int);
    
    sachd = new SAC_HEADER;
    sachd->rsacheader(_recName);

/*----------------------------------------------------------------------------*/
/* Assign some sac header info to RECORD values 
    Correct by time shift to free surface P at zero                            
    et (the time of last sample) is calculated from the adjusted _Bt, this way
    makes the zero lay on the P-wave peak                                      */
/*----------------------------------------------------------------------------*/
    _Bt = sachd->_b - _tsht;   
    _delta = sachd->_delta;
    _npts = sachd->_npts;
    et = (_npts - 1) * _delta + _Bt;
    
/*----------------------------------------------------------------------------*/
/* If et < _rfend, read [rf_start, et].
    The calculation is based on adjusted _Bt                                   */   
/*----------------------------------------------------------------------------*/
    if (et > _rfend) {      
        _nptsr = (int)((_rfend - _rfstart) / _delta) + 1;
    }   
    else {
        _nptsr = (int)((et - _rfstart) / _delta) + 1;
    }   
    _spt = (int) ((_rfstart - _Bt) / _delta);

/*----------------------------------------------------------------------------*/
/* Determine the number of points according to this RF sac's delta.
    We give more points for the size of "_data" array.
    They are initially set to zeros.
    Only read sac in [rf_start, rf_end] and assign to "_data"                     */    
/*----------------------------------------------------------------------------*/
    _dsize = _spt + _nptsr + 20;

    if (_data != NULL) {
        delete [] _data;
        _data = NULL;
    }   
    _data = new float[_dsize];
    if (_data == NULL) {
        fprintf(stderr, "Not enough memory for RF data!\n");
        fprintf(stderr, "Error in 'void RECORD::LoadData()'\n\n");
        exit(-1);
    }   
            
    zero(_data, _dsize);

    if (!rsac4(_recName, _data, _spt, _nptsr)) {
        fprintf(stderr, "\nError in LoadData(): Reading sac data for '%s'!!\n\n", _recName);
        exit(-1);
    }
            
    delete sachd;   
    sachd = NULL;
}   


/******************************************************************************/
/******************************************************************************/
void RECORD::FreeData()
{
    if (_data != NULL) {
        delete[] _data;
        _data = NULL;
    }   
}


/******************************************************************************/
/* Calcualte waveform amplitude at each phase arrival time
    t_0p1s, t_2p1s, and t_1p2s MUST have been calculated                       */
/******************************************************************************/
void RECORD::PeakAmp(HKS_INFO *hks)
{           
    a_0p1s = WaveScan(t_0p1s, hks->twl);
    a_2p1s = WaveScan(t_2p1s, hks->twl);
    a_1p2s = WaveScan(t_1p2s, hks->twl);        
}   


/*
void RECORD::PeakAmp(HKS_INFO *hks)
{
    int nzr;
    int izr;
    
    _data = new float[_npts];
        
    if (!rsac0(_recName, _data, _npts, _delta, _Bt)) {
        fprintf(stderr, "\nError in reading sac data for '%s'!!\n\n", _recName);
        exit(-1);
    }

    a_0p1s = new float[hks->nzr];
    a_2p1s = new float[hks->nzr];
    a_1p2s = new float[hks->nzr];
    
    for (izr = 0; izr < hks->nzr; izr++) {
        a_0p1s[izr] = WaveScan(t_0p1s[izr], hks->ntwl);
        a_2p1s[izr] = WaveScan(t_2p1s[izr], hks->ntwl);
        a_1p2s[izr] = WaveScan(t_1p2s[izr], hks->ntwl);     
    }
    
    delete [] _data;
}   
*/

        
/******************************************************************************/
/* Scan waveform till to "t", compute average amplitude within window "twl" 
    centered at "t".
    The search range is [_spt, _spt+_nptsr-1]                                  */
/******************************************************************************/
float RECORD::WaveScan(const float& t, const float& twl)
{
    float amp;
    float   rbuf;
    int ntwl;
    int ipt, i;
    int lft, rgt;
    
/*----------------------------------------------------------------------------*/
/* Set points of smoothing time window                                        */    
/*----------------------------------------------------------------------------*/
    rbuf = twl / _delta;
    if ((rbuf - (int)rbuf) > EPSON * _delta) {
        ntwl = (int)ceil(rbuf) + 1;
    }
    else{
        ntwl = (int)rbuf + 1;
    }

/*----------------------------------------------------------------------------*/
/* We can do this only when P-wave corresponding to point zero                */
/*----------------------------------------------------------------------------*/
    ipt = (int)((t - _Bt) / _delta);   
    lft = (ntwl - 1) / 2;
    rgt =  ntwl / 2;
        
    if (_data == NULL) {
        fprintf(stderr, "Error in using _data in float RECORD::WaveScan()!\n\n");
        exit(-1);
    }   
    
    amp = 0.0;
    for (i = ipt - lft; i <= ipt + rgt; i++) {
//      if (i >= _spt && i < (_spt+_nptsr)) {   
        if (i >= 0 && i < _dsize) { 
            amp += _data[i];
        }
    }   
    amp /= ((ipt+rgt) - (ipt-lft) + 1);

/*  amp = 0.0;
    icount = 0; 
    for (i = ipt - lft; i <= ipt + rgt; i++) {  
        if (i >= 0 && i < _npts) {  
            amp += _data[i];
            icount++;
        }
    }
    
    if (icount == 0) amp = 0.0;
    else amp /= icount; 
*/          
    return amp;     
}   

        
/******************************************************************************/
/* Calculate the first and the last time point that will be search
    t_0p1s = H * sqrt(1/vs/vs - p*p) - H * sqrt(1/vp/vp - p * p)
             = H * sqrt(k*k/vp/vp - p*p) - H * sqrt(1/vp/vp - p * p)
    This case uses the mhlb and rvlb, the uniform vp, and ray of each RF         
    t_1p2s = 2 * H * sqrt(1/vs/vs - p*p)
             = 2 * H * sqrt(k*k/vp/vp - p*p)
    This case uses the mhub and rvub, the uniform vp, and ray of each RF
    Usually, t_0p1s > 2.0 which means (H = 15km, k = 1.5) is impossible.
    So that set _rfstart = 2.0 s if it is less than 2.0 s.                     */
/******************************************************************************/
void RECORD::Ttps_range(HKS_INFO *hks)
{
    float p;
    float slow_p, slow_s;
        
    p = _ray * 180.0 / (EARTH_RADIUS - hks->mhlb) / M_PI;
    slow_p = sqrt(1.0 / hks->vp / hks->vp - p * p);
    slow_s = sqrt(hks->rvlb * hks->rvlb / hks->vp / hks->vp - p * p);
    _rfstart = hks->mhlb * (slow_s - slow_p);
    
    if (_rfstart < 2.0) _rfstart = 2.0;
        
    p = _ray * 180.0 / (EARTH_RADIUS - hks->mhub) / M_PI;
    slow_s = sqrt(hks->rvub * hks->rvub / hks->vp / hks->vp - p * p);
    _rfend = 2.0 * hks->mhub * slow_s;
}


/******************************************************************************/
/* Calculate each phase arrival time (t_0p1s etc) at izr_th H and k  
    combination                                                                */
/******************************************************************************/
void RECORD::Ttps(HKS_INFO *hks, const int& izr)
{
    float rv;
    float mh;
    float slow_p;
    float slow_s;
    float p;
    int ir;
    int iz;
    
    iz = izr / hks->nrv;
    ir = izr % hks->nrv;
        
    mh = hks->mhlb + iz * hks->mhdt;
    rv = hks->rvlb + ir * hks->rvdt;    
    
    p = _ray * 180.0 / (EARTH_RADIUS - mh) / M_PI;
    slow_p = sqrt(1.0 / hks->vp / hks->vp - p * p);
    _tp = mh * slow_p;
    slow_s = sqrt(rv * rv / hks->vp / hks->vp - p * p);
    _ts = mh * slow_s;

    t_0p1s = _ts - _tp;         
    t_2p1s = _ts + _tp; 
    t_1p2s = 2 * _ts;
}


/*
void RECORD::Ttps(HKS_INFO *hks, const int& izr)
{
    int ir;
    int iz;
    
    iz = izr / hks->nrv;
    ir = izr % hks->nrv;
        
    //mh = hks->mhlb + iz * hks->mhdt;
    //rv = hks->rvlb + ir * hks->rvdt;  
    
    if (_tp == NULL) {
        fprintf(stderr, "Error in using _tp in void RECORD::Ttps()!\n\n");
        exit(-1);
    }   
    if (_ts == NULL) {
        fprintf(stderr, "Error in using _ts in void RECORD::Ttps()!\n\n");
        exit(-1);
    }   
    t_0p1s = _ts[izr] - _tp[iz];            
    t_2p1s = _ts[izr] + _tp[iz]; 
    t_1p2s = 2 * _ts[izr];
}
*/


/******************************************************************************/
/* Calculate P and S parts of each phase arrival time (t_0p1s etc)  
    t_0p1s = H * sqrt(1/vs/vs - p*p) - H * sqrt(1/vp/vp - p * p)
             = H * sqrt(k*k/vp/vp - p*p) - H * sqrt(1/vp/vp - p * p)
    t_2p1s = H * sqrt(1/vs/vs - p*p) + H * sqrt(1/vp/vp - p * p)
             = H * sqrt(k*k/vp/vp - p*p) + H * sqrt(1/vp/vp - p * p)
    t_1p2s = 2 * H * sqrt(1/vs/vs - p*p)
             = 2 * H * sqrt(k*k/vp/vp - p*p)
    P part changes with iz (Moho depth).
    tp = H * sqrt(1/vp/vp - p * p)
    S part changes with izr (Moho depth and Vp/Vs ratio)
    ts = H * sqrt(k*k/vp/vp - p*p)                                             */
/******************************************************************************/
/*
void RECORD::TpTs(HKS_INFO *hks)
{
    float rv;
    float mh;
    float slow_p;
    float slow_s;
    float p;
    int ir;
    int iz;
    
    if (_tp != NULL) delete [] _tp;
    _tp = new float[hks->nmh];
    if (_tp == NULL) {
        fprintf(stderr, "Error in new float _tp in void RECORD::TpTs()!\n\n");
        exit(-1);
    }   
    if (_ts != NULL) delete [] _ts;
    _ts = new float[hks->nzr];
    if (_ts == NULL) {
        fprintf(stderr, "Error in new float _ts in void RECORD::TpTs()!\n\n");
        exit(-1);
    }   
    
    for (mh = hks->mhlb, iz = 0; iz < hks->nmh; iz++, mh += hks->mhdt) {
        
        p = _ray * 180.0 / (EARTH_RADIUS - mh) / M_PI;
        slow_p = sqrt(1 / hks->vp / hks->vp - p * p);
        _tp[iz] = mh * slow_p;
        
        for (rv = hks->rvlb, ir = 0; ir < hks->nrv; ir++, rv += hks->rvdt) {    
            slow_s = sqrt(rv * rv / hks->vp / hks->vp - p * p);
            _ts[iz * hks->nrv + ir] = mh * slow_s;
        }
    
    }   
}
*/


/******************************************************************************/
/* Calculate each phase arrival time (t_0p1s etc) at given crustal P-wave
    velocity vp, Moho depth mh, and Vp/Vs ratio rv                             */
/******************************************************************************/
void RECORD::Ttps_rvmh(const float& vp, const float& mh, const float& rv)
{
    float p;
    float slow_p, slow_s;
    float tp, ts;
        
    p = _ray * 180.0 / (EARTH_RADIUS-mh) / M_PI;
    slow_p = sqrt(1/vp/vp - p*p);
    tp = mh * slow_p;
        
    slow_s = sqrt(rv*rv/vp/vp - p*p);
    ts = mh * slow_s;

    t_0p1s = ts - tp;           
    t_2p1s = ts + tp; 
    t_1p2s = 2 * ts;
}

    
