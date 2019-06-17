#include "sac.h"

/*----------------------------------------------------------------------------*/
/*----------------------- class SAC ------------------------------------------*/
/*----------------------------------------------------------------------------*/
SAC_HEADER::SAC_HEADER()
{
	int i, j;
	
	_hdmtr = new SAC_HEADER_MTR;
	
	for (i = 0; i < 14; i++) {
		for (j = 0; j < 5; j++) {
			_hdmtr->headf[i][j] = -12345.0;
		}
	}
	for (i = 0; i < 8; i++) {
		for (j = 0; j < 5; j++) {
			_hdmtr->headi[i][j] = -12345;
		}
	}
	for (i = 0; i < 8; i++) {
		for (j = 0; j < 24; j++) {
			_hdmtr->headc[i][j] = ' ';
		}
	}
	_hdmtr->headc[0][0] = '-';
	_hdmtr->headc[0][1] = '1';
	_hdmtr->headc[0][2] = '2';
	_hdmtr->headc[0][3] = '3';
	_hdmtr->headc[0][4] = '4';
	_hdmtr->headc[0][5] = '5';
	_hdmtr->headc[0][6] = ' ';
	_hdmtr->headc[0][7] = ' ';
	
	_hdmtr->headc[0][8] = '-';
	_hdmtr->headc[0][9] = '1';
	_hdmtr->headc[0][10] = '2';
	_hdmtr->headc[0][11] = '3';
	_hdmtr->headc[0][12] = '4';
	_hdmtr->headc[0][13] = '5';
	_hdmtr->headc[0][14] = ' ';
	_hdmtr->headc[0][15] = ' ';
	_hdmtr->headc[0][16] = ' ';
	_hdmtr->headc[0][17] = ' ';
	_hdmtr->headc[0][18] = ' ';
	_hdmtr->headc[0][19] = ' ';
	_hdmtr->headc[0][20] = ' ';
	_hdmtr->headc[0][21] = ' ';
	_hdmtr->headc[0][22] = ' ';
	_hdmtr->headc[0][23] = ' ';
	
	for (i = 1; i < 8; i++) {
		_hdmtr->headc[i][0] = '-';
		_hdmtr->headc[i][1] = '1';
		_hdmtr->headc[i][2] = '2';
		_hdmtr->headc[i][3] = '3';
		_hdmtr->headc[i][4] = '4';
		_hdmtr->headc[i][5] = '5';
		_hdmtr->headc[i][6] = ' ';
		_hdmtr->headc[i][7] = ' ';
		
		_hdmtr->headc[i][8] = '-';
		_hdmtr->headc[i][9] = '1';
		_hdmtr->headc[i][10] = '2';
		_hdmtr->headc[i][11] = '3';
		_hdmtr->headc[i][12] = '4';
		_hdmtr->headc[i][13] = '5';
		_hdmtr->headc[i][14] = ' ';
		_hdmtr->headc[i][15] = ' ';

		_hdmtr->headc[i][16] = '-';
		_hdmtr->headc[i][17] = '1';
		_hdmtr->headc[i][18] = '2';
		_hdmtr->headc[i][19] = '3';
		_hdmtr->headc[i][20] = '4';
		_hdmtr->headc[i][21] = '5';
		_hdmtr->headc[i][22] = ' ';
		_hdmtr->headc[i][23] = ' ';

	}
	
	_b = -12345.0;
	_e = -12345.0;
	_o = -12345.0;
	_delta = -12345.0;
	_scale = -12345.0;
	_stla = -12345.0;
	_stlo = -12345.0;
	_stel = -12345.0;
	_stdp = -12345.0;
	_cmpaz = -12345.0;
	_cmpinc = -12345.0;
	_evla = -12345.0;
	_evlo = -12345.0;
	_evel = -12345.0;
	_evdp = -12345.0;
	_mag = -12345.0;
	_user0 = -12345.0;
	_user1 = -12345.0;
	_user2 = -12345.0;
	_user3 = -12345.0;
	_user4 = -12345.0;
	_user5 = -12345.0;
	_user6 = -12345.0;
	_user7 = -12345.0;
	_user8 = -12345.0;
	_user9 = -12345.0;
	_uusd1304 = -12345.0;
	_uusd1305 = -12345.0;
	
	_nzyear = -12345;
	_nzjday = -12345;
	_nzhour = -12345;
	_nzmin = -12345;
	_nzsec = -12345;
	_nzmsec = -12345;
	_nvhdr = -12345;
	_nevid = -12345;
	_norid = -12345;
	_nwfid = -12345;
	_iftype = -12345;
	_idep = -12345;
	_iztype = -12345;
	_ievtyp = -12345;
	_isynth = -12345;
	_leven = -12345;
	_lovrok = -12345;
	_lpspol = -12345;
	_lcalda = -12345;
	_imagtyp = -12345;
	_npts = -12345;
	
	_nvhdr = 6;
	_hdmtr->headi[1][1] = 6;

	_kstnm = new char[9];
	strcpy(_kstnm, "-12345");
	_kcmpnm = new char[9];
	strcpy(_kcmpnm, "-12345");
	_khole = new char[9];
	strcpy(_khole, "-12345");
	_knetwk = new char[9];
	strcpy(_knetwk, "-12345");
	_kdatrd = new char[9];
	strcpy(_kdatrd, "-12345");
	_kuser0 = new char[9];
	strcpy(_kuser0, "-12345");
	_kuser1 = new char[9];
	strcpy(_kuser1, "-12345");
	_kuser2 = new char[9];
	strcpy(_kuser2, "-12345");

}

/*----------------------------------------------------------------------------*/
SAC_HEADER::~SAC_HEADER()
{
	delete [] _kstnm;
	delete [] _khole;
	delete [] _knetwk;
	delete [] _kcmpnm;
	delete [] _kdatrd;
	delete [] _kuser0;
	delete [] _kuser1;
	delete [] _kuser2;
	delete _hdmtr;
	
	_kstnm = NULL;
	_khole = NULL;
	_knetwk = NULL;
	_kcmpnm = NULL;
	_kdatrd = NULL;
	_hdmtr = NULL;
	_kuser0 = NULL;	
	_kuser1 = NULL;	
	_kuser2 = NULL;	
}	
			
/*----------------------------------------------------------------------------*/
void SAC_HEADER::cpsachdr(SAC_HEADER* sachdr)
/*----------------------------------------------------------------------------*/
/* copy all non -12345 field from sachdr 
	original values are forced to changed if they are set 
	user defined and unused fields are not considered                          */
/*----------------------------------------------------------------------------*/
{	
	if (sachdr->_b != -12345.0) set_b(sachdr->_b);
	if (sachdr->_e != -12345.0) set_e(sachdr->_e);
	if (sachdr->_o != -12345.0) set_o(sachdr->_o);
	if (sachdr->_delta != -12345.0) set_delta(sachdr->_delta);
	if (sachdr->_scale != -12345.0) set_scale(sachdr->_scale);

	if (sachdr->_stla != -12345.0)
		set_sta(_kstnm, sachdr->_stla, _stlo, _stel, _stdp);		
	if (sachdr->_stlo != -12345.0)
		set_sta(_kstnm, _stla, sachdr->_stlo, _stel, _stdp); 
	if (sachdr->_stel != -12345.0)
		set_sta(_kstnm, _stla, _stlo, sachdr->_stel, _stdp);	
	if (sachdr->_stdp != -12345.0)		 
		set_sta(_kstnm, _stla, _stlo, _stel, sachdr->_stdp);
		 
	if (sachdr->_cmpaz != -12345.0) set_cmp(sachdr->_cmpaz, _cmpinc);	
	if (sachdr->_cmpinc != -12345.0) set_cmp(_cmpaz, sachdr->_cmpinc);

	if (sachdr->_evla != -12345.0)
		set_evt(sachdr->_evla, _evlo, _evel, _evdp, _mag, _imagtyp);	
	if (sachdr->_evlo != -12345.0)
		set_evt(_evla, sachdr->_evlo, _evel, _evdp, _mag, _imagtyp); 
	if (sachdr->_evel != -12345.0)
		set_evt(_evla, _evlo, sachdr->_evel, _evdp, _mag, _imagtyp);	
	if (sachdr->_evdp != -12345.0)	
		set_evt(_evla, _evlo, _evel, sachdr->_evdp, _mag, _imagtyp);		
	if (sachdr->_mag != -12345.0) set_mag(sachdr->_mag);
	if (sachdr->_imagtyp != -12345) set_imagtyp(sachdr->_imagtyp);
		

	if (sachdr->_nzyear != -12345 && sachdr->_nzjday != -12345 && 
		 sachdr->_nzhour != -12345 && sachdr->_nzmin != -12345 &&
		 sachdr->_nzsec != -12345 && sachdr->_nzmsec != -12345)
		set_reftime(sachdr->_nzyear, sachdr->_nzjday, sachdr->_nzhour, 
						sachdr->_nzmin, sachdr->_nzsec, sachdr->_nzmsec);

	if (sachdr->_nevid != -12345) set_css(sachdr->_nevid, _norid, _nwfid);
	if (sachdr->_norid != -12345) set_css(_nevid, sachdr->_norid, _nwfid);
	if (sachdr->_nwfid != -12345) set_css(_nevid, _norid, sachdr->_nwfid);
	 
	if (sachdr->_npts != -12345) set_npts(sachdr->_npts);
	if (sachdr->_iftype != -12345) set_iftype(sachdr->_iftype);
	if (sachdr->_idep != -12345) set_idep(sachdr->_idep);
	if (sachdr->_iztype != -12345) set_iztype(sachdr->_iztype);
	if (sachdr->_ievtyp != -12345)	set_ievtyp(sachdr->_ievtyp);
	if (sachdr->_isynth != -12345) set_isynth(sachdr->_isynth);
	if (sachdr->_leven != -12345) set_leven(sachdr->_leven);
	if (sachdr->_lpspol != -12345)	set_lpspol(sachdr->_lpspol);
	if (sachdr->_lovrok != -12345)	set_lovrok(sachdr->_lovrok);
	if (sachdr->_lcalda != -12345) set_lcalda(sachdr->_lcalda);

	if (strcmp(sachdr->_knetwk, "-12345")) set_knetwk(sachdr->_knetwk);
	if (strcmp(sachdr->_kstnm, "-12345")) set_kstnm(sachdr->_kstnm);
	if (strcmp(sachdr->_kcmpnm, "-12345")) set_kcmpnm(sachdr->_kcmpnm);
	if (strcmp(sachdr->_kdatrd, "-12345")) set_kdatrd(sachdr->_kdatrd);
	if (strcmp(sachdr->_khole, "-12345"))	set_khole(sachdr->_khole);

}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_b(float b)	
{
	_b = b;
	_hdmtr->headf[1][0] = b;
}
								
/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_e(float e)	
{
	_e = e;
	_hdmtr->headf[1][1] = e;
}
								
/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_o(float o)	
{
	_o = o;
	_hdmtr->headf[1][2] = o;
}
								
/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_delta(float delta)
{
	_delta = delta;
	_hdmtr->headf[0][0] = delta;
}	

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_scale(float scale)
{
	_scale = scale;
	_hdmtr->headf[0][3] = scale;
}	

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_npts(int npts)
{
	_npts = npts;
	_hdmtr->headi[1][4] = npts;
}	

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_reftime(int nzyear, int nzjday, 
										int nzhour, int nzmin, int nzsec, int nzmsec)
{
	_nzyear = nzyear;
	_nzjday = nzjday;
	_nzhour = nzhour;
	_nzmin = nzmin;
	_nzsec = nzsec;
	_nzmsec = nzmsec;
	
	_hdmtr->headi[0][0] = nzyear;
	_hdmtr->headi[0][1] = nzjday;
	_hdmtr->headi[0][2] = nzhour;
	_hdmtr->headi[0][3] = nzmin;
	_hdmtr->headi[0][4] = nzsec;
	_hdmtr->headi[1][0] = nzmsec;
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_sta(char *stnm, float stla, float stlo, float stel, float stdp) 
{
	strcpy(_kstnm, stnm);
	int lstnm = (strlen(_kstnm) <8 ) ? strlen(_kstnm) : 8;
	for (int j = 0; j < lstnm; j++) {
		_hdmtr->headc[0][j] = _kstnm[j];
	}	
	
	_stla = stla;
	_hdmtr->headf[6][1] = stla;
	_stlo = stlo;
	_hdmtr->headf[6][2] = stlo;	
	_stel = stel;
	_hdmtr->headf[6][3] = stel;
	_stdp = stdp;	
	_hdmtr->headf[6][4] = stdp;

}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_stdp(float stdp) 
{
	_stdp = stdp;	
	_hdmtr->headf[6][4] = stdp;

}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_evt(float evla, float evlo, float evel, float evdp, 
									float mag, int imagtyp)
{
	_evla = evla;
	_evlo = evlo;
	_evel = evel;
	_evdp = evdp;
	_mag = mag;
	if (imagtyp == IMB|| imagtyp == IMS || imagtyp == IML || imagtyp == IMW || 
			imagtyp == IMD || imagtyp == IMX) {
		_imagtyp = imagtyp;
	}
	else {
		_imagtyp = -12345;
	}
	
	_hdmtr->headf[7][0] = _evla;
	_hdmtr->headf[7][1] = _evlo;
	_hdmtr->headf[7][2] = _evel;
	_hdmtr->headf[7][3] = _evdp;
	_hdmtr->headf[7][4] = _mag;
	_hdmtr->headi[5][0] = _imagtyp;
	
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_mag(float mag)
{
	_mag = mag;
	_hdmtr->headf[7][4] = _mag;
	
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_imagtyp(int type)
{
	if (type == IMB|| type == IMS || type == IML || type == IMW || type == IMD ||
			type == IMX) {
		_imagtyp = type;
	}
	else {
		_imagtyp = -12345;
	}
	_hdmtr->headi[5][0] = _imagtyp;
	
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_cmp(float cmpaz, float cmpinc)
{
	_cmpaz = cmpaz;
	_cmpinc = cmpinc;
	_hdmtr->headf[11][2] = cmpaz;
	_hdmtr->headf[11][3] = cmpinc;
}	

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_css(int nevid, int norid, int nwfid)
{
	_nevid = nevid;
	_norid = norid;
	_nwfid = nwfid;
	_hdmtr->headi[1][2] = norid;
	_hdmtr->headi[1][3] = nevid;
	_hdmtr->headi[2][1] = nwfid;
	
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_iftype(int type)
{
	if (type == ITIME || type == IRLIM || type == IAMPH || type == IXY || type == IXYZ) {
		_iftype = type;
	}
	else {
		_iftype = -12345;
	}
	_hdmtr->headi[3][0] = _iftype;
	
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_idep(int idep)
{
	if (idep == IUNKN || idep == IDISP || idep == IVEL || idep == IVOLTS || idep == IACC) {
		_idep = idep;
	}
	else {
		_idep = -12345;
	}
	_hdmtr->headi[3][1] = _idep;
	
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_iztype(int type)
{
	if (type == IUNKN || type == IB || type == IDAY || type == IO || type == IA ||
			type == IT0 || type == IT1 || type == IT2 || type == IT3 || type == IT4 ||
			type == IT5 || type == IT6 || type == IT7 || type == IT8 || type == IT9) {
		_iztype = type;
	}
	else {
		_iztype = -12345;
	}
	_hdmtr->headi[3][2] = _iztype;
	
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_ievtyp(int type)
{
	if (type == IUNKN || type == INUCL || type == IPREN || type == IPOSTN || type == IQUAKE ||
			type == IPREQ || type == IPOSTQ || type == ICHEM || type == IQB || type == IQB1 ||
			type == IQB2 || type == IQBX || type == IQMT || type == IEQ || type == IEQ1 ||
			type == IEQ2 || type == IME || type == IEX || type == INU || type == INC || type == IO_ ||
			type == IL || type == IR || type == IT || type == IU || type == IOTHER) {
		_ievtyp = type;
	}
	else {
		_ievtyp = -12345;
	}
	_hdmtr->headi[4][2] = _ievtyp;
	
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_isynth(int type)
{
	_isynth = type;
	_hdmtr->headi[4][4] = type;
	
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_leven(int opt)
{
	if (opt == TRUE || opt == FALSE)
		_leven = opt;
	else
		_leven = -12345;
	_hdmtr->headi[7][0] = _leven;
}
		 
/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_lovrok(int opt)
{
	if (opt == TRUE || opt == FALSE)
		_lovrok = opt;
	else
		_lovrok = -12345;
	_hdmtr->headi[7][2] = _lovrok;
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_lpspol(int opt)
{
	if (opt == TRUE || opt == FALSE)
		_lpspol = opt;
	else
		_lpspol = -12345;
	_hdmtr->headi[7][1] = _lpspol;
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_lcalda(int opt)
{
	if (opt == TRUE || opt == FALSE)
		_lcalda = opt;
	else
		_lcalda = -12345;
	_hdmtr->headi[7][3] = _lcalda;
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_khole(char* khole)
{
	int i;
	int strleng = strlen(khole);	
	
	if (strleng == 0 || khole[0] == '\0') {
		_hdmtr->headc[1][0] = '-';
		_hdmtr->headc[1][1] = '1';
		_hdmtr->headc[1][2] = '2';
		_hdmtr->headc[1][3] = '3';
		_hdmtr->headc[1][4] = '4';
		_hdmtr->headc[1][5] = '5';
		_hdmtr->headc[1][6] = ' ';
		_hdmtr->headc[1][7] = ' ';
	}
	else if (strleng > 0 && strleng < 8) {
		for (i = 0; i < strleng; i++) {
			_hdmtr->headc[1][i] = khole[i];
			//printf("khole mtr =[%c]\n", _hdmtr[1][i]); 							
		}	
		for (i = strleng; i < 8; i++)
			_hdmtr->headc[1][i] = ' ';
		//printf("khole mtr =[%s]\n", _hdmtr->headc[1]); 							
	}		
	else {
		for (i = 0; i < 8; i++)
			_hdmtr->headc[1][i] = khole[i];
	}
	
	for (i = 0; i < 8; i++)	_khole[i] =	_hdmtr->headc[1][i];
	
	_khole[8] = '\0';
	//printf("khole =[%s]\n", _khole); 						

}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_knetwk(char* net)
{
	int i;
	int strleng = strlen(net);	
	
	if (strleng == 0 || net[0] == '\0') {
		_hdmtr->headc[7][0] = '-';
		_hdmtr->headc[7][1] = '1';
		_hdmtr->headc[7][2] = '2';
		_hdmtr->headc[7][3] = '3';
		_hdmtr->headc[7][4] = '4';
		_hdmtr->headc[7][5] = '5';
		_hdmtr->headc[7][6] = ' ';
		_hdmtr->headc[7][7] = ' ';
	}
	else if (strleng > 0 && strleng < 8) {
		for (i = 0; i < strleng; i++) {
			_hdmtr->headc[7][i] = net[i];
			//printf("khole mtr =[%c]\n", _hdmtr[1][i]); 							
		}	
		for (i = strleng; i < 8; i++)
			_hdmtr->headc[7][i] = ' ';
		//printf("khole mtr =[%s]\n", _hdmtr->headc[1]); 							
	}		
	else {
		for (i = 0; i < 8; i++)
			_hdmtr->headc[7][i] = net[i];
	}
	
	for (i = 0; i < 8; i++)	_knetwk[i] = _hdmtr->headc[7][i];
	
	_knetwk[8] = '\0';
	//printf("khole =[%s]\n", _khole); 						

}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_kstnm(char* stacode)
{
	int i;
	int strleng = strlen(stacode);	
	
	if (strleng == 0 || stacode[0] == '\0') {
		_hdmtr->headc[0][0] = '-';
		_hdmtr->headc[0][1] = '1';
		_hdmtr->headc[0][2] = '2';
		_hdmtr->headc[0][3] = '3';
		_hdmtr->headc[0][4] = '4';
		_hdmtr->headc[0][5] = '5';
		_hdmtr->headc[0][6] = ' ';
		_hdmtr->headc[0][7] = ' ';
	}
	else if (strleng > 0 && strleng < 8) {
		for (i = 0; i < strleng; i++) {
			_hdmtr->headc[0][i] = stacode[i];
			//printf("khole mtr =[%c]\n", _hdmtr[1][i]); 							
		}	
		for (i = strleng; i < 8; i++)
			_hdmtr->headc[0][i] = ' ';
		//printf("khole mtr =[%s]\n", _hdmtr->headc[1]); 							
	}		
	else {
		for (i = 0; i < 8; i++)
			_hdmtr->headc[0][i] = stacode[i];
	}
	
	for (i = 0; i < 8; i++)	_kstnm[i] =	_hdmtr->headc[0][i];
	
	_kstnm[8] = '\0';
	//printf("kshole =[%s]\n", _khole); 						

}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_kcmpnm(char* comp)
{
	int i;
	int strleng = strlen(comp);	
	
	if (strleng == 0 || comp[0] == '\0') {
		_hdmtr->headc[6][16] = '-';
		_hdmtr->headc[6][17] = '1';
		_hdmtr->headc[6][18] = '2';
		_hdmtr->headc[6][19] = '3';
		_hdmtr->headc[6][20] = '4';
		_hdmtr->headc[6][21] = '5';
		_hdmtr->headc[6][22] = ' ';
		_hdmtr->headc[6][23] = ' ';
	}
	else if (strleng > 0 && strleng < 8) {
		for (i = 0; i < strleng; i++) {
			_hdmtr->headc[6][i+16] = comp[i];
			//printf("khole mtr =[%c]\n", _hdmtr[1][i]); 							
		}	
		for (i = strleng; i < 8; i++)
			_hdmtr->headc[6][i+16] = ' ';
		//printf("khole mtr =[%s]\n", _hdmtr->headc[1]); 							
	}		
	else {
		for (i = 0; i < 8; i++)
			_hdmtr->headc[6][i+16] = comp[i];
	}
	
	for (i = 0; i < 8; i++)	_kcmpnm[i] = _hdmtr->headc[6][i+16];
	
	_kcmpnm[8] = '\0';
	//printf("kcmpnm =[%s]\n", _kcmpnm); 						

}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_kdatrd(char* kdatrd)
{
	int i;
	int strleng = strlen(kdatrd);	
	
	if (strleng == 0 || kdatrd[0] == '\0') {
		_hdmtr->headc[7][8] = '-';
		_hdmtr->headc[7][9] = '1';
		_hdmtr->headc[7][10] = '2';
		_hdmtr->headc[7][11] = '3';
		_hdmtr->headc[7][12] = '4';
		_hdmtr->headc[7][13] = '5';
		_hdmtr->headc[7][14] = ' ';
		_hdmtr->headc[7][15] = ' ';
	}
	else if (strleng > 0 && strleng < 8) {
		for (i = 0; i < strleng; i++) {
			_hdmtr->headc[7][i+8] = kdatrd[i];
			//printf("khole mtr =[%c]\n", _hdmtr[1][i]); 							
		}	
		for (i = strleng; i < 8; i++)
			_hdmtr->headc[7][i+8] = ' ';
		//printf("khole mtr =[%s]\n", _hdmtr->headc[1]); 							
	}		
	else {
		for (i = 0; i < 8; i++)
			_hdmtr->headc[7][i+8] = kdatrd[i];
	}
	
	for (i = 0; i < 8; i++)	_kdatrd[i] = _hdmtr->headc[7][i+8];
	
	_kdatrd[8] = '\0';
	//printf("khole =[%s]\n", _khole); 						

}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_userf(int i, float value)
{
	if (i == 0) {
		_hdmtr->headf[8][0] = value;
		_user0 = value;
	}	
	else if (i == 1) {
		_hdmtr->headf[8][1] = value;
		_user1 = value;
	}
	else if (i == 2) {
		_hdmtr->headf[8][2] = value;
		_user2 = value;
	}	
	else if (i == 3) {
		_hdmtr->headf[8][3] = value;
		_user3 = value;
	}	
	else if (i == 4) {
		_hdmtr->headf[8][4] = value;
		_user4 = value;
	}	
	else if (i == 5) {
		_hdmtr->headf[9][0] = value;
		_user5 = value;
	}	
	else if (i == 6) {
		_hdmtr->headf[9][1] = value;
		_user6 = value;
	}
	else if (i == 7) {
		_hdmtr->headf[9][2] = value;
		_user7 = value;
	}	
	else if (i == 8) {
		_hdmtr->headf[9][3] = value;
		_user8 = value;
	}
	else if (i == 9) {
		_hdmtr->headf[9][4] = value;
		_user9 = value;
	}
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_unusedf(int i, int j, float value)
{
	_hdmtr->headf[i-1][j-1] = value;
	if (i == 13 && j == 4) _uusd1304 = value;
	else if (i == 13 && j == 5) _uusd1305 = value;
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_unusedi(int i, int j, int value)
{
	_hdmtr->headi[i-1][j-1] = value;
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::set_userk(int index, char* value)
{
	int i = 0;
	int strleng = strlen(value);
	
	if (index == 0) {	
		
		if (strleng == 0 || value[0] == '\0') {
			_hdmtr->headc[5][16] = '-';
			_hdmtr->headc[5][17] = '1';
			_hdmtr->headc[5][18] = '2';
			_hdmtr->headc[5][19] = '3';
			_hdmtr->headc[5][20] = '4';
			_hdmtr->headc[5][21] = '5';
			_hdmtr->headc[5][22] = ' ';
			_hdmtr->headc[5][23] = ' ';
		}
		else if (strleng > 0 && strleng < 8) {
			for (i = 0; i < strleng; i++) {
				_hdmtr->headc[5][i+16] = value[i];
			}	
			for (i = strleng; i < 8; i++)
				_hdmtr->headc[5][i+16] = ' ';
		}		
		else {
			for (i = 0; i < 8; i++)
				_hdmtr->headc[5][i+16] = value[i];
		}
	
		for (i = 0; i < 8; i++)	_kuser0[i] = _hdmtr->headc[5][i+16];
		_kuser0[8] = '\0';

	}	
	else if (index == 1) {
	
		if (strleng == 0 || value[0] == '\0') {
			_hdmtr->headc[6][0] = '-';
			_hdmtr->headc[6][1] = '1';
			_hdmtr->headc[6][2] = '2';
			_hdmtr->headc[6][3] = '3';
			_hdmtr->headc[6][4] = '4';
			_hdmtr->headc[6][5] = '5';
			_hdmtr->headc[6][6] = ' ';
			_hdmtr->headc[6][7] = ' ';
		}
		else if (strleng > 0 && strleng < 8) {
			for (i = 0; i < strleng; i++) {
				_hdmtr->headc[6][i] = value[i];
			}	
			for (i = strleng; i < 8; i++)
				_hdmtr->headc[6][i] = ' ';
		}		
		else {
			for (i = 0; i < 8; i++)
				_hdmtr->headc[6][i] = value[i];
		}
	
		for (i = 0; i < 8; i++)	_kuser1[i] = _hdmtr->headc[6][i];	
		_kuser1[8] = '\0';
	
	}
	else if (index == 2) {

		if (strleng == 0 || value[0] == '\0') {
			_hdmtr->headc[6][8] = '-';
			_hdmtr->headc[6][9] = '1';
			_hdmtr->headc[6][10] = '2';
			_hdmtr->headc[6][11] = '3';
			_hdmtr->headc[6][12] = '4';
			_hdmtr->headc[6][13] = '5';
			_hdmtr->headc[6][14] = ' ';
			_hdmtr->headc[6][15] = ' ';
		}
		else if (strleng > 0 && strleng < 8) {
			for (i = 0; i < strleng; i++) {
				_hdmtr->headc[6][i+8] = value[i];
			}	
			for (i = strleng; i < 8; i++)
				_hdmtr->headc[6][i+8] = ' ';
		}		
		else {
			for (i = 0; i < 8; i++)
				_hdmtr->headc[6][i+8] = value[i];
		}
	
		for (i = 0; i < 8; i++)	_kuser2[i] = _hdmtr->headc[6][i+8];	
		_kuser2[8] = '\0';

	}	
}

/*----------------------------------------------------------------------------*/
void SAC_HEADER::rsacheader(char *kname)
{
	FILE* rsacp;
	int ndumy;
	int i, j;
	char cbuf[280];

	float float_swap(char cbuf[]);
	long long_swap(char cbuf[]);	
	long long_join(char cbuf[]);	 
	int int_swap(char cbuf[]);	
	int int_join(char cbuf[]);	 
	 
	rsacp = fopen(kname, "rb");
 	if (rsacp == NULL) {	
   	fprintf(stderr, "\nError in rsacheader: Cannot open sac header file %s!!\n\n",kname);
		fclose(rsacp);
		exit(0);
 	}
	
	if (fread(_hdmtr, sizeof(SAC_HEADER_MTR), 1, rsacp) != 1) {
	 	fprintf(stderr, "\nError in rsacheader: Reading sac header of file %s !! \n\n", kname);
  		fclose(rsacp);
		exit(0);
	}		

/*----------------------------------------------------------------------------*/
/* Test if need swap                                                          */
/*----------------------------------------------------------------------------*/
	rewind(rsacp);
	ndumy = fread(cbuf, 4 * FLOAT_ROWS * FLOAT_COLS * sizeof(char), 1, rsacp);
	ndumy = fread(cbuf, 4 * LONG_ROWS * LONG_COLS * sizeof(char), 1, rsacp);

	_nvhdr = (int) int_join(cbuf + 4 * 6);
	if (_nvhdr > 0 && _nvhdr < 10) {
		is_swap = 0;
	}
	else {	
		rewind(rsacp);

/*----------------------------------------------------------------------------*/
/* Need swap                                                                  */
/*----------------------------------------------------------------------------*/
		ndumy = fread(cbuf, 4 * FLOAT_ROWS * FLOAT_COLS * sizeof(char), 1, rsacp);
		for (i = 0; i < FLOAT_ROWS; i++) {
			for (j = 0; j < FLOAT_COLS; j++) {
				_hdmtr->headf[i][j] = float_swap(cbuf + 4 * (j + i * FLOAT_COLS));
			}
		}			
		ndumy = fread(cbuf, 4 * LONG_ROWS * LONG_COLS * sizeof(char), 1, rsacp);
		for (i = 0; i < LONG_ROWS; i++) {
			for (j = 0; j < LONG_COLS; j++) {
				_hdmtr->headi[i][j] = int_swap(cbuf + 4 * (j + i * LONG_COLS));
			}	
		}	
/*----------------------------------------------------------------------------*/
/* Chars do not need swap                                                     */
/*----------------------------------------------------------------------------*/
		is_swap = 1;
	}

	fclose(rsacp);
	
	_b = _hdmtr->headf[1][0];
	_e = _hdmtr->headf[1][1];
	_o = _hdmtr->headf[1][2];
	_delta = _hdmtr->headf[0][0];
	_scale = _hdmtr->headf[0][3];
	_npts = _hdmtr->headi[1][4];
	
	_nzyear = _hdmtr->headi[0][0];
	_nzjday = _hdmtr->headi[0][1];
	_nzhour = _hdmtr->headi[0][2];
	_nzmin = _hdmtr->headi[0][3];
	_nzsec = _hdmtr->headi[0][4];
	_nzmsec = _hdmtr->headi[1][0];

	_nvhdr = _hdmtr->headi[1][1];
	_norid = _hdmtr->headi[1][2];
	_nevid = _hdmtr->headi[1][3];
	_nwfid = _hdmtr->headi[2][1];
	_iftype = _hdmtr->headi[3][0];
	_idep = _hdmtr->headi[3][1];
	_iztype = _hdmtr->headi[3][2];
	_ievtyp = _hdmtr->headi[4][2];
	_isynth = _hdmtr->headi[4][4];
	_leven = _hdmtr->headi[7][0];
	_lovrok = _hdmtr->headi[7][2];
	_lcalda = _hdmtr->headi[7][3];

	_stla = _hdmtr->headf[6][1];
	_stlo = _hdmtr->headf[6][2];	
	_stel = _hdmtr->headf[6][3];
	_stdp = _hdmtr->headf[6][4];
	_cmpaz = _hdmtr->headf[11][2];
	_cmpinc = _hdmtr->headf[11][3];		

	_evla = _hdmtr->headf[7][0];
	_evlo = _hdmtr->headf[7][1];
	_evel = _hdmtr->headf[7][2];
	_evdp = _hdmtr->headf[7][3];
	_mag = _hdmtr->headf[7][4];
	_imagtyp = _hdmtr->headi[5][0];
	
	for (int j = 0; j < 8; j++) {
		_kstnm[j] = _hdmtr->headc[0][j];
		_kuser0[j] = _hdmtr->headc[5][j+16];
		_kuser1[j] = _hdmtr->headc[6][j];
		_kuser2[j] = _hdmtr->headc[6][j+8];
	}
	_kstnm[8] = '\0';
	_kuser0[8] = '\0';
	_kuser1[8] = '\0';
	_kuser2[8] = '\0';
	
/*----------------------------------------------------------------------------*/
/*	Add miscellaneous fields                                                   */
/*----------------------------------------------------------------------------*/
	_user0 = _hdmtr->headf[8][0];
	_user1 = _hdmtr->headf[8][1];
	_user2 = _hdmtr->headf[8][2];
	_user3 = _hdmtr->headf[8][3];
	_user4 = _hdmtr->headf[8][4];
	_user5 = _hdmtr->headf[9][0];
	_user6 = _hdmtr->headf[9][1];
	_user7 = _hdmtr->headf[9][2];
	_user8 = _hdmtr->headf[9][3];
	_user9 = _hdmtr->headf[9][4];
	
}

/*----------------------------------------------------------------------------*/
/* Write sac header into a file.
	The data can be added by other programs.                                   */		
/*----------------------------------------------------------------------------*/
void SAC_HEADER::wsacheader(char *kname)
{
	FILE* wsacp;

	wsacp = fopen(kname, "wb");
	if (wsacp == NULL) {
		fprintf(stderr, "\nError in wsacheader: Cannot open file %s when write sac header!!\n\n", kname);
		fclose(wsacp);
		exit(0);
	}

	// writing header
	if (fwrite(_hdmtr, sizeof(SAC_HEADER_MTR), 1, wsacp) != 1) {
		fprintf(stderr, "\nError in wsacheader: Writing sac header for file %s!! \n\n", kname);
		fclose(wsacp);
		exit(0);
	}
	  	
	fclose(wsacp);

}


/******************************************************************************/	
int rsac0(char *kname, float* data, int& npts, float& delta, float& b)
{
	FILE* rsacp;
	SAC_HEADER* sachd;
	int j;
	char cbuf[4];
	
	float float_swap(char cbuf[]);
	 
	sachd = new SAC_HEADER;
	sachd->rsacheader(kname);
	
	npts = sachd->_npts; 
	delta = sachd->_delta;
	b = sachd->_b;
	
	rsacp = fopen(kname, "rb");
 	if (rsacp == NULL) {
   	fprintf(stderr, "\nError in rsac0: Cannot open sac data file %s!!\n\n",kname);
		fclose(rsacp);
   	return 0;
 	}
		
	if (fread(sachd->_hdmtr, sizeof(SAC_HEADER_MTR), 1, rsacp) != 1) {
   	fprintf(stderr, "\nError in rsac0: Reading sac header of file %s !! \n\n", kname);
   	fclose(rsacp);
   	return 0;
 	}		
		
	
	if (sachd->is_swap == 1) {
		for (j = 0; j < npts; j++) {
			if (fread(cbuf, sizeof(char), 4, rsacp) != 4) {
	   		fprintf(stderr, "\nError in rsac0: Reading swaping sac data of file %s  != 4 at %d !! \n\n", kname, j);
   			fclose(rsacp);
   			return 0;
			}
			data[j] = (float) float_swap(cbuf);
		}
	}
	else {	
		if (fread(data, sizeof(float), npts, rsacp) != (unsigned int) npts) {
   		fprintf(stderr, "\nError in rsac0: Reading sac data of file %s  Npts = %d !! \n\n", kname, npts);
   		fclose(rsacp);
   		return 0;
		}
	}			

	fclose(rsacp);
	delete sachd;
	sachd = NULL;
	return 1;
}	

/*----------------------------------------------------------------------------*/
/*	Read part of sac data [spt, spt+nptsr-1], including point "spt".
	"spt" is accounted from the 1st sample of sac file (zero-shift).
	"spt" must be >= 0                              
	If "spt + nptsr" > sac file real length, read actrual length
	Return 1 if succesful.
	Note the difference of between rsac1, rsac2, rsac4:
	rsac1: The points [spt, spt+nptsr-1] in sac file is put in data[0, nptsr-1]
			 If spt+nptsr > npts(sac length), return error.
			 The return value is 1(succesful) or 0(error).
			 Used for read two sac files and process them, which require 
			 the same length of them.
	rsac2: The points [spt, spt+nptsr-1] in sac file is put in data[0, nptsr-1]
			 If spt+nptsr > npts(sac length), nptsr = npts - spt.
			 The return value is the actrually read number(succesful) or 0(error).
			 Used for cut sac file, or only process it as a new sac file.
	rsac4: The points [spt, spt+nptsr-1] in sac file is put in 
			 data[spt, spt+nptsr-1].
			 The other segments in "data" is unchanged.
			 Before calling, make sure zero "data" and assign enough "data" size.
			 If spt+nptsr > npts(sac length), nptsr = npts - spt.
			 The return value is the actrually read number(succesful) or 0(error).
			 Used for cut sac file, but keep the same time sequence of 
			 original sac file, like in H-k search                               */
/*----------------------------------------------------------------------------*/
int rsac1(char *kname, float* data, const int& spt, const int& nptsr)
{
	FILE* rsacp;
	SAC_HEADER* sachd;
	int ndumy;
	int npts;               // Points of data in SAC header 
	unsigned int nptsb;     // Returned points of fread
	int j;
	char cbuf[4];
	
	float float_swap(char cbuf[]);
	 
	if (spt < 0) {
   	fprintf(stderr, "\nError in rsac1: Starting point (= %d) < 0 !! \n\n", spt);
   	return 0;
 	}		

	sachd = new SAC_HEADER;
	sachd->rsacheader(kname);
	npts = sachd->_npts; 
	 
	if (spt >= npts) {
   	fprintf(stderr, "\nError in rsac1: Starting point (= %d) >= npts (= %d) !! \n\n", spt, npts);
   	return 0;
 	}		

	rsacp = fopen(kname, "rb");
 	if (rsacp == NULL) {
   	fprintf(stderr, "\nError in rsac1: Cannot open sac data file %s!!\n\n",kname);
		fclose(rsacp);
   	return 0;
 	}
	
	if (fread(sachd->_hdmtr, sizeof(SAC_HEADER_MTR), 1, rsacp) != 1) {
   	fprintf(stderr, "\nError in rsac1: Reading sac header of file %s !! \n\n", kname);
   	fclose(rsacp);
   	return 0;
 	}		
	
	if ((spt + nptsr) > npts) {
		fprintf(stderr, "\nError in rsac1 %s:\n", kname); 
		fprintf(stderr, "Start pts = %d, supposed read nptsr = %d\n", spt, nptsr);
   	fprintf(stderr, "Total (pts+nptsr) = %d > SAC data points = %d !! \n\n",
					 spt+nptsr, npts);
   	fclose(rsacp);					 
   	return 0;
 	}		
				
	if (fseek(rsacp, spt*sizeof(float), SEEK_CUR) != 0) {
		fprintf(stderr, "\nError in rsac1: Seek sac data in file %s!! \n\n", kname);
	}
	
	if (sachd->is_swap == 1) {
		for (j = 0; j < nptsr; j++) {
			ndumy = fread(cbuf, sizeof(char), 4, rsacp);
			data[j] = (float) float_swap(cbuf);
		}	
	}
	else {
		nptsb = fread(data, sizeof(float), nptsr, rsacp);
		if (nptsb != (unsigned int) nptsr) {
			fprintf(stderr, "\nError in rsac1 %s:\n", kname); 
			fprintf(stderr, "Start pts = %d, supposed read pts = %d, sac pts = %d\n", spt, nptsr, npts);
			fprintf(stderr, "Actrually reading sac pts = %d !! \n\n", nptsb);
			fclose(rsacp);
			return(0);
		}		
	}
		
	fclose(rsacp);
	delete sachd;
	sachd = NULL;
	return 1;
}	

/*----------------------------------------------------------------------------*/
/*	Read part of sac data [spt, spt+nptsr-1], including point "spt".
	"spt" is accounted from the 1st sample of sac file (zero-shift).
	"spt" must be >= 0                              
	If "spt + nptsr" > sac file real length, read actrual length
	Return the actrually read number of data
	Note the difference of between rsac1, rsac2, rsac4:
	rsac1: The points [spt, spt+nptsr-1] in sac file is put in data[0, nptsr-1]
			 If spt+nptsr > npts(sac length), return error.
			 The return value is 1(succesful) or 0(error).
			 Used for read two sac files and process them, which require 
			 the same length of them.
	rsac2: The points [spt, spt+nptsr-1] in sac file is put in data[0, nptsr-1]
			 If spt+nptsr > npts(sac length), nptsr = npts - spt.
			 The return value is the actrually read number(succesful) or 0(error).
			 Used for cut sac file, or only process it as a new sac file.
	rsac4: The points [spt, spt+nptsr-1] in sac file is put in 
			 data[spt, spt+nptsr-1].
			 The other segments in "data" is unchanged.
			 Before calling, make sure zero "data" and assign enough "data" size.
			 If spt+nptsr > npts(sac length), nptsr = npts - spt.
			 The return value is the actrually read number(succesful) or 0(error).
			 Used for cut sac file, but keep the same time sequence of 
			 original sac file, like in H-k search                               */
/*----------------------------------------------------------------------------*/
int rsac2(char *kname, float* data, const int& spt, const int& nptsr)
{
	FILE* rsacp;
	SAC_HEADER* sachd;
	int ndumy;
	int npts;               // Points of data in SAC header 
	int nptsa;              // Points of data should be actually read
	unsigned int nptsb;     // Returned points of fread
	int j;
	char cbuf[4];
	
	float float_swap(char cbuf[]);
	 
	if (spt < 0) {
   	fprintf(stderr, "\nError in rsac2: Starting point (= %d) < 0 !! \n\n", spt);
   	return 0;
 	}		

	sachd = new SAC_HEADER;
	sachd->rsacheader(kname);
	npts = sachd->_npts; 

	if (spt >= npts) {
   	fprintf(stderr, "\nError in rsac2: Starting point (= %d) >= npts (= %d) !! \n\n", spt, npts);
   	return 0;
 	}		
	 
	rsacp = fopen(kname, "rb");
 	if (rsacp == NULL) {
   	fprintf(stderr, "\nError in rsac2: Cannot open sac data file %s!!\n\n",kname);
		fclose(rsacp);
   	return 0;
 	}
	
	if (fread(sachd->_hdmtr, sizeof(SAC_HEADER_MTR), 1, rsacp) != 1) {
   	fprintf(stderr, "\nError in rsac2: Reading sac header of file %s !! \n\n", kname);
   	fclose(rsacp);
   	return 0;
 	}		
	
/*----------------------------------------------------------------------------*/
/* Calculate the actrually read number of points: nptsa                       */
/*----------------------------------------------------------------------------*/
	if ((spt + nptsr) > npts) {
		nptsa = npts - spt;
	}	
	else {
		nptsa = nptsr;
	}
						
	if (fseek(rsacp, spt*sizeof(float), SEEK_CUR) != 0) {
		fprintf(stderr, "\nError in rsac2: Seek sac data in file %s!! \n\n", kname);
	}
	
	if (sachd->is_swap == 1) {
		for (j = 0; j < nptsa; j++) {
			ndumy = fread(cbuf, sizeof(char), 4, rsacp);
			data[j] = (float) float_swap(cbuf);
		}	
	}
	else {
		nptsb = fread(data, sizeof(float), nptsa, rsacp);
		if (nptsb != (unsigned int) nptsa) {
			fprintf(stderr, "\nError in rsac2 %s:\n", kname); 
			fprintf(stderr, "Start pts = %d, supposed read pts = %d, sac pts = %d\n", spt, nptsr, npts);
			fprintf(stderr, "Calculated actrually reading sac pts = %d !! \n\n", nptsa);
			fprintf(stderr, "Actrually reading sac pts = %d !! \n\n", nptsb);
			fclose(rsacp);
			return 0;
		}		
	}
		
	fclose(rsacp);
	delete sachd;
	sachd = NULL;
	return nptsa;
}	

/*----------------------------------------------------------------------------*/
/*	Read part of sac data [spt, spt+nptsr-1], including point "spt".
	"spt" is accounted from the 1st sample of sac file (zero-shift).
	"spt" must be >= 0                              
	If "spt + nptsr" > sac file real length, read actrual length
	Return the actrually read number of data
	Note the difference of between rsac1, rsac2, rsac4:
	rsac1: The points [spt, spt+nptsr-1] in sac file is put in data[0, nptsr-1]
			 If spt+nptsr > npts(sac length), return error.
			 The return value is 1(succesful) or 0(error).
			 Used for read two sac files and process them, which require 
			 the same length of them.
	rsac2: The points [spt, spt+nptsr-1] in sac file is put in data[0, nptsr-1]
			 If spt+nptsr > npts(sac length), nptsr = npts - spt.
			 The return value is the actrually read number(succesful) or 0(error).
			 Used for cut sac file, or only process it as a new sac file.
	rsac4: The points [spt, spt+nptsr-1] in sac file is put in 
			 data[spt, spt+nptsr-1].
			 The other segments in "data" is unchanged.
			 Before calling, make sure zero "data" and assign enough "data" size.
			 If spt+nptsr > npts(sac length), nptsr = npts - spt.
			 The return value is the actrually read number(succesful) or 0(error).
			 Used for cut sac file, but keep the same time sequence of 
			 original sac file, like in H-k search                               */
/*----------------------------------------------------------------------------*/
int rsac4(char *kname, float* data, const int& spt, const int& nptsr)
{
	FILE* rsacp;
	SAC_HEADER* sachd;
	int ndumy;
	int npts;               // Points of data in SAC header 
	int nptsa;              // Points of data should be actually read
	unsigned int nptsb;     // Returned points of fread
	int j;
	char cbuf[4];
	
	float float_swap(char cbuf[]);
	 
	if (spt < 0) {
   	fprintf(stderr, "\nError in rsac4: Starting point (= %d) < 0 !! \n\n", spt);
   	return 0;
 	}		

	sachd = new SAC_HEADER;
	sachd->rsacheader(kname);
	npts = sachd->_npts; 

	if (spt >= npts) {
   	fprintf(stderr, "\nError in rsac4: Starting point (= %d) >= npts (= %d) !! \n\n", spt, npts);
   	return 0;
 	}		
	 
	rsacp = fopen(kname, "rb");
 	if (rsacp == NULL) {
   	fprintf(stderr, "\nError in rsac4: Cannot open sac data file %s!!\n\n",kname);
		fclose(rsacp);
   	return 0;
 	}
	
	if (fread(sachd->_hdmtr, sizeof(SAC_HEADER_MTR), 1, rsacp) != 1) {
   	fprintf(stderr, "\nError in rsac4: Reading sac header of file %s !! \n\n", kname);
   	fclose(rsacp);
   	return 0;
 	}		
	
/*----------------------------------------------------------------------------*/
/* Calculate the actrually read number of points: nptsa                       */
/*----------------------------------------------------------------------------*/
	if ((spt + nptsr) > npts) {
		nptsa = npts - spt;
	}	
	else {
		nptsa = nptsr;
	}
						
	if (fseek(rsacp, spt*sizeof(float), SEEK_CUR) != 0) {
		fprintf(stderr, "\nError in rsac4: Seek sac data in file %s!! \n\n", kname);
	}
	
	if (sachd->is_swap == 1) {
		for (j = 0; j < nptsa; j++) {
			ndumy = fread(cbuf, sizeof(char), 4, rsacp);
			data[j+spt] = (float) float_swap(cbuf);
		}	
	}
	else {
		nptsb = fread(data+spt, sizeof(float), nptsa, rsacp);
		if (nptsb != (unsigned int) nptsa) {
			fprintf(stderr, "\nError in rsac4 %s:\n", kname); 
			fprintf(stderr, "Start pts = %d, supposed read pts = %d, sac pts = %d\n", spt, nptsr, npts);
			fprintf(stderr, "Calculated actrually reading sac pts = %d !! \n\n", nptsa);
			fprintf(stderr, "Actrually reading sac pts = %d !! \n\n", nptsb);
			fclose(rsacp);
			return 0;
		}		
	}
		
	fclose(rsacp);
	delete sachd;
	sachd = NULL;
	return nptsa;
}	

/*----------------------------------------------------------------------------*/
int wsac0(char *kname, float* data, SAC_HEADER* sachd)
{
	FILE* wsacp;

	wsacp = fopen(kname, "wb");
	if (wsacp == NULL) {
		fprintf(stderr, "\nError in wsac0: Cannot open file %s when write sac file!!\n\n", kname);
		fclose(wsacp);
		return(0);
	}

	// writing header
	if (fwrite(sachd->_hdmtr, sizeof(SAC_HEADER_MTR), 1, wsacp) != 1) {
		fprintf(stderr, "\nError in wsac0: Writing sac header for file %s!! \n\n", kname);
		fclose(wsacp);
		return(0);
	}
	  
	// writing data
	if (fwrite(data, sizeof(float), sachd->_npts, wsacp) != (unsigned)(sachd->_npts)) {
   	fprintf(stderr, "\nError in wsac0: Writing sac data %s!! \n\n", kname);
		fclose(wsacp);
   	return(0);
	}
	/*for (int i = 0; i < sachd->_npts; i++) {
		if (fwrite(&data[i], sizeof(float), 1, wsacp) != 1) {
   		fprintf(stderr, "\nError in wsac0: Writing %d-th sac data %s!! \n\n", i, kname);
			fclose(wsacp);
   		return(0);
		}
	}*/
	
	fclose(wsacp);
	return (1);	

}


/*----------------------------------------------------------------------------*/
int wsaclarge(char *kname, double* data, SAC_HEADER* sachd)
{
	FILE* wsacp;
	float fdata;

	wsacp = fopen(kname, "wb");
	if (wsacp == NULL) {
		fprintf(stderr, "\nError in wsaclarge: Cannot open file %s when write sac file!!\n\n", kname);
		fclose(wsacp);
		return(0);
	}

	// writing header
	if (fwrite(sachd->_hdmtr, sizeof(SAC_HEADER_MTR), 1, wsacp) != 1) {
		fprintf(stderr, "\nError in wsaclarge: Writing sac header for file %s!! \n\n", kname);
		fclose(wsacp);
		return(0);
	}
	  
	// writing data
	for (int i = 0; i < sachd->_npts; i++) {
		fdata = (float)data[i];
		if (fwrite(&fdata, sizeof(float), 1, wsacp) != 1) {
   		fprintf(stderr, "\nError in wsaclarge: Writing %d-th sac data %s!! \n\n", i, kname);
			fclose(wsacp);
   		return(0);
		}
	}
	
	fclose(wsacp);
	return (1);	

}


/*===========================================================================*/
/* Swap long int between big-endien and little-endien       Youlin Chen      */
/*===========================================================================*/
long long_swap(char cbuf[])
{
	union {
		char cval[4];
		long lval;
	} l_union;

	l_union.cval[3] = cbuf[0];
	l_union.cval[2] = cbuf[1];
	l_union.cval[1] = cbuf[2];
	l_union.cval[0] = cbuf[3];
	return(l_union.lval);
}

/*===========================================================================*/
/* Swap int between big-endien and little-endien            Youlin Chen      */
/*===========================================================================*/
int int_swap(char cbuf[])
{
	union {
		char cval[4];
		int lval;
	} l_union;

	l_union.cval[3] = cbuf[0];
	l_union.cval[2] = cbuf[1];
	l_union.cval[1] = cbuf[2];
	l_union.cval[0] = cbuf[3];
	return(l_union.lval);
}

/*===========================================================================*/
/* Swap float between big-endien and little-endien       Youlin Chen         */
/*===========================================================================*/
float float_swap(char cbuf[])
{
	union {
		char cval[4];
		float fval;
	} f_union;

	f_union.cval[3] = cbuf[0];
	f_union.cval[2] = cbuf[1];
	f_union.cval[1] = cbuf[2];
	f_union.cval[0] = cbuf[3];
	return(f_union.fval);
}

/*===========================================================================*/
/* Convert char to long int                              Youlin Chen         */
/*===========================================================================*/
long long_join(char cbuf[])
{
	union {
		char cval[4];
		long lval;
	} l_union;

	l_union.cval[3] = cbuf[3];
	l_union.cval[2] = cbuf[2];
	l_union.cval[1] = cbuf[1];
	l_union.cval[0] = cbuf[0];
	return(l_union.lval);
}

/*===========================================================================*/
/* Convert char to int                                   Youlin Chen         */
/*===========================================================================*/
int int_join(char cbuf[])
{
	union {
		char cval[4];
		int lval;
	} l_union;

	l_union.cval[3] = cbuf[3];
	l_union.cval[2] = cbuf[2];
	l_union.cval[1] = cbuf[1];
	l_union.cval[0] = cbuf[0];
	return(l_union.lval);
}


