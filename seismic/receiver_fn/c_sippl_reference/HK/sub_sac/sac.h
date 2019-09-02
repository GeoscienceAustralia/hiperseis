#ifndef SAC_H
#define SAC_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cstdio> 

#include <cmath>

/*--- define enumerated header value -----------*/
#define ITIME  1  
#define IRLIM  2  
#define IAMPH  3  
#define IXY  4  
#define IUNKN  5  
#define IDISP  6  
#define IVEL  7  
#define IACC  8  
#define IB  9  
#define IDAY  10  
#define IO  11  
#define IA  12  
#define IT0  13  
#define IT1  14  
#define IT2  15  
#define IT3  16  
#define IT4  17  
#define IT5  18  
#define IT6  19  
#define IT7  20  
#define IT8  21  
#define IT9  22  
#define IRADNV  23  
#define ITANNV  24  
#define IRADEV  25  
#define ITANEV  26  
#define INORTH  27  
#define IEAST  28  
#define IHORZA  29  
#define IDOWN  30  
#define IUP  31  
#define ILLLBB  32  
#define IWWSN1  33  
#define IWWSN2  34  
#define IHGLP  35  
#define ISRO  36  
#define INUCL  37  
#define IPREN  38  
#define IPOSTN  39  
#define IQUAKE  40  
#define IPREQ  41  
#define IPOSTQ  42  
#define ICHEM  43  
#define IOTHER  44  
#define IGOOD  45  
#define IGLCH  46  
#define IDROP  47  
#define ILOWSN  48  
#define IRLDTA  49  
#define IVOLTS  50  
#define IXYZ  51  
#define IMB  52  
#define IMS  53  
#define IML  54  
#define IMW  55  
#define IMD  56  
#define IMX  57  
#define INEIC  58  
#define IPDE  59  
#define IISC  60  
#define IREB  61  
#define IUSGS  62  
#define IBRK  63  
#define ICALTECH  64  
#define ILLNL  65  
#define IEVLOC  66  
#define IJSOP  67  
#define IUSER  68  
#define IUNKNOWN  69  
#define IQB  70  
#define IQB1  71  
#define IQB2  72  
#define IQBX  73  
#define IQMT  74  
#define IEQ  75  
#define IEQ1  76  
#define IEQ2  77  
#define IME  78  
#define IEX  79  
#define INU  80  
#define INC  81  
#define IO_  82  
#define IL  83  
#define IR  84  
#define IT  85  
#define IU  86  

#define TRUE 1
#define FALSE 0

#define FLOAT_ROWS 14
#define FLOAT_COLS  5
#define LONG_ROWS   8
#define LONG_COLS   5
#define CHAR_ROWS   8
#define CHAR_COLS   3
#define CHAR_BYTES  8

typedef struct { 
  float headf[FLOAT_ROWS][FLOAT_COLS];
  int headi[LONG_ROWS][LONG_COLS]; 
  char headc[CHAR_ROWS][CHAR_COLS * CHAR_BYTES];
} 
SAC_HEADER_MTR;

class SAC_HEADER { 
  public:
		SAC_HEADER();
		~SAC_HEADER();
		void set_b(float b);
		void set_e(float e);
		void set_o(float o);
		void set_delta(float delta);
		void set_npts(int npts);
		void set_scale(float scale);
		void set_reftime(int nzyear, int nzjday, int nzhour, int nzmin, int nzsec, int nzmsec);
		void set_sta(char *stnm, float stla, float stlo, float stel, float stdp); 
		void set_stdp(float stdp);
		void set_evt(float evla, float evlo, float evel, float evdp, float mag, int imagtyp);
		void set_mag(float);
		void set_imagtyp(int);
		void set_cmp(float cmpaz, float cmpinc);
		void set_css(int nevid, int norid, int nwfid);
		void set_iftype(int);
		void set_idep(int);
		void set_iztype(int);
		void set_ievtyp(int);
		void set_isynth(int);
		void set_leven(int);
		void set_lpspol(int);
		void set_lovrok(int);
		void set_lcalda(int);
		void set_unusedf(int, int, float); /*index count from 1 */
		void set_unusedi(int, int, int);
		void set_userf(int, float);   /*int for 0-9 */
		void set_knetwk(char* net);
		void set_kstnm(char* stacode);
		void set_kcmpnm(char* comp);
		void set_kdatrd(char* kdatrd);
		void set_khole(char*);
		void set_userk(int, char*);  /*int for 0-2 */
		void cpsachdr(SAC_HEADER* sachdr);
		
		int _nzyear;
		int _nzjday;
		int _nzhour;
		int _nzmin;
		int _nzsec;
		int _nzmsec;
		int _nvhdr;
		int _nevid;
		int _norid;
		int _nwfid;
		int _iftype;
		int _idep;
		int _iztype;
		int _ievtyp;
		int _isynth;
		int _leven;
		int _lovrok;
		int _lpspol;
		int _lcalda;
		int _imagtyp;
		int _npts;
		
		float _b;
		float _e;
		float _o;
		float _delta;
		float _scale;
		float _stla;
		float _stlo;
		float _stel;
		float _stdp;
		float _cmpaz;
		float _cmpinc;
		float _evla;
		float _evlo;
		float _evel;
		float _evdp;
		float _mag;
		
		float _user0;
		float _user1;
		float _user2;
		float _user3;
		float _user4;
		float _user5;
		float _user6;
		float _user7;
		float _user8;
		float _user9;
		float _uusd1304;
		float _uusd1305;

		char* _kstnm;
		char* _khole;
		char* _knetwk;
		char* _kcmpnm;
		char* _kdatrd;
		char* _kuser0;
		char* _kuser1;
		char* _kuser2;
		
		void rsacheader(char *kname);
		void wsacheader(char *kname);
		friend int rsac0(char *kname, float* data, int& npts, float& delta, float& b);
		friend int rsac1(char *kname, float* data, const int& spt, const int& nptsr);
		friend int rsac2(char *kname, float* data, const int& spt, const int& nptsr);
		friend int rsac4(char *kname, float* data, const int& spt, const int& nptsr);
		friend int wsac0(char *kname, float* data, SAC_HEADER* sachd);
		friend int wsaclarge(char *kname, double* data, SAC_HEADER* sachd);

	private:
		SAC_HEADER_MTR* _hdmtr;
		int is_swap;	
}; 

#endif
