#ifndef _UTIL_TIME_H_
#define _UTIL_TIME_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <cstring>
#include <cstdlib>
#include <cctype>
#include <cstdio> 

#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <time.h>

#define MONTHS_ONE_YEAR 12

class TIME {
	public:
		TIME();
		void setT_sm(int yr, int mo, int da, int hr, int mi, float secm);
		void setT_ms(int yr, int mo, int da, int hr, int mi, int sec, int msec);
		void setT_jdms(int yr, int jday, int hr, int mi, int sec, int msec);
		void cpTime(TIME& t);
		int jday2mday();
		int mday2jday();
		void syear(int);
		void smon(int);
		void sday(int);
		void sjday(int);
		void shour(int);
		void smin(int);
		void ssec(int);
		void smsec(int);
		void ssecm(float);		
		int gyear() const {return _year; }
		int gmon() const {return _mon; }
		int gday()  const {return _day; }
		int gjday()  const {return _jday; }
		int ghour()  const {return _hour; }
		int gmin()  const {return _min; }
		int gsec() const {return _sec; }
		int gmsec()  const {return _msec; }
		float gsecm()  const {return _secm; }
		
		friend int time_increment(TIME& t1, double diff, TIME& t2);
		friend double time_dif(TIME& t2, TIME& t1);

	private:	
		int _year;
		int _mon;
		int _day;
		int _jday;
		int _hour; 
		int _min; 
		int _sec; 
		int _msec;
		float _secm; 
		
};

#endif
