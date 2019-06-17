#include "util_time.h"

/******************************************************************************/
TIME::TIME()
{
	_year = -12345;
	_mon = -12345;
	_day = -12345;
	_jday = -12345;
	_hour = -12345; 
	_min = -12345; 
	_sec = -12345; 
	_msec = -12345;
	_secm = -12345.0; 
}


/******************************************************************************/
void TIME::syear(int year)
{
	_year = year;
}
	
void TIME::smon(int mon)
{
	_mon = mon;
}
	
void TIME::sday(int day)
{
	_day = day;
}
	
void TIME::sjday(int jday)
{
	_jday = jday;
}
	
void TIME::shour(int hour)
{
	_hour = hour;
}
	
void TIME::smin(int min)
{
	_min = min;
}
	
void TIME::ssec(int sec)
{
	_sec = sec;
}
	
void TIME::smsec(int msec)
{
	_msec = msec;
}
	
void TIME::ssecm(float secm)
{
	_secm = secm;		
}


/******************************************************************************/
void TIME::setT_sm(int yr, int mo, int da, int hr, int mi, float secm)
{
	_year = yr;
	_mon = mo;
	_day = da;
	mday2jday();
	_hour = hr;
	_min = mi;
	_secm = secm;
	_sec = (int)secm;
	_msec = (int)((secm - (int)secm) * 1000);
}	

	
void TIME::setT_ms(int yr, int mo, int da, int hr, int mi, int sec, int msec)
{
	_year = yr;
	_mon = mo;
	_day = da;
	mday2jday();
	_hour = hr;
	_min = mi;
	_sec = sec;
	_msec = msec;
	_secm = (float)sec + (float)msec/1000.0;
}	

	
void TIME::setT_jdms(int yr, int jday, int hr, int mi, int sec, int msec)
{
	_year = yr;
	_jday = jday;
	jday2mday();
	_hour = hr;
	_min = mi;
	_sec = sec;
	_msec = msec;
	_secm = (float)sec + (float)msec/1000.0;
}	


void TIME::cpTime(TIME& t)
{
	_year = t.gyear();
	_mon = t.gmon();
	_day = t.gday();
	_hour = t.ghour();
	_min = t.gmin();
	_secm = t.gsecm();
	_sec = t.gsec();
	_msec = t.gmsec();
}

	
/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* From jday to mon and day:
	require jday != 0, refill mon and day field in TIME.                       */
/*----------------------------------------------------------------------------*/
int TIME::jday2mday()
{
	static int days[13] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
	static int days_leap[13] = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}; 
	int i;

	if (_year > 0 && _jday > 0) {
		if ( (_year % 4 == 0 && _year % 100 != 0) || (_year % 400 == 0) ) {
			if (_jday > 366) {
				printf("\nError julian day (%i) in year %i\n\n", _jday, _year);
				return -1;
			}
			else {
				for (i = 0; i < MONTHS_ONE_YEAR; i++) {
					if (_jday > days_leap[i] && _jday <= days_leap[i+1]) {
						_mon = i+1;
						_day = _jday - days_leap[i];
						break;
					}
				}
			}
		}
		else {
			if (_jday > 365) {
				printf("\nError julian day (%i) in year %i\n\n", _jday, _year);
				return -1;
			}
			else {
				for (i = 0; i < MONTHS_ONE_YEAR; i++) {
					if (_jday > days[i] && _jday <= days[i+1]) {
						_mon = i+1;
						_day = _jday - days[i];
						break;
					}
				}
			}
		}
	}
	else {
		printf("\nInvalid julian day (%i) in year %i\n\n", _jday, _year);
		return -1;
	}
	
	return 0;
}


/*----------------------------------------------------------------------------*/
/* From mon day to jday:
	require mon, day != 0, refill jday in TIME.                                */
/*----------------------------------------------------------------------------*/
int TIME::mday2jday()
{
	static int days[13] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
	static int days_leap[13] = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}; 

	if (_year > 0 && _mon > 0 && _day > 0) {
		if ( (_year % 4 == 0 && _year % 100 != 0) || (_year % 400 == 0) ) {
			_jday = days_leap[_mon-1] + _day;
		}
		else {
			_jday = days[_mon-1] + _day;
		}
	}
	else {
		printf("\nInvalid julian day (%i) in year %i\n\n", _jday, _year);
		return -1;
	}
	
	return 0;
}	
	
	

/******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Time increment and decrement:
	Calculates the t2 = t1 + diff(sec)
	Required fields in t1 are year, month, jday, hour, min and secm.
	If all fields are set in t1, their consistency is not checked.
	If jday is not set, calculate from year, month and day. 
	Return all fields in t2.  
	Leap years are accounted for.
	If diff > 0, t2 is a time after t1.
	If diff < 0, t2 is a time before t1.                                       */
/*----------------------------------------------------------------------------*/
int time_increment(TIME& t1, double diff, TIME& t2)
{
	static int days[13] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
	static int days_leap[13] = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}; 
	double t1sec;
	double t2sec;
	double secpart;
	int i;
	int more_loop;
	int in_mins;
	int in_hours;
	int in_days;
	int yr;

/*----------------------------------------------------------------------------*/
/* Process from mon and day to jday                                           */
/*----------------------------------------------------------------------------*/
	if (t1._jday == -12345) {
		t1.mday2jday();
	}

/*----------------------------------------------------------------------------*/
/* Calculate 2nd time after diff sec                                          */	
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
/* t1sec is number of seconds from t1.year/1/1 00:00:00.000 to t1             */
/*----------------------------------------------------------------------------*/
	if (t1._secm == -12345.0) {
		t1._secm = t1._sec + t1._msec/1000.0;
	}
	t1sec = ((t1._jday * 24.0 + t1._hour) * 60.0 + t1._min) * 60.0 + t1._secm; 
			  
	yr = t1._year;
  
/*----------------------------------------------------------------------------*/
/* t2sec is number of seconds from t1.yr/1/1 00:00:00.000 to t2               */
/*----------------------------------------------------------------------------*/
	do {					    
		t2sec = t1sec + diff;
		if (t2sec <= 0.0) {
/*----------------------------------------------------------------------------*/
/* rebuild t1sec and yr                                                       */
/*----------------------------------------------------------------------------*/
			yr = yr - 1;
			if ( (yr % 4 == 0 && yr % 100 != 0) || (yr % 400 == 0) ) {
				t1sec = t1sec + 366*24*60.0*60.0;
			}
			else {
				t1sec = t1sec + 365*24*60.0*60.0;
			}
			more_loop = 1;		 
		}
		else {
			more_loop = 0;
		}
	}	while(more_loop);

	in_mins = (int) (t2sec/60.0);
	secpart = fmod(t2sec, 60.0);
	t2._secm = (float) secpart;
	t2._sec = (int) secpart;
	t2._msec = (int) ((secpart - (int) secpart) * 1000);
  
	in_hours = in_mins / 60;
	t2._min = in_mins % 60;

	in_days = in_hours / 24;
	t2._hour = in_hours % 24;
  
	more_loop = 1;
	while(more_loop) {
		if ( (yr % 4 == 0 && yr % 100 != 0) || (yr % 400 == 0) ) {
			if (in_days > 366) {
				in_days = in_days - 366;
				yr++;
			}
			else {
				more_loop = 0;
				t2._jday = in_days;
				t2._year = yr;
				for (i = 0; i < MONTHS_ONE_YEAR; i++) {
					if (t2._jday > days_leap[i] && t2._jday <= days_leap[i+1]) {
						t2._mon = i+1;
						t2._day = t2._jday - days_leap[i];
						break;
					}
				}
			}
		}
		else {
			if (in_days > 365) {
				in_days = in_days - 365;
				yr++;
			}
			else {
				more_loop = 0;
				t2._jday = in_days;
				t2._year = yr;
				for (i = 0; i < MONTHS_ONE_YEAR; i++) {
					if (t2._jday > days[i] && t2._jday <= days[i+1]) {
						t2._mon = i+1;
						t2._day = t2._jday - days[i];
						break;
					}
				}
			}
		}
	} // while

	return 0;
}

/*----------------------------------------------------------------------------*/
/* Time difference 
	Calculates the time difference (t2 - t1)
	Required fields of t2 and t1 are year, mon, day, hour, min and secm.  
	If all fields are set in t1 and t2, their consistency is not checked.
	If day is not set, calculate from year, jday. 
	Leap years are accounted for.
	If t2 < t1, return negative time difference
	The time difference, in second, is returned in the double precision.       */
/*----------------------------------------------------------------------------*/
double time_dif(TIME& t2, TIME& t1)
{
	static int days[13] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
	static int days_leap[13] = {0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366}; 
	double t1sec;
	double t2sec;
	int min_yr;
	int max_yr;
	double diff;
	double yr_diff;
	int yr;
	
	if (t1._day == -12345 && t1._mon == -12345) t1.jday2mday();
	if (t2._day == -12345 && t2._mon == -12345) t2.jday2mday();
	if (t1._secm == -12345.0) t1._secm = t1._sec + t1._msec/1000.0;
	if (t2._secm == -12345.0) t2._secm = t2._sec + t2._msec/1000.0;
	
/*----------------------------------------------------------------------------*/
/* t1sec is number of seconds from t1.year/1/1 00:00:00.000 to t1             */
/*----------------------------------------------------------------------------*/
	if ( (t1._year % 4 == 0 && t1._year % 100 != 0) || (t1._year % 400 == 0) ) {
		t1sec = (( (days_leap[t1._mon-1]+t1._day-1) * 24.0 + t1._hour) * 60.0 + t1._min) * 60.0 
					+ t1._secm;  
	}
	else {
		t1sec = (( (days[t1._mon-1]+t1._day-1) * 24.0 + t1._hour) * 60.0 + t1._min) * 60.0 
					+ t1._secm;  
	}
	
/*----------------------------------------------------------------------------*/
/* t2sec is number of seconds from t2.year/1/1 00:00:00.000 to t2             */  
/*----------------------------------------------------------------------------*/
	if ( (t2._year % 4 == 0 && t2._year % 100 != 0) || (t2._year % 400 == 0) ) {
		t2sec = (( (days_leap[t2._mon-1]+t2._day-1) * 24.0 + t2._hour) * 60.0 + t2._min) * 60.0 
					+ t2._secm;  
	}
	else {
		t2sec = (( (days[t2._mon-1]+t2._day-1) * 24.0 + t2._hour) * 60.0 + t2._min) * 60.0 
					+ t2._secm;  
	}

/*----------------------------------------------------------------------------*/
/* Calculate the time in seconds from t1.year/1/1 00:00:00.000 to t2.year/1/1 00:00:00.000
	Consider the leap years
	Consider the case that t2 < t1                                             */ 
/*----------------------------------------------------------------------------*/
	if (t1._year < t2._year) {
		min_yr = t1._year;
		max_yr = t2._year;
	}
	else {
		min_yr = t2._year;
		max_yr = t1._year;
	}

	yr_diff = 0.0;
	for (yr = min_yr; yr < max_yr; yr++) {
		if ( (yr % 4 == 0 && yr % 100 != 0) || (yr % 400 == 0) ) {
			yr_diff = yr_diff + 366 * 24.0 * 60.0 * 60.0;
		}
		else {
			yr_diff = yr_diff + 365 * 24.0 * 60.0 * 60.0;
		}
	}

	if (t1._year < t2._year) {
		diff = yr_diff + t2sec - t1sec;
	}
	else if (t1._year > t2._year) {
		diff = -(yr_diff + t1sec - t2sec);
	}	
	else {
		diff = t2sec - t1sec;
	}

	return diff;

}


