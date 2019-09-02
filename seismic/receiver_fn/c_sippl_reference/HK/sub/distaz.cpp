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

#ifndef DR
#define DR 57.295780
#endif

#ifndef RE
#define RE 6371.227
#endif

#ifndef GEOCEN
#define GEOCEN 0.993305458
#endif

#ifndef PERD
#define PERD 111.18
#endif

/******************************************************************************/
int distaz(const float& epla, const float& eplo, const float& sla, const float& slo,
				float& del, float& dist, float& azi, float& baz)
/*----------------------------------------------------------------------------*/
/* A subroutine to calculate great circle distances and azimuth between 
	two points on the earth'h surface. 
	The earth is assummed to be an ellipsoid of revolution.
   This routine is from a program by lynn peseckis, 1979.

   input parameters :

   epla, eplo: Latitude and longitude of first point on earth's surface. 
               North latitude and east longitude is positive, 
               South latitude and west longitude is negative.
   sla, slo:   Latitude and longitude of second point on earth's surface.

   returned parameters :

   del ...............distance in degrees between two points.
   dist ..............distance in kilometers between two points.
   az ................azimuth from first point to second.
   baz ...............back azimuth from second point to first.                */
/*----------------------------------------------------------------------------*/
{
	double c;
	double bp;
	double ap;
	double gamma;
	double cbp;
	double sbp;
	double cap;
	double sap;
	double abtem;
	double altem;
	double baztem;
	double dlon;

/*----------------------------------------------------------------------------*/
/* begin calculation of distance                                              */
/*----------------------------------------------------------------------------*/
	if (fabs(epla) - 90.0 == 0.0) {
		bp = 90.0 - epla;
		ap = 90.0 - (atan(GEOCEN * tan(sla / DR))) * DR;
		if (epla <= 0.0) {
			azi = 0.0;
			baz = 180.0;
			del = fabs(ap - bp);
			dist = del * M_PI / 180.0 * RE;
			return 0;
		}
		else {
			azi = 180.0;
			baz = 0.0;
			del = fabs(ap - bp);
			dist = del * M_PI / 180.0 * RE;
			return 0;
		}
	}	 	
	else {
		bp = 90.0 - (atan(GEOCEN * tan(epla / DR))) * DR;
	}
	
	if (fabs(sla) - 90.0 == 0.0) {
		ap = 90.0 - sla;		 
		if (sla <= 0.0) {
			azi = 180.0;
			baz = 0.0;
			del = fabs(ap - bp);
			dist = del * M_PI / 180.0 * RE;
			return 0;
		}
		else {
			azi = 0.0;
			baz = 180.0;
			del = fabs(ap - bp);
			dist = del * M_PI / 180.0 * RE;
			return 0;
		}
	}
	else {
		ap = 90.0 - (atan(GEOCEN * tan(sla / DR))) * DR;
	}
	
	dlon = fabs(eplo-slo);
//	printf("here!!!!\n");
	if (dlon - 0.001 <= 0.0) {
		if (epla - sla <= 0.0) {
			azi = 0.0;
			baz = 180.0;
			del = fabs(ap - bp);
			dist = del * M_PI / 180.0 * RE;
			return 0;
		}
		else {
			azi = 180.0;
			baz = 0.0;
			del = fabs(ap - bp);
			dist = del * M_PI / 180.0 * RE;
			return 0;
		}
	}		 
	else {		 			 
		gamma = (slo - eplo) / DR;
		cbp = cos(bp / DR);
		sbp = sin(bp / DR);
		cap = cos(ap / DR);
		sap = sin(ap / DR);
		abtem = cbp * cap + sbp * sap * cos(gamma);
		c = acos(abtem);
		del = c * DR;
		altem = (cap - abtem * cbp) / (sin(c) * sbp);
		azi = acos(altem) * DR;
		baztem = (cbp - abtem * cap) / (sin(c) * sap);
		baz = acos(baztem) * DR;
		if (sin(gamma) < 0.0) {
			azi = 360.0 - azi;
		}		
		else {
			baz = 360.0 - baz;
		}
		dist = del * M_PI / 180.0 * RE;       
			return 0;
	}
	
}
