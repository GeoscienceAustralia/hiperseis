"""
Description:
    Implements a class that provides sunrise and sunset times, given lon, lat coordinates.
    Based on original implementation at:
    https://stackoverflow.com/questions/19615350/calculate-sunrise-and-sunset-times-for-a-given-gps-coordinate-within-postgresql

References:

CreationDate:   02/12/24
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     02/12/24   RH
"""

import math
from obspy.core import Stream, UTCDateTime

DAY_NIGHT_SECONDS = 86400

class Sun:
    def __init__(self, lon, lat):
        self.lon = lon
        self.lat = lat
    # end func

    def get_sunrise_time(self, day, month, year):
        return self.calc_sun_time(day, month, year, True)
    # end func

    def get_sunset_time(self, day, month, year):
        return self.calc_sun_time(day, month, year, False)
    # end func

    def calc_sun_time(self, day, month, year, isRiseTime, zenith=90.8):
        longitude = self.lon
        latitude = self.lat

        TO_RAD = math.pi / 180

        # 1. first calculate the day of the year
        N1 = math.floor(275 * month / 9)
        N2 = math.floor((month + 9) / 12)
        N3 = (1 + math.floor((year - 4 * math.floor(year / 4) + 2) / 3))
        N = N1 - (N2 * N3) + day - 30

        # 2. convert the longitude to hour value and calculate an approximate time
        lngHour = longitude / 15

        if isRiseTime:
            t = N + ((6 - lngHour) / 24)
        else:  # sunset
            t = N + ((18 - lngHour) / 24)

        # 3. calculate the Sun's mean anomaly
        M = (0.9856 * t) - 3.289

        # 4. calculate the Sun's true longitude
        L = M + (1.916 * math.sin(TO_RAD * M)) + (0.020 * math.sin(TO_RAD * 2 * M)) + 282.634
        L = self.force_range(L, 360)  # NOTE: L adjusted into the range [0,360)

        # 5a. calculate the Sun's right ascension

        RA = (1 / TO_RAD) * math.atan(0.91764 * math.tan(TO_RAD * L))
        RA = self.force_range(RA, 360)  # NOTE: RA adjusted into the range [0,360)

        # 5b. right ascension value needs to be in the same quadrant as L
        Lquadrant = (math.floor(L / 90)) * 90
        RAquadrant = (math.floor(RA / 90)) * 90
        RA = RA + (Lquadrant - RAquadrant)

        # 5c. right ascension value needs to be converted into hours
        RA = RA / 15

        # 6. calculate the Sun's declination
        sinDec = 0.39782 * math.sin(TO_RAD * L)
        cosDec = math.cos(math.asin(sinDec))

        # 7a. calculate the Sun's local hour angle
        cosH = (math.cos(TO_RAD * zenith) - (sinDec * math.sin(TO_RAD * latitude))) / (
                    cosDec * math.cos(TO_RAD * latitude))

        if cosH > 1:
            return None

        if cosH < -1:
            return None

        # 7b. finish calculating H and convert into hours

        if isRiseTime:
            H = 360 - (1 / TO_RAD) * math.acos(cosH)
        else:  # setting
            H = (1 / TO_RAD) * math.acos(cosH)

        H = H / 15

        # 8. calculate local mean time of rising/setting
        T = H + RA - (0.06571 * t) - 6.622

        # 9. adjust back to UTC
        UT = T - lngHour
        UT = self.force_range(UT, 24)  # UTC time in decimal format (e.g. 23.23)

        # 10. Return
        return UTCDateTime(year, month, day) + UT * 3600
    # end func

    def force_range(self, v, max):
        # force v to be >= 0 and < max
        if v < 0:
            return v + max
        elif v >= max:
            return v - max

        return v
    # end func
# end class

def day_night_coverage(day_st: Stream, lon, lat):
    """
    Computes overlap of input samples from a 24 hr stream with the extents of day and night
    on that date.
    @param day_st:
    @param lon:
    @param lat:
    @return: coverage of samples over daytime and nighttime as a tuple of fractions
    """
    def get_overlap(a, b):
        return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    # end func

    s = Sun(lon, lat)

    sunrise_ts = None
    sunset_ts = None
    date = None
    day = month = year = None
    day_coverage_seconds = night_coverage_seconds = 0
    day_length_seconds = night_length_seconds = 0
    for i, tr in enumerate(day_st):
        t_st = tr.stats.starttime
        t_et = tr.stats.endtime
        if (date is None):
            date = t_st.date
            day, month, year = date.day, date.month, date.year
            sunrise = s.get_sunrise_time(day, month, year)
            sunset = s.get_sunset_time(day, month, year)

            if (sunrise > sunset):  # night
                night_length_seconds = sunrise - sunset
                day_length_seconds = DAY_NIGHT_SECONDS - night_length_seconds
            else:  # day
                day_length_seconds = sunset - sunrise
                night_length_seconds = DAY_NIGHT_SECONDS - day_length_seconds
            # end if
        # end if

        assert date == t_st.date == t_et.date, \
            'Input stream crosses day boundary ({}, {}, {}, {}). Aborting..'. \
                format(day_st, date, t_st.date, t_et.date)

        t_len = t_et - t_st
        if (sunrise > sunset):  # night
            overlap = get_overlap([sunset, sunrise],
                                  [t_st, t_et])
            # print('night olap {}'.format(overlap))

            night_coverage_seconds += overlap
            day_coverage_seconds += t_len - overlap
        else:  # day
            overlap = get_overlap([sunrise, sunset],
                                  [t_st, t_et])
            # print('day olap {}'.format(overlap))

            day_coverage_seconds += overlap
            night_coverage_seconds += t_len - overlap
            # end if
    # end for

    return day_coverage_seconds / day_length_seconds, \
           night_coverage_seconds / night_length_seconds
# end func
