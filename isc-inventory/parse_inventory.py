from __future__ import absolute_import
import struct
import os
from obspy import UTCDateTime
isc_format1 = '4s12s10s8s20s20s'
isc_format2 = '4s12s10s8s20s'
isc_format3 = '5s12s10s8s20s20s'
isc_format4 = '5s12s10s8s20s'
five_char_isc_format1 = '5s11s10s8s20s20s'
five_char_isc_format2 = '5s11s10s8s20s'

"""
The first section is the data from the (ISC-)EHB station list showning station 
code, lat, lon, elev (m) and then the first and last dates the coordinates 
were used. No last date implies ongoing use.

The second section is our map from the code to FDSN where available, showing 
sta, network, location, channel and open date, close date and if we have it 
the sampling rates.
"""

isc_file_1 = os.path.join('isc-inventory', 'ehb.stn')
isc_file_2 = os.path.join('isc-inventory', 'iscehb.stn')

with open(isc_file_2, 'r') as sta:
    stations = 0
    for l, line in enumerate(sta):
        line = line.rstrip()
        line_length = len(line)
        if line_length in [74, 54, 75, 55]:
            if line_length == 74:
                isc_format = isc_format1
            elif line_length == 54:
                isc_format = isc_format2
            elif line_length == 75:
                isc_format = isc_format3
            else:  # line_length == 55
                isc_format = isc_format4
            sta_tuple = struct.unpack(isc_format, line)

            # exception handling due to 5 char station names
            if sta_tuple[0].strip():
                stations += 1
                print('found station {}'.format(stations),
                      'station_code: ', sta_tuple[0],
                      'sta details: ', sta_tuple)
                try:
                    if len(sta_tuple) == 5:
                        lat, lon, ele, start_date = sta_tuple[1:5]
                        end_date = False
                    elif len(sta_tuple) == 6:
                        lat, lon, ele, start_date, end_date = sta_tuple[1:6]

                    _ = float(lat)
                    _ = float(lon)
                    _ = float(ele)
                    start_date = UTCDateTime(start_date)
                    print('This station start date: ', start_date)
                    if end_date:
                        end_date = UTCDateTime(end_date)
                        print('This station end date: ', end_date)
                except ValueError:
                    print(sta_tuple)
                    if line_length == 54:
                        sta_tuple = struct.unpack(five_char_isc_format2, line)
                        print(sta_tuple)
                        # raise
                    elif line_length == 74:
                        sta_tuple = struct.unpack(five_char_isc_format1, line)
                        print(sta_tuple)
                    else:
                        raise
        else:
            # make sure we know all other formats before we ignore them
            assert len(line) == 80 or len(line) == 86 or len(line) == 92

