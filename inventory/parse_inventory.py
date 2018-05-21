from __future__ import absolute_import, print_function
import csv
import struct
from os.path import join, dirname
import logging
from collections import namedtuple
from obspy import read_inventory, UTCDateTime
import pandas as pd

log = logging.getLogger(__name__)

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

PASSIVE = dirname(dirname(__file__))

sc3_inventory = join(PASSIVE, 'inventory', 'stations.csv')
isc_file_1 = join(PASSIVE, 'inventory', 'ehb.stn')
isc_file_2 = join(PASSIVE, 'inventory', 'iscehb.stn')

Station = namedtuple('Station', 'station_code, latitude, longitude, '
                                'elevation, network_code')
def convert_inventory_to_df(inv):
    df = pd.DataFrame(columns=['station_code', 'latitude', 'longitude', 'elevation', 'network_code'])
    if inv and len(inv.networks) and len(inv.networks[0].stations):
        for i, st in enumerate(inv.networks[0].stations):
            df.loc[i] = [st.code, st.latitude, st.longitude, st.elevation, inv.networks[0].code]
    return df

def gather_isc_stations():
    log.info('Gathering isc and engdahl stations')
    sta_dict = {}
    # don't change file order, ehb.stn should be older stations list and we
    # keep the more modern ones in the dict
    for f in [isc_file_1, isc_file_2]:
        sta_dict = _read_sta_file(f, sta_dict)
    log.info('Gathered all stations')
    return sta_dict


def _read_sta_file(f, sta_dict):
    log.info('Gathering stations from {file}'.format(file=f))
    with open(f, 'r') as sta:
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
                sta_tuple = struct.unpack(isc_format, line.encode())

                # only select if station name is not empty string
                if sta_tuple[0].strip():
                    sta = sta_tuple[0].strip()
                    stations += 1
                    log.debug('found station {}, station_code: {}, '
                              'sta details: {}'.format(stations, sta_tuple[0],
                                                       sta_tuple))
                    try:
                        if len(sta_tuple) == 5:
                            lat, lon, ele, start_date = sta_tuple[1:5]
                            end_date = False
                        elif len(sta_tuple) == 6:
                            lat, lon, ele, start_date, end_date = \
                                sta_tuple[1:6]
                        else:  # no other formats should sneak in
                            raise Exception('We must have tuple of length'
                                            '5 or 6')
                        # check for conversion problems imply we need to use
                        # alternative formats
                        _ = float(lat.decode())
                        _ = float(lon.decode())
                        _ = float(ele.decode())
                        start_date = UTCDateTime(start_date.decode())
                        log.debug(
                            'Station: {sta} start date: {start_date}'.format(
                                sta=sta, start_date=start_date))
                        if end_date:
                            end_date = UTCDateTime(end_date)
                            log.debug(
                                'Station: {sta} end date: {end_date}'.format(
                                    sta=sta, end_date=end_date))
                    # exception handling due to 5 char station names
                    except ValueError:
                        log.debug(sta_tuple)
                        if line_length == 54:
                            sta_tuple = struct.unpack(five_char_isc_format2,
                                                      line.encode())
                            log.debug(sta_tuple)
                        elif line_length == 74:
                            sta_tuple = struct.unpack(five_char_isc_format1,
                                                      line.encode())
                            log.debug(sta_tuple)
                        else:
                            raise
                        sta = sta_tuple[0].strip()

                    if sta in sta_dict:
                        log.warning('Station {} already exists in the '
                                    'stations dict. '
                                    'Overwriting with most '
                                    'recent entry'.format(sta))
                    sta_dict[sta.decode()] = Station(
                        station_code=sta,
                        latitude=float(sta_tuple[1]),
                        longitude=float(sta_tuple[2]),
                        elevation=float(sta_tuple[3]),
                        network_code='dummy')

            else:
                # make sure we know all other formats before we ignore them
                assert line_length == 80 or line_length == 86 or \
                       line_length == 92
    return sta_dict


def read_all_stations():
    stations = read_stations(sc3_inventory)
    isc_stations = gather_isc_stations()
    stations.update(isc_stations)
    return stations


def read_stations(station_file):
    """
    Read station location from a csv file.
    :param station_file: str
        csv stations file handle passed in by click
    :return: stations_dict: dict
        dict of stations indexed by station_code for quick lookup
    """
    log.info('Reading seiscomp3 exported stations file')
    stations_dict = {}
    with open(station_file, 'r') as csv_file:
        reader = csv.reader(csv_file)
        next(reader)  # skip header
        for sta in reader:
            stations_dict[sta[0]] = Station(
                sta[0], float(sta[1]), float(sta[2]), float(sta[3]), sta[4]
            )
        log.info('Done reading seiscomp3 station files')
        return stations_dict


if __name__ == "__main__":
    #stations = gather_isc_stations()
    stations =  read_all_stations()

    print (len(stations.keys()))
#    inv = read_inventory('/home/ubuntu/7W.xml')
#    df = convert_inventory_to_df(inv)
#    df.to_csv('/home/ubuntu/7W.csv')
