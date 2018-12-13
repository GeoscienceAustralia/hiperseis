#!/usr/bin/env python
'''Creates database of stations from .STN files which are not in IRIS database,
   curates the data using heuristic rules, and exports new stations to station.xml
'''
import os
import sys

import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', -1)
pd.set_option('display.width', 240)
import cStringIO as sio

try:
    import tqdm
    show_progress = True
except:
    show_progress = False

USE_PICKLE = True

from obspy import read_inventory
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from collections import defaultdict

station_columns = ['Latitude', 'Longitude', 'Elevation', 'StationStart', 'StationEnd', 'ChannelCode', 'ChannelStart', 'ChannelEnd']

def read_eng(fname):
    ''' We read Engdahl file in this function which parses the following format of fixed width formatted columns:

     AAI   Ambon             BMG, Indonesia, IA-Ne              -3.6870  128.1945      0.0   2005001  2286324  I
     AAII                                                       -3.6871  128.1940      0.0   2005001  2286324  I
     AAK   Ala Archa         Kyrgyzstan                         42.6390   74.4940      0.0   2005001  2286324  I
     ABJI                                                       -7.7957  114.2342      0.0   2005001  2286324  I
     APSI                                                       -0.9108  121.6487      0.0   2005001  2286324  I
     AS01  Alice Springs Arra                                  -23.6647  133.9508      0.0   2005001  2286324  I
    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
              10        20        30        40        50        60        70        80

    Each Station Code (first column) might NOT be unique, and network codes are missing here, so duplicates will
    simply be noted and handled later.
    '''
    colspec_0 = ((0, 6), (59, 67), (68, 77), (78, 86))
    col_names = ['StationCode', 'Latitude', 'Longitude', 'Elevation']
    data_frame = pd.read_fwf(fname, colspecs=colspec_0, names=col_names, index_col=0)
    # Populate missing data
    # data_frame[['Latitude', 'Longitude', 'Elevation']] = np.nan
    data_frame['StationStart'] = pd.NaT
    data_frame['StationEnd'] = pd.NaT
    data_frame['ChannelStart'] = pd.NaT
    data_frame['ChannelEnd'] = pd.NaT
    data_frame['ChannelCode'] = ''
    data_frame = data_frame[station_columns]
    # Assumed network code for files of this format.
    data_frame = pd.concat([data_frame], keys=['GE'], names=['NetworkCode'])
    # Compute number of duplicates
    num_dupes = len(data_frame) - len(data_frame.index.unique())
    print("{0}: {1} stations found with {2} duplicates".format(fname, len(data_frame), num_dupes))
    return data_frame

def read_isc(fname):
    # Engdahl file is imported and stored, now lets find appropriate network codes based on coordinates
    # first we read catalogue supplied by ISC that inherited Engdahl work
    '''Read ISC station inventory supplied by ISC, having such format:

    109C     32.8892 -117.1100     0.0 2006-06-01 04:11:18 2008-01-04 01:26:30
    109C     32.8882 -117.1050   150.0 2008-01-04 01:26:30
                 FDSN 109C   TA -- BHZ 2004-05-04 23:00:00 2005-03-03 23:59:59
                 FDSN 109C   TA -- LHZ 2004-05-04 23:00:00 2005-03-03 23:59:59
                 FDSN 109C   TA -- BHZ 2005-04-11 00:00:00 2006-01-25 22:31:10
                 FDSN 109C   TA -- LHZ 2005-04-11 00:00:00 2006-01-25 22:31:10
    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
              10        20        30        40        50        60        70        80

    The following policies are applied in the merging of duplicate data within a given file
    when read by this function:
    * The first level of the index is the Network Code.
    * For given station matched by station code and lat/long, a non-zero elevation, if present,
      overwrites a 0.0 elevation.
    * Where multiple Network Codes are present for a given station header, the network codes
      are preserved as the first level of the table index. Within each network code for a 
      given station, we tabulate the channel stages.
    * Where Station Code on a channel data row does not match header Station Code, the mismatch
      is logged, assumed to be an error, and ignored since the ISC data format in this file 
      implies the station code should be the same as the header row.

    Useful reference for adding level of multiindex: https://stackoverflow.com/questions/14744068/prepend-a-level-to-a-pandas-multiindex
    '''

    def reportStationCount(df):
        num_unique_networks = len(df.index.levels[0])
        num_unique_stations = len(df.index.levels[1])
        print("{0}: {1} unique network codes, {2} unique station codes".format(fname, num_unique_networks, num_unique_stations))

    if USE_PICKLE:
        pkl_name = fname + ".pkl"
        if os.path.exists(pkl_name): 
            print("Reading cached " + fname)
            with open(pkl_name, 'rb') as f:
                import cPickle as pkl
                df_all = pkl.load(f)
                reportStationCount(df_all)
                return df_all

    print("Parsing " + fname)
    header_colspec = ((0, 5), (7, 16), (17, 26), (27, 34), (35, 54,), (55, 74))
    header_cols = ['StationCode', 'Latitude', 'Longitude', 'Elevation', 'StationStart', 'StationEnd']
    channel_colspec = ((13, 17), (17, 22), (23, 27), (30, 34), (35, 54), (55, 74))
    channel_cols = ['FDSN', 'StationCode', 'NetworkCode', 'ChannelCode', 'ChannelStart', 'ChannelEnd']
    df_list = []
    # Due to irregular data format, read one row at a time from file. Consider replacing with a pre-processing pass 
    # so that more than one row can be read at a time.
    if show_progress:
        pbar = tqdm.tqdm(total=os.path.getsize(fname), ascii=True)
    with open(fname, "r", buffering=64*1024) as f:
        line = f.readline()
        if show_progress:
            pbar.update(len(line))
        hdr = None
        while line.strip():
            channels = []
            while 'FDSN' in line:
                # Read channel rows
                line_input = sio.StringIO(line)
                ch_data = pd.read_fwf(line_input, colspecs=channel_colspec, names=channel_cols, index_col=1, nrows=1, parse_dates=[4, 5])
                assert ch_data.iloc[0]['FDSN'] == 'FDSN'
                channels.append(ch_data)
                line = f.readline()
                if show_progress:
                    pbar.update(len(line))
            if hdr is not None:
                # Combine header and channel data. 
                if channels:
                    ch_all = pd.concat(channels, sort=False)
                    netcode = ch_all.iloc[0]['NetworkCode']
                    ch_all.drop('FDSN', axis=1, inplace=True)
                    # TODO: Add warning here if there is more than one network code. Only the first will be respected.
                    ch_all.drop('NetworkCode', axis=1, inplace=True)
                    network_df = hdr.join(ch_all)
                    # Make sure column ordering is consistent
                    network_df = network_df[station_columns]
                    # Add a first level index of network code
                    network_df = pd.concat([network_df], keys=[netcode], names=['NetworkCode'])
                else:
                    netcode = 'IR'
                    hdr['ChannelStart'] = pd.NaT
                    hdr['ChannelEnd'] = pd.NaT
                    hdr['ChannelCode'] = 'SHZ'
                    hdr = hdr[station_columns]
                    network_df = pd.concat([hdr], keys=[netcode], names=['NetworkCode'])
                df_list.append(network_df)
                hdr = None
            # Read header row
            line_input = sio.StringIO(line)
            hdr = pd.read_fwf(line_input, colspecs=header_colspec, names=header_cols, index_col=0, nrows=1, parse_dates=[4, 5])
            line = f.readline()
            if show_progress:
                pbar.update(len(line))
    if show_progress:
        pbar.close()
    print("Concatenating records...")
    df_all = pd.concat(df_list, sort=False)
    if USE_PICKLE:
        with open(fname + ".pkl", "wb") as f:
            import cPickle as pkl
            pkl.dump(df_all, f, pkl.HIGHEST_PROTOCOL)
    reportStationCount(df_all)
    return df_all

def main(argv):
    # Read IRIS station database. Stations from STN files which exist here will 
    IRIS_all_file = "IRIS-ALL.xml"
    print("Reading " + IRIS_all_file)
    if USE_PICKLE:
        IRIS_all_pkl_file = IRIS_all_file + ".pkl"
        if os.path.exists(IRIS_all_pkl_file): # doesn't work reliably on EC2 CentOS instance for some reason
            with open(IRIS_all_pkl_file, 'rb') as f:
                import cPickle as pkl
                inv = pkl.load(f)
        else:
            with open(IRIS_all_file, 'r', buffering=1024*1024) as f:
                inv = read_inventory(f)
            with open(IRIS_all_pkl_file, 'wb') as f:
                import cPickle as pkl
                pkl.dump(inv, f, pkl.HIGHEST_PROTOCOL)
    else:
        with open(IRIS_all_file, 'r', buffering=1024*1024) as f:
            inv = read_inventory(f)

    # Read station database from ad-hoc formats
    ehb_data_bmg = read_eng('BMG.STN')
    ehb_data_isc = read_eng('ISC.STN')
    ehb_eng = pd.concat([ehb_data_bmg, ehb_data_isc], sort=False)
    ehb_eng.sort_values(['NetworkCode', 'StationCode'], inplace=True)

    isc1 = read_isc('ehb.stn')
    isc2 = read_isc('iscehb.stn')
    # isc1 = read_isc(os.path.join('test', 'ehb_test.stn'))
    # isc2 = read_isc(os.path.join('test', 'iscehb_test.stn'))
    ehb_isc = pd.concat([isc1, isc2], sort=False)
    ehb_isc.sort_values(['NetworkCode', 'StationCode'], inplace=True)

    pass

if __name__ == "__main__":
    main(sys.argv)
