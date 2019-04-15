#!/usr/bin/env python
import os.path
import matplotlib.pyplot as plt
import numpy as np
from obspy import read_inventory, read_events, UTCDateTime as UTC
from obspy.clients.fdsn import Client
from rf import RFStream
from rf import iter_event_data
# from rf.imaging import plot_profile_map
# from rf.profile import profile
from tqdm import tqdm

import pyasdf
from pyasdf import ASDFDataSet
from obspy.core import Stream


def get_events(lonlat, starttime, endtime):
    cat_file = 'DATA/catalog' + str(starttime).replace(' ', '-') + '-' + str(endtime).replace(' ', '-') + '.xml'
    if os.path.exists(cat_file):
        catalog = read_events(cat_file)
    else:
        client = Client('ISC')
        kwargs = {'starttime': starttime, 'endtime': endtime, 
                  'latitude': lonlat[1], 'longitude': lonlat[0],
                  # we can use distances from 15 degrees, see Levin et al.
                  'minradius': 15, 'maxradius': 90,
                  'minmagnitude': 5.5, 'maxmagnitude': 6.5}
        print("Following parameters for earthquake extraction will be used:")
        print('starttime', starttime, 'endtime', endtime, 'latitude', lonlat[1], 'longitude', 
              lonlat[0], 'minradius: 15', 'maxradius: 90', 'minmagnitude: 5.5', 'maxmagnitude : 6.5')
        catalog = client.get_events(**kwargs)
        catalog.write(cat_file, 'QUAKEML')
        print("Catalog loaded")
        print("Run this process again using qsub")
        exit(0)
    return catalog


def custom_get_waveforms(network, station, location, channel, starttime,
                         endtime, quality=None, minimumlength=None,
                         longestonly=None, filename=None, attach_response=False,
                         **kwargs):
    with pyasdf.ASDFDataSet('/g/data/ha3/Passive/_ANU/7X(2009-2011)/ASDF/7X(2009-2011).h5', mode='r') as asdfDataSet:
        st = Stream()
        # ignoring channel for now as all the 7D network waveforms have only BH? channels
        filteredList = [i for i in asdfDataSet.waveforms[network + '.' + station].list() if
                        'raw_recording' in i and
                        UTC(i.split("__")[1]) < starttime and
                        UTC(i.split("__")[2]) > endtime]
        for t in filteredList:
            st += asdfDataSet.waveforms[network + '.' + station][t]
        return st


# ---+----------Main---------------------------------

if __name__ == '__main__':

    # we use centre of Australia to calculate radius and gather events from 15 to 90 degrees
    lonlat = [133.88, -23.69]

    # Change parameters below
    data = os.path.join('DATA', '')
    invfile = data + '7X-inventory.xml'
    datafile = data + '7X-event_waveforms_for_rf.h5'

    start_time = '2009-12-01 00:00:00'
    end_time = '2011-04-01 00:00:00'
    inventory = read_inventory(invfile)

    # ----------------- End ----------------------

    catalog = get_events(lonlat, UTC(start_time), UTC(end_time))

    stream = RFStream()
    with tqdm() as pbar:
        for s in iter_event_data(catalog, inventory, custom_get_waveforms, pbar=pbar):
            for trace in s:
                stream.extend(s)

    stream.write(datafile, 'H5')
