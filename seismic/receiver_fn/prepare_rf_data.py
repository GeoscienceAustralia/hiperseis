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

import click
import pyasdf
from pyasdf import ASDFDataSet
from obspy.core import Stream


def get_events(lonlat, starttime, endtime, cat_file):
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


def custom_get_waveforms(waveform_datafile, network, station, location, channel, starttime,
                         endtime, quality=None, minimumlength=None,
                         longestonly=None, filename=None, attach_response=False,
                         **kwargs):
    with pyasdf.ASDFDataSet(waveform_datafile, mode='r') as asdfDataSet:
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

@click.command()
@click.option('--rf-trace-datafile', type=str, required=True, help='Path to output file, e.g. "7G_event_waveforms_for_rf.h5"')
def main(rf_trace_datafile):

    # rf_trace_datafile = os.path.join(os.path.split(__file__)[0], 'DATA', '7G_event_waveforms_for_rf.h5')
    assert not os.path.exists(rf_trace_datafile), "Won't delete existing file{}, remove manually.".format(rf_trace_datafile)

    # we use centre of Australia to calculate radius and gather events from 15 to 90 degrees
    lonlat = [133.88, -23.69]

    # Change parameters below
    invfile = "/g/data/ha3/Passive/SHARED_DATA/Inventory/networks/network_7G.xml"
    inventory = read_inventory(invfile)

    waveform_datafile = "/g/data/ha3/Passive/_ANU/7G(2013-2015)/ASDF/7G(2013-2015).h5"
    start_time = UTC('2014-01-01T00:00:06')
    end_time = UTC('2016-02-09T21:04:29')

    # ----------------- End ----------------------
    def fname_time(utc_datetime):
        return utc_datetime.strftime("%Y-%m-%dT%H%M%S")

    cat_file = os.path.join(os.path.split(__file__)[0], 'DATA', 'catalog_7G_' + fname_time(start_time) + '-' + fname_time(end_time) + '.xml')
    catalog = get_events(lonlat, start_time, end_time, cat_file)

    # Form closure to allow waveform source file to be derived from a setting (or command line input)
    def closure_get_waveforms(network, station, location, channel, starttime, endtime):
        return custom_get_waveforms(waveform_datafile, network, station, location, channel, starttime, endtime)

    stream = RFStream()
    with tqdm() as pbar:
        for s in iter_event_data(catalog, inventory, closure_get_waveforms, pbar=pbar):
            for trace in s:
                stream.extend(s)

    stream.write(rf_trace_datafile, 'H5')


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
