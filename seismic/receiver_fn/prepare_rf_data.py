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


def get_events(lonlat, starttime, endtime, cat_file, early_exit=True):
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
        if early_exit:
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
@click.option('--inventory-file', type=click.File('r'), required=True,
              help='Path to input inventory file corresponding to waveform file, '
              'e.g. "/g/data/ha3/Passive/_ANU/7X\(2009-2011\)/ASDF/7X\(2009-2011\)_ASDF.xml".')
@click.option('--waveform-file', type=click.File('rb'), required=True,
              help='Path to input h5 waveform file from which to extract traces for RF analysis, '
              'e.g. "/g/data/ha3/Passive/_ANU/7X\(2009-2011\)/ASDF/7X\(2009-2011\).h5".')
@click.option('--event-catalog-file', type=click.Path(dir_okay=False, writable=True), required=True,
              help='Path to event catalog file, e.g. "catalog_7X_for_rf.xml". '
              'If file already exists, it will be loaded, otherwise it will be created.')
@click.option('--rf-trace-datafile', type=click.File('wb'), required=True,
              help='Path to output file, e.g. "7X_event_waveforms_for_rf.h5"')
@click.option('--start-time', type=str, required=True,
              help='Start datetime in ISO 8601 format, e.g. "2009-06-16T03:42:00"')
@click.option('--end-time', type=str, required=True,
              help='End datetime in ISO 8601 format, e.g. "2011-04-01T23:18:49"')
def main(inventory_file, waveform_file, event_catalog_file, rf_trace_datafile, start_time, end_time):

    assert not os.path.exists(rf_trace_datafile.name), \
        "Won't delete existing file {}, remove manually.".format(rf_trace_datafile.name)

    # we use centre of Australia to calculate radius and gather events from 15 to 90 degrees
    lonlat = [133.88, -23.69]

    inventory = read_inventory(inventory_file)
    waveform_datafile = waveform_file
    start_time = UTC(start_time)
    end_time = UTC(end_time)

    exit_after_catalog = False
    catalog = get_events(lonlat, start_time, end_time, event_catalog_file, exit_after_catalog)

    # Form closure to allow waveform source file to be derived from a setting (or command line input)
    def closure_get_waveforms(network, station, location, channel, starttime, endtime):
        return custom_get_waveforms(waveform_datafile, network, station, location, channel, starttime, endtime)

    stream = RFStream()
    with tqdm() as pbar:
        for s in iter_event_data(catalog, inventory, closure_get_waveforms, pbar=pbar):
            for trace in s:
                stream.extend(s)

    print(stream)
    stream.write(rf_trace_datafile, 'H5')


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
