#!/usr/bin/env python
#
# Simple helper script to generate small data files for validation of Python setup on Raijin

import sys

from obspy import UTCDateTime, Stream
from obspy.clients.fdsn.client import Client

import pyasdf

NETWORK = 'AU'
CHANNEL = 'BHZ'
TAG = 'raw_recording'
TIME_RANGE = ("2011-03-11T00:00:00Z", "2011-03-12T01:00:00Z")

def generateStationTestData(sta):

    time_range = (UTCDateTime(TIME_RANGE[0]),
                  UTCDateTime(TIME_RANGE[1]))

    client = Client("IRIS")
    inv = client.get_stations(network=NETWORK, station=sta, channel=CHANNEL, 
                              starttime=time_range[0], endtime=time_range[1], level='channel')
    print(inv)

    traces = client.get_waveforms(network=NETWORK, station=sta, channel=CHANNEL, location='*', 
                                  starttime=time_range[0], endtime=time_range[1])
    print(traces)

    outfile = 'test_data_' + sta + '.h5'
    asdf_out = pyasdf.ASDFDataSet(outfile, mode='w')
    asdf_out.add_stationxml(inv)
    asdf_out.add_waveforms(traces, TAG)

    print("Saved data to " + outfile)


if __name__ == "__main__":
    generateStationTestData('ARMA')
    generateStationTestData('CMSA')
    generateStationTestData('QLP')

