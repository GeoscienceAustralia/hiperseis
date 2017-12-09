import os.path

import matplotlib.pyplot as plt
import numpy as np
from obspy import read_inventory, read_events, UTCDateTime as UTC
from obspy.clients.fdsn import Client
from rf import read_rf, RFStream
from rf import get_profile_boxes, iter_event_data, IterMultipleComponents
from rf.imaging import plot_profile_map
from rf.profile import profile
from tqdm import tqdm

import pyasdf
from pyasdf import ASDFDataSet
from obspy.core import Stream

def custom_get_waveforms(network, station, location, channel, starttime,
                  endtime, quality=None, minimumlength=None,
                  longestonly=None, filename=None, attach_response=False,
                  **kwargs):
    with pyasdf.ASDFDataSet('/g/data/ha3/Passive/_ANU/7D(2012-2013)/ASDF/7D(2012-2013).h5', mode='r') as asdfDataSet:
        st = Stream()
        #ignoring channel for now as all the 7D network waveforms have only BH? channels
        filteredList = [i for i in asdfDataSet.waveforms[network+'.'+station].list() if
                        'raw_recording' in i and
                        UTC(i.split("__")[1]) < starttime and
                        UTC(i.split("__")[2]) > endtime]
        for t in filteredList:
            st += asdfDataSet.waveforms[network+'.'+station][t]
        return st

data = os.path.join('data', '')
invfile = data + '7D-inventory.xml'
catfile = data + '7D-catalog.xml'
datafile = data + '7D-rf_profile_data.h5'
rffile = data + '7D-rf_profile_rfs.h5'
profilefile = data + '7D-rf_profile_profile.h5'

inventory = read_inventory(invfile)
catalog = read_events(catfile)

stream = RFStream()
with tqdm() as pbar:
    for s in iter_event_data(catalog, inventory, custom_get_waveforms, pbar=pbar):
        stream.extend(s)

stream.write(datafile, 'H5')

