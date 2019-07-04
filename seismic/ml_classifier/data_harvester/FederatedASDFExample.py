
# coding: utf-8

# In[1]:


import datetime
import glob, os, sys
from os.path import join, exists
import obspy
from obspy.core import Stream, UTCDateTime
from obspy import read, Trace
from obspy.signal.cross_correlation import xcorr
import numpy as np
import matplotlib
import logging
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
import math
from collections import defaultdict
from netCDF4 import Dataset

from obspy.signal.detrend import simple, spline
from obspy.signal.filter import bandpass

from ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

from obspy import UTCDateTime, read_events, read_inventory
from obspy.taup.taup_geo import calc_dist
from obspy.clients.iris import Client as IrisClient
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.signal.trigger import trigger_onset, z_detect, classic_sta_lta, recursive_sta_lta, ar_pick
from obspy.signal.rotate import rotate_ne_rt
from obspy.core.event import Pick, CreationInfo, WaveformStreamID, ResourceIdentifier, Arrival, Event,Origin, Arrival, OriginQuality, Magnitude, Comment


# In[3]:


fds = FederatedASDFDataSet('/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt', variant='db', use_json_db=True,
                           logger=None)


# In[5]:


stations = fds.get_stations('2009-05-17T00:00:00', '2009-05-18T00:00:00', station='QLP')

print(stations)


# In[7]:


s = fds.get_waveforms('AU', 'QLP', '', 'BHE', 
                      '2011-03-15T00:00:00', '2011-03-16T00:00:00',
                      automerge=True, trace_count_threshold=10)
print(s)
if(len(s)):s.plot(outfile='/home/118/rlt118/meme.png')

