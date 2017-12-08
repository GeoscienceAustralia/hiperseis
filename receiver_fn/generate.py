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

from asdf_util import ASDFUtil

data = os.path.join('data', '')
invfile = data + '7D-inventory.xml'
catfile = data + '7D-catalog.xml'
datafile = data + '7D-rf_profile_data.h5'
rffile = data + '7D-rf_profile_rfs.h5'
profilefile = data + '7D-rf_profile_profile.h5'

inventory = read_inventory(invfile)
catalog = read_events(catfile)

stream = RFStream()
asdfUtil = ASDFUtil('/g/data/ha3/Passive/_ANU/7D(2012-2013)/ASDF/7D(2012-2013).h5')
with tqdm() as pbar:
    for s in iter_event_data(catalog, inventory, asdfUtil.get_waveforms, pbar=pbar):
        stream.extend(s)

stream.write(datafile, 'H5')
