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

data = read_rf('data/7X-rf_profile_data-15deg.h5', 'H5')
data_cleaned = read_rf('data/7X-rf_profile_data-15deg-out.h5', 'H5')
inc_set = list(set([tr.stats.inclination for tr in data_cleaned]))
data_filtered = RFStream([tr for tr in data if tr.stats.inclination in inc_set and tr.stats.station not in ['MIJ2', 'MIL2']])
stream = RFStream()
for stream3c in tqdm(IterMultipleComponents(data_filtered, 'onset', 3)):
    stream3c.filter('bandpass', freqmin=0.1, freqmax=1)
    stream3c.trim2(-25, 75, 'onset')
    if len(stream3c) != 3:
        continue
    stream3c.rf()
    stream3c.moveout()
    stream.extend(stream3c)
stream.write('data/7X-rf_profile_rfs_0.1Hz_1Hz', 'H5')
ppoints = stream.ppoints(70)
boxes = get_profile_boxes((-18.4, 139.1), 135, np.linspace(0, 440, 80), width=500)
pstream = profile(stream, boxes)
pstream.write('data/7X-rf_profile_profile-15deg_0.1Hz_1Hz_70km_ppoints.h5', 'H5')
