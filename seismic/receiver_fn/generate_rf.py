#!/usr/bin/env python

# AG
# import os.path
# import matplotlib.pyplot as plt
# import numpy as np

from rf import read_rf, RFStream
from rf import IterMultipleComponents
# from rf.imaging import plot_profile_map
# from rf.profile import profile
from tqdm import tqdm

# pylint: disable=invalid-name

RESAMPLE_RATE_HZ = 100
FILTER_BAND_HZ = (0.01, 15.0)
TRIM_START_TIME_SEC = -25.0
TRIM_END_TIME_SEC = 75.0
EXCLUDE_STATION_CODE = ['MIJ2', 'MIL2']

data = read_rf('/g/data/ha3/am7399/shared/OA_event_waveforms_for_rf_20171001T120000-20171015T120000.h5', 'H5')

# exclude bad stations
inc_set = list(set([tr.stats.inclination for tr in data]))
data_filtered = RFStream([tr for tr in data if tr.stats.inclination in inc_set and tr.stats.station
                          not in EXCLUDE_STATION_CODE])

stream = RFStream()
for stream3c in tqdm(IterMultipleComponents(data, 'onset', 3)):
    stream3c.detrend('linear').resample(RESAMPLE_RATE_HZ)
    stream3c.taper(FILTER_BAND_HZ[0])
    stream3c.filter('bandpass', freqmin=FILTER_BAND_HZ[0], freqmax=FILTER_BAND_HZ[1])
    if len(stream3c) != 3:
        continue

    try:
        stream3c.rf()
    except ValueError as e:
        print("ERROR: Failed on stream:\n{}".format(stream3c))
        print(e)
        print("(continuing from error...)")
        continue

    stream3c.trim2(TRIM_START_TIME_SEC, TRIM_END_TIME_SEC, 'onset')
    stream.extend(stream3c)

stream.write('/g/data/ha3/am7399/shared/OA_event_waveforms_for_rf_20171001T120000-20171015T120000_LQT.h5', 'H5')
