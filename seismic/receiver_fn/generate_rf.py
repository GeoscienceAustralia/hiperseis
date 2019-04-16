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

data = read_rf('DATA/7X-event_waveforms_for_rf.h5', 'H5')

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

    # Workaround to preserve H5 header information. (TODO: check if this is fixed in latest obspy,
    # document which version(s) this bug applies to.)
    a1 = stream3c[0].stats['asdf']
    a2 = stream3c[1].stats['asdf']
    a3 = stream3c[2].stats['asdf']
    stream3c[0].stats['asdf'] = []
    stream3c[1].stats['asdf'] = []
    stream3c[2].stats['asdf'] = []

    stream3c.rf()
    stream3c[0].stats['asdf'] = a1
    stream3c[1].stats['asdf'] = a2
    stream3c[2].stats['asdf'] = a3
    stream3c.trim2(TRIM_START_TIME_SEC, TRIM_END_TIME_SEC, 'onset')
    stream.extend(stream3c)

stream.write('DATA/7X-rf_qlt', 'H5')
