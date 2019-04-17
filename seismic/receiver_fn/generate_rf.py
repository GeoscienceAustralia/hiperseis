#!/usr/bin/env python

# AG
# import os.path
# import matplotlib.pyplot as plt
import numpy as np

from rf import read_rf, RFStream
from rf import IterMultipleComponents
# from rf.imaging import plot_profile_map
# from rf.profile import profile
from tqdm import tqdm

# pylint: disable=invalid-name

RESAMPLE_RATE_HZ = 100
FILTER_BAND_HZ = (0.03, 0.80)
TAPER_LIMIT = 0.01
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
    stream3c.detrend('linear').interpolate(RESAMPLE_RATE_HZ)
    stream3c.taper(TAPER_LIMIT)
    stream3c.filter('bandpass', freqmin=FILTER_BAND_HZ[0], freqmax=FILTER_BAND_HZ[1], corners=2, zerophase=True)
    if len(stream3c) != 3:
        continue

    try:
        # LQT receiver functions are default
        # stream3c.rf()
        # ZRT receiver functions must be specified
        stream3c.rf(rotate='NE->RT')
        # stream3c.rf(rotate='NE->RT',deconvolve='freq',gauss=2.0)
    except ValueError as e:
        print("ERROR: Failed on stream:\n{}".format(stream3c))
        print(e)
        print("(continuing from error...)")
        continue

    # Note the parameters of gaussian pulse and its width where
    # Value of "a" | Frequency (hz) at which G(f) = 0.1 |  Approximate Pulse Width (s)
    # 10                      4.8                                0.50
    # 5                       2.4                                0.75
    # 2.5                     1.2                                1.00
    # 1.25                    0.6                                1.50
    # 1.0                     0.5                                1.67 (5/3)
    # 0.625                   0.3                                2.10
    # 0.5                     0.24                               2.36
    # 0.4                     0.2                                2.64
    # 0.2                     0.1                                3.73

    amax = {'amax': np.amax(stream3c[0].data)}
    stream3c[0].stats.update(amax)
#   stream3c[0].filter('bandpass', freqmin=0.03, freqmax=1.00, corners=2, zerophase=True)
#   stream3c[0].data = stream3c[0].data*(amax['amax']/np.amax(stream3c[0].data))

    amax = {'amax': np.amax(stream3c[1].data)}
    stream3c[1].stats.update(amax)
#   stream3c[1].filter('bandpass', freqmin=0.03, freqmax=1.00, corners=2, zerophase=True)
#   stream3c[1].data = stream3c[0].data*(amax['amax']/np.amax(stream3c[0].data))  # Is this supposed to be zero index here?

    amax = {'amax': np.amax(stream3c[2].data)}
    stream3c[2].stats.update(amax)
#   stream3c[2].filter('bandpass', freqmin=0.03, freqmax=1.00, corners=2, zerophase=True)
#   stream3c[2].data = stream3c[0].data*(amax['amax']/np.amax(stream3c[0].data))  # Is this supposed to be zero index here?

    stream3c.trim2(TRIM_START_TIME_SEC, TRIM_END_TIME_SEC, 'onset')
    stream.extend(stream3c)


for i, rf in enumerate(stream):
    event_id = {'event_id': i}
    for tr in rf:
        tr.stats.update(event_id)

stream.write('/g/data/ha3/am7399/shared/OA_event_waveforms_for_rf_20171001T120000-20171015T120000_LQT.h5', 'H5')
