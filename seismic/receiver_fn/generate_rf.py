# AG
import os.path
# import matplotlib.pyplot as plt
# import numpy as np

from rf import read_rf, RFStream
from rf import IterMultipleComponents
# from rf.imaging import plot_profile_map
# from rf.profile import profile
from tqdm import tqdm

data = read_rf('DATA/7X-event_waveforms_for_rf.h5', 'H5')

# exclude bad stations
inc_set = list(set([tr.stats.inclination for tr in data]))
data_filtered = RFStream([tr for tr in data if tr.stats.inclination in inc_set and tr.stats.station not in ['MIJ2', 'MIL2']])

stream = RFStream()
for stream3c in tqdm(IterMultipleComponents(data, 'onset', 3)):
    stream3c.detrend('linear').resample(100)
    stream3c.taper(0.01)
    stream3c.filter('bandpass', freqmin=0.01, freqmax=15)
    if len(stream3c) != 3:
        continue
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
    stream3c.trim2(-25, 75, 'onset')
    stream.extend(stream3c)

stream.write('DATA/7X-rf_qlt', 'H5')
