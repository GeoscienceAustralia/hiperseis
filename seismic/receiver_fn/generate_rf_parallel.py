from rf import read_rf, RFStream
from rf import IterMultipleComponents
from joblib import Parallel, delayed
import numpy as np


def do_rf(stream3c):

    stream3c=stream3c.detrend('linear').interpolate(100)
    stream3c=stream3c.taper(0.01)
    stream3c=stream3c.filter('bandpass', freqmin=0.03, freqmax=0.80, corners=2, zerophase=True)
    if len(stream3c) != 3:
       return RFStream()
    a1 = stream3c[0].stats['asdf']
    a2 = stream3c[1].stats['asdf']
    a3 = stream3c[2].stats['asdf']
    stream3c[0].stats['asdf'] = []
    stream3c[1].stats['asdf'] = []
    stream3c[2].stats['asdf'] = []

    # LQT receiver functions are default
#   stream3c.rf()
    # ZRT receiver functions must be specified
    stream3c.rf(rotate='NE->RT')
#   stream3c.rf(rotate='NE->RT',deconvolve='freq',gauss=2.0)
    '''
    Note the parameters of gaussian pulse and its width where

    Value of "a" | Frequency (hz) at which G(f) = 0.1 |  Approximate Pulse Width (s)

    10                      4.8                                0.50
    5                       2.4                                0.75
    2.5                     1.2                                1.00
    1.25                    0.6                                1.50
    1.0                     0.5                                1.67 (5/3)
    0.625                   0.3                                2.10
    0.5                     0.24                               2.36
    0.4                     0.2                                2.64
    0.2                     0.1                                3.73
    '''

    amax={'amax':np.amax(stream3c[0].data)}
    stream3c[0].stats['asdf'] = a1
    stream3c[0].stats.update(amax)
#   stream3c[0].filter('bandpass', freqmin=0.03, freqmax=1.00, corners=2, zerophase=True)
#   stream3c[0].data=stream3c[0].data*(amax['amax']/np.amax(stream3c[0].data))

    amax={'amax':np.amax(stream3c[1].data)}
    stream3c[1].stats['asdf'] = a2
    stream3c[1].stats.update(amax)
#   stream3c[1].filter('bandpass', freqmin=0.03, freqmax=1.00, corners=2, zerophase=True)
#   stream3c[1].data=stream3c[0].data*(amax['amax']/np.amax(stream3c[0].data))

    amax={'amax':np.amax(stream3c[2].data)}
    stream3c[2].stats['asdf'] = a3
    stream3c[2].stats.update(amax)
#   stream3c[2].filter('bandpass', freqmin=0.03, freqmax=1.00, corners=2, zerophase=True)
#   stream3c[2].data=stream3c[0].data*(amax['amax']/np.amax(stream3c[0].data))

    stream3c.trim2(-25, 75, 'onset')
#   print np.max(stream3c[0].data),np.max(stream3c[1].data),np.max(stream3c[2].data)
    return stream3c

print "Lets start the show..."
#data = read_rf('DATA/7X-event_waveforms_for_rf.h5', 'H5')
data = read_rf('DATA/7X-MA12.h5','H5')
print "Data in..."

'''
# we can exclude bad stations
inc_set = list(set([tr.stats.inclination for tr in data]))
data_filtered = RFStream([tr for tr in data if tr.stats.inclination in inc_set and tr.stats.station not in ['MIJ2', 'MIL2']])
'''

stream = RFStream()

rf_streams = Parallel(n_jobs=-1, verbose=1)(map(delayed(do_rf), IterMultipleComponents(data, 'onset', 3)))

for i, rf in enumerate(rf_streams):
    event_id = {'event_id': 0}
    event_id['event_id'] = i
    for tr in rf:
        tr.stats.update(event_id)
    stream.extend(rf)

stream.write('DATA/7X-MA12-rf_zrt', 'H5')
print "No worries, mate..."
