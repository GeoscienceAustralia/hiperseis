from rf import read_rf, RFStream
from rf import IterMultipleComponents
from joblib import Parallel, delayed

def do_rf(stream3c):

    stream3c.detrend('linear').resample(100)
    stream3c.taper(0.01)
    stream3c.filter('bandpass', freqmin=0.01, freqmax=15)
    if len(stream3c) != 3:
       return RFStream()
    a1=stream3c[0].stats['asdf']
    a2=stream3c[1].stats['asdf']
    a3=stream3c[2].stats['asdf']
    stream3c[0].stats['asdf']=[]
    stream3c[1].stats['asdf']=[]
    stream3c[2].stats['asdf']=[]

    # LQT receiver functions are default
#   stream3c.rf()
    # ZRT receiver functions must be specified
    stream3c.rf(rotate='NE->RT')

    stream3c[0].stats['asdf']=a1
    stream3c[1].stats['asdf']=a2
    stream3c[2].stats['asdf']=a3
    stream3c.trim2(-25, 75, 'onset')
    return stream3c

print "Lets start the show..."
data = read_rf('DATA/7X-event_waveforms_for_rf.h5', 'H5')
print "Data in..."

'''
# we can exclude bad stations
inc_set = list(set([tr.stats.inclination for tr in data]))
data_filtered = RFStream([tr for tr in data if tr.stats.inclination in inc_set and tr.stats.station not in ['MIJ2', 'MIL2']])
'''

stream = RFStream()

rf_streams=Parallel(n_jobs=30,verbose=1)(map(delayed(do_rf),IterMultipleComponents(data, 'onset', 3)))

for i,rf in enumerate(rf_streams):
    event_id={'eventid':0}
    event_id['eventid']=i
    rf.stats.update(event_id)
    stream.extend(rf)

stream.write('DATA/7X-rf_zrt', 'H5')
print "No worries, mate..."
