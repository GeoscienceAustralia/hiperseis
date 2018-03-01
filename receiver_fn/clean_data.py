import sys
import os
import numpy as np
from obspy.taup import TauPyModel
from obspy.signal.trigger import recursive_sta_lta, trigger_onset
from rf import read_rf, RFStream
# make changes to extract_picks to import this instead of repeating them here
#sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'seismic', 'pickers_integration'))
#from extract_picks import clean_trace, find_best_bounds

model = TauPyModel(model='iasp91')
pphase_search_margin=10

def clean_trace(tr, t1, t2, freqmin=1.0, freqmax=4.9):
    # add bandpass filtering sensitive to the distance d between source and receiver
    tr.trim(t1, t2)
    tr.detrend('linear')
    tr.taper(0.01)
    tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax)

def find_best_bounds(cft, samp_rate):
    bestupper = 2
    bestlower = 0.75
    max_margin = bestupper - bestlower
    leasttrigs = len(trigger_onset(cft, 2, 0.75, max_len=(60*samp_rate), max_len_delete=True))
    for upper in np.linspace(2, 10, 12):
        for lower in [0.875, 0.75, 0.625, 0.5, 0.375]:
            t = trigger_onset(cft, upper, lower, max_len=(60*samp_rate), max_len_delete=True)
            if len(t) > 0 and (upper - lower) > max_margin and len(t) <= leasttrigs:
                leasttrigs = len(t)
                max_margin = upper - lower
                bestupper = upper
                bestlower = lower
    return bestupper, bestlower

def extract_filter_params(trace):
    samp_rate = trace.stats.sampling_rate
    trigs = []
    p_arrivals = model.get_travel_times_geo(trace.stats.event_depth, trace.stats.event_latitude,
                     trace.stats.event_longitude, trace.stats.station_latitude, trace.stats.station_longitude, phase_list=['P',])
    mean_parrival = np.mean([parr.time for parr in p_arrivals])
    for band in [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9)]:
        tr_copy = trace.copy()
        clean_trace(tr_copy, tr_copy.stats.starttime, tr_copy.stats.endtime, freqmin=band[0], freqmax=band[1])
        cft = recursive_sta_lta(tr_copy.data, int(5*samp_rate), int(20*samp_rate))
        upper, lower = find_best_bounds(cft, samp_rate)
        trigs.extend([(onset, tr_copy.stats.channel, upper-lower, band, upper) for onset in trigger_onset(cft, upper, lower, max_len=(60*tr_copy.stats.sampling_rate), max_len_delete=True)])
    trigs = [trig for trig in trigs if abs(trace.stats.event_time + mean_parrival - trace.stats.starttime - (trig[0][0]/samp_rate)) < pphase_search_margin]
    if len(trigs) > 0:
        mintrigdiff = abs(trace.stats.event_time + mean_parrival - trace.stats.starttime - trigs[0][0][0])
        besttrig = trigs[0][0][0]
        best_trig_margin = trigs[0][2]
        best_band = trigs[0][3]
        best_upper = trigs[0][4]
        for trig in trigs:
            if abs(trace.stats.event_time + mean_parrival - trace.stats.starttime - trig[0][0]) <= mintrigdiff and \
                trig[2] >= best_trig_margin and trig[4] >= best_upper:
                mintrigdiff = abs(trace.stats.event_time + mean_parrival - trace.stats.starttime - trig[0][0])
                besttrig = trig[0][0]
                best_trig_margin = trig[2]
                best_band = trig[3]
                best_upper = trig[4]
        return (best_band, best_trig_margin, best_upper)
    else:
        print ('Something went wrong ... ')
    return None

in_streamfile = 'data/7X-rf_profile_data-15deg.h5'
out_streamfile = 'data/7X-rf_profile_data-15deg-out.h5'
st = read_rf(in_streamfile, 'H5')
stz = [tr for tr in st if tr.stats.channel.endswith('Z')]
stn = [tr for tr in st if tr.stats.channel.endswith('N')]
ste = [tr for tr in st if tr.stats.channel.endswith('E')]

output_from_z = []
for trz, trn, tre in zip(stz, stn, ste):
    band, margin, upper = extract_filter_params(trz)
    for tr in [trz, trn, tre]:
        clean_trace(tr, tr.stats.starttime, tr.stats.endtime, freqmin=band[0], freqmax=band[1])
        print('Plot the cleaned trace here')

st_out = stz + stn + ste
st_out_rf = RFStream(st_out)
st_out_rf.write(out_streamfile, 'H5') 
