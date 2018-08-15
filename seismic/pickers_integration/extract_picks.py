# Copyright
'''
This script accomplishes the task of picking the p-phase and s-phase arrival times
for waveform data sets that have been ingested as miniseed data into a Seiscomp3 (SC3
for short) server SDS archive.
'''

import click
import os
import numpy as np
import itertools
from obspy import UTCDateTime, read_events, read_inventory
from obspy.taup.taup_geo import calc_dist
from obspy.clients.iris import Client as IrisClient
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.signal.trigger import trigger_onset, z_detect, classic_sta_lta, recursive_sta_lta, ar_pick
from obspy.signal.rotate import rotate_ne_rt
from obspy.core.event import Pick, CreationInfo, WaveformStreamID, ResourceIdentifier, Arrival, Event,\
    Origin, Arrival, OriginQuality, Magnitude, Comment
import glob
import time

def s_phase_exists(evt, pick):
    '''
    Check if the given event contains the s-phase for a given
    p-phase pick.

    :type evt: class:`~obspy.core.event.event.Event`
    :param evt: event object to search for s-pick in.
    :type pick: obspy.core.event.Pick
    :param pick: the p-pick for which the s-pick needs to be checked.

    :rtype: bool
    :return: True if s-pick exist otherwise False
    '''
    exists=False
    filtered_picks = [p for p in evt.picks if p.waveform_id.station_code == pick.waveform_id.station_code]
    for p in filtered_picks:
        if p.phase_hint=='S':
            exists=True
            break
    return exists

def calc_distance(rec, src):
    '''
    Given the source and receiver location, calculate distance.
    Refer: https://github.com/obspy/obspy/blob/master/obspy/taup/taup_geo.py
    '''
    return calc_dist(src.latitude, src.longitude, rec.latitude, rec.longitude, 6371, 0)

def get_backazimuth(stalat, stalon, evtlat, evtlon):
    '''
    get the backazimuth for a given station-event pair

    :type stalat: float
    :param stalat: station latitude in degrees
    :type stalon: float
    :param stalon: station longitude in degrees
    :type evtlat: float
    :param evtlat: event/origin latitude in degrees
    :type evtlon: float
    :param evtlon: event/origin latitude in degrees
    :rtype: dict
    :return: Dictionary containing values for azimuth, backazimuth and distance.
    '''
    az = None
    client = IrisClient()
    try:
        az = client.distaz(stalat=stalat, stalon=stalon, evtlat=evtlat, evtlon=evtlon)
    except:
        print(('IRISClient.distaz failed with error => '+str(e)))
    return az

def createPickObject(net, sta, cha, time, backaz, phasehint, res=0.0, wt=1.0, comments_data=None):
    '''
    create the Pick object from its elements i.e. network, station, channel, time, backazimuth,
    phasehint(phase-type), residual, weight, text

    :type net: str
    :param net: Select the network code on which this pick was identified.
    :type sta: str
    :param sta: Select the station code on which this pick was identified.
    :type cha: str
    :param cha: Select the channel code on which this pick was identified.
    :type time: class:`~obspy.core.utcdatetime.UTCDateTime`
    :param time: The UTC timestamp for pick arrival
    :type backaz: float
    :param backaz: The angle from the event to the station.
    :type phasehint: str
    :param phasehint: Tentative phase identification as specified by the picker. e.g. 'P', 'S', etc
    :type res: float, optional
    :param res: Residual between observed and expected arrival time
    :type wt: float, optional
    :param float: Weight of the arrival time for computation of the associated Origin
    :type comments_data: list of :class:`~obspy.core.event.base.Comment`, optional
    :param comments_data: Additional comments.
    :rtype: (class:`~obspy.core.event.origin.Pick`, float, float)
    :return: a 3-tuple of the newly created Pick object, the residual and the weight that will potentially
        to be part of the bigger Event object
    '''
    if phasehint == 'S':
        if cha.endswith('N') or cha.endswith('1'):
            cha = cha[0:-1]+'R'
        elif cha.endswith('E') or cha.endswith('2'):
            cha = cha[0:-1]+'T'
    comments=[]
    if comments_data:
        comments.append(Comment(text='band = '+str(comments_data[0]), force_resource_id=False))
        comments.append(Comment(text='upper = '+str(comments_data[1]), force_resource_id=False))
        comments.append(Comment(text='margin = '+str(comments_data[2]), force_resource_id=False))
        comments.append(Comment(text='snr = '+str(comments_data[3]), force_resource_id=False))
    return (Pick(resource_id=ResourceIdentifier(id='smi:bilby2008.picker.ga.gov.au/pick/dummytrail'),
                 time=time,
                 waveform_id=WaveformStreamID(network_code=net, station_code=sta, channel_code=cha),
                 methodID=ResourceIdentifier('obspy/rec-sta-lta/aic'),
                 backazimuth=backaz,
                 phase_hint=phasehint,
                 evaluation_mode='automatic',
                 comments=comments,
                 creation_info=CreationInfo(author='niket',
                                 creation_time=UTCDateTime(),
                                 agency_id='niket-ga-picker')),
            res,
            wt)

def clean_trace(tr, trim_starttime, trim_endtime, freqmin=1.0, freqmax=4.9, zerophase=False):
    '''
    Clean the given trace for further processing (in-place)

    :type tr: class`~obspy.core.trace.Trace`
    :param tr: The given trace that needs to be cleaned
    :type trim_starttime: class`~obspy.core.utcdatetime.UTCDateTime`
    :param trim_starttime: The start of the trim window
    :type trim_endtime: class`~obspy.core.utcdatetime.UTCDateTime`
    :param trim_endtime: The end of the trim window
    :type freqmin: float, optional
    :param freqmin: Pass band low corner frequency
    :type freqmax: float, optional
    :param freqmax: Pass band high corner frequency
    :type zerophase: bool
    :param zerophase: If True, apply filter once forwards and once backwards.
    '''
    # add bandpass filtering sensitive to the distance d between source and receiver
    tr.trim(trim_starttime, trim_endtime)
    tr.detrend('linear')
    tr.taper(0.01)
    # added try catch after script aborted due to the error below:
    # ValueError: Selected corner frequency is above Nyquist.
    try:
        tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax, zerophase=zerophase)
    except:
        pass

def find_best_bounds(cft, samp_rate):
    '''
    Refer https://docs.obspy.org/tutorial/code_snippets/trigger_tutorial.html and
    https://docs.obspy.org/packages/autogen/obspy.signal.trigger.trigger_onset.html
    Iterate between different upper (thres1) and lower (thres2) threshold values
    to arrive at the range that indicates the most definitive pick

    :type cft: class:`numpy.ndarray`, dtype=float64
    :param cft: the characteristic function representing the ratio of the STA and LTA
    :type samp_rate: float
    :param samp_rate: the sampling rate of the waveform data on which STA/LTA was done
    :rtype: (float, float)
    :return: the tuple (upper_threshold, lower_threshold) that represents best thresholds
             that result in a definitive pick.
    '''
    bestupper = 1.5
    bestlower = 0.75
    max_margin = bestupper - bestlower
    leasttrigs = len(trigger_onset(cft, 1.5, 0.75, max_len=(60*samp_rate), max_len_delete=True))
    for upper in np.linspace(1.5, 8, 20):
        for lower in [0.875, 0.75, 0.625, 0.5, 0.375, 0.25]:
            t = trigger_onset(cft, upper, lower, max_len=(60*samp_rate), max_len_delete=True)
            if len(t) > 0 and (upper - lower) > max_margin and len(t) <= leasttrigs:
                leasttrigs = len(t)
                max_margin = upper - lower
                bestupper = upper
                bestlower = lower
    return bestupper, bestlower

def pick_phase(network, station, prefor, phase='P', p_Pick=None):
    '''
    Pick the given station for a given phase around the given preferred origin

    :type network: class`~obspy.core.inventory.network.Network`
    :param network: the network to which the given station belongs
    :type station: class`~obspy.core.inventory.station.Station`
    :param station: the station that needs to be picked
    :type prefor: class`~obspy.core.event.origin.Origin`
    :param prefor: the origin whose time needs to be taken as reference when
        looking for the given phase arrival
    :type phase: str, optional
    :param phase: the phase that needs to be picked e.g. 'P', 'S', etc
    :type p_Pick: class`~obspy.core.event.origin.Pick`
    :param p_Pick: the p-phase arrival, if known, when searching for s-phase
    :rtype: class`~obspy.core.event.origin.Pick`
    :return: the extracted Pick object after identifying the phase arrival
    '''
    return_pick = None
    snr = 0.0
    #p_bands = [(1, 6), (0.3, 2.3), (0.5, 2.5), (0.8, 2.8), (1, 3), (2, 4), (3, 5), (4, 6), (0.3, 1.3), (0.5, 1.5), (0.8, 1.8), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1.5, 2.5), (2.5, 3.5), (3.5, 4.5), (4.5, 5.5)]
    #s_bands = [(0.5, 2), (1, 2), (0.2, 1.5), (0.3, 2.0), (0.3, 1), (0.3, 0.7), (0.05, 0.3), (0.05, 1)]
    p_bands = [(0.5, 2.0), (0.3, 4.0)]
    s_bands = [(0.05, 0.2), (0.05, 1.0)]
    lookback = plookback if phase=='P' else slookback
    lookahead = plookahead if phase=='P' else slookahead
    arrivals = model.get_travel_times_geo(prefor.depth/1000, prefor.latitude,
             prefor.longitude, station.latitude, station.longitude, phase_list=("ttbasic",))
    target_arrivals = [ar for ar in arrivals if ar.phase.name.lower() == phase.lower()]
    target_arrivals.sort(key=lambda x: x.time)
    mean_target_arrival = np.mean([arr.time for arr in target_arrivals])
    distance = calc_distance(station, prefor)
    az = get_backazimuth(stalat=station.latitude, stalon=station.longitude, evtlat=prefor.latitude, evtlon=prefor.longitude)
    if not az:
        return None
    print(('station=> ' + str(station.code) + ' distance =>'+str(distance)))
    traces = []
    if len(target_arrivals) > 0:
        try:
            trim_starttime = prefor.time+target_arrivals[0].time-lookback
            trim_endtime = prefor.time+target_arrivals[-1].time+lookahead
            stz = client.get_waveforms(network.code, station.code, '*', 'BHZ,SHZ,HHZ', trim_starttime, trim_endtime)
            if stz and stz[0]:
                if trim_starttime < stz[0].stats.starttime:
                    trim_starttime = stz[0].stats.starttime
                if trim_endtime > stz[0].stats.endtime:
                    trim_endtime = stz[0].stats.endtime
            if not isclose(stz[0].stats.sampling_rate, 20.0):
                stz.resample(20.0)
            stz_raw = stz.copy()
            traces.append(stz[0])
            mult_const = 10
            if phase=='S':
                stn = client.get_waveforms(network.code, station.code, '*', 'BHN,SHN,HHN,BH1,SH1', trim_starttime, trim_endtime)
                ste = client.get_waveforms(network.code, station.code, '*', 'BHE,SHE,HHZ,BH2,SH2', trim_starttime, trim_endtime)
                if not isclose(stn[0].stats.sampling_rate, 20.0) or not isclose(ste[0].stats.sampling_rate, 20.0):
                    stn.resample(20.0)
                    ste.resample(20.0)
                stn.trim(trim_starttime, trim_endtime)
                ste.trim(trim_starttime, trim_endtime)
                r_t = rotate_ne_rt(stn[0].data, ste[0].data, az['backazimuth'])
                stn[0].data = r_t[0]
                ste[0].data = r_t[1]
                stn_raw = stn.copy()
                ste_raw = ste.copy()
                traces.append(stn[0])
                traces.append(ste[0])
                mult_const = 20
        except Exception as e:
            print(('FDSN client could not retrieve waveform for network='+str(network.code)+' and sta='+station.code+' around time='+str(prefor.time+mean_target_arrival)))
            print((str(e)))
            return None
        if phase == 'S':
            if stz == None or stn == None or ste == None or stz[0] == None or stn[0] == None or ste[0] == None or \
                stz[0].stats.sampling_rate != stn[0].stats.sampling_rate or \
                ste[0].stats.sampling_rate != stn[0].stats.sampling_rate:
                return None
        if stz and stz[0]:
            samp_rate = stz[0].stats.sampling_rate
            if phase == 'S' and p_Pick and p_Pick.time:
                if trim_starttime < p_Pick.time:
                    trim_starttime = p_Pick.time + 5
            trigs= []
            aic_arr = []
            bands = p_bands if phase == 'P' else s_bands
            for trace, band in itertools.product(traces, bands):
                tr_copy = trace.copy()
                if len(set(tr_copy.data)) == 1:
                    continue
                clean_trace(tr_copy, trim_starttime, trim_endtime, band[0], band[1])
                print(('band[0] => ' +str(band[0]) + ' band[1] => ' + str(band[1])))
                cft = recursive_sta_lta(tr_copy.data, int(5*samp_rate), int(mult_const*samp_rate))
                upper, lower = find_best_bounds(cft, tr_copy.stats.sampling_rate)
                trigs.extend([(onset, tr_copy.stats.channel, upper-lower, band, upper) for onset in trigger_onset(cft, upper, lower, max_len=(60*tr_copy.stats.sampling_rate), max_len_delete=True)])
            search_margin = (pphase_search_margin if distance < 10 else pphase_search_margin+5) if phase == 'P' else sphase_search_margin
            margin_threshold = 1.0 if phase == 'P' else 2.0
            trigs = [t for t in trigs if abs(prefor.time + mean_target_arrival - trim_starttime - (t[0][0]/samp_rate)) < search_margin]
            if len(trigs) > 0:
                mintrigdiff = abs(prefor.time + mean_target_arrival - trim_starttime - (trigs[0][0][0]/samp_rate))
                besttrig = trigs[0][0][0]
                best_cha = trigs[0][1]
                best_margin = trigs[0][2]
                best_band = trigs[0][3]
                best_upper = trigs[0][4]
                for trig in trigs:
                    if abs(prefor.time + mean_target_arrival - trim_starttime - (trig[0][0]/samp_rate)) <= mintrigdiff and \
                        trig[2] >= best_margin:
                        mintrigdiff = abs(prefor.time + mean_target_arrival - trim_starttime - (trig[0][0]/samp_rate))
                        besttrig = trig[0][0]
                        best_cha = trig[1]
                        best_margin = trig[2]
                        best_band = trig[3]
                        best_upper = trig[4]

                # add SNR to comments
                comments_data=(best_band, best_upper, best_margin)
                if best_upper > 1.5 and best_margin > margin_threshold:
                    tr_copy = stz[0].copy() if best_cha.endswith('Z') else (stn[0].copy() if best_cha.endswith('N') else ste[0].copy())
                    clean_trace(tr_copy, trim_starttime, trim_endtime, best_band[0], best_band[1])
                    aic, aic_deriv, aic_deriv_cleaned = calc_aic(tr_copy, phase)
                    aic[0:int(10*samp_rate)]=aic[np.argmax(aic)]
                    aic[-int(10*samp_rate):-1]=aic[np.argmax(aic)]
                    pick_index = np.argmin(aic)
                    aic_deriv_cleaned[0:int(10*samp_rate)]=0
                    aic_deriv_cleaned[-int(10*samp_rate):-1]=0
                    pick_index_deriv = np.argmax(aic_deriv_cleaned)
                    # internal plotting functions to access the variables
                    def plot_zerophase_compare():
                        import matplotlib.pyplot as plt
                        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(30,15))
                        start = 40 if phase == 'P' else (70 if distance > 10 else 50)
                        end = 80 if phase == 'P' else 160
                        tr_copy1=tr.copy()
                        tr_copy2=tr.copy()
                        tr_copy3=tr.copy()
                        tr_copy1.trim(trim_starttime+start, trim_endtime-end)
                        clean_trace(tr_copy2, trim_starttime+start, trim_endtime-end, best_band[0], best_band[1])
                        clean_trace(tr_copy3, trim_starttime+start, trim_endtime-end, best_band[0], best_band[1], zerophase=True)
                        axes.plot(tr_copy1, color='grey')
                        axes.text(0, int(min(tr_copy1.data)), 'RAW, sampling_rate='+str(tr_copy1.stats.sampling_rate), fontsize=12, color='grey')
                        axes.plot(tr_copy2, color='green')
                        axes.text(0, int(min(tr_copy2.data)+20), 'single pass filter, band='+str(best_band), fontsize=12, color='green')
                        axes.plot(tr_copy3, color='blue')
                        axes.text(0, int(min(tr_copy3.data)+40), 'two pass filter(zero phase is True), band='+str(best_band), fontsize=12, color='blue')
                        axes.legend()
                        plt.tight_layout()
                        fig.savefig('zerophase_'+network.code+'_'+station.code+'_'+tr.stats.channel+'_'+phase+'.png')
                        plt.close('all')

                    def plot_aic(snr, theo_trig):
                        import matplotlib.pyplot as plt
                        fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(30,15))
                        tr.trim(trim_starttime, trim_endtime)
                        axes[0].plot(tr.data, color='grey')
                        replacement = 'T' if tr.stats.channel.endswith('E') else ('R' if tr.stats.channel.endswith('N') else 'Z')
                        chan = tr.stats.channel[0:2]+replacement
                        axes[0].text(0, int(min(tr.data)), network.code+' '+station.code+' '+chan+' '+phase+' RAW Distance='+str(distance), fontsize=12)
                        axes[1].plot(list(range(theo_trig)), tr_copy.data[:theo_trig], color='blue')
                        #axes[1].plot(range(theo_trig, theo_trig+20), tr_copy.data[theo_trig:theo_trig+20], color='black')
                        for xc in range(theo_trig, theo_trig+10):
                            axes[1].axvline(x=xc, color='black')
                        axes[1].plot(list(range(theo_trig+10, tr_copy.stats.npts)), tr_copy.data[theo_trig+10:tr_copy.stats.npts], color='blue')
                        axes[1].text(0, int(min(tr_copy.data)), 'Filtered band='+str(best_band)+' SNR='+str(snr)+' sampling_rate='+str(tr.stats.sampling_rate), fontsize=12)
                        axes[1].text(theo_trig-100, int(min(tr_copy.data)), 'theoretical tt', fontsize=12)
                        axes[2].plot(aic, color='green')
                        axes[2].text(0, int(min(aic)), 'AIC', fontsize=12)
                        axes[3].plot(aic_deriv, color='red')
                        axes[3].text(0, int(min(aic_deriv)), 'AIC derivative', fontsize=12)
                        axes[4].plot(aic_deriv_cleaned, color='magenta')
                        axes[4].text(0, int(min(aic_deriv_cleaned)), 'AIC global minimum. residual = '+str((np.argmax(aic_deriv_cleaned)-theo_trig)/samp_rate), fontsize=12)
                        plt.tight_layout()
                        fig.savefig(network.code+'_'+station.code+'_'+tr.stats.channel+'_'+phase+'.png')
                        plt.close('all')
                    def plot_onsets(deriv=False):
                        import matplotlib.pyplot as plt
                        fig, axes = plt.subplots(nrows=len(s_bands)+1, ncols=3)
                        tr_n = stn[0].copy()
                        tr_e = ste[0].copy()
                        tr_z = stz[0].copy()
                        tr_n.trim(trim_starttime, trim_endtime)
                        tr_e.trim(trim_starttime, trim_endtime)
                        tr_z.trim(trim_starttime, trim_endtime)
                        axes[0, 0].plot(tr_n.data, color='grey')
                        axes[0, 0].text(int(tr_n.stats.npts*0.5), int(max(tr_n.data)*0.7), network.code+' '+station.code+' R RAW')
                        axes[0, 1].plot(tr_e.data, color='grey')
                        axes[0, 1].text(int(tr_n.stats.npts*0.5), int(max(tr_n.data)*0.7), network.code+' '+station.code+' T RAW')
                        axes[0, 2].plot(tr_z.data, color='grey')
                        axes[0, 2].text(int(tr_z.stats.npts*0.5), int(max(tr_z.data)*0.7), network.code+' '+station.code+' Z RAW')
                        for ind, band in enumerate(s_bands):
                            for col_index, comp in enumerate(['r', 't', 'z']):
                                tr_copy = stn[0].copy() if comp=='r' else (ste[0].copy() if comp=='t' else stz[0].copy())
                                tr_copy.trim(trim_starttime, trim_endtime)
                                clean_trace(tr_copy, trim_starttime, trim_endtime, band[0], band[1])
                                aic, aic_deriv, aic_deriv_cleaned = calc_aic(tr_copy, phase)
                                aic[0:int(10*samp_rate)]=aic[np.argmax(aic)]
                                aic[-int(10*samp_rate):-1]=aic[np.argmax(aic)]
                                pick_index = np.argmin(aic)
                                aic_deriv_cleaned[0:int(10*samp_rate)]=0
                                aic_deriv_cleaned[-int(10*samp_rate):-1]=0
                                pick_index_deriv = np.argmax(aic_deriv_cleaned)
                                index = pick_index_deriv if deriv else pick_index
                                snr_plot = calc_snr(tr_copy, trim_starttime + (index/samp_rate))
                                axes[ind+1, col_index].plot(list(range(int(index))), tr_copy.data[:int(index)], color='blue')
                                axes[ind+1, col_index].plot(list(range(int(index), int(index)+20)), tr_copy.data[int(index):int(index)+20], color='red')
                                axes[ind+1, col_index].plot(list(range(int(index)+20, tr_copy.stats.npts)), tr_copy.data[int(index)+20:tr_copy.stats.npts], color='blue')
                                axes[ind+1, col_index].text(int(tr_copy.stats.npts*0.5), int(max(tr_copy.data)*0.7), network.code+' '+station.code+' '+comp.upper()+' ('+str(band[0])+'Hz - '+str(band[1])+'Hz) SNR='+str(int(snr_plot))+' pick_index = '+str(index))
                        print('Break here while debugging and run plt.show()')

                    theoretical_trig = int((prefor.time + mean_target_arrival - trim_starttime)*samp_rate)
                    if theoretical_trig > 0 and abs(pick_index_deriv - theoretical_trig)/samp_rate < search_margin:
                        snr = calc_snr(tr_copy, trim_starttime + (pick_index_deriv/samp_rate))
                        #plot_onsets(deriv=True)
                        #plot_aic(snr, theoretical_trig)
                        res = trim_starttime + (pick_index_deriv/samp_rate) - prefor.time - mean_target_arrival
                        comments_data=comments_data+(snr,)
                        return_pick = createPickObject(network.code, station.code, best_cha, trim_starttime+(pick_index_deriv/samp_rate), az['backazimuth'] if az else None, phase, res, comments_data=comments_data)
                        print((phase+'-pick added'))
                    else:
                        print(('theo_trig => '+str(theoretical_trig)+'pick_index_deriv => ' + str(pick_index_deriv) + ' besttrig => ' + str(besttrig) + '. Investigate waveforms!'))

    if snr > 0 and (snr < 1.0 or snr > snr_max):
        print(('The calculated SNR => ' + str(snr) + '. Discarding this pick!'))
        return None
    return return_pick

def calc_snr(tr, time):
    '''
    Calculate the Signal-to-Noise ratio of the give trace around a given time

    :type tr: class`~obspy.core.trace.Trace`
    :param tr: the given trace for which the SNR is required
    :type time: class:`~obspy.core.utcdatetime.UTCDateTime`
    :param time: the timestamp in the input trace around with the SNR is required
    :rtype: float
    :return: the SNR
    '''
    if not tr or not time:
        print('Either of trace or the pick time are None.')
        return -1
    if time < tr.stats.starttime or time > tr.stats.endtime:
        print('The pick time lies outside the trace times')
        return -1
    tr_left = tr.copy()
    tr_right = tr.copy()
    tr_left.trim(time - 5, time)
    tr_right.trim(time, time + 5)
    return np.std(tr_right.data)/np.std(tr_left.data)

def isclose(a, b, rel_tol=1e-06, abs_tol=0.0):
    '''
    method to check if 2 float values are equal with a degree of tolerance

    refer: https://stackoverflow.com/questions/5595425/what-is-the-best-way-to-compare-floats-for-almost-equality-in-python
    '''
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def calc_aic(tr, phase):
    '''
    Calculate the AIC output for the given trace for a given phase type

    :type tr: class`~obspy.core.trace.Trace`
    :param tr: the given trace for which AIC needs to be calculated
    :type phase: str
    :param phase: the phase for which AIC is required
    :rtype: (np.array of float64, np.array of float64, np.array of float64)
    :return: a 3-tuple of:
        1. the AIC output
        2. the first derivative of AIC output
        3. the second derivative of AIC output
    '''
    npts = tr.stats.npts
    data = tr.data
    margin = int(tr.stats.sampling_rate*5) if phase=='P' else int(tr.stats.sampling_rate*10)
    aic = np.zeros(npts)
    for k in range(npts-2,0,-1):
        a = k*np.log10(np.std(data[:k])**2)+(npts-k-1)*np.log10(np.std(data[k:])**2)
        if a == -float('inf'):
            a = aic[k+1]
        aic[k] = a
    aic[0] = aic[1]
    aic[-1] = aic[-2]

    aic_deriv = []
    aic_deriv_cleaned = []
    for i in range(npts-1):
        aic_deriv.append(aic[i+1] - aic[i])
        if i < margin or i >= (npts - margin):
            aic_deriv_cleaned.append(min(aic_deriv))
        else:
            if aic[i - margin] < aic[i] < aic[i + margin] or \
                aic[i + margin] < aic[i] < aic[i - margin]:
                aic_deriv_cleaned.append(min(aic_deriv))
            elif aic[i] > aic[i - margin] and aic[i] > aic[i + margin]:
                aic_deriv_cleaned.append(min(aic_deriv))
            else:
                aic_deriv_cleaned.append(aic[i+1] - aic[i])

    return np.array(aic), np.array(aic_deriv), np.array(aic_deriv_cleaned)

def pick_p_phase(prefor, network=None, station=None, pick=None):
    '''
    A p-phase specific wrapper function around the generic pick_phase function

    :type network: class`~obspy.core.inventory.network.Network`
    :param network: the network to which the given station belongs
    :type prefor: class`~obspy.core.event.origin.Origin`
    :param prefor: the origin whose time needs to be taken as reference when
        looking for the given p-phase arrival
    :type station: class`~obspy.core.inventory.station.Station`, optional
    :param station: the station that needs to be picked
    :type pick: class`~obspy.core.event.origin.Pick`, optional
    :param phase: the pick whose station needs to be picked
    :rtype: ((class`~obspy.core.event.origin.Pick`, float, float), class`Network`, class`Station`)
    :return: 3-tuple of the extracted p-Pick tuple, the residual and the weight
        after identifying the phase arrival
    '''
    p_pick = None
    if station and pick:
        print('Only one of station or pick needs to be passed in; not both.')
        return None
    if station and network:
        p_pick = pick_phase(network, station, prefor)
        if p_pick:
            return (p_pick, network, station)
    elif pick:
        try:
            pick_st = client.get_stations(network=pick.waveform_id.network_code, station=pick.waveform_id.station_code)
        except:
            print(('FDSN client could not retrieve station info for sta='+pick.waveform_id.station_code+' and network='+pick.waveform_id.network_code))
            return None
        if not pick_st:
            print(('FDSN client returned an empty response for station query for station='+pick.waveform_id.station_code+' and network='+pick.waveform_id.network_code))
            return None
        p_pick = pick_phase(pick_st.networks[0], pick_st.networks[0].stations[0], prefor)
        if p_pick:
            return (p_pick, pick_st.networks[0], pick_st.networks[0].stations[0])
    else:
        print('Both the passed in station as well as the passed in pick are None')
        return None

def pick_s_phase(p_phase, net_st_list, prefor):
    '''
    A s-phase specific wrapper function around the generic pick_phase function

    :type p_phase: class`~obspy.core.event.origin.Pick`
    :param p_phase: the Pick object that was successfully created during p-phase picking
        for which the s-phase needs to be not extracted
    :type net_st_list: list of tuples of class`~~obspy.core.inventory.network.Network` and
        class`~~obspy.core.inventory.station.Station`
    :param net_st_list: list of tuples of network and station that needs to be searched for
        the network and station combination that corresponds to the p_phase Pick object
    :type prefor: class`~obspy.core.event.origin.Origin`
    :param prefor: the origin whose time needs to be taken as reference when
        looking for the given p-phase arrival
    :rtype: ((class`~obspy.core.event.origin.Pick`, float, float), class`Network`, class`Station`)
    :return: 3-tuple of the extracted s-Pick tuple, the network and the station after
        identifying the phase arrival
    '''
    if not p_phase or not p_phase[0] or len(net_st_list) < 1:
        return None
    ret_net_st_lst = [(n, s) for (n, s) in net_st_list if s.code == p_phase[0].waveform_id.station_code]
    if len(ret_net_st_lst) == 0:
        return None
    net, st = ret_net_st_lst[0]
    s_pick = pick_phase(net, st, prefor, phase='S', p_Pick=p_phase[0])
    if s_pick:
        return (s_pick, net, st)
    else:
        return None

def pick_phases(event, inventory=None):
    '''
    Pick the relevant stations for a given event for p-phase as well as s-phase arrivals
    :type event: class`~obspy.core.event.event.Event`
    :param event: the Event object which contains the preferred origin for which we want
        relevant stations picked
    :type inventory: class`~obspy.core.inventory.inventory.Inventory`, optional
    :param inventory: the Inventory object that contains the stations to be picked for this
        event. If not passed in, the stations that are part of the p-phase picks in the event
        object are picked for p-phase arrivals and the successful stations ofr s-phase
        arrivals.
    :rtype: (list of obspy.core.event.origin.Pick, list of obspy.core.event.origin.Pick,
        list of obspy.core.inventory.station.Station)
    :return: tuple of 3 items:
        1. list of the picked p-phase Pick objects
        2. list of the picked s-phase Pick objects
        3. list of Station objects relevant to this event's picking
    '''
    if not event:
        return None, None
    prefor = event.preferred_origin() or event.origins[0]

    p_phases = []
    p_stations = []
    p_pick = None
    if inventory:
        for net in inventory.networks:
            for st in net.stations:
                p_stations.append((net, st))
                p_pick = pick_phase(net, st, prefor)
                if p_pick:
                    p_phases.append(p_pick)
    else:
        for pick in [p for p in event.picks if p.phase_hint=='P']:
            try:
                pick_st = client.get_stations(network=pick.waveform_id.network_code, station=pick.waveform_id.station_code)
            except:
                print ('FDSN client could not retrieve station info for sta='+pick.waveform_id.station_code+' and network='+pick.waveform_id.network_code)
                continue
            if not pick_st:
                print ('FDSN client returned an empty response for station query for station='+pick.waveform_id.station_code+' and network='+pick.waveform_id.network_code)
                continue
            p_stations.append((pick_st.networks[0], pick_st.networks[0].stations[0]))
            p_pick = pick_phase(pick_st.networks[0], pick_st.networks[0].stations[0], prefor)
            if p_pick:
                p_phases.append(p_pick)
    s_phases = []
    s_pick = None
    net_st_list = []
    if inventory:
        for n in inventory.networks:
            for s in n.stations:
                net_st_list.append((n, s))
    for ppick in p_phases:
        if inventory:
            net, st = [(n, s) for (n, s) in net_st_list if s.code == ppick[0].waveform_id.station_code][0]
            s_pick = pick_phase(net, st, prefor, phase='S', p_Pick=ppick[0])
        else:
            s_stn = [s for s in p_stations if s[1].code == ppick[0].waveform_id.station_code][0]
            s_pick = pick_phase(s_stn[0], s_stn[1], prefor, phase='S', p_Pick=ppick[0])
        if s_pick:
            s_phases.append(s_pick)
    return p_phases, s_phases, p_stations

def pick_phases_parallel(event, inventory=None):
    '''
    Pick the relevant stations for a given event for p-phase as well as s-phase arrivals.
    Run parallel processes to perform the picking, for both p-phases as well as s-phases.

    :type event: class`~obspy.core.event.event.Event`
    :param event: the Event object which contains the preferred origin for which we want
        relevant stations picked
    :type inventory: class`~obspy.core.inventory.inventory.Inventory`, optional
    :param inventory: the Inventory object that contains the stations to be picked for this
        event. If not passed in, the stations that are part of the p-phase picks in the event
        object are picked for p-phase arrivals and the successful stations ofr s-phase
        arrivals.
    :rtype: (list of obspy.core.event.origin.Pick, list of obspy.core.event.origin.Pick,
        list of obspy.core.inventory.station.Station)
    :return: tuple of 3 items:
        1. list of the picked p-phase Pick objects
        2. list of the picked s-phase Pick objects
        3. list of Station objects relevant to this event's picking
    '''
    if not event:
        return None, None
    prefor = event.preferred_origin() or event.origins[0]

    p_phases = []
    p_stations = []
    nets = []
    stas = []
    p_pick = None
    requests = []
    count = 0
    pick_station_tuple_arr = []
    start_time = time.time()
    if inventory:
        for net in inventory.networks:
            for st in net.stations:
                requests.append(comm.isend((prefor, net, st, None), dest=count%(size-1)+1))
                count += 1
    else:
        for pick in event.picks:
            if pick.phase == 'P':
                requests.append(comm.isend((prefor, None, None, pick), dest=count%(size-1)+1))
                count += 1
    MPI.Request.Waitall(requests)
    for i in range(count):
        pick_station_tuple_arr.append(comm.recv(source=i%count+1))
    print(("--- time taken to run p-picking in parallel %s seconds ---" % (time.time() - start_time)))
    pick_station_tuple_arr = [x for x in pick_station_tuple_arr if x]
    if len(pick_station_tuple_arr) > 0:
        p_phases, nets, stas = list(zip(*pick_station_tuple_arr))
        p_phases = list(p_phases)
        nets = list(nets)
        stas = list(stas)
        p_stations = list(zip(nets, stas))

    s_phases = []
    s_pick = None
    count = 0
    requests = []
    if inventory:
        for n in inventory.networks:
            for s in n.stations:
                p_stations.append((n, s))
    pick_station_tuple_arr = []
    start_time = time.time()
    for ppick in p_phases:
        requests.append(comm.isend((ppick, p_stations, prefor), dest=count%(size-1)+1))
        count += 1
    MPI.Request.Waitall(requests)
    for i in range(count):
        pick_station_tuple_arr.append(comm.recv(source=i%count+1))
    print(("--- time taken to run s-picking in parallel %s seconds ---" % (time.time() - start_time)))
    pick_station_tuple_arr = [x for x in pick_station_tuple_arr if x]
    if len(pick_station_tuple_arr) > 0:
        s_phases, nets, stas = list(zip(*pick_station_tuple_arr))
        s_phases = list(s_phases)

    return p_phases, s_phases, p_stations

def add_picks_to_event(event, picks, stations):
    '''
    Add the given picks to the given Event object (in-place)

    :type event: class`~obspy.core.event.event.Event`
    :param event: The event object to which the picks need to be added
    :type picks: list of obspy.core.event.origin.Pick
    :param picks: the list of p-phase Pick objects we want added to the Event object
    :type stations: list of obspy.core.inventory.station.Station
    :param stations: list of station objects for this event for ready reference
    '''
    for pick in picks:
        st = [s for s in stations if s[1].code == pick[0].waveform_id.station_code][0][1]
        stalat = st.latitude
        stalon = st.longitude
        prefor = event.preferred_origin() or event.origins[0]
        evtlat = prefor.latitude
        evtlon = prefor.longitude
        az = get_backazimuth(stalat=stalat, stalon=stalon, evtlat=evtlat, evtlon=evtlon)
        arr = Arrival(resource_id=ResourceIdentifier(id=pick[0].resource_id.id+"#"),
                      pick_id=pick[0].resource_id,
                      phase=pick[0].phase_hint,
                      azimuth=az['azimuth'],
                      distance=az['distance'],
                      time_residual=pick[1],
                      time_weight=pick[2],
                      earth_model_id=ResourceIdentifier('quakeml:niket.ga.gov.au/earthmodel/iasp91'),
                      creation_info=CreationInfo(author='niket',
                                creation_time=UTCDateTime(),
                                agency_id='niket-ga-picker'))
        event.picks.append(pick[0])
        event.preferred_origin().arrivals.append(arr)

def add_ordered_pick_ids(picks):
    '''
    Order the resource ids of the pick objects that were returned by the parallel running jobs
    as part of pick_phases_parallel function. This is done in-place.

    :type picks: list of obspy.core.event.origin.Pick
    :param picks: the list of Pick objects whose resource ids need to be numbered as per current
        running pick count
    '''
    global pick_Count

    for p in picks:
        temp_id = p[0].resource_id.id
        temp_id = temp_id.replace('dummytrail', str(pick_Count))
        p[0].resource_id = ResourceIdentifier(id=temp_id)
        pick_Count += 1

def createEventObject(evt, p_picks, s_picks, stations):
    '''
    create the Event object from its elements i.e. the preferred origin, p-phase picks and s-phase picks.

    :type evt: class`~obspy.core.event.event.Event`
    :param evt: the reference event from which the preferred origin is used to create the new event
    :type p_picks: list of obspy.core.event.origin.Pick
    :param p_picks: the list of p-phase Pick objects we want added to the Event object being created
    :type s_picks: list of obspy.core.event.origin.Pick
    :param s_picks: the list of s-phase Pick objects we want added to the Event object being created
    :type stations: list of obspy.core.inventory.station.Station
    :param stations: list of station objects for this event for ready reference
    :rtype: class`~obspy.core.event.event.Event`
    :return: the newly created Event object containing the p-phase and s-phase picks
    '''
    global event_Count
    global origin_Count

    if len(p_picks) == 0:
        return None

    add_ordered_pick_ids(p_picks + s_picks)

    prefor = evt.preferred_origin() or evt.origins[0]
    magPresent = True
    if evt.preferred_magnitude() or len(evt.magnitudes) > 0:
        prefmag = evt.preferred_magnitude() or evt.magnitudes[0]
    else:
        magPresent = False
    creation_info = CreationInfo(author='niket_obspy_aic_picker',
                                 creation_time=UTCDateTime(),
                                 agency_uri=ResourceIdentifier(id='smi:niket.ga.gov.au/ga-picker'),
                                 agency_id='ga')

    origin = Origin(resource_id=ResourceIdentifier(id='smi:bilby2008.ga.gov.au/origin/'+str(origin_Count)),
                    time=prefor.time,
                    longitude=prefor.longitude,
                    latitude=prefor.latitude,
                    depth=prefor.depth,
                    method_id=ResourceIdentifier(id='niket_obspy_aic_picker'),
                    earth_model_id=ResourceIdentifier(id='iasp91'),
                    quality=OriginQuality(associated_phase_count=len(p_picks)+len(s_picks),
                            used_phase_count=len(p_picks)+len(s_picks)),
                    evaluation_mode='automatic',
                    creation_info=creation_info)

    arrList = []
    if magPresent:
        magnitude = Magnitude(resource_id=ResourceIdentifier(id='smi:bilby2008.ga.gov.au/origin/'+str(origin_Count)+'#netMag.Mb'),
                              mag=prefmag.mag,
                              magnitude_type=prefmag.magnitude_type,
                              origin_id=ResourceIdentifier(id='smi:bilby2008.ga.gov.au/origin/'+str(origin_Count)),
                              creation_info=creation_info)
    origin_Count += 1


    event = Event(resource_id=ResourceIdentifier(id='smi:bilby2008.picker.ga.gov.au/event/'+str(event_Count)),
                  creation_info=creation_info,
                  event_type='earthquake')
    event_Count += 1

    event.origins = [origin, ]
    if magPresent:
        event.magnitudes = [magnitude, ]
        event.preferred_magnitude_id = magnitude.resource_id

    event.preferred_origin_id = origin.resource_id

    add_picks_to_event(event, p_picks+s_picks, stations)

    return event

def dedup(inventory):
    '''
    remove duplicate stations, if present, from an inventory object (in-place)

    :type inventory: class`~obspy.core.inventory.inventory.Inventory`
    :param inventory: the given inventory object
    '''
    net_sta_list_uniq = []
    for net in inventory:
        ind_list = []
        for index, sta in enumerate(net):
            if net.code+'.'+sta.code in net_sta_list_uniq:
                ind_list.append(index)
            else:
                net_sta_list_uniq.append(net.code+'.'+sta.code)
        ind_list.sort(reverse=True)
        for ind in ind_list:
            del net.stations[ind]

@click.command()
@click.argument('sc3_host')
@click.argument('year')
@click.option('--input_folder', default=os.getcwd(), help='The input folder where ISC event catalog xml input files are kept.')
@click.option('--output_folder', default=os.getcwd(), help='The output folder where generated xml files are generated.')
@click.option('--parallel', is_flag=True, default=False, help='Needs to be run in conjunction with MPI.(DO NOT USE. NOT TESTED.)')
def main(sc3_host, year, input_folder, output_folder, parallel):
    '''
    This script accomplishes the task of picking the p-phase and s-phase arrival times for given set of events and waveform data sets that have been ingested as miniseed data into a Seiscomp3 (SC3 for short) server SDS archive.

Arguments:\n
    SC3_HOST	: the SC3 FDSN server where the permanent station waveform data is hosted\n
    YEAR	: the year to which the input event catalog corresponds
    '''
    global comm
    global size
    global rank
    global client
    if parallel:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
    client = Client('http://'+sc3_host)
    #evtfiles = glob.glob('/home/ubuntu/engdahl/*.xml')
    evtfiles = glob.glob(os.path.join(input_folder, '*.xml'))
    #outdir = '/home/ubuntu/bilby-out-final/2008'
    #inv = read_inventory('/home/ubuntu/7W_dummy_resp.xml')
    if rank and rank > 0:
        prefor_net_st = None
        while True:
            time.sleep(0.5)
            prefor_net_st_pick = comm.recv(source = 0)
            pick_station_tuple = None
            if len(prefor_net_st_pick) == 4 and isinstance(prefor_net_st_pick[0], Origin):
                pick_station_tuple = pick_p_phase(prefor_net_st_pick[0],
                                                  prefor_net_st_pick[1],
                                                  prefor_net_st_pick[2],
                                                  prefor_net_st_pick[3])
            elif len(prefor_net_st_pick) == 3 and isinstance(prefor_net_st_pick[0][0], Pick):
                pick_station_tuple = pick_s_phase(prefor_net_st_pick[0],
                                                  prefor_net_st_pick[1],
                                                  prefor_net_st_pick[2])
            else:
                print('The passed in arguments match neither for p-picking nor s-picking')
            req = comm.isend(pick_station_tuple, dest = 0)
            req.wait()
    else:
        try:
            inv = client.get_stations(starttime=UTCDateTime(year+'-01-01T00:00:00.000000Z'), endtime=UTCDateTime(year+'-12-31T23:59:59.000000Z'), latitude=-22.5, longitude=127.5, maxradius=45, level='channel')
            print(('len(inv)=>'+str(len(inv))))
        except:
            print(('FDSN client could not retrieve inventory for start year='+year+' and end year='+year))
            sys.exit(1)

        dedup(inv)
        for f in evtfiles:
            evts = read_events(f)
            if evts:
                count = 1
                for evt in evts:
                    if evt:
                        prefor = evt.preferred_origin() or evt.origins[0]
                        if prefor.depth >= 0 and prefor.time > UTCDateTime(year+'-01-01T00:00:00.000000Z') and prefor.time < UTCDateTime(year+'-12-31T23:59:59.000000Z'):
                            print(('Processing event => ' + str(evt)))
                            if parallel:
                                p_picks, s_picks, stations = pick_phases_parallel(evt, inv)
                            else:
                                p_picks, s_picks, stations = pick_phases(evt, inv)
                            evt_out = createEventObject(evt, p_picks, s_picks, stations)
                            if evt_out:
                                #evt_out.write(outdir+'/'+os.path.splitext(os.path.basename(f))[0]+'-'+str(count)+os.path.splitext(os.path.basename(f))[1], format='SC3ML')
                                evt_out.write(output_folder+'/'+os.path.splitext(os.path.basename(f))[0]+'-'+str(count)+os.path.splitext(os.path.basename(f))[1], format='SC3ML')
                    count += 1

if __name__ == '__main__':
    rank = None
    size = None
    comm = None
    client = None
    # the hostname/ipaddress of the Seiscomp3 machine on which the
    # waveforms to be picked are hosted underneath a FDSN web server
    #sc3_host_ip='13.211.209.88'
    #sc3_host='test-elb-521962885.ap-southeast-2.elb.amazonaws.com'
    #client = Client('http://'+sc3_host_ip+':8081')
    #client = Client('http://'+sc3_host)
    # the travel time velocity reference model object to be used
    # for fetching theoretical p and s phase arrivals
    model = TauPyModel(model='iasp91')
    # self explanatory variables. the p and s phase look-back and
    # look-ahead time windows for respective phases around the
    # theoretical onset time
    plookback=50
    plookahead=100
    slookback=100
    slookahead=200
    # the margin to allow for the p and s phase arrivals pick time
    pphase_search_margin=10
    sphase_search_margin=25

    # maximum allowed value for SNR
    snr_max=3000

    # some global variables (Seiscomp3 ingestion requirement) to have
    # unique resource identifiers for picks, origins and events that
    # will comprise the generated xml files
    pick_Count = 1
    origin_Count = 1
    event_Count = 1
    main()

