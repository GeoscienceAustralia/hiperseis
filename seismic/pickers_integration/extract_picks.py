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

sc3_host_ip='13.236.167.35'
client = Client('http://'+sc3_host_ip+':8081')
model = TauPyModel(model='iasp91')
plookback=50
plookahead=100
slookback=100
slookahead=200
pphase_search_margin=10
sphase_search_margin=25
global pick_Count
pick_Count = 1
global event_Count
event_Count = 1
global origin_Count
origin_Count = 1

def s_phase_exists(evt, pick):
    exists=False
    filtered_picks = [p for p in evt.picks if p.waveform_id.station_code == pick.waveform_id.station_code]
    for p in filtered_picks:
        if p.phase_hint=='S':
            exists=True
            break
    return exists

def calc_distance(rec, src):
    return calc_dist(src.latitude, src.longitude, rec.latitude, rec.longitude, 6371, 0)

def get_backazimuth(stalat, stalon, evtlat, evtlon):
    client = IrisClient()
    return client.distaz(stalat=stalat, stalon=stalon, evtlat=evtlat, evtlon=evtlon)

def createPickObject(net, sta, cha, time, backaz, phasehint, res=0.0, wt=1.0, comments_data=None):
    global pick_Count
    count = pick_Count
    pick_Count += 1
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
    return (Pick(resource_id=ResourceIdentifier(id='smi:bilby2008.picker.ga.gov.au/pick/'+str(count)),
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

def clean_trace(tr, t1, t2, freqmin=1.0, freqmax=4.9):
    # add bandpass filtering sensitive to the distance d between source and receiver
    tr.trim(t1, t2)
    tr.detrend('linear')
    tr.taper(0.01)
    # added try catch after script aborted due to the error below:
    # ValueError: Selected corner frequency is above Nyquist.
    try:
        tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax)
    except:
        pass

def find_best_bounds(cft, samp_rate):
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
    return_pick = None
    snr = 0.0
    #p_bands = [(1, 6), (0.3, 2.3), (0.5, 2.5), (0.8, 2.8), (1, 3), (2, 4), (3, 5), (4, 6), (0.3, 1.3), (0.5, 1.5), (0.8, 1.8), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1.5, 2.5), (2.5, 3.5), (3.5, 4.5), (4.5, 5.5)]
    #s_bands = [(0.5, 2), (1, 2), (0.2, 1.5), (0.3, 2.0), (0.3, 1), (0.3, 0.7), (0.05, 0.3), (0.05, 1)]
    p_bands = [(0.5, 2.0), (0.8, 2.3)]
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
    print ('station=> ' + str(station.code) + ' distance =>'+str(distance))
    traces = []
    if len(target_arrivals) > 0:
        try:
            trim_starttime = prefor.time+target_arrivals[0].time-lookback
            trim_endtime = prefor.time+target_arrivals[-1].time+lookahead
            stz = client.get_waveforms(network.code, station.code, '*', 'BHZ,SHZ,HHZ', trim_starttime, trim_endtime)
            stz_raw = stz.copy()
            traces.append(stz[0])
            mult_const = 10
            if phase=='S':
                stn = client.get_waveforms(network.code, station.code, '*', 'BHN,SHN,HHN,BH1,SH1', trim_starttime, trim_endtime)
                ste = client.get_waveforms(network.code, station.code, '*', 'BHE,SHE,HHZ,BH2,SH2', trim_starttime, trim_endtime)
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
        except Exception, e:
            print ('FDSN client could not retrieve waveform for sta='+station.code+' around time='+str(prefor.time+mean_target_arrival))
            print(str(e))
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
                print('band[0] => ' +str(band[0]) + ' band[1] => ' + str(band[1]))
                cft = recursive_sta_lta(tr_copy.data, int(5*samp_rate), int(mult_const*samp_rate))
                upper, lower = find_best_bounds(cft, tr_copy.stats.sampling_rate)
                trigs.extend([(onset, tr_copy.stats.channel, upper-lower, band, upper) for onset in trigger_onset(cft, upper, lower, max_len=(60*tr_copy.stats.sampling_rate), max_len_delete=True)])
            search_margin = pphase_search_margin if phase == 'P' else sphase_search_margin
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
                    aic, aic_deriv, aic_deriv_cleaned = calc_aic(tr_copy)
                    aic[0:int(10*samp_rate)]=aic[np.argmax(aic)]
                    aic[-int(10*samp_rate):-1]=aic[np.argmax(aic)]
                    pick_index = np.argmin(aic)
                    aic_deriv_cleaned[0:int(10*samp_rate)]=0
                    aic_deriv_cleaned[-int(10*samp_rate):-1]=0
                    pick_index_deriv = np.argmax(aic_deriv_cleaned)
                    # internal plotting function to access the variables
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
                                aic, aic_deriv, aic_deriv_cleaned = calc_aic(tr_copy)
                                aic[0:int(10*samp_rate)]=aic[np.argmax(aic)]
                                aic[-int(10*samp_rate):-1]=aic[np.argmax(aic)]
                                pick_index = np.argmin(aic)
                                aic_deriv_cleaned[0:int(10*samp_rate)]=0
                                aic_deriv_cleaned[-int(10*samp_rate):-1]=0
                                pick_index_deriv = np.argmax(aic_deriv_cleaned)
                                index = pick_index_deriv if deriv else pick_index
                                snr_plot = calc_snr(tr_copy, trim_starttime + (index/samp_rate))
                                axes[ind+1, col_index].plot(range(int(index)), tr_copy.data[:int(index)], color='blue')
                                axes[ind+1, col_index].plot(range(int(index), int(index)+20), tr_copy.data[int(index):int(index)+20], color='red')
                                axes[ind+1, col_index].plot(range(int(index)+20, tr_copy.stats.npts), tr_copy.data[int(index)+20:tr_copy.stats.npts], color='blue')
                                axes[ind+1, col_index].text(int(tr_copy.stats.npts*0.5), int(max(tr_copy.data)*0.7), network.code+' '+station.code+' '+comp.upper()+' ('+str(band[0])+'Hz - '+str(band[1])+'Hz) SNR='+str(int(snr_plot))+' pick_index = '+str(index))
                        print('Break here while debugging and run plt.show()')

                    theoretical_trig = int((prefor.time + mean_target_arrival - trim_starttime)*samp_rate)
                    if theoretical_trig > 0 and abs(pick_index_deriv - theoretical_trig)/samp_rate < search_margin/2:
                        snr = calc_snr(tr_copy, trim_starttime + (pick_index_deriv/samp_rate))
                        #plot_onsets(deriv=True)
                        res = trim_starttime + (pick_index_deriv/samp_rate) - prefor.time - mean_target_arrival
                        comments_data=comments_data+(snr,)
                        return_pick = createPickObject(network.code, station.code, best_cha, trim_starttime+(pick_index_deriv/samp_rate), az['backazimuth'] if az else None, phase, res, comments_data=comments_data)
                        print(phase+'-pick added')
                    else:
                        print('theo_trig => '+str(theoretical_trig)+'pick_index_deriv => ' + str(pick_index_deriv) + ' besttrig => ' + str(besttrig) + '. Investigate waveforms!')

    if snr > 0 and snr < 1.5:
        print('The calculated SNR => ' + str(snr) + '. Discarding this pick!')
        return None
    return return_pick

def calc_snr(tr, time):
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
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def calc_aic(tr):
    npts = tr.stats.npts
    data = tr.data
    margin = int(tr.stats.sampling_rate*5)
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
            aic_deriv_cleaned.append(0)
        else:
            if aic[i - margin] < aic[i] < aic[i + margin] or \
                aic[i + margin] < aic[i] < aic[i - margin]:
                aic_deriv_cleaned.append(0)
            else:
                aic_deriv_cleaned.append(aic[i+1] - aic[i])

    return np.array(aic), np.array(aic_deriv), np.array(aic_deriv_cleaned)

def pick_phases(event, inventory=None):
    if not event:
        return None, None
    prefor = event.preferred_origin() or event.origins[0]

    p_phases = []
    p_stations = []
    p_pick = None
    if inventory:
        for st in inventory.networks[0].stations:
            p_stations.append((inventory.networks[0], st))
            p_pick = pick_phase(inventory.networks[0], st, prefor)
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
    for ppick in p_phases:
        if inventory:
            st = [s for s in inventory.networks[0].stations if s.code == ppick[0].waveform_id.station_code][0]
            s_pick = pick_phase(inventory.networks[0], st, prefor, phase='S', p_Pick=ppick[0])
        else:
            s_stn = [s for s in p_stations if s[1].code == ppick[0].waveform_id.station_code][0]
            s_pick = pick_phase(s_stn[0], s_stn[1], prefor, phase='S', p_Pick=ppick[0])
        if s_pick:
            s_phases.append(s_pick)
    return p_phases, s_phases, p_stations

def add_picks_to_event(event, picks, stations):
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

def createEventObject(evt, p_picks, s_picks, stations):
    global event_Count
    global origin_Count

    if len(p_picks) == 0:
        return None

    ev_count = event_Count
    or_count = origin_Count
    event_Count += 1
    origin_Count += 1
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

    origin = Origin(resource_id=ResourceIdentifier(id='smi:bilby2008.ga.gov.au/origin/'+str(or_count)),
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
        magnitude = Magnitude(resource_id=ResourceIdentifier(id='smi:bilby2008.ga.gov.au/origin/'+str(or_count)+'#netMag.Mb'),
                              mag=prefmag.mag,
                              magnitude_type=prefmag.magnitude_type,
                              origin_id=ResourceIdentifier(id='smi:bilby2008.ga.gov.au/origin/'+str(or_count)),
                              creation_info=creation_info)

    event = Event(resource_id=ResourceIdentifier(id='smi:bilby2008.picker.ga.gov.au/event/'+str(ev_count)),
                  creation_info=creation_info,
                  event_type='earthquake')

    event.origins = [origin, ]
    if magPresent:
        event.magnitudes = [magnitude, ]
        event.preferred_magnitude_id = magnitude.resource_id

    event.preferred_origin_id = origin.resource_id

    add_picks_to_event(event, p_picks+s_picks, stations)

    return event

evtfiles = glob.glob('/home/ubuntu/engdahl/*.xml')
outdir = '/home/ubuntu/bilby-out-final/2008'
#inv = read_inventory('/home/ubuntu/7W_dummy_resp.xml')
for f in evtfiles:
    evts = read_events(f)
    if evts:
        count = 1
        for evt in evts:
            if evt:
                prefor = evt.preferred_origin() or evt.origins[0]
                if prefor.depth >= 0 and prefor.time > UTCDateTime('2008-08-01T00:00:00.000000Z') and prefor.time < UTCDateTime('2008-12-31T23:59:59.000000Z'):
                    print('Processing event => ' + str(evt))
                    #p_picks, s_picks, stations = pick_phases(evt, inv)
                    p_picks, s_picks, stations = pick_phases(evt)
                    evt_out = createEventObject(evt, p_picks, s_picks, stations)
                    if evt_out:
                        evt_out.write(outdir+'/'+os.path.splitext(os.path.basename(f))[0]+'-'+str(count)+os.path.splitext(os.path.basename(f))[1], format='SC3ML')
            count += 1


