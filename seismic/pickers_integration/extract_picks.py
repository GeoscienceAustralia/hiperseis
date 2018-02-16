import os
import numpy as np
from obspy import UTCDateTime, read_events
from obspy.clients.fdsn import Client
from obspy.clients.iris import Client as IrisClient
from obspy.taup import TauPyModel
from obspy.taup.taup_geo import calc_dist
from obspy.signal.trigger import ar_pick, recursive_sta_lta, trigger_onset
from obspy.core.event import Pick, CreationInfo, WaveformStreamID, ResourceIdentifier, Arrival
from multiprocessing import Pool
import glob

sc3_host_ip='54.206.31.42'
#sc3_host_ip='test-elb-521962885.ap-southeast-2.elb.amazonaws.com'
client = Client('http://'+sc3_host_ip+':8081')
#client = Client('http://'+sc3_host_ip)
model = TauPyModel(model='iasp91')
plookback=50
plookahead=100
slookback=100
slookahead=200
pphase_search_margin=10
sphase_search_margin=25
global p_pickCount
p_pickCount = 1
global s_pickCount
s_pickCount = 1

def s_phase_exists(evt, pick):
    exists=False
    filtered_picks = [p for p in evt.picks if p.waveform_id.station_code == pick.waveform_id.station_code]
    for p in filtered_picks:
        if p.phase_hint=='S':
            exists=True
            break
    return exists

def clean_trace(tr, t1, t2, freqmin=1.0, freqmax=4.9):
    # add bandpass filtering sensitive to the distance d between source and receiver
    tr.trim(t1, t2)
    tr.detrend('linear')
    tr.taper(0.01)
    tr.filter('bandpass', freqmin=freqmin, freqmax=freqmax)

def createPickObject(net, sta, cha, time, backaz, phasehint, res=0.0, wt=1.0):
    global p_pickCount
    global s_pickCount
    if phasehint == 'P':
        count = p_pickCount
        p_pickCount += 1
    else:
        count = s_pickCount
        s_pickCount += 1
    return (Pick(resource_id=ResourceIdentifier(id='smi:niket.picker.ga.gov.au/pick/'+str(count)),
                 time=time,
                 waveform_id=WaveformStreamID(network_code=net, station_code=sta, channel_code=cha),
                 methodID=ResourceIdentifier('obspy/stalta/arpicker'),
                 backazimuth=backaz,
                 phase_hint=phasehint,
                 evaluation_mode='automatic',
                 creation_info=CreationInfo(author='niket',
                                 creation_time=UTCDateTime(),
                                 agency_id='niket-ga-picker')),
            res,
            wt)

def createArrivalObject(event, pick):
    pass

def get_backazimuth(stalat, stalon, evtlat, evtlon):
    client = IrisClient()
    return client.distaz(stalat=stalat, stalon=stalon, evtlat=evtlat, evtlon=evtlon)

def calc_distance(rec, src):
    return calc_dist(src.latitude, src.longitude, rec.latitude, rec.longitude, 6371, 0)

def find_best_bounds(cft, samp_rate):
    bestupper = 1.5
    bestlower = 0.75
    max_margin = bestupper - bestlower
    leasttrigs = len(trigger_onset(cft, 1.5, 0.75, max_len=(60*samp_rate), max_len_delete=True))
    for upper in np.linspace(1.5, 5, 15):
        for lower in [0.875, 0.75, 0.625, 0.5, 0.375]:
            t = trigger_onset(cft, upper, lower, max_len=(60*samp_rate), max_len_delete=True)
            if len(t) > 0 and (upper - lower) > max_margin and len(t) <= leasttrigs:
                leasttrigs = len(t)
                max_margin = upper - lower
                bestupper = upper
                bestlower = lower
    return bestupper, bestlower

def pick_p_s_phases(network, station, prefor, pickP=False, available_ptime=None):
    ppick = spick = None
    arrivals = model.get_travel_times_geo(prefor.depth/1000, prefor.latitude,
                 prefor.longitude, station.latitude, station.longitude, phase_list=("ttbasic",))
    parrivals = [ar for ar in arrivals if ar.phase.name.lower() == 'p']
    sarrivals = [ar for ar in arrivals if ar.phase.name.lower() == 's']
    parrivals.sort(key=lambda x: x.time)
    sarrivals.sort(key=lambda x: x.time)
    mean_sarrival = np.mean([sarr.time for sarr in sarrivals])
    mean_parrival = np.mean([parr.time for parr in parrivals])
    if len(sarrivals) > 0 and len(parrivals) > 0:
        try:
            theoretical_ptime = prefor.time+mean_parrival
            trim_starttime = prefor.time+parrivals[0].time-plookback
            if not pickP and available_ptime:
                theoretical_ptime = available_ptime
                trim_starttime = available_ptime - plookback
            trim_endtime = prefor.time+sarrivals[-1].time+slookahead
            stz = client.get_waveforms(network.code, station.code, '*', 'BHZ,SHZ,HHZ', trim_starttime, trim_endtime)
            stn = client.get_waveforms(network.code, station.code, '*', 'BHN,SHN,HHN,BH1,SH1', trim_starttime, trim_endtime)
            ste = client.get_waveforms(network.code, station.code, '*', 'BHE,SHE,HHZ,BH2,SH2', trim_starttime, trim_endtime)
            stz_raw = stz.copy()
            stn_raw = stn.copy()
            ste_raw = ste.copy()
        except Exception, e:
            print ('FDSN client could not retrieve waveform for net='+network.code+', sta='+station.code+' around time='+str(prefor.time+arrivals[0].time))
            print(str(e))
            return None, None
        if stz and stn and ste and stz[0].stats.sampling_rate == stn[0].stats.sampling_rate and \
                                   stn[0].stats.sampling_rate == ste[0].stats.sampling_rate:
            samp_rate, f1, f2, lta_p, sta_p, lta_s, sta_s, m_p, m_s, l_p, l_s = \
                stz[0].stats.sampling_rate, 1.0, 4.9, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2
            distance = calc_distance(prefor, station)
            for tr in [stz[0], stn[0], ste[0]]:
                clean_trace(tr, trim_starttime, trim_endtime)
            ptime, stime = ar_pick(stz[0].data, stn[0].data, ste[0].data, samp_rate, f1, f2,
                                   lta_p, sta_p, lta_s, sta_s, m_p, m_s, l_p, l_s)
            if ptime > 0 and stime > 0 and stime > ptime and abs(trim_starttime + ptime - theoretical_ptime) < pphase_search_margin and \
                abs(prefor.time + mean_sarrival - trim_starttime - stime) < sphase_search_margin:
                az = get_backazimuth(stalat=station.latitude, stalon=station.longitude, evtlat=prefor.latitude, evtlon=prefor.longitude)
                p_res = trim_starttime + ptime - theoretical_ptime
                s_res = trim_starttime + stime - prefor.time - mean_sarrival
                ppick = createPickObject(network.code, station.code, stz[0].stats.channel, trim_starttime+ptime, az['backazimuth'] if az else None, 'P', p_res)
                spick = createPickObject(network.code, station.code, stn[0].stats.channel, trim_starttime+stime, az['backazimuth'] if az else None, 'S', s_res)
            elif (not pickP and available_ptime) or \
                 (pickP and ptime > 0 and abs(trim_starttime + ptime - available_ptime) < pphase_search_margin):
                if not pickP and available_ptime:
                    ptime = available_ptime - trim_starttime
                az = get_backazimuth(stalat=station.latitude, stalon=station.longitude, evtlat=prefor.latitude, evtlon=prefor.longitude)
                p_res = trim_starttime + ptime - theoretical_ptime
                ppick = createPickObject(network.code, station.code, stz[0].stats.channel, trim_starttime+ptime, az['backazimuth'] if az else None, 'P', p_res)
                trim_starttime = prefor.time+sarrivals[0].time-slookback
                trim_endtime = prefor.time+sarrivals[-1].time+slookahead
                stz_for_s = client.get_waveforms(network.code, station.code, '*', 'BHZ,SHZ,HHZ', trim_starttime, trim_endtime)
                stn_for_s = client.get_waveforms(network.code, station.code, '*', 'BHN,SHN,HHN,BH1,SH1', trim_starttime, trim_endtime)
                ste_for_s = client.get_waveforms(network.code, station.code, '*', 'BHE,SHE,HHZ,BH2,SH2', trim_starttime, trim_endtime)
                trigs= []
                if trim_starttime < available_ptime:
                    trim_starttime = available_ptime + 10
                for tr in [stz_for_s[0], stn_for_s[0], ste_for_s[0]]:
                    # replace steps in this loop with a function that runs on 4 cores in parallel
                    # with different filter windows and returns best trigger for this trace
                    for band in [(0.5, 3), (0.5, 2.5), (1, 3), (0.5, 2), (1, 2.5), (1.5, 3), (0.3, 2.0), (0.3, 1), (0.3, 0.7)]:
                        tr_copy = tr.copy()
                        clean_trace(tr_copy, trim_starttime, trim_endtime, band[0], band[1])
                        cft = recursive_sta_lta(tr_copy.data, int(5*samp_rate), int(20*samp_rate))
                        upper, lower = find_best_bounds(cft, tr_copy.stats.sampling_rate)
                        trigs.extend([(onset, tr_copy.stats.channel, upper-lower, band) for onset in trigger_onset(cft, upper, lower, max_len=(60*tr_copy.stats.sampling_rate), max_len_delete=True)])
                trigs = [trig for trig in trigs if abs(prefor.time + mean_sarrival - trim_starttime - (trig[0][0]/samp_rate)) < sphase_search_margin]
                if len(trigs) > 0:
                    mintrigdiff = abs(prefor.time + mean_sarrival - trim_starttime - trigs[0][0][0])
                    besttrig = trigs[0][0][0]
                    best_cha = trigs[0][1]
                    best_trig_margin = trigs[0][2]
                    best_band = trigs[0][3]
                    for trig in trigs:
                        if abs(prefor.time + mean_sarrival - trim_starttime - trig[0][0]) < mintrigdiff and \
                            trig[2] > best_trig_margin:
                            mintrigdiff = abs(prefor.time + mean_sarrival - trim_starttime - trig[0][0])
                            besttrig = trig[0][0]
                            best_cha = trig[1]
                            best_trig_margin = trig[2]
                            best_band = trig[3]
                    if best_trig_margin > 2:
                        s_res = trim_starttime + (besttrig/samp_rate) - prefor.time - mean_sarrival
                        spick = createPickObject(network.code, station.code, best_cha, trim_starttime+(besttrig/samp_rate), az['backazimuth'] if az else None, 'S', s_res)
            elif pickP and ptime <= 0:
                # this needs to be thought through and implemented, dummy behavior for now
                ppick = None
                spick = None
            if pickP:
                return ppick, spick
            else:
                return None, spick
    return None, None

def find_trigger_parallel(trace, starttime, endtime):
    pass

def pick_new_phases(event):
    prefor = event.preferred_origin() or event.origins[0]
    inv = client.get_stations(latitude=prefor.latitude, longitude=prefor.longitude, minradius=0, maxradius=90)
    picklist = []
    for n in inv.networks:
        for s in n.stations:
            if s.code not in [p.waveform_id.station_code for p in event.picks]:
                for p in pick_p_s_phases(n, s, prefor, pickP=True):
                    if p:
                        picklist.append(p)
    return picklist

def pick_sphase(event):
    if not event:
        return
    if len(event.picks) < 30:
        return
    s_phases = []
    for p in event.picks:
        net = p.waveform_id.network_code
        sta = p.waveform_id.station_code
        cha = p.waveform_id.channel_code
        loc = p.waveform_id.location_code
        t = p.time
        if p.phase_hint=='P':
            try:
                sts = client.get_stations(network=net, station=sta)
            except Exception, e:
                print ('FDSN client could not retrieve station details for net='+net+', sta='+sta)
                print(str(e))
                continue
            if not s_phase_exists(event, p):
                prefor = event.preferred_origin() or event.origins[0]
                s_phases.extend([p for p in pick_p_s_phases(sts.networks[0], sts.networks[0].stations[0], prefor, available_ptime=p.time) if p])
    return s_phases

def add_picks_to_event(event, picks):
    for pick in picks:
        net = pick[0].waveform_id.network_code
        sta = pick[0].waveform_id.station_code
        try:
            sts = client.get_stations(network=net, station=sta)
        except Exception, e:
            print ('FDSN client could not retrieve station details for net='+net+', sta='+sta)
            print(str(e))
            continue
        stalat = sts.networks[0].stations[0].latitude
        stalon = sts.networks[0].stations[0].longitude
        prefor = evts[0].preferred_origin() or evts[0].origins[0]
        evtlat = event.preferred_origin().latitude
        evtlon = event.preferred_origin().longitude
        az = get_backazimuth(stalat=stalat, stalon=stalon, evtlat=evtlat, evtlon=evtlon)
        arr = Arrival(resource_id=ResourceIdentifier(id=pick[0].resource_id.id+"#"),
                      pick_id=pick[0].resource_id,
                      phase='S',
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

evtfiles = glob.glob('event_xmls_sc3ml/976*.xml')
outdir = 'event_xmls_sc3ml_out'
for f in evtfiles:
    evts = read_events(f)
    if evts and evts[0]:
        prefor = evts[0].preferred_origin() or evts[0].origins[0]
        if prefor.time > UTCDateTime('2015-03-01T00:00:00.000000Z') and prefor.time < UTCDateTime('2015-04-01T00:00:00.000000Z'):
            s_picks = pick_sphase(evts[0])
#            new_picks = pick_new_phases(evts[0])
            if s_picks and len(s_picks) > 0:
                add_picks_to_event(evts[0], s_picks)
                evts[0].write(os.path.join(outdir, os.path.basename(f)), format='SC3ML')
