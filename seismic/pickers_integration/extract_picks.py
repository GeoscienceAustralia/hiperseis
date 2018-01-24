import numpy as np
from obspy import UTCDateTime, read_events
from obspy.clients.fdsn import Client
from obspy.clients.iris import Client as IrisClient
from obspy.taup import TauPyModel
from obspy.taup.taup_geo import calc_dist
from obspy.signal.trigger import ar_pick, classic_sta_lta, trigger_onset
from obspy.core.event import Pick, CreationInfo, WaveformStreamID, ResourceIdentifier
import glob

sc3_host_ip='54.79.27.220'
client = Client('http://'+sc3_host_ip+':8081')
model = TauPyModel(model='iasp91')
lookback=50
lookahead=50
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

def clean_trace(tr, t1, t2, d):
    # add bandpass filtering sensitive to the distance d between source and receiver
    tr.trim(t1, t2)
    tr.detrend('linear')
    tr.taper(0)
    tr.filter('bandpass', freqmin=1.0, freqmax=4.9)

def createPickObject(net, sta, cha, time, backaz, phasehint):
    global p_pickCount
    global s_pickCount
    if phasehint == 'P':
        count = p_pickCount
        p_pickCount += 1
    else:
        count = s_pickCount
        s_pickCount += 1
    return Pick(resource_id=ResourceIdentifier(id='smi:niket.picker.ga.gov.au/pick/'+str(count)),
                time=time,
                waveform_id=WaveformStreamID(network_code=net, station_code=sta, channel_code=cha),
                methodID=ResourceIdentifier('obspy/stalta/arpicker'),
                backazimuth=backaz,
                phase_hint=phasehint,
                evaluation_mode='automatic',
                creation_info=CreationInfo(author='niket',
                                creation_time=UTCDateTime(),
                                agency_id='niket-ga-picker'))

def get_backazimuth(stalat, stalon, evtlat, evtlon):
    client = IrisClient()
    return client.distaz(stalat=stalat, stalon=stalon, evtlat=evtlat, evtlon=evtlon)

def calc_distance(rec, src):
    return calc_dist(src.latitude, src.longitude, rec.latitude, rec.longitude, 6371, 0)

def find_best_bounds(cft):
    bestupper = 1.5
    bestlower = 0.75
    leasttrigs = len(trigger_onset(cft, 1.5, 0.75))
    for upper in np.linspace(1.5, 2.5, 0.25):
        for lower in [1.0, 0.75, 0.5]:
            t = trigger_onset(cft, upper, lower)
            if len(t) < leasttrigs and len(t) > 0:
                leasttrigs = len(t)
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
    if len(sarrivals) > 0 and len(parrivals) > 0:
        try:
            stz = client.get_waveforms(network.code, station.code, '*', 'BHZ,SHZ,HHZ', prefor.time+parrivals[0].time-lookback, prefor.time+sarrivals[-1].time+lookahead)
            stn = client.get_waveforms(network.code, station.code, '*', 'BHN,SHN,HHN,BH1,SH1', prefor.time+parrivals[0].time-lookback, prefor.time+sarrivals[-1].time+lookahead)
            ste = client.get_waveforms(network.code, station.code, '*', 'BHE,SHE,HHZ,BH2,SH2', prefor.time+parrivals[0].time-lookback, prefor.time+sarrivals[-1].time+lookahead)
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
            for st in [stz, stn, ste]:
                distance = calc_distance(prefor, station)
                clean_trace(st[0], prefor.time+parrivals[0].time-lookback, prefor.time+sarrivals[-1].time+lookahead, distance)
            ptime, stime = ar_pick(stz[0].data, stn[0].data, ste[0].data, samp_rate, f1, f2,
                                   lta_p, sta_p, lta_s, sta_s, m_p, m_s, l_p, l_s)
            if ptime > 0 and stime > 0 and stime > ptime and abs(lookback - ptime) < 10 and \
                abs(sarrivals[0].time - parrivals[0].time + lookback - stime) < 20:
                az = get_backazimuth(stalat=station.latitude, stalon=station.longitude, evtlat=prefor.latitude, evtlon=prefor.longitude)
                ppick = createPickObject(network.code, station.code, stz[0].stats.channel, prefor.time+arrivals[0].time-lookback+ptime, az['backazimuth'] if az else None, 'P')
                spick = createPickObject(network.code, station.code, stn[0].stats.channel, prefor.time+arrivals[0].time-lookback+stime, az['backazimuth'] if az else None, 'S')
            elif (pickP and ptime > 0 and abs(parrivals[0].time - arrivals[0].time + lookback - ptime) < 10) or \
                 (not pickP and available_ptime):
                if not pickP and available_ptime:
                    ptime = available_ptime
                az = get_backazimuth(stalat=station.latitude, stalon=station.longitude, evtlat=prefor.latitude, evtlon=prefor.longitude)
                ppick = createPickObject(network.code, station.code, stz[0].stats.channel, prefor.time+arrivals[0].time-lookback+ptime, az['backazimuth'] if az else None, 'P')
                trigs= []
                for tr in [stn[0], ste[0]]:
                    tr.trim(starttime=(prefor.time+parrivals[0].time-lookback+ptime))
                    tr.taper(0)
                    cft = classic_sta_lta(tr.data, int(5*samp_rate), int(20*samp_rate))
                    upper, lower = find_best_bounds(cft)
                    trigs.extend(trigger_onset(cft, upper, lower))
                trigs = [trig for trig in trigs if abs(sarrivals[0].time-arrivals[0].time+lookback-(trig[0]/samp_rate)) < 20]
                if len(trigs) > 0:
                    mintrigdiff = abs(sarrivals[0].time*samp_rate-trigs[0][0])
                    besttrig = trigs[0][0]
                    for trig in trigs:
                        if abs(sarrivals[0].time*samp_rate-trig[0]) < mintrigdiff:
                            mintrigdiff = abs(sarrivals[0].time*samp_rate-trig[0])
                            besttrig = trig[0]
                    spick = createPickObject(network.code, station.code, stn[0].stats.channel, prefor.time+arrivals[0].time-lookback+(besttrig/samp_rate), az['backazimuth'] if az else None, 'S')
            elif pickP and ptime <= 0:
                # this needs to be thought through and implemented, dummy behavior for now
                ppick = None
                spick = None
            if pickP:
                return ppick, spick
            else:
                return None, spick
    return None

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

evtfiles = glob.glob('event_sc3mls/*.xml')
for f in evtfiles:
    evts = read_events(f)
    if evts and evts[0]:
        prefor = evts[0].preferred_origin() or evts[0].origins[0]
        if prefor.time > UTCDateTime('2015-03-01T00:00:00.000000Z') and prefor.time < UTCDateTime('2015-04-01T00:00:00.000000Z'):
#            s_picks = pick_sphase(evts[0])
            new_picks = pick_new_phases(evts[0])

