from obspy import UTCDateTime, read_events
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.signal.trigger import ar_pick, classic_sta_lta, trigger_onset
import glob

sc3_host_ip='54.79.27.220'
trigOnThreshold=1.5
trigOffThreshold=0.5
client = Client('http://'+sc3_host_ip+':8081')

def s_phase_exists(evt, pick):
    exists=False
    filtered_picks = [p for p in evt.picks if p.waveform_id.station_code == pick.waveform_id.station_code]
    for p in filtered_picks:
        if p.phase_hint=='S':
            exists=True
            break
    return exists

def clean_trace(tr, t1, t2):
    tr.trim(t1, t2)
    tr.detrend('linear')
    tr.taper(0.05)
    tr.filter('bandpass', freqmin=2.0, freqmax=4.9)

def pick_sphase(event):
    model = TauPyModel(model='iasp91')

    if not event:
        return
    if len(event.picks) < 30:
        return
    for p in event.picks:
        net = p.waveform_id.network_code
        sta = p.waveform_id.station_code
        cha = p.waveform_id.channel_code
        loc = p.waveform_id.location_code
        t = p.time
        phasehint = p.phase_hint
        if phasehint=='P':
            try:
                sts = client.get_stations(network=net, station=sta)
                if not s_phase_exists(event, p):
                    prefor = event.preferred_origin() or event.origins[0]
                    arrivals = model.get_travel_times_geo(prefor.depth/1000, prefor.latitude,
                                 prefor.longitude, sts.networks[0].stations[0].latitude,
                                 sts.networks[0].stations[0].longitude, phase_list=("ttbasic",))
                    sarrivals = [ar for ar in arrivals if ar.phase.name.lower() == 's']
                    arrivals.sort(key=lambda x: x.time)
                    if len(sarrivals) > 0:
                        stz = client.get_waveforms(net, sta, '--', 'BHZ', prefor.time+arrivals[0].time-100, prefor.time+arrivals[-1].time+150)
                        stn = client.get_waveforms(net, sta, '--', 'BHN', prefor.time+arrivals[0].time-100, prefor.time+arrivals[-1].time+150)
                        ste = client.get_waveforms(net, sta, '--', 'BHE', prefor.time+arrivals[0].time-100, prefor.time+arrivals[-1].time+150)
                        if stz and stn and ste and stz[0].stats.sampling_rate == stn[0].stats.sampling_rate and \
                                                   stn[0].stats.sampling_rate == ste[0].stats.sampling_rate:
                            samp_rate, f1, f2, lta_p, sta_p, lta_s, sta_s, m_p, m_s, l_p, l_s = \
                                stz[0].stats.sampling_rate, 2.0, 4.9, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2
                            for st in [stz, stn, ste]:
                                clean_trace(st[0], prefor.time+arrivals[0].time-100, prefor.time+arrivals[-1].time+150)
                            ptime, stime = ar_pick(stz[0].data, stn[0].data, ste[0].data, samp_rate, f1, f2,
                                                   lta_p, sta_p, lta_s, sta_s, m_p, m_s, l_p, l_s)
                            if stime > ptime:
                                print ('ar_picker discovered ptime =>'+ptime+' and stime =>'+stime)
                            elif ptime > 95:
                                for tr in [stn[0], ste[0]]:
                                    cft = classic_sta_lta(tr.data, int(5*samp_rate), int(10*samp_rate))
                                    trigs = trigger_onset(cft, 1.5, 0.75)
                                    print [(tr.stats.channel, tr.stats.sampling_rate, abs(sarrivals[0].time*samp_rate-trig[0])) for trig in trigs if abs(sarrivals[0].time*samp_rate-trig[0]) < 1000]
            except:
                continue


evtfiles = glob.glob('event_sc3mls/*.xml')
for f in evtfiles:
    evts = read_events(f)
    if evts and evts[0]:
        prefor = evts[0].preferred_origin() or evts[0].origins[0]
        if prefor.time > UTCDateTime('2015-03-01T00:00:00.000000Z') and prefor.time < UTCDateTime('2015-04-01T00:00:00.000000Z'):
            pick_sphase(evts[0])
