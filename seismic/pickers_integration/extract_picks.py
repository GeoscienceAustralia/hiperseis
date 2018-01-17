from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core import read_events, read_inventory
from obspy.taup import TauPyModel
from obspy.signal.trigger import classic_sta_lta, trigger_onset

sc3_host_ip='54.79.27.220'

client = Client('http://'+sc3_host_ip+':8081')

inventory = read_inventory('inventory/inventory.xml')

def s_phase_exists(evt, pick):
    exists=False
    for p in evt.picks:
        if p.waveform_id.network_code==pick.waveform_id.network_code and
           p.waveform_id.station_code==pick.waveform_id.station_code and
           p.waveform_id.channel_code==pick.waveform_id.channel_code and
           p.phase_hint=='S':
            exists=True
            break
    return exists

def process_event(event):
    model = TauPyModel(model='iasp91')

    if not event:
        return
    if len(event.picks) < 30:
        return
    pcount = 0
    for p in event.picks:
        pcount += 1
        net = p.waveform_id.network_code
        sta = p.waveform_id.station_code
        cha = p.waveform_id.channel_code
        loc = p.waveform_id.location_code
        t = p.time
        phasehint = p.phase_hint
        if phasehint=='P':
            sts = client.get_stations(network=net, station=sta)
            if not s_phase_exists(event, p) and sts and sts.networks[0] and sts.networks[0].stations[0]:
                prefor = event.prederred_origin() or event.origins[0]
                arrivals = model.get_travel_times_geo(prefor.depth/1000, prefor.latitude,
                             prefor.longitude, sts.networks[0].stations[0].latitude,
                             sts.networks[0].stations[0].longitude, phase_list=("S",)):
                if arrivals and len(arrivals):
                    st = client.get_waveforms(net, sta, '--', 'BH?', prefor.time+arrival[0].time-100, prefor.time+arrival[-1].time+150)
                    if st and len(st):
                        st = [s if !s[0].stats.channel.endswith('Z') for s in st]
                        for tr in st:
                            tr.trim(prefor.time+arrival[0].time-100, prefor.time+arrival[-1].time+150)
                            df = tr.stats.sampling_rate
                            cft = classic_sta_lta(tr.data, int(5*df), int(10*df))
                            trigger_onset(cft, trigOnThreshold, trigOffThreshold)

def process_trace(tr, t1, t2):
    tr.trim(t1, t2)
    tr.detrend('linear')
    tr.taper(0.05)
    tr.filter('bandpass', freqmin=2.0, freqmax=9.0)
    tr.plot()
