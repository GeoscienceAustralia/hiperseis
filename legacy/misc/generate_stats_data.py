import glob
from obspy import read_events
from obspy.geodetics.base import degrees2kilometers
import csv

evtfiles = glob.glob('/home/ubuntu/isc-out-final/2010_iloc_phase_defined/*.xml')
with open('/home/ubuntu/p_stats_in_2010-m.csv', 'w') as p_out, open('/home/ubuntu/s_stats_in_2010-m.csv', 'w') as s_out:
    p_writer = csv.writer(p_out)
    s_writer = csv.writer(s_out)
    for f in evtfiles:
        evts = read_events(f)
        if evts and evts[0]:
            evt = evts[0]
            for i in range(len(evt.picks)):
                for j in range(len(evt.preferred_origin().arrivals)):
                    if evt.preferred_origin().arrivals[j].pick_id == evt.picks[i].resource_id.id:
                        dist_deg = evt.preferred_origin().arrivals[j].distance
                        dist_km = degrees2kilometers(dist_deg)
                        tt = evt.picks[i].time - evt.preferred_origin().time
                        res = evt.preferred_origin().arrivals[j].time_residual
                        if evt.preferred_origin().arrivals[j].phase == 'P':
                            p_writer.writerow([dist_km/tt, dist_deg, res])
                        elif evt.preferred_origin().arrivals[j].phase == 'S':
                            s_writer.writerow([dist_km/tt, dist_deg, res])
                        else:
                            print "New phase uncovered: " +str(evt.preferred_origin().arrivals[j].phase)

