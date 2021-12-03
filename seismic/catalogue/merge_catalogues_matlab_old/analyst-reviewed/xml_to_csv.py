from obspy import read_events
import glob

evt_files = glob.glob('*.xml')
csv_out_file = 'analyst_events.csv'
with open(csv_out_file, 'w') as csv_out:
    for evt_file in evt_files:
        evts = read_events(evt_file)
        for evt in evts:
            prefor = evt.preferred_origin() if evt.preferred_origin() else evt.origins[0]
            magPresent = True
            if evt.preferred_magnitude() or len(evt.magnitudes) > 0:
                prefmag = evt.preferred_magnitude() or evt.magnitudes[0]
            else:
                magPresent = False
            csv_out.write(str(prefor.latitude) + ',' + str(prefor.longitude) + ',' + str(prefor.depth/1000) + ',' + (str(prefmag.mag) if magPresent else '') + ',' + str(prefor.time) + '\n')

