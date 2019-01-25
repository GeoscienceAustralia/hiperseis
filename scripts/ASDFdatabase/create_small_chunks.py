import os
import numpy as np
from obspy import read_inventory, Stream
from obspy import UTCDateTime as UTC
from pyasdf import ASDFDataSet
import datetime

base_dir = '/g/data/ha3/Passive/_ANU/7B(1993-1995)'
asdf_file = os.path.join(base_dir, 'ASDF', '7B(1993-1995).h5')
out_dir = os.path.join(base_dir, 'small_mseed_DATA')

inv = read_inventory(os.path.join(base_dir, '7B.xml'))
asdf = ASDFDataSet(asdf_file, mode='r')

def create_chunk(trace, st_time, end_time, sta):
    trace.trim(starttime=st_time, endtime=end_time)
    st_out = Stream(traces=[trace, ])
    dest_dir = os.path.join(out_dir, str(trace.stats.starttime.timetuple().tm_year), str(tr.stats.starttime.timetuple().tm_yday))
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    st_out.write(os.path.join(dest_dir, sta.code+'_'+trace.stats.channel+'_'+str(st_time)+'_'+str(end_time)+'.ms'), format='MSEED')

for sta in inv.networks[0].stations:
    if asdf.waveforms.__contains__(inv.networks[0].code+'.'+sta.code):
        for i in asdf.waveforms[inv.networks[0].code+'.'+sta.code].list():
            if i.endswith('raw_recording'):
                start_time = UTC(i.split("__")[1])
                st = asdf.waveforms[inv.networks[0].code+'.'+sta.code][i]
                medn = np.median(st[0].data)
                while (abs(st[0].data[np.argmax(st[0].data)]) > 1e8 or abs(st[0].data[np.argmin(st[0].data)]) > 1e8):
                    if abs(st[0].data[np.argmax(st[0].data)]) > 1e8:
                        st[0].data[np.argmax(st[0].data)] = abs(medn) if st[0].data[np.argmax(st[0].data)] > 0 else -abs(medn)
                    if abs(st[0].data[np.argmin(st[0].data)]) > 1e8:
                        st[0].data[np.argmin(st[0].data)] = abs(medn) if st[0].data[np.argmin(st[0].data)] > 0 else -abs(medn)
                while (start_time+86400<UTC(i.split("__")[2])):
                    tr = st[0].copy()
                    create_chunk(tr, start_time, start_time+86400, sta)
                    start_time += 86400
                if start_time < UTC(i.split("__")[2]):
                    tr=st[0].copy()
                    create_chunk(tr, start_time, UTC(i.split("__")[2]), sta)

