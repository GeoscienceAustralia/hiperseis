import sys
from obspy import read_inventory, Stream
from rf import read_rf
from eqcorrscan.utils.clustering import cluster
from eqcorrscan.utils.stacking import linstack, PWS_stack
import pandas as pd

stream = read_rf('data/7X-rf_profile_rfs-cleaned.h5', 'H5')
inv = read_inventory('data/7X-inventory.xml')

def convert_ms_to_ascii(msfile_path, outfile_folder):
    msf = read(msfile_path)
    df = pd.DataFrame(columns=['offset_from_onset', 'value'])
    offset=-5.0
    for i in range(len(msf[0].data)):
        df.loc[i] = [offset, msf[0].data[i]]
        offset += msf[0].stats.delta
    df.to_csv(outfile_folder+'/'+os.path.splitext(os.path.basename(msfile_path))[0]+'.dat', sep=' ', index=False, header=False)

for st in inv.networks[0].stations:
    stream_list = [(Stream(traces=[tr, ]), i) for i, tr in enumerate(stream.select(component='Q', station=st.code).trim2(-5, 20, 'onset'))]
    try:
        groups = cluster(template_list=stream_list, show=False, corr_thresh=0.3, cores=4)
        if len(groups) < len(stream_list):
            group_max = groups[0]
            for g in groups:
                if len(g) > len(group_max):
                    group_max = g
            group_streams = [st_tuple[0] for st_tuple in group_max]
            stack = linstack(streams=group_streams)
            stack.write('data/'+stack[0].stats.network+'.'+stack[0].stats.station+'..'+stack[0].stats.channel+'-'+str(len(group_max))+'-linstack.ms', format='MSEED')
            stack = PWS_stack(streams=group_streams)
            stack.write('data/'+stack[0].stats.network+'.'+stack[0].stats.station+'..'+stack[0].stats.channel+'-'+str(len(group_max))+'-pwsstack.ms', format='MSEED')
    except:
        print "Error: ", sys.exc_info()[0]
