import os
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
    convert_ms_to_ascii(msf, outfile_folder)

def convert_stream_to_ascii(stream, outfile_folder, group_len):
    df = pd.DataFrame(columns=['offset_from_onset', 'value'])
    offset=-5.0
    for i in range(len(stream[0].data)):
        df.loc[i] = [offset, stream[0].data[i]]
        offset += stream[0].stats.delta
    df.to_csv(outfile_folder+'/'+stream[0].stats.network+'.'+stream[0].stats.station+'..'+stream[0].stats.channel+'.'+str(group_len)+'.dat', sep=' ', index=False, header=False)

for st in inv.networks[0].stations:
    stream_list_l = [(Stream(traces=[tr, ]), i) for i, tr in enumerate(stream.select(component='L', station=st.code).trim2(-5, 20, 'onset'))]
    stream_list = [(Stream(traces=[tr, ]), i) for i, tr in enumerate(stream.select(component='Q', station=st.code).trim2(-5, 20, 'onset'))]
    try:
        groups = cluster(template_list=stream_list, show=False, corr_thresh=0.3, cores=4)
        if len(groups) < len(stream_list):
            group_max = groups[0]
            for g in groups:
                if len(g) > len(group_max):
                    group_max = g
            group_max_l = Stream(traces=[stream_list_l[g[1]][0] for g in group_max])
            group_streams = [st_tuple[0] for st_tuple in group_max]
            group_streams_l = [Stream(traces=[st_tuple[0],]) for st_tuple in group_max_l]
            stack = PWS_stack(streams=group_streams)
            stack_l = PWS_stack(streams=group_streams_l)
            st_l_q = Stream(traces=[stack_l[0], stack[0]])
            st_l_q.normalize(global_max=True)
            st_q = Stream(traces=[st_l_q[1]])
            st_q.decimate(factor=2)
            convert_stream_to_ascii(st_q, 'data', len(group_max))
    except:
        print "Error: ", sys.exc_info()[0]
