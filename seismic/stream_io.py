#!/usr/bin/env python
"""Helper functions for seismic stream IO.
"""

import os
import itertools

import obspy
# import obspyh5
from obspy.taup import TauPyModel
from obspy.io.sac.sactrace import SACTrace

from seismic.receiver_fn.rf_util import KM_PER_DEG

def sac2hdf5(src_folder, basenames, channels, dest_h5_file, tt_model_id='iasp91'):
    """
    Convert collection of SAC files from a folder into a singel HDF5 stream file.

    :param src_folder:
    :param basenames:
    :param channels:
    :param dest_h5_file:
    :param tt_model_id:
    :return:
    """
    tt_model = TauPyModel(tt_model_id)
    traces = []
    for basename, channel in itertools.product(basenames, channels):
        fname = os.path.join(src_folder, '.'.join([basename, channel]))
        channel_stream = obspy.read(fname, 'sac')
        tr = channel_stream[0]
        event_depth_km = tr.stats.sac['evdp']
        dist_deg = tr.stats.sac['dist'] / KM_PER_DEG
        arrivals = tt_model.get_travel_times(event_depth_km, dist_deg, ('P',))
        arrival = arrivals[0]
        inc = arrival.incident_angle
        slowness = arrival.ray_param_sec_degree
        src_dic = tr.stats.sac
        sac_tr = SACTrace(nzyear=src_dic['nzyear'], nzjday=src_dic['nzjday'], nzhour=src_dic['nzhour'],
                          nzmin=src_dic['nzmin'], nzsec=src_dic['nzsec'], nzmsec=src_dic['nzmsec'])
        onset = sac_tr.reftime
        if 'nevid' in src_dic:
            event_id = src_dic['nevid']
        else:
            event_id = basename
        # end if
        stats = {'distance': dist_deg, 'back_azimuth': src_dic['baz'], 'inclination': inc,
                 'onset': onset, 'slowness': slowness, 'phase': 'P', 'tt_model': tt_model_id,
                 'event_id': event_id}
        tr.stats.update(stats)
        traces.append(tr)
    # end for

    stream_all = obspy.Stream(traces)
    stream_all.write(dest_h5_file, 'H5')

    return os.path.isfile(dest_h5_file)
# end func
