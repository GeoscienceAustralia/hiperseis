#!/usr/bin/env python
"""
Entry point for making synthetic seismograms.
"""

import os

from seismic.synthetics.backends import backend_tws, backend_syngine
from seismic.stream_io import write_h5_event_stream


def synthesize_dataset(method, output_file, net, sta, src_latlon, fs, time_window, **kwargs):
    """
    User function for creating a synthetic seismogram dataset of obspy streams in HDF5 format.
    Datasets generated can be loaded into class NetworkEventDataset.

    :param method: 'propmatrix' or 'syngine'
    :type method: str
    :param output_file: Destination file in which to write resultant streams in HDF5 format.
    :type output_file: str
    :param net: Network code of receiver
    :type net: str
    :param sta: Station code of receiver
    :type sta: str
    :param src_latlon: Iterable of source (lat, lon) locations
    :type src_latlon: iterable of pairs
    :param fs: Sampling rate
    :type fs: float
    :param time_window: Time window about onset. First value should be < 0, second should be > 0
    :type time_window: tuple(float, float)
    :return: Whether the dataset was successfully created.
    :rtype: bool
    """
    if method == 'propmatrix':
        backend = backend_tws.synthesizer()
    elif method == 'syngine':
        backend = backend_syngine.synthesizer()
    else:
        assert False, 'Method {} not supported'.format(method)
    # end if
    synthesizer = backend(**kwargs)
    synth_streams = synthesizer.synthesize(src_latlon, fs, time_window)
    for tr in synth_streams:
        tr.stats.network = net
        tr.stats.station = sta
    # end for
    # Use mode='w' to write brand new file.
    write_h5_event_stream(output_file, synth_streams, mode='w', ignore=('mseed',))
    return os.path.isfile(output_file)
# end func
