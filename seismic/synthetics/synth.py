#!/usr/bin/env python
"""
Entry point for making synthetic seismograms.
"""

import os

import numpy as np

from seismic.synthetics.backends import backend_tws, backend_syngine
from seismic.model_properties import LayerProps


def synthesize_dataset(method, output_file, net, sta, src_latlon, fs, time_window, **kwargs):
    """
    User function for creating a synthetic seismogram dataset of obspy streams in HDF5 format.
    Datasets generated can be loaded into class NetworkEventDataset.

    :param method: 'propmatrix' or 'syngine'
    :param output_file: Destination file in which to write resultant streams in HDF5 format.
    :param net: Network code of receiver
    :param sta: Station code of receiver
    :param src_latlon: Iterable of source event coordinates
    :param fs: Sampling rate
    :param time_window: Time window about onset. First value should be < 0, second should be > 0
    :return: Whether the dataset was successfully created.
    """
    backend = None
    if method == 'propmatrix':
        backend = backend_tws.synthesizer()
    elif method == 'syngine':
        backend = backend_syngine.synthesizer()
    else:
        assert False, 'Not supported'
    # end if
    synthesizer = backend(**kwargs)
    synth_streams = synthesizer.synthesize(src_latlon, fs, time_window)
    for tr in synth_streams:
        tr.stats.network = net
        tr.stats.station = sta
    # end for
    synth_streams.write(output_file, 'H5', ignore=('mseed',))
    return os.path.isfile(output_file)
# end func


if __name__ == '__main__':
    # Example using propagator matrix method
    crust = LayerProps(6.4, 4.2, 2.7, 35.0)
    mantle = LayerProps(8.2, 6.8, 3.3, np.nan)
    src_latlon = 2 * (np.random.random((5, 2)) - 0.5) + np.array([30, 160])
    fs = 100.0
    time_window = (-20, 40)
    generator_args = {
        'station_latlon': (-20, 140),
        'layerprops': [crust, mantle]
    }
    synthesize_dataset('propmatrix', 'test_prop_synth.h5', 'SY', 'AAA', src_latlon, fs, time_window, **generator_args)

    # Example using Syngine service
    generator_args = {
        'station_latlon': (-20, 140)
    }
    synthesize_dataset('syngine', 'test_syng_synth.h5', 'SY', 'AAA', src_latlon, fs, time_window, **generator_args)

# end if
