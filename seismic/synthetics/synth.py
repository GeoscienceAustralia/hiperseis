#!/usr/bin/env python
"""
Entry point for making synthetic seismograms.
"""

import numpy as np

from seismic.synthetics.backends import backend_tws, backend_syngine
from seismic.model_properties import LayerProps


def synthesize_dataset(method, output_file, src_latlon, fs, time_window, **kwargs):
    """

    :param method: 'propmatrix' or 'syngine'
    :param output_file: Destination file in which to write resultant streams in HDF5 format.
    :return:
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
    synth_streams.write(output_file, 'H5')
# end func


if __name__ == '__main__':
    crust = LayerProps(6.4, 4.2, 2.7, 35.0)
    mantle = LayerProps(8.2, 6.8, 3.3, np.nan)
    src_latlon = 2 * (np.random.random((5, 2)) - 0.5) + np.array([30, 160])
    fs = 100.0
    time_window = (-20, 40)
    generator_args = {
        'station_latlon': (-20, 140),
        'layerprops': [crust, mantle]
    }
    synthesize_dataset('propmatrix', 'test_prop_synth.h5', src_latlon, fs, time_window, **generator_args)
    pass
# end if
