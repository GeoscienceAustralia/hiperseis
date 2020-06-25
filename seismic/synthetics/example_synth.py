#!/usr/bin/env python
"""
Example usage of synth module.
"""

import numpy as np

from seismic.model_properties import LayerProps
from seismic.synthetics.synth import synthesize_dataset


def example_propmatrix():
    # Example using propagator matrix method
    crust = LayerProps(6.4, 3.8, 2.7, 38.0)
    mantle = LayerProps(8.2, 6.8, 3.3, np.nan)
    src_latlon = 20 * (np.random.random((5, 2)) - 0.5) + np.array([30, 160])
    fs = 100.0
    time_window = (-50, 150)
    generator_args = {
        'station_latlon': (-20, 140),
        'layerprops': [crust, mantle]
    }
    synthesize_dataset('propmatrix', 'test_prop_synth.h5', 'SY', 'AAA', src_latlon, fs, time_window, **generator_args)
# end func


def example_syngine():
    # Example using Syngine service
    src_latlon = 20 * (np.random.random((5, 2)) - 0.5) + np.array([30, 160])
    fs = 100.0
    time_window = (-50, 150)
    generator_args = {
        'station_latlon': (-20, 140)
    }
    synthesize_dataset('syngine', 'test_syng_synth.h5', 'SY', 'AAA', src_latlon, fs, time_window, **generator_args)
# end func


if __name__ == '__main__':
    example_propmatrix()
    example_syngine()
# end if
