
import numpy as np

from seismic.synthetics.synth import synthesize_dataset
from seismic.model_properties import LayerProps

mantle = LayerProps(8.0, 4.5, 3.3, np.nan)
sediment = LayerProps(2.5, 1.0, 2.0, 1.0)
upper_crust = LayerProps(5.6, 5.6/1.65, 2.4, 10.0)
lower_crust = LayerProps(6.4, 6.4/1.7, 2.7, 30.0)

num_events = 10
reference_locus = np.array([15, 140])
src_latlon = 5*(np.random.randn(num_events, 2) - 0.5) + reference_locus
fs = 100.0
time_window = (-20, 60)
generator_args = {
    'station_latlon': (-20, 137),
    'layerprops': [sediment, upper_crust, lower_crust, mantle]
}
synthesize_dataset('propmatrix', 'synth_events_3L.h5', 'SY', 'OAA', src_latlon, fs, time_window, **generator_args)

