import os

import numpy as np
import rf

from seismic.synthetics.synth import synthesize_dataset
from seismic.model_properties import LayerProps
from seismic.receiver_fn.generate_rf import event_waveforms_to_rf

mantle = LayerProps(8.0, 4.5, 3.3, np.nan)
crust = LayerProps(6.4, 6.4/1.7, 2.7, 35.0)

src_latlon = ((35, 140),)
fs = 100.0
time_window = (-20, 60)
generator_args = {
    'station_latlon': (-25, 140),
    'layerprops': [crust, mantle]
}
raw_file = 'synth_events_1L.h5'
synthesize_dataset('propmatrix', raw_file, 'SY', 'AAA', src_latlon, fs, time_window, **generator_args)

# Generate receiver function
rf_file = 'synth_rf_1L.h5'
resample_rate = 6.25  # Hz
event_waveforms_to_rf(raw_file, rf_file, resample_rate, taper_limit=0.05, filter_band=(0.02, 1.0),
                      trim_start_time=-5.0 - 275*resample_rate, trim_end_time=150,
                      rotation_type='ZRT', deconv_domain='iter')

# Read in RF and convert to text format
rf_dat_file = os.path.splitext(rf_file)[0] + '.dat'
rf_data = rf.read_rf(rf_file, format='h5')
rf_data.trim2(-5.0, 25.0, reftime='onset', nearest_sample=True)
rf_data = rf_data.select(component='R')[0]
times = rf_data.times() - (rf_data.stats.onset - rf_data.stats.starttime)
with open(rf_dat_file, 'w') as f:
    for t, d in zip(times, rf_data.data):
        f.write('{:2.2f}, {:.8f}\n'.format(t, d))
# end with
