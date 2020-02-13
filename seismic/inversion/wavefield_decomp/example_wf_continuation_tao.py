#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy import stats

from seismic.inversion.wavefield_decomp.network_event_dataset import NetworkEventDataset
from seismic.inversion.wavefield_decomp.wavefield_continuation_tao import WfContinuationSuFluxComputer
from seismic.inversion.wavefield_decomp.model_properties import LayerProps
from seismic.stream_quality_filter import curate_stream3c
from seismic.receiver_fn.rf_util import compute_vertical_snr


src_file = (r"/g/data/ha3/am7399/shared/OA_RF_analysis/" +
            r"OA_event_waveforms_for_rf_20170911T000036-20181128T230620_rev8.h5")
data_all = NetworkEventDataset(src_file, network='OA', station='BT23', location='0M')

# Time window of original data to use for processing. All traces must have at least this extent
# about the onset time.
TIME_WINDOW = (-20, 50)
# Narrower time window used for integration of energy flux
FLUX_WINDOW = (-10, 20)
# Cut window for selecting central wavelet
CUT_WINDOW = (-5, 30)

#------------------------------------------------------------------------------
# Apply windowinf, filtering and QC to loaded dataset before passing to Tao's algorithm.


def stream_snr_compute(stream):
    stream.taper(0.05)
    compute_vertical_snr(stream)
# end func

def amplitude_nominal(stream, max_amplitude):
    return ((np.max(np.abs(stream[0].data)) <= max_amplitude) and
            (np.max(np.abs(stream[1].data)) <= max_amplitude) and
            (np.max(np.abs(stream[2].data)) <= max_amplitude))
# end func

# Trim streams to time window
data_all.apply(lambda stream:
               stream.trim(stream[0].stats.onset + TIME_WINDOW[0], stream[0].stats.onset + TIME_WINDOW[1]))

# Apply curation to streams prior to rotation
data_all.curate(lambda _, evid, stream: curate_stream3c(evid, stream))

# Rotate to ZRT coordinates
data_all.apply(lambda stream: stream.rotate('NE->RT'))

# Detrend the traces
data_all.apply(lambda stream: stream.detrend('linear'))

# Run high pass filter to remove high amplitude, low freq noise, if present.
f_min = 0.05
data_all.apply(lambda stream: stream.filter('highpass', freq=f_min, corners=2, zerophase=True))

# Compute SNR of Z component to use as a quality metric
data_all.apply(stream_snr_compute)

# Filter by SNR
data_all.curate(lambda _1, _2, stream: stream[0].stats.snr_prior >= 3.0)

# It does not make sense to filter by similarity, since these are raw waveforms, not RFs,
# and the waveform will be dominated by the source waveform which differs for each event.

# Filter streams with incorrect number of traces
discard = []
for sta, ev_db in data_all.by_station():
    num_pts = np.array([tr.stats.npts for st in ev_db.values() for tr in st])
    expected_pts = stats.mode(num_pts)[0][0]
    for evid, stream in ev_db.items():
        if ((stream[0].stats.npts != expected_pts) or
            (stream[1].stats.npts != expected_pts) or
            (stream[2].stats.npts != expected_pts)):
            discard.append(sta, evid)
        # end if
    # end for
# end for
data_all.prune(discard)

# Filter streams with spuriously high amplitude
MAX_AMP = 10000
data_all.curate(lambda _1, _2, stream: amplitude_nominal(stream, MAX_AMP))


#------------------------------------------------------------------------------
# Pass cleaned up data set for test station to flux computer class.
data_OA = data_all.station('BT23')
flux_comp = WfContinuationSuFluxComputer(data_OA, TIME_WINDOW, CUT_WINDOW)

# Define bulk properties of mantle (lowermost half-space)
mantle_props = LayerProps(vp=8.0, vs=4.5, rho=3.3, thickness=np.Infinity)

# Assumed crust property constants
Vp_c = 6.4
rho_c = 2.7

# Assumed sediment property constants
# Vp_s = 2.1
# rho_s = 1.97

crust_props = LayerProps(Vp_c, 3.7, rho_c, 35)

single_layer_model = [crust_props]

energy, energy_per_event, wf_mantle = flux_comp(mantle_props, single_layer_model)
print(energy)
