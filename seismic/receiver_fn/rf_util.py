#!/usr/bin/env python
"""Utility functions to help with RF processing and analysis.
"""

from collections import defaultdict
import logging
import copy

import numpy as np
from scipy import signal
from scipy.signal import hilbert

import rf

# pylint: disable=invalid-name

KM_PER_DEG = 111.1949

logging.basicConfig()


def phase_weights(stream):
    """Phase weighting takes all the traces in a stream and computes a weighting for each sample in the
    stream between 0 and 1. The weighting represents how consistent is the phase angle of the signal
    at the same point in the time series across all streams.

    If phase weighting to accentuate higher multiples than Ps, then moveout should be applied first
    before calling this function.

    See https://doi.org/10.1111/j.1365-246X.1997.tb05664.x

    :param stream: Stream containing one or more traces from which phase coherence weightings will be generated.
    :type stream: obspy.core.stream.Stream
    :return: Array of normalized weighting factors with same length as traces in stream.
    :rtype: numpy.array
    """
    traces = np.array([tr.data for tr in stream])
    # Hilbert transform to separate complex amplitude from complex phase.
    analytic = hilbert(traces)
    # The complex part of the hilbert transform contains sine of the phase angle.
    # numpy.angle extracts the angle in radians from the imaginary component.
    angle = np.angle(analytic)
    # Using just the phase angle component (throwing away the amplitude) generate complex number
    # representing just the phase angle of the signal at each sample.
    i_phase = np.exp(1j * angle)
    # Compute the mean of all the complex phases. If they're random, due to summing in the complex
    # domain they will tend to cancel out. If they're not random, they will tend to sum coherently and
    # generated a large stacked complex amplitude.
    tphase = np.abs(np.mean(i_phase, axis=0))
    # Return normalized result against max amplitude, so that the most coherent part of the signal
    # has a scaling of 1.
    return tphase/np.max(tphase)


def find_rf_group_ids(stream):
    """For the given stream, which is expected to have an rf_group attribute in its traces' metadata, determine
    the unique set of group ids that the traces contain.

    :param stream: Stream containing traces with rf_group ids associated with them.
    :type stream: obspy.core.trace.Trace
    :return: Set of rf_group ids found in the traces
    :rtype: set(int)
    """
    # AttributeError may be raised here if rf_group attribute does not exist in the stats.
    group_ids = set((trace.stats.rf_group for trace in stream))
    return group_ids


def read_h5_rf(src_file, network=None, station=None, root='/waveforms'):
    """Helper function to load data from hdf5 file generated by rf library or script `rf_quality_filter.py`

    :param src_file: [description]
    :type src_file: [type]
    :param network: [description], defaults to None
    :type network: [type], optional
    :param station: [description], defaults to None
    :type station: [type], optional
    :param root: [description], defaults to '/waveforms'
    :type root: str, optional
    :return: [description]
    :rtype: [type]
    """
    logger = logging.getLogger(__name__)
    if (network is None and station is not None) or (network is not None and station is None):
        logger.warning("network and station should both be specified - IGNORING incomplete specification")
    elif network and station:
        group = root + '/{}.{}.0M'.format(network.upper(), station.upper())
    else:
        group = root
    # end if

    rf_data = rf.read_rf(src_file, format='h5', group=group)
    return rf_data


def rf_to_dict(rf_data):
    """Convert RF data loaded from function read_h5_rf() into a dict format for easier addressing
    of selected station and channel RF traces.

    :param rf_data: [description]
    :type rf_data: [type]
    :return: [description]
    :rtype: [type]
    """
    db = defaultdict(lambda: defaultdict(list))
    for s in rf_data:
        _, sta, _, cha = s.id.split('.')
        db[sta][cha].append(s)
    return db


def signed_nth_root(arr, order):
    """As per DOI https://doi.org/10.1038/217533a0.
    Muirhead, K.J. "Eliminating False Alarms when detecting Seismic Events Automatically"

    :param arr: [description]
    :type arr: [type]
    :param order: [description]
    :type order: [type]
    :return: [description]
    :rtype: [type]
    """
    if order == 1:
        return arr
    else:
        return np.sign(arr)*np.power(np.abs(arr), 1.0/order)


def signed_nth_power(arr, order):
    """As per DOI https://doi.org/10.1038/217533a0.
    Muirhead, K.J. "Eliminating False Alarms when detecting Seismic Events Automatically"

    :param arr: [description]
    :type arr: [type]
    :param order: [description]
    :type order: [type]
    :return: [description]
    :rtype: [type]
    """
    if order == 1:
        return arr
    else:
        return np.sign(arr)*np.power(np.abs(arr), order)


def filter_station_streams(db_station, freq_band=(None, None)):
    """Perform frequency filtering on all channels' traces for a given station.
    Returns a copy of db_station with streams containing filtered results.
    """
    db_station_filt = defaultdict(list)
    for i, (ch, streams) in enumerate(db_station.items()):
        if ch == 'size' or i >= 3:
            continue
        for s in streams:
            stream_filt = s.copy()
            if freq_band[0] is None and freq_band[1] is not None:
                stream_filt.filter('lowpass', zerophase=True, corners=2, freq=freq_band[1])
            elif freq_band[0] is not None and freq_band[1] is None:
                stream_filt.filter('highpass', zerophase=True, corners=2, freq=freq_band[0])
            elif freq_band[0] is not None and freq_band[1] is not None:
                stream_filt.filter('bandpass', zerophase=True, corners=2, freqmin=freq_band[0],
                                   freqmax=freq_band[1])
            # end if
            db_station_filt[ch].append(stream_filt)
        # end for
    # end for

    return db_station_filt


def filter_station_to_mean_signal(db_station, min_correlation=1.0):
    """Filter out streams which are not 'close enough' to the mean signal,
       based on simple correlation score.  The term "correlation" here really
       just means a similarity dot product (projection of individual trace onto
       the mean).
    """
    # Compute mean signals of channels in station
    mean_rfs = []
    for i, (ch, streams) in enumerate(db_station.items()):
        if ch == 'size' or i >= 3:
            continue
        for j, s in enumerate(streams):
            if j == 0:
                data_mean = copy.deepcopy(s.data)
            else:
                data_mean += s.data
            # end if
        # end for
        data_mean /= np.max(data_mean)
        mean_rfs.append(data_mean)
    # end for

    # Filter out signals that do not meet minimum coherence with mean signal for each channel
    db_station_filt = defaultdict(list)

    corrs = []
    for i, (ch, streams) in enumerate(db_station.items()):
        for j, s in enumerate(streams):
            corr = np.dot(s.data, mean_rfs[i])/np.dot(mean_rfs[i], mean_rfs[i])
            if corr >= min_correlation:
                db_station_filt[ch].append(s)
            corrs.append(corr)
        # end for
    # end for

    return db_station_filt, corrs


def compute_extra_rf_stats(db_station):
    """Compute extra statistics for each trace and add it to the RFTrace.stats structure.

    :param db_station: Dictionary with list of traces per channel for a given station.
    :type db_station: dict({str, list(RFTrace)})
    """
    for _ch, traces in db_station.items():
        for tr in traces:
            rms_amp = np.sqrt(np.mean(np.square(tr.data)))
            cplx_amp = np.abs(hilbert(tr.data))
            mean_cplx_amp = np.mean(cplx_amp)
            amp_20pc = np.percentile(cplx_amp, 20)
            amp_80pc = np.percentile(cplx_amp, 80)
            tr.stats.rms_amp = rms_amp
            tr.stats.mean_cplx_amp = mean_cplx_amp
            tr.stats.amp_20pc = amp_20pc
            tr.stats.amp_80pc = amp_80pc
        # end for
    # end for


def compute_vertical_snr(src_stream):
    """Compute the SNR of the Z component (Z before rotation or deconvolution)
    including the onset pulse (key 'snr_prior'). Stored results in metadata of input stream traces.

    Some authors compute this prior SNR on signal after rotation but before deconvolution, however
    that doesn't make sense for LQT rotation where the optimal rotation will result in the least
    energy in the L component. For simplicity we compute it on Z-component only which is a reasonable
    estimate for teleseismic events.

    :param src_stream: Seismic traces before rotation of raw stream.
    :type src_stream: rf.RFStream
    """
    logger = logging.getLogger(__name__)

    src_stream = src_stream.select(component='Z')

    # Compute max envelope amplitude from onset onwards relative to max envelope before onset.
    PRIOR_PICK_SIGNAL_WINDOW = (-5.0, 25.0)
    PRIOR_NOISE_SIGNAL_WINDOW = (None, -5.0)
    pick_signal = src_stream.slice2(*PRIOR_PICK_SIGNAL_WINDOW, reftime='onset')
    pick_signal = pick_signal.taper(0.5, max_length=0.5)
    pick_signal = np.array([tr.data for tr in pick_signal])
    if len(pick_signal.shape) == 1:
        pick_signal = pick_signal.reshape(1, -1)
    # Compute envelope of all traces
    pick_signal = np.absolute(signal.hilbert(pick_signal, axis=1))

    noise = src_stream.slice2(*PRIOR_NOISE_SIGNAL_WINDOW, reftime='onset')
    # Taper the slices so that the result is not overly affected by the phase of the signal at the ends.
    noise = noise.taper(0.5, max_length=0.5)
    noise = np.array([tr.data for tr in noise])
    if len(noise.shape) == 1:
        noise = noise.reshape(1, -1)
    noise = np.absolute(signal.hilbert(noise, axis=1))

    if pick_signal.shape[0] != noise.shape[0]:
        logger.error("Shape inconsistency between noise and signal slices: {}[0] != {}[0]"
                     .format(pick_signal.shape, noise.shape))
        md_dict = {'snr_prior': np.nan}
        for tr in src_stream:
            tr.stats.update(md_dict)
        # end for
    else:
        snr_prior = np.max(pick_signal, axis=1) / np.max(noise, axis=1)
        for i, tr in enumerate(src_stream):
            md_dict = {'snr_prior': snr_prior[i]}
            tr.stats.update(md_dict)
        # end for
    # end if


def compute_rf_snr(rf_stream):
    """Compute signal to noise (S/N) ratio of the RF itself about the onset pulse (key 'snr').

    In the LQT rotation case when rotation is working ideally, the onset pulse of the rotated transverse
    signals should be minimal, and a large pulse at t = 0 indicates lack of effective rotation of coordinate
    system, so for 'snr' we use a long time window after onset pulse, deliberately excluding the onset
    pulse, to maximize contribution to the SNR from the multiples after the onset pulse.

    :param rf_stream: R or Q component of Receiver Function
    :type rf_stream: rf.RFStream
    :return: SNR for each trace in the input stream
    :rtype: numpy.array
    """
    logger = logging.getLogger(__name__)

    PICK_SIGNAL_WINDOW = (1.0, 20.0)  # Consider tying the start of this window to 2x the minimum period present
    NOISE_SIGNAL_WINDOW = (None, -2.0)

    # Take everything up to 2 sec before onset as noise signal.
    noise = rf_stream.slice2(*NOISE_SIGNAL_WINDOW, reftime='onset')
    # Taper the slices so that the RMS is not overly affected by the phase of the signal at the ends.
    noise = noise.taper(0.5, max_length=0.5)
    noise = np.array([tr.data for tr in noise])
    if len(noise.shape) == 1:
        noise = noise.reshape(1, -1)

    # The time window from 1 sec before to 2 sec after onset as the RF P signal
    pick_signal = rf_stream.slice2(*PICK_SIGNAL_WINDOW, reftime='onset')
    pick_signal = pick_signal.taper(0.5, max_length=0.5)
    pick_signal = np.array([tr.data for tr in pick_signal])
    if len(pick_signal.shape) == 1:
        pick_signal = pick_signal.reshape(1, -1)

    if pick_signal.shape[0] != noise.shape[0]:
        logger.error("Shape inconsistency between noise and signal slices: {}[0] != {}[0]"
                     .format(pick_signal.shape, noise.shape))
        md_dict = {'snr': np.nan}
        for tr in rf_stream:
            tr.stats.update(md_dict)
        # end for
    else:
        snr = np.sqrt(np.mean(np.square(pick_signal), axis=1) / np.mean(np.square(noise), axis=1))
        for i, tr in enumerate(rf_stream):
            md_dict = {'snr': snr[i]}
            tr.stats.update(md_dict)
        # end for
    # end if
