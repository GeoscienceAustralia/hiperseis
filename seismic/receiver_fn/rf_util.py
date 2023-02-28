#!/usr/bin/env python
"""Utility functions to help with RF processing and analysis.
"""

from collections import defaultdict
import logging
import copy
import os, re
import numpy as np
from scipy import signal
from scipy.signal import hilbert, correlate

import obspy
import rf
import h5py

from seismic.stream_processing import assert_homogenous_stream
from seismic.receiver_fn.rf_network_dict import NetworkRFDict

# pylint: disable=invalid-name, logging-format-interpolation

logging.basicConfig()

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
# end func

def trim_hdf_keys(hdf_key_list:[str], networks_string:str, stations_string:str) -> [str]:
    """
    Trims a list of hdf_keys, filtering out unwanted networks and stations.
    @param hdf_key_list:
    @param networks_string: a space-separated list of networks. '*' includes all.
    @param stations_string: a space-separated list of stations or a text file
                            with station names in each row, w/wo location codes.
                            '*' includes all.
    @return: trimmed list
    """

    network_list = []
    station_list = []

    if(networks_string=='*'):
        network_list = []
    else:
        network_list = re.findall('\S+', networks_string)
        assert len(network_list), 'Invalid network list. Aborting..'
    # end if

    if(stations_string=='*'):
        station_list = []
    else:
        stations = []
        if(os.path.isfile(stations_string)):
            for iline, line in enumerate(open(stations_string, 'r').readlines()):
                if(not line.strip()): continue
                else: stations.append(line)
            # end for

            stations_string = ' '.join(stations)
        # end if

        station_list = re.findall('\S+', stations_string)
        assert len(station_list), 'Invalid station list. Aborting..'
    # end if

    # sanity checks
    for net in network_list:
        if net not in [hdf_key.split('.')[0] for hdf_key in hdf_key_list]:
            assert 0, 'Network {} not found in input dataset. Aborting..'.format(net)
    # end for

    for sta in station_list:
        if sta.split('.')[0] not in [hdf_key.split('.')[1] for hdf_key in hdf_key_list]:
            assert 0, 'Station {} not found in input dataset. Aborting..'.format(sta)
    # end for

    net_subset = [] # filter networks
    if(len(network_list)):
        for hdf_key in hdf_key_list:
            net, sta, loc = hdf_key.split('.')

            if(net in network_list): net_subset.append(hdf_key)
        # end for
    else:
        net_subset = hdf_key_list.copy()
    # end if

    sta_subset = [] # filter stations
    if(len(station_list)):
        for hdf_key in net_subset:
            net, sta, loc = hdf_key.split('.')

            if (sta in station_list): sta_subset.append(hdf_key)
            if ('.'.join([sta, loc]) in station_list): sta_subset.append(hdf_key)
        # end for
    else:
        sta_subset = net_subset.copy()
    # end if

    return sta_subset
# end func

def phase_weights(stream):
    """Phase weighting takes all the traces in a stream and computes a weighting for each sample in the
    stream between 0 and 1. The weighting represents how consistent is the phase angle of the signal
    at the same point in the time series across all streams.

    If phase weighting to accentuate higher multiples than Ps, then moveout should be applied first
    before calling this function.

    See https://doi.org/10.1111/j.1365-246X.1997.tb05664.x

    Note: this function should not be applied to streams with mixed components.

    :param stream: Stream containing one or more traces from which phase coherence weightings will be generated.
    :type stream: Iterable container of obspy.Trace
    :return: Array of normalized weighting factors with same length as traces in stream.
    :rtype: numpy.array
    """
    assert_homogenous_stream(stream, phase_weights.__name__)

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
# end func


def find_rf_group_ids(stream):
    """For the given stream, which is expected to have an rf_group attribute in its traces' metadata, determine
    the unique set of group ids that the traces contain.

    :param stream: Stream containing traces with rf_group ids associated with them.
    :type stream: obspy.Stream
    :return: Set of rf_group ids found in the traces
    :rtype: set(int)
    """
    # AttributeError may be raised here if rf_group attribute does not exist in the stats.
    group_ids = set((trace.stats.rf_group for trace in stream))
    return group_ids
# end func


def read_h5_rf(src_file, network=None, station=None, loc='', root='/waveforms'):
    """Helper function to load data from hdf5 file generated by rf library or script `rf_quality_filter.py`.
    For faster loading time, a particular network and station may be specified.

    :param src_file: File from which to load data
    :type src_file: str or Path
    :param network: Specific network to load, defaults to None
    :type network: str, optional
    :param station: Specific station to load, defaults to None
    :type station: str, optional
    :param root: Root path in hdf5 file where to start looking for data, defaults to '/waveforms'
    :type root: str, optional
    :return: All the loaded data in a rf.RFStream container.
    :rtype: rf.RFStream
    """
    logger = logging.getLogger(__name__)
    if (network is None and station is not None) or (network is not None and station is None):
        logger.warning("network and station should both be specified - IGNORING incomplete specification")
        group = root
    elif network and station:
        group = root + '/{}.{}.{}'.format(network.upper(), station.upper(), loc.upper())
    else:
        group = root
    # end if

    rf_data = rf.read_rf(src_file, format='h5', group=group)
    return rf_data
# end func


def rf_to_dict(rf_data):
    """Convert RF data loaded from function read_h5_rf() into a dict format for easier addressing
    of selected station and channel RF traces.

    :param rf_data: RFStream data
    :type rf_data: rf.RFStream
    :return: Nested dicts to find traces by station then channel code, with attached metadata.
    :rtype: seismic.receiver_fn.rf_network_dict.NetworkRFDict
    """
    return NetworkRFDict(rf_data)
# end func


def signed_nth_root(arr, order):
    """As per DOI https://doi.org/10.1038/217533a0.
    Muirhead, K.J. "Eliminating False Alarms when detecting Seismic Events Automatically"

    :param arr: Compute n-th root of input array, preserving sign of original data.
    :type arr: numpy.array
    :param order: Order of the root to compute
    :type order: float or int
    :return: Input array raised to 1/nth power.
    :rtype: numpy.array
    """
    if order == 1:
        return arr.copy()
    else:
        return np.sign(arr)*np.power(np.abs(arr), 1.0/order)
# end func


def signed_nth_power(arr, order):
    """As per DOI https://doi.org/10.1038/217533a0.
    Muirhead, K.J. "Eliminating False Alarms when detecting Seismic Events Automatically"

    :param arr: Compute n-th power of input array, preserving sign of original data.
    :type arr: numpy.array
    :param order: Order of the power to compute
    :type order: float or int
    :return: Input array raised to nth power.
    :rtype: numpy.array
    """
    if order == 1:
        return arr
    else:
        return np.sign(arr)*np.power(np.abs(arr), order)
# end func


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
# end func


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
# end func


def compute_extra_rf_stats(stream):
    """Compute extra statistics for each trace and add it to the RFTrace.stats structure.

    :param stream: RFStream to augment with additional metadata statistics.
    :type stream: rf.RFStream
    """
    for tr in stream:
        # RMS amplitude
        rms_amp = np.sqrt(np.mean(np.square(tr.data)))
        log10_rms_amp = np.log10(rms_amp)
        # Complex amplitude
        cplx_amp = np.abs(hilbert(tr.data))
        log10_cplx_amp = np.log10(cplx_amp)
        mean_cplx_amp = np.mean(cplx_amp)
        mean_log10_cplx_amp = np.mean(log10_cplx_amp)
        # Histogram percentiles
        log10_amp_20pc = np.percentile(log10_cplx_amp, 20)
        log10_amp_80pc = np.percentile(log10_cplx_amp, 80)
        # Store in stats metadata
        tr.stats.rms_amp = rms_amp
        tr.stats.log10_rms_amp = log10_rms_amp
        tr.stats.mean_cplx_amp = mean_cplx_amp
        tr.stats.mean_log10_cplx_amp = mean_log10_cplx_amp
        tr.stats.log10_amp_20pc = log10_amp_20pc
        tr.stats.log10_amp_80pc = log10_amp_80pc
    # end for
# end func


def compute_vertical_snr(src_stream):
    """Compute the SNR of the Z component (Z before deconvolution)
    including the onset pulse (key 'snr_prior'). Stores results in metadata of input stream traces.
    This SNR is a ratio of max envelopes.

    Some authors compute this prior SNR on signal after rotation but before deconvolution, however
    that doesn't make sense for LQT rotation where the optimal rotation will result in the least
    energy in the L component. For simplicity we compute it on Z-component only which is a reasonable
    estimate for teleseismic events.

    :param src_stream: Seismic traces before RF deconvolution of raw stream.
    :type src_stream: rf.RFStream or obspy.Stream
    """
    logger = logging.getLogger(__name__)

    if isinstance(src_stream, rf.RFStream):
        slice2 = lambda s, w: s.slice2(*w, reftime='onset')
    elif isinstance(src_stream, obspy.Stream):
        slice2 = lambda s, w: s.slice(w[0] if w[0] is None else s[0].stats.onset - w[0],
                                      w[1] if w[1] is None else s[0].stats.onset + w[1])
    else:
        assert False, "NYI"
    # end if

    def _set_nan_snr(stream):
        md_dict = {'snr_prior': np.nan}
        for tr in stream:
            tr.stats.update(md_dict)
        # end for
    # end func

    src_stream = src_stream.select(component='Z')

    # Compute max envelope amplitude from onset onwards relative to max envelope before onset.
    PRIOR_PICK_SIGNAL_WINDOW = (-5.0, 25.0)
    PRIOR_NOISE_SIGNAL_WINDOW = (None, -5.0)
    pick_signal = slice2(src_stream.copy(), PRIOR_PICK_SIGNAL_WINDOW)
    pick_signal = pick_signal.taper(0.5, max_length=0.5)
    pick_signal = np.array([tr.data for tr in pick_signal])
    if len(pick_signal.shape) == 1:
        pick_signal = pick_signal.reshape(1, -1)
    # Compute envelope of all traces
    if not np.any(pick_signal):
        _set_nan_snr(src_stream)
        return
    # end if
    pick_signal = np.absolute(signal.hilbert(pick_signal, axis=1))

    noise = slice2(src_stream.copy(), PRIOR_NOISE_SIGNAL_WINDOW)
    # Taper the slices so that the result is not overly affected by the phase of the signal at the ends.
    noise = noise.taper(0.5, max_length=0.5)
    noise = np.array([tr.data for tr in noise])
    if len(noise.shape) == 1:
        noise = noise.reshape(1, -1)
    if not np.any(noise):
        _set_nan_snr(src_stream)
        return
    # end if
    noise = np.absolute(signal.hilbert(noise, axis=1))

    if pick_signal.shape[0] != noise.shape[0]:
        logger.error("Shape inconsistency between noise and signal slices: {}[0] != {}[0]"
                     .format(pick_signal.shape, noise.shape))
        _set_nan_snr(src_stream)
    else:
        snr_prior = np.max(pick_signal, axis=1) / np.max(noise, axis=1)
        for i, tr in enumerate(src_stream):
            md_dict = {'snr_prior': snr_prior[i]}
            tr.stats.update(md_dict)
        # end for
    # end if
# end func


def compute_rf_snr(rf_stream):
    """Compute signal to noise (S/N) ratio of the RF itself about the onset pulse (key 'snr').
    This SNR is a ratio of RMS amplitudes. Stores results in metadata of input stream traces.

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

    PICK_SIGNAL_WINDOW = (1.0, 25.0)  # Consider tying the start of this window to 2x the minimum period present
    NOISE_SIGNAL_WINDOW = (None, -2.0)

    # Take everything up to 2 sec before onset as noise signal.
    noise = rf_stream.copy().slice2(*NOISE_SIGNAL_WINDOW, reftime='onset')
    # Taper the slices so that the RMS is not overly affected by the phase of the signal at the ends.
    noise = noise.taper(0.5, max_length=0.5)
    noise = np.array([tr.data for tr in noise])
    if len(noise.shape) == 1:
        noise = noise.reshape(1, -1)

    # The time window from 1 sec before to 2 sec after onset as the RF P signal
    pick_signal = rf_stream.copy().slice2(*PICK_SIGNAL_WINDOW, reftime='onset')
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
# end func


def choose_rf_source_channel(rf_type, db_station):
    """Choose source channel for RF analysis.

    :param rf_type: The RF rotation type, should be either 'ZRT' or 'LQT'
    :type rf_type: str
    :param db_station: Dict of traces for a given station keyed by channel code.
    :type db_station: dict(str, list(rf.RFTrace))
    :return: Channel code of the primary RF source channel
    :rtype: str
    """
    if rf_type[0:3].upper() == 'ZRT':
        prospective_channels = ['HHR', 'BHR', 'EHR', 'SHR', '**R']
        fallback = 'R'
    elif rf_type[0:3].upper() == 'LQT':
        prospective_channels = ['HHQ', 'BHQ', 'EHQ', 'SHQ', '**Q']
        fallback = 'Q'
    else:
        prospective_channels = []
        fallback = None
    # end if
    best_channel = None
    for c in prospective_channels:
        if c in db_station:
            best_channel = c
            break
        # end if
    # end for
    if best_channel is None:
        for c in db_station.keys():
            if c[-1] == fallback:
                best_channel = c
                break
            # end if
        # end for
    # end if
    return best_channel
# end func


def label_rf_quality_simple_amplitude(rf_type, traces, snr_cutoff=2.0, rms_amp_cutoff=0.2, max_amp_cutoff=1.0):
    """Add RF quality label for a collection of RFs based on simple amplitude criteria computed by
    quality filter script.  Adds quality label in-place.

    :param rf_type: The RF rotation type, should be either 'ZRT' or 'LQT'
    :type rf_type: str
    :param traces: Iterable collection of rf.RFTrace
    :type traces: Iterable collection of rf.RFTrace
    :param snr_cutoff: Minimum signal SNR, defaults to 2.0
    :type snr_cutoff: float, optional
    :param rms_amp_cutoff: Maximum accepted RMS amplitude of signal, defaults to 0.2
    :type rms_amp_cutoff: float, optional
    :param max_amp_cutoff: Maximum accepted amplitude of signal, defaults to 1.0
    :type max_amp_cutoff: float, optional
    """
    def _amplitude_metric_good(tr):
        if not (tr.stats.get('snr') and tr.stats.get('log10_amp_max') and tr.stats.get('log10_amp_rms')):
            return False
        return tr.stats.snr >= snr_cutoff and \
               tr.stats.log10_amp_max <= np.log10(max_amp_cutoff) and \
               tr.stats.log10_amp_rms <= np.log10(rms_amp_cutoff)
    # end func

    # Simple SNR and amplitude based filtering criteria matching formula from Babak in Matlab code, PLUS
    # additional requirement in the case of ZRT rotation that the peak occurs nearby to onset time.
    peak_location_tolerance_sec = 2.0
    if rf_type[0:3].upper() == 'ZRT':
        for tr in traces:
            times_rel = tr.times() - (tr.stats.onset - tr.stats.starttime)
            if (_amplitude_metric_good(tr) and
                    np.min(np.abs(times_rel[np.argwhere(tr.data == np.max(tr.data))])) <= peak_location_tolerance_sec):
                tr.stats.predicted_quality = 'a'
            else:
                tr.stats.predicted_quality = 'b'
            # end if
        # end for
    else:  # LQT
        for tr in traces:
            if _amplitude_metric_good(tr):
                tr.stats.predicted_quality = 'a'
            else:
                tr.stats.predicted_quality = 'b'
            # end if
        # end for
    # end if
# end func


def filter_crosscorr_coeff(rf_stream, time_window=(-2, 25), threshold_cc=0.70, min_fraction=0.15, apply_moveout=False):
    """For each trace in the stream, compute its correlation coefficient with the other traces.
    Return only traces matching cross correlation coefficient criteria based on C.Sippl (2016)
    [see http://dx.doi.org/10.1016/j.tecto.2016.03.031]

    :param rf_stream: Stream of RF traces to filter, should be **for a single component of a single station**
    :type rf_stream: rf.RFStream
    :param time_window: Time window to filter by, defaults to (-2, 25)
    :type time_window: tuple, optional
    :param threshold_cc: Threshold cross-correlation coefficient, defaults to 0.70.
        Denoted Xi in Sippl, who used value 0.80.
    :type threshold_cc: float, optional
    :param min_fraction: Minimum fraction of coefficients above threshold_cc, defaults to 0.15.
        Denoted tau in Sippl, who used value 0.15.
    :type min_fraction: float, optional
    :param apply_moveout: Whether to apply moveout correction to Ps phase prior to computing
        correlation coefficients.
    :type apply_moveout: bool
    :return: Filtered stream of RF traces
    :rtype: rf.RFStream
    """
    assert_homogenous_stream(rf_stream, filter_crosscorr_coeff.__name__)

    # Early exit if we don't have enough traces for similarity filtering to be meaningful.
    if len(rf_stream) < 3:
        return rf_stream
    # end if

    # Trim good RFs to time range so that subsequent cross-correlation computations relate to the
    # relevant period around and after onset.
    data_cc = rf_stream.copy().trim2(*time_window, reftime='onset')
    if not data_cc:
        return data_cc
    # end if

    # Apply optional moveout
    if apply_moveout:
        data_cc.moveout()
    # end if
    # Gather all RFs into a single array for efficient computation of correlation coefficients
    # between all traces
    data_array = np.array([tr.data for tr in data_cc])
    # Compute cross-correlation coefficients. cc matrix will be symmetric.
    # Each row of cc indicates the degree of correlation between each other trace.
    cc = np.corrcoef(data_array)
    # Determine mask of which traces meet the similarity filtering criteria
    fraction_above_threshold = np.sum(cc >= threshold_cc, axis=1)/len(data_cc)
    keep_trace_mask = (fraction_above_threshold >= min_fraction)
    kept_data = rf.RFStream([tr for i, tr in enumerate(rf_stream) if keep_trace_mask[i]])
    return kept_data
# end func

def filter_invalid_radial_component(rf_stream):
    """
    Filter out invalid radial RFs with amplitudes > 1 or troughs around onset time
    :param rf_stream: Stream of RF traces to filter, should be **for a single component of a single station**
    :type rf_stream: rf.RFStream
    :return: Filtered RF stream
    :rtype: rf.RFStream
    """
    if(len(rf_stream) == 0): return rf_stream

    assert rf_stream[0].stats.channel[-1] in ('R', 'Q'), 'Invalid cahnnel; must be either R or Q. Aborting..'
    assert_homogenous_stream(rf_stream, filter_invalid_radial_component.__name__)

    rf_stream_out = []
    for trc in rf_stream:
        if ((np.max(trc.data) <= 1.0) and \
            (np.max(trc.data) > -np.min(trc.data))):
            rf_stream_out.append(trc)
        # end if
    # end for

    return rf.RFStream(rf_stream_out)
# end func

def filter_by_distance(rf_stream, min_dist, max_dist):
    """
    Discard RFs that fall outside the distance range(min_dist, max_dist)
    @param rf_stream: RFStream
    @param min_dist: minimum angular distance
    @param max_dist: maximum angular distance
    @return: trimmed RFStream
    """

    rf_stream_out = []
    for trc in rf_stream:
        if(trc.stats.distance >= min_dist and trc.stats.distance <= max_dist):
            rf_stream_out.append(trc)
        # end if
    # end for

    return rf.RFStream(rf_stream_out)
# end func
