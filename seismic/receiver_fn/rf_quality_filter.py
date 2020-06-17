#!/usr/bin/env python
"""Filter out invalid RFs, and compute a variety of quality metrics for the
remaining RFs. These metrics can be combined in various ways downstream to
perform different kinds of filtering. Quality metrics are stored in the stats
of each trace.
"""

from builtins import range  # pylint: disable=redefined-builtin

import logging
import itertools
import traceback
from multiprocessing import Process, Manager

import numpy as np
import click
import h5py

from scipy import signal
from scipy import stats
# from fastdtw import fastdtw

from sklearn.cluster import DBSCAN
from joblib import Parallel, delayed
import rf

from seismic.receiver_fn.rf_process_io import async_write
from seismic.receiver_fn.rf_h5_file_station_iterator import IterRfH5StationEvents
from seismic.receiver_fn import rf_util

logging.basicConfig()

# pylint: disable=invalid-name, logging-format-interpolation


def _crossSpectrum(x, y, num_subsegs=32):
    # -------------------Remove mean-------------------
    # nperseg chosen arbitrary based on 126 samples RF signal, experiment to get best results
    nperseg = int(np.floor(x.size/num_subsegs))
    cross = np.zeros(nperseg, dtype='complex128')
    max_ind = int(nperseg * np.floor(x.size / nperseg))
    for ind in range(0, max_ind, nperseg):
        xp = x[ind: ind + nperseg]
        yp = y[ind: ind + nperseg]
        xp = xp - np.mean(xp)
        yp = yp - np.mean(xp)
        # Do FFT
        cfx = np.fft.fft(xp)
        cfy = np.fft.fft(yp)
        # Get cross spectrum
        cross += cfx.conj()*cfy
    # end for
    freq = np.fft.fftfreq(nperseg)
    return cross, freq


def _coh(y, y2):
    # This subroutine determines a coherence level between two signals on normilised frequency
    p11, freq = _crossSpectrum(y, y)
    p22, freq = _crossSpectrum(y2, y2)
    p12, freq = _crossSpectrum(y, y2)
    # coherence
    part1 = np.divide(np.abs(p12)**2, p11.real,
                      out=np.zeros_like(np.abs(p12)**2), where=p11.real != 0)
    coh = np.divide(part1, p22.real, out=np.zeros_like(part1), where=p22.real != 0)

    return freq[freq > 0], coh[freq > 0]


def rf_group_by_similarity(swipe, similarity_eps):
    """Cluster waveforms by similarity

    :param swipe: Numpy array of RF rowwise
    :type swipe: numpy.array
    :param similarity_eps: Tolerance on similarity between traced to be considered in the same group.
    :type similarity_eps: float
    :return: Index of the group for each trace. -1 if no group is found for a given trace.
    :rtype: numpy.array
    """

    def _compare_pairs_rmsdist(data):
        rms_dist = np.sqrt(np.nanmean(np.square(data[0] - data[1])))
        return rms_dist
    # end func

    distance = list(map(_compare_pairs_rmsdist, itertools.combinations(swipe, 2)))
    index = list((i, j) for ((i, _), (j, _)) in itertools.combinations(enumerate(swipe), 2))

    # First check that distance between points
    index = np.array(index)
    distance = np.array(distance)
    matrix = np.zeros((np.amax(index) + 1, 1 + np.amax(index))) + np.amax(distance)

    matrix[index[:, 0], index[:, 1]] = distance[:]
    clustering = DBSCAN(eps=similarity_eps, min_samples=5, metric='precomputed', n_jobs=-3).fit_predict(matrix)

    return clustering


def compute_max_coherence(orig_stream, f1, f2):
    """ Finding coherence between two signals in frequency domain
        f1 and f2 - normalised min and max frequencies, f1 < f2 <= ~0.5
        returns array of indexes for coherent traces with median

        Suggested minimum level of coherence for good results: 0.6
    """
    assert 0.0 <= f1 < f2 <= 1.0

    # Copy stream since moveout modifies self
    stream = orig_stream.copy()
    stream.moveout()

    # Clip off the noise - take everything after 2 sec before onset as seismic event signal.
    SIGNAL_WINDOW = (-2.0, None)
    pick_signal = stream.slice2(*SIGNAL_WINDOW, reftime='onset')
    # Taper at ends to ensure that end discontinuities do not inject spectral energy at high frequencies
    pick_signal = pick_signal.taper(0.5, max_length=1.0)

    data = np.array([tr.data for tr in pick_signal])
    median = np.median(data, axis=0)
    assert len(median) == data.shape[1]

    max_coh = []
    # TODO: Vectorize this loop
    for i in range(data.shape[0]):
        f, c = _coh(median, data[i, :])
        max_coh.append(np.amax(c[(f >= f1) & (f <= f2)]))
    # end for

    return np.array(max_coh)


# def knive(swipe, k_level1, sn_level2):
#     ind = np.ones((swipe.shape[0],), bool)
#     dev = []
#     ch = []
#     knive = []
#     sn = []
#     pulse_ind = np.max(t[t < 0])-1.
#
#     for i in range(swipe.shape[0]):
#         ch.append(np.amax(coh(average, swipe[i, :])))
#         dev.append(np.sum((swipe[i, :]-average)**2)/(swipe.shape[0]-1))
#         ind[i] = False
#         knive.append(np.std(swipe[ind]))
#         sn.append(np.std(swipe[i, t > pulse_ind]) /
#                   np.std(swipe[i, t < pulse_ind]))
#
#     knive = np.array(knive)
#     sn = np.array(sn)
#     return knive < k_level, sn > sn_level


def spectral_entropy(stream):
    """Compute the spectral entropy of a trace

    :param trace: Single channel seismic trace
    :type trace: rf.RFTrace
    :return: Spectral entropy of the trace waveform
    :rtype: float
    """
    data = np.array([tr.data for tr in stream])
    _, psd = signal.periodogram(data, detrend='linear')
    psd = psd.T/np.sum(psd, axis=1)
    psd = psd.T
    entropy = -np.sum(psd*np.log(psd), axis=1)
    return entropy


def get_rf_stream_components(stream):
    """Identify the RF component types and return them.

    :param stream: Stream containing mixed RF components.
    :type stream: rf.RFStream
    :return: (RF component type, primary RF component (R or Q), transverse RF component (T), source component (Z or L))
    :rtype: (str, rf.RFStream, rf.RFStream, rf.RFStream)
    """

    DEFAULT_RF_TYPE = 'LQT-Q'
    rf_type = DEFAULT_RF_TYPE

    primary_stream = stream.select(component='Q')

    if primary_stream:
        transverse_stream = stream.select(component='T')
        source_stream = stream.select(component='L')
    else:
        rf_type = 'ZRT-R'
        primary_stream = stream.select(component='R')
        if primary_stream:
            transverse_stream = stream.select(component='T')
            source_stream = stream.select(component='Z')
            # Only return Z traces which are labelled as type rf, as we shouldn't return raw trace data here.
            source_stream = source_stream.__class__(
                [tr for tr in source_stream if tr.stats.get('type') is not None and tr.stats['type'] == 'rf'])
        else:
            return None, None, None, None
        # end if
    # end if
    return rf_type, primary_stream, transverse_stream, source_stream


def compute_rf_quality_metrics(station_id, station_stream3c, similarity_eps):
    """Top level function for adding quality metrics to trace metadata.

    :param station_id: Station ID
    :type station_id: str
    :param station_stream3c: 3-channel stream
    :type station_stream3c: list(rf.RFStream) with 3 components
    :param similarity_eps: Distance threshold used for DBSCAN clustering
    :type similarity_eps: float
    :return Triplet of RF streams with Z, R or Q, and T components with populated
        quality metrics. Otherwise return None in case of failure.
    """

    logger = logging.getLogger(__name__)

    # Filter out traces with NaNs - simplifies downstream code so that can don't have to worry about NaNs.
    # We use the fact that traces are bundled into 3-channel triplets here to discard all or none of the related
    # channels for an event.
    nonan_streams = []
    for stream in station_stream3c:
        skip_stream = False
        for tr in stream:
            if tr.stats.type == 'rf' and np.any(np.isnan(tr.data)):
                logger.warning("NaN data found in station {} trace\n{}\n, skipping!".format(station_id, tr))
                skip_stream = True
                break
        # end for
        if skip_stream:
            continue
        nonan_streams.append(stream)
    # end for
    if len(nonan_streams) < len(station_stream3c):
        num_supplied = len(station_stream3c)
        num_discarded = num_supplied - len(nonan_streams)
        logger.info("Discarded {}/{} events from station {} due to NaNs in at least one channel"
                    .format(num_discarded, num_supplied, station_id))
    # end if

    # Early exit if nothing left
    if not nonan_streams:
        logger.warning("nonan_streams empty after filtering out nan traces! {}. Skipping station {}"
                       .format(nonan_streams, station_id))
        return None
    # end if

    # Flatten the traces into a single RFStream for subsequent processing
    rf_streams = rf.RFStream([tr for stream in nonan_streams for tr in stream if tr.stats.type == 'rf'])

    # Subsequent functions process the data in bulk square matrices, so it is essential all traces are the same length.
    # If not, processing will fail due to incompatible data structure. So here we filter out traces that do not have
    # the expected length. Expected length is assumed to be the most common length amongst all the traces.
    num_traces_before = len(rf_streams)
    all_trace_lens = np.array([len(tr) for tr in rf_streams])
    expected_len, _ = stats.mode(all_trace_lens, axis=None)
    expected_len = expected_len[0]
    if expected_len <= 1:
        logger.warning("Cannot compute quality metrics on trace length {} <= 1! Skipping station {}"
                       .format(expected_len, station_id))
        return None
    # end if
    keep_traces = []
    for tr in rf_streams:
        if len(tr) != expected_len:
            logger.error("Trace {} of station {} has inconsistent sample length {} (expected {}), discarding!"
                         .format(tr, station_id, len(tr), expected_len))
        else:
            keep_traces.append(tr)
        # end if
    # end for

    streams = rf.RFStream(keep_traces)
    num_traces_after = len(streams)
    if num_traces_after < num_traces_before:
        num_discarded = num_traces_before - num_traces_after
        logger.warning("Discarded {}/{} traces due to inconsistent trace length"
                       .format(num_discarded, num_traces_before))
    # end if

    # Extract RF type, the primary polarized component and transverse component (ignore source stream)
    rf_type, p_stream, t_stream, z_stream = get_rf_stream_components(streams)
    if rf_type is None:
        logger.error("Unrecognized RF type for station {}. File might not be RF file!".format(station_id))
        return None
    # end if

    # Note that we only compute quality metrics on the p_stream. The filtering of t_stream traces should match
    # the filtering of p_stream traces, so t_stream does not need independent metrics.

    # Compute S/N ratios for primary component RFs
    rf_util.compute_rf_snr(p_stream)

    # Compute spectral entropy for primary component RFs
    sp_entropy = spectral_entropy(p_stream)
    for i, tr in enumerate(p_stream):
        md_dict = {'entropy': sp_entropy[i]}
        tr.stats.update(md_dict)
    # end for

    # Compute log10 of amplitude metrics, as these are more useful than straight amplitudes for quality classifier
    for tr in p_stream:
        tr.stats['log10_amp_max'] = np.log10(tr.stats['amp_max'])
        tr.stats['log10_amp_rms'] = np.log10(tr.stats['amp_rms'])
        tr.stats['log10_z_amp_max'] = np.log10(tr.stats['z_amp_max'])
        tr.stats['log10_z_amp_rms'] = np.log10(tr.stats['z_amp_rms'])
    # end for

    # Define time windows relative to onset for computing statistical ratios
    EVENT_SIGNAL_WINDOW = (-5.0, 25.0)
    NOISE_SIGNAL_WINDOW = (None, -5.0)
    event_signal = p_stream.copy().slice2(*EVENT_SIGNAL_WINDOW, reftime='onset').taper(0.5, max_length=1.0)
    noise_signal = p_stream.copy().slice2(*NOISE_SIGNAL_WINDOW, reftime='onset').taper(0.5, max_length=1.0)
    rf_util.compute_extra_rf_stats(event_signal)
    rf_util.compute_extra_rf_stats(noise_signal)
    for _i, _tr in enumerate(p_stream):
        _tr.stats['delta_mean_log10_cplx_amp'] = (event_signal[_i].stats.mean_log10_cplx_amp -
                                                  noise_signal[_i].stats.mean_log10_cplx_amp)
        _tr.stats['delta_log10_amp_20pc'] = (event_signal[_i].stats.log10_amp_20pc -
                                             noise_signal[_i].stats.log10_amp_20pc)
        _tr.stats['delta_log10_amp_80pc'] = (event_signal[_i].stats.log10_amp_80pc -
                                             noise_signal[_i].stats.log10_amp_80pc)
        _tr.stats['delta_log10_rms_amp'] = event_signal[_i].stats.log10_rms_amp - noise_signal[_i].stats.log10_rms_amp
    # end for

    # Compute ratios of spectral histogram statistics
    noise_data = np.array([tr.data for tr in noise_signal])
    event_data = np.array([tr.data for tr in event_signal])
    noise_bins, noise_power = signal.welch(noise_data, detrend='linear')
    event_bins, event_power = signal.welch(event_data, detrend='linear')
    # Compute moments of the frequency distribution. Only use lower frequency bands up to 1/4 Nyquist.
    noise_bins = noise_bins[0:32]
    noise_power = noise_power[:, 0:32]
    event_bins = event_bins[0:32]
    event_power = event_power[:, 0:32]
    noise_m0 = np.sum(noise_power, axis=1)
    event_m0 = np.sum(event_power, axis=1)
    spectral_m0_ratio = np.log10(event_m0/noise_m0)
    noise_m1 = np.sum(noise_power*noise_bins, axis=1)
    event_m1 = np.sum(event_power*event_bins, axis=1)
    spectral_m1_ratio = np.log10(event_m1/noise_m1)
    noise_m2 = np.sum(noise_power*noise_bins**2, axis=1)
    event_m2 = np.sum(event_power*event_bins**2, axis=1)
    spectral_m2_ratio = np.log10(event_m2/noise_m2)
    for i, tr in enumerate(p_stream):
        md_dict = {
            'm0_delta': event_m0[i] - noise_m0[i],
            'm1_delta': event_m1[i] - noise_m1[i],
            'm2_delta': event_m2[i] - noise_m2[i],
            'm0_ratio': spectral_m0_ratio[i],
            'm1_ratio': spectral_m1_ratio[i],
            'm2_ratio': spectral_m2_ratio[i]
        }
        tr.stats.update(md_dict)
    # end for

    # Compute coherence metric within targeted normalized frequency band.
    # Note that settings here are relative to the sampling rate. If the sampling
    # rate changes and you want the same absolute frequency range to be used for
    # coherence, then these settings need to be updated.
    fn_low = 0.15
    fn_high = 0.3
    max_coherence = compute_max_coherence(p_stream, fn_low, fn_high)
    for i, tr in enumerate(p_stream):
        md_dict = {'max_coherence': max_coherence[i]}
        tr.stats.update(md_dict)
    # end for

    # TODO: Compute phase weighting vector per station per 2D (back_azimuth, distance) bin

    # Perform clustering for all traces in a station, and assign group IDs.
    # This will be super expensive when there are a lot of events, as the distance calculation grows as N^2.
    clustering_stream = p_stream.copy()
    clustering_stream = clustering_stream.trim2(-5.0, 25.0, 'onset')
    swipe = np.array([tr.data for tr in clustering_stream])
    if swipe.shape[0] > 1:
        ind = rf_group_by_similarity(swipe, similarity_eps)
    else:
        ind = np.array([0])
    # end if
    num_groups = np.amax(ind)
    logger.info("Station {}: detected {} clusters".format(station_id, num_groups))
    # Apply group
    for i, tr in enumerate(p_stream):
        md_dict = {'rf_group': ind[i] if ind[i] >= 0 else None}
        tr.stats.update(md_dict)
    # end for

    # TODO: Research techniques for grouping waveforms using singular value decomposition (SVD), possibly of
    # the complex waveform (from Hilbert transform) to determine the primary phase and amplitude components.
    # High similarity to the strongest eigenvectors might indicate waves in the primary group (group 0 in DBSCAN)
    # without the N^2 computational cost of DBSCAN.

    return (z_stream, p_stream, t_stream)
# end func


def rf_quality_metrics_queue(oqueue, station_id, station_stream3c, similarity_eps, drop_z=True):
    """Produce RF quality metrics in a stream and queue the QC'd components for
    downstream processing.

    :param oqueue: Output queue where filtered streams are queued
    :type oqueue: queue or multiprocessing.Manager.Queue
    :param station_id: Station ID
    :type station_id: str
    :param station_stream3c: 3-channel stream
    :type station_stream3c: list(rf.RFStream) with 3 components
    :param similarity_eps: Distance threshold used for DBSCAN clustering
    :type similarity_eps: float
    """
    streams_qual = compute_rf_quality_metrics(station_id, station_stream3c, similarity_eps)
    if streams_qual is not None:
        z_stream, p_stream, t_stream = streams_qual
        if drop_z:
            stream_qual = rf.RFStream([tr for doublet in zip(p_stream, t_stream)
                                       for tr in doublet])
        else:
            stream_qual = rf.RFStream([tr for triplet in zip(z_stream, p_stream, t_stream)
                                       for tr in triplet])
        # end if
        oqueue.put(stream_qual)
    # end if
# end func


@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(dir_okay=False), required=True)
@click.option('--temp-dir', type=click.Path(dir_okay=True), help="Temporary directory to use for best performance")
@click.option('--parallel/--no-parallel', default=True, show_default=True, help="Use parallel execution")
def main(input_file, output_file, temp_dir=None, parallel=True):
    """ This module filters RFs according to input options and then computes some quality metrics on each RF.
        This enables different downstream approaches to selecting and filtering for good quality RFs.

        The stats attribute of each RF is populated with these quality metrics. In addition, a new root group
        is added to the hdf5 file containing a Pandas DataFrame that tabulates the attributes of each trace
        to allow easy event filtering in the downstream workflow.

    Available methods:
    1. rf_group_by_similarity - grouping method based on calculation of euclidean distances and clustering by
       similarity ( aca machine learning approach)
    2. TODO: coherence - finding the coherent signals (in frequency domain) relative to median. Consequently, moveout
       should be applied to use this technique
    3. TODO knive - analysing the change of RMS relative to median. Noisy stations will give higher input. Moveout
       should be applied to use this technique
    4. S/N ratio
    5. Spectral entropy
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    similarity_eps = 0.05

    # Set up asynchronous buffered writing of results to file
    mgr = Manager()
    with h5py.File(input_file, mode='r') as h5f:
        config_str = h5f.attrs['metadata'] if 'metadata' in h5f.attrs else ''
    write_queue = mgr.Queue()
    output_thread = Process(target=async_write, args=(write_queue, output_file, 20, config_str))
    output_thread.daemon = True
    output_thread.start()

    logger.info("Processing source file {}".format(input_file))
    if parallel:
        logger.info("Parallel processing")
        Parallel(n_jobs=-3, verbose=5, max_nbytes='16M', temp_folder=temp_dir) \
            (delayed(rf_quality_metrics_queue)(write_queue, station_id, station_stream3c, similarity_eps)
             for station_id, station_stream3c in IterRfH5StationEvents(input_file))
    else:
        logger.info("Serial processing")
        for station_id, station_stream3c in IterRfH5StationEvents(input_file):
            try:
                rf_quality_metrics_queue(write_queue, station_id, station_stream3c, similarity_eps)
            except (ValueError, AssertionError) as e:
                traceback.print_exc()
                logger.error("Unhandled exception occurred in rf_quality_metrics_queue for station {}. "
                             "Data will be omitted for this station!\nError:\n{}".format(station_id, str(e)))
            # end try
        # end for
    # end if

    # Signal completion
    logger.info("Finishing...")
    write_queue.put(None)
    write_queue.join()

    logger.info("rf_quality_filter SUCCESS!")
# end func


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
