#!/usr/bin/env python

from builtins import range

import logging
import itertools
import traceback
from multiprocessing import Process, Manager

import numpy as np
import pandas as pd
import click

from scipy import signal
from scipy import stats
# from fastdtw import fastdtw
# from matplotlib.pyplot import plot, show, figure, ylim, xlabel, ylabel, legend, subplot2grid, GridSpec
# import matplotlib.pyplot as plt

from sklearn.cluster import DBSCAN
from joblib import Parallel, delayed
import obspy
import rf

from seismic.receiver_fn.rf_process_io import async_write
from seismic.receiver_fn.rf_h5_file_station_iterator import IterRfH5StationEvents
from seismic.receiver_fn.rf_util import compute_rf_snr

# from tqdm import tqdm

logging.basicConfig()

# pylint: disable=invalid-name, logging-format-interpolation

MIN_RESAMPLING_RATE_HZ = 8.0


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


def rf_group_by_similarity(swipe, similarity_eps, temp_dir=None):
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

    distance = Parallel(n_jobs=-3, verbose=5, max_nbytes='16M', temp_folder=temp_dir)(
        map(delayed(_compare_pairs_rmsdist), itertools.combinations(swipe, 2)))
    # distance = list(map(_compare_pairs_rmsdist, itertools.combinations(swipe, 2)))
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
        swipe - matrix with  waveforms organised rowwise
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
    # Taper aggressively to ensure that end discontinuities do not inject spectral energy at high frequencies
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

    :param stream: [description]
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


def rf_quality_metrics(oqueue, station_id, station_stream3c, similarity_eps, temp_dir=None, filter_band=None):

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
        return pd.DataFrame()
    # end if

    # Flatten the traces into a single RFStream for subsequent processing
    rf_streams = rf.RFStream([tr for stream in nonan_streams for tr in stream if tr.stats.type == 'rf'])
    raw_streams = obspy.core.stream.Stream(
        [tr for stream in nonan_streams for tr in stream if tr.stats.type == 'raw_resampled'])

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
        return pd.DataFrame()
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

    # Extract RF type and the primary polarized component (ignore other component streams)
    rf_type, p_stream, _, _ = get_rf_stream_components(streams)
    if rf_type is None:
        logger.error("Unrecognized RF type for station {}. File might not be RF file!".format(station_id))
        return None
    # end if

    # Filter p_stream prior to computing quality metrics and downsample to Nyquist freq (on the assumption that the
    # highest freq component is ~2*freq_max).
    if filter_band is not None:
        freq_min = filter_band[0]
        freq_max = filter_band[1]
        p_stream = p_stream.filter('bandpass', freqmin=freq_min, freqmax=freq_max, corners=2, zerophase=True)
        sampling_rate = max(4*freq_max, MIN_RESAMPLING_RATE_HZ)
        if sampling_rate < p_stream[0].stats.sampling_rate:
            p_stream = p_stream.interpolate(sampling_rate=sampling_rate, method='lanczos', a=20)
        # end if
    # end if

    # Compute S/N ratios for RFs
    compute_rf_snr(p_stream)

    # Compute spectral entropy for all traces
    sp_entropy = spectral_entropy(p_stream)
    for i, tr in enumerate(p_stream):
        md_dict = {'entropy': sp_entropy[i]}
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
        ind = rf_group_by_similarity(swipe, similarity_eps, temp_dir=temp_dir)
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
    # High similarity to the strongest eigenvectors indicates waves in the primary group (group 0 in DBSCAN)
    # without the N^2 computational cost.

    # Output resulting station streams. Only keeping the primary RF stream since the others are not needed
    # for RF analysis. Merge RF streams with raw waveforms for those events that made it through filtering.
    event_ids = {tr.stats.event_id: tr for tr in p_stream}
    for tr in raw_streams:
        if tr.stats.event_id in event_ids:
            # If it is the Z component, copy its SNR to the RF stream for convenient access later.
            if tr.stats.channel[-1].upper() == 'Z':
                if 'snr_prior' in tr.stats:
                    md_dict = {'snr_prior': tr.stats['snr_prior']}
                else:
                    md_dict = {'snr_prior': np.nan}
                # end if
                p_trace = event_ids[tr.stats.event_id]
                p_trace.stats.update(md_dict)
            # end if
            p_stream.append(tr)
        # end if
    # end for

    oqueue.put(p_stream)

    # Return Pandas DataFrame of tabulated metadata for this station
    return pd.DataFrame()


# -------------Main---------------------------------

@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(dir_okay=False), required=True)
@click.option('--temp-dir', type=click.Path(dir_okay=True), help="Temporary directory to use for best performance")
@click.option('--filter-band', type=(float, float), help="Filter pass band edges (Hz). Leave empty for no filtering.")
def main(input_file, output_file, temp_dir=None, filter_band=None):
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
    write_queue = mgr.Queue()
    output_thread = Process(target=async_write, args=(write_queue, output_file, 10))
    output_thread.daemon = True
    output_thread.start()

    logger.info("Processing source file {}".format(input_file))
    # Process in serial, delegate parallelization to station data processor function
    metadata_list = []
    for station_id, station_stream3c in IterRfH5StationEvents(input_file):
        try:
            md_table = rf_quality_metrics(write_queue, station_id, station_stream3c, similarity_eps, temp_dir=temp_dir,
                                          filter_band=filter_band)
            metadata_list.append(md_table)
        except (ValueError, AssertionError) as e:
            traceback.print_exc()
            logger.error("Unhandled exception occurred in rf_quality_metrics for station {}. "
                         "Data will be omitted for this station!\nError:\n{}".format(station_id, str(e)))
        # end try
    # end for

    # TODO: Filter out None items before concat, if required.
    all_metadata = pd.concat(metadata_list)

    # Signal completion
    logger.info("Finishing...")
    write_queue.put(None)
    write_queue.join()

    # TODO: Write Pandas DataFrame of metadata to output file

    logger.info("rf_quality_filter SUCCESS!")
# end func


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
