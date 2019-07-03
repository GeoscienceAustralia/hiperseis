#!/usr/bin/env python

from builtins import range

import logging

import itertools as iter
from copy import deepcopy

import numpy as np
import click
from scipy.spatial.distance import euclidean
from scipy import signal
from fastdtw import fastdtw
# from matplotlib.pyplot import plot, show, figure, ylim, xlabel, ylabel, legend, subplot2grid, GridSpec
# import matplotlib.pyplot as plt

from sklearn.cluster import DBSCAN

import rf
from joblib import Parallel, delayed

# from tqdm import tqdm

logging.basicConfig()

# pylint: disable=invalid-name


def _crossSpectrum(x, y):

    # -------------------Remove mean-------------------
    # nperseg chosen arbitrary based on 126 samples RF signal, experiment to get best results
    nperseg = x.size/20
    cross = np.zeros(nperseg, dtype='complex128')
    for ind in range(x.size / nperseg):

        xp = x[ind * nperseg: (ind + 1)*nperseg]
        yp = y[ind * nperseg: (ind + 1)*nperseg]
        xp = xp - np.mean(xp)
        yp = yp - np.mean(xp)

    # Do FFT
        cfx = np.fft.fft(xp)
        cfy = np.fft.fft(yp)

    # Get cross spectrum
        cross += cfx.conj()*cfy
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
    coh = np.divide(part1, p22.real, out=np.zeros_like(
        part1), where=p22.real != 0)

#   plot( freq[freq > 0], coh[freq > 0])
#   show()
#   return coh[freq > 0]

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
    # Helper function used to compute similarity distance between a pair of waveforms.
    # See https://pypi.org/project/fastdtw/
    def _compare_pairs(data):
        distance, _ = fastdtw(data[0], data[1], dist=euclidean)
        return distance
    # end func

    # Pre-compute distance metric that will be used for DBSCAN clustering.
    distance = Parallel(n_jobs=-1, verbose=1)(map(delayed(_compare_pairs), iter.combinations(swipe, 2)))
    index = list((i, j) for ((i, _), (j, _)) in iter.combinations(enumerate(swipe), 2))

    # First check that distance between points
    index = np.array(index)
    distance = np.array(distance)
    matrix = np.zeros((np.amax(index) + 1, 1 + np.amax(index))) + np.amax(distance)

    matrix[index[:, 0], index[:, 1]] = distance[:]
    clustering = DBSCAN(eps=similarity_eps, min_samples=2, metric='precomputed').fit_predict(matrix)

    return clustering


def coherence(swipe, level, f1, f2):
    """ Finding coherence between two signals in frequency domain
        swipe - matrix with  waveforms orginised rowwise
        level  - minimum level of coherence (>0.6) for good results
        f1 and f2 - normalised min and max frequencies
        returns array of indexes for coherent traces with median
    """
    # level - minimum coherence > 0.6 for good results, f2 <0.5 for RF
    median = np.median(swipe, axis=0)

    index = []
    for i in range(swipe.shape[0]):
        f, c = _coh(median, swipe[i, :])
        if np.amax(c[(f > f1) & (f < f2)]) > level:
            index.append(True)
        else:
            index.append(False)

    return np.array(index)


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


def remove_small_s2n(stream, min_snr):
    """Filter out streams with low signal to noise (S/N) ratio.

    :param stream: [description]
    :type stream: [type]
    :param ratio: [description]
    :type ratio: [type]
    :return: [description]
    :rtype: [type]
    """
    # Take this segment as noise signal (TODO: document this more clearly)
    noise = stream.slice2(-5, -2, 'onset')
    noise = np.array([tr.data for tr in noise])
    # Take this segment as RF Ps (?) conversion signal (TODO: document this more clearly)
    pick_signal = stream.slice2(-1, 2, 'onset')
    pick_signal = np.array([tr.data for tr in pick_signal])
    assert pick_signal.shape[0] == noise.shape[0]

    # Compute RMS vectorized across whole dataset
    rms = np.sqrt(np.mean(np.square(pick_signal), axis=1)) / np.sqrt(np.mean(np.square(noise), axis=1))
    # Augmenting with RMS of envelope (complex amplitude) rather than pure
    # wave amplitude, to capture the energy in the time rate of change.
    pick_env = np.absolute(signal.hilbert(pick_signal, axis=1))
    noise_env = np.absolute(signal.hilbert(noise, axis=1))
    rms_env = np.sqrt(np.mean(np.square(pick_env), axis=1)) / np.sqrt(np.mean(np.square(noise_env), axis=1))

    rms_mask = (rms >= min_snr) | (rms_env >= min_snr)
    newstream = rf.RFStream([s for i, s in enumerate(stream) if rms_mask[i]])

    return newstream


def spectral_entropy(trace):
    """Compute the spectral entropy of a trace

    :param trace: Single channel seismic trace
    :type trace: rf.RFTrace
    :return: Spectral entropy of the trace waveform
    :rtype: float
    """
    fs = trace.stats.sampling_rate
    _, psd = signal.periodogram(trace.data, fs, detrend='linear')
    psd = psd/np.sum(psd)
    entropy = -np.sum(psd*np.log(psd))
    return entropy


def get_rf_stream_components(stream):
    """Identify the RF component types and return them.

    :param stream: [description]
    :type stream: rf.RFStream
    :return: (RF component type, primary component, alt component 1, alt component 2)
    :rtype: (str, rf.RFStream, rf.RFStream, rf.RFStream)
    """

    DEFAULT_RF_TYPE = 'LQT-Q'
    rf_type = DEFAULT_RF_TYPE

    primary_stream = stream.select(component='Q')

    if primary_stream:
        sec_stream = stream.select(component='T')
        tert_stream = stream.select(component='L')
    else:
        rf_type = 'ZRT-R'
        primary_stream = stream.select(component='R')
        if primary_stream:
            sec_stream = stream.select(component='T')
            tert_stream = stream.select(component='Z')
        else:
            return None, None, None, None
        # end if
    # end if
    return rf_type, primary_stream, sec_stream, tert_stream


# -------------Main---------------------------------

@click.command()
@click.argument('input-file', type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument('output-file', type=click.Path(dir_okay=False), required=True)
def main(input_file, output_file):
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

    # Settings
    similarity_eps = 3.0
    # similarity_eps = 8.0

    min_snr = 5.0

    print("Reading the input file...")
    # Input file
    stream = rf.read_rf(input_file, 'H5')
    print("Reading is done...")
    # output file naming --> look at the end of the code

    # Pull out streams for analysis. FIXME: Variable naming here seems to assume the type will be 'ZRT-R'
    rf_type, prim_stream, t_stream, z_stream = get_rf_stream_components(stream)
    if rf_type is None:
        logger.error("Could not identify RF components from input stream")
        return
    # end if

    # first lets just remove plainly bad data
    logger.info("Number of traces before S/N cut out is: %d", len(prim_stream))
    prim_stream = remove_small_s2n(prim_stream, min_snr)
    logger.info("Number of traces after S/N cut out is: %d", len(prim_stream))


    # we have to decimate here otherwise clustering method wouldn't perform well. 5Hz sampling
    rf_stream = prim_stream.copy()
    # Filter specified below is only for data analysis and not applied to output data
    rf_stream = rf_stream.filter('bandpass', freqmin=0.05, freqmax=0.7).interpolate(sampling_rate=5.0)

    # original file will be interpolated to 100Hz
    prim_stream = prim_stream.trim2(starttime=-5, endtime=60, reftime='onset')

    # here we collect station names but maybe ID is more appropriate in case of having the same
    # station names in different deployments
    station_list = set()
    for stn_stream in prim_stream:
        station_list.add(stn_stream.stats.station.encode('utf-8'))

    print("Gathered ", len(station_list), " stations")

    # here we go with the main loop over stations
    filtered_stream = rf.RFStream()

    for i, stn_code in enumerate(station_list):
        station_code = stn_code.decode('utf-8')
        print("Station ", station_code, " ", i + 1, " of ", len(station_list))
        stn_stream = rf_stream.select(station=station_code).copy()

        # we choose short RF to simplify and speed up the processing
        stn_stream = stn_stream.trim2(-5, 20, 'onset')

        # but keep original traces as they are to use them at the end
        orig_p_stream = prim_stream.select(station=station_code)

        swipe = []
        # original_swipe = []

        for trace in stn_stream:
            swipe.append(trace.data)
        # for trace in orig_p_stream:
        #     original_swipe.append(trace.data)

        swipe = np.array(swipe)
        # original_swipe = np.array(original_swipe)

        print("Processing ", swipe.shape[0], " events")
        # we use clustering technique to find similar signals
        if swipe.shape[0] > 1:
            ind = rf_group_by_similarity(swipe, similarity_eps)
        else:
            ind = np.array([0])
            print(station_code, ' has only one trace')

        num_group = np.amax(ind)
        print("Number of detected groups: ", num_group + 1)

        # we have group indexes for each good quality RF trace and apply grouping to original RF traces for stacking
        for k in range(num_group + 1):
            # average can use weights and mean can work on masked arrays
            # stacked = np.average(original_swipe[ind == k, :], axis=0)
            # we choose only traces that belong to detected group
            rf_group = {'rf_group': k}
            for j, orig_trace in enumerate(orig_p_stream):
                if ind[j] == k:
                    orig_trace.stats.update(rf_group)
                    filtered_stream.append(orig_trace)
                    for tr in t_stream:
                        if tr.stats.event_id == orig_trace.stats.event_id:
                            tr.stats.update(rf_group)
                            filtered_stream.append(tr)
                    # end for
                    for tr in z_stream:
                        if tr.stats.event_id == orig_trace.stats.event_id:
                            tr.stats.update(rf_group)
                            filtered_stream.append(tr)
                    # end for
            # end for

            # # or here we can make a trick - coherent signal comes from different directions or segments.
            # # Therefore we assign stacked RF back to its original azimuths and angle of incidence.
            # for j in range(len(orig_p_stream)):
            #     if ind[j]==k:
            #         # here we replace original data by stacked rays. However original RFs with assigned
            #         # groups can be used as well and stacked later using migration image
            #         # this option can be more favourable to highlight small signals.
            #         # Comment out one line below to avoid stacking
            #         orig_p_stream[j].data=stacked.copy()
            #         filtered_stream.append(orig_p_stream[j])

        # end for
    # end for

    # # Some plots if required
    # ppoints = filtered_stream.ppoints(70)
    # boxes = rf.get_profile_boxes((-18.4, 139.1), 135, np.linspace(0, 440, 80), width=500)
    # pstream = rf.profile(filtered_stream, boxes)
    # pstream.plot_profile(scale=1.5,top='hist')
    # plt.show()

    # Output file
    filtered_stream.write(output_file, 'H5')
# end func


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
