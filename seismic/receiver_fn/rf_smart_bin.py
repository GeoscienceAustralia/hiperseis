#!/usr/bin/env python

from builtins import range

import logging

import numpy as np
import click
from scipy.spatial.distance import euclidean
from fastdtw import fastdtw
import itertools as iter
from matplotlib.pyplot import plot, show, figure, ylim, xlabel, ylabel, legend, subplot2grid, GridSpec
from copy import deepcopy
from obspy.core.utcdatetime import UTCDateTime

from sklearn.cluster import DBSCAN

# Here are the libraries to deal with RFSTREAM, it uses obspy classes for event and station.
# from rf import RFStream
import rf
from obspy.core.event.event import Event
from obspy.core.inventory.station import Station
from joblib import Parallel, delayed
import matplotlib.pyplot as plt

from rf.profile import profile
# from tqdm import tqdm
from rf.imaging import plot_profile_map
from rf import get_profile_boxes, iter_event_data, IterMultipleComponents

logging.basicConfig()

# pylint: disable=invalid-name

def compare_pairs(data):
    distance, _ = fastdtw(data[0], data[1], dist=euclidean)
    return distance


def crossSpectrum(x, y):

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


def coh(y, y2):
    # This subroutine determines a coherence level between two signals on normilised frequency
    p11, freq = crossSpectrum(y, y)
    p22, freq = crossSpectrum(y2, y2)
    p12, freq = crossSpectrum(y, y2)
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
    distance = Parallel(n_jobs=-1, verbose=1)(map(delayed(compare_pairs), iter.combinations(swipe, 2)))
    index = list((i, j) for ((i, _), (j, _)) in iter.combinations(enumerate(swipe), 2))

    # First check that distance between points
    index = np.array(index)
    distance = np.array(distance)
    matrix = np.zeros((np.amax(index) + 1, 1 + np.amax(index))) + np.amax(distance)

    matrix[index[:, 0], index[:, 1]] = distance[:]
    clustering = DBSCAN(eps=similarity_eps, min_samples=2, metric='precomputed').fit(matrix)

    return clustering.labels_


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
        f, c = coh(median, swipe[i, :])
        if np.amax(c[f > f1 & f < f2]) > level:
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

#     for i in range(swipe.shape[0]):
#         ch.append(np.amax(coh(average, swipe[i, :])))
#         dev.append(np.sum((swipe[i, :]-average)**2)/(swipe.shape[0]-1))
#         ind[i] = False
#         knive.append(np.std(swipe[ind]))
#         sn.append(np.std(swipe[i, t > pulse_ind]) /
#                   np.std(swipe[i, t < pulse_ind]))

#     knive = np.array(knive)
#     sn = np.array(sn)
#     return knive < k_level, sn > sn_level


def remove_small_s2n(stream, ratio, teleseismic_cutout):
    """Filter out streams with low signal to noise (S/N) ratio.

    :param stream: [description]
    :type stream: [type]
    :param ratio: [description]
    :type ratio: [type]
    :param teleseismic_cutout: [description]
    :type teleseismic_cutout: [type]
    :return: [description]
    :rtype: [type]
    """
    # Take this segment as noise signal (TODO: document this more clearly)
    noise = stream.slice2(-5, -2, 'onset')
    # Take this segment as RF Ps (?) conversion signal (TODO: document this more clearly)
    signal = stream.slice2(-1, 2, 'onset')
    newstream = rf.RFStream()
    for i in range(len(stream)):
        rms = np.sqrt(
            np.mean(np.square(signal[i].data)))/np.sqrt(np.mean(np.square(noise[i].data)))
        # Why do we confound the teleseismic distance filtering here with S/N filtering?
        if rms > ratio and stream[i].stats.distance > teleseismic_cutout:
            newstream.append(stream[i])
    return newstream


def get_rf_stream_components(stream):

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
    """ @package rf_smart_bin
    This code contains different approaches to select good quality RFs.
    Currently there are three methods
    1. rf_group_by_similarity - grouping method based on calculation of euclidean distances and clustering by
       similarity ( aca machine learning approach)
    2. coherence - finding the coherent signals (in frequency domain) relative to median. Consequently, moveout
       should be applied to use this technique
    3. knive - analysing the change of RMS relative to median. Noisy stations will give higher input. Moveout
       should be applied to use this technique
    """
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    # Settings
#    similarity_eps = 3.0
    similarity_eps = 8.0

    min_snr = 1.5
    # teleseismic cut out (35 degrees) for remove_small_s2n
    min_teleseismic_dist = 35.0

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
    prim_stream = remove_small_s2n(prim_stream, min_snr, min_teleseismic_dist)
    logger.info("Number of traces after S/N cut out is: %d", len(prim_stream))


    # we have to decimate here otherwise clustering method wouldn't perform well. 5Hz sampling
    rf_stream = prim_stream.copy()
    # Filter specified below is only for data analysis and not applied to output data
    rf_stream = rf_stream.filter('bandpass', freqmin=0.05, freqmax=0.7).interpolate(sampling_rate=5.0)

    # original file will be interpolated to 100Hz
    prim_stream = prim_stream.trim2(starttime=-5, endtime=60, reftime='onset')

    # here we collect station names but maybe ID is more appropriate in case of having the same
    # station names in different deployments
    station_list = []
    for i in range(len(rf_stream)):
        station_list.append(rf_stream[i].stats.station.encode('utf-8'))

    station_list = np.unique(np.array(station_list))
    print("Gathered ", len(station_list), " stations")

    # here we go with the main loop over stations
    out_file = rf.RFStream()

    for i in range(station_list.shape[0]):
        station_code = station_list[i].decode('utf-8')
        print("Station ", station_code, i+1, " of ", station_list.shape[0])
        traces = rf_stream.select(station=station_code).copy()

        # we choose short RF to simplify and speed up the processing
        traces = traces.trim2(-5, 20, 'onset')

        # but keep original traces as they are to use them at the end
        original_traces = prim_stream.select(station=station_code)

        swipe = []
        original_swipe = []

        for trace in traces:
            swipe.append(trace.data)
        for trace in original_traces:
            original_swipe.append(trace.data)

        swipe = np.array(swipe)
        original_swipe = np.array(original_swipe)

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
            for j in range(len(original_traces)):
                if ind[j] == k:
                    original_traces[j].stats.update(rf_group)
                    out_file.append(original_traces[j])
                    for tr in t_stream:
                        if tr.stats.event_id == original_traces[j].stats.event_id:
                            tr.stats.update(rf_group)
                            out_file.append(tr)
                    # end for
                    for tr in z_stream:
                        if tr.stats.event_id == original_traces[j].stats.event_id:
                            tr.stats.update(rf_group)
                            out_file.append(tr)
                    # end for
            # end for

            # # or here we can make a trick - coherent signal comes from different directions or segments.
            # # Therefore we assign stacked RF back to its original azimuths and angle of incidence.
            # for j in range(len(original_traces)):
            #     if ind[j]==k:
            #         # here we replace original data by stacked rays. However original RFs with assigned
            #         # groups can be used as well and stacked later using migration image
            #         # this option can be more favourable to highlight small signals.
            #         # Comment out one line below to avoid stacking
            #         original_traces[j].data=stacked.copy()
            #         out_file.append(original_traces[j])

        # end for
    # end for

    # # Some plots if required
    # ppoints = out_file.ppoints(70)
    # boxes = get_profile_boxes((-18.4, 139.1), 135, np.linspace(0, 440, 80), width=500)
    # pstream = profile(out_file, boxes)
    # pstream.plot_profile(scale=1.5,top='hist')
    # plt.show()

    # Output file
    out_file.write(output_file, 'H5')
# end func


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
