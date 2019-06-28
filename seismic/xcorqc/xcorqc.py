#!/usr/bin/env python
"""
Description:
    Cross-correlation functionality
References:

CreationDate:   29/06/17
Developer:      laurence.davies@ga.gov.au

Revision History:
    LastUpdate:     29/06/17   LD       First commit of xcor code.
    LastUpdate:     13/07/17   LD       Fixed xcor filtering issue when traces have different sample rates.
    LastUpdate:     11/08/17   RH       Implement ASDF-based cross-correlation workflow
    LastUpdate:     11/07/18   RH       Implemented parallel cross-correlator
    LastUpdate:     19/07/18   RH       Implemented cross-correlation approaches described in Habel et al. 2018

    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
import datetime
import glob
from os.path import join, exists
import logging
import math
from collections import defaultdict

import numpy as np
import scipy
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

import obspy
from obspy.core import Stream, UTCDateTime, Stats
from obspy import read, Trace
from obspy.signal.cross_correlation import xcorr
from obspy.signal.detrend import simple, spline
from obspy.signal.filter import bandpass
from scipy import signal

from seismic.xcorqc.fft import *
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from netCDF4 import Dataset
from functools import reduce

logging.basicConfig()


def setup_logger(name, log_file, level=logging.INFO):
    """
    Function to setup a logger; adapted from stackoverflow
    """
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler = logging.FileHandler(log_file, mode='w')
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger
# end func

def zeropad(tr, padlen):
    assert (tr.shape[0] < padlen)
    padded = np.zeros(padlen)
    padded[0:tr.shape[0]] = tr
    return padded
# end func

def zeropad_ba(tr, padlen):
    assert (tr.shape[0] < padlen)
    padded = np.zeros(padlen, dtype=np.complex_)
    s = int((padlen - tr.shape[0]) / 2)
    padded[s:(s + tr.shape[0])] = scipy.fftpack.fftshift(tr)
    return scipy.fftpack.ifftshift(padded)
# end func

def taper(tr, taperlen):
    tr[0:taperlen] *= 0.5 * (1 + np.cos(np.linspace(-math.pi, 0, taperlen)))
    tr[-taperlen:] *= 0.5 * (1 + np.cos(np.linspace(0, math.pi, taperlen)))
    return tr
# end func

def whiten(x, sr, fmin=None, fmax=None):

    if not fmin: fmin = 0
    if not fmax: fmax = sr/2.

    # apply hanning window
    l = x.shape[0]
    h = signal.hann(l)
    xf = rfft(h*x)

    # magnitude
    mag = np.abs(xf)

    # phase
    phase = np.unwrap(np.angle(xf))

    f = np.fft.rfftfreq(l, 1./float(sr))
    j = np.where((f >= fmin) & (f <= fmax))[0]

    mag[j] = 1.
    xnew = irfft(mag*np.exp(1j*phase))

    return xnew
# end func

def xcorr2(tr1, tr2, sta1_inv=None, sta2_inv=None,
           window_seconds=3600, window_overlap=0.1, interval_seconds=86400,
           resample_rate=None, flo=None, fhi=None,
           clip_to_2std=False, whitening=False,
           one_bit_normalize=False, envelope_normalize=False,
           verbose=1, logger=None):

    assert window_overlap <= 0.5, 'Overlap must be < 0.5'

    sr1 = tr1.stats.sampling_rate
    sr2 = tr2.stats.sampling_rate
    sr1_orig = sr1
    sr2_orig = sr2
    tr1_d_all = tr1.data  # refstn
    tr2_d_all = tr2.data
    lentr1_all = tr1_d_all.shape[0]
    lentr2_all = tr2_d_all.shape[0]
    window_samples_1 = window_seconds * sr1
    window_samples_2 = window_seconds * sr2
    interval_samples_1 = interval_seconds * sr1
    interval_samples_2 = interval_seconds * sr2
    itr1s = 0
    itr2s = 0
    sr = 0
    resll = []

    if(resample_rate):
        sr1 = resample_rate
        sr2 = resample_rate
    # end if
    sr = max(sr1, sr2)
    xcorlen = int(2 * window_seconds * sr - 1)
    fftlen = 2 ** (int(np.log2(xcorlen)) + 1)

    intervalCount = 0
    windowsPerInterval = []  # Stores the number of windows processed per interval
    intervalStartSeconds = []
    intervalEndSeconds = []
    while itr1s < lentr1_all and itr2s < lentr2_all:
        itr1e = min(lentr1_all, itr1s + interval_samples_1)
        itr2e = min(lentr2_all, itr2s + interval_samples_2)

        windowCount = 0
        wtr1s = int(itr1s)
        wtr2s = int(itr2s)
        resl = []

        while wtr1s < itr1e and wtr2s < itr2e:
            wtr1e = int(min(itr1e, wtr1s + window_samples_1))
            wtr2e = int(min(itr2e, wtr2s + window_samples_2))

            # Discard small windows
            if wtr1e - wtr1s < window_samples_1 or wtr2e - wtr2s < window_samples_2:
                wtr1s = wtr1e
                wtr2s = wtr2e
                continue
            # end if

            # Discard windows with masked regions, i.e. with gaps
            if (not (np.ma.is_masked(tr1_d_all[wtr1s:wtr1e])
                     or np.ma.is_masked(tr2_d_all[wtr2s:wtr2e]))):

                tr1_d = np.array(tr1_d_all[wtr1s:wtr1e], dtype=np.float32)
                tr2_d = np.array(tr2_d_all[wtr2s:wtr2e], dtype=np.float32)

                # resample
                if(resample_rate):
                    tr1_d = Trace(data=tr1_d,
                                  header=Stats(header={'sampling_rate':sr1_orig,
                                                       'npts':window_samples_1})).resample(resample_rate).data
                    tr2_d = Trace(data=tr2_d,
                                  header=Stats(header={'sampling_rate':sr2_orig,
                                                       'npts':window_samples_2})).resample(resample_rate).data
                # end if

                # remove response
                if(sta1_inv):
                    tr1 = Trace(data=tr1_d,
                                header=Stats(header={'sampling_rate':resample_rate if resample_rate else sr1_orig,
                                                     'npts':len(tr1_d),
                                                     'network':tr1.stats.network,
                                                     'station': tr1.stats.station,
                                                     'location': tr1.stats.location,
                                                     'channel': tr1.stats.channel,
                                                     'starttime': tr1.stats.starttime + wtr1s,
                                                     'endtime': tr1.stats.starttime + wtr1e}))
                    try:
                        tr1.remove_response(inventory=sta1_inv, water_level=50)
                    except Exception as e:
                        print (e)
                    # end try
                    tr1_d = tr1.data
                # end if

                if(sta2_inv):
                    tr2 = Trace(data=tr2_d,
                                header=Stats(header={'sampling_rate':resample_rate if resample_rate else sr2_orig,
                                                     'npts':len(tr2_d),
                                                     'network': tr2.stats.network,
                                                     'station': tr2.stats.station,
                                                     'location': tr2.stats.location,
                                                     'channel': tr2.stats.channel,
                                                     'starttime': tr2.stats.starttime + wtr2s,
                                                     'endtime': tr2.stats.starttime + wtr2e}))
                    try:
                        tr2.remove_response(inventory=sta2_inv, water_level=50)
                    except Exception as e:
                        print (e)
                    # end try
                    tr2_d = tr2.data
                # end if

                # detrend
                tr1_d = spline(tr1_d, 2, 1000)
                tr2_d = spline(tr2_d, 2, 1000)

                # zero-mean
                tr1_d -= np.mean(tr1_d)
                tr2_d -= np.mean(tr2_d)

                # taper
                if(window_overlap>0):
                    tr1_d = taper(tr1_d, int(window_overlap*tr1_d.shape[0]))
                    tr2_d = taper(tr2_d, int(window_overlap*tr2_d.shape[0]))
                # end if

                # apply zero-phase band-pass
                if(flo and fhi):
                    tr1_d = bandpass(tr1_d, flo, fhi, sr1, corners=6, zerophase=True)
                    tr2_d = bandpass(tr2_d, flo, fhi, sr2, corners=6, zerophase=True)
                # end if

                # clip to +/- 2*std
                if(clip_to_2std):
                    std_tr1 = np.std(tr1_d)
                    std_tr2 = np.std(tr2_d)
                    clip_indices_tr1 = np.fabs(tr1_d) > 2 * std_tr1
                    clip_indices_tr2 = np.fabs(tr2_d) > 2 * std_tr2

                    tr1_d[clip_indices_tr1] = 2 * std_tr1 * np.sign(tr1_d[clip_indices_tr1])
                    tr2_d[clip_indices_tr2] = 2 * std_tr2 * np.sign(tr2_d[clip_indices_tr2])
                # end if

                # spectral whitening
                if(whitening):
                    tr1_d = whiten(tr1_d, sr1)
                    tr2_d = whiten(tr2_d, sr2)
                # end if

                # 1-bit normalization
                if(one_bit_normalize):
                    tr1_d = np.sign(tr1_d)
                    tr2_d = np.sign(tr2_d)
                # end if

                # Apply Rhys Hawkins-style default time domain normalization
                if (clip_to_2std == 0 and one_bit_normalize == 0):
                    # 0-mean
                    tr1_d -= np.mean(tr1_d)
                    tr2_d -= np.mean(tr2_d)

                    # unit-std
                    tr1_d /= np.std(tr1_d)
                    tr2_d /= np.std(tr2_d)
                # end if

                if (sr1 < sr2):
                    fftlen2 = fftlen
                    fftlen1 = int((fftlen2 * 1.0 * sr1) / sr)
                    rf = zeropad_ba(fftn(zeropad(tr1_d, fftlen1), shape=[fftlen1]), fftlen2) * fftn(
                        zeropad(ndflip(tr2_d), fftlen2), shape=[fftlen2])
                elif (sr1 > sr2):
                    fftlen1 = fftlen
                    fftlen2 = int((fftlen1 * 1.0 * sr2) / sr)
                    rf = fftn(zeropad(tr1_d, fftlen1), shape=[fftlen1]) * zeropad_ba(
                        fftn(zeropad(ndflip(tr2_d), fftlen2), shape=[fftlen2]), fftlen1)
                else:
                    rf = fftn(zeropad(tr1_d, fftlen), shape=[fftlen]) * fftn(zeropad(ndflip(tr2_d), fftlen),
                                                                             shape=[fftlen])
                # end if

                resl.append(rf)
                windowCount += 1
            # end if

            wtr1s += int(window_samples_1 - window_samples_1 * window_overlap)
            wtr2s += int(window_samples_2 - window_samples_2 * window_overlap)
        # end while (windows within interval)

        if (verbose > 1):
            if(logger): logger.info('\tProcessed %d windows in interval %d' % (windowCount, intervalCount))
        # end fi

        if(np.fabs(itr1e - itr1s) < sr1_orig):
            itr1s = itr1e
            itr2s = itr2e
            continue
        # end if

        intervalStartSeconds.append(itr1s/sr1_orig)
        intervalEndSeconds.append(itr1e/sr1_orig)
        itr1s = itr1e
        itr2s = itr2e
        intervalCount += 1

        # Append an array of zeros if no windows were processed for the current interval
        if (windowCount == 0):
            resl.append(np.zeros(fftlen))
            if (verbose == 1):
                if(logger): logger.info('\tWarning: No windows processed due to gaps in data in current interval')
            # end if
        # end if

        windowsPerInterval.append(windowCount)

        mean = reduce((lambda tx, ty: tx + ty), resl) / len(resl)

        if (envelope_normalize):
            step = np.sign(np.fft.fftfreq(fftlen, 1.0 / sr))
            mean = mean + step * mean  # compute analytic
        # end if

        mean = ifftn(mean)

        if(envelope_normalize):
            # Compute magnitude of mean
            mean = np.abs(mean)
            normFactor = np.max(mean)

            # mean can be 0 for a null result
            if(normFactor > 0):
                mean /= normFactor
            # end if
        # end if

        resll.append(mean[:xcorlen])
    # end while (iteration over intervals)

    if (len(resll)):
        return np.array(resll), np.array(windowsPerInterval), \
               np.array(intervalStartSeconds, dtype='i8'), \
               np.array(intervalEndSeconds, dtype='i8')
    else:
        return None, None, None, None
    # end if
# end func

def dropBogusTraces(st, sampling_rate_cutoff=1):
    badTraces = [tr for tr in st if tr.stats.sampling_rate < sampling_rate_cutoff]

    for tr in badTraces: st.remove(tr)
# end func

def IntervalStackXCorr(refds, tempds,
                       start_time, end_time,
                       ref_net_sta, temp_net_sta,
                       ref_sta_inv, temp_sta_inv,
                       ref_cha,
                       temp_cha,
                       resample_rate=None,
                       buffer_seconds=864000, interval_seconds=86400,
                       window_seconds=3600, flo=None, fhi=None,
                       clip_to_2std=False, whitening=False,
                       one_bit_normalize=False, envelope_normalize=False,
                       ensemble_stack=False,
                       outputPath='/tmp', verbose=1, tracking_tag=''):
    """
    This function rolls through two ASDF data sets, over a given time-range and cross-correlates
    waveforms from all possible station-pairs from the two data sets. To allow efficient, random
    data access asdf data sources, an instance of a SeisDB object, instantiated from
    the corresponding Json database is passed in (tempds_db) -- although this parameter is not
    mandatory, data-access from large ASDF files will be slow without it.

    Station-ids to be processed from the two data-sources can be specified as lists of strings,
    while wildcards can be used to process all stations. Data is fetched from the sources in chunks
    to limit memory usage and data-windows with gaps are discarded.

    Cross-correlation results are written out for each station-pair, in the specified folder, as
    NETCDF4 files. Panoply (https://www.giss.nasa.gov/tools/panoply/), already installed on the
    NCI VDIs can be used to interrogate these results.

    :type refds: FederatedASDFDataSet
    :param refds: FederatedASDFDataSet containing reference-station data
    :type tempds: FederatedASDFDataSet
    :param tempds: FederatedASDFDataSet containing temporary-stations data
    :type start_time: UTCDateTime
    :param: start_time: Start-time (UTCDateTime format) for data to be used in cross-correlation
    :type end_time: UTCDateTime
    :param: end_time: End-time (UTCDateTime format) for data to be used in cross-correlation
    :type ref_net_sta: str
    :param ref_net_sta: Network.Station for the reference Dataset.
    :type temp_net_sta: str
    :param temp_net_sta: Network.Station for the temporary Dataset.
    :type ref_sta_inv: Inventory
    :param ref_sta_inv: Inventory containing instrument response for station
    :type temp_sta_inv: Inventory
    :param temp_sta_inv: Inventory containing instrument response for station
    :type ref_cha: str
    :param ref_cha: Channel name for the reference Dataset
    :type temp_cha: str
    :param temp_cha: Channel name for the temporary Dataset
    :type resample_rate: int
    :param resample_rate: Resampling rate (Hz). Applies to both data-sets
    :type buffer_seconds: int
    :param buffer_seconds: The amount of data to be fetched per call from the ASDFDataSets, because
                           we may not be able to fetch all the data (from start_time to end_time) at
                           once. The default is set to 10 days and should be a multiple of
                           interval_seconds.
    :type interval_seconds: int
    :param interval_seconds: The interval in seconds, over which cross-correlation windows are
                             stacked. Default is 1 day.
    :type window_seconds: int
    :param window_seconds: Length of cross-correlation window in seconds. Default is 1 hr.
    :type flo: float
    :param flo: Lower frequency for Butterworth bandpass filter
    :type fhi: float
    :param fhi: Upper frequency for Butterworth bandpass filter
    :type clip_to_2std: bool
    :param clip_to_2std: Clip data in each window to +/- 2 standard deviations
    :type whitening: bool
    :param whitening: Apply spectral whitening
    :type one_bit_normalize: bool
    :param one_bit_normalize: Apply one-bit normalization to data in each window
    :type envelope_normalize: bool
    :param envelope_normalize: Envelope via Hilbert transforms and normalize
    :type ensemble_stack: bool
    :param ensemble_stack: Outputs a single CC function stacked over all data for a given station-pair
    :type verbose: int
    :param verbose: Verbosity of printouts. Default is 1; maximum is 3.
    :type tracking_tag: str
    :param tracking_tag: File tag to be added to output file names so runtime settings can be tracked
    :type outputPath: str
    :param outputPath: Folder to write results to
    :return: 1: 1d np.array with time samples spanning [-window_samples:window_samples]
             2: A dictionary of 2d np.arrays containing cross-correlation results for each station-pair.
                Rows in each 2d array represent number of interval_seconds processed and columns
                represent stacked samples of length window_seconds.
             3: A dictionary of 1d np.arrays containing number of windows processed, within each
                interval_seconds period, for each station-pair. These Window-counts could be helpful
                in assessing robustness of results.
    """
    #######################################
    # check consistency of parameterization
    #######################################
    if(whitening and (ref_sta_inv or temp_sta_inv)):
        raise RuntimeError('Mutually exclusive parameterization: specify either spectral whitening or '
                           'instrument response removal')
    # end if

    if(clip_to_2std and one_bit_normalize):
        raise RuntimeError('Mutually exclusive parameterization: clip_to_2std and one-bit-normalizations'
                           'together is redundant')
    # end if

    # setup logger
    stationPair = '%s.%s'%(ref_net_sta, temp_net_sta)
    fn = os.path.join(outputPath, '%s.log'%(stationPair if not tracking_tag else '.'.join([stationPair, tracking_tag])))
    logger = setup_logger('%s.%s'%(ref_net_sta, temp_net_sta), fn)

    #######################################
    # Initialize variables for main loop
    #######################################
    startTime = UTCDateTime(start_time)
    endTime = UTCDateTime(end_time)

    cTime = startTime

    xcorrResultsDict = defaultdict(list)  # Results dictionary indexed by station-pair string
    windowCountResultsDict = defaultdict(list)  # Window-count dictionary indexed by station-pair string
    intervalStartTimesDict = defaultdict(list)
    intervalEndTimesDict = defaultdict(list)
    while (cTime < endTime):
        cStep = buffer_seconds

        if (cTime + cStep > endTime):
            cStep = endTime - cTime

        logger.info('====Time range  [%s - %s]====' % (str(cTime), str(cTime + cStep)))
        logger.info('Fetching data for station %s..' % (ref_net_sta))

        refSt = None
        try:
            rnc, rsc = ref_net_sta.split('.')
            stations = refds.get_stations(cTime, cTime + cStep, network=rnc, station=rsc)
            for codes in stations:
                if(ref_cha != codes[3]): continue
                refSt = refds.get_waveforms(codes[0], codes[1], codes[2], codes[3], cTime,
                                            cTime + cStep, automerge=False, trace_count_threshold=200)

                if (len(refSt) == 0): continue
                dropBogusTraces(refSt)

                if (verbose > 2):
                    logger.debug('\t\tData Gaps:')
                    refSt.print_gaps() # output sent to stdout; fix this
                    print ("\n")
                # end if
                # Merge reference station data. Note that we don't want to fill gaps; the
                # default merge() operation creates masked numpy arrays, which we can use
                # to detect and ignore windows that have gaps in their data.
                try:
                    refSt.merge()
                except:
                    logger.warning('\tFailed to merge traces..')
                    refSt = None
                # end try
            # end for
        except Exception as e:
            print(e)
            logger.warning('\tError encountered while fetching data. Skipping along..')

        if (refSt is None):
            logger.info('Failed to fetch data..')
            cTime += cStep
            continue
        elif (len(refSt) == 0):
            logger.info('Data source exhausted. Skipping time interval [%s - %s]' % (str(cTime), str(cTime + cStep)))
            cTime += cStep
            continue
        else:
            pass
            # print refSt

        logger.info('\tFetching data for station %s..' % (temp_net_sta))

        tempSt = None
        try:
            tnc, tsc = temp_net_sta.split('.')
            stations = tempds.get_stations(cTime, cTime + cStep, network=tnc, station=tsc)

            for codes in stations:
                if(temp_cha != codes[3]): continue
                tempSt = tempds.get_waveforms(codes[0], codes[1], codes[2], codes[3], cTime,
                                            cTime + cStep, automerge=False, trace_count_threshold=200)

                if (len(tempSt) == 0): continue
                dropBogusTraces(tempSt)

                if (verbose > 2):
                    logger.debug('\t\tData Gaps:')
                    tempSt.print_gaps() # output sent to stdout; fix this
                    print ("\n")
                # end if
                # Merge reference station data. Note that we don't want to fill gaps; the
                # default merge() operation creates masked numpy arrays, which we can use
                # to detect and ignore windows that have gaps in their data.
                try:
                    tempSt.merge()
                except:
                    logger.warning('\tFailed to merge traces..')
                    tempSt = None
                # end try
            # end if
        except Exception as e:
            print(e)
            logger.warning('\tError encountered while fetching data. Skipping along..')
        # end try

        if (tempSt is None):
            logger.info('Failed to fetch data..')
            cTime += cStep
            continue
        elif (len(tempSt) == 0):
            logger.info('Data source exhausted. Skipping time interval [%s - %s]' % (str(cTime), str(cTime + cStep)))
            cTime += cStep
            continue
        else:
            pass
            #print tempSt

        if (verbose > 2):
            logger.debug('\t\tData Gaps:')
            tempSt.print_gaps() # output sent to stdout; fix this
            print ("\n")

        logger.info('\tCross-correlating station-pair: %s' % (stationPair))
        xcl, winsPerInterval, \
        intervalStartSeconds, intervalEndSeconds = \
            xcorr2(refSt[0], tempSt[0], ref_sta_inv, temp_sta_inv,
                   window_seconds=window_seconds,
                   interval_seconds=interval_seconds,
                   resample_rate=resample_rate,
                   flo=flo, fhi=fhi,
                   clip_to_2std=clip_to_2std, whitening=whitening,
                   one_bit_normalize=one_bit_normalize,
                   envelope_normalize=envelope_normalize,
                   verbose=verbose, logger=logger)

        # Continue if no results were returned due to data-gaps
        if (xcl is None):
            logger.warning("\t\tWarning: no cross-correlation results returned for station-pair %s, " %
                  (stationPair) + " due to gaps in data.")
            cTime += cStep
            continue
        # end if

        xcorrResultsDict[stationPair].append(xcl)
        windowCountResultsDict[stationPair].append(winsPerInterval)

        intervalStartTimesDict[stationPair].append(cTime.timestamp + intervalStartSeconds)
        intervalEndTimesDict[stationPair].append(cTime.timestamp + intervalEndSeconds)

        cTime += cStep
    # wend (loop over time range)

    x = None
    skippedCount = 0
    # Concatenate results
    for k in list(xcorrResultsDict.keys()):
        combinedXcorrResults = None
        combinedWindowCountResults = None
        combinedIntervalStartTimes = None
        combinedIntervalEndTimes = None
        for i in np.arange(len(xcorrResultsDict[k])):
            if (i == 0):
                combinedXcorrResults = xcorrResultsDict[k][0]
                combinedWindowCountResults = windowCountResultsDict[k][0]
                combinedIntervalStartTimes = intervalStartTimesDict[k][0]
                combinedIntervalEndTimes = intervalEndTimesDict[k][0]

                # Generate time samples (only needs to be done once)
                if (x is None):
                    x = np.linspace(-window_seconds, window_seconds,
                                    xcorrResultsDict[k][0].shape[1])
                # end if

                if (ensemble_stack):
                    if (combinedXcorrResults.shape[0]>1):
                        combinedXcorrResults = np.expand_dims(np.sum(combinedXcorrResults,
                                                                     axis=0), axis=0)
                    # end if
                # end if
            else:
                if (combinedXcorrResults.shape[1] == xcorrResultsDict[k][i].shape[1]):
                    if(ensemble_stack):
                        if(xcorrResultsDict[k][i].shape[0]>1):
                            combinedXcorrResults += np.expand_dims(np.sum(xcorrResultsDict[k][i],
                                                                          axis=0), axis=0)
                        else:
                            combinedXcorrResults += xcorrResultsDict[k][i]
                        # end if
                    else:
                        combinedXcorrResults = np.concatenate((combinedXcorrResults,
                                                               xcorrResultsDict[k][i]))
                    # end if
                else:
                    if(ensemble_stack):
                        pass
                    else:
                        combinedXcorrResults = np.concatenate((combinedXcorrResults,
                                                               np.zeros((xcorrResultsDict[k][i].shape[0],
                                                                         combinedXcorrResults.shape[1]))))
                    # end if
                    logger.warning("\t\tVariable sample rates detected. Current station-pair: %s"%(k))
                # end if
                combinedWindowCountResults = np.concatenate((combinedWindowCountResults,
                                                             windowCountResultsDict[k][i]))
                combinedIntervalStartTimes = np.concatenate((combinedIntervalStartTimes,
                                                             intervalStartTimesDict[k][i]))
                combinedIntervalEndTimes = np.concatenate((combinedIntervalEndTimes,
                                                           intervalEndTimesDict[k][i]))
            # end if
        # end for

        # Replace lists with combined results
        xcorrResultsDict[k] = combinedXcorrResults
        windowCountResultsDict[k] = combinedWindowCountResults
        intervalStartTimesDict[k] = combinedIntervalStartTimes
        intervalEndTimesDict[k] = combinedIntervalEndTimes
    # end for

    # Save Results
    for i, k in enumerate(list(xcorrResultsDict.keys())):
        fn = os.path.join(outputPath, '%s.nc'%(k if not tracking_tag else '.'.join([k, tracking_tag])))

        root_grp = Dataset(fn, 'w', format='NETCDF4')
        root_grp.description = 'Cross-correlation results for station-pair: %s' % (k)

        # Dimensions
        root_grp.createDimension('lag', xcorrResultsDict[k].shape[1])
        lag = root_grp.createVariable('lag', 'f4', ('lag',))

        if(ensemble_stack):
            nsw = root_grp.createVariable('NumStackedWindows', 'i8')
            ist = root_grp.createVariable('IntervalStartTime', 'i8')
            iet = root_grp.createVariable('IntervalEndTime', 'i8')
            xc  = root_grp.createVariable('xcorr', 'f4', ('lag',))

            totalIntervalCount = int(np.sum(windowCountResultsDict[k] > 0))
            totalWindowCount = int(np.sum(windowCountResultsDict[k]))
            nsw[:] = totalWindowCount
            ist[:] = int(np.min(intervalStartTimesDict[k]))
            iet[:] = int(np.min(intervalEndTimesDict[k]))
            xc[:] = xcorrResultsDict[k] / float(totalIntervalCount)
        else:
            root_grp.createDimension('interval', xcorrResultsDict[k].shape[0])
            # Variables
            interval = root_grp.createVariable('interval', 'f4', ('interval',))
            nsw = root_grp.createVariable('NumStackedWindows', 'f4', ('interval',))
            ist = root_grp.createVariable('IntervalStartTimes', 'i8', ('interval',))
            iet = root_grp.createVariable('IntervalEndTimes', 'i8', ('interval',))
            xc = root_grp.createVariable('xcorr', 'f4', ('interval', 'lag',))

            # Populate variables
            interval[:] = np.arange(xcorrResultsDict[k].shape[0])
            nsw[:] = windowCountResultsDict[k]
            ist[:] = intervalStartTimesDict[k]
            iet[:] = intervalEndTimesDict[k]
            xc[:, :] = xcorrResultsDict[k]
        # end if

        lag[:] = x

        root_grp.close()
    # end for

    return x, xcorrResultsDict, windowCountResultsDict
# end func

