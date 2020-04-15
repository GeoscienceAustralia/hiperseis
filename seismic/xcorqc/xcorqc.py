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
import logging
import math
from collections import defaultdict

import numpy as np
import scipy

from obspy.core import Stream, UTCDateTime, Stats
from obspy import read, Trace
from obspy.signal.cross_correlation import xcorr
from obspy.signal.detrend import simple, spline
from obspy.signal.filter import bandpass, highpass, lowpass
from obspy.geodetics.base import gps2dist_azimuth
from scipy import signal

from seismic.xcorqc.fft import *
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.xcorqc.utils import drop_bogus_traces, get_stream
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

    logger = logging.getLogger(name+log_file)
    logger.setLevel(level)
    logger.addHandler(handler)
    logger.propagate = False
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


def whiten(a, sampling_rate, window_freq=0):
    """
    Applies spectral whitening to trace samples. When window_freq=0, all frequency bins are normalized
    by their amplitudes, i.e. all frequency bins end up with an amplitude of 1. When window_freq is
    nonzero, a smoothed amplitude spectrum (smoothing window length is as computed below) is used
    to normalize the frequency bins.

    :param a: trace samples
    :param sampling_rate: sampling rate
    :param window_freq: smoothing window length (Hz)
    return: spectrally whitened samples
    """
    # frequency step
    npts = a.shape[0]
    deltaf = sampling_rate / npts

    ffta = np.fft.rfft(a)

    # smooth amplitude spectrum
    halfwindow = int(round(window_freq / deltaf / 2.0))

    if halfwindow > 0:
        # moving average
        weight = np.convolve(np.abs(ffta), np.ones(halfwindow * 2 + 1) / (halfwindow * 2 + 1), mode='same')
    else:
        weight = np.abs(ffta)

    ss = ffta / weight

    a = np.fft.irfft(ss)

    return a
# end func


def xcorr2(tr1, tr2, sta1_inv=None, sta2_inv=None,
           instrument_response_output='vel', water_level=50.,
           window_seconds=3600, window_overlap=0.1, window_buffer_length=0,
           interval_seconds=86400, taper_length=0.05, resample_rate=None,
           flo=None, fhi=None, clip_to_2std=False, whitening=False,
           whitening_window_frequency=0, one_bit_normalize=False, envelope_normalize=False,
           verbose=1, logger=None):

    # Length of window_buffer in seconds
    window_buffer_seconds = window_buffer_length * window_seconds
    adjusted_taper_length = taper_length
    if window_buffer_seconds:
        # adjust taper length
        adjusted_taper_length = taper_length / (1. + window_buffer_length * 2.)
    # end if

    sr1 = tr1.stats.sampling_rate
    sr2 = tr2.stats.sampling_rate
    sr1_orig = sr1
    sr2_orig = sr2
    tr1_d_all = tr1.data  # refstn
    tr2_d_all = tr2.data
    lentr1_all = tr1_d_all.shape[0]
    lentr2_all = tr2_d_all.shape[0]
    window_samples_1 = (window_seconds + 2*window_buffer_seconds) * sr1
    window_samples_2 = (window_seconds + 2*window_buffer_seconds) * sr2
    interval_samples_1 = interval_seconds * sr1
    interval_samples_2 = interval_seconds * sr2
    # sr = 0
    resll = []

    # set day-aligned start-indices
    maxStartTime = max(tr1.stats.starttime, tr2.stats.starttime)
    dayAlignedStartTime = UTCDateTime(year=maxStartTime.year, month=maxStartTime.month,
                                      day=maxStartTime.day)
    itr1s = (dayAlignedStartTime - tr1.stats.starttime) * sr1
    itr2s = (dayAlignedStartTime - tr2.stats.starttime) * sr2

    if resample_rate:
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

        while (itr1s < 0) or (itr2s < 0):
            itr1s += (window_samples_1 - 2*window_buffer_seconds*sr1_orig) - \
                     (window_samples_1 - 2*window_buffer_seconds*sr1_orig) * window_overlap
            itr2s += (window_samples_2 - 2*window_buffer_seconds*sr2_orig) - \
                     (window_samples_2 - 2*window_buffer_seconds*sr2_orig) * window_overlap
        # end while

        if (np.fabs(itr1e - itr1s) < sr1_orig) or (np.fabs(itr2e - itr2s) < sr2_orig):
            itr1s = itr1e
            itr2s = itr2e
            continue
        # end if

        if (tr1.stats.starttime + itr1s / sr1_orig != tr2.stats.starttime + itr2s / sr2_orig):
            if logger:
                logger.warning('Detected misaligned traces..')

        windowCount = 0
        wtr1s = int(itr1s)
        wtr2s = int(itr2s)
        resl = []

        while wtr1s < itr1e and wtr2s < itr2e:
            wtr1e = int(min(itr1e, wtr1s + window_samples_1))
            wtr2e = int(min(itr2e, wtr2s + window_samples_2))

            # Discard small windows
            if ((wtr1e - wtr1s < window_samples_1) or (wtr2e - wtr2s < window_samples_2) or
                (wtr1e - wtr1s < sr1_orig) or (wtr2e - wtr2s < sr2_orig)):
                wtr1s = int(np.ceil(itr1e))
                wtr2s = int(np.ceil(itr2e))
                continue
            # end if

            # Discard windows with masked regions, i.e. with gaps or windows that are all zeros
            if (not (np.ma.is_masked(tr1_d_all[wtr1s:wtr1e]) or
                     np.ma.is_masked(tr2_d_all[wtr2s:wtr2e]) or
                     np.sum(tr1_d_all[wtr1s:wtr1e]) == 0 or
                     np.sum(tr2_d_all[wtr2s:wtr2e]) == 0)):

                # logger.info('%s, %s' % (tr1.stats.starttime + wtr1s / 200., tr1.stats.starttime + wtr1e / sr1_orig))
                # logger.info('%s, %s' % (tr2.stats.starttime + wtr2s / 200., tr2.stats.starttime + wtr2e / sr2_orig))

                tr1_d = np.array(tr1_d_all[wtr1s:wtr1e], dtype=np.float32)
                tr2_d = np.array(tr2_d_all[wtr2s:wtr2e], dtype=np.float32)

                # STEP 1: detrend
                tr1_d = signal.detrend(tr1_d)
                tr2_d = signal.detrend(tr2_d)

                # STEP 2: demean
                tr1_d -= np.mean(tr1_d)
                tr2_d -= np.mean(tr2_d)

                # STEP 3: remove response
                if sta1_inv:
                    resp_tr1 = Trace(data=tr1_d,
                                     header=Stats(header={'sampling_rate': sr1_orig,
                                                          'npts': len(tr1_d),
                                                          'network': tr1.stats.network,
                                                          'station': tr1.stats.station,
                                                          'location': tr1.stats.location,
                                                          'channel': tr1.stats.channel,
                                                          'starttime': tr1.stats.starttime + float(wtr1s)/sr1_orig,
                                                          'endtime': tr1.stats.starttime + float(wtr1e)/sr1_orig}))
                    try:
                        resp_tr1.remove_response(inventory=sta1_inv, output=instrument_response_output.upper(),
                                                 water_level=water_level)
                    except Exception as e:
                        logger.error(str(e))
                    # end try

                    tr1_d = resp_tr1.data
                # end if

                # remove response
                if sta2_inv:
                    resp_tr2 = Trace(data=tr2_d,
                                     header=Stats(header={'sampling_rate': sr2_orig,
                                                          'npts': len(tr2_d),
                                                          'network': tr2.stats.network,
                                                          'station': tr2.stats.station,
                                                          'location': tr2.stats.location,
                                                          'channel': tr2.stats.channel,
                                                          'starttime': tr2.stats.starttime + float(wtr2s)/sr2_orig,
                                                          'endtime': tr2.stats.starttime + float(wtr2e)/sr2_orig}))
                    try:
                        resp_tr2.remove_response(inventory=sta2_inv, output=instrument_response_output.upper(),
                                                 water_level=water_level)
                    except Exception as e:
                        logger.error(str(e))
                    # end try

                    tr2_d = resp_tr2.data
                # end if

                # STEPS 4, 5: resample after lowpass @ resample_rate/2 Hz
                if resample_rate:
                    tr1_d = lowpass(tr1_d, resample_rate/2., sr1_orig, corners=2, zerophase=True)
                    tr2_d = lowpass(tr2_d, resample_rate/2., sr2_orig, corners=2, zerophase=True)

                    tr1_d = Trace(data=tr1_d,
                                  header=Stats(header={'sampling_rate': sr1_orig,
                                                       'npts': window_samples_1})).resample(resample_rate,
                                                                                           no_filter=True).data
                    tr2_d = Trace(data=tr2_d,
                                  header=Stats(header={'sampling_rate': sr2_orig,
                                                       'npts': window_samples_2})).resample(resample_rate,
                                                                                           no_filter=True).data
                # end if

                # STEP 6: Bandpass
                if flo and fhi:
                    tr1_d = bandpass(tr1_d, flo, fhi, sr1, corners=2, zerophase=True)
                    tr2_d = bandpass(tr2_d, flo, fhi, sr2, corners=2, zerophase=True)
                # end if

                # STEP 7: time-domain normalization
                # clip to +/- 2*std
                if clip_to_2std:
                    std_tr1 = np.std(tr1_d)
                    std_tr2 = np.std(tr2_d)
                    clip_indices_tr1 = np.fabs(tr1_d) > 2 * std_tr1
                    clip_indices_tr2 = np.fabs(tr2_d) > 2 * std_tr2

                    tr1_d[clip_indices_tr1] = 2 * std_tr1 * np.sign(tr1_d[clip_indices_tr1])
                    tr2_d[clip_indices_tr2] = 2 * std_tr2 * np.sign(tr2_d[clip_indices_tr2])
                # end if

                # 1-bit normalization
                if one_bit_normalize:
                    tr1_d = np.sign(tr1_d)
                    tr2_d = np.sign(tr2_d)
                # end if

                # Apply Rhys Hawkins-style default time domain normalization
                if (clip_to_2std == 0) and (one_bit_normalize == 0):
                    # 0-mean
                    tr1_d -= np.mean(tr1_d)
                    tr2_d -= np.mean(tr2_d)

                    # unit-std
                    tr1_d /= np.std(tr1_d)
                    tr2_d /= np.std(tr2_d)
                # end if

                # STEP 8: taper
                if adjusted_taper_length > 0:
                    tr1_d = taper(tr1_d, int(np.round(adjusted_taper_length*tr1_d.shape[0])))
                    tr2_d = taper(tr2_d, int(np.round(adjusted_taper_length*tr2_d.shape[0])))
                # end if

                # STEP 9: spectral whitening
                if whitening:
                    tr1_d = whiten(tr1_d, sr1, window_freq=whitening_window_frequency)
                    tr2_d = whiten(tr2_d, sr2, window_freq=whitening_window_frequency)

                    # STEP 10: taper
                    if adjusted_taper_length > 0:
                        tr1_d = taper(tr1_d, int(np.round(adjusted_taper_length*tr1_d.shape[0])))
                        tr2_d = taper(tr2_d, int(np.round(adjusted_taper_length*tr2_d.shape[0])))
                    # end if
                # end if

                # STEP 11: Final bandpass
                # apply zero-phase bandpass
                if flo and fhi:
                    tr1_d = bandpass(tr1_d, flo, fhi, sr1, corners=2, zerophase=True)
                    tr2_d = bandpass(tr2_d, flo, fhi, sr2, corners=2, zerophase=True)
                # end if

                if window_buffer_seconds:
                    # extract window of interest from buffered window
                    tr1_d = tr1_d[int(window_buffer_seconds*sr1):-int(window_buffer_seconds*sr1)]
                    tr2_d = tr2_d[int(window_buffer_seconds*sr2):-int(window_buffer_seconds*sr2)]
                # end if

                # cross-correlate waveforms
                if sr1 < sr2:
                    fftlen2 = fftlen
                    fftlen1 = int((fftlen2 * 1.0 * sr1) / sr)
                    rf = zeropad_ba(fftn(zeropad(tr1_d, fftlen1), shape=[fftlen1]), fftlen2) * fftn(
                        zeropad(ndflip(tr2_d), fftlen2), shape=[fftlen2])
                elif sr1 > sr2:
                    fftlen1 = fftlen
                    fftlen2 = int((fftlen1 * 1.0 * sr2) / sr)
                    rf = fftn(zeropad(tr1_d, fftlen1), shape=[fftlen1]) * zeropad_ba(
                        fftn(zeropad(ndflip(tr2_d), fftlen2), shape=[fftlen2]), fftlen1)
                else:
                    rf = fftn(zeropad(tr1_d, fftlen), shape=[fftlen]) * fftn(zeropad(ndflip(tr2_d), fftlen),
                                                                             shape=[fftlen])
                # end if

                if not np.isnan(rf).any():
                    resl.append(rf)
                    windowCount += 1
                # end if
            # end if

            wtr1s += int((window_samples_1 - 2*window_buffer_seconds*sr1_orig) -
                         (window_samples_1 - 2*window_buffer_seconds*sr1_orig) * window_overlap)
            wtr2s += int((window_samples_2 - 2*window_buffer_seconds*sr2_orig) -
                         (window_samples_2 - 2*window_buffer_seconds*sr2_orig) * window_overlap)
        # end while (windows within interval)

        if verbose > 1:
            if logger:
                logger.info('\tProcessed %d windows in interval %d' % (windowCount, intervalCount))
        # end if

        intervalStartSeconds.append(itr1s/sr1_orig + tr1.stats.starttime.timestamp)
        intervalEndSeconds.append(itr1e/sr1_orig + tr1.stats.starttime.timestamp)
        itr1s = itr1e
        itr2s = itr2e
        intervalCount += 1

        # Append an array of zeros if no windows were processed for the current interval
        if windowCount == 0:
            resl.append(np.zeros(fftlen))
            if verbose > 1:
                if logger:
                    logger.info('\tWarning: No windows processed due to gaps in data in current interval')
            # end if
        # end if

        windowsPerInterval.append(windowCount)

        if windowCount > 0:
            mean = reduce((lambda tx, ty: tx + ty), resl) / float(windowCount)
        else:
            mean = reduce((lambda tx, ty: tx + ty), resl)
        # end if

        if envelope_normalize:
            step = np.sign(np.fft.fftfreq(fftlen, 1.0 / sr))
            mean = mean + step * mean  # compute analytic
        # end if

        mean = ifftn(mean)

        if envelope_normalize:
            # Compute magnitude of mean
            mean = np.abs(mean)
            normFactor = np.max(mean)

            # mean can be 0 for a null result
            if normFactor > 0:
                mean /= normFactor
            # end if
        # end if

        resll.append(mean[:xcorlen])
    # end while (iteration over intervals)

    if len(resll):
        return np.array(resll), np.array(windowsPerInterval), \
               np.array(intervalStartSeconds, dtype='i8'), \
               np.array(intervalEndSeconds, dtype='i8'), \
               sr
    else:
        return None, None, None, None, sr
    # end if
# end func


def IntervalStackXCorr(refds, tempds,
                       start_time, end_time,
                       ref_net_sta, temp_net_sta,
                       ref_sta_inv, temp_sta_inv,
                       instrument_response_output,
                       water_level,
                       ref_cha,
                       temp_cha,
                       baz_ref_net_sta,
                       baz_temp_net_sta,
                       resample_rate=None,
                       taper_length=0.05,
                       buffer_seconds=864000, interval_seconds=86400,
                       window_seconds=3600, window_overlap=0.1, window_buffer_length=0,
                       flo=None, fhi=None,
                       clip_to_2std=False, whitening=False, whitening_window_frequency=0,
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
    :type instrument_response_output: str
    :param instrument_response_output: Output of instrument response correction; can be either 'vel' or 'disp'
    :type water_level: float
    :param water_level: Water-level used during instrument response correction
    :type ref_cha: str
    :param ref_cha: Channel name for the reference Dataset
    :type temp_cha: str
    :param temp_cha: Channel name for the temporary Dataset
    :type baz_ref_net_sta: float
    :param baz_ref_net_sta: Back-azimuth of ref station from temp station in degrees
    :type baz_temp_net_sta: float
    :param baz_temp_net_sta: Back-azimuth of temp station from ref station in degrees
    :type resample_rate: float
    :param resample_rate: Resampling rate (Hz). Applies to both data-sets
    :type taper_length: float
    :param taper_length: Taper length as a fraction of window length
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
    :type window_overlap: float
    :param window_overlap: Window overlap fraction. Default is 0.1.
    :type window_buffer_length: float
    :param window_buffer_length: Buffer length as a fraction of 'window-seconds' around actual data windows of
                                 interest. This helps exclude effects of tapering and other edge artefacts from
                                 data windows before cross-correlation. Default is 0
    :type flo: float
    :param flo: Lower frequency for Butterworth bandpass filter
    :type fhi: float
    :param fhi: Upper frequency for Butterworth bandpass filter
    :type clip_to_2std: bool
    :param clip_to_2std: Clip data in each window to +/- 2 standard deviations
    :type whitening: bool
    :param whitening: Apply spectral whitening
    :type whitening_window_frequency: float
    :param whitening_window_frequency: Window frequency (Hz) used to determine length of averaging window
                                       for smoothing spectral amplitude
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
    :return: 1: 1d np.array with time samples spanning [-window_samples+dt:window_samples-dt]
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
    if resample_rate and fhi:
        if resample_rate < 2*fhi:
            raise RuntimeError('Resample-rate should be >= 2*fmax')

    if clip_to_2std and one_bit_normalize:
        raise RuntimeError('Mutually exclusive parameterization: clip_to_2std and one-bit-normalizations'
                           'together is redundant')
    # end if

    # setup logger
    stationPair = '%s.%s' % (ref_net_sta, temp_net_sta)
    fn = os.path.join(outputPath, '%s.log' % (stationPair if not tracking_tag else
                                              '.'.join([stationPair, tracking_tag])))
    logger = setup_logger('%s.%s' % (ref_net_sta, temp_net_sta), fn)

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
    sr = 0
    while cTime < endTime:
        cStep = buffer_seconds

        if (cTime + cStep) > endTime:
            cStep = endTime - cTime

        logger.info('====Time range  [%s - %s]====' % (str(cTime), str(cTime + cStep)))
        logger.info('Fetching data for station %s..' % ref_net_sta)

        refSt = None
        try:
            rnc, rsc = ref_net_sta.split('.')
            refSt = get_stream(refds, rnc, rsc, ref_cha, cTime, cTime + cStep, baz=baz_ref_net_sta,
                               logger=logger, verbose=verbose)
        except Exception as e:
            logger.error('\t'+str(e))
            logger.warning('\tError encountered while fetching data. Skipping along..')

        if refSt is None:
            logger.info('Failed to fetch data..')
            cTime += cStep
            continue
        elif len(refSt) == 0:
            logger.info('Data source exhausted. Skipping time interval [%s - %s]' % (str(cTime), str(cTime + cStep)))
            cTime += cStep
            continue
        else:
            pass
            # print refSt
        # end if

        logger.info('\tFetching data for station %s..' % temp_net_sta)

        tempSt = None
        try:
            tnc, tsc = temp_net_sta.split('.')
            tempSt = get_stream(tempds, tnc, tsc, temp_cha, cTime, cTime + cStep, baz=baz_temp_net_sta,
                                logger=logger, verbose=verbose)
        except Exception as e:
            logger.error('\t'+str(e))
            logger.warning('\tError encountered while fetching data. Skipping along..')
        # end try

        if tempSt is None:
            logger.info('Failed to fetch data..')
            cTime += cStep
            continue
        elif len(tempSt) == 0:
            logger.info('Data source exhausted. Skipping time interval [%s - %s]' % (str(cTime), str(cTime + cStep)))
            cTime += cStep
            continue
        else:
            pass
            # print tempSt
        # end if

        if verbose > 2:
            logger.debug('\t\tData Gaps:')
            tempSt.print_gaps() # output sent to stdout; fix this
            print("\n")

        logger.info('\tCross-correlating station-pair: %s' % stationPair)
        xcl, winsPerInterval, \
        intervalStartSeconds, intervalEndSeconds, sr = \
            xcorr2(refSt[0], tempSt[0], ref_sta_inv, temp_sta_inv,
                   instrument_response_output=instrument_response_output,
                   water_level=water_level,
                   window_seconds=window_seconds,
                   window_overlap=window_overlap,
                   window_buffer_length=window_buffer_length,
                   interval_seconds=interval_seconds,
                   resample_rate=resample_rate,
                   taper_length=taper_length,
                   flo=flo, fhi=fhi,
                   clip_to_2std=clip_to_2std,
                   whitening=whitening,
                   whitening_window_frequency=whitening_window_frequency,
                   one_bit_normalize=one_bit_normalize,
                   envelope_normalize=envelope_normalize,
                   verbose=verbose, logger=logger)

        # Continue if no results were returned due to data-gaps
        if xcl is None:
            logger.warning("\t\tWarning: no cross-correlation results returned for station-pair %s, " %
                  stationPair + " due to gaps in data.")
            cTime += cStep
            continue
        # end if

        xcorrResultsDict[stationPair].append(xcl)
        windowCountResultsDict[stationPair].append(winsPerInterval)

        intervalStartTimesDict[stationPair].append(intervalStartSeconds)
        intervalEndTimesDict[stationPair].append(intervalEndSeconds)

        cTime += cStep
    # wend (loop over time range)

    x = None
    # skippedCount = 0
    # Concatenate results
    for k in list(xcorrResultsDict.keys()):
        combinedXcorrResults = None
        combinedWindowCountResults = None
        combinedIntervalStartTimes = None
        combinedIntervalEndTimes = None
        for i in np.arange(len(xcorrResultsDict[k])):
            if i == 0:
                combinedXcorrResults = xcorrResultsDict[k][0]
                combinedWindowCountResults = windowCountResultsDict[k][0]
                combinedIntervalStartTimes = intervalStartTimesDict[k][0]
                combinedIntervalEndTimes = intervalEndTimesDict[k][0]

                # Generate time samples (only needs to be done once)
                if x is None:
                    dt = 1./sr
                    x = np.linspace(-window_seconds + dt, window_seconds - dt,
                                    xcorrResultsDict[k][0].shape[1])
                # end if

                if ensemble_stack:
                    if combinedXcorrResults.shape[0] > 1:
                        combinedXcorrResults = np.expand_dims(np.sum(combinedXcorrResults,
                                                                     axis=0), axis=0)
                    # end if
                # end if
            else:
                if combinedXcorrResults.shape[1] == xcorrResultsDict[k][i].shape[1]:
                    if ensemble_stack:
                        if xcorrResultsDict[k][i].shape[0] > 1:
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
                    if ensemble_stack:
                        pass
                    else:
                        combinedXcorrResults = np.concatenate((combinedXcorrResults,
                                                               np.zeros((xcorrResultsDict[k][i].shape[0],
                                                                         combinedXcorrResults.shape[1]))))
                    # end if
                    logger.warning("\t\tVariable sample rates detected. Current station-pair: %s" % k)
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
        fn = os.path.join(outputPath, '%s.nc' % (k if not tracking_tag else '.'.join([k, tracking_tag])))

        root_grp = Dataset(fn, 'w', format='NETCDF4')
        root_grp.description = 'Cross-correlation results for station-pair: %s' % k

        # Dimensions
        root_grp.createDimension('lag', xcorrResultsDict[k].shape[1])
        root_grp.createDimension('nchar', 10)

        lag = root_grp.createVariable('lag', 'f4', ('lag',))

        # Add metadata
        lon1 = root_grp.createVariable('Lon1', 'f4')
        lat1 = root_grp.createVariable('Lat1', 'f4')
        lon2 = root_grp.createVariable('Lon2', 'f4')
        lat2 = root_grp.createVariable('Lat2', 'f4')
        distance = root_grp.createVariable('Distance', 'f4')

        ref_sta_coords = refds.unique_coordinates[ref_net_sta]
        temp_sta_coords = tempds.unique_coordinates[temp_net_sta]
        lon1[:] = ref_sta_coords[0] if len(ref_sta_coords) == 2 else -999
        lat1[:] = ref_sta_coords[1] if len(ref_sta_coords) == 2 else -999
        lon2[:] = temp_sta_coords[0] if len(temp_sta_coords) == 2 else -999
        lat2[:] = temp_sta_coords[1] if len(temp_sta_coords) == 2 else -999
        if np.min([v != -999 for v in [lon1[:], lat1[:], lon2[:], lat2[:]]]):
            distance[:], _, _ = gps2dist_azimuth(lat1[:], lon1[:], lat2[:], lon2[:])
        # end if

        # Add data
        if ensemble_stack:
            nsw = root_grp.createVariable('NumStackedWindows', 'i8')
            avgnsw = root_grp.createVariable('AvgNumStackedWindowsPerInterval', 'f4')
            ist = root_grp.createVariable('IntervalStartTime', 'i8')
            iet = root_grp.createVariable('IntervalEndTime', 'i8')
            xc = root_grp.createVariable('xcorr', 'f4', ('lag',))

            totalIntervalCount = int(np.sum(windowCountResultsDict[k] > 0))
            totalWindowCount = int(np.sum(windowCountResultsDict[k]))
            nsw[:] = totalWindowCount
            avgnsw[:] = np.mean(windowCountResultsDict[k][windowCountResultsDict[k]>0])
            ist[:] = int(np.min(intervalStartTimesDict[k]))
            iet[:] = int(np.max(intervalEndTimesDict[k]))
            if totalIntervalCount > 0:
                xc[:] = xcorrResultsDict[k] / float(totalIntervalCount)
            else:
                xc[:] = xcorrResultsDict[k]
            # end if
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

        # Add and populate a new group for parameters used
        pg = root_grp.createGroup('Parameters')

        params = {'corr_chans': '%s.%s' % (ref_cha, temp_cha),
                  'instr_corr_applied_1': 1 if ref_sta_inv else 0,
                  'instr_corr_applied_2': 1 if temp_sta_inv else 0,
                  'instr_corr_output': instrument_response_output,
                  'instr_corr_water_level_db': water_level,
                  'resample_rate': resample_rate if resample_rate else -999,
                  'taper_length': taper_length,
                  'buffer_seconds': buffer_seconds,
                  'interval_seconds': interval_seconds,
                  'window_seconds': window_seconds,
                  'window_overlap': window_overlap,
                  'window_buffer_length': window_buffer_length,
                  'bandpass_fmin': flo if flo else -999,
                  'bandpass_fmax': fhi if fhi else -999,
                  'clip_to_2std': int(clip_to_2std),
                  'one_bit_normalize': int(one_bit_normalize),
                  'zero_mean_1std_normalize': int(clip_to_2std is False and one_bit_normalize is False),
                  'spectral_whitening': int(whitening),
                  'envelope_normalize': int(envelope_normalize),
                  'ensemble_stack': int(ensemble_stack)}

        if whitening:
            params['whitening_window_frequency'] = whitening_window_frequency

        for _k, _v in params.items():
            setattr(pg, _k, _v)
        # end for

        root_grp.close()
    # end for

    return x, xcorrResultsDict, windowCountResultsDict
# end func
