import datetime
import glob, os
from os.path import join, exists
import obspy
from obspy.core import Stream, UTCDateTime
from obspy import read, Trace
from obspy.signal.cross_correlation import xcorr
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy
import math
from fft import *
from collections import defaultdict
from netCDF4 import Dataset


# Greatest common divisor of more than 2 numbers.
def gcd(*numbers):
    """Return the greatest common divisor of the given integers"""
    from fractions import gcd
    return reduce(gcd, numbers)


# Assuming numbers are positive integers.

def lcm(*numbers):
    """Return lowest common multiple."""

    def lcm(a, b):
        return (a * b) // gcd(a, b)

    return reduce(lcm, numbers, 1)


def butterworthBPF(spectrum, samprate, f_lo, f_hi):
    f = np.fft.fftfreq(spectrum.shape[0], 1.0 / samprate)
    p = 4  # poles
    gain_lpf = 1.0 / (1 + (f / f_hi) ** (2 * p))
    gain_hpf = 1.0 / (1 + (f / f_lo) ** (2 * p))
    gain_bpf = gain_lpf - gain_hpf
    spectrum *= gain_bpf
    return spectrum


def zeropad(tr, padlen):
    assert (tr.shape[0] < padlen)
    padded = np.zeros(padlen)
    padded[0:tr.shape[0]] = tr
    return padded


def zeropad_ba(tr, padlen):
    assert (tr.shape[0] < padlen)
    padded = np.zeros(padlen, dtype=np.complex_)
    s = int((padlen - tr.shape[0]) / 2)
    padded[s:(s + tr.shape[0])] = scipy.fftpack.fftshift(tr)
    return scipy.fftpack.ifftshift(padded)


def zeropad_ba_old(tr, padlen):
    assert (tr.shape[0] < padlen)
    padded = np.zeros(padlen, dtype=np.complex_)
    s = int((padlen - tr.shape[0]) / 2)
    padded[s:(s + tr.shape[0])] = tr
    return padded


def taperDetrend(tr, taperlen):
    """ First detrend, then tapers signal"""
    tr -= np.median(tr)
    tr[0:taperlen] *= 0.5 * (1 + np.cos(np.linspace(-math.pi, 0, taperlen)))
    tr[-taperlen:] *= 0.5 * (1 + np.cos(np.linspace(0, math.pi, taperlen)))
    return tr


def xcorr2(tr1, tr2, window_seconds=3600, interval_seconds=86400,
           flo=0.9, fhi=1.1, verbose=1):
    sr1 = tr1.stats.sampling_rate
    sr2 = tr2.stats.sampling_rate
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
    resll = []

    intervalCount = 0
    windowsPerInterval = []  # Stores the number of windows processed per interval

    while itr1s < lentr1_all and itr2s < lentr2_all:
        itr1e = min(lentr1_all, itr1s + interval_samples_1)
        itr2e = min(lentr2_all, itr2s + interval_samples_2)
        sr = max(sr1, sr2)
        xcorlen = int(2 * window_seconds * sr - 1)
        fftlen = 2 ** (int(np.log2(xcorlen)) + 1)

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

                tr1_d = taperDetrend(tr1_d, int(sr1 / flo))
                tr2_d = taperDetrend(tr2_d, int(sr2 / flo))
                if (sr1 < sr2):
                    fftlen2 = fftlen
                    fftlen1 = int((fftlen2 * 1.0 * sr1) / sr)
                    outdims2 = np.array([fftlen2])
                    outdims1 = np.array([fftlen1])
                    rf = zeropad_ba(fftn(zeropad(tr1_d, fftlen1), shape=outdims1), fftlen2) * fftn(
                        zeropad(ndflip(tr2_d), fftlen2), shape=outdims2)
                elif (sr1 > sr2):
                    fftlen1 = fftlen
                    fftlen2 = int((fftlen1 * 1.0 * sr2) / sr)
                    outdims2 = np.array([fftlen2])
                    outdims1 = np.array([fftlen1])
                    rf = fftn(zeropad(tr1_d, fftlen1), shape=outdims1) * zeropad_ba(
                        fftn(zeropad(ndflip(tr2_d), fftlen2), shape=outdims2), fftlen1)
                else:
                    fftlen = 2 ** (int(np.log2(2 * window_samples_1 - 1)) + 1)
                    outdims = np.array([fftlen])
                    rf = fftn(zeropad(tr1_d, fftlen), shape=outdims) * fftn(zeropad(ndflip(tr2_d), fftlen),
                                                                            shape=outdims)
                # end if

                resl.append(rf)
                windowCount += 1
            # end if
            wtr1s = wtr1e
            wtr2s = wtr2e
        # end while (windows within interval)

        if (verbose > 1):
            print '\t\t\tProcessed %d windows in interval %d' % (windowCount, intervalCount)
        # end fi

        itr1s = itr1e
        itr2s = itr2e
        intervalCount += 1

        # Append an array of zeros if no windows were processed for the current interval
        if (windowCount == 0):
            resl.append(np.zeros(fftlen))
            if (verbose == 1):
                printf('\t\t\tWarning: No windows processed due to gaps in data in current interval')
            # end if
        # end if

        windowsPerInterval.append(windowCount)

        step = np.sign(np.fft.fftfreq(fftlen, 1.0 / sr))
        mean = reduce((lambda tx, ty: tx + ty), resl) / len(resl)
        mean = butterworthBPF(mean, sr, flo, fhi)
        mean = mean + step * mean  # compute analytic
        mean = ifftn(mean)

        # Compute magnitude of mean
        mean = np.abs(mean)

        # mask out the zero line for the max
        zerolower = int(xcorlen / 2 - sr / flo)
        zeroupper = int(xcorlen / 2 + sr / flo)
        maskmean = mean.view(np.ma.MaskedArray)
        maskmean[zerolower:zeroupper] = np.ma.masked

        #normFactor = np.ma.max(maskmean)
        #normFactor = 3.0 * np.ma.std(maskmean)  # three sigma
        normFactor = np.max(mean)

        # mean can be 0 for a null result
        if(normFactor > 0):
            mean /= normFactor
        #end if

        resll.append(mean[:xcorlen])
    # end while (iteration over intervals)

    if (len(resll)):
        return np.array(resll), np.array(windowsPerInterval)
    else:
        return None, None
    # end if
# end func

def IntervalStackXCorr(refds, tempds, tempds_db,
                       start_time, end_time,
                       ref_station_ids=None, temp_station_ids=None,
                       channel_wildcard='*Z',
                       refst_dec_factor=None, tempst_dec_factor=None,
                       buffer_seconds=864000, interval_seconds=86400,
                       window_seconds=3600, flo=0.9, fhi=1.1,
                       outputPath='/tmp', verbose=1):
    """
    :type refds: ASDFDataSet
    :param refds: ASDFDataSet containing reference-station data
    :type tempds: ASDFDataSet
    :param tempds: ASDFDataSet containing temporary-stations data
    :type tempds_db: SeisDB
    :param tempds_db: Json database corresponding to the temporary-stations ASDFDataSet.
    :type start_time: str
    :param: start_time: Start-time (UTCDateTime format) for data to be used in cross-correlation
    :type end_time: str
    :param: end_time: End-time (UTCDateTime format) for data to be used in cross-correlation
    :type ref_station_ids: List of strings
    :param ref_station_ids: List of reference station names to be used from the ASDFDataSet. Use
                            ['*'] for all stations available.
    :type temp_station_ids: List of strings
    :param temp_station_ids: List of temporary station names to be used from the ASDFDataSet.  Use
                             ['*'] for all stations available.
    :type channel_wildcard: str
    :param channel_wildcard: Wildcard to be used to filter channels. Default is '*Z'
    :type refst_dec_factor: int
    :param refst_dec_factor: Decimation factor for reference station data
    :type tempst_dec_factor: int
    :param tempst_dec_factor: Decimation factor for temporary stations data
    :type buffer_seconds: int
    :param buffer_seconds: The amount of data to be fetched per call from the ASDFDataSets, becuase
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
    :type verbose: int
    :param verbose: Verbosity of printouts. Default is 1; maximum is 3.
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

    if (ref_station_ids == ['*']):
        # Fetch station names if a wildcard was passed
        ref_station_ids = [x.split('.')[1] for x in refds.waveforms.list()]
    if (temp_station_ids == ['*']):
        # Fetch station names if a wildcard was passed
        temp_station_ids = [x.split('.')[1] for x in tempds.waveforms.list()]

    startTime = UTCDateTime(start_time)
    endTime = UTCDateTime(end_time)

    cTime = startTime

    xcorrResultsDict = defaultdict(list)  # Results dictionary indexed by station-pair string
    windowCountResultsDict = defaultdict(list)  # Window-count dictionary indexed by station-pair string
    while (cTime < endTime):
        cStep = buffer_seconds

        if (cTime + cStep > endTime):
            cStep = endTime - cTime

        print('=== Time Range: [%s - %s] ===' % (str(cTime), str(cTime + cStep)))
        for refStId in ref_station_ids:
            print('\tFetching data for Ref. Station %s..' % (refStId))
            refSt = refds.get_waveforms("*", refStId, "*", channel_wildcard, cTime,
                                        cTime + cStep, '*')
            if (refSt is None):
                print('\t\tFailed to fetch Ref. Station data. Skipping time interval..')
                continue
            elif (len(refSt) == 0):
                print('\t\tRef. Station data source exhausted. Skipping time interval..')
                continue
            # end for

            # Decimate traces if a decimation factor for reference station data has been
            # specified.
            if (refst_dec_factor is not None):
                for tr in refSt:
                    tr.decimate(refst_dec_factor)

            if (verbose > 2):
                print('\t\tData Gaps:')
                refSt.print_gaps()
                print ("\n")
            # end if
            # Merge reference station data. Note that we don't want to fill gaps; the
            # default merge() operation creates masked numpy arrays, which we can use
            # to detect and ignore windows that have gaps in their data.
            refSt.merge()

            for tempStId in temp_station_ids:
                stationPair = refStId + '.' + tempStId

                print('\tFetching data for Tmp. Station %s..' % (tempStId))
                tempSt = None

                if(tempds_db is not None):
                    # Fetching data from a large ASDF file, as is the case with temporary
                    # stations data, is significantly faster using a Json database.
                    tempSt = tempds_db.fetchDataByTime(tempds, tempStId, channel_wildcard,
                                                       cTime.timestamp,
                                                       (cTime + cStep).timestamp,
                                                       decimation_factor=tempst_dec_factor)
                else:
                    # Otherwise fetch data using the pyasdf interface.
                    tempSt = tempds.get_waveforms("*", tempStId, "*", channel_wildcard, cTime,
                                                  cTime + cStep, '*')

                    # Decimate traces
                    for tr in tempSt:
                        tr.decimate(tempst_dec_factor)

                    # Merge traces, preserving data gaps (see above).
                    tempSt.merge()
                #end if

                if (tempSt is None):
                    print('\t\tFailed to fetch Tmp. Station data. Skipping station..')
                    continue
                elif (len(tempSt) == 0):
                    print('\t\tTmp. Station data source exhausted. Skipping time interval..')
                    continue
                # end for

                if (verbose > 2):
                    print('\t\t\tData Gaps:')
                    tempSt.print_gaps()
                    print ("\n")
                # end if

                print('\t\tCross-correlating station-pair: %s' % (stationPair))
                xcl, winsPerInterval = xcorr2(refSt[0], tempSt[0], window_seconds=window_seconds,
                                              interval_seconds=interval_seconds,
                                              flo=flo, fhi=fhi, verbose=verbose)

                # Continue if no results were returned due to data-gaps
                if (xcl is None):
                    print("\t\tWarning: no cross-correlation results returned for station-pair %s, " %
                          (stationPair) + " due to gaps in data.")
                    continue
                # end if

                xcorrResultsDict[stationPair].append(xcl)
                windowCountResultsDict[stationPair].append(winsPerInterval)
                # end for (loop over temporary stations)
        # end for (loop over reference stations)

        cTime += cStep
        print "\n"
    # wend (loop over time range)

    print('=== Collating Results ===')
    x = None
    # Concatenate results
    for k in xcorrResultsDict.keys():
        combinedXcorrResults = None
        combinedWindowCountResults = None

        for i in np.arange(len(xcorrResultsDict[k])):
            if (i == 0):
                combinedXcorrResults = xcorrResultsDict[k][0]
                combinedWindowCountResults = windowCountResultsDict[k][0]

                # Generate time samples (only needs to be done once)
                if (x is None):
                    x = np.linspace(-window_seconds, window_seconds,
                                    xcorrResultsDict[k][0].shape[1])
                    # end if
            else:
                combinedXcorrResults = np.concatenate((combinedXcorrResults,
                                                       xcorrResultsDict[k][i]))
                combinedWindowCountResults = np.concatenate((combinedWindowCountResults,
                                                             windowCountResultsDict[k][i]))
                # end if
        # end for

        # Replace lists with combined results
        xcorrResultsDict[k] = combinedXcorrResults
        windowCountResultsDict[k] = combinedWindowCountResults
    # end for

    # Save Results
    for i, k in enumerate(xcorrResultsDict.keys()):
        fn = os.path.join(outputPath, '%s.nc'%(k))
        print('\t Writing: %s'%(fn))

        root_grp = Dataset(fn, 'w', format='NETCDF4')
        root_grp.description = 'Cross-correlation results for station-pair: %s' % (k)

        # Dimensions
        root_grp.createDimension('time', xcorrResultsDict[k].shape[0])
        root_grp.createDimension('lag', xcorrResultsDict[k].shape[1])

        # Variables
        time = root_grp.createVariable('time', 'f4', ('time',))
        lag = root_grp.createVariable('lag', 'f4', ('lag',))
        nsw = root_grp.createVariable('NumStackedWindows', 'f4', ('time',))
        xc = root_grp.createVariable('xcorr', 'f4', ('time', 'lag',))

        # Populate variables
        time[:] = np.arange(xcorrResultsDict[k].shape[0])
        lag[:] = x
        nsw[:] = windowCountResultsDict[k]
        xc[:, :] = xcorrResultsDict[k][::-1, :]  # Flipping rows
        root_grp.close()
    # end for

    return x, xcorrResultsDict, windowCountResultsDict
# end func

'''
# Leaving the earlier implementation and associated tests as is for the time being.
# This should be cleaned up and tests implemented in another ticket. 

def IntervalStackXCorr(refstn, st,
                       window_seconds=3600,
                       interval_seconds=86400, flo=0.9, fhi=1.1):
    comp_list = []
    a_list = []
    b_list = []
    ref_list = []
    stn_list = []
    xcorr_list = []
    xcorr_x_list = []
    for tr_1 in refstn:
        station_1 = tr_1.stats.station

        for tr_2 in st:
            station_2 = tr_2.stats.station

            comp_list.append(station_2 + '_' + station_1)

            xcorr_func = xcorr2(tr_1, tr_2, window_seconds, interval_seconds, flo, fhi)

            tr_1_len = tr_1.stats.endtime.timestamp - tr_1.stats.starttime.timestamp
            tr_2_len = tr_2.stats.endtime.timestamp - tr_2.stats.starttime.timestamp
            xcorr_len = tr_2_len - tr_1_len
            xcorr_start = tr_2.stats.starttime.timestamp

            X = np.linspace(-window_seconds, window_seconds, xcorr_func.shape[1])
            ref_list.append(tr_1)
            stn_list.append(tr_2)
            xcorr_list.append(xcorr_func)
            xcorr_x_list.append(X)
    return [xcorr_list, xcorr_x_list, comp_list]


def saveXCorr(xcorr_list, xcorr_x_list, xcor_output_dir, figname):
    os.chdir(xcor_output_dir)
    for i, l in enumerate(xcorr_list):
        np.savetxt(figname + '.xcor', xcorr_list[i])


def saveXCorrPlot(xcorr_list, xcorr_x_list, plot_output_dir, figname, comp_list):
    nplots = len(comp_list)
    if nplots == 0:
        return
    fig, ax = plt.subplots(1, nplots, squeeze=False)
    fig.set_figheight(23.39)  # 11.69)
    fig.set_figwidth(16.53)  # 8.27)

    for i, comp_xcorr in enumerate(comp_list):
        ax[0, i].pcolormesh(xcorr_x_list[i], np.arange(xcorr_list[i].shape[0] + 1), xcorr_list[i], vmin=0, vmax=1,
                            cmap='jet')  # , label=comp_xcorr, lw=1.2)
        ax[0, i].set_ylabel('Day')
        ax[0, i].set_xlabel('Lag (s)')
        ax[0, i].set_title(comp_xcorr)
        ax[0, i].invert_yaxis()

    os.chdir(plot_output_dir)
    plt.savefig(figname + '.png', dpi=300)
    plt.clf()

if __name__ == "__main__":
    sdr = '/g/data/ha3/Passive/Ref/'
    fn1 = 'AUQLP__BHZ__-14-05-23-00-00-00.msd'
    fn1b = 'AUQLP__BHZ__-14-05-24-00-00-00.msd'
    fn2 = 'AUINKA_BHZ__-14-05-23-00-00-00.msd'
    fn2b = 'AUINKA_BHZ__-14-05-24-00-00-00.msd'
    st1 = read(sdr + fn1)
    st1 += read(sdr + fn1b)
    st1.merge()
    st2 = read(sdr + fn2)
    st2 += read(sdr + fn2b)
    st2.merge()
    st2.resample(80.0, no_filter=True)
    print(st1)
    print(st2)
    ylist, x, comp_list = IntervalStackXCorr(st1, st2)
    print ylist
    print x
    print comp_list
    st1.resample(80.0, no_filter=True)
    ylist2, x2, comp_list2 = IntervalStackXCorr(st1, st2)
    print ylist2
    print x2
    print comp_list2
    # confirm the error is zero
    # stdyratio = np.std(ylist2[0]/ylist[0])
    maxdiff = np.max(np.abs(ylist2[0] - ylist[0]))
    saveXCorrPlot(ylist, x, '/g/data/ha3/', 'outxcor_dsr4', comp_list)
    saveXCorr(ylist, x2, '/g/data/ha3/', 'outxcor_8040')
    saveXCorr(ylist2, x2, '/g/data/ha3/', 'outxcor_8080')
    # print "Std y ratio = " + str(stdyratio)
    print "Max diff = " + str(maxdiff)
    # plt.subplot(211)
    # plt.plot(ylist[0][179999:180001])
    # plt.subplot(212)
    # plt.plot(ylist2[0][179999:180001])
    # plt.show()
    assert maxdiff < 0.05, "The xcor after resampling is not the same"
'''