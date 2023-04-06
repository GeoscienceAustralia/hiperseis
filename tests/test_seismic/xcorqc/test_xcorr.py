#!/bin/env python
"""
Description:
    Tests cross-correlation functionality
    .......
    .......
   
References:
 
CreationDate:   5/20/19
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     5/20/19   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from obspy.core import Trace, Stats, UTCDateTime
from seismic.xcorqc.xcorqc import taper, whiten, zeropad_ba, xcorr2
import pytest
import numpy as np
import itertools
import os
import tempfile

@pytest.fixture(params=[100, 500, 1000])
def trace_length(request):
    return request.param

@pytest.fixture(params=[1, 10, 40])
def sampling_rate(request):
    return request.param

@pytest.fixture(params=[1, 10, 30, 50])
def other_sampling_rate(request):
    return request.param

@pytest.fixture(params=[2, 2.25, 4, 4.5])
def pad_length_factor(request):
    return request.param

@pytest.fixture(params=[0.1, 0.2, 0.3,  0.5])
def taper_fraction(request):
    return request.param


def test_whiten(trace_length, sampling_rate):
    trc = np.zeros(trace_length)

    result = whiten(trc, float(sampling_rate))

    # Whitening a null trace produces nans
    assert np.isnan(np.nanmean(result))
# end func


def test_zeropad_ba(trace_length, pad_length_factor):
    trc = np.random.random(trace_length)

    fft = np.fft.fft(trc)
    fft_pad = zeropad_ba(fft, int(trace_length*pad_length_factor))

    # Check Parseval's theorem holds
    assert np.allclose(np.sum(trc ** 2), np.sum(np.abs(fft) ** 2) / float(trace_length))
    assert np.allclose(np.sum(trc ** 2), np.sum(np.abs(fft_pad) ** 2) / float(trace_length))
# end func


def test_taper(trace_length, taper_fraction):
    taper_length = int(trace_length * taper_fraction)
    trc = np.random.random(trace_length)

    trc_tapered = taper(trc.copy(), taper_length)

    trc_untapered = trc_tapered / taper(np.ones(trace_length), taper_length)

    # The untapered trace should equal the original trace except at indices 0 and N-1,
    # which are both nans.
    assert np.allclose(trc[1:int(trace_length)-1], trc_untapered[1:int(trace_length)-1])
# end func


def test_xcorr(trace_length, sampling_rate, other_sampling_rate):

    # Testing x-correlation of traces with different sampling rates

    tr1 = Trace(data=np.sin(np.linspace(0, 20*np.pi, trace_length*sampling_rate)),
                header=Stats(header={'sampling_rate': sampling_rate,
                                     'npts': trace_length*sampling_rate,
                                     'network': 'AU',
                                     'station': 'A'}))
    tr2 = Trace(data=np.cos(np.linspace(0, 200*np.pi, trace_length*other_sampling_rate)),
                header=Stats(header={'sampling_rate': other_sampling_rate,
                                     'npts': trace_length*other_sampling_rate,
                                     'network': 'AU',
                                     'station': 'B'}))

    arr, _, _, _, _, _, _ = xcorr2(tr1, tr2, window_seconds=trace_length/2, interval_seconds=trace_length)

    # Mean of the x-correlation function should be close to zero
    assert np.allclose(np.abs(np.mean(arr)), 0, atol=1e-2)
# end func


def generate_test_data(data_length_seconds:int,
                       starttime1:UTCDateTime,
                       starttime2:UTCDateTime):
    tr1 = Trace(data=np.random.random(data_length_seconds))
    tr2 = Trace(data=np.random.random(data_length_seconds))

    tr1.stats.starttime = starttime1
    tr2.stats.starttime = starttime2

    tr1.stats.network='AU'
    tr1.stats.station='X'
    tr1.stats.location=''
    tr1.stats.channel='BHZ'

    tr2.stats.network='AU'
    tr2.stats.station='Y'
    tr2.stats.location=''
    tr2.stats.channel='BHZ'

    return tr1, tr2
# end func

def test_window_counts():
    """
    Tests whether interval/window counts, with stacking disabled, are as expected
    """

    def expected_lines_iter():
        path = os.path.dirname(os.path.abspath(__file__))
        efn = os.path.join(path, 'data/expected/window_counts.txt')
        elines = open(efn, 'r').readlines()
        for eline in elines:
            yield eline
        # end for
    # end for
    elines = expected_lines_iter()

    DATA_LEN = [86400, 172800] # 1, 2 days
    WINDOW_SECONDS = [3600, 7200] # 1, 2 hrs
    READ_AHEAD_WINDOWS = [5, 10]
    OVERLAP = [0, 0.1, 0.75] # 0, 10, 75 %
    WINDOW_BUFFER_LENGTH = [0.1, 0.2] # 10, 20 %, at each end
    START_TIME_DELTAS = [[0, 0], [-200, 0], [0, -200],
                         [100, 200]]

    output_path = str(tempfile.mkdtemp(suffix='_test'))
    #output_path = '/tmp/'

    ofn = os.path.join(output_path, 'window_counts.txt')
    ofh = open(ofn, 'w+')
    all_output = []
    i = 1
    for dlen, raw, wsec, olap, wbl, delta in itertools.product(DATA_LEN,
                                                               READ_AHEAD_WINDOWS,
                                                               WINDOW_SECONDS,
                                                               OVERLAP,
                                                               WINDOW_BUFFER_LENGTH,
                                                               START_TIME_DELTAS):
        olines = []
        tr1, tr2 = generate_test_data(dlen,
                                      UTCDateTime(0)+delta[0],
                                      UTCDateTime(0)+delta[0])

        isec = wsec * (1 - olap) * raw + wsec * wbl * 2 + olap * wsec * 2

        xcorr, winPerInterval, \
        iStart, iEnd, wStart, \
        wEnd, sr = xcorr2(tr1,
                          tr2,
                          window_seconds=wsec,
                          interval_seconds=isec,
                          window_overlap=olap,
                          window_buffer_length=wbl,
                          apply_stacking=True)

        olines.append("============[test: {}]============\n".format(i))
        olines.append('Params: raw: {}, wsec: {}, olap: {}, wbl: {}\n'.format(raw, wsec, olap, wbl))
        olines.append('X-corr shape: {}\n'.format(xcorr.shape))
        olines.append('Windows per interval: {}\n'.format(winPerInterval))
        olines.append("Interval start- and end-times:\n".format(i))
        for s, e in zip(iStart, iEnd):
            olines.append('{} - {}\n'.format(UTCDateTime(s + tr1.stats.starttime.timestamp),
                                             UTCDateTime(e + tr1.stats.starttime.timestamp)))
        olines.append('\n')
        olines.append("Window start- and end-times:\n".format(i))
        for s, e in zip(wStart, wEnd):
            olines.append('{} - {}\n'.format(UTCDateTime(s + tr1.stats.starttime.timestamp),
                                             UTCDateTime(e + tr1.stats.starttime.timestamp)))
        olines.append('\n')

        for oline in olines: all_output.append(oline)

        if(1):
            for oline in olines:
                eline = next(elines)

                if(oline != eline):
                    assert 0, 'Output: {} does not match expected: {}'.format(oline, eline)
                # end if
            # end for
        # end if

        i += 1
    # end for

    for line in all_output: ofh.write(line)
    ofh.close()
# end func

def test_stacking_window_counts():
    """
    Tests whether interval/window counts are as expected
    """

    def expected_lines_iter():
        path = os.path.dirname(os.path.abspath(__file__))
        efn = os.path.join(path, 'data/expected/window_counts_stacked.txt')
        elines = open(efn, 'r').readlines()
        for eline in elines:
            yield eline
        # end for
    # end for
    elines = expected_lines_iter()

    DATA_LEN = [86400, 172800] # 1, 2 days
    INTERVAL_SECONDS = [43200, 86400] # 12, 24 hrs
    WINDOW_SECONDS = [3600, 7200] # 1, 2 hrs
    OVERLAP = [0, 0.1, 0.75] # 0, 10, 75 %
    WINDOW_BUFFER_LENGTH = [0.1, 0.2] # 10, 20 %, at each end
    START_TIME_DELTAS = [[0, 0], [-200, 0], [0, -200],
                         [100, 200]]

    output_path = str(tempfile.mkdtemp(suffix='_test'))
    #output_path = '/tmp/'

    ofn = os.path.join(output_path, 'window_counts_stacked.txt')
    ofh = open(ofn, 'w+')
    all_output = []
    i = 1
    for dlen, isec, wsec, olap, wbl, delta in itertools.product(DATA_LEN,
                                                                INTERVAL_SECONDS,
                                                                WINDOW_SECONDS,
                                                                OVERLAP,
                                                                WINDOW_BUFFER_LENGTH,
                                                                START_TIME_DELTAS):
        olines = []
        tr1, tr2 = generate_test_data(dlen,
                                      UTCDateTime(0)+delta[0],
                                      UTCDateTime(0)+delta[0])

        xcorr, winPerInterval, \
        iStart, iEnd, wStart, \
        wEnd, sr = xcorr2(tr1,
                          tr2,
                          window_seconds=wsec,
                          interval_seconds=isec,
                          window_overlap=olap,
                          window_buffer_length=wbl,
                          apply_stacking=True)

        olines.append("============[test: {}]============\n".format(i))
        olines.append('Params: wsec: {}, isec: {}, olap: {}, wbl: {}\n'.format(wsec, isec, olap, wbl))
        olines.append('X-corr shape: {}\n'.format(xcorr.shape))
        olines.append('Windows per interval: {}\n'.format(winPerInterval))
        olines.append("Interval start- and end-times:\n".format(i))
        for s, e in zip(iStart, iEnd):
            olines.append('{} - {}\n'.format(UTCDateTime(s + tr1.stats.starttime.timestamp),
                                             UTCDateTime(e + tr1.stats.starttime.timestamp)))
        olines.append('\n')
        olines.append("Window start- and end-times:\n".format(i))
        for s, e in zip(wStart, wEnd):
            olines.append('{} - {}\n'.format(UTCDateTime(s + tr1.stats.starttime.timestamp),
                                             UTCDateTime(e + tr1.stats.starttime.timestamp)))
        olines.append('\n')

        for oline in olines: all_output.append(oline)

        if(1):
            for oline in olines:
                eline = next(elines)

                if(oline != eline):
                    assert 0, 'Output: {} does not match expected: {}'.format(oline, eline)
                # end if
            # end for
        # end if

        i += 1
    # end for

    for line in all_output: ofh.write(line)
    ofh.close()
# end func


if __name__=="__main__":
    test_stacking_window_counts()
