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

from obspy.core import Trace, Stats
from seismic.xcorqc.xcorqc import taper, whiten, zeropad_ba, xcorr2
import pytest
import numpy as np


@pytest.fixture(params=[100, 2000, 5000])
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

    result = whiten(trc, sampling_rate)

    # Whitening a null trace elevates the mean amplitude to unity. The phase
    # being zero, an inverse fft should produce a trace with the first
    # element being 1, followed by zeros.
    assert np.allclose(np.mean(result), 1./float(trace_length), atol=1e-4)
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

    tr1 = Trace(data=np.ones(trace_length*sampling_rate),
                header=Stats(header={'sampling_rate': sampling_rate,
                                     'npts': trace_length*sampling_rate,
                                     'network': 'AU',
                                     'station': 'A'}))
    tr2 = Trace(data=np.ones(trace_length*other_sampling_rate),
                header=Stats(header={'sampling_rate': other_sampling_rate,
                                     'npts': trace_length*other_sampling_rate,
                                     'network': 'AU',
                                     'station': 'B'}))

    arr, _, _, _ = xcorr2(tr1, tr2, window_seconds=trace_length/2, interval_seconds=trace_length)

    # Mean of the x-correlation function should be close to zero
    assert np.allclose(np.mean(arr), 0+0j, atol=1e-3)
# end func
