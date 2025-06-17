#!/usr/bin/env python
"""Unit testing for stream processing functions
"""

import itertools

import numpy as np
import obspy
from unittest.mock import MagicMock, patch

from seismic.stream_processing import zne_order, zrt_order, zerophase_resample


def test_trace_ordering():
    test_stream = obspy.Stream([obspy.Trace(np.random.rand(20)) for _ in range(4)])

    # Test ZNE ordering
    ordered = ('BHZ', 'BHN', 'BHE', 'BHY')
    for perm in itertools.permutations(ordered):
        for i, tr in enumerate(test_stream):
            tr.stats.channel = perm[i]
        # end for
        test_stream.traces.sort(key=zne_order)
        assert tuple(tr.stats.channel for tr in test_stream) == ordered
    # end for

    # Test ZRT ordering
    ordered = ('BHZ', 'BHR', 'BHT', 'BHY')
    for perm in itertools.permutations(ordered):
        for i, tr in enumerate(test_stream):
            tr.stats.channel = perm[i]
        # end for
        test_stream.traces.sort(key=zrt_order)
        assert tuple(tr.stats.channel for tr in test_stream) == ordered
    # end for

# end func

def test_zerophase_resampling_with_invalid_types():
    # Test invalid item, not Stream or Trace
    try:
        zerophase_resample(123, 10)
    except TypeError:
        pass
    else:
        raise AssertionError("Expected TypeError for invalid item type")

@patch('seismic.stream_processing.lowpass')
def test_zerophase_resampling_success(mocked_lowpass, obspy_stats):
    # Test trace gets resampled and lowpass is called if resample_hz < sampling_rate
    mocked_resample = MagicMock(spec=obspy.Trace.resample)
    mock_trace = MagicMock(spec=obspy.Trace, data=np.array([1,2,3,4]), stats=obspy_stats, resample=mocked_resample)

    zerophase_resample(mock_trace, 1)

    mocked_resample.assert_called()
    mocked_lowpass.assert_called()
