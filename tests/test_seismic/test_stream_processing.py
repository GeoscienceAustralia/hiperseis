#!/usr/bin/env python
"""Unit testing for stream processing functions
"""

import itertools

import numpy as np
import obspy
from unittest.mock import MagicMock

from seismic.stream_processing import zne_order, zrt_order, zerophase_resample


def test_trace_ordering():
    test_stream = obspy.Stream([obspy.Trace(np.random.rand(20)) for _ in range(3)])

    # Test ZNE ordering
    ordered = ('BHZ', 'BHN', 'BHE')
    for perm in itertools.permutations(ordered):
        for i, tr in enumerate(test_stream):
            tr.stats.channel = perm[i]
        # end for
        test_stream.traces.sort(key=zne_order)
        assert tuple(tr.stats.channel for tr in test_stream) == ordered
    # end for

    # Test ZRT ordering
    ordered = ('BHZ', 'BHR', 'BHT')
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

def test_zerophase_resampling_success(obspy_stats):
    # Test resampling trace
    mocked_resample = MagicMock(spec=obspy.Trace.resample)
    mock_trace = MagicMock(spec=obspy.Trace, data=np.arange(4), stats=obspy_stats, resample=mocked_resample)

    zerophase_resample(mock_trace, 10)

    mocked_resample.assert_called()
