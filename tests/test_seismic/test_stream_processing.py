#!/usr/bin/env python
"""Unit testing for stream processing functions
"""

import itertools

import numpy as np
import obspy

from seismic.stream_processing import zne_order, zrt_order


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


if __name__ == "__main__":
    test_trace_ordering()
# end if
