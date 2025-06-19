#!/usr/bin/env python
"""Unit testing for seismic stream quality management.
"""

import numpy as np
import obspy

from seismic.stream_quality_filter import curate_stream3c


def test_curate_stream():
    npts = 300
    starttime = obspy.UTCDateTime('2020-02-05T15:00:00')
    stats = obspy.core.Stats({'sampling_rate': 10, 'npts': npts, 'inclination': 22,
                              'starttime': starttime, 'onset': starttime + 5})
    clean_stream = obspy.Stream()
    for channel in ['BHZ', 'BHN', 'BHE']:
        stats.update({'channel': channel})
        data = np.random.rand(npts)
        trace = obspy.Trace(data, stats)
        assert trace.stats.endtime > trace.stats.starttime
        clean_stream += trace
    # end for

    # Test that clean stream passes checks
    test_id = 'test'
    test_stream = clean_stream.copy()
    assert curate_stream3c(test_id, test_stream)

    # Test invalid inclination
    test_stream[0].stats.inclination = np.NaN
    assert not curate_stream3c(test_id, test_stream)

    # Test not enough streams
    test_stream = clean_stream.copy()
    test_stream.remove(test_stream[0])
    assert not curate_stream3c(test_id, test_stream)

    # Test making trace time spans consistent
    test_stream = clean_stream.copy()
    assert test_stream[2].stats.starttime == test_stream[1].stats.starttime == test_stream[0].stats.starttime
    assert test_stream[2].stats.endtime == test_stream[1].stats.endtime == test_stream[0].stats.endtime
    test_stream[0].trim(test_stream[0].stats.starttime, test_stream[0].stats.endtime - 5.0)
    assert test_stream[2].stats.endtime == test_stream[1].stats.endtime
    assert test_stream[1].stats.endtime != test_stream[0].stats.endtime
    assert curate_stream3c(test_id, test_stream)
    assert test_stream[2].stats.starttime == test_stream[1].stats.starttime == test_stream[0].stats.starttime
    assert test_stream[2].stats.endtime == test_stream[1].stats.endtime == test_stream[0].stats.endtime
    # Repeat for inconsistent starttimes
    test_stream[0].trim(test_stream[0].stats.starttime + 2.0, test_stream[0].stats.endtime)
    assert test_stream[2].stats.starttime == test_stream[1].stats.starttime
    assert test_stream[1].stats.starttime != test_stream[0].stats.starttime
    assert curate_stream3c(test_id, test_stream)
    assert test_stream[2].stats.starttime == test_stream[1].stats.starttime == test_stream[0].stats.starttime
    assert test_stream[2].stats.endtime == test_stream[1].stats.endtime == test_stream[0].stats.endtime

    # Test nan in channel data
    for i in range(3):
        test_stream = clean_stream.copy()
        assert curate_stream3c(test_id, test_stream)
        test_stream[i].data[10 + i] = np.NaN
        assert not curate_stream3c(test_id, test_stream)
    # end for

    # Test invariant channel data
    for i in range(3):
        test_stream = clean_stream.copy()
        assert curate_stream3c(test_id, test_stream)
        test_stream[i].data[:] = float(i)
        assert not curate_stream3c(test_id, test_stream)
    # end for

# end func


if __name__ == "__main__":
    test_curate_stream()
# end if
