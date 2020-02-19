#!/usr/bin/env python
"""Unit testing for class NetworkEventDataset.
"""

import os
# import itertools

import numpy as np

import obspy

from seismic.inversion.wavefield_decomp.network_event_dataset import NetworkEventDataset


def _mock_test_stream():
    network = 'AU'
    stations = ('TE01', 'TE02', 'TE03', 'TE04')
    channels = ('HHZ', 'HHN', 'HHE')
    # n_stations = len(stations)
    np.random.seed(19460706)
    ref_time = obspy.UTCDateTime.now() - 30*24*3600
    events = tuple(('0' + str(i), ref_time + i*24*3600) for i in range(10))
    test_stream = obspy.Stream()
    for station in np.random.permutation(stations):
        sta_events = np.random.permutation(events)[:np.random.randint(7, 11)]
        for event_id, event_time in sta_events:
            # Base signal bias on event id just for testing purposes,
            # so we can tell from event roughly what amplitude should be.
            trace_bias = int(event_id)
            station_event_time = event_time + np.random.uniform(-60, 60)
            for channel in np.random.permutation(channels):
                test_trace = obspy.Trace(np.random.randn(100) + trace_bias)
                test_trace.stats.update(
                    {'network': network,
                     'station': station,
                     'channel': channel,
                     'starttime': station_event_time,
                     'event_id': event_id}
                )
                assert test_trace.stats.endtime > test_trace.stats.starttime
                test_stream += test_trace
            # end for
        # end for
    # end for
    return test_stream
# end func


def _mock_test_file(test_stream, tmp_path):
    outfile = os.path.join(tmp_path, 'test_stream.h5')
    test_stream.write(outfile, 'h5')
    assert os.path.isfile(outfile)
    return outfile
# end func


def test_network_event_dataset(tmpdir):
    # Create test datasets, one as obspy.Stream and other as file.
    test_stream = _mock_test_stream()
    test_file = _mock_test_file(test_stream, tmpdir)

    def do_test_ned(test_ned):
        # Print class instance to exercise __repr__

        # Check len of class instance returns cumulative sum of events for all stations

        # Test that loaded streams are in ZNE order

        # Test that event ids are ordered

        # Test that station keys are ordered

        # Test class iteration looks over stations first then event per station

        # Test accessing all events for given station by station code

        # Test accessing all stations for given event by event id

        # Test pruning of streams

        # Test dataset curation

        # Test apply function

        # Test iteration by station code

        # Test iteration by event id

        pass
    # end func

    do_test_ned(NetworkEventDataset(test_stream))
    do_test_ned(NetworkEventDataset(test_file))

# end func

if __name__ == "__main__":
    import tempfile
    test_network_event_dataset(tempfile.mkdtemp())
# end if
