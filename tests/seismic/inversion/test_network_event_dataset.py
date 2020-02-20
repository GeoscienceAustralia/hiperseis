#!/usr/bin/env python
"""Unit testing for class NetworkEventDataset.
"""

import os

import numpy as np
import obspy

from seismic.inversion.wavefield_decomp.network_event_dataset import NetworkEventDataset


def _mock_test_stream():
    network = 'AU'
    stations = ('TE01', 'TE02', 'TE03', 'TE04')
    channels = ('HHZ', 'HHN', 'HHE')
    np.random.seed(19460706)
    ref_time = obspy.UTCDateTime.now() - 30*24*3600
    events = tuple(('0' + str(i), ref_time + i*24*3600) for i in range(10))
    test_stream = obspy.Stream()
    num_streams = 0
    for station in np.random.permutation(stations):
        # Randomly ignore up to a few events for each station
        sta_events = np.random.permutation(events)[:np.random.randint(7, 11)]
        for event_id, event_time in sta_events:
            # Base signal bias on event id just for testing purposes,
            # so we can tell from event roughly what amplitude should be.
            trace_bias = int(event_id)
            station_event_time = event_time + np.random.uniform(-60, 60)
            for channel in np.random.permutation(channels):
                test_trace = obspy.Trace(np.random.randn(100) + trace_bias)
                # Note: we need str conversions here because the use of numpy.random.permutations
                # changes the type of strings from str to numpy.string_, which is not compatible
                # with obspyh5.
                test_trace.stats.update(
                    {'network': network,
                     'station': str(station),
                     'channel': str(channel),
                     'starttime': str(station_event_time),
                     'event_id': event_id}
                )
                assert test_trace.stats.endtime > test_trace.stats.starttime
                test_stream += test_trace
            # end for
            num_streams += 1
        # end for
    # end for
    return test_stream, num_streams
# end func


def _mock_test_file(test_stream, tmp_path):
    # Cast to string needed here because some versions of pytest pass a LocalPath
    # instance for the tmpdir fixture which won't implicitly convert to string.
    outfile = os.path.join(str(tmp_path), 'test_stream.h5')
    test_stream.write(outfile, 'h5')
    assert os.path.isfile(outfile)
    return outfile
# end func


def test_network_event_dataset(tmpdir):
    # Create test datasets, one as obspy.Stream and other as file.
    test_stream, num_streams = _mock_test_stream()
    test_file = _mock_test_file(test_stream, tmpdir)

    def do_test_ned(test_ned):
        # Tests expected behaviour and attributed of class NetworkEventDataset.

        # Print class instance to exercise __repr__
        print(test_ned)

        # Check len of class instance returns cumulative sum of events for all stations
        assert len(test_ned) == num_streams

        # Test that station keys are ordered
        assert list(test_ned.db_sta.keys()) == sorted(test_ned.db_sta.keys())

        # Test that event ids are ordered
        assert list(test_ned.db_evid.keys()) == sorted(test_ned.db_evid.keys())

        # Test iteration by station code
        expected_stations = ['TE01', 'TE02', 'TE03', 'TE04']
        assert [sta for sta, _ in test_ned.by_station()] == expected_stations

        # Test iteration by event id
        expected_events = ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09']
        assert [evt for evt, _ in test_ned.by_event()] == expected_events

        # Test that instance identities of streams are same regardless of access method
        stream_ids_by_sta = set(id(stream) for _, ev_db in test_ned.by_station() for stream in ev_db.values())
        assert len(stream_ids_by_sta) == num_streams
        stream_ids_by_evt = set(id(stream) for _, st_db in test_ned.by_event() for stream in st_db.values())
        assert stream_ids_by_sta == stream_ids_by_evt

        # Test class iteration loops over stations first then event per station
        # Test that loaded streams are in ZNE order
        last_evid = None
        sta_next_count = 0
        last_sta = None
        for sta, evid, stream in test_ned:
            # Event id should change every iteration
            assert evid != last_evid
            last_evid = evid
            # Count how often station changes
            sta_next_count += int(sta != last_sta)
            last_sta = sta
            # Check ordering
            assert len(stream) == 3
            assert stream[0].stats.channel[-1] == 'Z'
            assert stream[1].stats.channel[-1] == 'N'
            assert stream[2].stats.channel[-1] == 'E'
        # end for
        assert sta_next_count == len(expected_stations)

        # Test accessing all events for given station by station code
        num_by_station = 0
        for sta in expected_stations:
            ev_db = test_ned.station(sta)
            assert np.all([evid in expected_events for evid in ev_db.keys()])
            num_by_station += len(ev_db)
        # end for
        assert num_by_station == num_streams
        assert test_ned.station('NON_EXISTENT') == None

        # Test accessing all stations for given event by event id
        num_by_event = 0
        for evid in expected_events:
            sta_db = test_ned.event(evid)
            assert np.all([sta in expected_stations for sta in sta_db.keys()])
            num_by_event += len(sta_db)
        # end for
        assert num_by_event == num_streams
        assert test_ned.event('NON_EXISTENT') == None

        # Test pruning of streams
        discard = (('TE03', '07'), ('TE01', '05'), ('TE02', '00'), ('TE04', '03'))
        test_ned.prune(discard)
        assert len(test_ned) == num_streams - 4

        # Test dataset curation
        def discard_low_amp(_1, _2, stream):
            return np.mean(stream[0].data) > 4.5
        # end if
        assert not np.all([stream[0].data.mean() > 4.5 for _1, _2, stream in test_ned])
        test_ned.curate(discard_low_amp)
        assert len(test_ned) < num_streams//2
        assert np.all([stream[0].data.mean() > 4.5 for _1, _2, stream in test_ned])

        # Test apply function
        def filter_traces(stream):
            for i in range(3):
                stream[i].data = stream[i].data.mean() * np.ones_like(stream[i].data)
        # end func
        assert np.all([tr.data.std() > 0 for _1, _2, stream in test_ned for tr in stream])
        test_ned.apply(filter_traces)
        assert np.all([np.isclose(tr.data.std(), 0) for _1, _2, stream in test_ned for tr in stream])

    # end func

    do_test_ned(NetworkEventDataset(test_stream))
    do_test_ned(NetworkEventDataset(test_file))

# end func

if __name__ == "__main__":
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        test_network_event_dataset(tmpdir)
# end if
