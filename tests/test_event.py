import os
from obspy import read_events


def test_obspy_read_sc3ml(test_dir):
    cat = read_events(os.path.join(test_dir, 'mocks', 'events',
                                   'ga2017qxlpiu.xml'),
                      format='SC3ML')
    assert len(cat.events) == 1
    event = cat.events[0]
    assert len(event.picks) == 33
    assert len(event.amplitudes) == 34
    assert len(event.station_magnitudes) == 34
    assert len(event.magnitudes) == 17
    assert len(event.origins) == 3
