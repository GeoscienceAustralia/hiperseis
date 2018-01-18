import os
import pytest
from iloc_rstt.iloc_event import ILocCatalog, DELTA, stations


@pytest.fixture
def analyst_event(test_dir):
    return os.path.join(test_dir, 'mocks', 'events',
                        'analyst_event_samples', '772009.xml')


def test_stations_in_range(analyst_event):
    cat = ILocCatalog(analyst_event)
    assert len(cat.events) == 1  # there is only one event in the xml
    cat.update()
    assert len(cat.iloc_events) == 1  # should create only one iloc event

    for e_old, e_new in zip(cat.events, cat.iloc_events):
        max_dist, orig_stas = e_new._farthest_station_dist()
        assert all(e_new.new_stations[DELTA] < 1.2*max_dist)
        orig_stations = stations[stations['station_code'].isin(orig_stas)]
        assert all(orig_stations[DELTA] <= 1.0*max_dist)

