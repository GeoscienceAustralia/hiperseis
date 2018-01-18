import pytest
from iloc_rstt.iloc_event import ILocCatalog, DELTA, stations

@pytest.fixture(params=[
    pytest.lazy_fixture('xml'),
    pytest.lazy_fixture('analyst_event')
])
def one_event(request):
    return request.param


def test_stations_in_range(one_event):
    cat = ILocCatalog(one_event)
    assert len(cat.events) == 1  # there is only one event in the xml
    cat.update()
    assert len(cat.iloc_events) == 1  # should create only one iloc event

    for e_old, e_new in zip(cat.events, cat.iloc_events):
        max_dist, orig_stas = e_new._farthest_station_dist()
        assert all(e_new.new_stations[DELTA] < 1.2*max_dist)
        orig_stations = stations[stations['station_code'].isin(orig_stas)]
        assert all(orig_stations[DELTA] <= 1.0*max_dist)
