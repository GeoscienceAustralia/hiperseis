import pytest
from obspy.geodetics import locations2degrees
from iloc_rstt.iloc_event import ILocCatalog, DELTA, stations, stations_dict

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
        assert all(e_new.new_stations[DELTA] < 1.2*max_dist)  # 120% is default
        orig_stations = stations[stations['station_code'].isin(orig_stas)]
        assert all(orig_stations[DELTA] <= 1.0*max_dist)

        origin = e_old.preferred_origin()
        for a in origin.arrivals:
            sta = a.pick_id.get_referred_object().waveform_id.station_code
            assert sta in orig_stas
            assert sta not in e_new.new_stations['station_code']
            station = stations_dict[sta]
            dis = locations2degrees(station.latitude, station.longitude,
                                    origin.latitude, origin.longitude)
            assert dis <= max_dist
