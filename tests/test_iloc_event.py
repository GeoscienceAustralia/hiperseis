import os
import pytest
from obspy.geodetics import locations2degrees
from obspy import read_events
from iloc_rstt.iloc_event import ILocCatalog, DELTA, stations, stations_dict


@pytest.fixture(params=[
    pytest.lazy_fixture('xml'),
    pytest.lazy_fixture('analyst_event')
], scope='module')
def one_event(request):
    return request.param


def test_stations_in_range(one_event):
    catalog = ILocCatalog(one_event)
    assert len(catalog.orig_events) == 1  # there is only one event in the xml
    catalog.update()
    assert len(catalog.events) == 1  # should create only one iloc event

    for e_old, e_new in zip(catalog.orig_events, catalog.events):
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


def test_iloc_catalog_write(random_filename, analyst_event):
    orig_cat = read_events(analyst_event)
    catalog = ILocCatalog(analyst_event)
    xml = random_filename(ext='.xml')
    catalog.update()
    catalog.write(xml, format='SC3ML', creation_info=catalog.creation_info)
    assert os.path.exists(xml)
    new_cat = read_events(xml, format='SC3ML')
    assert len(new_cat.events) == 1
    ev = new_cat.events[0]
    assert len(ev.picks) <= len(orig_cat.events[0].picks)
    assert len(ev.magnitudes) <= len(orig_cat.events[0].magnitudes)
