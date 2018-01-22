import os
import pytest
from obspy.geodetics import locations2degrees
from obspy import read_events, UTCDateTime
from iloc_rstt.iloc_event import ILocCatalog, DELTA, stations, stations_dict,\
    ILocEvent


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
        max_dist, sta_phs_dict = e_new._farthest_station_dist()
        assert all(e_new.all_stations[DELTA] < 1.2*max_dist)  # 120% is default
        orig_stations = stations[stations['station_code'].isin(sta_phs_dict)]
        assert all(orig_stations[DELTA] <= 1.0*max_dist)

        origin = e_old.preferred_origin()
        for a in origin.arrivals:
            sta = a.pick_id.get_referred_object().waveform_id.station_code
            assert sta in sta_phs_dict
            assert sta not in e_new.all_stations['station_code']
            station = stations_dict[sta]
            dis = locations2degrees(station.latitude, station.longitude,
                                    origin.latitude, origin.longitude)
            assert dis <= max_dist
            assert a.phase in sta_phs_dict[sta]
            pick = a.pick_id.get_referred_object()
            assert pick.time == UTCDateTime(sta_phs_dict[sta][a.phase]['time'])
            assert pick.waveform_id.channel_code == \
                   sta_phs_dict[sta][a.phase]['channel']


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


def test_write(random_filename, analyst_event):
    orig_cat = read_events(analyst_event)
    catalog = ILocCatalog(analyst_event)
    xml = random_filename(ext='.xml')
    catalog.update()
    ev = catalog.events[0]
    num_picks_b4 = len(ev.event.picks)
    ev.add_dummy_picks()
    print(len(ev.event.picks))
    ev.event.write(xml, format='SC3ML', creation_info=catalog.creation_info)
    assert os.path.exists(xml)
    new_cat = read_events(xml, format='SC3ML')
    assert num_picks_b4 + 5 == len(new_cat.events[0].picks)
