from __future__ import absolute_import, print_function
import os
from subprocess import check_call
import pytest
import random
import string
from obspy.geodetics import locations2degrees
from obspy import read_events, UTCDateTime
from obspy.core.event import (ResourceIdentifier, WaveformStreamID, Pick,
                              Arrival)
from iloc_rstt.iloc_event import ILocCatalog, DELTA, stations, stations_dict,\
    ILocEvent, STATION_CODE, NETWORK_CODE, DBFLAG

try:
    cmd = 'seiscomp check'.split()
    check_call(cmd)
    SC3 = True
except OSError:
    SC3 = False


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


class ILocEventDummy(ILocEvent):
    def __init__(self, event):
        super(ILocEventDummy, self).__init__(
            event=event
        )

    def add_dummy_picks(self):
        """
        Don't call this function for a real job
        """
        import random
        for i, r in enumerate(self.all_stations.iterrows()):
            d = r[1].to_dict()  # pandas series to dict
            res = ResourceIdentifier('custom_pick_{}'.format(i))

            wav_id = WaveformStreamID(station_code=d[STATION_CODE],
                                      network_code=d[NETWORK_CODE],
                                      channel_code='BHN')
            # randomly choose a time from the available picks
            p = Pick(resource_id=res, waveform_id=wav_id,
                     phase_hint='S', time=random.choice(self.picks).time)
            self.picks.append(p)  # add to picks list
            a = Arrival(pick_id=res, phase='S')
            self._pref_origin.arrivals.append(a)
            if i == 4:
                break


class ILocCatalogDummy(ILocCatalog):
    def __init__(self, event_xml):
        super(ILocCatalogDummy, self).__init__(
            event_xml=event_xml
        )

    def update(self):
        for e in self.orig_events:
            iloc_ev = ILocEventDummy(e)()
            self.events.append(iloc_ev)


def test_iloc_catalog_write(random_filename, analyst_event):
    orig_cat = read_events(analyst_event)
    orig_evt = orig_cat.events[0]
    catalog = ILocCatalogDummy(analyst_event)
    xml = random_filename(ext='.xml')
    catalog.update()
    ev = catalog.events[0]
    ev.add_dummy_picks()
    catalog.write(xml, format='SC3ML', creation_info=catalog.creation_info)
    assert os.path.exists(xml)
    new_cat = read_events(xml, format='SC3ML')
    assert len(new_cat.events) == 1
    new_ev = new_cat.events[0]
    assert len(new_ev.picks) >= len(orig_cat.events[0].picks)

    assert len(new_ev.magnitudes) >= len(orig_cat.events[0].magnitudes)
    new_cat = read_events(xml, format='SC3ML')
    new_evt = new_cat.events[0]
    assert len(orig_evt.picks) + 5 == len(new_evt.picks)
    assert len(new_evt.origins) == len(orig_evt.origins)
    assert len(new_evt.origins[0].arrivals) == len(new_evt.picks)


@pytest.mark.skipif(not SC3, reason='Skipped as seiscomp3 is not installed')
def test_sc3_db_write(random_filename, analyst_event, test_dir):

    orig_cat = read_events(analyst_event)
    orig_evt = orig_cat.events[0]
    catalog = ILocCatalogDummy(analyst_event)
    catalog.update()
    ev = catalog.events[0]
    ev.add_dummy_picks()
    ev = catalog.events[0]

    sqldb = 'sqlite3://' + random_filename(ext='.db')  # sqlite db

    check_call(['bash', os.path.join(test_dir, 'setupdb.sh'),
                sqldb,
                random_filename(ext='.xml'),  # inventory xml
                random_filename(ext='.xml')  # config xml
                ])

    catalog.insert_into_sc3db(dbflag=sqldb, plugins='dbsqlite3')

    cmd = 'scevtls -d'.split()
    cmd.append(DBFLAG)
    cmd += ('| grep {}'.format(ev._event.resource_id)).split()
    check_call(cmd)
