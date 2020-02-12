#!/usr/bin/env python
# coding: utf-8
"""Class encapsulating a collection of event waveforms for stations of one network.
"""

import os
from collections import defaultdict, OrderedDict

import obspy

from seismic.stream_io import read_h5_stream


class NetworkEventDataset:
    """Collection of 3-channel ZNE streams with traces aligned to a fixed time window about
     seismic P-wave arrival events, for a given network.

     Two indexes are provided. One indexes hierarchically by station code and
     event ID, yielding a 3-channel ZNE stream per event, so that you can easily gather all
     traces for a given station by iterating over events.

     The other index indexes hierarchically by event ID and station code, yielding a
     3-channel ZNE stream per station. Using this index you can easily gather all traces
     for a given event across multiple stations.
    """
    def __init__(self, stream_src, network=None, station=None, location=''):
        """Initialize from data source (file or obspy.Stream). Traces are COPIED into
        the dataset in order to leave input object intact, since many obspy functions
        mutate traces in-place.
        """
        if isinstance(stream_src, obspy.Stream):
            net = network
            sta = station
            loc = location or None
            if net or sta or loc:
                data_src = stream_src.select(net, sta, loc)
            else:
                data_src = stream_src
            # end if
        elif os.path.isfile(stream_src):
            data_src = read_h5_stream(stream_src, network, station, location)
        else:
            assert False, "Unknown data source {}".format(type(stream_src))
        # end if

        # Data in data_src collects all traces together under a single Stream object.
        # In order to get control over data slicing and traceability in processing, we
        # break it down into one Stream per ZNE channel triplet of a given event.
        self.network = network
        self.db_sta = defaultdict(lambda: defaultdict(obspy.Stream))
        self.db_evid = defaultdict(dict)
        for tr in data_src:
            net, sta, _, _ = tr.id.split('.')
            if self.network:
                assert net == self.network
            else:
                self.network = net
            # end if
            # Create single copy of the trace to be shared by both dicts.
            dupe_trace = tr.copy()
            self.db_sta[sta][tr.stats.event_id].append(dupe_trace)

        # end for
        # Index same obspy.Stream instances in event dict. This way, any changes
        # to a given event stream will be seen by both indexes.
        for sta, ev_db in self.db_sta.items():
            for evid, stream in ev_db.items():
                self.db_evid[evid][sta] = stream
            # end for
        # end for

        # Sort each stream into ZNE order.
        for ev_db in self.db_sta.values():
            for stream in ev_db.values():
                stream.sort(keys=['channel'], reverse=True)
            # end for
        # end for

        # Order station keys and event ids to ensure consistent ordering in slicing
        # results. Aids reporting and traceability.
        self.db_sta = OrderedDict(sorted([(sta, OrderedDict(sorted(sta_db.items(), key=lambda k: k[0])))
                                          for sta, sta_db in self.db_sta.items()], key=lambda k: k[0]))
        self.db_evid = OrderedDict(sorted([(evid, OrderedDict(sorted(ev_db.items(), key=lambda k: k[0])))
                                           for evid, ev_db in self.db_evid.items()], key=lambda k: k[0]))

    # end func


    def by_station(self):
        """
        Iterate over station sub-dictionaries
        :return: Iterable over the stations, each element consisting of pair containing
            (station code, event dict).
        """
        return iter(self.db_st.items())
    # end func


    def by_event(self):
        """
        Iterate over event sub-dictionaries
        :return: Iterable over the discrete events, each element consisting of pair containing
            (event id, station dict).
        """
        return iter(self.db_evid.items())
    # end func


    def prune(self, items):
        """
        Remove a given sequence of (station, evnet) pairs from the dataset.

        :param items: Iterable of (station, event) pairs
        :return: None
        """
        for station, event_id in items:
            self.db_sta[station].pop(event_id)
            self.db_evid[event_id].pop(station)
        # end for
    # end func

# end class


if __name__ == "__main__":
    # src_file = (r"/g/data/ha3/am7399/shared/OA_RF_analysis/" +
    #             r"OA_event_waveforms_for_rf_20170911T000036-20181128T230620_rev8.h5")
    # test_db = NetworkEventDataset(src_file, network='OA', station='BT23', location='0M')
    pass
# end if
