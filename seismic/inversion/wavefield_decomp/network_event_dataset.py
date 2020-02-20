#!/usr/bin/env python
# coding: utf-8
"""Class encapsulating a collection of event waveforms for stations of one network.
"""

import os
import sys

from sortedcontainers import SortedDict
import obspy

from seismic.stream_io import read_h5_stream
from seismic.receiver_fn.rf_util import zne_order, zrt_order


class NetworkEventDataset:
    """Collection of 3-channel ZNE streams with traces aligned to a fixed time window about
     seismic P-wave arrival events, for a given network.

     Two indexes are provided. One indexes hierarchically by station code and
     event ID, yielding a 3-channel ZNE stream per event, so that you can easily gather all
     traces for a given station by iterating over events.

     The other index indexes hierarchically by event ID and station code, yielding a
     3-channel ZNE stream per station. Using this index you can easily gather all traces
     for a given event across multiple stations.

     Preferably each input trace will already have an 'event_id' attribute in its stats. If
     not, an event ID will be invented based on station identifiers and time window.
    """
    def __init__(self, stream_src, network=None, station=None, location='', ordering='ZNE'):
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

        self.network = network

        # Data in data_src collects all traces together under a single Stream object.
        # In order to get control over data slicing and traceability in processing, we
        # break it down into one Stream per ZNE channel triplet of a given event.
        self.db_sta = SortedDict()
        for tr in data_src:
            net, sta, loc, _ = tr.id.split('.')
            if self.network:
                assert net == self.network
            else:
                self.network = net
            # end if
            # Create single copy of the trace to be shared by both dicts.
            dupe_trace = tr.copy()
            try:
                event_id = tr.stats.event_id
            except AttributeError:
                event_id = '.'.join([net, sta, loc, '_'.join([str(tr.stats.starttime), str(tr.stats.endtime)])])
            # end try
            self.db_sta.setdefault(sta, SortedDict()).setdefault(event_id, obspy.Stream()).append(dupe_trace)
        # end for

        # Index same obspy.Stream instances in event dict. This way, any changes
        # to a given event stream will be seen by both indexes.
        self.db_evid = SortedDict()
        for sta, ev_db in self.db_sta.items():
            for evid, stream in ev_db.items():
                self.db_evid.setdefault(evid, SortedDict())[sta] = stream
            # end for
        # end for

        # Sort each stream into specific order.
        if ordering.upper() == 'ZNE':
            ordinal = zne_order
        elif ordering.upper() == 'ZRT':
            ordinal = zrt_order
        else:
            ordinal = None
        # end if

        if ordinal is not None:
            self.apply(lambda x: x.traces.sort(key=ordinal))
        # end if

    # end func

    def __iter__(self):
        # Flat iterator. Loops over self.db_sta depth first and returns tuple of keys and matching stream.
        # for sta, ev_db in self.db_sta.items():
        #     for evid, stream in ev_db.items():
        #         yield (sta, evid, stream)
        #     # end for
        # # end for
        return ((sta, evid, stream) for sta, ev_db in self.db_sta.items() for evid, stream in ev_db.items())
    # end if

    def __len__(self):
        # Returns number of streams
        return sum((len(x) for x in self.db_sta.values()))
    # end func

    def __repr__(self):
        # Displays summary string for all streams
        return '\n'.join((evid + ', ' + str(stream) for _, evid, stream in iter(self)))
    # end func

    def station(self, station_code):
        """
        Accessor for events for a given station.
        :param station_code: Station to get
        :return: Event index for station, if station is found
        :rtype: SortedDict
        """
        return self.db_sta.get(station_code)
    # end func

    def event(self, event_id):
        """
        Accessor for stations for a given event.
        :param event_id: ID of event to look up
        :return: Station index for given event, if event ID is found
        :rtype: SortedDict
        """
        return self.db_evid.get(event_id)
    # end func

    def curate(self, curator):
        """
        Curate the dataset according to a callable curator. Curator call signature takes station code,
        event id and stream as input, and returns boolean whether to keep Stream or not.

        :param curator: Function or callable delegate to adjudicate whether to keep each given stream.
        :type curator: Callable
        :return: None
        """
        # Only need to loop over one db, since they both reference the same underlying Stream instances.
        PY2 = (sys.version_info[0] == 2)

        if PY2:
            from itertools import ifilterfalse as filterfalse
        else:
            from itertools import filterfalse
        # end if

        discard_items = [(x[0], x[1]) for x in filterfalse(lambda rec: curator(*rec), iter(self))]

        self.prune(discard_items)
    # end func

    def apply(self, callable):
        # Apply a callable across all streams. Use to apply uniform processing steps
        # to the whole dataset.
        for _1, _2, stream in iter(self):
            callable(stream)
    # end func

    def by_station(self):
        """
        Iterate over station sub-dictionaries
        :return: Iterable over the stations, each element consisting of pair containing
            (station code, event dict).
        """
        return iter(self.db_sta.items())
    # end func

    def by_event(self):
        """
        Iterate over event sub-dictionaries
        :return: Iterable over the discrete events, each element consisting of pair containing
            (event id, station dict).
        """
        return iter(self.db_evid.items())
    # end func

    def prune(self, items, cull=True):
        """
        Remove a given sequence of (station, event) pairs from the dataset.

        :param items: Iterable of (station, event) pairs
        :return: None
        """
        for station, event_id in items:
            self.db_sta[station].pop(event_id)
            self.db_evid[event_id].pop(station)
            if cull:
                if not self.db_sta[station]:
                    self.db_sta.pop(station)
                # end if
                if not self.db_evid[event_id]:
                    self.db_evid.pop(event_id)
                # end if
            # end if
        # end for

    # end func

# end class


if __name__ == "__main__":
    # from seismic.stream_quality_filter import curate_stream3c
    #
    # src_file = (r"/g/data/ha3/am7399/shared/OA_RF_analysis/" +
    #             r"OA_event_waveforms_for_rf_20170911T000036-20181128T230620_rev8.h5")
    # test_db = NetworkEventDataset(src_file, network='OA', station='BT23', location='0M')
    # print(len(test_db))
    # print(test_db)
    # for sta, _ in test_db.by_station():
    #     print(sta)
    # for evid, _ in test_db.by_event():
    #     print(evid)
    # test_db.prune((('BT23', 'smi:ISC/evid=615353118'),))
    # for sta, evid, stream in test_db:
    #     print(sta, evid)
    # test_db.curate(lambda *x: curate_stream3c(x[1], x[2]))
    # test_db.apply(lambda x: x.rotate('NE->RT'))
    # print(len(test_db))
    # print(test_db)
    pass
# end if
