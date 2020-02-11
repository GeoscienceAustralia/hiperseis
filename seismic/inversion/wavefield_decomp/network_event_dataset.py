#!/usr/bin/env python
# coding: utf-8
"""Class encapsulating a collection of event waveforms for stations of one network.
"""

import os
from collections import defaultdict, OrderedDict

import obspy

from seismic.stream_io import read_h5_stream
from seismic.receiver_fn.rf_util import zne_order


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
        self.db_evid = defaultdict(lambda: defaultdict(obspy.Stream))
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
            self.db_evid[tr.stats.event_id][sta].append(dupe_trace)
        # end for

        # Sort each stream into ZNE order.
        for ev_db in self.db_sta.values():
            for stream in ev_db.values():
                stream.sort(keys=['channel'], reverse=True)
            # end for
        # end for
        for st_db in self.db_evid.values():
            for stream in st_db.values():
                stream.sort(keys=['channel'], reverse=True)
            # end for
        # end for

        # Order station keys and event ids to ensure consistent ordering in slicing
        # results. Aids reporting and traceability.
        self.db_sta = OrderedDict(sorted([(sta, OrderedDict(sorted(sta_db.items(), key=lambda k: k[0]))) for sta, sta_db in self.db_sta], key=lambda k: k[0]))
        self.db_evid = OrderedDict(sorted([(evid, OrderedDict(sorted(ev_db.items(), key=lambda k: k[0]))) for evid, ev_db in self.db_evid], key=lambda k: k[0]))

    # end func


    def by_station(self):
        return iter(self.db_st.items())
    # end func


    def by_event(self):
        return iter(self.db_evid.items())
    # end func

# end class
