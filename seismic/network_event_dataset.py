#!/usr/bin/env python
# coding: utf-8
"""Class encapsulating a collection of event waveforms for stations of one network.
"""

import os
import sys

from sortedcontainers import SortedDict
import obspy

from seismic.stream_io import read_h5_stream, write_h5_event_stream
from seismic.stream_processing import zne_order, zrt_order


# pylint: disable=invalid-name


class NetworkEventDataset:
    """Collection of 3-channel ZNE streams with traces aligned to a fixed time window about
     seismic P-wave arrival events, *for a given network*.

     Two indexes are provided. One indexes hierarchically by station code and
     event ID, yielding a 3-channel ZNE stream per event, so that you can easily gather all
     traces for a given station by iterating over events.

     The other index indexes hierarchically by event ID and station code, yielding a
     3-channel ZNE stream per station. Using this index you can easily gather all traces
     for a given event across multiple stations.

     Preferably each input trace will already have an 'event_id' attribute in its stats. If
     not, an event ID will be invented based on station identifiers and time window.
    """
    def __init__(self, stream_src, network=None, station=None, location='', ordering='ZNE', root='/waveforms'):
        """
        Initialize from data source (file or obspy.Stream). Traces are COPIED into
        the dataset in order to leave input object intact, since many obspy functions
        mutate traces in-place.

        All streams in the input data source stream_src are expected to belong to the same network.
        This is checked as the data is ingested. A discrepant network code is an error condition.

        :param stream_src: Source of input streams. May be a file name or an Obspy Stream
        :type stream_src: str, pathlib.Path or obspy.Stream
        :param network: Network code of streams to load. If stream_src is an Obspy Stream, the \
            streams will be filtered to match this network code.
        :type network: str
        :param station: Station code of streams to load. If stream_src is an Obspy Stream, the \
            streams will be filtered to match this station code.
        :type station: str
        :param location: [OPTIONAL] Location code of streams to load. Leave as default (empty string) \
            if location code is empty in the data source.
        :type location: str
        :param ordering: Channel ordering to be applied to the data after loading. The channel labelling \
            must be consistent with the requested ordering - rotation to the coordinate system implied \
            by the ordering is *NOT* applied.
        :type ordering: str
        :raises AssertionError: If discrepant network code is found in input data
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
            data_src = read_h5_stream(stream_src, network, station, location, root=root)
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
        """
        Flat iterator. Loops over self.db_sta depth first and returns tuple of keys and matching stream.
        Equivalent to::

        ```Python
          for sta, ev_db in self.db_sta.items():
              for evid, stream in ev_db.items():
                  yield (sta, evid, stream)
        ```
        """
        return ((sta, evid, stream) for sta, ev_db in self.db_sta.items() for evid, stream in ev_db.items())
    # end if

    def __len__(self):
        """Returns number of streams"""
        return sum((len(x) for x in self.db_sta.values()))
    # end func

    def __repr__(self):
        """Displays summary string for all streams"""
        return '\n'.join((evid + ', ' + str(stream) for _, evid, stream in iter(self)))
    # end func

    def num_stations(self):
        """
        Get number of stations in the dataset.

        :return: Number of stations
        :rtype: int
        """
        return len(self.db_sta)
    # end func

    def station(self, station_code):
        """
        Accessor for events for a given station.

        :param station_code: Station to get
        :type station_code: str
        :return: Event index for station, if station is found
        :rtype: SortedDict
        """
        return self.db_sta.get(station_code)
    # end func

    def num_events(self):
        """
        Get number of events in the dataset.

        :return: Number of events
        :rtype: int
        """
        return len(self.db_evid)
    # end func

    def event(self, event_id):
        """
        Accessor for stations for a given event.

        :param event_id: ID of event to look up
        :type event_id: str
        :return: Station index for given event, if event ID is found, otherwise None
        :rtype: SortedDict or NoneType
        """
        return self.db_evid.get(event_id)
    # end func

    def curate(self, curator):
        """
        Curate the dataset according to a callable curator. Modifies collection in-place to remove
        streams that do not satisfy the curation criteria of the callable.
        Curator call signature must be consitent with::

            callable(station_code, event_id, stream) -> bool

        The callable returns a boolean indicating whether to keep the Stream or not.

        :param curator: Function or callable delegate to adjudicate whether to keep each given stream.
        :type curator: Callable
        :return: None
        """
        # Only need to loop over one db, since they both reference the same underlying Stream instances.
        PY2 = (sys.version_info[0] == 2)

        if PY2:
            from itertools import ifilterfalse as filterfalse  # pylint: disable=no-name-in-module, import-outside-toplevel
        else:
            from itertools import filterfalse  # pylint: disable=import-outside-toplevel
        # end if

        discard_items = [(x[0], x[1]) for x in filterfalse(lambda rec: curator(*rec), iter(self))]

        self.prune(discard_items)
    # end func

    def apply(self, _callable):
        """Apply a callable across all streams. Use to apply uniform processing steps to the whole dataset.

        :param _callable: Callable object that takes an obspy Stream as input and applies itself to that Stream. \
            Expect that stream may be mutated in-place by the callable.
        :type _callable: Any Callable compatible with the call signature.
        :return: None
        """
        for _1, _2, stream in iter(self):
            _callable(stream)
    # end func

    def by_station(self):
        """
        Iterate over station sub-dictionaries.

        :return: Iterable over the stations, each element consisting of pair containing \
            (station code, event dict).
        :rtype: Iterable(tuple)
        """
        return iter(self.db_sta.items())
    # end func

    def by_event(self):
        """
        Iterate over event sub-dictionaries.

        :return: Iterable over the discrete events, each element consisting of pair containing \
            (event id, station dict).
        :rtype: Iterable(tuple)
        """
        return iter(self.db_evid.items())
    # end func

    def prune(self, items, cull=True):
        """
        Remove a given sequence of (station, event) pairs from the dataset.

        :param items: Iterable of (station, event) pairs
        :type items: Iterable(tuple)
        :param cull: If True, then empty entries in the top level index will be removed.
        :type cull: boolean
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

    def write(self, output_h5_filename, index_format='event'):
        """
        Write event dataset back out to HDF5 file.

        :param output_h5_filename: Output file name
        :type output_h5_filename: str or path
        :param index_format: Format to use for index. Must be 'event' (default) or 'standard' (obspy default)
        :type index_format: str
        :return: True if file was written
        :rtype: boolean
        """
        assert not os.path.exists(output_h5_filename), 'Output file already exists'
        if index_format not in ['event', 'standard']:
            raise ValueError('Index format %s not supported' % index_format)
        # end if
        all_stream = obspy.Stream()
        for sta, evid, stream in iter(self):
            all_stream += stream
        # end for
        if index_format == 'event':
            write_h5_event_stream(output_h5_filename, all_stream, mode='w')
        elif index_format == 'standard':
            all_stream.write(output_h5_filename, format='H5', mode='w')
        # end if
        return os.path.isfile(output_h5_filename)
    # end func

# end class
