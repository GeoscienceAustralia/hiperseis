"""Utility functions to help with reading RFs from H5 files.
"""

from collections import defaultdict
import logging
import os, re

import h5py
import numpy as np
import rf
from obspyh5 import dataset2trace
from rf import RFStream

from seismic.misc import print_exception
from seismic.misc_p import parallel_abort
from seismic.stream_io import get_obspyh5_index
from seismic.stream_processing import zne_order

logging.basicConfig()

def trim_hdf_keys(hdf_key_list:[str], networks_string:str, stations_string:str) -> [str]:
    """
    Trims a list of hdf_keys, filtering out unwanted networks and stations.
    @param hdf_key_list:
    @param networks_string: a space-separated list of networks. '*' includes all.
    @param stations_string: a space-separated list of stations or a text file
                            with station names in each row, w/wo location codes.
                            '*' includes all.
    @return: trimmed list
    """

    network_list = []
    station_list = []

    if(networks_string=='*'):
        network_list = []
    else:
        network_list = re.findall('\S+', networks_string)
        assert len(network_list), 'Invalid network list. Aborting..'
    # end if

    if(stations_string=='*'):
        station_list = []
    else:
        stations = []
        if(os.path.isfile(stations_string)):
            for iline, line in enumerate(open(stations_string, 'r').readlines()):
                if(not line.strip()): continue
                else: stations.append(line)
            # end for

            stations_string = ' '.join(stations)
        # end if

        station_list = re.findall('\S+', stations_string)
        assert len(station_list), 'Invalid station list. Aborting..'
    # end if

    # sanity checks
    for net in network_list:
        if net not in [hdf_key.split('.')[0] for hdf_key in hdf_key_list]:
            assert 0, 'Network {} not found in input dataset. Aborting..'.format(net)
    # end for

    for sta in station_list:
        if sta.split('.')[0] not in [hdf_key.split('.')[1] for hdf_key in hdf_key_list]:
            assert 0, 'Station {} not found in input dataset. Aborting..'.format(sta)
    # end for

    net_subset = [] # filter networks
    if(len(network_list)):
        for hdf_key in hdf_key_list:
            net, sta, loc = hdf_key.split('.')

            if(net in network_list): net_subset.append(hdf_key)
        # end for
    else:
        net_subset = hdf_key_list.copy()
    # end if

    sta_subset = [] # filter stations
    if(len(station_list)):
        for hdf_key in net_subset:
            net, sta, loc = hdf_key.split('.')

            if (sta in station_list): sta_subset.append(hdf_key)
            if ('.'.join([sta, loc]) in station_list): sta_subset.append(hdf_key)
        # end for
    else:
        sta_subset = net_subset.copy()
    # end if

    return sta_subset
# end func

class MergeInfo:
    @staticmethod
    def _validate_nsl(nsl):
        net, sta, loc = None
        try:
            net, sta, loc = nsl.split('.')
        except Exception as e:
            print_exception(e)
            parallel_abort('Invalid triplet: {}. Aborting..'.format(nsl), logger=self.log)
        # end try
    # end func

    def __init__(self, merge_fn, rf_file, log=None):
        self.merge_dict = defaultdict(list)
        self.log = log
        try:
            lines = open(merge_fn, 'r').readlines()
            lines = [line.strip() for line in lines]
            for line in lines:
                src, dest = line.split(',')
                src, dest = src.strip(), dest.strip()
                for item in (src, dest): MergeInfo._validate_nsl(item)
                self.merge_dict[dest].append(src)
            # end for
        except Exception as e:
            print_exception(e)
            parallel_abort('Failed to read {}'.format(merge_fn), logger=self.log)
        # end try

        # ensure all stations listed in merge_stations_file are available
        if (len(self.merge_dict)):
            proc_hdfkeys = get_obspyh5_index(rf_file, seeds_only=True)
            proc_hdfkeys_set = set(proc_hdfkeys)
            not_found = None
            for k, v in self.merge_dict.items():
                if k not in proc_hdfkeys_set:
                    not_found = k
                    break
                # end if

                for item in v:
                    if item not in proc_hdfkeys_set:
                        not_found = item
                        break
                    # end if
                # end for
                if (not_found is not None): break
            # end for
            if (not_found is not None):
                parallel_abort('Station {} listed in {} is not available in {}. Aborting..', logger=self.log)
            # end if
        # end if
    # end func

    def is_merged(self, nsl):
        MergeInfo._validate_nsl(nsl)
        result = False
        for k, v in self.merge_dict.items():
            if(nsl in v):
                result = True
                break
            # end if
        # end for

        return result
    # end func

    def get_surrogates(self, nsl):
        MergeInfo._validate_nsl(nsl)

        return self.merge_dict[nsl]
    # end func
# end func

def read_h5_rf(src_file, network=None, station=None, loc='', root='/waveforms'):
    """Helper function to load data from hdf5 file generated by rf library or script `rf_quality_filter.py`.
    For faster loading time, a particular network and station may be specified.

    :param src_file: File from which to load data
    :type src_file: str or Path
    :param network: Specific network to load, defaults to None
    :type network: str, optional
    :param station: Specific station to load, defaults to None
    :type station: str, optional
    :param root: Root path in hdf5 file where to start looking for data, defaults to '/waveforms'
    :type root: str, optional
    :return: All the loaded data in a rf.RFStream container.
    :rtype: rf.RFStream
    """
    logger = logging.getLogger(__name__)
    if (network is None and station is not None) or (network is not None and station is None):
        logger.warning("network and station should both be specified - IGNORING incomplete specification")
        group = root
    elif network and station:
        group = root + '/{}.{}.{}'.format(network.upper(), station.upper(), loc.upper())
    else:
        group = root
    # end if

    rf_data = rf.read_rf(src_file, format='h5', group=group)
    return rf_data
# end func

class IterRfH5FileEvents(object):
    """Helper class to iterate over events in h5 file generated by extract_event_traces.py and pass
       them to RF generator. This class avoids having to load the whole file up front via obspy which
       is slow and not scalable.

       Yields (station_id, event_id, event_time, stream)

       Due to the expected hierarchy structure of the input H5 file, yielded event traces are grouped
       by station ID.

       Data yielded per event can easily be hundreds of kB in size, depending on the length of the
       event traces. rf library defaults to window of (-50, 150) sec about the P wave arrival time.
    """

    def __init__(self, h5_filename, memmap=False, channel_pattern=None):
        """
        Initializer

        :param h5_filename: Name of file containing event seismograms in HDF5 format, indexed in \
            seismic.stream_io.EVENTIO_H5INDEX format
        :type h5_filename: str or pathlib.Path
        :param memmap: If True, memmap the open file. Can improve tractability of handling very large files.
        :type memmap: bool
        :param channel_pattern: [OPTIONAL] Ordered list of preferred channels, e.g. 'HH*,BH*'. Channel mask to use
            to select channels returned.
        :type channel_pattern: str
        """
        self.h5_filename = h5_filename
        self.num_components = 3
        self.memmap_input = memmap
        self.channel_pattern = channel_pattern

    def _open_source_file(self):
        if self.memmap_input:
            try:
                return h5py.File(self.h5_filename, 'r', driver='core', backing_store=False)
            except OSError as e:
                logger = logging.getLogger(__name__)
                logger.error("Failure to memmap input file with error:\n{}\nReverting to default driver."
                             .format(str(e)))
        # end if
        return h5py.File(self.h5_filename, 'r')
        # end if

    def __iter__(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        logger.info("Scanning jobs metadata from file {}".format(self.h5_filename))
        with self._open_source_file() as f:
            wf_data = f['waveforms']
            num_stations = len(wf_data)
            count = 0
            event_count = 0
            create_event_id = False
            first_loop = True
            for station_id in wf_data:
                count += 1
                logger.info("Station {} {}/{}".format(station_id, count, num_stations))
                station_data = wf_data[station_id]
                for event_time in station_data:
                    event_traces = station_data[event_time]
                    if not event_traces:
                        continue

                    if first_loop:
                        first_loop = False
                        tmp = list(event_traces.keys())[0]
                        create_event_id = ('event_id' not in event_traces[tmp].attrs)

                    traces = []
                    for trace_id in event_traces:
                        trace = dataset2trace(event_traces[trace_id])
                        traces.append(trace)

                    stream = RFStream(traces=traces)

                    if len(stream) != self.num_components and self.channel_pattern is not None:
                        for ch_mask in self.channel_pattern.split(','):
                            _stream = stream.select(channel=ch_mask)
                            logging.info("Tried channel mask {}, got {} channels".format(ch_mask, len(_stream)))
                            if len(_stream) == self.num_components:
                                stream = _stream
                                break
                        # end for
                    # end if

                    if len(stream) != self.num_components:
                        logging.warning("Incorrect number of traces ({}) for stn {} event {}, skipping"
                                        .format(len(stream), station_id, event_time))
                        continue
                    # end if

                    # Force order of traces to ZNE ordering.
                    stream.traces = sorted(stream.traces, key=zne_order)
                    # Strongly assert expected ordering of traces. This must be respected so that
                    # RF normalization works properly.
                    assert stream.traces[0].stats.channel[-1] == 'Z'
                    assert stream.traces[1].stats.channel[-1] == 'N'
                    assert stream.traces[2].stats.channel[-1] == 'E'

                    event_count += 1
                    if create_event_id:
                        event_id = event_count
                    else:
                        event_id = traces[0].stats.event_id
                        assert np.all([(tr.stats.event_id == event_id) for tr in traces])
                    # end if

                    yield station_id, event_id, event_time, stream
                # end for
            # end for
        # end with
        logger.info("Yielded {} event traces to process".format(event_count))
    # end func
# end class

class IterRfH5StationEvents(object):
    """Helper class to iterate over stations in h5 file generated by extract_event_traces.py and pass
       them to RF generator. This class avoids having to load the whole file up front via obspy which
       is slow and not scalable.

       This class yields 3-channel traces per station per event.

       Data yielded per station can easily be many MB in size.
    """

    def __init__(self, h5_filename, memmap=False):
        """Initializer

        :param h5_filename: Name of file containing event seismograms in HDF5 format, indexed in \
            seismic.stream_io.EVENTIO_H5INDEX format
        :type h5_filename: str or pathlib.Path
        :param memmap: If True, memmap the open file. Can improve tractability of handling very large files.
        :type memmap: bool
        """
        self.h5_filename = h5_filename
        # self.num_components = 3
        self.memmap_input = memmap

    def _open_source_file(self):
        if self.memmap_input:
            try:
                return h5py.File(self.h5_filename, 'r', driver='core', backing_store=False)
            except OSError as e:
                logger = logging.getLogger(__name__)
                logger.error("Failure to memmap input file with error:\n{}\nReverting to default driver."
                             .format(str(e)))
        # end if
        return h5py.File(self.h5_filename, 'r')
        # end if

    def __iter__(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        logger.info("Scanning station groups from file {}".format(self.h5_filename))
        with self._open_source_file() as f:
            wf_data = f['waveforms']
            num_stations = len(wf_data)
            count = 0
            event_count = 0

            for station_id in wf_data:
                count += 1
                logger.info("Station {} {}/{}".format(station_id, count, num_stations))
                station_data = wf_data[station_id]
                station_stream3c = []
                for event_time in station_data:
                    event_traces = station_data[event_time]
                    # if len(event_traces) != self.num_components:
                    #     logging.warning("Incorrect number of traces ({}) for stn {} event {}, skipping"
                    #                     .format(len(event_traces), station_id, event_time))
                    #     continue
                    traces = []

                    for trace_id in event_traces:
                        trace = dataset2trace(event_traces[trace_id])
                        traces.append(trace)

                    event_count += 1
                    station_stream3c.append(RFStream(traces=traces).sort())
                # end for

                # Yield the results with 3-channel trace triplets grouped together in RFStream instances.
                yield station_id, station_stream3c
            # end for
        # end with
        logger.info("Yielded {} event traces to process".format(event_count))
    # end func
# end class
