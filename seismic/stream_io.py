#!/usr/bin/env python
"""Helper functions for seismic stream IO.
"""

import os
import itertools
import logging
import numpy as np
from sortedcontainers import SortedDict, SortedSet, SortedList
import obspy
from obspy.taup import TauPyModel
from obspy.io.sac.sactrace import SACTrace
import h5py
import obspyh5
from obspyh5 import dataset2trace, is_obspyh5, trace2group
from os.path import splitext

from seismic.units_utils import KM_PER_DEG
from rf.rfstream import rfstats, obj2stats
from rf.util import _get_stations

# pylint: disable=invalid-name


logging.basicConfig()


# Source event centric indexing format
EVENTIO_TF = '.datetime:%Y-%m-%dT%H:%M:%S'
EVENTIO_H5INDEX = (
    'waveforms/{network}.{station}.{location}/{event_time%s}/' % EVENTIO_TF +
    '{channel}_{starttime%s}_{endtime%s}' % (EVENTIO_TF, EVENTIO_TF)
    )

def safe_iter_event_data(events, inventory, get_waveforms, use_rfstats=True, phase='P',
                    request_window=None, pad=10, pbar=None, **kwargs):
    """
    Return iterator yielding three component streams per station and event.

    :param events: list of events or `~obspy.core.event.Catalog` instance
    :param inventory: `~obspy.core.inventory.inventory.Inventory` instance
        with station and channel information
    :param get_waveforms: Function returning the data. It has to take the
        arguments network, station, location, channel, starttime, endtime.
    :param phase: Considered phase, e.g. 'P', 'S', 'PP'
    :type request_window: tuple (start, end)
    :param request_window: requested time window around the onset of the phase
    :param float pad: add specified time in seconds to request window and
       trim afterwards again
    :param pbar: tqdm_ instance for displaying a progressbar
    :param kwargs: all other kwargs are passed to `~rf.rfstream.rfstats()`

    :return: three component streams with raw data

    Example usage with progressbar::

        from tqdm import tqdm
        from rf.util import iter_event_data
        with tqdm() as t:
            for stream3c in iter_event_data(*args, pbar=t):
                do_something(stream3c)

    .. _tqdm: https://pypi.python.org/pypi/tqdm
    """
    from rf.rfstream import rfstats, RFStream
    method = phase[-1].upper()
    if request_window is None:
        request_window = (-50, 150) if method == 'P' else (-100, 50)
    stations = _get_stations(inventory)
    if pbar is not None:
        pbar.total = len(events) * len(stations)
    for event, seedid in itertools.product(events, stations):
        if pbar is not None:
            pbar.update(1)
        origin_time = (event.preferred_origin() or event.origins[0])['time']
        try:
            args = (seedid[:-1] + stations[seedid], origin_time)
            coords = inventory.get_coordinates(*args)
        except Exception:  # station not available at that time
            continue

        stats = None
        if(use_rfstats):
            try:
                stats = rfstats(station=coords, event=event, phase=phase, **kwargs)
            except Exception as ex:
                from warnings import warn
                warn('Error "%s" in rfstats call for event %s, station %s.'
                     % (ex, event.resource_id, seedid))
                continue
            if not stats:
                continue
        else:
            stats = obj2stats(event, coords)
        # end if

        net, sta, loc, cha = seedid.split('.')

        if(use_rfstats):
            starttime = stats.onset + request_window[0]
            endtime = stats.onset + request_window[1]
        else:
            starttime = origin_time + request_window[0]
            endtime = origin_time + request_window[1]
        # end if

        kws = {'network': net, 'station': sta, 'location': loc,
               'channel': cha, 'starttime': starttime - pad,
               'endtime': endtime + pad}
        try:
            stream = get_waveforms(**kws)
            stream.trim(starttime, endtime)
            stream.merge()
        except Exception:  # no data available
            continue

        if len(stream) != 3:
            from warnings import warn
            warn('Need 3 component seismograms. %d components '
                 'detected for event %s, station %s.'
                 % (len(stream), event.resource_id, seedid))
            continue
        if any(isinstance(tr.data, np.ma.masked_array) for tr in stream):
            from warnings import warn
            warn('Gaps or overlaps detected for event %s, station %s.'
                 % (event.resource_id, seedid))
            continue

        for tr in stream:
            tr.stats.update(stats)
        # end for

        yield RFStream(stream)
# end func

def read_h5_stream(src_file, network=None, station=None, loc='', root='/waveforms'):
    """Helper function to load stream data from hdf5 file saved by obspyh5 HDF5 file IO.
    Typically the source file is generated using `extract_event_traces.py` script.
    For faster loading time, a particular network and station may be specified.

    :param src_file: File from which to load data
    :type src_file: str or Path
    :param network: Specific network to load, defaults to None
    :type network: str, optional
    :param station: Specific station to load, defaults to None
    :type station: str, optional
    :param root: Root path in hdf5 file where to start looking for data, defaults to '/waveforms'
    :type root: str, optional
    :return: All the loaded data in an obspy Stream.
    :rtype: obspy.Stream
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

    stream = obspy.read(src_file, format='h5', group=group)
    return stream
# end func


def get_obspyh5_index(src_file, seeds_only=False):
    """Scrape the index (only) from an obspyh5 file.

    :param src_file: Name of file to extract index from
    :type src_file: str or pathlib.Path
    :param seeds_only: If True, only get the seed IDs of the traces. Otherwise (default), get full index.
    :type seeds_only: bool
    :return: Sorted dictionary with index of waveforms in the file
    :rtype: sortedcontainers.SortedDict
    """
    # We use SortedSet rather than SortedList for the inner container
    # because it allows set operations that are useful, such as finding
    # event IDs for events which are common across multiple stations.
    assert is_obspyh5(src_file), '{} is not an obspyh5 file'.format(src_file)
    index = SortedDict()
    with h5py.File(src_file, mode='r') as h5f:
        root = h5f['/waveforms']
        if seeds_only:
            return SortedList(root.keys())
        # end if
        for seedid, station_grp in root.items():
            for event_grp in station_grp.values():
                evid = None
                for channel in event_grp.values():
                    evid = channel.attrs['event_id']
                    break
                # end for
                if evid:
                    index.setdefault(seedid, SortedSet()).add(evid)
                # end if
            # end for
        # end for
    # end with
    return index
# end func


def iter_h5_stream(src_file, headonly=False):
    """
    Iterate over hdf5 file containing streams in obspyh5 format.

    :param src_file: Path to file to read
    :type src_file: str or pathlib.Path
    :param headonly: Only read trace stats, do not read actual time series data
    :type headonly: bool
    :yield: obspy.Stream containing traces for a single seismic event.
    """
    assert is_obspyh5(src_file), '{} is not an obspyh5 file'.format(src_file)
    logger = logging.getLogger(__name__)
    fname = os.path.split(src_file)[-1]
    with h5py.File(src_file, mode='r') as h5f:
        root = h5f['/waveforms']
        for seedid, station_grp in root.items():
            logger.info('{}: Group {}'.format(fname, seedid))
            num_events = len(station_grp)
            for i, (_src_event_time, event_grp) in enumerate(station_grp.items()):
                traces = []
                for _trace_id, channel in event_grp.items():
                    traces.append(dataset2trace(channel, headonly=headonly))
                # end for
                evid = traces[0].stats.event_id
                for tr in traces[1:]:
                    assert tr.stats.event_id == evid
                # end for
                logger.info('Event {} ({}/{})'.format(evid, i + 1, num_events))
                yield seedid, evid, obspy.Stream(traces)
            # end for
        # end for
    # end with
# end func

def write_h5_event_stream(dest_h5_file, stream, index=EVENTIO_H5INDEX, mode='a', ignore=()):
    """
    Write stream to HDF5 file in event indexed format using obspy.

    :param dest_h5_file: File in which to write the stream.
    :type dest_h5_file: str or pathlib.Path
    :param stream: The stream to write
    :type stream: obspy.Stream
    :param mode: Write mode, such as 'w' or 'a'. Use 'a' to iteratively write multiple streams to one file.
    :type mode: str
    :param ignore: List of headers to ignore when writing attributes to group. Passed on directly to obspyh5.writeh5
    :type ignore: Any iterable of str
    """
    # Lock down format of obspyh5 node layout in HDF5 to ensure compatibility with
    # custom iterators.
    assert mode.lower() != 'r', 'Write mode cannot be \'r\''
    prior_index = obspyh5._INDEX
    obspyh5.set_index(index)
    stream.write(dest_h5_file, 'H5', mode=mode, ignore=ignore)
    obspyh5.set_index(prior_index)
# end func

def sac2hdf5(src_folder, basenames, channels, dest_h5_file, tt_model_id='iasp91'):
    """
    Convert collection of SAC files from a folder into a single HDF5 stream file.

    :param src_folder: Path to folder containing SAC files
    :type src_folder: str or Path
    :param basenames: List of base filenames (file name excluding extension) to load.
    :type basenames: list of str
    :param channels: List of channels to load. For each base filename in basenames,
        there is expected to be a file with each channel as the filename extension.
    :type channels: List of str
    :param dest_h5_file: Path to output file. Will be created, or overwritten if
        already exists.
    :type dest_h5_file: str or Path
    :param tt_model_id: Which travel time earth model to use for synthesizing
        trace metadata. Must be known to obspy.taup.TauPyModel
    :type tt_model_id: str
    :return: None
    """
    tt_model = TauPyModel(tt_model_id)
    traces = []
    for basename, channel in itertools.product(basenames, channels):
        fname = os.path.join(src_folder, '.'.join([basename, channel]))
        channel_stream = obspy.read(fname, 'sac')
        tr = channel_stream[0]
        event_depth_km = tr.stats.sac['evdp']
        dist_deg = tr.stats.sac['dist'] / KM_PER_DEG
        arrivals = tt_model.get_travel_times(event_depth_km, dist_deg, ('P',))
        arrival = arrivals[0]
        inc = arrival.incident_angle
        slowness = arrival.ray_param_sec_degree
        src_dic = tr.stats.sac
        sac_tr = SACTrace(nzyear=src_dic['nzyear'], nzjday=src_dic['nzjday'], nzhour=src_dic['nzhour'],
                          nzmin=src_dic['nzmin'], nzsec=src_dic['nzsec'], nzmsec=src_dic['nzmsec'])
        onset = sac_tr.reftime
        if 'nevid' in src_dic:
            event_id = src_dic['nevid']
        else:
            event_id = basename
        # end if
        stats = {'distance': dist_deg, 'back_azimuth': src_dic['baz'], 'inclination': inc,
                 'onset': onset, 'slowness': slowness, 'phase': 'P', 'tt_model': tt_model_id,
                 'event_id': event_id}
        tr.stats.update(stats)
        traces.append(tr)
    # end for

    stream_all = obspy.Stream(traces)
    stream_all.write(dest_h5_file, 'H5')
# end func
