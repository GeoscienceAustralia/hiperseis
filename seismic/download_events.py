#!/usr/bin/env python
"""Use waveform database and station inventory to extract raw traces for all seismic events within a given
magnitude and time range.
"""

import os.path
import logging
from mpi4py import MPI

import warnings

warnings.simplefilter("ignore", UserWarning)
# pylint: disable=wrong-import-position
import urllib3
import re
import numpy as np
import obspy
from obspy import read_inventory, read_events, UTCDateTime as UTC
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy.core import Stream, Trace, UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
from rf import iter_event_data
from tqdm import tqdm
from seismic.misc import setup_logger
import click

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.stream_processing import zne_order
from seismic.stream_io import safe_iter_event_data, write_h5_event_stream
import obspy.core.util.version
from obspy.core.inventory import Inventory
from obspy.taup import TauPyModel

from PhasePApy.phasepapy.phasepicker import aicdpicker
from seismic.pick_harvester.utils import Event, Origin, Magnitude
from seismic.pick_harvester.pick import extract_p, extract_s
from seismic.stream_processing import zerophase_resample

from collections import defaultdict

logging.basicConfig()

# pylint: disable=invalid-name, logging-format-interpolation

SW_MAX_DEPTH = 150  # km


def get_events(lonlat, starttime, endtime, cat_file, distance_range, magnitude_range):
    """Load event catalog (if available) or create event catalog from FDSN server.

    :param lonlat: (Longitude, latitude) of reference location for finding events
    :type lonlat: tuple(float, float)
    :param starttime: Start time of period in which to query events
    :type starttime: obspy.UTCDateTime or str in UTC datetime format
    :param endtime: End time of period in which to query events
    :type endtime: obspy.UTCDateTime or str in UTC datetime format
    :param cat_file: File containing event catalog, or file name in which to store event catalog
    :type cat_file: str or Path
    :param distance_range: Range of distances over which to query seismic events
    :type distance_range: tuple(float, float)
    :param magnitude_range: Range of event magnitudes over which to query seismic events.
    :type magnitude_range: tuple(float, float)
    :param early_exit: If True, exit as soon as new catalog has been generated, defaults to True
    :type early_exit: bool, optional
    :return: Event catalog
    :rtype: obspy.core.event.catalog.Catalog
    """
    log = setup_logger(__name__)

    min_magnitude = magnitude_range[0]
    max_magnitude = magnitude_range[1]
    client = Client('ISC')
    kwargs = {'starttime': starttime, 'endtime': endtime,
              'latitude': lonlat[1], 'longitude': lonlat[0],
              'minradius': distance_range[0], 'maxradius': distance_range[1],
              'minmagnitude': min_magnitude, 'maxmagnitude': max_magnitude}

    log.info("Following parameters will be used for earthquake event query:\n{}".format(kwargs))
    catalog = client.get_events(**kwargs)
    log.info("Catalog loaded from FDSN server")

    # Filter catalog before saving
    catalog = _filter_catalog_events(catalog)

    log.info("Creating catalog file: {}".format(cat_file))
    catalog.write(cat_file, 'QUAKEML')

    return catalog
# end func

def _filter_catalog_events(catalog):
    """Filter catalog with fixed filter criteria.

    :param catalog: Seismic event catalog
    :type catalog: obspy.core.event.catalog.Catalog
    :return: Filtered event catalog
    :rtype: obspy.core.event.catalog.Catalog
    """
    log = logging.getLogger(__name__)

    def _earthquake_event_filter(event):
        return event.get('event_type') == 'earthquake'

    # Type filter
    accepted_events = [e for e in catalog if _earthquake_event_filter(e)]
    catalog = obspy.core.event.catalog.Catalog(accepted_events)

    # Filter out events with missing magnitude or depth
    n_before = len(catalog)
    catalog = catalog.filter("magnitude > 0.0", "depth > 0.0")
    n_after = len(catalog)
    if n_after < n_before:
        log.info("Removed {} events from catalog with invalid magnitude or depth values".format(n_before - n_after))

    # Filter for standard error on travel time residuals
    n_before = len(catalog)
    catalog = catalog.filter("standard_error <= 5.0")
    n_after = len(catalog)
    if n_after < n_before:
        log.info("Removed {} events from catalog with high travel time residuals".format(n_before - n_after))

    return catalog
# end func


def asdf_get_waveforms(asdf_dataset, network, station, location, channel, starttime,
                       endtime):
    """Custom waveform getter function to retrieve waveforms from FederatedASDFDataSet.

    :param asdf_dataset: Instance of FederatedASDFDataSet to query
    :type asdf_dataset: seismic.ASDFdatabase.FederatedASDFDataSet
    :param network: Network code
    :type network: str
    :param station: Station code
    :type station: str
    :param location: Location code
    :type location: str
    :param channel: Channel code
    :type channel: str
    :param starttime: Start time of the waveform query
    :type starttime: str in UTC datetime format
    :param endtime: End time of the waveform query
    :type endtime: str in UTC datetime format
    :return: Stream containing channel traces
    :rtype: obspy.Stream of obspy.Traces
    """
    st = Stream()
    matching_stations = asdf_dataset.get_stations(starttime, endtime, network=network, station=station,
                                                  location=location)
    if matching_stations:
        channel = channel.replace('?', '.')  # replace greedy matching by single-character matching
        ch_matcher = re.compile(channel)
        for net, sta, loc, cha, _, _, _ in matching_stations:
            if ch_matcher.match(cha):
                st += asdf_dataset.get_waveforms(net, sta, loc, cha, starttime, endtime,
                                                 trace_count_threshold=50)
            # end if
        # end for
    # end if
    if st:
        try:
            st = Stream([tr for tr in st if tr.stats.asdf.tag == 'raw_recording'])
        except AttributeError:
            log = logging.getLogger(__name__)
            log.error("ASDF tag not found in Trace stats")
        # end try
    # end if
    return st


# end func

def timestamp_filename(fname, t0, t1):
    """Append pair of timestamps (start and end time) to file name in format that is
       compatible with filesystem file naming.

    :param fname: File name
    :type fname: str or path
    :param t0: first timestamp
    :type t0: obspy.UTCDateTime
    :param t1: second timestamp
    :type t1: obspy.UTCDateTime
    """
    t0_str = t0.strftime("%Y%m%dT%H%M%S")
    t1_str = t1.strftime("%Y%m%dT%H%M%S")
    bname, ext = os.path.splitext(fname)
    bname += ("_" + t0_str + "-" + t1_str)
    return bname + ext


# end func


def is_url(resource_path):
    """Convenience function to check if a given resource path is a valid URL

    :param resource_path: Path to test for URL-ness
    :type resource_path: str
    :return: True if input is a valid URL, False otherwise
    :rtype: bool
    """
    str_parsed = urllib3.util.url.parse_url(resource_path)
    return str_parsed.scheme and str_parsed.netloc


# end func

def trim_inventory(inventory, network_list, station_list):
    """
    Function to trim inventory with a given list of networks and stations.
    Note that duplicate station-names across different networks are not
    currently handled.

    :param inventory: obspy inventory
    :param network_list: a space-separated list of networks
    :param stations_list: a space-separated list of stations
    """

    log = logging.getLogger(__name__)

    if (network_list == '*'):
        network_list = []
    else:
        network_list = re.findall('\S+', network_list)
        assert len(network_list), 'Invalid network list. Aborting..'
    # end if

    if (station_list == '*'):
        station_list = []
    else:
        station_list = re.findall('\S+', station_list)
        assert len(station_list), 'Invalid station list. Aborting..'
    # end if

    if (len(network_list)):
        subset_inv = Inventory(networks=[], source=obspy.core.util.version.read_release_version())
        for net in network_list:
            subset_inv += inventory.select(network=net)
        # end for
        inventory = subset_inv
    # end if

    if (len(station_list)):
        subset_inv = Inventory(networks=[], source=obspy.core.util.version.read_release_version())
        for sta in station_list:
            subset_inv += inventory.select(station=sta)
        # end for
        inventory = subset_inv
    # end if

    net_codes = set()
    sta_codes = set()
    for net in inventory.networks:
        net_codes.add(net.code)
        for sta in net.stations:
            sta_codes.add(sta.code)
        # end for
    # end for

    if (len(sta_codes) == 0):
        log.error('Inventory is empty! Aborting..')
        exit(0)
    # end if

    log.info('Using %d networks: (%s)' % (len(net_codes), ', '.join(net_codes)))
    log.info('and %d stations: (%s)' % (len(sta_codes), ', '.join(sta_codes)))

    return inventory
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('data-source',
                type=click.Path(exists=True))
@click.option('--inventory-file', type=click.Path(exists=True, dir_okay=False), required=False, default=None,
              help=r'Optional path to input inventory file corresponding to waveform source provided through, '
                   r'--waveform-database. Note that this parameter is required only when the waveform source is '
                   r'not a definition file for a FederatedASDFDataSet, in which case, the relevant inventory '
                   r'is extracted internally.')
@click.option('--network-list', default='*', help='A space-separated list of networks (within quotes) to process.',
              type=str,
              show_default=True)
@click.option('--station-list', default='*', help='A space-separated list of stations (within quotes) to process.',
              type=str,
              show_default=True)
@click.option('--event-catalog-file', type=click.Path(dir_okay=False, writable=True), required=True,
              help='Path to output event catalog file, e.g. "catalog_7X_for_rf.xml".')
@click.option('--start-time', type=str, default='', show_default=True,
              help='Start datetime in ISO 8601 format, e.g. "2009-06-16T03:42:00". '
                   'If empty, will be inferred from the inventory file.')
@click.option('--end-time', type=str, default='', show_default=True,
              help='End datetime in ISO 8601 format, e.g. "2011-04-01T23:18:49". '
                   'If empty, will be inferred from the inventory file.')
@click.option('--event-distance-range', type=(float, float), default=(0, 180.0), show_default=True,
              help='Range of teleseismic distances (in degrees) to download events for')
@click.option('--magnitude-range', type=(float, float), default=(5.5, 10.0), show_default=True,
              help='Range of seismic event magnitudes to sample from the event catalog for P/S arrivals.')
def main(data_source, inventory_file, network_list, station_list, event_catalog_file,
         start_time, end_time, event_distance_range, magnitude_range):
    """
    DATA_SOURCE: Text file containing paths to ASDF files.
    """

    log = setup_logger('__func__')

    asdf_dataset = FederatedASDFDataSet(data_source)
    inventory = asdf_dataset.get_inventory()

    log.info("Loaded inventory..")

    inventory = trim_inventory(inventory, network_list=network_list, station_list=station_list)

    lonlat = None
    # Compute reference lonlat from the inventory.
    channels = inventory.get_contents()['channels']
    lonlat_coords = []
    for ch in channels:
        coords = inventory.get_coordinates(ch)
        lonlat_coords.append((coords['longitude'], coords['latitude']))
    lonlat_coords = np.array(lonlat_coords)
    lonlat = np.mean(lonlat_coords, axis=0)
    log.info("Inferred reference coordinates {}".format(lonlat))

    # If start and end time not provided, infer from date range of inventory.
    if not start_time:
        start_time = inventory[0].start_date
        for net in inventory:
            start_time = min(start_time, net.start_date)
        log.info("Inferred start time {}".format(start_time))
    # end if
    if not end_time:
        end_time = inventory[0].end_date
        if end_time is None:
            end_time = UTC.now()
        for net in inventory:
            end_time = max(end_time, net.end_date)
        log.info("Inferred end time {}".format(end_time))
    # end if

    start_time = UTC(start_time)
    end_time = UTC(end_time)

    catalog = get_events(lonlat, start_time, end_time, event_catalog_file, event_distance_range,
                         magnitude_range)
# end main

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
# end if