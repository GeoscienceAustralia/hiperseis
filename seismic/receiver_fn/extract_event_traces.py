#!/usr/bin/env python
"""Use waveform database and station inventory to extract raw traces for all seismic events within a given
magnitude and time range.
"""

import os.path
import logging
import re

import warnings
warnings.simplefilter("ignore", UserWarning)
import urllib3

import numpy as np
import obspy
from obspy import read_inventory, read_events, UTCDateTime as UTC
from obspy.clients.fdsn import Client
from obspy.core import Stream
from rf import iter_event_data
from tqdm import tqdm
import click

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.receiver_fn.rf_util import zne_order

logging.basicConfig()

# pylint: disable=invalid-name, logging-format-interpolation


def get_events(lonlat, starttime, endtime, cat_file, distance_range, magnitude_range, early_exit=True):
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
    log = logging.getLogger(__name__)

    # If file needs to be generated, then this function requires internet access.
    if os.path.exists(cat_file):
        # FIXME: This is a bad design - it will read events from a file generated using completely
        # different settings just because the file is there!
        log.warning("Loading catalog from file {} irrespective of command line options!!!".format(cat_file))
        log.info("Using catalog file: {}".format(cat_file))
        catalog = read_events(cat_file)
    else:
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

        log.info("Creating catalog file: {}".format(cat_file))
        catalog.write(cat_file, 'QUAKEML')

        if early_exit:
            print("Run this process again using qsub")
            exit(0)
        # end if
    # end if

    # Filter catalog before saving
    catalog = _filter_catalog_events(catalog)

    return catalog


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


def custom_get_waveforms(asdf_dataset, network, station, location, channel, starttime,
                         endtime, **kwargs):
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
    """
    st = Stream()
    matching_stations = asdf_dataset.get_stations(starttime, endtime, network=network, station=station,
                                                  location=location)
    if matching_stations:
        ch_matcher = re.compile(channel)
        for net, sta, loc, cha, _, _ in matching_stations:
            if ch_matcher.match(cha):
                st += asdf_dataset.get_waveforms(net, sta, loc, cha, starttime, endtime)
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


def _get_existing_index(rf_trace_datafile):
    import h5py
    try:
        rf_file_readonly = h5py.File(rf_trace_datafile, mode='r')
        idx = rf_file_readonly['waveforms']
        existing_index = dict((id, set(d for d in idx[id])) for id in idx)
    except Exception:
        existing_index = None
    return existing_index


def is_url(resource_path):
    """Convenience function to check if a given resource path is a valid URL

    :param resource_path: Path to test for URL-ness
    :type resource_path: str
    :return: True if input is a valid URL, False otherwise
    :rtype: bool
    """
    str_parsed = urllib3.util.url.parse_url(resource_path)
    return str_parsed.scheme and str_parsed.netloc


# ---+----------Main---------------------------------

@click.command()
@click.option('--inventory-file', type=click.Path(exists=True, dir_okay=False), required=True,
              help=r'Path to input inventory file corresponding to waveform file, '
              r'e.g. "/g/data/ha3/Passive/_ANU/7X\(2009-2011\)/ASDF/7X\(2009-2011\)_ASDF.xml".')
@click.option('--waveform-database', type=str, required=True,
              help=r'Location of waveform source database from which to extract traces. May be a recognized service '
                   r'provider from obspy.clients.fdsn.header.URL_MAPPINGS (e.g. "ISC"), an actual URL '
                   r'(e.g. "http://auspass.edu.au") or a file path. If detected as a URL, the obspy client '
                   r'get_waveform function will be used to retrieve waveforms from web service. Otherwise, if detected '
                   r'as a valid file path, then it must be the path to a definition file for a FederatedASDFDataSet, '
                   r'e.g. "/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt".')
@click.option('--event-catalog-file', type=click.Path(dir_okay=False, writable=True), required=True,
              help='Path to event catalog file, e.g. "catalog_7X_for_rf.xml". '
              'If file already exists, it will be loaded, otherwise it will be created by querying the ISC web '
              'service. Note that for traceability, start and end times will be appended to file name.')
@click.option('--rf-trace-datafile', type=click.Path(dir_okay=False, writable=True), required=True,
              help='Path to output file, e.g. "7X_event_waveforms_for_rf.h5". '
              'Note that for traceability, start and end times will be appended to file name. '
              'If file already exists, pre-existing records will be preserved, not overwritten.')
@click.option('--start-time', type=str, default='', show_default=True,
              help='Start datetime in ISO 8601 format, e.g. "2009-06-16T03:42:00". '
                   'If empty, will be inferred from the inventory file.')
@click.option('--end-time', type=str, default='', show_default=True,
              help='End datetime in ISO 8601 format, e.g. "2011-04-01T23:18:49". '
                   'If empty, will be inferred from the inventory file.')
@click.option('--taup-model', type=str, default='iasp91', show_default=True,
              help='Theoretical tau-p Earth model to use for Trace stats computation. Other possibilities, '
                   'such as ak135, are documented here: https://docs.obspy.org/packages/obspy.taup.html')
@click.option('--distance-range', type=(float, float), default=(30.0, 90.0), show_default=True,
              help='Range of teleseismic distances (in degrees) to sample relative to the mean lat,lon location')
@click.option('--magnitude-range', type=(float, float), default=(5.5, 7.0), show_default=True,
              help='Range of seismic event magnitudes to sample from the event catalog.')
def main(inventory_file, waveform_database, event_catalog_file, rf_trace_datafile, start_time, end_time, taup_model,
         distance_range, magnitude_range):

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    waveform_db_is_web = is_url(waveform_database) or waveform_database in obspy.clients.fdsn.header.URL_MAPPINGS
    if not waveform_db_is_web:
        assert os.path.exists(waveform_database), "Cannot find waveform database file {}".format(waveform_database)
    log.info("Using waveform data source: {}".format(waveform_database))

    assert not os.path.exists(rf_trace_datafile), \
        "Won't delete existing file {}, remove manually.".format(rf_trace_datafile)

    min_dist_deg = distance_range[0]
    max_dist_deg = distance_range[1]
    min_mag = magnitude_range[0]
    max_mag = magnitude_range[1]

    inventory = read_inventory(inventory_file)
    log.info("Loaded inventory {}".format(inventory_file))

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
    event_catalog_file = timestamp_filename(event_catalog_file, start_time, end_time)
    rf_trace_datafile = timestamp_filename(rf_trace_datafile, start_time, end_time)
    log.info("Traces will be written to: {}".format(rf_trace_datafile))

    exit_after_catalog = False
    catalog = get_events(lonlat, start_time, end_time, event_catalog_file, (min_dist_deg, max_dist_deg),
                         (min_mag, max_mag), exit_after_catalog)

    if waveform_db_is_web:
        existing_index = None
        log.info("Use fresh query results from web")
        client = Client(waveform_database)
        waveform_getter = client.get_waveforms
    else:
        # Form closure to allow waveform source file to be derived from a setting (or command line input)
        asdf_dataset = FederatedASDFDataSet(waveform_database, logger=log)
        def closure_get_waveforms(network, station, location, channel, starttime, endtime):
            return custom_get_waveforms(asdf_dataset, network, station, location, channel, starttime, endtime)
        existing_index = _get_existing_index(rf_trace_datafile)
        if existing_index is not None:
            log.warning("Resuming extraction using existing index from file {}".format(rf_trace_datafile))
        waveform_getter = closure_get_waveforms
    # end if

    with tqdm(smoothing=0) as pbar:
        stream_count = 0
        for s in iter_event_data(catalog, inventory, waveform_getter, tt_model=taup_model, pbar=pbar):
            # Write traces to output file in append mode so that arbitrarily large file
            # can be processed. If the file already exists, then existing streams will
            # be overwritten rather than duplicated.
            # Check first if rotation for unaligned *H1, *H2 channels to *HN, *HE is required.
            if s.select(component='1') and s.select(component='2'):
                s.rotate('->ZNE', inventory=inventory)
            # end if
            # Order the traces in ZNE ordering. This is required so that normalization
            # can be specified in terms of an integer index, i.e. the default of 0 in rf
            # library will normalize against the Z component.
            s.traces = sorted(s.traces, key=zne_order)
            # Assert the ordering of traces in the stream is ZNE.
            assert s.traces[0].stats.channel[-1] == 'Z'
            assert s.traces[1].stats.channel[-1] == 'N'
            assert s.traces[2].stats.channel[-1] == 'E'
            # Loop over ZNE traces
            for tr in s:
                grp_id = '.'.join(tr.id.split('.')[0:3])
                event_time = str(tr.meta.event_time)[0:19]
                pbar.set_description("{} -- {}".format(grp_id, event_time))
                if existing_index is not None:
                    # Skip records that already exist in the file to speed up generation
                    if grp_id in existing_index and event_time in existing_index[grp_id]:
                        pbar.write("Skipping {} -- {} already exists in output file".format(grp_id, event_time))
                        continue
                    else:
                        # Use don't override mode just in case our hand-crafted index is faulty
                        stream_count += 1
                        tr.write(rf_trace_datafile, 'H5', mode='a', override='dont')
                else:
                    stream_count += 1
                    tr.write(rf_trace_datafile, 'H5', mode='a')
            # end for
        # end for
        if stream_count == 0:
            log.warning("No traces found!")
        else:
            log.info("Wrote {} new stream to output file".format(stream_count))
        # end if
    # end with


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
