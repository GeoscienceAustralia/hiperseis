#!/usr/bin/env python

import os.path
import numpy as np
import logging
import re

import warnings
warnings.simplefilter("ignore", UserWarning)

from obspy import read_inventory, read_events, UTCDateTime as UTC
from obspy.clients.fdsn import Client
from rf import iter_event_data
from tqdm import tqdm

import click
from obspy.core import Stream

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet


def get_events(lonlat, starttime, endtime, cat_file, distance_range, magnitude_range, early_exit=True):
    """Load (if available) or create event catalog from FDSN server.

    :param lonlat: [description]
    :type lonlat: [type]
    :param starttime: [description]
    :type starttime: [type]
    :param endtime: [description]
    :type endtime: [type]
    :param cat_file: [description]
    :type cat_file: [type]
    :param distance_range: [description]
    :type distance_range: [type]
    :param magnitude_range: [description]
    :type magnitude_range: [type]
    :param early_exit: [description], defaults to True
    :type early_exit: bool, optional
    :return: [description]
    :rtype: [type]
    """
    log = logging.getLogger(__name__)

    # If file needs to be generated, then this function requires internet access.
    if os.path.exists(cat_file):
        # This is a bad model - it will read events from completely different settings just because the file is there!
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

        # Filter catalog before saving
        catalog = _filter_catalog_events(catalog)

        catalog.write(cat_file, 'QUAKEML')

        if early_exit:
            print("Run this process again using qsub")
            exit(0)
        # end if
    # end if

    return catalog


def _filter_catalog_events(catalog):
    """Filter catalog with fixed filter criteria.

    :param catalog: [description]
    :type catalog: [type]
    :return: [description]
    :rtype: [type]
    """

    def _known_event_filter(event):
        return event['event_type_certainty'] == 'known' and event['event_type'] == 'earthquake'

    # types filter
    accepted_events = [e for e in catalog if _known_event_filter(e)]
    filtered_catalog = obspy.core.event.catalog.Catalog(accepted_events)

    # inbuilt filter for standard error on travel time residuals
    filtered_catalog = filtered_catalog.filter("standard_error <= 5.0")

    return filtered_catalog


def custom_get_waveforms(asdf_dataset, network, station, location, channel, starttime,
                         endtime, quality=None, minimumlength=None,
                         longestonly=None, filename=None, attach_response=False,
                         **kwargs):
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


def get_existing_index(rf_trace_datafile):
    import h5py
    try:
        rf_file_readonly = h5py.File(rf_trace_datafile, mode='r')
        idx = rf_file_readonly['waveforms']
        existing_index = dict((id, set(d for d in idx[id])) for id in idx)
    except Exception:
        existing_index = None
    return existing_index

# ---+----------Main---------------------------------

@click.command()
@click.option('--inventory-file', type=click.Path(exists=True, dir_okay=False), required=True,
              help=r'Path to input inventory file corresponding to waveform file, '
              r'e.g. "/g/data/ha3/Passive/_ANU/7X\(2009-2011\)/ASDF/7X\(2009-2011\)_ASDF.xml".')
@click.option('--waveform-database', type=click.Path(exists=True, dir_okay=False), required=True,
              help=r'Path to waveform database definition file for FederatedASDFDataSet from which to extract traces '
                   r'for RF analysis, e.g. "/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt"')
@click.option('--event-catalog-file', type=click.Path(dir_okay=False, writable=True), required=True,
              help='Path to event catalog file, e.g. "catalog_7X_for_rf.xml". '
              'If file already exists, it will be loaded, otherwise it will be created by querying the ISC web '
              'service. Note that for traceability, start and end times will be appended to file name.')
@click.option('--rf-trace-datafile', type=click.Path(dir_okay=False, writable=True), required=True,
              help='Path to output file, e.g. "7X_event_waveforms_for_rf.h5". '
              'Note that for traceability, start and end times will be appended to file name.')
@click.option('--start-time', type=str, default='', show_default=True,
              help='Start datetime in ISO 8601 format, e.g. "2009-06-16T03:42:00". '
                   'If empty, will be inferred from the inventory file.')
@click.option('--end-time', type=str, default='', show_default=True,
              help='End datetime in ISO 8601 format, e.g. "2011-04-01T23:18:49". '
                   'If empty, will be inferred from the inventory file.')
@click.option('--distance-range', type=(float, float), default=(30.0, 90.0), show_default=True,
              help='Range of teleseismic distances (in degrees) to sample relative to the mean lat,lon location')
@click.option('--magnitude-range', type=(float, float), default=(5.0, 7.0), show_default=True,
              help='Range of seismic event magnitudes to sample from the event catalog.')
def main(inventory_file, waveform_database, event_catalog_file, rf_trace_datafile, start_time, end_time,
         distance_range, magnitude_range):

    assert not os.path.exists(rf_trace_datafile), \
        "Won't delete existing file {}, remove manually.".format(rf_trace_datafile)

    log = logging.getLogger(__name__)

    min_dist_deg = distance_range[0]
    max_dist_deg = distance_range[1]
    min_mag = magnitude_range[0]
    max_mag = magnitude_range[1]

    inventory = read_inventory(inventory_file)

    # Compute reference lonlat from the inventory.
    channels = inventory.get_contents()['channels']
    lonlat_coords = []
    for ch in channels:
        coords = inventory.get_coordinates(ch)
        lonlat_coords.append((coords['longitude'], coords['latitude']))
    lonlat_coords = np.array(lonlat_coords)
    lonlat = np.mean(lonlat_coords, axis=0)

    # If start and end time not provided, infer from date range of inventory.
    if not start_time:
        start_time = inventory[0].start_date
        for net in inventory:
            start_time = min(start_time, net.start_date)
    if not end_time:
        end_time = inventory[0].end_date
        for net in inventory:
            end_time = max(start_time, net.end_date)

    start_time = UTC(start_time)
    end_time = UTC(end_time)
    event_catalog_file = timestamp_filename(event_catalog_file, start_time, end_time)
    rf_trace_datafile = timestamp_filename(rf_trace_datafile, start_time, end_time)

    exit_after_catalog = False
    catalog = get_events(lonlat, start_time, end_time, event_catalog_file, (min_dist_deg, max_dist_deg),
                         (min_mag, max_mag), exit_after_catalog)
    # Filter out events with missing magnitude or depth
    n_before = len(catalog)
    catalog = catalog.filter("magnitude >= 0.0", "depth >= 0.0")
    n_after = len(catalog)
    if n_after < n_before:
        log.info("Removed {} events from catalog with invalid magnitude or depth values".format(n_before - n_after))

    # TODO: This can probably be sped up a lot by splitting event catalog across N processors

    # Form closure to allow waveform source file to be derived from a setting (or command line input)
    asdf_dataset = FederatedASDFDataSet(waveform_database, logger=log)
    def closure_get_waveforms(network, station, location, channel, starttime, endtime):
        return custom_get_waveforms(asdf_dataset, network, station, location, channel, starttime, endtime)

    existing_index = get_existing_index(rf_trace_datafile)

    with tqdm(smoothing=0) as pbar:
        for s in iter_event_data(catalog, inventory, closure_get_waveforms, pbar=pbar):
            # Write traces to output file in append mode so that arbitrarily large file
            # can be processed. If the file already exists, then existing streams will
            # be overwritten rather than duplicated.
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
                        tr.write(rf_trace_datafile, 'H5', mode='a', override='dont')
                else:
                    tr.write(rf_trace_datafile, 'H5', mode='a')
            # end for
        # end for
    # end with


if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
