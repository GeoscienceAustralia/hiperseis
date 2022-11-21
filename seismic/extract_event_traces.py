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
from obspy.core import Stream
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees
from rf import iter_event_data
from tqdm import tqdm
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

SW_MAX_DEPTH = 150 #km

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
        min_magnitude = magnitude_range[0]
        max_magnitude = magnitude_range[1]
        
        # For HPC systems with no internet access, the catalog file must be pre-generated
        log.warning("Loading catalog from file {} irrespective of command line options!!!".format(cat_file))
        log.info("Using catalog file: {}".format(cat_file))
        catalog = read_events(cat_file)
        
        # While events are downloaded, the ISC magnitude filter does not efectively cull earthquakes outside the requested 
        # magnitude-range. A secondary magnitude filter is therefore added here to cull events outside the magnitude-range.
        catalog = catalog.filter("time > {}".format(str(starttime)), "time < {}".format(str(endtime)),
                                 "magnitude >= {}".format(min_magnitude), "magnitude <= {}".format(max_magnitude))
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
        ch_matcher = re.compile(channel)
        for net, sta, loc, cha, _, _, _ in matching_stations:
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

    if(network_list=='*'):
        network_list = []
    else:
        network_list = re.findall('\S+', network_list)
        assert len(network_list), 'Invalid network list. Aborting..'
    # end if

    if(station_list=='*'):
        station_list = []
    else:
        station_list = re.findall('\S+', station_list)
        assert len(station_list), 'Invalid station list. Aborting..'
    # end if

    if(len(network_list)):
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

    if(len(sta_codes) == 0):
        log.error('Inventory is empty! Aborting..')
        exit(0)
    # end if

    log.info('Using %d networks: (%s)'%(len(net_codes), ', '.join(net_codes)))
    log.info('and %d stations: (%s)'%(len(sta_codes), ', '.join(sta_codes)))

    return inventory
#end func

class Picker():
    counter = 0 # static variable to count number of calls made

    def __init__(self, taup_model_name):
        self._taup_model = TauPyModel(model=taup_model_name)
        self._picker_list_p = []
        self._picker_list_s = []

        sigmalist = np.arange(8, 3, -1)
        for sigma in sigmalist:
            picker_p = aicdpicker.AICDPicker(t_ma=10, nsigma=sigma, t_up=1, nr_len=5,
                                             nr_coeff=2, pol_len=10, pol_coeff=10, uncert_coeff=3)
            picker_s = aicdpicker.AICDPicker(t_ma=15, nsigma=sigma, t_up=1, nr_len=5,
                                             nr_coeff=2, pol_len=10, pol_coeff=10, uncert_coeff=3)

            self._picker_list_p.append(picker_p)
            self._picker_list_s.append(picker_s)
        # end for
    # end func

    def pick(self, ztrace, ntrace, etrace, phase='P'):
        slope_ratio = -1
        arrival_time = UTC(-1)

        # construct a named array for event meta-data, as expected in extract_[p/s]
        event_fields = {'names': ['source', 'event_id', 'origin_ts', 'mag', 'lon', 'lat', 'depth_km'],
                        'formats': ['S10', 'i4', 'f8', 'f4', 'f4', 'f4', 'f4']}
        events = np.array([('', 0, ztrace.stats.event_time.timestamp,
                            ztrace.stats.event_magnitude, ztrace.stats.event_longitude,
                            ztrace.stats.event_latitude, ztrace.stats.event_depth)], dtype=event_fields)

        result = None
        if(phase == 'P'):
            result = extract_p(self._taup_model, self._picker_list_p, events[0], ztrace.stats.station_longitude,
                               ztrace.stats.station_latitude, Stream(ztrace), margin=5)
        elif(phase == 'S'):
            result = extract_s(self._taup_model, self._picker_list_s, events[0], ztrace.stats.station_longitude,
                               ztrace.stats.station_latitude, Stream(ntrace), Stream(etrace),
                               ntrace.stats.back_azimuth, margin=10)
        else:
            assert 0, 'Unknown phase: {}. Must be "P" or "S"'.format(phase)
        # end if

        if(result):
            picklist, residuallist, snrlist, _, _ = result
            best_pick_idx = np.argmax(snrlist[:, -1]) # hightest slope-ratio quality-estimate

            arrival_time = picklist[best_pick_idx]
            slope_ratio = snrlist[best_pick_idx, -1]
        # end if

        ztrace.stats.update({'arrival_time': arrival_time, 'slope_ratio': slope_ratio})
        ntrace.stats.update({'arrival_time': arrival_time, 'slope_ratio': slope_ratio})
        etrace.stats.update({'arrival_time': arrival_time, 'slope_ratio': slope_ratio})
    # end func
# end class

def sw_catalog(catalog, min_mag):
    # Trim catalog for surface waves, removing proximal events, as done in function catclean in:
    # https://github.com/jbrussell/DLOPy_v1.0/blob/master/pysave/locfuns.py
    def close(x1, x2):
        if (np.fabs(x1-x2) < 0.8):
            return True
        else:
            return False
        # end if
    #end func

    stime = np.array([e.preferred_origin().time.timestamp for e in catalog])
    lon = np.array([e.preferred_origin().longitude for e in catalog])
    lat = np.array([e.preferred_origin().latitude for e in catalog])
    depth = np.array([e.preferred_origin().depth for e in catalog])
    mag = np.array([float(e.magnitudes[0].mag) for e in catalog])

    repeated = []
    for i in np.arange((len(stime))):
        for j in np.arange((len(stime))):
            if(i==j): continue

            if (np.fabs(stime[j]-stime[i]) < 60*15 and
                    close(lat[j], lat[i]) and
                    close(lon[j], lon[i]) and
                    np.fabs(mag[j] - mag[i]) < 0.3):
                repeated.append(j)
            # end if
        # end for
    # end for
    repeated = set(repeated)

    out_cat = Catalog()
    for i, e in enumerate(catalog):
        if(e.magnitudes[0].mag < min_mag): continue
        if(e.preferred_origin().depth/1e3 > SW_MAX_DEPTH): continue

        if(i not in repeated):
            out_cat.append(e)
        # end if
    # end for

    return out_cat
# end func

def extract_data(catalog, inventory, waveform_getter, event_trace_datafile,
                 resample_hz, tt_model='iasp91', wave='P', pad=10,
                 dry_run=True):

    assert wave in ['P', 'S', 'SW'], 'Only P, S and SW (surface wave) is supported. Aborting..'

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)

    # descriptions
    descs = {'P': 'P-wave', 'S': 'S-wave', 'SW': 'Surface-wave'}

    # initialize event time-window dict
    request_window = defaultdict(tuple) # seconds
    request_window['P'] = (-70, 150)
    request_window['S'] = (-100, 150)
    request_window['SW'] = (-70, 4*60*60)

    # initialize phase-map dict
    phase_map = defaultdict(str) # seconds
    phase_map['P'] = 'P'
    phase_map['S'] = 'S'
    # for surface-waves we use the default phase (P), but internally, safe_iter_event_data
    # extracts data around event origin time
    phase_map['SW'] = 'P'

    # initialize event distance-range dict
    distance_range = defaultdict(tuple) # arc degrees
    distance_range['P'] = (30, 90)
    distance_range['S'] = (30, 90)
    distance_range['SW'] = (5, 175)

    # initialize dict that indicates whether rfstats should be generated
    rfstats_map = defaultdict(bool) # seconds
    rfstats_map['P'] = True
    rfstats_map['S'] = True
    rfstats_map['SW'] = False # for surface-waves we don't need rfstats

    # instantiate arrival-picker
    picker = Picker(taup_model_name=tt_model)

    # initialize trace-data organization scheme
    # P-waveforms are stored under root group 'waveforms' for backward compatibility
    tf = '.datetime:%Y-%m-%dT%H:%M:%S'
    h5_index = 'waveforms/{wave_type}/{network}.{station}.{location}/{event_time%s}/' % tf + \
                         '{channel}_{starttime%s}_{endtime%s}' % (tf, tf)

    # Initialize MPI
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    #################################################
    # data extraction is parallelized over stations
    #################################################
    nsl_dict = None
    if(rank==0):
        nsl_dict = []
        for i in np.arange(nproc): nsl_dict.append(defaultdict(list)) # net.sta.loc -> cha
        
        temp_dict = defaultdict(list)
        for item in inventory.get_contents()['channels']:
            tokens = item.split('.')
            temp_dict['.'.join(tokens[:3])] = tokens[-1]
        # end for

        njob = len(temp_dict) # total number of stations
        # Add made up entries to ensure MPI-barrier calls are balanced across all
        # processors
        nbogus = np.int(np.ceil(njob/float(nproc)))*nproc - njob
        for i in np.arange(nbogus): temp_dict['%i.%i.%i'%(i, i, i)] = '-1'

        cproc = 0
        for k, v in temp_dict.items():
            nsl_dict[cproc][k] = v
            cproc = (cproc + 1)%nproc
        # end for
    # end if

    nsl_dict = comm.scatter(nsl_dict, root=0)

    for nsl, cha in nsl_dict.items():
        if(cha == '-1'):
            # Nothing to do for made up entries, which exist for the sole purpose of balancing
            # MPI-barrier calls across all processors
            for irank in np.arange(nproc):
                comm.Barrier()
            # end for
        else:
            net, sta, loc = nsl.split('.')

            curr_inv = inventory.select(network=net, station=sta, location=loc)

            coord = curr_inv.get_coordinates(nsl + '.' + cha)
            sta_lon, sta_lat = coord['longitude'], coord['latitude']

            if(dry_run): continue

            stream_count = 0
            sta_stream = Stream()

            pbar=tqdm(desc=descs[wave], total=len(catalog) * len(curr_inv))
            for s in safe_iter_event_data(catalog, curr_inv, waveform_getter,
                                          use_rfstats=rfstats_map[wave],
                                          phase=phase_map[wave],
                                          tt_model=tt_model, pbar=None,#pbar,
                                          request_window=request_window[wave],
                                          pad=pad,
                                          dist_range=distance_range[wave]):
                # Write traces to output file in append mode so that arbitrarily large file
                # can be processed. If the file already exists, then existing streams will
                # be overwritten rather than duplicated.
                # Check first if rotation for unaligned *H1, *H2 channels to *HN, *HE is required.
                if not s:
                    continue
                # end if
                if s.select(component='1') and s.select(component='2'):
                    try:
                        s.rotate('->ZNE', inventory=inventory)
                    except Exception as e:
                        log.error('Unable to rotate to ZNE with error:\n{}'.format(str(e)))
                        continue
                    # end try
                # end if
                # Order the traces in ZNE ordering. This is required so that normalization
                # can be specified in terms of an integer index, i.e. the default of 0 in rf
                # library will normalize against the Z component.
                s.traces = sorted(s.traces, key=zne_order)
                # Assert the ordering of traces in the stream is ZNE.
                assert s[0].stats.channel[-1] == 'Z'
                assert s[1].stats.channel[-1] == 'N'
                assert s[2].stats.channel[-1] == 'E'

                # don't pick for surface-waves
                if(rfstats_map[wave]): picker.pick(s[0], s[1], s[2], phase=phase_map[wave])

                # Iterator returns rf.RFStream. Write traces from obspy.Stream to decouple from RFStream.
                grp_id = '.'.join(s.traces[0].id.split('.')[0:3])
                event_time = str(s.traces[0].meta.event_time)[0:19]

                out_stream = obspy.Stream([tr for tr in s])
                assert out_stream[0].stats.channel[-1] == 'Z'
                assert out_stream[1].stats.channel[-1] == 'N'
                assert out_stream[2].stats.channel[-1] == 'E'

                # resample after lowpass @ resample_rate / 2 Hz
                for tr in out_stream:
                    zerophase_resample(tr, resample_hz)
                   
                    tr.stats.update({'wave_type':wave})
                # end for

                sta_stream += out_stream
                stream_count += 1

                pbar.set_description("[{}] {} | {}".format(descs[wave], grp_id, event_time))
                pbar.update()
            # end for
            pbar.close()

            for irank in np.arange(nproc):
                if(irank == rank):
                    if(len(sta_stream)):
                        write_h5_event_stream(event_trace_datafile, sta_stream, index=h5_index, mode='a')
                    # end if
                # end if
                comm.Barrier()
            # end for

            if stream_count == 0:
                log.warning("{}: No traces found!".format(nsl))
            else:
                log.info("{}: Wrote {} {} streams to output file".format(nsl, stream_count, descs[wave]))
            # end if
        # end if
    # end for
# end func

# ---+----------Main---------------------------------

@click.command()
@click.option('--inventory-file', type=click.Path(exists=True, dir_okay=False), required=False, default=None,
              help=r'Optional path to input inventory file corresponding to waveform source provided through, '
                   r'--waveform-database. Note that this parameter is required only when the waveform source is '
                   r'not a definition file for a FederatedASDFDataSet, in which case, the relevant inventory '
                   r'is extracted internally.')
@click.option('--network-list', default='*', help='A space-separated list of networks (within quotes) to process.', type=str,
              show_default=True)
@click.option('--station-list', default='*', help='A space-separated list of stations (within quotes) to process.', type=str,
              show_default=True)
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
@click.option('--event-trace-datafile', type=click.Path(dir_okay=False, writable=True), required=True,
              help='Path to output file, e.g. "7X_event_waveforms.h5". '
                   'Note that for traceability, start and end datetimes will be appended to file name.')
@click.option('--start-time', type=str, default='', show_default=True,
              help='Start datetime in ISO 8601 format, e.g. "2009-06-16T03:42:00". '
                   'If empty, will be inferred from the inventory file.')
@click.option('--end-time', type=str, default='', show_default=True,
              help='End datetime in ISO 8601 format, e.g. "2011-04-01T23:18:49". '
                   'If empty, will be inferred from the inventory file.')
@click.option('--taup-model', type=str, default='iasp91', show_default=True,
              help='Theoretical tau-p Earth model to use for Trace stats computation. Other possibilities, '
                   'such as ak135, are documented here: https://docs.obspy.org/packages/obspy.taup.html')
@click.option('--distance-range', type=(float, float), default=(0, 180.0), show_default=True,
              help='Range of teleseismic distances (in degrees) to sample relative to the mean lat,lon location')
@click.option('--magnitude-range', type=(float, float), default=(5.5, 10.0), show_default=True,
              help='Range of seismic event magnitudes to sample from the event catalog for P/S arrivals.')
@click.option('--sw-magnitude-range', type=(float, float), default=(6.0, 10.0), show_default=True,
              help='Range of seismic event magnitudes to sample from the event catalog for surface waves.')
@click.option('--catalog-only', is_flag=True, default=False, show_default=True,
              help='If set, only generate catalog file and exit. Used for preparing '
                   'input file on HPC systems with no internet access.')
@click.option('--resample-hz', type=float, default=10, show_default=True,
              help='Resampling frequency (default 10 Hz) for output P/S traces')
@click.option('--sw-resample-hz', type=float, default=2, show_default=True,
              help='Resampling frequency (default 2 Hz) for surface waves')
@click.option('--p-data', is_flag=True, default=False, show_default=True,
              help='Extracts waveform data around P-arrival')
@click.option('--s-data', is_flag=True, default=False, show_default=True,
              help='Extracts waveform data around S-arrival')
@click.option('--sw-data', is_flag=True, default=False, show_default=True,
              help='Extracts waveform data around surface-wave arrival')
@click.option('--dry-run', is_flag=True, default=False, show_default=True,
              help='Reports events available to each station, by wave-type and exits without outputting any data. '
                   'Has no effect on --catalog-only mode.')
def main(inventory_file, network_list, station_list, waveform_database, event_catalog_file, event_trace_datafile,
         start_time, end_time, taup_model, distance_range, magnitude_range, sw_magnitude_range, catalog_only,
         resample_hz, sw_resample_hz, p_data, s_data, sw_data, dry_run):
    
    # Initialize MPI
    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    log = logging.getLogger('extract_event_traces')
    log.setLevel(logging.INFO)

    # sanity check
    owave_types = defaultdict(bool)
    if(not(p_data or s_data or sw_data) and not catalog_only):
        assert 0, 'At least one from [--p-data, --s-data, --sw-data] must be specified. Aborting'
    else:
        owave_types['P'] = p_data
        owave_types['S'] = s_data
        owave_types['SW'] = sw_data
    # end if

    inventory = None
    asdf_dataset = None
    pad = 10 # nominal padding for waveforms in seconds
    waveform_db_is_web = is_url(waveform_database) or waveform_database in obspy.clients.fdsn.header.URL_MAPPINGS
    if not waveform_db_is_web:
        assert os.path.exists(waveform_database), "Cannot find waveform database file {}".format(waveform_database)
        asdf_dataset = FederatedASDFDataSet(waveform_database)
        inventory = asdf_dataset.get_inventory()

        #################################################
        # Check if GPS clock-corrections are being applied
        # A large padding is used to allow for time-shifts
        # from clock-correction
        #################################################
        if (asdf_dataset.corrections_enabled()): pad = 3600
    else:
        assert inventory_file, 'Must provide inventory file if using a URL or an obspy client as waveform source'
        inventory = read_inventory(inventory_file)
        log.info("Loaded inventory {}".format(inventory_file))
    # end if

    inventory = trim_inventory(inventory, network_list=network_list, station_list=station_list)

    if(rank == 0): log.info("Using waveform data source: {}".format(waveform_database))

    min_dist_deg = distance_range[0]
    max_dist_deg = distance_range[1]
    min_mag = magnitude_range[0]
    max_mag = magnitude_range[1]
    
    lonlat = None
    if(rank == 0):
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
        if not os.path.exists(event_catalog_file):
            event_catalog_file = timestamp_filename(event_catalog_file, start_time, end_time)
        event_trace_datafile = timestamp_filename(event_trace_datafile, start_time, end_time)
        assert not os.path.exists(event_trace_datafile), \
            "Output file {} already exists, please remove!".format(event_trace_datafile)
        log.info("Traces will be written to: {}".format(event_trace_datafile))
    # end if
    lonlat = comm.bcast(lonlat, root=0)
    start_time = comm.bcast(start_time, root=0)
    end_time = comm.bcast(end_time, root=0)
    event_trace_datafile = comm.bcast(event_trace_datafile, root=0)

    exit_after_catalog = catalog_only
    catalog = get_events(lonlat, start_time, end_time, event_catalog_file, (min_dist_deg, max_dist_deg),
                         (min_mag, max_mag), exit_after_catalog)

    if waveform_db_is_web:
        log.info("Use fresh query results from web")
        client = Client(waveform_database)
        waveform_getter = client.get_waveforms
    else:
        # Form closure to allow waveform source file to be derived from a setting (or command line input)
        def closure_get_waveforms(network, station, location, channel, starttime, endtime):
            return asdf_get_waveforms(asdf_dataset, network, station, location, channel, starttime, endtime)
        waveform_getter = closure_get_waveforms
    # end if

    for wave, flag in owave_types.items():
        if(not flag): continue

        curr_catalog = catalog
        curr_resample_hz = resample_hz

        if(wave == 'SW'):
            # remove proximal events for surface-wave output
            curr_catalog = sw_catalog(catalog, sw_magnitude_range[0])
            curr_resample_hz = sw_resample_hz
        # end if

        extract_data(curr_catalog, inventory, waveform_getter, event_trace_datafile,
                     resample_hz=curr_resample_hz, tt_model=taup_model, wave=wave,
                     pad=pad,
                     dry_run=dry_run)
    # end for
    del asdf_dataset
    
    if(rank == 0):
        print("Finishing...")
        print("extract_event_traces SUCCESS!")
    # end if
# end main

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
# end if
