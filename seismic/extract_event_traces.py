#!/usr/bin/env python
"""Use waveform database and station inventory to extract raw traces for all seismic events within a given
magnitude and time range.
"""

import os.path
from mpi4py import MPI
import re
import numpy as np
import obspy
from obspy.core.event import Catalog
from obspy.core import Stream, Trace, UTCDateTime
import click
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.stream_processing import zne_order
from seismic.stream_io import safe_iter_event_data, write_h5_event_stream
import obspy.core.util.version
from obspy.core.inventory import Inventory
from obspy.taup import TauPyModel

from PhasePApy.phasepapy.phasepicker import aicdpicker
from seismic.pick_harvester.pick import extract_p, extract_s
from seismic.stream_processing import zerophase_resample
from seismic.catalogue.gcmt import GCMTCatalog
from collections import defaultdict
from seismic.misc import setup_logger, print_exception
from seismic.misc_p import parallel_abort
from pandas import DataFrame
import warnings
warnings.simplefilter("ignore", UserWarning)

SW_MAX_DEPTH = 150 #km
# descriptions
DESCS = {'P': 'P-wave', 'S': 'S-wave', 'SW': 'Surface-wave'}

def asdf_get_waveforms(ds:FederatedASDFDataSet, network, station, location, channel, starttime,
                       endtime):
    """Custom waveform getter function to retrieve waveforms from FederatedASDFDataSet.

    :param ds: Instance of FederatedASDFDataSet to query
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
    matching_stations = ds.get_stations(starttime, endtime, network=network, station=station,
                                                  location=location)
    if matching_stations:
        channel = channel.replace('?', '.') # replace greedy matching by single-character matching
        ch_matcher = re.compile(channel)
        for net, sta, loc, cha, _, _, _ in matching_stations:
            if ch_matcher.match(cha):
                st += ds.get_waveforms(net, sta, loc, cha, starttime, endtime,
                                                 trace_count_threshold=50)
            # end if
        # end for
    # end if

    return st
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
        arrival_time = UTCDateTime(-1)

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

def extract_data(recording_timespan_getter, waveform_getter,
                 catalog, inventory, event_trace_datafile, log_folder,
                 wave, request_window, time_range, distance_range, magnitude_range,
                 depth_range, min_areal_separation_km, resample_hz, tt_model='iasp91', pad=10,
                 dry_run=True):
    assert wave in ['P', 'S', 'SW'], 'Only P, S and SW (surface wave) is supported. Aborting..'

    # initialize phase-map dict
    phase_map = defaultdict(str) # seconds
    phase_map['P'] = 'P'
    phase_map['S'] = 'S'
    # for surface-waves we use the default phase (P), but internally, safe_iter_event_data
    # extracts data around event origin time
    phase_map['SW'] = 'P'

    # initialize dict that indicates whether rfstats should be generated
    rfstats_map = defaultdict(bool) # seconds
    rfstats_map['P'] = True
    rfstats_map['S'] = True
    rfstats_map['SW'] = False # for surface-waves we don't need rfstats

    # instantiate arrival-picker
    picker = Picker(taup_model_name=tt_model)

    # initialize trace-data organization scheme
    # P-waveforms are stored under root group 'waveforms' for backward compatibility
    tf = '.datetime:%Y-%m-%dT%H:%M:%S.%f'
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
            if(dry_run): continue
            # Nothing to do for made up entries, which exist for the sole purpose of balancing
            # MPI-barrier calls across all processors
            for irank in np.arange(nproc):
                comm.Barrier()
            # end for
        else:
            log_fn = os.path.join(log_folder, '{}.{}.log'.format(nsl, wave))
            log = setup_logger('__func__', log_fn)
            net, sta, loc = nsl.split('.')

            curr_inv = inventory.select(network=net, station=sta, location=loc)

            coord = curr_inv.get_coordinates(nsl + '.' + cha)
            sta_lon, sta_lat = coord['longitude'], coord['latitude']

            # set start- and end-times
            st, et = recording_timespan_getter(network=net, station=sta, location=loc)
            if(time_range[0] is None):
                time_range[0] = st
            else:
                time_range[0] = UTCDateTime(time_range[0])
                if(time_range[0] < st): time_range[0] = st
            # end if
            if(time_range[1] is None):
                time_range[1] = et
            else:
                time_range[1] = UTCDateTime(time_range[1])
                if(time_range[1] > et): time_range[1] = et
            # end if

            # tailor catalog for current station
            log.info(f"""Pruning catalog for:
\tlocation: {[sta_lon, sta_lat]} 
\ttime range: [{time_range[0]} -- {time_range[1]}] 
\tdistance range: [{distance_range[0]} -- {distance_range[1]}] deg 
\tmagnitude range: [{magnitude_range[0]} -- {magnitude_range[1]}] 
\tdepth range: [{depth_range[0]} -- {depth_range[1]}] km
\tminimum areal separation: {min_areal_separation_km} km """)
            curr_cat = catalog.prune(time_range, sta_lon, sta_lat,
                                     distance_range, magnitude_range,
                                     depth_range, min_areal_separation_km).to_obspy_catalog()
            log.info('A total of {} events retained in catalog.\n'.format(len(curr_cat)))

            log.info('Extracting data windows [{} -- {}] s around events..\n'.format(*request_window))

            if(dry_run): continue # nothing more to do for dry-runs

            stream_count = 0
            sta_stream = Stream()
            status = defaultdict(int)
            log.info('Data extraction stats:\n')
            for s in safe_iter_event_data(curr_cat, curr_inv, waveform_getter,
                                          use_rfstats=rfstats_map[wave],
                                          phase=phase_map[wave],
                                          tt_model=tt_model, pbar=None,
                                          request_window=request_window,
                                          pad=pad, status=status, log=log):
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
                resample_failed = False
                for tr in out_stream:
                    tr.detrend()
                    tr.taper(max_percentage=0.05, max_length=5)

                    try:
                        zerophase_resample(tr, resample_hz)
                    except Exception as e:
                        log.warn('Resampling failed for trace: {}, with exception: {}. Moving along..'.\
                                 format(tr.stats, e))
                        resample_failed = True
                        break
                    # end try
                    tr.stats.update({'wave_type':wave})
                # end for
                if(resample_failed): continue

                sta_stream += out_stream
                stream_count += 1
            # end for

            for irank in np.arange(nproc):
                if(irank == rank):
                    if(len(sta_stream)):
                        write_h5_event_stream(event_trace_datafile, sta_stream, index=h5_index, mode='a')
                    else:
                        t = Trace(data=np.array([0]),
                                  header={'network': net, 'station': sta,
                                          'location': loc, 'channel': 'XXX',
                                          'wave_type': wave,
                                          'station_longitude': sta_lon,
                                          'station_latitude': sta_lat,
                                          'event_time': UTCDateTime(0)})
                        write_h5_event_stream(event_trace_datafile, Stream([t]), index=h5_index, mode='a')
                    # end if
                # end if
                comm.Barrier()
            # end for

            log.info('\nSummary:\n', extra={'simple': True})
            for k, v in status.items():
                log.info('{}: good data found for {}/{} events.'.format \
                         (k, v, len(curr_cat)), extra={'simple': True})
            # end for

            warn_str = \
            " No {} traces found for {}! Added a null trace.".format(DESCS[wave], nsl)
            if stream_count == 0: log.warning(warn_str)
        # end if
    # end for
# end func

# ---+----------Main---------------------------------
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'], show_default=True)
@click.command()
@click.argument('data-source',
                type=click.Path(exists=True))
@click.option('--network-list', default='*', help='A space-separated list of networks (within quotes) to process.', type=str,
              show_default=True)
@click.option('--station-list', default='*', help='A space-separated list of stations (within quotes) to process.', type=str,
              show_default=True)
@click.option('--gcmt-catalog-file', type=click.Path(dir_okay=False), required=True,
              help='Path to gcmt catalog file. ')
@click.option('--output-file', type=click.Path(dir_okay=False, writable=True), required=True,
              help='Path to output file, e.g. "7X_event_waveforms.h5".')
@click.option('--log-folder', type=click.Path(dir_okay=True, file_okay=False, writable=True), required=True,
              help='Path to folder in which log files are to be output.')
@click.option('--start-time', type=str, default=None, show_default=True,
              help='Start datetime in ISO 8601 format, e.g. "2009-06-16T03:42:00". '
                   'If empty, will be inferred from the inventory file.')
@click.option('--end-time', type=str, default=None, show_default=True,
              help='End datetime in ISO 8601 format, e.g. "2011-04-01T23:18:49". '
                   'If empty, will be inferred from the inventory file.')
@click.option('--p-data', is_flag=True, default=False, show_default=True,
              help='Extracts waveform data around P-arrival')
@click.option('--s-data', is_flag=True, default=False, show_default=True,
              help='Extracts waveform data around S-arrival')
@click.option('--sw-data', is_flag=True, default=False, show_default=True,
              help='Extracts waveform data around surface-wave arrival')
@click.option('--p-magnitude-range', type=(float, float), default=(5.5, 10.0), show_default=True,
              help='Range of seismic event magnitudes to sample from the event catalog for P arrivals.')
@click.option('--s-magnitude-range', type=(float, float), default=(5.5, 10.0), show_default=True,
              help='Range of seismic event magnitudes to sample from the event catalog for S arrivals.')
@click.option('--sw-magnitude-range', type=(float, float), default=(6.0, 10.0), show_default=True,
              help='Range of seismic event magnitudes to sample from the event catalog for surface waves.')
@click.option('--p-data-window', type=(int, int), default=(-70, 150), show_default=True,
              help='Time window for waveform data around P-arrivals to extract. Has no effect without '
                   '--p-data')
@click.option('--s-data-window', type=(int, int), default=(-100, 150), show_default=True,
              help='Time window for waveform data around S-arrivals to extract. Has no effect without '
                   '--s-data')
@click.option('--sw-data-window', type=(int, int), default=(-70, 4*60*60), show_default=True,
              help='Time window for waveform data around SW-arrivals to extract. Has no effect without '
                   '--sw-data')
@click.option('--p-distance-range', type=(int, int), default=(30, 90), show_default=True,
              help='Range of epicentral distances (in degrees) for which P-arrival data at a given station '
                   'are to be fetched. Has no effect without --p-data')
@click.option('--s-distance-range', type=(int, int), default=(55, 85), show_default=True,
              help='Range of epicentral distances (in degrees) for which S-arrival data at a given station '
                   'are to be fetched. Has no effect without --s-data')
@click.option('--sw-distance-range', type=(int, int), default=(5, 175), show_default=True,
              help='Range of epicentral distances (in degrees) for which SW-arrival data at a given station '
                   'are to be fetched. Has no effect without --sw-data')
@click.option('--p-resample-hz', type=float, default=10, show_default=True,
              help='Resampling frequency (default 10 Hz) for output P traces')
@click.option('--s-resample-hz', type=float, default=10, show_default=True,
              help='Resampling frequency (default 10 Hz) for output S traces')
@click.option('--sw-resample-hz', type=float, default=2, show_default=True,
              help='Resampling frequency (default 2 Hz) for surface waves')
@click.option('--taup-model', type=str, default='iasp91', show_default=True,
              help='Theoretical tau-p Earth model to use for Trace stats computation. Other possibilities, '
                   'such as ak135, are documented here: https://docs.obspy.org/packages/obspy.taup.html')
@click.option('--dry-run', is_flag=True, default=False, show_default=True,
              help='Reports events available to each station, by wave-type and exits without outputting any data. ')
def main(data_source, network_list, station_list, gcmt_catalog_file, output_file, log_folder,
         start_time, end_time,
         p_data, s_data, sw_data,
         p_magnitude_range, s_magnitude_range, sw_magnitude_range,
         p_data_window, s_data_window, sw_data_window,
         p_distance_range, s_distance_range, sw_distance_range,
         p_resample_hz, s_resample_hz, sw_resample_hz,
         taup_model, dry_run):
    """
    DATA_SOURCE: Text file containing paths to ASDF files.
    """
    
    # Initialize MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    output_fn_base = os.path.splitext(os.path.basename(output_file))[0]
    log = setup_logger('__func__', os.path.join(log_folder, output_fn_base + '.log'))

    # sanity check
    owave_types = defaultdict(bool)
    if(not(p_data or s_data or sw_data)):
        assert 0, 'At least one from [--p-data, --s-data, --sw-data] must be specified. Aborting'
    else:
        owave_types['P'] = p_data
        owave_types['S'] = s_data
        owave_types['SW'] = sw_data
    # end if

    # initialize event magnitude range dict
    magnitude_range = defaultdict(tuple)
    magnitude_range['P'] = p_magnitude_range
    magnitude_range['S'] = s_magnitude_range
    magnitude_range['SW'] = sw_magnitude_range

    # initialize event data window dict
    request_window = defaultdict(tuple) # seconds
    request_window['P'] = p_data_window
    request_window['S'] = s_data_window
    request_window['SW'] = sw_data_window

    # initialize event distance range dict
    distance_range = defaultdict(tuple) # arc degrees
    distance_range['P'] = p_distance_range
    distance_range['S'] = s_distance_range
    distance_range['SW'] = sw_distance_range

    # initialize resampling dict
    resample_hz = defaultdict(tuple) # arc degrees
    resample_hz['P'] = p_resample_hz
    resample_hz['S'] = s_resample_hz
    resample_hz['SW'] = sw_resample_hz

    # initialize depth range dict
    depth_range = defaultdict(tuple) # arc degrees
    depth_range['P'] = [0, np.finfo('f4').max]
    depth_range['S'] = [0, np.finfo('f4').max]
    depth_range['SW'] = [0, 150] # max depth of 150 km

    # initialize areal event separation dict
    areal_separation_km = defaultdict(tuple) # arc degrees
    areal_separation_km['P'] = 0
    areal_separation_km['S'] = 0
    areal_separation_km['SW'] = 100 # events within 15 minutes of each other, of similar magnitude,
                                    # should be separated by at least 100 km

    inventory = None
    pad = 10 # nominal padding for waveforms in seconds
    fds = FederatedASDFDataSet(data_source)

    #################################################
    # Check if GPS clock-corrections are being applied
    # A large padding is used to allow for time-shifts
    # from clock-correction
    #################################################
    if (fds.corrections_enabled()): pad = 3600

    # trim inventory based on inputs
    log.info('Loading inventory...')
    inventory = fds.get_inventory()

    log.info('Trimming inventory...')
    inventory = trim_inventory(inventory, network_list=network_list, station_list=station_list)
    netsta_df = DataFrame(columns=['net.sta', 'lon', 'lat'])
    netsta_count = 0
    for net in inventory.networks:
        nc = net.code
        for sta in net.stations:
            sc = sta.code
            netsta = '{}.{}'.format(nc, sc)
            netsta_df.loc[netsta_count] = [netsta, *fds.unique_coordinates[netsta]]
            netsta_count += 1
        # end for
    # end for

    if(len(netsta_df) == 0):
        log.error('Inventory is empty! Aborting..')
        parallel_abort('')
    else:
        netsta_df.index += 1
        log.info('Inventory contains a total of {} stations: \n {}\n'.format(netsta_count,
                                                                             netsta_df.to_string()))
    # end if
    log.info('Loading GCMT catalog: {}..'.format(gcmt_catalog_file))
    catalog = GCMTCatalog(gcmt_catalog_file)

    if(rank == 0):
        assert not os.path.exists(output_file), \
            "Output file {} already exists, please remove!".format(output_file)
        log.info("Traces will be written to: {}\n".format(output_file))
    # end if

    # define closures for getting recording timespans and waveforms
    def recording_timespan_getter(network, station, location):
        return fds.get_recording_timespan(network=network, station=station, location=location)
    # end func

    def waveform_getter(network, station, location, channel, starttime, endtime):
        return asdf_get_waveforms(fds, network, station, location, channel, starttime, endtime)
    # end func

    for wave, flag in owave_types.items():
        if(not flag): continue

        if(rank == 0): log.info('Processing {} events..'.format(DESCS[wave]))
        extract_data(recording_timespan_getter, waveform_getter, catalog,
                     inventory, output_file, log_folder,
                     wave, request_window[wave], [start_time, end_time],
                     distance_range[wave], magnitude_range[wave],
                     depth_range[wave], areal_separation_km[wave],
                     resample_hz[wave], tt_model=taup_model, pad=pad, dry_run=dry_run)
    # end for

    del fds
    if(rank == 0):
        log.info("extract_event_traces SUCCESS!")
    # end if
# end main

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
# end if
