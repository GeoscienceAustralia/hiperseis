#!/bin/env python
"""
Description:
    Generates cross-correlations for data from staion-pairs in parallel

References:

CreationDate:   11/07/18

Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     11/07/18   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
import glob
from collections import defaultdict
from math import sqrt

from ordered_set import OrderedSet as set
import numpy as np
from scipy.spatial import cKDTree
import random
import click
import re
from mpi4py import MPI
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from obspy import UTCDateTime, read_inventory, Inventory
from obspy.geodetics.base import gps2dist_azimuth

from seismic.xcorqc.xcorqc import IntervalStackXCorr
from seismic.xcorqc.utils import ProgressTracker, getStationInventory, rtp2xyz, split_list
from seismic.misc import get_git_revision_hash

class Dataset:
    def __init__(self, asdf_file_name, netsta_list='*'):

        self._data_path = asdf_file_name
        self._earth_radius = 6371  # km

        self.fds = FederatedASDFDataSet(asdf_file_name)
        # Gather station metadata
        netsta_list_subset = set(netsta_list.split(' ')) if netsta_list != '*' else netsta_list
        self.netsta_list = []
        self.metadata = defaultdict(list)

        rtps = []
        for netsta in list(self.fds.unique_coordinates.keys()):
            if(netsta_list_subset != '*'):
                if netsta not in netsta_list_subset:
                    continue

            self.netsta_list.append(netsta)
            self.metadata[netsta] = self.fds.unique_coordinates[netsta]

            rtps.append([self._earth_radius,
                         np.radians(90 - self.metadata[netsta][1]),
                         np.radians(self.metadata[netsta][0])])
        # end for

        if(len(rtps) == 0):
            assert 0, 'No station-pairs found due to missing stations. Aborting..'
        # end if

        rtps = np.array(rtps)
        xyzs = rtp2xyz(rtps[:, 0], rtps[:, 1], rtps[:, 2])

        self._tree = cKDTree(xyzs)
        self._cart_location = defaultdict(list)
        for i, ns in enumerate(self.netsta_list):
            self._cart_location[ns] = xyzs[i, :]
        # end for
    # end func

    def get_closest_stations(self, netsta, other_dataset, nn=1):
        assert isinstance(netsta, str), 'station_name must be a string'
        assert isinstance(other_dataset, Dataset), 'other_dataset must be an instance of Dataset'
        netsta = netsta.upper()

        assert netsta in self.netsta_list, '%s not found'%(netsta)

        d, l = other_dataset._tree.query(self._cart_location[netsta], nn)

        if isinstance(l, int):
            l = np.array([l])

        l = l[l<len(other_dataset.netsta_list)]

        if isinstance(l, int):
            l = np.array([l])

        assert len(l), 'No stations found..'

        return list(np.array(other_dataset.netsta_list)[l])
    # end func

    def get_unique_station_pairs(self, other_dataset, nn=1):
        pairs = set()
        for st1 in self.netsta_list:
            st2list = None
            if (nn != -1):
                if self == other_dataset:
                    st2list = set(self.get_closest_stations(st1, other_dataset, nn=nn + 1))
                    if st1 in st2list:
                        st2list.remove(st1)
                    st2list = list(st2list)
                else:
                    st2list = self.get_closest_stations(st1, other_dataset, nn=nn)
            else:
                st2list = other_dataset.netsta_list
            # end if

            for st2 in st2list:
                pairs.add((st1, st2))
            # end for
        # end for

        pairs_subset = set()
        for item in pairs:
            if(item[0] == item[1]): continue

            dup_item = (item[1], item[0])
            if(dup_item not in pairs_subset and item not in pairs_subset):
                pairs_subset.add(item)
            # end if
        # end if

        return list(pairs_subset)
    # end func
# end class

def read_location_preferences(location_preferences_fn):
    result = defaultdict(lambda: None)

    if(location_preferences_fn):
        pref_list = open(location_preferences_fn, 'r').readlines()

        for pref in pref_list:
            pref = pref.strip()
            if (len(pref)):
                try:
                    netsta, loc = pref.split()
                    net, sta = netsta.split('.')

                    result[netsta] = loc
                except Exception as e:
                    print(str(e))
                    assert 0, 'Error parsing: {}'.format(pref)
                # end try
            # end if
        # end for
    # end if
    return result
# end func

def process(data_source1, data_source2, output_path,
            interval_seconds, window_seconds, window_overlap, window_buffer_length,
            read_ahead_window_seconds, resample_rate=None, taper_length=0.05, nearest_neighbours=1,
            fmin=None, fmax=None, netsta_list1='*', netsta_list2='*', pairs_to_compute=None,
            start_time='1970-01-01T00:00:00', end_time='2100-01-01T00:00:00',
            instrument_response_inventory=None, instrument_response_output='vel', water_level=50,
            clip_to_2std=False, whitening=False, whitening_window_frequency=0,
            one_bit_normalize=False, location_preferences=None,
            ds1_zchan=None, ds1_nchan=None, ds1_echan=None,
            ds2_zchan=None, ds2_nchan=None, ds2_echan=None, corr_chan=None,
            envelope_normalize=False, ensemble_stack=False, apply_stacking=True,
            restart=False, dry_run=False, no_tracking_tag=False, scratch_folder=None):
    """
    :param data_source1: Text file containing paths to ASDF files
    :param data_source2: Text file containing paths to ASDF files
    :param output_path: Output folder
    :param interval_seconds: Length of time window (s) over which to compute cross-correlations; e.g. 86400 for 1 day
    :param window_seconds: Length of stacking window (s); e.g 3600 for an hour. interval_seconds must be a multiple of \
                    window_seconds; no stacking is performed if they are of the same size.
    """
    if(os.path.exists(netsta_list1)):
        netsta_list1 = ' '.join(open(netsta_list1).readlines()).replace('\n', ' ').strip()
    # end if
    if(os.path.exists(netsta_list2)):
        netsta_list2 = ' '.join(open(netsta_list2).readlines()).replace('\n', ' ').strip()
    # end if

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    ds1 = Dataset(data_source1, netsta_list1)
    ds2 = Dataset(data_source2, netsta_list2)

    proc_stations = []
    location_preferences_dict = None
    time_tag = None
    git_hash = ''
    if (rank == 0):
        # get git-hash
        git_hash = get_git_revision_hash()

        # Register time tag with high resolution, since queued jobs can readily
        # commence around the same time.
        if(no_tracking_tag):
            time_tag = None
        else:
            time_tag = UTCDateTime.now().strftime("%y-%m-%d.T%H.%M.%S.%f")

        def outputConfigParameters():
            # output config parameters
            fn = 'correlator.%s.cfg' % (time_tag) if time_tag else 'correlator.cfg'
            fn = os.path.join(output_path, fn)

            f = open(fn, 'w+')
            f.write('Parameters Values:\n\n')
            f.write('%35s\t\t\t: %s\n' % ('DATA_SOURCE1', data_source1))
            f.write('%35s\t\t\t: %s\n' % ('DATA_SOURCE2', data_source2))
            f.write('%35s\t\t\t: %s\n' % ('OUTPUT_PATH', output_path))
            f.write('%35s\t\t\t: %s\n\n' % ('WINDOW_SECONDS', window_seconds))
            f.write('%35s\t\t\t: %s\n\n' % ('WINDOW_OVERLAP', window_overlap))
            f.write('%35s\t\t\t: %s\n\n' % ('READ_AHEAD_WINDOW_SECONDS', read_ahead_window_seconds))

            f.write('%35s\t\t\t: %s\n' % ('--window-buffer-length', window_buffer_length))
            f.write('%35s\t\t\t: %s\n' % ('--resample-rate', resample_rate))
            f.write('%35s\t\t\t: %s\n' % ('--taper-length', taper_length))
            f.write('%35s\t\t\t: %s\n' % ('--nearest-neighbours', nearest_neighbours))
            f.write('%35s\t\t\t: %s\n' % ('--fmin', fmin))
            f.write('%35s\t\t\t: %s\n' % ('--fmax', fmax))
            f.write('%35s\t\t\t: %s\n' % ('--station-names1', netsta_list1))
            f.write('%35s\t\t\t: %s\n' % ('--station-names2', netsta_list2))
            f.write('%35s\t\t\t: %s\n' % ('--start-time', start_time))
            f.write('%35s\t\t\t: %s\n' % ('--end-time', end_time))
            f.write('%35s\t\t\t: %s\n' % ('--instrument-response-inventory', instrument_response_inventory))
            f.write('%35s\t\t\t: %s\n' % ('--instrument-response-output', instrument_response_output))
            f.write('%35s\t\t\t: %s\n' % ('--corr-chan', corr_chan))
            f.write('%35s\t\t\t: %s\n' % ('--water-level', water_level))
            f.write('%35s\t\t\t: %s\n' % ('--clip-to-2std', clip_to_2std))
            f.write('%35s\t\t\t: %s\n' % ('--one-bit-normalize', one_bit_normalize))
            f.write('%35s\t\t\t: %s\n' % ('--envelope-normalize', envelope_normalize))
            f.write('%35s\t\t\t: %s\n' % ('--whitening', whitening))
            if(whitening):
                f.write('%35s\t\t\t: %s\n' % ('--whitening-window-frequency', whitening_window_frequency))
            f.write('%35s\t\t\t: %s\n' % ('--ensemble-stack', ensemble_stack))
            f.write('%35s\t\t\t: %s\n' % ('--restart', 'TRUE' if restart else 'FALSE'))
            f.write('%35s\t\t\t: %s\n' % ('--no-tracking-tag', 'TRUE' if no_tracking_tag else 'FALSE'))
            f.write('%35s\t\t\t: %s\n' % ('--scratch-folder', scratch_folder))

            f.close()
        # end func

        def cull_pairs(pairs, keep_list_fn):
            result = set()
            pairs_set = set()

            for pair in pairs:
                pairs_set.add('%s.%s'%(pair[0], pair[1]))
            # end for

            keep_list = open(keep_list_fn, 'r').readlines()
            for keep_pair in keep_list:
                keep_pair = keep_pair.strip()
                if(len(keep_pair)):
                    try:
                        knet1, ksta1, knet2, ksta2 = keep_pair.split('.')

                        keep_pair_alt = '%s.%s.%s.%s'%(knet2, ksta2, knet1, ksta1)

                        if(keep_pair in pairs_set or keep_pair_alt in pairs_set):
                            result.add(('%s.%s'%(knet1, ksta1), '%s.%s'%(knet2, ksta2)))
                    except Exception as e:
                        print(str(e))
                        assert 0, 'Error parsing: {}'.format(keep_pair)
                    # end try
                # end if
            # end for

            return list(result)
        # end func

        outputConfigParameters()

        pairs = ds1.get_unique_station_pairs(ds2, nn=nearest_neighbours)
        if(pairs_to_compute):
            # only keep pairs provided in the text file, given they exist in the data-sets
            pairs = cull_pairs(pairs, pairs_to_compute)
        # end if

        # print out station-pairs for dry runs
        if(dry_run):
            print('Computing %d station-pairs: '%(len(pairs)))
            for pair in pairs:
                print('.'.join(pair))
            # end for
        # end if

        random.Random(nproc).shuffle(pairs) # using nproc as seed so that shuffle produces the same
                                            # ordering when jobs are restarted.
        proc_stations = split_list(pairs, npartitions=nproc)
    # end if

    if(dry_run):
        # nothing more to do for dry runs
        return
    # end if

    location_preferences_dict = read_location_preferences(location_preferences)
    # broadcast workload to all procs
    proc_stations = comm.bcast(proc_stations, root=0)
    time_tag = comm.bcast(time_tag, root=0)
    git_hash = comm.bcast(git_hash, root=0)

    # read inventory
    inv = None
    stationInvCache = defaultdict(list)
    if(instrument_response_inventory):
        try:
            inv = read_inventory(instrument_response_inventory)
        except Exception as e:
            print (e)
        # end try
    # end if

    # Progress tracker
    progTracker = ProgressTracker(output_folder=output_path, restart_mode=restart)

    startTime = UTCDateTime(start_time)
    endTime = UTCDateTime(end_time)
    for pair in proc_stations[rank]:
        netsta1, netsta2 = pair

        if (progTracker.increment()):
            pass
        else:
            print (('Found results for station-pair: %s.%s. Moving along..'%(netsta1, netsta2)))
            continue
        # end if

        netsta1inv, stationInvCache = getStationInventory(inv, stationInvCache, netsta1, location_preferences_dict)
        netsta2inv, stationInvCache = getStationInventory(inv, stationInvCache, netsta2, location_preferences_dict)

        def evaluate_channels(cha1, cha2):
            result = []
            for netsta, cha, ds in zip((netsta1, netsta2), (cha1, cha2), (ds1, ds2)):
                if('*' not in cha1):
                    result.append(cha)
                else:
                    cha = cha.replace('*', '.*')  # hack to capture simple regex comparisons

                    net, sta = netsta.split('.')
                    stations = ds.fds.get_stations(start_time, end_time, net, sta)
                    for item in stations:
                        if(re.match(cha, item[3])):
                            result.append(item[3])
                            break
                        # end if
                    # end if
                # end if
            # end for

            return result
        # end func

        corr_chans = []
        if   (corr_chan == 'z'): corr_chans = evaluate_channels(ds1_zchan, ds2_zchan)
        elif (corr_chan == 'n'): corr_chans = evaluate_channels(ds1_nchan, ds2_nchan)
        elif (corr_chan == 'e'): corr_chans = evaluate_channels(ds1_echan, ds2_echan)
        elif (corr_chan == 't'): corr_chans = ['00T', '00T']
        else: raise ValueError('Invalid corr-chan')

        if(len(corr_chans)<2):
            print(('Either required channels are not found for station %s or %s, '
                   'or no overlapping data exists..')%(netsta1, netsta2))
            continue
        # end if

        baz_netsta1 = None
        baz_netsta2 = None
        if(corr_chan == 't'):
            try:
                sta1_lon, sta1_lat = ds1.fds.unique_coordinates[netsta1]
                sta2_lon, sta2_lat = ds2.fds.unique_coordinates[netsta2]
                _, baz_netsta2, baz_netsta1 = gps2dist_azimuth(sta1_lat, sta1_lon, sta2_lat, sta2_lon)
            except Exception as e:
                print (e)
                print (('Failed to compute back-azimuth for station-pairs; skipping %s.%s; '%(netsta1, netsta2)))
                continue
            # end try
        # end if

        IntervalStackXCorr(ds1.fds, ds2.fds, startTime,
                           endTime, netsta1, netsta2, netsta1inv, netsta2inv,
                           instrument_response_output, water_level,
                           corr_chans[0], corr_chans[1],
                           baz_netsta1, baz_netsta2, location_preferences_dict,
                           resample_rate, taper_length, read_ahead_window_seconds,
                           interval_seconds, window_seconds, window_overlap,
                           window_buffer_length, fmin, fmax, clip_to_2std, whitening,
                           whitening_window_frequency, one_bit_normalize, envelope_normalize,
                           ensemble_stack, apply_stacking, output_path, 2, time_tag,
                           scratch_folder, git_hash)
    # end for
# end func


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'], show_default=True)

@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('data-source1',
                type=click.Path('r'))
@click.argument('data-source2',
                type=click.Path('r'))
@click.argument('output-path', required=True,
                type=click.Path(exists=True))
@click.argument('window-seconds', required=True,
                type=int)
@click.argument('window-overlap', required=True,
                type=float)
@click.argument('read-ahead-windows', required=True,
                type=int)
@click.option('--stacking-interval-seconds', default=None, type=int,
              help="Interval, in seconds, to stack data windows over; default is None, "
                   "which outputs all cross-correlation windows, without performing any "
                   "stacking")
@click.option('--window-buffer-length', default=0, type=float,
              help="Length of buffer as a decimal percentage of 'window-seconds', at each end, around "
                    "actual data windows of interest. This helps exclude effects of tapering and "
                    "other edge artefacts from data windows before computing cross-correlations")
@click.option('--resample-rate', default=None, type=float,
              help="Resampling rate (Hz); applies to both datasets. Default is None, in which case, "
                   "data with a lower sampling rate is upsampled before computing cross-correlations")
@click.option('--taper-length', default=0.05, type=click.FloatRange(0, 0.5),
              help="Taper length as a decimal percentage of 'window-seconds', at each end; default 0.05")
@click.option('--nearest-neighbours', type=int, default=-1,
              help="Number of nearest neighbouring stations in data-source-2"
                   " to correlate against a given station in data-source-1. If"
                   " set to -1, correlations for a cartesian product of all stations"
                   " in both data-sets are produced -- note, this is computationally"
                   " expensive.")
@click.option('--fmin', default=None, type=float, help="Lowest frequency for bandpass filter; default is None")
@click.option('--fmax', default=None, type=float, help="Highest frequency for bandpass filter; default is None")
@click.option('--station-names1', default='*', type=str,
              help="Either station name(s) (space-delimited) or a text file containing NET.STA entries in each line to "
                   "process in data-source-1; default is '*', which processes all available stations.")
@click.option('--station-names2', default='*', type=str,
              help="Either station name(s) (space-delimited) or a text file containing NET.STA entries in each line to "
                   "process in data-source-2; default is '*', which processes all available stations.")
@click.option('--pairs-to-compute', default=None, type=click.Path('r'),
              help="Text file containing station pairs (NET.STA.NET.STA) for which cross-correlations are to be computed."
                   "Note that this parameter is intended as a way to restrict the number of computations to only the "
                   "station-pairs listed in the text-file.")
@click.option('--start-time', default='1970-01-01',
              type=str,
              help="Date and time (in UTC format) to start from")
@click.option('--end-time', default='2100-01-01',
              type=str,
              help="Date and time (in UTC format) to stop at")
@click.option('--instrument-response-inventory', default=None,
              type=click.Path('r'),
              help="FDSNxml inventory containing instrument response information. Note that when this parameter is provided, "
                   "instrument response corrections are automatically applied for matching stations with response "
                   "information.")
@click.option('--instrument-response-output',
              type=click.Choice(['vel', 'disp']),
              default='vel', help="Output of instrument response correction; must be either 'vel' (default) for velocity"
                                  " or 'disp' for displacement. Note, this parameter has no effect if instrument response"
                                  " correction is not performed.")
@click.option('--water-level', type=float, default=50.,
              help="Water-level in dB to limit amplification during instrument response correction "
                   "to a certain cut-off value. Note, this parameter has no effect if instrument "
                   "response correction is not performed.")
@click.option('--clip-to-2std', is_flag=True,
              help="Clip data in each window to +/- 2 standard deviations. Note that the default time-domain normalization "
                   "is N(0,1), i.e. 0-mean and unit variance")
@click.option('--whitening', is_flag=True,
              help="Apply spectral whitening")
@click.option('--whitening-window-frequency', type=float, default=0,
              help="Window frequency (Hz) determines the half-window length (of averaging window) used for smoothing weights "
                   "that scale the spectral amplitudes of the waveform being spectrally whitened. The default value of 0 "
                   "implies no smoothing of weights. Note that this parameter has no effect unless whitening is activated with "
                   "'--whitening'")
@click.option('--one-bit-normalize', is_flag=True,
              help="Apply one-bit normalization to data in each window.  Note that the default time-domain normalization "
                   "is N(0,1), i.e. 0-mean and unit variance")
@click.option('--location-preferences', default=None,
              type=click.Path('r'),
              help="A space-separated two-columned text file containing location code preferences for "
                   "stations in the form: 'NET.STA LOC'. Note that location code preferences need not "
                   "be provided for all stations -- the default is None. This approach allows "
                   "preferential selection of data from particular location codes for stations that "
                   "feature data from multiple location codes")
@click.option('--ds1-zchan', default='BHZ',
              type=str,
              help="Name of z-channel for data-source-1. This parameter and the five following are required to "
                   "specify channel names for the stations being cross-correlated. Simple wildcards, e.g. '*Z', are "
                   "also supported -- this is particularly useful when cross-correlating pairs of short-period and "
                   "broadband channels")
@click.option('--ds1-nchan', default='BHN',
              type=str,
              help="Name of n-channel for data-source-1")
@click.option('--ds1-echan', default='BHE',
              type=str,
              help="Name of e-channel for data-source-1")
@click.option('--ds2-zchan', default='BHZ',
              type=str,
              help="Name of z-channel for data-source-2")
@click.option('--ds2-nchan', default='BHN',
              type=str,
              help="Name of n-channel for data-source-2")
@click.option('--ds2-echan', default='BHE',
              type=str,
              help="Name of e-channel for data-source-2")
@click.option('--corr-chan', type=click.Choice(['z', 'n', 'e', 't']), default='z',
              help="Channels to be cross-correlated. Default is 'z' and 't' is for the transverse component"
                   ", rotated through Obspy NE->RT functionality, based on back-azimuth between station-pairs")
@click.option('--envelope-normalize', is_flag=True,
              help="Envelope (via Hilbert transform) and normalize cross-correlations. "
                    "This procedure is useful for clock-error detection "
                    "workflows")
@click.option('--ensemble-stack', is_flag=True, help="Outputs a single ensemble cross-correlation function over all data "
                                                     "for a given station-pair. In other words, stacks over "
                                                     "'interval-seconds' are in turn stacked to produce a "
                                                     "single cross-correlation function")
@click.option('--restart', default=False, is_flag=True, help='Restart job')
@click.option('--dry-run', default=False, is_flag=True, help='Dry run for printing out station-pairs and '
                                                             'additional stats.')
@click.option('--no-tracking-tag', default=False, is_flag=True, help='Do not tag output file names with a time-tag')
@click.option('--scratch-folder', default=None, help="Scratch folder for large jobs (e.g. $PBS_JOBFS on the NCI); "
                                                     "default is to use the standard temp folder")
def main(data_source1, data_source2, output_path, window_seconds, window_overlap, read_ahead_windows,
         stacking_interval_seconds, window_buffer_length, resample_rate, taper_length, nearest_neighbours,
         fmin, fmax, station_names1, station_names2, pairs_to_compute, start_time, end_time,
         instrument_response_inventory, instrument_response_output, water_level, clip_to_2std,
         whitening, whitening_window_frequency, one_bit_normalize, location_preferences,
         ds1_zchan, ds1_nchan, ds1_echan, ds2_zchan, ds2_nchan, ds2_echan, corr_chan, envelope_normalize,
         ensemble_stack, restart, dry_run, no_tracking_tag, scratch_folder):
    """
    DATA_SOURCE1: Text file containing paths to ASDF files \n
    DATA_SOURCE2: Text file containing paths to ASDF files \n
    OUTPUT_PATH: Output folder \n
    WINDOW_SECONDS: Length of time window (s); e.g. 3600 for an hour\n
    WINDOW_OVERLAP: Window overlap as a decimal percentage of WINDOW_SECONDS; e.g. 0.1 for 10% overlap\n
    READ_AHEAD_WINDOWS: An integer indicating the number of data windows that are to be read ahead, in
                        order to improve IO efficiency. Typically, a couple of days worth of data read
                        in at once has proven efficient for large runs; e.g. a value of 72 would roughly
                        fetch 3 days worth of data for WINDOW_SECONDS=3600. Note that the actual amount of
                        data fetched will depend on WINDOW_OVERLAP and --window-buffer-length\n
    """

    if(resample_rate): resample_rate = float(resample_rate)
    if(fmin): fmin = float(fmin)
    if(fmax): fmax = float(fmax)

    # sanity checks
    if(window_seconds <= 0): raise ValueError('WINDOW_SECONDS must be > 0')
    if(read_ahead_windows <= 0): raise ValueError('READ_AHEAD_WINDOWS must be > 0')
    if(window_buffer_length > 0.5): raise ValueError('--window-buffer-length must be <= 0.5')
    if(stacking_interval_seconds is not None):
        if(stacking_interval_seconds < window_seconds):
            raise ValueError('Invalid value for --stacking-interval-seconds, must be > WINDOW_SECONDS')
        if (stacking_interval_seconds > window_seconds * read_ahead_windows):
            raise ValueError("""Invalid value for --stacking-interval-seconds, \
                                must be < WINDOW_SECONDS*READ_AHEAD_WINDOWS""")
    # end if


    #######################################################
    # Compute amount of data to be read in in each IO call,
    # based on overlap and length of buffer around data
    # windows, if specified
    #######################################################
    interval_seconds = None
    read_ahead_window_seconds = None
    apply_stacking = None
    if(stacking_interval_seconds is None):
        if(ensemble_stack):
            raise ValueError('--ensemble-stack is only applicable with --stacking-interval-seconds. Aborting..')
        if(envelope_normalize):
            raise ValueError('--envelope-normalize is only applicable with --stacking-interval-seconds. Aborting..')

        read_ahead_window_seconds = window_seconds * (1 - window_overlap) * read_ahead_windows + \
                                    window_seconds * window_buffer_length * 2 + \
                                    window_overlap * window_seconds * 2
        interval_seconds = read_ahead_window_seconds
        apply_stacking = False
        #print(read_ahead_window_seconds)
    else:
        read_ahead_window_seconds = window_seconds * read_ahead_windows
        interval_seconds = stacking_interval_seconds
        apply_stacking = True
    # end if

    process(data_source1, data_source2, output_path, interval_seconds, window_seconds, window_overlap,
            window_buffer_length, read_ahead_window_seconds, resample_rate, taper_length, nearest_neighbours,
            fmin, fmax, station_names1, station_names2, pairs_to_compute, start_time, end_time,
            instrument_response_inventory, instrument_response_output, water_level, clip_to_2std, whitening,
            whitening_window_frequency, one_bit_normalize, location_preferences, ds1_zchan, ds1_nchan,
            ds1_echan, ds2_zchan, ds2_nchan, ds2_echan, corr_chan, envelope_normalize, ensemble_stack,
            apply_stacking, restart, dry_run, no_tracking_tag, scratch_folder)
# end func

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
# end if
