"""
Description:
    Implements the earthquake cross-correlation workflow

References:

CreationDate:   30/06/23
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     30/06/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os
import glob
from collections import defaultdict
from math import sqrt

from ordered_set import OrderedSet as set
import numpy as np
import random
import click
import re
from mpi4py import MPI
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from obspy import UTCDateTime, read_inventory, Inventory
from obspy.geodetics.base import gps2dist_azimuth

from seismic.xcorqc.xcorqc import IntervalStackXCorr
from seismic.xcorqc.utils import getStationInventory, read_location_preferences, Dataset, get_stream
from seismic.misc import get_git_revision_hash, rtp2xyz, split_list
from seismic.misc_p import ProgressTracker
from seismic.eqcorr.utils import GCMTCatalog
from obspy.core import Stream
from pandas import DataFrame
from pandas.core.series import Series

def process(data_source1, data_source2, gcmt_catalog, output_path,
            window_seconds, resample_rate=None, taper_length=0.05, nearest_neighbours=1,
            fmin=None, fmax=None, netsta_list1='*', netsta_list2='*', pairs_to_compute=None,
            instrument_response_inventory=None, instrument_response_output='vel', water_level=50,
            clip_to_2std=False, whitening=False, whitening_window_frequency=0,
            one_bit_normalize=False, location_preferences=None,
            ds1_zchan=None, ds1_nchan=None, ds1_echan=None,
            ds2_zchan=None, ds2_nchan=None, ds2_echan=None, corr_chan=None,
            restart=False, dry_run=False):

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
    git_hash = ''
    if (rank == 0):
        # get git-hash
        git_hash = get_git_revision_hash()

        def outputConfigParameters():
            # output config parameters
            fn = 'eq_correlator.cfg'
            fn = os.path.join(output_path, fn)

            f = open(fn, 'w+')
            f.write('Parameters Values:\n\n')
            f.write('%35s\t\t\t: %s\n' % ('DATA_SOURCE1', data_source1))
            f.write('%35s\t\t\t: %s\n' % ('DATA_SOURCE2', data_source2))
            f.write('%35s\t\t\t: %s\n' % ('GCMT_CATALOG', gcmt_catalog))
            f.write('%35s\t\t\t: %s\n' % ('OUTPUT_PATH', output_path))
            f.write('%35s\t\t\t: %s\n\n' % ('WINDOW_SECONDS', window_seconds))

            f.write('%35s\t\t\t: %s\n' % ('--resample-rate', resample_rate))
            f.write('%35s\t\t\t: %s\n' % ('--taper-length', taper_length))
            f.write('%35s\t\t\t: %s\n' % ('--nearest-neighbours', nearest_neighbours))
            f.write('%35s\t\t\t: %s\n' % ('--fmin', fmin))
            f.write('%35s\t\t\t: %s\n' % ('--fmax', fmax))
            f.write('%35s\t\t\t: %s\n' % ('--station-names1', netsta_list1))
            f.write('%35s\t\t\t: %s\n' % ('--station-names2', netsta_list2))
            f.write('%35s\t\t\t: %s\n' % ('--instrument-response-inventory', instrument_response_inventory))
            f.write('%35s\t\t\t: %s\n' % ('--instrument-response-output', instrument_response_output))
            f.write('%35s\t\t\t: %s\n' % ('--corr-chan', corr_chan))
            f.write('%35s\t\t\t: %s\n' % ('--water-level', water_level))
            f.write('%35s\t\t\t: %s\n' % ('--clip-to-2std', clip_to_2std))
            f.write('%35s\t\t\t: %s\n' % ('--one-bit-normalize', one_bit_normalize))
            f.write('%35s\t\t\t: %s\n' % ('--whitening', whitening))
            if(whitening):
                f.write('%35s\t\t\t: %s\n' % ('--whitening-window-frequency', whitening_window_frequency))
            f.write('%35s\t\t\t: %s\n' % ('--restart', 'TRUE' if restart else 'FALSE'))

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

        pairs = ds1.get_unique_station_pairs(ds2, nn=nearest_neighbours, require_overlap=False)
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

    # instantiate GCMTCatalog
    gcat = GCMTCatalog(gcmt_catalog)

    start_time = UTCDateTime('1900-01-01')
    end_time = UTCDateTime('2100-01-01')
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

        #######################################################
        # Evaluate channels
        #######################################################
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
        else: raise ValueError('Invalid corr-chan')

        if(len(corr_chans)<2):
            print('Required channels are not found for station %s or %s..'%(netsta1, netsta2))
            continue
        # end if

        #######################################################
        # Fetch data and compute cross-correlation
        #######################################################
        lon1, lat1 = ds1.fds.unique_coordinates[netsta1]
        lon2, lat2 = ds1.fds.unique_coordinates[netsta2]

        compatible_events = gcat.get_compatible_events(lon1, lat1, lon2, lat2, min_magnitude=5.8)
        print(len(compatible_events))
        for k, v in sorted(compatible_events.items(), key=lambda kv: kv[1]):
            net1, sta1 = netsta1.split('.')
            net2, sta2 = netsta2.split('.')
            cha1, cha2 = corr_chans
            e1 = gcat.cat.iloc[k[0]]
            e2 = gcat.cat.iloc[k[1]]

            e1s1 = e1s2 = e2s1 = e2s2 = None
            try:
                e1s1 = get_stream(ds1.fds, net1, sta1, cha1,
                                  e1['EventOrigintim'],
                                  e1['EventOrigintim'] + window_seconds,
                                  location_preferences_dict)
                e1s2 = get_stream(ds2.fds, net2, sta2, cha2,
                                  e1['EventOrigintim'],
                                  e1['EventOrigintim'] + window_seconds,
                                  location_preferences_dict)
                e2s1 = get_stream(ds1.fds, net1, sta1, cha1,
                                  e2['EventOrigintim'],
                                  e2['EventOrigintim'] + window_seconds,
                                  location_preferences_dict)
                e2s2 = get_stream(ds2.fds, net2, sta2, cha2,
                                  e2['EventOrigintim'],
                                  e2['EventOrigintim'] + window_seconds,
                                  location_preferences_dict)
            except Exception as e:
                print(e)
                pass
            # end try

            if(e1s1 and e1s2):
                process_event_pair(e1, e1, e1s1, e1s2, window_seconds)
            if(e2s1 and e2s2):
                process_event_pair(e2, e2, e2s1, e2s2, window_seconds)
            if(e1s1 and e2s2):
                process_event_pair(e1, e2, e1s1, e2s2, window_seconds)
            if(e2s1 and e1s2):
                process_event_pair(e2, e1, e2s1, e1s2, window_seconds)
        # end for
    # end for
# end func

def process_event_pair(e1:Series, e2:Series, s1:Stream, s2:Stream,
                       window_seconds:float):
    """
    @param e1: a row in a GCMTCatalog instance
    @param e2: a row in a GCMTCatalog instance
    @param s1: stream from station 1
    @param s2: stream for station 2
    @return:
    """

    def compute_cc():
        def align_traces(a, b):
            tr1 = a.copy()
            tr2 = b.copy()

            clen = np.max([tr1.data.shape[0], tr2.data.shape[0]])
            if(tr1.data.shape[0]<clen): tr1.data = zeropad(tr1.data, clen)
            if(tr2.data.shape[0]<clen): tr2.data = zeropad(tr2.data, clen)

            tr1.stats.starttime=UTCDateTime("1970-01-01")
            tr2.stats.starttime=UTCDateTime("1970-01-01")

            return tr1, tr2, int(np.floor(clen*tr1.stats.delta))
        # end func

        def pad_traces1(a, b):
            def zeropad_se(tr, padlen):
                assert (tr.shape[0] < padlen)
                padded = np.zeros(padlen)
                s = int((padlen - tr.shape[0]) / 2)
                padded[s:(s + tr.shape[0])] = tr
                return padded
            # end func

            tr1 = a.copy()
            tr2 = b.copy()

            clen = np.max([tr1.shape[0], tr2.shape[0]])
            if(tr1.shape[0]<clen): tr1 = zeropad_se(tr1, clen)
            if(tr2.shape[0]<clen): tr2 = zeropad_se(tr2, clen)

            return tr1, tr2
        # end func

        def sw_trace_window(tr, slon, slat, elon, elat):
            tr = tr.copy()
            MINVEL = 2.5 #km/s
            MAXVEL = 7.5 #km/s
            _, _, sta_dist = geod.inv(slon, slat,
                                      elon, elat)
            sr = tr.stats.sampling_rate
            tw_s, tw_e = sta_dist/1e3/MAXVEL, sta_dist/1e3/MINVEL
            tw = tukey(int((tw_e-tw_s)*sr), alpha=0.25)

            w = np.zeros(len(tr.data), dtype='f4')
            w[int(tw_s*sr) : int(tw_s*sr+len(tw))] = tw

            tr.data = tr.data * w
            return tr, w
        # end func

        eq1 = cat.iloc[params[i,0]]
        eq2 = cat.iloc[params[i,1]]

        tr1, tr2, window_seconds = pad_traces(meek_stream[i*2], qlp_stream[i*2])

        tr1w, w1 = sw_trace_window(tr1, *sta_coords_dict['MEEK'], eq1['lon'], eq1['lat'])
        tr2w, w2 = sw_trace_window(tr2, *sta_coords_dict['QLP'], eq1['lon'], eq1['lat'])

        cf1, _, _, _, sr =xcorr2(tr1w,
                                 tr2w,
                                 window_seconds=window_seconds,
                                 interval_seconds=window_seconds, whitening=True,
                                 whitening_window_frequency=0.02,
                                 flo=flo, fhi=fhi, resample_rate=4, verbose=10)
        dt = 1./float(sr)
        x1 = np.linspace(-window_seconds + dt, window_seconds - dt, len(cf1[0]))

        tr1, tr2, window_seconds = pad_traces(meek_stream[i*2], qlp_stream[i*2+1])

        tr1w, w1 = sw_trace_window(tr1, *sta_coords_dict['MEEK'], eq1['lon'], eq1['lat'])
        tr2w, w2 = sw_trace_window(tr2, *sta_coords_dict['QLP'], eq2['lon'], eq2['lat'])

        cf2, _, _, _, sr =xcorr2(tr1w,
                                 tr2w,
                                 window_seconds=window_seconds,
                                 interval_seconds=window_seconds, whitening=True,
                                 whitening_window_frequency=0.02,
                                 flo=flo, fhi=fhi, resample_rate=4, verbose=10)
        x2 = np.linspace(-window_seconds + dt, window_seconds - dt, len(cf2[0]))

        ci = params[i,0]
        cj = params[i,1]
        print(ci,cj)
        mt_angle = np.degrees(np.arccos(np.min([np.dot(cat.iloc[ci,4:10],
                                                       cat.iloc[cj,4:10]) / \
                   (np.linalg.norm(cat.iloc[ci,4:10]) *
                    np.linalg.norm(cat.iloc[cj,4:10])), 1.])))

        # compute lag at peak
        lag = x2[np.argmax(np.fabs(cf2))]

        # compute waveform similarity
        #print(cf1[0].shape, cf2[0].shape)
        a, b = pad_traces1(cf1[0], cf2[0])
        #plt.plot(a, c='r')
        #plt.plot(b, c='g')
        #a /= np.linalg.norm(a)
        #b /= np.linalg.norm(b)
        rms_fit = 1-np.sqrt(np.sum((a-b)**2)/np.sum(a**2))

        # compute lag offset
        cft1 = Trace(cf1[0], header={'sampling_rate':4})
        cft1.filter('bandpass', freqmin=.005, freqmax=0.05, zerophase=True)

        cft2 = Trace(cf2[0], header={'sampling_rate':4})
        cft2.filter('bandpass', freqmin=.005, freqmax=0.05, zerophase=True)

        c = np.correlate(cft1.data, cft2.data, "full")
        x3 = np.arange(-len(c)//2+1, len(c)//2-1)*dt
        lag_offset = np.fabs(x3[np.argmax(np.fabs(c))])

        return tr1, tr2, w1, w2, x1, x2, cft1, cft2, lag, lag_offset, mt_angle, rms_fit
    # end func

    if ((s1[0].stats.endtime - s1[0].stats.starttime) >= window_seconds and \
        (s2[0].stats.endtime - s2[0].stats.starttime) >= window_seconds):
        print(s1[0].stats.station, UTCDateTime(e1['EventOrigintim']), '----',
              s2[0].stats.station, UTCDateTime(e2['EventOrigintim']),
              'angle: {}'.format(GCMTCatalog.get_mt_angle(e1, e2)))
    # end func
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'], show_default=True)
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('data-source1',
                type=click.Path('r'))
@click.argument('data-source2',
                type=click.Path('r'))
@click.argument('gcmt-catalog',
                type=click.Path('r'))
@click.argument('output-path', required=True,
                type=click.Path(exists=True))
@click.argument('window-seconds', required=True,
                type=int)
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
@click.option('--restart', default=False, is_flag=True, help='Restart job')
@click.option('--dry-run', default=False, is_flag=True, help='Dry run for printing out station-pairs and '
                                                             'additional stats.')
def main(data_source1, data_source2, gcmt_catalog, output_path, window_seconds,
         resample_rate, taper_length, nearest_neighbours, fmin, fmax, station_names1, station_names2,
         pairs_to_compute, instrument_response_inventory, instrument_response_output, water_level,
         clip_to_2std, whitening, whitening_window_frequency, one_bit_normalize, location_preferences,
         ds1_zchan, ds1_nchan, ds1_echan, ds2_zchan, ds2_nchan, ds2_echan, corr_chan, restart, dry_run):
    """
    DATA_SOURCE1: Text file containing paths to ASDF files \n
    DATA_SOURCE2: Text file containing paths to ASDF files \n
    GCMT_CATALOG: GCMT catalog file name in text format. The expected columns (space-separated) are:
    lonp lon lat dep mrr mtt mff mrt mrf mtf exp EventOrigintim DC CLVD VOL Mw str1 dip1 rak1 \
    str2 dip2 rak2 plunP azP plunT azT Hlonp Hlon Hlat Hdep Ctime Chdur MwG base Bst Bco Bmp Bper \
    Sst Sco Smp Sper Mst Mco Mmp Mper\n
    OUTPUT_PATH: Output folder \n
    WINDOW_SECONDS: Length of time window (s); e.g. 3600 for an hour\n
    """

    if(resample_rate): resample_rate = float(resample_rate)
    if(fmin): fmin = float(fmin)
    if(fmax): fmax = float(fmax)

    process(data_source1, data_source2, gcmt_catalog, output_path, window_seconds,
            resample_rate, taper_length, nearest_neighbours,
            fmin, fmax, station_names1, station_names2, pairs_to_compute,
            instrument_response_inventory, instrument_response_output, water_level, clip_to_2std, whitening,
            whitening_window_frequency, one_bit_normalize, location_preferences, ds1_zchan, ds1_nchan,
            ds1_echan, ds2_zchan, ds2_nchan, ds2_echan, corr_chan, restart, dry_run)
# end func

if __name__ == '__main__':
    main()  # pylint: disable=no-value-for-parameter
# end if