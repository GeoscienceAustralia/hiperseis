#!/bin/env python
"""
Description:
    Generates a waveform-analytics report on raw data in mseed format

References:

CreationDate:   18/01/2024
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     07/02/23   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

import os, sys

is_windows = sys.platform.startswith('win')
is_osx = sys.platform.startswith('darwin')
is_linux = sys.platform.startswith('linux')

if(is_windows): # suppress warnings on windows
    import warnings
    warnings.filterwarnings('ignore', \
        message='Valid PROJ data directory not found')
# end if

import numpy as np
from obspy import UTCDateTime
import click
from collections import defaultdict
import matplotlib
from seismic.ASDFdatabase.utils import MseedIndex
from seismic.ASDFdatabase.utils import MAX_DATE, MIN_DATE
import psutil
import multiprocess
from multiprocess import Manager, freeze_support
from seismic.ASDFdatabase.analytics.sun import day_night_coverage, DAY_NIGHT_SECONDS
from seismic.ASDFdatabase.analytics.station_analytics import StationAnalytics
from seismic.ASDFdatabase.analytics.utils import ProgressTracker, get_response

if(is_windows | is_osx): matplotlib.use('TKAgg')
else: matplotlib.use('Agg')

def get_cpu_count(nproc):
    result = nproc
    cpu_count = psutil.cpu_count(logical=False)
    if (nproc > cpu_count or nproc == -1):
        result = cpu_count
    # end if
    return result
# end func

def select_channel(mi:MseedIndex, sd:UTCDateTime, ed:UTCDateTime)->dict(list([])):
    meta_list = mi.get_stations(sd, ed)
    cov_dict = mi.get_channel_coverages()
    channel_meta_dict = defaultdict(list)
    channel_cov_dict = defaultdict(float)

    for meta in meta_list:
        cha = meta[3]
        nslc = '.'.join(meta[:4])
        channel_meta_dict[cha].append(meta)
        channel_cov_dict[cha] += cov_dict[nslc]
    # end for

    if(len(meta_list) == 0):
        raise RuntimeError('No stations found between {} -- {}. Aborting..'.format(sd, ed))
    elif(len(meta_list) == 1):
        return channel_meta_dict
    else:
        answer_list = [i + 1 for i in np.arange(len(channel_meta_dict))]
        selection = []
        print('\n############################')
        print('# Multiple channels found: #')
        print('############################\n')
        sorted_channels = sorted(channel_meta_dict.keys())
        for i, cha in enumerate(sorted_channels):
            print('{}) {} ({:2.3f} days)'.format(answer_list[i],
                                                cha, channel_cov_dict[cha]/DAY_NIGHT_SECONDS))
        # end for

        while (not (np.all([s in answer_list for s in selection]) if len(selection) else False)):
            try:
                selection = list(input("""\n*** Please select desired channel *** \n(select multiple channels with entries separated by commas) \n(-1 to quit): """).split(','))
                selection = [int(s.strip()) for s in selection]

                if (-1 in selection): break
                if (not np.all([s in answer_list for s in selection])): print('Unexpected input(s) found.')
            except:
                print('Invalid selection. ')
                selection = []
            # end try
        # wend

        if(-1 in selection): exit(0)

        result = dict()
        for s in selection:
            channel = sorted_channels[s-1]
            channel_meta_list = channel_meta_dict[channel]

            print('\nData found under the following Network-Station-Location list will be analysed for channel {}:\n'.format(channel))
            for channel_meta in channel_meta_list:
                nslc = '.'.join(channel_meta)
                nsl = '.'.join(channel_meta[:3])
                print('{} ({:2.3f} days)'.format(nsl, cov_dict[nslc]/DAY_NIGHT_SECONDS))
            # end for

            result[channel] = channel_meta_list
        # end for

        return result
    # end if
# end func

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('mseed-folder', required=True,
                type=click.Path(exists=True))
@click.argument('mseed-pattern', required=True,
                type=str)
@click.argument('instrument-response', required=True,
                type=click.Path(exists=True))
@click.argument('output-folder', required=True,
                type=click.Path(exists=True, file_okay=False, dir_okay=True))
@click.option('--start-date', type=str, default=None, show_default=True,
              help="Start date in UTC format for processing data")
@click.option('--end-date', type=str, default=None, show_default=True,
              help="End date in UTC format for processing data")
@click.option('--nproc', type=int, default=-1, show_default=True,
              help="Number of parallel processes to use. Default is to use all available cores.")
def process_mseed(mseed_folder, mseed_pattern, instrument_response,
                  output_folder, start_date, end_date, nproc):
    """
    MSEED_FOLDER: Path to folder containing mseed files\n
    MSEED_PATTERN: File pattern to be used to capture files pertaining to specific channels.
                   Note that pattern must be specified within quotes. \n
    INSTRUMENT_RESPONSE: Path to instrument response in .resp format\n
    OUTPUT_FOLDER: Path to output folder\n
    """
    try:
        start_date = UTCDateTime(start_date) if start_date else None
        end_date   = UTCDateTime(end_date) if end_date else None
    except Exception as e:
        print(str(e))
        raise RuntimeError('Invalid start- or end-dates')
    # end try

    if(start_date and end_date and ((end_date.date - start_date.date).days<=0)):
        raise RuntimeError('Invalid start- and end-dates. Aborting..')
    # end if

    # instantiate MseedIndex
    print('Inspecting mseed files in {}..'.format(mseed_folder))
    mseed_index = MseedIndex(mseed_folder, mseed_pattern)

    # define bound getter methods
    def get_waveforms_func(net, sta, loc, cha, st, et):
        return mseed_index.get_waveforms(net, sta, loc, cha, st, et)
    # end func

    def get_time_range_func(net, sta, loc, cha):
        return mseed_index.get_time_range(net, sta, loc, cha)
    # end func

    sd = MIN_DATE if start_date is None else start_date
    ed = MAX_DATE if end_date is None else end_date

    nproc = get_cpu_count(nproc)

    print('Loading response..')
    resp = get_response(instrument_response)
    if(resp is not None):
        print('Found response: {}'.format(resp))
        sys.stdout.flush()
    else:
        raise(RuntimeError('No instrument response found. Aborting..'))
    # end if

    # Recordings often feature garbled network, station and location codes.
    # We only consider unique channel codes as meaningful and as such, get a list of
    # net, sta, loc, cha combinations that share the same channel code. We also
    # account for the possibility of each combination of net, sta, loc, cha to have
    # recordings under different sampling rates
    meta_dict = select_channel(mseed_index, sd, ed)
    sr_dict = mseed_index.get_channel_sampling_rates()

    # instantiate progress tracker
    manager = Manager()
    prog_tracker = ProgressTracker(manager)
    sa = StationAnalytics(get_time_range_func, get_waveforms_func,
                          resp, output_folder, prog_tracker, sd, ed, nproc)

    for cha, meta_list in meta_dict.items():

        print("""
        \n*** Processing data for channel {} ***\n
        """.format(cha))
        for meta in meta_list:
            net, sta, loc, cha = meta

            nslc = '.'.join(meta)
            sampling_rate = sr_dict[nslc]

            sa.analyse_data(net, sta, loc, cha, sampling_rate)
        # end for

        ofn = os.path.join(output_folder, '{}.pdf'.format(cha))
        sa.process_results(cha, ofn)
        print('\nDone..')
    # end for

# end func

if __name__ == "__main__":
    # add support for process-based multiprocessing for a Windows .exe
    if(is_windows): freeze_support()
    if(is_osx): multiprocess.set_start_method('spawn')

    process_mseed()
# end func
