#!/usr/env python
"""
Description:
    Reads waveform data from a FederatedASDFDataSet and generates data-quality plots into a multi-page PDF file

    For each station, the program read-in all the waveform data over the period specified (usually over a  year),
    and then perform statistics analysis over the waveform data for subsequent plotting.

    The first page of the PDF shows the stations on a basemap with marker colored according to the number of usable days
    ("good" data available in the day, not gap, not all zeros) And  the color is determined by c=mapper.to_rgba(days):
    Relatively speaking, black/blue/cold-ish colors indicate higher percentage of good data in the station.
    And red/yellow/warm-ish colors indicate higher percentage of problematic data.

    The following pages of the PDF will be plotting the daily-averaged seismic wave data
    with gaps and all-zeros days are shaded.

CreationDate:   19/09/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     19/09/2019   RH
    LastUpdate:     27/05/2020   FZ     Refactoring and docs run examples etc.

Todo:
    The script currently have a low pylint score 4.3/10.
    Need to refactor to  comply with Python standard pep-8: add required docstrings.
    use smaller function to modularize better.
    Consider to generate an additional page to executive summarise info for user.
"""

import gc
import logging
from collections import defaultdict

import click
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.dates import DateFormatter, AutoDateLocator
from mpi4py import MPI
import cartopy.crs as ccrs
from obspy import UTCDateTime
from ordered_set import OrderedSet as set
from tqdm import tqdm

from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.misc import split_list

logging.basicConfig()

def setup_logger(name, log_file, level=logging.INFO):
    """
    Function to setup a logger; adapted from stackoverflow
    """
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler = logging.FileHandler(log_file, mode='w')
    handler.setFormatter(formatter)

    logger = logging.getLogger(name + log_file)
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger


# end func


def process_data(rank, fds, stations, start_time, end_time):
    """

    :param rank: processor rank in parallel runs
    :param fds: FederatedASDFDataSer instance
    :param stations: list containing tuples (net, sta, loc, cha, lon, lat)
    :param start_time: start-time in UTCDateTime format
    :param end_time: end-time in UTCDateTime format
    :return: a list containing UTCDateTimes marking the start of each day
             and their corresponding means in each row
    """
    results = []

    day = 3600 * 24
    # Discard aux channels not in the following list
    chaSet = {'BH1', 'BH2', 'BHE', 'BHN', 'BHZ', 'BNZ', 'EHE', 'EHN', 'EHZ', 'HHE', 'HHN',
              'HHZ', 'LHE', 'LHN', 'LHZ', 'SHE', 'SHN', 'SHZ'}
    for s in tqdm(stations, desc='Rank: %d' % (rank)):
        # discard aux channels
        if (s[3] not in chaSet):
            results.append(([], []))
            continue
        # end if

        st, et = fds.get_global_time_range(s[0], s[1])

        if (start_time > st): st = start_time
        if (end_time < et): et = end_time

        # align st to day
        st = UTCDateTime(year=st.year, month=st.month, day=st.day)

        ct = st
        times = []
        means = []
        while (ct < et):
            times.append(ct)

            debug = False
            if not debug:
                stream = fds.get_waveforms(s[0], s[1], s[2], s[3],
                                           ct, ct + day,
                                           trace_count_threshold=200)
                if len(stream):
                    try:
                        tr = stream.merge()
                        means.append(np.mean(tr[0].data))
                    except:
                        # empty or streams that cannot be merged are marked by nans
                        means.append(np.nan)
                else:
                    means.append(np.nan)
                # end if
            else:
                # Used for debug purposes only
                if np.int_(np.round(np.random.random())):
                    means.append(np.random.random())
                else:
                    means.append(np.nan)
                # end if
            # end if

            ct += day
        # end while

        results.append([times, means])
    # end for

    return results


# end func


def plot_results(stations, results, output_basename):
    # collate indices for each channel for each station
    assert len(stations) == len(results)
    groupIndices = defaultdict(list)
    for i in np.arange(len(results)):
        groupIndices['%s.%s' % (stations[i][0], stations[i][1])].append(i)
    # end for

    # gather number of days of usable data for each station
    usableStationDays = defaultdict(int)
    maxUsableDays = -1e32
    minUsableDays = 1e32
    for k, v in groupIndices.items():
        for i, index in enumerate(v):
            x, means = results[index]

            means = np.array(means)
            days = np.sum(~np.isnan(means) & np.bool_(means != 0))
            if usableStationDays[k] < days:
                usableStationDays[k] = days
                maxUsableDays = max(maxUsableDays, days)
                minUsableDays = min(minUsableDays, days)
            # end if
        # end for
    # end for

    # Plot station map
    pdf = PdfPages('%s.pdf' % output_basename)

    fig = plt.figure(figsize=(20, 30))
    ax1 = fig.add_axes([0.05, 0.05, 0.9, 0.7], projection=ccrs.PlateCarree())

    minLon = 1e32
    maxLon = -1e32
    minLat = 1e32
    maxLat = -1e32
    for s in stations:
        lon, lat = s[4], s[5]
        if lon < 0:
            lon += 360
        # end if
        minLon = min(minLon, lon)
        maxLon = max(maxLon, lon)
        minLat = min(minLat, lat)
        maxLat = max(maxLat, lat)
    # end for

    minLon -= 1
    maxLon += 1
    minLat -= 1
    maxLat += 1

    ax1.set_extent([minLon, maxLon, minLat, maxLat], crs=ccrs.PlateCarree())
    # draw coastlines.
    ax1.coastlines('50m')

    # draw grid
    ax1.gridlines(draw_labels=True,
                  linewidth=1, color='gray',
                  alpha=0.5, linestyle='--')

    # plot stations
    norm = matplotlib.colors.Normalize(vmin=minUsableDays, vmax=maxUsableDays, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.jet_r)
    plotted = set()
    for s in stations:
        if s[1] in plotted:
            continue
        else:
            plotted.add(s[1])
        # end if

        lon, lat = s[4], s[5]

        days = usableStationDays['%s.%s' % (s[0], s[1])]
        ax1.scatter(lon, lat, s=400, transform=ccrs.PlateCarree(), marker='v', c=mapper.to_rgba(days),
                    edgecolor='none', label='%s: %d' % (s[1], days))
        ax1.annotate(s[1], xy=(lon + 0.05, lat + 0.05),
                     xycoords=ccrs.PlateCarree()._as_mpl_transform(ax1),
                     fontsize=22)
    # end for

    ax1.set_title("Network Name: %s" % s[0], fontsize=30, y=1.05)
    ax1.legend(prop={'size': 16}, loc=(0.2, 1.3),
               ncol=5, fancybox=True, title='No. of Usable Days',
               title_fontsize=16)
    #fig.colorbar(mapper, orientation='horizontal', label='Days')
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    # Plot results
    for k, v in groupIndices.items():
        axesCount = 0
        for i in v:
            assert (k == '%s.%s' % (stations[i][0], stations[i][1]))
            # only need axes for non-null results
            a, b = results[i]
            if len(a) and len(b):
                axesCount += 1
            # end if
        # end for
        fig, axes = plt.subplots(axesCount, sharex=True)
        fig.set_size_inches(20, 15)

        axes = np.atleast_1d(axes)

        if len(axes):
            axesIdx = 0
            for i, index in enumerate(v):
                try:
                    x, means = results[index]

                    if len(x) and len(means):
                        x = [a.matplotlib_date for a in x]
                        d = np.array(means)

                        if len(d):
                            d[0] = np.nanmedian(d)
                        # end if

                        dnorm = d
                        dnormmin = np.nanmin(dnorm)
                        dnormmax = np.nanmax(dnorm)

                        axes[axesIdx].scatter(x, dnorm, marker='.', s=20)
                        axes[axesIdx].plot(x, dnorm, c='k', label='24 hr mean\n'
                                                                  'Gaps indicate no-data', lw=2, alpha=0.7)
                        axes[axesIdx].grid(axis='x', linestyle=':', alpha=0.3)

                        axes[axesIdx].fill_between(x, dnormmax * np.int_(d == 0), dnormmin * np.int_(d == 0),
                                                   where=dnormmax * np.int_(d == 0) - dnormmin * np.int_(d == 0) > 0,
                                                   color='r', alpha=0.5, label='All 0 Samples')

                        axes[axesIdx].fill_between(x, dnormmax * np.int_(np.isnan(d)), dnormmin * np.int_(np.isnan(d)),
                                                   where=dnormmax * np.int_(np.isnan(d)) - dnormmin * np.int_(
                                                       np.isnan(d)) > 1,
                                                   color='b', alpha=0.5, label='No Data')

                        axes[axesIdx].xaxis.set_major_locator(AutoDateLocator())
                        axes[axesIdx].xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

                        for tick in axes[axesIdx].get_xticklabels():
                            tick.set_rotation(45)
                        # end for
                        axes[axesIdx].legend(loc='upper right', prop={'size': 12})
                        axes[axesIdx].tick_params(axis='both', labelsize=16)
                        stn = stations[index]
                        axes[axesIdx].set_title('Channel %s.%s' % (stn[2], stn[3]),
                                                fontsize=18, y=0.95, va='top')
                        axes[axesIdx].set_xlim(xmin=min(x), xmax=max(x))
                        axes[axesIdx].set_ylim(ymin=dnormmin, ymax=dnormmax)
                        axes[axesIdx].set_ylabel('Ampl.', fontsize=16)

                        axesIdx += 1
                    # end if
                except:
                    # plotting fails when each axes contain <2 values; just move on in those instances
                    logging.warning('Plotting failed on station %s' % k)
                # end try
            # end for
            axes[-1].set_xlabel('Days', fontsize=16)
        # end if

        plt.suptitle('%s Data Availability (~%d days)' % (k, usableStationDays[k]),
                     y=0.96, fontsize=20)
        pdf.savefig()
        plt.close()
        gc.collect()

    # end for

    pdf.close()
# end func


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('asdf-source', required=True,
                type=click.Path(exists=True))
@click.argument('start-time', required=True,
                type=str)
@click.argument('end-time', required=True,
                type=str)
@click.argument('net', required=True,
                type=str)
@click.argument('sta', required=True,
                type=str)
@click.argument('cha', required=True,
                type=str)
@click.argument('output-basename', required=True,
                type=str)
def process(asdf_source, start_time, end_time, net, sta, cha, output_basename):
    """
    ASDF_SOURCE: Text file containing a list of paths to ASDF files\n
    START_TIME: Start time in UTCDateTime format\n
    END_TIME: End time in UTCDateTime format\n
    NET: Network name\n
    STA: Station name ('*' for all stations; note that * must be in quotation marks)\n
    CHA: Channel name ('*' for all channels; note that * must be in quotation marks) \n
    OUTPUT_BASENAME: Basename of output file

    Example usage:
    mpirun -np 112 python plot_data_quality.py asdf_files.txt 1980:01:01 2020:01:01 OA '*' '*' data_quality.oa
    """

    start_time = UTCDateTime(start_time)
    end_time = UTCDateTime(end_time)
    if (sta == '*'): sta = None
    if (cha == '*'): cha = None

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()

    l = setup_logger(name=output_basename, log_file='%s.log' % output_basename)
    fds = FederatedASDFDataSet(asdf_source, logger=l)

    stations = []
    if rank == 0:
        stations = fds.get_stations(start_time, end_time, network=net, station=sta, channel=cha)

        stations = split_list(sorted(stations), nproc)
    # end if

    stations = comm.bcast(stations, root=0)
    results = process_data(rank, fds, sorted(stations[rank]), start_time, end_time)

    results = comm.gather(results, root=0)
    if rank == 0:
        results = [item for sublist in results for item in sublist]  # flatten sublists for each proc
        stations = [item for sublist in stations for item in sublist]  # flatten sublists for each proc
        plot_results(stations, results, output_basename)
    # end if


# end func


if (__name__ == '__main__'):
    """
    Example usage:
    mpirun -np 112 python plot_data_quality.py asdf_files.txt 1980:01:01 2020:01:01 OA '*' '*' data_quality.oa
    
    How to run without MPI in background 

    Example commandline:
    (hiperseispy37) fxz547@vdi-n16 /g/data/ha3/fxz547/Githubz/hiperseis (develop)
    $ nohup python seismic/ASDFdatabase/plot_data_quality.py /g/data/ha3/GASeisDataArchive/DevSpace/test_asdf_files.txt 2019:01:01 2019:03:01 'AU' '*' '*' test_plot_data_quality_AU &

    """

    process()
