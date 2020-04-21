#!/usr/env python
"""
Description:
    Reads waveform data from a FederatedASDFDataSet and generates data-quality
    plots
References:

CreationDate:   19/09/19
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     19/09/19   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from mpi4py import MPI
from collections import defaultdict
import logging
import gc

from ordered_set import OrderedSet as set
import numpy as np
from obspy import UTCDateTime
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

import click
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.dates import DateFormatter, AutoDateLocator
import matplotlib
import matplotlib.cm as cm
from tqdm import tqdm


logging.basicConfig()


def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
# end func


def setup_logger(name, log_file, level=logging.INFO):
    """
    Function to setup a logger; adapted from stackoverflow
    """
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler = logging.FileHandler(log_file, mode='w')
    handler.setFormatter(formatter)

    logger = logging.getLogger(name+log_file)
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

    day = 3600*24
    # Discard aux channels not in the following list
    chaSet = {'BH1','BH2','BHE','BHN','BHZ','BNZ','EHE','EHN','EHZ','HHE','HHN',
              'HHZ','LHE','LHN','LHZ','SHE','SHN','SHZ'}
    for s in tqdm(stations, desc='Rank: %d'%(rank)):
        # discard aux channels
        if (s[3] not in chaSet):
            results.append(([], []))
            continue
        # end if

        st, et = fds.get_global_time_range(s[0], s[1])

        if(start_time > st): st = start_time
        if(end_time < et): et = end_time

        # align st to day
        st = UTCDateTime(year=st.year, month=st.month, day=st.day)

        ct = st
        times = []
        means = []
        while (ct < et):
            times.append(ct)

            if(True):
                stream = fds.get_waveforms(s[0], s[1], s[2], s[3],
                                          ct, ct + day,
                                          trace_count_threshold=200)
                if (len(stream)):
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
                if(np.int_(np.round(np.random.random()))):
                    means.append(np.random.random())
                else:
                    means.append(np.nan)
                # end if
            # end if

            ct += day
            #break
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
        groupIndices['%s.%s'%(stations[i][0], stations[i][1])].append(i)
    # end for

    # gather number of days of usable data for each station
    usableStationDays = defaultdict(int)
    maxUsableDays = -1e32
    minUsableDays = 1e32
    for k,v in groupIndices.items():
        for i, index in enumerate(v):
            x, means = results[index]

            means = np.array(means)
            days = np.sum(~np.isnan(means) & np.bool_(means!=0))
            if(usableStationDays[k] < days):
                usableStationDays[k] = days

                if(maxUsableDays < days): maxUsableDays = days
                if(minUsableDays > days): minUsableDays = days
            # end if
        # end for
    # end for

    # Plot station map
    pdf = PdfPages('%s.pdf'%(output_basename))

    fig = plt.figure(figsize=(20, 30))
    ax1 = fig.add_axes([0.05, 0.05, 0.9, 0.7])
    ax2 = fig.add_axes([0.05, 0.7, 0.9, 0.3])
    ax2.set_visible(False)

    minLon = 1e32
    maxLon = -1e32
    minLat = 1e32
    maxLat = -1e32
    for s in stations:
        lon, lat = s[4], s[5]

        if(lon<0): lon += 360

        if (minLon > lon): minLon = lon
        if (maxLon < lon): maxLon = lon
        if (minLat > lat): minLat = lat
        if (maxLat < lat): maxLat = lat
    # end for

    minLon -= 1
    maxLon += 1
    minLat -= 1
    maxLat += 1

    m = Basemap(ax=ax1, projection='merc',
                resolution='i', llcrnrlat=minLat, urcrnrlat=maxLat,
                llcrnrlon=minLon, urcrnrlon=maxLon,
                lat_0=(minLat + maxLat) / 2., lon_0=(minLon + maxLon) / 2.)
    # draw coastlines.
    m.drawcoastlines()

    # draw grid
    parallels = np.linspace(np.around(minLat / 5) * 5 - 5, np.around(maxLat / 5) * 5 + 5, 6)
    m.drawparallels(parallels, labels=[True, True, False, False])
    meridians = np.linspace(np.around(minLon / 5) * 5 - 5, np.around(maxLon / 5) * 5 + 5, 6)
    m.drawmeridians(meridians, labels=[False, False, True, True])

    # plot stations
    norm = matplotlib.colors.Normalize(vmin=minUsableDays, vmax=maxUsableDays, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.jet_r)
    plotted = set()
    for s in stations:
        if(s[1] in plotted): continue
        else: plotted.add(s[1])

        lon, lat = s[4], s[5]

        px, py = m(lon, lat)
        pxl, pyl = m(lon, lat - 0.1)
        days = usableStationDays['%s.%s' % (s[0], s[1])]
        m.scatter(px, py, 50, marker='v',
                  c=mapper.to_rgba(days),
                  edgecolor='none', label='%s: %d'%(s[1], days))
        ax1.annotate(s[1], xy=(px, py), fontsize=5)
    # end for

    fig.axes[0].set_title("Network Name: %s"%s[0], fontsize=20, y=1.05)
    fig.axes[0].legend(prop={'size':5}, bbox_to_anchor=(0.2, 1.3),
                       ncol=5, fancybox=True, title='No. of Usable Days')

    pdf.savefig()

    # Plot results
    for k,v in groupIndices.items():
        axesCount = 0
        for i in v:

            assert (k == '%s.%s'%(stations[i][0], stations[i][1]))

            # only need axes for non-null results
            a, b = results[i]
            if(len(a) and len(b)): axesCount += 1
        # end for
        fig, axes = plt.subplots(axesCount, sharex=True)
        fig.set_size_inches(20, 15)

        axes = np.atleast_1d(axes)

        if(len(axes)):
            axesIdx = 0
            for i, index in enumerate(v):
                try:
                    x, means = results[index]

                    if(len(x) and len(means)):
                        x = [a.matplotlib_date for a in x]
                        d = np.array(means)

                        if(len(d)): d[0] = np.nanmedian(d)

                        # print k
                        # for val in d: print val
                        #dnorm = 2 * ((d - np.nanmin(d)) / (np.nanmax(d) - np.nanmin(d))) - 1
                        dnorm = d
                        dnormmin = np.nanmin(dnorm)
                        dnormmax = np.nanmax(dnorm)

                        axes[axesIdx].scatter(x, dnorm, marker='.')
                        axes[axesIdx].plot(x, dnorm, c='k', label='Mean %s over 24 Hrs\n'
                                                            'Gaps indicate no-data' % stations[index][3], lw=2)

                        axes[axesIdx].fill_between(x, dnormmax * np.int_(d == 0), dnormmin * np.int_(d == 0),
                                             where=dnormmax * np.int_(d == 0) - dnormmin * np.int_(d == 0) > 0,
                                             color='r', alpha=0.5, label='All 0 Samples')

                        axes[axesIdx].fill_between(x, dnormmax * np.int_(np.isnan(d)), dnormmin * np.int_(np.isnan(d)),
                                             where=dnormmax * np.int_(np.isnan(d)) - dnormmin * np.int_(np.isnan(d)) > 1,
                                             color='b', alpha=0.5, label='No Data')

                        axes[axesIdx].xaxis.set_major_locator(AutoDateLocator())
                        axes[axesIdx].xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

                        for tick in axes[axesIdx].get_xticklabels():
                            tick.set_rotation(45)
                        # end for
                        #axes[axesIdx].set_ylim(-1.5, 1.5)
                        axes[axesIdx].legend(loc='upper right', prop={'size':5})
                        axes[axesIdx].tick_params(axis='both', labelsize=14)

                        axesIdx += 1
                    # end if
                except:
                    # plotting fails when each axes contain <2 values; just move on in those instances
                    pass
                # end try
            # end for
            axes[-1].set_xlabel('Days', fontsize=14)
        # end if

        plt.suptitle('%s Data Availability (~%d days)'%(k, usableStationDays[k]), fontsize=20)
        pdf.savefig()
        gc.collect()

        #break
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

    l = setup_logger(name=output_basename, log_file='%s.log'%output_basename)
    fds = FederatedASDFDataSet(asdf_source, logger=l)

    stations = []
    if(rank == 0):
        stations = fds.get_stations(start_time, end_time, network=net, station=sta, channel=cha)

        stations = split_list(sorted(stations), nproc)
    # end if

    stations = comm.bcast(stations, root=0)
    results = process_data(rank, fds, sorted(stations[rank]), start_time, end_time)

    results = comm.gather(results, root=0)
    if (rank == 0):
        results = [item for sublist in results for item in sublist] # flatten sublists for each proc
        stations = [item for sublist in stations for item in sublist]  # flatten sublists for each proc
        plot_results(stations, results, output_basename)
    # end if
# end func


if (__name__=='__main__'):
    '''
    Example usage:
    mpirun -np 112 python plot_data_quality.py asdf_files.txt 1980:01:01 2020:01:01 OA '*' '*' data_quality.oa
    '''

    process()
# end if
