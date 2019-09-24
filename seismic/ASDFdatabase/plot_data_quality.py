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
import glob, os, sys
from collections import defaultdict
import numpy as np
from obspy import Stream, Trace, UTCDateTime
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet

import click, logging
import matplotlib.pyplot as plt
import glob
from mpl_toolkits.basemap import Basemap
from descartes import PolygonPatch
from shapely.geometry import Polygon
from matplotlib.patches import Ellipse
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.dates import DateFormatter, \
    AutoDateLocator, \
    HourLocator, \
    MinuteLocator, \
    epoch2num
from matplotlib.ticker import ScalarFormatter, FuncFormatter

import gc

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

def process_data(fds, stations):
    results = []

    day = 3600*24
    for s in stations:
        st, et = fds.get_global_time_range(s[0], s[1])

        # align st to day
        st = UTCDateTime(year=st.year, month=st.month, day=st.day)

        ct = st
        times = []
        means = []
        while (ct < et):
            times.append(ct)
            stream = fds.get_waveforms(s[0], s[1], s[2], s[3],
                                      ct, ct + day,
                                      trace_count_threshold=200)
            if (len(stream)):
                try:
                    tr = stream.merge()
                    means.append(np.mean(tr[0].data))
                except:
                    means.append(np.nan)
            else:
                means.append(np.nan)
            # end if
            ct += day
            #break
        # end while
        results.append([times, means])
    # end for

    return results
# end func

def plot_results(stations, results, output_basename):
    pdf = PdfPages('%s.pdf'%(output_basename))

    def taper(tr, taperlen):
        tr[0:taperlen] *= 0.5 * (1 + np.cos(np.linspace(-np.pi, 0, taperlen)))
        tr[-taperlen:] *= 0.5 * (1 + np.cos(np.linspace(0, np.pi, taperlen)))
        return tr
    # end func

    def drawBBox(minLon, minLat, maxLon, maxLat, bm, **kwargs):
        bblons = np.array([minLon, maxLon, maxLon, minLon, minLon])
        bblats = np.array([minLat, minLat, maxLat, maxLat, minLat])

        x, y = bm(bblons, bblats)
        xy = zip(x, y)
        poly = Polygon(xy)
        bm.ax.add_patch(PolygonPatch(poly, **kwargs))
    # end func

    fig = plt.figure(figsize=(10, 10))

    minLon = 1e32
    maxLon = -1e32
    minLat = 1e32
    maxLat = -1e32
    for s in stations:
        lon, lat = s[4], s[5]

        if (minLon > lon): minLon = lon
        if (maxLon < lon): maxLon = lon
        if (minLat > lat): minLat = lat
        if (maxLat < lat): maxLat = lat
    # end for

    minLon -= 1
    maxLon += 1
    minLat -= 1
    maxLat +=1

    m = Basemap(width=800000, height=800000, projection='lcc',
                resolution='l', lat_1=minLat, lat_2=maxLat,
                lat_0=(minLat + maxLat) / 2., lon_0=(minLon + maxLon) / 2.)
    # draw coastlines.
    m.drawcoastlines()

    # draw grid
    parallels = np.linspace(np.floor(minLat) - 5, np.ceil(maxLat) + 5, \
                            int((np.ceil(maxLat) + 5) - (np.floor(minLat) - 5)) + 1)
    m.drawparallels(parallels, labels=[True, True, False, False])
    meridians = np.linspace(np.floor(minLon) - 5, np.ceil(maxLon) + 5, \
                            int((np.ceil(maxLon) + 5) - (np.floor(minLon) - 5)) + 1)
    m.drawmeridians(meridians, labels=[False, False, True, True])

    # plot stations
    for s in stations:
        lon, lat = s[4], s[5]

        px, py = m(lon, lat)
        pxl, pyl = m(lon, lat - 0.1)
        m.scatter(px, py, 50, marker='v', c='g', edgecolor='none')
        plt.annotate(s[1], xy=(px, py), fontsize=5)
    # end for

    insetAx = fig.add_axes([0.75, 0.75, 0.125, 0.125])
    mInset = Basemap(resolution='c',  # c, l, i, h, f or None
                     ax=insetAx,
                     projection='merc',
                     lat_0=-20, lon_0=132,
                     llcrnrlon=110, llcrnrlat=-40, urcrnrlon=155, urcrnrlat=-10)
    # mInset.drawcoastlines()
    mInset.fillcontinents(color='lightgray')
    mInset.drawstates(color="grey")

    drawBBox(minLon, minLat, maxLon, maxLat, mInset, fill='True', facecolor='k')

    fig.axes[0].set_title("Network Name: %s"%s[0], fontsize=20, y=1.05)
    fig.axes[0].legend()
    plt.legend()

    pdf.savefig()

    # Plot results
    assert len(stations) == len(results)
    groupIndices = defaultdict(list)
    for i in np.arange(len(results)):
        groupIndices['%s.%s'%(stations[i][0], stations[i][1])].append(i)
    # end for

    for k,v in groupIndices.iteritems():
        fig, axes = plt.subplots(len(v), sharex=True)
        fig.set_size_inches(20, 15)

        for i, index in enumerate(v):
            x, means = results[index]
            x = [a.matplotlib_date for a in x]
            d = np.array(means)

            if(len(d)): d[0] = np.nanmedian(d)

            # print k
            # for val in d: print val
            #dnorm = 2 * ((d - np.nanmin(d)) / (np.nanmax(d) - np.nanmin(d))) - 1
            dnorm = d
            dnormmin = np.nanmin(dnorm)
            dnormmax = np.nanmax(dnorm)

            axes[i].scatter(x, dnorm, marker='.')
            axes[i].plot(x, dnorm, c='k', label='Mean %s over 24 Hrs '
                                                '\n(normalized between [-1,1])\n'
                                                'Gaps indicate no-data' % stations[index][3], lw=2)

            try:
                axes[i].fill_between(x, dnormmax * np.int_(d == 0), dnormmin * np.int_(d == 0),
                                     where=dnormmax * np.int_(d == 0) - dnormmin * np.int_(d == 0) > 0,
                                     color='r', alpha=0.5, label='All 0 Samples')

                axes[i].fill_between(x, dnormmax * np.int_(np.isnan(d)), dnormmin * np.int_(np.isnan(d)),
                                     where=dnormmax * np.int_(np.isnan(d)) - dnormmin * np.int_(np.isnan(d)) > 1,
                                     color='b', alpha=0.5, label='No Data')

                axes[i].xaxis.set_major_locator(AutoDateLocator())
                axes[i].xaxis.set_major_formatter(DateFormatter('%Y-%m-%d'))

                for tick in axes[i].get_xticklabels():
                    tick.set_rotation(45)
                # end for
                #axes[i].set_ylim(-1.5, 1.5)
                axes[i].legend()
                axes[i].tick_params(axis='both', labelsize=14)
            except:
                pass
            # end try
        # end for
        axes[-1].set_xlabel('Days', fontsize=14)
        plt.suptitle('%s Data Availability'%(k), fontsize=20)
        pdf.savefig()
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
    STA: Station name (* for all stations)\n
    CHA: Channel name (* for all channels) \n
    OUTPUT_BASENAME: Basename of output file
    """

    start_time = UTCDateTime(start_time)
    end_time = UTCDateTime(end_time)
    if (sta == '*'): sta = None
    if (cha == '*'): cha = None

    comm = MPI.COMM_WORLD
    nproc = comm.Get_size()
    rank = comm.Get_rank()
    proc_workload = None

    l = setup_logger(name=output_basename, log_file='%s.log'%output_basename)
    fds = FederatedASDFDataSet(asdf_source, logger=None)

    stations = []
    if(rank == 0):
        stations = fds.get_stations(start_time, end_time, network=net, station=sta, channel=cha)

        stations = split_list(stations, nproc)
    # end if

    stations = comm.bcast(stations, root=0)
    results = process_data(fds, stations[rank])

    results = comm.gather(results, root=0)
    if (rank == 0):
        results = [item for sublist in results for item in sublist] # flatten sublists for each proc
        stations = [item for sublist in stations for item in sublist]  # flatten sublists for each proc
        plot_results(stations, results, output_basename)
    # end if
# end func

if (__name__=='__main__'):
    process()
# end if