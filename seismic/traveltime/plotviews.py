from __future__ import print_function, absolute_import

import logging
import os
from collections import namedtuple
from math import asin

import matplotlib
import numpy as np
import pandas as pd

# Force matplotlib to not use any Xwindows backend
matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.lines import Line2D
from mpl_toolkits.basemap import Basemap
import click
from seismic import pslog

DPI = asin(1.0) / 90.0
R2D = 90. / asin(1.)
FLOAT_FORMAT = '%.4f'

log = logging.getLogger(__name__)

SOURCE_LATITUDE = 'source_latitude'
SOURCE_LONGITUDE = 'source_longitude'
STATION_LATITUDE = 'station_latitude'
STATION_LONGITUDE = 'station_longitude'
STATION_CODE = 'station_code'
FREQUENCY = 'no_of_summary_rays'

column_names = ['source_block', 'station_block',
                'residual', 'event_number',
                SOURCE_LONGITUDE, SOURCE_LATITUDE,
                'source_depth', STATION_LONGITUDE, STATION_LATITUDE,
                'observed_tt', 'locations2degrees', STATION_CODE, 'SNR', 'P_or_S']

# since we have Basemap in the virtualenv, let's just use that :)
ANZ = Basemap(llcrnrlon=100.0, llcrnrlat=-50.0, urcrnrlon=190.0, urcrnrlat=0.0) # original small region.
ANZ = Basemap(llcrnrlon=50.0, llcrnrlat=-70.0, urcrnrlon=200.0, urcrnrlat=30.0)

Region = namedtuple('Region', 'upperlat, bottomlat, leftlon, rightlon')


@click.group()
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cli(verbosity):
    pslog.configure(verbosity)


@cli.command()
@click.argument('arrivals_file', type=click.File(mode='r'))
@click.argument('region', type=str,
                metavar="str 'upperlat, bottomlat, leftlon, rightlon'")
def plot(arrivals_file, region):  # pragma: no cover
    """
    This command will output a `sources_in_region.png`, which will show all the
    sources inside the `region` specified by the region string which can be
    specified like '0 -50.0 100 190'. It will also output
    a`stations_in_region.png` showing the stations where arrivals were
    recorded within the `region`.
    The `cluster plot` command further outputs a
    `sources_and_stations_in_region.png` which should all sources and
    stations in the same plot that is within `region`.

    Output file from each other `cluster` `gather`, `sort`, `match` and `zone`
    can be visualised using the `cluster plot` command.
    """
    region = [float(s) for s in region.split()]
    reg = Region(*region)

    arrivals = pd.read_csv(arrivals_file, header=None, names=column_names,
                           sep=' ')
    arr_file_base = os.path.splitext(arrivals_file.name)[0]
    source = _source_or_stations_in_region(
        arrivals, reg, SOURCE_LATITUDE, SOURCE_LONGITUDE,
        'sources_in_region_{}.png'.format(arr_file_base))

    station = _source_or_stations_in_region(
        arrivals, reg, STATION_LATITUDE, STATION_LONGITUDE,
        'stations_in_region_{}.png'.format(arr_file_base))

    # sources and stations both in region
    sources_and_stations = arrivals[source & station]

    fig = plt.figure()

    _plot_on_map(sources_and_stations,
                 SOURCE_LONGITUDE, SOURCE_LATITUDE,
                 marker='*', color='r')
    _plot_on_map(sources_and_stations,
                 STATION_LONGITUDE, STATION_LATITUDE,
                 marker='^', color='b')

    plt.title('Sources and stations in \n region {}'.format(region))
    # plt.xlabel('Longitude')
    # plt.ylabel('Latitude')
    fig.savefig('sources_and_stations_in_region_{}.png'.format(arr_file_base))

    # rays originating and terminating in region
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for i, arr in enumerate(sources_and_stations.iterrows()):
        dat = arr[1]
        ax.add_line(Line2D([dat[SOURCE_LONGITUDE], dat[STATION_LONGITUDE]],
                           [dat[SOURCE_LATITUDE], dat[STATION_LATITUDE]],
                           color='b', zorder=i))
    ANZ.drawcoastlines(linewidth=2.0, color='k',
                       zorder=sources_and_stations.shape[0] + 1)

    # ax.set_xlim(reg.leftlon - 5, reg.rightlon + 5)
    # ax.set_ylim(reg.bottomlat - 5, reg.upperlat + 5)
    _draw_paras_merids(ANZ)
    plt.title('Ray paths in \n region {}'.format(region))
    # plt.xlabel('Longitude')
    # plt.ylabel('Latitude')
    fig.savefig('rays_in_region_{}.png'.format(arr_file_base))


def _plot_on_map(sources_and_stations, lon_str, lat_str,
                 marker, color):  # pragma: no cover
    lons = sources_and_stations[lon_str]
    lats = sources_and_stations[lat_str]
    x, y = ANZ(lons, lats)
    ANZ.scatter(x, y, marker=marker, color=color)
    ANZ.drawcoastlines(linewidth=2.0, color='k')
    _draw_paras_merids(ANZ)


def _source_or_stations_in_region(arrivals, region, lat_str, lon_str,
                                  fig_name):  # pragma: no cover
    condition = (
            (arrivals[lat_str] <= region.upperlat)
            &
            (arrivals[lat_str] >= region.bottomlat)
            &
            (arrivals[lon_str] <= region.rightlon)
            &
            (arrivals[lon_str] >= region.leftlon)
    )

    sources_in_region = arrivals[condition]

    _plot_figure(fig_name, lat_str, lon_str, sources_in_region)

    return condition


def _plot_figure(fig_name, lat_str, lon_str,
                 sources_in_region):  # pragma: no cover
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_aspect('equal')
    _plot_on_map(sources_in_region, lon_str, lat_str, marker='*', color='b')
    plt.title(fig_name.split('.')[0])
    # plt.xlabel('Longitude (degrees)')
    # plt.ylabel('Latitude (degrees)')
    fig.savefig(fig_name)


def _draw_paras_merids(m):
    """
    :param m: Basemap instance
    """
    # draw parallels and meridians.
    # label parallels on right and top
    # meridians on bottom and left
    parallels = np.arange(-60., 0, 10.)
    # labels = [left,right,top,bottom]
    m.drawparallels(parallels, labels=[False, True, True, False])
    meridians = np.arange(10., 351., 20.)
    m.drawmeridians(meridians, labels=[True, False, False, True])


@cli.command()
@click.argument('arrivals_file', type=click.File(mode='r'))
@click.argument('region', type=str,
                metavar="str 'upperlat, bottomlat, leftlon, rightlon'")
def plotcmd(arrivals_file, region):
    """
    Another plot command
    """
    from seismic import mpiops
    log.info('Begin plotcmd {}'.format(mpiops.rank))
    log.debug("The input arrival file is %s", arrivals_file)

    return


# ================= Quick Testings of the functions ====================
# cd  passive-seismic/
# export ELLIPCORR=/g/data1a/ha3/fxz547/Githubz/passive-seismic/ellip-corr/

# How to run this script standalone?
# $ fxz547@vdi-n2 /g/data/ha3/fxz547/travel_time_tomography/CSV_NewFormatAug10/FZ01-pst-cluster2/testrun3
# $ python /g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/plotviews.py plot sorted_S.csv '0 -50.0 100 190'
# ======================================================================
if __name__ == "__main__":
    cli()
