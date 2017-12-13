"""
Clustering of events and station for 3d inversion input files.
"""
from __future__ import print_function, absolute_import
import os
from os.path import dirname, join
import random
import logging
import csv
from collections import namedtuple
import fnmatch
from math import asin, sin, acos, sqrt
import numpy as np
import pandas as pd
from matplotlib import pylab as plt
import click
from obspy import read_events
from obspy.geodetics import locations2degrees
from obspy.geodetics.base import WGS84_A as RADIUS
from seismic import pslog
from seismic import mpiops
from inventory.parse_inventory import gather_isc_stations, Station

DPI = asin(1.0)/90.0
R2D = 90./asin(1.)
FLOAT_FORMAT = '%.4f'

log = logging.getLogger(__name__)

SOURCE_LATITUDE = 'source_latitude'
SOURCE_LONGITUDE = 'source_longitude'
STATION_LATITUDE = 'station_latitude'
STATION_LONGITUDE = 'station_longitude'

column_names = ['source_block', 'station_block',
                'residual', 'event_number',
                SOURCE_LONGITUDE, SOURCE_LATITUDE,
                'source_depth', STATION_LONGITUDE, STATION_LATITUDE,
                'observed_tt', 'locations2degrees', 'P_or_S']


PASSIVE = dirname(dirname(dirname(__file__)))
station_metadata = join(PASSIVE, 'inventory', 'stations.csv')
Region = namedtuple('Region', 'upperlat, bottomlat, leftlon, rightlon')


class Grid:
    def __init__(self, nx, ny, dz):
        self.nx = nx
        self.ny = ny
        self.dx = 360.0/nx
        self.dy = 180.0/ny
        self.dz = dz


@click.group()
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cli(verbosity):
    pslog.configure(verbosity)


def recursive_glob(dirname, ext='*.xml'):
    """
    source: https://stackoverflow.com/a/2186565/3321542
    """
    matches = []
    for root, dirnames, filenames in os.walk(dirname):
        for filename in fnmatch.filter(filenames, ext):
            matches.append(os.path.join(root, filename))
    return matches


@cli.command()
@click.argument('events_dir',
                type=click.Path(exists=True, file_okay=True, dir_okay=True,
                                writable=False, readable=True,
                                resolve_path=True))
@click.option('-o', '--output_file',
              type=str, default='outfile',
              help='output arrivals file basename')
@click.option('-x', '--nx', type=int, default=1440,
              help='number of segments from 0 to 360 degrees for longitude')
@click.option('-y', '--ny', type=int, default=720,
              help='number of segments from 0 to 180 degrees for latitude')
@click.option('-z', '--dz', type=float, default=10000.0,
              help='unit segment length of depth in meters')
@click.option('-w', '--wave_type',
              type=click.Choice(['P S', 'Pn Sn', 'Pg Sg', 'p s']),
              default='P S',
              help='Wave type pair to generate inversion inputs')
def gather(events_dir, output_file, nx, ny, dz, wave_type):
    """
    Gather all source-station block pairs for all events in a directory.
    """
    log.info("Gathering all arrivals")

    if os.path.isfile(events_dir):
        event_xmls = [events_dir]
    else:
        event_xmls = recursive_glob(events_dir, ext='*.xml')

    grid = Grid(nx=nx, ny=ny, dz=dz)

    # generate the stations dict
    stations = mpiops.run_once(_read_all_stations)

    process_many_events(event_xmls, grid, stations, wave_type, output_file)

    log.info('Gathered all arrivals in process {}'.format(mpiops.rank))

    mpiops.comm.barrier()

    if mpiops.rank == 0:
        log.info('Now joining all arrivals')
        for t in wave_type.split() + ['missing_stations']:
            _gather_all(output_file, t)


def _gather_all(output_file, s_type):

    final_s_file = output_file + '_' + s_type + '.csv'
    s_arrs = []
    for r in range(mpiops.size):
        s_file = output_file + '_' + s_type + '_{}.csv'.format(r)
        if os.stat(s_file).st_size:
            s_arrs.append(pd.read_csv(s_file, header=None))
        os.remove(s_file)

    if len(s_arrs):
        final_s_df = pd.concat(s_arrs)
        final_s_df.to_csv(final_s_file, header=False, index=False)
    else:
        with open(final_s_file, 'w') as sf:  # just create empty file
            pass


def _read_all_stations():
    stations = read_stations(station_metadata)
    isc_stations = gather_isc_stations()
    stations.update(isc_stations)
    return stations


class ArrivalWriter:
    """
    Convenience class for writing arrival data
    """

    def __init__(self, rank, wave_type, output_file):
        p_type, s_type = wave_type.split()
        p_file = output_file + '_' + p_type + '_{}.csv'.format(rank)
        s_file = output_file + '_' + s_type + '_{}.csv'.format(rank)
        st_file = output_file + '_missing_stations_{}.csv'.format(rank)

        self.p_handle = open(p_file, 'w')
        self.s_handle = open(s_file, 'w')
        self.st_handle = open(st_file, 'w')
        self.p_writer = csv.writer(self.p_handle)
        self.s_writer = csv.writer(self.s_handle)
        self.st_writer = csv.writer(self.st_handle)

    def write(self, cluster_info):
        log.info("Writing cluster info to output file in process {}".format(
            mpiops.rank))

        p_arr, s_arr, missing_stations = cluster_info
        for p in p_arr:
            self.p_writer.writerow(p)
        for s in s_arr:
            self.s_writer.writerow(s)
        for st in missing_stations:
            self.st_writer.writerow([st])

    def close(self):
        if mpiops.rank == 0:
            self.p_handle.close()
            self.s_handle.close()
            self.st_handle.close()


def process_many_events(event_xmls, grid, stations, wave_type, output_file,
                        seed=1):
    total_events = len(event_xmls)

    # when event xmls are of unequal complexity, this shuffle helps
    # distribute the workload evenly amongst processes
    random.seed(seed)
    random.shuffle(event_xmls)
    p_event_xmls = mpiops.array_split(event_xmls, mpiops.rank)

    log.info('Processing {} events of total {} using process {}'.format(
        len(p_event_xmls), total_events, mpiops.rank))

    arrival_writer = ArrivalWriter(mpiops.rank, wave_type=wave_type,
                                   output_file=output_file)

    for i, xml in enumerate(p_event_xmls):
        if xml is not None:
            p_arr = []
            s_arr = []
            missing_stations = []
            log.info('Reading event file {xml}: {i} of {files} in process'
                     ' {process}'.format(i=i+1, files=len(p_event_xmls),
                                         xml=os.path.basename(xml),
                                         process=mpiops.rank))
            # one event xml could contain multiple events
            for e in read_events(xml).events:
                p_arr_t, s_arr_t, m_st = process_event(e, stations, grid,
                                                       wave_type)
                p_arr += p_arr_t
                s_arr += s_arr_t
                missing_stations += m_st
                log.debug('processed event {e} from {xml}'.format(
                    e=e.resource_id, xml=xml))

            arrival_writer.write([p_arr, s_arr, missing_stations])

    log.info('Read all events in process {}'.format(mpiops.rank))
    arrival_writer.close()


def process_event(event, stations, grid, wave_type):
    """
    :param event: obspy.core.event.Event class instance
    :param stations: dict
        stations dict
    :param grid: Grid class instance
    :param wave_type: str
        Wave type pair to generate inversion inputs. See `gather` function.
    """
    p_type, s_type = wave_type.split()

    # use preferred origin timestamp as the event number
    # if preferred origin is not populated, use the first origin timestamp
    origin = event.preferred_origin() or event.origins[0]
    ev_number = int(origin.time.timestamp)

    p_arrivals = []
    s_arrivals = []
    missing_stations = []

    # other event parameters we need
    ev_latitude = origin.latitude
    ev_longitude = origin.longitude
    ev_depth = origin.depth

    if ev_latitude is None or ev_longitude is None or ev_depth is None:
        return p_arrivals, s_arrivals, missing_stations

    event_block = _find_block(grid, ev_latitude, ev_longitude, z=ev_depth)

    for arr in origin.arrivals:
        sta_code = arr.pick_id.get_referred_object(
        ).waveform_id.station_code

        # ignore arrivals not in stations dict, workaround for now for
        # ENGDAHL/ISC events
        # TODO: remove this condition once all ISC/ENGDAHL stations are
        # available
        # Actually it does not hurt retaining this if condition. In case,
        # a station comes in which is not in the dict, the data prep will
        # still work
        # Note some stations are still missing even after taking into account
        #  of all seiscomp3 stations, ISC and ENGDAHL stations
        if sta_code not in stations:
            log.warning('Station {} not found in inventory'.format(sta_code))
            missing_stations.append(str(sta_code))
            continue
        sta = stations[sta_code]

        degrees_to_source = locations2degrees(ev_latitude, ev_longitude,
                                              float(sta.latitude),
                                              float(sta.longitude))

        # ignore stations more than 90 degrees from source
        if degrees_to_source > 90.0:
            # log.info('Ignored this station arrival as distance from source '
            #          'is {} degrees'.format(degrees_to_source))
            continue

        # TODO: use station.elevation information
        station_block = _find_block(grid,
                                    float(sta.latitude), float(sta.longitude),
                                    z=0.0)

        # phase_type == 1 if P and 2 if S
        if arr.phase in wave_type.split():

            t_list = [event_block, station_block, arr.time_residual,
                      ev_number, ev_longitude, ev_latitude, ev_depth,
                      sta.longitude, sta.latitude,
                      (arr.pick_id.get_referred_object().time.timestamp -
                       origin.time.timestamp), degrees_to_source]

            p_arrivals.append(t_list + [1]) if arr.phase == p_type else \
                s_arrivals.append(t_list + [2])
        else:  # ignore the other phases
            pass
    return p_arrivals, s_arrivals, missing_stations


def _find_block(grid, lat, lon, z):
    y = 90. - lat
    x = lon if lon > 0 else lon + 360.0
    i = round(x / grid.dx) + 1
    j = round(y / grid.dy) + 1
    k = round(z / grid.dz) + 1
    block_number = (k - 1) * grid.nx * grid.ny + (j - 1) * grid.nx + i
    return int(block_number)


def read_stations(station_file):
    """
    Read station location from a csv file.
    :param station_file: str
        csv stations file handle passed in by click
    :return: stations_dict: dict
        dict of stations indexed by station_code for quick lookup
    """
    log.info('Reading seiscomp3 exported stations file')
    stations_dict = {}
    with open(station_file, 'r') as csv_file:
        reader = csv.reader(csv_file)
        next(reader)  # skip header
        for station in map(Station._make, reader):
            stations_dict[station.station_code] = station
        log.info('Done reading seiscomp3 station files')
        return stations_dict


@cli.command()
@click.argument('output_file',
                type=click.File(mode='r'))
@click.argument('residual_cutoff', type=float)
@click.option('-s', '--sorted_file',
              type=click.File(mode='w'), default='sorted.csv',
              help='output sorted and filter file.')
def sort(output_file, sorted_file, residual_cutoff):
    """
    Sort and filter the arrivals.

    Sort based on the source and station block number.
    There are two stages of filtering:
    1. Filter based on the time residual
    2. Filter based on median of observed travel time.

    If there are multiple source and station block combinations, we keep the
    row corresponding to the median observed travel time (observed_tt).

    :param output_file: output file from the gather stage
    :param sorted_file: str, optional
        optional sorted output file path. Default: sorted.csv.
    :param residual_cutoff: float
        residual seconds above which arrivals are rejected.
    :return: None
    """

    log.info('Filtering arrivals.')

    cluster_data = pd.read_csv(output_file, header=None,
                               names=column_names)
    cluster_data = cluster_data[abs(cluster_data['residual'])
                                < residual_cutoff]

    # groupby sorts by default
    # cluster_data.sort_values(by=['source_block', 'station_block'],
    #                          inplace=True)

    log.info('Sorting arrivals.')

    # groupby automatically sorts
    med = cluster_data.groupby(by=['source_block',
                                   'station_block'])[
        'observed_tt'].quantile(q=.5, interpolation='lower').reset_index()

    final_df = pd.merge(cluster_data, med, how='right',
                        on=['source_block', 'station_block', 'observed_tt'],
                        sort=True,
                        right_index=True)

    # Confirmed: drop_duplicates required due to possibly duplicated picks in
    #  the original engdahl events
    # refer: https://github.com/GeoscienceAustralia/passive-seismic/issues/51
    # The subset is specified as we have some stations that are very close?
    final_df.drop_duplicates(subset=['source_block', 'station_block',
                                     'event_number', 'source_longitude',
                                     'source_latitude', 'source_depth'],
                             keep='first',
                             inplace=True)

    final_df.to_csv(sorted_file, header=True, index=False)


@cli.command()
@click.argument('p_file', type=click.File(mode='r'))
@click.argument('s_file', type=click.File(mode='r'))
@click.option('-p', '--matched_p_file',
              type=click.File(mode='w'), default='matched_p.csv',
              help='output matched p file.')
@click.option('-s', '--matched_s_file',
              type=click.File(mode='w'), default='matched_s.csv',
              help='output matched s file.')
def match(p_file, s_file, matched_p_file, matched_s_file):
    """
    Match source and station blocks and output files with matched source and
    station blocks.

    :param p_file: str
        path to sorted P arrivals
    :param s_file: str
        path to sorted S arrivals
    :param matched_p_file: str, optional
        output p arrivals file with matched p and s arrivals source-block
        combinations
    :param matched_s_file: str, optional
        output s arrivals file with matched p and s arrivals source-block
        combinations

    :return:None
    """

    log.info('Matching p and s arrivals')

    p_arr = pd.read_csv(p_file)
    s_arr = pd.read_csv(s_file)

    blocks = pd.merge(p_arr[['source_block', 'station_block']],
                      s_arr[['source_block', 'station_block']],
                      how='inner',
                      on=['source_block', 'station_block'])
    matched_P = pd.merge(p_arr, blocks, how='inner',
                         on=['source_block', 'station_block'])[column_names]
    matched_S = pd.merge(s_arr, blocks, how='inner',
                         on=['source_block', 'station_block'])[column_names]
    matched_P.to_csv(matched_p_file, index=False, header=False)
    matched_S.to_csv(matched_s_file, index=False, header=False)


@cli.command()
@click.argument('region', type=str,
                metavar="str 'upperlat, bottomlat, leftlon, rightlon'")
@click.argument('matched_file', click.File(mode='r'),
                metavar='cluster_matched_file')
@click.option('-r', '--region_file', type=click.File('w'),
              default='region.csv',
              help='region file name.')
@click.option('-g', '--global_file', type=click.File('w'),
              default='global.csv',
              help='global file name.')
@click.option('-s', '--grid_size', type=float, default=0.0,
              help='grid size in degrees in the region. If grid size is '
                   'provided, cross region file will be created.')
@click.option('-c', '--cross_region_file', type=click.File('w'),
              default='cross_region.csv',
              help='cross region file name.')
def zone(region, matched_file, region_file, global_file, cross_region_file,
         grid_size):
    """
    `zone'ing the arrivals into three regions.
    """
    log.info('Calculating zones')
    region = [float(s) for s in region.split()]
    region = Region(*region)

    matched = pd.read_csv(matched_file, header=None, names=column_names)

    # convert longitude to co-longitude
    matched['source_longitude'] = matched['source_longitude'] % 360
    matched['station_longitude'] = matched['station_longitude'] % 360

    # if either the source or the station or both are inside region
    # else global, unless we want a cross region

    _in_region(region, matched, region_file=region_file,
               global_file=global_file, grid_size=grid_size,
               cross_region_file=cross_region_file)


@cli.command()
@click.argument('arrivals_file', type=click.File(mode='r'))
@click.argument('region', type=str,
                metavar="str 'upperlat, bottomlat, leftlon, rightlon'")
def plot(arrivals_file, region):
    """
    This command will output a `sources_in_region.png`, which will show all the
    sources inside the `region` specified by the region string which can be
    specified like '0 -50.0 100 160'. It will also output
    a`stations_in_region.png` showing the stations where arrivals were
    recorded within the `region`.
    The `cluster plot` command further outputs a
    `sources_and_stations_in_region.png` which should all sources and
    stations in the same plot that is within `region`.

    Note that `cluster plot` only accepts the final `zone` output files.
    """
    region = [float(s) for s in region.split()]
    region_tuple = Region(*region)
    arrivals = pd.read_csv(arrivals_file, header=None, names=column_names,
                           sep=' ')
    source = _source_or_stations_in_region(
        arrivals, region_tuple, SOURCE_LATITUDE, SOURCE_LONGITUDE,
        'sources_in_region.png')

    station = _source_or_stations_in_region(
        arrivals, region_tuple, STATION_LATITUDE, STATION_LONGITUDE,
        'stations_in_region.png')

    # sources and stations both in region
    sources_and_stations = arrivals[source & station]

    fig = plt.figure()
    plt.plot(sources_and_stations[SOURCE_LONGITUDE],
             sources_and_stations[SOURCE_LATITUDE], 'r*')
    plt.plot(sources_and_stations[STATION_LONGITUDE],
             sources_and_stations[STATION_LATITUDE], 'b^')
    plt.title('Sources and stations in \n region {}'.format(region))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig.savefig('sources_and_stations_in_region.png')
    fig.clear()

    # rays
    ax = fig.add_subplot(111)
    from matplotlib.lines import Line2D

    for arr in sources_and_stations.iterrows():
        dat = arr[1]
        ax.add_line(Line2D([dat[SOURCE_LONGITUDE], dat[STATION_LONGITUDE]],
                           [dat[SOURCE_LATITUDE], dat[STATION_LATITUDE]]))

    ax.set_xlim(region_tuple.leftlon - 5, region_tuple.rightlon + 5)
    ax.set_ylim(region_tuple.bottomlat - 5, region_tuple.upperlat + 5)
    plt.title('Ray paths in \n region {}'.format(region))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    fig.savefig('rays_in_region.png')


def _source_or_stations_in_region(arrivals, region, lat_str, lon_str,
                                  fig_name):
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


def _plot_figure(fig_name, lat_str, lon_str, sources_in_region):
    fig = plt.figure()
    plt.plot(sources_in_region[lon_str],
             sources_in_region[lat_str], '*')
    plt.title(fig_name.split('.')[0])
    plt.xlabel('Longitude (degrees)')
    plt.ylabel('Latitude (degrees)')
    fig.savefig(fig_name)


def _in_region(region, df, region_file, global_file, grid_size,
               cross_region_file):

    if grid_size > 0.0:
        df = _intersect_region(df, region, grid_size)

    # row indices of all in region arrivals
    df_region = df[
            (
                (
                    (region.leftlon < df['source_longitude']) &
                    (df['source_longitude'] < region.rightlon)
                )
                &
                (
                    (region.bottomlat < df['source_latitude']) &
                    (df['source_latitude'] < region.upperlat)
                )
            )
            |
            (
                (
                    (region.leftlon < df['station_longitude']) &
                    (df['station_longitude'] < region.rightlon)
                )
                &
                (
                    (region.bottomlat < df['station_latitude']) &
                    (df['station_latitude'] < region.upperlat)
                )
            )
    ][column_names]

    # dataframe excluding in region arrivals
    df_ex_region = df.iloc[df.index.difference(df_region.index)]

    if grid_size > 0.0:
        # cross region is in ex-region and cross-region==True
        df_ex_region[
            df_ex_region['cross_region'] == True][column_names].to_csv(
            cross_region_file, index=False, header=False,
            sep=' ', float_format=FLOAT_FORMAT)
        global_df = df_ex_region[df_ex_region['cross_region'] == False][
                    column_names]
    else:
        global_df = df_ex_region

    # Global region contain the remaining arrivals
    global_df.to_csv(global_file, index=False, header=False,
                     sep=' ', float_format=FLOAT_FORMAT)

    df_region.to_csv(region_file, index=False, header=False,
                     sep=' ', float_format=FLOAT_FORMAT)


def _intersect_region(df, region, grid_size):
    """
    Strategy to compute cross region: Intersect/cross region is computed first
    which will contain the `region`. The final intersect region will be
    be subtracted from the `region`.
    """

    pe = df['source_latitude']
    ps = df['station_latitude']
    re = df['source_longitude']
    rs = df['station_longitude']
    delta = df['locations2degrees']

    # operations on pd.Series
    nms = (delta/grid_size).astype(int)
    ar = pe*DPI
    ast = ps*DPI
    br = re*DPI
    bs = rs*DPI

    x1 = RADIUS*np.sin(ar)*np.cos(br)
    y1 = RADIUS*np.sin(ar)*np.sin(br)
    z1 = RADIUS*np.cos(ar)
    x2 = RADIUS*np.sin(ast)*np.cos(bs)
    y2 = RADIUS*np.sin(ast)*np.sin(bs)
    z2 = RADIUS*np.cos(ast)
    dx = (x2-x1)/nms
    dy = (y2-y1)/nms
    dz = (z2-z1)/nms

    in_cross = []

    # TODO: vectorize this loop
    for i, n in enumerate(nms):
        in_cross.append(_in_cross_region(dx[i], dy[i], dz[i], n, region, x1[i],
                        y1[i], z1[i]))
    df['cross_region'] = pd.Series(in_cross)
    return df


def _in_cross_region(dx, dy, dz, nms, region, x1, y1, z1):

    # TODO: vectorize this loop
    # TODO: tests for cross region
    for j in range(nms):

        x = x1 + dx * j
        y = y1 + dy * j
        z = z1 + dz * j
        r = sqrt(x ** 2 + y ** 2 + z ** 2)
        acosa = z / r
        if acosa < -1.:
            acosa = -1.

        if acosa > 1:
            acosa = 1.

        lat = acos(acosa) * R2D

        acosa = (x / r) / sin(lat * DPI)

        if acosa < -1.:
            acosa = -1.

        if acosa > 1.:
            acosa = 1.

        lon = acos(acosa) * R2D

        if y < 0.0:
            lon = 360.0 - lon

        if (lon > region.leftlon) and (lon < region.rightlon):
            if (lat > region.bottomlat) and (lat < region.upperlat):
                if (RADIUS - r) < 1000.0:
                    return True
    return False
