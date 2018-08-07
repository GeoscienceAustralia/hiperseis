"""
Clustering of events and station for 3d inversion input files.
"""
from __future__ import print_function, absolute_import

import csv
import fnmatch
import logging
import os
import random
from collections import namedtuple
from math import asin, sin, acos, sqrt
import matplotlib
import numpy as np
import pandas as pd

# Force matplotlib to not use any Xwindows backend
matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.lines import Line2D
from mpl_toolkits.basemap import Basemap
import click
from obspy import read_events
from obspy.geodetics import locations2degrees, gps2dist_azimuth
from obspy.geodetics.base import WGS84_A as RADIUS
from seismic import pslog
from seismic import mpiops
import ellipcorr
from inventory.parse_inventory import read_all_stations

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
ANZ = Basemap(llcrnrlon=100.0, llcrnrlat=-50.0,
              urcrnrlon=190.0, urcrnrlat=0.0)

Region = namedtuple('Region', 'upperlat, bottomlat, leftlon, rightlon')

PARAM_FILE_FORMAT = '''
    An example parameter file is provided in raytracer/params/param2x2. 
    A typical param file should have the following format:
    
    Global dataset following parameters:
    72 36 16 5. 5.
          0.
        110.
        280.
        410.
        660.
        840.
       1020.
       1250.
       1400.
       1600.
       1850.
       2050.
       2250.
       2450.
       2600.
       2750.
       2889.
       
    where 72 is number of cells in Lon, 36 number of cells in Lat, 
    16 is number of layers, and 5 is size of the cell
    
    For local parameterisation:
     
     
    100.  190.  -54.  0.  45 27  22
          0.
         35.
         70.
        110.
        160.
        210.
        260.
        310.
        360.
        410.
        460.
        510.
        560.
        610.
        660.
        710.
        810.
        910.
       1010.
       1110.
       1250.
       1400.
       1600.
     
    where 100, 190 are minLon and maxLon, '-54' and '0' are minlat and maxlat, 
    45 number of cells in Lon, 27 is the number of cells in lat, 22 
    is the number of layers
   '''


class Grid:
    """
    Encaptulate 3D cell properties in Uniform grid
    """

    def __init__(self, nx=360, ny=180, dz=10000):
        """
        constructor for Grid
        :param nx: integer number of cells in longitude (360)
        :param ny: lattitude cells (180)
        :param dz: depth in meteres (default=10KM)
        """
        self.nx = nx
        self.ny = ny
        self.dx = 360.0 / nx
        self.dy = 180.0 / ny
        self.dz = dz

    def find_block_number(self, lat, lon, z):
        """
        find the 3D-block number in this grid
        :param lat: lattitude
        :param lon: longitude
        :param z: elevation
        :return: int block number
        """
        y = 90. - lat
        x = lon % 360
        i = round(x / self.dx) + 1
        j = round(y / self.dy) + 1
        k = round(z / self.dz) + 1
        block_number = (k - 1) * self.nx * self.ny + (j - 1) * self.nx + i

        return int(block_number)

class Grid2:
    """
    A non-uniform Grid Model for source->events rays clustering and sorting.
    according to the inversion Fortran code model (file param1x1)
    72 36 16 5. 5.  (5x5degree global + a list of 1+16 depths)

    100. 190. -54. 0. 90 54  23  (1x1 degree in ANZ region)  + a list of 1+23 depths

    """

    def __init__(self, nx=360, ny=180, dz=10000):
        """
        constructor for a non-uniform Grid
        :param nx: integer number of cells in longitude (360)
        :param ny: lattitude cells (180)
        :param dz: depth in meteres (def=10000
        """

        # Region bounadry definition
        self.LON=(100.0, 190.0)
        self.LAT=(-54.0, 0.0)

        # region grid size in lat-long plane
        self.nx = 4 * 90  # 1/4 degree cell size in LONG
        self.ny = 4 * 54  # 1/4 degree cell size in LAT
        self.dx = (self.LON[1]-self.LON[0]) / self.nx  # 1/4 degree
        self.dy = (self.LAT[1]-self.LAT[0]) / self.ny  # 1/4 degree

        # Region's depths in KM according to Fortran param1x1
        self.rdepth=(0, 10, 35, 70, 110, 160, 210, 260, 310, 360, 410, 460, 510, 560,
                     610, 660, 710, 810, 910, 1010, 1110, 1250, 1400, 1600)

        # Region Depth in meters
        self.rmeters=1000*np.array(self.rdepth)
        self.dz=5000  # inside the ANZ zone 5KM uniform cellsize in depth


        # global grid model params (outside the region)
        self.gnx = 288
        self.gny = 144
        self.gdx = 5 * self.dx  # =360/self.gnx;  5 times as big as the regional grid size
        self.gdy = 5 * self.dy  # =180/self.gny;  5 times as big as the regional grid size

        self.gdepth=(0, 110, 280, 410, 660, 840, 1020, 1250, 1400, 1600, 1850, 2050, 2250, 2450, 2600, 2750, 2889)
        # Region Depth in meters
        self.gmeters = 1000 * np.array(self.gdepth)
        self.gdz=20000  # outside the ANZ zone 20KM uniform cellsize in depth


        self.REGION_MAX_BN = self._get_max_region_block_number()


        return

    def _get_max_region_block_number(self):
        """
        Compute the MAX BLOCK NUMBER for a point in the region box, with a fine grid,
        assuming the event's max depth is 2000KM
        :return: int max_nb
        """

        i = round(90.0 / self.dx) + 1
        j = round(54.0 / self.dy) + 1

        k = round(2000000.0 / self.dz) + 1  # assume 2000KM max depth
        bn = (k - 1) * self.nx * self.ny + (j - 1) * self.nx + i

        log.info("The REG_MAX_BN = %s", bn)

        return bn

    def is_point_in_region(self, lat, lon):
        """
        test if the event  or station point in the region box?
        :param lat:
        :param lon:
        :return: T/F
        """
        x = lon % 360  # convert lon into x which must be in [0,360)
        #y = 90. - lat  # convert lat into y which will be in [0,180)

        if (lat <= self.LAT[1] and lat >= self.LAT[0] and x <= self.LON[1] and x >= self.LON[0]):
            log.debug("Point in the region = (%s, %s)", lat, lon)
            return True
        else:
            log.debug("Point not in the region, = (%s, %s)", lat, lon)
            return False


    def find_block_number(self, lat, lon, z):
        """
        find the index-number of an events/station in a non-uniform grid
        each spatial point (lat, lon, z) is mapped to a uniq block_number.
        :param lat: latitude
        :param lon: longitude
        :param z: depth
        :return: int block number
        """
        x = lon % 360  # convert lon into x which must be in [0,360)
        y = 90. - lat  # convert lat into y which will be in [0,180)

        if self.is_point_in_region(lat,lon) is True:

            i = round(x / self.dx) + 1
            j = round(y / self.dy) + 1

            k = round(z / self.dz) + 1
            block_number = (k - 1) * self.nx * self.ny + (j - 1) * self.nx + i

        else:

            i = round(x / self.gdx) + 1
            j = round(y / self.gdy) + 1

            k = round(z / self.gdz) + 1

            g_block_number = (k - 1) * self.gnx * self.gny + (j - 1) * self.gnx + i

            block_number = g_block_number + self.REGION_MAX_BN  # offset by the region max block number


        return int(block_number)


@click.group()
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cli(verbosity):

    pslog.configure(verbosity)

def recursive_glob(dirname, ext='*.xml'):
    """
    Under the dirname recursively find all files with extension ext.
    Return a list of the full-path to the files of interest.
    See: https://stackoverflow.com/a/2186565/3321542
    :param dirname: a single dir OR a list of dirs.
    :param ext: eg, ".xml"
    :return: a list of path2files
    """

    if isinstance(dirname, (list,)): # the input argument is a list of directory
        filelist = []
        for adir in dirname:
            filelist.extend(recursive_glob(adir))
        return filelist
    else: # input variable is a single dir
        matches = []
        for root, dirnames, filenames in os.walk(dirname):
            for filename in fnmatch.filter(filenames, ext):
                matches.append(os.path.join(root, filename))
        return matches


def get_paths_from_csv(csvfile):
    """
    Parse a text/csv file to extract a list of paths, where events xml files are stored, to be gathered.
    :param csvfile: csv file
    :return: list_of_paths
    """

    paths = []  # ["/g/data/ha3/events_xmls_sc3ml/", "/g/data/ha3/fxz547/travel_time_tomography/new_events20180516/"]

    with open(csvfile) as csvf:
        reader = csv.reader(csvf)

        for row in reader:
            print(', '.join(row))

            paths.extend(row)

    return paths


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
    Gather all source-station block pairs for all events in a list of directory.
    """
    log.info("Gathering all arrivals")

    if os.path.isfile(events_dir):  # is a text csv file containing multiple dirs.
        event_dirs = get_paths_from_csv(events_dir)
        event_xmls = recursive_glob(event_dirs)
    elif os.path.isdir(events_dir):
        event_xmls = recursive_glob(events_dir, ext='*.xml')
    else:
        event_xmls = None
        raise Exception("Invalid Input Events Dir")

    # grid = Grid(nx=nx, ny=ny, dz=dz)  # original uniform grid
    grid = Grid2( )  # use a non-uniform Grid construction

    # generate the stations dict
    stations = mpiops.run_once(read_all_stations)

    process_many_events(event_xmls, grid, stations, wave_type, output_file)

    log.info('Gathered all arrivals in process {}'.format(mpiops.rank))

    mpiops.comm.barrier()

    if mpiops.rank == 0:
        log.info('Now joining all arrivals')
        for t in wave_type.split() + ['missing_stations',
                                      'participating_stations']:
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


class ArrivalWriter:
    """
    Convenience class for writing arrival data
    """

    def __init__(self, rank, wave_type, output_file):
        p_type, s_type = wave_type.split()
        p_file = output_file + '_' + p_type + '_{}.csv'.format(rank)
        s_file = output_file + '_' + s_type + '_{}.csv'.format(rank)
        miss_st_file = output_file + '_missing_stations_{}.csv'.format(rank)
        st_file = output_file + '_participating_stations_{}.csv'.format(rank)

        self.p_handle = open(p_file, 'w')
        self.s_handle = open(s_file, 'w')
        self.miss_st_handle = open(miss_st_file, 'w')
        self.st_handle = open(st_file, 'w')

        self.p_writer = csv.writer(self.p_handle)
        self.s_writer = csv.writer(self.s_handle)
        self.missing_st_writer = csv.writer(self.miss_st_handle)
        self.st_writer = csv.writer(self.st_handle)

    def write(self, cluster_info):
        log.info("Writing cluster info to output file in process {}".format(
            mpiops.rank))

        p_arr, s_arr, missing_stations, arr_stations = cluster_info
        for p in p_arr:
            self.p_writer.writerow(p)
        for s in s_arr:
            self.s_writer.writerow(s)

        for st in missing_stations:
            self.missing_st_writer.writerow([st])

        for st in arr_stations:
            self.st_writer.writerow([st])

    def close(self):
        if mpiops.rank == 0:
            self.p_handle.close()
            self.s_handle.close()
            self.miss_st_handle.close()
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
    process_event_counter = 0

    for i, xml in enumerate(p_event_xmls):
        if xml is not None:
            p_arr = []
            s_arr = []
            missing_stations = []
            arriving_stations = []
            log.info('Reading event file {xml}: {i} of {files} in process'
                     ' {process}'.format(i=i + 1, files=len(p_event_xmls),
                                         xml=os.path.basename(xml),
                                         process=mpiops.rank))
            # one event xml could contain multiple events
            try:
                for e in read_events(xml).events:
                    process_event_counter += 1
                    p_arr_t, s_arr_t, m_st, a_st = process_event(
                        e, stations, grid, wave_type, process_event_counter)
                    p_arr += p_arr_t
                    s_arr += s_arr_t
                    missing_stations += m_st
                    arriving_stations += a_st

                    log.debug('processed event {e} from {xml}'.format(
                        e=e.resource_id, xml=xml))
            except ValueError as e:
                log.warning('ValueError in processing event {}'.format(xml))
                log.warning(e)
            except Exception as e:
                log.warning('Unknown Exception in '
                            'processing event {}'.format(xml))
                log.warning(e)

            arrival_writer.write([p_arr, s_arr, missing_stations,
                                  arriving_stations])

    log.info('Read all events in process {}'.format(mpiops.rank))
    arrival_writer.close()


def process_event(event, stations, grid, wave_type, counter):
    """
    :param event: obspy.core.event.Event class instance
    :param stations: dict
        stations dict
    :param grid: Grid class instance
    :param wave_type: str
        Wave type pair to generate inversion inputs. See `gather` function.
    :param counter: int
        event counter in this process
    """
    p_type, s_type = wave_type.split()

    # use preferred origin timestamp as the event number
    # if preferred origin is not populated, use the first origin timestamp
    origin = event.preferred_origin() or event.origins[0]

    # TODO: longer term uncomment the following line when fortran TT code is
    # updated to use longer integers
    # ev_number = int(origin.time.timestamp)

    # TODO: delete this definition of ev_number when fortran code can use
    # longer integers

    # the following definition ensures a 8 digit event number that is also
    # unique
    assert counter < 100000, 'Counter must be less than 100000'
    ev_number = int(str(counter) + '{0:0=3d}'.format(mpiops.rank))

    p_arrivals = []
    s_arrivals = []
    missing_stations = []
    arrival_staions = []

    # other event parameters we need
    ev_latitude = origin.latitude
    ev_longitude = origin.longitude
    ev_depth = origin.depth

    if ev_latitude is None or ev_longitude is None or ev_depth is None:
        return p_arrivals, s_arrivals, missing_stations, arrival_staions

    event_block = grid.find_block_number(ev_latitude, ev_longitude, z=ev_depth)

    for arr in origin.arrivals:

        snr_value = getSNR(arr)
        log.debug("Arrival Pick SNR value: %s", snr_value)

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

        log.debug("events and station latlong: %s, %s, %s, %s", ev_latitude, ev_longitude,
                                              sta.latitude, sta.longitude)

        try:
            degrees_to_source = locations2degrees(ev_latitude, ev_longitude,
                                              sta.latitude,
                                              sta.longitude)
        except Exception as e:
            log.warning("location to degree error %s", e )

        log.debug("location to degree= %s", degrees_to_source)
        # ignore stations more than 90 degrees from source
        if degrees_to_source > 90.0:
            # log.info('Ignored this station arrival as distance from source '
            #          'is {} degrees'.format(degrees_to_source))
            continue

        # TODO: use station.elevation information
        station_block = grid.find_block_number(sta.latitude, sta.longitude, z=0.0)

        if arr.phase in wave_type.split():
            log.debug("Began ellipticity_corr ")

            azim_v = gps2dist_azimuth(ev_latitude, ev_longitude, sta.latitude, sta.longitude)[1]

            log.debug("Check input params to ellipticity_corr = %s, %s, %s, %s, %s", arr.phase, degrees_to_source, ev_depth, 90-ev_latitude, azim_v )

            ellipticity_corr = ellipcorr.ellipticity_corr(
                phase=arr.phase,
                edist=degrees_to_source,
                edepth=ev_depth / 1000.0,
                # TODO: check co-latitude definition
                # no `ecolat` bounds check in fortran ellipcorr subroutine
                # no `origin.latitude` bounds check in obspy
                ecolat=90 - ev_latitude,  # conversion to co-latitude
                azim= azim_v
            )

            log.debug("ellipticity_corr = %s", ellipticity_corr)
            
            t_list = [event_block, station_block, arr.time_residual,
                      ev_number, ev_longitude, ev_latitude, ev_depth,
                      sta.longitude, sta.latitude,
                      (arr.pick_id.get_referred_object().time.timestamp -
                       origin.time.timestamp) + ellipticity_corr,
                      degrees_to_source,
                      sta_code, snr_value]
            arrival_staions.append(sta_code)
            p_arrivals.append(t_list + [1]) if arr.phase == p_type else \
                s_arrivals.append(t_list + [2])
        else:  # ignore the other phases
            pass
    return p_arrivals, s_arrivals, missing_stations, arrival_staions

# FZ: moved this function to the grid class
# def _find_block(grid, lat, lon, z):
#     y = 90. - lat
#     x = lon % 360
#     i = round(x / grid.dx) + 1
#     j = round(y / grid.dy) + 1
#     k = round(z / grid.dz) + 1
#     block_number = (k - 1) * grid.nx * grid.ny + (j - 1) * grid.nx + i
#     return int(block_number)

def getSNR(arrival):
    """
    From the arrival get the SNR value.
    This algorithm depend on how the snr value is coded in the xml file
    :param arrival:
    :return: a float SNR value
    """
    snr_v = arrival.pick_id.get_referred_object().comments[3]  # Comment(text='snr = 10.7157568852')

    snrlist = str(snr_v).split("snr =")
    snrv = snrlist[-1][:-2]  # the last item of the split, trimming two chars ')

    return float(snrv)

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

    cmdline usage:
    cluster sort outfile_P.csv 5. -s sorted_P.csv
    cluster sort outfile_S.csv 10. -s sorted_S.csv


    :param output_file: output file from the gather stage (eg, outfile_P.csv)
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
    cluster_data['source_depth'] = cluster_data['source_depth'] / 1000.0  # scale meter to KM
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
                        sort=True, right_index=True)

    # Confirmed: drop_duplicates required due to possibly duplicated picks in
    #  the original engdahl events
    # refer: https://github.com/GeoscienceAustralia/passive-seismic/issues/51
    # The subset is specified as we have some stations that are very close?
    final_df.drop_duplicates(subset=['source_block', 'station_block',
                                     'event_number', SOURCE_LONGITUDE,
                                     SOURCE_LATITUDE, 'source_depth'],
                             keep='first', inplace=True)
    # FZ: this drop_duplicate will still keep many identical duplicated rows with only event_number different

    final_df.to_csv(sorted_file, header=False, index=False, sep=' ')



@cli.command()
@click.argument('output_file',
                type=click.File(mode='r'))
@click.argument('residual_cutoff', type=float)
@click.option('-s', '--sorted_file',
              type=click.File(mode='w'), default='sorted2.csv',
              help='output sorted and filter file.')
def sort2(output_file, sorted_file, residual_cutoff):
    """
    Sort and filter the arrivals.

    Sort based on the source and station block number.
    There are two stages of filtering:
    1. Filter based on the time residual
    2. Filter based on best Signal_to_Noise-Ratio seismic wave: If there are multiple source and station block combinations, we keep the
    row corresponding to the highest SNR value

    cmdline usage:
    cluster sort outfile_P.csv 5. -s sorted_P.csv
    cluster sort outfile_S.csv 10. -s sorted_S.csv


    :param output_file: output file from the gather stage (eg, outfile_P.csv)
    :param sorted_file: str, optional
        optional sorted output file path. Default: sorted.csv.
    :param residual_cutoff: float
        residual seconds above which arrivals are rejected.
    :return: pandas_df
    """

    log.info('Filtering arrivals.')

    cluster_data = pd.read_csv(output_file, header=None,
                               names=column_names)

    cluster_data = cluster_data[abs(cluster_data['residual'])
                                < residual_cutoff]

    cluster_data['source_depth'] = cluster_data['source_depth'] / 1000.0  # convert to KM?

    # groupby sorts by default
    # cluster_data.sort_values(by=['source_block', 'station_block'],
    #                          inplace=True)

    log.info('Sorting arrivals.')

    # groupby automatically sorts
    #     med = cluster_data.groupby(
    #         by=['source_block', 'station_block']
    #         )['observed_tt'].quantile(q=.5, interpolation='lower').reset_index() # use a seq index:0,1,2,....

    #  SNR value max
    med = cluster_data.groupby(by=['source_block', 'station_block'])['SNR'].max().reset_index()  # use a seq index:0,1,2,.

    # med dataframe has three columns:  [source_block, station_block ,observed_tt]

    final_df = pd.merge(cluster_data, med, how='right',
                        on=['source_block', 'station_block', 'SNR'],
                        sort=True,
                        right_index=True)

    # Confirmed: drop_duplicates required due to possibly duplicated picks in
    #  the original engdahl events
    # refer: https://github.com/GeoscienceAustralia/passive-seismic/issues/51
    # The subset is specified as we have some stations that are very close?
    final_df.drop_duplicates(subset=['source_block', 'station_block'],
                             keep='first', inplace=True)

    final_df.to_csv(sorted_file, header=False, index=False, sep=' ')

    return final_df


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

    p_arr = pd.read_csv(p_file, header=None, names=column_names, sep=' ')
    s_arr = pd.read_csv(s_file, header=None, names=column_names, sep=' ')

    blocks = pd.merge(p_arr[['source_block', 'station_block']],
                      s_arr[['source_block', 'station_block']],
                      how='inner',
                      on=['source_block', 'station_block'])
    matched_P = pd.merge(p_arr, blocks, how='inner',
                         on=['source_block', 'station_block'])[column_names]
    matched_S = pd.merge(s_arr, blocks, how='inner',
                         on=['source_block', 'station_block'])[column_names]
    matched_P.to_csv(matched_p_file, index=False, header=False, sep=' ')
    matched_S.to_csv(matched_s_file, index=False, header=False, sep=' ')


@cli.command()
@click.option('-z', '--region', type=str, default='',
              metavar="str 'upperlat, bottomlat, leftlon, rightlon'")
@click.option('-p', '--parameter_file', type=str, default='',
              metavar="inversion_parameter_file")
@click.argument('matched_file', click.File(mode='r'),
                metavar='cluster_matched_or_sorted_file')
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
@click.option('-t', '--stats', type=bool, default=True,
              help='Calculate station stats switch.')
@click.option('-j', '--reject_stations_file', type=click.File('r'),
              default=None, help='Calculate station stats switch.')
def zone(region, parameter_file, matched_file, region_file, global_file,
         cross_region_file, grid_size, stats, reject_stations_file):
    """
    `zone'ing the arrivals into three regions.
    Note: Arrivals don't have to be `match`ed for `zone`ing. Sorted P/p and S/s
    arrivals can also be used for `zone`ing.
    """

    log.info('Calculating zones')

    region = Region(*_get_region_string(parameter_file, region))

    matched = pd.read_csv(matched_file, header=None, names=column_names,
                          sep=' ')
    df_region, global_df, x_region_df = _in_region(region, matched, grid_size)

    if reject_stations_file is not None:
        reject_stations = pd.read_csv(reject_stations_file, header=None,
                                      names=[STATION_CODE])
        reject_stations_set = set(reject_stations[STATION_CODE].values)
        r_rows = [False if (x in reject_stations_set) else True for x
                  in df_region[STATION_CODE]]
        df_region = df_region[r_rows]
        g_rows = [False if (x in reject_stations_set) else True for x
                  in global_df[STATION_CODE]]
        global_df = global_df[g_rows]

    if stats:
        for df, fname in zip(
                [matched, df_region, global_df],
                [matched_file, region_file.name, global_file.name]):
            _write_stats(df, fname)

    # exclude station_code for final output files
    column_names.remove(STATION_CODE)
    column_names.remove('SNR')


    global_df[column_names].to_csv(global_file, index=False, header=False,
                                   sep=' ', float_format=FLOAT_FORMAT)

    df_region[column_names].to_csv(region_file, index=False, header=False,
                                   sep=' ', float_format=FLOAT_FORMAT)

    if x_region_df.shape[0]:  # create only if non empty df is returned
        x_region_df[column_names].to_csv(
            cross_region_file, index=False, header=False,
            sep=' ', float_format=FLOAT_FORMAT)


def _write_stats(df, original_file):
    matched_stats_file = os.path.splitext(original_file)[0] + '_stats.csv'
    with open(matched_stats_file, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow([STATION_CODE, STATION_LONGITUDE, STATION_LATITUDE,
                         FREQUENCY])
        for sta, grp in df.groupby(STATION_CODE):
            writer.writerow([sta,
                             grp.iloc[0][STATION_LONGITUDE],
                             grp.iloc[0][STATION_LATITUDE],
                             grp.shape[0]])


def _get_region_string(parameter_file, region):
    if not (parameter_file or region):
        raise ValueError('One of parameter file or region string need to be '
                         'supplied')

    if parameter_file and region:
        log.info('Parameter file will be used for zoning and region string '
                 'will be ignored')

    if parameter_file:
        error_str = 'Check param file format. \n ' \
                    'Here is some help: {}'.format(PARAM_FILE_FORMAT)
        try:
            return _parse_parameter_file(parameter_file)
        except ValueError:
            raise ValueError(error_str)
        except IndexError:
            raise IndexError(error_str)
        except Exception:
            raise Exception('Some unknown parsing error.\n {}'.format(
                error_str))

    else:
        return [float(s) for s in region.split()]


def _parse_parameter_file(param_file):
    with open(param_file, 'r') as f:
        lines = [l.strip() for l in f]

    global_params = [int(float(l)) for l in lines[0].split()]

    local_parms = [float(l) for l in lines[global_params[2] + 2].split()]

    return [local_parms[3], local_parms[2], local_parms[0], local_parms[1]]


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


def _in_region(region, df, grid_size):
    # convert longitude to co-longitude
    df[SOURCE_LONGITUDE] = df[SOURCE_LONGITUDE] % 360
    df[STATION_LONGITUDE] = df[STATION_LONGITUDE] % 360

    if grid_size > 0.0:
        df = _intersect_region(df, region, grid_size)

    # if either the source or the station or both are inside region
    # else global, unless we want a cross region

    # row indices of all in region arrivals
    df_region = df[
        (
            (
                (region.leftlon < df[SOURCE_LONGITUDE]) &
                (df[SOURCE_LONGITUDE] < region.rightlon)
            )
            &
            (
                (region.bottomlat < df[SOURCE_LATITUDE]) &
                (df[SOURCE_LATITUDE] < region.upperlat)
            )
        )
        |
        (
            (
                (region.leftlon < df[STATION_LONGITUDE]) &
                (df[STATION_LONGITUDE] < region.rightlon)
            )
            &
            (
                (region.bottomlat < df[STATION_LATITUDE]) &
                (df[STATION_LATITUDE] < region.upperlat)
            )
        )
        ]

    # dataframe excluding in region arrivals
    df_ex_region = df.iloc[df.index.difference(df_region.index)]

    if grid_size > 0.0:
        # cross region is in ex-region and cross-region==True
        x_region_df = df_ex_region[
            df_ex_region['cross_region'] == True]

        # Global region contain the remaining arrivals
        global_df = df_ex_region[df_ex_region['cross_region'] == False]
        return df_region, global_df, x_region_df
    else:
        global_df = df_ex_region
        return df_region, global_df, pd.DataFrame()


def _intersect_region(df, region, grid_size):  # pragma: no cover
    """
    Strategy to compute cross region: Intersect/cross region is computed first
    which will contain the `region`. The final intersect region will be
    be subtracted from the `region`.
    """

    pe = df[SOURCE_LATITUDE]
    ps = df[STATION_LATITUDE]
    re = df[SOURCE_LONGITUDE]
    rs = df[STATION_LONGITUDE]
    delta = df['locations2degrees']

    # operations on pd.Series
    nms = (delta / grid_size).astype(int)
    ar = pe * DPI
    ast = ps * DPI
    br = re * DPI
    bs = rs * DPI

    x1 = RADIUS * np.sin(ar) * np.cos(br)
    y1 = RADIUS * np.sin(ar) * np.sin(br)
    z1 = RADIUS * np.cos(ar)
    x2 = RADIUS * np.sin(ast) * np.cos(bs)
    y2 = RADIUS * np.sin(ast) * np.sin(bs)
    z2 = RADIUS * np.cos(ast)
    dx = (x2 - x1) / nms
    dy = (y2 - y1) / nms
    dz = (z2 - z1) / nms

    in_cross = []

    # TODO: vectorize this loop
    for i, n in enumerate(nms):
        in_cross.append(_in_cross_region(dx[i], dy[i], dz[i], n, region, x1[i],
                                         y1[i], z1[i]))
    df['cross_region'] = pd.Series(in_cross)
    return df


def _in_cross_region(dx, dy, dz, nms, region, x1, y1, z1):  # pragma: no cover

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


@cli.command()
@click.argument('arrivals_file', type=click.File(mode='r'))
@click.argument('region', type=str,
                metavar="str 'upperlat, bottomlat, leftlon, rightlon'")
def plotcmd(arrivals_file, region):
    """
    Another plot command
    """
    log.info('Begin plotcmd {}'.format(mpiops.rank))
    log.debug("The input arrival file is %s", arrivals_file)

    return

# ================= Quick Testings of the functions ====================
# cd  passive-seismic/
# python seismic/cluster/cluster.py
#
# Usage: cluster.py [OPTIONS] COMMAND [ARGS]...
#
# Options:
#   -v, --verbosity [DEBUG|INFO|WARNING|ERROR]
#                                   Level of logging
#   --help                          Show this message and exit.
#
# Commands:
#   gather  Gather all source-station block pairs for all...
#   match   Match source and station blocks and output...
#   plot    This command will output a...
#   sort    Sort and filter the arrivals.
#   zone    `zone'ing the arrivals into three regions.
# to show help on subcommands:
# $ python cluster/cluster.py gather --help

# How to run this script standalone?
# $ python seismic/cluster/cluster.py plot region_P_1deg.csv '0 -50.0 100 190'
# $ python seismic/cluster/cluster.py -v DEBUG plotcmd region_P_1deg.csv "ab cd"

# $ export ELLIPCORR=/g/data1a/ha3/fxz547/Githubz/passive-seismic/ellip-corr/
# $ python cluster/cluster.py gather -o out /g/data/ha3/events_xmls_test
# $ python ../seismic/cluster/cluster.py gather -o All2dirs /g/data/ha3/fxz547/travel_time_tomography/new_events20180516.run2/events_paths.csv &> ALL_2dirs.log &

# python ../seismic/cluster/cluster.py -v DEBUG gather /g/data/ha3/fxz547/Githubz/passive-seismic/testdata/
# python ../seismic/cluster/cluster.py -v DEBUG gather /g/data/ha3/niket/mtIsa_rf/events_for_fei
# ======================================================================
if __name__ == "__main__":

    cli()
