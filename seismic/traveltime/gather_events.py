"""
Parse multiple events xml files to gather all seismic events-arrivals
(Refactored From the original seismic.cluster.cluster.py)
output CSV files containing the following columns:
['source_block', 'station_block', 'residual', 'event_number',
SOURCE_LONGITUDE, SOURCE_LATITUDE, 'source_depth', STATION_LONGITUDE, STATION_LATITUDE,
'observed_tt', 'calculated_locations2degrees', xml_distance, STATION_CODE, 'SNR', 'P_or_S']

How to Run:
export ELLIPCORR=/g/data1a/ha3/fxz547/Githubz/passive-seismic/ellip-corr/
cd  passive-seismic/tempworks
# python ../seismic/traveltime/gather_events.py  -v DEBUG gather /g/data/ha3/fxz547/Githubz/passive-seismic/testdata/
# python ../seismic/traveltime/gather_events.py  -v DEBUG gather /g/data/ha3/fxz547/Githubz/passive-seismic/some_events_xml/
#  OR
#<fzhang@ubuntu16qos>/Softlab/Githubz/passive-seismic/tempworks (master)
# $ python2 ../seismic/traveltime/gather_events.py  -v DEBUG gather ./events_xmls_test

"""

from __future__ import print_function, absolute_import

import csv
import ellipcorr
import fnmatch
import logging
import os
import random
from math import asin

import click
import pandas as pd
from obspy import read_events
from obspy.geodetics import locations2degrees, gps2dist_azimuth

from inventory.parse_inventory import read_all_stations
from seismic import mpiops
from seismic import pslog

# Only If this gather_events will do computation of the block numbers in a grid
# from seismic.traveltime.cluster_grid import Grid2

DPI = asin(1.0) / 90.0
R2D = 90. / asin(1.)
FLOAT_FORMAT = '%.4f'

log = logging.getLogger(__name__)


# SOURCE_LATITUDE = 'source_latitude'
# SOURCE_LONGITUDE = 'source_longitude'
# STATION_LATITUDE = 'station_latitude'
# STATION_LONGITUDE = 'station_longitude'
# STATION_CODE = 'station_code'
# FREQUENCY = 'no_of_summary_rays'
#
# column_names = ['source_block', 'station_block',
#                 'residual', 'event_number',
#                 SOURCE_LONGITUDE, SOURCE_LATITUDE,
#                 'source_depth', STATION_LONGITUDE, STATION_LATITUDE,
#                 'observed_tt', 'locations2degrees', STATION_CODE, 'SNR', 'P_or_S']
#
# # since we have Basemap in the virtualenv, let's just use that :)
# # ANZ = Basemap(llcrnrlon=100.0, llcrnrlat=-50.0,  urcrnrlon=190.0, urcrnrlat=0.0)
#
# Region = namedtuple('Region', 'upperlat, bottomlat, leftlon, rightlon')


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

    if isinstance(dirname, (list,)):  # the input argument is a list of directory
        filelist = []
        for adir in dirname:
            filelist.extend(recursive_glob(adir))
        return filelist
    else:  # input variable is a single dir
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

        for arow in reader:
            row = arow.strip()  # remove the write spaces \t \n  to have correct linux path

            if row.startswith("#"):
                pass  # comment line
            else:
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
    Gather all source-station arrivals for all events xml files in a list of directory.
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

    # grid definition is needed for clustering purpose, not in gather events
    # grid = Grid(nx=nx, ny=ny, dz=dz)  # original uniform grid
    grid = None  # Grid2()  # use a non-uniform Grid construction

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
    :param grid: can be None; Grid class instance
    :param wave_type: str
        Wave type to generate inversion inputs. See `gather` function.
    :param counter: int
        event counter in this process
    """
    p_type, s_type = wave_type.split()  # ('P', 'S')

    # use preferred origin timestamp as the event number
    # if preferred origin is not populated, use the first origin timestamp
    origin = event.preferred_origin() or event.origins[0]

    # Use the following line when fortran TT code is able to use longer integers
    ev_number = int(origin.time.timestamp)

    # the following definition ensures a 8 digit event number that is also unique
    # delete this definition of ev_number when fortran code can use longer integers
    # assert counter < 100000, 'Counter must be less than 100000'
    # ev_number = int(str(counter) + '{0:0=3d}'.format(mpiops.rank))

    # what columns will be writing out in the gather process??
    col_names = ['source_block', 'station_block', 'residual', 'event_number',
                 'source_longitude', 'source_latitude', 'source_depth',
                 'station_longitude', 'station_latitude', 'observed_tt',
                 'ARRIVAL_TIME', 'ORIGIN_TIME', 'ELLIPTICITY_CORR',  # line for debug
                 'locations2degrees',  'ARRIVAL_DISTANCE','station_code', 'SNR',
                 'P_or_S']

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

    if grid is None:
        event_block = 0
    else:
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
            log.warning("location to degree error %s", e)

        log.debug("location to degree= %s", degrees_to_source)
        # ignore stations more than 90 degrees from source
        if degrees_to_source > 90.0:
            # log.info('Ignored this station arrival as distance from source '
            #          'is {} degrees'.format(degrees_to_source))
            continue

        if grid is None:
            station_block = -1
        else:
            # TODO: use station.elevation information
            station_block = grid.find_block_number(sta.latitude, sta.longitude, z=0.0)

        if arr.phase in wave_type.split():
            log.debug("Began ellipticity_corr ")

            azim_v = gps2dist_azimuth(ev_latitude, ev_longitude, sta.latitude, sta.longitude)[1]

            log.debug("Check input params to ellipticity_corr = %s, %s, %s, %s, %s", arr.phase, degrees_to_source,
                      ev_depth, 90 - ev_latitude, azim_v)

            ellipticity_corr = ellipcorr.ellipticity_corr(
                phase=arr.phase,
                edist=degrees_to_source,
                edepth=ev_depth / 1000.0,
                # TODO: check co-latitude definition
                # no `ecolat` bounds check in fortran ellipcorr subroutine
                # no `origin.latitude` bounds check in obspy
                ecolat=90 - ev_latitude,  # conversion to co-latitude
                azim=azim_v
            )

            log.debug("ellipticity_corr = %s", ellipticity_corr)

            t_list = [event_block, station_block, arr.time_residual, ev_number,
                      ev_longitude, ev_latitude, ev_depth,
                      sta.longitude, sta.latitude,
                      (arr.pick_id.get_referred_object().time.timestamp - origin.time.timestamp) + ellipticity_corr,
                      arr.pick_id.get_referred_object().time, origin.time, ellipticity_corr,  # line for debug
                      degrees_to_source, arr.distance,
                      sta_code, snr_value]
            arrival_staions.append(sta_code)
            p_arrivals.append(t_list + [1]) if arr.phase == p_type else \
                s_arrivals.append(t_list + [2])
        else:  # ignore the other phases
            pass
    return p_arrivals, s_arrivals, missing_stations, arrival_staions


# FZ has moved this function to the grid class
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

    try:
        snr_v = arrival.pick_id.get_referred_object().comments[3]  # Comment(text='snr = 10.7157568852')

        snrlist = str(snr_v).split("snr =")
        snrv = snrlist[-1][:-2]  # the last item of the split, trimming two chars ')
    except Exception as ex:
        log.warning("could not get SNR value becuase %s. Use artificial snrv = 0.0", str(ex))
        snrv = 0.0

    return float(snrv)


# ================= Quick Testings of the functions ====================
# $ export ELLIPCORR=/g/data1a/ha3/fxz547/Githubz/passive-seismic/ellip-corr/

# python ../seismic/traveltime/gather_events.py  [OPTIONS] COMMAND [ARGS]...
#
# Options:
#   -v, --verbosity [DEBUG|INFO|WARNING|ERROR]
#                                   Level of logging
#   --help                          Show this message and exit.
#
# Commands:
#   gather  Gather all source-station block pairs for all...

# to show help on subcommands:
# python seismic/traveltime/gather_events.py  gather --help

# test runs to gather many events xml files in folders
# cd  passive-seismic/tempworks
# python ../seismic/traveltime/gather_events.py  -v DEBUG gather /g/data/ha3/fxz547/Githubz/passive-seismic/testdata/
# python ../seismic/traveltime/gather_events.py  -v DEBUG gather /g/data/ha3/niket/mtIsa_rf/events_for_fei
# python seismic/traveltime/gather_events.py  gather -o outfile /g/data/ha3/events_xmls_test
# python ../seismic/traveltime/gather_events.py  gather -o All2dirs /g/data/ha3/fxz547/travel_time_tomography/new_events20180516.run2/events_paths.csv &> ALL_2dirs.log &

# ======================================================================
if __name__ == "__main__":
    cli()
