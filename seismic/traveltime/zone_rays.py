from __future__ import print_function, absolute_import

import csv
import logging
import os
from collections import namedtuple
from math import asin, sin, acos, sqrt

import click
import numpy as np
import pandas as pd
from obspy.geodetics.base import WGS84_A as RADIUS

from seismic.traveltime import pslog

DPI = asin(1.0) / 90.0
R2D = 90. / asin(1.)
FLOAT_FORMAT = '%.4f'

logging.basicConfig()
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


Region = namedtuple('Region', 'upperlat, bottomlat, leftlon, rightlon')


@click.group()
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cli(verbosity):
    pslog.configure(verbosity)


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
                    'Here is some help: {}'.format("PARAM_FILE_FORMAT")
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



# ================= Quick Testings of the functions ====================
# cd  passive-seismic/
# export ELLIPCORR=/g/data1a/ha3/fxz547/Githubz/passive-seismic/ellip-corr/

# How to run this script standalone?
#$ fxz547@vdi-n2 /g/data/ha3/fxz547/travel_time_tomography/CSV_NewFormatAug10/FZ01-pst-cluster2/testrun3
#$ python /g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/zone_rays.py  -v DEBUG zone sorted_S.csv -z '0 -50.0 100 190' -r sorted_region_S.csv -g sorted_global_S.csv
# ======================================================================
if __name__ == "__main__":
    cli()
