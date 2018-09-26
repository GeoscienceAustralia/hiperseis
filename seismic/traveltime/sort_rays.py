"""
Read in event-arrival seismic rays and sort them according to a Grid discretization, and select the best quality data
"""

from __future__ import print_function, absolute_import

import logging
from collections import namedtuple
from math import asin

import click
import pandas as pd

from seismic import pslog
from seismic.traveltime.cluster_grid import Grid, Grid2

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

Region = namedtuple('Region', 'upperlat, bottomlat, leftlon, rightlon')


@click.group()
@click.option('-v', '--verbosity',
              type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR']),
              default='INFO', help='Level of logging')
def cli(verbosity):
    pslog.configure(verbosity)


# @cli.command()
# @click.argument('output_file',
#                 type=click.File(mode='r'))
# @click.argument('residual_cutoff', type=float)
# @click.option('-s', '--sorted_file',
#               type=click.File(mode='w'), default='sorted.csv',
#               help='output sorted and filter file.')
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

    input file header:
    col_names=['source_block', 'station_block', 'residual', 'event_number',
            'source_longitude','source_latitude','source_depth',
            'station_longitude','station_latitude', 'observed_tt', 'locations2degrees', 'station_code','SNR', 'P_or_S']

    :param output_file: output file from the gather stage (eg, outfile_P.csv)
    :param sorted_file: str, optional
        optional sorted output file path. Default: sorted.csv.
    :param residual_cutoff: float
        residual seconds above which arrivals are rejected.
    :return: None
    """

    log.info('Reading in and Filtering arrivals.')

    #cluster_data = pd.read_csv(output_file, header=None,  names=column_names) # original fixed header
    cluster_data = pd.read_csv(output_file,  header='infer') # if input file has correct header line

    cluster_data = cluster_data[abs(cluster_data['residual']) < residual_cutoff]

    # groupby sorts by default
    # cluster_data.sort_values(by=['source_block', 'station_block'],
    #                          inplace=True)

    log.info('Apply a grid model to cluster the rays.')
    # mygrid = Grid(nx=36, ny=18, dz=10000) # use a new grid model
    mygrid = Grid2(ndis=2) # use a new grid model: default ndis=2

    # Re-define the source_block and station_block number according to the mygrid model
    cluster_data['source_block'] = cluster_data.apply(
        lambda x: mygrid.find_block_number(x.source_latitude,x.source_longitude, x.source_depth)[0], axis=1)

    cluster_data['station_block'] = cluster_data.apply(
        lambda x: mygrid.find_block_number(x.station_latitude, x.station_longitude, 0.0)[0], axis=1)


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
    # final_df.drop_duplicates(subset=['source_block', 'station_block',
    #                                  'event_number', SOURCE_LONGITUDE,
    #                                  SOURCE_LATITUDE, 'source_depth'],
    #                          keep='first', inplace=True)
    # FZ: this drop_duplicate will still keep many identical duplicated rows with only event_number different

    # use the following to keep only unique  prim_key: ['source_block', 'station_block']
    #final_df.drop_duplicates(subset=['source_block', 'station_block'],keep='first', inplace=True)
    # keep all station_code. Note that some near-by stations may be cluster into one station_block_number
    final_df.drop_duplicates(subset=['source_block', 'station_block', 'station_code'],keep='first', inplace=True)

    final_df['source_depth'] = final_df['source_depth'] / 1000.0  # scale meter to KM before wrting to csv
    final_df.to_csv(sorted_file, header=False, index=False, sep=' ')

    return final_df


# @cli.command()
# @click.argument('output_file',
#                 type=click.File(mode='r'))
# @click.argument('residual_cutoff', type=float)
# @click.option('-s', '--sorted_file',
#               type=click.File(mode='w'), default='sorted2.csv',
#               help='output sorted and filter file.')
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
    med = cluster_data.groupby(by=['source_block', 'station_block'])[
        'SNR'].max().reset_index()  # use a seq index:0,1,2,.

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

def translate_csv(in_csvfile, out_csvfile):
    """
    Read in a csv file, re-grid each row according to a new Grid model.
    Write into another csv file with re-calculated block_numbers and six new columns of Grid cell centers
    :param in_csvfile: path to an input csv file
    :param out_csvfile: path to an output csv file
    :return: out_csvfile
    """

    log.info('Reading in an input CSV files.')

    incsv = pd.read_csv(in_csvfile, header=None, names=column_names) # tweek if needed according to input file format

    # groupby sorts by default
    # cluster_data.sort_values(by=['source_block', 'station_block'],
    #                          inplace=True)

    log.info('Apply a grid model to cluster the rays.')
    # mygrid = Grid(nx=36, ny=18, dz=10000) # use a new grid model
    mygrid = Grid2() # use a new grid model

    # Re-define the source_block and station_block number according to the mygrid model
    incsv['source_block'] = incsv.apply(
        lambda x: mygrid.find_block_number(x.source_latitude,x.source_longitude, x.source_depth)[0], axis=1)


    incsv['source_xc'] = incsv.apply(lambda x: mygrid.find_block_number(x.source_latitude, x.source_longitude, x.source_depth)[1], axis=1)
    incsv['source_yc'] = incsv.apply(lambda x: mygrid.find_block_number(x.source_latitude, x.source_longitude, x.source_depth)[2], axis=1)
    incsv['source_zc'] = incsv.apply(lambda x: mygrid.find_block_number(x.source_latitude, x.source_longitude, x.source_depth)[3], axis=1)/1000.0 # KM

    incsv['station_block'] = incsv.apply(
        lambda x: mygrid.find_block_number(x.station_latitude, x.station_longitude, 0.0)[0], axis=1)

    incsv['station_xc'] = incsv.apply(lambda x: mygrid.find_block_number(x.station_latitude, x.station_longitude, 0.0)[1], axis=1)
    incsv['station_yc'] = incsv.apply(lambda x: mygrid.find_block_number(x.station_latitude, x.station_longitude, 0.0)[2], axis=1)
    incsv['station_zc'] = incsv.apply(lambda x: mygrid.find_block_number(x.station_latitude, x.station_longitude, 0.0)[3], axis=1)/1000.0

    incsv['source_depth'] = incsv['source_depth'] / 1000.0  # ?? scale meter to KM before wrting to csv

    # write out csv file as you want.
    #incsv.to_csv(out_csvfile, header=False, index=False, sep=' ')
    incsv.to_csv(out_csvfile, header=True, index=False, sep=',')

    return out_csvfile


# ================= Quick Testings of the functions ====================
# cd  passive-seismic/
# export ELLIPCORR=/g/data1a/ha3/fxz547/Githubz/passive-seismic/ellip-corr/

# How to run this script standalone?
#$ fxz547@vdi-n2 /g/data/ha3/fxz547/travel_time_tomography/CSV_NewFormatAug10/FZ01-pst-cluster2/testrun3
#$ python /g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/sort_rays.py sort2 S_out.csv 10. -s sorted_S.csv
# ======================================================================
if __name__ == "__main__":
    # cli()   # This is click API command entry point

    import sys
    inf = sys.argv[1]
    outf= sys.argv[2]
    translate_csv(inf, outf)
