"""
Description:
    Read in event-arrival seismic rays and sort them according to a discretization model of Earth,
    The output will feed into an inversion program.

Developer:
    fei.zhang@ga.gov.au
"""
from __future__ import print_function, absolute_import

import json
import logging
import sys
from collections import namedtuple
from math import asin

import click
import ellipcorr
import pandas as pd
from obspy.geodetics import gps2dist_azimuth, locations2degrees

from seismic.traveltime import pslog
from seismic.traveltime.cluster_grid import Grid2

DPI = asin(1.0) / 90.0
R2D = 90. / asin(1.)
FLOAT_FORMAT = '%.4f'

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

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
    """
    CLI group logging config
    :param verbosity:
    :return:
    """
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
    col_names=['source_block', 'station_block', 'residual', 'event_number','source_longitude','source_latitude','source_depth',
            'station_longitude','station_latitude', 'observed_tt', 'locations2degrees', 'station_code','SNR', 'P_or_S']

    :param output_file: output file from the gather stage (eg, outfile_P.csv)
    :param sorted_file: str, optional
        optional sorted output file path. Default: sorted.csv.
    :param residual_cutoff: float residual seconds above which arrivals are rejected.
    :return: None
    """

    log.info('Reading in and Filtering arrivals.')

    # cluster_data = pd.read_csv(output_file, header=None,  names=column_names) # original fixed header
    cluster_data = pd.read_csv(output_file, header='infer')  # if input file has correct header line

    cluster_data = cluster_data[abs(cluster_data['residual']) < residual_cutoff]

    # groupby sorts by default
    # cluster_data.sort_values(by=['source_block', 'station_block'],
    #                          inplace=True)

    log.info('Apply a grid model to cluster the rays.')
    # mygrid = Grid(nx=36, ny=18, dz=10000) # use a new grid model
    mygrid = Grid2(ndis=2)  # use a new grid model: default ndis=2

    # Re-define the source_block and station_block number according to the mygrid model
    cluster_data['source_block'] = cluster_data.apply(
        lambda x: mygrid.find_block_number(x.source_latitude, x.source_longitude, x.source_depth)[0], axis=1)

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
    # final_df.drop_duplicates(subset=['source_block', 'station_block'],keep='first', inplace=True)
    # keep all station_code. Note that some near-by stations may be cluster into one station_block_number
    # final_df.drop_duplicates(subset=['source_block', 'station_block', 'station_code'],keep='first', inplace=True)
    final_df.drop_duplicates(subset=['source_block', 'station_block'], keep='first', inplace=True)

    # make sure the originDepth is in KM
    final_df['source_depth'] = final_df['source_depth'] / 1000.0  # convert meter into KM before writing output csv
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

    cluster_data['source_depth'] = cluster_data['source_depth'] / 1000.0  # convert to KM for output file

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

    incsv = pd.read_csv(in_csvfile, header=None, names=column_names)  # tweek if needed according to input file format

    # groupby sorts by default
    # cluster_data.sort_values(by=['source_block', 'station_block'],
    #                          inplace=True)

    log.info('Apply a grid model to cluster the rays.')
    # mygrid = Grid(nx=36, ny=18, dz=10000) # use a new grid model
    mygrid = Grid2()  # use a new grid model

    # Re-define the source_block and station_block number according to the mygrid model
    incsv['source_block'] = incsv.apply(
        lambda x: mygrid.find_block_number(x.source_latitude, x.source_longitude, x.source_depth)[0], axis=1)

    incsv['source_xc'] = incsv.apply(
        lambda x: mygrid.find_block_number(x.source_latitude, x.source_longitude, x.source_depth)[1], axis=1)
    incsv['source_yc'] = incsv.apply(
        lambda x: mygrid.find_block_number(x.source_latitude, x.source_longitude, x.source_depth)[2], axis=1)
    incsv['source_zc'] = incsv.apply(
        lambda x: mygrid.find_block_number(x.source_latitude, x.source_longitude, x.source_depth)[3],
        axis=1) / 1000.0  # KM

    incsv['station_block'] = incsv.apply(
        lambda x: mygrid.find_block_number(x.station_latitude, x.station_longitude, 0.0)[0], axis=1)

    incsv['station_xc'] = incsv.apply(
        lambda x: mygrid.find_block_number(x.station_latitude, x.station_longitude, 0.0)[1], axis=1)
    incsv['station_yc'] = incsv.apply(
        lambda x: mygrid.find_block_number(x.station_latitude, x.station_longitude, 0.0)[2], axis=1)
    incsv['station_zc'] = incsv.apply(
        lambda x: mygrid.find_block_number(x.station_latitude, x.station_longitude, 0.0)[3], axis=1) / 1000.0

    incsv['source_depth'] = incsv['source_depth'] / 1000.0  # ?? scale meter to KM before wrting to csv

    # write out csv file as you want.
    # incsv.to_csv(out_csvfile, header=False, index=False, sep=' ')
    incsv.to_csv(out_csvfile, header=True, index=False, sep=',')

    return out_csvfile


def compute_ellipticity_corr(arrival_phase, ev_latitude, ev_longitude, ev_depth_km, sta_latitude, sta_longitude,
                             degrees_to_source):
    myazim = gps2dist_azimuth(ev_latitude, ev_longitude, sta_latitude, sta_longitude)[1] # [1] shall be taken
    # see https://docs.obspy.org/_modules/obspy/geodetics/base.html#gps2dist_azimuth
    # this function returns 3 values (Great_circle_distance_in_m, azimuth_A->B_in_degrees, azimuth_B->A_in degrees)

    log.debug("Check input params to ellipticity_corr = %s, %s, %s, %s, %s", arrival_phase, degrees_to_source,
              ev_depth_km, 90 - ev_latitude, myazim)

    ellipticity_corr = ellipcorr.ellipticity_corr(
        phase=arrival_phase,
        edist=degrees_to_source,
        edepth=ev_depth_km,
        ecolat=90 - ev_latitude,  # conversion to co-latitude
        azim=myazim
    )

    log.debug("ellipticity_corr = %s", ellipticity_corr)

    return ellipticity_corr


def sort_csv_in_grid(inputcsv, outputcsv, phase, mygrid, column_name_map):
    """ Read in a csv file, re-grid each row according to a given Grid model.
    Write into output csv file with re-calculated block_numbers re-named columns

    :param inputcsv: path to an input csv file
    :param outputcsv:
    :return:
    """

    if phase.upper() == 'P':
        residual_cutoff = 5.0
    elif phase.upper() == 'S':
        residual_cutoff = 10.0
    else:
        print("!!! Phase must be either P or S !!!")
        raise Exception("Phase < %s >  Not Recognised! " % (phase))

    log.info('Reading in the input CSV files.')

    csv_data = pd.read_csv(inputcsv, sep='\s+', header='infer')
    # , names=column_names)  # tweek if needed according to input file format

    print(csv_data.head(3))

    log.info("renaming columns according to the mapping dict %s" % str(column_name_map))

    csv_data.rename(columns=column_name_map, inplace=True)

    # Thing below is based on the assumption that the csv_data pdf has the following columns:
    #  ['source_block', 'station_block',  'tt_residual', 'event_number',   'source_lon', 'source_lat',
    # 'source_depth_km', 'station_lon',  'station_lat',  'observed_tt', 'locations2degrees',  'station_code',  'p_or_s']

    log.info('Select useful rows/rays by applying filters.')

    csv_data = csv_data[abs(csv_data['tt_residual']) < residual_cutoff]
    # remove shallow event originDepth<0.01KM?

    log.info('Apply a grid model to discretize the rays, before sorting clustering.')

    # # Re-define the source_block and station_block number according to the mygrid model
    # originLon
    # originLat
    # originDepthKm
    # net...
    # stationLon
    # stationLat
    csv_data['source_block'] = csv_data.apply(
        lambda x: mygrid.find_block_number(x.source_lat, x.source_lon, 1000 * (x.source_depth_km))[0], axis=1)

    csv_data['station_block'] = csv_data.apply(
        lambda x: mygrid.find_block_number(x.station_lat, x.station_lon, 0.0)[0], axis=1)

    csv_data['observed_tt'] = csv_data.pickTimestamp - csv_data.originTimestamp

    # cluster_data.to_csv(outputcsv+"_debug.CSV", header=True, index=False, sep=',')

    log.info('Begin Sorting arrivals.')

    # groupby automatically sorts
    med = csv_data.groupby(by=['source_block', 'station_block'])[
        'observed_tt'].quantile(q=.5, interpolation='lower').reset_index()

    final_df = pd.merge(csv_data, med, how='right',
                        on=['source_block', 'station_block', 'observed_tt'],
                        sort=True, right_index=True)

    # use the following to keep only unique  prim_key: ['source_block', 'station_block']
    # final_df.drop_duplicates(subset=['source_block', 'station_block'],keep='first', inplace=True)
    # Note that some near-by stations may be cluster into one station_block_number, if want to keep stations try
    # final_df.drop_duplicates(subset=['source_block', 'station_block', 'station_code'],keep='first', inplace=True)

    final_df.drop_duplicates(subset=['source_block', 'station_block'], keep='first', inplace=True)

    # elliptic correction to the  observed_travel_time;

    final_df['locations_to_degrees'] = final_df.apply(lambda x: locations2degrees( x.source_lat, x.source_lon,
                                                                x.station_lat, x.station_lon), axis=1)
    final_df['my_azim'] = final_df.apply(lambda x: gps2dist_azimuth( x.source_lat, x.source_lon,
                                                                x.station_lat, x.station_lon)[1], axis=1)

    final_df['my_bazim'] = final_df.apply(lambda x: gps2dist_azimuth(x.source_lat, x.source_lon,
                                                                    x.station_lat, x.station_lon)[2], axis=1)
    final_df['ellipticity_corr'] = final_df.apply(lambda x:
                                          compute_ellipticity_corr(phase, x.source_lat, x.source_lon,x.source_depth_km,
                                                                   x.station_lat, x.station_lon, x.distance),
                                          axis=1)

    final_df['observed_tt'] = final_df.observed_tt + final_df.ellipticity_corr




    # make sure the originDepth/source_depth is in KM for required by inversion program

    final_df.to_csv(outputcsv, header=True, index=False, sep=',')  # use comma separator,

    # inpdf.to_csv(outputcsv, header=True, index=False, sep=' ')   # mismatch columns in space-delimited csv file as the NaN => empty space !

    if phase =='P':
        final_df['P_or_S'] =1
    elif phase =='S':
        final_df['P_or_S'] = 2
    else:
        raise Exception("Phase must be P or S !!!")

    final_df['event_number'] = final_df.apply(lambda x: int(x.originTimestamp), axis=1)

    # the following values are required for inversion program. the event_number defined as int(originTimestamp)
    # the columns must be in the order:
    required_columns = ['source_block', 'station_block', 'tt_residual', 'event_number',
                 'source_lon', 'source_lat', 'source_depth_km',
                 'station_lon', 'station_lat', 'observed_tt', 'locations_to_degrees',  'P_or_S']


    pdf4inv = final_df[required_columns]

    inv_txt = "%s_inv.txt" % outputcsv
    pdf4inv.to_csv(inv_txt,header=False, index=False, sep=' ', float_format='%.6f')  # space delimitted txt file

    return outputcsv


def get_columns_dict(jsonfile):
    with open(jsonfile) as json_file:
        col_names_dict = json.load(json_file)

    json_file_dump = jsonfile + ".dump"

    with open(json_file_dump, 'w') as fp:
        json.dump(col_names_dict, fp, indent=4)

    return col_names_dict


# ================= Quick Testings of the functions ====================
# cd  passive-seismic/
# export ELLIPCORR=/g/data1a/ha3/fxz547/Githubz/passive-seismic/ellip-corr/

# How to run this script standalone?
# python /g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/sort_rays.py
# /g/data/ha3/rakib/seismic/pst/tests/results/p_arrivals_mag_4_and_above.txt p_arrivals_mag_4_and_above_sorted1x1.csv P
# $  python /g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/sort_rays.py
# /g/data/ha3/rakib/seismic/pst/tests/results/s_arrivals_mag_4_and_above.txt s_arrivals_mag_4_and_above_sorted1x1.csv S parfile column_json
# $ fxz547@vdi-n2 /g/data/ha3/fxz547/travel_time_tomography/CSV_NewFormatAug10/FZ01-pst-cluster2/testrun3
# $ python /g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/sort_rays.py sort2 S_out.csv 10. -s sorted_S.csv
# ======================================================================
if __name__ == "__main__":
    # cli()   # This is click API command entry point

    inf = sys.argv[1]
    outf = sys.argv[2]

    phase = sys.argv[3]  # P or S => residual_cutoff=5 for Pwave. residual_cutoff=10 for Swave
    in_param_file = sys.argv[4]  # '/g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/param1x1'

    if len(sys.argv) > 5:
        col_json_file = sys.argv[5]
        columns_dict = get_columns_dict(col_json_file)

    # define a grid for clustering the rays
    # mygrid = Grid2(param_file='/g/data/ha3/fxz547/Githubz/passive-seismic/seismic/traveltime/param2x2')
    mygrid = Grid2(param_file=in_param_file)

    sort_csv_in_grid(inf, outf, phase, mygrid,
                     columns_dict)  # residual_cutoff=5 for Pwave. residual_cutoff=10 for Swave
