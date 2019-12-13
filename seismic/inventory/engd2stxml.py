#!/usr/bin/env python
"""
==================================================
Engdahl and ISC STN data conversion to FDSN station XML
==================================================

Creates database of stations from .STN files,
curates the data using heuristic rules, and exports new stations to FDSN
station XML format network with non-empty, nominal instrument response data.

Cleanup steps applied:
* Removes "blacklisted" networks that add little value and cause problems due to station code conflicts.
* Add default station dates where missing.
* Make "future" station end dates consistent to max Pandas timestamp.
* Remove records with illegal station codes.
* Remove duplicate station records.
* Merge overlapping channel dates for given NET.STAT.CHAN to a single epoch.
"""

# pylint: disable=too-many-locals, invalid-name

from __future__ import division

import os
import sys
import argparse
import datetime

from ordered_set import OrderedSet as set
import numpy as np
import scipy as sp
import pandas as pd
import requests

import obspy
from obspy import read_inventory
from seismic.inventory.pdconvert import dataframe_to_fdsn_station_xml
from seismic.inventory.table_format import (TABLE_SCHEMA, TABLE_COLUMNS, PANDAS_MAX_TIMESTAMP,
                                            DEFAULT_START_TIMESTAMP, DEFAULT_END_TIMESTAMP)
from seismic.inventory.inventory_util import NOMINAL_EARTH_RADIUS_KM, SORT_ORDERING, extract_unique_sensors_responses

PY2 = sys.version_info[0] < 3

if PY2:
    import cStringIO as sio  # pylint: disable=import-error
    import cPickle as pkl  # pylint: disable=import-error
else:
    import io as sio  # pylint: disable=ungrouped-imports
    import pickle as pkl

print("Using Python version {0}.{1}.{2}".format(*sys.version_info))
print("Using obspy version {}".format(obspy.__version__))

try:
    import tqdm
    show_progress = True
except ImportError:
    show_progress = False
    print("Run 'pip install tqdm' to see progress bar.")

# Script requires numpy >= 1.15.4. The source of the error is not yet identified, but has been
# demonstrated on multiple platforms with lower versions of numpy.
vparts = np.version.version.split('.', 2)
(major, minor, maint) = [int(x) for x in vparts]
if major < 1 or (major == 1 and minor < 14):  # pragma: no cover
    print("Not supported error: Requires numpy >= 1.14.2, found numpy {0}".format(".".join(vparts)))
    sys.exit(1)
else:
    print("Using numpy {0}".format(".".join(vparts)))

# Pandas table display options to reduce aggressiveness of truncation. Due to size of data sometimes we
# need to see more details in the table.
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', -1)
pd.set_option('display.width', 240)

# Global constants. Assumption of spherical earth model. This is quite a weak assumption, but
# the obspy function locations2degrees() uses spherical earth model too, so no difference to
# distance method used here.
DIST_TOLERANCE_KM = 2.0
DIST_TOLERANCE_RAD = DIST_TOLERANCE_KM / NOMINAL_EARTH_RADIUS_KM
COSINE_DIST_TOLERANCE = np.cos(DIST_TOLERANCE_RAD)

# List of networks to remove outright. See ticket PST-340.
# CI: too many station code conflicts with other global networks
# 7B, 7D, 7F, 7G, 7W, 7X: these are Australian deployments and GA has more comprehensive and
#     accurate data on these networks than is available in the .STN files.
BLACKLISTED_NETWORKS = ("CI", "7B", "7D", "7F", "7G", "7W", "7X")

# Timestamp to be added to output file names, so that each run generates unique log files.
rt_timestamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")


def read_eng(fname):
    """
    Read Engdahl STN file having the following format of fixed width formatted columns:

    :: AAI   Ambon             BMG, Indonesia, IA-Ne              -3.6870  128.1945      0.0   2005001  2286324  I
    :: AAII                                                       -3.6871  128.1940      0.0   2005001  2286324  I
    :: AAK   Ala Archa         Kyrgyzstan                         42.6390   74.4940      0.0   2005001  2286324  I
    :: ABJI                                                       -7.7957  114.2342      0.0   2005001  2286324  I
    :: APSI                                                       -0.9108  121.6487      0.0   2005001  2286324  I
    :: AS01  Alice Springs Arra                                  -23.6647  133.9508      0.0   2005001  2286324  I
    :: 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    ::          10        20        30        40        50        60        70        80

    Each Station Code (first column) might NOT be unique, and network codes are missing here, so all records are
    placed under 'GE' network.

    :param fname: STN file name to load
    :type fname: str
    :return: Pandas Dataframe containing the loaded data in column order of TABLE_COLUMNS.
    :rtype: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    """
    # Intervals of column numbers delineating input fields.
    colspec = ((0, 6), (59, 67), (68, 77), (78, 86))
    col_names = ['StationCode', 'Latitude', 'Longitude', 'Elevation']
    data_frame = pd.read_fwf(fname, colspecs=colspec, names=col_names, dtype=TABLE_SCHEMA)
    # Assumed network code for files of this format.
    data_frame['NetworkCode'] = 'GE'
    # Populate missing data.
    data_frame['StationStart'] = pd.NaT
    data_frame['StationEnd'] = pd.NaT
    data_frame['ChannelStart'] = pd.NaT
    data_frame['ChannelEnd'] = pd.NaT
    # Default channel code.
    data_frame['ChannelCode'] = 'BHZ'
    # Sort columns into preferred order
    data_frame = data_frame[list(TABLE_COLUMNS)]
    # Compute and report number of duplicates
    num_dupes = len(data_frame) - len(data_frame['StationCode'].unique())
    print("{0}: {1} stations found with {2} duplicates".format(fname, len(data_frame), num_dupes))
    return data_frame


def reportUnpickleFail(filename):  # pragma: no cover
    """
    Standard failure report message when trying to unpickle file.

    :param filename: The name of the file that failed to unpickle.
    :type filename: str
    """
    print("PKL LOAD FAILED: {} file incompatible or corrupt, please delete. "
          "Falling back to full parse.".format(filename))


def read_isc(fname, use_pickle=False):
    """
    Read ISC station inventory having such format and convert to Pandas DataFrame:

    :: 109C     32.8892 -117.1100     0.0 2006-06-01 04:11:18 2008-01-04 01:26:30
    :: 109C     32.8882 -117.1050   150.0 2008-01-04 01:26:30
    ::              FDSN 109C   TA -- BHZ 2004-05-04 23:00:00 2005-03-03 23:59:59
    ::              FDSN 109C   TA -- LHZ 2004-05-04 23:00:00 2005-03-03 23:59:59
    ::              FDSN 109C   TA -- BHZ 2005-04-11 00:00:00 2006-01-25 22:31:10
    ::              FDSN 109C   TA -- LHZ 2005-04-11 00:00:00 2006-01-25 22:31:10
    :: 0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
    ::           10        20        30        40        50        60        70        80

    The lines starting with a station code are HEADER rows, and provide the station coordinates.
    The idented lines starting with FDSN provide distinct station, network and channel data for the
    given station location.

    :param fname: STN file name to load
    :type fname: str
    :return: Pandas Dataframe containing the loaded data in column order of TABLE_COLUMNS.
    :rtype: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    """
    # Intervals of column numbers delineating input fields.
    header_colspec = ((0, 5), (7, 16), (17, 26), (27, 34), (35, 54,), (55, 74))
    header_cols = ['StationCode', 'Latitude', 'Longitude', 'Elevation', 'StationStart', 'StationEnd']
    channel_colspec = ((13, 17), (18, 23), (24, 27), (31, 34), (35, 54), (55, 74))
    channel_cols = ['FDSN', 'StationCode', 'NetworkCode', 'ChannelCode', 'ChannelStart', 'ChannelEnd']

    # Timestamps in source data which present far future. Will be replaced by max supported Pandas timestamp.
    ISC_INVALID_TIMESTAMPS = ["2500-01-01 00:00:00",
                              "2500-12-31 23:59:59",
                              "2599-01-01 00:00:00",
                              "2599-12-31 00:00:00",
                              "2599-12-31 23:59:59",
                              "2999-12-31 23:59:59",
                              "5000-01-01 00:00:00"]

    # Nested helper function
    def reportStationCount(df):
        """
        Convenience function to report on number of unique network and station codes in the dataframe.

        :param df: Dataframe to report on
        :type df: pandas.DataFrame
        """
        num_unique_networks = len(df['NetworkCode'].unique())
        num_unique_stations = len(df['StationCode'].unique())
        print("{0}: {1} unique network codes, {2} unique station codes".format(fname, num_unique_networks,
                                                                               num_unique_stations))

    if use_pickle:  # pragma: no cover
        pkl_name = fname + ".pkl"
        if os.path.exists(pkl_name):
            print("Reading cached " + fname)
            try:
                with open(pkl_name, 'rb') as f:
                    df_all = pkl.load(f)
                    reportStationCount(df_all)
                    return df_all
            except Exception:
                reportUnpickleFail(pkl_name)

    print("Parsing " + fname)
    df_list = []
    # Due to irregular data format, read one row at a time from file. Consider replacing with a pre-processing pass
    # so that more than one row can be read at a time.
    if show_progress:
        pbar = tqdm.tqdm(total=os.path.getsize(fname), ascii=True)
    with open(fname, "r", buffering=64 * 1024) as f:
        line = f.readline()
        if show_progress:
            pbar.update(len(line))
        hdr = None
        while line.strip():
            channels = []
            while 'FDSN' in line:
                # Read channel rows.
                # Substitute max timestamp from source data with the lower value limited by Pandas.
                for ts_unsupported in ISC_INVALID_TIMESTAMPS:
                    line = line.replace(ts_unsupported, PANDAS_MAX_TIMESTAMP)
                line_input = sio.StringIO(line)
                ch_data = pd.read_fwf(line_input, colspecs=channel_colspec, names=channel_cols, nrows=1,
                                      dtype=TABLE_SCHEMA, na_filter=False, parse_dates=[4, 5])
                assert ch_data.iloc[0]['FDSN'] == 'FDSN'
                channels.append(ch_data)
                line = f.readline()
                if show_progress:
                    pbar.update(len(line))
            if hdr is not None:
                # Always store header data as a station record.
                hdr['NetworkCode'] = 'IR'     # pylint: disable=unsupported-assignment-operation
                hdr['ChannelStart'] = pd.NaT  # pylint: disable=unsupported-assignment-operation
                hdr['ChannelEnd'] = pd.NaT    # pylint: disable=unsupported-assignment-operation
                hdr['ChannelCode'] = 'BHZ'    # pylint: disable=unsupported-assignment-operation
                # Standardize column ordering
                hdr = hdr[list(TABLE_COLUMNS)]  # pylint: disable=unsubscriptable-object
                if channels:
                    # If channel data is also present, store it too.
                    ch_all = pd.concat(channels, sort=False)
                    ch_all.drop('FDSN', axis=1, inplace=True)
                    # Set the station date range to at least encompass the channels it contains.
                    st_min = np.array([hdr['StationStart'].min(), ch_all['ChannelStart'].min()]).min()
                    st_max = np.array([hdr['StationEnd'].max(), ch_all['ChannelEnd'].max()]).max()
                    hdr['StationStart'] = st_min
                    hdr['StationEnd'] = st_max
                    # Assign common fields to the channel rows.
                    ch_all[['Latitude', 'Longitude', 'Elevation', 'StationStart', 'StationEnd']] = \
                        hdr[['Latitude', 'Longitude', 'Elevation', 'StationStart', 'StationEnd']]
                    # Make sure column ordering is consistent
                    network_df = ch_all[list(TABLE_COLUMNS)]
                    df_list.append(network_df)
                df_list.append(hdr)
                hdr = None
            # Read header row
            if line:
                line_input = sio.StringIO(line)
                hdr = pd.read_fwf(line_input, colspecs=header_colspec, names=header_cols, nrows=1, dtype=TABLE_SCHEMA,
                                  parse_dates=[4, 5])
                line = f.readline()
            else:
                hdr = None
            if show_progress:
                pbar.update(len(line))
    if show_progress:
        pbar.close()
    print("Concatenating records...")
    df_all = pd.concat(df_list, sort=False)
    if use_pickle:  # pragma: no cover
        with open(fname + ".pkl", "wb") as f:
            pkl.dump(df_all, f, pkl.HIGHEST_PROTOCOL)
    reportStationCount(df_all)
    return df_all


def remove_blacklisted(df):
    """
    Remove network codes that are explicitly blacklisted due to QA issues or undesirable overlap
    with trusted FDSN station codes.

    :param df: Dataframe of initially loaded data from STN files
    :type df: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    """
    for badnet in BLACKLISTED_NETWORKS:
        df = df[df["NetworkCode"] != badnet]
    return df


def remove_illegal_stationNames(df):
    """
    Remove records for station names that do not conform to expected naming convention.
    Such names can cause problems in downstream station, in particular names with asterisk.

    :param df: Dataframe containing station records from which illegal station codes should be
        removed (modified in-place)
    :type df: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    """
    import re
    pattern = re.compile(r"^[a-zA-Z0-9]{1}[\w\-]{1,4}$")
    removal_index = []
    for (netcode, statcode), data in df.groupby(['NetworkCode', 'StationCode']):
        # assert isinstance(statcode, str)
        if not pattern.match(statcode):
            print("UNSUPPORTED Station Code: {0}.{1}".format(netcode, statcode))
            removal_index.extend(data.index.tolist())

    if removal_index:
        df.drop(removal_index, inplace=True)


def latlong_to_cosinedistance(latlong_deg_set1, latlong_deg_set2):
    """
    Compute the approximate cosine distance between each station of 2 sets.

    Each set is specified as a numpy column vector of [latitude, longitude] positions in degrees.

    This function performs an outer product and will produce matrix of size N0 x N1, where
    N0 is the number of rows in latlong_deg_set1 and N1 is the number of rows in latlong_deg_set2.

    Returns np.ndarray containing cosines of angles between each pair of stations from
    the input arguments.
    If input is 1D, convert to 2D for consistency of matrix orientations.

    :param latlong_deg_set1: First set of numpy column vector of [latitude, longitude] positions in
        degrees
    :type latlong_deg_set1: np.ndarray
    :param latlong_deg_set2: Second set of numpy column vector of [latitude, longitude] positions in
        degrees
    :type latlong_deg_set2: np.ndarray
    :return: Array containing cosines of angles between each pair of stations from the input
        arguments.
    :rtype: np.ndarray
    """
    if len(latlong_deg_set1.shape) == 1:
        latlong_deg_set1 = np.reshape(latlong_deg_set1, (1, -1))
    if len(latlong_deg_set2.shape) == 1:
        latlong_deg_set2 = np.reshape(latlong_deg_set2, (1, -1))

    set1_latlong_rad = np.deg2rad(latlong_deg_set1)  # pylint: disable=assignment-from-no-return
    set2_latlong_rad = np.deg2rad(latlong_deg_set2)  # pylint: disable=assignment-from-no-return

    set1_polar = np.column_stack((
        np.sin(set1_latlong_rad[:, 0]) * np.cos(set1_latlong_rad[:, 1]),
        np.sin(set1_latlong_rad[:, 0]) * np.sin(set1_latlong_rad[:, 1]),
        np.cos(set1_latlong_rad[:, 0])))

    set2_polar = np.column_stack((
        np.sin(set2_latlong_rad[:, 0]) * np.cos(set2_latlong_rad[:, 1]),
        np.sin(set2_latlong_rad[:, 0]) * np.sin(set2_latlong_rad[:, 1]),
        np.cos(set2_latlong_rad[:, 0]))).T

    cosine_dist = np.dot(set1_polar, set2_polar)

    # Collapse result to minimum number of dimensions necessary
    result = np.squeeze(cosine_dist)
    if np.isscalar(result):
        result = np.array([result], ndmin=1)
    return result


def compute_neighboring_station_matrix(df):
    """
    Compute sparse matrix representing index of neighboring stations.

    Ordering of matrix corresponds to ordering of Dataframe df, which is expected to be sequential
    integer indexed. For a given station index i, then the non-zero off-diagonal entries in row i
    of the returned matrix indicate the indices of adjacent, nearby stations.

    :param df: Dataframe containing station records.
    :type df: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    :return: Sparse binary matrix having non-zero values at indices of neighboring stations.
    :rtype: scipy.sparse.csr_matrix
    """
    self_latlong = df[["Latitude", "Longitude"]].values
    # In order to keep calculation tractable for potentially large matrices without resort to
    # on-disk memmapped arrays, we split the second operand into parts and compute parts of the
    # result, then recombine them to get the final (sparse) result.
    sparse_cos_dist = []
    num_splits = max(self_latlong.shape[0] // 2000, 1)
    for m in np.array_split(self_latlong, num_splits):
        partial_result = latlong_to_cosinedistance(self_latlong, m)
        partial_result = sp.sparse.csr_matrix(partial_result >= COSINE_DIST_TOLERANCE)
        sparse_cos_dist.append(partial_result)
    cos_dist = sp.sparse.hstack(sparse_cos_dist, "csr")
    return cos_dist


def remove_duplicate_stations(df, neighbor_matrix):
    """
    Remove stations which are identified as duplicates:
    * Removes duplicated stations in df based on station code and locality of lat/long coordinates.
    * Removes duplicated stations based on codes and channel data matching, IRRESPECTIVE of locality

    :param df: Dataframe containing station records. Is modified during processing.
    :type df: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    :param neighbor_matrix: Sparse binary matrix having non-zero values at indices of neighboring
        stations.
    :type neighbor_matrix: scipy.sparse.csr_matrix
    :return: Dataframe containing station records with identified duplicates removed.
    :rtype: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    """
    assert len(df) == neighbor_matrix.shape[0]
    assert neighbor_matrix.shape[0] == neighbor_matrix.shape[1]
    num_stations = len(df)
    # Firstly, remove stations by nearness to other stations with matching codes and channel data
    removal_rows = set()
    matching_criteria = ["NetworkCode", "StationCode", "ChannelCode", "ChannelStart", "ChannelEnd"]
    print("  LOCATION duplicates...")
    with open("LOG_LOCATION_DUPES_" + rt_timestamp + ".txt", 'w') as log:
        if show_progress:
            pbar = tqdm.tqdm(total=len(df), ascii=True)
        for i in range(num_stations):
            if show_progress:
                pbar.update()
            if i in removal_rows:
                continue
            row = neighbor_matrix.getrow(i)
            neighbors = row.nonzero()[1]
            # Only consider upper diagonal so that we don't doubly register duplicates
            neighbors = neighbors[neighbors > i]
            if len(neighbors) < 1:
                continue
            key = df.loc[i, matching_criteria]
            # Check which of the nearby stations match network and station code. We only remove rows if these match,
            # otherwise just raise warning.
            attrs_match = np.array([((k[1] == key) | (k[1].isna() & key.isna()))
                                    for k in df.loc[neighbors, matching_criteria].iterrows()])
            duplicate_mask = np.all(attrs_match, axis=1)
            if np.any(duplicate_mask):
                duplicate_index = neighbors[duplicate_mask]
                log.write("WARNING: Duplicates of\n{0}\nare "
                          "being removed:\n{1}\n----\n".format(df.loc[[i]], df.loc[duplicate_index]))
                removal_rows.update(duplicate_index)
        if show_progress:
            pbar.close()
    removal_rows = np.array(sorted(list(removal_rows)))
    if removal_rows.size > 0:
        print("Removing following {0} duplicates due to identical network, station and "
              "channel data:\n{1}".format(len(removal_rows), df.loc[removal_rows]))
        df.drop(removal_rows, inplace=True)

    # Secondly, remove stations with same network and station code, but which are further away than the
    # threshold distance and have no distinguishing channel data. We deliberately exclude station start
    # and end dates from consideration here, as these are extended during file read to cover range of
    # contained channels, and therefore might not match in code dupes.
    matching_criteria = ["ChannelCode", "ChannelStart", "ChannelEnd"]
    removal_index = set()
    print("  CODE duplicates...")
    with open("LOG_CODE_DUPES_" + rt_timestamp + ".txt", 'w') as log:
        if show_progress:
            pbar = tqdm.tqdm(total=len(df), ascii=True)
        for _, data in df.groupby(['NetworkCode', 'StationCode']):
            if show_progress:
                pbar.update(len(data))
            if len(data) <= 1:
                continue
            for row_index, channel in data.iterrows():
                if row_index in removal_index:
                    continue
                key = channel[matching_criteria]
                # Consider a likely duplicate if all matching criteria are same as the key.
                # Note that NA fields will compare False even if both are NA, which is what we want here
                # since we don't want to treat records with same codes as duplicates if the matching_criteria
                # are NA, as this removes records that are obviously not duplicates.
                duplicate_mask = (data[matching_criteria] == key)
                index_mask = np.all(duplicate_mask, axis=1) & (data.index > row_index)
                duplicate_index = data.index[index_mask]
                if not duplicate_index.empty:
                    log.write("WARNING: Apparent duplicates of\n{0}\nare being removed:\n{1}\n"
                              "----\n".format(data.loc[[row_index]], data.loc[duplicate_index]))
                    removal_index.update(duplicate_index.tolist())
        if show_progress:
            pbar.close()
    removal_index = np.array(sorted(list(removal_index)))
    if removal_index.size > 0:
        print("Removing following {0} duplicates due to undifferentiated network and station "
              "codes:\n{1}".format(len(removal_index), df.loc[removal_index]))
        df.drop(removal_index, inplace=True)

    return df


def populate_default_station_dates(df):
    """
    Replace all missing station start and end dates with default values.
    Replace all missing channel start and end dates with their corresponding station/end dates.

    :param df: Dataframe in which to fill in missing station start and end dates.
    :type df: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    """
    # Do it for both station dates AND channel dates, as seiscomp3 will treat empty dates as
    # overlapping other time intervals and discard records.
    st_isna_start_mask = df.StationStart.isna()
    st_isna_end_mask = df.StationEnd.isna()
    df.loc[st_isna_start_mask, 'StationStart'] = DEFAULT_START_TIMESTAMP
    df.loc[st_isna_end_mask, 'StationEnd'] = DEFAULT_END_TIMESTAMP
    assert not np.any(df.StationStart.isna())
    assert not np.any(df.StationEnd.isna())

    ch_isna_start_mask = df.ChannelStart.isna()
    ch_isna_end_mask = df.ChannelEnd.isna()
    df.loc[ch_isna_start_mask, 'ChannelStart'] = df.StationStart[ch_isna_start_mask]
    df.loc[ch_isna_end_mask, 'ChannelEnd'] = df.StationEnd[ch_isna_end_mask]
    assert not np.any(df.ChannelStart.isna())
    assert not np.any(df.ChannelEnd.isna())


def merge_overlapping_channel_epochs(df):
    """
    Removed overlapping time intervals for a given network.station.channel, as this
    needs to be unique.

    This function expects the input DataFrame to have a sequential integer index.

    :param df: Dataframe of station records in which to merge overlapping channel dates
    :type df: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    """
    if show_progress:
        pbar = tqdm.tqdm(total=len(df), ascii=True)
    removal_index = []
    for _, data in df.groupby(['NetworkCode', 'StationCode', 'ChannelCode']):
        if show_progress:
            pbar.update(len(data))
        if len(data) <= 1:
            continue
        data.loc[:, 'CumMaxEnd'] = data.ChannelEnd.cummax().shift(1)
        data.loc[:, 'overlap'] = (data.loc[:, 'ChannelStart'] < data.loc[:, 'CumMaxEnd']).astype(int)
        intervals = [(data.index[0], data.index[0])]
        for idx, row in data.iloc[1:].iterrows():
            if row.overlap > 0:
                intervals[-1] = (intervals[-1][0], idx)
            else:
                intervals.append((idx, idx))
        intervals = [i for i in intervals if i[1] > i[0]]
        for first, last in intervals:
            ch_start = data.loc[first:last, 'ChannelStart'].min()
            ch_end = data.loc[first:last, 'ChannelEnd'].max()
            df.loc[first, ['ChannelStart', 'ChannelEnd']] = np.array([ch_start, ch_end])
            removal_index.extend(range(first + 1, last + 1))

    if show_progress:
        pbar.close()

    if removal_index:
        df.drop(removal_index, inplace=True)


def cleanup_database(df):
    """
    Main cleanup function encompassing the sequential data cleanup steps.

    :param df: Dataframe of station records to clean up
    :type df: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    :return: Cleaned up dataframe of station records
    :rtype: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    """
    # Returns cleaned up df.

    print("Removing stations with illegal station code...")
    num_before = len(df)
    remove_illegal_stationNames(df)
    df.reset_index(drop=True, inplace=True)
    if len(df) < num_before:
        print("Removed {0}/{1} stations because their station codes are not "
              "compliant".format(num_before - len(df), num_before))

    print("Cleaning up station duplicates...")
    num_before = len(df)
    neighbor_matrix = compute_neighboring_station_matrix(df)
    df = remove_duplicate_stations(df, neighbor_matrix)
    df.reset_index(drop=True, inplace=True)
    if len(df) < num_before:
        print("Removed {0}/{1} stations flagged as duplicates".format(num_before - len(df), num_before))

    print("Filling in missing station dates with defaults...")
    populate_default_station_dates(df)

    print("Merging overlapping channel epochs...")
    num_before = len(df)
    merge_overlapping_channel_epochs(df)
    df.reset_index(drop=True, inplace=True)
    if len(df) < num_before:
        print("Merged {0}/{1} channel records with overlapping epochs".format(num_before - len(df), num_before))

    return df


def write_portable_inventory(df, fname):
    """
    Write the final database to re-usable file formats.

    :param df: Dataframe containing complete inventory of station records.
    :type df: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    :param fname: Base filename to export to. Output filename will be appended with
        file format extensions.
    :type fname: str
    """
    df.to_csv(fname + ".csv", index=False)
    df.to_hdf(fname + ".h5", mode='w', key='inventory')
    # Human readable, fixed width column format.
    with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 1000):
        inv_str = str(df)
        with open(fname + ".txt", "w") as f:
            f.write(inv_str)


def main(iris_xml_file, stations_folder, output_directory, test_mode=False):
    """
    Main entry point.

    :param iris_xml_file: Name of IRIS xml file to load (generated by script update_iris_inventory.py)
    :type iris_xml_file: str or pathlib.Path
    :param stations_folder: Path to folder containing STN files to process
    :type stations_folder: str or pathlib.Path
    """
    if test_mode:
        use_pickle = False
    else:  # pragma: no cover
        use_pickle = True

    # Read station database from ad-hoc formats
    ehb_data_bmg = read_eng(os.path.join(stations_folder, 'BMG.STN'))
    ehb_data_isc = read_eng(os.path.join(stations_folder, 'ISC.STN'))
    ehb = pd.concat([ehb_data_bmg, ehb_data_isc], sort=False)

    isc1 = read_isc(os.path.join(stations_folder, 'ehb.stn'), use_pickle)
    isc2 = read_isc(os.path.join(stations_folder, 'iscehb.stn'), use_pickle)
    isc = pd.concat([isc1, isc2], sort=False)

    db = pd.concat([ehb, isc], sort=False)

    print("Removing blacklisted networks...")
    db = remove_blacklisted(db)

    # Include date columns in sort so that NaT values sink to the bottom. This means when duplicates are removed,
    # the record with the least NaT values will be favored to be kept.
    db.sort_values(SORT_ORDERING, inplace=True)
    db.reset_index(drop=True, inplace=True)

    # Perform cleanup on each database
    db = cleanup_database(db)

    # Read IRIS station database.
    print("Reading " + iris_xml_file)
    if PY2:
        import io
        with io.open(iris_xml_file, mode='r', encoding='utf-8') as f:
            iris_inv = read_inventory(f)
    else:
        with open(iris_xml_file, mode='r', encoding='utf-8') as f:
            iris_inv = read_inventory(f)

    # Extract nominal sensor and response data from inventory, indexed by channel code.
    nominal_instruments = extract_unique_sensors_responses(iris_inv, requests,
                                                           blacklisted_networks=BLACKLISTED_NETWORKS,
                                                           test_mode=test_mode)

    # Write whole database to FDSN station xml file
    outfile_base = os.path.join(output_directory, "INVENTORY_" + rt_timestamp)
    dataframe_to_fdsn_station_xml(db, nominal_instruments, outfile_base + ".xml")

    # Write whole database in portable csv and hdf5 formats
    write_portable_inventory(db, outfile_base)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--iris", help="Path to IRIS station xml database file.")
    parser.add_argument("-s", "--stations-path", help="Path to folder containing STN files.")
    parser.add_argument("-o", "--output-path", help="Path to folder in which output inventory file(s) should be saved.")
    args = parser.parse_args()
    assert args.iris is not None
    assert args.stations_path is not None
    assert args.output_path is not None
    iris_file = args.iris.strip()
    stations_path = args.stations_path.strip()
    output_path = args.output_path.strip()
    print("Using IRIS source " + iris_file)
    print("Using STN files from " + stations_path)
    print("Output path is " + output_path)

    main(iris_file, stations_path, output_path)
