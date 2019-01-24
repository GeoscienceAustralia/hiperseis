#!/usr/bin/env python
"""Creates database of stations from .STN files which are not in IRIS database,
   curates the data using heuristic rules, and exports new stations to station.xml
"""
from __future__ import division

import os
import sys
import argparse

import numpy as np
import scipy as sp
import datetime
import pandas as pd
from obspy import read_inventory
from obspy.geodetics.base import locations2degrees
from pdconvert import pd2Network
from plotting import saveNetworkLocalPlots
from table_format import TABLE_SCHEMA, TABLE_COLUMNS, PANDAS_MAX_TIMESTAMP, DEFAULT_START_TIMESTAMP, DEFAULT_END_TIMESTAMP

if sys.version_info[0] < 3:
    import cStringIO as sio
    import pathlib2 as pathlib
    import cPickle as pkl
else:
    import io as sio
    import pathlib
    import pickle as pkl

try:
    import tqdm
    show_progress = True
except:
    show_progress = False
    print("Run 'pip install tqdm' to see progress bar.")

# Script requires numpy >= 1.15.4. The source of the error is not yet identified, but has been
# demonstrated on multiple platforms with lower versions of numpy.
vparts = np.version.version.split('.', 2)
(major, minor, maint) = [int(x) for x in vparts]
if major < 1 or (major == 1 and minor < 14):
    print("Not supported error: Requires numpy >= 1.14.2, found numpy {0}".format(".".join(vparts)))
    sys.exit(1)
else:
    print("Using numpy {0}".format(".".join(vparts)))

USE_PICKLE = True

TEST_MODE = False

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', -1)
pd.set_option('display.width', 240)

NOMINAL_EARTH_RADIUS_KM = 6378.1370
DIST_TOLERANCE_KM = 2.0
DIST_TOLERANCE_RAD = DIST_TOLERANCE_KM / NOMINAL_EARTH_RADIUS_KM
COSINE_DIST_TOLERANCE = np.cos(DIST_TOLERANCE_RAD)

# See https://pandas.pydata.org/pandas-docs/stable/timeseries.html#timestamp-limitations


rt_timestamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")


def read_eng(fname):
    """Read Engdahl file having the following format of fixed width formatted columns:

     AAI   Ambon             BMG, Indonesia, IA-Ne              -3.6870  128.1945      0.0   2005001  2286324  I
     AAII                                                       -3.6871  128.1940      0.0   2005001  2286324  I
     AAK   Ala Archa         Kyrgyzstan                         42.6390   74.4940      0.0   2005001  2286324  I
     ABJI                                                       -7.7957  114.2342      0.0   2005001  2286324  I
     APSI                                                       -0.9108  121.6487      0.0   2005001  2286324  I
     AS01  Alice Springs Arra                                  -23.6647  133.9508      0.0   2005001  2286324  I
    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
              10        20        30        40        50        60        70        80

    Each Station Code (first column) might NOT be unique, and network codes are missing here, so duplicates will
    simply be noted and handled later.

    Returns Pandas Dataframe containing the loaded data in column order of TABLE_COLUMNS.
    """
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


# @profile
def read_isc(fname):
    """Read ISC station inventory supplied by ISC that inherited Engdahl work, having such format:

    109C     32.8892 -117.1100     0.0 2006-06-01 04:11:18 2008-01-04 01:26:30
    109C     32.8882 -117.1050   150.0 2008-01-04 01:26:30
                 FDSN 109C   TA -- BHZ 2004-05-04 23:00:00 2005-03-03 23:59:59
                 FDSN 109C   TA -- LHZ 2004-05-04 23:00:00 2005-03-03 23:59:59
                 FDSN 109C   TA -- BHZ 2005-04-11 00:00:00 2006-01-25 22:31:10
                 FDSN 109C   TA -- LHZ 2005-04-11 00:00:00 2006-01-25 22:31:10
    0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
              10        20        30        40        50        60        70        80

    The lines starting with a station code are HEADER rows, and provide the station coordinates.
    The idented lines starting with FDSN provide distinct station, network and channel data for the given station location.

    Returns Pandas Dataframe containing the loaded data in column order of TABLE_COLUMNS.
    """
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
        """Convenience function to report on number of unique network and station codes in the dataframe."""
        num_unique_networks = len(df['NetworkCode'].unique())
        num_unique_stations = len(df['StationCode'].unique())
        print("{0}: {1} unique network codes, {2} unique station codes".format(fname, num_unique_networks, num_unique_stations))

    if USE_PICKLE:
        pkl_name = fname + ".pkl"
        if os.path.exists(pkl_name):
            print("Reading cached " + fname)
            with open(pkl_name, 'rb') as f:
                df_all = pkl.load(f)
                reportStationCount(df_all)
                return df_all

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
                hdr['NetworkCode'] = 'IR'
                hdr['ChannelStart'] = pd.NaT
                hdr['ChannelEnd'] = pd.NaT
                hdr['ChannelCode'] = 'BHZ'
                # Standardize column ordering
                hdr = hdr[list(TABLE_COLUMNS)]
                if channels:
                    # If channel data is also present, store it too.
                    ch_all = pd.concat(channels, sort=False)
                    ch_all.drop('FDSN', axis=1, inplace=True)
                    # Set the station date range to at least encompass the channels it contains.
                    st_min = min(hdr['StationStart'].min(), ch_all['ChannelStart'].min())
                    st_max = max(hdr['StationEnd'].max(), ch_all['ChannelEnd'].max())
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
            line_input = sio.StringIO(line)
            hdr = pd.read_fwf(line_input, colspecs=header_colspec, names=header_cols, nrows=1, dtype=TABLE_SCHEMA, parse_dates=[4, 5])
            line = f.readline()
            if show_progress:
                pbar.update(len(line))
    if show_progress:
        pbar.close()
    print("Concatenating records...")
    df_all = pd.concat(df_list, sort=False)
    if USE_PICKLE:
        with open(fname + ".pkl", "wb") as f:
            pkl.dump(df_all, f, pkl.HIGHEST_PROTOCOL)
    reportStationCount(df_all)
    return df_all


# @profile
def removeIllegalStationNames(df):
    """Remove records for station names that do not conform to expected naming convention.

       Such names can cause problems in downstream processing, in particular names with asterisk.
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


# @profile
def latLong2CosineDistance(latlong_deg_set1, latlong_deg_set2):
    """Compute the approximate cosine distance between each station of 2 sets.

       Each set is specified as a numpy column vector of [latitude, longitude] positions in degrees.

       This function performs an outer product and will produce matrix of size N0 x N1, where
       N0 is the number of rows in latlong_deg_set1 and N1 is the number of rows in latlong_deg_set2.

       Returns np.ndarray containing cosines of angles between each pair of stations from
       the input arguments.
    """
    # If input is 1D, convert to 2D for consistency of matrix orientations.
    if len(latlong_deg_set1.shape) == 1:
        latlong_deg_set1 = np.reshape(latlong_deg_set1, (1, -1))
    if len(latlong_deg_set2.shape) == 1:
        latlong_deg_set2 = np.reshape(latlong_deg_set2, (1, -1))

    set1_latlong_rad = np.deg2rad(latlong_deg_set1)
    set2_latlong_rad = np.deg2rad(latlong_deg_set2)

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


# @profile
def removeIrisDuplicates(df, iris_inv):
    """Remove stations which duplicate records in IRIS database.

       Remove stations from df which are present in the global IRIS inventory, in cases
       where the station codes match and the distance from the IRIS station is within threshold.

       Mutates df in-place.
    """
    if show_progress:
        pbar = tqdm.tqdm(total=len(df), ascii=True)
    removal_index = []
    with open("LOG_IRIS_DUPES_" + rt_timestamp + ".txt", 'w') as log:
        for (netcode, statcode), data in df.groupby(['NetworkCode', 'StationCode']):
            iris_query = iris_inv.select(network=netcode, station=statcode, channel="*HZ")
            if len(iris_query) <= 0:
                # No IRIS record matching this station
                if show_progress:
                    pbar.update(len(data))
                continue
            # Pull out matching stations. Since some station codes have asterisk, which is interpreted as a wildcard
            # by the obspy query, we need to filter against matching exact statcode.
            matching_stations = [s for n in iris_query.networks for s in n.stations if s.code == statcode and n.code == netcode]
            iris_station0 = matching_stations[0]
            # Check that the set of stations from IRIS are themselves within the distance tolerance of one another.
            iris_stations_dist = [np.deg2rad(locations2degrees(iris_station0.latitude, iris_station0.longitude, s.latitude, s.longitude)) * NOMINAL_EARTH_RADIUS_KM
                                  for s in matching_stations]
            iris_stations_dist = np.array(iris_stations_dist)
            within_tolerance_mask = (iris_stations_dist < DIST_TOLERANCE_KM)
            if not np.all(within_tolerance_mask):
                log.write("WARNING: Not all IRIS stations localized within distance tolerance for station code {0}.{1}. "
                          "Distances(km) = {2}\n".format(netcode, statcode, iris_stations_dist[~within_tolerance_mask]))
            # Compute cosine distances between this group's set of stations and the IRIS station locagtion.
            ref_latlong = np.array([iris_station0.latitude, iris_station0.longitude])
            stations_latlong = data[["Latitude", "Longitude"]].values
            distfunc = lambda r: np.deg2rad(locations2degrees(ref_latlong[0], ref_latlong[1], r[0], r[1])) * NOMINAL_EARTH_RADIUS_KM  # noqa
            surface_dist = np.apply_along_axis(distfunc, 1, stations_latlong)
            # assert isinstance(surface_dist, np.ndarray)
            if not surface_dist.shape:
                surface_dist = np.reshape(surface_dist, (1,))
            mask = (surface_dist < DIST_TOLERANCE_KM)
            if np.isscalar(mask):
                mask = np.array([mask], ndmin=1)
            duplicate_index = np.array(data.index.tolist())[mask]
            if len(duplicate_index) < len(data):
                kept_station_distances = surface_dist[(~mask)]
                log.write("WARNING: Some ISC stations outside distance tolerance of IRIS location for station {0}.{1}, not dropping. "
                          "(Possible issues with station date ranges?) Distances(km) = {2}\n".format(netcode, statcode, kept_station_distances))
            removal_index.extend(duplicate_index.tolist())
            if show_progress:
                pbar.update(len(data))
    if show_progress:
        pbar.close()

    if removal_index:
        df.drop(removal_index, inplace=True)


# @profile
def computeNeighboringStationMatrix(df):
    """Compute sparse matrix representing index of neighboring stations.

       Ordering of matrix corresponds to ordering of Dataframe df, which is expected to be sequential
       integer indexed. For a given station index i, then the non-zero off-diagonal entries in row i
       of the returned matrix indicate the indices of adjacent, nearby stations.

       Returns sparse binary matrix having non-zero values at indices of neighboring stations.
    """
    self_latlong = df[["Latitude", "Longitude"]].values
    # In order to keep calculation tractable for potentially large matrices without resort to
    # on-disk memmapped arrays, we split the second operand into parts and compute parts of the
    # result, then recombine them to get the final (sparse) result.
    sparse_cos_dist = []
    num_splits = max(self_latlong.shape[0] // 2000, 1)
    for m in np.array_split(self_latlong, num_splits):
        partial_result = latLong2CosineDistance(self_latlong, m)
        partial_result = sp.sparse.csr_matrix(partial_result >= COSINE_DIST_TOLERANCE)
        sparse_cos_dist.append(partial_result)
    cos_dist = sp.sparse.hstack(sparse_cos_dist, "csr")
    return cos_dist


# @profile
def removeDuplicateStations(df, neighbor_matrix):
    """Remove stations which are identified as duplicates.

       Remove duplicated stations in df based on station code and locality of lat/long coordinates.
       Remove duplicated station based on codes and channel data, IRRESPECTIVE of locality.

       Mutates df in-place.
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
            attrs_match = np.array([((k[1] == key) | (k[1].isna() & key.isna())) for k in df.loc[neighbors, matching_criteria].iterrows()])
            duplicate_mask = np.all(attrs_match, axis=1)
            if np.any(duplicate_mask):
                duplicate_index = neighbors[duplicate_mask]
                log.write("WARNING: Duplicates of\n{0}\nare being removed:\n{1}\n----\n".format(df.loc[[i]], df.loc[duplicate_index]))
                removal_rows.update(duplicate_index)
        if show_progress:
            pbar.close()
    removal_rows = np.array(sorted(list(removal_rows)))
    if removal_rows.size > 0:
        print("Removing following {0} duplicates due to identical network, station and channel data:\n{1}".format(len(removal_rows), df.loc[removal_rows]))
        df.drop(removal_rows, inplace=True)

    # Secondly, remove stations with same network and station code, but which are further away than the threshold distance
    # and have no distinguishing channel data. We deliberately exclude station start and end dates from consideration here,
    # as these are extended during file read to cover range of contained channels, and therefore might not match in code dupes.
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
                # Note that NA fields will compare False even if both are NA, which is what we want here since we don't want to treat
                # records with same codes as duplicates if the matching_criteria are NA, as this removes records that are obviously
                # not duplicates.
                duplicate_mask = (data[matching_criteria] == key)
                index_mask = np.all(duplicate_mask, axis=1) & (data.index > row_index)
                duplicate_index = data.index[index_mask]
                if not duplicate_index.empty:
                    log.write("WARNING: Apparent duplicates of\n{0}\nare being removed:\n{1}\n----\n".format(data.loc[[row_index]], data.loc[duplicate_index]))
                    removal_index.update(duplicate_index.tolist())
        if show_progress:
            pbar.close()
    removal_index = np.array(sorted(list(removal_index)))
    if removal_index.size > 0:
        print("Removing following {0} duplicates due to undifferentiated network and station codes:\n{1}".format(len(removal_index), df.loc[removal_index]))
        df.drop(removal_index, inplace=True)

    return df


def populateDefaultStationDates(df):
    """Replace all missing station start and end dates with default values.
    """
    df.StationStart[df.StationStart.isna()] = DEFAULT_START_TIMESTAMP
    df.StationEnd[df.StationEnd.isna()] = DEFAULT_END_TIMESTAMP


# @profile
def cleanupDatabase(df, iris_inv):
    """Clean up the dataframe df.

       Returns cleaned up df.
    """

    print("Removing stations with illegal station code...")
    num_before = len(df)
    removeIllegalStationNames(df)
    df.reset_index(drop=True, inplace=True)
    if len(df) < num_before:
        print("Removed {0}/{1} stations because their station codes are not compliant".format(num_before - len(df), num_before))

    print("Removing stations which replicate IRIS...")
    num_before = len(df)
    removeIrisDuplicates(df, iris_inv)
    df.reset_index(drop=True, inplace=True)
    if len(df) < num_before:
        print("Removed {0}/{1} stations because they exist in IRIS".format(num_before - len(df), num_before))

    print("Cleaning up station duplicates...")
    num_before = len(df)
    neighbor_matrix = computeNeighboringStationMatrix(df)
    df = removeDuplicateStations(df, neighbor_matrix)
    df.reset_index(drop=True, inplace=True)
    if len(df) < num_before:
        print("Removed {0}/{1} stations flagged as duplicates".format(num_before - len(df), num_before))

    print("Filling in missing station dates with defaults...")
    populateDefaultStationDates(df)

    return df


def exportStationXml(df, output_folder, filename_base):
    """Export the dataset in df to Station XML format.

       Given a dataframe containing network and station codes grouped under networks, for each
       network create and obspy inventory object and export to stationxml file. Write an overall
       list of stations based on global inventory.
    """
    from obspy.core.inventory import Inventory, Network, Station, Channel, Site

    pathlib.Path(output_folder).mkdir(exist_ok=True)
    print("Exporting stations to folder {0}".format(output_folder))

    if show_progress:
        pbar = tqdm.tqdm(total=len(df), ascii=True)
        progressor = pbar.update
    else:
        progressor = None

    global_inventory = Inventory(networks=[], source='EHB')
    for netcode, data in df.groupby('NetworkCode'):
        net = pd2Network(netcode, data, progressor)
        net_inv = Inventory(networks=[net], source=global_inventory.source)
        global_inventory.networks.append(net)
        fname = "{0}{1}.xml".format(filename_base, netcode)
        try:
            net_inv.write(os.path.join(output_folder, fname), format="stationxml", validate=True)
        except Exception as e:
            print(e)
            print("FAILED writing file {0} for network {1}, continuing".format(fname, netcode))
            continue
    if show_progress:
        pbar.close()

    # Write global inventory text file in FDSN stationtxt inventory format.
    global_inventory.write("station.txt", format="stationtxt")


def writeFinalInventory(df, fname):
    """Write the final database to re-usable file formats."""
    df.to_csv(fname + ".csv", index=False)
    df.to_hdf(fname + ".h5", mode='w', key='inventory')


def exportNetworkPlots(df, plot_folder):
    if show_progress:
        pbar = tqdm.tqdm(total=len(df), ascii=True)
    try:
        saveNetworkLocalPlots(df, plot_folder, pbar.update)
        if show_progress:
            pbar.close()
    except:
        if show_progress:
            pbar.close()
        raise


def main(iris_xml_file):
    # Read station database from ad-hoc formats
    if TEST_MODE:
        ehb_data_bmg = read_eng(os.path.join('test', 'BMG_test.STN'))
        ehb_data_isc = read_eng(os.path.join('test', 'ISC_test.STN'))
    else:
        ehb_data_bmg = read_eng('BMG.STN')
        ehb_data_isc = read_eng('ISC.STN')
    ehb = pd.concat([ehb_data_bmg, ehb_data_isc], sort=False)

    if TEST_MODE:
        isc1 = read_isc(os.path.join('test', 'ehb_test.stn'))
        isc2 = read_isc(os.path.join('test', 'iscehb_test.stn'))
    else:
        isc1 = read_isc('ehb.stn')
        isc2 = read_isc('iscehb.stn')
    isc = pd.concat([isc1, isc2], sort=False)

    db = pd.concat([ehb, isc], sort=False)
    # Include date columns in sort so that NaT values sink to the bottom. This means when duplicates are removed,
    # the record with the least NaT values will be favored to be kept.
    db.sort_values(['NetworkCode', 'StationCode', 'StationStart', 'StationEnd', 'ChannelCode', 'ChannelStart', 'ChannelEnd'], inplace=True)
    db.reset_index(drop=True, inplace=True)

    # Read IRIS station database.
    print("Reading " + iris_xml_file)
    if False:  # TODO: only keep this code if proven that pickled station xml file loads faster...
        IRIS_all_pkl_file = iris_xml_file + ".pkl"
        if os.path.exists(IRIS_all_pkl_file):
            with open(IRIS_all_pkl_file, 'rb') as f:
                iris_inv = pkl.load(f)
        else:
            with open(iris_xml_file, 'r') as f:
                iris_inv = read_inventory(f)
            with open(IRIS_all_pkl_file, 'wb') as f:
                pkl.dump(iris_inv, f, pkl.HIGHEST_PROTOCOL)
    else:
        with open(iris_xml_file, 'r') as f:
            iris_inv = read_inventory(f)

    # Perform cleanup on each database
    db = cleanupDatabase(db, iris_inv)

    if TEST_MODE:
        output_folder = "output_test"
    else:
        output_folder = "output"

    exportStationXml(db, output_folder, "network_")

    writeFinalInventory(db, "INVENTORY_" + rt_timestamp)

    plot_folder = "plots"
    print("Exporting network plots to folder {0}".format(plot_folder))
    exportNetworkPlots(db, plot_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--iris", help="Path to IRIS station xml database file to use to exclude station codes from STN sources.")
    args = parser.parse_args()
    if args.iris is None:
        parser.print_help()
        sys.exit(0)
    else:
        iris_xml_file = args.iris.strip()
    print("Using IRIS source " + iris_xml_file)

    # import cProfile as prof
    # statsfile = 'perfstats.stat'
    # # prof.run('main(iris_xml_file)', statsfile)
    # prof.run('read_isc(os.path.join(\'test\', \'iscehb_test.stn\'))', statsfile)
    # import pstats
    # p = pstats.Stats(statsfile)
    # p.sort_stats('tottime').print_stats(50)
    main(iris_xml_file)
