#!/usr/bin/env python
"""Merge a master inventory (e.g. from IRIS) with custom inventories.
"""

import os
import sys
import datetime
import pickle as pkl

import click
import numpy as np
import pandas as pd
import obspy
from obspy import read_inventory

from seismic.inventory.pdconvert import inventory2Dataframe
from obspy.geodetics.base import locations2degrees

from seismic.inventory.inventory_util import NOMINAL_EARTH_RADIUS_KM, SORT_ORDERING

try:
    import tqdm
    show_progress = True
except ImportError:
    show_progress = False
    print("Run 'pip install tqdm' to see progress bar.")

# Pandas table display options to reduce aggressiveness of truncation. Due to size of data sometimes we
# need to see more details in the table.
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_colwidth', -1)
pd.set_option('display.width', 240)

PY2 = sys.version_info[0] < 3

# Timestamp to be added to output file names, so that each run generates unique log files.
rt_timestamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")

# Distance within which IRIS stations will be considered at same notional location as other stations
# that match on other metadata axes.
DIST_TOLERANCE_KM = 5.0

# def removeIrisDuplicates(db_other, db_iris):
#     """
#     Remove stations which duplicate records in IRIS database.
#     """
#     if show_progress:
#         pbar = tqdm.tqdm(total=len(db_other), ascii=True)
#     removal_index = []
#     with open("LOG_IRIS_DUPES_" + rt_timestamp + ".txt", 'w') as log:
#         for (netcode, statcode), data in db_other.groupby(['NetworkCode', 'StationCode']):
#             iris_query = iris_inv.select(network=netcode, station=statcode, channel="*HZ")
#             if len(iris_query) <= 0:
#                 # No IRIS record matching this station
#                 if show_progress:
#                     pbar.update(len(data))
#                 continue
#             # Pull out matching stations. Since some station codes have asterisk, which is interpreted as a wildcard
#             # by the obspy query, we need to filter against matching exact statcode.
#             matching_stations = [s for n in iris_query.networks for s in n.stations if s.code == statcode and n.code == netcode]
#             iris_station0 = matching_stations[0]
#             # Check that the set of stations from IRIS are themselves within the distance tolerance of one another.
#             iris_stations_dist = [np.deg2rad(locations2degrees(iris_station0.latitude, iris_station0.longitude, s.latitude, s.longitude)) * NOMINAL_EARTH_RADIUS_KM
#                                   for s in matching_stations]
#             iris_stations_dist = np.array(iris_stations_dist)
#             within_tolerance_mask = (iris_stations_dist < DIST_TOLERANCE_KM)
#             if not np.all(within_tolerance_mask):
#                 log.write("WARNING: Not all IRIS stations localized within distance tolerance for station code {0}.{1}. "
#                           "Distances(km) = {2}\n".format(netcode, statcode, iris_stations_dist[~within_tolerance_mask]))
#             # Compute cosine distances between this group's set of stations and the IRIS station location.
#             ref_latlong = np.array([iris_station0.latitude, iris_station0.longitude])
#             stations_latlong = data[["Latitude", "Longitude"]].values
#             distfunc = lambda r: np.deg2rad(locations2degrees(ref_latlong[0], ref_latlong[1], r[0], r[1])) * NOMINAL_EARTH_RADIUS_KM  # noqa
#             surface_dist = np.apply_along_axis(distfunc, 1, stations_latlong)
#             # assert isinstance(surface_dist, np.ndarray)
#             if not surface_dist.shape:
#                 surface_dist = np.reshape(surface_dist, (1,))
#             mask = (surface_dist < DIST_TOLERANCE_KM)
#             if np.isscalar(mask):
#                 mask = np.array([mask], ndmin=1)
#             duplicate_index = np.array(data.index.tolist())[mask]
#             if len(duplicate_index) < len(data):
#                 # Disable false positive from pylint on following line
#                 kept_station_distances = surface_dist[(~mask)]  # pylint: disable=invalid-unary-operand-type
#                 log.write("WARNING: Some ISC stations outside distance tolerance of IRIS location for station {0}.{1}, not dropping. "
#                           "(Possible issues with station date ranges?) Distances(km) = {2}\n".format(netcode, statcode, kept_station_distances))
#             removal_index.extend(duplicate_index.tolist())
#             if show_progress:
#                 pbar.update(len(data))
#     if show_progress:
#         pbar.close()

#     if removal_index:
#         df.drop(removal_index, inplace=True)


def prune_iris_duplicates(db_other, db_iris):
    # Station codes and channel codes are reliable, but network codes and station/channel dates are not.
    # Station lat/lon locations are approximately reliable.
    # Therefore the merge method here matches records with matching station code, channel code and approximate
    # location. For each match, we examine station dates, and discard any records of db_other that overlap
    # with db_iris.
    #
    # db_iris is not changed by this function.

    def earth_distance(row):
        lat0 = row['Latitude_left']
        long0 = row['Longitude_left']
        lat1 = row['Latitude_right']
        long1 = row['Longitude_right']
        return np.deg2rad(locations2degrees(lat0, long0, lat1, long1)) * NOMINAL_EARTH_RADIUS_KM

    db_other['Index'] = db_other.index
    db_iris['Index'] = db_iris.index
    df_m = db_other.reset_index().merge(db_iris, on=['StationCode', 'ChannelCode'], suffixes=('_left', '_right')
                                        ).set_index('index')
    dist = df_m.apply(earth_distance, axis=1)
    df_m['Distance'] = dist
    # Reorder columns for easier visual comparison
    df_m = df_m[['NetworkCode_left', 'NetworkCode_right', 'StationCode', 'ChannelCode',
                 'Latitude_left', 'Latitude_right', 'Longitude_left', 'Longitude_right', 'Distance',
                 'Elevation_left', 'Elevation_right',
                 'StationStart_left', 'StationStart_right', 'StationEnd_left', 'StationEnd_right',
                 'ChannelStart_left', 'ChannelStart_right', 'ChannelEnd_left', 'ChannelEnd_right',
                 'Index_left', 'Index_right']]

    df_m = df_m[(df_m['Distance'] <= DIST_TOLERANCE_KM)]

    # # Make sure that for those records that we keep from db_other, the network codes match those from IRIS
    # # (irrespective of dates).
    # db_other[(db_other.index == df_m['Index_left']), 'NetworkCode'] = df_m['NetworkCode_right']

    # Identify which index values in df_m represent overlapping channel dates of IRIS and other data records.
    # The dates that overlap are the records that we DON'T want to keep.
    overlapping_dates = ((df_m['ChannelStart_left'] < df_m['ChannelEnd_right']) &
                         (df_m['ChannelEnd_left'] > df_m['ChannelStart_right']))
    index_to_keep = df_m.loc[~overlapping_dates, 'Index_left'].unique()

    db_pruned = db_other.loc[index_to_keep]

    return db_pruned

def load_FDSN_station_xml(inventory_file):  # pylint: disable=invalid-name
    """Load a FDSN stationxml file and convert to Pandas dataframe

    :param inventory_file: [description]
    :type inventory_file: str or path
    """
    if PY2:
        import io
        with io.open(inventory_file, mode='r', encoding='utf-8') as f:
            obspy_inv = read_inventory(f)
    else:
        with open(inventory_file, mode='r', encoding='utf-8') as f:
            obspy_inv = read_inventory(f)

    db = inventory2Dataframe(obspy_inv)
    return db


@click.command()
@click.option('--iris-inv', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True,
              help='IRIS inventory file in FDSN stationxml format. These records take precedence over custom files.')
@click.option('--extra-inv', type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True), required=True,
              help='Extra inventory file in FDSN stationxml format to merge with the IRIS inventory.')
@click.option('--output-file', type=click.Path(exists=False, dir_okay=False), required=True,
              help='Name of file to output final result to.')
def main(iris_inv, extra_inv, output_file):

    # Load IRIS inventory
    if os.path.exists(os.path.splitext(iris_inv)[0] + ".pkl"):
        with open(os.path.splitext(iris_inv)[0] + ".pkl", 'rb') as f:
            db_iris = pkl.load(f)
    else:
        db_iris = load_FDSN_station_xml(iris_inv)
        with open(os.path.splitext(iris_inv)[0] + ".pkl", 'wb') as f:
            pkl.dump(db_iris, f, pkl.HIGHEST_PROTOCOL)

    # Load custom inventory
    if os.path.exists(os.path.splitext(extra_inv)[0] + ".pkl"):
        with open(os.path.splitext(extra_inv)[0] + ".pkl", 'rb') as f:
            db_other = pkl.load(f)
    else:
        db_other = load_FDSN_station_xml(extra_inv)
        with open(os.path.splitext(extra_inv)[0] + ".pkl", 'wb') as f:
            pkl.dump(db_other, f, pkl.HIGHEST_PROTOCOL)

    print("Merging {} IRIS records with {} custom records".format(len(db_iris), len(db_other)))
    num_before = len(db_other)
    # removeIrisDuplicates(db_other, db_iris)
    db_other = prune_iris_duplicates(db_other, db_iris)
    if len(db_other) < num_before:
        print("Removed {0}/{1} stations because they exist in IRIS".format(num_before - len(db_other), num_before))

    # Merge db_iris and db_other, then sort by network, station, channel, etc...
    db_merged = pd.concat([db_iris, db_other], sort=False)
    db_merged.sort_values(SORT_ORDERING, inplace=True)
    db_merged.reset_index(drop=True, inplace=True)

    # Save merged to new index


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
