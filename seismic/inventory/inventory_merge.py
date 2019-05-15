#!/usr/bin/env python
"""Merge a master inventory (e.g. from IRIS) with custom inventories.
"""

# pylint: disable=too-many-locals

import os
import datetime
import pickle as pkl

import click
import numpy as np
import pandas as pd
from obspy.core.inventory import Network, Station
from obspy.geodetics.base import locations2degrees

from seismic.inventory.pdconvert import inventory2Dataframe
from seismic.inventory.inventory_split import split_inventory_by_network
from seismic.inventory.inventory_util import load_station_xml, NOMINAL_EARTH_RADIUS_KM, SORT_ORDERING

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

# Timestamp to be added to output file names, so that each run generates unique log files.
rt_timestamp = datetime.datetime.now().strftime("%Y%m%dT%H%M%S")

# Distance within which IRIS stations will be considered at same notional location as other stations
# that match on other metadata axes.
DIST_TOLERANCE_KM = 5.0  # pylint: disable=invalid-name


def prune_iris_duplicates(db_other, db_iris):
    """Prune the records out of db_other which duplicate records that exist in db_iris.

    Filtering concepts:
    Station codes and channel codes are reliable, but network codes and station/channel dates are not.
    Station lat/lon locations are approximately reliable. Therefore the merge method here matches records
    with matching station code, channel code and approximate location. For each match, we examine station dates,
    and discard any records of db_other that overlap with db_iris.

    db_iris is not changed by this function.

    :param db_other: DataFrame containing metadata records of custom inventory
    :type db_other: pandas.DataFrame
    :param db_iris: DataFrame containing metadata records of IRIS inventory
    :type db_iris: pandas.DataFrame
    :return: DataFrame containing records from db_other which are unique and do not appear in IRIS.
    :rtype: pandas.DataFrame
    """

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


@click.command()
@click.option('--iris-inv', type=click.Path(exists=True, dir_okay=False, readable=True), required=True,
              help='IRIS inventory file in FDSN stationxml format. These records take precedence over custom files.')
@click.option('--custom-inv', type=click.Path(exists=True, dir_okay=False, readable=True), required=True,
              help='Custom inventory file in FDSN stationxml format to merge with the IRIS inventory.')
@click.option('--output-file', type=click.Path(exists=False, dir_okay=False), required=True,
              help='Name of file to output final result to.')
def main(iris_inv, custom_inv, output_file, split_output_folder=None):
    """Merge an IRIS inventory with a custom inventory, filtering out any custom inventory records
       that duplicate IRIS records.

    :param iris_inv: Station XML file from which to load IRIS inventory
    :type iris_inv: str or Path to file
    :param custom_inv: Station XML file with custom records to merge with IRIS inventory
    :type custom_inv: str or Path to file
    :param output_file: File name of output file which will contain merged IRIS and custom inventory.
    :type output_file: str or Path to file
    """

    # Load IRIS inventory
    print("Loading {}...".format(iris_inv))
    if os.path.exists(os.path.splitext(iris_inv)[0] + ".pkl"):
        with open(os.path.splitext(iris_inv)[0] + ".pkl", 'rb') as f:
            inv_iris = pkl.load(f)
            db_iris = pkl.load(f)
    else:
        inv_iris = load_station_xml(iris_inv)
        # and convert to Pandas dataframe
        db_iris = inventory2Dataframe(inv_iris)
        with open(os.path.splitext(iris_inv)[0] + ".pkl", 'wb') as f:
            pkl.dump(inv_iris, f, pkl.HIGHEST_PROTOCOL)
            pkl.dump(db_iris, f, pkl.HIGHEST_PROTOCOL)

    # Load custom inventory
    print("Loading {}...".format(custom_inv))
    if os.path.exists(os.path.splitext(custom_inv)[0] + ".pkl"):
        with open(os.path.splitext(custom_inv)[0] + ".pkl", 'rb') as f:
            inv_other = pkl.load(f)
            db_other = pkl.load(f)
    else:
        inv_other = load_station_xml(custom_inv)
        # and convert to Pandas dataframe
        db_other = inventory2Dataframe(inv_other)
        with open(os.path.splitext(custom_inv)[0] + ".pkl", 'wb') as f:
            pkl.dump(inv_other, f, pkl.HIGHEST_PROTOCOL)
            pkl.dump(db_other, f, pkl.HIGHEST_PROTOCOL)

    print("Merging {} IRIS records with {} custom records...".format(len(db_iris), len(db_other)))
    num_before = len(db_other)
    # removeIrisDuplicates(db_other, db_iris)
    db_other = prune_iris_duplicates(db_other, db_iris)
    db_other.sort_values(SORT_ORDERING, inplace=True)
    db_other.reset_index(drop=True, inplace=True)
    if len(db_other) < num_before:
        print("Removed {0}/{1} stations because they exist in IRIS".format(num_before - len(db_other), num_before))
    print("{} custom records remaining".format(len(db_other)))

    # Merge inv_other into inv_iris, only keeping the records of inv_other that are present in db_other
    inv_merged = inv_iris  # Note: this aliases inv_iris, it does not make a copy.
    print("Filtering {} records and merging with IRIS...".format(custom_inv))
    num_added = 0
    if show_progress:
        num_entries = sum(len(station.channels) for network in inv_other.networks for station in network.stations)
        pbar = tqdm.tqdm(total=num_entries, ascii=True)
        pbar.set_description("Matched {}/{}".format(num_added, len(db_other)))
    for network in inv_other.networks:
        # Duplicate network data, but keep stations empty
        net = Network(network.code, stations=[], description=network.description, comments=network.comments,
                      start_date=network.start_date, end_date=network.end_date)
        add_network = False
        for station in network.stations:
            if show_progress:
                pbar.update(len(station.channels))
            # Duplicate station data, but keep channels empty
            sta = Station(station.code, station.latitude, station.longitude, station.elevation, channels=[],
                          site=station.site, creation_date=station.creation_date,
                          termination_date=station.termination_date, description=station.description,
                          comments=station.comments, start_date=station.start_date, end_date=station.end_date)
            for channel in station.channels:
                lat = channel.latitude if channel.latitude else station.latitude
                lon = channel.longitude if channel.longitude else station.longitude
                ele = channel.elevation if channel.elevation else station.elevation
                sta_start = np.datetime64(station.start_date)
                cha_start = np.datetime64(channel.start_date)
                mask = ((db_other['NetworkCode'] == network.code) &
                        (db_other['StationCode'] == station.code) &
                        (db_other['ChannelCode'] == channel.code) &
                        (db_other['Latitude'] == lat) &
                        (db_other['Longitude'] == lon) &
                        (db_other['Elevation'] == ele) &
                        (db_other['StationStart'] == sta_start) &
                        (db_other['ChannelStart'] == cha_start)
                       )
                if np.any(mask):
                    db_match = db_other[mask]
                    assert len(db_match) == 1, 'Found multiple matches, expected only one for {}'.format(db_match)
                    add_network = True
                    num_added += 1
                    if show_progress:
                        pbar.set_description("Matched {}/{}".format(num_added, len(db_other)))
                    sta.channels.append(channel)
            # end for
            if add_network:
                net.stations.append(sta)
        # end for
        if add_network:
            inv_merged += net
    # end for
    if show_progress:
        pbar.close()

    print("Added {} custom records to IRIS inventory".format(num_added))

    # Write merged inventory text file in FDSN stationxml inventory format.
    print("Writing merged inventory to {}".format(output_file))
    inv_merged.write(output_file, format="stationxml")
    inv_merged.write(os.path.splitext(output_file)[0] + ".txt", format="stationtxt")

    if split_output_folder is not None:
        print("Writing separate network files to folder {}".format(split_output_folder))
        split_inventory_by_network(inv_merged, split_output_folder)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
