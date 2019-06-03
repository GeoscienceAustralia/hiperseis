#!/usr/bin/env python
"""Functions to read picks from custom csv file into Pandas DataFrame format and related
   data filtering utilities.
"""

import datetime

from collections import OrderedDict

import numpy as np
import pandas as pd

import obspy

# pylint: disable=invalid-name

PICKS_TABLE_SCHEMA = OrderedDict((
    ('#eventID', object),
    ('originTimestamp', np.float64),
    ('mag', np.float64),
    ('originLon', np.float64),
    ('originLat', np.float64),
    ('originDepthKm', np.float64),
    ('net', object),
    ('sta', object),
    ('cha', object),
    ('pickTimestamp', np.float64),
    ('phase', object),
    ('stationLon', np.float64),
    ('stationLat', np.float64),
    ('az', np.float64),
    ('baz', np.float64),
    ('distance', np.float64),
    ('ttResidual', np.float64),
    ('snr', np.float64),
    ('qualityMeasureCWT', np.float64),
    ('domFreq', np.float64),
    ('qualityMeasureSlope', np.float64),
    ('bandIndex', np.int64),
    ('nSigma', np.int64)))

PICKS_TABLE_COLUMNS = PICKS_TABLE_SCHEMA.keys()


def read_picks_ensemble(csv_file):
    """Read picks from CSV file using PICKS_TABLE_SCHEMA and return as a Pandas DataFrame.

    :param csv_file: Input file containing picks ensemble
    :type csv_file: str or path
    """
    df_raw_picks = pd.read_csv(csv_file, ' ', header=0, dtype=PICKS_TABLE_SCHEMA)
    return df_raw_picks


def get_network_stations(df, netcode):
    """
    Get the unique station codes belonging to a given network from a Pandas DataFrame

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param netcode: Network code for which station codes will be returned
    :type netcode: str
    :return: Sorted list of station code strings extracted from df
    :rtype: list(str)
    """
    return sorted(df[df['net'] == netcode]['sta'].unique().tolist())


def get_network_location_mean(df, netcode):
    """
    Get the mean station latitude and longitude coordinates for all stations in a given network.

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param netcode: Network code for which mean coordinates will be returned
    :type netcode: str
    :return: Mean (latitude, longitude) coordinates of stations in the network
    :rtype: tuple(float, float)
    """
    # If a network.station pair appears more than once, we have to make sure we only count it once in the mean.
    # Since each network.station pair might not have the same coordinates, we compute the location of a given
    # network.station as the mean of its records.
    mean_lat = 0.0
    mean_lon = 0.0
    count = 0
    df_net = df.loc[(df['net'].str.upper() == netcode.upper()), ['sta', 'stationLat', 'stationLon']]
    for _, r in df_net.groupby('sta'):
        mean_lat += r['stationLat'].mean()
        mean_lon += r['stationLon'].mean()
        count += 1
    if count:
        mean_lat = mean_lat/count
        mean_lon = mean_lon/count
        return (mean_lat, mean_lon)
    else:
        return (np.nan, np.nan)


def get_network_date_range(df, netcode):
    """
    Get the date range of pick events in df for a given network code

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param netcode: Network code whose pick event dates min and max will be returned
    :type netcode: str
    :return: Min and max dates of picks for given network
    :rtype: tuple(obspy.UTCDateTime, obspy.UTCDateTime)
    """
    mask = (df['net'].str.upper() == netcode.upper())
    if np.any(mask):
        df_net = df.loc[mask]
        min_date = df_net['originTimestamp'].min()
        max_date = df_net['originTimestamp'].max()
        return (obspy.UTCDateTime(min_date), obspy.UTCDateTime(max_date))
    else:
        return (None, None)


def get_station_date_range(df, netcode, statcode):
    """
    Get the date range of pick events in df for a given network and station code

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param netcode: Network code
    :type netcode: str
    :param statcode: Station code
    :type statcode: str
    :return: Min and max dates of picks for given network and station
    :rtype: tuple(obspy.UTCDateTime, obspy.UTCDateTime)
    """
    net_mask = (df['net'].str.upper() == netcode.upper())
    sta_mask = (df['sta'].str.upper() == statcode.upper())
    mask = (net_mask & sta_mask)
    if np.any(mask):
        df_netsta = df.loc[mask]
        min_date = df_netsta['originTimestamp'].min()
        max_date = df_netsta['originTimestamp'].max()
        return (obspy.UTCDateTime(min_date), obspy.UTCDateTime(max_date))
    else:
        return (None, None)


def compute_matching_network_mask(df, net_dict):
    """Compute the mask for df of network codes and station codes that match the sequence
    specified in net_dict.

    This function differs from `pandas.DataFrame.isin()` function in that it only matches net.sta pairs from
    the same index in the 'net' and 'sta' lists from net_dict, whereas `isin()` matches net.sta pairs from
    arbitrary different positions in the lists.

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param net_dict: Lists of corresponding network and station codes to match pariwise against columns 'net and 'sta'.
    :type net_dict: dict of corresponding network and station codes under keys 'net' and 'sta' respectively.
    :return: Row mask for df haev value True in rows whose ['net', 'sta'] columns match one of the ('net', 'sta')
        pairs generated from net_dict.
    :rtype: numpy.array(bool)
    """
    if len(net_dict['net']) != len(net_dict['sta']):
        raise ValueError('Expect net and sta list to be same length')
    _df = df[['net', 'sta']]
    mask = np.array([False]*len(_df))
    for n, s in zip(net_dict['net'], net_dict['sta']):
        mask = (mask | ((_df['net'] == n) & (_df['sta'] == s)))
    return mask


def get_overlapping_date_range(df, network_1, network_2):
    """
    Get the range of dates for which pick events in df from network_1 set of stations overlap with
    any picks from network_2 set of stations.

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param network_1: Network and station codes for first network
    :type network_1: dict of corresponding network and station codes under keys 'net' and 'sta'
    :param network_2: Network and station codes for second network
    :type network_2: dict of corresponding network and station codes under keys 'net' and 'sta'
    :return: Start and end dates of datetime range during which pick events overlap for
        network_1 and network_2
    :rtype: tuple(obspy.UTCDateTime, obspy.UTCDateTime)
    """
    # TODO: Split the first part out into a separate function for determining event IDs of events in common.

    if not network_1 or not network_2:
        return (None, None)

    if 'net' not in network_1 or 'net' not in network_2:
        raise KeyError('Expect \'net\' in dict keys')
    if 'sta' not in network_1 or 'sta' not in network_2:
        raise KeyError('Expect \'sta\' in dict keys')
    # TODO: Change data structures so that such a mismatch as this is not possible.
    if (len(network_1['net']) != len(network_1['sta'])) or (len(network_2['net']) != len(network_2['sta'])):
        raise ValueError('Expect net and sta list to be same length')

    # Find sets of records that match the two input station sets.
    mask_1 = compute_matching_network_mask(df, network_1)
    mask_2 = compute_matching_network_mask(df, network_2)
    mask = (mask_1 | mask_2)
    if not np.any(mask):
        return (None, None)

    # Find which event IDs are common to both sets of records
    ev_1 = set(df.loc[mask_1, "#eventID"].values)
    ev_2 = set(df.loc[mask_2, "#eventID"].values)
    common_events = ev_1 & ev_2
    event_mask = df['#eventID'].isin(list(common_events))
    # Combine the common events with the networks of interest
    event_and_net_mask = (event_mask & (mask_1 | mask_2))
    df_masked = df[event_and_net_mask]

    # Extract the min and max timestamps
    if not df_masked.empty:
        return (obspy.UTCDateTime(df_masked['originTimestamp'].min()),
                obspy.UTCDateTime(df_masked['originTimestamp'].max()))

    return (None, None)


def generate_large_events_catalog(df_picks, min_magnitude=8.0, min_record_count=400, label_historical_events=True):
    """Use input picks dataset to identify dates of large seismic events

    :param df_picks: Dataframe of picks data following PICKS_TABLE_SCHEMA
    :type df_picks: pandas.DataFrame
    :param min_magnitude: Minimum seismic magnitude to be considered 'large' event, defaults to 8.0
    :type min_magnitude: float, optional
    :param min_record_count: Minimum number of records per day required with magnitude >= min_magnitude to be included
        in the event catalog, defaults to 400
    :type min_record_count: int, optional
    :param label_historical_events: Whether to populate fixed list of known historical events, defaults to True
    :type label_historical_events: bool, optional
    :return: Dataframe of large seismic events indexed by date.
    :rtype: pandas.DataFrame
    """
    df_largemag = df_picks[df_picks['mag'] >= min_magnitude]
    temp = df_largemag['originTimestamp'].transform(datetime.datetime.utcfromtimestamp)
    temp = temp.transform(lambda x: x.strftime("%Y-%m-%d"))
    df_largemag = df_largemag.assign(day=temp)
    df_largemag = df_largemag.sort_values(['day', 'originTimestamp'])

    day_largemag_count = [(day, len(df_day)) for day, df_day in df_largemag.groupby('day')]
    dates, counts = zip(*day_largemag_count)
    largemag_dict = {'date': dates, 'counts': counts, 'name': ['unknown']*len(counts)}
    largemag_events_df = pd.DataFrame(largemag_dict, columns=['date', 'counts', 'name'])

    significant_events = largemag_events_df[largemag_events_df['counts'] >= min_record_count]
    significant_events = significant_events.set_index('date')

    if label_historical_events:
        significant_events.loc['2001-06-23', 'name'] = '2001 South Peru Earthquake'
        significant_events.loc['2001-11-14', 'name'] = '2001 Kunlun earthquake'
        significant_events.loc['2002-11-03', 'name'] = '2002 Denali earthquake'
        significant_events.loc['2003-09-25', 'name'] = '2003 Tokachi-Oki earthquake'
        significant_events.loc['2004-12-26', 'name'] = '2004 Indian Ocean earthquake and tsunami'
        significant_events.loc['2005-03-28', 'name'] = '2005 Nias-Simeulue earthquake'
        significant_events.loc['2009-09-29', 'name'] = '2009 Samoa earthquake and tsunami'
        significant_events.loc['2010-02-27', 'name'] = '2010 Chile earthquake'
        significant_events.loc['2011-03-11', 'name'] = '2011 Tohoku earthquake and tsunami'
        significant_events.loc['2012-04-11', 'name'] = '2012 Indian Ocean earthquakes'
        significant_events.loc['2013-02-06', 'name'] = '2013 Solomon Islands earthquakes'
        significant_events.loc['2013-09-24', 'name'] = '2013 Balochistan earthquakes'
        significant_events.loc['2014-04-01', 'name'] = '2014 Iquique earthquake'
        significant_events.loc['2015-09-16', 'name'] = '2015 Illapel earthquake'
        significant_events.loc['2016-08-24', 'name'] = '2016 Myanmar earthquake'

    return significant_events
