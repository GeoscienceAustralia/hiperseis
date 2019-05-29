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
    mean_lat = df[df['net'] == netcode]['stationLat'].mean()
    mean_lon = df[df['net'] == netcode]['stationLon'].mean()
    return (mean_lat, mean_lon)


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
    mask = (df['net'] == netcode)
    df_net = df.loc[mask]
    min_date = df_net['originTimestamp'].min()
    max_date = df_net['originTimestamp'].max()
    return (obspy.UTCDateTime(min_date), obspy.UTCDateTime(max_date))


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
    mask = (df['net'] == netcode)
    df_net = df.loc[mask]
    mask = (df_net['sta'] == statcode)
    df_sta = df_net.loc[mask]
    min_date = df_sta['originTimestamp'].min()
    max_date = df_sta['originTimestamp'].max()
    return (obspy.UTCDateTime(min_date), obspy.UTCDateTime(max_date))


def get_overlapping_date_range(df, ref_network, target_network):
    """
    Get the range of dates for which pick events in df from ref_network overlap with
    any picks from target_network

    :param df: Pandas dataframe loaded from a pick event ensemble
    :type df: pandas.DataFrame
    :param ref_network: Network and station codes for reference network
    :type ref_network: dict of corresponding network and station codes under keys 'net' and 'sta'
    :param target_network: Network and station codes for target network
    :type target_network: dict of corresponding network and station codes under keys 'net' and 'sta'
    :return: Start and end dates of datetime range during which pick events overlap for
        ref_network and target_network
    :rtype: tuple(obspy.UTCDateTime, obspy.UTCDateTime)
    """
    mask_ref = df[list(ref_network)].isin(ref_network).all(axis=1)
    mask_targ = df[list(target_network)].isin(target_network).all(axis=1)
    mask = mask_ref | mask_targ
    if not np.any(mask):
        return (None, None)
    df_nets = df.loc[mask]
    keep_events = [e for e, d in df_nets.groupby('#eventID')
                   if np.any(d[list(ref_network)].isin(ref_network).all(axis=1)) and
                   np.any(d[list(target_network)].isin(target_network).all(axis=1))]
    event_mask = df_nets['#eventID'].isin(keep_events)
    df_nets = df_nets[event_mask]
    return (obspy.UTCDateTime(df_nets['originTimestamp'].min()), obspy.UTCDateTime(df_nets['originTimestamp'].max()))


def generate_large_events_catalog(df_picks, min_magnitude=8.0, label_historical_events=True):
    """Use input picks dataset to identify dates of large seismic events

    :param df_picks: Dataframe of picks data following PICKS_TABLE_SCHEMA
    :type df_picks: pandas.DataFrame
    :param min_magnitude: Minimum seismic magnitude to be considered 'large' event, defaults to 8.0
    :type min_magnitude: float, optional
    :param label_historical_events: Whether to populate fixed list of known historical events, defaults to True
    :type label_historical_events: bool, optional
    :return: Dataframe of large seismic events indexed by date.
    :rtype: pandas.DataFrame
    """
    df_largemag = df_picks[df_picks['mag'] >= min_magnitude]
    df_largemag['day'] = df_largemag['originTimestamp'].transform(datetime.datetime.utcfromtimestamp)\
        .transform(lambda x: x.strftime("%Y-%m-%d"))
    df_largemag = df_largemag.sort_values(['day', 'originTimestamp'])

    day_largemag_count = [(day, len(df_day)) for day, df_day in df_largemag.groupby('day')]
    dates, counts = zip(*day_largemag_count)
    largemag_dict = {'date': dates, 'counts': counts, 'name': ['unknown']*len(counts)}
    largemag_events_df = pd.DataFrame(largemag_dict, columns=['date', 'counts', 'name'])

    anonymous_event_count_threshold = 400
    significant_events = largemag_events_df[largemag_events_df['counts'] >= anonymous_event_count_threshold]
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
