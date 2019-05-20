#!/usr/bin/env python
"""
Helper functions for converting between Pandas dataframe and FDSN Inventory,
Network, Station and Channel objects.
"""

# pylint: disable=too-many-locals

from collections import defaultdict
import numpy as np
import pandas as pd

from obspy.core import utcdatetime
from obspy.core.inventory import Network, Station, Channel, Site

from seismic.inventory.table_format import TABLE_COLUMNS, PANDAS_MAX_TIMESTAMP


def pd2Station(statcode, station_df, instrument_register=None):
    """
    Convert Pandas dataframe with unique station code to obspy Station object.

    :param statcode: Station code
    :type statcode: str
    :param station_df: Dataframe containing records for a single station code.
    :type station_df: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    :param instrument_register: Dictionary of nominal instrument responses indexed by channel code, defaults to None
    :param instrument_register: dict of {str, Instrument(obspy.core.inventory.util.Equipment,
        obspy.core.inventory.response.Response)}, optional
    :return: Station object containing the station information from the dataframe
    :rtype: obspy.core.inventory.station.Station
    """
    station_data = station_df.iloc[0]
    st_start = station_data['StationStart']
    assert pd.notnull(st_start)
    st_start = utcdatetime.UTCDateTime(st_start)
    st_end = station_data['StationEnd']
    assert pd.notnull(st_end)
    st_end = utcdatetime.UTCDateTime(st_end)
    station = Station(statcode,
                      station_data['Latitude'],
                      station_data['Longitude'],
                      station_data['Elevation'],
                      start_date=st_start, creation_date=st_start,
                      end_date=st_end, termination_date=st_end,
                      site=Site(name=' '))
    for _, d in station_df.iterrows():
        ch_start = d['ChannelStart']
        ch_start = utcdatetime.UTCDateTime(ch_start) if not pd.isnull(ch_start) else None
        ch_end = d['ChannelEnd']
        ch_end = utcdatetime.UTCDateTime(ch_end) if not pd.isnull(ch_end) else None
        ch_code = d['ChannelCode']
        instrument = instrument_register[ch_code]
        if instrument is not None:
            sensor = instrument.sensor
            response = instrument.response
        elif 'LAST_RESORT' in instrument_register:
            last_resort = instrument_register['LAST_RESORT']
            sensor = last_resort.sensor
            response = last_resort.response
        else:
            sensor = None
            response = None
        cha = Channel(ch_code, '', float(d['Latitude']), float(d['Longitude']), float(d['Elevation']),
                      depth=0.0, azimuth=0.0, dip=-90.0, sample_rate=0.0, clock_drift_in_seconds_per_sample=0.0,
                      start_date=ch_start, end_date=ch_end, sensor=sensor, response=response)
        station.channels.append(cha)
    return station


def pd2Network(netcode, network_df, instrument_register, progressor=None):
    """
    Convert Pandas dataframe with unique network code to obspy Network object.

    :param netcode: Network code
    :type netcode: str
    :param network_df: Dataframe containing records for a single network code.
    :type network_df: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    :param instrument_register: Dictionary of nominal instrument responses indexed by channel code, defaults to None
    :param instrument_register: dict of {str, Instrument(obspy.core.inventory.util.Equipment,
        obspy.core.inventory.response.Response)}, optional
    :param progressor: Progress bar functor to receive progress updates, defaults to None
    :param progressor: Callable object receiving incremental update on progress, optional
    :return: Network object containing the network information from the dataframe
    :rtype: obspy.core.inventory.network.Network
    """
    net = Network(netcode, stations=[], description=' ')
    for statcode, ch_data in network_df.groupby('StationCode'):
        station = pd2Station(statcode, ch_data, instrument_register)
        net.stations.append(station)
        if progressor:
            progressor(len(ch_data))
    return net


def inventory2Dataframe(inv_object, show_progress=True):
    """
    Convert a obspy Inventory object to a Pandas Dataframe.

    :param inv_object: Obspy inventory object to convert to dataframe
    :type inv_object: obspy.core.inventory.inventory.Inventory
    :param show_progress: Whether to use a progress bar, defaults to True
    :param show_progress: bool, optional
    :return: Pandas Dataframe with sequential integer index and sorted by [NetworkCode, StationCode].
        Only populates entries for non-empty channels.
    :rtype: pandas.DataFrame conforming to table_format.TABLE_SCHEMA
    """
    if show_progress:
        import tqdm
        num_entries = sum(len(station.channels) for network in inv_object.networks for station in network.stations)
        pbar = tqdm.tqdm(total=num_entries, ascii=True)

    max_end_timestamp = np.datetime64(str(PANDAS_MAX_TIMESTAMP), 's')
    d = defaultdict(list)
    for network in inv_object.networks:
        for station in network.stations:
            if show_progress:
                pbar.update(len(station.channels))
            for channel in station.channels:
                d['NetworkCode'].append(network.code)
                d['StationCode'].append(station.code)
                lat = channel.latitude if channel.latitude else station.latitude
                lon = channel.longitude if channel.longitude else station.longitude
                ele = channel.elevation if channel.elevation else station.elevation
                d['Latitude'].append(lat)
                d['Longitude'].append(lon)
                d['Elevation'].append(ele)
                d['StationStart'].append(np.datetime64(station.start_date))
                d['StationEnd'].append(min(np.datetime64(station.end_date, 's'), max_end_timestamp))
                d['ChannelCode'].append(channel.code)
                d['ChannelStart'].append(np.datetime64(channel.start_date))
                d['ChannelEnd'].append(min(np.datetime64(channel.end_date, 's'), max_end_timestamp))
    if show_progress:
        pbar.close()
    df = pd.DataFrame.from_dict(d)
    df = df[list(TABLE_COLUMNS)]
    df.sort_values(['NetworkCode', 'StationCode'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df
