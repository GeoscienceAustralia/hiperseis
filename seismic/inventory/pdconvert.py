#!/usr/bin/env python
"""Helper functions for converting between Pandas dataframe and FDSN Inventory,
   Network, Station and Channel objects.
"""

from collections import defaultdict
import numpy as np
import pandas as pd

from obspy import read_inventory
from obspy.core import utcdatetime
from obspy.core.inventory import Inventory, Network, Station, Channel, Site

from table_format import TABLE_SCHEMA, TABLE_COLUMNS, PANDAS_MAX_TIMESTAMP


def pd2Station(statcode, station_df, instrument_register=None):
    """Convert Pandas dataframe with unique station code to FDSN Station object."""
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
        else:
            sensor = None
            response = None
        cha = Channel(ch_code, '', float(d['Latitude']), float(d['Longitude']), float(d['Elevation']),
                      depth=0.0, azimuth=0.0, dip=-90.0,
                      start_date=ch_start, end_date=ch_end,
                      sensor=sensor, response=response)
        station.channels.append(cha)
    return station


def pd2Network(netcode, network_df, instrument_register, progressor=None):
    """Convert Pandas dataframe with unique network code to FDSN Network object."""
    net = Network(netcode, stations=[], description=' ')
    for statcode, ch_data in network_df.groupby('StationCode'):
        station = pd2Station(statcode, ch_data, instrument_register)
        net.stations.append(station)
        if progressor:
            progressor(len(ch_data))
    return net


def inventory2Dataframe(inv_object, show_progress=True):
    """Convert a obspy Inventory object to a Pandas Dataframe.

       Returns a Pandas Dataframe with sequential integer index and sorted by [NetworkCode, StationCode].
       Only populates entries for non-empty channels.
    """

    if show_progress:
        import tqdm
        num_entries = sum(len(station.channels) for network in inv_object.networks for station in network.stations)
        pbar = tqdm.tqdm(total=num_entries, ascii=True)

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
                d['StationEnd'].append(np.datetime64(station.end_date))
                d['ChannelCode'].append(channel.code)
                d['ChannelStart'].append(np.datetime64(channel.start_date))
                d['ChannelEnd'].append(np.datetime64(channel.end_date))
    if show_progress:
        pbar.close()
    df = pd.DataFrame.from_dict(d)
    df = df[list(TABLE_COLUMNS)]
    df.sort_values(['NetworkCode', 'StationCode'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df
