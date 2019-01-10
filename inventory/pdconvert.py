#!/usr/bin/env python
"""Helper functions for converting Pandas dataframe to FDSN Inventory, Network, Station and Channel objects.
"""

import pandas as pd

from obspy.core import utcdatetime
from obspy.core.inventory import Inventory, Network, Station, Channel, Site

DEFAULT_START_TIMESTAMP = pd.Timestamp("1964-1-1 00:00:00")
DEFAULT_END_TIMESTAMP = pd.Timestamp.max


def pd2Station(statcode, station_df):
    """Convert Pandas dataframe with unique station code to FDSN Station object."""
    station_data = station_df.iloc[0]
    st_start = station_data['StationStart']
    st_start = utcdatetime.UTCDateTime(st_start) if not pd.isnull(st_start) else DEFAULT_START_TIMESTAMP
    st_end = station_data['StationEnd']
    st_end = utcdatetime.UTCDateTime(st_end) if not pd.isnull(st_end) else DEFAULT_END_TIMESTAMP
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
        cha = Channel(d['ChannelCode'], '', float(d['Latitude']), float(d['Longitude']), float(d['Elevation']),
                      depth=0.0, azimuth=0.0, dip=-90.0,
                      start_date=ch_start,
                      end_date=ch_end)
        station.channels.append(cha)
    return station


def pd2Network(netcode, network_df, progressor=None):
    """Convert Pandas dataframe with unique network code to FDSN Network object."""
    net = Network(netcode, stations=[], description=' ')
    for statcode, ch_data in network_df.groupby('StationCode'):
        station = pd2Station(statcode, ch_data)
        net.stations.append(station)
        if progressor:
            progressor(len(ch_data))
    return net
