#!/usr/bin/env python
"""Utility functions and constants shared amongst inventory management modules.
"""

import sys

from obspy import read_inventory

PY2 = sys.version_info[0] < 3  # pylint: disable=invalid-name

NOMINAL_EARTH_RADIUS_KM = 6378.1370  # pylint: disable=invalid-name

SORT_ORDERING = ['NetworkCode', 'StationCode', 'StationStart', 'StationEnd',  # pylint: disable=invalid-name
                 'ChannelCode', 'ChannelStart', 'ChannelEnd']


def load_station_xml(inventory_file):
    """Load a stationxml file

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

    return obspy_inv
