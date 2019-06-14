#!/usr/bin/env python
"""
Utility functions supporting plotting for cross-correlation visualizations.
"""

import math
import datetime
import numpy as np

from shapely.geometry import Polygon
from descartes import PolygonPatch


def distance(origin, destination):
    """
    Compute the distance in km between origin coordinates and destination coordinates.
    The coordinates are (latitude, longitude) couplets in units of degrees.

    :param origin: Coordinates of origin point
    :type origin: tuple(float, float)
    :param destination: Coordinates of destination point
    :type destination: tuple(float, float)
    :return: Epicentral distance between origin and destination in kilometres
    :rtype: float
    """
    lat1, lon1 = origin
    lat2, lon2 = destination
    radius = 6371  # km

    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) * math.sin(dlat / 2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon / 2) * math.sin(dlon / 2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = radius * c

    return d


def drawBBox(min_lon, min_lat, max_lon, max_lat, base_map, **kwargs):
    """
    Draw bounding box on a basemap

    :param min_lon: Minimum longitude
    :type min_lon: float
    :param min_lat: Minimum latitude
    :type min_lat: float
    :param max_lon: Maximum longitude
    :type max_lon: float
    :param max_lat: Maximum latitude
    :type max_lat: float
    :param base_map: Basemap on which to draw the bounding box
    :type base_map: mpl_toolkits.basemap.Basemap
    """
    bblons = np.array([min_lon, max_lon, max_lon, min_lon, min_lon])
    bblats = np.array([min_lat, min_lat, max_lat, max_lat, min_lat])

    x, y = base_map(bblons, bblats)
    xy = zip(x, y)
    poly = Polygon(xy)
    base_map.ax.add_patch(PolygonPatch(poly, **kwargs))


def timestamps_to_plottable_datetimes(time_series):
    """
    Convert a series of float (or equivalent) timestamp values to matplotlib plottable datetimes.

    :param time_series: Series of timestamps
    :type time_series: iterable container
    :return: Equivalent series of plottable timestamps
    :rtype: numpy.array('datetime64[ms]') with millisecond resolution
    """
    plt_times = np.array([datetime.datetime.utcfromtimestamp(v)
                          for v in time_series]).astype('datetime64[ms]')
    return plt_times
