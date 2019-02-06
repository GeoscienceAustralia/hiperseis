#!/usr/bin/env python
from __future__ import print_function

import sys
import datetime
# textwrap is used to get dynamic level of indentation, depending on the parentage of the object
# being printed.
import textwrap
from obspy import UTCDateTime

if sys.version_info[0] < 3:
    indentprint = lambda x: x
else:
    indentprint = lambda x: textwrap.indent(x, "  ")  # pylint: disable=no-member


class Origin:
    """
    Container for seismic event origin (location) within the earth
    """

    def __init__(self, utctime, lat, lon, depthkm):
        self.utctime = utctime
        self.lat = lat
        self.lon = lon
        self.depthkm = depthkm
        self.magnitude_list = []
        self.arrivals = []

    def location(self):
        """
        Get the location attribute as (lat, long, depth). Lat and long are in degrees, depth is in km.

        :return: Location of the seismic event origin
        :rtype: tuple(float, float, float)
        """
        return (self.lat, self.lon, self.depthkm)

    def epicenter(self):
        """
        Get the epicenter attribute as (lat, long). Lat and long are in degrees.

        :return: Location of the seismic event epicenter
        :rtype: tuple(double, double)
        """
        return (self.lat, self.lon)

    def __str__(self):
        return "Time: {0}\nLocation: ({1}, {2}, {3})".format(self.utctime, self.lat, self.lon, self.depthkm)


class Event:
    """
    Container for a seismic event with high level attributes.
    """

    def __init__(self):
        self.public_id = None
        self.preferred_origin = None
        self.preferred_magnitude = None
        self.origin_list = []

    def __str__(self):
        return "ID: {0}\nPreferred origin:\n{1}\nPreferred Magnitude:\n{2}".format(
               self.public_id,
               indentprint(str(self.preferred_origin)),
               indentprint(str(self.preferred_magnitude)))


class Magnitude:
    """Seismic event magnitude"""
    def __init__(self, mag, mag_type):
        self.magnitude_value = mag
        self.magnitude_type = mag_type

    def __str__(self):
        return "Magnitude: {0} ({1})".format(self.magnitude_value, self.magnitude_type)


class Arrival:
    """Arrival of seismic event signal at other location"""
    def __init__(self, net, sta, loc, cha, lon, lat, elev, phase, utctime, distance):
        self.net = net
        self.sta = sta
        self.loc = loc
        self.cha = cha
        self.lon = lon
        self.lat = lat
        self.elev = elev
        self.phase = phase
        self.utctime = utctime
        self.distance = distance  # degrees

    def __str__(self):
        return "{}.{} ({}, {}, {}), {}: {} -> {}, {}".\
            format(self.net, self.sta, self.lat, self.lon, self.elev, self.cha, self.phase, self.utctime, self.distance)


if __name__ == "__main__":
    # Test code by instantiating each type.
    utctime = UTCDateTime()
    arr = Arrival('GE', 'KUM', '', 'BHZ', 50.0, -5.0, 0.0, 'P', utctime, 55.0)
    mag = Magnitude(7.2, "mb")
    event_time = utctime - datetime.timedelta(minutes=25.0)
    origin = Origin(event_time, 20.0, 80.0, 60.0)
    event = Event()
    event.public_id = 'test'
    event.preferred_origin = origin
    event.preferred_magnitude = mag
    print(event)
    print(arr)
