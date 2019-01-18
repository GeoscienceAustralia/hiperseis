#!/usr/bin/env python

import textwrap


class Origin:
    """Seismic event origin"""
    def __init__(self, utctime, lat, lon, depthkm):
        self.utctime = utctime
        self.lat = lat
        self.lon = lon
        self.depthkm = depthkm
        self.magnitude_list = []
        self.arrivals = []

    def location(self):
        return (self.lat, self.lon, self.depthkm)

    def epicenter(self):
        return (self.lat, self.lon)

    def __str__(self):
        return "Time: {0}\nLocation: ({1}, {2}, {3})".format(self.utctime, self.lat, self.lon, self.depthkm)


class Event:
    """Seismic event"""
    def __init__(self):
        self.public_id = None
        self.preferred_origin = None
        self.preferred_magnitude = None
        self.origin_list = []

    def __str__(self):
        return "ID: {0}\nPreferred origin:\n{1}\nPreferred Magnitude:\n{2}".format(
               self.public_id,
               textwrap.indent(str(self.preferred_origin), "  "),
               textwrap.indent(str(self.preferred_magnitude), "  "))


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
        self.distance = distance # degrees

    def __str__(self):
        return "{}.{} ({}, {}, {}), {}: {} -> {}, {}".format(
                self.net, self.sta, self.lat, self.lon, self.elev, self.cha, self.phase, self.utctime, self.distance
        )


