#from mpi4py import MPI
from os.path import join, exists
import os, glob, fnmatch, sys
import re
from collections import defaultdict
from obspy import UTCDateTime
import numpy as np
import math

from math import radians, cos, sin, asin, sqrt
import numpy as np
import scipy
from scipy.spatial import cKDTree
from random import shuffle
import bisect
import textwrap
from tqdm.auto import tqdm


def recursive_glob(treeroot, pattern):
    results = []
    for base, dirs, files in os.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results


class Origin:
    def __init__(self, utctime, lat, lon, depthkm):
        self.utctime = utctime
        self.lat = lat
        self.lon = lon
        self.depthkm = depthkm
        self.magnitude_list = []
        self.arrivals = []

    def location(self):
        return (lat, lon, depthkm)

    def __str__(self):
        return "Time: {0}\nLocation: ({1}, {2}, {3})".format(self.utctime, self.lat, self.lon, self.depthkm)


class Event:
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
    def __init__(self, mag, mag_type):
        self.magnitude_value = mag
        self.magnitude_type = mag_type

    def __str__(self):
        return "Magnitude: {0} ({1})".format(self.magnitude_value, self.magnitude_type)


class Arrival:
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
        self.distance = distance

    def __str__(self):
        return "{}.{} ({}, {}, {}), {}: {} -> {}, {}".format(
                self.net, self.sta, self.lat, self.lon, self.elev, self.cha, self.phase, self.utctime, self.distance
        )


class CatalogCSV:
    def __init__(self, event_folder):
        self.event_folder = event_folder
        # self.comm = MPI.COMM_WORLD
        # self.nproc = self.comm.Get_size()
        # self.rank = self.comm.Get_rank()

        self.event_folder = event_folder

        # retrieve list of all csv files
        self.csv_files = sorted(recursive_glob(self.event_folder, '*.csv'))
        self._load_events()

    def _load_events(self):
        event_dict = defaultdict(list)
        station_dict = defaultdict(list)

        event_id = None
        if True: # (self.rank==0):
            for ifn, fn in enumerate(self.csv_files):
                print('Reading %s' % (fn))

                progress = tqdm(total=os.path.getsize(fn), desc=fn)

                for line in open(fn, 'r'):
                    progress.update(len(line))
                    if (line[0] == '#'):
                        try:
                            event_id, event = self._parse_event_header(line)
                            event_dict[event_id] = event
                        except:
                            event_id = None
                            continue
                    else:
                        assert event_id is not None
                        try:
                            arrival = self._parse_arrival(line)
                            event_dict[event_id].preferred_origin.arrivals.append(arrival)
                            bisect.insort(station_dict[arrival.sta], event_id)
                        except:
                            continue

                    # end if
                # end for
                progress.close()
            # end for
        # end if

        self.event_dict = event_dict
        self.station_dict = station_dict
    # end func


    def _parse_event_header(self, line):
        items = line.split(',')
        vals = list(map(float, items[1:]))

        year = int(vals[0])
        month = int(vals[1])
        day = int(vals[2])
        hour = int(vals[3] if vals[3] >=0 else 0)
        minute = int(vals[4] if vals[4] >=0 else 0)
        second = vals[5] if vals[5] >=0 else 0

        lon = vals[6]
        lat = vals[7]
        if (lon < -180 or lon > 180):
            raise Exception
        if (lat < -90 or lat > 90):
            raise Exception

        depth = vals[8] if vals[8] >=0 else 0

        mb = vals[10]
        ms = vals[11]
        mi = vals[12]
        mw = vals[13]
        mag = 0
        magtype='mw'
        if mw > 0:
            mag = mw
            magtype='mw'
        elif ms > 0:
            mag = ms
            magtype = 'ms'
        elif mb > 0:
            mag = mb
            magtype = 'mb'
        elif mi > 0:
            mag = mi
            magtype = 'mi'

        eventid = int(items[-1].strip())

        utctime = None
        utctime = UTCDateTime(year, month, day, hour, minute, second)

        origin = Origin(utctime, lat, lon, depth)
        event = Event()
        event.public_id = eventid
        event.preferred_origin = origin
        event.preferred_magnitude = Magnitude(mag, magtype)

        return eventid, event
            

    def _parse_arrival(self, line):
        items = line.split(',')
        vals = list(map(float, items[8:]))
        year = int(vals[0])
        month = int(vals[1])
        day = int(vals[2])
        hour = int(vals[3])
        minute = int(vals[4])
        second = vals[5]

        utctime = UTCDateTime(year, month, day, hour, minute, second)

        try:
            lon = float(items[4])
        except:
            lon = 0

        try:
            lat = float(items[5])
        except:
            lat = 0

        try:
            elev = float(items[6])
        except:
            elev = 0

        distance = vals[-1]
        netcode = items[3].strip()
        statcode = items[0].strip()
        arrival = Arrival(netcode, statcode, items[2].strip(), items[1].strip(),
                          lon, lat, elev, items[7].strip(), utctime, distance)
        return arrival


    def get_events(self):
        return self.event_dict.values()

# end class

#event_src_folder = '/g/data/ha3/am7399/temp'
#cat = CatalogCSV(event_src_folder)
