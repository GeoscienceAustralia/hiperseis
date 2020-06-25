#!/usr/bin/env python
"""
Lightweight reader for CSV seismic event catalogs, indexing the found events by event ID and station.

This was adapted for the speical use case of distance-to-event QA checks performed for ticket PST-340.
"""

import os
import fnmatch
from collections import defaultdict
import random as rnd

from obspy import UTCDateTime
from tqdm.auto import tqdm

from seismic.inventory.dataio.event_attrs import Origin, Event, Magnitude, Arrival


def recursive_glob(treeroot, pattern):
    """
    Generate a complete list of files matching pattern under the root of a directory hierarchy.

    :param treeroot: Path to the root of the directory tree.
    :type treeroot: str or pathlib.Path
    :param pattern: File name pattern to match, e.g. "\*.csv"
    :type pattern: str
    :return: List of paths to the files matching the pattern, qualified relative to treeroot
    :rtype: list(str)
    """
    results = []
    for base, dirs, files in os.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results


class CatalogCSV:
    """
    Lightweight parser for seismic event catalog.

    Catalog is format as follows::

       #EHB, 2005, 09, 16, 07, 28, 39.001,  126.93300,    4.18700,    2.90000,    28, 4.50, -999.00, -999.00, -999.00,       1, 134.3000, 1
       FITZ, BHZ,  ,  ,  ,  ,  , P , 2005, 09, 16, 07, 33, 37.00,  22.180
       WR0 , BHZ,  ,  ,  ,  ,  , P , 2005, 09, 16, 07, 34, 06.00,  25.130
       KUM , BHZ,  ,  ,  ,  ,  , P , 2005, 09, 16, 07, 35, 02.00,  26.220
       MEEK, BHZ,  ,  ,  ,  ,  , P , 2005, 09, 16, 07, 35, 04.00,  31.680
       FORT, BHZ,  ,  ,  ,  ,  , P , 2005, 09, 16, 07, 35, 32.00,  34.780
       STKA, BHZ,  ,  ,  ,  ,  , P , 2005, 09, 16, 07, 36, 04.00,  38.480
       KSM , BHZ,  ,  ,  ,  ,  , Pn, 2005, 09, 16, 07, 32, 39.00,  16.820
       KAKA, BHZ,  ,  ,  ,  ,  , P , 2005, 09, 16, 07, 33, 04.00,  17.660
       FITZ, BHZ,  ,  ,  ,  ,  , P , 2005, 09, 16, 07, 33, 36.00,  22.180
       ...

    The header row, which starts with '#', indicates the number of phases listed in the subsequent lines of arrival data
    (28 in this example).
    """

    def __init__(self, event_folder, sampling_factor=1.0):
        self.event_folder = event_folder
        self.sampling_factor = sampling_factor  # Proportion of events that get sampled
        # self.comm = MPI.COMM_WORLD
        # self.nproc = self.comm.Get_size()
        # self.rank = self.comm.Get_rank()

        # retrieve list of all csv files
        self.csv_files = sorted(recursive_glob(self.event_folder, '*.csv'))
        self._load_events()

    def _load_events(self):
        event_dict = {}
        station_dict = defaultdict(lambda: defaultdict(list))

        event_id = None
        if True:  # (self.rank==0):
            for ifn, fn in enumerate(self.csv_files):
                print('Reading %s' % (fn))

                progress = tqdm(total=os.path.getsize(fn), desc=fn)

                for line in open(fn, 'r'):
                    progress.update(len(line))
                    if line[0] == '#':
                        try:
                            if rnd.random() < self.sampling_factor:
                                event_id, event = self._parse_event_header(line)
                                event_dict[event_id] = event
                            else:
                                event_id = None
                        except:
                            event_id = None
                            continue
                    elif event_id is not None:
                        try:
                            arrival = self._parse_arrival(line)
                            event_dict[event_id].preferred_origin.arrivals.append(arrival)
                            station_dict[arrival.sta][arrival.phase].append((event_id, arrival.distance))
                        except:
                            continue
                progress.close()

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
        magtype = 'mw'
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

        eventid = int(items[15].strip())

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


if __name__ == "__main__":
    # Test path for standalone execution.
    self_path = os.path.dirname(os.path.abspath(__file__))
    event_src_folder = os.path.join(self_path, 'test')
    cat = CatalogCSV(event_src_folder)
    for e in cat.get_events():
        print(e)
