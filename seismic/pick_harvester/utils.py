from mpi4py import MPI
from lxml import etree as ET
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

def recursive_glob(treeroot, pattern):
    results = []
    for base, dirs, files in os.walk(treeroot):
        goodfiles = fnmatch.filter(files, pattern)
        results.extend(os.path.join(base, f) for f in goodfiles)
    return results
# end func

def split_list(lst, npartitions):
    result = []
    for i in np.arange(npartitions):
        result.append([])
    # end for
    count = 0
    for iproc in np.arange(npartitions):
        for i in np.arange(np.divide(len(lst), npartitions)):
            result[iproc].append(lst[count])
            count += 1
    # end for
    for iproc in np.arange(np.mod(len(lst), npartitions)):
        result[iproc].append(lst[count])
        count += 1
    # end for

    return result
# end func

class Origin:
    def __init__(self, utctime, lat, lon, depthkm):
        self.utctime = utctime
        self.lat = lat
        self.lon = lon
        self.depthkm = depthkm
        self.magnitude_list = []
    # end func
# end class

class Event:
    def __init__(self):
        self.public_id = None
        self.preferred_origin = None
        self.preferred_magnitude = None
        self.origin_list = []
    # end func
# end class

class Magnitude:
    def __init__(self, mag, mag_type):
        self.magnitude_value = mag
        self.magnitude_type = mag_type
    # end func
# end class

class EventParser:
    def __init__(self, xml_filename):
        self.xml_filename = xml_filename
        self.events  = []
        self.origins = []
        self.magnitudes = []

        for xmlevent, elem in ET.iterparse(xml_filename):
            tag = ET.QName(elem.tag).localname
            if (tag.upper() == 'EVENT'):
                self.events.append(elem)
            elif (tag.upper() == 'ORIGIN'):
                self.origins.append(elem)
            elif (tag.upper() == 'MAGNITUDE'):
                self.magnitudes.append(elem)
            # end if
        # end for

        self.origin_dict = defaultdict(list)
        self.magnitude_dict = defaultdict(list)

        # parse magnitudes
        for m in self.magnitudes:
            pid = None
            try:
                pid = int(re.split("[/=]", m.get('publicID'))[-1])
            except:
                pass
            # end try
            if(pid): self.magnitude_dict[pid] = self.parseMagnitude(m)
        # end for

        # parse origins
        for o in self.origins:
            pid = int(re.split("[/=]", o.get('publicID'))[-1])
            self.origin_dict[pid] = self.parseOrigin(o)
            # end for
        assert len(self.origins) == len(self.origin_dict), 'Invalid originIDs found..'

        self.event_dict = defaultdict(list)
        # parse events
        for e in self.events:
            pid = int(re.split("[/=]", e.get('publicID'))[-1])
            self.event_dict[pid] = self.parseEvent(e)
        # end for
        assert len(self.events) == len(self.event_dict), 'Invalid eventIDs found..'

    # end func

    def parseEvent(self, e):
        resultEvent = Event()
        origList = []

        resultEvent.public_id = e.get('publicID')
        for element in e:
            if ('preferredOriginID' in element.tag):
                oid = int(re.split("[/=]", element.text)[-1])
                resultEvent.preferred_origin = self.origin_dict[oid]
            elif('preferredMagnitudeID' in element.tag):
                mid = int(re.split("[/=]", element.text)[-1])
                resultEvent.preferred_magnitude = self.magnitude_dict[mid]
            elif ('origin' in element.tag and \
                  'Reference' not in element.tag):
                oid = int(re.split("[/=]", element.get('publicID'))[-1])
                origList.append(self.origin_dict[oid])
            # end if
        # end for
        resultEvent.origin_list = origList
        return resultEvent
    # end func

    def parseOrigin(self, o):
        utctime = None
        lat = None
        lon = None
        depthkm = None
        pref_mag = None
        magList = []
        for element in o:
            if ('time' in element.tag):
                for subelement in element:
                    if ('value' in subelement.tag):
                        utctime = UTCDateTime(subelement.text)
                        # end if
                        # end for
            elif ('latitude' in element.tag):
                for subelement in element:
                    if ('value' in subelement.tag):
                        lat = float(subelement.text)
                    # end if
                # end for
            # end if
            elif ('longitude' in element.tag):
                for subelement in element:
                    if ('value' in subelement.tag):
                        lon = float(subelement.text)
                    # end if
                # end for
            # end if
            elif ('depth' in element.tag):
                for subelement in element:
                    if ('value' in subelement.tag):
                        depthkm = float(subelement.text)
                    # end if
                # end for
            elif ('magnitude' in element.tag and \
                          'station' not in element.tag):
                mid = None
                try:
                    mid = int(re.split("[/=]", element.get('publicID'))[-1])
                except:
                    pass
                # end try
                if(mid): magList.append(self.magnitude_dict[mid])
            # end if
        # end for

        assert None not in [utctime, lat, lon], 'Failed to find required values for Origin'

        o = Origin(utctime, lat, lon, depthkm)
        o.magnitude_list = magList
        return o
    # end func

    def parseMagnitude(self, m):
        mag_value = None
        mag_type = None
        for element in m:
            if ('magnitude' in element.tag):
                for subelement in element:
                    if ('value' in subelement.tag):
                        mag_value = float(subelement.text)
                    # end if
                # end for
            elif ('type' in element.tag):
                mag_type = element.text
            # end if
        # end for

        return Magnitude(mag_value, mag_type)
    # end func

    def getEvents(self):
        return self.event_dict
    # end func
# end class

class Catalog():
    def __init__(self, event_folder):
        self.event_folder = event_folder
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        # retrieve list of all event xml files
        self.xml_files = recursive_glob(self.event_folder, '*.xml')
        shuffle(self.xml_files) # shuffle to avoid concentration of large files
        self.proc_workload = split_list(self.xml_files, self.nproc)

        # broadcast workload to all procs
        self.proc_workload = self.comm.bcast(self.proc_workload, root=0)

        print 'Rank %d: processing %d files' % (self.rank, len(self.proc_workload[self.rank]))

        self._load_events()
    # end func

    def _load_events(self):
        eventList = []
        poTimestamps = []
        for ifn, fn in enumerate(self.proc_workload[self.rank]):
            es = EventParser(fn).getEvents()
            for i, (eid, e) in enumerate(es.iteritems()):
                po = e.preferred_origin
                if (not (po.depthkm >= 0)): continue

                eventList.append(e)
                poTimestamps.append(po.utctime.timestamp)
            # end for
        # end for

        eventList = self.comm.allgather(eventList)
        poTimestamps = self.comm.allgather(poTimestamps)

        allEventList = []
        allPOTimestamps = []
        for iproc in np.arange(self.nproc):
            for i, e in enumerate(eventList[iproc]):
                allEventList.append(e)
                allPOTimestamps.append(poTimestamps[iproc][i])
            # end for
        # end for

        self.allPOTimestamps = np.array(allPOTimestamps)
        self.allEventList = allEventList

        if (self.rank == 0):
            print 'Collected %d event origins' % (len(self.allEventList))

            hasPM = 0
            hasMultipleMags = 0
            for e in self.allEventList:
                o = e.preferred_origin
                if (e.preferred_magnitude): hasPM += 1
                if (len(o.magnitude_list)): hasMultipleMags += 1
            # end for

            print '%d preferred origins have a preferred magnitude' % (hasPM)
            print '%d preferred origins have at least one magnitude' % (hasMultipleMags)
        # end if
    # end func

    def get_events(self):
        return self.allEventList
    # end func

    def get_preferred_origin_timestamps(self):
        return self.allPOTimestamps
    # end func
# end class

class CatalogCSV:
    def __init__(self, event_folder):
        self.event_folder = event_folder
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        self.event_folder = event_folder

        # retrieve list of all csv files
        self.csv_files = sorted(recursive_glob(self.event_folder, '*.csv'))
        self._load_events()
    # end func

    def _load_events(self):
        eventList = []
        poTimestamps = []

        if(self.rank==0):
            for ifn, fn in enumerate(self.csv_files):
                print 'Reading %s' % (fn)

                for line in open(fn, 'r'):
                    if(line[0]=='#'):
                        items = line.split(',')
                        vals = map(float, items[1:])

                        year = int(vals[0])
                        month = int(vals[1])
                        day = int(vals[2])
                        hour = int(vals[3] if vals[3] >=0 else 0)
                        minute = int(vals[4] if vals[4] >=0 else 0)
                        second = vals[5] if vals[5] >=0 else 0

                        lon = vals[6]
                        lat = vals[7]
                        if (lon <- 180 or lon > 180): continue
                        if (lat <- 90 or lat > 90): continue

                        depth = vals[8] if vals[8] >=0 else 0

                        mb = vals[10]
                        ms = vals[11]
                        mi = vals[12]
                        mw = vals[13]
                        mag = 0
                        magtype='mw'
                        if(mw>0):
                            mag = mw
                            magtype='mw'
                        elif(ms>0):
                            mag = ms
                            magtype = 'ms'
                        elif(mb>0):
                            mag = mb
                            magtype = 'mb'
                        elif(mi>0):
                            mag = mi
                            magtype = 'mi'
                        # end if

                        eventid = vals[-1]

                        utctime = None
                        try:
                            utctime = UTCDateTime(year, month, day,
                                                  hour, minute, second)
                        except:
                            continue
                        # end try

                        origin = Origin(utctime, lat, lon, depth)
                        event = Event()
                        event.public_id = eventid
                        event.preferred_origin = origin
                        event.preferred_magnitude = Magnitude(mag, magtype)

                        eventList.append(event)
                        poTimestamps.append(origin.utctime.timestamp)
                    # end if
                # end for
            # end for
        # end if

        self.allEventList = self.comm.bcast(eventList, root=0)
        self.allPOTimestamps = self.comm.bcast(poTimestamps, root=0)
        self.allPOTimestamps = np.array(self.allPOTimestamps)
    # end func

    def get_events(self):
        return self.allEventList
    # end func

    def get_preferred_origin_timestamps(self):
        return self.allPOTimestamps
    # end func
# end class