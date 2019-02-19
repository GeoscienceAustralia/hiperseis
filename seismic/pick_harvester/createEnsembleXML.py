#!/bin/env python
"""
Description:
    This script was initially written for inserting new picks into the ISC catalogue.
    We now use a unified csv catalogue (that Babak has prepared) and this script merges existing
    picks with those picked by our parallel picker and creates self-consistent SC3ML files 
    to be ingested into Seiscomp3.
    
CreationDate:   20/11/18
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     20/11/18   RH

"""

import click
import os, glob, fnmatch, sys
from random import shuffle
from obspy import UTCDateTime, read_events, read_inventory
from obspy.taup.taup_geo import calc_dist
from obspy.clients.iris import Client as IrisClient
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
from obspy.signal.trigger import trigger_onset, z_detect, classic_sta_lta, recursive_sta_lta, ar_pick
from obspy.signal.rotate import rotate_ne_rt
from obspy.core.event import Pick as OPick, \
    CreationInfo as OCreationInfo, \
    WaveformStreamID as OWaveformStreamID, \
    ResourceIdentifier as OResourceIdentifier, \
    Arrival as OArrical, Event as OEvent,\
     Origin as OOrigin, Arrival as OArrival, \
    OriginQuality as OOriginQuality, Magnitude as OMagnitude, \
    Comment as OComment, Catalog as OCatalog
import numpy as np
import scipy
from scipy.spatial import cKDTree
from mpi4py import MPI
from collections import defaultdict 
from pykml import parser
import uuid, copy
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees

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
        self.arrival_list = []
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
    # end func
# end class

class FDSNInv:
    def rtp2xyz(self, r, theta, phi):
        if(type(r)==float):
            xout = [0,0,0]
            rst = r * np.sin(theta);
            xout[0] = rst * np.cos(phi)
            xout[1] = rst * np.sin(phi)
            xout[2] = r * np.cos(theta)            
        else:
            xout = np.zeros((r.shape[0], 3))
            rst = r * np.sin(theta);
            xout[:, 0] = rst * np.cos(phi)
            xout[:, 1] = rst * np.sin(phi)
            xout[:, 2] = r * np.cos(theta)
        # end if
        return xout
    # end func        
    
    def __init__(self, fn, host=None, port=None):        
        def tree():
            def the_tree():
                return defaultdict(the_tree)
            # end func
            return the_tree()
        # end func
        
        if(host and port):
            self.host = host
            self.port = port
            try:
                self.client = Client('http://'+self.host+':%d'%self.port)

                if not self.client: raise Exception('Connection failed..')
            except Exception as e:
                print e
                print 'Failed to connect to client. Aborting..'
            # end try    
            self.inv = self.client.get_stations(level='station')
        else:
            self.fn = fn
            self.inv = read_inventory(fn)
        # end if

        self.t = tree()
        for n in self.inv.networks:
            for s in n.stations:
                self.t[n.code][s.code] = [s.longitude, s.latitude, s.elevation]        
            # end for
        # end for
        
        self.nsList = []
        self.nsCoordsList = []
        for nk, n in self.t.iteritems():
            for sk, s in n.iteritems():
                self.nsList.append('%s.%s'%(nk, sk))
                self.nsCoordsList.append(s)
            #end for
        # end for
        self.nsList = np.array(self.nsList)
        self.nsCoordsList = np.array(self.nsCoordsList)   
        
        ncoords = self.nsCoordsList.shape[0]
        self.xyz = self.rtp2xyz(6371e3*np.ones(ncoords),
                                np.radians(90-self.nsCoordsList[:,1]),
                                np.radians(self.nsCoordsList[:,0]))
                                
        self.kdtree = cKDTree(self.xyz)
    # end func
    
    def getClosestStations(self, lon, lat, maxdist=1e3):
        xyz = self.rtp2xyz(6371e3, 
                           np.radians(90-lat),
                           np.radians(lon))
        
        indices = self.kdtree.query_ball_point(xyz, maxdist)
        if(len(indices)):
            return list(self.nsList[indices]), self.nsCoordsList[indices]
        else:
            return None
        # end func
    # end func
# end class    

class Catalog():
    def __init__(self, isc_coords_file, fdsn_inventory, our_picks, event_folder, output_path):
        
        self.event_folder = event_folder
        self.output_path = output_path
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()
        self.counter = 0
        
        self.isc_coords_file = isc_coords_file
        k = parser.fromstring(open(self.isc_coords_file).read())

        ns = k.nsmap.values()[0]
        elist=k.findall('.//{%s}Placemark'%(ns))
        self.isc_coords_dict = defaultdict(list)

        for e in elist:
            self.isc_coords_dict[e['name']] = map(float, (str(e['Point']['coordinates'])).split(','))
        # end for


        # retrieve list of all event xml files
        self.event_folder = event_folder

        # retrieve list of all csv files
        self.csv_files = sorted(recursive_glob(self.event_folder, '*.csv'))

        self.fdsn_inventory = fdsn_inventory
        self.our_picks = our_picks
        self._load_events()
    # end func

    def _load_events_helper(self):
        eventList = []
        poTimestamps = []
        lines = []

        for ifn, fn in enumerate(self.csv_files):
            print 'Reading %s' % (fn)
            for iline, line in enumerate(open(fn, 'r').readlines()):
                lines.append(line)
            # end for
        # end for

        if(self.rank==0):
            for iline, line in enumerate(lines):
                evtline = ''
                if(line[0]=='#'):
                    items = line.split(',')
                    vals = map(float, items[1:])

                    year = int(vals[0])
                    month = int(vals[1])
                    day = int(vals[2])
                    hour = int(vals[3])
                    minute = int(vals[4])
                    second = vals[5]

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
                    event.public_id = int(eventid)
                    event.preferred_origin = origin
                    event.preferred_magnitude = Magnitude(mag, magtype)

                    eventList.append(event)
                    poTimestamps.append(origin.utctime.timestamp)
                else:
                    eventList[-1].preferred_origin.arrival_list.append(iline)
                # end if

                #if(iline%1000==0): print iline
            # end for
        # end if

        eventList = split_list(eventList, self.nproc)
        poTimestamps = split_list(poTimestamps, self.nproc)
        # broadcast workload to all procs
        eventList = self.comm.bcast(eventList, root=0)
        poTimestamps = self.comm.bcast(poTimestamps, root=0)

        print 'Processing %d events on rank %d'%(len(eventList[self.rank]), self.rank)
        for e in eventList[self.rank]:
            lineIndices = copy.copy(e.preferred_origin.arrival_list)
            e.preferred_origin.arrival_list = []
            for lineIndex in lineIndices:
                items = lines[lineIndex].split(',')
                vals = map(float, items[8:])

                year = int(vals[0])
                month = int(vals[1])
                day = int(vals[2])
                hour = int(vals[3])
                minute = int(vals[4])
                second = vals[5]

                utctime = None
                try:
                    utctime = UTCDateTime(year, month, day,
                                          hour, minute, second)
                except:
                    continue
                # end try

                try: lon = float(items[4])
                except: lon = 0
                try: lat = float(items[5])
                except: lat = 0
                try: elev = float(items[6])
                except: elev = 0
                
                distance = vals[-1]
                a = Arrival(items[3].strip(), items[0].strip(), items[2].strip(), items[1].strip(),
                            lon, lat, elev,
                            items[7].strip(), utctime, distance)
                e.preferred_origin.arrival_list.append(a)
            # end for
        # end for
        
        self.eventList = eventList[self.rank]
        self.poTimestamps = np.array(poTimestamps[self.rank])
    # end func
    
    def get_id(self):
        self.counter += 1
        return str(self.counter) + 'r%d'%self.rank
    # end func

    def _load_events(self):
        self._load_events_helper()
        cache = {}
        notFound = defaultdict(int)     
        oEvents = []
        missingStations = defaultdict(int)
        for e in self.eventList:
            if(e.preferred_origin and len(e.preferred_origin.arrival_list)):
                cullList = []
                for a in e.preferred_origin.arrival_list:                   
                    if(len(a.net)): continue
                    
                    seedid = '%s.%s.%s.%s'%(a.net, a.sta, a.loc, a.cha)
                    newCode = None
                    if(seedid not in cache):
                        sc = a.sta
                        lonlat = self.isc_coords_dict[sc]
                        if(len(lonlat)==0): 
                            cullList.append(a)
                            continue
                        # end if

                        r = self.fdsn_inventory.getClosestStations(lonlat[0], lonlat[1], maxdist=1e3)
                        #if(a.sta=='KUM'): print a.net, a.sta, a.loc, a.cha, r
                        if(not r):
                            notFound[sc]+=1
                        else:
                            for cr in r[0]:
                                c = cr.split('.')[0]
                                newCode = c
                            # end for
                        # end if

                        if(newCode):
                            cache[seedid] = newCode
                        # end if
                    else:
                        newCode = cache[seedid]
                    # end if

                    if(newCode):
                        #print a.net, newCode
                        a.net = newCode

                        sc = self.fdsn_inventory.t[a.net][a.sta]
                        if(type(sc)==defaultdict): 
                            cullList.append(a)
                            continue
                        # end if
                        da = gps2dist_azimuth(e.preferred_origin.lat, 
                                              e.preferred_origin.lon, 
                                              sc[1], sc[0])
                        dist = kilometers2degrees(da[0]/1e3)
                        if(np.fabs(a.distance-dist)>0.5): 
                            cullList.append(a)
                        # end if
                    # end if
                # end for
                for c in cullList: e.preferred_origin.arrival_list.remove(c)
            # end if
            
            # Create obspy event object         
            ci = OCreationInfo(author='GA', creation_time=UTCDateTime(),
                               agency_id='GA-iteration-1')   
            oid = self.get_id()
            origin = OOrigin(resource_id=OResourceIdentifier(id=oid),
                             time=UTCDateTime(e.preferred_origin.utctime),
                             longitude=e.preferred_origin.lon,
                             latitude=e.preferred_origin.lat,
                             depth=e.preferred_origin.depthkm*1e3,
                             method_id=OResourceIdentifier(id='unknown'),
                             earth_model_id=OResourceIdentifier(id='iasp91'),
                             evaluation_mode='automatic',
                             creation_info=ci)
            magnitude = OMagnitude(resource_id=OResourceIdentifier(id=self.get_id()),
                                   mag=e.preferred_magnitude.magnitude_value,
                                   magnitude_type=e.preferred_magnitude.magnitude_type,
                                   origin_id=OResourceIdentifier(id=oid),
                                   creation_info=ci)
            event = OEvent(resource_id=OResourceIdentifier(id=self.get_id()),
                           creation_info=ci,
                           event_type='earthquake')
            event.origins = [origin]
            event.magnitudes = [magnitude]
            event.preferred_magnitude_id = magnitude.resource_id
            event.preferred_origin_id = origin.resource_id
            
            # Insert old picks
            for a in e.preferred_origin.arrival_list:
                if(type(self.fdsn_inventory.t[a.net][a.sta]) == defaultdict):
                    missingStations[a.net+'.'+a.sta] += 1
                    continue
                # end if
                oldPick  = OPick( resource_id=OResourceIdentifier(id=self.get_id()),
                                  time=UTCDateTime(a.utctime),
                                  waveform_id=OWaveformStreamID(network_code=a.net,
                                                                station_code=a.sta,
                                                                channel_code=a.cha),
                                 methodID=OResourceIdentifier('unknown'),
                                 phase_hint=a.phase,
                                 evaluation_mode='automatic',
                                 creation_info=ci )

                oldArr = OArrival( resource_id=OResourceIdentifier(id=oldPick.resource_id.id+"#"),
                                   pick_id=oldPick.resource_id,
                                   phase=oldPick.phase_hint,
                                   distance=a.distance,
                                   earth_model_id=OResourceIdentifier('quakeml:ga.gov.au/earthmodel/iasp91'),
                                   creation_info=ci )
                
                event.picks.append(oldPick)
                event.preferred_origin().arrivals.append(oldArr)
            # end for

            # Insert our picks
            opList = self.our_picks.picks[e.public_id]
            if(len(opList)):
                for op in opList:
                    if(type(self.fdsn_inventory.t[op[1]][op[2]]) == defaultdict):
                        missingStations[op[1]+'.'+op[2]] += 1
                        continue
                    # end if
                    newPick  = OPick( resource_id=OResourceIdentifier(id=self.get_id()),
                                      time=UTCDateTime(op[0]),
                                      waveform_id=OWaveformStreamID(network_code=op[1],
                                                                    station_code=op[2],
                                                                    channel_code=op[3]),
                                     methodID=OResourceIdentifier('phasepapy/aicd'),
                                     backazimuth=op[-1],
                                     phase_hint=op[4],
                                     evaluation_mode='automatic',
                                     comments=op[6],
                                     creation_info=ci )

                    newArr = OArrival( resource_id=OResourceIdentifier(id=newPick.resource_id.id+"#"),
                                       pick_id=newPick.resource_id,
                                       phase=newPick.phase_hint,
                                       azimuth=op[-2],
                                       distance=op[-3],
                                       time_residual=op[5],
                                       time_weight=1.,
                                       earth_model_id=OResourceIdentifier('quakeml:ga.gov.au/earthmodel/iasp91'),
                                       creation_info=ci )
                    event.picks.append(newPick)
                    event.preferred_origin().arrivals.append(newArr)
                # end for
            # end if

            quality= OOriginQuality(associated_phase_count= len(e.preferred_origin.arrival_list) + 
                                                             len(self.our_picks.picks[e.public_id]),
                                    used_phase_count=len(e.preferred_origin.arrival_list) + 
                                                     len(self.our_picks.picks[e.public_id]))            
            event.preferred_origin().quality = quality
            oEvents.append(event)
        # end for // loop over e

        #print notFound
        print self.rank, missingStations

        cat = OCatalog(events=oEvents)        
        ofn = self.output_path + '/%d.xml'%(self.rank)
        cat.write(ofn, format='SC3ML')
    # end func
# end class

class OurPicks:
    def __init__(self, fnList, phaseList):
        self.fnList = fnList
        self.picks = defaultdict(list)
        
        for fn, phase in zip(self.fnList, phaseList):
            for iline, line in enumerate(open(fn, 'r')):
                if(iline == 0): continue

                items = line.split()
                """
                0 eventID
                1 originTimestamp
                2 mag
                3 originLon
                4 originLat
                5 originDepthKm
                6 net
                7 sta
                8 cha
                9 pickTimestamp
                10 stationLon
                11 stationLat
                12 az
                13 baz
                14 distance
                15 ttResidual
                16 snr
                17 qualityMeasureCWT
                18 domFreq
                19 qualityMeasureSlope
                20 bandIndex
                21 nSigma
                """

                t = float(items[9])
                d = float(items[14])
                baz = float(items[13])
                az = float(items[12])
                snr  = float(items[16])
                quality_measure_cwt = float(items[17])
                dom_freq = float(items[18])
                quality_measure_slope = float(items[19])
                bi  = int(items[20])
                nsigma = int(items[21])

                mag  = float(items[2])
                elon = float(items[3])
                elat = float(items[4])
                ed = float(items[5])
                slon = float(items[10])
                slat = float(items[11])
                residual = float(items[15])
                net = items[6]
                sta = items[7]
                cha = items[8]
                
                auxData = []
                auxData.append(OComment(text = 'phasepapy_snr = ' + str(snr) +
                                               ', quality_measure_cwt = ' + str(quality_measure_cwt) +
                                               ', dom_freq = ' + str(dom_freq) +
                                               ', quality_measure_slope = ' + str(quality_measure_slope) +
                                               ', band_index = ' + str(bi) +
                                               ', nsigma = ' + str(nsigma),
                                        force_resource_id=False))
                self.picks[int(float(items[0]))].append([t, net, sta, cha, phase, residual, auxData, elat, elon, slat, slon,
                                             d, az, baz])
            # end for
            print 'Read %s'%fn
        # end for
    # end func
# end class

if __name__=="__main__":
    PATH = '/g/data/ha3/Passive/Events/Unified/'
    INV = '/g/data1a/ha3/rakib/seismic/pst/tests/output/sc3Inventory.xml'
    ISC_COORDS_FILE = '/g/data1a/ha3/rakib/seismic/pst/tests/output/stations.kml'
    OFPATH = '/g/data1a/ha3/rakib/seismic/pst/tests/output/updatedEvents/'

    fi = FDSNInv(INV)
    if(0):
        op = OurPicks(['/g/data/ha3/rakib/seismic/pst/tests/results/p_arrivals_mag_4_and_above.txt',
                       '/g/data/ha3/rakib/seismic/pst/tests/results/s_arrivals_mag_4_and_above.txt',
                       '/g/data/ha3/rakib/seismic/pst/tests/results/temp/p_arrivals.txt',
                       '/g/data/ha3/rakib/seismic/pst/tests/results/temp/s_arrivals.txt'],
                       ['P', 'S', 'P', 'S'])
    else:
        op = OurPicks(['/g/data/ha3/rakib/seismic/pst/tests/output/pands_all_take12/p_combined.txt',
                       '/g/data/ha3/rakib/seismic/pst/tests/output/pands_all_take12/s_combined.txt'],
                       ['P', 'S'])
    # end if
    c = Catalog(ISC_COORDS_FILE, fi, op, PATH, OFPATH)
# end if

