from mpi4py import MPI
from ordered_set import OrderedSet as set
import numpy as np
from obspy import UTCDateTime
from scipy.spatial import cKDTree
import pyproj
import gc

from seismic.pick_harvester.utils import split_list

class ParametricData:
    def __init__(self, csv_catalog, auto_pick_files=[], auto_pick_phases=[],
                 events_only=False, phase_list='P Pg Pb P* Pn S Sg Sb S* Sn'):
        self.csv_catalog = csv_catalog
        self.auto_pick_files = auto_pick_files
        self.auto_pick_phases = auto_pick_phases
        self.events_only = events_only
        self.phase_list = set(map(str.strip, phase_list.split())) if len(phase_list) else set()
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()
        self.events = []
        self.arrivals = []
        self.kdtree = None
        self._lonlatalt2xyz = None

        # sanity check
        for item in self.phase_list:
            if(',' in item or not len(item)): raise ValueError('Invalid phase {} in phase_list'.format(item))
        # end for

        self.event_fields = {'names': ['source', 'event_id', 'origin_ts', 'mag', 'lon', 'lat', 'depth_km'],
                             'formats': ['S10', 'i4', 'f8', 'f4', 'f4', 'f4', 'f4']}
        self.arrival_fields = {'names': ['event_id', 'net', 'sta', 'loc', 'cha', 'lon', 'lat', 'elev_m', 'phase', 'arrival_ts', 'quality_measure_slope'],
                               'formats': ['i4', 'S10', 'S10', 'S10', 'S10', 'f4', 'f4', 'f4', 'S10', 'f8', 'f4']}

        # load events and arrivals
        self.events, self.arrivals = self._load_catalog()

        # initialize kd-tree
        self._initialize_kdtree()

        print('rank {}: {} arrivals'.format(self.rank, len(self.arrivals)))
        #if(self.rank == 0): print(self.arrivals)
    # end func

    def has_arrival(self, event_id, network, station, phase_list=''):
        if(',' in phase_list): raise ValueError('A space-separated list of phases expected')
        if(len(phase_list)): phase_list = set(map(str.strip, phase_list.split()))

        event_filter = self.arrivals['event_id'] == event_id
        m_nets = self.arrivals[event_filter]['net']
        m_stas = self.arrivals[event_filter]['sta']
        m_phases = self.arrivals[event_filter]['phase']

        #print(m_nets, m_stas, m_phases)
        result = False
        for net, sta, phase in zip(m_nets, m_stas, m_phases):
            net, sta, phase = net.decode(), sta.decode(), phase.decode()
            if(network == net and station == sta):
                if(len(phase_list)):
                    if(phase in phase_list):
                        result = True
                        break
                    # end if
                else:
                    result = True
                    break
                # end if
            # end if
        # end for
        return result
    # end func

    def _initialize_kdtree(self, ellipsoidal_distance=True):
        ER = 6371e3 #km

        elons = self.events['lon']
        elats = self.events['lat']
        ealts = -self.events['depth_km'] * 1e3
        xyz = None
        if(ellipsoidal_distance):
            transformer = pyproj.Transformer.from_crs(
                            {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
                            {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'})
            self._lonlatalt2xyz = lambda lon, lat, alt: np.vstack(transformer.transform(lon, lat, alt,
                                                                                        radians=False)).T
        else:
            def rtp2xyz(r, theta, phi):
                xout = np.zeros((r.shape[0], 3))
                rst = r * np.sin(theta)
                xout[:, 0] = rst * np.cos(phi)
                xout[:, 1] = rst * np.sin(phi)
                xout[:, 2] = r * np.cos(theta)
                return xout
            # end func

            self._lonlatalt2xyz = lambda lon, lat, alt: rtp2xyz(np.atleast_1d(ER + alt),
                                                                np.atleast_1d(np.radians(90 - lat)),
                                                                np.atleast_1d(np.radians(lon)))
        # end if
        xyz = self._lonlatalt2xyz(elons, elats, ealts)
        self.kdtree = cKDTree(xyz)

        if(self.rank == 0):
            pass
            #print(xyz)
        # end if
    # end func

    def _gather_arrivals(self, rank_arrivals):
        # gather arrivals from all ranks on all ranks
        arrival_counts = np.array(self.comm.allgather(len(rank_arrivals)))

        nelem = np.sum(arrival_counts)
        displacements = np.zeros(self.nproc)
        displacements[1:] = np.cumsum(arrival_counts[:-1])
        global_arrivals = np.empty(nelem, dtype=self.arrival_fields)
        type_map = {'i4': MPI.INT, 'f4': MPI.FLOAT, 'f8': MPI.DOUBLE}
        for name, dtype in zip(self.arrival_fields['names'], self.arrival_fields['formats']):
            if(dtype in ['i4', 'f4', 'f8']):
                temp = np.zeros(nelem, dtype=dtype)
                self.comm.Allgatherv(np.array(rank_arrivals[name]),
                                     [temp, arrival_counts,
                                      displacements, type_map[dtype]])
                global_arrivals[name][:] = temp
            else:
                length = 10
                temp = np.empty(nelem*length, dtype='b')
                self.comm.Allgatherv(np.array(rank_arrivals[name]).tobytes(),
                                     [temp, arrival_counts*length,
                                      displacements*length, MPI.CHAR])
                global_arrivals[name] = np.frombuffer(temp, dtype='S10')
            # end if
        # end for
        return global_arrivals
    # end func

    def _load_catalog(self):
        events = []
        arrival_ids = []
        arrivals = []

        if(self.rank==0):
            print ('Loading events from {}'.format(self.csv_catalog))

            event_id = -1
            for iline, line in enumerate(open(self.csv_catalog, 'r')):
                if(line[0]=='#'):
                    """
                      0      1    2   3   4   5   6    7     8       9      10   11  12  13  14   15      16        17
                    source, YYYY, MM, DD, hh, mm, ss, elon, elat, edepthkm, nph, mb, ms, ml, mw, ev_id, azim_gap, int_id
                    """
                    items = line.split(',')

                    year, month, day, hour, minute = map(int, items[1:6])
                    second, lon, lat, depth_km = map(float, items[6:10])
                    mb, ms, ml, mw = map(float, items[11:15])

                    if(lon < -180 or lon > 180): raise ValueError('Invalid origin lon on line {}'.format(iline))
                    if(lat < -90 or lat > 90): raise ValueError('Invalid origin lat on line {}'.format(iline))
                    if(depth_km < 0): raise ValueError('Invalid origin depth on line {}'.format(iline))

                    mag = 0
                    if(mw>0): mag = mw
                    elif(ms>0): mag = ms
                    elif(mb>0): mag = mb
                    elif(ml>0): mag = ml

                    event_id = int(items[-1])
                    otime = None
                    try: otime = UTCDateTime(year, month, day, hour, minute, second)
                    except: raise ValueError('Invalid date parameters on line {}'.format(iline))

                    events.append((items[0][1:], event_id, np.float64(otime.timestamp), mag, lon, lat, depth_km))
                else:
                    if(self.events_only): continue
                    else: arrival_ids.append((iline, event_id)) # line number and event_id for arrivals
                # end if
            # end for

            # convert events to a structured array
            events = np.array(events, dtype=self.event_fields)

            arrival_ids = split_list(arrival_ids, self.nproc)
        # end if (load events on rank 0)

        # broadcast events to all ranks
        events = self.comm.bcast(events, root=0)

        if(not self.events_only):
            # process arrivals on all ranks
            arrival_ids = self.comm.scatter(arrival_ids, root=0)
            file_content = open(self.csv_catalog).readlines()

            print ('Loading {} arrivals from {} on rank {}'.format(len(arrival_ids), self.csv_catalog, self.rank))
            for iline, event_id in arrival_ids:
                """
                 0    1    2    3    4    5    6        7      8    9  10  11  12  13     14        15
                sta, cha, loc, net, lon, lat, elev_m, phase, YYYY, MM, DD, hh, mm, ss, ang_dist, event_id
                """
                items = file_content[iline].split(',')

                sta, cha, loc, net = map(str.strip, items[0:4])
                lon = lat = elev_m = 0
                try: lon, lat, elev_m = map(float, items[4:7])
                except: continue

                phase = items[7].strip()
                # skip unwanted arrivals
                if(phase not in self.phase_list): continue

                year, month, day, hour, minute = map(int, items[8:13])
                second = float(items[13])

                if(lon < -180 or lon > 180): raise ValueError('Invalid origin lon on line {}'.format(iline))
                if(lat < -90 or lat > 90): raise ValueError('Invalid origin lat on line {}'.format(iline))

                atime = None
                try: atime = UTCDateTime(year, month, day, hour, minute, second)
                except: raise ValueError('Invalid date parameters on line {}'.format(iline))

                arrivals.append((event_id, net, sta, loc, cha, lon, lat, elev_m, phase,
                                 np.float64(atime.timestamp), -1)) # -1 for quality_measure_slope
            # end for
            file_content = None # free up memory
            gc.collect()

            # convert arrivals to a structured array
            arrivals = np.array(arrivals, dtype=self.arrival_fields)

            # load auto_pick files
            for auto_pick_fn, phase in zip(self.auto_pick_files, self.auto_pick_phases):
                arrivals = np.hstack([arrivals, self._load_automatic_picks(auto_pick_fn, phase)])
            # end for

            arrivals = self._gather_arrivals(arrivals)
        # end if

        return events, arrivals
    # end func

    def _load_automatic_picks(self, fn, phase):
        arrival_ids = []
        if (self.rank == 0):
            for iline, line in enumerate(open(fn, 'r')):
                if(not iline): continue # skip header
                arrival_ids.append(iline)
            # end for
            arrival_ids = split_list(arrival_ids, self.nproc)
        # end if
        arrival_ids = self.comm.scatter(arrival_ids, root=0)

        print('Loading {} automatic arrivals from {} on rank {}'.format(len(arrival_ids), fn, self.rank))

        arrivals = []
        file_content = open(fn).readlines()
        for iline in arrival_ids:
            """
            #eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp stationLon stationLat stationElev_m az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma
              0             1         2      3         4            5       6   7   8        9           10         11           12      13 14     15        16      17        18            19          20              21       22  
            """
            items = file_content[iline].split()
            event_id = int(items[0])
            loc = ''
            net, sta, cha = map(str.strip, items[6:9])
            arrival_time = np.float64(items[9])
            lon, lat, elev_m = map(float, items[10:13])
            quality_measure_slope = float(items[20])

            arrivals.append((event_id, net, sta, loc, cha, lon, lat, elev_m, phase, arrival_time, quality_measure_slope))
        # end for
        arrivals = np.array(arrivals, dtype=self.arrival_fields)

        return arrivals
    # end func
# end class

