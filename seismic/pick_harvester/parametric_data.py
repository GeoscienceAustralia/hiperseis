"""
Description:
    Reads parametric data from a catalog and optionally provided automatic picks

References:

CreationDate:   10/06/22
Developer:      rakib.hassan@ga.gov.au

Revision History:
    LastUpdate:     10/06/22   RH
    LastUpdate:     dd/mm/yyyy  Who     Optional description
"""

from ordered_set import OrderedSet as set
import numpy as np
from obspy import UTCDateTime
from collections import defaultdict
from mpi4py import MPI
import uuid
import os, shutil
from itertools import combinations
from pyproj import Geod
from seismic.pick_harvester.utils import split_list
import h5py

class ParametricData:
    def __init__(self, csv_catalog, auto_pick_files=[], auto_pick_phases=[],
                 events_only=False, phase_list='P Pg Pb Pn S Sg Sb Sn', temp_dir='./'):
        """
        :param csv_catalog: path to catalog file in csv format
        :param auto_pick_files: optional list of files containing automatic arrivals. Each file
                                should contain arrivals of the same phase (e.g. P and S arrivals
                                output from pick.py)
        :param auto_pick_phases: optional list of phases corresponding to each file provided in
                                 auto_pick_files
        :param events_only: whether only events should be loaded
        :param phase_list: a space-separated list of phases to be read from the catalogue -- all
                           other phases not in this list are discarded
        :param temp_dir: path to a temporary folder to be used for syncing data across processors.
                         Note that this temporary folder must be accessible by all MPI ranks, e.g.
                         within a project folder on the NCI.
        """
        self.EARTH_RADIUS_KM = 6371.
        self.DEG2KM = np.pi/180 * self.EARTH_RADIUS_KM
        self.STATION_DIST_M = 1e3 # distance in metres within which neighbouring stations are coalesced
        self.csv_catalog = csv_catalog
        self.auto_pick_files = auto_pick_files
        self.auto_pick_phases = auto_pick_phases
        self.events_only = events_only
        self.phase_list = set(map(str.strip, phase_list.split())) if len(phase_list) else set()
        self.p_phases = np.array([item for item in self.phase_list if item[0]=='P'])
        self.s_phases = np.array([item for item in self.phase_list if item[0]=='S'])
        self._temp_dir = temp_dir
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()
        self._geod = Geod(a=180/np.pi, f=0)
        self.events = []
        self.event_id_to_idx = None
        self.arrivals = []
        self.local_events_indices = None
        self.local_arrivals_indices = None

        # create temp_dir
        self._make_temp_dir()

        # sanity check
        for item in self.phase_list:
            if (',' in item or not len(item)): raise ValueError('Invalid phase {} in phase_list'.format(item))
        # end for

        self.event_fields = {'names': ['source', 'event_id', 'origin_ts', 'mag', 'lon', 'lat', 'depth_km'],
                             'formats': ['S10', 'i4', 'f8', 'f4', 'f4', 'f4', 'f4']}
        self.arrival_fields = {
            'names': ['event_id', 'net', 'sta', 'loc', 'cha', 'lon', 'lat', 'elev_m', 'phase', 'arrival_ts',
                      'quality_measure_slope'],
            'formats': ['i4', 'S10', 'S10', 'S10', 'S10', 'f4', 'f4', 'f4', 'S10', 'f8', 'f4']}

        # load events and arrivals
        if(1):
            self.events, self.arrivals = self._load_catalog()
            if(self.rank == 0):
                np.save('events.npy', self.events)
                np.save('arrivals.npy', self.arrivals)
            # end if
        else:
            self.events = np.load('events.npy')
            self.arrivals= np.load('arrivals.npy')
        # end if

        np.random.seed(0)
        np.random.shuffle(self.arrivals)
        np.random.shuffle(self.events)

        # create a map to translate event-id to array index
        self.event_id_to_idx = np.ones(np.max(self.events['event_id']) + 1, dtype='i4') * -1
        for i in np.arange(len(self.events)): self.event_id_to_idx[self.events['event_id'][i]] = i

        if(self.rank == 0): print('Loaded catalogue with {} events and {} arrivals..'.format(len(self.events),
                                                                                             len(self.arrivals)))
        # coalesce network codes
        self.local_events_indices = np.array(split_list(np.arange(len(self.events)), self.nproc)[self.rank], dtype='i4')
        self.local_arrivals_indices = np.array(split_list(np.arange(len(self.arrivals)), self.nproc)[self.rank], dtype='i4')
        self._coalesce_network_codes()

        # print(self.arrivals)
    # end func

    def has_arrival(self, event_id, network, station, station_lon, station_lat, phase_list=''):
        if (',' in phase_list): raise ValueError('A space-separated list of phases expected')
        if (len(phase_list)): phase_list = set(map(str.strip, phase_list.split()))

        event_imask = self.arrivals['event_id'] == event_id
        m_nets = self.arrivals[event_imask]['net']
        m_stas = self.arrivals[event_imask]['sta']
        m_phases = self.arrivals[event_imask]['phase']
        m_lons = self.arrivals[event_imask]['lon']
        m_lats = self.arrivals[event_imask]['lat']

        # print(m_nets, m_stas, m_phases)
        result = False
        for net, sta, phase, lon, lat in zip(m_nets, m_stas, m_phases, m_lons, m_lats):
            net, sta, phase = net.decode(), sta.decode(), phase.decode()
            if (network == net and station == sta):
                if (len(phase_list)):
                    if (phase in phase_list):
                        result = True
                        break
                    # end if
                else:
                    result = True
                    break
                # end if
            elif (station == sta):
                # use station coordinates to account for cases when network codes
                # may have been mangled in _coalesce_network_codes
                _, _, dist = self._geod.inv(lon, lat, station_lon, station_lat)
                if(dist < self.STATION_DIST_M):
                    if (len(phase_list)):
                        if (phase in phase_list):
                            result = True

                            #print(net, sta, [lon, lat], network, station, [station_lon, station_lat])
                            break
                        # end if
                    else:
                        result = True
                        break
                    # end if
                # end if
            # end if
        # end for
        return result

    # end func

    def __del__(self):
        self.comm.Barrier()
        if (self.rank == 0):
            print('Removing temp-dir: {}..'.format(self._temp_dir))
            shutil.rmtree(self._temp_dir)
        # end if
    # end func

    def _coalesce_network_codes(self):
        if(self.rank == 0): print('Coalescing network codes..')

        if(0):
            self.arrivals['net'] = np.load('coalesced_net.npy')
            return
        # end if

        iter_count = 0
        while(1):
            dupdict = defaultdict(set)
            coordsdict = defaultdict(list)
            for arrival in self.arrivals:
                if(arrival['net'] not in dupdict[arrival['sta']]):
                    dupdict[arrival['sta']].add(arrival['net'])
                    coordsdict[arrival['net'] + b'.' + arrival['sta']] = [arrival['lon'], arrival['lat']]
                # end if
            # end for

            swap_map = defaultdict(str)
            swap_fail_map = defaultdict(list)
            for sta, nets in dupdict.items():
                if(len(nets)>1):
                    nets = sorted(list(nets), reverse=True)
                    for net1, net2 in combinations(nets, r=2):
                        lon1, lat1 = coordsdict[net1 + b'.' + sta]
                        lon2, lat2 = coordsdict[net2 + b'.' + sta]

                        _, _, dist = self._geod.inv(lon1, lat1, lon2, lat2)
                        dist *= self.DEG2KM * 1e3 #m
                        if(dist < self.STATION_DIST_M):
                            swap_map[(net1, sta)] = net2
                        else:
                            swap_fail_map[(net1, sta)] = [net2, dist]
                        # end if
                    # end for
                # end if
            # end for

            sum = 0
            local_net_codes = np.array(self.arrivals['net'][self.local_arrivals_indices])
            local_sta_codes = np.array(self.arrivals['sta'][self.local_arrivals_indices])
            for (net1, sta), net2 in swap_map.items():
                net_sta_match = (local_net_codes == net1) & (local_sta_codes == sta)
                sum += np.sum(net_sta_match)

                local_net_codes[net_sta_match] = net2
            # end for

            self.arrivals['net'] = self._sync_var(local_net_codes)

            iter_count += 1
            if(len(swap_map) == 0): break
        # wend
        if(self.rank==0): np.save('coalesced_net.npy', self.arrivals['net'])
    # end func

    def _make_temp_dir(self):
        if(self.rank == 0):
            self._temp_dir = os.path.join(self._temp_dir, str(uuid.uuid4()))
            os.makedirs(self._temp_dir, exist_ok=True)
        # end if
        self.comm.Barrier()
        self._temp_dir = self.comm.bcast(self._temp_dir, root=0)
    # end func

    def _load_catalog(self):
        events = []
        arrival_ids = []
        arrivals = []

        if (self.rank == 0):
            print('Loading events from {}'.format(self.csv_catalog))
            event_id = -1
            for iline, line in enumerate(open(self.csv_catalog, 'r')):
                if (line[0] == '#'):
                    """
                      0      1    2   3   4   5   6    7     8       9      10   11  12  13  14   15      16        17
                    source, YYYY, MM, DD, hh, mm, ss, elon, elat, edepthkm, nph, mb, ms, ml, mw, ev_id, azim_gap, int_id
                    """
                    items = line.split(',')

                    year, month, day, hour, minute = map(int, items[1:6])
                    second, lon, lat, depth_km = map(float, items[6:10])
                    mb, ms, ml, mw = map(float, items[11:15])

                    if (lon < -180 or lon > 180):
                        print('Invalid origin lon on line {}. Moving along..'.format(iline + 1))
                        continue
                    # end if
                    if (lat < -90 or lat > 90):
                        print('Invalid origin lat on line {}. Moving along..'.format(iline + 1))
                        continue
                    # end if
                    if (depth_km < 0):
                        print('Invalid origin depth on line {}. Moving along..'.format(iline + 1))
                        continue
                    # end if

                    mag = 0
                    if (mw > 0):
                        mag = mw
                    elif (ms > 0):
                        mag = ms
                    elif (mb > 0):
                        mag = mb
                    elif (ml > 0):
                        mag = ml

                    event_id = int(items[-1])
                    otime = None
                    try:
                        otime = UTCDateTime(year, month, day, hour, minute, second)
                    except:
                        raise ValueError('Invalid date parameters on line {}'.format(iline))

                    events.append((items[0][1:], event_id, np.float64(otime.timestamp), mag, lon, lat, depth_km))
                else:
                    if (self.events_only):
                        continue
                    else:
                        arrival_ids.append((iline, event_id))  # line number and event_id for arrivals
                # end if
            # end for

            # convert events to a structured array
            events = np.array(events, dtype=self.event_fields)
            np.save(os.path.join(self._temp_dir, 'events.npy'), events)

            if (not self.events_only):
                print('Loading {} arrivals from {}'.format(len(arrival_ids), self.csv_catalog))
                idx = 0
                for f_iline, line in enumerate(open(self.csv_catalog, 'r')):
                    if (idx >= len(arrival_ids)): break
                    iline, event_id = arrival_ids[idx]
                    if (f_iline != iline):
                        continue
                    else:
                        idx = idx + 1
                    # end if
                    """
                     0    1    2    3    4    5    6        7      8    9  10  11  12  13     14        15
                    sta, cha, loc, net, lon, lat, elev_m, phase, YYYY, MM, DD, hh, mm, ss, ang_dist, event_id
                    """
                    items = line.split(',')

                    sta, cha, loc, net = map(str.strip, items[0:4])
                    lon = lat = elev_m = 0
                    try:
                        lon, lat, elev_m = map(float, items[4:7])
                    except:
                        continue

                    phase = items[7].strip()
                    # skip unwanted arrivals
                    if (phase not in self.phase_list): continue

                    year, month, day, hour, minute = map(int, items[8:13])
                    second = float(items[13])

                    if (lon < -180 or lon > 180): raise ValueError('Invalid origin lon on line {}'.format(iline))
                    if (lat < -90 or lat > 90): raise ValueError('Invalid origin lat on line {}'.format(iline))

                    atime = None
                    try:
                        atime = UTCDateTime(year, month, day, hour, minute, second)
                    except:
                        raise ValueError('Invalid date parameters on line {}'.format(iline))

                    arrivals.append((event_id, net, sta, loc, cha, lon, lat, elev_m, phase,
                                     np.float64(atime.timestamp), -1))  # -1 for quality_measure_slope
                # end for

                # convert arrivals to a structured array
                arrivals = np.array(arrivals, dtype=self.arrival_fields)

                # load auto_pick files
                for auto_pick_fn, phase in zip(self.auto_pick_files, self.auto_pick_phases):
                    arrivals = np.hstack([arrivals, self._load_automatic_picks(auto_pick_fn, phase)])
                # end for
                assert np.all(arrivals['event_id'] >= 0), 'Invalid event-ids found in arrivals'
            # end if
            np.save(os.path.join(self._temp_dir, 'arrivals.npy'), arrivals)
        # end if
        self.comm.Barrier()

        if (self.rank > 0):
            events = np.load(os.path.join(self._temp_dir, 'events.npy'))
            if (not self.events_only): arrivals = np.load(os.path.join(self._temp_dir, 'arrivals.npy'))
        # end if

        return events, arrivals
    # end func

    def _load_automatic_picks(self, fn, phase):
        print('Loading automatic arrivals from {}'.format(fn))

        arrivals = []
        for iline, line in enumerate(open(fn, 'r')):
            if (iline == 0): continue
            """
            #eventID originTimestamp mag originLon originLat originDepthKm net sta cha pickTimestamp stationLon stationLat stationElev_m az baz distance ttResidual snr qualityMeasureCWT domFreq qualityMeasureSlope bandIndex nSigma
              0             1         2      3         4            5       6   7   8        9           10         11           12      13 14     15        16      17        18            19          20              21       22  
            """
            items = line.split()
            event_id = int(items[0])
            loc = ''
            net, sta, cha = map(str.strip, items[6:9])
            arrival_time = np.float64(items[9])
            lon, lat, elev_m = map(float, items[10:13])
            quality_measure_slope = float(items[20])

            arrivals.append(
                (event_id, net, sta, loc, cha, lon, lat, elev_m, phase, arrival_time, quality_measure_slope))
        # end for
        arrivals = np.array(arrivals, dtype=self.arrival_fields)

        return arrivals
    # end func

    def _sync_var(self, rank_values):
        # sync variable across ranks
        counts = np.array(self.comm.allgather(len(rank_values)), dtype='i4')
        nelem = np.sum(counts)
        dtype = rank_values.dtype
        displacements = np.zeros(self.nproc, dtype='i4')
        displacements[1:] = np.cumsum(counts[:-1])
        global_values = np.zeros(nelem, dtype=dtype)

        fn = self._temp_dir + '/sync.h5'

        for irank in np.arange(self.nproc):
            if(self.rank == irank):
                hf = h5py.File(fn, 'a')
                dset = hf.create_dataset("%d" % (self.rank), rank_values.shape, dtype=rank_values.dtype)
                dset[:] = rank_values
                hf.close()
            # end if
            self.comm.Barrier()
        # end for

        hf = h5py.File(fn, 'r')
        for irank in np.arange(self.nproc):
            b, e = displacements[irank], displacements[irank]+counts[irank]
            global_values[b:e] = hf['{}'.format(irank)][:]
        # end for
        hf.close()

        self.comm.Barrier()
        if(self.rank == 0): os.remove(fn)

        return global_values
    # end func
# end class

if __name__ == "__main__":
    if (0):
        pd = ParametricData('./small_merge_catalogues_output.csv',
                            auto_pick_files=['small_p_combined.txt', 'small_s_combined.txt'],
                            auto_pick_phases=['P', 'S'],
                            events_only=False)
    else:
        pd = ParametricData('./merge_catalogues_output.csv',
                            auto_pick_files=['old_p_combined.txt', 'old_s_combined.txt'],
                            auto_pick_phases=['P', 'S'],
                            events_only=False)
    # end if
# end if
