from ordered_set import OrderedSet as set
import numpy as np
from obspy import UTCDateTime
from collections import defaultdict
from seismic.pick_harvester.utils import split_list
from mpi4py import MPI
import uuid
import os, shutil


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
        self.event_id_to_idx = None
        self.arrivals = []
        self.source_enum = defaultdict(int)
        self.arrival_source = None

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
        self.events, self.arrivals = self._load_catalog()

        # create a map to translate event-id to array index
        self.event_id_to_idx = np.ones(np.max(self.events['event_id']) + 1, dtype='i4') * -1
        for i in np.arange(len(self.events)): self.event_id_to_idx[self.events['event_id'][i]] = i

        # label arrivals by event-source
        if (not self.events_only): self._label_arrivals()

        print(
            'Completed loading catalogue with {} events and {} arrivals..'.format(len(self.events), len(self.arrivals)))
        # print(self.arrivals)

    # end func

    def has_arrival(self, event_id, network, station, phase_list=''):
        if (',' in phase_list): raise ValueError('A space-separated list of phases expected')
        if (len(phase_list)): phase_list = set(map(str.strip, phase_list.split()))

        event_filter = self.arrivals['event_id'] == event_id
        m_nets = self.arrivals[event_filter]['net']
        m_stas = self.arrivals[event_filter]['sta']
        m_phases = self.arrivals[event_filter]['phase']

        # print(m_nets, m_stas, m_phases)
        result = False
        for net, sta, phase in zip(m_nets, m_stas, m_phases):
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
            # end if
        # end for
        return result

    # end func

    def _label_arrivals(self):
        self.arrival_source = np.zeros(len(self.arrivals), dtype='i4')

        sources = set(self.events['source'])
        for i, source in enumerate(sources): self.source_enum[source] = i + 1

        print('Labelling arrivals by event source {}..'.format(self.source_enum.items()))

        arrival_source = np.zeros(len(self.arrivals), dtype='i4')
        for source in sources:
            enum = self.source_enum[source]
            sids = np.argwhere(self.events['source'] == source).flatten()
            source_eids = set(list(self.events['event_id'][sids]))

            for i, eid in enumerate(self.arrivals['event_id']):
                if (eid in source_eids): arrival_source[i] = enum
            # end for
        # end for

        self.arrival_source = arrival_source

        assert np.all(self.arrival_source > 0), 'Arrivals found with no corresponding event-ids..'

        if (0):
            print('arrival_source {}'.format(self.arrival_source))
            for i in np.arange(len(self.arrivals)):
                print(self.arrivals[i])
            # end for
        # end if

    # end func

    def _load_catalog(self):
        events = []
        arrival_ids = []
        arrivals = []

        tempdir = None
        if (self.rank == 0):
            tempdir = os.path.join(os.path.dirname(self.csv_catalog), str(uuid.uuid4()))
            os.makedirs(tempdir, exist_ok=True)

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
            np.save(os.path.join(tempdir, 'events.npy'), events)

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
            np.save(os.path.join(tempdir, 'arrivals.npy'), arrivals)
        # end if
        self.comm.Barrier()

        tempdir = self.comm.bcast(tempdir, root=0)
        if (self.rank > 0):
            events = np.load(os.path.join(tempdir, 'events.npy'))
            if (not self.events_only): arrivals = np.load(os.path.join(tempdir, 'arrivals.npy'))
        # end if
        self.comm.Barrier()
        if (self.rank == 0): shutil.rmtree(tempdir)

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
