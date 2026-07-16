from obspy.core import UTCDateTime
import numpy as np
from rtree import index
from glob import glob
from collections import defaultdict
from obspy import read
from obspy.core import Stream, Trace
import os
from tqdm import tqdm
from ordered_set import OrderedSet as set
from obspy import Inventory
import obspy
import copy

MAX_DATE = UTCDateTime(4102444800.0) #2100-01-01
MIN_DATE = UTCDateTime(-2208988800.0) #1900-01-01

UVW2ENZ = (1 / np.sqrt(6)) * np.array([
    [2.0, -1.0, -1.0],
    [0.0,  np.sqrt(3), -np.sqrt(3)],
    [np.sqrt(2), np.sqrt(2), np.sqrt(2)],
])

def galperin_to_enz(u, v, w):
    """
    Transform Galperin U, V, W components to orthonormal X, Y, Z.

    Parameters: u, v, w : array-like Galperin components.
    Returns
    -------
    e, n, z : ndarray Orthogonal Cartesian components.
    """
    ts_uvw = np.array([u, v, w])

    e, n, z = UVW2ENZ @ ts_uvw

    return e, n, z
# end func

def get_file_signature(file_path):
    """Return a dictionary of key file attributes to track changes, including absolute path."""
    abs_path = os.path.abspath(file_path)

    if not os.path.exists(abs_path):
        raise FileNotFoundError("File not found: {}".format(abs_path))
    # end if

    stat = os.stat(abs_path)

    return {
        'abs_path': abs_path,
        'st_size': stat.st_size,
        'st_mtime': stat.st_mtime,
        'st_ctime': stat.st_ctime,
        'st_ino': stat.st_ino,
        'st_dev': stat.st_dev,
    }
# end func

def cleanse_inventory(iinv: Inventory) -> Inventory:
    oinv = iinv.copy()

    for net in oinv.networks:
        net.comments = []
        for sta in net.stations:
            sta.comments = []
            for cha in sta.channels:
                cha.comments = []
            # end for
        # end for
    # end for

    return oinv
# end func

class InventoryAggregator:
    def __init__(self):
        tree = lambda: defaultdict(tree)
        self.net_dict = tree()
        self.sta_dict = tree()
        self.cha_dict = tree()
    # end func

    def append(self, inv: Inventory):
        for net in inv.networks:
            nc = net.code

            if (type(self.net_dict[nc]) == defaultdict):
                onet = copy.deepcopy(net)
                onet.stations = []
                self.net_dict[nc] = onet
            else:
                # update start/end dates for existing networks
                onet = self.net_dict[nc]

                if(onet.start_date and net.start_date):
                    if (onet.start_date > net.start_date): onet.start_date = net.start_date
                elif(onet.start_date is None and net.start_date):
                    onet.start_date = net.start_date
                # end if

                if(onet.end_date and net.end_date):
                    if(onet.end_date < net.end_date): onet.end_date = net.end_date
                # end if
            # end if

            for sta in net.stations:
                sc = sta.code

                if (type(self.sta_dict[nc][sc]) == defaultdict):
                    osta = copy.deepcopy(sta)
                    osta.channels = []
                    self.sta_dict[nc][sc] = osta
                else:
                    # update start/end dates for existing stations
                    osta = self.sta_dict[nc][sc]

                    if(osta.start_date and sta.start_date):
                        if (osta.start_date > sta.start_date): osta.start_date = sta.start_date
                    elif(osta.start_date is None and sta.start_date):
                        osta.start_date = sta.start_date
                    # end if

                    if(osta.end_date and sta.end_date):
                        if(osta.end_date < sta.end_date): osta.end_date = sta.end_date
                    # end if
                # end if

                for cha in sta.channels:
                    cc = cha.code
                    lc = cha.location_code

                    # set responses to None
                    try:
                        cc.response = None
                    except:
                        pass

                    if (type(self.cha_dict[nc][sc][lc][cc]) == defaultdict):
                        cha_copy = copy.deepcopy(cha)

                        """
                        Channel start- and end-dates in the inventory do not reflect actual waveform
                        data holdings. We therefore set the start- and end-times to None, so usable 
                        data is not lost e.g. when rotating waveform data for which corresponding 
                        metadata for the correct timeframes do not exist.
                        """
                        cha_copy.start_date = None
                        cha_copy.end_date = None
                        self.cha_dict[nc][sc][lc][cc] = cha_copy
                    # end if
                # end for
            # end for
        # end for
    # end func

    def summarize(self):
        oinv = Inventory(networks=[],
                         source=obspy.core.util.version.read_release_version())

        for nc in self.net_dict.keys():
            net = self.net_dict[nc]

            for sc in self.sta_dict[nc].keys():
                sta = self.sta_dict[nc][sc]

                for lc in self.cha_dict[nc][sc].keys():
                    for cc in self.cha_dict[nc][sc][lc].keys():
                        cha = self.cha_dict[nc][sc][lc][cc]

                        sta.channels.append(cha)
                    # end for
                # end for
                net.stations.append(sta)
            # end for

            # ensure network has a failsafe start-date
            if(net.start_date is None): net.start_date = UTCDateTime('1970-01-01')
            oinv.networks.append(net)
        # end for

        return oinv
    # end func
# end class

class MseedIndex:
    class StreamCache:
        def __init__(self):
            self.streams = defaultdict(list)
            self.read_times = defaultdict(list)
        # end func

        def get(self, fn):
            if (fn in self.streams.keys()):
                # print('found stream..')
                return self.streams[fn]
            else:
                # print('reading stream..')
                result = self._add(fn)
                self._cleanup()

                return result
            # end if
        # end func

        def flush(self):
            self.streams = defaultdict(list)
            self.read_times = defaultdict(list)
        # end func

        def _cleanup(self):
            MAX_STREAMS = 5
            # before = len(self.streams)
            while (len(self.streams) > MAX_STREAMS):
                time_key = sorted(self.read_times.keys())[0]
                file_key = self.read_times[time_key]

                self.read_times.pop(time_key)
                self.streams.pop(file_key)
            # wend
            # after = len(self.streams)
            # if(before > after): print('cleaned up {} streams..'.format(before-after))
        # end func

        def _add(self, fn):
            try:
                self.streams[fn] = read(fn)
                self.read_times[UTCDateTime.now().timestamp] = fn
                return self.streams[fn]
            except Exception as e:
                print("Failed to read {} with error {}. Moving along..".format(fn, e))
            # end try
        # end func
    # end class

    def __init__(self, mseed_folder, pattern):
        self.mseed_folder = mseed_folder
        self.tree = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
        self.stream_cache = MseedIndex.StreamCache()
        self.mseed_files = np.array(sorted(glob(os.path.join(self.mseed_folder, pattern))))
        self.coverage_dict = defaultdict(float) # dict keyed by nslc, with coverage as values

        fc = len(self.mseed_files)
        if(fc > 0):
            print('Found {} files:'.format(fc))
            print(os.path.basename(self.mseed_files[0]))
            print('..')
            print('..')
        else:
            raise RuntimeError('No mseed files found with pattern {}. Aborting..'.format(pattern))
        # end if

        print('Reading metadata from mseed files..')
        self.meta_list = []
        for i, mseed_file in enumerate(tqdm(self.mseed_files)):
            st = None
            try:
                st = read(mseed_file, headonly=True)
            except Exception as e:
                print("Failed to read {} with error {}. Moving along..".format(mseed_file, e))
                continue
            # end try

            for tr in st:
                nc, sc, lc, cc, st, et = \
                    tr.stats.network, tr.stats.station, tr.stats.location, \
                        tr.stats.channel, tr.stats.starttime.timestamp, \
                        tr.stats.endtime.timestamp

                # skip bogus traces
                if(nc == sc == lc == cc == ''): continue
                self.meta_list.append([i, nc, sc, lc, cc, st, et])

                nslc = '.'.join((nc, sc, lc, cc))
                # store coverage for unique channels and sampling rates
                cov = float(et - st) # coverage in seconds
                self.coverage_dict[nslc] += cov
            # end for
            # if (i > 0): break
        # end for

        print('\nCreating metadata index for a total of {} traces found..'.format(len(self.meta_list)))

        for row in tqdm(self.meta_list):
            idx, nc, sc, lc, cc, st, et = row

            if (type(self.tree[nc][sc][lc][cc]) != index.Index):
                self.tree[nc][sc][lc][cc] = index.Index()
            # end if
            self.tree[nc][sc][lc][cc].insert(idx, (st, 1, et, 1))
        # end for
    # end func

    def get_channel_coverages(self) -> defaultdict:
        return self.coverage_dict
    # end func

    def __getstate__(self):
        #print('pickling..')
        return self.__dict__
    # end func

    def __setstate__(self, d):
        #print('unpickling..')
        self.__dict__ = d

        # recreate tree, because serialization/deserialization across
        # processes do not preserve hashed objects
        self.tree = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
        for row in self.meta_list:
            idx, nc, sc, lc, cc, st, et = row

            if (type(self.tree[nc][sc][lc][cc]) != index.Index):
                self.tree[nc][sc][lc][cc] = index.Index()
            # end if
            self.tree[nc][sc][lc][cc].insert(idx, (st, 1, et, 1))
        # end for
    # end func

    def flush_cache(self):
        self.stream_cache.flush()
    # end func

    def get_waveforms(self, net, sta, loc, cha, st: UTCDateTime, et: UTCDateTime):
        epsilon = 1e-5
        st_ts = st.timestamp + epsilon
        et_ts = et.timestamp - epsilon

        result = Stream([])
        try:
            target_index = self.tree[net][sta][loc][cha]

            if (type(target_index) == index.Index):
                file_indices = np.array(list(target_index.intersection((st_ts, 1, et_ts, 1))), dtype='i4')

                # since file names are repeated for multiple traces, we need a unique set
                for mfile in set(self.mseed_files[file_indices]):
                    temp_stream = self.stream_cache.get(mfile).select(network=net,
                                                                      station=sta,
                                                                      location=loc,
                                                                      channel=cha)
                    result += temp_stream.slice(st, et, nearest_sample=False).copy()
                # end for
            else:
                print('empty index')
            # end if
        except Exception as e:
            print('error in mseedindex', str(e))
        # end try

        return result
    # end func

    def get_stations(self, st: UTCDateTime, et: UTCDateTime, net=None, sta=None, loc=None, cha=None):
        epsilon = 1e-5
        st_ts = st.timestamp + epsilon
        et_ts = et.timestamp - epsilon

        _net = _sta = _loc = _cha = None

        if (net == None):
            _net = self.tree.keys()
        else:
            _net = [net]

        result = []
        for nc in _net:
            if (sta == None):
                _sta = self.tree[nc].keys()
            else:
                _sta = [sta]

            for sc in _sta:
                if (loc == None):
                    _loc = self.tree[nc][sc].keys()
                else:
                    _loc = [loc]

                for lc in _loc:
                    if (cha == None):
                        _cha = self.tree[nc][sc][lc].keys()
                    else:
                        _cha = [cha]

                    for cc in _cha:
                        target_index = self.tree[nc][sc][lc][cc]
                        if (type(target_index) == index.Index):
                            entries = list(target_index.intersection((st_ts, 1, et_ts, 1)))

                            if (len(entries)): result.append((nc, sc, lc, cc))
                        # end if
                    # end for
                # end for
            # end for
        # end for

        return result
    # end func

    def get_time_range(self, net, sta, loc, cha):

        target_index = self.tree[net][sta][loc][cha]

        if (type(target_index) == index.Index):
            bounds = target_index.bounds
            return UTCDateTime(bounds[0]), UTCDateTime(bounds[2])
        # end if

        return MAX_DATE, MIN_DATE
    # end func
# end func

if __name__=="__main__":
    msi = MseedIndex('/g/data/ha3/ac5759/semi-perm-iris', '*AXCOZ*HHZ*j168.mseed')

    print(msi.tree['AU'].keys())
    print(msi.tree['AU']['AXCOZ'].keys())
    r = msi.get_waveforms('AU', 'AXCOZ', '00', 'HHZ', UTCDateTime(2022, 1, 3, 22, 9, 26), UTCDateTime(2023, 1, 4, 0, 0))
    print(r)
    print(msi.get_time_range('AU', 'AXCOZ', '00', 'HHZ'))
# end if
