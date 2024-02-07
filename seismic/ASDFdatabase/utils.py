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
from seismic.misc import split_list
from obspy import Inventory

MAX_DATE = UTCDateTime(4102444800.0)
MIN_DATE = UTCDateTime(-2208988800.0)

def remove_comments(iinv: Inventory) -> Inventory:
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
            # end for
            # if (i > 0): break
        # end for

        print('\nCreating metadata index for {} mseed files..'.format(len(self.meta_list)))

        for row in tqdm(self.meta_list):
            idx, nc, sc, lc, cc, st, et = row

            if (type(self.tree[nc][sc][lc][cc]) != index.Index):
                self.tree[nc][sc][lc][cc] = index.Index()
            # end if
            self.tree[nc][sc][lc][cc].insert(idx, (st, 1, et, 1))
            # end for
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
