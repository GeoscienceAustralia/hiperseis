from obspy.core import UTCDateTime
import numpy as np
from rtree import index
from glob import glob
from mpi4py import MPI
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
    def __init__(self, mseed_folder, pattern):
        self.mseed_folder = mseed_folder
        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()
        self.tree = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))

        work_load = None
        offsets = None
        if (self.rank == 0):
            self.mseed_files = np.array(sorted(glob(os.path.join(self.mseed_folder, pattern))))

            #self.mseed_files = self.mseed_files[:1000]

            work_load = split_list(self.mseed_files, self.nproc)
            counts = np.array([len(item) for item in work_load])
            offsets = np.append(0, np.cumsum(counts[:-1]))
        # end if
        self.local_mseed_files = self.comm.scatter(work_load)
        self.offsets = self.comm.bcast(offsets, root=0)

        print('Reading metadata from mseed files..')
        meta_list = []
        for i, mseed_file in enumerate(tqdm(self.local_mseed_files, desc='Rank {}'.format(self.rank))):
            st = None
            try:
                st = read(mseed_file, headonly=True)
            except:
                print("Failed to read {}. Moving along..". format(mseed_file))
                continue
            # end try

            for tr in st:
                nc, sc, lc, cc, st, et = \
                    tr.stats.network, tr.stats.station, tr.stats.location, \
                    tr.stats.channel, tr.stats.starttime.timestamp, \
                    tr.stats.endtime.timestamp

                meta_list.append([self.offsets[self.rank] + i, nc, sc, lc, cc, st, et])
            # end for
            #if (i > 0): break
        # end for

        meta_list = self.comm.gather(meta_list, root=0)
        if (self.rank == 0):
            meta_list = [item for ritem in meta_list for item in ritem]  # flatten list of lists

            print('\nCreating metadata index for {} traces..'.format(len(meta_list)))

            for row in tqdm(meta_list):
                idx, nc, sc, lc, cc, st, et = row

                if (type(self.tree[nc][sc][lc][cc]) != index.Index):
                    self.tree[nc][sc][lc][cc] = index.Index()
                # end if
                self.tree[nc][sc][lc][cc].insert(idx, (st, 1, et, 1))
            # end for
        # end if
    # end func

    def get_waveforms(self, net, sta, loc, cha, st: UTCDateTime, et: UTCDateTime):
        assert self.rank == 0, 'This function is only accessible from Rank 0. Aborting..'

        epsilon = 1e-5
        st_ts = st.timestamp + epsilon
        et_ts = et.timestamp - epsilon

        result = Stream([])
        try:
            target_index = self.tree[net][sta][loc][cha]

            if(type(target_index) == index.Index):
                file_indices = np.array(list(target_index.intersection((st_ts, 1, et_ts, 1))), dtype='i4')

                # since file names are repeated for multiple traces, we need a unique set
                for mfile in set(self.mseed_files[file_indices]):
                    result += read(mfile).select(network=net,
                                                 station=sta,
                                                 location=loc,
                                                 channel=cha).slice(st, et, nearest_sample=False).copy()
                # end for
            # end if
        except Exception as e:
            print(str(e))
        # end try

        return result
    # end func

    def get_stations(self, st:UTCDateTime, et:UTCDateTime, net=None, sta=None, loc=None, cha=None):
        assert self.rank == 0, 'This function is only accessible from Rank 0. Aborting..'

        epsilon = 1e-5
        st_ts = st.timestamp + epsilon
        et_ts = et.timestamp - epsilon

        _net = _sta = _loc = _cha = None

        if(net==None): _net = self.tree.keys()
        else: _net = [net]

        result = []
        for nc in _net:
            if (sta==None): _sta = self.tree[nc].keys()
            else: _sta = [sta]

            for sc in _sta:
                if (loc==None): _loc = self.tree[nc][sc].keys()
                else: _loc = [loc]

                for lc in _loc:
                    if (cha==None): _cha = self.tree[nc][sc][lc].keys()
                    else: _cha= [cha]

                    for cc in _cha:
                        target_index = self.tree[nc][sc][lc][cc]
                        if (type(target_index) == index.Index):
                            entries = list(target_index.intersection((st_ts, 1, et_ts, 1)))

                            if(len(entries)): result.append((nc, sc, lc, cc))
                        # end if
                    # end for
                # end for
            # end for
        # end for

        return result
    # end func

    def get_time_range(self, net, sta, loc, cha):
        assert self.rank == 0, 'This function is only accessible from Rank 0. Aborting..'

        target_index = self.tree[net][sta][loc][cha]

        if (type(target_index) == index.Index):
            bounds = target_index.bounds
            return UTCDateTime(bounds[0]), UTCDateTime(bounds[2])
        # end if

        return MAX_DATE, MIN_DATE
    # end func
# end func

if __name__=="__main__":
    msi = MseedIndex('/g/data/ha3/ac5759/semi-perm-iris', '*.mseed')

    if(msi.rank == 0):
        print(msi.tree['AU'].keys())
        print(msi.tree['AU']['AXCOZ'].keys())
        r = msi.get_waveforms('AU', 'AXCOZ', '00', 'HHZ', UTCDateTime(2021, 1, 3, 22, 9, 26), UTCDateTime(2021, 1, 4, 0, 0))
        print(r)
        print(msi.get_time_range('AU', 'AXCOZ', '00', 'HHN'))
    # end if
# end if
