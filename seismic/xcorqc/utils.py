from mpi4py import MPI
import os
import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict
from obspy import UTCDateTime, read_inventory, Inventory, Stream
from obspy.geodetics.base import gps2dist_azimuth
from tempfile import SpooledTemporaryFile

# define utility functions
def rtp2xyz(r, theta, phi):
    xout = np.zeros((r.shape[0], 3))
    rst = r * np.sin(theta)
    xout[:, 0] = rst * np.cos(phi)
    xout[:, 1] = rst * np.sin(phi)
    xout[:, 2] = r * np.cos(theta)
    return xout
# end func

def xyz2rtp(x, y, z):
    rout = np.zeros((x.shape[0], 3))
    tmp1 = x * x + y * y
    tmp2 = tmp1 + z * z
    rout[0] = np.sqrt(tmp2)
    rout[1] = np.arctan2(sqrt(tmp1), z)
    rout[2] = np.arctan2(y, x)
    return rout
# end func

def getStationInventory(master_inventory, inventory_cache, netsta):
    netstaInv = None
    if (master_inventory):
        if (inventory_cache is None): inventory_cache = defaultdict(list)
        net, sta = netsta.split('.')

        if (isinstance(inventory_cache[netsta], Inventory)):
            netstaInv = inventory_cache[netsta]
        else:
            inv = master_inventory.select(network=net, station=sta)
            if(len(inv.networks)):
                inventory_cache[netsta] = inv
                netstaInv = inv
            # end if
        # end if
    # end if

    return netstaInv, inventory_cache
# end func

def split_list(lst, npartitions):
    k, m = divmod(len(lst), npartitions)
    return [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(npartitions)]
# end func

def drop_bogus_traces(st, sampling_rate_cutoff=1):
    """
    Removes spurious traces with suspect sampling rates.
    :param st: Obspy Stream
    :param sampling_rate_cutoff: sampling rate threshold
    :return: Input stream is updated inplace
    """
    badTraces = [tr for tr in st if tr.stats.sampling_rate < sampling_rate_cutoff]

    for tr in badTraces: st.remove(tr)
# end func

def _get_stream_00T(fds, net, sta, cha, start_time, end_time,
                      baz=None, trace_count_threshold=200,
                      logger=None, verbose=1):

    stations = fds.get_stations(start_time, end_time, network=net, station=sta)

    stations_nch = [s for s in stations if 'N' == s[3][-1].upper() or '1' == s[3][-1]]  # only N channels
    stations_ech = [s for s in stations if 'E' == s[3][-1].upper() or '2' == s[3][-1]]  # only E channels

    stt = Stream()
    if (len(stations_nch) > 0 and len(stations_nch) == len(stations_ech)):
        for codesn, codese in zip(stations_nch, stations_ech):

            stn = fds.get_waveforms(codesn[0], codesn[1], codesn[2], codesn[3],
                                    start_time,
                                    end_time,
                                    trace_count_threshold=trace_count_threshold)
            ste = fds.get_waveforms(codese[0], codese[1], codese[2], codese[3],
                                    start_time,
                                    end_time,
                                    trace_count_threshold=trace_count_threshold)

            if (len(stn) == 0): continue
            if (len(ste) == 0): continue
            
            drop_bogus_traces(stn)
            drop_bogus_traces(ste)

            # Merge station data. Note that we don't want to fill gaps; the
            # default merge() operation creates masked numpy arrays, which we can use
            # to detect and ignore windows that have gaps in their data.
            try:
                stn.merge()
                ste.merge()

                max_start_time = np.max([stn[0].stats.starttime, ste[0].stats.starttime])
                min_end_time   = np.min([stn[0].stats.endtime, ste[0].stats.endtime])

                stn = stn.slice(starttime = max_start_time, endtime= min_end_time)
                ste = ste.slice(starttime=max_start_time, endtime=min_end_time)
            except Exception as e:
                if logger: logger.warning('\tFailed to merge traces..')
                st = None
                raise
            # end try

            bazr = np.radians(baz)
            tdata = -ste[0].data * np.cos(bazr) + stn[0].data * np.sin(bazr)

            stt_curr = ste.copy()
            stt_curr[0].data = tdata
            #stt_curr[0].stats.channel = '00T'

            stt += stt_curr
        # end for
    # end if

    return stt
# end func

def get_stream(fds, net, sta, cha, start_time, end_time,
               baz=None, trace_count_threshold=200,
               logger=None, verbose=1):

    if (cha == '00T'): return _get_stream_00T(fds, net, sta, cha, start_time, end_time,
                                              baz=baz, trace_count_threshold=trace_count_threshold,
                                              logger=logger, verbose=verbose)
    st = Stream()
    stations = fds.get_stations(start_time, end_time, network=net, station=sta)
    for codes in stations:
        if (cha != codes[3]): continue
        st += fds.get_waveforms(codes[0], codes[1], codes[2], codes[3], start_time,
                               end_time, trace_count_threshold=trace_count_threshold)
    # end for

    drop_bogus_traces(st)

    if (verbose > 2):
        if logger: logger.debug('\t\tData Gaps:')
        st.print_gaps()  # output sent to stdout; fix this
        print ("\n")
    # end if

    # Merge station data. Note that we don't want to fill gaps; the
    # default merge() operation creates masked numpy arrays, which we can use
    # to detect and ignore windows that have gaps in their data.
    try:
        st.merge()
    except Exception as e:
        if logger: logger.warning('\tFailed to merge traces..')
        st = None
        raise
    # end try

    return st
# end func

class ProgressTracker:
    def __init__(self, output_folder, restart_mode=False):
        self.output_folder = output_folder
        self.restart_mode = restart_mode

        self.comm = MPI.COMM_WORLD
        self.nproc = self.comm.Get_size()
        self.rank = self.comm.Get_rank()

        self.prev_progress = 0 # progress from a previous run
        self.progress = 0
        self.proc_fn = os.path.join(output_folder, 'prog.%d.txt' % (self.rank))

        if(self.restart_mode):
            if(not os.path.exists(self.proc_fn)):
                raise Exception('Progress file (%s) not found'%(self.proc_fn))
            # end if

            self.prev_progress = int(open(self.proc_fn).read())
        # end if
    # end func

    def increment(self):
        self.progress += 1
        if(self.restart_mode and (self.prev_progress > 0) and (self.progress < self.prev_progress)):
            return False
        else:
            tmpfn = self.proc_fn + '.tmp'
            f = open(tmpfn, 'w+')
            f.write(str(self.progress))
            f.close()
            os.rename(tmpfn, self.proc_fn)

            return True
        # end if
    # end func
# end class

class SpooledXcorrResults:
    """
    Spooled storage for cross-correlations. Stacked cross-correlations computed were previously
    gathered in memory, before being written to netCDF4 files at the end. Because we now need to
    output all cross-correlations, unstacked, the memory requirements have jumped by a factor of
    ~26. We now write the CCs to a spooled storage as they are being computed.
    """
    def __init__(self, ncols, dtype=np.float32, max_size_mb=2048, prefix='', dir=None):
        self._prefix = prefix
        self._ncols = ncols
        self._nrows = 0
        self._dtype = dtype
        self._max_size_mb = max_size_mb

        self._file = SpooledTemporaryFile(prefix = self._prefix, mode = 'w+b', max_size = max_size_mb * 1024**2, dir=dir)
    # end func
    
    @property
    def ncols(self):
        return self._ncols
    # end func
    
    @property
    def nrows(self):
        return self._nrows
    # end func

    def write_row(self, row_array):
        assert(row_array.dtype == self._dtype) 
        assert(len(row_array.shape) == 1)
        assert(row_array.shape[0] == self._ncols)

        self._file.write(row_array.data)
        
        self._nrows += 1
    # end func

    def read_row(self, row_idx):
        if(row_idx < self._nrows):
            seek_loc = row_idx * self._ncols * np.dtype(self._dtype).itemsize
            self._file.seek(seek_loc)

            nbytes = self._ncols * np.dtype(self._dtype).itemsize
            row = np.frombuffer(self._file.read(nbytes), dtype=self._dtype)
        
            return row
        else:
            return None
    # end func

    def __del__(self):
        self._file.close()
    # end func
# end class

