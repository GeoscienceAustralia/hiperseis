from mpi4py import MPI
import os
import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict
from obspy import UTCDateTime, read_inventory, Inventory, Stream
from obspy.geodetics.base import gps2dist_azimuth
from tempfile import SpooledTemporaryFile
from scipy.interpolate import interp1d

def getStationInventory(master_inventory, inventory_cache, netsta, location_preferences_dict):
    netstaInv = None
    if (master_inventory):
        if (inventory_cache is None): inventory_cache = defaultdict(list)
        net, sta = netsta.split('.')

        if (isinstance(inventory_cache[netsta], Inventory)):
            netstaInv = inventory_cache[netsta]
        else:
            inv = master_inventory.select(network=net, station=sta, location=location_preferences_dict[netsta])
            if(len(inv.networks)):
                inventory_cache[netsta] = inv
                netstaInv = inv
            # end if
        # end if
    # end if

    return netstaInv, inventory_cache
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

def fill_gaps(data, dt, max_gap_seconds=3):
    """
    Fills gaps <= max_gap_seconds in length via linear interpolation
    """

    if (not np.ma.is_masked(data)): return data

    gaps = np.ma.clump_masked(data)

    for gap in gaps:
        s, e = gap.start, gap.stop
        if ((e - s) * dt > max_gap_seconds):
            return data
    # end for

    valIndices = np.ma.where(~data.mask)[0]

    if (valIndices.shape[0] < 2): return data

    tOld = valIndices * dt
    tNew = np.arange(data.shape[0]) * dt

    io = interp1d(tOld, data[valIndices], bounds_error=False, fill_value='extrapolate')

    result = io(tNew)
    return result
# end func

def _get_stream_00T(fds, net, sta, cha, start_time, end_time, location_preferences_dict,
                      baz=None, trace_count_threshold=200,
                      logger=None, verbose=1):

    netsta = net + '.' + sta
    stations = fds.get_stations(start_time, end_time, network=net, station=sta,
                                location=location_preferences_dict[netsta])

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

def get_stream(fds, net, sta, cha, start_time, end_time, location_preferences_dict,
               baz=None, trace_count_threshold=200,
               logger=None, verbose=1):

    if (cha == '00T'): return _get_stream_00T(fds, net, sta, cha, start_time, end_time, location_preferences_dict,
                                              baz=baz, trace_count_threshold=trace_count_threshold,
                                              logger=logger, verbose=verbose)
    st = Stream()
    netsta = net + '.' + sta
    stations = fds.get_stations(start_time, end_time, network=net, station=sta,
                                location=location_preferences_dict[netsta])
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

class SpooledMatrix:
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

    def get_matrix(self):
        mat = np.empty([])
        if(self.nrows>0 and self.ncols>0):
            mat = np.zeros((self.nrows, self.ncols), dtype=self._dtype)
            for i in np.arange(self.nrows): mat[i, :] = self.read_row(i)
        # enf if

        return mat
    # end if

    def __del__(self):
        self._file.close()
    # end func
# end class

