from mpi4py import MPI
import os
import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict
from obspy import UTCDateTime, read_inventory, Inventory, Stream
from obspy.geodetics.base import gps2dist_azimuth
from tempfile import SpooledTemporaryFile
from scipy.interpolate import interp1d
from seismic.ASDFdatabase.FederatedASDFDataSet import FederatedASDFDataSet
from seismic.misc import rtp2xyz
from seismic.misc import get_git_revision_hash, rtp2xyz, split_list
import os, psutil

class Dataset:
    def __init__(self, asdf_file_name, netsta_list='*'):

        self._data_path = asdf_file_name
        self._earth_radius = 6371  # km

        self.fds = FederatedASDFDataSet(asdf_file_name)
        # Gather station metadata
        netsta_list_subset = set(netsta_list.split(' ')) if netsta_list != '*' else netsta_list
        self.netsta_list = []
        self.metadata = defaultdict(list)

        rtps = []
        for netsta in list(self.fds.unique_coordinates.keys()):
            if(netsta_list_subset != '*'):
                if netsta not in netsta_list_subset:
                    continue

            self.netsta_list.append(netsta)
            self.metadata[netsta] = self.fds.unique_coordinates[netsta]

            rtps.append([self._earth_radius,
                         np.radians(90 - self.metadata[netsta][1]),
                         np.radians(self.metadata[netsta][0])])
        # end for

        if(len(rtps) == 0):
            assert 0, 'No station-pairs found due to missing stations. Aborting..'
        # end if

        rtps = np.array(rtps)
        xyzs = rtp2xyz(rtps[:, 0], rtps[:, 1], rtps[:, 2])

        self._tree = cKDTree(xyzs)
        self._cart_location = defaultdict(list)
        for i, ns in enumerate(self.netsta_list):
            self._cart_location[ns] = xyzs[i, :]
        # end for
    # end func

    def get_closest_stations(self, netsta, other_dataset, nn=1):
        assert isinstance(netsta, str), 'station_name must be a string'
        assert isinstance(other_dataset, Dataset), 'other_dataset must be an instance of Dataset'
        netsta = netsta.upper()

        assert netsta in self.netsta_list, '%s not found'%(netsta)

        d, l = other_dataset._tree.query(self._cart_location[netsta], nn)

        if isinstance(l, int):
            l = np.array([l])

        l = l[l<len(other_dataset.netsta_list)]

        if isinstance(l, int):
            l = np.array([l])

        assert len(l), 'No stations found..'

        return list(np.array(other_dataset.netsta_list)[l])
    # end func

    def get_unique_station_pairs(self, other_dataset, nn=1, require_overlap=True,
                                 min_distance_km=None, max_distance_km=None):
        pairs = set()
        for ns1 in self.netsta_list:
            ns2list = None
            if (nn != -1):
                if self == other_dataset:
                    ns2list = set(self.get_closest_stations(ns1, other_dataset, nn=nn + 1))
                    if ns1 in ns2list:
                        ns2list.remove(ns1)
                    ns2list = list(ns2list)
                else:
                    ns2list = self.get_closest_stations(ns1, other_dataset, nn=nn)
            else:
                ns2list = other_dataset.netsta_list
            # end if

            for ns2 in ns2list:
                pairs.add((ns1, ns2))
            # end for
        # end for

        pairs_subset = set()
        for item in pairs:
            if(item[0] == item[1]): continue

            dup_item = (item[1], item[0])
            if(dup_item not in pairs_subset and item not in pairs_subset):
                pairs_subset.add(item)
            # end if
        # end if

        result_pairs = pairs_subset

        # cull pairs based on temporal overlap if specified
        # (note: gaps are not considered)
        if(require_overlap):
            overlapped_pairs = set()
            range_cache = defaultdict(tuple)
            for ns1, ns2 in pairs_subset:
                st1 = et1 = st2 = et2 = None

                if(ns1 in range_cache.keys()):
                    st1, et1 = range_cache[ns1]
                else:
                    net1, sta1 = ns1.split('.')
                    st1, et1 = self.fds.get_global_time_range(net1, sta1)
                    range_cache[ns1] = (st1, et1)
                # end if

                if(ns2 in range_cache.keys()):
                    st2, et2 = range_cache[ns2]
                else:
                    net2, sta2 = ns2.split('.')
                    st2, et2 = other_dataset.fds.get_global_time_range(net2, sta2)
                    range_cache[ns2] = (st2, et2)
                # end if

                # check for temporal overlap
                if((st1 <= et2) and (st2 <= et1)):
                    overlapped_pairs.add((ns1, ns2))
                # end if
            # end for
            #for k,v in range_cache.items(): print(k, v)

            result_pairs = overlapped_pairs
        # end if

        # filter by distance between pairs
        if((min_distance_km is not None) or (max_distance_km is not None)):
            dist_filtered_pairs = set()
            for ns1, ns2 in result_pairs:
                lon1, lat1 = self.metadata[ns1]
                lon2, lat2 = other_dataset.metadata[ns2]

                dist, _, _ = gps2dist_azimuth(lat1, lon1, lat2, lon2)
                dist /= 1e3 # km

                if(min_distance_km is not None):
                    if(dist < min_distance_km): continue
                if(max_distance_km is not None):
                    if(dist > max_distance_km): continue

                dist_filtered_pairs.add((ns1, ns2))
                #print('.'.join([ns1, ns2]), lon1, lat1, lon2, lat2, dist)
            # end for

            result_pairs = dist_filtered_pairs
        # end if

        return list(result_pairs)
    # end func
# end class

def read_location_preferences(location_preferences_fn):
    result = defaultdict(lambda: None)

    if(location_preferences_fn):
        pref_list = open(location_preferences_fn, 'r').readlines()

        for pref in pref_list:
            pref = pref.strip()
            if (len(pref)):
                try:
                    netsta, loc = pref.split()
                    net, sta = netsta.split('.')

                    result[netsta] = loc
                except Exception as e:
                    print(str(e))
                    assert 0, 'Error parsing: {}'.format(pref)
                # end try
            # end if
        # end for
    # end if
    return result
# end func

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

class MemoryTracker:
    def __init__(self, burnin_steps=7, outlier_factor=10, logger=None):
        assert burnin_steps > 1, 'Burnin-steps must be > 1'

        self.burnin_steps = burnin_steps
        self.outlier_factor = outlier_factor
        self.logger = logger

        self.usage_mb = np.zeros(self.burnin_steps)
        self.usage_mb[-1] = psutil.Process(os.getpid()).memory_info().rss / (1024 ** 2)
    # end func

    def update(self):
        self.usage_mb[:-1] = self.usage_mb[1:]
        self.usage_mb[-1] = psutil.Process(os.getpid()).memory_info().rss / (1024 ** 2)

        if(np.sum(self.usage_mb > 0) < self.burnin_steps): return

        # track outliers based on scaled median-absolute-deviations
        med = np.median(self.usage_mb)
        d = np.abs(self.usage_mb - med)
        s = d / med if med else np.zeros(len(d))

        outlier_indices = np.where(s >= self.outlier_factor)[0]
        for oi in outlier_indices:
            # only report anomalously large outliers
            if(self.usage_mb[oi] > med):
                msg = 'Anomalous memory usage detected: median({:3.2f} MB), current({:3.2f} MB)'.\
                    format(med, self.usage_mb[oi])
                if(self.logger):
                    self.logger.warning(msg)
                else:
                    print('Warning: {}'.format(msg))
                # end if
            # end if
        # end for
    # end func
# end class
