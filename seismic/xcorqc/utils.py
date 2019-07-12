from mpi4py import MPI
import os
import numpy as np
from scipy.spatial import cKDTree

from obspy import UTCDateTime, read_inventory, Inventory
from obspy.geodetics.base import gps2dist_azimuth

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
                      baz=None, automerge=False, trace_count_threshold=200,
                      logger=None, verbose=1):

    stations = fds.get_stations(start_time, end_time, network=net, station=sta)

    stations_nch = [s for s in stations if 'N' == s[3][-1].upper() or '1' == s[3][-1]]  # only N channels
    stations_ech = [s for s in stations if 'E' == s[3][-1].upper() or '2' == s[3][-1]]  # only E channels

    stt = None
    if (len(stations_nch) > 0 and len(stations_nch) == len(stations_ech)):
        for codesn, codese in zip(stations_nch, stations_ech):

            stn = fds.get_waveforms(codesn[0], codesn[1], codesn[2], codesn[3],
                                    start_time,
                                    end_time,
                                    automerge=automerge,
                                    trace_count_threshold=trace_count_threshold)
            ste = fds.get_waveforms(codese[0], codese[1], codese[2], codese[3],
                                    start_time,
                                    end_time,
                                    automerge=True,
                                    trace_count_threshold=trace_count_threshold)

            if (len(stn) == 0): continue
            if (len(ste) == 0): continue

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
            tdata = - ste[0].data * np.cos(bazr) + stn[0].data * np.sin(bazr)

            stt = ste.copy()
            stt[0].data = tdata
            #stt[0].stats.channel = '00T'
        # end for
    # end if

    if (stt and len(stt)):
        drop_bogus_traces(stt)
    # end if

    return stt
# end func

def get_stream(fds, net, sta, cha, start_time, end_time,
               baz=None, automerge=False, trace_count_threshold=200,
               logger=None, verbose=1):

    if (cha == '00T'): return _get_stream_00T(fds, net, sta, cha, start_time, end_time,
                                              baz=baz, automerge=automerge,
                                              trace_count_threshold=trace_count_threshold,
                                              logger=logger, verbose=verbose)
    st = None
    stations = fds.get_stations(start_time, end_time, network=net, station=sta)
    for codes in stations:
        if (cha != codes[3]): continue
        st = fds.get_waveforms(codes[0], codes[1], codes[2], codes[3], start_time,
                               end_time, automerge=automerge, trace_count_threshold=trace_count_threshold)

        if (len(st) == 0): continue
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
    # end for

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
        if(self.restart_mode and (self.prev_progress > 0) and (self.progress <= self.prev_progress)):
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