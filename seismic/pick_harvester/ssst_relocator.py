from ordered_set import OrderedSet as set
import numpy as np
from obspy import UTCDateTime
from scipy.spatial import cKDTree
import pyproj
from pyproj import Geod
from collections import defaultdict
from seismic.pick_harvester.utils import split_list

from seismic.pick_harvester.parametric_data import ParametricData
import importlib
import sys, os
from os.path import abspath, dirname

from seismic.pick_harvester.travel_time import TTInterpolator
from seismic.pick_harvester.ellipticity import Ellipticity
import traceback
import psutil
import h5py
from itertools import combinations

class SSSTRelocator(ParametricData):
    def __init__(self, csv_catalog, auto_pick_files=[], auto_pick_phases=[],
                 events_only=False, phase_list='P Pg Pb Pn S Sg Sb Sn', temp_dir='./'):
        super(SSSTRelocator, self).__init__(csv_catalog, auto_pick_files, auto_pick_phases,
                                            events_only, phase_list, temp_dir)
        self.EARTH_RADIUS = 6371. #km
        self.DEG2KM = np.pi/180 * self.EARTH_RADIUS
        self.kdtree = None
        self._lonlatalt2xyz = None
        self.local_events_indices = np.array(split_list(np.arange(len(self.events)), self.nproc)[self.rank], dtype='i4')
        self.local_arrivals_indices = np.array(split_list(np.arange(len(self.arrivals)), self.nproc)[self.rank], dtype='i4')
        self.ecdist = None
        self.edepth_km = None
        self.elat = None
        self.azim = None
        self.elev_corr = None
        self.ellip_corr = None
        self.vsurf = None
        self._tti = TTInterpolator()
        self._ellipticity = Ellipticity()
        self._geod = Geod(a=180/np.pi, f=0)

        self._initialize_kdtree()

        # initialize surface velocity array
        vs_surf = 3.46        # Sg velocity km/s for elevation corrections
        vp_surf = 5.8         # Pg velocity km/s for elevation corrections
        pidx = self.is_P()
        sidx = self.is_S()
        self.vsurf = np.zeros(len(self.arrivals), dtype='f4')
        self.vsurf[pidx] = vp_surf
        self.vsurf[sidx] = vs_surf
    # end func

    def is_P(self):
        result = self.arrivals['phase'].astype('S1') == b'P'
        return result
    # end func

    def is_S(self):
        result = self.arrivals['phase'].astype('S1') == b'S'
        return result
    # end func

    def coalesce_network_codes(self):
        self.arrivals['net'] = np.load('coalesced_net.npy')
        return

        MDIST = 500 # maximum distance in metres
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
                        if(dist < MDIST):
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
        #if(self.rank==0): np.save('coalesced_net.npy', self.arrivals['net'])
    # end func

    def compute_essentials(self):
        elon = self.events['lon'][self.event_id_to_idx[self.arrivals['event_id']]]
        elat = self.events['lat'][self.event_id_to_idx[self.arrivals['event_id']]]
        slon = self.arrivals['lon']
        slat = self.arrivals['lat']

        self.edepth_km = self.events['depth_km'][self.event_id_to_idx[self.arrivals['event_id']]]
        self.elat = self.events['lat'][self.event_id_to_idx[self.arrivals['event_id']]]

        azim, _, ecdist = self._geod.inv(elon[self.local_arrivals_indices],
                                         elat[self.local_arrivals_indices],
                                         slon[self.local_arrivals_indices],
                                         slat[self.local_arrivals_indices])
        self.azim = self._sync_var(np.array(azim.astype('f4')))
        self.ecdist = self._sync_var(np.array(ecdist.astype('f4')))
    # end func

    def compute_elevation_correction(self):
        vsurf = self.vsurf[self.local_arrivals_indices]
        elev_km = self.arrivals['elev_m'][self.local_arrivals_indices] / 1e3

        dtdd = self._tti.get_dtdd(self.arrivals['phase'][self.local_arrivals_indices],
                                  self.ecdist[self.local_arrivals_indices],
                                  self.edepth_km[self.local_arrivals_indices])
        elev_corr = vsurf * (dtdd / self.DEG2KM)
        elev_corr = np.power(elev_corr, 2.)
        elev_corr[elev_corr > 1.] = 1./elev_corr[elev_corr > 1.]
        elev_corr = np.sqrt(1. - elev_corr)
        elev_corr = elev_corr * elev_km / vsurf

        elev_corr = np.array(elev_corr.astype('f4'))
        self.elev_corr = self._sync_var(elev_corr)
    # end func

    def compute_ellipticity_correction(self):
        ellip_corr = self._ellipticity.get_correction(self.arrivals['phase'][self.local_arrivals_indices],
                                                      self.ecdist[self.local_arrivals_indices],
                                                      self.edepth_km[self.local_arrivals_indices],
                                                      self.elat[self.local_arrivals_indices],
                                                      self.azim[self.local_arrivals_indices])
        ellip_corr = np.array(ellip_corr.astype('f4'))
        self.ellip_corr = self._sync_var(ellip_corr)
    # end func

    def _initialize_kdtree(self, ellipsoidal_distance=False):
        ER = self.EARTH_RADIUS * 1e3 #m

        elons = self.events['lon']
        elats = self.events['lat']
        ealts = -self.events['depth_km'] * 1e3 #m
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
    if(0):
        sr = SSSTRelocator('./small_merge_catalogues_output.csv',
                           auto_pick_files=['small_p_combined.txt', 'small_s_combined.txt'],
                           auto_pick_phases=['P', 'S'],
                           events_only=False)
        if(sr.rank==0): print('computing essentials..')
        sr.compute_essentials()
        if(sr.rank==0): print('computing elevation corrections..')
        sr.compute_elevation_correction()
        if(sr.rank==0): print('computing ellipticity corrections..')
        sr.compute_ellipticity_correction()
        if(sr.rank == 0): print(sr.azim, sr.ecdist, sr.elev_corr, sr.ellip_corr)

        test_elev_corr = sr._sync_var(np.array(split_list(sr.elev_corr, sr.nproc)[sr.rank]))
        #test_elev_corr = sr._sync_var_h5(np.array(split_list(sr.elev_corr, sr.nproc)[sr.rank]))

        try:
            assert np.all(sr.elev_corr == test_elev_corr), 'Sync-var failed on rank {}'.format(sr.rank)
        except:
            np.savetxt('{}.synced.txt'.format(sr.rank), test_elev_corr)
            np.savetxt('{}.ellip_corr.txt'.format(sr.rank), sr.elev_corr)
        # end try
    else:
        sr = SSSTRelocator('./merge_catalogues_output.csv',
                           auto_pick_files=['p_combined.txt', 's_combined.txt'],
                           auto_pick_phases=['P', 'S'],
                           events_only=False)

        if(sr.rank==0): print('coalescing network codes..')
        sr.coalesce_network_codes()
        if(sr.rank==0): print('computing essentials..')
        sr.compute_essentials()
        if(sr.rank==0): print('computing elevation corrections..')
        sr.compute_elevation_correction()
        if(sr.rank==0): print('computing ellipticity corrections..')
        sr.compute_ellipticity_correction()

        print('Rank {}: memory used: {}'.format(sr.rank,
                                                round(psutil.Process().memory_info().rss / 1024. / 1024., 2)))

        if(sr.rank == 0): print(sr.azim, sr.ecdist, sr.elev_corr, sr.ellip_corr)

        for name, dtype in zip(sr.arrival_fields['names'], sr.arrival_fields['formats']):
            test_var = sr._sync_var(sr.arrivals[name][sr.local_arrivals_indices])

            assert np.all(test_var == sr.arrivals[name]), 'Sync-{} failed on rank {}'.format(name, sr.rank)
        # end for
    # end if
# end if
