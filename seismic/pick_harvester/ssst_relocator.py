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
        self.vsurf = None
        self.ett = None
        self.residual = None
        self._source_enum = defaultdict(int)
        self._arrival_source = None
        self._tti = TTInterpolator()
        self._ellipticity = Ellipticity()
        self._geod = Geod(a=180/np.pi, f=0)

        self._label_arrivals()
        self._coalesce_network_codes()
        self._initialize_kdtree()

        # compute empirical tt (need only be computed once)
        self.ett = self.arrivals['arrival_ts'] - \
                   self.events['origin_ts'][self.event_id_to_idx[self.arrivals['event_id']]]

        # initialize surface velocity array (only needs to be computed once, since
        # elevations are fixed and P/S arrivals are not interchanged)
        vs_surf = 3.46        # Sg velocity km/s for elevation corrections
        vp_surf = 5.8         # Pg velocity km/s for elevation corrections
        pidx = self.is_P(self.arrivals['phase'])
        sidx = self.is_S(self.arrivals['phase'])
        self.vsurf = np.zeros(len(self.arrivals), dtype='f4')
        self.vsurf[pidx] = vp_surf
        self.vsurf[sidx] = vs_surf
    # end func

    def is_P(self, phase):
        result = phase.astype('S1') == b'P'
        return result
    # end func

    def is_S(self, phase):
        result = phase.astype('S1') == b'S'
        return result
    # end func

    def _compute_residual_helper(self, phase, mask=None):
        if(mask is None): mask = np.ones(len(phase), dtype='i4')

        # event indices for local arrivals
        indices = self.event_id_to_idx[self.arrivals['event_id'][self.local_arrivals_indices]]
        elon = self.events['lon'][indices][mask]
        elat = self.events['lat'][indices][mask]
        edepth_km = self.events['depth_km'][indices][mask]

        slon = self.arrivals['lon'][self.local_arrivals_indices][mask]
        slat = self.arrivals['lat'][self.local_arrivals_indices][mask]
        vsurf = self.vsurf[self.local_arrivals_indices][mask]
        elev_km = self.arrivals['elev_m'][self.local_arrivals_indices][mask] / 1e3

        azim, _, ecdist = self._geod.inv(elon, elat, slon, slat)

        elev_corr = self.compute_elevation_correction(phase, vsurf, elev_km, ecdist, edepth_km)
        ellip_corr = self.compute_ellipticity_correction(phase, ecdist, edepth_km, elat, azim)
        ptt = self._tti.get_tt(phase, ecdist, edepth_km) + elev_corr + ellip_corr

        residual = self.ett[self.local_arrivals_indices][mask] - ptt.astype('f8')

        if(self.rank == 0): print(elev_corr, ellip_corr, residual)

        return residual
    # end func

    def compute_residual(self):
        phase = self.arrivals['phase'][self.local_arrivals_indices]

        local_residual = self._compute_residual_helper(phase)

        self.residual = self._sync_var(local_residual)
    # end func

    def redefine_phases(self):
        phase = self.arrivals[self.local_arrivals_indices]
        is_p = self.is_P(phase)
        is_s = self.is_S(phase)
        p_indices = np.argwhere(is_p)
        s_indices = np.argwhere(is_s)

        p_residuals = np.zeros((len(p_indices), len(self.p_phases)))
        s_residuals = np.zeros((len(s_indices), len(self.s_phases)))

        test_p_phase = np.zeros(len(p_indices), self.arrivals['phase'].dtype)
        for ip, p in enumerate(self.p_phases):
            test_p_phase[:] = p
            p_residuals[:, ip] = self._compute_residual_helper(test_p_phase)
        # end for

        test_s_phase = np.zeros(len(s_indices), self.arrivals['phase'].dtype)
        for ip, p in enumerate(self.s_phases):
            test_s_phase[:] = p
            s_residuals[:, ip] = self._compute_residual_helper(test_s_phase)
        # end for

        new_phase = np.zeros(len(self.local_arrivals_indices), dtype=self.arrivals['phase'].dtype)
        new_phase[is_p] = b'Px'
        new_phase[is_s] = b'Sx'

        new_phase[p_indices] = self.p_phases[np.argmin(p_residuals, axis=1)]
        new_phase[s_indices] = self.s_phases[np.argmin(s_residuals, axis=1)]
    # end func

    def compute_elevation_correction(self, phase, vsurf, elev_km, ecdist, edepth_km):
        dtdd = self._tti.get_dtdd(phase, ecdist, edepth_km)
        elev_corr = vsurf * (dtdd / self.DEG2KM)
        elev_corr = np.power(elev_corr, 2.)
        elev_corr[elev_corr > 1.] = 1./elev_corr[elev_corr > 1.]
        elev_corr = np.sqrt(1. - elev_corr)
        elev_corr = elev_corr * elev_km / vsurf

        elev_corr = np.array(elev_corr.astype('f4'))

        return elev_corr
    # end func

    def compute_ellipticity_correction(self, phase, ecdist, edepth_km, elat, azim):
        ellip_corr = self._ellipticity.get_correction(phase, ecdist, edepth_km, elat, azim)
        ellip_corr = np.array(ellip_corr.astype('f4'))

        return ellip_corr
    # end func

    def compute_ssst_correction(self):
        pass
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

    def _label_arrivals(self):
        sources = set(self.events['source'])
        for i, source in enumerate(sources): self._source_enum[source] = i + 1

        if(self.rank == 0): print('Labelling arrivals by event source {}..'.format(self._source_enum.items()))

        arrival_source = np.zeros(self.local_arrivals_indices.shape, dtype='i4')
        local_event_ids = self.arrivals['event_id'][self.local_arrivals_indices]
        for source in sources:
            enum = self._source_enum[source]
            sids = np.argwhere(self.events['source'] == source).flatten()
            source_eids = self.events['event_id'][sids]

            arrival_source[np.isin(local_event_ids, source_eids)] = enum
        # end for

        assert np.all(arrival_source > 0), 'Arrivals found with no corresponding event-ids..'

        self._arrival_source = self._sync_var(arrival_source)

        # create source-type attributes marking arrivals
        for source in sources:
            setattr(self, 'is_{}'.format(source.decode()), self._arrival_source == self._source_enum[source])
        # end for

        # attribute marking automatic picks
        setattr(self, 'is_AUTO', self.arrivals['quality_measure_slope'] > -1)
    # end func

    def _coalesce_network_codes(self):
        self.arrivals['net'] = np.load('coalesced_net.npy')
        return

        MDIST = 1e3 # maximum distance in metres
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
        if(self.rank==0): np.save('coalesced_net.npy', self.arrivals['net'])
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
    else:
        sr = SSSTRelocator('./merge_catalogues_output.csv',
                           auto_pick_files=['p_combined.txt', 's_combined.txt'],
                           auto_pick_phases=['P', 'S'],
                           events_only=False)

        if(sr.rank==0): print('computing residuals..')
        sr.compute_residual()

        print('Rank {}: memory used: {}'.format(sr.rank,
                                                round(psutil.Process().memory_info().rss / 1024. / 1024., 2)))

        if(0):
            for name, dtype in zip(sr.arrival_fields['names'], sr.arrival_fields['formats']):
                test_var = sr._sync_var(sr.arrivals[name][sr.local_arrivals_indices])

                assert np.all(test_var == sr.arrivals[name]), 'Sync-{} failed on rank {}'.format(name, sr.rank)
            # end for
        # end if
    # end if
# end if
