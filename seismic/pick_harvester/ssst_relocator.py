from ordered_set import OrderedSet as set
import numpy as np
from obspy import UTCDateTime
from scipy.spatial import cKDTree
import pyproj
from collections import defaultdict
from seismic.pick_harvester.utils import split_list

from seismic.pick_harvester.parametric_data import ParametricData
import importlib
import sys, os
from tqdm import tqdm
from seismic.pick_harvester.travel_time import TTInterpolator
from seismic.pick_harvester.ellipticity import Ellipticity
import traceback
import psutil
import h5py

class SSSTRelocator(ParametricData):
    def __init__(self, csv_catalog, auto_pick_files=[], auto_pick_phases=[],
                 events_only=False, phase_list='P Pg Pb Pn S Sg Sb Sn', temp_dir='./'):
        super(SSSTRelocator, self).__init__(csv_catalog, auto_pick_files, auto_pick_phases,
                                            events_only, phase_list, temp_dir)
        self.BALL_RADIUS_KM = 11 #~0.1 arc-degrees
        self._lonlatalt2xyz = None
        self.vsurf = None
        self.residual = None
        self._source_enum = defaultdict(int)
        self._arrival_source = None
        self._tti = TTInterpolator()
        self._ellipticity = Ellipticity()
        self.station_id = None
        self.ssst_tcorr = None

        np.random.seed(0)
        np.random.shuffle(self.arrivals)

        self._assign_station_ids()
        self._label_arrivals()
        self._initialize_spatial_functors()

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

    def _assign_station_ids(self):
        if(self.rank == 0): print('Assigning station IDs..')

        netsta_ids = None
        if(self.rank == 0):
            nets = self.arrivals['net']
            stas = self.arrivals['sta']

            netstas = np.unique(np.char.add(np.char.add(nets, b'.'), stas))
            netsta_ids = {}
            for sid, netsta in enumerate(netstas):
                netsta_ids[netsta] = sid
            # end for
        # end if
        netsta_ids = self.comm.bcast(netsta_ids)

        station_id = np.zeros(len(self.local_arrivals_indices), dtype='i4')
        local_nets = self.arrivals['net'][self.local_arrivals_indices]
        local_stas = self.arrivals['sta'][self.local_arrivals_indices]
        local_netstas = np.char.add(np.char.add(local_nets, b'.'), local_stas)
        for i, local_netsta in enumerate(local_netstas):
            station_id[i] = netsta_ids[local_netsta]
        # end for

        self.station_id = self._sync_var(station_id)
    # end func

    def _compute_residual_helper(self, phase, imask=None):
        """
        Operates on local data

        @param phase: dim(local_arrival_indices)
        @param imask: dim(local_arrival_indices)
        @return:
        """
        if (imask is None): imask = np.ones(len(phase), dtype='?')

        # event indices for local arrivals
        indices = self.event_id_to_idx[self.arrivals['event_id'][self.local_arrivals_indices]]
        elon = self.events['lon'][indices][imask]
        elat = self.events['lat'][indices][imask]
        edepth_km = self.events['depth_km'][indices][imask]

        slon = self.arrivals['lon'][self.local_arrivals_indices][imask]
        slat = self.arrivals['lat'][self.local_arrivals_indices][imask]
        vsurf = self.vsurf[self.local_arrivals_indices][imask]
        elev_km = self.arrivals['elev_m'][self.local_arrivals_indices][imask] / 1e3

        azim, _, ecdist = self._geod.inv(elon, elat, slon, slat)

        elev_corr = self.compute_elevation_correction(phase, vsurf, elev_km, ecdist, edepth_km)
        ellip_corr = self.compute_ellipticity_correction(phase, ecdist, edepth_km, elat, azim)
        ptt = self._tti.get_tt(phase, ecdist, edepth_km) + elev_corr + ellip_corr

        # compute empirical tt
        ett = self.arrivals['arrival_ts'][self.local_arrivals_indices][imask] - \
              self.events['origin_ts'][indices][imask]

        residual = ett - ptt.astype('f8')

        if (self.rank == 0): print(elev_corr, ellip_corr, residual)

        return residual
    # end func

    def compute_residual(self):
        phase = self.arrivals['phase'][self.local_arrivals_indices]

        local_residual = self._compute_residual_helper(phase)

        self.residual = self._sync_var(local_residual)
    # end func

    def redefine_phases(self, imask=None):
        """
        Operates on local data, followed by a global sync

        @param imask: dim(local_arrival_indices)
        """
        if(imask is None): imask = np.ones(len(self.local_arrivals_indices), dtype='?')

        phase = self.arrivals['phase'][self.local_arrivals_indices]
        is_p = self.is_P(phase) & imask
        is_s = self.is_S(phase) & imask
        p_indices = np.argwhere(is_p).flatten()
        s_indices = np.argwhere(is_s).flatten()

        p_residuals = np.zeros((len(p_indices), len(self.p_phases)))
        s_residuals = np.zeros((len(s_indices), len(self.s_phases)))

        test_p_phase = np.zeros(len(p_indices), self.arrivals['phase'].dtype)
        for ip, p in enumerate(self.p_phases):
            test_p_phase[:] = p
            p_residuals[:, ip] = self._compute_residual_helper(test_p_phase, imask=is_p)
        # end for

        test_s_phase = np.zeros(len(s_indices), self.arrivals['phase'].dtype)
        for ip, p in enumerate(self.s_phases):
            test_s_phase[:] = p
            s_residuals[:, ip] = self._compute_residual_helper(test_s_phase, imask=is_s)
        # end for

        local_new_phase = np.array(self.arrivals['phase'][self.local_arrivals_indices])
        local_new_residual = np.array(self.residual[self.local_arrivals_indices])

        local_new_phase[p_indices] = b'Px'
        local_new_phase[s_indices] = b'Sx'

        local_new_phase[p_indices] = self.p_phases[np.argmin(np.fabs(p_residuals), axis=1)]
        local_new_phase[s_indices] = self.s_phases[np.argmin(np.fabs(s_residuals), axis=1)]

        p_argmin_indices = np.argmin(np.fabs(p_residuals), axis=1)
        s_argmin_indices = np.argmin(np.fabs(s_residuals), axis=1)
        local_new_residual[p_indices] = np.take_along_axis(p_residuals, p_argmin_indices[:, None], axis=1).flatten()
        local_new_residual[s_indices] = np.take_along_axis(s_residuals, s_argmin_indices[:, None], axis=1).flatten()

        self.arrivals['phase'] = self._sync_var(local_new_phase)
        self.residual = self._sync_var(local_new_residual)
    # end func

    def compute_elevation_correction(self, phase, vsurf, elev_km, ecdist, edepth_km):
        """
        Operates on input parameters without respect to local/global indices

        @param phase: dim(any)
        @param vsurf:  dim(any)
        @param elev_km:  dim(any)
        @param ecdist:  dim(any)
        @param edepth_km:  dim(any)
        @return:
        """
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
        """
        Operates on input parameters without respect to local/global indices

        @param phase:  dim(any)
        @param ecdist:  dim(any)
        @param edepth_km:  dim(any)
        @param elat:  dim(any)
        @param azim:  dim(any)
        @return:
        """
        ellip_corr = self._ellipticity.get_correction(phase, ecdist, edepth_km, elat, azim)
        ellip_corr = np.array(ellip_corr.astype('f4'))

        return ellip_corr
    # end func

    def compute_ssst_correction(self, imask=None):
        """
        Operates on global data
        @param imask: dim(len(arrivals))
        """
        if(imask is None): imask = np.ones(len(self.arrivals), dtype='?')

        stas = np.unique(self.station_id)
        local_sids = split_list(stas, self.nproc)[self.rank]
        local_corrs = np.zeros(len(self.arrivals), dtype='f4')

        for sid in tqdm(local_sids, desc='Rank {}: '.format(self.rank)):
            sta_arrival_indices = np.argwhere((self.station_id == sid) & (imask)).flatten()
            sta_arrivals = self.arrivals[sta_arrival_indices]

            phases = set(sta_arrivals['phase'])
            phase_trees = {}
            phase_residuals = {}
            phase_arrival_indices = {}
            for p in phases:
                pimask = sta_arrivals['phase'] == p

                event_indices = self.event_id_to_idx[sta_arrivals['event_id'][pimask]]
                elons = self.events['lon'][event_indices]
                elats = self.events['lat'][event_indices]
                ealts = -self.events['depth_km'][event_indices] * 1e3  # m

                xyz = self._lonlatalt2xyz(elons, elats, ealts)
                phase_trees[p] = cKDTree(xyz)
                phase_residuals[p] = self.residual[sta_arrival_indices][pimask]
                phase_arrival_indices[p] = sta_arrival_indices[pimask]
            # end for

            for i, sta_arrival in enumerate(sta_arrivals):
                idx = self.event_id_to_idx[sta_arrival['event_id']]
                elon = self.events['lon'][idx]
                elat = self.events['lat'][idx]
                ealt = -self.events['depth_km'][idx] * 1e3  # m
                cxyz = self._lonlatalt2xyz(elon, elat, ealt).flatten()
                indices = phase_trees[sta_arrival['phase']].query_ball_point(cxyz,
                                                                             self.BALL_RADIUS_KM * 1e3)
                corr = np.median(phase_residuals[sta_arrival['phase']][indices])
                local_corrs[sta_arrival_indices[i]] = corr
            # end for
        # end for

        self.ssst_tcorr = self._sum(local_corrs)
    # end func

    def _initialize_spatial_functors(self, ellipsoidal_distance=False):
        if(self.rank == 0): print('Initializing spatial functors..')

        ER = self.EARTH_RADIUS_KM * 1e3 #m

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

    def _sum(self, rank_values):
        # sum of local values across ranks
        counts = np.array(self.comm.allgather(len(rank_values)), dtype='i4')
        assert(np.all(counts[0]==counts))

        nelem = counts[0]
        dtype = rank_values.dtype
        global_values = np.zeros(nelem, dtype=dtype)

        fn = self._temp_dir + '/sum.h5'

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
            global_values[:] += hf['{}'.format(irank)][:]
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

        if(sr.rank == 0): old_phase = np.array(sr.arrivals['phase'])

        if(sr.rank==0): print('Computing residuals..')
        sr.compute_residual()

        if(sr.rank==0): print('Redefining phases..')
        sr.redefine_phases(imask=~sr.is_ISC[sr.local_arrivals_indices] & sr.is_AUTO[sr.local_arrivals_indices])

        if(sr.rank==0): print('Computing SSST corrections..')
        sr.compute_ssst_correction(imask=~sr.is_ISC & sr.is_AUTO)
        if(sr.rank==0):
            np.save('test_arrivals.npy', sr.arrivals)
            np.save('ssst_tcorr.npy', sr.ssst_tcorr)
        # end if

        if(sr.rank == 0):
            filt = ~sr.is_ISC & sr.is_AUTO
            # filt = np.ones(len(sr.is_ISC), dtype='?')
            for p in np.concatenate((sr.p_phases, sr.s_phases, np.array(['Px', 'Sx'], dtype='U2'))):
                p = p.encode()
                before = np.sum(old_phase[filt] == p)
                after = np.sum(sr.arrivals['phase'][filt] == p)
                change = 0
                if (before): change = (before - after) / float(before) * 100
                print('Auto phase {}: before: {} after: {} change: {}%'.format(p, before, after, change))
            # end for
        # end if

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
